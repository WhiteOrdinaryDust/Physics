'use strict';
/*
 * particle_probe_v3.js
 * ---------------------------------------------------------------------------
 * STEP 1 (repaired): a detector whose controls are BOTH clean.
 *
 * Previous 1D detector failed its NEGATIVE control (rule 30 random reported
 * "particles") because it matched tiles at every position on the whole ring,
 * which manufactures coincidences. It also tracked connected runs, which is the
 * wrong granularity: a rule-110 glider is not a single run.
 *
 * THE FIX (background-relative detection)
 *   1. Find the equilibrium background phase of the rule empirically (do not
 *      assume it).
 *   2. Call a cell "excited" if it differs from that background phase.
 *   3. A particle is a maximal connected run of excited cells.
 *   4. Match those structures over period p and require constant non-zero drift.
 *
 * This mirrors the Life detector (connected components), which passes its
 * positive control (glider -> p4:v(1,1)) and negative control (blinker -> none).
 *
 * NOTE ON VELOCITY LIMITS
 *   In a radius-1 CA, drift over period p is bounded by |v*p| <= p, so the period
 *   window IS the velocity window. Neither detector needs an arbitrary v cutoff.
 * ---------------------------------------------------------------------------
 */

// ============================================================ ECA core
function ruleTable(rule) {
  const t = new Uint8Array(8);
  for (let n = 0; n < 8; n++) t[n] = (rule >> n) & 1;
  return t;
}

function step(row, W, tbl) {
  const out = new Uint8Array(W);
  for (let i = 0; i < W; i++) {
    out[i] = tbl[(row[(i - 1 + W) % W] << 2) | (row[i] << 1) | row[(i + 1) % W]];
  }
  return out;
}

function runRows(rule, init, W, T) {
  const tbl = ruleTable(rule);
  const rows = new Uint8Array((T + 1) * W);
  let cur = Uint8Array.from(init);
  rows.set(cur, 0);
  for (let t = 0; t < T; t++) {
    cur = step(cur, W, tbl);
    rows.set(cur, (t + 1) * W);
  }
  return rows;
}

// ---- step 1: find the background, empirically, WITH A CREDIBILITY GATE
//
// A background is only meaningful if a large fraction of cells actually agree
// with it. Without this gate, a chaotic rule (rule 30) "finds" a period-14 phase
// at 50% agreement, which marks EVERY cell as excited and makes the detector
// hallucinate structure in pure noise. That was the negative-control failure.
const BG_MIN_AGREEMENT = 0.90;

function findBackground(rows, W, T, maxP) {
  let best = null;
  for (let P = 1; P <= maxP; P++) {
    for (let p = 0; p < P; p++) {
      let agree = 0, total = 0;
      for (let t = P; t <= T; t++) {
        const base = t * W, prev = (t - P) * W;
        for (let i = 0; i < W; i++) {
          if (rows[prev + ((i + p) % W)] === rows[base + i]) agree++;
          total++;
        }
      }
      const frac = agree / total;
      if (!best || frac > best.frac) best = { P, p, frac };
    }
  }
  best.credible = best.frac >= BG_MIN_AGREEMENT;
  return best;
}

// ---- step 2-3: excited runs at time t
function excitedRuns(rows, W, t, bg) {
  const { P, p } = bg;
  const refT = t - P;
  if (refT < 0) return [];
  const base = t * W, prev = refT * W;
  const exc = new Uint8Array(W);
  for (let i = 0; i < W; i++) {
    exc[i] = (rows[base + i] !== rows[prev + ((i + p) % W)]) ? 1 : 0;
  }
  // maximal runs, wrap-aware
  let z = -1;
  for (let i = 0; i < W; i++) if (exc[i] === 0) { z = i; break; }
  if (z === -1) return [{ start: 0, len: W, pattern: null }];
  const out = [];
  let i = 0;
  while (i < W) {
    const idx = (z + i) % W;
    if (exc[idx] === 1) {
      let len = 0;
      while (len < W && exc[(z + i + len) % W] === 1) len++;
      const start = (z + i) % W;
      let pattern = '';
      for (let j = 0; j < len; j++) pattern += rows[base + ((start + j) % W)];
      out.push({ start, len, pattern });
      i += len;
    } else i++;
  }
  return out;
}

// ---- step 4: match structures over period p with constant non-zero drift
function detect(rows, W, T, bg, minLen, maxPeriod, minRepeats) {
  // no credible background => this rule has nothing to be a defect IN
  if (!bg.credible) return new Map();
  const history = [];   // array of Map: pattern -> start
  const hits = new Map();
  for (let t = 0; t <= T; t++) {
    const frame = new Map();
    for (const r of excitedRuns(rows, W, t, bg)) {
      if (r.len < minLen) continue;
      const key = `${r.len}:${r.pattern}`;
      if (!frame.has(key)) frame.set(key, r.start);
    }
    history.push(frame);
    if (history.length > maxPeriod + 1) history.shift();

    for (const [key, x] of frame) {
      for (let p = 1; p <= maxPeriod; p++) {
        const old = history[history.length - 1 - p];
        if (!old) continue;
        const x0 = old.get(key);
        if (x0 === undefined) continue;
        let dx = x - x0;
        if (dx > W / 2) dx -= W;
        if (dx < -W / 2) dx += W;
        if (dx === 0) continue;                 // must move
        if (Math.abs(dx) > p) continue;         // light-speed bound
        const label = `${key.split(':')[0]}:p${p}:v${dx}`;
        hits.set(label, (hits.get(label) || 0) + 1);
      }
    }
  }
  return hits;
}

// ============================================================ tests
const W = 400, T = 400;

console.log('=== 1D detector, background-relative, both controls ===\n');

// ---- rule 110, seeded glider
{
  const ether = [0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0];
  const init = new Uint8Array(W);
  for (let i = 0; i < W; i++) init[i] = ether[i % ether.length];
  const g = [1, 1, 1, 0, 1, 1, 1];
  for (let j = 0; j < g.length; j++) init[(150 + j) % W] ^= g[j];
  const rows = runRows(110, init, W, T);
  const bg = findBackground(rows, W, T, 20);
  console.log(`rule 110 seeded glider`);
  console.log(`  background found: period ${bg.P}, phase ${bg.p}, agreement ${(bg.frac * 100).toFixed(2)}%`);
  const hits = detect(rows, W, T, bg, 3, 6, 1);
  const ent = [...hits.entries()].sort((a, b) => b[1] - a[1]);
  if (ent.length === 0) console.log('  *** FAILED: no drifting structure ***');
  else {
    for (const [k, v] of ent.slice(0, 8)) console.log(`    ${k.padEnd(22)} ${v}`);
    // ---- visual confirmation: show the defect at three times
    console.log('  visual check (reference row vs actual row; "|" marks differing cells):');
    const { P, p } = bg;
    for (const dt of [0, 40, 80]) {
      const t = 60 + dt;
      if (t > T) continue;
      const base = t * W, prev = (t - P) * W;
      let first = -1;
      for (let i = 0; i < W; i++) {
        if (rows[base + i] !== rows[prev + ((i + p) % W)]) { first = i; break; }
      }
      if (first < 0) { console.log(`    t=${t}: no defect`); continue; }
      const lo = (first - 6 + W) % W;
      let ref = '', act = '', mark = '';
      for (let j = 0; j < 26; j++) {
        const i = (lo + j) % W;
        ref += rows[prev + ((i + p) % W)];
        act += rows[base + i];
        mark += (rows[base + i] !== rows[prev + ((i + p) % W)]) ? '|' : ' ';
      }
      console.log(`    t=${t}  ref ${ref}`);
      console.log(`          act ${act}`);
      console.log(`              ${mark}   (defect starts near index ${first})`);
    }
  }
  console.log('');
}

// ---- rule 110, random (mixed: gliders arise from defects)
{
  let s = 999 >>> 0;
  const rnd = () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
  const init = new Uint8Array(W);
  for (let i = 0; i < W; i++) init[i] = rnd() < 0.5 ? 1 : 0;
  const rows = runRows(110, init, W, T);
  const bg = findBackground(rows, W, T, 20);
  console.log(`rule 110 random start`);
  console.log(`  background found: period ${bg.P}, phase ${bg.p}, agreement ${(bg.frac * 100).toFixed(2)}%`);
  const hits = detect(rows, W, T, bg, 3, 6, 1);
  const ent = [...hits.entries()].sort((a, b) => b[1] - a[1]);
  if (ent.length === 0) console.log('  no drifting structures');
  else for (const [k, v] of ent.slice(0, 8)) console.log(`    ${k.padEnd(22)} ${v}`);
  console.log('');
}

// ---- rule 30 random: NEGATIVE control, must be clean
{
  let s = 12345 >>> 0;
  const rnd = () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
  const init = new Uint8Array(W);
  for (let i = 0; i < W; i++) init[i] = rnd() < 0.5 ? 1 : 0;
  const rows = runRows(30, init, W, T);
  const bg = findBackground(rows, W, T, 20);
  console.log(`rule 30 random (NEGATIVE CONTROL)`);
  console.log(`  background found: period ${bg.P}, phase ${bg.p}, agreement ${(bg.frac * 100).toFixed(2)}%`);
  const hits = detect(rows, W, T, bg, 3, 6, 1);
  const ent = [...hits.entries()].sort((a, b) => b[1] - a[1]);
  if (ent.length === 0) console.log('  PASS: no drifting structures');
  else for (const [k, v] of ent.slice(0, 8)) console.log(`    ${k.padEnd(22)} ${v}`);
  console.log('');
}

// ---- rule 184 (traffic): diagnose whether the background phase itself drifts
{
  const W2 = 400, T2 = 400;
  let s = 777 >>> 0;
  const rnd = () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
  const init = new Uint8Array(W2);
  for (let i = 0; i < W2; i++) init[i] = (i % 17 === 0) ? 1 : 0;
  const rows = runRows(184, init, W2, T2);
  const bg = findBackground(rows, W2, T2, 20);
  console.log(`rule 184 (traffic, known to have particles)`);
  console.log(`  background: period ${bg.P}, phase ${bg.p}, agreement ${(bg.frac * 100).toFixed(2)}%, credible=${bg.credible}`);

  // Diagnostic: is the background phase STABLE, or does it drift?
  // Compare the best phase at two different time windows.
  const bestPhaseIn = (t0, t1) => {
    let best = null;
    for (let P = 1; P <= 20; P++) {
      for (let p = 0; p < P; p++) {
        let agree = 0, total = 0;
        for (let t = Math.max(P, t0); t <= t1; t++) {
          const base = t * W2, prev = (t - P) * W2;
          for (let i = 0; i < W2; i++) {
            if (rows[prev + ((i + p) % W2)] === rows[base + i]) agree++;
            total++;
          }
        }
        const frac = agree / total;
        if (!best || frac > best.frac) best = { P, p, frac };
      }
    }
    return best;
  };
  const early = bestPhaseIn(20, 150);
  const late = bestPhaseIn(250, 380);
  console.log(`  phase in t=20..150  : P=${early.P} p=${early.p} agree=${(early.frac * 100).toFixed(2)}%`);
  console.log(`  phase in t=250..380 : P=${late.P} p=${late.p} agree=${(late.frac * 100).toFixed(2)}%`);
  console.log(`  -> phase ${early.p === late.p && early.P === late.P ? 'STABLE' : 'DRIFTS'}: global-phase detection is ${
    early.p === late.p && early.P === late.P ? 'applicable' : 'NOT applicable (the background itself moves)'}`);
}
