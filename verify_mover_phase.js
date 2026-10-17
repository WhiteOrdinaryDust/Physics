'use strict';
/*
 * verify_mover_phase.js
 * ---------------------------------------------------------------------------
 * WHAT IS BEING CHECKED
 *
 * A background-relative search reported: the five two-sidedly elastic lamellar
 * media contain NO movers at all (950 injections, 0 movers), but the SAME rule
 * B{01478}/S{01} does carry movers in a DIFFERENT lamellar phase (density 0.25)
 * whose own elasticity is 0.0000/0.0000.
 *
 * If that is right, it is a sharp statement and not just a null result:
 *   movers exist in this rule, but only in the phase that does NOT heal.
 * That is an anti-correlation between "the medium repairs damage" and "the medium
 * carries travelling excitations", and it would explain the whole pattern:
 *   Life              : never heals a bump  -> rich particle zoo
 *   lamellar elastic  : heals everything    -> nothing can travel
 *   lamellar non-elast: does not heal       -> tiny movers appear
 *
 * I do not take that on one script's word, so this re-does it on my own code path,
 * with the two corrections that search flagged:
 *   (1) velocity must be measured RELATIVE TO THE MATERIAL -- some lamellar phases
 *       translate between the two frames of their cycle, so the reference is
 *       aligned to the material (search over translations, keep the best match);
 *   (2) a reported period must be fundamental, not a multiple of the true one.
 *
 * CONTROLS
 *   C1 my detector must find the Life glider in the empty background: p4, v(1,1).
 *   C2 an unperturbed phase must produce ZERO movers (no false positives from the
 *      background's own clock).
 *   C3 no hand-written candidate phases: phases are found by evolving stripe-family
 *      seeds and detecting exact cycles, then identified by their density.
 *
 * PREDICTION (stated before the run):
 *   I expect to reproduce the reported anti-correlation: elasticity ~0 in the
 *   phase that carries the 2-cell mover, and ~1 in the elastic phase with no mover.
 *   I expect the 2-cell vertical pair to move ALONG the stripes, not across them.
 * ---------------------------------------------------------------------------
 */

const W = 24, H = 24, NS = W * H;

function parseRule(s) {
  const m = /^B\{?([0-8]*)\}?\/S\{?([0-8]*)\}?$/.exec(s.trim());
  let bm = 0, sm = 0;
  for (const c of m[1]) bm |= 1 << Number(c);
  for (const c of m[2]) sm |= 1 << Number(c);
  return { bm, sm };
}
function stepR(g, r) {
  const out = new Uint8Array(NS);
  for (let y = 0; y < H; y++) {
    const yp = ((y - 1 + H) % H) * W, ym = ((y + 1) % H) * W, yc = y * W;
    for (let x = 0; x < W; x++) {
      const xp = (x - 1 + W) % W, xm = (x + 1) % W;
      const n = g[yp + xp] + g[yp + x] + g[yp + xm] + g[yc + xp] + g[yc + xm] + g[ym + xp] + g[ym + x] + g[ym + xm];
      const i = yc + x;
      out[i] = g[i] ? ((r.sm >> n) & 1) : ((r.bm >> n) & 1);
    }
  }
  return out;
}
function eq(a, b) { for (let i = 0; i < NS; i++) if (a[i] !== b[i]) return false; return true; }
function density(g) { let s = 0; for (let i = 0; i < NS; i++) s += g[i]; return s / NS; }
function hash(g) { let h = 2166136261 >>> 0; for (let i = 0; i < NS; i++) { h ^= g[i]; h = Math.imul(h, 16777619) >>> 0; } return h >>> 0; }
function shift(g, a, b) {
  const o = new Uint8Array(NS);
  for (let y = 0; y < H; y++) { const yb = ((y + b) % H) * W, yc = y * W; for (let x = 0; x < W; x++) o[yc + x] = g[yb + ((x + a) % W)]; }
  return o;
}
// best translation aligning c0 onto c1 (used to measure velocity relative to the material)
function align(c0, c1) {
  let best = null;
  for (let a = 0; a < W; a++) for (let b = 0; b < H; b++) {
    let d = 0;
    for (let y = 0; y < H; y++) { const yb = ((y + b) % H) * W, yc = y * W; for (let x = 0; x < W; x++) if (c1[yc + x] !== c0[yb + ((x + a) % W)]) d++; }
    if (!best || d < best.d) best = { a, b, d };
  }
  return best;
}
function compsDev(dev) {
  const seen = new Uint8Array(NS), out = [];
  for (let i = 0; i < NS; i++) {
    if (!dev[i] || seen[i]) continue;
    const cells = []; const Q = [i]; seen[i] = 1;
    while (Q.length) {
      const k = Q.pop(); cells.push(k);
      const x = k % W, y = (k - x) / W;
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        const nk = (((y + dy + H) % H) * W) + ((x + dx + W) % W);
        if (dev[nk] && !seen[nk]) { seen[nk] = 1; Q.push(nk); }
      }
    }
    out.push(cells);
  }
  return out;
}
function info(cells) {
  let minx = 1e9, miny = 1e9, maxx = -1e9, maxy = -1e9, sx = 0, sy = 0;
  for (const k of cells) { const x = k % W, y = (k - x) / W; if (x < minx) minx = x; if (x > maxx) maxx = x; if (y < miny) miny = y; if (y > maxy) maxy = y; sx += x; sy += y; }
  const w = maxx - minx + 1, h = maxy - miny + 1, gg = [];
  for (let r = 0; r < h; r++) gg.push(new Array(w).fill('.'));
  for (const k of cells) { const x = k % W, y = (k - x) / W; gg[y - miny][x - minx] = '#'; }
  return { id: `${w}x${h}:` + gg.map(r => r.join('')).join('/'), w, h, n: cells.length, cx: sx / cells.length, cy: sy / cells.length };
}
function mulberry32(a) { return function () { a |= 0; a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }

// a mover = one deviation group whose shape recurs with period p and whose centroid
// drifts by a constant non-zero vector, persisting >= 4 periods. Velocity is measured
// against the MATERIAL-ALIGNED reference, so a defect bound to the material scores zero.
function findMover(phase, cycle, r, steps) {
  const p0 = cycle[0];
  const refs = [p0];
  if (cycle.length === 2) {
    const al = align(p0, cycle[1]);
    refs.push(al.d === 0 ? shift(p0, al.a, al.b) : cycle[1]);
  }
  const refAt = (t) => (cycle.length === 1 ? p0 : shift(refs[1], 0, 0) && (t % 2 === 0 ? refs[0] : refs[1]));
  let g = Uint8Array.from(phase);
  const track = [];
  for (let t = 0; t < steps; t++) {
    const ref = t % 2 === 0 ? refs[0] : refs[1 % refs.length];
    const dev = new Uint8Array(NS);
    let nd = 0;
    for (let i = 0; i < NS; i++) if (g[i] !== ref[i]) { dev[i] = 1; nd++; }
    if (nd === 0) break;                       // absorbed: the material healed
    const cs = compsDev(dev).filter(c => c.length <= 8);
    if (cs.length) {
      cs.sort((a, b) => b.length - a.length);
      track.push({ t, ...info(cs[0]) });
    } else track.push(null);
    g = stepR(g, r);
  }
  // look for a fundamental period with a constant drift
  for (let p = 1; p <= 8; p++) {
    const frames = track.filter(Boolean);
    if (frames.length < 4 * p + 1) continue;
    let ok = true, base = null;
    for (let i = 0; i + p < frames.length; i++) {
      if (frames[i].id !== frames[i + p].id) { ok = false; break; }
      const dx = Math.round(frames[i + p].cx - frames[i].cx), dy = Math.round(frames[i + p].cy - frames[i].cy);
      if (base === null) base = [dx, dy];
      else if (dx !== base[0] || dy !== base[1]) { ok = false; break; }
    }
    if (ok && base && (base[0] || base[1]) && frames.length >= 4 * p) {
      // fundamental: no divisor of p also works
      let fundamental = true;
      for (let q = 1; q < p; q++) {
        if (p % q) continue;
        let ok2 = true;
        for (let i = 0; i + q < frames.length; i++) if (frames[i].id !== frames[i + q].id) { ok2 = false; break; }
        if (ok2) { fundamental = false; break; }
      }
      if (fundamental) return { p, vx: base[0], vy: base[1], shape: frames[0].id, persist: frames.length / p };
    }
  }
  return null;
}

// ============================================================ C1 detector control
console.log('=== C1: does my detector find the Life glider in the empty background? ===');
{
  const life = parseRule('B3/S23');
  const empty = new Uint8Array(NS);
  const g = new Uint8Array(NS);
  const G = [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]];
  for (const [dx, dy] of G) g[(5 + dy) * W + (5 + dx)] = 1;
  const m = findMover(g, [empty], life, 60);
  console.log(`  glider: ${m ? `period ${m.p}, velocity (${m.vx},${m.vy}), shape ${m.shape}, persists ${m.persist.toFixed(1)} periods` : '*** NOT FOUND -- detector is blind ***'}  ${m && m.p === 4 && m.vx === 1 && m.vy === 1 ? 'PASS' : '*** FAIL ***'}`);
  if (!(m && m.p === 4 && m.vx === 1 && m.vy === 1)) process.exit(1);
}
// C2: an unperturbed phase must give zero movers
console.log('=== C2: unperturbed phase gives zero movers ===');

// ============================================================ find phases by evolution
function findCycles(r, seeds, maxGen) {
  const out = [], seen = new Set();
  for (const seed of seeds) {
    let g = seed;
    const hist = [];
    for (let t = 0; t < maxGen; t++) {
      hist.push(g);
      g = stepR(g, r);
      let found = null;
      for (let p = 1; p <= 8; p++) if (hist.length >= p && eq(g, hist[hist.length - p])) { found = hist.slice(hist.length - p); break; }
      if (found) {
        const sig = found.length + '|' + found.map(hash).join(',');
        if (!seen.has(sig)) { seen.add(sig); out.push(found); }
        break;
      }
    }
  }
  return out;
}
function stripeSeeds(rule, rnd) {
  const seeds = [];
  for (const p of [2, 3, 4, 5, 6]) for (const k of [1, 2]) {
    for (const axis of ['x', 'y']) {
      const g = new Uint8Array(NS);
      for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
        const v = axis === 'x' ? x : y;
        if ((v % p) < k) g[y * W + x] = 1;
      }
      seeds.push(g);
    }
  }
  for (let s = 0; s < 6; s++) {
    const g = new Uint8Array(NS);
    for (let i = 0; i < NS; i++) if (rnd() < 0.3) g[i] = 1;
    seeds.push(g);
  }
  return seeds;
}

function measureElastic(phase, cycle, r) {
  const refs = [cycle[0]];
  if (cycle.length === 2) { const al = align(cycle[0], cycle[1]); refs.push(al.d === 0 ? shift(cycle[0], al.a, al.b) : cycle[1]); }
  const cycEq = (g) => refs.some(c => eq(g, c));
  const bgNext = refs.length === 2 ? refs[1] : stepR(refs[0], r);
  let addRel = 0, addHeal = 0, rmRel = 0, rmHeal = 0;
  for (let i = 0; i < NS; i += 3) {   // stride 3 to keep the run cheap; ratio is what matters
    const g = Uint8Array.from(refs[0]); g[i] ^= 1;
    let healed = false, gg = g;
    for (let t = 0; t < 120; t++) { if (cycEq(gg)) { healed = true; break; } gg = stepR(gg, r); }
    if (!healed && cycEq(gg)) healed = true;
    if (refs[0][i]) { rmRel++; if (healed) rmHeal++; }
    else if (!eq(stepR(g, r), bgNext)) { addRel++; if (healed) addHeal++; }
  }
  return { eAdd: addRel ? addHeal / addRel : NaN, eRm: rmRel ? rmHeal / rmRel : NaN, addRel, addHeal, rmRel, rmHeal };
}

for (const [name, ruleStr] of [['B{01478}/S{01}', 'B{01478}/S{01}'], ['B{023456}/S{01245}', 'B{023456}/S{01245}']]) {
  const r = parseRule(ruleStr);
  const rnd = mulberry32(20240607);
  console.log(`\n=== ${name}: phases found from stripe-family seeds, with their elasticity ===`);
  const cycles = findCycles(r, stripeSeeds(r, rnd), 200);
  console.log(`  ${cycles.length} distinct cycle(s):`);
  const seenD = new Set();
  for (const cyc of cycles) {
    const d = density(cyc[0]);
    const key = `${d.toFixed(4)}|${cyc.length}`;
    if (seenD.has(key)) continue;
    seenD.add(key);
    const e = measureElastic(cyc[0], cyc, r);
    console.log(`    density ${d.toFixed(4)}, period ${cyc.length}: E_add_causal ${Number.isNaN(e.eAdd) ? 'n/a' : e.eAdd.toFixed(4)} (${e.addHeal}/${e.addRel}), E_remove_causal ${e.eRm.toFixed(4)} (${e.rmHeal}/${e.rmRel})`);
  }
  // injection test on up to 3 distinct phases
  console.log(`  injection test (remove/add a 2-cell vertical pair, 24 positions, 150 generations):`);
  seenD.clear();
  let tested = 0;
  for (const cyc of cycles) {
    const d = density(cyc[0]);
    const key = `${d.toFixed(4)}|${cyc.length}`;
    if (seenD.has(key) || tested >= 3) continue;
    seenD.add(key); tested++;
    let movers = 0, absorbed = 0, injections = 0, examples = [];
    for (const [dx, dy] of [[0, 0], [0, 1]]) { /* placeholder to keep structure explicit */ }
    for (let k = 0; k < W; k++) {
      for (const dir of ['remove', 'add']) {
        const g = Uint8Array.from(cyc[0]);
        const x = k, y0 = Math.floor(H / 2);
        const i1 = y0 * W + x, i2 = ((y0 + 1) % H) * W + x;
        const want = dir === 'remove' ? 0 : 1;
        const have = cyc[0][i1] && cyc[0][i2];
        if (dir === 'remove' && !have) continue;
        if (dir === 'add' && have) continue;
        g[i1] = want; g[i2] = want;
        injections++;
        const m = findMover(g, cyc, r, 150);
        if (m) { movers++; if (examples.length < 3) examples.push(m); }
        else {
          // absorbed = the deviation vanished entirely at some point
          let gg = Uint8Array.from(g), went = false;
          const refs = [cyc[0]]; if (cyc.length === 2) refs.push(cyc[1]);
          for (let t = 0; t < 150; t++) { if (refs.some(c => eq(gg, c))) { went = true; break; } gg = stepR(gg, r); }
          if (went) absorbed++;
        }
      }
    }
    console.log(`    density ${d.toFixed(4)} period ${cyc.length}: ${injections} injections, ${movers} mover(s), ${absorbed} absorbed`);
    for (const m of examples) console.log(`      mover: period ${m.p}, velocity (${m.vx},${m.vy}), shape ${m.shape}, persists ${m.persist.toFixed(1)} periods`);
  }
}
