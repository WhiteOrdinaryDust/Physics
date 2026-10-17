'use strict';
/*
 * elastic_blob_impact.js
 * ---------------------------------------------------------------------------
 * THE TEST THAT DECIDES WHETHER THE NEW ELASTIC MATTER IS USEFUL.
 *
 * crystalline_solid_search.js found 5 rules whose background is a crystal AND
 * heals both a one-cell hole and a one-cell bump with real samples (48-96 flips).
 * That is the first verified two-sided-elastic matter in this programme, and it is
 * what a screen needs in order to be bumped without recording the bump.
 *
 * But a PARTICLE is not one cell. The LWSS is 9-12 cells; a glider is 5. A medium
 * that repairs a one-cell scratch may still be destroyed by a particle-sized blow,
 * and if so the whole point is lost: we would again have a detector, not a screen.
 *
 * So this script repeats the elasticity measurement with blobs of k = 1, 2, 3, 5
 * cells, added or removed, and reports how the heal rate falls with k.
 *
 * PREDICTION (stated before the run):
 *   - Life block lattice (control): heal rate 0 for additions at every k (already
 *     known: a bump is fatal), 1.0 for k=1 removals, and it should degrade for
 *     larger removals.
 *   - the lamellar elastic rules: I expect k=1 and k=2 to heal at or near 1.0,
 *     and k=5 (particle-sized) to heal FAR less often, because a five-cell hole in
 *     a four-cell-thick slab removes more than a full cross-section of the slab.
 *   If that holds, the honest consequence is: these media are elastic to point
 *   defects but NOT to particle impacts, so a particle that touches the screen
 *   still leaves permanent damage. Elasticity of a point defect is not enough.
 *   If instead k=5 also heals near 1.0, we have a screen that can be hit by a
 *   particle-sized object and recover -- and that is the missing piece of the
 *   whole programme.
 * ---------------------------------------------------------------------------
 */

const W = 24, H = 24, NS = W * H;
const T = 100;

function stepR(g, bm, sm) {
  const out = new Uint8Array(NS);
  for (let y = 0; y < H; y++) {
    const yp = ((y - 1 + H) % H) * W, ym = ((y + 1) % H) * W, yc = y * W;
    for (let x = 0; x < W; x++) {
      const xp = (x - 1 + W) % W, xm = (x + 1) % W;
      const n = g[yp + xp] + g[yp + x] + g[yp + xm] + g[yc + xp] + g[yc + xm] + g[ym + xp] + g[ym + x] + g[ym + xm];
      const i = yc + x;
      out[i] = g[i] ? ((sm >> n) & 1) : ((bm >> n) & 1);
    }
  }
  return out;
}
function eq(a, b) { for (let i = 0; i < NS; i++) if (a[i] !== b[i]) return false; return true; }
function density(g) { let s = 0; for (let i = 0; i < NS; i++) s += g[i]; return s / NS; }
function mulberry32(a) { return function () { a |= 0; a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }

function parseRule(s) {
  const m = /^B\{?([0-8]*)\}?\/S\{?([0-8]*)\}?$/.exec(s.trim());
  let bm = 0, sm = 0;
  for (const c of m[1]) bm |= 1 << Number(c);
  for (const c of m[2]) sm |= 1 << Number(c);
  return { bm, sm };
}

// converge to the background cycle. The seed is deterministic, and deliberately the
// same kind of seed that crystalline_solid_search.js used (LCG at density 0.5 plus a
// 2x2 block at the origin), because these lamellar rules do NOT reach their cycle
// from an arbitrary random soup within a reasonable time.
function converge(r, seedTag) {
  let s = (987654321 + seedTag * 715827883) % 2147483647;
  let g = new Uint8Array(NS);
  for (let i = 0; i < NS; i++) { s = (s * 1103515245 + 12345) & 0x7fffffff; g[i] = ((s >> 16) & 1) ? 1 : 0; }
  g[0] = 1; g[1] = 1; g[W] = 1; g[W + 1] = 1;
  const hist = [];
  for (let t = 0; t < 800; t++) {
    hist.push(g);
    g = stepR(g, r.bm, r.sm);
    for (let p = 1; p <= 8; p++) {
      if (hist.length >= p && eq(g, hist[hist.length - p])) return { cycle: hist.slice(hist.length - p), period: p };
    }
  }
  return null;
}

// all connected k-cell blobs inside a given set of allowed cells, sampled
function blobsOfSize(allowed, k, limit, rnd) {
  const idx = [...allowed];
  const out = [];
  const seen = new Set();
  const key = (arr) => arr.slice().sort((a, b) => a - b).join(',');
  let guard = 0;
  while (out.length < limit && guard++ < limit * 40) {
    // grow a connected blob from a random start
    const start = idx[Math.floor(rnd() * idx.length)];
    const blob = [start];
    const inBlob = new Set(blob);
    while (blob.length < k) {
      const cand = [];
      for (const c of blob) {
        const x = c % W, y = (c - x) / W;
        for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
          const nk = (((y + dy + H) % H) * W) + ((x + dx + W) % W);
          if (!inBlob.has(nk) && allowed.has(nk)) cand.push(nk);
        }
      }
      if (!cand.length) break;
      const pick = cand[Math.floor(rnd() * cand.length)];
      blob.push(pick); inBlob.add(pick);
    }
    if (blob.length !== k) continue;
    const kk = key(blob);
    if (seen.has(kk)) continue;
    seen.add(kk); out.push(blob);
  }
  return out;
}

function measureRule(name, ruleStr, seedTag, fixedBg, truePeriod) {
  const r = parseRule(ruleStr);
  const bg = fixedBg ? { cycle: [fixedBg], period: truePeriod || 1 } : converge(r, seedTag);
  if (!bg) { console.log(`  ${name}: no background converged`); return; }
  const cycle = bg.cycle;
  const bgNext = stepR(cycle[0], r.bm, r.sm);
  const cycEq = (g) => cycle.some(c => eq(g, c));
  const heals = (start) => {
    let g = start;
    for (let t = 0; t < T; t++) { if (cycEq(g)) return true; g = stepR(g, r.bm, r.sm); }
    return cycEq(g);
  };
  const live = new Set(), dead = new Set();
  for (let i = 0; i < NS; i++) (cycle[0][i] ? live : dead).add(i);
  const rnd = mulberry32(4242);

  console.log(`  ${name}  (${ruleStr})  period ${bg.period}, density ${density(cycle[0]).toFixed(4)}, live ${live.size}`);
  console.log('    k    direction   relevant   inert   healed   heal rate');
  for (const k of [1, 2, 3, 5]) {
    for (const dir of ['remove', 'add']) {
      const pool = dir === 'remove' ? live : dead;
      const limit = k === 1 ? 400 : k === 2 ? 250 : k === 3 ? 200 : 150;
      const blobs = k === 1 ? [...pool].slice(0, limit).map(i => [i]) : blobsOfSize(pool, k, limit, rnd);
      let healed = 0, inert = 0, rel = 0;
      for (const blob of blobs) {
        const g = Uint8Array.from(cycle[0]);
        for (const i of blob) g[i] = dir === 'remove' ? 0 : 1;
        // additions must pass the causal filter: a cell that is about to be born
        // anyway is not a perturbation at all, and counting it inflates the rate
        if (dir === 'add' && eq(stepR(g, r.bm, r.sm), bgNext)) { inert++; continue; }
        rel++;
        if (heals(g)) healed++;
      }
      const rate = rel ? healed / rel : NaN;
      console.log(`    ${String(k).padStart(2)}   ${dir.padEnd(10)}  ${String(rel).padStart(8)}   ${String(inert).padStart(5)}   ${String(healed).padStart(6)}   ${Number.isNaN(rate) ? 'n/a' : rate.toFixed(4)}`);
    }
  }
  console.log('');
}

console.log('=== CONTROL: Life block lattice, constructed directly (must show: additions never heal) ===');
{
  const lat = new Uint8Array(NS);
  for (let j = 0; j < H; j += 4) for (let i = 0; i < W; i += 4) {
    lat[j * W + i] = 1; lat[j * W + i + 1] = 1; lat[(j + 1) * W + i] = 1; lat[(j + 1) * W + i + 1] = 1;
  }
  measureRule('Life block lattice (spacing 4)', 'B3/S23', 0, lat);
}

// The five lamellar backgrounds transcribed from the pictures printed by
// crystalline_solid_search.js. Transcribing a machine-found background is only safe
// if it is then verified to be a genuine member of a short cycle -- so every one of
// them is checked (period and density) before it is used, and a wrong transcription
// shows up immediately as a failed check instead of silently poisoning the table.
function makeBg(pred) {
  const g = new Uint8Array(NS);
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) g[y * W + x] = pred(x, y) ? 1 : 0;
  return g;
}
const CAND = [
  ['lamellar, bands 4 thick / 2 gap', 'B{023456}/S{01245}', (x, y) => (y % 6) >= 1 && (y % 6) <= 4, 0.6667],
  ['lamellar, stripes', 'B{023456}/S{012456}', (x, y) => (x % 6) >= 1 && (x % 6) <= 4, 0.6667],
  ['lamellar, sparse stripes', 'B{06}/S{28}', (x, y) => (x % 6) === 0 || (x % 6) === 5, 0.3333],
  ['lamellar, sparse stripes +8', 'B{068}/S{28}', (x, y) => (x % 6) === 0 || (x % 6) === 5, 0.3333],
  ['single slab 4 thick', 'B{01478}/S{01}', (x, y) => (y % 12) >= 1 && (y % 12) <= 4, 0.3333],
];
console.log('=== verification that each transcribed background is a genuine cycle member ===');
const verified = [];
for (const [name, ruleStr, pred, expDens] of CAND) {
  const r = parseRule(ruleStr);
  const bg = makeBg(pred);
  const d = density(bg);
  const s1 = stepR(bg, r.bm, r.sm);
  const s2 = stepR(s1, r.bm, r.sm);
  let period = 0;
  if (eq(s1, bg)) period = 1; else if (eq(s2, bg)) period = 2;
  const okD = Math.abs(d - expDens) < 1e-9;
  const okP = period > 0;
  console.log(`  ${name.padEnd(32)} ${ruleStr.padEnd(20)} density ${d.toFixed(4)} ${okD ? '(matches the search)' : `*** expected ${expDens} ***`}  period ${period || 'NOT ON A CYCLE'}  ${okP ? 'PASS' : '*** FAIL ***'}`);
  if (okP) verified.push([name, ruleStr, bg, period]);
  else console.log(`     -> skipped, the transcription is wrong`);
}
console.log('  (period 2 means the background blinks: its state alternates between two phases)');
console.log('');
console.log('=== blob impact test on the 5 verified two-sidedly elastic crystalline rules ===');
for (const [name, ruleStr, bg, period] of verified) measureRule(name, ruleStr, 0, bg, period);
console.log('  Read the table as: does the medium repair a hole/bump of k cells?');
console.log('  A particle (glider = 5 cells, LWSS = 9-12) is k >= 5. If k=5 heal rates collapse,');
console.log('  these media repair scratches but not impacts, and a particle still writes a');
console.log('  permanent record -- which is exactly what we were trying to escape.');
