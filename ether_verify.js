'use strict';
/*
 * ether_verify.js
 * ---------------------------------------------------------------------------
 * Verify the rule-110 ether's actual time period and phase velocity.
 *
 * The timescale scan measured a dominant temporal period of 7 for the ether,
 * while the ether's spatial period is 14 and I had expected a time period of 4
 * (that is the phase period of the GLIDER). I could not reconcile 7 by
 * reasoning, so this run measures it directly instead of arguing about it:
 *
 *   1. evolve the ether from a clean initial condition
 *   2. extract its 14-cell repeating pattern at each step
 *   3. find the smallest T with pattern(t+T) == pattern(t)      -> time period
 *   4. record the spatial shift over that T                     -> phase shift
 *   5. phase velocity = shift / T
 *
 * Also measure the glider's speed in the same run for comparison, and report the
 * ratio (the "medium speed vs matter speed" quantity).
 * ---------------------------------------------------------------------------
 */

const W = 560;                    // multiple of 14, room for a long run
const TMAX = 200;
const ETHER = [0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0];
const SPATIAL = ETHER.length;     // 14

function step110(row, w) {
  const out = new Uint8Array(w);
  for (let i = 0; i < w; i++) {
    const idx = (row[(i - 1 + w) % w] << 2) | (row[i] << 1) | row[(i + 1) % w];
    out[i] = (110 >> idx) & 1;
  }
  return out;
}

// the 14-cell pattern starting at position p
function tile(row, p, len) {
  let s = '';
  for (let i = 0; i < len; i++) s += row[(p + i) % row.length];
  return s;
}
function minRotation(s) {
  let best = s;
  for (let i = 1; i < s.length; i++) {
    const r = s.slice(i) + s.slice(0, i);
    if (r < best) best = r;
  }
  return best;
}

console.log('=== rule 110 ether: measured period and phase velocity ===\n');

// ---- 1. clean ether
let row = new Uint8Array(W);
for (let i = 0; i < W; i++) row[i] = ETHER[i % SPATIAL];

const canon = [];
for (let t = 0; t <= TMAX; t++) {
  canon.push(minRotation(tile(row, 0, SPATIAL)));
  row = step110(row, W);
}
console.log('first 10 canonical 14-cell patterns:');
for (let t = 0; t < 10; t++) console.log(`  t=${String(t).padStart(3)}  ${canon[t]}`);

// ---- 2. time period of the pattern
let period = -1;
for (let T = 1; T <= 60; T++) {
  let ok = true;
  for (let t = 0; t + T <= TMAX; t++) if (canon[t] !== canon[t + T]) { ok = false; break; }
  if (ok) { period = T; break; }
}
console.log(`\ncanonical pattern time period: ${period === -1 ? 'not found within 60' : period}`);

// ---- 3. phase shift over one period (position of the pattern start)
function findPhase(row, pat) {
  for (let p = 0; p < row.length; p++) if (tile(row, p, SPATIAL) === pat) return p;
  return -1;
}
// rebuild from t=0
row = new Uint8Array(W);
for (let i = 0; i < W; i++) row[i] = ETHER[i % SPATIAL];
const snaps = [];
for (let t = 0; t <= TMAX; t++) { snaps.push(Uint8Array.from(row)); row = step110(row, W); }

if (period > 0) {
  const p0 = findPhase(snaps[0], canon[0]);
  const pT = findPhase(snaps[period], canon[0]);
  let shift = pT - p0;
  console.log(`\npattern start at t=0    : ${p0}`);
  console.log(`pattern start at t=${String(period).padStart(2)}   : ${pT}`);
  console.log(`spatial shift over T    : ${shift} cells`);
  console.log(`phase velocity          : ${(shift / period).toFixed(4)} cells/step`);
  console.log(`substrate max speed     : 1 cell/step (radius-1 rule)`);
}

// ---- 4. glider speed in the same substrate, for the ratio
console.log('\n=== glider speed in the same substrate ===');
{
  // seed one glider onto the ether and follow its bounding box
  let r2 = new Uint8Array(W);
  for (let i = 0; i < W; i++) r2[i] = ETHER[i % SPATIAL];
  const g = [1, 1, 1, 0, 1, 1, 1];
  for (let j = 0; j < g.length; j++) r2[(200 + j) % W] ^= g[j];
  // find the defect: cells differing from the pure ether phase
  const trace = [];
  for (let t = 0; t <= 120; t++) {
    // locate the largest block of cells that disagree with the ether of some phase
    let bestP = 0, bestScore = -1;
    for (let p = 0; p < SPATIAL; p++) {
      let sc = 0;
      for (let i = 0; i < W; i++) if (r2[i] === ETHER[(i + p) % SPATIAL]) sc++;
      if (sc > bestScore) { bestScore = sc; bestP = p; }
    }
    // defect positions relative to that best phase
    const bad = [];
    for (let i = 0; i < W; i++) if (r2[i] !== ETHER[(i + bestP) % SPATIAL]) bad.push(i);
    const centroid = bad.length ? bad.reduce((a, b) => a + b, 0) / bad.length : -1;
    trace.push({ t, n: bad.length, centroid, phase: bestP });
    r2 = step110(r2, W);
  }
  const clean = trace.filter(x => x.n > 0 && x.n < 40);
  console.log('  selected samples (defect position over time):');
  for (const x of clean.filter((_, i) => i % 10 === 0).slice(0, 10)) {
    console.log(`    t=${String(x.t).padStart(3)}  defect cells=${String(x.n).padStart(3)}  centroid=${x.centroid.toFixed(1)}`);
  }
  // speed from first and last clean samples
  if (clean.length >= 2) {
    const a = clean[0], b = clean[clean.length - 1];
    let dc = b.centroid - a.centroid;
    if (dc > W / 2) dc -= W; if (dc < -W / 2) dc += W;
    console.log(`\n  glider centroid moved ${dc.toFixed(1)} cells in ${b.t - a.t} steps`);
    console.log(`  glider speed = ${(dc / (b.t - a.t)).toFixed(4)} cells/step`);
  }
}

console.log('\n=== reading ===');
console.log('  The substrate clock is 1 step and its max signal speed is 1 cell/step.');
console.log('  If the ether pattern period really is 7 with a 7-cell shift, phase');
console.log('  velocity is 1 cell/step -- the same as the substrate maximum, which would');
console.log('  mean the ether is itself a ballistic structure, not a slower medium.');
console.log('  A ratio c1/v close to 1 would mean the two-layer separation is ABSENT.');
