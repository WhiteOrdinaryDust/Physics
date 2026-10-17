'use strict';
/*
 * additivity_test.js
 * ---------------------------------------------------------------------------
 * STEP 7: is the cascade amplification ADDITIVE?
 *
 * Why this decides whether the "reactivity spectrum" is a real law or an artifact:
 *
 *   The previous round found that amplification A takes discrete values
 *   (5.6, 2.3, 13.2, 31.9, ...) depending sharply on impact geometry. If two
 *   independent cascades simply ADD (combined peak = sum of peaks), then A
 *   behaves like an extensive conserved-ish quantity -- something that could
 *   act as a genuine constraint on the particle sector. If instead they
 *   interfere, quench, or saturate, then A is just a per-collision geometric
 *   accident with no compositional law, and it cannot carry the weight we want.
 *
 * METHOD (the control matters as much as the test):
 *   1. measure each island ALONE                 -> peak_i
 *   2. run two islands SEPARATED far apart       -> peak_sep
 *      - this must be measured first, because if the torus is too small the two
 *        cascades collide with each other and the comparison is meaningless
 *      - prediction if additive: peak_sep ~ peak_1 + peak_2 - 10 (each island
 *        needs its own 10-cell seed, and the seeds are counted once)
 *   3. run the two islands OVERLAPPING            -> peak_over
 *      - the interesting case: do they add, quench, or interfere?
 * ---------------------------------------------------------------------------
 */

const SIZE = 600;
function wrapCell(x, y) { return ((y % SIZE) + SIZE) % SIZE * SIZE + (((x % SIZE) + SIZE) % SIZE); }

function step(set) {
  const counts = new Map();
  for (const k of set) {
    const x = k % SIZE, y = (k - x) / SIZE;
    for (let dy = -1; dy <= 1; dy++) {
      for (let dx = -1; dx <= 1; dx++) {
        if (dx === 0 && dy === 0) continue;
        const nk = wrapCell(x + dx, y + dy);
        counts.set(nk, (counts.get(nk) || 0) + 1);
      }
    }
  }
  const out = new Set();
  for (const [k, c] of counts) if (c === 3 || (c === 2 && set.has(k))) out.add(k);
  return out;
}

const GLIDER0 = [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]];
function rotCells(turns) {
  return GLIDER0.map(([dx, dy]) => {
    let a = dx, b = dy;
    for (let i = 0; i < turns; i++) { const na = -b, nb = a; a = na; b = nb; }
    return [a, b];
  });
}
function gliderAt(x, y, turns) {
  const s = new Set();
  for (const [dx, dy] of rotCells(turns)) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }

/*
 * One island = the two-glider pair (heading, ox, oy) placed at (bx, by).
 * Returns the pair of glider sets.
 */
function island(bx, by, heading, ox, oy) {
  return [gliderAt(bx, by, 0), gliderAt(bx + ox, by + oy, heading)];
}

function measure(cellsList, steps) {
  let s = cellsList.length === 1 ? cellsList[0] : union(...cellsList);
  const initMass = s.size;
  let peak = initMass, peakT = 0;
  for (let t = 1; t <= steps; t++) {
    s = step(s);
    if (s.size > peak) { peak = s.size; peakT = t; }
    if (s.size === 0) break;
  }
  return { initMass, peak, peakT, final: s.size, A: peak / initMass };
}

const STEPS = 2500;

// ---- candidates from the reactivity sweep, with their known A
const ISLANDS = [
  { name: 'h1(6,1)', heading: 1, ox: 6, oy: 1, A_known: 31.9 },
  { name: 'h2(0,1)', heading: 2, ox: 0, oy: 1, A_known: 22.8 },
  { name: 'h2(4,1)', heading: 2, ox: 4, oy: 1, A_known: 14.0 },
  { name: 'h3(1,0)', heading: 3, ox: 1, oy: 0, A_known: 12.8 },
  { name: 'h1(3,5)', heading: 1, ox: 3, oy: 5, A_known: 13.2 },
];

console.log('=== additivity test for cascade amplification ===');
console.log(`torus ${SIZE}x${SIZE}, ${STEPS} steps\n`);

console.log('--- step 1: each island ALONE (baseline) ---');
const base = new Map();
for (const isl of ISLANDS) {
  const pair = island(150, 150, isl.heading, isl.ox, isl.oy);
  const r = measure(pair, STEPS);
  base.set(isl.name, r);
  console.log(`  ${isl.name.padEnd(10)} init=${r.initMass}  peak=${String(r.peak).padStart(4)}@t${String(r.peakT).padStart(5)}  final=${String(r.final).padStart(4)}  A=${r.A.toFixed(1)}  (sweep said ${isl.A_known})`);
}

// ---- step 2: pairs of islands, SEPARATED far apart
console.log('\n--- step 2: two islands SEPARATED (must not interact) ---');
const separation = 140;
console.log(`  islands placed ${separation} cells apart; a glider travels at most ${Math.floor(STEPS / 4)} cells in ${STEPS} steps`);
const SEP_PAIRS = [
  ['h1(6,1)', 'h2(0,1)'],
  ['h1(6,1)', 'h3(1,0)'],
  ['h2(4,1)', 'h3(1,0)'],
  ['h2(0,1)', 'h2(4,1)'],
];
const sepResults = [];
for (const [n1, n2] of SEP_PAIRS) {
  const i1 = ISLANDS.find(z => z.name === n1), i2 = ISLANDS.find(z => z.name === n2);
  const p1 = island(120, 120, i1.heading, i1.ox, i1.oy);
  const p2 = island(120 + separation, 120 + separation, i2.heading, i2.ox, i2.oy);
  const r = measure([...p1, ...p2], STEPS);
  const b1 = base.get(n1), b2 = base.get(n2);
  const predicted = b1.peak + b2.peak - 20;   // seeds counted twice
  sepResults.push({ n1, n2, r, predicted });
  console.log(`  ${n1} + ${n2}:  init=${r.initMass} peak=${r.peak} final=${r.final} A=${r.A.toFixed(1)}`);
  console.log(`      additive prediction (peak1+peak2-20) = ${predicted}   ratio observed/predicted = ${(r.peak / predicted).toFixed(3)}`);
}

// ---- step 3: same pairs but OVERLAPPING (collide with each other)
console.log('\n--- step 3: the same island pairs OVERLAPPED (cascades interact) ---');
for (const [n1, n2] of SEP_PAIRS) {
  const i1 = ISLANDS.find(z => z.name === n1), i2 = ISLANDS.find(z => z.name === n2);
  const p1 = island(150, 150, i1.heading, i1.ox, i1.oy);
  // place the second island's seed right on top of the first collision site
  const p2 = island(150, 150, i2.heading, i2.ox, i2.oy);
  const r = measure([...p1, ...p2], STEPS);
  const b1 = base.get(n1), b2 = base.get(n2);
  const predicted = b1.peak + b2.peak - 20;
  console.log(`  ${n1} + ${n2}:  init=${r.initMass} peak=${r.peak} final=${r.final} A=${r.A.toFixed(1)}  (sep. would predict ${predicted})`);
}

// ---- verdict
console.log('\n=== verdict ===');
const ratios = sepResults.map(x => x.r.peak / x.predicted);
const meanR = ratios.reduce((a, b) => a + b, 0) / ratios.length;
const spread = Math.max(...ratios) - Math.min(...ratios);
console.log(`  separated: observed/predicted ratios = ${ratios.map(r => r.toFixed(3)).join(', ')}`);
console.log(`  mean ${meanR.toFixed(3)}, spread ${spread.toFixed(3)}`);
console.log('');
if (spread < 0.15) {
  console.log('  => ADDITIVE. Independent cascades combine by simple addition.');
  console.log('     Amplification behaves as an extensive quantity, which is what a');
  console.log('     compositional law would require.');
} else {
  console.log('  => NOT cleanly additive: the separated cascades do not combine as');
  console.log('     peak1 + peak2. Either they still interact, or the quantity is not');
  console.log('     extensive. Either way it cannot be treated as a conserved charge');
  console.log('     without further work.');
}
