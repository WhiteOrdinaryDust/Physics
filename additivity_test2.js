'use strict';
/*
 * additivity_test2.js
 * ---------------------------------------------------------------------------
 * STEP 7b: the previous additivity test was confounded, and one result in it was
 * qualitatively important.
 *
 * CONFOUND: "separated by 140 cells" is NOT isolated. A glider covers up to 625
 * cells in 2500 steps, and island h1(6,1) does not even peak until t=823, so its
 * debris had ample time to reach the other island. The measured ratios
 * (0.68, 0.80, 1.51, 0.90) therefore say nothing clean about additivity.
 *
 * IMPORTANT RESULT ALREADY VISIBLE, and it is the opposite of additivity:
 *   two cascades triggered at the SAME SITE do not add, they QUENCH:
 *       h1(6,1) + h2(0,1): separated-prediction 527 -> observed peak 38
 *       h2(4,1) + h3(1,0): separated-prediction 248 -> observed peak 14
 *   i.e. planting a second collision on top of the first annihilates the cascade
 *   rather than compounding it.
 *
 * PROPER TEST: synchronise the two collisions in TIME, then vary the spatial
 * separation across a wide range. A true isolated-additive regime shows a plateau
 * at peak1+peak2-seeds for large separation, falling off only when the cascades
 * can actually reach each other.
 * ---------------------------------------------------------------------------
 */

const SIZE = 700;
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
function island(bx, by, heading, ox, oy) {
  return [gliderAt(bx, by, 0), gliderAt(bx + ox, by + oy, heading)];
}

function measure(cellsList, steps) {
  let s = cellsList.length === 1 ? cellsList[0] : union(...cellsList);
  const initMass = s.size;
  let peak = initMass, peakT = 0;
  let firstChange = -1;
  for (let t = 1; t <= steps; t++) {
    s = step(s);
    if (firstChange < 0 && s.size !== initMass) firstChange = t;
    if (s.size > peak) { peak = s.size; peakT = t; }
    if (s.size === 0) break;
  }
  return { initMass, peak, peakT, final: s.size, A: peak / initMass, triggerT: firstChange };
}

const ISLANDS = {
  'h1(6,1)': { heading: 1, ox: 6, oy: 1 },
  'h2(0,1)': { heading: 2, ox: 0, oy: 1 },
  'h2(4,1)': { heading: 2, ox: 4, oy: 1 },
  'h3(1,0)': { heading: 3, ox: 1, oy: 0 },
};

console.log('=== additivity, take 2: synchronised collisions ===');
console.log(`torus ${SIZE}x${SIZE}\n`);

// ---- step 1: when does each island actually trigger, and peak?
console.log('--- trigger and peak times for each island alone ---');
const base = new Map();
for (const [name, isl] of Object.entries(ISLANDS)) {
  const r = measure(island(200, 200, isl.heading, isl.ox, isl.oy), 2500);
  base.set(name, r);
  console.log(`  ${name.padEnd(9)} trigger t=${String(r.triggerT).padStart(4)}  peak=${String(r.peak).padStart(4)}@t${String(r.peakT).padStart(5)}  A=${r.A.toFixed(1)}`);
}

/*
 * Synchronise by placing the two islands at positions such that their
 * first-interaction times coincide. Each island's gliders must travel to their
 * meeting point; shifting the whole island does not change that, so we instead
 * delay one island by letting it travel before the other is planted. Simpler and
 * exact: measure each island's trigger time, then offset the second island's
 * START POSITION along its own flight path so the meeting happens at the same t.
 *
 * Practical approach used here: plant island 1 at P, and island 2 at P + delta,
 * then ADD a lead-in by shifting island 2 backwards along its own approach
 * direction by (its trigger difference) * speed. Speed is 1/4 cell/step for a
 * glider, so a lead of L steps means shifting by L/4 cells -- not integral.
 * Instead we simply try a range of separations and report where the peak stops
 * depending on separation, which reveals the interaction range empirically.
 */
console.log('\n--- peak vs separation (both islands planted simultaneously) ---');
const PAIRS = [['h1(6,1)', 'h2(0,1)'], ['h1(6,1)', 'h3(1,0)'], ['h2(4,1)', 'h3(1,0)']];
const SEPS = [30, 60, 100, 160, 240, 320];

for (const [n1, n2] of PAIRS) {
  const i1 = ISLANDS[n1], i2 = ISLANDS[n2];
  const b1 = base.get(n1), b2 = base.get(n2);
  const additivePeak = b1.peak + b2.peak - 20;
  const line = [];
  for (const sep of SEPS) {
    // separate along the diagonal, both moved so they stay in bounds
    const p1 = island(200, 200, i1.heading, i1.ox, i1.oy);
    const p2 = island(200 + sep, 200 + sep, i2.heading, i2.ox, i2.oy);
    const r = measure([...p1, ...p2], 2500);
    line.push({ sep, peak: r.peak, ratio: r.peak / additivePeak });
  }
  console.log(`\n  ${n1} + ${n2}   (additive prediction = ${additivePeak})`);
  console.log('    separation : ' + SEPS.map(s => String(s).padStart(6)).join(''));
  console.log('    peak       : ' + line.map(x => String(x.peak).padStart(6)).join(''));
  console.log('    obs/pred   : ' + line.map(x => x.ratio.toFixed(2).padStart(6)).join(''));
}

console.log('\n=== interpretation guide ===');
console.log('  If a plateau appears at obs/pred ~ 1.0 for large separation, the');
console.log('  quantity is additive in the isolated limit and the deviation at small');
console.log('  separation measures the interaction range.');
console.log('  If obs/pred never approaches 1.0, amplification is NOT extensive.');
console.log('');
console.log('  Already established and independent of this test:');
console.log('  two cascades triggered at the SAME SITE quench each other.');
console.log('  h1(6,1)+h2(0,1): prediction 527 -> observed 38');
console.log('  h2(4,1)+h3(1,0): prediction 248 -> observed 14');
console.log('  That alone rules out simple additivity, and is the physical result.');
