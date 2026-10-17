'use strict';
/*
 * saturation_test.js
 * ---------------------------------------------------------------------------
 * STEP 8: does the cascade size SATURATE with the number of seed islands?
 *
 * Motivation: in the previous test, two cascades did not add, and -- more telling
 * -- the combined peak barely changed when the second island was moved from 60 to
 * 320 cells away. A distant second cascade made no difference at all.
 *
 * That is the signature of SATURATION rather than additivity. This test measures
 * it directly: plant n identical cascade islands, well separated, and plot peak
 * mass against n.
 *
 *   peak ~ n            -> extensive (additive) quantity
 *   peak ~ constant     -> saturating; the cascade size is capped by something
 *                          internal to the dynamics, not by the number of seeds
 * ---------------------------------------------------------------------------
 */

const SIZE = 900;               // large enough to hold several separate cascades
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
  let s = union(...cellsList);
  const initMass = s.size;
  let peak = initMass, peakT = 0;
  const trace = [];
  for (let t = 1; t <= steps; t++) {
    s = step(s);
    if (s.size > peak) { peak = s.size; peakT = t; }
    if (t % 250 === 0) trace.push(s.size);
    if (s.size === 0) break;
  }
  return { initMass, peak, peakT, final: s.size, trace };
}

const SPEC = { heading: 1, ox: 6, oy: 1 };      // h1(6,1), the strongest cascade
const STEPS = 3000;
const SPACING = 200;                             // far apart; glider covers <= 750 cells

console.log('=== saturation test: n identical cascade islands ===');
console.log(`torus ${SIZE}x${SIZE}, island = h1(6,1), spacing ${SPACING} cells, ${STEPS} steps\n`);
console.log('   n   initMass   peak@t      final   peak/init   mass at t=500,1000,...,3000');

const rows = [];
for (let n = 1; n <= 8; n++) {
  // place islands on a grid so they are pairwise separated
  const cells = [];
  for (let i = 0; i < n; i++) {
    const col = i % 4, row = Math.floor(i / 4);
    const bx = 100 + col * SPACING, by = 100 + row * SPACING;
    cells.push(...island(bx, by, SPEC.heading, SPEC.ox, SPEC.oy));
  }
  const r = measure(cells, STEPS);
  rows.push({ n, ...r });
  console.log(`  ${String(n).padStart(3)}  ${String(r.initMass).padStart(8)}   ${String(r.peak).padStart(4)}@${String(r.peakT).padStart(4)}  ${String(r.final).padStart(8)}   ${(r.peak / r.initMass).toFixed(1).padStart(8)}   ${r.trace.join(', ')}`);
}

// ---- linearity check
console.log('\n=== linearity ===');
const p1 = rows[0].peak;
console.log(`  single-island peak = ${p1}`);
console.log('   n   observed peak   n x single   ratio');
for (const r of rows) {
  const pred = r.n * p1;
  console.log(`  ${String(r.n).padStart(3)}   ${String(r.peak).padStart(11)}   ${String(pred).padStart(10)}   ${(r.peak / pred).toFixed(3)}`);
}

const ratios = rows.map(r => r.peak / (r.n * p1));
const lastRatios = ratios.slice(3);
const meanLast = lastRatios.reduce((a, b) => a + b, 0) / lastRatios.length;

console.log('');
if (meanLast > 0.85) {
  console.log('  => ROUGHLY EXTENSIVE: peak grows with n, so amplification adds across');
  console.log('     well-separated seeds. It can be treated as an extensive quantity.');
} else if (meanLast < 0.5) {
  console.log('  => SATURATES: peak does NOT grow proportionally with seed count.');
  console.log('     The total cascade size is capped by the dynamics itself, not by how');
  console.log('     much seed material is supplied. There is a maximum eruption size.');
} else {
  console.log(`  => SUB-EXTENSIVE: peak grows slower than n (ratio ~${meanLast.toFixed(2)}).`);
  console.log('     Partial saturation: more seeds help, but with diminishing returns.');
}
