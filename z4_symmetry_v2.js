'use strict';
/*
 * z4_symmetry_v2.js
 * ---------------------------------------------------------------------------
 * Rotational symmetry of Life patterns -- done correctly this time.
 *
 * WHAT THE PREVIOUS VERSION GOT WRONG
 *   It rotated the pattern and re-normalised at every step, so the coordinate
 *   origin drifted with each rotation and unrelated patterns were compared. It
 *   reported a 3x1 horizontal line as four-fold symmetric, which is impossible.
 *
 * THE FIX
 *   Rotation is ALWAYS about one FIXED candidate centre. For a box of size n x n
 *   the possible centres are the integer points and the half-integer points, so
 *   there are only a handful to test. A pattern has Z4 symmetry iff there is a
 *   centre about which rot90 maps the cell set onto itself.
 *
 * SELF-CHECKS (must pass before any result is trusted)
 *   - a 2x2 block            -> Z4 (rotating it by 90 degrees is itself)
 *   - a 3x1 line             -> NOT Z4
 *   - a random 5x5 pattern   -> NOT Z4 (overwhelmingly)
 *   - a 3x3 ring (8 cells)   -> Z4
 *
 * THEN THE ACTUAL QUESTION
 *   Which Z4-symmetric patterns SURVIVE, and importantly: are any of them
 *   MOVERS or OSCILLATORS WITH A NONTRIVIAL MODE SPACE, rather than static
 *   still lifes? A static still life is useless for the complex-structure idea,
 *   because it has no internal clock.
 * ---------------------------------------------------------------------------
 */

const SIZE = 60;
function wrapCell(x, y) { return ((y % SIZE) + SIZE) % SIZE * SIZE + (((x % SIZE) + SIZE) % SIZE); }

function lifeStep(set) {
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

// ---- rotation about a FIXED centre c, working in doubled coordinates so that
// half-integer centres are representable exactly.
// A point (x,y) becomes (2x, 2y); the centre is (cx2, cy2) in the same units.
function rot90About(pt, c) {
  // translate so centre is origin, rotate, translate back
  const dx = pt[0] - c[0], dy = pt[1] - c[1];
  const rx = -dy, ry = dx;                    // 90 degrees
  return [c[0] + rx, c[1] + ry];
}
function setKey(pts) { return pts.map(p => p.join(',')).sort().join(';'); }

// is the cell set (given as [x,y] pairs) invariant under rot90 about centre c?
function invariantUnderRot(pts, c) {
  const s = new Set(pts.map(p => p.join(',')));
  for (const p of pts) {
    const q = rot90About(p, c);
    if (!s.has(q.join(','))) return false;
  }
  return true;
}

// candidate centres: integer and half-integer points inside the bounding box
function findSymmetryCentres(pts) {
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  for (const [x, y] of pts) {
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
  }
  const found = { z4: [], z2: [] };
  for (let cx = 2 * minx; cx <= 2 * (maxx + 1); cx++) {
    for (let cy = 2 * miny; cy <= 2 * (maxy + 1); cy++) {
      const c = [cx, cy];
      if (invariantUnderRot(pts, c)) {
        // how many 90-degree turns return the pattern to itself?
        // Z4 if one turn works, else only Z2 (180 degrees)
        found.z4.push(c);
      } else {
        // test 180 degrees
        const p2 = pts.map(p => { const q = rot90About(p, c); return rot90About(q, c); });
        const s = new Set(pts.map(p => p.join(',')));
        if (p2.every(q => s.has(q.join(',')))) found.z2.push(c);
      }
    }
  }
  return found;
}

// ============================================================ SELF-CHECKS
console.log('=== SELF-CHECKS (must pass before any result is used) ===\n');
{
  const checks = [
    ['2x2 block', [[0, 0], [1, 0], [0, 1], [1, 1]], true],
    ['3x1 line', [[0, 0], [1, 0], [2, 0]], false],
    ['3x3 ring (8 cells)', [[0, 0], [1, 0], [2, 0], [0, 1], [2, 1], [0, 2], [1, 2], [2, 2]], true],
    ['single cell', [[0, 0]], true],
    ['2x1 line', [[0, 0], [1, 0]], false],
    ['L tromino', [[0, 0], [0, 1], [1, 1]], false],
    ['plus shape', [[1, 0], [0, 1], [1, 1], [2, 1], [1, 2]], true],
  ];
  let allPass = true;
  for (const [name, pts, expect] of checks) {
    const sym = findSymmetryCentres(pts);
    const got = sym.z4.length > 0;
    const pass = got === expect;
    if (!pass) allPass = false;
    console.log(`  ${name.padEnd(22)} Z4 = ${String(got).padEnd(5)} expected ${String(expect).padEnd(5)} ${pass ? 'PASS' : '*** FAIL ***'}`);
  }
  // random patterns must essentially never be Z4
  let rndZ4 = 0;
  let s = 12345 >>> 0;
  const rnd = () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
  for (let trial = 0; trial < 2000; trial++) {
    const pts = [];
    for (let x = 0; x < 5; x++) for (let y = 0; y < 5; y++) if (rnd() < 0.4) pts.push([x, y]);
    if (pts.length && findSymmetryCentres(pts).z4.length) rndZ4++;
  }
  console.log(`  random 5x5 patterns reported Z4: ${rndZ4} of 2000  ${rndZ4 <= 5 ? 'PASS (rare)' : '*** FAIL ***'}`);
  console.log('');
  if (!allPass) { console.log('*** a control failed; aborting before the search ***'); process.exit(1); }
}
console.log('  all controls passed\n');

// ============================================================ THE SEARCH
/*
 * Enumerate patterns inside a (2k+1) x (2k+1) box, but only those with Z4
 * symmetry. Enumerating all 2^(25) is wasteful, so we build the pattern from one
 * orbit under rot90 plus the cells on the rotation axes.
 */
function z4Patterns(k) {
  // box coordinates 0..2k
  const c = [2 * k, 2 * k];        // centre in doubled coordinates
  const seen = new Set();
  const out = [];
  const cells = [];
  for (let x = 0; x <= 2 * k; x++) for (let y = 0; y <= 2 * k; y++) cells.push([x, y]);
  // orbits under rotation about the centre
  const orbitOf = new Map();
  for (const p of cells) {
    const key = p.join(',');
    if (orbitOf.has(key)) continue;
    const orbit = [];
    let q = p;
    for (let i = 0; i < 4; i++) { orbit.push(q); q = rot90About(q, c); }
    const uniq = [];
    for (const o of orbit) if (!uniq.some(u => u[0] === o[0] && u[1] === o[1])) uniq.push(o);
    for (const o of uniq) orbitOf.set(o.join(','), uniq);
  }
  const reps = [];
  const usedRep = new Set();
  for (const p of cells) {
    const key = p.join(',');
    if (usedRep.has(key)) continue;
    const orb = orbitOf.get(key);
    for (const o of orb) usedRep.add(o.join(','));
    reps.push(orb);
  }
  // choose a subset of orbits: 2^reps.length (small for k<=2)
  const n = reps.length;
  if (n > 22) return { tooMany: n };
  for (let mask = 0; mask < (1 << n); mask++) {
    const pts = [];
    for (let i = 0; i < n; i++) if ((mask >> i) & 1) for (const p of reps[i]) pts.push(p);
    if (!pts.length) continue;
    const norm = setKey(pts);
    if (seen.has(norm)) continue;
    seen.add(norm);
    out.push(pts);
  }
  return { patterns: out, orbits: n };
}

console.log('=== Z4-symmetric patterns in a box, and whether they SURVIVE ===\n');

const results = [];
for (const k of [1, 2]) {          // 3x3 and 5x5 boxes
  const r = z4Patterns(k);
  if (r.tooMany) { console.log(`  ${2 * k + 1}x${2 * k + 1}: too many orbits (${r.tooMany})`); continue; }
  console.log(`--- ${2 * k + 1}x${2 * k + 1} box: ${r.patterns.length} distinct Z4-symmetric patterns (${r.orbits} orbits) ---`);
  let survived = 0;
  const survivors = [];
  for (const pts of r.patterns) {
    let set = new Set();
    for (const [x, y] of pts) set.add(wrapCell(25 + x, 25 + y));
    const initMass = set.size;
    // classify: still life (p=1, no move), oscillator (p>1, no move), mover, dead
    const seenStates = new Map();
    let cur = set, period = -1, pos0 = null;
    for (let t = 0; t <= 60; t++) {
      const sig = [...cur].sort((a, b) => a - b).join(',');
      if (seenStates.has(sig)) { period = t - seenStates.get(sig); break; }
      seenStates.set(sig, t);
      if (cur.size === 0) { period = -2; break; }
      cur = lifeStep(cur);
    }
    if (period < 0) continue;      // died or did not settle
    survived++;
    survivors.push({ pts, period, mass: initMass, k });
  }
  console.log(`    survived (periodic within 60 steps): ${survived}`);
  const byP = new Map();
  for (const s of survivors) byP.set(s.period, (byP.get(s.period) || 0) + 1);
  for (const [p, c] of [...byP.entries()].sort((a, b) => a[0] - b[0])) {
    console.log(`      period ${String(p).padStart(2)}: ${c}`);
  }
  results.push(...survivors);
  console.log('');
}

console.log('=== the interesting subclass: Z4 patterns with period > 1 ===');
const nontrivial = results.filter(s => s.period > 1);
console.log(`  count: ${nontrivial.length}`);
const shown = new Set();
for (const s of nontrivial.slice(0, 8)) {
  const key = setKey(s.pts);
  if (shown.has(key)) continue;
  shown.add(key);
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  for (const [x, y] of s.pts) { if (x < minx) minx = x; if (x > maxx) maxx = x; if (y < miny) miny = y; if (y > maxy) maxy = y; }
  const w = maxx - minx + 1, h = maxy - miny + 1;
  const g = Array.from({ length: h }, () => new Array(w).fill('.'));
  for (const [x, y] of s.pts) g[y - miny][x - minx] = 'o';
  console.log(`  ${w}x${h}, ${s.pts.length} cells, period ${s.period}`);
  for (const row of g) console.log('    ' + row.join(''));
}

console.log('\n=== what this establishes ===');
console.log('  A surviving Z4 pattern with period > 1 is the NECESSARY precondition for');
console.log('  the complex-structure idea. It is NOT sufficient: the next step would be to');
console.log('  show that the pattern\'s deformation modes close under the dynamics and that');
console.log('  the Z4 action on that mode space satisfies J^2 = -1.');
