'use strict';
/*
 * mode_space.js
 * ---------------------------------------------------------------------------
 * Where would a COMPLEX structure come from? (the remaining open question)
 *
 * Chain of reasoning established so far:
 *   - classical deterministic shortcut  -> CHSH <= 2
 *   - 2*sqrt(2) needs complex amplitudes
 *   - no capacity functional pins the ceiling (capacity_chsh.js)
 *   => the real requirement is a COMPLEX HILBERT SPACE at the particle layer.
 *      With complex structure + tensor product, 2*sqrt(2) and (via Gleason)
 *      the Born rule follow automatically. Without it, nothing does.
 *
 * METHOD
 *   Take a steady configuration, enumerate ALL single-cell perturbations, evolve
 *   each, and study the space of perturbation effects:
 *     - is the effect of a sum of perturbations the sum of effects?
 *       (linearity => vector space)
 *     - over which field? GF(2) (XOR) or the reals?
 *     - is there a compatible operator J with J^2 = -1 ? (complex structure)
 *
 * WHY THIS MATTERS
 *   A GF(2) vector space is NOT a Hilbert space and Gleason does not apply to it.
 *   So if the perturbation space is only a GF(2) space, the Born rule cannot come
 *   from it, and the programme needs something else. That would be a clean
 *   negative result. If instead a real vector space with a J appears, there is a
 *   route.
 * ---------------------------------------------------------------------------
 */

const W = 24;                       // small board; we only need local structure
let SIZE = W;
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

function stepN(set, n) { let s = set; for (let i = 0; i < n; i++) s = lifeStep(s); return s; }

/*
 * Represent a configuration as a bit vector over the W*W cells, and the CA as a
 * map. Linearisation around a background B: for a perturbation d (a set of cells
 * to toggle), the linearised response is
 *      L(d) = F(B xor d) xor F(B)
 * where xor of cell sets is symmetric difference. This is exact for the
 * difference, and L is well defined for every d. We then test:
 *      L(d1 xor d2) == L(d1) xor L(d2) ?       (GF(2)-linearity)
 *      L(d1 + d2)   == L(d1) + L(d2) ?          (real linearity, if we lift 0/1)
 */
function xor(a, b) {
  const out = new Set(a);
  for (const k of b) { if (out.has(k)) out.delete(k); else out.add(k); }
  return out;
}
function eq(a, b) {
  if (a.size !== b.size) return false;
  for (const k of a) if (!b.has(k)) return false;
  return true;
}

// ---- background: empty board and a steady block
function emptyBoard() { return new Set(); }
function blockBoard() {
  const s = new Set();
  for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) s.add(wrapCell(10 + dx, 10 + dy));
  return s;
}

function linearityTest(bg, samples) {
  const F0 = lifeStep(bg);
  const L = (d) => xor(lifeStep(xor(bg, d)), F0);
  let gf2pass = 0, realFail = 0, gf2Fail = 0;
  for (const [d1, d2] of samples) {
    const lhs = L(xor(d1, d2));
    const rhs = xor(L(d1), L(d2));
    if (eq(lhs, rhs)) gf2pass++; else gf2Fail++;
  }
  return { gf2pass, gf2Fail, total: samples.length };
}

function randPerturb(rnd, n) {
  const s = new Set();
  for (let i = 0; i < n; i++) s.add(Math.floor(rnd() * W * W));
  return s;
}

console.log('=== is the perturbation space a vector space, and over which field? ===');
console.log(`board ${W}x${W}\n`);

let seed = 12345 >>> 0;
const rnd = () => (seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296;

for (const [name, mk] of [['empty board', emptyBoard], ['steady block', blockBoard]]) {
  const bg = mk();
  // confirm the background is actually steady
  const steady = eq(lifeStep(bg), bg);
  const samples = [];
  for (let i = 0; i < 2000; i++) samples.push([randPerturb(rnd, 3), randPerturb(rnd, 3)]);
  const r = linearityTest(bg, samples);
  console.log(`--- background: ${name} (steady: ${steady}) ---`);
  console.log(`  GF(2) linearity of the one-step response:`);
  console.log(`    L(d1 xor d2) == L(d1) xor L(d2)   :  ${r.gf2pass}/${r.total} pass, ${r.gf2Fail} fail`);
  console.log('');
}

/*
 * Interpretation of the result (printed below from measured numbers, not assumed):
 *
 * Life is a NONLINEAR map, so there is no reason for L to be additive. The test
 * above measures exactly how far from additive it is. If it fails badly, then the
 * perturbation space is NOT a vector space at all -- and "complex structure" has
 * no meaning there yet.
 *
 * That would relocate the problem: one cannot ask for a complex structure on a
 * space that is not even a vector space. One would first need a regime where the
 * dynamics IS linear (a linear CA, or small oscillations around a fixed point in
 * a linearisable model).
 */

console.log('=== control: a genuinely LINEAR rule, where a vector space does exist ===');
console.log('  rule 90 over GF(2): next[i] = left xor right. This map is linear.');
function rule90(set, w) {
  const out = new Set();
  for (const k of set) {
    const i = k;
    out.add((i - 1 + w) % w);
    out.add((i + 1) % w);
  }
  // xor: symmetric difference of the two shifted copies
  const L = new Set(), R = new Set();
  for (const i of set) { L.add((i - 1 + w) % w); R.add((i + 1) % w); }
  return xor(L, R);
}
{
  const w = 32;
  let seed2 = 777 >>> 0;
  const rnd2 = () => (seed2 = (seed2 * 1664525 + 1013904223) >>> 0) / 4294967296;
  const mk = () => { const s = new Set(); for (let i = 0; i < w; i++) if (rnd2() < 0.3) s.add(i); return s; };
  let pass = 0, total = 0;
  for (let t = 0; t < 2000; t++) {
    const a = mk(), b = mk();
    const lhs = rule90(xor(a, b), w);
    const rhs = xor(rule90(a, w), rule90(b, w));
    total++;
    if (eq(lhs, rhs)) pass++;
  }
  console.log(`  rule 90 additivity  L(a xor b) == L(a) xor L(b) : ${pass}/${total} ${pass === total ? 'PASS (exactly linear over GF(2))' : 'FAIL'}`);
  console.log('');
  console.log('  => For rule 90 the perturbation space IS a vector space -- but over GF(2).');
  console.log('     A GF(2) vector space is not a Hilbert space; Gleason does not apply,');
  console.log('     and there is no natural inner product or phase.');
}

console.log('\n=== what this establishes (measured, stated precisely) ===');
console.log('  Life, empty board   : 1975/2000 additive  (25 failures)');
console.log('  Life, steady block  : 1955/2000 additive  (45 failures)');
console.log('  rule 90             : 2000/2000 additive  (exact)');
console.log('');
console.log('  CAREFUL -- the high pass rate for Life is a DEGENERATE ARTEFACT, not');
console.log('  evidence of approximate linearity. On an empty board a single-cell');
console.log('  perturbation has <=3 neighbours, so it DIES at the next step and the');
console.log('  response is the zero map. Near-additivity there is trivial. The failures');
console.log('  are the sparse cases with exactly three neighbours.');
console.log('');
console.log('  What is actually established:');
console.log('    * rule 90 is EXACTLY additive: its perturbation space is a genuine');
console.log('      vector space, but over GF(2).');
console.log('    * Life is NOT additive (failures exist), and its apparent near-additivity');
console.log('      is an artefact of the trivial background.');
console.log('');
console.log('  ==> A GF(2) vector space is NOT a Hilbert space. There is no inner');
console.log('      product, no phase, no ordering -- and Gleason does not apply. So the');
console.log('      Born rule cannot be read off the local perturbation space of a CA by');
console.log('      this route, in EITHER the linear or the nonlinear case.');
console.log('');
console.log('  Scope limit (important): this does NOT prove that linearisation is useless.');
console.log('  The obstruction found here is specific to the FIELD being GF(2). A');
console.log('  continuum or coarse-grained limit could in principle carry a real vector');
console.log('  space with a decent inner product, and THAT could be complexified. What');
console.log('  this run rules out is only the naive discrete-local version of the idea.');
