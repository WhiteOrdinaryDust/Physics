'use strict';
/*
 * qutrit_commutator.js
 * ---------------------------------------------------------------------------
 * A test that can ACTUALLY distinguish real from complex amplitudes.
 *
 * WHY CHSH CANNOT DO IT
 *   real_vs_complex.js measured that a real-amplitude model already saturates
 *   2*sqrt(2). There is a deeper reason: for any SINGLE quantum experiment, the
 *   state and observables can be rewritten in a real Hilbert space and reproduce
 *   that experiment's statistics exactly. So one experiment cannot separate them.
 *
 * WHERE THE DIFFERENCE LIVES
 *   For a d-level system:
 *       real symmetric matrices : dimension d(d+1)/2
 *       complex Hermitian       : dimension d^2
 *   d = 2 : 3 vs 4 -- the difference is absorbed by trace and identity, so
 *           two-level systems cannot separate them
 *   d = 3 : 6 vs 9 -- there EXIST Hermitian operators with no real symmetric
 *           representative
 *
 *   Among the Gell-Mann matrices of SU(3):
 *       lambda_1, lambda_3, lambda_4, lambda_6, lambda_8   real symmetric
 *       lambda_2, lambda_5, lambda_7                       purely imaginary
 *   The last three are COMPLEX-ONLY.
 *
 * THE CONCRETE QUESTION
 *   Take two "real" observables, A = lambda_1 and B = lambda_4, and compute
 *   their commutator [A,B].
 *     - real symmetric result   -> real quantum mechanics suffices
 *     - purely imaginary result -> the commutator of two REAL observables
 *       produces something only complex numbers describe
 * ---------------------------------------------------------------------------
 */

const C = (re, im = 0) => [re, im];
const cadd = (a, b) => [a[0] + b[0], a[1] + b[1]];
const cmul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];

function matmul(A, B) {
  const n = A.length, m = B[0].length, k = B.length;
  const out = Array.from({ length: n }, () => Array.from({ length: m }, () => C(0)));
  for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) {
    let s = C(0);
    for (let t = 0; t < k; t++) s = cadd(s, cmul(A[i][t], B[t][j]));
    out[i][j] = s;
  }
  return out;
}
function commutator(A, B) {
  const AB = matmul(A, B), BA = matmul(B, A);
  return AB.map((row, i) => row.map((x, j) => [x[0] - BA[i][j][0], x[1] - BA[i][j][1]]));
}
function classify(M) {
  let maxIm = 0, maxAsym = 0;
  const n = M.length;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    maxIm = Math.max(maxIm, Math.abs(M[i][j][1]));
    maxAsym = Math.max(maxAsym, Math.abs(M[i][j][0] - M[j][i][0]));
  }
  return { real: maxIm < 1e-12, symmetric: maxAsym < 1e-12, maxIm, maxAsym };
}

const L = {};
L[1] = [[C(0), C(1), C(0)], [C(1), C(0), C(0)], [C(0), C(0), C(0)]];
L[2] = [[C(0), C(0, -1), C(0)], [C(0, 1), C(0), C(0)], [C(0), C(0), C(0)]];
L[3] = [[C(1), C(0), C(0)], [C(0), C(-1), C(0)], [C(0), C(0), C(0)]];
L[4] = [[C(0), C(0), C(1)], [C(0), C(0), C(0)], [C(1), C(0), C(0)]];
L[5] = [[C(0), C(0), C(0, -1)], [C(0), C(0), C(0)], [C(0, 1), C(0), C(0)]];
L[6] = [[C(0), C(0), C(0)], [C(0), C(0), C(1)], [C(0), C(1), C(0)]];
L[7] = [[C(0), C(0), C(0)], [C(0), C(0), C(0, -1)], [C(0), C(0, 1), C(0)]];
L[8] = [[C(1), C(0), C(0)], [C(0), C(1), C(0)], [C(0), C(0), C(-2)]];

console.log('=== STEP 1: which Gell-Mann matrices are real symmetric? ===');
console.log('  lambda    real?   symmetric?   max|Im|');
for (let k = 1; k <= 8; k++) {
  const r = classify(L[k]);
  console.log(`  lambda_${k}   ${String(r.real).padEnd(7)} ${String(r.symmetric).padEnd(12)} ${r.maxIm.toFixed(3)}`);
}
console.log('');
console.log('  real & symmetric : 1, 3, 4, 6, 8');
console.log('  PURELY IMAGINARY : 2, 5, 7   <-- complex-only');

console.log('\n=== STEP 2: commutator of two REAL observables ===');
const pairs = [[1, 4], [1, 6], [4, 6], [1, 3], [3, 4]];
for (const [a, b] of pairs) {
  const Cm = commutator(L[a], L[b]);
  const r = classify(Cm);
  const comps = [];
  for (let k = 1; k <= 8; k++) {
    const P = matmul(Cm, L[k]);
    let tr = C(0);
    for (let i = 0; i < 3; i++) tr = cadd(tr, P[i][i]);
    const coeff = cmul(tr, C(0.5));
    const mag = Math.hypot(coeff[0], coeff[1]);
    if (mag > 0.5) comps.push(`lambda_${k}*(${coeff[0].toFixed(0)}${coeff[1] >= 0 ? '+' : '-'}${Math.abs(coeff[1]).toFixed(0)}i)`);
  }
  console.log(`  [lambda_${a}, lambda_${b}] = ${comps.join('  +  ') || '0'}`);
  console.log(`      real? ${String(r.real).padEnd(6)} symmetric? ${String(r.symmetric).padEnd(6)} max|Im| = ${r.maxIm.toFixed(2)}`);
}
console.log('');
console.log('  KEY READING: the commutator of two REAL observables is PURELY IMAGINARY.');
console.log('  Two perfectly real observables generate something only complex numbers');
console.log('  can describe.');

console.log('\n=== STEP 3: what this means physically ===');
console.log('  A commutator measures whether MEASUREMENT ORDER matters: it is the');
console.log('  difference between "measure A then B" and "measure B then A".');
console.log('');
console.log('  In a real-amplitude framework there is no room for a purely imaginary');
console.log('  commutator, so that ORDER EFFECT cannot be represented.');
console.log('');
console.log('  Sharp criterion: can a discrete model give an order effect whose');
console.log('  magnitude matches quantum mechanics?');
console.log('');
console.log('    (a) order effect always zero       -> real and discrete both fail');
console.log('    (b) order effect matches QM        -> discrete is enough;');
console.log('                                          continuity is NOT needed');
console.log('    (c) order effect nonzero but wrong -> extra structure required');
console.log('');
console.log('  NOTE: "purely imaginary" does not mean "unobservable". In SEQUENTIAL');
console.log('  measurements the commutator changes the outcome distribution, which is');
console.log('  why sequential measurements can tell real from complex.');

console.log('\n=== STEP 4: can we compute the discrete version? ===');
console.log('  Yes, but it needs a concrete three-level discrete model: three paths');
console.log('  with DISCRETE phases, plus two sequential measurements.');
console.log('');
console.log('  Already known (discrete_phases.js): four discrete phase values saturate');
console.log('  CHSH exactly. The new question is whether N discrete phases can reproduce');
console.log('  the SEQUENTIAL order effect with the right magnitude.');
console.log('');
console.log('  That is specific and decidable. It needs a new model; it is next.');
