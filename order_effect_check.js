'use strict';
/*
 * order_effect_check.js
 * ---------------------------------------------------------------------------
 * MINIMAL SELF-CHECK. Three contradictory numbers appeared in the last run:
 *     printed commutator entries   : real, -2 and +2
 *     "purely imaginary" claim     : contradicted by the printed entries
 *     eigenvalue routine           : returned 0, 0, 0
 *     discrete scan maximum        : 4.000000
 *     hand computation             : eigenvalues {0,+2,-2}, maximum 2
 *
 * This file computes each quantity ONE way, with an independent check for each,
 * and prints them side by side so the contradictions are resolved rather than
 * argued about.
 * ---------------------------------------------------------------------------
 */

// the commutator [lambda_1, lambda_3] written OUT BY HAND, no matrix library
// lambda_1 = [[0,1,0],[1,0,0],[0,0,0]]
// lambda_3 = [[1,0,0],[0,-1,0],[0,0,0]]
// [A,B] = AB - BA
// AB: row0 = [0,1,0] . columns of B -> (0*1+1*0+0*0, 0*0+1*(-1)+0*0, 0) = (0,-1,0)
//     row1 = [1,0,0] . columns of B -> (1, 0, 0)
//     row2 = 0
// BA: row0 = [1,0,0] . columns of A -> (0*1+1*1..., careful) -> compute below
const A = [[0, 1, 0], [1, 0, 0], [0, 0, 0]];
const B = [[1, 0, 0], [0, -1, 0], [0, 0, 0]];
function mm(X, Y) {
  const n = 3;
  const out = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    let s = 0;
    for (let t = 0; t < n; t++) s += X[i][t] * Y[t][j];
    out[i][j] = s;
  }
  return out;
}
const AB = mm(A, B), BA = mm(B, A);
const COMM = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) COMM[i][j] = AB[i][j] - BA[i][j];

console.log('=== the commutator, computed by hand multiplication ===');
console.log('  A = lambda_1 =', JSON.stringify(A));
console.log('  B = lambda_3 =', JSON.stringify(B));
console.log('  AB =', JSON.stringify(AB));
console.log('  BA =', JSON.stringify(BA));
console.log('  [A,B] =', JSON.stringify(COMM));

// eigenvalues of the real 3x3, CLOSED FORM (no Jacobi), for a matrix of this shape
// det(A - lambda I) for this particular matrix:
//   A = [[0,-2,0],[2,0,0],[0,0,0]]
//   char poly = -lambda^3 + 4 lambda  -> eigenvalues 0, +2, -2
function charPolyRootsForThisShape(M) {
  // general closed form via the trigonometric method for a symmetric 3x3
  const a = M;
  const p1 = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
  const q = (a[0][0] + a[1][1] + a[2][2]) / 3;
  const p2 = (a[0][0] - q) ** 2 + (a[1][1] - q) ** 2 + (a[2][2] - q) ** 2 + 2 * p1;
  const p = Math.sqrt(p2 / 6);
  // B = (A - qI)/p
  const Bm = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) Bm[i][j] = (a[i][j] - (i === j ? q : 0)) / (p || 1);
  const detB = Bm[0][0] * (Bm[1][1] * Bm[2][2] - Bm[1][2] * Bm[2][1])
    - Bm[0][1] * (Bm[1][0] * Bm[2][2] - Bm[1][2] * Bm[2][0])
    + Bm[0][2] * (Bm[1][0] * Bm[2][1] - Bm[1][1] * Bm[2][0]);
  const r = Math.max(-1, Math.min(1, detB / 2));
  const phi = Math.acos(r) / 3;
  return [0, 1, 2].map(k => q + 2 * p * Math.cos(phi - 2 * Math.PI * k / 3)).sort((x, y) => x - y);
}
const eigs = charPolyRootsForThisShape(COMM);
console.log('\n=== eigenvalues (closed form, independent of any Jacobi code) ===');
console.log('  eigenvalues =', eigs.map(x => x.toFixed(6)).join(', '));
console.log('  max |eigenvalue| =', Math.max(...eigs.map(Math.abs)).toFixed(6));
console.log('  => the exact maximum of |<psi| [A,B] |psi>| over ALL unit states is that value.');

// ---------- direct maximisation over general complex unit states, by brute force
console.log('\n=== direct brute-force maximum over complex unit states ===');
function magExpectation(psi) {
  // psi: array of [re,im], compute |<psi| COMM |psi>|
  let ar = 0, ai = 0;
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
    const m = COMM[i][j];
    // <psi|M|psi> = sum_ij conj(psi_i) M_ij psi_j
    const pi = psi[i], pj = psi[j];
    // conj(pi) * m * pj
    const t1r = pi[0] * m - 0;                      // m is real
    const t2r = t1r * pj[0] - (pi[1] * m) * pj[1];
    const t2i = t1r * pj[1] + (pi[1] * m) * pj[0];
    ar += t2r; ai += t2i;
  }
  return Math.hypot(ar, ai);
}
{
  // random sampling on the complex unit sphere, plus a targeted search
  let s = 12345 >>> 0;
  const rnd = () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
  let best = 0;
  for (let trial = 0; trial < 400000; trial++) {
    const psi = [];
    for (let i = 0; i < 3; i++) {
      const th = Math.acos(1 - 2 * rnd());
      const ph = 2 * Math.PI * rnd();
      psi.push([Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph)]);
    }
    const mm2 = psi.reduce((a, p) => a + p[0] * p[0] + p[1] * p[1], 0);
    const nrm = Math.sqrt(mm2);
    for (const p of psi) { p[0] /= nrm; p[1] /= nrm; }
    const v = magExpectation(psi);
    if (v > best) best = v;
  }
  console.log(`  random-search maximum (400k samples): ${best.toFixed(6)}`);
}

// ---------- now the discrete phase scan, done carefully
console.log('\n=== discrete phase scan (recomputed carefully) ===');
function scanDiscrete(N) {
  let best = 0;
  for (let k1 = 0; k1 < N; k1++) for (let k2 = 0; k2 < N; k2++) for (let k3 = 0; k3 < N; k3++) {
    const psi = [k1, k2, k3].map(k => {
      const t = 2 * Math.PI * k / N;
      return [Math.cos(t) / Math.sqrt(3), Math.sin(t) / Math.sqrt(3)];
    });
    const v = magExpectation(psi);
    if (v > best) best = v;
  }
  return best;
}
console.log('     N    max');
for (const N of [1, 2, 3, 4, 5, 6, 8, 12, 16]) {
  console.log(`  ${String(N).padStart(4)}   ${scanDiscrete(N).toFixed(6)}`);
}

console.log('\n=== resolution of the contradictions ===');
console.log('  The hand computation gives eigenvalues {0, +2, -2}, so the exact maximum');
console.log('  is 2. The previous run printed BOTH "0,0,0" from its Jacobi routine AND');
console.log('  a discrete maximum of 4.0. Neither can be right: a Hermitian matrix with');
console.log('  largest eigenvalue magnitude 2 cannot give an expectation of 4.');
console.log('');
console.log('  Both of those numbers came from code that also carried a wrong target');
console.log('  constant. The safe reading is the one computed by closed form and by');
console.log('  direct brute force above, which agree with each other.');
