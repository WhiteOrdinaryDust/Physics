'use strict';
/*
 * order_effect_final.js
 * ---------------------------------------------------------------------------
 * FINAL, with a reference implementation used as a control.
 *
 * The saga so far on this one quantity has produced four different numbers:
 *     wrong target constant (2/sqrt3)      -> script error
 *     Jacobi eigenvalues 0,0,0             -> routine error
 *     discrete scan maximum 4.0            -> double-counting error
 *     hand computation 2                   -> correct so far as it goes
 *     brute force + discrete scan 0.0      -> magExpectation error
 *
 * Rather than patch again, this file builds the expectation value from a plain
 * complex matrix-vector product and VALIDATES it against cases whose answers are
 * known by inspection before using it for anything.
 * ---------------------------------------------------------------------------
 */

// ---------- complex helpers, written once
const cz = (r, i = 0) => ({ r, i });
const cmul = (a, b) => cz(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r);
const cadd = (a, b) => cz(a.r + b.r, a.i + b.i);
const cconj = (a) => cz(a.r, -a.i);

// M: 3x3 of {r,i};  v: array of 3 {r,i}  ->  M v
function matvec(M, v) {
  const out = [];
  for (let i = 0; i < 3; i++) {
    let s = cz(0);
    for (let j = 0; j < 3; j++) s = cadd(s, cmul(M[i][j], v[j]));
    out.push(s);
  }
  return out;
}
// <u|v>
function inner(u, v) {
  let s = cz(0);
  for (let i = 0; i < 3; i++) s = cadd(s, cmul(cconj(u[i]), v[i]));
  return s;
}
function norm(v) { return Math.sqrt(v.reduce((a, x) => a + x.r * x.r + x.i * x.i, 0)); }
function normalize(v) { const n = norm(v); return v.map(x => cz(x.r / n, x.i / n)); }

// expectation <v|M|v>
function expectation(M, v) { return inner(v, matvec(M, v)); }

// ---------- the operator, built by explicit matrix products
const L1 = [[cz(0), cz(1), cz(0)], [cz(1), cz(0), cz(0)], [cz(0), cz(0), cz(0)]];
const L3 = [[cz(1), cz(0), cz(0)], [cz(0), cz(-1), cz(0)], [cz(0), cz(0), cz(0)]];
function matmul(A, B) {
  const out = [];
  for (let i = 0; i < 3; i++) {
    out.push([]);
    for (let j = 0; j < 3; j++) {
      let s = cz(0);
      for (let t = 0; t < 3; t++) s = cadd(s, cmul(A[i][t], B[t][j]));
      out[i].push(s);
    }
  }
  return out;
}
const AB = matmul(L1, L3), BA = matmul(L3, L1);
const COMM = AB.map((row, i) => row.map((x, j) => cz(x.r - BA[i][j].r, x.i - BA[i][j].i)));

// ============================================================ CONTROLS FIRST
console.log('=== CONTROLS: the expectation routine must pass these ===');
{
  const e0 = expectation(COMM, [cz(1), cz(0), cz(0)]);
  console.log(`  <0|COMM|0> = ${e0.r.toFixed(6)} + ${e0.i.toFixed(6)}i   expected 0  ${Math.abs(e0.r) < 1e-12 && Math.abs(e0.i) < 1e-12 ? 'PASS' : '*** FAIL ***'}`);
  const e1 = expectation(COMM, [cz(0), cz(1), cz(0)]);
  console.log(`  <1|COMM|1> = ${e1.r.toFixed(6)} + ${e1.i.toFixed(6)}i   expected 0  ${Math.abs(e1.r) < 1e-12 && Math.abs(e1.i) < 1e-12 ? 'PASS' : '*** FAIL ***'}`);
  // |+> = (|0>+|1>)/sqrt2  ->  <+|COMM|+> = (0 + 2 + (-2) + 0)/2 = 0
  const plus = normalize([cz(1), cz(1), cz(0)]);
  const ep = expectation(COMM, plus);
  console.log(`  <+|COMM|+> = ${ep.r.toFixed(6)} + ${ep.i.toFixed(6)}i   expected 0  ${Math.abs(ep.r) < 1e-12 && Math.abs(ep.i) < 1e-12 ? 'PASS' : '*** FAIL ***'}`);
  // |x> = (|0> + i|1>)/sqrt2 -> <x|COMM|x>: compute by hand
  //   COMM|x> : row0 = -2 * (i/sqrt2), row1 = 2*(1/sqrt2), row2 = 0
  //   <x| = (1/sqrt2)(<0| - i<1|)
  //   <x|COMM|x> = (1/sqrt2)[ conj(1/sqrt2) * (-2 i/sqrt2) + conj(-i/sqrt2) * (2/sqrt2) ]
  //              = (1/2)[ -2i + (i)(2) ] = (1/2)[ -2i + 2i ] = 0
  const x1 = normalize([cz(1), cz(0, 1), cz(0)]);
  const ex = expectation(COMM, x1);
  console.log(`  <x|COMM|x> = ${ex.r.toFixed(6)} + ${ex.i.toFixed(6)}i   expected 0  ${Math.abs(ex.r) < 1e-12 && Math.abs(ex.i) < 1e-12 ? 'PASS' : '*** FAIL ***'}`);
  // |y> = (|0> - |1>)/sqrt2 rotated: take (i, 1, 0)/sqrt2
  //   COMM|y>: row0 = -2*(1/sqrt2) = -sqrt2 ; row1 = 2*(i/sqrt2)
  //   <y| = (1/sqrt2)(-i, 1, 0)
  //   <y|COMM|y> = (1/2)[ (-i)(-sqrt2) + (1)(2i/sqrt2)·... ] careful below
  const y1 = normalize([cz(0, 1), cz(1), cz(0)]);
  const ey = expectation(COMM, y1);
  // hand: <y|COMM|y> = (1/2)[ conj(i)*(-2*1) + conj(1)*(2*i) ]/1 ... do numerically
  console.log(`  <y|COMM|y> = ${ey.r.toFixed(6)} + ${ey.i.toFixed(6)}i   (no closed form asserted here)`);
  // The decisive control: an eigenvector of COMM must give an eigenvalue.
  // COMM has eigenvalues {0,0,2} with eigenvector for 2 being (|0>+i|1>)/sqrt2? check:
  //   COMM (1,i,0)/sqrt2 = (row0: -2*i/sqrt2 = -i*sqrt2, row1: 2*1/sqrt2 = sqrt2)
  //   is that 2 * (1,i,0)/sqrt2 = (sqrt2, i*sqrt2)?  row0 = -i sqrt2 vs i sqrt2 -> signs differ
  // so instead: eigenvector for +2 is (i,1,0)/sqrt2? COMM(i,1,0)/sqrt2 = ( -2/sqrt2, 2i/sqrt2 ) = (-sqrt2, i sqrt2)
  // and 2*(i,1,0)/sqrt2 = (i sqrt2, sqrt2). Not equal either. The eigenvector is (-i,1,0)/sqrt2:
  const eig = normalize([cz(0, -1), cz(1), cz(0)]);
  const ee = expectation(COMM, eig);
  console.log(`  <eig|COMM|eig> = ${ee.r.toFixed(6)} + ${ee.i.toFixed(6)}i   expected 2  ${Math.abs(ee.r - 2) < 1e-9 ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ============================================================ now the scans
console.log('=== exact maximum over ALL complex unit states (direct search) ===');
{
  let s = 20240607 >>> 0;
  const rnd = () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
  let best = 0, bestPsi = null;
  for (let trial = 0; trial < 300000; trial++) {
    const v = normalize([0, 1, 2].map(() => cz(rnd() - 0.5, rnd() - 0.5)));
    const e = expectation(COMM, v);
    const mag = Math.hypot(e.r, e.i);
    if (mag > best) { best = mag; bestPsi = v; }
  }
  console.log(`  random search (300k): max |<COMM>| = ${best.toFixed(6)}`);
  console.log(`  (exact answer is 2, from eigenvalues 0,0,2)`);
}

console.log('\n=== discrete phase scan (N-th roots of unity) ===');
console.log('     N    max |<COMM>|    fraction of 2');
for (const N of [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24, 32]) {
  let best = 0;
  for (let k1 = 0; k1 < N; k1++) for (let k2 = 0; k2 < N; k2++) for (let k3 = 0; k3 < N; k3++) {
    const v = [k1, k2, k3].map(k => {
      const t = 2 * Math.PI * k / N;
      return cz(Math.cos(t) / Math.sqrt(3), Math.sin(t) / Math.sqrt(3));
    });
    const e = expectation(COMM, v);
    const mag = Math.hypot(e.r, e.i);
    if (mag > best) best = mag;
  }
  console.log(`  ${String(N).padStart(4)}    ${best.toFixed(6)}          ${(best / 2 * 100).toFixed(2)}%`);
}

console.log('\n=== what this says ===');
console.log('  The operator [lambda_1, lambda_3] is REAL and has eigenvalues {0,0,2},');
console.log('  so the largest possible order effect is 2.');
console.log('  N=1 and N=2 (phases only 0 or pi) are the REAL-amplitude case.');
console.log('  Read the table above to see which N reaches the value 2.');
