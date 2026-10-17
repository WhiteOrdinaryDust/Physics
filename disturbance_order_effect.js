'use strict';
/*
 * disturbance_order_effect.js
 * ---------------------------------------------------------------------------
 * CAN PURE CLASSICAL DISTURBANCE PRODUCE THE ORDER EFFECT?
 *
 * The user's mechanism, stated precisely: an internal measurement necessarily
 * changes the system (S_n depends on S_{n-1}, and measuring changes S_{n-1}), so
 * the later evolution differs. Therefore measurement ORDER matters.
 *
 * That mechanism is real. The question is whether it is ENOUGH.
 *
 * WHY THE OBVIOUS TEST FAILS
 *   Compare "measure A then B" with "measure B then A" and look at the joint
 *   distribution. This gives a nonzero difference, but in a deterministic world
 *   the difference VANISHES in the marginal of either measurement:
 *        sum_b  P(a,b | A then B)  =  P(a)  =  sum_b P(a,b | B then A)
 *   because summing over the first outcome erases whatever the disturbance did.
 *   So marginal-based order effects cannot distinguish anything here.
 *
 * THE QUANTITY THAT DOES MATTER
 *   The order-effect OPERATOR, i.e. the commutator C = [A,B], and its
 *   expectation value. That is what the earlier files computed, and it is what
 *   quantum mechanics predicts to have magnitude up to 2 for the pair
 *   A = lambda_1, B = lambda_3.
 *
 * THE DECISIVE QUESTION
 *   Classical disturbance acts by projection / mixing: it produces REAL
 *   processes. A real process cannot generate the purely imaginary structure
 *   that the commutator needs. So the prediction to be tested is:
 *
 *     PREDICTION (written before running):
 *       the largest value achievable by a disturbance-type (real, positivity
 *       preserving) process is EXACTLY HALF the quantum maximum.
 *       For our pair: 1.0 instead of 2.0.
 *
 *   If that holds, then disturbance is NECESSARY but NOT SUFFICIENT: the missing
 *   ingredient is exactly the coherent complex structure, i.e. INTERFERENCE.
 * ---------------------------------------------------------------------------
 */

const cz = (r, i = 0) => ({ r, i });
const cmul = (a, b) => cz(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r);
const cadd = (a, b) => cz(a.r + b.r, a.i + b.i);
const cconj = a => cz(a.r, -a.i);

function matmul(A, B) {
  const n = A.length, m = B[0].length, k = B.length;
  const out = [];
  for (let i = 0; i < n; i++) {
    out.push([]);
    for (let j = 0; j < m; j++) {
      let s = cz(0);
      for (let t = 0; t < k; t++) s = cadd(s, cmul(A[i][t], B[t][j]));
      out[i].push(s);
    }
  }
  return out;
}
function matvec(M, v) {
  return v.map((_, i) => {
    let s = cz(0);
    for (let j = 0; j < v.length; j++) s = cadd(s, cmul(M[i][j], v[j]));
    return s;
  });
}
function inner(u, v) {
  let s = cz(0);
  for (let i = 0; i < u.length; i++) s = cadd(s, cmul(cconj(u[i]), v[i]));
  return s;
}
function expect(M, v) { const e = inner(v, matvec(M, v)); return e; }
function norm(v) { return Math.sqrt(v.reduce((a, x) => a + x.r * x.r + x.i * x.i, 0)); }
function normalize(v) { const n = norm(v); return v.map(x => cz(x.r / n, x.i / n)); }

const L1 = [[cz(0), cz(1), cz(0)], [cz(1), cz(0), cz(0)], [cz(0), cz(0), cz(0)]];
const L3 = [[cz(1), cz(0), cz(0)], [cz(0), cz(-1), cz(0)], [cz(0), cz(0), cz(0)]];
const AB = matmul(L1, L3), BA = matmul(L3, L1);
const C = AB.map((row, i) => row.map((x, j) => cz(x.r - BA[i][j].r, x.i - BA[i][j].i)));

console.log('=== controls ===');
{
  const cases = [
    ['|0>', [cz(1), cz(0), cz(0)], 0],
    ['(i,1,0)/sqrt2  [the optimum]', normalize([cz(0, 1), cz(1), cz(0)]), 2],
    ['(1,1,0)/sqrt2  [real state]', normalize([cz(1), cz(1), cz(0)]), 0],
  ];
  let ok = true;
  for (const [nm, v, want] of cases) {
    const e = expect(C, v);
    const mag = Math.hypot(e.r, e.i);
    const pass = Math.abs(mag - want) < 1e-9;
    if (!pass) ok = false;
    console.log(`  ${nm.padEnd(32)} |<C>| = ${mag.toFixed(6)}  want ${want}  ${pass ? 'PASS' : '*** FAIL ***'}`);
  }
  if (!ok) { console.log('*** controls failed, aborting ***'); process.exit(1); }
}
console.log('  controls passed\n');

// ============================================================ the test
/*
 * A "disturbance-type" (classical) process maps a state to another state by
 * projection and mixing ONLY. Concretely, we model the reachable states of a
 * classical disturbance model as the real, non-negative mixtures:
 *
 *     rho_classical  =  sum_k p_k |s_k><s_k| ,
 *     p_k >= 0, sum p_k = 1, |s_k> a REAL unit vector in the computational basis
 *
 * i.e. distributions over real basis-like states. This is the state space of a
 * classical stochastic system written in quantum notation.
 *
 * We compute the largest |Tr(rho C)| over that set and compare with the quantum
 * maximum 2.
 */
console.log('=== TEST: largest |<C>| achievable by a classical (real, positive) ensemble ===');

// The classical ensemble: probability vectors over the three real basis states
// |1>, |2>, |3> (and their real superposition is NOT available classically, since
// a classical system cannot be in a coherent superposition -- only in a mixture).
// So the extreme points are the three basis states themselves.
console.log('\n  extreme points (the "pure" classical states) and their |<C>|:');
let classicalMax = 0;
for (let k = 0; k < 3; k++) {
  const v = [0, 1, 2].map(j => cz(j === k ? 1 : 0));
  const e = expect(C, v);
  const mag = Math.hypot(e.r, e.i);
  classicalMax = Math.max(classicalMax, mag);
  console.log(`    |${k}>   ->  |<C>| = ${mag.toFixed(6)}`);
}

// mixtures can only interpolate between extreme points (|.<C>.| is convex in rho),
// so the maximum over mixtures is again attained at an extreme point.
console.log(`\n  maximum over all classical mixtures = ${classicalMax.toFixed(6)}`);

// Now allow the classical model the best case it can have: a REAL state that is a
// real superposition (this is the most generous "classical + real amplitude" case,
// which is what a real-valued wave function would give).
console.log('\n  most generous real-amplitude case (real superposition allowed):');
let realMax = 0, bestReal = null;
{
  // search real unit vectors on the 2-sphere
  const N = 400;
  for (let i = 0; i <= N; i++) {
    const th = Math.PI * i / N;
    for (let j = 0; j <= N; j++) {
      const ph = 2 * Math.PI * j / N;
      const v = [cz(Math.sin(th) * Math.cos(ph)), cz(Math.sin(th) * Math.sin(ph)), cz(Math.cos(th))];
      const e = expect(C, v);
      const mag = Math.hypot(e.r, e.i);
      if (mag > realMax) { realMax = mag; bestReal = v; }
    }
  }
}
console.log(`    maximum over ALL real unit vectors = ${realMax.toFixed(6)}`);

console.log('\n=== RESULT ===');
console.log(`  classical (mixtures of basis states) : ${classicalMax.toFixed(6)}`);
console.log(`  real superpositions allowed          : ${realMax.toFixed(6)}`);
console.log(`  complex (quantum) maximum            : 2.000000  (from eigenvalues 0,0,2)`);
console.log('');
console.log('  PREDICTION STATED BEFORE THE RUN: the real/classical maximum is EXACTLY');
console.log('  HALF the quantum maximum, i.e. 1.0.');
console.log(`  OBSERVED: real maximum = ${realMax.toFixed(6)}`);

console.log('\n=== interpretation ===');
if (Math.abs(realMax - 1) < 1e-6) {
  console.log('  The prediction holds exactly: real processes reach 1.0, complex reach 2.0.');
  console.log('  => Disturbance is necessary but NOT sufficient. Its reach is exactly half.');
  console.log('  => The missing factor of 2 is the INTERFERENCE term.');
} else if (realMax > 1.5) {
  console.log('  The prediction FAILED: real processes reach much more than half.');
  console.log('  => My reasoning about "real processes cannot do it" was too strong.');
} else {
  console.log(`  The real maximum is ${realMax.toFixed(6)}; compare with the predicted 1.0.`);
  console.log('  Read the number rather than the narrative.');
}
