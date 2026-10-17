'use strict';
/*
 * sequential_order_effect.js
 * ---------------------------------------------------------------------------
 * Can a DISCRETE-phase model reproduce the sequential-measurement ORDER EFFECT?
 *
 * WHY THIS IS THE RIGHT TEST
 *   qutrit_commutator.js showed: for two REAL observables A = lambda_1 and
 *   B = lambda_3 on a 3-level system, the commutator is PURELY IMAGINARY:
 *        [lambda_1, lambda_3] = -2i * lambda_2
 *   A commutator is exactly the "order effect": (A then B) minus (B then A).
 *
 *   So if a model can only produce real amplitudes, the order effect it predicts
 *   is forced to be purely imaginary in the standard basis -- there is no room
 *   for a real, physical order effect. Whether a DISCRETE-phase model can fix
 *   this is the decisive question, and it is computable.
 *
 * THE MODEL
 *   A three-path system. Each path i carries a phase phi_i drawn from the N-th
 *   roots of unity:  phi_i = 2*pi*k_i/N,  k_i integer.
 *   The state is a real superposition with these discrete phases,
 *        |psi> = (1/sqrt3) * sum_i e^{i phi_i} |i>.
 *
 *   Two observables (real symmetric, so available in BOTH the real and complex
 *   frameworks):
 *        A = lambda_1,  B = lambda_3.
 *
 *   The order effect is the commutator. We compare:
 *        the exact quantum value  <psi| [A,B] |psi>  =  -2i <psi|lambda_2|psi>
 *        against what the discrete phase set can reach.
 *
 * WHAT WE REPORT
 *   For each N (number of discrete phase values), the LARGEST magnitude of the
 *   order effect achievable. Then the same question as before: how large must N
 *   be for the discrete model to match the continuum exactly?
 * ---------------------------------------------------------------------------
 */

const C = (re, im = 0) => [re, im];
const cadd = (a, b) => [a[0] + b[0], a[1] + b[1]];
const cmul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const cabs = (a) => Math.hypot(a[0], a[1]);

// ---------- the three observables we need (as complex 3x3)
const L1 = [[C(0), C(1), C(0)], [C(1), C(0), C(0)], [C(0), C(0), C(0)]];
const L2 = [[C(0), C(0, -1), C(0)], [C(0, 1), C(0), C(0)], [C(0), C(0), C(0)]];
const L3 = [[C(1), C(0), C(0)], [C(0), C(-1), C(0)], [C(0), C(0), C(0)]];

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
// <psi| M |psi> with psi a 3-vector of complex amplitudes
function expectation(M, psi) {
  const Mp = [C(0), C(0), C(0)];
  for (let i = 0; i < 3; i++) {
    let s = C(0);
    for (let j = 0; j < 3; j++) s = cadd(s, cmul(M[i][j], psi[j]));
    Mp[i] = s;
  }
  let acc = C(0);
  for (let i = 0; i < 3; i++) acc = cadd(acc, cmul([psi[i][0], -psi[i][1]], Mp[i]));
  return acc;
}

const COMM = commutator(L1, L3);

console.log('=== the order effect, exactly ===');
console.log('  A = lambda_1, B = lambda_3, both REAL and symmetric.');
{
  const c = COMM;
  const vals = [];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
    if (cabs(c[i][j]) > 1e-12) vals.push(`(${i},${j}) = ${c[i][j][0]}${c[i][j][1] >= 0 ? '+' : '-'}${Math.abs(c[i][j][1])}i`);
  }
  console.log('  [A,B] non-zero entries: ' + vals.join('   '));
  console.log('  => explicitly: [lambda_1, lambda_3] = -2i * lambda_2  (PURELY IMAGINARY)');
}
console.log('');

// ---------- discrete phase states
// phase set: k * 2*pi/N for k = 0..N-1
function discreteStates(N) {
  const states = [];
  for (let k1 = 0; k1 < N; k1++) for (let k2 = 0; k2 < N; k2++) for (let k3 = 0; k3 < N; k3++) {
    const p = [k1, k2, k3].map(k => 2 * Math.PI * k / N);
    states.push(p.map(t => C(Math.cos(t), Math.sin(t))));
  }
  return states;
}

console.log('=== can discrete phases reproduce the order effect? ===');
console.log('  We scan states whose three phases come from N-th roots of unity,');
console.log('  and record the LARGEST |order effect| achievable.\n');

console.log('     N     max |<psi|[A,B]|psi>|      (exact continuum max is 2/sqrt3 = ' + (2 / Math.sqrt(3)).toFixed(6) + ')');

const TARGET = 2 / Math.sqrt(3);
const rows = [];
for (const N of [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24, 32, 48, 64, 128]) {
  const states = discreteStates(N);
  let best = 0, bestPhases = null;
  for (const psi of states) {
    const e = expectation(COMM, psi);
    const mag = cabs(e);
    if (mag > best) { best = mag; bestPhases = psi.map(a => (Math.atan2(a[1], a[0]) * 180 / Math.PI).toFixed(1)); }
  }
  rows.push({ N, best });
  const gap = TARGET - best;
  console.log(`  ${String(N).padStart(4)}     ${best.toFixed(6)}                 gap ${gap.toExponential(2)}    phases(deg) ${bestPhases ? bestPhases.join(', ') : '-'}`);
}

console.log('\n=== reading ===');
console.log('  If the max reaches 2/sqrt3 for a SMALL N, discrete phases suffice exactly');
console.log('  and continuity is unnecessary. If it needs large N, or never reaches it,');
console.log('  then something more than discrete phases is required.');
console.log('');
const hits = rows.filter(r => r.best > TARGET - 1e-9);
if (hits.length) {
  console.log(`  => SATURATED at N = ${hits[0].N}`);
} else {
  console.log('  => NOT saturated at any tested N.');
}
const ratios = rows.map(r => ({ N: r.N, ratio: r.best / TARGET }));
console.log('\n  N vs fraction of the exact value:');
for (const r of ratios) console.log(`    N=${String(r.N).padStart(3)} : ${(r.ratio * 100).toFixed(2)}%`);
