'use strict';
/*
 * finite_state_sufficiency.js
 * ---------------------------------------------------------------------------
 * THE CLEAN QUESTION: is a FINITE set of states enough?
 *
 * Where we are. The order-effect test gave, for a model with three EQUAL
 * amplitudes and N discrete phases:
 *        N = 3, 6            ->  1.154701  (= 2/sqrt3)
 *        N = 4, 8, 16, 24, 32 ->  1.333333  (= 4/3)
 *        exact quantum max    ->  2
 * So the discrete-phase model saturates at 2/3 of the quantum value, and raising
 * N does not help.
 *
 * The diagnosis: the optimum requires the state (i, 1, 0)/sqrt2, which has UNEQUAL
 * amplitudes (the third is zero). The model forced all three amplitudes equal, so
 * it could never reach it.
 *
 * THEREFORE the limitation is not "discreteness". It is WHICH STATES ARE
 * AVAILABLE. This file makes that precise with a control the opposite way: add
 * the single optimum state to a finite set and check that the maximum is now
 * reached exactly.
 *
 * If one extra state saturates it, then "continuous vs discrete" was the wrong
 * axis all along: what matters is whether the reachable set of states contains
 * the right structure.
 * ---------------------------------------------------------------------------
 */

const cz = (r, i = 0) => ({ r, i });
const cmul = (a, b) => cz(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r);
const cadd = (a, b) => cz(a.r + b.r, a.i + b.i);
const cconj = (a) => cz(a.r, -a.i);

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
function matvec(M, v) {
  return [0, 1, 2].map(i => {
    let s = cz(0);
    for (let j = 0; j < 3; j++) s = cadd(s, cmul(M[i][j], v[j]));
    return s;
  });
}
function inner(u, v) {
  let s = cz(0);
  for (let i = 0; i < 3; i++) s = cadd(s, cmul(cconj(u[i]), v[i]));
  return s;
}
function norm(v) { return Math.sqrt(v.reduce((a, x) => a + x.r * x.r + x.i * x.i, 0)); }
function normalize(v) { const n = norm(v); return v.map(x => cz(x.r / n, x.i / n)); }
function expect(M, v) { const e = inner(v, matvec(M, v)); return Math.hypot(e.r, e.i); }

const L1 = [[cz(0), cz(1), cz(0)], [cz(1), cz(0), cz(0)], [cz(0), cz(0), cz(0)]];
const L3 = [[cz(1), cz(0), cz(0)], [cz(0), cz(-1), cz(0)], [cz(0), cz(0), cz(0)]];
const AB = matmul(L1, L3), BA = matmul(L3, L1);
const COMM = AB.map((row, i) => row.map((x, j) => cz(x.r - BA[i][j].r, x.i - BA[i][j].i)));

console.log('=== controls (must pass before anything is used) ===');
{
  const c1 = expect(COMM, [cz(1), cz(0), cz(0)]);
  const c2 = expect(COMM, [cz(0), cz(1), cz(0)]);
  console.log(`  |<0|COMM|0>| = ${c1.toFixed(6)}  expected 0   ${c1 < 1e-12 ? 'PASS' : '*** FAIL ***'}`);
  console.log(`  |<1|COMM|1>| = ${c2.toFixed(6)}  expected 0   ${c2 < 1e-12 ? 'PASS' : '*** FAIL ***'}`);
  const opt = normalize([cz(0, 1), cz(1), cz(0)]);
  const c3 = expect(COMM, opt);
  console.log(`  |<opt|COMM|opt>| = ${c3.toFixed(6)}  expected 2   ${Math.abs(c3 - 2) < 1e-9 ? 'PASS' : '*** FAIL ***'}`);
  const opt2 = normalize([cz(0, -1), cz(1), cz(0)]);
  const c4 = expect(COMM, opt2);
  console.log(`  |<opt'|COMM|opt'>| = ${c4.toFixed(6)}  expected 2  ${Math.abs(c4 - 2) < 1e-9 ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ---------- the two state sets we compare
// SET A: equal amplitudes, N discrete phases  (what the "discrete phase" model gives)
function setEqualDiscrete(N) {
  const out = [];
  for (let k1 = 0; k1 < N; k1++) for (let k2 = 0; k2 < N; k2++) for (let k3 = 0; k3 < N; k3++) {
    out.push([k1, k2, k3].map(k => {
      const t = 2 * Math.PI * k / N;
      return cz(Math.cos(t) / Math.sqrt(3), Math.sin(t) / Math.sqrt(3));
    }));
  }
  return out;
}
// SET B: the SAME discrete phase values, but with amplitude freedom
//   i.e. we allow any (a1,a2,a3) real non-negative with a1^2+a2^2+a3^2=1,
//   sampled on a coarse grid, and phases from the N-th roots of unity.
function setFreeAmplitudeDiscrete(N, steps) {
  const out = [];
  const phase = (k) => cz(Math.cos(2 * Math.PI * k / N), Math.sin(2 * Math.PI * k / N));
  for (let i = 0; i <= steps; i++) {
    for (let j = 0; i + j <= steps; j++) {
      const k = steps - i - j;
      const a = [i / steps, j / steps, k / steps];
      const nrm = Math.hypot(...a);
      if (nrm < 1e-12) continue;
      for (let p1 = 0; p1 < N; p1++) for (let p2 = 0; p2 < N; p2++) for (let p3 = 0; p3 < N; p3++) {
        out.push([
          cmul(cz(a[0] / nrm), phase(p1)),
          cmul(cz(a[1] / nrm), phase(p2)),
          cmul(cz(a[2] / nrm), phase(p3)),
        ]);
      }
    }
  }
  return out;
}

console.log('=== SET A: equal amplitudes, N discrete phases (the previous model) ===');
console.log('     N    max |<COMM>|    fraction of 2');
for (const N of [3, 4, 6, 8, 12, 16]) {
  let best = 0;
  for (const v of setEqualDiscrete(N)) best = Math.max(best, expect(COMM, v));
  console.log(`  ${String(N).padStart(4)}    ${best.toFixed(6)}          ${(best / 2 * 100).toFixed(2)}%`);
}

console.log('\n=== SET B: SAME discrete phases, but amplitudes allowed to vary ===');
console.log('     N  steps    max |<COMM>|    fraction of 2');
for (const N of [3, 4, 6]) {
  for (const steps of [4, 8, 16]) {
    let best = 0;
    for (const v of setFreeAmplitudeDiscrete(N, steps)) best = Math.max(best, expect(COMM, v));
    console.log(`  ${String(N).padStart(4)}  ${String(steps).padStart(5)}    ${best.toFixed(6)}          ${(best / 2 * 100).toFixed(2)}%`);
  }
}

console.log('\n=== the decisive control: add ONLY the single optimum state ===');
{
  // start from the equal-amplitude discrete set, then add the one optimum state
  const N = 4;
  const base = setEqualDiscrete(N);
  let bestBase = 0;
  for (const v of base) bestBase = Math.max(bestBase, expect(COMM, v));
  const opt = normalize([cz(0, 1), cz(1), cz(0)]);
  let bestWithOpt = 0;
  for (const v of base) bestWithOpt = Math.max(bestWithOpt, expect(COMM, v));
  bestWithOpt = Math.max(bestWithOpt, expect(COMM, opt));
  console.log(`  equal-amplitude discrete set (N=4) alone : ${bestBase.toFixed(6)}  (${(bestBase / 2 * 100).toFixed(1)}%)`);
  console.log(`  the SAME set plus ONE extra state        : ${bestWithOpt.toFixed(6)}  (${(bestWithOpt / 2 * 100).toFixed(1)}%)`);
  console.log('  => a single additional state saturates the bound exactly.');
}

console.log('\n=== what this establishes ===');
console.log('  1. The earlier limitation was NOT discreteness. Raising N inside the');
console.log('     equal-amplitude model does not help at all (it stalls at 4/3 = 66.7%).');
console.log('  2. Allowing the AMPLITUDES to vary, with the very same discrete phases,');
console.log('     lets the bound be approached much more closely.');
console.log('  3. ONE extra state - the optimum itself - saturates it exactly.');
console.log('');
console.log('  So the axis that matters is neither "continuous vs discrete" nor "how many');
console.log('  phases". It is WHETHER THE REACHABLE SET OF STATES CONTAINS THE RIGHT');
console.log('  STRUCTURE: two amplitudes with a relative phase of 90 degrees.');
console.log('');
console.log('  That is a statement about the AVAILABLE STATE SPACE, and it is exactly');
console.log('  what a complex vector space with arbitrary amplitudes supplies.');
