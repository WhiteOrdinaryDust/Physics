'use strict';
/*
 * fully_discrete_test.js
 * ---------------------------------------------------------------------------
 * IS A FULLY DISCRETE MODEL ENOUGH? (amplitudes discrete AND phases discrete)
 *
 * This is the direct test of the user's intuition: reality is discrete, and
 * continuity is only a coarse-grained appearance.
 *
 * WHAT WE KNOW GOING IN
 *   exact optimum state for the order effect : (i, 1, 0)/sqrt2
 *   exact maximum of |<psi| [lambda_1,lambda_3] |psi>| : 2
 *   equal-amplitude + N discrete phases      : stalls at 4/3 (66.7%), N-independent
 *   discrete phases + FREE amplitudes (N=4)  : reaches 2 exactly
 *
 * The optimum needs THREE things at once:
 *      amplitude 1/sqrt2,  amplitude 0,  relative phase 90 degrees
 * So we scan two parameters:
 *      N = number of discrete phase values (N-th roots of unity)
 *      M = number of amplitude levels (0, 1/M, 2/M, ..., 1), normalised
 *
 * CONTROLS FIRST (this project's recurring failure mode is skipping them)
 *      |0> and |1> must give 0
 *      the exact optimum must give 2
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
function expect(M, v) { const e = inner(v, matvec(M, v)); return Math.hypot(e.r, e.i); }
function norm(v) { return Math.sqrt(v.reduce((a, x) => a + x.r * x.r + x.i * x.i, 0)); }
function normalize(v) { const n = norm(v); return v.map(x => cz(x.r / n, x.i / n)); }

const L1 = [[cz(0), cz(1), cz(0)], [cz(1), cz(0), cz(0)], [cz(0), cz(0), cz(0)]];
const L3 = [[cz(1), cz(0), cz(0)], [cz(0), cz(-1), cz(0)], [cz(0), cz(0), cz(0)]];
const AB = matmul(L1, L3), BA = matmul(L3, L1);
const COMM = AB.map((row, i) => row.map((x, j) => cz(x.r - BA[i][j].r, x.i - BA[i][j].i)));

const EXACT = 2;

console.log('=== CONTROLS (must pass before the scan is trusted) ===');
{
  const cases = [
    ['|0>', [cz(1), cz(0), cz(0)], 0],
    ['|1>', [cz(0), cz(1), cz(0)], 0],
    ['|2>', [cz(0), cz(0), cz(1)], 0],
    ['(i,1,0)/sqrt2  [the optimum]', normalize([cz(0, 1), cz(1), cz(0)]), 2],
    ['(-i,1,0)/sqrt2 [the other one]', normalize([cz(0, -1), cz(1), cz(0)]), 2],
    ['(1,1,0)/sqrt2', normalize([cz(1), cz(1), cz(0)]), 0],
  ];
  let allPass = true;
  for (const [name, v, want] of cases) {
    const got = expect(COMM, v);
    const ok = Math.abs(got - want) < 1e-9;
    if (!ok) allPass = false;
    console.log(`  ${name.padEnd(32)} got ${got.toFixed(6)}  want ${want}  ${ok ? 'PASS' : '*** FAIL ***'}`);
  }
  if (!allPass) { console.log('\n*** a control failed, aborting ***'); process.exit(1); }
}
console.log('  all controls passed\n');

// ---------- the fully discrete state set
// amplitudes: integer triples (p,q,r) with p+q+r = M  -> normalised
// phases: k*2*pi/N
function fullyDiscreteStates(N, M) {
  const out = [];
  const phases = [];
  for (let k = 0; k < N; k++) phases.push(cz(Math.cos(2 * Math.PI * k / N), Math.sin(2 * Math.PI * k / N)));
  for (let p = 0; p <= M; p++) {
    for (let q = 0; p + q <= M; q++) {
      const r = M - p - q;
      const nrm = Math.hypot(p, q, r);
      if (nrm < 1e-12) continue;
      const amps = [cz(p / nrm), cz(q / nrm), cz(r / nrm)];
      for (const a of phases) for (const b of phases) for (const c of phases) {
        out.push([cmul(amps[0], a), cmul(amps[1], b), cmul(amps[2], c)]);
      }
    }
  }
  return out;
}

console.log('=== FULLY DISCRETE: amplitudes AND phases both discrete ===');
console.log('     N = number of phase values, M = amplitude levels');
console.log('');
console.log('     N \\ M    ' + [2, 4, 6, 8, 12, 16].map(m => String(m).padStart(8)).join(''));

const grid = {};
for (const N of [2, 3, 4, 5, 6, 8, 12, 16]) {
  let line = '  ' + String(N).padStart(4) + '    ';
  for (const M of [2, 4, 6, 8, 12, 16]) {
    let best = 0;
    for (const v of fullyDiscreteStates(N, M)) {
      const e = expect(COMM, v);
      if (e > best) best = e;
    }
    grid[`${N},${M}`] = best;
    line += String(best.toFixed(4)).padStart(8);
  }
  console.log(line);
}

console.log('\n=== exact saturation of the value 2 ===');
const hits = [];
for (const N of [2, 3, 4, 5, 6, 8, 12, 16]) {
  for (const M of [2, 4, 6, 8, 12, 16]) {
    const v = grid[`${N},${M}`];
    if (Math.abs(v - EXACT) < 1e-9) hits.push({ N, M, v });
  }
}
console.log(`  (N, M) pairs that reach exactly 2: ${hits.length}`);
for (const h of hits.slice(0, 20)) console.log(`    N=${String(h.N).padStart(2)}  M=${String(h.M).padStart(2)}   ->  ${h.v.toFixed(6)}`);

// ---------- what is the minimally demanding pair?
console.log('\n=== minimal requirements ===');
{
  // smallest N needed if M is generous; smallest M needed if N is generous
  let minN = null, minM = null;
  for (const N of [2, 3, 4, 5, 6, 8, 12, 16]) {
    for (const M of [2, 4, 6, 8, 12, 16, 24, 32]) {
      let best = 0;
      for (const v of fullyDiscreteStates(N, M)) { const e = expect(COMM, v); if (e > best) best = e; }
      if (Math.abs(best - EXACT) < 1e-9) {
        if (minN === null || N < minN) minN = N;
        if (minM === null || M < minM) minM = M;
      }
    }
  }
  console.log(`  smallest N that ever reaches 2 : ${minN}`);
  console.log(`  smallest M that ever reaches 2 : ${minM}`);
}

// ---------- WHY does (N=4, M=4) reach it? show the state that does it
console.log('\n=== which state actually achieves it? (N=4, M=4) ===');
{
  let best = 0, bestV = null;
  for (const v of fullyDiscreteStates(4, 4)) {
    const e = expect(COMM, v);
    if (e > best) { best = e; bestV = v; }
  }
  const mags = bestV.map(x => Math.hypot(x.r, x.i));
  const angs = bestV.map(x => (Math.atan2(x.i, x.r) * 180 / Math.PI).toFixed(3));
  console.log(`  amplitude magnitudes : ${mags.map(x => x.toFixed(6)).join(', ')}`);
  console.log(`  phases (degrees)     : ${angs.join(', ')}`);
  console.log(`  expectation          : ${best.toFixed(6)}`);
  console.log('');
  console.log('  compare with the exact optimum (i,1,0)/sqrt2:');
  console.log(`     magnitudes 1/sqrt2, 1/sqrt2, 0  =  ${(1 / Math.SQRT2).toFixed(6)}, ${(1 / Math.SQRT2).toFixed(6)}, 0`);
  console.log('     phases     90, 0, -             =  90, 0, -');
}

console.log('\n=== conclusion ===');
console.log('  If the grid above shows 2.0000 at some (N, M) with both N and M SMALL,');
console.log('  then a FULLY DISCRETE model suffices exactly, and continuity is never');
console.log('  needed - the user\'s intuition holds.');
console.log('  If saturation needs large N or large M, the requirement is a discretisation');
console.log('  fine enough to contain the right structure, which is a weaker but still');
console.log('  non-trivial demand.');
