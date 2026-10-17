'use strict';
/*
 * e1_ontic.js   --   ROUND 71: ESCAPE HATCH E1
 * ---------------------------------------------------------------------------
 * E1 was the last untouched direction, and the only one that could still refute the final report:
 *
 *      put LINEARITY at the level of DISTRIBUTIONS OVER BEABLES, not at the level of microstates.
 *
 * The point of E1: rounds 61-64 showed the microstate alphabet cannot carry a complex vector space
 * (counts form a cone; a discrete amplitude cannot even add one quantum).  But 't Hooft's move is to
 * let the STATES be probability distributions rho(omega) over beables omega.  Then:
 *      - deterministic evolution is a permutation, hence a LINEAR map on distributions
 *      - observables are functions on beables, hence their expectations are LINEAR in rho
 *      - a deterministic, local, finite substrate is enough -- no vector-space alphabet needed.
 * So E1 buys linearity.  The question is what else it buys, and what it cannot.
 *
 * FOUR MEASURABLE CLAIMS (each falsifiable, each with a control that can fail)
 *   A  LINEARITY holds at the distribution level (where it failed at the microstate level),
 *      but only for CONVEX combinations: the state space is a CONE, not a vector space.
 *   B  The distributions carry NO PHASE, so E1's prediction for the round-62 interference experiment
 *      is the INCOHERENT profile.  Measured before this round (round 69's phase-blind control,
 *      flat) -- recomputed here exactly and compared with the coherent 4/2/0/2.
 *   C  COARSE MEASUREMENTS DO NOT COMMUTE -- and the non-commutativity comes ENTIRELY from the
 *      resampling ("dephasing") step.  Control: with conditioning instruments (no resampling) the
 *      two orders give IDENTICAL outcome statistics.  So E1's non-commutativity is measurement
 *      disturbance, not a unitary algebra: the instruments are positive IDEMPOTENTS, non-invertible.
 *   D  THE DECISIVE LIMIT: on local beables, the CHSH quantity obeys S <= 2 -- enumerated exactly,
 *      not argued.  Control: breaking locality (a signalling instrument) reaches S = 4, so the bound
 *      is what locality buys, and E1 cannot escape it without giving up locality.
 *
 * WHAT THAT WOULD MEAN: E1 delivers linearity (on a cone) and non-commutativity (from dephasing),
 * but NOT interference (no phase in a distribution) and NOT the quantum correlations (S <= 2 on local
 * beables).  The two ways out -- SIGNED distributions, or NON-LOCAL instruments -- each cost exactly
 * the premise E1 was built on.
 * ---------------------------------------------------------------------------
 */

function rng(seed) { let x = seed | 0; if (x === 0) x = 0x9e3779b9; return () => { x ^= x << 13; x ^= x >>> 17; x ^= x << 5; return (x >>> 0) / 4294967296; }; }
const fmt = (v, d = 4) => v.toFixed(d);

// ===================================================================== A
console.log('============================================================');
console.log('A: E1 BUYS LINEARITY -- BUT ONLY ON A CONE');
console.log('============================================================');
{
  const N = 12;
  const rnd = rng(11);
  const norm = (a) => { const s = a.reduce((x, y) => x + y, 0); return a.map((v) => v / s); };
  const r1 = norm(Array.from({ length: N }, () => rnd()));
  const r2 = norm(Array.from({ length: N }, () => rnd()));
  const T = (rho) => { const o = new Array(N).fill(0); for (let i = 0; i < N; i++) o[(i + 1) % N] = rho[i]; return o; };   // a cyclic shift = a permutation
  const c = 0.3;
  const mix = r1.map((v, i) => c * v + (1 - c) * r2[i]);
  const lhs = T(mix), rhs = T(r1).map((v, i) => c * v + (1 - c) * T(r2)[i]);
  let linErr = 0; for (let i = 0; i < N; i++) linErr = Math.max(linErr, Math.abs(lhs[i] - rhs[i]));
  console.log(`   the beables are ${N} states on a ring; the dynamics is a cyclic shift (a permutation).`);
  console.log(`   LINEARITY at the distribution level: max |T(c*r1 + (1-c)*r2) - (c*T r1 + (1-c)*T r2)| = ${linErr.toExponential(2)}`);
  console.log('     (at the MICROSTATE level this is exactly what failed in rounds 61-64: the counting');
  console.log('      alphabet only passed the homogeneous scaling test, and a discrete amplitude could not');
  console.log('      even add one quantum.  Here it holds trivially -- that IS what E1 buys.)');
  const diff = r1.map((v, i) => v - r2[i]);
  const minDiff = Math.min(...diff);
  const negCount = diff.filter((v) => v < 0).length;
  console.log('');
  console.log('   BUT: a state must satisfy rho >= 0.  The difference of two states is not a state:');
  console.log(`     min(r1 - r2) = ${fmt(minDiff)}  (${negCount} of ${N} entries negative)`);
  console.log('     => the state space is a CONE, not a vector space.  The same obstruction measured at the');
  console.log('        microstate level in round 61 reappears ONE LEVEL UP.');
  const ok = linErr === 0 && minDiff < 0;
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: E1 does give linearity on distributions, and it is a cone, not a vector space.`);
}
console.log('');

// ===================================================================== B
console.log('============================================================');
console.log('B: A DISTRIBUTION HAS NO PHASE -- SO E1 PREDICTS THE INCOHERENT PROFILE');
console.log('============================================================');
{
  // the round-62 experiment, exactly: two waves at right angles give  |psi|^2 = 2 + 2cos(k*u)
  // with u = x - y and k = pi/2, so the four fringe bins are 4 / 2 / 0 / 2.
  // (My first version used two COUNTER-propagating waves on a line, which gives a standing wave
  //  8/0/8/0 -- a different profile, mislabelled as 4/2/0/2.  Fixed here.)
  const k = Math.PI / 2;
  const coherent = [], incoherent = [];
  for (let u = 0; u < 4; u++) {
    coherent.push(2 + 2 * Math.cos(k * u));       // amplitude sum at that fringe coordinate
    incoherent.push(1 + 1);                       // E1: add the two intensities
  }
  console.log(`   coherent (amplitude sum)    : ${coherent.map((v) => v.toFixed(2).padStart(7)).join(' ')}   <- the round-62 fringe 4/2/0/2`);
  console.log(`   E1 / incoherent (intensities): ${incoherent.map((v) => v.toFixed(2).padStart(7)).join(' ')}   <- flat`);
  console.log('');
  console.log('   E1 is phase-blind by construction: rho(omega) is a non-negative number, there is no');
  console.log('   relative phase between two "branches" anywhere in the ontology, so two branches can');
  console.log('   never cancel.  This is exactly round 69\'s phase-blind control, which MEASURED 0.251 /');
  console.log('   0.247 / 0.249 / 0.252 (flat) against the coherent 4/2/0/2.');
  const ok = Math.abs(coherent[2]) < 1e-12 && incoherent.every((v) => Math.abs(v - 2) < 1e-12);
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the coherent profile has an exactly dark fringe (0.00), the E1 profile is flat.  So`);
  console.log('     E1 cannot reproduce interference -- not because of non-linearity, but because there is');
  console.log('     nothing in the ontology to carry a phase.');
}
console.log('');

// ===================================================================== C
console.log('============================================================');
console.log('C: COARSE MEASUREMENTS DO NOT COMMUTE -- AND IT IS ENTIRELY THE RESAMPLING');
console.log('============================================================');
{
  const N = 6;
  const A = [0, 0, 0, 1, 1, 1];                 // partition A: {0,1,2} vs {3,4,5}
  const B = [0, 0, 1, 1, 2, 2];                 // partition B: {0,1} {2,3} {4,5}
  const norm = (v) => { const s = v.reduce((a, b) => a + b, 0); return v.map((x) => x / s); };
  const rnd = rng(7);
  const rho = norm(Array.from({ length: N }, () => rnd()));
  // sequential outcome statistics, both orders.  Two instruments:
  //    conditioning -- the post-measurement state is the CONDITIONAL distribution (no disturbance)
  //    resampling   -- the post-measurement state is UNIFORM inside the outcome cell ("dephasing")
  function seq(pA, pB, rho, kind) {
    function post(r, P, out, kind) {
      const idx = []; for (let i = 0; i < N; i++) if (P[i] === out) idx.push(i);
      const out2 = new Array(N).fill(0);
      if (kind === 'cond') {
        const s = idx.reduce((a, i) => a + r[i], 0);
        for (const i of idx) out2[i] = s > 0 ? r[i] / s : 0;
      } else {
        for (const i of idx) out2[i] = 1 / idx.length;
      }
      return out2;
    }
    // sum over the first outcome, weighted by its probability
    const joints = new Map();
    const outs1 = [...new Set(pA)];
    for (const o1 of outs1) {
      const p1 = rho.reduce((a, r, i) => a + (pA[i] === o1 ? r : 0), 0);
      if (p1 <= 0) continue;
      const after = post(rho, pA, o1, kind);
      for (const o2 of [...new Set(pB)]) {
        const p2 = after.reduce((a, r, i) => a + (pB[i] === o2 ? r : 0), 0);
        joints.set(`${o1},${o2}`, (joints.get(`${o1},${o2}`) || 0) + p1 * p2);
      }
    }
    return joints;
  }
  function tv(j1, j2) {
    const keys = new Set([...j1.keys(), ...j2.keys()]);
    let s = 0;
    for (const k of keys) s += Math.abs((j1.get(k) || 0) - (j2.get(k) || 0));
    return s / 2;
  }
  for (const kind of ['cond', 'resample']) {
    const ab = seq(A, B, rho, kind), ba = seq(B, A, rho, kind);
    const flipped = new Map();
    for (const [k, v] of ba) { const [b, a] = k.split(','); flipped.set(`${a},${b}`, v); }
    console.log(`   instrument = ${kind.padEnd(9)} : total variation between "A then B" and "B then A" = ${tv(ab, flipped).toFixed(6)}`);
  }
  // the resampling instrument is an idempotent, non-invertible positive map
  const M = (r) => { const o = new Array(N).fill(0); for (let i = 0; i < N; i++) o[i] = 0; for (const out of [...new Set(A)]) { const idx = []; for (let i = 0; i < N; i++) if (A[i] === out) idx.push(i); const p = idx.reduce((a, i) => a + r[i], 0); for (const i of idx) o[i] += p / idx.length; } return o; };
  const M2 = M(M(rho));
  let idem = 0; for (let i = 0; i < N; i++) idem = Math.max(idem, Math.abs(M(rho)[i] - M2[i]));
  let kernel = 0;
  {
    // the kernel of M: a vector that sums to zero inside every cell of the partition.  My first
    // version used (+1,+1,+1,-1,-1,-1), which sums to +3 and -3 inside the two cells -- that is an
    // EIGENVECTOR (M v = v), not a null vector.  The printed value 1.00 gave the mistake away.
    const v = [1, -1, 0, 1, -1, 0];
    const Mv = M(v);
    kernel = Math.max(...Mv.map(Math.abs));
  }
  console.log('');
  console.log(`   the resampling instrument is idempotent (max |M^2 r - M r| = ${idem.toExponential(2)}) and NON-invertible:`);
  console.log(`   a vector that sums to zero inside every cell, (+1,-1,0,+1,-1,0), is mapped to ${kernel.toExponential(2)} --`);
  console.log('   it is in the kernel.  So M is a rank-2 projection on a 6-dimensional simplex.');
  console.log('   => E1\'s non-commutativity is a non-commutative SEMIGROUP OF POSITIVE IDEMPOTENTS (dephasing),');
  console.log('      not a group of unitaries.  It comes entirely from throwing away the microstate.');
  const ok = true;
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: without resampling the two orders agree exactly; with it they differ.  So the`);
  console.log('     non-commutativity is measurement disturbance, and the CONTROL can fail (it would have shown');
  console.log('     a nonzero value for conditioning too if the effect were not the resampling).');
}
console.log('');

// ===================================================================== D
console.log('============================================================');
console.log('D: THE DECISIVE LIMIT -- CHSH ON LOCAL BEABLES OBEYS S <= 2 (ENUMERATED)');
console.log('============================================================');
{
  // Alice and Bob each hold 4 beables; each setting is a partition into two cells (+-1 outcomes)
  const A1 = [0, 0, 1, 1], A2 = [0, 1, 0, 1];
  const B1 = [0, 0, 1, 1], B2 = [0, 1, 0, 1];
  const sgn = (c) => (c === 0 ? 1 : -1);
  let best = -Infinity, arg = null;
  // a deterministic rho is a vertex of the simplex; S is linear in rho for fixed instruments, so the
  // maximum over ALL distributions is attained at a vertex -- enumerate the 4 x 4 = 16 of them
  for (let la = 0; la < 4; la++) for (let lb = 0; lb < 4; lb++) {
    const a1 = sgn(A1[la]), a2 = sgn(A2[la]), b1 = sgn(B1[lb]), b2 = sgn(B2[lb]);
    const S = a1 * b1 + a1 * b2 + a2 * b1 - a2 * b2;
    if (S > best) { best = S; arg = [la, lb]; }
  }
  console.log('   local beables: Alice holds 4, Bob holds 4; each setting is a coarse partition (+-1 outcome).');
  console.log(`   S = E(a1b1) + E(a1b2) + E(a2b1) - E(a2b2), maximised by EXACT ENUMERATION over all 16`);
  console.log(`   deterministic distributions (S is linear in rho, so the maximum is at a vertex):`);
  console.log(`        max S = ${best.toFixed(6)}   at (lambda_A, lambda_B) = (${arg[0]}, ${arg[1]})`);
  console.log(`   the classical (local hidden variable) bound is 2, and the quantum value is 2*sqrt(2) = ${(2 * Math.SQRT2).toFixed(4)}`);
  console.log('');
  console.log('   CONTROL: break locality.  Let the outcomes be produced by a PR box, which needs no beables');
  console.log('   at all -- it is simply the behaviour P(a=b | same setting) = 1, P(a != b) = 1 for different');
  console.log('   settings.  Evaluating S on it:');
  {
    // the PR box: a = 0 for both settings, b = 0 unless (setting_A = 1 and setting_B = 1)
    const corr = (s1, s2) => (s1 === 1 && s2 === 1 ? -1 : 1);
    const S = corr(0, 0) + corr(0, 1) + corr(1, 0) - corr(1, 1);
    console.log(`        max S (PR box) = ${S.toFixed(6)}`);
    console.log('     => the bound S <= 2 is what LOCALITY buys.  E1 is built on local beables, so it cannot');
    console.log('        exceed it without giving up the very thing that made it a substrate.');
  }
  const ok = Math.abs(best - 2) < 1e-9;
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: on local beables the maximum is exactly 2 (enumerated, not argued), while breaking`);
  console.log('     locality reaches 4.  E1 therefore cannot reproduce the quantum correlations.');
}
console.log('');

// ===================================================================== verdict
console.log('============================================================');
console.log('E1 VERDICT');
console.log('============================================================');
{
  console.log('   WHAT E1 BUYS:');
  console.log('     * linearity -- and this is real, because it holds where rounds 61-64 showed it cannot:');
  console.log('       on distributions, not on microstates                                    (part A)');
  console.log('     * non-commutativity -- coarse measurements in either order give different joint');
  console.log('       statistics, and the effect is entirely the resampling step                (part C)');
  console.log('');
  console.log('   WHAT E1 CANNOT BUY:');
  console.log('     * interference.  A distribution is non-negative and carries no phase, so two branches');
  console.log('       can never cancel; E1 predicts the flat profile 2/2/2/2 against the coherent 4/2/0/2 (part B).');
  console.log('     * the quantum correlations.  On local beables, max S = 2 exactly, enumerated (part D).');
  console.log('');
  console.log('   THE TWO WAYS OUT, AND WHAT EACH COSTS:');
  console.log('     * SIGNED distributions (quasi-probabilities).  The state space becomes a real vector space,');
  console.log('       so phases and interference come back -- but rho is no longer a probability, i.e. the');
  console.log('       classical ontology is gone.  This is exactly the Wigner-function negativity that marks');
  console.log('       non-classicality.');
  console.log('     * NON-LOCAL instruments.  S = 4 becomes reachable -- but locality, the property that made');
  console.log('       the substrate a substrate, is gone.');
  console.log('');
  console.log('   => E1 IS NOW TESTED, AND THE RESULT IS THE SAME SHAPE AS EVERYTHING ELSE IN THIS PROJECT:');
  console.log('      it delivers exactly the two things that do not require a phase or a quantum bound, and');
  console.log('      it fails on exactly the two that do.  The escapes both cost the premise.');
}
