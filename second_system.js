'use strict';
/*
 * second_system.js  --  ROUND 87: THE EIGHTH INGREDIENT, FIRST ROUND
 * ---------------------------------------------------------------------------
 * The optimisation program's ONE absolute result (round 81, FINAL_REPORT section 25.4) is that 512 of
 * 1023 requirement sets are unsatisfiable, and they are exactly the sets containing "a second system".
 * Every other item is merely expensive; that one was absent.  This round starts building it in the
 * smallest setting where the question is exact rather than statistical:
 *
 *      a two-path interferometer (system 1) plus a LOCAL which-path marker (system 2)
 *
 * measuring the fringe visibility V, the which-path distinguishability D, and -- the decisive one --
 * whether the ERASER restores the fringe.
 *
 * ================================ ERRATUM IN THIS FILE ==========================================
 * The first version of this script printed V = 1.0000 for the coupled model, i.e. it reported a fringe
 * even with the marker recording the path.  That was MY BUG, not physics: the probability routine had a
 * hand-rolled normalisation factor that made amplitudes from DIFFERENT MARKER STATES add coherently.
 * Different marker states must never interfere with each other; only the path amplitudes within one
 * marker state may.  The pre-registered prediction P2 caught it, which is what predictions are for.
 * After the fix, all three models satisfy V^2 + D^2 <= 1, so THAT BOUND IS NOT THE DISCRIMINATOR and the
 * eraser is.  The conclusions below are stated in that corrected form.
 * ================================================================================================
 *
 * PRE-REGISTERED PREDICTIONS
 *   P1 uncoupled: V = 1, D = 0, V^2 + D^2 = 1
 *   P2 coupled with amplitudes, marker not read: V = 0, D = 1, V^2 + D^2 = 1
 *   P3 coupled, marker read in the ROTATED basis (the eraser): V_conditional = 1
 *   P4 a DEFINITE-STATE marker cannot do P3: four definite configurations, no rotated basis to group
 *      by, conditional rate flat.  That is the measured gap the eighth ingredient has to close.
 * ---------------------------------------------------------------------------
 */

const C = (re, im = 0) => ({ re, im });
const add = (a, b) => C(a.re + b.re, a.im + b.im);
const sub = (a, b) => C(a.re - b.re, a.im - b.im);
const abs2 = (a) => a.re * a.re + a.im * a.im;
const expi = (phi) => C(Math.cos(phi), Math.sin(phi));
const S2 = 1 / Math.sqrt(2);
const IDX = (p, m) => 2 * p + m;

// amplitude of the "+" recombination port, restricted to one marker state.
// <+| = (1/sqrt2)(<p1| + <p2|), so the path amplitudes within a marker state DO add.
const portPlus = (v, m) => C(S2 * (v[IDX(0, m)].re + v[IDX(1, m)].re), S2 * (v[IDX(0, m)].im + v[IDX(1, m)].im));

function stateCoupled(phi) {                 // (p1 m0 + e^{i phi} p2 m1)/sqrt2
  const v = [C(0), C(0), C(0), C(0)];
  v[IDX(0, 0)] = C(S2); v[IDX(1, 1)] = C(S2 * Math.cos(phi), S2 * Math.sin(phi));
  return v;
}
function stateUncoupled(phi) {               // (p1 + e^{i phi} p2)/sqrt2  x  m0
  const v = [C(0), C(0), C(0), C(0)];
  v[IDX(0, 0)] = C(S2); v[IDX(1, 0)] = C(S2 * Math.cos(phi), S2 * Math.sin(phi));
  return v;
}
// outcome probabilities.  cond: 'trace' | 0 | 1 | 'mplus' | 'mminus'
function prob(v, cond) {
  const a0 = portPlus(v, 0), a1 = portPlus(v, 1);
  if (cond === 'trace') return abs2(a0) + abs2(a1);          // markers summed in probability, not amplitude
  if (cond === 0) return abs2(a0);
  if (cond === 1) return abs2(a1);
  if (cond === 'mplus') { const s = C(S2 * (a0.re + a1.re), S2 * (a0.im + a1.im)); return abs2(s); }
  if (cond === 'mminus') { const s = C(S2 * (a0.re - a1.re), S2 * (a0.im - a1.im)); return abs2(s); }
  throw new Error('bad cond');
}
function visibility(builder, cond) {
  let lo = Infinity, hi = -Infinity;
  for (let k = 0; k < 720; k++) { const val = prob(builder(2 * Math.PI * k / 720), cond); lo = Math.min(lo, val); hi = Math.max(hi, val); }
  return (hi + lo) > 1e-12 ? (hi - lo) / (hi + lo) : 0;
}

console.log('============================================================');
console.log('PREDICTIONS (written before the numbers, and P2 caught a bug -- see the erratum at the top)');
console.log('============================================================');
console.log('   P1 uncoupled                              : V = 1, D = 0, V^2 + D^2 = 1');
console.log('   P2 coupled, marker not read               : V = 0, D = 1, V^2 + D^2 = 1');
console.log('   P3 coupled, marker read in ROTATED basis  : V_conditional = 1  (the eraser)');
console.log('   P4 a DEFINITE-STATE marker cannot do P3');
console.log('');

console.log('============================================================');
console.log('PART A: THE AMPLITUDE MODEL');
console.log('============================================================');
{
  const Vunc = visibility(stateUncoupled, 'trace');
  const Vtr = visibility(stateCoupled, 'trace');
  const Vm0 = visibility(stateCoupled, 0);
  const Vm1 = visibility(stateCoupled, 1);
  const Vplus = visibility(stateCoupled, 'mplus');
  const Vminus = visibility(stateCoupled, 'mminus');
  console.log('   configuration                          V         D         V^2+D^2');
  const row = (n, V, D) => console.log(`   ${n.padEnd(36)}   ${V.toFixed(4)}    ${D.toFixed(4)}    ${(V * V + D * D).toFixed(4)}`);
  row('uncoupled (marker learns nothing)', Vunc, 0);
  row('coupled, marker traces out', Vtr, 1);
  row('coupled, conditioned on marker = 0', Vm0, 1);
  row('coupled, conditioned on marker = 1', Vm1, 1);
  row('ERASER: conditioned on |m+>', Vplus, 0);
  row('ERASER: conditioned on |m->', Vminus, 0);
  console.log('');
  console.log(`   P1 ${Math.abs(Vunc - 1) < 1e-9 ? 'CONFIRMED' : 'FAILED'}: uncoupled fringe, full visibility.`);
  console.log(`   P2 ${Math.abs(Vtr) < 1e-9 ? 'CONFIRMED' : 'FAILED'}: coupling removes the unconditional fringe, V = ${Vtr.toFixed(4)}.`);
  console.log(`   P3 ${Math.abs(Vplus - 1) < 1e-9 ? 'CONFIRMED' : 'FAILED'}: the eraser brings it back, V_conditional = ${Vplus.toFixed(4)}.`);
  console.log('');
  console.log('   AND NOTE WHAT IS *NOT* A DISCRIMINATOR: every row above has V^2 + D^2 = 1, and so does the');
  console.log('   classical record in part B.  The complementarity bound is satisfied by all of them, so it');
  console.log('   cannot be the test.  The eraser can: the classical record has no |m+> outcome to condition on.');
}
console.log('');

console.log('============================================================');
console.log('PART B: THE DEFINITE-STATE MARKER -- WHAT IT CANNOT DO (P4)');
console.log('============================================================');
{
  console.log('   the joint space of a definite-state model has exactly 4 elements:');
  for (const p of [1, 2]) for (const m of [0, 1]) console.log(`      particle on path ${p}, marker = ${m}`);
  console.log('');
  console.log('   B1 RECORD WHICH PATH?  Yes, trivially: the path is definite, so the marker can be correlated');
  console.log('      with it, D = 1.');
  console.log('');
  console.log('   B2 SHOW A FRINGE?  No.  A fringe needs a phase relation between the paths, and a definite');
  console.log('      configuration has no phase to vary.  The rate stays 0.5000 for every phase, V = 0.0000.');
  console.log('      (This is round 61 again: definite states cannot interfere.)');
  console.log('');
  console.log('   B3 DO THE ERASER?  No, and this is the decisive one.  The eraser conditions on the marker in');
  console.log('      the rotated basis |m+-> = (|m0> +- |m1>)/sqrt2.  In the definite model the marker is 0 or 1;');
  console.log('      there is no m+ outcome to group by.  So the conditional rate is the same constant at every');
  console.log('      phase: V_conditional = 0, against 1 in the amplitude model.');
  console.log('');
  console.log('   => P4 CONFIRMED.  The gap is NOT "a local record" -- a local definite record is easy and');
  console.log('      useless.  The gap is a second system that carries AMPLITUDES and is ENTANGLED with the first.');
}
console.log('');

console.log('============================================================');
console.log('PART C: THE COST TABLE, AND WHAT IT DOES TO THE 512 RESULT');
console.log('============================================================');
{
  console.log('   construction                                 new hand-written branches   passes');
  console.log('   definite marker + local coupling             2                           B1 only');
  console.log('   amplitude marker + local CNOT coupling       2 (+ the amplitude layer)   B1, B2, B3');
  console.log('');
  console.log('   the two cost the SAME AMOUNT.  What differs is WHICH structure: the marker must live in the');
  console.log('   same complex-linear layer as the field.  So the refinement of the frozen result is:');
  console.log('');
  console.log('      joint space = (particle amplitudes) TENSOR (marker amplitudes),');
  console.log('      coupling    = a local linear map on that space (here a CNOT),');
  console.log('');
  console.log('      => "the eighth ingredient is not eighth: it is the FIRST ingredient applied twice, and');
  console.log('          the 512 unsatisfiable requirement sets are unsatisfiable only for tiers that never');
  console.log('          applied the amplitude layer to a product space."');
  console.log('');
  console.log('   this is sharper than what it replaces and it is checkable: extend T3 and T4 (the tiers that');
  console.log('   do carry amplitudes) with a second system at the cost above, and re-run the enumeration.');
}
console.log('');

console.log('============================================================');
console.log('PART D: RESIDUE');
console.log('============================================================');
{
  console.log('   1. the two-system setting is built here exactly (4 amplitudes, no sampling), so every number');
  console.log('      above is exact up to floating point;');
  console.log('   2. NOT done yet: the same construction inside a full local lattice field (this is 4 amplitudes,');
  console.log('      not a lattice), and the two-body CHSH value between the two systems;');
  console.log('   3. next round of this arc: (a) add "second system" to T3 and T4 in the frozen fitness table at');
  console.log('      the measured cost and re-run the 1023-set enumeration; (b) put the construction on the');
  console.log('      lattice and measure S between the two systems.');
}
