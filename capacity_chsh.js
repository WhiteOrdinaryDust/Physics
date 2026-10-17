'use strict';
/*
 * capacity_chsh.js
 * ---------------------------------------------------------------------------
 * STEP: derive what a "shortcut capacity" constraint implies for CHSH.
 *
 * BACKGROUND (established earlier with shortcut_chsh.js, exhaustive over 16.7M
 * classical models):
 *   - a shortcut made of classical shared randomness caps CHSH at 2
 *   - reaching 2*sqrt(2) requires complex amplitudes, i.e. a two-level quantum
 *     system shared between the parties
 *   - NO-SIGNALLING ALONE DOES NOT CAP IT: the PR box is no-signalling with
 *     CHSH = 4. So a further principle is required.
 *
 * THE QUESTION HERE
 *   Candidate principle: the shortcut has FINITE CAPACITY. It can only carry a
 *   bounded amount of correlation. Does some natural normalisation of "capacity"
 *   pin the maximum at exactly 2*sqrt(2), rather than at 2 or at 4?
 *
 * TOOL
 *   Horodecki criterion: for a two-qubit state with correlation matrix T, the
 *   largest CHSH value is
 *        S_max = 2 * sqrt( t1^2 + t2^2 )
 *   where t1 >= t2 are the two largest singular values of T.
 *   For the pure state cos(a)|00> + sin(a)|11> one finds
 *        t1 = 2 cos a sin a = sin 2a,  t2 = sin 2a
 *   so S_max = 2 * sqrt(2) * sin(2a), maximised (= 2*sqrt(2)) at a = pi/4.
 *
 *   The capacity parameter enters by deciding WHICH states are allowed. We test
 *   several candidate normalisations below.
 * ---------------------------------------------------------------------------
 */

// ---- the pure state cos(a)|00> + sin(a)|11>
// Concurrence C = 2 cos a sin a = sin 2a. Schmidt rank 2 for 0 < a < pi/2.
function tValues(a) {
  const C = Math.sin(2 * a);          // concurrence
  return { t1: C, t2: C, concurrence: C };
}
function sMax(a) {
  const { t1, t2 } = tValues(a);
  return 2 * Math.sqrt(t1 * t1 + t2 * t2);
}

console.log('=== CHSH maximum for the two-qubit family cos(a)|00> + sin(a)|11> ===');
console.log('   a(deg)   concurrence   S_max      S_max/2sqrt2');
for (let deg = 5; deg <= 45; deg += 5) {
  const a = deg * Math.PI / 180;
  const C = Math.sin(2 * a);
  console.log(`  ${String(deg).padStart(6)}   ${C.toFixed(4).padStart(10)}   ${sMax(a).toFixed(4).padStart(7)}   ${(sMax(a) / (2 * Math.SQRT2)).toFixed(4).padStart(11)}`);
}
console.log(`\n  maximum over a: ${Math.max(...Array.from({length:91},(_,i)=>sMax((i+1)*Math.PI/180))).toFixed(6)}`);
console.log(`  2*sqrt(2)     : ${(2 * Math.SQRT2).toFixed(6)}`);

/*
 * Now the substantive part: WHICH capacity normalisation is physical?
 *
 * A "shortcut" between nodes A and B is a coupling. Its capacity should measure
 * how much correlation it can carry. Three natural candidates:
 *
 *   N1. CAPACITY = SCHMIDT RANK ONLY.
 *       "The shortcut can carry at most one bit of shared information."
 *       Allowed: all states of Schmidt rank <= 2 (any a).
 *       => S_max reaches 2*sqrt(2). Too permissive: it also allows a = pi/4 for
 *          free, so it caps nothing below Tsirelson. USELESS as a cap.
 *
 *   N2. CAPACITY = CONCURRENCE (or entanglement entropy), bounded by the flow
 *       through a link of bounded degree.
 *       If a link can carry at most concurrence c, then sin(2a) <= c and
 *          S_max(c) = 2*sqrt(2)*c.
 *       This is a straight line from 2 (at c = 1/sqrt2) to 2sqrt2 (at c = 1).
 *       => a knob, not a derivation. Any value in [2, 2sqrt2] is reachable.
 *
 *   N3. CAPACITY = A GLOBAL CONSERVATION LAW (monogamy), summed over partners:
 *          sum_B C(A:B)^2 <= 1     (the CKW constraint for a qubit A)
 *       With one partner this allows C = 1, hence 2sqrt2. With two partners it
 *       limits the split, but does not by itself lower the single-partner max.
 *
 * The honest conclusion is drawn in the output below.
 */

console.log('\n=== test of the candidate capacity normalisations ===');

console.log('\nN1. capacity = Schmidt rank <= 2 (one bit of shared information)');
console.log('    allowed concurrence range: (0, 1]');
console.log(`    implied max CHSH: ${sMax(Math.PI / 4).toFixed(4)}  = 2*sqrt(2)`);
console.log('    -> CAPS NOTHING. Rank alone permits the maximally entangled state,');
console.log('       so it cannot explain why 2*sqrt(2) is the ceiling rather than 4.');

console.log('\nN2. capacity = concurrence bounded by c');
console.log('    implied max CHSH = 2*sqrt(2)*c');
for (const c of [0.5, 0.7071, 0.8, 0.9, 1.0]) {
  console.log(`      c=${c.toFixed(4)} -> S_max=${(2 * Math.SQRT2 * c).toFixed(4)}`);
}
console.log('    -> A FREE KNOB. Every value in [2, 2sqrt2] is reachable by choosing c.');
console.log('       This reproduces the spectrum but does not pin the ceiling.');

console.log('\nN3. capacity = monogamy sum over partners');
console.log('    sum_B C(A:B)^2 <= 1');
console.log('    one partner: C=1 allowed -> S_max = 2sqrt2 (no further restriction)');
console.log('    two partners: C(A:B)^2 + C(A:C)^2 <= 1 -> each partner gets less');
console.log('    -> constrains DISTRIBUTION across partners, not the two-party ceiling.');

console.log('\n=== the decisive numerical question ===');
console.log('  Is there any capacity functional that (a) is bounded by the local');
console.log('  structure, and (b) forces S_max <= 2sqrt2 with equality attainable?');
console.log('');
console.log('  S_max = 2*sqrt(t1^2+t2^2) and for pure two-qubit states t1=t2=C.');
console.log('  So S_max = 2sqrt2 * C. Therefore:');
console.log('     S_max <= 2sqrt2   <=>   C <= 1');
console.log('  and C <= 1 is automatic for ANY two-qubit state (concurrence is');
console.log('  normalised to at most 1).');
console.log('');
console.log('  ==> S_max <= 2sqrt2 is therefore EQUIVALENT to C <= 1, and for a');
console.log('      two-qubit state C <= 1 is automatic -- it is the NORMALISATION of');
console.log('      concurrence, not a physical restriction.');
console.log('');
console.log('  So this is not an explanation: it restates "the resource is a qubit".');
console.log('  The real content would have to be a substrate argument for WHY the');
console.log('  resource is two-level. See the dimension control below before drawing');
console.log('  any conclusion from this line of reasoning.');

// ---- does a higher-dimensional resource beat 2sqrt2?
console.log('\n=== control: dimension and the Tsirelson bound ===');
console.log('  CAREFUL -- the obvious guess is wrong.');
console.log('  For CHSH specifically, the Tsirelson bound 2*sqrt(2) holds for ALL');
console.log('  dimensions: a d-level (qudit) resource does NOT beat 2*sqrt(2) on CHSH.');
console.log('  So on this inequality, "capacity = 1 qubit" is NOT what produces the cap.');
console.log('');
console.log('  Dimension dependence shows up in OTHER Bell inequalities, where the');
console.log('  maximum violation does grow with d (e.g. the CGLMP family). I am not');
console.log('  quoting specific CGLMP values here because I have not verified them.');
console.log('');
console.log('  ==> Consequence for the argument: the identification');
console.log('      "2sqrt2 ceiling  <=>  two-level shortcut"');
console.log('      is NOT established. Both halves are separately true');
console.log('      (S_max = 2sqrt2*C for two-qubit states; Tsirelson holds in all');
console.log('      dimensions) but they do not combine into the claimed equivalence.');
console.log('');
console.log('  What the two-qubit calculation DOES establish rigorously:');
console.log('    * S_max = 2*sqrt(2) * C  exactly, for the pure two-qubit family');
console.log('    * therefore a concurrence bound c gives S_max = 2sqrt(2)*c');
console.log('    * and "no cap below Tsirelson" is equivalent to "C = 1 is allowed"');
console.log('  That is a clean parametrisation, but it is a knob (N2), not a');
console.log('  derivation of the ceiling.');
