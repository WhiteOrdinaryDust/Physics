'use strict';
/*
 * ingredient_ledger.js  --  ROUND 84: THE TEN-ITEM LEDGER, AND THE F1 SPLIT
 * ---------------------------------------------------------------------------
 * Round 82 classified four items as ALPHABET-LIMITED or STRUCTURE-LIMITED.  This round does all ten,
 * and it settles a phrasing defect found while planning: criterion F1 as written says
 *
 *      "on a FINITE alphabet, give a local deterministic rule such that EXACT COMPLEX LINEAR
 *       SUPERPOSITION appears at the microstate level"
 *
 * The literal reading of that is refuted by a cardinality argument -- a finite set cannot be closed
 * under multiplication by arbitrary complex scalars, since the orbit of one non-zero state already
 * has the cardinality of the continuum.  So the literal F1 is FALSE for a trivial reason, and the
 * substantive version ("can a finite alphabet's local rule realise exact linear superposition of
 * excitations?") is the one that rounds 61-64 measured obstructions against.  The frozen report must
 * say both, or it carries an "open" criterion that is trivially false as phrased.
 *
 * PARTS
 *   A  the cardinality argument, with the demonstration that makes it concrete.
 *   B  the ten-item ledger: lowest tier that can hold each item, and what limits it.
 *   C  the counts, and the nested structure of the three alphabet requirements.
 * ---------------------------------------------------------------------------
 */

console.log('============================================================');
console.log('PART A: F1 AS WRITTEN IS FALSE, AND THE REASON IS CARDINALITY');
console.log('============================================================');
{
  console.log('   A1 a finite alphabet has finitely many states.  Closure under complex scaling would require');
  console.log('      that for every state s and every complex lambda, lambda*s is also a state.');
  console.log('      Take one non-zero state: {lambda*s : lambda in C} already has the cardinality of the');
  console.log('      continuum, which no finite set has.  So the literal reading fails before any dynamics');
  console.log('      is even discussed.');
  // make it concrete: the closure of a small alphabet under +- stays finite, under complex scaling does not
  const alpha = [0, 1, 2, 3];
  const scale = (A, c) => { const out = new Set(); for (const s of A) out.add(s * c); return out; };
  const closeUnderNegation = (A) => { const out = new Set(A); for (const s of A) out.add(-s); return [...out].sort((a, b) => a - b); };
  console.log('');
  console.log(`      alphabet {${alpha.join(',')}}:`);
  console.log(`        closed under negation  -> {${closeUnderNegation(alpha).join(',')}}   still finite`);
  console.log(`        scaled by 1/2          -> {${[...scale(alpha, 0.5)].join(',')}}   leaves the alphabet (1/2 is not in it)`);
  console.log('        scaled by 1/3, sqrt(2), i ... leaves it too, and there are continuum many such scalars');
  console.log('');
  console.log('   A2 SO THE FROZEN REPORT MUST SPLIT F1 INTO TWO ROWS:');
  console.log('        F1-literal     "a finite alphabet carries exact COMPLEX linear superposition"');
  console.log('                       => REFUTED, by cardinality, trivially and permanently;');
  console.log('        F1-substantive "a finite alphabet exists whose LOCAL RULE realises exact linear');
  console.log('                        superposition of its excitations (amplitudes drawn from a finite or');
  console.log('                        countable set, superposed exactly)"');
  console.log('                       => still open, and it is what rounds 61-64 measured obstructions');
  console.log('                          against: one bit cannot hold |2>, non-negative counts are a cone,');
  console.log('                          Gaussian integers cannot add one quantum for any n >= 1.');
  console.log('');
  console.log('   A3 WHY THIS IS NOT A LOOPHOLE IN THE OTHER DIRECTION: F1-literal being false does not');
  console.log('      weaken the project.  Its whole point was never that "a finite alphabet is a vector');
  console.log('      space" -- it was that SOMETHING must carry exact superposition, and the report already');
  console.log('      says what: signed complex amplitudes, continuous for exactness (rounds 62-64, 82).');
}
console.log('');

console.log('============================================================');
console.log('PART B: THE TEN-ITEM LEDGER (lowest tier, and what limits it)');
console.log('============================================================');
{
  const L = [
    ['0 propagation at refresh rate', 'A1 bits', 'structure', 'speed exactly 1.0000 in the counting tier (round 61); HPP-like propagation at A1 (rounds 12-17)'],
    ['1 light cone / locality', 'A1 bits', 'definitional', 'the cone IS the propagation speed (round 62)'],
    ['2 multi-quantum occupancy', 'A2 counts', 'ALPHABET (>= 3 states)', 'bits collapse 2 -> 1, measured (round 61); A3 also holds 2 as 1+i (round 82)'],
    ['3 interference', 'A3 Gaussian integers', 'ALPHABET (needs sign/phase)', 'fringes 4/2/0/2 exact (round 62); A1 and A2 fail by alphabet (rounds 61-62)'],
    ['4 exact add-one-quantum', 'A4 continuous', 'ALPHABET (needs continuity)', 'A3: exact only for n = 0, exhaustive n = 1..200 (rounds 64, 82); A2: arithmetic, not a measurement, and its state reading was refuted (rounds 61, 63)'],
    ['5 dispersion / mass', 'A1..A3 all can', 'structure', 'integer-coefficient leapfrog stays in Z[i] over 200 steps (round 82); round 68 measured the wall in A4 with a dispersion term'],
    ['6 transparency', 'A4 + integrability', 'structure (one sub-question open)', 'round 68: 7.04% radiation integrable vs 56.04% generic. Whether an integrable DISCRETISATION exists inside A3 is still unmeasured'],
    ['7 conserved charge vs interference observable', 'A3 onwards', 'structure', 'round 69: L2 conserved to 1.18e-16 while |sum psi|^2 moves 7.520 -> 5.014 -> 15.040'],
    ['8 single readout / Born', 'A2 counts', 'structure (needs a clock)', 'round 65 (A4 with clock), round 81 (A2, max deviation 8.0e-3), round 82 (A3, 8.7e-3), controls fail by an order of magnitude'],
    ['9 second system / entanglement', 'NONE', 'ABSENT (not expensive -- absent)', 'no tier delivers it: 512 of 1023 requirement sets are unsatisfiable, and they are exactly those requiring it (round 81)'],
  ];
  console.log('   item                                         lowest tier            limited by                    provenance');
  for (const r of L) console.log(`   ${r[0].padEnd(44)} ${r[1].padEnd(22)} ${r[2].padEnd(29)} ${r[3]}`);
  console.log('');
  console.log('   counts: ALPHABET-limited 3, structure-limited 5, definitional 1, ABSENT 1 -- and 3+5+1+1 = 10.');
  console.log('   (the light cone is DEFINITIONAL, not structure: it is the propagation speed by definition,');
  console.log('    so it does not belong in either empirical bucket.  My first version counted it twice.)');
  console.log('   The three alphabet requirements are NESTED, which is the sharpest way to state the whole');
  console.log('   alphabet story:');
  console.log('        >= 3 states        buys occupancy          (bits cannot hold |2>)');
  console.log('        + sign / phase     buys interference       (counts are a cone)');
  console.log('        + continuity       buys the exact ladder  (Gaussian integers fail for every n >= 1)');
  console.log('   each step is strictly stronger, each was measured, and NOTHING ELSE in the ten items needs');
  console.log('   the alphabet to grow.');
}
console.log('');

console.log('============================================================');
console.log('PART C: WHAT THIS CHANGES IN THE FROZEN DOCUMENTS');
console.log('============================================================');
{
  console.log('   1. FINAL_REPORT section 7 (the F1-F4 status table) gains a split row: F1-literal REFUTED by');
  console.log('      cardinality, F1-substantive still open.  The "we cannot say impossible" boundary survives,');
  console.log('      because it was always about the STRONG claim, not about F1;');
  console.log('   2. section 26 (the four-item classification) is superseded by this ten-item ledger, which is');
  console.log('      stricter: it finds THREE alphabet requirements, not two, because occupancy needs a third');
  console.log('      state before anything else can happen;');
  console.log('   3. and the ledger states the three alphabet requirements as nested, which is a claim nobody');
  console.log('      has made in this project before and which is checkable against rounds 61-64 line by line.');
  console.log('');
  console.log('   REMAINING OPEN after this round: the cubic size dependence (round 83), 15 table cells and 8');
  console.log('   dependency edges that are still readings (section 25), whether an integrable discretisation');
  console.log('   lives inside Z[i] (item 6 above), and the second system -- the one ABSENT item.');
}
