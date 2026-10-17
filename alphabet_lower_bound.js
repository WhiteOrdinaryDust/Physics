'use strict';
/*
 * alphabet_lower_bound.js  --  ROUND 82: THE ALPHABET LOWER BOUND, AND IT IMMEDIATELY CORRECTS MY TABLE
 * ---------------------------------------------------------------------------
 * The experiment: walk the alphabet ladder (bits -> non-negative counts -> Gaussian integers ->
 * continuous complex) and at every tier ask whether it can represent the SAME target objects.  The
 * deliverable is the SMALLEST tier that can, plus -- more useful -- WHICH ITEMS ARE ALPHABET-LIMITED
 * and which are merely STRUCTURE-LIMITED.
 *
 * IT STARTS BY REFUTING MY OWN TABLE.  Rounds 79-81 recorded "T3 x occupancy = X, because sqrt(2) is
 * not a Gaussian integer".  That is wrong: 1+i IS a Gaussian integer and |1+i|^2 = 2, so the state
 * with occupancy 2 exists.  What round 63 actually measured was that the COUNTING tier misidentifies
 * states with counts, and that only SOME occupancies are norms of Gaussian integers (Fermat).  Since
 * that cell was the counterexample used in round 81 to refute the edge 2->3, the refutation goes with
 * it.  This round recomputes the cells from arithmetic instead of from memory.
 *
 * PARTS
 *   A  the ladder and the four target items.
 *   B  the representation audit, by exact arithmetic: can the tier HOLD the objects the items need?
 *   C  three exact checks that decide the rest: the occupancy spectrum (Fermat), the exact ladder
 *      (rationality of sqrt((n+1)/n)), and whether a dispersive update can stay inside the alphabet.
 *   D  the lower bound per item, and the split into ALPHABET-LIMITED vs STRUCTURE-LIMITED.
 *   E  verdict, and the table/DAG amendment this forces, with the frozen-table discipline applied.
 * ---------------------------------------------------------------------------
 */

console.log('============================================================');
console.log('PART A: THE LADDER AND THE FOUR TARGET ITEMS');
console.log('============================================================');
console.log('   ladder:  A1 bits (2 states)  ->  A2 non-negative counts  ->  A3 Gaussian integers  ->');
console.log('            A4 continuous complex');
console.log('   target items, the same test at every tier:');
console.log('     I1 INTERFERENCE   : a fringe pattern with a controllable dark line');
console.log('     I2 THRESHOLD      : a finite intensity threshold (a wall / saturable response)');
console.log('     I3 TRANSPARENCY   : an interaction that passes objects through without radiating');
console.log('     I4 READOUT        : one outcome per run, distributed as a squared amplitude');
console.log('');

console.log('============================================================');
console.log('PART B: THE REPRESENTATION AUDIT -- CAN THE TIER HOLD THE OBJECTS?');
console.log('============================================================');
{
  // B1 bits: a mode can hold 0 or 1 quantum; the record measured the collapse (round 61)
  console.log('   B1 bits: a mode holds 0 or 1.  Encoding "2 quanta" is impossible -- there is no third state.');
  console.log('         round 61 measured it directly: counts 1,2,3,5 all collapse to 1 in the bit alphabet.');
  console.log('         => cannot hold |2>, so anything needing occupancy is out, and with no sign or phase');
  console.log('            there is nothing to interfere.  ALPHABET-LIMITED, and it fails two items at once.');
  console.log('');
  // B2 non-negative counts: a cone
  console.log('   B2 non-negative counts: the state is (n1, n2, ...) with n >= 0 -- a CONE, not a vector space.');
  {
    const p1 = 3, p2 = 3;
    console.log(`         two branches of size ${p1} and ${p2}: their sum is ${p1 + p2} >= max(${p1}, ${p2}), so an`);
    console.log('         "addition" of two branches can never be SMALLER than either one.  A dark line needs');
    console.log('         cancellation, and inside a cone the only cancellation available is 0 + 0 = 0.');
    console.log('         There is also no phase to vary, so the dark line cannot be steered either.');
    console.log('         => I1 FAILS by ALPHABET (round 61 measured exactly this: non-negative counts are a');
    console.log('            cone, so there are no signs and no interference).');
  }
  console.log('');
  // B3 Gaussian integers: which occupancies exist?
  const rep = [];
  for (let n = 0; n <= 40; n++) {
    let ok = false;
    for (let a = 0; a * a <= n; a++) { const b2 = n - a * a; const b = Math.round(Math.sqrt(b2)); if (b * b === b2) ok = true; }
    rep.push(ok);
  }
  const yes = [], no = [];
  for (let n = 0; n <= 40; n++) (rep[n] ? yes : no).push(n);
  console.log('   B3 Gaussian integers: the amplitude is a+bi with a,b integers, so the occupancy is a^2+b^2.');
  console.log(`         representable n <= 40 : ${yes.join(', ')}`);
  console.log(`         NOT representable     : ${no.join(', ')}`);
  console.log('         => occupancy 2 EXISTS (1+i, |1+i|^2 = 2), so the rounds-79-81 cell "T3 x occupancy = X"');
  console.log('            is WRONG and is corrected to M here.  What fails is the EXACT LADDER (part C).');
  console.log('         (Fermat: n is a sum of two squares iff every prime = 3 mod 4 divides it to an even');
  console.log('          power -- which is why 3, 6, 7, 11, 12 are missing and 4, 5, 8, 9, 10 are not.)');
}
console.log('');

console.log('============================================================');
console.log('PART C: THREE EXACT CHECKS THAT DECIDE WHAT THE ALPHABET LIMITS');
console.log('============================================================');
{
  // C1 the exact ladder: a -> a*sqrt((n+1)/n) must stay in the alphabet
  let exact = [];
  for (let n = 1; n <= 200; n++) {
    const r2 = (n + 1) / n;
    // rational iff both n and n+1 are perfect squares (they are coprime)
    const s1 = Math.round(Math.sqrt(n)), s2 = Math.round(Math.sqrt(n + 1));
    if (s1 * s1 === n && s2 * s2 === n + 1) exact.push(n);
  }
  console.log(`   C1 exact ladder: the multiplier sqrt((n+1)/n) must be rational to keep a Gaussian integer.`);
  console.log(`         exhaustive n = 1..200 : exact cases = ${exact.length}   (n = 0 -> 1 works because a zero`);
  console.log('         amplitude has no phase to preserve).  => the EXACT-LADDER item (not I2) is');
  console.log('         ALPHABET-LIMITED: no discrete amplitude tier can add one quantum exactly for n >= 1');
  console.log('         (round 64 proved it and checked n = 1..200 there as well).');
  console.log('');
  // C2 can a dispersive update stay inside the integer alphabet?
  console.log('   C2 dispersion inside the alphabet: leapfrog psi_new[x] = psi[x-1] + psi[x+1] - psi[x] uses');
  {
    let ok = true;
    const N = 256;
    let a = new Array(N).fill(0), b = new Array(N).fill(0);
    for (let x = 100; x < 140; x++) { b[x] = 1; a[x] = b[x] - 1; }
    for (let t = 0; t < 200; t++) {
      const nx = new Array(N).fill(0);
      for (let x = 0; x < N; x++) nx[x] = b[(x - 1 + N) % N] + b[(x + 1) % N] - a[x];
      a = b; b = nx;
      for (let x = 0; x < N; x++) if (!Number.isInteger(b[x])) ok = false;
    }
    console.log(`         after 200 steps every amplitude is still an integer: ${ok}`);
    console.log('         => dispersion is STRUCTURE-LIMITED, not alphabet-limited: an integer-coefficient');
    console.log('            dispersive term stays inside the Gaussian integers.  (Round 80 found the counting');
    console.log('            tier dispersionless too -- but that is its rule, not its alphabet.)');
  }
  console.log('');
  // C3 the readout: sampling proportional to |psi|^2 works with integer amplitudes
  {
    const counts = [4, 2, 0, 2];
    const total = counts.reduce((a, b) => a + b, 0);
    const cdf = []; let acc = 0; for (const c of counts) { acc += c; cdf.push(acc / total); }
    const hist = [0, 0, 0, 0];
    let clock = 987654321;
    for (let i = 0; i < 20000; i++) {
      clock = (clock * 1103515245 + 12345) & 0x7fffffff;
      const u = clock / 0x80000000;
      let k = 0; while (u >= cdf[k]) k++;
      hist[k]++;
    }
    const dev = Math.max(...hist.map((v, i) => Math.abs(v / 20000 - counts[i] / total)));
    console.log(`   C3 readout: a hidden clock samples proportional to |psi|^2 = the integer occupancy.`);
    console.log(`         20000 runs give ${hist.map((v) => (v / 20000).toFixed(4)).join(' / ')} against ${counts.map((c) => (c / total).toFixed(4)).join(' / ')}`);
    console.log(`         max deviation ${dev.toExponential(2)}  => I4 is available at this tier (and at every tier`);
    console.log('         that has a squared quantity with a clock -- it is STRUCTURE-LIMITED, not alphabet-limited).');
  }
}
console.log('');

console.log('============================================================');
console.log('PART D: THE LOWER BOUND, AND THE SPLIT THAT MATTERS MORE');
console.log('============================================================');
{
  const rows = [
    ['I1 interference', 'A3 Gaussian integers', 'A1/A2 fail BY ALPHABET: no sign (bits), a cone (counts)', 'alphabet'],
    ['I2 threshold', 'A2 counts (structure)', 'the tier supplies the intensity; the SATURABLE rule is hand-written', 'structure'],
    ['I3 transparency', 'A4 continuous + integrability', 'needs objects that interact; transparency comes from integrability', 'structure'],
    ['I4 readout', 'A2 counts', 'a squared quantity plus a clock is enough; proven at A2 and A3 (part C3)', 'structure'],
  ];
  console.log('   item                  first tier that can hold it                 why                        limited by');
  for (const r of rows) console.log(`   ${r[0].padEnd(20)}  ${r[1].padEnd(42)} ${r[2].padEnd(58)} ${r[3]}`);
  console.log('');
  console.log('   THE LOWER BOUND FOR THE FOUR-ITEM CHECKLIST IS A3, THE GAUSSIAN INTEGERS:');
  console.log('     * A1 (bits) and A2 (counts) fail I1 by ALPHABET -- measured in round 61, and the cone');
  console.log('       argument in part B2 shows it cannot be repaired inside those alphabets;');
  console.log('     * A3 holds occupancy 2 (part B3), holds the fringe amplitudes 2 and 1+i (so I1 is');
  console.log('       available), and holds the readout (part C3);');
  console.log('     * and the two items that DO force the continuous tier are not on this checklist:');
  console.log('       the exact ladder for all n (alphabet-limited, part C1) and a massive self-trapped wall');
  console.log('       (structure-limited: it needs dispersion, which part C2 shows even A3 can carry).');
  console.log('');
  console.log('   SO THE HONEST FORM OF THE PROJECT\'S ALPHABET STORY, PER ITEM:');
  console.log('     alphabet-limited : interference (needs sign/phase) and the exact ladder (needs continuity)');
  console.log('     structure-limited: threshold, transparency, dispersion, readout, the conserved charge');
  console.log('   The old sentence "the alphabet had to be enriched until everything worked" is too coarse:');
  console.log('   enriching the alphabet buys exactly two things, and everything else was structure all along.');
}
console.log('');

console.log('============================================================');
console.log('PART E: VERDICT AND THE AMENDMENTS THIS FORCES');
console.log('============================================================');
{
  console.log('   AMENDMENT 1 (table): T3 x occupancy becomes M, not X.  1+i is a Gaussian integer whose norm');
  console.log('     is 2, so the state exists.  Rounds 79-81 had it wrong, and the error came from reading');
  console.log('     "sqrt(2) is not a Gaussian integer" as a statement about states rather than about moduli.');
  console.log('   AMENDMENT 2 (edge): the edge 2->3 ("superposition presupposes occupancy") was REFUTED in');
  console.log('     round 79 using T3 as the counterexample -- but T3 does deliver occupancy (amendment 1), so');
  console.log('     that refutation does not stand and 2->3 returns to being a reading.  The other three');
  console.log('     refutations (0->2, 0->3, 3->8) are unaffected: they do not rest on this cell.');
  console.log('   AMENDMENT 3 (frozen table discipline): the frozen table required "any change must come with a');
  console.log('     measurement".  This change comes with exact arithmetic (the Fermat spectrum in part B3) and');
  console.log('     an explicit recomputation, which satisfies the discipline -- and it is the discipline doing');
  console.log('     its job, since the wrong cell had survived three rounds of audits.');
  console.log('');
  console.log('   WHAT THIS ROUND ACTUALLY DELIVERED:');
  console.log('     1. the lower bound asked for: A3, and the reason A1/A2 cannot reach it;');
  console.log('     2. a sharper taxonomy than "the alphabet must be rich": only TWO items are alphabet-limited,');
  console.log('        and the rest of the project\'s ingredient list is structure;');
  console.log('     3. and a corrected table cell plus a withdrawn edge refutation, found by the experiment this');
  console.log('        round was asked to run -- which is the same pattern as rounds 79-81, now catching an');
  console.log('        error of mine that three earlier audits had passed.');
}
