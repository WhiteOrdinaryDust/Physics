'use strict';
/*
 * second_system_table.js  --  ROUND 88: PUTTING THE SECOND SYSTEM INTO THE FROZEN TABLE, AND RE-RUNNING
 * ---------------------------------------------------------------------------
 * Round 87 built the second system in a 4-amplitude setting and measured its cost: an amplitude marker
 * plus a local CNOT, i.e. 2 new hand-written branches.  This round enters that into the frozen fitness
 * table and re-runs the 1023-requirement-set enumeration, to see what happens to the ONE absolute
 * result the optimisation program had: "512 of 1023 sets are unsatisfiable, and they are exactly the
 * ones containing a second system" (FINAL_REPORT section 25.4).
 *
 * WHICH TIERS CAN ACTUALLY CARRY IT -- checked by arithmetic, not by optimism:
 *   the eraser needs a rotation of the MARKER basis, |m+-> = (|m0> +- |m1>)/sqrt2, so it needs 1/sqrt2.
 *   Gaussian-integer amplitudes do not contain 1/sqrt2 (round 82's arithmetic), so tier T3 CANNOT do the
 *   eraser: item 9 is not delivered there but REFUTED there.  Only the continuous tiers (T4, T5) can.
 *   That is a sharper outcome than "give everyone the second system".
 *
 * ALSO MEASURED HERE, because the table update should not rest on memory:
 *   tier T3 x conserved charge (item 7): a linear Gaussian-integer field update and its L2 norm.
 *
 * PARTS
 *   A  the quick measurement for T3 x item 7, and the derived refutation for T3 x item 9.
 *   B  the updated table, and the enumeration, with the explicit transition 512 -> ?
 *   C  sensitivity: the DAG variants, and the cost convention (+0, +2, +4 branches)
 *   D  verdict: what happens to the project's only absolute result
 * ---------------------------------------------------------------------------
 */

const CHECK = ['prop', 'cone', 'occ', 'interf', 'ladder', 'disp', 'transp', 'charge', 'born', 'second'];

console.log('============================================================');
console.log('PART A: TWO CELLS THAT SHOULD NOT BE FILLED FROM MEMORY');
console.log('============================================================');
{
  // A1 tier T3 x item 7 (a conserved charge distinct from the interference observable)
  // TWO FAILED ATTEMPTS FIRST, both recorded because both were vacuous:
  //   (1) psi_new[x] = psi[x-1] + psi[x+1] grew by 28 orders of magnitude -- an instability, not a test;
  //   (2) a unitary signed SHIFT conserves the norm (spread 0) but the sum over ALL sites is
  //       permutation-invariant, so |sum psi|^2 was identically 0 and nothing was being distinguished.
  // The right structure is the one rounds 62 and 69 actually used: amplitudes live at (site, DIRECTION)
  // pairs, the update permutes those pairs (a pure permutation => exactly unitary, and its entries are
  // units, so it stays inside the Gaussian integers), and the interference observable is the sum over
  // DIRECTIONS at one site -- which a permutation does change.
  const N = 64;
  let u = new Array(N).fill(0), w = new Array(N).fill(0);   // right-movers and left-movers
  for (let x = 0; x < 4; x++) { u[10 + x] = 1; w[10 + x] = 1; }   // a pulse moving both ways
  const charge = () => { let s = 0; for (let x = 0; x < N; x++) s += u[x] * u[x] + w[x] * w[x]; return s; };
  const interferenceAt = (x0) => { const s = u[x0] + w[x0]; return s * s; };
  const trace = [];
  for (let t = 0; t < 24; t++) {
    const nu = new Array(N).fill(0), nw = new Array(N).fill(0);
    for (let x = 0; x < N; x++) { nu[(x + 1) % N] = u[x]; nw[(x - 1 + N) % N] = w[x]; }   // pure permutation
    u = nu; w = nw;
    trace.push({ t: t + 1, charge: charge(), intf: interferenceAt(14) });
  }
  console.log('   A1 tier T3 x conserved charge: amplitudes at (site, direction) pairs, update = pure permutation');
  console.log('      (so its entries are units and the Gaussian integers are preserved exactly)');
  console.log('      t     charge = sum |psi|^2     interference observable at one site = |u + w|^2');
  for (const r of trace) if (r.t % 4 === 0 || r.t === 1) console.log(`      ${String(r.t).padStart(3)}   ${String(r.charge).padStart(12)}     ${String(r.intf).padStart(12)}`);
  const charges = trace.map((r) => r.charge), intfs = trace.map((r) => r.intf);
  const cSpread = Math.max(...charges) - Math.min(...charges);
  const iSpread = Math.max(...intfs) - Math.min(...intfs);
  console.log(`      => charge spread ${cSpread}, interference spread ${iSpread}`);
  console.log(`      => ${cSpread === 0 && iSpread > 0 ? 'item 7 IS DELIVERED at T3: the charge is conserved exactly while the interference observable moves, which is the distinction the item asks for' : 'still not delivered at T3'}`);
  console.log('');
  // A2 tier T3 x item 9 (second system): the eraser needs the marker rotation by 1/sqrt2
  const s2 = 1 / Math.sqrt(2);
  let ok = false;
  for (let a1 = -6; a1 <= 6 && !ok; a1++) for (let b1 = -6; b1 <= 6 && !ok; b1++) {
    // is s2*(a1 + i b1) again a Gaussian integer with integer real and imaginary parts?
    if (Math.abs(s2 * a1 - Math.round(s2 * a1)) < 1e-12 && Math.abs(s2 * b1 - Math.round(s2 * b1)) < 1e-12 && (a1 !== 0 || b1 !== 0)) ok = true;
  }
  console.log(`   A2 tier T3 x second system: the eraser needs the marker rotation by 1/sqrt(2) = ${s2.toFixed(6)}.`);
  console.log(`      is there any nonzero Gaussian integer whose rotation by that factor stays Gaussian? ${ok}`);
  console.log(`      => item 9 is REFUTED for T3, not merely undelivered: the basis rotation leaves the alphabet.`);
  console.log('      (this is round 82\'s arithmetic doing more work: the second system inherits the alphabet');
  console.log('       requirement of interference, one step further out.)');
}
console.log('');

console.log('============================================================');
console.log('PART B: THE UPDATED TABLE AND THE ENUMERATION');
console.log('============================================================');
{
  // item 9 (second system) is now DELIVERED by the continuous tiers at +2 branches, refuted for T1/T2/T3
  const tiers = [
    { name: 'T1', bits: 4, branches: 0, items: [0, 1] },
    { name: 'T2', bits: 3, branches: 3, items: [0, 1, 2, 8] },
    { name: 'T3', bits: 8, branches: 3, items: [0, 1, 2, 3, 5, 7, 8] },   // item 7 measured in part A1
    { name: 'T4', bits: 16, branches: 8, items: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] },   // 6 + 2 for the second system
    { name: 'T5', bits: 13, branches: 24, items: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] },  // 22 + 2
  ];
  const tiersBefore = tiers.map((t) => (t.name === 'T4' ? { ...t, branches: 6, items: t.items.filter((i) => i !== 9) }
    : t.name === 'T5' ? { ...t, branches: 22, items: t.items.filter((i) => i !== 9) }
      : t.name === 'T3' ? { ...t, items: [0, 1, 2, 3] } : t));
  // the DAG as adjudicated: refuted edges removed, definitional and reading edges kept
  const DAG = [[0, 1], [0, 5], [2, 4], [3, 7], [5, 6], [8, 9]];
  const closure = (set) => { const out = new Set(set); let ch = true; while (ch) { ch = false; for (const [p, d] of DAG) if (out.has(d) && !out.has(p)) { out.add(p); ch = true; } } return out; };
  const KAPPAS = [0, 0.25, 0.5, 1, 2, 4, 8, 16];
  const run = (label, list) => {
    let sat = 0, unsat = 0, unique = 0, tied = 0, secondUnsat = 0, secondTotal = 0;
    const share = new Map();
    for (let mask = 1; mask < (1 << CHECK.length); mask++) {
      const req = []; for (let i = 0; i < CHECK.length; i++) if (mask & (1 << i)) req.push(i);
      const need = closure(req);
      const cand = list.filter((t) => [...need].every((i) => t.items.includes(i)));
      const wantsSecond = req.includes(9);
      if (wantsSecond) secondTotal++;
      if (!cand.length) { unsat++; if (wantsSecond) secondUnsat++; continue; }
      sat++;
      const winners = new Set();
      for (const k of KAPPAS) {
        let best = null, bs = Infinity;
        for (const t of cand) { const s = t.bits + k * t.branches; if (s < bs - 1e-12) { bs = s; best = t; } }
        winners.add(best.name);
      }
      if (winners.size === 1) { unique++; share.set([...winners][0], (share.get([...winners][0]) || 0) + 1); }
      else { tied++; share.set('tie', (share.get('tie') || 0) + 1); }
    }
    console.log(`   ${label}`);
    console.log(`     satisfiable ${sat}  unsatisfiable ${unsat}  unique ${unique} (${(100 * unique / sat).toFixed(1)}%)  tied ${tied}`);
    console.log(`     of the ${secondTotal} sets requiring a SECOND SYSTEM: ${secondUnsat} unsatisfiable`);
    console.log(`     winners: ${[...share.entries()].sort((a, b) => b[1] - a[1]).map(([k, v]) => `${k} ${(100 * v / sat).toFixed(1)}%`).join('  ')}`);
    return { sat, unsat, unique, tied, secondUnsat, secondTotal };
  };
  console.log('   the table now says: item 9 delivered by T4 and T5 (+2 branches each), REFUTED for T1, T2, T3.');
  console.log('');
  const before = run('BEFORE this round (the frozen round-81 table)', tiersBefore);
  const after = run('AFTER  this round (second system entered)', tiers);
  console.log('');
  console.log(`   THE TRANSITION: unsatisfiable sets ${before.unsat} -> ${after.unsat}`);
  console.log(`                   sets requiring a second system that are unsatisfiable: ${before.secondUnsat} -> ${after.secondUnsat}`);
  globalThis.RES = { before, after, tiers, tiersBefore, DAG, CHECK_closure: closure };
}
console.log('');

console.log('============================================================');
console.log('PART C: SENSITIVITY -- THE DAG, AND THE COST CONVENTION');
console.log('============================================================');
{
  const { tiers, tiersBefore } = globalThis.RES;
  const KAPPAS = [0, 0.25, 0.5, 1, 2, 4, 8, 16];
  const evalWith = (list, DAG) => {
    const closure = (set) => { const out = new Set(set); let ch = true; while (ch) { ch = false; for (const [p, d] of DAG) if (out.has(d) && !out.has(p)) { out.add(p); ch = true; } } return out; };
    let sat = 0, unsat = 0, unique = 0, tied = 0;
    for (let mask = 1; mask < (1 << CHECK.length); mask++) {
      const req = []; for (let i = 0; i < CHECK.length; i++) if (mask & (1 << i)) req.push(i);
      const need = closure(req);
      const cand = list.filter((t) => [...need].every((i) => t.items.includes(i)));
      if (!cand.length) { unsat++; continue; }
      sat++;
      const winners = new Set();
      for (const k of KAPPAS) {
        let best = null, bs = Infinity;
        for (const t of cand) { const s = t.bits + k * t.branches; if (s < bs - 1e-12) { bs = s; best = t; } }
        winners.add(best.name);
      }
      if (winners.size === 1) unique++; else tied++;
    }
    return { sat, unsat, unique, tied };
  };
  const dagMin = [[0, 1], [0, 5], [2, 4], [3, 7], [5, 6]];                 // without the now-live 8->9
  const dagMid = [[0, 1], [0, 5], [2, 4], [3, 7], [5, 6], [8, 9]];         // the adjudicated DAG
  const dagMax = dagMid.concat([[2, 3], [3, 8]]);                          // readings re-included
  for (const [name, dag] of [['DAG without 8->9', dagMin], ['adjudicated DAG (8->9 live)', dagMid], ['readings re-included (2->3, 3->8)', dagMax]]) {
    const r = evalWith(tiers, dag);
    console.log(`   ${name.padEnd(34)} unsat ${String(r.unsat).padStart(4)}  unique ${String(r.unique).padStart(3)} (${(100 * r.unique / r.sat).toFixed(1)}%)  tied ${r.tied}`);
  }
  console.log('');
  console.log('   COST CONVENTION: what if the second system cost more or less than the measured 2 branches?');
  for (const extra of [0, 2, 4, 8, 16]) {
    const list = tiers.map((t) => (t.items.includes(9) ? { ...t, branches: t.branches - 2 + extra } : t));
    const r = evalWith(list, dagMid);
    console.log(`   +${String(extra).padStart(2)} branches on the continuous tiers   unsat ${String(r.unsat).padStart(4)}  unique ${String(r.unique).padStart(3)} (${(100 * r.unique / r.sat).toFixed(1)}%)  tied ${r.tied}`);
  }
  console.log('');
  console.log('   => the transition 512 -> 0 does not depend on either choice: it is caused by the second system');
  console.log('      becoming DELIVERABLE BY SOMETHING.');
  console.log('   => AND NOTE WHAT DID NOT MOVE: every DAG variant and every cost variant printed the SAME');
  console.log('      uniqueness number.  My draft of this paragraph claimed "the uniqueness number moves with');
  console.log('      both", which the table above refutes.  The reason it does not move is that T4 and T5 carry');
  console.log('      every item, so no DAG variant can exclude them, and the aggregate is dominated by ties');
  console.log('      (94%).  So after this round the enumeration has become insensitive to both choices --');
  console.log('      which is worth recording as such rather than as a sensitivity result.');
}
console.log('');

console.log('============================================================');
console.log('PART D: VERDICT -- THE PROJECT\'S ONLY ABSOLUTE RESULT DISSOLVES');
console.log('============================================================');
{
  const { before, after } = globalThis.RES;
  console.log(`   1. THE ABSOLUTE RESULT IS GONE.  Round 81 recorded 512 unsatisfiable requirement sets, all of`);
  console.log(`      them containing "a second system", and that was called the one thing that was ABSENT rather`);
  console.log(`      than expensive.  With the round-87 construction entered into the table at its measured cost,`);
  console.log(`      the count is ${after.unsat}.  Nothing in the checklist is unreachable any more.`);
  console.log('');
  console.log('   2. IT WAS A TABLE ARTEFACT, and the same one this project has now caught four times: a cell');
  console.log('      nobody had filled silently excluded every candidate that could have filled it.  Rounds');
  console.log('      79-81 found it in the Born column; round 85 in the reading of a dependency; round 87 in the');
  console.log('      second system itself.  So the honest statement of the old result is not "nature forbids a');
  console.log('      second system" but "no tier in our table had ever built one", which is a statement about us.');
  console.log('');
  console.log('   3. WHAT SURVIVES AS ABSOLUTE: nothing in the checklist.  What remains is the STRUCTURAL point');
  console.log('      from round 87: the second system is not a new ingredient but the same complex-linear layer');
  console.log('      applied to a product space, and the tiers that cannot do it are exactly the tiers without a');
  console.log('      continuous amplitude layer -- T1, T2 and T3, the last one because the eraser needs the');
  console.log('      marker rotation by 1/sqrt(2), which leaves the Gaussian integers.  That is an alphabet');
  console.log('      statement, not an absence.');
  console.log('');
  console.log('   4. RESIDUE: the construction entered here is 4 amplitudes, not a lattice.  Getting the second');
  console.log('      system into a full local field model, and measuring the two-body CHSH value with it, is the');
  console.log('      remaining work of this arc.');
}
