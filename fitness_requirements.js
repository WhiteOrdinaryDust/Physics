'use strict';
/*
 * fitness_requirements.js  --  ROUND 78: REQUIREMENTS FIRST, COST SECOND
 * ---------------------------------------------------------------------------
 * Round 77 showed the defect plainly: scoring "number of checklist items delivered" treated the items
 * as interchangeable units, so the Gaussian-integer tier (which delivers INTERFERENCE) came out
 * dominated by the counting tier (which delivers OCCUPANCY).  No physicist would make that trade.
 *
 * Engineering does not rank by counting interchangeable features; it states REQUIREMENTS and then
 * minimises cost.  This round does exactly that, and it asks three questions:
 *
 *   PART A  the items are not independent: build the PREREQUISITE structure (a DAG) with provenance
 *           for every edge, so "delivered" can be closed under requirements.
 *   PART B  requirement-based selection: for every subset of requirements (closed under prerequisites)
 *           find the cheapest tier.  Report how many requirement sets are satisfiable at all, which
 *           are NOT, and -- the interesting number -- how often the answer is UNIQUE and robust to the
 *           cost weights.  This is the honest form of "optimisation necessity": uniqueness relative to
 *           a stated requirement, not uniqueness in the absolute.
 *   PART C  data-driven item weights: weight each item by how many other items it ENABLES (its
 *           out-degree in the DAG), re-rank, and check whether the round-77 anomaly disappears.
 *   PART D  what survives every choice made here: the robust statements, including the absolute one.
 *   PART E  verdict.
 *
 * STATED UP FRONT, because it is the point: the DAG edges are this project's own claims about what
 * depends on what, each with a round number.  A different dependency structure would move the numbers;
 * PART D therefore reports only what holds across the variants tested.
 * ---------------------------------------------------------------------------
 */

const CHECK = [
  'propagation at the refresh rate',      // 0
  'light cone / locality',               // 1
  'multi-quantum occupancy (|2>)',       // 2
  'linear superposition / interference', // 3
  'exact add-one-quantum for all n',     // 4
  'dispersion / mass (a wall)',          // 5
  'transparency (integrable collisions)',// 6
  'conserved charge vs interference obs',// 7
  'single outcome / Born statistics',    // 8
  'second system (entanglement)',        // 9
];
// prerequisite edges: [prerequisite, dependent, provenance]
const EDGES = [
  [0, 1, 'a light cone presupposes propagation (round 62)'],
  [0, 2, 'occupancy presupposes a carrier (rounds 61-64)'],
  [0, 3, 'superposition presupposes a field to superpose (round 62)'],
  [0, 5, 'mass/dispersion presupposes propagation (round 68)'],
  [2, 4, 'the exact ladder extends occupancy (round 64)'],
  [5, 6, 'transparency presupposes colliding objects (round 68)'],
  [3, 7, 'a conserved charge vs an interference observable needs both (round 69)'],
  [3, 8, 'Born statistics needs amplitudes (round 63)'],
  [8, 9, 'a second system is what complementarity/entanglement requires (rounds 65-66)'],
];
const TIERS = [
  { name: 'T1 two-state lattice', bits: 4, branches: 0, items: [0, 1] },
  { name: 'T2 counting (M=8)', bits: 3, branches: 2, items: [0, 1, 2, 4] },
  { name: 'T3 Gaussian integers', bits: 8, branches: 3, items: [0, 1, 3] },
  { name: 'T4 continuous complex', bits: 16, branches: 6, items: [0, 1, 2, 3, 4, 7] },
  { name: 'T5 full model', bits: 13, branches: 22, items: [0, 1, 2, 3, 4, 5, 6, 7, 8] },
];

console.log('============================================================');
console.log('PART A: THE ITEMS ARE NOT INDEPENDENT -- THE PREREQUISITE DAG');
console.log('============================================================');
{
  console.log('   edge                                    provenance');
  for (const [p, d, why] of EDGES) console.log(`   ${CHECK[p].padEnd(38)} -> ${CHECK[d].padEnd(36)} ${why}`);
  const closure = (set) => {
    const out = new Set(set);
    let changed = true;
    while (changed) {
      changed = false;
      for (const [p, d] of EDGES) if (out.has(d) && !out.has(p)) { out.add(p); changed = true; }
    }
    return out;
  };
  const depth = new Array(CHECK.length).fill(0);
  for (let it = 0; it < 20; it++) {
    let ch = false;
    for (const [p, d] of EDGES) if (depth[d] < depth[p] + 1) { depth[d] = depth[p] + 1; ch = true; }
    if (!ch) break;
  }
  const enables = new Array(CHECK.length).fill(0);
  for (const [p] of EDGES) enables[p]++;
  console.log('');
  console.log('   item                                       depth   enables (out-degree)');
  CHECK.forEach((c, i) => console.log(`   ${c.padEnd(42)} ${String(depth[i]).padStart(4)}   ${String(enables[i]).padStart(6)}`));
  console.log('');
  console.log('   the out-degree is a data-driven measure of how FOUNDATIONAL an item is: an item that many');
  console.log('   others depend on is not interchangeable with one that nothing depends on.');
  // what each tier actually delivers after closing under prerequisites
  console.log('');
  console.log('   tier                          delivered (raw)   delivered (closed under prerequisites)');
  for (const t of TIERS) {
    const cl = closure(t.items);
    t.closure = [...cl].sort((a, b) => a - b);
    console.log(`   ${t.name.padEnd(28)}  ${String(t.items.length).padStart(6)}            ${String(cl.size).padStart(6)}   [${t.closure.join(',')}]`);
  }
  global.CLOSURE = closure; global.DEPTH = depth; global.ENABLES = enables;
}
console.log('');

console.log('============================================================');
console.log('PART B: REQUIREMENTS FIRST, COST SECOND (the honest form of "optimisation necessity")');
console.log('============================================================');
{
  const closure = global.CLOSURE;
  const cost = (t, kappa) => t.bits + kappa * t.branches;
  const KAPPAS = [0, 0.25, 0.5, 1, 2, 4, 8, 16];
  let sat = 0, unsat = 0, unique = 0, tied = 0;
  const unsatWhy = new Map();
  const examples = [];
  for (let mask = 0; mask < (1 << CHECK.length); mask++) {
    if (mask === 0) continue;
    const req = [];
    for (let i = 0; i < CHECK.length; i++) if (mask & (1 << i)) req.push(i);
    const need = closure(req);
    // a tier qualifies iff it DELIVERS EVERY required item (my first version had this inverted, asking
    // whether the tier needs nothing beyond the requirement -- which produced a one-off count error
    // and 255 unsatisfiable sets with "no single item explaining it")
    const cand = TIERS.filter((t) => [...need].every((i) => t.items.includes(i)));
    if (cand.length === 0) {
      unsat++;
      const missing = [...need].filter((i) => TIERS.every((t) => !t.items.includes(i)));
      const key = missing.map((i) => CHECK[i]).join(' + ') || '(no single item explains it)';
      unsatWhy.set(key, (unsatWhy.get(key) || 0) + 1);
      continue;
    }
    sat++;
    // weight-robust uniqueness: the same tier wins for EVERY kappa
    const winners = new Set();
    for (const k of KAPPAS) {
      let best = null, bs = Infinity;
      for (const t of cand) { const s = cost(t, k); if (s < bs - 1e-12) { bs = s; best = t; } }
      winners.add(best.name);
    }
    if (winners.size === 1) unique++; else tied++;
    if (examples.length < 6 && req.length <= 2) examples.push({ req: req.map((i) => CHECK[i]), need: [...need].map((i) => CHECK[i]), winners: [...winners] });
  }
  console.log(`   all ${(1 << CHECK.length) - 1} non-empty requirement sets tested (each closed under prerequisites)`);
  console.log('');
  console.log(`     satisfiable                      : ${sat}`);
  console.log(`     NOT satisfiable by ANY tier      : ${unsat}  (${(100 * unsat / (sat + unsat)).toFixed(1)}%)`);
  console.log(`     of the satisfiable ones, UNIQUE answer robust to the cost weights : ${unique} (${(100 * unique / sat).toFixed(1)}%)`);
  console.log(`                              tied between tiers                      : ${tied} (${(100 * tied / sat).toFixed(1)}%)`);
  console.log('');
  console.log('   why the unsatisfiable ones fail (counted by the missing item):');
  for (const [k, v] of [...unsatWhy.entries()].sort((a, b) => b[1] - a[1])) console.log(`     ${String(v).padStart(4)} sets require : ${k}`);
  console.log('');
  console.log('   worked examples (requirement -> cheapest tier, unique or tied):');
  for (const e of examples) {
    console.log(`     require [${e.req.join(' + ')}]`);
    console.log(`       closure needs [${e.need.join(' + ')}]`);
    console.log(`       winner(s) across all cost weights: ${e.winners.join(' / ')}`);
  }
  console.log('');
  console.log(`   => ${unique} of ${sat} satisfiable requirement sets (${(100 * unique / sat).toFixed(1)}%) have a UNIQUE cheapest tier that is`);
  console.log(`      robust to the cost weights; the other ${tied} are ties, and every tie is a WEAK requirement set (e.g.`);
  console.log('      "propagation" alone is satisfied by all five tiers, so the cheapest flips between T1 and T2 as');
  console.log('      the branch cost is weighted).  Round 77\'s non-uniqueness came from NOT stating requirements.');
}
console.log('');

console.log('============================================================');
console.log('PART C: DATA-DRIVEN ITEM WEIGHTS -- DOES THE ROUND-77 ANOMALY GO AWAY?');
console.log('============================================================');
{
  const enables = global.ENABLES;
  const w = CHECK.map((_, i) => 1 + enables[i]);
  const scoreW = (t) => t.items.reduce((s, i) => s + w[i], 0);
  console.log('   weight of an item = 1 + (how many items depend on it)  [out-degree + 1]');
  console.log('');
  console.log(`   ${'item'.padEnd(42)} weight`);
  CHECK.forEach((c, i) => console.log(`   ${c.padEnd(42)} ${String(w[i]).padStart(5)}`));
  console.log('');
  console.log('   tier                          raw count   weighted score   weighted score per bit of cost');
  const rows = TIERS.map((t) => ({ t, raw: t.items.length, ws: scoreW(t), eff: scoreW(t) / t.bits }));
  for (const r of rows) console.log(`   ${r.t.name.padEnd(28)}  ${String(r.raw).padStart(9)}   ${String(r.ws).padStart(14)}   ${r.eff.toFixed(3).padStart(10)}`);
  console.log('');
  const rankW = rows.slice().sort((a, b) => b.eff - a.eff);
  console.log('   ranking by weighted score per bit : ' + rankW.map((r) => r.t.name.split(' ')[0]).join(' > '));
  const t2 = rows.find((r) => r.t.name.startsWith('T2')), t3 = rows.find((r) => r.t.name.startsWith('T3'));
  console.log('');
  console.log(`   the round-77 anomaly was T3 (interference) dominated by T2 (occupancy).  With foundationality`);
  console.log(`   weights: T3 weighted score ${t3.ws} at ${t3.t.bits} bits (efficiency ${t3.eff.toFixed(3)}) versus T2 ${t2.ws} at ${t2.t.bits} bits (${t2.eff.toFixed(3)}).`);
  console.log(`   => ${t3.eff > t2.eff ? 'the anomaly is GONE: interference now counts for more than occupancy, as it should.' : 'the anomaly REMAINS: even weighted by what they enable, occupancy outscores interference here.'}`);
}
console.log('');

console.log('============================================================');
console.log('PART D: WHAT SURVIVES EVERY CHOICE MADE IN THIS ROUND');
console.log('============================================================');
{
  // variant DAGs: does the headline change if an edge is added or removed?
  const variants = {
    'base': EDGES,
    'no 3->8 (Born independent of superposition)': EDGES.filter(([p, d]) => !(p === 3 && d === 8)),
    'extra 2->3 (occupancy needed for superposition)': EDGES.concat([[2, 3, 'variant']]),
  };
  console.log('   checking the headline claims against dependency-structure variants:');
  console.log('');
  console.log('     variant                                          sets needing item 10   sets with no tier');
  for (const [label, edges] of Object.entries(variants)) {
    const closure = (set) => { const out = new Set(set); let ch = true; while (ch) { ch = false; for (const [p, d] of edges) if (out.has(d) && !out.has(p)) { out.add(p); ch = true; } } return out; };
    let need10 = 0, none = 0, tested = 0;
    for (let mask = 1; mask < (1 << CHECK.length); mask++) {
      const req = []; for (let i = 0; i < CHECK.length; i++) if (mask & (1 << i)) req.push(i);
      const need = closure(req); tested++;
      if (need.has(9)) need10++;
      if (TIERS.every((t) => ![...need].every((i) => t.items.includes(i)))) none++;
    }
    console.log(`     ${label.padEnd(48)} ${String(need10).padStart(10)} / ${tested}      ${String(none).padStart(6)}`);
  }
  console.log('');
  console.log('   STRONGER ROBUSTNESS TEST: do the WINNERS for the worked examples move when the dependency');
  console.log('   structure changes?  (This is the part of the analysis the DAG could actually break.)');
  const exampleReqs = [[0], [1], [2], [3], [5], [3, 8]];
  console.log('');
  console.log('     requirement(s)                      base                                      no 3->8                                    extra 2->3');
  for (const req of exampleReqs) {
    const cells = [];
    for (const edges of Object.values(variants)) {
      const closure = (set) => { const out = new Set(set); let ch = true; while (ch) { ch = false; for (const [p, d] of edges) if (out.has(d) && !out.has(p)) { out.add(p); ch = true; } } return out; };
      const need = closure(req);
      const cand = TIERS.filter((t) => [...need].every((i) => t.items.includes(i)));
      const winners = new Set();
      for (const k of [0, 0.25, 0.5, 1, 2, 4, 8, 16]) {
        let best = null, bs = Infinity;
        for (const t of cand) { const s = t.bits + k * t.branches; if (s < bs - 1e-12) { bs = s; best = t; } }
        if (best) winners.add(best.name.trim().split(' ')[0]);
      }
      cells.push([...winners].join('/') || 'none');
    }
    console.log(`     ${req.map((i) => CHECK[i].slice(0, 16)).join(' + ').padEnd(38)} ${cells.map((c) => c.padEnd(42)).join('')}`);
  }
  console.log('');
  console.log('   ROBUST ACROSS ALL VARIANTS:');
  console.log('     * every requirement set that includes the SECOND SYSTEM is unsatisfiable -- no tier in this');
  console.log('       project delivers it at any cost.  This is the one ABSOLUTE statement in the whole');
  console.log('       optimisation exercise, and it is the same fact as FINAL_REPORT section 4.5\'s missing');
  console.log('       eighth ingredient, now expressed as a requirement with no solution;');
  console.log('     * the cheap tiers deliver the first items cheaply and the last items are ruinous;');
  console.log('     * uniqueness appears only relative to a stated requirement set.');
  console.log('');
  console.log('   NOT ROBUST (moves with the dependency structure and the weights, so it is a choice, not a fact):');
  console.log('     * exactly which tier wins for a given requirement set when there is a tie;');
  console.log('     * and here is the sharpest instance, in the table above: requiring INTERFERENCE names T3 as');
  console.log('       the cheapest tier under the base structure and under "Born independent of superposition",');
  console.log('       but names T4/T5 as soon as the structure also says "superposition presupposes occupancy".');
  console.log('       All three variants are defensible readings of this project\'s own results, and they');
  console.log('       disagree about who provides the one item that motivated this whole round.  So even the');
  console.log('       most load-bearing requirement is only as stable as the dependency claim under it.');
}
console.log('');

console.log('============================================================');
console.log('PART E: VERDICT');
console.log('============================================================');
{
  console.log('   1. TWO ROUTES WERE TRIED and they behaved differently -- this is the main result of the round:');
  console.log('        * WEIGHTING the items by how foundational they are (part C) does NOT repair round 77\'s');
  console.log('          anomaly: T3 (interference) still loses to T2 (occupancy) on weighted score per bit,');
  console.log('          because the anomaly is about COST (8 bits against 3), not about counting.  No additive');
  console.log('          weighting of a checklist can see that interference is worth more than occupancy;');
  console.log('        * REQUIRING first and costing second (part B) DOES repair it: 464 of 511 satisfiable');
  console.log('          requirement sets (90.8%) get a unique cheapest tier, and that tier is the same for');
  console.log('          every cost weighting tried.  As soon as interference is required, T2 is excluded.');
  console.log('   2. what remains a choice is the REQUIREMENT SET (and the dependency structure).  So');
  console.log('      "optimisation necessity" is real but RELATIVE: it says "given what you need, this is the');
  console.log('      cheapest way", not "this is the substrate the world must have";');
  console.log('   3. and the one absolute result: the requirement "a second system / entanglement" has NO');
  console.log('      solution anywhere in this landscape -- exactly the 512 unsatisfiable requirement sets, and');
  console.log('      the same fact rounds 65-66 found from the other direction.  Every other capability is');
  console.log('      merely expensive; that one is absent.');
  console.log('');
  console.log('   => relative to round 77, the program gains a well-posed question and loses a hope:');
  console.log('      "rank by taste" is replaced by "state requirements, then minimise cost" (a good method),');
  console.log('      while "optimisation will single out the substrate" survives only after you have already');
  console.log('      said what you want -- which is the part physics was supposed to tell you.');
}
