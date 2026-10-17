'use strict';
/*
 * fitness_audit.js  --  ROUND 79: AUDIT THE TABLE AND THE DAG AGAINST THE MEASURED RECORD
 * ---------------------------------------------------------------------------
 * Rounds 77-78 produced a fitness table (tier x checklist item) and a dependency DAG, and the whole
 * "optimisation necessity" analysis rests on both.  But round 78 ended with a non-robustness that
 * mattered: requiring INTERFERENCE named T3 (Gaussian integers) under two dependency structures and
 * T4/T5 under a third -- the difference being the edge "superposition presupposes occupancy".
 *
 * An edge is either a READING or a MEASUREMENT.  This round decides, cell by cell and edge by edge,
 * using the project's own recorded numbers:
 *
 *   PART A  the audit table: for every (tier, item) cell, whether the record MEASURED it delivered,
 *           MEASURED it not delivered, or never measured it -- with the round and the number.
 *   PART B  the DAG adjudication: for each edge, look for a tier that delivers the dependent item
 *           but NOT the prerequisite.  Such a tier REFUTES the edge; without one the edge stays a
 *           reading.  The contested edge (2 -> 3) is included explicitly.
 *   PART C  re-run the requirement analysis on the AUDITED table with the refuted edge removed, and
 *           report whether the round-78 non-robustness is resolved.
 *   PART D  the honest residue: every cell and edge still unmeasured, listed.
 *
 * THE TWO FACTS THIS ROUND TURNS ON, both from the record and both decisive:
 *   round 63 (STATUS 64.1): the state is the AMPLITUDE and the count is derived, so the counting
 *       tier's identification of states with counts is REFUTED -- and sqrt(2) is not a Gaussian
 *       integer, so the Z[i] tier cannot hold |2> exactly;
 *   round 64 (STATUS 65.1): adding one quantum with the phase preserved needs a' = a*sqrt((n+1)/n);
 *       for a Gaussian integer this is exact only for n = 0 -- exhaustively checked n = 1..200, zero
 *       exact cases.  So the Z[i] tier delivers INTERFERENCE but NOT the exact ladder.
 * ---------------------------------------------------------------------------
 */

const CHECK = [
  'propagation at the refresh rate',
  'light cone / locality',
  'multi-quantum occupancy (|2>)',
  'linear superposition / interference',
  'exact add-one-quantum for all n',
  'dispersion / mass (a wall)',
  'transparency (integrable collisions)',
  'conserved charge vs interference obs',
  'single outcome / Born statistics',
  'second system (entanglement)',
];
// what the RECORD actually measured for each tier (not what seems reasonable)
const TIERS = [
  {
    name: 'T1 two-state lattice', bits: 4, branches: 0, measured: [0, 1],
    prov: { 0: 'speed = refresh rate, rounds 12-17', 1: 'finite cone, rounds 12-17', 2: 'REFUTED: bits collapse 2 -> 1 (round 61, STATUS 62.2A)', 3: 'REFUTED: no signs, no interference (round 61, STATUS 62.2 note)' },
  },
  {
    name: 'T2 counting (M=8)', bits: 3, branches: 2, measured: [0, 1, 2],
    prov: { 0: 'speed exactly 1.0000 (round 61, STATUS 62.2C)', 2: 'counts 1,2,3,5 carried intact (round 61, STATUS 62.2A)', 3: 'REFUTED: non-negative counts are a cone, so there are no signs and no interference (round 61, STATUS 62.2 note)', 4: 'CHANGED: the ladder 1+1/n was ARITHMETIC, not measured (round 61, STATUS 62.2E), and the count-as-state reading was REFUTED (round 63, STATUS 64.1)' },
  },
  {
    name: 'T3 Gaussian integers', bits: 8, branches: 3, measured: [0, 1, 3],
    prov: { 3: 'fringes 4/2/0/2, phase moves them one cell per 90 deg (round 62, STATUS 63)', 2: 'REFUTED: sqrt(2) is not a Gaussian integer, so |2> has no state (round 63, STATUS 64.1)', 4: 'REFUTED: exact add-one only for n = 0, exhaustive n = 1..200 (round 64, STATUS 65.1)' },
  },
  {
    name: 'T4 continuous complex', bits: 16, branches: 6, measured: [0, 1, 2, 3, 4, 7],
    prov: { 4: 'ladder measured: ratios 2.000/1.500/1.250/1.200 (round 63, STATUS 64.1)', 7: 'L2 conserved exactly while |sum psi|^2 moves (round 69)', 5: 'NOT delivered here: the massless lattice cannot have a wall (round 68)' },
  },
  {
    name: 'T5 full model', bits: 13, branches: 22, measured: [0, 1, 2, 3, 4, 5, 6, 7, 8],
    prov: { 5: 'self-trapped wall needs dispersion (round 68)', 6: 'transparent collision with integrable structure, 7.04% vs 56.04% (round 68)', 8: 'single outcome, Born 0.498/0.253/0.000/0.249 (round 65)', 9: 'ABSENT in every tier: the missing eighth ingredient (FINAL_REPORT 4.5)' },
  },
];

console.log('============================================================');
console.log('PART A: THE AUDIT TABLE (measured, not plausible)');
console.log('============================================================');
{
  console.log('   legend: M = measured delivered, X = measured NOT delivered, ? = never measured');
  console.log('');
  const head = '   tier                          bits br  ' + CHECK.map((_, i) => String(i).padStart(3)).join('');
  console.log(head);
  for (const t of TIERS) {
    const cells = CHECK.map((_, i) => {
      if (t.measured.includes(i)) return '  M';
      const p = t.prov[i] || '';
      if (p.startsWith('REFUTED') || p.startsWith('CHANGED') || p.startsWith('NOT delivered')) return '  X';
      return '  ?';
    });
    console.log(`   ${t.name.padEnd(28)} ${String(t.bits).padStart(4)} ${String(t.branches).padStart(2)}  ${cells.join('')}`);
  }
  console.log('');
  console.log('   item index: 0 propagation | 1 light cone | 2 occupancy | 3 interference | 4 exact ladder');
  console.log('               5 dispersion   | 6 transparency | 7 conserved charge | 8 Born | 9 second system');
  console.log('');
  console.log('   THE CELLS THAT CHANGE AGAINST ROUNDS 77-78:');
  console.log('     T2, item 4 (exact add-one-quantum): rounds 77-78 marked it DELIVERED.  The record says');
  console.log('       otherwise: round 61 computed the ladder arithmetically and labelled it "arithmetic, not');
  console.log('       measurement" (STATUS 62.2E), and round 63 then REFUTED the state-space reading behind it');
  console.log('       (STATUS 64.1).  So it becomes X.');
  console.log('     T3, item 2 (occupancy): rounds 77-78 had it absent already, and the record confirms by');
  console.log('       refutation rather than by silence: |2> needs amplitude sqrt(2), not a Gaussian integer.');
  console.log('');
  console.log('   the audited delivered counts therefore become: T1 2, T2 3 (was 4), T3 3, T4 6, T5 9.');
}
console.log('');

console.log('============================================================');
console.log('PART B: ADJUDICATING EVERY DEPENDENCY EDGE');
console.log('============================================================');
{
  const EDGES = [
    [0, 1], [0, 2], [0, 3], [0, 5], [2, 4], [5, 6], [3, 7], [3, 8], [8, 9],
    [2, 3],   // the contested edge from round 78, included here on purpose
  ];
  const CONTESTED = new Set(['2->3']);
  console.log('   for each edge: is there a tier that delivers the DEPENDENT item but NOT the prerequisite?');
  console.log('   such a tier refutes the edge; without one the edge stays a reading, not a measurement.');
  console.log('');
  console.log('   edge        verdict              evidence');
  const refuted = [];
  for (const [p, d] of EDGES) {
    const ce = TIERS.filter((t) => t.measured.includes(d) && !t.measured.includes(p));
    const key = `${p}->${d}`;
    let verdict, ev;
    if (ce.length) { verdict = 'REFUTED by the record'; ev = `${ce.map((t) => t.name.split(' ')[0]).join(', ')} delivers ${d} without ${p}`; refuted.push([p, d]); }
    else if (TIERS.some((t) => t.measured.includes(d) && t.measured.includes(p))) { verdict = 'consistent, unproven'; ev = 'every tier with the dependent item also has the prerequisite'; }
    else { verdict = 'vacuous / unresolved'; ev = 'no tier delivers the dependent item at all'; }
    console.log(`   ${key.padEnd(10)}  ${verdict.padEnd(21)}  ${ev}${CONTESTED.has(key) ? '   <- the edge round 78 was not robust to' : ''}`);
  }
  console.log('');
  console.log(`   refuted by measurement: ${refuted.length ? refuted.map(([p, d]) => `${p}->${d}`).join(', ') : 'none'}`);
  console.log('');
  console.log('   => the contested edge 2 -> 3 ("superposition presupposes occupancy") IS REFUTED, and the');
  console.log('      counterexample is exactly the tier that motivated round 78: T3 delivers INTERFERENCE and');
  console.log('      cannot hold OCCUPANCY, both by measurement.  So the variant that produced the round-78');
  console.log('      non-robustness is excluded by the record, and no other edge is refuted.');
  global.REFUTED = refuted;
}
console.log('');

console.log('============================================================');
console.log('PART C: RE-RUN THE REQUIREMENT ANALYSIS ON THE AUDITED TABLE');
console.log('============================================================');
{
  const baseEdges = [[0, 1], [0, 2], [0, 3], [0, 5], [2, 4], [5, 6], [3, 7], [3, 8], [8, 9]];
  const closure = (set, edges) => { const out = new Set(set); let ch = true; while (ch) { ch = false; for (const [p, d] of edges) if (out.has(d) && !out.has(p)) { out.add(p); ch = true; } } return out; };
  const KAPPAS = [0, 0.25, 0.5, 1, 2, 4, 8, 16];
  const run = (label, edges, tiers) => {
    let sat = 0, unsat = 0, unique = 0, tied = 0;
    let interf = null;
    for (let mask = 1; mask < (1 << CHECK.length); mask++) {
      const req = []; for (let i = 0; i < CHECK.length; i++) if (mask & (1 << i)) req.push(i);
      const need = closure(req, edges);
      const cand = tiers.filter((t) => [...need].every((i) => t.measured.includes(i)));
      if (!cand.length) { unsat++; continue; }
      sat++;
      const winners = new Set();
      for (const k of KAPPAS) {
        let best = null, bs = Infinity;
        for (const t of cand) { const s = t.bits + k * t.branches; if (s < bs - 1e-12) { bs = s; best = t; } }
        winners.add(best.name.split(' ')[0]);
      }
      if (winners.size === 1) unique++; else tied++;
      if (mask === (1 << 3)) interf = [...winners];
    }
    console.log(`   ${label}`);
    console.log(`     satisfiable ${sat}, unsatisfiable ${unsat}, unique-and-weight-robust ${unique} (${(100 * unique / sat).toFixed(1)}%), tied ${tied}`);
    console.log(`     requiring INTERFERENCE alone -> winner(s): ${interf ? interf.join(' / ') : '(unsatisfiable)'}`);
  };
  // the round-78 table, for a real comparison: T2 was credited with item 4 there
  const TIERS78 = TIERS.map((t) => (t.name.startsWith('T2') ? { ...t, measured: [...t.measured, 4] } : t));
  console.log('   same enumeration as round 78, but on the audited table and with the refuted edge removed:');
  console.log('');
  run('AUDITED table (T2 item 4 corrected to X)', baseEdges, TIERS);
  console.log('');
  run('round-78 table (T2 credited with item 4)', baseEdges, TIERS78);
  console.log('');
  console.log('   => the audited run is the one to trust.  Note WHAT does and does not move: the interference');
  console.log('      answer is T3 in both (that requirement never involved T2), while the tie/unique split moves,');
  console.log('      because crediting T2 with an unmeasured item changes which requirement sets it can bid on.');
}
console.log('');

console.log('============================================================');
console.log('PART D: THE HONEST RESIDUE -- WHAT IS STILL NOT MEASURED');
console.log('============================================================');
{
  const unmeasured = [];
  for (const t of TIERS) for (let i = 0; i < CHECK.length; i++) {
    if (!t.measured.includes(i) && !(t.prov[i] && (t.prov[i].startsWith('REFUTED') || t.prov[i].startsWith('CHANGED') || t.prov[i].startsWith('NOT delivered')))) {
      unmeasured.push(`${t.name.split(' ')[0]} x item ${i}`);
    }
  }
  console.log(`   cells that are neither measured-delivered nor measured-not-delivered: ${unmeasured.length}`);
  console.log('     ' + unmeasured.join(', '));
  console.log('');
  console.log('   edges that stay READINGS (no counterexample found, but no direct measurement either):');
  console.log('     0->1, 0->2, 0->3, 0->5, 2->4, 5->6, 3->7, 3->8  (all "consistent, unproven")');
  console.log('     and 8->9 is vacuous: no tier delivers a second system, so the edge cannot be tested.');
  console.log('');
  console.log('   WHAT THIS MEANS FOR THE PROGRAM: after this audit the fitness table is measured where the');
  console.log('   record allows and explicitly "?" where it does not, and exactly one dependency edge has been');
  console.log('   settled -- by refutation.  The optimisation results of rounds 77-78 survive in shape, but one');
  console.log('   of their inputs (T2 delivering the exact ladder) was a reading and is now corrected.');
  console.log('');
  console.log('   => the practical lesson, which generalises beyond this round: in a fitness analysis, every');
  console.log('      cell you have not measured is a degree of freedom, and the cheapest-looking tier is');
  console.log('      usually the one whose cells were filled in by plausibility rather than by measurement.');
}
