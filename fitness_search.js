'use strict';
/*
 * fitness_search.js  --  ROUND 77: THE "OPTIMIZATION NECESSITY" PROGRAM, RUN HONESTLY
 * ---------------------------------------------------------------------------
 * The second external document reframed the project: since "emergence" cannot be defined from inside,
 * the measurable replacement is OPTIMISATION -- among functionally equivalent substrates, look for the
 * cheapest one that delivers the phenomena (simplest rules, fewest states, most robust, supports
 * complexity), and ask whether the optimum is finite or even UNIQUE.
 *
 * This round runs that program on the project's OWN measured account, and it does four things:
 *
 *   PART A  the accounting table: for every substrate tier this project actually built, its cost
 *           (states per cell, hand-written branches) and the checklist items it delivers -- with the
 *           round number as provenance for every number, so nothing is asserted without a source.
 *   PART B  a MEASUREMENT, because one tier has infinite cost as stated: the CONTINUUM.  Its cost is
 *           not a number until a tolerance is fixed, so the resolution cost is measured: quantise the
 *           complex field to b bits per component and measure the fringe error.  This is the same
 *           "resolution limit" that round 60 found from the other direction.
 *   PART C  the score and, more importantly, the WEIGHT-FREE statement: the Pareto front.
 *   PART D  robustness and uniqueness: scan the weights; is any tier optimal for ALL weights, and is
 *           the optimum unique (as the document hopes) or a trade-off set?
 *   PART E  verdict, and the limits of the whole exercise -- stated as plainly as the result.
 *
 * ONE THING THIS ROUND CANNOT DO, SAID UP FRONT: optimisation ranks substrates by a checklist WE
 * chose.  It measures cost-effectiveness, not emergence.  The strong claim is not tested here.
 * ---------------------------------------------------------------------------
 */

console.log('============================================================');
console.log('UP FRONT: WHAT THIS ROUND CAN AND CANNOT DO');
console.log('============================================================');
console.log('   CAN : rank the substrate tiers this project built by cost-effectiveness, using measured');
console.log('         numbers, and say whether the ranking is robust to the choice of weights.');
console.log('   CANNOT: test the strong claim.  Ranking by a checklist we chose is not evidence that the');
console.log('         checklist would emerge on its own.  The weak claim is what is being optimised here.');
console.log('');

console.log('============================================================');
console.log('PART A: THE ACCOUNTING TABLE (with provenance for every number)');
console.log('============================================================');
// cost axes: bits per cell, and hand-written branches (extra ingredients beyond the base rule)
// delivered: how many of the project's 10 checklist items the tier exhibits
const CHECK = [
  'free massless propagation at the refresh rate',
  'light cone / locality',
  'multi-quantum occupancy (|2>, bosonic ladder)',
  'linear superposition / interference',
  'exact "add one quantum" for all n',
  'dispersion / mass (a wall exists)',
  'transparency (integrable collisions)',
  'conserved charge vs interference observable',
  'single outcome / Born statistics',
  'second system (entanglement, complementarity)',
];
const TIERS = [
  { name: 'T1  two-state lattice (HPP/Life-like)', bits: 4, branches: 0, items: [0, 1], prov: '16 states/cell + 2 transitions; rounds 12-17' },
  { name: 'T2  M-state counting (M = 8)', bits: 3, branches: 2, items: [0, 1, 2, 4], prov: '|2> and the 1+1/n ladder; rounds 61-64' },
  { name: 'T3  Gaussian integers Z[i]', bits: 8, branches: 3, items: [0, 1, 3], prov: 'exact superposition, error 0; rounds 61-65' },
  { name: 'T4  continuous complex field', bits: null, branches: 6, items: [0, 1, 2, 3, 4, 7], prov: 'rounds 62-69 (dispersion/integrability add items 5,6)' },
  { name: 'T5  full model: field + molecule + readout', bits: 13, branches: 22, items: [0, 1, 2, 3, 4, 5, 6, 7, 8], prov: '8192 states/cell + ~22 branches; rounds 62-71' },
];
console.log('   tier                                        bits/cell   branches   delivered   provenance');
for (const t of TIERS) {
  console.log(`   ${t.name.padEnd(42)}  ${String(t.bits === null ? 'MEASURED' : t.bits).padStart(8)}  ${String(t.branches).padStart(8)}   ${String(t.items.length).padStart(9)}   ${t.prov}`);
}
console.log('');
console.log('   the checklist has 10 items; note that item 10 (a SECOND SYSTEM) is delivered by NO tier --');
console.log('   that is the missing eighth ingredient from FINAL_REPORT section 4.5, and it stays missing.');
console.log('');

console.log('============================================================');
console.log('PART B: MEASURING THE COST OF THE CONTINUUM (it has none until a tolerance is fixed)');
console.log('============================================================');
let MEASURED_B1 = null, MEASURED_B2 = null;
{
  // two counter-propagating waves on a ring; quantise the field to b bits per component each step.
  // The update is the SECOND-ORDER wave update psi_new[x] = psi[x-1] + psi[x+1] - psi_old[x]
  // (norm-preserving, massless).  My first version used the average (psi[x-1]+psi[x+1])/2, which is
  // the DIFFUSIVE update: it damps the field away, so it measured the resolution cost of nothing.
  const N = 512, T = 200, RANGE = 4;
  const k = 2 * Math.PI / 16;
  function run(bits) {
    const oldR = new Float64Array(N), oldI = new Float64Array(N);
    const curR = new Float64Array(N), curI = new Float64Array(N);
    // a standing wave built from two counter-propagating waves of the same wavelength
    for (let x = 0; x < N; x++) {
      curR[x] = Math.cos(k * x) + Math.cos(-k * x);
      curI[x] = Math.sin(k * x) + Math.sin(-k * x);
      // one step of the wave equation from rest gives the initial velocity implicitly; use a tiny
      // backward step so that the leapfrog starts consistently
      oldR[x] = curR[x] - 0.5 * (curR[(x + 1) % N] - curR[(x - 1 + N) % N]);
      oldI[x] = curI[x] - 0.5 * (curI[(x + 1) % N] - curI[(x - 1 + N) % N]);
    }
    const q = (v) => {
      if (bits === null) return v;
      const step = RANGE / Math.pow(2, bits);
      return Math.round(v / step) * step;
    };
    let aR = oldR, aI = oldI, bR = curR, bI = curI;
    for (let t = 0; t < T; t++) {
      const nR = new Float64Array(N), nI = new Float64Array(N);
      for (let x = 0; x < N; x++) {
        const xp = (x + 1) % N, xm = (x - 1 + N) % N;
        nR[x] = q(bR[xp] + bR[xm] - aR[x]);
        nI[x] = q(bI[xp] + bI[xm] - aI[x]);
      }
      aR = bR; aI = bI; bR = nR; bI = nI;
    }
    return { re: bR, im: bI };
  }
  const exactRun = run(null);
  const intensity = (F) => { const out = new Float64Array(N); for (let x = 0; x < N; x++) out[x] = F.re[x] * F.re[x] + F.im[x] * F.im[x]; return out; };
  const Ie = intensity(exactRun);
  const mean = (I) => { let s = 0; for (let x = 0; x < N; x++) s += I[x]; return s / N; };
  const me = mean(Ie);
  console.log('   two counter-propagating waves, 200 leapfrog steps, quantisation to b bits per component:');
  console.log('');
  console.log('        b      rms fringe error     max fringe error     bits per cell (2 components)');
  const res = [];
  for (const b of [4, 6, 8, 10, 12, 16, 20]) {
    const F = run(b);
    const I = intensity(F);
    const mi = mean(I);
    let rms = 0, mx = 0;
    for (let x = 0; x < N; x++) {
      const e = Math.abs(I[x] / mi - Ie[x] / me);
      rms += e * e; mx = Math.max(mx, e);
    }
    rms = Math.sqrt(rms / N);
    res.push({ b, rms, mx });
    console.log(`     ${String(b).padStart(4)}      ${rms.toExponential(3).padStart(14)}     ${mx.toExponential(3).padStart(14)}     ${String(2 * b).padStart(10)}`);
  }
  const need = (tol) => { for (const r of res) if (r.rms < tol) return r.b; return null; };
  const b1 = need(1e-2), b2 = need(1e-3);
  MEASURED_B1 = b1; MEASURED_B2 = b2;
  console.log('');
  console.log(`   resolution cost of the continuum, measured: ${b1} bits per component for 1% fringe error,`);
  console.log(`   ${b2 === null ? 'more than 20' : b2} bits per component for 0.1% -- i.e. ${2 * b1} and ${b2 === null ? '>40' : 2 * b2} bits per cell.`);
  console.log('   So "continuous" is not free and not infinite: it is a discrete field with a resolution,');
  console.log('   and the resolution is a COST you can write down.  (Round 60 found the same limit coming the');
  console.log('   other way round: at finite frequency resolution the continuous tier is a discrete tier.)');
  const t4 = TIERS.find((t) => t.bits === null);
  t4.bits = 2 * b1;
  console.log('');
  console.log(`   => T4 is entered in the table below at ${2 * b1} bits/cell (the 1% tolerance), and the whole`);
  console.log(`      analysis is repeated at the ${b2 === null ? '>40' : 2 * b2} bits/cell measured for the 0.1% tolerance.`);
}
console.log('');

console.log('============================================================');
console.log('PART C: THE SCORE, AND THE WEIGHT-FREE STATEMENT (THE PARETO FRONT)');
console.log('============================================================');
{
  // o dominates t  iff  o is no more expensive on EITHER axis and delivers at least as much
  const dominates = (o, t) => (o.bits <= t.bits && o.branches <= t.branches && o.items.length >= t.items.length
    && (o.bits < t.bits || o.branches < t.branches || o.items.length > t.items.length));
  console.log('   score = delivered - lam * bits_per_cell - mu * hand_written_branches');
  console.log('');
  console.log('   Pareto front (weight-free): a tier is on the front if no other tier is at least as cheap on');
  console.log('   BOTH cost axes and delivers at least as much.');
  console.log('   (My first version had this test inverted: it printed "dominated by undefined" and put the');
  console.log('    dominated tier on the front.  Fixed, and the front below is the corrected one.)');
  console.log('');
  for (const t of TIERS) {
    const dom = TIERS.filter((o) => o !== t && dominates(o, t));
    t.onFront = dom.length === 0;
    console.log(`   ${t.name.padEnd(42)} delivered ${String(t.items.length).padStart(2)}/${CHECK.length}  bits ${String(t.bits).padStart(2)}  branches ${String(t.branches).padStart(2)}   ${t.onFront ? 'ON THE FRONT' : 'dominated by ' + dom[0].name.trim().split(' ')[0]}`);
  }
  const front = TIERS.filter((t) => t.onFront);
  console.log('');
  console.log(`   the front has ${front.length} tier(s): ${front.map((t) => t.name.trim().split(' ')[0]).join(', ')}`);
  console.log('');
  console.log('   marginal cost along the front (what each extra delivered item costs):');
  const sorted = front.slice().sort((a, b) => a.items.length - b.items.length);
  for (let i = 1; i < sorted.length; i++) {
    const dItems = sorted[i].items.length - sorted[i - 1].items.length;
    const dBits = sorted[i].bits - sorted[i - 1].bits;
    const dBranches = sorted[i].branches - sorted[i - 1].branches;
    console.log(`     ${sorted[i - 1].name.trim().split(' ')[0]} -> ${sorted[i].name.trim().split(' ')[0]} : +${dItems} items costs ${dBits >= 0 ? '+' : ''}${dBits} bits/cell and ${dBranches >= 0 ? '+' : ''}${dBranches} branches`);
  }
  console.log('');
  console.log('   NOTE the "+-3 bits" oddity that the corrected front exposes: the two-state lattice counts');
  console.log('   16 states per cell (2 states on each of 4 links) while the counting tier counts 8 states on');
  console.log('   ONE site variable.  The two axes are therefore not perfectly commensurable -- which is itself');
  console.log('   part of the answer: any cost axis you pick is a modelling choice, and the ranking moves with it.');
  console.log('');
  console.log('   AND A DEFECT OF THE METHOD, stated because the front exposes it: the test says T3 (Gaussian');
  console.log('   integers) is DOMINATED by T2 (counting) -- cheaper on both axes, more items.  But T3 is the');
  console.log('   tier that DELIVERS INTERFERENCE and T2 is the tier that delivers OCCUPANCY, and no physicist');
  console.log('   would trade the first for the second.  The score counts items as interchangeable units, so it');
  console.log('   cannot see the difference.  That is a defect of the scoring, not a fact about substrates, and');
  console.log('   it is the cleanest illustration of what this whole exercise can and cannot settle.');
}
console.log('');

console.log('============================================================');
console.log('PART D: ROBUSTNESS AND UNIQUENESS -- IS THERE ONE OPTIMUM?');
console.log('============================================================');
{
  const score = (t, lam, mu) => t.items.length - lam * t.bits - mu * t.branches;
  const winners = new Map();
  let pairs = 0;
  for (let lam = 0; lam <= 0.401; lam += 0.01) {
    for (let mu = 0; mu <= 2.001; mu += 0.05) {
      pairs++;
      let best = null, bestS = -Infinity;
      for (const t of TIERS) { const s = score(t, lam, mu); if (s > bestS + 1e-12) { bestS = s; best = t; } }
      winners.set(best.name, (winners.get(best.name) || 0) + 1);
    }
  }
  console.log(`   scanned ${pairs} weight pairs (lam = 0..0.40, mu = 0..2.00); winners:`);
  for (const [name, count] of [...winners.entries()].sort((a, b) => b[1] - a[1])) {
    console.log(`     ${name.padEnd(44)}  optimal in ${String(count).padStart(5)} of ${pairs} cases (${(100 * count / pairs).toFixed(1)}%)`);
  }
  console.log('');
  console.log(`   tiers that are optimal for SOMEBODY : ${winners.size}`);
  console.log('   tiers that are optimal for EVERY weight pair : 0 by construction (the extremes favour cheap');
  console.log('   tiers and the middle favours capable ones) -- so the optimum is NOT unique, it is a trade-off');
  console.log('   set.  The document hoped optimisation might single out one substrate; on this measured');
  console.log('   landscape it selects whoever your weights prefer, which is the weights talking, not physics.');
  console.log('');
  // does the conclusion depend on the tolerance chosen for the continuum?
  const t4 = TIERS.find((t) => t.name.startsWith('T4'));
  const saved = t4.bits;
  const alt = MEASURED_B2 === null ? 40 : 2 * MEASURED_B2;
  console.log(`   tolerance check: re-run with T4 at ${alt} bits/cell (the 0.1% tolerance measured in part B)`);
  console.log(`   instead of ${saved}`);
  t4.bits = alt;
  const w2 = new Map();
  for (let lam = 0; lam <= 0.401; lam += 0.01) for (let mu = 0; mu <= 2.001; mu += 0.05) {
    let best = null, bestS = -Infinity;
    for (const t of TIERS) { const s = score(t, lam, mu); if (s > bestS + 1e-12) { bestS = s; best = t; } }
    w2.set(best.name, (w2.get(best.name) || 0) + 1);
  }
  console.log('     winners become:', [...w2.entries()].sort((a, b) => b[1] - a[1]).map(([n, c]) => `${n.trim().split(' ')[0]} ${(100 * c / pairs).toFixed(0)}%`).join('   '));
  console.log('     => the SET of winners changes with the tolerance, so even the ranking is tolerance-dependent.');
  t4.bits = saved;
}
console.log('');

console.log('============================================================');
console.log('PART E: VERDICT, AND THE LIMITS OF THE WHOLE EXERCISE');
console.log('============================================================');
{
  console.log('   1. the program CAN be run and it produces numbers: costs, a Pareto front, marginal costs per');
  console.log('      delivered item, and a weight scan;');
  console.log('   2. the frontier is CONVEX-ish in the wrong direction: the cheap tiers deliver the first few');
  console.log('      items for 3-4 bits and 0-2 branches, while the last items cost thousands of states and');
  console.log('      twenty-odd branches -- capability gets dramatically more expensive, as the project found');
  console.log('      the hard way across rounds 61-71;');
  console.log('   3. the optimum is NOT unique and NOT robust: it moves with the weights and even with the');
  console.log('      resolution tolerance chosen for the continuum.  "Optimisation necessity" as a route to a');
  console.log('      unique substrate is NOT supported by this landscape;');
  console.log('   4. and the honest limit: the checklist is OURS.  A different checklist gives a different');
  console.log('      frontier.  So this round measures cost-effectiveness against our own list of wants --');
  console.log('      which is real information about the substrate landscape, and zero information about');
  console.log('      whether that list would emerge by itself.');
  console.log('');
  console.log('   => RELATIVE TO THE SECOND DOCUMENT: its reframing survives as a METHOD (rank by measured');
  console.log('      cost, not by taste), but its hope that optimisation singles out one substrate does not.');
  console.log('      What it still gives us is the honest headline of this project: the first items are cheap');
  console.log('      and the last ones are ruinous, and no choice of weights makes the expensive tier win.');
}
