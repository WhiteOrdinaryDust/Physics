'use strict';
/*
 * entanglement_e4.js   --   ROUND 72: ESCAPE HATCH E4
 * ---------------------------------------------------------------------------
 * THE USER'S PROPOSAL, MADE TESTABLE:
 *
 *   "entanglement should also be a hidden-variable matter -- only the hidden variable lives between
 *    the OBSERVER and ONE END of the observed entangled body."
 *
 * That is a precise statement about Bell's derivation.  Bell needs FOUR assumptions:
 *      (i)   outcomes are determined by the beables lambda
 *      (ii)  LOCALITY: Alice's outcome does not depend on Bob's setting
 *      (iii) MEASUREMENT INDEPENDENCE: p(lambda | a,b) = p(lambda)
 *      (iv)  no conspiracy / "free" settings
 * The proposal denies (iii), and denies it in a specific way: lambda is correlated with the observer's
 * own setting at one end.  Our own rounds 12-17 recorded a version of this ("dropping statistical
 * independence relaxes the bound from 2 to 4, not to 2*sqrt(2)"), so this round re-tests that with an
 * explicit, no-signalling construction -- and the earlier verdict turns out to need REFINING.
 *
 * THE MODEL CLASS (smallest interesting one, and no-signalling BY CONSTRUCTION)
 *      lambda = (r, s):  r in {+1,-1} is a SHARED random sign, uniform;  s in {0,1} is the
 *      "settings memory" carried by the source.
 *      responses:  A = r * alpha(a,s),  B = r * beta(b,s),   alpha, beta in {+1,-1}
 *      measurement dependence: p(s | a,b) is ARBITRARY.
 * Because r flips BOTH outcomes, every marginal is exactly 1/2 whatever p(s|a,b) is:
 *      p(A=+|a,b) = p(B=+|a,b) = 1/2   =>   NO SIGNALLING, automatically.
 * So this class is local + deterministic + no-signalling, and the ONLY thing it gives up is (iii).
 *
 * PARTS
 *   A  CONTROL: restore measurement independence (p(s|a,b) = p(s)) => the class IS the local
 *      hidden-variable model => max S = 2.000000.  (Enumerated exactly: S is linear in p, so the
 *      maximum sits at a vertex, p in {0,1}.)
 *   B  Drop it: the range opens.  An explicit PR-box construction reaches S = 4.000000.
 *   C  THE QUESTION: can this class reproduce the QUANTUM tuple exactly, S = 2*sqrt(2) = 2.828427,
 *      while staying no-signalling?  (Enumerate all 256 response tables and check per setting pair
 *      whether the target correlation lies in the reachable interval.)
 *   D  THE COST, measured: the fitted model contains one free distribution PER SETTING PAIR, i.e.
 *      4 tuned numbers for 4 targets -- and it has NO rule for an unforeseen setting pair.  Quantify
 *      that: show that the quantum model predicts every setting from ONE parameter, and that the
 *      simplest extension of the fitted model to a fifth setting deviates from -cos.
 * ---------------------------------------------------------------------------
 */

const TGT = (ta, tb) => -Math.cos(((ta - tb) * Math.PI) / 180);
const PAIRS = [[0, 45], [0, -45], [90, 45], [90, -45]];        // a in {0,90}, b in {45,-45}
const target = PAIRS.map(([a, b]) => TGT(a, b));
const S_of = (E) => Math.abs(E[0] + E[1] + E[2] - E[3]);

// ---- the model class: alpha[a][s], beta[b][s] in {+1,-1}; p[s] per setting pair
function responses(t) {                                        // t = 0..255 -> alpha, beta tables
  const alpha = [[1, 1], [1, 1]], beta = [[1, 1], [1, 1]];
  for (let i = 0; i < 4; i++) alpha[i < 2 ? 0 : 1][i % 2] = ((t >> i) & 1) ? 1 : -1;
  for (let i = 0; i < 4; i++) beta[i < 2 ? 0 : 1][i % 2] = ((t >> (i + 4)) & 1) ? 1 : -1;
  return { alpha, beta };
}
function E_of(alpha, beta, p) {                                 // p[ip][s] for the 4 setting pairs
  return PAIRS.map((_, ip) => {
    const [ai, bi] = [Math.floor(ip / 2), ip % 2];              // a index, b index
    return p[ip][0] * alpha[ai][0] * beta[bi][0] + p[ip][1] * alpha[ai][1] * beta[bi][1];
  });
}

console.log('============================================================');
console.log('A: CONTROL -- WITH MEASUREMENT INDEPENDENCE THE BOUND IS 2');
console.log('============================================================');
{
  let best = -Infinity, arg = null;
  for (let t = 0; t < 256; t++) {
    const { alpha, beta } = responses(t);
    // measurement independence: one shared p(s) for every setting pair
    for (const ps of [[1, 0], [0, 1]]) {
      const p = PAIRS.map(() => ps);
      const S = S_of(E_of(alpha, beta, p));
      if (S > best) { best = S; arg = { t, ps }; }
    }
  }
  console.log(`   p(s|a,b) = p(s) for all pairs (that IS assumption iii), S is linear in p so the maximum`);
  console.log(`   sits at a vertex p in {0,1}:  max S = ${best.toFixed(6)}`);
  console.log(`   Bell's classical bound is 2, quantum is 2*sqrt(2) = ${(2 * Math.SQRT2).toFixed(6)}`);
  console.log('');
  console.log(`  => ${Math.abs(best - 2) < 1e-9 ? 'PASS' : '*** FAIL ***'}: with measurement independence this class cannot beat 2 -- the control behaves.`);
}
console.log('');

console.log('============================================================');
console.log('B: DROP IT -- AND AN EXPLICIT PR BOX FALLS OUT (S = 4)');
console.log('============================================================');
{
  // alpha(a,s) = +1 ; beta(0,s) = +1 ; beta(1,0) = +1, beta(1,1) = -1 ; p(0|a,b=1) = a ? 0 : 1
  const alpha = [[1, 1], [1, 1]], beta = [[1, 1], [1, -1]];
  const p = [
    [1, 0],   // (a=0,b=0): E = +1
    [1, 0],   // (a=0,b=1): E = +1     <- p depends on a here: the measurement dependence
    [1, 0],   // (a=1,b=0): E = +1
    [0, 1],   // (a=1,b=1): E = -1
  ];
  const E = E_of(alpha, beta, p);
  console.log(`   E(a,b) = ${E.map((v) => v.toFixed(3)).join(', ')}`);
  console.log(`   S = |E00 + E01 + E10 - E11| = ${S_of(E).toFixed(6)}   (the PR box: the no-signalling maximum)`);
  console.log('   every marginal is 1/2 for every (a,b) because the shared sign r flips both outcomes:');
  console.log('   => this violates Bell WITHOUT signalling.  What it gives up is (iii), nothing else.');
  console.log('');
  console.log(`  => ${Math.abs(S_of(E) - 4) < 1e-9 ? 'PASS' : '*** FAIL ***'}: measurement dependence alone reaches 4 -- so it is NOT pinned at 2*sqrt(2).`);
}
console.log('');

console.log('============================================================');
console.log('C: THE QUESTION -- CAN THE QUANTUM TUPLE 2*sqrt(2) BE REACHED EXACTLY, NO-SIGNALLING?');
console.log('============================================================');
{
  console.log(`   targets E(a,b) = -cos(a-b):  ${target.map((v) => v.toFixed(6)).join(', ')}`);
  console.log(`   target S = ${S_of(target).toFixed(6)}`);
  console.log('');
  let hits = [];
  for (let t = 0; t < 256; t++) {
    const { alpha, beta } = responses(t);
    // per pair, p is free in the simplex: reachable interval is [min(q0,q1), max(q0,q1)]
    const ok = PAIRS.every((_, ip) => {
      const [ai, bi] = [Math.floor(ip / 2), ip % 2];
      const q0 = alpha[ai][0] * beta[bi][0], q1 = alpha[ai][1] * beta[bi][1];
      return target[ip] >= Math.min(q0, q1) - 1e-12 && target[ip] <= Math.max(q0, q1) + 1e-12;
    });
    if (ok) {
      const p = PAIRS.map((_, ip) => {
        const [ai, bi] = [Math.floor(ip / 2), ip % 2];
        const q0 = alpha[ai][0] * beta[bi][0], q1 = alpha[ai][1] * beta[bi][1];
        if (Math.abs(q0 - q1) < 1e-12) return [1, 0];
        const x = (target[ip] - q1) / (q0 - q1);
        return [x, 1 - x];
      });
      hits.push({ t, alpha, beta, p, E: E_of(alpha, beta, p) });
    }
  }
  console.log(`   response tables that can hit ALL FOUR targets exactly: ${hits.length} of 256`);
  if (hits.length) {
    const h = hits[0];
    console.log(`   example table ${h.t}:  alpha = ${JSON.stringify(h.alpha)}  beta = ${JSON.stringify(h.beta)}`);
    console.log(`   the required MEASUREMENT-DEPENDENT distributions p(s|a,b):`);
    PAIRS.forEach(([a, b], ip) => console.log(`     (a=${String(a).padStart(2)}°, b=${String(b).padStart(3)}°) : p(s=0) = ${h.p[ip][0].toFixed(6)}`));
    console.log(`   achieved E = ${h.E.map((v) => v.toFixed(6)).join(', ')}`);
    console.log(`   achieved S = ${S_of(h.E).toFixed(6)}`);
    console.log('');
    console.log(`  => ${Math.abs(S_of(h.E) - 2 * Math.SQRT2) < 1e-9 ? 'PASS' : '*** FAIL ***'}: THE QUANTUM VALUE IS REACHED EXACTLY by a model that is local, deterministic`);
    console.log('     and NO-SIGNALLING -- it gives up only measurement independence.');
    console.log('     ⇒ THE USER\'S PROPOSAL WORKS.  Our rounds 12-17 sentence "superdeterminism is');
    console.log('       necessary but not SUFFICIENT" was too strong: it IS sufficient to reproduce the value,');
    console.log('       in this explicitly no-signalling class.  What it is not sufficient for is part D.');
  } else {
    console.log('  => not reachable in this class; the earlier verdict stands unchanged.');
  }
}
console.log('');

console.log('============================================================');
console.log('D: THE COST, MEASURED -- IT RETRODICTS, IT DOES NOT PREDICT');
console.log('============================================================');
{
  // the fitted model carries one free distribution per setting pair: 4 tuned numbers for 4 targets
  console.log('   the model fitted in part C contains p(s|a,b) for the FOUR setting pairs it was fitted to:');
  console.log('     tuned parameters : 4 (one distribution per pair)  +  8 response signs');
  console.log('     quantities fitted: 4 correlations  =>  zero predictive slack');
  console.log('');
  console.log('   a model that PREDICTS would fix its parameters ONCE and then compute any setting pair.');
  console.log('   The quantum singlet does exactly that: one parameter (the state |psi>) gives -cos(a-b) for');
  console.log('   EVERY pair.  Its parameter count is independent of how many settings you test:');
  const cols = [4, 8, 16, 64, 256];
  console.log('      settings tested : ' + cols.map((n) => String(n).padStart(5)).join(''));
  console.log('      quantum model   : ' + cols.map(() => String(1).padStart(5)).join('') + '   <- one parameter, independent of the number of settings');
  console.log('      fitted E4 model : ' + cols.map((n) => String(n).padStart(5)).join('') + '   <- one tuning per setting pair');
  console.log('');
  console.log('   and the extension is not free: take the simplest rule for an UNFORESEEN pair, "reuse the');
  console.log('   nearest tuned pair", and ask what the fitted model says at (a=22.5°, b=0°):');
  {
    const TGT2 = (ta, tb) => -Math.cos(((ta - tb) * Math.PI) / 180);
    const quantum = TGT2(22.5, 0);
    // nearest tuned pair to (22.5, 0): among the four, (0,45) or (0,-45) -- both at 45 deg away
    const nearest = [TGT2(0, 45), TGT2(0, -45)];              // the model has no entry for 22.5°
    console.log(`      quantum prediction at (22.5°, 0°) : ${quantum.toFixed(6)}`);
    console.log(`      nearest tuned entries             : ${nearest.map((v) => v.toFixed(6)).join(' , ')}`);
    console.log(`      deviation of that extension       : ${Math.abs(quantum - nearest[0]).toFixed(6)}  (and the choice of which neighbour is arbitrary)`);
  }
  console.log('');
  console.log('   ⇒ the fitted model has NO answer at settings it was not tuned for; any answer requires');
  console.log('     MORE tuning.  That is the whole cost of E4, and it is not a technicality:');
  console.log('     a derived mechanism would give the correlation p(lambda|a,b) from the substrate\'s rules,');
  console.log('     with a parameter count that does not grow with the number of settings tested.');
}
console.log('');

console.log('============================================================');
console.log('E4 VERDICT');
console.log('============================================================');
{
  console.log('   WHAT E4 BUYS (measured):');
  console.log('     * with measurement independence restored, this class stops at S = 2 exactly   (part A)');
  console.log('     * dropping it -- while keeping locality, determinism and no-signalling -- reaches');
  console.log('       S = 4 (an explicit PR box)                                                  (part B)');
  console.log('     * and it can sit EXACTLY at the quantum value S = 2*sqrt(2) = 2.828427         (part C)');
  console.log('   ⇒ the proposal is correct about WHERE the escape is, and it is enough to reproduce the value.');
  console.log('');
  console.log('   WHAT IT COSTS (measured):');
  console.log('     * the source must carry a distribution p(s|a,b) that is correlated with the observers\'');
  console.log('       OWN settings -- that is the "hidden variable between the observer and the entangled body"');
  console.log('     * fitted, that is one free distribution per setting pair, so the parameter count grows with');
  console.log('       the number of settings tested, while the model to be reproduced uses ONE parameter  (part D)');
  console.log('     * the model has no prediction for an unforeseen setting pair                       (part D)');
  console.log('');
  console.log('   SO THE HONEST STATEMENT IS:');
  console.log('     "Entanglement is a hidden-variable matter, and the hidden variable lives between the');
  console.log('      observer and one end" -- YES, and it reproduces 2*sqrt(2) exactly without signalling.');
  console.log('      The open task is not to ALLOW that correlation but to DERIVE it from the substrate,');
  console.log('      with a parameter count that does not grow with the number of settings tested.');
  console.log('     (This also corrects our own FINAL_REPORT 14.4: a Bell test does NOT exclude hidden-variable');
  console.log('      substrates -- it excludes them only if measurement independence is also assumed.)');
}
