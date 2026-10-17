'use strict';
/*
 * e4_anatomy.js  --  WHAT THE ROUND-72 MODEL ACTUALLY IS (an anatomy, for the record)
 * ---------------------------------------------------------------------------
 * Round 72 is the round that landed on 2*sqrt(2) EXACTLY (S = 2.828427).  This script dissects it:
 * how many solutions there are, what the simplest one looks like, how no-signalling comes for free,
 * and where the "2*sqrt(2) rather than 4" comes from.  Nothing new is claimed here; it is the same
 * model, taken apart.
 *
 * THE DESIGN, IN TWO MOVES
 *   move 1 (the free lunch):   lambda = (r, s) with r in {+1,-1} UNIFORM and SHARED, and the responses
 *                              are A = r*alpha(a,s), B = r*beta(b,s).  Because the SAME r multiplies
 *                              both outcomes, every marginal is exactly 1/2 whatever else happens:
 *                                   p(A=+|a,b) = sum_s p(s|a,b) * (1/2) = 1/2
 *                              => NO-SIGNALLING BY CONSTRUCTION, not by tuning.
 *   move 2 (the price):        let p(s|a,b) depend on BOTH settings (measurement dependence).  Then
 *                                   E(a,b) = p*q0 + (1-p)*q1,   q0,q1 = alpha*beta in {+1,-1}
 *                              and any target E in [-1,1] is reachable by  p = (1-E)/2  (or (1+E)/2).
 *                              The correlation moves OUT of the response functions and INTO the
 *                              source's bookkeeping about the settings.
 * ---------------------------------------------------------------------------
 */

const PAIRS = [[0, 45], [0, -45], [90, 45], [90, -45]];
const TGT = (a, b) => -Math.cos(((a - b) * Math.PI) / 180);
const target = PAIRS.map(([a, b]) => TGT(a, b));
const S_of = (E) => Math.abs(E[0] + E[1] + E[2] - E[3]);
const responses = (t) => {
  const alpha = [[1, 1], [1, 1]], beta = [[1, 1], [1, 1]];
  for (let i = 0; i < 4; i++) alpha[i < 2 ? 0 : 1][i % 2] = ((t >> i) & 1) ? 1 : -1;
  for (let i = 0; i < 4; i++) beta[i < 2 ? 0 : 1][i % 2] = ((t >> (i + 4)) & 1) ? 1 : -1;
  return { alpha, beta };
};

console.log('============================================================');
console.log('1. WHAT IS BEING ASKED OF THE MODEL');
console.log('============================================================');
console.log('   four requirements at once:');
console.log('     (1) LOCAL      : Alice\'s outcome may not depend on Bob\'s setting');
console.log('     (2) DETERMINISTIC given lambda');
console.log('     (3) NO-SIGNALLING : each marginal independent of the other setting');
console.log('     (4) VIOLATE BELL, and land on the QUANTUM value 2*sqrt(2) = ' + (2 * Math.SQRT2).toFixed(6));
console.log('   a measurement-INDEPENDENT local model cannot do (4): it stops at 2.  So one of Bell\'s other');
console.log('   assumptions has to go -- the proposal is measurement independence, p(lambda|a,b) = p(lambda).');
console.log('');

console.log('============================================================');
console.log('2. THE FUNCTIONAL FORM, AND WHY ANY TARGET IS REACHABLE');
console.log('============================================================');
console.log('   E(a,b) = sum_s p(s|a,b) * alpha(a,s) * beta(b,s)');
console.log('   with s in {0,1} this is  E = p*q0 + (1-p)*q1,  q_s = alpha(a,s)*beta(b,s) in {+1,-1}');
console.log('   => the reachable interval for that setting pair is [min(q0,q1), max(q0,q1)]');
console.log('   => and p is then forced: p = (E - q1)/(q0 - q1)');
console.log('');

console.log('============================================================');
console.log('3. HOW MANY SOLUTIONS? AND WHAT DOES THE SIMPLEST ONE LOOK LIKE?');
console.log('============================================================');
{
  const hits = [];
  for (let t = 0; t < 256; t++) {
    const { alpha, beta } = responses(t);
    const ok = PAIRS.every((_, ip) => {
      const [ai, bi] = [Math.floor(ip / 2), ip % 2];
      const q0 = alpha[ai][0] * beta[bi][0], q1 = alpha[ai][1] * beta[bi][1];
      return target[ip] >= Math.min(q0, q1) - 1e-12 && target[ip] <= Math.max(q0, q1) + 1e-12;
    });
    if (ok) hits.push({ t, alpha, beta });
  }
  const settingFree = hits.filter((h) => h.alpha[0][0] === h.alpha[1][0] && h.alpha[0][1] === h.alpha[1][1] && h.beta[0][0] === h.beta[1][0] && h.beta[0][1] === h.beta[1][1]);
  console.log(`   response tables that can hit all four targets : ${hits.length} of 256`);
  console.log(`   ... of which the RESPONSES DO NOT DEPEND ON THE SETTINGS AT ALL : ${settingFree.length}`);
  console.log('   (in those solutions the whole correlation is carried by p(s|a,b): the source\'s bookkeeping)');
  console.log('');
  // pick the simplest: alpha(s=0)=+1, alpha(s=1)=-1, beta(s)= -1
  const alpha = [[1, -1], [1, -1]], beta = [[-1, -1], [-1, -1]];
  console.log('   THE SIMPLEST SOLUTION (all physics in the source, none in the responses):');
  console.log('        A = r * alpha(s),   alpha(s=0) = +1, alpha(s=1) = -1      (independent of a)');
  console.log('        B = r * beta(s),    beta(s)  = -1                          (independent of b)');
  console.log(`        so q0 = alpha(0)*beta(0) = ${alpha[0][0] * beta[0][0]}  and  q1 = alpha[0][1] * beta[0][1] = ${alpha[0][1] * beta[0][1]}`);
  console.log('        E = p*(-1) + (1-p)*(+1) = 1 - 2p    =>    p(s=0|a,b) = (1 - E)/2');
  console.log('');
  console.log('     setting pair      target E(a,b)     required p(s=0)      achieved E          deviation');
  let worst = 0;
  PAIRS.forEach(([a, b], ip) => {
    const p = (1 - target[ip]) / 2;
    const E = p * (-1) + (1 - p) * (+1);
    worst = Math.max(worst, Math.abs(E - target[ip]));
    console.log(`     (a=${String(a).padStart(2)}°,b=${String(b).padStart(3)}°)     ${target[ip].toFixed(6).padStart(11)}     ${p.toFixed(6).padStart(12)}     ${E.toFixed(6).padStart(12)}     ${Math.abs(E - target[ip]).toExponential(1)}`);
  });
  const E = PAIRS.map((_, ip) => { const p = (1 - target[ip]) / 2; return p * (-1) + (1 - p) * (+1); });
  console.log('');
  console.log(`     S = ${S_of(E).toFixed(6)}    (2*sqrt(2) = ${(2 * Math.SQRT2).toFixed(6)})     worst deviation ${worst.toExponential(1)}`);
  console.log('');
  console.log('   WHY IT IS EXACT: the required p values are (1 -+ 1/sqrt(2))/2 = 0.853553 and 0.146447, and');
  console.log('   a two-branch mixture realises ANY value in [-1,1] exactly -- there is no approximation step.');
  console.log('');
  console.log('   AND NOTE: p(s=0|a,b) = (1 - E)/2 = (1 + cos D)/2 = cos^2(D/2).  The source\'s bias for one bit');
  console.log('   IS cos^2(D/2).  That is the whole conspiracy.');
}
console.log('');

console.log('============================================================');
console.log('4. NO-SIGNALLING, CHECKED EXPLICITLY (NOT ASSERTED)');
console.log('============================================================');
{
  const alpha = [1, -1], beta = [-1, -1];
  console.log('   build the joint distribution: p(a,b,s,r) = (1/4)*(1/4)*p(s|a,b)*(1/2), outcomes A = r*alpha(s), B = r*beta(s)');
  console.log('');
  console.log('     setting pair     p(A=+1|a,b)     p(B=+1|a,b)');
  PAIRS.forEach(([a, b], ip) => {
    const p = (1 - target[ip]) / 2;
    let pA = 0, pB = 0;
    for (const s of [0, 1]) {
      const ps = s === 0 ? p : 1 - p;
      for (const r of [1, -1]) {
        const w = ps * 0.5;
        if (r * alpha[s] === +1) pA += w;
        if (r * beta[s] === +1) pB += w;
      }
    }
    console.log(`     (a=${String(a).padStart(2)}°,b=${String(b).padStart(3)}°)      ${pA.toFixed(12)}      ${pB.toFixed(12)}`);
  });
  console.log('');
  console.log('   every marginal is exactly 1/2: no-signalling holds for every pair, and it holds because the');
  console.log('   shared sign r multiplies BOTH outcomes -- it cancels out of every marginal.');
}
console.log('');

console.log('============================================================');
console.log('5. SO WHY 2*sqrt(2) AND NOT 4?');
console.log('============================================================');
{
  const PR = [1, 1, 1, -1];
  console.log('   the SAME machine can realise the PR box:');
  console.log(`     target E = ${PR.join(', ')}  =>  p(s=0|a,b) = ${PR.map((e) => ((1 - e) / 2).toFixed(3)).join(', ')}`);
  console.log(`     S = ${S_of(PR).toFixed(6)}`);
  console.log('');
  console.log('   => the machine reaches any target tuple in [-1,1]^4 that its response table can cover.');
  console.log('      2*sqrt(2) is simply the S value OF THE TUPLE WE ASKED FOR (-cos of the four optimal');
  console.log('      angles).  The model does not KNOW 2*sqrt(2); it is TOLD the four target correlations.');
  console.log('      That is exactly why round 72 could only retrodict: one tuning per setting pair.');
  console.log('');
  console.log('   => and that is what round 73 fixed: one FIXED density, no per-pair tuning, so it predicts');
  console.log('      angles it was never given (at 2.828698 on an 8192-bin grid, exact in the continuum).');
}
console.log('');

console.log('============================================================');
console.log('6. THE CONTROL: PUT MEASUREMENT INDEPENDENCE BACK');
console.log('============================================================');
{
  let best = -Infinity;
  for (let t = 0; t < 256; t++) {
    const { alpha, beta } = responses(t);
    for (const ps of [[1, 0], [0, 1]]) {
      const E = PAIRS.map((_, ip) => { const [ai, bi] = [Math.floor(ip / 2), ip % 2]; return ps[0] * alpha[ai][0] * beta[bi][0] + ps[1] * alpha[ai][1] * beta[bi][1]; });
      best = Math.max(best, S_of(E));
    }
  }
  console.log(`   with p(s|a,b) = p(s) (the same distribution for every pair) the maximum is S = ${best.toFixed(6)}`);
  console.log('   => the free lunch is move 1; the violation is bought entirely by move 2.');
  console.log('');
  console.log('   IN ONE SENTENCE: the model is local and deterministic, its randomness is a shared sign that');
  console.log('   cancels from every marginal, and ALL of the Bell violation lives in the source knowing the');
  console.log('   observers\' settings.  No signal travels between the wings; the correlation was in the');
  console.log('   initial conditions.');
}
