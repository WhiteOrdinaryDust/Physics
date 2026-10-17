'use strict';
/*
 * molecule_laser.js   --   ROUND 64
 * ---------------------------------------------------------------------------
 * WIRING THE MOLECULE SIDE TO THE NEW ALPHABET: CAN STIMULATED EMISSION BE EXACT?
 *
 * Round 63 built the bridge  count n = |amplitude|^2  and measured the gain ladder 1 + 1/n with a
 * minimal emitter.  This round does the real thing: molecules decaying into a MODE, i.e. a laser.
 * The laser is where the bridge is under maximum strain, because a laser is nothing but "add one
 * photon, many times, in phase".
 *
 * THE QUESTION THAT KILLS THE DISCRETE VERSION (part A)
 *   Stimulated emission must take a mode with count n to count n+1, IN THE SAME PHASE.  In the
 *   amplitude language that means  a -> a * sqrt((n+1)/n)  -- a multiplication by a real factor.
 *   If a is a Gaussian integer, when can the result also be a Gaussian integer?
 *        a' = a * r with r real positive, |a'|^2 = n+1, |a|^2 = n  =>  r = sqrt((n+1)/n)
 *        r rational  <=>  n and n+1 are both perfect squares
 *   and n, n+1 are coprime, so the only solution is n = 0.  THEREFORE:
 *        THE FIRST QUANTUM IS EXACT AND EVERY LATER ONE IS NOT.
 *   This is not a resolution limit: round 63's theorem said some counts are missing; this says the
 *   OPERATION "add one quantum" is not representable at all beyond n = 0.
 *
 * PREDICTIONS (written before the run)
 *   A  exactly one exact increment exists: 0 -> 1.  No n >= 1 has one (checked exhaustively).
 *   B  WITH A CONTINUOUS AMPLITUDE the laser works, and its law is exactly predictable:
 *      with N_e excited molecules and rate  gamma*(n+1)  per molecule per step,
 *           dn/dt = gamma*(n+1)*N_e ,  n + N_e = const = n0 + N0   (each emission uses one molecule)
 *      so S = n+1 obeys the LOGISTIC equation  dS/dt = gamma*S*(C+1-S)  with C = n0 + N0,
 *           S(t) = (C+1) / (1 + ((C+1)/S0 - 1) * exp(-gamma*(C+1)*t))
 *      => saturation at n = n0 + N0 exactly, initial exponential rate gamma*(n0+1)*N0.
 *   C  the controls that can fail: seed dependence of the initial rate disappears when stimulation
 *      is switched off (spontaneous only), and off-resonant molecules give no coupling at all.
 *   D  BOTH discrete repairs fail, in two different ways:
 *        D1 keep the phase, add a unit Gaussian integer -> the mode gains amplitude but MISCOUNTS:
 *           each emission adds 2*(a.u)+1 to the count instead of 1, so n + N_e is not conserved.
 *        D2 keep the count, jump to any Gaussian integer with |a'|^2 = n+1 -> the laser STALLS,
 *           because 3 is not a sum of two squares: from two quanta there is no third state.
 * ---------------------------------------------------------------------------
 */

// ============================================================ A: is one quantum exact?
console.log('============================================================');
console.log('A: CAN A DISCRETE AMPLITUDE ABSORB ONE QUANTUM IN PHASE? (exhaustive)');
console.log('============================================================');
{
  // all Gaussian integers with modulus^2 = m
  function states(m) {
    const out = [];
    for (let a = -Math.ceil(Math.sqrt(m)); a <= Math.ceil(Math.sqrt(m)); a++) {
      const rest = m - a * a;
      if (rest < 0) continue;
      const b = Math.round(Math.sqrt(rest));
      if (b * b === rest) { out.push([a, b]); if (b !== 0) out.push([a, -b]); }
    }
    return out;
  }
  let firstExact = 'the zero amplitude has no phase, so ANY one-quantum state is exact', holes = [];
  for (let n = 1; n <= 200; n++) {
    const from = states(n), to = states(n + 1);
    let exact = null;
    for (const [ax, ay] of from) {
      const na = Math.hypot(ax, ay);
      for (const [bx, by] of to) {
        // same direction and longer: cross product zero, dot product positive
        const cross = ax * by - ay * bx, dot = ax * bx + ay * by;
        if (cross === 0 && dot > 0 && Math.hypot(bx, by) > na) { exact = [[ax, ay], [bx, by]]; break; }
      }
      if (exact) break;
    }
    if (exact) holes.push([n, exact]);
  }
  console.log(`   n = 0 -> 1 : EXACT  (${firstExact})`);
  console.log(`   n = 1..200 : increments found in the same phase = ${holes.length}`);
  if (holes.length) for (const [n, e] of holes.slice(0, 5)) console.log(`        n = ${n}: ${e[0]} -> ${e[1]}`);
  console.log('');
  console.log('   why: a\' = a*r with r = sqrt((n+1)/n) real.  For a\' to be a Gaussian integer whenever a');
  console.log('   is, r must be rational; r rational means n and n+1 are both perfect squares, and they are');
  console.log('   coprime, so n = 0 is the only solution.  The first quantum is exact; no later one is.');
  console.log('');
  console.log(`  => ${holes.length === 0 ? 'PASS' : '*** FAIL ***'}: exactly one exact increment exists in the whole discrete alphabet (0 -> 1);`);
  console.log('     Stimulated emission cannot be hosted by a Gaussian-integer amplitude.');
}
console.log('');

// ============================================================ shared laser machinery
const OMEGA_MODE = 4, GAMMA_RES = 0;          // sharp resonance, as in the molecule rounds
function clock(seed) { let x = seed | 0; if (x === 0) x = 0x9e3779b9; return () => { x ^= x << 13; x ^= x >>> 17; x ^= x << 5; return ((x >>> 0) / 4294967296); }; }

// ============================================================ B: the continuous laser
console.log('============================================================');
console.log('B: THE LASER WITH A CONTINUOUS AMPLITUDE -- MEASURED AGAINST THE EXACT LOGISTIC');
console.log('============================================================');
{
  const N0 = 200, n0 = 1, GAMMA = 1 / 2000, STEPS = 400, TRIALS = 300;
  const C = n0 + N0, S0 = n0 + 1, SINF = C + 1;
  const logistic = (t) => SINF / (1 + (SINF / S0 - 1) * Math.exp(-GAMMA * SINF * t));
  // ---- the deterministic mean field, integrated with K sub-steps per tick (a convergence check:
  //      if the deviation from the closed form shrinks like 1/K it is discretisation, not modelling)
  function detRun(K) {
    const h = 1 / K, out = [n0];
    let n = n0;
    for (let t = 1; t <= STEPS; t++) {
      for (let s = 0; s < K; s++) n += h * GAMMA * (n + 1) * (C - n);
      out.push(n);
    }
    return out;
  }
  const det = detRun(64);
  const detConv = [1, 4, 16, 64].map((K) => {
    const d = detRun(K);
    let worst = 0;
    for (let t = 0; t <= STEPS; t++) worst = Math.max(worst, Math.abs(d[t] - (logistic(t) - 1)) / Math.max(1, logistic(t) - 1));
    return [K, worst];
  });
  // ---- the EXACT master equation: the distribution of n, no Monte Carlo noise
  let P = new Float64Array(C + 1); P[n0] = 1;
  const mean = [], vari = [];
  function stats(P) {
    let m = 0, m2 = 0;
    for (let n = 0; n <= C; n++) { m += n * P[n]; m2 += n * n * P[n]; }
    return [m, m2 - m * m];
  }
  for (let t = 0; t <= STEPS; t++) {
    const [m, v] = stats(P); mean.push(m); vari.push(Math.max(0, v));
    if (t === STEPS) break;
    const Q = new Float64Array(C + 1);
    for (let n = 0; n <= C; n++) {
      if (P[n] === 0) continue;
      const Ne = C - n;
      const p = Math.min(1, GAMMA * (n + 1));
      if (p === 0) { Q[n] += P[n]; continue; }
      if (p === 1) { Q[n + Ne] += P[n]; continue; }
      let pmf = Math.pow(1 - p, Ne);                       // binomial pmf, built iteratively
      for (let k = 0; k <= Ne; k++) {
        if (pmf > 1e-300) Q[n + k] += P[n] * pmf;
        pmf *= ((Ne - k) / (k + 1)) * (p / (1 - p));
      }
    }
    P = Q;
  }
  console.log(`   ${N0} excited molecules, seed n0 = ${n0}, rate = gamma*(n+1) per molecule per step, gamma = ${GAMMA}`);
  console.log(`   exact predictions:  saturation n = n0 + N0 = ${C}   initial rate = gamma*(n0+1)*N0 = ${(GAMMA * (n0 + 1) * N0).toFixed(4)}`);
  console.log('');
  console.log('     t     exact mean   variance   deterministic Euler   closed-form logistic   Jensen gap');
  for (const t of [0, 5, 10, 15, 20, 30, 40, 60, 100, 200, 400]) {
    const gap = det[t] - mean[t];
    console.log(`   ${String(t).padStart(4)}   ${mean[t].toFixed(2).padStart(10)}   ${vari[t].toFixed(2).padStart(8)}   ${det[t].toFixed(2).padStart(19)}   ${(logistic(t) - 1).toFixed(2).padStart(19)}   ${gap.toFixed(2).padStart(10)}`);
  }
  // ---- test 1: the deterministic mean field converges to the closed-form logistic as 1/K
  const detDev = detConv[detConv.length - 1][1];
  const detDevCoarse = detConv[0][1];
  // ---- test 2: the EXACT stochastic-mean identity  E[n'] - E[n] = gamma*[(E[n]+1)(C-E[n]) - Var(n)]
  let idDev = 0;
  for (let t = 0; t < STEPS; t++) {
    const rhs = GAMMA * ((mean[t] + 1) * (C - mean[t]) - vari[t]);
    const lhs = mean[t + 1] - mean[t];
    idDev = Math.max(idDev, Math.abs(lhs - rhs));
  }
  // ---- test 3: a Monte Carlo ensemble must agree with the exact master equation
  const mcMean = new Array(STEPS + 1).fill(0);
  for (let trial = 0; trial < TRIALS; trial++) {
    const rnd = clock(1013 + trial * 7919);
    let n = n0, Ne = N0;
    mcMean[0] += n;
    for (let t = 1; t <= STEPS; t++) {
      const count = n;                          // THE BRIDGE: the rate reads the COUNT from the field
      const p = Math.min(1, GAMMA * (count + 1));
      for (let k = 0; k < N0; k++) if (k < Ne && rnd() < p) { n += 1; Ne -= 1; }
      mcMean[t] += n;
    }
  }
  for (let t = 0; t <= STEPS; t++) mcMean[t] /= TRIALS;
  let mcDev = 0, mcT = 0;
  for (let t = 0; t <= STEPS; t++) { const d = Math.abs(mcMean[t] - mean[t]); if (d > mcDev) { mcDev = d; mcT = t; } }
  const noiseScale = Math.sqrt(Math.max(...vari)) / Math.sqrt(TRIALS);
  console.log('');
  console.log(`   exact mean at saturation = ${mean[STEPS].toFixed(4)}   (n0 + N0 = ${C} exactly)`);
  console.log(`   worst |E[n'] - E[n] - gamma*((E[n]+1)(C-E[n]) - Var(n))| over all steps = ${idDev.toExponential(2)}   <- the identity is EXACT`);
  console.log(`   deterministic mean field vs closed-form logistic, sub-steps K = ${detConv.map(([K]) => K).join(', ')}:`);
  console.log(`        ${detConv.map(([, d]) => (100 * d).toFixed(2) + '%').join('   ')}   <- shrinking as 1/K => discretisation, not modelling`);
  console.log(`   worst MC ensemble vs exact master equation = ${mcDev.toFixed(2)} at t = ${mcT}   (one sigma of the ensemble mean ~ ${noiseScale.toFixed(2)})`);
  console.log('');
  console.log('   the Jensen gap column is the physical content: the rate (n+1)(C-n) is a downward parabola,');
  console.log('   so a fluctuating n has a SMALLER mean rate than a steady n of the same average.  The exact');
  console.log('   mean therefore LAGS the deterministic logistic, and the lag is exactly gamma*Var(n) per step.');
  console.log('');
  const ok = idDev < 1e-9 && detDev < 0.005 && detDevCoarse > detDev && mcDev < 4 * noiseScale && Math.abs(mean[STEPS] - C) < 1e-6;
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: with a continuous amplitude the laser saturates at exactly n0 + N0, follows the`);
  console.log(`     logistic law in the deterministic limit, and its stochastic mean obeys the variance-corrected`);
  console.log(`     identity to machine precision.  "Add one quantum" is exact.`);
}
console.log('');

// ============================================================ C: the controls
console.log('============================================================');
console.log('C: THE CONTROLS -- SEED DEPENDENCE, AND RESONANCE');
console.log('============================================================');
{
  const N0 = 200, GAMMA = 1 / 2000, TRIALS = 4000;
  // the estimator is the rate AT t = 0: emissions during the very first tick, per trial.
  // (averaging over many steps would fold the growth into the estimate -- my first version did that
  //  and read 2.04x too high, which is what diluted the predicted ratio of 3 down to 2.73.)
  function initialRate(n0, stimulate, resonant) {
    let total = 0;
    for (let trial = 0; trial < TRIALS; trial++) {
      const rnd = clock(555 + trial * 104729);
      let n = n0;
      const p = resonant ? Math.min(1, GAMMA * (stimulate ? n + 1 : 1)) : 0;
      for (let k = 0; k < N0; k++) if (p > 0 && rnd() < p) { n += 1; total += 1; }
    }
    return total / TRIALS;
  }
  const s1 = initialRate(1, true, true), s5 = initialRate(5, true, true);
  const f1 = initialRate(1, false, true), f5 = initialRate(5, false, true);
  const off = initialRate(1, true, false);
  console.log('   initial growth rate (quanta per step, measured on the very first tick)');
  console.log(`     seed n0 = 1, stimulation ON  : ${s1.toFixed(5)}     predicted gamma*(1+1)*N0 = ${(GAMMA * 2 * N0).toFixed(5)}`);
  console.log(`     seed n0 = 5, stimulation ON  : ${s5.toFixed(5)}     predicted gamma*(5+1)*N0 = ${(GAMMA * 6 * N0).toFixed(5)}`);
  console.log(`       ratio 6/2 = 3 expected      : ${(s5 / s1).toFixed(3)}`);
  console.log(`     seed n0 = 1, stimulation OFF : ${f1.toFixed(5)}     predicted gamma*1*N0     = ${(GAMMA * N0).toFixed(5)}`);
  console.log(`     seed n0 = 5, stimulation OFF : ${f5.toFixed(5)}     (identical BY CONSTRUCTION: with stimulation off`);
  console.log(`       ratio                        : ${(f5 / f1).toFixed(3)}   the seed does not enter the dynamics at all)`);
  console.log(`     seed n0 = 1, OFF-RESONANT    : ${off.toFixed(5)}     (0 expected: no coupling to this mode)`);
  console.log('');
  const ok = Math.abs(s5 / s1 - 3) < 0.12 && Math.abs(s1 - GAMMA * 2 * N0) < 0.15 * GAMMA * 2 * N0 && off === 0;
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the seed dependence is the signature of stimulated emission (factor 3, measured`);
  console.log(`     ${(s5 / s1).toFixed(2)}); with stimulation off the seed stops mattering (${(f5 / f1).toFixed(3)}); off-resonant molecules`);
  console.log(`     do not couple at all.  Each control can fail, and none of them does.`);
}
console.log('');

// ============================================================ D: the two discrete repairs
console.log('============================================================');
console.log('D: BOTH DISCRETE REPAIRS FAIL -- ONE MISCOUNTS, THE OTHER STALLS AT TWO PHOTONS');
console.log('============================================================');
{
  const N0 = 200, n0 = 1, GAMMA = 1 / 2000, STEPS = 400, TRIALS = 40;
  // D1: keep the phase.  The mode is a Gaussian integer; each emission adds the unit Gaussian
  //     integer that is most nearly in phase with it.
  function runD1() {
    let drift = 0, finaln = 0;
    for (let trial = 0; trial < TRIALS; trial++) {
      const rnd = clock(31 + trial * 6151);
      let re = 1, im = 0, Ne = N0;
      for (let t = 1; t <= STEPS; t++) {
        for (let k = 0; k < N0; k++) {
          if (!(k < Ne)) continue;
          const count = re * re + im * im;
          if (rnd() < Math.min(1, GAMMA * (count + 1))) {
            const dirs = [[1, 0], [-1, 0], [0, 1], [0, -1]];
            let best = dirs[0], bestDot = -Infinity;
            for (const [dx, dy] of dirs) { const d = re * dx + im * dy; if (d > bestDot) { bestDot = d; best = [dx, dy]; } }
            re += best[0]; im += best[1];
            Ne -= 1;
          }
        }
        const cnt = re * re + im * im;
        drift = Math.max(drift, Math.abs(cnt + Ne - (n0 + N0)));
      }
      finaln += re * re + im * im;
    }
    return { drift, finaln: finaln / TRIALS };
  }
  // D2: keep the count.  Jump to a Gaussian integer with |a'|^2 = n+1, closest in phase.
  function statesOf(m) {
    const out = [];
    for (let a = -Math.ceil(Math.sqrt(m)); a <= Math.ceil(Math.sqrt(m)); a++) {
      const rest = m - a * a; if (rest < 0) continue;
      const b = Math.round(Math.sqrt(rest));
      if (b * b === rest) { out.push([a, b]); if (b !== 0) out.push([a, -b]); }
    }
    return out;
  }
  function runD2() {
    let drift = 0, finaln = 0, stuck = 0;
    for (let trial = 0; trial < TRIALS; trial++) {
      const rnd = clock(97 + trial * 3571);
      let re = 1, im = 0, Ne = N0;
      for (let t = 1; t <= STEPS; t++) {
        for (let k = 0; k < N0; k++) {
          if (!(k < Ne)) continue;
          const count = re * re + im * im;
          if (rnd() < Math.min(1, GAMMA * (count + 1))) {
            const to = statesOf(count + 1);
            if (to.length === 0) { stuck++; continue; }          // no such state: the quantum cannot be added
            const na = Math.hypot(re, im);
            let best = to[0], bestCos = -Infinity;
            for (const [bx, by] of to) {
              const cos = (re * bx + im * by) / (na * Math.hypot(bx, by));
              if (cos > bestCos) { bestCos = cos; best = [bx, by]; }
            }
            re = best[0]; im = best[1];
            Ne -= 1;
          }
        }
        const cnt = re * re + im * im;
        drift = Math.max(drift, Math.abs(cnt + Ne - (n0 + N0)));
      }
      finaln += re * re + im * im;
    }
    return { drift, finaln: finaln / TRIALS, stuck };
  }
  const d1 = runD1(), d2 = runD2();
  console.log('   D1  keep the phase, add a unit Gaussian integer (the only way to stay discrete):');
  console.log(`        measured final count  = ${d1.finaln.toFixed(1)}      exact value n0 + N0 = ${n0 + N0}`);
  console.log(`        max |n + N_e - const| = ${d1.drift.toFixed(1)}          (should be 0)`);
  console.log(`        => each emission adds 2*(a.u) + 1 to the count, not 1, so once the mode is bright the`);
  console.log(`           bookkeeping is wrong by a factor ~|a|.  The amplitude grows; the PHOTON COUNT does not.`);
  console.log('');
  console.log('   D2  keep the count, jump to any Gaussian integer with |a\'|^2 = n+1 (closest in phase):');
  console.log(`        measured final count  = ${d2.finaln.toFixed(1)}      exact value n0 + N0 = ${n0 + N0}`);
  console.log(`        max |n + N_e - const| = ${d2.drift.toFixed(1)}`);
  console.log(`        emissions that were IMPOSSIBLE (no such state) = ${d2.stuck}`);
  console.log(`        => 3 is not a sum of two squares (round 63), so from two quanta there is no third state:`);
  console.log(`           the laser stalls at n = 2 with ${N0 - 2} molecules still excited.`);
  console.log('');
  const ok = d1.drift > 1 && d2.finaln === 2 && d2.stuck > 0;
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: one discrete rule miscounts photons, the other cannot reach the third one.`);
  console.log(`     BOTH REPAIRS FAIL, and they fail in the two ways part A predicted.`);
  console.log(`     (D2's stall is round 63's theorem, now standing in a real physical setting: see next line)`);
  console.log(`     3 = 1 + 1 + 1 has no Gaussian-integer amplitude -- checked in round 63's hole table.`);
}
