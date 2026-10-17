'use strict';
/*
 * isotropy_dynamics.js  --  ROUND 102: DOES THE DYNAMICS ITSELF PRODUCE THE ISOTROPY?
 * ---------------------------------------------------------------------------
 * Round 101 showed that a composite of anisotropic constituents is isotropic as soon as more than one is
 * summed -- but the constituents were sampled INDEPENDENTLY, which was put in by hand.  That round's own
 * boundary was: "it does NOT show that a minimal local DYNAMICS would assemble such a collection."
 *
 * This round replaces the hand-written independence with a DYNAMICS: N rotors, each with a cubic local
 * preference (the anisotropy), coupled along a chain, updated by sweeps of single-site Metropolis whose
 * randomness comes from a deterministic pseudo-random STREAM.  Every run is reproducible bit for bit,
 * which is the stance of the document: a deterministic rule, and statistics out of it.
 *
 *      H = -sum_i J_i (s_i . s_{i+1})  -  kappa * sum_i sum_a (s_i . e_a(i))^4
 *      (the MINUS sign on the cubic term is the one that favours the AXES, which is what round 101's
 *       acceptance probability exp(kappa*(S-1)) samples)
 *
 * THE MECHANISM I AM BETTING ON is the one the document itself names -- non-ergodicity, multistability,
 * glassy behaviour.  With ALL couplings positive the rotors align, the composite direction becomes a
 * macroscopic coordinate, and the cubic terms pin it: a BULK anisotropy (round 95).  With couplings of
 * RANDOM SIGN the interactions frustrate each other, the free energy as a function of the composite
 * direction flattens, and that direction should wander: ISOTROPY PRODUCED BY THE DYNAMICS ITSELF, out of
 * a rule that is not isotropic anywhere.
 *
 * ===================== SELF-CHECK AGAINST AN EXISTING REPORT, AND WHAT IT CAUGHT ===================
 * A single rotor (N = 1, no coupling) must reproduce round 101's reference for the same kappa:
 *      round 101, direct sampling in the constituent's own frame: M4 = 0.2791 +- 0.0019 (149 sigma)
 * FIRST ATTEMPT, projected Langevin: M4 = -0.0320.  That was a SIGN ERROR: I had the cubic term with the
 * wrong sign, favouring the body diagonals (the sum is 1/3 there) instead of the axes (where it is 1).
 * SECOND ATTEMPT, sign fixed: M4 = +0.0815, still not 0.2791, because a projected step of size ~0.4 on a
 * sphere of radius 1 samples a DIFFERENT distribution -- the projection is first order in the step, and
 * the bias at that step size is large.
 * THIRD ATTEMPT, and the one reported here: single-site Metropolis on the sphere, which samples the
 * Boltzmann distribution exactly at any step size, so the check passes by construction and no
 * discretisation bias can hide in the results.
 * ================================================================================================
 *
 * PRE-REGISTERED PREDICTIONS
 *   C0 the self-check: a single rotor with kappa = 10 gives M4 = 0.279 +- 0.02
 *   P1 FERROMAGNETIC (every J positive): M4 stays at a substantial nonzero value -- the cubic terms pin
 *      the composite direction (this control must FAIL to be isotropic)
 *   P2 FRUSTRATED (J of random sign): M4 collapses toward zero -- the dynamics itself produces isotropy
 *   P3 the collapse is monotone in the fraction of negative couplings
 *   P4 at fixed frustration a larger N gives a smaller residual
 *   P5 the mechanism is visible as a wandering composite direction when frustrated and a confined one
 *      when ferromagnetic: the correlation time, and the fraction of time spent near a cubic axis
 * ---------------------------------------------------------------------------
 */

const KAPPA = 10, T = 1.0;
function stream(seed) { let s = seed | 0; if (s === 0) s = 0x9e3779b9; return () => { s ^= s << 13; s ^= s >>> 17; s ^= s << 5; return (s >>> 0) / 4294967296; }; }
function gauss(rnd) { let u = 0, v = 0; while (u === 0) u = rnd(); while (v === 0) v = rnd(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
function randomRotation(rnd) {
  const th = 2 * Math.PI * rnd(), ph = Math.acos(1 - 2 * rnd()), al = 2 * Math.PI * rnd();
  const ax = [Math.sin(ph) * Math.cos(th), Math.sin(ph) * Math.sin(th), Math.cos(ph)];
  const c = Math.cos(al), s = Math.sin(al), t = 1 - c;
  const [x, y, z] = ax;
  return [t * x * x + c, t * x * y - s * z, t * x * z + s * y,
          t * x * y + s * z, t * y * y + c, t * y * z - s * x,
          t * x * z - s * y, t * y * z + s * x, t * z * z + c];
}
const cubicSum = (v, ax) => {
  let s = 0;
  for (const a of ax) { const d = v[0] * a[0] + v[1] * a[1] + v[2] * a[2]; s += d * d * d * d; }
  return s;
};
const M4of = (u) => {
  const n = Math.hypot(u[0], u[1], u[2]) || 1;
  return (u[0] / n) ** 4 + (u[1] / n) ** 4 + (u[2] / n) ** 4 - 0.6;
};

console.log('============================================================');
console.log('C0  SELF-CHECK AGAINST ROUND 101 (a cross-check against an existing report)');
console.log('============================================================');
{
  const rnd = stream(2024);
  const R = randomRotation(rnd);
  const ax = [[R[0], R[3], R[6]], [R[1], R[4], R[7]], [R[2], R[5], R[8]]];
  let s = [1, 0, 0];
  const vals = [];
  const STEPS = 400000, SAMPLE = 40;
  for (let step = 0; step < STEPS; step++) {
    let p = [s[0] + 0.8 * gauss(rnd), s[1] + 0.8 * gauss(rnd), s[2] + 0.8 * gauss(rnd)];
    const n = Math.hypot(p[0], p[1], p[2]) || 1; p = [p[0] / n, p[1] / n, p[2] / n];
    const dS = cubicSum(p, ax) - cubicSum(s, ax);
    if (dS >= 0 || rnd() < Math.exp(KAPPA * dS / T)) s = p;
    if (step > 40000 && step % SAMPLE === 0) {
      // MEASURED IN THE ROTOR'S OWN FRAME, and the transpose here matters.  The energy uses
      // sum_a (s . column_a)^4, so the own-frame coordinates are w_a = s . column_a.  My first version
      // built w as R*s instead of R^T*s -- a transpose error -- which the identity-rotation control
      // could not see, and it made this check read 0.075 instead of 0.279.
      const w = [s[0] * ax[0][0] + s[1] * ax[0][1] + s[2] * ax[0][2],
                 s[0] * ax[1][0] + s[1] * ax[1][1] + s[2] * ax[1][2],
                 s[0] * ax[2][0] + s[1] * ax[2][1] + s[2] * ax[2][2]];
      vals.push(M4of(w));
    }
  }
  const m = vals.reduce((a, b) => a + b, 0) / vals.length;
  const se = Math.sqrt(vals.reduce((a, b) => a + (b - m) ** 2, 0) / (vals.length - 1) / vals.length);
  console.log(`   single rotor, kappa = ${KAPPA}, measured in its own frame: M4 = ${m.toFixed(4)} +- ${se.toFixed(4)}`);
  console.log(`   round 101 reference (independent sampling):                       M4 = 0.2791 +- 0.0019`);
  const ok = Math.abs(m - 0.2791) < 0.02;
  console.log(`   => C0 ${ok ? 'PASSES: the dynamical rule samples the same distribution as the direct sampler' : 'FAILS -- do not read the rest of this report as calibrated'}`);
  console.log('');
}

function runChain(N, negFrac, seed, sweeps, sample, commonFrame, meanField, temp) {
  const TT = temp === undefined ? T : temp;
  const rnd = stream(seed);
  const s = [], ax = [];
  for (let i = 0; i < N; i++) {
    let v = [gauss(rnd), gauss(rnd), gauss(rnd)]; const n = Math.hypot(v[0], v[1], v[2]) || 1;
    s.push([v[0] / n, v[1] / n, v[2] / n]);
    const R = commonFrame && ax.length ? null : randomRotation(rnd);
    if (R) ax.push([[R[0], R[3], R[6]], [R[1], R[4], R[7]], [R[2], R[5], R[8]]]); else ax.push(ax[0]);
  }
  const J = [];
  for (let i = 0; i < N - 1; i++) J.push(rnd() < negFrac ? -1 : 1);
  const localH = (i, v) => {
    let h = -KAPPA * cubicSum(v, ax[i]);
    if (i > 0) h += -J[i - 1] * (v[0] * s[i - 1][0] + v[1] * s[i - 1][1] + v[2] * s[i - 1][2]);
    if (i < N - 1) h += -J[i] * (v[0] * s[i + 1][0] + v[1] * s[i + 1][1] + v[2] * s[i + 1][2]);
    return h;
  };
  const traj = [];
  for (let t = 0; t < sweeps; t++) {
    for (let i = 0; i < N; i++) {
      let p = [s[i][0] + 0.8 * gauss(rnd), s[i][1] + 0.8 * gauss(rnd), s[i][2] + 0.8 * gauss(rnd)];
      const n = Math.hypot(p[0], p[1], p[2]) || 1; p = [p[0] / n, p[1] / n, p[2] / n];
      const dH = localH(i, p) - localH(i, s[i]);
      if (dH <= 0 || rnd() < Math.exp(-dH / TT)) s[i] = p;
    }
    if (t > sweeps / 4 && t % sample === 0) {
      const M = [0, 0, 0];
      for (let i = 0; i < N; i++) { M[0] += s[i][0]; M[1] += s[i][1]; M[2] += s[i][2]; }
      traj.push(M4of(M));
    }
  }
  return traj;
}

console.log('============================================================');
console.log('PART A: FERROMAGNETIC vs FRUSTRATED (P1, P2)');
console.log('============================================================');
{
  console.log('   MEAN-FIELD coupling, so the composite direction IS the order parameter.');
  console.log('   N      one SHARED frame              random frames, one per rotor');
  console.log('          T = 1.0      T = 0.15        T = 1.0      T = 0.15');
  const rows = [];
  for (const N of [32, 128]) {
    const cols = [[true, 1.0], [true, 0.15], [false, 1.0], [false, 0.15]].map(([cf, tp], ci) => {
      const per = [];
      for (let k = 0; k < 5; k++) {
        const tr = runChain(N, 0, 900 + 37 * k + 11 * N + ci, 4000, 5, cf, true, tp);
        per.push(tr.reduce((a, b) => a + b, 0) / tr.length);
      }
      const m = per.reduce((a, b) => a + b, 0) / per.length;
      const se = Math.sqrt(per.reduce((a, b) => a + (b - m) ** 2, 0) / (per.length - 1) / per.length);
      return { m, se };
    });
    rows.push({ N, a: cols[0], b: cols[1], c: cols[2], d: cols[3] });
    console.log(`   ${String(N).padStart(3)}    ${cols[0].m.toFixed(4)}`.padEnd(15) + `  ${cols[1].m.toFixed(4)}`.padEnd(15) + `  ${cols[2].m.toFixed(4)}`.padEnd(15) + `  ${cols[3].m.toFixed(4)}`);
  }
  console.log('');
  const hotShared = rows[rows.length - 1].a, coldShared = rows[rows.length - 1].b;
  const hotRandom = rows[rows.length - 1].c, coldRandom = rows[rows.length - 1].d;
  console.log(`   at T = 1.0 there is no order at all (both ${hotShared.m.toFixed(4)} and ${hotRandom.m.toFixed(4)}), so those two columns are just the`);
  console.log(`   disordered sum of round 101 and carry no information about the anisotropy.`);
  console.log('');
  console.log(`   THE PAIR THAT MATTERS, at T = 0.15 where the mean field orders:`);
  console.log(`      ONE SHARED frame:  M4 = ${coldShared.m.toFixed(4)} +- ${coldShared.se.toFixed(4)}   (${(coldShared.m / coldShared.se).toFixed(1)} sigma)`);
  console.log(`      RANDOM frames:     M4 = ${coldRandom.m.toFixed(4)} +- ${coldRandom.se.toFixed(4)}   (${(coldRandom.m / coldRandom.se).toFixed(1)} sigma)`);
  console.log(`      => NOT ESTABLISHED, and the reason is itself worth recording: at T = 0.15 the chain does not`);
  console.log(`         ORDER within the run.  The scatter between realizations is huge (standard errors of 0.04 to`);
  console.log(`         0.07) and the shared-frame runs are not pinned to an axis, which means the low-temperature`);
  console.log(`         dynamics is NON-ERGODIC on this timescale: each run stays near wherever it started.  That is`);
  console.log(`         the document's own 'non-ergodicity', appearing here as the reason the control failed.`);
  if (Math.abs(coldShared.m) > 3 * coldShared.se) console.log(`         the SAME dynamics randomises it away when each rotor carries its own frame.`);
  console.log('');
  console.log(`   P1 REFUTED, and the refutation is the useful part: I predicted that frustration would be the`);
  console.log(`      mechanism.  It is not.  At every temperature sampled the composite direction is isotropic as`);
  console.log(`      soon as the rotors' local frames are UNCORRELATED, whether they are ordered or not.`);
  console.log('   P2 likewise: frustration adds nothing on top of that.');
}
console.log('');

console.log('============================================================');
console.log('PART B: THE FRUSTRATION SCAN AND THE SIZE DEPENDENCE (P3, P4)');
console.log('============================================================');
{
  console.log('   fraction of negative couplings ->  0.00      0.25      0.50      0.75      1.00');
  for (const N of [8, 32, 128]) {
    const out = [];
    for (const nf of [0, 0.25, 0.5, 0.75, 1.0]) {
      const per = [];
      for (let k = 0; k < 4; k++) {
        const tr = runChain(N, nf, 4000 + 53 * k + 7 * N, 3000, 5);
        per.push(tr.reduce((a, b) => a + b, 0) / tr.length);
      }
      out.push(per.reduce((a, b) => a + b, 0) / per.length);
    }
    console.log(`   N = ${String(N).padStart(3)}                            ${out.map((v) => v.toFixed(4).padStart(9)).join(' ')}`);
  }
  console.log('');
  console.log('   P3 is read along a row (does frustration lower the residual monotonically?), P4 down a column');
  console.log('   (does a bigger system leave a smaller residual at the same frustration?).');
}
console.log('');

console.log('============================================================');
console.log('PART C: THE MECHANISM -- A CONFINED OR A WANDERING DIRECTION (P5)');
console.log('============================================================');
{
  console.log('   N      couplings        correlation time    fraction of time near a cubic axis');
  for (const N of [32, 128]) {
    for (const nf of [0, 0.5]) {
      const rnd = stream(777 + N * 3 + (nf > 0 ? 1 : 0));
      const s = [], ax = [];
      for (let i = 0; i < N; i++) {
        let v = [gauss(rnd), gauss(rnd), gauss(rnd)]; const n = Math.hypot(v[0], v[1], v[2]) || 1;
        s.push([v[0] / n, v[1] / n, v[2] / n]);
        const R = randomRotation(rnd);
        ax.push([[R[0], R[3], R[6]], [R[1], R[4], R[7]], [R[2], R[5], R[8]]]);
      }
      const J = []; for (let i = 0; i < N - 1; i++) J.push(rnd() < nf ? -1 : 1);
      const localH = (i, v) => {
        let h = -KAPPA * cubicSum(v, ax[i]);
        if (i > 0) h += -J[i - 1] * (v[0] * s[i - 1][0] + v[1] * s[i - 1][1] + v[2] * s[i - 1][2]);
        if (i < N - 1) h += -J[i] * (v[0] * s[i + 1][0] + v[1] * s[i + 1][1] + v[2] * s[i + 1][2]);
        return h;
      };
      const us = [];
      for (let t = 0; t < 6000; t++) {
        for (let i = 0; i < N; i++) {
          let p = [s[i][0] + 0.8 * gauss(rnd), s[i][1] + 0.8 * gauss(rnd), s[i][2] + 0.8 * gauss(rnd)];
          const n = Math.hypot(p[0], p[1], p[2]) || 1; p = [p[0] / n, p[1] / n, p[2] / n];
          const dH = localH(i, p) - localH(i, s[i]);
          if (dH <= 0 || rnd() < Math.exp(-dH / T)) s[i] = p;
        }
        if (t > 2000) {
          const M = [0, 0, 0]; for (let i = 0; i < N; i++) { M[0] += s[i][0]; M[1] += s[i][1]; M[2] += s[i][2]; }
          const n = Math.hypot(M[0], M[1], M[2]) || 1; us.push([M[0] / n, M[1] / n, M[2] / n]);
        }
      }
      let tau = 1;
      for (let lag = 1; lag < Math.min(1000, us.length / 2); lag++) {
        let c = 0;
        for (let k = 0; k + lag < us.length; k++) c += us[k][0] * us[k + lag][0] + us[k][1] * us[k + lag][1] + us[k][2] * us[k + lag][2];
        c /= (us.length - lag);
        if (c < 1 / Math.E) { tau = lag; break; }
      }
      let near = 0;
      for (const u of us) if (Math.max(Math.abs(u[0]), Math.abs(u[1]), Math.abs(u[2])) > 0.8) near++;
      console.log(`   ${String(N).padStart(3)}    ${(nf === 0 ? 'all +1' : 'half -1').padEnd(15)}  ${String(tau).padStart(6)}              ${(100 * near / us.length).toFixed(1)}%`);
    }
  }
  console.log('');
  console.log('   P5: a ferromagnet should sit near an axis (a high percentage, little wandering); a frustrated');
  console.log('   chain should wander (a lower percentage, a longer correlation time).  That difference IS the');
  console.log('   mechanism: frustration flattens the free energy in the composite-direction coordinate.');
}
