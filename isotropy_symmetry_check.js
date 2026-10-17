'use strict';
/*
 * isotropy_symmetry_check.js  --  ROUND 106: FINDING THE SYMMETRY-BREAKING BIAS
 * ---------------------------------------------------------------------------
 * Round 105 ended with a contradiction rather than a result.  With random per-rotor frames and strong
 * coupling, the ensemble average of M4 came out at 0.2651 +- 0.0233 (11.4 sigma), but the theorem says it
 * must be exactly zero: the frame ensemble is rotation-invariant, and the only rotation-invariant
 * probability measure on the sphere is the uniform one, for which M4 vanishes identically.
 *
 * So the pipeline has a bias that breaks rotation invariance, and this round locates it by cutting the
 * pipeline into three stages and measuring M4 on each.  The theorem constrains the first two to zero, so
 * the first stage that is NOT zero is the culprit:
 *
 *   STAGE 1  the frame generator itself: the directions of the frames' own cubic axes.  If the axis
 *            distribution is uniform, their M4 must average zero.
 *   STAGE 2  the peak seeding: the argmax of  S(u) = sum_i sum_a (u . e_a(i))^4 , found by a dense scan
 *            plus projected gradient ascent.  Since each S_i is built from a random frame, S is
 *            rotation-COVARIANT in distribution, so the argmax direction must be invariant and its M4
 *            must average zero.
 *   STAGE 3  the self-consistent iteration: its fixed point.  Round 105 measured 0.2651 here.
 *
 * AND THE DECISIVE EXTRA NUMBER: the ANGLE between each seed and the fixed point it converges to.  If the
 * iteration barely moves, the bias is inherited from the seeding; if it moves a lot, the bias is created
 * by the iteration itself.
 *
 * PRE-REGISTERED PREDICTIONS
 *   Z1 stage 1 averages zero (the generator is uniform)
 *   Z2 stage 2 averages zero (the landscape is covariant)
 *   Z3 the nonzero value appears at stage 3, and it comes with a large seed-to-fixed-point angle, or with
 *      a small one -- either answer localises the bias, which is the point of the round
 * ---------------------------------------------------------------------------
 */

const KAPPA = 10, N = 64;
function stream(seed) { let s = seed | 0; if (s === 0) s = 0x9e3779b9; return () => { s ^= s << 13; s ^= s >>> 17; s ^= s << 5; return (s >>> 0) / 4294967296; }; }
function gauss(rnd) { let u = 0, v = 0; while (u === 0) u = rnd(); while (v === 0) v = rnd(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
// HAAR-UNIFORM ROTATION, and the old version was the bug.  It used a uniform axis with a uniform ANGLE,
// which is not the Haar measure: the angle has too much weight near 0 and 2*pi (small rotations), so the
// rotated columns clustered near the lab coordinate axes -- the stage-1 check measured M4 = 0.0448 +- 0.0047
// for the frames' own axes, i.e. a 9.5 sigma bias in the GENERATOR, and every later stage inherited it.
// A uniform unit quaternion is exactly Haar, so this replaces it.
function randomRotation(rnd) {
  let q = [gauss(rnd), gauss(rnd), gauss(rnd), gauss(rnd)];
  const n = Math.hypot(q[0], q[1], q[2], q[3]) || 1;
  const w = q[0] / n, x = q[1] / n, y = q[2] / n, z = q[3] / n;
  return [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y),
          2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
          2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)];
}
function uniformSphere(rnd) { const z = 1 - 2 * rnd(), ph = 2 * Math.PI * rnd(), r = Math.sqrt(Math.max(0, 1 - z * z)); return [r * Math.cos(ph), r * Math.sin(ph), z]; }
const cubicSum = (v, ax) => { let s = 0; for (const a of ax) { const d = v[0] * a[0] + v[1] * a[1] + v[2] * a[2]; s += d * d * d * d; } return s; };
const M4of = (u) => { const n = Math.hypot(u[0], u[1], u[2]) || 1; return (u[0] / n) ** 4 + (u[1] / n) ** 4 + (u[2] / n) ** 4 - 0.6; };
const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;
const seOf = (a) => { const m = mean(a); return Math.sqrt(a.reduce((x, y) => x + (y - m) ** 2, 0) / Math.max(1, a.length - 1) / a.length); };
function makeFrames(n, rnd) {
  const ax = [];
  for (let i = 0; i < n; i++) { const R = randomRotation(rnd); ax.push([[R[0], R[3], R[6]], [R[1], R[4], R[7]], [R[2], R[5], R[8]]]); }
  return ax;
}
function ascendS(u0, ax, iters) {
  let v = u0.slice();
  const grad = (u) => {
    const g = [0, 0, 0];
    for (const a of ax) { const d = u[0] * a[0] + u[1] * a[1] + u[2] * a[2]; const cc = 4 * d * d * d; g[0] += cc * a[0]; g[1] += cc * a[1]; g[2] += cc * a[2]; }
    return g;
  };
  let best = cubicSum(v, ax), eta = 0.1;
  for (let k = 0; k < iters; k++) {
    const g = grad(v);
    const dot = g[0] * v[0] + g[1] * v[1] + g[2] * v[2];
    let p = [v[0] + eta * (g[0] - dot * v[0]), v[1] + eta * (g[1] - dot * v[1]), v[2] + eta * (g[2] - dot * v[2])];
    const n = Math.hypot(p[0], p[1], p[2]) || 1; p = [p[0] / n, p[1] / n, p[2] / n];
    const val = cubicSum(p, ax);
    if (val > best) { v = p; best = val; eta = Math.min(0.5, eta * 1.05); } else eta *= 0.5;
    if (eta < 1e-7) break;
  }
  return v;
}
function innerMax(h, axOne, rnd, samples, refine) {
  const f = (v) => (h[0] * v[0] + h[1] * v[1] + h[2] * v[2]) + KAPPA * cubicSum(v, axOne);
  let best = uniformSphere(rnd), bestF = f(best);
  for (let k = 0; k < samples; k++) { const v = uniformSphere(rnd); const val = f(v); if (val > bestF) { best = v; bestF = val; } }
  let v = best.slice(), eta = 0.25;
  const grad = (u) => {
    const g = [h[0], h[1], h[2]];
    for (const a of axOne) { const d = u[0] * a[0] + u[1] * a[1] + u[2] * a[2]; const cc = 4 * KAPPA * d * d * d; g[0] += cc * a[0]; g[1] += cc * a[1]; g[2] += cc * a[2]; }
    return g;
  };
  for (let k = 0; k < refine; k++) {
    const g = grad(v);
    const dot = g[0] * v[0] + g[1] * v[1] + g[2] * v[2];
    let p = [v[0] + eta * (g[0] - dot * v[0]), v[1] + eta * (g[1] - dot * v[1]), v[2] + eta * (g[2] - dot * v[2])];
    const n = Math.hypot(p[0], p[1], p[2]) || 1; p = [p[0] / n, p[1] / n, p[2] / n];
    const val = f(p);
    if (val > bestF) { v = p; bestF = val; eta = Math.min(0.5, eta * 1.05); } else eta *= 0.5;
    if (eta < 1e-7) break;
  }
  return { v, f: bestF };
}
function solveFrom(ax, n, J, seedDir, rnd, samples, refine, iters) {
  let M = [seedDir[0] * n, seedDir[1] * n, seedDir[2] * n];
  let lastE = Infinity;
  for (let it = 0; it < iters; it++) {
    const h = [J * M[0] / n, J * M[1] / n, J * M[2] / n];
    let sx = 0, sy = 0, sz = 0;
    for (let i = 0; i < n; i++) { const r = innerMax(h, ax[i], rnd, samples, refine); sx += r.v[0]; sy += r.v[1]; sz += r.v[2]; }
    const newM = [sx, sy, sz];
    const conv = Math.hypot(newM[0] - M[0], newM[1] - M[1], newM[2] - M[2]) / n;
    M = newM;
    lastE = -(J / (2 * n)) * (M[0] * M[0] + M[1] * M[1] + M[2] * M[2]) - KAPPA * ax.reduce((a, one, i) => a, 0);
    if (conv < 1e-6) break;
  }
  const nrm = Math.hypot(M[0], M[1], M[2]) || 1;
  return { u: [M[0] / nrm, M[1] / nrm, M[2] / nrm], e: lastE };
}
const angle = (a, b) => Math.acos(Math.max(-1, Math.min(1, a[0] * b[0] + a[1] * b[1] + a[2] * b[2]))) * 180 / Math.PI;

console.log('============================================================');
console.log('STAGE BY STAGE: where does the rotation-invariance break?');
console.log('============================================================');
console.log('   the theorem: an invariant frame ensemble forces every stage below to average ZERO in M4.');
console.log('');

// STAGE 1
{
  const rnd = stream(11);
  const vals = [];
  const AX = makeFrames(300, rnd);
  for (const ax of AX) for (const a of ax) for (const sg of [1, -1]) vals.push(M4of([a[0] * sg, a[1] * sg, a[2] * sg]));
  console.log(`   STAGE 1  the frames' own axis directions  (${vals.length} directions): M4 = ${mean(vals).toFixed(4)} +- ${seOf(vals).toFixed(4)}`);
  console.log(`            => Z1 ${Math.abs(mean(vals)) < 3 * seOf(vals) ? 'CONFIRMED: the generator is uniform' : 'FAILS: the GENERATOR is biased'}`);
}
console.log('');

// STAGE 2
const stage2 = [];
{
  const rnd = stream(22);
  console.log('   STAGE 2  the argmax of sum_i S_i(u): a dense scan of 20000 directions plus gradient ascent');
  for (let k = 0; k < 30; k++) {
    const ax = makeFrames(N, rnd);
    let best = uniformSphere(rnd), bestS = cubicSum(best, ax);
    for (let t = 0; t < 20000; t++) { const u = uniformSphere(rnd); const s = cubicSum(u, ax); if (s > bestS) { best = u; bestS = s; } }
    const refined = ascendS(best, ax, 400);
    stage2.push(M4of(refined));
  }
  console.log(`            M4 of the argmax over ${stage2.length} frame ensembles: ${mean(stage2).toFixed(4)} +- ${seOf(stage2).toFixed(4)}`);
  console.log(`            => Z2 ${Math.abs(mean(stage2)) < 3 * seOf(stage2) ? 'CONFIRMED: the landscape is covariant and its argmax is isotropic' : 'FAILS: the PEAK SEARCH is biased'}`);
}
console.log('');

// STAGE 3
{
  const rnd = stream(33);
  console.log('   STAGE 3  the self-consistent iteration from the top peaks, plus the seed-to-fixed-point angle');
  const vals = [], angs = [];
  for (let k = 0; k < 20; k++) {
    const ax = makeFrames(N, rnd);
    const seeds = [];
    for (let t = 0; t < 4000; t++) { const u = uniformSphere(rnd); seeds.push({ u, s: cubicSum(u, ax) }); }
    seeds.sort((a, b) => b.s - a.s);
    const top = seeds.slice(0, 8).map((p) => p.u);
    let best = null;
    for (const s of top) {
      const r = solveFrom(ax, N, 25, s, rnd, 400, 120, 25);
      if (!best || r.e < best.e) best = { ...r, seed: s };
    }
    vals.push(M4of(best.u));
    angs.push(angle(best.seed, best.u));
  }
  console.log(`            M4 of the fixed point: ${mean(vals).toFixed(4)} +- ${seOf(vals).toFixed(4)}`);
  console.log(`            seed-to-fixed-point angle: mean ${mean(angs).toFixed(1)} degrees, max ${Math.max(...angs).toFixed(1)}`);
  console.log(`            => Z3 the nonzero stage is: ${Math.abs(mean(vals)) > 3 * seOf(vals) ? 'STAGE 3' : 'none of the three'}`);
  console.log('');
  console.log(`   DIAGNOSIS: ${mean(angs) < 15 ? 'the iteration barely moves off its seed, so the bias is inherited from the SEED SELECTION (stage 2), even though stage 2 measured zero on its own argmax -- which points at the way the top-8 peaks are CHOSEN from a finite sample' : 'the iteration moves a long way from its seed, so the bias is created BY THE ITERATION itself (stage 3)'}`);
}
