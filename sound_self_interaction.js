'use strict';
/*
 * sound_self_interaction.js   (third version -- the first two designs were broken)
 *
 * WHAT WAS WRONG
 *   v1 passed the A-run in as both pulses (no B-only run was ever generated), so the
 *      "coupling" it measured was just the difference between the two pulses. It read 100% of
 *      the pulse amplitude before the pulses could possibly have met.
 *   v2 tried to test linearity by comparing amplitude A with 2A. Its own control exposed the
 *      flaw: with collisions OFF the ratio must be EXACTLY 2 and it came out 2.76% off. The
 *      reason is structural, not numerical: the extra particles added by a probability bump are
 *      the set {u < p+A}, and {u < p+2A} is NOT twice {u < p+A} cell by cell. So "the shot noise
 *      cancels exactly" was false for that design.
 *
 * THE DESIGN THAT IS EXACT
 *   Use two perturbations with DISJOINT spatial support. Write the occupancy test as
 *   u < p_base + d1(x) + d2(x). Where the supports are disjoint, this equals
 *   (u < p_base+d1) OR (u < p_base+d2), so
 *        cfg(both) = cfg(base) + extra1 + extra2   AS SETS, exactly, for the same u-field.
 *   Therefore
 *        R(t) = d1(t) + d2(t) - d12(t)
 *   is IDENTICALLY ZERO at t=0, and it stays exactly zero until the two perturbations'
 *   influence cones meet. Any growth of R after that is the mode-mode coupling, measured with
 *   NO noise floor at all -- not a small number, but zero.
 *
 *   This is the same signature the Life glider test had: nonlinearity is not a property of the
 *   medium everywhere, it switches on at contact.
 *
 * PREDICTIONS (before this run)
 *   U1  the two smooth, long-wavelength bumps launch pulses on the SOUND branch, travelling at
 *       about 0.7 cells/step (not 1.0).
 *   U2  CONTROL with collisions OFF: R(t) == 0 exactly, at every time, to machine precision.
 *       (Free streaming cannot couple two disjoint perturbations.)
 *   U3  with collisions ON: R(t) is exactly 0 until the fronts meet, then grows. The fronts
 *       start 64 cells apart and close at ~1.4 cells/step, so the meeting time is t ~ 45.
 *   U4  the coupling is weak: after the crossing, R stays a small fraction of the pulse height,
 *       and the outgoing pulses keep most of their height (>90%).
 * ---------------------------------------------------------------------------
 */

const W = 256, H = 16, N = W * H;
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
const NB = (1 << 0) | (1 << 2), EW = (1 << 1) | (1 << 3);
const PBASE = 0.25, T = 200, SEEDS = 12, AMP = 0.05;
const B1 = [32, 96], B2 = [160, 224];        // disjoint supports
const MEET = Math.round((B2[0] - B1[1]) / (2 * 0.707));

function collide(g) { for (let i = 0; i < N; i++) { const v = g[i]; if (v === NB) g[i] = EW; else if (v === EW) g[i] = NB; } }
function propagate(g) {
  const o = new Uint8Array(N);
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
    const v = g[y * W + x];
    if (!v) continue;
    for (let d = 0; d < 4; d++) if ((v >> d) & 1) o[((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W)] |= 1 << d;
  }
  return o;
}
function step(g, doCollide) { const h = Uint8Array.from(g); if (doCollide) collide(h); return propagate(h); }
function mulberry32(a) { return function () { a |= 0; a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
// smooth long-wavelength bump: sin^2 over its window, zero outside
function bumpShape(x, a, b) {
  if (x < a || x > b) return 0;
  const s = Math.sin(Math.PI * (x - a) / (b - a));
  return s * s;
}
function makeGas(u, use1, use2, Amp) {
  const g = new Uint8Array(N);
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
    let p = PBASE;
    if (use1) p += Amp * bumpShape(x, B1[0], B1[1]);
    if (use2) p += Amp * bumpShape(x, B2[0], B2[1]);
    let v = 0;
    for (let d = 0; d < 4; d++) if (u[((y * W + x) * 4) + d] < p) v |= 1 << d;
    g[y * W + x] = v;
  }
  return g;
}
function profile(g) {
  const p = new Float64Array(W);
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) { let v = g[y * W + x], c = 0; while (v) { c += v & 1; v >>= 1; } p[x] += c; }
  for (let x = 0; x < W; x++) p[x] /= H;
  return p;
}
function runAll(u, doCollide, Amp) {
  let gb = makeGas(u, false, false, Amp), g1 = makeGas(u, true, false, Amp), g2 = makeGas(u, false, true, Amp), g12 = makeGas(u, true, true, Amp);
  const out = {};
  for (const key of ['b', '1', '2', '12']) out[key] = [];
  for (let t = 0; t <= T; t++) {
    const pb = profile(gb), p1 = profile(g1), p2 = profile(g2), p12 = profile(g12);
    const d1 = new Float64Array(W), d2 = new Float64Array(W), d12 = new Float64Array(W);
    for (let x = 0; x < W; x++) { d1[x] = p1[x] - pb[x]; d2[x] = p2[x] - pb[x]; d12[x] = p12[x] - pb[x]; }
    out.b.push(pb); out['1'].push(d1); out['2'].push(d2); out['12'].push(d12);
    gb = step(gb, doCollide); g1 = step(g1, doCollide); g2 = step(g2, doCollide); g12 = step(g12, doCollide);
  }
  return out;
}
const US = [];
for (let s = 0; s < SEEDS; s++) { const rnd = mulberry32(31337 + s * 104729); const u = new Float64Array(N * 4); for (let i = 0; i < N * 4; i++) u[i] = rnd(); US.push(u); }

// ============================================================ exactness at t=0
console.log('============================================================');
console.log('EXACTNESS CHECK (this is the control that v1 and v2 both failed)');
console.log('============================================================');
{
  const run = US.map(u => runAll(u, true, AMP));
  let worst = 0;
  for (let s = 0; s < SEEDS; s++) for (let x = 0; x < W; x++) {
    const r = run[s]['1'][0][x] + run[s]['2'][0][x] - run[s]['12'][0][x];
    if (Math.abs(r) > worst) worst = Math.abs(r);
  }
  console.log(`  worst |R| at t=0 over all seeds and cells: ${worst}   ${worst === 0 ? 'PASS (the three runs are exactly additive as particle sets)' : '*** FAIL ***'}`);
  // sanity on the perturbation size: a bump of Amp per direction should raise the density by 4*Amp at its peak
  let mxSeed = 0;
  for (let x = 0; x < W; x++) mxSeed = Math.max(mxSeed, Math.abs(run[0]['1'][0][x]));
  const mean = new Float64Array(W);
  for (const r of run) for (let x = 0; x < W; x++) mean[x] += r['1'][0][x] / run.length;
  const mxMean = Math.max(...mean.map(Math.abs));
  console.log(`  perturbation size check: expected peak ${(4 * AMP).toFixed(3)} particles/cell; measured single seed ${mxSeed.toFixed(4)}, seed average ${mxMean.toFixed(4)}`);
  console.log(`    => ${Math.abs(mxMean - 4 * AMP) < 0.05 ? 'the bump does what it says' : '*** the bump is larger than intended: the perturbations are NOT in the small-amplitude regime ***'}`);
}

// ============================================================ U1 pulse speed
console.log('');
console.log('============================================================');
console.log('U1: do the bumps launch SOUND pulses (about 0.7 cells/step)?');
console.log('============================================================');
const ON = US.map(u => runAll(u, true, AMP));
const OFF = US.map(u => runAll(u, false, AMP));
function meanField(runs, key, t) { const p = new Float64Array(W); for (const r of runs) for (let x = 0; x < W; x++) p[x] += r[key][t][x] / runs.length; return p; }
{
  console.log('  position of the largest positive deviation in a window that follows the right-going wave:');
  console.log('     t    peak x   speed from previous row');
  let prevX = null, prevT = null;
  const speeds = [];
  for (const t of [20, 40, 60, 80, 100, 120]) {
    const p = meanField(ON, '1', t);
    let bx = 96, bv = -1e9;
    for (let x = 96; x < 200; x++) if (p[x] > bv) { bv = p[x]; bx = x; }
    const sp = (prevX === null) ? '' : ((bx - prevX) / (t - prevT)).toFixed(3);
    if (prevX !== null) speeds.push((bx - prevX) / (t - prevT));
    console.log(`    ${String(t).padStart(3)}    ${String(bx).padStart(4)}     ${sp}`);
    prevX = bx; prevT = t;
  }
  const meanSpeed = speeds.reduce((a, b) => a + b, 0) / speeds.length;
  console.log(`    mean speed over those intervals: ${meanSpeed.toFixed(3)}`);
  console.log(`    => ${Math.abs(meanSpeed - 0.707) < 0.2 ? 'SOUND branch (about 0.7), as predicted' : (meanSpeed > 0.85 ? 'BALLISTIC (about 1.0)' : 'in between: the pulse is not a clean sound mode')}`);
}
console.log('');

// ============================================================ U2/U3/U4 coupling
console.log('============================================================');
console.log('U2/U3: R(t) = d1 + d2 - d12, the mode-mode coupling');
console.log('============================================================');
function residualSeries(runs) {
  const out = [];
  for (let t = 0; t <= T; t++) {
    const R = new Float64Array(W);
    for (const r of runs) for (let x = 0; x < W; x++) R[x] += (r['1'][t][x] + r['2'][t][x] - r['12'][t][x]) / runs.length;
    let mx = 0;
    for (let x = 0; x < W; x++) if (Math.abs(R[x]) > mx) mx = Math.abs(R[x]);
    out.push(mx);
  }
  return out;
}
const Roff = residualSeries(OFF), Ron = residualSeries(ON);
{
  console.log(`  the fronts are predicted to meet at t ~ ${MEET}`);
  console.log('     t     |R| collisions OFF     |R| collisions ON    pulse height (single bump)');
  for (const t of [10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200]) {
    const ph = Math.max(...meanField(ON, '1', t).map(Math.abs));
    console.log(`    ${String(t).padStart(4)}        ${Roff[t].toExponential(2).padStart(10)}          ${Ron[t].toExponential(2).padStart(10)}            ${ph.toFixed(4)}`);
  }
  const offMax = Math.max(...Roff);
  console.log(`  U2 control: max |R| with collisions OFF = ${offMax.toExponential(3)}  ${offMax < 1e-9 ? 'PASS: identically zero, as free streaming requires' : '*** FAIL ***'}`);
  const before = Math.max(...Ron.slice(0, MEET - 5));
  const after = Math.max(...Ron.slice(MEET + 20));
  console.log(`  U3 on the meeting: max |R| before t=${MEET - 5} is ${before.toExponential(3)}; after t=${MEET + 20} it is ${after.toExponential(3)}`);
  console.log(`    => ${before < 1e-9 && after > 1e-6 ? 'the coupling is EXACTLY ZERO until the fronts meet and then switches on -- the same signature as the glider test' : (before > 1e-9 ? 'nonzero before the meeting: something is wrong' : 'no growth after the meeting')}`);
  const ph = Math.max(...meanField(ON, '1', 60).map(Math.abs));
  console.log(`  U4 coupling relative to the pulse height: ${(100 * after / ph).toFixed(2)}% after the crossing`);
}
// pass-through fidelity: outgoing pulse height with and without the partner
{
  const p1 = meanField(ON, '1', 200), p12 = meanField(ON, '12', 200);
  const h1 = Math.max(...p1.map(Math.abs));
  let h12 = 0, h12b = 0;
  for (let x = 224; x < W; x++) h12 = Math.max(h12, Math.abs(p12[x]));
  for (let x = 0; x < 96; x++) h12b = Math.max(h12b, Math.abs(p12[x]));
  console.log('');
  console.log(`  pass-through: at t=200 the pulse from bump 1 alone has height ${h1.toFixed(4)};`);
  console.log(`    with the partner present the outgoing heights are ${h12.toFixed(4)} (right side) and ${h12b.toFixed(4)} (left side)`);
}
console.log('');
console.log('============================================================');
console.log('THE DECIDING NUMBER: how does the coupling scale with amplitude?');
console.log('============================================================');
console.log('  (a weakly nonlinear medium must show the coupling falling faster than the signal)');
console.log('    amplitude   peak pulse   peak |R| after crossing   |R| / pulse');
const AMPS = [0.0125, 0.025, 0.05];
const rows = [];
for (const Amp of AMPS) {
  const runs = US.map(u => runAll(u, true, Amp));
  let pulse = 0;
  for (const t of [40, 60, 80]) pulse = Math.max(pulse, ...meanField(runs, '1', t).map(Math.abs));
  let Rmax = 0;
  for (let t = MEET + 20; t <= T; t++) {
    const R = new Float64Array(W);
    for (const r of runs) for (let x = 0; x < W; x++) R[x] += (r['1'][t][x] + r['2'][t][x] - r['12'][t][x]) / runs.length;
    for (let x = 0; x < W; x++) Rmax = Math.max(Rmax, Math.abs(R[x]));
  }
  rows.push({ Amp, pulse, Rmax, frac: Rmax / pulse });
  console.log(`    ${Amp.toFixed(4)}      ${pulse.toFixed(4)}      ${Rmax.toFixed(4)}                  ${(100 * Rmax / pulse).toFixed(1)}%`);
}
{
  // fit the exponent: R ~ Amp^p
  const a = rows[0], b = rows[rows.length - 1];
  const p = Math.log(b.Rmax / a.Rmax) / Math.log(b.Amp / a.Amp);
  console.log(`  fitted exponent p in  |R| ~ amplitude^p :  p = ${p.toFixed(2)}`);
  console.log(`    (p = 1 means the coupling is as strong as the signal and the medium is NOT weakly`);
  console.log(`     nonlinear at any amplitude -- a messenger that always interacts with itself;`);
  console.log(`     p = 2 means a second-order hydrodynamic nonlinearity that dies away at small`);
  console.log(`     amplitude -- i.e. a medium that can carry an approximately free photon.)`);
}
console.log('');
console.log('============================================================');
console.log('WHAT THIS SETTLES');
console.log('============================================================');
console.log('  R(t) is exactly zero while the two sound pulses are apart -- not small, ZERO -- and it');
console.log('  is exactly zero at ALL times when collisions are switched off. So whatever appears in');
console.log('  R after the meeting is genuine mode-mode coupling, with no noise floor to interpret.');
console.log('  That is the number that decides whether this collective mode can serve as a messenger:');
console.log('  zero would make it an ideal photon, comparable to the pulse height would make it a');
console.log('  soliton, and the truth is the value printed above.');
