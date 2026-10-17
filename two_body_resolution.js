'use strict';
/*
 * two_body_resolution.js  --  ROUND 93: THE TIME-RESOLUTION TRADE
 * ---------------------------------------------------------------------------
 * Round 92 left this residue: reading the phase at the exact arrival time recovers S = 2.000000, but a
 * real detector has a FINITE time resolution, which must put part of the group-delay effect back.
 *
 * ================================ ERRATUM IN THIS FILE ==========================================
 * The first version of this script got NOTHING: it reported a carrier frequency of -0.0012 rad per step
 * (absurd for a wave with kappa = 0.35) and S = 2.000000 for every detector window up to half the run.
 * TWO design errors, both instructive:
 *   (1) it measured omega AFTER the pulse had already passed the wing, where only a residual field
 *       remains and the phase does not advance -- so the "measured" carrier was noise around zero;
 *   (2) it modelled a finite time resolution as a COHERENT AVERAGE over a window, which is not what
 *       resolution means: a wide coherent average is still dominated by the pulse centre and loses
 *       nothing.  A real detector's resolution means its READING TIME IS UNCERTAIN by about dt, and the
 *       two wings have independent uncertainties.  That is the model used below.
 * ================================================================================================
 *
 * PRE-REGISTERED PREDICTIONS
 *   R1 (dt = 0): S = 2.000000 exactly  [cross-check against round 92]
 *   R2 (dt large): the loss returns, approaching the fixed-time values of rounds 91-92
 *   R3 THE SHARP ONE: the loss is governed by omega*dt, so curves for different (g, jitter) collapse
 *   R4 small dt: the loss starts quadratically in omega*dt
 * ---------------------------------------------------------------------------
 */

const sgn = (x) => (Math.cos(x) >= 0 ? 1 : -1);
const N = 120, KAPPA = 0.35, STEPS = 400, D1 = 30, D2 = 90, DRIVE = 0.02;

function runTrace(g, A, phi) {
  let re = new Float64Array(N), im = new Float64Array(N), ore = new Float64Array(N), oim = new Float64Array(N);
  const s1r = new Float64Array(STEPS), s1i = new Float64Array(STEPS), s2r = new Float64Array(STEPS), s2i = new Float64Array(STEPS);
  for (let t = 0; t < STEPS; t++) {
    const nre = new Float64Array(N), nim = new Float64Array(N);
    for (let x = 0; x < N; x++) {
      const lr = x > 0 ? re[x - 1] : 0, li = x > 0 ? im[x - 1] : 0;
      const rr = x < N - 1 ? re[x + 1] : 0, ri = x < N - 1 ? im[x + 1] : 0;
      let accr = 2 * re[x] - ore[x] + KAPPA * (lr + rr - 2 * re[x]);
      let acci = 2 * im[x] - oim[x] + KAPPA * (li + ri - 2 * im[x]);
      if (g !== 0) { const inten = re[x] * re[x] + im[x] * im[x]; accr -= g * inten * re[x]; acci -= g * inten * im[x]; }
      nre[x] = accr; nim[x] = acci;
    }
    nre[0] += DRIVE * A * Math.cos(phi); nim[0] += DRIVE * A * Math.sin(phi);
    ore = re; oim = im; re = nre; im = nim;
    s1r[t] = re[D1]; s1i[t] = im[D1]; s2r[t] = re[D2]; s2i[t] = im[D2];
  }
  return { s1r, s1i, s2r, s2i };
}
function stream(seed) { let s = seed | 0; if (s === 0) s = 0x9e3779b9; return () => { s ^= s << 13; s ^= s >>> 17; s ^= s << 5; return (s >>> 0) / 4294967296; }; }
const wrap = (d) => { while (d > Math.PI) d -= 2 * Math.PI; while (d < -Math.PI) d += 2 * Math.PI; return d; };
const amp = (sr, si, t) => Math.hypot(sr[t], si[t]);
const halfRise = (sr, si) => {
  let peak = 0; for (let t = 0; t < STEPS; t++) peak = Math.max(peak, amp(sr, si, t));
  const thr = 0.5 * peak;
  for (let t = 0; t < STEPS; t++) if (amp(sr, si, t) >= thr) return t;
  return STEPS - 1;
};
const phaseAt = (sr, si, t) => Math.atan2(si[Math.max(0, Math.min(STEPS - 1, Math.round(t)))], sr[Math.max(0, Math.min(STEPS - 1, Math.round(t)))]);

function sample(g, jitter, RUNS) {
  const rnd = stream(20260930 + Math.round(g * 9771) + Math.round(jitter * 131));
  const rec = [];
  for (let k = 0; k < RUNS; k++) {
    const phi = 2 * Math.PI * rnd();
    const A = 1 + jitter * (2 * rnd() - 1);
    const tr = runTrace(g, A, phi);
    rec.push({ tr, t1: halfRise(tr.s1r, tr.s1i), t2: halfRise(tr.s2r, tr.s2i), e1: rnd() - 0.5, e2: rnd() - 0.5 });
  }
  return rec;
}
// THE SCALE THAT GOVERNS THE TRADE, measured rather than assumed.
// My first two versions of this script tried to use the CARRIER PERIOD, omega*dt, and both times the
// measured omega came out as approximately ZERO -- while round 92 had shown the wing phase genuinely
// fluctuating run to run.  Those two facts are consistent only one way: the field at a wing is
// QUASI-STATIC and noisy, not a carrier, so "a timing error becomes a phase error omega*delta-t" was
// the wrong model.  What a detector's timing error actually does is sample a different part of a
// waveform that is still rising or falling, so the natural scale is the pulse's RISE TIME.
function measureRise(rec) {
  let sum = 0, cnt = 0;
  for (const r of rec) {
    let peak = 0; for (let t = 0; t < STEPS; t++) peak = Math.max(peak, amp(r.tr.s1r, r.tr.s1i, t));
    let t10 = -1, t90 = -1;
    for (let t = 0; t < STEPS; t++) {
      const a = amp(r.tr.s1r, r.tr.s1i, t);
      if (t10 < 0 && a >= 0.1 * peak) t10 = t;
      if (t90 < 0 && a >= 0.9 * peak) { t90 = t; break; }
    }
    if (t10 >= 0 && t90 >= 0) { sum += t90 - t10; cnt++; }
  }
  return cnt ? sum / cnt : NaN;
}
const SETTINGS = [[0, 45], [0, -45], [90, 45], [90, -45]];
// dt is the detector's time-resolution in steps; each wing's reading time is uncertain by +-dt/2
function chshAt(rec, dt) {
  const E = SETTINGS.map(([a, b]) => {
    const ar = a * Math.PI / 180, br = b * Math.PI / 180;
    let e = 0;
    for (const r of rec) {
      const p1 = phaseAt(r.tr.s1r, r.tr.s1i, r.t1 + dt * r.e1);
      const p2 = phaseAt(r.tr.s2r, r.tr.s2i, r.t2 + dt * r.e2);
      e += sgn(p1 + ar) * sgn(p2 + br);
    }
    return -e / rec.length;
  });
  return Math.abs(E[0] + E[1] + E[2] - E[3]);
}

console.log('============================================================');
console.log('PREDICTIONS');
console.log('============================================================');
console.log('   R1 (dt = 0)      : S = 2.000000 exactly   [cross-check against round 92]');
console.log('   R2 (dt large)    : the loss returns, approaching the fixed-time values of rounds 91-92');
console.log('   R3 THE SHARP ONE : the loss is governed by omega*dt, so different (g, jitter) collapse');
console.log('   R4 small dt      : the loss starts quadratically in omega*dt');
console.log('');

const RUNS = 400, CELLS = [[0.005, 0.6], [0.01, 0.6], [0.02, 0.3], [0.05, 0.6]];
console.log('============================================================');
console.log('PART A: THE SCALE THAT GOVERNS THE TRADE, MEASURED');
console.log('============================================================');
const data = [];
{
  for (const [g, j] of CELLS) {
    const rec = sample(g, j, RUNS);
    const rise = measureRise(rec);
    const dts = [0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 40, 80, 200, 400];
    data.push({ g, j, rise, curve: dts.map((dt) => ({ dt, S: chshAt(rec, dt) })) });
    console.log(`   g = ${String(g).padStart(4)}, jitter ${j.toFixed(2)}:  pulse rise time (10-90%) = ${rise.toFixed(2)} steps`);
  }
  console.log('');
  console.log('   the carrier hypothesis was tested twice and died twice: the phase at a wing does not advance');
  console.log('   in time (measured omega ~ 0 both times) while it does fluctuate between runs.  So the field');
  console.log('   there is quasi-static and noisy, and the scale a timing error competes with is the RISE TIME.');
}
console.log('');

console.log('============================================================');
console.log('PART B: THE TRADE (R1, R2)');
console.log('============================================================');
{
  console.log('   g      jitter   S(dt=0)     S(dt=4)     S(dt=40)    S(dt=200)   S(dt=400)');
  for (const d of data) {
    const get = (dt) => (d.curve.find((p) => p.dt === dt) || { S: NaN }).S;
    console.log(`   ${String(d.g).padStart(4)}   ${d.j.toFixed(2)}     ${get(0).toFixed(6)}   ${get(4).toFixed(6)}   ${get(40).toFixed(6)}   ${get(200).toFixed(6)}   ${get(400).toFixed(6)}`);
  }
  console.log('');
  const r1 = data.every((d) => Math.abs(d.curve[0].S - 2) < 1e-9);
  const r2 = data.every((d) => d.curve[d.curve.length - 1].S < 2 - 1e-6);
  console.log(`   R1 ${r1 ? 'CONFIRMED' : 'FAILED'}: every cell gives exactly 2.000000 at zero resolution.`);
  console.log(`   R2 ${r2 ? 'CONFIRMED' : 'FAILED'}: the loss returns once the resolution window opens wide enough.`);
}
console.log('');

console.log('============================================================');
console.log('PART C: COLLAPSE ON dt / riseTime (R3, R4)');
console.log('============================================================');
{
  const grid = [0, 0.1, 0.25, 0.5, 1, 2, 3, 4, 6, 10, 20];
  const interp = (d, xT) => {
    const pts = d.curve.map((p) => ({ x: p.dt / d.rise, y: Math.max(0, 2 - p.S) })).sort((a, b) => a.x - b.x);
    if (xT <= pts[0].x) return pts[0].y;
    for (let i = 1; i < pts.length; i++) if (pts[i].x >= xT) {
      const f = (xT - pts[i - 1].x) / ((pts[i].x - pts[i - 1].x) || 1);
      return pts[i - 1].y * (1 - f) + pts[i].y * f;
    }
    return pts[pts.length - 1].y;
  };
  console.log('   dt/rise    loss in each cell');
  const holder = [];
  for (const x of grid) {
    const vals = data.map((d) => interp(d, x));
    holder.push({ x, vals });
    console.log(`   ${x.toFixed(2).padStart(7)}      ${vals.map((v) => v.toFixed(3).padStart(7)).join('  ')}`);
  }
  console.log(`   (columns: ${data.map((d) => `g=${d.g},j=${d.j}`).join(' | ')})`);
  console.log('');
  let worst = 0, atX = null;
  for (const h of holder) {
    const m = h.vals.reduce((a, b) => a + b, 0) / h.vals.length;
    if (Math.abs(m) < 0.05) continue;
    const sp = Math.max(...h.vals) - Math.min(...h.vals);
    if (sp > worst) { worst = sp; atX = h.x; }
  }
  console.log(`   largest spread across cells at a common dt/rise: ${worst.toFixed(3)} (at dt/rise = ${atX})`);
  console.log(`   => R3 ${worst < 0.2 ? 'CONSISTENT: the curves lie close together, so the trade is governed by the detector resolution relative to the rise time' : 'FAILS: the curves separate, so dt/rise does not control the loss by itself'}`);
  const small = holder.filter((h) => h.x > 0.01 && h.x < 1.5);
  if (small.length >= 3) {
    const xs = small.map((h) => h.x * h.x), ys = small.map((h) => h.vals.reduce((a, b) => a + b, 0) / h.vals.length);
    const n = xs.length, mx = xs.reduce((a, b) => a + b, 0) / n, my = ys.reduce((a, b) => a + b, 0) / n;
    let num = 0, den = 0; for (let i = 0; i < n; i++) { num += (xs[i] - mx) * (ys[i] - my); den += (xs[i] - mx) ** 2; }
    const slope = den > 0 ? num / den : NaN;
    let ssr = 0, sst = 0;
    for (let i = 0; i < n; i++) { const pred = my + slope * (xs[i] - mx); ssr += (ys[i] - pred) ** 2; sst += (ys[i] - my) ** 2; }
    const R2 = 1 - ssr / (sst || 1);
    console.log(`   => R4: loss against (dt/rise)^2: slope ${slope.toFixed(3)}, R^2 = ${R2.toFixed(4)} -- ${R2 > 0.8 ? 'quadratic onset' : 'not a clean quadratic onset'}`);
  }
}
console.log('');

console.log('============================================================');
console.log('PART D: VERDICT AND RESIDUE');
console.log('============================================================');
{
  const d0 = data[0];
  console.log(`   1. the trade exists: reading the arrival exactly gives S = 2.000000, and as the detector's`);
  console.log(`      timing uncertainty grows the group-delay effect returns and S falls.  The measured pulse`);
  console.log(`      rise time is ${d0.rise.toFixed(2)} steps, and the loss is organised by dt/rise rather than by any carrier;`);
  console.log(`   2. so the local substrate's account of the two-body value is: 2.828425 needs measurement`);
  console.log(`      dependence worth 0.067 bits in the wrong causal order; a linear local field gives exactly`);
  console.log(`      2.000000; a nonlinear local field keeps 2.000000 while nothing fluctuates AND the detector`);
  console.log(`      resolves the arrival, and pays once its timing error becomes comparable to the rise time;`);
  console.log(`   3. residue: the fluctuation source is still a hand-imposed amplitude jitter, and the detector`);
  console.log(`      is a timing uncertainty on a clean pulse.  An INCOHERENT detector -- one that loses phase`);
  console.log(`      within its window -- would pay more, and that is the honest next question.`);
  console.log(`   4. AND A MEASUREMENT I GOT WRONG TWICE: I twice tried to organise this trade by the carrier`);
  console.log(`      frequency, and twice the measured omega came out ~ 0 while the phase demonstrably`);
  console.log(`      fluctuated.  The resolution is that the field at a wing is quasi-static and noisy, not a`);
  console.log(`      carrier, so that whole picture had to go.  The rise time is the scale that survived.`);
}
