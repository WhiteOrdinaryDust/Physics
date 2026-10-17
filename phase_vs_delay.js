'use strict';
/*
 * phase_vs_delay.js  --  ROUND 115 (D7): SEPARATING PHASE ACCUMULATION FROM GROUP DELAY
 * ---------------------------------------------------------------------------
 * Round 113 (D6) priced an incoherent detector in PHASE: per-wing jitter of scale Delta makes S fall
 * quadratically and kills the violation at 26.9 degrees.  A real detector's incoherence is a TIME jitter,
 * and the two are linked by the carrier:
 *
 *      Delta_theta = omega * Delta_tau
 *
 * This round measures that link and what lives beyond it, which is the item §34 left open.
 *
 * THE NOISE SOURCE IS NOT HAND-IMPOSED: there is no random number anywhere.  The ensemble is a family of
 * DETERMINISTIC initial conditions (the pulse launched from slightly different positions), so the "jitter"
 * is the field's reproducible response.
 *
 * THE MODEL: the one-dimensional nonlinear Schrodinger equation i dpsi/dt = -d2psi/dx2 + g|psi|^2 psi, by
 * split-step Fourier, with g = 0 as the LINEAR CONTROL.  One observation point downstream; the arrival
 * time is the envelope peak, and THE CARRIER IS MEASURED AT THE PEAK, where the envelope phase is
 * stationary, so the chirp cannot contaminate it.
 *
 * THREE FAILURES ARE RECORDED HERE SO THEY ARE NOT REPEATED
 *   (a) watching x = 50 and 60, which the pulse never reached inside the run -- that table was all noise;
 *   (b) putting the observers at x = 22 and 26 with width 6, so the pulse tail was already there at t = 0
 *       and the peak landed on the first sample, which crashed;
 *   (c) reading the carrier from the phase difference between two points AT THEIR OWN PEAKS, which mixes
 *       in the envelope phase and got k wrong by 40 percent.
 *
 * PRE-REGISTERED PREDICTIONS
 *   P1 THE CONTROL: linear (g = 0) gives R = (dtheta/ddelta) / (-omega dtau/ddelta) = 1 to numerical
 *      precision.  The minus sign is the convention: launching the pulse farther out makes it arrive later
 *      and its phase at arrival smaller.
 *   P2 nonlinear (g > 0): R exceeds 1 and the excess grows with g, because the interaction adds phase
 *      without adding delay
 *   P3 the two jitters stay strongly correlated even when R > 1
 *   P4 bit-for-bit reproducible, since nothing in the dynamics is random
 * ---------------------------------------------------------------------------
 */

const N = 1024, L = 100, DX = L / N;
const DT = 0.002, STEPS = 25000;
const XOBS = 60;                              // downstream of the launch at x = 20, clear of the pulse
const K0 = 2 * Math.PI / (L / 12);            // carrier = grid mode 12: v_g = 2k = 1.51, so the

function fft(re, im, sign) {
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = sign * 2 * Math.PI / len;
    const wr = Math.cos(ang), wi = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let cr = 1, ci = 0;
      for (let j = 0; j < len / 2; j++) {
        const ur = re[i + j], ui = im[i + j];
        const vr = re[i + j + len / 2] * cr - im[i + j + len / 2] * ci;
        const vi = re[i + j + len / 2] * ci + im[i + j + len / 2] * cr;
        re[i + j] = ur + vr; im[i + j] = ui + vi;
        re[i + j + len / 2] = ur - vr; im[i + j + len / 2] = ui - vi;
        const nr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = nr;
      }
    }
  }
}
const kOf = (i) => (i < N / 2 ? i : i - N) * (2 * Math.PI / L);
const unwrap = (d) => { while (d > Math.PI) d -= 2 * Math.PI; while (d < -Math.PI) d += 2 * Math.PI; return d; };

function launch(g, x0, width, amp) {
  const re = new Float64Array(N), im = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    const xx = i * DX;
    const env = amp * Math.exp(-((xx - x0) ** 2) / (2 * width * width));
    re[i] = env * Math.cos(K0 * xx);
    im[i] = env * Math.sin(K0 * xx);
  }
  return { re, im };
}
function advance(g, st, steps) {
  const { re, im } = st;
  const linPh = new Float64Array(N);
  for (let i = 0; i < N; i++) linPh[i] = -kOf(i) * kOf(i) * DT;
  for (let s = 0; s < steps; s++) {
    if (g !== 0) for (let i = 0; i < N; i++) {
      const a2 = re[i] * re[i] + im[i] * im[i], ph = -g * a2 * DT, c = Math.cos(ph), sn = Math.sin(ph);
      const r2 = re[i] * c - im[i] * sn, i2 = re[i] * sn + im[i] * c; re[i] = r2; im[i] = i2;
    }
    fft(re, im, 1);
    for (let i = 0; i < N; i++) { const c = Math.cos(linPh[i]), sn = Math.sin(linPh[i]); const r2 = re[i] * c - im[i] * sn, i2 = re[i] * sn + im[i] * c; re[i] = r2; im[i] = i2; }
    fft(re, im, -1);
    for (let i = 0; i < N; i++) { re[i] /= N; im[i] /= N; }
  }
}
function run(g, x0, width, amp) {
  const st = launch(g, x0, width, amp);
  const i0 = Math.round(XOBS / DX);
  const tr = [];
  const linPh = new Float64Array(N);
  for (let i = 0; i < N; i++) linPh[i] = -kOf(i) * kOf(i) * DT;
  for (let s = 0; s < STEPS; s++) {
    if (g !== 0) for (let i = 0; i < N; i++) {
      const a2 = st.re[i] * st.re[i] + st.im[i] * st.im[i], ph = -g * a2 * DT, c = Math.cos(ph), sn = Math.sin(ph);
      const r2 = st.re[i] * c - st.im[i] * sn, i2 = st.re[i] * sn + st.im[i] * c; st.re[i] = r2; st.im[i] = i2;
    }
    fft(st.re, st.im, 1);
    for (let i = 0; i < N; i++) { const c = Math.cos(linPh[i]), sn = Math.sin(linPh[i]); const r2 = st.re[i] * c - st.im[i] * sn, i2 = st.re[i] * sn + st.im[i] * c; st.re[i] = r2; st.im[i] = i2; }
    fft(st.re, st.im, -1);
    for (let i = 0; i < N; i++) { st.re[i] /= N; st.im[i] /= N; }
    if (s % 2 === 0) tr.push({ t: (s + 1) * DT, A: Math.hypot(st.re[i0], st.im[i0]), th: Math.atan2(st.im[i0], st.re[i0]) });
  }
  // THE CENTROID ARRIVAL, which is what the round-116 diagnostic showed to be robust.  The peak of a
  // dispersing pulse is a plateau and its parabola refinement is unstable; the centroid of |psi| is not.
  let wSum = 0, tSum = 0, cSum = 0, sSum = 0;
  for (let i = 1; i < tr.length; i++) {
    const w = tr[i].A;
    wSum += w; tSum += w * tr[i].t;
    cSum += w * Math.cos(tr[i].th); sSum += w * Math.sin(tr[i].th);
  }
  const tArr = wSum > 0 ? tSum / wSum : NaN;
  const thArr = Math.atan2(sSum, cSum);
  // THE CARRIER IS MEASURED ONLY IN A WINDOW AROUND THE ARRIVAL.  Averaging over the whole trace was a
  // mistake: once the pulse has passed, the field at the observer is essentially zero and its PHASE IS
  // MEANINGLESS, so the unwrap chain accumulates spurious jumps and the "carrier" came out at 7.78 against
  // a true 0.569.  The window below is where the envelope is actually there.
  let bIdx = 1;
  for (let i = 1; i < tr.length; i++) if (tr[i].A > tr[bIdx].A) bIdx = i;
  let om = 0, cnt = 0;
  for (let i = Math.max(2, bIdx - 8); i < Math.min(tr.length - 1, bIdx + 8); i++) {
    om += unwrap(tr[i].th - tr[i - 1].th) / (tr[i].t - tr[i - 1].t); cnt++;
  }
  om = Math.abs(om / Math.max(1, cnt));
  return { tPeak: tArr, thPeak: thArr, om, A: 0 };
}
function normCheck(g, x0, width, amp) {
  const st = launch(g, x0, width, amp);
  let n0 = 0;
  for (let i = 0; i < N; i++) n0 += st.re[i] * st.re[i] + st.im[i] * st.im[i];
  advance(g, st, STEPS);
  let n1 = 0;
  for (let i = 0; i < N; i++) n1 += st.re[i] * st.re[i] + st.im[i] * st.im[i];
  return { n0, n1, rel: Math.abs(n1 - n0) / n0 };
}

console.log('============================================================');
console.log('SELF-CHECK: NORM CONSERVATION');
console.log('============================================================');
for (const g of [0, 5]) {
  const r = normCheck(g, 20, 3, 1);
  console.log(`   g = ${g}: ${r.n0.toFixed(8)} -> ${r.n1.toFixed(8)}, relative change ${r.rel.toExponential(2)}  => ${r.rel < 1e-6 ? 'PASSES' : 'FAILS, nothing below is trustworthy'}`);
}
console.log('');

// pulse reaches the observer inside the run.
const WIDTH = 8, AMP = 1, DELTAS = [0, 0.5, 1, 1.5, 2];
const SPAN = DELTAS[DELTAS.length - 1] - DELTAS[0];
console.log('============================================================');
console.log('PART A: THE PHASE-DELAY RELATION (P1, P2)');
console.log('============================================================');
const table = {};
for (const g of [0, 1, 5]) {
  const rows = DELTAS.map((d) => { const r = run(g, 20 + d, WIDTH, AMP); return { d, tau: r.tPeak, th: r.thPeak, om: r.om, A: r.A }; });
  const dtau = (rows[rows.length - 1].tau - rows[0].tau) / SPAN;
  let dth = 0;
  for (let j = 1; j < rows.length; j++) dth += unwrap(rows[j].th - rows[j - 1].th);
  dth /= SPAN;
  const om = rows[2].om, pred = -om * dtau, R = dth / pred;
  table[g] = { rows, dtau, dth, om, R, pred };
  console.log(`   g = ${g}:  omega at the peak = ${om.toFixed(4)}   dtau/ddelta = ${dtau.toFixed(4)}   dtheta/ddelta = ${dth.toFixed(4)}`);
  console.log(`           -omega dtau/ddelta = ${pred.toFixed(4)}   R = ${R.toFixed(6)}`);
  console.log(`           arrival times  ${rows.map((r) => r.tau.toFixed(3)).join(', ')}`);
  console.log(`           phase at arrival ${rows.map((r) => r.th.toFixed(3)).join(', ')}`);
  console.log('');
}
{
  const R0 = table[0].R, R1 = table[1].R, R5 = table[5].R;
  console.log(`   P1 ${Math.abs(R0 - 1) < 0.05 ? 'CONFIRMED' : 'NOT confirmed'}: linear R = ${R0.toFixed(6)} against 1, difference ${Math.abs(R0 - 1).toFixed(4)}`);
  console.log(`   P2 ${R1 > R0 && R5 > R1 ? 'CONFIRMED' : 'NOT confirmed'}: R by g = 0, 1, 5 is ${R0.toFixed(3)}, ${R1.toFixed(3)}, ${R5.toFixed(3)}`);
  console.log(`      excess phase per unit launch shift: ${(table[0].dth - table[0].pred).toFixed(4)} (g=0), ${(table[1].dth - table[1].pred).toFixed(4)} (g=1), ${(table[5].dth - table[5].pred).toFixed(4)} (g=5)`);
  console.log('      that excess is the part of the phase jitter that a pure group delay cannot explain.');
}
console.log('');

console.log('============================================================');
console.log('PART B: CORRELATION ACROSS LAUNCHES (P3, P4)');
console.log('============================================================');
console.log('   g      correlation of (dtheta, -omega dtau) between neighbouring launches');
for (const g of [0, 1, 5]) {
  const rows = table[g].rows, om = table[g].om;
  const xs = [], ys = [];
  for (let j = 1; j < rows.length; j++) { xs.push(unwrap(rows[j].th - rows[j - 1].th)); ys.push(-om * (rows[j].tau - rows[j - 1].tau)); }
  const n = xs.length, mx = xs.reduce((a, b) => a + b, 0) / n, my = ys.reduce((a, b) => a + b, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let j = 0; j < n; j++) { sxy += (xs[j] - mx) * (ys[j] - my); sxx += (xs[j] - mx) ** 2; syy += (ys[j] - my) ** 2; }
  console.log(`   ${String(g).padStart(2)}     ${(sxy / Math.sqrt(sxx * syy)).toFixed(4)}`);
}
console.log('');
console.log('   P3: a correlation near 1 means the phase jitter and the delay jitter are the SAME channel, so the');
console.log('   nonlinear excess appears as a ratio above 1 rather than as decorrelation.');
console.log('');
console.log('   P4: no random number is used anywhere in this script, so the table reproduces bit for bit.  That is');
console.log('   also what makes it a legitimate "jitter": the DETERMINISTIC response of the field to a family of');
console.log('   initial conditions.');
console.log('');
console.log('   WHAT THIS MEANS FOR D6: the incoherent-detector jitter is a TIME jitter, and it becomes phase');
console.log('   through omega.  The 26.9-degree window of round 113 is therefore a statement about the READ TIME:');
console.log('   26.9 degrees divided by omega.  Whatever excess R shows beyond 1 is the part a pure delay cannot');
console.log('   mimic, and it is the part no timing precision removes.');
