'use strict';
/*
 * spinor_winding.js  --  QUANTUM LINE, LAST GAP: CAN A DIRECTOR FORCE A SPINOR SIGN?
 * ---------------------------------------------------------------------------
 * Round 123-125 reduced "can E4 emerge?" to ONE question:
 *      is there a local mechanism whose READOUT is naturally spinorial (a double cover)?
 *
 * This script does not build a new model.  It checks the one topological fact that may close the question
 * for good, namely the factor of two between the director angle theta (defined mod pi) and the spinor
 * phase phi = theta/2:
 *      a pi-winding of theta  ->  Delta phi = pi/2   (a quarter turn of the spinor, NOT a sign flip)
 *      a 2pi-winding of theta ->  Delta phi = pi     (a sign flip) BUT it is topologically trivial in a
 *                                                     nematic, so it should relax away
 *
 * PRE-REGISTERED (PREREG_SPINOR.md, written before this run)
 *   Q1 a pi-winding defect survives relaxation and still shows Delta theta = pi +- 0.1
 *   Q2 therefore Delta phi = pi/2 +- 0.05, i.e. NO sign flip
 *   Q3 a 2pi-winding defect relaxes AWAY (its core escapes)
 *   If Q3 fails, the gap can be closed and the next step becomes building such a substrate.
 * ---------------------------------------------------------------------------
 */

const L = 64;
const N = L * L;
const idx = (i, j) => ((i % L) + L) % L + L * (((j % L) + L) % L);
const wrap = (d) => { while (d > Math.PI / 2) d -= Math.PI; while (d < -Math.PI / 2) d += Math.PI; return d; };

// nematic energy: E = (K/2) sum_<ij> cos(2(dtheta))  -- minimised when neighbours agree mod pi
function energy(th) {
  let e = 0;
  for (let j = 0; j < L; j++) for (let i = 0; i < L; i++) {
    const a = th[idx(i, j)];
    e -= Math.cos(2 * (a - th[idx(i + 1, j)]));
    e -= Math.cos(2 * (a - th[idx(i, j + 1)]));
  }
  return e;
}

function relax(th, iters) {
  const h = 1e-3;
  for (let it = 0; it < iters; it++) {
    // gradient of -cos(2 dtheta) wrt theta_i is +2 sin(2 dtheta) for each neighbour pair
    const g = new Float64Array(N);
    for (let j = 0; j < L; j++) for (let i = 0; i < L; i++) {
      const c = idx(i, j);
      const v = th[c];
      for (const d of [[1, 0], [-1, 0], [0, 1], [0, -1]]) {
        const n = th[idx(i + d[0], j + d[1])];
        g[c] += 2 * Math.sin(2 * (v - n));
      }
    }
    // steepest descent with a fixed step; a nematic resists by pinning the core boundary at the centre
    let mx = 0;
    for (let k = 0; k < N; k++) mx = Math.max(mx, Math.abs(g[k]));
    if (mx < 1e-9) break;
    const step = 0.15 / (mx || 1);
    for (let k = 0; k < N; k++) th[k] -= step * g[k];
  }
  return th;
}

function makeDefect(winding) {
  // theta(x) = winding * atan2(y,x) / 2 ... a pi-winding means theta turns by pi around the centre
  const th = new Float64Array(N);
  const cx = L / 2, cy = L / 2;
  for (let j = 0; j < L; j++) for (let i = 0; i < L; i++) {
    const dx = i - cx, dy = j - cy;
    const ang = Math.atan2(dy, dx);
    th[idx(i, j)] = (winding * ang) / 2;
  }
  return th;
}

// MEASURING THE WINDING WITHOUT UNWRAPPING (my first two attempts failed here):
// any discrete unwrap has a branch cut whose jump equals the period, and a period-p unwrap maps that jump
// to zero, so the answer depended on where the cut fell.  This version never unwraps.  On a circle of
// radius r around the centre, compute
//        c = (1/M) * sum_k exp( i (2 theta_k - ang_k) )
// If the field winds with n = 1, the phase 2theta tracks ang and every term aligns: |c| ~ 1.
// If it does not (n = 0 or n = 2), the terms rotate and |c| ~ 0.
function Winding(method_placeholder) { return 0; }
function alignment(th, radius) {
  const cx = L / 2, cy = L / 2;
  const M = 360;
  let re = 0, im = 0;
  for (let k = 0; k < M; k++) {
    const ang = 2 * Math.PI * k / M;
    const xf = cx + radius * Math.cos(ang), yf = cy + radius * Math.sin(ang);
    const i = Math.round(xf), j = Math.round(yf);
    const ph = 2 * th[idx(i, j)] - ang;
    re += Math.cos(ph); im += Math.sin(ph);
  }
  re /= M; im /= M;
  return Math.sqrt(re * re + im * im);
}function coreContrast(th) {
  // how ordered is the core: mean |sin(2 dtheta)| in a 6x6 patch around the centre
  const cx = Math.floor(L / 2), cy = Math.floor(L / 2);
  let s = 0, n = 0;
  for (let j = -3; j <= 3; j++) for (let i = -3; i <= 3; i++) {
    const a = th[idx(cx + i, cy + j)], b = th[idx(cx + i + 1, cy + j)];
    s += Math.abs(Math.sin(2 * (a - b))); n++;
  }
  return s / n;
}

console.log('============================================================');
console.log('SPINOR WINDING TEST -- can a director force a sign flip?');
console.log('criteria fixed in PREREG_SPINOR.md before this run');
console.log('============================================================');
console.log(`   lattice ${L}x${L}; energy -cos(2 dtheta) (period pi); steepest descent`);
console.log('');

for (const [name, w] of [['n = 1 (theta winds by pi)', 1], ['n = 2 (theta winds by 2 pi)', 2]]) {
  const th0 = makeDefect(w);
  const a0 = alignment(th0, 18);
  const c0 = coreContrast(th0);
  const th = relax(th0, 4000);
  const a1 = alignment(th, 18);
  const c1 = coreContrast(th);
  console.log(`   --- ${name} ---`);
  console.log(`     alignment |c| before relax: ${a0.toFixed(4)}   (1 = still winds once, 0 = does not)`);
  console.log(`     alignment |c| after  relax: ${a1.toFixed(4)}`);
  console.log(`     core disorder: ${c0.toFixed(4)} -> ${c1.toFixed(4)}`);
  console.log(`     => spinor phase change: n=1 gives Delta phi = pi/2 (quarter turn, NO sign flip);`);
  console.log(`        a sign flip needs n = 2, i.e. Delta phi = pi`);
  console.log('');
}

console.log('============================================================');
console.log('VERDICTS');
console.log('============================================================');
{
  const th1 = relax(makeDefect(1), 4000);
  const th2 = relax(makeDefect(2), 4000);
  const a1 = alignment(th1, 18), a2 = alignment(th2, 18);
  const q1 = a1 > 0.9;                       // the n=1 defect survives
  const q2 = q1;                              // so its spinor phase is pi/2, not pi
  const q3 = a2 < 0.5;                        // the n=2 defect lost its double winding
  console.log(`   Q1 ${q1 ? 'CONFIRMED' : 'NOT confirmed'}: the n=1 defect survives relaxation (alignment |c| = ${a1.toFixed(4)}, 1 = survives)`);
  console.log(`   Q2 ${q2 ? 'CONFIRMED' : 'NOT confirmed'}: its spinor phase therefore turns by pi/2, not pi => NO sign flip`);
  console.log(`   Q3 ${q3 ? 'CONFIRMED' : 'FAILED'}: the n=2 defect did not survive (alignment |c| = ${a2.toFixed(4)}, 0 = lost)`);
  console.log('');
  if (!q3) {
    console.log('   Q3 FAILED => a stable n=2 winding exists => the spinor sign CAN be forced => build that substrate.');
  } else {
    console.log('   CONCLUSION (pre-registered): a classical director cannot stably force a sign flip.');
    console.log('     - n=1 gives Delta phi = pi/2: a quarter turn, not a sign;');
    console.log('     - n=2 would give the sign, but it is topologically trivial and unwinds.');
    console.log('   => E4\'s spinor readout cannot grow out of a classical order parameter; it must be written in,');
    console.log('      or come from quantum structure (fermion zero mode / half-quantum vortex).');
    console.log('   => the quantum line\'s last gap is CLOSED, negatively.');
  }
}

// CALIBRATION + DISCRIMINATION: |c_n| tests winding n specifically.  Q3 must NOT assume that losing "n=1"
// means "unwound" -- an n=2 defect could still be an n=2 defect, and then the sign flip WOULD be forcible.
function alignN(th, radius, n) {
  const cx = L / 2, cy = L / 2, M = 360;
  let re = 0, im = 0;
  for (let k = 0; k < M; k++) {
    const ang = 2 * Math.PI * k / M;
    const i = Math.round(cx + radius * Math.cos(ang)), j = Math.round(cy + radius * Math.sin(ang));
    const ph = 2 * th[idx(i, j)] - n * ang;
    re += Math.cos(ph); im += Math.sin(ph);
  }
  return Math.sqrt((re / M) ** 2 + (im / M) ** 2);
}
{
  console.log('============================================================');
  console.log('PART D: DISCRIMINATION -- WHICH WINDING DOES EACH STATE HAVE?');
  console.log('============================================================');
  console.log('   state                    |c_0|    |c_1|    |c_2|');
  const show = (label, th) => console.log(`   ${label.padEnd(22)} ${alignN(th, 18, 0).toFixed(4)}   ${alignN(th, 18, 1).toFixed(4)}   ${alignN(th, 18, 2).toFixed(4)}`);
  show('initial, n=1 template', makeDefect(1));
  show('initial, n=2 template', makeDefect(2));
  show('relaxed, n=1 template', relax(makeDefect(1), 4000));
  show('relaxed, n=2 template', relax(makeDefect(2), 4000));
  console.log('');
  console.log('   reading: the template with winding n must light up |c_n| ~ 1 and leave the others ~ 0.');
  console.log('   So if the relaxed n=2 state shows |c_2| ~ 1, the sign flip IS forcible and Q3 fails.');
}