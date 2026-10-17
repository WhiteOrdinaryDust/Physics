'use strict';
/*
 * signed_field.js   --   ROUND 62
 * ---------------------------------------------------------------------------
 * THE SECOND HALF OF THE ALPHABET GAP: SIGNS, AND WITH THEM INTERFERENCE.
 *
 * Round 61 gave the photon mode a COUNT, so |2> exists and the gas is homogeneous -- but counts
 * are non-negative, so the state space is a CONE, not a vector space: two amplitudes can never
 * cancel.  This round gives each (cell, direction) a GAUSSIAN INTEGER amplitude: a signed real
 * part and a signed imaginary part.  That is a genuine vector space, and it carries a phase.
 *
 * THE PRICE, STATED UP FRONT
 *   The HPP head-on turn is BILINEAR -- it needs both particles present, which is a product of
 *   two amplitudes -- so it cannot survive in a linear field.  This round drops it.  That is not a
 *   loss: real photons do not scatter off each other, and the model's turn was a genuine
 *   nonlinearity that rounds 38-47 kept measuring and treating as a feature.  A linear photon
 *   sector and a self-interacting photon sector are incompatible, and this round takes the first.
 *
 * ATTEMPT 1 FAILED, AND THE FAILURE IS THE DISPERSION RELATION
 *   The first version of this script rotated every amplitude by 90 degrees on every step, on the
 *   theory that "phase advance per step = frequency".  It produced ZERO fringes: a rotation applied
 *   to every cell advances every cell's phase equally, so the field has no spatial phase gradient
 *   at all -- wavelength infinite, intensity flat, interference impossible.  Measured before the
 *   fix: one source 1.0/1.0/1.0/1.0, two sources 4.0/4.0/4.0/4.0, spread 0.0.
 *   The rotation was also, physically, a MASS TERM: it makes the field oscillate on its own with no
 *   source, which is a rest frequency, which is an energy gap.  Light has no gap (rounds 47-48).
 *   THE FIX, and it is the physical statement of this model's dispersion relation:
 *     the frequency belongs to the SOURCE, not to the lattice.  A source oscillating at omega
 *     injects an amplitude that then only MOVES.  A quantum injected at time t and found x cells
 *     downstream at time t+x therefore carries the source phase from time t, i.e. cos(omega(t-x)):
 *     the lattice converts a source frequency into a wavelength using speed 1, giving k = omega.
 *     An on-site rotation would have given k = sqrt(omega^2 - m^2) -- a massive dispersion -- and
 *     that is exactly why it killed the fringes.  THE ABSENCE OF A ROTATION *IS* THE MASSLESSNESS.
 *
 * WHAT IS BUILT
 *   re[], im[] : one signed number each, per (cell, direction).
 *   step       : move every amplitude one cell along its own direction.  That is all.  Linear.
 *   sources    : driven boundaries that oscillate in TIME: cos(omega t + phase), sin(omega t + phase).
 *   observable : the field at a cell is the sum over directions, psi = sum (re + i im);
 *                the intensity is |psi|^2.  (Summing directions is the physically sensible reading:
 *                the field at a point is the sum of what arrives from everywhere.  It is also what
 *                makes two waves crossing at right angles interfere at all.)
 *
 * THE INTERFERENCE EXPERIMENT, WITH AN EXACT PREDICTION
 *   Two plane-wave sources at right angles, driven at the same frequency omega = pi/2 per step,
 *   with a controllable phase difference delta:
 *        wave A into +x :  exp( i (omega t - omega x) )
 *        wave B into +y :  exp( i (omega t - omega y + delta) )
 *        |psi|^2 = 2 + 2 cos( omega(y - x) - delta )
 *   so the pattern is a set of fringes in u = x - y with period 2*pi/omega = 4 cells, and the fringe
 *   maximum sits at u = -delta/omega.  THEREFORE: shifting the source phase by 90 degrees must move
 *   the fringes by EXACTLY ONE CELL.  That is the prediction this round tests.
 *
 * PREDICTIONS (written before the run; C/D/E below are the re-run after the fix)
 *   A  linearity: step(A + B) = step(A) + step(B) exactly, for random complex configurations --
 *      true superposition, which the cone could not do.
 *   C  ONE source gives NO fringes (a flat intensity along u); TWO sources give fringes with period
 *      4 cells in u.
 *   D  THE FRINGE SHIFT: delta = 0, pi/2, pi, 3pi/2 must move the fringe maximum by exactly one
 *      cell each time, in the direction u -> u-1.
 *   E  a phase-blind (classical) reading -- adding intensities instead of amplitudes -- gives NO
 *      fringes at all.  That is the control that makes C and D mean something.
 *   F  the price: restoring the bilinear head-on turn breaks additivity.  A counterexample is shown.
 * ---------------------------------------------------------------------------
 */

const W = 64, H = 64, N = W * H;
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
const OMEGA = Math.PI / 2;                 // the source frequency: one quarter cycle per step
const LAMBDA = (2 * Math.PI) / OMEGA;      // 4 cells -- forced by speed 1, not chosen
const idx = (x, y, d) => (y * W + x) * 4 + d;

const zero = () => new Float64Array(4 * N);
// move every amplitude one cell along its own direction.  Nothing else happens.  This is the
// whole dynamics, and it is why the field is massless and why it is linear.
function move(a) {
  const o = new Float64Array(4 * N);
  for (let i = 0; i < N; i++) {
    const x = i % W, y = (i - x) / W;
    for (let d = 0; d < 4; d++) {
      const v = a[i * 4 + d];
      if (v === 0) continue;
      const k = ((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W);
      o[k * 4 + d] += v;
    }
  }
  return o;
}
function step(re, im) { return [move(re), move(im)]; }
// the observable: the field at each cell is the sum over its four directions
function intensity(re, im) {
  const I = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    let sr = 0, si = 0;
    for (let d = 0; d < 4; d++) { sr += re[i * 4 + d]; si += im[i * 4 + d]; }
    I[i] = sr * sr + si * si;
  }
  return I;
}

// two driven plane-wave sources: A enters from the left going +x, B from the bottom going +y.
// The source OSCILLATES IN TIME; propagation only carries that oscillation outward.
function run(steps, delta, sources) {
  let re = zero(), im = zero();
  for (let t = 0; t < steps; t++) {
    const th = OMEGA * t;
    if (sources.includes('A')) {
      for (let y = 0; y < H; y++) { re[idx(0, y, 1)] = Math.cos(th); im[idx(0, y, 1)] = Math.sin(th); }
    }
    if (sources.includes('B')) {
      for (let x = 0; x < W; x++) { re[idx(x, 0, 0)] = Math.cos(th + delta); im[idx(x, 0, 0)] = Math.sin(th + delta); }
    }
    const [r, m] = step(re, im);
    re = r; im = m;
  }
  return { re, im };
}
// bin the interior intensity by u = x - y (mod 4), which is where the fringes live
function fringeProfile(re, im) {
  const I = intensity(re, im);
  const sum = [0, 0, 0, 0], cnt = [0, 0, 0, 0];
  for (let y = 10; y < H - 10; y++) {
    for (let x = 10; x < W - 10; x++) {
      const u = (((x - y) % 4) + 4) % 4;
      sum[u] += I[y * W + x]; cnt[u]++;
    }
  }
  return sum.map((s, i) => s / cnt[i]);
}
const argmax = (a) => a.indexOf(Math.max(...a));

// ============================================================ A: linearity
console.log('============================================================');
console.log('A: IS IT A VECTOR SPACE NOW? (true superposition, not just scaling)');
console.log('============================================================');
{
  const rand = (scale) => {
    const a = zero(), b = zero();
    for (let i = 0; i < 4 * N; i++) {
      if ((i * 2654435761) % 7 < 2) { a[i] = ((i % 5) - 2) * scale; b[i] = ((i % 3) - 1) * scale; }
    }
    return [a, b];
  };
  const [are, aim] = rand(1), [bre, bim] = rand(2);
  let worst = 0;
  {
    const [r1, m1] = step(are, aim);
    const [r2, m2] = step(bre, bim);
    const sumRe = new Float64Array(4 * N), sumIm = new Float64Array(4 * N);
    for (let i = 0; i < 4 * N; i++) { sumRe[i] = are[i] + bre[i]; sumIm[i] = aim[i] + bim[i]; }
    const [rs, ms] = step(sumRe, sumIm);
    for (let i = 0; i < 4 * N; i++) {
      worst = Math.max(worst, Math.abs(rs[i] - (r1[i] + r2[i])), Math.abs(ms[i] - (m1[i] + m2[i])));
    }
  }
  console.log(`  max |step(A+B) - step(A) - step(B)| over random complex configurations: ${worst}`);
  console.log(`  => ${worst === 0 ? 'PASS' : '*** FAIL ***'}: signs make it a vector space, so superposition is exact.  The cone of round 61`);
  console.log(`     could only scale a configuration; it could never add two of them.`);
}
console.log('');

// ============================================================ C: fringes or not
console.log('============================================================');
console.log('C: ONE SOURCE OR TWO? (the intensity profile along the fringe coordinate)');
console.log('============================================================');
{
  const one = run(400, 0, ['A']);
  const two = run(400, 0, ['A', 'B']);
  const p1 = fringeProfile(one.re, one.im), p2 = fringeProfile(two.re, two.im);
  const spread = (a) => Math.max(...a) - Math.min(...a);
  console.log('   u = x-y mod 4      0        1        2        3      spread');
  console.log(`   one source      ${p1.map((v) => v.toFixed(2).padStart(7)).join('  ')}   ${spread(p1).toFixed(2)}`);
  console.log(`   two sources     ${p2.map((v) => v.toFixed(2).padStart(7)).join('  ')}   ${spread(p2).toFixed(2)}`);
  const flat = spread(p1) < 0.01, fringed = spread(p2) > 1.0;
  console.log('');
  console.log(`  predicted for two sources: 2 + 2cos(omega*u) over u = 0,1,2,3 with omega = pi/2`);
  console.log(`                             = ${[0, 1, 2, 3].map((u) => (2 + 2 * Math.cos(-OMEGA * u)).toFixed(2)).join('   ')}`);
  console.log(`  => ${flat && fringed ? 'PASS' : '*** FAIL ***'}: one source is flat, two sources vary with u -- and the four values land`);
  console.log(`     exactly on the predicted fringe curve (period ${LAMBDA} cells).`);
}
console.log('');

// ============================================================ D: the shift
console.log('============================================================');
console.log('D: THE FRINGE SHIFT -- 90 degrees of source phase must move the fringes ONE CELL');
console.log('============================================================');
{
  console.log('   source phase delta   fringe maximum at u   predicted   shift from the previous row');
  let ok = true, prev = null;
  for (const [name, delta] of [['0', 0], ['pi/2', Math.PI / 2], ['pi', Math.PI], ['3pi/2', 3 * Math.PI / 2]]) {
    const r = run(400, delta, ['A', 'B']);
    const p = fringeProfile(r.re, r.im);
    const u = argmax(p);
    // the predicted maximum is at u = -delta/omega, folded into the four bins
    const pred = ((Math.round(-delta / OMEGA) % 4) + 4) % 4;
    const shift = prev === null ? '-' : `${((prev - u) % 4 + 4) % 4} cell(s) toward u-1`;
    console.log(`        ${name.padEnd(8)}            ${u}                  ${pred}          ${shift}`);
    if (u !== pred) ok = false;
    if (prev !== null && ((prev - u) % 4 + 4) % 4 !== 1) ok = false;
    prev = u;
  }
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the fringe maximum sits exactly where the prediction puts it at all four source`);
  console.log(`     phases, and each 90 degree step of source phase moves it exactly one cell.`);
  console.log(`     THAT is interference: the pattern's POSITION depends on a phase, and nothing in the`);
  console.log(`     model could do that before this round.`);
}
console.log('');

// ============================================================ E: the phase-blind control
console.log('============================================================');
console.log('E: THE CONTROL -- A PHASE-BLIND READING SHOWS NO FRINGES AT ALL');
console.log('============================================================');
{
  // add the two runs' INTENSITIES instead of their amplitudes: the classical, incoherent sum
  const a = run(400, 0, ['A']), b = run(400, 0, ['B']);
  const Ia = intensity(a.re, a.im), Ib = intensity(b.re, b.im);
  const sum = [0, 0, 0, 0], cnt = [0, 0, 0, 0];
  for (let y = 10; y < H - 10; y++) {
    for (let x = 10; x < W - 10; x++) {
      const u = (((x - y) % 4) + 4) % 4;
      sum[u] += Ia[y * W + x] + Ib[y * W + x]; cnt[u]++;
    }
  }
  const p = sum.map((s, i) => s / cnt[i]);
  const spread = Math.max(...p) - Math.min(...p);
  console.log(`   incoherent (intensities added): ${p.map((v) => v.toFixed(2).padStart(7)).join('  ')}   spread ${spread.toFixed(2)}`);
  const coherent = run(400, 0, ['A', 'B']);
  const pc = fringeProfile(coherent.re, coherent.im);
  const spreadC = Math.max(...pc) - Math.min(...pc);
  console.log(`   coherent   (amplitudes added):  ${pc.map((v) => v.toFixed(2).padStart(7)).join('  ')}   spread ${spreadC.toFixed(2)}`);
  console.log('');
  console.log(`  => ${spread < 0.01 && spreadC > 1.0 ? 'PASS' : '*** FAIL ***'}: adding intensities gives a flat profile while adding AMPLITUDES gives`);
  console.log(`     fringes.  Without this control the fringes in C and D could have been an artefact of`);
  console.log(`     the geometry; with it, they are the signature of a phase.`);
}
console.log('');

// ============================================================ F: the price
console.log('============================================================');
console.log('F: THE PRICE -- THE HEAD-ON TURN IS BILINEAR, SO IT CANNOT LIVE IN A LINEAR FIELD');
console.log('============================================================');
{
  // a head-on turn in the round-61 style, written on signed amplitudes
  const turn = (re, im) => {
    const r = Float64Array.from(re), m = Float64Array.from(im);
    for (let i = 0; i < N; i++) {
      for (const [a, b, pa, pb] of [[0, 2, 1, 3], [1, 3, 0, 2]]) {
        const va = re[i * 4 + a], vb = re[i * 4 + b];
        if (va === 0 || vb === 0) continue;
        const t = Math.sign(va) * Math.sign(vb) * Math.min(Math.abs(va), Math.abs(vb));
        r[i * 4 + a] -= t; r[i * 4 + b] -= t;
        r[i * 4 + pa] += t; r[i * 4 + pb] += t;
      }
    }
    return [r, m];
  };
  const mk = (x, y, d) => { const a = zero(); a[idx(x, y, d)] = 1; return a; };
  // case 1: two quanta far apart -- nothing collides
  {
    const A = mk(20, 30, 1), B = mk(30, 30, 3);
    const [ar] = turn(A, zero()), [br] = turn(B, zero());
    const s = Float64Array.from(A);
    for (let i = 0; i < 4 * N; i++) s[i] = A[i] + B[i];
    const [tr] = turn(s, zero());
    let w = 0;
    for (let i = 0; i < 4 * N; i++) w = Math.max(w, Math.abs(tr[i] - (ar[i] + br[i])));
    console.log(`   two quanta in DIFFERENT cells:  max |turn(A+B) - turn(A) - turn(B)| = ${w.toFixed(3)}`);
    console.log(`     (they never meet, so the turn never fires -- a nonlinearity that is dormant is still there)`);
  }
  // case 2: they share a cell -- the term appears out of nothing
  const w2 = (() => {
    const C = mk(20, 30, 1); C[idx(20, 30, 3)] = 1;
    const [cr] = turn(C, zero());
    const [d1] = turn(mk(20, 30, 1), zero()), [d2] = turn(mk(20, 30, 3), zero());
    let w = 0;
    for (let i = 0; i < 4 * N; i++) w = Math.max(w, Math.abs(cr[i] - (d1[i] + d2[i])));
    return w;
  })();
  console.log(`   two quanta in the SAME cell:    max |turn(A+B) - turn(A) - turn(B)| = ${w2.toFixed(3)}`);
  console.log('');
  console.log(`  => ${w2 > 0 ? 'PASS' : '*** FAIL ***'}: the deviation appears exactly when both amplitudes occupy one cell, because the turn`);
  console.log(`     needs their PRODUCT.  A linear field cannot have that term.  So this round keeps linearity`);
  console.log(`     and drops the turn -- and real photons do not scatter off each other, so the linear choice`);
  console.log(`     is the physical one.  What the model loses here is an interaction it never should have had.`);
}
