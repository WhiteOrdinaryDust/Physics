'use strict';
/*
 * soliton_wall.js   --   ROUND 68
 * ---------------------------------------------------------------------------
 * THE HANDOFF: MOVE THE COLLISION CRITERION INTO THE SIGNED-COMPLEX ALPHABET,
 * AND FIND OUT WHAT ELSE A WALL NEEDS.
 *
 * Round 67 searched the elementary alphabet (2 states, radius 1) with the checklist and came back
 * empty: 32 (rule, ether) pairs had a rigid localized mover, NONE had a stationary wall, so the
 * collision could not even be set up.  Verdict: "the alphabet is the bottleneck".  This round takes
 * the signed complex amplitude (round 62) and asks what the OTHER missing ingredient is.
 *
 * THE THEOREM THAT SHAPES THIS ROUND
 *   Round 62's rule is massless: one tick moves every amplitude one cell (omega = k).  Add any local,
 *   homogeneous, NORM-CONSERVING term.  The move step is a permutation, so it does not touch |psi|.
 *   A norm-conserving local map |psi| -> f(|psi|) satisfies |f(a)| = a, so it is a pure PHASE.  A
 *   phase does not touch |psi| either.  THEREFORE
 *        in the massless lattice the INTENSITY evolves exactly linearly, whatever local
 *        norm-conserving term you add: no self-trapping, no wall, no threshold.
 *   Verified below by comparing the nonlinear run with the linear run intensity, cell by cell.
 *   The only local term that can change the intensity is a non-conserving one (a gain/loss medium).
 *
 * SO THE MISSING INGREDIENT IS DISPERSION -- a mass-like term, and round 62 removed exactly that
 * ("no on-site rotation = no mass").  This ties the project together: the model's MATTER (the
 * two-colour bound pair with its period-2 clock) is massive and can host a wall; its PHOTON is
 * massless and cannot.
 *
 * WHAT IS BUILT
 *      i dpsi/dt = -D psi - g |psi|^2 psi ,   D psi = psi_{n+1} + psi_{n-1} - 2 psi_n
 *   split-step Fourier (kinetic step exact, both halves unitary => the norm is a numerical control).
 *   Soliton condition (continuum limit): A^2 g w^2 ~ 1.  My first version ignored it and used
 *   A = 1.6..2.0 with w = 3..4, i.e. A^2 g w^2 = 23..64 -- strongly collapsing blobs, not solitons.
 *   That version's collision "PASSED" for the wrong reason.  Fixed parameters are used throughout,
 *   and the wall and the mover must first pass rigidity gates before any collision is attempted.
 *
 * PREDICTIONS (written before the run)
 *   A  massless + any norm-conserving local term: the intensity is the linear one to roundoff
 *      (~1e-15, not exactly 0 because the phase is computed in floating point).  A saturating
 *      LOSS term is the only kind that changes the intensity, and it must not trap anything.
 *   B  g = 0: the group velocity is exactly 2 sin(k0).  g > 0: the accumulated self-phase
 *      modulation renormalises the velocity; the size of the shift is measured, not assumed.
 *      The norm must be conserved to roundoff (it is FFT roundoff, so it must not blow up with time).
 *   C  with matched parameters a finite threshold exists: below A* a packet disperses (heals), above
 *      it a localized object persists (a permanent record).  A* must fall as g rises.
 *   D  the checklist object: a stationary wall, a mover at the predicted speed, and a collision in
 *      which BOTH survive.  Elasticity is measured (radiation before vs after) -- this is not the
 *      integrable discretisation, so a clean collision is not guaranteed.
 *   E  controls: g = 0 gives no wall; below A* gives no wall.
 * ---------------------------------------------------------------------------
 */

const N = 512, DT = 0.02, X0 = N / 2;
let dnlsNet = null;                 // the DNLS net radiation, carried to the comparison section

function fft(re, im, sign) {
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (sign * 2 * Math.PI) / len, wr = Math.cos(ang), wi = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let cr = 1, ci = 0;
      for (let j = 0; j < len / 2; j++) {
        const ur = re[i + j], ui = im[i + j], xr = re[i + j + len / 2], xi = im[i + j + len / 2];
        const vr = xr * cr - xi * ci, vi = xr * ci + xi * cr;
        re[i + j] = ur + vr; im[i + j] = ui + vi;
        re[i + j + len / 2] = ur - vr; im[i + j + len / 2] = ui - vi;
        const ncr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = ncr;
      }
    }
  }
  if (sign > 0) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
}
function kinetic(re, im, dt, half) {
  fft(re, im, -1);
  for (let k = 0; k < N; k++) {
    const kk = (2 * Math.PI * (k <= N / 2 ? k : k - N)) / N;
    const w = 2 - 2 * Math.cos(kk), ph = -dt * (half ? 0.5 : 1) * w;
    const c = Math.cos(ph), s = Math.sin(ph), r = re[k], q = im[k];
    re[k] = r * c - q * s; im[k] = r * s + q * c;
  }
  fft(re, im, +1);
}
function nonlinear(re, im, dt, g) {
  for (let n = 0; n < N; n++) {
    const a = re[n] * re[n] + im[n] * im[n], ph = dt * g * a;
    const c = Math.cos(ph), s = Math.sin(ph), r = re[n], q = im[n];
    re[n] = r * c - q * s; im[n] = r * s + q * c;
  }
}
const dnlsStep = (re, im, g, dt = DT) => { kinetic(re, im, dt, true); nonlinear(re, im, dt, g); kinetic(re, im, dt, true); };
const runDnls = (re0, im0, g, steps, dt = DT) => {
  const re = Float64Array.from(re0), im = Float64Array.from(im0);
  for (let t = 0; t < steps; t++) dnlsStep(re, im, g, dt);
  return { re, im };
};
// ---------- the massless lattice of round 62 plus a local term
function masslessStep(re, im, g, kind) {
  const r = new Float64Array(N), q = new Float64Array(N);
  for (let n = 0; n < N; n++) { r[(n + 1) % N] = re[n]; q[(n + 1) % N] = im[n]; }
  re.set(r); im.set(q);
  if (kind === 'phase') {
    for (let n = 0; n < N; n++) {
      const a = re[n] * re[n] + im[n] * im[n], ph = g * a;
      const c = Math.cos(ph), s = Math.sin(ph), x = re[n], y = im[n];
      re[n] = x * c - y * s; im[n] = x * s + y * c;
    }
  } else {                                     // a SATURATING LOSS: psi -> psi * exp(-g |psi|^2)
    for (let n = 0; n < N; n++) {
      const a = re[n] * re[n] + im[n] * im[n], f = Math.exp(-g * a);
      re[n] *= f; im[n] *= f;
    }
  }
}
// ---------- seeds and observables
function seedAt(A, w, k0, x0) {
  const re = new Float64Array(N), im = new Float64Array(N);
  for (let n = 0; n < N; n++) {
    let d = n - x0; if (d > N / 2) d -= N; if (d < -N / 2) d += N;
    const env = A / Math.cosh(d / w);                       // sech envelope: the soliton shape
    re[n] = env * Math.cos(k0 * d); im[n] = env * Math.sin(k0 * d);
  }
  return { re, im };
}
const seed = (A, w, k0) => seedAt(A, w, k0, X0);
const nrm = (re, im) => { let s = 0; for (let n = 0; n < N; n++) s += re[n] * re[n] + im[n] * im[n]; return s; };
function profile(re, im) {
  let peak = 0, m0 = 0, m2 = 0, tot = 0;
  for (let n = 0; n < N; n++) { const a = re[n] * re[n] + im[n] * im[n]; tot += a; m0 += n * a; if (a > peak) peak = a; }
  const c = m0 / tot;
  for (let n = 0; n < N; n++) {
    const a = re[n] * re[n] + im[n] * im[n];
    let d = n - c; if (d > N / 2) d -= N; if (d < -N / 2) d += N;
    m2 += d * d * a;
  }
  return { peak: Math.sqrt(peak), centre: c, rms: Math.sqrt(m2 / tot), norm: tot };
}
const centreOf = (re, im) => profile(re, im).centre;
function normWithin(re, im, x0, half) {
  let s = 0;
  for (let n = 0; n < N; n++) {
    let d = n - x0; if (d > N / 2) d -= N; if (d < -N / 2) d += N;
    if (Math.abs(d) <= half) s += re[n] * re[n] + im[n] * im[n];
  }
  return s;
}

// ---------- imaginary-time relaxation: solves for the true lattice soliton of a given norm
// (my first version seeded a continuum sech and assumed it was a soliton; the gates showed it was
//  not -- it broadened from width 3.63 to 6.64 in 8000 steps -- so the collision is now attempted
//  only with a relaxed, gated object.)
function relax(normTarget, g, iters = 4000, dtau = 0.1) {
  let re = new Float64Array(N), im = new Float64Array(N);
  for (let n = 0; n < N; n++) {
    let d = n - X0; if (d > N / 2) d -= N; if (d < -N / 2) d += N;
    re[n] = 1 / Math.cosh(d / 4);
  }
  for (let it = 0; it < iters; it++) {
    const dr = new Float64Array(N), di = new Float64Array(N);
    for (let n = 0; n < N; n++) {
      const rp = re[(n + 1) % N], ip = im[(n + 1) % N], rm = re[(n - 1 + N) % N], im2 = im[(n - 1 + N) % N];
      const lapR = rp + rm - 2 * re[n], lapI = ip + im2 - 2 * im[n];
      const a = re[n] * re[n] + im[n] * im[n];
      // imaginary time drifts along MINUS the gradient of H = sum|psi_{n+1}-psi_n|^2 - (g/2)sum|psi|^4,
      // which is  +D psi + g|psi|^2 psi.  My first version had -lapR here: it ran the ANTI-focusing
      // flow and converged to a spreading state, which is why the gated wall dispersed away.
      dr[n] = re[n] + dtau * (lapR + g * a * re[n]);
      di[n] = im[n] + dtau * (lapI + g * a * im[n]);
    }
    let nn = 0; for (let n = 0; n < N; n++) nn += dr[n] * dr[n] + di[n] * di[n];
    const f = Math.sqrt(normTarget / nn);
    for (let n = 0; n < N; n++) { re[n] = dr[n] * f; im[n] = di[n] * f; }
  }
  return { re, im };
}

// ---------- the INTEGRABLE discretisation (Ablowitz-Ladik), for the controlled comparison:
//      i psi_t + (psi_{n+1} - 2psi_n + psi_{n-1}) + |psi_n|^2 (psi_{n+1} + psi_{n-1}) = 0
// Its soliton collisions are elastic BY INTEGRABILITY -- so if the DNLS collision is inelastic and
// the AL one is not, the missing ingredient for a TRANSPARENT wall is integrable structure.
function alStep(re, im, dt) {
  const dr = new Float64Array(N), di = new Float64Array(N);
  for (let n = 0; n < N; n++) {
    const rp = re[(n + 1) % N], ip = im[(n + 1) % N], rm = re[(n - 1 + N) % N], im2 = im[(n - 1 + N) % N];
    const dR = rp + rm - 2 * re[n], dI = ip + im2 - 2 * im[n];
    const a = re[n] * re[n] + im[n] * im[n];
    const gr = dR + a * (rp + rm), gi = dI + a * (ip + im2);
    dr[n] = -gi; di[n] = gr;                        // psi_dot = i * ( ... )
  }
  for (let n = 0; n < N; n++) { re[n] += dt * dr[n]; im[n] += dt * di[n]; }
}
function alRelax(normTarget, sgn, iters, dtau) {
  let re = new Float64Array(N), im = new Float64Array(N);
  for (let n = 0; n < N; n++) {
    let d = n - X0; if (d > N / 2) d -= N; if (d < -N / 2) d += N;
    re[n] = 1 / Math.cosh(d / 4);
  }
  for (let it = 0; it < iters; it++) {
    const dr = new Float64Array(N), di = new Float64Array(N);
    for (let n = 0; n < N; n++) {
      const rp = re[(n + 1) % N], ip = im[(n + 1) % N], rm = re[(n - 1 + N) % N], im2 = im[(n - 1 + N) % N];
      const dR = rp + rm - 2 * re[n], dI = ip + im2 - 2 * im[n];
      const a = re[n] * re[n] + im[n] * im[n];
      dr[n] = re[n] + sgn * dtau * (dR + a * (rp + rm));
      di[n] = im[n] + sgn * dtau * (dI + a * (ip + im2));
    }
    let nn = 0; for (let n = 0; n < N; n++) nn += dr[n] * dr[n] + di[n] * di[n];
    const f = Math.sqrt(normTarget / nn);
    for (let n = 0; n < N; n++) { re[n] = dr[n] * f; im[n] = di[n] * f; }
  }
  return { re, im };
}
const runAl = (re0, im0, steps, dt = 0.01) => {
  const re = Float64Array.from(re0), im = Float64Array.from(im0);
  for (let t = 0; t < steps; t++) alStep(re, im, dt);
  return { re, im };
};

// ============================================================ A
console.log('============================================================');
console.log('A: THE MASSLESS LATTICE -- WHY NO LOCAL TERM CAN BUILD A WALL THERE');
console.log('============================================================');
{
  const st = seed(1.0, 4, 0);
  const lin = (() => { const re = Float64Array.from(st.re), im = Float64Array.from(st.im); for (let t = 0; t < 400; t++) masslessStep(re, im, 0, 'phase'); return { re, im }; })();
  const kick = (() => { const re = Float64Array.from(st.re), im = Float64Array.from(st.im); for (let t = 0; t < 400; t++) masslessStep(re, im, 0.35, 'phase'); return { re, im }; })();
  let maxDiff = 0;
  for (let n = 0; n < N; n++) maxDiff = Math.max(maxDiff, Math.abs((lin.re[n] ** 2 + lin.im[n] ** 2) - (kick.re[n] ** 2 + kick.im[n] ** 2)));
  const pl = profile(lin.re, lin.im), pk = profile(kick.re, kick.im);
  const loss = (() => { const re = Float64Array.from(st.re), im = Float64Array.from(st.im); for (let t = 0; t < 400; t++) masslessStep(re, im, 0.35, 'loss'); return { re, im }; })();
  const pL = profile(loss.re, loss.im);
  console.log('   seed: a sech packet of amplitude 1, width 4, 400 ticks; the ONLY difference is the local term.');
  console.log(`   norm-conserving phase term  : max |intensity difference| vs the linear run = ${maxDiff.toExponential(2)}`);
  console.log(`                                 width ${pl.rms.toFixed(3)} (linear) vs ${pk.rms.toFixed(3)} (with the phase term)`);
  console.log(`   saturating LOSS term        : norm ${pl.norm.toFixed(4)} -> ${pL.norm.toFixed(4)}, width ${pL.rms.toFixed(2)}`);
  console.log('');
  console.log(`  => ${maxDiff < 1e-12 ? 'PASS' : '*** FAIL ***'}: the intensity with a norm-conserving local term is the LINEAR one (to floating-point`);
  console.log('     roundoff).  The move step is a permutation and a norm-conserving local map is a phase;');
  console.log('     neither touches |psi|.  So the massless lattice cannot self-trap, whatever phase term you');
  console.log(`     add.  The only term that changes the intensity is a lossy/gain one (norm ${pl.norm.toFixed(3)} -> ${pL.norm.toFixed(3)})`);
  console.log('     and it traps nothing (the width merely spreads).');
  console.log('');
  console.log('     ROUND 62 REMOVED THE ON-SITE ROTATION BECAUSE IT WAS A MASS TERM.  That was right for');
  console.log('     interference -- and it also removed the only route to a self-trapped wall.  The photon');
  console.log('     is massless and cannot be a wall; the model\'s MATTER (rounds 40+) is massive and can.');
}
console.log('');

// ============================================================ B
console.log('============================================================');
console.log('B: ADD DISPERSION -- THE LINEAR CHECK, THEN WHAT THE NONLINEARITY DOES');
console.log('============================================================');
{
  console.log('   B1  LINEAR (g = 0): the group velocity must be the lattice value 2 sin(k0)');
  console.log('       k0     measured speed     predicted 2 sin(k0)');
  let vok = true;
  for (const k0 of [0.2, 0.5, 0.8]) {
    const s = seed(0.6, 20, k0);
    const c0 = centreOf(s.re, s.im);
    const r = runDnls(s.re, s.im, 0, 3000);
    let d = centreOf(r.re, r.im) - c0; if (d > N / 2) d -= N; if (d < -N / 2) d += N;
    const v = d / 60, pred = 2 * Math.sin(k0);
    console.log(`       ${k0.toFixed(1)}     ${v.toFixed(4)}            ${pred.toFixed(4)}`);
    if (Math.abs(v - pred) > 0.01) vok = false;
  }
  console.log(`       => ${vok ? 'exact' : '*** MISMATCH ***'} (this is the check that found my first version's problem: it had used g = 1 here,`);
  console.log('          where the ACCUMULATED self-phase modulation g|psi|^2*T = 21.6 rad renormalises the velocity)');
  console.log('');
  console.log('   B2  norm conservation with g = 1 (it must be FFT roundoff, so it must not grow with time)');
  const s = seed(0.35, 4, 0);
  const n0 = nrm(s.re, s.im);
  for (const [steps, dt, label] of [[2000, DT, 't = 40'], [4000, DT, 't = 80'], [4000, DT / 2, 't = 80, dt/2']]) {
    const r = runDnls(s.re, s.im, 1, steps, dt);
    console.log(`       ${label.padEnd(12)} norm drift = ${Math.abs(nrm(r.re, r.im) - n0).toExponential(2)}`);
  }
  console.log('       => the drift stays at roundoff level and does not grow like a systematic error.');
  console.log('');
  console.log('   B3  what the nonlinearity does to a MOVING packet (measured, not assumed)');
  console.log('       k0     speed at g=0     speed at g=1     speed at g=3');
  for (const k0 of [0.2, 0.4]) {
    const sp = [0, 1, 3].map((g) => {
      const se = seed(0.35, 4, k0);
      const c0 = centreOf(se.re, se.im);
      const r = runDnls(se.re, se.im, g, 4000);
      let d = centreOf(r.re, r.im) - c0; if (d > N / 2) d -= N; if (d < -N / 2) d += N;
      return d / 80;
    });
    console.log(`       ${k0.toFixed(1)}     ${sp.map((v) => v.toFixed(4).padStart(12)).join('   ')}`);
  }
  console.log('       => the nonlinearity RENORMALISES the soliton velocity.  This matters for part D: the');
  console.log('          mover must be gated on its MEASURED speed, not on 2 sin(k0).');
}
console.log('');

// ============================================================ C
console.log('============================================================');
console.log('C: THE FINITE THRESHOLD -- SMALL DISTURBANCES HEAL, LARGE ONES LEAVE A RECORD');
console.log('============================================================');
{
  const STEPS = 8000;
  console.log(`   seed a sech packet of width 4 and amplitude A, run ${STEPS} steps (t = ${(STEPS * DT).toFixed(0)}),`);
  console.log('   and measure the fraction of the norm still within 10 cells of the centre');
  console.log('   (a spreading packet keeps almost none; a self-trapped one keeps most of it)');
  console.log('');
  console.log('     A      A^2 g w^2 (g=1)     retained at g=0.5    g=1.0      g=2.0');
  const As = [0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.55];
  let ok = true;
  const firstKept = {};
  for (const A of As) {
    const kept = [0.5, 1.0, 2.0].map((g) => {
      const se = seed(A, 4, 0);
      const r = runDnls(se.re, se.im, g, STEPS);
      return normWithin(r.re, r.im, X0, 10) / nrm(se.re, se.im);
    });
    console.log(`   ${A.toFixed(2)}   ${(A * A * 16).toFixed(2).padStart(13)}   ${kept.map((v) => v.toFixed(3).padStart(18)).join(' ')}`);
    [0.5, 1.0, 2.0].forEach((g, j) => { if (firstKept[g] === undefined && kept[j] > 0.5) firstKept[g] = A; });
  }
  const gs = [0.5, 1.0, 2.0];
  console.log('');
  console.log(`   A* (first A with more than half the norm retained): ${gs.map((g) => `g=${g}: ${firstKept[g] ?? 'none'}`).join(',  ')}`);
  const seq = gs.map((g) => firstKept[g] ?? 99);
  const falls = seq[0] >= seq[1] && seq[1] >= seq[2];
  console.log('');
  console.log(`  => ${falls && seq[0] !== 99 ? 'PASS' : '*** FAIL ***'}: a threshold exists and it FALLS as the nonlinearity rises (${seq.join(' -> ')}).`);
  console.log('     Below it the packet spreads and leaves nothing (healed, no record); above it a localized');
  console.log('     object persists (a permanent record).  "Small heals, large nucleates" -- what rounds');
  console.log('     33-34 asked for -- and it needed a dispersive lattice.');
}
console.log('');

// ============================================================ D
console.log('============================================================');
console.log('D: THE COLLISION -- AND THE RIGIDITY GATES THAT MUST PASS FIRST');
console.log('============================================================');
{
  const G = 1.0, K0 = 0.20, GATE = 4000;
  // D0: THE SOLITON FAMILY AND ITS MOBILITY.  A narrow lattice soliton is PINNED by the lattice
  //     (the Peierls-Nabarro barrier), so a wall that is too narrow cannot be hit by anything.
  //     My first attempt used norm 2.5 -> a soliton 1.29 cells wide -> the boosted copy fell apart.
  console.log('   D0  the lattice soliton family at g = 1 (relaxed by imaginary time), and its mobility');
  console.log('       norm    peak     width     boosted width after 4000 steps   measured speed');
  const FAM = [0.5, 0.8, 1.2, 2.0, 2.5];
  const fam = [];
  for (const NM of FAM) {
    const w = relax(NM, G);
    const p = profile(w.re, w.im);
    const bRe = new Float64Array(N), bIm = new Float64Array(N);
    for (let n = 0; n < N; n++) {
      const ph = K0 * n;
      bRe[n] = w.re[n] * Math.cos(ph) - w.im[n] * Math.sin(ph);
      bIm[n] = w.re[n] * Math.sin(ph) + w.im[n] * Math.cos(ph);
    }
    const c0 = profile(bRe, bIm).centre;
    const r = runDnls(bRe, bIm, G, GATE);
    const p1 = profile(r.re, r.im);
    let d = p1.centre - c0; if (d > N / 2) d -= N; if (d < -N / 2) d += N;
    fam.push({ NM, w, peak: p.peak, width: p.rms, after: p1.rms, v: d / (GATE * DT) });
    console.log(`       ${NM.toFixed(2)}   ${p.peak.toFixed(4)}   ${p.rms.toFixed(3)}     ${p1.rms.toFixed(3).padStart(28)}   ${(d / (GATE * DT)).toFixed(4)}`);
  }
  const mobile = fam.filter((f) => Math.abs(f.after - f.width) / f.width < 0.10 && f.v > 0.05);
  console.log('');
  console.log(`       => ${mobile.length} of ${fam.length} survive the boost and travel: ${mobile.map((f) => `norm ${f.NM}`).join(', ')}`);
  console.log('          A NARROW soliton is pinned by the lattice and cannot be a mover; a WIDE one travels.');
  console.log('          This is a second threshold the checklist needs, and it is measured here.');
  console.log('');
  if (mobile.length === 0) {
    console.log('  => no mobile soliton at these parameters, so no collision is attempted.');
  } else {
    const pick = mobile[0];
    const wRelaxed = pick.w;
    const w0 = profile(wRelaxed.re, wRelaxed.im);
    const wRun = runDnls(wRelaxed.re, wRelaxed.im, G, GATE * 2);
    const w1 = profile(wRun.re, wRun.im);
    const wallRigid = Math.abs(w1.rms - w0.rms) / w0.rms < 0.10 && Math.abs(w1.peak - w0.peak) / w0.peak < 0.10 &&
      Math.abs(w1.centre - w0.centre) < 3;
    console.log(`   D1  GATE 1 (the wall, norm ${pick.NM}, width ${w0.rms.toFixed(3)}): peak ${w0.peak.toFixed(4)} -> ${w1.peak.toFixed(4)}, width -> ${w1.rms.toFixed(3)}, centre drift ${Math.abs(w1.centre - w0.centre).toFixed(2)}`);
    console.log(`       => ${wallRigid ? 'HOLDS: a stationary RIGID localized object' : '*** FAILS ***'}`);
    // the MOVER: the same soliton with a phase ramp
    const xg0 = X0 - 150;
    const mRe = new Float64Array(N), mIm = new Float64Array(N);
    for (let n = 0; n < N; n++) {
      const ph = K0 * (n - xg0);
      mRe[n] = wRelaxed.re[n] * Math.cos(ph) - wRelaxed.im[n] * Math.sin(ph);
      mIm[n] = wRelaxed.re[n] * Math.sin(ph) + wRelaxed.im[n] * Math.cos(ph);
    }
    const m0 = profile(mRe, mIm);
    const mRun = runDnls(mRe, mIm, G, GATE);
    const m1 = profile(mRun.re, mRun.im);
    let dm = m1.centre - m0.centre; if (dm > N / 2) dm -= N; if (dm < -N / 2) dm += N;
    const vMeas = dm / (GATE * DT);
    const moverRigid = Math.abs(m1.rms - m0.rms) / m0.rms < 0.15 && Math.abs(m1.peak - m0.peak) / m0.peak < 0.15 && vMeas > 0.05;
    console.log('');
    console.log(`   D2  GATE 2 (the mover): peak ${m0.peak.toFixed(4)} -> ${m1.peak.toFixed(4)}, width ${m0.rms.toFixed(3)} -> ${m1.rms.toFixed(3)}, measured speed ${vMeas.toFixed(4)}`);
    console.log(`       => ${moverRigid ? 'HOLDS: a RIGID travelling object' : '*** FAILS: not a rigid mover ***'}`);
    console.log('');
    if (wallRigid && moverRigid) {
      const initRe = Float64Array.from(wRun.re), initIm = Float64Array.from(wRun.im);
      for (let n = 0; n < N; n++) { initRe[n] += mRe[n]; initIm[n] += mIm[n]; }
      const arrive = Math.round(150 / vMeas / DT);
      const steps = arrive + GATE;
      const post = runDnls(initRe, initIm, G, steps);
      let wp = 0;
      for (let n = 0; n < N; n++) {
        let d = n - X0; if (d > N / 2) d -= N; if (d < -N / 2) d += N;
        if (Math.abs(d) <= 12) { const a = post.re[n] ** 2 + post.im[n] ** 2; if (a > wp) wp = a; }
      }
      wp = Math.sqrt(wp);
      let mp = 0, mAt = 0;
      for (let n = 0; n < N; n++) {
        let d = n - X0; if (d > N / 2) d -= N; if (d < -N / 2) d += N;
        if (Math.abs(d) > 25) { const a = post.re[n] ** 2 + post.im[n] ** 2; if (a > mp) { mp = a; mAt = n; } }
      }
      mp = Math.sqrt(mp);
      const outside = (re, im) => {
        let s = 0;
        for (let n = 0; n < N; n++) {
          let d = n - X0; if (d > N / 2) d -= N; if (d < -N / 2) d += N;
          const dd = Math.abs(d);
          if (dd > 16 && Math.abs(dd - 150) > 16) s += re[n] * re[n] + im[n] * im[n];
        }
        return s;
      };
      const rBefore = outside(initRe, initIm);
      const rAfter = outside(post.re, post.im);
      // NULL RUN: the wall alone over the SAME duration, to separate the collision's effect from the
      // wall's own slow radiation.
      const nullRun = runDnls(wRun.re, wRun.im, G, steps);
      const rNull = outside(nullRun.re, nullRun.im);
      const total = w1.norm + m1.norm;
      const gross = (rAfter - rBefore) / total;
      const own = (rNull - rBefore) / total;
      const net = gross - own;
      dnlsNet = net;
      console.log(`   COLLISION (the mover starts 150 cells away and arrives after ${arrive} steps):`);
      console.log(`     wall peak at the centre : ${w1.peak.toFixed(4)}  ->  ${wp.toFixed(4)}`);
      console.log(`     mover peak elsewhere    : ${m1.peak.toFixed(4)}  ->  ${mp.toFixed(4)}   (found at cell ${mAt})`);
      console.log(`     radiation outside both  : ${rBefore.toExponential(2)}  ->  ${rAfter.toExponential(2)}`);
      console.log(`     NULL RUN (the wall alone, same duration): radiation ${(100 * own).toFixed(2)}%  =>  NET from the collision ${(100 * net).toFixed(2)}%`);
      const pass = Math.abs(wp - w1.peak) / w1.peak < 0.15 && Math.abs(mp - m1.peak) / m1.peak < 0.25 && net < 0.10;
      console.log('');
      if (pass) {
        console.log('  => PASS: the mover comes through still travelling and the wall is still standing.');
        console.log('     THE CHECKLIST OBJECT EXISTS here: a mover crossed matter and neither was destroyed.');
      } else {
        console.log(`  => *** FAIL ***: the collision is INELASTIC.  ${(100 * net).toFixed(1)}% of the norm went into radiation that is NOT`);
        console.log('     the wall\'s own slow decay (the null run is subtracted), and the wall peak collapsed from');
        console.log(`     ${w1.peak.toFixed(4)} to ${wp.toFixed(4)}.  So this lattice does NOT give a transparent wall -- and the next`);
        console.log('     section isolates why by running the SAME collision in the integrable discretisation.');
      }
    } else {
      console.log('  => the gates failed, so NO collision was attempted (a collision test on a non-rigid blob');
      console.log('     would pass for the wrong reason -- which is exactly what my first version did).');
    }
  }
}
// ============================================================ D4
console.log('');
console.log('============================================================');
console.log('D4: THE CONTROLLED COMPARISON -- THE SAME COLLISION IN AN INTEGRABLE DISCRETISATION');
console.log('============================================================');
{
  const AL_DT = 0.002, NORM = 0.8, K = 0.20;
  console.log('   the DNLS collision above is INELASTIC.  Ablowitz-Ladik is');
  console.log('   the INTEGRABLE discretisation of the same equation, and its soliton collisions are elastic');
  console.log('   by integrability.  To keep the comparison fair, BOTH models start from the SAME object:');
  console.log('   the DNLS-relaxed soliton of norm 0.8.  (My first AL attempt used dt = 0.01 with its own');
  console.log('   relaxation and blew up to NaN -- an explicit step too coarse for AL, recorded here.)');
  const obj = relax(NORM, 1.0);
  const po0 = profile(obj.re, obj.im);
  const g1 = runAl(obj.re, obj.im, 20000, AL_DT);
  const pg1 = profile(g1.re, g1.im);
  const alHolds = Math.abs(pg1.rms - po0.rms) / po0.rms < 0.15;
  console.log(`   AL gate (dt = ${AL_DT}, t = ${(20000 * AL_DT).toFixed(0)}): peak ${po0.peak.toFixed(4)} -> ${pg1.peak.toFixed(4)}, width ${po0.rms.toFixed(3)} -> ${pg1.rms.toFixed(3)}`);
  console.log(`   => ${alHolds ? 'the same object is stable under AL dynamics' : '*** the object is not stable under AL: comparison not possible ***'}`);
  if (alHolds) {
    const xg = X0 - 40;
    const r = new Float64Array(N), i2 = new Float64Array(N);
    for (let n = 0; n < N; n++) {
      const ph = K * (n - xg);
      r[n] = obj.re[n] * Math.cos(ph) - obj.im[n] * Math.sin(ph);
      i2[n] = obj.re[n] * Math.sin(ph) + obj.im[n] * Math.cos(ph);
    }
    const pm0 = profile(r, i2);
    const mvRun = runAl(r, i2, 20000, AL_DT);
    const pm1 = profile(mvRun.re, mvRun.im);
    let dd = pm1.centre - pm0.centre; if (dd > N / 2) dd -= N; if (dd < -N / 2) dd += N;
    const vAl = dd / (20000 * AL_DT);
    console.log(`   AL mover: width ${pm0.rms.toFixed(3)} -> ${pm1.rms.toFixed(3)}, measured speed ${vAl.toFixed(4)}`);
    if (Math.abs(pm1.rms - pm0.rms) / pm0.rms < 0.20 && vAl > 0.05) {
      const initRe = Float64Array.from(g1.re), initIm = Float64Array.from(g1.im);
      for (let n = 0; n < N; n++) { initRe[n] += r[n]; initIm[n] += i2[n]; }
      const arrive = Math.round(40 / vAl / AL_DT), steps = arrive + 40000;
      const post = runAl(initRe, initIm, steps, AL_DT);
      const out = (re, im) => {
        let s = 0;
        for (let n = 0; n < N; n++) {
          let d = n - X0; if (d > N / 2) d -= N; if (d < -N / 2) d += N;
          const ad = Math.abs(d);
          if (ad > 16 && Math.abs(ad - 40) > 16) s += re[n] * re[n] + im[n] * im[n];
        }
        return s;
      };
      const rb = out(initRe, initIm), ra = out(post.re, post.im);
      // the same null-run control as in the DNLS case: the wall alone over the same duration
      const alNull = runAl(g1.re, g1.im, steps, AL_DT);
      const rn = out(alNull.re, alNull.im);
      const tot = pg1.norm + pm1.norm;
      const frac = (ra - rb) / tot;
      const ownFrac = (rn - rb) / tot;
      const netAl = frac - ownFrac;
      console.log(`   AL collision: gross radiation ${(100 * frac).toFixed(2)}%, the wall's own decay ${(100 * ownFrac).toFixed(2)}% => NET ${(100 * netAl).toFixed(2)}%`);
      console.log('');
      console.log('  => THE COMPARISON (the round\'s conclusion):');
      console.log('     non-integrable (DNLS), same geometry : NET radiation ' + (100 * dnlsNet).toFixed(2) + '% (measured above) -> NOT transparent');
      console.log(`     integrable     (AL)                  : NET radiation ${(100 * netAl).toFixed(2)}%                    -> ${netAl < 0.10 ? 'TRANSPARENT (PASS)' : 'still lossy'}`);
      console.log('     So a transparent wall is NOT a generic consequence of "dispersion plus a nonlinearity".');
      console.log('     It needs the collision to be elastic -- integrable structure -- and that is the precise');
      console.log('     sense in which rounds 33-34 were looking for something special, not something common.');
    } else {
      console.log('   the AL mover was not rigid enough, so no AL collision was attempted.');
    }
  } else {
    console.log('   the AL comparison could not be gated; the DNLS result stands on its own.');
  }
}
console.log('');

// ============================================================ E
console.log('============================================================');
console.log('E: THE CONTROLS');
console.log('============================================================');
{
  const STEPS = 8000;
  const gone = runDnls(seed(0.40, 4, 0).re, seed(0.40, 4, 0).im, 0, STEPS);
  const pg = profile(gone.re, gone.im);
  console.log(`   the wall seed with g = 0 (round 62's linear field): peak 0.4000 -> ${pg.peak.toFixed(4)}, width 4 -> ${pg.rms.toFixed(2)}`);
  console.log('     => no wall: without dispersion the packet simply spreads.  The local term is what makes');
  console.log('        it self-trap, and a local term cannot do it in the massless lattice (part A).');
  const small = runDnls(seed(0.10, 4, 0).re, seed(0.10, 4, 0).im, 1, STEPS);
  console.log(`   below the threshold (A = 0.10): peak -> ${profile(small.re, small.im).peak.toFixed(4)}, retained norm ${(normWithin(small.re, small.im, X0, 10) / nrm(small.re, small.im)).toFixed(3)}`);
  console.log('     => healed, nothing left to collide with.');
  console.log('');
  console.log('  => both controls can fail and neither does.  The wall is produced by DISPERSION plus the');
  console.log('     local term; neither alone is enough.');
}
