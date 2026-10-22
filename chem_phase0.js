'use strict';
/*
 * chem_phase0.js  --  CHEMISTRY LINE, PHASE 0 (see PREREG_CHEMISTRY.md, written BEFORE this file)
 * ---------------------------------------------------------------------------
 * THE QUESTION (and only this): in a phi^4 substrate, can two "particles" collide and form a
 * REPRODUCIBLE bound product?  Not "is there chemistry" -- that is Phase 1, run only if this passes.
 *
 * THE MODEL is the RELATIVISTIC phi^4 (Klein-Gordon), not the gradient dynamics:
 *      phi_tt = phi_xx - (phi^3 - phi)
 * The gradient version has no inertia, so a kink only drifts and never collides -- which is why the
 * gradient substrates of rounds 97-125 could not have been chemistry substrates at all.
 *
 * THE PARTICLES are kinks (Q = +-1/2).  A kink and an antikink attract, so they always meet; the
 * interesting question is what they do afterwards, and whether the answer is a FUNCTION of (v, d).
 *
 * PRE-REGISTERED CRITERIA (fixed in PREREG_CHEMISTRY.md before running)
 *   H1 reproducibility: perturb the initial state by 1e-9 and re-run.
 *        flip rate <= 10%  -> the outcome table is a function of (v,d)  -> go to Phase 1
 *        flip rate >  30%  -> the scattering is chaotic                  -> PHASE 0 FAILS, stop
 *   H2 bound products: fraction of "bind" outcomes
 *        >= 20% -> bound products exist | < 5% -> stop | in between -> grey, report as partial
 *   OUTCOME CLASSES, fixed in advance:
 *        pass        |Q_final| > 0.4                      (the pair came back out as two kinks)
 *        bind        Q_final ~ 0 and separation < R_bind = 12 at the end
 *        radiated    Q_final ~ 0 and separation >= R_bind (annihilated into radiation)
 *   MY BET, recorded so it can be wrong: H1 FAILS (phi^4 kink-antikink scattering is known to be
 *   resonant/chaotic), and H2 is partially satisfied (low-v bions exist but are not reproducible).
 * ---------------------------------------------------------------------------
 */

const N = 1024, DX = 0.25, DT = 0.05;
const L = (N - 1) * DX;
const T_END = 300;
const STEPS = Math.round(T_END / DT);
const R_BIND = 12.0;
const E_KINK = 2 * Math.SQRT2 / 3;          // energy of one kink at rest

const x = (i) => i * DX;
const SQ2 = Math.SQRT2;

// a kink (rising -1 -> +1) boosted to velocity v>0, at position x0 at t=0
function kinkMoving(xx, x0, v, sign) {
  // sign = +1 for kink, -1 for antikink
  const g = 1 / Math.sqrt(1 - v * v);
  const arg = g * (xx - x0) / SQ2;
  const th = Math.tanh(arg);
  const sech2 = 1 - th * th;
  const phi = sign * th;
  const phit = sign * (-v * g / SQ2) * sech2;   // d/dt of the boosted profile
  return { phi, phit };
}

function init(v, d, eps) {
  const x1 = 0.35 * L;
  const x2 = x1 + d;
  const phi = new Float64Array(N), phiOld = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    const xx = x(i);
    const k = kinkMoving(xx, x1, +v, +1);       // kink moving right
    const a = kinkMoving(xx, x2, -v, -1);       // antikink moving left
    let p = k.phi * a.phi;                       // product form: -1 | +1 | -1
    let pt = k.phit * a.phi + k.phi * a.phit;
    // the 1e-9 perturbation for the reproducibility test
    p += eps * Math.sin(7.0 * xx);
    phi[i] = p;
    phiOld[i] = p - DT * pt;                     // phi at t = -DT
  }
  phi[0] = -1; phi[N - 1] = -1;
  phiOld[0] = -1; phiOld[N - 1] = -1;
  return { phi, phiOld };
}

function energy(phi, phiOld) {
  let e = 0;
  for (let i = 1; i < N - 1; i++) {
    const phit = (phi[i] - phiOld[i]) / DT;
    const phix = (phi[i + 1] - phi[i - 1]) / (2 * DX);
    const V = 0.25 * (phi[i] * phi[i] - 1) ** 2;
    e += (0.5 * phit * phit + 0.5 * phix * phix + V) * DX;
  }
  return e;
}

function separation(phi) {
  // the two steepest-gradient positions; returns their distance (or Infinity if fewer than two)
  const g = [];
  for (let i = 2; i < N - 2; i++) g.push({ i, s: Math.abs(phi[i + 1] - phi[i - 1]) });
  g.sort((p, q) => q.s - p.s);
  if (g.length < 2) return Infinity;
  const a = g[0].i, b = g[1].i;
  if (Math.abs(a - b) < 4) return 0;
  return Math.abs(a - b) * DX;
}

// COUNT THE KINKS PROPERLY.  My first classifier used Q = (phi[N-3]-phi[2])/2, but BOTH boundaries are
// pinned to -1 in this configuration, so Q is identically zero and the "pass" branch could never fire --
// which is exactly why the first table reported pass = 0.  A kink is a local maximum of |phi_x| above a
// threshold; counting those peaks distinguishes "two kinks came out" from "one lump" from "vacuum".
// THIRD ATTEMPT AT THE CLASSIFIER (attempts 1 and 2 were both wrong: attempt 1 used a boundary Q that is
// identically zero because both ends are pinned, attempt 2 counted |phi_x| peaks and missed merged bions).
// This one clusters the ENERGY DENSITY, which is the standard way: a kink is a localised energy peak, a
// merged bion is one broader peak, the vacuum has none.  Threshold and minimum width are fixed here.
function lumps(phi, phiOld) {
  const eps = new Float64Array(N);
  for (let i = 1; i < N - 1; i++) {
    const phit = (phi[i] - phiOld[i]) / DT;
    const phix = (phi[i + 1] - phi[i - 1]) / (2 * DX);
    eps[i] = 0.5 * phit * phit + 0.5 * phix * phix + 0.25 * (phi[i] * phi[i] - 1) ** 2;
  }
  let maxE = 0;
  for (let i = 1; i < N - 1; i++) maxE = Math.max(maxE, eps[i]);
  if (maxE < 0.05) return { n: 0, centers: [], maxE, frac: 0 };
  const thr = 0.12 * maxE;                    // fixed in advance
  const out = [];
  let i = 1;
  while (i < N - 1) {
    if (eps[i] > thr) {
      let j = i, sum = 0, wsum = 0, e = 0;
      while (j < N - 1 && eps[j] > thr) { sum += eps[j] * j; wsum += eps[j]; e += eps[j]; j++; }
      if (j - i >= 3) out.push({ c: (sum / wsum) * DX, e: e * DX, w: (j - i) * DX });
      i = j;
    } else i++;
  }
  let tot = 0;
  for (let k = 1; k < N - 1; k++) tot += eps[k] * DX;
  return { n: out.length, centers: out.map((o) => o.c), widths: out.map((o) => o.w), maxE, frac: tot };
}

function run(v, d, eps, record) {
  const { phi, phiOld } = init(v, d, eps);
  let p = phi, o = phiOld;
  let n = new Float64Array(N);
  const sTrace = [];
  for (let s = 0; s < STEPS; s++) {
    for (let i = 1; i < N - 1; i++) {
      const lap = (p[i + 1] - 2 * p[i] + p[i - 1]) / (DX * DX);
      const acc = lap - (p[i] ** 3 - p[i]);
      n[i] = 2 * p[i] - o[i] + DT * DT * acc;
    }
    n[0] = -1; n[N - 1] = -1;
    const t = o; o = p; p = n; n = t;
    if (record && s % 100 === 0) sTrace.push({ t: (s + 1) * DT, s: separation(p) });
  }
  const Q = (p[N - 3] - p[2]) / 2;
  const E = energy(p, o);
  const sEnd = separation(p);
  const lm = lumps(p, o);
  const nK = lm.n;
  const radiated = E - 2 * E_KINK;
  let cls;
  if (nK === 0) cls = 'vacuum';
  else if (nK === 1) cls = 'bind';                    // one lump: a bion (oscillating bound state)
  else {
    // two or more lumps: are they moving apart?  compare the spread at the end with a snapshot earlier
    cls = sEnd >= R_BIND ? 'pass' : 'bind';
  }
  return { Q, E, sEnd, nK, nLumps: nK, radiated, cls, sTrace, lumpC: lm.centers };
}

console.log('============================================================');
console.log('CHEMISTRY PHASE 0 -- phi^4 (Klein-Gordon) kink-antikink scattering');
console.log('criteria and thresholds were fixed in PREREG_CHEMISTRY.md BEFORE this run');
console.log('============================================================');
console.log(`   N=${N} DX=${DX} DT=${DT} T=${T_END} (${STEPS} steps) R_bind=${R_BIND}`);
console.log('');

const VS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40];
const DS = [6, 8, 10, 12, 16];

console.log('============================================================');
console.log('PART A: THE OUTCOME TABLE (clean start, eps = 0)');
console.log('============================================================');
console.log('   v \\ d     6        8        10       12       16');
const table = [];
for (const v of VS) {
  const row = [];
  for (const d of DS) {
    const r = run(v, d, 0, false);
    row.push(r);
  }
  table.push(row);
  console.log(`   ${v.toFixed(2)}    ` + row.map((r) => r.cls.padEnd(8)).join(' '));
}
console.log('');
const counts = { pass: 0, bind: 0, vacuum: 0 };
for (const row of table) for (const r of row) counts[r.cls]++;
const total = VS.length * DS.length;
console.log(`   totals over ${total} cells: pass ${counts.pass}, bind ${counts.bind}, vacuum ${counts.vacuum}`);
console.log('');

console.log('============================================================');
console.log('PART B: H1 -- REPRODUCIBILITY UNDER A 1e-9 PERTURBATION');
console.log('============================================================');
{
  let flips = 0, cells = 0;
  const flipList = [];
  for (let iv = 0; iv < VS.length; iv++) {
    for (let id = 0; id < DS.length; id++) {
      const v = VS[iv], d = DS[id];
      const base = table[iv][id].cls;
      const pert = run(v, d, 1e-9, false).cls;
      cells++;
      if (base !== pert) { flips++; flipList.push(`v=${v} d=${d}: ${base}->${pert}`); }
    }
  }
  const rate = flips / cells;
  console.log(`   cells: ${cells}, outcome flips under a 1e-9 perturbation: ${flips} (${(100 * rate).toFixed(1)}%)`);
  for (const f of flipList.slice(0, 8)) console.log('     ' + f);
  console.log('');
  console.log(`   H1 verdict (pre-registered): rate <= 10% -> function of (v,d); rate > 30% -> CHAOTIC`);
  console.log(`   => ${rate <= 0.10 ? 'PASSES: the table is a function of (v,d)' : rate > 0.30 ? 'FAILS: the scattering is chaotic' : 'GREY: intermediate (10-30%)'}`);
}
console.log('');

console.log('============================================================');
console.log('PART C: H2 -- DO BOUND PRODUCTS EXIST?');
console.log('============================================================');
{
  const frac = counts.bind / total;
  console.log(`   bind fraction over the table: ${(100 * frac).toFixed(1)}%`);
  console.log(`   H2 verdict (pre-registered): >= 20% exists | < 5% stop | in between grey`);
  console.log(`   => ${frac >= 0.20 ? 'PASSES: bound products exist' : frac < 0.05 ? 'FAILS: no bound products' : 'GREY: partial'}`);
}
console.log('');

console.log('============================================================');
console.log('PART D: DIAGNOSTICS (why the outcome is what it is)');
console.log('============================================================');
{
  console.log('   v      d     nK   E_total    E - 2*E_kink   separation   class');
  for (const v of [0.05, 0.20, 0.40]) {
    for (const d of [6, 12]) {
      const r = run(v, d, 0, false);
      console.log(`   ${v.toFixed(2)}   ${String(d).padStart(2)}    ${String(r.nK).padStart(2)}   ${r.E.toFixed(4).padStart(8)}   ${r.radiated.toFixed(4).padStart(10)}   ${String(r.sEnd.toFixed(2)).padStart(8)}    ${r.cls}`);
    }
  }
  console.log('');
  console.log('   readings: nK = number of kink peaks at the end (2 = both survived, 1 = merged, 0 = vacuum);');
  console.log('   E_total far above 2*E_kink means a lot of the energy went into radiation;');
  console.log('   separation near zero at the end means the pair is sitting on top of each other.');
}
