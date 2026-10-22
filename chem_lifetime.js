'use strict';
/*
 * chem_lifetime.js  --  CHEMISTRY LINE, PHASE 1-A: IS THE "BIND" REAL?  (see PREREG_CHEMISTRY_P1.md)
 * ---------------------------------------------------------------------------
 * Phase 0 found that 80% of the (v,d) cells end with ONE energy lump at T = 300.  But "still one lump"
 * is not "bound" -- it may only be "radiating slowly".  This script re-runs the 28 bind cells to TEN
 * times that length and asks whether the lump is still there and still holds its energy.
 *
 * The lesson being applied is from the diffusion round (D5): a snapshot is not an equilibrium.
 *
 * PRE-REGISTERED CRITERIA (fixed before the run)
 *   genuine bion   : at T = 3000, one lump and in-lump energy fraction f >= 0.80
 *   slow decay     : at T = 3000, zero lumps, or f < 0.50
 *   grey           : anything else
 *   two representative cells run to T = 9000 (three times longer) as a long-time check
 *   MY BETS: low v (0.05-0.15) survive at >= 80%; high v (0.25-0.40) at <= 50%; overall 50-75%
 * ---------------------------------------------------------------------------
 */

const N = 1024, DX = 0.25, DT = 0.05;
const L = (N - 1) * DX;
const SQ2 = Math.SQRT2;
const E_KINK = 2 * Math.SQRT2 / 3;

const x = (i) => i * DX;

function kinkMoving(xx, x0, v, sign) {
  const g = 1 / Math.sqrt(1 - v * v);
  const arg = g * (xx - x0) / SQ2;
  const th = Math.tanh(arg);
  const sech2 = 1 - th * th;
  return { phi: sign * th, phit: sign * (-v * g / SQ2) * sech2 };
}

function init(v, d) {
  const x1 = 0.35 * L, x2 = x1 + d;
  const phi = new Float64Array(N), phiOld = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    const xx = x(i);
    const k = kinkMoving(xx, x1, +v, +1), a = kinkMoving(xx, x2, -v, -1);
    const p = k.phi * a.phi, pt = k.phit * a.phi + k.phi * a.phit;
    phi[i] = p; phiOld[i] = p - DT * pt;
  }
  phi[0] = -1; phi[N - 1] = -1; phiOld[0] = -1; phiOld[N - 1] = -1;
  return { phi, phiOld };
}

function lumps(phi, phiOld) {
  const eps = new Float64Array(N);
  for (let i = 1; i < N - 1; i++) {
    const phit = (phi[i] - phiOld[i]) / DT;
    const phix = (phi[i + 1] - phi[i - 1]) / (2 * DX);
    eps[i] = 0.5 * phit * phit + 0.5 * phix * phix + 0.25 * (phi[i] * phi[i] - 1) ** 2;
  }
  let maxE = 0, tot = 0;
  for (let i = 1; i < N - 1; i++) { maxE = Math.max(maxE, eps[i]); tot += eps[i] * DX; }
  if (maxE < 0.05) return { n: 0, inLump: 0, tot, maxE };
  const thr = 0.12 * maxE;
  let n = 0, inLump = 0;
  let i = 1;
  while (i < N - 1) {
    if (eps[i] > thr) {
      let j = i, e = 0;
      while (j < N - 1 && eps[j] > thr) { e += eps[j] * DX; j++; }
      if (j - i >= 3) { n++; inLump += e; }
      i = j;
    } else i++;
  }
  return { n, inLump, tot, maxE };
}

function runTo(v, d, T) {
  const st = init(v, d);
  let p = st.phi, o = st.phiOld;
  let n = new Float64Array(N);
  const steps = Math.round(T / DT);
  const trace = [];
  for (let s = 0; s < steps; s++) {
    for (let i = 1; i < N - 1; i++) {
      const lap = (p[i + 1] - 2 * p[i] + p[i - 1]) / (DX * DX);
      n[i] = 2 * p[i] - o[i] + DT * DT * (lap - (p[i] ** 3 - p[i]));
    }
    n[0] = -1; n[N - 1] = -1;
    const t = o; o = p; p = n; n = t;
    if ((s + 1) % 200 === 0) {
      const lm = lumps(p, o);
      trace.push({ t: (s + 1) * DT, n: lm.n, f: lm.tot > 0 ? lm.inLump / lm.tot : 0, tot: lm.tot, maxE: lm.maxE });
    }
  }
  const last = trace[trace.length - 1];
  let cls;
  if (last.n === 0 || last.f < 0.50) cls = 'slow-decay';
  else if (last.n === 1 && last.f >= 0.80) cls = 'bion';
  else cls = 'grey';
  return { last, trace, cls };
}

// the 28 cells Phase 0 called "bind" (recomputed by hand from chem_phase0_report.txt)
const BIND = [];
for (const v of [0.05, 0.10, 0.15]) for (const d of [6, 8, 10, 12, 16]) BIND.push([v, d]);
for (const d of [8, 10, 12, 16]) BIND.push([0.20, d]);
for (const d of [6, 8, 10, 12, 16]) BIND.push([0.25, d]);
for (const d of [6, 10]) BIND.push([0.30, d]);
for (const d of [8, 10]) BIND.push([0.40, d]);

console.log('============================================================');
console.log('CHEMISTRY PHASE 1-A -- LIFETIME TEST OF THE 28 "BIND" CELLS');
console.log('criteria fixed in PREREG_CHEMISTRY_P1.md before this run');
console.log('============================================================');
console.log(`   cells: ${BIND.length}; T = 3000 (10x Phase 0); long check T = 9000 on two cells`);
console.log('');

console.log('============================================================');
console.log('PART A: THE LIFETIME TABLE');
console.log('============================================================');
console.log('   v      d     lumps(T=3000)   in-lump fraction   total E     verdict');
const res = [];
for (const [v, d] of BIND) {
  const r = runTo(v, d, 3000);
  res.push({ v, d, ...r });
  console.log(`   ${v.toFixed(2)}   ${String(d).padStart(2)}        ${String(r.last.n).padStart(2)}              ${r.last.f.toFixed(3).padStart(6)}         ${r.last.tot.toFixed(4)}    ${r.cls}`);
}
console.log('');
const cnt = { bion: 0, 'slow-decay': 0, grey: 0 };
for (const r of res) cnt[r.cls]++;
console.log(`   verdicts: bion ${cnt.bion}, slow-decay ${cnt['slow-decay']}, grey ${cnt.grey} (of ${res.length})`);
const survive = cnt.bion / res.length;
console.log('');

console.log('============================================================');
console.log('PART B: MY BETS VS THE RESULT');
console.log('============================================================');
{
  const lowV = res.filter((r) => r.v <= 0.15);
  const highV = res.filter((r) => r.v >= 0.25);
  const lowRate = lowV.filter((r) => r.cls === 'bion').length / lowV.length;
  const highRate = highV.filter((r) => r.cls === 'bion').length / highV.length;
  console.log(`   low  v (0.05-0.15): bion ${(100 * lowRate).toFixed(0)}%   [my bet: >= 80%]   => ${lowRate >= 0.8 ? 'BET WON' : 'BET LOST'}`);
  console.log(`   high v (0.25-0.40): bion ${(100 * highRate).toFixed(0)}%   [my bet: <= 50%]   => ${highRate <= 0.5 ? 'BET WON' : 'BET LOST'}`);
  console.log(`   overall           : ${(100 * survive).toFixed(0)}%       [my bet: 50-75%]  => ${survive >= 0.5 && survive <= 0.75 ? 'BET WON' : 'BET LOST'}`);
  console.log('');
  console.log(`   Phase 0 said 80% bound; after a 10x longer run it is ${(100 * survive).toFixed(0)}%.`);
  console.log(`   => ${survive >= 0.5 ? 'PASSES the gate: stable species exist, Phase 1-B can proceed' : 'FAILS the gate: no stable species, the chemistry line stops here'}`);
}
console.log('');

console.log('============================================================');
console.log('PART C: LONG-TIME CHECK (T = 9000, three times longer)');
console.log('============================================================');
for (const [v, d] of [[0.05, 6], [0.25, 10]]) {
  const r = runTo(v, d, 9000);
  const f0 = r.trace[Math.floor(r.trace.length / 3)];
  console.log(`   v=${v} d=${d}: at t=${f0.t}: n=${f0.n}, f=${f0.f.toFixed(3)} | at t=${r.last.t}: n=${r.last.n}, f=${r.last.f.toFixed(3)} => ${r.cls}`);
}
console.log('');
console.log('   this is the check that a snapshot cannot give: does the lump keep shrinking after T = 3000?');
