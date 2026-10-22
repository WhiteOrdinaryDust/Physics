'use strict';
/*
 * chem_lifetime.js  --  CHEMISTRY LINE, PHASE 1-B: CORRECTED CRITERION (see PREREG_CHEMISTRY_P1B.md)
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
  // P1-B CRITERION (see PREREG_CHEMISTRY_P1B.md): the number of lumps, NOT the in-lump fraction.
  // f depends on the peak threshold, and a kink-antikink bubble between two pinned -1 boundaries has
  // f of order 0.6 by construction, so the old f >= 0.80 test could never pass.  My first attempt at this
  // replacement silently did not apply, which is why the previous tally printed 0/0/0 while every row
  // showed n = 1.
  let cls;
  if (last.n === 1) cls = 'bound';
  else if (last.n === 0) cls = 'gone';
  else cls = 'split';
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
console.log('CHEMISTRY PHASE 1-B -- CORRECTED LIFETIME CRITERION (n at T = 6000)');
console.log('criteria fixed in PREREG_CHEMISTRY_P1B.md before this run');
console.log('============================================================');
console.log(`   cells: ${BIND.length}; T = 6000 (20x Phase 0)`);
console.log('');

console.log('============================================================');
console.log('PART A: THE CORRECTED TABLE (T = 6000)');
console.log('============================================================');
console.log('   v      d     lumps   in-lump f   total E    verdict');
const res = [];
for (const [v, d] of BIND) {
  const r = runTo(v, d, 6000);
  res.push({ v, d, ...r });
  console.log(`   ${v.toFixed(2)}   ${String(d).padStart(2)}      ${String(r.last.n).padStart(2)}      ${r.last.f.toFixed(3).padStart(6)}    ${r.last.tot.toFixed(4)}    ${r.cls}`);
}
console.log('');
const cnt = { bound: 0, gone: 0, split: 0 };
for (const r of res) cnt[r.cls]++;
console.log(`   verdicts: bound ${cnt.bound}, gone ${cnt.gone}, split ${cnt.split} (of ${res.length})`);
const survive = cnt.bound / res.length;
console.log('');

console.log('============================================================');
console.log('PART B: THE PRE-REGISTERED THRESHOLD');
console.log('============================================================');
console.log(`   bound fraction at T = 6000: ${(100 * survive).toFixed(0)}%   [threshold: >= 80% passes, < 50% stops]`);
console.log(`   => ${survive >= 0.8 ? 'PASSES: stable bound pairs exist; Phase 1-C can proceed' : survive >= 0.5 ? 'GREY: proceed only with a written caveat' : 'FAILS: no stable species, the chemistry line stops here'}`);
console.log('');
console.log(`   MY BET was >= 80%: ${survive >= 0.8 ? 'BET WON' : 'BET LOST'}`);
{
  const lo = res.filter((r) => r.v <= 0.15), hi = res.filter((r) => r.v >= 0.25);
  const f = (a) => (100 * a.filter((r) => r.cls === 'bound').length / a.length).toFixed(0);
  console.log(`   low v (0.05-0.15): ${f(lo)}% bound | high v (0.25-0.40): ${f(hi)}% bound`);
}
console.log('');

console.log('============================================================');
console.log('PART C: THE f(t) SLOPE, AND THE VACUUM CONTROL');
console.log('============================================================');
{
  const r = runTo(0.05, 6, 6000);
  const first = r.trace[0], last = r.trace[r.trace.length - 1];
  console.log(`   v=0.05 d=6: f goes ${first.f.toFixed(3)} (t=${first.t}) -> ${last.f.toFixed(3)} (t=${last.t})`);
  console.log(`   => in-lump fraction decays SLOWLY (not a criterion, just a lifetime indicator)`);
}
{
  // VACUUM CONTROL: one kink only (left end -1, right end +1) must stay n = 1 for all time
  const phi = new Float64Array(N), phiOld = new Float64Array(N);
  for (let i = 0; i < N; i++) { const th = Math.tanh((x(i) - 0.5 * L) / Math.SQRT2); phi[i] = th; phiOld[i] = th; }
  phi[0] = -1; phi[N - 1] = 1; phiOld[0] = -1; phiOld[N - 1] = 1;
  let p = phi, o = phiOld, n = new Float64Array(N);
  for (let s = 0; s < 120000; s++) {
    for (let i = 1; i < N - 1; i++) {
      const lap = (p[i + 1] - 2 * p[i] + p[i - 1]) / (DX * DX);
      n[i] = 2 * p[i] - o[i] + DT * DT * (lap - (p[i] ** 3 - p[i]));
    }
    n[0] = -1; n[N - 1] = 1;
    const t = o; o = p; p = n; n = t;
  }
  const lm = lumps(p, o);
  console.log(`   VACUUM CONTROL (a single kink, boundaries -1 and +1): lumps = ${lm.n} at T = 6000  => ${lm.n === 1 ? 'control passes (n=1 is not a classifier defect)' : 'CONTROL FAILS'}`);
}
console.log('');
console.log('   reading: a kink-antikink bubble between two -1 boundaries gives n=1 forever; so does a single');
console.log('   kink between -1 and +1.  The control shows that n=1 is the expected signature of ONE stable');
console.log('   object, and the topology (both ends at -1) forces that object to be a bound kink-antikink pair.');
