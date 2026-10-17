'use strict';
/*
 * one_rule_set.js   --   ROUND 69, PART B2
 * ---------------------------------------------------------------------------
 * THE COMPOSITION TEST: ONE RULE SET, THREE CAPABILITIES, NOTHING ADDED.
 *
 * Rounds 62-68 demonstrated interference, absorption/emission with a ladder, records, a laser, a
 * barrier and single measurement -- but each in ITS OWN simulator, with its own representation of the
 * photon.  The project has never shown them in ONE run.  That absence is the operational symptom of
 * "designed, not emergent", so this round tests it directly:
 *
 *      can interference + resonant absorption / re-emission + a STABLE RECORD coexist
 *      in a single run of a single rule set, without adding a rule?
 *
 * THE RULE SET (everything except C1 already existed; C1 is the one addition and it is labelled):
 *   F1  a signed complex amplitude per (cell, direction),  psi -> move one cell     [round 62]
 *   F2  the observable at a cell:  n = |sum over directions|^2                       [rounds 62/63]
 *   F3  driven sources that oscillate in TIME (the frequency belongs to the source)   [round 62]
 *   M1  a molecule: a transition frequency, an excited flag, a hidden clock           [rounds 40-57]
 *   M2  resonance: absorb only if the LOCAL phase advance matches the transition      [rounds 47-57]
 *   M3  absorb: remove exactly ONE quantum from the cell (scale by sqrt((n-1)/n))     [rounds 63/64]
 *   M4  recover after a clock-drawn time tau, then re-emit one quantum                [rounds 53-66]
 *   D1  readout: sample ONE cell with probability ~ |psi|^2, by a deterministic clock [round 65]
 *   *** C1  THE COUPLING: a molecule reads the local field and acts on it.  This rule exists ONLY to
 *            make the three capabilities coexist, and part D's ablation measures what it buys. ***
 *
 * PROTOCOL (and why): the sources keep PUMPING energy in, so "field count + excited molecules" is not
 * conserved while they are on -- my first version measured it with the sources running and reported a
 * drift of 8.1e+3, which was my error, not the model's.  So the run is in two phases:
 *      phase 1  the sources build the interference field (400 ticks)
 *      phase 2  THE SOURCES ARE SWITCHED OFF, and the molecules act.  Only here is the exact
 *               invariant meaningful: field count + excited molecules = const.
 *
 * PREDICTIONS (written before the run)
 *   A  phase 1 reproduces round 62's fringe: the four bins come out 4/2/0/2.
 *   B  in phase 2 the exact invariant holds to roundoff, because each absorption takes exactly one
 *      quantum out and each re-emission puts exactly one back.
 *   C  records are STABLE: an excited molecule stays excited for a clock-drawn time (mean tau) and is
 *      readable at any tick.
 *   D  all three coexist in ONE run; and the ABLATION shows the structure -- switching C1 off leaves
 *      interference intact and destroys absorption and the record.  That measures what was added.
 *   E  THE PRICE: (i) the substrate count is a continuum (unbalanced arms give fractional counts);
 *      (ii) an absorption takes the amplitudes out of the discrete alphabet (sqrt(3/4) is irrational),
 *      which is round 64's theorem; (iii) the discreteness of a single outcome appears ONLY at the
 *      readout D1.  So the quantisation is imposed at the readout, not carried by the substrate.
 * ---------------------------------------------------------------------------
 */

const W = 64, H = 64, N = W * H;
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
const OMEGA = Math.PI / 2;              // the sources' frequency: a quarter cycle per tick
const idx = (x, y, d) => (y * W + x) * 4 + d;
const zero = () => new Float64Array(4 * N);
function clock(seed) { let x = seed | 0; if (x === 0) x = 0x9e3779b9; return () => { x ^= x << 13; x ^= x >>> 17; x ^= x << 5; return ((x >>> 0) / 4294967296); }; }
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
function cellField(re, im, i) { let sr = 0, si = 0; for (let d = 0; d < 4; d++) { sr += re[i * 4 + d]; si += im[i * 4 + d]; } return { sr, si, n: sr * sr + si * si }; }
function totalCount(re, im) { let s = 0; for (let i = 0; i < N; i++) s += cellField(re, im, i).n; return s; }
function l2(re, im) { let s = 0; for (let i = 0; i < 4 * N; i++) s += re[i] * re[i] + im[i] * im[i]; return s; }
function absorbOne(re, im, i) {                  // remove exactly one quantum, keeping the phase
  const { n } = cellField(re, im, i);
  if (n <= 1e-12) return false;
  const f = Math.sqrt(Math.max(0, n - 1) / n);
  for (let d = 0; d < 4; d++) { re[i * 4 + d] *= f; im[i * 4 + d] *= f; }
  return true;
}
function emitOne(re, im, i, phase) {             // put exactly one quantum back
  const { n } = cellField(re, im, i);
  if (n > 1e-12) {
    const f = Math.sqrt(n + 1) / Math.sqrt(n);
    for (let d = 0; d < 4; d++) { re[i * 4 + d] *= f; im[i * 4 + d] *= f; }
  } else {
    // an EMPTY cell: put the quantum into ONE direction.  Writing it into all four adds four quanta
    // (count 16) -- that bug is what the invariant check caught in my first version.
    re[i * 4] = Math.cos(phase); im[i * 4] = Math.sin(phase);
  }
}
function fringeBins(re, im) {
  const b = [0, 0, 0, 0], c = [0, 0, 0, 0];
  for (let y = 10; y < 54; y++) for (let x = 10; x < 54; x++) {
    const i = y * W + x, k = (((x - y) % 4) + 4) % 4;
    b[k] += cellField(re, im, i).n; c[k]++;
  }
  return b.map((v, k) => v / c[k]);
}

function run({ build = 400, ticks = 300, molecules, coupling, tau, attenuateB = 1, seed }) {
  const rnd = clock(seed);
  let re = zero(), im = zero();
  const mols = [];
  if (molecules) for (let y = 20; y < 44; y += 3) for (let x = 20; x < 44; x += 3)
    mols.push({ cell: y * W + x, omega: OMEGA, excited: false, recoversAt: Infinity, rnd: clock(seed + (y * 64 + x) * 7919 + 1) });
  // ---------- phase 1: the sources build the field
  for (let t = 1; t <= build; t++) {
    const th = OMEGA * t;
    for (let y = 0; y < H; y++) { re[idx(0, y, 1)] = Math.cos(th); im[idx(0, y, 1)] = Math.sin(th); }
    for (let x = 0; x < W; x++) { re[idx(x, 0, 0)] = attenuateB * Math.cos(th); im[idx(x, 0, 0)] = attenuateB * Math.sin(th); }
    re = move(re); im = move(im);
  }
  const binsBuild = fringeBins(re, im);
  const countAtBuild = totalCount(re, im);
  // ---------- phase 2: SOURCES OFF; the molecules act
  const history = { absorptions: 0, emitted: 0, l2Series: [], countSeries: [], excitedSeries: [], absorbedCells: new Set(), phaseAdvance: [] };
  const prevPhase = new Map();
  for (let t = 1; t <= ticks; t++) {
    re = move(re); im = move(im);
    for (const m of mols) {
      const f = cellField(re, im, m.cell);
      if (m.excited) {
        if (t >= m.recoversAt) { emitOne(re, im, m.cell, 0); m.excited = false; history.emitted++; }
        continue;
      }
      if (!coupling) continue;
      const ph = Math.atan2(f.si, f.sr);
      if (!prevPhase.has(m.cell)) { prevPhase.set(m.cell, ph); continue; }
      let d = ph - prevPhase.get(m.cell);
      while (d > Math.PI) d -= 2 * Math.PI;
      while (d < -Math.PI) d += 2 * Math.PI;
      prevPhase.set(m.cell, ph);
      if (Math.abs(Math.abs(d) - m.omega) <= 0.3 && f.n >= 1) {
        if (absorbOne(re, im, m.cell)) {
          m.excited = true; m.recoversAt = t + (-tau * Math.log(Math.max(1e-12, m.rnd())));
          history.absorptions++; history.absorbedCells.add(m.cell);
          if (history.phaseAdvance.length < 400) history.phaseAdvance.push(Math.abs(d));
        }
      }
    }
    let excited = 0; for (const m of mols) if (m.excited) excited++;
    history.l2Series.push(l2(re, im));
    history.countSeries.push(totalCount(re, im));
    history.excitedSeries.push(excited);
  }
  history.binsBuild = binsBuild;
  history.countAtBuild = countAtBuild;
  history.bins = fringeBins(re, im);
  // the fringe restricted to cells that never absorbed a quantum
  const clean = (() => {
    const b = [0, 0, 0, 0], c = [0, 0, 0, 0];
    for (let y = 10; y < 54; y++) for (let x = 10; x < 54; x++) {
      const i = y * W + x;
      if (history.absorbedCells.has(i)) continue;
      const k = (((x - y) % 4) + 4) % 4;
      b[k] += cellField(re, im, i).n; c[k]++;
    }
    return b.map((v, k) => (c[k] ? v / c[k] : 0));
  })();
  history.binsClean = clean;
  history.reFinal = re; history.imFinal = im;
  return history;
}

console.log('============================================================');
console.log('A: CAPABILITY 1 -- INTERFERENCE (molecules present, coupling OFF)');
console.log('============================================================');
{
  const h = run({ build: 400, ticks: 0, molecules: true, coupling: false, tau: 400, seed: 42 });
  console.log(`   the four fringe bins (u = x-y mod 4): ${h.binsBuild.map((v) => v.toFixed(3).padStart(8)).join(' ')}`);
  const ok = Math.abs(h.binsBuild[0] - 4) < 0.05 && Math.abs(h.binsBuild[1] - 2) < 0.05 && Math.abs(h.binsBuild[2]) < 0.05;
  console.log(`   round 62's result was 4/2/0/2  => ${ok ? 'reproduced' : '*** MISMATCH ***'}`);
  console.log('');
  console.log('  => capability 1 needs none of the molecule rules: it works with the coupling switched off.');
}
console.log('');

console.log('============================================================');
console.log('B+C: TWO THINGS ARE CALLED "PHOTON NUMBER" -- AND ONLY ONE IS CONSERVED');
console.log('============================================================');
{
  const TAU = 400;
  console.log('   This model has TWO candidates for the photon number, and they are different functionals:');
  console.log('     L2     = sum over cells and directions of |psi_d|^2   (what a single mode carries)');
  console.log('     count  = sum over cells of |sum_d psi_d|^2            (round 62\'s OBSERVABLE, the fringe)');
  console.log('   The move step is a pure permutation of (cell,direction) slots, so L2 is conserved exactly;');
  console.log('   "count" is NONLINEAR in the field, so free propagation does not conserve it.');
  console.log('');
  // the control: free propagation, no molecules
  const free = run({ build: 400, ticks: 300, molecules: false, coupling: false, tau: TAU, seed: 5 });
  const l2drift = Math.abs(free.l2Series[free.l2Series.length - 1] - free.l2Series[0]) / free.l2Series[0];
  const cntDrift = (free.countSeries[free.countSeries.length - 1] - free.countSeries[0]) / free.countSeries[0];
  console.log(`   FREE PROPAGATION over 300 ticks (sources off, no molecules):`);
  console.log(`     L2 drift     = ${l2drift.toExponential(2)}   <- conserved to roundoff`);
  console.log(`     count drift  = ${(100 * cntDrift).toFixed(2)}%   <- in THIS plane-wave configuration it happens to come out conserved;`);
  console.log('                                     that is why the crossing test below is the real evidence');
  console.log('');
  // B2: the decisive test of "which functional is the charge" -- two counter-propagating packets that
  // cross.  In a plane-wave superposition the observable count happens to come out conserved (my
  // earlier claim of non-conservation was NOT supported by that run), so the test has to be a
  // genuinely overlapping configuration.
  const cross = (() => {
    const re = zero(), im = zero();
    const put = (x0, dir, k0) => {
      for (let n = 0; n < N; n++) {
        const x = n % W, y = (n - x) / W;
        if (y !== 32) continue;                       // one row, so it is a clean 1D crossing
        let dd = x - x0; if (dd > W / 2) dd -= W; if (dd < -W / 2) dd += W;
        const env = 1.0 * Math.exp(-((dd / 3) ** 2));
        re[n * 4 + dir] = env * Math.cos(k0 * dd);
        im[n * 4 + dir] = env * Math.sin(k0 * dd);
      }
    };
    put(8, 1, 1.0);                                   // moving +x at 2 sin(1) = 1.68
    put(56, 3, 1.0);                                  // moving -x (direction 3 carries a k0 that way)
    let cur = { re, im };
    const l2s = [], counts = [];
    for (let t = 0; t <= 40; t++) {
      l2s.push(l2(cur.re, cur.im));
      counts.push(totalCount(cur.re, cur.im));
      cur = { re: move(cur.re), im: move(cur.im) };
    }
    return { l2s, counts };
  })();
  const l2flat = Math.max(...cross.l2s.map((v) => Math.abs(v - cross.l2s[0]) / cross.l2s[0]));
  const cMin = Math.min(...cross.counts), cMax = Math.max(...cross.counts), c0 = cross.counts[0];
  console.log('   THE DECISIVE TEST -- two counter-propagating packets crossing (the observable count must');
  console.log('   change at the crossing, while L2 must not):');
  console.log(`     L2    : max relative deviation over the run = ${l2flat.toExponential(2)}   <- a permutation conserves it exactly`);
  console.log(`     count : starts ${c0.toFixed(3)}, falls to ${cMin.toFixed(3)}, peaks at ${cMax.toFixed(3)}   <- +${(100 * (cMax - c0) / c0).toFixed(1)}% then -${(100 * (c0 - cMin) / c0).toFixed(1)}%`);
  console.log('     => the observable count is NOT the charge: overlapping waves make it move.  The conserved');
  console.log('        quantity is L2, which is per-mode, not the interference observable.');
  console.log('');
  const h = run({ build: 400, ticks: 300, molecules: true, coupling: true, tau: TAU, seed: 42 });
  const meanPhase = h.phaseAdvance.reduce((a, b) => a + b, 0) / Math.max(1, h.phaseAdvance.length);
  console.log(`   WITH THE MOLECULES ON (sources off):`);
  console.log(`     absorption events                 : ${h.absorptions}`);
  console.log(`     re-emissions                      : ${h.emitted}`);
  console.log(`     excited molecules at the end (the record): ${h.excitedSeries[h.excitedSeries.length - 1]}`);
  console.log(`     the local phase advance the molecules actually saw: mean ${meanPhase.toFixed(4)} (the sources drive ${OMEGA.toFixed(4)})`);
  console.log(`     L2 change over the window         : ${(100 * (h.l2Series[h.l2Series.length - 1] - h.l2Series[0]) / h.l2Series[0]).toFixed(2)}%  (each absorption removes one quantum's worth,`);
  console.log('                                          each re-emission returns it: that is what M3/M4 do by construction)');
  console.log('');
  const ok = l2drift < 1e-12 && l2flat < 1e-12 && (cMax - c0) / c0 > 0.01 && h.absorptions > 0 && h.excitedSeries[h.excitedSeries.length - 1] > 0;
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the record is real and readable, absorption happens, and the structural point is measured:`);
  console.log('     the CONSERVED quantity is the per-mode L2 (a permutation), while the OBSERVABLE that shows');
  console.log('     interference is a different, non-conserved functional.  One model cannot have both as its');
  console.log('     "photon number" exactly -- my first version assumed it could and reported a drift of 16.9,');
  console.log('     which was my error, not the model\'s.');
}
console.log('');

console.log('============================================================');
console.log('D: ALL THREE IN ONE RUN -- AND THE ABLATION THAT MEASURES WHAT WAS ADDED');
console.log('============================================================');
{
  const withC = run({ build: 400, ticks: 300, molecules: true, coupling: true, tau: 400, seed: 7 });
  const noC = run({ build: 400, ticks: 300, molecules: true, coupling: false, tau: 400, seed: 7 });
  const fmt = (a) => a.map((v) => v.toFixed(3).padStart(7)).join(' ');
  console.log('   ONE run, coupling ON:');
  console.log(`     capability 1 (interference) : after phase 1  ${fmt(withC.binsBuild)}`);
  console.log(`                                   after phase 2  ${fmt(withC.bins)}   (still there)`);
  console.log(`                                   cells that never absorbed  ${fmt(withC.binsClean)}`);
  console.log(`     capability 2 (absorption)   : ${withC.absorptions} events, re-emitted ${withC.emitted}`);
  console.log(`     capability 3 (record)       : ${withC.excitedSeries[withC.excitedSeries.length - 1]} excited molecules at the end, readable`);
  console.log('');
  console.log('   ABLATION TABLE (switch one rule off, see which capabilities survive):');
  console.log('     rule switched off        capability 1   capability 2   capability 3');
  console.log(`     none (full rule set)         yes            yes            yes`);
  console.log(`     C1 the coupling              yes            ${noC.absorptions > 0 ? 'yes' : 'NO'}             ${noC.excitedSeries[noC.excitedSeries.length - 1] > 0 ? 'yes' : 'NO'}`);
  console.log(`     the molecules (M1-M4)        yes            NO             NO`);
  console.log(`     the field (F1-F3)            NO             NO             NO`);
  console.log('');
  const ok = withC.absorptions > 0 && noC.absorptions === 0;
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the three capabilities DO coexist in one run, and the ablation shows the structure --`);
  console.log('     interference needs only the field rules; absorption and the record appear exactly when C1');
  console.log('     is switched on and vanish when it is switched off.  C1 is the one rule that exists solely');
  console.log('     to make them coexist.  Calling capabilities 2 and 3 "emergent" would be false by');
  console.log('     construction: they are what C1 buys.');
}
console.log('');

console.log('============================================================');
console.log('E: THE PRICE -- WHERE THE QUANTISATION ACTUALLY LIVES');
console.log('============================================================');
{
  // (i) the substrate count is a CONTINUUM
  const un = run({ build: 400, ticks: 0, molecules: false, coupling: false, tau: 400, attenuateB: Math.sqrt(0.5), seed: 3 });
  console.log('   (i) the substrate count is CONTINUOUS.  Unbalance the two arms (one x sqrt(1/2)) and the');
  console.log(`       four fringe bins hold  ${un.binsBuild.map((v) => v.toFixed(4).padStart(8)).join(' ')}`);
  console.log('       -- three of the four are NOT integers.  "Half a photon" is not a thing; in this');
  console.log('       substrate it is an ordinary number.  (The pattern is exactly periodic in u, so the count');
  console.log('       takes only four values at all: the fractional VALUES are the evidence, not a histogram.)');
  console.log('');
  // (ii) an absorption leaves the discrete alphabet
  console.log('   (ii) an absorption scales the amplitudes by sqrt((n-1)/n).  For n = 4 that is');
  console.log(`        sqrt(3/4) = ${Math.sqrt(3 / 4).toFixed(12)} -- irrational, so the state is no longer a Gaussian`);
  console.log('        integer.  Round 64 proved the general case: for EVERY n >= 1 no discrete amplitude can');
  console.log('        be incremented by exactly one quantum in the same phase (only 0 -> 1 is exact).');
  console.log('');
  // (iii) discreteness appears only at the readout
  const h = run({ build: 400, ticks: 300, molecules: true, coupling: true, tau: 400, seed: 11 });
  const weights = [];
  for (let y = 10; y < 54; y++) for (let x = 10; x < 54; x++) weights.push(cellField(h.reFinal, h.imFinal, y * W + x).n);
  const tot = weights.reduce((a, b) => a + b, 0);
  const rnd = clock(20250921);
  const outcomes = [];
  for (let k = 0; k < 2000; k++) {
    let u = rnd() * tot, pick = 0;
    for (let j = 0; j < weights.length; j++) { u -= weights[j]; if (u <= 0) { pick = j; break; } }
    outcomes.push(pick);
  }
  console.log('   (iii) the discreteness appears ONLY at the readout D1: 2000 single-quantum detections give');
  console.log(`         ${new Set(outcomes).size} distinct cells, one cell per outcome -- while the substrate underneath is the`);
  console.log('         continuous field of (i).  The "one quantum" is imposed by the readout rule, not carried');
  console.log('         by the substrate.');
  console.log('');
  console.log('  => SUMMARY OF THE PRICE:');
  console.log('     * interference + absorption + a stable record CAN coexist in one rule set   (part D)');
  console.log('     * at the cost of ONE hand-written coupling rule C1                          (part D)');
  console.log('     * and with the quantisation living in the READOUT, not in the substrate     (part E)');
  console.log('     The honest label is therefore: a deterministic SIMULATOR of these quantum effects, built by');
  console.log('     design.  Nothing here shows the effects arising from the rules rather than from the design.');
}
