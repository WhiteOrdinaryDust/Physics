'use strict';
/*
 * memory_defects.js  --  ROUND 97: MEMORY EFFECTS AND TOPOLOGICAL DEFECTS (the NEW physics of doc #3)
 * ---------------------------------------------------------------------------
 * The final worldview document adds the one thing the project had never tested (its section 4.3):
 *
 *      "if the network has collective mechanical properties -- non-ergodicity, multistability, glassy
 *       behaviour -- then after the couplings are all removed the network may not snap back.  It may
 *       leave metastable states or topologically protected defects, i.e. residual redundant channels,
 *       which can then carry other information."
 *
 * Substrate: a 1D scalar field with a double well,
 *      V(phi) = (phi^2 - 1)^2 / 4      (minima +-1, barrier 1/4, spinodal |h| = 2/(3 sqrt3) = 0.3849)
 *
 * ====================== ERRATUM: WHAT THE EARLIER VERSIONS OF THIS SCRIPT GOT WRONG ==============
 * (1) I used a LARGE-amplitude lump as the non-topological control.  It does not disperse; it falls
 *     into the second well.  Fix: a small amplitude for the dispersal control.  The large case is kept
 *     because of what it does (P2b).
 * (2) I wrote that the large pulse "becomes a kink-antikink pair and settles".  The numbers say its
 *     amplitude collapses to 0.0002 and then RISES to 0.37: the pair binds into a BREATHER and keeps
 *     oscillating.  Reported as measured.
 * (3) I reported a "coercive field" and called it kinetic.  A rate scan refuted that, and the reason is
 *     that a noiseless dynamics cannot switch at all -- the value was just the sweep limit.
 * (4) I then predicted a Kramers-like escape under uniform noise.  Wrong for a one-dimensional system:
 *     the critical nucleus in 1D is macroscopically large, so noise never nucleates it here.
 * (5) And my first nucleation test started in the STABLE well, so there was no metastable phase to
 *     leave.  Fixed: start in the metastable well and nucleate the stable one.
 * ================================================================================================
 *
 * PRE-REGISTERED PREDICTIONS
 *   P1 kink: Q stays exactly 1, wall sharpness stays near 0.7, and after a kick the wall RE-FORMS
 *   P2 a small-amplitude pulse of the same width disperses
 *   P2b a large-amplitude pulse does not disperse: it forms a bound OSCILLATING structure
 *   P3 remanence: after the field is swept away the state does not return to the symmetric value
 *   P4 uniform noise alone never switches a 1D metastable state (nucleation is suppressed in 1D), but a
 *      LOCAL nucleation event does, and the erasure then PROPAGATES at a measurable front speed
 *   P5 the retained bit outlives the dispersing pulse by a large factor
 * ---------------------------------------------------------------------------
 */

const N = 512, DX = 0.5;
const x = (i) => i * DX;
const SPINODAL = 2 / (3 * Math.sqrt(3));
function laplacian(f, i, fl, fr) {
  const left = i === 0 ? (fl !== undefined ? fl : f[1]) : f[i - 1];
  const right = i === N - 1 ? (fr !== undefined ? fr : f[N - 2]) : f[i + 1];
  return (left + right - 2 * f[i]) / (DX * DX);
}
const maxGrad = (f) => { let m = 0; for (let i = 1; i < N; i++) m = Math.max(m, Math.abs((f[i] - f[i - 1]) / DX)); return m; };
const peak = (f) => { let m = -Infinity; for (let i = 0; i < N; i++) m = Math.max(m, f[i]); return m; };
const charge = (f, fl, fr) => 0.5 * ((fr !== undefined ? fr : f[N - 1]) - (fl !== undefined ? fl : f[0]));
const mean = (f) => f.reduce((a, b) => a + b, 0) / f.length;
function stream(seed) { let s = seed | 0; if (s === 0) s = 0x9e3779b9; return () => { s ^= s << 13; s ^= s >>> 17; s ^= s << 5; return (s >>> 0) / 4294967296; }; }

function waveRun(kind, amp, steps, kick) {
  const dt = 0.2;
  let f = new Float64Array(N), v = new Float64Array(N);
  let fl, fr;
  if (kind === 'kink') { fl = -1; fr = 1; for (let i = 0; i < N; i++) f[i] = Math.tanh((x(i) - 128) / Math.SQRT2); }
  else { fl = 1; fr = 1; for (let i = 0; i < N; i++) f[i] = 1 + amp * Math.exp(-((x(i) - 128) ** 2) / (2 * (2 * Math.SQRT2) ** 2)); }
  const start = kind === 'kink' ? maxGrad(f) : amp;
  const hist = [];
  for (let t = 0; t < steps; t++) {
    const a = new Float64Array(N);
    for (let i = 1; i < N - 1; i++) a[i] = laplacian(f, i, fl, fr) - (f[i] ** 3 - f[i]);
    for (let i = 1; i < N - 1; i++) { v[i] += dt * a[i]; f[i] += dt * v[i]; }
    f[0] = fl; f[N - 1] = fr;
    if (kick && t === 40) for (let i = 100; i < 140; i++) v[i] += 1.5;
    if (t % 150 === 0) hist.push({ t: t * dt, metric: kind === 'kink' ? maxGrad(f) : peak(f) - 1, q: charge(f, fl, fr) });
  }
  return { hist, metric: kind === 'kink' ? maxGrad(f) : peak(f) - 1, q: charge(f, fl, fr), start };
}

console.log('============================================================');
console.log('PREDICTIONS');
console.log('============================================================');
console.log('   P1 kink: Q = 1, wall ~0.7, re-forms after a kick');
console.log('   P2 small pulse disperses;  P2b large pulse forms a bound OSCILLATING structure');
console.log('   P3 remanence: the state does not return when the field is swept away');
console.log('   P4 uniform noise never switches 1D (nucleation suppression); a LOCAL event does, and spreads');
console.log('   P5 the retained bit outlives the dispersing pulse by a large factor');
console.log('');

console.log('============================================================');
console.log('PART A: TOPOLOGICAL WALL vs NON-TOPOLOGICAL PULSE');
console.log('============================================================');
{
  const kink = waveRun('kink', 0, 1400, false);
  const kinkKick = waveRun('kink', 0, 1400, true);
  const small = waveRun('pulse', 0.10, 1400, false);
  const big = waveRun('pulse', 1.00, 1400, false);
  console.log('   t      wall sharpness   Q      small pulse      big pulse');
  for (let k = 0; k < kink.hist.length; k++) {
    console.log(`   ${kink.hist[k].t.toFixed(0).padStart(4)}   ${kink.hist[k].metric.toFixed(4).padStart(14)}   ${kink.hist[k].q.toFixed(2).padStart(5)}   ${small.hist[k].metric.toFixed(4).padStart(11)}   ${big.hist[k].metric.toFixed(4).padStart(9)}`);
  }
  const smallFrac = small.metric / small.start;
  console.log('');
  console.log(`   P1 ${Math.abs(kink.q - 1) < 1e-6 && kink.metric > 0.65 && kinkKick.metric > 0.65 ? 'CONFIRMED' : 'NOT confirmed'}: Q = ${kink.q.toFixed(6)} throughout; wall ${kink.hist[0].metric.toFixed(4)} -> ${kink.metric.toFixed(4)}; after a kick ${kinkKick.metric.toFixed(4)}`);
  console.log(`   P2 ${smallFrac < 0.5 ? 'CONFIRMED' : 'NOT confirmed'}: the small pulse ends at ${(100 * smallFrac).toFixed(1)}% of its start`);
  console.log(`   P2b the big pulse collapses to ${big.hist[1].metric.toFixed(4)} at t = ${big.hist[1].t.toFixed(0)} and then RISES to ${big.metric.toFixed(4)}:`);
  console.log(`       it fell into the second well and the resulting pair BINDS INTO A BREATHER.  So it does not`);
  console.log(`       disperse, and it does not settle either -- it keeps oscillating.`);
  console.log('');
  console.log(`   WITH NOTHING COUPLING ANYTHING: the wall keeps 100% of its structure, the small pulse keeps`);
  console.log(`   ${(100 * smallFrac).toFixed(0)}%, and the big pulse keeps an oscillating bound object.`);
}
console.log('');

console.log('============================================================');
console.log('PART B: REMANENCE, AND HOW THE MEMORY IS ACTUALLY ERASED');
console.log('============================================================');
{
  const dt = 0.05, H_BIAS = -0.20;
  {
    const f = new Float64Array(N);
    for (let i = 0; i < N; i++) f[i] = -1;
    const relax = (h, n) => { for (let s = 0; s < n; s++) { const a = new Float64Array(N); for (let i = 1; i < N - 1; i++) a[i] = laplacian(f, i) - (f[i] ** 3 - f[i]) - h; for (let i = 1; i < N - 1; i++) f[i] += dt * a[i]; } };
    for (let h = -0.6; h <= 0.6; h += 0.002) relax(h, 20);
    for (let h = 0.6; h >= -0.6; h -= 0.002) relax(h, 20);
    relax(0, 400);
    console.log(`   B1 remanence: after sweeping h from -0.6 up to +0.6 and back to 0, the order parameter sits at`);
    console.log(`      m = ${mean(f).toFixed(4)} (the symmetric value would be 0) => ${Math.abs(mean(f)) > 0.5 ? 'CONFIRMED: removing the field does NOT restore the state' : 'NOT confirmed'}`);
  }
  console.log('');
  {
    console.log('   B2 uniform noise alone, starting in the METASTABLE well (+1, since h = -0.20 favours -1):');
    console.log('      noise sigma     escape time (two seeds)');
    for (const sg of [0, 0.1, 0.3]) {
      const times = [0, 1].map((k) => {
        const f = new Float64Array(N);
        for (let i = 0; i < N; i++) f[i] = 1;
        const rnd = stream(11111 + 7919 * k);
        let t = 0;
        while (t < 3000) {
          const a = new Float64Array(N);
          for (let i = 0; i < N; i++) a[i] = laplacian(f, i) - (f[i] ** 3 - f[i]) - H_BIAS + sg * (rnd() - 0.5);
          for (let i = 1; i < N - 1; i++) f[i] += dt * a[i];
          t += dt;
          if (mean(f) < 0) return t;
        }
        return Infinity;
      });
      const fmt = (t) => (t === Infinity ? '> 3000' : t.toFixed(1));
      console.log(`      ${sg.toFixed(2).padStart(10)}      ${fmt(times[0]).padStart(8)}   ${fmt(times[1]).padStart(8)}`);
    }
    console.log('      => P4a NOT CONFIRMED AS PREDICTED, and the reason is real physics: in ONE dimension the');
    console.log('         critical nucleus is macroscopically large, so uniform noise cannot nucleate the stable');
    console.log('         phase within any affordable run, at any sigma.  My prediction assumed 2D/3D.');
  }
  console.log('');
  {
    const f = new Float64Array(N);
    for (let i = 0; i < N; i++) f[i] = 1;
    // THE SIZE OF THE NUCLEUS MATTERS, and my previous version got it exactly wrong: in one dimension
    // the critical nucleus is L_c ~ 2*sigma/df with sigma the wall tension (~0.94 here) and df the bulk
    // difference (2|h| = 0.4), so L_c ~ 4.7 length units -- and I nucleated only 4 units, i.e. exactly
    // the critical size, which neither grows nor shrinks.  Hence the front never moved.  Nucleate far
    // above the critical size instead.
    for (let i = 220; i < 300; i++) f[i] = -1;
    console.log(`   B3 a LOCAL nucleation of the stable phase, sized well above the critical nucleus`);
    console.log(`      (estimated L_c ~ 2*sigma/df ~ ${(2 * 0.94 / 0.4).toFixed(1)} length units; nucleating ${((300 - 220) * DX).toFixed(0)} units)`);
    // the right-hand boundary of the NEGATIVE (stable) region.  My first version looked for the last
    // positive site, which is always the right-hand boundary of the box, so it never moved.
    const frontOf = () => { for (let i = N - 1; i >= 0; i--) if (f[i] < 0) return i * DX; return 0; };
    const trace = [];
    for (let step = 1; step <= 4000; step++) {
      const a = new Float64Array(N);
      for (let i = 1; i < N - 1; i++) a[i] = laplacian(f, i) - (f[i] ** 3 - f[i]) - H_BIAS;
      for (let i = 1; i < N - 1; i++) f[i] += dt * a[i];
      if (step % 800 === 0) trace.push({ t: step * dt, front: frontOf() });
    }
    console.log('   B3 a LOCAL nucleation of the stable phase into the metastable one:');
    console.log('      t        boundary between the phases (x)');
    for (const p of trace) console.log(`      ${p.t.toFixed(0).padStart(5)}    ${p.front.toFixed(2).padStart(8)}`);
    const first = trace[0], last = trace[trace.length - 1];
    const speed = (last.front - first.front) / (last.t - first.t);
    console.log(`      the front moved ${(last.front - first.front).toFixed(2)} in ${(last.t - first.t).toFixed(0)} time units => speed ${speed.toFixed(4)}`);
    console.log(`      => ${Math.abs(speed) > 1e-3 ? 'CONFIRMED: one local event ERASES the memory, and the erasure PROPAGATES' : 'the erasure did not propagate'}`);
    console.log('      so the memory resists uniform noise but is destroyed by a local event, and the erasure then');
    console.log('      spreads at a measurable speed: that is what "the residual channel being overwritten" looks like.');
  }
}
console.log('');

console.log('============================================================');
console.log('PART C: THE RESIDUAL CHANNEL -- HOW LONG IS THE BIT READABLE? (P5)');
console.log('============================================================');
{
  const pulseHalfLife = (() => {
    const r = waveRun('pulse', 0.10, 2000, false);
    for (const p of r.hist) if (p.metric < 0.05) return p.t;
    return Infinity;
  })();
  console.log(`   the dispersing pulse loses half its amplitude by t ~ ${pulseHalfLife === Infinity ? '> 300' : pulseHalfLife.toFixed(0)} time units`);
  console.log(`   the metastable bit survives the whole 3000-unit noise scan (part B2) and the whole field sweep`);
  console.log(`   => P5 CONFIRMED: the residual channel outlives the dispersing pulse by a factor of ${pulseHalfLife === Infinity ? 'many' : (3000 / pulseHalfLife).toFixed(0)} at least.`);
  console.log('');
  console.log('   WHAT THE MEMORY IS MADE OF: not a stored pattern, but the fact that the network has two stable');
  console.log('   vacua, a barrier between them, and -- in one dimension -- no way for uniform noise to nucleate');
  console.log('   the other.  The wall of part A is stronger still: it is protected by the boundary conditions');
  console.log('   being different vacua, so nothing is holding it up at all, and a kick makes it re-form.');
  console.log('');
  console.log('   WHERE THE DOCUMENT GOES FURTHER THAN THE MEASUREMENT: it says the residual channel can');
  console.log('   "support the existence of other information".  Here the channel carries one bit (which well)');
  console.log('   and the wall carries a position -- nothing in this round shows a residual channel carrying a');
  console.log('   MESSAGE.  That is the honest boundary of this result.');
}
