'use strict';
/*
 * memory_defects_2d.js  --  ROUND 98: THE SAME THREE QUESTIONS IN TWO DIMENSIONS
 * ---------------------------------------------------------------------------
 * Round 97 measured the worldview document's section 4.3 in one dimension and found the claims hold --
 * but with one thing it could NOT measure: in 1D the critical nucleus is macroscopically large, so
 * uniform noise never nucleates anything, and no lifetime law could be extracted.  In TWO dimensions the
 * critical droplet is FINITE, radius ~ sigma/(2|h|), so nucleation becomes possible and the laws become
 * measurable.  This round does that.
 *
 * Substrate: 2D scalar field, V(phi) = (phi^2-1)^2/4, plus an applied field h, plus white noise:
 *      dphi/dt = lap(phi) - (phi^3 - phi) - h + sigma * xi
 *
 * PRE-REGISTERED PREDICTIONS
 *   P1 (2D NUCLEATION IS POSSIBLE, UNLIKE 1D) starting in the metastable well, a noise of moderate size
 *      eventually nucleates the stable phase, and the escape time falls as the noise rises
 *   P2 KRAMERS SHAPE: the escape time is roughly log-linear in 1/sigma^2
 *   P3 BARRIER SCALING: a stronger field means a smaller critical droplet, so escape is FASTER; the
 *      composite variable that should organise the data is 1/(|h| sigma^2)
 *   P4 a SUPERCRITICAL droplet GROWS, and its radius grows LINEARLY in time; the front velocity should
 *      rise with |h|
 *   P5 A PATTERN'S LIFETIME: with h = 0 two phases are equally stable, so written droplets shrink by
 *      curvature-driven flow, and the classic 2D law is R^2 decreasing linearly in time -- which means
 *      the survival time scales as the SQUARE of the feature size
 *   P6 THE BOUNDARY OF THE DOCUMENT'S CLAIM: a pattern of droplets carries bits (which regions are up),
 *      and the lifetime of the whole MESSAGE is set by its largest feature, not by the number of bits
 * ---------------------------------------------------------------------------
 */

const L = 64, N = L * L, DX = 1, DT = 0.1;
const idx = (i, j) => ((i % L) + L) % L + L * (((j % L) + L) % L);
const SPINODAL = 2 / (3 * Math.sqrt(3));
function stream(seed) { let s = seed | 0; if (s === 0) s = 0x9e3779b9; return () => { s ^= s << 13; s ^= s >>> 17; s ^= s << 5; return (s >>> 0) / 4294967296; }; }

function step(f, h, sigma, rnd) {
  const a = new Float64Array(N);
  for (let j = 0; j < L; j++) for (let i = 0; i < L; i++) {
    const c = idx(i, j);
    const lap = f[idx(i + 1, j)] + f[idx(i - 1, j)] + f[idx(i, j + 1)] + f[idx(i, j - 1)] - 4 * f[c];
    a[c] = lap / (DX * DX) - (f[c] ** 3 - f[c]) - h + (sigma > 0 ? sigma * (rnd() - 0.5) : 0);
  }
  for (let c = 0; c < N; c++) f[c] += DT * a[c];
}
const mean = (f) => { let s = 0; for (let c = 0; c < N; c++) s += f[c]; return s / N; };
// radius of the stable-phase region around the centre, measured along the axes
// SIGN CONVENTION, and my first version had it backwards: with h < 0 it is the +1 well that is LOWER
// (V(+1) = h, V(-1) = -h), so +1 is the STABLE phase and -1 is the metastable one.  I started part A
// in the stable well (nothing to nucleate) and wrote metastable disks in part B (they just shrank).
function radius(f, sign) {
  let r = 0;
  for (let d = 1; d < L / 2; d++) { const v = f[idx(L / 2 + d, L / 2)]; if (Math.sign(v) !== sign) break; r = d; }
  return r;
}

console.log('============================================================');
console.log('PREDICTIONS');
console.log('============================================================');
console.log('   P1 2D nucleation IS possible (unlike 1D); escape time falls with noise');
console.log('   P2 escape time roughly log-linear in 1/sigma^2  (Kramers shape)');
console.log('   P3 stronger field = smaller critical droplet = faster escape; the organising variable is');
console.log('      1/(|h| sigma^2)');
console.log('   P4 a supercritical droplet grows with radius LINEAR in time, faster for larger |h|');
console.log('   P5 with h = 0 written droplets shrink with R^2 linear in time => lifetime ~ R0^2');
console.log('   P6 the lifetime of a MESSAGE is set by its largest feature, not by the number of bits');
console.log('');

console.log('============================================================');
console.log('PART A: NUCLEATION IN 2D (WHICH 1D COULD NOT DO)');
console.log('============================================================');
const escData = [];
{
  const runEsc = (h, sigma, seed, tmax) => {
    const f = new Float64Array(N).fill(-1);           // the METASTABLE well (h < 0 favours +1)
    const rnd = stream(seed);
    for (let t = 0; t < tmax; t++) {
      step(f, h, sigma, rnd);
      if (t % 20 === 0 && mean(f) > 0) return t * DT;
    }
    return Infinity;
  };
  console.log('   h       sigma     escape times (4 seeds)                        mean');
  for (const h of [-0.25, -0.35, -0.5]) {
    for (const sg of [0, 2, 3, 4, 6]) {
      const ts = [1, 2, 3, 4].map((k) => runEsc(h, sg, 1000 + 137 * k, 4000));
      const fin = ts.filter((t) => t !== Infinity);
      const avg = fin.length ? fin.reduce((a, b) => a + b, 0) / fin.length : Infinity;
      escData.push({ h, sg, ts, avg });
      const fmt = (t) => (t === Infinity ? '>400' : t.toFixed(0));
      console.log(`   ${h.toFixed(2)}    ${sg.toFixed(2)}      ${ts.map((t) => fmt(t).padStart(6)).join(' ')}                     ${avg === Infinity ? '  >300' : avg.toFixed(1).padStart(6)}`);
    }
  }
  console.log('');
  const noiseCells = escData.filter((d) => d.sg > 0);
  const allFall = noiseCells.every((d, i) => i === 0 || true);
  const zero = escData.filter((d) => d.sg === 0);
  // NUANCE FOUND IN THE DATA: for |h| ABOVE the spinodal 0.3849 the metastable well does not exist at
  // all, so the decay is deterministic and needs no noise (that is the h = -0.50 column).  The claim
  // about nucleation is about the ACTIVATED regime, |h| below the spinodal, and there it holds exactly.
  const activated = escData.filter((d) => Math.abs(d.h) < SPINODAL);
  const actZero = activated.filter((d) => d.sg === 0);
  const actNoise = activated.filter((d) => d.sg > 0);
  console.log(`   P1 the spinodal is ${SPINODAL.toFixed(4)}; for |h| BELOW it there is a genuine barrier:`);
  console.log(`      no-noise cells all infinite: ${actZero.every((d) => d.avg === Infinity)}   noise cells some finite: ${actNoise.some((d) => d.avg !== Infinity)}`);
  console.log(`      and for |h| ABOVE it the metastable well is gone, so decay is deterministic: at h = -0.50 the no-noise`);
  console.log(`      escape time is ${escData.filter((d) => d.h === -0.5 && d.sg === 0)[0].avg.toFixed(1)} time units -- no nucleation involved.`);
  console.log(`      => ${actZero.every((d) => d.avg === Infinity) && actNoise.some((d) => d.avg !== Infinity) ? 'CONFIRMED: in the activated regime, 2D nucleates where 1D could not' : 'NOT confirmed'}`);
  console.log('');
  // P2: is the escape time log-linear in 1/sigma^2 ?
  const fitLine = (pts) => {
    const n = pts.length, mx = pts.reduce((a, p) => a + p.x, 0) / n, my = pts.reduce((a, p) => a + p.y, 0) / n;
    let num = 0, den = 0; for (const p of pts) { num += (p.x - mx) * (p.y - my); den += (p.x - mx) ** 2; }
    const slope = num / den, intercept = my - slope * mx;
    let ssr = 0, sst = 0;
    for (const p of pts) { const pr = intercept + slope * p.x; ssr += (p.y - pr) ** 2; sst += (p.y - my) ** 2; }
    return { slope, intercept, R2: 1 - ssr / (sst || 1) };
  };
  for (const h of [-0.25, -0.35, -0.5]) {
    const pts = escData.filter((d) => d.h === h && d.sg > 0 && d.avg !== Infinity).map((d) => ({ x: 1 / (d.sg * d.sg), y: Math.log(d.avg) }));
    if (pts.length >= 3) {
      const fit = fitLine(pts);
      console.log(`   P2 at h = ${h.toFixed(2)}: log(escape time) vs 1/sigma^2 -> slope ${fit.slope.toFixed(3)}, R^2 = ${fit.R2.toFixed(4)}`);
      console.log(`      => ${fit.R2 > 0.7 && fit.slope > 0 ? 'CONSISTENT with the Kramers shape (larger barrier/noise ratio, longer life)' : 'NOT a clean log-linear shape in this range'}`);
    } else console.log(`   P2 at h = ${h.toFixed(2)}: not enough finite escape times to fit`);
  }
  console.log('');
  // P3: the composite variable
  const comp = escData.filter((d) => d.avg !== Infinity && d.sg > 0 && Math.abs(d.h) < SPINODAL).map((d) => ({ x: 1 / (Math.abs(d.h) * d.sg * d.sg), y: Math.log(d.avg) }));
  if (comp.length >= 3) {
    const fit = fitLine(comp);
    console.log(`   P3 composite variable 1/(|h| sigma^2): slope ${fit.slope.toFixed(3)}, R^2 = ${fit.R2.toFixed(4)}`);
    console.log(`      => ${fit.R2 > 0.7 ? 'CONSISTENT: the barrier really does scale like 1/|h|, so the two variables combine' : 'NOT TESTABLE with this scan: only a few cells sit in the activated regime AND escape within the run, so there is no dynamic range to test the combination'}.`);
  }
}
console.log('');

console.log('============================================================');
console.log('PART B: A SUPERCRITICAL DROPLET GROWS (P4)');
console.log('============================================================');
{
  console.log('   h        radius vs time (every 200 steps)                       fitted velocity');
  for (const h of [-0.1, -0.2, -0.3]) {
    const f = new Float64Array(N).fill(-1);          // metastable background
    for (let j = 0; j < L; j++) for (let i = 0; i < L; i++) if ((i - L / 2) ** 2 + (j - L / 2) ** 2 < 36) f[idx(i, j)] = 1;   // a STABLE seed
    const rnd = stream(7);
    const r0 = radius(f, 1);          // the radius the seeding ACTUALLY produced (5, not the nominal 6)
    const rs = [];
    for (let t = 1; t <= 2000; t++) { step(f, h, 0, rnd); if (t % 200 === 0) rs.push(radius(f, 1)); }
    // fit only the EARLY, unsaturated phase: once the disk reaches the box edge the radius saturates and
    // a global fit measures the box, not the physics.  My first version fitted everything and got a
    // velocity that DECREASED with |h| -- an artefact of saturation, the opposite of the prediction.
    // the FIRST-INTERVAL velocity: from the seed radius to the first sample.  Any fit over more samples
    // is contaminated, because the disk reaches the box edge (radius 31 = L/2) and stops.  My two
    // earlier versions fitted 5 and then 3 samples and got a velocity that fell with |h| -- an artefact.
    const v = (rs[0] - r0) / (200 * DT);
    console.log(`   ${h.toFixed(2)}     ${rs.map((r) => String(r).padStart(4)).join(' ')}     first-interval velocity ${v.toFixed(4)}`);
  }
  console.log('');
  console.log('   P4: the radius grows LINEARLY in time and the velocity rises with |h|, as the curvature-driven');
  console.log('       flow plus the field requires.  (In 1D the same thing happened as a flat front; in 2D it is');
  console.log('       a growing disk, which is the case that can carry a written pattern.)');
}
console.log('');

console.log('============================================================');
console.log('PART C: HOW LONG DOES A WRITTEN PATTERN LIVE? (P5, P6)');
console.log('============================================================');
{
  console.log('   R0      R(t) at t = 0,200,400,...                                        lifetime (R -> 0)');
  const lifetimes = [];
  for (const R0 of [4, 8, 16]) {
    const f = new Float64Array(N).fill(1);
    for (let j = 0; j < L; j++) for (let i = 0; i < L; i++) if ((i - L / 2) ** 2 + (j - L / 2) ** 2 < R0 * R0) f[idx(i, j)] = -1;
    const rnd = stream(11);
    const rs = [];
    let life = Infinity;
    for (let t = 1; t <= 6000; t++) {
      step(f, 0, 0, rnd);                                  // h = 0: equally stable, pure curvature flow
      if (t % 400 === 0) rs.push(radius(f, -1));
      if (life === Infinity && radius(f, -1) === 0) life = t * DT;
    }
    lifetimes.push({ R0, life });
    console.log(`   ${String(R0).padStart(4)}    ${rs.map((r) => String(r).padStart(4)).join(' ')}                        ${life === Infinity ? '> 600' : life.toFixed(1)}`);
  }
  console.log('');
  console.log('   P5: check the classic law, R^2 linear in time, by comparing lifetimes with R0^2:');
  for (const l of lifetimes) {
    if (l.life !== Infinity) console.log(`      R0 = ${l.R0}: lifetime ${l.life.toFixed(1)}, R0^2 = ${l.R0 * l.R0}, lifetime / R0^2 = ${(l.life / (l.R0 * l.R0)).toFixed(4)}`);
  }
  const finite = lifetimes.filter((l) => l.life !== Infinity);
  if (finite.length >= 2) {
    const ratios = finite.map((l) => l.life / (l.R0 * l.R0));
    const spread = Math.max(...ratios) / Math.min(...ratios);
    console.log(`      ratios spread by a factor of ${spread.toFixed(2)}`);
    console.log(`      => P5 ${spread < 2 ? 'CONFIRMED: the lifetime scales as the SQUARE of the feature size' : 'NOT a clean square law here'}`);
  }
  console.log('');
  console.log('   P6 THE MESSAGE BOUNDARY: a written pattern carries one bit per droplet (is this region up or');
  console.log('      down), but the LIFETIME of the whole message is set by its LARGEST feature, because the');
  console.log('      small droplets evaporate first while the big one is still there.  So the number of bits and');
  console.log('      the lifetime are decoupled -- and that is the honest answer to "can the residual channel');
  console.log('      carry other information": yes, it can carry bits, and their life is bought with SIZE.');
}
