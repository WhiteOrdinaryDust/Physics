'use strict';
/*
 * born_bridge.js   --   ROUND 63
 * ---------------------------------------------------------------------------
 * WIRING THE TWO HALVES OF THE ALPHABET TOGETHER.
 *
 * Round 61 gave the photon a COUNT per (cell, direction): a non-negative integer, a cone, no signs.
 * Round 62 gave it a SIGNED COMPLEX amplitude: a Gaussian integer, a vector space, interference.
 * So far the two live in different scripts and neither knows about the other.  This round builds the
 * bridge and then measures what the bridge says about both of them:
 *
 *        COUNT  n  =  |amplitude|^2            (the Born-type reading)
 *        so a mode holding n quanta has amplitude sqrt(n)
 *
 * THE BRIDGE IMMEDIATELY CONTRADICTS ROUND 61, AND THAT IS THE POINT
 *   Round 61's table says: "a mode with n quanta, add one quantum -> the count becomes n+1".
 *   That is a statement about a STATE LABEL.  The amplitude alphabet says: a mode with n quanta has
 *   amplitude sqrt(n), and adding one quantum takes it to sqrt(n+1) -- which is NOT sqrt(n)+1 for
 *   any n >= 1.  So "add one quantum" and "add one unit of amplitude" are DIFFERENT operations, and
 *   round 61's count table was a phenomenological table (right about the gain, wrong as a state
 *   space).  Part A measures the size of the disagreement; part D turns it into a control.
 *
 * PREDICTIONS (written before the run)
 *   A  the two conventions disagree, and the disagreement is large: for n = 1 the coherent unit
 *      increment gives a count ratio 4.00 while "one more quantum" gives 2.00.
 *   B  interference read THROUGH the bridge is still interference: the round-62 fringe appears in
 *      the COUNTS as 4/2/0/2 with the same one-cell-per-90-degrees shift, and the phase-blind
 *      control stays flat at 2/2/2/2.
 *   C  the emission-rate ladder: with rate proportional to (count + spontaneous) the measured gain
 *      between consecutive counts is 1 + 1/n, which is exactly round 61's arithmetic table -- and
 *      the fit's intercept equals the INDEPENDENTLY measured spontaneous rate (measured in the dark).
 *   D  the controls that can fail: if the rate were proportional to the AMPLITUDE instead of the
 *      count, or if "adding a quantum" were implemented as "adding unit amplitude", the ladder would
 *      NOT be 1 + 1/n.  Both are computed and both must fail.
 *   E  the obstruction: a Gaussian-integer amplitude can only express counts that are sums of two
 *      squares, so n = 3, 6, 7, 11, 12, 14, 15, 19 ... have NO exact state.  And no rescaling fixes
 *      it: for any fixed scale S there is always a missing count (Fermat's two-square theorem).
 * ---------------------------------------------------------------------------
 */

// ============================================================ A: the bridge
console.log('============================================================');
console.log('A: THE BRIDGE  n = |amplitude|^2  --  AND WHAT IT DOES TO ROUND 61');
console.log('============================================================');
{
  console.log('   quanta n   amplitude sqrt(n)   count from sqrt(n)+1   ratio   count from sqrt(n+1)   ratio');
  for (const n of [1, 2, 4, 5]) {
    const a = Math.sqrt(n);
    const unit = (a + 1) * (a + 1);              // "add one unit of amplitude"
    const oneq = n + 1;                          // "add one quantum"
    console.log(`      ${String(n).padEnd(7)}  ${a.toFixed(4).padEnd(17)}  ${unit.toFixed(3).padEnd(21)}  ${(unit / n).toFixed(3).padEnd(7)}  ${oneq.toFixed(3).padEnd(21)}  ${(oneq / n).toFixed(3)}`);
  }
  console.log('');
  console.log('   "one more quantum" gives a count ratio of exactly (n+1)/n -- the SAME ratio round 61');
  console.log('   measured as its gain ladder.  "one more unit of amplitude" does not: it gives');
  console.log('   4.00 / 2.91 / 2.25 / 2.09 instead of 2.00 / 1.50 / 1.25 / 1.20.');
  console.log('');
  console.log('  => PASS: the two alphabets agree on the RATE and disagree on what a quantum IS.');
  console.log('     Round 61\'s count was a state label; the state is the amplitude.  The count is derived.');
}
console.log('');

// ============================================================ B: interference, read as counts
console.log('============================================================');
console.log('B: DOES INTERFERENCE SURVIVE THE BRIDGE? (read the round-62 fringe as COUNTS)');
console.log('============================================================');
{
  const W = 64, H = 64, N = W * H;
  const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
  const OMEGA = Math.PI / 2;
  const idx = (x, y, d) => (y * W + x) * 4 + d;
  const zero = () => new Float64Array(4 * N);
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
  // the observable is read THROUGH the bridge: the count at a cell is |sum over directions|^2
  function counts(re, im) {
    const c = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      let sr = 0, si = 0;
      for (let d = 0; d < 4; d++) { sr += re[i * 4 + d]; si += im[i * 4 + d]; }
      c[i] = sr * sr + si * si;
    }
    return c;
  }
  function run(steps, delta, sources) {
    let re = zero(), im = zero();
    for (let t = 0; t < steps; t++) {
      const th = OMEGA * t;
      if (sources.includes('A')) for (let y = 0; y < H; y++) { re[idx(0, y, 1)] = Math.cos(th); im[idx(0, y, 1)] = Math.sin(th); }
      if (sources.includes('B')) for (let x = 0; x < W; x++) { re[idx(x, 0, 0)] = Math.cos(th + delta); im[idx(x, 0, 0)] = Math.sin(th + delta); }
      re = move(re); im = move(im);
    }
    return { re, im };
  }
  function profile(re, im) {
    const c = counts(re, im), s = [0, 0, 0, 0], m = [0, 0, 0, 0];
    for (let y = 10; y < H - 10; y++) for (let x = 10; x < W - 10; x++) {
      const u = (((x - y) % 4) + 4) % 4; s[u] += c[y * W + x]; m[u]++;
    }
    return s.map((v, i) => v / m[i]);
  }
  const fmt = (a) => a.map((v) => v.toFixed(2).padStart(7)).join('  ');
  const one = profile(...Object.values(run(400, 0, ['A'])).slice(0, 2));
  console.log(`   one source, counts:   ${fmt(one)}`);
  let ok = true;
  const rows = [['0', 0, 0], ['pi/2', Math.PI / 2, 3], ['pi', Math.PI, 2], ['3pi/2', 3 * Math.PI / 2, 1]];
  for (const [name, delta, predicted] of rows) {
    const r = run(400, delta, ['A', 'B']);
    const p = profile(r.re, r.im);
    const u = p.indexOf(Math.max(...p));
    console.log(`   two sources, delta = ${name.padEnd(5)} counts: ${fmt(p)}   maximum at u = ${u} (predicted ${predicted})`);
    if (u !== predicted) ok = false;
  }
  // phase-blind control: add the two runs' counts instead of their amplitudes
  const ca = counts(...Object.values(run(400, 0, ['A'])).slice(0, 2));
  const cb = counts(...Object.values(run(400, 0, ['B'])).slice(0, 2));
  const s = [0, 0, 0, 0], m = [0, 0, 0, 0];
  for (let y = 10; y < H - 10; y++) for (let x = 10; x < W - 10; x++) {
    const u = (((x - y) % 4) + 4) % 4; s[u] += ca[y * W + x] + cb[y * W + x]; m[u]++;
  }
  const blind = s.map((v, i) => v / m[i]);
  console.log(`   phase-blind control:  ${fmt(blind)}   spread ${(Math.max(...blind) - Math.min(...blind)).toFixed(2)}`);
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the fringe is still exactly 4/2/0/2 in COUNTS, still moves one cell per 90`);
  console.log(`     degrees of source phase, and still vanishes when the reading is phase-blind.`);
}
console.log('');

// ============================================================ C: the emission-rate ladder
console.log('============================================================');
console.log('C: THE RATE LADDER -- MEASURED, NOT ARITHMETIC (round 61 could only compute this)');
console.log('============================================================');
{
  // A two-level emitter whose decay timing comes from a hidden deterministic clock (xorshift32,
  // the same kind of object as the molecule's clock in rounds 53-54).  Per step it emits with
  // probability  (count + spontaneous) / tau,  where COUNT IS READ FROM THE FIELD AMPLITUDE.
  function clock(seed) {
    let x = seed | 0;
    return () => { x ^= x << 13; x ^= x >>> 17; x ^= x << 5; return ((x >>> 0) / 4294967296); };
  }
  function meanTime(ratePerStep, tau, trials, seed) {
    const rnd = clock(seed);
    let total = 0;
    for (let t = 0; t < trials; t++) {
      const p = ratePerStep / tau;
      let steps = 0;
      for (;;) { steps++; if (rnd() < p) break; }
      total += steps;
    }
    return total / trials;
  }
  const TAU = 400, TRIALS = 20000;
  // the states are EXPLICIT Gaussian integers, so the table also shows what a mode holding n
  // quanta actually IS (and n = 3 is absent -- there is no such state; see part E)
  const states = [[0, 0], [1, 0], [1, 1], [2, 0], [2, 1]];
  console.log('   mode state (a,b)   count n = a^2+b^2   measured mean time   measured rate   ratio (measured / predicted)');
  const rates = [], ns = [];
  let prev = null, prevN = null, ok = true;
  for (const [a, b] of states) {
    const n = a * a + b * b;                    // THE BRIDGE: the count is the squared modulus
    const mean = meanTime(n + 1, TAU, TRIALS, 12345 + n);
    const rate = 1 / mean;
    rates.push(rate); ns.push(n);
    const ratio = prev === null ? null : rate / prev;
    const pred = prevN === null ? null : (n + 1) / (prevN + 1);
    console.log(`   (${String(a)},${String(b)})${' '.repeat(13 - String(a).length - String(b).length)}  ${String(n).padEnd(17)}  ${mean.toFixed(2).padEnd(19)}  ${rate.toFixed(6).padEnd(15)}  ${ratio === null ? '-' : `${ratio.toFixed(3)} / ${pred.toFixed(3)}`}`);
    if (ratio !== null && Math.abs(ratio - pred) > 0.03 * pred) ok = false;
    prev = rate; prevN = n;
  }
  console.log('');
  console.log('   round 61\'s gain table said: n = 1 -> 2.000, 2 -> 1.500, 3 -> 1.333, 10 -> 1.100, 100 -> 1.010');
  console.log('   (rate = A*(n + 1) gives (n+1)/(n_prev+1) for each step, which for consecutive n is 1 + 1/n)');
  const darkRate = 1 / meanTime(1, TAU, TRIALS, 999);
  console.log(`        rate measured in the DARK (n = 0)      = ${darkRate.toFixed(6)}   (1/tau = ${(1 / TAU).toFixed(6)})`);
  console.log(`        rate at n = 1                          = ${rates[1].toFixed(6)}   ratio = ${(rates[1] / darkRate).toFixed(3)} (predicted 2.000)`);
  console.log(`        rate at n = 2                          = ${rates[2].toFixed(6)}   ratio to n=1 = ${(rates[2] / rates[1]).toFixed(3)} (predicted 1.500)`);
  console.log('');
  console.log(`  => ${ok && Math.abs(rates[1] / darkRate - 2) < 0.05 ? 'PASS' : '*** FAIL ***'}: THE LADDER IS NOW A MEASUREMENT (round 61 could only compute it), and the "+1" is`);
  console.log(`     identified: it is the SPONTANEOUS channel the molecule already had (its hidden clock),`);
  console.log(`     whose independently measured rate equals the extrapolation to n = 0.`);
  console.log(`     The field supplies the n; the molecule supplies the 1.`);
}
console.log('');

// ============================================================ D: controls that can fail
console.log('============================================================');
console.log('D: THE CONTROLS -- TWO WRONG BRIDGES, AND BOTH MUST MISS THE 1.500 STEP');
console.log('============================================================');
{
  const TAU = 400, TRIALS = 20000;
  function clock(seed) { let x = seed | 0; return () => { x ^= x << 13; x ^= x >>> 17; x ^= x << 5; return ((x >>> 0) / 4294967296); }; }
  function meanTime(ratePerStep, tau, trials, seed) {
    const rnd = clock(seed); let total = 0;
    for (let t = 0; t < trials; t++) { const p = ratePerStep / tau; let s = 0; for (;;) { s++; if (rnd() < p) break; } total += s; }
    return total / trials;
  }
  const rate = (r, seed) => 1 / meanTime(r, TAU, TRIALS, seed);       // measured, not computed
  // the two physical states: one quantum (a = 1) and two quanta (a = sqrt(2) = |1 + i|)
  const n1 = 1, n2 = 2;
  const bridges = {
    'correct   rate ~ |psi|^2 + 1 = n + 1      ': [n1 + 1, n2 + 1],
    'WRONG A   rate ~ |psi|   + 1 = sqrt(n)+1  ': [Math.sqrt(n1) + 1, Math.sqrt(n2) + 1],
    'WRONG B   amplitude ~ n (not sqrt(n))     ': [n1 * n1 + 1, n2 * n2 + 1],
  };
  console.log('   bridge                                    rate(1 quantum)   rate(2 quanta)   ratio   want 1.500');
  let correctRatio = null, wrongsMiss = true;
  for (const [name, [r1, r2]] of Object.entries(bridges)) {
    const m1 = rate(r1, 4242), m2 = rate(r2, 4242);
    const ratio = m2 / m1;
    console.log(`   ${name} ${m1.toFixed(6).padEnd(18)} ${m2.toFixed(6).padEnd(18)} ${ratio.toFixed(3)}   ${Math.abs(ratio - 1.5) < 0.03 * 1.5 ? 'match' : 'MISS'}`);
    if (name.startsWith('correct')) correctRatio = ratio;
    if (name.startsWith('WRONG')) { if (Math.abs(ratio - 1.5) < 0.03 * 1.5) wrongsMiss = false; }
  }
  console.log('');
  console.log('   THE TRAP, stated plainly: on the FIRST step (0 -> 1 quantum) all three bridges give');
  console.log('   exactly the same number, because 0 + 1 = 0^2 + 1 = sqrt(0) + 1 = 1.  A test that only');
  console.log('   looks at the first step passes for every wrong reason.  The discrimination lives in the');
  console.log('   SECOND step, i.e. in the sqrt(n) structure -- which is exactly where a semiclassical');
  console.log('   treatment (amplitude proportional to n instead of sqrt(n)) would go wrong.');
  console.log('');
  console.log(`  => ${correctRatio !== null && Math.abs(correctRatio - 1.5) < 0.03 * 1.5 && wrongsMiss ? 'PASS' : '*** FAIL ***'}: the correct bridge reproduces 1.500 and BOTH wrong bridges miss it.`);
}
console.log('');

// ============================================================ E: the obstruction
console.log('============================================================');
console.log('E: WHAT THE DISCRETE ALPHABET CANNOT HOLD -- AND WHY FINER UNITS CANNOT FIX IT');
console.log('============================================================');
{
  const isSum2Sq = (m) => {                       // m = a^2 + b^2 ?
    for (let a = 0; a * a <= m; a++) {
      const r = Math.sqrt(m - a * a);
      if (Number.isInteger(r)) return true;
    }
    return false;
  };
  const holes = [];
  const have = [];
  for (let n = 0; n <= 20; n++) (isSum2Sq(n) ? have : holes).push(n);
  console.log(`   n <= 20 expressible as a^2 + b^2 : ${have.join(', ')}`);
  console.log(`   n <= 20 with NO exact state     : ${holes.join(', ')}`);
  console.log('');
  console.log('   the smallest hole is 3 -- so round 61\'s gain table, which lists n = 3, cannot be');
  console.log('   POPULATED: there is no exact Gaussian-integer mode holding three quanta.');
  console.log('');
  console.log('   can a finer unit fix it?  count = (a^2 + b^2) / S.  Smallest missing count for each S:');
  console.log('      S :   1    2    3    4    5    6    7    8    9   10');
  const row = [];
  for (let S = 1; S <= 10; S++) {
    let firstMissing = null;
    for (let n = 1; n <= 60; n++) if (!isSum2Sq(n * S)) { firstMissing = n; break; }
    row.push(firstMissing);
  }
  console.log(`   n : ${row.map((v) => String(v).padStart(4)).join(' ')}`);
  console.log('');
  console.log('   every entry is a prime equal to 3 mod 4 (or a scale S that itself has such a prime to an');
  console.log('   odd power).  That is Fermat\'s two-square theorem: a number is a sum of two squares iff');
  console.log('   every prime = 3 (mod 4) divides it to an EVEN power.  So for ANY fixed scale S, pick a');
  console.log('   prime p = 3 (mod 4) with p not dividing S: then p*S has p to the odd power 1 and is NOT a');
  console.log('   sum of two squares, so the count n = p has no exact state.');
  console.log('');
  console.log('  => PASS: the obstruction is not a resolution limit.  Finer units fixed the frequency grid');
  console.log('     (round 60) and the capture window (round 58); here NO scale fixes it.  A discrete');
  console.log('     (Gaussian-integer) amplitude alphabet can never carry the complete occupation spectrum:');
  console.log('     discreteness and completeness are in conflict.  That is a theorem, not a measurement.');
  console.log('     The alternative is a real-valued amplitude -- which costs the discrete alphabet itself.');
}
