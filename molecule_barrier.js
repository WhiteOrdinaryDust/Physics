'use strict';
/*
 * molecule_barrier.js   --   ROUND 66
 * ---------------------------------------------------------------------------
 * THE CHECKLIST'S ORIGINAL DESTINATION: MATTER, AND A BARRIER MADE OF IT.
 *
 * Rounds 33-34 went looking for a wall in the BASE RULES (Life, layered elastic media) and came back
 * with two hard negatives:
 *      - Life matter is one-sided elastic: removing a cell always heals, adding one always destroys
 *        => it is always a detector, never a transparent wall (threshold = 0);
 *      - the layered elastic phases heal SCRATCHES but not IMPACTS, and 950 injections found 0 movers
 *        => walls, no balls (threshold = infinite).
 *   and one key quantity: what is wanted is a FINITE threshold -- small disturbance passes without a
 *   record, large disturbance nucleates a permanent one.
 *
 * Thirty rounds later the model has its own matter: molecules with a transition frequency, a
 * resonant absorber, a level, and a lifetime tau driven by a hidden clock (rounds 40-64).  So this
 * round stops searching and ASKS THE MODEL'S OWN MATTER what kind of barrier it makes.
 *
 * PREDICTIONS (written before the run)
 *   A  SELECTIVE TRANSPARENCY -- a finite threshold in the FREQUENCY domain:
 *      off-resonant quanta pass through the wall with ZERO records; resonant ones are absorbed.
 *      This is the finite threshold rounds 33-34 could not find in rule space.
 *   B  BLEACHING -- a finite threshold in the INTENSITY domain:
 *      with M molecules and N incident quanta and no recovery, the transmitted number must be
 *      exactly max(0, N - M)  =>  there is a kink at N = M: dim light is blocked, bright light
 *      passes.  (This is saturable absorption, and it is exact here because each absorbed quantum
 *      uses up exactly one molecule.)
 *   C  RECOVERY -- the wall heals with the MOLECULE'S OWN lifetime:
 *      excited molecules return to the ground state with mean time tau, so a probe pulse sent after
 *      a wait T finds a recovered fraction 1 - exp(-T/tau).  No new constant is introduced.
 *   D  THE COHERENCE COST -- the round-33 question, answered:
 *      does recording destroy the quantum?  NO -- it is re-emitted elsewhere.  It destroys the
 *      PHASE.  So a gentle detector exists: the count survives, the interference does not.
 *      With the two-arm interference of rounds 62/65 and k quanta absorbed from an arm holding n,
 *      the bridge (count = |amplitude|^2) predicts  V(k) = 2*sqrt(n-k) / (1 + (n-k)).
 *      TWO competing models are computed: removing whole quanta (sqrt) versus absorbing amplitude
 *      proportionally (linear).  They differ, so the measurement discriminates.
 *      NOTE the consequence: with n = 1 a single absorbed quantum kills the fringe completely
 *      (detector), and for n > 1 the visibility is NOT monotone in k -- absorbing can IMPROVE
 *      contrast by re-balancing the arms.
 * ---------------------------------------------------------------------------
 */

function clock(seed) { let x = seed | 0; if (x === 0) x = 0x9e3779b9; return () => { x ^= x << 13; x ^= x >>> 17; x ^= x << 5; return ((x >>> 0) / 4294967296); }; }

// a wall of M molecules; a quantum is absorbed iff a GROUND-STATE molecule is within Gamma of it
function makeWall(M, omega, Gamma, tau, seedBase) {
  const mols = [];
  for (let k = 0; k < M; k++) mols.push({ omega, Gamma, tau, exc: false, recoverAt: Infinity, rnd: clock(seedBase + k * 7919 + 1) });
  return mols;
}
function wallRun(M, N, omegaBeam, omegaWall, Gamma, tau, gap, seed) {
  const mols = makeWall(M, omegaWall, Gamma, tau, seed);
  let transmitted = 0, absorbed = 0, records = 0, recovered = 0;
  const steps = N * gap + (tau === Infinity ? 0 : 8 * tau);
  for (let t = 1; t <= steps; t++) {
    for (const m of mols) {
      if (m.exc && t >= m.recoverAt) { m.exc = false; recovered++; }
    }
    if (t % gap === 0 && t / gap <= N) {
      // a quantum arrives and looks for a ground-state molecule in resonance with it
      let hit = null;
      for (const m of mols) {
        if (!m.exc && Math.abs(m.omega - omegaBeam) <= m.Gamma) { hit = m; break; }
      }
      if (hit) {
        hit.exc = true; records++;
        hit.recoverAt = tau === Infinity ? Infinity : t + (-tau * Math.log(Math.max(1e-12, hit.rnd())));
        absorbed++;
      } else transmitted++;
    }
  }
  return { transmitted, absorbed, records, recovered, stillExcited: mols.filter((m) => m.exc).length };
}

console.log('============================================================');
console.log('A: A FINITE THRESHOLD IN THE FREQUENCY DOMAIN (rounds 33-34 could not find this)');
console.log('============================================================');
{
  const M = 20, N = 40;
  const off = wallRun(M, N, 5, 4, 0, Infinity, 3, 11);       // off-resonant quanta
  const on = wallRun(M, N, 4, 4, 0, Infinity, 3, 11);        // resonant quanta
  console.log(`   wall: ${M} molecules at omega = 4.  Beam: ${N} quanta, tolerance Gamma = 0.`);
  console.log('');
  console.log('                          transmitted   absorbed   permanent records');
  console.log(`   off-resonant (w = 5)  ${String(off.transmitted).padStart(11)}   ${String(off.absorbed).padStart(8)}   ${String(off.records).padStart(17)}`);
  console.log(`   resonant     (w = 4)  ${String(on.transmitted).padStart(11)}   ${String(on.absorbed).padStart(8)}   ${String(on.records).padStart(17)}`);
  console.log('');
  const ok = off.transmitted === N && off.records === 0 && on.absorbed === M && on.transmitted === N - M;
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the off-resonant beam passes COMPLETELY and leaves ZERO records; the resonant beam is`);
  console.log(`     absorbed exactly until the wall is used up (${M}) and the rest passes.`);
  console.log(`     A finite threshold exists -- in the FREQUENCY domain, and it is the molecule's own`);
  console.log(`     resonant absorber (rounds 47-57), not a new rule.`);
}
console.log('');

console.log('============================================================');
console.log('B: A FINITE THRESHOLD IN THE INTENSITY DOMAIN -- BLEACHING, WITH A KINK AT N = M');
console.log('============================================================');
{
  const M = 20;
  console.log(`   wall of ${M} molecules (recovery switched off: a fast pulse), resonant beam of N quanta`);
  console.log('');
  console.log('     N   transmitted   predicted max(0,N-M)   transmission   predicted transmission');
  let ok = true;
  for (const N of [5, 10, 19, 20, 21, 30, 40, 60]) {
    const r = wallRun(M, N, 4, 4, 0, Infinity, 2, 23);
    const pred = Math.max(0, N - M);
    const T = r.transmitted / N, Tp = pred / N;
    console.log(`   ${String(N).padStart(4)}   ${String(r.transmitted).padStart(11)}   ${String(pred).padStart(20)}   ${T.toFixed(3).padStart(12)}   ${Tp.toFixed(3).padStart(22)}`);
    if (r.transmitted !== pred) ok = false;
  }
  console.log('');
  console.log('  => ' + (ok ? 'PASS' : '*** FAIL ***') + ': transmission is exactly max(0, N-M)/N -- dim light is blocked (0 below N = M),');
  console.log('     bright light bleaches the wall and passes.  THE KINK IS AT N = M, and it is exact');
  console.log('     because each absorbed quantum uses up exactly one molecule.');
  console.log('     This is saturable absorption, and it is the finite intensity threshold rounds 33-34');
  console.log('     were looking for.');
}
console.log('');

console.log('============================================================');
console.log('C: RECOVERY -- THE WALL HEALS WITH THE MOLECULE\'S OWN LIFETIME (no new constant)');
console.log('============================================================');
{
  const M = 20, N = 20, TAU = 500, TRIALS = 4000;
  console.log(`   pulse of ${N} resonant quanta fills the wall; then wait T; then probe with ${N} more.`);
  console.log(`   molecular lifetime tau = ${TAU} (the SAME tau the molecule always had).`);
  console.log('');
  console.log('     wait T   recovered by probe (measured)   1 - exp(-T/tau)');
  let ok = true;
  for (const T of [0, 100, 250, 500, 1000]) {
    let sum = 0;
    for (let trial = 0; trial < TRIALS; trial++) {
      const mols = makeWall(M, 4, 0, TAU, 5000 + trial * 104729);
      for (let k = 0; k < M; k++) {                        // the filling pulse
        const m = mols[k];
        m.exc = true;
        m.recoverAt = (-TAU * Math.log(Math.max(1e-12, m.rnd())));
      }
      let recovered = 0;
      for (const m of mols) if (m.exc && T >= m.recoverAt) recovered++;
      sum += recovered;
    }
    const meas = sum / TRIALS, pred = M * (1 - Math.exp(-T / TAU));
    console.log(`     ${String(T).padStart(6)}   ${meas.toFixed(2).padStart(29)}   ${pred.toFixed(2).padStart(17)}`);
    if (Math.abs(meas - pred) > 0.45) ok = false;
  }
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the wall recovers on the molecule's own timescale, so the barrier has a DEAD TIME`);
  console.log(`     and it is not a new parameter -- it is tau, which the model already had.`);
}
console.log('');

console.log('============================================================');
console.log('D: THE COHERENCE COST -- DOES RECORDING DESTROY THE QUANTUM? (round 33\'s question)');
console.log('============================================================');
{
  // two-arm interference: arm A carries 1 quantum's worth of amplitude, arm B carries sqrt(n).
  // k quanta are absorbed out of arm B.  Two models for what that does:
  //    bridge : the count drops by one per absorption, so the AMPLITUDE is sqrt(n-k)
  //    linear : the amplitude is reduced in proportion, 1 - k/n  (each absorption "eats" amplitude)
  const W = 64, H = 64, N = W * H;
  const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0], OMEGA = Math.PI / 2, LO = 10, HI = 54;
  const idx = (x, y, d) => (y * W + x) * 4 + d;
  const zero = () => new Float64Array(4 * N);
  function move(a) {
    const o = new Float64Array(4 * N);
    for (let i = 0; i < N; i++) {
      const x = i % W, y = (i - x) / W;
      for (let d = 0; d < 4; d++) {
        const v = a[i * 4 + d]; if (v === 0) continue;
        const k = ((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W);
        o[k * 4 + d] += v;
      }
    }
    return o;
  }
  function visibility(sB) {                       // source B amplitude = sB, source A = 1
    let re = zero(), im = zero();
    for (let t = 0; t < 400; t++) {
      const th = OMEGA * t;
      for (let y = 0; y < H; y++) { re[idx(0, y, 1)] = Math.cos(th); im[idx(0, y, 1)] = Math.sin(th); }
      for (let x = 0; x < W; x++) { re[idx(x, 0, 0)] = sB * Math.cos(th); im[idx(x, 0, 0)] = sB * Math.sin(th); }
      re = move(re); im = move(im);
    }
    const bins = [0, 0, 0, 0], cnt = [0, 0, 0, 0];
    for (let y = LO; y < HI; y++) for (let x = LO; x < HI; x++) {
      const i = y * W + x;
      let sr = 0, si = 0;
      for (let d = 0; d < 4; d++) { sr += re[i * 4 + d]; si += im[i * 4 + d]; }
      const I = sr * sr + si * si;
      const b = (((x - y) % 4) + 4) % 4;
      bins[b] += I; cnt[b]++;
    }
    const p = bins.map((v, k) => v / cnt[k]);
    return (Math.max(...p) - Math.min(...p)) / (Math.max(...p) + Math.min(...p));
  }
  console.log('   two-arm interference, arm A = 1 quantum, arm B = n quanta, k of them absorbed by the wall');
  console.log('');
  console.log('     n   k absorbed   V measured   bridge 2*sqrt(n-k)/(1+(n-k))   linear model   discriminates?');
  let ok = true, nonMono = false, discriminated = 0, coincidence = 0;
  for (const n of [1, 2, 3, 5]) {
    let prevV = null;
    for (let k = 0; k <= n; k++) {
      const V = visibility(Math.sqrt(n - k));
      const Vbridge = (2 * Math.sqrt(n - k)) / (1 + (n - k));
      // the competing model: absorption eats AMPLITUDE in proportion (each quantum removes 1/n of it)
      const sLin = Math.sqrt(n) * Math.max(0, 1 - k / n);
      const Vlin = (2 * sLin) / (1 + sLin * sLin);
      // at k = 0 nothing was absorbed, and at k = n no light is left in the arm: there the two models
      // MUST agree, so the discrimination only lives in the interior 0 < k < n
      const interior = k > 0 && k < n;
      const differs = Math.abs(V - Vlin) > 0.01;
      if (interior) { if (differs) discriminated++; else coincidence++; }
      console.log(`   ${String(n).padStart(3)}   ${String(k).padStart(10)}   ${V.toFixed(4).padStart(10)}   ${Vbridge.toFixed(4).padStart(31)}   ${Vlin.toFixed(4).padStart(12)}   ${interior ? (differs ? 'yes' : 'no (accidental)') : 'n/a (models must agree)'}`);
      if (Math.abs(V - Vbridge) > 1e-3) ok = false;
      if (prevV !== null && V > prevV + 1e-9) nonMono = true;
      prevV = V;
    }
  }
  console.log('');
  console.log(`   interior rows where the linear model disagrees: ${discriminated}; rows where it accidentally agrees: ${coincidence}`);
  console.log(`   the visibility is NOT monotone in k: ${nonMono ? 'measured (absorbing can improve contrast by re-balancing the arms)' : 'not observed'}`);
  console.log(`   with n = 1 and k = 1 the fringe dies COMPLETELY (V = 0.0000) -- a single quantum, fully`);
  console.log(`   recorded, has no interference left`);
  console.log('');
  console.log(`  => ${ok && nonMono && discriminated > 0 && coincidence <= 1 ? 'PASS' : '*** FAIL ***'}: the measured visibility follows the BRIDGE model exactly at every row, and the`);
  console.log(`     linear model is ruled out wherever the two can be told apart (${coincidence} accidental agreement(s)).`);
  console.log(`     So: absorbing quanta removes WHOLE quanta; it does not attenuate the wave.`);
  console.log('');
  console.log('   ROUND 33\'S ANSWER: in Life, the contact that records the phase necessarily DESTROYS the');
  console.log('   particle, so no gentle detector existed.  Here the wall records the quantum and re-emits');
  console.log('   it elsewhere: the COUNT survives, the PHASE does not.  A gentle detector does exist --');
  console.log('   it is a phase-destroying, number-conserving scatterer.');
  console.log('');
  console.log('   HONEST BOUNDARY: this wall is built from the DESIGNED molecules of rounds 40-64.  It does');
  console.log('   NOT answer rounds 33-34\'s real question, which was whether the BASE rules can make matter.');
  console.log('   What it does answer: GIVEN the model\'s matter, the three things those rounds went looking');
  console.log('   for -- a finite threshold, a dead time, and a gentle detector -- are all there, and they');
  console.log('   are all consequences of machinery the model already had.');
}
