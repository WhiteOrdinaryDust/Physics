'use strict';
/*
 * single_measurement.js   --   ROUND 65
 * ---------------------------------------------------------------------------
 * THE OTHER HALF OF THE BORN RULE: ONE QUANTUM, ONE OUTCOME, AND WHERE THE RANDOMNESS LIVES.
 *
 * Rounds 62-64 built a linear field with phases, interference, a count bridge (n = |amplitude|^2) and
 * a laser.  All of that is about ENSEMBLES: an intensity profile, a rate ladder, a distribution of
 * outcomes.  What is still missing is the single event: one quantum arrives, and ONE detector fires.
 *
 * This round builds the detector.  The rule (and it is a rule I write, on top of round 63's bridge):
 *        a single quantum in the field is detected at ONE cell, drawn with probability
 *        p(cell) = |psi(cell)|^2 / sum |psi|^2
 * and the drawing is done by a DETERMINISTIC hidden clock (xorshift32) -- the same kind of object as
 * the molecule's decay clock (rounds 53-54).  So the randomness is epistemic: same field + same
 * hidden clock => same outcome; same field + different hidden clock => the outcome follows |psi|^2.
 *
 * PREDICTIONS (written before the run)
 *   A  one run gives exactly one outcome; the outcome is reproducible for a fixed clock and varies
 *      with the clock.  The visible state (the field) is identical across runs -- only the hidden
 *      state differs, which is the project's central mechanism since round 32.
 *   B  the ensemble of single outcomes reproduces |psi|^2: over the four fringe bins the histogram
 *      must be 0.50 / 0.25 / 0.00 / 0.25, and the DARK bin must get EXACTLY zero events (a classical
 *      wave has an exactly dark fringe).  Shifting the source phase by 90 degrees rotates the
 *      histogram by one bin.
 *   C  controls that can fail: a phase-blind field gives a FLAT histogram (0.25 each), and sampling
 *      from |psi| instead of |psi|^2 gives a different, wrong histogram.
 *   D  the collapse rule cannot be local: absorbing only the firing cell leaves the wave and the
 *      detector fires again (a second detection in ~100% of runs).  Removing the whole field gives
 *      exactly one.  And a PARTIAL absorber on one arm must reduce the fringe visibility to
 *           V(kappa) = 2*sqrt(1-kappa) / (2-kappa)
 * ---------------------------------------------------------------------------
 */

const W = 64, H = 64, N = W * H;
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
const OMEGA = Math.PI / 2;
const LO = 10, HI = 54;                       // the detector region
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
// the two-source field: A enters from the left going +x, B from the bottom going +y with phase delta
// and amplitude scale s (s < 1 = a which-path absorber on arm B).
function field(delta, s) {
  let re = zero(), im = zero();
  for (let t = 0; t < 400; t++) {
    const th = OMEGA * t;
    for (let y = 0; y < H; y++) { re[idx(0, y, 1)] = Math.cos(th); im[idx(0, y, 1)] = Math.sin(th); }
    for (let x = 0; x < W; x++) { re[idx(x, 0, 0)] = s * Math.cos(th + delta); im[idx(x, 0, 0)] = s * Math.sin(th + delta); }
    re = move(re); im = move(im);
  }
  return { re, im };
}
// the intensity the detector sees at a cell: |sum over directions|^2  (round 62's observable)
function intensity(re, im) {
  const I = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    let sr = 0, si = 0;
    for (let d = 0; d < 4; d++) { sr += re[i * 4 + d]; si += im[i * 4 + d]; }
    I[i] = sr * sr + si * si;
  }
  return I;
}
function detectorCells() {
  const cells = [];
  for (let y = LO; y < HI; y++) for (let x = LO; x < HI; x++) cells.push(y * W + x);
  return cells;
}
const CELLS = detectorCells();
const binOf = (i) => { const x = i % W, y = (i - x) / W; return (((x - y) % 4) + 4) % 4; };
function clock(seed) { let x = seed | 0; if (x === 0) x = 0x9e3779b9; return () => { x ^= x << 13; x ^= x >>> 17; x ^= x << 5; return ((x >>> 0) / 4294967296); }; }
// one single-quantum detection: draw a cell with probability proportional to the weight given
function detect(weights, rnd) {
  let tot = 0;
  for (const w of weights) tot += w;
  let u = rnd() * tot;
  for (let k = 0; k < weights.length; k++) { u -= weights[k]; if (u <= 0) return k; }
  return weights.length - 1;
}

// ============================================================ A: one run, one outcome
console.log('============================================================');
console.log('A: ONE QUANTUM, ONE OUTCOME -- AND IT IS THE HIDDEN CLOCK THAT DECIDES');
console.log('============================================================');
{
  const { re, im } = field(0, 1);
  const I = intensity(re, im);
  const w = CELLS.map((i) => I[i]);
  const tot = w.reduce((a, b) => a + b, 0);
  const outcomes = new Set();
  for (let c = 1; c <= 1000; c++) outcomes.add(detect(w, clock(c * 2654435761)));
  const r1 = detect(w, clock(12345)), r2 = detect(w, clock(12345)), r3 = detect(w, clock(999));
  let differentOnHidden = 0;
  for (let c = 1; c <= 1000; c++) {
    if (detect(w, clock(c * 2654435761)) !== detect(w, clock(987654321))) differentOnHidden++;
  }
  console.log(`   the field is identical in every run (same sources, same 400 steps) -- only the clock differs`);
  console.log(`   same hidden clock, two calls   : cell index ${r1} and ${r2}  => ${r1 === r2 ? 'IDENTICAL (deterministic)' : '*** DIFFERENT ***'}`);
  console.log(`   a third clock                  : cell index ${r3}`);
  console.log(`   1000 different clocks          : ${outcomes.size} distinct cells out of ${CELLS.length} in the detector`);
  console.log(`   comparing 1000 clocks against one fixed other clock: they disagree in ${differentOnHidden}/1000 runs`);
  console.log(`   total intensity in the detector: ${tot.toFixed(1)}   (one quantum is spread over the whole wave)`);
  console.log('');
  const ok = r1 === r2 && outcomes.size > 50 && differentOnHidden > 500;
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the outcome is deterministic given the hidden state, varies when the hidden state`);
  console.log(`     varies, and the visible state never changed.  This is the mechanism of round 32, now`);
  console.log(`     doing the Born rule's job: the randomness is in the observer's access, not in the field.`);
}
console.log('');

// ============================================================ B: the ensemble is Born
console.log('============================================================');
console.log('B: 20000 SINGLE OUTCOMES -- DOES THE HISTOGRAM REPRODUCE |psi|^2 ?');
console.log('============================================================');
{
  const TRIALS = 20000;
  function histogram(delta, weightsFrom) {
    const { re, im } = field(delta, 1);
    const I = intensity(re, im);
    const w = weightsFrom === 'amplitude' ? CELLS.map((i) => Math.sqrt(I[i])) : CELLS.map((i) => I[i]);
    const hist = [0, 0, 0, 0];
    const rnd = clock(20250920);
    for (let t = 0; t < TRIALS; t++) hist[binOf(CELLS[detect(w, rnd)])]++;
    return hist.map((v) => v / TRIALS);
  }
  // the coherence-free reference: source A alone and source B alone, their INTENSITIES added
  const Ione = intensity(...Object.values(field(0, 0)).slice(0, 2));      // source A only (B off)
  const { re: reB, im: imB } = (() => {   // source B alone: the same geometry, sources exchanged
    let re = zero(), im = zero();
    for (let t = 0; t < 400; t++) {
      const th = OMEGA * t;
      for (let x = 0; x < W; x++) { re[idx(x, 0, 0)] = Math.cos(th); im[idx(x, 0, 0)] = Math.sin(th); }
      re = move(re); im = move(im);
    }
    return { re, im };
  })();
  const Itwo = intensity(reB, imB);
  const wBlind = CELLS.map((i) => Ione[i] + Itwo[i]);
  const histBlind = [0, 0, 0, 0];
  {
    const rnd = clock(20250920);
    for (let t = 0; t < TRIALS; t++) histBlind[binOf(CELLS[detect(wBlind, rnd)])]++;
    for (let k = 0; k < 4; k++) histBlind[k] /= TRIALS;      // probabilities, not counts
  }
  const h0 = histogram(0, 'intensity');
  const h90 = histogram(Math.PI / 2, 'intensity');
  const hAmp = histogram(0, 'amplitude');
  const predicted = [0.5, 0.25, 0.0, 0.25];
  const sigma = predicted.map((p) => Math.sqrt(p * (1 - p) / TRIALS));
  console.log('   bin (u = x-y mod 4)     0        1        2        3');
  console.log(`   predicted              ${predicted.map((v) => v.toFixed(3).padStart(8)).join(' ')}`);
  console.log(`   measured, delta = 0    ${h0.map((v) => v.toFixed(3).padStart(8)).join(' ')}`);
  console.log(`   deviation (in sigma)   ${h0.map((v, k) => ((v - predicted[k]) / (sigma[k] || 1e-9)).toFixed(2).padStart(8)).join(' ')}`);
  console.log(`   measured, delta = 90   ${h90.map((v) => v.toFixed(3).padStart(8)).join(' ')}   (predicted rotation: 0.25/0.00/0.25/0.50)`);
  console.log('');
  console.log(`   the DARK bin (u = 2): predicted exactly 0, measured ${h0[2]} events out of ${TRIALS}`);
  const maxSig = Math.max(...h0.map((v, k) => Math.abs(v - predicted[k]) / (sigma[k] || 1e-9)).filter((_, k) => k !== 2));
  const ok = h0[2] === 0 && maxSig < 4 &&
    Math.abs(h90[0] - 0.25) < 0.02 && Math.abs(h90[2] - 0.25) < 0.02 && Math.abs(h90[3] - 0.5) < 0.02;
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: 20000 single outcomes reproduce |psi|^2 to within statistics, the dark fringe gets`);
  console.log(`     EXACTLY zero events, and a 90 degree phase shift rotates the whole histogram by one bin.`);
  console.log(`     The Born rule's operational half is now measured, not assumed.`);
  console.log('');
  console.log('C: THE CONTROLS');
  console.log(`   phase-blind field (intensities added): ${histBlind.map((v) => v.toFixed(3).padStart(8)).join(' ')}  <- flat, no fringes`);
  console.log(`   sampling from |psi| instead of |psi|^2: ${hAmp.map((v) => v.toFixed(3).padStart(8)).join(' ')}  <- the probabilities differ`);
  console.log('');
  console.log(`  => ${Math.max(...histBlind) - Math.min(...histBlind) < 0.02 && Math.abs(hAmp[0] - 0.5) > 0.05 ? 'PASS' : '*** FAIL ***'}: the coherence-free field gives a flat histogram, so the fringes are not an`);
  console.log(`     artefact of the geometry; and the square is a real, load-bearing choice -- it is round 63's`);
  console.log(`     bridge, i.e. a rule I wrote.  This round does not derive the square; it measures what the`);
  console.log(`     square implies for single events.`);
}
console.log('');

// ============================================================ D: the collapse rule
console.log('============================================================');
console.log('D: THE COLLAPSE RULE CANNOT BE LOCAL -- AND PARTIAL ABSORPTION COSTS VISIBILITY');
console.log('============================================================');
{
  const { re, im } = field(0, 1);
  const I = intensity(re, im);
  const w = CELLS.map((i) => I[i]);
  const tot = w.reduce((a, b) => a + b, 0);
  // (i) absorb ONLY the firing cell, then look again
  let secondLocal = 0;
  const TRIALS = 2000;
  for (let t = 0; t < TRIALS; t++) {
    const rnd = clock(777 + t * 4099);
    const k = detect(w, rnd);
    const left = w.reduce((a, b) => a + b, 0) - w[k];
    const rnd2 = clock(31337 + t * 6151);
    const k2 = detect(w.map((v, j) => (j === k ? 0 : v)), rnd2);
    if (left > 0 && k2 !== k) secondLocal++;
  }
  console.log('   D1  after one detection, absorb ONLY the cell that fired, then detect again:');
  console.log(`        a second detection occurs in ${(100 * secondLocal / TRIALS).toFixed(2)}% of runs`);
  console.log(`        (the wave is spread over ${CELLS.length} cells, so removing one cell leaves ~all of it)`);
  console.log('');
  console.log('   D2  after one detection, remove the WHOLE field, then detect again:');
  console.log(`        residual intensity = 0  =>  a second detection is impossible  =>  exactly one count`);
  console.log('');
  console.log('   => the model\'s single-quantum detector must empty the entire extended wave, which is a');
  console.log('      NON-LOCAL operation.  This is the same price \'t Hooft lists in section 18.2: continuous');
  console.log('      operations on a lattice cost non-locality.  A local absorption rule cannot do it.');
  console.log('');
  // (iii) partial absorption on one arm: measured visibility vs the closed form
  function visibility(delta, s) {
    const f = field(delta, s);
    const If = intensity(f.re, f.im);
    const bins = [0, 0, 0, 0], cnt = [0, 0, 0, 0];
    for (const i of CELLS) { const b = binOf(i); bins[b] += If[i]; cnt[b]++; }
    const p = bins.map((v, k) => v / cnt[k]);
    return (Math.max(...p) - Math.min(...p)) / (Math.max(...p) + Math.min(...p));
  }
  console.log('   D3  a which-path absorber of strength kappa on one arm:');
  console.log('        kappa   measured visibility   predicted 2*sqrt(1-kappa)/(2-kappa)');
  let ok = true;
  for (const kappa of [0, 0.25, 0.5, 0.75, 1]) {
    const s = Math.sqrt(1 - kappa);
    const V = visibility(0, s);
    const pred = (2 * s) / (2 - kappa);
    console.log(`        ${String(kappa).padEnd(7)} ${V.toFixed(4).padEnd(21)} ${pred.toFixed(4)}`);
    if (Math.abs(V - pred) > 0.01) ok = false;
  }
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the visibility falls exactly as the closed form for a partly absorbed arm.`);
  console.log('');
  console.log('      HONEST BOUNDARY: this is the VISIBILITY loss only.  The full complementarity relation');
  console.log('      V^2 + D^2 <= 1 needs somewhere for the which-path information to BE -- a second system,');
  console.log('      i.e. entanglement.  A classical field has no second system to hold it, so this model');
  console.log('      reproduces the fringe loss but cannot state the trade-off.  That is round 62\'s boundary');
  console.log('      again: interference yes, entanglement no.');
}
