'use strict';
/*
 * nematic_substrate.js  --  ROUND 124 (D8'''): WHAT E4 ACTUALLY NEEDS, AND WHETHER A DIRECTOR CAN SUPPLY IT
 * ---------------------------------------------------------------------------
 * ROUND 123 SAID something too strong, and I am correcting it here first.  It showed that a uniform direction
 * in three dimensions, read at half its polar angle, reproduces the E4 measure |sin 2 phi| / 4.  That is
 * true.  But the S = 2.8284 it reported came from TWO ingredients, and the report attributed the whole
 * effect to the measure.  The second ingredient is the (a+b)/2 shift in the effective angle:
 *
 *      theta = phi + (a+b)/2
 *
 * which makes each wing's outcome depend on BOTH settings -- measurement dependence, the 0.201786 bits of
 * round 108.  WITHOUT it the model is a standard local hidden-variable model, and then NO distribution can
 * exceed 2, by Bell's theorem.  This script measures that, and then asks the real question.
 *
 * THE REAL QUESTION: a nematic (director) field is orientation defined only up to a sign, so
 *   - the half-angle (spinor-like) reading is NATIVE to it, and
 *   - a nematic couples to a field through cos^2, i.e. with period pi, so its response to TWO applied fields
 *     is naturally governed by their AVERAGE.
 * If both are native, then the E4 mechanism has, for the first time, a classical substrate in which neither
 * ingredient is put in by hand.  This round checks the second one directly by minimising the director energy.
 *
 * PRE-REGISTERED PREDICTIONS
 *   P1 THE NECESSITY TABLE: of the four combinations {shift, no shift} x {uniform, E4 measure}, only
 *      (shift + E4 measure) exceeds 2.  The other three are at or below 2.
 *   P2 in particular the E4 measure WITHOUT the shift is at or below 2 -- Bell's theorem, and the correction
 *      of round 123's overstatement.
 *   P3 THE DIRECTOR: minimising F(theta) = -A[cos^2(theta - a) + cos^2(theta - b)] puts the director at the
 *      AVERAGE (a+b)/2 modulo pi/2, so the measurement-dependence structure is NATIVE to a nematic.
 *   P4 therefore a 3D nematic (half angle for the measure) that responds to both analysers (average for the
 *      shift) has both ingredients without hand-insertion.  THAT is the next substrate to build, and this
 *      round states the prediction for it so that it can fail.
 * ---------------------------------------------------------------------------
 */

const SQ2 = Math.SQRT2, PI = Math.PI;
const sgn = (x) => (x >= 0 ? 1 : -1);
const SETTINGS = [[0, PI / 4], [0, -PI / 4], [PI / 2, PI / 4], [PI / 2, -PI / 4]];
const SfromE = (e) => Math.abs(e[0] + e[1] + e[2] - e[3]);
function stream(seed) { let s = seed | 0; if (s === 0) s = 0x9e3779b9; return () => { s ^= s << 13; s ^= s >>> 17; s ^= s << 5; return (s >>> 0) / 4294967296; }; }
const uniQ = () => 1 / (2 * PI);
const e4Q = (phi) => Math.abs(Math.sin(2 * phi)) / 4;

// E(a,b) with a switch for the (a+b)/2 shift; no shift means lambda alone fixes the outcome => local model
function Eof(q, a, b, shift, nPhi) {
  const sc = shift ? (a + b) / 2 : 0;
  let acc = 0, w = 0;
  for (let i = 0; i < nPhi; i++) {
    const phi = (i + 0.5) * (2 * PI) / nPhi, wgt = q(phi) * (2 * PI) / nPhi;
    const th = phi + sc;
    acc += wgt * sgn(Math.cos(th - a)) * -sgn(Math.cos(th - b));
    w += wgt;
  }
  return acc / w;
}
const Sof = (q, shift, nPhi) => SfromE(SETTINGS.map(([a, b]) => Eof(q, a, b, shift, nPhi)));

console.log('============================================================');
console.log('PREDICTIONS');
console.log('============================================================');
console.log('   P1 only (shift + E4 measure) exceeds 2; the other three combinations do not');
console.log('   P2 the E4 measure WITHOUT the shift is at or below 2 (Bell), which corrects round 123');
console.log('   P3 a director minimising -A[cos^2(theta-a) + cos^2(theta-b)] sits at the AVERAGE of the two');
console.log('      settings, so the measurement-dependence structure is native to a nematic');
console.log('   P4 a 3D nematic therefore has BOTH ingredients natively -- stated as a prediction to be tested');
console.log('');

console.log('============================================================');
console.log('PART A: THE NECESSITY TABLE (P1, P2)');
console.log('============================================================');
{
  console.log('   measure            (a+b)/2 shift    S           verdict');
  const rows = [];
  for (const [name, q] of [['uniform', uniQ], ['E4 |sin 2phi|', e4Q]]) {
    for (const shift of [false, true]) {
      const S = Sof(q, shift, 4000);
      rows.push({ name, shift, S });
      const verdict = S > 2.01 ? 'EXCEEDS 2' : 'at or below 2';
      console.log(`   ${name.padEnd(18)} ${shift ? 'YES' : ' no'}             ${S.toFixed(4)}     ${verdict}`);
    }
  }
  console.log('');
  const only = rows.filter((r) => r.S > 2.01);
  console.log(`   P1 ${only.length === 1 && only[0].shift && only[0].name.startsWith('E4') ? 'CONFIRMED' : 'NOT confirmed'}: exactly one of the four exceeds 2, namely the E4 measure WITH the shift`);
  const noShiftE4 = rows.find((r) => r.name.startsWith('E4') && !r.shift);
  console.log(`   P2 ${noShiftE4.S <= 2.01 ? 'CONFIRMED' : 'NOT confirmed'}: the E4 measure without the shift gives S = ${noShiftE4.S.toFixed(4)}`);
  console.log('');
  console.log('   THE CORRECTION OF ROUND 123: that round reported S = 2.8284 for the half-angle reading and');
  console.log('   attributed it to the measure.  It comes from the measure AND the shift together.  Without the');
  console.log('   shift no distribution can help, which is Bell\'s theorem and was already measured in round 72.');
  console.log('   So the measure is necessary but nowhere near sufficient.');
}
console.log('');

console.log('============================================================');
console.log('PART B: THE DIRECTOR RESPONDS TO THE AVERAGE (P3)');
console.log('============================================================');
{
  // a nematic: F(theta) = -A [cos^2(theta - a) + cos^2(theta - b)], period pi in theta.
  // Algebraically F = const - A cos(2 theta - (a+b)) cos(a - b), so the minimum is at theta = (a+b)/2 mod pi
  // -- UNLESS cos(a-b) = 0, i.e. unless the two settings differ by 90 degrees, in which case the two fields
  // are the SAME director orientation and F is completely flat.  My first version compared angles without
  // excluding that degenerate case, so three rows looked like failures of the prediction when in fact the
  // answer was undefined there.  The degenerate rows are now reported as flat, and every other row is
  // compared against the analytic value.
  const A = 1.0;
  const F = (th, a, b) => -A * (Math.cos(th - a) ** 2 + Math.cos(th - b) ** 2);
  const minimise = (a, b) => {
    let best = 0, bestF = Infinity;   // MINIMISE the energy (my first version kept the maximum)
    for (let k = 0; k < 20000; k++) {
      const th = PI * k / 20000;
      const f = F(th, a, b);
      if (f < bestF) { bestF = f; best = th; }
    }
    return { th: best, spread: F(best, a, b) - F(best + PI / 2, a, b) };
  };
  const mod = (x, p) => ((x % p) + p) % p;
  console.log('   a (deg)   b (deg)   |a-b| mod 180   theta* mod pi   analytic (a+b)/2 mod pi   offset       note');
  let worst = 0, tested = 0;
  for (const [aDeg, bDeg] of [[0, 45], [0, 90], [30, 120], [-20, 70], [10, 10], [0, 180], [15, 75], [-40, 25], [5, 200]]) {
    const a = aDeg * PI / 180, b = bDeg * PI / 180;
    const r = minimise(a, b);
    const avg = (a + b) / 2;
    const diff = mod(Math.abs(aDeg - bDeg), 180);
    const degenerate = Math.abs(diff - 90) < 1e-9;
    // CORRECT CRITERION, from the algebra: F = const - A cos(2 theta - (a+b)) cos(a-b), so the minimiser
    // satisfies 2 theta* = (a+b) + [pi if cos(a-b) < 0] modulo 2 pi.  Comparing raw (a+b)/2 is wrong
    // because a director does not distinguish a from a+pi: for (a,b) = (5,200) the two fields are at 5
    // and 20 degrees, and the average is 12.5, which is exactly where the director went.
    const flip = Math.cos(a - b) < 0 ? PI : 0;
    let d2 = mod(2 * r.th - (a + b) - flip, 2 * PI);
    if (d2 > PI) d2 -= 2 * PI;
    const off = d2 / 2;
    const note = degenerate ? 'DEGENERATE: the two fields are the same director orientation, F is flat' : '';
    if (!degenerate) { worst = Math.max(worst, Math.abs(off)); tested++; }
    console.log(`   ${String(aDeg).padStart(5)}     ${String(bDeg).padStart(5)}        ${diff.toFixed(1).padStart(6)}      ${(mod(r.th, PI) * 180 / PI).toFixed(2).padStart(7)}          ${(mod(avg, PI) * 180 / PI).toFixed(2).padStart(7)}          ${(off * 180 / PI).toFixed(3).padStart(8)}    ${note}`);
  }
  console.log('');
  console.log(`   P3 ${worst < 0.02 ? 'CONFIRMED' : 'NOT confirmed'}: over the ${tested} non-degenerate rows the largest offset from (a+b)/2 is ${(worst * 180 / PI).toFixed(4)} degrees`);
  console.log('   (modulo pi, which is all a director defines).  So the AVERAGE structure is native to a nematic, and');
  console.log('   the energy period pi in theta is the half-angle nature stated in part C.');
}