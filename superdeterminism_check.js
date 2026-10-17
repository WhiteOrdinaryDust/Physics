'use strict';
/*
 * superdeterminism_check.js
 * ---------------------------------------------------------------------------
 * Does accepting SUPERDETERMINISM give you 2*sqrt(2)?
 *
 * THE CLAIM UNDER TEST
 *   "Observers and particles both emerge from the same deterministic substrate,
 *    so the initial condition determines everything downstream -- and therefore
 *    statistical independence fails, and therefore Bell's bound does not apply."
 *
 * WHAT BELL'S THEOREM ACTUALLY ASSUMES
 *        locality  +  determinism  +  STATISTICAL INDEPENDENCE   =>   CHSH <= 2
 *   Dropping statistical independence does remove the bound. So the first half of
 *   the claim is correct.
 *
 * THE QUESTION THIS FILE SETTLES
 *   Once statistical independence is dropped, what IS the maximum? Is it 2*sqrt(2)
 *   (Tsirelson, the quantum value) or 4 (the no-signalling / PR-box value)?
 *
 * MODEL
 *   A superdeterministic model: a hidden variable lambda, plus the measurement
 *   settings (A,B), whose joint statistics are NOT independent of lambda. For each
 *   (lambda, A, B) the outcomes (a,b) are deterministic. We search over all such
 *   deterministic assignments and report the largest achievable CHSH.
 *
 *   NOTE: no-signalling is NOT imposed, because in a superdeterministic model the
 *   settings are correlated with lambda by construction, so the operational
 *   distribution P(a,b|x,y) is not required to be no-signalling.
 * ---------------------------------------------------------------------------
 */

function chshFromTable(P) {
  // P[x][y][a][b]
  const E = [[0, 0], [0, 0]];
  for (let x = 0; x < 2; x++) for (let y = 0; y < 2; y++) {
    let same = 0, diff = 0;
    for (let a = 0; a < 2; a++) for (let b = 0; b < 2; b++) {
      if (a === b) same += P[x][y][a][b]; else diff += P[x][y][a][b];
    }
    E[x][y] = same - diff;
  }
  return E[0][0] + E[0][1] + E[1][0] - E[1][1];
}
function noSignallingDev(P) {
  let dev = 0;
  for (let x = 0; x < 2; x++) {
    let p0y0 = 0, p0y1 = 0;
    for (let b = 0; b < 2; b++) { p0y0 += P[x][0][0][b]; p0y1 += P[x][1][0][b]; }
    dev = Math.max(dev, Math.abs(p0y0 - p0y1));
  }
  for (let y = 0; y < 2; y++) {
    let q0x0 = 0, q0x1 = 0;
    for (let a = 0; a < 2; a++) { q0x0 += P[0][y][a][0]; q0x1 += P[1][y][a][0]; }
    dev = Math.max(dev, Math.abs(q0x0 - q0x1));
  }
  return dev;
}

console.log('=== does superdeterminism give 2*sqrt(2)? ===\n');

/*
 * EXHAUSTIVE SEARCH over superdeterministic models.
 *
 * lambda takes N values with arbitrary weights. For each lambda we enumerate the
 * 16 possible deterministic response tables (a,b) as a function of (x,y), stored
 * as a 4-bit mask. The model is fully specified by N such masks plus the weights.
 *
 * Because CHSH is linear in the weights, the maximum over weights is attained at
 * a vertex: all weight on ONE lambda. That is the pure deterministic case, which
 * obviously cannot exceed 2. The interesting regime is therefore a MIXTURE over
 * lambda with DIFFERENT response tables -- and crucially, we let the lambda
 * distribution depend on the settings, which is exactly the superdeterministic
 * move.
 *
 * Concretely we search: for each pair of settings (x,y), choose a distribution
 * over lambda (equivalently, choose for each (x,y) one of the 16 response masks
 * and a weight). This is the most general superdeterministic model at finite N.
 */
function searchSuperdeterministic(N) {
  // response mask m: bit index (x*2+y) gives (a,b) as (a<<1)|b
  const masks = [];
  for (let m = 0; m < 16; m++) {
    // table[x][y] = [a,b]
    const t = [[[0, 0], [0, 0]], [[0, 0], [0, 0]]];
    for (let x = 0; x < 2; x++) for (let y = 0; y < 2; y++) {
      const bit = (m >> (x * 2 + y)) & 1;
      // use a second bit for b: combine two masks, see below
      t[x][y] = [bit, 0];
    }
    masks.push(t);
  }
  // The full (a,b) needs 4 bits: one for each of the 4 (x,y) pairs for a, and
  // one for each for b. So a response table is a mask in 0..255.
  let best = -Infinity, bestDesc = null, bestNS = -Infinity;
  const totalMasks = 256;
  // choose one table per (x,y) -> 256^4 is too many; instead choose N tables with
  // weights, which for N=4 subsumes the per-(x,y) choice. Enumerate N=4 with
  // all-equal weights to keep it finite: 256^4 is 4.3e9, too many.
  // Reduce: we only need, for each (x,y), the ability to pick ANY (a,b) pair.
  // Since lambda can be chosen per (x,y) in the superdeterministic model, the
  // effective freedom is: for each (x,y) pick (a,b) in {0,1}^2, i.e. 4 options.
  // That gives 4^4 = 256 models -- and it is the maximal freedom available.
  for (let code = 0; code < 256; code++) {
    const P = [
      [[[0, 0], [0, 0]], [[0, 0], [0, 0]]],
      [[[0, 0], [0, 0]], [[0, 0], [0, 0]]],
    ];
    for (let xy = 0; xy < 4; xy++) {
      const x = xy >> 1, y = xy & 1;
      const choice = (code >> (xy * 2)) & 3;   // 0..3 = (a,b)
      const a = choice >> 1, b = choice & 1;
      P[x][y][a][b] += 1;
    }
    const S = chshFromTable(P);
    if (S > best) { best = S; bestDesc = code; }
    if (noSignallingDev(P) < 1e-12 && S > bestNS) bestNS = S;
  }
  return { best, bestNS, bestDesc };
}

console.log('--- exhaustive search over superdeterministic models ---');
console.log('  (settings may be correlated with the hidden variable; no-signalling NOT imposed)');
{
  const r = searchSuperdeterministic(4);
  console.log(`  max CHSH                                    : ${r.best.toFixed(6)}`);
  console.log(`  max CHSH among no-signalling ones           : ${r.bestNS.toFixed(6)}`);
  console.log(`  Bell bound (statistical independence)       : 2.000000`);
  console.log(`  Tsirelson bound (quantum mechanics)         : ${(2 * Math.SQRT2).toFixed(6)}`);
  console.log(`  PR box / no-signalling limit                : 4.000000`);
  console.log('');
  console.log(`  => dropping statistical independence lifts the bound from 2 to ${r.best.toFixed(0)},`);
  console.log('     NOT to 2*sqrt(2).');
}

// ---- now the quantum comparison, for contrast
console.log('\n--- for contrast: what gives exactly 2*sqrt(2)? ---');
{
  // |Phi+> correlation matrix is diag(1,-1,1); Horodecki gives 2*sqrt(2)
  const t1 = 1, t2 = 1;
  console.log(`  two-qubit maximally entangled state, correlation singular values 1,1`);
  console.log(`  Horodecki bound 2*sqrt(t1^2+t2^2) = ${(2 * Math.sqrt(t1 * t1 + t2 * t2)).toFixed(6)}`);
  console.log('  This value comes from the ALGEBRA of the tensor product and of local');
  console.log('  unitary rotations. There is no analogue of that algebra in a classical');
  console.log('  deterministic dynamics: classical configurations do not superpose with');
  console.log('  complex amplitudes, so there is nothing for a "local rotation" to act on.');
}

console.log('\n=== conclusion ===');
console.log('  1. CORRECT: in a substrate where observers and particles both emerge,');
console.log('     statistical independence fails automatically. Bell\'s bound therefore');
console.log('     does not apply, and CHSH > 2 becomes POSSIBLE. No extra assumption');
console.log('     is needed for this step, and it is not a conspiracy.');
console.log('');
console.log('  2. INCOMPLETE: removing statistical independence does not produce');
console.log('     2*sqrt(2). It removes the ceiling entirely, up to the no-signalling');
console.log('     limit of 4. The value 2*sqrt(2) is a separate, additional fact that');
console.log('     comes from complex amplitudes and tensor products.');
console.log('');
console.log('  So superdeterminism is NECESSARY for a deterministic substrate to violate');
console.log('  Bell, but NOT SUFFICIENT to reproduce quantum mechanics. The remaining');
console.log('  gap is the same complex-amplitude gap identified four times earlier in');
console.log('  this project -- superdeterminism removes the obstacle, it does not supply');
console.log('  the missing structure.');
console.log('');
console.log('=== what this changes, and what it does not ===');
console.log('  CHANGES  : the goal is now well posed. We no longer need to explain how a');
console.log('             deterministic substrate avoids Bell\'s bound -- it simply does not');
console.log('             obey the hypothesis. The task is to derive WHICH structure gives');
console.log('             2*sqrt(2) rather than 4.');
console.log('  DOES NOT : supply that structure. And a superdeterministic model is at');
console.log('             risk of unfalsifiability, because any mismatch can in principle');
console.log('             be blamed on a setting-system correlation. The safety valve is');
console.log('             that the model must PREDICT the value, not merely permit it.');
