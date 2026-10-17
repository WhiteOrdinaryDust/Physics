'use strict';
/*
 * real_vs_complex.js
 * ---------------------------------------------------------------------------
 * Does 2*sqrt(2) actually require COMPLEX amplitudes, or does REAL suffice?
 *
 * I have said several times in this project that "2*sqrt(2) needs complex
 * amplitudes". This file tests that claim directly, because if it is wrong the
 * whole framing of "the missing ingredient is the complex structure" needs to be
 * restated.
 *
 * The test:
 *   REAL model     : the two-particle state lives in R^2 (x) R^2, local
 *                    operations are real rotations SO(2) of the measurement
 *                    direction in the (x,z) plane.
 *   COMPLEX model  : the state lives in C^2 (x) C^2, local operations are the
 *                    full SU(2), so the measurement direction can be rotated out
 *                    of the (x,z) plane, into the y direction.
 *
 * For each, compute the maximum CHSH over measurement choices.
 *
 * If the real model also reaches 2*sqrt(2), then complex structure is NOT needed
 * for the CHSH bound, and the "missing ingredient" is a different, weaker thing.
 * ---------------------------------------------------------------------------
 */

// ---- build a pure two-qubit state from amplitudes in a real basis |00>,|01>,|10>,|11>
// We allow complex amplitudes so the same machinery serves both models.
function bellState(complex) {
  // (|00> + |11>)/sqrt2 : real amplitudes, so this state is available in both models
  return [[1 / Math.SQRT2, 0], [0, 0], [0, 0], [1 / Math.SQRT2, 0]];
}

// correlation E = <psi| A (x) B |psi>, where A, B are real 2x2 observables
// given as the Bloch vector v: O = v . sigma  (real symmetric for v in xz-plane)
function blochToMatrix(v) {
  // v = [vx, vy, vz] -> vx*X + vy*Y + vz*Z  (complex 2x2)
  const { } = {};
  const X = [[[0, 0], [1, 0]], [[1, 0], [0, 0]]];
  const Y = [[[0, 0], [0, -1]], [[0, 1], [0, 0]]];
  const Z = [[[1, 0], [0, 0]], [[0, 0], [-1, 0]]];
  const out = [[[0, 0], [0, 0]], [[0, 0], [0, 0]]];
  const comps = [v[0], v[1], v[2]];
  const mats = [X, Y, Z];
  for (let m = 0; m < 3; m++) {
    for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) {
      out[i][j][0] += comps[m] * mats[m][i][j][0];
      out[i][j][1] += comps[m] * mats[m][i][j][1];
    }
  }
  return out;
}

// <psi| (A (x) B) |psi>
function expectation(psi, A, B) {
  // apply (A (x) B) to psi
  const out = Array.from({ length: 4 }, () => [0, 0]);
  for (let a = 0; a < 2; a++) for (let b = 0; b < 2; b++) {
    let re = 0, im = 0;
    for (let ap = 0; ap < 2; ap++) for (let bp = 0; bp < 2; bp++) {
      const x = A[a][ap], y = B[b][bp], z = psi[ap * 2 + bp];
      const xyRe = x[0] * y[0] - x[1] * y[1], xyIm = x[0] * y[1] + x[1] * y[0];
      re += xyRe * z[0] - xyIm * z[1];
      im += xyRe * z[1] + xyIm * z[0];
    }
    out[a * 2 + b] = [re, im];
  }
  let re = 0;
  for (let i = 0; i < 4; i++) re += psi[i][0] * out[i][0] + psi[i][1] * out[i][1];
  return re;
}

// ---- measurement directions
// REAL model: directions restricted to the (x,z) plane, so vy = 0
//   A2_theta = cos(theta) Z + sin(theta) X
// COMPLEX model: the full sphere, so we can also use Y
function dirXZ(theta) { return [Math.sin(theta), 0, Math.cos(theta)]; }
function dirFull(theta, phi) {
  return [Math.sin(theta) * Math.cos(phi), Math.sin(theta) * Math.sin(phi), Math.cos(theta)];
}

function chshFor(psi, dirsA, dirsB) {
  const [A0, A1] = dirsA, [B0, B1] = dirsB;
  const E = (v, w) => expectation(psi, blochToMatrix(v), blochToMatrix(w));
  return E(A0, B0) + E(A0, B1) + E(A1, B0) - E(A1, B1);
}

console.log('=== does CHSH = 2*sqrt(2) require COMPLEX amplitudes? ===\n');
const psi = bellState(true);

/*
 * For |Phi+> the correlation of two measurement directions is exactly the cosine
 * of the angle between them:
 *        E(thetaA, thetaB) = cos(thetaA - thetaB)
 * for directions in the same plane, whether that plane is (x,z) [real model] or
 * any plane through the sphere [complex model]. We therefore evaluate CHSH
 * analytically instead of on a grid (the previous version nested four loops of
 * 720 points, i.e. ~2.7e11 evaluations, and timed out).
 */
function chshAngles(a0, a1, b0, b1) {
  const E = (u, v) => Math.cos(u - v);
  return E(a0, b0) + E(a0, b1) + E(a1, b0) - E(a1, b1);
}

console.log('--- REAL model: measurement directions confined to one plane ---');
{
  // coarse scan, then local refinement
  let best = -Infinity, bestA = null;
  const N = 180;
  for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) for (let k = 0; k < N; k++) {
    const S = chshAngles(i * Math.PI / N, j * Math.PI / N, k * Math.PI / N, 0);
    if (S > best) { best = S; bestA = [i * Math.PI / N, j * Math.PI / N, k * Math.PI / N, 0]; }
  }
  let cur = bestA.slice(), curS = best;
  for (let iter = 0; iter < 20000; iter++) {
    const idx = iter % 4;
    const step = 0.2 / (1 + iter / 500);
    for (const d of [step, -step]) {
      const trial = cur.slice();
      trial[idx] += d;
      const S = chshAngles(trial[0], trial[1], trial[2], trial[3]);
      if (S > curS) { curS = S; cur = trial; }
    }
  }
  const deg = (a) => (a * 180 / Math.PI);
  const norm = (a) => { let x = deg(a) % 360; if (x < 0) x += 360; return x; };
  console.log(`  max CHSH (real)        : ${curS.toFixed(6)}`);
  console.log(`  optimal angles (deg)   : ${cur.map(a => norm(a).toFixed(1)).join(', ')}`);
  // check the textbook choice explicitly
  const textbook = chshAngles(0, Math.PI / 2, Math.PI / 4, -Math.PI / 4);
  console.log(`  textbook choice A=0,90  B=45,-45 : ${textbook.toFixed(6)}`);
  console.log(`  2*sqrt(2)                        : ${(2 * Math.SQRT2).toFixed(6)}`);
  console.log('');
}

console.log('--- COMPLEX model: full SU(2) local rotations ---');
{
  // Restricting to a plane already saturates the bound, so allowing the full
  // sphere cannot do better. We state the fact and verify the same value.
  console.log(`  the same angle configuration is available in the complex model, so`);
  console.log(`  max CHSH (complex)     : ${chshAngles(0, Math.PI / 2, Math.PI / 4, -Math.PI / 4).toFixed(6)}`);
  console.log(`  (2*sqrt(2) is the universal Tsirelson bound, so neither model exceeds it)`);
  console.log('');
}

console.log('=== what this shows ===');
console.log('  Both models reach 2*sqrt(2) on this state. So the CHSH bound does NOT');
console.log('  distinguish real from complex amplitudes -- the claim I made earlier in');
console.log('  this project ("2*sqrt(2) needs complex amplitudes") is NOT supported by');
console.log('  this test and should be restated.');
console.log('');
console.log('  What IS true, and what the exhaustive classical search showed, is weaker:');
console.log('    - a CLASSICAL BIT caps CHSH at 2');
console.log('    - a continuous REAL amplitude already reaches 2*sqrt(2)');
console.log('  So the missing ingredient is not "complex versus real" but');
console.log('  "continuous amplitude structure versus a discrete slot".');
console.log('');
console.log('=== why this matters for the strategy ===');
console.log('  If a continuous real amplitude is enough for CHSH, then what the particle');
console.log('  needs is a CONTINUUM of local states with a rotation acting on them -- not');
console.log('  specifically the complex numbers. That is a weaker requirement, and a');
console.log('  lattice with a four-fold internal symmetry (Z4) is a candidate source of');
console.log('  exactly such a rotation, with the continuum arising as the coarse-grained');
console.log('  limit of many discrete phases.');
console.log('');
console.log('  CAUTION: CHSH is not the only quantum experiment. There are experiments');
console.log('  that do distinguish real from complex (multi-slit interference with three');
console.log('  or more paths, sequential measurements). I have NOT checked those here and');
console.log('  will not claim a conclusion about them.');
