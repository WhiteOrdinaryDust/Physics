'use strict';
/*
 * shortcut_chsh.js
 * ---------------------------------------------------------------------------
 * JUDGEMENT RUN for the architecture:
 *
 *   L0  substrate      : minimal deterministic local rules
 *   L1  particle layer : particles emerge from L0; NOT bits --
 *                        they carry complex attributes (phase/amplitude-like state)
 *   L2  world layer    : spatial relations emerge from particle couplings
 *   photon             : the minimal event of L1; travels geodesics
 *
 * The question:
 *   Can a CLASSICAL DETERMINISTIC SHORTCUT EDGE between two particles produce
 *   CHSH = 2*sqrt(2) while keeping no-signaling?
 *
 * Four shortcut types are computed and each is checked for no-signaling.
 * ---------------------------------------------------------------------------
 */

// ---------------------------------------------------------------- CHSH helpers
// correlations E(x,y) = P(a==b) - P(a!=b), x,y in {0,1}
// S = E00 + E01 + E10 - E11
const AOPTS = [[0, 0], [0, 1], [1, 0], [1, 1]]; // [a(x=0), a(x=1)]

function emptyJoint() {
  return [
    [[[0, 0], [0, 0]], [[0, 0], [0, 0]]],
    [[[0, 0], [0, 0]], [[0, 0], [0, 0]]],
  ];
}

function chshFromJoint(joint) {
  const E = [[0, 0], [0, 0]];
  for (let x = 0; x < 2; x++) {
    for (let y = 0; y < 2; y++) {
      let same = 0, diff = 0;
      for (let a = 0; a < 2; a++) {
        for (let b = 0; b < 2; b++) {
          const p = joint[x][y][a][b];
          if (a === b) same += p; else diff += p;
        }
      }
      E[x][y] = same - diff;
    }
  }
  return { E, S: E[0][0] + E[0][1] + E[1][0] - E[1][1] };
}

// max marginal shift: 0 means no-signalling, >0 means a usable superluminal channel
function signalingCheck(joint) {
  let maxDev = 0;
  // P(a=0|x) must not depend on y
  for (let x = 0; x < 2; x++) {
    let p0y0 = 0, p0y1 = 0;
    for (let b = 0; b < 2; b++) { p0y0 += joint[x][0][0][b]; p0y1 += joint[x][1][0][b]; }
    maxDev = Math.max(maxDev, Math.abs(p0y0 - p0y1));
  }
  // P(b=0|y) must not depend on x
  for (let y = 0; y < 2; y++) {
    let q0x0 = 0, q0x1 = 0;
    for (let a = 0; a < 2; a++) { q0x0 += joint[0][y][a][0]; q0x1 += joint[1][y][a][0]; }
    maxDev = Math.max(maxDev, Math.abs(q0x0 - q0x1));
  }
  return maxDev;
}

// ---------------------------------------------------------------- CASE 1
// Shortcut = shared randomness ("lambda"), the most general thing a classical
// deterministic edge can supply. lambda uniform on NL values:
//   P(a,b|x,y) = (1/NL) * sum_lambda [a=A(x,lambda)] [b=B(y,lambda)]
// For fixed lambda this is a PRODUCT distribution. A no-signalling mixture of
// product distributions IS a local-hidden-variable model, so CHSH <= 2 is forced.
// We enumerate exhaustively on a small lambda space to confirm it.
function case1_exhaustive(NL) {
  const nA = Math.pow(4, NL), nB = Math.pow(4, NL);
  const total = nA * nB;
  if (total > (1 << 24)) return { tooLarge: true, total };

  let bestAny = -Infinity, bestNS = -Infinity;
  let bestAnyMarginal = 0, bestNSMarginal = 0;
  const Acfg = new Uint8Array(NL), Bcfg = new Uint8Array(NL);

  for (let ia = 0; ia < nA; ia++) {
    let t = ia;
    for (let l = 0; l < NL; l++) { Acfg[l] = t & 3; t >>= 2; }
    for (let ib = 0; ib < nB; ib++) {
      let u = ib;
      for (let l = 0; l < NL; l++) { Bcfg[l] = u & 3; u >>= 2; }

      const joint = emptyJoint();
      for (let l = 0; l < NL; l++) {
        const A = AOPTS[Acfg[l]], B = AOPTS[Bcfg[l]];
        for (let x = 0; x < 2; x++) {
          for (let y = 0; y < 2; y++) joint[x][y][A[x]][B[y]] += 1 / NL;
        }
      }
      const { S } = chshFromJoint(joint);
      const dev = signalingCheck(joint);
      if (S > bestAny) { bestAny = S; bestAnyMarginal = dev; }
      if (dev < 1e-12 && S > bestNS) { bestNS = S; bestNSMarginal = dev; }
    }
  }
  return { bestAny, bestAnyMarginal, bestNS, bestNSMarginal, total };
}

// ---------------------------------------------------------------- CASE 2
// Shortcut that CAN signal: Bob's outcome is allowed to depend on Alice's setting.
// CHSH = 4 becomes reachable, but the marginals move with the other setting.
//
// Note on the old bug: the pattern b = (x XOR y) XOR 1 gives E = [-1,+1,+1,-1],
// i.e. S = -1+1+1-(-1) = +2, not 4. The maximal pattern needs E = [+1,+1,+1,-1],
// which is b = x AND y:
//        b=1
// x=0  [ 1  -1 ]
// x=1  [-1   1 ]
function case2_signalling() {
  const joint = emptyJoint();
  for (let x = 0; x < 2; x++) {
    for (let y = 0; y < 2; y++) {
      const a = 0;
      const b = (x === 1 && y === 1) ? 1 : 0;
      joint[x][y][a][b] += 1;
    }
  }
  return { joint, ...chshFromJoint(joint), dev: signalingCheck(joint) };
}

// ---------------------------------------------------------------- CASE 3
// Shortcut = a resource carrying COMPLEX attributes (quantum).
// |Phi+> = (|00> + |11>)/sqrt(2) with observables at angles theta:
//   E(x,y) = cos(thetaA_x - thetaB_y)
// Optimal angles give 2*sqrt(2) (Tsirelson).
function case3_exact(angles) {
  const E = [[0, 0], [0, 0]];
  for (let x = 0; x < 2; x++) {
    for (let y = 0; y < 2; y++) E[x][y] = Math.cos(angles[x] - angles[2 + y]);
  }
  return { E, S: E[0][0] + E[0][1] + E[1][0] - E[1][1] };
}

// Same thing sampled, to show it is a real operational model and not just a formula.
//
// Three failed attempts lived here, and the failure mode was always the same:
// substituting a hand-derived closed form for the complex inner product. The
// phases carry the sign structure, so a real-valued shortcut silently produces a
// state that is not even a valid measurement model (angle-dependent marginals).
// The fix is to build the basis amplitudes and take the inner product.
//
// |Phi+> = (|00> + |11>)/sqrt2
// Alice(+-) basis: |+_A> = (|0> + e^{+i tA}|1>)/sqrt2
// Bob(+-)   basis: |+_B> = (|0> + e^{-i tB}|1>)/sqrt2
function basisPlus(t) { return [1 / Math.SQRT2, Math.cos(t) / Math.SQRT2, Math.sin(t) / Math.SQRT2]; }
function basisMinus(t) { return [1 / Math.SQRT2, -Math.cos(t) / Math.SQRT2, -Math.sin(t) / Math.SQRT2]; }

// <Phi+| vA (x) vB> with v = [re0, re1, im1] (component 0 is purely real)
function overlapPhiPlus(vA, vB) {
  const re = vA[0] * vB[0] + (vA[1] * vB[1] - vA[2] * vB[2]);
  const im = vA[1] * vB[2] + vA[2] * vB[1];
  return [re / Math.SQRT2, im / Math.SQRT2];
}
function bornProbs(tA, tB) {
  const pA = basisPlus(tA), mA = basisMinus(tA);
  const pB = basisPlus(-tB), mB = basisMinus(-tB);
  const o = (vA, vB) => { const c = overlapPhiPlus(vA, vB); return c[0] * c[0] + c[1] * c[1]; };
  return { pp: o(pA, pB), pm: o(pA, mB), mp: o(mA, pB), mm: o(mA, mB) };
}

function case3_sampled(angles, trialsPerPair, settings, seed) {
  let s = seed >>> 0;
  const rnd = () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
  const joint = emptyJoint();
  for (const [x, y] of settings) {
    const b = bornProbs(angles[x], angles[2 + y]);
    const cum = [b.pp, b.pp + b.pm, b.pp + b.pm + b.mp, 1];
    for (let n = 0; n < trialsPerPair; n++) {
      const u = rnd();
      let idx = 0;
      while (idx < 3 && u > cum[idx]) idx++;
      const a = (idx === 0 || idx === 1) ? 0 : 1;
      const bb = (idx === 0 || idx === 2) ? 0 : 1;
      joint[x][y][a][bb] += 1;
    }
    for (let a = 0; a < 2; a++)
      for (let bb = 0; bb < 2; bb++) joint[x][y][a][bb] /= trialsPerPair;
  }
  return { joint, ...chshFromJoint(joint), dev: signalingCheck(joint) };
}

// ---------------------------------------------------------------- CASE 4
// PR box: no-signalling, CHSH = 4, NOT quantum. Shows that no-signaling alone
// does not cap CHSH at 2*sqrt(2) -- a further principle is needed.
function case4_prBox() {
  const joint = emptyJoint();
  const table = [
    [[[0, 0], [1, 1]], [[0, 0], [1, 1]]],
    [[[0, 0], [1, 1]], [[0, 1], [1, 0]]],
  ];
  for (let x = 0; x < 2; x++)
    for (let y = 0; y < 2; y++)
      for (const [a, b] of table[x][y]) joint[x][y][a][b] += 0.5;
  return { joint, ...chshFromJoint(joint), dev: signalingCheck(joint) };
}

// ---------------------------------------------------------------- run
console.log('=== shortcut / CHSH judgement run ===\n');

const NL = 6;
const c1 = case1_exhaustive(NL);
console.log('CASE 1: shortcut = shared randomness (classical deterministic edge)');
if (c1.tooLarge) {
  console.log(`  (lambda space too large: ${c1.total})`);
} else {
  console.log(`  lambda space ${NL}, exhaustive over ${c1.total.toLocaleString()} models`);
  console.log(`  max CHSH, any model            : ${c1.bestAny.toFixed(6)}   (marginal shift ${c1.bestAnyMarginal.toFixed(6)})`);
  console.log(`  max CHSH, no-signalling models : ${c1.bestNS.toFixed(6)}   (marginal shift ${c1.bestNSMarginal.toFixed(6)})`);
}
console.log('  -> the ceiling is 2, and the models that reach 2 mostly do not signal\n');

const c2 = case2_signalling();
console.log('CASE 2: shortcut allowed to SIGNAL (b depends on x)');
console.log(`  CHSH                            : ${c2.S.toFixed(6)}`);
console.log(`  max marginal shift              : ${c2.dev.toFixed(6)}`);
console.log('  -> 4 is reachable, but only by opening a superluminal channel\n');

const angles = [0, Math.PI / 2, Math.PI / 4, -Math.PI / 4];
const c3 = case3_exact(angles);
console.log('CASE 3: shortcut = resource with COMPLEX attributes (quantum)');
console.log(`  observables A = [0, pi/2], B = [pi/4, -pi/4]`);
for (let x = 0; x < 2; x++) {
  console.log(`    E[x=${x}] = [${c3.E[x].map(v => v.toFixed(6)).join(', ')}]`);
}
console.log(`  CHSH                            : ${c3.S.toFixed(6)}`);
console.log(`  2*sqrt(2)                       : ${(2 * Math.SQRT2).toFixed(6)}`);
const ALL = [[0, 0], [0, 1], [1, 0], [1, 1]];
const c3s = case3_sampled(angles, 200000, ALL, 7);
console.log(`  sampled check (800k trials)     : ${c3s.S.toFixed(4)}`);
const c3ns = case3_sampled(angles, 1000000, ALL, 11);
console.log(`  sampled no-signalling shift     : ${c3ns.dev.toExponential(2)}`);
console.log('  -> no-signalling AND above 2: requires complex amplitudes\n');

const c4 = case4_prBox();
console.log('CASE 4: PR box (non-signalling, stronger than quantum)');
console.log(`  CHSH                            : ${c4.S.toFixed(6)}`);
console.log(`  max marginal shift              : ${c4.dev.toFixed(6)}`);
console.log('  -> no-signalling alone does NOT cap CHSH at 2*sqrt(2)\n');

console.log('=== summary: what each shortcut type buys ===');
console.log('shortcut type                        CHSH     no-signalling?');
console.log('shared randomness (classical)        2.000    yes');
console.log('signalling channel                   4.000    NO');
console.log('complex attributes (quantum)         2.828    yes');
console.log('PR box                               4.000    yes (non-quantum)');

// ---------------------------------------------------------------- scan of the quantum family
// How much of [2, 2*sqrt(2)] can a complex-attribute resource cover?
console.log('\n=== CHSH vs measurement-angle spread (complex-attribute resource) ===');
console.log('  phi (B angles = +-phi)   CHSH');
for (let deg = 0; deg <= 45; deg += 5) {
  const phi = deg * Math.PI / 180;
  const S = case3_exact([0, Math.PI / 2, phi, -phi]).S;
  console.log(`  ${String(deg).padStart(3)} deg                  ${S.toFixed(6)}`);
}
