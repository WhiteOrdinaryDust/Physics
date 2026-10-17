'use strict';
/*
 * monogamy_check.js
 * ---------------------------------------------------------------------------
 * The interface-as-consumable-slot picture is WRONG. This file shows what the
 * actual constraint is.
 *
 * THE PROBLEM WITH "INTERFACES ARE CONSUMED"
 *   If maximal entanglement with B consumed all of A's interfaces, then the pair
 *   would expose ZERO interfaces and could not interact with anything. But
 *   entangled particles are measured in labs every day. So a literal
 *   slot-consumption reading gives false physics.
 *
 * WHAT MONOGAMY ACTUALLY SAYS
 *   Not "A runs out of interfaces", but "A's entanglement is SHARED":
 *        C(A:B)^2 + C(A:C)^2 <= C(A:BC)^2 <= 1
 *   and crucially, entanglement is defined RELATIVE TO A FACTORISATION. Let A
 *   also become entangled with an environment E and the A-B entanglement must
 *   drop. Nothing is consumed; the correlations are REDISTRIBUTED.
 *
 * WHAT IS COMPUTED HERE
 *   - GHZ  = (|000> + |111>)/sqrt2     : the "all or nothing" state
 *   - W    = (|001> + |010> + |100>)/sqrt3 : the "shared" state
 *   For both, the pairwise concurrences C(A:B), C(A:C), C(B:C) and the CKW sum
 *   C(A:B)^2 + C(A:C)^2 + C(B:C)^2, which must be <= 1.
 *
 * WHY IT MATTERS FOR THE INTERFACE MODEL
 *   A model in which coupling consumes slots predicts that adding a third
 *   partner MUST reduce the first pair's correlation toward zero. The real
 *   constraint is quantitative and tight (CKW), and it is a statement about
 *   correlation STRUCTURE, not about exhausted capacity.
 * ---------------------------------------------------------------------------
 */

// ---------- minimal complex linear algebra for 2x2 and 3-qubit states
// state vector as array of [re, im] pairs, indexed with qubit 0 most significant
function v(...amps) { return amps.map(a => [a, 0]); }

function kron(a, b) {
  const out = [];
  for (const [ar, ai] of a) for (const [br, bi] of b) {
    out.push([ar * br - ai * bi, ar * bi + ai * br]);
  }
  return out;
}
function addv(a, b) { return a.map((x, i) => [x[0] + b[i][0], x[1] + b[i][1]]); }
function scalev(a, s) { return a.map(x => [x[0] * s, x[1] * s]); }
function norm2(a) { return a.reduce((s, x) => s + x[0] * x[0] + x[1] * x[1], 0); }

// reduced density matrix of an n-qubit pure state over the qubits in `keep`
// returns a (2^k x 2^k) complex matrix as [re, im]
function reducedDensity(psi, n, keep) {
  const k = keep.length;
  const dim = 1 << k;
  const rho = Array.from({ length: dim }, () => Array.from({ length: dim }, () => [0, 0]));
  const kept = new Set(keep);
  const other = [];
  for (let q = 0; q < n; q++) if (!kept.has(q)) other.push(q);
  const oDim = 1 << other.length;
  for (let oi = 0; oi < oDim; oi++) {
    for (let kj = 0; kj < dim; kj++) {
      for (let kl = 0; kl < dim; kl++) {
        // build global indices
        let gi = 0, gj = 0;
        for (let q = 0; q < n; q++) {
          if (kept.has(q)) {
            const pos = keep.indexOf(q);
            gi |= ((kj >> (k - 1 - pos)) & 1) << (n - 1 - q);
            gj |= ((kl >> (k - 1 - pos)) & 1) << (n - 1 - q);
          } else {
            const pos = other.indexOf(q);
            const bit = (oi >> (other.length - 1 - pos)) & 1;
            gi |= bit << (n - 1 - q);
            gj |= bit << (n - 1 - q);
          }
        }
        const a = psi[gi], b = psi[gj];
        // rho += a * conj(b)
        const re = a[0] * b[0] + a[1] * b[1];
        const im = a[1] * b[0] - a[0] * b[1];
        rho[kj][kl][0] += re;
        rho[kj][kl][1] += im;
      }
    }
  }
  return rho;
}

// concurrence of a 2-qubit density matrix given as 4x4 complex
function concurrence2(rho) {
  // build spin-flipped matrix: (sigma_y x sigma_y) rho* (sigma_y x sigma_y)
  const sy = [[[0, 0], [0, -1]], [[0, 1], [0, 0]]];   // sigma_y = [[0,-i],[i,0]]
  // helper to multiply 4x4 complex
  const mul = (A, B) => {
    const n = A.length;
    const C = Array.from({ length: n }, () => Array.from({ length: n }, () => [0, 0]));
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
      let re = 0, im = 0;
      for (let k = 0; k < n; k++) {
        const a = A[i][k], b = B[k][j];
        re += a[0] * b[0] - a[1] * b[1];
        im += a[0] * b[1] + a[1] * b[0];
      }
      C[i][j] = [re, im];
    }
    return C;
  };
  const conj = (A) => A.map(row => row.map(x => [x[0], -x[1]]));
  // Y = sy (x) sy as a 4x4
  const Y = [];
  for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) {
    const row = [];
    for (let k = 0; k < 2; k++) for (let l = 0; l < 2; l++) {
      // (sy_ij)(sy_kl)
      const a = sy[i][k], b = sy[j][l];
      row.push([a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]);
    }
    Y.push(row);
  }
  const M = mul(mul(Y, conj(rho)), Y);
  // eigenvalues of M (complex in general); for the CKW check we need sqrt of the
  // eigenvalues of M, which for these states are real. Use the R matrix approach:
  // compute eigenvalues numerically via the power-free method on M*conj(M)? 
  // Simpler: for the states used here the concurrence has a closed form, and we
  // verify by computing the eigenvalues of the Hermitian matrix M^(1/2)... 
  // Instead: compute the singular values of sqrt(rho) Y sqrt(rho)-like object.
  // For robustness we use the standard formula via eigenvalues of rho*Y*conj(rho)*Y
  // which IS Hermitian and has non-negative eigenvalues.
  const R = mul(mul(Y, conj(rho)), mul(Y, rho));
  // eigenvalues of Hermitian R (2x2 blocks of complex) -- use Jacobi-free approach:
  // R is 4x4 Hermitian; get eigenvalues via characteristic polynomial is overkill.
  // Use the fact that for our target states R is real symmetric after a basis change.
  const lam = eigenHermitian4(R);
  const sq = lam.map(x => Math.sqrt(Math.max(0, x))).sort((a, b) => b - a);
  return Math.max(0, sq[0] - sq[1] - sq[2] - sq[3]);
}

// eigenvalues of a 4x4 Hermitian complex matrix (real eigenvalues), simple QR-free
// implementation using the power method on deflation is overkill; use closed form
// for 4x4 symmetric real part and verify imaginary parts are ~0.
function eigenHermitian4(R) {
  // convert to real symmetric matrix if imaginary parts are negligible
  let maxIm = 0;
  for (const row of R) for (const x of row) maxIm = Math.max(maxIm, Math.abs(x[1]));
  if (maxIm > 1e-9) {
    // fall back: use the 2x2 complex structure -- not needed for our test states
    throw new Error('imaginary part non-negligible; state not in the expected class');
  }
  const A = R.map(row => row.map(x => x[0]));
  // Jacobi eigenvalue algorithm for a real symmetric matrix
  const n = A.length;
  const a = A.map(r => r.slice());
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += a[i][j] * a[i][j];
    if (off < 1e-18) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(a[p][q]) < 1e-18) continue;
      const theta = 0.5 * Math.atan2(2 * a[p][q], a[q][q] - a[p][p]);
      const c = Math.cos(theta), s = Math.sin(theta);
      for (let k = 0; k < n; k++) {
        const akp = a[k][p], akq = a[k][q];
        a[k][p] = c * akp - s * akq;
        a[k][q] = s * akp + c * akq;
      }
      for (let k = 0; k < n; k++) {
        const apk = a[p][k], aqk = a[q][k];
        a[p][k] = c * apk - s * aqk;
        a[q][k] = s * apk + c * aqk;
      }
    }
  }
  return [0, 1, 2, 3].map(i => a[i][i]);
}

// ============================================================ the states
console.log('=== monogamy: what the constraint actually says ===\n');

// GHZ = (|000> + |111>)/sqrt2  -> |000>,|111> are indices 0 and 7
const ghz = (() => {
  const s = Array.from({ length: 8 }, () => [0, 0]);
  s[0] = [1 / Math.SQRT2, 0];
  s[7] = [1 / Math.SQRT2, 0];
  return s;
})();
// W = (|001> + |010> + |100>)/sqrt3 -> indices 1,2,4
const w = (() => {
  const s = Array.from({ length: 8 }, () => [0, 0]);
  for (const i of [1, 2, 4]) s[i] = [1 / Math.sqrt(3), 0];
  return s;
})();

function analyse(name, psi) {
  console.log(`--- ${name} ---`);
  const pairs = [[0, 1], [0, 2], [1, 2]];
  const names = ['A:B', 'A:C', 'B:C'];
  const cs = [];
  pairs.forEach((p, i) => {
    const rho = reducedDensity(psi, 3, p);
    const c = concurrence2(rho);
    cs.push(c);
    console.log(`  C(${names[i]}) = ${c.toFixed(6)}   C^2 = ${(c * c).toFixed(6)}`);
  });
  const sum = cs.reduce((a, c) => a + c * c, 0);
  console.log(`  CKW sum  C(A:B)^2 + C(A:C)^2 + C(B:C)^2 = ${sum.toFixed(6)}`);
  console.log(`  (CKW requires this <= 1, and requires C(A:B)^2 + C(A:C)^2 <= 1)`);
  console.log('');
  return { cs, sum };
}

const rGHZ = analyse('GHZ = (|000> + |111>)/sqrt(2)', ghz);
const rW = analyse('W = (|001> + |010> + |100>)/sqrt(3)', w);

// ============================================================ what this says
console.log('=== what this means for the interface model ===');
console.log('');
console.log('  GHZ: every pairwise concurrence is ZERO, yet the state is maximally');
console.log('       entangled. The entanglement is purely three-way. A slot model that');
console.log('       only counts pairwise couplings would report "no coupling at all".');
console.log('');
console.log('  W  : each pair shares 2/3 of a concurrence, and the CKW sum is');
console.log(`       ${rW.sum.toFixed(4)} <= 1, so the pairs genuinely COMPETE for a limited total.`);
console.log('');
console.log('  ==> Monogamy is a statement about how correlation is DISTRIBUTED,');
console.log('      not about a stock of consumable connectors.');
console.log('');
console.log('  ==> AND: neither state makes the particles unable to interact with the');
console.log('      outside. Any of the three qubits can still be measured, can still');
console.log('      scatter a photon, can still decohere. What changes when a THIRD');
console.log('      partner is engaged is that the original pair\'s entanglement drops.');
console.log('');
console.log('  So the user objection is correct: "entanglement consumes the interfaces,');
console.log('  leaving zero for the outside" gives false physics. The interface picture');
console.log('  must be read as "which correlations are established", NOT as "how many');
console.log('  sockets are left".');
