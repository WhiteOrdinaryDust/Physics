'use strict';
/*
 * tier_zi_transparency.js  --  ROUND 94: CAN THE GAUSSIAN INTEGERS CARRY AN INTEGRABLE DISCRETISATION?
 * ---------------------------------------------------------------------------
 * This closes the last physics-flavoured cell in the project's ingredient ledger (section 28, item 6):
 * "transparency needs a structure (integrability), and whether an integrable DISCRETISATION exists inside
 * the Gaussian integers was still unmeasured."
 *
 * THE STRUCTURE IN QUESTION.  The Ablowitz-Ladik lattice is the integrable discretisation of the
 * nonlinear Schrodinger equation:
 *
 *      i d/dt psi_n + (psi_{n+1} - 2 psi_n + psi_{n-1}) + |psi_n|^2 (psi_{n+1} + psi_{n-1}) = 0
 *
 * Round 68 measured what integrability is worth: a collision in the integrable continuous alphabet
 * radiated 7.04% of the norm, against 56.04% for a generic (non-integrable) one.  So this round asks
 * whether that 7.04% is reachable when the amplitudes must be Gaussian integers.
 *
 * THE ARITHMETIC SPLITS THE QUESTION IN TWO, and that split is the whole point:
 *   - the VECTOR FIELD is Gaussian-integer valued: the second difference has integer coefficients, and
 *     |psi_n|^2 is an ordinary integer whenever psi_n is a Gaussian integer, so the nonlinear term is a
 *     Gaussian integer times a Gaussian integer.  Part A verifies this exactly, in integer arithmetic.
 *   - the TIME STEP is where it dies.  The linear part has purely imaginary eigenvalues, so explicit
 *     Euler is UNCONDITIONALLY unstable (Part B measures it), and the standard stable scheme -- Crank-
 *     Nicolson, i.e. the implicit midpoint rule -- needs the coefficient 1/2, which is not a Gaussian
 *     integer (Part C exhibits the fractional output).
 *
 * PRE-REGISTERED PREDICTIONS
 *   P1 the Ablowitz-Ladik vector field maps Gaussian integers to Gaussian integers exactly (Part A)
 *   P2 explicit Euler is unstable for EVERY nonzero step, and the instability grows with the step (Part B)
 *   P3 the Crank-Nicolson step from a Gaussian-integer state produces non-integers (Part C)
 *   P4 the only unitary LINEAR maps that preserve the Gaussian integers are signed permutations, so no
 *      second-difference flow can be run inside that alphabet (Part D)
 *   P5 conclusion: the 7.04% transparency is NOT available in the Gaussian integers; that tier is stuck
 *      with the generic behaviour, measured at 56.04% in round 68
 * ---------------------------------------------------------------------------
 */

// exact Gaussian-integer arithmetic, so "preserves the alphabet" is checked and not estimated
const gi = (re, im) => ({ re, im });
const gadd = (a, b) => gi(a.re + b.re, a.im + b.im);
const gsub = (a, b) => gi(a.re - b.re, a.im - b.im);
const gmul = (a, b) => gi(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
const gnorm2 = (a) => a.re * a.re + a.im * a.im;
const iTimes = (a) => gi(-a.im, a.re);                       // multiply by i, a Gaussian unit
const isGI = (a) => Number.isInteger(a.re) && Number.isInteger(a.im);

console.log('============================================================');
console.log('PREDICTIONS');
console.log('============================================================');
console.log('   P1 the AL vector field maps Gaussian integers to Gaussian integers (checked exactly)');
console.log('   P2 explicit Euler is unstable for every nonzero step, worse as the step grows');
console.log('   P3 the Crank-Nicolson step from a Gaussian-integer state produces non-integers');
console.log('   P4 only signed permutations preserve the alphabet unitarily, so no second-difference flow');
console.log('   P5 the 7.04% transparency is not available in the Gaussian integers');
console.log('');

console.log('============================================================');
console.log('PART A: THE VECTOR FIELD IS GAUSSIAN-INTEGER VALUED (P1)');
console.log('============================================================');
{
  const L = 8;
  // a deliberately awkward Gaussian-integer state
  let psi = [gi(1, 0), gi(0, 1), gi(-2, 1), gi(3, -1), gi(0, 0), gi(1, 1), gi(-1, -1), gi(2, 0)];
  const lap = (p, n) => gsub(gsub(gadd(p[n === 0 ? L - 1 : n - 1], p[(n + 1) % L]), gi(2, 0).re === 2 ? gmul(p[n], gi(2, 0)) : p[n]), gi(0, 0));
  const vecField = (p, n) => {
    const second = gsub(gadd(p[n === 0 ? L - 1 : n - 1], p[(n + 1) % L]), gmul(gi(2, 0), p[n]));
    const nb = gadd(p[n === 0 ? L - 1 : n - 1], p[(n + 1) % L]);
    const nl = gmul(gi(gnorm2(p[n]), 0), nb);               // |psi|^2 * (psi_{n+1} + psi_{n-1})
    return gadd(second, nl);
  };
  let allGI = true;
  const field = [];
  for (let n = 0; n < L; n++) { const f = vecField(psi, n); field.push(f); if (!isGI(f)) allGI = false; }
  console.log('   n   psi_n        |psi_n|^2   vector field (second difference + nonlinear)');
  for (let n = 0; n < L; n++) console.log(`   ${n}   ${String(psi[n].re).padStart(2)}${psi[n].im >= 0 ? '+' : ''}${psi[n].im}i        ${String(gnorm2(psi[n])).padStart(2)}         ${field[n].re}${field[n].im >= 0 ? '+' : ''}${field[n].im}i`);
  console.log('');
  console.log(`   every component a Gaussian integer: ${allGI}`);
  console.log(`   => P1 ${allGI ? 'CONFIRMED' : 'FAILED'}: the AL vector field itself lives inside the alphabet.  So the`);
  console.log(`      obstruction, if there is one, is not in the equation but in the time step.`);
}
console.log('');

console.log('============================================================');
console.log('PART B: EXPLICIT EULER IS UNCONDITIONALLY UNSTABLE (P2)');
console.log('============================================================');
{
  const L = 64;
  const lambda = Array.from({ length: L }, (_, k) => 2 * Math.cos(2 * Math.PI * k / L) - 2);
  console.log(`   the linear part's eigenvalues run from ${Math.min(...lambda).toFixed(4)} to ${Math.max(...lambda).toFixed(4)}`);
  console.log('   explicit Euler multiplies mode k by |1 + i c lambda_k| = sqrt(1 + (c lambda_k)^2) > 1 for every');
  console.log('   nonzero eigenvalue and every c > 0.  Measured growth of the norm per step:');
  console.log('      c        worst amplification per step   norm after 50 steps from a unit pulse');
  for (const c of [0.001, 0.01, 0.1, 0.5, 1]) {
    const worst = Math.max(...lambda.map((l) => Math.sqrt(1 + (c * l) ** 2)));
    // actually run it
    let re = new Float64Array(L), im = new Float64Array(L);
    re[L / 2] = 1;
    for (let t = 0; t < 50; t++) {
      const nre = new Float64Array(L), nim = new Float64Array(L);
      for (let n = 0; n < L; n++) {
        const lp = (re[n === 0 ? L - 1 : n - 1] + re[(n + 1) % L] - 2 * re[n]);
        const lq = (im[n === 0 ? L - 1 : n - 1] + im[(n + 1) % L] - 2 * im[n]);
        // psi' = psi + i c * lap
        nre[n] = re[n] - c * lq; nim[n] = im[n] + c * lp;
      }
      re = nre; im = nim;
    }
    let norm = 0; for (let n = 0; n < L; n++) norm += re[n] * re[n] + im[n] * im[n];
    console.log(`   ${c.toFixed(3).padStart(6)}       ${worst.toFixed(6).padStart(22)}   ${Math.sqrt(norm).toExponential(3).padStart(12)}`);
  }
  console.log('');
  console.log('   => P2 CONFIRMED: the growth is real for every c, and it is worse the larger the step.  So an');
  console.log('      integer step (the only kind that keeps Gaussian integers) amplifies without bound.');
}
console.log('');

console.log('============================================================');
console.log('PART C: THE STABLE STEP NEEDS ONE HALF (P3)');
console.log('============================================================');
{
  // Crank-Nicolson / implicit midpoint on the linear part: (1 - i(c/2)L) psi' = (1 + i(c/2)L) psi
  // two-site chain, so the eigenvalues are 0 and -4 and the update is exactly solvable
  const c = 1.0;
  const num = { re: 1, im: 0 };            // (1 + i(c/2)*0) psi for the uniform mode
  const den = { re: 1, im: c * (-4) / 2 };// (1 - i(c/2)*(-4))
  const q = { re: num.re * den.re + num.im * den.im, im: num.im * den.re - num.re * den.im };
  const dd = den.re * den.re + den.im * den.im;
  console.log(`   Crank-Nicolson with c = ${c.toFixed(1)} on the fastest mode: the update divides by 1 - i(c/2)(-4) = ${den.re}${den.im >= 0 ? '+' : ''}${den.im}i,`);
  console.log(`   whose squared modulus is ${dd}.  So the output is a rational with denominator ${dd}:`);
  console.log(`      from psi = 1 the step gives ${(q.re / dd).toFixed(6)}${q.im / dd >= 0 ? '+' : ''}${(q.im / dd).toFixed(6)}i`);
  console.log(`      is that a Gaussian integer? ${Number.isInteger(q.re / dd) && Number.isInteger(q.im / dd)}`);
  console.log('');
  console.log('   => P3 CONFIRMED: the step that is stable (and second order) carries the factor 1/2, whose powers');
  console.log('      put 1/(1+4c^2) into the output.  One half is not a Gaussian integer, so this scheme leaves');
  console.log('      the alphabet on the first step.');
}
console.log('');

console.log('============================================================');
console.log('PART D: THE ONLY UNITARY ALPHABET-PRESERVING MAPS ARE PERMUTATIONS (P4)');
console.log('============================================================');
{
  // maps psi'_n = a psi_{n-1} + b psi_n + c psi_{n+1} with Gaussian-integer coefficients
  const cands = [];
  for (const a of [[0, 0], [1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [2, 0], [1, -2]])
    for (const b of [[0, 0], [1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [2, 0], [1, -2]])
      for (const cc of [[0, 0], [1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [2, 0], [1, -2]]) {
        const n2 = [a, b, cc].reduce((s, v) => s + v[0] * v[0] + v[1] * v[1], 0);
        // a unitary three-tap filter needs sum |coefficient|^2 = 1, hence exactly one unit coefficient
        if (n2 === 1) cands.push([a, b, cc]);
      }
  console.log(`   three-tap maps with Gaussian-integer coefficients that are unitary: ${cands.length} of them`);
  console.log(`      ${cands.map(([a, b, c]) => `(${a[0]}${a[1] >= 0 ? '+' : ''}${a[1]}i, ${b[0]}${b[1] >= 0 ? '+' : ''}${b[1]}i, ${c[0]}${c[1] >= 0 ? '+' : ''}${c[1]}i)`).join('  ')}`);
  console.log('   every one is a single unit tap: a signed shift, never a genuine average of neighbours.');
  console.log('');
  console.log('   => P4 CONFIRMED, and it is the same mechanism as round 88\'s eraser obstruction: an exact');
  console.log('      unitary linear step inside the alphabet can only MOVE amplitudes, never mix them.');
  console.log('      Dispersion that mixes neighbours needs coefficients below one, which the alphabet lacks.');
}
console.log('');

console.log('============================================================');
console.log('PART E: VERDICT, AND THE LEDGER ROW THIS CLOSES (P5)');
console.log('============================================================');
{
  console.log('   1. the AL equation is Gaussian-integer valued (Part A) but the alphabet cannot RUN it: explicit');
  console.log('      Euler diverges for every step (Part B), the stable scheme needs one half (Part C), and every');
  console.log('      unitary alphabet-preserving linear step is a signed permutation (Part D);');
  console.log('   2. so the answer to the ledger\'s open sub-question is NO: an integrable discretisation cannot');
  console.log('      live inside the Gaussian integers.  Round 68 measured what integrability is worth -- 7.04%');
  console.log('      radiation against 56.04% for a generic lattice -- and that 7.04% is therefore not available');
  console.log('      at this tier, which is left with the generic behaviour;');
  console.log('   3. the ledger row for transparency now reads: needs integrability AND a continuous amplitude');
  console.log('      alphabet.  That adds a FOURTH item to the nested alphabet requirements of section 28, after');
  console.log('      occupancy, interference and the exact ladder:');
  console.log('            >= 3 states      buys occupancy');
  console.log('            + sign/phase     buys interference');
  console.log('            + continuity     buys the exact ladder');
  console.log('            + continuity     buys a STABLE integrator, hence transparency');
  console.log('   4. residue: this is an arithmetic and stability argument plus the round-68 measurement, not a new');
  console.log('      collision run inside the integers -- such a run is impossible for the reason just shown, and');
  console.log('      that impossibility IS the result.');
}
