'use strict';
/*
 * e4_generic.js   --   ROUND 74: IS THE CONSPIRACY *DERIVED*, OR ONLY *ENGINEERED*?
 * ---------------------------------------------------------------------------
 * Round 73 wrote the conspiracy down exactly:  q(theta) = |sin 2*theta|/4  in closed form, zero free
 * parameters, no-signalling, reproducing -cos D for every angle.  But q was REVERSE-ENGINEERED from
 * the target.  The task named there: get the same q from a rule that never mentions -cos(D).
 *
 * THIS ROUND SIMPLIFIES THE QUESTION FIRST.  Fix the two rules to purely geometric ones:
 *
 *      (R1) BISECTOR RULE  : the source's absolute phase is  theta = (a+b)/2 + psi ,
 *                            where psi is the source's OWN internal phase (its hidden state).
 *                            "The source sits between the two devices and references its phase to
 *                            their bisector" -- no target involved, and no D involved.
 *      (R2) PROJECTION RULE: A = r*sign(cos(theta-a)),  B = -r*sign(cos(theta-b))
 *                            (a bit is the SIGN of a projection: Malus-type, geometric).
 *
 * With those two fixed, the whole correlation is a functional of ONE object -- the source's own
 * measure p(psi):
 *
 *      E(D) = - integral p(psi) * sign(cos(psi - D/2)) * sign(cos(psi + D/2)) dpsi
 *
 * So "deriving E4 from a substrate" reduces to a single, decidable question:
 *
 *      DOES A LOCAL RULE GENERATE THE MEASURE p(psi) = |sin 2*psi| / 4 ?
 *
 * and this script measures how special that measure is: against natural closed forms, against a
 * one-parameter deformation, and against the invariant measures of simple deterministic maps.
 *
 * PARTS
 *   A  the machinery: the engineered measure vs the natural ones (uniform, isotropic sphere angle,
 *      von Mises, cos^2, a cos(2psi) family).  Deviation from -cos, and S.  (The control: uniform
 *      must give the triangle wave and S = 2.)
 *   B  how special is it: admixtures with the uniform measure, and the |sin 2psi|^p family.
 *   C  A DERIVATION ATTEMPT: the invariant measures of simple deterministic circle maps
 *      psi -> psi + w + K*sin(m*psi) (a standard, local, parameter-free-in-spirit family).
 *      Best achievable deviation and S.
 *   D  THE HONEST CATCH, MEASURED: a map with the required invariant measure EXISTS in closed form
 *      (conjugate the doubling map by the required CDF F(psi) = (1-cos 2psi)/4).  But conjugating by
 *      a DIFFERENT F gives a DIFFERENT measure, i.e. the measure is whatever you put into F.
 *   E  no-target audit, and verdict.
 * ---------------------------------------------------------------------------
 */

const N = 4096;
const PS = Array.from({ length: N }, (_, i) => (2 * Math.PI * i) / N);
const sgn = (x) => (Math.cos(x) >= 0 ? 1 : -1);
const norm = (w) => { const s = w.reduce((a, b) => a + b, 0); return w.map((v) => v / s); };
const TGT = (d) => -Math.cos((d * Math.PI) / 180);

// the ONE functional.  NOTE: it takes the measure and the angle; the target never enters.
function E_of(p, dDeg) {
  const h = (dDeg * Math.PI) / 180 / 2;
  let e = 0;
  for (let i = 0; i < N; i++) e += p[i] * sgn(PS[i] - h) * sgn(PS[i] + h);
  return -e;
}
function report(p, label) {
  let worst = 0, at = null;
  for (let d = 0; d <= 180; d += 1) { const v = Math.abs(E_of(p, d) - TGT(d)); if (v > worst) { worst = v; at = d; } }
  // CHSH at the four optimal settings (0,45),(0,-45),(90,45),(90,-45); E is even in D, so this is |3E(45) - E(135)|
  const S4 = Math.abs(E_of(p, 45) + E_of(p, 45) + E_of(p, 45) - E_of(p, 135));
  return { label, worst, at, S: S4 };
}

const ENGINEERED = norm(PS.map((t) => Math.abs(Math.sin(2 * t)) / 4));

console.log('============================================================');
console.log('A: THE TWO RULES ARE FIXED AND GEOMETRIC; ONLY THE SOURCE MEASURE VARIES');
console.log('============================================================');
console.log('   (R1) bisector rule   : theta = (a+b)/2 + psi        <- the source references its phase to');
console.log('                                                          the bisector of the two devices');
console.log('   (R2) projection rule : A = r*sign(cos(theta-a)),  B = -r*sign(cos(theta-b))');
console.log('   => E(D) = - integral p(psi) sign(cos(psi-D/2)) sign(cos(psi+D/2)) dpsi');
console.log('   => the SHAPE of the correlation is a functional of the source measure p(psi) ALONE.');
console.log('   The target -cos(D) enters only in the comparison columns below, never in the rules.');
console.log('');
{
  const cases = [
    ['engineered  |sin 2psi|/4', ENGINEERED],
    ['uniform', norm(PS.map(() => 1))],
    ['isotropic sphere angle  |sin psi|/4', norm(PS.map((t) => Math.abs(Math.sin(t)) / 4))],
    ['cos^2 psi', norm(PS.map((t) => Math.cos(t) ** 2))],
    ['sin^2 psi', norm(PS.map((t) => Math.sin(t) ** 2))],
    ['(1+cos 2psi)/2pi', norm(PS.map((t) => 1 + Math.cos(2 * t)))],
    ['(1-cos 2psi)/2pi', norm(PS.map((t) => 1 - Math.cos(2 * t)))],
  ];
  for (const k of [0.5, 1, 2, 4, 8]) cases.push([`von Mises kappa=${k}`, norm(PS.map((t) => Math.exp(k * Math.cos(t))))]);
  console.log('   measure                                 max dev from -cos     at D     CHSH S');
  for (const [label, p] of cases) {
    const r = report(p, label);
    console.log(`   ${label.padEnd(38)}  ${r.worst.toExponential(3).padStart(12)}   ${String(r.at).padStart(5)}   ${r.S.toFixed(6)}`);
  }
  console.log('');
  console.log('   READ THIS TABLE CAREFULLY -- it says more than "only one works":');
  console.log('     * the engineered measure lands on 2.828427 EXACTLY (its 1.4e-3 deviation is the 4096-bin grid);');
  console.log('     * natural measures scatter S BOTH BELOW AND ABOVE 2*sqrt(2): from 1.0997 (sin^2) up to');
  console.log('       2.9003 (cos^2).  So measurement dependence alone does not pin the value at all;');
  console.log('     * and anything ABOVE 2.828427 is excluded by experiment (quantum mechanics never exceeds');
  console.log('       the Tsirelson bound), so cos^2 and the large-kappa von Mises measures are ruled out by');
  console.log('       that alone, quite apart from not matching -cos.');
  console.log('   => the requirement is a CONJUNCTION: reproduce -cos(D) for all D AND never exceed 2*sqrt(2).');
  console.log('      Of everything scanned here, only the engineered measure passes both.');
}
console.log('');

console.log('============================================================');
console.log('B: HOW SPECIAL IS IT?  (deformation and sensitivity)');
console.log('============================================================');
{
  const UNI = norm(PS.map(() => 1));
  console.log('   (B1) admixture:  p_eps = (1-eps)*engineered + eps*uniform');
  console.log('        eps      max dev from -cos      CHSH S');
  for (const eps of [0, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 1]) {
    const p = ENGINEERED.map((v, i) => (1 - eps) * v + eps * UNI[i]);
    const r = report(p, '');
    console.log(`        ${eps.toFixed(2)}     ${r.worst.toExponential(3).padStart(14)}      ${r.S.toFixed(6)}`);
  }
  console.log('');
  console.log('   (B2) family:  p ~ |sin 2psi|^p, normalised');
  console.log('        p        max dev from -cos      CHSH S');
  for (const pw of [0.5, 0.75, 0.9, 1.0, 1.1, 1.25, 1.5, 2]) {
    const p = norm(PS.map((t) => Math.pow(Math.abs(Math.sin(2 * t)), pw)));
    const r = report(p, '');
    console.log(`        ${pw.toFixed(2)}     ${r.worst.toExponential(3).padStart(14)}      ${r.S.toFixed(6)}`);
  }
  console.log('');
  console.log('   (B3) family:  p ~ 1 + A*cos(2*psi),  scan the amplitude A (2^0 .. 2^2)');
  console.log('        A        max dev from -cos      CHSH S');
  let bestA = { A: 0, S: -Infinity };
  for (const A of [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3]) {
    const p = norm(PS.map((t) => Math.max(0, 1 + A * Math.cos(2 * t))));
    const r = report(p, '');
    if (r.S > bestA.S) bestA = { A, S: r.S };
    console.log(`        ${A.toFixed(2)}     ${r.worst.toExponential(3).padStart(14)}      ${r.S.toFixed(6)}`);
  }
  console.log('');
  console.log('   => the required measure is a POINT in measure space: a 1% admixture of the uniform measure');
  console.log('      already shows, and the deviation grows linearly with it.');
  console.log('   => and note what the S column does: it is NOT maximised by the engineered measure');
  console.log('      (S = 2.828427 there, but 2.885 at p = 1.1, 3.0817 at p = 1.5, 2.9199 at A = 1.5).');
  console.log('      The engineered measure is special for matching -cos, not for making S large -- and the');
  console.log('      ones that overshoot 2.828427 are excluded by experiment.');
}
console.log('');

console.log('============================================================');
console.log('C: A DERIVATION ATTEMPT -- INVARIANT MEASURES OF SIMPLE DETERMINISTIC MAPS');
console.log('============================================================');
{
  console.log('   map family (standard circle map, local, one step):  psi -> psi + w + K*sin(m*psi)  mod 2pi');
  console.log('   measured by ORBIT AVERAGE (no histogram noise):  E(D) = - mean_n s(psi_n - D/2) s(psi_n + D/2)');
  console.log('');
  const M = 40000, TRANSIENT = 4000;
  function orbitE(psi0, w, K, m, dDeg) {
    const h = (dDeg * Math.PI) / 180 / 2;
    let psi = psi0, acc = 0, cnt = 0;
    for (let n = 0; n < TRANSIENT + M; n++) {
      psi = psi + w + K * Math.sin(m * psi);
      psi = ((psi % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
      if (n >= TRANSIENT) { acc += sgn(psi - h) * sgn(psi + h); cnt++; }
    }
    return -acc / cnt;
  }
  const DS = []; for (let d = 0; d <= 180; d += 5) DS.push(d);
  let best = { dev: Infinity, w: null, K: null, m: null, S: -Infinity };
  let bestS = { S: -Infinity, w: null, K: null, m: null, dev: Infinity };
  for (const m of [1, 2, 3]) for (const w of [0.1, 0.3, 0.5, 0.7, 1.0, Math.PI / 2, 2.0]) for (const K of [0.2, 0.5, 1.0, 2.0]) {
    let dev = 0;
    for (const d of DS) dev = Math.max(dev, Math.abs(orbitE(0.1234, w, K, m, d) - TGT(d)));
    const e = (d) => orbitE(0.1234, w, K, m, d);
    const S = Math.abs(e(-45) + e(45) + e(45) - e(135));
    if (dev < best.dev) best = { dev, w, K, m, S };
    if (S > bestS.S) bestS = { S, w, K, m, dev };
  }
  console.log(`   best match to -cos : dev = ${best.dev.toExponential(3)}   at w=${best.w.toFixed(3)}, K=${best.K}, m=${best.m}   (its S = ${best.S.toFixed(6)})`);
  console.log(`   best CHSH value    : S = ${bestS.S.toFixed(6)}   at w=${bestS.w.toFixed(3)}, K=${bestS.K}, m=${bestS.m}   (its dev = ${bestS.dev.toExponential(3)})`);
  console.log('');
  console.log('   NOTE WHAT THE SECOND LINE SAYS: the simple local maps are NOT too weak.  A near-rigid');
  console.log('   rotation (small K) makes psi ~ w*n, i.e. an almost deterministic phase, and that reaches');
  console.log('   S = 4.000000 -- the PR-box maximum, which EXPERIMENT EXCLUDES.  So simple maps fail in both');
  console.log('   directions at once: too strong for the bound, and nowhere near -cos in shape.');
  console.log('   (A sample, not a theorem: "searched and empty".)');
}
console.log('');

console.log('============================================================');
console.log('D: THE HONEST CATCH, MEASURED -- A MAP WITH THE REQUIRED MEASURE DOES EXIST, BY CONSTRUCTION');
console.log('============================================================');
{
  // Conjugate the doubling map by a CDF, and get the invariant measure from a proper Ulam transfer
  // operator: each bin's mass is spread over the images of SUBSAMPLES of the bin (using bin centres
  // alone collapses onto periodic orbits -- that was the bug in my previous two attempts).
  // Self-check first: the plain doubling map must give the UNIFORM density.
  const B = 512, SUB = 16, ITERS = 2000;
  const dv = (2 * Math.PI) / B;
  const PSB = Array.from({ length: B }, (_, i) => (i + 0.5) * dv);
  const engB = norm(PSB.map((t) => Math.abs(Math.sin(2 * t)) / 4));
  const cdfOfB = (p) => {
    const c = new Array(B + 1).fill(0);
    for (let i = 0; i < B; i++) c[i + 1] = c[i] + p[i];
    for (let i = 0; i <= B; i++) c[i] /= c[B];
    return c;
  };
  const cdfAtB = (c, psi) => {
    const x = (((psi % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI)) / (2 * Math.PI) * B;
    const i = Math.min(B - 1, Math.floor(x)), f = x - i;
    return c[i] + f * (c[i + 1] - c[i]);
  };
  const invCdfB = (c, u) => {
    u = ((u % 1) + 1) % 1;
    let lo = 0, hi = B;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (c[mid] < u) lo = mid + 1; else hi = mid; }
    const i = Math.max(1, lo), c0 = c[i - 1], c1 = c[i];
    const f = c1 > c0 ? (u - c0) / (c1 - c0) : 0;
    return (i - 1 + f) * dv;
  };
  const conjugate = (p) => { const c = cdfOfB(p); return (psi) => invCdfB(c, 2 * cdfAtB(c, psi)); };
  const doubling = (psi) => (2 * psi) % (2 * Math.PI);
  function invariantDensity(f) {
    const binOf = (x) => Math.min(B - 1, Math.floor((((x % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI)) / dv));
    const tgt = new Int32Array(B * SUB);
    for (let i = 0; i < B; i++) for (let k = 0; k < SUB; k++) tgt[i * SUB + k] = binOf(f(PSB[i] + (k + 0.5 - SUB / 2) * dv / SUB));
    let rho = new Float64Array(B).fill(1 / B);
    for (let it = 0; it < ITERS; it++) {
      const next = new Float64Array(B);
      const w = 1 / SUB;
      for (let i = 0; i < B; i++) { const m = rho[i] * w; if (m === 0) continue; for (let k = 0; k < SUB; k++) next[tgt[i * SUB + k]] += m; }
      let s = 0; for (let i = 0; i < B; i++) s += next[i];
      if (!(s > 0)) break;
      for (let i = 0; i < B; i++) next[i] /= s;
      rho = next;
    }
    return Array.from(rho);
  }
  const E_ofB = (p, dDeg) => { const h = (dDeg * Math.PI) / 180 / 2; let e = 0; for (let i = 0; i < B; i++) e += p[i] * sgn(PSB[i] - h) * sgn(PSB[i] + h); return -e; };
  const DS = []; for (let d = 0; d <= 180; d += 5) DS.push(d);
  const uniB = norm(PSB.map(() => 1));
  console.log(`   B = ${B} bins, ${SUB} subsamples per bin, ${ITERS} transfer-operator iterations.`);
  console.log('');
  console.log('   map                                  L1 to its intended density    max dev from -cos      S');
  const rows = [
    ['plain doubling map (self-check)', doubling, uniB],
    ['conjugate by the CDF of |sin 2psi|/4', conjugate(engB), engB],
    ['conjugate by the CDF of uniform', conjugate(uniB), uniB],
  ];
  for (const [label, f, intended] of rows) {
    const rho = invariantDensity(f);
    let l1 = 0; for (let i = 0; i < B; i++) l1 += Math.abs(rho[i] - intended[i]);
    l1 /= 2;
    let dev = 0; for (const d of DS) dev = Math.max(dev, Math.abs(E_ofB(rho, d) - TGT(d)));
    const e = (d) => E_ofB(rho, d);
    const S = Math.abs(e(-45) + e(45) + e(45) - e(135));
    console.log(`   ${label.padEnd(36)} ${l1.toExponential(3).padStart(18)}          ${dev.toExponential(3).padStart(12)}   ${S.toFixed(6)}`);
  }
  console.log('');
  console.log('   => the self-check passes (doubling gives the uniform density), so the machinery is sound, and');
  console.log('      the engineered conjugacy DOES give a deterministic map whose invariant measure is');
  console.log('      |sin 2psi|/4, reproducing -cos(D) and landing on 2*sqrt(2).');
  console.log('      BUT the measure is EXACTLY whatever goes into the conjugacy: choosing the CDF of');
  console.log('      |sin 2psi|/4 is the same reverse-engineering as round 73, one level down -- the density');
  console.log('      moved from "a distribution I wrote" to "a map I built to have that distribution".');
  console.log('      So a rule with the right measure EXISTS; what does not exist is an independently');
  console.log('      motivated reason to pick that rule rather than another.');
}
console.log('');

console.log('============================================================');
console.log('E: NO-TARGET AUDIT, AND VERDICT');
console.log('============================================================');
{
  console.log('   audit: the two rules of this round are');
  console.log('        theta = (a+b)/2 + psi                       <- contains a and b, but NOT their difference,');
  console.log('                                                       and no trigonometric function of it');
  console.log('        A = r*sign(cos(theta-a)), B = -r*sign(cos(theta-b))   <- purely local projections');
  console.log('   the target -cos(D) appears ONLY inside this script\'s comparison columns, never inside E_of().');
  console.log('   (check with:  grep -n "TGT" e4_generic.js   -> only in reporting lines)');
  console.log('');
  console.log('   VERDICT:');
  console.log('     1. with the two geometric rules fixed, E(D) is a FUNCTIONAL OF THE SOURCE MEASURE ALONE;');
  console.log('     2. passing the physical test means BOTH: reproduce -cos(D) for all D, and never exceed');
  console.log('        2*sqrt(2).  The engineered measure |sin 2psi|/4 does both;');
  console.log('     3. natural closed forms do neither reliably: they scatter S from 1.0997 up to 2.9003, so');
  console.log('        several of them (cos^2, large-kappa von Mises, p = 1.1 or 1.5, A = 1.5) OVERSHOOT');
  console.log('        Tsirelson and are excluded by experiment outright;');
  console.log('     4. a sampled family of simple local circle maps does not reproduce -cos, and reaches');
  console.log('        S = 4.000000 (the PR maximum) for near-rigid rotations -- wrong in both directions;');
  console.log('     5. BUT a deterministic map whose invariant measure IS |sin 2psi|/4 does exist and works:');
  console.log('        conjugate the doubling map by the CDF of |sin 2psi|/4 (verified here against a passing');
  console.log('        self-check).  So "a rule that generates the measure" EXISTS -- what is missing is any');
  console.log('        reason to choose THAT rule, since the measure is exactly whatever the conjugacy encodes.');
  console.log('');
  console.log('   SO THE HONEST STATE OF E4 AFTER ROUNDS 72-74:');
  console.log('     * the quantum correlations CAN be had from local, deterministic, no-signalling models,');
  console.log('       exactly and with zero free parameters .................... YES (rounds 72-73)');
  console.log('     * a DETERMINISTIC RULE whose invariant measure is the required one .... YES (part D),');
  console.log('       but it is built from the required CDF, so the target is still an input;');
  console.log('     * a rule that produces it for INDEPENDENT reasons (locality, symmetry, simplicity) ...');
  console.log('       NOT FOUND (this round: natural densities and simple circle maps both fail)');
  console.log('     => E4 is therefore like every other ingredient in this project: WRITABLE, AND');
  console.log('        GENERATABLE ONLY BY CONSTRUCTION -- never derived from something else.');
}
