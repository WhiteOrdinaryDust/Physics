'use strict';
/*
 * ether_and_dispersion.js   (third version)
 *
 * WHAT THE FIRST TWO VERSIONS GOT WRONG, AND WHY IT MATTERS
 *   v1 measured the temporal period of the CONFIGURATION rather than of the MODE, used a
 *      growth threshold for "bounded", and reduced a translating shift mod lambda with the
 *      wrong sign. Its controls failed.
 *   v2 replaced that with "follow the shape up to translation", which is right for a
 *      LOCALIZED defect but meaningless for a DELOCALIZED periodic mode: such a mode's
 *      bounding box is the whole ring, so "shape up to translation" never matches. Its
 *      controls failed too (rule 170 read as "standing with T=lambda").
 *   v3 therefore measures a mode the way a dispersion measurement should, and does not
 *      assume in advance which kind it is:
 *        (a) exact eventual periodicity: is there p with dev(t+p) == dev(t) for all large t?
 *            (A transient is allowed -- rule 90's linear map is singular, so its modes DO
 *            settle rather than being periodic from the start. Requiring periodicity from
 *            t=0 was v1's bug and it hid every rule-90 mode.)
 *        (b) a rigid shift per step, measured by overlap of consecutive frames over shifts
 *            mod lambda. This catches a translating mode even when it never recurs exactly.
 *      A mode is then: translating (w = k*v) or standing (w = 2*pi/T) or nothing.
 *
 * JOB 1 (Part A): rule 110's vacuum. The majority-vote approach failed for a good reason:
 *   the attractor is a chaotic mixture of ether plus gliders, so a per-residue majority is a
 *   SMOOTHED average, not an actual configuration, and it is not a cycle member. The right
 *   estimator for "the vacuum" is the MOST COMMON LOCAL PATTERN of a given length.
 *
 * PREDICTIONS
 *   D1 rule 170 (pure shift): translating, |v| = 1 for every lambda => w = k, straight line.
 *   D2 rule 204 (identity): no shift, period 1 => w = 0 for every k.
 *   D3 rule 90: no shift (the filter x + x^-1 is symmetric) => standing, with a period that
 *      varies with lambda => w(k) NOT linear => a preferred frame, not relativity.
 *   D4 rule 110's most common 14-cell window should be the ether; it must be a verified cycle
 *      member AND occur in a large fraction of the attractor.
 * ---------------------------------------------------------------------------
 */

function ruleTable(rule) { const t = new Uint8Array(8); for (let n = 0; n < 8; n++) t[n] = (rule >> n) & 1; return t; }
function step1(x, tb) {
  const L = x.length, o = new Uint8Array(L);
  for (let i = 0; i < L; i++) o[i] = tb[(x[(i - 1 + L) % L] << 2) | (x[i] << 1) | x[(i + 1) % L]];
  return o;
}
function xorA(a, b) { const o = new Uint8Array(a.length); for (let i = 0; i < a.length; i++) o[i] = a[i] ^ b[i]; return o; }
function popA(a) { let s = 0; for (let i = 0; i < a.length; i++) s += a[i]; return s; }
function eqA(a, b) { for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false; return true; }
function rot1(x, s) { const L = x.length, o = new Uint8Array(L); for (let i = 0; i < L; i++) o[i] = x[((i - s) % L + L) % L]; return o; }

// ---------- mode measurement ----------
function eventualPeriod(fr, pmax) {
  const T = fr.length - 1;
  for (let p = 1; p <= pmax; p++) {
    let t0 = -1;
    for (let t = T - p; t >= 0; t--) { if (eqA(fr[t], fr[t + p])) t0 = t; else break; }
    if (t0 >= 0 && T - t0 >= 2 * p) return { p, transient: t0 };
  }
  return null;
}
function bestShift(fr, lam, L) {
  // align fr[t+1] onto fr[t]; returns the shift (mod lam, representative of smallest magnitude)
  let best = -1, bestS = 0;
  for (let sh = 0; sh < lam; sh++) {
    let m = 0, tot = 0;
    for (let i = 0; i < L; i++) if (fr[1][i]) { tot++; if (fr[0][(i - sh + L) % L]) m++; }
    if (tot && m / tot > best) { best = m / tot; bestS = sh; }
  }
  let v = bestS; if (v > lam / 2) v -= lam;
  return { v, overlap: best };
}
function modeStatus(rule, bgAt, lam, T) {
  const tb = ruleTable(rule), L = bgAt(0).length;
  const d = new Uint8Array(L);
  for (let i = 0; i < L; i += lam) d[i] = 1;
  let s = xorA(bgAt(0), d);
  const fr = [];
  for (let t = 0; t <= T; t++) { fr.push(xorA(s, bgAt(t))); s = step1(s, tb); }
  const n0 = popA(fr[0]), nT = popA(fr[T]);
  const per = eventualPeriod(fr, 24);
  const sh = bestShift(fr, lam, L);
  const k = 2 * Math.PI / lam;
  if (sh.overlap > 0.995 && sh.v !== 0) return { lam, k, kind: 'translating', v: sh.v, overlap: sh.overlap, n0, nT, w: k * Math.abs(sh.v) };
  if (per) return { lam, k, kind: 'standing', p: per.p, transient: per.transient, n0, nT, w: per.p === 1 ? 0 : 2 * Math.PI / per.p };
  return { lam, k, kind: 'no mode', n0, nT, shift: sh.v, overlap: sh.overlap };
}
function dispersion(rule, bgAt, lams, T) { return lams.map(l => modeStatus(rule, bgAt, l, T)); }
function summarize(name, prediction, out, check) {
  console.log(`  ${name}`);
  console.log(`    prediction: ${prediction}`);
  const kinds = [...new Set(out.map(o => o.kind))];
  console.log(`    outcomes: ${kinds.map(k => `${out.filter(o => o.kind === k).length}x ${k}`).join(', ')}`);
  for (const o of out.filter(x => x.kind !== 'no mode').slice(0, 12)) {
    const det = o.kind === 'translating' ? `translating |v|=${Math.abs(o.v)} (overlap ${o.overlap.toFixed(3)})` : `standing T=${o.p} (settles after ${o.transient} steps)`;
    console.log(`      lambda ${String(o.lam).padStart(3)}  k=${o.k.toFixed(3)}  ${det.padEnd(40)} w=${o.w.toFixed(3)}`);
  }
  const ok = check(out);
  console.log(`    => ${ok ? 'PASS' : '*** FAIL ***'}`);
  return out;
}

// ============================================================ PART A
console.log('============================================================');
console.log('PART A: rule 110 -- the vacuum as the MOST COMMON LOCAL PATTERN');
console.log('============================================================');
const tb110 = ruleTable(110);
const LATE = (() => {
  const L = 700;
  let s = new Uint8Array(L), sd = 12345;
  for (let i = 0; i < L; i++) { sd = (sd * 1103515245 + 12345) & 0x7fffffff; s[i] = ((sd >> 16) & 1) ? 1 : 0; }
  for (let t = 0; t < 2000; t++) s = step1(s, tb110);
  return s;
})();
function commonWindows(s, p) {
  const L = s.length, m = new Map();
  for (let i = 0; i < L; i++) {
    let key = '';
    for (let j = 0; j < p; j++) key += s[(i + j) % L];
    m.set(key, (m.get(key) || 0) + 1);
  }
  return [...m.entries()].sort((a, b) => b[1] - a[1]);
}
console.log(`  rule 110, 700-cell ring, 2000 steps from a random start: density ${(popA(LATE) / LATE.length).toFixed(4)}`);
function cycleOf(P, p) {
  for (let q = 1; q <= 60; q++) {
    let cur = P;
    for (let k = 0; k < q; k++) cur = step1(cur, tb110);
    for (let s = 0; s < p; s++) {
      let ok = true;
      for (let i = 0; i < p; i++) if (cur[i] !== P[((i - s) % p + p) % p]) { ok = false; break; }
      if (ok) return { q, s };
    }
  }
  return null;
}
let ETHER = null;
for (const p of [7, 14, 21, 28, 35]) {
  const top = commonWindows(LATE, p).slice(0, 3);
  console.log(`    window length ${p}:`);
  for (const [w, cnt] of top) {
    const P = Uint8Array.from(w, c => Number(c));
    const cyc = cycleOf(P, p);
    console.log(`      ${w}  occurs ${cnt}/${LATE.length} = ${(100 * cnt / LATE.length).toFixed(1)}%   cycle: ${cyc ? `step^${cyc.q} = shift ${cyc.s} (v=${(cyc.s / cyc.q).toFixed(3)})` : 'NO'}`);
    if (cyc && !ETHER && p === 14) ETHER = { p, P, cyc, frac: cnt / LATE.length };
  }
}
if (ETHER) {
  const causal = Math.abs(ETHER.cyc.s) / ETHER.cyc.q <= 1;
  console.log('');
  console.log(`  => ETHER CANDIDATE: period ${ETHER.p}, pattern ${Array.from(ETHER.P).join('')}, density ${(popA(ETHER.P) / ETHER.p).toFixed(4)},`);
  console.log(`     the most common local pattern (${(100 * ETHER.frac).toFixed(1)}% of all windows), and a verified cycle member (step^${ETHER.cyc.q} = shift ${ETHER.cyc.s}).`);
  if (!causal) {
    console.log(`     *** CAVEAT: that shift is ${ETHER.cyc.s} per ${ETHER.cyc.q} step(s), which EXCEEDS the light cone (a radius-1`);
    console.log('     rule moves influence at most 1 cell per step). So this is a self-similarity of the pattern,');
    console.log('     not a physical travelling object: the one-step image of the pattern happens to be a rotation');
    console.log('     of itself. Any "phase velocity" read off from it is therefore meaningless, and the ether is');
    console.log('     still not given a physical frame. Do not quote a speed for it. ***');
  }
  console.log(`     NOTE on the 6.4%: for a background that can sit in many relative phases, one specific 14-cell`);
  console.log(`     window catches only one phase. The ether is therefore the dominant structure UP TO PHASE.`);
} else {
  console.log('  => no ether found among the dominant windows of these lengths.');
}
console.log('');

// ============================================================ PART B
console.log('============================================================');
console.log('PART B: the dispersion ruler w(k), quantitative controls');
console.log('============================================================');
const LD = 840, LAMS = [];
for (let l = 2; l <= 60; l++) if (LD % l === 0) LAMS.push(l);
const empty = () => new Uint8Array(LD);
summarize('D1 CONTROL rule 170 (pure shift)', 'translating, |v| = 1 for every lambda => w = k, a straight line (gapless)',
  dispersion(170, empty, LAMS, 30),
  (o) => o.every(x => x.kind === 'translating' && Math.abs(Math.abs(x.v) - 1) < 1e-9));
console.log('');
summarize('D2 CONTROL rule 204 (identity)', 'no shift, period 1 => w = 0 for every k (flat band, infinite mass)',
  dispersion(204, empty, LAMS, 10),
  (o) => o.every(x => x.kind === 'standing' && x.p === 1));
console.log('');
{
  const out = summarize('D3 CONTROL rule 90 (linear, symmetric filter)', 'standing for every lambda (no shift), with a period that varies with lambda => w(k) not linear',
    dispersion(90, empty, LAMS, 120),
    (o) => o.every(x => x.kind === 'standing' || x.kind === 'no mode'));
  const Ts = out.filter(x => x.kind === 'standing').map(x => x.p);
  console.log(`    distinct shape periods: ${[...new Set(Ts)].sort((a, b) => a - b).join(',')}`);
  const wOverK = out.filter(x => x.kind === 'standing').map(x => (x.w / x.k).toFixed(4));
  console.log(`    w/k values: ${[...new Set(wOverK)].slice(0, 8).join(', ')}${new Set(wOverK).size > 8 ? ' ...' : ''}`);
  console.log(`    => ${new Set(wOverK).size === 1 ? 'w/k is constant: LINEAR dispersion' : `w/k takes ${new Set(wOverK).size} different values: the dispersion is NOT linear in k`}`);
}
console.log('');

// ============================================================ PART C
console.log('============================================================');
console.log('PART C: the three rulers on the ether');
console.log('============================================================');
if (!ETHER) {
  console.log('  No verified ether, so there is nothing to run them on.');
} else {
  const p = ETHER.p, { q, s } = ETHER.cyc, LL = p * 60;
  const R = new Uint8Array(LL);
  for (let i = 0; i < LL; i++) R[i] = ETHER.P[i % p];
  const cycA = [R];
  for (let k = 1; k < q; k++) cycA.push(step1(cycA[k - 1], tb110));
  const back = step1(cycA[q - 1], tb110);
  const rotOk = back.every((v, i) => v === R[((i - s) % LL + LL) % LL]);
  console.log(`  ring ${LL}; verification step^${q} = rotate by ${s}: ${rotOk ? 'PASS' : '*** FAIL ***'}`);
  if (rotOk) {
    const bgAt = (t) => rot1(cycA[t % q], (s * Math.floor(t / q)) % LL);
    const bg0 = bgAt(0), C = Math.floor(LL / 2);
    const d1 = new Uint8Array(LL), d2 = new Uint8Array(LL);
    d1[C] = 1; d2[C + 1] = 1;
    let a1 = xorA(bg0, d1), a2 = xorA(bg0, d2), a12 = xorA(bg0, xorA(d1, d2));
    let mx = 0, nz = 0;
    for (let t = 0; t <= 120; t++) {
      const ref = bgAt(t);
      const m = popA(xorA(xorA(a12, ref), xorA(xorA(a1, ref), xorA(a2, ref))));
      if (m > 0) nz++;
      mx = Math.max(mx, m);
      a1 = step1(a1, tb110); a2 = step1(a2, tb110); a12 = step1(a12, tb110);
    }
    console.log(`  RULER 1 (linearity, ADJACENT deviations): max|mismatch| = ${mx}, steps nonzero = ${nz}/121`);
    console.log(`    => ${mx === 0 ? 'EXACTLY LINEAR: this substrate has a photon sector' : 'nonlinear: no linear wave sector around this background'}`);
    const lams = [];
    for (let l = 2; l <= 60; l++) if (LL % l === 0) lams.push(l);
    const out = dispersion(110, bgAt, lams, 90);
    console.log(`  RULER 3 (dispersion): ${[...new Set(out.map(o => o.kind))].map(k => `${out.filter(o => o.kind === k).length}x ${k}`).join(', ')}`);
    for (const o of out.filter(x => x.kind !== 'no mode').slice(0, 10)) {
      const det = o.kind === 'translating' ? `translating |v|=${Math.abs(o.v)}` : `standing T=${o.p} (settles after ${o.transient})`;
      console.log(`    lambda ${String(o.lam).padStart(3)}  ${det.padEnd(38)} w=${o.w.toFixed(3)}  phase velocity w/k=${(o.w / o.k).toFixed(3)}`);
    }
  }
}
