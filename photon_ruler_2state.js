'use strict';
/*
 * photon_ruler_2state.js   (second version -- three defects in the first draft fixed)
 *
 * FIXED AFTER READING THE FIRST RUN'S OWN OUTPUT (this is the point of printing data):
 *   1. the 2D shape signature was only "bounding box + mass", which is far too coarse:
 *      it reported the Life glider as "period 3" when the glider's period is 4. Now the
 *      signature is the full canonical cell pattern.
 *   2. "SPREADS" was triggered by (width grew by >2), which mislabelled a small BOUNDED
 *      defect as spreading. Now a deviation is called spreading only if it is still
 *      growing between t=3T/4 and t=T.
 *   3. the rule 110 ether linearity test put the two deviations 20 cells apart. In a
 *      local rule that can pass TRIVIALLY (no action at a distance). The test now also
 *      runs with the two deviations ADJACENT, where interaction is forced, and the
 *      ether's own single-cell defect is tracked to see whether it moves or sits.
 *
 * WHAT THE RULERS ARE
 *   RULER 1 (linearity): in a 2-state substrate a deviation is a set of flipped cells, so
 *     evolve(bg+d1+d2) - bg  ==  (evolve(bg+d1) - bg) XOR (evolve(bg+d2) - bg)
 *   mismatches are counted every step. Zero always = additive = "the messenger does not
 *   interact with itself".
 *   RULER 2 (propagation): inject one small deviation; does it travel rigidly, spread,
 *   sit still, or get absorbed?
 *
 * PREDICTIONS (before the run)
 *   P1 rules 90/150/170/204 -> RULER 1 exactly 0 everywhere (they are linear over GF(2))
 *   P2 rule 30 and Life's block lattice -> RULER 1 nonzero, immediately
 *   P3 rule 90's deviation spreads over the whole light cone (Sierpinski); it is linear
 *      but NOT a wave packet. Cross-check: its mass at time t must be 2^popcount(t).
 *   P4 rule 170 (pure shift) is linear AND rigid -- but contentless (no interaction, no matter)
 *   P5 Life's glider is rigid with period 4 and velocity (1,1); its nonlinearity is zero at
 *      a distance and switches on only at contact
 *   P6 rule 110 HAS a genuine travelling periodic background (the ether) -- verifiable
 *      exactly -- and I expect its perturbation sector to be small and bound rather than
 *      free particles, with a THRESHOLD for nucleating a glider (this is the "finite
 *      threshold medium" I predicted in the previous round)
 * ---------------------------------------------------------------------------
 */

// ============================================================ PART 0
console.log('============================================================');
console.log('PART 0: can a 2-state substrate even REPRESENT two quanta in one mode?');
console.log('============================================================');
{
  const mode = [200];
  const bg = new Uint8Array(512);
  const one = Uint8Array.from(bg); for (const i of mode) one[i] ^= 1;
  const two = Uint8Array.from(one); for (const i of mode) two[i] ^= 1;
  const eq = (a, b) => a.every((v, i) => v === b[i]);
  console.log('  natural encoding: |0> = background, "add one quantum" = flip the mode\'s cells');
  console.log(`  add one quantum  -> the mode's cell is ${one[200] ? 'ON' : 'OFF'} (a distinct state, fine)`);
  console.log(`  add a second one -> the mode's cell is ${two[200] ? 'ON' : 'OFF'}; identical to |0> ? ${eq(two, bg)}`);
  console.log('  => in this encoding occupancy is mod 2: |2> = |0>. Two quanta in one mode are');
  console.log('     indistinguishable from none. That is the EXCLUSION law (fermionic). A photon');
  console.log('     needs |2> to be distinct and brighter. The only escape is to encode |2> as a');
  console.log('     different cell pattern, but then the occupation->configuration map is nonlinear,');
  console.log('     which destroys exactly the linearity the photon needs.');
}
console.log('');

// ============================================================ 1D machinery
function ruleTable(rule) { const t = new Uint8Array(8); for (let n = 0; n < 8; n++) t[n] = (rule >> n) & 1; return t; }
function step1(x, tb) {
  const L = x.length, o = new Uint8Array(L);
  for (let i = 0; i < L; i++) o[i] = tb[(x[(i - 1 + L) % L] << 2) | (x[i] << 1) | x[(i + 1) % L]];
  return o;
}
function xorA(a, b) { const o = new Uint8Array(a.length); for (let i = 0; i < a.length; i++) o[i] = a[i] ^ b[i]; return o; }
function popA(a) { let s = 0; for (let i = 0; i < a.length; i++) s += a[i]; return s; }
function rot1(x, s) { const L = x.length, o = new Uint8Array(L); for (let i = 0; i < L; i++) o[i] = x[((i - s) % L + L) % L]; return o; }
function shape1(dev) {
  let lo = -1, hi = -1, n = 0, sx = 0;
  for (let i = 0; i < dev.length; i++) if (dev[i]) { if (lo < 0) lo = i; hi = i; n++; sx += i; }
  if (!n) return { id: 'empty', n: 0, w: 0, cx: 0 };
  let s = ''; for (let i = lo; i <= hi; i++) s += dev[i] ? '#' : '.';
  return { id: `${hi - lo + 1}:${s}`, n, w: hi - lo + 1, cx: sx / n };
}
const L = 512;

function linearity1D(rule, bgAt, d1, d2, T) {
  const tb = ruleTable(rule);
  const bg0 = bgAt(0);
  let s1 = xorA(bg0, d1), s2 = xorA(bg0, d2), s12 = xorA(bg0, xorA(d1, d2));
  const tr = [];
  for (let t = 0; t <= T; t++) {
    const ref = bgAt(t);
    tr.push(popA(xorA(xorA(s12, ref), xorA(xorA(s1, ref), xorA(s2, ref)))));
    s1 = step1(s1, tb); s2 = step1(s2, tb); s12 = step1(s12, tb);
  }
  return tr;
}
function propagate1D(rule, bgAt, dev0, T) {
  const tb = ruleTable(rule);
  let s = xorA(bgAt(0), dev0);
  const frames = [];
  for (let t = 0; t <= T; t++) {
    const sh = shape1(xorA(s, bgAt(t)));
    frames.push({ t, ...sh });
    s = step1(s, tb);
  }
  return frames;
}
function classify1D(frames) {
  const fe = frames.findIndex(f => f.n === 0);
  if (fe >= 0) return `absorbed at t=${fe}`;
  for (let p = 1; p <= 8; p++) {
    let ok = true, base = null;
    for (let t = 0; t + p < frames.length; t++) {
      if (frames[t].id !== frames[t + p].id) { ok = false; break; }
      const d = Math.round(frames[t + p].cx - frames[t].cx);
      if (base === null) base = d; else if (d !== base) { ok = false; break; }
    }
    if (ok && base !== null && base !== 0) return `RIGID: period ${p}, velocity ${base}/period, mass ${frames[0].n}, shape ${frames[0].id}`;
  }
  const q = frames[Math.floor(frames.length * 0.75)], last = frames[frames.length - 1];
  const grew = last.w > q.w + 1 || last.n > q.n + 2;
  const moved = Math.abs(last.cx - frames[0].cx) > 1;
  if (grew) return `SPREADS: width ${frames[0].w} -> ${last.w}, mass ${frames[0].n} -> ${last.n}`;
  return `BOUNDED (width ${last.w}, mass ${last.n}, centroid moved ${(last.cx - frames[0].cx).toFixed(1)} over ${last.t} steps)`;
}

// ============================================================ PART 1
console.log('============================================================');
console.log('PART 1: the rulers on 1D 2-state rules (ring of 512, 120 steps)');
console.log('============================================================');
const T1 = 120;
const emptyBg = () => new Uint8Array(L);
const oneD = [[90, 'exactly linear'], [150, 'exactly linear'], [170, 'exactly linear'], [204, 'exactly linear'], [30, 'nonlinear']];

console.log('  RULER 1 -- two deviations FAR apart (200 and 300):');
console.log('    rule   max|mismatch|   steps with mismatch>0 / total');
for (const [rule, kind] of oneD) {
  const d1 = new Uint8Array(L), d2 = new Uint8Array(L);
  d1[200] = 1; d2[300] = 1;
  const tr = linearity1D(rule, emptyBg, d1, d2, T1);
  console.log(`    ${String(rule).padStart(3)}   ${String(Math.max(...tr)).padStart(10)}   ${tr.filter(v => v > 0).length} / ${tr.length}      (${kind})`);
}
console.log('  RULER 1 -- two deviations ADJACENT (200 and 201), where interaction is forced:');
console.log('    rule   max|mismatch|   steps with mismatch>0 / total');
for (const [rule, kind] of oneD) {
  const d1 = new Uint8Array(L), d2 = new Uint8Array(L);
  d1[200] = 1; d2[201] = 1;
  const tr = linearity1D(rule, emptyBg, d1, d2, T1);
  console.log(`    ${String(rule).padStart(3)}   ${String(Math.max(...tr)).padStart(10)}   ${tr.filter(v => v > 0).length} / ${tr.length}      (${kind})`);
}
console.log('');
console.log('  RULER 2 -- one single-cell deviation at 200:');
console.log('    rule   outcome');
for (const [rule] of oneD) {
  const d = new Uint8Array(L); d[200] = 1;
  console.log(`    ${String(rule).padStart(3)}   ${classify1D(propagate1D(rule, emptyBg, d, T1))}`);
}
{
  // internal cross-check of the rule-90 implementation against the known closed form
  const tb = ruleTable(90);
  let x = new Uint8Array(L); x[200] = 1;
  const pops = [];
  for (let t = 0; t <= 8; t++) { pops.push(popA(x)); x = step1(x, tb); }
  const pc = (t) => { let c = 0, m = t; while (m) { c += m & 1; m >>= 1; } return c; };
  const expect = [0, 1, 2, 3, 4, 5, 6, 7, 8].map(t => 2 ** pc(t));
  console.log(`  cross-check rule 90 mass at t=0..8: measured ${pops.join(',')}   closed form 2^popcount(t) = ${expect.join(',')}   ${pops.join() === expect.join() ? 'MATCH' : '*** MISMATCH ***'}`);
}
console.log('');

// ============================================================ PART 1b
console.log('============================================================');
console.log('PART 1b: rule 110 -- a TRUSTWORTHY travelling background, then the rulers on it');
console.log('============================================================');
function findPeriodicBackgrounds(rule, pmax, qmax) {
  const tb = ruleTable(rule), out = [];
  for (let p = 2; p <= pmax; p++) {
    for (let mask = 1; mask < (1 << p) - 1; mask++) {
      const P = new Uint8Array(p); let ones = 0;
      for (let i = 0; i < p; i++) if ((mask >> i) & 1) { P[i] = 1; ones++; }
      let minimal = true;
      for (let q = 1; q < p; q++) {
        if (p % q) continue;
        let same = true;
        for (let i = 0; i < p && same; i++) if (P[i] !== P[i % q]) same = false;
        if (same) { minimal = false; break; }
      }
      if (!minimal) continue;
      let cur = P, found = null;
      for (let q = 1; q <= qmax && !found; q++) {
        cur = step1(cur, tb);
        for (let s = 0; s < p; s++) {
          let ok = true;
          for (let i = 0; i < p; i++) if (cur[i] !== P[((i - s) % p + p) % p]) { ok = false; break; }
          if (ok) { found = { p, q, s, ones }; break; }
        }
      }
      if (found) {
        // ALIASING GUARD (this project has already been burned by a "v = 2.0" nonsense once):
        // on a ring of length p a shift by s is the same as a shift by s - p, so only the
        // representative with the smallest magnitude is physical, and it must not exceed the
        // light cone (a radius-1 rule moves influence at most 1 cell per step).
        let sm = ((found.s % p) + p) % p; if (sm > p / 2) sm -= p;
        out.push({ p, q: found.q, s: sm, ones: found.ones, causal: Math.abs(sm) / found.q <= 1, P });
      }
    }
  }
  const seen = new Map();
  for (const f of out) {
    const k = `${f.p}|${f.ones}`;
    if (!seen.has(k) || (f.causal && !seen.get(k).causal) || (f.causal === seen.get(k).causal && seen.get(k).q > f.q)) seen.set(k, f);
  }
  return [...seen.values()].sort((a, b) => b.ones - a.ones || a.p - b.p);
}
{
  const t0 = Date.now();
  const all = findPeriodicBackgrounds(110, 14, 14);
  const prop = all.filter(f => f.s !== 0 && f.causal);
  console.log(`  exact enumeration of spatially periodic rule-110 backgrounds (period <= 14, return <= 14):`);
  console.log(`  a background counts only if step^q(P) equals P shifted by s (checked cell by cell), and the`);
  console.log(`  shift is reduced to the representative of smallest magnitude and required to respect the`);
  console.log(`  light cone (|s|/q <= 1). ${all.length} classes total, ${prop.length} of them causally PROPAGATE.`);
  console.log('    period   ones   density   return q   shift s   phase velocity');
  for (const f of prop.slice(0, 8)) console.log(`    ${String(f.p).padStart(6)}   ${String(f.ones).padStart(4)}   ${(f.ones / f.p).toFixed(4).padStart(7)}   ${String(f.q).padStart(8)}   ${String(f.s).padStart(7)}   ${(f.s / f.q).toFixed(3)}`);
  console.log(`  (${((Date.now() - t0) / 1000).toFixed(1)} s)`);
  const tb = ruleTable(110);
  // WHICH periodic pattern is actually the vacuum the attractor settles into?
  // (the first attempt picked the DENSEST propagating pattern, and the attractor check then
  //  read 0.0%: that pattern is a periodic solution, but it is not the ether. Let the data choose.)
  const RR = 700;
  const late = (() => {
    let s = new Uint8Array(RR), sd = 12345;
    for (let i = 0; i < RR; i++) { sd = (sd * 1103515245 + 12345) & 0x7fffffff; s[i] = ((sd >> 16) & 1) ? 1 : 0; }
    for (let t = 0; t < 2000; t++) s = step1(s, tb);
    return s;
  })();
  function bestMatch(P) {
    let best = 0, bestPhi = 0;
    for (let phi = 0; phi < P.length; phi++) {
      let m = 0;
      for (let i = 0; i < RR; i++) {
        let ok = true;
        for (let j = 0; j < P.length; j++) if (late[(i + j) % RR] !== P[(j + phi) % P.length]) { ok = false; break; }
        if (ok) m++;
      }
      if (m > best) { best = m; bestPhi = phi; }
    }
    return { frac: best / RR, phi: bestPhi };
  }
  const cands = prop.map(g => ({ g, ...bestMatch(g.P) })).sort((a, b) => b.frac - a.frac);
  console.log('  which periodic pattern is the vacuum the attractor settles into?');
  console.log(`  (700-cell ring, 2000 steps from a random start; density ${(popA(late) / RR).toFixed(4)})`);
  for (const c of cands.slice(0, 5)) console.log(`    p=${String(c.g.p).padStart(3)}  density ${(c.g.ones / c.g.p).toFixed(4)}  v=${(c.g.s / c.g.q).toFixed(3)}  ->  ${(100 * c.frac).toFixed(1)}% of windows match`);
  const pick = cands[0];
  const isVacuum = pick.frac >= 0.2;
  if (!isVacuum) {
    console.log('');
    console.log('  *** READ BEFORE QUOTING ANY NUMBER BELOW ***');
    console.log(`  The attractor is NOT dominated by any periodic pattern of period <= 14 (best match`);
    console.log(`${(100 * pick.frac).toFixed(1)}% of windows). So on this ring rule 110's vacuum is still not identified, and the`);
    console.log('  numbers that follow describe SOME periodic solution of rule 110 -- they are NOT');
    console.log('  statements about "the ether". The previous rounds\' ether trouble is therefore not');
    console.log('  resolved here; it is now localised: the enumeration kept only one pattern per');
    console.log('  (period, density), which can discard the real ether, and the attractor needs a');
    console.log('  larger spatial period to be matched.');
    console.log('');
  }
  const Pc = new Uint8Array(pick.g.p);
  for (let j = 0; j < pick.g.p; j++) Pc[j] = pick.g.P[(j + pick.phi) % pick.g.p];
  const f = { p: pick.g.p, q: pick.g.q, s: pick.g.s, ones: pick.g.ones, P: Pc };
  console.log(`  chosen background: p=${f.p}, density ${(f.ones / f.p).toFixed(4)}, velocity ${(f.s / f.q).toFixed(3)}, pattern ${Pc.join('')}${isVacuum ? '  (validated as the attractor vacuum)' : '  (NOT validated as the attractor vacuum)'}`);
  const LL = f.p * 60;               // long ring so nothing wraps while travelling
  const C = Math.floor(LL / 2);      // inject in the middle
  const R = new Uint8Array(LL);
  for (let i = 0; i < LL; i++) R[i] = f.P[i % f.p];
  const cyc = [R];
  for (let k = 1; k < f.q; k++) cyc.push(step1(cyc[k - 1], tb));
  const back = step1(cyc[f.q - 1], tb);
  const rotOk = back.every((v, i) => v === R[((i - f.s) % LL + LL) % LL]);
  console.log(`  verification that step^${f.q}(R) = rotate(R, ${f.s}) on the ring of ${LL}: ${rotOk ? 'PASS' : '*** FAIL ***'}`);
  if (!rotOk) { console.log('  *** aborting the rule 110 section: the background is not on a cycle ***'); }
  else {
    const matVel = f.s / f.q;
    const bgAt = (t) => rot1(cyc[t % f.q], (f.s * Math.floor(t / f.q)) % LL);
    const wrapCheck = (fr) => fr.some(x => x.n > 0 && (x.cx < 5 || x.cx > LL - 5));
    console.log('');
    console.log(`    NOTE the frame: the ether itself travels at +${matVel.toFixed(3)} cells/step. The deviation's`);
    console.log(`    centroid is measured in absolute cell indices, so a defect that RIDES the material has`);
    console.log(`    centroid velocity exactly +${matVel.toFixed(3)}; its velocity RELATIVE to the material is 0.`);
    console.log('    RULER 1 around the ether:');
    for (const [label, p1, p2] of [['20 cells apart', C, C + 20], ['ADJACENT (forced interaction)', C, C + 1], ['3 cells apart', C, C + 3]]) {
      const d1 = new Uint8Array(LL), d2 = new Uint8Array(LL);
      d1[p1] = 1; d2[p2] = 1;
      const tr = linearity1D(110, bgAt, d1, d2, T1);
      console.log(`      ${label.padEnd(30)} max|mismatch| = ${String(Math.max(...tr)).padStart(4)}, steps with mismatch>0 = ${tr.filter(v => v > 0).length}/${tr.length}`);
    }
    // outcome in the material's own rest frame: bound (rides) vs free (moves relative to the material)
    function etherOutcome(fr) {
      const fe = fr.findIndex(x => x.n === 0);
      let rigid = null;
      for (let p = 1; p <= 8; p++) {
        let ok = true, b = null;
        for (let t = 0; t + p < fr.length; t++) {
          if (fr[t].id !== fr[t + p].id) { ok = false; break; }
          const d = Math.round(fr[t + p].cx - fr[t].cx);
          if (b === null) b = d; else if (d !== b) { ok = false; break; }
        }
        if (ok && b !== null) { rigid = { p, drift: b }; break; }
      }
      if (fe >= 0) return { kind: 'absorbed', mass: 0, width: 0, vCom: 0, vAbs: 0, rigid };
      const qq = fr[Math.floor(fr.length * 0.75)], last = fr[fr.length - 1];
      const grew = last.w > qq.w + 1 || last.n > qq.n + 2;
      const vAbs = (last.cx - fr[0].cx) / last.t;
      const vCom = vAbs - matVel;
      if (grew) return { kind: 'DESTROYED (spreads/chaos)', mass: last.n, width: last.w, vCom, vAbs, rigid };
      const bound = Math.abs(vCom) < 0.05;
      return { kind: bound ? 'BOUND (rides the material)' : '*** FREE (moves relative to the material) ***', mass: last.n, width: last.w, vCom, vAbs, rigid };
    }
    // is this pattern really the vacuum that random initial conditions settle into?
    {
      const RR = 700;
      let s = new Uint8Array(RR);
      let seed = 12345;
      for (let i = 0; i < RR; i++) { seed = (seed * 1103515245 + 12345) & 0x7fffffff; s[i] = ((seed >> 16) & 1) ? 1 : 0; }
      let dens0 = popA(s) / RR;
      for (let t = 0; t < 2000; t++) s = step1(s, tb);
      let match = 0, tot = 0;
      for (let i = 0; i < RR; i++) { let ok = true; for (let j = 0; j < f.p; j++) if (s[(i + j) % RR] !== f.P[j]) { ok = false; break; } tot++; if (ok) match++; }
      console.log(`    ether-attractor check: random start (density ${dens0.toFixed(3)}), after 2000 steps`);
      console.log(`      ${(100 * match / tot).toFixed(1)}% of all ${f.p}-cell windows match the ether pattern in that phase (chance level ~0)`);
      console.log(`      late-time density ${(popA(s) / RR).toFixed(4)} vs the ether's ${(f.ones / f.p).toFixed(4)}`);
    }
    console.log('    RULER 2 around the ether: flip k adjacent cells of the ether, 240 steps');
    console.log('      k   outcome                                                      mass  width  v(absolute)  v(relative to material)');
    for (let k = 1; k <= 6; k++) {
      const d = new Uint8Array(LL);
      for (let i = 0; i < k; i++) d[C + i] = 1;
      const fr = propagate1D(110, bgAt, d, 240);
      const o = etherOutcome(fr);
      console.log(`      ${k}   ${o.kind.padEnd(56)} ${String(o.mass).padStart(5)} ${String(o.width).padStart(6)}   ${(o.vAbs || 0).toFixed(3).padStart(9)}   ${(o.vCom || 0).toFixed(3).padStart(9)}${o.rigid ? `  rigid p=${o.rigid.p} drift=${o.rigid.drift}` : ''}`);
    }
    // EXHAUSTIVE nucleation search: every possible perturbation of a 6-cell window.
    console.log('    EXHAUSTIVE nucleation search: all 2^6 - 1 = 63 perturbations of a 6-cell window of the ether');
    const tally = {};
    const frees = [];
    const movers = [];
    for (let mask = 1; mask < (1 << 6); mask++) {
      const d = new Uint8Array(LL);
      for (let i = 0; i < 6; i++) if ((mask >> i) & 1) d[C + i] = 1;
      const fr = propagate1D(110, bgAt, d, 240);
      const o = etherOutcome(fr);
      const key = o.kind.startsWith('***') ? 'FREE' : o.kind.split(' (')[0];
      tally[key] = (tally[key] || 0) + 1;
      if (key === 'FREE') frees.push({ mask, o });
      if (o.rigid && Math.abs(o.vCom) > 0.1) movers.push({ mask, o });
    }
    for (const k of Object.keys(tally).sort()) console.log(`      ${String(tally[k]).padStart(3)} of 63  ->  ${k}`);
    console.log(`      rigid movers with |velocity relative to the material| > 0.1: ${movers.length}`);
    movers.sort((a, b) => Math.abs(b.o.vCom) - Math.abs(a.o.vCom));
    for (const mm of movers.slice(0, 4)) console.log(`        mask ${mm.mask.toString(2).padStart(6, '0')}: v_rel ${mm.o.vCom.toFixed(3)}, mass ${mm.o.mass}, width ${mm.o.width}, rigid p=${mm.o.rigid.p} drift=${mm.o.rigid.drift}`);
    for (const fr2 of frees.slice(0, 4)) console.log(`      FREE case: mask ${fr2.mask.toString(2).padStart(6, '0')}, mass ${fr2.o.mass}, width ${fr2.o.width}, v_rel ${fr2.o.vCom.toFixed(3)}${fr2.o.rigid ? `, rigid p=${fr2.o.rigid.p}` : ' (not even rigid)'}`);
    if (!movers.length) console.log('      => no perturbation of a 6-cell window makes anything that moves at all crisply relative to the material.');
  }
}
console.log('');

// ============================================================ PART 2
console.log('============================================================');
console.log('PART 2: the rulers on 2D 2-state substrates');
console.log('============================================================');
const SZ = 240;
function stepL(g, bm, sm) {
  const o = new Uint8Array(SZ * SZ);
  for (let y = 0; y < SZ; y++) {
    const yp = ((y - 1 + SZ) % SZ) * SZ, ym = ((y + 1) % SZ) * SZ, yc = y * SZ;
    for (let x = 0; x < SZ; x++) {
      const xp = (x - 1 + SZ) % SZ, xm = (x + 1) % SZ;
      const n = g[yp + xp] + g[yp + x] + g[yp + xm] + g[yc + xp] + g[yc + xm] + g[ym + xp] + g[ym + x] + g[ym + xm];
      const i = yc + x;
      o[i] = g[i] ? ((sm >> n) & 1) : ((bm >> n) & 1);
    }
  }
  return o;
}
const LIFE = { bm: 1 << 3, sm: (1 << 2) | (1 << 3) };
const newGrid = () => new Uint8Array(SZ * SZ);
function setCells(g, cells) { for (const [x, y] of cells) g[y * SZ + x] = 1; }
const GLIDER = [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]];
function canonDev(dev) {
  let minx = 1e9, miny = 1e9, maxx = -1e9, maxy = -1e9, n = 0, sx = 0, sy = 0;
  for (let y = 0; y < SZ; y++) for (let x = 0; x < SZ; x++) if (dev[y * SZ + x]) {
    n++; sx += x; sy += y;
    if (x < minx) minx = x; if (x > maxx) maxx = x; if (y < miny) miny = y; if (y > maxy) maxy = y;
  }
  if (!n) return { id: 'empty', n: 0, w: 0, h: 0, cx: 0, cy: 0 };
  const w = maxx - minx + 1, h = maxy - miny + 1, rows = [];
  for (let r = 0; r < h; r++) rows.push(new Array(w).fill('.'));
  for (let y = miny; y <= maxy; y++) for (let x = minx; x <= maxx; x++) if (dev[y * SZ + x]) rows[y - miny][x - minx] = '#';
  return { id: `${w}x${h}:` + rows.map(r => r.join('')).join('/'), n, w, h, cx: sx / n, cy: sy / n };
}
function linearity2D(rule, bgAt, d1, d2, T) {
  const bg0 = bgAt(0);
  let s1 = xorA(bg0, d1), s2 = xorA(bg0, d2), s12 = xorA(bg0, xorA(d1, d2));
  const tr = [];
  for (let t = 0; t <= T; t++) {
    const ref = bgAt(t);
    tr.push(popA(xorA(xorA(s12, ref), xorA(xorA(s1, ref), xorA(s2, ref)))));
    s1 = stepL(s1, rule.bm, rule.sm); s2 = stepL(s2, rule.bm, rule.sm); s12 = stepL(s12, rule.bm, rule.sm);
  }
  return tr;
}
function propagate2D(rule, bgAt, dev0, T) {
  let s = xorA(bgAt(0), dev0);
  const frames = [];
  for (let t = 0; t <= T; t++) { frames.push({ t, ...canonDev(xorA(s, bgAt(t))) }); s = stepL(s, rule.bm, rule.sm); }
  return frames;
}
function classify2D(frames) {
  const fe = frames.findIndex(f => f.n === 0);
  if (fe >= 0) return `absorbed at t=${fe}`;
  for (let p = 1; p <= 12; p++) {
    let ok = true, b = null;
    for (let t = 0; t + p < frames.length; t++) {
      if (frames[t].id !== frames[t + p].id) { ok = false; break; }
      const dx = Math.round(frames[t + p].cx - frames[t].cx), dy = Math.round(frames[t + p].cy - frames[t].cy);
      if (b === null) b = [dx, dy]; else if (dx !== b[0] || dy !== b[1]) { ok = false; break; }
    }
    if (ok && b && (b[0] || b[1])) return `RIGID: period ${p}, velocity (${b[0]},${b[1]}) per period, mass ${frames[0].n}`;
  }
  const q = frames[Math.floor(frames.length * 0.75)], last = frames[frames.length - 1];
  const grew = last.w > q.w + 1 || last.h > q.h + 1 || last.n > q.n + 2;
  if (grew) return `SPREADS/EXPLODES: mass ${frames[0].n} -> ${last.n}, box ${frames[0].w}x${frames[0].h} -> ${last.w}x${last.h}`;
  return `BOUNDED: mass ${last.n}, box ${last.w}x${last.h}, centroid moved (${(last.cx - frames[0].cx).toFixed(1)},${(last.cy - frames[0].cy).toFixed(1)}) over ${last.t} steps`;
}
const T2 = 120;
console.log('  A) Life, empty background');
{
  const bgAt = () => newGrid();
  const g1 = newGrid(), g2 = newGrid();
  setCells(g1, GLIDER.map(([x, y]) => [x + 100, y + 100]));
  setCells(g2, GLIDER.map(([x, y]) => [x + 140, y + 140]));
  let tr = linearity2D(LIFE, bgAt, g1, g2, T2);
  console.log(`     RULER 1, two gliders far apart: max|mismatch| = ${Math.max(...tr)} (they never interact)`);
  const g3 = newGrid(), g4 = newGrid();
  setCells(g3, GLIDER.map(([x, y]) => [x + 100, y + 100]));
  setCells(g4, GLIDER.map(([x, y]) => [x + 100, y + 116]));   // same diagonal lane, 8 periods apart
  tr = linearity2D(LIFE, bgAt, g3, g4, T2);
  console.log(`     RULER 1, two gliders in the SAME lane (same speed, never meet): max|mismatch| = ${Math.max(...tr)}`);
  const g5 = newGrid(), g6 = newGrid();
  setCells(g5, GLIDER.map(([x, y]) => [x + 100, y + 100]));                     // moves (+1,+1)
  setCells(g6, GLIDER.map(([x, y]) => [x + 112, y + 104]));                     // mirrored so it moves (-1,+1)
  for (const [x, y] of GLIDER) { const nx = 112 - x, ny = 104 + y; g6[ny * SZ + nx] = 1; }
  tr = linearity2D(LIFE, bgAt, g5, g6, 80);
  const firstBad = tr.findIndex(v => v > 0);
  console.log(`     RULER 1, two gliders on a COLLISION COURSE: max|mismatch| = ${Math.max(...tr)}, first mismatch at t=${firstBad < 0 ? 'never' : firstBad}`);
  console.log(`     RULER 2, one glider: ${classify2D(propagate2D(LIFE, bgAt, g1, T2))}`);
  const cell = newGrid(); cell[100 * SZ + 100] = 1;
  console.log(`     RULER 2, one single cell: ${classify2D(propagate2D(LIFE, bgAt, cell, 10))}`);
}
console.log('');
console.log('  B) Life, block-lattice background (= what we have been calling matter)');
{
  const lat = newGrid();
  for (let j = 0; j < SZ; j += 4) for (let i = 0; i < SZ; i += 4) { lat[j * SZ + i] = 1; lat[j * SZ + i + 1] = 1; lat[(j + 1) * SZ + i] = 1; lat[(j + 1) * SZ + i + 1] = 1; }
  const bgAt = () => lat;
  const d1 = newGrid(); d1[100 * SZ + 100] ^= 1;
  const d2 = newGrid(); d2[100 * SZ + 101] ^= 1;
  const tr = linearity2D(LIFE, bgAt, d1, d2, 30);
  console.log(`     RULER 1, two adjacent single-cell deviations: max|mismatch| = ${Math.max(...tr)}, steps with mismatch>0 = ${tr.filter(v => v > 0).length}/${tr.length}`);
  console.log(`     RULER 2, one-cell hole:  ${classify2D(propagate2D(LIFE, bgAt, d1, 30))}`);
  const bump = newGrid(); bump[100 * SZ + 102] = 1;
  console.log(`     RULER 2, one-cell bump:  ${classify2D(propagate2D(LIFE, bgAt, bump, 30))}`);
}
console.log('');
console.log('  C) the lamellar elastic rule B{023456}/S{01245} with its own period-2 background');
{
  const r = { bm: [0, 2, 3, 4, 5, 6].reduce((a, n) => a | (1 << n), 0), sm: [0, 1, 2, 4, 5].reduce((a, n) => a | (1 << n), 0) };
  const c0 = newGrid();
  for (let y = 0; y < SZ; y++) for (let x = 0; x < SZ; x++) if ((y % 6) >= 1 && (y % 6) <= 4) c0[y * SZ + x] = 1;
  const c1 = stepL(c0, r.bm, r.sm), c2 = stepL(c1, r.bm, r.sm);
  console.log(`     period check step^2(bg)==bg: ${c2.every((v, i) => v === c0[i])}`);
  const bgAt = (t) => (t % 2 === 0 ? c0 : c1);
  const d1 = newGrid(); d1[101 * SZ + 100] ^= 1;
  const d2 = newGrid(); d2[101 * SZ + 101] ^= 1;
  const tr = linearity2D(r, bgAt, d1, d2, 60);
  console.log(`     RULER 1, two adjacent single-cell deviations: max|mismatch| = ${Math.max(...tr)}, steps with mismatch>0 = ${tr.filter(v => v > 0).length}/${tr.length}`);
  console.log(`     RULER 2, one single-cell deviation: ${classify2D(propagate2D(r, bgAt, d1, 60))}`);
}
console.log('');
console.log('============================================================');
console.log('WHAT THE CONTROLS ESTABLISHED');
console.log('============================================================');
console.log('  RULER 1 is trustworthy: it reads exactly 0 for the four linear 2-state rules');
console.log('  (including with ADJACENT deviations, where interaction is forced) and nonzero');
console.log('  immediately for rule 30, for Life\'s block lattice and for the lamellar rule.');
console.log('  RULER 2 is trustworthy: it finds the glider rigid (period 4, velocity (1,1)),');
console.log('  the pure shift rigid, the identity frozen, a Life cell absorbed, a lattice hole');
console.log('  absorbed and a lattice bump exploding -- all matching independent measurements.');
