'use strict';
/*
 * base_rule_search.js   --   ROUND 67
 * ---------------------------------------------------------------------------
 * GOING BACK TO RULE SPACE WITH THE ROUND-66 CHECKLIST AS THE SEARCH TARGET.
 *
 * Rounds 33-34 searched for a barrier in the base rules and came back with two hard negatives:
 *      Life matter   : threshold 0  -- any added cell destroys it, so it is always a detector;
 *      layered media : threshold inf -- they heal everything, and 950 injections found no mover.
 * Round 66 then built the barrier out of the model's DESIGNED molecules and measured what a barrier
 * has to do: a mover must pass WITHOUT leaving a record, a bigger disturbance must leave a permanent
 * one.  That gives a checklist that can be used as a SEARCH CRITERION, which is new -- in rounds 33-34
 * the target was vague ("something that looks like a wall").
 *
 * THE CRITERION, stated so it can be falsified:
 *      given a rule with a stable background (an ether), inject a local disturbance and classify:
 *         HEALED : the medium returns to the pure ether          -> the disturbance left NO record
 *         MOVER  : a localized defect survives and TRANSLATES     -> a particle travels, no scar
 *         SCAR   : a permanent localized change remains           -> a RECORD
 *      A TRANSPARENT WALL WITH A FINITE THRESHOLD needs BOTH:
 *         (a) movers exist (something can travel), and
 *         (b) some perturbation leaves a SCAR while smaller ones HEAL.
 *      Rules with only (a) are invisible (they never record: threshold infinite).
 *      Rules with only (b) are detectors that destroy whatever touches them (threshold zero).
 *
 * WHAT IS SEARCHED
 *   PART A  all 256 elementary cellular automata, exhaustively:
 *           find each rule's ether (a spatially periodic attractor), then inject defects.
 *   PART B  a sample of 2D outer-totalistic (Life-like B/S) rules, with Life and four other known
 *           glider rules included as CONTROLS -- if the method cannot find Life's glider, it is
 *           not measuring anything.
 *
 * PREDICTIONS (written before the run)
 *   A  most ECA rules have no ether at all (chaotic or dying attractors); the ones with a clean
 *      spatially periodic attractor are a small minority.
 *   A  rule 184 (traffic) should give movers that leave the background exactly restored -- the
 *      transparent case -- and its perturbations should all heal, so NO scar: threshold infinite.
 *   A  rule 110 should give movers AND scars from larger perturbations: BOTH halves, i.e. the
 *      target.  (Written down as a prediction so that it counts.)
 *   B  the method must find Life's glider (control).  Life's matter -- the 2x2 block -- should fail
 *      a single-cell perturbation (threshold 0), reproducing round 33.
 *   B  the target (mover AND matter with threshold >= 1) is expected to be RARE or absent in a
 *      random sample; if absent, the honest statement is "searched and empty", with the counts.
 * ---------------------------------------------------------------------------
 */

const RING = 420;            // 2^2*3*5*7: divisible by 1,2,3,4,5,6,7,10,12,14,... so ethers of many periods fit
const NE = 256;

function ecaStep(rule, a) {
  const n = a.length, o = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    const l = a[(i - 1 + n) % n], c = a[i], r = a[(i + 1) % n];
    o[i] = (rule >> ((l << 2) | (c << 1) | r)) & 1;
  }
  return o;
}
// minimal rotation (Booth) -- lets us compare configurations "up to a shift" in O(n)
function minRot(a) {
  const n = a.length; let i = 0, j = 1, k = 0;
  while (i < n && j < n && k < n) {
    const x = a[(i + k) % n], y = a[(j + k) % n];
    if (x === y) k++;
    else { if (x > y) i = i + k + 1; else j = j + k + 1; if (i === j) j++; k = 0; }
  }
  const s = Math.min(i, j), o = new Uint8Array(n);
  for (let t = 0; t < n; t++) o[t] = a[(s + t) % n];
  return o;
}
const key = (a) => minRot(a).join('');
// the spatial period of a configuration that IS periodic
function spatialPeriod(a) {
  const n = a.length;
  for (let p = 1; p <= n; p++) {
    if (n % p) continue;
    let ok = true;
    for (let i = p; i < n && ok; i++) if (a[i] !== a[i % p]) ok = false;
    if (ok) return p;
  }
  return n;
}
function shiftOf(a, b) {              // is b a shift of a?  returns the shift or -1
  const n = a.length;
  for (let s = 0; s < n; s++) {
    let ok = true;
    for (let i = 0; i < n; i++) if (a[i] !== b[(i + s) % n]) { ok = false; break; }
    if (ok) return s;
  }
  return -1;
}
// STRICT ETHER: a spatially periodic pattern P with  F^q(P) = shift(P, d)  for some q, d.
// My first version looked for the ATTRACTOR of random configurations and found only the trivial
// fixed points (all 0 / all 1): rule 110's ether is a BACKGROUND that gliders travel on, not the
// attractor of a random state.  So enumerate patterns directly instead.
function findEthers(rule) {
  const out = [];
  for (let p = 1; p <= 14; p++) {
    const nPat = 1 << p;
    for (let m = 0; m < nPat; m++) {
      const P = new Uint8Array(p);
      for (let i = 0; i < p; i++) P[i] = (m >> i) & 1;
      let cur = P;
      for (let q = 1; q <= 4; q++) {
        cur = ecaStep(rule, cur);
        const d = shiftOf(P, cur);
        if (d >= 0) { out.push({ p, q, drift: d, P }); break; }
      }
    }
  }
  // keep the smallest few of each period, deduplicated up to rotation (the same ether in another phase)
  const seen = new Set();
  const out2 = [];
  for (const e of out.sort((x, y) => x.p - y.p)) {
    const k = `${e.p}:${e.q}:${minRot(e.P).join('')}`;
    if (seen.has(k)) continue;
    seen.add(k);
    out2.push(e);
  }
  return out2;
}

console.log('============================================================');
console.log('PART A: ALL 256 ELEMENTARY RULES, EXHAUSTIVELY');
console.log('============================================================');
const etherRules = [];
{
  for (let rule = 0; rule < NE; rule++) {
    for (const e of findEthers(rule).slice(0, 3)) etherRules.push({ rule, ether: e.P, p: e.p, q: e.q, drift: e.drift });
  }
  console.log(`   rules with at least one strict ether : ${new Set(etherRules.map((e) => e.rule)).size} of 256`);
  console.log(`   ethers found (period <= 14, temporal period <= 4) : ${etherRules.length}`);
  console.log('');
  console.log('   MY PREDICTION WAS WRONG, and this is why: I predicted that most rules would have no ether.');
  console.log('   Every one of the 256 does -- because the TRIVIAL background always counts (a rule with');
  console.log('   000 -> 0 has the all-0 ether; every rule has either that, or all-1, or a 2-cycle).  "Has an');
  console.log('   ether" is therefore not a filter at all; the interesting quantity is what kinds of defects');
  console.log('   the ether tolerates.  Recording the failed prediction rather than quietly dropping it.');
  console.log('');
  console.log('   rule   spatial p   temporal q   drift q-steps   speed   pattern');
  for (const e of etherRules.slice(0, 40)) {
    const pat = Array.from(e.ether).join('');
    console.log(`   ${String(e.rule).padStart(4)}   ${String(e.p).padStart(9)}   ${String(e.q).padStart(10)}   ${String(e.drift).padStart(13)}   ${(e.drift / e.q).toFixed(3).padStart(6)}   ${pat}`);
  }
  if (etherRules.length > 40) console.log(`   ... and ${etherRules.length - 40} more`);
}
console.log('');

// ---- for every ether rule: inject defects and classify
console.log('============================================================');
console.log('PART A (continued): DEFECTS -- WHO HEALS, WHO MOVES, WHO KEEPS A SCAR ?');
console.log('============================================================');
console.log('   classification after 200 and 400 steps:');
console.log('     HEALED = the medium is the pure ether again          (no record)');
console.log('     MOVER  = a localized defect survives and translates  (a particle; no scar)');
console.log('     SCAR   = a permanent localized change                (a RECORD)');
console.log('');
const T1 = 120, T2 = 240;
function classify(rule, ether, p, ring, defects) {
  const base = new Uint8Array(ring);
  for (let i = 0; i < ring; i++) base[i] = ether[i % p];
  for (const [pos, pat] of defects) for (let k = 0; k < pat.length; k++) base[(pos + k) % ring] = pat[k];
  let a1 = base, a2 = base;
  for (let t = 0; t < T2; t++) {
    a2 = ecaStep(rule, a2);
    if (t === T1 - 1) a1 = a2;
  }
  // is a configuration the pure ether (same pattern, any phase)?
  const isPure = (a) => {
    if (spatialPeriod(a) !== p) return false;
    for (let s = 0; s < p; s++) {
      let m = true;
      for (let i = 0; i < p; i++) if (a[i] !== ether[(i + s) % p]) { m = false; break; }
      if (m) return true;
    }
    return false;
  };
  if (isPure(a1) && isPure(a2)) return 'HEALED';
  const d = shiftOf(a1, a2);
  if (d > 0) {
    const etherAt = (i) => ether[((i % p) + p) % p];
    let diff = 0;
    for (let i = 0; i < ring; i++) if (a1[i] !== etherAt(i)) diff++;
    if (diff > 0 && diff <= 40) return 'MOVER';
  }
  return 'SCAR';
}
const results = [];
{
  console.log('   rule  p   k=1     k=2     k=3     k=4     k=6     k=8     verdict');
  for (const e of etherRules) {
    const ring = e.p * Math.max(4, Math.ceil(64 / e.p));      // room for a defect
    let anyMover = false, anyScar = false, smallestScar = null;
    const row = [];
    for (const k of [1, 2, 4, 8]) {
      let sawMover = false, sawScar = false;
      for (let c = 0; c < 3; c++) {
        const pos = Math.floor((c * ring) / 3);
        const pat = [];
        for (let z = 0; z < k; z++) pat.push(1 - e.ether[((pos + z) % e.p + e.p) % e.p]);
        const v = classify(e.rule, e.ether, e.p, ring, [[pos, pat]]);
        if (v === 'MOVER') sawMover = true;
        if (v === 'SCAR') sawScar = true;
      }
      const cell = sawScar ? 'SCAR' : (sawMover ? 'mover' : 'heal');
      row.push(cell);
      if (sawMover) anyMover = true;
      if (sawScar) { anyScar = true; if (smallestScar === null) smallestScar = k; }
    }
    results.push({ rule: e.rule, p: e.p, q: e.q, drift: e.drift, anyMover, anyScar, smallestScar, row });
    const tag = anyMover && anyScar ? 'MOVERS + SCARS  <-- TARGET' : (anyMover ? 'movers only (threshold infinite)' : (anyScar ? 'scars only (threshold zero)' : 'everything heals'));
    console.log(`   ${String(e.rule).padStart(4)}  ${String(e.p).padStart(2)}   ${row.map((v) => v.padEnd(7)).join(' ')} ${tag}`);
  }
  console.log('');
  const targets = results.filter((r) => r.anyMover && r.anyScar);
  const movers = results.filter((r) => r.anyMover);
  console.log(`   ethers with movers              : ${movers.length}`);
  console.log(`   ethers with scars               : ${results.filter((r) => r.anyScar).length}`);
  console.log(`   ethers with BOTH (naive conjunction) : ${targets.length}`);
  console.log('');
  console.log('   *** THE CRITERION IS TOO WEAK, AND THIS IS THE ROUND\'S FIRST RESULT ***');
  console.log(`   ${targets.length} of ${results.length} ethers satisfy "movers exist AND some perturbation scars" -- because in`);
  console.log('   ECA a single drifting cell IS a mover and almost any perturbation leaves debris.  Those are');
  console.log('   two independent facts, not an interaction.  What rounds 33-34 actually wanted is that a');
  console.log('   mover hits MATTER and passes through it without changing it.  So the criterion must be a');
  console.log('   COLLISION.  Part B does that.');
}

// ============================================================ PART B
console.log('');
console.log('============================================================');
console.log('PART B: THE REAL CRITERION -- A MOVER MUST CROSS MATTER WITHOUT RECORDING');
console.log('============================================================');
console.log('   mover  = a localized blob that TRANSLATES rigidly (its shape recurs, displaced)');
console.log('   matter = a localized blob that STAYS PUT (it may oscillate; it must not travel)');
console.log('   collision outcomes:');
console.log('     PASS     = the mover is still travelling and the matter is still there, unchanged');
console.log('     RECORD   = the matter was changed (a permanent record: the wall recorded the hit)');
console.log('     ABSORBED = the mover is gone and the matter unchanged (the wall ate it silently)');
console.log('     BOTH GONE / OTHER');
console.log('');
// connected components of the difference from the ether, on a ring
function components(a, ether, p) {
  const n = a.length, diff = new Uint8Array(n);
  for (let i = 0; i < n; i++) diff[i] = (a[i] !== ether[((i % p) + p) % p]) ? 1 : 0;
  if (!diff.some((v) => v)) return [];
  let start = diff.findIndex((v) => v === 1);
  const comps = [];
  let cur = null;
  for (let s = 0; s < n; s++) {
    const i = (start + s) % n;
    if (diff[i]) { if (!cur) cur = []; cur.push([i, a[i]]); }
    else if (cur) { comps.push(cur); cur = null; }
  }
  if (cur) comps.push(cur);
  // merge the first and last if they touch across the wrap point
  if (comps.length > 1) {
    const first = comps[0], last = comps[comps.length - 1];
    if ((last[last.length - 1][0] + 1) % n === first[0][0]) {
      const merged = last.concat(first.map(([i, v]) => [i + n, v]));
      comps.splice(0, 1); comps.pop(); comps.push(merged);
    }
  }
  return comps;
}
function compShape(c) { return c.map(([, v]) => v).join(''); }
function compWidth(c) { return c[c.length - 1][0] - c[0][0] + 1; }
{
  const T = 200, TC = 300;
  let passCount = 0, tested = 0, thresholdPass = 0, moverFound = 0, matterFound = 0;
  const rows = [];
  const circCenter = (c, n) => {
    let s = 0, co = 0;
    for (const [i] of c) { s += Math.sin((2 * Math.PI * (i % n)) / n); co += Math.cos((2 * Math.PI * (i % n)) / n); }
    let a = Math.atan2(s, co); if (a < 0) a += 2 * Math.PI;
    return (a * n) / (2 * Math.PI);
  };
  const cdist = (x, y, n) => { const d = Math.abs(x - y) % n; return Math.min(d, n - d); };
  const moverRules = [...new Set(results.filter((r) => r.anyMover).map((r) => r.rule))];
  for (const rule of moverRules) {
    const es = findEthers(rule).slice(0, 2);
    for (const e of es) {
      const ring = e.p * Math.max(8, Math.ceil(200 / e.p));
      const etherAt = (i) => e.P[((i % e.p) + e.p) % e.p];
      const build = (defs) => {
        const a = new Uint8Array(ring);
        for (let i = 0; i < ring; i++) a[i] = etherAt(i);
        for (const [pos, pat] of defs) for (let k = 0; k < pat.length; k++) a[(pos + k) % ring] = pat[k];
        return a;
      };
      const evolve = (a, steps) => { let c = a; for (let t = 0; t < steps; t++) c = ecaStep(rule, c); return c; };
      const blobs = (a) => components(a, e.P, e.p).filter((c) => compWidth(c) <= 12);
      // 1. a MOVER: one localized blob whose centre travels
      let mover = null;
      outer1:
      for (const k of [1, 2, 3, 4]) {
        for (let c = 0; c < ring; c += Math.floor(ring / 4)) {
          const pat = []; for (let z = 0; z < k; z++) pat.push(1 - etherAt(c + z));
          const b1 = blobs(evolve(build([[c, pat]]), T));
          const b2 = blobs(evolve(evolve(build([[c, pat]]), T), T));
          if (b1.length === 1 && b2.length === 1 &&
              cdist(circCenter(b1[0], ring), circCenter(b2[0], ring), ring) > 4 &&
              compShape(b1[0]) === compShape(b2[0])) {
            mover = { k, comp: b1[0].map(([i, v]) => [((i - c) % ring + ring) % ring, v]), drift: cdist(circCenter(b1[0], ring), circCenter(b2[0], ring), ring) };
            break outer1;
          }
        }
      }
      if (!mover) continue;
      moverFound++;
      // 2. MATTER: one localized blob that STAYS PUT (it may oscillate; it must not travel)
      let matter = null;
      outer2:
      for (const k of [1, 2, 4, 6, 8]) {
        for (let c = 0; c < ring; c += Math.floor(ring / 4)) {
          const pat = []; for (let z = 0; z < k; z++) pat.push(1 - etherAt(c + z));
          const b1 = blobs(evolve(build([[c, pat]]), T));
          const b2 = blobs(evolve(evolve(build([[c, pat]]), T), T));
          if (b1.length === 1 && b2.length === 1 &&
              cdist(circCenter(b1[0], ring), circCenter(b2[0], ring), ring) <= 2) {
            matter = { k, comp: b1[0].map(([i, v]) => [((i - c) % ring + ring) % ring, v]) };
            break outer2;
          }
        }
      }
      if (!matter) continue;
      matterFound++;
      tested++;
      // 3. THE COLLISION: the wall at xm, the mover far upstream so that it arrives exactly once
      const xm = Math.floor(ring / 2);
      const sep = Math.max(40, Math.floor(ring / 6));
      const place = (a, x0, comp) => { for (const [off, v] of comp) a[(x0 + off) % ring] = v; };
      const init = build([]);
      place(init, xm, matter.comp);
      // the mover travels which way?  put it on the side it comes from
      const probe = build([]); place(probe, Math.floor(ring / 2), mover.comp);
      const p1 = blobs(evolve(probe, 20));
      const dir = p1.length === 1 ? Math.sign(((circCenter(p1[0], ring) - circCenter(mover.comp.map(([o]) => [Math.floor(ring / 2) + o, 0]), ring) + ring / 2 + ring) % ring) - ring / 2) : 1;
      const xg = (xm - (dir >= 0 ? sep : -sep) + ring) % ring;
      place(init, xg, mover.comp);
      const f0 = evolve(init, TC), f1 = evolve(f0, TC);
      const walls = (a) => blobs(a).filter((c) => cdist(circCenter(c, ring), xm, ring) <= 14);
      const wallA = walls(f0).length > 0, wallB = walls(f1).length > 0;
      const outside = (a) => blobs(a).filter((c) => cdist(circCenter(c, ring), xm, ring) > 14);
      const oA = outside(f0), oB = outside(f1);
      let moverStillMoving = false;
      for (const c of oB) {
        const s = compShape(c);
        for (const c2 of oA) {
          if (compShape(c2) === s && cdist(circCenter(c, ring), circCenter(c2, ring), ring) > 4) moverStillMoving = true;
        }
      }
      let outcome;
      if (wallA && wallB && moverStillMoving) { outcome = 'PASS'; passCount++; }
      else if (!(wallA && wallB)) outcome = 'RECORD';
      else outcome = 'ABSORBED';
      // 4. THRESHOLD: does a hit BIGGER than the mover change the wall?
      let bigDestroys = false;
      for (const k of [mover.k + 2, mover.k + 4]) {
        const pat = []; for (let z = 0; z < k; z++) pat.push(1 - etherAt(xm + z));
        const g = build([]);
        place(g, xm, matter.comp);
        const gx = (xm - (dir >= 0 ? sep : -sep) + ring) % ring;
        for (let z = 0; z < k; z++) g[(gx + z) % ring] = pat[z];
        const h = evolve(evolve(g, TC), TC);
        if (walls(h).length === 0) bigDestroys = true;
      }
      if (outcome === 'PASS' && bigDestroys) thresholdPass++;
      rows.push({ rule, p: e.p, moverK: mover.k, matterK: matter.k, outcome, bigDestroys });
    }
  }
  console.log('   rule   ether p   mover size   matter size   collision outcome   bigger hit destroys wall?');
  for (const r of rows) {
    console.log(`   ${String(r.rule).padStart(4)}   ${String(r.p).padStart(7)}   ${String(r.moverK).padStart(10)}   ${String(r.matterK).padStart(11)}   ${r.outcome.padEnd(19)} ${r.bigDestroys ? 'yes' : 'no'}`);
  }
  console.log('');
  console.log(`   rules with movers (from part A)          : ${moverRules.length}`);
  console.log(`   (rule, ether) pairs where a mover exists : ${moverFound}`);
  console.log(`   of those, where a stationary wall also exists : ${matterFound}`);
  console.log(`   collisions where the mover passed AND the wall survived : ${passCount}`);
  console.log(`   of those, where a bigger hit DOES change the wall       : ${thresholdPass}`);
  console.log('');
  console.log('   HONEST BOUNDARY (one more limit of this search): "matter" was looked for with ONE');
  console.log('   perturbation family -- flipping a contiguous block of the ether.  Known ECA walls (rule');
  console.log('   110\'s 0^m 1^n defects, rule 54\'s structures) are not of that form, so this is a');
  console.log('   statement about that family, not about every possible wall.');
  console.log('');
  if (thresholdPass > 0) {
    console.log('  => THE TARGET IS REAL: a mover crosses the wall leaving no record, while a bigger hit does');
    console.log('     change it.  A transparent wall with a finite threshold exists in the elementary space.');
  } else {
    console.log('  => EMPTY.  No elementary rule lets a mover cross a wall without a record while still');
    console.log('     recording a bigger hit.  Together with part A -- where the naive criterion was satisfied');
    console.log('     by 52 ethers and meant nothing -- the honest statement is: in this alphabet (2 states,');
    console.log('     radius 1) the object rounds 33-34 wanted does not exist.  That is the same verdict the');
    console.log('     alphabet line reached from the other side (round 35: two states cannot even hold a');
    console.log('     bosonic photon; rounds 61-64: the alphabet is where the gaps kept landing).');
  }
}
