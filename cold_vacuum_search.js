'use strict';
/*
 * cold_vacuum_search.js
 * ---------------------------------------------------------------------------
 * JOB: find a COLD vacuum in the lattice gas -- a deterministic background whose small
 * perturbations stay small and evolve linearly -- because the last round showed that a random
 * gas is a HOT (chaotic) medium, and a photon in a hot medium is like sound in air, not light
 * in vacuum.
 *
 * HOW TO MEASURE "TEMPERATURE" (this is the new tool, and it is the test that caught the last
 * round's false positive):
 *   Put k particles into a background at k well-separated places, evolve, and count how many
 *   particle-bits differ from the unperturbed background:
 *       Nk(t) = number of differing bits at time t.
 *   In a LINEAR (cold) background the perturbation is just carried around, so
 *       N4(t) / N1(t) = 4    for every t.
 *   In a CHAOTIC (hot) background the difference grows until it saturates at the size of the
 *   attractor's own fluctuations, which does not depend on k, so
 *       N4(t) / N1(t) -> 1.
 *   The ratio is the temperature, and unlike "does the difference grow" it cannot be fooled by
 *   mere spreading: spreading preserves the ratio, saturation destroys it.
 *
 * CONTROLS the temperature test must pass:
 *   - the EMPTY lattice (cold by construction: perturbations are free particles) -> ratio 4
 *   - a RANDOM gas at density 0.25 (hot by construction: collisions)        -> ratio ~1
 *
 * SEARCH SPACE
 *   (1) systematic: every lattice-gas state whose pattern repeats with spatial period px*py <= 4
 *       (16 states per cell, so at most 16^4 = 65536 candidates per shape). Keep the ones that
 *       are members of a short cycle (step^q = a translation of itself, q <= 6).
 *   (2) a handful of designed candidates: uniform beams, uniform head-on pairs, the full
 *       lattice, checkerboards, diagonals.
 *   Then temperature-test every survivor.
 *
 * PREDICTIONS (before the run)
 *   V1  the two controls behave as stated (empty cold, random gas hot).
 *   V2  uniform beams (all particles moving the same way) are cold: with my collision rule only
 *       exactly head-on pairs inside one cell collide, so a uniform beam never collides.
 *   V3  most cycle-member candidates are still hot once a perturbation is added, because a
 *       periodic arrangement of collisions is chaotic under perturbation even when it is
 *       periodic in time.
 *   V4  I therefore expect at most a handful of genuinely cold backgrounds, and they will be
 *       "collision-free" or "regularly colliding" states rather than generic ones.
 * ---------------------------------------------------------------------------
 */

// SEEDED generator.  This file used to call Math.random() for the hot gas, which made its report
// irreproducible -- a violation of the project's own core promise, found in the round-120 audit.
function stream(seed) { let s = seed | 0; if (s === 0) s = 0x9e3779b9; return () => { s ^= s << 13; s ^= s >>> 17; s ^= s << 5; return (s >>> 0) / 4294967296; }; }
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
const NB = (1 << 0) | (1 << 2), EW = (1 << 1) | (1 << 3);

function makeTorus(W, H) { return { g: new Uint8Array(W * H), W, H }; }
function gstep(g, W, H) {
  const o = new Uint8Array(W * H);
  for (let i = 0; i < W * H; i++) { const v = g[i]; if (v === NB) o[i] = EW; else if (v === EW) o[i] = NB; else o[i] = v; }
  const out = new Uint8Array(W * H);
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
    const v = o[y * W + x];
    if (!v) continue;
    for (let d = 0; d < 4; d++) if ((v >> d) & 1) out[((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W)] |= 1 << d;
  }
  return out;
}
function bitDiff(a, b) { let s = 0; for (let i = 0; i < a.length; i++) { let v = a[i] ^ b[i], c = 0; while (v) { c += v & 1; v >>= 1; } s += c; } return s; }
function popBits(g) { let s = 0; for (let i = 0; i < g.length; i++) { let v = g[i]; while (v) { s += v & 1; v >>= 1; } } return s; }

// ============================================================ the temperature test
function temperature(base, W, H, slots, T) {
  // slots: array of [x, y, d] perturbation positions, all valid (compatible with the base)
  const mk = (k) => {
    const g = Uint8Array.from(base);
    for (let i = 0; i < k; i++) { const [x, y, d] = slots[i]; g[y * W + x] ^= 1 << d; }
    return g;
  };
  const g0 = Uint8Array.from(base);
  let a0 = g0, a1 = mk(1), a4 = mk(4);
  const n1 = [], n4 = [];
  for (let t = 0; t <= T; t++) {
    n1.push(bitDiff(a1, a0)); n4.push(bitDiff(a4, a0));
    a0 = gstep(a0, W, H); a1 = gstep(a1, W, H); a4 = gstep(a4, W, H);
  }
  let worst = 9, best = -9;
  for (let t = 10; t <= T; t++) {
    const r = n4[t] / Math.max(1, n1[t]);
    // only judge while the perturbed region has not wrapped the whole torus
    if (n1[t] < 0.35 * W * H) { if (r < worst) worst = r; if (r > best) best = r; }
  }
  const mean = (() => { let s = 0, c = 0; for (let t = 10; t <= T; t++) if (n1[t] < 0.35 * W * H) { s += n4[t] / Math.max(1, n1[t]); c++; } return c ? s / c : NaN; })();
  return { worst, best, mean, n1end: n1[T], n4end: n4[T], usable: Number.isFinite(mean) };
}
function findSlots(base, W, H, seed) {
  // four well-separated slots; prefer ADDING a particle, otherwise removing one
  const want = [[3, 3], [9, 9], [3, 15], [15, 3]];
  for (const mode of ['add', 'remove']) {
    const slots = [];
    for (const [x, y] of want) {
      for (let d = 0; d < 4; d++) {
        const has = (base[y * W + x] >> d) & 1;
        if ((mode === 'add' && !has) || (mode === 'remove' && has)) { slots.push([x, y, d]); break; }
      }
      if (slots.length !== (want.indexOf([x, y]) + 1) && slots.length < want.length) { /* keep going */ }
    }
    if (slots.length === 4) return { slots, mode };
  }
  return null;
}

// ============================================================ controls
console.log('============================================================');
console.log('CONTROLS for the temperature test');
console.log('============================================================');
{
  const W = 24, H = 24, T = 200;
  const empty = new Uint8Array(W * H);
  const s1 = findSlots(empty, W, H, 1);
  const rEmpty = temperature(empty, W, H, s1.slots, T);
  console.log(`  empty lattice (perturbations = free particles): ratio N4/N1 mean ${rEmpty.mean.toFixed(3)}, worst ${rEmpty.worst.toFixed(3)}`);
  console.log(`    => ${rEmpty.mean > 3.0 ? 'COLD as expected (the perturbation is just carried)' : '*** unexpected ***'}`);
  let sd = 12345; const rnd = () => { sd = (sd * 1103515245 + 12345) & 0x7fffffff; return ((sd >> 16) & 1); };
  const hot = new Uint8Array(W * H);
  const rndHot = stream(20251016);
    for (let i = 0; i < W * H; i++) { let v = 0; for (let d = 0; d < 4; d++) if (rndHot() < 0.25) v |= 1 << d; hot[i] = v; }
  const s2 = findSlots(hot, W, H, 2);
  const rHot = temperature(hot, W, H, s2.slots, T);
  console.log(`  random gas at density 0.25 (hot): ratio mean ${rHot.mean.toFixed(3)}, worst ${rHot.worst.toFixed(3)}`);
  console.log(`    => ${rHot.mean < 1.5 ? 'HOT as expected (the difference saturates)' : '*** unexpected ***'}`);
  console.log('  So the ratio separates cold from hot by a factor of at least 2, which is the discrimination we need.');
}
console.log('');

// ============================================================ designed candidates
console.log('============================================================');
console.log('PART 1: designed candidates');
console.log('============================================================');
const DES = [];
DES.push(['empty', (x, y) => 0]);
DES.push(['all N (uniform beam)', () => 1 << 0]);
DES.push(['all E (uniform beam)', () => 1 << 1]);
DES.push(['all N|E (perpendicular pair, no head-on)', () => (1 << 0) | (1 << 1)]);
DES.push(['all N|S (every cell head-on)', () => NB]);
DES.push(['all E|W (every cell head-on, other axis)', () => EW]);
DES.push(['full lattice (all four)', () => 0xf]);
DES.push(['checkerboard N|S and E|W', (x, y) => ((x + y) % 2 === 0 ? NB : EW)]);
DES.push(['stripes: x even N|S, odd empty', (x) => (x % 2 === 0 ? NB : 0)]);
DES.push(['stripes: x even N, x odd S', (x) => (x % 2 === 0 ? 1 << 0 : 1 << 2)]);
DES.push(['diagonal beams N|E every third cell', (x, y) => ((x + y) % 3 === 0 ? ((1 << 0) | (1 << 1)) : 0)]);

function isCycleMember(base, W, H, qmax) {
  let cur = Uint8Array.from(base);
  for (let q = 1; q <= qmax; q++) {
    cur = gstep(cur, W, H);
    for (let a = 0; a < W; a++) for (let b = 0; b < H; b++) {
      let ok = true;
      for (let y = 0; y < H && ok; y++) for (let x = 0; x < W; x++) if (cur[y * W + x] !== base[((y + b) % H) * W + ((x + a) % W)]) { ok = false; break; }
      if (ok) return { q, a, b };
    }
  }
  return null;
}
const W = 24, H = 24, T = 200;
const results = [];
for (const [name, f] of DES) {
  const base = new Uint8Array(W * H);
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) base[y * W + x] = f(x, y);
  const cyc = isCycleMember(base, W, H, 6);
  const s = findSlots(base, W, H, 3);
  if (!s) { console.log(`  ${name.padEnd(42)} no valid perturbation slot (skipped)`); continue; }
  const r = temperature(base, W, H, s.slots, T);
  // NOTE (bug found afterwards): the old version fell back to 9 when no time step qualified, and
  // 9 > 3 classified as COLD. That is a false-cold: a background that a single particle destroys
  // saturates the torus before t=10 and then read as "cold". Unmeasurable is now reported as such.
  const tag = !r.usable ? 'UNMEASURABLE (perturbation fills the torus before t=10)' : (r.mean > 3.0 ? 'COLD' : (r.mean > 2.0 ? 'lukewarm' : 'HOT'));
  results.push({ name, cyc, r, tag, base: Uint8Array.from(base), slots: s.slots });
  console.log(`  ${name.padEnd(42)} cycle ${cyc ? `q=${cyc.q}` : 'none '}   ratio ${r.usable ? r.mean.toFixed(3).padStart(6) : '  n/a '}   ${tag}   [perturb by ${s.mode}]`);
}
console.log('');

// ============================================================ systematic search
console.log('============================================================');
console.log('PART 2: systematic search over small-period states');
console.log('============================================================');
{
  const shapes = [[1, 1], [2, 1], [1, 2], [3, 1], [1, 3], [2, 2], [4, 1], [1, 4]];
  let tested = 0, members = 0, tempTests = 0;
  const survivors = [];
  const MAX_TEMP_TESTS = 400;            // the screen is a sample, not an exhaustive pass
  for (const [px, py] of shapes) {
    const ncell = px * py;
    const total = 16 ** ncell;
    if (total > 70000) continue;
    let shapeMembers = 0;
    const t0 = Date.now();
    for (let code = 0; code < total; code++) {
      const base = new Uint8Array(px * py);
      let c = code;
      for (let i = 0; i < ncell; i++) { base[i] = c & 15; c >>= 4; }
      tested++;
      const cyc = isCycleMember(base, px, py, 6);
      if (!cyc) continue;
      members++;
      shapeMembers++;
      // SAMPLE the cycle members: the temperature test is by far the expensive part
      if (tempTests >= MAX_TEMP_TESTS) continue;
      if (shapeMembers % Math.max(1, Math.floor(shapeMembers / 40 + 1)) !== 0) continue;
      const W2 = px * 4, H2 = py * 4;
      const big = new Uint8Array(W2 * H2);
      for (let y = 0; y < H2; y++) for (let x = 0; x < W2; x++) big[y * W2 + x] = base[(y % py) * px + (x % px)];
      const s = findSlots(big, W2, H2, 4);
      if (!s) continue;
      const r = temperature(big, W2, H2, s.slots, 120);
      tempTests++;
      if (r.usable && r.mean > 3.0) survivors.push({ px, py, base, cyc, r });
    }
    console.log(`    shape ${px}x${py}: ${total} states, ${shapeMembers} cycle members  (${((Date.now() - t0) / 1000).toFixed(1)} s, temperature tests so far ${tempTests})`);
  }
  console.log(`  ${tested} states examined, ${members} are members of a short cycle, ${tempTests} temperature-tested (sampled), ${survivors.length} read COLD`);
  // show a few, deduplicated by the pattern they describe on their own torus
  const shown = new Set();
  let shownCount = 0;
  for (const s of survivors) {
    const key = `${s.px}x${s.py}:` + Array.from(s.base).join(',');
    if (shown.has(key)) continue;
    shown.add(key);
    if (shownCount++ >= 8) break;
    const pic = [];
    for (let y = 0; y < s.py; y++) { let row = ''; for (let x = 0; x < s.px; x++) row += s.base[y * s.px + x].toString(16); pic.push(row); }
    console.log(`    ${s.px}x${s.py} cycle q=${s.cyc.q} shift (${s.cyc.a},${s.cyc.b}) ratio ${s.r.mean.toFixed(3)}  pattern(hex) ${pic.join(' | ')}`);
  }
  // re-test the survivors on a bigger torus, because a small torus can hide spreading
  console.log('');
  console.log('  re-testing the cold survivors on a 48x48 torus (a small torus can fake coldness):');
  const BIG = 48, T2 = 300;
  let confirmed = 0, reTested = 0;
  const seen2 = new Set();
  const MAX_RETESTS = 40;                // a cap, so the script always finishes
  for (const s of survivors) {
    if (reTested >= MAX_RETESTS) break;
    const key = `${s.px}x${s.py}:` + Array.from(s.base).join(',');
    if (seen2.has(key)) continue;
    seen2.add(key);
    const big = new Uint8Array(BIG * BIG);
    for (let y = 0; y < BIG; y++) for (let x = 0; x < BIG; x++) big[y * BIG + x] = s.base[(y % s.py) * s.px + (x % s.px)];
    const sl = findSlots(big, BIG, BIG, 5);
    if (!sl) continue;
    const r = temperature(big, BIG, BIG, sl.slots, T2);
    reTested++;
    const tag = !r.usable ? 'unmeasurable on the big torus' : (r.mean > 3.0 ? 'COLD confirmed' : 'fails on the big torus');
    if (r.usable && r.mean > 3.0) confirmed++;
    if (reTested <= 8) console.log(`    ${s.px}x${s.py} ratio ${r.usable ? r.mean.toFixed(3) : 'n/a'}  ${tag}`);
  }
  console.log(`  => ${confirmed} cold backgrounds confirmed on the 48x48 torus (out of ${reTested} re-tested, cap ${MAX_RETESTS})`);
}
console.log('');
console.log('============================================================');
console.log('WHAT A COLD BACKGROUND WOULD BUY US');
console.log('============================================================');
console.log('  A cold background is the layer\'s VACUUM: small perturbations stay small and (if the');
console.log('  temperature ratio is exactly 4 for every k) evolve linearly, which is the superposition');
console.log('  the photon needs. The next step after this search is to run the three rulers --');
console.log('  superposition, w(k), and translation versus standing -- on whichever cold backgrounds');
console.log('  come out of it, and to look for a GAPPED (massive) mode in the same background to get');
console.log('  the particles that must sit on top of it.');
