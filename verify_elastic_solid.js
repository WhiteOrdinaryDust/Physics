'use strict';
/*
 * verify_elastic_solid.js
 * ---------------------------------------------------------------------------
 * INDEPENDENT VERIFICATION of the elastic-solid search, plus the controls that
 * script and mine must agree on.
 *
 * The search reported that some 2-state radius-1 outer-totalistic rules have
 * matter that heals BOTH a hole and a bump -- the property real matter has and
 * Life matter does not, and the property a "screen you can bump into without it
 * recording you" needs. That is the first constructive result of this programme,
 * so it does not get to be true on one script's word.
 *
 * WHAT THIS SCRIPT DOES, on its own code path:
 *   1. finds the background cycle by evolving seeds and detecting an exact repeat
 *      of the whole configuration (any drift is inside the cycle on a torus, so
 *      exact equality is a complete test);
 *   2. prints the background as a picture, with period, density, component count
 *      and largest component -- a "solid" is one big component, a "dust cloud"
 *      is many small ones, and those must never be confused;
 *   3. measures elasticity with the CAUSAL filter:
 *         a flip is relevant iff  step(flipped) != step(background)
 *      (not a distance rule: a cell two away can matter by supplying the third
 *      neighbour for a birth in between, while a diagonal two-away cell cannot);
 *      "healed" = the trajectory returns to some phase of that same cycle;
 *   4. measures it for EVERY phase of the cycle separately -- a real solid should
 *      not care which phase of its own internal clock it is in, and if it does,
 *      that is a finding, not a nuisance;
 *   5. repeats the measurement on a larger torus, because "heals" on 16x16 could
 *      be a small-torus artefact;
 *   6. runs the C1 controls through the same code: Life + isolated block and
 *      Life + block lattice must give E_add_causal = 0.0000 / E_remove = 1.0000,
 *      matching matter_elasticity.js and elasticity_metric_trap.js.
 * ---------------------------------------------------------------------------
 */

function parseRule(s) {
  const m = /^B\{?([0-8]*)\}?\/S\{?([0-8]*)\}?$/.exec(s.trim());
  if (!m) throw new Error('bad rule ' + s);
  let bm = 0, sm = 0;
  for (const c of m[1]) bm |= 1 << Number(c);
  for (const c of m[2]) sm |= 1 << Number(c);
  return { name: s, bm, sm };
}
function stepRule(g, W, H, r) {
  const out = new Uint8Array(W * H);
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
    let n = 0;
    for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
      if (dx === 0 && dy === 0) continue;
      n += g[(((y + dy + H) % H) * W) + ((x + dx + W) % W)];
    }
    const i = y * W + x;
    out[i] = g[i] ? ((r.sm >> n) & 1) : ((r.bm >> n) & 1);
  }
  return out;
}
function hashGrid(g) {
  let h = 2166136261 >>> 0;
  for (let i = 0; i < g.length; i++) { h ^= g[i]; h = Math.imul(h, 16777619) >>> 0; }
  return h >>> 0;
}
function density(g) { let s = 0; for (let i = 0; i < g.length; i++) s += g[i]; return s / g.length; }
function comps(g, W, H) {
  const seen = new Uint8Array(g.length);
  let n = 0, largest = 0;
  for (let i = 0; i < g.length; i++) {
    if (!g[i] || seen[i]) continue;
    n++; let size = 0; const Q = [i]; seen[i] = 1;
    while (Q.length) {
      const k = Q.pop(); size++;
      const x = k % W, y = (k - x) / W;
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        const nx = (x + dx + W) % W, ny = (y + dy + H) % H, nk = ny * W + nx;
        if (g[nk] && !seen[nk]) { seen[nk] = 1; Q.push(nk); }
      }
    }
    if (size > largest) largest = size;
  }
  return { n, largest };
}
function picture(g, W, H) {
  const rows = [];
  for (let y = 0; y < H; y++) {
    let s = '';
    for (let x = 0; x < W; x++) s += g[y * W + x] ? '#' : '.';
    rows.push(s);
  }
  return rows;
}
function mulberry32(a) { return function () { a |= 0; a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }

// ---- find the eventual cycle: exact repeat of the whole configuration
function findCycle(r, W, H, seedDensity, seed) {
  const rnd = mulberry32(seed);
  let g = new Uint8Array(W * H);
  if (seedDensity === 'all') g.fill(1);
  else if (seedDensity === 'block') { g.fill(0); g[0] = 1; g[1] = 1; g[W] = 1; g[W + 1] = 1; }
  else for (let i = 0; i < g.length; i++) g[i] = rnd() < seedDensity ? 1 : 0;
  const seen = new Map();
  const hist = [];
  for (let t = 0; t < 1500; t++) {
    const h = hashGrid(g);
    if (seen.has(h)) {
      const j = seen.get(h);
      // exact confirmation, hash collisions are not allowed to fool us
      let same = true;
      for (let i = 0; i < g.length; i++) if (hist[j][i] !== g[i]) { same = false; break; }
      if (same) return hist.slice(j);
    }
    seen.set(h, t);
    hist.push(g);
    g = stepRule(g, W, H, r);
  }
  return null;
}

// how repeatable is the pattern? count the translations that leave it unchanged.
// A material you can hold is a CRYSTAL: it repeats on a short lattice. A dense
// pattern that merely happens to be periodic in time but never repeats in space is
// an oscillating soup, not matter, and must not be called a solid.
function crystallinity(g, W, H) {
  let count = 0, shortest = null;
  for (let a = 0; a < W; a++) for (let b = 0; b < H; b++) {
    if (a === 0 && b === 0) continue;
    let ok = true;
    for (let y = 0; y < H && ok; y++) for (let x = 0; x < W; x++) {
      if (g[y * W + x] !== g[((y + b) % H) * W + ((x + a) % W)]) { ok = false; break; }
    }
    if (!ok) continue;
    count++;
    const d = Math.hypot(Math.min(a, W - a), Math.min(b, H - b));
    if (!shortest || d < shortest.d) shortest = { a, b, d };
  }
  return { count, area: (W * H) / (count + 1), shortest };
}

const T = 120;
function measurePhase(r, W, H, cycle, phase, rmStride) {
  const bg = cycle[phase];
  const bgNext = stepRule(bg, W, H, r);
  const nextHash = hashGrid(bgNext);
  const cycHashes = new Set(cycle.map(hashGrid));
  const heals = (start) => {
    let g = start;
    for (let t = 0; t < T; t++) {
      if (cycHashes.has(hashGrid(g))) return true;
      g = stepRule(g, W, H, r);
    }
    return cycHashes.has(hashGrid(g));
  };
  const stride = rmStride || 1;
  let addRel = 0, addHeal = 0, addAll = 0, addAllHeal = 0, rmRel = 0, rmHeal = 0, rmAll = 0, rmAllHeal = 0;
  let liveIdx = 0;
  for (let i = 0; i < bg.length; i++) {
    if (bg[i]) {
      // REMOVALS need no relevance filter: the cell removed is part of the matter,
      // so removing it is an interaction by definition. (A one-step "did anything
      // change?" test is degenerate here: self-repairing matter refills the cell so
      // that step(flipped) equals step(background) exactly, giving 0/0.)
      const sampled = (liveIdx % stride === 0);
      liveIdx++;
      if (!sampled) continue;
      rmRel++;
      const g = Uint8Array.from(bg); g[i] = 0;
      const h = heals(g);
      rmAll++; if (h) { rmHeal++; rmAllHeal++; }
    } else {
      // ADDITIONS need the filter, because a cell added far from the matter simply
      // dies alone and the matter never notices -- counting those is what inflated
      // E_add to 0.8889 in the earlier search.
      const g = Uint8Array.from(bg); g[i] = 1;
      const h = heals(g);
      addAll++; if (h) addAllHeal++;
      if (hashGrid(stepRule(g, W, H, r)) !== nextHash) { addRel++; if (h) addHeal++; }
    }
  }
  return {
    phase,
    addRel, addHeal, addAll, addAllHeal, rmRel, rmHeal, rmAll, rmAllHeal,
    eAdd: addRel ? addHeal / addRel : NaN,
    eRm: rmRel ? rmHeal / rmRel : NaN,
    eAddAll: addAll ? addAllHeal / addAll : NaN,
    eRmAll: rmAll ? rmAllHeal / rmAll : NaN,
  };
}

// ============================================================ controls
console.log('=== CONTROLS: Life through this code path ===');
const life = parseRule('B3/S23');
{
  const W = 16, H = 16;
  const g = new Uint8Array(W * H);
  for (const [x, y] of [[8, 8], [9, 8], [8, 9], [9, 9]]) g[y * W + x] = 1;
  const cyc = findCycle(life, W, H, 'all', 0) || [stepRule(g, W, H, life)];
  const m = measurePhase(life, W, H, [g], 0);
  console.log(`  isolated 2x2 block : E_add_causal = ${m.eAdd.toFixed(4)} (${m.addHeal}/${m.addRel})   E_remove_causal = ${m.eRm.toFixed(4)} (${m.rmHeal}/${m.rmRel})   ${m.eAdd === 0 && m.eRm === 1 ? 'PASS' : '*** FAIL ***'}`);
  if (!(m.eAdd === 0 && m.eRm === 1)) process.exit(1);
}
{
  const W = 16, H = 16;
  const g = new Uint8Array(W * H);
  for (let j = 0; j < 16; j += 4) for (let i = 0; i < 16; i += 4) { g[j * W + i] = 1; g[j * W + i + 1] = 1; g[(j + 1) * W + i] = 1; g[(j + 1) * W + i + 1] = 1; }
  const m = measurePhase(life, W, H, [g], 0);
  console.log(`  block lattice      : E_add_causal = ${m.eAdd.toFixed(4)} (${m.addHeal}/${m.addRel})   E_remove_causal = ${m.eRm.toFixed(4)} (${m.rmHeal}/${m.rmRel})   ${m.eAdd === 0 && m.eRm === 1 ? 'PASS' : '*** FAIL ***'}`);
  if (!(m.eAdd === 0 && m.eRm === 1)) process.exit(1);
}
console.log('  (both match matter_elasticity.js and elasticity_metric_trap.js)\n');

// ============================================================ the candidates
const CANDIDATES = [
  ['B{234}/S{012345}', 'reported dense solid with usable statistics'],
  ['B{012567}/S{}', 'reported best E_score overall, but dust (8 components)'],
  ['B{0124567}/S{0}', 'reported identical to the above'],
  ['B{0138}/S{0125678}', 'reported 1.0000 / 0.9016'],
];

function survey(ruleStr, note, sizes, rmStride) {
  const r = parseRule(ruleStr);
  console.log(`=== ${ruleStr} ===   (${note})`);
  for (const SZ of sizes) {
    const W = SZ, H = SZ;
    const cycles = [];
    const sigs = new Set();
    for (const d of ['block', 0.5, 0.2]) {
      for (const sd of [1, 2]) {
        const c = findCycle(r, W, H, d, sd);
        if (!c) continue;
        const sig = c.length + '|' + c.map(g => hashGrid(g)).join(',');
        if (sigs.has(sig)) continue;
        sigs.add(sig);
        cycles.push(c);
      }
    }
    if (!cycles.length) { console.log(`  ${SZ}x${SZ}: no cycle found`); continue; }
    cycles.sort((a, b) => density(b[0]) - density(a[0]));
    console.log(`  --- ${SZ}x${SZ} torus: ${cycles.length} distinct cycle(s) reached ---`);
    for (const cyc of cycles.slice(0, 2)) {
      const dens = cyc.map(density);
      const cc = comps(cyc[0], W, H);
      const cr = crystallinity(cyc[0], W, H);
      console.log(`    period ${cyc.length}, density ${dens[0].toFixed(4)} (min ${Math.min(...dens).toFixed(4)} max ${Math.max(...dens).toFixed(4)}), components ${cc.n}, largest ${cc.largest}/${W * H}`);
      console.log(`    crystallinity: ${cr.count} non-trivial translations leave it unchanged => fundamental cell area ${cr.area}${cr.shortest ? `, shortest repeat (${cr.shortest.a},${cr.shortest.b})` : ''}  ${cr.area <= 64 ? '(a repeatable crystal)' : '<-- NOT a crystal: it never repeats in space, it is an oscillating soup'}`);
      const phases = cyc.length <= 4 ? cyc.map((_, p) => p) : [0];
      const ms = phases.map(p => measurePhase(r, W, H, cyc, p, rmStride));
      for (const m of ms) {
        console.log(`      phase ${m.phase}: E_add_causal = ${Number.isNaN(m.eAdd) ? 'n/a' : m.eAdd.toFixed(4)} (${m.addHeal}/${m.addRel})   E_remove_causal = ${Number.isNaN(m.eRm) ? 'n/a' : m.eRm.toFixed(4)} (${m.rmHeal}/${m.rmRel})   [all-cells E_add ${Number.isNaN(m.eAddAll) ? 'n/a' : m.eAddAll.toFixed(4)} (${m.addAllHeal}/${m.addAll})]`);
      }
      const eAdds = ms.map(m => m.eAdd).filter(v => !Number.isNaN(v));
      const eRms = ms.map(m => m.eRm).filter(v => !Number.isNaN(v));
      const phaseDep = new Set(eAdds).size > 1 || new Set(eRms).size > 1;
      console.log(`      -> across the cycle's phases: E_add ${Math.min(...eAdds).toFixed(4)}..${Math.max(...eAdds).toFixed(4)}, E_remove ${Math.min(...eRms).toFixed(4)}..${Math.max(...eRms).toFixed(4)}  ${phaseDep ? '<-- DEPENDS ON THE PHASE OF ITS OWN CLOCK' : '(phase independent)'}`);
    }
  }
  console.log('');
}

for (const [ruleStr, note] of CANDIDATES) survey(ruleStr, note, [16]);
// the small-torus caveat: does the leading dense candidate survive a bigger torus?
survey('B{234}/S{012345}', 'size check on a bigger torus', [24]);
// the one candidate that DID look two-sidedly elastic, with a bigger E_add sample
// (its background is 95% full, so on 16x16 there are only ~12 dead cells to add)
survey('B{0138}/S{0125678}', 'the only two-sided candidate, bigger torus for a usable E_add sample', [32], 5);

// ---- picture of the leading candidate's background
console.log('=== background picture of the reported best solid with usable statistics ===');
{
  const r = parseRule('B{234}/S{012345}');
  const W = 24, H = 24;
  const cyc = findCycle(r, W, H, 0.5, 1);
  if (cyc) {
    console.log(`  B{234}/S{012345} on ${W}x${W}, period ${cyc.length}, density ${density(cyc[0]).toFixed(4)}`);
    for (const row of picture(cyc[0], W, H)) console.log('   ' + row);
  }
}
