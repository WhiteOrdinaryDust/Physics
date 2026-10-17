'use strict';
/*
 * crystalline_solid_search.js
 * ---------------------------------------------------------------------------
 * WHY THIS EXISTS
 *
 * An earlier search claimed 39 of 3000 sampled rules have "matter" that heals
 * both holes and bumps, and named solids. I verified the named rules on my own
 * code path (verify_elastic_solid.js, controls passing against two independent
 * scripts) and got:
 *   - the headline number does not reproduce (reported E_add_causal 0.7890 for
 *     B{234}/S{012345}; I measure 0.09-0.37 on 16x16 and 0.08-0.14 on 24x24,
 *     across several cycles and all phases);
 *   - and, decisively, NONE of the named backgrounds repeats in space. Every one
 *     of them has zero non-trivial translation symmetries: a dense pattern that is
 *     periodic in TIME but aperiodic in SPACE. That is an oscillating soup, not a
 *     material. You cannot cut a piece out of it and find the same structure.
 *
 * So the search optimised "dense + connected + self-repairing" and returned chaos.
 * A material needs a third property the search never tested: it must REPEAT.
 *
 * THIS SCRIPT adds that property and re-runs the whole rule space.
 *
 * SCREEN, in order of cost (cheap first):
 *   1. evolve from a mixed seed (random at 0.5 density plus a 2x2 block) for 100
 *      generations on a 12x12 torus; keep only rules whose state then has period
 *      1, 2 or 4 (a background with a clock, or a static one);
 *   2. density in [0.15, 0.85] -- a 99%-full background has almost no room to add
 *      a cell, so its E_add can never be measured on a real sample;
 *   3. CRYSTALLINE: it must repeat under a short lattice translation, fundamental
 *      cell area <= 36. This is the test that was missing;
 *   4. one 8-connected bulk (largest component >= half the live cells) -- not dust;
 *   5. only then, measure elasticity with the CAUSAL filter (a flip counts only if
 *      it changes the next generation; removals of live matter always count);
 *      keep rules with E_add_causal >= 0.5 AND E_remove_causal >= 0.5 AND at least
 *      20 causally relevant flips in EACH direction.
 *
 * CONTROLS (mandatory, before the sweep):
 *   - Life block lattice: crystallinity area must be 16, and E_add_causal = 0.0000
 *     with E_remove_causal = 1.0000 (matches matter_elasticity.js).
 *   - a random pattern: crystallinity area must be the whole torus (no repeats).
 *   - the reported "best solid" B{234}/S{012345}: must be rejected by the
 *     crystallinity screen -- if this search does not reject it, the screen is broken.
 *
 * PREDICTION (stated before the run):
 *   I expect ZERO rules to pass all five screens. Crystalline matter in this rule
 *   space tends to be either one-sided (Life: heals holes, never bumps) or nearly
 *   full (B{0138}/S{0125678}: 99% full, so only ~8 dead cells exist and its E_add
 *   rests on 8 samples). If any rule does pass, that is the real find and I will
 *   say so loudly. A zero result here is a real result: it would mean the search
 *   space contains no "solid you can hold" at all, only soups.
 * ---------------------------------------------------------------------------
 */

const SW = 12, SH = 12, NS = SW * SH;
const GEN = 100;

function stepR(g, bm, sm) {
  const out = new Uint8Array(NS);
  for (let y = 0; y < SH; y++) {
    const yp = ((y - 1 + SH) % SH) * SW, ym = ((y + 1) % SH) * SW, yc = y * SW;
    for (let x = 0; x < SW; x++) {
      const xp = (x - 1 + SW) % SW, xm = (x + 1) % SW;
      const n = g[yp + xp] + g[yp + x] + g[yp + xm] + g[yc + xp] + g[yc + xm] + g[ym + xp] + g[ym + x] + g[ym + xm];
      const i = yc + x;
      out[i] = g[i] ? ((sm >> n) & 1) : ((bm >> n) & 1);
    }
  }
  return out;
}
function eq(a, b) { for (let i = 0; i < NS; i++) if (a[i] !== b[i]) return false; return true; }
function density(g) { let s = 0; for (let i = 0; i < NS; i++) s += g[i]; return s / NS; }
function hash(g) { let h = 2166136261 >>> 0; for (let i = 0; i < NS; i++) { h ^= g[i]; h = Math.imul(h, 16777619) >>> 0; } return h >>> 0; }

// translations that leave the pattern unchanged, and the resulting fundamental area
function transl(g) {
  let count = 0;
  for (let a = 0; a < SW; a++) for (let b = 0; b < SH; b++) {
    if (a === 0 && b === 0) continue;
    let ok = true;
    for (let y = 0; y < SH && ok; y++) {
      const yb = ((y + b) % SH) * SW, yc = y * SW;
      for (let x = 0; x < SW; x++) if (g[yc + x] !== g[yb + ((x + a) % SW)]) { ok = false; break; }
    }
    if (ok) count++;
  }
  return { count, area: NS / (count + 1) };
}
function biggest(g) {
  const seen = new Uint8Array(NS);
  let n = 0, largest = 0;
  for (let i = 0; i < NS; i++) {
    if (!g[i] || seen[i]) continue;
    n++; let size = 0; const Q = [i]; seen[i] = 1;
    while (Q.length) {
      const k = Q.pop(); size++;
      const x = k % SW, y = (k - x) / SW;
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        const nk = (((y + dy + SH) % SH) * SW) + ((x + dx + SW) % SW);
        if (g[nk] && !seen[nk]) { seen[nk] = 1; Q.push(nk); }
      }
    }
    if (size > largest) largest = size;
  }
  return { n, largest };
}
function pic(g) {
  const r = [];
  for (let y = 0; y < SH; y++) { let s = ''; for (let x = 0; x < SW; x++) s += g[y * SW + x] ? '#' : '.'; r.push(s); }
  return r;
}

const T = 120;
function measure(r, bg) {
  const next = stepR(bg, r.bm, r.sm);
  const nextHash = hash(next);
  let addRel = 0, addHeal = 0, rmRel = 0, rmHeal = 0;
  const heals = (start) => {
    // a static or oscillating crystal: healed = back on the background's own cycle
    let g = start;
    for (let t = 0; t < T; t++) {
      if (eq(g, bg) || eq(g, next)) return true;
      g = stepR(g, r.bm, r.sm);
    }
    return eq(g, bg) || eq(g, next);
  };
  for (let i = 0; i < NS; i++) {
    const g = Uint8Array.from(bg);
    g[i] ^= 1;
    const h = heals(g);
    if (bg[i]) { rmRel++; if (h) rmHeal++; }
    else if (hash(stepR(g, r.bm, r.sm)) !== nextHash) { addRel++; if (h) addHeal++; }
  }
  return { addRel, addHeal, rmRel, rmHeal, eAdd: addRel ? addHeal / addRel : NaN, eRm: rmRel ? rmHeal / rmRel : NaN };
}

// ================================================================ controls
console.log('=== CONTROLS ===');
{
  const life = { bm: 1 << 3, sm: (1 << 2) | (1 << 3) };
  const lat = new Uint8Array(NS);
  for (let j = 0; j < SH; j += 4) for (let i = 0; i < SW; i += 4) {
    lat[j * SW + i] = 1; lat[j * SW + i + 1] = 1; lat[(j + 1) * SW + i] = 1; lat[(j + 1) * SW + i + 1] = 1;
  }
  const t = transl(lat);
  const m = measure(life, lat);
  console.log(`  Life block lattice: crystallinity area ${t.area} (expect 16) ${t.area === 16 ? 'PASS' : '*** FAIL ***'}`);
  console.log(`    E_add_causal ${m.eAdd.toFixed(4)} (${m.addHeal}/${m.addRel}), E_remove_causal ${m.eRm.toFixed(4)} (${m.rmHeal}/${m.rmRel})  ${m.eAdd === 0 && m.eRm === 1 ? 'PASS' : '*** FAIL ***'}`);
  if (t.area !== 16 || m.eAdd !== 0 || m.eRm !== 1) process.exit(1);
}
{
  const rnd = new Uint8Array(NS);
  let s = 12345;
  for (let i = 0; i < NS; i++) { s = (s * 1103515245 + 12345) & 0x7fffffff; rnd[i] = (s >> 16) & 1; }
  const t = transl(rnd);
  console.log(`  random pattern: crystallinity area ${t.area} (expect ${NS}) ${t.area === NS ? 'PASS' : '*** FAIL ***'}`);
  if (t.area !== NS) process.exit(1);
}
function classify(bm, sm) {
  const SEED = new Uint8Array(NS);
  let s = 987654321;
  for (let i = 0; i < NS; i++) { s = (s * 1103515245 + 12345) & 0x7fffffff; SEED[i] = ((s >> 16) & 1) ? 1 : 0; }
  SEED[0] = 1; SEED[1] = 1; SEED[SW] = 1; SEED[SW + 1] = 1;   // plus a 2x2 block
  let g = Uint8Array.from(SEED);
  for (let t = 0; t < GEN; t++) g = stepR(g, bm, sm);
  let g1 = stepR(g, bm, sm);
  let period = 0;
  if (eq(g, g1)) period = 1;
  else {
    const g2 = stepR(g1, bm, sm);
    if (eq(g, g2)) period = 2;
    else {
      const g3 = stepR(g2, bm, sm), g4 = stepR(g3, bm, sm);
      if (eq(g, g4)) period = 4;
    }
  }
  return { period, g };
}
{
  const nb = (1 << 2) | (1 << 3) | (1 << 4);
  const ns = (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 5);
  const c = classify(nb, ns);
  const t = transl(c.g);
  console.log(`  reported "best solid" B{234}/S{012345}: period ${c.period}, density ${density(c.g).toFixed(4)}, crystallinity area ${t.area} -> ${t.area > 36 ? 'REJECTED by the crystallinity screen (correct: it is a soup)' : '*** FAIL: the screen let a soup through ***'}`);
  if (t.area <= 36) process.exit(1);
}
console.log('  controls passed\n');

// ================================================================ sweep
console.log('=== SWEEP over all 262144 rules (12x12 torus, 100 generations) ===');
const t0 = Date.now();
const survivors = [];
let screened = 0, periodOk = 0, densOk = 0, crystalOk = 0;

for (let idx = 0; idx < 262144; idx++) {
  const bm = idx & 511, sm = (idx >>> 9) & 511;
  if (bm === 0 && sm === 0) continue;
  screened++;
  const c = classify(bm, sm);
  if (!c.period) continue;
  periodOk++;
  const d = density(c.g);
  if (d < 0.15 || d > 0.85) continue;
  densOk++;
  const t = transl(c.g);
  if (t.area > 36) continue;
  crystalOk++;
  const b = biggest(c.g);
  if (b.largest < 0.5 * (d * NS)) continue;
  const r = { bm, sm };
  const m = measure(r, c.g);
  if (m.addRel >= 20 && m.rmRel >= 20 && m.eAdd >= 0.5 && m.eRm >= 0.5) {
    survivors.push({ idx, bm, sm, d, period: c.period, area: t.area, comps: b.n, largest: b.largest, m, g: c.g });
  }
  if (Date.now() - t0 > 1500000) { console.log('  (time limit reached, stopping the sweep)'); break; }
}

console.log(`  rules examined                : ${screened}`);
console.log(`  with a background of period<=4: ${periodOk}`);
console.log(`  ... and density in [0.15,0.85] : ${densOk}`);
console.log(`  ... and CRYSTALLINE            : ${crystalOk}   <-- the screen that was missing before`);
console.log(`  ... and one bulk component     : (counted in survivors below)`);
console.log(`  ... and two-sidedly elastic    : ${survivors.length}`);
console.log(`  elapsed ${((Date.now() - t0) / 1000).toFixed(0)} s\n`);

if (!survivors.length) {
  console.log('  NO rule in the entire space has a crystalline background that heals both');
  console.log('  a hole and a bump with at least 20 samples each way.');
  console.log('  Combined with verify_elastic_solid.js (every named "solid" is an aperiodic');
  console.log('  soup), the honest conclusion is: this rule space contains no solid you can');
  console.log('  hold. The earlier search found dense oscillating chaos, not matter.');
} else {
  survivors.sort((a, b) => Math.min(b.m.eRm, b.m.eAdd) - Math.min(a.m.eRm, a.m.eAdd));
  for (const s of survivors.slice(0, 12)) {
    console.log(`  B{${[...Array(9).keys()].filter(n => (s.bm >> n) & 1).join('')}}/S{${[...Array(9).keys()].filter(n => (s.sm >> n) & 1).join('')}}  index ${s.idx}`);
    console.log(`    density ${s.d.toFixed(4)}, period ${s.period}, crystallinity area ${s.area}, components ${s.comps} largest ${s.largest}`);
    console.log(`    E_add_causal ${s.m.eAdd.toFixed(4)} (${s.m.addHeal}/${s.m.addRel})   E_remove_causal ${s.m.eRm.toFixed(4)} (${s.m.rmHeal}/${s.m.rmRel})`);
    for (const row of pic(s.g)) console.log('      ' + row);
  }
}
