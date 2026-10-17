'use strict';
/*
 * eca_ether_sweep.js
 * ---------------------------------------------------------------------------
 * Sweep all 256 elementary cellular automata for PROPAGATING periodic backgrounds.
 *
 * WHY
 *   The two substrates tried so far both failed to provide a "slow medium":
 *     Life       : no periodic background at all (ash is inert -- measured)
 *     rule 110   : has 133 distinct periodic backgrounds, but every one has
 *                  shift = 0 in canonical form, i.e. NONE of them propagates
 *   A two-layer structure needs a background whose phase ADVANCES in space
 *   (shift != 0 over its temporal period T), so that a slower clock exists.
 *
 * CRITERION
 *   A rule qualifies if there exists a periodic seed (period L <= LMAX) such that
 *        pattern(t + T) == pattern(t) rotated by k,  with k != 0
 *   and the measured phase velocity |k|/T does not exceed the light-speed limit
 *   (for a radius-1 rule, at most 1 cell per step).
 *
 * REUSE
 *   The self-checking machinery from ether_search_v2/ether_canonical:
 *     - L passed explicitly, defined once per call
 *     - canonical form so that T is the true period (no aliasing)
 *     - pattern must stay L-periodic throughout, else the seed is rejected
 *     - negative control: a random seed must NOT qualify
 *     - speed-limit check as an explicit assertion
 * ---------------------------------------------------------------------------
 */

function makeRule(rule) {
  const t = new Uint8Array(8);
  for (let n = 0; n < 8; n++) t[n] = (rule >> n) & 1;
  return t;
}
function step(row, w, t) {
  const out = new Uint8Array(w);
  for (let i = 0; i < w; i++) {
    out[i] = t[(row[(i - 1 + w) % w] << 2) | (row[i] << 1) | row[(i + 1) % w]];
  }
  return out;
}
function canonStr(bits) {
  let best = null;
  for (let k = 0; k < bits.length; k++) {
    const r = bits.slice(k).concat(bits.slice(0, k)).join('');
    if (best === null || r < best) best = r;
  }
  return best;
}

/*
 * Analyse one periodic seed under one rule.
 * Returns null unless the seed is a background with a WELL-DEFINED phase advance.
 */
function analyseSeed(bits, L, t, TMAX, MAXT) {
  const w = L * 6;
  if (w < L + 4) return null;
  let row = new Uint8Array(w);
  for (let i = 0; i < w; i++) row[i] = bits[i % L];
  const snaps = [];
  for (let a = 0; a <= TMAX; a++) {
    snaps.push(Uint8Array.from(row));
    row = step(row, w, t);
  }
  // the state must be a pure L-periodic pattern at every time
  const canonAt = [];
  for (let a = 0; a <= TMAX; a++) {
    const win = [];
    for (let i = 0; i < L; i++) win.push(snaps[a][i]);
    let periodic = true;
    for (let i = 0; i < 4 * L && i < w; i++) if (snaps[a][i] !== win[i % L]) { periodic = false; break; }
    if (!periodic) return null;
    canonAt.push(canonStr(win));
  }
  // true temporal period of the canonical form
  for (let T = 1; T <= MAXT; T++) {
    let ok = true;
    for (let a = 0; a + T <= TMAX; a++) if (canonAt[a] !== canonAt[a + T]) { ok = false; break; }
    if (!ok) continue;
    const A = canonAt[0], B = canonAt[T];
    let k = -1;
    for (let s = 0; s < L; s++) if (B.slice(s).concat(B.slice(0, s)) === A) { k = s; break; }
    if (k < 0) continue;
    // a displacement of k and of k-L are the same physical advance
    let kk = k;
    if (kk > L / 2) kk -= L;
    const v = kk / T;
    return { T, k: kk, v };
  }
  return null;
}

// ============================================================ negative control
console.log('=== NEGATIVE CONTROL (must fail to qualify) ===');
{
  const t = makeRule(110);
  const L = 13, bits = new Array(L).fill(0); bits[0] = 1;
  const r = analyseSeed(bits, L, t, 60, 20);
  console.log(`  rule 110, seed = single 1 : ${r ? '*** qualified -> SEARCH BROKEN ***' : 'did not qualify (correct)'}`);
  if (r) process.exit(1);
}
console.log('');

// ============================================================ the sweep
const LMAX = 12;
const TMAX = 60;
const MAXT = 30;

console.log(`=== sweeping all 256 elementary CA (L <= ${LMAX}, TMAX=${TMAX}) ===\n`);

const qualifying = [];
const t0 = Date.now();

for (let rule = 0; rule < 256; rule++) {
  const t = makeRule(rule);
  let best = null;
  let count = 0;
  for (let L = 2; L <= LMAX; L++) {
    const maxMask = 1 << L;
    for (let mask = 1; mask < maxMask; mask++) {
      const bits = [];
      for (let i = 0; i < L; i++) bits.push((mask >> i) & 1);
      const r = analyseSeed(bits, L, t, TMAX, MAXT);
      if (!r) continue;
      if (Math.abs(r.v) > 1 + 1e-9) continue;        // speed-limit assertion
      if (r.k === 0) continue;                        // NOT propagating: reject here
      count++;
      const ones = bits.reduce((a, b) => a + b, 0);
      if (!best || ones < best.ones || (ones === best.ones && Math.abs(r.v) > Math.abs(best.v))) {
        best = { L, bits: bits.join(''), ones, ...r };
      }
    }
  }
  if (best) {
    qualifying.push({ rule, count, ...best });
    console.log(`  rule ${String(rule).padStart(3)}  propagating backgrounds: ${String(count).padStart(5)}   sparsest: L=${best.L} ${best.bits} (ones=${best.ones}) T=${best.T} k=${best.k} v=${best.v.toFixed(4)}`);
  }
}

console.log(`\nsweep done in ${((Date.now() - t0) / 1000).toFixed(1)}s`);
console.log(`rules with propagating periodic backgrounds: ${qualifying.length} of 256`);

// ---- summary table
if (qualifying.length) {
  console.log('\n=== all qualifying rules, with their fastest propagating background ===');
  console.log('  rule   count   best v     L   ones   T   k');
  const sorted = [...qualifying].sort((a, b) => Math.abs(b.v) - Math.abs(a.v));
  for (const q of sorted) {
    console.log(`  ${String(q.rule).padStart(4)}   ${String(q.count).padStart(5)}   ${q.v.toFixed(4).padStart(7)}   ${String(q.L).padStart(2)}   ${String(q.ones).padStart(4)}   ${String(q.T).padStart(2)}   ${String(q.k).padStart(3)}`);
  }
  const nonZero = qualifying.filter(q => q.v !== 0);
  console.log(`\nof these, ${nonZero.length} have a strictly non-zero phase velocity (the others`);
  console.log('qualified through seeds whose representative shift was 0 after sign folding).');
  const fastest = sorted.filter(q => q.v !== 0)[0];
  if (fastest) {
    const ratio = 1 / Math.abs(fastest.v);
    console.log(`\nfastest propagating background: rule ${fastest.rule}, v = ${fastest.v.toFixed(4)} cells/step`);
    console.log(`  => c0 / c1 = ${ratio.toFixed(2)}  (substrate max speed 1 / phase speed |v|)`);
    console.log('  A ratio well above 1 is the "slow medium" a two-layer structure needs.');
  }
} else {
  console.log('\nNO elementary CA in L<=12 has a propagating periodic background.');
}
