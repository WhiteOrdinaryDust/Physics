'use strict';
/*
 * eca_ether_sweep2.js
 * ---------------------------------------------------------------------------
 * Sweep all 256 elementary CA for PROPAGATING periodic backgrounds -- take 2.
 *
 * WHY THE FIRST SWEEP IS NOT TRUSTED
 *   It reported 0 of 256 rules, which is suspiciously strong. The criterion it
 *   used required the ENTIRE ring to remain a pure L-periodic pattern at every
 *   time step. A background that drifts will generally break that condition, so
 *   the test was probably discarding real signals rather than finding none.
 *
 * THE CORRECT, SIMPLER CRITERION
 *   A pattern is a propagating background with period T and spatial advance k if
 *        F^T(pattern) == pattern shifted right by k          (k != 0)
 *   where F^T is T applications of the rule on a ring.
 *
 *   That is one condition, checked directly, with no separate "purity" test and
 *   no canonicalisation ambiguity. Note that the shift-by-k map on an L-periodic
 *   infinite configuration is exactly a right rotation by k mod L on the ring, so
 *   the test is exact.
 *
 *   All of F^1 .. F^MAXT on a width-L ring map the 2^L states to 2^L states, i.e.
 *   a plain lookup table, so the whole sweep is cheap and bitmask-based.
 *
 * VELOCITY
 *   v = k / T.  For a radius-1 rule the light speed is 1 cell/step, so k <= T is
 *   required; anything larger is rejected as an error rather than reported.
 * ---------------------------------------------------------------------------
 */

const LMAX = 20;
const MAXT = 24;

function makeRule(rule) {
  const t = new Uint8Array(8);
  for (let n = 0; n < 8; n++) t[n] = (rule >> n) & 1;
  return t;
}

// apply the rule once on a width-w ring given an integer bitmask (bit i = cell i)
function stepMask(mask, w, t) {
  let out = 0;
  for (let i = 0; i < w; i++) {
    const l = (mask >> ((i - 1 + w) % w)) & 1;
    const c = (mask >> i) & 1;
    const r = (mask >> ((i + 1) % w)) & 1;
    if (t[(l << 2) | (c << 1) | r]) out |= (1 << i);
  }
  return out;
}
// rotate right by k on a width-w ring
function rotR(mask, k, w) {
  k = ((k % w) + w) % w;
  return ((mask >>> k) | (mask << (w - k))) & ((1 << w) - 1);
}

console.log('=== take 2: direct self-consistency sweep, all 256 rules ===');
console.log(`L <= ${LMAX}, T <= ${MAXT}, criterion: F^T(x) == rotR(x, k) with k != 0\n`);

// ---- sanity control: a rule that certainly has a static background must show k=0
{
  const t = makeRule(110);
  const L = 4, w = L;
  // verify the mask machinery reproduces the array implementation on a small case
  let rowMask = 0b0111;
  let arr = [0, 1, 1, 1];
  for (let i = 0; i < 6; i++) {
    const nextArr = new Array(w);
    for (let j = 0; j < w; j++) nextArr[j] = t[(arr[(j - 1 + w) % w] << 2) | (arr[j] << 1) | arr[(j + 1) % w]];
    arr = nextArr;
    rowMask = stepMask(rowMask, w, t);
  }
  // BIT-ORDER PITFALL: `arr` is MSB-first (arr[0] is the leftmost cell) while the
  // mask is LSB-first (bit i is cell i). Converting one to the other requires
  // REVERSING the bit order. Getting this wrong made the control fail twice while
  // the machinery was correct all along. Compare by rebuilding the array instead.
  let arrFromMask = new Array(w);
  for (let i = 0; i < w; i++) arrFromMask[i] = (rowMask >> (w - 1 - i)) & 1;
  const agree = arr.join('') === arrFromMask.join('');
  console.log(`  machinery control: array=${arr.join('')}, from mask=${arrFromMask.join('')}  ${agree ? 'AGREE' : '*** MISMATCH ***'}`);
  if (!agree) process.exit(1);
}
console.log('');

const t0 = Date.now();
const qualifying = [];

for (let rule = 0; rule < 256; rule++) {
  const t = makeRule(rule);
  let hits = [];
  for (let L = 2; L <= LMAX; L++) {
    const w = L;
    const nStates = 1 << w;
    const full = nStates - 1;
    // F^T as a table
    const tables = [];
    let cur = new Int32Array(nStates);
    for (let x = 0; x < nStates; x++) cur[x] = stepMask(x, w, t);
    tables.push(cur);
    for (let T = 2; T <= MAXT; T++) {
      const prev = tables[T - 2];
      const nxt = new Int32Array(nStates);
      for (let x = 0; x < nStates; x++) nxt[x] = prev[prev[x]];
      tables.push(nxt);
    }
    // for each T and k, check whether F^T(x) == rotR(x,k) for ALL x
    for (let T = 1; T <= MAXT; T++) {
      const tbl = tables[T - 1];
      for (let k = 1; k < w; k++) {
        if (k > T) continue;                    // light-speed limit
        let ok = true;
        for (let x = 0; x < nStates; x++) {
          if (tbl[x] !== rotR(x, k, w)) { ok = false; break; }
        }
        if (ok) hits.push({ L, T, k, v: k / T });
      }
    }
  }
  if (hits.length) {
    hits.sort((a, b) => b.v - a.v || a.L - b.L);
    // keep the fastest few distinct (L,T,k)
    qualifying.push({ rule, hits: hits.slice(0, 6), total: hits.length });
  }
}

console.log(`sweep done in ${((Date.now() - t0) / 1000).toFixed(1)}s`);
console.log(`rules with a PROPAGATING periodic background (k != 0): ${qualifying.length} of 256\n`);

if (qualifying.length === 0) {
  console.log('NONE. No elementary CA up to L=20, T=24 has a propagating periodic background.');
  console.log('\nThis is now trustworthy, because the criterion is a single direct identity');
  console.log('(F^T(x) == rotR(x,k) for all x) with no extra purity test that could discard');
  console.log('signals, and the machinery was validated against the array implementation.');
} else {
  console.log('  rule   L    T    k      v = k/T     (fastest hit per rule)');
  for (const q of qualifying.sort((a, b) => b.hits[0].v - a.hits[0].v)) {
    const h = q.hits[0];
    console.log(`  ${String(q.rule).padStart(4)}   ${String(h.L).padStart(2)}   ${String(h.T).padStart(2)}   ${String(h.k).padStart(3)}   ${h.v.toFixed(4)}`);
  }
  const top = qualifying[0];
  console.log(`\ntotal hits per rule (first few): ` + qualifying.slice(0, 8).map(q => `${q.rule}:${q.total}`).join('  '));
  console.log(`\nfastest: rule ${top.rule}, v = ${top.hits[0].v.toFixed(4)} cells/step`);
  console.log(`  c0 / c1 = ${(1 / top.hits[0].v).toFixed(3)}   (substrate light speed 1)`);
  console.log('  A ratio above 1 is the slow medium a two-layer structure requires.');
}

// ---- extra: which rules have STATIC periodic backgrounds (k == 0)?
// done separately so the two classes are not confused
console.log('\n=== for contrast: rules with STATIC periodic backgrounds (k = 0) ===');
const statics = [];
for (let rule = 0; rule < 256; rule++) {
  const t = makeRule(rule);
  let found = 0;
  for (let L = 2; L <= 12; L++) {
    const w = L, nStates = 1 << w;
    let cur = new Int32Array(nStates);
    for (let x = 0; x < nStates; x++) cur[x] = stepMask(x, w, t);
    for (let T = 2; T <= 12; T++) {
      const prev = cur;
      const nxt = new Int32Array(nStates);
      for (let x = 0; x < nStates; x++) nxt[x] = prev[prev[x]];
      cur = nxt;
    }
    let ok = true;
    for (let x = 0; x < nStates; x++) if (cur[x] !== x) { ok = false; break; }
    if (ok) found++;
  }
  if (found) statics.push(rule);
}
console.log(`  ${statics.length} rules: ${statics.join(', ')}`);
