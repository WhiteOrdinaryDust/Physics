'use strict';
/*
 * eca_cycle_structure.js
 * ---------------------------------------------------------------------------
 * WHERE COULD A 90-DEGREE PHASE COME FROM?
 *
 * The remaining gap: what mechanism makes the substrate produce the structure
 * "two amplitudes with a relative phase of 90 degrees"?
 *
 * A concrete candidate. 90 degrees = a fourth root of unity = PERIOD 4.
 * In a cellular automaton the time evolution is a PERMUTATION of the
 * configurations, and a permutation forces complex eigenvalues according to its
 * cycle lengths:
 *
 *      cycle length 1  ->  eigenvalue 1
 *      cycle length 2  ->  +1, -1        (REAL ONLY: no phase at all)
 *      cycle length 3  ->  cube roots of unity      (phase 120 degrees)
 *      cycle length 4  ->  +i, -i                   (EXACTLY 90 degrees)
 *      cycle length 6  ->  sixth roots of unity     (60 degrees)
 *
 * So if the substrate contains a 4-cycle, the 90-degree phase is not inserted by
 * hand: it is forced by the mathematics of the permutation.
 *
 * THIS FILE computes, for every elementary CA rule, the CYCLE LENGTH
 * DISTRIBUTION of its time-evolution permutation on a small ring.
 *
 * CONTROLS FIRST (this project's recurring failure mode is skipping them):
 *      rule 51  : complement rule, F is an involution -> ALL cycles length 1 or 2
 *      rule 204 : identity                            -> all cycles length 1
 *      rule 0   : everything maps to the empty state  -> exactly one 1-cycle
 * ---------------------------------------------------------------------------
 */

function makeRule(rule) {
  const t = new Uint8Array(8);
  for (let n = 0; n < 8; n++) t[n] = (rule >> n) & 1;
  return t;
}
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

/*
 * Cycle length distribution of the functional graph x -> F(x) on w bits.
 * Every weakly connected component of a functional graph contains exactly one
 * cycle, so we find the cycles by walking with a visited marker.
 */
function cycleDistribution(rule, w) {
  const t = makeRule(rule);
  const n = 1 << w;
  const F = new Int32Array(n);
  for (let x = 0; x < n; x++) F[x] = stepMask(x, w, t);

  const color = new Uint8Array(n);   // 0 unseen, 1 on current path, 2 done
  const dist = new Map();
  for (let s = 0; s < n; s++) {
    if (color[s] !== 0) continue;
    // walk from s
    const path = [];
    let x = s;
    while (color[x] === 0) {
      color[x] = 1;
      path.push(x);
      x = F[x];
    }
    if (color[x] === 1) {
      // found a new cycle: the tail from the first occurrence of x
      const idx = path.indexOf(x);
      const len = path.length - idx;
      dist.set(len, (dist.get(len) || 0) + 1);
    }
    for (const p of path) color[p] = 2;
  }
  return dist;
}

// ============================================================ CONTROLS
console.log('=== CONTROLS (must pass before the scan is trusted) ===\n');
{
  const w = 8;
  const n = 1 << w;

  // rule 51: complement of each cell -> F(F(x)) = x, so all cycles have length 1 or 2
  {
    const d = cycleDistribution(51, w);
    const lens = [...d.keys()].sort((a, b) => a - b);
    const total = [...d.values()].reduce((a, b) => a + b, 0);
    const ok = lens.every(L => L === 1 || L === 2);
    // expected: how many 2-cycles? states with F(x) != x pair up: (n - fixed)/2
    const t = makeRule(51);
    let fixed = 0;
    for (let x = 0; x < n; x++) if (stepMask(x, w, t) === x) fixed++;
    const expected2 = (n - fixed) / 2;
    console.log(`  rule 51 (complement): cycle lengths present = ${lens.join(', ')}`);
    console.log(`      cycles of length 2: ${d.get(2) || 0}, expected ${expected2}   ${(d.get(2) || 0) === expected2 ? 'PASS' : '*** FAIL ***'}`);
    console.log(`      only lengths 1 and 2? ${ok ? 'PASS' : '*** FAIL ***'}`);
    if (!ok) process.exit(1);
  }

  // rule 204: identity
  {
    const d = cycleDistribution(204, w);
    console.log(`  rule 204 (identity): cycle lengths present = ${[...d.keys()].join(', ')}`);
    const ok = [...d.keys()].length === 1 && d.get(1) === n;
    console.log(`      all ${n} states are 1-cycles? ${ok ? 'PASS' : '*** FAIL ***'}`);
    if (!ok) process.exit(1);
  }

  // rule 0: everything maps to the empty state
  {
    const d = cycleDistribution(0, w);
    console.log(`  rule 0 (all die): cycle lengths present = ${[...d.keys()].join(', ')}`);
    const ok = [...d.keys()].length === 1 && d.get(1) === 1;
    console.log(`      exactly one 1-cycle? ${ok ? 'PASS' : '*** FAIL ***'}`);
    if (!ok) process.exit(1);
  }
}
console.log('\n  all controls passed\n');

// ============================================================ THE SCAN
const WIDTH = 10;
console.log(`=== cycle length distributions, all 256 elementary CA (ring width ${WIDTH}) ===\n`);

const results = [];
for (let rule = 0; rule < 256; rule++) {
  const d = cycleDistribution(rule, WIDTH);
  const lens = [...d.keys()].sort((a, b) => a - b);
  const maxLen = Math.max(...lens);
  results.push({ rule, dist: d, lens, maxLen });
}

// which rules have a 4-cycle?
const has4 = results.filter(r => r.dist.has(4));
console.log(`rules whose permutation contains a 4-CYCLE (=> eigenvalues +/- i, i.e. 90 deg): ${has4.length} of 256`);
console.log('');
console.log('  rule    cycle lengths present          number of 4-cycles');
for (const r of has4) {
  console.log(`  ${String(r.rule).padStart(4)}    ${r.lens.join(', ').padEnd(30)} ${r.dist.get(4)}`);
}

// longest cycle over all rules
console.log('\n=== maximum cycle length, per rule (top 20) ===');
const byMax = [...results].sort((a, b) => b.maxLen - a.maxLen);
console.log('  rule   max cycle length   lengths present');
for (const r of byMax.slice(0, 20)) {
  console.log(`  ${String(r.rule).padStart(4)}   ${String(r.maxLen).padStart(16)}   ${r.lens.join(', ')}`);
}

// how many rules have ONLY real eigenvalues (all cycles length 1 or 2)?
const onlyReal = results.filter(r => r.lens.every(L => L === 1 || L === 2));
console.log(`\n=== rules with NO phase at all (every cycle length 1 or 2): ${onlyReal.length} of 256 ===`);
console.log('  ' + onlyReal.map(r => r.rule).join(', '));

console.log('\n=== reading ===');
console.log('  A 4-cycle exists in the time-evolution permutation exactly when the rule\'s');
console.log('  dynamics contain a state that returns only after 4 steps and not earlier.');
console.log('  Its presence means the time-evolution operator has eigenvalues +/- i, i.e.');
console.log('  a 90-degree phase, WITHOUT anything being inserted by hand.');
console.log('');
console.log('  CAUTION, stated up front: having a 90-degree eigenvalue in the time');
console.log('  evolution is NOT the same as the observable quantity we computed earlier');
console.log('  needing it. This file only answers whether the RAW MATERIAL exists.');
console.log('  Connecting it to the observable is a separate step and is NOT done here.');
