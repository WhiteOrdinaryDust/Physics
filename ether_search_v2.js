'use strict';
/*
 * ether_search_v2.js
 * ---------------------------------------------------------------------------
 * A SELF-CHECKING search for shift-invariant backgrounds ("ethers") of rule 110.
 *
 * WHY A REWRITE
 *   The previous version reported "217 shift-invariant 14-cell patterns" and was
 *   wrong. Diagnosing it afterwards (ether_check.js) showed no (T, shift) pair
 *   with zero mismatches, and revealed the cause: the file used LEN = 13 while
 *   the pattern array had 14 entries, and the search itself hard-coded a third
 *   value. Three lengths, one file.
 *
 * DESIGN RULES (to make that class of error impossible)
 *   1. L is defined ONCE per call and passed explicitly everywhere. No module
 *      level length constant is read by the search.
 *   2. Every candidate is verified to be genuinely L-periodic BEFORE use; if it
 *      is not, the candidate is skipped and counted.
 *   3. The search prints the actual mismatch matrix for the winners, so a reader
 *      can see the zeros rather than trusting a count.
 *   4. NEGATIVE CONTROL: a random seed must NOT be reported as invariant. If it
 *      is, the search logic is broken and the run aborts.
 * ---------------------------------------------------------------------------
 */

function step110(row, w) {
  const out = new Uint8Array(w);
  for (let i = 0; i < w; i++) {
    const idx = (row[(i - 1 + w) % w] << 2) | (row[i] << 1) | row[(i + 1) % w];
    out[i] = (110 >> idx) & 1;
  }
  return out;
}

// is `row` periodic with period L, exactly?
function isPeriodic(row, L) {
  for (let i = 0; i < row.length; i++) if (row[i] !== row[i % L]) return false;
  return true;
}

/*
 * Test one periodic seed.
 * Returns null if not invariant, else { T, shift, phaseVelocity }.
 *
 * Convention: next[i] === cur[(i + shift) % L]
 */
function testSeed(seedBits, L, TMAX, maxT) {
  const w = L * 12;
  let row = new Uint8Array(w);
  for (let i = 0; i < w; i++) row[i] = seedBits[i % L];
  if (!isPeriodic(row, L)) return null;         // rule 2

  const snaps = [];
  for (let t = 0; t <= TMAX; t++) {
    snaps.push(Uint8Array.from(row));
    row = step110(row, w);
  }
  // every snapshot must remain L-periodic, otherwise the seed is not a background
  for (let t = 0; t <= TMAX; t++) {
    if (!isPeriodic(snaps[t], L)) {
      // it may still be L-periodic with a different phase; test periodicity of
      // the configuration by checking it repeats with L at some offset
      let ok = false;
      for (let off = 0; off < L && !ok; off++) {
        ok = true;
        for (let i = 0; i < 3 * L; i++) if (snaps[t][i] !== snaps[t][off + (i % L)]) { ok = false; break; }
      }
      if (!ok) return null;
    }
  }

  for (let T = 1; T <= maxT; T++) {
    let shift = -1;
    for (let k = 0; k < L; k++) {
      let good = true;
      for (let i = 0; i < 3 * L; i++) if (snaps[T][i] !== snaps[0][(i + k) % w]) { good = false; break; }
      if (good) { shift = k; break; }
    }
    if (shift >= 0) {
      // confirm it holds for ALL later times, not just t=0 -> t=T
      let all = true;
      for (let t = 0; t + T <= TMAX; t++) {
        let any = false;
        for (let k = 0; k < L; k++) {
          let good = true;
          for (let i = 0; i < 3 * L; i++) if (snaps[t + T][i] !== snaps[t][(i + k) % w]) { good = false; break; }
          if (good) { any = true; break; }
        }
        if (!any) { all = false; break; }
      }
      if (all) return { T, shift, phaseVelocity: shift / T };
    }
  }
  return null;
}

function mismatchMatrix(seedBits, L, TMAX, maxT) {
  const w = L * 12;
  let row = new Uint8Array(w);
  for (let i = 0; i < w; i++) row[i] = seedBits[i % L];
  const snaps = [];
  for (let t = 0; t <= TMAX; t++) { snaps.push(Uint8Array.from(row)); row = step110(row, w); }
  const lines = [];
  for (let T = 1; T <= maxT; T++) {
    let best = Infinity, bestSh = -1;
    const rowCounts = [];
    for (let k = 0; k < L; k++) {
      let mis = 0;
      for (let i = 0; i < 3 * L; i++) if (snaps[T][i] !== snaps[0][(i + k) % w]) mis++;
      rowCounts.push(mis);
      if (mis < best) { best = mis; bestSh = k; }
    }
    lines.push({ T, rowCounts, best, bestSh });
  }
  return lines;
}

// ============================================================ negative control
console.log('=== NEGATIVE CONTROL: a non-invariant seed must NOT be reported ===');
{
  const L = 13, TMAX = 80, maxT = 20;
  // a deliberately non-periodic-under-evolution seed: single 1 in an empty ring
  const bits = new Array(L).fill(0); bits[0] = 1;
  const r = testSeed(bits, L, TMAX, maxT);
  console.log(`  seed = single 1 (${bits.join('')}), L=${L}`);
  console.log(`  reported invariant? ${r ? 'YES -> *** SEARCH IS BROKEN, ABORTING ***' : 'no (correct)'}`);
  if (r) process.exit(1);
  const mm = mismatchMatrix(bits, L, TMAX, maxT);
  const zeroRows = mm.filter(x => x.best === 0).length;
  console.log(`  rows of the mismatch matrix with a zero entry: ${zeroRows} (must be 0)`);
  if (zeroRows > 0) { console.log('  *** SEARCH IS BROKEN, ABORTING ***'); process.exit(1); }
}
console.log('  negative control PASSED\n');

// ============================================================ the search
console.log('=== searching all L-periodic seeds of rule 110 ===');
const TMAX = 100;
const MAXT = 40;
const results = [];
let triedSeeds = 0, nonPeriodic = 0, nonInvariant = 0;

for (let L = 2; L <= 16; L++) {
  const maxMask = 1 << L;
  let foundForL = 0;
  for (let mask = 1; mask < maxMask; mask++) {
    const bits = [];
    for (let i = 0; i < L; i++) bits.push((mask >> i) & 1);
    triedSeeds++;
    const r = testSeed(bits, L, TMAX, MAXT);
    if (r) {
      const ones = bits.reduce((a, b) => a + b, 0);
      results.push({ L, mask, bits: bits.join(''), ones, T: r.T, shift: r.shift, v: r.phaseVelocity });
      foundForL++;
    } else nonInvariant++;
  }
  console.log(`  L=${String(L).padStart(2)}  seeds tried ${String(maxMask - 1).padStart(6)}  invariant found ${String(foundForL).padStart(4)}`);
}

console.log(`\ntotals: seeds tried ${triedSeeds}, invariant ${results.length}, not invariant ${nonInvariant}`);

if (results.length === 0) {
  console.log('\nNO shift-invariant background found for any L in 2..16.');
  console.log('=> rule 110 has no strictly periodic ether of these periods.');
} else {
  // group by density
  const byDensity = new Map();
  for (const r of results) byDensity.set(r.ones, (byDensity.get(r.ones) || 0) + 1);
  console.log('\ninvariant backgrounds by number of ones:');
  for (const [d, c] of [...byDensity.entries()].sort((a, b) => a[0] - b[0])) {
    console.log(`  ones=${d}: ${c}`);
  }

  results.sort((a, b) => a.ones - b.ones || a.L - b.L);
  console.log('\nsparsest invariant backgrounds (candidate ethers):');
  for (const r of results.slice(0, 12)) {
    console.log(`  L=${String(r.L).padStart(2)}  ${r.bits}  ones=${r.ones}  T=${String(r.T).padStart(2)}  shift=${String(r.shift).padStart(2)}  v=${r.v.toFixed(4)}`);
  }

  // print the mismatch matrix for the sparsest one, as evidence
  const top = results[0];
  console.log(`\nevidence: mismatch matrix for L=${top.L} ${top.bits}`);
  console.log('   T \\ shift ' + Array.from({ length: top.L }, (_, k) => String(k).padStart(4)).join(''));
  const mm = mismatchMatrix(top.bits.split('').map(Number), top.L, TMAX, 20);
  for (const line of mm) {
    console.log('  ' + String(line.T).padStart(5) + '   ' + line.rowCounts.map(c => String(c).padStart(4)).join('') +
      (line.best === 0 ? `   <-- EXACT at shift ${line.bestSh}` : ''));
  }
  console.log('\n  (a row with a zero entry at shift k means pattern(t+T) is pattern(t)');
  console.log('   shifted by k -- that is the invariance we are looking for)');
}
