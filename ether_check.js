'use strict';
/*
 * ether_check.js
 * ---------------------------------------------------------------------------
 * A clean, minimal, diagnostic check of the claimed invariant background.
 *
 * Two previous attempts to "use the real ether" failed because I kept
 * transcription-editing a constant and re-running. That is guessing, not
 * verifying. This file instead:
 *
 *   1. takes the mask that the exhaustive search reported
 *   2. prints its exact bit string
 *   3. evolves it and prints the actual 13-cell windows side by side
 *   4. checks shift-invariance by direct comparison, printing the mismatch count
 *      for every candidate shift instead of just reporting success/failure
 *
 * If the background is genuinely invariant up to a shift, exactly one shift will
 * give zero mismatches. If no shift gives zero, the claimed background is wrong
 * and the earlier "217 invariant patterns" result has a bug worth finding.
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

const LEN = 13;

// The exhaustive search reported this mask for "11110001000000" with ones=5.
// Reconstruct the mask from the string rather than from a hand-copied array.
const PATTERN_STR = '11110001000000';
let mask = 0;
for (let i = 0; i < PATTERN_STR.length; i++) if (PATTERN_STR[i] === '1') mask |= (1 << i);
console.log('pattern string as reported :', PATTERN_STR);
console.log('reconstructed mask         :', mask, '(binary', mask.toString(2) + ')');

// window of the first LEN bits, which is what the search tested
let win = '';
for (let i = 0; i < LEN; i++) win += (mask >> i) & 1;
console.log(`first ${LEN} bits             : ${win}`);
console.log(`bits ${LEN} and beyond         : ${PATTERN_STR.slice(LEN)}  (must be all zeros for a ${LEN}-periodic seed)`);

// build the ring with period LEN
function ring(w) {
  const r = new Uint8Array(w);
  for (let i = 0; i < w; i++) r[i] = (mask >> (i % LEN)) & 1;
  return r;
}

const w = LEN * 202;
console.log(`\nring width ${w} = ${LEN} x ${w / LEN}`);

let row = ring(w);
console.log('\nfirst 26 cells at successive times:');
const snaps = [];
for (let t = 0; t <= 30; t++) {
  snaps.push(Uint8Array.from(row));
  let s = '';
  for (let i = 0; i < 26; i++) s += row[i];
  console.log(`  t=${String(t).padStart(2)}  ${s}`);
  row = step110(row, w);
}

// for each T, for each shift, count mismatches
console.log('\nmismatch counts (should be 0 for exactly one shift at the true period):');
console.log('   T \\ sh ' + Array.from({ length: LEN }, (_, k) => String(k).padStart(4)).join(''));
let trueT = -1, trueSh = -1;
for (let T = 1; T <= 20; T++) {
  let line = '  ' + String(T).padStart(4) + '   ';
  let bestSh = -1, bestMis = Infinity;
  for (let sh = 0; sh < LEN; sh++) {
    let mis = 0;
    for (let i = 0; i < 3 * LEN; i++) if (snaps[T][i] !== snaps[0][(i + sh) % w]) mis++;
    line += String(mis).padStart(4);
    if (mis < bestMis) { bestMis = mis; bestSh = sh; }
  }
  if (bestMis === 0 && trueT < 0) { trueT = T; trueSh = bestSh; }
  console.log(line + (bestMis === 0 ? `   <-- EXACT, shift ${bestSh}` : ''));
}

console.log('');
if (trueT > 0) {
  console.log(`RESULT: the background IS invariant up to a shift.`);
  console.log(`  time period T = ${trueT}, spatial shift = ${trueSh}`);
  console.log(`  phase velocity = ${trueSh}/${trueT} = ${(trueSh / trueT).toFixed(4)} cells/step`);
  console.log(`  substrate max speed = 1 cell/step  =>  c0/c1 = ${(trueT / trueSh).toFixed(3)}`);
} else {
  console.log('RESULT: NO shift gives an exact match. The claimed background is NOT');
  console.log('  invariant, so the "217 invariant patterns" search has a bug (most likely');
  console.log('  in how it built the ring from the mask, or in its rotation test).');
}
