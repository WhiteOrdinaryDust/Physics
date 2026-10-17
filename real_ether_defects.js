'use strict';
/*
 * real_ether_defects.js
 * ---------------------------------------------------------------------------
 * Re-do rule-110 particle detection on the REAL ether.
 *
 * The previous round found that the 14-cell "ether" used in all earlier rule-110
 * work is not an invariant background at all (its density jumps 5 -> 10 at the
 * first step). So the earlier rule-110 particle results -- including the v = 2/3
 * glider speed -- were measured against a wrong background and are untrusted.
 *
 * HERE
 *   - use the empirically found invariant background (density 5, T=13, shift 2)
 *   - build a large ring of it
 *   - scan for DEFECTS: contiguous blocks where the ring differs from the ether
 *     phase expected at that time (ether(t+d) where d = 2t/13, interpolated)
 *   - track defect blocks over time and measure their drift -> matter speed
 *   - compare with the ether phase velocity 2/13 = 0.1538, giving the ratio
 *     c1 / v_matter that the two-layer picture needs
 * ---------------------------------------------------------------------------
 */

const W = 2520;                   // 2520 = 13 * 193.8... use a multiple of 13
const WW = 2522;                  // 2522 = 13 * 194
const TMAX = 400;

function step110(row, w) {
  const out = new Uint8Array(w);
  for (let i = 0; i < w; i++) {
    const idx = (row[(i - 1 + w) % w] << 2) | (row[i] << 1) | row[(i + 1) % w];
    out[i] = (110 >> idx) & 1;
  }
  return out;
}

// The invariant background found empirically is the bit string
//     11110001000000   (bit i printed as (mask>>i)&1, so bit 0 first)
// = bits 0..4 = 1, bits 5..13 = 0  -> five ones then nine zeros.
// NOTE: an earlier version of this file transcribed it wrongly as
// [1,1,1,1,0,0,0,1,0,0,0,0,0,0] (an extra 1 at index 7), which broke every
// consistency check below. Keep the transcription literal.
const PHASE = [1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
const LEN = 13;                   // its true period is 13, not 14
// fill a ring with the background
function etherRing(w) {
  const r = new Uint8Array(w);
  for (let i = 0; i < w; i++) r[i] = PHASE[i % LEN];
  return r;
}

console.log('=== rule 110 on the REAL ether (period 13, phase velocity 2/13) ===\n');

// ---- confirm the background is invariant up to a shift, and measure it
{
  const w = 2626;                 // 13 * 202
  let row = etherRing(w);
  const snaps = [];
  for (let t = 0; t <= 40; t++) { snaps.push(Uint8Array.from(row)); row = step110(row, w); }
  // find (T, shift)
  let found = null;
  for (let T = 1; T <= 30 && !found; T++) {
    let ok = true;
    for (let t = 0; t + T <= 40; t++) {
      // is snaps[t+T] a rotation of snaps[t]?
      let matched = false;
      for (let sh = 0; sh < LEN; sh++) {
        let good = true;
        // by definition of the shift convention below, pattern(t+T) equals
        // pattern(t) shifted by sh, i.e. next[i] == cur[(i + sh) % w]
        for (let i = 0; i < 3 * LEN; i++) {
          if (snaps[t + T][i] !== snaps[t][(i + sh) % w]) { good = false; break; }
        }
        if (good) { matched = true; break; }
      }
      if (!matched) { ok = false; break; }
    }
    if (ok) found = T;
  }
  console.log(`background invariance: period T = ${found === null ? 'not found' : found}`);
  if (found) {
    // measure shift
    let sh = -1;
    for (let k = 0; k < LEN; k++) {
      let good = true;
      for (let i = 0; i < 3 * LEN; i++) if (snaps[found][i] !== snaps[0][(i + k) % w]) { good = false; break; }
      if (good) { sh = k; break; }
    }
    console.log(`  shift over T = ${sh} cells  -> phase velocity = ${(sh / found).toFixed(4)} cells/step`);
  }
}

// ---- scan for defects: cells that differ from the ether at the expected phase
function defectBlocks(row, expected) {
  const bad = new Uint8Array(row.length);
  let any = false;
  for (let i = 0; i < row.length; i++) { bad[i] = row[i] !== expected[i] ? 1 : 0; if (bad[i]) any = true; }
  if (!any) return [];
  let z = -1;
  for (let i = 0; i < row.length; i++) if (!bad[i]) { z = i; break; }
  if (z === -1) return [{ start: 0, len: row.length }];
  const out = [];
  let i = 0;
  while (i < row.length) {
    const idx = (z + i) % row.length;
    if (bad[idx]) {
      let n = 0;
      while (n < row.length && bad[(z + i + n) % row.length]) n++;
      out.push({ start: (z + i) % row.length, len: n });
      i += n;
    } else i++;
  }
  return out;
}

/*
 * Seed a local defect into the ether and watch it. The seed is a small block that
 * differs from the background.
 */
const seeds = [
  { name: 'block of 3', cells: [1, 1, 1], at: 0 },
  { name: 'block of 5', cells: [1, 0, 1, 0, 1], at: 0 },
  { name: 'glider-like triplet', cells: [1, 1, 1, 0, 1, 1, 1], at: 0 },
  { name: 'single 1', cells: [1], at: 0 },
  { name: 'single flip', cells: [1], at: 4, flip: true },
];

console.log('\n--- seeding local defects into the real ether and tracking them ---');
for (const sd of seeds) {
  const w = 2002;                  // 13 * 154
  let row = etherRing(w);
  const pos = Math.floor(w / 2);
  for (let j = 0; j < sd.cells.length; j++) {
    const i = (pos + j) % w;
    if (sd.flip) row[i] ^= 1; else row[i] = sd.cells[j];
  }
  const track = [];
  for (let t = 0; t <= TMAX; t++) {
    // expected ether at time t: ether(t + shift(t)) with shift(t) = round(2t/13)
    const sh = Math.round(2 * t / 13);
    const expected = new Uint8Array(w);
    for (let i = 0; i < w; i++) expected[i] = PHASE[(i + sh) % LEN];
    const blocks = defectBlocks(row, expected);
    const total = blocks.reduce((a, b) => a + b.len, 0);
    const big = blocks.filter(b => b.len >= 2);
    const centroid = big.length ? big.reduce((a, b) => a + (b.start + b.len / 2) * b.len, 0) / big.reduce((a, b) => a + b.len, 0) : -1;
    track.push({ t, nBlocks: blocks.length, total, centroid });
    if (t === TMAX) break;
    row = step110(row, w);
  }
  const alive = track.filter(x => x.total > 0 && x.total < 60);
  console.log(`\n  seed: ${sd.name}`);
  console.log(`    defect-cell count over time: t=0:${track[0].total}  t=50:${track[50].total}  t=200:${track[200].total}  t=400:${track[400].total}`);
  if (alive.length >= 20) {
    const a = alive[0], b = alive[alive.length - 1];
    let dc = b.centroid - a.centroid;
    if (dc > w / 2) dc -= w; if (dc < -w / 2) dc += w;
    console.log(`    tracked from t=${a.t} to t=${b.t}: centroid moved ${dc.toFixed(1)} cells`);
    console.log(`    => drift speed = ${(dc / (b.t - a.t)).toFixed(4)} cells/step`);
    console.log(`    => ratio to ether phase velocity (0.1538): ${(Math.abs(dc / (b.t - a.t)) / 0.1538).toFixed(3)}`);
  } else {
    console.log('    defect did not persist in a trackable form');
  }
}

console.log('\n=== reference numbers ===');
console.log('  substrate clock            : 1 step');
console.log('  substrate max signal speed : 1 cell/step');
console.log('  ether phase velocity c1    : 2/13 = 0.1538 cells/step');
console.log('  => c0/c1 = 6.5');
