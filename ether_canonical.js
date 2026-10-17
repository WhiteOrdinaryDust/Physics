'use strict';
/*
 * ether_canonical.js
 * ---------------------------------------------------------------------------
 * Remove the shift ambiguity from the ether search, and print the candidates.
 *
 * PROBLEM WITH THE PREVIOUS OUTPUT
 *   It reported, for L=4, "T=1, shift=2, v=2.0" -- a phase velocity of 2 cells
 *   per step, which is impossible for a radius-1 rule (max 1 cell/step).
 *   The cause is an ALIASING AMBIGUITY: when the temporal period T is smaller
 *   than the spatial period L, several different (T, shift) pairs describe the
 *   same evolution, so "shift/T" is not well defined.
 *
 * FIX
 *   Work with the CANONICAL FORM (the lexicographically smallest rotation) of the
 *   pattern at each time. Then:
 *     - the temporal period T* is the smallest T with canon(t+T) == canon(t)
 *     - the spatial advance is the shift needed to map canon(t) to canon(t+T)
 *   Reporting T* removes the alias. Also require the measured velocity to be
 *   <= 1, and flag any candidate that violates that as an error rather than
 *   silently reporting it.
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

function canonStr(bits) {
  let best = null;
  for (let k = 0; k < bits.length; k++) {
    const r = bits.slice(k).concat(bits.slice(0, k)).join('');
    if (best === null || r < best) best = r;
  }
  return best;
}

/*
 * Analyse a periodic seed properly.
 * Returns null if not a background, else { T, shift, v, canon }.
 */
function analyse(bits, L, TMAX, maxT) {
  const w = L * 12;
  let row = new Uint8Array(w);
  for (let i = 0; i < w; i++) row[i] = bits[i % L];
  const snaps = [];
  for (let t = 0; t <= TMAX; t++) {
    snaps.push(Uint8Array.from(row));
    row = step110(row, w);
  }
  // canonical L-window at each time; if the state is not a pure L-periodic
  // pattern the seed is not a background
  const canonAt = [];
  for (let t = 0; t <= TMAX; t++) {
    const win = [];
    for (let i = 0; i < L; i++) win.push(snaps[t][i]);
    // verify the window repeats with period L over a long stretch
    let periodic = true;
    for (let i = 0; i < 3 * L; i++) if (snaps[t][i] !== win[i % L]) { periodic = false; break; }
    if (!periodic) return null;
    canonAt.push(canonStr(win));
  }
  // smallest temporal period of the canonical form, requiring it to hold for the
  // whole run, not just near the end
  for (let T = 1; T <= maxT; T++) {
    let ok = true;
    for (let t = 0; t + T <= TMAX; t++) if (canonAt[t] !== canonAt[t + T]) { ok = false; break; }
    if (!ok) continue;
    // spatial advance: shift k with canon(t) rotated by k giving canon(t+T)
    const a = canonAt[0], b = canonAt[T];
    let shift = -1;
    for (let k = 0; k < L; k++) {
      if (b.slice(k).concat(b.slice(0, k)) === a) { shift = k; break; }
    }
    if (shift < 0) continue;
    // a shift of k and a shift of k-L are the same physical displacement;
    // choose the representative of smallest magnitude
    let s = shift;
    if (s > L / 2) s -= L;
    const v = s / T;
    return { T, shift, signedShift: s, v, canon: a, alias: Math.abs(v) > 1 + 1e-9 };
  }
  return null;
}

console.log('=== canonical ether search (aliasing removed) ===\n');

const TMAX = 120, MAXT = 40;
const found = [];
let aliased = 0;

for (let L = 2; L <= 16; L++) {
  const maxMask = 1 << L;
  let n = 0;
  for (let mask = 1; mask < maxMask; mask++) {
    const bits = [];
    for (let i = 0; i < L; i++) bits.push((mask >> i) & 1);
    const r = analyse(bits, L, TMAX, MAXT);
    if (!r) continue;
    n++;
    if (r.alias) { aliased++; continue; }      // reject impossible velocities
    const ones = bits.reduce((a, b) => a + b, 0);
    found.push({ L, bits: bits.join(''), ones, ...r });
  }
  console.log(`  L=${String(L).padStart(2)}  found ${String(n).padStart(4)}  (aliased/rejected ${aliased} so far)`);
}

console.log(`\ntotal valid backgrounds: ${found.length}, rejected for velocity > 1: ${aliased}`);

// ---- de-duplicate by canonical form
const uniq = new Map();
for (const f of found) if (!uniq.has(f.canon)) uniq.set(f.canon, f);
console.log(`distinct canonical backgrounds: ${uniq.size}`);

// ---- the density-5 family (what I originally believed the ether to be)
console.log('\n=== density-5 backgrounds (the family I originally thought was the ether) ===');
const d5 = [...uniq.values()].filter(f => f.ones === 5).sort((a, b) => a.L - b.L || a.T - b.T);
console.log(`count: ${d5.length}`);
for (const f of d5.slice(0, 20)) {
  console.log(`  L=${String(f.L).padStart(2)}  canon=${f.canon}  T=${String(f.T).padStart(3)}  signedShift=${String(f.signedShift).padStart(3)}  v=${f.v.toFixed(4)}`);
}
if (d5.length === 0) console.log('  none found');

// ---- summary by density with velocity ranges
console.log('\n=== all densities: count, and velocity range ===');
const byD = new Map();
for (const f of uniq.values()) {
  if (!byD.has(f.ones)) byD.set(f.ones, []);
  byD.get(f.ones).push(f);
}
for (const [d, list] of [...byD.entries()].sort((a, b) => a[0] - b[0])) {
  const vs = list.map(x => x.v);
  const Ts = list.map(x => x.T);
  console.log(`  ones=${String(d).padStart(2)}  count=${String(list.length).padStart(4)}  T in [${Math.min(...Ts)}, ${Math.max(...Ts)}]  v in [${Math.min(...vs).toFixed(3)}, ${Math.max(...vs).toFixed(3)}]`);
}

// ---- global check: did any candidate exceed the speed limit?
const tooFast = [...uniq.values()].filter(f => Math.abs(f.v) > 1 + 1e-9);
console.log(`\nspeed-limit check: ${tooFast.length} backgrounds exceed 1 cell/step` +
  (tooFast.length ? '  *** contradiction with the radius-1 limit ***' : '  (all within the limit)'));
