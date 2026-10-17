'use strict';
/*
 * find_ether.js
 * ---------------------------------------------------------------------------
 * FIND the rule-110 ether empirically instead of trusting a hand-written pattern.
 *
 * The previous run settled it: the 14-cell pattern I had been using as "the
 * ether" is NOT invariant under rule 110. Its density jumps from 5 to 10 ones
 * between t=0 and t=1, and its canonical form shows no period within 60 steps.
 * Therefore the "period 7, covering 88.5%" result from the timescale scan was a
 * TRANSIENT of a wrong initial pattern, not a property of the ether. That result
 * is retracted.
 *
 * WHAT AN ETHER ACTUALLY IS
 *   A background that is invariant up to a shift, i.e. there exist T and d with
 *        pattern(t + T) = pattern(t) shifted by d
 *   Then T is the background's time period and d/T its phase velocity.
 *
 * We search for it: try every 14-cell pattern of small density that is invariant
 * up to shift over 100 steps, seeded on a ring.
 * ---------------------------------------------------------------------------
 */

const W = 280;                    // multiple of 14
const TMAX = 100;

function step110(row, w) {
  const out = new Uint8Array(w);
  for (let i = 0; i < w; i++) {
    const idx = (row[(i - 1 + w) % w] << 2) | (row[i] << 1) | row[(i + 1) % w];
    out[i] = (110 >> idx) & 1;
  }
  return out;
}

function isShiftInvariant(init, len, periods) {
  // returns {T, shift} if pattern(t+T) is a rotation of pattern(t) for all t
  let row = Uint8Array.from(init);
  const snaps = [];
  for (let t = 0; t <= TMAX; t++) { snaps.push(Uint8Array.from(row)); row = step110(row, W); }
  const pat = (t) => { let s = ''; for (let i = 0; i < len; i++) s += snaps[t][i]; return s; };
  const isRot = (a, b) => {
    if (a.length !== b.length) return false;
    return (b + b).includes(a);
  };
  for (const T of periods) {
    let ok = true;
    for (let t = 0; t + T <= TMAX; t++) if (!isRot(pat(t), pat(t + T))) { ok = false; break; }
    if (ok) {
      // measure the shift: find d with pattern(t+T) == pattern(t) rotated by d
      const a = pat(0), b = pat(T);
      let d = -1;
      for (let k = 0; k < len; k++) {
        const rot = b.slice(k) + b.slice(0, k);
        if (rot === a) { d = k; break; }
      }
      return { T, shift: d };
    }
  }
  return null;
}

console.log('=== searching for a rule-110 ether (shift-invariant background) ===');
console.log(`ring ${W}, checking all 14-cell patterns of density <= 8, T up to 28\n`);

// enumerate 14-bit patterns with small density
const len = 14;
const found = [];
for (let mask = 0; mask < (1 << len); mask++) {
  let bits = 0;
  for (let i = 0; i < len; i++) if ((mask >> i) & 1) bits++;
  if (bits === 0 || bits > 8) continue;
  const init = new Uint8Array(W);
  for (let i = 0; i < W; i++) init[i] = (mask >> (i % len)) & 1;
  const periods = [];
  for (let T = 1; T <= 28; T++) periods.push(T);
  const r = isShiftInvariant(init, len, periods);
  if (r) {
    let s = '';
    for (let i = 0; i < len; i++) s += (mask >> i) & 1;
    found.push({ pattern: s, density: bits, T: r.T, shift: r.shift, v: r.shift / r.T });
  }
}
found.sort((a, b) => a.density - b.density || a.T - b.T);
console.log(`shift-invariant 14-cell patterns found: ${found.length}`);
for (const f of found.slice(0, 25)) {
  console.log(`  ${f.pattern}  ones=${f.density}  T=${String(f.T).padStart(2)}  shift=${String(f.shift).padStart(2)}  phase velocity=${f.v.toFixed(4)}`);
}

// ---- density distribution of what was found
if (found.length) {
  const dens = new Map();
  for (const f of found) dens.set(f.density, (dens.get(f.density) || 0) + 1);
  console.log('\n  densities of shift-invariant patterns:');
  for (const [d, c] of [...dens.entries()].sort((a, b) => a[0] - b[0])) {
    console.log(`    density ${d}: ${c} patterns`);
  }
  const sparsest = found[0];
  console.log(`\n  sparsest invariant background: ${sparsest.pattern} (density ${sparsest.density})`);
  console.log(`  its time period ${sparsest.T}, phase shift ${sparsest.shift}, phase velocity ${sparsest.v.toFixed(4)} cells/step`);
} else {
  console.log('\n  NO shift-invariant 14-cell pattern found. The background of rule 110, if any,');
  console.log('  is not of period 14, or is not reachable from a strictly periodic seed.');
}
