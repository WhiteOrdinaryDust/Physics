'use strict';
/*
 * find_orthogonal_ship.js
 * ---------------------------------------------------------------------------
 * WHY: the LWSS generator in collision.js is BROKEN. Its own control prints
 *      "*** FAIL: expected 1 component, got 2 ***" -- the 9-cell pattern is two
 *      pieces, not a ship. (That failure was read as a pass in an earlier round;
 *      re-running collision.js shows it plainly.) Everything that needs a
 *      particle moving in a straight cardinal direction depends on this, and a
 *      particle moving perpendicular to a screen is the whole point of a slit
 *      experiment, so the ship has to be found, not remembered.
 *
 * METHOD: brute force. Enumerate every pattern in a 5x4 box (2^20 masks, which
 *   covers the LWSS bounding box), keep those of 7..13 cells, evolve 8
 *   generations, and keep the ones whose shape at t=4 and t=8 is the SAME shape
 *   displaced by a constant non-zero integer vector. No hand-written pattern,
 *   no literature lookup: the rule itself decides.
 *
 * OUTPUT: every distinct ship found, with its period, velocity, component count
 *   per phase (so a two-piece impostor is visible as such), and the cell list in
 *   a form that can be pasted into another script.
 * ---------------------------------------------------------------------------
 */

const STRIDE = 64, OFF = 16, S2 = STRIDE * STRIDE;
const cbuf = new Int32Array(S2), mark = new Int32Array(S2), lmark = new Int32Array(S2);
let stamp = 0, lstamp = 0;

function kOf(x, y) { return (y + OFF) * STRIDE + (x + OFF); }
function stepCells(cells) {
  stamp++; lstamp++;
  for (const k of cells) lmark[k] = lstamp;
  const touched = [];
  for (const k of cells) {
    const x = k % STRIDE, y = (k - x) / STRIDE;
    for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
      if (dx === 0 && dy === 0) continue;
      const nk = (y + dy) * STRIDE + (x + dx);
      if (mark[nk] !== stamp) { mark[nk] = stamp; cbuf[nk] = 0; touched.push(nk); }
      cbuf[nk]++;
    }
  }
  const out = [];
  for (const nk of touched) {
    const c = cbuf[nk];
    if (c === 3 || (c === 2 && lmark[nk] === lstamp)) out.push(nk);
  }
  return out;
}
function canonKeys(cells) {
  let minx = 1e9, miny = 1e9, maxx = -1e9, maxy = -1e9, sx = 0, sy = 0;
  for (const k of cells) {
    const x = k % STRIDE, y = (k - x) / STRIDE;
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
    sx += x; sy += y;
  }
  const w = maxx - minx + 1, h = maxy - miny + 1, g = [];
  for (let r = 0; r < h; r++) g.push(new Array(w).fill('.'));
  for (const k of cells) {
    const x = k % STRIDE, y = (k - x) / STRIDE;
    g[y - miny][x - minx] = 'o';
  }
  return { id: `${w}x${h}:` + g.map(r => r.join('')).join('/'), w, h, n: cells.length, cx: sx / cells.length, cy: sy / cells.length };
}
function compCount(cells) {
  const live = new Set(cells), seen = new Set();
  let n = 0;
  for (const s of live) {
    if (seen.has(s)) continue;
    n++; const Q = [s]; seen.add(s);
    while (Q.length) {
      const k = Q.pop(), x = k % STRIDE, y = (k - x) / STRIDE;
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        if (dx === 0 && dy === 0) continue;
        const nk = (y + dy) * STRIDE + (x + dx);
        if (live.has(nk) && !seen.has(nk)) { seen.add(nk); Q.push(nk); }
      }
    }
  }
  return n;
}

const W = 5, H = 4, BITS = W * H;
console.log(`=== brute force: patterns in a ${W}x${H} box that travel in a straight line ===`);
console.log(`  ${1 << BITS} masks to test\n`);

const found = new Map();
let pass = 0;
const t0 = Date.now();

for (let mask = 1; mask < (1 << BITS); mask++) {
  let m = mask, pc = 0;
  while (m) { m &= m - 1; pc++; }
  if (pc < 6 || pc > 13) continue;
  const cells = [];
  for (let b = 0; b < BITS; b++) if (mask & (1 << b)) cells.push(kOf(16 + (b % W), 16 + ((b / W) | 0)));
  const c0 = canonKeys(cells);
  const fr = [cells];
  let cur = cells, blewUp = false;
  for (let t = 0; t < 8; t++) {
    cur = stepCells(cur);
    if (cur.length > 40) { blewUp = true; break; }
    fr.push(cur);
  }
  if (blewUp) continue;
  const c4 = canonKeys(fr[4]), c8 = canonKeys(fr[8]);
  if (c0.id !== c4.id || c0.id !== c8.id) continue;
  const d1x = Math.round(c4.cx - c0.cx), d1y = Math.round(c4.cy - c0.cy);
  const d2x = Math.round(c8.cx - c4.cx), d2y = Math.round(c8.cy - c4.cy);
  if (d1x !== d2x || d1y !== d2y) continue;
  if (!((Math.abs(d1x) === 2 && d1y === 0) || (Math.abs(d1y) === 2 && d1x === 0))) continue;
  pass++;
  const key = `${c0.id}|${d1x},${d1y}`;
  if (!found.has(key)) {
    const comps = [0, 1, 2, 3].map(t => compCount(fr[t]));
    found.set(key, { id: c0.id, n: c0.n, w: c0.w, h: c0.h, vx: d1x, vy: d1y, comps, mask, cells: [...cells] });
  }
  if (Date.now() - t0 > 240000) { console.log('  (time limit reached, stopping the enumeration)'); break; }
}

console.log(`masks that are period-4 straight-line travellers: ${pass}`);
console.log(`distinct shapes: ${found.size}\n`);
const list = [...found.values()].sort((a, b) => a.n - b.n);
for (const f of list) {
  console.log(`--- ${f.n} cells, ${f.w}x${f.h}, period 4, velocity (${f.vx},${f.vy})`);
  console.log(`    pieces per phase t=0..3: ${f.comps.join(', ')}${f.comps.every(c => c === 1) ? '   (a genuine single particle)' : '   <-- NOT a single particle at every phase'}`);
  console.log(`    pattern:`);
  for (const row of f.id.split(':')[1].split('/')) console.log(`      ${row}`);
  const rel = [];
  const minx = Math.min(...f.cells.map(k => k % STRIDE)), miny = Math.min(...f.cells.map(k => (k - (k % STRIDE)) / STRIDE));
  for (const k of f.cells) { const x = k % STRIDE, y = (k - x) / STRIDE; rel.push([x - minx, y - miny]); }
  rel.sort((a, b) => a[1] - b[1] || a[0] - b[0]);
  console.log(`    cells: [${rel.map(p => `[${p[0]},${p[1]}]`).join(',')}]`);
}
console.log(`\nelapsed ${((Date.now() - t0) / 1000).toFixed(1)} s`);
