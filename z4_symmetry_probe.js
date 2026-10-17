'use strict';
/*
 * z4_symmetry_probe.js
 * ---------------------------------------------------------------------------
 * A SENTINEL CHECK, not a conclusion.
 *
 * NEW IDEA UNDER TEST
 *   A complex structure is a linear operator J with J^2 = -1. Such an operator
 *   typically comes from a 90-degree rotation. In a discrete system, a four-fold
 *   (Z4) INTERNAL SYMMETRY is exactly what supplies J^2 = -1.
 *
 *   So the question "where does complex structure come from" can be replaced by
 *   the checkable question: "is there a stable particle whose internal symmetry
 *   is Z4?" If yes, its deformation-mode space may carry a complex structure
 *   WITHOUT it being inserted by hand.
 *
 * THIS FILE ONLY ASKS THE CHEAP PRELIMINARY
 *   Do small patterns with four-fold symmetry exist in Life AT ALL?
 *   - enumerate small patterns, test 90-degree rotational symmetry
 *   - for those, check whether they survive as stable objects
 *   - also test reflection symmetry, for contrast
 *
 *   If NO four-fold symmetric small pattern survives, this route is dead in Life
 *   and there is no point building the mode-space machinery on top of it.
 * ---------------------------------------------------------------------------
 */

const SIZE = 40;
function wrapCell(x, y) { return ((y % SIZE) + SIZE) % SIZE * SIZE + (((x % SIZE) + SIZE) % SIZE); }

function lifeStep(set) {
  const counts = new Map();
  for (const k of set) {
    const x = k % SIZE, y = (k - x) / SIZE;
    for (let dy = -1; dy <= 1; dy++) {
      for (let dx = -1; dx <= 1; dx++) {
        if (dx === 0 && dy === 0) continue;
        const nk = wrapCell(x + dx, y + dy);
        counts.set(nk, (counts.get(nk) || 0) + 1);
      }
    }
  }
  const out = new Set();
  for (const [k, c] of counts) if (c === 3 || (c === 2 && set.has(k))) out.add(k);
  return out;
}

// ---- symmetry tests on a pattern given as a set of (x,y) offsets
function normalise(cells) {
  let minx = Infinity, miny = Infinity;
  for (const [x, y] of cells) { if (x < minx) minx = x; if (y < miny) miny = y; }
  return cells.map(([x, y]) => [x - minx, y - miny]).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
}
function same(a, b) {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i][0] !== b[i][0] || a[i][1] !== b[i][1]) return false;
  return true;
}
function bbox(cells) {
  let w = 0, h = 0;
  for (const [x, y] of cells) { if (x + 1 > w) w = x + 1; if (y + 1 > h) h = y + 1; }
  return [w, h];
}
function rot90(cells, w, h) {
  // rotate the w x h box by 90 degrees: (x,y) -> (h-1-y, x)
  return normalise(cells.map(([x, y]) => [h - 1 - y, x]));
}
function reflX(cells, w) { return normalise(cells.map(([x, y]) => [w - 1 - x, y])); }
function reflY(cells, h) { return normalise(cells.map(([x, y]) => [x, h - 1 - y])); }

function hasZ4(cells) {
  let [w, h] = bbox(cells);
  let c = normalise(cells);
  for (let i = 0; i < 4; i++) {
    c = rot90(c, w, h);
    [w, h] = bbox(c);
    if (i === 3) return same(c, normalise(cells));
  }
  return false;
}
function hasRefl(cells) {
  const [w, h] = bbox(cells);
  const c = normalise(cells);
  return same(normalise(reflX(cells, w)), c) || same(normalise(reflY(cells, h)), c);
}

// ---- enumerate all patterns inside a 3x3 and a 4x4 box, test symmetry and survival
function enumerateBox(n) {
  const results = [];
  const total = 1 << (n * n);
  for (let mask = 1; mask < total; mask++) {
    const cells = [];
    for (let i = 0; i < n * n; i++) if ((mask >> i) & 1) cells.push([i % n, Math.floor(i / n)]);
    if (!hasZ4(cells)) continue;
    // place it away from the edges and test survival
    let set = new Set();
    for (const [x, y] of cells) set.add(wrapCell(15 + x, 15 + y));
    const start = [...set].sort((a, b) => a - b).join(',');
    let cur = set, period = -1, died = false;
    const seen = new Map();
    for (let t = 0; t <= 40; t++) {
      const sig = [...cur].sort((a, b) => a - b).join(',');
      if (seen.has(sig)) { period = t - seen.get(sig); break; }
      seen.set(sig, t);
      if (cur.size === 0) { died = true; break; }
      cur = lifeStep(cur);
    }
    if (died) continue;
    if (period < 0) continue;
    results.push({ cells, n, mask, period, size: cells.length });
  }
  return results;
}

console.log('=== sentinel: are there Z4-symmetric patterns in Life at all? ===\n');

for (const n of [3, 4, 5]) {
  const t0 = Date.now();
  const r = enumerateBox(n);
  console.log(`--- ${n}x${n} box: ${r.length} four-fold symmetric patterns that survive ---`);
  // group by period
  const byP = new Map();
  for (const x of r) byP.set(x.period, (byP.get(x.period) || 0) + 1);
  for (const [p, c] of [...byP.entries()].sort((a, b) => a[0] - b[0])) {
    console.log(`    period ${String(p).padStart(2)}: ${c} patterns`);
  }
  // print the smallest few explicitly
  r.sort((a, b) => a.size - b.size);
  for (const x of r.slice(0, 4)) {
    const [w, h] = bbox(x.cells);
    console.log(`    example: ${w}x${h}, ${x.size} cells, period ${x.period}`);
    const g = Array.from({ length: h }, () => new Array(w).fill('.'));
    for (const [cx, cy] of x.cells) g[cy][cx] = 'o';
    for (const row of g) console.log('      ' + row.join(''));
  }
  console.log(`  (${((Date.now() - t0) / 1000).toFixed(1)}s)\n`);
}

console.log('=== the decisive question for the Z4 route ===');
console.log('  A four-fold symmetric pattern that SURVIVES is necessary for the idea,');
console.log('  but not sufficient. The next requirement would be that the pattern\'s');
console.log('  DEFORMATION MODES close under the dynamics and that the surviving');
console.log('  symmetry acts as J with J^2 = -1 on that mode space.');
console.log('');
console.log('  Additionally, a STATIC four-fold symmetric still life is useless here:');
console.log('  its mode space is trivial (it does not move, so there is no internal');
console.log('  clock). What the idea needs is a four-fold symmetric MOVER or a');
console.log('  four-fold symmetric OSCILLATOR with a nontrivial mode space.');
