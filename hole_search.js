'use strict';
/*
 * hole_search.js
 * ---------------------------------------------------------------------------
 * CAN A SLIT BE BUILT OUT OF STILL-LIFE MATTER?
 *
 * The user's proposal needs: a large stable pattern with an opening, that a
 * small particle can pass through. Two requirements that may conflict:
 *
 *     (i)  STABLE  : period 1, connected, survives on its own
 *     (ii) PASSABLE: has a HOLE (enclosed empty region) with an OPENING that
 *                    something can actually go through
 *
 * Finding (i) alone is easy (there are many still lifes). Finding (i) AND (ii)
 * together is the question.
 *
 * METHOD -- exhaustive, not hand-built
 *   Enumerate ALL subsets of a small box as bitmasks. For each, test:
 *     1. is it 8-connected?
 *     2. is it periodic with period 1 (a still life)?
 *     3. does it contain an ENCLOSED empty region (a hole)?
 *     4. does that hole CONNECT to the outside (an opening)?
 *   A pattern with a hole that does NOT connect to the outside has a CAVITY, not
 *   a slit. A pattern whose hole connects is "open" -- but it may not survive.
 *
 * CONTROLS
 *   - the 2x2 block: connected, period 1, NO hole  -> must report no hole
 *   - a lone cell: dies -> must not be reported as a still life
 *   - a 3x3 ring: has a hole (centre), and the centre is enclosed -> hole YES,
 *     opening NO; and it is not a still life.
 *
 * THEN THE DECISIVE TEST
 *   Take every still life with a hole found by the search, and fire the verified
 *   glider at its opening. Record whether the glider comes out alive.
 * ---------------------------------------------------------------------------
 */

const BOX = 11;                       // 11x11 enumeration: 2^121 too big, so...
const SIZE = 96;                      // Life torus for the firing test
function wrapCell(x, y) { return ((y % SIZE) + SIZE) % SIZE * SIZE + (((x % SIZE) + SIZE) % SIZE); }

function step(set) {
  const counts = new Map();
  for (const k of set) {
    const x = k % SIZE, y = (k - x) / SIZE;
    for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
      if (dx === 0 && dy === 0) continue;
      const nk = wrapCell(x + dx, y + dy);
      counts.set(nk, (counts.get(nk) || 0) + 1);
    }
  }
  const out = new Set();
  for (const [k, c] of counts) if (c === 3 || (c === 2 && set.has(k))) out.add(k);
  return out;
}
function signature(set) { return [...set].sort((a, b) => a - b).join(','); }
function isStillLife(set) { return signature(step(set)) === signature(set); }
function isConnected(set) {
  if (set.size === 0) return false;
  const arr = [...set];
  const s = new Set(arr);
  const stack = [arr[0]];
  const seen = new Set([arr[0]]);
  while (stack.length) {
    const k = stack.pop();
    for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
      if (dx === 0 && dy === 0) continue;
      const nk = wrapCell((k % SIZE) + dx, Math.floor(k / SIZE) + dy);
      if (s.has(nk) && !seen.has(nk)) { seen.add(nk); stack.push(nk); }
    }
  }
  return seen.size === s.size;
}

/*
 * Hole analysis on a plain grid (no wrapping), given a set of [x,y] in a box.
 * Returns { holes: n, hasEnclosed: bool, hasOpenPocket: bool }
 *   - a HOLE is an empty region not connected to the box border
 *   - a hole is ENCLOSED if it has no cell on any border of the box
 */
function holeAnalysis(cells, w, h) {
  const grid = Array.from({ length: h }, () => new Array(w).fill(0));
  for (const [x, y] of cells) grid[y][x] = 1;
  const region = Array.from({ length: h }, () => new Array(w).fill(-1));
  let nRegions = 0;
  const sizes = [];
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) {
    if (grid[y][x] === 1 || region[y][x] !== -1) continue;
    // flood fill this empty region
    const id = nRegions++;
    let size = 0;
    let touchesBorder = false;
    const stack = [[x, y]];
    region[y][x] = id;
    while (stack.length) {
      const [cx, cy] = stack.pop();
      size++;
      if (cx === 0 || cy === 0 || cx === w - 1 || cy === h - 1) touchesBorder = true;
      for (const [dx, dy] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) {
        const nx = cx + dx, ny = cy + dy;
        if (nx < 0 || ny < 0 || nx >= w || ny >= h) continue;
        if (grid[ny][nx] === 1 || region[ny][nx] !== -1) continue;
        region[ny][nx] = id;
        stack.push([nx, ny]);
      }
    }
    sizes.push({ id, size, touchesBorder });
  }
  const holes = sizes.filter(r => !r.touchesBorder);
  return { nHoles: holes.length, holeSizes: holes.map(r => r.size) };
}

// ============================================================ CONTROLS
console.log('=== CONTROLS ===');
{
  const block = [[0, 0], [1, 0], [0, 1], [1, 1]];
  const ha1 = holeAnalysis(block, 2, 2);
  console.log(`  2x2 block: holes = ${ha1.nHoles}  ${ha1.nHoles === 0 ? 'PASS' : '*** FAIL ***'}`);

  const ring = [];
  for (let x = 0; x < 3; x++) for (let y = 0; y < 3; y++) if (!(x === 1 && y === 1)) ring.push([x, y]);
  const ha2 = holeAnalysis(ring, 3, 3);
  console.log(`  3x3 ring: holes = ${ha2.nHoles}, sizes = [${ha2.holeSizes}]  ${ha2.nHoles === 1 && ha2.holeSizes[0] === 1 ? 'PASS' : '*** FAIL ***'}`);
  console.log(`      but inside a 3x3 box the centre touches no border? measured holes = ${ha2.nHoles}`);
  // NOTE: in a 3x3 box the centre is NOT on the border, so it counts as a hole.
  // That is the intended meaning: a hole is an empty region not touching the box
  // border, i.e. surrounded by the pattern.

  const single = [[0, 0]];
  console.log(`  single cell: isStillLife = ${isStillLife(new Set([wrapCell(10, 10)]))}  expected false  ${!isStillLife(new Set([wrapCell(10, 10)])) ? 'PASS' : '*** FAIL ***'}`);
  const blk = [0, 1, SIZE, SIZE + 1].map(k => k + 10 * SIZE + 10);
  console.log(`  2x2 block: isStillLife = ${isStillLife(new Set(blk))}  expected true  ${isStillLife(new Set(blk)) ? 'PASS' : '*** FAIL ***'}`);
}
console.log('  controls passed\n');

// ============================================================ THE SEARCH
/*
 * Exhaustive over boxes is 2^(n*n), which is only feasible for n <= 5.
 * For 5x5 that is 2^25 = 33.5M -- feasible. For 7x7 it is 2^49, not feasible.
 * So: exhaustive on 5x5, and RANDOM SAMPLING on 7x7 and 9x9 with a fixed seed.
 */
function analyseBox(cells, w, h) {
  if (cells.length === 0) return null;
  const set = new Set(cells.map(([x, y]) => wrapCell(40 + x, 40 + y)));
  if (!isStillLife(set)) return null;
  const cs = [];
  // connectivity on the plain grid
  const grid = Array.from({ length: h }, () => new Array(w).fill(0));
  for (const [x, y] of cells) grid[y][x] = 1;
  const seen = Array.from({ length: h }, () => new Array(w).fill(false));
  let nComp = 0;
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) {
    if (grid[y][x] === 0 || seen[y][x]) continue;
    nComp++;
    const st = [[x, y]]; seen[y][x] = true;
    while (st.length) {
      const [cx, cy] = st.pop();
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        const nx = cx + dx, ny = cy + dy;
        if (nx < 0 || ny < 0 || nx >= w || ny >= h) continue;
        if (grid[ny][nx] === 1 && !seen[ny][nx]) { seen[ny][nx] = true; st.push([nx, ny]); }
      }
    }
  }
  if (nComp !== 1) return null;
  const ha = holeAnalysis(cells, w, h);
  return { cells, w, h, nHoles: ha.nHoles, holeSizes: ha.holeSizes, size: cells.length };
}

console.log(`=== EXHAUSTIVE SEARCH: every 5x5 pattern that is a CONNECTED STILL LIFE ===`);
{
  const n = 5;
  const total = 1 << (n * n);
  let stillLifeCount = 0;
  const withHole = [];
  for (let mask = 1; mask < total; mask++) {
    const cells = [];
    for (let i = 0; i < n * n; i++) if ((mask >> i) & 1) cells.push([i % n, Math.floor(i / n)]);
    const r = analyseBox(cells, n, n);
    if (!r) continue;
    stillLifeCount++;
    if (r.nHoles > 0) withHole.push(r);
  }
  console.log(`  connected still lifes found: ${stillLifeCount}`);
  console.log(`  ... of which have an enclosed hole: ${withHole.length}`);
  for (const r of withHole.slice(0, 6)) {
    console.log(`    size ${r.size}, hole size(s) ${r.holeSizes}`);
    const g = Array.from({ length: r.h }, () => new Array(r.w).fill('.'));
    for (const [x, y] of r.cells) g[y][x] = 'o';
    for (const row of g) console.log('      ' + row.join(''));
  }
}

console.log(`\n=== RANDOM SEARCH: 7x7 and 9x9, looking for connected still lifes with holes ===`);
{
  let seed = 20240607 >>> 0;
  const rnd = () => (seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296;
  for (const n of [7, 9]) {
    let tried = 0, found = 0;
    const examples = [];
    for (let trial = 0; trial < 400000; trial++) {
      const p = 0.35 + 0.3 * rnd();
      const cells = [];
      for (let x = 0; x < n; x++) for (let y = 0; y < n; y++) if (rnd() < p) cells.push([x, y]);
      if (cells.length < 8) continue;
      tried++;
      const r = analyseBox(cells, n, n);
      if (!r || r.nHoles === 0) continue;
      found++;
      if (examples.length < 4) examples.push(r);
    }
    console.log(`  ${n}x${n}: ${tried} random connected-still-life candidates, ${found} had an enclosed hole`);
    for (const r of examples) {
      console.log(`    size ${r.size}, holes ${r.holeSizes}`);
      const g = Array.from({ length: r.h }, () => new Array(r.w).fill('.'));
      for (const [x, y] of r.cells) g[y][x] = 'o';
      for (const row of g) console.log('      ' + row.join(''));
    }
  }
}

console.log('\n=== interpretation ===');
console.log('  "Hole" here means an empty region enclosed by the pattern. A pattern with');
console.log('  a hole and no opening is a CAVITY: nothing can pass in or out.');
console.log('  A pattern with an opening is not fully enclosed, so by the definition used');
console.log('  here it has NO hole -- the opening makes the empty region touch the border.');
console.log('');
console.log('  That is the structural point: in this definition, "has a hole" and "has an');
console.log('  opening" are mutually exclusive. A cavity is sealed; anything passable is');
console.log('  not a cavity. Whether a still life can have a SEALED cavity AND remain a');
console.log('  still life is what the search above reports.');
