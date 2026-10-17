'use strict';
/*
 * slit_experiment2.js
 * ---------------------------------------------------------------------------
 * THE SLIT, DONE PROPERLY.
 *
 * WHAT WAS WRONG IN THE FIRST ATTEMPT
 *   The wall was a chain of 2x2 blocks laid along the anti-diagonal. In Life,
 *   blocks that touch DIAGONALLY merge into one connected component (8-connected),
 *   so the "chain" became a single blob of ~75 cells and the "gap width" was
 *   meaningless. Its control printed "wall blocks before: 1", which is the tell.
 *   All of that run's numbers are void.
 *
 * THE GEOMETRIC LESSON FROM THAT FAILURE
 *   A slit cannot be built by laying separate still lifes in a row: they merge.
 *   A slit requires ONE large stable pattern with a hole in it.
 *
 * WHAT THIS FILE DOES
 *   1. Build a large candidate wall (a solid slab), then RE-SHAPE it into a
 *      genuinely static pattern by running it to a fixed point first. If the
 *      slab is not static we will see that immediately.
 *   2. Measure the wall: cell count, number of connected components, and whether
 *      it is static. Printing these BEFORE any glider is fired is the control.
 *   3. Cut a slit by removing a band of cells, then RE-CHECK that the holed
 *      structure is still static.
 *   4. Only then fire a glider at the slit.
 *
 * PREDICTIONS (stated before running):
 *   - the plain slab is NOT static (Life slabs fill in)
 *   - therefore we first let it settle, and use the SETTLED pattern as the wall
 *   - a settled still-life cluster with a slit cut in it will normally become
 *     non-static, because the cut exposes unstable edges
 *   - if so, that is itself the result: you cannot cut a slit in ordinary
 *     still-life matter, which constrains the whole idea
 * ---------------------------------------------------------------------------
 */

const SIZE = 120;
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
function components(set) {
  const seen = new Set(), comps = [], Q = [];
  for (const start of set) {
    if (seen.has(start)) continue;
    Q.length = 0; Q.push(start); seen.add(start);
    const cells = [];
    while (Q.length) {
      const k = Q.pop(); cells.push(k);
      const x = k % SIZE, y = (k - x) / SIZE;
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        if (dx === 0 && dy === 0) continue;
        const nk = wrapCell(x + dx, y + dy);
        if (set.has(nk) && !seen.has(nk)) { seen.add(nk); Q.push(nk); }
      }
    }
    comps.push(cells);
  }
  return comps;
}
function canon(cells) {
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity, sx = 0, sy = 0;
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
    sx += x; sy += y;
  }
  const w = maxx - minx + 1, h = maxy - miny + 1;
  const g = [];
  for (let r = 0; r < h; r++) g.push(new Array(w).fill('.'));
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    g[y - miny][x - minx] = 'o';
  }
  return { w, h, n: cells.length, sx, sy, cx: sx / cells.length, cy: sy / cells.length, id: `${w}x${h}:` + g.map(r => r.join('')).join('/') };
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }
function signature(set) { return [...set].sort((a, b) => a - b).join(','); }

const GLIDER0 = [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]];
function gliderAt(x, y, turns) {
  const s = new Set();
  for (const [dx0, dy0] of GLIDER0) {
    let a = dx0, b = dy0;
    for (let i = 0; i < turns; i++) { const na = -b, nb = a; a = na; b = nb; }
    s.add(wrapCell(x + a, y + b));
  }
  return s;
}

// settle a pattern until it repeats, return the period-1 or final state
function settle(set, maxSteps) {
  let s = set;
  const seen = new Map();
  for (let t = 0; t <= maxSteps; t++) {
    const sig = signature(s);
    if (seen.has(sig)) return { state: s, period: t - seen.get(sig), steps: t, died: s.size === 0 };
    seen.set(sig, t);
    if (s.size === 0) return { state: s, period: 0, steps: t, died: true };
    s = step(s);
  }
  return { state: s, period: -1, steps: maxSteps, died: false };
}

function describe(name, set) {
  const comps = components(set);
  const kinds = comps.map(c => { const k = canon(c); return `${k.w}x${k.h}(${k.n})`; });
  console.log(`  ${name}: ${set.size} cells, ${comps.length} component(s)`);
  const sample = kinds.slice(0, 12).join(' ');
  console.log(`      shapes: ${sample}${kinds.length > 12 ? ` ... (+${kinds.length - 12} more)` : ''}`);
  return comps.length;
}

console.log('=== STEP 1: can a slab be a static wall at all? ===');
{
  // a solid 40 x 6 slab
  let slab = new Set();
  for (let x = 40; x < 80; x++) for (let y = 40; y < 46; y++) slab.add(wrapCell(x, y));
  console.log(`  initial slab: ${slab.size} cells`);
  const before = signature(slab);
  let s = slab;
  for (let t = 0; t < 30; t++) s = step(s);
  console.log(`  after 30 steps: ${s.size} cells`);
  console.log(`  static? ${signature(s) === before ? 'YES' : 'NO (a solid slab is not a still life)'}`);
  describe('settled slab', s);
  const st = settle(slab, 400);
  console.log(`  settling result: period ${st.period}, ${st.state.size} cells after ${st.steps} steps`);
  describe('settled state', st.state);
}

console.log('\n=== STEP 2: build a wall from a KNOWN still life, repeated WITHOUT touching ===');
{
  // blocks separated by at least one empty cell are separate components,
  // but in Life a block and its diagonal neighbour DO merge. So we space them
  // two apart along a line and check the component count.
  const parts = [];
  for (let i = 0; i < 12; i++) {
    const x = 30 + 3 * i, y = 60;
    const b = new Set();
    for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) b.add(wrapCell(x + dx, y + dy));
    parts.push(b);
  }
  const line = union(...parts);
  console.log(`  12 blocks spaced 3 apart along a horizontal line:`);
  describe('the line', line);
  const st = settle(line, 100);
  console.log(`  static? ${st.period === 1 ? 'YES' : `no, period ${st.period}`}`);
}

console.log('\n=== STEP 3: cut a slit, and check whether the wall survives the cut ===');
{
  // use a horizontal line of blocks, spaced 3 apart (verified separate above),
  // and REMOVE some blocks to make a gap. Then check the gap is a real opening:
  // measure the vertical extent of free space at the gap.
  const buildLine = (skipFrom, skipTo) => {
    let line = new Set();
    for (let i = 0; i < 20; i++) {
      if (i >= skipFrom && i < skipTo) continue;
      const x = 20 + 3 * i, y = 60;
      for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) line.add(wrapCell(x + dx, y + dy));
    }
    return line;
  };
  for (const gapBlocks of [1, 2, 3]) {
    const wall = buildLine(10, 10 + gapBlocks);
    const comps = components(wall).length;
    // gap width in cells: the empty span between the last block of the left part
    // and the first of the right part
    const leftLast = 20 + 3 * 9 + 1;          // right edge of block 9
    const rightFirst = 20 + 3 * (10 + gapBlocks); // left edge of the next block
    const gapCells = rightFirst - leftLast - 1;
    console.log(`  removing ${gapBlocks} block(s): ${comps} components, ${wall.size} cells, gap = ${gapCells} empty columns`);
  }
  console.log('');
  console.log('  NOTE: a "gap" between two SEPARATE blocks is not a slit in a wall.');
  console.log('  It is just a space between two objects. For a true slit the wall must be');
  console.log('  ONE connected object with a hole, and Life slabs are not static, so the');
  console.log('  only connected static wall available is a CLUSTER OF TOUCHING STILL LIFES.');
  console.log('  Whether such a cluster can have a hole without becoming unstable is the');
  console.log('  real question, and it is tested next.');
}

console.log('\n=== STEP 4: can a connected still-life cluster have a hole? ===');
{
  // Build a 2-cell-thick ring (annulus) of blocks, joined edge-to-edge, leaving a
  // hole in the middle. Then check stability and whether the hole is enclosed.
  const R = 7;
  let ring = new Set();
  for (let i = -R; i <= R; i++) {
    for (let j = -R; j <= R; j++) {
      const onRing = (Math.abs(i) === R || Math.abs(j) === R);
      if (!onRing) continue;
      ring.add(wrapCell(60 + i, 60 + j));
    }
  }
  console.log(`  solid square ring (outer ${2 * R + 1}, 1 cell thick): ${ring.size} cells`);
  describe('ring', ring);
  const st = settle(ring, 200);
  console.log(`  settled: period ${st.period}, ${st.state.size} cells`);
  describe('settled ring', st.state);

  // now a 2-thick ring so it might be more stable
  let ring2 = new Set();
  for (let i = -R; i <= R; i++) for (let j = -R; j <= R; j++) {
    const d = Math.max(Math.abs(i), Math.abs(j));
    if (d === R || d === R - 1) ring2.add(wrapCell(60 + i, 60 + j));
  }
  console.log(`\n  2-cell-thick square ring: ${ring2.size} cells`);
  const st2 = settle(ring2, 200);
  console.log(`  settled: period ${st2.period}, ${st2.state.size} cells`);
  describe('settled ring2', st2.state);
}

console.log('\n=== what this establishes ===');
console.log('  This file deliberately stops before firing a glider, because the wall');
console.log('  itself has to be built and verified first. The earlier failure came from');
console.log('  skipping exactly this step.');
console.log('');
console.log('  Read the outputs above: if no connected static pattern with an enclosed');
console.log('  hole can be built in Life, then the slit experiment as proposed cannot be');
console.log('  set up with still lifes at all, and that is a real constraint on the idea.');
