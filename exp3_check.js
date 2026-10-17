'use strict';
/*
 * exp3_check.js
 * ---------------------------------------------------------------------------
 * Why did EXPERIMENT 3 vary?
 *
 * Experiment 3 shifted the glider's start by k cells along its own direction of
 * travel, with a FIXED blinker. Mass over k = 0..7 was
 *     127, 96, 122, 109, 128, 153, 148, 146
 *
 * Naive reasoning says it must be identical: shifting by k cells shortens the
 * flight by 4k steps, and the blinker has period 2, so 4k is a multiple of 2 and
 * the blinker's phase at the collision should be unchanged.
 *
 * That reasoning has to be wrong somewhere. This file finds out where, by
 * measuring the actual arrival moment and the blinker's phase at that moment.
 *
 * METHOD
 *   For each k, record for every step t:
 *     - the number of components
 *     - the blinker's own canonical shape (identify the component that is the
 *       blinker, i.e. the one whose shape is 3x1 or 1x3 sitting at the fixed spot)
 *     - the first step at which the glider and blinker MERGE into one component
 *   Then check the blinker's phase at the merge step for each k.
 * ---------------------------------------------------------------------------
 */

const SIZE = 160;
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
  return { w, h, n: cells.length, minx, miny, cx: sx / cells.length, cy: sy / cells.length, id: `${w}x${h}:` + g.map(r => r.join('')).join('/') };
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }

const GLIDER0 = [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]];
function gliderAt(x, y) {
  const s = new Set();
  for (const [dx, dy] of GLIDER0) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function blinkerAt(cx, cy, horizontal) {
  const s = new Set();
  for (let d = -1; d <= 1; d++) {
    if (horizontal) s.add(wrapCell(cx + d, cy));
    else s.add(wrapCell(cx, cy + d));
  }
  return s;
}

const BLINKER_CX = 62, BLINKER_CY = 62;
const BLINKER_SHAPES = new Set();
// the blinker only ever has these two shapes
BLINKER_SHAPES.add('3x1:ooo');
BLINKER_SHAPES.add('1x3:o/o/o');

console.log('=== for each k, when does the glider reach the blinker, and what phase is it in? ===\n');
console.log('   k   mergeStep   blinker shape at merge   blinker shape 1 step before   final mass');

for (let k = 0; k < 8; k++) {
  const obstacle = blinkerAt(BLINKER_CX, BLINKER_CY, false);
  const start = gliderAt(32 + k, 32 + k);
  let s = union(obstacle, start);

  let mergeStep = -1;
  let shapeAtMerge = '-', shapeBefore = '-';
  let prevShape = '-';
  for (let t = 0; t <= 220; t++) {
    const comps = components(s);
    // find the component whose shape is a blinker shape AND whose position is near the blinker
    let blinkerShape = null;
    for (const c of comps) {
      const cc = canon(c);
      if (BLINKER_SHAPES.has(cc.id) && Math.abs(cc.cx - (BLINKER_CX + 1)) < 3 && Math.abs(cc.cy - (BLINKER_CY + 1)) < 3) {
        blinkerShape = cc.id;
        break;
      }
    }
    if (comps.length === 1 && mergeStep < 0) {
      mergeStep = t;
      shapeAtMerge = blinkerShape || '(merged, no separate blinker)';
      shapeBefore = prevShape;
    }
    prevShape = blinkerShape || prevShape;
    s = step(s);
  }

  // final mass after 200 steps
  let f = union(obstacle, start);
  for (let t = 0; t < 200; t++) f = step(f);

  console.log(`  ${String(k).padStart(2)}   ${String(mergeStep).padStart(9)}   ${String(shapeAtMerge).padEnd(24)} ${String(shapeBefore).padEnd(28)} ${f.size}`);
}

console.log('\n=== interpretation ===');
console.log('  The blinker has only two shapes: 3x1 (horizontal) and 1x3 (vertical).');
console.log('  If the "merge step" advances by 4 per unit k (as the flight shortening');
console.log('  predicts) and the blinker has period 2, then the blinker shape AT the');
console.log('  merge alternates only if the merge step changes by an ODD number.');
console.log('  Read the table: does the shape at merge ever change with k?');
