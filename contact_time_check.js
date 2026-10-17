'use strict';
/*
 * contact_time_check.js
 * ---------------------------------------------------------------------------
 * THE CONFOUND THAT phase_arrival_check.js EXPOSED, AND HOW TO SETTLE IT.
 *
 *   The four phases of the ship do NOT sit in the same bounding box at the same
 *   generation in general. Because the ship advances 2 cells per 4 generations
 *   and wobbles in y on the way, the boxes coincide only when t is a multiple of
 *   4 (measured from the launch). At t=112,117,... wait -- in the run, the boxes
 *   coincided at t = 112, 116, 120, 124 (every 4th generation) and differed in
 *   between.
 *
 *   So "the phase changed the outcome" in lwss_slit.js could be nothing but a
 *   sub-lattice position difference: the colliding particle may simply be sitting
 *   1 cell lower or higher and 1 cell further along when it first disturbs the
 *   screen. That is an impact parameter, not a hidden internal state.
 *
 * THE EXACT TEST (no bounding-box proxies, no hand-waving)
 *   The ship alone is a faithful description of the particle up to the instant
 *   the screen's evolution is first affected, because until then the two do not
 *   influence each other. So, step by step, compare
 *        step(ship U wall)      with      step(ship) U step(wall)
 *   The first generation t* where these differ is the exact instant the screen
 *   learns anything at all. At that instant we record
 *        - the ship's bounding box,
 *        - which screen cells were affected (that is the record appearing),
 *        - whether that cell is on the LOWER or the UPPER edge of the slit.
 *
 *   If all four phases have the same t* and the same box at t*, then the four
 *   runs are the same impact parameter and the difference between them is purely
 *   internal -- the hidden phase. If t* or the box differs, the phase effect is
 *   confounded with position and must be re-measured with the launch moved.
 * ---------------------------------------------------------------------------
 */

const SIZE = 600;
function wrapCell(x, y) { return ((((y % SIZE) + SIZE) % SIZE) * SIZE) + (((x % SIZE) + SIZE) % SIZE); }
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
function canon(cells) {
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
  }
  return { minx, miny, maxx, maxy };
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }
function setEq(a, b) { if (a.size !== b.size) return false; for (const k of a) if (!b.has(k)) return false; return true; }
function shipRight(x, y) {
  const s = new Set();
  for (const [dx, dy] of [[1, 0], [2, 0], [3, 0], [4, 0], [0, 1], [4, 1], [4, 2], [0, 3], [3, 3]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function normalizeTo(set, x0, y0) {
  const c = canon(set);
  const dx = x0 - c.minx, dy = y0 - c.miny;
  const out = new Set();
  for (const k of set) { const x = k % SIZE, y = (k - x) / SIZE; out.add(wrapCell(x + dx, y + dy)); }
  return out;
}
function blockAt(x, y) {
  const s = new Set();
  for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function wallBlocks(X0, Y0, nx, ny, s, jskip) {
  let m = new Set();
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
    if (j === jskip) continue;
    m = union(m, blockAt(X0 + i * s, Y0 + j * s));
  }
  return m;
}

const X0 = 250, Y0 = 150, NY = 9, JSKIP = 4, LAUNCH_X = 190;

function firstContact(wall, y0, phase) {
  let sh = shipRight(150, 150);
  for (let i = 0; i < phase; i++) sh = step(sh);
  sh = normalizeTo(sh, LAUNCH_X, y0);
  for (let t = 0; t < 200; t++) {
    const both = step(union(sh, wall));
    const free = union(step(sh), wall ? step(wall) : new Set());
    if (!setEq(both, free)) {
      const affected = [];
      for (const k of both) if (!free.has(k)) affected.push(k);
      for (const k of free) if (!both.has(k)) affected.push(k);
      return { t, box: canon(sh), affected, mass: sh.size };
    }
    sh = step(sh);
  }
  return null;
}

console.log('=== WHEN does the screen first learn something, and WHERE? ===');
console.log('  the screen is the spacing-6, 1-column wall; slit rows 170..179');
console.log('  lower edge block row occupies y=168,169;  upper edge block row y=180,181');
console.log('  (a ship cell is in range of the wall when it is within 1 of those rows)\n');
console.log('  thickness spacing  b    phase   t*    ship box at t*        mass  first screen cells touched');
for (const [nx, s, bs] of [[1, 6, [2, 3, 4]], [1, 4, [0, 1, 2]], [2, 6, [2, 3]], [3, 6, [2, 3]]]) {
  const wall = wallBlocks(X0, Y0, nx, NY, s, JSKIP);
  const gapLo = Y0 + s * (JSKIP - 1) + 2;
  for (const b of bs) {
    const rows = [];
    for (let phase = 0; phase < 4; phase++) {
      const r = firstContact(wall, gapLo + b, phase);
      if (!r) { rows.push(null); continue; }
      rows.push(r);
      const cells = r.affected.map(k => `(${k % SIZE},${Math.floor(k / SIZE)})`).slice(0, 6).join(' ');
      console.log(`  ${String(nx).padStart(9)} ${String(s).padStart(7)} ${String(b).padStart(4)}   ${String(phase).padStart(5)}  ${String(r.t).padStart(4)}    ${`${r.box.minx}..${r.box.maxx},${r.box.miny}..${r.box.maxy}`.padEnd(20)} ${String(r.mass).padStart(4)}  ${cells}`);
    }
    const ts = new Set(rows.filter(Boolean).map(r => r.t));
    const boxes = new Set(rows.filter(Boolean).map(r => `${r.box.minx}..${r.box.maxx},${r.box.miny}..${r.box.maxy}`));
    console.log(`      -> distinct t*: ${[...ts].join(',')}   distinct boxes: ${[...boxes].join(' | ')}`);
    console.log(`      -> ${ts.size === 1 && boxes.size === 1 ? 'SAME impact instant and SAME box for all four phases: the difference is internal' : 'DIFFERENT instant or box: the phases do not arrive at the same place at the same time'}`);
    console.log('');
  }
}
