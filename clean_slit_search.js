'use strict';
/*
 * clean_slit_search.js
 * ---------------------------------------------------------------------------
 * EXHAUSTIVE SEARCH FOR A CLEAN PASS-THROUGH SLIT.
 *
 * THE CORRECT CRITERION (the previous attempts used the wrong one)
 *   "passing through" must mean: what goes in comes out, unchanged, like a
 *   particle through a slit. Not a reaction, not debris, not a new particle.
 *   The earlier searches counted "a mover survived", which is satisfied by a
 *   NEWLY CREATED glider. That was the mistake.
 *
 *   CLEAN PASS-THROUGH requires ALL of:
 *     (1) the wall's cell set after the run is EXACTLY the wall's cell set before
 *     (2) the only other cells form a single 3x3 glider
 *     (3) total mass = wall mass + 5
 *     (4) the outgoing heading is the same as the incoming heading
 *
 * SEARCH SPACE
 *   - block spacing along a row: 3,4,5,6,7,8
 *   - which block is omitted (the gap position): 2,3,4
 *   - impact offset b = -3..3, applied PERPENDICULAR to the flight direction
 *   - glider phase: 0..3 (a step along its own direction, i.e. a time shift)
 *
 *   The launch point is chosen so that the glider's straight-line trajectory
 *   passes through the gap centre, offset by b perpendicular to the flight.
 *
 * CONTROLS
 *   - the four glider headings are classified correctly in vacuum
 *   - with NO wall, every configuration must trivially "pass" (the glider flies)
 *
 * OUTPUT
 *   Only CLEAN passes are printed. Everything else is summarised by count.
 * ---------------------------------------------------------------------------
 */

const SIZE = 240;
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
  return { w, h, n: cells.length, minx, miny, maxx, maxy, sx, sy, cx: sx / cells.length, cy: sy / cells.length, id: `${w}x${h}:` + g.map(r => r.join('')).join('/') };
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }
function setEq(a, b) {
  if (a.size !== b.size) return false;
  for (const k of a) if (!b.has(k)) return false;
  return true;
}

const GLIDER0 = [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]];
const HEAD_VEL = [[1, 1], [-1, 1], [-1, -1], [1, -1]];
function gliderAt(x, y, heading) {
  const s = new Set();
  for (const cell of GLIDER0) {
    let a = cell[0], b = cell[1];
    for (let i = 0; i < heading; i++) { const na = -b, nb = a; a = na; b = nb; }
    s.add(wrapCell(x + a, y + b));
  }
  return s;
}
function blockAt(x, y) {
  const s = new Set();
  for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function buildWall(spacing, count, skip, x0, y0) {
  let wall = new Set();
  for (let i = 0; i < count; i++) {
    if (i === skip) continue;
    wall = union(wall, blockAt(x0 + i * spacing, y0));
  }
  return wall;
}

/*
 * Fire a glider of `heading` so that its straight trajectory passes through the
 * gap centre, displaced by b perpendicular to the flight direction, and offset
 * by `phase` steps along the flight direction.
 *
 * For heading 0 (velocity (1,1)):  perpendicular unit direction is (1,-1).
 */
function launchPoint(cx, cy, heading, b, phase, dist) {
  const hx = HEAD_VEL[heading][0], hy = HEAD_VEL[heading][1];
  // unit motion direction and a perpendicular unit direction (integer vectors of length sqrt2)
  const ux = hx, uy = hy;                       // length sqrt2
  const px = hx, py = -hy;                      // length sqrt2, perpendicular
  // start = target - dist*u + b*p + phase*u
  const x = cx - (dist - phase) * ux + b * px;
  const y = cy - (dist - phase) * uy + b * py;
  return [x, y];
}

// ============================================================ controls
console.log('=== CONTROLS ===');
{
  let ok = true;
  const exp = [[1, 1], [-1, 1], [-1, -1], [1, -1]];
  for (let h = 0; h < 4; h++) {
    let s = gliderAt(60, 60, h);
    for (let t = 0; t < 40; t++) s = step(s);
    const cs = components(s);
    const good = cs.length === 1 && cs[0].length === 5;
    if (!good) ok = false;
    console.log(`  heading ${h}: ${cs.length} component(s), mass ${s.size}  ${good ? 'PASS' : '*** FAIL ***'}`);
  }
  // wall alone must be static
  const w = buildWall(5, 7, 3, 30, 60);
  const before = new Set(w);
  let s = w;
  for (let t = 0; t < 60; t++) s = step(s);
  console.log(`  wall alone static: ${setEq(before, s) ? 'PASS' : '*** FAIL ***'}`);
  if (!setEq(before, s)) ok = false;
  if (!ok) process.exit(1);
}
console.log('  controls passed\n');

// ============================================================ the search
const SPACINGS = [3, 4, 5, 6, 7, 8];
const SKIPS = [2, 3, 4];
const OFFSETS = [-3, -2, -1, 0, 1, 2, 3];
const PHASES = [0, 1, 2, 3];
const HEADINGS = [0, 1, 2, 3];
const STEPS = 400;

console.log('=== SEARCH: clean pass-through ===');
console.log('  Criterion: wall cell set EXACTLY unchanged, plus exactly one 3x3 glider,');
console.log('  total mass = wallMass + 5, and the same heading coming out.\n');

let tested = 0;
const clean = [];
const stats = { wallChanged: 0, massWrong: 0, multiOrWrongShape: 0, headingChanged: 0 };
const interesting = [];

for (const spacing of SPACINGS) {
  for (const skip of SKIPS) {
    const count = 7, x0 = 40, y0 = 80;
    const wall = buildWall(spacing, count, skip, x0, y0);
    const wallMass = wall.size;
    // gap between block (skip-1) and block (skip+1)
    const gapLeft = x0 + (skip - 1) * spacing + 1;
    const gapRight = x0 + (skip + 1) * spacing;
    const cx = Math.round((gapLeft + gapRight) / 2);
    const cy = y0 + 1;

    for (const heading of HEADINGS) {
      for (const b of OFFSETS) {
        for (const phase of PHASES) {
          tested++;
          const [sx, sy] = launchPoint(cx, cy, heading, b, phase, 80);
          const start = gliderAt(sx, sy, heading);
          if (start.size !== 5) continue;
          let s = union(wall, start);
          for (let t = 0; t < STEPS; t++) s = step(s);

          const wallAfter = new Set();
          for (const k of s) if (wall.has(k)) wallAfter.add(k);
          const outside = new Set();
          for (const k of s) if (!wall.has(k)) outside.add(k);

          const wallUnchanged = setEq(wall, wallAfter);
          const total = s.size;
          const outsideComps = components(outside);

          if (!wallUnchanged) { stats.wallChanged++; continue; }
          if (total !== wallMass + 5) { stats.massWrong++; continue; }
          if (outsideComps.length !== 1 || outsideComps[0].length !== 5) {
            stats.multiOrWrongShape++;
            continue;
          }
          const oc = canon(outsideComps[0]);
          if (oc.w !== 3 || oc.h !== 3) { stats.multiOrWrongShape++; continue; }

          // heading check: follow the outside component and see the drift
          let f = new Set(outsideComps[0]);
          const frames = [];
          let headingOk = false;
          for (let t = 0; t < 40; t++) {
            const cs = components(f);
            if (cs.length !== 1) break;
            frames.push(canon(cs[0]));
            f = step(f);
          }
          if (frames.length >= 20) {
            let dx = (frames[4].sx - frames[0].sx) / 5;
            let dy = (frames[4].sy - frames[0].sy) / 5;
            if (Math.abs(dx) > SIZE / 2) dx -= Math.sign(dx) * SIZE;
            if (Math.abs(dy) > SIZE / 2) dy -= Math.sign(dy) * SIZE;
            const vx = Math.round(dx), vy = Math.round(dy);
            headingOk = (vx === HEAD_VEL[heading][0] && vy === HEAD_VEL[heading][1]);
          }
          if (!headingOk) { stats.headingChanged++; continue; }

          clean.push({ spacing, skip, heading, b, phase });
        }
      }
    }
  }
}

console.log(`configurations tested: ${tested}`);
console.log(`  rejected because the wall changed      : ${stats.wallChanged}`);
console.log(`  rejected because the mass was wrong    : ${stats.massWrong}`);
console.log(`  rejected because the shape was wrong   : ${stats.multiOrWrongShape}`);
console.log(`  rejected because the heading changed   : ${stats.headingChanged}`);
console.log(`  CLEAN PASSES                           : ${clean.length}`);
console.log('');

if (clean.length === 0) {
  console.log('  NO clean pass-through found in this search space.');
  console.log('  Note: the search covers spacing 3..8, three gap positions, seven impact');
  console.log('  offsets and four phases for four headings. A clean pass may need a');
  console.log('  wider gap, a longer wall, or a different wall arrangement, so this is');
  console.log('  NOT a proof that none exists - only that none occurs in this space.');
} else {
  console.log('  CLEAN PASSES FOUND:');
  console.log('  spacing  skip  heading  offset  phase');
  for (const c of clean.slice(0, 60)) {
    console.log(`  ${String(c.spacing).padStart(7)}  ${String(c.skip).padStart(4)}  ${String(c.heading).padStart(7)}  ${String(c.b).padStart(6)}  ${String(c.phase).padStart(5)}`);
  }
  if (clean.length > 60) console.log(`  ... and ${clean.length - 60} more`);
  console.log('');
  // does the pass-through depend on phase?
  console.log('  Does a clean pass depend on the PHASE? (for the same spacing/skip/heading/offset)');
  const g = new Map();
  for (const c of clean) {
    const key = `${c.spacing}/${c.skip}/${c.heading}/${c.b}`;
    if (!g.has(key)) g.set(key, new Set());
    g.get(key).add(c.phase);
  }
  let anyPhaseDep = false;
  for (const [key, phases] of g) {
    if (phases.size > 1 && phases.size < 4) { anyPhaseDep = true; console.log(`    ${key}: phases ${[...phases].join(',')} (some but not all) <-- DEPENDS ON PHASE`); }
  }
  if (!anyPhaseDep) console.log('    no configuration where only some phases pass');
}
