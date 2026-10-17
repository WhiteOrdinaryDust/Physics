'use strict';
/*
 * dynamic_wall_test.js
 * ---------------------------------------------------------------------------
 * A DYNAMIC WALL.
 *
 * The user's requirement: the wall must itself run dynamically under the SAME
 * rule, not be a frozen arrangement. A still life is a legal Life configuration
 * but it has period 1, so it has no internal phase and therefore nothing for a
 * glider's phase to be (mis)matched against. A period-2 oscillator DOES have a
 * phase.
 *
 * DESIGN
 *   WALL   : a blinker, period 2, placed in the glider's path. Life has no
 *            background, so a blinker is the simplest genuinely dynamic barrier.
 *   CONTROL: a block, period 1, same position.
 *   GLIDER : launched along (+1,+1). Its starting position along its own
 *            direction of travel, taken mod 4, sets its phase at any later
 *            moment. Advancing the start by 1 cell reduces the flight time by
 *            4 steps, so the arrival PHASE on the wall shifts.
 *
 * WHAT IS MEASURED
 *   For each glider start offset k = 0..7, the terminal state of the whole
 *   system: final mass, number of components, movers, still lifes.
 *
 * PREDICTIONS, STATED BEFORE RUNNING
 *   1. CONTROL (block, period 1): the outcome is IDENTICAL for every k,
 *      because the obstacle has no phase to be matched.
 *   2. WALL (blinker, period 2): the outcome MAY depend on k. If it does, then
 *      two firings that look identical differ only by the hidden phase relation,
 *      which is the user's mechanism. If it does not, the mechanism fails here.
 *   3. The vacuum glider is identical for every k, after removing its position.
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
  return { w, h, n: cells.length, sx, sy, cx: sx / cells.length, cy: sy / cells.length, id: `${w}x${h}:` + g.map(r => r.join('')).join('/') };
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }
function signature(set) { return [...set].sort((a, b) => a - b).join(','); }

const GLIDER0 = [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]];
function gliderAt(x, y) {
  const s = new Set();
  for (const [dx, dy] of GLIDER0) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function blockAt(x, y) {
  const s = new Set();
  for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
// blinker, period 2. Horizontal phase by default; place its centre at (cx,cy).
function blinkerAt(cx, cy, horizontal) {
  const s = new Set();
  for (let d = -1; d <= 1; d++) {
    if (horizontal) s.add(wrapCell(cx + d, cy));
    else s.add(wrapCell(cx, cy + d));
  }
  return s;
}

/*
 * Run one experiment. The obstacle is given as a function of nothing (fixed in
 * space). The glider starts at (startX + k, startY + k), which both advances its
 * phase and shortens its flight by 4k steps.
 */
function run(obstacle, startX, startY, k, steps) {
  const start = gliderAt(startX + k, startY + k);
  let s = union(obstacle, start);
  const initMass = s.size;
  const history = [];
  for (let t = 0; t < steps; t++) {
    history.push(components(s).map(c => canon(c)));
    s = step(s);
  }
  const finalComps = history[history.length - 1];
  let movers = 0, stills = 0;
  for (const c of finalComps) {
    let moved = false;
    for (let t = history.length - 2; t >= 0 && !moved; t--) {
      for (const old of history[t]) {
        if (old.id === c.id && (Math.abs(old.cx - c.cx) > 1e-6 || Math.abs(old.cy - c.cy) > 1e-6)) { moved = true; break; }
      }
    }
    if (moved) movers++; else stills++;
  }
  return { initMass, finalMass: s.size, nComps: finalComps.length, movers, stills };
}

// ============================================================ CONTROLS
console.log('=== CONTROLS ===');
{
  // vacuum glider: same shape for every k
  let ok = true;
  const ids = [];
  for (let k = 0; k < 8; k++) {
    let s = gliderAt(30 + k, 30 + k);
    for (let t = 0; t < 60; t++) s = step(s);
    const cs = components(s);
    ids.push(cs.length === 1 ? canon(cs[0]).id : 'multi:' + cs.length);
  }
  console.log(`  vacuum glider over k=0..7: ${new Set(ids).size === 1 ? 'identical shape' : 'DIFFERS'}  ${new Set(ids).size === 1 ? 'PASS' : '*** FAIL ***'}`);
  if (new Set(ids).size !== 1) ok = false;

  // blinker must be period 2
  const b = blinkerAt(60, 60, true);
  const b1 = step(b), b2 = step(b1);
  console.log(`  blinker period: ${signature(b2) === signature(b) && signature(b1) !== signature(b) ? '2' : 'NOT 2'}  ${signature(b2) === signature(b) ? 'PASS' : '*** FAIL ***'}`);
  if (signature(b2) !== signature(b)) ok = false;
  console.log(`  blinker phases: ${signature(b) === signature(blinkerAt(60, 60, false)) ? '(vertical == horizontal after a step)' : ''}`);

  // block must be period 1
  const bl = blockAt(60, 60);
  console.log(`  block period: ${signature(step(bl)) === signature(bl) ? '1' : 'NOT 1'}  ${signature(step(bl)) === signature(bl) ? 'PASS' : '*** FAIL ***'}`);
  if (!ok) process.exit(1);
}
console.log('  controls passed\n');

// ============================================================ the two walls
function sweep(label, obstacle, startX, startY, steps) {
  console.log(`--- ${label} ---`);
  console.log('   k   finalMass  components  movers  stills');
  const rows = [];
  for (let k = 0; k < 12; k++) {
    const r = run(obstacle, startX, startY, k, steps);
    rows.push({ k, ...r });
    console.log(`  ${String(k).padStart(2)}  ${String(r.finalMass).padStart(9)}  ${String(r.nComps).padStart(10)}  ${String(r.movers).padStart(6)}  ${String(r.stills).padStart(6)}`);
  }
  const masses = rows.map(r => r.finalMass);
  const comps = rows.map(r => r.nComps);
  const movers = rows.map(r => r.movers);
  const stills = rows.map(r => r.stills);
  const varies = new Set(masses).size > 1 || new Set(comps).size > 1 ||
    new Set(movers).size > 1 || new Set(stills).size > 1;
  console.log(`  mass values: ${[...new Set(masses)].sort((a, b) => a - b).join(', ')}`);
  console.log(`  components : ${[...new Set(comps)].sort((a, b) => a - b).join(', ')}`);
  console.log(`  movers     : ${[...new Set(movers)].sort((a, b) => a - b).join(', ')}`);
  console.log(`  => OUTCOME ${varies ? 'VARIES with k' : 'IDENTICAL for all k'}`);
  console.log('');
  return varies;
}

console.log('=== EXPERIMENT: fixed obstacle, glider phase varied ===\n');

const blockVaries = sweep(
  'CONTROL WALL: a single BLOCK (period 1, no phase)',
  blockAt(62, 62), 32, 32, 220
);

const blinkerVaries = sweep(
  'DYNAMIC WALL: a single BLINKER (period 2, has a phase)',
  blinkerAt(62, 62, true), 32, 32, 220
);

const blinkerVaries2 = sweep(
  'DYNAMIC WALL: a single BLINKER, vertical phase',
  blinkerAt(62, 62, false), 32, 32, 220
);

console.log('=== VERDICTS ===');
console.log(`  block   (period 1): outcome ${blockVaries ? 'VARIES' : 'identical'}  <- predicted identical`);
console.log(`  blinker (period 2): outcome ${blinkerVaries ? 'VARIES' : 'identical'}`);
console.log(`  blinker, other phase: outcome ${blinkerVaries2 ? 'VARIES' : 'identical'}`);
console.log('');
if (!blockVaries && (blinkerVaries || blinkerVaries2)) {
  console.log('  => The user\'s mechanism is CONFIRMED in this setup: a period-1 obstacle');
  console.log('     gives phase-independent results, while a period-2 (dynamic) obstacle');
  console.log('     gives phase-DEPENDENT results. Identical-looking firings differ only');
  console.log('     by the hidden phase relation.');
} else if (!blockVaries && !blinkerVaries && !blinkerVaries2) {
  console.log('  => Neither wall shows phase dependence. The glider passes or is absorbed');
  console.log('     the same way regardless of the timing.');
} else if (blockVaries) {
  console.log('  => The CONTROL itself varies, which means the sweep is not isolating');
  console.log('     phase: something else changes with k. Treat all results as suspect.');
}
