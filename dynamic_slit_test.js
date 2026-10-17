'use strict';
/*
 * dynamic_slit_test.js
 * ---------------------------------------------------------------------------
 * A DYNAMIC SLIT: gap edges made of period-2 oscillators.
 *
 * WHY
 *   The previous search used a STATIC slit and found that the outgoing direction
 *   never depended on the launch offset. That is forced: shifting a glider along
 *   its own direction of travel is just a time translation, and a static wall
 *   gives the same trajectory. To make the outgoing direction depend on timing,
 *   the slit itself must be dynamic.
 *
 *   We already established (dynamic_wall_test.js) that a period-2 obstacle gives
 *   phase-dependent OUTCOMES while a period-1 obstacle does not. This file asks
 *   the sharper question: does the OUTGOING DIRECTION change?
 *
 * SETUP
 *   A row of blocks as before, but the two blocks flanking the gap are replaced
 *   by blinkers (period 2). Everything else is as in slit_search.js.
 *   The glider is launched along the diagonal; its launch offset k (in cells,
 *   along its own direction) advances the arrival time by 4k steps, which after
 *   mod 2 changes which blinker phase it meets.
 *
 * MEASUREMENT
 *   After the run, find the MOVING component and measure its displacement over
 *   one period, giving the outgoing direction among the four diagonals. Also
 *   record the total number of movers (a secondary "screen" reading).
 *
 * CONTROLS
 *   - the four glider headings must be classified correctly in vacuum
 *   - a blinker must have period 2
 *
 * PREDICTIONS, STATED BEFORE RUNNING
 *   - static slit (blocks): outgoing direction independent of k
 *   - dynamic slit (blinkers): outgoing direction MAY depend on k
 * ---------------------------------------------------------------------------
 */

const SIZE = 200;
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
function blinkerAt(cx, cy, horizontal) {
  const s = new Set();
  for (let d = -1; d <= 1; d++) {
    if (horizontal) s.add(wrapCell(cx + d, cy));
    else s.add(wrapCell(cx, cy + d));
  }
  return s;
}

/*
 * Follow a single component in isolation and report whether it moves and how.
 */
function mobility(seed) {
  let s = new Set(seed);
  const frames = [];
  for (let t = 0; t < 40; t++) {
    const cs = components(s);
    if (cs.length !== 1) return { mover: false, why: 'not single' };
    frames.push(canon(cs[0]));
    s = step(s);
  }
  for (let p = 1; p <= 12; p++) {
    let ok = true;
    for (let t = 0; t + p < frames.length; t++) if (frames[t].id !== frames[t + p].id) { ok = false; break; }
    if (!ok) continue;
    let dx = (frames[p].sx - frames[0].sx) / frames[0].n;
    let dy = (frames[p].sy - frames[0].sy) / frames[0].n;
    if (Math.abs(dx) > SIZE / 2) dx -= Math.sign(dx) * SIZE;
    if (Math.abs(dy) > SIZE / 2) dy -= Math.sign(dy) * SIZE;
    return { mover: Math.round(dx) !== 0 || Math.round(dy) !== 0, vx: Math.round(dx), vy: Math.round(dy), p };
  }
  return { mover: false, why: 'no recurrence' };
}

/*
 * Build a wall: a row of blocks, with the two blocks flanking the gap optionally
 * replaced by blinkers.
 */
function buildSlit(spacing, count, skip, x0, y0, dynamicEdges, edgePhase) {
  let wall = new Set();
  for (let i = 0; i < count; i++) {
    if (i === skip) continue;
    const bx = x0 + i * spacing;
    const isEdge = dynamicEdges && (i === skip - 1 || i === skip + 1);
    if (isEdge) {
      wall = union(wall, blinkerAt(bx, y0, edgePhase));
    } else {
      wall = union(wall, blockAt(bx, y0));
    }
  }
  return wall;
}

// ============================================================ controls
console.log('=== CONTROLS ===');
{
  let ok = true;
  const exp = [[1, 1], [-1, 1], [-1, -1], [1, -1]];
  for (let h = 0; h < 4; h++) {
    // rotate by hand for the control
    const s = new Set();
    for (const [dx0, dy0] of GLIDER0) {
      let a = dx0, b = dy0;
      for (let i = 0; i < h; i++) { const na = -b, nb = a; a = na; b = nb; }
      s.add(wrapCell(60 + a, 60 + b));
    }
    const m = mobility([...s]);
    const good = m.mover && m.p === 4;
    if (!good) ok = false;
    console.log(`  heading ${h}: period ${m.p}, drift (${m.vx},${m.vy})  ${good ? 'PASS' : '*** FAIL ***'}`);
  }
  const b = blinkerAt(60, 60, false);
  const b2 = step(step(b));
  console.log(`  blinker period 2: ${b2.size === b.size && step(b).size === b.size ? 'PASS' : '*** FAIL ***'}`);
  const w = buildSlit(5, 7, 3, 30, 60, false, true);
  console.log(`  static slit wall: ${components(w).length} components (expected 6)`);
  if (!ok) process.exit(1);
}
console.log('  controls passed\n');

// ============================================================ the experiment
/*
 * Fire the glider so that its trajectory goes through the gap, and vary k (the
 * launch offset along its own direction), which changes the arrival time.
 */
function runSlit(dynamicEdges, edgePhase, k, steps) {
  const spacing = 5, count = 7, skip = 3, x0 = 30, y0 = 60;
  const wall = buildSlit(spacing, count, skip, x0, y0, dynamicEdges, edgePhase);
  // gap between block (skip-1)=2 and block (skip+1)=4
  const gapLeft = x0 + 2 * spacing + 1;     // right edge of block 2 (block is 2 wide)
  const gapRight = x0 + 4 * spacing;        // left edge of block 4
  const cx = Math.floor((gapLeft + gapRight) / 2);
  const cy = y0 + 1;
  const dist = 70;
  const start = gliderAt(cx - dist + k, cy - dist + k);
  let s = union(wall, start);
  for (let t = 0; t < steps; t++) s = step(s);
  const comps = components(s);
  const movers = [];
  let stills = 0;
  for (const c of comps) {
    const m = mobility(c);
    if (m.mover) movers.push(`(${m.vx},${m.vy})`);
    else stills++;
  }
  return { finalMass: s.size, nComps: comps.length, movers, stills };
}

console.log('=== STATIC slit (control): does the outgoing direction depend on k? ===');
console.log('   k   finalMass  comps  movers            stills');
{
  const rows = [];
  for (let k = 0; k < 12; k++) {
    const r = runSlit(false, true, k, 320);
    rows.push(r);
    console.log(`  ${String(k).padStart(2)}   ${String(r.finalMass).padStart(9)}  ${String(r.nComps).padStart(5)}  ${r.movers.join(' ').padEnd(18)} ${r.stills}`);
  }
  const sig = rows.map(r => r.movers.join('+'));
  console.log(`  mover signatures: ${[...new Set(sig)].join('  |  ')}`);
  console.log(`  => ${new Set(sig).size === 1 ? 'INDEPENDENT of k' : 'DEPENDS on k'}`);
}

console.log('\n=== DYNAMIC slit (blinker edges): does the outgoing direction depend on k? ===');
console.log('   k   finalMass  comps  movers            stills');
{
  const rows = [];
  for (let k = 0; k < 12; k++) {
    const r = runSlit(true, true, k, 320);
    rows.push(r);
    console.log(`  ${String(k).padStart(2)}   ${String(r.finalMass).padStart(9)}  ${String(r.nComps).padStart(5)}  ${r.movers.join(' ').padEnd(18)} ${r.stills}`);
  }
  const sig = rows.map(r => r.movers.join('+'));
  console.log(`  mover signatures: ${[...new Set(sig)].join('  |  ')}`);
  console.log(`  => ${new Set(sig).size === 1 ? 'INDEPENDENT of k' : 'DEPENDS on k'}`);
}

console.log('\n=== summary of the screen reading ===');
console.log('  For an ideal slit experiment we would want: identical inputs, and an');
console.log('  outgoing DIRECTION that depends on the timing. Read the two verdicts.');
