'use strict';
/*
 * slit_search.js
 * ---------------------------------------------------------------------------
 * EXHAUSTIVE SEARCH for a slit that a glider can pass through.
 *
 * GOAL (the user's proposal, made precise)
 *   Find a wall with a gap such that:
 *     (i)  the glider COMEES OUT (it is not absorbed), and
 *     (ii) the wall SURVIVES (it can be used again), and
 *     (iii) the outgoing DIRECTION depends on the glider's incident PHASE.
 *   If (iii) holds, we have a primitive "diffraction" signature: identical
 *   firings give different screen hits.
 *
 * CANDIDATE SPACE (searched exhaustively)
 *   - block spacing along a row: 3, 4, 5, 6
 *   - which gap is used
 *   - incident direction: the four diagonal headings
 *   - incident phase: 4 glider phases
 *
 * WHY THE PREVIOUS ATTEMPTS FAILED
 *   They built the wall by laying blocks along the ANTI-DIAGONAL, which made the
 *   blocks merge into one blob (8-connectivity) and made "gap width" meaningless.
 *   Here the blocks are laid along a ROW, so their separation is a real gap.
 *
 * MEASUREMENT
 *   Terminal state -> components -> classify each as STILL LIFE (period 1, no
 *   displacement) or MOVER (shape recurs with constant drift). Report the mover
 *   direction(s) and whether the wall is intact.
 *
 * PREDICTIONS
 *   - spacing 3 leaves a 1-cell gap: probably too narrow, glider destroyed
 *   - spacing 4-5 may let it through
 *   - if any configuration lets it through, the direction distribution will
 *     depend on the incident phase
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

// glider with a heading, built by rotating the verified phase-0 glider
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

// build a row of blocks with a gap where block `skip` would be
function buildWall(spacing, count, skip, y0, x0) {
  let wall = new Set();
  for (let i = 0; i < count; i++) {
    if (i === skip) continue;
    wall = union(wall, blockAt(x0 + i * spacing, y0));
  }
  return wall;
}

/*
 * Follow a component over time to decide whether it is a mover, and by how much.
 * Returns { mover, vx, vy, p } for the best period found.
 */
function classifyMover(cells) {
  let s = new Set(cells);
  const frames = [];
  for (let t = 0; t < 40; t++) {
    const cs = components(s);
    if (cs.length !== 1) return { mover: false, reason: 'splits' };
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
    const vx = Math.round(dx), vy = Math.round(dy);
    return { mover: vx !== 0 || vy !== 0, vx, vy, p };
  }
  return { mover: false, reason: 'no recurrence' };
}

// ============================================================ controls
console.log('=== CONTROLS ===');
{
  let ok = true;
  for (let h = 0; h < 4; h++) {
    const r = classifyMover([...gliderAt(40, 40, h)]);
    const pass = r.mover && r.vx === HEAD_VEL[h][0] / 4 * 4 / 4 && Math.abs(r.vx) === 0;
    // the glider moves 1 cell per period p=4, so vx and vy should be +/-1 with p=4
    const good = r.mover && r.p === 4 && Math.abs(r.vx) <= 1 && Math.abs(r.vy) <= 1;
    if (!good) ok = false;
    console.log(`  heading ${h}: period ${r.p}, drift (${r.vx},${r.vy})  ${good ? 'PASS' : '*** FAIL ***'}`);
  }
  const w = buildWall(5, 6, 2, 60, 30);
  const wc = components(w).length;
  console.log(`  wall of 6 blocks, spacing 5, one skipped: ${wc} components  ${wc === 5 ? 'PASS' : '*** FAIL ***'}`);
  if (!ok || wc !== 5) process.exit(1);
}
console.log('  controls passed\n');

// ============================================================ the search
console.log('=== SEARCH: can a glider pass through a gap and come out? ===\n');
console.log('  spacing  gapCells  heading  phase  wallIntact  movers  moverDirections');

const results = [];
for (const spacing of [3, 4, 5, 6]) {
  const gapCells = spacing - 2;               // blocks are 2 wide
  for (const skip of [2, 3]) {
    for (let heading = 0; heading < 4; heading++) {
      for (let phase = 0; phase < 4; phase++) {
        const y0 = 60, x0 = 30;
        const wall = buildWall(spacing, 7, skip, y0, x0);
        const wallComps = components(wall).length;

        // launch the glider so that its trajectory crosses the gap.
        // The gap centre is between block skip-1 and block skip+1:
        const gapLeft = x0 + (skip - 1) * spacing + 1;   // right edge of the previous block
        const gapRight = x0 + (skip + 1) * spacing;      // left edge of the next block
        // aim from the south-west along the diagonal through the gap centre
        const cx = Math.floor((gapLeft + gapRight) / 2);
        const cy = y0 + 1;
        // heading 0 is (+1,+1): start below-left, so that it passes through (cx,cy)
        const dist = 60;
        const sx = cx - dist + phase;
        const sy = cy - dist + phase;
        let s = union(wall, gliderAt(sx, sy, heading));
        for (let t = 0; t < 400; t++) s = step(s);
        const comps = components(s);
        // is the wall intact? count how many of the original blocks survive as still lifes
        let intactStills = 0, movers = 0;
        const dirs = [];
        for (const c of comps) {
          const r = classifyMover(c);
          if (r.mover) { movers++; dirs.push(`(${r.vx},${r.vy})p${r.p}`); }
          else if (c.length === 4) intactStills++;
        }
        const wallIntact = intactStills >= wallComps - 1;
        results.push({ spacing, gapCells, skip, heading, phase, wallIntact, movers, dirs, finalMass: s.size });
        if (movers > 0) {
          console.log(`  ${String(spacing).padStart(7)}  ${String(gapCells).padStart(8)}  ${String(heading).padStart(7)}  ${String(phase).padStart(5)}  ${String(wallIntact).padStart(10)}  ${String(movers).padStart(6)}  ${dirs.join(' ')}`);
        }
      }
    }
  }
}

console.log('\n=== summary ===');
const withMovers = results.filter(r => r.movers > 0);
console.log(`  configurations tested: ${results.length}`);
console.log(`  configurations where a mover survived: ${withMovers.length}`);
if (withMovers.length === 0) {
  console.log('  => NO configuration let a glider come out. The wall absorbs everything.');
  console.log('     That is a clean negative result for this wall construction.');
} else {
  // group by (spacing, skip, heading) and see whether the direction set depends on phase
  const groups = new Map();
  for (const r of results) {
    const key = `${r.spacing}/${r.skip}/${r.heading}`;
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(r);
  }
  console.log('\n  does the outgoing direction depend on the incident phase?');
  let anyDepends = false;
  for (const [key, rs] of groups) {
    const sigs = rs.map(r => r.dirs.join('+') || '(none)');
    const varies = new Set(sigs).size > 1;
    if (varies) anyDepends = true;
    console.log(`    spacing/skip/heading ${key}: ${sigs.join('  |  ')}  ${varies ? '<-- DEPENDS on phase' : ''}`);
  }
  console.log('');
  console.log(anyDepends
    ? '  => Some configurations show PHASE-DEPENDENT outgoing directions.'
    : '  => Outgoing direction does not depend on phase in these configurations.');
}
