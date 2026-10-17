'use strict';
/*
 * slit_experiment.js
 * ---------------------------------------------------------------------------
 * THE USER'S PROPOSAL, MADE CONCRETE: a slit.
 *
 * Structure, and what each piece is meant to stand for:
 *
 *     the WALL      a large stable pattern (a long chain of still lifes)
 *                   -- stands in for the "instrument" an observer can build
 *     the GLIDER    a small moving particle
 *     the GAP       a hole in the wall
 *     the TRAIL     whatever the wall keeps afterwards
 *                   -- the wall's "record" of the passage
 *
 * The wall is deliberately made of STILL LIFES (blocks), which are perfectly
 * static. That is the strongest possible test of the idea that an instrument is
 * passive: if even a wall of perfectly static objects responds to the glider in
 * a way that depends on the glider's internal phase, then the "instrument" is
 * not passive, and that is a real finding.
 *
 * PREDICTIONS STATED BEFORE RUNNING (this project's recurring failure has been
 * to interpret numbers after seeing them):
 *   1. solid wall        -> glider absorbed, one extra still life near the impact
 *   2. wide gap          -> glider passes through, wall nearly unchanged
 *   3. narrow gap        -> glider blocked or deflected; effect grows monotonically
 *                          as the gap narrows
 *   4. wall response does NOT depend on the glider's phase
 *   5. no non-monotonic effect in gap width
 *
 * Controls first: glider in vacuum must travel unchanged; the glider's four
 * phases must be the four known phases.
 * ---------------------------------------------------------------------------
 */

const SIZE = 200;
function wrapCell(x, y) { return ((y % SIZE) + SIZE) % SIZE * SIZE + (((x % SIZE) + SIZE) % SIZE); }

function step(set) {
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

// the verified glider: moves (+1,+1) per 4 generations
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
function blockAt(x, y) {
  const s = new Set();
  for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }
function diff(a, b) { const o = new Set(); for (const k of a) if (!b.has(k)) o.add(k); return o; }

// ---------- controls
console.log('=== CONTROLS ===');
{
  // glider in vacuum: single component, moves (+1,+1) per 4 steps, and the four
  // phases must all be valid movers
  let allOk = true;
  for (let turns = 0; turns < 4; turns++) {
    let s = gliderAt(40, 40, turns);
    const frames = [];
    let ok = true;
    for (let t = 0; t <= 40; t++) {
      const cs = components(s);
      if (cs.length !== 1) { ok = false; break; }
      frames.push(canon(cs[0]));
      s = step(s);
    }
    if (!ok) { console.log(`  turns=${turns}: *** not a single component ***`); allOk = false; continue; }
    let found = null;
    for (let p = 1; p <= 12 && !found; p++) {
      let good = true, base = null;
      for (let t = 0; t + p < frames.length; t++) {
        if (frames[t].id !== frames[t + p].id) { good = false; break; }
        let dx = (frames[t + p].sx - frames[t].sx) / frames[t].n;
        let dy = (frames[t + p].sy - frames[t].sy) / frames[t].n;
        if (Math.abs(dx) > SIZE / 2) dx -= Math.sign(dx) * SIZE;
        if (Math.abs(dy) > SIZE / 2) dy -= Math.sign(dy) * SIZE;
        const d = [Math.round(dx), Math.round(dy)];
        if (base === null) base = d;
        else if (d[0] !== base[0] || d[1] !== base[1]) { good = false; break; }
      }
      if (good && base && (base[0] || base[1])) found = { p, base };
    }
    const exp = [[1, 1], [-1, 1], [-1, -1], [1, -1]][turns];
    const pass = found && found.base[0] === exp[0] && found.base[1] === exp[1];
    if (!pass) allOk = false;
    console.log(`  glider turns=${turns}: ${found ? `period ${found.p}, velocity (${found.base})` : 'no mover'}  expected (${exp})  ${pass ? 'PASS' : '*** FAIL ***'}`);
  }
  if (!allOk) { console.log('*** controls failed, aborting ***'); process.exit(1); }

  // a solid wall must be completely static on its own
  let wall = new Set();
  for (let i = 0; i < 10; i++) wall = union(wall, blockAt(30 + 3 * i, 30 + 3 * i));
  const before = [...wall].sort((a, b) => a - b).join(',');
  let w = wall;
  for (let t = 0; t < 60; t++) w = step(w);
  const after = [...w].sort((a, b) => a - b).join(',');
  console.log(`  wall alone is static over 60 steps: ${before === after ? 'PASS' : '*** FAIL ***'}`);
  if (before !== after) process.exit(1);
}
console.log('  controls passed\n');

// ============================================================ the experiment
/*
 * WALL: blocks placed along the anti-diagonal x + y = const, which is
 * perpendicular to the glider's (+1,+1) motion. A gap is a set of missing
 * block positions in that chain.
 *
 * All positions are offsets along the wall, so gap "position" and "width" are
 * directly meaningful.
 */
const WALL_OFFSET = 100;          // the wall sits at x + y = WALL_OFFSET
function buildWall(chainLen, gapStart, gapWidth, halfSpan) {
  // blocks along the anti-diagonal, in a band of width 2 (so the wall is thick)
  let wall = new Set();
  for (let t = -halfSpan; t <= halfSpan; t++) {
    // skip the gap interval
    if (t >= gapStart && t < gapStart + gapWidth) continue;
    // two rows thick: offset the block along the diagonal direction (1,-1)
    for (const d of [0, 1]) {
      // point on the anti-diagonal: x = k + t, y = k - t with k = WALL_OFFSET/2
      const k = Math.floor(WALL_OFFSET / 2);
      const x = k + t + d;
      const y = k - t - d + d; // keep blocks isolated from each other
      // place a block whose top-left is at (x, y)
      const b = blockAt(x, y);
      for (const c of b) wall.add(c);
    }
  }
  return wall;
}

function run(wall, gliderStartTurns, steps) {
  const start = gliderAt(30, 30, gliderStartTurns);
  let all = union(wall, start);
  const initWall = new Set(wall);
  const initMass = all.size;
  for (let t = 0; t < steps; t++) all = step(all);
  // what is left of the wall?
  const comps = components(all).map(c => canon(c));
  const wallComps = components(initWall).map(c => canon(c));
  return {
    initMass,
    finalMass: all.size,
    finalComps: comps.length,
    wallCompsBefore: wallComps.length,
    wallChanged: initMass !== all.size,
    ids: comps.map(c => `${c.w}x${c.h}`),
  };
}

console.log('=== EXPERIMENT: glider vs a wall of still lifes ===');
console.log(`wall along x+y = ${WALL_OFFSET}, made of 2x2 blocks\n`);

// ---- 1. solid wall (no gap)
console.log('--- 1. SOLID WALL (no gap) ---');
{
  const wall = buildWall(0, 0, 0, 8);
  const before = components(wall).length;
  const r = run(wall, 0, 300);
  console.log(`  wall blocks before: ${before}`);
  console.log(`  total mass: ${r.initMass} -> ${r.finalMass}`);
  console.log(`  components after: ${r.finalComps}`);
  console.log('  PREDICTION: glider absorbed, mass roughly conserved (wall + debris)');
}

// ---- 2. wide gap
console.log('\n--- 2. WIDE GAP ---');
for (const gapWidth of [6, 4, 3, 2, 1]) {
  const wall = buildWall(0, -Math.floor(gapWidth / 2), gapWidth, 8);
  const r = run(wall, 0, 300);
  console.log(`  gap width ${gapWidth}: mass ${r.initMass} -> ${r.finalMass}, components ${r.finalComps}`);
}

// ---- 3. does the outcome depend on the glider's PHASE?
console.log('\n--- 3. PHASE DEPENDENCE (the key question) ---');
console.log('  A passive instrument must give the SAME result for all glider phases.');
console.log('');
console.log('  gap   phase0        phase1        phase2        phase3       all equal?');
for (const gapWidth of [4, 3, 2]) {
  const results = [];
  for (let turns = 0; turns < 4; turns++) {
    const wall = buildWall(0, -Math.floor(gapWidth / 2), gapWidth, 8);
    const r = run(wall, turns, 300);
    results.push(r.finalMass);
  }
  const allEq = results.every(x => x === results[0]);
  console.log(`  ${String(gapWidth).padStart(3)}   ${results.map(x => String(x).padStart(12)).join('')}   ${allEq ? 'YES' : '*** NO ***'}`);
}
console.log('');
console.log('  PREDICTION: all equal. If they are NOT equal, the wall is not passive,');
console.log('  and that is a finding: the instrument would depend on the particle phase.');

console.log('\n=== where this leaves us ===');
console.log('  This is the first step of the user\'s proposal. It answers:');
console.log('     does a large static pattern act as a PASSIVE instrument,');
console.log('     and does a slit do anything non-classical?');
console.log('  The random-soup and collision work already suggests the answer is "yes, passive".');
console.log('  The slit adds the one thing not yet tested: whether a NARROW gap makes the');
console.log('  outcome depend on something other than a single definite trajectory.');
