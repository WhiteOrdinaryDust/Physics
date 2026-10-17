'use strict';
/*
 * hidden_phase_test2.js
 * ---------------------------------------------------------------------------
 * THE CORRECTED TEST.
 *
 * WHAT WAS WRONG IN THE PREVIOUS ATTEMPT
 *   It shifted the wall AND the glider together by k cells. But Life is
 *   translation invariant, so shifting the entire apparatus is just translating
 *   the same experiment. The outcome was necessarily identical, and the test
 *   could not have detected anything. Its "negative result" is void.
 *
 * THE CORRECT TEST
 *   Keep the WALL FIXED and move the GLIDER's starting position by k cells along
 *   its own direction of travel. Because the glider has period 4 with a
 *   displacement of 1 cell per period, the value of k modulo 4 IS the glider's
 *   internal phase at any later moment. So this isolates exactly the hidden
 *   parameter the user is pointing at.
 *
 *   Both the wall and the target point stay fixed, so any change in outcome is
 *   attributable to the glider's phase and nothing else.
 *
 * CONTROLS
 *   1. glider in vacuum, started at k = 0..7: outcome must be identical
 *      (translation invariance again, but now applied to the glider alone)
 *   2. wall alone: static
 *
 * SCENARIOS
 *   A. single block in the glider's path      (minimal obstacle)
 *   B. a short wall with a gap                (the slit-like case)
 *
 * PREDICTIONS, STATED BEFORE RUNNING
 *   - control 1 and 2 pass
 *   - scenario A: the outcome MAY vary with k
 *   - scenario B: ditto
 *   - if neither varies, then for these obstacles the hidden phase is irrelevant
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

/*
 * Run one experiment: a fixed obstacle set plus a glider started at
 * (startX + k, startY + k) whose direction is always (+1,+1).
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
  const moverIds = [];
  for (const c of finalComps) {
    let moved = false;
    for (let t = history.length - 2; t >= 0 && !moved; t--) {
      for (const old of history[t]) {
        if (old.id === c.id && (Math.abs(old.cx - c.cx) > 1e-6 || Math.abs(old.cy - c.cy) > 1e-6)) { moved = true; break; }
      }
    }
    if (moved) { movers++; moverIds.push(c.id.split(':')[0]); }
    else stills++;
  }
  return { initMass, finalMass: s.size, nComps: finalComps.length, movers, stills, moverIds: moverIds.sort() };
}

// ============================================================ CONTROLS
console.log('=== CONTROLS ===');
{
  console.log('--- control 1: glider alone, started at k ---');
  const ids = [];
  for (let k = 0; k < 8; k++) {
    let s = gliderAt(30 + k, 30 + k);
    for (let t = 0; t < 80; t++) s = step(s);
    ids.push(signature(step(s)) === signature(s) ? 'dead' : 'alive:' + components(s).map(c => canon(c).id.split(':')[0]).join('+'));
  }
  console.log(`  outcomes: ${ids.join('  ')}`);
  console.log(`  all identical? ${new Set(ids).size === 1 ? 'PASS' : '*** FAIL ***'}`);
  if (new Set(ids).size !== 1) process.exit(1);

  console.log('--- control 2: obstacle alone must be static ---');
  const b = blockAt(60, 60);
  const before = signature(b);
  let s = b;
  for (let t = 0; t < 40; t++) s = step(s);
  console.log(`  single block static? ${signature(s) === before ? 'PASS' : '*** FAIL ***'}`);
  if (signature(s) !== before) process.exit(1);
}
console.log('  controls passed\n');

// ============================================================ SCENARIO A
console.log('=== SCENARIO A: ONE BLOCK in the glider\'s path ===');
console.log('  The block is FIXED. Only the glider start is moved by k.\n');
{
  const obstacle = blockAt(60, 60);
  // the glider is launched from the southwest; it reaches x=60 around t = 4*(60-32) = 112
  const startX = 32, startY = 32;
  console.log('   k    initMass  finalMass  components  movers  stills   mover shapes');
  const rows = [];
  for (let k = 0; k < 12; k++) {
    const r = run(obstacle, startX, startY, k, 200);
    rows.push({ k, ...r });
    console.log(`  ${String(k).padStart(2)}   ${String(r.initMass).padStart(8)}  ${String(r.finalMass).padStart(9)}  ${String(r.nComps).padStart(10)}  ${String(r.movers).padStart(6)}  ${String(r.stills).padStart(6)}   ${r.moverIds.join(' ')}`);
  }
  console.log('');
  const masses = rows.map(r => r.finalMass);
  const comps = rows.map(r => r.nComps);
  const movers = rows.map(r => r.movers);
  console.log(`  final mass varies? ${new Set(masses).size > 1 ? 'YES  -> ' + [...new Set(masses)].sort((a,b)=>a-b).join(', ') : 'no'}`);
  console.log(`  components vary?  ${new Set(comps).size > 1 ? 'YES  -> ' + [...new Set(comps)].sort((a,b)=>a-b).join(', ') : 'no'}`);
  console.log(`  movers vary?      ${new Set(movers).size > 1 ? 'YES  -> ' + [...new Set(movers)].sort((a,b)=>a-b).join(', ') : 'no'}`);
  console.log('  periodicity in k:');
  for (const p of [1, 2, 3, 4, 6]) {
    let per = true;
    for (let i = 0; i + p < rows.length; i++) if (masses[i] !== masses[i + p]) { per = false; break; }
    console.log(`    period ${p}: ${per ? 'CONSISTENT' : 'no'}`);
  }
  const deps = new Set(masses).size > 1 || new Set(comps).size > 1 || new Set(movers).size > 1;
  console.log(`\n  VERDICT: hidden phase ${deps ? 'MATTERS' : 'does NOT matter'} for a single block.`);
}

// ============================================================ SCENARIO B
console.log('\n=== SCENARIO B: a FIXED wall with a gap ===');
{
  const wallGap = 6;         // skip block index 6
  const buildWall = () => {
    let wall = new Set();
    for (let i = 0; i < 12; i++) {
      if (i === wallGap) continue;
      wall = union(wall, blockAt(50 + 3 * i, 70));
    }
    return wall;
  };
  const obstacle = buildWall();
  console.log(`  wall: ${obstacle.size} cells, ${components(obstacle).length} components (blocks)`);
  console.log('   k    finalMass  components  movers  stills');
  const rows = [];
  for (let k = 0; k < 12; k++) {
    const r = run(obstacle, 28, 48, k, 220);
    rows.push({ k, ...r });
    console.log(`  ${String(k).padStart(2)}   ${String(r.finalMass).padStart(9)}  ${String(r.nComps).padStart(10)}  ${String(r.movers).padStart(6)}  ${String(r.stills).padStart(6)}`);
  }
  console.log('');
  const masses = rows.map(r => r.finalMass);
  const comps = rows.map(r => r.nComps);
  const movers = rows.map(r => r.movers);
  console.log(`  final mass varies? ${new Set(masses).size > 1 ? 'YES  -> ' + [...new Set(masses)].sort((a,b)=>a-b).join(', ') : 'no'}`);
  console.log(`  components vary?  ${new Set(comps).size > 1 ? 'YES  -> ' + [...new Set(comps)].sort((a,b)=>a-b).join(', ') : 'no'}`);
  console.log(`  movers vary?      ${new Set(movers).size > 1 ? 'YES  -> ' + [...new Set(movers)].sort((a,b)=>a-b).join(', ') : 'no'}`);
  const deps = new Set(masses).size > 1 || new Set(comps).size > 1 || new Set(movers).size > 1;
  console.log(`\n  VERDICT: hidden phase ${deps ? 'MATTERS' : 'does NOT matter'} for this wall.`);
}

console.log('\n=== summary ===');
console.log('  The question was: can two firings that LOOK the same give different');
console.log('  results because the glider phase differs? Read the two verdicts above.');
console.log('  Note that a single block is a PERFECTLY STATIC obstacle, so if the phase');
console.log('  matters even there, the effect cannot come from the obstacle\'s own');
console.log('  dynamics -- it must come from the glider\'s phase relative to the lattice.');
