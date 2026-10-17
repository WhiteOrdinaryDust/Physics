'use strict';
/*
 * hidden_phase_test.js
 * ---------------------------------------------------------------------------
 * THE USER'S PROPOSAL, AND A CLEAN WAY TO TEST IT.
 *
 * The idea: two firings that LOOK identical may have different inputs, because
 * there are hidden parameters:
 *     - the glider's internal PHASE (4 distinct shapes, all called "a glider")
 *     - the wall's own internal phase (it is a dynamic stable pattern too)
 *     - the background
 * So "the same experiment" can give different results. It is a combinatorial
 * matter of which phases line up.
 *
 * WHY A RIGID TRANSLATION IS THE CLEAN TEST
 *   In vacuum a glider moves (+1,+1) per 4 steps, and Life is translation
 *   invariant. So shifting the WHOLE setup (wall + glider together) by k cells
 *   changes nothing about the glider's own dynamics. What it DOES change is the
 *   glider's position modulo 4 relative to the wall, i.e. exactly the hidden
 *   relative phase.
 *
 *   Therefore: if the outcome depends on k, the dependence is not on geometry
 *   (geometry is identical for every k) but on the hidden phase. That is the
 *   user's mechanism, isolated and measurable.
 *
 * PRIOR EVIDENCE
 *   An earlier, badly-built run (slit_experiment.js) already hinted at phase
 *   dependence, but its wall geometry was broken (17 blocks had merged into one
 *   blob), so its numbers were void. This file rebuilds the wall correctly.
 *
 * WALL CONSTRUCTION (verified in slit_experiment2.js)
 *   2x2 blocks spaced 3 apart along a line stay SEPARATE and are perfectly
 *   static. That is the only wall shape we have verified.
 *
 * PREDICTIONS, STATED BEFORE RUNNING
 *   1. glider in vacuum: outcome identical for every shift k  (control)
 *   2. wall alone: static for every shift k                  (control)
 *   3. glider + wall: outcome VARIES with k
 *   4. the variation is PERIODIC in k with period 4 (the glider's period)
 * ---------------------------------------------------------------------------
 */

const SIZE = 140;
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
 * Baseline wall: 12 blocks spaced 3 apart along a horizontal line at y = 70,
 * x = 40, 43, 46, ...  With the gap made by SKIPPING block index 6.
 * Verified separate-and-static in slit_experiment2.js.
 */
function buildWall(shift, skipFrom, skipCount) {
  let wall = new Set();
  for (let i = 0; i < 12; i++) {
    if (i >= skipFrom && i < skipFrom + skipCount) continue;
    const x = 40 + 3 * i + shift;
    const y = 70 + shift;
    wall = union(wall, blockAt(x, y));
  }
  return wall;
}

// classify the terminal state: count components, find movers and their velocity
function classify(set, steps) {
  // record component signatures over time to detect movers
  let s = set;
  const history = [];
  for (let t = 0; t < steps; t++) {
    history.push(components(s).map(c => canon(c)));
    s = step(s);
  }
  const finalComps = history[history.length - 1];
  let movers = 0, stills = 0, others = 0;
  const moverIds = [];
  for (const c of finalComps) {
    // find whether this shape appeared earlier at a DIFFERENT position
    let moved = false;
    for (let t = history.length - 2; t >= 0 && !moved; t--) {
      for (const old of history[t]) {
        if (old.id === c.id && (Math.abs(old.cx - c.cx) > 1e-6 || Math.abs(old.cy - c.cy) > 1e-6)) {
          moved = true; break;
        }
      }
    }
    if (moved) { movers++; moverIds.push(c.id); }
    else stills++;
  }
  return { nComps: finalComps.length, movers, stills, mass: set.size };
}

// ============================================================ CONTROLS
console.log('=== CONTROLS ===');
{
  // control 1: glider alone, shifted by k -> identical outcome for every k
  console.log('--- control 1: glider in vacuum, shifted by k ---');
  let ok = true;
  const base = [];
  for (let k = 0; k < 8; k++) {
    let s = gliderAt(30 + k, 30 + k);
    for (let t = 0; t < 60; t++) s = step(s);
    const cs = components(s);
    base.push({ k, n: cs.length, mass: s.size, id: cs.length === 1 ? canon(cs[0]).id : 'multi' });
  }
  const ids = new Set(base.map(b => b.id));
  console.log(`  outcomes for k = 0..7: ${base.map(b => b.id.split(':')[0]).join(' ')}`);
  console.log(`  all identical? ${ids.size === 1 ? 'PASS' : '*** FAIL ***'}`);
  if (ids.size !== 1) ok = false;

  // control 2: wall alone, shifted, must be static
  console.log('--- control 2: wall alone, shifted by k, must be static ---');
  let ok2 = true;
  for (let k = 0; k < 4; k++) {
    const wall = buildWall(k, 6, 1);
    const before = signature(wall);
    let s = wall;
    for (let t = 0; t < 40; t++) s = step(s);
    const same = signature(s) === before;
    if (!same) ok2 = false;
    console.log(`  shift ${k}: ${wall.size} cells -> static? ${same ? 'YES' : 'NO'}`);
  }
  if (!ok2) ok = false;
  if (!ok) { console.log('*** a control failed, aborting ***'); process.exit(1); }
}
console.log('  controls passed\n');

// ============================================================ THE TEST
console.log('=== TEST: shift the WHOLE setup by k, and see whether the outcome changes ===');
console.log('  The glider is aimed at the gap. Geometry is identical for every k.');
console.log('  Only the glider-phase-relative-to-wall changes.\n');
console.log('   k    components   movers   stills   mass');
const rows = [];
for (let k = 0; k < 12; k++) {
  const wall = buildWall(k, 6, 1);
  // glider aimed at the gap: gap is where block index 6 would be, i.e. x = 40+18 = 58
  const gx = 58 + k - 30, gy = 70 + k - 30;
  const start = gliderAt(gx, gy);
  const all = union(wall, start);
  const initMass = all.size;
  let s = all;
  for (let t = 0; t < 120; t++) s = step(s);
  const r = classify(s, 40);
  rows.push({ k, ...r, initMass });
  console.log(`  ${String(k).padStart(2)}   ${String(r.nComps).padStart(10)}   ${String(r.movers).padStart(6)}   ${String(r.stills).padStart(6)}   ${String(r.mass).padStart(5)}`);
}

console.log('\n=== analysis ===');
const masses = rows.map(r => r.mass);
const comps = rows.map(r => r.nComps);
const movers = rows.map(r => r.movers);
const allSameMass = new Set(masses).size === 1;
const allSameComps = new Set(comps).size === 1;
const allSameMovers = new Set(movers).size === 1;
console.log(`  final mass   : ${allSameMass ? 'IDENTICAL for all k' : `VARIES: ${[...new Set(masses)].sort((a,b)=>a-b).join(', ')}`}`);
console.log(`  components   : ${allSameComps ? 'identical' : `varies: ${[...new Set(comps)].sort((a,b)=>a-b).join(', ')}`}`);
console.log(`  movers       : ${allSameMovers ? 'identical' : `varies: ${[...new Set(movers)].sort((a,b)=>a-b).join(', ')}`}`);

// periodicity in k
console.log('\n  is the pattern periodic in k?');
for (const p of [1, 2, 3, 4, 6]) {
  let per = true;
  for (let i = 0; i + p < rows.length; i++) if (masses[i] !== masses[i + p]) { per = false; break; }
  console.log(`    period ${p}: ${per ? 'CONSISTENT' : 'no'}`);
}

console.log('\n=== conclusion ===');
if (!allSameMass || !allSameComps || !allSameMovers) {
  console.log('  The outcome DEPENDS on the shift k, even though the geometry is identical.');
  console.log('  => The user\'s mechanism is confirmed: two firings that look the same are');
  console.log('     not the same, because the hidden relative phase differs.');
} else {
  console.log('  The outcome is IDENTICAL for every shift k.');
  console.log('  => For this wall and this glider, the hidden phase does NOT matter.');
  console.log('     That would be a negative result for the mechanism in this setup.');
}
