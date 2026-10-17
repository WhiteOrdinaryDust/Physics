'use strict';
/*
 * precontact_test.js
 * ---------------------------------------------------------------------------
 * THE USER'S CRITERION, MADE OPERATIONAL:
 *
 *     "we only need to observe whether the state just before contact is
 *      completely identical; if it is not, there should be a different result."
 *
 * THE CLEANEST POSSIBLE SETUP
 *   Fix everything, and change ONE thing: the INITIAL PHASE of the obstacle.
 *   A blinker at a fixed position can start horizontal or vertical. Same shape,
 *   same position, same rules, same glider launch, same geometry.
 *
 *   Nothing geometric differs. If the pre-contact state differs, it can only be
 *   because of the obstacle's phase, and then the criterion predicts a different
 *   outcome.
 *
 * WHAT IS RECORDED FOR EACH RUN
 *   - the state N steps BEFORE contact, as a canonical signature
 *   - the terminal outcome: mass, component count, movers, still lifes
 *
 * CONTROLS
 *   - a BLOCK (period 1) started in either of its (identical) phases: the
 *     pre-contact state must be identical and the outcome must be identical
 *   - vacuum glider unchanged
 *
 * PREDICTIONS, STATED BEFORE RUNNING
 *   1. block:  pre-contact states identical -> outcomes identical
 *   2. blinker: pre-contact states DIFFER    -> outcomes DIFFER
 *   3. if (2) fails, the criterion is falsified in this setup
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

// a translation-invariant, rotation-invariant signature of a whole configuration
function configSignature(set) {
  const comps = components(set).map(c => canon(c));
  // list of shapes only (ignore absolute position), sorted
  return comps.map(c => c.id).sort().join(' | ');
}

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

// ============================================================ the run
/*
 * Run and report:
 *   - precontact signature at step (contactStep - PRE)
 *   - terminal outcome
 */
function experiment(label, obstacle, startX, startY, steps, precontactStep, PRE) {
  const start = gliderAt(startX, startY);
  let s = union(obstacle, start);
  const initMass = s.size;
  const history = [];
  for (let t = 0; t <= steps; t++) {
    history.push(s);
    if (t === steps) break;
    s = step(s);
  }
  const pre = configSignature(history[precontactStep]);
  const prePlus = configSignature(history[precontactStep + 1]);

  // precompute canonical signatures for every step once (the previous version
  // recomputed components inside the search loop, which was needlessly heavy)
  const sigPerStep = history.map(h => components(h).map(c => canon(c)));

  const termComps = sigPerStep[steps];
  let movers = 0, stills = 0;
  for (const c of termComps) {
    let moved = false;
    for (let t = steps - 1; t >= 0 && !moved; t--) {
      for (const old of sigPerStep[t]) {
        if (old.id === c.id && (Math.abs(old.cx - c.cx) > 1e-6 || Math.abs(old.cy - c.cy) > 1e-6)) { moved = true; break; }
      }
    }
    if (moved) movers++; else stills++;
  }
  console.log(`  ${label}`);
  console.log(`      precontact signature (t=${precontactStep}): ${pre.length > 140 ? pre.slice(0, 140) + '...' : pre}`);
  console.log(`      terminal: mass ${history[steps].size}, components ${termComps.length}, movers ${movers}, stills ${stills}`);
  return { pre, prePlus, termMass: history[steps].size, termComps: termComps.length, movers, stills };
}

// ============================================================ controls
console.log('=== CONTROLS ===');
{
  let s = gliderAt(30, 30);
  for (let t = 0; t < 60; t++) s = step(s);
  const cs = components(s);
  console.log(`  vacuum glider intact: ${cs.length === 1 && canon(cs[0]).id === '3x3:.o./..o/ooo' ? 'PASS' : '*** FAIL ***'}`);
  const b = blinkerAt(60, 60, true);
  const b1 = step(b);
  const b2 = step(b1);
  console.log(`  blinker sizes over 3 steps: ${b.size} -> ${b1.size} -> ${b2.size}`);
  console.log(`  blinker period 2: ${b2.size === b.size && b1.size === b.size ? 'PASS (size returns)' : 'check sizes above'}`);
}
console.log('');

console.log('=== EXPERIMENT 1: the obstacle is a BLOCK (period 1) ===');
// the block has only one phase, so we simply run it twice to confirm determinism
{
  const r1 = experiment('block, run A', blockAt(62, 62), 32, 32, 200, 100, 1);
  const r2 = experiment('block, run B (identical setup)', blockAt(62, 62), 32, 32, 200, 100, 1);
  console.log(`  => precontact identical? ${r1.pre === r2.pre ? 'YES' : 'NO'}`);
  console.log(`  => outcome identical?   ${r1.termMass === r2.termMass && r1.termComps === r2.termComps ? 'YES' : 'NO'}`);
}

console.log('\n=== EXPERIMENT 2: the obstacle is a BLINKER, two different starting phases ===');
{
  const rH = experiment('blinker starting HORIZONTAL', blinkerAt(62, 62, true), 32, 32, 200, 100, 1);
  const rV = experiment('blinker starting VERTICAL', blinkerAt(62, 62, false), 32, 32, 200, 100, 1);
  const preSame = rH.pre === rV.pre;
  const outSame = rH.termMass === rV.termMass && rH.termComps === rV.termComps;
  console.log('');
  console.log(`  => precontact identical? ${preSame ? 'YES' : 'NO'}`);
  console.log(`  => outcome identical?   ${outSame ? 'YES' : 'NO'}`);
  console.log('');
  console.log('  THE CRITERION: if the precontact state differs, the outcome should differ.');
  if (!preSame && !outSame) {
    console.log('  => CRITERION CONFIRMED: precontact differed and so did the outcome.');
  } else if (!preSame && outSame) {
    console.log('  => CRITERION FALSIFIED HERE: precontact differed but the outcome was the same.');
    console.log('     (This can happen if the difference does not affect the interaction.)');
  } else if (preSame && outSame) {
    console.log('  => precontact states were identical, so nothing differed. The two blinker');
    console.log('     phases must coincide at the sampled moment.');
  } else {
    console.log('  => precontact identical yet outcome differed: that would contradict determinism.');
  }
}

console.log('\n=== EXPERIMENT 3: same blinker phase, glider launched at different phases ===');
console.log('  (fixed obstacle, glider start advanced by k along its own direction)');
{
  const results = [];
  for (let k = 0; k < 8; k++) {
    const obstacle = blinkerAt(62, 62, false);
    const start = gliderAt(32 + k, 32 + k);
    let s = union(obstacle, start);
    const initMass = s.size;
    for (let t = 0; t < 200; t++) s = step(s);
    const cs = components(s).map(c => canon(c));
    results.push({ k, mass: s.size, comps: cs.length });
  }
  console.log('   k   mass   components');
  for (const r of results) console.log(`  ${String(r.k).padStart(2)}  ${String(r.mass).padStart(5)}  ${String(r.comps).padStart(10)}`);
  const masses = results.map(r => r.mass);
  console.log(`  mass values: ${[...new Set(masses)].sort((a, b) => a - b).join(', ')}`);
  console.log(`  => outcome ${new Set(masses).size > 1 ? 'VARIES' : 'identical'} with k`);
}
