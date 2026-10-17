'use strict';
/*
 * pair_coupling.js
 * ---------------------------------------------------------------------------
 * STEP 10: separate "two particles" from "two CORRELATED particles", and measure
 * how sharply the response depends on arrival time.
 *
 * Previous round found:
 *   - ONE injected glider  : no effect at any injection time (final mass always
 *                            exactly 155 vs baseline 161)  -> zero coupling
 *   - a full ISLAND (a pair of gliders aimed at each other): peak ranges from 38
 *                            to 1229 depending ONLY on arrival time
 *   But the two cases had different initial mass (11 vs 20 cells), so the
 *   comparison was not clean. That is fixed here with three arms:
 *
 *   ARM A : one glider                     (11 cells)
 *   ARM B : two gliders flying IN PARALLEL (20 cells, no mutual collision)
 *   ARM C : two gliders aimed AT EACH OTHER(20 cells, mutually colliding)
 *
 *   A vs B isolates "is there a second particle at all"
 *   B vs C isolates "is the pair internally CORRELATED (on a collision course)"
 *
 * The B-vs-C comparison is the one that matters: if B behaves like A (no effect)
 * and only C couples, then what the cascade responds to is not the presence of
 * particles but the RELATION between them.
 * ---------------------------------------------------------------------------
 */

const SIZE = 700;
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

const GLIDER0 = [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]];
function rotCells(turns) {
  return GLIDER0.map(([dx, dy]) => {
    let a = dx, b = dy;
    for (let i = 0; i < turns; i++) { const na = -b, nb = a; a = na; b = nb; }
    return [a, b];
  });
}
function gliderAt(x, y, turns) {
  const s = new Set();
  for (const [dx, dy] of rotCells(turns)) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }

const H = 1, OX = 6, OY = 1;          // the h1(6,1) cascade
const BX = 300, BY = 300;
const STEPS = 2200;
const D = 60;                          // injection distance from the collision site

// A moving glider of heading `turns` placed at distance D from (BX,BY) along
// direction (1,1) will reach the site in D/1 * 1/ (speed 1/4) = 4D steps.
// We want arrival at (injectT + 4D), so start it at injectT.
function armA(injectT) { return [gliderAt(BX + D, BY + D, 2)]; }              // one glider, heading (-1,-1)
function armB(injectT) {                                                       // two gliders, parallel, heading (-1,-1)
  return [gliderAt(BX + D, BY + D, 2), gliderAt(BX + D + 30, BY + D, 2)];
}
function armC(injectT) {                                                       // two gliders aimed at each other
  // the h1(6,1) pair itself, shifted so that its mutual collision lands on the site
  return [gliderAt(BX + D, BY + D, 0), gliderAt(BX + D + OX, BY + D + OY, H)];
}

function runArm(armFn, injectT) {
  let cur = union(...armFn(0).map(() => new Set()));   // empty
  cur = union(gliderAt(BX, BY, 0), gliderAt(BX + OX, BY + OY, H));   // seed the cascade
  for (let t = 0; t < injectT; t++) cur = step(cur);
  cur = union(cur, ...armFn(injectT));
  const initMass = cur.size;
  let peak = initMass, peakT = injectT;
  for (let t = injectT; t < STEPS; t++) {
    cur = step(cur);
    if (cur.size > peak) { peak = cur.size; peakT = t; }
    if (cur.size === 0) break;
  }
  return { initMass, peak, peakT, final: cur.size };
}

console.log('=== pair coupling: one particle vs two vs two-correlated ===');
console.log(`torus ${SIZE}x${SIZE}; cascade h1(6,1) at (${BX},${BY}); ${STEPS} steps`);
console.log('injected from distance %d, so arrival time = injectT + %d\n', D, 4 * D);

// baseline
const base = runArm(() => [], 0);
console.log('--- baseline (cascade alone) ---');
console.log(`  peak ${base.peak}  final ${base.final}\n`);

const TIMES = [0, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 400, 600, 900];

console.log('=== ARM A: ONE glider injected ===');
const A = [];
for (const t of TIMES) { const r = runArm(armA, t); A.push({ t, ...r }); }
console.log('  injectT : ' + TIMES.map(t => String(t).padStart(6)).join(''));
console.log('  final   : ' + A.map(r => String(r.final).padStart(6)).join(''));
console.log('  peak    : ' + A.map(r => String(r.peak).padStart(6)).join(''));

console.log('\n=== ARM B: TWO gliders, parallel (no mutual collision) ===');
const B = [];
for (const t of TIMES) { const r = runArm(armB, t); B.push({ t, ...r }); }
console.log('  injectT : ' + TIMES.map(t => String(t).padStart(6)).join(''));
console.log('  final   : ' + B.map(r => String(r.final).padStart(6)).join(''));
console.log('  peak    : ' + B.map(r => String(r.peak).padStart(6)).join(''));

console.log('\n=== ARM C: TWO gliders aimed at each other (correlated pair) ===');
const C = [];
for (const t of TIMES) { const r = runArm(armC, t); C.push({ t, ...r }); }
console.log('  injectT : ' + TIMES.map(t => String(t).padStart(6)).join(''));
console.log('  final   : ' + C.map(r => String(r.final).padStart(6)).join(''));
console.log('  peak    : ' + C.map(r => String(r.peak).padStart(6)).join(''));

// ---- spread analysis
function spread(arr) {
  const fins = arr.map(r => r.final);
  const peaks = arr.map(r => r.peak);
  return {
    finMin: Math.min(...fins), finMax: Math.max(...fins),
    pkMin: Math.min(...peaks), pkMax: Math.max(...peaks),
  };
}
const sA = spread(A), sB = spread(B), sC = spread(C);

console.log('\n=== spread of outcomes across injection times ===');
console.log(`  ARM A (1 glider)      final ${sA.finMin}..${sA.finMax}   peak ${sA.pkMin}..${sA.pkMax}`);
console.log(`  ARM B (2 parallel)    final ${sB.finMin}..${sB.finMax}   peak ${sB.pkMin}..${sB.pkMax}`);
console.log(`  ARM C (2 correlated)  final ${sC.finMin}..${sC.finMax}   peak ${sC.pkMin}..${sC.pkMax}`);
console.log(`  baseline              final ${base.final}   peak ${base.peak}`);

console.log('\n=== sharpness of the response in ARM C (finest available resolution) ===');
const fine = C.filter(r => r.t <= 300);
for (const r of fine) {
  console.log(`  injectT=${String(r.t).padStart(4)}  final=${String(r.final).padStart(4)}  peak=${String(r.peak).padStart(5)}  (${(r.peak / base.peak).toFixed(2)}x baseline peak)`);
}

console.log('\n=== verdict ===');
const cSpread = sC.pkMax - sC.pkMin;
const bSpread = sB.pkMax - sB.pkMin;
const aSpread = sA.pkMax - sA.pkMin;
console.log(`  peak spread: A=${aSpread}  B=${bSpread}  C=${cSpread}`);
console.log('');
if (cSpread > 3 * Math.max(aSpread, bSpread, 1)) {
  console.log('  => CORRELATION IS WHAT COUPLES.');
  console.log('     Two particles flying in parallel behave like one particle;');
  console.log('     two particles on a collision course couple strongly and with a');
  console.log('     sharp dependence on arrival time.');
} else if (bSpread > 3 * Math.max(aSpread, 1)) {
  console.log('  => mere particle COUNT couples; internal correlation is not decisive.');
} else {
  console.log('  => A, B, C all behave similarly: the cascade is insensitive to what is');
  console.log('     injected, and the earlier large swings may have been mis-measured.');
}
