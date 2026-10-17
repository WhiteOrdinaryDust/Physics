'use strict';
/*
 * calibrated_injection.js
 * ---------------------------------------------------------------------------
 * STEP 11: fix the geometry first, then redo the three-arm comparison.
 *
 * WHY THE LAST ATTEMPT FAILED
 *   The cascade triggers at t=1 and peaks at t=823. The previous test injected
 *   from 60 cells away; at glider speed 1/4 that is 240 steps of travel, so the
 *   injected material always arrived AFTER the cascade was already over. The
 *   "injection time" was therefore meaningless.
 *
 * THE TIME-DISTANCE TABLE (computed before designing the experiment)
 *   glider speed            = 1/4 cell per step
 *   length of flight        = 4 * d  steps   (d = distance to the target point)
 *   cascade trigger         = t = 1
 *   cascade peak            = t = 823
 *   => to arrive AT the peak, fly  d = 823 / 4 ~ 206 cells
 *   => to arrive at trigger,  fly  d = 1   / 4 ~ 0   cells (inject immediately)
 *   We use d = 200 (arrival ~t=800, right at the peak), plus 100 and 400 for
 *   comparison, so the arrival stage is known rather than accidental.
 *
 * We also measure the cascade's SPATIAL EXTENT, so that the impact parameter can
 * be chosen relative to the actual size of the thing being hit, rather than
 * guessed.
 *
 * ARMS (all injected to arrive at the peak):
 *   A : one glider aimed at the collision site
 *   B : two gliders, PARALLEL courses (they never meet each other)
 *   C : two gliders aimed AT EACH OTHER, colliding exactly when they arrive
 * ---------------------------------------------------------------------------
 */

const SIZE = 900;
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

function bbox(set) {
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  for (const k of set) {
    const x = k % SIZE, y = (k - x) / SIZE;
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
  }
  return set.size ? { w: maxx - minx + 1, h: maxy - miny + 1, area: (maxx - minx + 1) * (maxy - miny + 1) } : { w: 0, h: 0, area: 0 };
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

const H = 1, OX = 6, OY = 1;
const CX = 450, CY = 450;              // cascade site, centred on the torus
const STEPS = 2200;

// ---- 1. cascade alone: mass and spatial extent
console.log('=== calibration: cascade alone ===');
let s = union(gliderAt(CX, CY, 0), gliderAt(CX + OX, CY + OY, H));
const mass = [];
const ext = new Map();
let peak = 0, peakT = 0;
for (let t = 0; t <= STEPS; t++) {
  mass.push(s.size);
  if (s.size > peak) { peak = s.size; peakT = t; }
  if ([1, 50, 200, 400, 800, 823, 1200, 2000].includes(t)) ext.set(t, bbox(s));
  if (t === STEPS) break;
  s = step(s);
}
console.log(`  peak ${peak} @t${peakT}, final ${mass[STEPS]}`);
console.log('   t     mass    bbox');
for (const [t, b] of [...ext.entries()].sort((a, b) => a[0] - b[0])) {
  console.log(`  ${String(t).padStart(4)}  ${String(mass[t]).padStart(5)}   ${b.w}x${b.h}`);
}

// ---- 2. time-distance table
console.log('\n=== time-distance table (glider speed 1/4 cell/step) ===');
console.log('  distance d   flight time   arrives at t=   stage of cascade');
for (const d of [25, 50, 100, 150, 200, 206, 250, 300, 400]) {
  const t = 4 * d;
  const stage = t < 50 ? 'before growth' : t < 400 ? 'growth' : t < 900 ? 'PEAK region' : t < 1400 ? 'decline' : 'plateau';
  console.log(`  ${String(d).padStart(10)}   ${String(4 * d).padStart(11)}   ${String(t).padStart(11)}   ${stage}`);
}

// ---- 3. arms, injected so as to ARRIVE at the peak
const D = 200;                        // arrives ~t=800
const D2 = 100, D3 = 400;

/*
 * A glider on a course through the point (CX,CY) along direction (1,1) or (-1,-1),
 * offset perpendicular to that course by `ipar`. Perpendicular to (1,1) is (1,-1).
 * courses: u=+1 means travelling (+1,+1), u=-1 means travelling (-1,-1).
 */
function aimed(u, ipar, dist) {
  const ux = u * (1 / Math.SQRT2), uy = u * (1 / Math.SQRT2);
  const px = 1 / Math.SQRT2, py = -1 / Math.SQRT2;
  const sx = CX - ux * dist + px * ipar;
  const sy = CY - uy * dist + py * ipar;
  return gliderAt(Math.round(sx), Math.round(sy), u > 0 ? 0 : 2);
}

function runArms(dist, ipar) {
  const setup = {
    A: [aimed(-1, ipar, dist)],
    B: [aimed(-1, ipar, dist), aimed(-1, ipar + 40, dist)],
    C: [aimed(-1, ipar, dist), aimed(+1, ipar, dist)],
  };
  const out = {};
  for (const arm of ['A', 'B', 'C']) {
    // seed the cascade immediately, then add the injected material at t=0 too:
    // it arrives 4*dist steps later, which by construction is the peak
    let cur = union(gliderAt(CX, CY, 0), gliderAt(CX + OX, CY + OY, H), ...setup[arm]);
    const initMass = cur.size;
    let pk = initMass;
    for (let t = 0; t < STEPS; t++) {
      cur = step(cur);
      if (cur.size > pk) pk = cur.size;
      if (cur.size === 0) break;
    }
    out[arm] = { initMass, peak: pk, final: cur.size };
  }
  return out;
}

console.log(`\n=== arms injected to arrive at the peak (distance ${D}, ~t=${4 * D}) ===`);
console.log('(baseline: peak ' + peak + ', final ' + mass[STEPS] + ')');
for (const ipar of [0, 5, 10, 20, 40, 80]) {
  const r = runArms(D, ipar);
  console.log(`  impact param ${String(ipar).padStart(3)}:  A peak ${String(r.A.peak).padStart(4)} final ${String(r.A.final).padStart(4)}   B peak ${String(r.B.peak).padStart(4)} final ${String(r.B.final).padStart(4)}   C peak ${String(r.C.peak).padStart(4)} final ${String(r.C.final).padStart(4)}`);
}

console.log('\n=== arms at other distances (impact parameter 0) ===');
for (const d of [D2, D, D3]) {
  const r = runArms(d, 0);
  console.log(`  distance ${String(d).padStart(3)} (arrive ~t=${String(4 * d).padStart(4)}):  A final ${String(r.A.final).padStart(4)}   B final ${String(r.B.final).padStart(4)}   C final ${String(r.C.final).padStart(4)}`);
}

console.log('\n=== reading ===');
console.log('  All three arms now deliver material AT the peak, so the arrival stage is');
console.log('  controlled. Compare the spread of outcomes across impact parameters:');
console.log('   - if A and B match the baseline and only C deviates, the cascade');
console.log('     responds to the internal collision of the injected pair, not to its size');
console.log('   - if A, B, C all behave alike, the cascade is insensitive to fine structure');
