'use strict';
/*
 * z4_dynamics.js
 * ---------------------------------------------------------------------------
 * THE DECIDING TEST for route C: can a Z4-symmetric pattern sweep a RICH set of
 * internal states, or only a handful?
 *
 * Background. Earlier we found (via the self-checks) that on the integer lattice a
 * 90-degree rotation about any lattice-consistent centre requires every cell to
 * lie on one of the two axes through that centre -- a 2x2 block is NOT Z4, a
 * cross is. So the Z4-symmetric shapes available are crosses.
 *
 * What route C needs. Complex (or at least continuous) amplitude structure is to
 * emerge from DISCRETE phases. For that to be credible, a Z4-symmetric particle
 * must sweep a rich internal orbit. If every Z4 pattern is a still life or a
 * period-2/3 blinker, then there are only 1-3 distinguishable internal states and
 * there is no continuum to emerge from. Route C dies.
 *
 * MEASURED
 *   for each Z4-symmetric seed: does it survive, what is its period, and how many
 *   DISTINCT internal states does it visit per period?
 * ---------------------------------------------------------------------------
 */

const SIZE = 80;
function wrapCell(x, y) { return ((y % SIZE) + SIZE) % SIZE * SIZE + (((x % SIZE) + SIZE) % SIZE); }

function lifeStep(set) {
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

// canonical signature of a configuration (translation-invariant), for period finding
function signature(set) {
  if (set.size === 0) return 'empty';
  const cells = [...set];
  let minx = Infinity, miny = Infinity;
  for (const k of cells) { const x = k % SIZE, y = (k - x) / SIZE; if (x < minx) minx = x; if (y < miny) miny = y; }
  // use the centroid's fractional part plus the shape, so translation does not
  // change the signature but internal phase does
  let sx = 0, sy = 0;
  const rel = [];
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    sx += x; sy += y;
    rel.push([x - minx, y - miny]);
  }
  rel.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  const shape = rel.map(p => p.join(',')).join(';');
  const n = cells.length;
  // fractional part of the centroid identifies the phase within a period when the
  // object translates; include it so a translating object is not "period 1"
  const fx = ((sx / n) % 1 + 1) % 1, fy = ((sy / n) % 1 + 1) % 1;
  return `${shape}|${fx.toFixed(4)},${fy.toFixed(4)}`;
}

// pure shape signature, ignoring phase (for still-life / oscillator detection)
function shapeOnly(set) {
  if (set.size === 0) return 'empty';
  const cells = [...set];
  let minx = Infinity, miny = Infinity;
  for (const k of cells) { const x = k % SIZE, y = (k - x) / SIZE; if (x < minx) minx = x; if (y < miny) miny = y; }
  const rel = cells.map(k => { const x = k % SIZE, y = (k - x) / SIZE; return [x - minx, y - miny]; });
  rel.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  return rel.map(p => p.join(',')).join(';');
}

// ---- Z4 cross seeds: arms of length a along +x,-x,+y,-y
function cross(arm) {
  const s = new Set();
  for (let d = -arm; d <= arm; d++) {
    s.add(wrapCell(40 + d, 40));
    s.add(wrapCell(40, 40 + d));
  }
  return s;
}
// thicker crosses: arms with width w
function fatCross(arm, width) {
  const s = new Set();
  for (let d = -arm; d <= arm; d++) {
    for (let o = -Math.floor(width / 2); o <= Math.floor(width / 2); o++) {
      s.add(wrapCell(40 + d, 40 + o));
      s.add(wrapCell(40 + o, 40 + d));
    }
  }
  return s;
}

function analyse(name, init, steps) {
  const initMass = init.size;
  const seenSig = new Map();
  const seenShape = new Map();
  let cur = init, period = -1, shapePeriod = -1, died = false;
  const states = [];
  for (let t = 0; t <= steps; t++) {
    const sig = signature(cur);
    const sh = shapeOnly(cur);
    states.push(sig);
    if (seenSig.has(sig)) { period = t - seenSig.get(sig); break; }
    seenSig.set(sig, t);
    if (seenShape.has(sh) && shapePeriod < 0) shapePeriod = t - seenShape.get(sh);
    seenShape.set(sh, t);
    if (cur.size === 0) { died = true; break; }
    cur = lifeStep(cur);
  }
  const distinct = new Set(states).size;
  const verdict = died ? 'DIES'
    : period === 1 ? 'still life'
      : period > 0 ? `period ${period}`
        : 'not settled';
  console.log(`  ${name.padEnd(26)} mass ${String(initMass).padStart(3)}  ${verdict.padEnd(14)} distinct internal states: ${distinct}`);
  return { name, initMass, period, distinct, died };
}

console.log('=== Z4-symmetric (cross) patterns in Life: how rich are their orbits? ===\n');

const results = [];
console.log('--- pure crosses, arm length 1..8 ---');
for (let a = 1; a <= 8; a++) {
  results.push(analyse(`cross arm=${a}`, cross(a), 400));
}
console.log('');

console.log('--- thicker crosses ---');
for (let a = 1; a <= 5; a++) {
  for (const w of [2, 3]) {
    results.push(analyse(`cross arm=${a} width=${w}`, fatCross(a, w), 400));
  }
}
console.log('');

console.log('--- crosses with a gap in the arms (arm = on/off pattern) ---');
for (let a = 3; a <= 6; a++) {
  // arm cells at distance 1 and a, skipping the middle
  const s = new Set();
  for (const d of [1, 2, a - 1, a]) { s.add(wrapCell(40 + d, 40)); s.add(wrapCell(40 - d, 40)); s.add(wrapCell(40, 40 + d)); s.add(wrapCell(40, 40 - d)); }
  results.push(analyse(`hollow cross span=${a}`, s, 400));
}
console.log('');

// ---- summary
console.log('=== summary: is there ANY Z4 pattern with a rich internal orbit? ===');
const alive = results.filter(r => !r.died);
console.log(`  patterns tested      : ${results.length}`);
console.log(`  that survive         : ${alive.length}`);
const byPeriod = new Map();
for (const r of alive) byPeriod.set(r.period, (byPeriod.get(r.period) || 0) + 1);
for (const [p, c] of [...byPeriod.entries()].sort((a, b) => a[0] - b[0])) {
  console.log(`    period ${String(p).padStart(3)}: ${c}`);
}
const best = alive.slice().sort((a, b) => (b.period || 0) - (a.period || 0))[0];
if (best) {
  console.log(`  longest period found : ${best.period} ("${best.name}", ${best.distinct} distinct states)`);
}
console.log('');
console.log('=== reading ===');
console.log('  Route C needs a Z4 pattern whose orbit sweeps MANY distinguishable');
console.log('  internal states, so that a continuous phase could emerge as a');
console.log('  coarse-grained limit. A maximum period of 1-4 gives at most a handful of');
console.log('  states, which is not a continuum and not enough for a cos(theta) law.');
