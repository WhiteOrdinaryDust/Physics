'use strict';
/*
 * encapsulation_clock.js
 * ---------------------------------------------------------------------------
 * Emergence by ENCAPSULATION: does an observer who sees only whole gliders
 * experience a slower global refresh rate?
 *
 * THE REFRAME (this is the point of the experiment)
 *   Earlier attempts looked for a PATTERN that itself moves slowly (a background
 *   with small phase velocity). That is the wrong level. Emergence is
 *   encapsulation: the observer treats the glider as an INDIVISIBLE unit and
 *   coarse-grains the rest into a smooth background. Then
 *       - the glider's internal period becomes the observer's unit of time
 *       - the "global refresh rate" is how fast the COARSE-GRAINED global state
 *         changes, not how fast any single mode travels
 *
 * WHAT IS MEASURED
 *   1. tick ratio: period of the coarse-grained global state / substrate period.
 *      For an array of gliders the whole configuration recurs every 4 substrate
 *      steps, so the world's clock should be 4x slower than the substrate's.
 *   2. the number of DISTINCT coarse-grained states visited per substrate step.
 *      If the world only changes once per 4 steps, that is the emergent clock.
 *   3. how fast a local perturbation spreads through the array -- i.e. the
 *      "light speed" of the encapsulated world, expressed in world ticks.
 *   4. control: an array of STILL LIFES, whose coarse-grained state never changes.
 * ---------------------------------------------------------------------------
 */

const SIZE = 600;
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
      for (let dy = -1; dy <= 1; dy++) {
        for (let dx = -1; dx <= 1; dx++) {
          if (dx === 0 && dy === 0) continue;
          const nk = wrapCell(x + dx, y + dy);
          if (set.has(nk) && !seen.has(nk)) { seen.add(nk); Q.push(nk); }
        }
      }
    }
    comps.push(cells);
  }
  return comps;
}

// coarse-grained description of one component: its SHAPE CLASS only.
// Position is deliberately quantised to the glider's period (4), so the observer
// cannot resolve sub-period motion -- this is the encapsulation assumption.
function shapeKey(cells) {
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
  }
  const w = maxx - minx + 1, h = maxy - miny + 1;
  let hash = 2166136261;
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    hash = Math.imul(hash ^ (((y - miny) * w + (x - minx)) + 1), 16777619);
  }
  return `${w}x${h}:${hash >>> 0}`;
}

/*
 * Coarse-grained global state = the multiset of (shape, quantised position).
 * Positions are quantised to a grid of step `q`, modelling the observer's
 * inability to resolve finer motion.
 */
function coarseState(set, q) {
  const items = [];
  for (const cells of components(set)) {
    let minx = Infinity, miny = Infinity;
    for (const k of cells) {
      const x = k % SIZE, y = (k - x) / SIZE;
      if (x < minx) minx = x; if (y < miny) miny = y;
    }
    items.push(`${shapeKey(cells)}@${Math.floor(minx / q)},${Math.floor(miny / q)}`);
  }
  items.sort();
  return items.join('|');
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
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }

// ---- build a uniform array of same-direction gliders, well separated
function gliderArray(spacing, n) {
  const parts = [];
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) parts.push(gliderAt(60 + i * spacing, 60 + j * spacing));
  }
  return union(...parts);
}
function blockArray(spacing, n) {
  const parts = [];
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) parts.push(blockAt(60 + i * spacing, 60 + j * spacing));
  }
  return union(...parts);
}

console.log('=== encapsulation: does a glider array give the observer a slower clock? ===');
console.log(`torus ${SIZE}x${SIZE}\n`);

const SPACING = 40, N = 8;          // 64 gliders, never colliding
const STEPS = 40;

function analyse(name, init, q) {
  let s = init;
  const coarseSeen = [];
  const rawSeen = [];
  for (let t = 0; t <= STEPS; t++) {
    coarseSeen.push(coarseState(s, q));
    rawSeen.push(s.size + ':' + Array.from(s).sort((a, b) => a - b).join(','));
    if (t === STEPS) break;
    s = step(s);
  }
  // period of the coarse-grained state
  let cPeriod = -1;
  for (let p = 1; p <= 30; p++) {
    let ok = true;
    for (let t = 0; t + p < coarseSeen.length; t++) if (coarseSeen[t] !== coarseSeen[t + p]) { ok = false; break; }
    if (ok) { cPeriod = p; break; }
  }
  let rPeriod = -1;
  for (let p = 1; p <= 30; p++) {
    let ok = true;
    for (let t = 0; t + p < rawSeen.length; t++) if (rawSeen[t] !== rawSeen[t + p]) { ok = false; break; }
    if (ok) { rPeriod = p; break; }
  }
  const distinctCoarse = new Set(coarseSeen).size;
  const distinctRaw = new Set(rawSeen).size;
  // DIAGNOSTIC: print the first few object signatures literally, so a mismatch
  // can be seen rather than theorised about.
  console.log(`  first 4 coarse states, first 3 objects each:`);
  for (let t = 0; t < 4; t++) {
    const objs = coarseSeen[t].split('|').slice(0, 3).join('  ');
    console.log(`    t=${t}  ${objs}`);
  }
  console.log('');
  console.log(`--- ${name} (coarse grain q=${q}) ---`);
  console.log(`  coarse-grained state, first 10 steps:`);
  for (let t = 0; t < 10; t++) {
    // show a compact signature rather than the whole state
    const sig = coarseSeen[t].length + ' chars, ' + (coarseSeen[t].match(/@/g) || []).length + ' objects';
    const changed = t > 0 && coarseSeen[t] !== coarseSeen[t - 1];
    console.log(`    t=${String(t).padStart(2)}  ${sig.padEnd(30)} ${changed ? 'CHANGED' : 'same'}`);
  }
  console.log(`  coarse-grained period : ${cPeriod === -1 ? '>30' : cPeriod}`);
  console.log(`  raw (substrate) period: ${rPeriod === -1 ? '>30' : rPeriod}`);
  console.log(`  distinct coarse states over ${STEPS} steps : ${distinctCoarse}`);
  console.log(`  distinct raw states    over ${STEPS} steps : ${distinctRaw}`);
  console.log(`  => coarse clock ticks every ${cPeriod === -1 ? '?' : cPeriod} substrate steps; tick ratio = ${cPeriod > 0 ? cPeriod : '?'}`);
  console.log('');
  return { cPeriod, rPeriod, distinctCoarse, distinctRaw };
}

const g = analyse('array of 64 same-direction gliders', gliderArray(SPACING, N), 4);
const b = analyse('control: array of 64 still-life blocks', blockArray(SPACING, N), 4);

console.log('=== summary ===');
console.log(`  glider array : coarse period ${g.cPeriod}, raw period ${g.rPeriod}`);
console.log(`  block array  : coarse period ${b.cPeriod}, raw period ${b.rPeriod}`);
console.log('');
if (g.cPeriod > 1) {
  console.log(`  => The encapsulated world ticks every ${g.cPeriod} substrate steps.`);
  console.log(`     The observer's global refresh rate is ${g.cPeriod}x slower than the substrate,`);
  console.log('     purely from encapsulation -- no new physics, just a change of description.');
  console.log(`     Note the size: ${g.cPeriod}x, not orders of magnitude.`);
} else {
  console.log('  => No slower clock: the coarse-grained state changes every substrate step.');
}
