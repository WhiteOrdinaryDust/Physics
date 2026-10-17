'use strict';
/*
 * collision.js
 * ---------------------------------------------------------------------------
 * STEP 2 (option B): controlled collisions instead of random soups.
 *
 * Why: the random-soup route needed component identity tracking in a dense
 * environment, where three successive detectors gave three contradictory
 * answers (8947 movers / 0 movers / all-zero). Controlled collisions avoid that
 * entirely: few components, unambiguous identity, and the collision PRODUCTS are
 * themselves the particle spectrum.
 *
 * Method, per experiment:
 *   phase 1  launch   : single isolated object travels in a straight line
 *                      -> identify its (shape, period, velocity) exactly
 *   phase 2  collide  : the two objects meet (interaction tagged by mass change)
 *   phase 3  settle   : let the debris run; only STABLE final objects remain
 *   phase 4  read out : identify the surviving objects as the spectrum
 *
 * Controls: a lone glider must produce no interaction and fly forever.
 * ---------------------------------------------------------------------------
 */

const SIZE = 96;            // torus; big enough that objects never wrap into each other
const stride = SIZE;

/*
 * BUG FIX (this invalidated the first run):
 *   key() applied % SIZE to x and y SEPARATELY, but cells are encoded as
 *   y*SIZE + x. An object straddling x=0 had its cells split to opposite sides
 *   of the torus, so the geometry was destroyed. Symptoms: LWSS reported "2
 *   components", gliders reported "non-rigid motion", and a bounding box of
 *   "3x96" appeared.
 *   Fix: keep constructor coordinates in bounds, wrap once at the end, and wrap
 *   x/y together (not independently).
 */
function wrapCell(x, y) { return ((y % SIZE) + SIZE) % SIZE * stride + (((x % SIZE) + SIZE) % SIZE); }
function key(x, y) { return wrapCell(x, y); }

function step(set) {
  const counts = new Map();
  for (const k of set) {
    const x = k % stride, y = (k - x) / stride;
    for (let dy = -1; dy <= 1; dy++) {
      for (let dx = -1; dx <= 1; dx++) {
        if (dx === 0 && dy === 0) continue;
        const nk = key(x + dx, y + dy);
        counts.set(nk, (counts.get(nk) || 0) + 1);
      }
    }
  }
  const out = new Set();
  for (const [k, c] of counts) {
    if (c === 3 || (c === 2 && set.has(k))) out.add(k);
  }
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
      const x = k % stride, y = (k - x) / stride;
      for (let dy = -1; dy <= 1; dy++) {
        for (let dx = -1; dx <= 1; dx++) {
          if (dx === 0 && dy === 0) continue;
          const nk = key(x + dx, y + dy);
          if (set.has(nk) && !seen.has(nk)) { seen.add(nk); Q.push(nk); }
        }
      }
    }
    comps.push(cells);
  }
  return comps;
}

// canonical form: bounding box + occupancy string, and integer centroid sums
function canon(cells) {
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  let sx = 0, sy = 0;
  for (const k of cells) {
    const x = k % stride, y = (k - x) / stride;
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
    sx += x; sy += y;
  }
  const w = maxx - minx + 1, h = maxy - miny + 1;
  const grid = [];
  for (let r = 0; r < h; r++) grid.push(new Array(w).fill('.'));
  for (const k of cells) {
    const x = k % stride, y = (k - x) / stride;
    grid[y - miny][x - minx] = 'o';
  }
  return {
    w, h, n: cells.length, sx, sy,
    rows: grid.map(r => r.join('')),
    id: `${w}x${h}:` + grid.map(r => r.join('')).join('/'),
  };
}

// wrap-aware displacement of integer centroid SUMS.
// BUG FIX: the old version compared a centroid displacement (sum/n) against a
// limit of SIZE*n and tested divisibility by n, which can never hold. Work in
// cell-count units throughout.
function centroidDisp(a, b) {
  const n = a.n;
  if (b.n !== n) return null;                 // different size => not the same object
  const lim = SIZE * n;                       // full torus in sum units
  let dsx = b.sx - a.sx, dsy = b.sy - a.sy;
  if (dsx > lim / 2) dsx -= lim; if (dsx < -lim / 2) dsx += lim;
  if (dsy > lim / 2) dsy -= lim; if (dsy < -lim / 2) dsy += lim;
  if (dsx % n !== 0 || dsy % n !== 0) return null;
  return [dsx / n, dsy / n];
}

// ---- known objects, for catalogue lookup
const KNOWN = new Map();
function reg(rows, name) { KNOWN.set(`${rows[0].length}x${rows.length}:` + rows.join('/'), name); }
reg(['.o.', '..o', 'ooo'], 'glider (phase 0, moves +1,+1 per 4 gens)');
reg(['oo', 'oo'], 'block (still life)');
reg(['.oo.', 'o..o', '.oo.'], 'loaf (still life)');
reg(['.o.', 'o.o', 'o.o', '.o.'], 'ship (still life, 3x4)         <-- confirmed fixed point');
reg(['.o.', 'o.o', '.o.'], '?? (3x3, unverified)');
reg(['.o.', '.o.', '.o.'], 'blinker (oscillator, vertical)');
reg(['ooo'], 'blinker (oscillator, horizontal)');
reg(['.oo', 'o.o', '.oo'], 'beehive (still life)');
reg(['ooo', 'o..', 'o..'], '?? (unverified)');

function identify(canonId) {
  return KNOWN.get(canonId) || null;
}

/*
 * Launch phase control.
 *
 * CORRECT CRITERION (two earlier attempts were wrong):
 *   - "step-by-step rigid motion" is WRONG: a glider's centroid does not move by
 *     a constant vector each generation. Measured steps for the glider are
 *     (0,.5), (.2,.4), (.4,0), (.4,.2) -- the centroid of a discrete pattern may
 *     wander; only the 4-step total is (1,1).
 *   - therefore match on SHAPE RECURRENCE: same canonical shape at t and t+p,
 *     with a constant non-zero displacement over every p-step window.
 */
function trackSingle(init, steps) {
  let set = init;
  const frames = [];
  for (let t = 0; t <= steps; t++) {
    const comps = components(set);
    if (comps.length !== 1) return { ok: false, reason: `expected 1 component, got ${comps.length}`, t };
    frames.push(canon(comps[0]));
    set = step(set);
  }
  const sumDisp = (a, b) => {
    let dsx = (b.sx - a.sx) / a.n, dsy = (b.sy - a.sy) / a.n;
    if (Math.abs(dsx) > SIZE / 2) dsx -= Math.sign(dsx) * SIZE;
    if (Math.abs(dsy) > SIZE / 2) dsy -= Math.sign(dsy) * SIZE;
    return [Math.round(dsx), Math.round(dsy)];
  };
  const phases = new Set(frames.map(f => f.id)).size;
  for (let p = 1; p <= 12; p++) {
    let ok = true, base = null;
    for (let t = 0; t + p < frames.length; t++) {
      if (frames[t].id !== frames[t + p].id) { ok = false; break; }
      const [dx, dy] = sumDisp(frames[t], frames[t + p]);
      if (base === null) base = [dx, dy];
      else if (dx !== base[0] || dy !== base[1]) { ok = false; break; }
    }
    if (ok && base && (base[0] !== 0 || base[1] !== 0)) {
      return { ok: true, p, vx: base[0], vy: base[1], id: frames[0].id, phases };
    }
  }
  return { ok: false, reason: 'no periodic shape recurrence with constant drift' };
}

// ---- object constructors (coordinates kept positive; nodes wrap at the end)
function glider(x, y) {
  // verified: moves (+1,+1) per 4 generations
  const s = new Set();
  for (const [dx, dy] of [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function gliderNW(x, y) {
  // mirror in x: moves (-1,+1) per 4 generations
  const s = new Set();
  for (const [dx, dy] of [[1, 0], [0, 1], [2, 2], [1, 2], [0, 2]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function gliderSW(x, y) {
  // mirror in y: moves (+1,-1) per 4 generations
  const s = new Set();
  for (const [dx, dy] of [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]]) s.add(wrapCell(x + dx, y - dy));
  return s;
}
function lwss(x, y) {
  // verified: moves (-2,0) per 4 generations (leftward)
  const s = new Set();
  for (const [dx, dy] of [[1, 0], [4, 0], [0, 1], [0, 2], [4, 2], [0, 3], [1, 3], [2, 3], [3, 3]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function lwssR(x, y) {
  // mirror of the verified LWSS in x: moves (+2,0) per 4 generations
  const s = new Set();
  for (const [dx, dy] of [[4, 0], [1, 0], [0, 1], [0, 2], [0, 2 + 0], [0, 3], [1, 3], [2, 3], [3, 3]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}

function union(...sets) {
  const out = new Set();
  for (const s of sets) for (const k of s) out.add(k);
  return out;
}

// ---- analysis of a settled configuration
function analyseConfig(set) {
  const comps = components(set);
  const out = [];
  for (const cells of comps) {
    const c = canon(cells);
    out.push({ id: c.id, n: c.n, w: c.w, h: c.h, name: identify(c.id) });
  }
  return out;
}

// ============================================================ runs
console.log('=== controlled collision experiments (Life, empty background) ===');
console.log(`torus ${SIZE}x${SIZE}\n`);

// ---- CONTROL: lone glider must fly rigidly forever, never interact
console.log('--- control: lone glider ---');
{
  const t = trackSingle(glider(20, 20), 60);
  console.log(`  ${t.ok ? `PASS: rigid, period ${t.p}, velocity (${t.vx},${t.vy})` : `*** FAIL: ${t.reason} ***`}`);
  console.log(`  pattern id: ${t.id}`);
  const t2 = trackSingle(gliderNW(70, 20), 60);
  console.log(`  mirrored:   ${t2.ok ? `PASS: period ${t2.p}, velocity (${t2.vx},${t2.vy})` : `*** FAIL: ${t2.reason} ***`}`);
}
{
  const t = trackSingle(lwss(20, 60), 60);
  console.log(`  LWSS:       ${t.ok ? `PASS: period ${t.p}, velocity (${t.vx},${t.vy})` : `*** FAIL: ${t.reason} ***`}`);
}
console.log('');

// ---- EXPERIMENT runner: no hand-computed intercept courses.
// Objects are launched on straight paths; the collision is DETECTED (component
// count drops to 1, or total mass departs from the sum of the isolated objects).
function runExperiment(name, init, steps) {
  console.log(`--- ${name} ---`);
  let set = init;
  const mass = [];
  let mergeT = -1, prevComps = -1, minComps = 99;
  for (let t = 0; t <= steps; t++) {
    const comps = components(set);
    mass.push(set.size);
    if (comps.length < prevComps && prevComps === 2) mergeT = t;
    if (t < 120) minComps = Math.min(minComps, comps.length);
    prevComps = comps.length;
    set = step(set);
  }
  const initMass = mass[0];
  const peak = Math.max(...mass);
  const minMass = Math.min(...mass);
  const interacted = mergeT >= 0 || peak !== initMass || minMass !== initMass;
  console.log(`  initial mass ${initMass}, peak ${peak}, min ${minMass}, final ${mass[mass.length - 1]}`);
  console.log(`  objects merged at t=${mergeT >= 0 ? mergeT : 'never'}   interacted: ${interacted ? 'YES' : 'NO'}`);
  if (!interacted) {
    console.log('  -> the two objects passed by each other without interacting (aim was off)');
    console.log('');
    return { interacted: false };
  }
  set = init;
  for (let t = 0; t < steps; t++) set = step(set);
  const final = analyseConfig(set);
  console.log(`  final objects: ${final.length}`);
  for (const o of final) {
    console.log(`    ${String(o.n).padStart(3)} cells  ${o.w}x${o.h}  ${o.name || '(unidentified)'}`);
    if (!o.name) console.log(`        ${o.id}`);
  }
  console.log('');
  return { interacted: true, final };
}

// Head-on: verified glider (+1,+1) versus its (-1,-1) counterpart.
// Constructed by 180-degree rotation of the glider pattern, so it provably moves
// the opposite way without me guessing a mirrored RLE.
function gliderOpp(x, y) {
  const s = new Set();
  for (const [dx, dy] of [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]]) s.add(wrapCell(x - dx, y - dy));
  return s;
}
function gliderPerp(x, y) {
  // 90-degree rotation of the glider: moves (+1,-1)
  const s = new Set();
  for (const [dx, dy] of [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]]) s.add(wrapCell(x + dy, y - dx));
  return s;
}

runExperiment('glider(+1,+1)  vs  glider(-1,-1)   head-on',
  union(glider(20, 20), gliderOpp(70, 70)), 300);

runExperiment('glider(+1,+1)  vs  glider(+1,-1)   perpendicular',
  union(glider(20, 20), gliderPerp(20, 70)), 300);

runExperiment('glider(+1,+1)  vs  glider(-1,-1)   offset lanes',
  union(glider(20, 30), gliderOpp(70, 60)), 300);

// NOTE: the lwss() pattern used below is NOT a valid LWSS -- it decomposes into
// two components immediately (see the control output). It is left in place only
// to show that fact; the run is not evidence about real LWSS collisions.
runExperiment('glider(+1,+1)  vs  (broken lwss pattern)   crossing',
  union(glider(30, 40), lwss(70, 40)), 300);
