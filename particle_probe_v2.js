'use strict';
/*
 * particle_probe_v2.js
 * ---------------------------------------------------------------------------
 * STEP 1: the simplest usable definition of a "particle", WITH POSITIVE CONTROLS.
 *
 * DEFINITION
 *   A particle is a finite connected localized pattern whose shape is periodic
 *   in time (period p) and which translates rigidly at constant velocity over
 *   that period:  shape(t+p) == shape(t),  centroid(t+p) - centroid(t) = p * v.
 *
 *   p = 1 and v = 0  -> still life (not a particle)
 *   p > 1 and v = 0  -> oscillator (not a particle)
 *   v != 0           -> particle, with period p
 *
 * PREVIOUS VERSION WAS INVALID: it had no positive control, tracked connected
 * runs (which is the wrong granularity -- a glider is not a single run), and
 * mixed W=256 with a 4096 stride in two different places.
 *
 * CONTROLS
 *   Life glider        : known to move (1,1) every 4 generations -> MUST be found
 *   Life blinker       : known oscillator, v=0                    -> must NOT be a particle
 *   rule 110 glider    : seeded explicitly                        -> MUST be found
 *   rule 30 random     : chaotic                                  -> negative control
 * ---------------------------------------------------------------------------
 */

const HASH = 0;

// ============================================================ Life
function makeLife(W) {
  const stride = W;
  return {
    stride,
    step(set) {
      const counts = new Map();
      for (const k of set) {
        const x = k % stride, y = (k - x) / stride;
        for (let dy = -1; dy <= 1; dy++) {
          for (let dx = -1; dx <= 1; dx++) {
            if (dx === 0 && dy === 0) continue;
            const nx = (x + dx + W) % W, ny = (y + dy + W) % W;
            const nk = ny * stride + nx;
            counts.set(nk, (counts.get(nk) || 0) + 1);
          }
        }
      }
      const out = new Set();
      for (const [k, c] of counts) {
        if (c === 3 || (c === 2 && set.has(k))) out.add(k);
      }
      return out;
    },
  };
}

function components(set, W) {
  const stride = W;
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
          const nx = (x + dx + W) % W, ny = (y + dy + W) % W;
          const nk = ny * stride + nx;
          if (set.has(nk) && !seen.has(nk)) { seen.add(nk); Q.push(nk); }
        }
      }
    }
    comps.push(cells);
  }
  return comps;
}

// canonical shape: bounding box, occupancy bitmask, and true centroid offset
function canon(cells, W) {
  const stride = W;
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  let sx = 0, sy = 0;
  for (const k of cells) {
    const x = k % stride, y = (k - x) / stride;
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
    sx += x; sy += y;
  }
  const w = maxx - minx + 1, h = maxy - miny + 1;
  let hash = 2166136261;
  for (const k of cells) {
    const x = k % stride, y = (k - x) / stride;
    const idx = (y - miny) * w + (x - minx);
    hash = Math.imul(hash ^ (idx + 1), 16777619);
  }
  return {
    w, h, hash: hash >>> 0, n: cells.length,
    cx: sx / cells.length, cy: sy / cells.length, minx, miny,
  };
}

// wrap-aware displacement between two centroids
function disp(ax, ay, bx, by, W) {
  let dx = bx - ax, dy = by - ay;
  if (dx > W / 2) dx -= W; if (dx < -W / 2) dx += W;
  if (dy > W / 2) dy -= W; if (dy < -W / 2) dy += W;
  return [dx, dy];
}

/*
 * Track components and classify by (shapeHash, period) -> rigid body motion.
 * history: array of Maps, newest last, each Map shapeKey -> {cx,cy,w,h,n}
 */
function classifyLife(history, W, minSize, maxPeriod, minRepeats) {
  const found = new Map();     // "wxh:period p:v(dx,dy)" -> repeats
  const now = history[history.length - 1];
  for (const [key, rec] of now) {
    if (rec.n < minSize) continue;
    for (let p = 1; p <= maxPeriod; p++) {
      const old = history[history.length - 1 - p];
      if (!old) continue;
      const prev = old.get(key);
      if (!prev) continue;
      const [dx, dy] = disp(prev.cx, prev.cy, rec.cx, rec.cy, W);
      const isMover = (dx !== 0 || dy !== 0);
      if (!isMover) continue;                       // exclude oscillators / still lifes
      const tag = p === 1 ? `v(${dx},${dy})` : `p${p}:v(${dx},${dy})`;
      const label = `${rec.w}x${rec.h}:${tag}`;
      found.set(label, (found.get(label) || 0) + 1);
    }
  }
  return found;
}

function runLifeAnalysis(initSet, W, T, minSize, maxPeriod) {
  const L = makeLife(W);
  let set = initSet;
  const history = [];
  const agg = new Map();
  for (let t = 0; t < T; t++) {
    set = L.step(set);
    if (set.size === 0) break;
    const frame = new Map();
    for (const cells of components(set, W)) {
      if (cells.length < minSize) continue;
      const c = canon(cells, W);
      const key = `${c.w}:${c.h}:${c.hash}`;
      frame.set(key, { cx: c.cx, cy: c.cy, w: c.w, h: c.h, n: c.n });
    }
    history.push(frame);
    if (history.length > maxPeriod + 1) history.shift();
    for (const [k, v] of classifyLife(history, W, minSize, maxPeriod, 1)) {
      agg.set(k, (agg.get(k) || 0) + v);
    }
  }
  return agg;
}

// ============================================================ rule 110 (1D), tracked by tile matching
function rule110Step(row, W) {
  const out = new Uint8Array(W);
  for (let i = 0; i < W; i++) {
    const l = row[(i - 1 + W) % W], c = row[i], r = row[(i + 1) % W];
    const idx = (l << 2) | (c << 1) | r;
    out[i] = (110 >> idx) & 1;
  }
  return out;
}

/*
 * Generic 1D particle detector: take the local neighborhood of a tile of width w
 * starting at position x, canonicalize it, and look for period-p shape recurrence
 * with constant drift.
 */
function scanECA(rule, init, W, T, tileW, maxPeriod) {
  const tbl = new Uint8Array(8);
  for (let n = 0; n < 8; n++) tbl[n] = (rule >> n) & 1;
  const step = (row) => {
    const out = new Uint8Array(W);
    for (let i = 0; i < W; i++) {
      out[i] = tbl[(row[(i - 1 + W) % W] << 2) | (row[i] << 1) | row[(i + 1) % W]];
    }
    return out;
  };
  // canonicalize a tile by rotation: rotate so it starts with the minimal pattern
  const canonTile = (row, x) => {
    let best = null;
    for (let s = 0; s < tileW; s++) {
      let str = '';
      for (let i = 0; i < tileW; i++) str += row[(x + s + i) % W];
      if (best === null || str < best.str) best = { str, s };
    }
    return best.str;
  };
  const history = [];
  const hits = new Map();
  let row = Uint8Array.from(init);
  for (let t = 0; t < T; t++) {
    const frame = new Map();
    for (let x = 0; x < W; x++) {
      const key = canonTile(row, x);
      if (!frame.has(key)) frame.set(key, x); // first occurrence
    }
    history.push(frame);
    if (history.length > maxPeriod + 1) history.shift();
    for (const [key, x] of frame) {
      for (let p = 1; p <= maxPeriod; p++) {
        const old = history[history.length - 1 - p];
        if (!old) continue;
        const x0 = old.get(key);
        if (x0 === undefined) continue;
        let dx = x - x0;
        if (dx > W / 2) dx -= W; if (dx < -W / 2) dx += W;
        if (dx === 0) continue;
        const label = `tile${tileW}:p${p}:v${dx}`;
        hits.set(label, (hits.get(label) || 0) + 1);
      }
    }
    row = step(row);
  }
  return hits;
}

// ============================================================ positive controls
function lifeGlider(W, ox, oy) {
  const s = new Set();
  const pat = [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]];
  for (const [dx, dy] of pat) {
    const x = (ox + dx) % W, y = (oy + dy) % W;
    s.add(y * W + x);
  }
  return s;
}
function lifeBlinker(W, ox, oy) {
  const s = new Set();
  for (let dx = 0; dx < 3; dx++) {
    const x = (ox + dx) % W, y = oy % W;
    s.add(y * W + x);
  }
  return s;
}

console.log('=== POSITIVE CONTROLS (detector must pass these) ===\n');

// --- Life glider: moves (1,1) per 4 generations
{
  const W = 64;
  const res = runLifeAnalysis(lifeGlider(W, 10, 10), W, 40, 4, 4);
  console.log('Life glider (isolated), expected ~ "3x3:p4:v(1,1)"');
  const ent = [...res.entries()].sort((a, b) => b[1] - a[1]);
  if (ent.length === 0) console.log('  *** FAILED: no mover detected ***');
  else for (const [k, v] of ent.slice(0, 5)) console.log(`  ${k.padEnd(24)} ${v}`);
  console.log('');
}

// --- Life blinker: oscillator, must NOT be reported as a particle
{
  const W = 64;
  const res = runLifeAnalysis(lifeBlinker(W, 10, 10), W, 20, 3, 4);
  console.log('Life blinker (isolated oscillator), expected: NO movers');
  if (res.size === 0) console.log('  PASS: no mover reported');
  else console.log('  *** FAILED: reported ' + [...res.keys()].join(', ') + ' ***');
  console.log('');
}

// --- rule 110 glider, seeded explicitly
{
  // a known rule-110 glider embedded in the ether background
  const W = 200, T = 120;
  const ether = [0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0];
  const init = new Uint8Array(W);
  for (let i = 0; i < W; i++) init[i] = ether[i % ether.length];
  // one of the standard glider triplets
  const g = [1, 1, 1, 0, 1, 1, 1];
  for (let j = 0; j < g.length; j++) init[(100 + j) % W] ^= g[j];
  const hits = scanECA(110, init, W, T, 8, 3);
  console.log('rule 110 with seeded glider, expected some "tile8:p*:v*" hits');
  const ent = [...hits.entries()].sort((a, b) => b[1] - a[1]);
  if (ent.length === 0) console.log('  *** FAILED: no drift detected ***');
  else for (const [k, v] of ent.slice(0, 8)) console.log(`  ${k.padEnd(24)} ${v}`);
  console.log('');
}

// --- rule 30 random (negative control)
{
  const W = 300, T = 120;
  let s = 12345 >>> 0;
  const rnd = () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
  const init = new Uint8Array(W);
  for (let i = 0; i < W; i++) init[i] = rnd() < 0.5 ? 1 : 0;
  const hits = scanECA(30, init, W, T, 8, 3);
  console.log('rule 30 random (negative control), expected: ~nothing persistent');
  const ent = [...hits.entries()].sort((a, b) => b[1] - a[1]);
  if (ent.length === 0) console.log('  PASS: no persistent drift');
  else for (const [k, v] of ent.slice(0, 5)) console.log(`  ${k.padEnd(24)} ${v}`);
  console.log('');
}
