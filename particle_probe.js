'use strict';
/*
 * particle_probe.js
 * ---------------------------------------------------------------------------
 * STEP 1: the simplest usable definition of a "particle".
 *
 * DEFINITION
 *   A particle is a finite connected localized pattern that translates rigidly
 *   (shape preserved) with constant velocity.
 *
 *   Canonical form: translate the pattern so its bounding box sits at the origin.
 *   Particle  <=>  canonical form is constant over k consecutive steps AND the
 *                  centroid displacement per step is constant.
 *   Oscillators (displacement 0) are counted separately, not as particles.
 *
 * WHY THIS DEFINITION
 *   - decidable, no substrate assumptions
 *   - it is the operational definition of "an object an internal observer could
 *     track": something whose shape survives, i.e. something with an identity
 *   - it excludes oscillators, which are 'field' rather than 'particles'
 *
 * CONTROLS
 *   rule 110 : known to have gliders            -> positive control
 *   rule 30  : chaotic                          -> negative control
 *   Life     : what emerges from random soups?  -> open question
 *
 * The question this answers:
 *   How RARE are particles, under the simplest definition, in a minimal rule?
 * ---------------------------------------------------------------------------
 */

// ============================================================ 1D elementary CA
function makeRuleTable(rule) {
  const t = new Uint8Array(8);
  for (let n = 0; n < 8; n++) t[n] = (rule >> n) & 1;
  return t;
}

function evolveECA(rule, init, steps, W) {
  const tbl = makeRuleTable(rule);
  const rows = new Uint8Array((steps + 1) * W);
  let cur = Uint8Array.from(init);
  rows.set(cur, 0);
  const nxt = new Uint8Array(W);
  for (let t = 0; t < steps; t++) {
    for (let i = 0; i < W; i++) {
      const l = cur[(i - 1 + W) % W], c = cur[i], r = cur[(i + 1) % W];
      nxt[i] = tbl[(l << 2) | (c << 1) | r];
    }
    cur.set(nxt);
    rows.set(cur, (t + 1) * W);
  }
  return rows;
}

// Split a 1D ring row into maximal runs of 1s. Returns [{start,len}].
function runs1D(row, W) {
  const out = [];
  // rotate so index 0 is a 0 if possible, to avoid splitting a wrapped run
  let z = -1;
  for (let i = 0; i < W; i++) if (row[i] === 0) { z = i; break; }
  if (z === -1) return [{ start: 0, len: W }]; // all ones
  let i = 0;
  while (i < W) {
    const idx = (z + i) % W;
    if (row[idx] === 1) {
      let len = 0;
      while (len < W && row[(z + i + len) % W] === 1) len++;
      out.push({ start: (z + i) % W, len });
      i += len;
    } else i++;
  }
  return out;
}

function scanECA(rule, seeds, W, T, minRun, minLife) {
  const found = new Map(); // key = len:dx -> count
  const runCounts = new Map();
  for (const seed of seeds) {
    let s = seed >>> 0;
    const rnd = () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
    const init = new Uint8Array(W);
    for (let i = 0; i < W; i++) init[i] = rnd() < 0.5 ? 1 : 0;
    const rows = evolveECA(rule, init, T, W);

    // track each run: (len, start). Record displacement per step.
    let prev = new Map(); // key "start:len" -> {len,start,age,dx}
    for (let t = 0; t <= T; t++) {
      const row = rows.subarray(t * W, (t + 1) * W);
      const rs = runs1D(row, W);
      const cur = new Map();
      for (const r of rs) {
        if (r.len < minRun) continue;
        // find a previous run of the same length that is adjacent (moved by dx)
        let best = null, bestDx = null;
        for (const dx of [-1, 0, 1]) {
          const k = ((r.start - dx) % W + W) % W + ':' + r.len;
          if (prev.has(k)) { best = prev.get(k); bestDx = dx; break; }
        }
        if (best && best.dx === bestDx) {
          best.age++;
          best.start = r.start;
          cur.set(r.start + ':' + r.len, best);
          if (best.age >= minLife) {
            const key = `len${r.len}:v${bestDx}`;
            found.set(key, (found.get(key) || 0) + 1);
          }
        } else {
          cur.set(r.start + ':' + r.len, { len: r.len, start: r.start, age: 1, dx: bestDx === null ? 0 : bestDx });
        }
      }
      prev = cur;
      for (const v of prev.values()) {
        const k = `len${v.len}`;
        runCounts.set(k, (runCounts.get(k) || 0) + 1);
      }
    }
  }
  return { found, runCounts };
}

// ============================================================ Life (2D)
function lifeStep(set) {
  const counts = new Map();
  for (const k of set) {
    const x = k % 4096, y = (k - x) / 4096;
    for (let dy = -1; dy <= 1; dy++) {
      for (let dx = -1; dx <= 1; dx++) {
        if (dx === 0 && dy === 0) continue;
        const nx = (x + dx + 4096) % 4096, ny = (y + dy + 4096) % 4096;
        const nk = ny * 4096 + nx;
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

// connected components (8-connected) of a set of cells
function components(set) {
  const seen = new Set();
  const comps = [];
  const Q = [];
  for (const start of set) {
    if (seen.has(start)) continue;
    Q.length = 0; Q.push(start); seen.add(start);
    const cells = [];
    while (Q.length) {
      const k = Q.pop();
      cells.push(k);
      const x = k % 4096, y = (k - x) / 4096;
      for (let dy = -1; dy <= 1; dy++) {
        for (let dx = -1; dx <= 1; dx++) {
          if (dx === 0 && dy === 0) continue;
          const nx = (x + dx + 4096) % 4096, ny = (y + dy + 4096) % 4096;
          const nk = ny * 4096 + nx;
          if (set.has(nk) && !seen.has(nk)) { seen.add(nk); Q.push(nk); }
        }
      }
    }
    comps.push(cells);
  }
  return comps;
}

// canonical form: bounding box + bitmask of occupancy, plus centroid offset
function canon(cells) {
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  let sx = 0, sy = 0;
  for (const k of cells) {
    const x = k % 4096, y = (k - x) / 4096;
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
    sx += x; sy += y;
  }
  const w = maxx - minx + 1, h = maxy - miny + 1;
  const bits = new Uint8Array(w * h);
  for (const k of cells) {
    const x = k % 4096, y = (k - x) / 4096;
    bits[(y - miny) * w + (x - minx)] = 1;
  }
  return { w, h, bits, minx, miny, cx: sx / cells.length, cy: sy / cells.length, n: cells.length };
}

function runLife(seed, W, T, minSize, minLife) {
  let s = seed >>> 0;
  const rnd = () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
  let set = new Set();
  for (let i = 0; i < W * W; i++) if (rnd() < 0.3) set.add(i);

  const found = new Map();   // "wxh:v(dx,dy)" -> shape-steps, for MOVERS only
  const stillFound = new Map();
  let prev = new Map();      // shapeKey -> {age, dx, dy, minx, miny}

  for (let t = 0; t < T; t++) {
    set = lifeStep(set);
    if (set.size === 0) break;
    const comps = components(set);
    const cur = new Map();
    for (const cells of comps) {
      if (cells.length < minSize) continue;
      const c = canon(cells);
      let hash = 0;
      for (let i = 0; i < c.bits.length; i++) hash = (hash * 31 + c.bits[i]) | 0;
      const shapeKey = `${c.w}:${c.h}:${hash}`;
      const rec = prev.get(shapeKey);
      if (rec) {
        let dx = c.minx - rec.minx, dy = c.miny - rec.miny;
        if (dx > 2048) dx -= 4096; if (dx < -2048) dx += 4096;
        if (dy > 2048) dy -= 4096; if (dy < -2048) dy += 4096;
        // a particle keeps its shape AND keeps the same non-zero velocity
        if (rec.dx === null || (rec.dx === dx && rec.dy === dy)) {
          const isMover = (dx !== 0 || dy !== 0);
          const ok = isMover;
          if (ok) {
            rec.age++;
            rec.dx = dx; rec.dy = dy;
            rec.minx = c.minx; rec.miny = c.miny;
            if (rec.age >= minLife) {
              const key = `${c.w}x${c.h}:v(${dx},${dy})`;
              found.set(key, (found.get(key) || 0) + 1);
            }
            cur.set(shapeKey, rec);
          } else {
            // stationary rigid shape = oscillator / still life, counted separately
            rec.age++;
            rec.dx = 0; rec.dy = 0;
            if (rec.age >= minLife) {
              const key = `${c.w}x${c.h}`;
              stillFound.set(key, (stillFound.get(key) || 0) + 1);
            }
            cur.set(shapeKey, rec);
          }
          continue;
        }
      }
      cur.set(shapeKey, { w: c.w, h: c.h, minx: c.minx, miny: c.miny, dx: null, dy: null, age: 1 });
    }
    prev = cur;
  }
  return { found, stillFound };
}

// ============================================================ run
const SEEDS = [1, 7, 42, 1234];

console.log('=== particle probe: simplest definition (translate rigidly, constant velocity) ===\n');

// ---- rule 110 (positive control) and rule 30 (negative control)
for (const rule of [110, 30]) {
  const W = 512, T = 1200;
  const { found, runCounts } = scanECA(rule, SEEDS, W, T, 3, 20);
  console.log(`--- ECA rule ${rule} (ring ${W}, ${T} steps, ${SEEDS.length} seeds) ---`);
  const ent = [...found.entries()].sort((a, b) => b[1] - a[1]).slice(0, 12);
  if (ent.length === 0) {
    console.log('  no translating persistent runs found');
  } else {
    console.log('  candidate particles (runs that translate, age>=20):');
    for (const [k, v] of ent) console.log(`    ${k.padEnd(18)} occurrences ${v}`);
  }
  const top = [...runCounts.entries()].sort((a, b) => b[1] - a[1]).slice(0, 6);
  console.log('  most common run lengths: ' + top.map(([k, v]) => `${k}:${v}`).join('  '));
  console.log('');
}

// ---- Life
console.log('--- Life (256x256 torus, random soup density 0.3) ---');
for (const seed of SEEDS) {
  const t0 = Date.now();
  const { found, stillFound } = runLife(seed, 256, 300, 4, 12);
  const movers = [...found.entries()].sort((a, b) => b[1] - a[1]);
  const still = [...stillFound.entries()].sort((a, b) => b[1] - a[1]);
  console.log(`  seed ${String(seed).padStart(4)}  (${((Date.now() - t0) / 1000).toFixed(1)}s)`);
  console.log(`    MOVERS (rigid + constant non-zero velocity): ${movers.length} distinct`);
  for (const [k, v] of movers.slice(0, 10)) console.log(`      ${k.padEnd(22)} ${v}`);
  console.log(`    stationary rigid shapes (oscillators/still lifes): ${still.length} distinct`);
  for (const [k, v] of still.slice(0, 8)) console.log(`      ${k.padEnd(22)} ${v}`);
}
