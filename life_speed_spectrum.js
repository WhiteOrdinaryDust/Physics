'use strict';
/*
 * life_speed_spectrum.js
 * ---------------------------------------------------------------------------
 * STEP 1b: make the Life "zero movers from random soup" claim solid.
 *
 * Why this is needed: in the previous run the Life detector only looked for
 * period-4 shape recurrence. That detects gliders (p=4, v=1/4) but would MISS
 * anything with a longer period. A zero result from a narrow window is weak
 * evidence. Here the window is widened to all periods 1..MAXP.
 *
 * Note on velocity limits: in Life, a pattern of period p can drift at most p
 * cells per period (light speed 1 cell/generation, radius-1 rule), and only
 * diagonal (1,1) or orthogonal (2,0) drifts are known to occur. So scanning
 * periods 1..MAXP covers every possible mover; no separate velocity cutoff is
 * needed or wanted.
 *
 * Outputs:
 *   - a histogram of the velocity spectrum of rigid localized structures
 *   - the fraction of structures that are movers vs oscillators
 *   - a known-mover sanity check (glider on a big torus, surrounded by nothing)
 * ---------------------------------------------------------------------------
 */

const W = 200;
const T = 600;
const MIN_SIZE = 4;
const MAXP = 20;

function lifeStep(set) {
  const counts = new Map();
  for (const k of set) {
    const x = k % W, y = (k - x) / W;
    for (let dy = -1; dy <= 1; dy++) {
      for (let dx = -1; dx <= 1; dx++) {
        if (dx === 0 && dy === 0) continue;
        const nk = ((y + dy + W) % W) * W + ((x + dx + W) % W);
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
      const x = k % W, y = (k - x) / W;
      for (let dy = -1; dy <= 1; dy++) {
        for (let dx = -1; dx <= 1; dx++) {
          if (dx === 0 && dy === 0) continue;
          const nk = ((y + dy + W) % W) * W + ((x + dx + W) % W);
          if (set.has(nk) && !seen.has(nk)) { seen.add(nk); Q.push(nk); }
        }
      }
    }
    comps.push(cells);
  }
  return comps;
}

function canon(cells) {
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  let sx = 0, sy = 0;
  for (const k of cells) {
    const x = k % W, y = (k - x) / W;
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
    sx += x; sy += y;
  }
  const w = maxx - minx + 1, h = maxy - miny + 1;
  let hash = 2166136261;
  for (const k of cells) {
    const x = k % W, y = (k - x) / W;
    hash = Math.imul(hash ^ (((y - miny) * w + (x - minx)) + 1), 16777619);
  }
  return { w, h, hash: hash >>> 0, n: cells.length, cx: sx / cells.length, cy: sy / cells.length };
}

function disp(ax, ay, bx, by) {
  let dx = bx - ax, dy = by - ay;
  if (dx > W / 2) dx -= W; if (dx < -W / 2) dx += W;
  if (dy > W / 2) dy -= W; if (dy < -W / 2) dy += W;
  return [dx, dy];
}

/*
 * Track shapes across time. Report, for every (shape, period) pair, the drift.
 * Classify:  v=0             -> oscillator / still life (not a particle)
 *            v != 0, |v|<=p  -> particle with velocity v/p
 */
function analyze(initSet, steps) {
  let set = initSet;
  const hist = [];        // newest last; each entry Map shapeKey -> {cx,cy,w,h,n}
  const movers = new Map();   // "wxh:v(dx,dy):p" -> count
  const stationary = new Map();
  let frames = 0;

  for (let t = 0; t < steps; t++) {
    set = lifeStep(set);
    if (set.size === 0) break;
    frames++;
    const frame = new Map();
    for (const cells of components(set)) {
      if (cells.length < MIN_SIZE) continue;
      const c = canon(cells);
      frame.set(`${c.w}:${c.h}:${c.hash}`, { cx: c.cx, cy: c.cy, w: c.w, h: c.h, n: c.n });
    }
    hist.push(frame);
    if (hist.length > MAXP + 1) hist.shift();

    for (const [key, rec] of frame) {
      for (let p = 1; p <= Math.min(MAXP, hist.length - 1); p++) {
        const old = hist[hist.length - 1 - p];
        const prev = old.get(key);
        if (!prev) continue;
        const [dx, dy] = disp(prev.cx, prev.cy, rec.cx, rec.cy);
        if (Math.abs(dx) > p || Math.abs(dy) > p) continue;   // light-speed bound
        if (dx === 0 && dy === 0) {
          if (p === 1) stationary.set(`${rec.w}x${rec.h}`, (stationary.get(`${rec.w}x${rec.h}`) || 0) + 1);
          continue;
        }
        const label = `${rec.w}x${rec.h}:p${p}:v(${dx},${dy})`;
        movers.set(label, (movers.get(label) || 0) + 1);
      }
    }
  }
  return { movers, stationary, frames };
}

function lifeGlider(ox, oy) {
  const s = new Set();
  for (const [dx, dy] of [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]]) {
    s.add(((oy + dy) % W) * W + ((ox + dx) % W));
  }
  return s;
}
function lifeLWSS(ox, oy) {
  // lightweight spaceship: travels 2 cells horizontally per 4 generations
  const s = new Set();
  for (const [dx, dy] of [[1, 0], [4, 0], [0, 1], [0, 2], [4, 2], [0, 3], [1, 3], [2, 3], [3, 3]]) {
    s.add(((oy + dy) % W) * W + ((ox + dx) % W));
  }
  return s;
}

console.log('=== Life velocity spectrum of rigid localized structures ===');
console.log(`grid ${W}x${W} torus, ${T} steps, min component size ${MIN_SIZE}, periods 1..${MAXP}\n`);

// ---- sanity checks on isolated known movers
console.log('--- sanity: isolated known movers (must be found) ---');
for (const [name, mk, expect] of [['glider', lifeGlider, 'p4:v(1,1)'], ['LWSS', lifeLWSS, 'p4:v(2,0)']]) {
  const { movers } = analyze(mk(50, 50), 60);
  const ent = [...movers.entries()].sort((a, b) => b[1] - a[1]).slice(0, 4);
  console.log(`  ${name} (expect ~${expect}):`);
  if (ent.length === 0) console.log('    *** NOT FOUND ***');
  else for (const [k, v] of ent) console.log(`    ${k.padEnd(24)} ${v}`);
}
console.log('');

// ---- random soups
console.log('--- random soups (density 0.3) ---');
const seeds = [1, 7, 42, 1234, 99, 2024];
let totalMoverLabels = 0, totalStationaryLabels = 0;
for (const seed of seeds) {
  let s = seed >>> 0;
  const rnd = () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
  const init = new Set();
  for (let i = 0; i < W * W; i++) if (rnd() < 0.3) init.add(i);
  const t0 = Date.now();
  const { movers, stationary, frames } = analyze(init, T);
  const mv = [...movers.entries()].sort((a, b) => b[1] - a[1]);
  const st = [...stationary.entries()].sort((a, b) => b[1] - a[1]);
  totalMoverLabels += mv.length;
  totalStationaryLabels += st.length;
  console.log(`  seed ${String(seed).padStart(5)} (${((Date.now() - t0) / 1000).toFixed(1)}s, ${frames} frames)`);
  console.log(`    movers: ${mv.length} distinct (shape,period,velocity) signatures`);
  for (const [k, v] of mv.slice(0, 6)) console.log(`      ${k.padEnd(26)} ${v}`);
  console.log(`    stationary: ${st.length} distinct shapes, top: ${st.slice(0, 5).map(([k, v]) => `${k}(${v})`).join(' ')}`);
}

console.log(`\n=== summary ===`);
console.log(`total distinct mover signatures across ${seeds.length} random soups: ${totalMoverLabels}`);
console.log(`total distinct stationary signatures:                        ${totalStationaryLabels}`);
console.log('Interpretation: a mover signature is a (shape, period, velocity) triple that');
console.log('recurred with rigid translation. Zero movers means the random soup produced');
console.log('no translating localized structures at all, only still lifes and oscillators.');
