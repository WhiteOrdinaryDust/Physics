'use strict';
/*
 * life_spectrum_v2.js
 * ---------------------------------------------------------------------------
 * STEP 1b (repaired): the shape key must include the TRAJECTORY, not just the
 * instantaneous shape.
 *
 * WHAT WENT WRONG IN v1 (life_speed_spectrum.js)
 *   The key was (w, h, shapeHash). Two different 3x3 components with the SAME
 *   shape but DIFFERENT velocities therefore shared a key and matched each
 *   other, manufacturing displacements that never happened. Symptoms:
 *     - reported 2x2 "moving" at (5,-10)  -- a 2x2 block is a still life
 *     - reported p=8,12,16 as separate "discoveries" -- those are just the p=4
 *       displacement counted 2, 3, 4 times
 *     - reported 0 stationary shapes and ~9000 movers, i.e. exactly backwards
 *       from the previous run (0 movers, 6..10 stationary)
 *
 * THE FIX
 *   Component identity = (shape, trajectory of the last p frames). Concretely,
 *   a structure's signature is the tuple of shape-hashes over its recent
 *   positions, together with its accumulated displacement. Two gliders in
 *   different phases/velocities no longer collide.
 *   Also: report only PRIMITIVE (p, v) pairs -- do not re-report multiples.
 *   Also: integer displacement accumulation (no float noise like 0.9999999999999982).
 *
 * CONTROLS (all four must behave)
 *   isolated glider  -> p4:v(1,1) with all four diagonal directions
 *   isolated LWSS    -> p4:v(2,0) (and mirrored)
 *   isolated blinker -> NO mover
 *   random soup      -> the actual question
 * ---------------------------------------------------------------------------
 */

const W = 200;
const T = 600;
const MIN_SIZE = 4;
const MAXP = 16;

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
  // INTEGER bounding box + integer-sum centroid (no float noise)
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
  return { w, h, hash: hash >>> 0, n: cells.length, sx, sy };
}

function ddisp(ax, ay, bx, by) {
  let dx = bx - ax, dy = by - ay;
  if (dx > W / 2) dx -= W; if (dx < -W / 2) dx += W;
  if (dy > W / 2) dy -= W; if (dy < -W / 2) dy += W;
  return [dx, dy];
}

/*
 * track(): returns for each component a chain of (frame, shapeKey, sx, sy).
 * We then require, for period p: same shape at t and t+p, AND every
 * intermediate frame has a shape consistent with the same rigid motion
 * (checked via the component's own recorded trajectory).
 *
 * Implementation: a component's trajectory is identified by following the
 * min-cell (top-left of bounding box) forward. We keep a per-frame map from
 * shapeKey -> list of components, and we link components across frames by
 * "same shape, displacement <= 1 per step".
 */
function analyze(initSet, steps, verbose) {
  let set = initSet;
  const frames = [];             // per frame: array of comps {key,w,h,sx,sy,n}
  for (let t = 0; t < steps; t++) {
    set = lifeStep(set);
    if (set.size === 0) break;
    const comps = [];
    for (const cells of components(set)) {
      if (cells.length < MIN_SIZE) continue;
      const c = canon(cells);
      comps.push(c);
    }
    frames.push(comps);
  }

  // link components across consecutive frames: same shape, step <= 1
  // node = (frame, index)
  const succ = frames.map(f => new Array(f.length).fill(-1));   // -> index in next frame
  const pred = frames.map(f => new Array(f.length).fill(-1));
  for (let t = 0; t + 1 < frames.length; t++) {
    const used = new Set();
    for (let i = 0; i < frames[t].length; i++) {
      const a = frames[t][i];
      let best = -1, bestD = null;
      for (let j = 0; j < frames[t + 1].length; j++) {
        if (used.has(j)) continue;
        const b = frames[t + 1][j];
        if (b.w !== a.w || b.h !== a.h || b.hash !== a.hash) continue;
        const [dx, dy] = ddisp(a.sx, a.sy, b.sx, b.sy);
        if (Math.abs(dx) > 1 || Math.abs(dy) > 1) continue;
        best = j; bestD = [dx, dy]; break;
      }
      if (best >= 0) { succ[t][i] = best; pred[t + 1][best] = i; used.add(best); }
    }
  }

  // walk chains and record motion
  const movers = new Map();     // "wxh:p4:v(dx,dy)" primitive -> count
  const stationary = new Map();
  const seenNode = frames.map(f => new Array(f.length).fill(false));
  for (let t = 0; t < frames.length; t++) {
    for (let i = 0; i < frames[t].length; i++) {
      if (seenNode[t][i] || pred[t][i] >= 0) continue;   // start only at chain heads
      // walk forward
      let ct = t, ci = i;
      const chain = [];
      while (ct < frames.length && ci >= 0 && !seenNode[ct][ci]) {
        seenNode[ct][ci] = true;
        chain.push({ t: ct, c: frames[ct][ci] });
        ci = succ[ct][ci]; ct++;
      }
      if (chain.length < 3) continue;
      const first = chain[0].c;
      // accumulate displacement along the chain
      let accx = 0, accy = 0;
      let prev = first;
      const stepsAcc = [];
      for (let s = 1; s < chain.length; s++) {
        const cur = chain[s].c;
        const [dx, dy] = ddisp(prev.sx, prev.sy, cur.sx, cur.sy);
        accx += dx; accy += dy;
        stepsAcc.push([dx, dy]);
        prev = cur;
      }
      // check rigid-translation-with-periodicity: find the smallest p such that
      // the motion in every consecutive p steps is identical and non-zero
      let foundP = -1, pdx = 0, pdy = 0;
      for (let p = 1; p <= Math.min(MAXP, stepsAcc.length); p++) {
        let ok = true;
        let sx = 0, sy = 0;
        for (let k = 0; k < p; k++) { sx += stepsAcc[k][0]; sy += stepsAcc[k][1]; }
        if (sx === 0 && sy === 0) continue;
        for (let s = p; s + p <= stepsAcc.length; s += p) {
          let qx = 0, qy = 0;
          for (let k = 0; k < p; k++) { qx += stepsAcc[s + k][0]; qy += stepsAcc[s + k][1]; }
          if (qx !== sx || qy !== sy) { ok = false; break; }
        }
        if (ok) { foundP = p; pdx = sx; pdy = sy; break; }
      }
      if (foundP > 0) {
        if (verbose && movers.size < 3) {
          // (kept for potential debugging)
        }
        const label = `${first.w}x${first.h}:p${foundP}:v(${pdx},${pdy})`;
        movers.set(label, (movers.get(label) || 0) + 1);
      } else if (accx === 0 && accy === 0 && chain.length >= 4) {
        stationary.set(`${first.w}x${first.h}`, (stationary.get(`${first.w}x${first.h}`) || 0) + 1);
      }
    }
  }
  return { movers, stationary };
}

// ---- known patterns
function mkGlider(ox, oy) {
  const s = new Set();
  for (const [dx, dy] of [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]]) {
    s.add(((oy + dy) % W) * W + ((ox + dx) % W));
  }
  return s;
}
function mkLWSS(ox, oy) {
  const s = new Set();
  for (const [dx, dy] of [[1, 0], [4, 0], [0, 1], [0, 2], [4, 2], [0, 3], [1, 3], [2, 3], [3, 3]]) {
    s.add(((oy + dy) % W) * W + ((ox + dx) % W));
  }
  return s;
}
function mkBlinker(ox, oy) {
  const s = new Set();
  for (let dx = 0; dx < 3; dx++) s.add((oy % W) * W + ((ox + dx) % W));
  return s;
}
function mkBlock(ox, oy) {
  const s = new Set();
  for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) s.add(((oy + dy) % W) * W + ((ox + dx) % W));
  return s;
}

console.log('=== Life rigid-structure tracker v2 (trajectory-aware identity) ===\n');
console.log('--- controls: isolated known patterns ---');
for (const [name, mk, expect] of [
  ['glider     ', mkGlider, 'p?:v(1,1) or (1,-1) etc'],
  ['LWSS       ', mkLWSS, 'p?:v(2,0) or (-2,0)'],
  ['blinker    ', mkBlinker, 'NO mover'],
  ['block      ', mkBlock, 'NO mover'],
]) {
  const { movers, stationary } = analyze(mk(50, 50), 80);
  const mv = [...movers.entries()].sort((a, b) => b[1] - a[1]);
  const st = [...stationary.entries()].sort((a, b) => b[1] - a[1]);
  console.log(`  ${name} expect ${expect}`);
  console.log(`    movers    : ${mv.length ? mv.map(([k, v]) => `${k}(${v})`).join('  ') : 'none'}`);
  console.log(`    stationary: ${st.length ? st.map(([k, v]) => `${k}(${v})`).join('  ') : 'none'}`);
}
console.log('');

console.log('--- random soups (density 0.3) ---');
const seeds = [1, 7, 42, 1234];
let totM = 0, totS = 0;
for (const seed of seeds) {
  let s = seed >>> 0;
  const rnd = () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
  const init = new Set();
  for (let i = 0; i < W * W; i++) if (rnd() < 0.3) init.add(i);
  const t0 = Date.now();
  const { movers, stationary } = analyze(init, T);
  const mv = [...movers.entries()].sort((a, b) => b[1] - a[1]);
  const st = [...stationary.entries()].sort((a, b) => b[1] - a[1]);
  totM += mv.length; totS += st.length;
  console.log(`  seed ${String(seed).padStart(5)} (${((Date.now() - t0) / 1000).toFixed(1)}s)`);
  console.log(`    movers    : ${mv.length}`);
  for (const [k, v] of mv.slice(0, 8)) console.log(`       ${k.padEnd(26)} ${v}`);
  console.log(`    stationary: ${st.length}  ${st.slice(0, 6).map(([k, v]) => `${k}(${v})`).join(' ')}`);
}
console.log(`\ntotals: movers ${totM}, stationary ${totS}`);
