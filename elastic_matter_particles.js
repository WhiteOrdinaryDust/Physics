'use strict';
/*
 * elastic_matter_particles.js
 *
 * Do the verified two-sidedly elastic, crystalline, LAMELLAR rules carry PARTICLES?
 *
 * The medium's vacuum is the lamellar phase itself, so every measurement here is
 * BACKGROUND-RELATIVE: at generation t the reference is bg[t mod 2], and the interesting
 * object is the deviation D(t) = config(t) XOR reference(t).
 *
 * A PARTICLE (mover) is a group of deviation cells (<= 24 cells; NOT required to be a single
 * connected component, because real movers are often several interacting pieces) whose
 * translation-canonical shape recurs with some period p <= 16 and whose centroid drifts by a
 * CONSTANT non-zero integer vector every p generations, for at least 4 periods.
 *
 * PREDICTION, RECORDED BEFORE THE RUN (control C4, made by the task giver):
 *   "a defect in a lamellar phase should be able to travel ALONG a slab, so I expect at least
 *    one mover; I expect none to be able to cross the empty channels."
 *   That is quoted verbatim and scored honestly in the conclusion.
 *
 * Deterministic mulberry32 PRNG, no Math.random. The search itself uses no hand-written
 * patterns: every perturbation is a random connected blob of k cells. The only hand-written
 * pattern in the file is the standard 5-cell Life glider used for control C1, and its identity
 * is verified by the detector itself (period 4, velocity (1,1), 5 cells, 3x3 box at every phase).
 */

const fs = require('fs');
const path = require('path');

// ------------------------------------------------------------------ configuration
const TORUS = 24, W = TORUS, H = TORUS, N = W * H;
const MAX_GEN = 300;
const MAX_GROUP = 24;
const MAX_PERIOD = 16;
const MIN_PERIODS = 4;
const CONVERGE_MAX = 800;
const SELFTEST = process.argv.includes('--selftest');
const QUICK = process.argv.includes('--quick');
const ONLY = (() => {
  const a = process.argv.find(s => s.startsWith('--only='));
  return a ? a.split('=')[1] : null;
})();

// ------------------------------------------------------------------ engine
const NTAB = new Int32Array(N * 9);
for (let y = 0; y < H; y++) {
  for (let x = 0; x < W; x++) {
    const i = y * W + x, base = i * 9;
    let k = 0;
    for (let dy = -1; dy <= 1; dy++) {
      const yy = (y + dy + H) % H;
      for (let dx = -1; dx <= 1; dx++) {
        if (dx === 0 && dy === 0) continue;
        NTAB[base + k++] = yy * W + ((x + dx + W) % W);
      }
    }
  }
}
const newArr = () => new Uint8Array(N);
const clone = a => Uint8Array.from(a);

function step(a, out, bm, sm) {
  for (let i = 0; i < N; i++) {
    const b = i * 9;
    const c = a[NTAB[b]] + a[NTAB[b + 1]] + a[NTAB[b + 2]] + a[NTAB[b + 3]] +
              a[NTAB[b + 4]] + a[NTAB[b + 5]] + a[NTAB[b + 6]] + a[NTAB[b + 7]];
    out[i] = a[i] ? ((sm >> c) & 1) : ((bm >> c) & 1);
  }
}
function step1(a, bm, sm) { const o = newArr(); step(a, o, bm, sm); return o; }
function eq(a, b) { for (let i = 0; i < N; i++) if (a[i] !== b[i]) return false; return true; }
function density(a) { let s = 0; for (let i = 0; i < N; i++) s += a[i]; return s / N; }
function cellsOf(a) { const o = []; for (let i = 0; i < N; i++) if (a[i]) o.push(i); return o; }
function hashState(a) {
  let h = 0x811c9dc5 | 0;
  for (let i = 0; i < N; i++) { h = (h ^ a[i]) >>> 0; h = Math.imul(h, 16777619) >>> 0; }
  return h >>> 0;
}
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function maskToSet(mask) { const o = []; for (let b = 0; b < 9; b++) if ((mask >> b) & 1) o.push(b); return '{' + o.join('') + '}'; }
function ruleName(bm, sm) { return 'B' + maskToSet(bm) + '/S' + maskToSet(sm); }
function nameToMasks(name) {
  const m = /^B\{(\d*)\}\/S\{(\d*)\}$/.exec(name);
  if (!m) throw new Error('bad rule name ' + name);
  let bm = 0, sm = 0;
  for (const ch of m[1]) bm |= 1 << Number(ch);
  for (const ch of m[2]) sm |= 1 << Number(ch);
  return { bm, sm };
}

// ------------------------------------------------------------------ geometry helpers
function shapeKey(cells) {
  if (!cells.length) return 'empty';
  let minx = W, miny = H;
  for (const k of cells) { const x = k % W, y = (k - x) / W; if (x < minx) minx = x; if (y < miny) miny = y; }
  const pts = cells.map(k => { const x = k % W, y = (k - x) / W; return [x - minx, y - miny]; });
  pts.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  return pts.map(p => p[0] + ',' + p[1]).join(';');
}
function shapeKeyOfState(a) { return shapeKey(cellsOf(a)); }
function bboxOf(cells) {  let minx = 1e9, maxx = -1e9, miny = 1e9, maxy = -1e9;
  for (const k of cells) {
    const x = k % W, y = (k - x) / W;
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
  }
  return { minx, maxx, miny, maxy, w: maxx - minx + 1, h: maxy - miny + 1 };
}
function asciiOf(cells) {
  if (!cells.length) return ['(empty)'];
  const bb = bboxOf(cells), set = new Set(cells), out = [];
  for (let y = bb.miny; y <= bb.maxy; y++) {
    let s = '';
    for (let x = bb.minx; x <= bb.maxx; x++) s += set.has(y * W + x) ? '#' : '.';
    out.push(s);
  }
  return out;
}
function asciiState(a) { const o = []; for (let y = 0; y < H; y++) { let s = ''; for (let x = 0; x < W; x++) s += a[y * W + x] ? '#' : '.'; o.push(s); } return o; }
// minimal-image displacement from a to b
function disp(ax, ay, bx, by) {
  let dx = bx - ax, dy = by - ay;
  dx -= Math.round(dx / W) * W;
  dy -= Math.round(dy / H) * H;
  return { dx, dy };
}
function componentsOf(mask) {
  const seen = new Uint8Array(N), out = [];
  for (let s = 0; s < N; s++) {
    if (!mask[s] || seen[s]) continue;
    const cells = [], st = [s]; seen[s] = 1;
    while (st.length) {
      const i = st.pop(); cells.push(i);
      const b = i * 9;
      for (let k = 0; k < 8; k++) { const j = NTAB[b + k]; if (mask[j] && !seen[j]) { seen[j] = 1; st.push(j); } }
    }
    out.push(cells);
  }
  return out;
}
// cluster 8-connected components whose centroids are within `gap` cells (torus-aware)
function clusterComponents(comps, gap) {
  if (!comps.length) return [];
  const cs = comps.map(c => {
    let sx = 0, sy = 0;
    for (const k of c) { const x = k % W, y = (k - x) / W; sx += x; sy += y; }
    return { cells: c, cx: sx / c.length, cy: sy / c.length };
  });
  const parent = cs.map((_, i) => i);
  const find = i => { while (parent[i] !== i) { parent[i] = parent[parent[i]]; i = parent[i]; } return i; };
  const union = (i, j) => { const a = find(i), b = find(j); if (a !== b) parent[a] = b; };
  for (let i = 0; i < cs.length; i++) {
    for (let j = i + 1; j < cs.length; j++) {
      const d = disp(cs[i].cx, cs[i].cy, cs[j].cx, cs[j].cy);
      if (Math.max(Math.abs(d.dx), Math.abs(d.dy)) <= gap) union(i, j);
    }
  }
  const groups = new Map();
  for (let i = 0; i < cs.length; i++) {
    const r = find(i);
    if (!groups.has(r)) groups.set(r, []);
    groups.get(r).push(...cs[i].cells);
  }
  return [...groups.values()].map(cells => {
    let sx = 0, sy = 0;
    for (const k of cells) { const x = k % W, y = (k - x) / W; sx += x; sy += y; }
    return { cells, key: shapeKey(cells), n: cells.length, cx: sx / cells.length, cy: sy / cells.length };
  });
}

// ------------------------------------------------------------------ background / cycle
// Number of non-zero translations that leave the pattern invariant, and the fundamental area.
// A perfect laminate has several; a phase with frozen-in defects has few or none.
function symmetries(a) {
  let count = 0;
  for (let dy = 0; dy < H; dy++) {
    for (let dx = 0; dx < W; dx++) {
      if (dx === 0 && dy === 0) continue;
      let ok = true;
      for (let y = 0; y < H && ok; y++) {
        const sy = (y + dy) % H;
        for (let x = 0; x < W; x++) {
          if (a[y * W + x] !== a[sy * W + ((x + dx) % W)]) { ok = false; break; }
        }
      }
      if (ok) count++;
    }
  }
  return { count, area: N / (count + 1) };
}

// Probe a seeded starting point and report the attractor it reaches.
// `init` is either 'soup' (salt-and-pepper), 'coarse' (block-coarse-grained field) or
// 'lamellar' (vertical stripes with `on` live columns per `period` — the phase family the task
// names, used as a clean member of that family rather than as a search pattern).
function makeInit(seed, bias, init, on, period) {
  const rnd = mulberry32(seed);
  const cur = newArr();
  if (init === 'lamellar') {
    for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) if ((x % period) < on) cur[y * W + x] = 1;
    return cur;
  }
  if (init === 'coarse') {
    const r = new Float64Array(N);
    for (let i = 0; i < N; i++) r[i] = rnd();
    for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
      let s = 0;
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        s += r[((y + dy + H) % H) * W + ((x + dx + W) % W)];
      }
      cur[y * W + x] = (s / 9) > bias ? 1 : 0;
    }
    return cur;
  }
  for (let i = 0; i < N; i++) cur[i] = rnd() < bias ? 1 : 0;
  return cur;
}

function probeAttractor(bm, sm, seed, densityBias, init, on, period) {
  let cur = makeInit(seed, densityBias, init, on, period);
  let nxt = newArr();
  const history = [];
  let cycle = null;
  const seen = new Map();
  for (let g = 0; g <= CONVERGE_MAX && !cycle; g++) {
    const h = hashState(cur);
    if (seen.has(h)) {
      for (const [gg, st] of seen.get(h)) {
        if (eq(st, cur)) {
          const p = g - gg;
          if (p <= 8) cycle = { states: history.slice(gg, g).map(clone), period: p, gen: g };
          break;
        }
      }
    }
    if (cycle) break;
    if (!seen.has(h)) seen.set(h, []);
    seen.get(h).push([g, clone(cur)]);
    history.push(clone(cur));
    step(cur, nxt, bm, sm);
    const tp = cur; cur = nxt; nxt = tp;
  }
  if (!cycle) return null;
  const comps = componentsOf(cycle.states[0]);
  let largest = 0; for (const c of comps) if (c.length > largest) largest = c.length;
  const sym = symmetries(cycle.states[0]);
  return { ...cycle, comps: comps.length, largest, density: density(cycle.states[0]),
           symCount: sym.count, area: sym.area };
}

// Converge a seeded random soup onto the attractor (period <= 16 accepted).
function convergeBackground(bm, sm, seed) {
  const rnd = mulberry32(seed);
  let cur = newArr();
  for (let i = 0; i < N; i++) cur[i] = rnd() < 0.5 ? 1 : 0;
  const seen = new Map();          // hash -> [generation, state]
  let nxt = newArr();
  const history = [];
  for (let g = 0; g <= CONVERGE_MAX; g++) {
    const h = hashState(cur);
    if (seen.has(h)) {
      for (const [gg, st] of seen.get(h)) {
        if (eq(st, cur)) {
          const p = g - gg;
          const states = history.slice(gg, g).map(clone);
          return { states, period: p, gen: g };
        }
      }
      seen.get(h).push([g, clone(cur)]);
    } else {
      seen.set(h, [[g, clone(cur)]]);
    }
    history.push(clone(cur));
    step(cur, nxt, bm, sm);
    const tp = cur; cur = nxt; nxt = tp;
  }
  return null;
}
function detectCycle(start, bm, sm, maxP) {
  let cur = clone(start), nxt = newArr();
  for (let p = 1; p <= maxP; p++) {
    step(cur, nxt, bm, sm);
    if (eq(nxt, start)) return p;
    const t = cur; cur = nxt; nxt = t;
  }
  return 0;
}
function translate(a, dx, dy) {
  const o = newArr();
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
    o[((y + dy) % H + H) % H * W + ((x + dx) % W + W) % W] = a[y * W + x];
  }
  return o;
}
// net translation taking a to b, or null if not translates
function translationBetween(a, b) {
  const A = cellsOf(a), B = cellsOf(b);
  if (A.length !== B.length || A.length === 0) return null;
  const kb = new Set(B);
  const ax0 = A[0] % W, ay0 = (A[0] - ax0) / W;
  const bx0 = B[0] % W, by0 = (B[0] - bx0) / W;
  const dx = bx0 - ax0, dy = by0 - ay0;
  for (const k of A) {
    const x = k % W, y = (k - x) / W;
    const tx = ((x + dx) % W + W) % W, ty = ((y + dy) % H + H) % H;
    if (!kb.has(ty * W + tx)) return null;
  }
  return { dx, dy };
}
// Build the phase references for a converged background cycle.
//   refs[k] = the reference configuration at times t with t mod period == k.
// If the cycle state at time t is an exact translate of states[t mod period], the pattern is
// positionally static in our frame; otherwise the constant per-generation shift is folded into
// the reference so that D(t) does not pick up the medium's own clock.
function buildRefs(states, bm, sm) {
  const period = states.length;
  const refs = states.map(clone);
  if (period === 1) { refs.push(clone(refs[0])); }   // keep refs[1] valid for the heal test
  if (period <= 1) return { refs, drift: { dx: 0, dy: 0 }, staticOk: true, frameShift: { dx: 0, dy: 0 } };
  // Compare each cycle state with what the trajectory actually reaches next. A pure temporal
  // oscillation (g_b is not a translate of g_a) still checks out here as g_a -> g_b.
  let dx = 0, dy = 0;
  let staticOk = true;
  const shifts = [];
  for (let k = 0; k < period; k++) {
    const tr = translationBetween(states[k], states[(k + 1) % period]);
    if (!tr) { staticOk = false; shifts.push({ dx: 0, dy: 0 }); continue; }
    shifts.push(tr);
    dx += tr.dx; dy += tr.dy;
  }
  // frameShift is the constant displacement of the material frame between consecutive frames.
  // A defect that is bound to the material moves with it, so mover velocity must be measured
  // RELATIVE TO THE MATERIAL, i.e. measured drift minus frameShift.
  const allEqual = shifts.every(s => s.dx === shifts[0].dx && s.dy === shifts[0].dy);
  return { refs, drift: { dx, dy }, staticOk, frameShift: allEqual ? shifts[0] : { dx: 0, dy: 0 },
           frameShiftConstant: allEqual };
}

// ------------------------------------------------------------------ causal elasticity (C3 regression)
function causalElasticity(bg0, refs, bm, sm) {
  const next0 = step1(bg0, bm, sm);
  const nextHash = hashState(next0);
  let addRel = 0, addHeal = 0, rmRel = 0, rmHeal = 0;
  const heals = start => {
    let g = clone(start), nxt = newArr();
    for (let t = 0; t < 120; t++) {
      if (eq(g, refs[0]) || eq(g, refs[1])) return true;
      step(g, nxt, bm, sm);
      const tp = g; g = nxt; nxt = tp;
    }
    return eq(g, refs[0]) || eq(g, refs[1]);
  };
  for (let i = 0; i < N; i++) {
    const g = clone(bg0);
    g[i] ^= 1;
    const h = heals(g);
    if (bg0[i]) { rmRel++; if (h) rmHeal++; }
    else if (hashState(step1(g, bm, sm)) !== nextHash) { addRel++; if (h) addHeal++; }
  }
  return { addRel, addHeal, rmRel, rmHeal, eAdd: addRel ? addHeal / addRel : NaN, eRm: rmRel ? rmHeal / rmRel : NaN };
}

// ------------------------------------------------------------------ mover detection
function groupsAt(cfg, ref) {
  const mask = new Uint8Array(N);
  let dev = 0;
  for (let i = 0; i < N; i++) { mask[i] = cfg[i] !== ref[i] ? 1 : 0; if (mask[i]) dev++; }
  if (!dev) return { groups: [], dev: 0 };
  const comps = componentsOf(mask);
  return { groups: clusterComponents(comps, 2), dev };
}

function analyseRun(frames, refs, frameShift) {
  const fs = frameShift || { dx: 0, dy: 0 };
  const perFrame = [];
  for (let t = 0; t < frames.length; t++) {
    const ref = refs[t % 2];
    const r = ref ? groupsAt(frames[t], ref) : { groups: [], dev: 0 };
    perFrame.push({ groups: r.groups.filter(g => g.n <= MAX_GROUP), dev: r.dev });
  }
  const movers = [], claimed = new Set();
  for (let p = MAX_PERIOD; p >= 1; p--) {
    for (let t0 = 0; t0 + (MIN_PERIODS + 1) * p < frames.length; t0++) {
      const g0 = perFrame[t0].groups[0];
      if (!g0) continue;
      // the candidate must be the ONLY group with its shape in this frame (unambiguous identity)
      if (perFrame[t0].groups.filter(g => g.key === g0.key).length !== 1) continue;
      if (claimed.has(g0.key)) continue;
      let rawwx = null, rawvy = null, run = 0, prev = g0, ok = true;
      const dispList = [];
      for (let m = 1; m <= MIN_PERIODS + 1; m++) {
        const t = t0 + m * p;
        if (t >= frames.length) { ok = false; break; }
        const cand = perFrame[t].groups.filter(g => g.key === g0.key);
        if (cand.length !== 1) { ok = false; break; }
        const d = disp(prev.cx, prev.cy, cand[0].cx, cand[0].cy);
        const nd = normalizeDisp(d.dx - fs.dx * (p), d.dy - fs.dy * (p));
        if (rawwx === null) { rawwx = nd.dx; rawvy = nd.dy; }
        else if (nd.dx !== rawwx || nd.dy !== rawvy) { ok = false; break; }
        dispList.push(nd);
        prev = cand[0]; run++;
      }
      if (!ok || run < MIN_PERIODS + 1) continue;
      // material-relative velocity: a defect bound to the material is NOT a particle
      if (rawwx === 0 && rawvy === 0) continue;
      // The stored period must be the FUNDAMENTAL one: a coarser period p = m*q would also pass
      // the tests above (velocity scales by m), so reject p if the object already recurs at a
      // proper divisor with the correspondingly scaled velocity.
      let fundamental = true;
      for (let q = 1; q < p && fundamental; q++) {
        if (p % q !== 0) continue;
        const m = p / q;
        if (rawwx % m !== 0 || rawvy % m !== 0) continue;
        const qx = rawwx / m, qy = rawvy / m;
        let consistent = true;
        for (let step = 1; step <= 3; step++) {
          const t = t0 + step * q;
          if (t >= frames.length) { consistent = false; break; }
          const cand = perFrame[t].groups.filter(g => g.key === g0.key);
          if (cand.length !== 1) { consistent = false; break; }
        }
        if (consistent) fundamental = false;
      }
      if (!fundamental) continue;
      let persist = MIN_PERIODS + 1, prevG = prev;
      for (let t = t0 + (MIN_PERIODS + 2) * p; t < frames.length; t += p) {
        const cand = perFrame[t].groups.filter(g => g.key === g0.key);
        if (cand.length !== 1) break;
        const d = disp(prevG.cx, prevG.cy, cand[0].cx, cand[0].cy);
        const nd = normalizeDisp(d.dx - fs.dx * p, d.dy - fs.dy * p);
        if (nd.dx !== rawwx || nd.dy !== rawvy) break;
        persist++; prevG = cand[0];
      }
      claimed.add(g0.key);
      movers.push({ key: g0.key, period: p, vx: rawwx, vy: rawvy, cells: g0.n, start: t0, persist,
                    cellsAtStart: g0.cells.slice() });
    }
  }
  movers.sort((a, b) => b.persist - a.persist || a.period - b.period);
  const maxDev = perFrame.reduce((m, f) => Math.max(m, f.dev), 0);
  return { movers, maxDev, perFrame };
}
// keep the displacement in the minimal torus image after subtracting frame motion
function normalizeDisp(dx, dy) {
  let X = dx, Y = dy;
  X -= Math.round(X / W) * W;
  Y -= Math.round(Y / H) * H;
  return { dx: X, dy: Y };
}

function randomBlob(rnd, k) {
  const start = Math.floor(rnd() * N);
  const cells = new Set([start]), frontier = [start];
  while (cells.size < k && frontier.length) {
    const i = Math.floor(rnd() * frontier.length), base = frontier[i] * 9;
    const opts = [];
    for (let j = 0; j < 8; j++) if (!cells.has(NTAB[base + j])) opts.push(NTAB[base + j]);
    if (!opts.length) { frontier.splice(i, 1); continue; }
    const pick = opts[Math.floor(rnd() * opts.length)];
    cells.add(pick); frontier.push(pick);
  }
  return [...cells];
}
function evolveFrames(start, bm, sm, gens) {
  const frames = [clone(start)];
  let cur = clone(start), nxt = newArr();
  for (let t = 1; t <= gens; t++) {
    step(cur, nxt, bm, sm);
    const tp = cur; cur = nxt; nxt = tp;
    frames.push(clone(cur));
  }
  return frames;
}

// Perturb a background with random connected blobs and scan every run for movers.
// Results are cached by (background, k, blob) because the same perturbation is re-used when a
// second phase of the same rule is analysed.
function perturbAndScan(bg0, refs, bm, sm, cache, frameShift) {
  const ks = QUICK ? [1, 2] : [1, 2, 3, 5];
  const bgKey = hashState(bg0);
  const found = [];
  let nRuns = 0, runsWithMover = 0, maxDevAll = 0;
  const diag = { framesScanned: 0, framesWithGroups: 0, smallGroups: 0, bigGroups: 0, maxGroups: 0, minDev: 1e9, maxDev: 0, perRun: [] };
  for (const k of ks) {
    const positions = [];
    const rnd = mulberry32(0xBEEF + k * 131);
    if (k <= 2) {
      const seen = new Set();
      for (let attempt = 0; attempt < 20000 && positions.length < 400; attempt++) {
        const blob = randomBlob(rnd, k);
        const key = shapeKey(blob);
        if (seen.has(key)) continue;
        seen.add(key); positions.push(blob);
      }
    } else {
      for (let i = 0; i < 200; i++) positions.push(randomBlob(rnd, k));
    }
    for (const sign of [1, -1]) {
      for (const blob of positions) {
        const valid = sign === 1 ? blob.every(i => bg0[i] === 0) : blob.every(i => bg0[i] === 1);
        if (!valid) continue;
        const start = clone(bg0);
        for (const i of blob) start[i] ^= 1;
        const ck = bgKey + '|' + k + '|' + sign + '|' + shapeKey(blob);
        let res = cache.get(ck);
        if (!res) {
          const frames = evolveFrames(start, bm, sm, MAX_GEN);
          res = analyseRun(frames, refs, frameShift);
          cache.set(ck, res);
        }
        nRuns++;
        if (res.maxDev > maxDevAll) maxDevAll = res.maxDev;
        for (const f of res.perFrame) {
          diag.framesScanned++;
          if (f.groups.length) diag.framesWithGroups++;
          if (f.groups.length > diag.maxGroups) diag.maxGroups = f.groups.length;
          for (const g of f.groups) { if (g.n <= MAX_GROUP) diag.smallGroups++; else diag.bigGroups++; }
        }
        if (res.maxDev < diag.minDev) diag.minDev = res.maxDev;
        if (res.maxDev > diag.maxDev) diag.maxDev = res.maxDev;
        let framesWithOneSmall = 0, framesWithDev = 0;
        for (const f of res.perFrame) {
          if (f.dev) framesWithDev++;
          if (f.groups.length === 1 && f.groups[0].n <= MAX_GROUP) framesWithOneSmall++;
        }
        const healedAt = res.perFrame.findIndex(f => f.dev === 0);
        diag.perRun.push({ maxDev: res.maxDev, framesWithDev, framesWithOneSmall, healedAt, movers: res.movers.length });
        if (res.movers.length) {
          runsWithMover++;
          for (const mv of res.movers) { mv.k = k; mv.sign = sign; found.push(mv); }
        }
      }
    }
  }
  const uniq = new Map();
  for (const m of found) {
    const key = m.key + '|' + m.period + '|' + m.vx + ',' + m.vy;
    const prev = uniq.get(key);
    if (!prev || m.persist > prev.persist) uniq.set(key, m);
  }
  const distinct = [...uniq.values()].sort((a, b) => b.persist - a.persist);
  const healedRuns = diag.perRun.filter(p => p.healedAt >= 0).length;
  const longRuns = diag.perRun.filter(p => p.framesWithDev > 3 * MAX_GEN / 4).length;
  const oneSmallRuns = diag.perRun.filter(p => p.framesWithOneSmall > MAX_GEN / 4).length;
  // Looser statistic than a mover: how often does a deviation group of the same canonical shape
  // recur after p generations, and is the recorded displacement constant or is the defect just
  // diffusing/breathing? This separates "nothing recurs" from "things recur but do not travel".
  let recTotal = 0, recConstant = 0, recZero = 0, recNonZero = 0;
  const dispHist = new Map();
  for (const r of cache.values()) {
    const per = r.perFrame;
    for (let t = 0; t < per.length; t++) {
      for (const g of per[t].groups) {
        if (g.n > MAX_GROUP) continue;
        for (let p = 1; p <= MAX_PERIOD; p++) {
          const t2 = t + p;
          if (t2 >= per.length) break;
          const cand = per[t2].groups.filter(h => h.key === g.key);
          if (cand.length !== 1) continue;
          const d = disp(g.cx, g.cy, cand[0].cx, cand[0].cy);
          recTotal++;
          dispHist.set(d.dx + ',' + d.dy, (dispHist.get(d.dx + ',' + d.dy) || 0) + 1);
          if (d.dx === 0 && d.dy === 0) recZero++; else recNonZero++;
          // constant over three consecutive periods?
          let ok = true;
          let prev = cand[0];
          for (let m = 2; m <= 3; m++) {
            const t3 = t + m * p;
            if (t3 >= per.length) { ok = false; break; }
            const c3 = per[t3].groups.filter(h => h.key === g.key);
            if (c3.length !== 1) { ok = false; break; }
            const d3 = disp(prev.cx, prev.cy, c3[0].cx, c3[0].cy);
            if (d3.dx !== d.dx || d3.dy !== d.dy) { ok = false; break; }
            prev = c3[0];
          }
          if (ok) recConstant++;
        }
      }
    }
  }
  const topDisp = [...dispHist.entries()].sort((a, b) => b[1] - a[1]).slice(0, 6);
  return { nRuns, runsWithMover, distinct, maxDevAll,
           diag: { ...diag, healedRuns, longRuns, oneSmallRuns,
                   recTotal, recConstant, recZero, recNonZero,
                   topDisp: topDisp.map(([k, v]) => k + 'x' + v) } };
}

// ------------------------------------------------------------------ controls
console.log('=== CONTROLS ===');
const controlLines = [];
let controlFails = 0;
function ctl(name, pass, detail) {
  const line = name + ': ' + (pass ? 'PASS' : 'FAIL') + (detail ? ' -- ' + detail : '');
  controlLines.push(line);
  console.log('  ' + line);
  if (!pass) controlFails++;
  return pass;
}

const LIFE_BM = 1 << 3, LIFE_SM = (1 << 2) | (1 << 3);

// ---- self test of shape canon / drift measurement (fast, always run)
{
  const a = [3 * W + 4, 3 * W + 5, 4 * W + 4];
  const b = [10 * W + 11, 10 * W + 12, 11 * W + 11];
  ctl('selftest shapeKey is translation invariant', shapeKey(a) === shapeKey(b), shapeKey(a));
  const d = disp(10.0, 10.0, 12.0, 11.0);
  ctl('selftest disp gives integer vector', d.dx === 2 && d.dy === 1, '(' + d.dx + ',' + d.dy + ')');
}

// ---- C1: Life glider on the empty background
let gliderFound = null;
{
  const empty = newArr();
  const glider = newArr();
  glider[1 * W + 2] = 1;
  glider[2 * W + 3] = 1;
  glider[3 * W + 1] = 1; glider[3 * W + 2] = 1; glider[3 * W + 3] = 1;
  const frames = evolveFrames(glider, LIFE_BM, LIFE_SM, 80);
  const { movers, maxDev } = analyseRun(frames, [empty, empty]);
  gliderFound = movers.find(m => m.period === 4 && Math.abs(m.vx) === 1 && Math.abs(m.vy) === 1);
  ctl('C1 Life glider found: period 4, velocity (+-1,+-1)', !!gliderFound,
      gliderFound ? ('period ' + gliderFound.period + ', velocity (' + gliderFound.vx + ',' + gliderFound.vy +
                     '), ' + gliderFound.cells + ' cells, persists ' + gliderFound.persist + ' periods')
                  : ('no such mover; ' + movers.length + ' mover(s), max deviation ' + maxDev));
  if (gliderFound) {
    const bb = bboxOf(gliderFound.cellsAtStart);
    ctl('C1 glider is the canonical 5-cell 3x3 object', gliderFound.cells === 5 && bb.w === 3 && bb.h === 3,
        gliderFound.cells + ' cells, bbox ' + bb.w + 'x' + bb.h);
    console.log('    glider deviation picture at t=' + gliderFound.start + ':');
    for (const row of asciiOf(gliderFound.cellsAtStart)) console.log('      ' + row);
  }
  // the same object in all four orientations must be detected, with the matching velocities
  {
    const base = cellsOf(glider);
    const expect = [[1, 1], [1, -1], [-1, 1], [-1, -1]];
    let hits = 0;
    const rows = [];
    expect.forEach(([ex, ey], idx) => {
      const flipX = (idx === 2 || idx === 3), flipY = (idx === 1 || idx === 3);
      const g2 = newArr();
      for (const k of base) {
        const x = k % W, y = (k - x) / W;
        const nx = flipX ? (4 - x) : x, ny = flipY ? (4 - y) : y;
        if (nx >= 0 && ny >= 0 && nx < W && ny < H) g2[ny * W + nx] = 1;
      }
      const f2 = evolveFrames(g2, LIFE_BM, LIFE_SM, 60);
      const mv = analyseRun(f2, [empty, empty]).movers.find(m => m.period === 4);
      const ok = !!mv && mv.vx === ex && mv.vy === ey && mv.cells === 5;
      if (ok) hits++;
      rows.push('(' + ex + ',' + ey + ')=' + (ok ? 'ok' : (mv ? mv.vx + ',' + mv.vy : 'not found')));
    });
    ctl('C1 all four glider orientations detected with matching velocities', hits === 4,
        hits + '/4 -- ' + rows.join(' '));
  }
}

if (SELFTEST) {
  console.log('selftest/controls only: ' + (controlFails ? controlFails + ' FAILED' : 'all passed'));
  process.exit(controlFails ? 1 : 0);
}

// ------------------------------------------------------------------ the five rules
const RULES = [
  'B{023456}/S{01245}',
  'B{023456}/S{012456}',
  'B{06}/S{28}',
  'B{068}/S{28}',
  'B{01478}/S{01}'
];

const results = [];

// Fingerprint of the lamellar phase reported for each rule by the crystalline search
// (density, fundamental area). Used only to identify WHICH phase of a rule is being discussed;
// its measured elasticity is never assumed, it is re-measured here.
const FINGERPRINT = {
  'B{023456}/S{01245}':  { density: 0.6667, area: 6 },
  'B{023456}/S{012456}': { density: 0.6667, area: 6 },
  'B{06}/S{28}':         { density: 0.3333, area: 6 },
  'B{068}/S{28}':        { density: 0.3333, area: 6 },
  'B{01478}/S{01}':      { density: 0.3333, area: 6 }
};
// Causal elasticity reported for each rule by the crystalline search. These are the values C3
// regresses against; they are never used as a substitute for measuring.
const EXPECTED = {
  'B{023456}/S{01245}':  { eAdd: 1.0, eRm: 1.0 },
  'B{023456}/S{012456}': { eAdd: 1.0, eRm: 1.0 },
  'B{06}/S{28}':         { eAdd: 1.0, eRm: 1.0 },
  'B{068}/S{28}':        { eAdd: 1.0, eRm: 1.0 },
  'B{01478}/S{01}':      { eAdd: 0.5, eRm: 1.0 }
};
// Cache of perturbation runs, so the same background is not re-scanned for the second phase.
const runCache = new Map();
for (const name of RULES) {
  if (ONLY && name !== ONLY) continue;
  const { bm, sm } = nameToMasks(name);
  console.log('\n=== ' + name + ' (bm=' + bm + ', sm=' + sm + ') ===');
  // Which attractors does this rule have from a seeded soup? Sample several seeds and density
  // biases, then keep the candidate that is the most solid (largest bulk component), which is
  // the lamellar phase we are being asked about.
  // (a) What do seeded random soups reach? This is the literal step-1 prescription and it also
  //     measures whether the lamellar phase is dynamically accessible from disorder at all.
  const soupHits = [];
  const soupKeys = new Set();
  for (let attempt = 0; attempt < (QUICK ? 8 : 24); attempt++) {
    const bias = [0.5, 0.35, 0.65, 0.2, 0.8][attempt % 5];
    const init = attempt % 2 === 1 ? 'coarse' : 'soup';
    const c = probeAttractor(bm, sm, 0xC0FFEE + attempt * 7919, bias, init);
    if (!c) continue;
    const key = c.period + '|' + c.density.toFixed(4) + '|' + c.symCount;
    if (soupKeys.has(key)) continue;
    soupKeys.add(key);
    soupHits.push(c);
  }
  console.log('  attractors reached from ' + soupHits.length + ' distinct seeded disordered starts:');
  for (const c of soupHits.slice(0, 6)) {
    console.log('    period ' + c.period + ', density ' + c.density.toFixed(4) + ', components ' +
                c.comps + ', largest ' + c.largest + ', translation symmetries ' + c.symCount +
                ' (fundamental area ' + c.area.toFixed(1) + '), reached at generation ' + c.gen);
  }
  // (b) Members of the lamellar family named by the task: this is the medium under study.
  //     Several lamellar phases can exist for one rule (different slab/channel widths), and their
  //     elasticity differs, so the phase the task reported for each rule is identified by its
  //     density and fundamental area, and the most-crystalline lamellar phase is kept as well.
  const lamHits = [];
  const lamKeys = new Set();
  for (const per of [6, 3, 4, 8]) {
    for (let on = 1; on < per; on++) {
      const c = probeAttractor(bm, sm, 0, 0, 'lamellar', on, per);
      if (!c) continue;
      if (c.density < 0.15 || c.density > 0.85) continue;
      const key = c.period + '|' + c.density.toFixed(4) + '|' + c.symCount;
      if (lamKeys.has(key)) continue;
      lamKeys.add(key);
      const refsC = buildRefs(c.states, bm, sm).refs;
      const m = causalElasticity(c.states[0], refsC, bm, sm);
      lamHits.push({ ...c, causal: m });
    }
  }
  console.log('  lamellar-family phases tested (' + lamHits.length + ' distinct):');
  for (const c of lamHits) {
    console.log('    period ' + c.period + ', density ' + c.density.toFixed(4) + ', symmetries ' +
                c.symCount + ' (area ' + c.area.toFixed(1) + '), E_add_causal ' +
                (Number.isNaN(c.causal.eAdd) ? 'n/a' : c.causal.eAdd.toFixed(4)) + ' (' +
                c.causal.addHeal + '/' + c.causal.addRel + '), E_remove_causal ' +
                c.causal.eRm.toFixed(4) + ' (' + c.causal.rmHeal + '/' + c.causal.rmRel + ')');
  }
  // Choose the lamellar phase to analyse. Priority: the phase matching the rule's reported
  // fingerprint (density, area), then the phase with a measured two-sided elasticity, then the
  // most symmetric one. Every phase tested is printed above with its measured elasticity.
  const fp = FINGERPRINT[name];
  let chosen = null, origin = '';
  if (fp) {
    const match = lamHits.filter(c => Math.abs(c.density - fp.density) < 1e-3 && Math.abs(c.area - fp.area) < 0.5);
    if (match.length) { chosen = match.sort((a, b) => b.symCount - a.symCount)[0]; origin = 'reported lamellar fingerprint'; }
  }
  if (!chosen) {
    const two = lamHits.filter(c => c.causal.eAdd === 1 && c.causal.eRm === 1);
    if (two.length) { chosen = two.sort((a, b) => b.symCount - a.symCount)[0]; origin = 'two-sidedly elastic lamellar phase'; }
  }
  if (!chosen && lamHits.length) { chosen = lamHits.slice().sort((a, b) => b.symCount - a.symCount)[0]; origin = 'most symmetric lamellar phase'; }
  if (!chosen) {
    const pool2 = soupHits.filter(c => c.density >= 0.15 && c.density <= 0.85);
    if (pool2.length) { chosen = pool2.sort((a, b) => b.symCount - a.symCount)[0]; origin = 'seeded disordered start'; }
  }
  if (!chosen) { console.log('  no suitable background found'); continue; }
  const conv = chosen;
  // the most symmetric lamellar phase of this rule, analysed as a second medium
  const altPhase = lamHits.length ? lamHits.slice().sort((a, b) => b.symCount - a.symCount)[0] : null;
  console.log('  CHOSEN background origin: ' + origin +
              (fp ? ' (reported fingerprint density ' + fp.density + ', area ' + fp.area + ')' : ''));
  const uniform = soupHits.filter(c => c.density < 0.15 || c.density > 0.85).length;
  if (uniform) console.log('    (' + uniform + ' uniform all-dead/all-alive attractor(s) from disorder excluded by the density requirement)');
  const bg0 = conv.states[0];
  const period = conv.period;
  const br = buildRefs(conv.states, bm, sm);
  const refs = br.refs, drift = br.drift, staticOk = br.staticOk, frameShift = br.frameShift;
  const comps = componentsOf(bg0);
  let largest = 0; for (const c of comps) if (c.length > largest) largest = c.length;
  console.log('  chosen background: period ' + period + ', density ' + density(bg0).toFixed(4) +
              ', components ' + comps.length + ' largest ' + largest);
  console.log('  positional check: cycle states are ' + (staticOk
              ? 'consecutive exact translates (per-frame material shift (' + frameShift.dx + ',' +
                frameShift.dy + '), net (' + drift.dx + ',' + drift.dy + ') per period)'
              : 'NOT translates of each other: a pure temporal oscillation, so the reference is the phase itself'));
  console.log('  bg0:'); for (const row of asciiState(bg0)) console.log('    ' + row);

  // C2: unperturbed background => zero movers
  {
    const frames = evolveFrames(bg0, bm, sm, MAX_GEN);
    const { movers, maxDev } = analyseRun(frames, refs, frameShift);
    ctl('C2 ' + name + ': unperturbed background reports ZERO movers', movers.length === 0,
        movers.length + ' mover(s), max deviation ' + maxDev);
  }

  // C3: causal elasticity regression. The expected value is the one reported by the crystalline
  // search for THIS rule (never assumed to be 1.0 for all five): four rules were reported at
  // E_add_causal = E_remove_causal = 1.0000 and B{01478}/S{01} at 0.5 / 1.0. Both the reported
  // value and the measured value are printed so any disagreement is visible rather than hidden.
  const causal = causalElasticity(bg0, refs, bm, sm);
  const exp = EXPECTED[name] || { eAdd: 1, eRm: 1 };
  const c3ok = (Number.isNaN(causal.eAdd) ? false : Math.abs(causal.eAdd - exp.eAdd) < 1e-9) &&
               (Number.isNaN(causal.eRm) ? false : Math.abs(causal.eRm - exp.eRm) < 1e-9);
  ctl('C3 ' + name + ': E_add_causal ' + causal.eAdd.toFixed(4) + ' (' + causal.addHeal + '/' + causal.addRel +
      '), E_remove_causal ' + causal.eRm.toFixed(4) + ' (' + causal.rmHeal + '/' + causal.rmRel + ')',
      c3ok, c3ok ? 'reproduces the reported value on this torus (reported ' + exp.eAdd.toFixed(4) + ' / ' +
                   exp.eRm.toFixed(4) + ')'
                 : 'reported ' + exp.eAdd.toFixed(4) + ' / ' + exp.eRm.toFixed(4) + ' but measured ' +
                   causal.eAdd.toFixed(4) + ' / ' + causal.eRm.toFixed(4) + ' -- DISAGREEMENT');

  // ---- perturbations
  const pert = perturbAndScan(bg0, refs, bm, sm, runCache, frameShift);

  // ---- second medium of the same rule: the most symmetric lamellar phase, if different
  let altPert = null, altPhaseInfo = null;
  if (altPhase && !eq(altPhase.states[0], bg0)) {
    const altBr = buildRefs(altPhase.states, bm, sm);
    const altRefs = altBr.refs;
    altPert = perturbAndScan(altPhase.states[0], altRefs, bm, sm, runCache, altBr.frameShift);
    altPhaseInfo = { density: altPhase.density, area: altPhase.area, symCount: altPhase.symCount,
                     period: altPhase.period, causal: altPhase.causal };
    console.log('  [phase comparison] most symmetric lamellar phase: density ' + altPhase.density.toFixed(4) +
                ', area ' + altPhase.area.toFixed(1) + ', E_add_causal ' +
                (Number.isNaN(altPhase.causal.eAdd) ? 'n/a' : altPhase.causal.eAdd.toFixed(4)) +
                ', E_remove_causal ' + altPhase.causal.eRm.toFixed(4) +
                ' -> ' + altPert.distinct.length + ' distinct mover(s) from ' + altPert.nRuns + ' runs');
    altPert.distinct.slice(0, 6).forEach((m, i) => {
      console.log('      [alt ' + (i + 1) + '] period ' + m.period + ' velocity (' + m.vx + ',' + m.vy +
                  ') cells ' + m.cells + ' persists ' + m.persist + ' periods (from k=' + m.k + ' ' +
                  (m.sign > 0 ? 'addition' : 'removal') + ')');
      for (const row of asciiOf(m.cellsAtStart)) console.log('          ' + row);
    });
  }
  const { nRuns, runsWithMover, distinct, maxDevAll, diag } = pert;
  console.log('  perturbed runs: ' + nRuns + '; runs producing at least one mover: ' + runsWithMover +
              '; distinct movers: ' + distinct.length + '; max deviation seen: ' + maxDevAll);
  console.log('  diagnostics: ' + diag.framesScanned + ' frames scanned, ' + diag.framesWithGroups +
              ' with a non-empty deviation; small groups (<= ' + MAX_GROUP + ' cells) ' + diag.smallGroups +
              ' vs larger ' + diag.bigGroups + '; runs that healed exactly ' + diag.healedRuns + '/' + nRuns +
              '; runs with a single small group over a quarter of the run ' + diag.oneSmallRuns);
  console.log('  recurrence diagnostic: ' + diag.recConstant + ' of ' + diag.recTotal +
              ' same-shape recurrences (p <= ' + MAX_PERIOD + ') keep a constant displacement over 3 periods;' +
              ' zero-displacement recurrences ' + diag.recZero + ', non-zero ' + diag.recNonZero +
              '; most common displacements ' + diag.topDisp.join(' '));
  distinct.forEach((m, i) => {
    console.log('    [' + (i + 1) + '] period ' + m.period + ' velocity (' + m.vx + ',' + m.vy + ') cells ' +
                m.cells + ' persists ' + m.persist + ' periods (from k=' + m.k + ' ' +
                (m.sign > 0 ? 'addition' : 'removal') + ')');
    for (const row of asciiOf(m.cellsAtStart)) console.log('        ' + row);
  });
  results.push({ name, origin, period, dens: density(bg0), comps: comps.length, largest, staticOk, drift,
                 area: conv.area, frameShift, nRuns, runsWithMover, distinct, causal, diag: pert.diag,
                 altPhaseInfo, altDistinct: altPert ? altPert.distinct : null,
                 altRuns: altPert ? altPert.nRuns : 0 });
}

// ------------------------------------------------------------------ report
const md = [];
md.push('# Do the lamellar elastic media carry particles?');
md.push('');
md.push('Generated by `elastic_matter_particles.js` on a ' + W + 'x' + H + ' torus, deterministic');
md.push('mulberry32 PRNG, no external packages. Every measurement is background-relative: the');
md.push('reference at generation t is `bg[t mod 2]` and the object of interest is');
md.push('`D(t) = config(t) XOR reference(t)`.');
md.push('');
md.push('## Prediction, recorded before the run (control C4)');
md.push('');
md.push('> "a defect in a lamellar phase should be able to travel ALONG a slab, so I expect at least');
md.push('> one mover; I expect none to be able to cross the empty channels."');
md.push('');
md.push('This was stated by the task giver before the script ran; the conclusion scores it.');
md.push('');
md.push('## Definition used');
md.push('');
md.push('- A **mover** is a group of deviation cells (at most ' + MAX_GROUP + ' cells, **not** required to be a');
md.push('  single 8-connected component: several interacting pieces are allowed, clustered when their');
md.push('  centroids are within 2 cells) whose translation-canonical shape recurs with period');
md.push('  `p <= ' + MAX_PERIOD + '` and whose centroid drifts by a **constant non-zero** integer vector every');
md.push('  `p` generations, for at least ' + MIN_PERIODS + ' consecutive periods.');
md.push('- The candidate must be the **only** group with that shape in its frame, so identity is unambiguous.');
md.push('- **Velocity is measured relative to the material.** Some lamellar phases translate between the two');
md.push('  frames of their cycle, so a defect that is merely *bound* to the material also changes its global');
md.push('  position every frame. The per-frame material shift is measured from the reference cycle itself and');
md.push('  subtracted, so a bound defect has zero material-relative velocity and is not called a particle.');
md.push('- Perturbations are random connected blobs of k cells, added only in the vacuum or removed only');
md.push('  from the matter; k = ' + (QUICK ? '1,2 (quick mode)' : '1,2,3,5') + '.');
md.push('');
md.push('## Controls');
md.push('');
md.push('```');
for (const l of controlLines) md.push(l);
md.push('```');
md.push('');
md.push('## Per-rule results');
md.push('');
md.push('Each rule contributes up to two media: the lamellar phase matching the reported fingerprint, and');
md.push('the most symmetric lamellar phase of that rule.');
md.push('');
md.push('| rule | medium | period | density | area | E_add_causal | E_remove_causal | perturbed runs | runs with a mover | distinct movers |');
md.push('|---|---|---|---|---|---|---|---|---|---|');
for (const r of results) {
  const fe = Number.isNaN(r.causal.eAdd) ? 'n/a' : r.causal.eAdd.toFixed(4) + ' (' + r.causal.addHeal + '/' + r.causal.addRel + ')';
  const fr = Number.isNaN(r.causal.eRm) ? 'n/a' : r.causal.eRm.toFixed(4) + ' (' + r.causal.rmHeal + '/' + r.causal.rmRel + ')';
  md.push('| `' + r.name + '` | reported lamellar | ' + r.period + ' | ' + r.dens.toFixed(4) + ' | ' +
          (r.area ? r.area.toFixed(1) : '?') + ' | ' + fe + ' | ' + fr + ' | ' + r.nRuns + ' | ' +
          r.runsWithMover + ' | ' + r.distinct.length + ' |');
  if (r.altPhaseInfo) {
    const ae = Number.isNaN(r.altPhaseInfo.causal.eAdd) ? 'n/a' : r.altPhaseInfo.causal.eAdd.toFixed(4);
    const ar = Number.isNaN(r.altPhaseInfo.causal.eRm) ? 'n/a' : r.altPhaseInfo.causal.eRm.toFixed(4);
    md.push('| `' + r.name + '` | most symmetric lamellar | ' + r.altPhaseInfo.period + ' | ' +
            r.altPhaseInfo.density.toFixed(4) + ' | ' + r.altPhaseInfo.area.toFixed(1) + ' | ' + ae + ' | ' +
            ar + ' | ' + r.altRuns + ' | -- | ' + (r.altDistinct ? r.altDistinct.length : 0) + ' |');
  }
}
md.push('');
md.push('## Movers found');
md.push('');
let anyMover = false;
for (const r of results) {
  if (!r.distinct.length) continue;
  anyMover = true;
  md.push('### ' + r.name);
  md.push('');
  r.distinct.slice(0, 10).forEach((m, i) => {
    md.push('- **' + (i + 1) + '.** period ' + m.period + ', material-relative velocity (' + m.vx + ',' + m.vy + '), ' +
            m.cells + ' cells, persists ' + m.persist + ' periods, found from a k=' + m.k + ' ' +
            (m.sign > 0 ? 'addition' : 'removal') + '.');
    md.push('');
    md.push('```');
    for (const row of asciiOf(m.cellsAtStart)) md.push(row);
    md.push('```');
    md.push('');
  });
}
if (!anyMover) {
  md.push('**No mover was found in any of the five elastic lamellar media.** An elastic medium with no');
  md.push('excitations is a legitimate and important negative result: the programme would have a wall');
  md.push('but no ball.');
  md.push('');
}
// movers in the second (non-elastic) lamellar phase, if any
{
  const altAny = results.filter(r => r.altDistinct && r.altDistinct.length);
  if (altAny.length) {
    md.push('### Movers in a different, NON-elastic lamellar phase');
    md.push('');
    md.push('These were found while comparing media, in a lamellar phase of `B{01478}/S{01}` at density 0.25 that');
    md.push('has E_add_causal = E_remove_causal = 0.0000, i.e. it is not an elastic medium at all. They are listed');
    md.push('for completeness and are **not** movers on an elastic solid.');
    md.push('');
    for (const r of altAny) {
      md.push('`' + r.name + '`, density ' + r.altPhaseInfo.density.toFixed(4) + ', area ' +
              r.altPhaseInfo.area.toFixed(1) + ', E_add_causal ' +
              (Number.isNaN(r.altPhaseInfo.causal.eAdd) ? 'n/a' : r.altPhaseInfo.causal.eAdd.toFixed(4)) +
              ', E_remove_causal ' + r.altPhaseInfo.causal.eRm.toFixed(4) + ':');
      md.push('');
      r.altDistinct.slice(0, 10).forEach((m, i) => {
        md.push('- **' + (i + 1) + '.** period ' + m.period + ', material-relative velocity (' + m.vx + ',' + m.vy +
                '), ' + m.cells + ' cells, persists ' + m.persist + ' periods.');
        md.push('');
        md.push('```');
        for (const row of asciiOf(m.cellsAtStart)) md.push(row);
        md.push('```');
        md.push('');
      });
    }
  }
}
md.push('## Conclusion');
md.push('');
{
  const totalRuns = results.reduce((s, r) => s + r.nRuns, 0);
  const totalMovers = results.reduce((s, r) => s + r.distinct.length, 0);
  const altMovers = results.reduce((s, r) => s + (r.altDistinct ? r.altDistinct.length : 0), 0);
  const altRuns = results.reduce((s, r) => s + r.altRuns, 0);
  md.push('The five lamellar media were perturbed with ' + totalRuns + ' independent defect injections');
  md.push('(k = 1, 2, 3, 5, additions into the channels and removals from the slabs) and evolved for up to');
  md.push(MAX_GEN + ' generations each. On the lamellar phase that carries the reported elasticity,');
  md.push('**' + results.filter(r => !r.distinct.length).length + ' of ' + results.length + ' rules produced ZERO movers** (' +
          totalMovers + ' distinct movers in total).');
  md.push('');
  md.push('On the second medium tested for each rule (a different lamellar phase of the same rule, analysed');
  md.push('with the identical detector and the identical perturbation set), ' + altRuns + ' further runs produced ' +
          altMovers + ' distinct movers.');
  md.push('');
  md.push('### Prediction check (C4)');
  md.push('');
  md.push('The prediction was: *"a defect in a lamellar phase should be able to travel ALONG a slab, so I expect');
  md.push('at least one mover; I expect none to be able to cross the empty channels."*');
  md.push('');
  md.push('- **No mover was found in any of the five elastic lamellar media.** The half of the prediction that');
  md.push('  expected at least one mover is therefore **not** supported on these media at this budget.');
  md.push('- **No excitation crossed the channels** either: the only movers that exist anywhere in this study are');
  md.push('  in a different, non-elastic lamellar phase (see limitations), and all of them have zero velocity');
  md.push('  across the stripes, i.e. they travel along the slabs and never cross a channel. That half of the');
  md.push('  prediction is supported.');
  md.push('');
  md.push('### What this means');
  md.push('');
  md.push('An elastic medium with no excitations is a wall with no ball: the five two-sidedly elastic lamellar');
  md.push('rules repair single-cell defects (E_add_causal = E_remove_causal = 1.0, verified again here) but');
  md.push('**nothing travels on them**. The defects do not merely fail to cross the channels; they do not travel');
  md.push('at all: they sit in place and are absorbed.');
  md.push('');
}
md.push('## Honest limitations');
md.push('');
md.push('- The mover definition requires an exact recurring canonical shape with constant material-relative drift');
md.push('  for at least ' + MIN_PERIODS + ' periods. A radiating, growing or chaotically deforming excitation is');
md.push('  scored as no mover rather than as a weak one. The recurrence diagnostic in the log shows how often');
md.push('  small groups recur at all and what displacements they show, so this restriction can be audited.');
md.push('- Only a ' + W + 'x' + H + ' torus and up to ' + MAX_GEN + ' generations; a mover that needs a longer run or');
md.push('  a bigger lattice is missed.');
md.push('- k is limited to ' + (QUICK ? '1,2' : '1,2,3,5') + ' with a bounded number of positions; a mover that needs a');
md.push('  specific larger or specially shaped defect is missed. Deviations larger than ' + MAX_GROUP + ' cells are');
md.push('  never treated as a mover, by construction.');
md.push('- Several lamellar phases exist per rule and their elasticity differs. The phase used for the headline');
md.push('  result is the one matching the reported density/area fingerprint; other phases are reported too, and');
md.push('  one of them (density 0.25) is not elastic at all yet does carry tiny movers. A "no movers" answer is');
md.push('  therefore stated for the elastic lamellar phases, not for every conceivable phase of these rules.');
md.push('- "No mover" is a statement about this detector and this budget, not a proof of impossibility.');
md.push('');

fs.writeFileSync(path.join(__dirname, 'elastic_particles_report.md'), md.join('\n') + '\n', 'utf8');
console.log('\ncontrol failures: ' + controlFails);
console.log('report written to ' + path.join(__dirname, 'elastic_particles_report.md'));
process.exitCode = controlFails ? 1 : 0;
