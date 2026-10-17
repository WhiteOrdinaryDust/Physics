'use strict';
/*
 * timescale_spectrum.js
 * ---------------------------------------------------------------------------
 * Emergence of a LOW-FREQUENCY TIME LAYER.
 *
 * Framing (this replaces the previous "find a massless mode" experiment, which
 * was looking for the wrong kind of object):
 *
 *   The photon is not a stable pattern. It is TRANSMISSION ITSELF, which an
 *   observer happens to register as a particle. So the thing to look for is not a
 *   particle called "photon" but a SECOND TIMESCALE:
 *
 *       the substrate changes every step (period 1 for raw cells),
 *       but how long must an observer wait before it can register a change?
 *
 *   A particle is a pattern whose shape recurs with period P. To resolve which
 *   PHASE such a pattern is in, you must wait a full period. So P is a hard
 *   lower bound on the particle-level clock. If large regions of the world are
 *   locked into a common period P >> 1, then a low-frequency layer has emerged,
 *   and P is its time unit.
 *
 * MEASUREMENTS
 *   1. per-cell temporal period  (smallest p with x(t+p) == x(t) over a window)
 *   2. the distribution of periods
 *   3. the SIZE of the largest connected spatial domain sharing one period
 *      -- a large domain means a coherent emergent clock, not local noise
 *   4. the largest domain with period >= 2 (a genuine low-frequency layer)
 * ---------------------------------------------------------------------------
 */

// ============================================================ Life
const LW = 200;
let SIZE = LW;
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

function lifePeriods(initSet, steps, sampleStride, maxP) {
  let s = new Set(initSet);
  const pts = [];
  for (let x = 0; x < LW; x += sampleStride) for (let y = 0; y < LW; y += sampleStride) pts.push(y * LW + x);
  const nPts = pts.length;
  const hist = [];
  for (const k of pts) hist.push(new Uint8Array(steps + 1));
  for (let t = 0; t <= steps; t++) {
    for (let i = 0; i < nPts; i++) hist[i][t] = s.has(pts[i]) ? 1 : 0;
    if (t === steps) break;
    s = lifeStep(s);
  }
  // period per point: smallest p such that the tail of the series is p-periodic
  const periods = new Int32Array(nPts);
  const WINDOW = Math.min(120, steps);
  for (let i = 0; i < nPts; i++) {
    const h = hist[i];
    let best = 0;
    for (let p = 1; p <= maxP; p++) {
      let ok = true;
      for (let t = steps - WINDOW; t + p <= steps; t++) {
        if (h[t] !== h[t + p]) { ok = false; break; }
      }
      if (ok) { best = p; break; }
    }
    periods[i] = best;               // 0 = not periodic within maxP
  }
  return { pts, periods, nPts };
}

// connected domains of equal period (4-neighbour on the sampled lattice)
function domainSizes(pts, periods, stride) {
  const idx = new Map();
  for (let i = 0; i < pts.length; i++) idx.set(pts[i], i);
  const seen = new Uint8Array(pts.length);
  const sizes = [];
  const step4 = [[1, 0], [-1, 0], [0, 1], [0, -1]];
  for (let i = 0; i < pts.length; i++) {
    if (seen[i]) continue;
    const p = periods[i];
    const stack = [i];
    seen[i] = 1;
    let size = 0;
    const members = [];
    while (stack.length) {
      const a = stack.pop();
      size++;
      members.push(a);
      const k = pts[a];
      const x = k % LW, y = (k - x) / LW;
      for (const [dx, dy] of step4) {
        const nx = x + dx * stride, ny = y + dy * stride;
        if (nx < 0 || ny < 0 || nx >= LW || ny >= LW) continue;
        const j = idx.get(ny * LW + nx);
        if (j === undefined || seen[j]) continue;
        if (periods[j] !== p) continue;
        seen[j] = 1;
        stack.push(j);
      }
    }
    sizes.push({ period: p, size, members });
  }
  return sizes;
}

console.log('=== does a low-frequency TIME LAYER emerge? ===\n');

// ---- Life
{
  const STEPS = 900, STRIDE = 4, MAXP = 12;
  let seedv = 42 >>> 0;
  const rnd = () => (seedv = (seedv * 1664525 + 1013904223) >>> 0) / 4294967296;
  const init = new Set();
  for (let i = 0; i < LW * LW; i++) if (rnd() < 0.3) init.add(i);
  console.log(`--- Life ${LW}x${LW}, density 0.3, ${STEPS} steps, sampling every ${STRIDE} cells ---`);
  const t0 = Date.now();
  const { pts, periods, nPts } = lifePeriods(init, STEPS, STRIDE, MAXP);
  console.log(`  sampled ${nPts} cells (${((Date.now() - t0) / 1000).toFixed(1)}s)`);

  const counts = new Map();
  for (const p of periods) counts.set(p, (counts.get(p) || 0) + 1);
  console.log('  period histogram (0 = not periodic within ' + MAXP + '):');
  for (const [p, c] of [...counts.entries()].sort((a, b) => a[0] - b[0])) {
    const label = p === 0 ? 'non-periodic' : `period ${p}`;
    console.log(`    ${label.padEnd(16)} ${String(c).padStart(6)}  (${(100 * c / nPts).toFixed(1)}%)`);
  }

  const doms = domainSizes(pts, periods, STRIDE);
  doms.sort((a, b) => b.size - a.size);
  console.log('  largest coherent domains (cells sharing one period):');
  for (const d of doms.slice(0, 8)) {
    console.log(`    period ${String(d.period).padStart(2)}  size ${String(d.size).padStart(6)} sampled cells  (~${d.size * STRIDE * STRIDE} real cells)`);
  }
  const bigP1 = doms.filter(d => d.period > 1).sort((a, b) => b.size - a.size)[0];
  console.log(`  largest domain with period >= 2: ${bigP1 ? `period ${bigP1.period}, size ${bigP1.size} sampled (~${bigP1.size * STRIDE * STRIDE} real cells)` : 'NONE'}`);
  console.log('');
}

// ============================================================ rule 110 ether
console.log('--- rule 110: the ether background (known period 4) ---');
{
  const W = 512, T = 1500, MAXP = 12;
  const ether = [0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0];
  const init = new Uint8Array(W);
  for (let i = 0; i < W; i++) init[i] = ether[i % ether.length];
  let row = Uint8Array.from(init);
  const hist = [];
  for (let i = 0; i < W; i++) hist.push(new Uint8Array(T + 1));
  for (let t = 0; t <= T; t++) {
    for (let i = 0; i < W; i++) hist[i][t] = row[i];
    if (t === T) break;
    const out = new Uint8Array(W);
    for (let i = 0; i < W; i++) {
      const idx = (row[(i - 1 + W) % W] << 2) | (row[i] << 1) | row[(i + 1) % W];
      out[i] = (110 >> idx) & 1;
    }
    row = out;
  }
  const WINDOW = 120;
  const periods = new Int32Array(W);
  for (let i = 0; i < W; i++) {
    const h = hist[i];
    for (let p = 1; p <= MAXP; p++) {
      let ok = true;
      for (let t = T - WINDOW; t + p <= T; t++) if (h[t] !== h[t + p]) { ok = false; break; }
      if (ok) { periods[i] = p; break; }
    }
  }
  const counts = new Map();
  for (const p of periods) counts.set(p, (counts.get(p) || 0) + 1);
  console.log('  period histogram:');
  for (const [p, c] of [...counts.entries()].sort((a, b) => a[0] - b[0])) {
    console.log(`    ${(p === 0 ? 'non-periodic' : 'period ' + p).padEnd(16)} ${String(c).padStart(6)}  (${(100 * c / W).toFixed(1)}%)`);
  }
  // spatial domains on the ring
  const seen = new Uint8Array(W);
  const doms = [];
  for (let i = 0; i < W; i++) {
    if (seen[i]) continue;
    const p = periods[i];
    let j = i, n = 0;
    while (!seen[j] && periods[j] === p) { seen[j] = 1; n++; j = (j + 1) % W; }
    doms.push({ period: p, size: n });
  }
  doms.sort((a, b) => b.size - a.size);
  console.log('  largest coherent domains:');
  for (const d of doms.slice(0, 5)) console.log(`    period ${String(d.period).padStart(2)}  size ${d.size} cells`);
}

console.log('\n=== how to read this ===');
console.log('  The substrate clock is 1 step, by construction. The question is whether');
console.log('  a SECOND, slower clock emerges and organises SPACE.');
console.log('  A large domain with period P >= 2 means a region of the world that an');
console.log('  internal observer cannot see change faster than once per P steps -> that');
console.log('  is an emergent low-frequency time layer, and P is its unit.');
