'use strict';
/*
 * quench_mechanism.js
 * ---------------------------------------------------------------------------
 * STEP 9: WHY does a second collision on the same site destroy the cascade?
 *
 * Observed (previous round): two cascade islands planted at the SAME site give a
 * peak of 38 instead of the ~527 that separated planting predicts. Adding seed
 * material KILLS the eruption.
 *
 * HYPOTHESIS
 *   The cascade is sustained by ACTIVE INTERMEDIATES -- gliders that keep
 *   re-colliding. If so, then anything that destroys those intermediates should
 *   quench the cascade, and the quench should depend on WHEN the disruption
 *   arrives:
 *     - a glider arriving before the cascade gets going   -> weak effect
 *     - a glider arriving AT the peak                     -> kills it
 *     - a glider arriving after the plateau has formed     -> no effect at all
 *
 * TEST
 *   Grow the h1(6,1) cascade alone and record its mass curve (baseline).
 *   Then repeat, injecting a single glider aimed at the collision site, at a
 *   range of injection times. Compare final mass and peak against the baseline.
 *
 * DIAGNOSTIC
 *   Also measure the number of active MOVERS in the cascade at several times.
 *   If the cascade really is maintained by movers, mover count should be high
 *   through the growth phase and collapse once the plateau forms.
 * ---------------------------------------------------------------------------
 */

const SIZE = 700;
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

function canon(cells) {
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity, sx = 0, sy = 0;
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
    sx += x; sy += y;
  }
  const w = maxx - minx + 1, h = maxy - miny + 1;
  const g = [];
  for (let r = 0; r < h; r++) g.push(new Array(w).fill('.'));
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    g[y - miny][x - minx] = 'o';
  }
  return { w, h, n: cells.length, sx, sy, cx: sx / cells.length, cy: sy / cells.length, id: `${w}x${h}:` + g.map(r => r.join('')).join('/') };
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
function island(bx, by, heading, ox, oy) {
  return [gliderAt(bx, by, 0), gliderAt(bx + ox, by + oy, heading)];
}

// ---- solo mover counting, for the diagnostic
function followComponent(seed, maxSteps) {
  const frames = [];
  let s = seed;
  for (let t = 0; t < maxSteps; t++) {
    const cs = components(s);
    if (cs.length !== 1) return { ok: false };
    frames.push(canon(cs[0]));
    if (s.size === 0) return { ok: false };
    s = step(s);
  }
  for (let p = 1; p <= 15; p++) {
    let ok = true;
    for (let t = 0; t + p < frames.length; t++) if (frames[t].id !== frames[t + p].id) { ok = false; break; }
    if (!ok) continue;
    const a = frames[0], b = frames[p];
    let dx = (b.sx - a.sx) / a.n, dy = (b.sy - a.sy) / a.n;
    if (Math.abs(dx) > SIZE / 2) dx -= Math.sign(dx) * SIZE;
    if (Math.abs(dy) > SIZE / 2) dy -= Math.sign(dy) * SIZE;
    return { ok: true, p, mover: Math.round(dx) !== 0 || Math.round(dy) !== 0 };
  }
  return { ok: false };
}

const H = 1, OX = 6, OY = 1;          // h1(6,1), the strongest cascade
const BX = 250, BY = 250;
const STEPS = 2000;

console.log('=== quench mechanism: does disrupting active intermediates kill the cascade? ===');
console.log(`torus ${SIZE}x${SIZE}, island h1(6,1) at (${BX},${BY}), ${STEPS} steps\n`);

// ---- baseline
console.log('--- baseline: cascade alone ---');
const baseCells = island(BX, BY, H, OX, OY);
let s = union(...baseCells);
const baseMass = [];
const baseMovers = new Map();
for (let t = 0; t <= STEPS; t++) {
  baseMass.push(s.size);
  if ([200, 400, 600, 800, 1000, 1400, 1800].includes(t)) {
    let m = 0, tot = 0;
    for (const cells of components(s)) { tot++; if (followComponent(new Set(cells), 40).mover) m++; }
    baseMovers.set(t, { movers: m, comps: tot });
  }
  if (t === STEPS) break;
  s = step(s);
}
let basePeak = 0, basePeakT = 0;
for (let t = 0; t < baseMass.length; t++) if (baseMass[t] > basePeak) { basePeak = baseMass[t]; basePeakT = t; }
console.log(`  peak ${basePeak} @t${basePeakT}, final ${baseMass[STEPS]}`);
console.log('  t     mass   movers/comps');
for (const [t, v] of [...baseMovers.entries()].sort((a, b) => a[0] - b[0])) {
  console.log(`  ${String(t).padStart(4)}  ${String(baseMass[t]).padStart(5)}   ${v.movers}/${v.comps}`);
}

// ---- injection sweep
console.log('\n--- a single glider injected toward the collision site, at various times ---');
console.log('  injectT   final   peak   verdict');
const injectResults = [];
for (const injectT of [0, 100, 200, 400, 600, 800, 1000, 1200, 1600]) {
  // start the cascade, then at injectT add one glider heading (-1,-1) at a
  // position that brings it into the cascade region a few tens of steps later
  let cur = union(...island(BX, BY, H, OX, OY));
  for (let t = 0; t < injectT; t++) cur = step(cur);
  const gx = BX + 60, gy = BY + 60;                 // heading toward the site
  cur = union(cur, gliderAt(gx, gy, 2));            // turns=2 -> velocity (-1,-1)
  let peak = cur.size;
  for (let t = injectT; t < STEPS; t++) {
    cur = step(cur);
    if (cur.size > peak) peak = cur.size;
  }
  const final = cur.size;
  const verdict = final < baseMass[STEPS] * 0.6 ? 'QUENCHED'
    : final > baseMass[STEPS] * 1.4 ? 'ENHANCED'
      : 'little effect';
  injectResults.push({ injectT, final, peak, verdict });
  console.log(`  ${String(injectT).padStart(7)}  ${String(final).padStart(5)}   ${String(peak).padStart(5)}   ${verdict}   (baseline final ${baseMass[STEPS]})`);
}

console.log('\n=== follow-up: inject a FULL ISLAND instead of one glider ===');
console.log('  (a single glider does nothing, so the quench needs more than one extra particle)');
console.log('  injectT   final   peak   verdict   (baseline final 161, peak 319)');
const islandResults = [];
for (const injectT of [0, 100, 200, 400, 600, 800, 1200]) {
  let cur = union(...island(BX, BY, H, OX, OY));
  for (let t = 0; t < injectT; t++) cur = step(cur);
  // plant the second island's own two gliders at the same site
  cur = union(cur, ...island(BX, BY, 2, 0, 1));   // h2(0,1) overlapping
  let peak = cur.size;
  for (let t = injectT; t < STEPS; t++) {
    cur = step(cur);
    if (cur.size > peak) peak = cur.size;
  }
  const final = cur.size;
  const verdict = final < 60 ? 'QUENCHED' : final > 200 ? 'ENHANCED' : 'partial';
  islandResults.push({ injectT, final, peak, verdict });
  console.log(`  ${String(injectT).padStart(7)}  ${String(final).padStart(5)}   ${String(peak).padStart(5)}   ${verdict}`);
}

console.log('\n=== interpretation ===');
const quenchedEarly = injectResults.filter(r => r.verdict === 'QUENCHED');
console.log(`  single-glider injections that quench : ${quenchedEarly.length} of ${injectResults.length}`);
console.log(`  full-island injections that quench   : ${islandResults.filter(r => r.verdict === 'QUENCHED').length} of ${islandResults.length}`);
console.log('');
console.log('  baseline final mass    : ' + baseMass[STEPS]);
console.log('  baseline peak          : ' + basePeak + ' @t' + basePeakT);
const peakMovers = baseMovers.get(basePeakT) || baseMovers.get(800);
console.log(`  mover count near the peak (t=${basePeakT}): ${peakMovers ? peakMovers.movers + '/' + peakMovers.comps : 'n/a'}`);
console.log('');
console.log('  Reading: if a SINGLE glider never quenches at any time, the cascade is');
console.log('  not held together by fragile intermediates; one extra particle is simply');
console.log('  absorbed. If a FULL ISLAND quenches at every time, the effect is');
console.log('  structural rather than kinetic.');
