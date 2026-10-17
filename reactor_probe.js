'use strict';
/*
 * reactor_probe.js
 * ---------------------------------------------------------------------------
 * STEP 5: is h1(6,1) a self-sustaining particle SOURCE, or just a few particles
 * bouncing around a closed torus?
 *
 * Context: the previous round found that most two-glider collisions either
 * annihilate completely or freeze into still lifes/oscillators, but a handful
 * produce movers. One configuration, h1(6,1), kept 4-7 movers and 22-38
 * components across every time window from t=200 to t=1200 with no sign of decay.
 *
 * That is suspicious for a reason worth stating precisely: the incoming pair
 * supplies exactly TWO gliders. If components keep being present at a rate far
 * above 2, something must be CREATING them. But there is a competing explanation
 * that must be ruled out: on a closed torus the two original gliders (plus any
 * debris they shed) could simply keep re-colliding forever without any net
 * creation.
 *
 * DISCRIMINATOR (total live-cell mass):
 *   - a closed re-collision system conserves or slowly DECAYS in mass
 *     (every collision of gliders with still lifes tends to consume structure)
 *   - a genuine source INCREASES mass, because it outputs new gliders
 * So we watch mass over a long run. Rising mass = creation. Falling = decay.
 * ---------------------------------------------------------------------------
 */

let SIZE = 600;
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

function followComponent(seed, maxSteps) {
  const frames = [];
  let s = seed;
  for (let t = 0; t < maxSteps; t++) {
    const cs = components(s);
    if (cs.length !== 1) return { ok: false, reason: 'split', t };
    frames.push(canon(cs[0]));
    if (s.size === 0) return { ok: false, reason: 'died', t };
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
    const vx = Math.round(dx), vy = Math.round(dy);
    return { ok: true, p, vx, vy, mover: vx !== 0 || vy !== 0, id: a.id };
  }
  return { ok: false, reason: 'no recurrence' };
}

// ============================================================ the experiment
console.log('=== is h1(6,1) a self-sustaining particle source? ===');
console.log(`torus ${SIZE}x${SIZE}\n`);

const STEPS = 6000;
const SAMPLE = 100;

let s = union(gliderAt(120, 120, 0), gliderAt(120 + 6, 120 + 1, 1));
console.log('t=0: two gliders, heading (1,1) and (-1,1), mass 10\n');
console.log('   t     mass   comps   movers   still   osc   other');
const series = [];
const t0 = Date.now();
for (let t = 1; t <= STEPS; t++) {
  s = step(s);
  if (t % SAMPLE === 0 || t === STEPS) {
    const comps = components(s);
    let movers = 0, still = 0, osc = 0, other = 0;
    for (const cells of comps) {
      const f = followComponent(new Set(cells), 60);
      if (f.ok && f.mover) movers++;
      else if (f.ok && f.p === 1) still++;
      else if (f.ok) osc++;
      else other++;
    }
    const rec = { t, mass: s.size, comps: comps.length, movers, still, osc, other };
    series.push(rec);
    console.log(`  ${String(t).padStart(5)}  ${String(rec.mass).padStart(6)}  ${String(rec.comps).padStart(6)}  ${String(movers).padStart(7)}  ${String(still).padStart(6)}  ${String(osc).padStart(5)}  ${String(other).padStart(6)}`);
  }
}
console.log(`\n(${((Date.now() - t0) / 1000).toFixed(1)}s for ${STEPS} steps)\n`);

// ---- trend analysis
const first = series.slice(0, Math.max(1, Math.floor(series.length / 4)));
const last = series.slice(-Math.max(1, Math.floor(series.length / 4)));
const avg = (arr, k) => arr.reduce((a, r) => a + r[k], 0) / arr.length;
const m1 = avg(first, 'mass'), m2 = avg(last, 'mass');
const c1 = avg(first, 'comps'), c2 = avg(last, 'comps');
const mv1 = avg(first, 'movers'), mv2 = avg(last, 'movers');

console.log('=== trend (first quarter vs last quarter) ===');
console.log(`  mass   : ${m1.toFixed(1)} -> ${m2.toFixed(1)}   (${m2 > m1 * 1.1 ? 'RISING' : m2 < m1 * 0.9 ? 'DECAYING' : 'stable'})`);
console.log(`  comps  : ${c1.toFixed(1)} -> ${c2.toFixed(1)}   (${c2 > c1 * 1.1 ? 'RISING' : c2 < c1 * 0.9 ? 'DECAYING' : 'stable'})`);
console.log(`  movers : ${mv1.toFixed(1)} -> ${mv2.toFixed(1)}   (${mv2 > mv1 * 1.1 ? 'RISING' : mv2 < mv1 * 0.9 ? 'DECAYING' : 'stable'})`);

console.log('');
if (m2 > m1 * 1.1) {
  console.log('  => MASS RISES: this is a SOURCE, not a closed re-collision loop.');
  console.log('     Something in this configuration keeps creating structure.');
} else if (m2 < m1 * 0.9) {
  console.log('  => MASS DECAYS: it is a slowly dying transient, not a steady source.');
  console.log('     The earlier "never settles" reading was just a long transient.');
} else {
  console.log('  => MASS IS STABLE: consistent with a closed system of particles');
  console.log('     re-colliding on the torus, i.e. NO net creation.');
}
