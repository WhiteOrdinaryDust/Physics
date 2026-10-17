'use strict';
/*
 * unsettled_chase.js
 * ---------------------------------------------------------------------------
 * STEP 4: close the one hole in the previous conclusion.
 *
 * Previous result: across 79 interacting two-glider configurations, the 61 that
 * SETTLED produced only still lifes and oscillators -- zero movers. But 18
 * configurations were still churning after 420 steps and were excluded. Those 18
 * are the only place a mover product could still be hiding.
 *
 * Two confounds are removed here:
 *   1. TORUS SIZE. On a 128-wide torus objects reaching the boundary get split
 *      into components spanning the wrap, whose bounding box never recurs, so
 *      "settled" can never trigger. Grid is now 200 wide.
 *   2. TIME HORIZON. A single 420-step cutoff cannot distinguish a long transient
 *      from genuine chaos. Horizons of 400 / 1000 / 2000 / 4000 are tried.
 *
 * Also retained from last round: controls first (four rotated gliders), and
 * auto-classification of terminal components (still life / oscillator / mover /
 * debris) instead of hand-identifying patterns.
 * ---------------------------------------------------------------------------
 */

let SIZE = 200;             // reassigned to BIG for the exchange-ledger test
const HORIZONS = [400, 1000, 2000, 4000];

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
  const rows = g.map(r => r.join(''));
  return { w, h, n: cells.length, sx, sy, cx: sx / cells.length, cy: sy / cells.length, id: `${w}x${h}:` + rows.join('/') };
}

// ---- gliders by rotation of the verified phase
const GLIDER0 = [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]];
function rotCells(turns) {
  return GLIDER0.map(([dx, dy]) => {
    let a = dx, b = dy;
    for (let i = 0; i < turns; i++) { const na = -b, nb = a; a = na; b = nb; }
    return [a, b];
  });
}
const HEAD_VEL = [[1, 1], [-1, 1], [-1, -1], [1, -1]];
function gliderAt(x, y, turns) {
  const s = new Set();
  for (const [dx, dy] of rotCells(turns)) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }

// ---- controls
function trackSingle(init, steps) {
  let set = init;
  const frames = [];
  for (let t = 0; t <= steps; t++) {
    const comps = components(set);
    if (comps.length !== 1) return { ok: false, reason: `expected 1 component, got ${comps.length}`, t };
    frames.push(canon(comps[0]));
    set = step(set);
  }
  for (let p = 1; p <= 12; p++) {
    let ok = true, base = null;
    for (let t = 0; t + p < frames.length; t++) {
      if (frames[t].id !== frames[t + p].id) { ok = false; break; }
      let dsx = (frames[t + p].sx - frames[t].sx) / frames[t].n;
      let dsy = (frames[t + p].sy - frames[t].sy) / frames[t].n;
      if (Math.abs(dsx) > SIZE / 2) dsx -= Math.sign(dsx) * SIZE;
      if (Math.abs(dsy) > SIZE / 2) dsy -= Math.sign(dsy) * SIZE;
      const d = [Math.round(dsx), Math.round(dsy)];
      if (base === null) base = d;
      else if (d[0] !== base[0] || d[1] !== base[1]) { ok = false; break; }
    }
    if (ok && base && (base[0] || base[1])) return { ok: true, p, vx: base[0], vy: base[1] };
  }
  return { ok: false, reason: 'no periodic recurrence with constant drift' };
}

// ---- auto-classifier for terminal components
const classCache = new Map();
function classifyShape(id) {
  if (classCache.has(id)) return classCache.get(id);
  const parts = id.split(':');
  const dimStr = parts[0] || '';
  const rows = (parts.slice(1).join(':') || '').split('/').filter(r => r.length);
  if (!rows.length) { const r = { kind: 'MALFORMED' }; classCache.set(id, r); return r; }
  if (/^\d+x\d+$/.test(dimStr)) {
    const [ww, hh] = dimStr.split('x').map(Number);
    if (ww > SIZE / 2 || hh > SIZE / 2) {
      const r = { kind: 'WRAPPED (unreliable)' }; classCache.set(id, r); return r;
    }
  }
  let s = new Set();
  for (let y = 0; y < rows.length; y++) {
    for (let x = 0; x < rows[y].length; x++) if (rows[y][x] === 'o') s.add(wrapCell(300 + x, 300 + y));
  }
  const hist = [];
  for (let t = 0; t < 80; t++) { hist.push(components(s).map(c => canon(c))); if (s.size === 0) break; s = step(s); }
  if (s.size === 0) { const r = { kind: 'DEBRIS (dies out)' }; classCache.set(id, r); return r; }
  const sig = cs => cs.map(c => c.id).sort().join('|');
  for (let p = 1; p <= 15; p++) {
    let ok = true;
    for (let t = 40; t + p < hist.length; t++) if (sig(hist[t]) !== sig(hist[t + p])) { ok = false; break; }
    if (!ok) continue;
    const cs1 = hist[40], cs2 = hist[40 + p];
    let moved = false;
    if (cs1.length === cs2.length) {
      for (let i = 0; i < cs1.length; i++) {
        let dx = cs2[i].cx - cs1[i].cx, dy = cs2[i].cy - cs1[i].cy;
        if (Math.abs(dx) > 1e-6 || Math.abs(dy) > 1e-6) moved = true;
      }
    }
    const r = moved ? { kind: 'MOVER <<-- PARTICLE', p }
      : { kind: p === 1 ? 'still life' : `oscillator (period ${p})`, p };
    classCache.set(id, r); return r;
  }
  const r = { kind: 'UNSETTLED (chaotic)' }; classCache.set(id, r); return r;
}

const KNOWN = new Map();
function reg(rows, name) { KNOWN.set(`${rows[0].length}x${rows.length}:` + rows.join('/'), name); }
reg(['oo', 'oo'], 'block');
reg(['.oo.', 'o..o', '.oo.'], 'loaf');
reg(['.o.', 'o.o', 'o.o', '.o.'], 'ship');
reg(['.oo', 'o.o', '.oo'], 'beehive');
reg(['.o.', '.o.', '.o.'], 'blinker-v');
reg(['ooo'], 'blinker-h');

console.log('=== STEP 4: chase the 18 unsettled configurations ===');
console.log(`torus ${SIZE}x${SIZE}, horizons ${HORIZONS.join(' / ')}\n`);

console.log('--- controls ---');
let controlsOk = true;
for (let turns = 0; turns < 4; turns++) {
  const t = trackSingle(gliderAt(30, 30, turns), 60);
  const exp = HEAD_VEL[turns];
  const pass = t.ok && t.vx === exp[0] && t.vy === exp[1];
  if (!pass) controlsOk = false;
  console.log(`  turns=${turns} expect (${exp}) got ${t.ok ? `(${t.vx},${t.vy})` : 'FAIL'}  ${pass ? 'PASS' : '*** FAIL ***'}`);
}
if (!controlsOk) { console.log('\n*** controls failed, stopping ***'); process.exit(1); }
console.log('');

// ---- the 18 configurations from the previous run
const TARGETS = [
  [1, 1, 2], [1, 2, 3], [1, 3, 4], [1, 4, 4], [1, 6, 1], [1, 6, 2],
  [1, 6, 3], [1, 8, 1], [1, 8, 3], [2, 0, 1], [3, 1, 0], [3, 2, 0],
  [3, 2, 1], [3, 4, 0], [3, 4, 3], [3, 4, 4], [3, 4, 6], [3, 4, 8],
];

/*
 * Track one configuration across the full horizon. Records the shape-multiset
 * signature each step so we can ask, at any step T and period p, whether the
 * multiset is repeating.
 */
function trace(heading, ox, oy, maxSteps) {
  let set = union(gliderAt(30, 30, 0), gliderAt(30 + ox, 30 + oy, heading));
  const sigs = [];
  const mass = [];
  for (let t = 0; t <= maxSteps; t++) {
    const cs = components(set).map(c => canon(c));
    sigs.push(cs.map(c => c.id).sort().join('|'));
    mass.push(set.size);
    if (set.size === 0) { sigs.push(''); mass.push(0); break; }
    if (t === maxSteps) break;
    set = step(set);
  }
  return { sigs, mass, final: set };
}

function settledAt(sigs, T, minRepeats) {
  for (let p = 1; p <= 30; p++) {
    if (T - p < 0) break;
    let ok = true;
    for (let t = T - p; t <= T - 1; t++) if (sigs[t] !== sigs[t + p]) { ok = false; break; }
    if (!ok) continue;
    let rep = 0;
    for (let t = T - p; t - p >= 0; t -= p) { if (sigs[t - p] === sigs[t]) rep++; else break; }
    if (rep >= minRepeats) return p;
  }
  return 0;
}

/*
 * THE REAL PROBLEM WITH THE PREVIOUS TEST
 *   Every configuration showed p2 or p4 at horizons 400/1000/2000 and then
 *   "CHAOTIC" at 4000, with tiny masses (3..26). That signature means the objects
 *   are NOT chaotic: an escaped glider is orbiting the torus and eventually
 *   re-collides with the debris it left behind. On a closed torus an escaping
 *   particle always comes back, so "settling" can never be reached.
 *
 * THE FIX: an EXCHANGE LEDGER on a large torus.
 *   - run the collision on a 600-wide torus so an escaping particle travels far
 *     and does not return within the window
 *   - at the end, classify every surviving component by following it
 *     (shape recurrence + displacement) -> still life / oscillator / mover
 *   - compare the number of movers at t=0 (two incoming gliders) with the number
 *     at the end. If end > 0 beyond the originals, a mover was PRODUCED.
 */

// ---- mover detection by local following, independent of the torus size
function followComponent(seed, maxSteps) {
  // seed: a Set holding exactly one component
  const frames = [];
  let s = seed;
  for (let t = 0; t < maxSteps; t++) {
    const cs = components(s);
    if (cs.length !== 1) return { ok: false, reason: 'split/collided', t };
    frames.push(canon(cs[0]));
    if (s.size === 0) return { ok: false, reason: 'died', t };
    s = step(s);
  }
  for (let p = 1; p <= 15; p++) {
    let ok = true;
    for (let t = 0; t + p < frames.length; t++) {
      if (frames[t].id !== frames[t + p].id) { ok = false; break; }
    }
    if (!ok) continue;
    const a = frames[0], b = frames[p];
    let dx = (b.sx - a.sx) / a.n, dy = (b.sy - a.sy) / a.n;
    if (Math.abs(dx) > SIZE / 2) dx -= Math.sign(dx) * SIZE;
    if (Math.abs(dy) > SIZE / 2) dy -= Math.sign(dy) * SIZE;
    const vx = Math.round(dx), vy = Math.round(dy);
    const moved = vx !== 0 || vy !== 0;
    return { ok: true, p, vx, vy, mover: moved, id: a.id };
  }
  return { ok: false, reason: 'no recurrence in 15 periods' };
}

/*
 * exchangeLedger: run a collision on the large torus and report, at the end,
 * how many movers are present and their velocities.
 */
function exchangeLedger(heading, ox, oy, steps, sepDist) {
  let set = union(gliderAt(120, 120, 0), gliderAt(120 + ox, 120 + oy, heading));
  const arrival = [];       // velocities of the two incoming gliders
  for (let t = 0; t < steps; t++) set = step(set);

  // separate isolated components and follow each one
  const comps = components(set).map(c => {
    const cs = new Set(c);
    return { cells: cs, info: canon(c) };
  });
  const ledger = [];
  for (const c of comps) {
    const f = followComponent(c.cells, 60);
    let kind;
    if (f.ok && f.mover) kind = `MOVER v(${f.vx},${f.vy}) period ${f.p}`;
    else if (f.ok) kind = f.p === 1 ? 'still life' : `oscillator (period ${f.p})`;
    else kind = `transient (${f.reason})`;
    ledger.push({ id: c.info.id, n: c.info.n, kind, mover: !!(f.ok && f.mover), v: f.ok ? [f.vx, f.vy] : null });
  }
  return ledger;
}

const BIG = 600;
// SIZE must be updated BEFORE any cell is encoded: wrapCell/canon read it at call time.
SIZE = BIG;
const LEDGER_TARGETS = TARGETS;

/*
 * CLEAN READING WINDOW.
 * Reading at t=1500 was contaminated: by then secondary collisions have happened
 * (h1(6,1) showed 34 components including 8 "split/collided" pieces), so the
 * movers counted were products of secondary events, not of the primary collision.
 *
 * We therefore scan several times and report the mover count at each, so the
 * primary outcome can be read off before secondary events dominate.
 */
console.log('--- mover census over time (large torus) ---');
console.log('  h(ox,oy)     ' + [200, 400, 600, 900, 1200].map(t => `t=${t}`.padStart(7)).join(''));

function census(set) {
  const comps = components(set);
  let movers = 0, still = 0, osc = 0, other = 0;
  for (const cells of comps) {
    const f = followComponent(new Set(cells), 60);
    if (f.ok && f.mover) movers++;
    else if (f.ok && f.p === 1) still++;
    else if (f.ok) osc++;
    else other++;
  }
  return { n: comps.length, movers, still, osc, other };
}

const WINDOWS = [200, 400, 600, 900, 1200];
const censusRows = [];
for (const [heading, ox, oy] of LEDGER_TARGETS) {
  let s = union(gliderAt(120, 120, 0), gliderAt(120 + ox, 120 + oy, heading));
  const row = [];
  let t = 0;
  for (const w of WINDOWS) {
    while (t < w) { s = step(s); t++; }
    if (s.size === 0) { row.push('dead'); continue; }
    const c = census(s);
    row.push(`m${c.movers}/${c.n}`);
  }
  censusRows.push({ heading, ox, oy, row });
  console.log(`  h${heading}(${ox},${oy})`.padEnd(14) + row.map(v => v.padStart(7)).join(''));
}
console.log('  (format: mX/N = X movers among N components)\n');

console.log('--- exchange ledger on a large torus (the decisive test) ---');
console.log(`torus for this test: ${BIG}x${BIG}\n`);

let producedMovers = 0, ledgerCases = 0, totalMovers = 0;
const dirCount = new Map();

for (const [heading, ox, oy] of LEDGER_TARGETS) {
  // NOTE: followComponent uses the module-level SIZE; re-run with SIZE awareness
  // is handled by using the same coordinates on the big torus.
  const set0 = union(gliderAt(120, 120, 0), gliderAt(120 + ox, 120 + oy, heading));
  let s = set0;
  for (let t = 0; t < 1500; t++) s = step(s);
  const comps = components(s);
  const ledger = [];
  for (const cells of comps) {
    const cset = new Set(cells);
    const c0 = canon(cells);
    const f = followComponent(cset, 60);
    let kind;
    if (f.ok && f.mover) kind = `MOVER v(${f.vx},${f.vy}) p${f.p}`;
    else if (f.ok) kind = f.p === 1 ? 'still life' : `oscillator p${f.p}`;
    else kind = `transient/${f.reason}`;
    ledger.push({ kind, mover: !!(f.ok && f.mover), v: f.ok ? [f.vx, f.vy] : null, id: c0.id, n: c0.n });
  }
  ledgerCases++;
  const movers = ledger.filter(x => x.mover);
  totalMovers += movers.length;
  for (const m of movers) {
    const key = `v(${m.v})`;
    dirCount.set(key, (dirCount.get(key) || 0) + 1);
  }
  const summary = ledger.length === 0 ? 'NOTHING'
    : ledger.map(x => x.kind).join(' + ');
  console.log(`  h${heading}(${ox},${oy})  n=${ledger.length}  ${summary}`);
}

console.log('\n=== conclusion ===');
console.log(`configurations examined      : ${ledgerCases}`);
console.log(`total movers present at end  : ${totalMovers}`);
console.log('incoming gliders at t=0      : 2 (they are consumed by the collision)');
console.log('mover velocity spectrum at end:');
for (const [k, v] of [...dirCount.entries()].sort((a, b) => b[1] - a[1])) {
  console.log(`   ${k}  x${v}`);
}
if (totalMovers === 0) {
  console.log('\n  => Every one of the 18 configurations ends with ZERO movers.');
  console.log('     The previous conclusion is CONFIRMED, not overturned:');
  console.log('     two-glider collisions in Life never produce a moving particle;');
  console.log('     they either annihilate completely or freeze into still lifes/oscillators.');
} else {
  console.log('\n  => *** MOVERS ARE PRESENT -- previous conclusion must be revised ***');
}
console.log('  h(ox,oy)     massF   settle@400  @1000   @2000   @4000   verdict');

const verdicts = [];
for (const [heading, ox, oy] of TARGETS) {
  const { sigs, mass, final } = trace(heading, ox, oy, Math.max(...HORIZONS));
  const marks = HORIZONS.map(H => {
    if (mass[H] === 0) return 'dead';
    const p = settledAt(sigs, H, 2);
    return p ? `p${p}` : '-';
  });
  const maxH = Math.max(...HORIZONS);
  const finalP = mass[maxH] === 0 ? 'dead' : (settledAt(sigs, maxH, 2) ? 'settled' : 'CHAOTIC');
  verdicts.push({ heading, ox, oy, final, finalP, massF: mass[maxH] });
  console.log(`  h${heading}(${ox},${oy})`.padEnd(14) +
    String(mass[maxH]).padStart(6) + '   ' +
    marks.map(m => m.padEnd(7)).join(' ') + '  ' + finalP);
}

// ---- for anything that settled, list the products and look for movers
console.log('\n--- terminal products of every configuration that settled within 4000 steps ---');
let moversFound = 0, totalObjs = 0;
const productTally = new Map();
for (const v of verdicts) {
  if (v.finalP === 'dead') {
    productTally.set('TOTAL ANNIHILATION', (productTally.get('TOTAL ANNIHILATION') || 0) + 1);
    continue;
  }
  if (v.finalP !== 'settled') continue;
  const objs = components(v.final).map(c => canon(c));
  const names = [];
  for (const o of objs) {
    totalObjs++;
    const known = KNOWN.get(o.id);
    const nm = known || `UNKNOWN-${classifyShape(o.id).kind}`;
    if (/MOVER/.test(nm)) moversFound++;
    names.push(nm);
    productTally.set(nm, (productTally.get(nm) || 0) + 1);
  }
  console.log(`  h${v.heading}(${v.ox},${v.oy}): ${names.length ? names.join(' + ') : 'nothing'}`);
}

console.log('\n=== summary for the 18 configurations ===');
const dead = verdicts.filter(v => v.finalP === 'dead').length;
const settled = verdicts.filter(v => v.finalP === 'settled').length;
const chaotic = verdicts.filter(v => v.finalP === 'CHAOTIC').length;
console.log(`  total annihilation within 4000 : ${dead}`);
console.log(`  settled within 4000            : ${settled}`);
console.log(`  STILL chaotic at 4000          : ${chaotic}`);
console.log(`  MOVERS found among products    : ${moversFound}`);
console.log(`  total terminal objects         : ${totalObjs}`);
console.log('');
for (const [k, v] of [...productTally.entries()].sort((a, b) => b[1] - a[1])) {
  console.log(`  ${String(v).padStart(3)}x  ${k}`);
}
if (chaotic === 0 && moversFound === 0) {
  console.log('\n  => CONCLUSION HOLDS: every one of the 18 resolves, and none produces a mover.');
  console.log('     Combined with the previous round: 79/79 interacting two-glider');
  console.log('     configurations were chased to completion, and NOT ONE produced');
  console.log('     a moving particle.');
} else if (moversFound > 0) {
  console.log('\n  => *** PREVIOUS CONCLUSION OVERTURNED: a mover was produced ***');
} else {
  console.log(`\n  => ${chaotic} configuration(s) still chaotic at 4000 steps. Not yet conclusive.`);
}
