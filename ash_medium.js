'use strict';
/*
 * ash_medium.js
 * ---------------------------------------------------------------------------
 * Is the still-life "ash" a MEDIUM that transmits, or just dead weight?
 *
 * Motivation. In rule 110 the candidate medium was an ether, and that thread is
 * now untrusted (the ether search had a bug). Life has no ether, but after a
 * cascade it freezes into a network of still lifes -- "ash". That ash is a
 * candidate medium: it is a stable background through which a disturbance might
 * propagate, and it needs no infrastructure we have not already verified.
 *
 * WHY THIS IS THE RIGHT QUESTION NOW
 *   A two-layer world needs a carrier whose timescale is slower than the
 *   substrate clock. If ash merely absorbs disturbances, there is no medium and
 *   no second layer. If ash amplifies or relays them, it IS a medium, and its
 *   relay time is the emergent clock.
 *
 * CONTROLS (all three run)
 *   vacuum  : glider through empty space -- must move (1,1) per 4 gens, unchanged
 *   sparse  : glider through a line of blocks
 *   lattice : glider through a denser grid of blocks
 *   single  : a lone block with no glider -- must simply die (no self-sustaining)
 * ---------------------------------------------------------------------------
 */

const SIZE = 400;
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
function gliderAt(x, y) {
  const s = new Set();
  for (const [dx, dy] of GLIDER0) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function blockAt(x, y) {
  const s = new Set();
  for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }

// follow one component in isolation: is it a mover, an oscillator, or dead?
function classify(seedSet) {
  const frames = [];
  let s = seedSet;
  for (let t = 0; t < 60; t++) {
    frames.push(components(s).map(c => canon(c)));
    if (s.size === 0) return { kind: 'DEBRIS (dies)' };
    s = step(s);
  }
  const sig = cs => cs.map(c => c.id).sort().join('|');
  for (let p = 1; p <= 15; p++) {
    let ok = true;
    for (let t = 0; t + p < frames.length; t++) if (sig(frames[t]) !== sig(frames[t + p])) { ok = false; break; }
    if (!ok) continue;
    const cs1 = frames[0], cs2 = frames[p];
    let moved = false;
    if (cs1.length === cs2.length) {
      for (let i = 0; i < cs1.length; i++) {
        if (Math.abs(cs2[i].cx - cs1[i].cx) > 1e-6 || Math.abs(cs2[i].cy - cs1[i].cy) > 1e-6) moved = true;
      }
    }
    return { kind: moved ? `MOVER p${p}` : (p === 1 ? 'still life' : `oscillator p${p}`), p };
  }
  return { kind: 'unsettled' };
}

function run(name, init, steps) {
  let s = init;
  const initMass = s.size;
  let peak = initMass, peakT = 0;
  const trace = [];
  for (let t = 1; t <= steps; t++) {
    s = step(s);
    if (s.size > peak) { peak = s.size; peakT = t; }
    if (t % 50 === 0) trace.push(s.size);
    if (s.size === 0) break;
  }
  const comps = components(s);
  const classes = comps.map(c => classify(new Set(c)).kind);
  const movers = classes.filter(k => /MOVER/.test(k)).length;
  console.log(`  ${name.padEnd(22)} init=${String(initMass).padStart(4)}  peak=${String(peak).padStart(4)}@t${String(peakT).padStart(3)}  final=${String(s.size).padStart(4)}  comps=${String(comps.length).padStart(3)}  movers=${movers}`);
  console.log(`      mass every 50 steps: ${trace.slice(0, 16).join(', ')}`);
  return { initMass, peak, peakT, final: s.size, comps: comps.length, movers, classes };
}

console.log('=== is the still-life ash a MEDIUM? ===');
console.log(`torus ${SIZE}x${SIZE}\n`);

const STEPS = 700;
const OX = 40, OY = 40;                 // glider start

// ---- control 1: glider in vacuum
console.log('--- control: glider in vacuum (must travel unchanged) ---');
run('glider, vacuum', gliderAt(OX, OY), STEPS);

// ---- control 2: a lone block, no glider (must die)
console.log('\n--- control: a lone block (must die; ash is not self-sustaining) ---');
run('single block', blockAt(200, 200), 60);
run('line of 5 blocks', union(blockAt(200, 200), blockAt(210, 200), blockAt(220, 200), blockAt(230, 200), blockAt(240, 200)), 60);

// ---- experiment: glider through a sparse line of blocks along its path
console.log('\n--- glider through a line of blocks on its diagonal ---');
{
  const parts = [gliderAt(OX, OY)];
  for (let i = 1; i <= 8; i++) parts.push(blockAt(OX + 20 * i, OY + 20 * i));
  run('sparse line (8 blocks)', union(...parts), STEPS);
}

// ---- experiment: glider through a dense lattice of blocks
console.log('\n--- glider through a dense lattice of blocks ---');
{
  const parts = [gliderAt(OX, OY)];
  for (let i = 1; i <= 6; i++) {
    for (let j = -2; j <= 2; j++) {
      parts.push(blockAt(OX + 20 * i + 6 * j, OY + 20 * i - 6 * j));
    }
  }
  run('dense lattice', union(...parts), STEPS);
}

// ---- experiment: glider through a random still-life field
console.log('\n--- glider through a random field of blocks ---');
{
  let seedv = 99 >>> 0;
  const rnd = () => (seedv = (seedv * 1664525 + 1013904223) >>> 0) / 4294967296;
  const parts = [gliderAt(OX, OY)];
  for (let i = 0; i < 60; i++) {
    const bx = OX + 10 + Math.floor(rnd() * 240);
    const by = OY + 10 + Math.floor(rnd() * 240);
    parts.push(blockAt(bx, by));
  }
  run('random block field', union(...parts), STEPS);
}

console.log('\n=== reading ===');
console.log('  vacuum   : glider should survive untouched -> baseline');
console.log('  ash      : if the glider is absorbed and the ash is unchanged, ash is');
console.log('             INERT (no medium, no second layer).');
console.log('             if new movers appear or the glider is deflected, the ash');
console.log('             RELAYS influence -> it is a medium.');
console.log('  lone block: must die, confirming ash is not self-sustaining on its own.');
