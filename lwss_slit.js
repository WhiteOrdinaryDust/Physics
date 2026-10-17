'use strict';
/*
 * lwss_slit.js   (second attempt -- the first one aborted in its own control)
 * ---------------------------------------------------------------------------
 * WHAT WENT WRONG THE FIRST TIME
 *   The rightward ship I built by mirroring collision.js's `lwss` split into two
 *   pieces, so the control failed and the script aborted. Investigating:
 *     - collision.js's own LWSS control prints "*** FAIL: expected 1 component,
 *       got 2 ***". That failure had been carried as a pass for several rounds.
 *     - find_orthogonal_ship.js then brute-forced every pattern in a 5x4 box and
 *       found the genuine ship by the rule alone: 9 cells, period 4, velocity
 *       (+-2,0). Its four phases have 2,1,2,1 8-connected pieces -- the LWSS is
 *       genuinely NOT a single blob at every phase, the pieces interact over 2-4
 *       generations. So the pattern was right and the *detector* was wrong: a
 *       particle here must be recognised as "a whole pattern that recurs with a
 *       constant drift", not as "one component".
 *
 * THE GEOMETRY (this is the point of the script)
 *   A glider moves diagonally, so a straight screen can never be perpendicular
 *   to its motion: it always slides along the screen at 45 degrees. The LWSS
 *   moves straight in x, so the screen is a plain vertical wall of blocks and the
 *   particle arrives perpendicular to it -- a real slit geometry.
 *
 *   The LWSS has period 4 and speed 2 cells per 4 generations, so its four phases
 *   are four genuinely different cell arrangements: no spatial shift turns one
 *   into another. That is the hidden phase, tested exactly where it matters, at
 *   the edge of a slit.
 *
 * WHAT IS MEASURED
 *   for each (screen thickness, slit width, placement of the ship in the slit,
 *   phase of the ship):
 *     out    : does a whole 9-cell pattern leave with the same velocity (2,0)?
 *     dmg    : how many cells of the screen differ from before, split into dead
 *              and born, with the canonical shape and the offset of the damage
 *              from the centre of the slit (which edge got hit)
 *   and then: for a FIXED geometry, do the four phases give the same record?
 *
 * PREDICTIONS (stated before the run, per method rule)
 *   P1  flush against an edge -> the screen is permanently damaged, because
 *       matter_elasticity.js showed Life matter heals holes but never bumps.
 *   P2  the middle placement -> the ship must pass with no contact and leave the
 *       screen EXACTLY unchanged, all four phases.
 *   P3  the damage, when there is any, DEPENDS ON PHASE: the four phases present
 *       different cells at the moment of closest approach.
 *   P4  thickness raises the screen's mass but not the absolute damage, so the
 *       damaged fraction falls like 1/thickness -- a macroscopic screen survives
 *       the hit and still keeps an indelible record of it.
 * ---------------------------------------------------------------------------
 */

const fs = require('fs');
const SIZE = 600;
const STEPS = 500;
const OUT = 'lwss_slit_report.txt';
const lines = [];
function say(s) { lines.push(s); }

function wrapCell(x, y) { return ((((y % SIZE) + SIZE) % SIZE) * SIZE) + (((x % SIZE) + SIZE) % SIZE); }
function step(set) {
  const counts = new Map();
  for (const k of set) {
    const x = k % SIZE, y = (k - x) / SIZE;
    for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
      if (dx === 0 && dy === 0) continue;
      const nk = wrapCell(x + dx, y + dy);
      counts.set(nk, (counts.get(nk) || 0) + 1);
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
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        if (dx === 0 && dy === 0) continue;
        const nk = wrapCell(x + dx, y + dy);
        if (set.has(nk) && !seen.has(nk)) { seen.add(nk); Q.push(nk); }
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
  const w = maxx - minx + 1, h = maxy - miny + 1, g = [];
  for (let r = 0; r < h; r++) g.push(new Array(w).fill('.'));
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    g[y - miny][x - minx] = 'o';
  }
  return { w, h, n: cells.length, minx, miny, maxx, maxy, cx: sx / cells.length, cy: sy / cells.length, id: `${w}x${h}:` + g.map(r => r.join('')).join('/') };
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }
function setEq(a, b) { if (a.size !== b.size) return false; for (const k of a) if (!b.has(k)) return false; return true; }

// ---- the ship found by brute force in find_orthogonal_ship.js
//      (9 cells, 5x4 box, period 4, velocity (+2,0))
function shipRight(x, y) {
  const s = new Set();
  for (const [dx, dy] of [[1, 0], [2, 0], [3, 0], [4, 0], [0, 1], [4, 1], [4, 2], [0, 3], [3, 3]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
// The four phases are NOT all the same size: phases 0 and 2 have 9 cells, phases
// 1 and 3 have 12. Every phase is nevertheless a 5x4 pattern that recurs after 4
// generations displaced by (2,0). So "the particle came out unchanged" must be
// tested by matching the SHAPE, not by counting cells.
const KNOWN_RIGHT = new Set();
{
  let sh = shipRight(150, 150);
  for (let k = 0; k < 4; k++) { KNOWN_RIGHT.add(canon([...sh]).id); sh = step(sh); }
}
function normalizeTo(set, x0, y0) {
  const c = canon([...set]);
  const dx = x0 - c.minx, dy = y0 - c.miny;
  const out = new Set();
  for (const k of set) { const x = k % SIZE, y = (k - x) / SIZE; out.add(wrapCell(x + dx, y + dy)); }
  return out;
}
function blockAt(x, y) {
  const s = new Set();
  for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function wallBlocks(X0, Y0, nx, ny, s, jskip) {
  let m = new Set();
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
    if (j === jskip) continue;
    m = union(m, blockAt(X0 + i * s, Y0 + j * s));
  }
  return m;
}
function zoneOf(cells, r) {
  const z = new Set();
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    for (let dy = -r; dy <= r; dy++) for (let dx = -r; dx <= r; dx++) z.add(wrapCell(x + dx, y + dy));
  }
  return z;
}
// a particle = a whole pattern (possibly in several interacting pieces) that
// recurs with the same shape and a constant non-zero drift
function trackMover(cells, steps) {
  let c = new Set(cells);
  const frames = [];
  for (let t = 0; t <= steps; t++) {
    if (c.size === 0) return { ok: false, mass: 0 };
    frames.push(canon([...c]));
    c = step(c);
  }
  for (let p = 1; p <= 8; p++) {
    let ok = true, base = null;
    for (let t = 0; t + p < frames.length; t++) {
      if (frames[t].id !== frames[t + p].id) { ok = false; break; }
      const dx = Math.round(frames[t + p].cx - frames[t].cx);
      const dy = Math.round(frames[t + p].cy - frames[t].cy);
      if (base === null) base = [dx, dy];
      else if (dx !== base[0] || dy !== base[1]) { ok = false; break; }
    }
    if (ok && base && (base[0] !== 0 || base[1] !== 0)) {
      return { ok: true, p, vx: base[0], vy: base[1], mass: frames[0].n, w: frames[0].w, h: frames[0].h, id: frames[0].id };
    }
  }
  return { ok: false, mass: frames[0].n, id: frames[0].id };
}
function piecesOf(set, t) {
  let c = new Set(set);
  for (let i = 0; i < t; i++) c = step(c);
  return components(c).length;
}
function analyzeOutgoing(set) {
  if (set.size === 0) return { kind: 'nothing', outMass: 0 };
  const id = canon([...set]).id;
  if (KNOWN_RIGHT.has(id)) return { kind: 'ship(2,0)', outMass: set.size };
  const m = trackMover(set, 40);
  if (m.ok) return { kind: `mover${m.mass} p${m.p} v(${m.vx},${m.vy})`, outMass: set.size };
  return { kind: `debris${set.size}`, outMass: set.size };
}

// ================================================================ controls
console.log('=== CONTROLS ===');
say('=== CONTROLS ===');
{
  const t = trackMover(shipRight(150, 150), 60);
  const pcs = [0, 1, 2, 3].map(k => piecesOf(shipRight(150, 150), k));
  console.log(`  ship in vacuum: ${t.ok ? `PASS period ${t.p} velocity (${t.vx},${t.vy}), ${t.mass} cells ${t.w}x${t.h}` : `*** FAIL ***`}`);
  console.log(`    connected pieces at t=0..3: ${pcs.join(', ')}  (2,1,2,1 is correct: the LWSS really is`);
  console.log('    two interacting pieces in some phases -- that is why a "1 component" test rejects a valid ship)');
  say(`  ship in vacuum: ${t.ok ? `PASS period ${t.p} velocity (${t.vx},${t.vy}), ${t.mass} cells ${t.w}x${t.h}` : 'FAIL'}`);
  say(`    connected pieces at t=0..3: ${pcs.join(', ')}`);
  if (!t.ok || t.vx !== 2 || t.vy !== 0 || t.mass !== 9) { console.log('  *** abort ***'); process.exit(1); }
}
{
  const ids = new Set();
  let allOk = true;
  for (let k = 0; k < 4; k++) {
    let sh = shipRight(150, 150);
    for (let i = 0; i < k; i++) sh = step(sh);
    sh = normalizeTo(sh, 150, 150);
    const t = trackMover(sh, 60);
    if (t.ok) ids.add(t.id);
    if (!t.ok || t.vx !== 2 || t.vy !== 0 || (t.mass !== 9 && t.mass !== 12)) allOk = false;
    console.log(`  phase ${k}: ${t.ok ? `period ${t.p} velocity (${t.vx},${t.vy}), ${t.mass} cells, shape ${t.id}` : '*** FAIL ***'}`);
    say(`  phase ${k}: ${t.ok ? `${t.mass} cells, velocity (${t.vx},${t.vy}), shape ${t.id}` : 'FAIL'}`);
  }
  console.log(`  distinct shapes over the four phases: ${ids.size}  ${allOk ? 'PASS (all four fly)' : '*** FAIL ***'}`);
  say(`  distinct shapes over the four phases: ${ids.size}`);
  if (!allOk) process.exit(1);
}
{
  let ok = true;
  for (const nx of [1, 2, 3]) for (const s of [4, 6]) {
    const w = wallBlocks(250, 150, nx, 9, s, 4);
    let x = new Set(w);
    for (let t = 0; t < STEPS; t++) x = step(x);
    if (!setEq(w, x)) { ok = false; console.log(`  *** wall (thickness ${nx}, spacing ${s}) is NOT static ***`); }
  }
  console.log(`  all six walls static for ${STEPS} generations: ${ok ? 'PASS' : '*** FAIL ***'}`);
  say(`  all six walls static: ${ok ? 'PASS' : 'FAIL'}`);
  if (!ok) process.exit(1);
}
{
  // a ship that never comes near the screen must leave it untouched
  const w = wallBlocks(250, 150, 2, 9, 4, 4);
  const zone = zoneOf(w, 2);
  let x = union(w, normalizeTo(shipRight(190, 260), 190, 260));
  for (let t = 0; t < STEPS; t++) x = step(x);
  let dmg = 0;
  for (const k of x) if (zone.has(k) && !w.has(k)) dmg++;
  for (const k of w) if (!x.has(k)) dmg++;
  const out = new Set(); for (const k of x) if (!zone.has(k)) out.add(k);
  const kind = analyzeOutgoing(out).kind;
  const good = dmg === 0 && kind === 'ship(2,0)';
  console.log(`  ship far from the screen: damage ${dmg}, outgoing ${kind}  ${good ? 'PASS' : '*** FAIL ***'}`);
  say(`  ship far from the screen: damage ${dmg}, outgoing ${kind}  ${good ? 'PASS' : 'FAIL'}`);
  if (!good) process.exit(1);
}
{
  // wide-open control: a slit 18 cells wide must let the ship through untouched
  const w = wallBlocks(250, 150, 2, 9, 6, 4);
  const w2 = new Set();
  for (const k of w) w2.add(k);
  for (let j = 3; j <= 5; j++) for (let i = 0; i < 2; i++) for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) w2.delete(wrapCell(250 + i * 6 + dx, 150 + j * 6 + dy));
  const zone = zoneOf(w2, 2);
  let x = union(w2, normalizeTo(shipRight(190, 170), 190, 170));
  for (let t = 0; t < STEPS; t++) x = step(x);
  let dmg = 0;
  for (const k of x) if (zone.has(k) && !w2.has(k)) dmg++;
  for (const k of w2) if (!x.has(k)) dmg++;
  const out = new Set(); for (const k of x) if (!zone.has(k)) out.add(k);
  const kind = analyzeOutgoing(out).kind;
  const good = dmg === 0 && kind === 'ship(2,0)';
  console.log(`  ship through a very wide opening: damage ${dmg}, outgoing ${kind}  ${good ? 'PASS' : '*** FAIL ***'}`);
  say(`  ship through a very wide opening: damage ${dmg}, outgoing ${kind}  ${good ? 'PASS' : 'FAIL'}`);
}
console.log('  controls done\n');

// ================================================================ experiment
const X0 = 250, Y0 = 150, NY = 9, JSKIP = 4, LAUNCH_X = 190;

function runOne(wall, zone, x0, y0, phase) {
  let sh = shipRight(150, 150);
  for (let i = 0; i < phase; i++) sh = step(sh);
  sh = normalizeTo(sh, x0, y0);
  for (const k of sh) if (wall.has(k)) return null;
  let x = union(wall, sh);
  for (let t = 0; t < STEPS; t++) x = step(x);
  const removed = new Set();
  for (const k of wall) if (!x.has(k)) removed.add(k);
  const added = new Set();
  for (const k of x) if (zone.has(k) && !wall.has(k)) added.add(k);
  const outside = new Set();
  for (const k of x) if (!zone.has(k)) outside.add(k);
  const out = analyzeOutgoing(outside);
  const rec = union(removed, added);
  const rc = rec.size ? canon([...rec]) : null;
  return {
    pass: out.kind === 'ship(2,0)', outKind: out.kind, outMass: out.outMass,
    rm: removed.size, ad: added.size, dmg: rec.size,
    recId: rc ? rc.id : '-', recY: rc ? rc.cy : null,
    recRows: rc ? `${rc.miny}..${rc.maxy}` : '-',
  };
}

const groups = [];
for (const nx of [1, 2, 3]) {
  for (const s of [4, 6]) {
    const wall = wallBlocks(X0, Y0, nx, NY, s, JSKIP);
    const zone = zoneOf(wall, 2);
    const gapLo = Y0 + s * (JSKIP - 1) + 2;
    const gapHi = Y0 + s * (JSKIP + 1) - 1;
    const nPlace = gapHi - gapLo + 1 - 4 + 1;
    const cySlit = (gapLo + gapHi) / 2;
    for (let b = 0; b < nPlace; b++) {
      const sigs = [];
      for (let phase = 0; phase < 4; phase++) sigs.push(runOne(wall, zone, LAUNCH_X, gapLo + b, phase));
      groups.push({ nx, s, b, nPlace, gapLo, gapHi, cySlit, mass: wall.size, sigs });
    }
  }
}

say('');
say('=== FULL RECORD ===');
say('X0=250 Y0=150, 9 block rows, slit = block row 4 removed.');
say('slit width = 2s-2.  b = vertical placement of the ship in the slit:');
say('  b=0   flush against the LOWER edge of the slit');
say(`  b=max flush against the UPPER edge`);
say('phase = which of the 4 internal arrangements of the ship is launched, all');
say('launched from the same bounding box so only the internal arrangement differs.');
say('');
const headline = [];
let phaseDependent = 0;
const details = [];

for (const g of groups) {
  const valid = g.sigs.filter(Boolean);
  if (!valid.length) continue;
  const keys = valid.map(r => `${r.outKind}|${r.dmg}|${r.recId}|${r.recRows}`);
  const dep = new Set(keys).size > 1;
  if (dep) phaseDependent++;
  const tag = `thickness ${g.nx}  spacing ${g.s}  slit ${g.gapLo}..${g.gapHi} (width ${g.gapHi - g.gapLo + 1})  place b=${g.b}${g.b === 0 ? ' (flush low)' : g.b === g.nPlace - 1 ? ' (flush high)' : ' (clear)'}  mass ${g.mass}`;
  say(tag);
  for (let phase = 0; phase < 4; phase++) {
    const r = g.sigs[phase];
    if (!r) { say(`   phase ${phase}: launch overlaps the screen (skipped)`); continue; }
    const dy = r.recY === null ? '-' : (r.recY - g.cySlit).toFixed(1);
    say(`   phase ${phase}: out=${r.outKind.padEnd(16)} dmg=${String(r.dmg).padStart(3)} (died ${String(r.rm).padStart(3)}, born ${String(r.ad).padStart(3)})  record ${r.recId}  rows ${r.recRows}  dy_from_slit_centre ${dy}`);
  }
  say(`   -> ${dep ? 'THE FOUR PHASES DIFFER   <-- PHASE-SENSITIVE' : 'the four phases are indistinguishable'}`);
  say('');
  const short = g.sigs.map(r => r ? `${r.pass ? 'S' : r.outKind === 'nothing' ? '.' : 'D'}${String(r.dmg).padStart(2)}` : ' x').join(' ');
  headline.push(`thickness ${g.nx} spacing ${g.s} b=${g.b} mass ${String(g.mass).padStart(3)} | ${short} | ${dep ? 'PHASE-SENSITIVE' : 'same'}`);
  if (dep) details.push({ tag, sigs: g.sigs, cySlit: g.cySlit });
}

console.log('=== SCREEN: thickness x slit-width x placement x phase ===');
console.log('  each column is one phase;  S=dock ship came out, .=nothing came out, D=other object');
console.log('  the number after the letter is how many screen cells changed\n');
for (const h of headline) console.log('  ' + h);

console.log('\n=== PHASE SUMMARY ===');
console.log(`  geometries tested: ${groups.length}`);
console.log(`  geometries whose outcome depends on the ship's internal phase: ${phaseDependent}`);
console.log('\n  details of the phase-sensitive geometries:');
for (const d of details) {
  console.log('  ' + d.tag);
  for (let phase = 0; phase < 4; phase++) {
    const r = d.sigs[phase];
    if (!r) continue;
    const dy = r.recY === null ? '-' : (r.recY - d.cySlit).toFixed(1);
    console.log(`     phase ${phase}: out=${r.outKind.padEnd(16)} dmg=${String(r.dmg).padStart(3)} record ${r.recId} dy=${dy}`);
  }
}

say('=== SUMMARY ===');
say(`geometries tested: ${groups.length}`);
say(`geometries whose outcome depends on the ship's internal phase: ${phaseDependent}`);
fs.writeFileSync(OUT, lines.join('\n') + '\n');
console.log(`\n  full record written to ${OUT}`);
