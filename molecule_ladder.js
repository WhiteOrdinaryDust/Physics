'use strict';
/*
 * molecule_ladder.js   --   ROUND 43 (rewritten: the previous version had a silent arithmetic bug)
 * ---------------------------------------------------------------------------
 * THE BUG THAT WAS FOUND (kept here as a permanent guard, see the first control)
 *   The grid was a Uint16Array, but the molecule's direction field sits at bits 16..17
 *   (DIRB = [1<<16, 1<<17]).  1<<17 = 131072 does not fit in 16 bits, so
 *        Uint16 = MBASE | (2 << 16)   was silently truncated to   256
 *   and dirOf() therefore ALWAYS returned 0.  The absorber did fire, but it re-emitted
 *   along direction 0 (+y) instead of the absorbed direction, and the molecule's
 *   bookkeeping momentum was always 0.  Every symptom in the last report follows from
 *   this one truncation.  Fix: 32-bit grids + a printed field-width control.
 *
 * TWO MORE FAULTS found while re-reading the phase order
 *   (a) the molecule emitted at its OLD cell and then moved to the next cell, while the
 *       emitted photon was propagated on top of it -- so molecule and photon ended the
 *       step in the SAME cell and the molecule re-absorbed its own photon.  Now the
 *       molecule moves one cell and the photon is emitted THERE, so the photon ends the
 *       step one cell ahead of the molecule.
 *   (b) the clock phase at absorption was left wherever the toggle put it, so whether the
 *       molecule decayed in the SAME step it absorbed (which makes absorption invisible --
 *       it looks like the photon passed through) or one step later was accidental.  Absorption
 *       now RESETS THE PHASE (clock := 0): the molecule never decays in the step it absorbs,
 *       and an absorbed photon is always held for exactly one step before it is re-emitted.
 *       The phase reset is also the physically sensible reading: absorption perturbs the
 *       internal state, so the molecule's own clock cannot be the thing that decides instantly.
 *
 * THE LADDER (unchanged in intent)
 *     LEVEL = how many photons the molecule holds (0..7)
 *     TICK  = its period-2 clock; on a tick it EMITS one photon and MOVES one cell
 *   so an excitation of N decays as a CASCADE of N photons, at one photon per tick, and
 *   the molecule advances exactly one cell per photon held.  Its speed is 1 cell per tick
 *   = 0.500 cells/step -- half the photon speed, which is what "it has a rest frame" means.
 *   That 0.500 is NOT put in by hand: it comes from the period-2 clock.
 *
 * MOMENTUM ACCOUNTING.  Free particle = (DX,DY); a molecule at level N = N*(DX,DY) along
 * its axis; ground molecule = 0.  I have to be honest about what this buys: under that
 * convention conservation is a BOOKKEEPING identity, not a discovery.  What it does buy is
 * that it catches algorithm faults of exactly the kind above (a wrong axis shows up as a
 * wrong total).  The PHYSICAL claims are the separate, non-tautological measurements:
 * the molecule's real displacement direction equals the absorbed photon's direction, the
 * displacement per photon is one cell, and the emitted photon travels at 1.000 along the axis.
 *
 * PREDICTIONS (written before the run)
 *   L1  controls: photon speed 1.000; ground molecule static with zero momentum; vacuum
 *       fixed; mode-0 memory still passes then reflects; field-width control shows Uint16
 *       truncating dir to 0 and Uint32 keeping it.
 *   L2  one photon absorbed at step 13; the boundary level is 1 at t=13; the tick comes at
 *       step 14, so at t=14 exactly 1 emission has happened; molecule displacement exactly
 *       (+1,0); emitted photon at speed 1.000 moving +x; total momentum (1,0) at EVERY step
 *       boundary AND at every half step (the state after the absorption phase).
 *   L3  the excited state is held for ONE step, so a second photon is absorbed only if it
 *       lands on the molecule on the very step after the first one -- a FINITE two-sided
 *       window.  Sweeping the second photon's arrival delay d (d<0 = it arrives first):
 *       |d| = 1 must reach level 2 just after the absorption phase, |d| >= 2 must not.
 *   L4  a molecule put at level 2 (hand-set control) emits exactly TWO photons, one per
 *       tick, and moves exactly 2 cells.
 *   L5  directional selectivity with everything else held fixed: a photon arriving at the
 *       SAME place at the SAME time is absorbed if it comes along the axis (+x) and passes
 *       through if it comes from another direction (+y).
 * ---------------------------------------------------------------------------
 */

const W = 160, H = 32, N = W * H;
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
const PBITS = 0xff;
const MBASE = 1 << 8;
const MODE1 = 1 << 9, MODE2 = 1 << 10;
const LVB = [1 << 11, 1 << 12, 1 << 13];        // level 0..7
const MCLK = 1 << 14;
const MREC = 1 << 15;
const DIRB = [1 << 16, 1 << 17];                // *** needs 18 bits of storage ***
const MARKER = ~PBITS;
function levOf(mol) { return ((mol >> 11) & 7); }
function dirOf(mol) { return ((mol >> 16) & 3); }
function modeOf(mol) { return ((mol >> 9) & 3); }
function withLev(mol, n) { return (mol & ~(LVB[0] | LVB[1] | LVB[2])) | ((n & 7) << 11); }
function withDir(mol, d) { return (mol & ~(DIRB[0] | DIRB[1])) | ((d & 3) << 16); }
function withMode(mol, m) { return (mol & ~(MODE1 | MODE2)) | ((m & 3) << 9); }
const TURN = new Map([[1 | 4, 2 | 8], [2 | 8, 1 | 4], [16 | 64, 32 | 128], [32 | 128, 16 | 64]]);
const BIND = new Map([[1 | 64, 0], [4 | 16, 1], [2 | 128, 2], [8 | 32, 3]]);
function reflect(parts) {
  let np = 0;
  for (let d = 0; d < 4; d++) {
    if ((parts >> d) & 1) np |= 1 << ((d + 2) % 4);
    if ((parts >> (d + 4)) & 1) np |= 1 << (((d + 2) % 4) + 4);
  }
  return np;
}
function singleDir(parts) { for (let d = 0; d < 4; d++) if ((parts >> d) & 1) return d; return -1; }

// ---------------------------------------------------------------- one local rule, three phases
function collide(g) {
  const o = new Uint32Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const parts = v & PBITS, mol = v & MARKER;
    if (mol) {
      const m = mol ^ MCLK;                                   // the clock always ticks
      if (modeOf(mol) === 1 && parts) {                       // MODE 1 = absorber
        const d = singleDir(parts);
        const lv = levOf(mol);
        const reset = m & ~MCLK;                              // absorption RESETS the internal phase
        if (lv === 0) { o[i] = withLev(withDir(reset & ~MREC, d < 0 ? 0 : d), 1); continue; }
        if (d === dirOf(mol) && lv < 7) { o[i] = withLev(reset, lv + 1); continue; }
        o[i] = m | parts; continue;                           // wrong axis: transparent
      }
      if (!parts) { o[i] = m; continue; }
      if (m & MREC) o[i] = (m & ~MREC) | reflect(parts);       // MODE 0 memory: READ
      else o[i] = (m | MREC) | parts;                          // MODE 0 memory: WRITE
      continue;
    }
    if (TURN.has(parts)) { o[i] = TURN.get(parts); continue; }
    if (BIND.has(parts)) { o[i] = MBASE | (BIND.get(parts) << 16); continue; }
    o[i] = parts;
  }
  return o;
}
// on a tick: the molecule moves ONE cell along its axis and emits ONE photon there
function moleculePhase(g) {
  const o = new Uint32Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i], mol = v & MARKER, parts = v & PBITS;
    if (!mol) { o[i] |= v; continue; }
    const ticked = (mol & MCLK) !== 0;
    const lv = levOf(mol);
    if (!ticked || lv === 0) { o[i] |= mol | parts; continue; }
    const d = dirOf(mol);
    const x = i % W, y = (i - x) / W;
    const k = ((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W);
    o[i] |= parts;                                            // any free particle stays behind
    o[k] |= withLev(mol, lv - 1) | (1 << d);                  // move one cell, emit one photon
  }
  return o;
}
function propagate(g) {
  const o = new Uint32Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i], keep = v & MARKER;
    if (keep) o[i] |= keep;
    if (!(v & PBITS)) continue;
    const x = i % W, y = (i - x) / W;
    for (let d = 0; d < 4; d++) {
      if ((v >> d) & 1) o[((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W)] |= 1 << d;
      if ((v >> (d + 4)) & 1) o[((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W)] |= 1 << (d + 4);
    }
  }
  return o;
}
function step(g) { return propagate(moleculePhase(collide(g))); }

// ---------------------------------------------------------------- measurement
function freeInfo(g) {
  let n = 0, px = 0, py = 0, ax = -1, ay = -1;
  for (let i = 0; i < N; i++) {
    const v = g[i] & PBITS;
    if (!v) continue;
    const x = i % W, y = (i - x) / W;
    for (let d = 0; d < 4; d++) {
      if ((v >> d) & 1) { n++; px += DX[d]; py += DY[d]; if (ax < 0) { ax = x; ay = y; } }
      if ((v >> (d + 4)) & 1) { n++; px += DX[d]; py += DY[d]; if (ax < 0) { ax = x; ay = y; } }
    }
  }
  return { n, px, py, ax, ay };
}
function molInfo(g) {
  for (let i = 0; i < N; i++) {
    const mol = g[i] & MARKER;
    if (!mol) continue;
    const lv = levOf(mol), d = dirOf(mol);
    return { cell: i, x: i % W, y: (i - i % W) / W, lv, dir: d, mode: modeOf(mol), px: lv * DX[d], py: lv * DY[d] };
  }
  return null;
}
function total(g) {
  const f = freeInfo(g), m = molInfo(g);
  return { px: f.px + (m ? m.px : 0), py: f.py + (m ? m.py : 0), f, m };
}
function newGrid() { return new Uint32Array(N); }
function withMolecule(mx, mode) {
  const g = newGrid();
  g[16 * W + (mx - 2)] = 1 << 1;              // colour A moving +x
  g[16 * W + (mx + 2)] = 1 << 7;              // colour B moving -x  -> they bind
  let a = Uint32Array.from(g);
  for (let t = 0; t < 3; t++) a = step(a);
  a[16 * W + mx] = withMode(a[16 * W + mx], mode);
  return a;
}
function clearPhotons(g) { for (let i = 0; i < N; i++) g[i] &= MARKER; return g; }

// ============================================================ controls
console.log('============================================================');
console.log('CONTROLS');
console.log('============================================================');
{
  const a = new Uint16Array(1), b = new Uint32Array(1);
  a[0] = MBASE | (2 << 16); b[0] = MBASE | (2 << 16);
  console.log(`  field-width control: Uint16 stores ${a[0]} -> dirOf ${(a[0] >> 16) & 3} (WRONG) | Uint32 stores ${b[0]} -> dirOf ${(b[0] >> 16) & 3} (right)`);
  console.log(`    => the grid must be 32-bit; this control exists so the bug cannot come back silently.`);
}
{
  const g = newGrid(); g[16 * W + 10] = 1 << 1;
  let a = g; for (let t = 0; t < 20; t++) a = step(a);
  const f = freeInfo(a);
  console.log(`  L1 photon alone: reached x=${f.ax} from x=10 in 20 steps => speed ${((f.ax - 10) / 20).toFixed(3)}  ${Math.abs((f.ax - 10) / 20 - 1) < 1e-9 ? 'PASS' : '*** FAIL ***'}`);
}
{
  let a = withMolecule(22, 0);
  const m0 = molInfo(a);
  for (let t = 0; t < 40; t++) a = step(a);
  const m1 = molInfo(a), tot = total(a);
  console.log(`  L1 ground molecule: x ${m0.x} -> ${m1.x}, level ${m1.lv}, momentum ${tot.px},${tot.py}  ${m0.x === m1.x && tot.px === 0 && tot.py === 0 ? 'PASS (static, no momentum)' : '*** FAIL ***'}`);
}
{
  let e = newGrid(); for (let t = 0; t < 30; t++) e = step(e);
  console.log(`  L1 empty vacuum fixed: ${freeInfo(e).n === 0 ? 'PASS' : '*** FAIL ***'}`);
}
{
  let g = withMolecule(22, 0);
  const shot = (board) => {
    board[16 * W + 10] = 1 << 1;
    let dirAfter = null, meet = -1;
    for (let k = 0; k < 40; k++) {
      board = step(board);
      let px = -1;
      for (let x = 0; x < W; x++) if (board[16 * W + x] & (1 << 1)) px = x;
      if (meet < 0 && px >= 0 && Math.abs(px - 22) <= 1 && k > 2) meet = k;
      if (meet >= 0 && k === meet + 6) { dirAfter = px - 22; break; }
    }
    clearPhotons(board);
    return { dirAfter, board };
  };
  const r1 = shot(g); g = r1.board;
  const r2 = shot(g); g = r2.board;
  console.log(`  L1 mode-0 memory still: shot1 ${r1.dirAfter > 0 ? 'PASSED' : 'reflected'}, shot2 ${r2.dirAfter > 0 ? 'passed' : 'REFLECTED'}  ${r1.dirAfter > 0 && r2.dirAfter <= 0 ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ============================================================ L2: one photon
console.log('============================================================');
console.log('L2: one photon in -- absorbed -- level 1 -- one cell of recoil -- one photon out');
console.log('============================================================');
{
  let g = withMolecule(22, 1);
  const x0 = molInfo(g).x;
  g[16 * W + 10] = 1 << 1;
  console.log('     t   分子x  lv(boundary)  lv(after absorb)  自由光子  free p   分子 p   TOTAL p');
  let bad = 0, badHalf = 0, maxBound = 0, maxHalf = 0, emitted = 0, prevN = 1, absorbT = -1;
  const trail = [];
  for (let t = 0; t <= 30; t++) {
    const T = total(g);
    const cT = total(collide(g));                 // half-step: state after the absorption phase
    if (!(T.px === 1 && T.py === 0)) bad++;
    if (!(cT.px === 1 && cT.py === 0)) badHalf++;
    if (T.m) { maxBound = Math.max(maxBound, T.m.lv); trail.push(T.m.x); }
    if (cT.m) maxHalf = Math.max(maxHalf, cT.m.lv);
    if (T.f.n === 0 && absorbT < 0) absorbT = t;
    if (T.f.n > prevN) emitted += T.f.n - prevN;
    prevN = T.f.n;
    if (t <= 18) console.log(`    ${String(t).padStart(3)}   ${T.m ? String(T.m.x).padStart(4) : '  - '}       ${T.m ? T.m.lv : '-'}              ${cT.m ? cT.m.lv : '-'}            ${String(T.f.n).padStart(3)}      ${String(T.f.px).padStart(3)}    ${String(T.m ? T.m.px : 0).padStart(4)}     ${T.px},${T.py}`);
    g = step(g);
  }
  const T = total(g), x1 = T.m.x;
  const p20 = (() => { let a = withMolecule(22, 1); a[16 * W + 10] = 1 << 1; for (let t = 0; t < 20; t++) a = step(a); return freeInfo(a); })();
  const p30 = (() => { let a = withMolecule(22, 1); a[16 * W + 10] = 1 << 1; for (let t = 0; t < 30; t++) a = step(a); return freeInfo(a); })();
  console.log('');
  console.log(`    level reached: at a step boundary ${maxBound}, just after the absorption phase ${maxHalf}   ${maxHalf === 1 && maxBound === 1 ? 'PASS (absorbed, level 1, and the level is VISIBLE)' : '*** FAIL ***'}`);
  console.log(`    absorbed at step ${absorbT}   ${absorbT === 13 ? 'PASS (as predicted)' : '*** FAIL (predicted 13) ***'}`);
  console.log(`    emissions counted: ${emitted}   ${emitted === 1 ? 'PASS (exactly one photon out)' : '*** FAIL ***'}`);
  console.log(`    molecule displacement: x ${x0} -> ${x1}  = ${x1 - x0} cell(s) to the +x   ${x1 - x0 === 1 ? 'PASS (recoil along the absorbed direction, 1 cell per photon = 0.500 cells/step)' : '*** FAIL ***'}`);
  console.log(`    emitted photon: at t=20 x=${p20.ax}, at t=30 x=${p30.ax} => speed ${((p30.ax - p20.ax) / 10).toFixed(3)}, direction p=(${p30.px},${p30.py})   ${Math.abs((p30.ax - p20.ax) / 10 - 1) < 1e-9 && p30.px === 1 && p30.py === 0 ? 'PASS' : '*** FAIL ***'}`);
  console.log(`    steps whose total momentum was not (1,0): boundary ${bad}, half-step ${badHalf}   ${bad === 0 && badHalf === 0 ? 'PASS: conserved at EVERY step' : '*** FAIL ***'}`);
  console.log(`    (honest note: under the bookkeeping convention conservation is an identity; what is NOT an identity here is the`);
  console.log(`     measured recoil direction and the emitted photon's speed.)`);
}
console.log('');

// ============================================================ L3: two photons, finite window
console.log('============================================================');
console.log('L3: two photons -- how long can the molecule stay excited?  sweep the launch offset');
console.log('============================================================');
{
  console.log('  photon 1 starts at x=10.  photon 2 starts at x=10-d, so it arrives d steps LATER');
  console.log('  (d>0 = behind photon 1, d<0 = ahead of it, both are sent along +x).');
  console.log('     d   max lv (boundary)  max lv (after absorb)  step it hit level 2  emissions  final free p');
  for (let d = -3; d <= 4; d++) {
    if (d === 0) continue;                       // same cell = the same photon, nothing to test
    let g = withMolecule(22, 1);
    g[16 * W + 10] = 1 << 1;
    g[16 * W + (10 - d)] = 1 << 1;
    let maxBound = 0, maxHalf = 0, emitted = 0, prevN = 2, half2At = -1;
    for (let t = 0; t <= 90; t++) {
      const T = total(g), cT = total(collide(g));
      if (T.m) maxBound = Math.max(maxBound, T.m.lv);
      if (cT.m) { maxHalf = Math.max(maxHalf, cT.m.lv); if (cT.m.lv === 2 && half2At < 0) half2At = t; }
      if (T.f.n > prevN) emitted += T.f.n - prevN;
      prevN = T.f.n;
      g = step(g);
    }
    const T = total(g);
    const want = Math.abs(d) === 1 ? 2 : 1;
    console.log(`    ${String(d).padStart(2)}          ${maxBound}                   ${maxHalf}                ${String(half2At).padStart(4)}             ${String(emitted).padStart(3)}       ${T.f.px},${T.f.py}   ${maxHalf === want ? 'as predicted' : '*** NOT AS PREDICTED (wanted ' + want + ') ***'}`);
  }
  console.log('  PREDICTION: the excited state is held for ONE step, so only |d| = 1 puts the second');
  console.log('  photon on the molecule while it is still excited => level 2 just after absorption.');
  console.log('  |d| >= 2 arrives after the molecule has already decayed => the level stays at 1.');
}
console.log('');

// ============================================================ L4: cascade from level 2
console.log('============================================================');
console.log('L4: a molecule at level 2 (hand-set control) -- a cascade of exactly TWO photons');
console.log('============================================================');
{
  let g = withMolecule(22, 1);
  const c = 16 * W + 22;
  g[c] = withLev(withDir(g[c], 1), 2);            // hand-set control: level 2, axis +x
  const x0 = molInfo(g).x;
  console.log('     t   分子x  lv   自由光子  free p   分子 p   TOTAL p');
  let emitted = 0, prevN = 0, maxBound = 0;
  for (let t = 0; t <= 14; t++) {
    const T = total(g);
    if (T.m) maxBound = Math.max(maxBound, T.m.lv);
    if (T.f.n > prevN) emitted += T.f.n - prevN;
    prevN = T.f.n;
    if (t <= 12) console.log(`    ${String(t).padStart(3)}   ${T.m ? String(T.m.x).padStart(4) : '  - '}    ${T.m ? T.m.lv : '-'}      ${String(T.f.n).padStart(3)}      ${String(T.f.px).padStart(3)}    ${String(T.m ? T.m.px : 0).padStart(4)}     ${T.px},${T.py}`);
    g = step(g);
  }
  const T = total(g);
  console.log(`    max level seen ${maxBound}, photons emitted ${emitted}, molecule moved ${T.m.x - x0} cells, final free momentum (${T.f.px},${T.f.py})`);
  console.log(`    => ${maxBound === 2 && emitted === 2 && T.m.x - x0 === 2 && T.f.px === 2 ? 'PASS: 2 photons out, 2 cells of recoil, one emission per tick' : '*** FAIL ***'}`);
}
console.log('');
// is the cascade a general ladder, or was level 2 a special case?
{
  console.log('L4b: the cascade for levels 1..4 (all hand-set controls).  Each level N must give');
  console.log('     exactly N photons out, N cells of recoil, and momentum (N,0) at EVERY step.');
  console.log('     level   photons out   cells moved   final free p   momentum (N,0) every step?   gaps between emissions');
  for (let L = 1; L <= 4; L++) {
    let g = withMolecule(22, 1);
    const c = 16 * W + 22;
    g[c] = withLev(withDir(g[c], 1), L);
    const x0 = molInfo(g).x;
    let emitted = 0, prevN = 0, ok = true;
    const emitAt = [];
    for (let t = 0; t <= 2 * L + 6; t++) {
      const T = total(g);
      if (!(T.px === L && T.py === 0)) ok = false;
      if (T.f.n > prevN) { emitted += T.f.n - prevN; for (let k = prevN; k < T.f.n; k++) emitAt.push(t); }
      prevN = T.f.n;
      g = step(g);
    }
    const T = total(g);
    const gaps = emitAt.slice(1).map((v, i) => v - emitAt[i]).join(',');
    console.log(`      ${L}          ${emitted}            ${T.m.x - x0}            ${T.f.px},${T.f.py}          ${ok ? 'yes' : '*** NO ***'}                     ${gaps || '-'}`);
  }
}
console.log('');

// ============================================================ L5: directional selectivity
console.log('============================================================');
console.log('L5: same place, same time, different direction');
console.log('============================================================');
{
  // photon 1 lands on the molecule at the end of step 12, is absorbed at step 13, is held
  // through the t=13 boundary, and the molecule ticks at step 14.  The test photon is timed
  // to be in the same cell DURING step 14's absorption phase, so the molecule is excited when
  // it arrives.  Everything is identical in the two runs except the test photon's direction.
  const run = (dir) => {
    let g = withMolecule(22, 1);
    g[16 * W + 10] = 1 << 1;                                    // the exciting photon, +x
    if (dir === 1) g[16 * W + 9] = 1 << 1;                      // test photon +x, in the cell at step 14
    else g[(16 - 13 + H) % H * W + 22] = 1 << 0;                // test photon +y from below, same arrival
    let half14 = 0, bound14 = null, sawYabove = false;
    for (let t = 0; t <= 20; t++) {
      if (t === 13) half14 = total(collide(g)).m.lv;            // collide at the t=13 boundary IS step 14's absorb phase
      if (t === 14) bound14 = total(g);
      g = step(g);
      for (let y = 0; y < H; y++) if ((g[y * W + 22] & (1 << 0)) && y < 16) sawYabove = true;
    }
    return { half14, bound14, sawYabove };
  };
  const ax = run(1), ay = run(0);
  console.log(`  test photon along the axis (+x): level just after step 14's absorb phase = ${ax.half14}; at t=14 level ${ax.bound14.m.lv}, free photons ${ax.bound14.f.n}, free p (${ax.bound14.f.px},${ax.bound14.f.py})   ${ax.half14 === 2 ? 'ABSORBED' : '*** not absorbed ***'}`);
  console.log(`  test photon from +y            : level just after step 14's absorb phase = ${ay.half14}; at t=14 level ${ay.bound14.m.lv}, free photons ${ay.bound14.f.n}, free p (${ay.bound14.f.px},${ay.bound14.f.py}); +y photon seen above the molecule: ${ay.sawYabove}   ${ay.half14 === 1 ? 'NOT absorbed' : '*** absorbed ***'}`);
  console.log(`  => ${ax.half14 === 2 && ay.half14 === 1 && ay.sawYabove ? 'PASS: same place, same time, different direction -> one absorbed, one not' : '*** FAIL ***'}`);
}
