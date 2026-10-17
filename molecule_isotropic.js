'use strict';
/*
 * molecule_isotropic.js   --   ROUND 44
 * ---------------------------------------------------------------------------
 * WHAT WAS STILL TIED TO THE MOLECULE'S SINGLE INTERNAL AXIS (round 43):
 *   (1) ABSORPTION was not isotropic: a photon was taken only if the molecule was in the
 *       ground state, or if its direction equalled the stored axis.  A photon from any other
 *       direction was simply transparent.  Consequence: only collinear photons could pile up,
 *       and level >= 3 was unreachable without hand-setting it.
 *   (2) EMISSION always left along that same stored axis, so the emitted direction was a
 *       copy of the absorbed direction -- there was no room for a hidden variable to matter.
 *
 * THE CHANGE.  The molecule's internal state is now three separate things:
 *     LEVEL   = how many photons it is holding                     (3 bits, 0..7)   energy
 *     (px,py) = net momentum it has absorbed minus emitted          (2 x 4-bit signed) momentum
 *     PHASE   = a free-running 2-bit internal phase, advanced once per tick
 * and the three phases of the local rule become:
 *     absorb : ANY photon in the cell, ANY direction, ANY colour -> LEVEL += 1, p += v(d).
 *              (if LEVEL is full the photon passes through untouched)
 *     emit   : on a tick with LEVEL > 0, emit ONE photon in the direction given by PHASE
 *              -- NOT by what was absorbed -- then LEVEL -= 1, p -= v(d), PHASE += 1.
 *     move   : on a tick the molecule moves one cell along the direction of its momentum
 *              (two axes alternate when the momentum is diagonal), using the momentum it had
 *              BEFORE the emission -- so a molecule that absorbed one photon moves one cell
 *              during its decay tick, exactly as in round 43.  A molecule whose momentum is
 *              zero does not move at all.
 *
 * WHAT THIS BUYS (and what the tests below are for):
 *   (a) absorption is isotropic, so photons from four different directions all pile up
 *       (LEVEL 4) with no hand-setting;
 *   (b) four cardinal photons in one cell carry ZERO net momentum, so such a molecule does
 *       not move at all -- the level is 4 and the recoil is 0.  Non-obvious, sharp, testable;
 *   (c) the emitted direction is a function of the HIDDEN PHASE, not of the absorbed photon:
 *       the same level and the same momentum give four different fates depending on the phase;
 *   (d) it can be seen with REAL photons: the same molecule excited by the same photon at two
 *       different times emits in two different directions, while everything visible at the
 *       moment of absorption was identical;
 *   (e) momentum is conserved exactly at every step (a molecule that emits keeps coasting --
 *       inertia), which is what makes the recoil bookkeeping non-trivial once directions mix.
 *
 * PREDICTIONS (written before the run)
 *   C   controls: encoding round-trips (all 256 px,py pairs); a freshly bound molecule has
 *       level 0 and p (0,0) and stays put for 40 steps; photon alone 1.000; vacuum fixed;
 *       mode-0 memory still passes then reflects.
 *   R2  the round-43 ladder still works through real photons: one photon absorbed at step 13,
 *       level 1 visible at t=13, exactly 1 emission, momentum (1,0) at every boundary AND
 *       every half step; two photons one step apart reach level 2 and give 2 emissions.
 *   R3a an already-excited molecule (level 1, p=(1,0)) absorbs a photon from ANY of the four
 *       directions (round 43 took only the axis one), the level goes to 2 each time and the
 *       momentum becomes (1,0)+v(d).
 *   R3b four cardinal photons in one cell: level 4, momentum (0,0), and the molecule does not
 *       move until it has emitted something.
 *   R4  hand-set level 1 with p=(1,0), phase 0..3: four different emitted directions
 *       (N,E,S,W), four different momenta afterwards (1,-1)/(0,0)/(1,1)/(2,0), four different
 *       paths, all four conserving momentum at every step.
 *   R5  real photons: injecting the same photon 2 cells further along makes it arrive 2 steps
 *       earlier, which is exactly one tick of internal phase => a DIFFERENT emitted direction.
 *   R6  inertia: the phase-3 case ends with level 0 and momentum (2,0) and keeps coasting.
 * ---------------------------------------------------------------------------
 */

const W = 160, H = 32, N = W * H;
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
const DIRNAME = ['N', 'E', 'S', 'W'];
const PBITS = 0xff;
const MBASE = 1 << 8;
const MODE1 = 1 << 9, MODE2 = 1 << 10;
const LVB = [1 << 11, 1 << 12, 1 << 13];        // LEVEL 0..7
const MCLK = 1 << 14;
const MREC = 1 << 15;
const PXS = 16, PYS = 20, PM = 15;              // two 4-bit signed fields (px+0, py+0)
const PHB = [1 << 24, 1 << 25];                 // PHASE 0..3
const AXISB = 1 << 26;                          // which axis the next diagonal step uses
const MARKER = ~PBITS;                          // marker bits never reach bit 31, so no sign traps
function levOf(mol) { return (mol >> 11) & 7; }
function modeOf(mol) { return (mol >> 9) & 3; }
function phOf(mol) { return (mol >> 24) & 3; }
function axisOf(mol) { return (mol >> 26) & 1; }
function sgn4(v) { return v > 7 ? v - 16 : v; }        // sign-extend a 4-bit field
function pxOf(mol) { return sgn4((mol >> PXS) & PM); }
function pyOf(mol) { return sgn4((mol >> PYS) & PM); }
function withLev(mol, n) { return (mol & ~(LVB[0] | LVB[1] | LVB[2])) | ((n & 7) << 11); }
function withPh(mol, p) { return (mol & ~(PHB[0] | PHB[1])) | ((p & 3) << 24); }
function withAxis(mol, a) { return (mol & ~AXISB) | ((a & 1) << 26); }
function withP(mol, px, py) { return (mol & ~((PM << PXS) | (PM << PYS))) | ((px & 15) << PXS) | ((py & 15) << PYS); }
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

// ---------------------------------------------------------------- one local rule, three phases
function collide(g) {
  const o = new Uint32Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const parts = v & PBITS, mol = v & MARKER;
    if (mol) {
      let m = mol ^ MCLK;                                     // the clock always ticks
      if (modeOf(mol) === 1 && parts) {                       // MODE 1 = isotropic absorber
        let lv = levOf(mol), px = pxOf(mol), py = pyOf(mol), rest = 0, took = 0;
        for (let d = 0; d < 4; d++) {
          for (let c = 0; c < 2; c++) {                       // both colours, every direction
            const bit = c ? (parts >> (d + 4)) & 1 : (parts >> d) & 1;
            if (!bit) continue;
            if (lv < 7) { lv++; px += DX[d]; py += DY[d]; took++; }
            else rest |= 1 << (d + (c ? 4 : 0));              // full: this photon passes through
          }
        }
        if (took) { o[i] = withP(withLev(m & ~MCLK, lv), px, py) | rest; continue; }  // absorption resets the phase
        o[i] = m | parts; continue;
      }
      if (!parts) { o[i] = m; continue; }
      if (m & MREC) o[i] = (m & ~MREC) | reflect(parts);       // MODE 0 memory: READ
      else o[i] = (m | MREC) | parts;                          // MODE 0 memory: WRITE
      continue;
    }
    if (TURN.has(parts)) { o[i] = TURN.get(parts); continue; }
    if (BIND.has(parts)) { o[i] = MBASE; continue; }           // a fresh molecule: p = (0,0)
    o[i] = parts;
  }
  return o;
}
// on a tick: the phase advances; if excited, emit one photon along the PHASE and recoil;
// then move one cell along the momentum direction (the momentum it had BEFORE emitting)
function moleculePhase(g) {
  const o = new Uint32Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i], mol = v & MARKER, parts = v & PBITS;
    if (!mol) { o[i] |= v; continue; }
    if (!(mol & MCLK)) { o[i] |= mol | parts; continue; }      // not a tick
    const ph = phOf(mol);
    let m = withPh(mol, (ph + 1) & 3);                        // the internal phase always runs
    let lv = levOf(mol), px = pxOf(mol), py = pyOf(mol);
    let emit = -1;
    if (lv > 0) { emit = ph; lv--; px -= DX[ph]; py -= DY[ph]; }
    // motion, using the momentum BEFORE the emission
    const mpx = px + (emit >= 0 ? DX[emit] : 0), mpy = py + (emit >= 0 ? DY[emit] : 0);
    let mx = 0, my = 0;
    if (mpx !== 0 && mpy !== 0) {
      if (axisOf(mol) === 0) mx = mpx > 0 ? 1 : -1; else my = mpy > 0 ? 1 : -1;
      m = withAxis(m, 1 - axisOf(mol));
    } else if (mpx !== 0) mx = mpx > 0 ? 1 : -1;
    else if (mpy !== 0) my = mpy > 0 ? 1 : -1;
    m = withP(withLev(m, lv), px, py);
    if (mx === 0 && my === 0) {
      o[i] |= m | parts;
      if (emit >= 0) o[i] |= 1 << emit;                       // emitted where it stands
      continue;
    }
    const dStep = mx !== 0 ? (mx > 0 ? 1 : 3) : (my > 0 ? 0 : 2);
    const x = i % W, y = (i - x) / W;
    const k = ((y + DY[dStep] + H) % H) * W + ((x + DX[dStep] + W) % W);
    o[i] |= parts;                                            // free particles stay behind
    o[k] |= m;
    if (emit >= 0) o[k] |= 1 << emit;
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
    const x = i % W;
    return { cell: i, x, y: (i - x) / W, lv: levOf(mol), px: pxOf(mol), py: pyOf(mol),
             ph: phOf(mol), axis: axisOf(mol), mode: modeOf(mol) };
  }
  return null;
}
function total(g) {
  const f = freeInfo(g), m = molInfo(g);
  return { px: f.px + (m ? m.px : 0), py: f.py + (m ? m.py : 0), f, m };
}
function freeDir(g) {
  const f = freeInfo(g);
  if (f.n !== 1) return -1;
  for (let d = 0; d < 4; d++) if (f.px === DX[d] && f.py === DY[d]) return d;
  return -2;
}
function newGrid() { return new Uint32Array(N); }
function withMolecule(mx, mode, ph0) {
  const g = newGrid();
  g[16 * W + (mx - 2)] = 1 << 1;              // colour A moving +x
  g[16 * W + (mx + 2)] = 1 << 7;              // colour B moving -x  -> they bind
  let a = Uint32Array.from(g);
  for (let t = 0; t < 3; t++) a = step(a);
  const c = 16 * W + mx;
  a[c] = withPh(withModeBits(a[c] & MARKER, mode), ph0);
  return a;
}
function withModeBits(mol, m) { return (mol & ~(MODE1 | MODE2)) | ((m & 3) << 9); }
// a hand-set control: an already-excited molecule with a chosen level, momentum and phase
function excited(mx, lv, px, py, ph, mode) {
  const g = withMolecule(mx, mode === undefined ? 1 : mode, ph);
  const c = 16 * W + mx;
  g[c] = withP(withLev(withPh(g[c] & MARKER, ph), lv), px, py);
  return g;
}
function clearPhotons(g) { for (let i = 0; i < N; i++) g[i] &= MARKER; return g; }

// ============================================================ controls
console.log('============================================================');
console.log('CONTROLS');
console.log('============================================================');
{
  let ok = true;
  for (let a = -8; a <= 7; a++) for (let b = -8; b <= 7; b++) {
    const m = withP(MBASE, a, b);
    if (pxOf(m) !== a || pyOf(m) !== b) ok = false;
  }
  const t = withP(withLev(withPh(MBASE, 2), 5), -3, 4);
  console.log(`  encoding: all 256 (px,py) pairs round-trip ${ok ? 'PASS' : '*** FAIL ***'}; level/phase/px,py = ${levOf(t)}/${phOf(t)}/${pxOf(t)},${pyOf(t)} ${levOf(t) === 5 && phOf(t) === 2 && pxOf(t) === -3 && pyOf(t) === 4 ? 'PASS' : '*** FAIL ***'}`);
}
{
  let a = withMolecule(22, 1, 0);
  const m0 = molInfo(a);
  let bad = 0;
  for (let t = 0; t < 40; t++) { a = step(a); const m = molInfo(a); if (Math.abs(m.x - 22) + Math.abs(m.y - 16) !== 0) bad++; }
  const m1 = molInfo(a), T = total(a);
  console.log(`  fresh molecule: level ${m0.lv}, p (${m0.px},${m0.py}) -> after 40 steps level ${m1.lv}, p (${m1.px},${m1.py}), x ${m1.x}, moved-on ${bad} steps   ${m0.lv === 0 && m0.px === 0 && m0.py === 0 && bad === 0 && T.px === 0 && T.py === 0 ? 'PASS (static, zero momentum)' : '*** FAIL ***'}`);
}
{
  const g = newGrid(); g[16 * W + 10] = 1 << 1;
  let a = g; for (let t = 0; t < 20; t++) a = step(a);
  const f = freeInfo(a);
  console.log(`  photon alone: x=10 -> ${f.ax} in 20 steps => speed ${((f.ax - 10) / 20).toFixed(3)}  ${Math.abs((f.ax - 10) / 20 - 1) < 1e-9 ? 'PASS' : '*** FAIL ***'}`);
}
{
  let e = newGrid(); for (let t = 0; t < 30; t++) e = step(e);
  console.log(`  empty vacuum fixed: ${freeInfo(e).n === 0 ? 'PASS' : '*** FAIL ***'}`);
}
{
  let g = withMolecule(22, 0, 0);
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
  console.log(`  mode-0 memory still: shot1 ${r1.dirAfter > 0 ? 'PASSED' : 'reflected'}, shot2 ${r2.dirAfter > 0 ? 'passed' : 'REFLECTED'}  ${r1.dirAfter > 0 && r2.dirAfter <= 0 ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ============================================================ R2: the ladder, through real photons
console.log('============================================================');
console.log('R2: the round-43 ladder must still work through real photons');
console.log('============================================================');
{
  let g = withMolecule(22, 1, 3);          // phase chosen so the decay tick lands on phase 1 (+x)
  g[16 * W + 10] = 1 << 1;
  let bad = 0, badHalf = 0, emitted = 0, prevN = 1, absorbT = -1, emitDir = -1, emitPh = -1;
  const x0 = molInfo(g).x;
  for (let t = 0; t <= 24; t++) {
    const T = total(g), cT = total(collide(g));
    if (!(T.px === 1 && T.py === 0)) bad++;
    if (!(cT.px === 1 && cT.py === 0)) badHalf++;
    if (T.f.n === 0 && absorbT < 0) absorbT = t;
    if (T.f.n > prevN) { emitted += T.f.n - prevN; if (emitDir < 0) { emitPh = T.m.ph; emitDir = freeDir(g); } }
    prevN = T.f.n;
    g = step(g);
  }
  const T = total(g);
  console.log(`  R2a one photon: absorbed at step ${absorbT} ${absorbT === 13 ? '(as predicted)' : '*** FAIL ***'}; emissions ${emitted} ${emitted === 1 ? 'PASS' : '*** FAIL ***'}`);
  const phBefore = (emitPh + 3) & 3;         // the phase that CHOSE the direction: the emission advances it
  console.log(`       emitted along ${DIRNAME[emitDir]}; the phase that chose it was ${phBefore} (it reads ${emitPh} after the step) ${emitDir === phBefore ? 'PASS: direction = phase' : '*** FAIL ***'}`);
  console.log(`       molecule x ${x0} -> ${T.m.x}, level ${T.m.lv}, p (${T.m.px},${T.m.py}); momentum not (1,0): boundary ${bad}, half-step ${badHalf}   ${T.m.x - x0 === 1 && T.m.px === 0 && T.m.py === 0 && bad === 0 && badHalf === 0 ? 'PASS (one cell of recoil, then it stops -- exactly round 43)' : '*** FAIL ***'}`);
}
{
  let g = withMolecule(22, 1, 3);
  g[16 * W + 10] = 1 << 1;
  g[16 * W + 11] = 1 << 1;                 // one step behind -> arrives 1 step later
  let maxHalf = 0, emitted = 0, prevN = 2, bad = 0;
  for (let t = 0; t <= 40; t++) {
    const T = total(g), cT = total(collide(g));
    if (!(T.px === 2 && T.py === 0)) bad++;
    if (cT.m) maxHalf = Math.max(maxHalf, cT.m.lv);
    if (T.f.n > prevN) emitted += T.f.n - prevN;
    prevN = T.f.n;
    g = step(g);
  }
  console.log(`  R2b two photons: highest level just after the absorption phase ${maxHalf} ${maxHalf === 2 ? 'PASS' : '*** FAIL ***'}; emissions ${emitted} ${emitted === 2 ? 'PASS' : '*** FAIL ***'}; momentum not (2,0) on ${bad} steps ${bad === 0 ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ============================================================ R3: absorption is isotropic now
console.log('============================================================');
console.log('R3: absorption is now isotropic (round 43 took only the axis direction)');
console.log('============================================================');
{
  console.log('  R3a an ALREADY-EXCITED molecule (level 1, p=(1,0)) meets one more photon:');
  console.log('      photon from   level after absorb   p after absorb   absorbed?');
  let okAll = true;
  for (let d = 0; d < 4; d++) {
    let g = excited(22, 1, 1, 0, 0, 1);
    const c = 16 * W + 22;
    g[c] |= 1 << d;                                        // one photon from direction d, in its cell
    const cT = total(collide(g));
    const got = cT.m.lv, px = cT.m.px, py = cT.m.py;
    const want = `${1 + DX[d]},${0 + DY[d]}`;
    const ok = got === 2 && px === 1 + DX[d] && py === 0 + DY[d];
    if (!ok) okAll = false;
    console.log(`         ${DIRNAME[d]}            ${got}                 (${px},${py})            ${ok ? 'yes (wanted ' + want + ')' : '*** NO (wanted ' + want + ') ***'}`);
  }
  console.log(`      => ${okAll ? 'PASS: all four directions are absorbed (round 43 absorbed only E)' : '*** FAIL ***'}`);
}
{
  console.log('');
  console.log('  R3b FOUR cardinal photons (N,E,S,W) in one cell:');
  let g = withMolecule(22, 1, 0);
  const c = 16 * W + 22;
  g[c] |= 0b1111;                                          // all four directions, colour A
  console.log('      t   level   p        cell(x,y)   free photons   TOTAL p');
  let emitted = 0, prevN = 4, bad = 0, maxLev = 0;      // 4 photons are on the board at t=0
  const cells = [];
  for (let t = 0; t <= 12; t++) {
    const T = total(g);
    if (!(T.px === 0 && T.py === 0)) bad++;
    if (T.m) { cells.push(`${T.m.x},${T.m.y}`); maxLev = Math.max(maxLev, T.m.lv); }
    if (T.f.n > prevN) emitted += T.f.n - prevN;
    prevN = T.f.n;
    if (t <= 10) console.log(`     ${String(t).padStart(3)}     ${T.m.lv}     (${String(T.m.px).padStart(2)},${String(T.m.py).padStart(2)})    ${String(T.m.x).padStart(3)},${String(T.m.y).padStart(2)}        ${T.f.n}          (${T.px},${T.py})`);
    g = step(g);
  }
  const T = total(g);
  console.log(`      highest level seen ${maxLev} ${maxLev === 4 ? 'PASS' : '*** FAIL ***'}, emissions ${emitted} ${emitted === 4 ? 'PASS' : '*** FAIL ***'}, final level ${T.m.lv}, final p (${T.m.px},${T.m.py}), free photons ${T.f.n} (free p (${T.f.px},${T.f.py}))`);
  console.log(`      momentum of the WHOLE thing never left (0,0): ${bad === 0 ? 'PASS' : '*** FAIL (' + bad + ' steps) ***'}; the four emitted photons together carry (${T.f.px},${T.f.py})`);
  console.log(`      => a molecule holding FOUR photons is not pushed anywhere (the four momenta cancel),`);
  console.log(`         so it does not move until it has emitted some of them.  Its cells: ${[...new Set(cells)].join(' -> ')}`);
}
console.log('');

// ============================================================ R4: the hidden phase decides
console.log('============================================================');
console.log('R4: ONE AND THE SAME STATE, FOUR DIFFERENT HIDDEN PHASES');
console.log('============================================================');
{
  console.log('  hand-set: level 1, momentum (1,0) (as if one +x photon had been absorbed), clock 0.');
  console.log('  Everything is identical except the 2-bit internal phase.');
  console.log('');
  console.log('   phase   emitted   p after emit   level at t=12   cell at t=12   momentum kept?');
  let allOk = true;
  const fates = [];
  for (let ph = 0; ph < 4; ph++) {
    let g = excited(22, 1, 1, 0, ph, 1);
    const want = `${1 - DX[ph]},${0 - DY[ph]}`;
    let emitDir = -1, emitPhBefore = -1, bad = 0, last = null;
    for (let t = 0; t <= 12; t++) {
      const T = total(g);
      if (!(T.px === 1 && T.py === 0)) bad++;
      if (t === 12) last = T;
      g = step(g);
      if (t === 0) { emitDir = freeDir(g); emitPhBefore = (molInfo(g).ph + 3) & 3; }
    }
    const T = last;
    const ok = emitDir === ph && emitDir === emitPhBefore && `${T.m.px},${T.m.py}` === want && bad === 0;
    if (!ok) allOk = false;
    fates.push(`(${T.m.x},${T.m.y})`);
    console.log(`     ${ph} (${DIRNAME[ph]})       ${DIRNAME[emitDir]}       (${String(T.m.px).padStart(2)},${String(T.m.py).padStart(2)}) [wanted ${want}]     ${T.m.lv}            ${T.m.x},${T.m.y}      ${bad === 0 ? 'yes' : '*** NO on ' + bad + ' steps ***'}`);
  }
  console.log('');
  console.log(`  => ${allOk ? 'PASS' : '*** FAIL ***'}: the same level and the same momentum, four different emitted`);
  console.log(`     directions and four different futures ${fates.join(' / ')} -- decided by the phase alone.`);
  console.log(`     (And the total momentum was (1,0) at every single step in all four runs.)`);
}
console.log('');

// ============================================================ R5: the hidden phase, with real photons
console.log('============================================================');
console.log('R5: the hidden phase shows up with REAL photons');
console.log('============================================================');
{
  console.log('  The same molecule, the same +x photon -- but launched 2 cells further along, so it');
  console.log('  arrives 2 steps earlier, which is exactly one tick of internal phase.');
  console.log('');
  console.log('   launched at   absorbed at step   phase that chose   emitted along   phase = direction?');
  const dirs = [];
  for (const x of [10, 12]) {
    let g = withMolecule(22, 1, 0);
    g[16 * W + x] = 1 << 1;
    let absorbT = -1, emitDir = -1, emitPh = -1, prevN = 1;
    for (let t = 0; t <= 30; t++) {
      const T = total(g);
      if (T.f.n === 0 && absorbT < 0) absorbT = t;
      if (T.f.n > prevN && emitDir < 0) { emitPh = (T.m.ph + 3) & 3; emitDir = freeDir(g); }
      prevN = T.f.n;
      g = step(g);
    }
    dirs.push(emitDir);
    console.log(`      x=${x}          ${String(absorbT).padStart(2)}                 ${emitPh} (${DIRNAME[emitPh]})            ${DIRNAME[emitDir]}            ${emitDir === emitPh ? 'yes' : '*** NO ***'}`);
  }
  console.log('');
  console.log(`  => ${dirs[0] !== dirs[1] && dirs[0] >= 0 ? 'PASS' : '*** FAIL ***'}: two runs that were identical in every visible respect at the moment of`);
  console.log(`     absorption ended with a photon leaving in DIFFERENT directions (${DIRNAME[dirs[0]]} vs ${DIRNAME[dirs[1]]}).`);
  console.log(`     The difference lives in the internal phase -- a hidden variable, not in anything an`);
  console.log(`     observer could have seen in the molecule's state (level 1, p=(1,0), same place, same time).`);
}
console.log('');

// ============================================================ R6: inertia
console.log('============================================================');
console.log('R6: inertia -- the recoil survives the emission');
console.log('============================================================');
{
  let g = excited(22, 1, 1, 0, 3, 1);        // phase 3 = emits -x, i.e. backwards
  console.log('   an excited molecule with p=(1,0) that emits BACKWARDS (-x): p must become (2,0)');
  console.log('   t    level   p       cell(x,y)');
  for (let t = 0; t <= 12; t++) {
    const T = total(g);
    if (t % 2 === 0 || t < 4) console.log(`   ${String(t).padStart(2)}     ${T.m.lv}    (${T.m.px},${T.m.py})   ${T.m.x},${T.m.y}`);
    g = step(g);
  }
  const T = total(g);
  console.log(`   after 12 steps: level ${T.m.lv}, p (${T.m.px},${T.m.py}), at x=${T.m.x} -- it is still moving with NO photons left.`);
  console.log(`   => ${T.m.lv === 0 && T.m.px === 2 && T.m.x > 25 ? 'PASS: level 0, momentum (2,0), still coasting at 1 cell per tick (0.500 cells/step)' : '*** FAIL ***'}`);
}
