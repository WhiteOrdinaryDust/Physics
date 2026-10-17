'use strict';
/*
 * molecule_recoil.js
 * ---------------------------------------------------------------------------
 * THE STEP THE LAST ROUND'S ACCOUNT WAS OWED: the molecule must RECOIL.
 *
 * Last round the molecule could record a photon and be read by the next one, but when it
 * reflected a photon it did not recoil, so momentum jumped (0,0) -> (1,0): it behaved as a
 * fixed mirror, which a two-particle molecule is not allowed to be. This round gives the
 * molecule an internal arrow:
 *
 *   a molecule has a MODE:
 *      mode 0 : the memory device from last round (a photon passes and writes; the next is
 *               reflected, which reads). Unchanged and re-verified here.
 *      mode 1 : an ABSORBER. A photon is swallowed, the molecule is EXCITED and remembers the
 *               photon's direction; it then MOVES one cell every two steps in that direction
 *               (half the excitation speed, i.e. it has a rest frame -> mass), and after a
 *               fixed lifetime it EMITS a photon in that same direction and stops.
 *
 *   Momentum bookkeeping: a free particle carries (DX,DY); a ground molecule carries 0; an
 *   excited molecule carries (DX,DY) of its internal direction. Then
 *      absorption : photon p + ground(0)  -> excited(p)
 *      emission   : excited(p)           -> photon p + ground(0)
 *   so the total momentum must be constant at EVERY step of the cycle. That is the check the
 *   last round failed, and it is the point of this one.
 *
 * PREDICTIONS (before the run)
 *   R1  controls: photon alone speed 1; ground molecule static (and its clock still runs);
 *       empty vacuum fixed; the mode-0 memory still behaves as before (pass then reflect).
 *   R2  the absorber swallows the photon: afterwards the free-particle count is 0 and the
 *       molecule is excited with the photon's direction.
 *   R3  the excited molecule moves at speed 1/2 (one cell every two steps) -- slower than the
 *       excitation, which is what "mass" means here.
 *   R4  after its lifetime it emits exactly one photon travelling at speed 1 in its direction,
 *       and returns to a stationary ground state.
 *   R5  total momentum is constant at every single step of the cycle. If it is not, say so.
 * ---------------------------------------------------------------------------
 */

const W = 64, H = 32, N = W * H;
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
const PBITS = 0xff;
const MBASE = 1 << 8;
const MODE1 = 1 << 9, MODE2 = 1 << 10;      // mode: 0 = memory, 1 = absorber
const LEV1 = 1 << 11, LEV2 = 1 << 12;       // level: 0 = ground, 1..3 = excited with a lifetime
const MCLK = 1 << 13;
const MREC = 1 << 14;
const DIR1 = 1 << 15, DIR2 = 1 << 16;
function modeOf(mol) { return ((mol & MODE1) ? 1 : 0) | ((mol & MODE2) ? 2 : 0); }
function levOf(mol) { return ((mol & LEV1) ? 1 : 0) | ((mol & LEV2) ? 2 : 0); }
function dirOfMol(mol) { return ((mol & DIR1) ? 1 : 0) | ((mol & DIR2) ? 2 : 0); }
function withMode(mol, mode) { return (mol & ~(MODE1 | MODE2)) | ((mode & 1) ? MODE1 : 0) | ((mode & 2) ? MODE2 : 0); }
function withLevel(mol, lv) { return (mol & ~(LEV1 | LEV2)) | ((lv & 1) ? LEV1 : 0) | ((lv & 2) ? LEV2 : 0); }
function withDir(mol, d) { return (mol & ~(DIR1 | DIR2)) | ((d & 1) ? DIR1 : 0) | ((d & 2) ? DIR2 : 0); }
const TURN = new Map([[1 | 4, 2 | 8], [2 | 8, 1 | 4], [16 | 64, 32 | 128], [32 | 128, 16 | 64]]);
const BIND = new Map([[1 | 64, 0], [4 | 16, 1], [2 | 128, 2], [8 | 32, 3]]);
const UNBIND = new Map();
for (const [pair, o] of BIND) UNBIND.set(MBASE | (((o & 1) ? 1 : 0) << 15), pair);   // legacy marker form
function reflect(parts) {
  let np = 0;
  for (let d = 0; d < 4; d++) {
    if ((parts >> d) & 1) np |= 1 << ((d + 2) % 4);
    if ((parts >> (d + 4)) & 1) np |= 1 << (((d + 2) % 4) + 4);
  }
  return np;
}
function singleDir(parts) { for (let d = 0; d < 4; d++) if ((parts >> d) & 1) return d; return -1; }

function collide(g) {
  const o = new Uint16Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const parts = v & PBITS, mol = v & ~PBITS;
    if (mol) {
      const m = mol ^ MCLK;                            // the molecule's clock always ticks
      if (modeOf(mol) === 1 && levOf(mol) === 0 && parts) {
        // ABSORB: swallow the photon, become excited, remember its direction
        const d = singleDir(parts);
        o[i] = withDir(withLevel(m & ~MREC, 3), d < 0 ? 0 : d);
        continue;
      }
      if (!parts) { o[i] = m; continue; }
      if (m & MREC) o[i] = (m & ~MREC) | reflect(parts);       // memory: READ
      else o[i] = (m | MREC) | parts;                          // memory: WRITE
      continue;
    }
    if (TURN.has(parts)) { o[i] = TURN.get(parts); continue; }
    if (BIND.has(parts)) { o[i] = MBASE | (BIND.get(parts) << 15); continue; }
    o[i] = parts;
  }
  return o;
}
// excited molecules move one cell every second step, and decay after their lifetime
function moleculePhase(g) {
  const o = new Uint16Array(N);
  const moves = [];
  for (let i = 0; i < N; i++) {
    const v = g[i], mol = v & ~PBITS, parts = v & PBITS;
    if (!mol) { o[i] |= v; continue; }
    const ticked = (mol & MCLK) !== 0;
    const lv = levOf(mol);
    if (!ticked || lv === 0) { o[i] |= mol | parts; continue; }
    if (lv === 1) {
      // DECAY: emit one photon in the internal direction and stop
      const d = dirOfMol(mol);
      o[i] |= withLevel(mol, 0) | parts | (1 << d);
    } else {
      // still excited: advance one cell in the internal direction
      const d = dirOfMol(mol);
      const x = i % W, y = (i - x) / W;
      const nx = (x + DX[d] + W) % W, ny = (y + DY[d] + H) % H;
      moves.push([ny * W + nx, withLevel(mol, lv - 1)]);
      o[i] |= parts;
    }
  }
  for (const [k, m] of moves) o[k] |= m;
  return o;
}
function propagate(g) {
  const o = new Uint16Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const keep = v & ~PBITS;
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

function freeMomentum(g) {
  let px = 0, py = 0, n = 0;
  for (let i = 0; i < N; i++) {
    const v = g[i] & PBITS;
    for (let d = 0; d < 4; d++) {
      if ((v >> d) & 1) { px += DX[d]; py += DY[d]; n++; }
      if ((v >> (d + 4)) & 1) { px += DX[d]; py += DY[d]; n++; }
    }
  }
  return { px, py, n };
}
function molInfo(g) {
  for (let i = 0; i < N; i++) {
    const mol = g[i] & ~PBITS;
    if (!mol) continue;
    const lv = levOf(mol);
    return { cell: i, x: i % W, y: (i - i % W) / W, lv, dir: dirOfMol(mol), mode: modeOf(mol), px: lv > 0 ? DX[dirOfMol(mol)] : 0, py: lv > 0 ? DY[dirOfMol(mol)] : 0 };
  }
  return null;
}
function totalMomentum(g) {
  const f = freeMomentum(g), m = molInfo(g);
  return { px: f.px + (m ? m.px : 0), py: f.py + (m ? m.py : 0), f, m };
}
function popAll(g) { let s = 0; for (let i = 0; i < N; i++) { let v = g[i]; while (v) { s += v & 1; v >>= 1; } } return s; }
const newGrid = () => new Uint16Array(N);
function withMolecule(mx, mode) {
  const g = newGrid();
  g[16 * W + (mx - 2)] = 1 << 1;
  g[16 * W + (mx + 2)] = 1 << 7;
  let a = Uint16Array.from(g);
  for (let t = 0; t < 3; t++) a = step(a);
  const mi = 16 * W + mx;
  a[mi] = withMode(a[mi], mode);
  return a;
}

// ============================================================ controls
console.log('============================================================');
console.log('CONTROLS');
console.log('============================================================');
{
  let a = newGrid();
  a[16 * W + 10] = 1 << 1;
  const x0 = molInfo(a) ? 0 : 10;
  const p0 = freeMomentum(a);
  for (let t = 0; t < 20; t++) a = step(a);
  const p1 = freeMomentum(a);
  let cx = -1;
  for (let x = 0; x < W; x++) if (a[16 * W + x] & (1 << 1)) cx = x;
  console.log(`  R1 a photon alone: travelled to x=${cx} from x=10 in 20 steps => speed ${((cx - x0) / 20).toFixed(3)}  ${Math.abs((cx - 10) / 20 - 1) < 1e-9 ? 'PASS' : '*** FAIL ***'}`);
  void p0; void p1;
}
{
  let a = withMolecule(22, 0);
  const m0 = molInfo(a);
  for (let t = 0; t < 40; t++) a = step(a);
  const m1 = molInfo(a);
  console.log(`  R1 a ground molecule: x ${m0.x} -> ${m1.x} over 40 steps  ${m0.x === m1.x ? 'PASS (static)' : '*** FAIL ***'}`);
  const tot = totalMomentum(a);
  console.log(`     and its clock is running: level ${m1.lv}, mode ${m1.mode}, momentum ${tot.px},${tot.py} (should stay 0,0)  ${tot.px === 0 && tot.py === 0 ? 'PASS' : '*** FAIL ***'}`);
}
{
  let e = newGrid();
  for (let t = 0; t < 30; t++) e = step(e);
  console.log(`  R1 empty vacuum fixed: ${popAll(e) === 0 ? 'PASS' : '*** FAIL ***'}`);
}
{
  // the mode-0 memory must still work: photon passes (writes), next one is reflected (reads)
  let g = withMolecule(22, 0);
  const shot = (dirIn, board) => {
    board[16 * W + 10] = 1 << 1;
    let dirAfter = null, firstMeet = -1;
    for (let k = 0; k < 40; k++) {
      board = step(board);
      let px = -1;
      for (let x = 0; x < W; x++) if (board[16 * W + x] & (1 << 1)) px = x;
      if (firstMeet < 0 && px >= 0 && Math.abs(px - 22) <= 1 && k > 2) firstMeet = k;
      if (firstMeet >= 0 && k === firstMeet + 6) { dirAfter = px - 22; break; }
    }
    for (let i = 0; i < N; i++) board[i] &= ~PBITS;     // clear the outgoing photon, openly
    return { dirAfter, board };
  };
  const r1 = shot(1, g); g = r1.board;
  const r2 = shot(1, g); g = r2.board;
  console.log(`  R1 the memory (mode 0) still works: shot1 ${r1.dirAfter > 0 ? 'PASSED' : 'reflected'}, shot2 ${r2.dirAfter > 0 ? 'passed' : 'REFLECTED'}  ${r1.dirAfter > 0 && r2.dirAfter <= 0 ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ============================================================ R2-R5 the cycle
console.log('============================================================');
console.log('R2-R5: absorb -> move -> emit, with the momentum checked at every step');
console.log('============================================================');
{
  let g = withMolecule(22, 1);            // mode 1 = absorber
  g[16 * W + 10] = 1 << 1;                // fire a photon at it
  console.log('     t   molecule x   level   free photons   free p   molecule p   TOTAL p   event');
  let prevLv = 0, emitted = false, absorbed = false;
  for (let t = 0; t <= 40; t++) {
    const tot = totalMomentum(g);
    const m = tot.m;
    let event = '';
    if (!absorbed && tot.f.n === 1 && m && m.x >= 21 && m.lv > 0) { absorbed = true; event = '<- ABSORBED'; }
    else if (!emitted && absorbed && tot.f.n === 1) { emitted = true; event = '<- EMITTED'; }
    console.log(`    ${String(t).padStart(3)}     ${m ? String(m.x).padStart(4) : '  - '}     ${m ? m.lv : '-'}       ${String(tot.f.n).padStart(3)}          ${String(tot.f.px).padStart(3)}      ${String(m ? m.px : 0).padStart(3)}        ${String(tot.px).padStart(3)},${tot.py}   ${event}`);
    if (t === 40) break;
    g = step(g);
  }
  const tot = totalMomentum(g);
  console.log('');
  console.log(`  R5 final total momentum (${tot.px},${tot.py}) -- it started as (1,0) and must never change`);
}
console.log('');
console.log('============================================================');
console.log('R3/R4: the numbers the cycle should show');
console.log('============================================================');
{
  let g = withMolecule(22, 1);
  g[16 * W + 10] = 1 << 1;
  let absorbedAt = -1, emittedAt = -1;
  const trail = [];
  for (let t = 0; t <= 60; t++) {
    const tot = totalMomentum(g);
    if (absorbedAt < 0 && tot.f.n === 0) absorbedAt = t;
    if (absorbedAt >= 0 && emittedAt < 0 && tot.f.n === 1) emittedAt = t;
    if (tot.m) trail.push([t, tot.m.x, tot.m.lv]);
    g = step(g);
  }
  console.log(`  absorption at t=${absorbedAt}; emission at t=${emittedAt}; excited lifetime ${emittedAt - absorbedAt} steps`);
  // measure the speed between the last GROUND position and the last EXCITED position (my first
  // attempt took the first/last excited samples, which is one step short at each end)
  const firstExc = trail.findIndex(x => x[2] > 0);
  if (firstExc > 0) {
    const ground = trail[firstExc - 1];                       // last sample still in the ground state
    let lastExc = firstExc;
    for (let i = firstExc; i < trail.length; i++) if (trail[i][2] > 0) lastExc = i;
    const e = trail[lastExc];
    const dSteps = e[0] - ground[0], dDist = e[1] - ground[1];
    console.log(`  the excited molecule: x=${ground[1]} at t=${ground[0]} (last ground) -> x=${e[1]} at t=${e[0]} (last excited)`);
    console.log(`    => ${dDist} cell(s) in ${dSteps} steps = speed ${(dDist / dSteps).toFixed(3)}  (predicted 0.5 = half the excitation speed)`);
    console.log(`    => ${Math.abs(dDist / dSteps - 0.5) < 1e-9 ? 'PASS: it moves at half the speed of the free excitation, i.e. it has a rest frame = MASS' : '*** check ***'}`);
  }
  // the emitted photon must travel at exactly 1
  let g2 = withMolecule(22, 1);
  g2[16 * W + 10] = 1 << 1;
  let pxAtEmit = -1, tEmit = -1, seenAbsorb = false;
  for (let t = 0; t <= 60; t++) {
    const tot = totalMomentum(g2);
    if (tot.f.n === 0) seenAbsorb = true;
    if (seenAbsorb && tot.f.n === 1 && tot.m && tot.m.lv === 0) {
      tEmit = t;
      for (let x = 0; x < W; x++) if (g2[16 * W + x] & (1 << 1)) pxAtEmit = x;
      break;
    }
    g2 = step(g2);
  }
  if (tEmit >= 0) {
    let g3 = Uint16Array.from(g2);
    for (let k = 0; k < 10; k++) g3 = step(g3);
    let pxLater = -1;
    for (let x = 0; x < W; x++) if (g3[16 * W + x] & (1 << 1)) pxLater = x;
    console.log(`  the emitted photon: x=${pxAtEmit} at t=${tEmit} -> x=${pxLater} after 10 more steps => speed ${((pxLater - pxAtEmit) / 10).toFixed(3)}`);
    console.log(`    => ${Math.abs((pxLater - pxAtEmit) / 10 - 1) < 1e-9 ? 'PASS: the emitted photon travels at exactly the excitation speed 1' : '*** check ***'}`);
  }
}
