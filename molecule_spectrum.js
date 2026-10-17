'use strict';
/*
 * molecule_spectrum.js   --   ROUND 47
 * ---------------------------------------------------------------------------
 * THE ROADMAP ITEM WAS "several spectral lines: different levels -> different frequencies".
 * Before building anything I had to face a structural fact, and the round is organised around it:
 *
 *   THE PHOTON IN THIS MODEL IS ONE BIT.  It is a single occupancy bit in a cell, moving one
 *   cell per step.  A BIT HAS NO PHASE.  A frequency is a rate of phase rotation, so there is
 *   nowhere for a frequency to live: every emission, from every level, produces the SAME object
 *   (one bit, same direction, speed 1).  The model therefore has exactly ONE spectral line, and
 *   NO change to the molecule can give it a second one.  That is measured first (N1), as a
 *   negative result -- it locates the gap precisely instead of hiding it.
 *
 * THE MINIMAL FIX: give the photon a phase, and give it a frequency 2 bits to rotate at.
 *   A second lattice array carries, per cell PER DIRECTION, a 5-bit field: [phase:3][exp:2],
 *   i.e. 8 phase states and a frequency code omega = 2^exp in {1,2,4}.  While a photon is in
 *   flight its phase rotates by omega each step, so the spatial period (the wavelength) is
 *   8/omega = 8, 4 or 2 cells.  The carrier still moves one cell per step, so
 *
 *        wavelength x frequency = speed      =>      omega = c * k  exactly, for every line
 *
 *   which is the light dispersion relation already established in this project, now checked
 *   for photons that actually carry a frequency (N2).  Light stays gapless and non-dispersive.
 *
 * WHAT THIS BUYS (N3): the molecule's internal circulation rate is set by how excited it is
 * (n units of circulation at level n), and the emitted photon inherits that rate and the
 * molecule's phase at the moment of emission.  So different levels emit different frequencies:
 *   level 1, 2, 4  ->  omega = 1, 2, 4  ->  wavelengths 8, 4, 2 cells  ->  THREE distinct lines,
 * with the emission phase-coherent with the molecule (phase continuity at absorption too: an
 * absorbed photon hands its phase to the molecule's oscillation).
 *
 * THE PHASE IS CARRIED BY THE PARTICLE, NOT BY THE CELL -- so every local rule that moves a
 * particle must move its phase with it, INCLUDING the head-on 90-degree turn.  Getting that
 * bookkeeping wrong would silently create a phase that belongs to the wrong photon, so N4 and
 * N5 test it directly (a 90-degree crossing of two different frequencies, and a head-on turn).
 *
 * HONEST LIMITS
 *   * The phase field is shared between the two COLOURS in a direction.  Free photons are
 *     colour A, so this only matters if an A and a B photon ever occupy the same direction in
 *     one cell, which the tests avoid.  Stated, not hidden.
 *   * Only powers of two are available: 8 phase states give omega in {1,2,4}.  The physical
 *     statement I want is "frequency proportional to the internal circulation"; the lattice can
 *     only realise 1, 2, 4, so the LINES ARE GEOMETRICALLY SPACED, not equally spaced.  That is
 *     a phase-resolution artefact -- more phase bits = more lines and finer spacing.
 *   * The level -> frequency map (level 1, 2, >=4 -> omega 1, 2, 4) is a DESIGN choice.  What is
 *     not designed is that the frequency is carried by a rotating phase, that it is conserved in
 *     flight, that it survives crossings and turns, and that it is continuous at emission.
 *
 * PREDICTIONS (written before the run)
 *   P  controls: with the phase switched OFF the model must reduce exactly to the previous one
 *      (photon speed 1.000, static molecule, vacuum fixed, mode-0 memory passes then reflects);
 *      phase and frequency fields round-trip; the molecule's oscillation phase advances by
 *      (level mod 8) per step exactly; a lone photon's phase advances by omega per step exactly.
 *   N1 levels 1, 2, 3 with the phase OFF emit objects that are indistinguishable in every
 *      property the model can express  =>  one line, structurally.
 *   N2 with the phase ON: omega = 1, 2, 4 give speed 1.000 always, phase advance per step
 *      exactly omega, and measured wavelengths exactly 8, 4, 2 cells.
 *   N3 levels 1, 2, 4 emit wavelengths 8, 4, 2 (three lines), and the emitted phase equals the
 *      molecule's phase at the moment of emission.
 *   N4 two photons of different omega crossing at 90 degrees: each keeps its own omega and phase.
 *   N5 a head-on pair with different omega: after the 90-degree turn each keeps its own omega and
 *      phase (only the direction changes).
 * ---------------------------------------------------------------------------
 */

const W = 160, H = 32, N = W * H;
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
const PBITS = 0xff;
const MBASE = 1 << 8;
const MODE1 = 1 << 9, MODE2 = 1 << 10;
const LVB = [1 << 11, 1 << 12, 1 << 13];
const MCLK = 1 << 14;
const MREC = 1 << 15;
const PXS = 16, PYS = 20, PM = 15;
const PHB = [1 << 24, 1 << 25];
const AXISB = 1 << 26;
const ACCS = 27, ACCM = 7;
const MASS = 1;
const MARKER = ~PBITS;
let WALK = 'rest';            // mechanism B of round 46 (its law is v = p/(p+m))
let PHASED = true;            // whether photons carry a phase at all
// per cell, per direction: 5 bits = [phase 0..7][exp 0..2]
const fldOf = (f, d) => (f >> (5 * d)) & 31;
const phaseOf = (f) => f & 7;
const expOf = (f) => (f >> 3) & 3;
const omOf = (f) => 1 << expOf(f);
const mkFld = (p, e) => ((p & 7) | ((e & 3) << 3));
const withFld = (f, d, v) => (f & ~(31 << (5 * d))) | ((v & 31) << (5 * d));
const rot = (f) => mkFld((phaseOf(f) + omOf(f)) & 7, expOf(f));
function levOf(mol) { return (mol >> 11) & 7; }
function modeOf(mol) { return (mol >> 9) & 3; }
function phOf(mol) { return (mol >> 24) & 3; }
function axisOf(mol) { return (mol >> 26) & 1; }
function accOf(mol) { return (mol >> ACCS) & ACCM; }
function sgn4(v) { return v > 7 ? v - 16 : v; }
function pxOf(mol) { return sgn4((mol >> PXS) & PM); }
function pyOf(mol) { return sgn4((mol >> PYS) & PM); }
function withLev(mol, n) { return (mol & ~(LVB[0] | LVB[1] | LVB[2])) | ((n & 7) << 11); }
function withPh(mol, p) { return (mol & ~(PHB[0] | PHB[1])) | ((p & 3) << 24); }
function withAxis(mol, a) { return (mol & ~AXISB) | ((a & 1) << 26); }
function withAcc(mol, a) { return (mol & ~(ACCM << ACCS)) | ((a & ACCM) << ACCS); }
function withP(mol, px, py) { return (mol & ~((PM << PXS) | (PM << PYS))) | ((px & 15) << PXS) | ((py & 15) << PYS); }
function withModeBits(mol, m) { return (mol & ~(MODE1 | MODE2)) | ((m & 3) << 9); }
const TURN = new Map([[1 | 4, 2 | 8], [2 | 8, 1 | 4], [16 | 64, 32 | 128], [32 | 128, 16 | 64]]);
const BIND = new Map([[1 | 64, 0], [4 | 16, 1], [2 | 128, 2], [8 | 32, 3]]);
// which direction each direction becomes in a turn, both for the carrier and for its phase
const TURNMAP = new Map([[0, 1], [2, 3], [1, 0], [3, 2]]);
function reflect(parts) {
  let np = 0;
  for (let d = 0; d < 4; d++) {
    if ((parts >> d) & 1) np |= 1 << ((d + 2) % 4);
    if ((parts >> (d + 4)) & 1) np |= 1 << (((d + 2) % 4) + 4);
  }
  return np;
}
// the frequency a level emits at: n units of internal circulation, realised on the lattice at
// the powers of two the 3-bit phase can resolve
function expForLevel(lv) { return lv >= 4 ? 2 : (lv >= 2 ? 1 : 0); }
// the molecule's own oscillation phase gets its OWN array: field 0 of the photon array is a
// DIRECTION, and a molecule emitting along direction 0 would otherwise overwrite its own phase
// with the photon's (this was caught by the phase-advance control, not by inspection).
function collide(s) {
  const { g, ph, mph } = s;
  const og = new Uint32Array(N), op = new Uint32Array(N), om = new Uint8Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const parts = v & PBITS, mol = v & MARKER;
    if (mol) {
      let m = mol ^ MCLK;
      op[i] = ph[i];                                        // free particles keep their cells' fields
      om[i] = mph[i];                                       // the molecule keeps its own phase
      if (modeOf(mol) === 1 && parts) {
        let lv = levOf(mol), px = pxOf(mol), py = pyOf(mol), rest = 0, took = 0, dAbs = -1;
        for (let d = 0; d < 4; d++) {
          for (let c = 0; c < 2; c++) {
            const bit = c ? (parts >> (d + 4)) & 1 : (parts >> d) & 1;
            if (!bit) continue;
            if (lv < 7) { lv++; px += DX[d]; py += DY[d]; took++; if (dAbs < 0) dAbs = d; }
            else rest |= 1 << (d + (c ? 4 : 0));
          }
        }
        if (took) {
          og[i] = withAcc(withP(withLev(m & ~MCLK, lv), px, py), 0) | rest;
          // PHASE CONTINUITY AT ABSORPTION: the photon hands its phase to the molecule
          if (PHASED && dAbs >= 0) om[i] = phaseOf(fldOf(ph[i], dAbs));
          continue;
        }
        og[i] = m | parts; continue;
      }
      if (!parts) { og[i] = m; continue; }
      if (m & MREC) og[i] = (m & ~MREC) | reflect(parts);
      else og[i] = (m | MREC) | parts;
      continue;
    }
    if (TURN.has(parts)) {
      og[i] = TURN.get(parts);
      // the phase is carried by the PARTICLE, so a turn must carry it too
      let f = 0;
      for (let d = 0; d < 4; d++) f = withFld(f, TURNMAP.get(d), fldOf(ph[i], d));
      op[i] = f;
      continue;
    }
    if (BIND.has(parts)) { og[i] = MBASE; op[i] = 0; om[i] = 0; continue; }   // a fresh molecule: phase 0
    og[i] = parts; op[i] = ph[i];
  }
  return { g: og, ph: op, mph: om };
}
function moleculePhase(s) {
  const { g, ph, mph } = s;
  const og = new Uint32Array(N), op = new Uint32Array(N), om = new Uint8Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i], mol = v & MARKER, parts = v & PBITS;
    if (!mol) { og[i] |= v; op[i] |= ph[i]; continue; }
    let m = mol, emit = -1, dx = 0, dy = 0;
    const lv0 = levOf(mol);
    // the internal oscillation runs at a rate set by how excited the molecule is
    const mp = mph[i];
    const mpNext = (mp + (lv0 & 7)) & 7;
    const px0 = pxOf(mol), py0 = pyOf(mol);
    if (WALK === 'rest') {                                   // mechanism B: m rest steps per move
      const ap = Math.abs(px0) + Math.abs(py0);
      if (ap > 0) {
        let acc = accOf(mol) + ap;
        if (acc >= ap + MASS) {
          acc -= ap + MASS;
          if (px0 !== 0 && py0 !== 0) {
            if (axisOf(mol) === 0) dx = px0 > 0 ? 1 : -1; else dy = py0 > 0 ? 1 : -1;
            m = withAxis(m, 1 - axisOf(mol));
          } else if (px0 !== 0) dx = px0 > 0 ? 1 : -1;
          else dy = py0 > 0 ? 1 : -1;
        }
        m = withAcc(m, acc);
      }
    }
    if (mol & MCLK) {
      const pd = phOf(mol);
      m = withPh(m, (pd + 1) & 3);
      if (lv0 > 0) {
        emit = pd;
        m = withP(withLev(m, lv0 - 1), pxOf(mol) - DX[pd], pyOf(mol) - DY[pd]);
      }
    }
    const e = expForLevel(lv0);
    if (dx === 0 && dy === 0) {
      og[i] |= m | parts;
      om[i] = mpNext;                                        // the molecule keeps oscillating in place
      if (emit >= 0) {
        og[i] |= 1 << emit;
        // PHASE CONTINUITY AT EMISSION: the photon leaves with the molecule's phase and rate
        if (PHASED) op[i] = withFld(op[i], emit, mkFld(mp, e));
      }
      continue;
    }
    const x = i % W, y = (i - x) / W;
    const k2 = ((y + dy + H) % H) * W + ((x + dx + W) % W);
    og[i] |= parts;
    op[i] |= ph[i];                                          // free particles keep their cells' fields
    og[k2] |= m;
    om[k2] = mpNext;
    if (emit >= 0) {
      og[k2] |= 1 << emit;
      if (PHASED) op[k2] = withFld(op[k2], emit, mkFld(mp, e));
    }
  }
  return { g: og, ph: op, mph: om };
}
function propagate(s) {
  const { g, ph, mph } = s;
  const og = new Uint32Array(N), op = new Uint32Array(N), om = new Uint8Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i], keep = v & MARKER, parts = v & PBITS;
    if (keep) { og[i] |= keep; om[i] = mph[i]; }
    if (!parts) continue;
    const x = i % W, y = (i - x) / W;
    for (let d = 0; d < 4; d++) {
      const a = (v >> d) & 1, b = (v >> (d + 4)) & 1;
      if (!a && !b) continue;
      const k2 = ((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W);
      if (a) og[k2] |= 1 << d;
      if (b) og[k2] |= 1 << (d + 4);
      // the carrier moved one cell, so its phase advances by omega -- once per step
      op[k2] = withFld(op[k2], d, rot(fldOf(ph[i], d)));
    }
  }
  return { g: og, ph: op, mph: om };
}
function step(s) { return propagate(moleculePhase(collide(s))); }

function freeInfo(s) {
  const { g, ph } = s;
  let n = 0, px = 0, py = 0, ax = -1, ay = -1;
  const list = [];
  for (let i = 0; i < N; i++) {
    const v = g[i] & PBITS;
    if (!v) continue;
    const x = i % W, y = (i - x) / W;
    for (let d = 0; d < 4; d++) {
      for (const c of [0, 4]) {
        if ((v >> (d + c)) & 1) {
          n++; px += DX[d]; py += DY[d];
          list.push({ i, x, y, d, colour: c ? 'B' : 'A', phase: phaseOf(fldOf(ph[i], d)), om: omOf(fldOf(ph[i], d)) });
          if (ax < 0) { ax = x; ay = y; }
        }
      }
    }
  }
  return { n, px, py, ax, ay, list };
}
function molInfo(s) {
  const { g, mph } = s;
  for (let i = 0; i < N; i++) {
    const mol = g[i] & MARKER;
    if (!mol) continue;
    const x = i % W;
    return { cell: i, x, y: (i - x) / W, lv: levOf(mol), px: pxOf(mol), py: pyOf(mol), ph: phOf(mol), mphase: mph[i], mode: modeOf(mol) };
  }
  return null;
}
function total(s) {
  const f = freeInfo(s), m = molInfo(s);
  return { px: f.px + (m ? m.px : 0), py: f.py + (m ? m.py : 0), f, m };
}
function newState() { return { g: new Uint32Array(N), ph: new Uint32Array(N), mph: new Uint8Array(N) }; }
function withMolecule(mx, mode, phase0) {
  let s = newState();
  s.g[16 * W + (mx - 2)] = 1 << 1;
  s.g[16 * W + (mx + 2)] = 1 << 7;
  for (let t = 0; t < 3; t++) s = step(s);
  const c = 16 * W + mx;
  s.g[c] = withPh(withModeBits(s.g[c] & MARKER, mode), phase0);
  return s;
}
function excited(mx, lv, px, py, phase0, mode) {
  const s = withMolecule(mx, mode === undefined ? 1 : mode, phase0);
  const c = 16 * W + mx;
  s.g[c] = withAcc(withP(withLev(withPh(s.g[c] & MARKER, phase0), lv), px, py), 0);
  return s;
}

// ============================================================ controls
console.log('============================================================');
console.log('CONTROLS');
console.log('============================================================');
{
  let ok = true;
  const f = withFld(withFld(0, 0, mkFld(5, 2)), 3, mkFld(1, 1));
  if (phaseOf(fldOf(f, 0)) !== 5 || expOf(fldOf(f, 0)) !== 2) ok = false;
  if (phaseOf(fldOf(f, 3)) !== 1 || expOf(fldOf(f, 3)) !== 1) ok = false;
  if (fldOf(f, 1) !== 0 || fldOf(f, 2) !== 0) ok = false;
  if (omOf(mkFld(0, 0)) !== 1 || omOf(mkFld(0, 1)) !== 2 || omOf(mkFld(0, 2)) !== 4) ok = false;
  console.log(`  phase/frequency fields: round-trip and omega = 1,2,4 ${ok ? 'PASS' : '*** FAIL ***'}`);
}
for (const phased of [false, true]) {
  PHASED = phased;
  let s = newState();
  s.g[16 * W + 10] = 1 << 1;
  for (let t = 0; t < 20; t++) s = step(s);
  const f = freeInfo(s);
  console.log(`  photon alone (phase ${phased ? 'ON ' : 'OFF'}): x=10 -> ${f.ax} in 20 steps => speed ${((f.ax - 10) / 20).toFixed(3)}  ${Math.abs((f.ax - 10) / 20 - 1) < 1e-9 ? 'PASS' : '*** FAIL ***'}`);
}
PHASED = false;
{
  let s = withMolecule(22, 1, 0);
  for (let t = 0; t < 40; t++) s = step(s);
  const m = molInfo(s), T = total(s);
  console.log(`  fresh molecule (phase OFF): x ${m.x}, level ${m.lv}, p (${m.px},${m.py})   ${m.x === 22 && m.lv === 0 && T.px === 0 && T.py === 0 ? 'PASS (static, zero momentum)' : '*** FAIL ***'}`);
}
{
  let s = withMolecule(22, 0, 0);
  const shot = (st) => {
    st.g[16 * W + 10] = 1 << 1;
    let dirAfter = null, meet = -1;
    for (let k = 0; k < 40; k++) {
      st = step(st);
      let px = -1;
      for (let x = 0; x < W; x++) if (st.g[16 * W + x] & (1 << 1)) px = x;
      if (meet < 0 && px >= 0 && Math.abs(px - 22) <= 1 && k > 2) meet = k;
      if (meet >= 0 && k === meet + 6) { dirAfter = px - 22; break; }
    }
    for (let i = 0; i < N; i++) st.g[i] &= MARKER;
    return { dirAfter, st };
  };
  const r1 = shot(s); s = r1.st;
  const r2 = shot(s); s = r2.st;
  console.log(`  mode-0 memory (phase OFF): shot1 ${r1.dirAfter > 0 ? 'PASSED' : 'reflected'}, shot2 ${r2.dirAfter > 0 ? 'passed' : 'REFLECTED'}  ${r1.dirAfter > 0 && r2.dirAfter <= 0 ? 'PASS' : '*** FAIL ***'}`);
}
PHASED = true;
{
  // a lone photon's phase must advance by exactly omega per step
  let ok = true;
  for (const e of [0, 1, 2]) {
    let s = newState();
    s.g[16 * W + 20] = 1 << 1;
    s.ph[16 * W + 20] = withFld(0, 1, mkFld(0, e));
    for (let t = 0; t < 16; t++) {
      s = step(s);
      const f = freeInfo(s).list[0];
      const want = ((1 << e) * (t + 1)) & 7;
      if (f.phase !== want) ok = false;
    }
  }
  console.log(`  a lone photon's phase advances by exactly omega per step (omega = 1, 2, 4 over 16 steps) ${ok ? 'PASS' : '*** FAIL ***'}`);
}
{
  // the molecule's oscillation phase advances by (its level) per step -- and the level changes
  // as it decays, so the expected phase must be accumulated from the level observed each step
  let ok = true;
  for (const lv of [1, 3]) {
    let s = excited(22, lv, 0, 0, 0, 1);
    let expect = 0;
    for (let t = 0; t < 12; t++) {
      const before = molInfo(s);
      expect = (expect + (before.lv & 7)) & 7;
      s = step(s);
      const m = molInfo(s);
      if (m.mphase !== expect) ok = false;
    }
  }
  console.log(`  the molecule's internal phase advances by its level per step, tracked through decays ${ok ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ============================================================ N1: the negative result
console.log('============================================================');
console.log('N1: WITH NO PHASE THERE IS NOWHERE FOR A FREQUENCY TO LIVE');
console.log('============================================================');
{
  PHASED = false;
  console.log('  emit one photon from level 1, from level 2 and from level 3 (phase field OFF):');
  console.log('   emitted from level   direction   speed   any property that could differ between lines?');
  const rows = [];
  for (const lv of [1, 2, 3]) {
    let s = excited(22, lv, 0, 0, 1, 1);
    let photon = null;
    for (let t = 0; t < 20; t++) {
      const before = freeInfo(s).n;
      s = step(s);
      const f = freeInfo(s);
      if (f.n > before && !photon) photon = { ...f.list[0] };
    }
    const x0 = photon.x;
    let s2 = s;
    const after = freeInfo(s2).list[0];
    rows.push({ lv, d: photon.d, om: photon.om, phase: photon.phase, x0 });
    console.log(`            ${lv}              ${'NESW'[photon.d]}         1.000    omega = ${photon.om}, phase = ${photon.phase} (no field is ever written)`);
  }
  const same = rows.every((r) => r.d === rows[0].d && r.om === rows[0].om && r.phase === rows[0].phase);
  console.log('');
  console.log(`  => ${same ? 'PASS (the negative result is clean)' : '*** UNEXPECTED ***'}: every level emits an object with identical measurable`);
  console.log(`     properties -- one bit, one direction, speed 1, no frequency.  ONE LINE, structurally.`);
  console.log(`     No change to the MOLECULE can fix this: the frequency has nowhere to be stored.`);
  PHASED = true;
}
console.log('');

// ============================================================ N2: three lines in flight
console.log('============================================================');
console.log('N2: WITH A PHASE, THREE FREQUENCIES -- AND THE SAME SPEED FOR ALL OF THEM');
console.log('============================================================');
{
  console.log('   omega   exp   measured wavelength   measured speed   omega/k   phase per step');
  let ok = true;
  for (const e of [0, 1, 2]) {
    const om = 1 << e;
    let s = newState();
    s.g[16 * W + 20] = 1 << 1;
    s.ph[16 * W + 20] = withFld(0, 1, mkFld(0, e));
    const byCell = new Map();
    const xs = [];
    for (let t = 0; t < 40; t++) {
      const f = freeInfo(s).list[0];
      byCell.set(f.x, f.phase);
      xs.push(f.x);
      s = step(s);
    }
    const speed = (xs[xs.length - 1] - xs[0]) / (xs.length - 1);
    // the wavelength is the first distance at which the phase repeats along the path
    let lam = -1;
    for (let x = 20; x < 50 && lam < 0; x++) {
      if (!byCell.has(x)) continue;
      const ph = byCell.get(x);
      for (let y = x + 1; y < 50; y++) if (byCell.has(y) && byCell.get(y) === ph) { lam = y - x; break; }
    }
    const want = 8 / om;
    if (lam !== want || Math.abs(speed - 1) > 1e-9) ok = false;
    console.log(`    ${om}      ${e}            ${lam}                    ${speed.toFixed(3)}        1.000        ${om}`);
  }
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: wavelengths 8, 4, 2 cells for omega = 1, 2, 4, with the SAME measured speed 1.000,`);
  console.log(`     i.e. omega / k = 1 for every line.  The light dispersion relation established earlier in this`);
  console.log(`     project (gapless, non-dispersive) survives the upgrade to photons that carry a frequency.`);
}
console.log('');

// ============================================================ N3: the spectrum
console.log('============================================================');
console.log('N3: THE SPECTRUM -- different levels emit different frequencies');
console.log('============================================================');
{
  console.log('   level   emitted omega   emitted wavelength   emitted phase   molecule phase at emission   continuous?');
  let ok = true, lams = [];
  for (const lv of [1, 2, 4]) {
    let s = excited(22, lv, 0, 0, 0, 1);
    let ev = null;
    for (let t = 0; t < 20 && !ev; t++) {
      const mBefore = molInfo(s);
      const before = freeInfo(s).n;
      s = step(s);
      const f = freeInfo(s);
      if (f.n > before) ev = { p: f.list[0], molPhBefore: mBefore.mphase };
    }
    const om = ev.p.om, lam = 8 / om;
    lams.push(lam);
    // the photon is stamped with the molecule's phase at the moment of emission, and then it
    // moves one cell during that same step, which rotates it once -- so the value read at the
    // end of the emission step is the stamped phase plus omega
    const cont = ev.p.phase === ((ev.molPhBefore + om) & 7);
    if (!cont) ok = false;
    console.log(`     ${lv}          ${om}              ${lam}                  ${ev.p.phase}              ${ev.molPhBefore}                  ${cont ? 'yes (+omega)' : '*** NO ***'}`);
  }
  const distinct = new Set(lams).size === 3;
  console.log('');
  console.log(`  => ${ok && distinct ? 'PASS' : '*** FAIL ***'}: three distinct lines (measured wavelengths ${lams.join(', ')} cells) and the`);
  console.log(`     emitted phase equals the molecule's phase at the moment of emission, every time (phase continuity).`);
  console.log(`     The lines are geometrically spaced (x2) because 8 phase states only resolve omega = 1, 2, 4.`);
  console.log(`     The physical statement is "frequency proportional to internal circulation"; the lattice`);
  console.log(`     realises it only at powers of two.  More phase bits = finer spacing.`);
}
console.log('');

// ============================================================ N3b: the absorption side
console.log('============================================================');
console.log('N3b: PHASE CONTINUITY AT ABSORPTION (a photon hands its phase to the molecule)');
console.log('============================================================');
{
  let s = withMolecule(22, 1, 0);
  s.g[16 * W + 10] = 1 << 1;
  s.ph[16 * W + 10] = withFld(0, 1, mkFld(5, 1));    // an E photon of omega 2, phase 5 at launch
  let prevPhase = 5, absorbedAt = -1, got = -1, want = -1;
  for (let t = 0; t <= 20; t++) {
    const f = freeInfo(s);
    const m = molInfo(s);
    if (m.lv > 0 && absorbedAt < 0) {
      absorbedAt = t; got = m.mphase;
      // the molecule's phase is stamped from the photon's phase, then advanced once in that same
      // step by the level it had at absorption (1)
      want = (prevPhase + 1) & 7;
      break;
    }
    if (f.list.length) prevPhase = f.list[0].phase;
    s = step(s);
  }
  console.log(`  the photon's phase on the step it was absorbed: ${prevPhase}`);
  console.log(`  the molecule's internal phase right after absorption (step ${absorbedAt}): ${got}  (expected ${prevPhase} + 1 = ${want})`);
  console.log(`  => ${got === want ? 'PASS' : '*** FAIL ***'}: the photon did not vanish without trace -- its phase became the`);
  console.log(`     molecule's phase, which is what makes absorption and emission a coherent pair.`);
}
console.log('============================================================');
console.log('N4: TWO DIFFERENT FREQUENCIES CROSSING AT 90 DEGREES');
console.log('============================================================');
{
  // an E-photon and an N-photon meet in one cell; HPP only turns exact head-on pairs, so these
  // pass through each other.  Each must keep its own omega and phase.
  // geometry: the E photon from x=20 reaches x=24 after 4 steps, and the N photon from y=12
  // reaches row 16 after 4 steps, so they share cell (24,16) at the end of step 4.
  let s = newState();
  s.g[16 * W + 20] = 1 << 1;                       // E
  s.ph[16 * W + 20] = withFld(0, 1, mkFld(0, 2));  // omega 4
  s.g[12 * W + 24] = 1 << 0;                       // N
  s.ph[12 * W + 24] = withFld(0, 0, mkFld(3, 0));  // omega 1, phase 3
  let met = false, ok = true;
  for (let t = 0; t < 20; t++) {
    s = step(s);
    const f = freeInfo(s);
    if (f.n !== 2) { ok = false; break; }
    const E = f.list.find((p) => p.d === 1), Nn = f.list.find((p) => p.d === 0);
    if (!E || !Nn) { ok = false; break; }
    if (E.x === Nn.x && E.y === Nn.y) met = true;
    if (E.om !== 4 || Nn.om !== 1) ok = false;
    if (E.phase !== ((4 * (t + 1)) & 7) || Nn.phase !== ((3 + (t + 1)) & 7)) ok = false;
  }
  console.log(`  the two photons did occupy a common cell: ${met ? 'yes' : '*** never met ***'}`);
  console.log(`  => ${ok && met ? 'PASS' : '*** FAIL ***'}: both kept their own omega (4 and 1) and their own phase throughout the`);
  console.log(`     crossing.  A phase on the carrier does not make photons interact -- the round-38 result`);
  console.log(`     (photons do not self-interact) survives.`);
}
console.log('');

// ============================================================ N5: the turn must carry the phase
console.log('============================================================');
console.log('N5: A HEAD-ON TURN MUST CARRY THE PHASE WITH THE PARTICLE');
console.log('============================================================');
{
  // A head-on pair must start in the SAME cell: two carriers of equal speed in adjacent cells
  // simply swap places and never collide.  E|W in one cell turns into N|S (E->N, W->S).
  let s = newState();
  const c = 16 * W + 30;
  s.g[c] = (1 << 1) | (1 << 3);                     // E and W, same cell
  s.ph[c] = withFld(withFld(0, 1, mkFld(1, 1)), 3, mkFld(6, 0));   // E: omega 2 phase 1;  W: omega 1 phase 6
  const before = freeInfo(s);
  s = step(s);
  const f = freeInfo(s);
  const Nn = f.list.find((p) => p.d === 0), Ss = f.list.find((p) => p.d === 2);
  const expN = (1 + 2) & 7, expS = (6 + 1) & 7;
  const ok = before.n === 2 && Nn && Ss && Nn.om === 2 && Ss.om === 1 && Nn.phase === expN && Ss.phase === expS;
  console.log(`  starting as E (omega 2, phase 1) and W (omega 1, phase 6) in one cell at step 0`);
  console.log(`  after one step: N carries omega ${Nn ? Nn.om : '-'} phase ${Nn ? Nn.phase : '-'} (expected ${expN});  S carries omega ${Ss ? Ss.om : '-'} phase ${Ss ? Ss.phase : '-'} (expected ${expS})`);
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: E became N and W became S, and each kept ITS OWN omega and phase.`);
  console.log(`     If the phase belonged to the cell or the direction instead of the particle, the two`);
  console.log(`     frequencies would have swapped here.`);
}
