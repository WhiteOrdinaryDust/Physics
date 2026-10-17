'use strict';
/*
 * molecule_beam.js   --   ROUND 50
 * ---------------------------------------------------------------------------
 * THE ROADMAP ITEM: a real light field (a continuous stream of photons, no hand-aiming), and a
 * quantitative recoil limit.
 *
 * TWO THINGS CHANGE FROM ROUND 49
 *   (1) THE PHASE REGISTER IS WIDENED to 8 bits (256 states) and SCALE to 32, which is the fix
 *       rounds 47-49 kept asking for.  The frequency a head-on photon needs in order to slow a
 *       molecule from p to p-1 is round(32*(1 - v(p))), which now resolves the momentum classes:
 *            p = 6   5    4   3   2   1   (and 0, which is out of range: 48)
 *            ->  5   5    6   8   11  16
 *       (p = 6 and p = 5 still share a class -- the top of the ladder is degenerate.)
 *   (2) THE LIGHT IS A BEAM, not aimed shots: every few steps a photon of one fixed frequency is
 *       injected into EVERY row at one edge, travelling -x.  Nothing is aimed at the molecule.
 *
 * WHAT A FIXED-FREQUENCY BEAM MUST DO (derived from E(p) = integral of v, nothing added)
 *   A molecule only absorbs while its momentum is in the beam's class, and each absorption drops
 *   it one unit, out of that class.  So a beam at frequency omega cools every faster molecule down
 *   to p*(omega) - 1 and leaves it there: the BEAM FREQUENCY PICKS THE FINAL MOMENTUM.  That is
 *   velocity-selective cooling, and it is why a fixed-frequency beam does not cool to zero.
 *   Predictions: omega = 6 -> settles near 3;  omega = 8 -> near 2;  omega = 11 -> near 1;
 *   omega = 16 -> near 0;  and a frequency nobody needs (13) must do NOTHING AT ALL.
 *   After settling, the molecule is not still: the emission recoil kicks it by -1, 0 or +1, so it
 *   fluctuates in a band around the class.  That band is the recoil limit, and it is measurable.
 *
 * HONEST NOTES
 *   * The beam saturates: photons travel in one direction, so a cell holds at most one of them --
 *     the field becomes a maximal plane wave rather than growing without bound.  Injected momentum
 *     is therefore an external source, and the conservation check accounts for it explicitly.
 *   * HPP turns exact head-on pairs, so a beam photon meeting an emitted photon can scatter.  That
 *     perturbs the light, not the molecule (the recoil is applied at emission), but it is a real
 *     limitation of this photon alphabet and is stated rather than hidden.
 *   * Absolute frequencies are still conventional (SCALE and the internal quantum are choices).
 * ---------------------------------------------------------------------------
 */

const W = 80, H = 12, N = W * H;
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
const MARKER = ~PBITS;
const FBITS = 8, FMASK = 255;                 // 8 bits per direction, both arrays
let MASS = 1;
const speed = (p) => (p <= 0 ? 0 : p / (p + MASS));
function energy(p) { let s = 0; for (let i = 1; i <= p; i++) s += speed(i); return s; }
const SCALE = 32, EINT = 1;
const clamp = (r) => Math.max(1, Math.min(255, r));
const norm = (px, py) => Math.abs(px) + Math.abs(py);
const rateFor = (e) => clamp(Math.round(SCALE * e));
const intRate = (px, py) => Math.max(0, Math.min(255, Math.round(SCALE * energy(norm(px, py)))));
const absRate = (px, py, k) => rateFor(EINT + energy(norm(px + DX[k], py + DY[k])) - energy(norm(px, py)));
const emitRate = (px, py, k) => rateFor(EINT + energy(norm(px, py)) - energy(norm(px - DX[k], py - DY[k])));
const fld = (f, d) => (f >>> (FBITS * d)) & FMASK;
const setFld = (f, d, v) => ((f & ~(FMASK << (FBITS * d))) | ((v & FMASK) << (FBITS * d))) >>> 0;
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
const TURNMAP = new Map([[0, 1], [2, 3], [1, 0], [3, 2]]);
const BIND = new Map([[1 | 64, 0], [4 | 16, 1], [2 | 128, 2], [8 | 32, 3]]);
function reflect(parts) {
  let np = 0;
  for (let d = 0; d < 4; d++) {
    if ((parts >> d) & 1) np |= 1 << ((d + 2) % 4);
    if ((parts >> (d + 4)) & 1) np |= 1 << (((d + 2) % 4) + 4);
  }
  return np;
}
function collide(s) {
  const { g, ph, om, mph } = s;
  const og = new Uint32Array(N), op = new Uint32Array(N), oo = new Uint32Array(N), omp = new Uint8Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const parts = v & PBITS, mol = v & MARKER;
    if (mol) {
      let m = mol ^ MCLK;
      op[i] = ph[i]; oo[i] = om[i]; omp[i] = mph[i];
      const mode = modeOf(mol);
      if ((mode === 1 || mode === 3) && parts) {
        const px0 = pxOf(mol), py0 = pyOf(mol);
        let lv = levOf(mol), px = px0, py = py0, rest = 0, took = 0, dAbs = -1, phAbs = 0, omAbs = 0;
        for (let d = 0; d < 4; d++) {
          for (let c = 0; c < 2; c++) {
            const bit = c ? (parts >> (d + 4)) & 1 : (parts >> d) & 1;
            if (!bit) continue;
            const want = absRate(px0, py0, d);
            if (mode === 3 && fld(om[i], d) !== want) { rest |= 1 << (d + (c ? 4 : 0)); continue; }
            if (lv < 7) {
              lv++; px += DX[d]; py += DY[d]; took++;
              if (dAbs < 0) { dAbs = d; phAbs = fld(ph[i], d); omAbs = fld(om[i], d); }
            } else rest |= 1 << (d + (c ? 4 : 0));
          }
        }
        if (took) {
          og[i] = withAcc(withP(withLev(m & ~MCLK, lv), px, py), 0) | rest;
          if (dAbs >= 0) { omp[i] = phAbs & 255; oo[i] = setFld(oo[i], 0, omAbs); }
          continue;
        }
        og[i] = m | rest; continue;
      }
      if (!parts) { og[i] = m; continue; }
      if (m & MREC) og[i] = (m & ~MREC) | reflect(parts);
      else og[i] = (m | MREC) | parts;
      continue;
    }
    if (TURN.has(parts)) {
      og[i] = TURN.get(parts);
      let f = 0, o = 0;
      for (let d = 0; d < 4; d++) {
        const dst = TURNMAP.get(d);
        f = setFld(f, dst, fld(ph[i], d));
        o = setFld(o, dst, fld(om[i], d));
      }
      op[i] = f; oo[i] = o;
      continue;
    }
    if (BIND.has(parts)) { og[i] = MBASE; continue; }
    og[i] = parts; op[i] = ph[i]; oo[i] = om[i];
  }
  return { g: og, ph: op, om: oo, mph: omp };
}
function moleculePhase(s) {
  const { g, ph, om, mph } = s;
  const og = new Uint32Array(N), op = new Uint32Array(N), oo = new Uint32Array(N), omp = new Uint8Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i], mol = v & MARKER, parts = v & PBITS;
    if (!mol) { og[i] |= v; op[i] |= ph[i]; oo[i] |= om[i]; continue; }
    let m = mol, emit = -1, dx = 0, dy = 0, rate = 0;
    const lv0 = levOf(mol);
    const px0 = pxOf(mol), py0 = pyOf(mol);
    const mp = mph[i] & 255;
    const mpNext = (mp + intRate(px0, py0)) & 255;
    {
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
        m = withP(withLev(m, lv0 - 1), px0 - DX[pd], py0 - DY[pd]);
      }
    }
    if (emit >= 0) rate = emitRate(px0, py0, emit);
    if (dx === 0 && dy === 0) {
      og[i] |= m | parts;
      omp[i] = mpNext;
      if (emit >= 0) {
        og[i] |= 1 << emit;
        op[i] = setFld(op[i], emit, mp); oo[i] = setFld(oo[i], emit, rate);
      }
      continue;
    }
    const x = i % W, y = (i - x) / W;
    const k2 = ((y + dy + H) % H) * W + ((x + dx + W) % W);
    og[i] |= parts;
    op[i] |= ph[i]; oo[i] |= om[i];
    og[k2] |= m;
    omp[k2] = mpNext;
    if (emit >= 0) {
      og[k2] |= 1 << emit;
      op[k2] = setFld(op[k2], emit, mp); oo[k2] = setFld(oo[k2], emit, rate);
    }
  }
  return { g: og, ph: op, om: oo, mph: omp };
}
function propagate(s) {
  const { g, ph, om, mph } = s;
  const og = new Uint32Array(N), op = new Uint32Array(N), oo = new Uint32Array(N), omp = new Uint8Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i], keep = v & MARKER, parts = v & PBITS;
    if (keep) { og[i] |= keep; omp[i] = mph[i]; }
    if (!parts) continue;
    const x = i % W, y = (i - x) / W;
    for (let d = 0; d < 4; d++) {
      const a = (v >> d) & 1, b = (v >> (d + 4)) & 1;
      if (!a && !b) continue;
      const k2 = ((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W);
      if (a) og[k2] |= 1 << d;
      if (b) og[k2] |= 1 << (d + 4);
      const f = fld(ph[i], d), rate = fld(om[i], d);
      op[k2] = setFld(op[k2], d, (f + rate) & 255);
      oo[k2] = setFld(oo[k2], d, rate);
    }
  }
  return { g: og, ph: op, om: oo, mph: omp };
}
// THE BEAM: every `period` steps, one photon of frequency `omega` travelling -x is placed on
// every row at column `W-2`.  Nothing is aimed: the molecule meets them by crossing the field.
const BEAM_X = W - 2;
function injectBeam(s, omega) {
  for (let y = 0; y < H; y++) {
    const c = y * W + BEAM_X;
    if (s.g[c] & MARKER) continue;                   // never inject into a molecule's cell
    s.g[c] |= 1 << 3;                                // direction 3 = W (-x)
    s.om[c] = setFld(s.om[c], 3, omega);
    s.ph[c] = setFld(s.ph[c], 3, 0);
  }
  return s;
}
function step(s) { return propagate(moleculePhase(collide(s))); }

function freeInfo(s) {
  const { g, om } = s;
  let n = 0, px = 0, py = 0;
  for (let i = 0; i < N; i++) {
    const v = g[i] & PBITS;
    if (!v) continue;
    for (let d = 0; d < 4; d++) {
      for (const c of [0, 4]) {
        if ((v >> (d + c)) & 1) { n++; px += DX[d]; py += DY[d]; }
      }
    }
  }
  return { n, px, py };
}
function molInfo(s) {
  const { g, mph } = s;
  for (let i = 0; i < N; i++) {
    const mol = g[i] & MARKER;
    if (!mol) continue;
    const x = i % W;
    return { cell: i, x, y: (i - x) / W, lv: levOf(mol), px: pxOf(mol), py: pyOf(mol),
             mphase: mph[i], mode: modeOf(mol), ph: phOf(mol) };
  }
  return null;
}
function newState() { return { g: new Uint32Array(N), ph: new Uint32Array(N), om: new Uint32Array(N), mph: new Uint8Array(N) }; }
function withMolecule(mx, mode, phase0) {
  let s = newState();
  s.g[6 * W + (mx - 2)] = 1 << 1;
  s.g[6 * W + (mx + 2)] = 1 << 7;
  for (let t = 0; t < 3; t++) s = step(s);
  const c = 6 * W + mx;
  s.g[c] = withPh(withModeBits(s.g[c] & MARKER, mode), phase0);
  return s;
}
function excited(mx, lv, px, py, phase0, mode) {
  const s = withMolecule(mx, mode === undefined ? 1 : mode, phase0);
  const c = 6 * W + mx;
  s.g[c] = withAcc(withP(withLev(withPh(s.g[c] & MARKER, phase0), lv), px, py), 0);
  return s;
}
// run the beam on a molecule and record what happens.  The beam frequency is a SCHEDULE of t,
// so a fixed frequency is schedule = () => omega and a chirp is a stepped schedule.
function beamRun(p0, schedule, steps, period) {
  let s = excited(20, 0, p0, 0, 0, 3);
  const trace = [], lvTrace = [];
  let injections = 0, abs = 0, emi = 0, prevLv = 0;
  for (let t = 1; t <= steps; t++) {
    if (t % period === 0) { injectBeam(s, schedule(t)); injections++; }
    s = step(s);
    const m = molInfo(s);
    if (m.lv > prevLv) abs += m.lv - prevLv;
    if (m.lv < prevLv) emi += prevLv - m.lv;
    prevLv = m.lv;
    trace.push(m.px);
    lvTrace.push(m.lv);
  }
  return { s, trace, lvTrace, injected: injections, abs, emi };
}

// ============================================================ controls
console.log('============================================================');
console.log('CONTROLS');
console.log('============================================================');
{
  const f = setFld(setFld(0, 0, 200), 3, 17);
  const ok = fld(f, 0) === 200 && fld(f, 3) === 17 && fld(f, 1) === 0 && fld(f, 2) === 0;
  console.log(`  8-bit phase/frequency fields round-trip (including the top byte) ${ok ? 'PASS' : '*** FAIL ***'}`);
}
{
  let ok = true;
  for (const r of [1, 32, 100, 200]) {
    let s = newState();
    s.g[6 * W + 20] = 1 << 1;
    s.om[6 * W + 20] = setFld(0, 1, r);
    for (let t = 0; t < 10; t++) s = step(s);
    let x = -1;
    for (let i = 0; i < N; i++) if (s.g[i] & (1 << 1)) { const xx = i % W; if (x < 0 || xx > x) x = xx; }
    if (x - 20 !== 10) ok = false;
  }
  console.log(`  photon speed 1.000 for frequencies 1, 32, 100, 200 ${ok ? 'PASS' : '*** FAIL ***'}`);
}
{
  let s = excited(20, 0, 4, 0, 0, 1);
  let bad = 0;
  for (let t = 0; t < 40; t++) { s = step(s); if (molInfo(s).px !== 4) bad++; }
  console.log(`  an un-illuminated molecule keeps its momentum for 40 steps ${bad === 0 ? 'PASS' : '*** FAIL ***'}`);
}
{
  let s = withMolecule(20, 0, 0);
  const shotMem = (st) => {
    const y = 6;
    st.g[y * W + 8] |= 1 << 1;
    let dirAfter = null, meet = -1;
    for (let k = 0; k < 30; k++) {
      st = step(st);
      let px = -1;
      for (let x = 0; x < W; x++) if (st.g[y * W + x] & (1 << 1)) px = x;
      if (meet < 0 && px >= 0 && Math.abs(px - 20) <= 1 && k > 2) meet = k;
      if (meet >= 0 && k === meet + 5) { dirAfter = px - 20; break; }
    }
    for (let i = 0; i < N; i++) st.g[i] &= MARKER;
    return { dirAfter, st };
  };
  const r1 = shotMem(s); s = r1.st;
  const r2 = shotMem(s); s = r2.st;
  console.log(`  mode-0 memory unchanged: shot1 ${r1.dirAfter > 0 ? 'PASSED' : 'reflected'}, shot2 ${r2.dirAfter > 0 ? 'passed' : 'REFLECTED'}  ${r1.dirAfter > 0 && r2.dirAfter <= 0 ? 'PASS' : '*** FAIL ***'}`);
}
{
  // momentum bookkeeping DURING the beam: the molecule's momentum plus the photons' momentum minus
  // what the beam has injected must be constant
  let s = excited(20, 0, 4, 0, 0, 3);
  let injected = 0, bad = 0;
  let base = null;
  for (let t = 1; t <= 60; t++) {
    if (t % 4 === 0) { injectBeam(s, 8); injected += H; }
    s = step(s);
    const f = freeInfo(s), m = molInfo(s);
    const q = m.px + f.px + injected;               // injected photons carry -1 each
    if (base === null) base = q;
    if (q !== base) bad++;
  }
  console.log(`  momentum accounting through the beam (molecule + photons - injected) constant on ${60 - bad}/60 steps ${bad === 0 ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ============================================================ B0: the classes
console.log('============================================================');
console.log('B0: THE MOMENTUM CLASSES (computed, before measuring)');
console.log('============================================================');
console.log('  frequency a head-on photon needs to slow a molecule from p to p-1: round(32*(1-v(p)))');
console.log('   p     v(p)      needed frequency     a beam at this frequency settles a molecule near');
const classes = {};
for (let p = 6; p >= 0; p--) {
  const need = p === 0 ? absRate(0, 0, 3) : rateFor(EINT + energy(p - 1) - energy(p));
  if (!classes[need]) classes[need] = [];
  classes[need].push(p);
}
for (let p = 6; p >= 0; p--) {
  const need = p === 0 ? absRate(0, 0, 3) : rateFor(EINT + energy(p - 1) - energy(p));
  const settle = Math.max(0, p - 1);
  console.log(`   ${p}   ${speed(p).toFixed(4)}          ${String(need).padStart(3)}                 ${p === 0 ? '-' : settle}${p === 0 ? '   (out of range: it would push a resting molecule the other way)' : ''}`);
}
console.log(`  degenerate classes: ${Object.entries(classes).filter(([, v]) => v.length > 1).map(([k, v]) => k + ' -> p ' + v.join('/')).join(', ') || 'none'}`);
console.log('');

// ============================================================ B1: velocity selectivity
console.log('============================================================');
console.log('B1: VELOCITY SELECTIVITY -- one beam, three molecules, only one interacts');
console.log('============================================================');
{
  // A single-frequency beam is resonant with exactly one momentum class, so putting molecules of
  // different speeds in the SAME field shows the selection directly.  (An earlier version of this
  // test used frequencies that were nobody's class and "measured" nothing happening: the physics
  // was right and the experiment was wrong.)
  const omega = 8;                                  // the class of p = 3
  console.log(`   beam at frequency ${omega} (resonant with the p = 3 class).  Three separate molecules:`);
  console.log('   starting px   absorptions   mean px after 200 steps   verdict');
  let ok = true;
  for (const p0 of [6, 3, 2]) {
    const r = beamRun(p0, () => omega, 200, 4);
    const tail = r.trace.slice(150);
    const mean = tail.reduce((a, b) => a + b, 0) / tail.length;
    const interacts = (p0 === 3);
    const good = interacts ? (r.abs > 0 && Math.abs(mean - 2) <= 1.2) : (r.abs === 0 && Math.abs(mean - p0) < 1e-9);
    if (!good) ok = false;
    console.log(`       ${p0}            ${String(r.abs).padStart(2)}            ${mean.toFixed(2)}                        ${good ? (interacts ? 'cooled, as predicted' : 'untouched, as predicted') : '*** OFF ***'}`);
  }
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the same light field cools the molecule that is IN its frequency class and`);
  console.log(`     is completely transparent to the others.  That is velocity selectivity -- a fixed`);
  console.log(`     frequency picks out a velocity, which is exactly why it cannot cool everything.`);
}
console.log('');

// ============================================================ B2: a fixed frequency stalls at the class edge
console.log('============================================================');
console.log('B2: A FIXED-FREQUENCY BEAM COOLS TO THE EDGE OF ITS CLASS AND STOPS');
console.log('============================================================');
{
  const omega = 5;                                  // the class of BOTH p = 6 and p = 5
  const r = beamRun(6, () => omega, 200, 4);
  const tail = r.trace.slice(150);
  const mean = tail.reduce((a, b) => a + b, 0) / tail.length;
  console.log(`  beam at frequency ${omega}; its class contains p = 6 and p = 5 and nothing lower.`);
  console.log(`  molecule starts at px = 6 -> absorptions ${r.abs}, emissions ${r.emi}, px over the last 10 steps ${r.trace.slice(-10).join(',')}`);
  console.log(`  => ${Math.abs(mean - 4) <= 1.2 && r.abs >= 1 ? 'PASS' : '*** FAIL ***'}: it is driven to the EDGE of the class (px = 4) and stops there, because the`);
  console.log(`     next step down needs frequency 6.  The absorption count is 1 or 2 depending on whether an`);
  console.log(`     emission recoil happens to help it out of the class -- the momentum is the reproducible`);
  console.log(`     quantity, not the number of absorptions.`);
}
console.log('');

// ============================================================ B3: the chirp walks it down
console.log('============================================================');
console.log('B3: A CHIRPED BEAM WALKS THE MOLECULE ALL THE WAY DOWN');
console.log('============================================================');
{
  const ramp = (list) => (t) => list[Math.min(list.length - 1, Math.floor((t - 1) / 25))];
  const full = beamRun(6, ramp([5, 6, 8, 11, 16]), 260, 4);
  const half = beamRun(6, ramp([5, 6, 8]), 200, 4);
  const tF = full.trace.slice(220), tH = half.trace.slice(160);
  const mF = tF.reduce((a, b) => a + b, 0) / tF.length;
  const mH = tH.reduce((a, b) => a + b, 0) / tH.length;
  console.log(`   schedule 5,6,8,11,16 (one step every 25) -> final mean px ${mF.toFixed(2)}, absorptions ${full.abs}`);
  console.log(`   schedule 5,6,8 (stops early)             -> final mean px ${mH.toFixed(2)}, absorptions ${half.abs}`);
  console.log(`  => ${mF < 1.0 && Math.abs(mH - 2) <= 1.2 ? 'PASS' : '*** FAIL ***'}: the chirp walks the molecule down the ladder, and the LAST frequency the beam`);
  console.log(`     reaches is what decides where it stops (16 -> near 0, stopping at 8 -> near 2).`);
  console.log(`     The molecule waits for the frequency to catch it rather than falling behind, so the chirp`);
  console.log(`     does not have to be slow -- the final frequency is the thing that matters.`);
}
console.log('');

// ============================================================ B4: the recoil limit
console.log('============================================================');
console.log('B4: THE RECOIL LIMIT -- the band it settles into, measured');
console.log('============================================================');
{
  // the ladder has to START at the molecule's own class, otherwise the field is simply transparent
  // (the mistake this test made in its first version: it opened at 11 while the molecule sat at 6)
  const r = beamRun(6, (t) => (t < 80 ? 5 : t < 140 ? 6 : t < 200 ? 8 : t < 260 ? 11 : 16), 460, 4);
  const stats = (a) => {
    const m = a.reduce((x, y) => x + y, 0) / a.length;
    return { m, sd: Math.sqrt(a.reduce((x, y) => x + (y - m) ** 2, 0) / a.length), mn: Math.min(...a), mx: Math.max(...a) };
  };
  const mid = stats(r.trace.slice(60, 300));      // while absorptions and emissions still alternate
  const tail = stats(r.trace.slice(310));         // after the cycling has stopped
  console.log(`  the beam walks the whole ladder (5,6,8,11,16) and then parks at 16, 460 steps:`);
  console.log(`   while it is still cycling (steps 60-300): mean px ${mid.m.toFixed(2)}, spread ${mid.sd.toFixed(2)}, range ${mid.mn}..${mid.mx}`);
  console.log(`   after it has settled     (steps 310-460): mean px ${tail.m.toFixed(2)}, spread ${tail.sd.toFixed(2)}, range ${tail.mn}..${tail.mx}`);
  console.log(`   events over the whole run: ${r.abs} absorptions, ${r.emi} emissions`);
  console.log('');
  console.log(`  => MY PREDICTION HERE WAS WRONG, and the measurement is more interesting than the guess.`);
  console.log(`     There IS a jitter while the molecule is cycling (${mid.sd.toFixed(2)} units, from the emission recoils),`);
  console.log(`     but once it is at rest it has nothing left to emit AND the beam is far off resonance: a`);
  console.log(`     resting molecule would need frequency 48 and the beam is at 16.  So the kicks stop with`);
  console.log(`     the cycling and px sits at exactly 0: this model reaches a PERFECT stop, not a recoil`);
  console.log(`     limit.  Real laser cooling keeps cycling even at rest because the light is close to`);
  console.log(`     resonance there -- the recoil limit lives in that near-resonant regime, and this coarse`);
  console.log(`     frequency ladder (classes 5, 6, 8, 11, 16, then 48) cannot express it.`);
}
console.log('');

// ============================================================ B5: a frequency nobody needs
console.log('============================================================');
console.log('B5: A FREQUENCY NO MOMENTUM CLASS NEEDS DOES NOTHING (now a real control)');
console.log('============================================================');
{
  const r = beamRun(6, () => 13, 160, 4);
  const tail = r.trace.slice(120);
  const mean = tail.reduce((a, b) => a + b, 0) / tail.length;
  console.log(`  beam at frequency 13 (the classes are ${Object.keys(classes).join(', ')}): absorptions ${r.abs}, mean px ${mean.toFixed(3)}`);
  console.log(`  => ${r.abs === 0 && Math.abs(mean - 6) < 1e-9 ? 'PASS' : '*** FAIL ***'}: 13 is nobody's class, so the field is transparent -- including in the co-moving`);
  console.log(`     direction, where the needed frequencies are much higher.  This control only means something`);
  console.log(`     because B1/B2/B3 show the beam DOES work at the right frequencies: an earlier version of`);
  console.log(`     this test "passed" while the beam was doing nothing at all.`);
}
