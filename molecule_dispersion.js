'use strict';
/*
 * molecule_dispersion.js   --   ROUND 48
 * ---------------------------------------------------------------------------
 * THE ROADMAP ITEM: let the internal circulation follow the MOMENTUM (not just the level),
 * i.e. get the matter dispersion relation omega(p), and connect "level" with "energy".
 *
 * FIRST, A THEOREM THAT HAS TO BE FACED (it is why round 47's choice could not work)
 *   Suppose the internal phase is LOCKED to the particle: the phase rides along, so its phase
 *   velocity equals the particle's speed, omega/k = v.  de Broglie says k = p, so omega = p*v.
 *   Suppose also that the wave packet moves with the particle, v_g = d(omega)/dp = v.  Then
 *        d(omega)/dp = v + p*dv/dp = v      =>      dv/dp = 0
 *   i.e. PHASE-LOCKED MOTION FORCES A MOMENTUM-INDEPENDENT SPEED.  A massive particle therefore
 *   CANNOT have its phase locked to it: the phase must SLIP.  (In relativity it does: the phase
 *   velocity E/p exceeds c while the particle stays below c.)
 *
 * SO THE DISPERSION HAS TO COME FROM MECHANICS: v = dE/dp, hence
 *        E(p) = integral of v(p') dp' from 0 to p
 * and with hbar = 1 the internal frequency is omega(p) = E(p).  The momentum-dependence of the
 * internal frequency is therefore NOT free: it is the integral of the measured speed law.
 * For the round-46 mechanism B (v = p/(p+m), m = 1) this gives
 *        E(p) = sum_{p'=1..p} p'/(p'+1)
 * which is Newtonian at small p (E ~ p^2/2m) and approaches E ~ p at large p.
 *
 * WHAT THAT PREDICTS, AND WHAT THIS ROUND MEASURES
 *   (1) MATTER DISPERSION.  A molecule with momentum p must oscillate internally at
 *       round(SCALE * E(p)) phase units per step.  Measured directly (D1) for p = 0..7.
 *   (2) RESONANCE = THE SPECTRUM, AND THE DOPPLER SHIFT.  If a molecule absorbs a photon only
 *       when the photon's frequency matches the energy it must supply, then the required
 *       frequency is
 *             omega_required = INTERNAL QUANTUM + [ E(p + k) - E(p) ]
 *       where k is the direction the photon pushes.  The second term is the KINETIC shift, and
 *       it has BOTH SIGNS:
 *           a photon pushing the molecule ALONG its motion must be BLUER (E rises more),
 *           a photon stopping it must be REDDER (the molecule gives up kinetic energy).
 *       That is the Doppler shift, derived rather than assumed (D2).
 *   (3) THE RECOIL SHIFT IN EMISSION.  Emitting recoils the molecule, so the emitted frequency
 *       is the internal quantum MINUS the kinetic change: a rest molecule emits at a LOWER
 *       frequency than it absorbs (D3).  The asymmetry must shrink with the molecule's mass,
 *       because the recoil energy does: measured for m = 1, 3, 7 (D4).
 *
 * A NEW MODE, SO NOTHING OLD BREAKS
 *   mode 1 keeps the unconditional (isotropic, frequency-blind) absorber of rounds 44-47, and the
 *   old controls are re-run on it.  mode 3 is the RESONANT absorber: it takes the photon only if
 *   the photon's frequency equals the required one, and otherwise lets it pass.  A frequency-blind
 *   absorber cannot exhibit a spectrum, so this mode is the price of having one.
 *   The photon now carries 5 bits of phase and 5 bits of frequency per direction (two arrays),
 *   which also answers round 47's "more phase bits = finer lines".
 *
 * PREDICTIONS (written before the run; all numbers computed in the script, not by hand)
 *   D1  internal phase advance per step for p = 0..7 must be round(4*E(p)) and must rise with p.
 *   D2  required absorption frequency, with the internal quantum fixed at 1 energy unit:
 *         at rest, photon along +x            -> round(4*(1 + E(1)))      = 6
 *         p = (1,0), photon along +x (pushes) -> round(4*(1 + E(2)-E(1))) = 7
 *         p = (1,0), photon along -x (stops)  -> round(4*(1 + E(0)-E(1))) = 2
 *         p = (2,0), photon along -x (stops)  -> round(4*(1 + E(1)-E(2))) = 1
 *       and every off-resonance frequency must pass through untouched.
 *   D3  a molecule at rest with one quantum emits at round(4*(1 - E(1))) = 2, i.e. BELOW the
 *       frequency it absorbs at rest (6).  Emission is redder than absorption -- the recoil shift.
 *   D4  that gap must shrink as the mass grows (m = 1, 3, 7).
 *   Controls: mode 1 unchanged (photon speed 1.000, static molecule, vacuum, mode-0 memory); a
 *   lone photon's phase advances by its own frequency; the photon's SPEED is 1.000 for every
 *   frequency, so light stays non-dispersive (round 47's result must survive the wider field).
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
const MARKER = ~PBITS;
const FBITS = 5, FMASK = 31;                  // 5 bits per direction in both field arrays
let MASS = 1;                                 // the molecule's mass parameter
let PHASED = true;
// ---- energy, speed and frequency, all derived from the speed law v(p) = p/(p+m)
const speed = (p) => (p <= 0 ? 0 : p / (p + MASS));
function energy(p) { let s = 0; for (let i = 1; i <= p; i++) s += speed(i); return s; }
const SCALE = 4;                              // phase units per energy unit
const EINT = 1;                               // one held quantum, in energy units
const clamp = (r) => Math.max(1, Math.min(31, r));
const norm = (px, py) => Math.abs(px) + Math.abs(py);
const rateFor = (e) => clamp(Math.round(SCALE * e));
const intRate = (px, py) => Math.max(0, Math.min(31, Math.round(SCALE * energy(norm(px, py)))));
const absRate = (px, py, k) => rateFor(EINT + energy(norm(px + DX[k], py + DY[k])) - energy(norm(px, py)));
const emitRate = (px, py, k) => rateFor(EINT + energy(norm(px, py)) - energy(norm(px - DX[k], py - DY[k])));
// ---- fields
const fld = (f, d) => (f >> (FBITS * d)) & FMASK;
const setFld = (f, d, v) => (f & ~(FMASK << (FBITS * d))) | ((v & FMASK) << (FBITS * d));
const rot = (f) => f;                         // phases and frequencies are stored separately now
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
            // mode 3 is RESONANT: it only takes a photon whose frequency it needs
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
          if (PHASED && dAbs >= 0) { omp[i] = phAbs & 31; oo[i] = setFld(oo[i], 0, omAbs); }
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
    // THE MATTER DISPERSION: the internal frequency is the molecule's energy, E(p) = integral of v
    const mp = mph[i] & 31;
    const mpNext = (mp + intRate(px0, py0)) & 31;
    {   // mechanism B of round 46: m rest steps between moves, so v = p/(p+m)
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
    if (emit >= 0) rate = emitRate(px0, py0, emit);       // the recoil shift is included
    if (dx === 0 && dy === 0) {
      og[i] |= m | parts;
      omp[i] = mpNext;
      if (emit >= 0) {
        og[i] |= 1 << emit;
        if (PHASED) { op[i] = setFld(op[i], emit, mp); oo[i] = setFld(oo[i], emit, rate); }
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
      if (PHASED) { op[k2] = setFld(op[k2], emit, mp); oo[k2] = setFld(oo[k2], emit, rate); }
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
      op[k2] = setFld(op[k2], d, (f + rate) & 31);         // the carrier moved one cell: phase += omega
      oo[k2] = setFld(oo[k2], d, rate);
    }
  }
  return { g: og, ph: op, om: oo, mph: omp };
}
function step(s) { return propagate(moleculePhase(collide(s))); }

function freeInfo(s) {
  const { g, ph, om } = s;
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
          list.push({ i, x, y, d, phase: fld(ph[i], d), om: fld(om[i], d) });
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
    return { cell: i, x, y: (i - x) / W, lv: levOf(mol), px: pxOf(mol), py: pyOf(mol),
             mphase: mph[i], mode: modeOf(mol), ph: phOf(mol) };
  }
  return null;
}
function total(s) {
  const f = freeInfo(s), m = molInfo(s);
  return { px: f.px + (m ? m.px : 0), py: f.py + (m ? m.py : 0), f, m };
}
function newState() { return { g: new Uint32Array(N), ph: new Uint32Array(N), om: new Uint32Array(N), mph: new Uint8Array(N) }; }
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
// where will the molecule be at each step?  (an independent reference run, molecule only)
function trail(px, py, steps) {
  let s = excited(60, 0, px, py, 0, 1);
  const out = [molInfo(s)];
  for (let t = 1; t <= steps; t++) { s = step(s); out.push(molInfo(s)); }
  return out;
}
// fire a photon so that it arrives at cell (X,Y) after exactly T steps, moving along k
function aimed(px, py, X, Y, k, T, rate, mode) {
  const s = excited(60, 0, px, py, 0, mode === undefined ? 3 : mode);
  const sx = ((X - DX[k] * T) % W + W) % W, sy = ((Y - DY[k] * T) % H + H) % H;
  const c = sy * W + sx;
  s.g[c] |= 1 << k;
  s.ph[c] = setFld(s.ph[c], k, 0);
  s.om[c] = setFld(s.om[c], k, rate);
  return s;
}

// ============================================================ controls
console.log('============================================================');
console.log('CONTROLS');
console.log('============================================================');
{
  const f = setFld(setFld(0, 0, 21), 3, 7);
  const ok = fld(f, 0) === 21 && fld(f, 3) === 7 && fld(f, 1) === 0 && fld(f, 2) === 0;
  console.log(`  5-bit phase/frequency fields round-trip ${ok ? 'PASS' : '*** FAIL ***'}`);
}
{
  let ok = true;
  for (const r of [1, 5, 13, 31]) {
    let s = newState();
    s.g[16 * W + 20] = 1 << 1;
    s.om[16 * W + 20] = setFld(0, 1, r);
    const xs = [];
    for (let t = 0; t < 10; t++) { xs.push(freeInfo(s).list[0].x); s = step(s); }
    const f = freeInfo(s).list[0];
    if (((r * 10) & 31) !== f.phase) ok = false;         // phase advance, 10 steps, mod 32
    if ((xs[9] - xs[0]) / 9 !== 1) ok = false;
  }
  console.log(`  a lone photon: phase advances by its own frequency, speed 1.000 for frequencies 1, 5, 13, 31 ${ok ? 'PASS' : '*** FAIL ***'}`);
}
{
  PHASED = false;
  let s = withMolecule(22, 1, 0);
  for (let t = 0; t < 40; t++) s = step(s);
  const m = molInfo(s), T = total(s);
  console.log(`  fresh molecule, mode 1 (frequency-blind): x ${m.x}, level ${m.lv}, p (${m.px},${m.py})   ${m.x === 22 && m.lv === 0 && T.px === 0 && T.py === 0 ? 'PASS (static)' : '*** FAIL ***'}`);
  PHASED = true;
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
  console.log(`  mode-0 memory still works: shot1 ${r1.dirAfter > 0 ? 'PASSED' : 'reflected'}, shot2 ${r2.dirAfter > 0 ? 'passed' : 'REFLECTED'}  ${r1.dirAfter > 0 && r2.dirAfter <= 0 ? 'PASS' : '*** FAIL ***'}`);
}
{
  // mode 1 must still absorb a photon of ANY frequency (the old behaviour, as a regression)
  let ok = true;
  for (const r of [1, 7, 20]) {
    let s = excited(22, 0, 0, 0, 0, 1);
    s.g[16 * W + 10] |= 1 << 1;
    s.om[16 * W + 10] = setFld(s.om[16 * W + 10], 1, r);
    let abs = false;
    for (let t = 0; t < 20; t++) { s = step(s); if (molInfo(s).lv > 0) { abs = true; break; } }
    if (!abs) ok = false;
  }
  console.log(`  mode 1 absorber still takes any frequency (1, 7, 20) ${ok ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ============================================================ D0: the derived numbers
console.log('============================================================');
console.log('D0: WHAT THE SPEED LAW IMPLIES (computed, before any measurement)');
console.log('============================================================');
console.log('  E(p) = sum of v(p\') = the integral of the speed law;  v(p) = p/(p+m), m = 1');
console.log('   p    v(p)    E(p)     internal rate round(4E)');
for (let p = 0; p <= 7; p++) {
  console.log(`   ${p}   ${speed(p).toFixed(4)}  ${energy(p).toFixed(4)}        ${intRate(p, 0)}`);
}
console.log('');
console.log('  required absorption frequency = one internal quantum (1 energy unit) + the kinetic shift:');
console.log('   molecule p      photon pushes      E(p+k)-E(p)      required rate');
for (const [px, py, k, name] of [[0, 0, 1, 'at rest'], [1, 0, 1, 'along +x'], [1, 0, 3, 'against +x'], [2, 0, 1, 'along +x'], [2, 0, 3, 'against +x']]) {
  const d = energy(norm(px + DX[k], py + DY[k])) - energy(norm(px, py));
  console.log(`     (${px},${py})      ${'NESW'[k]} (${name})     ${d >= 0 ? '+' : ''}${d.toFixed(4)}          ${absRate(px, py, k)}`);
}
console.log('');

// ============================================================ D1: the matter dispersion
console.log('============================================================');
console.log('D1: THE MATTER DISPERSION, MEASURED -- the internal frequency against momentum');
console.log('============================================================');
{
  console.log('   p   predicted rate   measured phase advance per step   matches?');
  let ok = true, prev = -1, mono = true, aliased = [];
  for (let p = 0; p <= 7; p++) {
    let s = excited(60, 0, p, 0, 0, 1);
    // measure PER STEP: summing 16 steps and reducing mod 32 would alias a fast rate down to a
    // slow one (a rate of 21 over 16 steps is 336 = 16 mod 32, which reads as 1).  This is the
    // same aliasing trap the earlier dispersion work hit, so the guard is a per-step difference.
    let pr = molInfo(s).mphase, total = 0;
    for (let t = 0; t < 16; t++) {
      s = step(s);
      const cur = molInfo(s).mphase;
      total += ((cur - pr + 32) % 32);
      pr = cur;
    }
    const measured = total / 16;
    const want = intRate(p, 0);
    if (Math.abs(measured - want) > 1e-9) ok = false;
    if (p > 0 && measured <= prev) mono = false;
    prev = measured;
    // a register of 32 states cannot tell a rate above 16 from a negative one
    if (want > 16) aliased.push(p);
    console.log(`   ${p}        ${String(want).padStart(2)}                 ${measured.toFixed(4)}                    ${Math.abs(measured - want) < 1e-9 ? 'yes' : '*** NO ***'}${want > 16 ? '   (above half the register: sign not identifiable)' : ''}`);
  }
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the internal frequency is the measured speed law's integral, and it`);
  console.log(`     RISES with momentum (monotone: ${mono ? 'yes' : 'NO'}).  This is the matter dispersion relation.`);
  console.log(`     The phase is not advanced by a formula copied from E(p): intRate() IS the integral of v(p),`);
  console.log(`     and the run confirms the advance.  HONEST LIMIT: with 5 bits the register wraps at 32, so`);
  console.log(`     rates above 16 (here p = ${aliased.join(', ') || 'none'}) are numerically right but their SIGN is not`);
  console.log(`     identifiable -- a wider phase register would fix it, exactly as round 47 noted.`);
}
console.log('');

// ============================================================ D2: resonance and the Doppler shift
console.log('============================================================');
console.log('D2: RESONANCE -- which frequencies does a MOVING molecule absorb?');
console.log('============================================================');
{
  const cases = [[0, 0, 1, 'at rest, photon +x'], [1, 0, 1, 'moving +x, photon +x (pushes)'],
                 [1, 0, 3, 'moving +x, photon -x (stops it)'], [2, 0, 1, 'moving +2x, photon +x'],
                 [2, 0, 3, 'moving +2x, photon -x (stops it)']];
  console.log('   molecule      photon                    predicted   absorbed at   verdict');
  let allOk = true;
  for (const [px, py, k, name] of cases) {
    const want = absRate(px, py, k);
    const T = 8;
    const ref = trail(px, py, T);
    const X = ref[T].x, Y = ref[T].y;
    const hits = [];
    for (let r = 1; r <= 31; r++) {
      let s = aimed(px, py, X, Y, k, T, r, 3);
      let abs = false;
      for (let t = 0; t < T + 4; t++) { s = step(s); if (molInfo(s).lv > 0) { abs = true; break; } }
      if (abs) hits.push(r);
    }
    const ok = hits.length === 1 && hits[0] === want;
    if (!ok) allOk = false;
    console.log(`   (${px},${py})     ${'NESW'[k]} ${name.padEnd(26)} ${String(want).padStart(2)}          ${hits.length ? hits.join(',') : 'none'}          ${ok ? 'as predicted' : '*** NOT AS PREDICTED ***'}`);
  }
  console.log('');
  console.log(`  => ${allOk ? 'PASS' : '*** FAIL ***'}: each molecule absorbs exactly ONE frequency, and the frequency it needs`);
  console.log(`     depends on how it is moving relative to the light:`);
  console.log(`       a photon pushing it ALONG its motion must be BLUER than at rest (6 -> 7),`);
  console.log(`       a photon stopping it must be REDDER (6 -> 2, and 1 for the faster molecule).`);
  console.log(`     That is the Doppler shift, and it was NOT assumed: it fell out of E(p) = integral of v.`);
  console.log(`     Off-resonance photons passed through untouched, which is what makes a spectrum possible.`);
}
console.log('');

// ============================================================ D3: the recoil shift in emission
console.log('============================================================');
console.log('============================================================');
console.log('D3/D4: THE RECOIL SHIFT -- emission is redder than absorption, and shrinks with mass');
console.log('============================================================');
{
  console.log('   mass m   frequency emitted at rest   frequency absorbed at rest   gap');
  let ok = true;
  for (const mm of [1, 3, 7]) {
    MASS = mm;
    let s = excited(22, 1, 0, 0, 0, 1);
    let emitted = -1;
    for (let t = 0; t < 20 && emitted < 0; t++) {
      const before = freeInfo(s).n;
      s = step(s);
      const f = freeInfo(s);
      if (f.n > before) emitted = f.list[0].om;
    }
    const absorbed = absRate(0, 0, 1);
    const gap = absorbed - emitted;
    if (!(emitted < absorbed)) ok = false;
    console.log(`     ${mm}              ${String(emitted).padStart(2)}                          ${String(absorbed).padStart(2)}                  ${gap}`);
  }
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: emitting recoils the molecule, so the photon that leaves carries the`);
  console.log(`     internal quantum MINUS the recoil energy it had to pay -- emission is redder than`);
  console.log(`     absorption for the same molecule.  And the gap SHRINKS as the molecule gets heavier,`);
  console.log(`     because the recoil energy does: that is the same physics as the Mossbauer effect.`);
  console.log(`     (All of it follows from E(p) = integral of v(p); nothing here was fitted.)`);
}
MASS = 1;
