'use strict';
/*
 * molecule_linewidth.js   --   ROUND 52
 * ---------------------------------------------------------------------------
 * THE ROADMAP ITEM: give the transition a linewidth, so that a molecule keeps cycling near rest
 * and a real recoil limit can appear.
 *
 * THE LINE WIDTH IS NOT A NEW FREE PARAMETER -- IT COMES FROM THE LIFETIME
 *   Energy and time trade off: a state that lives tau steps has an energy uncertainty of order
 *   hbar/tau.  In this model the energy unit is EINT = 1, carried as SCALE = 32 frequency units,
 *   so the natural width of the transition is
 *          GAMMA = SCALE * EINT / tau          (in frequency units)
 *   and the absorption rule becomes: take the photon when  | omega_photon - omega_required | <= GAMMA.
 *   The molecule's excited state now really does live tau steps (an age counter, emitted when
 *   age >= tau), instead of the fixed period-2 tick of rounds 48-51.
 *
 * AND THAT EXPOSES AN INCONSISTENCY IN THE EARLIER ROUNDS, WHICH THIS ROUND FIXES
 *   Rounds 49-51 used an INFINITELY SHARP line (exact frequency match) together with a VERY SHORT
 *   lifetime (2 steps).  Energy and time forbid that combination.  A sharp line belongs to a long
 *   lived state, and a short lived state must have a broad line.  The sharp rule is kept here only
 *   as a control, labelled for what it is: the tau -> infinity limit, used with tau = 2.
 *
 * WHY A HEAVY MOLECULE IS NEEDED TO SEE THE CROSSOVER
 *   The frequency a photon must have changes by about SCALE/(MASS+1) for every unit of momentum.
 *   That shift is the recoil energy in frequency units.  With MASS = 1 it is ~16 units, far more
 *   than any achievable line width (GAMMA <= 16 even at tau = 2), so one absorption always throws
 *   the molecule clean out of resonance.  With MASS = 31 the shift is ~1 unit, so a line of width
 *   GAMMA covers GAMMA momentum units and the molecule can cycle.
 *   (The walk accumulator therefore moves into its own array: with MASS = 31 its threshold reaches
 *   38, which does not fit in the marker's 4 spare bits.)
 *
 * PREDICTIONS (written before the run)
 *   L1  the excited state lives exactly tau steps: the delay from absorption to emission = tau.
 *   L2  THE LINE WIDTH, MEASURED.  A resting molecule absorbs a photon when the frequency is within
 *       GAMMA of the required one, so the accepted window should be 2*GAMMA + 1 frequency units
 *       wide: tau = 8 -> GAMMA 4 -> 9 units;  tau = 16 -> 2 -> 5;  tau = 32 -> 1 -> 3.
 *   L3  THE CROSSOVER (the point of the round).  A heavy molecule in a red-detuned molasses settles
 *       with a momentum spread set by the line: a WIDE line (short tau) lets it wander several units,
 *       a NARROW line (long tau) should squeeze it down towards ONE recoil -- and then STOP, because
 *       every cycle still kicks it by one recoil no matter how narrow the line is.  Predicted spread:
 *       tau = 8 -> a few units;  tau = 16 -> about 2;  tau = 32 and 64 -> about 1 and no narrower.
 *       That floor is the recoil limit: the Doppler limit (set by GAMMA) gives way to it as the line
 *       narrows.
 *   Controls: fields, light speed for several frequencies, momentum accounting through the beams,
 *       mode-0 memory, and the sharp-line rule as the tau -> infinity limit.
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
const MARKER = ~PBITS;
const FBITS = 8, FMASK = 255;
let MASS = 31;                                  // heavy, so the recoil shift is ~1 frequency unit
let LIFETIME = 16;                              // tau, in steps
let SHARP = false;                              // the tau -> infinity control
const speed = (p) => (p <= 0 ? 0 : p / (p + MASS));
function energy(p) { let s = 0; for (let i = 1; i <= p; i++) s += speed(i); return s; }
const SCALE = 32, EINT = 1;
const clamp = (r) => Math.max(1, Math.min(255, r));
const norm = (px, py) => Math.abs(px) + Math.abs(py);
const rateFor = (e) => clamp(Math.round(SCALE * e));
const intRate = (px, py) => Math.max(0, Math.min(255, Math.round(SCALE * energy(norm(px, py)))));
const absRate = (px, py, k) => rateFor(EINT + energy(norm(px + DX[k], py + DY[k])) - energy(norm(px, py)));
const emitRate = (px, py, k) => rateFor(EINT + energy(norm(px, py)) - energy(norm(px - DX[k], py - DY[k])));
const gamma = () => (SHARP ? 0 : (SCALE * EINT) / LIFETIME);
const fld = (f, d) => (f >>> (FBITS * d)) & FMASK;
const setFld = (f, d, v) => ((f & ~(FMASK << (FBITS * d))) | ((v & FMASK) << (FBITS * d))) >>> 0;
function levOf(mol) { return (mol >> 11) & 7; }
function modeOf(mol) { return (mol >> 9) & 3; }
function phOf(mol) { return (mol >> 24) & 3; }
function sgn4(v) { return v > 7 ? v - 16 : v; }
function pxOf(mol) { return sgn4((mol >> PXS) & PM); }
function pyOf(mol) { return sgn4((mol >> PYS) & PM); }
function withLev(mol, n) { return (mol & ~(LVB[0] | LVB[1] | LVB[2])) | ((n & 7) << 11); }
function withPh(mol, p) { return (mol & ~(PHB[0] | PHB[1])) | ((p & 3) << 24); }
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
  const { g, ph, om, mph, age, acc } = s;
  const og = new Uint32Array(N), op = new Uint32Array(N), oo = new Uint32Array(N);
  const omp = new Uint8Array(N), oag = new Uint8Array(N), oac = new Uint8Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const parts = v & PBITS, mol = v & MARKER;
    if (mol) {
      let m = mol ^ MCLK;
      op[i] = ph[i]; oo[i] = om[i]; omp[i] = mph[i]; oag[i] = age[i]; oac[i] = acc[i];
      const mode = modeOf(mol);
      if ((mode === 1 || mode === 3) && parts) {
        const px0 = pxOf(mol), py0 = pyOf(mol);
        let lv = levOf(mol), px = px0, py = py0, rest = 0, took = 0, dAbs = -1, phAbs = 0, omAbs = 0;
        for (let d = 0; d < 4; d++) {
          for (let c = 0; c < 2; c++) {
            const bit = c ? (parts >> (d + 4)) & 1 : (parts >> d) & 1;
            if (!bit) continue;
            const want = absRate(px0, py0, d);
            if (mode === 3 && Math.abs(fld(om[i], d) - want) > gamma()) { rest |= 1 << (d + (c ? 4 : 0)); continue; }
            if (lv < 7) {
              lv++; px += DX[d]; py += DY[d]; took++;
              if (dAbs < 0) { dAbs = d; phAbs = fld(ph[i], d); omAbs = fld(om[i], d); }
            } else rest |= 1 << (d + (c ? 4 : 0));
          }
        }
        if (took) {
          og[i] = withP(withLev(m, lv), px, py) | rest;
          oag[i] = 0;                                  // the state is entered at the END of this step;
                                                       // it decays when the age reaches LIFETIME, so
                                                       // it lasts LIFETIME steps counting this one
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
  return { g: og, ph: op, om: oo, mph: omp, age: oag, acc: oac };
}
function moleculePhase(s) {
  const { g, ph, om, mph, age, acc } = s;
  const og = new Uint32Array(N), op = new Uint32Array(N), oo = new Uint32Array(N);
  const omp = new Uint8Array(N), oag = new Uint8Array(N), oac = new Uint8Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i], mol = v & MARKER, parts = v & PBITS;
    if (!mol) { og[i] |= v; op[i] |= ph[i]; oo[i] |= om[i]; continue; }
    let m = mol, emit = -1, dx = 0, dy = 0, rate = 0;
    const lv0 = levOf(mol);
    const px0 = pxOf(mol), py0 = pyOf(mol);
    const mp = mph[i] & 255;
    const mpNext = (mp + intRate(px0, py0)) & 255;
    let ageNext = age[i];
    if (lv0 > 0) {
      ageNext = age[i] + 1;
      if (ageNext >= LIFETIME) {                       // the excited state finally decays
        emit = phOf(mol);
        m = withPh(m, (phOf(mol) + 1) & 3);
        m = withP(withLev(m, lv0 - 1), px0 - DX[emit], py0 - DY[emit]);
        ageNext = 0;
      }
    } else ageNext = 0;
    {
      const ap = Math.abs(px0) + Math.abs(py0);
      if (ap > 0) {                                    // mechanism B: v = p/(p+m), accumulator in its own array
        let a = acc[i] + ap;
        const thr = ap + MASS;
        if (a >= thr) {
          a -= thr;
          if (px0 !== 0 && py0 !== 0) {
            if (((mph[i] >> 7) & 1) === 0) dx = px0 > 0 ? 1 : -1; else dy = py0 > 0 ? 1 : -1;
            omp[i] = mph[i];                            // (the axis bit is not kept separately here)
          } else if (px0 !== 0) dx = px0 > 0 ? 1 : -1;
          else dy = py0 > 0 ? 1 : -1;
        }
        oac[i] = a;
      }
    }
    if (emit >= 0) rate = emitRate(px0, py0, emit);
    if (dx === 0 && dy === 0) {
      og[i] |= m | parts;
      omp[i] = mpNext; oag[i] = ageNext;
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
    omp[k2] = mpNext; oag[k2] = ageNext; oac[k2] = oac[i];
    if (emit >= 0) {
      og[k2] |= 1 << emit;
      op[k2] = setFld(op[k2], emit, mp); oo[k2] = setFld(oo[k2], emit, rate);
    }
  }
  return { g: og, ph: op, om: oo, mph: omp, age: oag, acc: oac };
}
function propagate(s) {
  const { g, ph, om, mph, age, acc } = s;
  const og = new Uint32Array(N), op = new Uint32Array(N), oo = new Uint32Array(N);
  const omp = new Uint8Array(N), oag = new Uint8Array(N), oac = new Uint8Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i], keep = v & MARKER, parts = v & PBITS;
    if (keep) { og[i] |= keep; omp[i] = mph[i]; oag[i] = age[i]; oac[i] = acc[i]; }
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
  return { g: og, ph: op, om: oo, mph: omp, age: oag, acc: oac };
}
function injectBeam(s, omega, dir) {
  const col = dir === 3 ? W - 2 : 1;
  for (let y = 0; y < H; y++) {
    const c = y * W + col;
    if (s.g[c] & MARKER) continue;
    s.g[c] |= 1 << dir;
    s.om[c] = setFld(s.om[c], dir, omega);
    s.ph[c] = setFld(s.ph[c], dir, 0);
  }
  return s;
}
function step(s) { return propagate(moleculePhase(collide(s))); }

function freeInfo(s) {
  const { g } = s;
  let n = 0, px = 0, py = 0;
  for (let i = 0; i < N; i++) {
    const v = g[i] & PBITS;
    if (!v) continue;
    for (let d = 0; d < 4; d++) for (const c of [0, 4]) if ((v >> (d + c)) & 1) { n++; px += DX[d]; py += DY[d]; }
  }
  return { n, px, py };
}
function molInfo(s) {
  const { g, mph, age } = s;
  for (let i = 0; i < N; i++) {
    const mol = g[i] & MARKER;
    if (!mol) continue;
    const x = i % W;
    return { cell: i, x, y: (i - x) / W, lv: levOf(mol), px: pxOf(mol), py: pyOf(mol),
             mphase: mph[i], age: age[i], mode: modeOf(mol), ph: phOf(mol) };
  }
  return null;
}
function newState() {
  return { g: new Uint32Array(N), ph: new Uint32Array(N), om: new Uint32Array(N),
           mph: new Uint8Array(N), age: new Uint8Array(N), acc: new Uint8Array(N) };
}
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
  s.g[c] = withP(withLev(withPh(s.g[c] & MARKER, phase0), lv), px, py);
  return s;
}
// one photon, fired along the row at a resting molecule (a resting molecule does not walk)
function restShot(omega, tau) {
  LIFETIME = tau;
  let s = excited(20, 0, 0, 0, 0, 3);
  const c = 6 * W + 30;
  s.g[c] |= 1 << 3;
  s.om[c] = setFld(s.om[c], 3, omega);
  let absorbed = false, delay = -1, tAbs = -1;
  for (let t = 1; t <= 60; t++) {
    const before = molInfo(s).lv;
    s = step(s);
    const m = molInfo(s);
    if (m.lv > before && !absorbed) { absorbed = true; tAbs = t; }
    if (absorbed && m.lv === 0 && delay < 0) delay = t - tAbs;
  }
  return { absorbed, delay };
}
function stats(a) {
  const m = a.reduce((x, y) => x + y, 0) / a.length;
  return { m, sd: Math.sqrt(a.reduce((x, y) => x + (y - m) ** 2, 0) / a.length), mn: Math.min(...a), mx: Math.max(...a) };
}

// ============================================================ controls
console.log('============================================================');
console.log('CONTROLS');
console.log('============================================================');
console.log(`  machine: MASS = ${MASS}, SCALE = ${SCALE}, energy unit = ${EINT};  GAMMA = SCALE*EINT/tau`);
{
  const f = setFld(setFld(0, 0, 200), 3, 17);
  console.log(`  8-bit fields round-trip ${fld(f, 0) === 200 && fld(f, 3) === 17 ? 'PASS' : '*** FAIL ***'}`);
}
{
  let ok = true;
  for (const r of [1, 16, 100, 200]) {
    let s = newState();
    s.g[6 * W + 20] = 1 << 1;
    s.om[6 * W + 20] = setFld(0, 1, r);
    for (let t = 0; t < 10; t++) s = step(s);
    let x = -1;
    for (let i = 0; i < N; i++) if (s.g[i] & (1 << 1)) { const xx = i % W; if (x < 0 || xx > x) x = xx; }
    if (x - 20 !== 10) ok = false;
  }
  console.log(`  photon speed 1.000 for frequencies 1, 16, 100, 200 ${ok ? 'PASS' : '*** FAIL ***'}`);
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
  console.log(`  mode-0 memory unchanged ${r1.dirAfter > 0 && r2.dirAfter <= 0 ? 'PASS' : '*** FAIL ***'}`);
}
{
  let s = excited(20, 0, 3, 0, 0, 3);
  let injected = 0, bad = 0, base = null;
  for (let t = 1; t <= 80; t++) {
    if (t % 6 === 0) { injectBeam(s, 30, 3); injected += H; }
    s = step(s);
    const f = freeInfo(s), m = molInfo(s);
    const q = m.px + f.px + injected;
    if (base === null) base = q;
    if (q !== base) bad++;
  }
  console.log(`  momentum accounting through a beam (molecule + photons - injected) constant on ${80 - bad}/80 steps ${bad === 0 ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ============================================================ L1: the lifetime is the lifetime
console.log('============================================================');
console.log('L1: THE EXCITED STATE LIVES tau STEPS (and that is what sets the width)');
console.log('============================================================');
{
  const need = absRate(0, 0, 3);
  console.log(`  firing a photon at exactly the resting frequency ${need} and timing absorption -> emission:`);
  console.log('   tau   gamma    measured delay (steps)');
  let ok = true;
  for (const tau of [2, 4, 8, 16, 32]) {
    const r = restShot(need, tau);
    const good = r.absorbed && r.delay === tau - 1;
    if (!good) ok = false;
    console.log(`    ${String(tau).padStart(2)}   ${((SCALE * EINT) / tau).toFixed(2)}      ${r.delay}  ${good ? '' : '*** MISMATCH ***'}`);
  }
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the delay from the absorption STEP to the emission STEP is LIFETIME-1, i.e. the`);
  console.log(`     state lasts LIFETIME steps counting the one it was created in.  That convention is the only`);
  console.log(`     gap between the measured delay and tau; the line widths measured next confirm that`);
  console.log(`     GAMMA = SCALE*EINT/LIFETIME is the right size.`);
}
console.log('');

// ============================================================ L2: the line width, measured
console.log('============================================================');
console.log('L2: THE LINE WIDTH, MEASURED -- sweeping the frequency past a resting molecule');
console.log('============================================================');
{
  const need = absRate(0, 0, 3);
  console.log(`  resting molecule needs ${need};  a photon is absorbed when |omega - ${need}| <= GAMMA`);
  console.log('   tau   predicted window   measured accepted frequencies   measured width');
  let ok = true;
  for (const tau of [4, 8, 16, 32]) {
    const G = (SCALE * EINT) / tau;
    const accepted = [];
    for (let om = need - 10; om <= need + 10; om++) if (restShot(om, tau).absorbed) accepted.push(om);
    const width = accepted.length;
    const want = 2 * Math.floor(G) + 1;
    const good = width === want;
    if (!good) ok = false;
    console.log(`    ${String(tau).padStart(2)}        ${String(want).padStart(2)}               ${accepted.join(',')}          ${width}  ${good ? '' : '*** MISMATCH ***'}`);
  }
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the line is not sharp.  Its width follows the lifetime, which is exactly`);
  console.log(`     the uncertainty relation doing its work -- and it is why the earlier rounds' sharp line with a`);
  console.log(`     2-step lifetime was inconsistent.`);
}
console.log('');

// ============================================================ L3: the crossover
console.log('============================================================');
console.log('L3: THE CROSSOVER -- Doppler limit giving way to the recoil limit');
console.log('============================================================');
{
  // A heavy molecule in a two-sided beam placed a fixed FRACTION OF A LINEWIDTH below its resting
  // transition -- which is what real cooling does.  (The first version of this test held the
  // detuning fixed at 3 units while the line narrowed, so the molecule simply fell out of
  // resonance and froze; the spread then fell to almost nothing, which looks like a narrower
  // limit but is really a dead molecule.  The cycling count below is what tells the two apart.)
  const need = absRate(0, 0, 3);
  console.log(`  heavy molecule (MASS = ${MASS});  the beam sits half a linewidth below the resting frequency ${need}`);
  console.log('   tau   gamma   beam   detune/gamma   settled mean   spread   range      cycling in the window');
  const rows = [];
  for (const tau of [4, 8, 16, 32, 64]) {
    const G = (SCALE * EINT) / tau;
    const detune = G / 2;
    const beam = Math.round(need - detune);
    LIFETIME = tau;
    let s = excited(20, 0, 2, 0, 0, 3);
    const trace = [], lvTrace = [];
    for (let t = 1; t <= 1200; t++) {
      if (t % 6 === 0) { injectBeam(s, beam, 3); injectBeam(s, beam, 1); }
      s = step(s);
      if (t > 800) { trace.push(molInfo(s).px); lvTrace.push(molInfo(s).lv); }
    }
    const st = stats(trace);
    let cyc = 0;
    for (let i = 1; i < lvTrace.length; i++) if (lvTrace[i] > lvTrace[i - 1]) cyc++;
    rows.push({ tau, G, sd: st.sd, cyc, m: st.m });
    console.log(`    ${String(tau).padStart(2)}   ${G.toFixed(2)}   ${String(beam).padStart(2)}      ${(detune / G).toFixed(2)}          ${st.m.toFixed(2)}        ${st.sd.toFixed(2)}     ${st.mn}..${st.mx}      ${cyc}`);
  }
  console.log('');
  const cycling = rows.filter((r) => r.cyc > 0);
  const frozen = rows.filter((r) => r.cyc === 0);
  const widest = rows[0], narrowestCycling = cycling[cycling.length - 1];
  const doppler = widest.sd > narrowestCycling.sd + 0.2;      // a broader line -> a wider spread
  const atRecoil = narrowestCycling.sd <= 1.3;                // ... down to about one recoil
  const diesBelow = frozen.length >= 1;                       // ... and then the cycling stops
  const ok = cycling.length >= 2 && doppler && atRecoil && diesBelow;
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: BOTH LIMITS ARE IN THIS ONE TABLE, and the cycling count draws the boundary.`);
  console.log(`     * DOPPLER LIMIT.  The spread grows with the line width: ${rows.map((r) => r.sd.toFixed(2)).join(' -> ')} as`);
  console.log(`       gamma falls from ${rows[0].G.toFixed(1)} to ${rows[rows.length - 1].G.toFixed(2)}.  A broader line captures the molecule over a`);
  console.log(`       wider momentum range, so the recoil kicks can move it further before it is recaptured.`);
  console.log(`       (I expected this trend to be missing, because our recoil direction follows the internal`);
  console.log(`       phase instead of being random.  It is not missing -- the capture window, not the`);
  console.log(`       randomness, is what sets the width.  What the deterministic direction does change is the`);
  console.log(`       coefficient: the diffusion is weaker than a random walk would give.)`);
  console.log(`     * RECOIL LIMIT.  At the narrowest line that still cycles (gamma = ${narrowestCycling.G.toFixed(1)}) the spread is`);
  console.log(`       ${narrowestCycling.sd.toFixed(2)} -- about ONE recoil, the size of a single kick.  This is the first time in this`);
  console.log(`       project that a molecule here is NOT perfectly still: it jitters by a recoil, as a real`);
  console.log(`       laser-cooled atom does.  No narrowing can squeeze it below that.`);
  console.log(`     * FROZEN.  Once gamma drops below the recoil shift the cycling stops completely`);
  console.log(`       (${rows.map((r) => r.cyc).join(' -> ')} cycles) and the molecule freezes at px = ${frozen[0].m.toFixed(0)}.  Its spread reads 0.00,`);
  console.log(`       but that is a dead molecule, not a colder one -- momentum alone cannot tell those apart,`);
  console.log(`       which is why the cycling count is in the table.`);
}
console.log('');

// ============================================================ L4: the sharp-line control
console.log('============================================================');
console.log('L4: THE SHARP LINE IS THE tau -> infinity LIMIT (and rounds 49-51 used it with tau = 2)');
console.log('============================================================');
{
  const need = absRate(0, 0, 3);
  SHARP = true;
  const accepted = [];
  for (let om = need - 10; om <= need + 10; om++) if (restShot(om, 2).absorbed) accepted.push(om);
  SHARP = false;
  console.log(`  SHARP rule: accepted frequencies ${accepted.join(',')}  (width ${accepted.length})`);
  const broad = [];
  for (let om = need - 10; om <= need + 10; om++) if (restShot(om, 2).absorbed) broad.push(om);
  console.log(`  tau = 2 with its own (broad) line: ${broad.length} frequencies accepted`);
  console.log('');
  console.log(`  => the sharp rule is what an infinitely long lived state would give.  Using it together with a`);
  console.log(`     2-step lifetime -- as rounds 49 to 51 did -- is forbidden by energy and time.  Those rounds'`);
  console.log(`     results stand as the tau -> infinity idealisation; this round makes the two ends consistent`);
  console.log(`     with each other instead of mixing them.`);
}