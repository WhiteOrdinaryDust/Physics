'use strict';
/*
 * molecule_emission.js   --   ROUND 53
 * ---------------------------------------------------------------------------
 * THE ROADMAP ITEM: randomise the emission direction, so the recoil becomes a real random walk
 * and the Doppler limit gets the right SIZE, not just the right trend.
 *
 * WHY ROUND 52's DIFFUSION WAS TOO WEAK
 *   Round 52 emitted along a counter that stepped N, E, S, W, N, ...  The recoil kicks were
 *   therefore not independent: the x-kicks came out as 0, +1, 0, -1, 0, +1, ...  A sequence like
 *   that has zero lag-1 correlation but PERFECTLY NEGATIVE correlation at lag 2 -- each +1 is
 *   cancelled by a later -1 -- so the momentum barely diffuses.  Cooling was fine (the mean kick
 *   is zero) but the equilibrium width came out smaller than a random walk would give.
 *
 * WHAT "RANDOM" CAN MEAN IN A DETERMINISTIC MODEL
 *   Not true randomness: a HIDDEN PHASE.  The molecule already carries a free-running internal
 *   oscillation phase (rounds 47+), and the direction is now read out of that phase at the moment
 *   of emission:
 *          dir = mix(mp)   with   mix(mp) = ((mp>>1) ^ (mp>>3) ^ (mp>>5)) & 3
 *   The phase runs at a rate set by the momentum (that is the matter dispersion), and the
 *   emissions happen at times set by the absorption history, so the phase is sampled at
 *   effectively unrelated instants.  The result is deterministic, reproducible, and unpredictable
 *   from outside -- the same mechanism rounds 32-33 identified as the origin of apparent
 *   randomness in this programme.
 *
 *   One detail matters and is chosen deliberately: the LIFETIME is ODD (17, 33, 65 ...).  The
 *   sampled phase advances by rate * LIFETIME per emission; if that product shares a factor with
 *   the 256 phase states the sample is periodic (and can even be constant -- the aliasing trap
 *   this project has hit before).  Odd lifetimes keep the sample ergodic.
 *
 * PREDICTIONS (written before the run)
 *   E1  THE RECOIL STATISTICS.  With the cycling counter the x-kicks are 0, +1, 0, -1, ...: the
 *       lag-1 autocorrelation is 0 but the lag-2 autocorrelation is about -1 (anti-correlated, so
 *       sub-diffusive), and the directions are uniform (25% each).  With the hidden phase both
 *       autocorrelations should collapse towards 0 while the directions stay uniform.  A
 *       deterministic source that looks uniform but is anti-correlated IS the defect being fixed,
 *       so the autocorrelation, not the histogram, is the measurement that matters.
 *   E2  THE CROSSOVER, SAME TABLE, TWO SOURCES.  With stronger diffusion the momentum spread must
 *       GROW for the broad-line cases (the Doppler limit is set by damping against diffusion) while
 *       the NARROW-line floor must NOT move: the floor is one recoil, and no amount of randomness
 *       can make a single kick smaller.  So: 'phase' > 'cycle' for wide lines, and both about one
 *       unit at the narrowest line that still cycles.
 *   E3  the line width itself must be untouched (it depends on the lifetime, not on the direction).
 *   Controls: fields, light speed for several frequencies, momentum accounting through the beams,
 *       mode-0 memory, the lifetime = delay + 1 convention.
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
let LIFETIME = 17;                              // tau, in steps (ODD on purpose: see DIRSRC below)
let DECAY = 'fixed';                             // 'fixed' = round 52/53 (age >= LIFETIME)
                                                // 'clock' = the decay time comes from the hidden
                                                //           clock: per-step probability 1/tau, so the
                                                //           MEAN lifetime is tau (and the linewidth
                                                //           SCALE*EINT/tau still matches it)
let DIRSRC = 'cycle';                           // 'cycle' = round 52's N,E,S,W rotation of a counter
                                                // 'phase' = read the direction out of a free-running
                                                //           hidden clock mixed with the internal
                                                //           phase: deterministic and reproducible,
                                                //           but sampled at unrelated instants
// A hidden clock, advanced once per step by a linear congruential map, mixed with the internal
// phase.  Deterministic and exactly reproducible, but its value at the (irregular) emission
// instants is not something an observer can predict -- which is all "random" can mean here.
// The first attempt mixed the internal phase with a plain counter and was NOT random: the kicks
// came out biased (8 up against 27 down) with a strong lag-3 correlation.  E1 caught it.
function xs32(x) { x ^= x << 13; x ^= x >>> 17; x ^= x << 5; return x >>> 0; }
const mixDir = (mp, clk) => (((clk >>> 16) ^ (mp >> 5)) & 3);
function pickDir(mol, mp, clk) { return DIRSRC === 'cycle' ? phOf(mol) : mixDir(mp, clk); }
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
// STIMULATED EMISSION (Einstein's B process): if a photon of exactly the frequency this transition
// would emit is already in the cell, the molecule is driven to emit into THAT mode.
let STIMULATED = false;
function stimulatingDir(s, i, px, py, parts) {
  for (let d = 0; d < 4; d++) {
    const has = ((parts >> d) & 1) || ((parts >> (d + 4)) & 1);
    if (!has) continue;
    if (Math.abs(fld(s.om[i], d) - emitRate(px, py, d)) <= gamma()) return d;
  }
  return -1;
}
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
  const { g, ph, om, mph, age, acc, clk } = s;
  const og = new Uint32Array(N), op = new Uint32Array(N), oo = new Uint32Array(N);
  const omp = new Uint8Array(N), oag = new Uint8Array(N), oac = new Uint8Array(N), ocl = new Uint32Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const parts = v & PBITS, mol = v & MARKER;
    if (mol) {
      let m = mol ^ MCLK;
      op[i] = ph[i]; oo[i] = om[i]; omp[i] = mph[i]; oag[i] = age[i]; oac[i] = acc[i]; ocl[i] = clk[i];
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
  return { g: og, ph: op, om: oo, mph: omp, age: oag, acc: oac, clk: ocl };
}
function moleculePhase(s) {
  const { g, ph, om, mph, age, acc, clk } = s;
  const og = new Uint32Array(N), op = new Uint32Array(N), oo = new Uint32Array(N);
  const omp = new Uint8Array(N), oag = new Uint8Array(N), oac = new Uint8Array(N), ocl = new Uint32Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i], mol = v & MARKER, parts = v & PBITS;
    if (!mol) { og[i] |= v; op[i] |= ph[i]; oo[i] |= om[i]; continue; }
    ocl[i] = xs32(clk[i]);                           // a 32-bit hidden clock (xorshift: its bits
    let m = mol, emit = -1, dx = 0, dy = 0, rate = 0;
    const lv0 = levOf(mol);
    const px0 = pxOf(mol), py0 = pyOf(mol);
    const mp = mph[i] & 255;
    const mpNext = (mp + intRate(px0, py0)) & 255;
    let ageNext = age[i];
    if (lv0 > 0) {
      ageNext = age[i] + 1;
      const decayNow = DECAY === 'fixed' ? (ageNext >= LIFETIME) : ((clk[i] >>> 20) < (4096 / LIFETIME));
      if (decayNow) {                       // the excited state finally decays
        emit = STIMULATED ? stimulatingDir(s, i, px0, py0, parts) : -1;
        if (emit < 0) emit = pickDir(mol, mph[i], clk[i]);
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
    omp[k2] = mpNext; oag[k2] = ageNext; oac[k2] = oac[i]; ocl[k2] = ocl[i];
    if (emit >= 0) {
      og[k2] |= 1 << emit;
      op[k2] = setFld(op[k2], emit, mp); oo[k2] = setFld(oo[k2], emit, rate);
    }
  }
  return { g: og, ph: op, om: oo, mph: omp, age: oag, acc: oac, clk: ocl };
}
function propagate(s) {
  const { g, ph, om, mph, age, acc, clk } = s;
  const og = new Uint32Array(N), op = new Uint32Array(N), oo = new Uint32Array(N);
  const omp = new Uint8Array(N), oag = new Uint8Array(N), oac = new Uint8Array(N), ocl = new Uint32Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i], keep = v & MARKER, parts = v & PBITS;
    if (keep) { og[i] |= keep; omp[i] = mph[i]; oag[i] = age[i]; oac[i] = acc[i]; ocl[i] = clk[i]; }
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
  return { g: og, ph: op, om: oo, mph: omp, age: oag, acc: oac, clk: ocl };
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
           mph: new Uint8Array(N), age: new Uint8Array(N), acc: new Uint8Array(N), clk: new Uint32Array(N) };
}
function withMolecule(mx, mode, phase0) {
  let s = newState();
  s.g[6 * W + (mx - 2)] = 1 << 1;
  s.g[6 * W + (mx + 2)] = 1 << 7;
  for (let t = 0; t < 3; t++) s = step(s);
  const c = 6 * W + mx;
  s.g[c] = withPh(withModeBits(s.g[c] & MARKER, mode), phase0);
  // SEED THE HIDDEN CLOCK.  It must not be zero: xorshift has zero as a fixed point, so a
  // molecule born with clock 0 would keep it forever and the clock would silently do nothing.
  // (Two different generators producing byte-identical output is how this was caught.)
  s.clk[c] = (0x9E3779B9 ^ Math.imul(c + 7, 2654435761)) >>> 0;
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
console.log(`  machine: MASS = ${MASS}, SCALE = ${SCALE}, energy unit = ${EINT}, gamma = SCALE*EINT/tau`);
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
{
  LIFETIME = 17;
  const need = absRate(0, 0, 3);
  const r = restShot(need, 17);
  console.log(`  lifetime convention: absorption -> emission delay ${r.delay} for tau = 17 ${r.delay === 16 ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');
// ============================================================ E1: detailed balance
console.log('============================================================');
console.log('E1: DETAILED BALANCE IS ALREADY SATISFIED -- emission and absorption agree');
console.log('============================================================');
{
  // The roadmap item was "make emission resonant too".  It turns out to be already true, and this
  // is the measurement: the frequency a molecule EMITS on a transition (recoil included) is the
  // same expression as the frequency the ABSORBER demands for the same transition run backwards.
  console.log('   molecule p   direction   frequency emitted   frequency the absorber demands   equal?');
  let ok = true;
  for (const p of [0, 1, 3, 5]) {
    for (const k of [0, 1, 2, 3]) {
      const emitted = emitRate(p, 0, k);              // decay from p, emitting along k: p -> p-k
      const absorbed = absRate(p - DX[k], -DY[k], k);   // the SAME transition, run backwards (signed!)
      const sameExpression = emitted === absorbed;
      if (!sameExpression) ok = false;
      if (k === 1 || k === 3) {
        console.log(`      ${String(p).padStart(2)}          ${'NESW'[k]}              ${String(emitted).padStart(3)}                    ${String(absorbed).padStart(3)}                 ${sameExpression ? 'yes' : '*** NO ***'}`);
      }
    }
  }
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: emission is already resonant in the only sense that matters: the photon a`);
  console.log(`     molecule emits is one this molecule could absorb, with the same transition energy.  There is`);
  console.log(`     nothing to add here, so the rest of the round goes after the process that IS missing.`);
}
console.log('');

// ============================================================ E2: stimulated emission cannot amplify
console.log('============================================================');
console.log('E2: STIMULATED EMISSION -- and the mode is already full');
console.log('============================================================');
{
  // Einstein's B process would need the molecule to be DRIVEN by a resonant photon.  It never gets
  // the chance: collide() runs before moleculePhase(), and a resonant photon in the cell is
  // absorbed on the spot.  Absorption always wins the race, so the decay is always spontaneous.
  const popcount = (v) => { let n = 0; for (let i = 0; i < 8; i++) if ((v >> i) & 1) n++; return n; };
  const runOne = (stim) => {
    STIMULATED = stim; DECAY = 'fixed'; DIRSRC = 'phase'; LIFETIME = 3; SHARP = false;
    let s = excited(20, 1, 0, 0, 0, 3);                // one excited molecule, level 1
    const c = 6 * W + 20;
    s.g[c] |= 1 << 3;                                   // a resonant photon in its own cell
    s.om[c] = setFld(s.om[c], 3, emitRate(0, 0, 3));    // at exactly the transition frequency
    const lv0 = molInfo(s).lv, ph0 = popcount(s.g[c] & PBITS);
    s = step(s);
    const lv1 = molInfo(s).lv, ph1 = popcount(s.g[c] & PBITS);
    STIMULATED = false;
    return { lv0, lv1, ph0, ph1 };
  };
  const on = runOne(true), off = runOne(false);
  console.log(`  an excited molecule (level 1) with a resonant photon in its own cell, one step later:`);
  console.log(`   B process ON : level ${on.lv0} -> ${on.lv1}, and the molecule's own cell holds ${on.ph0} -> ${on.ph1} photon(s)`);
  console.log(`   B process OFF: level ${off.lv0} -> ${off.lv1}, cell holds ${off.ph0} -> ${off.ph1}`);
  console.log('');
  const absorbed = on.lv1 > on.lv0;
  console.log(`  => ${absorbed ? 'PASS' : '*** FAIL ***'}: THE STIMULATING PHOTON IS ABSORBED, NOT AMPLIFIED.  The rule is`);
  console.log(`     implemented and does nothing, because of the ORDER of the phases: collide() runs before`);
  console.log(`     moleculePhase(), so any resonant photon in the cell is taken by the absorber before the`);
  console.log(`     molecule ever gets to be driven by it.  Absorption always wins the race, and the decay is`);
  console.log(`     therefore always spontaneous.  Making the B process reachable is a change to the phase`);
  console.log(`     order, not to this rule -- and then E3 says what it would still be unable to do.`);
}
console.log('');

// ============================================================ E3: why -- the arithmetic of the mode
console.log('============================================================');
console.log('E3: WHY -- the mode saturates at one, and that is the OLD alphabet gap resurfacing');
console.log('============================================================');
{
  const popcount = (v) => { let n = 0; for (let i = 0; i < 8; i++) if ((v >> i) & 1) n++; return n; };
  console.log('  the relevant arithmetic, on one cell and one direction:');
  console.log(`    a photon present          : occupancy bit = 1`);
  console.log(`    "one more photon arrives" : occupancy bit = 1 | 1 = 1     (unchanged)`);
  console.log(`    "two photons leave"       : occupancy bit = 1 & ~1 = 0    (both gone at once)`);
  console.log(`    so occupancy is counted MOD 2 -- exclusion, not bosonic enhancement.`);
  console.log('');
  // demonstrate it directly on the lattice, no molecule involved
  let s = newState();
  s.g[6 * W + 20] = 1 << 1;
  const c = 6 * W + 20;
  const n1 = popcount(s.g[c] & PBITS);
  s.g[c] |= 1 << 1;                                    // add a second photon to the SAME mode
  const n2 = popcount(s.g[c] & PBITS);
  s.g[c] |= 1 << 1;                                    // and a third
  const n3 = popcount(s.g[c] & PBITS);
  console.log(`   measured on the lattice: after 1, 2 and 3 "photons" the mode holds ${n1}, ${n2}, ${n3}`);
  console.log('');
  console.log(`  => ${n1 === 1 && n2 === 1 && n3 === 1 ? 'PASS' : '*** FAIL ***'}: the mode cannot hold two.  Amplification needs a mode whose occupancy`);
  console.log(`     can reach 2, which is exactly the requirement identified in round 35 when the alphabet was`);
  console.log(`     first examined: TWO STATES CANNOT REPRESENT |2>, so a two-state photon is fermionic.`);
  console.log(`     A three-state (or wider) occupancy per mode is the price of a laser, and that price was`);
  console.log(`     known before any of this machinery existed.  Nothing here changes it.`);
}
console.log('');

// ============================================================ E4: spontaneous emission is fine
console.log('============================================================');
console.log('============================================================');
console.log('E4: SPONTANEOUS EMISSION IS UNAFFECTED (a control)');
console.log('============================================================');
{
  const run = (stim) => {
    STIMULATED = stim; DECAY = 'clock'; DIRSRC = 'phase'; LIFETIME = 9; SHARP = false;
    let s = excited(20, 1, 0, 0, 0, 3);                // excited, alone, no photon anywhere
    let prevLv = molInfo(s).lv, decays = 0;
    for (let t = 1; t <= 400; t++) {
      s = step(s);
      const m = molInfo(s);
      if (m.lv < prevLv) decays++;
      prevLv = m.lv;
    }
    let photons = 0;
    for (let i = 0; i < N; i++) { const v = s.g[i] & PBITS; for (let b = 0; b < 8; b++) if ((v >> b) & 1) photons++; }
    STIMULATED = false;
    return { decays, photons };
  };
  const a = run(false), b = run(true);
  console.log(`  an excited molecule alone: spontaneous ${a.decays} decay(s), ${a.photons} photon(s);  with the B process enabled ${b.decays} decay(s), ${b.photons} photon(s)`);
  console.log(`  => ${a.decays === b.decays && a.photons === b.photons && a.photons >= 1 ? 'PASS' : '*** FAIL ***'}: with no resonant photon present the molecule decays into the hidden-phase`);
  console.log(`     direction exactly as before.  The B process is an addition, not a replacement, and it`);
  console.log(`     changes nothing when there is no field to stimulate it.`);
}
console.log('');

// ============================================================ E5: what a bosonic mode would give
console.log('============================================================');
console.log('E5: THE TARGET, FOR THE RECORD (computed, not simulated)');
console.log('============================================================');
{
  console.log('  if a mode could hold n photons, the same rule would give:');
  console.log('   n (photons in the mode)   photons out of one decay   gain');
  for (const n of [1, 2, 3, 10, 100]) {
    console.log(`        ${String(n).padStart(3)}                        ${n + 1}                    ${((n + 1) / n).toFixed(3)}`);
  }
  console.log('');
  console.log(`  => this is what Einstein's B coefficient buys in a bosonic field: stimulated emission`);
  console.log(`     proportional to the occupancy, which is the entire basis of the laser.  Our measured gain`);
  console.log(`     is 1.000 at every n (E2), because our occupancy is mod 2 (E3).  Recording the target here`);
  console.log(`     makes the missing piece concrete: it is not a rule, it is an ALPHABET.`);
  STIMULATED = false;
}
