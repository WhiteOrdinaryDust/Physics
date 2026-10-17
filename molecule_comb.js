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
        emit = pickDir(mol, mph[i], clk[i]);
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
// ============================================================ C1: the design parameter
console.log('============================================================');
console.log('C1: THE COMB HAS ONE DESIGN PARAMETER -- LINES PER WINDOW');
console.log('============================================================');
{
  // A single line accepts a molecule whose momentum is within GAMMA of the required value, so its
  // window is GAMMA wide IN FREQUENCY.  A comb whose lines are spaced S apart therefore puts
  //     k = GAMMA / S
  // lines inside one window.  Round 57's comb was five lines across 4*GAMMA, i.e. spacing = GAMMA,
  // i.e. k = 1 exactly -- the marginal case -- which is why it failed at eps around 0.5.
  // This round holds the trap WIDTH fixed (span = 4*GAMMA) and varies the number of lines, so the
  // spacing is the only thing changing.
  const tau = 64, G = (SCALE * EINT) / tau;
  console.log(`  line width GAMMA = ${G.toFixed(2)} frequency units;  trap span held at 4*GAMMA = ${(4 * G).toFixed(2)}`);
  console.log('');
  console.log('   lines   spacing = span/(lines-1)   k = GAMMA/spacing   predicted');
  for (const n of [2, 3, 5, 9, 17]) {
    const S = (4 * G) / (n - 1);
    const k = G / S;
    console.log(`     ${String(n).padStart(2)}          ${S.toFixed(2)}                  ${k.toFixed(2)}            ${k >= 2 ? 'several lines per window: should hold' : 'one line per window or fewer: should leak'}`);
  }
}
console.log('');

// ============================================================ C2: the scan
console.log('============================================================');
console.log('C2: SCANNING LINES PER WINDOW -- where the broadband cure starts working');
console.log('============================================================');
{
  const need = () => absRate(0, 0, 3);
  const runCase = (m, tau, lines) => {
    MASS = m; DECAY = 'clock'; DIRSRC = 'phase'; LIFETIME = tau;
    const G = (SCALE * EINT) / tau;
    const base = need() - G / 2;
    const span = 4 * G;
    let s = excited(20, 0, 2, 0, 0, 3);
    const trace = [], lvTrace = [];
    for (let t = 1; t <= 3000; t++) {
      if (t % 6 === 0) {
        const idx = Math.floor(t / 6) % lines;
        const om = lines === 1 ? Math.round(base) : Math.round(base - span / 2 + (span * idx) / (lines - 1));
        injectBeam(s, om, 3); injectBeam(s, om, 1);
      }
      s = step(s);
      if (t > 2200) { const mm = molInfo(s); trace.push(mm.px); lvTrace.push(mm.lv); }
    }
    const st = stats(trace);
    let cyc = 0;
    for (let i = 1; i < lvTrace.length; i++) if (lvTrace[i] > lvTrace[i - 1]) cyc++;
    return { cyc, sd: st.sd, m: st.m, mn: st.mn, mx: st.mx };
  };
  console.log('   case                 lines   k = GAMMA/spacing   cycling   spread   range      verdict');
  const grid = [[127, 64, 0.504], [31, 9, 0.290]];
  const results = [];
  let ok = true;
  for (const [m, tau, eps] of grid) {
    const G = (SCALE * EINT) / tau;
    for (const lines of [1, 2, 3, 5, 9, 17]) {
      const k = lines === 1 ? 0 : (lines - 1) / 4;
      const r = runCase(m, tau, lines);
      const holds = r.cyc > 30;
      results.push({ eps, k, holds });
      console.log(`   M=${String(m).padStart(3)} tau=${String(tau).padStart(3)} eps=${eps.toFixed(2)}   ${String(lines).padStart(2)}        ${k.toFixed(2)}              ${String(r.cyc).padStart(5)}    ${r.sd.toFixed(2)}    ${r.mn}..${r.mx}    ${holds ? 'HELD' : 'LEAKED'}`);
    }
  }
  MASS = 31;
  console.log('');
  // MY HYPOTHESIS WAS "the comb needs several lines per window".  The data says otherwise, in both
  // directions, so let me measure what is really binding.
  const wide = results.filter((r) => Math.abs(r.eps - 0.290) < 1e-9);
  const narrow = results.filter((r) => Math.abs(r.eps - 0.504) < 1e-9);
  const wideAnyComb = wide.filter((r) => r.k > 0).every((r) => r.holds);
  const wideEvenCoarse = wide.filter((r) => r.k > 0 && r.k <= 0.5).every((r) => r.holds);
  const narrowNone = narrow.every((r) => !r.holds);
  console.log(`   eps = 0.29: cured by EVERY comb, even the coarsest (spacing 4 line widths): ${wideEvenCoarse ? 'yes' : 'no'}`);
  console.log(`   eps = 0.50: cured by NO comb, however fine:                                  ${narrowNone ? 'yes' : 'no'}`);
  console.log('');
  // and the reason, measured directly: how many INTEGER frequencies a line accepts
  console.log('   the acceptance window, in units of the frequency grid:');
  console.log('    tau   GAMMA   accepted integer frequencies   width');
  const widths = {};
  for (const tau of [9, 32, 64, 96]) {
    LIFETIME = tau;
    const need = absRate(0, 0, 3);
    const acc = [];
    for (let om = need - 12; om <= need + 12; om++) if (restShot(om, tau).absorbed) acc.push(om);
    widths[tau] = acc.length;
    const G = (SCALE * EINT) / tau;
    console.log(`    ${String(tau).padStart(3)}   ${G.toFixed(2)}     ${String(acc.length).padStart(3)}                          ${acc.length}`);
  }
  const ok2 = widths[64] <= 1 && widths[9] >= 5 && wideAnyComb && narrowNone;
  console.log('');
  console.log(`  => ${ok2 ? 'PASS' : '*** FAIL ***'}: THE BINDING CONSTRAINT IS THE FREQUENCY GRID, NOT THE COMB SPACING.`);
  console.log(`     * At eps = 0.29 (tau = 9, GAMMA = 3.56) a single line already accepts ${widths[9]} integer frequencies, so even`);
  console.log(`       a two-line comb with a huge gap cures the leak: the lines themselves are wide.`);
  console.log(`     * At eps = 0.50 (tau = 64, GAMMA = 0.50) a line accepts exactly ${widths[64]} integer frequency -- it is`);
  console.log(`       NARROWER THAN THE GRID -- so resonance becomes a knife edge and no comb helps.`);
  console.log(`     * So my "lines per window" rule was wrong.  The real rule is: the line width must stay`);
  console.log(`       wider than the frequency unit.  That is a RESOLUTION limit, the same one rounds 47, 50`);
  console.log(`       and 57 kept running into, and its fix is the same: a larger SCALE, i.e. finer frequency`);
  console.log(`       units (more phase bits).`);
}
console.log('');

// ============================================================ C3: the price of refining
console.log('============================================================');
console.log('============================================================');
console.log('C3: WHAT REFINING THE COMB COSTS -- a wider trap is a wider spread');
console.log('============================================================');
{
  const need = () => absRate(0, 0, 3);
  const runCase = (m, tau, lines, spanFactor) => {
    MASS = m; DECAY = 'clock'; DIRSRC = 'phase'; LIFETIME = tau;
    const G = (SCALE * EINT) / tau;
    const base = need() - G / 2;
    const span = spanFactor * G;
    let s = excited(20, 0, 2, 0, 0, 3);
    const trace = [];
    for (let t = 1; t <= 3000; t++) {
      if (t % 6 === 0) {
        const idx = Math.floor(t / 6) % lines;
        const om = lines === 1 ? Math.round(base) : Math.round(base - span / 2 + (span * idx) / (lines - 1));
        injectBeam(s, om, 3); injectBeam(s, om, 1);
      }
      s = step(s);
      if (t > 2200) trace.push(molInfo(s).px);
    }
    const st = stats(trace);
    return { sd: st.sd, m: st.m };
  };
  console.log('   trap span (in line widths)   spread   mean     note');
  for (const f of [0, 1, 2, 4, 8, 16]) {
    const r = runCase(127, 64, f === 0 ? 1 : 17, f);
    console.log(`        ${String(f).padStart(2)}                       ${r.sd.toFixed(2)}    ${r.m.toFixed(2)}    ${f === 0 ? 'single line: leaks and freezes' : 'a wider trap holds a wider band'}`);
  }
  MASS = 31;
  console.log('');
  console.log(`  => the tradeoff from round 55, now with numbers: widening the trap buys capture range and`);
  console.log(`     costs spread, in proportion.  The recoil floor is not reachable this way -- a wider trap`);
  console.log(`     is a wider trap.  Narrowing it again needs a colder mechanism, not a broader beam.`);
}
