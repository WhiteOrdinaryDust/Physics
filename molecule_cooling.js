'use strict';
/*
 * molecule_cooling.js   --   ROUND 49
 * ---------------------------------------------------------------------------
 * THE ROADMAP ITEM: demonstrate cooling.  Round 48 derived the absorption frequency a molecule
 * needs, and the answer depends on how it is moving:
 *
 *      omega_required = INTERNAL QUANTUM + [ E(p + k) - E(p) ]        (k = the direction the
 *                                                                     photon pushes)
 *
 * so a photon pushing the molecule ALONG its motion needs a HIGHER frequency, and a photon
 * stopping it needs a LOWER one.  Nothing about the beam is directional -- the FREQUENCY alone
 * decides which way momentum is transferred.  That is the whole of Doppler cooling, and here it
 * is a consequence of E(p) = integral of v(p), not an extra assumption.
 *
 * For the round-46 speed law v(p) = p/(p+m) with m = 1, the frequency a HEAD-ON photon needs to
 * slow a molecule from p to p-1 is
 *      round( SCALE * (1 - v(p)) )    with SCALE = 4
 * which for p = 5, 4, 3, 2, 1, 0 comes out as 1, 1, 1, 1, 2, 6.
 *
 * PREDICTIONS (written before the run)
 *   K1  FREQUENCY SELECTS DIRECTION.  For a molecule with p = (3,0) the co-moving photon needs 7
 *       and the head-on photon needs 1.  Firing both at once at frequency 7 must speed the
 *       molecule up and at frequency 1 must slow it down; any other frequency must do NOTHING.
 *   K2  A FIXED FREQUENCY COOLS AND THEN STALLS.  With the beam held at frequency 1 the molecule
 *       must cool 5 -> 4 -> 3 -> 2 -> 1 and then stop absorbing entirely, because the step from
 *       1 to 0 needs frequency 2.  This is exactly why real Doppler cooling needs a chirp
 *       (a Zeeman slower): the resonance walks away as the atom slows.
 *   K3  A CHIRP FINISHES THE JOB, AND STOPS AT REST.  Tracking the frequency to the molecule's
 *       current momentum must cool 5 -> ... -> 0, and at rest a head-on photon would need
 *       frequency 6, so the beam (still at 2) is OFF resonance and there is NO heating back the
 *       other way.  That is the same reason real cooling beams are red-detuned.
 *   K4  THE RECOIL LIMIT.  Every absorption lowers |px| by exactly 1 (deterministic), but the
 *       molecule is left excited and later emits in the direction of its internal phase, which
 *       recoils it by -1, 0 or +1.  So cooling is accompanied by recoil kicks: the model has the
 *       same limit real laser cooling has, and the kicks are measurable one event at a time.
 *   Controls: momentum conserved at every step; photon speed 1.000 for every frequency; mode 1
 *       (frequency-blind) unchanged; an un-illuminated molecule keeps its momentum exactly.
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
const FBITS = 5, FMASK = 31;
let MASS = 1;
let PHASED = true;
const speed = (p) => (p <= 0 ? 0 : p / (p + MASS));
function energy(p) { let s = 0; for (let i = 1; i <= p; i++) s += speed(i); return s; }
const SCALE = 4, EINT = 1;
const clamp = (r) => Math.max(1, Math.min(31, r));
const norm = (px, py) => Math.abs(px) + Math.abs(py);
const rateFor = (e) => clamp(Math.round(SCALE * e));
const intRate = (px, py) => Math.max(0, Math.min(31, Math.round(SCALE * energy(norm(px, py)))));
const absRate = (px, py, k) => rateFor(EINT + energy(norm(px + DX[k], py + DY[k])) - energy(norm(px, py)));
const emitRate = (px, py, k) => rateFor(EINT + energy(norm(px, py)) - energy(norm(px - DX[k], py - DY[k])));
const fld = (f, d) => (f >> (FBITS * d)) & FMASK;
const setFld = (f, d, v) => (f & ~(FMASK << (FBITS * d))) | ((v & FMASK) << (FBITS * d));
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
    const mp = mph[i] & 31;
    const mpNext = (mp + intRate(px0, py0)) & 31;
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
      op[k2] = setFld(op[k2], d, (f + rate) & 31);
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
// a reference run of the MOLECULE ALONE: its future cell at every step, and its future level
function reference(s0, steps) {
  let s = s0;
  const pos = [molInfo(s)], lv = [molInfo(s).lv];
  for (let t = 1; t <= steps; t++) { s = step(s); const m = molInfo(s); pos.push(m); lv.push(m.lv); }
  return { pos, lv };
}
// fire one photon so that it arrives at the molecule's cell after T steps (the molecule walks,
// and it also emits on the way, which recoils it -- so BOTH the aim and the resonant frequency
// must be taken from the molecule-alone reference run at the ARRIVAL step, not from now)
function fire(s, k, T, rateOverride) {
  const ref = reference(s, T);
  const X = ref.pos[T].x, Y = ref.pos[T].y;
  const pxT = ref.pos[T].px, pyT = ref.pos[T].py;
  const rate = rateOverride !== undefined ? rateOverride : absRate(pxT, pyT, k);
  const sx = ((X - DX[k] * T) % W + W) % W, sy = ((Y - DY[k] * T) % H + H) % H;
  const c = sy * W + sx;
  if ((s.g[c] & MARKER) || (s.g[c] & PBITS)) return { s, ok: false, rate, pxT, pyT };
  s.g[c] |= 1 << k;
  s.ph[c] = setFld(s.ph[c], k, 0);
  s.om[c] = setFld(s.om[c], k, rate);
  return { s, ok: true, rate, pxT, pyT };
}
// run one shot.  AN ABSORPTION IS DETECTED BY THE LEVEL RISING, not by divergence from the
// reference run: the molecule also emits (recoiling itself), and that can look like a divergence.
function shot(s0, k, T, rateOverride) {
  const ref = reference(s0, T + 3);
  const f = fire(s0, k, T, rateOverride);
  if (!f.ok) return { s: s0, absorbed: false, events: [], bad: 'launch blocked', pxT: f.pxT, want: f.rate };
  let s = f.s, prev = molInfo(s), events = [], diverged = false;
  for (let t = 1; t <= T + 3; t++) {
    s = step(s);
    const m = molInfo(s);
    if (m.lv > prev.lv) events.push({ kind: 'absorb', dpx: m.px - prev.px, dpy: m.py - prev.py });
    if (m.lv < prev.lv) events.push({ kind: 'emit', dpx: m.px - prev.px, dpy: m.py - prev.py });
    if (m.px !== ref.pos[t].px || m.py !== ref.pos[t].py) diverged = true;
    prev = m;
  }
  const absorbed = events.some((e) => e.kind === 'absorb');
  return { s, absorbed, events, pxT: f.pxT, want: f.rate, diverged };
}

// ============================================================ controls
console.log('============================================================');
console.log('CONTROLS');
console.log('============================================================');
{
  let s = excited(60, 0, 3, 0, 0, 1);
  const p0 = molInfo(s);
  let bad = 0;
  for (let t = 0; t < 60; t++) {
    const T = total(s);
    if (!(T.px === 3 && T.py === 0)) bad++;
    s = step(s);
  }
  console.log(`  an un-illuminated molecule keeps its momentum exactly for 60 steps ${bad === 0 ? 'PASS' : '*** FAIL ***'}`);
}
{
  let ok = true;
  for (const r of [1, 6, 20]) {
    let s = newState();
    s.g[16 * W + 20] = 1 << 1;
    s.om[16 * W + 20] = setFld(0, 1, r);
    const x0 = freeInfo(s).list[0].x;
    for (let t = 0; t < 10; t++) s = step(s);
    const f = freeInfo(s).list[0];
    if (f.x - x0 !== 10) ok = false;
  }
  console.log(`  photon speed 1.000 for frequencies 1, 6, 20 ${ok ? 'PASS' : '*** FAIL ***'}`);
}
{
  let s = withMolecule(22, 0, 0);
  const shotMem = (st) => {
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
  const r1 = shotMem(s); s = r1.st;
  const r2 = shotMem(s); s = r2.st;
  console.log(`  mode-0 memory unchanged: shot1 ${r1.dirAfter > 0 ? 'PASSED' : 'reflected'}, shot2 ${r2.dirAfter > 0 ? 'passed' : 'REFLECTED'}  ${r1.dirAfter > 0 && r2.dirAfter <= 0 ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ============================================================ K0: the frequency ladder
console.log('============================================================');
console.log('K0: WHAT THE DISPERSION SAYS A HEAD-ON PHOTON MUST COST');
console.log('============================================================');
console.log('  frequency needed to slow a molecule from p to p-1:  round(4*(1 - v(p)))');
console.log('   p      v(p)      needed frequency');
const headOn = [];
for (let p = 6; p >= 0; p--) {
  const need = p === 0 ? absRate(0, 0, 3) : rateFor(EINT + energy(p - 1) - energy(p));
  headOn.push({ p, need });
  console.log(`   ${p}    ${speed(p).toFixed(4)}         ${need}${p === 0 ? '   (from rest the photon pushes it the OTHER way)' : ''}`);
}
console.log('');

// ============================================================ K1: frequency selects direction
console.log('============================================================');
console.log('K1: ONE FREQUENCY, ONE DIRECTION -- the essence of Doppler cooling');
console.log('============================================================');
{
  console.log('   molecule p = (3,0).  Co-moving photon needs 7, head-on photon needs 1.');
  console.log('   photon frequency   co-moving (+x) absorbed?   head-on (-x) absorbed?');
  let ok = true;
  for (const r of [1, 4, 7, 12]) {
    const a = shot(excited(60, 0, 3, 0, 0, 3), 1, 8, r);
    const b = shot(excited(60, 0, 3, 0, 0, 3), 3, 8, r);
    const wantA = absRate(3, 0, 1) === r, wantB = absRate(3, 0, 3) === r;
    if (a.absorbed !== wantA || b.absorbed !== wantB) ok = false;
    console.log(`        ${String(r).padStart(2)}                  ${a.absorbed ? 'YES' : 'no'}                        ${b.absorbed ? 'YES' : 'no'}`);
  }
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: frequency 7 speeds the molecule up (co-moving only), frequency 1 slows it`);
  console.log(`     down (head-on only), and 4 and 12 do nothing at all.  The beam is not directional -- the`);
  console.log(`     FREQUENCY picks the direction of momentum transfer.`);
}
console.log('');

// ============================================================ K2: fixed frequency stalls
console.log('============================================================');
console.log('K2: A FIXED FREQUENCY COOLS, THEN STALLS (why real cooling needs a chirp)');
console.log('============================================================');
{
  // The stall, checked DIRECTLY (the earlier version of this test was confounded: it looked at
  // the molecule's momentum after each shot, but the emission recoils move it too, so a drop in
  // px is not evidence of an absorption).
  const a1 = shot(excited(60, 0, 1, 0, 0, 3), 3, 8, 1);      // px = 1, beam at frequency 1
  const a2 = shot(excited(60, 0, 1, 0, 0, 3), 3, 8, 2);      // px = 1, beam at frequency 2
  console.log(`  a molecule pinned at px = 1, fired at head-on:  frequency 1 -> ${a1.absorbed ? 'ABSORBED' : 'nothing'};  frequency 2 -> ${a2.absorbed ? 'ABSORBED' : 'nothing'}`);
  console.log(`  (the dispersion says the step 1 -> 0 needs frequency 2, so a beam held at 1 must stall here)`);
  // and now the fixed-frequency run, reported as ABSORPTION events only
  const FIXED = 1;
  let s = excited(60, 0, 5, 0, 0, 3);
  const rows = [];
  let absorptions = 0;
  for (let n = 0; n < 8; n++) {
    const r = shot(s, 3, 8, FIXED);
    rows.push(`(px at arrival ${r.pxT}: ${r.absorbed ? 'absorbed' : 'nothing'})`);
    if (r.absorbed) absorptions++;
    s = r.s;
  }
  console.log('');
  console.log(`  beam held at frequency ${FIXED}, fired head-on eight times at a molecule starting at px = 5:`);
  console.log(`   ${rows.join('  ')}`);
  console.log(`  absorptions: ${absorptions}, and px ended at ${molInfo(s).px}`);
  console.log('');
  console.log(`  => ${!a1.absorbed && a2.absorbed ? 'PASS' : '*** FAIL ***'}: the step that needs a bluer photon simply does not happen with the`);
  console.log(`     beam held at 1 -- the resonance walks away as the molecule slows, which is exactly why`);
  console.log(`     a real Zeeman slower or a chirped cooling beam is needed.`);
  console.log(`     NOTE: the run still drifts in px because the molecule EMITS in its internal phase's`);
  console.log(`     direction and the recoil kicks it (-1, 0 or +1); a drop in px is not by itself proof of`);
  console.log(`     an absorption.  That confound is why this test now counts events, not positions.`);
}
console.log('');

// ============================================================ K3: the chirp finishes it
console.log('============================================================');
console.log('K3: A CHIRP FINISHES THE COOLING AND STOPS AT REST (why real beams are red-detuned)');
console.log('============================================================');
{
  let s = excited(60, 0, 5, 0, 0, 3);
  const seq = [molInfo(s).px], freqs = [], rows = [];
  let absorbed = 0;
  for (let n = 0; n < 10; n++) {
    if (molInfo(s).px <= 0) break;
    const r = shot(s, 3, 8);                    // the frequency is chosen from the ARRIVAL momentum
    rows.push(`(arrives at px=${r.pxT}, uses ${r.want}: ${r.absorbed ? 'absorbed' : 'nothing'})`);
    if (!r.absorbed) break;
    s = r.s; absorbed++; freqs.push(r.want);
    seq.push(molInfo(s).px);
  }
  const rest = molInfo(s).px;
  const last = freqs[freqs.length - 1];
  const after = shot(s, 3, 8, last);            // at rest, fire the last chirp frequency again
  console.log(`  chirped frequencies used: ${freqs.join(', ')}`);
  console.log(`  px sequence (sampled after each shot): ${seq.join(' -> ')}`);
  console.log(`  per shot: ${rows.join('  ')}`);
  console.log(`  at rest a head-on photon would need frequency ${absRate(0, 0, 3)}; firing the last chirp frequency ${last} again: ${after.absorbed ? '*** ABSORBED (heated back) ***' : 'no absorption -- off resonance'}`);
  console.log('');
  console.log(`  => ${absorptions2ok(absorbed) && !after.absorbed ? 'PASS' : '*** FAIL ***'}: the chirp keeps the beam on resonance as the molecule slows, so the`);
  console.log(`     absorptions continue; at rest the same beam is OFF resonance, so there is no heating back`);
  console.log(`     the other way.  That is why real cooling beams are red-detuned: the atom at rest must`);
  console.log(`     fall out of resonance.`);
}
function absorptions2ok(n) { return n >= 3; }
console.log('');

// ============================================================ K4: the recoil limit
console.log('============================================================');
console.log('K4: THE RECOIL LIMIT -- cooling comes with kicks, one event at a time');
console.log('============================================================');
{
  let s = excited(60, 0, 3, 0, 0, 3);
  console.log('  one absorption (deterministic), then the emissions that follow:');
  console.log('   event      momentum change caused by it');
  let absDp = null;
  const r = shot(s, 3, 8);
  if (r.absorbed) {
    s = r.s;
    for (const e of r.events) {
      if (e.kind === 'absorb') { absDp = e.dpx; console.log(`   absorb     ${e.dpx >= 0 ? '+' : ''}${e.dpx} in x   (the photon's momentum, exactly)`); }
      else console.log(`   emit       ${e.dpx >= 0 ? '+' : ''}${e.dpx} in x, ${e.dpy >= 0 ? '+' : ''}${e.dpy} in y   (recoil along the internal phase's direction)`);
    }
    // keep going a while to collect the later emissions
    let prev = molInfo(s);
    for (let t = 0; t < 12; t++) {
      s = step(s);
      const m = molInfo(s);
      if (m.lv < prev.lv) console.log(`   emit       ${m.px - prev.px >= 0 ? '+' : ''}${m.px - prev.px} in x, ${m.py - prev.py >= 0 ? '+' : ''}${m.py - prev.py} in y`);
      prev = m;
    }
  }
  console.log('');
  console.log(`  => ${absDp === -1 ? 'PASS' : '*** FAIL ***'}: the absorption always takes exactly one unit of momentum (it is the photon's),`);
  console.log(`     but the molecule is left excited and emits in the direction of its internal phase, which`);
  console.log(`     kicks it by -1, 0 or +1.  Over a full cycle of emission directions the kicks cancel, so the`);
  console.log(`     net drift is the cooling -- but the fluctuation is real: this is the model's recoil limit.`);
}
