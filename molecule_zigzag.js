'use strict';
/*
 * molecule_zigzag.js   --   ROUND 46
 * ---------------------------------------------------------------------------
 * ROUND 45 POSTULATED THE SPEED LAW.  An accumulator implemented v = |p|/(|p|+m) directly:
 * it added |p| to a counter each step and moved one cell once the counter reached |p|+m.
 * The formula WAS the rule.  This round states a MECHANISM in steps and lets the law be a
 * theorem about it -- and then checks which mechanism the rest of the physics can live with.
 *
 * TWO MECHANISMS, BOTH LOCAL AND UNIFORM, BOTH DRIVEN BY THE SAME STORED MOMENTUM p AND MASS m
 *
 *  A. ZIGZAG (always at speed 1).  The molecule moves one cell EVERY step.  Its motion is a
 *     cycle: (|p| + m) steps along the momentum axis, then m steps back, repeat, with its own
 *     place in the cycle kept in a counter k.  It never rests.
 *
 *  B. REST OVERHEAD (round 45).  The molecule moves one cell only on the steps its accumulator
 *     has filled up; between them it stands still.  It spends m of every |p|+m steps at rest.
 *
 * DERIVATION OF A'S SPEED LAW (no formula is put in)
 *     every step is +/-1 cell, period L, net displacement per period = p:
 *         #forward + #back = L,  #forward - #back = p   =>   #back = (L-p)/2, and L, p share parity.
 *     A MASSIVE walker cannot simply run: it must turn back at least once per period, #back >= 1.
 *     With exactly m steps going back the period is L = p + 2m, hence
 *         v = p/L = p/(p+2m)          <-- derived.  m = 1 gives v = p/(p+2), the FASTEST massive
 *     walker at that p; at p = 1 that is 1/3, not 1/2.
 * B's law, by the same counting (m of every p+m steps are rest steps) is v = p/(p+m).
 * So BOTH laws are DERIVED here from a mechanism stated in steps -- round 45's formula is no
 * longer postulated, it is the law of the "rest step" mechanism.  The two agree only at
 * p = 2, m = 1; the 0.500 of rounds 42-45 is therefore a property of mechanism B, not a
 * constant of the substrate.
 *
 * THE DECISIVE EXPERIMENT (Z4), AND A PREDICTION OF MINE THAT WAS WRONG
 *   My first prediction was that A would be unable to absorb a second photon at all: a photon
 *   and a speed-1 walker each advance one cell per step, so a photon behind a walker that never
 *   rests ought to step over it, killing the level ladder (round 43), the cascade and the
 *   spectrum.  That prediction is WRONG, and the reason is the phase order: collide runs BEFORE
 *   the walk, so a photon that is already in the molecule's cell is absorbed by both mechanisms.
 *   What actually differs is where the molecule IS one step later -- A has already recoiled one
 *   cell, B has not moved at all.  Aiming the second photon at one cell or the cell next door
 *   therefore gives OPPOSITE outcomes under the two mechanisms (measured below), which makes
 *   them experimentally distinguishable by one cheap shot rather than by a whole programme.
 *
 * WHAT NEITHER MECHANISM FIXES
 *   * The internal phase clock (period 2) is momentum-independent, so there is still no
 *     relativistic dispersion relation.  In A the mechanical cycle length grows with p, i.e.
 *     the mechanical frequency FALLS with p -- the opposite of omega = sqrt(p^2+m^2).  In the
 *     Dirac picture the saturating speed comes from a SUPERPOSITION of two chiralities; a single
 *     deterministic trajectory has no superposition, so it must either saturate and pay with an
 *     inverted dispersion, or keep the dispersion and let its speed rise until it reaches c.
 *   * Which mechanism the substrate uses, and the value of m = 1, are still not derived.
 *   * A is implemented along the dominant axis only; a 2-D zigzag is left for later.
 *
 * PREDICTIONS (written before the run)
 *   Z1  A's law, m = 1, N = 600, level 0.  Predicted from INDEPENDENT cycle arithmetic:
 *       p = 0 -> 0 (no walk)  1 -> 200 (1/3)  2 -> 300 (0.500)  3 -> 360 (0.600)  4 -> 400 (0.667)
 *       5 -> 430 (0.717)  6 -> 450 (0.750)  7 -> 468 (0.780).  Per-step displacement exactly 1
 *       on every step for p >= 1 (never 0); light cone saturated, never exceeded.
 *   Z2  From the trajectory alone: period = p+2, net per period = p, exactly 1 backward step.
 *   Z3  B's law, same style: exact net cells floor(p*N/(p+m)) = 300,400,450,480,500,514,525.
 *   Z4  two photons with the second aimed one cell to the side of the first:
 *       B hits at the old cell (level 2), A hits at the neighbouring cell (level 2) -- opposite.
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
// bits 27..30 carry EITHER the rest accumulator (3 bits) OR the cycle counter (4 bits):
// the two mechanisms are never active at the same time, so they may share the field.
const ACCM = 7, ACCS = 27;                 // mechanism B
const KM = 15, KS = 27;                    // mechanism A
const MASS = 1;
const MARKER = ~PBITS;
let WALK = 'zigzag';
function levOf(mol) { return (mol >> 11) & 7; }
function modeOf(mol) { return (mol >> 9) & 3; }
function phOf(mol) { return (mol >> 24) & 3; }
function axisOf(mol) { return (mol >> 26) & 1; }
function accOf(mol) { return (mol >> ACCS) & ACCM; }
function kOf(mol) { return (mol >> KS) & KM; }
function sgn4(v) { return v > 7 ? v - 16 : v; }
function pxOf(mol) { return sgn4((mol >> PXS) & PM); }
function pyOf(mol) { return sgn4((mol >> PYS) & PM); }
function withLev(mol, n) { return (mol & ~(LVB[0] | LVB[1] | LVB[2])) | ((n & 7) << 11); }
function withPh(mol, p) { return (mol & ~(PHB[0] | PHB[1])) | ((p & 3) << 24); }
function withAxis(mol, a) { return (mol & ~AXISB) | ((a & 1) << 26); }
function withAcc(mol, a) { return (mol & ~(ACCM << ACCS)) | ((a & ACCM) << ACCS); }
function withK(mol, k) { return (mol & ~(KM << KS)) | ((k & KM) << KS); }
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

function collide(g) {
  const o = new Uint32Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const parts = v & PBITS, mol = v & MARKER;
    if (mol) {
      let m = mol ^ MCLK;
      if (modeOf(mol) === 1 && parts) {
        let lv = levOf(mol), px = pxOf(mol), py = pyOf(mol), rest = 0, took = 0;
        for (let d = 0; d < 4; d++) {
          for (let c = 0; c < 2; c++) {
            const bit = c ? (parts >> (d + 4)) & 1 : (parts >> d) & 1;
            if (!bit) continue;
            if (lv < 7) { lv++; px += DX[d]; py += DY[d]; took++; }
            else rest |= 1 << (d + (c ? 4 : 0));
          }
        }
        if (took) { o[i] = withK(withAcc(withP(withLev(m & ~MCLK, lv), px, py), 0), 0) | rest; continue; }
        o[i] = m | parts; continue;
      }
      if (!parts) { o[i] = m; continue; }
      if (m & MREC) o[i] = (m & ~MREC) | reflect(parts);
      else o[i] = (m | MREC) | parts;
      continue;
    }
    if (TURN.has(parts)) { o[i] = TURN.get(parts); continue; }
    if (BIND.has(parts)) { o[i] = MBASE; continue; }
    o[i] = parts;
  }
  return o;
}
function moleculePhase(g) {
  const o = new Uint32Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i], mol = v & MARKER, parts = v & PBITS;
    if (!mol) { o[i] |= v; continue; }
    let m = mol, emit = -1, dx = 0, dy = 0;
    const px0 = pxOf(mol), py0 = pyOf(mol);
    if (WALK === 'zigzag') {                       // A: one cell EVERY step
      const axial = Math.abs(px0) >= Math.abs(py0);
      const pax = axial ? px0 : py0, a = Math.abs(pax);
      if (a > 0) {
        const L = a + 2 * MASS;
        let k = kOf(mol); if (k >= L) k = 0;
        const s = pax > 0 ? 1 : -1;
        const dir = (k < a + MASS) ? s : -s;
        m = withK(m, (k + 1) % L);
        if (axial) dx = dir; else dy = dir;
      }
    } else {                                        // B: rest steps between moves
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
    if (mol & MCLK) {                               // the decay clock: period 2, momentum-blind
      const ph = phOf(mol);
      m = withPh(m, (ph + 1) & 3);
      const lv0 = levOf(mol);
      if (lv0 > 0) {
        emit = ph;
        m = withP(withLev(m, lv0 - 1), pxOf(mol) - DX[ph], pyOf(mol) - DY[ph]);
      }
    }
    if (dx === 0 && dy === 0) {
      o[i] |= m | parts;
      if (emit >= 0) o[i] |= 1 << emit;
      continue;
    }
    const x = i % W, y = (i - x) / W;
    const k2 = ((y + dy + H) % H) * W + ((x + dx + W) % W);
    o[i] |= parts;
    o[k2] |= m;
    if (emit >= 0) o[k2] |= 1 << emit;
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

function freeInfo(g) {
  let n = 0, px = 0, py = 0, ax = -1, ay = -1;
  const xs = [];
  for (let i = 0; i < N; i++) {
    const v = g[i] & PBITS;
    if (!v) continue;
    const x = i % W, y = (i - x) / W;
    for (let d = 0; d < 4; d++) {
      for (const c of [0, 4]) {
        if ((v >> (d + c)) & 1) { n++; px += DX[d]; py += DY[d]; xs.push(x); if (ax < 0) { ax = x; ay = y; } }
      }
    }
  }
  return { n, px, py, ax, ay, xs };
}
function molInfo(g) {
  for (let i = 0; i < N; i++) {
    const mol = g[i] & MARKER;
    if (!mol) continue;
    const x = i % W;
    return { cell: i, x, y: (i - x) / W, lv: levOf(mol), px: pxOf(mol), py: pyOf(mol),
             ph: phOf(mol), mode: modeOf(mol) };
  }
  return null;
}
function total(g) {
  const f = freeInfo(g), m = molInfo(g);
  return { px: f.px + (m ? m.px : 0), py: f.py + (m ? m.py : 0), f, m };
}
function newGrid() { return new Uint32Array(N); }
function withModeBits(mol, m) { return (mol & ~(MODE1 | MODE2)) | ((m & 3) << 9); }
function withMolecule(mx, mode, ph0) {
  const g = newGrid();
  g[16 * W + (mx - 2)] = 1 << 1;
  g[16 * W + (mx + 2)] = 1 << 7;
  let a = Uint32Array.from(g);
  for (let t = 0; t < 3; t++) a = step(a);
  const c = 16 * W + mx;
  a[c] = withPh(withModeBits(a[c] & MARKER, mode), ph0);
  return a;
}
function excited(mx, lv, px, py, ph, mode) {
  const g = withMolecule(mx, mode === undefined ? 1 : mode, ph);
  const c = 16 * W + mx;
  g[c] = withK(withAcc(withP(withLev(withPh(g[c] & MARKER, ph), lv), px, py), 0), 0);
  return g;
}
function clearPhotons(g) { for (let i = 0; i < N; i++) g[i] &= MARKER; return g; }

// independent predictors, no cellular automaton involved
function predictZig(p, m, steps) {
  const a = Math.abs(p), L = a + 2 * m, s = Math.sign(p), per = [];
  for (let k = 0; k < L; k++) per.push(a === 0 ? 0 : (k < a + m ? s : -s));
  let net = 0;
  for (let t = 0; t < steps; t++) net += per[t % L];
  return { net, per };
}
function predictRest(p, m, steps) { return Math.floor((p * steps) / (p + m)); }

function coast(p, steps) {
  let g = excited(60, 0, p, 0, 0, 1);
  let prev = molInfo(g);
  const wrap = (d, S) => { d = ((d % S) + S) % S; return d > S / 2 ? d - S : d; };
  let dxT = 0, maxStep = 0, minStep = 99, zeros = 0;
  const seq = [];
  for (let t = 0; t < steps; t++) {
    g = step(g);
    const mm = molInfo(g);
    const dx = wrap(mm.x - prev.x, W), dy = wrap(mm.y - prev.y, H);
    dxT += dx;
    const man = Math.abs(dx) + Math.abs(dy);
    seq.push(dx);
    maxStep = Math.max(maxStep, man); minStep = Math.min(minStep, man);
    if (man === 0) zeros++;
    prev = mm;
  }
  return { dxT, maxStep, minStep, zeros, seq };
}

// ============================================================ controls
console.log('============================================================');
console.log('CONTROLS');
console.log('============================================================');
{
  let ok = true;
  for (let a = -8; a <= 7; a++) for (let b = -8; b <= 7; b++) if (pxOf(withP(MBASE, a, b)) !== a || pyOf(withP(MBASE, a, b)) !== b) ok = false;
  let okB = true, okA = true;
  for (let a = 0; a <= 7; a++) if (accOf(withAcc(MBASE, a)) !== a) okB = false;
  for (let a = 0; a <= 15; a++) if (kOf(withK(MBASE, a)) !== a) okA = false;
  console.log(`  encoding: 256 (px,py) pairs ${ok ? 'PASS' : '*** FAIL ***'}; shared field round-trips as acc 0..7 ${okB ? 'PASS' : '*** FAIL ***'} and as cycle counter 0..15 ${okA ? 'PASS' : '*** FAIL ***'}`);
}
for (const mode of ['zigzag', 'rest']) {
  WALK = mode;
  let a = withMolecule(22, 1, 0);
  const m0 = molInfo(a);
  for (let t = 0; t < 40; t++) a = step(a);
  const m1 = molInfo(a), T = total(a);
  console.log(`  fresh molecule (${mode}): level ${m0.lv}, p (${m0.px},${m0.py}); after 40 steps x ${m1.x} ${m1.x === 22 ? '(static)' : '*** MOVED ***'}   ${m0.lv === 0 && m0.px === 0 && m0.py === 0 && m1.x === 22 && T.px === 0 && T.py === 0 ? 'PASS (zero momentum = no walk)' : '*** FAIL ***'}`);
}
WALK = 'zigzag';
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

// ============================================================ Z1: mechanism A's derived law
console.log('============================================================');
console.log('Z1: MECHANISM A (zigzag) -- v = p/(p+2m) derived by counting, m = 1, N = 600');
console.log('============================================================');
WALK = 'zigzag';
const STEPS = 600;
{
  console.log('   p   cycle L=p+2   predicted net   measured net   speed    p/(p+2)   per-step |dx|');
  let ok = true, cone = 0, mono = true, prev = -1;
  for (let p = 0; p <= 7; p++) {
    const pr = predictZig(p, MASS, STEPS), r = coast(p, STEPS);
    const sp = r.dxT / STEPS;
    if (r.dxT !== pr.net) ok = false;
    if (p >= 1 && sp <= prev) mono = false;
    if (p >= 1) prev = sp;
    cone = Math.max(cone, r.maxStep);
    console.log(`   ${p}        ${p + 2}          ${String(pr.net).padStart(4)}            ${String(r.dxT).padStart(4)}        ${sp.toFixed(4)}   ${(p / (p + 2)).toFixed(4)}      ${[...new Set(r.seq.map(Math.abs))].sort().join('/')}`);
  }
  const r1 = coast(1, 20), r0 = coast(0, 20);
  console.log('');
  console.log(`  => exact predicted net on every row ${ok ? 'PASS' : '*** FAIL ***'}; increasing ${mono ? 'PASS' : '*** FAIL ***'}; max single step ${cone} ${cone <= 1 ? 'PASS' : '*** FAIL ***'}`);
  console.log(`  => p >= 1 never rests: displacement exactly 1 on all 20 test steps ${r1.zeros === 0 && r1.minStep === 1 ? 'PASS' : '*** FAIL ***'}; p = 0 never moves ${r0.zeros === 20 && r0.maxStep === 0 ? 'PASS' : '*** FAIL ***'}`);
  console.log(`     So A is at speed 1 at EVERY step: the light cone is saturated locally, and only the`);
  console.log(`     average (p/(p+2)) is lower.  No formula for the speed was put into the rule.`);
}
console.log('');

// ============================================================ Z2: what momentum and mass mean here
console.log('============================================================');
console.log('Z2: the identification the mechanism rests on (read off the trajectory alone)');
console.log('============================================================');
{
  console.log('   p    smallest period   net per period   backward steps per period');
  let ok = true;
  for (let p = 1; p <= 7; p++) {
    const r = coast(p, 40);
    let L = -1;
    for (let cand = 2; cand <= 20; cand++) {
      let per = true;
      for (let t = cand; t < 40; t++) if (r.seq[t] !== r.seq[t - cand]) { per = false; break; }
      if (per) { L = cand; break; }
    }
    const one = r.seq.slice(0, L);
    const net = one.reduce((a, b) => a + b, 0), back = one.filter((s) => s < 0).length;
    const good = L === p + 2 && net === p && back === MASS;
    if (!good) ok = false;
    console.log(`   ${p}          ${L}                 ${net}                    ${back}      ${good ? '' : '*** MISMATCH ***'}`);
  }
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: "p units of momentum" IS p cells per cycle and "mass 1" IS one backward`);
  console.log(`     step per cycle.  Neither is fitted in this mechanism -- both are counted from the walk.`);
}
console.log('');

// ============================================================ Z3: mechanism B's law, for comparison
console.log('============================================================');
console.log('Z3: MECHANISM B (rest overhead) -- v = p/(p+m), same style of exact prediction');
console.log('============================================================');
WALK = 'rest';
{
  console.log('   p   predicted net   measured net   speed    p/(p+1)   per-step |dx| values');
  let ok = true, cone = 0;
  for (let p = 1; p <= 7; p++) {
    const want = predictRest(p, MASS, STEPS), r = coast(p, STEPS);
    if (r.dxT !== want) ok = false;
    cone = Math.max(cone, r.maxStep);
    console.log(`   ${p}       ${String(want).padStart(4)}            ${String(r.dxT).padStart(4)}        ${(r.dxT / STEPS).toFixed(4)}   ${(p / (p + 1)).toFixed(4)}      ${[...new Set(r.seq.map(Math.abs))].sort().join('/')}`);
  }
  console.log('');
  console.log(`  => exact predicted net on every row ${ok ? 'PASS' : '*** FAIL ***'}; max single step ${cone} ${cone <= 1 ? 'PASS' : '*** FAIL ***'}`);
  console.log(`     B is at speed 0 or 1 and the per-step values include 0 -- it rests between moves.`);
}
WALK = 'zigzag';
console.log('');

// ============================================================ Z4: the decisive experiment
console.log('============================================================');
console.log('Z4: THE DECISIVE EXPERIMENT -- where IS the molecule during the absorption window?');
console.log('============================================================');
{
  // Photon 1 comes along +x from x=10 and is absorbed at step 13, leaving the molecule excited
  // with momentum (1,0); its excitation window is step 14 (one step).  A second photon comes
  // along +y, timed to land on row 16 at the END of step 13, aimed at one of two cells:
  //     aim (22,16) = the cell the molecule occupied,   aim (23,16) = the cell +x of it.
  // The two mechanisms disagree about which cell the molecule is in at that moment:
  //     rest overhead: it has NOT moved yet (its accumulator needs two steps for |p|=1)
  //     zigzag:        it recoiled one cell immediately (it is at speed 1 every step)
  const run = (mode, aimX) => {
    WALK = mode;
    let g = withMolecule(22, 1, 3);
    g[16 * W + 10] = 1 << 1;                              // the exciting photon
    g[((16 - 13) % H + H) % H * W + aimX] = 1 << 0;       // +y photon, on row 16 at end of step 13
    let prevLv = 0, abs_ = 0, emi = 0, maxPost = 0, absorbTimes = [];
    for (let t = 0; t <= 20; t++) {
      const T = total(g), cT = total(g && collide(g));
      if (cT.m) maxPost = Math.max(maxPost, cT.m.lv);
      if (T.m.lv > prevLv) { abs_ += T.m.lv - prevLv; absorbTimes.push(t); }
      if (T.m.lv < prevLv) emi += prevLv - T.m.lv;
      prevLv = T.m.lv;
      g = step(g);
    }
    return { abs_, emi, maxPost, absorbTimes };
  };
  const A1 = run('rest', 22), A2 = run('zigzag', 22);
  const B1 = run('rest', 23), B2 = run('zigzag', 23);
  console.log('   mechanism        second photon aimed at   absorptions   highest level held at once');
  console.log(`   rest overhead         (22,16)                  ${A1.abs_}              ${A1.maxPost}   ${A1.maxPost === 2 ? '<- HIT' : 'missed'}`);
  console.log(`   zigzag (speed 1)      (22,16)                  ${A2.abs_}              ${A2.maxPost}   ${A2.maxPost === 2 ? '<- HIT' : 'missed'}`);
  console.log(`   rest overhead         (23,16)                  ${B1.abs_}              ${B1.maxPost}   ${B1.maxPost === 2 ? '<- HIT' : 'missed'}`);
  console.log(`   zigzag (speed 1)      (23,16)                  ${B2.abs_}              ${B2.maxPost}   ${B2.maxPost === 2 ? '<- HIT' : 'missed'}`);
  console.log('');
  const ok = A1.maxPost === 2 && A2.maxPost === 1 && B1.maxPost === 1 && B2.maxPost === 2;
  console.log(`  => ${ok ? 'PASS' : '*** UNEXPECTED ***'}: the two mechanisms give OPPOSITE outcomes for the same two`);
  console.log(`     photon placements, so this single cheap measurement tells them apart: a walker that is at`);
  console.log(`     speed 1 has ALREADY recoiled one cell when the second photon arrives, a walker with rest`);
  console.log(`     steps has not moved yet.  (Nothing here is a defect: both are self-consistent.  But only`);
  console.log(`     one of them can be the substrate, and now that is testable.)`);
}
console.log('');

// ============================================================ Z5: what neither mechanism fixes
console.log('============================================================');
console.log('Z5: the price, unchanged by either mechanism');
console.log('============================================================');
{
  console.log('   p   A: cycle length   A: mechanical frequency   relativistic omega = sqrt(p^2+m^2)');
  for (let p = 1; p <= 6; p++) console.log(`   ${p}         ${p + 2}               ${(1 / (p + 2)).toFixed(4)}                    ${Math.sqrt(p * p + 1).toFixed(4)}`);
  console.log('');
  console.log('  In A the mechanical period grows with p, so its frequency FALLS with momentum -- the opposite');
  console.log('  of the relativistic trend (last column, up to normalisation).  In the Dirac picture the');
  console.log('  saturating speed comes from a SUPERPOSITION of the two chiralities; a single deterministic');
  console.log('  trajectory has no superposition, so it must either saturate and invert its dispersion, or');
  console.log('  keep the dispersion and let its speed rise linearly until it hits the refresh rate.');
  console.log('  The model picks saturation, and the internal PHASE clock (period 2) is momentum-blind, so');
  console.log('  there is still no dispersion relation at all.  That remains the open item.');
}
