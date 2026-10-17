'use strict';
/*
 * bound_state_search.js
 * ---------------------------------------------------------------------------
 * WHERE WE ARE
 *   The last rounds produced the photon half: in the HPP lattice gas the empty vacuum carries
 *   an excitation with speed exactly 1 (the lattice refresh rate), no dispersion and no loss,
 *   and two head-on photons evolve exactly linearly until the step they meet, where they depart
 *   from linearity by exactly 4 bits, permanently. What is missing is MATTER: a localized object
 *   with rest mass, i.e. something that can sit still or move slower than light and carry
 *   internal states.
 *
 * WHAT A PARTICLE IS, OPERATIONALLY (this is the definition the detector uses)
 *   Take a configuration, follow it, and look at the bounding box of its occupied cells:
 *     BOUND_STATIC  -- the box stays small AND the shape recurs with zero drift: rest mass.
 *     BOUND_MOVING  -- the box stays small AND the shape recurs with a constant non-zero drift:
 *                      a particle with mass moving slower than the excitation speed.
 *     FREE/SCATTER  -- the box grows out to the light cone: nothing bound.
 *     ANNIHILATED   -- the mass goes to zero.
 *   A rest state is what "mass" means here: a gapped mode, w -> w0 > 0 as k -> 0.
 *
 * DETECTOR CONTROLS (a detector without them is worthless -- this project has been burned)
 *   Life is the control substrate because its objects are known:
 *     a 2x2 block  must read BOUND_STATIC
 *     a glider     must read BOUND_MOVING with period 4 and velocity (1,1)
 *     one lone cell must read ANNIHILATED
 *   then, in the lattice gas, one free particle must read FREE/SCATTER.
 *
 * PREDICTIONS (before the run)
 *   B1  the controls behave as stated.
 *   B2  IN THE STANDARD HPP GAS THERE ARE NO BOUND STATES AT ALL, for a reason I can state in
 *       advance: the only interaction is a head-on pair turning 90 degrees, and that interaction
 *       is purely REPULSIVE -- it turns approach into separation. Two particles moving in the
 *       same direction never interact; two perpendicular ones pass through; two head-on ones
 *       scatter apart. No channel anywhere brings particles together, so nothing can bind.
 *       I therefore expect ZERO bound states among all two-particle configurations and among
 *       all small multi-particle clusters.
 *   B3  if B2 holds, the constructive conclusion is that matter needs an ATTRACTIVE CHANNEL,
 *       which means the collision rule must be extended -- and the alphabet has room for it:
 *       16 states per cell in the gas, 256 if a second species is added.
 * ---------------------------------------------------------------------------
 */

// ---------------- Life (2 states), used as the detector control substrate
const LS = 64;
function lifeStep(g, S) {
  const o = new Uint8Array(S * S);
  for (let y = 0; y < S; y++) for (let x = 0; x < S; x++) {
    let n = 0;
    for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
      if (!dx && !dy) continue;
      n += g[((y + dy + S) % S) * S + ((x + dx + S) % S)];
    }
    const i = y * S + x;
    o[i] = g[i] ? ((n === 2 || n === 3) ? 1 : 0) : (n === 3 ? 1 : 0);
  }
  return o;
}
// ---------------- HPP lattice gas (16 states per cell)
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
const NB = (1 << 0) | (1 << 2), EW = (1 << 1) | (1 << 3);
function gasStep(g, W, H) {
  const N = W * H, mid = new Uint8Array(N);
  for (let i = 0; i < N; i++) { const v = g[i]; mid[i] = (v === NB) ? EW : (v === EW ? NB : v); }
  const out = new Uint8Array(N);
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
    const v = mid[y * W + x];
    if (!v) continue;
    for (let d = 0; d < 4; d++) if ((v >> d) & 1) out[((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W)] |= 1 << d;
  }
  return out;
}
// ---------------- the detector
// v1 BUG (caught by its own control): "bounding box stays small and the shape recurs with a
// constant drift" is satisfied TRIVIALLY by a single free particle, and by any set of free
// particles flying in parallel -- their shape recurs because nothing happens. So boundedness
// alone cannot mean "bound".
// The real criterion has two parts:
//   (1) the object stays spatially bounded (it does not fly apart), AND
//   (2) its parts actually INTERACT: the evolution must differ from the evolution of the parts
//       taken separately. That difference R is measured exactly, the same way the photon-photon
//       test measured it. R == 0 means "several independent free particles", not a bound state.
function popBitsOf(g) { let s = 0; for (let i = 0; i < g.length; i++) { let v = g[i]; while (v) { s += v & 1; v >>= 1; } } return s; }
function splitBits(g) {
  const out = [];
  for (let i = 0; i < g.length; i++) for (let d = 0; d < 4; d++) if ((g[i] >> d) & 1) { const b = new Uint8Array(g.length); b[i] = 1 << d; out.push(b); }
  return out;
}
function bboxOf(g, W, H) {
  let lo = 1e9, hi = -1e9, lo2 = 1e9, hi2 = -1e9, n = 0, sx = 0, sy = 0;
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
    const v = g[y * W + x];
    if (!v) continue;
    let c = 0, w = v; while (w) { c += w & 1; w >>= 1; }
    n += c; sx += x * c; sy += y * c;
    if (x < lo) lo = x; if (x > hi) hi = x; if (y < lo2) lo2 = y; if (y > hi2) hi2 = y;
  }
  return { n, sx, sy, lo, hi, lo2, hi2, w: hi - lo + 1, h: hi2 - lo2 + 1 };
}
function classify(g0, stepFn, T, W, H, boxLimit) {
  const nParts = popBitsOf(g0);
  if (nParts <= 1) {
    // a single excitation cannot be bound to anything; just say whether it survives
    let g = Uint8Array.from(g0);
    for (let t = 0; t <= T; t++) { if (popBitsOf(g) === 0) return { kind: 'ANNIHILATED', t }; g = stepFn(g); }
    return { kind: 'SINGLE EXCITATION (free by definition)', n: nParts };
  }
  const parts = splitBits(g0);
  let g = Uint8Array.from(g0);
  let ps = parts.map(p => Uint8Array.from(p));
  const frames = [];
  let maxR = 0, firstR = -1;
  for (let t = 0; t <= T; t++) {
    // R = how much the real evolution differs from "each part flies on its own"
    const pred = new Uint8Array(W * H);
    for (const p of ps) for (let i = 0; i < W * H; i++) pred[i] ^= p[i];
    let R = 0;
    for (let i = 0; i < W * H; i++) { let v = g[i] ^ pred[i]; while (v) { R += v & 1; v >>= 1; } }
    if (R > 0 && firstR < 0) firstR = t;
    if (R > maxR) maxR = R;
    const b = bboxOf(g, W, H);
    if (!b.n) return { kind: 'ANNIHILATED', t };
    if (Math.max(b.w, b.h) > boxLimit) return { kind: 'FREE/SCATTER', t, Rmax: maxR };
    let sig = '';
    for (let y = b.lo2; y <= b.hi2; y++) { for (let x = b.lo; x <= b.hi; x++) sig += (g[y * W + x] ? '1' : '0') + ','; sig += '|'; }
    frames.push({ t, sig, w: b.w, h: b.h, n: b.n, cx: b.sx / b.n, cy: b.sy / b.n });
    g = stepFn(g);
    ps = ps.map(p => stepFn(p));
  }
  // stays bounded: now ask whether it is bound or just several independent free particles
  let rec = null;
  for (let p = 1; p <= 12 && !rec; p++) {
    let ok = true, base = null;
    for (let i = 0; i + p < frames.length; i++) {
      if (frames[i].sig !== frames[i + p].sig) { ok = false; break; }
      const dx = Math.round(frames[i + p].cx - frames[i].cx), dy = Math.round(frames[i + p].cy - frames[i].cy);
      if (base === null) base = [dx, dy]; else if (dx !== base[0] || dy !== base[1]) { ok = false; break; }
    }
    if (ok && base) rec = { p, vx: base[0], vy: base[1] };
  }
  if (maxR === 0) return { kind: 'INDEPENDENT FREE PARTICLES (bounded but never interact)', n: nParts, rec };
  if (rec && (rec.vx || rec.vy)) return { kind: 'BOUND_MOVING', p: rec.p, vx: rec.vx, vy: rec.vy, n: frames[0].n, Rmax: maxR, firstR };
  if (rec) return { kind: 'BOUND_STATIC', p: rec.p, n: frames[0].n, Rmax: maxR, firstR };
  return { kind: 'BOUND (irregular shape)', n: frames[0].n, w: frames[0].w, h: frames[0].h, Rmax: maxR, firstR };
}

// ============================================================ controls, on Life
console.log('============================================================');
console.log('DETECTOR CONTROLS (on Life, where the objects are known)');
console.log('============================================================');
{
  const block = new Uint8Array(LS * LS);
  for (const [x, y] of [[30, 30], [31, 30], [30, 31], [31, 31]]) block[y * LS + x] = 1;
  const r1 = classify(block, g => lifeStep(g, LS), 30, LS, LS, 6);
  console.log(`  a 2x2 block : ${r1.kind}${r1.p ? ` (period ${r1.p})` : ''}   ${r1.kind === 'BOUND_STATIC' ? 'PASS' : '*** FAIL ***'}`);
  const gl = new Uint8Array(LS * LS);
  for (const [dx, dy] of [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]]) gl[(30 + dy) * LS + (30 + dx)] = 1;
  const r2 = classify(gl, g => lifeStep(g, LS), 30, LS, LS, 6);
  console.log(`  a glider    : ${r2.kind}${r2.p ? ` (period ${r2.p}, velocity (${r2.vx},${r2.vy}))` : ''}   ${r2.kind === 'BOUND_MOVING' && r2.p === 4 && r2.vx === 1 && r2.vy === 1 ? 'PASS' : '*** FAIL ***'}`);
  const cell = new Uint8Array(LS * LS); cell[30 * LS + 30] = 1;
  const r3 = classify(cell, g => lifeStep(g, LS), 30, LS, LS, 6);
  console.log(`  a lone cell : ${r3.kind}   ${r3.kind === 'ANNIHILATED' ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');
console.log('  and in the lattice gas, the two cases the old detector confused:');
const GW = 48, GH = 48;
{
  const g = new Uint8Array(GW * GH);
  g[24 * GW + 24] = 1 << 1;
  const r = classify(g, x => gasStep(x, GW, GH), 20, GW, GH, 6);
  console.log(`  ONE free particle             : ${r.kind}   ${r.kind.startsWith('SINGLE') ? 'PASS (an excitation, not a bound state)' : '*** FAIL ***'}`);
  const g2 = new Uint8Array(GW * GH);
  g2[24 * GW + 24] = 1 << 1; g2[24 * GW + 28] = 1 << 1;    // two particles flying in parallel
  const r2 = classify(g2, x => gasStep(x, GW, GH), 20, GW, GH, 6);
  console.log(`  TWO parallel free particles   : ${r2.kind}   ${r2.kind.startsWith('INDEPENDENT') ? 'PASS (bounded but no interaction => not bound)' : '*** FAIL ***'}`);
  const g3 = new Uint8Array(GW * GH);
  g3[24 * GW + 24] = 1 << 1; g3[24 * GW + 26] = 1 << 3;    // head-on: they do interact
  const r3 = classify(g3, x => gasStep(x, GW, GH), 20, GW, GH, 6);
  console.log(`  TWO head-on particles         : ${r3.kind} (they scatter, so the box must grow)   ${r3.kind === 'FREE/SCATTER' ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ============================================================ B2: two-particle census
console.log('============================================================');
console.log('B2a: every two-particle configuration (relative displacement x directions)');
console.log('============================================================');
{
  const T = 20, LIM = 16;
  let tested = 0, bound = 0;
  const outcomes = {};
  for (let dx = -4; dx <= 4; dx++) for (let dy = -4; dy <= 4; dy++) {
    if (dx === 0 && dy === 0) continue;
    for (let d1 = 0; d1 < 4; d1++) for (let d2 = 0; d2 < 4; d2++) {
      const g = new Uint8Array(GW * GH);
      g[(24 + dy) * GW + (24 + dx)] |= 1 << d1;
      g[24 * GW + 24] |= 1 << d2;
      const r = classify(g, x => gasStep(x, GW, GH), T, GW, GH, LIM);
      tested++;
      outcomes[r.kind] = (outcomes[r.kind] || 0) + 1;
      if (r.kind.startsWith('BOUND')) bound++;
    }
  }
  console.log(`  ${tested} configurations tested, ${bound} read as BOUND`);
  for (const k of Object.keys(outcomes).sort()) console.log(`    ${String(outcomes[k]).padStart(5)}  ${k}`);
  console.log(`  => ${bound === 0 ? 'NO BOUND STATE among any two particles (as predicted: the interaction is purely repulsive)' : '*** BOUND STATES FOUND -- prediction wrong, investigate ***'}`);
}
console.log('');

// ============================================================ B2: small clusters
console.log('============================================================');
console.log('B2b: small clusters -- all configurations with up to 3 particles in a 3x3 box,');
console.log('     plus 20000 sampled 4-particle ones');
console.log('============================================================');
{
  // 3x3 box, 9 cells x 4 directions = 36 slots
  const slots = [];
  for (let y = 0; y < 3; y++) for (let x = 0; x < 3; x++) for (let d = 0; d < 4; d++) slots.push([x, y, d]);
  const NS = slots.length;
  const T = 16, LIM = 12;
  let tested = 0, bound = 0, earlyOut = 0;
  const outcomes = {};
  const runOne = (chosen) => {
    const g = new Uint8Array(GW * GH);
    for (const si of chosen) { const [x, y, d] = slots[si]; g[(24 + y) * GW + (24 + x)] |= 1 << d; }
    const r = classify(g, x => gasStep(x, GW, GH), T, GW, GH, LIM);
    outcomes[r.kind] = (outcomes[r.kind] || 0) + 1;
    if (r.kind.startsWith('BOUND')) bound++;
    return r;
  };
  // up to 3 particles: all combinations
  for (let i = 0; i < NS; i++) { tested++; runOne([i]); }
  for (let i = 0; i < NS; i++) for (let j = i + 1; j < NS; j++) { tested++; runOne([i, j]); }
  for (let i = 0; i < NS; i++) for (let j = i + 1; j < NS; j++) for (let k = j + 1; k < NS; k++) { tested++; runOne([i, j, k]); }
  // 4 particles: sampled
  let sd = 987654321;
  const rnd = () => { sd = (sd * 1103515245 + 12345) & 0x7fffffff; return sd / 0x7fffffff; };
  let sampled = 0, tries = 0;
  while (sampled < 20000 && tries < 400000) {
    tries++;
    const s = new Set();
    while (s.size < 4) s.add(Math.floor(rnd() * NS));
    const arr = [...s];
    tested++; sampled++;
    runOne(arr);
  }
  console.log(`  ${tested} clusters tested (${sampled} of them four-particle), ${bound} read as BOUND`);
  for (const k of Object.keys(outcomes).sort()) console.log(`    ${String(outcomes[k]).padStart(6)}  ${k}`);
  console.log(`  => ${bound === 0 ? 'NO BOUND CLUSTER AT ALL: this gas cannot hold anything together' : '*** BOUND CLUSTERS FOUND ***'}`);
}
console.log('');
console.log('============================================================');
console.log('WHAT THIS MEANS');
console.log('============================================================');
console.log('  If both censuses come back with zero bound states, then the conclusion is structural and');
console.log('  not a failure of search effort: the HPP interaction is a hard-core head-on turn, which');
console.log('  converts approach into separation, so no channel anywhere pulls particles together.');
console.log('  Photons exist (speed 1, gapless, lossless, exactly linear until contact) but matter');
console.log('  cannot: there is nothing to bind with.');
console.log('');
console.log('  The constructive next step is therefore to give the rule an ATTRACTIVE CHANNEL while');
console.log('  keeping reversibility and conservation. There is room in the alphabet: a second species');
console.log('  (two colours, 8 bits per cell = 256 states) allows a rule of the form');
console.log('      same colour head-on   -> turn 90 degrees   (scattering, as now)');
console.log('      opposite colour head-on -> form a composite that advances one cell every TWO steps');
console.log('  A composite that advances one cell per two steps travels at half the excitation speed,');
console.log('  which is exactly a massive particle: it has a rest frame, a speed below the refresh rate,');
console.log('  and (being a two-step object) an internal clock to carry states.');
