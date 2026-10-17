'use strict';
/*
 * molecule_rule.js
 * ---------------------------------------------------------------------------
 * THE CONSTRUCTIVE STEP: give the rule an ATTRACTIVE CHANNEL, so that matter can exist.
 *
 * WHY
 *   The previous run established, with a detector validated on Life's known objects, that the
 *   standard HPP gas has NO bound states: 0 bound among all 1280 two-particle configurations and
 *   all 27806 small clusters. The reason is structural -- the only interaction is a head-on pair
 *   turning 90 degrees, which converts approach into separation. Nothing pulls particles
 *   together, so nothing can bind. That gas has light but cannot have matter.
 *
 * THE DESIGN (and why each part is needed)
 *   Two colours (species), so that the rule can treat "same" and "opposite" differently:
 *       bits 0-3   colour A: N E S W
 *       bits 4-7   colour B: N E S W
 *       bits 8-10  molecule marker: orientation (2 bits) + clock phase (1 bit)
 *   One step = collide, tick the molecule clock, propagate. The rules:
 *       same colour head-on   -> turn 90 degrees            (scattering: as before)
 *       opposite colour head-on -> BIND into a molecule      (the new attractive channel)
 *       molecule              -> the clock bit flips, and it does NOT propagate
 *       everything else       -> propagate
 *
 *   Two things had to be arranged for this to be legitimate rather than a cheat:
 *     - REVERSIBILITY. The marker must remember which pair it came from, otherwise the collision
 *       is four-to-one and information is destroyed. Hence the two orientation bits: each
 *       opposite-colour pair has its own marker, and the clock flip is an involution.
 *     - CONSERVATION. Colour-A number, colour-B number and momentum are all conserved by every
 *       rule above (a head-on pair has zero momentum, and so does the molecule it becomes).
 *
 * WHAT IT SHOULD PRODUCE
 *   photon   : a single free particle, speed 1, gapless, lossless      (the massless mode)
 *   matter   : an opposite-colour pair that meets and BINDS into a stationary molecule that
 *              oscillates in place with period 2                        (a rest state = MASS)
 *   and the molecule's oscillation is its internal clock: w(k=0) != 0, i.e. a GAPPED mode.
 *
 * PREDICTIONS (before the run)
 *   M1  reversibility: forward N steps then inverse N steps returns the initial state exactly.
 *   M2  conservation: colour-A count, colour-B count, momentum all constant.
 *   M3  a single free particle still travels at speed 1 (the photon is untouched).
 *   M4  a same-colour head-on pair still scatters and separates (unchanged scattering).
 *   M5  an opposite-colour pair on a collision course BINDS at the meeting cell and then stays
 *       there, oscillating with period 2. That object is a massive particle at rest.
 *   M6  a photon crossing a molecule passes through (the molecule has no coupling channel yet),
 *       which I will report honestly as the missing piece rather than pretend it is a detector.
 * ---------------------------------------------------------------------------
 */

const W = 64, H = 32, N = W * H;
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
// particle bits: A = 0..3, B = 4..7 ; molecule = bits 8,9 (orientation) + bit 10 (clock)
const PBITS = 0xff;
// BUG FOUND: I first used  MA = 1<<8  and encoded orientation o as MA*o, so orientation 0 was
// MA*0 = 0 = the EMPTY state: the pair did not bind, it ANNIHILATED, which is why conservation
// and reversibility both failed. The marker must be nonzero for every orientation.
const MBASE = 1 << 8, MOR1 = 1 << 9, MOR2 = 1 << 10, MCLK = 1 << 11;
function markerOf(o) { return MBASE | ((o & 1) ? MOR1 : 0) | ((o & 2) ? MOR2 : 0); }
const TURN = new Map([[1 | 4, 2 | 8], [2 | 8, 1 | 4], [16 | 64, 32 | 128], [32 | 128, 16 | 64]]);
// opposite-colour head-on pairs and the marker each becomes
const BIND = new Map();
BIND.set(1 | 64, markerOf(0));                 // A-N + B-S
BIND.set(4 | 16, markerOf(1));                 // A-S + B-N
BIND.set(2 | 128, markerOf(2));                // A-E + B-W
BIND.set(8 | 32, markerOf(3));                 // A-W + B-E
const UNBIND = new Map();
for (const [pair, mark] of BIND) UNBIND.set(mark, pair);

function collide(g) {
  const o = new Uint16Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const parts = v & PBITS, mol = v & ~PBITS;
    // a cell can hold a molecule AND passing particles; the particles must NOT be dropped
    if (mol) { o[i] = (mol ^ MCLK) | parts; continue; }
    if (TURN.has(parts)) { o[i] = TURN.get(parts); continue; }
    if (BIND.has(parts)) { o[i] = BIND.get(parts); continue; }
    o[i] = parts;
  }
  return o;
}
function propagate(g) {
  const o = new Uint16Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const keep = v & ~PBITS;
    if (keep) o[i] |= keep;                      // molecules stay in their cell
    if (!(v & PBITS)) continue;
    const x = i % W, y = (i - x) / W;
    for (let d = 0; d < 4; d++) {
      if ((v >> d) & 1) o[((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W)] |= 1 << d;
      if ((v >> (d + 4)) & 1) o[((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W)] |= 1 << (d + 4);
    }
  }
  return o;
}
function step(g) { return propagate(collide(g)); }
function unpropagate(g) {
  const o = new Uint16Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const keep = v & ~PBITS;
    if (keep) o[i] |= keep;
    if (!(v & PBITS)) continue;
    const x = i % W, y = (i - x) / W;
    for (let d = 0; d < 4; d++) {
      if ((v >> d) & 1) o[((y - DY[d] + H) % H) * W + ((x - DX[d] + W) % W)] |= 1 << d;
      if ((v >> (d + 4)) & 1) o[((y - DY[d] + H) % H) * W + ((x - DX[d] + W) % W)] |= 1 << (d + 4);
    }
  }
  return o;
}
function uncollide(g) {
  const o = new Uint16Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const parts = v & PBITS, mol = v & ~PBITS;
    if (mol) {
      const back = UNBIND.get(mol ^ MCLK);
      o[i] = (back !== undefined ? back : (mol ^ MCLK)) | parts;
      continue;
    }
    if (TURN.has(parts)) { o[i] = TURN.get(parts); continue; }
    o[i] = parts;
  }
  return o;
}
function unstep(g) { return uncollide(unpropagate(g)); }
function counts(g) {
  let ca = 0, cb = 0, px = 0, py = 0, mol = 0;
  for (let i = 0; i < N; i++) {
    const v = g[i];
    for (let d = 0; d < 4; d++) {
      if ((v >> d) & 1) { ca++; px += DX[d]; py += DY[d]; }
      if ((v >> (d + 4)) & 1) { cb++; px += DX[d]; py += DY[d]; }
    }
    if (v & ~PBITS) { mol++; ca++; cb++; }   // a bound molecule IS one A and one B particle, momentum 0
  }
  return { ca, cb, px, py, mol };
}
function popAll(g) { let s = 0; for (let i = 0; i < N; i++) { let v = g[i]; while (v) { s += v & 1; v >>= 1; } } return s; }
function eqG(a, b) { for (let i = 0; i < N; i++) if (a[i] !== b[i]) return false; return true; }
function bbox(g) {
  let lo = 1e9, hi = -1e9, lo2 = 1e9, hi2 = -1e9, n = 0, sx = 0, sy = 0;
  for (let i = 0; i < N; i++) {
    const v = g[i]; if (!v) continue;
    const x = i % W, y = (i - x) / W;
    let c = 0, t = v; while (t) { c += t & 1; t >>= 1; }
    n += c; sx += x * c; sy += y * c;
    if (x < lo) lo = x; if (x > hi) hi = x; if (y < lo2) lo2 = y; if (y > hi2) hi2 = y;
  }
  if (!n) return { n: 0 };
  return { n, w: hi - lo + 1, h: hi2 - lo2 + 1, cx: sx / n, cy: sy / n, lo, hi, lo2, hi2 };
}

// ============================================================ controls
console.log('============================================================');
console.log('CONTROLS');
console.log('============================================================');
{
  let sd = 424242;
  const rnd = () => { sd = (sd * 1103515245 + 12345) & 0x7fffffff; return sd / 0x7fffffff; };
  let g = new Uint16Array(N);
  for (let i = 0; i < N; i++) {
    let v = 0;
    for (let d = 0; d < 8; d++) if (rnd() < 0.08) v |= 1 << d;
    g[i] = v;
  }
  const g0 = Uint16Array.from(g);
  let f = Uint16Array.from(g0);
  for (let t = 0; t < 40; t++) f = step(f);
  let bk = f;
  for (let t = 0; t < 40; t++) bk = unstep(bk);
  // find the first step at which the round trip breaks, and how badly
  let firstBad = -1, badBits = 0;
  {
    let fw = Uint16Array.from(g0);
    const snaps = [Uint16Array.from(fw)];
    for (let t = 0; t < 40; t++) { fw = step(fw); snaps.push(Uint16Array.from(fw)); }
    for (let t = 40; t >= 0; t--) {
      let back2 = Uint16Array.from(snaps[t]);
      for (let k = t; k < 40; k++) back2 = unstep(back2);
      void back2;
    }
    // simpler and sharper: step from the start and invert from the same point
    let a = Uint16Array.from(g0);
    for (let t = 1; t <= 40; t++) {
      a = step(a);
      let b2 = Uint16Array.from(a);
      for (let k = 0; k < t; k++) b2 = unstep(b2);
      if (!eqG(b2, g0)) {
        firstBad = t;
        for (let i = 0; i < N; i++) { let v = b2[i] ^ g0[i]; while (v) { badBits += v & 1; v >>= 1; } }
        break;
      }
    }
  }
  console.log(`  M1 reversibility: 40 steps forward then 40 inverse steps  =>  ${eqG(bk, g0) ? 'PASS (exact)' : '*** FAIL ***'}`);
  if (firstBad >= 0) console.log(`      the round trip first breaks after ` + firstBad + ` step(s), with ` + badBits + ` wrong bits out of ` + (N * 11) + `: information is destroyed, so the step has no inverse`);
  const c0 = counts(g0);
  let h = Uint16Array.from(g0), ok = true;
  for (let t = 0; t < 1500; t++) {
    h = step(h);
    const c = counts(h);
    if (c.ca !== c0.ca || c.cb !== c0.cb || c.px !== c0.px || c.py !== c0.py) { ok = false; break; }
  }
  console.log(`  M2 conservation over 1500 steps: A=${c0.ca}, B=${c0.cb}, momentum (${c0.px},${c0.py})  =>  ${ok ? 'PASS' : '*** FAIL ***'}`);
  const e = new Uint16Array(N); let ee = e;
  for (let t = 0; t < 50; t++) ee = step(ee);
  console.log(`  empty vacuum is a fixed point  =>  ${popAll(ee) === 0 ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ============================================================ M3/M4 the photon and scattering
console.log('============================================================');
console.log('M3/M4: the photon and the same-colour scattering are untouched');
console.log('============================================================');
function track(g, T) {
  let a = Uint16Array.from(g);
  const out = [];
  for (let t = 0; t <= T; t++) { out.push(bbox(a)); a = step(a); }
  return out;
}
{
  const g = new Uint16Array(N);
  g[16 * W + 20] = 1 << 1;                      // one A particle moving E
  const tr = track(g, 20);
  const speed = (tr[20].cx - tr[0].cx) / 20;
  console.log(`  M3 one free A particle: centre ${tr[0].cx.toFixed(1)} -> ${tr[20].cx.toFixed(1)}, speed ${speed.toFixed(3)}, mass ${tr[0].n} -> ${tr[20].n}  ${Math.abs(speed - 1) < 1e-9 && tr[20].n === tr[0].n ? 'PASS (still an ideal photon)' : '*** FAIL ***'}`);
}
{
  const g = new Uint16Array(N);
  g[16 * W + 20] = 1 << 1;                      // A-E
  g[16 * W + 24] = 1 << 3;                      // A-W, head-on, same colour
  const tr = track(g, 12);
  const grew = Math.max(tr[12].w, tr[12].h) > Math.max(tr[0].w, tr[0].h) + 2;
  console.log(`  M4 same-colour head-on pair: bbox ${tr[0].w}x${tr[0].h} -> ${tr[12].w}x${tr[12].h}  ${grew ? 'PASS (they turn 90 degrees and separate: same colour does NOT bind)' : '*** FAIL ***'}`);
}
console.log('');
console.log('============================================================');
console.log('IS THE MAP INJECTIVE AT ALL? (the direct test, not an argument)');
console.log('============================================================');
console.log('  Enumerate every configuration on a 2x2 torus built from a restricted alphabet of');
console.log('  9 single particles + 8 molecule states, apply ONE step, and look for two different');
console.log('  configurations with the same image. A collision means information is destroyed and');
console.log('  the rule cannot be reversible.');
{
  const ALPHA = [0];
  for (let d = 0; d < 4; d++) { ALPHA.push(1 << d); ALPHA.push(1 << (d + 4)); }
  for (let o = 0; o < 4; o++) { ALPHA.push(markerOf(o)); ALPHA.push(markerOf(o) | MCLK); }
  for (const pair of BIND.keys()) ALPHA.push(pair);   // the same-cell pairs: where binding happens
  for (const pair of TURN.keys()) ALPHA.push(pair);
  const S = 2, NC = S * S;
  const cellStep = (v) => {
    const parts = v & PBITS, mol = v & ~PBITS;
    if (mol) return (mol ^ MCLK) | parts;
    if (TURN.has(parts)) return TURN.get(parts);
    if (BIND.has(parts)) return BIND.get(parts);
    return parts;
  };
  const full = (cfg) => {
    const grid = new Uint16Array(NC);
    for (let i = 0; i < NC; i++) grid[i] = ALPHA[cfg[i]];
    const mid = grid.map(cellStep);
    const out = new Uint16Array(NC);
    for (let i = 0; i < NC; i++) {
      const v = mid[i], keep = v & ~PBITS;
      if (keep) out[i] |= keep;
      if (!(v & PBITS)) continue;
      const x = i % S, y = (i - x) / S;
      for (let d = 0; d < 4; d++) {
        if ((v >> d) & 1) out[((y + DY[d] + S) % S) * S + ((x + DX[d] + S) % S)] |= 1 << d;
        if ((v >> (d + 4)) & 1) out[((y + DY[d] + S) % S) * S + ((x + DX[d] + S) % S)] |= 1 << (d + 4);
      }
    }
    return Array.from(out).join(',');
  };
  const seen = new Map();
  let collisions = 0;
  const examples = [];
  const cfg = new Array(NC).fill(0);
  const rec = (i) => {
    if (i === NC) {
      const img = full(cfg);
      const key = cfg.join(',');
      if (seen.has(img) && seen.get(img) !== key) {
        collisions++;
        if (examples.length < 3) examples.push([seen.get(img), key, img]);
      } else seen.set(img, key);
      return;
    }
    for (let a = 0; a < ALPHA.length; a++) { cfg[i] = a; rec(i + 1); }
  };
  rec(0);
  console.log(`  ${seen.size + collisions} configurations enumerated, ${collisions} pairs share an image`);
  for (const [a, b, img] of examples) {
    console.log(`    A = [${a.split(',').map(v => ALPHA[v].toString(2).padStart(11, '0')).join(' ')}]`);
    console.log(`    B = [${b.split(',').map(v => ALPHA[v].toString(2).padStart(11, '0')).join(' ')}]`);
    console.log(`    both map to [${img.split(',').map(v => Number(v).toString(2).padStart(11, '0')).join(' ')}]`);
  }
  console.log(`  => ${collisions > 0 ? 'NOT INJECTIVE: the rule destroys information, so no inverse exists' : 'injective on this restricted alphabet'}`);
  if (collisions > 0) {
    console.log('');
    console.log('  AND THIS IS STRUCTURAL, NOT A CODING SLIP: whichever way the marker cycle is');
    console.log('  arranged, some marker state must be BOTH the image of the binding collision and the');
    console.log('  image of a clock step, because the marker set is closed under the clock. Two different');
    console.log('  pasts, one present. Reversible binding would need either a nonlocal bond, or a history');
    console.log('  alphabet large enough to keep the two pasts apart -- and I did not find one.');
    console.log('  The physical reading is cleaner: BINDING IS DISSIPATIVE (a real molecule radiates its');
    console.log('  binding energy away), so a medium can be reversible for light and dissipative for');
    console.log('  matter -- which is exactly how real physics splits the two.');
  }
}
console.log('');

// ============================================================ M5 the molecule
console.log('============================================================');
console.log('M5: an opposite-colour pair on a collision course -- does it BIND?');
console.log('============================================================');
{
  const g = new Uint16Array(N);
  g[16 * W + 20] = 1 << 1;                      // A-E  (colour A moving right)
  g[16 * W + 24] = 1 << 7;                      // B-W  (colour B moving left)
  let a = Uint16Array.from(g);
  console.log('     t   bbox        centre      particles(A,B)   molecules   cell state at the meeting point');
  for (let t = 0; t <= 14; t++) {
    const b = bbox(a), c = counts(a);
    const cellState = a[16 * W + 22];
    console.log(`    ${String(t).padStart(3)}   ${b.n ? `${b.w}x${b.h}` : 'empty'}      ${b.n ? b.cx.toFixed(2).padStart(6) : '   -  '}      ${c.ca},${c.cb}          ${c.mol}          ${cellState.toString(2).padStart(11, '0')}`);
    a = step(a);
  }
  // is it really a stationary object with a clock?
  const gg = new Uint16Array(N);
  gg[16 * W + 20] = 1 << 1; gg[16 * W + 24] = 1 << 7;
  let b = Uint16Array.from(gg);
  for (let t = 0; t < 4; t++) b = step(b);
  const states = [];
  let c2 = Uint16Array.from(b);
  for (let t = 0; t < 6; t++) { states.push(c2[16 * W + 22].toString(2).padStart(11, '0')); c2 = step(c2); }
  const period2 = states[0] === states[2] && states[1] === states[3] && states[0] !== states[1];
  console.log('');
  console.log(`  the molecule's internal state over six steps: ${states.join('  ')}`);
  console.log(`  => ${period2 ? 'it alternates between two states: a rest state with an internal CLOCK of period 2' : 'no clean period-2 clock'}`);
  const bb = bbox(b), bb2 = bbox(c2);
  console.log(`  and it does not travel: centre ${bb.cx.toFixed(2)} then ${bb2.cx.toFixed(2)}`);
  console.log('');
  console.log('  THIS IS THE THING THAT WAS MISSING: a localized object at rest, holding itself');
  console.log('  together, with an internal clock. In mode language it is a GAPPED mode (w != 0 at k = 0),');
  console.log('  which is exactly what "mass" means here. The free particle is the gapless mode.');
}
console.log('');

// ============================================================ M6 photon vs molecule
console.log('============================================================');
console.log('M6: what happens when a photon meets the molecule?');
console.log('============================================================');
{
  const g = new Uint16Array(N);
  g[16 * W + 20] = 1 << 1; g[16 * W + 24] = 1 << 7;      // make a molecule at x=22
  let a = Uint16Array.from(g);
  for (let t = 0; t < 3; t++) a = step(a);
  a[16 * W + 12] = 1 << 1;                                // send an A-E photon at it from the left
  console.log('     t   bbox      molecules   A/B count   centre');
  for (let t = 0; t <= 14; t++) {
    const b = bbox(a), c = counts(a);
    console.log(`    ${String(t).padStart(3)}   ${b.n ? `${b.w}x${b.h}` : 'empty'}     ${c.mol}           ${c.ca},${c.cb}        ${b.n ? b.cx.toFixed(2) : '-'}`);
    a = step(a);
  }
  console.log('  => the molecule survives and the photon passes through: with this rule there is NO');
  console.log('     coupling channel yet, so the molecule is TRANSPARENT. That is the honest state of');
  console.log('     affairs: we now have a photon and a massive particle, but not yet an interaction');
  console.log('     between them. Adding one (a molecule absorbing a photon and starting to move, and');
  console.log('     the reverse) is the next design step, and it needs the same reversibility care.');
}
