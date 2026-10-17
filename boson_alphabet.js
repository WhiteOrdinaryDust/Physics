'use strict';
/*
 * boson_alphabet.js   --   ROUND 61
 * ---------------------------------------------------------------------------
 * THE OLDEST GAP IN THIS PROGRAMME, AND THE FIRST HALF OF IT.
 *
 * Round 35 identified it: a photon mode in this model is ONE OCCUPANCY BIT, so occupancy is
 * counted modulo 2, so a mode cannot hold two quanta, so the alphabet cannot represent |2>.  That
 * is why the CHSH work needed amplitudes {0,1,2} and could not get them from two states, and why
 * round 56 measured a laser gain of exactly 1.000 where a bosonic mode should give 1 + n.
 *
 * WHAT THIS ROUND CHANGES: the photon sector stores a COUNT per (cell, direction) instead of a
 * bit.  Everything else about the gas is kept the same -- one cell per step, exact head-on pairs
 * turn -- and the turn rule is extended conservatively:
 *
 *      a cell holding counts m and n in opposite directions (d and d+2)
 *          turns  min(m,n)  quanta, one into each perpendicular direction,
 *          and leaves the surplus |m-n| travelling the way the surplus was going.
 *      total number:  2*min + |m-n| = m + n                      (conserved)
 *      total momentum: the turned pair cancels, the surplus keeps its sign     (conserved)
 *
 * WHAT THIS BUYS, AND WHAT IT DOES NOT
 *   * |2> EXISTS.  Two quanta can sit in one mode, and the mode is a bosonic ladder.
 *   * The gas is now EXACTLY HOMOGENEOUS: scaling a configuration by a positive integer scales
 *     its evolution by the same integer, collisions included (min and |m-n| are homogeneous).
 *     That is the "linearity ruler" of round 35, passed at integer level.
 *   * It is still a CONE, not a vector space: counts are non-negative, so there are no negative
 *     or complex amplitudes and therefore no interference.  The gap has MOVED -- from "no |2>"
 *     to "no signs" -- which is progress, not arrival.
 *   * Stimulated emission would now give gain 1 + 1/n instead of 1.000 (round 56's target).  The
 *     molecule side is NOT rewired in this round; the gain is computed, not measured, and that is
 *     stated plainly.
 *
 * PREDICTIONS (written before the run)
 *   A  coexistence: two quanta placed in one mode survive in the count alphabet and collapse to
 *      one in the bit alphabet (which is the same arithmetic as 1|1 = 1).
 *   B  linearity: for every configuration tried, evolving the doubled state equals twice the
 *      evolution of the single state, exactly, collisions included.
 *   C  light: a pulse of n quanta travels at exactly 1 cell per step with its width unchanged, for
 *      n = 1, 2, 3, 5 -- so gapless, non-dispersive light survives the alphabet change.
 *   D  conservation: total number and total momentum are conserved over a long collision-heavy run.
 *   E  the amplification target (computed): gain 1 + 1/n for counts, 1.000 for bits.
 * ---------------------------------------------------------------------------
 */

const W = 64, H = 24, N = W * H;
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
const DIRNAME = ['N', 'E', 'S', 'W'];
let MODE = 'count';                    // 'count' = the new alphabet, 'bit' = the old one
const idx = (x, y, d) => (y * W + x) * 4 + d;

function newGas() { return new Uint8Array(4 * N); }
function clone(g) { return Uint8Array.from(g); }
// the collision: exact head-on pairs turn min(m,n) each way, the surplus continues.
// BOTH pairs in a cell must be computed from the ORIGINAL counts and written once: doing the
// second pair after the first would zero the quanta the first pair had just turned into it
// (that mistake cost 2 quanta out of 799, and test D caught it).
function collide(g) {
  const o = clone(g);
  for (let i = 0; i < N; i++) {
    const a0 = g[i * 4 + 0], a1 = g[i * 4 + 1], a2 = g[i * 4 + 2], a3 = g[i * 4 + 3];
    let o0 = 0, o1 = 0, o2 = 0, o3 = 0;
    // the (N,S) pair turns into (E,W)
    if (a0 > 0 && a2 > 0) {
      const t = Math.min(a0, a2), s = Math.abs(a0 - a2);
      o1 += t; o3 += t;
      if (s > 0) { if (a0 > a2) o0 += s; else o2 += s; }
    } else { o0 += a0; o2 += a2; }
    // the (E,W) pair turns into (N,S)
    if (a1 > 0 && a3 > 0) {
      const t = Math.min(a1, a3), s = Math.abs(a1 - a3);
      o0 += t; o2 += t;
      if (s > 0) { if (a1 > a3) o1 += s; else o3 += s; }
    } else { o1 += a1; o3 += a3; }
    o[i * 4 + 0] = o0; o[i * 4 + 1] = o1; o[i * 4 + 2] = o2; o[i * 4 + 3] = o3;
  }
  return o;
}
// move every quantum one cell along its own direction
function propagate(g) {
  const o = new Uint8Array(4 * N);
  for (let i = 0; i < N; i++) {
    const x = i % W, y = (i - x) / W;
    for (let d = 0; d < 4; d++) {
      const c = g[i * 4 + d];
      if (!c) continue;
      const k = ((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W);
      o[k * 4 + d] += c;
    }
  }
  return o;
}
function step(g) {
  let a = collide(g);
  if (MODE === 'bit') { for (let i = 0; i < a.length; i++) if (a[i] > 1) a[i] = 1; }
  a = propagate(a);
  if (MODE === 'bit') { for (let i = 0; i < a.length; i++) if (a[i] > 1) a[i] = 1; }
  return a;
}
function totals(g) {
  let n = 0, px = 0, py = 0;
  for (let i = 0; i < N; i++) for (let d = 0; d < 4; d++) {
    const c = g[i * 4 + d];
    if (!c) continue;
    n += c; px += c * DX[d]; py += c * DY[d];
  }
  return { n, px, py };
}
// a pulse of `count` quanta per cell along +x, centred at x0 on row y0, `len` cells long
function pulse(count, len, x0, y0) {
  const g = newGas();
  for (let k = 0; k < len; k++) g[idx(x0 + k, y0, 1)] = count;
  return g;
}
function pulseStats(g, y0) {
  let n = 0, sx = 0, mn = 1e9, mx = -1e9;
  for (let x = 0; x < W; x++) {
    const c = g[idx(x, y0, 1)];
    if (!c) continue;
    n += c; sx += c * x; mn = Math.min(mn, x); mx = Math.max(mx, x);
  }
  return { n, centroid: n ? sx / n : NaN, width: n ? mx - mn + 1 : 0 };
}

// ============================================================ A: coexistence
console.log('============================================================');
console.log('A: CAN A MODE HOLD TWO QUANTA?');
console.log('============================================================');
{
  let ok = true;
  console.log('   quanta placed in one mode   count alphabet   bit alphabet (old)');
  for (const c of [1, 2, 3, 5]) {
    const g = newGas(); g[idx(10, 12, 1)] = c;
    MODE = 'count'; const a = step(g);
    MODE = 'bit'; const b = step(g);
    MODE = 'count';
    const ca = totals(a).n, cb = totals(b).n;
    const good = (c === 1 ? ca === 1 && cb === 1 : ca === c && cb === 1);
    if (!good) ok = false;
    console.log(`        ${c}                        ${String(ca).padStart(2)}               ${String(cb).padStart(2)}`);
  }
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the count alphabet carries the number through; the bit alphabet collapses`);
  console.log(`     everything to one, because "one more quantum" is 1|1 = 1 there.  That single line is the`);
  console.log(`     whole of round 35's alphabet argument, and it is now measured rather than argued.`);
}
console.log('');

// ============================================================ B: linearity
console.log('============================================================');
console.log('B: IS THE GAS NOW LINEAR? (the round-35 ruler, on a wider alphabet)');
console.log('============================================================');
{
  // build a few configurations, including ones that collide, and check step(2g) == 2*step(g)
  const cases = [];
  {
    const g = newGas();
    g[idx(20, 12, 1)] = 1; g[idx(30, 12, 3)] = 1;         // head-on pair, one quantum each
    cases.push(['single head-on pair', g]);
  }
  {
    const g = newGas();
    g[idx(20, 12, 1)] = 2; g[idx(30, 12, 3)] = 1;         // a head-on pair with a surplus
    cases.push(['pair with a surplus', g]);
  }
  {
    const g = newGas();
    for (let x = 10; x < 16; x++) g[idx(x, 12, 1)] = 1;   // a pulse
    for (let x = 30; x < 34; x++) g[idx(x, 12, 3)] = 1;   // meeting it
    cases.push(['pulse meets a pulse', g]);
  }
  let ok = true;
  console.log('   configuration            steps   max |step(2g) - 2*step(g)|');
  for (const [name, g] of cases) {
    let a = g, b = clone(g);
    for (let i = 0; i < b.length; i++) b[i] *= 2;
    let worst = 0;
    for (let t = 1; t <= 12; t++) {
      a = step(a); b = step(b);
      for (let i = 0; i < a.length; i++) worst = Math.max(worst, Math.abs(2 * a[i] - b[i]));
    }
    if (worst !== 0) ok = false;
    console.log(`   ${name.padEnd(24)} 12      ${worst}`);
  }
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: doubling the state doubles the evolution EXACTLY, collisions included --`);
  console.log(`     because min(m,n) and |m-n| are homogeneous.  Round 35 asked for "can it be superposed";`);
  console.log(`     for positive integer amplitudes the answer is now yes.`);
  console.log(`     HONEST LIMIT: counts are non-negative, so this is a CONE, not a vector space.  There are no`);
  console.log(`     negative or complex amplitudes, so there is still no interference.  The gap has moved from`);
  console.log(`     "no |2>" to "no signs" -- which is the next thing to attack, and it is a different thing.`);
}
console.log('');

// ============================================================ C: light
console.log('============================================================');
console.log('C: IS THE LIGHT STILL IDEAL? (speed exactly 1, no spreading, any amplitude)');
console.log('============================================================');
{
  console.log('   quanta per cell   speed (cells/step)   width at t=0 -> t=40   number at t=0 -> t=40');
  let ok = true;
  for (const c of [1, 2, 3, 5]) {
    let g = pulse(c, 6, 10, 12);
    const s0 = pulseStats(g, 12);
    for (let t = 0; t < 40; t++) g = step(g);
    const s1 = pulseStats(g, 12);
    const speed = (s1.centroid - s0.centroid) / 40;
    if (Math.abs(speed - 1) > 1e-9 || s1.width !== s0.width || s1.n !== s0.n) ok = false;
    console.log(`        ${c}              ${speed.toFixed(4)}              ${s0.width} -> ${s1.width}                ${s0.n} -> ${s1.n}`);
  }
  console.log('');
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: the speed is exactly the refresh rate and the pulse does not spread, at every`);
  console.log(`     amplitude.  Gapless, non-dispersive light survives the alphabet change, and the number is`);
  console.log(`     now a real amplitude rather than a bit.`);
}
console.log('');

// ============================================================ D: conservation
console.log('============================================================');
console.log('D: IS ANYTHING LEAKED? (a long, collision-heavy run)');
console.log('============================================================');
{
  let g = newGas();
  for (let i = 0; i < 400; i++) {
    const x = 1 + ((i * 37) % (W - 2)), y = 1 + ((i * 13) % (H - 2)), d = i % 4;
    const c = 1 + (i % 3);
    g[idx(x, y, d)] += c;
  }
  const t0 = totals(g);
  let ok = true;
  for (let t = 0; t < 300; t++) {
    g = step(g);
    const tt = totals(g);
    if (tt.n !== t0.n || tt.px !== t0.px || tt.py !== t0.py) { ok = false; break; }
  }
  const t1 = totals(g);
  console.log(`  start: ${t0.n} quanta, momentum (${t0.px},${t0.py})   after 300 steps: ${t1.n} quanta, momentum (${t1.px},${t1.py})`);
  console.log(`  => ${ok ? 'PASS' : '*** FAIL ***'}: number and both momenta are conserved at every step of a 400-quantum,`);
  console.log(`     collision-heavy run.  The conservative extension of the turn rule works.`);
}
console.log('');

// ============================================================ E: the amplification target
console.log('============================================================');
console.log('E: WHAT THE ALPHABET WOULD BUY A LASER (computed, not measured)');
console.log('============================================================');
{
  console.log('   quanta already in the mode   bit alphabet: gain   count alphabet: gain');
  for (const n of [1, 2, 3, 10, 100]) {
    console.log(`        ${String(n).padStart(3)}                       1.000               ${((n + 1) / n).toFixed(3)}`);
  }
  console.log('');
  console.log(`  => round 56 measured the left column (1.000 at every n, because the mode saturates at one)`);
  console.log(`     and named the right column as the target.  The target is now REPRESENTABLE: a mode holds`);
  console.log(`     n quanta, so a stimulated emission adds one and the gain is 1 + 1/n.`);
  console.log(`     NOT DONE HERE: the molecule is not rewired to this alphabet in this round, so this column`);
  console.log(`     is arithmetic, not a measurement.  That is the next step, and it is a wiring job rather`);
  console.log(`     than an alphabet question.`);
}
