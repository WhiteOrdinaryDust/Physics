'use strict';
/*
 * molecule_photon_coupling.js
 * ---------------------------------------------------------------------------
 * THE MISSING PIECE: a coupling channel between the photon and the molecule.
 *
 * Until now the molecule was TRANSPARENT: a photon went through and neither side changed.
 * Now it gets a one-bit memory with a WRITE and a READ:
 *
 *     record = 0   a photon arrives -> the photon passes UNDISTURBED and the record becomes 1
 *     record = 1   a photon arrives -> the photon is REFLECTED and the record becomes 0
 *     no photon    -> the record simply persists (the molecule's clock keeps ticking)
 *
 * So the first photon is recorded while paying nothing, and the second photon's fate reveals
 * what the first one did. That is a memory cell written by light and read by light -- the
 * "toggle" from the original functional-device list.
 *
 * WHY THIS MATTERS FOR THE PROGRAMME
 *   Round 34 tried to find a wall that a particle can touch without leaving a record, by
 *   brute-force searching cell patterns, and failed; the conclusion then was that a wall built
 *   from "cell matter" is always a detector. Here the same thing is reached the other way round,
 *   by DESIGN: the molecule records the passage while the messenger is left completely
 *   unperturbed. That is exactly the user's "construct the devices by function, do not search".
 *
 * WHY IT IS ALSO THE HONEST PLACE TO STOP AND LOOK
 *   The record rule as written makes the molecule a FIXED MIRROR: it reflects a photon without
 *   recoiling, so momentum is NOT conserved when the record is read. (An infinitely heavy mirror
 *   is legitimate physics; a two-particle molecule is not infinitely heavy.) Giving it recoil is
 *   the next design step and requires a molecule that can MOVE. I report this rather than hide it.
 *
 * PREDICTIONS (before the run)
 *   C1  controls: a photon alone still travels at speed 1; a lone molecule is still static with a
 *       period-2 clock; the empty vacuum is still a fixed point.
 *   C2  THE WRITE: a photon sent at a molecule with record 0 passes THROUGH unchanged (same speed,
 *       same direction, no delay longer than one cell) and flips the record to 1.
 *   C3  THE RECORD PERSISTS: left alone for many steps, the record bit does not change.
 *   C4  THE READ: the next photon is REFLECTED (it leaves in the direction it came from), and the
 *       record flips back to 0. So the molecule alternates pass / reflect / pass / reflect.
 *   C5  the record channel itself is invertible (its preimage is unique), so it adds no new
 *       irreversibility; the only information-destroying step remains the BINDING collision,
 *       which the previous round already proved non-injective with a concrete counterexample.
 * ---------------------------------------------------------------------------
 */

const W = 64, H = 32, N = W * H;
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
const PBITS = 0xff;
const MBASE = 1 << 8, MOR1 = 1 << 9, MOR2 = 1 << 10, MCLK = 1 << 11, MREC = 1 << 12;
function markerOf(o) { return MBASE | ((o & 1) ? MOR1 : 0) | ((o & 2) ? MOR2 : 0); }
const TURN = new Map([[1 | 4, 2 | 8], [2 | 8, 1 | 4], [16 | 64, 32 | 128], [32 | 128, 16 | 64]]);
const BIND = new Map();
BIND.set(1 | 64, markerOf(0));
BIND.set(4 | 16, markerOf(1));
BIND.set(2 | 128, markerOf(2));
BIND.set(8 | 32, markerOf(3));
const UNBIND = new Map();
for (const [pair, mark] of BIND) UNBIND.set(mark, pair);
function reflect(parts) {
  let np = 0;
  for (let d = 0; d < 4; d++) {
    if ((parts >> d) & 1) np |= 1 << ((d + 2) % 4);
    if ((parts >> (d + 4)) & 1) np |= 1 << (((d + 2) % 4) + 4);
  }
  return np;
}

function collide(g) {
  const o = new Uint16Array(N);
  for (let i = 0; i < N; i++) {
    const v = g[i];
    const parts = v & PBITS, mol = v & ~PBITS;
    if (mol) {
      const m = mol ^ MCLK;                          // the molecule's clock always ticks
      if (!parts) { o[i] = m; continue; }            // alone: the record just persists
      if (m & MREC) o[i] = (m & ~MREC) | reflect(parts);   // READ: reflect, clear the record
      else o[i] = (m | MREC) | parts;                      // WRITE: the photon passes untouched
      continue;
    }
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
      const m = mol ^ MCLK;                          // undo the clock tick
      if (!parts) { o[i] = m; continue; }
      if (m & MREC) o[i] = (m & ~MREC) | parts;              // post-write -> pre had no record
      else o[i] = (m | MREC) | reflect(parts);               // post-read  -> pre had the record
      continue;
    }
    if (TURN.has(parts)) { o[i] = TURN.get(parts); continue; }
    const back = UNBIND.get(parts);
    o[i] = (back !== undefined ? back : parts);
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
    if (v & ~PBITS) { mol++; ca++; cb++; }
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
  return { n, w: hi - lo + 1, h: hi2 - lo2 + 1, cx: sx / n, cy: sy / n };
}
// build a molecule at (mx,16) by colliding an A-E and a B-W two cells apart
function withMolecule(mx) {
  const g = new Uint16Array(N);
  g[16 * W + (mx - 2)] = 1 << 1;          // A-E
  g[16 * W + (mx + 2)] = 1 << 7;          // B-W
  let a = Uint16Array.from(g);
  for (let t = 0; t < 3; t++) a = step(a);
  return a;
}

// ============================================================ controls
console.log('============================================================');
console.log('CONTROLS');
console.log('============================================================');
{
  const g = new Uint16Array(N);
  g[16 * W + 10] = 1 << 1;
  let a = Uint16Array.from(g);
  const x0 = bbox(a).cx;
  for (let t = 0; t < 20; t++) a = step(a);
  const x1 = bbox(a).cx;
  console.log(`  C1 a photon alone: speed ${((x1 - x0) / 20).toFixed(3)}, mass ${popAll(a)} / 1  ${Math.abs((x1 - x0) / 20 - 1) < 1e-9 && popAll(a) === 1 ? 'PASS' : '*** FAIL ***'}`);
  const m = withMolecule(22);
  const b = bbox(m);
  const rec0 = m[16 * W + 22] & MREC;
  let mm = Uint16Array.from(m);
  for (let t = 0; t < 40; t++) mm = step(mm);
  const b2 = bbox(mm);
  console.log(`  C1 a lone molecule: centre ${b.cx.toFixed(2)} -> ${b2.cx.toFixed(2)} over 40 steps (must not drift)  ${Math.abs(b2.cx - b.cx) < 1e-9 ? 'PASS' : '*** FAIL ***'}`);
  console.log(`     its record after 40 steps alone: ${(mm[16 * W + 22] & MREC) ? 1 : 0} (was ${rec0 ? 1 : 0})  => C3 record persists: ${((mm[16 * W + 22] & MREC) ? 1 : 0) === (rec0 ? 1 : 0) ? 'PASS' : '*** FAIL ***'}`);
  const e = new Uint16Array(N); let ee = e;
  for (let t = 0; t < 30; t++) ee = step(ee);
  console.log(`  C1 empty vacuum fixed: ${popAll(ee) === 0 ? 'PASS' : '*** FAIL ***'}`);
}
console.log('');

// ============================================================ the write / read cycle
console.log('============================================================');
console.log('C2/C4: write with one photon, read with the next');
console.log('============================================================');
function runSequence(label) {
  let g = withMolecule(22);
  const MX = 22;
  let t = 0;
  console.log(`  ${label}`);
  console.log('    photon#  outcome                       record  photon direction  steps to travel 12 cells');
  for (let shot = 1; shot <= 3; shot++) {
    g[16 * W + 10] = 1 << 1;                     // fire an A-E photon from x=10
    let firstMeet = -1, dirAfter = null, traveled = 0, prevX = 10;
    for (let k = 0; k < 40; k++) {
      g = step(g); t++;
      const b = bbox(g);
      // find the free photon (the A particle that is not part of the molecule)
      let px = -1;
      for (let x = 0; x < W; x++) { const v = g[16 * W + x]; if (v & (1 << 1)) px = x; }
      if (firstMeet < 0 && px >= 0 && Math.abs(px - MX) <= 1 && k > 2) { firstMeet = k; }
      if (firstMeet >= 0 && k === firstMeet + 6) { dirAfter = px - MX; break; }
      void b; void traveled; void prevX;
    }
    const rec = (g[16 * W + MX] & MREC) ? 1 : 0;
    const outcome = dirAfter === null ? '(never reached the molecule)' : (dirAfter > 0 ? 'PASSED THROUGH (continues +x)' : 'REFLECTED (goes back -x)');
    console.log(`    ${String(shot).padStart(6)}  ${outcome.padEnd(30)}  ${rec}       ${dirAfter === null ? '-' : (dirAfter > 0 ? '+x' : '-x')}`);
    // clear the outgoing photon before the next shot: on a torus a reflected photon would come
    // round again and confuse the bookkeeping. This is a controlled intervention, stated openly.
    for (let i = 0; i < N; i++) g[i] &= ~PBITS;
  }
  const c = counts(g);
  console.log(`    final: molecule ${c.mol}, A=${c.ca}, B=${c.cb}, momentum (${c.px},${c.py})`);
  return g;
}
runSequence('three photons fired at the same molecule in turn:');
console.log('');
console.log('  READING: if the molecule were still transparent, every photon would pass and the');
console.log('  record would never change. If the memory works, the sequence must alternate.');
console.log('');

// ============================================================ does it really record?
console.log('============================================================');
console.log('the decisive comparison: with a photon vs without');
console.log('============================================================');
{
  let withP = withMolecule(22);
  withP[16 * W + 10] = 1 << 1;
  let noP = withMolecule(22);
  console.log('     t   record WITH a photon   record WITHOUT a photon');
  for (let t = 0; t <= 16; t++) {
    console.log(`    ${String(t).padStart(3)}        ${(withP[16 * W + 22] & MREC) ? 1 : 0}                        ${(noP[16 * W + 22] & MREC) ? 1 : 0}`);
    withP = step(withP); noP = step(noP);
  }
  console.log('  => the two columns differ exactly when the photon is at the molecule: that is the record.');
}
console.log('');

// ============================================================ does the coupling conserve momentum?
console.log('============================================================');
console.log('C5: what the coupling costs, and whether it adds irreversibility');
console.log('============================================================');
{
  let g = withMolecule(22);
  const before = counts(g);
  g[16 * W + 10] = 1 << 1;
  let recAt = -1;
  for (let t = 0; t < 60; t++) { g = step(g); if (recAt < 0 && (g[16 * W + 22] & MREC)) recAt = t; }
  const after = counts(g);
  console.log(`  momentum before the shot (${before.px},${before.py}), after ${after.px},${after.py} steps ${after.py - before.py === 0 ? '' : ''}(${after.px - before.px},${after.py - before.py})`);
  console.log(`    => ${after.px === before.px && after.py === before.py ? 'conserved' : 'NOT conserved: the molecule reflected a photon without recoiling, i.e. it acts as a FIXED MIRROR.'}`);
  console.log(`    (an infinitely heavy mirror may do that; a two-particle molecule may not. Recoil needs a`);
  console.log(`     molecule that can MOVE, which is the next design step.)`);
  // invertibility of the record channel itself
  const ALPHA = [0];
  for (let d = 0; d < 4; d++) { ALPHA.push(1 << d); ALPHA.push(1 << (d + 4)); }
  for (let o = 0; o < 4; o++) { ALPHA.push(markerOf(o)); ALPHA.push(markerOf(o) | MCLK); ALPHA.push(markerOf(o) | MREC); ALPHA.push(markerOf(o) | MCLK | MREC); }
  for (const pair of BIND.keys()) ALPHA.push(pair);
  for (const pair of TURN.keys()) ALPHA.push(pair);
  const S = 2, NC = S * S;
  const oneStep = (cfg) => {
    const grid = new Uint16Array(NC);
    for (let i = 0; i < NC; i++) grid[i] = ALPHA[cfg[i]];
    // collide (same logic as above, on a 2x2 torus)
    const mid = new Uint16Array(NC);
    for (let i = 0; i < NC; i++) {
      const v = grid[i], parts = v & PBITS, mol = v & ~PBITS;
      if (mol) {
        const m = mol ^ MCLK;
        if (!parts) mid[i] = m;
        else if (m & MREC) mid[i] = (m & ~MREC) | reflect(parts);
        else mid[i] = (m | MREC) | parts;
        continue;
      }
      if (TURN.has(parts)) { mid[i] = TURN.get(parts); continue; }
      if (BIND.has(parts)) { mid[i] = BIND.get(parts); continue; }
      mid[i] = parts;
    }
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
  const kinds = { binding: 0, other: 0 };
  const examples = [];
  const cfg = new Array(NC).fill(0);
  const rec = (i) => {
    if (i === NC) {
      const img = oneStep(cfg);
      const key = cfg.join(',');
      if (seen.has(img) && seen.get(img) !== key) {
        collisions++;
        const A = seen.get(img).split(',').map(Number), B = cfg.slice();
        const isBind = A.includes(ALPHA.indexOf(BIND.get(1 | 64))) || B.includes(ALPHA.indexOf(BIND.get(1 | 64)));
        // classify: does one of them contain a same-cell opposite-colour pair (a binding state)?
        const hasPair = (arr) => arr.some(k => BIND.has(ALPHA[k]));
        if (hasPair(A) || hasPair(B)) kinds.binding++; else kinds.other++;
        if (examples.length < 2) examples.push([A, B, img]);
      } else seen.set(img, key);
      return;
    }
    for (let a = 0; a < ALPHA.length; a++) { cfg[i] = a; rec(i + 1); }
  };
  rec(0);
  console.log('');
  console.log(`  injectivity of the NEW rule on a 2x2 torus (alphabet now includes the record states):`);
  console.log(`    ${seen.size + collisions} configurations, ${collisions} pairs share an image`);
  console.log(`    of those: ${kinds.binding} involve a same-cell pair (the binding collision), ${kinds.other} do not`);
  for (const [A, B, img] of examples) {
    console.log(`    example: [${A.map(k => ALPHA[k].toString(2).padStart(13, '0')).join(' ')}]`);
    console.log(`             [${B.map(k => ALPHA[k].toString(2).padStart(13, '0')).join(' ')}]  both -> [${img.split(',').map(v => Number(v).toString(2).padStart(13, '0')).join(' ')}]`);
  }
  console.log(`    => ${kinds.other === 0 ? 'the record channel adds NO new irreversibility: everything that collides is the binding step' : 'the record channel ALSO destroys information'}`);
}
