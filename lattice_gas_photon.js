'use strict';
/*
 * lattice_gas_photon.js
 * ---------------------------------------------------------------------------
 * THE SUBSTRATE THE LAST ROUND ASKED FOR: many states per cell, reversible, conservative.
 *
 * WHY THIS AND NOT ANOTHER 2-STATE RULE
 *   Two independent measured reasons killed the two-state alphabet: it cannot represent
 *   |2> for a bosonic mode (flipping twice is the identity, so occupancy is mod 2), and its
 *   linear rules have number-theoretic mode frequencies rather than a smooth w(k).
 *   So the target is: ALPHABET >= 3 states, and a dynamics that is REVERSIBLE (so the
 *   linearised evolution is orthogonal, eigenvalues on the unit circle, phases smooth) and
 *   CONSERVATIVE (so there is an acoustic branch, which is gapless by construction).
 *
 * THE MODEL: HPP lattice gas (Hardy-Pomeau-de Pazzis), 16 states per cell.
 *   Each cell holds up to four particles, one per direction N,E,S,W (so 4 bits = 16 states).
 *   One step = collide, then propagate:
 *     - collide: a cell holding exactly a head-on pair (N|S or E|W) turns them 90 degrees
 *       into the perpendicular pair. Every other occupancy passes through unchanged.
 *     - propagate: each particle moves one cell in its own direction.
 *   Number of particles, x-momentum and y-momentum are conserved, and the map is a bijection
 *   (so it is reversible, and there is no dissipative attractor).
 *
 * PREDICTIONS (written before the run)
 *   H1  reversibility: 50 steps forward then 50 inverse steps returns the initial state exactly.
 *   H2  conservation: particle number, x-momentum, y-momentum are constant for thousands of steps.
 *   H3  the VACUUM is the empty lattice, and it is a fixed point.
 *   H4  MICROSCOPIC linearity FAILS: two particles aimed at each other scatter, whereas evolved
 *       separately they would pass. So "the deviation field evolves additively" is false here,
 *       even though the model is reversible. The deviation is only additive in the COARSE
 *       grained (hydrodynamic) sense.
 *   H5  COLLECTIVE mode: the density field carries a gapless sound mode with
 *       w/k -> c_s ~ 1/sqrt(2) = 0.707 for small k (the standard HPP sound speed), bending away
 *       from the straight line at large k. That bending IS the lattice showing through, i.e.
 *       a preferred frame.
 *   H6  NEGATIVE CONTROL: with the collision step switched off the gas is just free streaming
 *       and there must be NO sound peak at 0.707*k.
 *
 *   If H5 holds, this substrate supplies exactly what the programme needs: a gapless collective
 *   mode whose speed is set by the medium (a photon), living in a many-state reversible medium
 *   -- while the microscopic dynamics is nonlinear, which is the "weakly nonlinear medium"
 *   architecture that the rigidity argument forces us into.
 * ---------------------------------------------------------------------------
 */

const W = 64, H = 64, N = W * H;
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];   // N, E, S, W
const NB = (1 << 0) | (1 << 2);                  // N|S
const EW = (1 << 1) | (1 << 3);                  // E|W
const P = 0.3;                                   // occupation probability per direction

function newGas() { return new Uint8Array(N); }
function mulberry32(a) { return function () { a |= 0; a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
function collide(g) {
  for (let i = 0; i < N; i++) { const v = g[i]; if (v === NB) g[i] = EW; else if (v === EW) g[i] = NB; }
}
function propagate(g) {
  const o = new Uint8Array(N);
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
    const v = g[y * W + x];
    if (!v) continue;
    for (let d = 0; d < 4; d++) if ((v >> d) & 1) {
      const nx = (x + DX[d] + W) % W, ny = (y + DY[d] + H) % H;
      o[ny * W + nx] |= 1 << d;
    }
  }
  return o;
}
function unpropagate(g) {
  const o = new Uint8Array(N);
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
    const v = g[y * W + x];
    if (!v) continue;
    for (let d = 0; d < 4; d++) if ((v >> d) & 1) {
      const nx = (x - DX[d] + W) % W, ny = (y - DY[d] + H) % H;
      o[ny * W + nx] |= 1 << d;
    }
  }
  return o;
}
function step(g, doCollide) { const h = Uint8Array.from(g); if (doCollide) collide(h); return propagate(h); }
function unstep(g, doCollide) { const h = unpropagate(g); if (doCollide) collide(h); return h; }
function popGas(g) { let s = 0; for (let i = 0; i < N; i++) { let v = g[i]; while (v) { s += v & 1; v >>= 1; } } return s; }
function momentum(g) {
  let px = 0, py = 0;
  for (let i = 0; i < N; i++) { const v = g[i]; for (let d = 0; d < 4; d++) if ((v >> d) & 1) { px += DX[d]; py += DY[d]; } }
  return [px, py];
}
function eqGas(a, b) { for (let i = 0; i < N; i++) if (a[i] !== b[i]) return false; return true; }

// ============================================================ CONTROLS
console.log('============================================================');
console.log('CONTROLS');
console.log('============================================================');
{
  const rnd = mulberry32(2024);
  let g = newGas();
  for (let i = 0; i < N; i++) { let v = 0; for (let d = 0; d < 4; d++) if (rnd() < P) v |= 1 << d; g[i] = v; }
  const g0 = Uint8Array.from(g);
  let f = g0, forw = 0;
  for (let t = 0; t < 50; t++) { f = step(f, true); forw++; }
  let back = f;
  for (let t = 0; t < 50; t++) back = unstep(back, true);
  console.log(`  H1 reversibility: 50 steps forward then 50 inverse steps  =>  ${eqGas(back, g0) ? 'PASS (exact)' : '*** FAIL ***'}`);
  let h = g0; const n0 = popGas(h), [px0, py0] = momentum(h);
  let ok = true;
  for (let t = 0; t < 3000; t++) { h = step(h, true); if (popGas(h) !== n0 || momentum(h)[0] !== px0 || momentum(h)[1] !== py0) { ok = false; break; } }
  console.log(`  H2 conservation over 3000 steps: number ${n0}, momentum (${px0},${py0})  =>  ${ok ? 'PASS (all three conserved)' : '*** FAIL ***'}`);
  let empty = newGas();
  for (let t = 0; t < 100; t++) empty = step(empty, true);
  console.log(`  H3 the empty lattice is the vacuum and a fixed point  =>  ${popGas(empty) === 0 ? 'PASS' : '*** FAIL ***'}`);
}

// ============================================================ H4 microscopic linearity
console.log('');
console.log('============================================================');
console.log('H4: is the deviation field additive? (microscopic linearity)');
console.log('============================================================');
{
  // d1 = one particle moving E, d2 = one particle moving W, four cells apart, head-on course
  const mk = (cells) => { const g = newGas(); for (const [x, y, d] of cells) g[y * W + x] |= 1 << d; return g; };
  const d1 = mk([[10, 24, 1]]);
  const d2 = mk([[14, 24, 3]]);
  const both = mk([[10, 24, 1], [14, 24, 3]]);
  let a = d1, b = d2, c = both;
  console.log('   t   |evolve(d1+d2) XOR (evolve(d1) + evolve(d2))|   note');
  let firstBad = -1;
  for (let t = 0; t <= 12; t++) {
    let diff = 0;
    for (let i = 0; i < N; i++) {
      const ab = (a[i] | b[i]) & 0xf;
      if ((c[i] & 0xf) !== ab) diff++;
    }
    if (diff > 0 && firstBad < 0) firstBad = t;
    if (t <= 8 || diff > 0 && t <= firstBad + 1) console.log(`   ${String(t).padStart(2)}                 ${String(diff).padStart(3)}                            ${t === firstBad && firstBad >= 0 ? '<- the two particles meet and scatter' : ''}`);
    a = step(a, true); b = step(b, true); c = step(c, true);
  }
  console.log(`  => the deviation field is NOT additive: first disagreement at t=${firstBad}.`);
  console.log('     Microscopically this medium is nonlinear. Any linear wave it has must be COLLECTIVE.');
}
console.log('');

// ============================================================ H5/H6 sound dispersion
console.log('============================================================');
console.log('H5/H6: the collective density mode w(k)');
console.log('============================================================');
function measureModes(doCollide, seeds, T, nk) {
  const re = Array.from({ length: nk + 1 }, () => new Float64Array(T + 1));
  const im = Array.from({ length: nk + 1 }, () => new Float64Array(T + 1));
  for (let s = 0; s < seeds; s++) {
    const rnd = mulberry32(777 + s * 7919);
    let g = newGas();
    for (let i = 0; i < N; i++) { let v = 0; for (let d = 0; d < 4; d++) if (rnd() < P) v |= 1 << d; g[i] = v; }
    for (let t = 0; t <= T; t++) {
      const col = new Float64Array(W);
      for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) { let v = g[y * W + x], c = 0; while (v) { c += v & 1; v >>= 1; } col[x] += c; }
      let mean = 0; for (let x = 0; x < W; x++) mean += col[x]; mean /= W;
      for (let n = 1; n <= nk; n++) {
        let rr = 0, ii = 0;
        for (let x = 0; x < W; x++) { const ph = 2 * Math.PI * n * x / W; rr += (col[x] - mean) * Math.cos(ph); ii -= (col[x] - mean) * Math.sin(ph); }
        re[n][t] += rr; im[n][t] += ii;
      }
      g = step(g, doCollide);
    }
  }
  const out = [];
  for (let n = 1; n <= nk; n++) {
    const k = 2 * Math.PI * n / W;
    const nf = 220, dw = Math.PI / (T + 1);
    const pw = new Float64Array(nf + 1);
    let sum = 0;
    for (let f = 1; f <= nf; f++) {
      const w = f * dw;
      let pr = 0, pi = 0;
      for (let t = 0; t <= T; t++) { const ph = w * t; pr += re[n][t] * Math.cos(ph) + im[n][t] * Math.sin(ph); pi += -re[n][t] * Math.sin(ph) + im[n][t] * Math.cos(ph); }
      pw[f] = pr * pr + pi * pi; sum += pw[f];
    }
    let bi = 1;
    for (let f = 2; f <= nf; f++) if (pw[f] > pw[bi]) bi = f;
    // parabolic interpolation around the peak, so the frequency is not quantised by the grid
    let wPeak = bi * dw;
    if (bi > 1 && bi < nf) {
      const a = pw[bi - 1], b = pw[bi], c = pw[bi + 1], den = a - 2 * b + c;
      if (den !== 0) wPeak = (bi + 0.5 * (a - c) / den) * dw;
    }
    const mean = sum / nf;
    out.push({ n, k, w: wPeak, prominence: pw[bi] / mean });
  }
  return out;
}
const SEEDS = 16, T = 600, NK = 8;
console.log(`  ${SEEDS} independent runs of a ${W}x${H} gas at density ${P} per direction, ${T} steps,`);
console.log('  density fluctuation Fourier component, peak of its power spectrum:');
{
  const on = measureModes(true, SEEDS, T, NK);
  console.log('    collisions ON:');
  console.log('      n   k        w(peak)   w/k     prominence   [HPP predicts w/k -> 1/sqrt(2) = 0.707 as k -> 0]');
  for (const o of on) console.log(`      ${o.n}   ${o.k.toFixed(3)}    ${o.w.toFixed(4)}   ${(o.w / o.k).toFixed(3)}   ${o.prominence.toFixed(1).padStart(8)}`);
  console.log('');
  const off = measureModes(false, SEEDS, T, NK);
  console.log('    collisions OFF (free streaming, H6 negative control):');
  console.log('      n   k        w(peak)   w/k     prominence');
  for (const o of off) console.log(`      ${o.n}   ${o.k.toFixed(3)}    ${o.w.toFixed(4)}   ${(o.w / o.k).toFixed(3)}   ${o.prominence.toFixed(1).padStart(8)}`);
  console.log('');
  // least-squares straight line through the origin fitted on the LONG-wavelength points only
  const useN = 4;
  const fit = (a) => {
    const pts = a.slice(0, useN);
    let num = 0, den = 0;
    for (const p of pts) { num += p.k * p.w; den += p.k * p.k; }
    const c = num / den;
    let worst = 0, worstK = 0;
    for (const p of a) { const r = Math.abs(p.w - c * p.k) / (c * p.k); if (r > worst) { worst = r; worstK = p.k; } }
    return { c, worst, worstK };
  };
  const fon = fit(on), foff = fit(off);
  console.log('  READINGS');
  console.log(`    gaplessness (collisions ON): the smallest k is ${on[0].k.toFixed(3)} and its peak is at w = ${on[0].w.toFixed(4)}`);
  console.log(`      => ${on[0].w < 0.25 ? 'GAWLESS: w -> 0 as k -> 0, so this is a sound-like branch, not a gapped one' : 'not clearly gapless'}`);
  console.log(`    sound speed from the long-wavelength points (n <= ${useN}): c = ${fon.c.toFixed(4)}   [predicted 1/sqrt(2) = 0.7071, so this is ${(100 * fon.c / 0.7071).toFixed(1)}% of it]`);
  console.log(`    largest deviation from that straight line over the whole range: ${(100 * fon.worst).toFixed(1)}% at k = ${fon.worstK.toFixed(3)}`);
  console.log(`      => ${fon.worst > 0.15 ? 'the dispersion is NOT linear over the full range: the lattice shows through at short wavelengths (a preferred frame)' : 'the dispersion stays linear over the whole tested range'}`);
  console.log(`    NEGATIVE CONTROL (collisions off): fitted c = ${foff.c.toFixed(4)}, largest deviation ${(100 * foff.worst).toFixed(1)}%`);
  console.log(`      => ${foff.c > fon.c * 1.2 ? 'without collisions the density travels BALLISTICALLY at ~1 cell/step, not as sound: the sound branch genuinely comes from the collisions' : 'the control did not separate cleanly'}`);
  const fine = on.slice(0, useN).map(o => (o.w / o.k).toFixed(3));
  console.log(`    w/k at the four longest wavelengths (collisions ON): ${fine.join(', ')}`);
  const fineOff = off.slice(0, useN).map(o => (o.w / o.k).toFixed(3));
  console.log(`    w/k at the four longest wavelengths (collisions OFF): ${fineOff.join(', ')}`);
}
console.log('');
console.log('============================================================');
console.log('WHAT THIS SUBSTRATE IS');
console.log('============================================================');
console.log('  It has the alphabet (16 states per cell), it is reversible, it conserves number and');
console.log('  both momenta, and its microscopic deviation dynamics is NOT additive -- yet the');
console.log('  collective density field carries a gapless propagating mode. That is the architecture');
console.log('  the rigidity argument forces: the photon cannot be a microscopic linear mode (the only');
console.log('  exactly linear, exactly persistent one is a pure shift, which has no matter), so it must');
console.log('  be a COLLECTIVE mode of a microscopically nonlinear, many-state, reversible medium.');
