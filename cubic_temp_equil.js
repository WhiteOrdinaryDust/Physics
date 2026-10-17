'use strict';
/*
 * cubic_temp_equil.js  --  ROUND 109 DIAGNOSTIC: why the cross-check failed, and how many sweeps it needs
 * ---------------------------------------------------------------------------
 * cubic_temperature.js failed its cross-check against round 95: at lambda = 2, T = 0.7, L = 12 it gave
 * B2 = 1.6162e-1 +- 4.95e-4 where round 95 recorded 1.0240e-1 +- 8.72e-3.  The statistical error is far too
 * small to explain that, so it is a systematic difference, and the first suspect is EQUILIBRATION: the
 * initial configuration is aligned along n_hat, so any transverse anisotropy has to GROW in, and my
 * measurement window (sweeps 1500 to 3000) may simply be too early.
 *
 * This script measures B2 against the number of sweeps, with the measurement window always the second
 * half, so the question "where does it settle?" is answered directly.  If it drifts down toward 0.1024 the
 * diagnosis is confirmed; if it settles somewhere else, then the difference from round 95 lies elsewhere
 * and that is what the next step must look at.
 *
 * PRE-REGISTERED PREDICTIONS
 *   Q1 B2 falls as the number of sweeps grows, i.e. the first run was not equilibrated
 *   Q2 it approaches round 95's 1.0240e-1
 *   Q3 if Q1 fails, the difference is NOT equilibration and the round-95 script must be re-read
 * ---------------------------------------------------------------------------
 */

const TH = (22.5 * Math.PI) / 180;
const NHAT = [Math.cos(TH), Math.sin(TH), 0];
const E2 = [0, 0, 1];
const E1 = (() => { const a = [E2[1] * NHAT[2] - E2[2] * NHAT[1], E2[2] * NHAT[0] - E2[0] * NHAT[2], E2[0] * NHAT[1] - E2[1] * NHAT[0]]; const n = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / n, a[1] / n, a[2] / n]; })();
function build(L) {
  const N = L * L * L;
  const idx = (x, y, z) => ((x % L) + L) % L + L * ((((y % L) + L) % L) + L * (((z % L) + L) % L));
  const nbrs = [];
  for (let x = 0; x < L; x++) for (let y = 0; y < L; y++) for (let z = 0; z < L; z++) {
    nbrs[idx(x, y, z)] = [idx(x + 1, y, z), idx(x - 1, y, z), idx(x, y + 1, z), idx(x, y - 1, z), idx(x, y, z + 1), idx(x, y, z - 1)];
  }
  return { L, N, idx, nbrs };
}
function rotMat(axis, ang) {
  const [ux, uy, uz] = axis, c = Math.cos(ang), sn = Math.sin(ang), t = 1 - c;
  return [t * ux * ux + c, t * ux * uy - sn * uz, t * ux * uz + sn * uy,
          t * ux * uy + sn * uz, t * uy * uy + c, t * uy * uz - sn * ux,
          t * ux * uz - sn * uy, t * uy * uz + sn * ux, t * uz * uz + c];
}
function applyMat(s, m) {
  const out = s.slice();
  for (let i = 0; i < s.length / 3; i++) {
    const a = s[3 * i], b = s[3 * i + 1], c = s[3 * i + 2];
    out[3 * i] = m[0] * a + m[1] * b + m[2] * c;
    out[3 * i + 1] = m[3] * a + m[4] * b + m[5] * c;
    out[3 * i + 2] = m[6] * a + m[7] * b + m[8] * c;
  }
  return out;
}
function gaugeFix(g, s, nhat) {
  let mx = 0, my = 0, mz = 0;
  for (let i = 0; i < g.N; i++) { mx += s[3 * i]; my += s[3 * i + 1]; mz += s[3 * i + 2]; }
  const m = Math.hypot(mx, my, mz); if (m < 1e-12) return s;
  const u = [mx / m, my / m, mz / m];
  const dot = Math.max(-1, Math.min(1, u[0] * nhat[0] + u[1] * nhat[1] + u[2] * nhat[2]));
  if (dot > 1 - 1e-14) return s;
  const ax = u[1] * nhat[2] - u[2] * nhat[1], ay = u[2] * nhat[0] - u[0] * nhat[2], az = u[0] * nhat[1] - u[1] * nhat[0];
  const an = Math.hypot(ax, ay, az); if (an < 1e-14) return s;
  return applyMat(s, rotMat([ax / an, ay / an, az / an], Math.acos(dot)));
}
function B2(g, s) {
  let num = 0, den = 0;
  for (let i = 0; i < g.N; i++) {
    const a = s[3 * i], b = s[3 * i + 1], c = s[3 * i + 2];
    const p1 = a * E1[0] + b * E1[1] + c * E1[2];
    const p2 = a * E2[0] + b * E2[1] + c * E2[2];
    num += p1 * p1 - p2 * p2; den += a * a + b * b + c * c;
  }
  return den > 0 ? num / den : 0;
}
function stream(seed) { let x = seed | 0; if (x === 0) x = 0x9e3779b9; return () => { x ^= x << 13; x ^= x >>> 17; x ^= x << 5; return (x >>> 0) / 4294967296; }; }
function runTrace(L, lambda, T, sweeps, seed) {
  const g = build(L), rnd = stream(seed);
  let s = new Float64Array(3 * g.N);
  for (let i = 0; i < g.N; i++) {
    const a = NHAT[0] + 0.6 * (rnd() - 0.5), b = NHAT[1] + 0.6 * (rnd() - 0.5), c = NHAT[2] + 0.6 * (rnd() - 0.5);
    const n = Math.hypot(a, b, c) || 1;
    s[3 * i] = a / n; s[3 * i + 1] = b / n; s[3 * i + 2] = c / n;
  }
  const beta = 1 / T, delta = 0.6;
  const trace = [];
  for (let sw = 1; sw <= sweeps; sw++) {
    for (let i = 0; i < g.N; i++) {
      let hx = 0, hy = 0, hz = 0;
      for (const j of g.nbrs[i]) { hx += s[3 * j]; hy += s[3 * j + 1]; hz += s[3 * j + 2]; }
      const ax = s[3 * i], ay = s[3 * i + 1], az = s[3 * i + 2];
      let px = ax + delta * (2 * rnd() - 1), py = ay + delta * (2 * rnd() - 1), pz = az + delta * (2 * rnd() - 1);
      const pn = Math.hypot(px, py, pz) || 1; px /= pn; py /= pn; pz /= pn;
      const dE = -(px - ax) * hx - (py - ay) * hy - (pz - az) * hz
        - lambda * ((px ** 4 + py ** 4 + pz ** 4) - (ax ** 4 + ay ** 4 + az ** 4));
      if (dE <= 0 || rnd() < Math.exp(-beta * dE)) { s[3 * i] = px; s[3 * i + 1] = py; s[3 * i + 2] = pz; }
    }
    s = gaugeFix(g, s, NHAT);
    if (sw % 500 === 0) trace.push({ sw, b2: B2(g, s) });
  }
  return trace;
}

console.log('============================================================');
console.log('PREDICTIONS');
console.log('============================================================');
console.log('   Q1 B2 falls as the number of sweeps grows (the first run was not equilibrated)');
console.log('   Q2 it approaches round 95\'s 1.0240e-1');
console.log('   Q3 if Q1 fails, the difference is not equilibration');
console.log('');
console.log('   B2 against the sweep number, lambda = 2, T = 0.7, L = 12 (round 95\'s own parameters):');
console.log('   sweeps      seed 11      seed 22      seed 33');
{
  const tr = [11, 22, 33].map((sd) => runTrace(12, 2, 0.7, 12000, sd));
  const n = tr[0].length;
  for (let k = 0; k < n; k++) {
    console.log(`   ${String(tr[0][k].sw).padStart(6)}      ${tr[0][k].b2.toFixed(4)}       ${tr[1][k].b2.toFixed(4)}       ${tr[2][k].b2.toFixed(4)}`);
  }
  console.log('');
  const last = [tr[0][n - 1].b2, tr[1][n - 1].b2, tr[2][n - 1].b2];
  const mean = last.reduce((a, b) => a + b, 0) / last.length;
  const half = [tr[0][Math.floor(n / 2)].b2, tr[1][Math.floor(n / 2)].b2, tr[2][Math.floor(n / 2)].b2];
  const meanHalf = half.reduce((a, b) => a + b, 0) / half.length;
  console.log(`   at 6000 sweeps the mean is ${meanHalf.toFixed(4)};  at 12000 sweeps it is ${mean.toFixed(4)}   (round 95: 0.1024)`);
  console.log(`   => Q1 ${Math.abs(mean - meanHalf) > 0.005 ? 'CONFIRMED: it is still drifting, so the first run was indeed unequilibrated' : 'NOT confirmed: it looks settled'}`);
  console.log(`   => Q2 ${Math.abs(mean - 0.1024) < 0.02 ? 'CONFIRMED: it settles at round 95\'s value' : 'NOT confirmed: it settles elsewhere, so the difference lies elsewhere'}`);
  console.log('');
  console.log('   ON THE DIRECTION, stated carefully because I got it wrong in the first draft of this note: the');
  console.log('   initial configuration is aligned along n_hat, so its transverse anisotropy starts near zero, and');
  console.log('   one would expect B2 to GROW in -- which would make a short run UNDERestimate, not overestimate.');
  console.log('   My run came out ABOVE round 95, so either the equilibrium value really is larger and round 95');
  console.log('   was the unequilibrated one, or the difference lies somewhere else entirely.  The trace above is');
  console.log('   what decides it; the reasoning is not allowed to.');
}
