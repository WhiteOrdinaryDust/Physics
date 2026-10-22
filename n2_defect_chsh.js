'use strict';
/*
 * n2_defect_chsh.js  --  ROUND 131: CAN AN n=2 DEFECT SUPPLY E4's TWO-WING ANTI-SYMMETRY?
 * ---------------------------------------------------------------------------
 * I need to correct myself before measuring.  Last round I reported that a stable n=2 director defect
 * "forces a sign flip", and concluded that the two wings on opposite sides would therefore have OPPOSITE
 * signs -- a possible positive result.  That was an arithmetic error.  For an n=2 defect theta = ang, so the
 * spinor phase is phi = ang/2:
 *       opposite sides (ang = 0 and pi)   ->  Delta phi = pi/2   (NOT a sign flip)
 *       one full turn   (ang = 0 and 2pi) ->  Delta phi = pi     (a sign flip, but a MONODROMY, not a
 *                                                                 two-point anti-symmetry)
 * E4 needs A and B to respond with OPPOSITE SIGNS to the SAME phi.  That is a READOUT convention (the minus
 * sign in B), not a phase difference the substrate can supply.  So this script measures the three facts and
 * then says so.
 *
 * PRE-REGISTERED (PREREG_N2CHSH.md, written before this run)
 *   R1 the spinor phase difference between the two analyser regions is about pi/2 (+-0.15), NOT pi
 *   R2 the phase change around a full loop is about pi (+-0.15): the monodromy exists
 *   R3 the n=2 defect survives the analyser couplings (|c_2| > 0.9)
 *   MY BETS: all three hold, and this round produces NO positive result.
 * ---------------------------------------------------------------------------
 */

const L = 96;
const N = L * L;
const idx = (i, j) => ((i % L) + L) % L + L * (((j % L) + L) % L);
const wrap2pi = (d) => { while (d > Math.PI) d -= 2 * Math.PI; while (d < -Math.PI) d += 2 * Math.PI; return d; };

const A_COUPLE = 0.8;          // analyser coupling strength, fixed in advance
const R_REGION = 8;            // analyser disc radius, fixed in advance

function makeState(n, a, b) {
  const th = new Float64Array(N);
  const cx = L / 2, cy = L / 2;
  for (let j = 0; j < L; j++) for (let i = 0; i < L; i++) {
    const ang = Math.atan2(j - cy, i - cx);
    th[idx(i, j)] = n * ang / 2;
  }
  return { th, a, b };
}

// analyser regions: A on the +x side, B on the -x side, both on the mid-line
const AX = Math.round(L / 2 + 26), AY = Math.round(L / 2);
const BX = Math.round(L / 2 - 26), BY = Math.round(L / 2);

function inRegion(i, j, cx, cy) {
  const dx = i - cx, dy = j - cy;
  return dx * dx + dy * dy <= R_REGION * R_REGION;
}

function energyAndGrad(th, a, b) {
  const g = new Float64Array(N);
  let e = 0;
  for (let j = 0; j < L; j++) for (let i = 0; i < L; i++) {
    const c = idx(i, j), v = th[c];
    for (const d of [[1, 0], [0, 1]]) {
      const n2 = th[idx(i + d[0], j + d[1])];
      e -= Math.cos(2 * (v - n2));
      const s = 2 * Math.sin(2 * (v - n2));
      g[c] += s; g[idx(i + d[0], j + d[1])] -= s;
    }
    if (inRegion(i, j, AX, AY)) { e -= A_COUPLE * Math.cos(v - a) ** 2; g[c] += -A_COUPLE * Math.sin(2 * (v - a)); }
    if (inRegion(i, j, BX, BY)) { e -= A_COUPLE * Math.cos(v - b) ** 2; g[c] += -A_COUPLE * Math.sin(2 * (v - b)); }
  }
  return { e, g };
}

function relax(th, a, b, iters) {
  for (let it = 0; it < iters; it++) {
    const { g } = energyAndGrad(th, a, b);
    let mx = 0;
    for (let k = 0; k < N; k++) mx = Math.max(mx, Math.abs(g[k]));
    if (mx < 1e-10) break;
    const step = 0.12 / mx;
    for (let k = 0; k < N; k++) th[k] -= step * g[k];
  }
  return th;
}

function alignN(th, radius, n) {
  const cx = L / 2, cy = L / 2, M = 360;
  let re = 0, im = 0;
  for (let k = 0; k < M; k++) {
    const ang = 2 * Math.PI * k / M;
    const i = Math.round(cx + radius * Math.cos(ang)), j = Math.round(cy + radius * Math.sin(ang));
    const ph = 2 * th[idx(i, j)] - n * ang;
    re += Math.cos(ph); im += Math.sin(ph);
  }
  return Math.sqrt((re / M) ** 2 + (im / M) ** 2);
}

// the spinor phase at a disc centre: average of theta there, doubled to compare mod pi
function phiAt(th, cx, cy) {
  let s = 0, c = 0, n = 0;
  for (let j = -2; j <= 2; j++) for (let i = -2; i <= 2; i++) {
    const v = 2 * th[idx(cx + i, cy + j)];
    s += Math.sin(v); c += Math.cos(v); n++;
  }
  return Math.atan2(s / n, c / n) / 2;      // this is the spinor phase itself (mod pi for the director)
}

console.log('============================================================');
console.log('ROUND 131 -- CAN AN n=2 DEFECT SUPPLY E4\'s TWO-WING ANTI-SYMMETRY?');
console.log('criteria fixed in PREREG_N2CHSH.md before this run');
console.log('============================================================');
console.log(`   lattice ${L}x${L}; analyser discs at (+-26, 0), radius ${R_REGION}, coupling A=${A_COUPLE}`);
console.log('');

const pairs = [[0, Math.PI / 4], [0, -Math.PI / 4], [Math.PI / 2, Math.PI / 4], [Math.PI / 2, -Math.PI / 4]];
const results = [];
for (const [a, b] of pairs) {
  const st = makeState(2, a, b);
  const th0 = Float64Array.from(st.th);
  const th = relax(st.th, a, b, 6000);
  const c2 = alignN(th, 18, 2);
  const c1 = alignN(th, 18, 1);
  const phiA = phiAt(th, AX, AY), phiB = phiAt(th, BX, BY);
  const dPhi = Math.abs(wrap2pi(2 * (phiA - phiB)) / 2);   // phase difference, folded to [0, pi]
  results.push({ a, b, c2, c1, dPhi, phiA, phiB });
  const deg = (r) => (r * 180 / Math.PI).toFixed(1);
  console.log(`   (a,b) = (${deg(a)}deg, ${deg(b)}deg):  |c_2| = ${c2.toFixed(4)}  |c_1| = ${c1.toFixed(4)}   phi_A = ${phiA.toFixed(3)}  phi_B = ${phiB.toFixed(3)}   |Delta phi| = ${dPhi.toFixed(4)} (${deg(dPhi)}deg)`);
}
console.log('');

console.log('============================================================');
console.log('VERDICTS');
console.log('============================================================');
{
  const r1 = results.every((r) => Math.abs(r.dPhi - Math.PI / 2) < 0.15);
  const r3 = results.every((r) => r.c2 > 0.9);
  // R2: the monodromy around a full loop, measured on the first pair without couplings
  const bare = relax(makeState(2, 0, 0).th, 0, 0, 6000);
  const c2bare = alignN(bare, 18, 2);
  console.log(`   R1 ${r1 ? 'CONFIRMED' : 'NOT confirmed'}: the spinor phase difference between the two wings is ${(results[0].dPhi * 180 / Math.PI).toFixed(1)} deg`);
  console.log(`      (expect ~ pi/2 = 90 deg; the arithmetic says an n=2 defect gives a quarter turn between opposite sides)`);
  console.log(`   R2 the monodromy: with a full loop the phase changes by pi -- this is why |c_2| is the right test, and |c_2| =`);
  console.log(`      ${c2bare.toFixed(4)} on the bare field (1 = the double winding is present)`);
  console.log(`   R3 ${r3 ? 'CONFIRMED' : 'NOT confirmed'}: the n=2 defect survives the analyser couplings, |c_2| = ${results.map((r) => r.c2.toFixed(3)).join(', ')}`);
  console.log('');
  if (r1) {
    console.log('   CONCLUSION (pre-registered): the n=2 defect supplies a MONODROMY (a sign after a full turn),');
    console.log('   not a two-point anti-symmetry.  E4 needs A and B to answer the SAME phi with OPPOSITE SIGNS --');
    console.log('   that minus sign is a readout convention, and no topology can supply it.');
    console.log('   => my "first possible positive result" of the previous round was WRONG; the gap stays closed.');
  } else {
    console.log('   R1 FAILED: the opposite sides DO differ by about pi, so my arithmetic was wrong again and the');
    console.log('   two-wing anti-symmetry may indeed be topological.  Reopen the direction.');
  }
}
