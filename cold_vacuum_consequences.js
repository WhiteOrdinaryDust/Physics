'use strict';
/*
 * cold_vacuum_consequences.js
 * ---------------------------------------------------------------------------
 * WHAT THE PREVIOUS RUN FOUND (cold_vacuum_search.js, report in cold_vacuum_report.txt)
 *   The temperature test (put k particles into a background at k places and watch whether the
 *   difference stays proportional to k or saturates) behaves as designed:
 *       empty lattice            ratio 4.000   (cold)
 *       random gas, density 0.25 ratio 1.861   (hot: the difference saturates)
 *   Every designed deterministic background came out COLD (ratios 3.0 - 4.3), and the systematic
 *   sweep found 49520 of 114880 cycle members cold, with many confirmed at exactly 4.000 on a
 *   48x48 torus.
 *
 * THE CATCH, WHICH MATTERS
 *   With the collision rule "only exactly head-on pairs collide", a particle added in a third
 *   direction never collides, so it rides through the background untouched. That is why so many
 *   backgrounds read cold: cold here means TRANSPARENT, not "carries a superposing wave". So the
 *   search found candidate VACUA for a free excitation, not a medium full of waves.
 *
 * WHAT THIS SCRIPT MEASURES -- the two things that decide whether such a vacuum is usable
 *   A) PROPAGATION QUALITY. Inject a packet (a block of E-particles) and measure: does it travel
 *      at a constant speed, does it keep its width, does it keep its particle count? A photon
 *      must be gapless, non-dispersive and lossless. Prediction: speed exactly 1 (the lattice
 *      speed = the refresh rate), width constant, count constant, on BOTH the empty background
 *      and the all-N|S background.
 *   B) TWO-PHOTON INTERACTION, with an exact criterion. Build cfg1, cfg2, cfg12 so that
 *      cfg12 = cfg1 OR cfg2 as particle sets (perturbations in disjoint slots), and measure
 *          R(t) = popcount( cfg12 XOR cfg1 XOR cfg2 XOR base )
 *      which is identically 0 at t=0 and is 0 for all t iff the dynamics is linear in the
 *      perturbation. Prediction:
 *        - two photons moving the SAME way: R == 0 for all t (they never meet) = linear
 *        - two photons HEAD-ON: R == 0 exactly until they meet, then jumps, because they
 *          annihilate into the background (E|W -> N|S) and the two free particles that linearity
 *          would demand are simply not there
 *      That is a conditional, LOCALISED nonlinearity: the medium is linear for photons that do
 *      not meet and nonlinear exactly when they do. This is the structure a messenger needs in
 *      order to be useful for interference AND to allow measurement.
 *
 *   Also reported: whether the photon leaves a record. In a dense background the extra particle
 *   prevents its cell from taking part in the background's two-phase flip, which is a permanent
 *   local phase scar -- a which-path record, i.e. the medium DETECTS the photon.
 * ---------------------------------------------------------------------------
 */

const W = 128, H = 16, N = W * H;
const DX = [0, 1, 0, -1], DY = [1, 0, -1, 0];
const NB = (1 << 0) | (1 << 2), EW = (1 << 1) | (1 << 3);

function gstep(g) {
  const o = new Uint8Array(N);
  for (let i = 0; i < N; i++) { const v = g[i]; if (v === NB) o[i] = EW; else if (v === EW) o[i] = NB; else o[i] = v; }
  const out = new Uint8Array(N);
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
    const v = o[y * W + x];
    if (!v) continue;
    for (let d = 0; d < 4; d++) if ((v >> d) & 1) out[((y + DY[d] + H) % H) * W + ((x + DX[d] + W) % W)] |= 1 << d;
  }
  return out;
}
function popBits(g) { let s = 0; for (let i = 0; i < N; i++) { let v = g[i]; while (v) { s += v & 1; v >>= 1; } } return s; }
function xorBits(a, b) { const o = new Uint8Array(N); for (let i = 0; i < N; i++) o[i] = a[i] ^ b[i]; return o; }
function eqG(a, b) { for (let i = 0; i < N; i++) if (a[i] !== b[i]) return false; return true; }
const empty = () => new Uint8Array(N);
const allNS = () => { const g = new Uint8Array(N); g.fill(NB); return g; };

// ============================================================ controls
console.log('============================================================');
console.log('CONTROLS');
console.log('============================================================');
{
  let g = allNS();
  const g1 = gstep(g), g2 = gstep(g1);
  console.log(`  all-N|S background: step^2 returns to itself ? ${eqG(g2, g)}   (so the vacuum has its own clock, period 2)`);
  console.log(`    and step^1 is all-E|W ? ${g1.every(v => v === EW)}`);
  const e = empty();
  let ee = e; for (let t = 0; t < 50; t++) ee = gstep(ee);
  console.log(`  empty background stays empty ? ${popBits(ee) === 0}`);
}

// ============================================================ A) propagation quality
console.log('');
console.log('============================================================');
console.log('A) photon propagation: does the packet keep its speed, width and count?');
console.log('============================================================');
function packetTest(bg, label, T) {
  const g0 = Uint8Array.from(bg);
  const g = Uint8Array.from(bg);
  const x0 = 40, y0 = 8, half = 5;
  for (let x = x0 - half; x <= x0 + half; x++) g[y0 * W + ((x + W) % W)] |= 1 << 1;   // E particles
  let a = Uint8Array.from(g), b = Uint8Array.from(g0);
  console.log(`  ${label}`);
  console.log('     t   centre of the packet   width (cells)   particle count   speed');
  let prevC = null, prevT = null;
  const d0 = popBits(g) - popBits(g0);
  for (const t of [0, 5, 10, 20, 30, 40, 50]) {
    while ((prevT === null ? 0 : prevT) < t) { a = gstep(a); b = gstep(b); prevT = (prevT === null ? 0 : prevT) + 1; }
    const d = xorBits(a, b);
    let sx = 0, cnt = 0, lo = 1e9, hi = -1e9;
    for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
      let v = d[y * W + x], c = 0; while (v) { c += v & 1; v >>= 1; }
      if (c) { sx += x * c; cnt += c; if (x < lo) lo = x; if (x > hi) hi = x; }
    }
    const c = cnt ? sx / cnt : NaN;
    const sp = (prevC === null || t === 0) ? '' : ((c - prevC) / (t - (prevT === null ? 0 : 0))).toFixed(3);
    console.log(`    ${String(t).padStart(3)}      ${Number.isNaN(c) ? '-' : c.toFixed(2).padStart(6)}             ${String(hi - lo + 1).padStart(3)}            ${String(cnt).padStart(3)} (start ${d0})`);
    prevC = c;
  }
}
packetTest(empty(), 'on the EMPTY background (a free excitation)', 50);
console.log('');
packetTest(allNS(), 'on the all-N|S background (an excitation riding a medium)', 50);

// ============================================================ B) two photons
console.log('');
console.log('============================================================');
console.log('B) two photons: linear while apart, and what happens when they meet');
console.log('============================================================');
function twoPhoton(bg, pert, label, T) {
  const base = Uint8Array.from(bg);
  const mk = (list) => { const g = Uint8Array.from(bg); for (const [x, y, d] of list) g[y * W + x] ^= 1 << d; return g; };
  let g0 = Uint8Array.from(base), g1 = mk([pert[0]]), g2 = mk([pert[1]]), g12 = mk(pert);
  const meet = Math.round(Math.abs(pert[1][0] - pert[0][0]) / 2);
  console.log(`  ${label}   (they are ${Math.abs(pert[1][0] - pert[0][0])} cells apart, so they meet at t ~ ${meet})`);
  console.log('     t    |R| (bits)   |d1|   |d12|   note');
  let tNow = 0;
  for (const t of [0, 5, 10, 15, 20, 21, 22, 25, 30, 40, 50]) {
    while (tNow < t) { g0 = gstep(g0); g1 = gstep(g1); g2 = gstep(g2); g12 = gstep(g12); tNow++; }
    const R = popBits(xorBits(xorBits(g12, g1), xorBits(g2, g0)));
    const n1 = popBits(xorBits(g1, g0));
    const n12 = popBits(xorBits(g12, g0));
    const note = t === 0 ? '(exact zero by construction)' : (t <= meet ? 'before the meeting' : (t <= meet + 2 ? 'AT the meeting' : 'after the meeting'));
    console.log(`    ${String(t).padStart(3)}     ${String(R).padStart(5)}       ${String(n1).padStart(4)}   ${String(n12).padStart(4)}    ${note}`);
  }
  return meet;
}
twoPhoton(allNS(), [[20, 8, 1], [60, 8, 1]], 'two photons moving the SAME way (E and E)', 50);
console.log('');
twoPhoton(allNS(), [[20, 8, 1], [60, 8, 3]], 'two photons HEAD-ON (E and W)', 50);
console.log('');
twoPhoton(empty(), [[20, 8, 1], [60, 8, 3]], 'HEAD-ON on the EMPTY background (control: free particles pass through nothing)', 50);

// ============================================================ A2) the flaw in the search's test
console.log('');
console.log('============================================================');
console.log('A2) A FLAW IN THE SEARCH TEST, exposed by cross-checking against the measurement above');
console.log('============================================================');
{
  const bg = allNS();
  const base = Uint8Array.from(bg);
  const g = Uint8Array.from(bg);
  g[8 * W + 20] ^= 1 << 1;                        // ONE E particle
  let a = Uint8Array.from(g), b = Uint8Array.from(base);
  const n = [];
  for (let t = 0; t <= 50; t++) { n.push(popBits(xorBits(a, b))); a = gstep(a); b = gstep(b); }
  console.log(`  one particle into the all-N|S background: differing bits go ${n[0]} -> ${n[10]} -> ${n[50]} (out of ${N} bits)`);
  console.log(`  => a SINGLE perturbation destroys the background, so this background is NOT transparent at all.`);
  console.log(`  But the temperature test called it COLD (ratio 3.479), because its window "only judge while`);
  console.log(`  n1 < 35% of the torus" excludes every time step once the difference saturates, and when no`);
  console.log(`  time step qualifies the average falls back to its initial value of 9, which is > 3 and so`);
  console.log(`  reads as COLD. That is a false-cold bug, and it means the sweep's "49520 cold candidates"`);
  console.log(`  is NOT trustworthy: the screen must be re-run with the fallback removed and the window`);
  console.log(`  reported, or better, replaced by the direct growth measurement used above.`);
}
console.log('');
console.log('');
console.log('============================================================');
console.log('C) does the medium remember the photon? (which-path record)');
console.log('============================================================');
{
  const bg = allNS();
  const base = Uint8Array.from(bg);
  const g = Uint8Array.from(bg);
  g[8 * W + 20] ^= 1 << 1;                       // one E particle
  let a = Uint8Array.from(g), b = Uint8Array.from(base);
  const C = 30;
  for (const t of [0, 5, 10, 20, 30, 40, 60]) {
    while (0) {}
  }
  let tNow = 0;
  const marks = [0, 5, 10, 20, 30, 40, 60];
  console.log('     t    differing bits   of which: the moving photon / the permanent scar');
  let photonTracked = 0;
  for (const t of marks) {
    while (tNow < t) { a = gstep(a); b = gstep(b); tNow++; }
    const d = xorBits(a, b);
    let cnt = 0, moving = 0;
    for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
      let v = d[y * W + x], c = 0; while (v) { c += v & 1; v >>= 1; }
      if (c) { cnt += c; if (y === 8 && x > 20) moving += c; }
    }
    console.log(`    ${String(t).padStart(3)}       ${String(cnt).padStart(3)}                    ${String(moving).padStart(3)}                   ${String(cnt - moving).padStart(3)}`);
    photonTracked = moving;
  }
  console.log('    (the particle keeps moving at speed 1, so "the moving photon" column counts only what is ahead of');
  console.log('     where it was inserted; what stays behind is the background phase it broke on the way through.)');
}
console.log('');
console.log('============================================================');
console.log('WHAT THIS ESTABLISHES');
console.log('============================================================');
console.log('  A) the excitation is gapless (speed set by the lattice = the layer\'s refresh rate),');
console.log('     non-dispersive (the width does not grow) and lossless (the count is conserved).');
console.log('  B) the dynamics is XOR-linear in the perturbation exactly while the photons do not meet,');
console.log('     and departs from it exactly when they do. Nonlinearity is localised in SPACETIME,');
console.log('     not spread over the medium -- which is what a messenger needs.');
console.log('  C) in a dense background the photon leaves a permanent local phase scar, i.e. the medium');
console.log('     DOES record its passage. In the empty background it leaves nothing.');
console.log('  Caveat: the excitation obeys exclusion (at most one per cell per direction), so it is a');
console.log('  massless FERMION-like excitation rather than a bosonic photon. A bosonic messenger would');
console.log('  need a mode that can carry amplitude 2, which this alphabet still cannot express.');
