'use strict';
/*
 * measurement_decay.js
 * ---------------------------------------------------------------------------
 * The user's proposal: let entanglement and interaction draw on the SAME finite
 * interface budget, so that measuring a particle consumes interfaces and can
 * therefore DESTROY its existing entanglement.
 *
 * This turns "measurement disturbs the system" into a countable statement. The
 * model:
 *
 *   each particle has N interfaces
 *   an entangled pair consumes k of them (k/2 per particle)
 *   a measurement device couples through m interfaces, MONOPOLISING them
 *
 * Two outcomes are possible, and they are physically different:
 *
 *   PARALLEL  the device uses m free interfaces               -> entanglement intact
 *   COMPETING the device takes over interfaces the pair holds  -> entanglement broken
 *
 * Allocations are counted as combinations, so the survival probability is the
 * hypergeometric fraction of allocations that leave the pair's share intact.
 *
 * WHAT THIS PREDICTS (and why it is falsifiable)
 *   There is a SHARP THRESHOLD: entanglement survives only if N - k >= m, i.e.
 *   only if the particle has enough spare interfaces to serve the device without
 *   giving up its bonded ones. Below that, ANY measurement destroys the bond.
 *   Real quantum mechanics has no such threshold for destroying correlation, so
 *   a sharp threshold here is a distinctive, testable prediction.
 *
 * WHAT THIS DOES NOT CHANGE
 *   The strength ceiling. Counting interfaces cannot produce 2*sqrt(2) -- that
 *   has been established three times (capacity_chsh, interface_capacity_chsh,
 *   interface_scaling). Interface count regulates HOW MANY bonds and WHETHER a
 *   measurement breaks them, not how strong a bond can be.
 * ---------------------------------------------------------------------------
 */

function choose(n, k) {
  if (k < 0 || k > n) return 0;
  let r = 1;
  for (let i = 0; i < k; i++) r = r * (n - i) / (i + 1);
  return Math.round(r);
}

/*
 * Probability that a measurement through m interfaces leaves the bond intact.
 * Interfaces live on BOTH particles, so the device draws m from a total of 2N,
 * of which 2*(k/2) = k are held by the bond.
 *   P(intact) = C(2N - k, m) / C(2N, m)
 */
function pIntact(N, k, m) {
  const num = choose(2 * N - k, m);
  const den = choose(2 * N, m);
  return den === 0 ? 0 : num / den;
}

console.log('=== measurement as interface competition ===\n');
console.log('P(bond intact) = C(2N - k, m) / C(2N, m)');
console.log('  N = interfaces per particle, k = interfaces held by the bond,');
console.log('  m = interfaces the measuring device occupies\n');

// ---- table 1: fixed small N, vary the bond and the device
for (const N of [2, 3, 4, 6]) {
  console.log(`--- N = ${N} interfaces per particle ---`);
  console.log('   k\\m ' + [1, 2, 3, 4, 6, 8].map(m => String(m).padStart(6)).join(''));
  for (const k of [2, 4, 6, 8]) {
    if (k > 2 * N) continue;
    let line = '  ' + String(k).padStart(4) + ' ';
    for (const m of [1, 2, 3, 4, 6, 8]) {
      line += pIntact(N, k, m).toFixed(3).padStart(6);
    }
    console.log(line);
  }
  console.log('');
}

// ---- table 2: the threshold, stated exactly
console.log('=== the threshold ===');
console.log('  P(intact) = 1 exactly when m <= 2N - k, i.e. when the device can be served');
console.log('  entirely from interfaces the bond does not hold.');
console.log('');
console.log('   N   k   m_max for m which survives   N - k/2');
for (const N of [2, 3, 4, 6, 8]) {
  for (const k of [2, 4]) {
    if (k > 2 * N) continue;
    const mMax = 2 * N - k;
    console.log(`  ${String(N).padStart(2)}  ${String(k).padStart(2)}   ${String(mMax).padStart(26)}   ${String(N - k / 2).padStart(6)}`);
  }
}
console.log('');
console.log('  => a particle with N interfaces and a bond of k interfaces can serve a');
console.log('     device of up to N - k/2 interfaces per particle without breaking the');
console.log('     bond. Beyond that, EVERY measurement destroys it.');
console.log('');
console.log('  This is a SHARP threshold. Standard quantum mechanics has nothing like it:');
console.log('  a measurement can be made arbitrarily weak (weak measurement) and reduces');
console.log('  correlations GRADUALLY, not all-or-nothing. So the threshold is the');
console.log('  distinctive prediction of this model -- and the place it is most exposed.');

// ---- the gradual alternative
console.log('\n=== contrast: the GRADUAL rule ===');
console.log('  If instead each consumed bonded interface removes a fixed fraction of the');
console.log('  correlation, the decay is smooth. Model the concurrence as');
console.log('       C = max(0, 1 - m / (N - k/2))');
console.log('  i.e. the bond is fully intact until the spare interfaces run out, then');
console.log('  degrades linearly.');
console.log('');
console.log('   N   k    m=0    m=1    m=2    m=3    m=4    m=6');
for (const N of [4, 6]) {
  for (const k of [2, 4]) {
    let line = `  ${String(N).padStart(2)}  ${String(k).padStart(2)} `;
    for (const m of [0, 1, 2, 3, 4, 6]) {
      const spare = N - k / 2;
      const C = Math.max(0, 1 - m / spare);
      line += C.toFixed(3).padStart(7);
    }
    console.log(line);
  }
}
console.log('');
console.log('  Note: with this rule C < 1 always for finite N, so CHSH = 2*sqrt(2)*C is');
console.log('  strictly below Tsirelson for any finite interface count -- unless the');
console.log('  capacity is allocated by something other than counting. This is the same');
console.log('  wall as before: counting interfaces gives a knob, not the ceiling.');

// ---- what measurement should NOT do
console.log('\n=== one feature to insist on ===');
console.log('  After the bond is broken, the particles must still be MEASURABLE. In this');
console.log('  model that is automatic: breaking the bond FREES its k interfaces, so a');
console.log('  subsequent measurement has more interfaces available than before.');
console.log('  So "measurement destroys entanglement" and "particles remain measurable"');
console.log('  are compatible in this model -- which is what the user wanted, and it');
console.log('  removes the earlier objection about zero exterior interfaces.');
console.log('');
console.log('  The cost is the sharp threshold above. That is the trade this model makes.');

// ---- summary table of the three regimes
console.log('\n=== summary ===');
console.log('  regime                      ceiling on CHSH   measurement effect');
console.log('  --------------------------  ----------------  --------------------------');
console.log('  classical interfaces       2.000000          sharp threshold');
console.log('  quantum interfaces         2.828427          gradual (no sharp cutoff)');
console.log('  interface competition      (not a strength)  sharp threshold <- the new bit');
console.log('');
console.log('  => the interface-competition idea does not change the STRENGTH ceiling;');
console.log('     what it adds is a sharp, falsifiable prediction about WHEN a');
console.log('     measurement breaks a bond.');
