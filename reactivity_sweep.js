'use strict';
/*
 * reactivity_sweep.js
 * ---------------------------------------------------------------------------
 * STEP 6: measure the "reactivity" of two-particle collisions.
 *
 * Previous round, one configuration (h1(6,1)) turned out to be a self-limiting
 * CASCADE: two gliders (mass 10) triggered 16-21x mass amplification and then
 * froze into a plateau of still lifes and oscillators held for 4000+ steps.
 *
 * That is exactly the shape of a chain reaction, so the natural quantity to
 * measure is its multiplication factor. Here we sweep all two-glider geometries
 * and record, per configuration:
 *
 *     A = peak mass / initial mass        ("amplification")
 *     frozen = whether the configuration ends as inert ash
 *
 * The question: is A a smooth function of the impact offset (continuous
 * reactivity), or is there a THRESHOLD -- most impacts do nothing and a few
 * erupt (bimodal, like a critical system)?
 *
 * Torus note: SIZE=600 and 2500 steps means a glider (speed 1/4) travels at most
 * 625 cells, so nothing wraps around and re-enters. That was the failure mode
 * that invalidated the earlier "settled" analysis.
 * ---------------------------------------------------------------------------
 */

const SIZE = 600;
function wrapCell(x, y) { return ((y % SIZE) + SIZE) % SIZE * SIZE + (((x % SIZE) + SIZE) % SIZE); }

function step(set) {
  const counts = new Map();
  for (const k of set) {
    const x = k % SIZE, y = (k - x) / SIZE;
    for (let dy = -1; dy <= 1; dy++) {
      for (let dx = -1; dx <= 1; dx++) {
        if (dx === 0 && dy === 0) continue;
        const nk = wrapCell(x + dx, y + dy);
        counts.set(nk, (counts.get(nk) || 0) + 1);
      }
    }
  }
  const out = new Set();
  for (const [k, c] of counts) if (c === 3 || (c === 2 && set.has(k))) out.add(k);
  return out;
}

const GLIDER0 = [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]];
function rotCells(turns) {
  return GLIDER0.map(([dx, dy]) => {
    let a = dx, b = dy;
    for (let i = 0; i < turns; i++) { const na = -b, nb = a; a = na; b = nb; }
    return [a, b];
  });
}
const HEAD_VEL = [[1, 1], [-1, 1], [-1, -1], [1, -1]];
function gliderAt(x, y, turns) {
  const s = new Set();
  for (const [dx, dy] of rotCells(turns)) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }

const STEPS = 2500;
const T_MEASURE = 2000;      // read mass after the cascade has had time to finish

/*
 * Amplification of one configuration.
 *   peak    = largest mass seen at any time (the eruption size)
 *   final   = mass at T_MEASURE (the frozen residue)
 *   A       = peak / initial
 * Configurations that never interact keep mass 10 throughout: A = 1.
 */
function amplify(heading, ox, oy) {
  const A0 = gliderAt(120, 120, 0);
  const B0 = gliderAt(120 + ox, 120 + oy, heading);
  const initMass = A0.size + B0.size;
  let s = union(A0, B0);
  if (s.size !== initMass) return { overlap: true };
  let peak = s.size, peakT = 0, changed = false;
  for (let t = 1; t <= STEPS; t++) {
    s = step(s);
    if (s.size > peak) { peak = s.size; peakT = t; }
    if (s.size !== initMass) changed = true;
    if (s.size === 0) return { initMass, peak, peakT, final: 0, A: 0, interacted: true, diedT: t };
  }
  return { initMass, peak, peakT, final: s.size, A: peak / initMass, interacted: changed };
}

console.log('=== reactivity sweep: amplification of two-glider collisions ===');
console.log(`torus ${SIZE}x${SIZE}, ${STEPS} steps, mass read at t=${T_MEASURE}`);
console.log('A = peak mass / initial mass (initial mass 10 for all)\n');

const OFFS = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12];
const rows = [];
for (const heading of [1, 2, 3]) {
  for (const ox of OFFS) {
    for (const oy of OFFS) {
      const r = amplify(heading, ox, oy);
      if (r.overlap) continue;
      rows.push({ heading, ox, oy, ...r });
    }
  }
}

// ---- amplification as a function of offset: a grid of numbers
console.log('=== amplification grid (rows oy, cols ox) ===');
for (const heading of [1, 2, 3]) {
  console.log(`\nheading turns=${heading}  B velocity (${HEAD_VEL[heading]})`);
  console.log('   oy\\ox ' + OFFS.map(o => String(o).padStart(5)).join(''));
  for (const oy of OFFS) {
    let line = '   ' + String(oy).padStart(5) + ' ';
    for (const ox of OFFS) {
      const r = rows.find(q => q.heading === heading && q.ox === ox && q.oy === oy);
      if (!r) { line += '    x'; continue; }
      const v = r.A === Infinity ? 0 : r.A;
      line += String(v.toFixed(1)).padStart(5);
    }
    console.log(line);
  }
}

// ---- distribution of A: continuous or bimodal?
console.log('\n=== distribution of amplification A ===');
const buckets = new Map();
for (const r of rows) {
  const a = r.A;
  const b = a === 0 ? '0 (annihilated)'
    : a === 1 ? '1.0 (no interaction)'
      : a < 1.5 ? '1.0-1.5 (barely)'
        : a < 3 ? '1.5-3 (mild)'
          : a < 6 ? '3-6 (moderate)'
            : a < 12 ? '6-12 (strong)'
              : a < 25 ? '12-25 (CASCADE)'
                : '>25 (huge cascade)';
  buckets.set(b, (buckets.get(b) || 0) + 1);
}
const maxA = Math.max(...rows.map(r => r.A));
console.log(`  configurations: ${rows.length}   max A = ${maxA.toFixed(1)}`);
for (const [k, v] of [...buckets.entries()]) {
  console.log(`  ${String(v).padStart(4)}x  ${k}`);
}

// ---- are the cascade configurations clustered in offset space?
console.log('\n=== where do the cascades live? (A >= 6) ===');
const big = rows.filter(r => r.A >= 6).sort((a, b) => b.A - a.A);
console.log(`  ${big.length} of ${rows.length} configurations amplify by 6x or more:`);
for (const r of big.slice(0, 30)) {
  console.log(`    h${r.heading}(${r.ox},${r.oy})  A=${r.A.toFixed(1)}  peak=${r.peak}@t${r.peakT}  final=${r.final}`);
}

// ---- threshold test: is reactivity smooth in ox for fixed heading/oy?
console.log('\n=== threshold test: A as a function of ox, at fixed oy ===');
for (const heading of [1, 2, 3]) {
  for (const oy of [0, 1, 2, 4]) {
    const line = OFFS.map(ox => {
      const r = rows.find(q => q.heading === heading && q.ox === ox && q.oy === oy);
      return r ? r.A.toFixed(1).padStart(6) : '     x';
    }).join('');
    console.log(`  h${heading} oy=${oy} :${line}`);
  }
}
console.log('  (columns are ox = ' + OFFS.join(', ') + ')');

// ---- how many interactions freeze completely?
console.log('\n=== outcome classes ===');
const dead = rows.filter(r => r.A === 0).length;
const inert = rows.filter(r => r.A === 1).length;
const finite = rows.filter(r => r.A > 1 && r.A < 6).length;
const cascade = rows.filter(r => r.A >= 6).length;
console.log(`  total annihilation        : ${dead}`);
console.log(`  no interaction            : ${inert}`);
console.log(`  bounded reaction (A<6)    : ${finite}`);
console.log(`  cascade (A>=6)            : ${cascade}`);

console.log('');
if (cascade > 0 && finite > 0 && big.length < rows.length * 0.25) {
  console.log('  => BIMODAL: most collisions are inert or bounded, a minority erupt.');
  console.log('     Reactivity is THRESHOLD-LIKE, not smoothly tunable.');
} else if (cascade > rows.length * 0.25) {
  console.log('  => Cascades are common; reactivity is not a rare threshold phenomenon here.');
} else {
  console.log('  => No cascades found at all in this sweep.');
}
