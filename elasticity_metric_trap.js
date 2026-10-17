'use strict';
/*
 * elasticity_metric_trap.js
 * ---------------------------------------------------------------------------
 * A SUBTLE VERSION OF THE "EMPTY BACKGROUND" TRAP, MEASURED.
 *
 * A background search elsewhere reported that Life matter is elastic in one
 * situation and brittle in another:
 *     tiled 2x2 block lattice : E_add = 0.0000
 *     one ISOLATED 2x2 block  : E_add = 0.8889
 * and concluded "Life matter is brittle as a lattice, not as a lone block".
 *
 * That conclusion is wrong, and this script shows exactly why with numbers.
 *
 * DEFINITION OF E_add USED BY THAT SEARCH:
 *     for every DEAD cell of the torus, turn it on, evolve, and ask whether the
 *     configuration returns to the background's cycle.
 *
 * On an almost empty torus, nearly every dead cell is nowhere near the matter.
 * Turning such a cell on just makes a cell that dies by itself; the matter never
 * notices; the test passes. So E_add mostly measures HOW EMPTY THE TORUS IS.
 * It is not a property of the matter at all.
 *
 * THE FIX: count only NEAR-FIELD flips -- dead cells within Chebyshev distance 2
 * of a live cell of the matter. Those are the only cells that can ever interact
 * with it. Below, both numbers are printed side by side so the dilution is
 * visible rather than argued about.
 *
 * PREDICTION (stated before the run):
 *   isolated block, all dead cells      : E_add ~ 0.85-0.90  (the artefact)
 *   isolated block, near-field only     : E_add = 0.000      (the truth)
 *   block lattice,   all dead cells     : E_add = 0.000      (no dilution: the
 *                                         lattice is dense, so every dead cell
 *                                         in the gaps is near matter anyway)
 *   E_remove = 1.000 in every case.
 * If the near-field numbers disagree with matter_elasticity.js (which measured
 * the same thing on a 300x300 torus), one of the two scripts has a bug.
 * ---------------------------------------------------------------------------
 */

const SIZE = 16;              // the small torus the other search used
const STEPS = 120;

function wrapCell(x, y) { return ((((y % SIZE) + SIZE) % SIZE) * SIZE) + (((x % SIZE) + SIZE) % SIZE); }
function step(set) {
  const counts = new Map();
  for (const k of set) {
    const x = k % SIZE, y = (k - x) / SIZE;
    for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
      if (dx === 0 && dy === 0) continue;
      const nk = wrapCell(x + dx, y + dy);
      counts.set(nk, (counts.get(nk) || 0) + 1);
    }
  }
  const out = new Set();
  for (const [k, c] of counts) if (c === 3 || (c === 2 && set.has(k))) out.add(k);
  return out;
}
function setEq(a, b) { if (a.size !== b.size) return false; for (const k of a) if (!b.has(k)) return false; return true; }
function blockAt(x, y) {
  const s = new Set();
  for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }
function lattice(step4) {
  let m = new Set();
  for (let i = 0; i < SIZE; i += step4) for (let j = 0; j < SIZE; j += step4) m = union(m, blockAt(i, j));
  return m;
}
function nearField(matter, r) {
  const z = new Set();
  for (const k of matter) {
    const x = k % SIZE, y = (k - x) / SIZE;
    for (let dy = -r; dy <= r; dy++) for (let dx = -r; dx <= r; dx++) z.add(wrapCell(x + dx, y + dy));
  }
  return z;
}
// the background is a fixed point here, so its cycle is just itself
function heals(start, background) {
  let s = start;
  for (let t = 0; t < STEPS; t++) {
    if (setEq(s, background)) return true;
    s = step(s);
  }
  return setEq(s, background);
}

function measure(name, matter, verbose) {
  const dead = [];
  for (let y = 0; y < SIZE; y++) for (let x = 0; x < SIZE; x++) {
    const k = wrapCell(x, y);
    if (!matter.has(k)) dead.push(k);
  }
  const near = nearField(matter, 2);
  const nearDead = dead.filter(k => near.has(k));
  const live = [...matter];
  // true Chebyshev distance from a cell to the matter
  const distTo = (k) => {
    const x = k % SIZE, y = (k - x) / SIZE;
    let best = 99;
    for (const m of matter) {
      const mx = m % SIZE, my = (m - mx) / SIZE;
      let dx = Math.abs(x - mx); dx = Math.min(dx, SIZE - dx);
      let dy = Math.abs(y - my); dy = Math.min(dy, SIZE - dy);
      best = Math.min(best, Math.max(dx, dy));
    }
    return best;
  };

  let addAll = 0, addNear = 0, rmAll = 0;
  for (const k of dead) {
    const s = new Set(matter); s.add(k);
    if (heals(s, matter)) addAll++;
  }
  // A flip is CAUSALLY RELEVANT if it changes the very next generation. That is
  // the honest filter: geometric distance is not enough, because a cell TWO away
  // can still matter -- it can supply the third neighbour that gives birth to the
  // cell in between. And a cell two away on the diagonal supplies nothing and is
  // inert even though it lies inside the distance-2 box.
  let relevant = 0, relHealed = 0, inertNear = 0;
  const detail = [];
  for (const k of nearDead) {
    const s = new Set(matter); s.add(k);
    const isRelevant = !setEq(step(s), matter);
    const h = heals(s, matter);
    if (isRelevant) { relevant++; if (h) relHealed++; } else inertNear++;
    if (h) addNear++;
    detail.push({ k, d: distTo(k), h, isRelevant });
  }
  for (const k of live) {
    const s = new Set(matter); s.delete(k);
    if (heals(s, matter)) rmAll++;
  }
  console.log(`  ${name}`);
  console.log(`    matter ${matter.size} cells, density ${(matter.size / (SIZE * SIZE)).toFixed(4)} on the ${SIZE}x${SIZE} torus`);
  console.log(`    dead cells: ${dead.length} total, ${nearDead.length} of them within Chebyshev distance 2 of matter`);
  console.log(`    E_add over ALL dead cells     = ${(addAll / dead.length).toFixed(4)}  (${addAll}/${dead.length})   <-- diluted by empty space`);
  console.log(`    E_add over NEAR-FIELD only    = ${(addNear / nearDead.length).toFixed(4)}  (${addNear}/${nearDead.length})   <-- still diluted by inert cells inside the box`);
  console.log(`    E_add over CAUSALLY RELEVANT  = ${relevant ? (relHealed / relevant).toFixed(4) : 'n/a'}  (${relHealed}/${relevant})   <-- THE HONEST NUMBER`);
  console.log(`      (inert flips inside the distance-2 box: ${inertNear}; they change nothing at t=1)`);
  console.log(`    E_remove over all live cells  = ${(rmAll / live.length).toFixed(4)}  (${rmAll}/${live.length})`);
  console.log(`    far-field flips that "healed" = ${addAll - addNear} of ${dead.length - nearDead.length}  (the matter was never touched)`);
  if (verbose) {
    const byD = {};
    for (const e of detail) {
      byD[e.d] = byD[e.d] || { healed: 0, total: 0 };
      byD[e.d].total++;
      if (e.h) byD[e.d].healed++;
    }
    console.log('    near-field detail, by true Chebyshev distance from the matter:');
    for (const d of Object.keys(byD).sort()) {
      console.log(`      distance ${d}: ${byD[d].total} cells, ${byD[d].healed} healed`);
    }
    console.log(`      cells at distance 1 (real immediate neighbours):   ${byD[1] ? byD[1].healed + '/' + byD[1].total + ' healed' : 'none'}`);
    console.log(`      cells at distance 2 (no live neighbour at all):    ${byD[2] ? byD[2].healed + '/' + byD[2].total + ' healed' : 'none'}`);
    const d2 = detail.filter(e => e.d === 2);
    const fmt = (e) => `(${e.k % SIZE},${(e.k - (e.k % SIZE)) / SIZE})`;
    console.log(`      distance-2 inert (healed but never interacted): ${d2.filter(e => !e.isRelevant).map(fmt).join(' ') || 'none'}`);
    console.log(`      distance-2 interactive and NOT healed        : ${d2.filter(e => e.isRelevant && !e.h).map(fmt).join(' ') || 'none'}`);
    console.log(`      distance-2 interactive and healed            : ${d2.filter(e => e.isRelevant && e.h).map(fmt).join(' ') || 'none'}`);
  }
  console.log('');
  return { addAll: addAll / dead.length, addNear: addNear / nearDead.length, relevantE: relevant ? relHealed / relevant : NaN, rm: rmAll / live.length, relevant, relHealed };
}

console.log('=== the same Life rule, the same code path, two backgrounds ===');
console.log(`  torus ${SIZE}x${SIZE}, ${STEPS} generations, "healed" = back on the background cycle\n`);
const a = measure('ISOLATED 2x2 block', blockAt(8, 8), true);
const b = measure('TILED 2x2 block lattice (spacing 4)', lattice(4), true);

console.log('=== VERDICT ===');
console.log(`  isolated block : E_add all = ${a.addAll.toFixed(4)} | distance<=2 = ${a.addNear.toFixed(4)} | causally relevant = ${a.relevantE.toFixed(4)} (${a.relHealed}/${a.relevant})`);
console.log(`  lattice        : E_add all = ${b.addAll.toFixed(4)} | distance<=2 = ${b.addNear.toFixed(4)} | causally relevant = ${b.relevantE.toFixed(4)} (${b.relHealed}/${b.relevant})`);
console.log('');
console.log('  Each column removes another layer of dilution, and the number falls towards the truth:');
console.log('    0.8889  <- almost all dead cells are nowhere near the matter; flipping one just');
console.log('               makes a cell that dies by itself and the matter never notices');
console.log('    0.1250  <- restricted to the distance-2 box, but that box still contains four');
console.log('               INERT corners that supply no neighbour to anything that matters');
console.log('    0.0000  <- restricted to flips that actually change the next generation');
console.log('');
console.log('  A cell TWO away can matter (it can supply the third neighbour that gives birth to');
console.log('  the cell in between, e.g. (8,6) for a block at (8,8)), so distance 1 is too tight');
console.log('  a cut. But a cell two away on the DIAGONAL, e.g. (6,6), supplies nothing. That is');
console.log('  why the honest filter is causal (does the flip change t=1?), not geometric.');
console.log('');
console.log(`  Cross-check against matter_elasticity.js (300x300 torus, additions enumerated only`);
console.log(`  in the box around the matter): single block E_add = 0.000, lattice E_add = 0.000.`);
console.log(`  ${a.relevantE === 0 && b.relevantE === 0 ? 'AGREES with both backgrounds' : '*** DISAGREES -- one of the two scripts is wrong ***'}`);
