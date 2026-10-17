'use strict';
/*
 * matter_elasticity.js
 * ---------------------------------------------------------------------------
 * WHY THIS SCRIPT EXISTS
 *
 *   A real screen is a macroscopic, very stable object: you can pick it up.
 *   The walls we have been building (a row of five 2x2 blocks) are neither
 *   large nor robust, so the earlier negative results -- "any contact damages
 *   the wall", "a clean pass requires that the glider never touches anything" --
 *   might be artefacts of a bad wall.
 *
 * WHAT "STABLE" ACTUALLY MEANS
 *
 *   A macroscopic object is not robust because it is big. It is robust because
 *   it is ELASTIC: a small local disturbance (a hole, a bump, a scratch) is
 *   repaired and the object returns to the same state. Elasticity is exactly
 *   what makes a real screen a NON-detector: a particle passing by perturbs it,
 *   the perturbation dies away, and no record is left behind.
 *
 *   So the question is quantitative and computable. For a piece of matter M,
 *   define its elasticity
 *
 *       E_remove = fraction of one-cell REMOVALS that are repaired
 *       E_add    = fraction of one-cell ADDITIONS that are repaired
 *
 *   with "repaired" meaning: after 200 generations the matter occupies exactly
 *   its original cells again. A real solid has E_remove ~ E_add ~ 1.
 *
 * PREDICTIONS (stated before the run, per method rule)
 *
 *   P1  Removing one cell from a 2x2 block IS repaired, in one step: the three
 *       remaining cells each keep exactly 2 neighbours and the empty corner has
 *       exactly 3, so the block is restored exactly.  Expect E_remove = 1.00
 *       for an isolated block AND for a block inside a lattice (blocks 4 apart
 *       do not interact).
 *   P2  Adding one cell next to a block is FATAL. The added cell and the block
 *       cell beside it both see 4 neighbours and die, and the block unravels.
 *       Expect E_add = 0.00.
 *   P3  Removing a whole block (a vacancy) is PERMANENT: the surrounding blocks
 *       cannot reach into the hole. Expect E_vacancy = 0.00.
 *
 *   So the prediction is that Life matter has a ONE-SIDED restoring force: it
 *   heals holes, never bumps, and never fills a vacancy. A real solid does all
 *   three. If that comes out, then "our wall is not stable enough" is not a
 *   matter of scale (make it bigger) but of KIND: Life matter is brittle, and
 *   no amount of extra blocks will give it the property that makes a real
 *   screen transparent-but-not-a-detector.
 * ---------------------------------------------------------------------------
 */

const SIZE = 300;
const STEPS = 200;

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
function components(set) {
  const seen = new Set(), comps = [], Q = [];
  for (const start of set) {
    if (seen.has(start)) continue;
    Q.length = 0; Q.push(start); seen.add(start);
    const cells = [];
    while (Q.length) {
      const k = Q.pop(); cells.push(k);
      const x = k % SIZE, y = (k - x) / SIZE;
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        if (dx === 0 && dy === 0) continue;
        const nk = wrapCell(x + dx, y + dy);
        if (set.has(nk) && !seen.has(nk)) { seen.add(nk); Q.push(nk); }
      }
    }
    comps.push(cells);
  }
  return comps;
}
function bboxOf(cells) {
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
  }
  return { minx, miny, maxx, maxy, w: maxx - minx + 1, h: maxy - miny + 1 };
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }
function blockAt(x, y) {
  const s = new Set();
  for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function lattice(X0, Y0, nx, ny, s) {
  let m = new Set();
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) m = union(m, blockAt(X0 + i * s, Y0 + j * s));
  return m;
}

// --------------------------------------------------------------- classify
function runCase(material, live, steps, bb) {
  let s = live;
  for (let t = 0; t < steps; t++) s = step(s);
  const inBox = new Set();
  for (const k of s) {
    const x = k % SIZE, y = (k - x) / SIZE;
    if (x >= bb.minx && x <= bb.maxx && y >= bb.miny && y <= bb.maxy) inBox.add(k);
  }
  let removed = 0, added = 0;
  for (const k of material) if (!inBox.has(k)) removed++;
  for (const k of inBox) if (!material.has(k)) added++;
  const outside = s.size - inBox.size;
  let kind;
  if (removed === 0 && added === 0) kind = outside === 0 ? 'RESTORED' : 'RESTORED+EMIT';
  else if (inBox.size === 0) kind = 'GONE';
  else kind = 'DAMAGED';
  return { kind, removed, added, outside };
}

// ================================================================ controls
console.log('=== CONTROLS ===');
{
  const m = lattice(150, 150, 3, 3, 4);
  const bb = bboxOf(m);
  const r = runCase(m, m, STEPS, bb);
  console.log(`  matter alone, 200 generations: ${r.kind}  ${r.kind === 'RESTORED' ? 'PASS' : '*** FAIL ***'}`);
  // sanity: a perturbation far away must be reported as an emission, not as damage
  let alien = union(m, blockAt(150, 250));
  const r2 = runCase(m, alien, STEPS, bb);
  // a lone block is a still life, so this must come out as RESTORED+EMIT (the
  // matter is intact, and the extra block simply stayed where it was)
  console.log(`  far-away extra block: ${r2.kind} (expected RESTORED+EMIT)`);
  if (r.kind !== 'RESTORED') { console.log('  *** control failed, aborting ***'); process.exit(1); }
}
console.log('');

// ================================================================ materials
const MATERIALS = [
  ['single block', lattice(150, 150, 1, 1, 4)],
  ['row of 3 blocks', lattice(146, 150, 3, 1, 4)],
  ['slab 3x3 blocks', lattice(146, 146, 3, 3, 4)],
  ['wall column 1x5', lattice(150, 140, 1, 5, 4)],
];

const summary = [];

for (const [name, material] of MATERIALS) {
  const bb = bboxOf(material);
  console.log(`=== ${name} ===  mass ${material.size}, bounding box ${bb.w}x${bb.h}`);

  // ---- what perturbations to try
  const removals1 = [...material];
  const additions1 = [];
  for (let y = bb.miny - 1; y <= bb.maxy + 1; y++) {
    for (let x = bb.minx - 1; x <= bb.maxx + 1; x++) {
      const k = wrapCell(x, y);
      if (!material.has(k)) additions1.push(k);
    }
  }
  // two-cell removals, only from within one block, sampled: every pair of cells
  // of the material that are neighbours
  const removals2 = [];
  for (const a of material) {
    const ax = a % SIZE, ay = (a - ax) / SIZE;
    for (const b of material) {
      if (b <= a) continue;
      const bx = b % SIZE, by = (b - bx) / SIZE;
      if (Math.abs(ax - bx) <= 1 && Math.abs(ay - by) <= 1) removals2.push([a, b]);
    }
  }
  // vacancies: remove one whole block at a time (only the blocks of the material)
  const vacancies = [];
  for (let i = 0; i < 8; i++) for (let j = 0; j < 8; j++) {
    const b0 = wrapCell(bb.minx + i * 4, bb.miny + j * 4);
    if (!material.has(b0)) continue;
    const blk = blockAt(bb.minx + i * 4, bb.miny + j * 4);
    let all = true;
    for (const k of blk) if (!material.has(k)) all = false;
    if (all) vacancies.push(blk);
  }

  const tally = {};
  const examples = {};
  const record = (cls, r, label) => {
    tally[cls] = tally[cls] || {};
    tally[cls][r.kind] = (tally[cls][r.kind] || 0) + 1;
    if (r.kind !== 'RESTORED' && !examples[cls]) {
      examples[cls] = `${label}: ${r.kind} (removed ${r.removed}, added ${r.added}, emitted ${r.outside})`;
    }
  };

  for (const k of removals1) {
    const live = new Set(material); live.delete(k);
    record('remove 1 cell', runCase(material, live, STEPS, bb), `cell (${k % SIZE},${Math.floor(k / SIZE)})`);
  }
  for (const [a, b] of removals2) {
    const live = new Set(material); live.delete(a); live.delete(b);
    record('remove 2 cells', runCase(material, live, STEPS, bb), 'pair');
  }
  for (const k of additions1) {
    const live = new Set(material); live.add(k);
    record('add 1 cell', runCase(material, live, STEPS, bb), `cell (${k % SIZE},${Math.floor(k / SIZE)})`);
  }
  for (const blk of vacancies) {
    const live = new Set(material); for (const k of blk) live.delete(k);
    record('remove 1 whole block', runCase(material, live, STEPS, bb), 'vacancy');
  }

  console.log('  perturbation class      cases   restored   restored+emit   damaged   gone');
  for (const cls of ['remove 1 cell', 'remove 2 cells', 'add 1 cell', 'remove 1 whole block']) {
    const t = tally[cls];
    if (!t) continue;
    const g = (k) => t[k] || 0;
    const total = Object.values(t).reduce((a, b) => a + b, 0);
    console.log(`    ${cls.padEnd(22)}${String(total).padStart(5)}${String(g('RESTORED')).padStart(11)}${String(g('RESTORED+EMIT')).padStart(16)}${String(g('DAMAGED')).padStart(10)}${String(g('GONE')).padStart(7)}`);
    if (examples[cls]) console.log(`      e.g. ${examples[cls]}`);
  }
  const eRemove = tally['remove 1 cell'] ? (tally['remove 1 cell']['RESTORED'] || 0) / removals1.length : NaN;
  const eAdd = tally['add 1 cell'] ? (tally['add 1 cell']['RESTORED'] || 0) / additions1.length : NaN;
  const eVac = vacancies.length ? (tally['remove 1 whole block']['RESTORED'] || 0) / vacancies.length : NaN;
  console.log(`  E_remove(1 cell) = ${eRemove.toFixed(3)}   E_add(1 cell) = ${eAdd.toFixed(3)}   E_vacancy = ${Number.isNaN(eVac) ? 'n/a' : eVac.toFixed(3)}`);
  console.log('');
  summary.push({ name, mass: material.size, eRemove, eAdd, eVac });
}

console.log('=== SUMMARY ===');
console.log('  matter              mass   E_remove   E_add   E_vacancy');
for (const s of summary) {
  console.log(`  ${s.name.padEnd(20)}${String(s.mass).padStart(4)}${s.eRemove.toFixed(3).padStart(11)}${s.eAdd.toFixed(3).padStart(8)}${(Number.isNaN(s.eVac) ? 'n/a' : s.eVac.toFixed(3)).padStart(12)}`);
}
console.log('');
console.log('  A real solid is elastic in BOTH directions: E_remove ~ 1 AND E_add ~ 1.');
console.log('  Interpretation of the numbers above decides whether "make the wall bigger"');
console.log('  can ever turn Life matter into a screen that a particle can touch and survive.');
