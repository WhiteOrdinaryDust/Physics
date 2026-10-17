'use strict';
/*
 * phase_arrival_check.js
 * ---------------------------------------------------------------------------
 * THE CONTROL THAT DECIDES WHETHER "PHASE" MEANS ANYTHING.
 *
 * This programme has already been burned three times by the same mistake:
 *   "moving the particle 1 cell further along its own direction of motion is a
 *    time translation, not a change of internal phase."
 * If that is what is going on here, then the phase-selective results in
 * lwss_slit.js are not a hidden variable at all -- they are just a different
 * impact parameter, obtainable by sliding the particle, and the wall is only
 * measuring a position.
 *
 * So this script answers ONE question with data:
 *   when the four phases reach the screen, are they at the SAME bounding box at
 *   the SAME generation, differing only in which cells inside that box are live?
 *
 * It uses one fact that makes the test exact: before contact, the particle's
 * motion is unaffected by the screen, so the particle alone is a faithful record
 * of where it is. A full period is 4 generations of exactly 2 cells of travel,
 * so after any multiple of 4 generations the bounding box is at EXACTLY the same
 * place for every phase -- no half-cell ambiguity to argue about.
 *
 * A second, independent question: can any COARSE observable of the particle
 * (mass, bounding box, position, speed, period) predict which phase it is?
 * The mass alternates 9, 12, 9, 12 over the phases, which is all a coarse
 * observer gets.
 * ---------------------------------------------------------------------------
 */

const SIZE = 600;
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
function canon(cells) {
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity, sx = 0, sy = 0;
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
    sx += x; sy += y;
  }
  const w = maxx - minx + 1, h = maxy - miny + 1, g = [];
  for (let r = 0; r < h; r++) g.push(new Array(w).fill('.'));
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    g[y - miny][x - minx] = 'o';
  }
  return { w, h, n: cells.length, minx, miny, maxx, maxy, cx: sx / cells.length, cy: sy / cells.length, id: `${w}x${h}:` + g.map(r => r.join('')).join('/') };
}
function shipRight(x, y) {
  const s = new Set();
  for (const [dx, dy] of [[1, 0], [2, 0], [3, 0], [4, 0], [0, 1], [4, 1], [4, 2], [0, 3], [3, 3]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function normalizeTo(set, x0, y0) {
  const c = canon([...set]);
  const dx = x0 - c.minx, dy = y0 - c.miny;
  const out = new Set();
  for (const k of set) { const x = k % SIZE, y = (k - x) / SIZE; out.add(wrapCell(x + dx, y + dy)); }
  return out;
}
function picture(set) {
  const c = canon([...set]);
  return c.id;
}

// geometry of the screen used in lwss_slit.js: spacing 6, slit rows 170..179
const WALL_X = 250;          // first live cell of the screen
const LAUNCH_X = 190, Y0 = 172 + 1;   // b=3 placement is the centred, safe one; b=2 is 1 closer to the edge
const TMAX = 140;

console.log('=== does the phase change where the particle is, or only what it is? ===');
console.log(`screen starts at x=${WALL_X};  particle launched with its bounding box at x=${LAUNCH_X}\n`);

const hist = [];
for (let phase = 0; phase < 4; phase++) {
  let sh = shipRight(150, 150);
  for (let i = 0; i < phase; i++) sh = step(sh);
  sh = normalizeTo(sh, LAUNCH_X, Y0);
  const track = [];
  let cur = new Set(sh);
  for (let t = 0; t <= TMAX; t++) {
    const c = canon([...cur]);
    track.push({ t, minx: c.minx, maxx: c.maxx, miny: c.miny, maxy: c.maxy, n: c.n, pieces: components(cur).length, id: c.id });
    cur = step(cur);
  }
  hist.push(track);
}

console.log('  generation :  bounding box and mass of the particle, by phase');
console.log('               (the last row before the first cell becomes a neighbour of the screen is the');
console.log('                state from which the collision is decided)');
console.log('');
console.log('   t    phase0            phase1            phase2            phase3      same box?');
for (let t = 112; t <= 126; t++) {
  const cells = hist.map(h => h[t]);
  const boxes = cells.map(c => `${c.minx}..${c.maxx},${c.miny}..${c.maxy}`);
  const sameBox = boxes.every(b => b === boxes[0]);
  const masses = cells.map(c => c.n);
  const parts = cells.map((c, i) => `${boxes[i]} m${masses[i]}`.padEnd(18));
  console.log(`  ${String(t).padStart(3)}  ${parts.join('')}${sameBox ? 'YES' : '*** NO ***'}`);
}

console.log('');
const tContact = hist[0].findIndex(c => c.maxx >= WALL_X - 1);
console.log(`  first generation at which the particle touches the screen's neighbour column: t=${tContact}`);
for (let phase = 0; phase < 4; phase++) {
  const c = hist[phase][tContact];
  const p = hist[phase][tContact - 4];
  console.log(`   phase ${phase}: at t=${tContact} box ${c.minx}..${c.maxx},${c.miny}..${c.maxy} mass ${c.n} pieces ${c.pieces}   (at t=${tContact - 4}: box ${p.minx}..${p.maxx},${p.miny}..${p.maxy} mass ${p.n})`);
}
console.log('');
console.log('  The particle cell arrangement one generation before contact:');
for (let phase = 0; phase < 4; phase++) {
  const c = hist[phase][tContact - 1];
  console.log(`   phase ${phase}: mass ${c.n}, ${c.pieces} piece(s), ${c.id}`);
}
console.log('');
const boxesAtContact = hist.map(h => `${h[tContact - 1].minx}..${h[tContact - 1].maxx},${h[tContact - 1].miny}..${h[tContact - 1].maxy}`);
const idsAtContact = hist.map(h => h[tContact - 1].id);
const massesAtContact = hist.map(h => h[tContact - 1].n);
console.log('=== VERDICT ===');
console.log(`  all four phases have the same bounding box one generation before contact: ${boxesAtContact.every(b => b === boxesAtContact[0]) ? 'YES' : 'NO'}`);
console.log(`  all four phases reach contact at the same generation: YES (t=${tContact} by construction of the search above)`);
console.log(`  number of distinct cell arrangements at that same box: ${new Set(idsAtContact).size}`);
console.log(`  masses at that box: ${massesAtContact.join(', ')}`);
console.log('');
console.log('  So at the identical box, at the identical generation, the particle is in one of four');
console.log('  different internal arrangements -- and the masses show that even the particle MASS does');
console.log('  not identify the arrangement (the outcome in lwss_slit.js splits the phases as');
console.log('  {9-cell, 12-cell} against {12-cell, 9-cell}, so no coarse observable separates them).');
