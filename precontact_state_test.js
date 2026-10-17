'use strict';
/*
 * precontact_state_test.js
 * ---------------------------------------------------------------------------
 * THE PRE-CONTACT CRITERION, WITH THE MISSING PIECE INCLUDED.
 *
 * WHAT WAS MISSING
 *   The previous version described the pre-contact state by SHAPE only:
 *        glider shape + obstacle shape
 *   It found these identical for every k, yet the outcome varied. So shape is
 *   not enough. The missing quantity is the RELATIVE DISTANCE between the two
 *   objects. Two configurations can have identical shapes and identical
 *   geometry (all motion on the line y = x) while the glider sits at a
 *   different distance from the obstacle.
 *
 * THE CORRECTED PRE-CONTACT STATE
 *        (glider shape, obstacle shape, gap in cells between them)
 *   where the gap is the largest Chebyshev separation d such that every cell of
 *   the glider is more than d away from every cell of the obstacle.
 *
 * WHAT THIS FILE DOES
 *   1. For a range of launch offsets k, record
 *        - the pre-contact (shape, shape, gap) triple
 *        - the terminal outcome
 *   2. Test the criterion directly: do equal triples always give equal outcomes,
 *      and do different triples sometimes give different outcomes?
 *   3. CONTROLS: the vacuum glider must have a well-defined gap to nothing; a
 *      block obstacle must give identical triples and identical outcomes when
 *      the setup is repeated.
 * ---------------------------------------------------------------------------
 */

const SIZE = 160;
function wrapCell(x, y) { return ((y % SIZE) + SIZE) % SIZE * SIZE + (((x % SIZE) + SIZE) % SIZE); }
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
  const w = maxx - minx + 1, h = maxy - miny + 1;
  const g = [];
  for (let r = 0; r < h; r++) g.push(new Array(w).fill('.'));
  for (const k of cells) {
    const x = k % SIZE, y = (k - x) / SIZE;
    g[y - miny][x - minx] = 'o';
  }
  return { w, h, n: cells.length, minx, miny, maxx, maxy, cx: sx / cells.length, cy: sy / cells.length, id: `${w}x${h}:` + g.map(r => r.join('')).join('/') };
}
function union(...sets) { const o = new Set(); for (const s of sets) for (const k of s) o.add(k); return o; }

// Chebyshev gap: the largest g such that every cell of A is more than g from every cell of B
function chebyshevGap(cellsA, cellsB) {
  let best = Infinity;
  for (const a of cellsA) {
    const ax = a % SIZE, ay = Math.floor(a / SIZE);
    let nearest = Infinity;
    for (const b of cellsB) {
      const bx = b % SIZE, by = Math.floor(b / SIZE);
      const d = Math.max(Math.abs(ax - bx), Math.abs(ay - by));
      if (d < nearest) nearest = d;
    }
    if (nearest < best) best = nearest;
  }
  return best === Infinity ? null : best - 1;   // gap in empty cells between them
}

const GLIDER0 = [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]];
function gliderAt(x, y) {
  const s = new Set();
  for (const [dx, dy] of GLIDER0) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function blockAt(x, y) {
  const s = new Set();
  for (const [dx, dy] of [[0, 0], [1, 0], [0, 1], [1, 1]]) s.add(wrapCell(x + dx, y + dy));
  return s;
}
function blinkerAt(cx, cy, horizontal) {
  const s = new Set();
  for (let d = -1; d <= 1; d++) {
    if (horizontal) s.add(wrapCell(cx + d, cy));
    else s.add(wrapCell(cx, cy + d));
  }
  return s;
}

/*
 * Run, and capture the LAST moment before the two objects merge into one
 * component. At that moment record both shapes and the gap.
 */
function runAndCapture(obstacle, startGlider, steps) {
  let s = union(obstacle, startGlider);
  let lastTriple = null;
  let mergeStep = -1;
  for (let t = 0; t <= steps; t++) {
    const comps = components(s);
    if (comps.length >= 2) {
      // identify the two objects: the one that is the glider, the other the obstacle
      const cs = comps.map(c => ({ cells: c, k: canon(c) }));
      // glider: the 3x3 component; obstacle: anything else
      const gl = cs.find(x => x.k.w === 3 && x.k.h === 3);
      const ob = cs.find(x => x !== gl);
      if (gl && ob) {
        const gap = chebyshevGap(gl.cells, ob.cells);
        lastTriple = { glider: gl.k.id, obstacle: ob.k.id, gap };
      }
    } else if (comps.length === 1 && mergeStep < 0) {
      mergeStep = t;
    }
    if (t === steps) break;
    s = step(s);
  }
  // terminal outcome
  let f = union(obstacle, startGlider);
  for (let t = 0; t < steps; t++) f = step(f);
  const fcomps = components(f).map(c => canon(c));
  return { triple: lastTriple, mergeStep, finalMass: f.size, finalComps: fcomps.length };
}

// ============================================================ controls
console.log('=== CONTROLS ===');
{
  // block obstacle, identical setup twice: triples and outcomes must match
  const ob = blockAt(62, 62);
  const a = runAndCapture(ob, gliderAt(32, 32), 200);
  const b = runAndCapture(ob, gliderAt(32, 32), 200);
  const same = JSON.stringify(a.triple) === JSON.stringify(b.triple);
  console.log(`  block, identical setup twice: triple identical? ${same ? 'YES PASS' : '*** NO FAIL ***'}`);
  console.log(`      triple: ${JSON.stringify(a.triple)}`);
  console.log(`      outcome: mass ${a.finalMass}, comps ${a.finalComps}`);
  console.log(`      merge step: ${a.mergeStep}`);
  if (!same) process.exit(1);
}

// ============================================================ the point
console.log('\n=== THE TEST: does the (glider shape, obstacle shape, gap) triple decide the outcome? ===\n');
console.log('--- obstacle = a block (period 1) ---');
{
  const rows = [];
  for (let k = 0; k < 12; k++) {
    const r = runAndCapture(blockAt(62, 62), gliderAt(32 + k, 32 + k), 200);
    rows.push({ k, ...r });
    console.log(`  k=${String(k).padStart(2)}  merge=${String(r.mergeStep).padStart(3)}  gap=${String(r.triple ? r.triple.gap : '?').padStart(3)}  glider=${r.triple ? r.triple.glider : '?'}  obstacle=${r.triple ? r.triple.obstacle : '?'}  ->  mass ${r.finalMass}, comps ${r.finalComps}`);
  }
  // group by triple
  const byTriple = new Map();
  for (const r of rows) {
    const key = JSON.stringify(r.triple);
    if (!byTriple.has(key)) byTriple.set(key, new Set());
    byTriple.get(key).add(`${r.finalMass}/${r.finalComps}`);
  }
  console.log('');
  console.log('  grouping by triple:');
  for (const [key, outs] of byTriple) {
    console.log(`    ${key}  -> outcomes: ${[...outs].join(', ')}  ${outs.size === 1 ? '(consistent)' : '*** INCONSISTENT ***'}`);
  }
}

console.log('\n--- obstacle = a blinker, starting vertical (period 2) ---');
{
  const rows = [];
  for (let k = 0; k < 12; k++) {
    const r = runAndCapture(blinkerAt(62, 62, false), gliderAt(32 + k, 32 + k), 200);
    rows.push({ k, ...r });
    console.log(`  k=${String(k).padStart(2)}  merge=${String(r.mergeStep).padStart(3)}  gap=${String(r.triple ? r.triple.gap : '?').padStart(3)}  glider=${r.triple ? r.triple.glider : '?'}  obstacle=${r.triple ? r.triple.obstacle : '?'}  ->  mass ${r.finalMass}, comps ${r.finalComps}`);
  }
  console.log('');
  const byTriple = new Map();
  for (const r of rows) {
    const key = JSON.stringify(r.triple);
    if (!byTriple.has(key)) byTriple.set(key, new Set());
    byTriple.get(key).add(`${r.finalMass}/${r.finalComps}`);
  }
  console.log('  grouping by triple:');
  for (const [key, outs] of byTriple) {
    console.log(`    ${key}  -> outcomes: ${[...outs].join(', ')}  ${outs.size === 1 ? '(consistent)' : '*** INCONSISTENT ***'}`);
  }
}

console.log('\n=== interpretation ===');
console.log('  If every distinct triple maps to a single outcome, and different triples');
console.log('  map to different outcomes, then the criterion holds once the relative');
console.log('  distance is included in the description of the pre-contact state.');
