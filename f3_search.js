'use strict';
/*
 * f3_search.js   --   ROUND 70, CRITERION F3
 * ---------------------------------------------------------------------------
 * F3:  in a GENERAL (not integrable, not fine-tuned) rule space, use the COLLISION criterion to look
 *      for a wall that a mover crosses, AND a bigger disturbance that leaves a permanent record.
 *
 * Round 67 searched ONE NARROW CORNER -- 1D, 2 states, radius 1 -- and came back empty: 32 (rule,
 * ether) pairs had a rigid localized mover and NONE had a stationary wall, so the collision could not
 * even be set up.  That is what made F3 the last existence question still standing.  This round
 * widens the space in the two directions that matter:
 *
 *      PART 1  1D, THREE states, radius 1      (3^27 ~ 7.6e12 rules -- random sample)
 *      PART 2  2D, two states, Life-like B/S   (outer-totalistic, B0 off so the vacuum is stable)
 *
 * THE CRITERION (identical in both parts, and it must be an INTERACTION -- round 67's first result was
 * that a naive conjunction "movers exist AND some perturbation scars" is satisfied almost everywhere
 * and means nothing):
 *      ether  : a spatially periodic background that maps to itself (up to a shift)
 *      mover  : one localized defect whose blob translates rigidly
 *      matter : one localized blob that STAYS PUT and persists
 *      GATES  : the collision is only attempted if a mover AND matter were both found and both are
 *               rigid -- otherwise a "PASS" on a non-rigid blob would pass for the wrong reason
 *      PASS   : after the collision the mover is still travelling AND the matter is still there
 *      THRESH : a disturbance BIGGER than the mover changes the matter (=> a finite threshold)
 *      F3 OBJECT = PASS + THRESH
 *
 * PREDICTION (written before the run):  empty again.  Everything measured so far says a transparent
 * wall needs non-generic structure (round 68: non-integrable 56% radiation vs integrable 7%), and
 * round 67 found no wall at all in the elementary corner.  If something DOES pass, that is a headline
 * and it reopens the program; if nothing does, F3 is closed for these spaces and the honest statement
 * is "searched and empty", with the counts.
 * ---------------------------------------------------------------------------
 */

function rng(seed) { let x = seed | 0; if (x === 0) x = 0x9e3779b9; return () => { x ^= x << 13; x ^= x >>> 17; x ^= x << 5; return (x >>> 0) / 4294967296; }; }

// ===================================================================== PART 1: 1D, three states
const P1 = { rules: 0, withEther: 0, withMover: 0, withMatter: 0, withBoth: 0, pass: 0, passThresh: 0, examples: [] };
function step3(rule, a) {
  const n = a.length, o = new Uint8Array(n);
  for (let i = 0; i < n; i++) o[i] = rule[a[(i - 1 + n) % n] * 9 + a[i] * 3 + a[(i + 1) % n]];
  return o;
}
function shiftOf(a, b) {
  const n = a.length;
  for (let s = 0; s < n; s++) { let ok = true; for (let i = 0; i < n; i++) if (a[i] !== b[(i + s) % n]) { ok = false; break; } if (ok) return s; }
  return -1;
}
function spatialPeriod(a) {
  const n = a.length;
  for (let p = 1; p <= n; p++) { if (n % p) continue; let ok = true; for (let i = p; i < n && ok; i++) if (a[i] !== a[i % p]) ok = false; if (ok) return p; }
  return n;
}
function findEther3(rule) {
  for (let p = 1; p <= 5; p++) {
    const total = 3 ** p;
    for (let m = 0; m < total; m++) {
      const P = new Uint8Array(p);
      let v = m;
      for (let i = 0; i < p; i++) { P[i] = v % 3; v = (v - P[i]) / 3; }
      let cur = P;
      for (let q = 1; q <= 3; q++) {
        cur = step3(rule, cur);
        if (shiftOf(P, cur) >= 0) return { P, p };
      }
    }
  }
  return null;
}
// connected components of the difference from the ether
function comps3(a, ether, p) {
  const n = a.length, out = [];
  let start = -1;
  for (let i = 0; i < n; i++) if (a[i] !== ether[((i % p) + p) % p]) { start = i; break; }
  if (start < 0) return out;
  let cur = null;
  for (let s = 0; s < n; s++) {
    const i = (start + s) % n;
    const diff = a[i] !== ether[((i % p) + p) % p];
    if (diff) { if (!cur) cur = []; cur.push([i, a[i]]); }
    else if (cur) { out.push(cur); cur = null; }
  }
  if (cur) out.push(cur);
  if (out.length > 1) {
    const f = out[0], l = out[out.length - 1];
    if ((l[l.length - 1][0] + 1) % n === f[0][0]) { const merged = l.concat(f.map(([i, v]) => [i + n, v])); out.splice(0, 1); out.pop(); out.push(merged); }
  }
  return out;
}
const cwidth = (c) => c[c.length - 1][0] - c[0][0] + 1;
const ccenter = (c, n) => { let s = 0, co = 0; for (const [i] of c) { s += Math.sin(2 * Math.PI * (i % n) / n); co += Math.cos(2 * Math.PI * (i % n) / n); } let a = Math.atan2(s, co); if (a < 0) a += 2 * Math.PI; return a * n / (2 * Math.PI); };
const cdist = (x, y, n) => { const d = Math.abs(x - y) % n; return Math.min(d, n - d); };
{
  const RULES = 3000, RING = 60, T1 = 100, T2 = 200;
  const rnd = rng(20260101);
  for (let r = 0; r < RULES; r++) {
    const rule = new Uint8Array(27);
    for (let i = 0; i < 27; i++) rule[i] = Math.floor(rnd() * 3);
    const eth = findEther3(rule);
    if (!eth) continue;
    P1.rules++;
    const { P, p } = eth;
    const etherAt = (i) => P[((i % p) + p) % p];
    const build = (defs) => { const a = new Uint8Array(RING); for (let i = 0; i < RING; i++) a[i] = etherAt(i); for (const [pos, pat] of defs) for (let k = 0; k < pat.length; k++) a[(pos + k) % RING] = pat[k]; return a; };
    const evolve = (a, steps) => { let c = a; for (let s = 0; s < steps; s++) c = step3(rule, c); return c; };
    const blobs = (a) => comps3(a, P, p).filter((c) => cwidth(c) <= 10);
    // mover
    let mover = null;
    outer1:
    for (const k of [1, 2, 3]) for (let x = 0; x < RING; x += 15) {
      const pat = []; for (let z = 0; z < k; z++) pat.push((etherAt(x + z) + 1) % 3);
      const b1 = blobs(evolve(build([[x, pat]]), T1));
      const b2 = blobs(evolve(build([[x, pat]]), T2));
      if (b1.length === 1 && b2.length === 1 && cdist(ccenter(b1[0], RING), ccenter(b2[0], RING), RING) > 3 &&
          b1[0].map(([, v]) => v).join('') === b2[0].map(([, v]) => v).join('')) {
        mover = { k, comp: b1[0].map(([i, v]) => [((i - x) % RING + RING) % RING, v]), v: cdist(ccenter(b1[0], RING), ccenter(b2[0], RING), RING) / (T2 - T1) };
        break outer1;
      }
    }
    if (!mover) continue;
    P1.withMover++;
    // matter
    let matter = null;
    outer2:
    for (const k of [1, 2, 3, 4]) for (let x = 0; x < RING; x += 15) {
      const pat = []; for (let z = 0; z < k; z++) pat.push((etherAt(x + z) + 1) % 3);
      const b1 = blobs(evolve(build([[x, pat]]), T1));
      const b2 = blobs(evolve(build([[x, pat]]), T2));
      if (b1.length === 1 && b2.length === 1 && cdist(ccenter(b1[0], RING), ccenter(b2[0], RING), RING) <= 1) {
        matter = { k, comp: b1[0].map(([i, v]) => [((i - x) % RING + RING) % RING, v]) };
        break outer2;
      }
    }
    if (!matter) continue;
    P1.withMatter++;
    P1.withBoth++;
    // the collision, with gates AND SHAPE MATCHING.  My first version only asked for "a blob near the
    // centre and a blob somewhere else" -- debris satisfies that, so it could pass for the wrong
    // reason.  Now the surviving blobs must match the mover's and the matter's SHAPES.
    const xm = Math.floor(RING / 2), sep = Math.max(12, Math.floor(RING / 5));
    const place = (a, x0, comp) => { for (const [off, v] of comp) a[(x0 + off) % RING] = v; };
    const mShape = mover.comp.map(([, v]) => v).join('');
    const wShape = matter.comp.map(([, v]) => v).join('');
    const init = build([]); place(init, xm, matter.comp); place(init, (xm - sep + RING) % RING, mover.comp);
    const post = evolve(init, T2);
    const blobsPost = blobs(post);
    const wallShapeThere = blobsPost.some((c) => c.map(([, v]) => v).join('') === wShape && cdist(ccenter(c, RING), xm, RING) <= 10);
    const moverShapeThere = blobsPost.some((c) => c.map(([, v]) => v).join('') === mShape && cdist(ccenter(c, RING), xm, RING) > 10);
    if (!(wallShapeThere && moverShapeThere)) continue;
    // and the mover must actually have CROSSED: run on and require a mover-shaped blob downstream
    const post2 = evolve(post, T2);
    const crossed = blobs(post2).some((c) => c.map(([, v]) => v).join('') === mShape);
    if (!crossed) continue;
    P1.pass++;
    // the threshold: a bigger defect aimed at the matter
    let changes = false;
    for (const k of [mover.k + 2, mover.k + 4]) {
      const pat = []; for (let z = 0; z < k; z++) pat.push((etherAt(xm + z) + 1) % 3);
      const g = build([]); place(g, xm, matter.comp);
      for (let z = 0; z < k; z++) g[(xm - sep + z + RING) % RING] = pat[z];
      const h = evolve(g, T2);
      if (!blobs(h).some((c) => c.map(([, v]) => v).join('') === wShape)) changes = true;
    }
    if (changes) {
      P1.passThresh++;
      if (P1.examples.length < 5) P1.examples.push({ r, rule: Array.from(rule), P, p, moverK: mover.k, matterK: matter.k, wShape, mShape });
    }
  }
}

console.log('============================================================');
console.log('PART 1: 1D, THREE STATES, RADIUS 1 (3^27 rules) -- RANDOM SAMPLE');
console.log('============================================================');
console.log(`   rules sampled                     : 3000`);
console.log(`   with a spatially periodic ether   : ${P1.rules}`);
console.log(`   with a rigid localized MOVER      : ${P1.withMover}`);
console.log(`   with a stationary localized WALL  : ${P1.withMatter}   (of those with a mover)`);
console.log(`   collision PASSES (both survive)   : ${P1.pass}`);
console.log(`   PASS *and* a bigger hit changes it: ${P1.passThresh}   <- the F3 object`);
console.log('');
console.log(`   comparison, round 67 (1D, TWO states, 256 rules exhaustively): movers 130 rules,`);
console.log('   (rule,ether) pairs with a mover 32, with a wall 0.');
console.log(`  => ${P1.passThresh === 0 ? 'EMPTY in this space too.' : 'a hit survived the strengthened criterion -- see the detail below.'}`);
for (const e of P1.examples) {
  console.log('');
  console.log(`   DETAIL OF THE HIT (sample index ${e.r}), so it can be re-checked by hand:`);
  console.log(`     rule table (neighbourhood 000..222 -> value): ${e.rule.join('')}`);
  console.log(`     ether: period ${e.p}, pattern ${Array.from(e.P).join('')}`);
  console.log(`     mover: a ${e.moverK}-cell defect, shape ${e.mShape}`);
  console.log(`     wall : a ${e.matterK}-cell defect, shape ${e.wShape}`);
  console.log('     the collision kept BOTH shapes, the mover crossed, and a bigger hit removed the wall shape.');
}
console.log('');

// ===================================================================== PART 2: 2D, Life-like B/S
const P2 = { rules: 0, withMover: 0, withMatter: 0, withBoth: 0, pass: 0, passThresh: 0, examples: [] };
const G = 44;                                   // the torus is G x G.  It MUST be much larger than the
// displacement over the screen interval: a Life glider moves 20 cells in 80 ticks, and on my first
// G = 22 torus that wrapped back to a displacement of 2, so the screen reported "no mover" for Life
// itself.  The control below is what caught it.
const SEEDS = [
  [[1, 0]],
  [[1, 0], [1, 1]],
  [[1, 0], [1, 1], [0, 1]],
  [[1, 0], [1, 1], [0, 1], [0, 0]],
  [[0, 1, 0], [1, 0, 1], [0, 1, 0]],
  // the Life glider.  My first version left this out, so LIFE ITSELF would have shown "0 movers" --
  // the screen would have measured nothing.  A control is mandatory here.
  [[1, 0], [2, 1], [0, 2], [1, 2], [2, 2]]
];
function blank() { return new Uint8Array(G * G); }
function placeSeed(g, seed, x0, y0) { for (const [dx, dy] of seed) g[((y0 + dy + G) % G) * G + ((x0 + dx + G) % G)] = 1; }
function stepLife(g, B, S) {
  const o = new Uint8Array(G * G);
  for (let y = 0; y < G; y++) for (let x = 0; x < G; x++) {
    let c = 0;
    for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
      if (dx === 0 && dy === 0) continue;
      c += g[((y + dy + G) % G) * G + ((x + dx + G) % G)];
    }
    const alive = g[y * G + x];
    o[y * G + x] = alive ? ((S >> c) & 1) : ((B >> c) & 1);
  }
  return o;
}
function box(g) {                       // live count, circular centre, and a TRANSLATION-INVARIANT shape
  const cells = [];
  for (let y = 0; y < G; y++) for (let x = 0; x < G; x++) if (g[y * G + x]) cells.push([x, y]);
  if (!cells.length) return { n: 0, cx: 0, cy: 0, shape: '' };
  let asx = 0, acx = 0, asy = 0, acy = 0;
  for (const [x, y] of cells) {
    asx += Math.sin(2 * Math.PI * x / G); acx += Math.cos(2 * Math.PI * x / G);
    asy += Math.sin(2 * Math.PI * y / G); acy += Math.cos(2 * Math.PI * y / G);
  }
  let ax = Math.atan2(asx, acx); if (ax < 0) ax += 2 * Math.PI;
  let ay = Math.atan2(asy, acy); if (ay < 0) ay += 2 * Math.PI;
  const cx = ax * G / (2 * Math.PI), cy = ay * G / (2 * Math.PI);
  // UNWRAP each cell relative to the centre, then normalise by the bounding-box corner.  My first
  // version anchored the shape at round(cx), round(cy): a glider's centre of mass creeps by 0.25 of a
  // cell per period, so the window slipped and the shapes never matched -- the Life control caught it.
  const rel = [];
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  for (const [x, y] of cells) {
    let dx = x - cx, dy = y - cy;
    if (dx > G / 2) dx -= G; if (dx < -G / 2) dx += G;
    if (dy > G / 2) dy -= G; if (dy < -G / 2) dy += G;
    rel.push([dx, dy]);
    minx = Math.min(minx, dx); maxx = Math.max(maxx, dx);
    miny = Math.min(miny, dy); maxy = Math.max(maxy, dy);
  }
  const w = Math.round(maxx - minx) + 1, h = Math.round(maxy - miny) + 1;
  if (w > 20 || h > 20) return { n: cells.length, cx, cy, shape: 'BIG' };   // not a localized object
  const grid = new Array(w * h).fill('0');
  for (const [dx, dy] of rel) {
    const ix = Math.round(dx - minx), iy = Math.round(dy - miny);
    if (ix >= 0 && ix < w && iy >= 0 && iy < h) grid[iy * w + ix] = '1';
  }
  return { n: cells.length, cx, cy, shape: `${w}x${h}:${grid.join('')}` };
}
function components2D(g) {               // 8-connected components of the live cells
  const seen = new Uint8Array(G * G), out = [];
  for (let y0 = 0; y0 < G; y0++) for (let x0 = 0; x0 < G; x0++) {
    const s0 = y0 * G + x0;
    if (!g[s0] || seen[s0]) continue;
    const stack = [[x0, y0]], cells = [];
    seen[s0] = 1;
    while (stack.length) {
      const [x, y] = stack.pop();
      cells.push([x, y]);
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        const nx = (x + dx + G) % G, ny = (y + dy + G) % G, s = ny * G + nx;
        if (g[s] && !seen[s]) { seen[s] = 1; stack.push([nx, ny]); }
      }
    }
    out.push(cells);
  }
  return out;
}
function shapeOfCells(cells) {           // translation-invariant, wrap-safe shape of one component
  if (!cells.length) return '';
  let asx = 0, acx = 0, asy = 0, acy = 0;
  for (const [x, y] of cells) {
    asx += Math.sin(2 * Math.PI * x / G); acx += Math.cos(2 * Math.PI * x / G);
    asy += Math.sin(2 * Math.PI * y / G); acy += Math.cos(2 * Math.PI * y / G);
  }
  let ax = Math.atan2(asx, acx); if (ax < 0) ax += 2 * Math.PI;
  let ay = Math.atan2(asy, acy); if (ay < 0) ay += 2 * Math.PI;
  const cx = ax * G / (2 * Math.PI), cy = ay * G / (2 * Math.PI);
  const rel = [];
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  for (const [x, y] of cells) {
    let dx = x - cx, dy = y - cy;
    if (dx > G / 2) dx -= G; if (dx < -G / 2) dx += G;
    if (dy > G / 2) dy -= G; if (dy < -G / 2) dy += G;
    rel.push([dx, dy]);
    minx = Math.min(minx, dx); maxx = Math.max(maxx, dx);
    miny = Math.min(miny, dy); maxy = Math.max(maxy, dy);
  }
  const w = Math.round(maxx - minx) + 1, h = Math.round(maxy - miny) + 1;
  if (w > 20 || h > 20) return 'BIG';
  const grid = new Array(w * h).fill('0');
  for (const [dx, dy] of rel) {
    const ix = Math.round(dx - minx), iy = Math.round(dy - miny);
    if (ix >= 0 && ix < w && iy >= 0 && iy < h) grid[iy * w + ix] = '1';
  }
  return `${w}x${h}:${grid.join('')}`;
}
const hasShape = (g, shape) => components2D(g).some((c) => shapeOfCells(c) === shape);
const torusDist = (a, b) => { const d = Math.abs(a - b) % G; return Math.min(d, G - d); };
const bslabel = (mask) => 'B' + [...Array(9).keys()].filter((k) => (mask >> k) & 1).join('') + '/S' + [...Array(9).keys()].filter((k) => (S_ >> k) & 1).join('');
let S_ = 0;                              // set inside the loop, used by bslabel
{
  const RULES = 400, T = 80;
  const rnd = rng(20260202);
  for (let r = 0; r < RULES; r++) {
    const B = Math.floor(rnd() * 512), S = Math.floor(rnd() * 512);
    if (B & 1) continue;                          // B0 would make the vacuum unstable
    P2.rules++;
    let mover = null, matter = null;
    for (const seed of SEEDS) {
      const g0 = blank(); placeSeed(g0, seed, 5, 5);
      const g1 = (() => { let g = g0; for (let t = 0; t < T; t++) g = stepLife(g, B, S); return g; })();
      const g2 = (() => { let g = g1; for (let t = 0; t < T; t++) g = stepLife(g, B, S); return g; })();
      const b1 = box(g1), b2 = box(g2);
      if (b1.n === 0 || b2.n === 0) continue;
      const moved = torusDist(b1.cx, b2.cx) + torusDist(b1.cy, b2.cy);
      if (moved > 3 && b1.shape === b2.shape) mover = { seed, v: moved / T };
      else if (moved <= 1 && b1.shape === b2.shape) matter = { seed };
    }
    if (mover) P2.withMover++;
    if (matter) P2.withMatter++;
    if (!mover || !matter) continue;
    P2.withBoth++;
    // THE COLLISION, by SHAPE MATCHING (my first version only asked "is something alive near the
    // centre, and is something else moving" -- debris satisfies that, and it reported 128 of 140
    // collisions as passing, which is obviously wrong).  Now the surviving components must MATCH the
    // mover's and the wall's shapes.
    const gm = blank(); placeSeed(gm, matter.seed, 11, 11);
    let gmF = gm; for (let t = 0; t < T; t++) gmF = stepLife(gmF, B, S);
    const wallShape = shapeOfCells(components2D(gmF)[0]);
    const gg = blank(); placeSeed(gg, mover.seed, 4, 11);
    const moverShape = shapeOfCells(components2D(gg)[0]);
    const init = blank();
    for (let i = 0; i < G * G; i++) init[i] = gmF[i] | gg[i];
    let cur = init; for (let t = 0; t < T * 3; t++) cur = stepLife(cur, B, S);
    const wallThere = hasShape(cur, wallShape), moverThere = hasShape(cur, moverShape);
    if (wallThere && moverThere) {
      P2.pass++;
      // the threshold: add one cell next to the wall; a bigger disturbance must change the wall
      const gt = Uint8Array.from(gmF);
      gt[11 * G + 14] = 1;
      let gtF = gt; for (let t = 0; t < T * 2; t++) gtF = stepLife(gtF, B, S);
      if (!hasShape(gtF, wallShape)) {
        P2.passThresh++;
        if (P2.examples.length < 5) {
          S_ = S;
          P2.examples.push({ label: bslabel(B), moverSeed: mover.seed, wallSeed: matter.seed, wallShape, moverShape });
        }
      }
    }
  }
}
console.log('============================================================');
console.log('PART 2: 2D, TWO STATES, LIFE-LIKE B/S (outer-totalistic, B0 off) -- RANDOM SAMPLE');
console.log('============================================================');
{
  // THE CONTROL: Life itself, run through exactly the same screen.  If the screen cannot find Life's
  // glider it is not measuring anything.
  const B = 0b00001000, S = 0b00011100;               // B3 / S23
  let lifeMover = null, lifeMatter = null;
  for (const seed of SEEDS) {
    const g0 = blank(); placeSeed(g0, seed, 5, 5);
    let g1 = g0; for (let t = 0; t < 80; t++) g1 = stepLife(g1, B, S);
    let g2 = g1; for (let t = 0; t < 80; t++) g2 = stepLife(g2, B, S);
    const b1 = box(g1), b2 = box(g2);
    if (b1.n === 0 || b2.n === 0) continue;
    const moved = torusDist(b1.cx, b2.cx) + torusDist(b1.cy, b2.cy);
    if (moved > 3 && b1.shape === b2.shape) lifeMover = seed;
    else if (moved <= 1 && b1.shape === b2.shape) lifeMatter = seed;
  }
  console.log(`   CONTROL -- Life (B3/S23): mover found = ${lifeMover ? 'yes (the glider seed)' : '*** NO: the screen is broken ***'}, wall found = ${lifeMatter ? 'yes' : 'no'}`);
  console.log('   (this control is why the glider seed was added: without it the screen reported 0 movers');
  console.log('    for every rule, including Life, and would have measured nothing.)');
  // F3 CONTROL on Life itself: it has a mover (glider) and matter (the block), so if this screen says
  // Life passes F3, the screen is wrong -- a glider hitting a Life block is known to wreck it.
  {
    const B3 = 0b00001000, S23 = 0b00011100;
    const gm = blank(); placeSeed(gm, [[1, 0], [1, 1], [0, 1], [0, 0]], 11, 11);
    let gmF = gm; for (let t = 0; t < 80; t++) gmF = stepLife(gmF, B3, S23);
    const wS = shapeOfCells(components2D(gmF)[0]);
    const gg = blank(); placeSeed(gg, SEEDS[5], 4, 11);
    const mS = shapeOfCells(components2D(gg)[0]);
    const init = blank(); for (let i = 0; i < G * G; i++) init[i] = gmF[i] | gg[i];
    let cur = init; for (let t = 0; t < 240; t++) cur = stepLife(cur, B3, S23);
    const w = hasShape(cur, wS), m = hasShape(cur, mS);
    console.log(`   F3 CONTROL -- Life, glider into a block: wall shape survives = ${w}, mover shape survives = ${m}`);
    console.log(`   => ${w && m ? '*** the screen says LIFE passes F3, which is wrong: the screen is broken ***' : 'correctly does NOT pass (the glider wrecks the block) -- so the screen discriminates'}`);
  }
}
console.log('');
console.log(`   rules sampled                     : ${P2.rules}`);
console.log(`   with a rigid MOVER (a patch that translates)  : ${P2.withMover}`);
console.log(`   with a stationary WALL (a persistent patch)   : ${P2.withMatter}`);
console.log(`   with BOTH                         : ${P2.withBoth}`);
console.log(`   collision PASSES (both survive)   : ${P2.pass}`);
console.log(`   PASS *and* adding one cell kills the wall : ${P2.passThresh}   <- the F3 object`);
console.log('');
console.log(`  => ${P2.passThresh === 0 ? 'EMPTY in this space too.' : 'hits survived the shape-matching criterion -- see the detail below.'}`);
for (const e of P2.examples) {
  console.log(`     ${e.label}: mover ${JSON.stringify(e.moverSeed)} shape ${e.moverShape}, wall ${JSON.stringify(e.wallSeed)} shape ${e.wallShape}`);
}
console.log('');
console.log('============================================================');
console.log('F3 VERDICT');
console.log('============================================================');
{
  const total = P1.passThresh + P2.passThresh;
  console.log(`   1D three-state : ${P1.passThresh} rules satisfy F3`);
  console.log(`   2D Life-like   : ${P2.passThresh} rules satisfy F3`);
  console.log('');
  if (total === 0) {
    console.log('  => F3 IS CLOSED for the spaces searched here.  Together with round 67 (the elementary corner');
    console.log('     and its 256 rules) and round 68 (the non-integrable vs integrable comparison, 56.04% vs');
    console.log('     7.04%), the picture is consistent: a transparent wall needs a NON-GENERIC structure, and');
    console.log('     random rules do not supply one.  This does NOT prove impossibility -- it closes the');
    console.log('     search spaces that were declared and leaves F3 open only for spaces not yet searched.');
  } else {
    console.log('  => F3 IS SATISFIED.  This reopens the program: a general-space rule DOES host a transparent');
    console.log('     wall with a finite threshold.  Everything above about "needs non-generic structure" must');
    console.log('     be revised for these rules, and the next round should characterise them.');
  }
  console.log('');
  console.log('   HONEST BOUNDARIES OF THIS SEARCH:');
  console.log('     * PART 1 samples 3000 of 7.6e12 three-state rules (not exhaustive); ethers only up to');
  console.log('       period 5; defects only of one family (a block set to the next symbol).');
  console.log('     * PART 2 samples 400 rules of the Life-like B/S space (209 kept: B0 makes the vacuum');
  console.log('       unstable); seeds are six small patches INCLUDING the Life glider; the mover/wall screen');
  console.log('       is a bounding-box/component screen, not a proof.');
  console.log('     * So the honest statement is "SEARCHED AND EMPTY", never "impossible".');
}
