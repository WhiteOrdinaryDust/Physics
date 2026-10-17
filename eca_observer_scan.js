'use strict';
/*
 * eca_observer_scan.js
 * ---------------------------------------------------------------------------
 * STEP 1 of the "internal observer" program.
 *
 * Ontology assumed (nothing more):
 *   - a deterministic, local, finite-speed cellular automaton
 *   - a global configuration exists at every time step (the "outside" view)
 *   - an "observer" is NOT assumed. It is a candidate window to be screened.
 *
 * A candidate observer is a space-time window:
 *   sensor   = L adjacent cells            (its finite interface)
 *   memory   = the block's last k states   (its finite memory)
 *   target   = the same block's next state (what it is coupled to)
 *
 * For every (rule, L, k) we estimate, with a Miller-Madow bias correction:
 *   complexity   C = H(W_next)                  bits of "what the world does to me"
 *   self-pred    I(W_k ; W_next)                raw predictive power about itself
 *   self-density = I(W_k ; W_next) / L          per-cell predictive budget
 *   info-closure = [H(W_next|W_k) - H(W_next|W_1)] / H(W_next)
 *                  how much the k-step memory beats a memoryless observer
 *   adv-density  = [I(W_k;W_next) - I(W_k;W_k(+D))]/L
 *                  self-model advantage over modelling a causally-shielded region
 *   observerScore= self-density * (clamped complexity) * (clamped closure)
 *
 * CONTROLS:
 *   R30  random      -> chaotic / statistically uniform (negative control)
 *   R110 glider     -> coherent moving structures  (positive control)
 *   R110 random
 *   R90  random
 *
 * The question the instrument answers:
 *   Is there ANY peak structure in (L,k) space, or is the landscape flat?
 * ---------------------------------------------------------------------------
 */

const RULES = [30, 90, 110, 184];
const MODES = { R30: 'random', R90: 'random', R110: 'glider', R184: 'random' };
const SEEDS = [12345, 999331];
const W = 256;          // ring width
const T = 2000;         // time steps
const D = 8;            // displacement for the causally-shielded control region
const LS = [2, 3, 5];   // L=8 costs 2^32 state space per pass; left out of the default run
const KS = [1, 2, 3];
const NEIGH = { 30: 1, 90: 1, 110: 1, 184: 1 };

// ---------------------------------------------------------------- PRNG (mulberry32)
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// ---------------------------------------------------------------- rule table
function makeRuleTable(rule) {
  const tbl = new Uint8Array(8);
  for (let n = 0; n < 8; n++) tbl[n] = (rule >> n) & 1;
  return tbl;
}

// ---------------------------------------------------------------- initial states
function initRandom(w, rnd, fill) {
  const s = new Uint8Array(w);
  for (let i = 0; i < w; i++) s[i] = rnd() < fill ? 1 : 0;
  return s;
}

// Rule 110 background "ether" (periodic) + a few gliders on top.
function initGlider110(w, rnd) {
  let s = initRandom(w, rnd, 0.5);
  // quiescent-ish background: mostly the 110 ether phase
  const ether = [0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0];
  for (let i = 0; i < w; i++) s[i] = ether[i % ether.length];
  // inject left/right moving gliders (the standard 110 glider triplets)
  const gliders = [
    [1, 1, 1, 0, 1], [0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1],
    [0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0],
  ];
  const anchors = [
    Math.floor(w * 0.25), Math.floor(w * 0.5), Math.floor(w * 0.75),
  ];
  for (let g = 0; g < gliders.length; g++) {
    const p = gliders[g], a = anchors[g];
    for (let j = 0; j < p.length; j++) s[(a + j) % w] ^= p[j];
  }
  return s;
}

// ---------------------------------------------------------------- evolution
// Returns rows padded to (W+2P) so block lookups need no modulo in the hot loop.
function evolve(rule, init, steps, P) {
  const tbl = makeRuleTable(rule);
  const Wp = W + 2 * P;
  const rows = new Uint8Array((steps + 1) * Wp);
  let cur = Uint8Array.from(init);
  const writeRow = (t, state) => {
    const off = t * Wp;
    for (let i = 0; i < W; i++) rows[off + P + i] = state[i];
    for (let i = 0; i < P; i++) {
      rows[off + i] = state[(W - P + i) % W];
      rows[off + P + W + i] = state[i % W];
    }
  };
  writeRow(0, cur);
  const nxt = new Uint8Array(W);
  // the already-padded row for time t is exactly what the rule needs:
  // index P+i in the padded row == cell i, so the triple (i-1,i,i+1) is at (P+i-1 .. P+i+1)
  for (let t = 0; t < steps; t++) {
    const base = t * Wp;
    for (let i = 0; i < W; i++) {
      const idx = (rows[base + P + i - 1] << 2) | (rows[base + P + i] << 1) | rows[base + P + i + 1];
      nxt[i] = tbl[idx];
    }
    cur.set(nxt);
    writeRow(t + 1, cur);
  }
  return rows;
}

// ---------------------------------------------------------------- analysis
// Build log2 lookup: logLUT[c] = log2(c/n), 1 <= c <= n. Table size n ~ 1e6 floats = 8MB.
let LOG_N = -1, LOG_LUT = null;
function ensureLogLUT(n) {
  if (LOG_N === n && LOG_LUT) return;
  const t = new Float64Array(n + 1);
  const invN = 1 / n;
  for (let c = 1; c <= n; c++) t[c] = Math.log2(c * invN);
  LOG_N = n; LOG_LUT = t;
}

// H from a raw count histogram. Uses the LUT so there is no transcendental in the hot loop.
function entropyCounts(arr, S, n) {
  const t = LOG_LUT;
  let h = 0, nz = 0;
  for (let i = 0; i < S; i++) {
    const c = arr[i];
    if (c > 0) { h -= (c / n) * t[c]; nz++; }
  }
  return h + (nz - 1) / (2 * n * Math.LN2);
}

// H from a probability distribution (already normalised to sum 1)
function entropyProbs(p, S) {
  let h = 0;
  for (let i = 0; i < S; i++) {
    const v = p[i];
    if (v > 0) h -= v * Math.log2(v);
  }
  return h < 0 ? 0 : h;
}

// Histogram of (L-cell block) -> (next L-cell block), used as the memoryless baseline.
function buildJoint1(rows, L) {
  const P = 8;
  const Wp = W + 2 * P;
  const SB = 1 << L;
  const j = new Int32Array(SB * SB);
  for (let t = 0; t < T; t++) {
    const baseA = t * Wp;
    const baseB = (t + 1) * Wp;
    for (let p = 0; p < W; p++) {
      let w1 = 0, wn = 0;
      for (let i = 0; i < L; i++) {
        w1 = (w1 << 1) | rows[baseA + p + 4 + i];
        wn = (wn << 1) | rows[baseB + p + 4 + i];
      }
      j[w1 * SB + wn]++;
    }
  }
  return j;
}

function analyze(rows, L, k, joint1) {
  const P = L * k + 4;
  const Wp = W + 2 * P;
  const n = W * T;
  const S = 1 << (L * k);
  const SB = 1 << L;
  ensureLogLUT(n);

  const joint = new Int32Array(S * SB);
  const enval = new Int32Array(S);

  // nested layout: enval[w*SB + e] counts history w with environment symbol e
  const b = new Uint8Array(L);
  const bNext = new Uint8Array(L);
  const bEnv = new Uint8Array(L);

  for (let t = 0; t < T; t++) {
    const baseB = (t + 1) * Wp;
    for (let p = 0; p < W; p++) {
      for (let i = 0; i < L; i++) {
        b[i] = 0;
        bNext[i] = rows[baseB + p + P + i];
        bEnv[i] = rows[baseB + ((p + D) % W) + P + i];
      }
      for (let j = k - 1; j >= 0; j--) {
        const base = (t + j) * Wp;
        for (let i = 0; i < L; i++) b[i] = (b[i] << 1) | rows[base + p + P + i];
      }
      let wh = 0, wn = 0, we = 0;
      for (let i = 0; i < L; i++) {
        wh = (wh << 1) | b[i];
        wn = (wn << 1) | bNext[i];
        we = (we << 1) | bEnv[i];
      }
      joint[wh * SB + wn]++;
      enval[wh * SB + we]++;
    }
  }

  // scratch buffer reused for conditional distributions (single allocation, size 2^L)
  const scratch = new Float64Array(SB);

  // ---- H(next | history), p(next), and the history marginal
  const pn = new Float64Array(SB);
  const ph = new Float64Array(S);
  let Hhn = 0;
  for (let w = 0; w < S; w++) {
    const o = w * SB;
    let tot = 0;
    for (let x = 0; x < SB; x++) tot += joint[o + x];
    if (tot === 0) continue;
    const inv = 1 / tot;
    const pw = tot / n;
    ph[w] = pw;
    for (let x = 0; x < SB; x++) {
      const c = joint[o + x];
      pn[x] += c / n;
      scratch[x] = c * inv;
    }
    Hhn += pw * entropyProbs(scratch, SB);
  }
  const Hh = entropyCounts(ph, S, n);
  const Hn = entropyProbs(pn, SB);

  // ---- H(next | single-step history)
  let Hn_given_1 = 0;
  for (let w = 0; w < SB; w++) {
    const o = w * SB;
    let tot = 0;
    for (let x = 0; x < SB; x++) tot += joint1[o + x];
    if (tot === 0) continue;
    const inv = 1 / tot;
    for (let x = 0; x < SB; x++) scratch[x] = joint1[o + x] * inv;
    Hn_given_1 += (tot / n) * entropyProbs(scratch, SB);
  }

  // ---- H(environment | history) and p(environment)
  const pe = new Float64Array(SB);
  let Hhe = 0;
  for (let w = 0; w < S; w++) {
    const o = w * SB;
    let tot = 0;
    for (let x = 0; x < SB; x++) tot += enval[o + x];
    if (tot === 0) continue;
    const inv = 1 / tot;
    for (let x = 0; x < SB; x++) {
      const c = enval[o + x];
      pe[x] += c / n;
      scratch[x] = c * inv;
    }
    Hhe += ph[w] * entropyProbs(scratch, SB);
  }
  const He = entropyProbs(pe, SB);

  const Hmax = L * k;
  const Iself = Math.max(0, Hh + Hn - Hhn);
  const Icross = Math.max(0, Hh + He - Hhe);
  const closure = Hn > 1e-9
    ? Math.max(0, Hn_given_1 - Hhn) / Hn
    : 0;
  const complexityRatio = Hmax > 0 ? Hn / Hmax : 0;
  const selfDensity = Iself / L;
  const advDensity = (Iself - Icross) / L;
  const complexityClamp = Math.max(0, Math.min(1, complexityRatio / 0.25));
  const closureClamp = Math.max(0, Math.min(1, closure / 0.10));
  const observerScore = selfDensity * complexityClamp * closureClamp;

  return {
    Hh, Hn, He, Hhn, Iself, Icross, closure, complexityRatio,
    selfDensity, advDensity, observerScore, S,
  };
}

// ---------------------------------------------------------------- incremental flush
// Written after every (rule, seed) pair so a killed run still leaves usable data.
const fs = require('fs');
const path = require('path');
const outDir = __dirname;
const CSV_HEADER = 'rule,mode,seed,L,k,density,Hh,Hn,He,Hhn,Iself,Icross,closure,complexityRatio,selfDensity,advDensity,observerScore';
function flushResults() {
  let csv = CSV_HEADER + '\n';
  for (const r of results) {
    csv += [r.rule, r.mode, r.seed, r.L, r.k, r.density, r.Hh, r.Hn, r.He, r.Hhn,
      r.Iself, r.Icross, r.closure, r.complexityRatio, r.selfDensity, r.advDensity,
      r.observerScore].join(',') + '\n';
  }
  fs.writeFileSync(path.join(outDir, 'results.csv'), csv);
  fs.writeFileSync(path.join(outDir, 'results.json'), JSON.stringify(results, null, 1));
}

// ---------------------------------------------------------------- main
const t0 = Date.now();
const results = [];
console.log('=== ECA internal-observer scan ===');
console.log(`grid W=${W} T=${T}  rules=${RULES.join(',')}  L=${LS.join(',')}  k=${KS.join(',')}`);
console.log('');

for (const rule of RULES) {
  for (const seed of SEEDS) {
    const rnd = mulberry32(seed);
    const mode = MODES['R' + rule] || 'random';
    const init = mode === 'glider' ? initGlider110(W, rnd) : initRandom(W, rnd, 0.5);
    const Pmax = Math.max(...LS) * Math.max(...KS) + 4;
    const rows = evolve(rule, init, T, Pmax);
    // density sanity check on the whole run
    let ones = 0;
    for (let i = 0; i < rows.length; i++) ones += rows[i];
    const density = ones / rows.length;

    for (const L of LS) {
      const joint1 = buildJoint1(rows, L);
      for (const k of KS) {
        const r = analyze(rows, L, k, joint1);
        results.push({ rule, seed, mode, L, k, density, ...r });
      }
    }
    console.log(`done  rule=${rule} mode=${mode} seed=${seed} density=${density.toFixed(4)}`);
    flushResults();
  }
}

// ---------------------------------------------------------------- report
function fmt(x, d = 4) { return x.toFixed(d).padStart(9); }

console.log('\n=== per-rule summary (mean over L,k,seeds) ===');
console.log('rule  mode     L  k   Hn     Iself   Icross  closure  selfDens  advDens  score');
for (const rule of RULES) {
  const mode = MODES['R' + rule];
  for (const L of LS) {
    for (const k of KS) {
      const g = results.filter(r => r.rule === rule && r.L === L && r.k === k);
      const avg = (f) => g.reduce((s, r) => s + r[f], 0) / g.length;
      console.log(
        String(rule).padEnd(5), mode.padEnd(8), String(L).padStart(2), String(k).padStart(3),
        fmt(avg('Hn')), fmt(avg('Iself')), fmt(avg('Icross')), fmt(avg('closure')),
        fmt(avg('selfDensity')), fmt(avg('advDensity')), fmt(avg('observerScore'))
      );
    }
  }
}

console.log('\n=== top 15 candidates by observerScore ===');
console.log('rule  mode     L  k   seed    Hn     Iself   closure  selfDens  score');
const top = [...results].sort((a, b) => b.observerScore - a.observerScore).slice(0, 15);
for (const r of top) {
  console.log(
    String(r.rule).padEnd(5), r.mode.padEnd(8), String(r.L).padStart(2), String(r.k).padStart(3),
    String(r.seed).padStart(8), fmt(r.Hn), fmt(r.Iself), fmt(r.closure),
    fmt(r.selfDensity), fmt(r.observerScore)
  );
}

console.log('\n=== does memory length k change closure? (mean over L,seeds) ===');
for (const rule of RULES) {
  const line = [];
  for (const k of KS) {
    const g = results.filter(r => r.rule === rule && r.k === k);
    line.push(`k=${k}:${(g.reduce((s, r) => s + r.closure, 0) / g.length).toFixed(4)}`);
  }
  console.log(`rule ${String(rule).padEnd(4)} ${line.join('  ')}`);
}

console.log('\n=== asymmetry check: self vs cross predictive density ===');
for (const rule of RULES) {
  const g = results.filter(r => r.rule === rule);
  const s = g.reduce((a, r) => a + r.selfDensity, 0) / g.length;
  const c = g.reduce((a, r) => a + r.Icross / r.L, 0) / g.length;
  console.log(`rule ${String(rule).padEnd(4)} selfDens=${s.toFixed(5)} crossDens=${c.toFixed(5)} ratio=${(s / c).toFixed(4)}`);
}

console.log(`\nelapsed ${((Date.now() - t0) / 1000).toFixed(1)}s`);

// ---------------------------------------------------------------- dump
flushResults();
console.log(`wrote ${path.join(outDir, 'results.csv')} and results.json`);
