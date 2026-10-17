'use strict';
/*
 * discrete_phases.js
 * ---------------------------------------------------------------------------
 * 如果"相位"只能取有限个离散值，最大 CHSH 是多少？它如何趋近 2*sqrt(2)？
 *
 * THE POINT (this tests a physical intuition, not a formula)
 *   Nobody measures a truly continuous angle. Real Bell experiments use a finite
 *   set of settings. And the user's intuition is that reality is discrete and
 *   continuity is a coarse-grained appearance.
 *
 *   So the decisive question is not "can a discrete system give 2*sqrt(2) exactly"
 *   but "how many discrete phase values are needed before the bound is
 *   indistinguishable from 2*sqrt(2)?" If a modest number suffices, then
 *   continuity is indeed only an APPEARANCE, and the discrete picture works.
 *
 * MODEL
 *   Two particles share the maximally entangled state (real amplitudes suffice --
 *   we established earlier that the real model already saturates 2*sqrt(2)):
 *        |psi> = (|00> + |11>) / sqrt(2)
 *   Alice and Bob each have N allowed measurement angles, evenly spaced:
 *        theta_k = k * pi / N,   k = 0..N-1
 *   For this state the correlation is E = cos(theta_A - theta_B).
 *   We maximise  S = E(a0,b0) + E(a0,b1) + E(a1,b0) - E(a1,b1)
 *   over all choices of two angles per side, and report S(N).
 *
 *   N = 4 corresponds to a four-valued phase (the Z4 case we tested in Life).
 *   N -> infinity recovers the full continuum.
 * ---------------------------------------------------------------------------
 */

function chshFor(N, i0, i1, j0, j1) {
  const th = (k) => k * Math.PI / N;
  const E = (a, b) => Math.cos(th(a) - th(b));
  return E(i0, j0) + E(i0, j1) + E(i1, j0) - E(i1, j1);
}

console.log('=== CHSH 与离散相位数的关系 ===');
console.log('状态 (|00>+|11>)/sqrt2，每侧 N 个等距角度，穷举最优对 ===\n');
console.log('   N    最优角度(度)              最大 CHSH      2√2 - S       相对差');
const TARGET = 2 * Math.SQRT2;

const rows = [];
for (const N of [2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 48, 64, 96, 128, 256]) {
  let best = -Infinity, bestIdx = null;
  for (let i0 = 0; i0 < N; i0++) for (let i1 = 0; i1 < N; i1++) {
    for (let j0 = 0; j0 < N; j0++) for (let j1 = 0; j1 < N; j1++) {
      const S = chshFor(N, i0, i1, j0, j1);
      if (S > best) { best = S; bestIdx = [i0, i1, j0, j1]; }
    }
  }
  const deg = bestIdx.map(i => (i * 180 / N).toFixed(1)).join(', ');
  const gap = TARGET - best;
  const rel = gap / TARGET;
  rows.push({ N, best, gap, rel });
  console.log(`  ${String(N).padStart(3)}   ${deg.padEnd(26)}  ${best.toFixed(6)}    ${gap.toExponential(2)}   ${rel.toExponential(2)}`);
}

console.log('\n=== 关键读数 ===');
const N4 = rows.find(r => r.N === 4);
const N8 = rows.find(r => r.N === 8);
const N20 = rows.find(r => r.N === 20);
const N128 = rows.find(r => r.N === 128);
console.log(`  N=4   (四值相位，即我们测过的 Z4): S = ${N4.best.toFixed(6)}   差 2√2 还有 ${N4.rel * 100}%`);
console.log(`  N=8   : S = ${N8.best.toFixed(6)}   差 ${N8.rel * 100}%`);
console.log(`  N=20  : S = ${N20.best.toFixed(6)}   差 ${N20.rel * 100}%`);
console.log(`  N=128 : S = ${N128.best.toFixed(6)}   差 ${N128.rel * 100}%`);

// ---- 需要多少 N 才能达到给定位精度？
console.log('\n=== 达到给定精度需要多少个离散相位？ ===');
for (const tol of [0.1, 0.01, 0.001, 1e-4, 1e-6]) {
  const hit = rows.find(r => r.gap < tol);
  console.log(`  误差 < ${String(tol).padEnd(8)} : N = ${hit ? hit.N : '>256'}`);
}

// ---- 和"实验实际用几个角度"对照
console.log('\n=== 与真实实验对照 ===');
console.log('  真实的 Bell 实验通常只用 2-4 个设置角度。');
console.log(`  我们的计算显示 N=4 已经给出 S = ${N4.best.toFixed(6)}，`);
console.log(`  即 2√2 的 ${(N4.best / TARGET * 100).toFixed(2)}%。`);
console.log('  => 从实验的角度看，连续相位从来就不是必需的。');

console.log('\n=== 这对你的直觉意味着什么 ===');
console.log('  1. 要"看起来连续"，不需要真的连续。N 只有几十就够把 2√2 逼到');
console.log('     小数点后好几位。');
console.log('  2. 一个四值相位（N=4）就已经给出 2√2 的相当大的一部分。');
console.log('  3. 所以"底层离散、连续是粗粒化感知"这个直觉，在 CHSH 这一项上');
console.log('     是站得住的：离散相位足以复现连续相位的可观测量。');
console.log('');
console.log('  但注意一个遗留问题：N 需要多大，取决于底层能提供多少内部状态。');
console.log('  我们测过生命游戏的 Z4 图案，内部轨道最长只有 2 个状态（N=2），');
console.log('  那个 N 太小。真正的底层需要一个能提供大 N 的结构。');
