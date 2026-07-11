# mm2-fast optimizations + Apple Silicon tuning

This fork ports the [mm2-fast](https://github.com/bwa-mem2/mm2-fast) optimizations
(Kalikar et al., *Nature Comput. Sci.* 2022) onto minimap2 v2.31 and adds tuning
for Apple Silicon / AArch64.

**Byte-identity (important):**
- The **Apple Silicon tuning** (`arm_tune=1` + modern sse2neon) and the **AVX2/AVX-512
  alignment** kernels are byte-identical to stock v2.31, verified by `cmp` on PAF and
  SAM, single/multi-thread, across presets and on real human HiFi reads → GRCh38 (the
  AVX kernels additionally validated over 120k random cases).
- The **vectorized chaining** (`vchain`/`neonchain`) is byte-identical on small/simple
  inputs but **NOT on all real data**: it must drop the `max_skip` heuristic (a
  loop-carried dependency), computing the `max_skip=∞` result. On real HiFi this changes
  a small fraction of records (~0.6% of lines: mostly secondary-chain scores). This is
  the same limitation as upstream mm2-fast. Kept off by default; enable only if you
  accept mm2-fast's `max_skip=∞` chaining semantics.

## Build options

| Command | What it enables | Target |
|---|---|---|
| `make` | Stock minimap2 (x86 SSE dispatch) | x86 |
| `make aarch64=1` | ARM/NEON build (modern sse2neon shim) | Apple/AArch64 |
| `make aarch64=1 arm_tune=1` | + `-O3 -mcpu=native -flto` | Apple/AArch64 |
| `make aarch64=1 arm_tune=1 vchain=1` | + vectorized DP chaining (SIMDe→NEON) | Apple/AArch64 |
| `make avx=1` | + AVX2 two-piece extension kernel | x86 (AVX2) |
| `make avx512=1` | + AVX2 **and** AVX-512 kernels | x86, **gcc/icc only** |
| `make vchain=1` | + vectorized DP chaining (native AVX2/AVX-512) | x86 |

Notes:
- `vchain=1` needs a C++ compiler (`$(CXX)`, default `c++`) and links the C++
  runtime. On ARM it uses **SIMDe** (`lib/simde`) to lower AVX2 intrinsics to NEON.
- `avx512=1` requires **gcc or icc**: the mm2-fast AVX-512 kernel uses non-constant
  shuffle immediates that clang rejects. `avx=1` (AVX2) builds cleanly on clang and gcc.
- `arm_tune=1` uses `-mcpu=native` (ties the binary to the build host's uarch).

## What each optimization does, and where it actually helps

### 1. Vectorized DP chaining (mm2-fast, `vchain=1`)
The chaining inner loop is replaced by the mm2-fast Structure-of-Arrays,
`max_skip`-free vectorized kernel (`contrib/parallel_chaining_v2_22.h`, wrapped by
`mm2fast_chain.cpp`). Only engaged for the regime where it is provably identical to
scalar `comput_sc()`: single-segment, non-cDNA, equal ref/query gap limits;
everything else falls back to the scalar loop.

- **x86**: native AVX2 (8-wide) / AVX-512 (16-wide) — a real speedup (paper: up to 3.1x
  chaining vs default).
- **Apple M-series**: **regresses**. Left **off by default on ARM**. Measured on real
  HiFi (chr1 index, 5k reads, single-thread CPU time, chaining-dominated PAF):

  | chaining kernel | time | vs scalar |
  |---|---|---|
  | scalar (`max_skip=25`) | 9.35s | 1.00x |
  | vchain (SIMDe→NEON) | 13.00s | 0.72x |
  | `neonchain` (hand-written native 4-wide NEON) | 12.83s | 0.73x |

  Root cause: the vectorized kernel must drop `max_skip` (a loop-carried dependency),
  so it does the full `max_skip=∞` work (~2.6x more than the scalar default). NEON's
  4-wide int32 width cannot recoup that. A hand-written native NEON kernel
  (`neon_chain.c`, `neonchain=1`) was implemented to remove SIMDe's 8→4 AVX2 lowering
  overhead — it is byte-identical to the SIMDe kernel (verified) but only ~1% faster
  than it (clang already lowers SIMDe well); the bottleneck is algorithmic, not SIMD
  width. **Conclusion: vectorized chaining does not help on Apple Silicon.** Both
  kernels are kept as opt-in for x86 (where wide AVX-512 does win) and for reference.

### 2. AVX2 / AVX-512 two-piece extension alignment (mm2-fast, `avx=1`/`avx512=1`)
`ksw2_extd2_avx.c` widens the default dual-affine extension DP from 128-bit SSE to
256/512-bit. Runtime CPU dispatch (`ksw2_dispatch.c`) routes to it. Reconciled with
the v2.31 SSE kernel (the `mte_q = r - en0` fix) and validated **byte-identical to the
SSE kernel over 120,000 random cases** (portably, via SIMDe emulation).

- **x86**: real speedup (paper: 1.8–2.2x on the alignment stage with AVX-512).
- **Apple M-series**: no benefit (NEON caps at 128-bit = the existing SSE width), so it
  is not built on ARM.

### 3. Apple Silicon tuning (`arm_tune=1` + modern sse2neon)
The dominant cost for map-ont/map-hifi is the ksw2 alignment (NEON) path. Two changes
help it on Apple Silicon:
- Replaced the legacy vendored sse2neon shim with the modern
  [DLTcollab sse2neon](https://github.com/DLTcollab/sse2neon) (`sse2neon/sse2neon.h`;
  the old shim is kept as `sse2neon/emmintrin_legacy.h`).
- `arm_tune=1` adds `-O3 -mcpu=native -flto`.

## Measured results (Apple M4 Max)

Byte-identical throughout. Real human HiFi (2153 GIAB HG002 reads → GRCh38):

| Build | map-hifi mapping time | speedup |
|---|---|---|
| stock `aarch64=1` (`-O2`, legacy shim) | baseline | 1.00x |
| `aarch64=1 arm_tune=1` (modern shim, `-O3 -mcpu=native -flto`) | faster | ~1.06x |

(Small-genome presets showed up to ~1.12–1.15x; on real HiFi the isolated mapping
speedup is ~1.06x. This is a genuine but modest Apple-Silicon win — far from mm2-fast's
1.5–1.8x, which came from x86 AVX-512 that has no NEON equivalent.)
