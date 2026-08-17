# CLAUDE.md — 1D_solver/microbenchmark

WENO5-Z reconstruction isolated from the solver, so the FP32/FP64
simultaneous-execution question can be measured without SLAU flux,
Runge-Kutta, or primitive-decode noise. Solver-wide concerns live in
`1D_solver/CLAUDE.md`.

| file | role |
|---|---|
| `weno_micro.f90` | kernels + driver; **self-timing** via CUDA events |
| `run_weno_nsys.sh` | sweep driver, writes `summary.csv` and a ranked table (self-timed only — its per-row nsys capture and `parse_weno_nsys.py` were removed 2026-08-17; the name survives for the historical findings docs) |
| `CMakeLists.txt` | standalone `weno_micro` target, own `CASE_GPU_CC` handling |

Findings: `../report/rtx4060_weno_micro_split.md`.

## Run it

```bash
bash 1D_solver/microbenchmark/run_weno_nsys.sh \
  --out 1D_solver/report/weno_micro_a100 --gpu-cc 80 --repeat 3
```

Or drive the binary directly — `<mode> <nx> <nrepeat> <nlaunch>`:

```bash
cd 1D_solver/microbenchmark
cmake -B build -DCMAKE_Fortran_COMPILER=nvfortran -DCASE_GPU_CC=89
cmake --build build -j
./build/weno_micro var_fp32_seq 4194304 1 10
```

**Always pass `--gpu-cc` / `-DCASE_GPU_CC` explicitly.** It defaults to
`ccnative`, and a stale cached value builds for the wrong architecture, which
fails at *launch*, not at build.

## Stencil width: the `7` / `9` mode suffix

Unsuffixed modes are WENO5-Z; `…7` is WENO7-Z and `…9` is WENO9-Z. Only the
five **core co-issue modes** have wider versions — `var_fp64_seq`,
`var_fp32_seq`, `weight_poly32_seq`, `weight_poly32_serial`,
`weight_poly32_warp` — plus the balanced-split additions: `w32_poly64_*` and
`w64_only_*` exist at all three widths, while `w32mix*` and every `var{3,4,5}_*`
mode are WENO9-only (that is where the A100 signal lives; WENO5 sits near the
memory/issue floor). The old shared-memory/barrier layout variants stay
WENO5-only, since their ranking is already settled.

The WENO7/9 device routines (and their `_right` twins) are **generated**:

```bash
python3 ../report/check_weno_order.py --gen-micro 4   # WENO7-Z split form
python3 ../report/check_weno_order.py --gen-micro 5   # WENO9-Z
```

`--gen-micro` emits the *split* form (separate normalised weights and candidate
polynomials) rather than the fused `result(vf)` form `--gen` emits for the
solver — folding the weight normalisation into a final division the way the
solver does would move FP64 work across the very boundary this benchmark
measures. `micro_split_test()` in that script order-verifies the split form on
**both biases** at 5.00 / 7.02 / 8.95.

**WENO5-Z stays hand-written on purpose.** Its betas use the textbook
sum-of-two-squares form, which is a WENO5 special case — the published WENO7/9
betas are dense quadratic forms. Regenerating r=3 would swap 2 squares for 6
products and change WENO5's instruction count, invalidating every prior WENO5
timing.

### What widening actually showed

Measured at `nx=1048576`, RTX 4060:

| | WENO5 | WENO7 | WENO9 |
|---|---:|---:|---:|
| `var_fp64_seq` (FP64 baseline) | 5411 µs | 9807 µs (1.81×) | 16080 µs (2.97×) |
| `var_fp32_seq` (weights → FP32) | 432 µs | 529 µs | 715 µs |
| **demotion speedup** | **12.5×** | **18.5×** | **22.5×** |
| `weight_poly32_serial` | 5389 | 9083 | 14188 |
| `weight_poly32_warp` | 5373 | 9078 | 14809 |
| **serial → warp (co-issue)** | **1.003×** | **1.001×** | **0.958×** |

`ncu` confirms the mechanism: `smsp__inst_executed_pipe_fp64.sum` scales
1.00 / 1.81 / 2.76 and time scales identically, matching the 1.00 / 1.82 / 2.79
predicted from the solver's per-face SASS counts.

**Widening the stencil does not rescue warp co-issue.** The controlled
`serial → warp` pair is worth nothing at WENO5 and WENO7, and at WENO9 the warp
split is **4.4% slower** — its shared memory reaches 46 KB, which caps occupancy
at one block per SM. What *does* grow with the stencil is the precision
demotion: 12.5× → 18.5× → 22.5×. This is the same conclusion the solver reached
by a different route — relocating work is free, removing it is what pays.

## The mode families

The question the whole directory exists to answer is *which half of WENO goes
to FP32*. The families answer it differently:

| family | what goes to FP32 | note |
|---|---|---|
| `var_*` | the WENO-Z **weights** (smoothness indicators + ratios) for selected variables | `var_fp32_*` demotes all of ρ/u/p; `var_rho64`/`u64`/`p64` keep one in FP64 |
| `weight_poly32_*` | the **candidate polynomials**, weights stay FP64 | the split originally tried |
| `w32_poly64_*` | the **weights** (~88% of the arithmetic), polynomials + combine stay FP64 | the reverse split — the A100 pipe-balance family (1:2 wants N_FP32 = 2·N_FP64). Order-verified: `check_weno_order.py --w32-split` shows the design order down to an ~2e-8 L1 floor, because weight errors multiply O(h^r) candidate differences |
| `w32mix{h,1}_poly64_seq9` | all but the first 1 / 2 of the 6 var-bias **weight calls** | balance tuning between `weight_poly32` (0 of 6) and `w32_poly64` (6 of 6) |
| `var{3,4,5}_{fp64,fp32,k*}_seq9` / `_k*_warp9` | all but the first k **variables**, at nv = 3/4/5 (ρ,u,p / ρ,u,v,p / ρ,u,v,w,p) | variable-granular balance; nv=5 models the 3D solver. `_warp9` needs **no shared memory and no barrier** — each role owns complete reconstructions and stores its own columns. One parameterized kernel pair (`weno_varsplit_{seq,warp}9(n,nrepeat,nv,k,…)`); `var3_fp64/fp32_seq9` reproduce `var_fp64/fp32_seq9` bit-for-bit, quantifying the runtime-trip-count codegen delta |
| `w64_only_seq{,7,9}` | nothing — the FP64 **weights stream alone**, no polynomials | ablation: `t(weight_poly32_seqX) − t(w64_only_seqX)` isolates whether the FP32 half of `weight_poly32_seqX` rides free (ILP co-issue) or serializes. Output is a fixed weight combination, so its checksum is not comparable to anything else |

The weights carry ~96% of the FP64-pipe work and every division, so `var_*` vs
`weight_poly32_*` are not close: demoting the weights is worth ~24× on Ada,
demoting the polynomials ~1.1×. See `../report/rtx4060_weno_micro_split.md`.
On A100 the interesting question inverts — time tracks the FP64-pipe count, so
the FP32 work already hides under the FP64 stream, and the balanced splits
exist to load *both* pipes near their ratio.

The `w32_poly64_serial/warp` pair fixes two known unfairnesses of the older
pair: only the real(4) weights cross shared memory (9216/12288/15360 B at
5/7/9, vs 27648–46080 B), and the checksum acc comes from registers instead of
re-reading `out()` from global.

Within `weight_poly32_*`, the suffixes are a controlled ladder:

| suffix | what it changes |
|---|---|
| `_seq` | one thread does both pieces; 128 threads, no shared memory |
| `_serial` | 256 threads + shared memory + barrier, **both** pieces still on the lower warp role |
| `_warp` | identical layout to `_serial`, but roles split FP64/FP32 |
| `_serial_oncebar` / `_warp_oncebar` | barrier hoisted out of the inner loop |
| `_wsmem_*` | only weights in shared memory; FP32 polynomials stay in registers |
| `_wsmem_tile2_*` | two faces per role, amortising one barrier |
| `_halfwarp_serial` / `_halfwarp_shfl` | shared-memory-free `__shfl_xor` diagnostic |

**`_serial` → `_warp` is the only clean isolation of co-issue**: identical
thread count, shared memory and barrier count, only the role assignment
differs. Compare against `_seq` separately — it has a different thread budget,
so it is not a controlled pair.

`_halfwarp_*` splits lanes *within* one warp, so it is lane-divergent and is a
diagnostic, not the cross-warp strategy.

## Harness rules

- **`nrepeat` is pinned to 1** and `run_weno_nsys.sh` rejects `--nrepeat` /
  `--nrepeat-list` outright. At `nrepeat>1` nvfortran hoists the inner loop for
  any mode that writes only private registers (`weight_poly32_seq`,
  `halfwarp_serial`: ~2–4 µs per extra repeat) but not for modes that touch
  shared memory or cross a barrier (~700–1400 µs). The rows stop being
  comparable, silently.
- **Rounds are interleaved** (round-major, not mode-major) because GPU clocks
  cannot be locked without root here. Compare `time_min_us` across rounds.
- **The binary times itself** with CUDA events over `nlaunch` launches and
  prints `time_min_us` / `time_avg_us`. No profiler is involved — and that is
  deliberate: `ncu` and `nsys` have both hung repeatedly on this hardware while
  the unprofiled binary ran fine. The sweep script's per-row `--nsys` capture
  was removed 2026-08-17 (it now rejects the flag); profile a single mode by
  hand if a timeline is ever needed.
- **Correctness gate on every row.** The output buffer is poisoned before
  launch, the checksum covers the *whole* array, and a non-finite value is an
  `error stop`, not a printed warning. `summary.csv` records `checksum_all` and
  `nonfinite`; a nonzero `nonfinite` invalidates the timing.
- **`checksum_all` is a launch/NaN guard, NOT a correctness check.**
  `init_input` builds a piecewise-linear field (a ramp plus one step), and every
  WENO candidate is exact on linear data, so the reconstruction is insensitive
  to the weights almost everywhere. Measured: giving `v^+` the *wrong* optimal
  weights — the third-order bug below — moves the result by at most **4.4e-16**
  (one ulp) and the whole-array sum by 3e-11 out of 1.27e6, far under the
  printed resolution. Anything that changes *which* WENO you are computing must
  be checked with `../report/check_weno_order.py`, not with this checksum.
  The FP32-weights split has its own gate there: `--w32-split` runs the
  float32-weights × float64-polynomials numpy twin. Checksums are comparable
  only within one family at one width — 7/9 halo offsets, var4/var5 column
  counts and `w64_only`'s weights-only output all shift the sum legitimately,
  while `w32_poly64_{seq,serial,warp}` at one width and `var{N}_k*_{seq,warp}9`
  at one (nv,k) must match **bit-for-bit**.
- **`maxregcount:128`**, raised from the solver's 96. At 96 the WENO9-Z kernels
  pinned to the cap and spilled 64–84 B, injecting LDL/STL into the instruction
  stream being measured. At 128 they settle at their natural 114 registers with
  zero spill, and WENO5/7 are unaffected (61/64/60 either way), so prior WENO5
  timings stay comparable. `w32_poly64_{serial,warp}9` sit at exactly 128 with
  **zero spill** (their `p(30)` real(8) polynomial array lives across the
  barrier); re-check ptxinfo after any change that adds live values there.

## Three traps that each cost a debugging session

- **Right-biased reconstruction needs MIRRORED optimal weights.** Until
  2026-08-17 `weights64`/`weights32` used `d = 0.1/0.6/0.3` for *both* sides
  while `poly64_right` used the mirrored candidate ordering, making `v^+` third
  order instead of fifth. Instruction counts are identical, so **every timing
  result from before the fix still stands** — but the values were wrong, and the
  checksum could not see it (above). There are now explicit
  `weights{5,7,9}_{64,32}_{left,right}` routines; the `_right` set carries `d`
  reversed. The same bug existed in the solver's `weno5z_right`.


- **Stale `CASE_GPU_CC` → `cudaErrorInvalidPtx` (218).** A cached CC that does
  not match the GPU builds cleanly and then fails at every launch, because the
  `lto` device code cannot be JIT'd across architectures. This was invisible
  until the driver grew a `cudaGetLastError()` check — `cudaDeviceSynchronize()`
  alone returns 0 when the launch itself was rejected, so the untouched buffer
  was reported as a valid `checksum=0.0`. `CMakeLists.txt` now echoes the target
  architecture and `run_weno_nsys.sh` copies it into `progress.log`.
- **nvfortran 24.7 miscompiles nested `value` dummies.** A device routine with
  `intent(in), value` dummies, called from *another* device routine that also
  has `value` dummies, returns NaN at `-O2`/`-fast`; `-O0` is correct and
  `-Kieee`/`-Mnofprelaxed`/`-Mnoflushz`/`-Mnovect` do not help. One nesting
  level is fine, which is why `weight_poly32_*` always worked while every
  `var_*` mode returned NaN — including `var_fp64_seq`, which contains no FP32
  at all. **Do not put `value` on the real(8) dummies in `weno_micro.f90`.**
  (`integer, value` on the `attributes(global)` kernels is required and fine.)
  That NaN was for a while misattributed to FP32 overflow in the WENO-Z
  `tau5/(β+ε)` ratio; the `ratio_cap` guard in `weights32` is still correct and
  still needed, but it was not the cause.

## Keeping the arithmetic in step with the solver

`weno_micro.f90` carries its own copy of the WENO-Z math, so it can drift from
`src/calc_weno.f90`. It currently mirrors the division reduction described in
`../report/rtx4060_weno_division_reduction.md`:

- `weights64` — batched ratio inversion **and** one reciprocal for `1/s`: 6 → 2
  divisions.
- `weights32` — normalisation batched only: 6 → 4. The ratio batch would
  underflow, since in FP32 `eps³ = 1e-60` flushes to zero and `1/0 → Inf` makes
  the ratios `0*Inf = NaN` on a smooth plateau.
- `poly64_*` / `poly32_*` — the `1/6` folded into the coefficients: 3 → 0.

The weights/poly **split is preserved on purpose** — folding the `1/6` into the
final normalisation the way the solver does would move FP64 work across the
very boundary being measured.
