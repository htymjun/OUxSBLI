# CLAUDE.md — 1D_solver/microbenchmark

WENO5-Z reconstruction isolated from the solver, so the FP32/FP64
simultaneous-execution question can be measured without SLAU flux,
Runge-Kutta, or primitive-decode noise. Solver-wide concerns live in
`1D_solver/CLAUDE.md`.

| file | role |
|---|---|
| `weno_micro.f90` | kernels + driver; **self-timing** via CUDA events |
| `run_weno_nsys.sh` | sweep driver, writes `summary.csv` and a ranked table |
| `parse_weno_nsys.py` | nsys CSV/SQLite parser, only used with `--nsys` |
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

## The two mode families

The question the whole directory exists to answer is *which half of WENO goes
to FP32*. The two families answer it differently:

| family | what goes to FP32 | note |
|---|---|---|
| `var_*` | the WENO-Z **weights** (smoothness indicators + ratios) for selected variables | `var_fp32_*` demotes all of ρ/u/p; `var_rho64`/`u64`/`p64` keep one in FP64 |
| `weight_poly32_*` | the **candidate polynomials**, weights stay FP64 | the split originally tried |

The weights carry ~96% of the FP64-pipe work and every division, so the two
families are not close: demoting the weights is worth ~24× on Ada, demoting the
polynomials ~1.1×. See `../report/rtx4060_weno_micro_split.md`.

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
  prints `time_min_us` / `time_avg_us`. No profiler is needed — and that is
  deliberate: `ncu` and `nsys` have both hung repeatedly on this hardware while
  the unprofiled binary ran fine. `--nsys` adds a timeline; if it fails the row
  still carries valid timing.
- **Correctness gate on every row.** The output buffer is poisoned before
  launch, the checksum covers the *whole* array, and a non-finite value is an
  `error stop`, not a printed warning. `summary.csv` records `checksum_all` and
  `nonfinite`; a nonzero `nonfinite` invalidates the timing.

## Two traps that each cost a debugging session

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
