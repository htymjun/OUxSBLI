# Benchmark E/F/G vs fused dR convective kernels, then remove the old generation

## Context

The convective flux kernels currently exist in **two generations** in `3D_solver/src/`:

- **OLD**: `calc_{keep,slau,roe,hybrid}_{x,y,z}[_in]` write materialized face-flux
  arrays `E(nx-1,ny-2,nz-2,5)`, `F`, `G`; then `calc_dR_from_EFG`
  ([calc_steps.f90:13](3D_solver/src/calc_steps.f90#L13)) derives the flux
  divergence `dR` from them.
- **NEW (fused)**: `calc_{...}_{dx,dy,dz}[_in]` write/accumulate the divergence
  **directly into `dR`**, never materializing E/F/G (warp-shuffle for x,
  atomicAdd for y/z).

The fused generation was added by the recent "kernel split increased the compute
throughput thanks to nice memory access" work, but the old generation was left in
place and is still compiled. The goal: **measure whether fused dR is actually
faster** (ETGV and NSTGV at nx=ny=nz=128), and if so **delete the old
generation** — which also frees the E/F/G allocations.

### What exploration established (verified by direct reads)

- **No shipped 3D case executes the old kernels.** All 8 (DHIT, ETGV, IVST, KHI,
  NSTGV, SBLI, STZ, TBL) take the fused path. Only compile-time references remain.
- Three *buildable* combos still route to the old path and must be guarded before
  deletion: `NS/LES + Roe`, `NS/LES + VISC_ORDER==4`, `COMMZ + NS/LES`.
- The Euler compatibility `#:else` at
  [calc_flux_base.f90.fypp:323-329](3D_solver/src/calc_flux_base.f90.fypp#L323)
  is **already unreachable** — its `#:if SCHEME in ('KEEP','SLAU','Hybrid','Roe')`
  tuple is the entire legal SCHEME domain per `src/mod_constant.f90.fypp`.
- **`3D_solver_mixed` is a real live consumer**: its `CMakeLists.txt:63-64`
  compiles `3D_solver/src/calc_slau_kernel[_internal].f90.fypp` and
  `src/calc_flux_base_mixed.f90.fypp:71-73,134` calls `calc_slau_{x,y,z}_in`
  directly. Per decision below, those specific kernels are **retained**.
  `2D_solver`, `3D_solver_curv`, `3D_solver_fltflt` have their own independent
  copies and are unaffected.
- E+F+G at 128³ = **~230 MiB (~35% of the ~650 MiB total device footprint)**,
  allocated unconditionally at
  [preprocess.f90.fypp:36](3D_solver/src/preprocess.f90.fypp#L36).

### Decisions taken (confirmed with user)

1. **Removal scope**: guard the three unused combos with `#:stop`, then delete the
   old generation. (Not the fuller "convert remaining paths first" option.)
2. **`3D_solver_mixed`**: keep `calc_slau_{x,y,z}[_in]` + their modules for it,
   with an explanatory comment; leave the mixed fork untouched.

## Phase 1 — Benchmark (do this first; Phase 2 is gated on the result)

### 1a. Add a temporary `LEGACY_EFG` A/B toggle

In [calc_flux_base.f90.fypp](3D_solver/src/calc_flux_base.f90.fypp), using the
existing default idiom already at lines 1-4 (`#:if not defined('VISC_ORDER')`),
so **no case `config.fypp` needs editing**:

```
#:if not defined('LEGACY_EFG')
  #:set LEGACY_EFG = False
#:endif
```

- **Euler branch (line ~303)**: append `and not LEGACY_EFG` to the fast-path
  condition. This makes the existing (currently dead) compat branch
  `calc_conv` + `calc_dR_from_EFG` reachable. ETGV is `SCHEME='KEEP'`, so
  neither leg runs a shock sensor → **already a clean convective-only A/B**.
- **NS/LES branch (line ~390)**: append `and not LEGACY_EFG`, and add a new
  `#:elif <same condition> and LEGACY_EFG` branch that isolates **only** the
  convective change:
  1. LES-only `calc_mut` + `set_bc_mut` (identical to fast path)
  2. sensor **identical to fast path** — `calc_Ducros_div` for `VISC_ORDER>2`
  3. old `calc_${scheme}_{x,y,z}[_in]` → E/F/G (reuse the existing
     `args_t`/`args_sensor` fypp vars, same as `calc_conv` builds them)
  4. `calc_dR_from_EFG` — fresh `dR` write
  5. fused `calc_dEv/dFv/dGv${VISC_ORDER}$[_in]` — accumulate onto `dR`

  **Do NOT reuse `calc_conv` here**: it bundles plain `calc_Ducros`
  ([line 158](3D_solver/src/calc_flux_base.f90.fypp#L158)) whereas the fast path
  uses `calc_Ducros_div` ([line 419](3D_solver/src/calc_flux_base.f90.fypp#L419)),
  which would confound the comparison with a different sensor kernel.

  Ordering is safe: `calc_dR_from_EFG` declares `dR` `intent(out)` and writes
  every point of `(nx-2,ny-2,nz-2,5)`, so it is a complete fresh write before the
  viscous accumulators (separate serialized launches on the default stream).

### 1b. Harness

Throwaway Python script **in the scratchpad, not the repo**.

**⚠ CORRECTION — `Case(**kwargs)` CANNOT carry a new fypp flag.** The first
attempt passed `LEGACY_EFG=False/True` as a `Case` kwarg. `Case._setup()`
(`ouxsbli/case.py:110-152`) routes a param to `config.fypp` only if it already
matches a `#:set NAME` line there, **otherwise to `mod_globals.f90`** — where an
unknown Fortran name is silently discarded. Net effect: the flag reached *neither*
file, all four builds compiled the **fused** path, and the whole v1 run was
fused-vs-fused. All v1 numbers were discarded.

Correct approach (v2, `scratchpad/bench2.py`): call `Case._setup()` for the
nx/ny/nz/nt/np patching and env detection, then **append `#:set LEGACY_EFG = …`
to the workdir's `config.fypp` directly**, then invoke cmake/make manually (not
`Case.build()`, which would re-run `_setup()` and wipe the injected line). cmake
needs explicit `-DCMAKE_Fortran_COMPILER=…/nvfortran` and
`-DMPI_Fortran_COMPILER=…/mpif90` or it falls back to GNU `f95` and fails to find
MPI.

**Mandatory guard:** after each build, assert on the generated
`build/calc_flux_base.f90` that the intended generation is actually wired in
(fused ⇒ contains `call calc_<scheme>_dx` and NO `call calc_dR_from_EFG`;
legacy ⇒ the converse). This is what would have caught the v1 no-op immediately,
and no timing number should be trusted without it.

- Separate workdir per (case × leg) so builds don't clobber each other.
- `np=1` → one I/O at the very end; `nt=200` **pinned explicitly** because
  NSTGV's `nt` is otherwise derived from a physical end time
  ([NSTGV/mod_globals.f90:57](3D_solver/NSTGV/mod_globals.f90#L57)) and would
  otherwise vary. Both legs must run an identical step count.
- 2 MPI ranks (matches both cases' `calc.sh`); only rank 0 computes, rank 1 is the
  I/O partner, so there is no GPU contention.
- 128³ ≈ 650 MiB of 8 GiB — no OOM risk.

### 1c. Metrics — two independent measurements

The built-in `calculation time:` print is unusable for this
([3D_solver/src/main.f90:67-80](3D_solver/src/main.f90#L67)): it is `cpu_time`,
truncated to **whole seconds**, and includes allocation + all I/O.

1. **Wall clock**: `time.perf_counter()` around the `mpirun`; 3 repeats;
   report min and median.
2. **Per-kernel GPU time** (the number that actually isolates the change):
   `nsys profile -t cuda -f true -o <tag> mpirun -n 2 ./build/a.out`, then
   `nsys stats --report cuda_gpu_kern_sum`. Compare
   `calc_slau_dx_in+dy_in+dz_in` (fused) vs
   `calc_slau_x_in+y_in+z_in + calc_dr_from_efg` (legacy). Also confirm the
   viscous and step kernels are **unchanged** between legs — that proves the
   isolation actually worked. (Reuses the existing `profile.sh` pattern; prefer
   `nsys` over `ncu --set full`, which replays every kernel and is far slower.)

### 1d. Correctness gate (before trusting any timing)

Both legs must produce matching `data/kinetic_energy.d` and `data/entropy.d` to
~roundoff. If they diverge, the legacy leg is not a fair functional equivalent and
the comparison is void — investigate before proceeding.

### 1e. Decision rule

Proceed to Phase 2 **only if fused dR is faster on both cases by both metrics**.
Report the measured numbers either way, including the case where it is not faster.

### 1f. Report deliverable — `3D_solver/nsys_ncu/benchmark_dr_vs_efg.md`

Write the results to a **committed markdown report**, co-located with the existing
profiling captures since `CLAUDE.md`'s Notice section already directs readers to
`3D_solver/nsys_ncu/` for performance data (easy to move to `docs/` or add to
`mkdocs.yml`'s nav later if it should be published). Contents:

- **Verdict up front** — one line: is fused dR faster, and by how much.
- **Setup**: exact grid (128³), pinned `nt`/`np`, rank count, GPU (RTX 4060 Laptop,
  8 GiB), compiler (`nvfortran` 24.7), and the exact `mpirun` / `nsys` command lines
  used, so the numbers are reproducible.
- **What was compared, precisely**: that the A/B is *convective-only* — same sensor
  kernel, same viscous generation in both legs — and how that was enforced
  (`LEGACY_EFG`), so a future reader does not mistake it for a whole-solver diff.
- **Wall-clock table**: per case × leg, min and median of 3 repeats, plus % delta.
- **Per-kernel GPU-time table** from `nsys stats --report cuda_gpu_kern_sum`:
  the convective kernels in each leg (`calc_slau_dx_in` … vs
  `calc_slau_x_in` … + `calc_dr_from_efg`), their totals, and the delta. Include
  the unchanged viscous/step kernel rows as **evidence the isolation held**.
- **Correctness evidence**: the `kinetic_energy.d` / `entropy.d` agreement from 1d.
- **Memory**: the E/F/G footprint that removal frees (~230 MiB at 128³, ~35%).
- **Caveats**, honestly stated: e.g. 128 vs 129 for periodic TGV, `cpu_time`-based
  in-solver timer not used and why, any run-to-run variance observed, and — if the
  result is mixed or negative — say so plainly rather than shading it positive.

Keep raw artifacts alongside it (the `.nsys-rep` / `nsys stats` CSV output) if they
are a reasonable size; otherwise record the commands to regenerate them.

## Phase 2 — Remove the old generation (gated on Phase 1)

### Pre-flight check — DONE, results below

- **`calc_dR_from_EFG` is referenced only inside `3D_solver/src/`**
  (`calc_time_dev.f90.fypp` ×7, `calc_flux_base.f90.fypp`, and its own definition
  in `calc_steps.f90:13`). Although `3D_solver/src/calc_steps.f90` *is* compiled by
  `3D_solver_fltflt` and `3D_solver_mixed`, neither calls it, and `2D_solver` has
  its own `2D_solver/src/calc_steps.f90`. → **safe to delete.**
- **`3D_solver/src/preprocess.f90.fypp` IS shared** — it is in the `_BASELINE_FYPP`
  list of both `3D_solver_fltflt/CMakeLists.txt:67-71` and
  `3D_solver_mixed/CMakeLists.txt:59-63`, both read from `../../3D_solver/src/`.
  Changing `allocate_device_mem`'s signature therefore reaches both forks.
- **⚠ Both forks are ALREADY STALE against that signature — pre-existing, not
  caused by this work.** The current signature is
  `(..., T, mu, mut, qc2, ux, vy, wz, E, F, G)`, but:
  - `3D_solver_mixed/src/calc_time_dev_mixed.f90.fypp:82-83` calls
    `(..., mut, qc2, E, F, G, inv_dx_r4, inv_dy_r4, inv_dz_r4)` — no `ux,vy,wz`,
    plus three r4 args the signature does not have;
  - `3D_solver_fltflt/src/calc_time_dev_ff.f90.fypp:62,157` calls
    `(..., mut, qc2, E, F, G)` — no `ux,vy,wz`;
  - `preprocess.f90.fypp:39`'s `#:if MIXED_FORK` branch allocates
    `xix_r4/etay_r4/zetaz_r4`, which are **not declared in that subroutine**
    (only in `pre_calc`, line 78).

  These are fallout from the (uncommitted) `ux/vy/wz` div-u refactor that predates
  this task. **Consequence for the plan:** "confirm `3D_solver_mixed` /
  `3D_solver_fltflt` still build" is NOT a usable regression gate — confirm the
  pre-existing failure first, then require only that my change does not make it
  *worse*. Report this to the user rather than silently fixing unrelated forks.

### Guards (reuse the `#:stop` idiom already in `src/mod_constant.f90.fypp`)

Add clear `#:stop` messages for the three combos now unsupported:
`NS/LES + Roe`, `NS/LES + VISC_ORDER==4`, `COMMZ + NS/LES`.

### Delete

- `calc_keep_kernel[_internal].f90.fypp` — `calc_keep_{x,y,z}[_in]`, `calc_keep_z_in_koff`
- `calc_hybrid_kernel[_internal].f90.fypp` — `calc_hybrid_{x,y,z}[_in]`, `calc_hybrid_z_in_koff`
- `calc_roe_kernel[_internal].f90.fypp` — `calc_roe_{x,y,z}[_in]`
  (this **supersedes** last session's bug-#3 fix: the buggy legacy Roe order
  dispatch disappears entirely — note it in the commit message)
- `calc_flux_base.f90.fypp` — `calc_conv`, `calc_conv_EF`, `calc_conv_G_koff`,
  both compat `#:else` branches, the `calc_EFG_halo` NS/LES branch, and the
  `LEGACY_EFG` scaffolding from Phase 1
- `calc_steps.f90` — `calc_dR_from_EFG`
- `preprocess.f90.fypp` — the E/F/G allocation and their `allocate_device_mem` params
- `calc_time_dev.f90.fypp` — E/F/G declarations, `calc_EFG`/`calc_EFG_halo`
  argument lists, deallocates, and the COMMZ `calc_dR_from_EFG` calls under
  `#:if VISC != 'Euler'` (all 6 RungeKutta variants)
- Verify `calc_slau_z_in_koff` is now unreferenced (its only caller was
  `calc_conv_G_koff`; COMMZ Euler uses fused `calc_conv_dG_koff`) and delete it too

### Retain (with an explanatory comment)

`calc_slau_{x,y,z}[_in]` and the `calc_slau_kernel[_internal]` modules — called
directly by `3D_solver_mixed`, not by any `3D_solver` path.

### Docs

- Update `CLAUDE.md`: drop the now-wrong E/F/G descriptions, document the three
  newly unsupported config combos, record the ~230 MiB device-memory reduction, and
  add a pointer to the benchmark report so it is discoverable from the Notice
  section that already indexes `3D_solver/nsys_ncu/`.
- Append a short **"Post-removal re-measurement"** section to
  `3D_solver/nsys_ncu/benchmark_dr_vs_efg.md` with the verification-step 4 numbers
  (no regression + actual memory drop), so the report reflects the shipped state
  rather than only the pre-removal experiment.

## Verification

1. Rebuild **all 8** 3D cases; also confirm `3D_solver_mixed`, `3D_solver_fltflt`,
   `3D_solver_curv`, and a 2D case still build (mixed is the one with a real
   coupling).
2. `pytest ouxsbli/tests/` — `test_etgv`, `test_st`, `test_evc`, `test_os`,
   `test_corn`. `test_etgv` is the key physics gate (it passed after the previous
   session's RK4 change and must still pass).
3. Confirm each of the three `#:stop` guards fires on a deliberately-configured
   throwaway case, and that all 8 shipped configs still preprocess cleanly.
4. Re-run the 128³ ETGV/NSTGV benchmark post-removal: confirm no perf regression
   and record the actual memory drop (`nvidia-smi` peak or `cudaMemGetInfo`), then
   append those numbers to the report per the Docs step above.
5. Clean up all throwaway case dirs and scratch workdirs; final `git status` review.
   The only new file that should remain is
   `3D_solver/nsys_ncu/benchmark_dr_vs_efg.md` (plus any retained raw profiling
   artifacts).

## Notes / risks

- Phase 1 edits `3D_solver/src/calc_flux_base.f90.fypp`, which is shared by every
  3D case — the `LEGACY_EFG` default of `False` keeps all existing builds on
  today's code path, but rebuild at least one case per branch shape to confirm.
- The fused y/z kernels use `atomicAdd`; if contention makes them *slower* at 128³
  than the E/F/G path, report that rather than forcing the removal — the whole
  point of Phase 1 is that the answer is measured, not assumed.
- 128 (not 129) is used as requested; for these periodic TGV cases that changes
  exact periodicity slightly but is irrelevant to a performance comparison, and
  both legs use the same grid.
