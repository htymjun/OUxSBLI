# Finish the sliding-window (slice-window) refactor of `calc_steps.f90`

## Context

`3D_solver/src/calc_steps.f90` currently has an **uncommitted** working-tree change that already implements the "slice window" idea: `calc_step1`, `calc_step`, `calc_step2_3`, and `calc_step4` were converted from fully 3D-thread-indexed kernels (`i,j,k` all derived from `blockIdx`/`threadIdx`) to 2D-thread-indexed kernels where each thread owns one `(i,j)` column and internally loops `do k = 1, nz-2`, keeping only two z-planes of the z-flux `G` in registers (`G_before(5)`, `G_after(5)`) instead of re-reading `G(:,i,k,j)` and `G(:,i,k+1,j)` from global memory every iteration. This halves redundant global loads of `G` and follows naturally from the just-committed `7726a30 "transpose G"`, which reordered `G` to `G(5,nx-2,nz-1,ny-2)` so that consecutive-`k` reads are memory-adjacent.

The kernel bodies are correct and complete (confirmed by reading the file and tracing that `E`/`F`/`G` are fully populated in global memory — even under the `COMMZ` overlapped-halo path — before any `calc_step*` launch, so reading arbitrary `k` per thread is safe). **What's missing is that the launch configuration was never updated to match.** `calc_step1/step/step2_3/step4` are still launched as `<<<blocks,threads>>>` in `3D_solver/src/calc_time_dev.f90.fypp`, where `blocks` is a 3D grid computed by `set_block_3D` in `src/set_coordinate.f90:41` as:

```fortran
blocks = dim3((nx-2+threads%x-1)/threads%x, (ny-2+threads%y-1)/threads%y, (nz-2+threads%z-1)/threads%z)
```

Since `threads%z` is `1` in every case's `mod_globals.f90`, this still launches `nz-2` blocks in the z-direction (e.g. 511 for a 513³ grid). The new kernel bodies never reference `blockIdx%z` — every one of those `nz-2` blocks would redundantly perform the *entire* z-sweep for the same `(i,j)` columns. That's not just wasted work, it defeats the entire point of the optimization (turns an O(N) launch into an O(N²)-ish one). This must be fixed for the refactor to actually pay off.

Critically, `blocks` **cannot simply be redefined as 2D**, because it's also used, unchanged, by two genuinely 3D-indexed elementwise kernels dispatched from `3D_solver/src/calc_flux_base.f90.fypp`: `calc_Ducros` (`3D_solver/src/calc_hybrid.f90:11`, uses `blockIdx%z` directly as `k`) and `calc_mut` (`3D_solver/src/calc_les.f90:96`, same pattern). Those must keep the full 3D grid. So a new, dedicated grid variable is needed for the step kernels only.

## Plan

**1. Add a new 2D grid variable `blocksStep`, used only by the RK step kernels.**

- `src/set_coordinate.f90`, `set_block_3D` (lines 30–42): add `blocksStep` as a new `intent(out) :: type(dim3)` argument, computed as:
  ```fortran
  blocksStep = dim3((nx-2+threads%x-1)/threads%x, (ny-2+threads%y-1)/threads%y, 1)
  ```
  (mirrors the existing 2D pattern already used for `blocksE/blocksF` in `set_block_2D`, just with `threads` instead of `threadsE/threadsF`). Leave `blocks` itself untouched — it keeps its current 3D meaning for `calc_Ducros`/`calc_mut`.

- `3D_solver/<CASE>/mod_globals.f90` for all 8 cases (`DHIT, ETGV, IVST, KHI, NSTGV, SBLI, STZ, TBL`): add `blocksStep` to the existing `type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks` declaration line.

- `3D_solver/src/main.f90`: add `blocksStep` to the `use mod_globals, only : ...` list (lines 4-6) and to the `call set_block_3D(...)` argument list (lines 30-31).

- `3D_solver/src/calc_time_dev.f90.fypp`: add `blocksStep` to the `use mod_globals, only : ...` list (line 9), then change every `calc_step1<<<blocks,threads>>>`, `calc_step<<<blocks,threads>>>`, `calc_step2_3<<<blocks,threads>>>`, `calc_step4<<<blocks,threads>>>` launch to `<<<blocksStep,threads>>>` — this covers all 6 `RungeKutta` variants (RK3/RK4 × plain/RESCALE/COMMZ), 24 call sites total (lines 88,94,99,171,181,191,269,274,279,284,354,364,374,384,473,484,496,582,593,604,615). `calc_flux_base.f90.fypp`'s `calc_Ducros`/`calc_mut` launches keep `<<<blocks,threads>>>` unchanged.

**2. Small cleanup in `3D_solver/src/calc_steps.f90` itself.**

The in-progress diff added `use mod_globals, only : dt, threads` to the module header, but `threads` is never referenced anywhere in this file (kernel *bodies* don't touch `dim3` launch-config variables — those only matter at the call site). Remove `threads` from that import; keep `dt` as-is (pre-existing, unrelated to this task).

**3. Verification.**

- Build a plain case, a RESCALE case, and a COMMZ case to exercise all 6 fypp branches in `calc_time_dev.f90.fypp` and confirm `blocksStep` wiring compiles cleanly:
  ```bash
  cd 3D_solver/NSTGV && cmake -B build && cmake --build build -j   # RK3, plain
  cd 3D_solver/SBLI  && cmake -B build && cmake --build build -j   # RESCALE
  cd 3D_solver/STZ   && cmake -B build && cmake --build build -j   # COMMZ
  ```
  Spot-check the generated `build/main.f90` / `build/calc_time_dev.f90` to confirm `blocksStep` is threaded through correctly (per CLAUDE.md guidance to always check fypp-generated output).
- Run one case (`mpirun -n 2 ./a.out` in `3D_solver/NSTGV/build`) and compare against a pre-change baseline (stash the changes, rerun) — since this is purely a memory-access/launch-config optimization, results should be numerically identical to before.
- Run `pytest ouxsbli/tests/` if any of the tested cases (`ETGV` for `test_etgv.py`) go through this code path, to catch any regression automatically.
