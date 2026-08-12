# Remove pipeline-async loads, reorder host Q layout, fypp-ify output precision

## Context

Three follow-up cleanups to `2D_solver` (and, for the third item, the shared
`src/print.f90` used by `3D_solver` and `3D_solver_curv` too), found while
auditing the codebase after the earlier SoA/fma port:

1. **`pipelineMemcpyAsync`** (from the compiler-provided `wmma` module) is used
   in `2D_solver/src/load_smem_visc2.f90` and `load_smem_visc4.f90` to
   asynchronously copy `Q_2`/`Q_3` into shared memory. `3D_solver/src/load_smem_visc2.f90`
   does the equivalent load with a **plain synchronous assignment** — no pipeline
   API at all. The user wants the 2D versions simplified to match (async-copy
   removed, direct assignment kept). This is the only place in `2D_solver` using
   the pipeline API; `3D_solver/src/calc_steps_smem.f90` and
   `3D_solver_curv/src/load_smem_visc2_curv.f90` also use it but are explicitly
   out of scope (not `2D_solver`).

2. **Host-side `Q` layout**: after the earlier SoA pass, `2D_solver`'s *device*
   arrays became `Q_1..Q_4(nx,ny)`, but host-only arrays (used only in `set.f90`,
   `set_init_common.f90`, `main.f90`, `preprocess.f90.fypp`'s `pre_calc`, the
   `_init` cyclic-BC variants in `set_bc_common.f90`, and `src/print.f90`) were
   deliberately left as `Q(nx,4,ny)` (component in the *middle* dimension),
   mirroring how 3D also left its host `Q` alone during its own SoA pass — except
   3D's host `Q` was *already* `Q(nx,ny,nz,5)` (component **last**), which is the
   convention 2D should now adopt too: `Q(nx,ny,4)`, indexed `Q(i,j,l)`.

3. **Output precision (`io` kind in `src/print.f90`)** is set via a bare C
   preprocessor line `#define io 4` (resolved by nvfortran's `-Mpreprocess`, not
   fypp) with zero conditional logic — it's hardcoded to single precision for
   *every* case, in both solvers. `test_evc.py`'s own docstring says as much:
   *"You should edit `/src/print.f90` file. io must be 8 because half precision
   (default) is not enough to check grid convergence."* This is exactly why
   `test_evc.py` currently fails even on unmodified code: single-precision VTK
   output has too little dynamic range to resolve a 4th/6th-order convergence
   trend (I confirmed in the prior session that the *pre-existing*, unmodified
   codebase reproduces the identical failing L2 numbers to 15 significant
   digits — same root cause). Converting `io` to an fypp-driven `config.fypp`
   value and having `test_evc.py` request `io=8` through `Case(...)` (instead of
   hand-editing `print.f90`) is expected to make the test measure a real
   convergence order instead of a precision floor.

## Change 1 — Remove `pipelineMemcpyAsync` from `2D_solver`

- `2D_solver/src/load_smem_visc2.f90`: in `load_smem_visc2_x`/`load_smem_visc2_y`,
  replace the `call pipelineMemcpyAsync(u(idx_l), Q_2(i,j))` /
  `call pipelineMemcpyAsync(v(idx_l), Q_3(i,j))` pair (inside the `if` branch)
  with plain `u(idx_l) = Q_2(i,j)`; `v(idx_l) = Q_3(i,j)`; drop the
  `call pipelineCommit()` / `call pipelineWaitPrior(0)` calls (keep
  `call syncthreads()`); drop `use wmma`.
- `2D_solver/src/load_smem_visc4.f90`: same pattern in `load_smem_visc4_x`/`load_smem_visc4_y`,
  but note the `pipelineCommit()` currently sits *between* the memcpy loop and a
  second loop that computes `uy`/`vy` (or `ux`/`vx`) directly from `Q_2`/`Q_3` in
  DRAM (independent of the shared-memory copy) — with synchronous assignment
  that second loop can simply follow the first with no intervening commit/wait,
  ending in one `call syncthreads()` before both are used. Drop `use wmma`.
- Match `3D_solver/src/load_smem_visc2.f90`'s style exactly (no pipeline calls,
  direct assignment, single trailing `syncthreads()`) as the reference pattern.

## Change 2 — Host `Q(nx,4,ny)` → `Q(nx,ny,4)` across `2D_solver`

Reorder every host (non-`device`) `Q`/`QJ` array from component-middle to
component-last, updating both the declaration and every `Q(i,l,j)`-style index
to `Q(i,j,l)`. Device SoA arrays (`Q_1..Q_4`, `QJ_1..QJ_4`) are untouched — this
is host-only.

Files (mechanical index-reorder in each):
- `2D_solver/{OS,ST,EVC,DSL,BL,BL0,SBLI}/set.f90` — each `set_init` subroutine's
  `Q(nx,4,ny)` declaration + all `Q(i,l,j) = ...` assignments.
- `2D_solver/src/set_init_common.f90` — `set_init_tbl`'s `Q(nx,4,ny)` + assignments
  (including the `Q(:,l,1)` boundary-fixup slices and the `Q(2,4,2)`-style scalar
  reads).
- `2D_solver/src/main.f90` — the `Q(nx,dimension+2,ny)` allocation, the
  post-loop `Q(i,m,j) = Jacobian(i,j) * Q(i,m,j)` rescale triple-loop, and the
  binary `read(10) Q` / `write(10) Q` restart calls (the on-disk layout changes
  with this — flag it, but `RESTART=False` on every current 2D case so no
  existing `recal/Q*.dat` is at risk).
- `2D_solver/src/preprocess.f90.fypp` — `pre_calc`'s `Q(nx,4,ny)` dummy arg and
  the divide-by-Jacobian loop (`Q(i,l,j) = Q(i,l,j) / Jacobian_cpu(i,j)`); the
  `QJ_1 = Q(:,1,:)`-style slice reads become `QJ_1 = Q(:,:,1)` etc.
- `2D_solver/src/calc_time_dev.f90.fypp` — both `RungeKutta` (RK3/RK4) `Q(nx,4,ny)`
  dummy-arg declarations (pass-through only, no index math on `Q` itself).
- `src/print.f90` (shared with 3D — only the 2D-facing procedures change; 3D's
  `Q(nx,ny,nz,5)` procedures are already component-last and untouched):
  `make_1d_for_print2`, `send_recv_for_print_even2` (+ its `Q(:,l,:) = QJ_l` slice
  writes), `send_recv_for_print_odd2`, `print0`.
- `2D_solver/src/set_bc_common.f90` — the CPU-only `_init` variants
  (`set_bc_cyclic2_init`, `set_bc_cyclic4_init`, `set_bc_cyclic6_init`), all
  `Q(a,:,b)`-style slice assignments become `Q(a,b,:)`.
- `2D_solver/BL0/extract.f90` — standalone post-processing tool reading the same
  binary `Q` file `main.f90` writes; update its `Q(nx,4,ny)` declaration and the
  `Q(nx1:nx2,l,j)`-style sums to keep it in sync with the new on-disk layout.

## Change 3 — fypp-ify output precision (`io` in `src/print.f90`)

- Convert `src/print.f90` → `src/print.f90.fypp`: replace the leading
  `#define io 4` with
  ```
  #:include 'config.fypp'
  #:if not defined('OUTPUT_PRECISION')
    #:set OUTPUT_PRECISION = 4
  #:endif
  ...
  integer, parameter :: io = ${OUTPUT_PRECISION}$
  ```
  placed as a module-level parameter (the rest of the file — the runtime
  `if (io == 4) ... else ...` branches selecting `MPI_REAL4`/`Float32` vs
  `MPI_REAL8`/`Float64` — needs no change; those aren't C macros, they're
  ordinary Fortran conditionals on a parameter and stay as-is per the user's
  literal ask ("use fypp instead of C-style macro" targets the *definition*,
  not every downstream branch).
- Add `#:set OUTPUT_PRECISION = 4` to every case's `config.fypp` (7× `2D_solver`,
  8× `3D_solver`, `3D_solver_curv/{NACA,CORN}`) so the default (single
  precision, current behavior) is explicit and discoverable alongside `ORDER`/`TVD`/etc.,
  matching the existing `#:if not defined(...)` fallback pattern used for
  `VISC_ORDER` — cases that don't set it still default to 4 via the fallback
  above, so this is a documentation/discoverability addition, not a required edit.
- Wire `print` into the fypp-generation step (alongside `mod_constant`/`calc_muscl`
  in `_PARENT_FYPP`) instead of `_BASE`'s static-file list, in all three
  CMakeLists that reference `print.f90`: `2D_solver/CMakeLists.txt`,
  `3D_solver/CMakeLists.txt`, `3D_solver_curv/CMakeLists.txt`. It stays a real
  `module print` (unlike `calc_scheme_math.f90.fypp`, which is include-only) so
  it's still added to `_GEN`/`add_executable` — no need for the include-only
  custom-target trick used for `calc_scheme_math`.
- Add `"output_precision": "OUTPUT_PRECISION"` to `ouxsbli/case.py`'s `_ALIAS`
  table so `Case(..., output_precision=8)` works.
- Update `ouxsbli/tests/test_evc.py`: pass `output_precision=8` in the `Case(...)`
  call in `_run_evc`, and remove the now-stale docstring note about manually
  editing `print.f90`.

## Verification

- After Change 1: rebuild all 7 `2D_solver` cases; rerun `test_os.py`/`test_st.py`
  (`-m integration`) — expect bit-identical pass (synchronization-strategy-only
  change, no math difference).
- After Change 2: rebuild all 7 cases again; rerun `test_os.py`/`test_st.py`;
  smoke-run `BL` and `SBLI` a few seconds (as done previously) and check VTK
  output for NaNs/sane ranges, since these two exercise `set_init_common.f90`'s
  `set_init_tbl` most heavily.
- After Change 3: rebuild all 7 `2D_solver` + all 8 `3D_solver` + `3D_solver_curv/{NACA,CORN}`
  cases (print.f90.fypp now generated everywhere) to catch any case where the
  fypp expansion breaks; run the **full** `pytest ouxsbli/tests/` (`test_os`,
  `test_st`, `test_evc`, `test_etgv`, `test_corn`) — `test_evc.py` is the key
  signal: with real double-precision output it should now measure genuine
  convergence order. If it still fails to meet the ≥3.6/≥5.7 thresholds even at
  `io=8`, that indicates a separate, real numerical issue worth reporting
  rather than a blocking failure of this task (the task here is fixing the
  precision *mechanism*, not guaranteeing the scheme's order).
- Spot-check the fypp-generated `print.f90` in at least one 2D and one 3D case's
  `build/` directory for correct `io` resolution before trusting the rebuild,
  per the project's existing fypp-template convention.

## Additional optimization/refactoring observations (report only, not implemented here)

While making these changes I'll keep an eye out for anything else worth flagging
(e.g. the `BL/set.f90` local `Qp(4,ny)`/`Qi(5,nyi)` helper arrays use a third,
different component-first layout, unrelated to `Q(nx,4,ny)` but inconsistent
with the rest of the codebase; whether `calc_div`'s pattern could extend to
`VISC_ORDER=6` if ever needed; anything else surfaced by re-reading these files
closely) and summarize them at the end for you to consider — no code changes
beyond the three items above without asking first.
