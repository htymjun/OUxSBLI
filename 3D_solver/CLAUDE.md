# CLAUDE.md — 3D_solver/

This file documents the 3D Cartesian solver. See the repository root `CLAUDE.md` for project overview, shared configuration reference (`config.fypp` variables, `mod_globals.f90`/`mod_constant.f90`), and tooling rules.

## fypp-preprocessed source files

| Template (in `src/` or `3D_solver/src/`) | Purpose |
|------------------------------------------|---------|
| `src/mod_constant.f90.fypp` | Kind-dispatch Fortran parameters from config values |
| `3D_solver/src/calc_flux_base.f90.fypp` | Top-level flux dispatcher (convective + viscous) |
| `3D_solver/src/calc_time_dev.f90.fypp` | RK time-stepping orchestration |
| `3D_solver/src/calc_keep_kernel.f90.fypp` | KEEP convective kernel |
| `3D_solver/src/calc_slau_kernel.f90.fypp` | SLAU convective kernel |
| `3D_solver/src/calc_roe_kernel.f90.fypp` | Roe convective kernel |
| `3D_solver/src/calc_hybrid_kernel.f90.fypp` | Hybrid KEEP↔SLAU kernel |
| `3D_solver/src/calc_keep_kernel_internal.f90.fypp` | Interior-only KEEP variant |
| `3D_solver/src/calc_slau_kernel_internal.f90.fypp` | Interior-only SLAU variant |
| `3D_solver/src/calc_roe_kernel_internal.f90.fypp` | Interior-only Roe variant |
| `3D_solver/src/calc_hybrid_kernel_internal.f90.fypp` | Interior-only Hybrid variant |
| `3D_solver/src/calc_visc2.f90.fypp` | 2nd-order viscous kernels |
| `3D_solver/src/calc_visc_high.f90.fypp` | 4th/6th-order viscous kernels (`VISC_ORDER`-controlled), boundary-aware (interior/fallback branch per point) |
| `3D_solver/src/calc_visc_high_internal.f90.fypp` | Interior-only 4th/6th-order viscous variant (used where a direction has no BC, i.e. periodic) |
| `3D_solver/src/calc_visc_cent.f90.fypp` | Shared `interp`/`diff`/`calc_tau_straight`/`calc_tau_cross` helpers, `#:include`d into both `calc_visc_high[_internal].f90.fypp` |
| `3D_solver/src/load_smem_visc_cent.f90.fypp` | Shared-memory loaders for `calc_visc_high[_internal]`'s `calc_Ev/Fv/Gv` |
| `3D_solver/src/preprocess.f90.fypp` | Device memory allocation helpers |

## Architecture

Source is split across directories; CMake/vpath merges them at build time:

- `src/` — shared templates and utilities: `mod_constant.f90.fypp`, `main.f90`, `cpu_gpu_mpi.f90`, `calc_muscl.f90`, `calc_physical_quantities.f90`, `print.f90`, `set_coordinate.f90`, `set_compressible_bl.f90`
- `3D_solver/src/` — solver kernels (mostly `.f90.fypp` templates): convective/viscous scheme modules, time-integration, preprocessing
- `3D_solver/<CASE>/` — per-case configuration: `mod_globals.f90` (grid/physical params, thread blocks), `set.f90` (grid/IC/BC), `config.fypp` (build-time scheme/method flags), `CMakeLists.txt`, `calc.sh`

### Data Flow (Cartesian)

```
main.f90
  └─ set_coordinate()        # grid metrics & Jacobians
  └─ set() [set.f90]         # grid, IC, BC (case-specific)
  └─ calc_time_dev()         # time-loop entry point  [fypp-generated]
        └─ calc_flux_base()  # dispatches to convective + viscous kernels  [fypp-generated]
        └─ calc_steps()      # RK stage update (Q += dt * RHS)
        └─ calc_para()       # MPI ghost-cell exchange (x) or z-halo (COMMZ)
        └─ print()           # VTK output
```

## Cases

| Case | Description |
|------|-------------|
| BL   | Quasi-2D laminar flat-plate boundary layer (extrudes 2D_solver/BL through a thin periodic z; validated vs. Blasius) |
| DHIT | Decaying homogeneous isotropic turbulence |
| ETGV | Supersonic Taylor-Green vortex (entropy-preserving) |
| EVC  | Quasi-2D Euler vortex convection (extrudes 2D_solver/EVC through a thin periodic z; grid-convergence study) |
| IVST | Inviscid vortex / smooth test case |
| KHI  | Kelvin-Helmholtz instability |
| NSTGV| NS Taylor-Green vortex |
| OS   | Quasi-2D oblique shock + wall reflection (extrudes 2D_solver/OS through a thin periodic z; validated vs. Rankine-Hugoniot) |
| SBLI | Shock-boundary layer interaction (RESCALE=True) |
| STZ  | z-direction halo exchange validation (COMMZ=True, BC_Z=True) |
| TBL  | Supersonic turbulent boundary layer, M=2.5 adiabatic wall (NS, KEEP, ORDER=6, RESCALE=True); validated against Guarini et al. (2000) — see the TBL note below |

BL, EVC, and OS are "quasi-2D" ports of their `2D_solver/` namesakes: the same
real physics extruded uniformly through a thin, periodic spanwise `z`
(`nz=9`, fixed — see below), which exercises 3D-only code paths (G-flux,
`calc_div_wz`, 3D grid metrics, the 3D VTK writer) that the 2D solver
structurally cannot reach, while still checking against the same closed-form
answer (Blasius, Rankine-Hugoniot, grid convergence) as the 2D case.

**`nz=9`, not the theoretical 6th-order-stencil minimum of 7.** Direct
testing on all three cases showed `nz=7` (1 real interior z-plane) trips a
GPU-kernel edge case that corrupts the domain within a single RK step, even
for a perfectly z-uniform field — reproduced independent of scheme and
viscosity. `nz=9` (3 interior planes) and `nz=13` (7 interior planes) both run
clean; `nz=9` was kept as the cheaper of the two, odd so that `nz/2` lands
exactly on the middle interior plane. Root cause not isolated further (would
require debugging the shared-memory tiling in `calc_slau_kernel.f90.fypp`/
`calc_keep_kernel.f90.fypp` for a case at the tiling's degenerate minimum
size); flagged here as a known limitation for any future `3D_solver/` case
that considers a similarly thin `nz`.

**BL's inlet is an active Dirichlet BC, unlike 2D_solver/BL's "leave i=1
frozen, no inlet BC" idiom** (`3D_solver/BL/set.f90`). Direct testing showed
that with `BC_X=True` and `i=1` left untouched, the 3D boundary-aware x
kernel corrupts columns `i=2..7` (exactly `ORDER` columns in) into NaN within
a single RK step — reproduced with both SLAU and KEEP, with and without
viscosity, and at every `nx` tried, so it is not scheme-, physics-, or
resolution-dependent. No existing 3D case had ever combined `BC_X=True` with
an untouched `i=1`: SBLI and TBL (the only other `BC_X=True` 3D cases) both
actively rewrite their leftmost columns every step via rescaling. Re-asserting
the same freestream state at `i=1` every step is numerically a no-op relative
to leaving it frozen, and avoids whatever this untested combination trips in
the boundary-aware kernel. OS was unaffected from the start: 2D_solver/OS's
own `set_bc` already actively rewrites `i=1` every step (a Dirichlet inflow),
so its 3D port inherited a working pattern without needing this fix.

### TBL: sizing the case against a reference Reynolds number

`blt` in `TBL/mod_globals.f90` is the **rescaling setpoint for inlet δ99**, and
every domain length is a multiple of it, so it is the single knob that sets the
Reynolds number. For the case's freestream (M=2.5, 100 kPa / 295 K total,
Re_unit = 9.87×10⁶ /m) the useful conversions are

```
Re_theta   = 751 * delta99[mm]          # from the profile shape, theta/delta99 = 0.0762
dRe_theta  = 13.9 per mm of x           # = Cf/2 * Re_unit at Cf = 0.00282
delta_nu   = 7.53 um                    # wall unit
```

so `blt = 1.9 mm` gives Re_θ ≈ 1430 at the inlet, growing to ≈ 1665 at the
outlet and crossing Guarini's 1577 near x = 0.63·Lx. Derive these from
`ouxsbli.analysis.tbl_stats.reference_station()` rather than by hand — it builds
the paper's own composite profile on the actual solver grid.

Do **not** size the domain from the composite δ used in the reference's
tables: δ99 is about 20% smaller (277 wall units versus 335 at Re_θ = 1577), so
using the composite value silently makes the box too small in δ99 units.

Changing `blt` invalidates a restart, because `Lx`, `Ly` and `Lz` all scale with
it and the grid changes. Calibrate before the long run, not after.

### TBL: running it

Two stages, switched by `stage.sh`, which patches `endT`/`np`/`step_offset` in
`mod_globals.f90` and `RESTART` in `config.fypp` and rebuilds:

```bash
cd 3D_solver/TBL
./calc.sh A            # spin-up,  60k steps = 12.1 flow-throughs -> recal/
./calc.sh B            # sampling, 120k steps = 24.2 flow-throughs, 101 snapshots
qsub -v STAGE=A job_miyabi.sh      # or job_tsubame.sh on TSUBAME
```

**Two MPI ranks, one GPU** — `main.f90` pairs them, so the even rank computes
and the odd rank runs the rescaling on the host and writes VTK. Unlike SWTBLI
this cannot be raised to 4; rescaling supports only the single `rerank` pair.

Stage A must run first: it carries the laminar-plus-noise initial condition
through transition, and only its `recal/Q00001.dat` + `recal/Qm.dat` let stage B
start from a stationary layer. Average **stage B only**.

`export FC=nvfortran` before any fresh `cmake -B build`. CMake's Fortran search
does not know `nvfortran`, so it silently picks gfortran and then fails at
`find_package(MPI)`; the case `CMakeLists.txt` cannot fix this because its
`project()` has already run by the time the shared file is included.

Statistics are computed offline from the snapshots by
`ouxsbli/analysis/tbl_stats.py`, which reads `.vtr` through
`ouxsbli/analysis/vtr_raw.py` — numpy only, no `vtk` package — so the ~13 GB of
stage B output can be reduced to a ~12 MB moment file on the compute node
before anything is copied back:

```bash
python -m ouxsbli.analysis.tbl_stats accumulate data --out tbl_acc.npz --halves
python -m ouxsbli.analysis.tbl_stats report tbl_acc.npz
```

## Convective Schemes (dispatched from `calc_flux_base.f90.fypp`)

| Scheme | Files | Best for |
|--------|-------|----------|
| KEEP (Kinetic Energy & Entropy Preserving) | `calc_keep_kernel.f90.fypp`, `calc_keep_kernel_internal.f90.fypp` | Smooth vortical flows (TGV, KHI) |
| SLAU (Simple Low-dissipation AUSM) | `calc_slau_kernel.f90.fypp`, `calc_slau_kernel_internal.f90.fypp` | Compressible turbulence with shocks (SBLI, TBL) |
| Roe | `calc_roe_kernel.f90.fypp`, `calc_roe_kernel_internal.f90.fypp` | Supersonic/hypersonic discontinuities |
| Hybrid (KEEP↔SLAU via Ducros sensor) | `calc_hybrid_kernel.f90.fypp`, `calc_hybrid_kernel_internal.f90.fypp` | Mixed smooth/shocked regions |

Viscous discretization: `calc_visc2.f90.fypp` (Gaitonde & Visbal 2nd-order, `VISC_ORDER=2`) or `calc_visc_high.f90.fypp` + `calc_visc_high_internal.f90.fypp` (4th/6th-order, `VISC_ORDER=4` or `6`).

## Cell-Center Velocity Gradients (`ux`, `vy`, `wz`)

For `VISC_ORDER > 2` (NS/LES), `calc_div.f90` precomputes the three diagonal
velocity gradients `ux=∂u/∂x`, `vy=∂v/∂y`, `wz=∂w/∂z` at cell centers once per
RK stage (`calc_div_${VISC_ORDER}$_in`, called from `calc_flux_base.f90.fypp`
right after `calc_quantities_T_3D`) and stores them to DRAM. `calc_Ev/Fv/Gv`
(in both `calc_visc_high.f90.fypp` and `calc_visc_high_internal.f90.fypp`)
then read the two gradients they don't own directly from these arrays instead
of each direction's kernel recomputing the other two independently in shared
memory — e.g. `calc_Ev` needs `vy,wz` (not `ux`, which it derives locally from
its own `u` shared-memory tile); it reads them straight from DRAM rather than
via `load_smem_visc_cent`, which now only loads the two cross-derivatives
unique to that direction (e.g. `uy,uz` for `calc_Ev`). Non-contiguous
(y/z-sweep) reads use scalar `calc_tau_straight_s`/`interp2_${N}$_s` helpers
in `calc_visc_cent.f90.fypp` (same pattern as the existing `mu`/`T` scalar
reads) rather than gathering into a local array.
Each of `ux,vy,wz` only needs a margin in its own direction (see
`calc_div_4_in`/`calc_div_6_in`'s independent per-component gating), which is
why one un-specialized `calc_div` kernel serves every consumer regardless of
that direction's `BC_X/Y/Z` setting. `calc_div_2` is not currently wired to
any consumer (`calc_visc2.f90.fypp` computes its trace terms directly from
`Q` via face differences, not a cell-centered gradient) — it's kept for a
possible future reuse by `calc_hybrid`'s `calc_Ducros` shock sensor at 2nd
order (see the `calc_Ducros` note in Notice below).
**Known limitation:** only wired into `calc_flux_base.f90.fypp`'s non-`COMMZ`
path; no current case combines `COMMZ=True` with `VISC in ('NS','LES')` at
`VISC_ORDER>2` (the only `COMMZ=True` case, STZ, is Euler-only), so
`calc_EFG_halo` still uses the old per-kernel recompute path.

## Adding a New Test Case

1. `cp -r 3D_solver/NSTGV 3D_solver/MYCASE`
2. Edit `mod_globals.f90` — grid size, physical parameters, thread block sizes
3. Edit `set.f90` — grid generation, initial conditions, boundary condition calls
4. Edit `config.fypp` — VISC, SCHEME, ORDER, BC_X/Y/Z, etc.
5. Edit `calc.sh` — MPI rank count and runtime args
6. `cmake -B build && cmake --build build -j`

## Notice
* From an occupancy perspective, the subroutines invoked within `calc_flux_base.f90` should not be executed on separate streams.
* `calc_flux_base.f90` and `calc_steps.f90` are the main bottleneck. You should optimize them.
* `calc_flux_base.f90` calls `calc_*_kernel.f90`, `calc_*_kernel_internal.f90`, `calc_visc*.f90`, and `calc_div.f90` (runs once per RK stage, right after `calc_quantities_T_3D`, for `VISC_ORDER>2`). They are the main bottleneck.
* Roe scheme is not used. KEEP, SLAU, Hybrid schemes should be optimized.
* STZ case exists specifically to validate the z-direction halo exchange (COMMZ=True, BC_Z=True); more validation and 6-point stencil support are still required.
* `RESTART=True` + `RESCALE=True` used to abort on startup: `pre_rescale`
  (`preprocess.f90.fypp`) requires `recal/Qm.dat`, but `write_Qm` — the only
  routine that wrote it — was never called from anywhere. `rescale_recv_send`
  now checkpoints the profile through `write_Qm_restart` every 1000 timesteps
  once rescaling has engaged, which is what makes a staged TBL run resumable.
  Note the resulting semantics: on restart `id_recal` is true, so `calc_mean` is
  skipped and `Qm` stays **frozen** at the restored profile for the whole run.
* Two known quirks in `calc_rescale.f90`, both benign, both left alone: the
  `jup` loop that sets `u99` has no `exit`, so `u99` is effectively
  `0.99*Um(ny)` (the freestream value you want anyway); and the `step` index
  passed from `calc_time_dev.f90.fypp` is `np*(t2-1)+t1` rather than
  `nt*(t2-1)+t1`, which makes the time column of `data/rescaling.d` a sawtooth
  but does not affect the running mean (`calc_mean` keeps its own counter).
* `calc_hybrid.f90`'s `calc_Ducros` shock sensor (used by SLAU/Roe/Hybrid schemes) recomputes its own full velocity-gradient tensor via always-2nd-order central differences, independent of `calc_div`'s higher-order `ux,vy,wz` — reusing `calc_div`'s output there would change the sensor's numerical values (mixed-order dilatation vs. vorticity) and needs a fallback for the edge band `calc_div` doesn't cover, so it hasn't been done; flagged here as a known optimization candidate for SLAU/Hybrid NS/LES cases (currently: SBLI).
* Nsight Compute profiling captures live in `3D_solver/nsys_ncu/*.ncu-rep` (open with `ncu --import <file> --print-summary per-kernel`, or the `.csv` triage exports); use these before speculating about kernel performance.
* `preprocess.f90.fypp`'s `pre_calc` (the initial t=0 snapshot send) used to hardcode `real(4)`/`MPI_REAL4` for its flux-flat buffers regardless of `OUTPUT_PRECISION`, unlike the rest of the print pipeline (`print.f90.fypp`'s `send_recv_for_print_even/odd3`), which already dispatched on the `io` kind parameter correctly. This was never exercised because no case had set `OUTPUT_PRECISION=8` before EVC/BL/OS's tests needed it (float32 can't resolve a 6th-order convergence trend or tight Cf tolerances); fixed to use `real(io)` and branch on `io` like the rest of the pipeline. No existing case's behavior changes (`io` still resolves to 4 for all of them).
