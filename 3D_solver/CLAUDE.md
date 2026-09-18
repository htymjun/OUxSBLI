# CLAUDE.md — 3D_solver/

This file documents the 3D Cartesian solver. See the repository root `CLAUDE.md` for project overview, shared configuration reference (`config.fypp` variables, `mod_globals.f90`/`mod_constant.f90`), and tooling rules.

## fypp-preprocessed source files

| Template (in `src/` or `3D_solver/src/`) | Purpose |
|------------------------------------------|---------|
| `src/mod_constant.f90.fypp` | Kind-dispatch Fortran parameters from config values |
| `3D_solver/src/calc_flux_base.f90.fypp` | Top-level flux dispatcher (convective + viscous) |
| `3D_solver/src/calc_time_dev.f90.fypp` | RK time-stepping orchestration |
| `3D_solver/src/calc_para.f90.fypp` | MPI halo exchange: x-direction pack/unpack and the COMMZ z-halo `start_exchange_z`/`finish_exchange_z` (`MPI='CPU'` host staging or `MPI='GPU'` CUDA-aware) |
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
| SWLBLI | Oblique-shock / laminar flat-plate boundary-layer interaction, M=2.15 (NS, Hybrid, ORDER=6, BC_X=BC_Y=True, periodic z); the only NS case that can run with `COMMZ=True` (`MPI='GPU'`) — see the COMMZ section below. `BC_FORCING=True`: its `set_bc` takes `(x, z, phi_l_gpu, phi_m_gpu, t_now)` for a blowing/suction strip (currently commented out inside `set_bc`) |
| SWTBLI | Oblique-shock / turbulent boundary-layer interaction (rescaling inflow, RESCALE branch of `calc_time_dev`) |
| TBL  | Supersonic turbulent boundary layer, M=2.5 adiabatic wall (NS, KEEP, ORDER=6, RESCALE=True); validated against Guarini et al. (2000) — see the TBL notes below |

BL, EVC, and OS are "quasi-2D" ports of their `2D_solver/` namesakes: the same real physics extruded uniformly through a thin, periodic spanwise `z` (`nz=9`, fixed), which exercises 3D-only code paths (G-flux, `calc_div_wz`, 3D grid metrics, the 3D VTK writer) that the 2D solver structurally cannot reach, while still checking against the same closed-form answer (Blasius, Rankine-Hugoniot, grid convergence). Two constraints found by direct testing, both reproduced independent of scheme, viscosity and resolution:

* **`nz=9`, not the theoretical 6th-order-stencil minimum of 7.** `nz=7` (1 real interior z-plane) trips a GPU-kernel edge case that corrupts the domain within a single RK step, even for a perfectly z-uniform field. `nz=9` (3 interior planes) and `nz=13` both run clean; `nz=9` was kept as the cheaper, and is odd so `nz/2` lands exactly on the middle interior plane. Root cause not isolated (would require debugging the shared-memory tiling in `calc_slau_kernel.f90.fypp`/`calc_keep_kernel.f90.fypp` at the tiling's degenerate minimum size) — flagged as a known limitation for any future case that considers a similarly thin `nz`.
* **BL's inlet is an active Dirichlet BC**, unlike 2D_solver/BL's "leave `i=1` frozen, no inlet BC" idiom (`3D_solver/BL/set.f90`). With `BC_X=True` and `i=1` left untouched, the 3D boundary-aware x kernel corrupts columns `i=2..7` (exactly `ORDER` columns in) into NaN within a single RK step. No other 3D case had ever combined `BC_X=True` with an untouched `i=1`: SBLI and TBL actively rewrite their leftmost columns every step via rescaling, and OS inherited a working Dirichlet inflow from 2D_solver/OS. Re-asserting the same freestream state at `i=1` every step is numerically a no-op relative to leaving it frozen.

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

so `blt = 1.9 mm` gives Re_θ ≈ 1430 at the inlet, growing to ≈ 1665 at the outlet
and crossing Guarini's 1577 near x = 0.63·Lx. Derive these from
`ouxsbli.analysis.tbl_stats.reference_station()` rather than by hand — it builds the
paper's own composite profile on the actual solver grid. Do **not** size the domain
from the composite δ used in the reference's tables: δ99 is about 20% smaller (277
wall units versus 335 at Re_θ = 1577), so the composite value silently makes the box
too small in δ99 units. Changing `blt` invalidates a restart, because `Lx`, `Ly` and
`Lz` all scale with it and the grid changes — calibrate before the long run.

### TBL: running it

Two stages, switched by `stage.sh`, which patches `endT`/`np`/`step_offset` in
`mod_globals.f90` and `RESTART` in `config.fypp` and rebuilds:

```bash
cd 3D_solver/TBL
./calc.sh A            # spin-up,  60k steps = 12.1 flow-throughs -> recal/
./calc.sh B            # sampling, 120k steps = 24.2 flow-throughs, 101 snapshots
qsub -v STAGE=A job_miyabi.sh      # or job_tsubame.sh on TSUBAME
```

* **Two MPI ranks, one GPU** — `main.f90` pairs them, so the even rank computes and
  the odd rank runs the rescaling on the host and writes VTK. Unlike SWTBLI this
  cannot be raised to 4; rescaling supports only the single `rerank` pair.
* Stage A must run first: it carries the laminar-plus-noise initial condition through
  transition, and only its `recal/Q00001.dat` + `recal/Qm.dat` let stage B start from
  a stationary layer. Average **stage B only**.
* `export FC=nvfortran` before any fresh `cmake -B build`. CMake's Fortran search does
  not know `nvfortran`, so it silently picks gfortran and then fails at
  `find_package(MPI)`; the case `CMakeLists.txt` cannot fix this because its
  `project()` has already run by the time the shared file is included.
* Statistics are computed offline from the snapshots by `ouxsbli/analysis/tbl_stats.py`,
  which reads `.vtr` through `ouxsbli/analysis/vtr_raw.py` — numpy only, no `vtk`
  package — so the ~13 GB of stage B output can be reduced to a ~12 MB moment file on
  the compute node before anything is copied back:

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

## KEEP_TVD: giving the KEEP flux a limiter

`SCHEME='KEEP'` is a central, kinetic-energy- and entropy-preserving flux with no
numerical dissipation, and `calc_keep_kernel[_internal].f90.fypp` never reference
`id_tvd` — so `TVD` in `config.fypp` reaches SLAU/Roe/Hybrid only and is **silently
inert for KEEP**. On a flow with a discontinuity nothing damps the 2-cell mode and the
solution rings; that is what "KEEP oscillates on STZ's Sod tube but SLAU does not"
means, and it is the scheme behaving as designed, not a bug.

`KEEP_TVD = True` (default `False`) makes the face flux

```
F = F_KEEP  -  1/2 |lambda| (U_R - U_L),     lambda = max(|u_n| + c)_{L,R}
```

with `U_L`/`U_R` reconstructed by the **same** `delta4`/`delta6` MUSCL chain SLAU and
Hybrid already use, so `TVD` picks the limiter for KEEP exactly as it does for them.
The algebra lives in `3D_solver/src/calc_keep_tvd.f90.fypp`, which is `#:include`d into
both KEEP kernel modules (raw text, like `calc_visc_cent.f90.fypp` — so it is **not** in
CMake's `_SHARED_FYPP` and CMake does not track it as a dependency; touch a kernel
template or wipe `build/` after editing it).

Measured on the host with the real limiter chain (`ORDER=6`, `TVD='tvd'`):

| field at the face | `\|D(1)\|` |
|---|---|
| smooth `rho = 1 + 0.1 sin(2*pi*z)`, dz = 5.0e-2 | 6.0e-4 |
| same, dz = 1.25e-2 | 3.6e-9 |
| same, dz = 3.1e-3 | 3.3e-12 |
| Sod jump (1 -> 0.125, p 1 -> 0.1) | 5.2e-1 |

i.e. once the wave is resolved the added dissipation converges at ~O(dz^5) — below the
6th-order convective truncation error — while a genuine discontinuity gets the full
O(1) upwind term. KEEP is still no longer *exactly* entropy preserving, so do not switch
this on for TGV/DHIT-class runs that depend on that property; `ouxsbli/tests/test_etgv.py`
is the check that would catch it.

Constraints, enforced by `$:error` in `calc_flux_base.f90.fypp`:

* `KEEP_TVD` requires `SCHEME='KEEP'`.
* `TVD='hybrid'` is rejected: the threshold limiter needs the Ducros sensor, which is not
  passed to the KEEP kernels. Use `TVD='tvd'`, or `SCHEME='Hybrid'`, which already
  switches KEEP<->SLAU on that sensor and leaves KEEP untouched where it is quiet.
* `KEEP_TVD=True` disables the fused convective+viscous kernel (`fused_conv_visc`), so
  `calc_keep_visc_kernel.f90.fypp` did not need a second copy of the dissipation. Costs
  one extra launch per direction on the no-BC/no-COMMZ cases that would have fused.

With `KEEP_TVD` unset or `False` every existing case generates **byte-identical**
Fortran — verified across all 12 `3D_solver` cases.

## Cell-Center Velocity Gradients (`ux`, `vy`, `wz`)

For `VISC_ORDER > 2` (NS/LES), `calc_div.f90` precomputes the three diagonal velocity
gradients `ux=∂u/∂x`, `vy=∂v/∂y`, `wz=∂w/∂z` at cell centers once per RK stage
(`calc_div_${VISC_ORDER}$_in`, called from `calc_flux_base.f90.fypp` right after
`calc_quantities_T_3D`) and stores them to DRAM. `calc_Ev/Fv/Gv` (in both
`calc_visc_high.f90.fypp` and `calc_visc_high_internal.f90.fypp`) then read the two
gradients they don't own straight from those arrays instead of recomputing them in
shared memory — e.g. `calc_Ev` needs `vy,wz` (not `ux`, which it derives locally from
its own `u` tile), so `load_smem_visc_cent` now only loads the two cross-derivatives
unique to that direction (`uy,uz` for `calc_Ev`). Non-contiguous (y/z-sweep) reads use
the scalar `calc_tau_straight_s`/`interp2_${N}$_s` helpers in `calc_visc_cent.f90.fypp`
(same pattern as the existing `mu`/`T` scalar reads) rather than gathering into a local
array. Each of `ux,vy,wz` only needs a margin in its own direction (see
`calc_div_4_in`/`calc_div_6_in`'s independent per-component gating), which is why one
un-specialized `calc_div` kernel serves every consumer regardless of `BC_X/Y/Z`.
`calc_div_2` is not wired to any consumer (`calc_visc2.f90.fypp` computes its trace
terms directly from `Q` via face differences) — it's kept for possible reuse by
`calc_hybrid`'s `calc_Ducros` sensor at 2nd order (see Notice). Every
`calc_div_*_{4,6}_in` kernel takes a z range `(k_lo, k_hi)`: the non-COMMZ path passes
the whole slab (`1, nz` for `ux`/`vy`; `io_v+2, nz-io_v-1` for `wz`), the COMMZ path
splits it around the halo exchange as described next.

## COMMZ: z-halo exchange and what may run before it lands

`COMMZ=True` decomposes z across the even (compute) ranks, one slab per rank
with `OV = ORDER//2` ghost planes on each z end (`OVV = VISC_ORDER//2` for the
viscous stencil; equal to `OV` in every current case). Per RK stage
`calc_time_dev` does: `start_exchange_z` (pack, non-blocking MPI) →
`calc_EFG` (stencil-safe interior only) → `finish_exchange_z` (wait, unpack)
→ `calc_EFG_halo` (everything whose stencil reaches a ghost plane) → RK
update → `set_bc`. The ranges are derived from the stencils, not tuned, and
are listed at the top of `calc_flux_base.f90.fypp`:

| quantity | stencil (planes) | before the exchange | after the exchange (`calc_EFG_halo`) |
|---|---|---|---|
| primitives, `mu`, `T`, `ux`, `vy` | k | `[OV+1, nz-OV]` | `[1, OV]`, `[nz-OV+1, nz]` |
| Ducros sensor | k±1 | whole slab | whole slab again |
| `wz` | k±OVV | `[OV+1+OVV, nz-OV-OVV]` | `[OV+1, OV+OVV]`, `[nz-OV-OVV+1, nz-OV]` |
| E/F convective (plane k) | k±1 via the sensor | `[OV+2, nz-OV-1]` | `[2, OV+1]`, `[nz-OV, nz-1]` |
| G convective (face k) | k-io..k+io+1 and sensor | `[GLO, nz-GHI]` (= `[2OV, nz-2OV]` for ORDER 4/6) | `[1, GLO-1]`, `[nz-GHI+1, nz-1]` |
| E/F viscous (plane k) | k±OVV (`uz`/`vz`, `wz`) | `[OV+1+OVV, nz-OV-OVV]` | `[2, OV+OVV]`, `[nz-OV-OVV+1, nz-1]` |
| G viscous (face k) | k-OVV+1..k+OVV | `[OV+OVV, nz-OV-OVV]` | `[1, OV+OVV-1]`, `[nz-OV-OVV+1, nz-1]` |

Rules that make this correct, and that any edit must preserve:

* **Convective before viscous, always.** The convective kernels *assign* the flux
  array and the viscous kernels *increment* it, so within a stage every face must
  receive its convective value before its viscous increment. The pre-fix code
  violated this: `calc_EFG` subtracted the z viscous flux on faces `[OV+1, nz-OV-1]`
  and `calc_EFG_halo` then re-assigned the convective flux on `[OV, 2OV-1]`, silently
  deleting τ_zz on faces `OV+1..2OV-1` (3 of the 4 periodic faces for the nz=10 SWLBLI
  slab). Since τ_zz = −(2/3)μ(u_x+v_y) ≠ 0 for z-uniform flow, the missing faces acted
  as a deterministic ρw source: after 1000 steps the pre-fix `COMMZ=True` binary had
  max|w| ≈ 0.24 m/s and Δp ≈ 150 Pa between z planes of a nominally 2D SBLI, while
  `COMMZ=False` stays 2D to the last bit. That was the "COMMZ=True transitions,
  COMMZ=False does not" report.
* **The `_z_in_koff` convective kernels clamp to `[io+1, nz-io-1]` like their non-koff
  twins.** `calc_conv_G_koff` is called with `k_lo=1` and `k_hi=nz-1`, i.e. over faces
  whose `k-io .. k+io+1` stencil is not inside the slab. The koff kernels used to test
  only `k_lo <= k <= k_hi`, so those faces were written from shared-memory slots the
  loader never filled (and, for SLAU/Hybrid, from out-of-bounds `T(i,j,k-io:k+io+1)`
  reads). Only ghost planes consume those faces and `set_bc_cyclic_z` overwrites them
  every stage, so it never corrupted the solution -- but it made the ghost planes
  run-to-run nondeterministic, which is exactly what the `COMMZ=True` vs `COMMZ=False`
  bit-comparison relies on being clean.
* **The high-z range of every pair starts at `max(..., low_end+1)`** so thin slabs
  (nz-2·OV small) never double-count; empty ranges are skipped (`if (hi >= lo)`, and
  the `calc_div` kernels return early).
* **The E/F x/y kernels have `_koff` twins** (`calc_<scheme>_x[_in]_koff`, generated
  only when `COMMZ` is set, so non-COMMZ generated code is byte-identical) that sweep
  `[k_lo, k_hi]`; `calc_conv_EF_koff` / `calc_conv_G_koff` size the z grid from the
  range (`zblocks`).
* **With one compute rank (2 MPI ranks) the exchange is a self-copy** equal to
  `set_bc_cyclic_z`, so `COMMZ=True` must reproduce `COMMZ=False` bit for bit; that
  comparison (VTK and `recal/Q*.dat`) is the regression test for this path.
* `set_bc_cyclic_z` stays in SWLBLI's `set_bc` under COMMZ on purpose: it is a finite
  placeholder for the ghost planes between stages (nothing reads them before
  `finish_exchange_z`), and `calc_time_dev` repeats the exchange right before
  `send_recv_for_print_even` so VTK/restart ghost planes hold the neighbours' data.
* `MPI='GPU'` hands device buffers to MPI; the `cudaDeviceSynchronize` in
  `start_exchange_z` is **required** (flatten kernels vs. MPI's own stream/IPC/RDMA
  reads, and the previous stage's `reconstruct_z_*` vs. the new `MPI_Irecv`). Removing
  it makes results run-to-run nondeterministic.
* **An MPI stack that cannot really move device memory corrupts the halo without
  reporting anything**, so `check_exchange_z` (`calc_para.f90.fypp`) runs one pattern
  exchange before the time loop and aborts with a diagnosis instead of letting the run
  drift into NaN. Not hypothetical: on the WSL2 dev box (one GPU, HPC-X 2.19,
  `opal_built_with_cuda_support=true`, UCX cuda transports present) a 4-rank SWLBLI run
  used to turn the whole field into NaN within one output interval, and a 90-line
  standalone probe doing nothing but this exchange still returns **every received
  element zero** at 4 ranks (the receive buffer goes -1 → 0, so the fill kernel ran and
  MPI wrote the zeros). Two ranks always work because `rank_lo == rank_hi == self`.
  Forcing the transport fixes the standalone case:

  | 4-rank device-pointer exchange (standalone probe) | result |
  |---|---|
  | default UCX transport selection | **corrupted (all zeros)** |
  | `UCX_MEMTYPE_CACHE=n` | **corrupted** |
  | `UCX_TLS=self,sm,cuda_copy` | correct |
  | `--mca pml ob1 --mca btl self,vader,smcuda --mca btl_smcuda_use_cuda_ipc 0` | correct |

  The solver's own 4-rank `MPI='GPU'` runs became clean once GPU selection was made
  node-local (below) — 3/3 runs finite, z-uniform to the last bit and bit-identical to
  the 2-rank `COMMZ=False` reference — even though the standalone probe, which calls
  `cudaSetDevice` on every rank including the I/O ranks, still fails. So the
  device-pointer path on this host is at best fragile and sensitive to how the CUDA
  contexts are set up. A cluster can hit the same class of failure for its own reasons
  (UCX_TLS misconfigured, `nvidia_peermem` not loaded, ranks in containers without IPC
  access), and the failure is silent wrong data rather than an error, which is why the
  check lives in the code. `MPI='CPU'` (the fypp default) sidesteps the question and is
  verified at 2 and 4 ranks.
* Related, and independent of any MPI stack: GPU selection is **node-local**
  (`main.f90` splits `MPI_COMM_TYPE_SHARED` and uses the rank within the node, modulo
  the device count). The old global `mygpu = myrank / 2` is correct only on one node —
  across several nodes the higher ranks index devices their node does not have, and
  `cudaSetDevice` returns success for an out-of-range device (measured: it returned 0
  for device 1 on a single-GPU host, and the rank then ran with a context that made its
  pointers unusable to MPI), so `check_gpu` range-checks it too.
* Not implemented / not verified: `COMMZ` with `VISC_ORDER=2` (fypp error), `LES` +
  `COMMZ` (`calc_mut` still runs before the exchange), `RK=4` + `COMMZ` (compiles at
  both `MPI='CPU'` and `MPI='GPU'`, never run), `BC_FORCING` + `RESCALE`.
* SWLBLI's `set_grid` follows STZ: `dz = Lz / (N_compute·(nz-6))`, plane `k=4` of the
  first slab at `z=-Lz/2`, so `Lz` is the global periodic span and `nz` counts the 3+3
  ghost planes.
* **`BC_Z=True` is inoperative under `COMMZ`, and STZ relies on it.** `start_exchange_z`
  wraps unconditionally (`rank_lo = mod(myrank-2+nranks, nranks)`), so the first and
  last compute ranks exchange with each other; `calc_conv_G_koff` likewise always uses
  the `_in` z kernel. STZ's `set_bc` does write zero-gradient z ghost planes, but only
  *after* `calc_step`, and the next stage's `finish_exchange_z` overwrites them before
  any ghost-dependent flux is evaluated — so those values are never read and the run is
  effectively z-periodic. For STZ's Sod tube (`set_init`: left state on `myrank < 2`,
  right state otherwise) that adds a **second** discontinuity at the periodic seam, and
  the shock/contact/expansion it launches travel inward from both domain ends. Any
  `COMMZ` case that needs a real physical z boundary must gate the wrap on `BC_Z` and
  move the BC to just after `finish_exchange_z`.
* **STZ cannot be compared against a single-rank reference as written.** `set_init`
  selects its state from `myrank`, not from `z`, so at 2 ranks (one compute rank) the
  whole domain is the left state and there is no shock tube at all. To use STZ as a
  halo regression test, give it a smooth z-dependent initial condition (e.g.
  `rho = 1 + 0.1*sin(2*pi*z/Lz)`) written in terms of `z(k)`, which is already
  globally offset by `iz_offset`.
* **`TVD` is a no-op for `SCHEME='KEEP'`.** The KEEP kernel templates never reference
  `id_tvd` and `KEEP2/4/6` carry no limiter, so `TVD='tvd'`/`'hybrid'` only affects
  SLAU/Roe/Hybrid. KEEP on a shock tube therefore runs fully unlimited.

### SWLBLI: the blowing/suction strip (`BC_FORCING`)

`BC_FORCING=True` only changes the `set_bc` *signature* (`(x, z, phi_l_gpu, phi_m_gpu, t_now)`); the strip itself is still commented out inside `SWLBLI/set.f90`. The phases are drawn once in `calc_time_dev` and `MPI_BCAST` so every rank sees the same `phi_l`/`phi_m`, and `t_now` is the step index, so all RK stages of a step share one forcing time (RK3 and RK4 alike). Four things to check before switching the strip on — all found by reading the code, none tested:

* `h_t` sums `sin(beta_force*t + 2π φ_m)` over `m = m_min..m_max` **without a `dble(m)` factor**, unlike `g_z` which does carry `dble(l)`; as written every temporal mode has the same frequency and the sum collapses to a single sine.
* `beta_force = 75000` is commented `(Hz)` but enters as `sin(β t)`, i.e. 75000 rad/s ≈ 11.9 kHz — multiply by 2π if Hz was meant.
* Spanwise modes `l = 12..24` need a span that resolves them: at the current `Lz = 0.5·blt` with `nz-6` planes per rank their wavelengths are far below `dz`.
* `Z_l`, `T_m`, `xs_gpu`, `zs_gpu` are rebuilt and copied host→device inside `set_bc`, i.e. once per RK stage, although only `t_now` changes between calls.

## Adding a New Test Case

1. `cp -r 3D_solver/NSTGV 3D_solver/MYCASE`
2. Edit `mod_globals.f90` — grid size, physical parameters, thread block sizes
3. Edit `set.f90` — grid generation, initial conditions, boundary condition calls
4. Edit `config.fypp` — VISC, SCHEME, ORDER, BC_X/Y/Z, etc.
5. Edit `calc.sh` — MPI rank count and runtime args
6. `cmake -B build && cmake --build build -j`

## Notice
* From an occupancy perspective, the subroutines invoked within `calc_flux_base.f90` should not be executed on separate streams.
* `calc_flux_base.f90` and `calc_steps.f90` are the main bottleneck. You should optimize them. `calc_flux_base.f90` calls `calc_*_kernel.f90`, `calc_*_kernel_internal.f90`, `calc_visc*.f90`, and `calc_div.f90` (runs once per RK stage, right after `calc_quantities_T_3D`, for `VISC_ORDER>2`).
* Roe scheme is not used. KEEP, SLAU, Hybrid schemes should be optimized.
* STZ case exists specifically to validate the z-direction halo exchange (COMMZ=True, BC_Z=True, Euler, 2nd order). The 6th-order NS path is exercised by SWLBLI with `COMMZ=True`; see the COMMZ section above.
* `calc_time_dev.f90.fypp`'s RK3 and RK4 loops (with and without COMMZ) call `set_bc` with the extra `(x, z, phi_l_gpu, phi_m_gpu, t_now)` arguments only when the case's `config.fypp` sets `BC_FORCING = True` (SWLBLI); every other case keeps the plain signature. Adding a case-specific argument to `set_bc` again must go through such a flag, or every other 3D case stops compiling.
* `RESTART=True` + `RESCALE=True` used to abort on startup: `pre_rescale`
  (`preprocess.f90.fypp`) requires `recal/Qm.dat`, but `write_Qm` — the only routine
  that wrote it — was never called from anywhere. `rescale_recv_send` now checkpoints
  the profile through `write_Qm_restart` every 1000 timesteps once rescaling has
  engaged, which is what makes a staged TBL run resumable. Note the resulting
  semantics: on restart `id_recal` is true, so `calc_mean` is skipped and `Qm` stays
  **frozen** at the restored profile for the whole run.
* Two known quirks in `calc_rescale.f90`, both benign, both left alone: the `jup` loop
  that sets `u99` has no `exit`, so `u99` is effectively `0.99*Um(ny)` (the freestream
  value you want anyway); and the `step` index passed from `calc_time_dev.f90.fypp` is
  `np*(t2-1)+t1` rather than `nt*(t2-1)+t1`, which makes the time column of
  `data/rescaling.d` a sawtooth but does not affect the running mean (`calc_mean` keeps
  its own counter).
* `calc_hybrid.f90`'s `calc_Ducros` shock sensor (used by SLAU/Roe/Hybrid schemes) recomputes its own full velocity-gradient tensor via always-2nd-order central differences, independent of `calc_div`'s higher-order `ux,vy,wz` — reusing `calc_div`'s output there would change the sensor's numerical values (mixed-order dilatation vs. vorticity) and needs a fallback for the edge band `calc_div` doesn't cover, so it hasn't been done; flagged as a known optimization candidate for SLAU/Hybrid NS/LES cases (currently: SBLI).
* Nsight Compute profiling captures live in `3D_solver/nsys_ncu/*.ncu-rep` (open with `ncu --import <file> --print-summary per-kernel`, or the `.csv` triage exports); use these before speculating about kernel performance.
* `preprocess.f90.fypp`'s `pre_calc` (the initial t=0 snapshot send) used to hardcode `real(4)`/`MPI_REAL4` for its flux-flat buffers regardless of `OUTPUT_PRECISION`, unlike `print.f90.fypp`'s `send_recv_for_print_even/odd3`, which already dispatched on the `io` kind parameter. This was never exercised because no case had set `OUTPUT_PRECISION=8` before EVC/BL/OS's tests needed it (float32 can't resolve a 6th-order convergence trend or tight Cf tolerances); fixed to use `real(io)` and branch on `io` like the rest of the pipeline. No existing case's behavior changes (`io` still resolves to 4 for all of them).
