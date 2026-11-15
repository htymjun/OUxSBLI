# CLAUDE.md — 3D_solver/

This file documents the 3D Cartesian solver. See the repository root `CLAUDE.md` for project overview, shared configuration reference (`config.fypp` variables, `mod_globals.f90`/`mod_constant.f90`), and tooling rules.

## fypp-preprocessed source files

| Template (in `src/` or `3D_solver/src/`) | Purpose |
|------------------------------------------|---------|
| `src/mod_constant.f90.fypp` | Kind-dispatch Fortran parameters from config values |
| `src/calc_weno.f90.fypp` | WENO5-Z face reconstruction (FP64, W32P64), called by `calc_muscl`'s `delta6` when `RECON='WENO'` |
| `3D_solver/src/calc_flux_base.f90.fypp` | Top-level flux dispatcher (convective + viscous) |
| `3D_solver/src/calc_time_dev.f90.fypp` | RK time-stepping orchestration |
| `3D_solver/src/calc_para.f90.fypp` | MPI halo exchange: x-direction pack/unpack and the COMMZ z-halo `start_exchange_z`/`finish_exchange_z` (`MPI='CPU'` host staging or `MPI='GPU'` CUDA-aware) |
| `3D_solver/src/calc_keep_kernel.f90.fypp` | KEEP convective kernel |
| `3D_solver/src/calc_slau_kernel.f90.fypp` | SLAU convective kernel |
| `3D_solver/src/calc_hybrid_kernel.f90.fypp` | Hybrid KEEP↔SLAU kernel |
| `3D_solver/src/calc_keep_kernel_internal.f90.fypp` | Interior-only KEEP variant |
| `3D_solver/src/calc_slau_kernel_internal.f90.fypp` | Interior-only SLAU variant |
| `3D_solver/src/calc_hybrid_kernel_internal.f90.fypp` | Interior-only Hybrid variant |
| `3D_solver/src/calc_visc2.f90.fypp` | 2nd-order viscous kernels |
| `3D_solver/src/calc_visc_high.f90.fypp` | 4th/6th-order viscous kernels (`VISC_ORDER`-controlled), boundary-aware (interior/fallback branch per point) |
| `3D_solver/src/calc_visc_high_internal.f90.fypp` | Interior-only 4th/6th-order viscous variant (used where a direction has no BC, i.e. periodic) |
| `3D_solver/src/calc_visc_cent.f90.fypp` | Shared `interp`/`diff`/`calc_tau_straight`/`calc_tau_cross` helpers, `#:include`d into both `calc_visc_high[_internal].f90.fypp` |
| `3D_solver/src/load_smem_visc_cent.f90.fypp` | Shared-memory loaders for `calc_visc_high[_internal]`'s `calc_Ev/Fv/Gv` |
| `3D_solver/src/preprocess.f90.fypp` | Device memory allocation helpers |

## Architecture

Source is split across directories; CMake/vpath merges them at build time:

- `src/` — shared templates and utilities: `mod_constant.f90.fypp`, `sbli.f90` (SWTBLI's two-domain main; the single-domain `main.f90` is in `3D_solver/src/`), `cpu_gpu_mpi.f90`, `calc_muscl.f90.fypp`, `calc_physical_quantities.f90`, `print.f90.fypp`, `set_coordinate.f90`, `set_compressible_bl.f90`
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
| STZ  | z-direction halo exchange validation: viscous Sod shock tube along z, same physical state as 2D_solver/ST (NS, KEEP, ORDER=4, COMMZ=True, `MPI='GPU'`, BC_Z=True) |
| SWLBLI | Oblique-shock / laminar flat-plate boundary-layer interaction, M=2.15, β=32.72°, Xsh=100 mm (NS, SLAU, ORDER=6, BC_X=BC_Y=True, periodic z). Runs with `RESTART=True`: the initial field is the converged 2D solution of `2D_solver/SWLBLI`, extruded in z by that case's `convert_2d_to_3d[_2].f90` (see `2D_solver/CLAUDE.md`), so a fresh clone needs that 2D run first. It is the only NS case that can run with `COMMZ=True` (`MPI='GPU'`) — see the COMMZ section below. `BC_FORCING=True`: its `set_bc` takes `(x, z, phi_l_gpu, phi_m_gpu, t_now)` for a blowing/suction strip (currently commented out inside `set_bc`) |
| SWTBLI | Incident-shock / turbulent boundary-layer interaction with large separation, M=2.28, θ=12° (grid generation 2; generation 1 was IUSTI-like θ=9.5°), Re_θ≈500 (NS, SLAU-family + WENO5-Z `W32P64`, ORDER=6, RESCALE=True). Two coupled domains run by `src/sbli.f90` on 4 ranks: a recycling-rescaling TBL generator (ranks 0/1) feeds an SBLI domain (ranks 2/3). Staged with `stage.sh`/`run_stage.sh`; see the SWTBLI section below |
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

### SWTBLI: two coupled domains

* **Layout** (`SWTBLI/mod_globals.f90`, lengths in δ = inlet δ99 = `blt` = 1.0724 mm).
  Grid generation 2 (current; stages G0, C1, C2 N, P N):
  * TBL generator: 11δ × 4.07δ, dx0 = 0.0403δ (nx1=274, ny1=128), recycling plane
    `nre2` at 9.95δ, `Qm` averaged over 9.35–9.95δ.
  * SBLI domain (1397×455×114, 72 M points, 20 GB on the GPU): inlet at the TBL plane
    `nre2-5` (9.75δ), inviscid impingement `x_imp` 31δ (`ximp_in`) downstream, uniform
    to `x_imp+xdown` (22δ), then an 80-cell outlet buffer (+3 %/cell); 18.1δ high,
    uniform dy=dx up to `y_str`=16δ, `gy_top`=+3 %/cell above.
  * Generation 1 (A1..B0): dx0 = 0.0604δ, 773×218, 12.2δ high, `ximp_in`=25,
    `xdown`=18, uniform up to 9.5δ and +5 %/cell above. Two problems in its
    schlieren/p'rms led to generation 2: grid-scale waves trailing the reflected
    shock (worst where it entered the coarse +5 %/cell layer), and the wave reflected
    at the top by the reflected shock coming back into the recovering boundary layer
    (inviscid landing estimate x_imp+28δ, inside the domain). With 18δ and θ=12° the
    estimate is beyond the outlet (`report` → `top_reflection`). The longer inflow
    (`ximp_in` 31) keeps the top entry of the incident shock, `x_top` = x_imp − 24.3δ
    at θ=12°, 6.7δ downstream of the inlet.
  * dz (4δ/108), `dyw0` and `dt` are the same in both generations: `dt` is set by
    the acoustic CFL of the first wall cell (x CFL is ~0.2 on the finer grid).
  * Both share dx, the y-grid function and the span, so the SBLI inlet columns
    1..3 are the TBL planes `nre2-5..nre2-3` verbatim. Rank 0 sends them through
    host memory (`set.f90`), not as device pointers.
* **Incident shock**: enters through the SBLI top at `x_top = x_imp − Ly2/tan β1`.
  `theta_deg` (0 = no shock) is ramped 4 → 7 → 9.5 over the A2 stages (generation 1)
  and 9.5 → 11 → 12 over G0/C1/C2 (generation 2). The top
  boundary, each part found necessary by looking at schlieren GIFs:
  * Upstream of `x_top` it is a Riemann-invariant condition on the nominal state 1. The
    outer flow there sits ~2 % above p0 (displacement of the growing layer). Do not
    rescale state 1 to that level: the measured level then feeds back into itself.
  * At `x_top` it blends (tanh, ~4 cells) into state 2. States 2/3 are rescaled
    isentropically by f = running mean of the top pressure upstream of `x_top` / p0
    (≈1.02). With nominal states the mismatch radiated weak waves behind the incident
    shock.
  * From `x_top+3δ`, while the running-mean v there is still below v2/2, it is a
    Riemann condition on state 2.
  * Once the reflected/separation shock has crossed the top, the boundary copies row
    ny−1 at the foot of the right-running characteristic, `x − dy/tan(μ+φ)` from the
    running-mean local state. Two other conditions were tried and both reflected the
    exiting shock more strongly (peak dp/dx·δ/p0 at y=9δ: 0.59 and 0.92 against 0.33):
    the former SWTBLI switch to state 3 decided cell by cell from the local v, and a
    Riemann condition on the running-mean local state.
* **SBLI inlet above the TBL domain** (`j ≥ ny1`): a running time mean (τ=10 δ/u0) of
  the spanwise-mean TBL row `ny1−1`. Copying the instantaneous spanwise mean injected
  column-uniform density pulses that convected through the domain as vertical fronts.
* **Scheme**: SLAU (HR-SLAU2) + WENO5-Z `W32P64`. `RIEMANN='LF'` was rejected: its
  first-order wall face dissipates tangential momentum at ~c·ρu, several times the
  viscous wall stress, and the wall shear collapsed within ~10 δ/U (TBL and SWTBLI alike).
* **Stages, generation 2** (`./run_stage.sh STAGE [N]`, or several in a row with
  `nohup ./chain.sh G0 C1 "C2 1" &`; archived in `runs/<stage>/`):
  * `regrid.py runs/B0 mod_globals.f90 recal_gen2_init` builds the initial field
    from the end of B0: linear (x, y) interpolation, x aligned at x_imp; the 6δ
    longer inflow is a linear stretch of the old one (factor 1.46, continuous, no
    seam); above the old top the measured three-state inviscid field with straight
    shocks, blended over 1δ. `Qm.dat` is interpolated in y. It prints Cf/p_w of the
    source and the result side by side.
  * T2 (100-step smoke test), G0 (θ=9.5°, 50 δ/U, flushes the interpolation),
    C1 (θ=11°, 150 δ/U), C2 N (θ=12°, 300 δ/U segments until x_sh and x_0 pass
    `stationarity`), P N (θ=12°, 75 δ/U statistics segments, VTK every 2.5 δ/U).
  * Cost on the RTX PRO 6000 pair: 0.48 s/step (SBLI rank is the bottleneck), i.e.
    3.6 h per 150 δ/U.
  * `auto_gen2.sh` (start with `CHAIN_PID=<chain.sh pid>` to wait for a running chain)
    continues C2 N until `track.json` reports both x_sh and x_0 stationary (at most 4),
    then P 1..4 and `wall report` of the last C2 + P segments against B0 into `runs/final/`.
  * Raising θ moves `x_top` downstream; for the first ~30 δ/U the upstream top-pressure
    window still sees the old incident shock, so f starts high (1.089 at C1) and relaxes
    with τ = 10 δ/U.
  * `post_stage.sh` (called by `run_stage.sh`) runs turb/track/report/wall/gif,
    then keeps only the last 3 snapshots in `data/` (an SBLI snapshot is ~1.5 GB) and
    the `recal/` of the last two generation-2 stages. Stage order: `runs/STAGES_GEN2`.
    `track --prev` continues the previous stage's `track.npz`, so the stationarity
    test spans stages whose VTK files are gone; `report`/`track --ref` take the
    no-shock baseline from `runs/A1b/analysis/ref_rows.npz` (`cache-ref`),
    interpolated in physical x onto any grid generation.
* **Stages, generation 1** (parameters in `runs/<stage>/mod_globals.f90`; `stage.sh`
  refuses them now; their `data/` was deleted after the analysis, `analysis/`,
  `gif/`, `recal/` remain; B0 keeps its `data/` as the θ=9.5° reference):
  * A1: 2δ span, no shock, from a composite turbulent profile with
    streak/vortex/noise perturbations. A1b continues it, with `Qm` re-seeded at the
    recycling plane, to bring the inlet back to δ99 = blt (Re_θ ≈ 500).
  * A2a/A2b: θ=4°/7°.
  * A2c N: θ=9.5°, 300 δ/U segments until the separation shock is
    statistically stationary.
  * `convert_2d_to_4d.py`: tiles the 2δ restart into 4δ and adds 1 % noise. The
    noise is needed because the half-span recycling shift maps a Lz/2-periodic
    field onto itself.
  * B0: 4δ, short decorrelation. Its end state is the initial field of the production run.
* **Analysis**: `python -m ouxsbli.analysis.swtbli_stats {turb|track|report|zcorr|gif|wall|cache-ref} runs/<stage>`.
  It reads `data/1` (TBL) and `data/3` (SBLI) with `vtr_raw` and covers the turbulence check,
  shock tracking with the stationarity test, the boundary diagnostics and the GIFs.
  The grid is rebuilt from `mod_globals.f90` by `make_x`/`make_y` (twins of `set.f90`).
  * `wall accumulate RUNDIR` → `analysis/wall_acc.npz` (sums of span-mean p_w, p_w²,
    τ_w and mean fields); `wall report NPZ... [--compare NPZ] [--ref-csv CSV]` merges
    segments and plots/writes Cf and p_w/p1 against Van Driest II (upstream), the
    free-interaction theory (F=4.22 separation, 6.0 plateau; Zukoski 1+M/2), the
    inviscid p3, other runs, and digitised literature data in the IUSTI
    normalisation x* = (x−x_imp)/L_int, p* = (p_w−p1)/(p3−p1). A plateau is reported
    only when dp_w/dx has an interior minimum inside the bubble.
  * `report` → `post_reflected`: grid-scale residual (deviation from a ≥0.3δ running
    mean, the same physical width on both grids) of mean p and p'rms 0.5–2.5δ behind
    the reflected shock, over the inviscid jump p3−p2; `top_reflection`: inviscid
    landing point of the wave reflected at the top. Before the generation-2 change
    `P["y_str"]` was parsed without the `*blt` factor, so the old
    `prms_ratio_top_to_9p5delta` (0.94) compared the top row with itself; the
    corrected `prms_ratio_top_to_ystr` is 0.31 for B0.
* `src/sbli.f90` (the two-domain main) was brought up to date: it now uses
  `set_block_3D`/`set_Jacobian_xy3`, node-local GPU selection (TBL pair on GPU0,
  SBLI pair on GPU1 when two are visible), and requires exactly 4 ranks.

## Convective Schemes (dispatched from `calc_flux_base.f90.fypp`)

| Scheme | Files | Best for |
|--------|-------|----------|
| KEEP (Kinetic Energy & Entropy Preserving) | `calc_keep_kernel.f90.fypp`, `calc_keep_kernel_internal.f90.fypp` | Smooth vortical flows (TGV, KHI) |
| SLAU (Simple Low-dissipation AUSM) | `calc_slau_kernel.f90.fypp`, `calc_slau_kernel_internal.f90.fypp` | Compressible turbulence with shocks (SWLBLI, SWTBLI) |
| Hybrid (KEEP↔SLAU via Ducros sensor) | `calc_hybrid_kernel.f90.fypp`, `calc_hybrid_kernel_internal.f90.fypp` | Mixed smooth/shocked regions |

Viscous discretization: `calc_visc2.f90.fypp` (Gaitonde & Visbal 2nd-order, `VISC_ORDER=2`) or `calc_visc_high.f90.fypp` + `calc_visc_high_internal.f90.fypp` (4th/6th-order, `VISC_ORDER=4` or `6`).

### WENO5-Z reconstruction and the LF flux (`RECON`, `WENO_PREC`, `RIEMANN`)

Shared by `2D_solver` and `3D_solver`; not wired into `3D_solver_curv`.

**Where the code lives.** No convective kernel was edited. `delta6(sensor, a(6), al, ar)` in
`src/calc_muscl.f90.fypp` already has WENO5's geometry: the face lies between `a(3)` and `a(4)`,
`al` comes from `a(1..5)` and `ar` from `a(2..6)`. With `RECON='WENO'`, its body calls
`weno5z_fp64` or `weno5z_w32p64` from `src/calc_weno.f90.fypp`. Every SLAU-family kernel that
calls `delta6` therefore reconstructs with WENO:
- `calc_slau_kernel[_internal]`;
- the fused `calc_slau_visc_kernel`;
- the SLAU branch of `calc_hybrid_kernel[_internal]`.

`RIEMANN='LF'` makes `id_slau` an `integer(8)`. That selects `LFflux`, a local Lax-Friedrichs
(Rusanov) flux in `src/calc_scheme_math.f90.fypp`, through the existing `interface SLAU` blocks.
`HR` (the wiggle detector) is ignored, and `c` is not scaled by a metric because the Cartesian
`Normal_*` are unit vectors. That is also why `LFflux` is deliberately absent from the
`3D_solver_curv` interfaces: `RIEMANN='LF'` fails to compile there instead of silently using a
non-unit normal.

`mod_constant.f90.fypp` stops the build when:
- `RECON='WENO'` is used without `ORDER=6`;
- `RECON='WENO'` or `RIEMANN='LF'` is used with `SCHEME='KEEP'`.

**Near walls (`BC_*=True`).** WENO applies only where `delta6` was used. Faces 2 and n−2 stay
MUSCL3rd via `delta4`, so `TVD` still matters there. Faces 1 and n−1 stay first order.

**W32P64 is not WENO_ILP's arithmetic verbatim.** The verbatim port narrows the stencil values to
FP32 and combines with `Σ w_k p_k`. That has an absolute error floor of about 6e-8·|v|, because
the separately rounded FP32 weights do not sum to 1. On exact cell averages it stalled at ~1e-7
(O(1) data) and 7.5e-3 (1e5 + 1e-2·sin). On 2D EVC it gave L2 ≈ 6e-8 at every grid. The
shipped version makes two changes:
- it builds the betas from **FP64 first differences narrowed to FP32**;
- it combines as **`p0 + w1(p1−p0) + w2(p2−p0)`**.

With both, it matches FP64 to the printed digits, at 5th order on both biases and on EVC. The
`calc_weno.f90.fypp` header has the measurements and the SASS counts. On the 2D `calc_slau_x_in`,
FP64 WENO takes 4096 instructions, 36 CALLs; W32P64 takes 2504, the same 8 CALLs as MUSCL.

FP32 rules carried over from WENO_ILP:
- the ratio cap;
- unbatched FP32 ratios;
- mirrored right-bias linear weights;
- the approximate reciprocal for FP32 only.

`__fdividef` compiles to `MUFU.RCP + FMUL` without range guards, because `-fast` already makes
FP32 flush-to-zero. Do not add a fast-math flag.

**Expected order on smooth flow.** The SLAU/LF kernels reconstruct primitive *point* values and
then evaluate a nonlinear flux. That leaves an O(h²) error whatever the reconstruction. On 2D
EVC at n=256, FP64 WENO + SLAU gives 3.9528e-10 and MUSCL + SLAU gives 3.9526e-10. Measured
orders over 64/128/256:

| configuration | order |
|---|---|
| WENO + LF | 3.47 |
| WENO + SLAU | 3.18 |
| MUSCL + SLAU | 1.79 |

The reconstruction alone is 5th order: see `ouxsbli/tests/test_weno.py`, which gates that
together with FP64 WENO + LF on EVC and the OS shock.

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
  ghost planes. Observation on STZ (not changed): it sets `BC_Z=True` yet
  `start_exchange_z` wraps the outer ranks periodically and `calc_conv_G_koff` always
  uses the `_in` z kernel.

### SWLBLI: the blowing/suction strip (`BC_FORCING`)

`BC_FORCING=True` only changes the `set_bc` *signature* (`(x, z, phi_l_gpu, phi_m_gpu, t_now)`); the strip itself is still commented out inside `SWLBLI/set.f90`. The phases are drawn once in `calc_time_dev` and `MPI_BCAST` so every rank sees the same `phi_l`/`phi_m`, and `t_now` is the step index, so all RK stages of a step share one forcing time (RK3 and RK4 alike). Four things to check before switching the strip on — all found by reading the code, none tested:

* `h_t` sums `sin(beta_force*t + 2π φ_m)` over `m = m_min..m_max` **without a `dble(m)` factor**, unlike `g_z` which does carry `dble(l)`; as written every temporal mode has the same frequency and the sum collapses to a single sine.
* `beta_force = 75000` is commented `(Hz)` but enters as `sin(β t)`, i.e. 75000 rad/s ≈ 11.9 kHz — multiply by 2π if Hz was meant.
* Spanwise modes `l = 12..24` need a span that resolves them: at the current `Lz = 4·blt` with `nz-6 = 134` planes per compute rank, `l = 24` has a wavelength of `Lz/24 ≈ 5.6·dz` on one compute rank — resolved, but only marginally for a 6th-order stencil.
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
* KEEP, SLAU, Hybrid schemes should be optimized (the unused Roe kernels were removed).
* STZ case exists specifically to validate the z-direction halo exchange (COMMZ=True, BC_Z=True; currently NS, KEEP, 4th order, `MPI='GPU'`). The 6th-order NS path is exercised by SWLBLI with `COMMZ=True`; see the COMMZ section above.
* `calc_time_dev.f90.fypp`'s RK3 and RK4 loops (with and without COMMZ) call `set_bc` with the extra `(x, z, phi_l_gpu, phi_m_gpu, t_now)` arguments only when the case's `config.fypp` sets `BC_FORCING = True` (SWLBLI); every other case keeps the plain signature. Adding a case-specific argument to `set_bc` again must go through such a flag, or every other 3D case stops compiling.
* `RESTART=True` + `RESCALE=True` used to abort on startup: `pre_rescale`
  (`preprocess.f90.fypp`) requires `recal/Qm.dat`, but `write_Qm` — the only routine
  that wrote it — was never called from anywhere (it has since been removed). `rescale_recv_send` now checkpoints
  the profile through `write_Qm_restart` every 1000 timesteps once rescaling has
  engaged, which is what makes a staged TBL run resumable. On restart the running
  mean now **resumes** from the checkpoint with the weight of 2000 steps (not the
  true sample count of the previous stage): `pre_rescale` loads `recal/Qm.dat` into
  the compute rank's device `Qm`, `calc_time_dev` starts the counter `step` at 2001 (the
  restored profile weighs 2000 steps), and `calc_mean`/`rescale_recv_send` no longer
  skip on `id_recal`. It used to stay frozen at the restored profile, which makes the
  inlet δ drift exponentially over a long multi-segment run: the outer scaling maps
  the recycled plane with the frozen δ99 of `Qm`, so every recycling pass multiplies
  the inlet δ by (growth from inlet to recycling plane) × (target δ / frozen δ).
  Related: `Qm` should be averaged at the recycling plane (`nre1` close to `nre2`), not
  over a wide window upstream of it, or the same ratio exceeds 1 and the inlet δ grows
  a few per cent per pass even without restarts (seen in SWTBLI; TBL still uses
  0.5–0.9 nx and has not been re-checked).
* Two known quirks in `calc_rescale.f90`, both benign, both left alone: the `jup` loop
  that sets `u99` has no `exit`, so `u99` is effectively `0.99*Um(ny)` (the freestream
  value you want anyway); and the `step` index passed from `calc_time_dev.f90.fypp` is
  `np*(t2-1)+t1` rather than `nt*(t2-1)+t1`, which makes the time column of
  `data/rescaling.d` a sawtooth but does not affect the running mean (`calc_mean` keeps
  its own counter).
* `calc_hybrid.f90`'s `calc_Ducros` shock sensor (used by SLAU/Hybrid schemes) recomputes its own full velocity-gradient tensor via always-2nd-order central differences, independent of `calc_div`'s higher-order `ux,vy,wz` — reusing `calc_div`'s output there would change the sensor's numerical values (mixed-order dilatation vs. vorticity) and needs a fallback for the edge band `calc_div` doesn't cover, so it hasn't been done; flagged as a known optimization candidate for SLAU/Hybrid NS/LES cases.
* Nsight Compute profiling captures live in `3D_solver/nsys_ncu/*.ncu-rep` (open with `ncu --import <file> --print-summary per-kernel`, or the `.csv` triage exports); use these before speculating about kernel performance.
* `preprocess.f90.fypp`'s `pre_calc` (the initial t=0 snapshot send) used to hardcode `real(4)`/`MPI_REAL4` for its flux-flat buffers regardless of `OUTPUT_PRECISION`, unlike `print.f90.fypp`'s `send_recv_for_print_even/odd3`, which already dispatched on the `io` kind parameter. This was never exercised because no case had set `OUTPUT_PRECISION=8` before EVC/BL/OS's tests needed it (float32 can't resolve a 6th-order convergence trend or tight Cf tolerances); fixed to use `real(io)` and branch on `io` like the rest of the pipeline. No existing case's behavior changes (`io` still resolves to 4 for all of them).
