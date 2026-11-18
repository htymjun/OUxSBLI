# CLAUDE.md — 2D_solver/

This file documents the 2D solver. See the repository root `CLAUDE.md` for project overview, shared configuration reference, and tooling rules.

## Source Layout

A standalone 2D solver sharing the same convective/viscous kernels as the 3D Cartesian solver. Source layout mirrors the 3D structure:

- `2D_solver/src/` — shared 2D utilities (main, grid, BCs)
- `2D_solver/<CASE>/` — per-case config: `mod_globals.f90`, `set.f90`, `config.fypp`, `CMakeLists.txt`, `calc.sh`

The 2D solver uses CMake (like 3D); `2D_solver/src/calc_flux_base.f90.fypp` is preprocessed by the per-case CMakeLists.txt using the same `fypp -I<case-dir>` pattern.

## Cases

| Case | Description |
|------|-------------|
| BL   | Laminar flat-plate boundary layer (M=0.1, validated vs. Blasius) |
| DSL  | Double shear layer |
| EVC  | Euler vortex convection (grid-convergence study) |
| OS   | 2D oblique shock (M=2, θ=8°, SLAU, Euler) |
| SBLI | Oblique shock / laminar BL interaction (M=2.15, β=30.8°, Re_Xsh≈1e5) |
| ST   | Sod shock tube |
| SWLBLI | 2D precursor of `3D_solver/SWLBLI`: oblique shock / laminar BL interaction (M=2.15, β=32.72°, p_tot=100 kPa, Xsh=100 mm, 2700×480, SLAU) — see below |

## WENO5-Z and Lax-Friedrichs

`RECON='WENO'`, `WENO_PREC='FP64'|'W32P64'` and `RIEMANN='SLAU'|'LF'` work here exactly as in
3D. The SLAU-family kernels call the shared `delta6` / `interface SLAU`, so no 2D-specific code
exists. See the WENO section of `3D_solver/CLAUDE.md`.

`EVC` and `OS` carry the three lines in their `config.fypp`, so `Case()` can patch them.

## BL and SBLI: Shared Flat-Plate Setup

BL and SBLI share the same flat-plate setup: a uniform freestream over a
symmetry wall upstream of the leading edge and a no-slip wall downstream of it,
with the inlet column left frozen at its `set_init` state (the interior-only RK
kernels never touch i=1, and neither case defines an inlet BC). In SBLI that
frozen column is also the oblique-shock generator — above the incident-shock
trace it holds the analytic post-shock state. Both place the leading edge at
x = 0 via the `i_LE` constant in `mod_globals.f90`, so wall quantities can be
compared against theory and reference data without an origin offset.

Their **top boundaries differ, and must**. SBLI uses a Riemann-invariant far
field, which holds the incoming invariant fixed and therefore ties the boundary
pressure to the wall-normal velocity through the acoustic impedance:
`(p_b - p0)/q_inf = 2*v_b/(u0*M0)`. That is a non-reflection property for
transients, *not* a statement that p → p0 in steady state, and the `1/M0` makes
it unusable subsonically — at BL's M=0.1 it converted the displacement-induced
`v_b` into a spurious 0.07·q favorable pressure gradient along the plate,
accelerating the edge flow 3% and inflating Cf by >10% at the trailing edge. BL
therefore imposes `p = p0` at the top and takes everything else from the
outgoing characteristics (entropy and tangential momentum upwinded on the sign
of v). Do not "unify" the two BCs; at SBLI's M=2.15 there is no 1/M0
amplification and the Riemann form is fine there.

For the same reason BL's **outlet** imposes `p = p0` and extrapolates only
density and momentum. Extrapolating all four conservatives — which every other
2D case still does — supplies zero conditions where subsonic outflow needs
exactly one, so the exit pressure floats; in BL it settled 0.008·q below p0 and
dragged a favorable gradient ~25 mm back up the plate. A zero-incidence flat
plate has dp/dx = 0, so p = p0 is the correct single condition. Anything
subsonic added to `2D_solver/` should copy BL's two boundaries, not OS/ST's.

A closing caution on this case: **do not switch BL to `SCHEME='KEEP'`.** It
looks like the right choice (no shocks, and no upwind dissipation to thicken a
laminar layer) but it diverges at t ≈ 2.5 ms — with `TVD='none'` nothing damps
the mode seeded by the abrupt slip→no-slip switch at the leading edge. Lowering
dt does not help: at dt = 3e-9 (CFL_y = 0.011) it fails at the same *physical*
time as at 2e-8, so the growth is time-step-independent and this is the scheme,
not a CFL violation. SLAU is stable throughout. SLAU's dissipation was never
what inflated Cf here; the boundary conditions were.

Once the BCs stopped coupling pressure to a tall domain, the domain no longer
needed to be tall: `Ly` dropped 30→15 mm and `ny` 97→49 in the same ratio (so
`dy_wall`, and thus near-wall accuracy where Cf is measured, is unchanged), and
`dt` rose to the largest value that stayed clean in short probes with margin
(4e-8, CFL_y≈0.58; probes ran clean to at least CFL_y≈1.16 before a 2x safety
factor was applied). Time/step is dominated by fixed launch overhead (measured
~446 µs of ~613 µs at this grid), so cutting `endT` and raising `dt` — not
shrinking the grid — is what actually saves wall-clock time; the production run
went from ~26 min to ~2.5 min this way with no change to `SCHEME`, `ORDER`,
`VISC_ORDER`, or `BC_X`/`BC_Y`, i.e. the same boundary-aware (non-`_in`) kernels
still run.

Digitized reference data for SBLI (Moro et al., Degrez et al., Vila-Perez et
al.) lives in `2D_solver/SBLI/ref/`, together with the `Cf_all.py` / `Cp_all.py`
overlay scripts. Note that `.gitignore` excludes `2D_solver/*/data*` and `*.dat`
globally; `ref/` is kept tracked by an explicit `!2D_solver/*/ref/**` negation.

## SWLBLI: 2D precursor for the 3D warm start

`2D_solver/SWLBLI` is the 2D run whose converged field initialises
`3D_solver/SWLBLI` (which is built with `RESTART=True`). It started as a copy of
SBLI and keeps the flat-plate setup above, but it is **not** the validated SBLI
case — different shock angle, Reynolds number, domain and top boundary — so it has
no reference data and no test. SBLI stays at the Moro et al. conditions that
`ouxsbli/tests/test_sbli.py` checks; do not move SWLBLI parameters back into it.

Differences from SBLI that matter:

* The oblique shock enters through the **top** boundary instead of the frozen
  inlet column (whose post-shock override is commented out, so `i=1` holds plain
  freestream). `set_bc` splits the top row by `x + Ly/tan(beta)` against `Xsh`:
  upstream of the shock trace a Riemann far field with the *freestream* exterior
  state; on the trace (and `x < 0.09` m, hard-coded) the analytic post-shock state
  as a Dirichlet value; downstream zeroth-order extrapolation.
* tanh wall clustering `s = 3.0` (SBLI: 1.6).

Warm-start workflow:

1. Build and run this case; its final state is written to `recal/Q00001.dat`
   (physical conservatives — `main.f90` multiplies `Q/J` back by the Jacobian — as a
   raw stream of shape `(nx, ny, 4)`).
2. Compile and run `convert_2d_to_3d.f90` (spanwise sine disturbance in `v`) or
   `convert_2d_to_3d_2.f90` (round-off-level random noise) in this directory. Both
   copy every z plane, insert `rho*w = 0`, and write one `recal/Q_3d0000N.dat`
   per compute rank. Their `nx`, `ny`, `nz_new`, `n_compute_ranks` (and, in
   `convert_2d_to_3d.f90`, also `Lz` and the grid/flow constants) are hard-coded and **must match
   `3D_solver/SWLBLI/mod_globals.f90`** — the file is a raw dump with no header, so
   a mismatch is not detected, it just scrambles the field. As committed,
   `convert_2d_to_3d_2.f90` has `nz_new = 35` while the 3D case has `nz = 140`.
3. Copy `recal/Q_3d0000N.dat` to `3D_solver/SWLBLI/recal/Q0000N.dat` (the name the
   3D `main.f90` reads) and run the 3D case.

## Precision modes (`PREC`)

`PREC = 'FP64'` (default), `'DF'` or `'FP32'` in `config.fypp` (or `Case(prec=...)`).
`src/prec.fypp` defines the two types every template uses:

| | FP64 | DF | FP32 |
|---|---|---|---|
| `RT` (state, fluxes, primitives) | real(8) | `type(fltflt)` | real(4) |
| `WT` (metrics, coefficients, constants) | real(8) | real(8) | real(4) |

`set.f90` and `mod_globals.f90` are not fypp templates; they see the same choice through cpp:
`REAL_T` (= `RT`), `REAL_W` (= `WT`), `PREC_DF`, `PREC_FP32` (set in `2D_solver/CMakeLists.txt`).
FP32 also needs `sp = 4` in `mod_globals.f90` (the `#ifdef PREC_FP32` lines every case carries);
`mod_constant` refuses to compile otherwise.

Rules for code that runs on the device (they keep FP64 bit-identical and FP32 free of FP64 instructions):

* No `!$cuf kernel do` loop may touch a DF value beyond a plain copy: nvfortran 25.5 cannot compile
  fltflt arithmetic, `dble()` or a defined assignment there. Boundary conditions are therefore
  hand-written `attributes(global)` kernels (every case's `set_bc`), and `calc_div` calls the
  device subroutine `div4` from its CUF loop. The cyclic copies in `set_bc_common` stay CUF loops
  (DF passes the state as 8-byte words through `IGNORE_TKR`).
* No real(8) constant or literal in FP32 arithmetic. Templates wrap device code in the `fp` block of
  `src/prec.fypp` (turns `1.d0` into `1.e0` when `PREC = 'FP32'`). The boundary kernels in `set.f90` declare
  their constants as `REAL_W, parameter` locals (`rho0 = g_rho0`, `half = 0.5d0`, ...) and
  leave mixed real(8) expressions only on the host, for example SWLBLI's top-boundary zones,
  which are computed once as index bounds `i_sh` and `i_09`.

Accuracy against FP64 (40 steps, `prec2d` runs; `ouxsbli/tests/test_prec_2d.py`): see
`docs/precision_modes.md`. HR-SLAU2 configurations differ at 1e-9..1e-3 under DF and FP32 for the
same reason as in 3D (the 0/1 wiggle switch); compare them with `SLAU_VARIANT = 'SLAU'`.

## Constant viscosity, cell sources, raw state files

Three `config.fypp` switches for periodic runs (all off by default; the generated code of every
existing case is unchanged with them off):

* `VISCOSITY = 'Const'`: `calc_quantities_T_2D` sets `mu = mu0` instead of Sutherland's law. The
  case's `mod_globals.f90` defines `real(8), parameter :: mu0`.
* `SOURCE = True` (needs `VISC = 'NS'`, `BC_X = BC_Y = False`): `src/calc_source.f90.fypp`.
  S_m = rho (f - fric u), S_E = u.S_m - q_cool with q_cool = cool (<rho e> - p0 / (gamma - 1)).
  `calc_source_term` runs after every `calc_EF` (it uses the primitives of that stage) and stores
  `Sr = dt / J * S` in the layout of `Rs`; `calc_R` subtracts it with the stage weight of the fluxes
  (`coef` / 1 / `coef3`). q_cool is uniform and is evaluated once per step, at the first stage, from
  the interior mean pressure (block partial sums on the device, final sum on the host in real(8)).
  The case provides `set_force(nx, ny, x, y, fx, fy)` in module `set` (body force per unit mass,
  real(8), all cells) and `fric`, `cool`, `p0` in `mod_globals.f90`.
* `OUTPUT_RAW = True`: at step 0 and every output, the even rank writes `data/Q#####.dat` (float64
  `Q(N1,N2,4)` of the interior cells, conservative, x fastest, no ghost layers) and
  `data/Q#####.dat.txt` (`t=... step=... N=N1`); no VTK is written. `read_raw` in
  `src/io_raw.f90.fypp` reads such a file into `Q(nx,ny,4)` and fills the ghost layers periodically,
  so a case's `set_init` can restart from it with another `ORDER`.

Checked (64^2, 40 steps): friction decay and cooling relaxation of uniform states match the exact
RK solution to 6e-15 in FP64 and DF; a forced NS run agrees between DF and FP64 to 3e-14 (rho) and
2e-14 relative (rho E); the FP32 SASS has no FP64 instruction.

## Thermal fluctuations (`FLUCT`)

`FLUCT = True` (needs `VISC = 'NS'`, `VISCOSITY = 'Const'`, `VISC_ORDER = 2`, `RK = 3`, periodic x and y)
adds the Landau-Lifshitz stochastic stress and heat flux on the faces (`src/calc_fluct.f90.fypp`,
after `calc_EF`, before the source term) and, with `SOURCE = True`, the friction noise in
`calc_source_term`. The case defines `kB` (real(8), Boltzmann constant in code units) and `seed`
in `mod_globals.f90`. The noise covariance is that of the linearised 2nd-order viscous operator
(discrete fluctuation-dissipation); see the header of `calc_fluct.f90.fypp`.

`src/mod_rng.f90.fypp`: splitmix64-hashed counters (stream key, field, periodic cell/face index),
Box-Muller in FP32, two streams per step combined as `xi_A + beta_s xi_B` for the three RK3 stages
(Delong et al. 2013). Samples and amplitudes are FP32 in every `PREC`, so FP64, DF and FP32 runs
draw identical noise. `rng_step0` is the global index of the run's first step (set it in
`set_init` when restarting). Compiler notes (nvfortran 25.5): a kernel cannot take a derived-type
`value` argument (the stage noise goes as `keyA, keyB, beta`); `gauss` / `stage_gauss` take their
arguments by reference (with `value` dummies the inlined `stage_gauss` returned one sample for every
index at -O2); the 24-bit integers go through `integer(4)` variables before `real()` (otherwise
`I2F.S64` in the SASS).

Checked: the 64-bit integer sequence equals an independent numpy implementation bit for bit and the
samples agree to FP32 rounding (7e-7) with the right moments; a gas at rest (64^2, 4e4 steps) reaches
equipartition, mode-averaged (k >= 5) |u|^2, |v|^2, |T|^2, |rho|^2 at 0.997-1.007 of
kB T / (rho dV), kB T^2 / (rho cv dV), rho^2 kB T / (p dV) in FP64, DF and FP32; DF with noise
agrees with FP64 to 3e-14 after 40 steps; the FP32 SASS has no FP64 instruction.

## Adding a New Test Case

1. `cp -r 2D_solver/OS 2D_solver/MYCASE`
2. Edit `mod_globals.f90` — grid size, physical parameters
3. Edit `set.f90` — grid generation, initial conditions, boundary condition calls
4. Edit `config.fypp` — VISC, SCHEME, ORDER, BC_X/Y/Z, etc.
5. Edit `calc.sh` — runtime args
6. `cmake -B build && cmake --build build -j`
