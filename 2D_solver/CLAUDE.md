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

## Adding a New Test Case

1. `cp -r 2D_solver/OS 2D_solver/MYCASE`
2. Edit `mod_globals.f90` — grid size, physical parameters
3. Edit `set.f90` — grid generation, initial conditions, boundary condition calls
4. Edit `config.fypp` — VISC, SCHEME, ORDER, BC_X/Y/Z, etc.
5. Edit `calc.sh` — runtime args
6. `cmake -B build && cmake --build build -j`
