# Quasi-2D 3D-solver analogs of test_bl.py, test_evc.py, test_os.py

## Context

`test_bl.py`, `test_evc.py`, and `test_os.py` validate the solver against exact
analytical answers (Blasius, vortex-returns-to-IC, Rankine-Hugoniot), but only
through `2D_solver/`. The 3D Cartesian solver (`3D_solver/`) shares most of its
convective/viscous kernels with the 2D solver but also has code the 2D solver
structurally cannot exercise: the G-flux (z-direction), `calc_div_wz`, 3D grid
metrics/Jacobians, 3D `calc_para` x-halo exchange, and the 3D VTK writer. None
of the existing 3D tests (`test_etgv.py`, `test_corn.py` for the curvilinear
tree) validate a wall-bounded viscous case or a shock reflection against a
closed-form answer in 3D. The goal here is to close that gap cheaply: extrude
each 2D case's real physics uniformly through a thin, periodic spanwise `z`,
run it through the full 3D solver, and check that (a) it still matches the
same analytical answer and (b) the field comes out genuinely uniform in `z` —
which is itself a direct test of the 3D-only code paths above.

This pattern already exists once in the repo, just not under `3D_solver/`:
`3D_solver_curv/CORN` + `test_corn.py` do exactly this for a curvilinear
oblique-shock case (`nz=4`, pick one interior z-slice, reuse the 2D analytical
helpers). `3D_solver/SBLI` independently shows the exact mechanism needed to
combine real x/y boundary conditions with a periodic-only z: apply real x/y
BCs, then call `set_bc_cyclic_z` (`3D_solver/src/set_bc_common.f90:331-349`).
`3D_solver/IVST` (`ny=nz=7`) and `3D_solver_curv/NACA` (`nz=17`, comment
"quasi-2D spanwise") confirm thin-periodic-transverse is an established idiom.
`ouxsbli/tests/utils/vtk_reader.py` and `ouxsbli/analysis/wall.py` already
handle real 3D fields with no `nz=1` assumption — `wall_coeffs_from_vtr`
already slices `[0,:,:]` before doing 2D math, so it works unmodified on a
`(nz,ny,nx)` array as long as the field is genuinely uniform in z.

Three new cases: `3D_solver/BL`, `3D_solver/EVC`, `3D_solver/OS`, each a
"quasi-2D" extrusion of its 2D counterpart, plus three new test files.

## Shared design decisions (apply to all three cases)

- **`nz = 7`, fixed, never patched by a test.** `set_grid_cyclic6_3D`
  (`src/set_coordinate.f90:284-337`) and `set_bc_cyclic_z` need `nz = 2·ghost+1
  = 7` minimum for a 6th-order periodic z-stencil (all three 2D sources use
  `ORDER=6`). Verified by hand-tracing `set_bc_cyclic_z`'s ghost-copy indices at
  `nz=7`: `QJ(1)=QJ(2), QJ(2)=QJ(3), QJ(3)=QJ(4)[interior], QJ(5)=QJ(4),
  QJ(6)=QJ(5), QJ(7)=QJ(6)` — some of these read *other ghosts*, not the
  interior plane directly, but by induction this is harmless: if the field is
  already uniform in z before the call (true from `t=0` since the IC is
  k-independent), every value being copied is already the same number, so it
  stays uniform after. This only holds because the field is genuinely
  z-invariant; it would not be a correct periodic BC for a real varying field
  at this `nz`.
- **`nranks=2`** (1 GPU compute rank + 1 writer rank), matching
  `3D_solver/{ETGV,IVST,NSTGV}/calc.sh`. `nz=7` is far too thin to need
  z-domain-decomposition (`COMMZ=True`, `nranks=4`, à la STZ/SBLI) — that would
  be solving a different problem than intended here.
- **`mod_globals.f90` must declare `nre1, nre2, rerank, start_rescale`** even
  with `RESCALE=False` — `calc_rescale.f90` is compiled unconditionally into
  every 3D case and does `use mod_globals, only : nre1, nre2, rerank, ...,
  start_rescale` at file scope. Set them to the inert values `IVST`/`ETGV` use:
  `nre1=1, nre2=nx, rerank=0, start_rescale=0`. Omitting them is a link error,
  not a runtime issue.
- **`Q` in `set_init`/`set_bc` is `Q(nx,ny,nz,5)`, component *last*** — verified
  against `main.f90`'s call and `IVST`/`ETGV`/`SBLI`'s own `set.f90`. This is
  *not* the `Q(nx,5,ny,nz)` shape CLAUDE.md's "Conservative Variable Layout"
  section describes — that shape is the internal packed device buffer used
  inside `calc_steps_smem.f90`/`preprocess.f90.fypp`, never a case's own arrays.
  Getting this backwards silently scrambles the IC.
- **`set_bc` must not branch on `myrank`** (unlike SBLI's two-block dispatch) —
  with `nranks=2`, `set_bc` only ever runs on the compute rank; `myrank` is an
  unused dummy argument, as in `IVST`/`ETGV`/`NSTGV`.
- **`CMakeLists.txt`/`calc.sh`**: copy from a simple single-block case
  (`3D_solver/NSTGV` for BL/OS's real-x/y-BC pattern, `3D_solver/ETGV` for
  EVC's fully-periodic pattern) — **not** `3D_solver/SBLI`, whose
  `CMakeLists.txt` overrides `CASE_MAIN_SOURCE` for its two-block driver, which
  none of these three need. `calc.sh` is `mpiexec -n 2 ./build/a.out`.

## `3D_solver/BL` — quasi-2D laminar flat-plate boundary layer

Ports `2D_solver/BL` (M=0.1, `p=p0` top/outlet BC, `SLAU`, `ORDER=6`,
`BC_X=True, BC_Y=True`) plus `BC_Z=False`.

**`mod_globals.f90`**: copy every physical/geometric constant from
`2D_solver/BL/mod_globals.f90` verbatim (`R, gamma, M0, p_tot, T0, Pr, p0, rho0,
u0, rf, Taw, i_LE = 2*nx/12+1`). Add `Lz = 1.d0*blt`, `nz = 7`. `Lz`'s absolute
value is physically inert (flow is z-uniform, `w≡0`) but should not be shrunk
to microns: the z-direction still carries an acoustic eigenvalue `±c0≈340m/s`
in the SLAU flux, so `CFL_z = c0·dt/Lz` — at `Lz=1mm, dt=4e-8s` this is
`≈0.0136`, far looser than the y-direction's `≈0.58`, which is the intended
relationship (z must not become the binding stability constraint). Thread
blocks: `threadsE=dim3(32,1,1)`, `threadsF=dim3(32,4,1)`,
`threadsG=dim3(32,1,4)` and matching `Ev/Fv/Gv`, `threads=dim3(32,4,1)` —
copied from `3D_solver/SBLI`/`NSTGV`'s NS/`BC_X=True` precedent.

Grid size / cost is **unmeasured** — the 2D BL cost-reduction work established
that 2D per-step time is launch-overhead dominated, but the 3D port adds
whole-new kernel launches (`calc_slau_z_in`, `calc_div_wz`) with no 2D
analogue, so 3D per-step cost cannot be extrapolated from the 2D number. Do
this in two passes, mirroring exactly how the 2D BL cost pass proceeded:
1. **Correctness pass**: `nx=65, ny=25, nz=7`, `endT` cut to ~5·`np`·`dt`
   worth of steps. Confirms build/run/VTK-output/z-uniformity before spending
   time on the real grid.
2. **Physics pass**: `nx=257, ny=49` (2D BL's validated values — same
   `dy_wall`, same Cf-window point count that `test_bl.py`'s tolerances were
   tuned against), `nz=7`, `dt=4.d-8`, `endT=1.d-2`. Measure wall-clock; if too
   slow for a routine (non-`slow`-marked) test, apply the 2D BL playbook in
   order: raise `dt` toward the largest value with margin (the y-direction CFL
   analysis is unaffected by `nz`), shorten `endT` once the quasi-steady check
   converges sooner, and only as a last resort shrink `ny`/`Ly` together (as
   2D BL's own cost pass did) to hold `dy_wall` fixed.

**`set.f90`** — `set_grid`: x uniform + `i_LE` placement and y tanh-stretch
(`s=2.4`) copied verbatim from 2D BL; z uniform, matching
`set_grid_cyclic6_3D`'s own convention (`dz1=Lz/(nz-6)`, `z(k)=dz1*(k-4)`, so
`k=4` is the sole interior center). `set_init`: uniform freestream + wall row
zeroed for `i>=i_LE`, identical to 2D BL, looped over `k=1,nz` (IC is
k-independent by construction — this is what makes the z-uniformity check
meaningful). `set_bc`: copy the 2D BL outlet/top/wall math exactly (`p=p0`
outlet, `p=p0` top with entropy/tangential-momentum upwinding on `sign(v)`,
symmetry-then-no-slip wall at `i_LE`), adding a `k` loop (`!$cuf kernel do(2)`
over `(i,k)` or `(j,k)` in place of 2D's `do(1)` over `i`/`j` alone), then one
call to `set_bc_cyclic_z(nx, ny, nz, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)` at the end.
`Jacobian` stays 2D (`(nx,ny)`), identical at every `k`, so no BC math itself
changes — only the loop nesting.

**`config.fypp`**: `VISC='NS'`, `SCHEME='SLAU'` (KEEP diverges on this case
regardless of `dt` — do not use it, same reason as 2D BL), `ORDER=6`,
`TVD='none'`, `BC_X=True`, `BC_Y=True`, `BC_Z=False`, `RK=3`. Leave
`VISC_ORDER` unset (defaults to `ORDER=6`) — do not override to 4 "to widen the
z-interior band"; per CLAUDE.md's Notice, `VISC_ORDER=4` for an NS build is an
under-tested branch and should not be adopted opportunistically here.

**Known, understood limitation — not a bug, must shape the z-uniformity
check.** With `nz=7`, `calc_div.f90`'s `calc_div_ux_6_in`/`calc_div_vy_6_in`
only populate `k=4` (the sole plane satisfying `io_v+2<=k<=nz-io_v-1` for
`VISC_ORDER=6`). `calc_visc_high.f90.fypp`'s `calc_Ev6`/`calc_Fv6` handle this
safely (their `else`/fallback branch at other `k` recomputes the x/y viscous
terms directly from `Q` via a different, lower-order local formula rather than
reading the invalid-there precomputed arrays) — so there is no memory-safety
issue — but it means **the x/y viscous flux is computed by a genuinely
different numerical formula at `k=4` than at `k∈{2,3,5,6}`**. Since BL's x/y
viscous terms (`ux`, `vy`) are NOT physically zero (real BL development/
profile), this is expected to produce a small, benign, non-zero difference
between k-planes — unlike the pure z-derivative terms (`wz`, cross-derivatives
of `w`), which stay exactly zero at every k regardless of formula, since `w≡0`
identically. Consequence: **use `k=4` (Fortran, `k_py=3`) as the canonical
plane for every Cf/profile comparison against Blasius** — it is the one plane
that receives the same `VISC_ORDER=6` treatment 2D BL was validated against —
and treat the z-uniformity check for this case as a *bounded, not
near-machine-precision* check (see below).

## `3D_solver/EVC` — quasi-2D Euler vortex convection (grid convergence)

Ports `2D_solver/EVC` (fully periodic, `BC_X=False, BC_Y=False`) plus
`BC_Z=False` — the simplest of the three since it needs none of BL/OS's real
boundary-condition porting.

**`mod_globals.f90`**: copy every constant from `2D_solver/EVC/mod_globals.f90`
verbatim (`Lx=Ly=0.1, gamma=1.4, Pr=0.71, R=287.15, M0=0.05, beta=1/50,
theta=0, Rc=0.005, p0=1e5, T0=300, u0=M0·sqrt(gamma·R·T0), rho0=p0/(R·T0),
CFL=0.03, dt=CFL·Lx/((nx-1)·u0)`). Confirmed the period is `T = 2·Lx/u0` (two
convective periods) — `test_evc.py`'s own docstring says "one full period,
`T=Lx/u0`", which is wrong against the actual `mod_globals.f90` formula; carry
the *real* formula (`2·Lx/u0`) into the 3D case and don't propagate the
docstring's error. Add `Lz` (inert — no z-velocity, no viscous kernels at all
since `VISC='Euler'`; any positive value works, e.g. `Lz=0.01`), `nz=7`.
Thread blocks: `dim3(32,1,1)/(32,4,1)/(32,1,4)` for E/F/G, matching
`3D_solver/ETGV`.

**`set.f90`**: `set_grid` is one call, `set_grid_cyclic(id_accuracy, nx, ny,
nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)` — no hand-written grid code needed.
`set_init`: copy 2D EVC's vortex superposition formula
(`ex=exp(-0.5·((x-xc)²+(y-yc)²)/Rc²)`, `T=T0-0.5·(u0·beta)²/Cp·ex²`,
`rho=rho0·(T/T0)^(1/(gamma-1))`, swirl `du,dv`, rotate by `theta`) verbatim,
looped over `k=1,nz` with `Q(i,j,k,4)=0` (no spanwise velocity) — the formula
never references `z(k)`, so it produces bit-identical `Q` at every `k` by
construction. `set_bc` is one call, `set_bc_cyclic(id_accuracy, nx, ny, nz,
Q_1..Q_5)`.

**`config.fypp`**: `VISC='Euler'`, `SCHEME='KEEP'`, `ORDER=6`, `TVD='none'`,
`RK=4`, `BC_X=False, BC_Y=False, BC_Z=False`, `OUTPUT_PRECISION=8` (needed for
the same reason 2D EVC's test needs it — single precision can't resolve a
4th/6th-order convergence trend). `SCHEME`/`ORDER` here are just build
defaults; the test overrides both per-run via `Case(scheme=..., accuracy=...)`
exactly as `test_evc.py` already does.

**Ghost-trim interaction with `nz=7` at each tested `accuracy`**: verified from
`set_grid_cyclic{2,4,6}_3D`'s ghost widths — `accuracy=6` gives 1 real z-plane
at `nz=7` (`k=4` only), `accuracy=4` gives 3 (`k=3,4,5`), `accuracy=2` gives 5.
**`nz=7` works unmodified for every `accuracy` the test uses (4 and 6) — do
not vary `nz` per accuracy and do not let `Case()` patch it.** When trimming
ghosts for the L2 error (same `accuracy//2` trim `test_evc.py` already applies
to x/y), apply the identical trim to the z axis: `Rho[g:-g, g:-g, g:-g]` where
`g=accuracy//2` (array order `(Nz,Ny,Nx)` from `getScalar`). This always
leaves `>=1` z-plane; take `rho1[0]` after trimming (don't average over z — it
adds nothing since all real planes are bit-identical by construction here,
`VISC='Euler'` means no viscous fallback-formula effect exists to average
away).

## `3D_solver/OS` — quasi-2D oblique shock + wall reflection

Ports `2D_solver/OS` (M=2, incident shock β=37.2°, reflects off a slip wall at
y=0, `Euler`, `SLAU`, `ORDER=6`, `BC_X=True, BC_Y=True`) plus `BC_Z=False`.
Structurally like BL (real x/y BC + `set_bc_cyclic_z`) but simpler physics
(`Euler` — no viscous kernels, no fallback-formula nuance).

**`mod_globals.f90`**: copy every constant from `2D_solver/OS/mod_globals.f90`
verbatim (`R=287.03, gamma=1.4, Pr=0.72, M0=2.0, p_tot=100e3, T_tot=295,
beta=π·37.2/180, Lx=5·blt, Ly=2·blt, nx=257, ny=129, dt=3e-9, endT=0.1e-3,
np=1`, and the full incident-shock Rankine-Hugoniot chain giving `rho2, p2,
ux, uy` — this is the Fortran twin of `oblique_shock.py`'s functions). Add
`Lz=1.d0*blt`, `nz=7`. This case is Euler-only and already cheap in 2D
(`endT=0.1e-3` is tiny) — a 3D cost-reduction pass is not expected to be
necessary, but still measure it rather than assume.

**`set.f90`** — `set_grid`: plain uniform x/y (copied from 2D OS, no
stretching), z uniform matching `set_grid_cyclic6_3D`'s convention (same as
BL). `set_init`: uniform pre-shock freestream everywhere, then the
post-incident-shock state `(rho2, rho2·ux, rho2·uy, p2)` Dirichlet-injected on
the top row (`j=ny`) for `i >= 0.1·nx` onward — copied from 2D OS, looped over
`k`. `set_bc`: top row re-applies the pre/post-shock Dirichlet split every
step, bottom wall (`j=1`) is a slip/reflecting wall (`QJ_3(i,1,k) =
-QJ_3(i,2,k)`), left (`i=1`) Dirichlet freestream inlet, right (`i=nx`)
extrapolated outlet — all copied from 2D OS with a `k` loop added, then one
call to `set_bc_cyclic_z` at the end.

**`config.fypp`**: `VISC='Euler'`, `SCHEME='SLAU'`, `ORDER=6`, `TVD='hybrid'`,
`BC_X=True, BC_Y=True, BC_Z=False`, `RK=3`.

**No viscous-fallback nuance here** (unlike BL): `VISC='Euler'` means no
`calc_visc_high`/`calc_div` kernels exist in this build at all — only
`calc_slau_x/y` (along-sweep, no z-cross-term dependency) and `calc_Ducros`
(fixed 2nd-order formula, no interior/fallback split) touch the field. Every
k-plane should receive genuinely bit-identical IEEE-754 arithmetic every RK
stage.

## Test files

New files `ouxsbli/tests/test_bl_3d.py`, `test_evc_3d.py`, `test_os_3d.py`,
mirroring their 2D counterparts' structure/tolerances closely, with these
differences:

- **Workdirs**: `tmp/bl3d`, `tmp/evc3d_{scheme}{accuracy}_{nx}`, `tmp/os3d` —
  distinct from the 2D tests' `tmp/bl`, `tmp/evc_{scheme}{accuracy}_{nx}`,
  `tmp/os` (already present on disk) so `case.build()`'s `rmtree`/`copytree`
  can't collide and `pytest -k 3d` can isolate the new tests.
- **Canonical z-slice**: pick `k_py = nz//2 = 3` (Fortran `k=4`, the sole true
  6th-order-treated interior plane) for every physics comparison (Cf, velocity
  profile, oblique-shock region means) — reuse `wall_coeffs_from_vtr`/
  `edge_state`/`fprime`/the `oblique_shock.py` helpers completely unmodified,
  just index `rho[3,:,:]` etc. instead of `rho[0,:,:]`, mirroring
  `test_corn.py`'s existing precedent of using an interior index (there,
  `rho[1,...]`), not `k=0`.
- **New assertion in every file — the actual point of this work — z-uniformity
  of the raw field.** Two tiers, justified by the mechanism above:
  - **OS and EVC (Euler, no viscous kernels)**: tight, near-machine-precision
    check across *all* k, both ghost and interior:
    `assert np.abs(field[k]-field[k_ref]).max()/np.abs(field[k_ref]).mean() <
    1e-10` for every `k`. Use `output_precision=8` in `Case(...)` so the
    check isn't muddied by the float32 write cast.
  - **BL (NS, viscous fallback-formula effect)**: two-part —
    (a) planes sharing the same fallback formula (`k_py∈{1,2}` vs `{4,5}` in
    0-based, i.e. Fortran `k=2,3` and `k=5,6`) should still agree to
    `~1e-10`, since they share one code path; (b) the true-interior plane
    (`k_py=3`) vs. a fallback plane may differ by more — assert this
    difference is small **relative to the physics tolerances**, not
    near-machine-precision: start with a placeholder (e.g. `< 0.1 ×
    CF_MEAN_TOL`), measure the actual value once the case runs, and tighten
    to ~1.5× measured, exactly the convention `test_bl.py` already uses for
    `CF_MEAN_TOL`/`EDGE_ACCEL_TOL`. Document *why* this tier is looser
    (viscous fallback formula, not a bug) directly in the assertion's failure
    message so a future reader doesn't mistake it for a regression.
- **EVC's ghost trim** extends to the z axis (`g=accuracy//2` on all three
  axes of the `(Nz,Ny,Nx)` array from `getScalar`), then take the first
  (only, at `accuracy=6`) remaining z-plane — see the case section above.
- Reuse tolerances from the 2D tests as the starting point for the
  non-z-uniformity assertions (`CF_MEAN_TOL=0.018` etc. for BL,
  `PRE_RTOL=POST_RTOL=0.03` for OS, `min_order>=3.6/5.7` for EVC's KEEP4/KEEP6)
  and adjust only if the quasi-2D port's k=4-plane measurement actually
  requires it — it shouldn't, since k=4 receives identical treatment to the
  validated 2D case.
- Mark `@pytest.mark.integration` on all three (existing convention); add
  `@pytest.mark.slow` only if a measured build+run exceeds roughly the
  `test_sbli.py` threshold (~8 min) — EVC in particular runs 2 schemes × 3
  grid levels × build+run, so measure before deciding.

## Implementation order

Build and verify one case fully (source → standalone `calc.sh` run → Python
test) before starting the next, since they're independent but each carries
real Fortran-authoring risk:

1. **EVC first** — simplest (fully periodic, reuses generic `set_grid_cyclic`/
   `set_bc_cyclic` wholesale, no hand-written BC math, no viscous nuance) —
   shakes out the CMakeLists/calc.sh/Case()-workflow mechanics cheaply.
2. **BL second** — highest risk (hand-written wall/top/outlet BC math ported
   into 3D, the viscous z-fallback nuance, the two-phase cost-tuning pass).
3. **OS third** — same real-BC structural pattern as BL but simpler physics
   (Euler, no fallback nuance, already-cheap 2D timings) — should go quickly
   once BL's pattern is proven out.

## Verification

For each case, in order:

1. **Structural build**: `cd 3D_solver/<CASE> && cmake -B build && cmake
   --build build -j`. Fix any compile error before proceeding (expect the
   `nre1/nre2/rerank/start_rescale` boilerplate to be the most likely miss).
2. **Kernel-dispatch check**: grep the generated `build/calc_flux_base.f90` to
   confirm `BC_X/Y/Z` took effect as intended — BL/OS should show
   `calc_slau_x`/`calc_slau_y` (no `_in` suffix) and `calc_slau_z_in` (with
   `_in`); EVC should show `calc_keep_x_in`/`calc_keep_y_in`/`calc_keep_z_in`
   (all `_in`, matching what 2D EVC already exercises, now in 3D).
3. **Short correctness run** via `calc.sh` (`mpiexec -n 2`) at cheap
   settings (BL: the `nx=65,ny=25` correctness pass described above; OS/EVC:
   a handful of steps): confirm no NaN, VTK output opens, and the z-uniformity
   check passes at whatever tier applies — this is the cheap gate to clear
   before spending time on the full physics run.
4. **Full physics run + z-uniformity + analytical comparison**, first
   standalone via `calc.sh`/manual Python (to see actual numbers and tune any
   placeholder tolerance), then wired into the new pytest file.
5. **`pytest ouxsbli/tests/test_bl_3d.py ouxsbli/tests/test_evc_3d.py
   ouxsbli/tests/test_os_3d.py -v`** — all passing, with real (not placeholder)
   tolerances recorded with their measured values in comments, matching the
   documentation convention already established in `test_bl.py`.
6. Update `CLAUDE.md`'s 3D case table and "Python Test Suite" table to list
   the three new cases/tests.
