# Refactor 2D BL & SBLI cases + pytest validation

## Context

Two 2D cases were just added/reworked: `2D_solver/BL` (untracked; replaces the old two-stage BL0→`Qin.dat`→BL workflow with a self-contained **M=0.1** flat plate that grows its own boundary layer from a leading edge) and `2D_solver/SBLI` (**M=2.15**, β=30.8° oblique shock impinging on a laminar flat-plate BL at Xsh=80 mm, Re_Xsh≈1e5 — the Degrez benchmark). Goals: clean up both cases, and add pytest validation — BL vs incompressible Blasius (Cf and u-profile; user-confirmed), SBLI vs Degrez reference Cp/Cf with a **reduced config** (user-confirmed, ~20–40 min budget).

Blocking problems found in exploration:
1. **SBLI does not compile**: its `set.f90` uses SoA (`Q(nx,4,ny)`, single packed `QJ`) but `2D_solver/src/main.f90:32,62` and `calc_time_dev.f90.fypp:90,96,102,182-199` use AoS with four split `QJ_1..QJ_4(nx,ny)` arrays. BL's `set.f90` is the working reference pattern.
2. **BL leading edge is at x≈−20.3 mm** (wall starts at `i > 2*nx/12`, `x(1)=-0.04`), while `Cf.py` assumes LE at x=0 → Re_x wrong by 20 mm.
3. **Reference data is gitignored**: `2D_solver/SBLI/data_1/` matches `.gitignore:18` (`2D_solver/*/data*`) and the files also match `*.dat` (line 2).
4. **Dead code**: BL's `mod_globals.f90:46-61` oblique-shock block produces garbage at M=0.1 (neutralized by hand-overrides); unused imports/locals in both `set.f90`s; `2D_solver/src/set_bc_tbl_sbli.f90` is never compiled (absent from CMake); `set_init_common.f90`/`set_compressible_bl.f90` are pulled in via `CASE_EXTRA_SOURCES` but no longer called.
5. Plotting scripts (`BL/Cf.py`, `SBLI/data_1/Cf_all.py`, `Cp_all.py`) import external, non-repo `myvtk`; repo equivalent exists at `ouxsbli/tests/utils/vtk_reader.py` (same API).
6. Committed run settings are prohibitive for tests (BL: 3.3M steps with dt ~50× below CFL; SBLI: 5M steps on 564×513). The `Case` harness (`ouxsbli/case.py`) patches `dt/endT/np/nx/ny` per test.

## Stage 1 — Make SBLI compile (AoS conversion)

**`2D_solver/SBLI/mod_globals.f90`** — add named geometry constants (plain `parameter ::` so the pytest patcher can rewrite them):
```fortran
! flat-plate geometry
real(8), parameter :: x_in = -16.d0 * blt   ! x(1); shock geometry depends on this
real(8), parameter :: Xsh  = 80.d0 * blt    ! inviscid shock impingement point
integer, parameter :: i_LE = nint(-x_in * dble(nx-1) / Lx) + 1  ! first no-slip point; =52 at nx=564 (matches old nx/11+1); LE stays ~x=0 for any nx
```

**`2D_solver/SBLI/set.f90`** — rewrite to BL's working convention, preserving physics exactly:
- Header: `use mod_globals, only : gamma, rho0, u0, p0, rho2, p2, ux, uy, beta, x_in, Xsh, i_LE`; `use mod_constant, only : gamma_1, over_gamma_1`. Drop `use set_bc_common` / `use set_init_common` and unused symbols.
- `set_grid`: `x(1) = x_in`; delete commented-out clustering block (:35-44), unused `ny_b`, `dy1`.
- `set_init`: `Q(nx,4,ny)` → `Q(nx,ny,4)` (indices `Q(i,l,j)` → `Q(i,j,l)`); shock-injection condition `(ys(j)/0.08d0) > 1.2d0*dtan(beta)` → `ys(j) > (Xsh - x_in) * dtan(beta)` (mathematically identical: 1.2 = (0.08+0.016)/0.08); wall stamp loop `do i = nx/11+1, nx` → `do i = i_LE, nx`. Comment that the inlet column i=1 is never updated (interior-only kernels, no inlet BC) — it is the frozen shock generator.
- `set_bc`: signature → `set_bc(myrank, nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4)` with four `(nx,ny)` device arrays (match `calc_time_dev.f90.fypp` call sites). Outlet: flatten `do j/do l` nest to BL's 4-assignment form. Top Riemann loop: mechanical `QJ(i,l,·)` → `QJ_l(i,·)`; **preserve exactly** the entropy switch (`vb>=0`: interior entropy/`ub=uin`; else `p2/rho2**gamma`, `ub=ux`) and `Rm = uy - 2*c_ext*over_gamma_1`; replace the full-width parens in the :105 comment. Bottom: `if (i <= nx/11)` → `if (i < i_LE)` (identical set). Delete unused `l, ireq, ierr, istat`. Add trailing newline.

**`2D_solver/SBLI/CMakeLists.txt`**: `CASE_EXTRA_SOURCES` → `""` (drop `set_compressible_bl.f90` + `set_init_common.f90`, both uncalled), matching OS/ST/EVC/DSL layout.

**Verify**: `cmake -B build && cmake --build build -j` (first-ever compile of this case). Short smoke run (`mkdir -p data && mpirun -n 2 ./build/a.out`, kill after Q00001.vtr): fields finite, ρ>0, inlet column above y=57.2 mm holds the post-shock state, shock trace forming.

## Stage 2 — BL cleanup + leading edge at x=0

**`2D_solver/BL/mod_globals.f90`**:
- Delete the dead oblique-shock block (lines 46-61: `beta, Ms, Ms2, theta, T2, p2, rho2, u1, v1, a1, u2, v2, u_magnitude, ux, uy`). Keep `rf`/`Taw`. Keep `nx, ny, M0, dt, endT, np, nt` as plain parameter lines (`test_patcher.py` patches `nx`/`M0` in this file).
- Add: `integer, parameter :: i_LE = 2 * nx / 12 + 1  ! first no-slip wall point; set_grid puts the LE at x = 0`.
- Drop stale trailing comments (`!176!20`, `!2.15d0 !2`, …).

**`2D_solver/BL/set.f90`**:
- Header: `use mod_globals, only : gamma, rho0, u0, p0, i_LE`; `use mod_constant, only : gamma_1, over_gamma, over_gamma_1`. Drop `use set_bc_common`/`use set_init_common` and unused symbols.
- `set_grid`: `x(1) = -0.04d0` → `x(1) = -dble(i_LE - 1) * dx1` (LE lands exactly at x=0; whole-grid shift only — flow field is index-based and bit-identical). Name the stretch factor `real(8), parameter :: s = 1.6d0`. Delete unused `ny_b`, `dy1`.
- `set_init`: `do i = 2*nx/12 + 1, nx` → `do i = i_LE, nx`.
- `set_bc`: `if (i <= 2*nx/12)` → `if (i < i_LE)`; delete dead `T = Taw - rf*u0**2/(2*Cp)` (:70), unused `ub, l, ireq, ierr, istat`. Comment the intentional frozen-inlet (no i=1 BC) behavior.

**`2D_solver/BL/CMakeLists.txt`**: `CASE_EXTRA_SOURCES` → `""`.

**Verify**: `mv data data_prev`, clean rebuild, run 1–2 outputs; new `Q00000.vtr` fields bitwise-equal to `data_prev/Q00000.vtr`, x-coords shifted by +0.0203125 exactly, y identical. `pytest ouxsbli/tests/test_patcher.py` still passes.

## Stage 3 — Dead shared files + calc.sh

- `git rm 2D_solver/src/set_bc_tbl_sbli.f90 2D_solver/src/set_init_common.f90` — verified: only referenced by the `use` lines and `CASE_EXTRA_SOURCES` entries removed in Stages 1–2. Keep BCs **case-local** (the two cases' top BCs are genuinely different physics; the dead module matches neither; repo convention is BCs in each case's `set.f90`). Keep repo-root `src/set_compressible_bl.f90` (3D solver uses it).
- Both `calc.sh` (byte-identical): move the `cp set.f90 mod_globals.f90 data/` provenance snapshot **before** the backgrounded `mpirun` and add `config.fypp` to it; `mkdir -p data`.

**Verify**: clean rebuild of BL, SBLI, and one untouched case (`2D_solver/OS`).

## Stage 4 — Reference data, .gitignore, plotting scripts

- `.gitignore`: append `!2D_solver/*/ref/**` (needed: global `*.dat` at line 2 catches ref files; `ref/` itself matches no dir rule, so negation works).
- `mv 2D_solver/SBLI/data_1/ → 2D_solver/SBLI/ref/` (files currently untracked): 6 `data_*_C{p,f}.dat` files (rename `Vila-Pérez` → ASCII `Vila-Perez`) + `Cf_all.py`, `Cp_all.py`. `git add`, confirm `git check-ignore` returns nothing.
- Refactor the three plotting scripts (kept as user-facing tools): replace `from myvtk import …` with `sys.path.insert(0, repo_root)` + imports from `ouxsbli.tests.utils.vtk_reader` and the new `ouxsbli.analysis` modules (Stage 5); inline the tiny `myParams` rcParams block; default snapshot = `latest_vtr(<case>/data)` with optional argv override (fixes stale `./Q03000.vtr`/`./Q00800.vtr`); ref-data paths relative to `Path(__file__).parent`. **`BL/Cf.py` only**: delete dead Eckert block, drop the meaningless `Xsh_mm=80` normalizer (plot vs x in mm), fix `X_MAX_MM` to the actual domain.

## Stage 5 — Shared Python analysis utils (single source of truth for tests + scripts)

**`ouxsbli/analysis/wall.py`** (new; numpy-only, lazy VTK import; lifted from `BL/Cf.py` / `Cf_all.py` / `Cp_all.py`):
- `sutherland_mu(T)` (constants mu_ref=1.716e-5, T_ref=273.2, S=111.0 — match `src/mod_constant.f90.fypp:32-33`)
- `one_sided_deriv_3pt(...)`, `wall_dudy(y, u2d)` (u=0 at wall + u[1],u[2], non-uniform 3-pt)
- `compute_cf_cp(y, rho2d, u2d, p2d, R, p_inf, q_inf) -> (cf, cp)` — T_w=p_w/(ρ_w R), μ_w Sutherland, τ_w=μ_w·du/dy|w, cf=τ_w/q∞, cp=(p_w−p∞)/q∞
- `wall_coeffs_from_vtr(path, R, p_inf, q_inf) -> (x, cf, cp)`
- `find_zero_crossings(x, f)`, `load_reference(path, x_scale=1.0)` (loadtxt, sort by col 0, scale x)

**`ouxsbli/analysis/blasius.py`** (new): 29-point `ETA_TAB`/`FP_TAB` transcribed from `src/set_compressible_bl.f90:6-16` (solver's own ground truth; interp error ≤~2.5e-3); `fprime(eta)`, `eta(y, x, u_inf, nu_inf)`, `cf_blasius(x, rho_inf, u_inf, mu_inf) = 0.664/sqrt(Re_x)`. No scipy.

## Stage 6 — `ouxsbli/tests/test_bl.py`

Module-scoped fixture builds+runs once; `pytestmark = pytest.mark.integration`.
```python
Case(source="2D_solver/BL", workdir="./tmp/bl",
     dt=5e-8,           # CFL_y ≈ 0.18 (committed 3e-9 was ≈0.011)
     endT=1e-2,         # ≈2.8 flow-throughs (Lx/u0 = 3.5 ms)
     np=10,             # nt auto-derives to 20000 (do NOT pass nt)
     output_precision=8)
case.run(nranks=2)      # 2D print path requires the even/odd rank pair
```
- Grid preconditions: `y[0]==0`, `min(|x|) < 1e-9` (LE on-grid at 0), `x[-1] >= 0.079`.
- **Cf vs Blasius**: window x ∈ [32, 72] mm (Re_x 1.8e4–4.1e4; avoids LE singularity and outlet buffer). Assert mean|cf/cf_Blasius − 1| < 0.06, max < 0.12, cf > 0 and monotonically decreasing.
- **u-profile**: column nearest x=56 mm (Re_x≈3.2e4); η = y·√(u0/(ν0·x)); ≥8 points in 0<η≤6; max|u/u0 − f′(η)| ≤ 0.04, RMS ≤ 0.02; max|v|/u0 < 0.02.
- **Steadiness guard**: window-mean Cf change between last two snapshots < 1% (run before comparisons; fail with the numbers if unconverged).
- Est. runtime: ~200k steps ≈ 2 min run (measured ~1.9–2k steps/s on the RTX 4060 Laptop) + build ≈ **3–5 min total**.

## Stage 7 — `ouxsbli/tests/test_sbli.py`

`pytestmark = [pytest.mark.integration, pytest.mark.slow]`; module-scoped fixture.
```python
Case(source="2D_solver/SBLI", workdir="./tmp/sbli",
     nx=276, ny=257,    # dx=0.64 mm exactly → LE exactly at x=0 (i_LE=26); y(2)≈82 µm
     dt=3e-8,           # CFL ≈ 0.22 combined; viscous limit ≫
     endT=6e-3,         # ≈25 flow-throughs
     np=10,             # nt auto-derives to 20000
     output_precision=8)
```
Reference: **Moro's data** (user decision — newer than Degrez's): `REPO_ROOT/2D_solver/SBLI/ref/data_Moro_{Cp,Cf}.dat`, loaded with `x_scale=1.25` (digitizer x/L → X/Xsh; the scripts apply it to all datasets). Derived at runtime from ref data (not hardcoded); for orientation: separation X/Xsh≈0.72, reattachment ≈1.25, Cf_min≈−8.4e-4 near X/Xsh≈1.04, plateau Cp≈0.164–0.168. (Degrez/Vila-Perez files stay in `ref/` for the plotting scripts.)
- **Cp vs Moro**: interp sim Cp onto ref points in X/Xsh ∈ [0.4, 1.75]; normalize by plateau (max ref Cp = 0.168) since near-zero Cp makes relative error ill-conditioned: mean err < 0.10, max < 0.20.
- **Cp plateau**: mean over X/Xsh ∈ [1.4, 1.8] vs ref mean, `assert_close_relative(rtol=0.12)`.
- **Cf separation bubble** (window X/Xsh ∈ [0.3, 1.5]): min(Cf) < −1e-4; zero-crossing locations within ±0.15 X/Xsh of the ref crossings (via shared `find_zero_crossings`).
- **Steadiness guard**: max|ΔCp|/0.174 < 0.02 and crossings move < 0.05 X/Xsh between last two snapshots.
- Est. runtime: 200k steps on 71k points ≈ 8–30 min (marker `slow` isolates it).
- Fallback ladder if the bubble isn't captured on the reduced grid: (i) ny→321; (ii) drop crossing asserts, keep Cp quantitative + concrete qualitative Cf check (local min in [0.6,1.3] with <25% of local Blasius Cf, within ±0.15 of the Moro ref min at ≈1.04).

## Stage 8 — Harness/config/doc updates

- `ouxsbli/case.py:120`: extend copytree `ignore_patterns` with `"data", "recal", "nohup.out"` — otherwise BL's 26 stale `Q*.vtr` get copied into the workdir and `latest_vtr()` silently reads the old run.
- `ouxsbli/tests/conftest.py`: register `slow` marker (`deselect with -m "not slow"`).
- `pyproject.toml`: add `scipy` to dependencies (pre-existing undeclared dep of `utils/oblique_shock.py`).
- `CLAUDE.md`: 2D case table — BL is now "Laminar flat-plate boundary layer (M=0.1, validated vs Blasius)"; add `test_bl.py`/`test_sbli.py` rows to the test table.

## Verification (end-to-end)

1. After each Fortran stage: clean rebuild of the touched case; BL Stage-2 bitwise-identity check vs `data_prev`.
2. `python -m py_compile` on all three refactored plotting scripts; run each against real data (no `myvtk`).
3. `pytest ouxsbli/tests/test_patcher.py` (no GPU) after mod_globals edits.
4. `pytest ouxsbli/tests/test_bl.py -v` (~5 min) — iterate tolerances only if failures are marginal and systematic (document virtual-origin fit fallback if needed).
5. `pytest ouxsbli/tests/test_sbli.py -v -m slow` (~10–30 min); apply the fallback ladder if the bubble is too shallow.
6. `pytest ouxsbli/tests/ -m "not slow"` to confirm nothing else broke; `git status` shows `ref/` tracked, no `data/`/`build/` leakage.

## Commit order

(1) SBLI AoS fix + constants; (2) BL cleanup + LE shift (combines with the already-staged BL/BL0 deletions); (3) dead shared files + CMake + calc.sh; (4) ref-data move + .gitignore + script refactor; (5) analysis utils + tests + conftest/pyproject/case.py + CLAUDE.md.
