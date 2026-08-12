# Apply `fma()` and CUDA Fortran cache-hint intrinsics across the hot kernels

## Context

The user is manually hardening the GPU hot path for precision/performance by replacing
chained multiply-adds with the `fma()` intrinsic (single-rounding fused multiply-add) and
wants to extend the same treatment to CUDA Fortran cache-hint load/store intrinsics
(`__ldca`, `__ldcg`, `__ldcs`, `__ldlu`, `__ldcv`, `__stwb`, `__stcg`, `__stcs`, `__stwt`).

Git status shows this work already in progress and staged in
[3D_solver/src/calc_slau_3d.f90](3D_solver/src/calc_slau_3d.f90) and
[3D_solver/src/calc_visc_cent.f90.fypp](3D_solver/src/calc_visc_cent.f90.fypp) — both convert
symmetric-stencil weighted sums into pairwise sum/diff temporaries followed by a chain of
`fma()` calls, with the original expression kept as a comment above for traceability. No
cache-hint intrinsics exist anywhere in the repo yet — this will be a new pattern.

Repo-wide search confirms: no existing `fma`-wrapping fypp macro, no prior cache-hint usage,
and `load_smem_visc2.f90` / `load_smem_visc_cent.f90.fypp` / `load_smem_visc4.f90` /
`load_smem_visc_me4_base.f90` already stage shared-memory tiles via
`pipelineMemcpyAsync`/`pipelineWaitPrior` — those are a different (already-optimal) mechanism
and must not be touched.

Scope (per user's choice): 3D solver **and** 2D solver hot paths, for KEEP / SLAU / Hybrid
schemes. **Roe is explicitly excluded** (CLAUDE.md: "Roe scheme is not used... should [not]
be optimized").

Verification (per user's choice): build-only — `cmake --build` must succeed for representative
cases after each file's changes. No run/VTK-diff required.

## The established FMA pattern to follow

From the already-modified `calc_visc_cent.f90.fypp`, the house style is:
1. Keep the original expression as a `!`-commented line directly above.
2. Factor repeated/symmetric stencil points into named sum/diff temporaries first
   (independent of `fma`).
3. Chain `fma(coef, term, accumulator)` calls, innermost-first, ending in the final variable.
4. Wrap in a `block ... end block` when new temporaries are introduced, matching existing
   scoping style.

Example already committed (`calc_visc_cent.f90.fypp:6-12`):
```fortran
block
  real(8) s1, s2
  !ai = 0.0625d0 * (-a(1) + 9.d0 * (a(2) + a(3)) - a(4))
  s1 = a(2) + a(3)
  s2 = a(1) + a(4)
  ai = 0.0625d0 * fma(9.d0, s1, -s2)
end block
```
Apply this same shape everywhere below — don't invent a different convention.

## FMA targets

### 3D solver (`3D_solver/src/`)

- **`calc_steps.f90`** (NEW — highest call frequency in the codebase; runs every RK
  sub-stage, every interior point, every case regardless of scheme):
  - `calc_R` (lines ~27-37): the y- and z-direction accumulation lines
    `R(l) = R(l) + dtdzdx * (-F(l,...) + F(l,...))` and
    `R(l) = R(l) + dtdxdy * (-G(l,...) + G(l,...))` are literal `a*b+c` — convert to
    `R(l) = fma(dtdzdx, F(l,i,j+1,k)-F(l,i,j,k), R(l))` etc.
  - `calc_step` (line 103-104): `Q2(...) = Q(...) - coef1*R(l)` → `fma(-coef1, R(l), Q(...))`;
    `Rs(...) = Rs(...) + coef2*R(l)` → `fma(coef2, R(l), Rs(...))`.
  - `calc_step2_3` (line 143): `(coef1*Qin + coef2*Qout - R(l)) * coef4_inv` → nested
    `fma(coef1, Qin, fma(coef2, Qout, -R(l))) * coef4_inv`.
  - `calc_step4` (line 180): `Q(...) - R(l)*one_sixth` → `fma(-one_sixth, R(l), Q(...))`.
  - `calc_step1` (line 68) is a plain subtraction, no multiply — leave as-is.

- **`calc_slau_3d.f90`** (finish what's already started — same file, currently only the
  `pres` blocks were converted):
  - `SLAU1`/`HRSLAU2` `mass` line (`mass = 0.25d0 * (rho1*(un1+Vtp) + rho2*(un2-Vtm) -
    x*dp*over_c)`) — 3-product sum, nested-fma candidate.
  - `F2/F3/F4` lines (`mass1*u1 + mass2*u2 + pres*Norm(2)`, same for v/w) in both `SLAU1`
    and `HRSLAU2` — 3-term product sum, convert to nested fma.

- **`calc_keep_3d.f90`** (KEEP4/KEEP6 — active scheme for NSTGV, `ORDER=6`):
  - `KEEP4`: `F(1) = one_third*RV1 - (RV2+RV3)*one_24`; `pres` (4-term weighted sum);
    `F(2)/F(3)/F(4)` (4-product sums); `ene` internal+kinetic energy accumulation; `F(5)`
    pressure-diffusion term.
  - `KEEP6`: same shapes with 6 terms — `F(1)`, `pres`, `ruu/ruv/ruw` (6-product sums),
    `ene`, `F(5)` — largest per-call FMA payoff in the 3D convective path since KEEP6 runs
    at `ORDER=6`.

- **`calc_visc_me4_base.f90`** (active viscous stencil for NSTGV, `VISC_STENCIL='ME4Base'`):
  - `flux4` (`(-a1+26*a2-a3)*one_24`) — same shape as the already-converted `interp4`/`diff4`
    in `calc_visc_cent.f90.fypp`; factor `s=a1+a3`, `ans = one_24*fma(26.d0, a2, -s)`.
  - `calc_tau_straight`/`_LES`, `calc_tau_cross`/`_LES`: the `tmp1/tmp2/tmp3` stencil
    combinations and the final `t11`/`t12`/`ut11`/`vt12` `flux4`-shaped combinations — mirror
    the `calc_visc_cent.f90.fypp` treatment already done for the `N==4` branch.

- **`calc_visc_high.f90.fypp`** / **`calc_visc_high_internal.f90.fypp`** (6th-order viscous,
  used when `VISC_ORDER=6`): same `flux4`/`tau`-style weighted-sum patterns as
  `calc_visc_me4_base.f90` but wider stencils — apply the same nested-fma treatment.

- **`calc_visc2.f90.fypp`** (2nd-order viscous fallback): `muy`/`mvy`/`mux`/`mvx`
  (3-term product sums for face-averaged viscosity × velocity-difference), and
  `viscous_work` (2-product sum) — nested-fma candidates.

- **`calc_hybrid.f90`** (Ducros sensor, used when `SCHEME='Hybrid'`, e.g. TBL case):
  - `div = dudx + dvdy + dwdz` → nested fma of the three direction-derivative products.
  - `rot(1)**2+rot(2)**2+rot(3)**2` inside the `fd(...)` sensor formula → nested fma of
    three squares.

- **`src/calc_muscl.f90`** (`MUSCL4thnonTVD`, used by all `ORDER=6` reconstructions):
  `al = a2 + (-0.4d0*d1+2.2d0*d2+4.8d0*d3-0.6d0*d4)*one_twelfth` — 4-term weighted-diff sum.

### 2D solver (`2D_solver/src/`) — mirror the 3D treatment

- **`calc_steps.f90`** — same `calc_R`/`calc_step*` shapes as the 3D file above.
- **`calc_slau_2d.f90`** — port the exact two-`fma` `pres` refactor already done in
  `calc_slau_3d.f90` (same formula, not yet converted here), plus the `mass` and `F2/F3`
  (2D has no w-component) product sums.
- **`calc_keep_2d.f90`** — `KEEP4`/`KEEP6` `ruu/ruv`, `ene`, pressure-diffusion terms
  (2D drops the w/`RV*w` terms but is otherwise structurally identical to `calc_keep_3d.f90`).
- **`calc_visc_me4_base.f90`** (2D copy) — same `flux4`/`tau` treatment as the 3D file.
- **`calc_visc4.f90.fypp`** / **`calc_visc4_internal.f90.fypp`** — same weighted-sum shapes
  as 3D's `calc_visc_high*.f90.fypp`.
- **`calc_visc2.f90.fypp`** — same `muy`/`mvy`/`mux`/`mvx` treatment as 3D.
- **`calc_hybrid.f90`** (2D copy) — same Ducros `div`/`rot` treatment.

Skip entirely: `calc_roe_3d.f90`, `calc_roe_kernel(_internal).f90.fypp` (3D and 2D) — Roe is
unused per CLAUDE.md.

## Cache-hint intrinsic targets

CUDA Fortran exposes these as device-callable intrinsics taking a global-memory array
element reference: `value = __ldca(E(1,i,j,k))` for loads, `call __stcs(Q2(i,l,j,k), value)`
for stores. Apply only to genuine, un-staged global-memory traffic — do not touch anything
already staged through `pipelineMemcpyAsync` in the `load_smem_*` files.

- **`calc_steps.f90`** (highest impact — every RK stage, whole domain):
  - `calc_R`: `E(l,i,j,k)`/`E(l,i+1,j,k)`, `F(l,...)`, `G(l,...)` reads — un-staged, and each
    element is re-read by the neighboring thread on its other side → **`__ldca`**.
  - `dtdxdy(i,j)`, `dtdydz(j,k)`, `dtdzdx(i,k)` — 2D metric arrays broadcast across the third
    block dimension (shared by many threads) → **`__ldca`**.
  - `Q2(i+1,l,j+1,k+1) = ...` / `Qout(...) = ...` writes — write-once per launch, not re-read
    until the *next* kernel → **`__stcs`** (streaming store, avoid evicting cache needed by
    this kernel's own E/F/G reads).
  - `Rs(i,l,j,k)` read-modify-write accumulator in `calc_step`/`calc_step4` — genuinely reused
    across stages → leave as default (**`__stwb`**), do not mark streaming.

- **`calc_visc_high_internal.f90.fypp`** / 2D **`calc_visc4_internal.f90.fypp`**:
  - Un-staged wide-stencil reads of `mu`/`mut`/`T`/`qc2` (only `u,v,w` derivatives are staged
    via `load_smem_*`) — overlapping windows between adjacent threads → **`__ldca`**.
  - `E(l,i,j-1,k-1) = E(l,...) - txx` style read-modify-write of the flux array written
    earlier by the convective kernel, touched once per thread here → **`__ldcg`** for the
    read, **`__stcg`** for the write.

- **`calc_visc2.f90.fypp`** (3D and 2D): same two patterns as above at 2nd-order-stencil
  scale — un-staged `mu`/`mut`/`qc2`/cross-plane `Q` component reads → **`__ldca`**;
  `E/F/G(...) = ... - txx` read-modify-write → **`__ldcg`** read / **`__stcg`** write.

- **`calc_slau_kernel(_internal).f90.fypp`**, **`calc_hybrid_kernel(_internal).f90.fypp`**,
  **`calc_keep_kernel(_internal).f90.fypp`** (3D and 2D):
  - `sensor(i,j,k)`/`sensor(i+1,j,k)` (and y/z) reads, un-staged, single precision, reused by
    the neighboring thread → **`__ldca`**.
  - `E/F/G(:,i,j-1,k-1) = KEEP(...)`/`SLAU(...)`/`Hybrid(...)` flux-array stores — write-once,
    not re-read by this kernel → **`__stcs`**.

Do not modify: any array access inside `load_smem_visc2.f90`, `load_smem_visc_cent.f90.fypp`,
`load_smem_visc4.f90`, `load_smem_visc_me4_base.f90` (already optimal via async pipeline), or
`calc_flux_base.f90.fypp` (pure kernel-launch orchestration, no array math).

## Verification

Build-only, per the user's choice — after modifying each file (or logical group of files),
run `cmake --build build -j` for at least one case that actually exercises the changed code
path, using the compiler already present in this environment
(`/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/compilers/bin/nvfortran`, GPU visible via
`nvidia-smi`):

- 3D KEEP + ME4Base viscous (exercises `calc_steps.f90`, `calc_keep_3d.f90`,
  `calc_visc_me4_base.f90`, `calc_muscl.f90`): `cd 3D_solver/NSTGV && cmake -B build &&
  cmake --build build -j`.
- 3D SLAU + Hybrid + 6th-order viscous (exercises `calc_slau_3d.f90`, `calc_hybrid.f90`,
  `calc_visc_high*.f90.fypp`): build the SBLI or TBL case.
- 2D equivalents: build one 2D case exercising KEEP (e.g. `DSL`/`EVC`) and one exercising
  SLAU (`OS` or `SBLI`) to compile `2D_solver/src/*`.

A clean `nvfortran`/`fypp` build for each representative case is the pass criterion — no
functional run or VTK comparison is required for this pass. Flag any fypp template where a
generated-output check (`build/*.f90`) reveals the `fma`/cache-hint edit landed in the wrong
`#:if` branch for a given config combination.
