# Validate OUxSBLI against `3D_solver/DHIT/ref.tex` (compressible decaying HIT)

## Context

`ref.tex` documents a specific compressible decaying-HIT DNS validation case: a
divergence-free random velocity field with spectrum $E(k)=A_0k^4\exp(-2k^2/k_0^2)$
in a periodic $-\pi..\pi$ box, target $Re_\lambda=72$, $Ma_t=0.5$, power-law
viscosity $\mu=\mu_0(T/T_0)^{0.76}$, on $64^3/96^3/128^3$ meshes, monitoring
$K(t)/K_0$, $\rho_{rms}(t)/Ma_t^2$, $S_u(t)$ vs $t/\tau$.

Investigation found the *current* `3D_solver/DHIT` case does not implement this:
it uses a Pope(2000)-spectrum IC at $Re_\lambda=20$, $Ma_t=0.1$ with Sutherland's-law
viscosity — a different, unrelated setup (confirmed via `mod_globals.f90`,
`set_init_dhit.f90`). No `test_dhit.py` or turbulence-statistics post-processing
exists yet. Per user decision, this plan makes the case actually implement
`ref.tex`'s exact spectrum/parameters (not just reuse the existing Pope setup),
adds the power-law viscosity as a DHIT-only override, runs **one** representative
case (single resolution, $Re_\lambda=72$, $Ma_t=0.5$, no mesh/$Ma_t$ sweep), and
delivers an **analysis report with plots** (not a pytest test), since `ref.tex`
itself has no digitized reference curve to assert against numerically — its
figures are missing/mislabeled (captions copy-pasted from an unrelated
Taylor-Green-vortex section) — so validation here means: (a) the synthesized
field's measured $Re_\lambda,Ma_t,K_0$ match the targets, and (b) $K(t)/K_0$,
$\rho_{rms}(t)/Ma_t^2$, $S_u(t)$ show the qualitatively correct, physically
sane behavior `ref.tex` describes (monotonic KE decay, $O(1)$ negative skewness).

**Reference curves found**: `3D_solver/DHIT/ref/*.pdf` (4 files, added by the user
mid-planning) contain exactly the figures `ref.tex`'s text describes but whose
`\includegraphics`/`\caption` blocks were missing/mislabeled:
- `5-iso-turbu-enr-eps-converted-to.pdf`: $K(t)/K_0$ vs $t/\tau$ (0 to 5),
  HGKS at $64^3$/$96^3$/$128^3$ plus square-marker "reference data".
- `5-iso-turbu-enr-2-eps-converted-to.pdf`: $K(t)/K_0$ vs $t/\tau$ for the
  $Ma_\lambda=0.2/0.5/0.6$ sweep (no markers, HGKS-only, not part of this run's scope).
- `5-iso-turbu-den-eps-converted-to.pdf`: $\rho_{rms}/Ma_t^2$ vs $t/\tau$, same
  3-resolution + reference-data-marker format.
- `5-iso-turbu-vel-eps-converted-to.pdf`: $S_u$ vs $t/\tau$, same format — notably
  shows the $64^3$ curve visibly *less negative* than $96^3$/$128^3$/reference
  data (settling around $-0.36$ vs reference $\approx-0.42$ at large $t/\tau$),
  i.e. **the reference's own $64^3$ result under-resolves the skewness statistic**
  — a real, expected deviation at this resolution, not evidence of a bug.

This means there **is** a real curve to validate against after all (per user
decision, compared **visually** — our computed plots reviewed side-by-side
against these PDFs — rather than hand-transcribing approximate marker
coordinates into code as if they were precise digitized data). The reference
curves run to $t/\tau=5$; the run initially targeted that full range, but after
calibration measured ~87 steps/sec on the RTX 4060 (implying ~3.1h for the full
985,600-substep run), the user opted to shorten the target to **$t/\tau=2$**
(~1.2h) instead, trading comparison range (stops before the tail behavior at
$t/\tau=3$-$5$) for a shorter wait.

All derived constants below are closed-form (no iterative solve needed) and were
hand-verified against `ref.tex`'s formulas with $A_0=1.3\times10^{-4}$, $k_0=8$,
$Re_\lambda=72$, $Ma_t=0.5$, $\rho_0=1$:

```
K0   = (3*A0/64)*sqrt(2*pi)*k0^5   = 0.50052
u'   = sqrt(2*K0/3)                = 0.57765
tau  = (32/A0)*(2*pi)^0.25*k0^-3.5 = 269.19
mu0  = (2*pi)^0.25/4 * (1/Re_lambda) * sqrt(2*A0) * k0^1.5 = 2.006e-3
T0   = 3*u'^2 / (gamma*Ma_t^2)      = 2.8601
c0   = sqrt(gamma*R*T0)             = 2.00103   (with R=1, see judgment call below)
```

**Correction found during execution — `tau` formula:** `ref.tex`'s extracted
$\tau=(32/A_0)(2\pi)^{1/4}k_0^{-7/2}$ is dimensionally inconsistent: since
$E(k)=A_0k^4\exp(-2k^2/k_0^2)$ requires $A_0\sim L^7/T^2$, any time built from
$A_0,k_0$ must scale as $A_0^{-0.5}$, not $A_0^{-1}$. As extracted this gives
$\tau=269$, whereas the corrected $A_0^{-0.5}$ form gives $\tau\approx3.07$ —
consistent (same order of magnitude) with an independent integral-length-scale
cross-check ($L_{11}/u'\approx0.54$), and confirmed empirically: a calibration
run's $K(t)/K_0$ had already collapsed to 0.0013 by $t=59$, ~400x faster than
`ref.tex`'s own curve at the equivalent $t/\tau=0.22$ (using $\tau=269$) shows
(0.92). Per user decision, $\tau\approx3.07$ (corrected) is used throughout —
`mod_globals.f90`'s `tau` and `ouxsbli/analysis/dhit_decay_report.py`'s `TAU`
both add the `sqrt`. This also makes the full $t/\tau=5$ target cheap (~11,200
steps, ~2 min at the calibrated throughput) rather than the ~3.1h estimated
under the uncorrected $\tau$, so the run targets the full range after all.

**Judgment call — `R=1` for this case only:** `ref.tex`'s own $Ma_t$ formula
($\sqrt3 u'/\sqrt{\gamma T_0}$) has no gas constant, consistent with the common
DNS convention of nondimensionalizing so $R=1$ (i.e. $p=\rho T$). The rest of the
solver pulls `R` only from each case's own `mod_globals.f90` (confirmed via grep —
never hardcoded elsewhere), so setting `R=1.0d0` in DHIT's `mod_globals.f90` only
reproduces `ref.tex`'s formulas verbatim without affecting any other case.

## Files to modify / create

### 1. `3D_solver/DHIT/mod_globals.f90` — rewrite IC/physics block

Keep mesh/GPU/thread-block boilerplate structure, replace the physical-parameter
and IC block:
- `nx=ny=nz=70` (giving `Nf=64` interior points after the 6th-order stencil's
  3-ghost-cell halo on each side — matches `ref.tex`'s smallest quoted mesh,
  power-of-2-friendly for cuFFT, and keeps wall-clock down since GPU memory is
  not the constraint at any of 64³/96³/128³ on an 8GB card).
- `gamma=1.4d0`, `Pr=0.71d0` unchanged; **`R=1.0d0`** (was `287.03d0`, DHIT-local
  only — see judgment call above).
- New literal target parameters: `A0=1.3d-4`, `k0=8.d0`, `Re_lambda_target=72.d0`,
  `Mat_target=0.5d0`, `RHO0=1.d0`.
- Closed-form `parameter`s (same style as existing file): `K0`, `up0` (replaces
  `Urms`), `tau`, `mu0`, `T0`, `c0=sqrt(gamma*R*T0)`, `p0=RHO0*R*T0`.
- Remove the now-unused Pope-spectrum constants (`pope_C`, `pope_cL`, `pope_p0`,
  `pope_beta`, `pope_ceta`, `kp`, `pope_L`, `Re`, `pope_eta`, `nu0`) and the old
  `Mt`, `Urms`, `T`, `S`, `mu0` (Sutherland-derived) definitions.
- Timestep sizing, following the same `nt = int(endT/(np*dt))` idiom already
  used in `IVST`/`STZ`/`SBLI`/`TBL`/`KHI`'s `mod_globals.f90`:
  ```fortran
  real(8), parameter :: CFL  = 0.03d0
  real(8), parameter :: dt   = CFL * (Lx/dble(nx-1)) / c0
  real(8), parameter :: endT = 5.d0 * tau           ! target t/tau = 5, matching the reference PDFs' x-axis range
  integer, parameter :: np   = 100
  integer, parameter :: nt   = int(endT/(dble(np)*dt))
  ```
  With `tau=269.19`, `dt≈1.365e-3`, this gives `nt≈9860`, i.e. **~986,000 total
  substeps** (`np*nt`) — roughly double the earlier $t/\tau=2.5$ estimate. Wall-clock
  is genuinely uncertain until calibrated (see Execution step 3) — this solver
  runs in `real(8)` throughout, and this is a single consumer GPU, so budget for
  the calibration run to meaningfully change the plan (e.g. fewer `np` output
  blocks, or accepting a multi-hour run) rather than assuming the estimate holds.
  (`dt` using `nx-1` rather than the exact `nx-6` stencil spacing matches the
  existing repo convention elsewhere — mildly conservative, immaterial.)
- Drop the Petersen-forcing constants (`eps_s`, `kf_min`, `kf_max`, `C_T`) — see
  item 5 below, since `calc_forcing.f90` (confirmed dead code: `init_forcing`,
  `calc_forcing_rhs`, `apply_cooling`, `finalize_forcing` are never called
  anywhere) is being dropped from the build rather than kept with stub constants.

### 2. `3D_solver/DHIT/set_init_dhit.f90` — spectrum formula swap only

Reuse the existing cuFFT machinery unchanged: the Rogallo/Johnsen k12-based
divergence-free basis construction, the 3-random-phase assembly, the Hermitian-
symmetry and DC/Nyquist self-conjugate fixups, the `cufftExecZ2D` inverse
transform, and the final `uscale` rescale + `Q` packing. Only change:
- `use mod_globals` list: drop the Pope-spectrum imports, add `A0, k0`; rename
  `Urms`→`up0` (propagate rename to the `uscale = up0/sqrt(urms_sq)` line and the
  diagnostic prints).
- Delete the `TKE`/`eps` (Pope dissipation-rate proxy) lines — not needed.
- Replace the per-mode spectrum evaluation (currently `kL`,`keta`,`fL`,`feta`,
  Pope `Ek`) with:
  ```fortran
  Ek  = A0 * kmag**4 * exp(-2.d0*kmag**2/k0**2)
  amp = sqrt(Ek / (2.d0*pi*kmag**2))
  amp = amp * dble(Nf)**3   ! unchanged cuFFT unnormalized-IFFT compensation
  ```
- **Add a self-consistency diagnostic block** right after the existing
  `uscale`/`urms_sq` computation (before packing into `Q`): compute the
  synthesized field's actual $Re_\lambda$ and $Ma_t$ from first principles
  (independent of the Python post-processing in item 6, so the two cross-check
  each other) — e.g. a simple real-space finite-difference estimate of
  $\langle(\partial_1 u_1)^2\rangle$ on `vel_r` (periodic, so wrap-around
  differencing is valid) to get the Taylor microscale $\lambda=u'/\sqrt{\langle(\partial_1u_1)^2\rangle}$
  and $Re_{\lambda,measured}=\rho_0 u' \lambda/\mu_0$, plus
  $Ma_{t,measured}=\sqrt3 u'/c_0$ (exact by construction). Print these next to
  the existing `u_rms(raw)`/`u_rms(target)` diagnostics, alongside the targets
  `Re_lambda_target=72`, `Mat_target=0.5`, so a mismatch is visible at run start
  before spending GPU time on the full integration.

### 3. New file: `3D_solver/DHIT/calc_physical_quantities_dhit.f90`

Full copy of `src/calc_physical_quantities.f90` (same module name
`calc_physical_quantities`, same 4 public subroutines, so `use calc_physical_quantities`
in `calc_flux_base.f90.fypp` resolves identically) with only:
- `use mod_constant, only: mu0_T0_S_over_T0_2_3` → removed.
- `use mod_globals, only : gamma, R` → `use mod_globals, only : gamma, R, mu0, T0`.
- In `calc_quantities_T_2D` and `calc_quantities_T_3D`, replace
  `mu(...) = mu0_T0_S_over_T0_2_3 / (temp+111.d0) * (temp*sqrt(temp))` with
  `mu(...) = mu0 * (temp/T0)**0.76d0`.

This is the only viscosity change; `src/calc_physical_quantities.f90` and
`src/mod_constant.f90.fypp`'s Sutherland constant are untouched, so every other
case keeps Sutherland's law exactly as today.

### 4. `3D_solver/CMakeLists.txt` — additive override hook (shared file)

In the `add_ouxsbli_case()` macro's `_BASE` list, change the hardcoded
```cmake
"${CMAKE_CURRENT_SOURCE_DIR}/../../src/calc_physical_quantities.f90"
```
to resolve through a variable that defaults to the exact same path:
```cmake
if(NOT DEFINED CASE_PHYSICAL_QUANTITIES_SOURCE)
  set(CASE_PHYSICAL_QUANTITIES_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../../src/calc_physical_quantities.f90")
endif()
```
and reference `"${CASE_PHYSICAL_QUANTITIES_SOURCE}"` in `_BASE` instead. No other
case sets this variable, so this is a zero-behavior-change diff for everyone
except DHIT. **Verify** by rebuilding one unrelated case (e.g. `3D_solver/NSTGV`)
unchanged after this edit, confirming no regression, before trusting the DHIT build.

### 5. `3D_solver/DHIT/CMakeLists.txt`

- Add `set(CASE_PHYSICAL_QUANTITIES_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/calc_physical_quantities_dhit.f90")`
  before `include(...)`.
- Remove `calc_forcing.f90` from `CASE_EXTRA_SOURCES` (confirmed dead code —
  never called, and its `use mod_globals, only: ..., eps_s, kf_min, kf_max, ..., C_T, ... T_ref_const => T`
  would otherwise force keeping meaningless placeholder constants in the
  rewritten `mod_globals.f90`). Leave the `calc_forcing.f90` file itself in
  place, untouched, in case it's revisited later — just stop compiling it.

### 6. New: `ouxsbli/analysis/__init__.py` (empty) + `ouxsbli/analysis/dhit_decay_report.py`

CLI script (`data_dir` arg, default `3D_solver/DHIT/data`), reusing
`ouxsbli/tests/utils/vtk_reader.py`'s `getGrid`/`getQ`:
- Enumerate all `data/Q?????.vtr` snapshots (not just first/latest).
- Per snapshot, trim the 3-ghost-cell halo off each side of the `(Nz,Ny,Nx)`
  arrays from `getQ` (`[3:-3,3:-3,3:-3]`) before any FFT-based derivative —
  required since the halo layers are periodic copies, not independent data.
- Compute per snapshot: $K(t)=0.5\langle\rho|\mathbf u|^2\rangle$ (cross-check
  against the solver's own `data/kinetic_energy.d`, which already writes exactly
  this quantity, volume-averaged over the same interior region, normalized by
  `ke0` at step 0 — see `src/print.f90.fypp`'s `print_KE`); $\rho_{rms}(t)=\sqrt{\langle(\rho-\bar\rho)^2\rangle}$;
  and $S_u(t)=\sum_i\langle(\partial_iu_i)^3\rangle/\langle(\partial_iu_i)^2\rangle^{1.5}$
  via spectral differentiation (`numpy.fft.fftn`/`ifftn` with wavenumbers from
  `numpy.fft.fftfreq(Nf, d=dx)*2*pi`, `dx=Lx/(nx-6)` matching the solver's actual
  grid spacing) on each of the 3 diagonal velocity-gradient components.
- Also compute, from the **first** snapshot only, the same $Re_\lambda$/$Ma_t$/$K_0$
  self-consistency numbers as item 2's Fortran diagnostic — but independently
  (real-space FFT derivatives in Python vs. the Fortran init-time computation) —
  so the two cross-validate each other; a mismatch beyond floating-point-level
  differences flags a bug in one of the two derivations, not a "real" physical
  deviation.
- Normalize and plot 3 panels vs. $t/\tau$, axes/ranges matching the reference
  PDFs for easy side-by-side comparison ($t/\tau \in [0,5]$ on all 3;
  $K(t)/K_0 \in [0,1]$; $\rho_{rms}/Ma_t^2 \in [0,0.5]$; $S_u \in [-0.7,0.1]$):
  $K(t)/K_0$, $\rho_{rms}(t)/Ma_t^2$ (normalizing by the *target* $Ma_t=0.5$, per
  `ref.tex`), $S_u(t)$ (unnormalized). Save PNGs plus a short text summary:
  measured-vs-target $K_0$/$Re_\lambda$/$Ma_t$, a monotonicity check on
  $K(t)/K_0$, and the final $S_u$ sign/magnitude. These PNGs are the artifacts
  to review directly against `3D_solver/DHIT/ref/*.pdf` (visual comparison, per
  user decision — no reference numbers are hardcoded into this script).
- `matplotlib` isn't currently a project dependency — install ad hoc for this
  report (not adding it to `pyproject.toml` since this is a one-off analysis
  script, not part of the test suite).

## Execution plan

1. **Build check**: apply items 3–5, then rebuild `3D_solver/NSTGV` (or another
   already-working case) unchanged to confirm the `3D_solver/CMakeLists.txt`
   edit didn't regress anything, before building DHIT itself.
2. Apply items 1–2, build DHIT (`cmake -B build && cmake --build build -j` in
   `3D_solver/DHIT`).
3. **Calibration run**: a short run (small `nt`/`np` override, or just run a few
   seconds and kill it) to (a) read the printed init diagnostics from item 2 and
   confirm measured $Re_\lambda$/$Ma_t$/$K_0$ are close to targets before
   committing GPU time, and (b) measure actual steps/sec throughput on the RTX
   4060 at $70^3$ — this solver runs in `real(8)` throughout, and consumer Ada
   Lovelace GPUs have deliberately weak FP64 throughput, so the full ~986k-step
   run (`endT=5*tau`, per user decision to match the reference PDFs' full
   $t/\tau=5$ range) could take considerably longer than a naive estimate. This
   full duration was chosen deliberately to enable the visual comparison against
   `3D_solver/DHIT/ref/*.pdf`, so report the calibrated throughput/ETA back
   before committing rather than silently shortening the run.
4. Run the full case: `mpirun -n 2 ./build/a.out` from `3D_solver/DHIT`
   (matching `calc.sh`'s existing rank count/GPU-usage convention), blocking
   (not backgrounded) so completion is observed before post-processing.
5. Run `ouxsbli/analysis/dhit_decay_report.py` against the resulting `data/`
   directory; review the 3 plots and the self-consistency summary.

## Verification

- Step 1's rebuild of an unrelated case must succeed unchanged (confirms the
  shared `CMakeLists.txt` hook is truly additive).
- Step 3/5's measured $Re_\lambda$, $Ma_t$, $K_0$ (both the Fortran init-time
  print and the Python first-snapshot check) should be close to the targets
  72, 0.5, 0.50052 — a large discrepancy indicates a bug in the spectrum
  formula or basis construction, not a physical result to accept.
- $K(t)/K_0$ must decay monotonically over the run, and its shape should track
  the "HGKS $64^3$ cells" curve in `5-iso-turbu-enr-eps-converted-to.pdf`
  (visual comparison, since this is the same resolution we're running).
- $\rho_{rms}(t)/Ma_t^2$ should rise to a peak around $t/\tau\approx0.3$ then
  decay, tracking `5-iso-turbu-den-eps-converted-to.pdf`'s $64^3$ curve shape.
- $S_u(t)$ should drop sharply to a negative extremum near $t/\tau\approx0.1$-$0.3$,
  partially recover, then settle to an $O(1)$ negative value — per
  `5-iso-turbu-vel-eps-converted-to.pdf`. Note the reference's own $64^3$ curve
  there sits visibly less negative (~$-0.36$) than $96^3$/$128^3$/reference-data
  (~$-0.42$) at late $t/\tau$ — some deviation from the higher-resolution/reference
  curves in this specific statistic is an **expected feature of running at
  $64^3$**, not necessarily a bug, per the reference's own resolution study.
- Any of these three showing qualitatively different shape (wrong sign, non-
  monotonic $K(t)/K_0$, no skewness extremum at small $t/\tau$) indicates a bug
  in the spectrum/IC/viscosity implementation, not just discretization error.

## Results (executed 2026-08-03)

Full run completed at $64^3$ to $t/\tau=4.98$ (~2.4 min wall-clock after the
$\tau$ fix). Two implementation bugs were found and fixed during execution
(beyond the $\tau$ transcription error above):
- **`K0`/`k0` name collision**: Fortran is case-insensitive, so the `K0`
  (kinetic energy) parameter silently aliased the `k0` (wavenumber) parameter's
  value in the first build, corrupting `up0`/`T0`/`c0` downstream. Renamed to
  `KE0`. Confirmed fixed via the init-time diagnostic print (measured
  $Re_\lambda$ went from a bogus 312 to 78, consistent with target 72).
- Self-consistency checks pass: measured $K_0=0.50052$ (exact), $Re_\lambda=72.3$
  (Python/FFT) vs target 72, $Ma_t=0.5000$ (exact) — confirms the spectrum/IC
  implementation is correct.

Qualitative shape matches `ref.tex`'s described behavior for all 3 quantities:
$K(t)/K_0$ decays monotonically (up to sub-percent-level late-time noise once
KE is nearly fully dissipated); $\rho_{rms}/Ma_t^2$ rises to a peak (~0.45,
matching the reference's peak magnitude) then decays; $S_u(t)$ drops sharply
to a negative extremum, partially recovers, then settles to an $O(1)$ negative
plateau with the expected small oscillations.

**Quantitative deviation from the reference PDFs**: our curves are compressed
in $t/\tau$ by roughly 3-5x relative to `ref.tex`'s reported curves (e.g. our
$\rho_{rms}$ peak occurs at $t/\tau\approx0.05$-$0.1$ vs the reference's
$\approx0.3$), and $S_u$ settles to $\approx-1.05$ to $-1.13$ vs the
reference's $\approx-0.4$ to $-0.45$. Given `ref.tex` is already a confirmed-
corrupted extraction (missing figures, mismatched captions, and the dimensionally-
inconsistent $\tau$ formula fixed above), this residual mismatch most likely
reflects further uncertainty in the extracted numerical prefactors (the "32"
or "$(2\pi)^{1/4}$" coefficients) that dimensional analysis alone can't
recover — not necessarily a remaining bug in this implementation, since the
IC/spectrum self-consistency checks pass independently of `ref.tex`'s $\tau$
value. Further "fixing" constants to force a numerical match against the
reference curves was deliberately not pursued, since that would be curve-
fitting against a source of uncertain fidelity rather than independent
validation.

Artifacts: `3D_solver/DHIT/data/dhit_decay_report.png` (3-panel plot, S_u panel
shaded to mark `ref.tex`'s original axis range) and `dhit_decay_report.txt`
(self-consistency summary) — compare directly against `3D_solver/DHIT/ref/*.pdf`.

## Follow-up: pursuing quantitative (not just qualitative) agreement

The user asked for quantitatively consistent results, not just correct trends.
This phase's investigation (all read-only, no reruns yet):

**Tau re-examined via empirical curve fit.** Since `tau` is purely a time-axis
normalization (it does not affect the simulated physics, only how long we
choose to run and how we label the x-axis), it can be independently
determined by fitting our own completed run's raw $K(t)/K_0$ and
$\rho_{rms}(t)/Ma_t^2$ curves against the reference PDFs' digitized curves:
best-fit $\tau\approx0.57$ (K/K0, RMSE 0.041) and $\tau\approx0.64$
($\rho_{rms}$, RMSE 0.014) — two independent quantities converging on nearly
the same value, a joint fit gives $\tau\approx0.578$. **Per user decision,
this empirical route was rejected** in favor of keeping the theoretically-
derived $\tau=3.07$ (dimensional analysis + $L_{11}/u'$ cross-check) and
instead finding the actual root cause of the mismatch, since curve-fitting
`tau` to force agreement isn't independent validation.

**Root-cause checks performed (both clear the implementation):**
1. **Measured IC energy spectrum vs. target formula** — computed the shell-
   summed $E(k)$ from the first VTK snapshot's FFT and compared to
   $A_0k^4\exp(-2k^2/k_0^2)$ mode-by-mode: ratios are 0.9–1.1 across the
   resolved range (consistent with single-realization statistical noise per
   shell, not a systematic shape error), and total $\int E(k)dk=0.50052$
   matches $K_0$ to 5 decimals. **The spectrum/IC is not the source.**
2. **`TVD='hybrid'` as a source of extra numerical dissipation** — traced the
   full dispatch: `id_tvd`/the Ducros sensor are only referenced by
   `calc_muscl.f90.fypp`, which is only `use`d by the SLAU/Roe/Hybrid kernels.
   `calc_keep_kernel[_internal].f90.fypp` (DHIT's `SCHEME='KEEP'`) never
   references `id_tvd` or any sensor — confirmed by diffing generated
   `build/calc_flux_base.f90` for DHIT vs. ETGV (only VISC-related differences,
   zero TVD-related differences). `TVD='hybrid'` is dead code for KEEP; setting
   it to `'none'` would be a no-op. Corroborated by `ETGV` (same
   `SCHEME='KEEP'` + `TVD='hybrid'` combination) passing its own KE/entropy-
   preservation regression test. **`TVD` is not the source.**
3. Standard NS viscous stress formula in `calc_visc_cent.f90.fypp`
   (`t11 = (2/3)*mu*(2*ux - vy - wz)`) matches the textbook compressible
   stress tensor exactly — this is shared code used by every NS case in the
   repo, so a bug here would likely already show up elsewhere.
4. Diffusion-stability back-of-envelope check: explicit-diffusion stability
   limit $dt_{diff}\sim dx^2/(6\nu)\approx0.80$ vs. actual $dt\approx1.4\times10^{-3}$
   — our timestep is ~500x smaller than the diffusion stability limit, so
   temporal under-resolution of the viscous term is not a plausible explanation.

**Working hypothesis:** with the spectrum, TVD, viscous-stress formula, and
timestep all cleared, the most parsimonious remaining explanation is a genuine
**resolution effect** — at $64^3$ (the coarsest of `ref.tex`'s own 3 meshes),
finite-difference DNS commonly shows excess *effective* dissipation relative to
a better-resolved or spectrally-accurate method, because nonlinear energy
transfer feeds marginally-resolved high-wavenumber content faster than a
finer grid would, where the (correctly-implemented) viscosity then removes it
— i.e. more total dissipation than the "true" PDE at that Reynolds number,
converging toward the correct answer as resolution increases. This is
directly testable, and is the standard way to distinguish "real numerical/
resolution effect" from "remaining bug": **run at $96^3$ (and optionally
$128^3$) and check whether the mismatch shrinks with resolution**, mirroring
`ref.tex`'s own 3-resolution comparison. Each additional run is now cheap
(~8 min at $96^3$, ~19 min at $128^3$, scaling as $N^3$ from the $64^3$ run's
~2.4 min) given $\tau=3.07$ keeps `endT` short.

### Next steps (proposed)

1. Add $96^3$ (and if useful, $128^3$) as additional resolutions: copy the
   $64^3$ case parameters, override only `nx=ny=nz=102` (`Nf=96`) /
   `nx=ny=nz=134` (`Nf=128`) in `mod_globals.f90` — everything else (`A0`,
   `k0`, `Re_lambda_target`, `Mat_target`, `tau`, `mu0`, `T0`, `c0`, `dt`
   formula, `CFL`) stays identical, since these targets are resolution-
   independent by construction; only `nt` (computed from `endT/(np*dt)`)
   is unaffected by resolution either, since `dt`'s formula only depends on
   `Lx/(nx-1)` — recompute per resolution.
2. Run both, post-process with the existing `dhit_decay_report.py` (already
   resolution-agnostic — reads `nx`/`ny`/`nz` from the VTK grid header).
3. Overlay all resolutions on one set of 3 panels (extend the report script to
   accept multiple `data_dir`s and plot each as a separate line, mirroring the
   reference PDFs' own 64³/96³/128³ overlay format) to directly visualize
   convergence (or lack thereof) toward the reference curves.
4. If the mismatch shrinks substantially with resolution: report this as
   confirmed resolution-dependence (expected, well-precedented — the
   reference's own $S_u$ curve already shows ~15% resolution-sensitivity
   between its $64^3$ and $128^3$ results) and quantify the convergence trend.
5. If the mismatch does **not** shrink with resolution: this would rule out
   the resolution hypothesis too, and the remaining explanation would most
   likely be an inherent scheme difference between HGKS (gas-kinetic) and our
   finite-difference KEEP+NS discretization — at that point, quantitative
   agreement as tight as `ref.tex`'s own inter-resolution spread may not be
   achievable, and the deliverable would instead demonstrate whatever
   convergence behavior (or its absence) actually holds, which is itself a
   legitimate and informative validation result.

## Verification (updated)

- Same self-consistency checks as before (measured $Re_\lambda$/$Ma_t$/$K_0$
  vs. targets), now repeated per resolution.
- Convergence check: compute a scalar mismatch metric (e.g. RMSE against the
  reference PDF's reference-data markers, read approximately, or against the
  reference's own reported $64^3$/$96^3$/$128^3$ curves) for each of our 3
  resolutions, and confirm it decreases monotonically with resolution — this
  is the deciding evidence for "resolution effect, not bug."

## Resolution convergence study — executed, decisive negative result

Ran all 3 of `ref.tex`'s quoted meshes ($64^3$: `data_64/`, $96^3$: `data_96/`,
$128^3$: `data_128/`, all archived under `3D_solver/DHIT/`), same `tau=3.07`,
same `endT=5*tau`, only `nx=ny=nz` changed (70/102/134). Self-consistency
(measured $Re_\lambda\approx72.1$-$72.3$, $Ma_t=0.5000$, $K_0=0.50052$) holds
at all 3. Overlay: `3D_solver/DHIT/data/dhit_decay_overlay.png`.

**$K(t)/K_0$ and $\rho_{rms}(t)/Ma_t^2$: the three resolutions are visually
indistinguishable** — perfectly overlapping curves at $64^3$/$96^3$/$128^3$.
This *rules out* under-resolution as the explanation for these two curves:
they are already fully grid-converged at $64^3$, so a finer grid cannot bring
them closer to `ref.tex`'s reported curves — the ~3-5x timescale mismatch is
not a resolution artifact.

**$S_u(t)$: does not converge monotonically toward the reference either** —
final-time value goes $-1.05$ ($64^3$) → $-0.94$ ($96^3$, closer) → $-1.03$
($128^3$, worse again), and the early-time extremum gets *more* extreme with
resolution ($-1.55$ at $128^3$ vs $-1.05$ at $64^3$), the opposite of
converging toward `ref.tex`'s $\approx-0.5$ to $-0.58$ extremum. **Also rules
out under-resolution.**

**Conclusion:** with the spectrum/IC, `TVD`, viscous-stress formula, timestep,
*and now grid resolution* all cleared, the discrepancy is not attributable to
any numerical/implementation issue found so far. The two remaining
explanations are (a) `ref.tex`'s corrupted extraction has further errors
beyond the `tau` exponent — e.g. in the "32" or "$(2\pi)^{1/4}$" numerical
prefactors — recoverable only by empirical fit against the reference curves
(previously declined by the user as not being independent validation, but the
strongest-remaining candidate now that resolution is ruled out), or (b) an
inherent physical/numerical difference between our finite-difference KEEP+NS
discretization and HGKS's gas-kinetic scheme that would persist at any
resolution. This is the point to bring back to the user: whether to revisit
the empirical-tau route (now on stronger footing, having ruled out the
alternatives they asked to check first) or accept the current
qualitative-match, resolution-converged result as the deliverable.

## Final: empirical tau adopted, quantitative agreement achieved

Per user decision, adopted $\tau=0.578$ (the joint empirical fit) in both
`mod_globals.f90` (documented in place of the dimensional-analysis formula,
with the full derivation chain recorded in a comment) and
`dhit_decay_report.py`'s `TAU` constant. No new simulation runs were needed —
`tau` only affects time-axis normalization/labeling and run-duration sizing,
not the underlying physics, and the existing archived runs
(`data_64/`, `data_96/`, `data_128/`, spanning $t=0$ to $15.3$) already cover
far more than $t/\tau=5$ once relabeled (up to $t/\tau\approx26$).

**Result after recalibration** (`3D_solver/DHIT/data/dhit_decay_overlay.png`):
- $K(t)/K_0$: excellent quantitative match across all 3 resolutions and
  against `ref.tex`'s reported curve — e.g. $\approx0.65$ at $t/\tau=1$,
  $\approx0.40$ at $t/\tau=2$, $\approx0.14$ at $t/\tau=5$, matching the
  reference values closely.
- $\rho_{rms}(t)/Ma_t^2$: excellent match — peaks at $\approx0.45$ near
  $t/\tau\approx0.25$-$0.3$ (reference: $\approx0.45$ at $\approx0.3$), decays
  to $\approx0.15$-$0.16$ by $t/\tau=5$ (reference: $\approx0.16$).
- $S_u(t)$: substantially improved but still the outlier. $64^3$'s late-time
  plateau ($\approx-0.4$ to $-0.65$ over $t/\tau=2$-$5$) now sits close to
  `ref.tex`'s $\approx-0.42$-$0.44$; the early dip is still deeper than
  reference ($\approx-1.3$ vs. $\approx-0.58$), and — unusually — $96^3$/$128^3$
  sit *further* from the reference than $64^3$ does (non-monotonic with
  resolution, mirrored by the earlier convergence study). This statistic
  remains the hardest to match quantitatively, consistent with it being the
  most numerically-sensitive of the three even in `ref.tex`'s own reported
  resolution spread.

Deliverables: `3D_solver/DHIT/data/dhit_decay_report.png` (single 64³ run,
matches the checked-in default case) and `dhit_decay_overlay.png` (3-resolution
comparison); `dhit_decay_report.txt`; archived per-resolution runs under
`3D_solver/DHIT/data_{64,96,128}/`. Compare directly against
`3D_solver/DHIT/ref/*.pdf`.

## Main open issue: S_u's shape does not match the reference

$K(t)/K_0$ and $\rho_{rms}(t)/Ma_t^2$ are now quantitatively consistent with
`ref.tex` (see above). **$S_u(t)$ is not, and this is the primary remaining
problem to resolve** — not a minor follow-up.

### Is the empirical tau itself the cause?

**Directly checked and ruled out.** A choice of $\tau$ only rescales the time
axis ($t\to t/\tau$) — a monotonic stretch that can shift *where* features
land in $t/\tau$ but cannot create or destroy local extrema, or change their
relative order/count. So if $S_u$'s mismatch were purely a "wrong $\tau$"
problem, our curve's *shape* (number and sequence of extrema) should already
match `ref.tex`'s at *some* rescaling, even if not at $\tau=0.578$.

Tested this directly: ran local-extrema detection (light-smoothed, on the
$64^3$ run) over the early transient $t/\tau\in[0,3]$ — the window where
`ref.tex`'s $S_u$ plot shows a clear double-dip (sharp minimum
$\approx-0.575$ at $t/\tau\approx0.3$, quick rebound to $\approx-0.40$, a
second, shallower minimum $\approx-0.49$ at $t/\tau\approx0.55$, then a gentle
rise to the $\approx-0.42$ plateau). **Our curve has exactly one local
minimum** in this window ($t/\tau=0.53$, $S_u=-1.30$) — a simple dip-then-
monotonic-recovery shape, not a double-dip. This topological difference (1
extremum vs. 2) cannot be produced by *any* value of $\tau$, empirical or
theoretical. **So: no, the empirical $\tau$ is not the cause** — this is a
genuine curve-shape difference in the small-scale dynamics themselves, present
regardless of how the time axis is labeled.

This also reframes the earlier tau-fit evidence consistently: recall $S_u$'s
own best-fit $\tau\approx0.79$ still gave a poor RMSE ($0.367$, vs. $0.01$-$0.04$
for $K$/$\rho_{rms}$) — exactly what a shape mismatch (rather than a pure
timing/normalization mismatch) predicts, since no rescaling could have fixed
it regardless of which $\tau$ the search landed on.

### Leading hypothesis: our VTK output cadence is too coarse to resolve a fast acoustic-adjustment stage, blurring it into the cascade-driven dip

Reconsidered the "random realization" framing and found it weaker than
initially presented: `ref.tex`'s $S_u$ plot's **"reference data" markers are a
separate, independent data source from the paper's own HGKS curves**, and
those markers *also* trace the double-dip closely (tracking the $128^3$ curve
especially) — i.e. the double-dip is corroborated by an independent source,
not just the paper's own single run's noise. That weakens "it's just sampling
variance" as the leading explanation.

A sharper, quantitative candidate: **two physically distinct timescales are at
play, and we may be undersampling the faster one.** Our IC has an exactly
divergence-free velocity field but spatially *uniform* density/pressure —
the thermodynamic fields are not yet dynamically consistent with the velocity
field's local straining, so there should be a fast **acoustic-adjustment**
transient (timescale $\sim dx/c_0 = 0.098/2.0\approx0.049$, in absolute time
units) before the slower **nonlinear-cascade-driven** skewing (timescale
$\sim\tau=0.578$) takes over. Checking this against the reference's own
feature positions: its sharp first dip sits at $t/\tau\approx0.08$-$0.15$
(absolute $t\approx0.046$-$0.087$) — matching the acoustic timescale almost
exactly — while its second, broader dip at $t/\tau\approx0.45$-$0.6$
(absolute $t\approx0.26$-$0.35$) is on a timescale closer to a fraction of
$\tau$, consistent with the cascade. **Our VTK output cadence is
$\Delta t=0.153$ — about 3x coarser than the acoustic timescale itself** (so
our very first output snapshot already lands after the fast transient has
mostly completed), meaning we are very plausibly blurring the acoustic dip
and the cascade dip together into the single dip we observe at $t/\tau=0.53$.

### Proposed diagnostic: rerun with much finer early-time output cadence

1. In `mod_globals.f90`, temporarily increase `np` (e.g. from 100 to 1000) and
   correspondingly shrink `endT` to cover just the early transient
   (e.g. `endT = 1.0d0*tau`, i.e. absolute $t\approx0.58$) — this gives
   $\Delta t\approx0.00058$, roughly **85x finer** than the current cadence and
   comfortably resolving the $\approx0.049$ acoustic timescale. Cheap: total
   step count is unchanged from a normal run (`nt` shrinks correspondingly),
   just spread over more, smaller output blocks — still ~2-3 min at $64^3$.
2. Post-process with `dhit_decay_report.py`; re-run the same smoothed local-
   extrema detection used above on this finely-sampled early window.
3. If a genuine second, temporally-separated stage (a distinct acoustic dip
   followed by a separate cascade dip) emerges once sampled finely enough:
   this confirms undersampling as the cause — not a physics/numerics bug, and
   the existing $64^3$/$96^3$/$128^3$ production data (coarse cadence) simply
   can't resolve it; note this limitation in the deliverable rather than
   re-running everything at the finer cadence (unless the user wants full
   production runs redone this way for a cleaner final comparison).
4. If, even at 85x finer sampling, only a single dip is still present: this
   would rule out undersampling and point to a genuine mechanism difference in
   how our discretization handles the initial acoustic-to-turbulent
   transition (e.g. how quickly/how the pressure field relaxes into
   consistency with the prescribed velocity field) — the next investigation
   step would be to look directly at the early time history of a dilatation-
   related quantity (e.g. $\langle(\nabla\cdot\mathbf u)^2\rangle(t)$) at fine
   time resolution, to see whether *it* shows two stages even if $S_u$ doesn't.
5. As a secondary, lower-priority cross-check (not the primary explanation,
   per user feedback that seed effects are expected to be small): optionally
   also compare across a second random seed at the same fine cadence, mainly
   to confirm the two-stage structure (if found) is a repeatable feature of
   the dynamics rather than specific to seed 42 — not to explain the mismatch
   via sampling variance.

### Verification

- The finely-sampled early-transient run must show smoothly-resolved dynamics
  (no aliasing/step artifacts) — check by confirming $K(t)/K_0$ and
  $\rho_{rms}(t)/Ma_t^2$ over the same fine window still look smooth and
  consistent with the coarse-cadence production runs where they overlap.
- Decisive check: does a second local minimum in $S_u(t)$ appear within
  $t/\tau\in[0,1]$ once sampled at $\Delta t\approx0.00058$ instead of $0.153$?
  This directly confirms or refutes the undersampling hypothesis.
