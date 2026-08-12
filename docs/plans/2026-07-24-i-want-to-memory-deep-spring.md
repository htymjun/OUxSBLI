# SoA Memory-Layout Refactor for `3D_solver/` and `3D_solver_curv/`

## Context

The conservative-variable array is currently `Q(nx,5,ny,nz)` (device) — the variable
dimension (rho, rho*u, rho*v, rho*w, rho*E) sits in the *second* position, so each
physical variable is strided by `nx` elements rather than living in its own
contiguous array. This hurts GPU memory coalescing in every convective/viscous
kernel, which reads all 5 components per cell. The goal is true SoA: separate device
arrays per variable, not merely a reordered single array (a transpose alone doesn't
change coalescing behavior for a genuinely split access pattern).

Scope is `3D_solver/` and `3D_solver_curv/` only. `2D_solver/` must not change, even
though it shares some top-level `src/` files with the 3D solvers.

## Confirmed scope decisions

1. **Split into 5 separate device arrays** (真 SoA): the conservative storage
   (`QJ`, `QJ2`, `QJs`) and the primitive/reconstructed buffer that feeds every
   convective/viscous kernel (currently allocated as `ruvwp`, but almost always
   received as a dummy argument literally named `Q` in kernel signatures).
2. **Transpose only** (single array, variable index moved to the last/slowest
   position, NOT split): flux arrays `E`/`F`/`G`, the RK residual accumulator `Rs`,
   and host-side `Q` (IC generation in `set.f90`, restart binary I/O in `main.f90`
   / `sbli.f90` / `main_curv.f90.fypp`).
3. MPI ghost-cell exchange (`calc_para.f90`, `calc_para_curv.f90`) keeps **one
   combined message per exchange direction** — pack from the 5 split source arrays
   into the same single 1‑D buffer, rather than sending 5 separate messages.
4. Restart-file (`.dat`) on-disk byte layout will change as an accepted consequence
   of transposing host `Q`. No migration shim — old restart files simply won't be
   binary-compatible after this change.

## Naming convention (apply uniformly in both solvers)

Split arrays: append `_1`..`_5` to the original name, each dimensioned `(nx,ny,nz)`,
in the fixed order `[rho, rho*u, rho*v, rho*w, rho*E]`:
- `QJ` → `QJ_1..QJ_5`; `QJ2` → `QJ2_1..QJ2_5`; `QJs` (Cartesian RK4 only) → `QJs_1..QJs_5`
- primitive/`ruvwp` buffer → `Q_1..Q_5` (this is what nearly every kernel signature's
  dummy arg `Q` becomes — since that's the array actually read in every hot loop).
  **Trim, don't uniformly pass all 5**: a kernel's signature should only take the
  component arrays it actually reads, not all 5 for uniformity. In particular,
  every viscous-related kernel never touches pressure (component 5) — confirmed:
  `calc_visc2.f90.fypp`/`calc_visc_high*`/`load_smem_visc*`/`calc_visc2_curv.f90`/
  `load_smem_visc2_curv.f90` read only `Q_2,Q_3,Q_4` (u,v,w, for stress tensors);
  `calc_hybrid.f90`/`calc_hybrid_curv.f90`'s Ducros sensor reads only `Q_2,Q_3,Q_4`;
  `calc_les.f90::calc_mut` reads `Q_1,Q_2,Q_3,Q_4` (needs rho too, for `mut`) but
  not `Q_5`. Only the convective kernels (KEEP/SLAU/Roe/Hybrid flux evaluation) and
  `calc_quantities_3D`/`_T_3D` genuinely need all 5. Trimming these signatures is a
  real (small) benefit of the split — fewer registers/dummy args per kernel launch
  — and should be done per-kernel, verified by checking which components each
  kernel body actually reads before finalizing its signature.
- `calc_steps.f90`/`calc_steps_curv.f90` local dummy names `Q2`, `Qin`, `Qout` →
  `Q2_1..Q2_5`, `Qin_1..Qin_5`, `Qout_1..Qout_5`
- `calc_rescale.f90` flattened buffers `Qm(ny*5)`, `Qre(ny*(nz-6)*5)` → `Qm_1..Qm_5`
  (each `ny`), `Qre_1..Qre_5` (each `ny*(nz-6)`), dropping the `5*(j-1)+l` index
  arithmetic entirely. `Qm_cpu`, `Qre_cpu` mirrors get the same `_1.._5` treatment.

Transpose-only arrays keep their original name, just reorder dims so the variable
index is last: `E(5,nx-1,ny-2,nz-2)` → `E(nx-1,ny-2,nz-2,5)` (same for `F`,`G`);
`Rs(nx-2,5,ny-2,nz-2)` → `Rs(nx-2,ny-2,nz-2,5)`; host `Q(nx,5,ny,nz)` →
`Q(nx,ny,nz,5)`.

**Key mechanical rule**: any `do l=1,5` loop (or Fortran colon-slice across the
variable dimension, e.g. `Q(1,:,j,k) = Q(nx-1,:,j,k)`) that touches a **split**
array must be manually unrolled into 5 explicit statements — there's no way to loop
over 5 differently-named arrays. A loop touching only **transpose-only** arrays
just needs its index reordered (`Rs(i,j,k,l)` instead of `Rs(i,l,j,k)`) — no
unrolling required. If a loop's body mixes both (common in the RK step kernels),
the whole loop must be unrolled, since the split-array part forces it.

## Verified corrections from investigation (apply these, not the naive assumption)

- `ruvwp` is currently allocated `(5,nx,ny,nz)` in `preprocess.f90.fypp` but every
  kernel consumes it as `(nx,5,ny,nz)` — a latent shape mismatch papered over by
  Fortran sequence association. The split naturally fixes this (each `Q_n` is
  simply `(nx,ny,nz)`, unambiguous).
- `set_bc_common.f90`'s `_init` routines (`set_bc_cyclic2/4/6_init`) run on the
  **host**, pre-split, single-array `Q` — confirmed no `device` attribute. Their
  colon-slice copies (`Q(1,:,j,k) = Q(nx-1,:,j,k)`, `Q(1:2,:,j,k) = Q(nx-3:nx-2,:,j,k)`,
  etc., ~21 sites across the three routines) need only a dimension reorder
  (`Q(1,j,k,:) = Q(nx-1,j,k,:)`), **not** per-component unrolling.
- The **GPU** variants in the same file (`set_bc_cyclic2/4/6`, `set_bc_cyclic_z`)
  operate on the split conservative array and do use real `do l=1,5` loops (many,
  large bodies for the 4th/6th-order corner cases) — these **do** need full
  unrolling.
- `set_bc_common.f90` is the literal same file compiled into both `3D_solver` and
  `3D_solver_curv` (confirmed via both CMakeLists), but curv has **no call sites**
  into it (`NACA/set.f90` and `CORN/set.f90` only call grid/metric setup, never
  `set_bc_cyclic*`). So this file only needs edits driven by `3D_solver`'s call
  sites; curv just needs it to keep compiling.
- In `3D_solver_curv`, `set_bc`'s `Q` dummy argument is confirmed device and bound
  to `QJ`/`QJ2` at the call site in `calc_time_dev_curv.f90.fypp` — i.e. it's the
  split conservative array, not the primitive buffer.
- `3D_solver/src/calc_steps_slice.f90` and `calc_steps_smem.f90` are confirmed dead
  code (not in any case's CMake `_BASE` list, would conflict with `calc_steps.f90`
  if compiled) — do not touch.
- `3D_solver_curv/src/main_curv.f90` (non-`.fypp`) is a stale, out-of-date tracked
  file not used by the CMake build (which regenerates from `main_curv.f90.fypp`) —
  do not touch; just make sure verification builds actually regenerate from the
  `.fypp` template.

## Implementation order

Work bottom-up (leaf kernels → dispatcher → orchestration → host) so most files are
mechanical, independent edits, and the few genuinely hard files (`calc_steps.f90`,
`calc_para.f90`, `calc_rescale.f90`, per-case `set_bc`) are done last against
already-finalized call signatures.

### Step 0 — shared coordination point
`src/calc_physical_quantities.f90`'s `calc_quantities_3D`/`calc_quantities_T_3D`
(NOT `calc_quantities_2D`/`calc_quantities_T_2D`, which stay untouched for
2D_solver) are called from both solvers. Change their signature once:
`(nx,ny,nz,Jacobian,QJ_1,QJ_2,QJ_3,QJ_4,QJ_5,Q_1,Q_2,Q_3,Q_4,Q_5,T,[mu,]k_lo,k_hi)`.
Do this first; both solvers' `calc_flux_base*.f90.fypp` call sites adapt to it.

### Step 1 — leaf kernel files (mutually independent, mechanical per-component rename)
Cartesian: `calc_keep_kernel.f90.fypp` (+`_internal`), `calc_slau_kernel.f90.fypp`
(+`_internal`), `calc_roe_kernel.f90.fypp` (+`_internal`), `calc_hybrid_kernel.f90.fypp`
(+`_internal`), `calc_hybrid.f90` (`calc_Ducros`), `calc_les.f90` (`calc_mut`),
`calc_visc2.f90.fypp`, `calc_visc_high.f90.fypp`/`_internal`, `load_smem_visc2.f90`,
`load_smem_visc_cent.f90.fypp`, `load_smem_visc_me4_base.f90`.
Curv equivalents: `calc_keep_kernel_curv.f90`, `calc_slau_kernel_curv.f90`,
`calc_hybrid_kernel_curv.f90`, `calc_hybrid_curv.f90`, `calc_visc2_curv.f90`,
`load_smem_visc2_curv.f90`.

Pattern in all of these: `rho(idx)=Q(i,1,j,k)` → `rho(idx)=Q_1(i,j,k)` (already
scalar, no loop); flux writes like `E(:,i,j-1,k-1) = KEEP2(...)*S` become
`E(i,j-1,k-1,:) = KEEP2(...)*S` (just move the colon — E/F/G going var-first→var-last
means these vectorized stores stay valid array-section assignments, no manual
unrolling needed). The `include`d small physics functions (`calc_keep_3d.f90`,
`calc_slau_3d.f90`, `calc_roe_3d.f90`, `calc_visc_cent.f90.fypp`,
`calc_visc_me4_base.f90`) operate only on small local vectors, never touch
`Q`/`E`/`F`/`G` directly — leave untouched.

Apply the trimming rule from the naming convention here: `calc_visc2.f90.fypp`,
`calc_visc_high.f90.fypp`/`_internal`, `load_smem_visc2.f90`,
`load_smem_visc_cent.f90.fypp`, `load_smem_visc_me4_base.f90`,
`calc_visc2_curv.f90`, `load_smem_visc2_curv.f90` should only declare `Q_2,Q_3,Q_4`
dummy args (no `Q_1`/`Q_5` — pressure and density are never read in these files,
verified during exploration). `calc_hybrid.f90`/`calc_hybrid_curv.f90`'s
`calc_Ducros`/`calc_Ducros_curv` similarly take only `Q_2,Q_3,Q_4`. `calc_les.f90`'s
`calc_mut` takes `Q_1,Q_2,Q_3,Q_4` (needs rho for the `mut` scaling) but not `Q_5`.
The convective kernels (`calc_keep_kernel*`, `calc_slau_kernel*`, `calc_roe_kernel*`,
`calc_hybrid_kernel*`) and `calc_quantities_3D`/`_T_3D` (Step 0) take all 5, since
flux evaluation genuinely needs the full state.

### Step 2 — dispatchers
`calc_flux_base.f90.fypp` (Cartesian) / `calc_flux_base_curv.f90.fypp` (curv): thread
the split argument lists through `calc_conv`/`calc_conv_curv`, `calc_EFG`/`calc_EFG_curv`,
`calc_EFG_halo` (COMMZ path), and the `calc_quantities_3D`/`_T_3D` calls from Step 0.
Pure argument-list fan-out, no loops of its own.

### Step 3 — allocation/preprocessing
`preprocess.f90.fypp` / `preprocess_curv.f90.fypp`: split `ruvwp` allocation into
`Q_1..Q_5(nx,ny,nz)`; `E`/`F`/`G` allocations transpose-only. In `pre_calc`/
`pre_calc_curv`: the Jacobian-normalization loop on host `Q` is transpose-only
(reorder, no unroll); the `QJ = Q` whole-array host→device copy becomes 5 explicit
statements `QJ_1 = Q(:,:,:,1)`, …, `QJ_5 = Q(:,:,:,5)` — because host `Q` is
transposed to put the variable index last, `Q(:,:,:,m)` is a contiguous host block
matching contiguous device `QJ_m`, so this is a clean, efficient per-component
transfer, not a strided one. `pre_rescale` (Cartesian only): split `Qre`/`Qm` into
`_1..5`.

### Step 4 — print/VTK output
`src/print.f90` (3D-specific routines only: `make_1d_for_print3`,
`send_recv_for_print_even3`/`_odd3` — do NOT touch the 2D-specific routines in the
same file) and `src/print_curv.f90` (curv-only file): rename `QJ` reads to
`QJ_1..QJ_5`; `Q = QJ` whole-array device→host copies become 5 explicit per-component
copies, same contiguity argument as Step 3.

### Step 5 — RK time-stepping (highest-effort files, real unrolling)
`3D_solver/src/calc_steps.f90` / `3D_solver_curv/src/calc_steps_curv.f90`:
- `calc_R`/`calc_R_curv`: `E`/`F`/`G` reads are already scalar per-component, just
  reorder the index (`E(i,j,k,1)` instead of `E(1,i,j,k)`) — no unroll.
- `calc_step1`(`_curv`), `calc_step2_3`(`_curv`): unroll the `do l=1,5` loop that
  writes the split `Q2`/`Qout` array (5 explicit statements each).
- `calc_step`(`_curv`), `calc_step4`(`_curv`) (RK4 path — dead in all current cases
  except Cartesian ETGV uses RK4; curv's NACA/CORN are RK3 so this path is dead
  there but must still compile/be logically correct): loop touches both split
  `Q`/`Q2` and transpose-only `Rs` — unroll the whole loop, reindexing `Rs` to
  `Rs(i,j,k,l)` inline.

### Step 6 — MPI ghost-cell exchange
`3D_solver/src/calc_para.f90` / `3D_solver_curv/src/calc_para_curv.f90`: every
`flatten*`/`reconstruct*` routine (12 loops in Cartesian: `flatten`, `flatten_left`,
`flatten_right`, `flatten_rescale`, `reconstruct`, `reconstruct_left`,
`reconstruct_right`, `reconstruct_sbli_inlet`, `flatten_z_lo/hi`,
`reconstruct_z_lo/hi`; 4 loops in curv: `flatten_z_lo/hi`, `reconstruct_z_lo/hi`,
currently dead since `COMMZ=False` in all curv cases but must still be correct)
must have its `do l=1,5` loop unrolled into 5 explicit statements sourcing from
`QJ_1..QJ_5`, writing into the **same single 1‑D buffer** at component-block
offsets (e.g. component `m` occupies offset `(m-1)*nx*ny*overlap` within the slab
buffer) — collapsing a `do(4)` cuf-kernel to `do(3)` plus 5 literal statements.
The actual `MPI_SENDRECV`/`cudaMemcpy` calls on the flat buffers are unaffected
(same wire format/message count).

### Step 7 — rescale (Cartesian SBLI/TBL only, `calc_rescale.f90`)
Split `Qm`/`Qre` per the naming convention; remove the `5*(j-1)+l` /
`ny*5*(k-1)+5*(j-1)+l` index arithmetic throughout `calc_mean`, `copy`,
`step_rescale`, `write_Qm`, `rescale_recv_send`. For `set_rescale` (host-only,
currently the most arithmetic-heavy routine): reshape the incoming `Qre_1..5`
dummy args to explicit `(ny,nz)` shape locally — same underlying contiguous
memory, but turns most of the flat-offset arithmetic into plain 2‑D per-component
array accesses. For the one MPI send/recv of `Qre`/`Qm` in `step_rescale`/
`rescale_recv_send` (keeping "one combined message" per the scope decision): add a
small pack/unpack step copying `Qre_1..5`/`Qm_1..5` into a scratch flat buffer
immediately before/after the `CPUGPU_MPI_SEND`/`_RECV` call, rather than changing
the generic MPI wrapper in `cpu_gpu_mpi.f90` (which must stay untouched — it only
ever takes one flat buffer and is shared with 2D_solver).

### Step 8 — boundary conditions
- `set_bc_common.f90` (Cartesian call sites only, per corrections above): GPU
  `set_bc_cyclic2/4/6`, `set_bc_cyclic_z` — unroll every `do l=1,5` loop (largest
  volume of repetitive-but-mechanical unrolling in the whole refactor, especially
  the 4th/6th-order corner cases). CPU `_init` variants — reorder only, no unroll.
  `set_bc_mut_common` — untouched (operates on `mut`/`qc2`, not `Q`).
- `set_bc_tbl_sbli.f90` — already fully unrolled per literal index; trivial rename.
- Per-case `set.f90` (all 8 Cartesian cases + NACA/CORN): `set_init`/IC generation
  is host, transpose-only (reorder). `set_bc` is device, split — most bodies are
  already per-component literal statements (trivial rename), except:
  - `IVST/set.f90::set_bc` has 2 genuine `do l=1,5` loops around a column-cache
    (`Qc` also needs splitting into `Qc_1..5`).
  - `SBLI/set.f90::set_bc`/`set_bc_Gaussian` have 5 loops total, one reading `Qre`
    via flat-index arithmetic (coordinate with Step 7's `Qre` split).
  - `STZ/set.f90::set_bc` has 2 loops (z-lo/z-hi zero-gradient fill).
  - `TBL/set.f90::set_bc` has 3 loops, one reading `Qre`.
  - `NACA/set.f90` and `CORN/set.f90::set_bc` — nearly all per-literal-index
    (trivial), except the z-periodic full-slab copy `Q(:,:,:,1) = Q(:,:,:,nz-1)` /
    `Q(:,:,:,nz) = Q(:,:,:,2)` which, since this `Q` is the split device array,
    expands to 5 pairs of per-component slab copies.
  - `DHIT`, `ETGV`, `KHI`, `NSTGV` — no loops of their own; just expand the call
    into `set_bc_cyclic` to pass 5 arrays.

### Step 9 — top-level orchestration and host I/O
`3D_solver/src/calc_time_dev.f90.fypp` (6 fypp-conditional `RungeKutta` variants)
and `3D_solver_curv/src/calc_time_dev_curv.f90.fypp` (2 variants): update every
allocation and call site to the finalized split signatures — no new logic, pure
argument-list fan-out, done last since every callee signature is now stable.
`main.f90`, `src/sbli.f90`, `main_curv.f90.fypp`: host `Q(nx,5,ny,nz)` →
`Q(nx,ny,nz,5)`; restart read/write stays whole-array (byte layout changes,
accepted); the pre-write Jacobian-rescale loop is transpose-only (reorder, no
unroll). `set_init_common.f90::set_init_tbl` — host, transpose-only reorder.

## Files explicitly out of scope — do not touch

- `2D_solver/` entirely, including its own `set_bc_common.f90` (a distinct file
  from the 3D one) and the 2D-specific routines inside shared files
  (`calc_quantities_2D`/`_T_2D` in `calc_physical_quantities.f90`;
  `make_1d_for_print2`, `send_recv_for_print_even2`/`_odd2`, `print0`,
  `print_vtk_2D` in `print.f90`).
- `src/cpu_gpu_mpi.f90` — generic flat-buffer MPI wrapper, shared with 2D_solver,
  layout-agnostic; no changes needed.
- `src/calc_muscl.f90.fypp` — shared with 2D_solver; operates only on small local
  stencil arrays already extracted by the caller, never touches `Q` directly.
- `src/set_coordinate.f90`, `src/set_compressible_bl.f90` — grid metrics /
  Blasius-profile generators; confirmed no `Q` access.
- `3D_solver/src/calc_steps_slice.f90`, `calc_steps_smem.f90` — dead code, not
  compiled into any case.
- `3D_solver_curv/src/main_curv.f90` (non-fypp stale duplicate) — not part of the
  build graph.
- `calc_keep_3d.f90`, `calc_slau_3d.f90`, `calc_roe_3d.f90`, `calc_visc_cent.f90.fypp`,
  `calc_visc_me4_base.f90` (as `include`d physics bodies) — operate on small local
  vectors only.
- `mod_constant.f90.fypp`'s `Normal_x/y/z(5)` one-hot vectors — unrelated to `Q`'s
  storage layout, used only as function arguments.

## Verification

1. **Build all 8 Cartesian cases** (`DHIT`, `ETGV`, `IVST`, `KHI`, `NSTGV`, `SBLI`,
   `STZ`, `TBL`) plus **both curv cases** (`NACA`, `CORN`) from clean build
   directories (`cmake -B build && cmake --build build -j`), in this rough order to
   surface signature bugs early on the simplest configs first: `IVST` (Euler, RK3,
   SLAU, no BCs beyond periodic) → `ETGV` (Euler, RK4 — only case using RK4, so
   uniquely exercises `calc_step`/`calc_step4`) → `KHI` (Euler, KEEP) → `DHIT`/`NSTGV`
   (NS, viscous) → `TBL` (LES, Hybrid scheme) → `SBLI` (NS, RESCALE=True — the
   highest-risk file, `calc_rescale.f90`) → `STZ` (COMMZ=True — exercises the
   z-halo pack/unpack path) → `NACA`/`CORN`.
2. **Automated regression tests** — `pytest ouxsbli/tests/test_etgv.py` (only
   pytest test targeting `3D_solver/`; checks entropy conservation and KE bounds
   over 20000 steps — very sensitive to any component-index mixup from the
   unrolling) and `pytest ouxsbli/tests/test_corn.py` (only test targeting
   `3D_solver_curv/`; checks post-shock pressure/density ratio against
   oblique-shock theory at 1% tolerance). Both should pass with no tolerance
   relaxation, since this is a pure relayout with no intended arithmetic change —
   any discrepancy beyond floating-point reassociation noise indicates a
   component-index bug in one of the manually-unrolled loops or a mismatched
   argument order at a call site.
3. **Manual smoke runs** for paths the two pytest tests don't reach: `SBLI` (live
   `RESCALE` + `Qre` inlet path), `STZ` with `nranks>=4` (live COMMZ z-halo
   exchange), `TBL` (LES + Hybrid scheme), `NACA` (no dedicated pytest test at all
   — flag this as a pre-existing gap, not something to newly fix here). Check
   `data/entropy.d`/`data/kinetic_energy.d` for boundedness/no-NaN, and do a
   restart round-trip (`RESTART=True` after a run) to confirm the transposed host
   `Q` layout is self-consistent even though the on-disk format changed.
4. **Roe scheme**: no case's `config.fypp` sets `SCHEME='Roe'`, so it's compiled
   everywhere but never runtime-exercised. As a one-off manual check (not
   committed), temporarily set `SCHEME='Roe'` in one case (e.g. `IVST`), confirm it
   builds and runs a few steps without NaN, then revert.
