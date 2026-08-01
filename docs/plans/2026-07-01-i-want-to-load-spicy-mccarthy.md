# Sliding-window G-flux loading in a new `calc_steps_slice.f90`

## Context

The prior shared-memory experiment (`3D_solver/src/calc_steps_smem.f90`, staged but not wired into the build — you manage swapping it in manually) stages the whole per-block G tile into shared memory via a cooperative load before computing the z-flux divergence. You measured it as "as fast as" the naive direct-global-read baseline (`calc_steps.f90`). That file stays as-is (kept for reference/comparison) — this next experiment goes in a **new** file, `3D_solver/src/calc_steps_slice.f90`.

The reason: the per-block tile only lets threads reuse G *within* a block. For every case except ETGV, `threads%z == 1` (the CUDA block is only 1 thread deep in z), so there is zero in-block reuse — every G plane still gets read from global memory exactly twice (once as its own block's "current" plane, once as the neighboring z-block's halo), identical to the naive baseline's redundancy. Even ETGV's `threads%z=2` only saves a fraction of that.

The fix you asked for is the classic GPU stencil "z-sweep"/sliding-window pattern: collapse the z-tiling entirely (`gridDim%z=1`), have each `(i,j)` thread walk the **entire** z-column (`k=1..nz-2`) in an internal loop, and keep only the 2 planes it actually needs (`G_before`, `G_after`) alive at once — sliding the window forward by exactly one new global read per iteration. This reads every G plane from global memory exactly once, total, regardless of grid size, eliminating the block-boundary redundancy that made the block-tiled version a no-op.

You flagged the risk correctly: naively holding `G_before(5)`/`G_after(5)` as plain local (register) arrays across a loop that can run 500+ iterations (e.g. NSTGV, nz=513) risks the compiler blowing up register usage / hurting occupancy. The fix is to declare them `shared`, indexed per-thread by `(it,jt)` — this is **not** cross-thread cooperative loading (each thread only ever touches its own `(it,jt)` slot), it's a deliberate manual "spill to shared memory" to keep them out of the register file regardless of loop length. Because of that — no thread ever reads another thread's slot — this design needs **no `syncthreads()` at all**, and the early-return-before-a-barrier hazard from the block-tiled design doesn't apply here: `if (i > nx-2 .or. j > ny-2) return` up front is completely safe.

## Change

### 1. `3D_solver/src/calc_steps_slice.f90` — new file

A fresh `module calc_steps` (same module name as `calc_steps.f90`/`calc_steps_smem.f90`, matching the existing convention of swapping the physical file in manually — no CMake changes, per your instruction). No `load_smem_G`, `load_smem_G_cooperative`, the `Gs(5,threads%x,threads%z+1,threads%y)` tile, the `valid`/`syncthreads()` dance, or `use wmma` (no `pipelineMemcpyAsync`/`pipelineCommit`/`pipelineWaitPrior` — plain synchronous loads replace them, since there's no cooperative load to pipeline against anymore, and no async overlap to arrange).

The `calc_R_EF`/`calc_R_G` split in `calc_steps_smem.f90` existed specifically to let `calc_R_EF` run while the async G-tile pipeline copy was in flight (`pipelineWaitPrior` sat between them). That rationale is gone here — there's no async pipeline left to overlap against — so go back to a single combined subroutine (matching the original, pre-split `calc_R` shape), just with `G` replaced by the two window planes:
```fortran
!$dir inline
attributes(device) subroutine calc_R(nx, ny, nz, i, j, k, dtdxdy, dtdydz, dtdzdx, E, F, G_before, G_after, R)
  integer, intent(in), value              :: nx, ny, nz
  integer, intent(in), value              :: i, j, k
  real(8), intent(in), value              :: dtdxdy, dtdydz, dtdzdx
  real(8), intent(in), device, contiguous :: E(5,nx-1,ny-2,nz-2)
  real(8), intent(in), device, contiguous :: F(5,nx-2,ny-1,nz-2)
  real(8), intent(in), contiguous         :: G_before(5), G_after(5)
  real(8), intent(out), contiguous        :: R(5)
  ! x direction
  R(1) = dtdydz * (-E(1,i,j,k) + E(1,i+1,j,k))
  R(2) = dtdydz * (-E(2,i,j,k) + E(2,i+1,j,k))
  R(3) = dtdydz * (-E(3,i,j,k) + E(3,i+1,j,k))
  R(4) = dtdydz * (-E(4,i,j,k) + E(4,i+1,j,k))
  R(5) = dtdydz * (-E(5,i,j,k) + E(5,i+1,j,k))
  ! y direction
  R(1) = R(1) + dtdzdx * (-F(1,i,j,k) + F(1,i,j+1,k))
  R(2) = R(2) + dtdzdx * (-F(2,i,j,k) + F(2,i,j+1,k))
  R(3) = R(3) + dtdzdx * (-F(3,i,j,k) + F(3,i,j+1,k))
  R(4) = R(4) + dtdzdx * (-F(4,i,j,k) + F(4,i,j+1,k))
  R(5) = R(5) + dtdzdx * (-F(5,i,j,k) + F(5,i,j+1,k))
  ! z direction (sliding window instead of global G)
  R(1) = R(1) + dtdxdy * (-G_before(1) + G_after(1))
  R(2) = R(2) + dtdxdy * (-G_before(2) + G_after(2))
  R(3) = R(3) + dtdxdy * (-G_before(3) + G_after(3))
  R(4) = R(4) + dtdxdy * (-G_before(4) + G_after(4))
  R(5) = R(5) + dtdxdy * (-G_before(5) + G_after(5))
end subroutine calc_R
```

Restructure each of the four kernels (`calc_step1`, `calc_step`, `calc_step2_3`, `calc_step4`) the same way — shown for `calc_step1`, apply identically to the other three (only the per-stage update math inside the loop body differs, exactly as it already differs today):
```fortran
attributes(global) subroutine calc_step1(nx, ny, nz, coef, dtdxdy, dtdydz, dtdzdx, E, F, G, Q, Q2)
  ... dummy args unchanged (E, F, G, Q, Q2 still plain global 4D arrays) ...
  real(8), shared :: G_before(5,threads%x,threads%y)
  real(8), shared :: G_after(5,threads%x,threads%y)
  real(8) R(5), coef_dtdxdy, coef_dtdydz, coef_dtdzdx
  integer i, j, k, l, it, jt
  it = threadIdx%x
  jt = threadIdx%y
  i  = (blockIdx%x-1)*blockDim%x + it
  j  = (blockIdx%y-1)*blockDim%y + jt
  if (i > nx-2 .or. j > ny-2) return
  G_before(1,it,jt) = G(1,i,j,1)
  G_before(2,it,jt) = G(2,i,j,1)
  G_before(3,it,jt) = G(3,i,j,1)
  G_before(4,it,jt) = G(4,i,j,1)
  G_before(5,it,jt) = G(5,i,j,1)
  do k = 1, nz-2
    G_after(1,it,jt) = G(1,i,j,k+1)
    G_after(2,it,jt) = G(2,i,j,k+1)
    G_after(3,it,jt) = G(3,i,j,k+1)
    G_after(4,it,jt) = G(4,i,j,k+1)
    G_after(5,it,jt) = G(5,i,j,k+1)
    coef_dtdxdy = coef * dtdxdy(i,j)
    coef_dtdydz = coef * dtdydz(j,k)
    coef_dtdzdx = coef * dtdzdx(i,k)
    call calc_R(nx, ny, nz, i, j, k, coef_dtdxdy, coef_dtdydz, coef_dtdzdx, E, F, G_before(:,it,jt), G_after(:,it,jt), R)
    do l = 1, 5
      Q2(i+1,l,j+1,k+1) = Q(i+1,l,j+1,k+1) - R(l)
    enddo
    G_before(1,it,jt) = G_after(1,it,jt)
    G_before(2,it,jt) = G_after(2,it,jt)
    G_before(3,it,jt) = G_after(3,it,jt)
    G_before(4,it,jt) = G_after(4,it,jt)
    G_before(5,it,jt) = G_after(5,it,jt)
  enddo
end subroutine calc_step1
```
`G_before(:,it,jt)`/`G_after(:,it,jt)` are contiguous 5-element sections (component is the leading/fastest dimension), matching the `contiguous` dummies in `calc_R`. `E`/`F` handling inside `calc_R` is otherwise unchanged from today — still a direct global read per iteration.

`calc_step`/`calc_step4` additionally index their `Rs(nx-2,5,ny-2,nz-2)` accumulator by the loop's `k`, same as they already do; `calc_step2_3` is the same pattern with its 3-coefficient blend.

### 2. `3D_solver/src/calc_time_dev.f90.fypp` — collapse the launch grid's z-dimension

All 21 occurrences of `<<<blocks,threads>>>` in this file are exclusively the `calc_step1`/`calc_step`/`calc_step2_3`/`calc_step4` launches (verified — no other kernel type uses this exact literal pair). Since each thread now walks the whole z-column itself, the launch must become 2D:
```fortran
<<<dim3(blocks%x,blocks%y,1),dim3(threads%x,threads%y,1)>>>
```
`blocks%x`/`blocks%y`/`threads%x`/`threads%y` are already correctly sized for the x,y tiling (from `set_coordinate.f90`'s `set_block_3D`); only the z-component is overridden to 1 here, inline, so no changes are needed to `mod_globals.f90` or `set_coordinate.f90` in any case directory.

**Caveat worth flagging, not fixing now**: collapsing `gridDim%z` shrinks the total block count for small grids. ETGV (66³, `threads=(32,4,2)`) goes from `2×16×32=1024` blocks down to `2×16=32` — likely too few to fully occupy a modern GPU, so ETGV specifically may trade memory-traffic reduction for reduced parallelism. Larger cases (e.g. NSTGV, 513³) keep thousands of blocks even after collapsing z, so shouldn't be affected. Worth checking profiler occupancy on both a small and large case rather than assuming the traffic win dominates everywhere.

## Verification

Since `calc_steps_slice.f90` isn't wired into the build (by your choice — no CMake changes), test it by temporarily swapping it in for `calc_steps.f90`, then restoring via git afterward:
```bash
cp 3D_solver/src/calc_steps.f90 /tmp/calc_steps_baseline.f90.bak   # your manual process; not committing anything
cp 3D_solver/src/calc_steps_slice.f90 3D_solver/src/calc_steps.f90
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/comm_libs/openmpi4/bin:/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/compilers/bin:$PATH
cd 3D_solver/ETGV && rm -rf build && FC=mpif90 cmake -B build && cmake --build build -j
cd build && mkdir -p data && mpirun -n 2 ./a.out
tail data/entropy.d data/kinetic_energy.d   # compare against the ETGV run already verified this session (entropy ~1e-4, KE decaying 0.02->0.0152)
```
Then restore the baseline (`git checkout -- 3D_solver/src/calc_steps.f90`) and repeat on a larger case (e.g. NSTGV, if it fits this machine's 8GB GPU at reduced grid size, or IVST as a quick smoke test) to sanity-check the occupancy caveat above doesn't cause a regression there. Compare `nsys` output (`profile.sh` in each case dir) between the naive baseline and this version to confirm reduced global memory traffic on the G array is actually showing up, not just "didn't crash."
