# WENO FP32/FP64 microbenchmark

This directory isolates WENO5-Z reconstruction from the 1D SLAU solver so the
FP32/FP64 simultaneous-execution question can be measured without SLAU flux,
Runge-Kutta, or primitive-variable setup noise.

The important comparisons are:

| Mode pair | What it tests |
|---|---|
| `var_*_seq` vs `var_*_warp` | Physical-variable split. Selected variables use FP64 WENO; the rest use FP32 WENO. This mirrors the earlier SLAU-WENO split, but removes the SLAU flux. |
| `weight_poly32_seq` vs `weight_poly32_warp` | Expression-block split. WENO-Z nonlinear weights use FP64; candidate polynomials use FP32. |
| `weight_poly32_serial` vs `weight_poly32_warp` | Same 256-thread/shared-memory layout. `serial` computes FP64 weights and FP32 polynomials sequentially on the lower warp role, while `warp` computes them concurrently on lower/upper roles. This is the cleanest test of FP32/FP64 warp-level overlap. |
| `weight_poly32_serial_oncebar` vs `weight_poly32_warp_oncebar` | Same full-shared-memory layout as above, but the block barrier is moved outside the inner repeat loop. Use `--nrepeat 2`, `4`, or `8` to see whether barrier fixed cost hides FP32/FP64 overlap. |
| `weight_poly32_wsmem_serial` vs `weight_poly32_wsmem_warp` | Reduced shared-memory layout. Only FP64 WENO weights are stored in shared memory; FP32 polynomials remain in registers across the barrier. This lowers static shared memory from the full weights+polynomials layout. |
| `weight_poly32_wsmem_tile2_serial` vs `weight_poly32_wsmem_tile2_warp` | Reduced shared-memory layout with two faces per thread role. This amortizes one barrier over two WENO reconstructions. |
| `weight_poly32_halfwarp_serial` vs `weight_poly32_halfwarp_shfl` | Shared-memory-free diagnostic using half-warp lane split and `__shfl_xor`. This removes smem/barrier cost, but because one warp is lane-divergent it is not the main cross-warp simultaneous-issue strategy. |

Run the default WENO weight/poly study on A100 with:

```bash
cd /path/to/OUxSBLI
bash 1D_solver/microbenchmark/run_weno_nsys.sh \
  --out 1D_solver/report/weno_micro_a100 \
  --gpu-cc 80 \
  --repeat 3
```

By default this runs the FP32/FP64 overlap and smem-cost modes with
`--nrepeat-list "1 2 4 8"`. The older `var_*` physical-variable split modes are
still available through `--modes`, but are not included in the default run.

For a single-repeat quick smoke run:

```bash
bash 1D_solver/microbenchmark/run_weno_nsys.sh \
  --out 1D_solver/report/weno_micro_a100_smoke \
  --gpu-cc 80 \
  --repeat 1 \
  --nrepeat-list "1"
```

The script writes `summary.csv`, raw `nsys-rep`, SQLite exports, and
`nsys stats` CSV files under `--out`.
