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

Run on A100 with:

```bash
cd /path/to/OUxSBLI
bash 1D_solver/microbenchmark/run_weno_nsys.sh \
  --out 1D_solver/report/weno_micro_a100 \
  --gpu-cc 80 \
  --repeat 3
```

The script writes `summary.csv`, raw `nsys-rep`, SQLite exports, and
`nsys stats` CSV files under `--out`.
