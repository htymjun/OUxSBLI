"""Cost of the KEEP precision modes on ETGV (all-periodic Euler).

For each mode and KEEP order: build 3D_solver/ETGV at n^3 through ouxsbli.Case,
run a few RK4 steps under nsys, and report
  * GPU time per step of the convective KEEP kernels and of all kernels,
  * a SASS census of the x-direction KEEP kernel: FP64-pipe instructions
    (D*, F2F/I2F/F2I involving F64, as in diffFP32's sass_census) and FP32 ones,
  * the census ratio to FP64 under a throughput model of an FP64:FP32 = 1:k
    machine (cost ~ n_fp64 + n_fp32 / k), k = 64 here, k = 2 on GH200/H100.

  python -m ouxsbli.analysis.bench_prec --n 256 --steps 20 --out bench.md
"""
import argparse
import csv
import io
import os
import re
import subprocess
from pathlib import Path

from ouxsbli import Case

MODES = {"FP64": {}, "DF": {"prec": "df"}, "MIXED": {"keep_prec": "mixed"}}
FP64_OP = re.compile(r"\b(DADD|DMUL|DFMA|DSETP|DMNMX)\b|\b[FI]2[FI]\.[^ ;]*64")
FP32_OP = re.compile(r"\b(FADD|FMUL|FFMA|FSETP|FMNMX|FSEL|MUFU)\b")


def census(exe, kernel):
    sass = subprocess.run(["cuobjdump", "-sass", str(exe)], capture_output=True, text=True).stdout
    body = ""
    for part in sass.split("Function : ")[1:]:
        if part.split("\n", 1)[0].strip().startswith(kernel):
            body = part
            break
    return len(FP64_OP.findall(body)), len(FP32_OP.findall(body))


def kernel_times(rep):
    out = subprocess.run(["nsys", "stats", "-r", "cuda_gpu_kern_sum", "--format", "csv", rep],
                         capture_output=True, text=True).stdout
    rows = list(csv.DictReader(io.StringIO(out[out.index("Time (%)") if "Time (%)" in out else 0:])))
    return {r["Name"]: float(r["Total Time (ns)"]) for r in rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=256)
    ap.add_argument("--steps", type=int, default=20)
    ap.add_argument("--orders", default="2,6")
    ap.add_argument("--modes", default=",".join(MODES))
    ap.add_argument("--work", default="./tmp/bench_prec")
    ap.add_argument("--out", default="bench_prec.md")
    a = ap.parse_args()
    lines = [f"ETGV {a.n}^3, RK4, {a.steps} steps, 2 ranks; GPU time per step [ms]\n",
             "| KEEP | mode | KEEP kernels | all kernels | KEEP / FP64 | FP64 ops | FP32 ops | 1:64 model | 1:2 model |",
             "|---|---|---|---|---|---|---|---|---|"]
    for order in map(int, a.orders.split(",")):
        ref = None
        for mode in a.modes.split(","):
            wd = Path(a.work) / f"{mode.lower()}_{order}"
            case = Case(source="3D_solver/ETGV", workdir=str(wd), visc="euler", scheme="keep",
                        accuracy=order, rk="rk4", nx=a.n + 2, ny=a.n + 2, nz=a.n + 2,
                        nt=a.steps, np=1, dt=1.e-3, **MODES[mode])
            case.build()
            (wd / "data").mkdir(exist_ok=True)
            rep = str(wd / "prof")
            subprocess.run(["nsys", "profile", "-t", "cuda", "--force-overwrite", "true", "-o", rep,
                            "mpirun", "-n", "2", "./build/a.out"], cwd=wd, check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            t = kernel_times(rep + ".nsys-rep")
            keep = sum(v for k, v in t.items() if "calc_keep" in k) / a.steps / 1e6
            total = sum(t.values()) / a.steps / 1e6
            n64, n32 = census(wd / "build" / "a.out", "calc_keep_kernel_internal_calc_keep_x_in_")
            if ref is None:
                ref = (keep, n64, n32)
            m64 = (n64 + n32 / 64) / (ref[1] + ref[2] / 64)
            m2 = (n64 + n32 / 2) / (ref[1] + ref[2] / 2)
            lines.append(f"| {order} | {mode} | {keep:.3f} | {total:.3f} | {keep / ref[0]:.2f} "
                         f"| {n64} | {n32} | {m64:.2f} | {m2:.2f} |")
            print(lines[-1], flush=True)
    Path(a.out).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
