"""How much of the fluid that carried the initial disturbance is still in the
domain / in the separation bubble of a 3D_solver/SWLBLI run with TRACER=True.

The solver (3D_solver/src/tracer.f90) seeds Lagrangian particles at grid nodes
of the disturbed boundary layer and writes ``data/data_trace/P######.dat`` every
``trace_nstep`` steps. Each particle stands for the fixed fluid mass
``m_p = rho(x_0) dV_0`` of the node it was seeded at (``dV_0`` = node control
volume times the seeding strides in x, y; in z the period over the seeded
planes), so sums of ``m_p`` are Lagrangian estimates
of a passive scalar ``c`` that was 1 in the seeded region at t = 0:

  f_mass(t)    tagged mass still in the domain / tagged mass at t = 0
  f_dist(t)    the same with the weight m_p * shape(y_0)**2, shape(y) =
               exp(-((y - y_c)/y_w)**2) of convert_2d_to_3d.f90, i.e. the share
               of the initial disturbance energy whose fluid is still inside
  c_bub(t)     tagged mass inside the bubble / total mass of the bubble, the
               mean of c over the bubble (needs the flow snapshot at that time)
  f_dist_bub(t)  disturbance-weighted tagged fluid inside the bubble / the whole
               disturbance weight at t = 0 (share of the initial disturbance
               energy whose fluid is in the bubble)

The bubble is the region <u>_z(x, y) < --u-cut (0: reversed flow) of the flow
snapshot Q#####.vtr nearest in time; a particle is inside when <u>_z,
bilinearly interpolated to its (x, y), is below the cut. The wall row j = 1
carries no particles and is left out of the bubble mass. Only the rows
y <= max(seed y) are read from each 3.6 GB snapshot (np.memmap).

Usage
-----
The file needs only numpy (and matplotlib for the figure). Copy it into the
solver output directory (Q#####.vtr there, particle files in data_trace/) and
run it there::

    cd data
    python tracer_residual.py                     # all defaults
    python tracer_residual.py --dt-flow 1e-4      # flow snapshot interval nt*dt [s]
    python tracer_residual.py --yc 0.1e-3 --yw 0.03e-3 --u-cut 0 --out res_eps1e-6
    python tracer_residual.py -h                  # this text and the defaults

or from the repository root::

    python -m ouxsbli.analysis.tracer_residual path/to/data

Arguments:

  data        output directory with Q#####.vtr and data_trace/P######.dat
              (default: the current directory)
  --dt-flow   time between flow snapshots, nt*dt of mod_globals.f90 [s]
              (default 25000*4e-9 = 1e-4); Q#####.vtr is matched to a particle
              file when their times agree to 1 % of this
  --yc, --yw  y_center, y_width of the disturbance in convert_2d_to_3d.f90 [m]
              (default 0.1e-3, 0.03e-3); only f_dist and f_dist_bub use them
  --u-cut     the bubble is <u>_z below this [m/s] (default 0: reversed flow)
  --out       prefix of the output <out>.dat and <out>.png, written in the
              current directory (default tracer_residual)

Needed: data_trace/P000000.dat and Q00000.vtr (t = 0); without the later
Q#####.vtr the bubble columns are NaN. Columns of <out>.dat:
t n_alive f_mass f_dist c_bub f_dist_bub bub_top.
"""

import argparse
import glob
import os
import re

import numpy as np

NG = 3  # ghost planes on each side in z
_DTYPE = {"Float32": np.dtype("<f4"), "Float64": np.dtype("<f8")}


# _header/_arrays are copies of ouxsbli/analysis/vtr_raw.py (see there for the
# file layout), so that this script also runs on its own inside data/
def _header(path, limit=65536):
    with open(path, "rb") as f:
        head = f.read(limit)
    marker = head.find(b"<AppendedData")
    if marker < 0:
        raise ValueError(f"{path}: no <AppendedData> block in the first {limit} bytes")
    start = head.find(b"_", marker)
    if start < 0:
        raise ValueError(f"{path}: no '_' payload marker")
    return head[:marker].decode("ascii", "replace"), start + 1


def _arrays(xml):
    """[(name, dtype, ncomp, offset)] in file order; coordinates have name None."""
    out = []
    for tag in re.findall(r"<DataArray\b[^>]*/>", xml):
        t = re.search(r'type="([^"]+)"', tag)
        o = re.search(r'offset="\s*(\d+)"', tag)
        if not t or not o:
            continue
        name = re.search(r'Name="([^"]+)"', tag)
        ncomp = re.search(r'NumberOfComponents="\s*(\d+)"', tag)
        out.append((name.group(1) if name else None,
                    _DTYPE[t.group(1)],
                    int(ncomp.group(1)) if ncomp else 1,
                    int(o.group(1))))
    return out


def read_particles(path):
    """Fields of one data/data_trace/P######.dat (see the header of tracer.f90)."""
    with open(path, "rb") as f:
        n = int(np.fromfile(f, "<i4", 1)[0])
        step = int(np.fromfile(f, "<i8", 1)[0])
        t = float(np.fromfile(f, "<f8", 1)[0])
        x, y, z = (np.fromfile(f, "<f4", n) for _ in range(3))
        status = np.fromfile(f, "<i4", n)
        t_exit = np.fromfile(f, "<f4", n)
    return dict(n=n, step=step, t=t, x=x, y=y, z=z, status=status, t_exit=t_exit)


def _open_vtr(path):
    xml, base = _header(path)
    ext = re.search(r'WholeExtent="\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', xml)
    i0, i1, j0, j1, k0, k1 = (int(v) for v in ext.groups())
    shape = (k1 - k0 + 1, j1 - j0 + 1, i1 - i0 + 1)
    coords, fields = [], {}
    for name, dt, ncomp, off in _arrays(xml):
        if name is None:
            coords.append(np.fromfile(path, dtype=dt, count=shape[2 - len(coords)],
                                      offset=base + off + 4))
        else:
            fields[name] = np.memmap(path, dtype=dt, mode="r", offset=base + off + 4,
                                     shape=shape + ((ncomp,) if ncomp > 1 else ()))
    return coords, fields


def span_mean(path, jmax):
    """x, y[:jmax] and the spanwise means of rho and u over the unique planes."""
    (x, y, _), f = _open_vtr(path)
    nz = f["rho"].shape[0]
    rho = np.asarray(f["rho"][NG:nz - NG, :jmax, :]).mean(axis=0, dtype=np.float64)
    u = np.asarray(f["velocity"][NG:nz - NG, :jmax, :, 0]).mean(axis=0, dtype=np.float64)
    return x.astype(np.float64), y[:jmax].astype(np.float64), rho, u


def node_widths(a):
    """control-volume widths of the nodes a (half cells at the ends)."""
    d = np.diff(a)
    return 0.5 * (np.r_[d[:1], d] + np.r_[d, d[-1:]])


def bilinear(x, y, field, xs, ys):
    """field(y, x) interpolated to the points (xs, ys), clamped to the grid."""
    i = np.clip(np.searchsorted(x, xs) - 1, 0, x.size - 2)
    j = np.clip(np.searchsorted(y, ys) - 1, 0, y.size - 2)
    fx = np.clip((xs - x[i]) / (x[i + 1] - x[i]), 0.0, 1.0)
    fy = np.clip((ys - y[j]) / (y[j + 1] - y[j]), 0.0, 1.0)
    return ((1 - fy) * ((1 - fx) * field[j, i] + fx * field[j, i + 1])
            + fy * ((1 - fx) * field[j + 1, i] + fx * field[j + 1, i + 1]))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 epilog=__doc__[__doc__.index("Usage\n"):],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data", nargs="?", default=".",
                    help="solver output directory with Q#####.vtr and data_trace/P######.dat"
                         " (default: the current directory)")
    ap.add_argument("--dt-flow", type=float, default=25000 * 4.0e-9,
                    help="time between flow snapshots [s] (nt*dt of mod_globals.f90)")
    ap.add_argument("--yc", type=float, default=0.1e-3, help="y_center of the disturbance [m]")
    ap.add_argument("--yw", type=float, default=0.03e-3, help="y_width of the disturbance [m]")
    ap.add_argument("--u-cut", type=float, default=0.0, help="bubble = <u>_z below this [m/s]")
    ap.add_argument("--out", default="tracer_residual", help="prefix of the .dat/.png output")
    a = ap.parse_args(argv)

    pfiles = sorted(glob.glob(os.path.join(a.data, "data_trace", "P[0-9]*.dat")))
    if not pfiles:
        raise SystemExit(f"no P######.dat in {os.path.join(a.data, 'data_trace')}")
    p0 = read_particles(pfiles[0])
    if p0["step"] != 0:
        raise SystemExit(f"{pfiles[0]} is at step {p0['step']}; the seeding file P000000.dat is needed")

    # particle masses from the t = 0 flow snapshot (seeds sit exactly on nodes)
    q0 = os.path.join(a.data, "Q00000.vtr")
    (xg, yg, zg), _ = _open_vtr(q0)
    jmax = int(np.searchsorted(yg, p0["y"].max())) + 2
    x, y, rho0, _ = span_mean(q0, jmax)
    i0 = np.searchsorted(x, p0["x"]); j0 = np.searchsorted(y, p0["y"])
    si = np.diff(np.unique(i0)).min() if np.unique(i0).size > 1 else 1
    zs = np.unique(p0["z"])
    dz = float(zg[NG + 1] - zg[NG])
    sk = int(round((zs[1] - zs[0]) / dz)) if zs.size > 1 else 1
    lz = (zg.size - 2 * NG) * dz
    # z weight lz / (number of seeded planes), not sk * dz: the nz-6 unique planes need
    # not be a multiple of sk (134 planes, sk = 4 -> 34 * 4 = 136, +1.5 % in c_bub)
    m_p = rho0[j0, i0] * node_widths(x)[i0] * si * node_widths(y)[j0] * lz / zs.size
    w_p = m_p * np.exp(-((p0["y"] - a.yc) / a.yw) ** 2) ** 2
    print(f"{p0['n']} particles, strides si={si} sk={sk}, seeded y <= {p0['y'].max()*1e3:.3f} mm,"
          f" reading rows j < {jmax} of the flow snapshots")

    area = node_widths(x)[None, :] * node_widths(y)[:, None]
    rows = []
    for pf in pfiles:
        p = read_particles(pf)
        alive = p["status"] == 0
        row = dict(t=p["t"], n_alive=int(alive.sum()),
                   f_mass=m_p[alive].sum() / m_p.sum(), f_dist=w_p[alive].sum() / w_p.sum(),
                   c_bub=np.nan, f_dist_bub=np.nan, bub_top=np.nan)
        m = int(round(p["t"] / a.dt_flow))
        qf = os.path.join(a.data, f"Q{m:05d}.vtr")
        if abs(p["t"] - m * a.dt_flow) < 0.01 * a.dt_flow and os.path.exists(qf):
            _, _, rho, u = span_mean(qf, jmax)
            bub = u < a.u_cut
            bub[0] = False
            inside = alive & (bilinear(x, y, u, p["x"].astype(np.float64), p["y"].astype(np.float64)) < a.u_cut)
            mass_bub = lz * (rho * area)[bub].sum()
            if bub.any():
                row["c_bub"] = m_p[inside].sum() / mass_bub
                row["bub_top"] = y[np.nonzero(bub.any(axis=1))[0].max()]
                if bub[-3:].any():
                    print(f"WARNING t={p['t']*1e3:.3f} ms: the bubble reaches the top of the seeded"
                          " region, c_bub is underestimated")
            row["f_dist_bub"] = w_p[inside].sum() / w_p.sum()
        rows.append(row)
        print(f"t={row['t']*1e3:8.4f} ms  alive {row['n_alive']:9d}  f_mass {row['f_mass']:.4e}"
              f"  f_dist {row['f_dist']:.4e}  c_bub {row['c_bub']:.4e}  f_dist_bub {row['f_dist_bub']:.4e}")

    keys = ("t", "n_alive", "f_mass", "f_dist", "c_bub", "f_dist_bub", "bub_top")
    np.savetxt(a.out + ".dat", np.array([[r[k] for k in keys] for r in rows]),
               header=" ".join(keys) + f"   (bubble: <u>_z < {a.u_cut:g} m/s; y_c={a.yc:g} y_w={a.yw:g} m)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = np.array([r["t"] for r in rows]) * 1e3
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for k, lab in (("f_mass", "tagged mass in domain"), ("f_dist", "disturbance-weighted, domain"),
                   ("c_bub", "tagged mass fraction of bubble"), ("f_dist_bub", "disturbance-weighted, bubble")):
        v = np.array([r[k] for r in rows])
        ok = np.isfinite(v) & (v > 0)
        ax.semilogy(t[ok], v[ok], "-o", ms=3, label=lab)
    ax.axhline(1e-2, color="gray", lw=0.8, ls="--")
    ax.set_xlabel("t [ms]")
    ax.set_ylabel("fraction remaining")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(a.out + ".png", dpi=150)
    print("wrote", a.out + ".dat", a.out + ".png")


if __name__ == "__main__":
    main()
