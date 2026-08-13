"""Dependency-free reader for the solver's .vtr snapshots (numpy only).

``ouxsbli/tests/utils/vtk_reader.py`` needs the ``vtk`` package, which is often
absent on a compute cluster -- exactly where you want to reduce 13 GB of
snapshots down to a megabyte of moments before copying anything back. The
files are plain VTK XML RectilinearGrid with raw appended binary
(``src/print.f90.fypp``), so parsing them directly needs nothing but numpy.

Layout written by ``print_xml``: a text header, then ``<AppendedData
encoding="raw">``, then ``_`` and, for each array in turn, a 4-byte little
endian count followed by the data. Note the count includes its own 4 bytes,
which is not VTK's usual convention -- so seek by each DataArray's declared
``offset`` attribute and skip 4, rather than trusting the counts to chain.

Point data is ``rho``, ``p`` and a 3-component ``velocity``, in C order
(k, j, i), i.e. shaped ``(nz, ny, nx)`` to match ``vtk_reader.getQ``.
"""

import re

import numpy as np

_DTYPE = {"Float32": np.dtype("<f4"), "Float64": np.dtype("<f8")}


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


def read_vtr(path):
    """Return ``(x, y, z, {name: array})`` with arrays shaped (nz, ny, nx[, ncomp])."""
    xml, base = _header(path)
    ext = re.search(r'WholeExtent="\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', xml)
    if not ext:
        raise ValueError(f"{path}: no WholeExtent")
    i0, i1, j0, j1, k0, k1 = (int(v) for v in ext.groups())
    nx, ny, nz = i1 - i0 + 1, j1 - j0 + 1, k1 - k0 + 1

    arrays = _arrays(xml)
    coord_n = (nx, ny, nz)
    coords, fields, ncoord = [], {}, 0
    with open(path, "rb") as f:
        for name, dt, ncomp, off in arrays:
            f.seek(base + off + 4)  # +4 skips this block's own count
            if name is None:
                n = coord_n[ncoord]
                coords.append(np.frombuffer(f.read(n * dt.itemsize), dtype=dt))
                ncoord += 1
            else:
                n = nx * ny * nz * ncomp
                a = np.frombuffer(f.read(n * dt.itemsize), dtype=dt)
                fields[name] = a.reshape((nz, ny, nx, ncomp) if ncomp > 1 else (nz, ny, nx))
    if ncoord != 3:
        raise ValueError(f"{path}: expected 3 coordinate arrays, found {ncoord}")
    return coords[0], coords[1], coords[2], fields


def getGrid(path):
    """``(ni, nj, nk, x, y, z)`` -- signature-compatible with vtk_reader.getGrid."""
    x, y, z, _ = read_vtr(path)
    return len(x), len(y), len(z), x, y, z


def getQ(path, Nx=None, Ny=None, Nz=None, reader=None):
    """``(rho, u, v, w, p)`` shaped (nz, ny, nx) -- as vtk_reader.getQ.

    The size arguments are accepted for signature compatibility and checked
    against the file rather than trusted.
    """
    _, _, _, f = read_vtr(path)
    rho, p, vel = f["rho"], f["p"], f["velocity"]
    if Nx is not None and rho.shape != (Nz, Ny, Nx):
        raise ValueError(f"{path}: shape {rho.shape} != requested {(Nz, Ny, Nx)}")
    return rho, vel[..., 0], vel[..., 1], vel[..., 2], p
