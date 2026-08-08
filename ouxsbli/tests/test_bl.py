"""
BL: laminar flat-plate boundary layer (2D_solver/BL).

A uniform M=0.1 freestream develops a boundary layer over a flat plate whose
leading edge sits at x=0 (symmetry BC upstream of it, no-slip downstream).
At this Mach number compressibility is negligible, so both the skin-friction
distribution and the streamwise velocity profile must collapse onto the
incompressible Blasius solution.

The comparison assertions (Blasius Cf/profile match, quasi-steady guard,
freestream-not-accelerated guard) are shared with test_bl_3d.py -- see
utils/bl_common.py for the full rationale, including the acoustic-impedance
story behind why test_freestream_is_not_accelerated exists at all.

Production's own defaults (2D_solver/BL/mod_globals.f90) are now cheap enough
(nx=257, ny=49, Ly=15mm, dt=4e-8, endT=1e-2 -> ~2.5 min) that this test runs the
case unmodified apart from output cadence and precision.
"""
import pathlib

import pytest

from ouxsbli import Case
from .utils.vtk_reader import extract_number
from .utils.bl_common import (  # noqa: F401 -- imported for pytest collection
    test_grid_places_leading_edge_at_origin,
    test_flow_is_quasi_steady,
    test_freestream_is_not_accelerated,
    test_cf_matches_blasius,
    test_velocity_profile_matches_blasius,
)

pytestmark = pytest.mark.integration


@pytest.fixture(scope="module")
def bl_run():
    """Build and run the BL case once, returning the last two snapshots and k_ref=0."""
    workdir = "./tmp/bl"

    # nx/ny/Ly/dt/endT are left at their (now cheap) production defaults -- see
    # 2D_solver/BL/mod_globals.f90. Only the snapshot cadence and precision are
    # testing-specific: np=10 gives fewer, larger-precision outputs than
    # production's np=100, at the same dt this makes nt derive to 25000 (do not
    # pass nt explicitly).
    case = Case(
        source="2D_solver/BL",
        workdir=workdir,
        np=10,
        output_precision=8,
    )
    case.build()
    case.run(nranks=2)

    data_dir = pathlib.Path(workdir) / "data"
    snaps = sorted(data_dir.glob("Q*.vtr"), key=lambda p: extract_number(p.name))
    assert len(snaps) >= 2, f"expected at least 2 snapshots in {data_dir}, got {len(snaps)}"
    return snaps[-2], snaps[-1], 0
