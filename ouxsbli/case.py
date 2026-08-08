"""
Case — lifecycle manager for a single OUxSBLI simulation run.

Usage::

    from ouxsbli import Case

    case = Case(
        source  = "3D_solver/NSTGV",
        workdir = "/tmp/run_01",
        Re      = 800.0,
        nx      = 65,
        scheme  = "KEEP",
        visc    = "NS",
        accuracy = 6,
    )
    case.build()
    case.run(nranks=2)
"""
import os
import shutil
import subprocess
import pathlib
from typing import Any
from .patcher import patch, _is_macro_line
from .ic import ic_to_params, bc_to_params

# ---------------------------------------------------------------------------
# Alias table: friendly Case() kwarg → fypp macro name in config.fypp
# ---------------------------------------------------------------------------
_ALIAS = {
    "scheme":     "SCHEME",
    "visc":       "VISC",
    "accuracy":   "ORDER",
    "visc_order": "VISC_ORDER",
    "tvd":        "TVD",
    "slau":       "SLAU_VARIANT",
    "rescale":    "RESCALE",
    "recal":      "RESTART",
    "rk":         "RK",
    "gpumpi":     "GPUMPI",
    "bc_x":       "BC_X",
    "bc_y":       "BC_Y",
    "commz":      "COMMZ",
    "output_precision": "OUTPUT_PRECISION",
}

# ---------------------------------------------------------------------------
# Value normalisation for config.fypp string parameters
# fypp comparisons are case-sensitive: 'Euler' ≠ 'EULER', 'none' ≠ 'NONE'.
# ---------------------------------------------------------------------------
_VALUE_NORMALIZE: dict[str, dict[str, str]] = {
    "SCHEME":       {"keep": "KEEP", "slau": "SLAU", "hybrid": "Hybrid"},
    "VISC":         {"euler": "Euler", "ns": "NS", "les": "LES"},
    "TVD":          {"none": "none", "tvd": "tvd", "hybrid": "hybrid",
                     "minmod": "tvd", "muscl4": "hybrid"},  # backward-compat aliases
    "SLAU_VARIANT": {"slau": "SLAU", "hrslau2": "HRSLAU2"},
}

# ---------------------------------------------------------------------------
# RK stage-count normalisation  (human-friendly strings → int)
# ---------------------------------------------------------------------------
_RK_NORMALIZE: dict[str, int] = {
    "tvd_rk3": 3, "rk3": 3,
    "rk4": 4, "classical_rk4": 4,
}


class Case:
    """Manage one simulation run (copy source -> patch config.fypp + mod_globals.f90 -> cmake/make -> run)."""

    def __init__(self, source: str, workdir: str, **params: Any) -> None:
        """
        Parameters
        ----------
        source:
            Path to the original case directory (e.g. ``"3D_solver/NSTGV"``).
            This directory is never modified.
        workdir:
            Path for the new working directory that will be created.
        **params:
            Parameter overrides.  Recognised aliases (e.g. ``scheme``,
            ``visc``, ``accuracy``) are expanded to the fypp macro names used
            in ``config.fypp``.  All other parameters (``nx``, ``ny``, ``Re``,
            etc.) are forwarded to ``mod_globals.f90``.
        """
        self.repo_root = pathlib.Path(__file__).resolve().parent.parent
        self.source  = (self.repo_root / source).resolve()
        self.workdir = pathlib.Path(workdir).resolve()
        self.build_dir = self.workdir / "build"
        # Expand friendly aliases; do NOT blanket-uppercase values —
        # config.fypp comparisons are case-sensitive (e.g. 'Euler' ≠ 'EULER').
        # Expand IC/BC specs into flat parameter dicts before alias expansion.
        raw = dict(params)
        ic_spec = raw.pop("ic", None)
        bc_spec = raw.pop("bc", None)
        if ic_spec is not None:
            raw.update(ic_to_params(ic_spec))
        if bc_spec is not None:
            raw.update(bc_to_params(bc_spec))

        self.params: dict[str, Any] = {
            _ALIAS.get(k.lower(), k): v
            for k, v in raw.items()
        }
        self._built = False

    # ------------------------------------------------------------------
    # Setup (copy + patch)
    # ------------------------------------------------------------------

    def _setup(self) -> None:
        """Copy source directory to workdir, then patch config.fypp and mod_globals.f90."""
        if self.workdir.exists():
            shutil.rmtree(self.workdir)

        shutil.copytree(
            self.source,
            self.workdir,
            # data/recal are excluded so snapshots from a previous run of the
            # source case can't be mistaken for this run's output.
            ignore=shutil.ignore_patterns(
                "build", "CMakeCache.txt", "CMakeFiles", "*.cmake",
                "data", "data_*", "recal", "nohup.out",
            ),
        )

        config_path  = self.workdir / "config.fypp"
        globals_path = self.workdir / "mod_globals.f90"

        if not config_path.exists():
            raise FileNotFoundError(f"config.fypp not found in {self.source}")
        if not globals_path.exists():
            raise FileNotFoundError(f"mod_globals.f90 not found in {self.source}")

        config_text  = config_path.read_text()
        globals_text = globals_path.read_text()
        config_lines = config_text.splitlines()

        config_params:  dict[str, Any] = {}
        globals_params: dict[str, Any] = {}

        for k, v in self.params.items():
            # Normalise RK string → integer stage count
            if k == "RK" and isinstance(v, str):
                v = _RK_NORMALIZE.get(v.lower(), v)

            if any(_is_macro_line(line, k) for line in config_lines):
                # param lives in config.fypp — apply fypp-correct casing for strings
                if isinstance(v, str) and k in _VALUE_NORMALIZE:
                    v = _VALUE_NORMALIZE[k].get(v.lower(), v)
                config_params[k] = v
            else:
                # param is not in config.fypp → route to mod_globals.f90
                globals_params[k] = v

        config_path.write_text(patch(config_text, config_params))
        globals_path.write_text(patch(globals_text, globals_params))

    # ------------------------------------------------------------------
    # Build
    # ------------------------------------------------------------------

    def build(self) -> None:
        """Apply patches, generate build system via CMake, and compile."""
        self._setup()

        self.build_dir.mkdir(exist_ok=True)
        # Set the CMake src dir to "original case dir".
        # This ensures that relative paths to ../CMakeLists.txt and src/ work correctly.
        cmake_cmd = [
            "cmake",
            str(self.source),  # This line is important
            f"-DCASE_DIR={self.workdir}"
        ]

        self._run_cmd(cmake_cmd, cwd=self.build_dir)
        self._run_cmd(["make"], cwd=self.build_dir)

        self._built = True

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(self, nranks: int = 2) -> None:
        """Launch the simulation with mpirun.

        Parameters
        ----------
        nranks:
            Number of MPI ranks (default 2).
        """
        if not self._built:
            raise RuntimeError("Must call build() before run()")

        data_dir = self.workdir / "data"
        data_dir.mkdir(exist_ok=True)

        # NOTE: Depending on CMakeLists.txt, the executable may be placed in build_dir
        # Output artifacts remain cleanly in workdir.
        exe_path = self.build_dir / "a.out"
        if not exe_path.exists():
            # Fallback path if CMake installs the executable to the workdir root instead
            exe_path = self.workdir / "a.out"

        self._run_cmd(["mpirun", "-n", str(nranks), str(exe_path)], cwd=self.workdir)

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _build_env(self) -> dict:
        """Return os.environ copy with NVIDIA HPC SDK paths prepended."""
        env = os.environ.copy()
        nvfortran = shutil.which("nvfortran")
        if nvfortran is None:
            hpcsdk_root = pathlib.Path("/opt/nvidia/hpc_sdk")
            candidates = sorted(
                hpcsdk_root.glob("Linux_x86_64/*/compilers/bin/nvfortran"),
                reverse=True,
            )
            if candidates:
                nvfortran = str(candidates[0])
        if nvfortran is None:
            raise RuntimeError(
                "NVIDIA HPC SDK not found. Install it from "
                "https://developer.nvidia.com/hpc-sdk or add nvfortran to PATH."
            )
        compiler_bin = pathlib.Path(nvfortran).parent
        sdk_root = compiler_bin.parent.parent
        extra = [str(compiler_bin)]
        for mpi_subpath in ("comm_libs/hpcx/bin", "comm_libs/openmpi4/bin"):
            mpi_bin = sdk_root / mpi_subpath
            if mpi_bin.exists():
                extra.append(str(mpi_bin))
                env.setdefault("FC", str(mpi_bin / "mpif90"))
                break
        env["PATH"] = ":".join(extra) + ":" + env.get("PATH", "")
        return env

    def _run_cmd(self, cmd: list[str], cwd: pathlib.Path) -> None:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=True,
            text=True,
            env=self._build_env(),
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"Command {cmd} failed (exit {result.returncode}):\n"
                f"STDOUT:\n{result.stdout}\n"
                f"STDERR:\n{result.stderr}"
            )

    def __repr__(self) -> str:
        return (
            f"Case(source={self.source.name!r}, "
            f"workdir={str(self.workdir)!r}, "
            f"params={self.params})"
        )
