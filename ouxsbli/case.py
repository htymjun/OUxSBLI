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

``source`` is a case directory of the source tree (``"2D_solver/EVC"``), an
absolute path, or a path relative to the current directory; the source tree is
found by :mod:`ouxsbli.paths` (``OUXSBLI_ROOT`` or the package location), so a
script can run from any directory. An existing case directory is used in place
with ``Case.from_dir(path)``.
"""
import os
import re
import shutil
import subprocess
import pathlib
import warnings
from typing import Any
from .patcher import patch, _is_macro_line
from .ic import ic_to_params, bc_to_params
from .paths import repo_root, resolve_case

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
    "recon":      "RECON",
    "weno_prec":  "WENO_PREC",
    "riemann":    "RIEMANN",
    "prec":       "PREC",
    "keep_prec":  "KEEP_PREC",
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
    "RECON":        {"muscl": "MUSCL", "weno": "WENO"},
    "WENO_PREC":    {"fp64": "FP64", "w32p64": "W32P64", "mixed": "MIXED"},
    "PREC":         {"fp64": "FP64", "df": "DF"},
    "KEEP_PREC":    {"fp64": "FP64", "mixed": "MIXED"},
    "RIEMANN":      {"slau": "SLAU", "lf": "LF"},
}

_CONFIG_KEYS = set(_ALIAS.values()) | set(_VALUE_NORMALIZE)

# ---------------------------------------------------------------------------
# RK stage-count normalisation  (human-friendly strings → int)
# ---------------------------------------------------------------------------
_RK_NORMALIZE: dict[str, int] = {
    "tvd_rk3": 3, "rk3": 3,
    "rk4": 4, "classical_rk4": 4,
}


# Oldest NVIDIA HPC SDK that builds OUxSBLI correctly (see docs/installation.md)
MIN_HPCSDK = (25, 5)


def _sdk_version(nvfortran) -> tuple:
    """(major, minor) of the HPC SDK that owns *nvfortran*, from its install path
    (.../<version>/compilers/bin/nvfortran); (0, 0) if the path has no version."""
    m = re.search(r"/(\d+)\.(\d+)/compilers/bin/nvfortran$", str(pathlib.Path(nvfortran).resolve()))
    return (int(m.group(1)), int(m.group(2))) if m else (0, 0)


class Case:
    """Manage one simulation run (copy source -> patch config.fypp + mod_globals.f90 -> cmake/make -> run)."""

    def __init__(self, source: str, workdir: str, verbose: bool = False, **params: Any) -> None:
        """
        Parameters
        ----------
        source:
            The original case directory: relative to the source tree
            (``"3D_solver/NSTGV"``), absolute, or relative to the current
            directory. This directory is never modified.
        workdir:
            Path for the new working directory that will be created
            (relative paths are taken from the current directory).
        verbose:
            Show the output of cmake, make and mpirun instead of capturing it.
        **params:
            Parameter overrides.  Recognised aliases (e.g. ``scheme``,
            ``visc``, ``accuracy``) are expanded to the fypp macro names used
            in ``config.fypp``.  All other parameters (``nx``, ``ny``, ``Re``,
            etc.) are forwarded to ``mod_globals.f90``.
        """
        self.repo_root = repo_root()
        self.source  = resolve_case(source)
        self.workdir = pathlib.Path(workdir).expanduser().resolve()
        self.build_dir = self.workdir / "build"
        self.verbose = verbose
        self._inplace = False
        # Expand friendly aliases; do NOT blanket-uppercase values —
        # config.fypp comparisons are case-sensitive (e.g. 'Euler' ≠ 'EULER').
        # Expand IC/BC specs into flat parameter dicts before alias expansion.
        raw = dict(params)
        # OUXSBLI_CONFIG="prec=df,keep_prec=mixed" adds overrides to every Case,
        # so an existing test runs unchanged in another precision mode
        for kv in filter(None, os.environ.get("OUXSBLI_CONFIG", "").split(",")):
            k, v = kv.split("=", 1)
            raw[k.strip()] = v.strip()
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

    @classmethod
    def from_dir(cls, path: str = ".", verbose: bool = False) -> "Case":
        """Use the existing case directory *path* in place (no copy, no patching)."""
        case = cls.__new__(cls)
        case.repo_root = repo_root()
        case.source = case.workdir = resolve_case(pathlib.Path(path).expanduser().resolve())
        case.build_dir = case.workdir / "build"
        case.verbose = verbose
        case._inplace = True
        case.params = {}
        case._built = (case.build_dir / "a.out").exists()
        return case

    # ------------------------------------------------------------------
    # Setup (copy + patch)
    # ------------------------------------------------------------------

    def _setup(self) -> None:
        """Copy source directory to workdir, then patch config.fypp and mod_globals.f90."""
        if self._inplace:
            return
        if self.workdir == self.source:
            raise ValueError(f"workdir is the source case directory ({self.source}); use Case.from_dir()")
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

            # A known fypp macro goes to config.fypp even when the case's file does not
            # set it yet (e.g. VISC_ORDER, PREC): patch() appends the missing #:set line.
            if k in _CONFIG_KEYS or any(_is_macro_line(line, k) for line in config_lines):
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
        # The working directory is a complete case (its own CMakeLists.txt, config.fypp,
        # mod_globals.f90, set.f90); the shared sources come from OUXSBLI_ROOT.
        cmake_cmd = [
            "cmake",
            "-S", str(self.workdir),
            "-B", str(self.build_dir),
            f"-DOUXSBLI_ROOT={self.repo_root}",
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
        """Return os.environ copy with NVIDIA HPC SDK paths prepended.

        OUxSBLI needs HPC SDK 25.5 or newer. If the nvfortran on PATH is older
        (or absent), the newest installed SDK under the same root (or
        /opt/nvidia/hpc_sdk) is used instead; without one, this raises.
        """
        env = os.environ.copy()
        on_path = shutil.which("nvfortran")
        roots = [pathlib.Path("/opt/nvidia/hpc_sdk/Linux_x86_64")]
        if on_path is not None:
            # <root>/<version>/compilers/bin/nvfortran
            roots.insert(0, pathlib.Path(on_path).resolve().parents[3])
        nvfortran = on_path if on_path and _sdk_version(on_path) >= MIN_HPCSDK else None
        if nvfortran is None:
            found = sorted(
                (c for r in dict.fromkeys(roots) for c in r.glob("*/compilers/bin/nvfortran")
                 if _sdk_version(c) >= MIN_HPCSDK),
                key=_sdk_version,
            )
            if found:
                nvfortran = str(found[-1])
                if on_path is not None:
                    warnings.warn(
                        f"nvfortran on PATH ({on_path}) is older than HPC SDK "
                        f"{MIN_HPCSDK[0]}.{MIN_HPCSDK[1]}; using {nvfortran}", stacklevel=3)
        if nvfortran is None:
            raise RuntimeError(
                f"NVIDIA HPC SDK {MIN_HPCSDK[0]}.{MIN_HPCSDK[1]} or newer not found"
                + (f" (nvfortran on PATH: {on_path})" if on_path else "") + ". Install it from "
                "https://developer.nvidia.com/hpc-sdk or add its nvfortran to PATH."
            )
        compiler_bin = pathlib.Path(nvfortran).parent
        sdk_root = compiler_bin.parent.parent
        extra = [str(compiler_bin)]
        for mpi_subpath in ("comm_libs/hpcx/bin", "comm_libs/openmpi4/bin"):
            mpi_bin = sdk_root / mpi_subpath
            if mpi_bin.exists():
                extra.append(str(mpi_bin))
                if on_path != nvfortran:
                    env["FC"] = str(mpi_bin / "mpif90")
                else:
                    env.setdefault("FC", str(mpi_bin / "mpif90"))
                break
        env["PATH"] = ":".join(extra) + ":" + env.get("PATH", "")
        return env

    def _run_cmd(self, cmd: list[str], cwd: pathlib.Path) -> None:
        if self.verbose:
            result = subprocess.run(cmd, cwd=cwd, env=self._build_env())
            if result.returncode != 0:
                raise RuntimeError(f"Command {cmd} failed (exit {result.returncode})")
            return
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
