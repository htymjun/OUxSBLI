"""
Locate the OUxSBLI source tree and case directories.

The solver sources (``2D_solver/``, ``3D_solver/``, ``3D_solver_curv/``, ``src/``,
``extern/``) are not part of the installed Python package, so they are found at
run time:

1. ``OUXSBLI_ROOT`` if it is set (selects a checkout explicitly, e.g. one of
   several worktrees);
2. otherwise the directory that contains this package (an editable install,
   ``pip install -e .``, or a script run from the checkout).
"""
import os
import pathlib

_MARKERS = ("2D_solver", "3D_solver", "src")
SOLVERS = ("2D_solver", "3D_solver", "3D_solver_curv")


def _is_root(path: pathlib.Path) -> bool:
    return all((path / m).is_dir() for m in _MARKERS)


def root_origin() -> str:
    """How :func:`repo_root` finds the tree: ``'OUXSBLI_ROOT'`` or ``'package'``."""
    return "OUXSBLI_ROOT" if os.environ.get("OUXSBLI_ROOT") else "package"


def repo_root() -> pathlib.Path:
    """Return the root of the OUxSBLI source tree."""
    env = os.environ.get("OUXSBLI_ROOT")
    if env:
        path = pathlib.Path(env).expanduser().resolve()
        if not _is_root(path):
            raise RuntimeError(
                f"OUXSBLI_ROOT={env} is not an OUxSBLI checkout "
                f"(expected the directories {', '.join(_MARKERS)} in it)."
            )
        return path
    path = pathlib.Path(__file__).resolve().parent.parent
    if _is_root(path):
        return path
    raise RuntimeError(
        f"The OUxSBLI solver sources were not found next to the Python package ({path}).\n"
        "Install the package in editable mode from a checkout (pip install -e <OUxSBLI>) "
        "or set OUXSBLI_ROOT to the checkout."
    )


def resolve_case(source) -> pathlib.Path:
    """Return the case directory named by *source*.

    An absolute path is used as is. A relative path is looked up in the source
    tree first (``"2D_solver/EVC"``) and then relative to the current directory
    (a case directory of your own, outside the tree).
    """
    path = pathlib.Path(source).expanduser()
    if path.is_absolute():
        candidates = [path]
    else:
        candidates = [repo_root() / path, pathlib.Path.cwd() / path]
    for c in candidates:
        if (c / "config.fypp").is_file():
            return c.resolve()
    tried = ", ".join(str(c) for c in candidates)
    raise FileNotFoundError(f"No case directory (with a config.fypp) for {source!r}; tried {tried}")


def template_cases() -> list[str]:
    """Case directories shipped in the source tree, as paths relative to the root."""
    root = repo_root()
    return sorted(
        str(p.parent.relative_to(root))
        for s in SOLVERS
        for p in (root / s).glob("*/config.fypp")
    )
