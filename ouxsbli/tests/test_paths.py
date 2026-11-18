"""
OUxSBLI from any directory: source-tree discovery (`ouxsbli.paths`), case
resolution, the `ouxsbli` command line, and a case directory outside the tree
that builds on its own (`cmake -S <case> -B <case>/build`, OUXSBLI_ROOT).
"""
import os
import pathlib
import re
import shutil
import subprocess
import sys

import pytest

from ouxsbli import Case
from ouxsbli.cli import main
from ouxsbli.paths import repo_root, resolve_case, template_cases
from .conftest import REPO_ROOT


def test_root_from_package_and_env(monkeypatch, tmp_path):
    monkeypatch.delenv("OUXSBLI_ROOT", raising=False)
    assert repo_root() == REPO_ROOT
    monkeypatch.setenv("OUXSBLI_ROOT", str(REPO_ROOT))
    assert repo_root() == REPO_ROOT
    monkeypatch.setenv("OUXSBLI_ROOT", str(tmp_path))  # not a checkout
    with pytest.raises(RuntimeError, match="OUXSBLI_ROOT"):
        repo_root()


def test_resolve_case(monkeypatch, tmp_path):
    monkeypatch.delenv("OUXSBLI_ROOT", raising=False)
    monkeypatch.chdir(tmp_path)  # nothing here: tree-relative and absolute names still resolve
    evc = REPO_ROOT / "2D_solver" / "EVC"
    assert resolve_case("2D_solver/EVC") == evc
    assert resolve_case(str(evc)) == evc
    # a case of your own, named relative to the current directory
    shutil.copytree(evc, tmp_path / "mine" / "EVC")
    assert resolve_case("mine/EVC") == tmp_path / "mine" / "EVC"
    with pytest.raises(FileNotFoundError):
        resolve_case("mine/NOPE")
    assert "2D_solver/EVC" in template_cases() and "3D_solver/NSTGV" in template_cases()


def test_cli_init_patches_a_copy(monkeypatch, tmp_path, capsys):
    monkeypatch.delenv("OUXSBLI_ROOT", raising=False)
    monkeypatch.delenv("OUXSBLI_CONFIG", raising=False)
    monkeypatch.chdir(tmp_path)
    assert main(["init", "2D_solver/EVC", "evc", "nx=66", "accuracy=4", "prec=df"]) == 0
    assert "integer, parameter :: nx = 66" in (tmp_path / "evc" / "mod_globals.f90").read_text()
    cfg = (tmp_path / "evc" / "config.fypp").read_text()
    assert re.search(r"#:set\s+ORDER\s*=\s*4", cfg) and re.search(r"#:set\s+PREC\s*=\s*'DF'", cfg)
    assert (tmp_path / "evc" / "CMakeLists.txt").exists() and (tmp_path / "evc" / "set.f90").exists()
    assert main(["root"]) == 0 and str(REPO_ROOT) in capsys.readouterr().out
    assert main(["build", "no_such_dir"]) == 1


def test_installed_package_finds_the_tree_from_elsewhere(tmp_path):
    """`import ouxsbli` in a fresh interpreter started outside the checkout."""
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "OUXSBLI_ROOT")}
    r = subprocess.run([sys.executable, "-c", "import ouxsbli; print(ouxsbli.repo_root())"],
                       cwd=tmp_path, env=env, capture_output=True, text=True)
    if r.returncode != 0 and "No module named" in r.stderr:
        pytest.skip("ouxsbli is not installed (pip install -e .)")
    assert r.returncode == 0, r.stderr
    root = pathlib.Path(r.stdout.strip())
    assert (root / "2D_solver").is_dir() and (root / "src").is_dir()


@pytest.mark.integration
def test_case_outside_the_tree_builds_and_runs(monkeypatch, tmp_path):
    monkeypatch.delenv("OUXSBLI_CONFIG", raising=False)
    monkeypatch.chdir(tmp_path)
    case = Case("2D_solver/EVC", "evc", nx=66, ny=66, nt=5, np=1)
    case._setup()
    # plain CMake in the copied directory, the tree named by the environment only
    env = case._build_env()
    env["OUXSBLI_ROOT"] = str(REPO_ROOT)
    for cmd in (["cmake", "-S", ".", "-B", "build"], ["cmake", "--build", "build"]):
        r = subprocess.run(cmd, cwd=tmp_path / "evc", env=env, capture_output=True, text=True)
        assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    assert (tmp_path / "evc" / "build" / "a.out").exists()
    Case.from_dir("evc").run(nranks=2)
    assert list((tmp_path / "evc" / "data").glob("*.vtr"))
