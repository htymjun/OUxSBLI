"""Shared pytest fixtures and markers for OUxSBLI integration tests."""
import pathlib
import pytest


REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "integration: marks tests that build and run the CUDA Fortran solver "
        "(requires NVIDIA HPC SDK and a CUDA-capable GPU)",
    )
    config.addinivalue_line(
        "markers",
        "slow: long-running integration case (~10-30 min); "
        'deselect with -m "not slow"',
    )


@pytest.fixture(scope="session")
def repo_root():
    return REPO_ROOT


def assert_close_relative(computed, reference, rtol, label=""):
    """Assert that computed is within rtol of reference (relative error)."""
    rel_err = abs(computed - reference) / abs(reference)
    assert rel_err < rtol, (
        f"{label + ': ' if label else ''}"
        f"rel_err={rel_err:.3%}, computed={computed:.4g}, ref={reference:.4g}"
    )

