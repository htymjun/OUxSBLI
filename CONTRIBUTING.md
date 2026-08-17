# Contributing to OUxSBLI

Thank you for contributing to OUxSBLI.

OUxSBLI is a GPU-accelerated CFD solver for compressible flows. Contributions
should prioritize numerical correctness, reproducibility, maintainability,
and performance.

## Where to start

Before opening a new Issue or Discussion, search existing topics.

- **Bug Report** — incorrect behavior, crashes, build failures, or unexpected numerical results.
- **Feature Request** — a concrete feature or capability that should be added.
- **Validation Issue** — numerical validation problems or new validation cases.
- **Performance Ideas** — optimization ideas, including speculative ideas that are not yet ready for implementation.
- **Q&A / General Discussions** — questions and broader technical discussions.

As a rule of thumb:

> Issues describe work that may need to be done. Discussions are for ideas and conversations that are not yet committed work.

## Development workflow

Create a focused branch from `stable`:

```text
stable
  |
  +-- feature/<name>
  +-- fix/<name>
  +-- perf/<name>
  +-- validation/<name>
  +-- docs/<name>
```

Examples:

```text
feature/waLE-model
fix/slau-boundary-condition
perf/flux-kernel
validation/blasius
docs/installation
```

Open a Pull Request against `stable` when the change is ready for review.

## Repository-specific development

OUxSBLI uses CUDA Fortran, MPI, CMake, and fypp.

Before changing solver code, read the relevant `README.md` and `CLAUDE.md` documentation for the
subtree you are modifying.

Important rules include:

- Edit `.f90.fypp` templates rather than generated `.f90` files in `build/`.
- Modify a case's `config.fypp` for compile-time numerical configuration.
- Verify generated Fortran when changing fypp conditions.
- Run relevant Python tests and numerical validation when applicable.
- Do not commit generated build artifacts.

## Testing and validation

For the Python test suite:

```bash
pytest ouxsbli/tests/
```

For solver cases, follow the build and run instructions in the corresponding
solver documentation.

Numerical changes should be validated against an analytical solution,
reference dataset, experiment, or established benchmark whenever practical.

## Performance changes

Performance claims should be supported by measurements whenever possible.

Include:

- GPU model
- grid size
- MPI rank count
- compiler / NVIDIA HPC SDK version
- CUDA version
- runtime before and after
- relevant Nsight Systems / Nsight Compute observations

A speedup should not be considered sufficient evidence by itself if the
optimization changes numerical behavior.

## Pull Requests

Keep each PR focused on one logical change.

A PR should explain:

1. What changed.
2. Why it changed.
3. How it was implemented.
4. How it was validated.
5. What numerical or performance impact it has.

Use the Pull Request template.

## AI-assisted development

AI coding assistants may be used for development.

Contributors remain responsible for:

- correctness
- numerical validity
- performance claims
- scientific interpretation
- reproducibility
- licensing

AI-generated code must be reviewed and validated by the contributor.

In particular, passing compilation or unit tests does not establish that a
CFD implementation is physically or numerically correct.

## Commit messages

Prefer concise, descriptive commit messages.

Examples:

```text
fix: correct wall boundary condition
feat: add WALE SGS model
perf: fuse HRSLAU2 flux kernels
test: add Blasius validation
docs: update GPU build instructions
```

## Review philosophy

Reviews should focus on:

- correctness
- numerical stability
- physical validity
- reproducibility
- performance
- maintainability

Constructive technical discussion is encouraged.
