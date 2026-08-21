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

Create a focused branch from `stable` for ordinary development:

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
The `stable` branch should remain the stable CUDA Fortran baseline.

CUDA C porting work is handled separately:

```text
cuda-c-stable
  |
  +-- cuda-c-<name>
```

All branches for the `Porting to CUDA C` milestone must use the
`cuda-c-*` naming pattern and open Pull Requests against `cuda-c-stable`, not
`stable`. The CUDA C branch family must preserve the existing CUDA Fortran
implementation while adding CUDA C paths alongside it. Duplicated build/test
surfaces for CUDA Fortran and CUDA C are acceptable during the migration.

For AI-assisted issue-driven development, use the OUxSBLI AI workflow in
`docs/ai_development_workflow.md`. In short: issues hold the scientific and
numerical decisions, AI assistants may propose and implement within the
approved scope, and the human maintainer reviews unresolved choices before
commit.

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

For issue-driven AI work:

- Keep human-owned scientific and numerical choices in the issue.
- Ask the AI for a plan before implementation when the change touches solver
  behavior, validation tolerances, performance claims, or public API behavior.
- Do not let the AI silently choose final boundary conditions, model constants,
  default schemes, validation tolerances, or benchmark acceptance criteria.
- For CUDA C migration work, split the milestone into a foundation issue and
  child issues for fused kernels, build-system changes, validation reports,
  translated kernel families, benchmarks, and documentation.
- Use `cuda-c-*` branches and target `cuda-c-stable` for CUDA C migration PRs;
  do not merge migration work directly into `stable`.
- Keep the existing CUDA Fortran implementation available on `cuda-c-stable`
  while CUDA C paths are introduced.
- CUDA Fortran performance PRs merged into `stable` may trigger a CUDA C
  carry-over issue when they carry a performance-style label, a `perf:` commit,
  a `perf:` PR title, or a `perf/*` branch name. The follow-up issue evaluates
  usefulness for `cuda-c-stable`; it does not imply automatic porting.
- The CUDA C milestone starts from the current CUDA Fortran `stable` baseline.
  The historical `CUDA_C_JAX_Python` branch is out of scope; reusing code or
  structure from it requires a decision recorded in an issue first.
- Include human-readable validation evidence in issues and PRs; pytest output
  alone is not sufficient for solver-behavior changes.
- Review the final diff and validation evidence before committing.

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
