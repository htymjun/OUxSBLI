# AI-assisted issue workflow

This workflow is for issue-driven development with AI assistants in OUxSBLI.
It is intentionally more gated than a small single-kernel research prototype:
OUxSBLI has multiple solvers, fypp-generated Fortran, GPU/MPI execution paths,
and numerical validation requirements that cannot be delegated entirely to an
assistant.

## Why OUxSBLI needs explicit gates

Issue-driven AI development can be lightweight when an issue maps to a compact
algorithmic change, a narrow benchmark, and a single validation surface.
OUxSBLI needs a stricter workflow because a change can cross several axes at
once:

- solver tree: `2D_solver/`, `3D_solver/`, `3D_solver_curv/`, or `ouxsbli/`
- generated code: `.f90.fypp` templates produce case-local build outputs
- physics: Euler, Navier-Stokes, LES, wall models, shocks, curvilinear grids
- numerics: scheme, order, limiter, boundary condition, time integration
- execution: GPU kernels, MPI decomposition, output precision, VTK I/O
- evidence: compilation, analytical validation, reference data, benchmarks

Therefore, the issue is not just a task ticket. It is the place where the
human owner records the scientific and numerical decisions that the AI must
treat as constraints.

## What is automated and what is not

Be explicit about this, because every gate below is otherwise easy to mistake
for something CI enforces.

The repository has exactly two GitHub Actions workflows:

- `.github/workflows/deploy-docs.yml` — publishes MkDocs on pushes to `stable`
- `.github/workflows/cuda-c-followup-from-perf.yml` — opens a CUDA C carry-over
  issue when a performance PR is merged into `stable`

There is **no build CI, no pytest CI, and no GPU runner**. Nothing in this
repository compiles the solver, runs a case, or checks a tolerance
automatically. Every build result, numerical comparison, and benchmark number in
an issue or PR is produced by hand on a developer machine and pasted in as
evidence. Consequently:

- validation claims must name the machine they came from
- every performance or numerical claim must record GPU model, HPC SDK version,
  MPI rank count, grid size, `dt`, and end time
- "CI is green" is never an argument in this project, because there is no CI to
  be green

## Issue states

Use these states as labels. The label is the source of truth; a heading inside
the issue body is a convenience, not the state.

| State | Label | Owner | Meaning |
|---|---|---|---|
| Triage | `triage` | Human | Decide whether this is a bug, feature, validation, performance idea, or discussion. |
| Needs spec | `needs-spec` | Human | Scientific/numerical choices are still open. AI may analyze, but must not implement a final design. |
| Ready for AI plan | `ready-for-ai-plan` | Human | The required decisions are recorded. AI may propose an implementation plan. |
| Plan review | `plan-review` | Human | AI has proposed a plan; human accepts, edits, or rejects it. |
| AI implementation | `ai-implementation` | AI | AI implements only the approved scope. |
| Needs human decision | `needs-human-decision` | Human | Implementation exposed a choice not covered by the issue. Work pauses at that boundary. |
| Validation review | `validation-review` | Human + AI | Build/test/validation results are reviewed against acceptance criteria. |
| Ready to commit | `ready-to-commit` | Human | Human has approved the final diff and commit message. |

State transitions are a human action. An AI assistant proposes a transition in
an issue comment ("scope implemented, requesting `validation-review`") and must
not edit labels, milestones, assignees, or issue state itself.

The labels above do not exist in a fresh clone of this repository; see
[Repository bootstrap](#repository-bootstrap).

## Repository bootstrap

The policies in this document reference branches, labels, and automation that
must exist before the workflow can actually run. Do this once, in order.

1. **Create the CUDA C integration branch.**

   ```bash
   git fetch origin
   git switch -c cuda-c-stable origin/stable
   git push -u origin cuda-c-stable
   ```

2. **Create the labels.** A fresh repository has only the nine GitHub defaults,
   and GitHub silently ignores unknown labels in issue forms, so unbacked labels
   make the state model above a no-op.

   ```bash
   gh label create triage               --color ededed --description "Needs classification"
   gh label create needs-spec           --color fbca04 --description "Scientific/numerical choices still open"
   gh label create ready-for-ai-plan    --color 0e8a16 --description "AI may propose an implementation plan"
   gh label create plan-review          --color 1d76db --description "AI plan awaiting human review"
   gh label create ai-implementation    --color 5319e7 --description "AI is implementing the approved scope"
   gh label create needs-human-decision --color d93f0b --description "Paused on a human decision"
   gh label create validation-review    --color c5def5 --description "Validation evidence under review"
   gh label create ready-to-commit      --color 0e8a16 --description "Diff and commit message approved"
   gh label create ai-assisted          --color bfd4f2 --description "Implemented mostly by an AI assistant"
   gh label create cuda-c               --color 76b900 --description "CUDA C porting milestone"
   gh label create performance          --color d4c5f9 --description "Performance / optimization work"
   gh label create validation           --color c2e0c6 --description "Numerical validation work"
   ```

3. **Confirm the milestone.** The milestone is titled `Porting to CUDA C`.

   ```bash
   gh api repos/htymjun/OUxSBLI/milestones --jq '.[] | "\(.number) \(.title) [\(.state)]"'
   ```

4. **Merge the carry-over automation into `stable`.** GitHub runs a
   `pull_request_target` workflow from the *base* branch of the PR, so
   `cuda-c-followup-from-perf.yml` does nothing while it lives only on a feature
   branch. Until it is on `stable`, carry-over issues must be opened by hand with
   the CUDA C Performance Carry-over template.

5. **Open the foundation issue** and attach it to `Porting to CUDA C`.

## Milestone structure

Use milestones for multi-issue efforts whose value only appears after several
pieces land. For example, `Porting to CUDA C` should not be a single issue:
it is a milestone containing a foundation issue plus small executable child
issues.

Recommended milestone shape:

- **Foundation issue**: defines the minimum runnable architecture and stays open
  until its child issues are complete or the scope is deliberately redefined.
- **Implementation issues**: each changes one coherent code path, kernel family,
  build-system surface, or API boundary.
- **Validation issues**: each owns one visible validation target with commands,
  outputs, tolerances, and plots/tables that a human can inspect.
- **Performance issues**: each owns one benchmark surface, profiling method, and
  before/after comparison.
- **Documentation issues**: record decisions, migration notes, and known
  limitations that should not live only in chat logs.

For a large milestone, the foundation issue should list its child issues and
state which ones block downstream work. Child issues should link back to the
foundation issue and name what they unblock.

## CUDA C porting milestone

For `Porting to CUDA C`, use issue sizes that can be reviewed independently.
Do not make the first issue "rewrite the solver in CUDA C." A good issue should
usually fit one of these scopes:

- expose or preserve one build path without changing numerical behavior
- translate one narrow kernel family or helper layer
- fuse one existing hot path while preserving the current validation result
- add one validation report that produces human-readable evidence
- benchmark one case, grid, GPU, compiler, and rank configuration
- document one migration decision or unresolved risk

### Starting point and prior work

The milestone starts from the current CUDA Fortran `stable` baseline.

The `CUDA_C_JAX_Python` branch contains earlier experiments (`CUDA_solver/` with
`flux.cu`, `visc.cu`, `rk.cu`, `wrapper.cpp`, and a `jax_solver/` tree, last
touched 2025-07-07). That branch is **historical and out of scope for this
milestone**: it is not the baseline, it is not maintained against current
`stable`, and its layout does not define the target architecture. Copying code
or structure from it is an explicit human decision that must be recorded in an
issue before the code appears in a PR — an AI assistant must not pull from that
branch on its own initiative.

### Branch and merge policy

- Keep `stable` as the stable CUDA Fortran baseline.
- Use `cuda-c-stable` as the integration branch for the CUDA C porting
  milestone.
- Name every CUDA C porting branch `cuda-c-*`.
- Open CUDA C porting PRs against `cuda-c-stable`, not `stable`.
- Preserve the existing CUDA Fortran implementation on `cuda-c-stable`; CUDA C
  paths should be introduced alongside it until the human maintainer explicitly
  approves a different policy.
- Duplicated build, pytest, validation, and benchmark surfaces for CUDA Fortran
  and CUDA C are acceptable during the migration.

Nothing enforces this policy automatically. It is checked by the reviewer using
the CUDA C section of the pull-request template.

### Numerical equivalence policy

A language port raises a question that a tolerance number alone does not answer,
so the foundation issue must record the equivalence standard before any kernel is
translated:

- **What is compared**: which case, grid, `dt`, step count, and which quantities
  (conserved variables, kinetic energy, dissipation rate, wall quantities).
- **Against what**: the CUDA Fortran baseline at a named commit, not "the
  previous run".
- **How strict**: bitwise identity or a stated tolerance per quantity.

Bitwise agreement between `nvfortran -fast` and `nvcc` should not be assumed as
the default target. FMA contraction, reduction order, fast-math transformations,
and math-library implementations differ between the two compilers, so an
otherwise correct translation can disagree in the last bits and a
bitwise-or-nothing rule would reject it. If the maintainer does want bitwise
agreement, the issue must say so and accept the resulting constraints on
optimization flags.

Every translated-kernel issue cites the recorded standard. An AI assistant must
not invent, widen, or reinterpret a tolerance; that is a
[human decision gate](#human-decision-gates).

### Performance carry-over policy

- CUDA Fortran remains an active source of optimization ideas because some
  contributors may be more effective there than in CUDA C.
- When a PR into `stable` contains a useful CUDA Fortran performance change,
  the project should evaluate whether the idea belongs on `cuda-c-stable`.
- The evaluation is not an automatic port. It is a human-reviewed CUDA C
  follow-up issue that can decide to port now, defer, or reject.
- `cuda-c-followup-from-perf.yml` watches PRs merged into `stable` and creates
  the follow-up issue when it sees any of: a performance-style label
  (`perf`, `performance`, `optimization`), a conventional `perf:` commit, a
  `perf:` PR title, or a `perf/*` head branch. It attaches the issue to
  `Porting to CUDA C` when that milestone is open, and labels it `cuda-c`,
  `performance`, `ai-assisted`, `triage`.
- The follow-up issue must compare the CUDA Fortran baseline and CUDA C path
  before any performance claim is accepted.

Known limitations of that automation:

- It only sees **pull requests merged into `stable`**. A direct push to `stable`
  produces no follow-up issue.
- It must be present on `stable` to run at all (`pull_request_target` uses the
  base branch's copy of the workflow).
- It attaches a milestone only when an open milestone matches
  `Porting to CUDA C` (compared case-insensitively); otherwise it warns in the
  run log and creates the issue without one.
- Missing labels do not block issue creation: the workflow retries without
  labels and warns, so the labels must then be applied by hand.
- **It detects nothing unless a signal is added on purpose.** None of the last
  twelve pull requests merged into `stable` carried any of the four signals —
  including one titled "Apply SoA and FMA optimization to 2D solver", which is
  exactly the kind of change this automation exists to catch. Branch names in
  this project are topic names (`SoA`, `fused`, `mixed`, `fypp`), titles are not
  conventional commits, and PRs are usually merged without labels. Adding the
  `performance` label before merging is the one habit that makes the carry-over
  path work.

The signals are meant to be steered deliberately. To force an evaluation, add
the `performance` label. To skip one, keep the branch name, PR title, and commit
prefixes free of the `perf` markers, and open a carry-over issue by hand later if
the idea turns out to matter after all.

### Suggested initial issue graph

Work the rows roughly in order; each row's blocking relationship is listed.

| # | Issue | Purpose | Blocks |
|---|---|---|---|
| 0 | Milestone bootstrap | Create `cuda-c-stable`, create the labels, confirm the milestone, and merge the carry-over workflow into `stable`. | Every other row |
| 1 | Foundation: CUDA C porting architecture | Decide the migration boundary, build layout, `cuda-c-stable` policy, interoperability policy, numerical equivalence standard, and first runnable target. | All CUDA C child issues |
| 2 | Fused kernel baseline | Define which current Fortran kernels should be fused before translation, and validate that fusion preserves results. | CUDA C translation of that hot path |
| 3 | Minimal CUDA C build path | Add the smallest CMake/NVHPC/NVCC path that can compile and run an isolated CUDA C kernel or bridge. | Translated kernel work |
| 4 | Validation report format | Add a report section or script output that humans can inspect without reading only pytest logs. | Numerical acceptance of ported kernels |
| 5 | First translated kernel family | Translate one approved kernel family after the fusion boundary is clear. | Later broad translation |
| 6 | Performance benchmark | Compare CUDA Fortran and CUDA C for the same case/configuration with recorded hardware and compiler details. | Performance claims |
| 7 | CUDA Fortran performance carry-over | Evaluate merged CUDA Fortran performance PRs for possible CUDA C adoption. | Optional CUDA C optimization issues |
| 8 | Milestone exit review | Check the exit criteria below and decide merge-back, continuation, or termination. | Closing the milestone |

The exact first fused-kernel target is a human decision. AI may profile and
recommend candidates, but the issue should record the selected target before
implementation.

### Milestone exit, merge-back, and termination

The milestone needs a definition of done, or `cuda-c-stable` becomes a permanent
second baseline by default. The exit review issue (row 8) decides between three
outcomes, and the maintainer owns the decision.

**Merge back into `stable`** requires all of:

- every case that the milestone claimed as in-scope runs on the CUDA C path
- the validation surfaces named in the validation issues pass at the recorded
  equivalence standard, with evidence linked from the exit issue
- a performance comparison against the CUDA Fortran baseline on the same
  hardware, meeting the threshold the maintainer recorded in the foundation
  issue — a slower port is a valid reason not to merge back
- `README.md`, `CONTRIBUTING.md`, `CLAUDE.md`, and the docs describe the build
  and run procedure that survives the merge
- an explicit decision about the CUDA Fortran path: it stays until the
  maintainer says otherwise, and removing it is a separate issue with its own
  review

**Continue** — scope was too large for one milestone; the exit issue records
what landed and opens a successor milestone.

**Defer or terminate** — closing the milestone without merging back is a
legitimate outcome. If the port is abandoned, the exit issue must record why
(performance, maintainability, effort, tooling), and what happens to
`cuda-c-stable`: kept as a reference branch or deleted. That record is the point
of the issue; an abandoned effort with no written reason invites a repeat.

## Validation sections

Do not rely only on pytest output for AI-driven numerical work. Pytest is useful
for automation, but reviewers also need evidence they can read and reason about.

Every solver-behavior issue should include a validation section with:

- commands run
- case path and configuration
- grid size, time step, end time, rank count, and GPU/compiler when relevant
- reference solution, dataset, or baseline commit
- tolerance or acceptance metric
- human-readable result: table, short log excerpt, plot path, VTK snapshot, or
  benchmark summary
- known limitations and unresolved discrepancies

For CUDA C migration, validation should explicitly compare the new path against
the current CUDA Fortran path before performance is interpreted. It is fine for
the migration branch to run both CUDA Fortran and CUDA C pytest/validation
surfaces in parallel.

### Evidence surfaces in this repository

Name the surface you used, so evidence stays comparable across issues instead of
each issue inventing its own.

- **Build**: `cmake -B build && cmake --build build -j` in the case directory,
  for the representative cases the issue names (for example
  `3D_solver/NSTGV`, `3D_solver/ETGV`, `2D_solver/OS`, `3D_solver_curv/NACA`).
  State which cases were built and which were not.
- **Generated Fortran**: when a `.f90.fypp` template changes, inspect the
  generated `.f90` in `<CASE>/build/` for every affected config combination.
  The root `CLAUDE.md` records the specific trap here: no current case builds
  NS/LES with `VISC_ORDER=4`, so that branch is only exercised by overriding
  `VISC_ORDER` in a case's `config.fypp` on purpose.
- **Python tests**: `pytest ouxsbli/tests/`, or the specific files that cover
  the change — `test_etgv.py`, `test_nstgv_ke_eps.py`, `test_evc.py`,
  `test_evc_3d.py`, `test_os.py`, `test_os_3d.py`, `test_bl.py`,
  `test_bl_3d.py`, `test_sbli.py`, `test_corn.py`, `test_st.py`,
  `test_patcher.py`. Quote the summary line, and say which tests were skipped
  and why.
- **Solver runs**: the case path, `config.fypp` values, `mod_globals.f90` grid
  and `dt`, and the rank count. 2D output needs two ranks (the even rank
  computes, the odd rank writes VTK), so `mpirun -n 1` produces no files.
- **Performance**: `profile.sh` (nsys/ncu) with reports kept under `ncu/`, plus
  the kernel times being compared. Never compare numbers from different GPUs or
  SDK versions without saying so.
- **Session records**: longer investigations belong in `docs/plans/`, which is
  excluded from the published docs, so it is the right place for raw logs and
  intermediate reasoning that a reader of the issue does not need.

## Workflow

1. Open or update an issue.

   Use the normal bug, feature, or validation templates for straightforward
   work. Use the AI specification template when the issue is expected to be
   implemented mostly by an AI assistant, and the CUDA C Porting Task template
   for milestone work.

2. Record human-owned decisions before implementation.

   The issue must state the intended solver tree, case(s), scheme/order choices,
   boundary-condition policy, validation target, and acceptance criteria when
   they matter. Unknowns should be marked explicitly as open decisions.

3. Ask the AI for an implementation plan.

   The AI should read the relevant `CLAUDE.md` files, inspect existing patterns,
   and propose a scoped plan. For OUxSBLI, the plan must identify affected
   templates and at least one generated-output check when fypp conditions are
   changed.

4. Human reviews the plan.

   The human can approve the plan, modify it, or narrow it. The AI must not
   silently choose final scientific or numerical behavior when the issue leaves
   a meaningful choice open.

5. AI implements the approved scope.

   The AI may make local engineering decisions that do not change the issue's
   scientific intent, such as matching existing helper patterns or choosing a
   focused test location. If a new physical model, boundary condition, default
   scheme, tolerance, benchmark metric, or public API behavior must be chosen,
   the AI should stop and ask.

6. AI reports validation evidence.

   The report should distinguish compilation, Python tests, short solver runs,
   analytical/reference validation, and performance measurements. A passing
   build is not enough evidence for numerical changes.

7. Human reviews the diff before commit.

   The AI may suggest a commit message, but the final decision to commit belongs
   to the human unless the human has explicitly delegated that action for the
   current issue.

## Human decision gates

The AI must ask for a human decision before committing to any of these:

- a new numerical method, model constant, limiter, sensor, or coefficient
- a default scheme, order, TVD mode, RK mode, precision, or restart behavior
- boundary-condition behavior or extrapolation policy
- validation tolerance, benchmark metric, or reference dataset selection
- the numerical equivalence standard for a ported kernel, or any change to it
- accepting a known numerical discrepancy as tolerable
- reducing a test grid, end time, or tolerance to save runtime
- changing MPI rank count, decomposition direction, or GPU assignment policy
- changing public Python API behavior
- importing code or structure from the historical `CUDA_C_JAX_Python` branch
- deleting tracked data, reference files, cases, or documentation
- changing labels, milestones, or issue state

The AI may decide without asking when the choice is purely local and already
constrained by existing project style, for example:

- using `rg` to inspect code
- editing `.f90.fypp` rather than generated build output
- adding a focused regression test next to similar tests
- reusing an existing helper, case pattern, or documentation structure
- formatting a file consistently with nearby files

## AI implementation checklist

Before editing:

- Link the issue and summarize the approved scope.
- Read the relevant root and subtree `CLAUDE.md` files.
- Identify generated files that must be inspected after fypp changes.
- Identify which tests or solver runs are practical on the current machine.
- For milestone work, confirm the branch is `cuda-c-*` and based on
  `cuda-c-stable`.

Before review:

- Show the changed files.
- Explain numerical and performance impact separately.
- Report exact commands run and whether they passed.
- Report commands that could not be run and why.
- Name the hardware and compiler behind any timing or tolerance claim.
- Call out every open human decision.

Before commit:

- Confirm the issue is in `ready-to-commit`.
- Confirm generated build artifacts are not staged.
- Confirm no unrelated user changes were reverted.
- Use a concise commit message that names the issue's logical change.

## Recommended issue shape

The issue templates in `.github/ISSUE_TEMPLATE/` already encode this; the list
is here for issues opened without a template.

- Goal
- Non-goals
- Affected solver tree and cases
- Human-owned decisions
- Open questions
- Proposed validation
- Acceptance criteria
- Notes for AI

This keeps the issue small enough for an assistant to execute, while keeping
the scientific judgment with the human maintainer.
