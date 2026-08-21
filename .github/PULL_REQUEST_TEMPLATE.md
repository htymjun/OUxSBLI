## Summary

Describe what this PR changes.

## Related Issue

Closes #

For AI-assisted work, link the AI specification issue and summarize any human
decisions that constrained the implementation.

For `Porting to CUDA C` milestone work, the PR branch should be `cuda-c-*` and
the base branch should be `cuda-c-stable`.

## Type of change

- [ ] Bug fix
- [ ] New feature
- [ ] Numerical method
- [ ] Performance optimization
- [ ] Refactoring
- [ ] Documentation
- [ ] Test / validation
- [ ] Other

## Scientific / technical motivation

Why is this change needed?

## Implementation

Briefly describe the implementation.

## Validation

Describe how the change was verified.

- [ ] Existing tests pass
- [ ] New tests added
- [ ] Numerical validation performed
- [ ] Reference solution compared
- [ ] Benchmark performed
- [ ] Not applicable

### Results

Provide human-readable results, plots, benchmark numbers, error values, VTK
snapshot paths, or report paths. Do not rely only on pytest pass/fail output
for solver-behavior changes. There is no build or test CI in this repository, so
state the machine behind every result: GPU model, HPC SDK version, MPI rank
count, grid size, `dt`, and end time.

## Performance impact

For performance-related changes:

| Metric | Before | After |
|---|---:|---:|
| Runtime | | |
| GPU utilization | | |
| Memory usage | | |
| MPI communication | | |

If this is a CUDA Fortran performance PR targeting `stable` and the optimization
should be evaluated for carry-over to `cuda-c-stable`, give it at least one
carry-over signal: a performance-style label, a `perf:` commit, a `perf:` PR
title, or a `perf/*` branch name.

## Numerical impact

- [ ] No known numerical impact
- [ ] Accuracy changed
- [ ] Stability changed
- [ ] Reproducibility may change
- [ ] Not applicable

## AI-assisted development

- [ ] Not AI-assisted
- [ ] AI plan was reviewed before implementation
- [ ] Human-owned scientific / numerical decisions came from the issue
- [ ] Validation section includes human-readable evidence, not only pytest output
- [ ] AI-raised open decisions were resolved by a human
- [ ] Final diff was reviewed before commit

## Compatibility

- [ ] Existing cases remain compatible
- [ ] Configuration changes required
- [ ] Breaking change

## CUDA C Porting

- [ ] Not CUDA C porting work
- [ ] Branch name is `cuda-c-*`
- [ ] Base branch is `cuda-c-stable`
- [ ] Existing CUDA Fortran path remains available
- [ ] CUDA Fortran and CUDA C validation surfaces are both documented where applicable

## Checklist

- [ ] Code follows the existing project structure
- [ ] Documentation updated if necessary
- [ ] Tests/validation performed
- [ ] No generated build files committed
- [ ] PR is focused on one logical change
