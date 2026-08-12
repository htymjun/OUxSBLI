# Split CLAUDE.md into a hierarchical set of memory files

## Context

The root `CLAUDE.md` (434 lines / ~26KB) has grown to cover four largely
independent subsystems — the 3D Cartesian solver, the curvilinear solver, the
2D solver, and the Python test suite — in one flat file. Every session loads
all of it regardless of which subsystem is actually being touched, and the 2D
solver section in particular has accumulated a large amount of hard-won,
case-specific debugging narrative (the BL/SBLI boundary-condition post-mortem)
that has little relevance when someone is working in `3D_solver/`.

Claude Code supports **nested memory files**: a `CLAUDE.md` placed inside a
subdirectory is loaded automatically once Claude reads/edits a file in that
subtree, in addition to the root file. Splitting along the existing directory
boundaries — `3D_solver/`, `3D_solver_curv/`, `2D_solver/`, `ouxsbli/` — lets
each subsystem's detail load only when relevant, while the root file shrinks
to the cross-cutting material (project identity, shared config reference,
tooling rules) plus a short map pointing to the nested files.

This is a **reorganization, not a rewrite**: all prose must move verbatim
(including the specific numbers — CFL values, timings, tolerances — in the
2D solver lessons-learned section). Nothing gets summarized or dropped, only
relocated, and a small number of cross-references need updating because the
text they point to moves to a different file.

Confirmed via an Explore pass: no nested `CLAUDE.md` files exist yet, the
target directories (`3D_solver/`, `3D_solver_curv/`, `2D_solver/`, `ouxsbli/`)
all exist with the case lists the current CLAUDE.md already describes, and
`.claude/settings.json` has no memory-loading overrides that would affect this.
The repo-root `docs/` (mkdocs site) and `ouxsbli/tests/utils/` etc. are
unrelated to this change.

## Files to create/modify

- **`CLAUDE.md`** (root) — trimmed to cross-cutting content + a documentation map
- **`3D_solver/CLAUDE.md`** (new) — 3D Cartesian architecture, fypp template table, cases, schemes, cell-center gradients, kernel-performance Notice items
- **`3D_solver_curv/CLAUDE.md`** (new) — curvilinear architecture and data flow
- **`2D_solver/CLAUDE.md`** (new) — 2D cases table + the full BL/SBLI boundary-condition/tuning narrative
- **`ouxsbli/CLAUDE.md`** (new) — Python test suite, analytical helpers, post-processing utilities

`tutorials/` stays undocumented-by-nested-file — its section (L389-396) is a
few lines and doesn't justify a fifth file; it stays at root as-is.

## Section-by-section mapping (current line numbers)

| Current section (lines) | Destination | Notes |
|---|---|---|
| L1-7 Title, intro, "What This Project Is" | root | unchanged |
| L9-46 Build & Run | root, trimmed | keep compiler requirement + output format; keep one command block per solver but shortened, since this is the first thing a new session needs — do not reduce to just prose |
| L48-73 Fypp "How it works" (config.fypp concept + example) | root | generic across all three solvers, stays as-is |
| L75-96 fypp-preprocessed source file table (3D-specific) | `3D_solver/CLAUDE.md` | verbatim table |
| L97 "2D solver uses CMake... `2D_solver/src/calc_flux_base.f90.fypp`..." | `2D_solver/CLAUDE.md` | one line, pairs naturally with 2D's own build section |
| L99-109 "Changing the scheme or method" | root, generalized | replace the 3D-only example with a solver-agnostic one-liner ("edit config.fypp in the case dir, rebuild — same pattern for 2D_solver/<CASE> and 3D_solver_curv/<CASE>") |
| L111-117 Architecture intro + 3D directory-layout bullets | `3D_solver/CLAUDE.md` | bullets name only 3D paths |
| L119-132 Curvilinear Solver subsection | `3D_solver_curv/CLAUDE.md` | verbatim |
| L134-145 Data Flow (Cartesian) | `3D_solver/CLAUDE.md` | verbatim |
| L147-159 Data Flow (Curvilinear) | `3D_solver_curv/CLAUDE.md` | verbatim, grouped with Curvilinear Solver |
| L161-233 2D Solver subsection (cases table + full BL/SBLI narrative) | `2D_solver/CLAUDE.md` | verbatim, this is the bulk of the file size reduction |
| L235-246 3D Cartesian Cases table | `3D_solver/CLAUDE.md` | verbatim |
| L248-257 Convective Schemes table | `3D_solver/CLAUDE.md` | verbatim |
| L259-286 Cell-Center Velocity Gradients | `3D_solver/CLAUDE.md` | verbatim; internal "see calc_Ducros note in Notice below" reference stays valid since both live in the same new file |
| L288-331 Configuration (config.fypp / mod_globals.f90 / mod_constant.f90) | root | cross-cutting reference table, stays |
| L333-335 Conservative Variable Layout | root | stays |
| L337-341 MPI Decomposition | root | stays |
| L343-387 Python Test Suite | `ouxsbli/CLAUDE.md` | verbatim; root keeps a 2-line pointer + the `pytest ouxsbli/tests/` command |
| L389-396 Tutorials | root | stays (too small to split out) |
| L398-406 "Adding a New Test Case: For 3D Cartesian" | `3D_solver/CLAUDE.md` | verbatim |
| L408-414 "Adding a New Test Case: For 2D cases" | `2D_solver/CLAUDE.md` | verbatim |
| L416 Notice bullet 1 (occupancy/streams, calc_flux_base) | `3D_solver/CLAUDE.md` | |
| L417-419 Notice bullets 2-3 (calc_flux_base/calc_steps bottleneck, calc_div chain) | `3D_solver/CLAUDE.md` | |
| L420 Notice bullet 4 (Roe unused, optimize KEEP/SLAU/Hybrid) | `3D_solver/CLAUDE.md` | |
| L421 Notice bullet 5 (id_accuracy) | root | generic, pairs with Configuration section |
| L422 Notice bullet 6 (no contiguous/shared on shared-mem args) | root | generic coding rule, applies to any kernel work |
| L423 Notice bullet 7 (cpu_gpu_mpi.f90 future change) | root | generic heads-up on a shared file |
| L424 Notice bullet 8 (STZ validation status) | `3D_solver/CLAUDE.md` | STZ is a 3D case |
| L425 Notice bullet 9 (fypp template verification caution, VISC_ORDER=4) | root | generic caution about `.f90.fypp` templates broadly, pairs with Fypp Preprocessing section |
| L426 Notice bullet 10 (calc_Ducros / calc_hybrid, SBLI) | `3D_solver/CLAUDE.md` | references calc_div, which is 3D-specific |
| L427 Notice bullet 11 (ncu profiling captures path) | `3D_solver/CLAUDE.md` | path is `3D_solver/nsys_ncu/` |
| L429-433 Strict Tooling Rules | root | unchanged, global |

## Cross-references that need fixing

- `ouxsbli/CLAUDE.md` (moved Python Test Suite section) currently says "see
  the cost note under '2D Solver' above" (orig L365) and "described under '2D
  Solver' above" (orig L376) — both need to become explicit pointers to
  `2D_solver/CLAUDE.md` since that section is no longer "above" in the same
  file.
- Root's new documentation map (below) is itself the fix for the general
  "where did X move to" problem for a human skimming the file.

## New root CLAUDE.md structure

1. Title + intro (unchanged)
2. **New: "Documentation Map"** section right after the intro, one line per
   nested file naming what it covers, e.g.:
   - `3D_solver/CLAUDE.md` — 3D Cartesian build/fypp templates, cases, convective schemes, cell-center gradient optimization, kernel-performance notes
   - `3D_solver_curv/CLAUDE.md` — curvilinear O-grid architecture, data flow, NACA/CORN cases
   - `2D_solver/CLAUDE.md` — 2D cases, and the BL/SBLI shared flat-plate setup with the boundary-condition/tuning lessons learned
   - `ouxsbli/CLAUDE.md` — Python test suite, analytical helpers, post-processing utilities
3. Build & Run (trimmed)
4. Fypp Preprocessing System (generic parts only)
5. Architecture (short paragraph + pointers, no per-solver bullets)
6. Configuration (unchanged)
7. Conservative Variable Layout (unchanged)
8. MPI Decomposition (unchanged)
9. Python Test Suite → 2-line pointer + pytest command only
10. Tutorials (unchanged)
11. Adding a New Test Case → pointer line only ("see 3D_solver/CLAUDE.md or 2D_solver/CLAUDE.md")
12. Notice (only the 4 generic bullets)
13. Strict Tooling Rules (unchanged)

Each new nested file opens with one line orienting a reader who opens it
directly, e.g. "This documents the 3D Cartesian solver; see the repository
root CLAUDE.md for project overview, shared configuration reference, and
tooling rules."

## Verification

- After the split, `wc -w` the root file and confirm it's substantially
  shorter (rough target: under 40% of the original 26KB).
- `wc -w` across root + all 4 nested files and confirm the total word count
  is close to the original (allowing for the small amount of new pointer/map
  text added) — this catches accidental content loss.
- `grep -n "above\|below"` across the new files to catch any remaining
  stale same-file cross-references that now point across files.
- Manually diff the moved 2D-solver narrative (BL/SBLI section) against the
  original line-for-line — this is the section where losing a specific number
  (a CFL value, a timing, a tolerance) would silently erase debugging history.
- Confirm no other file in the repo references `CLAUDE.md` by path in a way
  that assumes single-file content (quick `grep -rn "CLAUDE.md" --include=*.md --include=*.py --include=*.sh .` excluding `.git`).
