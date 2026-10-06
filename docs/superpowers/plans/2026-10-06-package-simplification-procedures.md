# Package simplification — procedures every task follows

Read this whole file before starting any task of
`docs/superpowers/plans/2026-10-06-package-simplification.md`. It holds the commands, the
rules and the per-kind procedures; your task brief names which procedure applies.

Paths:

- Repo: `/home/denis.linux/.julia/dev/EnzymeRates`, branch `package-simplification`.
- Spec: `docs/superpowers/specs/2026-10-06-package-simplification-design.md`.
- Item catalog: `docs/superpowers/specs/2026-10-06-package-simplification-items.md`. Every item
  ID in a task (`types/T5`, `enum/S1`, `time/TT02`, `loc/TY1`, …) has a section there headed
  `### <id> — <title>` (source items) or `#### <id> — <title>` (test items). Find it with
  `grep -n '^#\+ <id> —' <catalog>`. The section is the detailed requirement: what to change,
  what it must not change, the evidence, and the test impact. Line numbers in the catalog refer
  to main at `e59ebd6`; earlier tasks move code, so re-locate by function name.
- Harness: `.superpowers/harness/` (git-ignored).

## Rules

1. **One Julia process at a time.** The machine has 7.7 GB and no swap. Before starting Julia,
   run `pgrep -af julia | grep -v languageserver`; if a test or harness process from an earlier
   step is still running, wait for it or kill it. Never run two Julia processes at once.
2. **Never push, never switch branches, never rewrite history** (no rebase, no amend of earlier
   commits, no reset of commits). Commit on `package-simplification` only.
3. **Do not dispatch subagents.** Review comes from the controller.
4. **Four invariants stay intact:** the `rate_equation` performance contract
   (`test_rate_equation_performance`: 0 allocations, < 120 ns), Canonical Step Form (step
   direction, step order and group order canonicalized in the `Step`, `Mechanism` and
   `AllostericMechanism` constructors), the `name(p, m)` chokepoint (all `Parameter → Symbol`
   rendering goes through it), and Pass-2 Wegscheider absorption
   (`_build_kinetic_rename_map` / `_build_wegscheider_rename_map`).
5. **CLAUDE.md applies in full.** In particular: every file keeps its two `ABOUTME:` lines;
   comments describe the code as it is (no "new", "old", "moved", "refactored", "unified"); never
   remove a comment unless it is false; delete a deleted function's docstring with it; 92-char
   lines, 4-space indents, match the surrounding style; TDD for every bug fix and new behaviour;
   test output stays pristine; never delete a failing test to make the suite pass.
6. **Docs move with the code.** When a change renames or deletes something a page under
   `docs/src/` names, update that page in the same commit (`grep -rn '<name>' docs/src`).
7. **Tests that test only deleted code go with it** (approved items loc/T25 and loc/SRC1). Any
   other test deletion must be an approved catalog item.
8. **Behaviour-preserving means the harness passes.** If `harness.jl check` reports a
   difference in a commit the catalog marks "Behaviour change: none", the change is wrong: fix
   it. Never re-record the snapshot to make a behaviour-preserving commit pass.
9. **Commit messages** use the repo style (imperative subject, a body that says what and why,
   group-prefixed item IDs) and end with these two lines, verbatim:

   ```
   Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
   Claude-Session: https://claude.ai/code/session_015DG9FSn9ZSbZhzeSDkqqnS
   ```

10. **Report** to the report file your dispatch names: what you changed per item ID, net line
    counts (`git diff --stat`), the commands you ran with their pass/fail summary, harness
    result, and anything you could not do. Return only status, commits, a one-line test summary
    and concerns.

## Commands

Focused run of one test file (warm, ~1–5 min; skips Aqua and JET; a helper defined in another
test file may be undefined here, an artifact of the focused run):

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project --heap-size-hint=2500M -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")' 2>&1 | tail -40
```

Full suite (10–20 min, longer than the 10-minute foreground limit). Start it with the Bash
tool's `run_in_background: true`:

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])' > /tmp/enzymerates-suite.log 2>&1; echo "SUITE_EXIT $?" >> /tmp/enzymerates-suite.log
```

then poll in bounded stretches (each call returns within 9 minutes; repeat until it prints the
exit line):

```bash
timeout 540 bash -c 'until grep -q "SUITE_EXIT" /tmp/enzymerates-suite.log; do sleep 20; done'; tail -25 /tmp/enzymerates-suite.log
```

Pass means `SUITE_EXIT 0`, `Testing EnzymeRates tests passed`, and no unexpected warnings or
errors in the log. Before committing, skim the log for warnings your change introduced.

Harness (5–8 min; run the same way, in the background, then poll for `HARNESS_EXIT`):

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project --heap-size-hint=2500M .superpowers/harness/harness.jl check > /tmp/enzymerates-harness.log 2>&1; echo "HARNESS_EXIT $?" >> /tmp/enzymerates-harness.log
```

`check` exits 0 and prints `HARNESS OK` when every recorded file matches; otherwise it prints
the first differing lines per file and exits 1. `record` (same command with `record`) rewrites
the snapshot; only procedure B uses it.

Timed suite (procedure T): `julia --project --heap-size-hint=2500M .superpowers/harness/timed_runtests.jl > /tmp/enzymerates-timed.log 2>&1`
in the background, then poll for `TOTAL`.

## Procedure R — behaviour-preserving refactor or deletion

1. Read every catalog section your task names, then the code it touches.
2. Make the change, one item (or one tightly coupled item group) at a time.
3. After each item: run the focused test files that cover the touched code.
4. When the task's items are done: run the harness `check` → must print `HARNESS OK`. Run the
   full suite → must pass.
5. Commit. One commit per item group the task names (the task says where to split).

## Procedure B — approved behaviour change or bug fix

1. Read the catalog section and the code.
2. TDD: write the test that pins the new behaviour; run it and watch it fail for the expected
   reason. (For a pure removal of an internal method, the failing check is a test asserting the
   method is gone or the call site that changes.)
3. Implement; run the focused tests until green.
4. Run the harness `check`. Confirm every reported difference is exactly the approved change and
   nothing else (write the diff summary into your report). Then run
   `harness.jl record` to accept it.
5. Run the full suite; re-pin any test that pinned the old behaviour, and only those.
6. Commit alone, one commit per behaviour change. The body states what changed for users and how
   many recorded mechanisms or strings moved.

## Procedure S — test-speed or test-size change (test files only)

1. Read the catalog sections and the testsets they name.
2. Edit the test. Deletions must be the approved catalog item and keep every assertion the
   catalog says is unique to the deleted test (move it to the subsuming test where the catalog
   says so).
3. Run the focused test file before and after and record both times from the `Test Summary`
   line in your report.
4. Run the full suite; commit. No harness run is needed when `src/` is untouched.

## Procedure T — timing checkpoint

Run the timed suite, then paste into your report the `FILE` lines, the `TOTAL` line, and the 15
slowest testsets from the verbose summary (`grep -E '^ +[^|]+\|' log | sort -t'|' -k3`).
