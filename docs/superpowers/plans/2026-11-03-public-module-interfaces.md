# Public Module Interfaces Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Declare each module's interface with Julia's `public`, raise the minimum Julia to 1.11, and revise `.claude/CLAUDE.md` so that simplifying the package is the default.

**Architecture:** Each `const API` tuple becomes a `public` statement listing the same names, and the existing `Meta.parseall` scanner reads that statement, so all four boundary rules, the docs tables and the diagram keep their meaning. The rule text comes verbatim from the spec.

**Tech Stack:** Julia 1.11+, `Test`, Aqua, JET, Codecov.

**Spec:** `docs/superpowers/specs/2026-11-03-public-module-interfaces-design.md`

## Global Constraints

- `Project.toml`: `julia = "1.11"`. No new dependency anywhere, ExplicitImports and coverage tools included.
- Each module's public names equal its old `API` list. No interface name is added or removed.
- No compile budget or timing bound is loosened. A failure that is not a real bug goes to Denis.
- `CLAUDE.md` and `codecov.yml` text is copied verbatim from the spec. The rule lists are closed: report any other conflict to Denis.
- 92-character lines, `_testhelper_` prefix on test helpers, one Julia process at a time.
- Full suite: `timeout 1200 julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`, in the background.
- Focused run: `julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_boundaries.jl")'`
- Never push, open a PR or merge. `git add` named files only; `.claude/settings.local.json` stays out of every commit.

## Review Focus

- A typo in a `public` statement: Julia accepts a name that does not exist, so a test must fail on it (Task 1).
- A `public` statement split in two within a module: the scanner must read both, not the last (Task 1).
- A module with no `public` statement: the scanner must fail loudly and not report an empty interface (Task 1).
- `public` listing a macro: `@enzyme_mechanism` must match the `Symbol("@enzyme_mechanism")` the scanner finds at call sites (Task 1).
- The new minimum itself: nothing else runs Julia 1.11, so one suite run must (Task 2).

---

### Task 1: `public` replaces `const API`

**Files:**
- Modify: `src/{Mechanisms,Derivation,Fitting,Enumeration,DSL,Data}/<Module>.jl`, `Project.toml:42`
- Modify: `test/module_map.jl:1,82-98,123`, `test/test_boundaries.jl`, `docs/src/developer.md:15-18,116,314-319`

**Interfaces:**
- Produces: `_testhelper_api(root, mod::Symbol) -> Set{Symbol}`, the same signature and result as today, read from `public` statements.

- [ ] **Step 1: Move the fixtures to `public` and add the new cases**

In `test/test_boundaries.jl`, rewrite every fixture `const API = (:a, :b)` as `public a, b` (lines 150-155, 173, 207-210, 230-231, 276, 280). `Symbol("@enzyme_mechanism")` becomes `@enzyme_mechanism`. In the testset "module_sources and api", add:

```julia
_testhelper_with_tree(
    "Fitting/Fitting.jl" => "public fit\npublic refit") do root
    @test _testhelper_api(root, :Fitting) == Set([:fit, :refit])
end
_testhelper_with_tree("Fitting/Fitting.jl" => "fit() = 1") do root
    @test_throws ErrorException _testhelper_api(root, :Fitting)
end
```

Add one testset after "module boundaries":

```julia
@testset "every public name is defined" begin
    for m in _testhelper_MODULES, n in _testhelper_api(_testhelper_SRC, m)
        @test isdefined(getfield(EnzymeRates, m), n)
    end
end
```

- [ ] **Step 2: Run to verify it fails**

Run the focused command. Expected: FAIL, with `declares no` errors from `_testhelper_api` on the rewritten fixtures.

- [ ] **Step 3: Teach the scanner to read `public`**

In `test/module_map.jl`, delete `_testhelper_api_name`. `_testhelper_api` walks the file and unions the arguments of every `ex.head === :public` node into a `Set{Symbol}`; their arguments are plain symbols, macros as `Symbol("@name")`. With no such node it errors with `"$path has no public statement"`.

- [ ] **Step 4: Run to verify the fixtures pass and `src/` fails**

Run the focused command. Expected: the fixture testsets pass; those reading `src/` error with `has no public statement`.

- [ ] **Step 5: Replace each `const API` and raise the minimum**

In each `src/<Module>/<Module>.jl`, replace the tuple with one `public` statement of the same names, wrapped at 92 characters (`public parse_mechanism, @enzyme_reaction, @enzyme_mechanism, @allosteric_mechanism`). In `Project.toml`, set `julia = "1.11"`.

- [ ] **Step 6: Update the wording**

- `test/test_boundaries.jl`: the `ABOUTME` and testset names say "public names". "every cross-module name is in its API list" becomes "every cross-module name is public in its module"; "rule 1 finds a name outside the API list…" and "rule 3 finds underscore names in an API list…" follow. The docstring at line 56 follows.
- `test/module_map.jl`: the `ABOUTME` and the docstrings at lines 123 and 135 follow.
- `docs/src/developer.md`: line 15-16 becomes "lists its interface in a `public` statement"; line 18 and line 116 say "public names"; rule 1 becomes "that is not public in that module"; rule 3 becomes "an underscore name in a `public` statement or used across modules". Lines 319-320 become: "To use another module's function, it must be public there, under a name without an underscore."

- [ ] **Step 7: Run to verify everything passes**

Run the focused command. Expected: PASS, including "the developer docs' module map matches src/" and "Mechanisms API is the data model", neither of which changes. A table change here would mean a name changed.

- [ ] **Step 8: Prove the names did not change**

Run: `git diff -U0 -- src | grep -E '^[-+](const API|public|    )'` and compare the removed and added names by eye. Expected: the same names per module.

- [ ] **Step 9: Run the full suite and commit**

Run the full-suite command. Expected: all pass, with no warnings in the log. Commit `Project.toml`, the six module files, both test files and `docs/src/developer.md` as "Modules declare their interfaces with public; the minimum Julia is 1.11".

### Task 2: Run the suite on the new minimum

**Files:**
- Modify: `test/test_compile_budget.jl:12-14`, `ROADMAP.md:72-84,234`

- [ ] **Step 1: Run on Julia 1.11 in a clean worktree**

A separate worktree keeps the 1.13 `Manifest.toml` intact.

```bash
juliaup add 1.11
git worktree add "$SCRATCHPAD/wt11" HEAD && cd "$SCRATCHPAD/wt11"
timeout 1800 /usr/bin/time -v julia +1.11 --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])' 2>&1 | tee "$SCRATCHPAD/full_1_11.log"
```

Expected: all pass. Record the Julia patch version and "Maximum resident set size". For a failure, find the root cause: fix a real bug in its own commit with a test, and report anything else to Denis without loosening a bound. If the OOM killer ends the run, report that with the peak.

- [ ] **Step 2: Measure the init trace count on 1.11**

In the worktree, print the count that `test/test_compile_budget.jl` compares with `INIT_TRACE_BUDGET`, by a temporary print that is never committed.

- [ ] **Step 3: Update the budget comment**

Replace "and 87 on Julia 1.10" with the 1.11 figure, and adjust the `baseline 2026-10-09: 57-87` range to the figures that remain. `INIT_TRACE_BUDGET` stays 200.

- [ ] **Step 4: Update `ROADMAP.md`**

If 1.11 peaked under 4 GB, delete "On Julia 1.10 the full suite can peak above 4 GB". Otherwise rewrite its measurements for 1.11, keeping the ter-ter garbage explanation. In "Test the compat bounds", write `julia = "1.11"`.

- [ ] **Step 5: Clean up and commit**

Remove the worktree with `git worktree remove`. Run `test/test_compile_budget.jl` as a focused run. Expected: PASS. Commit both files as "Measure the suite on Julia 1.11".

### Task 3: Codecov target and housekeeping

**Files:**
- Create: `codecov.yml`
- Modify: `.gitignore`, `ROADMAP.md` (Ideas), `src/Mechanisms/names.jl:1`

- [ ] **Step 1: Add `codecov.yml`**

Use the YAML block from the spec, under "Elsewhere", verbatim.

- [ ] **Step 2: Validate it**

Run: `curl -s --data-binary @codecov.yml https://codecov.io/validate`. Expected: `Valid!`.

- [ ] **Step 3: Add the rest**

- `.gitignore`: add `.claude/settings.local.json`.
- `ROADMAP.md`, under Ideas: "Review `identify_rate_equation`'s keyword arguments" (at least 18, such as `eq_complexity_filter = 337`), which become internal constants unless users have a real reason to change them; and "Make `DSL.parse_mechanism` private if nothing needs it", a public name no other module in `src/` uses.
- `src/Mechanisms/names.jl`: add a two-line `ABOUTME` header saying what the file holds, read from its contents.

- [ ] **Step 4: Commit**

Run: `git status --short`. Expected: no `.claude/settings.local.json`. Commit the four files as "A fixed Codecov patch target; ignore local settings; roadmap ideas".

### Task 4: Rewrite `CLAUDE.md`

**Files:**
- Modify: `.claude/CLAUDE.md`

Denis approved this wording by approving the spec. Run this task in the main session, not in a subagent.

- [ ] **Step 1: Apply "Added"**

Copy each block quote from the spec's "Added" section to the place its heading names, verbatim. In "Module boundaries (guard)", the three-paragraph quote and its enforcement paragraph replace the whole current paragraph.

- [ ] **Step 2: Apply "Changed" and "Deleted"**

Apply each table row and each bullet verbatim. The sections "Issue tracking" and "Learning and Memory Management" go entirely, along with both journal bullets in "Our relationship" and the sentence "Commit your journal entries."

- [ ] **Step 3: Check what remains**

Run: `grep -n -i 'journal\|TodoWrite\|const API\|SMALLEST\|PROVE they\|park on\|about 26 ns\|Article' .claude/CLAUDE.md`. Expected: no output.

Check that every backticked name and path in "Instructions for this repository" still exists in the tree. Report any that does not, without editing it.

- [ ] **Step 4: Commit**

Commit `.claude/CLAUDE.md` as "CLAUDE.md: simplification by default, deep modules, and rules on approval, dispatch, comments, docs and names".

### Task 5: Final verification

**Files:**
- Delete: the spec and this plan, after the final review.

- [ ] **Step 1: Run the full suite on the final tree**

Run the full-suite command. Expected: all pass, with no warnings in the log.

- [ ] **Step 2: Build the docs**

Run: `julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl")'`, then `git checkout docs/Project.toml`. Expected: exit 0, with only the deployment-environment warning.

- [ ] **Step 3: Whole-branch review, then delete the working documents**

After the review and any fixes, `git rm` the spec and this plan, commit, and stop. The PR opens after #99 merges and the branch rebases onto `main`; both are Denis's calls. The Codecov check on that PR confirms the 95% target.
