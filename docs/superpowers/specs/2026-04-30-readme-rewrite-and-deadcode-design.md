# README rewrite and cheap dead-doc cleanup

## Goal

Bring the user-facing documentation in line with the package's current state.
The existing `README.md` and `SPEC.md` predate several major refactors
(mechanism-types refactor, allosteric-state refactor, mechanism-enumeration
refactor, `identify_rate_equation` implementation) and now describe APIs that
no longer exist alongside features that already work. Replace the README
with a focused 5-section introduction targeting the package's core
functionality, delete documentation that has become obsolete, and add a
lightweight test that keeps the README's code blocks runnable in CI.

Source-level dead code in `src/` and `test/` is **out of scope for this
phase** — it requires its own analysis pass once the README has clarified
what the public surface actually is. That cleanup gets a separate spec.

## What's wrong with the current README

A non-exhaustive list:

- Calls `IdentifyRateEquationProblem`, `IdentifyRateEquationResults`, and
  `identify_rate_equation` "not yet implemented" — they're implemented in
  `src/identify_rate_equation.jl` and tested in
  `test/test_identify_rate_equation.jl`.
- DSL examples for `@enzyme_mechanism` use a nested `species: begin … end`
  block; the current DSL is flat (`substrates:`, `products:`, `enzymes:` at
  top level).
- DSL examples for allosteric mechanisms reuse `@enzyme_mechanism` with a
  `metabolites:` block; the allosteric macro is now a separate
  `@allosteric_mechanism`.
- Type signatures shown as
  `EnzymeMechanism{Species,Reactions,EquilibriumSteps,ParamConstraints}`;
  current type is `EnzymeMechanism{Metabolites, Reactions}` — kinetic groups
  encode what used to be separate type parameters.
- Documents an `enumerate_mechanisms` function and an "8-stage pipeline"
  that no longer exist; the current API is composable
  `init_mechanisms` / `expand_mechanisms` / `dedup!` building blocks.
- References a `graph(m)` accessor that no longer exists.
- "Known Limitations" section warns about slow `rate_equation` compilation
  on large mechanisms; the SS-step cap and factored allosteric expansion
  have made this a non-issue. If a user hits it now, that's a bug, not a
  documented limitation.
- `SPEC.md` overlaps heavily and adds its own staleness (a "Migration from
  Current Format (NOT YET DONE)" section describing a completed migration).

## Out of scope

- **Source-level dead code in `src/`.** Phase 2 — separate spec, separate
  effort. Includes possibly-unused helpers, accessors, exported symbols
  that are no longer covered by tests, and vestigial branches from before
  the type-parameter refactor.
- **Test scaffolding cleanup.** `test/outcmaes/`, leftover diagnostic test
  specs, etc. Phase 2.
- **`CLAUDE.md` updates.** Internal-facing, not user-facing; can drift one
  release cycle longer without harm. Will get caught when its assertions
  start blocking work.
- **Building a documentation site (Documenter.jl).** Decision recorded:
  defer until the package needs more than a single README. When that
  happens, switch.

## New README structure

Five sections plus a one-line install snippet, in this order.

### Section 1 — Brief description

Two short paragraphs.

- Paragraph 1: what the package does. Identify the best enzyme rate
  equation from kinetic data by enumerating biochemically valid
  mechanisms, fitting each to data, and selecting the simplest mechanism
  that adequately describes the data via cross-validation.
- Paragraph 2: differentiator. First-class support for MWC allostery and
  steady-state / rapid-equilibrium hybrid mechanisms, with automatic
  Haldane / Wegscheider thermodynamic constraints derived from the
  mechanism's cycle structure.

Followed by:

```julia
# README-SKIP-IN-TEST
using Pkg
Pkg.add(url="https://github.com/DenisTitovLab/EnzymeRates.jl")
```

### Section 2 — Define a mechanism, derive its rate equation, simulate, fit

Running example: a uni-uni reaction `S ⇌ P` catalyzed by an MWC homodimer
in which substrate, product, and an allosteric activator all bind only in
the R conformation (`:OnlyR` everywhere). This example shows MWC machinery
end-to-end while keeping the rate equation small enough to print on a
single page.

The section walks through:

1. Define the mechanism with `@allosteric_mechanism`. Catalytic site has
   `[E,S] ⇌ [ES]` (`::OnlyR`), `[ES] <--> [EP]` (SS catalytic interconversion),
   `[EP] ⇌ [E,P]` (`::OnlyR`). One regulatory site at oligomeric multiplicity 2
   with ligand `A::OnlyR`.
2. Inspect: `parameters(m)`, `rate_equation_string(m)`, single
   `rate_equation(m, concs, params)` evaluation.
3. Simulate data from `rate_equation` on a concentration grid (S × A
   varying, P at two levels including 0; ~50 points across multiple
   `group` values), with multiplicative log-normal noise. `Random.seed!`
   for reproducibility.
4. Fit recovered parameters with `FittingProblem` + `fit_rate_equation`,
   show recovered values next to the true values.

Concrete deliverable: every code block in this section runs in a single
REPL session in under a few seconds (excluding any compile-time penalty
on first call).

### Section 3 — Define an `EnzymeReaction`, recover the mechanism via `identify_rate_equation`

Same chemistry, now defined as a *reaction* (no mechanism specified):

```julia
rxn = @enzyme_reaction begin
    substrates: S
    products:   P
    allosteric_regulators: A
    oligomeric_state: 2
end
```

Reuse `data` from section 2. Construct
`IdentifyRateEquationProblem(rxn, data; Keq=…)`. Call
`identify_rate_equation(prob; optimizer=…, max_param_count=…, …)`. The
call block is marked `# README-SKIP-IN-TEST` so CI does not run it
(`test_identify_rate_equation.jl` already exercises this path with
reduced settings).

Show:

- `results.best` — the recovered mechanism
- `rate_equation_string(results.best)` — confirms it matches the section-2
  mechanism (or describes the equivalent recovered form)
- A `head(results.cv_results)` snippet — first few rows so the reader sees
  CV scores and fitted parameter columns

### Section 4 — How rate-equation derivation works (biochemist intuition)

Three subsections, ~1 paragraph each. No math derivations beyond pointing
at the printed rate equation; biochemist-friendly prose.

- **Steady-state vs rapid equilibrium.** `<-->` denotes a QSSA step
  (King-Altman applies; kf and kr are independent parameters). `⇌` denotes
  a rapid-equilibrium step (only the binding constant K matters; kf and kr
  collapse to one parameter). A typical mechanism mixes both, and the
  package handles the mixed Cha-style derivation automatically.
  `parameters(m)` reflects this — each RE step contributes one K, each SS
  step contributes a kf and a kr.
- **Haldane and Wegscheider relationships.** When the mechanism contains
  thermodynamic cycles, the rate constants around each cycle must satisfy
  the equilibrium constant. The package detects these cycles
  automatically, declares some k's as dependent, and computes them from
  the rest plus a user-supplied `Keq`. The user fits *independent* k's;
  dependent k's are derived. Identifiability is structural, and
  `structural_identifiability_deficit(m)` reports it.
- **MWC allostery (R/T conformations).** Two-state model: enzyme exists
  in R (active) and T (inactive) conformations with conformational
  equilibrium `L = [T]/[R]`. Each kinetic group and each regulator can
  independently be `:OnlyR`, `:OnlyT`, `:EqualRT`, or `:NonequalRT`. The
  full rate equation is the sum of R-state and T-state contributions
  weighted by the partition function, raised to oligomeric power `n`.
  Refer back to the printed rate equation in section 2 as a worked
  example of these factors.

### Section 5 — How `identify_rate_equation` and mechanism enumeration work

Three subsections, biochemist intuition only.

- **Enumeration as composable building blocks.** `init_mechanisms(reaction)`
  produces all biochemically minimal mechanisms by combining catalytic
  topologies with dead-end inhibition subsets, merging steps that bind the
  same metabolite into shared kinetic groups. `expand_mechanisms(specs,
  reaction)` applies single-move expansions (RE→SS conversion, splitting a
  kinetic group, adding a dead-end regulator, becoming allosteric, changing
  an allosteric state), keyed by estimated parameter count. `dedup!(cache)`
  canonicalizes specs and removes duplicates. The enumeration is grounded
  in chemical reasoning, not blind combinatorics: a step is "elementary"
  only if it changes one site by one event with atom balance preserved.
- **Beam search across parameter counts.** Fit all init mechanisms → keep
  the top fraction by training loss → expand survivors to the next
  parameter-count level → fit, dedup, keep top → repeat until no
  improvement or `max_param_count` reached. Beam width balances coverage
  against runtime.
- **Model selection by leave-one-group-out CV.** Top-N fitted candidates
  per parameter count are LOOCV'd. The "best" mechanism is the one with
  minimum training loss at the parameter count whose CV score is lowest —
  the simplest mechanism that generalizes. The `group` column in the data
  defines folds (one experiment per group, sharing a single `E_total`).

## Runnability test

A new test file `test/test_readme_runs.jl`:

```julia
# ABOUTME: Extracts ```julia blocks from README.md and runs them in one
# ABOUTME: REPL session, skipping blocks tagged with # README-SKIP-IN-TEST.

using Test
using EnzymeRates
using Random

@testset "README runs" begin
    md = read(joinpath(@__DIR__, "..", "README.md"), String)
    blocks = String[]
    for m in eachmatch(r"```julia\n(.*?)\n```"s, md)
        block = m.captures[1]
        startswith(strip(block), "# README-SKIP-IN-TEST") && continue
        push!(blocks, block)
    end
    @test !isempty(blocks)

    script = join(blocks, "\n\n")
    sandbox = Module()
    Core.eval(sandbox, :(using EnzymeRates, Random))
    Core.eval(sandbox, Meta.parse("begin\n$script\nend"))
end
```

Wired into `test/runtests.jl` alongside other test files.

### Convention

- All code blocks use the ` ```julia ` fence (preserves syntax highlighting
  on GitHub).
- A block is skipped from extraction if its **first non-blank line** is
  `# README-SKIP-IN-TEST`.
- All extracted blocks share one anonymous `Module()` — module-level state
  (e.g., `m`, `data`, `params`) defined in earlier blocks is visible to
  later blocks. Snippets do **not** need to be individually runnable.

### Blocks expected to carry the skip marker

- The `Pkg.add(...)` installation snippet.
- The `identify_rate_equation(prob; …)` call in section 3 (multi-minute
  search; `test_identify_rate_equation.jl` already exercises the same
  path with reduced settings).

### Known regex limitations (acceptable as written)

- The regex `r"```julia\n(.*?)\n```"s` does not handle nested triple
  backticks. The README plan does not include any. If we ever introduce
  some, the regex needs upgrading.
- Code blocks indented inside list items (Markdown allows fenced blocks
  indented by 4 spaces inside lists) won't match. The README plan keeps
  all code blocks at the document's top level.

If a block fails, the test reports a stack trace pointing into the
concatenated script. Acceptable for a single-author single-file README.
If locating blame becomes painful, add per-block boundary markers to the
concatenated script.

## Cheap doc cleanup

Files to delete after a grep pass confirms nothing in `src/`, `test/`,
`Project.toml`, or `.github/` references them.

| File | Reason |
|------|--------|
| `SPEC.md` | Superseded by the new README; contains stale "(NOT YET DONE)" migration note. |
| `CODE_SIMPLIFICATION_PROMPT.md` | Repo-root one-shot prompt artifact for an earlier refactor. |
| `PLAN_IMPLEMENTATION_PROMPT.md` | Repo-root planning artifact. |
| `PLAN_RESS_DEDUP.md` | Planning artifact for a completed refactor. |
| `ralph.sh` | Developer scratch script (agent runner); verify before deleting. |
| `.ralph-logs/` | Empty log dir for `ralph.sh`; orphaned once the script is gone. Untracked but exists on disk. |
| `scripts/verify_counts.py` | Topology-count verification used during the enumeration refactor; no longer referenced. Delete the empty `scripts/` directory if this is its only file. |

Plus inline cleanup:

- Remove `(NOT YET DONE)` and `(not yet implemented)` markers anywhere
  they appear (anything still calling out `IdentifyRateEquationProblem` /
  `identify_rate_equation` as unimplemented).

Files left alone in this phase:

- `CLAUDE.md` — internal-facing; separate update pass.
- `docs/superpowers/specs/` and `docs/superpowers/plans/` — historical
  specs and plans; preserved.
- Anything in `src/` or `test/` — phase 2.
- `test/outcmaes/` — runtime output from CMA-ES, already in `.gitignore`,
  auto-regenerates each test run. Cosmetic to delete.
- `.CondaPkg/` — local conda env built by `OptimizationPyCMA` test dep,
  already in `.gitignore`. Deletion forces a multi-MB redownload with no
  permanent benefit.

### Verification

Before each deletion:

1. `git grep <basename>` and `git grep <full-path>` across the repo.
2. If references appear in a non-historical file (anything other than
   another `docs/superpowers/...` artifact or a previous spec referencing
   it by name), surface the reference and decide before deleting.
3. `git log --oneline -- <file>` to confirm the file isn't actively being
   modified.

## Acceptance criteria

The phase is done when:

- `README.md` exactly matches the 5-section structure above. No "Known
  Limitations", no API reference table, no `(NOT YET DONE)` text, no
  references to APIs that don't exist.
- `Pkg.test("EnzymeRates")` passes with the new `test_readme_runs.jl`
  exercising every non-skipped README code block in one shared module
  without errors.
- Files in the cleanup table are deleted; no dangling references remain.
- One commit (or a small number of focused commits) on the working branch
  containing all of the above. The commit messages cite this spec.

## Open questions deferred to implementation

- Exact optimizer choice for `fit_rate_equation` in section 2 — pick
  whichever gives a fast, deterministic-with-`Random.seed!` recovery for
  the example mechanism. Likely `LBFGSB` or `PyCMAOpt`.
- Exact numerical settings for `identify_rate_equation` in section 3 —
  whatever makes recovery succeed for the simulated data. Block is
  skip-marked, so CI runtime is not a constraint; only the *shown* call
  needs to be a sensible default a reader could realistically use.
- Whether to leave the example data table inline in the README or move it
  to a small helper that builds the simulated data from `rate_equation`.
  Inline is more readable; helper avoids 50 lines of `(group=…, S=…, …)`.
  Decide while writing — preference for inline if it fits.
