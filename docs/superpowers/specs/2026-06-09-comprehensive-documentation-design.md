# Comprehensive Documentation — Design

**Date:** 2026-06-09
**Branch:** `docs-comprehensive`
**Status:** Design approved; ready for implementation plan.

## Goal

Build a Documenter.jl documentation site for EnzymeRates.jl that faithfully
describes how the package works. The site becomes the canonical reference and
absorbs the user-facing and architecture content now scattered across the
README and CLAUDE.md. After the migration, the README shrinks to a landing
page and CLAUDE.md shrinks to agent-behavioral rules plus pointers into the
docs.

Three pillars structure the content: deriving rate equations, fitting them to
data, and identifying the best equation. Each pillar opens with a runnable
tutorial and continues with concept pages.

## Scope

One design (this document) and one phased implementation plan. The work spans
Documenter scaffolding, doctest and citation tooling, the three content
pillars, a developer/architecture page, an automated API page, a roadmap, and
the README/CLAUDE.md trim. The phasing keeps each step independently
reviewable and green.

## 1. Tooling and infrastructure

**Stack:** Documenter.jl + DocumenterCitations.jl.

**`docs/` subproject:** a `docs/Project.toml` that `dev`s EnzymeRates and adds
Documenter, DocumenterCitations, and the example dependencies (the optimizer
and any plotting package the tutorials need). These stay out of the main
`Project.toml` so the Aqua stale-dependency check stays green.

**`docs/make.jl`:**

```julia
makedocs(;
    sitename = "EnzymeRates.jl",
    doctest = true,
    checkdocs = :exports,
    plugins = [CitationBibliography("docs/src/refs.bib"; style = :authoryear)],
    format = Documenter.HTML(;
        canonical = "https://DenisTitovLab.github.io/EnzymeRates.jl",
        edit_link = "main",
    ),
)
deploydocs(; repo = "github.com/DenisTitovLab/EnzymeRates.jl.git", devbranch = "main")
```

`doctest = true` runs and checks the `jldoctest` blocks in docstrings.
`checkdocs = :exports` fails the build when an exported name lacks a
docstring. `devbranch = "main"` is set explicitly — DataDrivenEnzymeRateEqs.jl
leaves it commented, which defaults it to `"master"` and silently blocks the
`dev/` deploy on a `main`-default repo. We avoid that bug.

**CI:** add a `docs` job to the existing `.github/workflows/CI.yml` (not a
separate `Documentation.yml`), matching DataDrivenEnzymeRateEqs.jl and the
PkgTemplates default. The job uses `julia-actions/julia-docdeploy@v1`, reuses
the workflow-level `on:` block (push `main` + tags, `pull_request`,
`workflow_dispatch`), and grants `permissions: { contents: write }`.

**Auth:** pass both `DOCUMENTER_KEY` and `GITHUB_TOKEN` as env on the docdeploy
step, with `permissions: contents: write`. Denis has installed `DOCUMENTER_KEY`
(a write-access deploy key plus the matching Actions secret, both verified
present 2026-06-10), so Documenter authenticates via the SSH deploy key and
falls back to `GITHUB_TOKEN` if the key ever fails — the deploy succeeds either
way. The SSH key also makes **tagged-release** docs work: `TagBot.yml` already
passes `ssh: ${{ secrets.DOCUMENTER_KEY }}`, so a tag triggers the docs workflow
and builds the `stable/`/`vX.Y.Z/` versions. Denis flips the repo Pages source
to the `gh-pages` branch once the first CI deploy creates it (this serving step
is independent of the auth method).

**Build vs deploy:** `makedocs` builds on every trigger and turns a PR check
red on a broken page or failing doctest. `deploydocs` self-gates: it pushes to
gh-pages only on a push to `main` (into `dev/`) or a tag (into `vX.Y.Z/`, then
repoints `stable/`). PRs verify; they never publish.

## 2. Site structure

```
Home (index.md)            — what it is, install, ~10-line quickstart, links
Getting Started            — end-to-end example: define → derive → fit → identify (fast)

Deriving rate equations
├─ Rate equations from textbooks   — rate_equation_string tutorial
├─ Rapid equilibrium vs steady state
├─ The Cha / King–Altman algorithm — division-free, per-segment derivation
├─ Thermodynamic constraints       — Haldane / Wegscheider, parameter reduction
├─ Ping-pong mechanisms            — residual-on-:E representation
├─ Dead-end inhibitor binding      — the all-RE assumption
└─ MWC allostery                   — RE ligand binding & A/I transitions, 4-tag taxonomy

Fitting rate equations
├─ Fitting tutorial & data format
├─ Normalized vs absolute rate     — scale_k_to_kcat, kcat rescaling
└─ Loss & optimizers               — Optimization.jl choices (BYO optimizer)

Identifying the best rate equation
├─ Identify tutorial               — the fast width-1-beam example
├─ Model selection                 — LOOCV, paired 1-SE rule, permutation test
└─ The enumeration engine          — init mechanisms + six expansion moves

Developer / Architecture           — curated internals
API Reference                      — @autodocs (Private=false) + @index
Roadmap                            — KNF, identifiability, iso models, plotting, outlier detection
References                         — @bibliography
```

Each pillar opens with a runnable tutorial, then concept pages — a "do it,
then understand it" path. Getting Started stands apart from the pillars as a
single scannable arc for newcomers; it reuses the fast identify example so
nothing slow lands on the front path. Developer, API, Roadmap, and References
are reference material, not the learning path. Concept pages cross-link into
the relevant API entry with `[name](@ref)`.

Two settled choices: Getting Started is its own page (Home stays thin), and
concept pages are split fine-grained (one per topic) so each `@example`
sandbox stays small and pages deep-link cleanly.

## 3. Doctest and example strategy

**`jldoctest`** (run and output-checked) — only in docstrings of the
deterministic exported functions: `rate_equation_string`, `parameters`,
`metabolites`, `rate_equation`. These guard against equation and format drift.
Their output is byte-stable by construction; the audit confirmed
`rate_equation_string` produces identical bytes across separate cold Julia
processes.

**`@example`** (run, output captured, not checked) — every tutorial page, and
any block showing `fit_rate_equation` or `identify_rate_equation` output. Fit
results depend on a random multi-start and cannot be pinned; only their keys
and shapes are stable.

Doctest authoring rules from the audit:

- Show `m isa EnzymeMechanism`, never `typeof(m)` (the `Sig` type parameter is
  a huge unreadable string).
- Prefer `print(rate_equation_string(m))` over the bare return value to avoid
  quote and escape noise.
- For counts that vary with enumeration changes (mechanism totals), assert
  invariants in prose or `@example`, not pinned numbers.
- `_select_best_n_params` and `_onesided_permutation_p` are pure and seedable,
  so they *can* be `jldoctest`ed on hand-built inputs even though the full
  `identify_rate_equation` cannot.

**Caveat to record:** the `jldoctest` blocks assume `rate_equation_string`
output stays byte-stable across dependency bumps. The canonicalization
invariants make it deterministic for fixed code, but a dependency that
reorders term printing would force a `doctest = :fix` regeneration. Low risk,
easily managed.

## 4. The fast `identify_rate_equation` tutorial

The full search runs ~1 hour and cannot execute in CI. The tutorial instead
runs a real but fast search. These constraints are load-bearing, not
incidental:

- **Noiseless simulated data** and a **mechanism with no degenerate analogs**,
  so the width-1 greedy beam deterministically recovers the generating
  mechanism and the prose can claim recovery honestly.
- Beam arguments `loss_rel_threshold = 1.0`, `loss_abs_threshold = 0.0`,
  `min_beam_width = 1` collapse the beam to one survivor per parameter-count
  level. Verified against `_select_beam` (`src/identify_rate_equation.jl:317`):
  `cutoff = loss_rel_threshold * best + loss_abs_threshold = best`, so only the
  lowest-loss candidate clears it.
- A **small reaction** (uni-uni → few `init_mechanisms`) and a low
  **`max_param_count`** keep the initial level cheap.

The block runs in seconds as an `@example`. The prose describes the full
hour-long production search alongside, noting the default wider beam explores
more candidates.

## 5. Citations

DocumenterCitations.jl, **author-year** style, with inline `[Key](@cite)`
citations and a `@bibliography` block on the References page. The starting
reference set:

| Concept | Citation |
|---|---|
| Rate-equation derivation (RE + mixed SS) | Cha 1968, *J Biol Chem* |
| King–Altman pattern method | King & Altman 1956, *J Phys Chem* |
| MWC allostery | Monod, Wyman & Changeux 1965, *J Mol Biol* |
| Enzyme-kinetics reference text | Segel 1975, *Enzyme Kinetics* |
| Thermodynamic constraints | Haldane 1930 / Wegscheider 1901 |
| Model selection (1-SE rule) | Hastie, Tibshirani & Friedman, *ESL* |
| Optimizer | Hansen, CMA-ES |

More references can be added later. Gathering exact DOIs and correct BibTeX is
a verification step in Phase 0 — a small workflow fetches and confirms each
entry rather than guessing.

## 6. API page

A consolidated **API Reference** page built with `@autodocs` (`Modules =
[EnzymeRates]`, `Private = false`), topped with `@index`. Adding a new export
later surfaces it with no edit to the page. `checkdocs = :exports` makes the
build fail on any undocumented export — the real guarantee that the page stays
current.

**Docstring gaps** (the only two exports that fail `checkdocs = :exports`):

| Export | Problem | Fix |
|---|---|---|
| `EnzymeReaction` | `#`-comment, no docstring (`src/types.jl:272-276`) | Write a type docstring; convert and expand the comment. |
| `metabolites` | Docstring detached by an intervening `#`-comment (`src/types.jl:1167-1176`) | Move the comment below the function signature so the docstring attaches. |

Several attached docstrings are terse and need expansion before they read as
documentation — `rate_equation_string` most of all (it documents neither the
`mode` argument, the multi-line format, nor the Full-mode restriction).
Expansion is content work, not a `checkdocs` failure.

## 7. CLAUDE.md trim and README disposition

**README** → thin landing page: badges, a one-paragraph description, install,
a ~10-line quickstart, and links into the docs. Retire
`test/test_readme_runs.jl`; its coverage moves into Documenter doctests, a net
increase. The README's stale prose disappears with the trim rather than being
corrected in place.

**CLAUDE.md** → splits by audience. When a block migrates, rewrite it
correctly from the audit (approved: fix stale references as part of this
pass). Dual-nature blocks keep a one-line guard in CLAUDE.md and move their
explanation to the docs.

### Migration table

**Stay in a slimmed CLAUDE.md (agent-behavioral rules):**
Foundational rules, Our relationship, Proactiveness, Designing software, TDD,
Writing code, Naming, Code Comments, Version Control, Testing, Issue tracking,
Systematic Debugging, Learning and Memory Management, Commands, Workflow, Code
Style, and the `rate_equation` "runtime perf is non-negotiable" block
(cross-linked from the Developer page).

**Move to the Developer page (maintainer-only):**
the `name(p, m)` parameter-naming chokepoint and its AST-walker test (correct
the reference: enforcement lives in `test/test_types.jl:1577-1644`, not a
`test_chokepoint.jl`); Canonical Step Form function names
(`_canonical_group_order!`, `_step_canonical_key`, `_entry_kind`,
`_canonical_iso_direction`); the `EnzymeMechanism{Sig}` / `_sig_of` /
`Mechanism(em)` lift and the concrete-vs-singleton split; the
`AllostericEnzymeMechanism` three-type-parameter rationale; Source Layout; the
compile-budget / `@generated` Known Issues note.

**Move to how-it-works pages (user concepts), corrected:**
Allosteric state taxonomy (standardize on A/I); Regulator representation and
`@enzyme_reaction` DSL roles (rewrite — the current block is almost entirely
stale); Mechanism enumeration building blocks (six moves, flat return);
Vmax/kcat normalization (correct the function names); Catalytic topology
constraints C1–C10 (chemistry is user-facing, function names dev-only);
Haldane/Wegscheider (add the exact-integer null space, the Haldane-vs-
Wegscheider distinction, and the division-free derivation story).

The governing principle: **CLAUDE.md keeps the rules for how Claude works in
this repo; the docs own what the package does and how it is architected.**

## 8. Staleness ledger (fix during migration)

The audit found references in CLAUDE.md and the README that no longer match
the code. Fix each as its content migrates. The high-impact items:

**Phantom or renamed names in CLAUDE.md:**

- `_is_ss_rate_constant` → `_ss_rate_constant_names` (classifies by `Parameter`
  subtype, `src/rate_eq_derivation.jl:621`), not "lowercase k + digit".
- `_kcat_components` → `_kcat_groups_from_polys`
  (`src/rate_eq_derivation.jl:650`).
- `test/test_chokepoint.jl` → does not exist; enforcement is
  `test/test_types.jl:1577-1644`.
- `regulator_roles` → does not exist. `regulators(rxn::EnzymeReaction)` returns
  `Vector{RegulatorMults}` (`src/types.jl:339`); bare Symbols come from the
  `regulators(em::EnzymeMechanism)` overload (`src/types.jl:1162`).
- `@enzyme_reaction` role labels (`regulators:`, `:unknown`, `:dead_end`) →
  wrong. `dead_end_inhibitors:` and `competitive_inhibitors:` both map to
  `:competitive`; roles are `:competitive` and `:allosteric` only
  (`src/dsl.jl:39-44`).
- `has_residual` in `backtrack!` → the parameter is
  `pingpong_intermediate::Bool` (`src/mechanism_enumeration.jl:328`);
  `has_residual` is a `Species` accessor.
- The "covalent intermediate always on `:E`" claim holds for the enumerator,
  not the DSL (which allows separate labels such as `F`). Scope the claim.

**README drift:**

- `rate_equation_string` "prints" → it **returns a String** (R:66, R:175,
  R:235).
- Model selection (R:285–295) describes plain argmin-CV; the real rule is the
  smallest bucket below `n_min` passing both the paired 1-SE test
  (`se_threshold = 1.0`) and the one-sided permutation test
  (`perm_p_threshold = 0.16`). The CV score is the **mean of log** per-fold
  losses (`src/identify_rate_equation.jl:867`).
- MWC names `K_R`/`K_T` → real names use A/I tokens: `K_A_S_E`, `K_I_S_E`;
  `:EqualAI` renders with no token. Catalytic `:OnlyI` is a hard constructor
  error.
- Enumeration: `dedup!` → `_dedup_flat!` (= `unique!`); canonicalization is in
  the constructor, not dedup; `expand_mechanisms` returns a **flat** vector;
  there are **six** moves (the README lists five, omitting
  `_expand_add_allosteric_regulator`); `_expand_re_to_ss` flips a whole kinetic
  group atomically, not one step.

The audit's per-topic briefs (15 of them) and the cross-cutting critic report
are the authoritative source for content; reference them while writing each
page.

## 9. Coverage gaps to add

Behaviors worth documenting that the outline omits:

1. `IdentifyRateEquationResults` struct accessors (`best`, `cv_results`
   schema, frontier).
2. The CSV / `progress.log` artifacts `identify_rate_equation` writes
   (`save_dir` is mandatory; cluster-visible, resumable).
3. `max_param_count` — caps actual fitted params, bounds search depth.
4. `@enzyme_reaction` atom-bracket syntax (`S[C]`, `A[C1H1]`) — load-bearing
   for ping-pong and multi-substrate reactions.
5. The `::Inh` role tag — bind a real metabolite in `CompetitiveInhibitor`
   role while keeping its name so `concs.NAME` drives it.
6. `FitFailure` and loud-failure semantics — exception text captured,
   `_loocv` raises on a non-finite fold, CSV rows carry `retcode`/`error`.
7. `rate_equation` as a first-class exported call with its 0-allocation /
   sub-100 ns contract.
8. `oligomeric_state:` vs `allowed_catalytic_multiplicities:` vs
   `catalytic_multiplicity:` — three related knobs and a known footgun.

## 10. Resolved decisions and defaults

Approved this session:

- **CLAUDE.md staleness:** fix stale references as content migrates.
- **Enumeration framing:** present moves as "mostly +1 parameter," but state
  that a Haldane/Wegscheider constraint can absorb the parameter (net **+0**)
  and that changing a steady-state group from `EqualAI` to `NonequalAI` adds
  **+2** (independent `kf` and `kr` per state). This is why the search re-fits
  and buckets by actual fitted-param count rather than assuming +1. Mark the
  exact per-move deltas **verify-during-writing**.
- **MWC terminology:** standardize on **A/I** (active/inactive), matching the
  parameter names users type; state the R≡A, T≡I correspondence once.

Defaults, overridable at spec review:

- Textbook tutorial uses the **RE form** (clean `1 + S/K_S_E + P/K_P_E`
  denominator), which reads like the literature.
- Provide a **symbol-mapping note** (package `K_S_E`, `k_ES_to_EP` ↔ textbook
  `Km`, `kcat`, `kf`, `kr`).
- Show the allosteric Haldane-RHS rendering quirk **as-is**, annotated as a
  known rendering quirk (bare active-state base names on the RHS while `params`
  lists A/I-suffixed names — confirmed real behavior, do not "fix").
- Present the **relative** path (`scale_k_to_kcat::Real`, default) as primary
  and the **absolute** turnover path (`nothing`) as advanced.
- Bless **`PyCMAOpt()`** (`OptimizationPyCMA`) with
  `BBO_adaptive_de_rand_1_bin_radiuslimited` as a tested alternative; state
  that the optimizer is **bring-your-own** (the base package depends only on
  `Optimization`).
- Treat `rescale_parameter_values` as the public kcat helper; `_kcat_forward`
  and friends stay internal.
- Add a type docstring for `AllostericEnzymeMechanism`.

## 11. Phasing

Each phase ends green (tests pass, docs build passes).

- **Phase 0 — Infra scaffold.** `docs/` subproject; `make.jl`; `refs.bib`
  (DOIs gathered via a verification workflow); `docs` job in `CI.yml`
  (`DOCUMENTER_KEY` + `GITHUB_TOKEN`); a minimal Home plus empty section stubs.
  After the first CI deploy creates `gh-pages`, Denis sets the repo Pages source
  to that branch. **Exit:** docs build green locally and in CI; first gh-pages
  deploy eyeballed.
- **Phase 1 — Docstring gaps + API page.** Write the `EnzymeReaction`
  docstring; fix the `metabolites` detachment; expand the terse deterministic-
  function docstrings with their `jldoctest` blocks; build the consolidated API
  Reference page. **Exit:** `checkdocs = :exports` passes.
- **Phase 2 — Derivation pillar.** The seven derivation pages, sourced from the
  audit briefs; `jldoctest` where deterministic.
- **Phase 3 — Fitting pillar.** Tutorial + data format; normalized vs absolute;
  loss + optimizers (BYO note).
- **Phase 4 — Identify pillar.** Fast-example tutorial; model selection (1-SE +
  permutation); the enumeration engine (six moves + actual-count framing); plus
  coverage gaps (results struct, CSV artifacts, `max_param_count`).
- **Phase 5 — Developer page + Getting Started + Roadmap.** Migrate maintainer
  internals (corrected per audit); write the end-to-end Getting Started; write
  the Roadmap.
- **Phase 6 — Trim and slim.** Apply the CLAUDE.md migration table (cut, fix
  stale references, leave guard lines); slim the README to a landing page;
  retire `test/test_readme_runs.jl`. **Exit:** full test suite green, docs
  build green, all doctests pass.

Phase 0 proves the CI/deploy/doctest pipeline before any content exists. The
trim is last so the README and CLAUDE.md point only at pages that already
exist.

## Verification

- `make.jl` builds locally without warnings.
- The `docs` CI job is green on a PR and deploys on merge.
- `checkdocs = :exports` passes (every export documented).
- Every `jldoctest` passes; every `@example` runs.
- The fast identify tutorial runs in seconds and recovers its generating
  mechanism.
- The full test suite stays green after the README/CLAUDE.md trim.
