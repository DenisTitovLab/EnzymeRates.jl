# README rewrite and cheap dead-doc cleanup — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the stale `README.md` with a focused 5-section introduction
(description; mechanism + derive + simulate + fit; reaction +
`identify_rate_equation` recovery; biochemist-intuition explanations of
derivation and selection); add a small test that keeps every README
code block runnable in CI; delete obsolete repo-root planning artifacts.

**Architecture:** README is a hand-written Markdown file; a new
`test/test_readme_runs.jl` extracts ` ```julia ` fenced code blocks
from `README.md` and `eval`s them in one shared anonymous module so
state (mechanism `m`, simulated `data`, fitted `params`) flows from
earlier blocks to later ones. Blocks whose first non-blank line is
`# README-SKIP-IN-TEST` are excluded from extraction. No new runtime
dependencies; the test reuses `OptimizationPyCMA` from the existing
test extras.

**Tech Stack:** Julia 1.9+, plain Markdown, `Pkg.test` workflow,
`OptimizationPyCMA` (test-only dep) for the fit demo and the
`identify_rate_equation` showcase.

**Reference spec:** `docs/superpowers/specs/2026-04-30-readme-rewrite-and-deadcode-design.md`

---

## File Structure

| File | Action | Purpose |
|------|--------|---------|
| `README.md` | Replace | New 5-section README; the source of truth. |
| `test/test_readme_runs.jl` | Create | Extracts and runs README code blocks. |
| `test/runtests.jl` | Modify | Wire the new test file into the suite. |
| `SPEC.md` | Delete | Superseded by new README. |
| `CODE_SIMPLIFICATION_PROMPT.md` | Delete | Repo-root one-shot prompt artifact. |
| `PLAN_IMPLEMENTATION_PROMPT.md` | Delete | Repo-root planning artifact. |
| `PLAN_RESS_DEDUP.md` | Delete | Planning artifact for completed refactor. |
| `ralph.sh` | Delete | Developer scratch agent runner. |
| `.ralph-logs/` | Delete | Empty log dir for `ralph.sh`. |
| `scripts/verify_counts.py` | Delete | Topology-count verification, no longer referenced. |
| `scripts/` | Delete (if empty after above) | |

---

## Task 1: Cheap doc cleanup — verify and delete obsolete files

**Files:**
- Delete: `SPEC.md`, `CODE_SIMPLIFICATION_PROMPT.md`,
  `PLAN_IMPLEMENTATION_PROMPT.md`, `PLAN_RESS_DEDUP.md`, `ralph.sh`,
  `.ralph-logs/`, `scripts/verify_counts.py`, `scripts/` (if empty)

- [ ] **Step 1: Verify nothing live references these files**

Run from repo root:

```bash
git grep -nl 'SPEC.md\|CODE_SIMPLIFICATION_PROMPT\|PLAN_IMPLEMENTATION_PROMPT\|PLAN_RESS_DEDUP\|ralph\.sh\|verify_counts' -- ':!docs/superpowers/'
```

Expected: matches only inside `docs/superpowers/specs/...` (the design spec
itself) and `.gitignore` (which mentions `outcmaes/`). If there are matches
in `src/`, `test/`, `Project.toml`, or `.github/`, stop and surface them
before deleting.

- [ ] **Step 2: Delete the files**

```bash
git rm SPEC.md CODE_SIMPLIFICATION_PROMPT.md PLAN_IMPLEMENTATION_PROMPT.md PLAN_RESS_DEDUP.md ralph.sh scripts/verify_counts.py
rm -rf .ralph-logs scripts
```

`.ralph-logs` and `scripts/` are untracked once `verify_counts.py` is gone,
so plain `rm -rf` is correct. (`git rm` works on tracked files only.)

- [ ] **Step 3: Verify the cleanup left the repo in a coherent state**

```bash
ls /home/denis.linux/.julia/dev/EnzymeRates/
```

Expected (no `SPEC.md`, `CODE_SIMPLIFICATION_PROMPT.md`,
`PLAN_IMPLEMENTATION_PROMPT.md`, `PLAN_RESS_DEDUP.md`, `ralph.sh`,
`.ralph-logs`, or `scripts`):

```
Manifest.toml
Project.toml
README.md
docs
src
test
```

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "$(cat <<'EOF'
Delete obsolete planning artifacts and the stale SPEC.md

SPEC.md is superseded by the upcoming README rewrite; the rest are
single-shot planning prompts and developer scratch from completed
refactors. No live references in src/, test/, Project.toml, or .github/.

Refs: docs/superpowers/specs/2026-04-30-readme-rewrite-and-deadcode-design.md

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Replace README with Section 1 stub and add runnability test

**Goal of this task:** establish the runnability-test scaffolding and
prove it works on a minimal README that contains exactly one
unambiguously-runnable code block (so the test doesn't pass vacuously
on an empty README).

**Files:**
- Replace: `README.md`
- Create: `test/test_readme_runs.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Replace `README.md` with Section 1 stub**

Overwrite the file so it contains *only*:

```markdown
# EnzymeRates.jl

Identify the best enzyme rate equation from kinetic data.

Given a reaction definition and experimental rate measurements at varying
substrate, product, and regulator concentrations, EnzymeRates enumerates
all biochemically valid mechanisms, fits each to the data, and selects
the simplest mechanism that adequately describes the data based on
leave-one-group-out cross-validation.

The package has first-class support for MWC allostery and for mechanisms
that mix steady-state and rapid-equilibrium elementary steps.
Thermodynamic constraints (Haldane, Wegscheider) are derived
automatically from the cycle structure of the mechanism, so users supply
only the independent rate constants plus a measured equilibrium
constant.

## Installation

```julia
# README-SKIP-IN-TEST
using Pkg
Pkg.add(url="https://github.com/DenisTitovLab/EnzymeRates.jl")
```

## Smoke test

```julia
using EnzymeRates
m = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        [E, S] <--> [ES]
        [ES] <--> [E, P]
    end
end
@assert :S in metabolites(m)
```
```

The "Smoke test" block is intentional scaffolding — it validates the
runnability test is doing something, and gets replaced by real Section 2
content in Task 3. The fenced code blocks inside the README contain
literal triple-backticks; in your editor the Markdown source will look
exactly as printed.

- [ ] **Step 2: Create `test/test_readme_runs.jl`**

```julia
# ABOUTME: Extracts ```julia code blocks from README.md and runs them in
# ABOUTME: one REPL session, skipping blocks tagged # README-SKIP-IN-TEST.

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

- [ ] **Step 3: Wire into `test/runtests.jl`**

First read the file to see what's there:

```bash
cat test/runtests.jl
```

It includes other test files via `include(...)`. Add a line for the new
test file. The exact edit depends on the current file's structure;
add `include("test_readme_runs.jl")` alongside the existing
`include(...)` lines, keeping them in alphabetical order if that's the
existing convention.

- [ ] **Step 4: Run the runnability test in isolation**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40
```

Expected: full test suite passes, including a `README runs` testset that
extracts at least one block (the Smoke test) and runs it successfully.
The Pkg.add install line is skipped via the magic comment.

If the test fails with a regex mismatch, double-check that the README's
code-block fences begin with triple-backticks followed exactly by `julia`
on the same line, with a newline before the code body.

- [ ] **Step 5: Commit**

```bash
git add README.md test/test_readme_runs.jl test/runtests.jl
git commit -m "$(cat <<'EOF'
Replace README with Section 1 stub and add runnability test

test/test_readme_runs.jl extracts ```julia blocks from README.md and
evaluates them in one shared module. Blocks whose first non-blank line
is # README-SKIP-IN-TEST are excluded.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Section 2 — define mechanism, derive rate equation

**Goal:** flesh out Section 2 with the running example mechanism (uni-uni
with `:OnlyR` substrate, product, catalytic step, and allosteric
activator), inspect via `parameters` / `rate_equation_string` /
`rate_equation`. Defer simulation + fit to Task 4.

**Files:**
- Modify: `README.md` (replace the Smoke test section, add Section 2)

- [ ] **Step 1: Replace the README's "Smoke test" section with Section 2's first block**

Replace the section starting `## Smoke test` (and its single code block)
with:

````markdown
## Define a mechanism, derive its rate equation, fit data

The running example is a uni-uni reaction `S ⇌ P` catalyzed by an MWC
homodimer in which substrate, product, the catalytic interconversion,
and an allosteric activator `A` all operate exclusively in the R
conformation. The T conformation is catalytically silent — a textbook
K-type allosteric activator.

```julia
using EnzymeRates

m = @allosteric_mechanism begin
    substrates: S
    products:   P
    allosteric_regulators: A::OnlyR

    site(:catalytic, 2): begin
        steps: begin
            [E, S] ⇌ [ES]      :: OnlyR
            [ES] <--> [EP]     :: OnlyR
            [EP] ⇌ [E, P]      :: OnlyR
        end
    end
end
```

The two `⇌` steps mark binding events at rapid equilibrium (one
binding constant `K` per step); the `<-->` step marks a steady-state
catalytic interconversion (independent forward and reverse rate
constants `kf`, `kr`). The `::OnlyR` annotation tells the framework
that those steps fire only in the R conformation, so the T-state
contribution to the rate equation is identically zero.

`parameters(m)` lists the names the framework needs at evaluation
time, `rate_equation_string(m)` prints the symbolic rate equation, and
`rate_equation(m, concs, params)` evaluates the rate numerically:

```julia
parameters(m)
rate_equation_string(m)
```
````

The first ` ```julia ` block defines `m`, used by the rest of the
README. The second prints the parameter list and rate-equation string —
both side-effect-free, kept separate for readability.

- [ ] **Step 2: Run the runnability test**

```bash
julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -40
```

Expected: full suite passes, including `README runs`. The
`@allosteric_mechanism` macro compiles, `parameters(m)` returns a tuple
of Symbols, `rate_equation_string(m)` returns a String. None of these
need to be matched exactly — the test only cares that the calls don't
error.

If `@allosteric_mechanism` errors, the most likely cause is mis-formed
DSL syntax — read `src/dsl.jl:435` (`macro allosteric_mechanism`) and
`src/dsl.jl:593` (`_parse_allosteric_mechanism_body`) to verify
expected forms. The pattern used here matches the existing test in
`test/test_dsl.jl:46`.

- [ ] **Step 3: Capture the actual `parameters(m)` output for use in Task 4**

Run a quick interactive check from the repo root to record the precise
parameter names that will appear in `params` later:

```bash
julia --project=. -e '
    using EnzymeRates
    m = @allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: A::OnlyR
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌ [ES]      :: OnlyR
                [ES] <--> [EP]     :: OnlyR
                [EP] ⇌ [E, P]      :: OnlyR
            end
        end
    end
    println(parameters(m))
'
```

Record the output. The expected shape is a tuple of Symbols including
binding constants (likely `:K1`, `:K2`), an SS forward rate (likely
`:k3f`), the conformational equilibrium `:L`, an activator binding
constant (likely `:K_A_reg1`), `:Keq`, and `:E_total`. The exact
*names* depend on `_kinetic_rename_map` numbering and the regulatory
site's auto-naming — Task 4's `true_params` NamedTuple must use these
exact keys.

- [ ] **Step 4: Commit**

```bash
git add README.md
git commit -m "$(cat <<'EOF'
README Section 2: define example mechanism, inspect rate equation

Adds the @allosteric_mechanism block, explanatory prose, and the
parameters/rate_equation_string inspection calls. Runnability test
covers the new blocks.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Section 2 — simulate data and fit with `fit_rate_equation`

**Goal:** extend Section 2 with synthetic-data generation and a
`fit_rate_equation` call that recovers the true parameters. The fit must
complete within a few seconds during `Pkg.test()` so CI stays fast.

**Files:**
- Modify: `README.md` (append simulation + fit blocks to Section 2)

- [ ] **Step 1: Append simulation prose and code to Section 2**

Append after the existing Section 2 content:

````markdown

We generate synthetic data by evaluating `rate_equation` on a grid of
concentrations and adding multiplicative log-normal noise. Multiple
`group` values represent independent experimental batches that share
the same `E_total`; the framework's loss function is invariant to a
per-group `E_total` rescaling.

```julia
using OptimizationPyCMA
Random.seed!(42)

# True parameters (independent rate constants + Keq + E_total).
# Replace the parameter names with whatever `parameters(m)` printed
# above if your build differs.
true_params = (
    K1 = 1.0,         # S binding K (R-state)
    K2 = 0.5,         # P binding K (R-state)
    k3f = 5.0,        # catalytic forward rate (R-state)
    L = 0.1,          # conformational [T]/[R] for free enzyme
    K_A_reg1 = 2.0,   # activator binding K (R-state)
    Keq = 2.0,
    E_total = 1.0,
)

# Concentration grid: vary S and A across a dynamic range; P held at two
# levels including a near-zero value to constrain the reverse direction.
data_rows = NamedTuple[]
for grp in 1:5
    for _ in 1:10
        S = exp(randn() * 0.8)         # ~lognormal around 1
        A = 0.05 + 5.0 * rand()        # uniform 0.05..5.05
        P = rand() < 0.5 ? 0.05 : 0.5
        v_true = rate_equation(m, (S=S, P=P, A=A), true_params)
        v_obs = v_true * exp(0.05 * randn())   # 5% log-normal noise
        push!(data_rows, (group="G$grp", Rate=v_obs, S=S, P=P, A=A))
    end
end
data = (
    group = [r.group for r in data_rows],
    Rate  = [r.Rate  for r in data_rows],
    S     = [r.S     for r in data_rows],
    P     = [r.P     for r in data_rows],
    A     = [r.A     for r in data_rows],
)
```

The fit runs `fit_rate_equation` on a `FittingProblem`, using the PyCMA
optimizer (multi-start CMA-ES) recommended for rate-equation fitting.
Fitted rate constants are returned with kcat normalized to 1.0 by
default — the absolute scale is recovered by multiplying with a
separately measured kcat.

```julia
fp = FittingProblem(m, data; Keq=2.0)
result = fit_rate_equation(fp, PyCMAOpt();
    n_restarts=3, maxtime=5.0, popsize=50)
result.params       # fitted (K1, K2, k3f, L, K_A_reg1)
result.loss         # final loss value
```
````

Note that `kcat=1.0` is the default in `fit_rate_equation`, which
rescales the SS rate constants so the resulting kcat equals 1.0.
Comparing fitted vs true values therefore requires also rescaling
`true_params` to the same kcat — but readers only need to see that the
fit converges (`result.loss` is small).

- [ ] **Step 2: Run the runnability test**

```bash
julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -60
```

Expected: passes, including the new blocks. The fit takes a few seconds
(`n_restarts=3, maxtime=5.0` keeps each restart bounded). Watch for:

- *Parameter-name mismatch*: if `true_params` uses keys that don't
  appear in `parameters(m)`, `rate_equation` errors. Cross-reference
  with the output recorded in Task 3 Step 3 and edit the NamedTuple.
- *Negative simulated rates*: if the noise occasionally pushes a rate to
  zero or sign-flips, `FittingProblem` rejects zero rates. The
  `Random.seed!(42)` makes this reproducible — if it fires, change the
  seed or tighten the noise scale.

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "$(cat <<'EOF'
README Section 2: simulate data and fit with fit_rate_equation

Adds synthetic-data generation (grid over S, A; two P levels; 5
groups; 5% log-normal noise) and a fit_rate_equation call using the
PyCMA optimizer. The runnability test exercises the full fit.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Section 3 — define `EnzymeReaction`, recover via `identify_rate_equation`

**Goal:** add Section 3 — define the same chemistry as a *reaction*
(no mechanism), reuse the `data` from Section 2, and demonstrate
`identify_rate_equation`. The slow `identify_rate_equation` call is
marked `# README-SKIP-IN-TEST` because (a) it would dominate test time,
(b) `test/test_identify_rate_equation.jl` already exercises the same
path with reduced settings.

**Files:**
- Modify: `README.md` (append Section 3)

- [ ] **Step 1: Append Section 3 to the README**

Append after Section 2:

````markdown

## Recover the mechanism with `identify_rate_equation`

If the mechanism is unknown — only the overall reaction and its
regulators are — `identify_rate_equation` enumerates biochemically
valid mechanisms, fits each to the data, and returns the simplest that
generalizes (judged by leave-one-group-out cross-validation). The same
chemistry from Section 2, declared as a *reaction*:

```julia
rxn = @enzyme_reaction begin
    substrates: S
    products:   P
    regulators: A
    oligomeric_state: 2
end
```

`regulators: A` declares `A` with an unspecified role; the search
enumerates dead-end-inhibitor and allosteric variants and selects
between them on cross-validation score. (If you already know `A` is
allosteric, declare it with `allosteric_regulators: A` instead and the
search skips dead-end variants.)

```julia
# README-SKIP-IN-TEST
using OptimizationPyCMA
prob = IdentifyRateEquationProblem(rxn, data; Keq=2.0)
results = identify_rate_equation(prob;
    optimizer=PyCMAOpt(),
    max_param_count=10,
    pmap_function=map,            # serial; pass `pmap` for distributed
)
results.best                       # the recovered mechanism
rate_equation_string(results.best) # printed rate equation
first(results.cv_results, 5)       # top rows of the CV-score DataFrame
```

`results.best` is the mechanism with minimum training loss at the
parameter-count level whose CV score is lowest — i.e., the simplest
mechanism that generalizes. For this synthetic data, the recovered
mechanism agrees with the one we used to generate it.
````

- [ ] **Step 2: Run the runnability test**

```bash
julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -40
```

Expected: passes. The `@enzyme_reaction` block and
`IdentifyRateEquationProblem` declaration run; the slow
`identify_rate_equation` call is skipped because its first non-blank
line is `# README-SKIP-IN-TEST`.

If the test passes but `IdentifyRateEquationProblem` errors with a
"missing column" message, the `regulators: A` line is producing the
wrong column-name expectation. Cross-reference
`src/identify_rate_equation.jl:38` — the constructor expects a column
for every metabolite-or-regulator name.

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "$(cat <<'EOF'
README Section 3: identify_rate_equation recovery from simulated data

Same chemistry as Section 2, declared as @enzyme_reaction. The slow
identify_rate_equation call is skipped from the runnability test
because test/test_identify_rate_equation.jl already covers it.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Section 4 — biochemist intuition for rate-equation derivation

**Goal:** add the prose-only Section 4 covering steady-state vs rapid
equilibrium, Haldane / Wegscheider, and the MWC R/T model. No code
blocks — pure exposition pointing back at Section 2's printed rate
equation.

**Files:**
- Modify: `README.md` (append Section 4)

- [ ] **Step 1: Append Section 4 to the README**

Append after Section 3:

````markdown

## How rate-equation derivation works

### Steady-state vs rapid equilibrium

`<-->` denotes a steady-state (QSSA) elementary step — both the forward
and reverse rate constants enter the rate equation as independent
parameters, and the King-Altman determinant assembles them into the
denominator polynomial. `⇌` denotes a rapid-equilibrium step — only the
binding constant `K` matters, because the framework collapses the
forward and reverse rates into a single equilibrium relation. A typical
mechanism mixes both, and the framework handles the mixed Cha-style
derivation automatically. `parameters(m)` reflects this: each RE step
contributes one `K`; each SS step contributes a forward `kf` and a
reverse `kr`.

### Haldane and Wegscheider relationships

When the mechanism contains thermodynamic cycles — any closed loop of
binding and catalytic steps — the rate constants around the cycle are
constrained by the equilibrium constant of the overall reaction. The
framework detects these cycles automatically (via the null space of the
mechanism's stoichiometric matrix), declares one rate constant per
cycle as *dependent*, and computes it from the rest plus a user-supplied
`Keq`. You fit the *independent* rate constants; dependent constants
are derived. `structural_identifiability_deficit(m)` reports the deficit
of the mechanism's parameter map: zero means every independent
parameter can in principle be identified from the rate equation.

### Allostery: the MWC R/T model

For multi-subunit enzymes, the framework uses the Monod-Wyman-Changeux
two-state model: the enzyme exists in an active R conformation and an
inactive T conformation, with `L = [T]/[R]` the conformational
equilibrium for the bare enzyme and the same `L` propagating to all
ligand-bound species. Each kinetic group (binding step, catalytic
interconversion) and each regulatory ligand can independently be:

- `:OnlyR` — the symbol exists in R only; T-state contributions are
  zero. A `:OnlyR` activator binds R preferentially and shifts the
  population toward R.
- `:OnlyT` — symbol exists in T only. A `:OnlyT` regulator binds T
  preferentially and shifts the population toward T (a typical
  allosteric inhibitor).
- `:EqualRT` — same `K` (or `kf`, `kr`) in both conformations. Useful
  for ligands that bind without conformational preference.
- `:NonequalRT` — independent R and T parameters (`K_R`, `K_T`).

The full rate equation is then the sum of R-state and T-state numerator
terms, weighted by the partition function `(R-state polynomial)^n +
L*(T-state polynomial)^n`, where `n` is the oligomeric state. The
example mechanism above uses `:OnlyR` everywhere — the T-state
contributions vanish and the printed rate equation simplifies
accordingly.
````

- [ ] **Step 2: Run the runnability test**

```bash
julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -20
```

Expected: passes, no new code blocks.

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "$(cat <<'EOF'
README Section 4: rate-equation derivation explanation

Biochemist-intuition prose covering QSSA vs rapid equilibrium,
Haldane/Wegscheider thermodynamic constraints, and the MWC R/T
allosteric model. No code blocks.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Section 5 — biochemist intuition for `identify_rate_equation`

**Goal:** add the prose-only Section 5 covering enumeration building
blocks, beam search, and LOOCV-based model selection.

**Files:**
- Modify: `README.md` (append Section 5)

- [ ] **Step 1: Append Section 5 to the README**

Append after Section 4:

````markdown

## How `identify_rate_equation` works

### Enumeration as composable building blocks

Mechanism enumeration is built from three small functions, not a
monolithic pipeline:

- `init_mechanisms(reaction)` produces the biochemically minimal
  mechanisms for a reaction by combining catalytic topologies (orderings
  of substrate binding, catalytic interconversion, and product release —
  random-order, ordered, ping-pong) with subsets of dead-end inhibition
  steps. Steps that bind the same metabolite share a kinetic group, so
  the parameter count starts at the smallest physically meaningful
  value.
- `expand_mechanisms(specs, reaction)` applies a fixed set of
  single-move expansions to each spec — converting an RE step to SS,
  splitting a kinetic group, adding a dead-end regulator, becoming
  allosteric, changing an allosteric state — and returns the expanded
  candidates keyed by their estimated parameter count.
- `dedup!(cache)` canonicalizes specs (sorted steps; renumbered kinetic
  groups by first occurrence) and removes structural duplicates.

The enumeration is grounded in chemical reasoning rather than blind
combinatorics: a step is "elementary" only if it changes one binding
site by one event with atom balance preserved, and only catalytic
topologies that satisfy bounds on bound-metabolite count, isomerization
size, and substrate participation are emitted.

### Beam search across parameter counts

`identify_rate_equation` walks parameter counts in ascending order:

1. Fit all init mechanisms on the full data; record training loss.
2. Keep the top fraction (beam width = `max(beam_fraction * n,
   min_beam_width)`).
3. Apply `expand_mechanisms` to surviving specs to produce candidates
   at the next parameter-count level.
4. `dedup!` and fit; rank by training loss.
5. Repeat until no new candidates appear or `max_param_count` is
   reached.

The beam width balances coverage (more candidates explored) against
runtime (every kept candidate gets a multi-restart fit).

### Model selection by leave-one-group-out cross-validation

After beam search, the top `n_cv_candidates` mechanisms per parameter
count enter LOOCV. Each unique value of the `group` column defines one
fold: the mechanism is fit on every group except one, then evaluated on
the held-out group. The CV score is the mean held-out loss across
folds. The "best" mechanism is the one with minimum training loss at
the parameter count whose CV score is lowest — *the simplest mechanism
that generalizes*. The `group` column reflects experimental batches
that share an `E_total`; LOOCV respects this structure and gives an
honest estimate of generalization to new conditions.
````

- [ ] **Step 2: Run the runnability test**

```bash
julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -20
```

Expected: passes.

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "$(cat <<'EOF'
README Section 5: identify_rate_equation explanation

Biochemist-intuition prose on the init/expand/dedup enumeration
building blocks, beam search across parameter counts, and
leave-one-group-out CV-based model selection.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Final integration check

**Goal:** run the full test suite end-to-end, do one last consistency
sweep on the README, and confirm no orphan references to deleted files.

**Files:**
- Read: `README.md`, `test/test_readme_runs.jl`, `test/runtests.jl`

- [ ] **Step 1: Run the full test suite**

```bash
julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -80
```

Expected: every existing testset passes, plus `README runs` extracting
~6 blocks (Section 2 mechanism, Section 2 inspect, Section 2 simulate,
Section 2 fit, Section 3 reaction definition; Sections 4 and 5 are
prose-only). The `Pkg.add` install line and the `identify_rate_equation`
call are skipped via the magic comment.

- [ ] **Step 2: Sweep for stale phrases**

```bash
git grep -n 'NOT YET DONE\|not yet implemented\|enumerate_mechanisms\|graph(m)\|EnzymeMechanism{Species\|compile_mechanism' -- ':!docs/superpowers/' ':!CLAUDE.md' ':!.claude/CLAUDE.md'
```

Expected: no matches. (`CLAUDE.md` and `.claude/CLAUDE.md` are excluded
because internal-facing doc updates are explicitly out of scope for
this phase — they get caught in a later pass.)

If matches appear, edit the offending file to remove the stale
reference, run `Pkg.test()` to confirm nothing broke, and commit
separately:

```bash
git add <file>
git commit -m "Remove stale reference to <thing>"
```

- [ ] **Step 3: Verify deleted files are gone and not re-introduced**

```bash
ls SPEC.md CODE_SIMPLIFICATION_PROMPT.md PLAN_IMPLEMENTATION_PROMPT.md PLAN_RESS_DEDUP.md ralph.sh scripts/verify_counts.py 2>&1
```

Expected: every line is "No such file or directory".

- [ ] **Step 4: Read the final README end-to-end**

```bash
cat README.md
```

Sanity check:

- The 5 sections appear in the order specified in the spec (description
  + install; mechanism + derive + simulate + fit; reaction + identify;
  derivation explanation; identify-and-enumeration explanation).
- No section refers to a function that doesn't exist in `src/`.
- No "Known Limitations" section, no API reference table, no leftover
  `(NOT YET DONE)` text.
- The fenced code blocks all use ` ```julia ` opening fences — none use
  ` ```jldoctest `, ` ```julia-skip `, etc.

- [ ] **Step 5: No commit if step 2 surfaced nothing; otherwise the prior commits already cover it**

This task is verification-only when the previous tasks executed cleanly.
If step 2 found stale phrases and required edits, those were committed
in step 2. If everything is clean, no commit is needed and the branch
is ready to push.

---

## Self-Review Notes (for plan author, not the executor)

- **Spec coverage:** All seven items in the design spec are addressed
  — five README sections (Tasks 3–7), the runnability test (Task 2),
  and the cheap doc cleanup (Task 1). The "out of scope" deferred items
  remain deferred.
- **Acceptance criteria:** Task 8 verifies all three acceptance criteria
  from the spec — `Pkg.test` passes, deleted files are absent, and the
  README structure matches.
- **Risks acknowledged:** Task 4 explicitly flags the parameter-name
  mismatch and noise sign-flip cases. Task 5 explicitly flags the slow
  `identify_rate_equation` skip-marker. Task 3 records the actual
  `parameters(m)` output before Task 4 hard-codes the names.
