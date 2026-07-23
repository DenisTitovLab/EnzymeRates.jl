# Enumeration combinatorics doc, shared_catalytic_site, raw-loss LOOCV — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `shared_catalytic_site` constraint to `EnzymeReaction`, remove the redundant outer `log` from LOOCV aggregation, and document how the enumeration count explodes.

**Architecture:** Three independent changes on branch `enum-combinatorics-shared-catalytic-site-loocv`. Task 1 removes the outer log from LOOCV model selection. Tasks 2–4 add the `shared_catalytic_site` field, DSL label, and enumeration filter. Task 5 documents the enumeration combinatorics. Each task ends with its own tests and commit.

**Tech Stack:** Julia, Documenter.jl, Test stdlib, DataFrames.

## Global Constraints

- 92-character line length, 4-space indentation.
- Names must describe what code does; after removing the log, any name containing "log" that no longer logs is false and must be renamed.
- TDD: write the failing test, run it red, implement, run it green, commit.
- `rate_equation` performance contract is untouched by every task here.
- Design spec: `docs/superpowers/specs/2026-07-21-enumeration-combinatorics-shared-catalytic-site-loocv-loss-design.md`.
- Run tests with `julia --project=. -e 'using Pkg; Pkg.test()'` (cold — pays precompilation + JIT each time). Individual test files are NOT standalone-includable: they rely on `test/runtests.jl`'s `using Test, EnzymeRates, …` preamble and shared fixtures. Always run the full suite for the red/green cycle and before every commit.

---

## Task 1: Raw-loss LOOCV (Change 3)

Remove the outer `log` from `cv_score`, the paired mean/SE, and the permutation-test input, so model selection compares raw per-fold losses. The per-fold loss is already a mean squared log-ratio; the outer log was a double log.

**Files:**
- Modify: `src/identify_rate_equation.jl` (lines 1080–1090, 1108–1117, 1124, 1189, 1206–1214, 1265–1267, and docstrings 152–160, 187–205, 1036–1070)
- Modify: `docs/src/identify/model_selection.md` (lines 75–90, 94, 99–101, 142–146)
- Test: `test/test_identify_rate_equation.jl` (lines 428–502, 671–880)

**Interfaces:**
- Produces: `_select_best_n_params` returns `diagnostics` whose per-bucket NamedTuple field is renamed `mean_log_loss_diff` → `mean_loss_diff` (fields now `(mean_loss_diff, se_paired, permutation_p)`). `cv_results` column renamed identically. `_cv_fold_loss` returns the raw finite fold loss (no eps floor).

- [ ] **Step 1: Recalibrate the fold-score fixtures and rename the diagnostic field in the tests (write the failing tests)**

In `test/test_identify_rate_equation.jl`, the `_select_best_n_params` fixtures author fold scores in log space via `exp.([...])` so the code's internal `log` recovers them. With the log removed, drop the `exp.(` wrapper (and its closing `)`) from every `cv_fold_scores` entry in the four testsets (`paired SE math`, `AND-combiner truth table`, `edge cases`, and the single/tie fixtures). The authored numbers now ARE the raw fold scores, so every expected diff/SE/mean value stays the same. Concretely:

- In each `cv_fold_scores = [ exp.([...]), exp.([...]) ]`, rewrite to `[ [...], [...] ]`.
- Rename every `mean_log_loss_diff` to `mean_loss_diff` (fixtures at lines 693, 696, 717, 735, 760, 798, 801 and any others — grep the file).
- Line 736–737 `mean(log.(exp.(...)) .- log.(exp.(...)))` becomes `mean([0.11, 0.12, 0.115, 0.115] .- [0.09, 0.13, 0.10, 0.12])`.
- Line 697–699 comment about "FP roundoff through exp/log makes std(diffs) ≈ 1e-17": with exp/log gone, `std([0.5,0.5,0.5,0.5,0.5,0.5]) == 0.0` exactly. Change the assertion `isapprox(d5.se_paired, 0.0; atol = 1e-10)` to `d5.se_paired == 0.0` and replace the stale comment with `# uniform diffs → std is exactly 0`.
- In the `_cv_fold_loss` testsets, remove the eps-floor assumption: line 456 `all(s -> s >= eps(Float64), scores)` → `all(s -> s >= 0.0, scores)`; line 501 `one >= eps(Float64) && isfinite(one)` → `one >= 0.0 && isfinite(one)`. Rename the two testset titles that say "floored at eps" / "floored + finite" (lines 428, 483) to "…: per-fold scores, finite" and "…: one fold, finite".

- [ ] **Step 2: Run the tests to verify they fail**

Run: `julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -40` (or run just `test/test_identify_rate_equation.jl` through the runner).
Expected: FAIL — the code still emits `mean_log_loss_diff` (KeyError/`type NamedTuple has no field mean_loss_diff`) and still applies `log`, so the recalibrated numeric expectations mismatch.

- [ ] **Step 3: Remove the outer log and rename in `src/identify_rate_equation.jl`**

In `_select_best_n_params` (line 1072+):
- Line 1080–1082: rename `log_scores` → `fold_scores` and drop the log:
  ```julia
      fold_scores = Dict(
          row.n_params => row.cv_fold_scores
          for row in eachrow(reps))
      means = Dict(n => mean(fs) for (n, fs) in fold_scores)
  ```
- Line 1083 `log_means` → `means` (done above); line 1086 `argmin(n -> (log_means[n], n), keys(log_means))` → `argmin(n -> (means[n], n), keys(means))`; line 1087 `length(log_scores[n_min])` → `length(fold_scores[n_min])`.
- Line 1084–1085 comment "Tie-break on log-mean" → "Tie-break on mean".
- Lines 1089–1096 NamedTuple field `mean_log_loss_diff` → `mean_loss_diff` (declaration and the `diagnostics[n_min]` literal).
- Lines 1099–1108: `keys(log_means)` → `keys(means)` (lines 1099, 1101); `ls = log_scores[n]` → `fs = fold_scores[n]` (line 1103) and its uses (1104 `length(ls)`, 1108 `diffs = ls .- log_scores[n_min]`) → `fs` / `fold_scores[n_min]`.
- Line 1114 `mean_log_loss_diff = md` → `mean_loss_diff = md`.
- Line 1124 `d.mean_log_loss_diff` → `d.mean_loss_diff`.

In `_cv_model_selection`:
- Line 1189: `cv_df.cv_score = [mean(log.(v)) for v in cv_df.cv_fold_scores]` → `cv_df.cv_score = [mean(v) for v in cv_df.cv_fold_scores]`.
- Line 1208: `for fld in (:mean_log_loss_diff, :se_paired, :permutation_p)` → `for fld in (:mean_loss_diff, :se_paired, :permutation_p)`.

In `_cv_fold_loss` (line 1244+):
- Delete the eps floor (lines 1265–1267). The function's last expression becomes `test_loss` (the value already asserted finite at 1262). Remove the two-line comment "Floor at eps so log(score) is finite…". Update the function docstring line 1240 "return the test loss floored at `eps(Float64)`" → "return the finite test loss".

Docstrings:
- Lines 152–160 (`se_threshold`, `perm_p_threshold` kwargs): "mean paired log-loss difference" → "mean paired loss difference".
- Lines 187–205 (Model selection section): "lowest mean log-fold-loss" → "lowest mean fold-loss"; "paired log-loss differences" → "paired fold-loss differences"; "log-fold-loss" → "fold-loss"; column name `mean_log_loss_diff` → `mean_loss_diff` (line 203).
- Lines 1036–1070 (`_select_best_n_params` docstring): "log-transformed per-fold LOOCV scores" → "per-fold LOOCV scores"; the `mean_log_loss_diff` field name (1045) → `mean_loss_diff`; "lowest mean log-fold-loss" (1044, 1051) → "lowest mean fold-loss"; "`n_min`'s log-folds" (1051) → "`n_min`'s folds"; "mean(`log_scores`)" (1060) → "mean(`fold_scores`)"; "on `log_scores[n_min]`" wording aligned to `fold_scores`.

Verify nothing else references the old names:
Run: `grep -rn "mean_log_loss_diff\|log_scores\|log_means\|log-fold-loss\|log per-fold\|log.(v)\|floored at eps\|log(score)" src/ | grep -v "log_abs_rates"`
Expected: no matches (the inner `loss!` log in `fitting.jl` and `log_abs_rates` are unrelated and must remain).

- [ ] **Step 4: Update `docs/src/identify/model_selection.md`**

- Lines 75–78: the fold-score description says it is "floored at `eps(Float64)` so the log stays finite…". Rewrite: "For each fold the mechanism is fit on all other groups and scored on the held-out group. Leaving one group out this way estimates how well the mechanism predicts a new experiment…" (drop the floor/log clause).
- Lines 82–90: replace the cv_score block:
  ```markdown
  The CV score for a mechanism is the **mean of the per-fold losses**:

  ```julia
  cv_score = mean(fold_scores)
  ```

  Each fold loss is already a mean squared log-ratio (see the loss definition),
  so averaging the fold losses directly is the natural aggregate — no further
  transform is applied.
  ```
- Line 94: "lowest mean log CV score" → "lowest mean CV score".
- Lines 99–101: "mean of the paired log-fold-loss differences" → "mean of the paired fold-loss differences".
- Line 142: table row `cv_score` "Mean of log per-fold losses (lower is better)." → "Mean of per-fold losses (lower is better)."
- Line 143: table row rename `mean_log_loss_diff` → `mean_loss_diff`, description "Mean paired fold-loss difference vs the `n_min` bucket's representative. `0.0` for `n_min`."

- [ ] **Step 5: Run the tests to verify they pass**

Run: `julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -40`
Expected: PASS — full suite green, test output pristine.

- [ ] **Step 6: Commit**

```bash
git add src/identify_rate_equation.jl docs/src/identify/model_selection.md test/test_identify_rate_equation.jl
git commit -m "LOOCV model selection: aggregate raw per-fold loss, not log-loss

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01B8LHSsU7VvCZ7Tr6ZzXF5N"
```

---

## Task 2: `shared_catalytic_site` field on `EnzymeReaction` (Change 2a)

Add the field, constructor keyword, validation, `==`/`hash`, and accessor.

**Files:**
- Modify: `src/types.jl` (struct 329–390, `==`/`hash` 401–408, accessors 392–399, struct docstring above 329)
- Test: `test/test_types.jl` (new testset near the existing `EnzymeReaction` construction tests, ~line 870)

**Interfaces:**
- Produces: field `shared_catalytic_site::Vector{Tuple{Symbol,Symbol}}`, keyword `shared_catalytic_site` (default `Tuple{Symbol,Symbol}[]`) on the `EnzymeReaction` inner constructor, and accessor `shared_catalytic_site(r::EnzymeReaction)::Vector{Tuple{Symbol,Symbol}}`. Stored normalized to `(substrate, product)` order, deduplicated, sorted. Consumed by Task 3 (DSL) and Task 4 (enforcement).

- [ ] **Step 1: Write the failing constructor tests**

Add to `test/test_types.jl` (follow the existing `ReactantAtoms`/`EnzymeReaction` construction pattern at lines 873–903):

```julia
@testset "EnzymeReaction shared_catalytic_site" begin
    A = EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:A), [:C => 1])
    B = EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:B), [:C => 1])
    P = EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])
    Q = EnzymeRates.ReactantAtoms(EnzymeRates.Product(:Q), [:C => 1])
    noregs = EnzymeRates.RegulatorMults[]

    # Accepts a valid (substrate, product) pair; stored normalized + sorted.
    r = EnzymeRates.EnzymeReaction([A, B, P, Q], noregs, Int[1];
        shared_catalytic_site = [(:A, :P)])
    @test EnzymeRates.shared_catalytic_site(r) == [(:A, :P)]

    # Normalizes product-first input to (substrate, product).
    r2 = EnzymeRates.EnzymeReaction([A, B, P, Q], noregs, Int[1];
        shared_catalytic_site = [(:P, :A)])
    @test EnzymeRates.shared_catalytic_site(r2) == [(:A, :P)]

    # Sorts and dedups multiple pairs.
    r3 = EnzymeRates.EnzymeReaction([A, B, P, Q], noregs, Int[1];
        shared_catalytic_site = [(:B, :Q), (:A, :P)])
    @test EnzymeRates.shared_catalytic_site(r3) == [(:A, :P), (:B, :Q)]

    # Default is empty.
    @test EnzymeRates.shared_catalytic_site(
        EnzymeRates.EnzymeReaction([A, P], noregs, Int[1])) ==
        Tuple{Symbol,Symbol}[]

    # Rejections.
    @test_throws ErrorException EnzymeRates.EnzymeReaction(  # unknown name
        [A, B, P, Q], noregs, Int[1]; shared_catalytic_site = [(:A, :Z)])
    @test_throws ErrorException EnzymeRates.EnzymeReaction(  # two substrates
        [A, B, P, Q], noregs, Int[1]; shared_catalytic_site = [(:A, :B)])
    @test_throws ErrorException EnzymeRates.EnzymeReaction(  # two products
        [A, B, P, Q], noregs, Int[1]; shared_catalytic_site = [(:P, :Q)])
    @test_throws ErrorException EnzymeRates.EnzymeReaction(  # duplicate pair
        [A, B, P, Q], noregs, Int[1]; shared_catalytic_site = [(:A, :P), (:P, :A)])

    # A metabolite that is BOTH a substrate and a competitive inhibitor may
    # still be the substrate side of a shared pair (validation keys on the
    # substrate/product role, not the inhibitor role).
    regA = EnzymeRates.RegulatorMults(EnzymeRates.CompetitiveInhibitor(:A), [1])
    r_dual = EnzymeRates.EnzymeReaction([A, B, P, Q], [regA], Int[1];
        shared_catalytic_site = [(:A, :P)])
    @test EnzymeRates.shared_catalytic_site(r_dual) == [(:A, :P)]

    # == / hash reflect the field.
    @test r != EnzymeRates.EnzymeReaction([A, B, P, Q], noregs, Int[1])
    @test hash(r) != hash(EnzymeRates.EnzymeReaction([A, B, P, Q], noregs, Int[1]))
    @test r == EnzymeRates.EnzymeReaction([A, B, P, Q], noregs, Int[1];
        shared_catalytic_site = [(:A, :P)])
end
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -30` (test files are not standalone-includable — they rely on `runtests.jl`'s `using` preamble, so run the full suite).
Expected: FAIL — `shared_catalytic_site` is not a keyword / accessor yet (`MethodError`/`UndefVarError`).

- [ ] **Step 3: Add the field, keyword, validation, accessor, `==`/`hash`**

In `src/types.jl`, extend the struct and inner constructor (329–390):

```julia
struct EnzymeReaction
    reactants::Vector{ReactantAtoms}
    regulators::Vector{RegulatorMults}
    allowed_catalytic_multiplicities::Vector{Int}
    shared_catalytic_site::Vector{Tuple{Symbol,Symbol}}

    function EnzymeReaction(reactants::Vector{ReactantAtoms},
                            regulators::Vector{RegulatorMults},
                            allowed_catalytic_multiplicities::Vector{Int};
                            shared_catalytic_site = Tuple{Symbol,Symbol}[])
```

After `sub_names`/`prod_names` are built (line 351–352), add validation that normalizes each pair to `(substrate, product)`:

```julia
        sub_set  = Set(sub_names)
        prod_set = Set(prod_names)
        normalized_shared = Tuple{Symbol,Symbol}[]
        for pair in shared_catalytic_site
            length(pair) == 2 || error(
                "EnzymeReaction: each shared_catalytic_site entry must name " *
                "exactly two metabolites, got $pair")
            a, b = pair[1], pair[2]
            if a in sub_set && b in prod_set
                push!(normalized_shared, (a, b))
            elseif b in sub_set && a in prod_set
                push!(normalized_shared, (b, a))
            else
                error("EnzymeReaction: shared_catalytic_site pair $pair must " *
                      "name one declared substrate and one declared product")
            end
        end
        sorted_shared = sort(unique(normalized_shared))
        length(sorted_shared) == length(normalized_shared) || error(
            "EnzymeReaction: duplicate shared_catalytic_site pair")
```

Change the final `new(...)` (line 388) to pass the field:

```julia
        new(sorted_reactants, sorted_regulators, sorted_mults, sorted_shared)
```

Add the accessor beside the others (after line 395):

```julia
shared_catalytic_site(r::EnzymeReaction) = r.shared_catalytic_site
```

Extend `==` (line 401–403) with `&& a.shared_catalytic_site == b.shared_catalytic_site`, and fold the field into `Base.hash` (line 404–407) as the outermost wrapper:

```julia
Base.hash(r::EnzymeReaction, h::UInt) =
    hash(r.shared_catalytic_site,
         hash(r.allowed_catalytic_multiplicities,
              hash(r.regulators,
                   hash(r.reactants, hash(:EnzymeReaction, h)))))
```

Update the `EnzymeReaction` struct docstring (above line 329) to document the `shared_catalytic_site` keyword: one sentence — a list of `(substrate, product)` pairs that bind the same catalytic site and so are never both bound to it at once.

- [ ] **Step 4: Run to verify pass**

Run the full suite: `julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: PASS — the new testset green and no regression (the `==`/`hash` change and new field must not break existing reaction tests, which construct with three positional args).

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Add shared_catalytic_site field + validation to EnzymeReaction

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01B8LHSsU7VvCZ7Tr6ZzXF5N"
```

---

## Task 3: `shared_catalytic_site:` DSL label (Change 2b)

Parse the label in `@enzyme_reaction` and pass the pairs to the constructor.

**Files:**
- Modify: `src/dsl.jl` (`_VALID_REACTION_LABELS` 41–46, `_parse_reaction_block` 59–101, macro emit 31–39, new helpers near 100)
- Test: `test/test_dsl.jl` (new testset)

**Interfaces:**
- Consumes: the `EnzymeReaction(...; shared_catalytic_site=...)` keyword from Task 2.
- Produces: `@enzyme_reaction` accepts `shared_catalytic_site: (sub, prod), ...`, emitting a reaction whose `shared_catalytic_site(r)` equals the constructor-keyword result.

- [ ] **Step 1: Write the failing DSL tests**

Add to `test/test_dsl.jl`:

```julia
@testset "@enzyme_reaction shared_catalytic_site" begin
    rxn = @enzyme_reaction begin
        substrates: A[C], B[C]
        products:   P[C], Q[C]
        shared_catalytic_site: (A, P), (B, Q)
    end
    @test EnzymeRates.shared_catalytic_site(rxn) == [(:A, :P), (:B, :Q)]

    # Single pair.
    rxn1 = @enzyme_reaction begin
        substrates: A[C], B[C]
        products:   P[C], Q[C]
        shared_catalytic_site: (P, A)
    end
    @test EnzymeRates.shared_catalytic_site(rxn1) == [(:A, :P)]

    # Malformed: bare (unparenthesized) names error.
    @test_throws LoadError @eval @enzyme_reaction begin
        substrates: A[C]
        products:   P[C]
        shared_catalytic_site: A, P
    end

    # Malformed: three-name tuple errors.
    @test_throws LoadError @eval @enzyme_reaction begin
        substrates: A[C], B[C]
        products:   P[C], Q[C]
        shared_catalytic_site: (A, P, Q)
    end
end
```

(Macro-expansion errors surface as `LoadError` from `@eval`; if the local convention wraps differently, match the existing malformed-DSL tests in `test/test_dsl.jl`.)

- [ ] **Step 2: Run to verify failure**

Run: `julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: FAIL — `unknown label shared_catalytic_site:`.

- [ ] **Step 3: Implement the DSL label**

In `src/dsl.jl`:

Add `:shared_catalytic_site` to `_VALID_REACTION_LABELS` (line 41–46).

In `_parse_reaction_block` (59–101): add an accumulator `shared = Tuple{Symbol,Symbol}[]` beside `subs/prods/regs` (line 62–65), a branch in the `for` loop:

```julia
        elseif label === :shared_catalytic_site
            append!(shared, _parse_shared_site_pairs(values))
```

and return `shared` in the NamedTuple: `(; subs, prods, regs, mults, shared)`.

Add the helper (near `_parse_atom_bracket_entries`):

```julia
"""
Parse `shared_catalytic_site:` entries. Each value is a two-name tuple
`(sub, prod)`; roles are resolved by the `EnzymeReaction` constructor.
Returns `Vector{Tuple{Symbol,Symbol}}`.
"""
function _parse_shared_site_pairs(values)
    out = Tuple{Symbol,Symbol}[]
    for v in values
        v isa Expr && v.head === :tuple && length(v.args) == 2 &&
            all(a -> a isa Symbol, v.args) || error(
                "@enzyme_reaction `shared_catalytic_site:`: each entry must be " *
                "a `(substrate, product)` pair of two names; got $v.")
        push!(out, (v.args[1]::Symbol, v.args[2]::Symbol))
    end
    out
end

"""Build the `Vector{Tuple{Symbol,Symbol}}` `Expr` for shared-site pairs."""
function _build_shared_site_expr(shared)
    entries = [:(($(QuoteNode(a)), $(QuoteNode(b)))) for (a, b) in shared]
    :(Tuple{Symbol,Symbol}[$(entries...)])
end
```

In the macro (31–39), build and pass the keyword:

```julia
macro enzyme_reaction(block)
    parsed = _parse_reaction_block(block)
    reactants_expr  = _build_reactants_expr(parsed.subs, parsed.prods)
    mults_expr      = _build_catalytic_mults_expr(parsed.mults)
    regulators_expr = _build_regulators_expr(parsed.regs, parsed.mults)
    shared_expr     = _build_shared_site_expr(parsed.shared)
    return esc(:(EnzymeRates.EnzymeReaction(
        $reactants_expr, $regulators_expr, $mults_expr;
        shared_catalytic_site = $shared_expr,
    )))
end
```

Add a bullet to the `@enzyme_reaction` docstring (4–30) documenting `shared_catalytic_site:` with the `(sub, prod)` syntax.

- [ ] **Step 4: Run to verify pass**

Run: `julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/dsl.jl test/test_dsl.jl
git commit -m "Add shared_catalytic_site: label to @enzyme_reaction DSL

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01B8LHSsU7VvCZ7Tr6ZzXF5N"
```

---

## Task 4: Enforce `shared_catalytic_site` in enumeration (Change 2c)

Filter the substrate/product competition patterns so no dead-end form ever binds both members of a declared pair.

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_expand_substrate_product_dead_ends`, after line 975)
- Test: `test/test_mechanism_enumeration.jl` (new testset)

**Interfaces:**
- Consumes: `shared_catalytic_site(reaction)` from Task 2.
- Produces: `init_mechanisms(reaction)` yields no `Mechanism` whose catalytic-site species binds both members of a declared pair; the count strictly drops versus the unconstrained reaction.

- [ ] **Step 1: Write the failing enforcement tests**

Add to `test/test_mechanism_enumeration.jl`:

```julia
@testset "shared_catalytic_site prunes catalytic-site co-occupancy" begin
    rxn = @enzyme_reaction begin
        substrates: A[C], B[C]
        products:   P[C], Q[C]
        shared_catalytic_site: (A, P)
    end
    _reactant_names(sp) = Set(EnzymeRates.name(b)
        for b in EnzymeRates.bound(sp) if b isa EnzymeRates.Reactant)
    binds_both(sp) = :A in _reactant_names(sp) && :P in _reactant_names(sp)
    @test !any(binds_both(sp)
        for m in EnzymeRates.init_mechanisms(rxn)
        for g in EnzymeRates.steps(m) for s in g
        for sp in (EnzymeRates.from_species(s), EnzymeRates.to_species(s)))
end

@testset "shared_catalytic_site strictly reduces mechanism count" begin
    base = @enzyme_reaction begin
        substrates: A[C], B[C]
        products:   P[C], Q[C]
    end
    constrained = @enzyme_reaction begin
        substrates: A[C], B[C]
        products:   P[C], Q[C]
        shared_catalytic_site: (A, P)
    end
    n_base = length(unique!(collect(EnzymeRates.init_mechanisms(base))))
    n_con  = length(unique!(collect(EnzymeRates.init_mechanisms(constrained))))
    @test n_con < n_base
end
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: FAIL — without the filter, `init_mechanisms` builds an `E·A·P` dead-end form, so `binds_both` is true and `n_con == n_base`.

- [ ] **Step 3: Add the pattern filter**

In `_expand_substrate_product_dead_ends` (`src/mechanism_enumeration.jl`), directly after `patterns = _competition_patterns(sub_names, prod_names)` (line 974–975):

```julia
    # A declared shared catalytic site forbids its (substrate, product) pair
    # from co-occupying the catalytic site, so keep only competition patterns
    # whose forbidden-edge set contains every declared pair. The complete
    # bipartite pattern contains all edges and always survives, so the list is
    # never empty.
    shared = shared_catalytic_site(reaction)
    if !isempty(shared)
        patterns = filter(
            pat -> all(edge -> edge in pat, shared), patterns)
    end
```

- [ ] **Step 4: Run to verify pass**

Run: `julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: PASS — full suite green.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Enforce shared_catalytic_site by pruning competition patterns

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01B8LHSsU7VvCZ7Tr6ZzXF5N"
```

---

## Task 5: Enumeration combinatorics doc subsection (Change 1)

Append a subsection to the enumeration page deriving the candidate count as a product of combinatorial factors, plugged into a random bi-bi with two regulators, closing with the filters (including `shared_catalytic_site`).

**Files:**
- Modify: `docs/src/identify/enumeration_engine.md` (append after line 181)
- Scratch (not committed): a Julia snippet under the session scratchpad to verify factor formulas.

**Interfaces:** none (prose only; no runnable block ships).

- [ ] **Step 1: Verify the factor formulas against real counts (offline, not committed)**

Write a scratch script and run it to anchor the numbers:

```julia
using EnzymeRates, Combinatorics  # Combinatorics only if convenient; else hand-count
# Substrate/product competition-pattern count = bipartite graphs on
# n_s × n_p with every row/column degree ≥ 1 (inclusion–exclusion):
N(ns, np) = sum((-1)^(i+j) * binomial(ns, i) * binomial(np, j) *
                2.0^((ns - i) * (np - j)) for i in 0:ns, j in 0:np)
@show N(1, 1)          # expect 1
@show N(2, 1)          # expect 3
@show N(2, 2)          # expect 7  ← the bi-bi dead-end factor

# Cross-check the engine agrees for the dead-end/topology factors on small cases:
for (subs, prods) in [(["S"],["P"]), (["A","B"],["P","Q"])]
    rxn = # build via @enzyme_reaction with 1-C atoms per reactant, balanced
    @show length(unique!(collect(EnzymeRates.init_mechanisms(rxn))))
end
```

Record the observed `init_mechanisms` counts and confirm the derived per-stage factors reproduce them. These numbers back the prose; the script itself is not committed.

- [ ] **Step 2: Write the subsection**

Append to `docs/src/identify/enumeration_engine.md`:

```markdown
## How many mechanisms? The combinatorics of enumeration

The candidate count is a **product** of independent combinatorial choices, so it
grows explosively with the number of substrates, products, and regulators. Each
factor below is a function of the reaction's counts — `n_s` substrates, `n_p`
products, `n_i` competitive inhibitors, `n_a` allosteric regulators.

- **Catalytic binding topologies.** Ordered mechanisms fix a binding sequence;
  the random-order topology adds one more. The number of distinct sequences
  grows with `n_s!` and `n_p!`.
- **Substrate/product dead-end forms.** A dead-end form binds a mix of
  substrates and products at the catalytic site. Which pairs may co-occupy is a
  bipartite competition graph on `n_s × n_p` in which every reactant has degree
  at least one. The count is
  `N(n_s, n_p) = Σᵢⱼ (−1)^{i+j} C(n_s,i) C(n_p,j) 2^{(n_s−i)(n_p−j)}` —
  `N(1,1)=1`, `N(2,1)=3`, `N(2,2)=7`. This is exactly the choice
  `shared_catalytic_site` prunes.
- **The seven expansion moves**, each multiplying the count: rapid-
  equilibrium/steady-state flips (up to `2^{groups}`), kinetic-group splits,
  competitive-inhibitor sites (a second competition count over `n_i` inhibitors),
  promotion to allosteric, allosteric ligands and their `:OnlyA`/`:EqualAI`/
  `:NonequalAI` states, and regulatory-site merges.

**A worked example — random-order bi-bi with two regulators** (`n_s = n_p = 2`,
one competitive inhibitor, one allosteric regulator). Multiplying the binding
topologies, the `N(2,2)=7` dead-end forms, and the per-move factors lands in the
**tens of thousands** of candidate mechanisms. Structural deduplication and
validity pruning trim the raw product, but the growth *rate* is the point: one
more reactant or regulator multiplies the universe again.

Fitting even a fraction of that is infeasible, which is why the search is
filtered rather than exhaustive. The beam keeps only the promising candidates at
each parameter count; `eq_complexity_filter` drops equations whose
denominator is too dense to be identifiable; required-regulator seeding skips
the partially-regulated lower shelf; and `shared_catalytic_site` removes
mechanisms a chemist already knows are wrong. See
[Best mechanism selection](@ref) for the filters that bound the search.
```

Replace "tens of thousands" with the actual order of magnitude from Step 1 if it differs materially; keep the framing honest (a calculated estimate, not a measured enumeration count).

- [ ] **Step 3: Verify cross-references resolve**

Confirm each `@ref`/target used exists as a heading or docstring anchor:
Run: `grep -rn "Best mechanism selection\|eq_complexity_filter\|shared_catalytic_site" docs/src/`
Expected: `## The cross-validation selection rule` page titled "Best mechanism selection" exists in `make.jl`; `eq_complexity_filter` is a documented kwarg. If a target is not a valid Documenter anchor, use plain prose or the correct `@ref`. Best-effort: if the docs environment is instantiated, run `julia --project=docs docs/make.jl` and confirm no cross-reference warnings for the new section; otherwise rely on the grep check.

- [ ] **Step 4: Commit**

```bash
git add docs/src/identify/enumeration_engine.md
git commit -m "Doc: enumeration combinatorics subsection motivating filters

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01B8LHSsU7VvCZ7Tr6ZzXF5N"
```

---

## Final verification

- [ ] Run the full suite once more: `julia --project=. -e 'using Pkg; Pkg.test()'` — all green, output pristine.
- [ ] `grep -rn "competing_reactants\|mean_log_loss_diff" src/ test/ docs/src/` returns nothing.
- [ ] Review the branch as a whole (superpowers:requesting-code-review) before finishing (superpowers:finishing-a-development-branch).
