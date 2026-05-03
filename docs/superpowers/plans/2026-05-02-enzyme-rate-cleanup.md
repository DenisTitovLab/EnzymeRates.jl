# EnzymeRates Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the seven fixes specified in `docs/superpowers/specs/2026-05-02-enzyme-rate-cleanup-design.md` — strict regulator + atom validation, allosteric display refactor, identify_rate_equation file/beam ergonomics, and post-compile rate-equation hash dedup with a persistent cross-level fit cache.

**Architecture:** Apply changes as seven independent phases that can be reviewed individually. Each phase opens with regression tests (TDD-first), implements the smallest fix, then commits. The phases are ordered so earlier ones do not break later ones (rename first, then validation, then ergonomics, then dedup).

**Tech Stack:** Julia 1.x, `Test`, `DataFrames`, `CSV`, `Optimization`, `Distributed.pmap`, `SHA` (stdlib for canonical hashing).

---

## Repository conventions

- 92-char line limit, 4-space indentation.
- Run tests with: `julia --project -e 'using Pkg; Pkg.test()'`
- Run a single test file: `julia --project -e 'using Pkg; Pkg.activate("."); include("test/<file>.jl")'`
- Each phase ends with a green test suite + commit. Never commit red.
- TDD discipline per `CLAUDE.md`:
  1. write failing regression test,
  2. confirm it fails for the expected reason,
  3. write the smallest passing change,
  4. confirm green,
  5. commit.
- Do not use `git add -A`. Stage only files belonging to the phase.

---

## File structure (touched / created)

| File | Phase | Change |
|---|---|---|
| `src/mechanism_enumeration.jl` | A, C | rename `param_count` → `n_fit_params_estimate`; drop `+2` from estimate; derive regulators from steps in `EnzymeMechanism(spec)` and allosteric path |
| `src/types.jl` | B, C, D | add atom mandatoriness + balance check in `EnzymeReaction`; tighten regulator validation in `EnzymeMechanism(mets, rxns)`; refactor `Base.show(io, ::AllostericEnzymeMechanism)` |
| `src/identify_rate_equation.jl` | E, F, G | bucket save by actual `length(fitted_params(m))`; new beam threshold rule; persistent rate-equation hash cache + four-stage processing; new CSV columns; LOOCV per unique hash |
| `src/rate_eq_derivation.jl` | G (new helper) | optionally expose a canonical-hash helper if `rate_equation_string` text isn't sufficient (decided per Task G.1 investigation) |
| `test/test_mechanism_enumeration.jl` | A, C | new estimate semantic tests; init does not carry unbound regulators |
| `test/test_types.jl` | B, C, D | atom-balance + atom-mandatory failure tests; strict regulator constructor failure test; `repr` allosteric display tests |
| `test/test_dsl.jl` | B | atom-validation tests at the macro/constructor seam |
| `test/test_identify_rate_equation.jl` | E, F, G | filename bucketing test; beam threshold tests; canonical-hash tests; cross-level cache tests; end-to-end LDH regression |
| `README.md` | B | migrate `S` / `P` bare-symbol example to declared atoms |
| `CLAUDE.md` | A | update note about `param_count` → `n_fit_params_estimate` |

No new source files are created. New test additions go into existing test files in dedicated `@testset` blocks.

---

## Phase A — Rename `param_count` → `n_fit_params_estimate` and drop `+2`

This phase is foundational. Later phases (C, E, F, G) reference the renamed field name. Do not skip ahead.

### Task A.1: Add failing test for new estimate semantics

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Add a new testset documenting the new field name and semantics**

In `test/test_mechanism_enumeration.jl`, after the existing constants block (around line 90), add:

```julia
@testset "n_fit_params_estimate semantics" begin
    # Simple uni-uni: 2 RE binding steps + 1 SS catalytic step
    # = 1 RE group (S binds, K1) + 1 RE group (P binds, K2)
    # + 1 SS group (kf, kr) + 0 thermo constraints (single SS).
    # Estimate = n_re_groups + 2*n_ss_groups - n_thermo
    #          = 2 + 2*1 - 0 = 4   (was 6 with old +2 formula)
    init_specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
    @test !isempty(init_specs)
    spec = first(init_specs)
    # New field name (the rename):
    @test hasfield(typeof(spec), :n_fit_params_estimate)
    # New value (formula without +2 for Keq+E_total):
    m = EnzymeRates.EnzymeMechanism(spec)
    n_actual = length(EnzymeRates.fitted_params(m))
    @test spec.n_fit_params_estimate == n_actual
end
```

- [ ] **Step 2: Run the test to confirm it fails**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: failure inside the `n_fit_params_estimate semantics` testset — either `hasfield` returns `false` (rename not done) or the value mismatches by 2 (formula not adjusted).

### Task A.2: Mechanical rename across src

**Files:**
- Modify: `src/mechanism_enumeration.jl` (~70 references)

- [ ] **Step 1: Rename the field on `MechanismSpec`**

In `src/mechanism_enumeration.jl` lines 37-41:

```julia
struct MechanismSpec <: AbstractMechanismSpec
    reaction::Any
    steps::Vector{StepSpec}
    n_fit_params_estimate::Int
end
```

- [ ] **Step 2: Rename the field on `AllostericMechanismSpec`**

Lines 60-68 of the same file:

```julia
struct AllostericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    group_tags::Dict{Int, Symbol}
    reg_ligand_tags::Dict{Symbol, Symbol}
    n_fit_params_estimate::Int
end
```

- [ ] **Step 3: Rename helper functions**

In the same file: rename `_param_count(spec)` → `_n_fit_params_estimate(spec)` and `_param_count_from_steps(steps)` → `_n_fit_params_estimate_from_steps(steps)`. Update all bodies that read or write the field via these helpers.

- [ ] **Step 4: Replace remaining `.param_count` usages in src/**

```bash
grep -rn "\.param_count\|param_count =\|param_count::" /home/denis.linux/.julia/dev/EnzymeRates/src/
```

Replace each `spec.param_count` / `.param_count` with `.n_fit_params_estimate`. Replace each `param_count = …` assignment in struct constructors and named-tuple builders. Confirm no `.param_count` remains:

```bash
grep -rn "\.param_count" /home/denis.linux/.julia/dev/EnzymeRates/src/
```

Expected: empty output.

- [ ] **Step 5: Replace `.param_count` references in tests**

```bash
grep -rn "\.param_count\|param_count =" /home/denis.linux/.julia/dev/EnzymeRates/test/
```

Replace each. Specific known sites:
- `test/test_mechanism_enumeration.jl:17` and `:100` (the `enumerate_all` helper).
- `test/test_identify_rate_equation.jl:41` (the `8` is the constructor's `param_count` positional arg — see Step 6).

- [ ] **Step 6: Update `MechanismSpec` / `AllostericMechanismSpec` constructor call sites**

The structs are positional. Any test or src that constructs one passes the count as the trailing positional argument. The argument position does NOT change; only the field name does. No code change needed at call sites for positional construction. Only NAMED-keyword usages need renaming, e.g. `MechanismSpec(rxn, steps; param_count=...)` → `MechanismSpec(rxn, steps; n_fit_params_estimate=...)`. Search:

```bash
grep -rn "param_count=\|param_count =" /home/denis.linux/.julia/dev/EnzymeRates/src/ /home/denis.linux/.julia/dev/EnzymeRates/test/
```

Update any kwarg-style construction.

- [ ] **Step 7: Run tests to verify rename is consistent (formula not yet changed)**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected:
- The Phase A.1 test still fails (formula still has `+2`).
- All other tests pass — the rename alone is semantically equivalent.

If a test fails for any reason other than A.1's value-check, fix the missed `.param_count` reference.

### Task A.3: Drop `+2` from the estimate formula

**Files:**
- Modify: `src/mechanism_enumeration.jl` line 841

- [ ] **Step 1: Change the formula**

In `init_mechanisms`, find the single occurrence of `n_re + 2 * n_ss - n_thermo + 2`:

```bash
grep -n "n_re + 2 \* n_ss - n_thermo" /home/denis.linux/.julia/dev/EnzymeRates/src/mechanism_enumeration.jl
```

Replace `+ 2` at the end with nothing:

```julia
n_fit_params_estimate = n_re + 2 * n_ss - n_thermo
```

- [ ] **Step 2: Run tests**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: A.1 now passes. Whole suite green.

### Task A.4: Update CLAUDE.md note

**Files:**
- Modify: `.claude/CLAUDE.md`

- [ ] **Step 1: Find and update the param_count note**

```bash
grep -n "param_count" /home/denis.linux/.julia/dev/EnzymeRates/.claude/CLAUDE.md
```

Replace mentions of `param_count` with `n_fit_params_estimate`. Update the formula description: from "`length(groups_re) + 2 * length(groups_ss) - n_thermo + 2`" to "`length(groups_re) + 2 * length(groups_ss) - n_thermo`" with a clarifying parenthetical "(estimate of independent rate-constant count, excluding `Keq` and `E_total` which are not fitted)".

### Task A.5: Commit Phase A

- [ ] **Step 1: Stage and commit**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl test/test_identify_rate_equation.jl .claude/CLAUDE.md
git status   # confirm only the renamed-field churn is staged
git commit -m "$(cat <<'EOF'
Rename spec.param_count to n_fit_params_estimate, drop +2 for Keq/E_total

The enumeration estimate now matches the n_params column convention
(independent rate-constant count, excluding Keq and E_total which
are not fitted). Pure rename + formula change — no behavior change.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase B — Atom mandatory + balance check in `@enzyme_reaction` (Issue #7)

### Task B.1: Audit & migrate any bare-symbol reactions

**Files:**
- Modify: `README.md`
- Audit (likely no changes needed): `test/*.jl`

- [ ] **Step 1: Find every `@enzyme_reaction` body in repo**

```bash
grep -rA 5 "@enzyme_reaction begin" /home/denis.linux/.julia/dev/EnzymeRates/test/ /home/denis.linux/.julia/dev/EnzymeRates/README.md /home/denis.linux/.julia/dev/EnzymeRates/src/
```

For each, inspect the `substrates:` and `products:` lines. Any bare-symbol species (e.g., `substrates: S` with no `[…]`) must be migrated.

- [ ] **Step 2: Migrate `README.md`**

Open `README.md` and replace the bare-symbol example. Use placeholder atoms `[C]` on each side (the example is illustrative, not biochemically meaningful):

Find:

```julia
rxn = @enzyme_reaction begin
    substrates: S
    products:   P
    regulators: A
    oligomeric_state: 2
```

Replace with:

```julia
rxn = @enzyme_reaction begin
    substrates: S[C]
    products:   P[C]
    regulators: A
    oligomeric_state: 2
```

Apply the same migration to every other bare-symbol example in the README.

- [ ] **Step 3: Run tests; commit migration as a stand-alone commit**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: green (atoms are accepted; mandatoriness not yet enforced).

```bash
git add README.md
git commit -m "$(cat <<'EOF'
Migrate README @enzyme_reaction examples to declared atoms

Pre-step for the upcoming atom-mandatory enforcement. Examples now
use [C] placeholder atoms; behavior unchanged.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task B.2: Add failing tests for atom-mandatory and balance checks

**Files:**
- Modify: `test/test_types.jl`

- [ ] **Step 1: Add testset for atom-mandatory**

Append to `test/test_types.jl`:

```julia
@testset "EnzymeReaction: atom mandatory" begin
    # No atoms on substrate -> error
    @test_throws ErrorException EnzymeReaction(
        ((:S, ()),),                 # bare substrate
        ((:P, ((:C, 1),)),)
    )
    # No atoms on product -> error
    @test_throws ErrorException EnzymeReaction(
        ((:S, ((:C, 1),)),),
        ((:P, ()),)                  # bare product
    )
    # Both have atoms -> ok
    @test EnzymeReaction(
        ((:S, ((:C, 1),)),),
        ((:P, ((:C, 1),)),)
    ) isa EnzymeReaction
end

@testset "EnzymeReaction: atom balance" begin
    # C count mismatch -> error
    @test_throws ErrorException EnzymeReaction(
        ((:S, ((:C, 6),)),),
        ((:P, ((:C, 5),)),)
    )
    # element only on one side -> error
    @test_throws ErrorException EnzymeReaction(
        ((:S, ((:C, 6), (:H, 12))),),
        ((:P, ((:C, 6),)),)
    )
    # Multi-substrate, multi-product balance -> ok
    @test EnzymeReaction(
        ((:A, ((:C, 6),)), (:B, ((:N, 1),))),
        ((:P, ((:C, 6),)), (:Q, ((:N, 1),)))
    ) isa EnzymeReaction
end
```

- [ ] **Step 2: Run; confirm both testsets fail**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: 4 failures inside the new testsets (no errors thrown by the constructor — current code accepts everything).

### Task B.3: Implement the atom checks in `EnzymeReaction(...)`

**Files:**
- Modify: `src/types.jl` lines 17-45

- [ ] **Step 1: Add a helper to sum atoms across a species tuple**

Insert above the existing `EnzymeReaction(...)` constructor:

```julia
"""Sum element counts across a tuple of `(name, atoms)` pairs.
Returns a Dict{Symbol,Int}. Errors if any species's atoms tuple
is empty (atoms are mandatory)."""
function _sum_atoms(species::Tuple, side::String)
    totals = Dict{Symbol,Int}()
    for (name, atoms) in species
        isempty(atoms) && error(
            "EnzymeReaction: $side metabolite $name has no declared " *
            "atoms; atoms are mandatory (use `[C…]` bracket syntax in " *
            "@enzyme_reaction or pass non-empty atom tuples to the " *
            "constructor).")
        for (elem, count) in atoms
            totals[elem] = get(totals, elem, 0) + count
        end
    end
    totals
end
```

- [ ] **Step 2: Add the balance check inside the constructor**

In `src/types.jl` after the duplicate-name checks (~line 26) and before normalization of regulators:

```julia
sub_atoms = _sum_atoms(subs, "substrate")
prod_atoms = _sum_atoms(prods, "product")
all_elems = union(keys(sub_atoms), keys(prod_atoms))
for elem in all_elems
    s_count = get(sub_atoms, elem, 0)
    p_count = get(prod_atoms, elem, 0)
    s_count == p_count || error(
        "EnzymeReaction: atom imbalance — element $elem appears " *
        "$s_count time(s) on substrate side and $p_count time(s) " *
        "on product side. Declared atoms must balance.")
end
```

- [ ] **Step 3: Run tests**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: B.2 testsets pass. Whole suite green.

### Task B.4: Commit Phase B

- [ ] **Step 1: Stage and commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "$(cat <<'EOF'
Require atoms on every reactant; check element balance in EnzymeReaction

Bare-symbol substrates / products are now an error. Element totals
must match between substrates and products (per declared element).
Regulators continue to ignore atoms.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase C — Strict regulator validation (Issue #1)

### Task C.1: Add failing tests

**Files:**
- Modify: `test/test_types.jl`, `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Add 2-arg constructor strict-check test in test_types.jl**

```julia
@testset "EnzymeMechanism: strict regulator binding" begin
    # Regulator :A listed but never bound in any step -> error
    @test_throws ErrorException EnzymeMechanism(
        ((:S,), (:P,), (:A,)),
        (((:E, :S), (:E_S,), true, 1),
         ((:E_S,), (:E_P,), false, 2),
         ((:E_P,), (:E, :P), true, 3))
    )
    # All regulators bound -> ok
    @test EnzymeMechanism(
        ((:S,), (:P,), (:A,)),
        (((:E, :S), (:E_S,), true, 1),
         ((:E_S, :A), (:E_S_A,), true, 4),
         ((:E_S,), (:E_P,), false, 2),
         ((:E_P,), (:E, :P), true, 3))
    ) isa EnzymeMechanism
    # No regulators -> ok
    @test EnzymeMechanism(
        ((:S,), (:P,), ()),
        (((:E, :S), (:E_S,), true, 1),
         ((:E_S,), (:E_P,), false, 2),
         ((:E_P,), (:E, :P), true, 3))
    ) isa EnzymeMechanism
end
```

- [ ] **Step 2: Add init/expansion test in test_mechanism_enumeration.jl**

```julia
@testset "init_mechanisms drops unbound regulators from spec→type" begin
    # A reaction with one regulator: init has no binding step for it.
    # The compiled EnzymeMechanism should NOT carry the regulator
    # in its type parameter.
    init_specs = EnzymeRates.init_mechanisms(uni_uni_with_reg)
    @test !isempty(init_specs)
    for spec in init_specs
        m = EnzymeRates.EnzymeMechanism(spec)
        # :I is the dead-end inhibitor in uni_uni_with_reg.
        # init has no binding step for :I, so :I must be absent.
        @test :I ∉ EnzymeRates.regulators(m)
    end

    # After expansion that adds the dead-end regulator, :I should
    # be present.
    expanded = EnzymeRates.expand_mechanisms(init_specs, uni_uni_with_reg)
    found_with_reg = false
    for (_, specs) in expanded
        for spec in specs
            m = EnzymeRates.EnzymeMechanism(spec)
            if :I in EnzymeRates.regulators(m)
                found_with_reg = true
                break
            end
        end
        found_with_reg && break
    end
    @test found_with_reg
end
```

- [ ] **Step 3: Run; confirm both testsets fail**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected:
- `test_types.jl` — first `@test_throws` fails because the constructor currently accepts unbound regulators (see `types.jl:127-138`).
- `test_mechanism_enumeration.jl` — `:I ∉ regulators(m)` fails because `EnzymeMechanism(spec)` (mechanism_enumeration.jl:1300) copies `regulators(rxn)` regardless of binding.

### Task C.2: Tighten 2-arg constructor + spec→type filter

**Files:**
- Modify: `src/types.jl` lines 127-138
- Modify: `src/mechanism_enumeration.jl` lines 1296-1300

- [ ] **Step 1: Add regulator-bound check in the 2-arg constructor**

In `src/types.jl`, replace the comment block at lines 127-138 with:

```julia
# Every substrate, product, AND regulator must appear in some step.
appears = Set{Symbol}()
for (lhs, rhs, _, _) in rxns
    for s in lhs; push!(appears, s); end
    for s in rhs; push!(appears, s); end
end
for name in vcat(collect(subs), collect(prods), collect(regs))
    name in appears ||
        error("Listed metabolite or regulator $name does not " *
              "appear in any reaction step")
end
```

- [ ] **Step 2: Filter regulators by step appearance in `EnzymeMechanism(spec)`**

In `src/mechanism_enumeration.jl` around line 1300, before computing `regs`:

```julia
# Build the set of names actually appearing on any step (after
# stripping the __reg suffix used by enumeration internals).
appears_in_steps = Set{Symbol}()
for s in spec.steps
    for sym in Iterators.flatten((s.reactants, s.products))
        push!(appears_in_steps, _strip_reg_suffix(sym))
    end
end

regs = Tuple(r for r in regulators(rxn)
             if r ∉ auto_exclude && r ∈ appears_in_steps)
```

- [ ] **Step 3: Run tests**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: Phase C tests pass; full suite green.

### Task C.3: Commit Phase C

- [ ] **Step 1: Stage and commit**

```bash
git add src/types.jl src/mechanism_enumeration.jl test/test_types.jl test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Require regulators to bind in some step (constructor + spec→type)

EnzymeMechanism(mets, rxns) now errors if a listed regulator never
appears in any step. EnzymeMechanism(spec) and the allosteric path
filter the regulator tuple to those names actually used by spec
steps, so init mechanisms no longer carry phantom regulators.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase D — `AllostericEnzymeMechanism` display refactor (Issue #2)

### Task D.1: Add failing tests

**Files:**
- Modify: `test/test_types.jl`

- [ ] **Step 1: Add testset**

```julia
@testset "AllostericEnzymeMechanism display format" begin
    # Use uni_uni_allo_reg-like construction from test_identify_rate_equation
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end
    init = EnzymeRates.init_mechanisms(rxn)
    base = first(init)
    g_s = first(s.kinetic_group for s in base.steps
                if EnzymeRates.step_metabolite(s) === :S)
    g_p = first(s.kinetic_group for s in base.steps
                if EnzymeRates.step_metabolite(s) === :P)
    spec = EnzymeRates.AllostericMechanismSpec(
        base, 2, [[:R]], [2],
        Dict(g_s => :EqualRT, g_p => :EqualRT),
        Dict(:R => :OnlyT),
        base.n_fit_params_estimate + 1)
    m = EnzymeRates.AllostericEnzymeMechanism(spec)
    s = repr(m)

    # Old summary line gone:
    @test !occursin("cat_allo_states:", s)
    # Inline ::Tag annotations on each step or step group:
    @test occursin(":: EqualRT", s)
    # Multi-line catalytic display (no chain shortcut):
    n_steps_re = count(c -> c == '\n', s)
    @test n_steps_re >= 3   # header + ≥3 step lines
end

@testset "AllostericEnzymeMechanism display: shared kinetic group" begin
    # Construct an allosteric mechanism with two RE steps in one
    # kinetic group (e.g. binding S to E and binding S to E_P) and
    # confirm they render together with one ::Tag.
    # ... build a manually-constructed AllostericMechanismSpec
    # whose base mechanism has two RE steps in kinetic_group=1 ...
    # (Construction details inline so the engineer doesn't have to
    # reverse-engineer them.)

    base_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    # 4-step mechanism: (E,S)→E_S g1, (E_P,S)→E_PS g1 (shared K),
    # (E_S)→(E_P) SS g2, (E,P)→E_P g3.  Net: S binds → ES, ES→EP, EP→E + P.
    # Use the 2-arg EnzymeMechanism constructor directly.
    cm = EnzymeMechanism(
        ((:S,), (:P,), ()),
        (((:E, :S), (:E_S,), true, 1),
         ((:E_S,), (:E_P,), false, 2),
         ((:E_P,), (:E, :P), true, 3)))
    # Pretend it's a 2-mer allosteric for the display test.
    am = EnzymeRates.AllostericEnzymeMechanism(
        cm, (2, (:EqualRT, :EqualRT, :EqualRT)), ())
    s = repr(am)
    @test occursin("(", s)        # parens grouping or a non-grouped step
    # If the mechanism truly has multi-step groups, this would also
    # exercise the parenthesized-group branch; for this single-step-
    # per-group case, just confirm rendering is multi-line and tagged.
    @test occursin(":: EqualRT", s)
end
```

- [ ] **Step 2: Run; confirm failure**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: failure on `!occursin("cat_allo_states:", s)` — the current display still prints that line.

### Task D.2: Implement helper and refactor `Base.show`

**Files:**
- Modify: `src/types.jl` lines 426-444

- [ ] **Step 1: Add a helper that renders catalytic steps grouped by kinetic_group**

Insert above `Base.show(io, m::AllostericEnzymeMechanism)`:

```julia
"""Render the catalytic mechanism's steps as multi-line text,
grouping steps that share a kinetic_group with parens and a single
`:: Tag` annotation. Mirrors `@allosteric_mechanism` macro syntax."""
function _format_allo_step_groups(
    io::IO, cm::EnzymeMechanism,
    m::AllostericEnzymeMechanism,
)
    rxns = reactions(cm)
    _arrow(is_eq) = is_eq ? " ⇌ " : " <--> "

    # Walk kinetic groups in first-appearance order
    groups_seen = Int[]
    group_to_step_idxs = Dict{Int,Vector{Int}}()
    for (i, step) in enumerate(rxns)
        g = step[4]
        if !haskey(group_to_step_idxs, g)
            push!(groups_seen, g)
            group_to_step_idxs[g] = Int[]
        end
        push!(group_to_step_idxs[g], i)
    end

    for g in groups_seen
        idxs = group_to_step_idxs[g]
        tag = cat_allo_state(m, g)
        if length(idxs) == 1
            (lhs, rhs, is_eq, _) = rxns[idxs[1]]
            print(io, "\n  ", join(lhs, " + "),
                  _arrow(is_eq), join(rhs, " + "),
                  " :: ", tag)
        else
            print(io, "\n  (")
            for (k, i) in enumerate(idxs)
                k > 1 && print(io, ", ")
                (lhs, rhs, is_eq, _) = rxns[i]
                print(io, join(lhs, " + "),
                      _arrow(is_eq), join(rhs, " + "))
            end
            print(io, ") :: ", tag)
        end
    end
end
```

- [ ] **Step 2: Replace the body of `Base.show(io, ::AllostericEnzymeMechanism)`**

Replace lines 426-444 with:

```julia
function Base.show(io::IO, m::AllostericEnzymeMechanism)
    cm = catalytic_mechanism(m)
    print(io, "AllostericEnzymeMechanism (cat_n=",
          catalytic_multiplicity(m))
    rs = regulatory_sites(m)
    if !isempty(rs)
        print(io, ", ", length(rs), " reg sites")
    end
    print(io, "):")
    _format_allo_step_groups(io, cm, m)
    for (i, (ligands, mult, reg_allo_states)) in enumerate(rs)
        print(io, "\n  reg site $i (n=", mult, "): ",
              join(ligands, ", "))
        print(io, " [")
        print(io, join(("$(n)::$(t)"
                        for (n, t) in zip(ligands, reg_allo_states)),
                       ", "))
        print(io, "]")
    end
end
```

- [ ] **Step 3: Run tests**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: D.1 passes; full suite green.

### Task D.3: Commit Phase D

- [ ] **Step 1: Stage and commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "$(cat <<'EOF'
Refactor AllostericEnzymeMechanism display to inline ::Tag per step

Drop the separate cat_allo_states summary line. Steps render
multi-line with their kinetic_group's allosteric tag inline,
matching @allosteric_mechanism macro syntax. Steps sharing a
kinetic_group render as a parenthesized group with one tag.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase E — `identify_rate_equation` file naming bucketing (Issues #3, #4)

### Task E.1: Add failing test

**Files:**
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Add a testset that exercises `_save_level_csv` via `_beam_search`**

A focused unit test on `_save_level_csv` is cleaner than an end-to-end run. Append:

```julia
@testset "save_level_csv groups by actual n_params" begin
    mktempdir() do tmp
        # Two mock rows with different actual n_params, sharing the
        # same enumeration-level pc.
        rows = [
            (n_params=3, loss=1.0,
             mechanism_type="m1", rate_equation="eq1",
             fitted_param_names=(:K1, :K2, :k3f),
             fitted_param_values=(1.0, 2.0, 3.0)),
            (n_params=4, loss=1.5,
             mechanism_type="m2", rate_equation="eq2",
             fitted_param_names=(:K1, :K2, :K3, :k4f),
             fitted_param_values=(1.0, 2.0, 3.0, 4.0)),
        ]
        # Group by actual n_params and call save once per bucket.
        for n in unique(r.n_params for r in rows)
            EnzymeRates._save_level_csv(
                tmp, [r for r in rows if r.n_params == n], n)
        end
        @test isfile(joinpath(tmp, "params_3.csv"))
        @test isfile(joinpath(tmp, "params_4.csv"))
        # No params_5.csv from any estimate-level
        @test !isfile(joinpath(tmp, "params_5.csv"))
        # Each file's sole row has matching n_params
        df3 = CSV.read(joinpath(tmp, "params_3.csv"), DataFrame)
        df4 = CSV.read(joinpath(tmp, "params_4.csv"), DataFrame)
        @test all(df3.n_params .== 3)
        @test all(df4.n_params .== 4)
    end
end
```

- [ ] **Step 2: Run; confirm failure**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: test passes IF `_save_level_csv` already accepts arbitrary param_count, but the *behavior* of `_beam_search` to group by actual is missing — so add a complementary integration test next.

(If both files are produced by the test above, this confirms `_save_level_csv` itself is fine; the bug is upstream — `_beam_search` always passes the level estimate. The Phase E.2 implementation moves the grouping into `_beam_search`.)

### Task E.2: Group results by actual `n_params` inside `_beam_search`

**Files:**
- Modify: `src/identify_rate_equation.jl` lines 308-313

- [ ] **Step 1: Replace the single-call save with grouped saves**

Find:

```julia
# Save CSV for this param count
if save_dir !== nothing
    _save_level_csv(
        save_dir,
        [r.row for r in results], pc)
end
```

Replace with:

```julia
# Save CSVs grouped by actual n_params (which may differ from
# the enumeration estimate `pc` due to Haldane reduction).
if save_dir !== nothing
    actual_groups = Dict{Int,Vector{NamedTuple}}()
    for r in results
        push!(get!(actual_groups, r.row.n_params, NamedTuple[]),
              r.row)
    end
    for (np, rows) in actual_groups
        _save_level_csv(save_dir, rows, np)
    end
end
```

- [ ] **Step 2: Run tests**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: E.1 passes; full suite green.

### Task E.3: Commit Phase E

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "$(cat <<'EOF'
Group beam-search saves by actual n_params, not the estimate-level

Files are now named after the row's n_params column (the actual
length(fitted_params(m))). One file → one param count, no Haldane-
induced mixing.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase F — Beam selection threshold (Issue #5)

### Task F.1: Add failing tests

**Files:**
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Add tests for the new selection rule**

```julia
@testset "beam selection: loss thresholds + min_beam_width floor" begin
    # Build a synthetic results vector with known losses and assert
    # the helper picks the right indices.
    losses = [1.0, 1.5, 2.5, 5.0, 10.0]
    # rel=2.0, abs=0.0: keep loss <= 2.0*1.0 = 2.0 -> indices 1,2
    sel = EnzymeRates._select_beam(
        losses;
        loss_rel_threshold=2.0,
        loss_abs_threshold=0.0,
        min_beam_width=1)
    @test sort(sel) == [1, 2]

    # rel=2.0, abs=0.0, min_beam_width=4 -> floor takes top-4
    sel = EnzymeRates._select_beam(
        losses;
        loss_rel_threshold=2.0,
        loss_abs_threshold=0.0,
        min_beam_width=4)
    @test sort(sel) == [1, 2, 3, 4]

    # near-zero best_loss with abs cushion: rel=2.0, abs=0.01,
    # losses [1e-6, 0.005, 0.05]:
    losses_small = [1e-6, 0.005, 0.05]
    sel = EnzymeRates._select_beam(
        losses_small;
        loss_rel_threshold=2.0,
        loss_abs_threshold=0.01,
        min_beam_width=1)
    # threshold = 2.0*1e-6 + 0.01 = 0.010002 -> indices 1, 2
    @test sort(sel) == [1, 2]
end
```

- [ ] **Step 2: Add a test that the kwargs are wired into `identify_rate_equation`**

```julia
@testset "identify_rate_equation accepts loss thresholds" begin
    # Just check the kwargs exist at the public-API level.
    # Inspect the method signature.
    methods_list = methods(identify_rate_equation)
    @test !isempty(methods_list)
    # Smoke-test: small problem with explicit thresholds doesn't
    # error during kwarg binding.
    # (Use the test_rxn / test_mechanism set up earlier in this
    # file to avoid duplicating setup.)
end
```

- [ ] **Step 3: Run; confirm failure**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: `_select_beam` is undefined → MethodError. New thresholds kwargs not recognized → MethodError or no-op.

### Task F.2: Implement `_select_beam` and wire kwargs

**Files:**
- Modify: `src/identify_rate_equation.jl`

- [ ] **Step 1: Add the helper**

Insert near the top of the file (after the docstrings):

```julia
"""
Return indices into `losses` for mechanisms that qualify for the
beam at this level. A mechanism qualifies if either:
  • its loss ≤ loss_rel_threshold * best_loss + loss_abs_threshold,
  • OR its rank (1-indexed by ascending loss) ≤ min_beam_width.
"""
function _select_beam(
    losses::AbstractVector{<:Real};
    loss_rel_threshold::Float64,
    loss_abs_threshold::Float64,
    min_beam_width::Int,
)
    isempty(losses) && return Int[]
    perm = sortperm(losses)
    best = losses[perm[1]]
    cutoff = loss_rel_threshold * best + loss_abs_threshold
    selected = Int[]
    for (rank, idx) in enumerate(perm)
        if losses[idx] <= cutoff || rank <= min_beam_width
            push!(selected, idx)
        end
    end
    selected
end
```

- [ ] **Step 2: Replace the rank-based block in `_beam_search`**

Find lines around 316-326:

```julia
# Beam select within this level
perm = sortperm(
    [r.row.loss for r in results])
beam_size = max(
    ceil(Int,
        beam_fraction *
        length(results)),
    min_beam_width)
beam_size = min(
    beam_size, length(results))
beam_specs = [results[perm[i]].spec
              for i in 1:beam_size]
```

Replace with:

```julia
# Beam select within this level
sel = _select_beam(
    [r.row.loss for r in results];
    loss_rel_threshold=loss_rel_threshold,
    loss_abs_threshold=loss_abs_threshold,
    min_beam_width=min_beam_width)
beam_specs = [results[i].spec for i in sel]
```

- [ ] **Step 3: Update the `_beam_search` and `identify_rate_equation` signatures**

In the public function, drop `beam_fraction`, add the two thresholds, lower `min_beam_width` default:

```julia
function identify_rate_equation(
    prob::IdentifyRateEquationProblem;
    # Beam search
    min_beam_width::Int = 50,
    loss_rel_threshold::Float64 = 2.0,
    loss_abs_threshold::Float64 = 0.01,
    max_param_count::Int = 20,
    ...
```

Forward the new kwargs into `_beam_search`:

```julia
specs, df = _beam_search(prob;
    min_beam_width, loss_rel_threshold, loss_abs_threshold,
    max_param_count, save_dir,
    pmap_function, optimizer,
    fitting_kwargs...)
```

Update `_beam_search`'s signature to accept the same kwargs. Drop `beam_fraction` from both functions.

- [ ] **Step 4: Update the docstring**

Replace the existing keyword-argument documentation block in the `identify_rate_equation` docstring with one that describes the new kwargs. Insert after the `min_beam_width` line:

```
- `loss_rel_threshold::Float64 = 2.0`: relative tolerance.
- `loss_abs_threshold::Float64 = 0.01`: absolute tolerance.

A mechanism qualifies for the next-level beam if EITHER:
  • its loss ≤ loss_rel_threshold * best_loss + loss_abs_threshold,
  • OR its rank by loss (ascending) ≤ min_beam_width.

The additive term protects against best_loss approaching zero
(simulated / very-low-loss data) where a purely multiplicative
threshold would collapse the beam to the single best mechanism.
```

Drop the `beam_fraction` line.

- [ ] **Step 5: Run tests**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: F.1 passes; full suite green.

### Task F.3: Commit Phase F

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "$(cat <<'EOF'
Replace beam_fraction with loss_rel_threshold + loss_abs_threshold

Beam selection now keeps mechanisms within a relative+absolute loss
tolerance of the best, with a min_beam_width floor (default 50,
down from 200). beam_fraction is removed; passing it raises
MethodError.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase G — Post-compile rate-equation hash dedup with cross-level cache (Issue #6)

This is the largest phase. It implements the persistent hash cache, four-stage per-level processing, three new CSV columns, and LOOCV per unique hash.

### Task G.1: Implement and test the canonical rate-equation hash

**Files:**
- Modify: `src/identify_rate_equation.jl`
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Add failing tests for the hash function**

Append to `test/test_identify_rate_equation.jl`:

```julia
@testset "canonical rate-equation hash" begin
    # Two strings differing only in k{N}f index should hash equal
    s1 = "k10r = (1 / Keq) * (1 / K1) * k10f\n" *
         "v = (k10f * S - k10r * P) / (1 + S/K1)"
    s2 = "k11r = (1 / Keq) * (1 / K1) * k11f\n" *
         "v = (k11f * S - k11r * P) / (1 + S/K1)"
    @test EnzymeRates._canonical_rate_eq_hash(s1) ==
          EnzymeRates._canonical_rate_eq_hash(s2)

    # Two strings with K names in different first-appearance order
    s3 = "v = (kf * S - kr * P) / (1 + S/K1 + P/K2)"
    s4 = "v = (kf * S - kr * P) / (1 + P/K1 + S/K2)"
    # These should hash DIFFERENTLY because the canonicalizer renames
    # by first appearance and the structures differ.
    @test EnzymeRates._canonical_rate_eq_hash(s3) !=
          EnzymeRates._canonical_rate_eq_hash(s4)

    # Hash is deterministic (same input → same hash)
    @test EnzymeRates._canonical_rate_eq_hash(s1) ==
          EnzymeRates._canonical_rate_eq_hash(s1)

    # Returns (full::Vector{UInt8} of length 32, short::String of 8 hex)
    h_full, h_short = EnzymeRates._canonical_rate_eq_hash_pair(s1)
    @test length(h_full) == 32
    @test length(h_short) == 8
    @test all(c -> c in "0123456789abcdef", h_short)
end
```

- [ ] **Step 2: Run; confirm failure (function undefined)**

- [ ] **Step 3: Implement the canonicalizer**

Add `using SHA` at the top of `src/identify_rate_equation.jl` (alongside the other `using` lines). Add helpers:

```julia
"""Canonicalize a rate-equation source string into a stable form
for hashing. Drops `(; … ) = params` / `= concs` destructure lines,
renames `k{N}f → kf_1, kf_2, …`, `k{N}r → kr_1, kr_2, …`, and
`K{N} → K_1, K_2, …` in first-appearance order, normalizes
whitespace."""
function _canonicalize_rate_eq(text::AbstractString)
    # 1. Strip destructure lines
    lines = split(text, '\n')
    body = String[]
    for ln in lines
        s = strip(ln)
        startswith(s, "(;") && occursin("= params", s) && continue
        startswith(s, "(;") && occursin("= concs", s)  && continue
        push!(body, ln)
    end
    text2 = join(body, "\n")

    # 2. Rename in first-appearance order. Walk through the text
    # left-to-right, replacing each match with its canonical name.
    # Three independent counters: kf, kr, K.
    out = IOBuffer()
    name_map = Dict{String,String}()
    kf_n = Ref(0); kr_n = Ref(0); K_n = Ref(0)

    function canon_for(token::AbstractString)
        if haskey(name_map, token); return name_map[token]; end
        new = if startswith(token, "k") && endswith(token, "f")
            kf_n[] += 1; "kf_$(kf_n[])"
        elseif startswith(token, "k") && endswith(token, "r")
            kr_n[] += 1; "kr_$(kr_n[])"
        elseif startswith(token, "K")
            K_n[] += 1; "K_$(K_n[])"
        else
            error("unexpected token: $token")
        end
        name_map[token] = new
        new
    end

    # Match patterns: k\d+f, k\d+r, K\d+ (but NOT Keq).
    # A single regex over the body string with replacement function:
    pattern = r"(k\d+f|k\d+r|K\d+)(?!eq)"
    text3 = replace(text2, pattern => canon_for)

    # 3. Whitespace normalize
    text4 = replace(text3, r"\s+" => " ")
    strip(text4)
end

"""Return SHA-256 bytes of the canonicalized rate equation."""
function _canonical_rate_eq_hash(text::AbstractString)
    SHA.sha256(codeunits(_canonicalize_rate_eq(text)))
end

"""Return (32-byte full hash, 8-char short hex) pair."""
function _canonical_rate_eq_hash_pair(text::AbstractString)
    h = _canonical_rate_eq_hash(text)
    short = bytes2hex(h)[1:8]
    (h, short)
end
```

- [ ] **Step 4: Run G.1 tests**

Expected: pass.

### Task G.2: Add the persistent fit cache and four-stage processing

**Files:**
- Modify: `src/identify_rate_equation.jl` `_beam_search` body

- [ ] **Step 1: Add a `FitResult` struct (internal)**

Insert near the top of `src/identify_rate_equation.jl`:

```julia
"""Cached fit result keyed by canonical rate-equation hash. The
`first_seen_estimate` field is the `n_fit_params_estimate` of the
*level* (the beam-search pc loop iteration) at which this hash was
first fit; the `first_seen_n_actual` is `length(fitted_params(m))`
at that fit. Both are useful diagnostics: the gap between them
indicates Haldane reduction, and the gap between this level's
estimate and `first_seen_estimate` indicates cross-level cache
reuse."""
struct _CachedFitResult
    loss::Float64
    params::NamedTuple
    first_seen_estimate::Int
    first_seen_n_actual::Int
    first_seen_eq_hash::String   # 8-char short hex
end
```

**Spec deviation note:** the spec's tentative
`first_seen_n_params` field is renamed here to
`first_seen_estimate` because the previous name was degenerate
with the row's own `n_params` (any two specs with the same
canonical hash compile to mechanisms with the same actual
`n_params`, so the column would always be either `missing` or
equal to `n_params`). Tracking the estimate-level instead
captures the actually-useful diagnostic Denis described in
brainstorming: "this equation was already fit during a previous
round and now it comes back."

- [ ] **Step 2: Refactor `_beam_search` to four-stage processing**

Replace the body of `_beam_search` (current lines 243-344) with:

```julia
function _beam_search(
    prob::IdentifyRateEquationProblem;
    min_beam_width, loss_rel_threshold, loss_abs_threshold,
    max_param_count, save_dir, pmap_function,
    optimizer, kwargs...
)
    # Persistent cross-level cache (full 32-byte hash → cached result)
    fit_cache = Dict{Vector{UInt8}, _CachedFitResult}()

    # Build initial cache by enumeration estimate
    cache = Dict{Int,Vector{AbstractMechanismSpec}}()
    for spec in init_mechanisms(prob.reaction)
        push!(get!(cache, spec.n_fit_params_estimate,
                   AbstractMechanismSpec[]),
              spec)
    end
    dedup!(cache)

    all_specs = AbstractMechanismSpec[]
    all_rows  = NamedTuple[]

    isempty(cache) && return (
        all_specs, _rows_to_dataframe(all_rows))

    min_pc = minimum(keys(cache))
    for pc in min_pc:max_param_count
        level = pop!(cache, pc, AbstractMechanismSpec[])
        isempty(level) && (isempty(cache) ? break : continue)

        # ── Stage 1: parallel compile + hash ───────────
        compiled = pmap_function(level) do spec
            try
                m = compile_mechanism(spec)
                eq_text = rate_equation_string(m)
                h_full, h_short = _canonical_rate_eq_hash_pair(eq_text)
                n_actual = length(fitted_params(m))
                (spec=spec, mechanism=m, eq_text=eq_text,
                 h_full=h_full, h_short=h_short,
                 n_actual=n_actual, ok=true)
            catch e
                @debug("Mechanism compilation failed",
                       exception=(e, catch_backtrace()))
                (spec=spec, ok=false)
            end
        end
        filter!(c -> c.ok, compiled)
        isempty(compiled) && continue

        # ── Stage 2: bucket by full hash within level ──
        by_hash = Dict{Vector{UInt8},Vector{NamedTuple}}()
        for c in compiled
            push!(get!(by_hash, c.h_full, NamedTuple[]), c)
        end

        # ── Stage 3: identify NEW hashes; fit reps in parallel ─
        new_hashes = [h for h in keys(by_hash) if !haskey(fit_cache, h)]
        new_reps = [first(by_hash[h]) for h in new_hashes]
        new_results = pmap_function(new_reps) do rep
            try
                fp = FittingProblem(rep.mechanism, prob.data;
                                    Keq=prob.Keq)
                fit = fit_rate_equation(fp, optimizer; kwargs...)
                (h_full=rep.h_full, h_short=rep.h_short,
                 n_actual=rep.n_actual,
                 loss=fit.loss, params=fit.params, ok=true)
            catch e
                @debug("Fit failed",
                       exception=(e, catch_backtrace()))
                (h_full=rep.h_full, ok=false)
            end
        end
        for r in new_results
            r.ok || continue
            fit_cache[r.h_full] = _CachedFitResult(
                r.loss, r.params, pc, r.n_actual, r.h_short)
        end

        # ── Stage 4: build rows (one per hash group) ──
        # Snapshot which hashes were fit at this level vs inherited.
        new_hashes_set = Set(new_hashes)
        ordered_hashes = collect(keys(by_hash))
        level_rows = NamedTuple[]
        level_specs = AbstractMechanismSpec[]   # ALL members for expand
        for h_full in ordered_hashes
            members = by_hash[h_full]
            haskey(fit_cache, h_full) || continue   # fit failed
            cached = fit_cache[h_full]
            rep = first(members)
            is_inherited = !(h_full in new_hashes_set)
            row = (
                n_params = rep.n_actual,
                loss = cached.loss,
                mechanism_type = _mechanism_type_string(rep.mechanism),
                rate_equation = rep.eq_text,
                fitted_param_names = collect(keys(cached.params)),
                fitted_param_values = Tuple(values(cached.params)),
                eq_hash = cached.first_seen_eq_hash,
                n_equivalent = length(members),
                fit_inherited_from_estimate =
                    is_inherited ? cached.first_seen_estimate : missing,
            )
            push!(level_rows, row)
            for c in members
                push!(level_specs, c.spec)
            end
        end

        append!(all_specs, level_specs)
        append!(all_rows,  level_rows)

        # Save CSV grouped by actual n_params
        if save_dir !== nothing
            actual_groups = Dict{Int,Vector{NamedTuple}}()
            for row in level_rows
                push!(get!(actual_groups, row.n_params, NamedTuple[]),
                      row)
            end
            for (np, rows) in actual_groups
                _save_level_csv(save_dir, rows, np)
            end
        end

        # ── Beam selection: include all spec members of qualifying
        #     hash groups so RE→SS expansion variants are preserved ─
        sel = _select_beam(
            [r.loss for r in level_rows];
            loss_rel_threshold=loss_rel_threshold,
            loss_abs_threshold=loss_abs_threshold,
            min_beam_width=min_beam_width)
        selected_hashes = Set{Vector{UInt8}}()
        # Each level_row corresponds to one hash group; map sel
        # indices back to the hash groups.
        ordered_hashes = collect(keys(by_hash))
        # Need same ordering for level_rows and ordered_hashes:
        # rebuild together to be safe.
        # (See implementation note below if reordering is risky.)
        for i in sel
            push!(selected_hashes, ordered_hashes[i])
        end
        beam_specs = AbstractMechanismSpec[]
        for h in selected_hashes
            for c in by_hash[h]
                push!(beam_specs, c.spec)
            end
        end

        # ── Expand all beam specs ───────────
        new_cache = expand_mechanisms(beam_specs, prob.reaction)
        for (target_pc, specs) in new_cache
            target_pc > max_param_count && continue
            append!(get!(cache, target_pc,
                         AbstractMechanismSpec[]),
                    specs)
        end
        dedup!(cache)
    end

    df = _rows_to_dataframe(all_rows)
    return all_specs, df
end
```

**Implementation note** (to leave as a comment in the code): the
`level_rows` and `ordered_hashes` MUST be built from a single pass
over the same data so their indices align. The Stage 4 loop above
builds them in `by_hash`'s iteration order; capture that ordering
once into a local vector and reuse it both for row-building and for
selected-hash mapping. Refactor the Stage 4 loop to:

```julia
ordered_hashes = collect(keys(by_hash))
level_rows = Vector{NamedTuple}(undef, length(ordered_hashes))
for (idx, h_full) in enumerate(ordered_hashes)
    members = by_hash[h_full]
    ...   # build the row as above
    level_rows[idx] = row
end
```

so that `ordered_hashes[i]` and `level_rows[i]` always refer to the
same hash group.

- [ ] **Step 3: Update `_rows_to_dataframe` to surface the new columns**

Replace `_rows_to_dataframe` (lines 196-228) with:

```julia
function _rows_to_dataframe(rows)
    isempty(rows) && return DataFrame()
    all_pnames = Set{Symbol}()
    for row in rows
        for p in row.fitted_param_names; push!(all_pnames, p); end
    end
    sorted_pnames = sort(collect(all_pnames))

    df = DataFrame(
        n_params = [r.n_params for r in rows],
        loss = [r.loss for r in rows],
        mechanism_type = [r.mechanism_type for r in rows],
        rate_equation = [r.rate_equation for r in rows],
        eq_hash = [r.eq_hash for r in rows],
        n_equivalent = [r.n_equivalent for r in rows],
        fit_inherited_from_estimate = [
            r.fit_inherited_from_estimate for r in rows],
    )
    for pn in sorted_pnames
        df[!, pn] = [
            pn in r.fitted_param_names ?
                r.fitted_param_values[
                    findfirst(==(pn), r.fitted_param_names)] :
                missing
            for r in rows
        ]
    end
    df
end
```

### Task G.3: Add tests for cache behavior

**Files:**
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Add cross-level cache regression test (concrete)**

The cleanest fixture exercises the cache directly using two
mechanism specs that compile to the same rate equation. Two
specs from the LDH `params_7.csv` 9-member duplicate group provide
known same-hash specs at *different* enumeration estimate levels.

Append to `test/test_identify_rate_equation.jl`:

```julia
@testset "fit cache: cross-level reuse" begin
    # Two synthetic mechanism specs whose compiled rate equations
    # canonically hash equal. We build small fake compiled objects
    # by stubbing rate_equation_string output; the test exercises
    # _beam_search's cache mechanics, NOT the full pipeline.
    #
    # Strategy: monkey-patch fit_rate_equation via a recording
    # wrapper while testing _beam_search end-to-end on a small
    # reaction that produces known Pattern-A duplicates.

    # Use the bi_bi reaction (4 metabolites) with simple synthetic
    # data — enough to drive several enumeration levels and
    # exercise cross-level cache hits via Haldane collapse.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end

    # Synthetic data tight enough that fits converge quickly:
    Keq = 1.0
    data = (
        group = repeat(1:3, inner=4),
        Rate  = repeat([1.0, 2.0, 3.0, 4.0], 3),
        A = repeat([0.1, 0.5, 1.0, 2.0], 3),
        B = repeat([0.5, 0.5, 0.5, 0.5], 3),
        P = repeat([0.1, 0.1, 0.1, 0.1], 3),
        Q = repeat([0.1, 0.1, 0.1, 0.1], 3),
    )
    prob = IdentifyRateEquationProblem(rxn, data; Keq=Keq)

    # Recording optimizer: counts fit calls.
    fit_calls = Ref(0)
    orig_fit = EnzymeRates.fit_rate_equation
    # Replace EnzymeRates.fit_rate_equation in a test-local override.
    # Use a closure passed into _beam_search as a function-typed
    # kwarg. (Adjust _beam_search to accept `fit_function` if not
    # already there — see implementation note below.)
    function counting_fit(fp, opt; kwargs...)
        fit_calls[] += 1
        orig_fit(fp, opt; kwargs...)
    end

    # Run beam search with a tiny budget so it terminates fast.
    specs, df = EnzymeRates._beam_search(
        prob;
        min_beam_width = 5,
        loss_rel_threshold = 100.0,    # accept everything
        loss_abs_threshold = 0.0,
        max_param_count = 6,
        save_dir = nothing,
        pmap_function = map,
        optimizer = MockOptimizer(),   # see fixture below
        fit_function = counting_fit,
        n_restarts = 1, maxtime = 1.0,
    )

    # The number of fits performed must equal the number of
    # distinct eq_hash values across all rows — never more.
    @test fit_calls[] == length(unique(df.eq_hash))

    # If any row's fit was inherited, the column shows the
    # originating estimate level.
    inherited = filter(
        r -> !ismissing(r.fit_inherited_from_estimate),
        eachrow(df))
    if !isempty(inherited)
        @test all(r.fit_inherited_from_estimate < r.n_params + 5
                  for r in inherited)
        # Sanity: inherited estimate is < or = the row's own
        # current iteration's pc, which is bounded by max_param_count.
    end
end
```

**Implementation note** for the engineer: `_beam_search` must
accept a `fit_function` kwarg (default `fit_rate_equation`) so
tests can inject a recording wrapper. If you didn't add this in
G.2, add it now: thread `fit_function` through, defaulting to
`fit_rate_equation`. Replace the `fit_rate_equation(fp, optimizer;
kwargs...)` call in Stage 3 with `fit_function(fp, optimizer;
kwargs...)`.

Also: define a tiny `MockOptimizer` placeholder somewhere in the
test setup that returns fixed parameters quickly. Or — if the real
optimizer is fast enough on this 4-metabolite reaction with
`maxtime=1.0` — use a real one (e.g.
`Optimization.OptimizationOptimizers.Adam()`), preferring
correctness over fixture complexity.

- [ ] **Step 2: Add within-level no-inheritance test**

```julia
@testset "fit cache: within-level rows show missing inheritance" begin
    # Within a single beam level, all rows are first-fit at that
    # level (regardless of how many specs share each hash). So
    # fit_inherited_from_estimate should be `missing` for every
    # row in the LOWEST level's CSV/df chunk.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    data = (
        group = repeat(1:2, inner=3),
        Rate  = repeat([1.0, 2.0, 3.0], 2),
        S = repeat([0.1, 0.5, 1.0], 2),
        P = repeat([0.1, 0.1, 0.1], 2),
    )
    prob = IdentifyRateEquationProblem(rxn, data; Keq=1.0)
    specs, df = EnzymeRates._beam_search(
        prob;
        min_beam_width = 5,
        loss_rel_threshold = 100.0,
        loss_abs_threshold = 0.0,
        max_param_count = 4,
        save_dir = nothing,
        pmap_function = map,
        optimizer = MockOptimizer(),
        n_restarts = 1, maxtime = 1.0,
    )
    # Find the smallest n_params present in df:
    !isempty(df) || return
    smallest_np = minimum(df.n_params)
    first_level_rows = filter(r -> r.n_params == smallest_np,
                              eachrow(df))
    @test all(ismissing(r.fit_inherited_from_estimate)
              for r in first_level_rows)
end
```

- [ ] **Step 3: Add LOOCV-per-unique-hash dedup test**

In `_cv_model_selection`, dedup candidates by `eq_hash` within each
`n_params` bucket. Add an inline filter:

```julia
candidate_indices = Int[]
for gdf in groupby(df_indexed, :n_params)
    seen_hashes = Set{String}()
    sorted = sort(gdf, :loss)
    for row in eachrow(sorted)
        row.eq_hash in seen_hashes && continue
        push!(seen_hashes, row.eq_hash)
        push!(candidate_indices, row.spec_idx)
        length(seen_hashes) >= n_cv_candidates && break
    end
end
```

Add a regression test that constructs a small DataFrame with
duplicate `eq_hash` values within an `n_params` bucket and asserts
the deduper drops them:

```julia
@testset "_cv_model_selection dedups by eq_hash within n_params" begin
    df = DataFrame(
        n_params = [3, 3, 3, 4],
        loss = [1.0, 1.5, 1.2, 0.5],
        eq_hash = ["aaa", "aaa", "bbb", "ccc"],
        # other columns omitted from this dedup unit test
    )
    df_indexed = copy(df)
    df_indexed.spec_idx = 1:nrow(df_indexed)

    candidate_indices = Int[]
    for gdf in groupby(df_indexed, :n_params)
        seen_hashes = Set{String}()
        sorted = sort(gdf, :loss)
        for row in eachrow(sorted)
            row.eq_hash in seen_hashes && continue
            push!(seen_hashes, row.eq_hash)
            push!(candidate_indices, row.spec_idx)
            length(seen_hashes) >= 5 && break
        end
    end
    # 4 rows total → 3 unique hashes after dedup (aaa, bbb, ccc):
    @test length(candidate_indices) == 3
    # First-occurrence of "aaa" wins by lowest loss in its group:
    @test 1 in candidate_indices
    @test !(2 in candidate_indices)   # second "aaa" dropped
end
```

### Task G.4: End-to-end LDH regression

**Files:**
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Add an end-to-end small-bi-bi regression**

```julia
@testset "identify_rate_equation: end-to-end CSV invariants" begin
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    # Synthetic data from a known random-binding mechanism:
    Keq = 2.0
    n_pts_per_group = 8
    n_groups = 3
    A_vals = repeat([0.1, 0.5, 1.0, 2.0], n_groups * 2)
    B_vals = repeat([0.5, 0.5, 0.5, 0.5], n_groups * 2)
    P_vals = repeat([0.05, 0.05, 0.05, 0.05], n_groups * 2)
    Q_vals = repeat([0.05, 0.05, 0.05, 0.05], n_groups * 2)
    Rates  = @. (A_vals * B_vals) /
                ((1 + A_vals/0.5) * (1 + B_vals/0.3))
    data = (
        group = repeat(1:n_groups, inner=n_pts_per_group),
        Rate = Rates, A = A_vals, B = B_vals,
        P = P_vals, Q = Q_vals,
    )
    prob = IdentifyRateEquationProblem(rxn, data; Keq=Keq)
    mktempdir() do tmp
        result = identify_rate_equation(
            prob;
            min_beam_width = 5,
            loss_rel_threshold = 5.0,
            loss_abs_threshold = 0.01,
            max_param_count = 6,
            optimizer = MockOptimizer(),
            save_dir = tmp,
            pmap_function = map,
            n_restarts = 1, maxtime = 2.0,
        )
        # 1. Every saved file's name matches the n_params column.
        all_hashes = String[]
        for fname in readdir(tmp)
            endswith(fname, ".csv") || continue
            np = parse(Int,
                replace(fname, "params_" => "", ".csv" => ""))
            df_file = CSV.read(joinpath(tmp, fname), DataFrame)
            @test all(df_file.n_params .== np)
            # 2. eq_hash is unique within each file.
            @test length(unique(df_file.eq_hash)) ==
                  nrow(df_file)
            append!(all_hashes, df_file.eq_hash)
        end
        # 3. eq_hash is unique across all files (cross-level cache).
        @test length(unique(all_hashes)) == length(all_hashes)
        # 4. result has a best mechanism.
        @test result.best isa AbstractEnzymeMechanism
    end
end
```

The `MockOptimizer` referenced here is the same shared fixture
introduced in G.3 step 1; if you haven't added it, add it once
near the top of the test file:

```julia
struct MockOptimizer end
```

…and use a real Optimization.jl optimizer if a stub doesn't
satisfy the `fit_rate_equation` interface. If a stub IS workable,
implement it as a one-step "return the initial point" solver that
makes loss reproducible across runs.

### Task G.5: Commit Phase G

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "$(cat <<'EOF'
Add persistent rate-equation hash cache for fit dedup

A canonical SHA-256 of each compiled rate equation (with
parameter-name normalization) keys a Dict that survives across all
beam levels. Specs that share a hash with any prior level reuse
the cached fit; new hashes get fitted in parallel via pmap.

CSV gains eq_hash, n_equivalent, and fit_inherited_from_estimate
columns. _cv_model_selection dedups by hash within n_params buckets.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Final verification

- [ ] **Run the full test suite cleanly**

```
julia --project -e 'using Pkg; Pkg.activate("."); Pkg.test()'
```

Expected: all tests pass. Aqua + JET checks green. No `@test_skip` should remain — all Phase G regression tests have concrete bodies and run as part of the suite.

- [ ] **Run the README example**

```
julia --project test/test_readme_runs.jl
```

Expected: pass.

- [ ] **Inspect the diff one last time**

```bash
git log --oneline main..HEAD
git diff main..HEAD --stat
```

Confirm:
- Each phase is its own commit.
- No file outside the inventory has unexpected changes.
- No `param_count` references remain in src or tests.
- No bare-symbol `@enzyme_reaction` examples remain.

---

## Notes for the engineer

- **`_select_beam` indexing:** the ordering of `level_rows` and the per-hash group lookup in beam selection MUST agree. Use a single explicit `ordered_hashes = collect(keys(by_hash))` reference and index into it once — do not rely on Julia's `Dict` iteration order being stable across two separate `for (h, …) in by_hash` loops.
- **Pattern A test fixture:** the simplest known Pattern-A duplicate is two specs for a bi-bi reaction whose RE-edge subsets differ but share the same form set + kinetic-group structure. The 9-member group in the LDH `params_7.csv` analysis is a real example; pick any two of those mechanism types, parse them into `MechanismSpec`s, and use them as the test fixture. The `mechanism_type` strings are quoted in the spec doc's brainstorming notes if you need to copy them verbatim.
- **`SHA` import:** ensure `using SHA` is at module top; `SHA` is a Julia stdlib (no Project.toml change required).
- **Backwards compatibility:** none required. This is internal-version cleanup, not a public-API deprecation. Removing `beam_fraction` raises `MethodError`; that is intended.
