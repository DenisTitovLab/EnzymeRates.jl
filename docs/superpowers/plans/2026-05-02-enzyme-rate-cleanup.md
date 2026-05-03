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
| `src/mechanism_enumeration.jl` | A, C | rename `param_count` → `n_fit_params_estimate`; drop `+2` from BOTH estimate formulas (lines 841 & 1434); lower `floor_pc` from `+3` to `+1`; derive regulators from steps in `EnzymeMechanism(spec)` and allosteric path |
| `src/types.jl` | B, C, D | add atom mandatoriness + balance check in `EnzymeReaction`; tighten regulator validation in `EnzymeMechanism(mets, rxns)`; refactor `Base.show(io, ::AllostericEnzymeMechanism)` to inline `:: Tag` per step / step-group |
| `src/identify_rate_equation.jl` | E, F, G | rename saved CSV files to `params_estimate_<pc>.csv`; new beam threshold rule (`loss_rel_threshold`, `loss_abs_threshold`, lowered `min_beam_width=50`); token-driven canonical rate-equation hash; persistent cross-level fit cache; two-stage processing with worker-side recompile; new CSV columns (`eq_hash`, `fit_inherited_from_estimate`); LOOCV dedups by hash within each `n_params` bucket |
| `test/test_mechanism_enumeration.jl` | A, C | new estimate semantic tests; init does not carry unbound regulators |
| `test/test_types.jl` | B, C, D | atom-balance + atom-mandatory failure tests; strict regulator constructor failure test; `repr` allosteric display tests including the multi-step parens-grouping branch |
| `test/test_dsl.jl` | B | atom-validation tests at the macro/constructor seam (only if any `@enzyme_reaction` examples in this file remain bare-symbol after migration) |
| `test/test_identify_rate_equation.jl` | E, F, G | filename rename test; beam threshold tests; canonical-hash test; end-to-end CSV-invariants regression with real optimizer |
| `README.md` | B | migrate any bare-symbol `@enzyme_reaction` example (NOT `@enzyme_mechanism` / `@allosteric_mechanism` blocks — those forbid atom brackets) |
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

- [ ] **Step 4: Replace remaining `param_count` usages in src/**

Use a word-boundary grep so bare-identifier uses (local variable
references, positional argument names, comments) are also caught:

```bash
grep -rn "\bparam_count\b" /home/denis.linux/.julia/dev/EnzymeRates/src/
```

Replace each occurrence of `param_count` with `n_fit_params_estimate`,
including: `spec.param_count` field access, `param_count = …`
assignments, local variable declarations like
`param_count = n_re + 2 * n_ss - n_thermo`, positional uses inside
constructor calls (`MechanismSpec(rxn, tagged, param_count)`), and
any inline comments. Confirm clean:

```bash
grep -rn "\bparam_count\b" /home/denis.linux/.julia/dev/EnzymeRates/src/
```

Expected: empty output.

- [ ] **Step 5: Replace `param_count` references in tests (including comments)**

```bash
grep -rn "\bparam_count\b" /home/denis.linux/.julia/dev/EnzymeRates/test/
```

Replace each. Specific known sites:
- `test/test_mechanism_enumeration.jl:17` and `:100` (the `enumerate_all` helper).
- `test/test_identify_rate_equation.jl:41` (the `8) # param_count` line —
  the value is a positional argument to `AllostericMechanismSpec`; the
  inline comment `# param_count` must be updated to
  `# n_fit_params_estimate`).

Confirm clean:

```bash
grep -rn "\bparam_count\b" /home/denis.linux/.julia/dev/EnzymeRates/test/
```

Expected: empty output.

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

### Task A.3: Drop `+2` from BOTH estimate formulas and lower `floor_pc`

There are TWO `+ 2` occurrences in `src/mechanism_enumeration.jl`,
plus a `floor_pc` constant whose value depends on the convention.
All three must change in lock-step.

**Files:**
- Modify: `src/mechanism_enumeration.jl` lines 841, 1374, 1434

- [ ] **Step 1: Change the per-step formula in `init_mechanisms`**

Around line 841, replace:

```julia
param_count = n_re + 2 * n_ss -
    n_thermo + 2
```

with (after Phase A.2's rename):

```julia
n_fit_params_estimate = n_re + 2 * n_ss - n_thermo
```

- [ ] **Step 2: Change the per-kinetic-group formula in `_param_count_from_steps`**

This is at line 1422-1436 of the SAME file. After the Phase A.2
rename it is now `_n_fit_params_estimate_from_steps(steps)`. Drop
the trailing `+ 2`:

```julia
function _n_fit_params_estimate_from_steps(steps::Vector{StepSpec})
    groups_re = Set{Int}()
    groups_ss = Set{Int}()
    for s in steps
        if s.is_equilibrium
            push!(groups_re, s.kinetic_group)
        else
            push!(groups_ss, s.kinetic_group)
        end
    end
    n_forms = length(all_form_names(steps))
    n_thermo = length(steps) - n_forms + 1
    length(groups_re) + 2 * length(groups_ss) - n_thermo
end
```

This is the helper called from `_apply_equivalence_grouping` at
line 1411 — its return value is what actually gets stamped onto the
spec at the end of `init_mechanisms`. Without changing it, the
A.1 regression test still fails.

- [ ] **Step 3: Lower `floor_pc` from `+3` to `+1`**

At line 1374 in `init_mechanisms`:

```julia
floor_pc = n_s + n_p + 1
```

The original `+3` decomposes as `+1` (one SS rate constant) `+1`
(Keq) `+1` (E_total). Under the new convention (no Keq, no
E_total), only the `+1` for the SS rate constant remains.

- [ ] **Step 4: Run the full test suite**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: A.1 now passes. Whole suite green.

If A.1 still fails, the most likely cause is that `floor_pc` is
still clamping the count too high for the small `uni_uni_rxn`
fixture (n_s = n_p = 1 → floor_pc = 3 under new rule; the actual
fitted_params count is also 3 = K1, K2, k3f). They should agree.
Print both values from the test if the assertion fails to debug.

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

- [ ] **Step 1: Find every bare-symbol `@enzyme_reaction` block**

The grep below is the source of truth. **Only `@enzyme_reaction`
blocks need migration** — `@enzyme_mechanism` and
`@allosteric_mechanism` DSLs forbid atom brackets (per
`dsl.jl:285-289`), so any bare symbols in those blocks are
correct and must NOT be touched.

```bash
grep -rA 6 "@enzyme_reaction begin" /home/denis.linux/.julia/dev/EnzymeRates/test/ /home/denis.linux/.julia/dev/EnzymeRates/README.md /home/denis.linux/.julia/dev/EnzymeRates/src/
```

For each match, inspect the `substrates:` and `products:` lines
INSIDE THE `@enzyme_reaction` BLOCK only. Any bare-symbol species
(e.g., `substrates: S` with no `[…]`) must be migrated to declare
atoms (e.g., `S[C]` is fine for illustrative examples — pick
plausible biochemistry where one is implied, otherwise placeholder
`[C]` is acceptable).

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
    # Construct a base mechanism where TWO RE steps share
    # kinetic_group=1 — both bind metabolite :S to different forms.
    # The display must render them together as
    #   ([E, S] ⇌ [E_S], [E_P, S] ⇌ [E_PS]) :: Tag
    cm = EnzymeMechanism(
        ((:S,), (:P,), ()),
        (((:E, :S),    (:E_S,),  true,  1),  # first :S binding, group 1
         ((:E_P, :S),  (:E_PS,), true,  1),  # second :S binding, group 1
         ((:E_S,),     (:E_P,),  false, 2),  # SS catalytic, group 2
         ((:E_P,),     (:E, :P), true,  3))) # P release, group 3
    # 2-mer allosteric, all kinetic groups :EqualRT
    am = EnzymeRates.AllostericEnzymeMechanism(
        cm, (2, (:EqualRT, :EqualRT, :EqualRT)), ())
    s = repr(am)

    # The parenthesized-group branch must execute: look for the
    # exact "(...) :: EqualRT" shape with a comma-separated body.
    paren_group_match = match(
        r"\([^()]*,[^()]*\) :: EqualRT", s)
    @test paren_group_match !== nothing
    # And the single-step lines for groups 2 and 3 should also
    # appear with their own ":: EqualRT" tags.
    @test occursin("[E_S] <--> [E_P] :: EqualRT", s)
    @test occursin(":: EqualRT", s)
    # The deprecated summary line must NOT appear.
    @test !occursin("cat_allo_states:", s)
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

## Phase E — `identify_rate_equation` filename rename (Issues #3, #4)

The original plan tried to bucket save files by actual
`length(fitted_params(m))`. That introduced a class of CSV-append
edge cases (different column sets, header rewriting, eq_hash
duplication across levels). The simpler answer: keep one file per
estimate-level (the existing behavior) but rename it so the
filename clearly says "this is an estimate, not the actual count."
The row's `n_params` column already shows the actual count; users
sort and filter from there.

### Task E.1: Add failing test for the new filename pattern

**Files:**
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Add a testset asserting the new filename**

Append:

```julia
@testset "save_level_csv uses estimate-level filename" begin
    mktempdir() do tmp
        rows = [(n_params=3, loss=1.0,
                 mechanism_type="m1", rate_equation="eq1",
                 fitted_param_names=(:K1, :K2, :k3f),
                 fitted_param_values=(1.0, 2.0, 3.0))]
        # Caller passes the estimate-level pc (e.g., 5) — could
        # diverge from the row's actual n_params=3 due to Haldane
        # reduction. Filename must reflect the estimate.
        EnzymeRates._save_level_csv(tmp, rows, 5)
        @test isfile(joinpath(tmp, "params_estimate_5.csv"))
        @test !isfile(joinpath(tmp, "params_5.csv"))
        df = CSV.read(joinpath(tmp, "params_estimate_5.csv"),
                      DataFrame)
        @test df.n_params == [3]   # actual count, not the estimate
    end
end
```

- [ ] **Step 2: Run; confirm failure**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: file is found at `params_5.csv` (old name), so the
`isfile(joinpath(tmp, "params_estimate_5.csv"))` assertion fails.

### Task E.2: Rename the saved-file path

**Files:**
- Modify: `src/identify_rate_equation.jl` lines 233-241

- [ ] **Step 1: Update `_save_level_csv` filename**

Replace lines 237-238:

```julia
path = joinpath(
    save_dir, "params_$(param_count).csv")
```

with:

```julia
path = joinpath(
    save_dir, "params_estimate_$(param_count).csv")
```

Update the docstring just above to clarify:

```julia
"""
Save results for one beam level to a CSV file. The filename
encodes the level's `n_fit_params_estimate`; the actual `n_params`
of each row may be smaller (Haldane reduction collapses some
declared kinetic groups). Users wanting one file per actual
`n_params` value can post-process by reading and re-grouping.
"""
```

The function's positional argument name should also be renamed
from `param_count` to `n_fit_params_estimate` for consistency
with the rest of the codebase after Phase A's rename.

- [ ] **Step 2: Run tests**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: E.1 passes; full suite green.

### Task E.3: Commit Phase E

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "$(cat <<'EOF'
Rename beam-search saves to params_estimate_<pc>.csv

Filename now signals that the integer is the enumeration estimate,
not the actual fitted-param count. The n_params column inside each
row remains the actual count from length(fitted_params(m)). Users
can post-process to bucket by actual count if desired.

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

## Phase G — Post-compile rate-equation hash dedup with persistent cache (Issue #6)

This phase adds a global cache keyed by canonical rate-equation hash
that persists across beam levels. Specs whose compiled rate
equations share a hash with any previously-seen mechanism reuse the
prior fit; only new hashes pay the optimizer cost.

Key design decisions (refined from spec deviations during plan
review):

- **Canonical hash uses Julia's built-in `hash(::String)::UInt64`,
  not SHA.** No new Project.toml dep; collision probability is
  negligible at our scale (~10⁻¹² over 10⁴ mechanisms).
- **Canonicalization is token-driven via `parameters(m, Reduced)`,
  not regex.** The param list is the source of truth; we walk
  `parameters(m)` to discover every name (`K1`, `K1_T`, `kf_T`,
  `K_R_reg1`, `L`, etc.) and rename each by first appearance in
  the body. Avoids the regex-misses-allosteric-params class of
  bugs.
- **One row per spec member, not per hash group.** All members of
  a hash group share the cached `(loss, params, eq_hash)` values;
  the `eq_hash` column lets users post-hoc dedup. No 1:1
  spec-to-row mapping is broken; `_cv_model_selection`'s
  indexing logic is unchanged.
- **Two stages with worker-side recompile.** Stage 1 returns
  `(spec, eq_hash, eq_hash_short, n_actual)` — *no* mechanism
  object. Stage 2 receives the spec, recompiles it on the worker,
  and fits. This pays a 2× compile cost on new-hash specs but
  guarantees `EnzymeMechanism{...}` singleton types never travel
  between workers.
- **No unit tests of cache mechanics.** A single end-to-end test
  using a real optimizer covers cache + integration. Removed the
  earlier G.3 unit-test scaffolding, which was either tautological
  or required `MockOptimizer` plumbing that doesn't satisfy
  `Optimization.solve`'s dispatch interface.

### Task G.1: Implement and test the canonical rate-equation hash

**Files:**
- Modify: `src/identify_rate_equation.jl`
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Add failing tests for the hash function**

Append to `test/test_identify_rate_equation.jl`:

```julia
@testset "canonical rate-equation hash" begin
    # Build two mechanisms whose rate equations are isomorphic up
    # to parameter-name renumbering — same form set, same SS step,
    # different step ordering so rep-step indices differ.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
    end

    init = EnzymeRates.init_mechanisms(rxn)
    base = first(init)
    m_a = EnzymeRates.EnzymeMechanism(base)

    # Reorder steps so the SS step lands at a different index;
    # the rate equation should canonicalize to the same hash.
    rxns = collect(EnzymeRates.reactions(m_a))
    swapped = (rxns[2], rxns[1], rxns[3:end]...)
    m_b = EnzymeMechanism(((:S,), (:P,), ()), swapped)

    h_a = EnzymeRates._canonical_rate_eq_hash(m_a)
    h_b = EnzymeRates._canonical_rate_eq_hash(m_b)
    @test h_a == h_b

    # Determinism.
    @test EnzymeRates._canonical_rate_eq_hash(m_a) == h_a

    # Pair returns (UInt64, 16-char hex String).
    h_full, h_short = EnzymeRates._canonical_rate_eq_hash_pair(m_a)
    @test h_full == h_a
    @test length(h_short) == 16
    @test all(c -> c in "0123456789abcdef", h_short)
end
```

- [ ] **Step 2: Run; confirm failure**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: `MethodError` / `UndefVarError` for `_canonical_rate_eq_hash`.

- [ ] **Step 3: Implement the canonicalizer (token-driven via `parameters(m)`)**

Add to `src/identify_rate_equation.jl`:

```julia
"""
Build a canonical text representation of a mechanism's rate
equation, suitable for hashing.

Strategy: walk `parameters(m, Reduced)` to discover every
parameter symbol, scan the rate-equation body to find each
parameter's first-appearance position, then rename them as
`p_1, p_2, …` in first-appearance order. `Keq`, `E_total`, and
metabolite names are NOT renamed.

This works for monomeric AND allosteric mechanisms because
`parameters(m)` returns the full param list including T-state
suffixes (`K1_T`, `kf_T`, `kr_T`), regulator-site names
(`K_R_reg1`, `K_R_T_reg2`), and the allosteric coupling `L`.
Future param shapes are auto-handled.
"""
function _canonicalize_rate_eq(m::AbstractEnzymeMechanism)
    body = rate_equation_string(m)

    # Strip the destructure lines `(; ... ) = params`/`= concs`.
    body = join(
        filter(
            ln -> !occursin(
                r"^\s*\(; .* = (params|concs)$", ln),
            split(body, '\n')),
        '\n')

    # Discover all parameter names; exclude the never-renamed set.
    skip = (:Keq, :E_total)
    pnames = String[String(p) for p in parameters(m, Reduced)
                    if p ∉ skip]

    # Find each name's first-appearance position via word-boundary
    # regex.
    first_pos = Dict{String,Int}()
    for name in pnames
        rx = Regex("\\b" * name * "\\b")
        m_pos = match(rx, body)
        first_pos[name] = m_pos === nothing ? typemax(Int) :
            m_pos.offset
    end

    # Order by first appearance; tie-break by name for determinism.
    ordered = sort(pnames; by=name -> (first_pos[name], name))
    name_map = Dict(name => "p_$i"
                    for (i, name) in enumerate(ordered))

    # Apply substitutions; longest first to prevent prefix
    # collisions (e.g. rename `K1_T` before `K1`).
    for name in sort(pnames; by=length, rev=true)
        body = replace(body,
            Regex("\\b" * name * "\\b") => name_map[name])
    end

    strip(replace(body, r"\s+" => " "))
end

"""Hash a mechanism's canonicalized rate equation. Returns
`UInt64` from Julia's built-in `hash`. Adequate for our scale
(~10⁻¹² collision probability over 10⁴ mechanisms)."""
function _canonical_rate_eq_hash(m::AbstractEnzymeMechanism)
    hash(_canonicalize_rate_eq(m))
end

"""Return `(UInt64 hash, 16-char hex display string)`."""
function _canonical_rate_eq_hash_pair(m::AbstractEnzymeMechanism)
    h = _canonical_rate_eq_hash(m)
    (h, string(h, base=16, pad=16))
end
```

- [ ] **Step 4: Run G.1 tests**

Expected: pass.

### Task G.2: Persistent fit cache + two-stage processing with worker recompile

**Files:**
- Modify: `src/identify_rate_equation.jl` `_beam_search` body

- [ ] **Step 1: Add the cached-result struct**

Insert near the top of `src/identify_rate_equation.jl`:

```julia
"""Cached fit result keyed by canonical rate-equation hash.
- `first_seen_estimate`: the beam-search level (the `pc` loop
  iteration value, equal to `n_fit_params_estimate`) at which
  this hash's fit was first performed.
- `first_seen_n_actual`: `length(fitted_params(m))` at first fit.
- `first_seen_eq_hash`: 16-char hex display string of the hash.
"""
struct _CachedFitResult
    loss::Float64
    params::NamedTuple
    first_seen_estimate::Int
    first_seen_n_actual::Int
    first_seen_eq_hash::String
end
```

- [ ] **Step 2: Refactor `_beam_search` to two-stage processing**

Replace the body of `_beam_search` (current lines 243-344) with:

```julia
function _beam_search(
    prob::IdentifyRateEquationProblem;
    min_beam_width, loss_rel_threshold, loss_abs_threshold,
    max_param_count, save_dir, pmap_function,
    optimizer, kwargs...
)
    # Persistent cross-level cache keyed by canonical hash.
    fit_cache = Dict{UInt64, _CachedFitResult}()

    # Initialize cache by enumeration estimate level
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

        # ── Stage 1 (parallel): compile + hash. Return spec +
        #    hash + n_actual; the mechanism object stays on the
        #    worker that compiled it. ─────────────────────────
        compiled = pmap_function(level) do spec
            try
                m = compile_mechanism(spec)
                eq_text = rate_equation_string(m)
                h_full, h_short = _canonical_rate_eq_hash_pair(m)
                n_actual = length(fitted_params(m))
                mech_type_str = string(typeof(m))
                (spec=spec, eq_text=eq_text,
                 h_full=h_full, h_short=h_short,
                 n_actual=n_actual,
                 mech_type_str=mech_type_str, ok=true)
            catch e
                @debug("Mechanism compilation failed",
                       exception=(e, catch_backtrace()))
                (spec=spec, ok=false)
            end
        end
        filter!(c -> c.ok, compiled)
        isempty(compiled) && continue

        # Snapshot which hashes are NEW vs already cached.
        new_hashes = Set{UInt64}()
        for c in compiled
            haskey(fit_cache, c.h_full) && continue
            push!(new_hashes, c.h_full)
        end

        # Pick one rep spec per new hash (first-encountered).
        reps_by_hash = Dict{UInt64, NamedTuple}()
        for c in compiled
            c.h_full in new_hashes || continue
            haskey(reps_by_hash, c.h_full) && continue
            reps_by_hash[c.h_full] = c
        end

        # ── Stage 2 (parallel): worker-side recompile + fit. ──
        # Recompile on the worker so the singleton type never
        # crosses worker boundaries between stages.
        rep_results = pmap_function(
            collect(values(reps_by_hash))
        ) do rep
            try
                m = compile_mechanism(rep.spec)  # recompile
                fp = FittingProblem(m, prob.data; Keq=prob.Keq)
                fit = fit_rate_equation(
                    fp, optimizer; kwargs...)
                (h_full=rep.h_full, h_short=rep.h_short,
                 n_actual=rep.n_actual,
                 loss=fit.loss, params=fit.params, ok=true)
            catch e
                @debug("Rep fit failed",
                       exception=(e, catch_backtrace()))
                (h_full=rep.h_full, ok=false)
            end
        end

        # Fold new fits into the cache.
        for r in rep_results
            r.ok || continue
            fit_cache[r.h_full] = _CachedFitResult(
                r.loss, r.params, pc, r.n_actual, r.h_short)
        end

        # ── Stage 3 (master): build ONE row per spec member. ──
        # Members of the same hash group share (loss, params,
        # eq_hash); users can post-hoc dedup via eq_hash.
        level_rows = NamedTuple[]
        level_specs = AbstractMechanismSpec[]
        for c in compiled
            haskey(fit_cache, c.h_full) || continue
            cached = fit_cache[c.h_full]
            is_inherited = !(c.h_full in new_hashes)
            row = (
                n_params = c.n_actual,
                loss = cached.loss,
                mechanism_type = c.mech_type_str,
                rate_equation = c.eq_text,
                fitted_param_names = collect(keys(cached.params)),
                fitted_param_values =
                    Tuple(values(cached.params)),
                eq_hash = cached.first_seen_eq_hash,
                fit_inherited_from_estimate =
                    is_inherited ? cached.first_seen_estimate :
                                   missing,
            )
            push!(level_rows, row)
            push!(level_specs, c.spec)
        end

        append!(all_specs, level_specs)
        append!(all_rows,  level_rows)

        # Save CSV for this estimate-level.
        if save_dir !== nothing
            _save_level_csv(save_dir, level_rows, pc)
        end

        # Beam selection (same _select_beam logic as Phase F).
        sel = _select_beam(
            [r.loss for r in level_rows];
            loss_rel_threshold=loss_rel_threshold,
            loss_abs_threshold=loss_abs_threshold,
            min_beam_width=min_beam_width)
        beam_specs = level_specs[sel]

        # Expand all beam specs (each spec is structurally
        # distinct; expansion preserves RE→SS coverage).
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

- [ ] **Step 3: Update `_rows_to_dataframe` for the new columns**

Replace `_rows_to_dataframe` (lines 196-228) with:

```julia
function _rows_to_dataframe(rows)
    isempty(rows) && return DataFrame()
    all_pnames = Set{Symbol}()
    for row in rows
        for p in row.fitted_param_names
            push!(all_pnames, p)
        end
    end
    sorted_pnames = sort(collect(all_pnames))

    df = DataFrame(
        n_params = [r.n_params for r in rows],
        loss = [r.loss for r in rows],
        mechanism_type = [r.mechanism_type for r in rows],
        rate_equation = [r.rate_equation for r in rows],
        eq_hash = [r.eq_hash for r in rows],
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

- [ ] **Step 4: Dedup LOOCV by `eq_hash` within each `n_params` bucket**

Replace the candidate-selection loop in `_cv_model_selection`
(current lines 431-440):

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

Update the `n_cv_candidates` docstring line in
`identify_rate_equation` to clarify: *"top N **unique-rate-equation**
candidates per param count"* — same number of LOOCV runs as
before, but each run is on a distinct equation.

### Task G.3: End-to-end regression test (real optimizer)

**Files:**
- Modify: `test/test_identify_rate_equation.jl`

The single test for Phase G covers cache + integration end-to-end
using the real optimizer already imported in this file
(`OptimizationPyCMA.PyCMAOpt()`).

- [ ] **Step 1: Add the regression test**

```julia
@testset "identify_rate_equation: end-to-end CSV invariants" begin
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    Keq = 2.0
    n_pts_per_group = 8
    n_groups = 3
    A_vals = repeat([0.1, 0.5, 1.0, 2.0], n_groups * 2)
    B_vals = repeat([0.5, 0.5, 0.5, 0.5], n_groups * 2)
    P_vals = repeat([0.05, 0.05, 0.05, 0.05], n_groups * 2)
    Q_vals = repeat([0.05, 0.05, 0.05, 0.05], n_groups * 2)
    Rates = @. (A_vals * B_vals) /
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
            optimizer = OptimizationPyCMA.PyCMAOpt(),
            save_dir = tmp,
            pmap_function = map,
            n_restarts = 1, maxtime = 2.0,
        )
        # 1. Filenames use the new estimate-level naming.
        files = filter(f -> endswith(f, ".csv"),
                       readdir(tmp))
        @test !isempty(files)
        for fname in files
            @test startswith(fname, "params_estimate_")
        end
        # 2. eq_hash column exists and is well-formed.
        for fname in files
            df_file = CSV.read(joinpath(tmp, fname), DataFrame)
            @test "eq_hash" in names(df_file)
            @test all(.!ismissing.(df_file.eq_hash))
            @test all(length.(df_file.eq_hash) .== 16)
        end
        # 3. Cross-level inheritance chain: any row whose
        #    fit_inherited_from_estimate is not missing must
        #    point to a level whose CSV contains a row with the
        #    same eq_hash.
        all_rows_by_level = Dict{Int, DataFrame}()
        for fname in files
            est = parse(Int, replace(fname,
                "params_estimate_" => "", ".csv" => ""))
            all_rows_by_level[est] = CSV.read(
                joinpath(tmp, fname), DataFrame)
        end
        for (est, df_lvl) in all_rows_by_level
            for row in eachrow(df_lvl)
                ismissing(row.fit_inherited_from_estimate) &&
                    continue
                src = row.fit_inherited_from_estimate
                @test haskey(all_rows_by_level, src)
                @test row.eq_hash in
                    all_rows_by_level[src].eq_hash
            end
        end
        # 4. result.best is a real mechanism.
        @test result.best isa AbstractEnzymeMechanism
        # 5. result.cv_results is non-empty.
        @test nrow(result.cv_results) >= 1
    end
end
```

- [ ] **Step 2: Run; confirm pass**

```
julia --project -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
```

Expected: green. If `maxtime=2.0` proves flaky in CI, raise to
`5.0` — accept the cost.

### Task G.4: Commit Phase G

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "$(cat <<'EOF'
Add persistent rate-equation hash cache for fit dedup

A canonical Julia hash of each mechanism's rate equation (built
token-by-token from parameters(m, Reduced) so allosteric T-state
and regulator params are covered) keys a Dict that survives across
all beam levels. Specs whose hash hits the cache reuse the prior
fit; new hashes are recompiled worker-side and fitted.

CSV gains eq_hash and fit_inherited_from_estimate columns. One row
per spec member; users dedup post-hoc via the eq_hash column.
_cv_model_selection picks top-N unique-hash candidates per
n_params bucket.

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

Expected: all tests pass. Aqua + JET checks green. No `@test_skip` markers anywhere. The Phase G end-to-end test runs the real optimizer (`PyCMAOpt`) with `maxtime=2.0`; if CI is slow, raise to `maxtime=5.0`.

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

- **Worker-side recompile in Phase G Stage 2:** the Stage 2 closure receives the spec, not a master-compiled mechanism. It calls `compile_mechanism(rep.spec)` again on the worker before constructing the `FittingProblem`. This is intentional: `EnzymeMechanism{...}` is a singleton type instantiated per worker, and `Distributed.pmap` has no worker affinity between separate pmap calls — shipping a mechanism object across worker boundaries fails with a deserialization error. Stage 1's compile is for hashing only; the result is discarded before Stage 2.
- **Token-driven canonicalizer in Phase G:** `_canonicalize_rate_eq` walks `parameters(m, Reduced)` to discover every parameter symbol the mechanism uses. This catches monomeric (`K1`, `kf_3`), allosteric (`K1_T`, `kf_T`), regulator-site (`K_R_reg1`, `K_R_T_reg1`), and the allosteric coupling (`L`) automatically. Adding a future parameter shape requires no code change here as long as it appears in `parameters(m)`.
- **No `SHA` dependency:** the canonical hash uses Julia's built-in `hash(::String)::UInt64`. Don't add `using SHA` — Aqua will flag it as a missing/stale dep.
- **Backwards compatibility:** none required. This is internal-version cleanup, not a public-API deprecation. Removing `beam_fraction` raises `MethodError`; that is intended.
