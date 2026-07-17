# `:OnlyA` Haldane Validator Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reject allosteric mechanisms whose `:OnlyA` tags make a catalytic cycle's Haldane relation unsatisfiable, and stop the enumerator from generating them.

**Architecture:** A pure validator (`_onlya_haldane_violation`) reuses the existing cycle machinery: build a check graph containing every catalytic group except `:OnlyA` *chemical* groups, assemble its Haldane cycles with `_assemble_constraints`, and inspect the sign pattern of the `:OnlyA`-binding K columns in each cycle row. All-same-sign means the cycle's Haldane cannot be satisfied. The two moves that write `:OnlyA` (`_expand_to_allosteric`, `_expand_promote_catalytic_to_onlya`) close over the promotions the validator forces, so the enumerator only emits valid mechanisms. The validator is wired into the `AllostericMechanism` constructor last, after every in-tree fixture is corrected.

**Tech Stack:** Julia 1.12, EnzymeRates.jl. Tests via `Pkg.test()` / `TestEnv`.

**Spec:** `docs/superpowers/specs/2026-07-15-onlya-haldane-validator-design.md`. Read it before Task 1.

**Branch:** `onlya-haldane-validator` (already created, spec committed at `334c22c`).

## Global Constraints

- 92-character line length, 4-space indentation.
- All new files start with two `# ABOUTME: ` lines.
- `rate_equation` stays allocation-free and under 120 ns. This change touches
  construction and derivation only; never edit the `rate_equation` codegen path.
- All `Parameter → Symbol` rendering goes through the `name(p, m)` chokepoint. Do
  not write `Symbol("K…")`/`Symbol("k…")` literals — the AST-walker test at
  `test/test_types.jl:1577-1644` fails the build on any.
- **The validator must never call `_state_allo_mechanism`.** That function
  constructs an `AllostericMechanism` (`rate_eq_derivation.jl:1232`), so calling
  it from the constructor recurses forever. Build the check graph as a plain
  `Mechanism`.
- **The validator must never substitute a numeric `Inf` or compare infinities.**
  Work with the sign pattern of exponents only. Comparing infinities discards the
  K-system degree of freedom and rejects valid mechanisms.
- Run the full suite before the final commit: `julia --project -e 'using Pkg; Pkg.test()'`.
- Commit after each task. Never skip a pre-commit hook.

---

## Background the implementer needs

An `:OnlyA` binding asserts `K_I = ∞`. Both conformations share `Keq`, so writing
`K_I = K_A/ε` with `ε → 0⁺` and dividing the two Haldane relations gives:

```
k_I_f/k_I_r  =  (k_A_f/k_A_r) · ∏ε_p / ∏ε_s
```

A non-`:OnlyA` catalytic tag asserts that ratio is a finite nonzero constant,
which requires `∏ε_p/∏ε_s` to be `O(1)`. Each `ε` is independent, so writing
`ε_i = t^{c_i}` (`c_i > 0`, `t → 0⁺`) makes the monomial `t^{Σ a_i c_i}`. That is
`O(1)` exactly when the exponents `a_i` carry **both signs** — then the `c_i` can
be chosen to cancel. All-same-sign forces the monomial to `0` or `∞`, and the
only escape is `k_I = 0`, which only an `:OnlyA` catalytic tag expresses.

`_assemble_constraints` already yields these exponents. Its docstring
(`thermodynamic_constr_for_rate_eq_derivation.jl:315-333`) says: "`A` is the
constraint matrix (rows = independent Wegscheider/Haldane cycles, columns =
parameters in `all_params` order)" and "Binding K's are Kd in the polynomial while
cycle products use 1/Kd, so binding-K column entries carry a sign flip on top of
the cycle incidence."

So the algorithm is:

1. Build a check `Mechanism` from every catalytic group **except** `:OnlyA`
   chemical (iso) groups. Dropping those is what encodes the `k_I = 0` escape: a
   cycle through an `:OnlyA` chemical step does not exist, so it produces no row
   and cannot be violated.
2. `_assemble_constraints` it.
3. For each row, collect `sign(A[i, c])` over the columns `c` belonging to
   `:OnlyA` **binding** groups, ignoring zeros.
4. Empty (no `:OnlyA` binding in this cycle) → fine. Both signs → balanced
   K-system → fine. One sign → **violation**.

Bindings that complete no cycle — competitive inhibitors, dead ends, regulator
sites — never appear in a row, so they are excluded automatically. This is why the
check must not key on metabolite names: LDH's `Lactate`/`NAD` competitive
inhibitors share names with its products.

---

### Task 1: The Haldane validator

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl` (append after `_assemble_constraints`)
- Test: `test/test_types.jl` (new testset at end of file)

**Interfaces:**
- Consumes: `_assemble_constraints(mech, rename; step_params, all_params, is_i_state)`
  from `src/thermodynamic_constr_for_rate_eq_derivation.jl:334`; `Mechanism`,
  `Step`, `is_binding`, `is_iso`, `_flat_steps`, `_step_parameters`, `name`.
- Produces: `_onlya_haldane_violation(rxn::EnzymeReaction,
  cat_steps::Vector{Vector{Step}}, cat_allo_states::Vector{Symbol})
  → Union{Nothing, String}`. Returns `nothing` when valid, else a message.
  Tasks 2, 3, and 7 call it.

- [ ] **Step 1: Write the failing test**

Append to `test/test_types.jl`:

```julia
@testset "OnlyA Haldane validator" begin
    ER = EnzymeRates
    # Uni-uni S -> P. Tags: (S binding, chemical step, P binding).
    function uni(s_tag, cat_tag, p_tag)
        m = @allosteric_mechanism begin
            substrates: S
            products:   P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
            end
        end
        am = ER.AllostericMechanism(m)
        tags = copy(ER.cat_allo_states(am))
        # find each group by its representative step
        for (g, grp) in enumerate(ER.steps(am))
            bm = ER.bound_metabolite(grp[1])
            tags[g] = bm === nothing ? cat_tag :
                      ER.name(bm) === :S ? s_tag : p_tag
        end
        (ER.reaction(am), ER.steps(am), tags)
    end

    # no :OnlyA anywhere -> valid
    @test ER._onlya_haldane_violation(uni(:EqualAI, :EqualAI, :EqualAI)...) === nothing
    # :OnlyA on the substrate only, catalysis :EqualAI -> VIOLATION
    @test ER._onlya_haldane_violation(uni(:OnlyA, :EqualAI, :EqualAI)...) isa String
    # :OnlyA on the product only, catalysis :EqualAI -> VIOLATION
    @test ER._onlya_haldane_violation(uni(:EqualAI, :EqualAI, :OnlyA)...) isa String
    # :OnlyA on the substrate, catalysis :OnlyA -> the k_I = 0 escape -> valid
    @test ER._onlya_haldane_violation(uni(:OnlyA, :OnlyA, :EqualAI)...) === nothing
    # balanced: :OnlyA on both sides, catalysis :EqualAI -> valid
    @test ER._onlya_haldane_violation(uni(:OnlyA, :EqualAI, :OnlyA)...) === nothing
    # V-system: :OnlyA chemical step only -> valid
    @test ER._onlya_haldane_violation(uni(:EqualAI, :OnlyA, :EqualAI)...) === nothing
    # :NonequalAI catalysis is also a finite-nonzero assertion -> same verdict
    @test ER._onlya_haldane_violation(uni(:OnlyA, :NonequalAI, :EqualAI)...) isa String
end
```

- [ ] **Step 2: Run it and confirm it fails**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_types.jl")' 2>&1 | tail -20
```

Expected: `UndefVarError: _onlya_haldane_violation`.

- [ ] **Step 3: Implement**

Append to `src/thermodynamic_constr_for_rate_eq_derivation.jl`:

```julia
"""
    _onlya_haldane_violation(rxn, cat_steps, cat_allo_states)
        → Union{Nothing, String}

Return `nothing` when every catalytic thermodynamic cycle's Haldane relation
stays satisfiable under the `K_I = K_A/ε`, `ε → 0⁺` limit an `:OnlyA` binding
asserts; otherwise return a message naming the offending cycle.

A cycle's Haldane carries `∏ε_p/∏ε_s` on the inactive side. The `ε` are
independent, so that monomial can be held finite exactly when its exponents
carry both signs. All-same-sign drives it to `0` or `∞`, and only `k_I = 0` —
an `:OnlyA` catalytic tag — absorbs that.

The check graph drops `:OnlyA` chemical groups, so a cycle running through one
never appears and never reports a violation: that is the `k_I = 0` escape.
Bindings completing no cycle (competitive inhibitors, dead ends, regulator
sites) never enter a row and take no part.

Builds a plain `Mechanism`; it must not call `_state_allo_mechanism`, which
would construct an `AllostericMechanism` and recurse.
"""
function _onlya_haldane_violation(rxn::EnzymeReaction,
                                  cat_steps::Vector{Vector{Step}},
                                  cat_allo_states::Vector{Symbol})
    keep = [g for g in eachindex(cat_steps)
            if !(cat_allo_states[g] === :OnlyA && is_iso(cat_steps[g][1]))]
    isempty(keep) && return nothing
    onlyA_steps = Set{Step}()
    for g in eachindex(cat_steps)
        cat_allo_states[g] === :OnlyA && is_binding(cat_steps[g][1]) &&
            union!(onlyA_steps, cat_steps[g])
    end
    isempty(onlyA_steps) && return nothing
    cm = Mechanism(rxn, [copy(cat_steps[g]) for g in keep])
    sp = _step_parameters(cm)
    A, _, columns, _ = _assemble_constraints(cm, Dict{Symbol, Symbol}();
                                             step_params = sp)
    sym_col = Dict(c => i for (i, c) in enumerate(columns))
    onlyA_cols = Set{Int}()
    for (j, (s, _)) in enumerate(_flat_steps(cm))
        s in onlyA_steps || continue
        sym = name(sp[j][1], cm)
        haskey(sym_col, sym) && push!(onlyA_cols, sym_col[sym])
    end
    isempty(onlyA_cols) && return nothing
    for i in axes(A, 1)
        signs = Set{Int}()
        for c in onlyA_cols
            A[i, c] == 0 || push!(signs, A[i, c] > 0 ? 1 : -1)
        end
        isempty(signs) && continue
        length(signs) == 1 || continue
        offenders = sort!([string(columns[c]) for c in onlyA_cols if A[i, c] != 0])
        return "an :OnlyA binding ($(join(offenders, ", "))) leaves a " *
               "catalytic cycle's Haldane relation unsatisfiable: the inactive " *
               "conformation cannot run that cycle at a finite nonzero rate. " *
               "Tag the cycle's chemical step :OnlyA, or tag an opposing " *
               "binding :OnlyA so the affinities diverge together."
    end
    nothing
end
```

- [ ] **Step 4: Run the test and confirm it passes**

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_types.jl")' 2>&1 | tail -20
```

Expected: all 7 assertions pass.

- [ ] **Step 5: Validate against the real HPC corpus**

The session cached 200 LDH and 359 PFKP real mechanism types. Run:

```bash
julia --project -e '
using EnzymeRates; const ER = EnzymeRates
S = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/fb96a560-d2b8-4694-a474-ff2f47968c52/scratchpad"
for f in ("ldh_err_types.txt", "ldh_ok_types.txt", "pfkp_err_types.txt", "pfkp_ok_types.txt")
    lines = readlines(joinpath(S, f)); bad = 0
    for l in lines
        am = ER.AllostericMechanism(Core.eval(ER, Meta.parse(String(l)))())
        ER._onlya_haldane_violation(ER.reaction(am), ER.steps(am),
                                    ER.cat_allo_states(am)) === nothing || (bad += 1)
    end
    println("$f: $bad / $(length(lines)) flagged invalid")
end'
```

Expected, matching the spec's measured census:
`ldh_err_types.txt: 71 / 71`, `ldh_ok_types.txt: 129 / 129`,
`pfkp_err_types.txt: 104 / 118`, `pfkp_ok_types.txt: 172 / 241`.

If the numbers differ, the validator is wrong — do not proceed. The 14 PFKP
errors it does not flag are the documented out-of-scope V-system defect.

- [ ] **Step 6: Commit**

```bash
git add src/thermodynamic_constr_for_rate_eq_derivation.jl test/test_types.jl
git commit -m "Add :OnlyA Haldane validator (not yet wired in)"
```

---

### Task 2: Close `_expand_to_allosteric` over forced promotions

**Files:**
- Modify: `src/mechanism_enumeration.jl:1772-1800` (`_expand_to_allosteric`)
- Test: `test/test_mechanism_enumeration.jl` (new testset)

**Interfaces:**
- Consumes: `_onlya_haldane_violation` from Task 1.
- Produces: `_valid_onlya_completions(rxn, cat_steps, tags) → Vector{Vector{Symbol}}`
  — every minimal tag vector that promotes additional `:EqualAI` groups to
  `:OnlyA` until the Haldane holds. Task 3 calls it too.

- [ ] **Step 1: Write the failing test**

Append to `test/test_mechanism_enumeration.jl`:

```julia
@testset "to_allosteric emits only Haldane-valid mechanisms" begin
    ER = EnzymeRates
    rxn = EnzymeReaction(
        substrates = [:S => (C = 6,)], products = [:P => (C = 6,)],
        Keq = 1.0, allowed_catalytic_multiplicities = (2,))
    m = @enzyme_mechanism begin
        substrates: S
        products:   P
        catalytic_steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    kids = ER._expand_to_allosteric(ER.Mechanism(m), rxn)
    @test !isempty(kids)
    for am in kids
        @test ER._onlya_haldane_violation(
            ER.reaction(am), ER.steps(am), ER.cat_allo_states(am)) === nothing
    end
    # every child still carries at least one :OnlyA — the move's purpose
    @test all(am -> :OnlyA in ER.cat_allo_states(am), kids)
end
```

If the `EnzymeReaction`/`@enzyme_mechanism` construction above does not match this
repo's current API, copy the exact construction used by the existing
`biuni_seed()` at `test/test_mechanism_enumeration.jl:4058-4070` and adapt the
metabolite names — do not invent an API.

- [ ] **Step 2: Run it and confirm it fails**

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_mechanism_enumeration.jl")' 2>&1 | tail -25
```

Expected: FAIL — children with an `:OnlyA` binding and `:EqualAI` catalysis are
flagged.

- [ ] **Step 3: Implement the completion search**

Add to `src/mechanism_enumeration.jl`, above `_expand_to_allosteric`:

```julia
"""
    _valid_onlya_completions(rxn, cat_steps, tags) → Vector{Vector{Symbol}}

Every minimal completion of `tags` that satisfies the Haldane relation, found by
promoting additional `:EqualAI` groups to `:OnlyA`. Returns `[tags]` unchanged
when `tags` already holds. Searches promotion subsets by increasing size and
returns every valid vector at the first size that yields one, so both
completions of a one-sided `:OnlyA` — promote the chemical step, or promote an
opposing binding — are emitted. Returns an empty vector when no completion
exists.
"""
function _valid_onlya_completions(rxn::EnzymeReaction,
                                  cat_steps::Vector{Vector{Step}},
                                  tags::Vector{Symbol})
    _onlya_haldane_violation(rxn, cat_steps, tags) === nothing && return [copy(tags)]
    cand = [g for g in eachindex(tags) if tags[g] === :EqualAI]
    for k in 1:length(cand)
        found = Vector{Symbol}[]
        for combo in Combinatorics.combinations(cand, k)
            t = copy(tags)
            for g in combo
                t[g] = :OnlyA
            end
            _onlya_haldane_violation(rxn, cat_steps, t) === nothing && push!(found, t)
        end
        isempty(found) || return found
    end
    Vector{Symbol}[]
end
```

`Combinatorics` is not a dependency. Do not add one — write the subset loop by
hand instead:

```julia
    for k in 1:length(cand)
        found = Vector{Symbol}[]
        _each_subset(cand, k) do combo
            t = copy(tags)
            for g in combo
                t[g] = :OnlyA
            end
            _onlya_haldane_violation(rxn, cat_steps, t) === nothing && push!(found, t)
        end
        isempty(found) || return found
    end
```

with a small local helper:

```julia
# Call `f` on every size-`k` subset of `xs`.
function _each_subset(f, xs::Vector{Int}, k::Int)
    n = length(xs)
    idx = collect(1:k)
    k > n && return
    while true
        f([xs[i] for i in idx])
        i = k
        while i ≥ 1 && idx[i] == n - k + i
            i -= 1
        end
        i == 0 && return
        idx[i] += 1
        for j in (i + 1):k
            idx[j] = idx[j - 1] + 1
        end
    end
end
```

- [ ] **Step 4: Wire it into `_expand_to_allosteric`**

Replace the binding branch (`src/mechanism_enumeration.jl:1786-1795`). The
current code is:

```julia
            else
                push!(results, AllostericMechanism(
                    reaction(m), copy(steps(m)), new_tags, cn, RegulatorySite[]))
            end
```

Replace with:

```julia
            else
                for t in _valid_onlya_completions(rxn, steps(m), new_tags)
                    push!(results, AllostericMechanism(
                        reaction(m), copy(steps(m)), t, cn, RegulatorySite[]))
                end
            end
```

Leave the `is_iso` branch alone — a `:OnlyA` chemical group with all-`:EqualAI`
bindings is a V-system and already valid.

- [ ] **Step 5: Run the test and confirm it passes**

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_mechanism_enumeration.jl")' 2>&1 | tail -25
```

Expected: the new testset passes. Other testsets in this file may now fail on
child counts — that is Task 6's job. Record which fail; do not fix them here.

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Close _expand_to_allosteric over Haldane-forced promotions"
```

---

### Task 3: Close `_expand_promote_catalytic_to_onlya` and fix its docstring

**Files:**
- Modify: `src/mechanism_enumeration.jl:1969-1990`
- Test: `test/test_mechanism_enumeration.jl`

**Interfaces:**
- Consumes: `_valid_onlya_completions` from Task 2.

- [ ] **Step 1: Write the failing test**

```julia
@testset "promote_catalytic_to_onlya emits only Haldane-valid mechanisms" begin
    ER = EnzymeRates
    for am in ER._expand_to_allosteric(ER.Mechanism(biuni_seed_mech()), biuni_seed_rxn())
        for kid in ER._expand_promote_catalytic_to_onlya(am)
            @test ER._onlya_haldane_violation(
                ER.reaction(kid), ER.steps(kid), ER.cat_allo_states(kid)) === nothing
        end
    end
end
```

Use whatever seed helpers Task 6 leaves in place; if `biuni_seed()` still returns
a single mechanism, adapt to it. The assertion — every child is valid — is the
contract; the seed plumbing is incidental.

- [ ] **Step 2: Run it and confirm it fails**

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_mechanism_enumeration.jl")' 2>&1 | tail -25
```

- [ ] **Step 3: Implement**

Replace the body of `_expand_promote_catalytic_to_onlya(am::AllostericMechanism)`:

```julia
function _expand_promote_catalytic_to_onlya(am::AllostericMechanism)
    results = AllostericMechanism[]
    for g in 1:length(cat_allo_states(am))
        cat_allo_states(am)[g] == :EqualAI || continue
        new_states = copy(cat_allo_states(am))
        new_states[g] = :OnlyA
        for t in _valid_onlya_completions(reaction(am), steps(am), new_states)
            push!(results, _with_cat_allo_states(am, t))
        end
    end
    unique!(results)
    results
end
```

`unique!` is needed: two different seed promotions can close to the same tag
vector.

- [ ] **Step 4: Rewrite the docstring**

Replace the docstring at `src/mechanism_enumeration.jl:1969-1978`:

```julia
"""
    _expand_promote_catalytic_to_onlya(am::AllostericMechanism)
        → Vector{AllostericMechanism}

Catalytic-state move. For each catalytic kinetic group tagged `:EqualAI`, emit
the variants with that group set to `:OnlyA` — binding (K-type) and
iso/catalytic (V-type) groups alike — closed over whatever further promotions
the Haldane relation forces (`_valid_onlya_completions`). Promoting a binding on
one side of the reaction drives `∏ε_p/∏ε_s` to `0` or `∞`, so the inactive cycle
must either lose its chemical step or gain an opposing `:OnlyA` binding; both
completions are emitted.

The parameter count varies across children. Breaking the inactive cycle removes
Wegscheider constraints, so parameters the solver derived become parameters the
fit must find: a 6-parameter LDH parent yields 6-, 7-, and 8-parameter children.
The catalytic steps, multiplicity, regulatory sites, and every other tag pass
through unchanged.
"""
```

- [ ] **Step 5: Run the test and confirm it passes**

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_mechanism_enumeration.jl")' 2>&1 | tail -25
```

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Close promote_catalytic_to_onlya; correct its Δ0 docstring"
```

---

### Task 4: Rewrite the three `:OnlyA` ground-truth gates

**Files:**
- Modify: `test/allosteric_ground_truth.jl`

**Interfaces:** none produced.

All three `:OnlyA` gates encode the invalid combination: an `:OnlyA` binding with
`:EqualAI` catalysis, plus a phantom inactive species reachable only by inactive
catalysis. Read the whole file first.

- [ ] **Step 1: Confirm the gates currently fail the validator**

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/allosteric_ground_truth.jl")' 2>&1 | tail -20
```

Record the current pass counts before changing anything.

- [ ] **Step 2: Rewrite gate 1 (`uni_onlyA_flux`)**

In the mass-action network: delete the `:ES_I` species, both of its edges, and
its entry in `cat_edges`. In the paired `@allosteric_mechanism`, retag the
chemical step `:: OnlyA`. Update the `fitted_params` assertion — the catalysis
parameter renames `k_ES_to_EP` → `k_A_ES_to_EP`.

- [ ] **Step 3: Rewrite gates 2 and 3 the same way**

Gate 2 (`multi_onlyA_flux`) drops `:EAB_I`. Gate 3
(`metab_dfree_onlyA_flux`) drops `:ES_I` and `:ESB_I`. Retag each chemical step
`:: OnlyA` and update each `fitted_params` assertion.

- [ ] **Step 4: Add a kcat assertion to each gate**

`_kcat_forward` has no coverage on this class at all — that is how 1354 codegen
errors reached the cluster. Add to each gate, using that gate's own parameter
NamedTuple:

```julia
    @test isfinite(EnzymeRates._kcat_forward(m, prm))
```

- [ ] **Step 5: Run and confirm green**

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/allosteric_ground_truth.jl")' 2>&1 | tail -20
```

Expected: every gate passes, agreeing with its ground truth to `rtol = 1e-4` or
better. The session measured `≤ 4.5e-7` for the rewritten forms. If a gate
disagrees by more than that, stop — the rewrite is wrong, not the tolerance.

Expect gates 1 and 2 to now take the raw normalization branch rather than the
cross-weight branch. That is correct and expected: their `d_free_I` becomes `1`.
Gate 3 keeps the cross-weight branch covered.

- [ ] **Step 6: Commit**

```bash
git add test/allosteric_ground_truth.jl
git commit -m "Rewrite :OnlyA ground-truth gates to the thermodynamically valid form"
```

---

### Task 5: Fix the two invalid `MECHANISM_TEST_SPECS`

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Modify: `test/reference/allosteric_golden_reference.txt`

- [ ] **Step 1: Locate them**

```bash
grep -n "LDH i-state NonequalAI 6-group\|LDH i-state NonequalAI 5-group" \
  test/mechanism_definitions_for_test_enzyme_derivation.jl
```

Both carry an `:OnlyA` binding with a non-`:OnlyA` chemical step.

- [ ] **Step 2: Retag each chemical step `:: OnlyA`**

- [ ] **Step 3: Recompute the expected counts**

`expected_n_wegscheider_constraints` and `expected_n_independent_params` both
move — breaking the inactive cycle removes constraints, so the independent count
goes **up**. Measure, do not guess:

```bash
julia --project -e 'using TestEnv; TestEnv.activate();
include("test/mechanism_definitions_for_test_enzyme_derivation.jl");
# print fitted_params length for the two edited specs'
```

Write the measured values into the spec entries.

- [ ] **Step 4: Regenerate the goldens**

Find the golden-regeneration entry point:

```bash
grep -rn "allosteric_golden_reference" test/*.jl | head
```

Regenerate, then `git diff` the golden file. **Only the two edited specs'
`REDUCED_STRING`, `PARAMS_FULL`, and `PARAMS_REDUCED` blocks may change.** If
any other block moves, stop and investigate — the change has leaked.

- [ ] **Step 5: Run the derivation suite**

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_rate_eq_derivation.jl")' 2>&1 | tail -20
```

The performance gate (`test_rate_equation_performance`, 0 allocations and under
120 ns) must stay green. If it does not, stop and report — do not adjust the
threshold.

- [ ] **Step 6: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/reference/allosteric_golden_reference.txt
git commit -m "Retag the two invalid LDH i-state specs; regenerate goldens"
```

---

### Task 6: Fix `biuni_seed()` and the enumeration expectations

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

`biuni_seed()` (`:4058-4070`) is itself an instance of the invalid combination:
`E + A ⇌ E(A) :: OnlyA` with `E(A,B) <--> E(P) :: EqualAI`. It seeds four
promote-move testsets and will not construct once Task 7 lands.

- [ ] **Step 1: Retag `biuni_seed()`'s chemical step `:: OnlyA`**

- [ ] **Step 2: Run the file and collect every failure**

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_mechanism_enumeration.jl")' 2>&1 | tail -40
```

- [ ] **Step 3: Update each child-count assertion to the measured value**

Both moves now emit more children per parent. For each failing count, print the
actual children and their tags, confirm each is Haldane-valid and structurally
distinct, then write the measured count in. **Do not update a count without
looking at what the children are** — a wrong count that "looks plausible" is how
this class of bug got here.

- [ ] **Step 4: Drop or rescope the Δ0 test**

The Δ0 testset at `:4134-4141` asserts a property that is false in general — it
passes only because `biuni_seed()` happens to satisfy it. The real LDH type
promotes 6 parameters to 6, 7, 7, and 8. Delete the testset. Its docstring claim
was already corrected in Task 3.

The "every promotion changes the rate equation" testset at `:4143` may also fail:
in the balanced region every catalytic tag renders the same equation. If it does,
delete it and note why in the commit message — that degeneracy is documented in
the spec as known and unaddressed.

- [ ] **Step 5: Run and confirm green**

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_mechanism_enumeration.jl")' 2>&1 | tail -20
```

- [ ] **Step 6: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Fix biuni_seed and enumeration expectations for the closed moves"
```

---

### Task 7: Wire the validator into the constructor

**Files:**
- Modify: `src/types.jl:607-655` (the `AllostericMechanism` inner constructor)
- Test: `test/test_types.jl`

This is the flip. Every in-tree fixture must already be valid (Tasks 4-6).

- [ ] **Step 1: Write the failing test**

```julia
@testset "AllostericMechanism rejects an unsatisfiable Haldane" begin
    @test_throws ErrorException @allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: EqualAI
            E(P) ⇌ E + P      :: EqualAI
        end
    end
end
```

- [ ] **Step 2: Run it and confirm it fails**

Expected: no exception thrown.

- [ ] **Step 3: Implement**

In the inner constructor, immediately before `new(...)` — after canonicalization,
so the check sees canonical groups:

```julia
        violation = _onlya_haldane_violation(reaction, cat_steps, cat_allo_states)
        violation === nothing ||
            error("AllostericMechanism: $violation")
```

- [ ] **Step 4: Run the full suite**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40
```

This is the moment the change either holds or does not. Every failure is either a
fixture Tasks 4-6 missed or a validator bug. **Fix the fixture or the validator —
never weaken the check to make a test pass.**

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Reject thermodynamically impossible :OnlyA mechanisms at construction"
```

---

### Task 8: Documentation

**Files:**
- Modify: `docs/src/deriving/mwc_allostery.md`

- [ ] **Step 1: Find the rule's prose**

Lines 207-209 already say "its catalytic step should be `:OnlyA` too, since a
state that cannot bind the substrate cannot catalyze." Change *should* to *must*,
and say the constructor now enforces it.

- [ ] **Step 2: State the balanced alternative**

Add: `:OnlyA` on both a substrate-side and a product-side binding is the
thermodynamically consistent exclusive K-system — the affinities diverge together
and their ratio stays free — so it keeps `:EqualAI` catalysis legal.

- [ ] **Step 3: Do not touch the `:NonequalAI` collapse example**

Lines 155-179 stay as they are. Rejecting a collapsed `:NonequalAI` is a
follow-up spec, and that `@example` still runs correctly today.

- [ ] **Step 4: Verify every `@example` block still executes**

```bash
julia --project -e 'using Pkg; Pkg.build("EnzymeRates")' 2>&1 | tail -5
grep -c "@example" docs/src/deriving/mwc_allostery.md
```

Run each `@example` block's code manually. Any block constructing an `:OnlyA`
binding with non-`:OnlyA` catalysis now throws and must be corrected.

- [ ] **Step 5: Commit**

```bash
git add docs/src/deriving/mwc_allostery.md
git commit -m "Document the :OnlyA catalytic rule as enforced"
```

---

## Final verification

- [ ] Full suite green: `julia --project -e 'using Pkg; Pkg.test()'`
- [ ] Performance gate green (0 allocations, under 120 ns).
- [ ] The Task 1 Step 5 corpus check still reports 71/71, 129/129, 104/118, 172/241.
- [ ] `git diff main --stat` touches only: `src/types.jl`,
      `src/thermodynamic_constr_for_rate_eq_derivation.jl`,
      `src/mechanism_enumeration.jl`, the four test files, one golden, and
      `docs/src/deriving/mwc_allostery.md`.
