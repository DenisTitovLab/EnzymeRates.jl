# Enumeration Review Round 3 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One definition of chemistry in the enumeration engine: a chemistry step is an isomerization (`is_iso`), which is how the enumerator writes every mechanism; the `_is_chemistry` predicate goes away, the move tests write their fixtures the enumerator's way, and the assumption is guarded at the engine's entry and stated in the docs.

**Architecture:** The derivation keeps accepting a step that folds chemistry into a release (`E(A) <--> F + P`, Segel's ping-pong, and the residual form on the ping-pong docs page): `rate_equation` lifts macro mechanisms through the same `Mechanism` type the engine uses, so no constructor check is possible. The engine instead relies on the enumerator, which only ever emits chemistry as an isomerization to a product-bound form followed by a release step, and `expand_mechanisms` rejects a parent with a binding step that changes the covalent residual. A binding step may change conformation (Denis: a conformation change is not a chemical reaction and can accompany binding); only the residual is chemistry. Test fixtures for the moves follow the enumerator's form, which becomes rule 4 of the enumeration-test rules.

**Tech Stack:** Julia 1.12, EnzymeRates internals, `Test`, TestEnv for focused runs.

**Spec:** Denis's round-2 review discussion of 2026-09-16 (questions 1, 2, 4 and the conformation comment) on top of `docs/superpowers/plans/2026-09-15-enumeration-review-round-2.md`. Question 5 (combinatorics page) is out of scope for this branch.

## Global Constraints

- 92-character lines, 4-space indentation, match surrounding style; every new test inline with the mechanism macro; move tests assert exact child sets; helpers prefixed `_testhelper_`; new testsets go after `end # top-level testset` in `test/test_mechanism_enumeration.jl` (CLAUDE.md "Enumeration-engine tests").
- Never edit an expectation to match output. If a move disagrees with an expected child set below, print both sets (as `EnzymeRates.steps.(...)`) and stop; the plan's derivation is then wrong and Denis decides.
- No comment may describe what the code used to do or that anything changed.
- One Julia process at a time; focused runs in the foreground with a 20-minute timeout; no background monitors; the full suite runs detached with the heap hint (below).
- Commit after each task with the attribution trailer from the session's system reminder.

## Running tests

Focused run of one file (about 4 minutes for the enumeration file; the `UndefVarError: random_reduced_params` in `init division-freeness (bi_bi_pp)` is a known cross-file artifact):

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")'
```

Full suite, detached, then poll `kill -0 <pid>` in bounded loops of at most 9.5 minutes:

```bash
nohup julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])' > "$CLAUDE_JOB_DIR/tmp/full_suite.log" 2>&1 &
```

## The enumerator's form of the docs ping-pong

The ping-pong docs page (`docs/src/deriving/ping_pong.md`) writes chemistry folded into the release. The enumerator writes the same mechanism as an isomerization to a product-bound residual form, then a release. Rendered form names are `E`, `EA`, `EP_res_+A_-P`, `E_res_+A_-P`, `EB_res_+A_-P`, `EQ`.

```julia
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
        end
    end)
```

Six singleton groups: A-binding (RE), chemistry 1 (SS iso), P-release (RE), B-binding (RE), chemistry 2 (SS iso), Q-release (RE). The reaction with an inhibitor:

```julia
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        competitive_inhibitors: I
    end
```

---

### Task 1: Move fixtures in the enumerator's form

**Files:**
- Test: `test/test_mechanism_enumeration.jl` (four sites, all after `end # top-level testset`)

**Interfaces:**
- Consumes: `EnzymeRates._flux_carrying_groups`, `_context_bipartitions`, `_expand_re_to_ss`, `_expand_add_dead_end_regulator`, `_add_competitive_inhibitor` (unchanged).
- Produces: fixtures that Task 2 keeps green when `_is_chemistry` is deleted.

Locate each site by the quoted line, not by number.

- [ ] **Step 1: `_flux_carrying_groups` ping-pong assertion.** In the testset that ends with `@test all(EnzymeRates._flux_carrying_groups(pp))`, replace the comment and the fixture:

```julia
    # Ping-pong: every group lies on the single catalytic cycle, so every
    # group is flux-carrying.
    pp = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
        end
    end)
    @test all(EnzymeRates._flux_carrying_groups(pp))
```

- [ ] **Step 2: residual-context fixture in `@testset "_context_bipartitions"`.** The fixture bound to `res` keeps its three-step Q group (its first step, `E + Q ⇌ E(Q)`, is also the release step of the second half-reaction) and writes the chemistry the enumerator's way:

```julia
    res = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + Q ⇌ E(Q), E(; residual = A - P) + Q ⇌ E(Q; residual = A - P),
             E(B; residual = A - P) + Q ⇌ E(B, Q; residual = A - P))
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
        end
    end)
```

The assertions that follow (`length(rbps) == 2`, division by B then by residual) do not change: the Q group's context forms are the same three forms.

- [ ] **Step 3: flip test.** Replace the whole `@testset "_expand_re_to_ss: ping-pong with combined chemistry-release steps"` with:

```julia
    @testset "_expand_re_to_ss: ping-pong" begin
        # The rapid-equilibrium subgraph is two trees, {E, E(A), E(Q)} and
        # {E(P; res), E(; res), E(B; res)}, so every RE group is a bridge and
        # flipping any one of them alone raises the segment count; no pair is
        # minimal.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(P; residual = A - P)
                E(P; residual = A - P) ⇌ E(; residual = A - P) + P
                E(; residual = A - P) + B ⇌ E(B; residual = A - P)
                E(B; residual = A - P) <--> E(Q)
                E(Q) ⇌ E + Q
            end
        end)
        flipA = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(P; residual = A - P)
                E(P; residual = A - P) ⇌ E(; residual = A - P) + P
                E(; residual = A - P) + B ⇌ E(B; residual = A - P)
                E(B; residual = A - P) <--> E(Q)
                E(Q) ⇌ E + Q
            end
        end)
        flipP = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(P; residual = A - P)
                E(P; residual = A - P) <--> E(; residual = A - P) + P
                E(; residual = A - P) + B ⇌ E(B; residual = A - P)
                E(B; residual = A - P) <--> E(Q)
                E(Q) ⇌ E + Q
            end
        end)
        flipB = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(P; residual = A - P)
                E(P; residual = A - P) ⇌ E(; residual = A - P) + P
                E(; residual = A - P) + B <--> E(B; residual = A - P)
                E(B; residual = A - P) <--> E(Q)
                E(Q) ⇌ E + Q
            end
        end)
        flipQ = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(P; residual = A - P)
                E(P; residual = A - P) ⇌ E(; residual = A - P) + P
                E(; residual = A - P) + B ⇌ E(B; residual = A - P)
                E(B; residual = A - P) <--> E(Q)
                E(Q) <--> E + Q
            end
        end)
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 4
        @test Set(kids) == Set([flipA, flipP, flipB, flipQ])
    end
```

- [ ] **Step 4: regulator exact-children test.** Replace the whole `@testset "_expand_add_dead_end_regulator: inhibitor mirrors one half-reaction"`. The plain enumerator-form ping-pong gives only three children and no mirror, because the inhibitor targets the forms that bind a competing ligand and no two of them are adjacent. Two dead-end bindings make the second half-reaction's forms targets: `E(Q) + A ⇌ E(A, Q)` (A binds before Q has left) and `E(B; res) + P ⇌ E(B, P; res)` (P still bound when B binds). Derivation, with A → {E, E(Q)}, B → {E(; res)}, P → {E(; res), E(B; res)}, Q → {E} as the forms that bind each ligand, and a form already carrying a competing ligand never a target:

| competes with | targets after exclusion | mirrored groups |
|---|---|---|
| A, P | E, E(Q), E(; res), E(B; res) | B-binding, chemistry 2, Q-release |
| A, Q | E | none |
| A; P, Q | E, E(; res), E(B; res) | B-binding |
| B, P | E(; res) | none |
| B, Q (also B; P,Q and A,B; Q and A,B; P,Q) | E, E(; res) | none |
| A, B; P | E, E(Q), E(; res) | Q-release |

```julia
@testset "_expand_add_dead_end_regulator: inhibitor mirrors one half-reaction" begin
    # The inhibitor binds the forms that bind a competing ligand, never a form
    # already carrying one, so a competing substrate's binding step is never
    # mirrored and no inhibitor-bound branch completes the net reaction. With A
    # binding E(Q) and P binding E(B; res) as dead ends, the pattern "compete
    # with A and P" puts I on E, E(Q), E(; res) and E(B; res): the second
    # half-reaction (B binds, chemistry, Q leaves) runs with I bound, and A
    # cannot bind E(I). Six patterns give distinct targets.
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        competitive_inhibitors: I
    end
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
        end
    end)
    mech(block) = EnzymeRates.Mechanism(block)
    second_half = mech(@enzyme_mechanism begin      # I competes with A and P
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            (E(; residual = A - P) + B ⇌ E(B; residual = A - P),
             E(I::Inh; residual = A - P) + B ⇌ E(B, I::Inh; residual = A - P))
            (E(B; residual = A - P) <--> E(Q),
             E(B, I::Inh; residual = A - P) <--> E(I::Inh, Q))
            (E(Q) ⇌ E + Q, E(I::Inh, Q) ⇌ E(I::Inh) + Q)
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            (E + I ⇌ E(I::Inh), E(Q) + I ⇌ E(I::Inh, Q),
             E(; residual = A - P) + I ⇌ E(I::Inh; residual = A - P),
             E(B; residual = A - P) + I ⇌ E(B, I::Inh; residual = A - P))
        end
    end)
    only_E = mech(@enzyme_mechanism begin           # I competes with A and Q
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            E + I ⇌ E(I::Inh)
        end
    end)
    B_binding = mech(@enzyme_mechanism begin        # I competes with A, P and Q
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            (E(; residual = A - P) + B ⇌ E(B; residual = A - P),
             E(I::Inh; residual = A - P) + B ⇌ E(B, I::Inh; residual = A - P))
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            (E + I ⇌ E(I::Inh),
             E(; residual = A - P) + I ⇌ E(I::Inh; residual = A - P),
             E(B; residual = A - P) + I ⇌ E(B, I::Inh; residual = A - P))
        end
    end)
    only_Eres = mech(@enzyme_mechanism begin        # I competes with B and P
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            E(; residual = A - P) + I ⇌ E(I::Inh; residual = A - P)
        end
    end)
    E_and_Eres = mech(@enzyme_mechanism begin       # I competes with B and Q
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            (E + I ⇌ E(I::Inh), E(; residual = A - P) + I ⇌ E(I::Inh; residual = A - P))
        end
    end)
    Q_release = mech(@enzyme_mechanism begin        # I competes with A, B and P
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            (E(Q) ⇌ E + Q, E(I::Inh, Q) ⇌ E(I::Inh) + Q)
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            (E + I ⇌ E(I::Inh), E(Q) + I ⇌ E(I::Inh, Q),
             E(; residual = A - P) + I ⇌ E(I::Inh; residual = A - P))
        end
    end)
    expected = [second_half, only_E, B_binding, only_Eres, E_and_Eres, Q_release]
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(kids) == 6
    @test Set(EnzymeRates.steps.(kids)) == Set(EnzymeRates.steps.(expected))
    @test all(c -> EnzymeRates.reaction(c) ==
                   EnzymeRates._add_competitive_inhibitor(rxn, :I), kids)
end
```

Group order inside each expected mechanism does not matter (the constructor canonicalizes); the mirror steps must sit in the same group as their counterpart, as written. If the macro rejects the bound-ligand spelling `E(I::Inh, Q)`, use the order the existing fixtures in this file use for an inhibitor next to another ligand (`E(B, I::Inh; residual = A - P)` is one) and report it.

- [ ] **Step 5: Run the file.** Focused run of `test/test_mechanism_enumeration.jl`. Expected: the four edited testsets pass (`_flux_carrying_groups`, `_context_bipartitions` 32/32, the flip testset 2/2, the regulator testset 3/3). Any child-set mismatch: print and stop.

- [ ] **Step 6: Commit.**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Write ping-pong move fixtures the enumerator's way"
```

---

### Task 2: Chemistry is the isomerization step

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_is_chemistry` and its docstring; `_flux_carrying_groups`)
- Test: `test/test_mechanism_enumeration.jl` (`@testset "_expand_add_dead_end_regulator: inhibitor never runs the net reaction"`)

**Interfaces:**
- Consumes: `is_iso(s::Step)` from `src/types.jl`.
- Produces: `_flux_carrying_groups(m)` keyed on `is_iso`; no `_is_chemistry` anywhere.

- [ ] **Step 1: Delete `_is_chemistry`** (the docstring block starting `    _is_chemistry(s::Step) -> Bool` through the function's `end`).

- [ ] **Step 2: `_flux_carrying_groups`.** In its docstring replace `containing a chemistry step (one `_is_chemistry` recognizes), i.e.` with `containing a chemistry step (an isomerization), i.e.`. In its body replace `push!(edge_group, g); push!(edge_is_chemistry, _is_chemistry(s))` with `push!(edge_group, g); push!(edge_is_chemistry, is_iso(s))`.

- [ ] **Step 3: The aggregate test.** In the testset named above replace `any(EnzymeRates._is_chemistry, branch) && (n_half += 1)` with `any(EnzymeRates.is_iso, branch) && (n_half += 1)`.

- [ ] **Step 4: Confirm nothing else refers to the predicate.**

```bash
grep -rn "_is_chemistry" src test docs/src .claude
```

Expected: no output.

- [ ] **Step 5: Run the file.** Focused run of `test/test_mechanism_enumeration.jl`; all testsets outside the giant block pass as in Task 1, plus `inhibitor never runs the net reaction` 352/352 (the `n_half > 0` assertion still holds because the seeds' chemistry steps are isomerizations).

- [ ] **Step 6: Commit.**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Take the isomerization step as the chemistry step in the flip move"
```

---

### Task 3: Guard the engine's assumption at `expand_mechanisms`

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`expand_mechanisms`, plus one helper next to `_assert_atom_conserving`)
- Test: `test/test_mechanism_enumeration.jl` (new testset after `@testset "expand_mechanisms: ter-ter random-order seed within budget"`)

**Interfaces:**
- Consumes: `is_binding`, `residual`, `from_species`, `to_species`, `name` (all in `src/types.jl`); `steps(m)` for both mechanism types.
- Produces: `_assert_chemistry_is_iso(m)`, called on every parent in `expand_mechanisms`.

Why here and not the `Mechanism` constructor: `rate_equation` lifts every macro mechanism through `Mechanism(em)` (`src/rate_eq_derivation.jl`, the `mech = Mechanism(em)` line), and the Segel ping-pong in `test/mechanism_definitions_for_test_enzyme_derivation.jl` and the docs page both fold chemistry into a release step on purpose. The derivation supports that form; only the moves assume the enumerator's form.

- [ ] **Step 1: Write the failing test.**

```julia
@testset "expand_mechanisms rejects chemistry folded into a release step" begin
    # The moves take the isomerization step as the chemistry step, which is
    # how the enumerator writes every mechanism. A mechanism written for the
    # derivation with chemistry folded into a release (the ping-pong docs
    # page) is not a valid parent.
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    folded = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end)
    err = try
        EnzymeRates.expand_mechanisms([folded], rxn); nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("changes the covalent residual", sprint(showerror, err))
    # A binding step may change conformation; only the residual is chemistry.
    conf = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ Estar(A)
            Estar(A) + B ⇌ Estar(A, B)
            Estar(A, B) <--> Estar(P, Q)
            Estar(P, Q) ⇌ Estar(P) + Q
            Estar(P) ⇌ E + P
        end
    end)
    @test EnzymeRates.expand_mechanisms([conf], rxn) isa Vector
end
```

- [ ] **Step 2: Run to verify it fails.** Focused run; expected: the first `@test err isa ErrorException` fails (`expand_mechanisms` returns children today) and the second assertion errors on `showerror(nothing)`.

- [ ] **Step 3: Implement.** Next to `_assert_atom_conserving` in `src/mechanism_enumeration.jl` add:

```julia
"""
    _assert_chemistry_is_iso(m)

The moves take the isomerization step as the chemistry step, which is how the
enumerator writes every mechanism: chemistry isomerizes to a product-bound form
and the release is its own step. A binding step may change the enzyme's
conformation but never its covalent residual; a mechanism written for the
derivation with chemistry folded into a release step is not a valid parent.
"""
function _assert_chemistry_is_iso(m::Union{Mechanism, AllostericMechanism})
    for group in steps(m), s in group
        is_binding(s) && residual(from_species(s)) != residual(to_species(s)) &&
            error("binding step $(name(from_species(s))) ⇌ " *
                  "$(name(to_species(s))) changes the covalent residual; the " *
                  "moves need the chemistry as an isomerization and the " *
                  "release as its own step")
    end
end
```

and in `expand_mechanisms` change the parent loop to:

```julia
    for m in mechs
        _assert_chemistry_is_iso(m)
        _add_expansions_mech!(result, m, rxn)
    end
```

- [ ] **Step 4: Run to verify it passes.** Focused run; the new testset 3/3; the budget testset and every other testset unchanged.

- [ ] **Step 5: Commit.**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Reject a parent whose chemistry is folded into a release step"
```

---

### Task 4: Split children from a conformation context and from a residual context

**Files:**
- Test: `test/test_mechanism_enumeration.jl` (two testsets after `@testset "_expand_split_kinetic_group (context bipartitions)"`)

**Interfaces:**
- Consumes: `EnzymeRates._expand_split_kinetic_group(m)` (exact children), `EnzymeRates._context_bipartitions`.

Round 2 tested the two context families only through `_context_bipartitions`; the 102-children pin did not move because no seed gains a child from them. These fixtures do gain: a dead-end binding constant is not tied by any cycle, so giving it its own kinetic group adds an independent parameter.

- [ ] **Step 1: Conformation context.**

```julia
@testset "_expand_split_kinetic_group: by conformation" begin
    # S binds both conformations in one kinetic group; Estar(S) is a dead end
    # whose dissociation constant no cycle ties, so dividing the group by
    # conformation frees one parameter. The other groups are singletons.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E + S ⇌ E(S), Estar + S ⇌ Estar(S))
            E(S) <--> Estar(P)
            Estar(P) ⇌ Estar + P
            Estar <--> E
        end
    end)
    by_conformation = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            Estar + S ⇌ Estar(S)
            E(S) <--> Estar(P)
            Estar(P) ⇌ Estar + P
            Estar <--> E
        end
    end)
    kids = EnzymeRates._expand_split_kinetic_group(m)
    @test length(kids) == 1
    @test Set(kids) == Set([by_conformation])
end
```

- [ ] **Step 2: Residual context.** Uses the `res` fixture of Task 1 Step 2 (Q binds E, the modified enzyme, and the modified enzyme with B bound; the two residual-bearing Q complexes are dead ends). The Q group admits two divisions, by B and by residual; each alone frees a dead-end constant, so each is a minimal gaining set and the pair is not minimal.

```julia
@testset "_expand_split_kinetic_group: by residual" begin
    # Q binds free E on the cycle and the two modified forms as dead ends. The
    # Q group divides by B ({E, E(; res)} | {E(B; res)}) and by residual
    # ({E} | {E(; res), E(B; res)}); each division frees a dead-end constant
    # on its own, so each is a child and the pair is not minimal.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + Q ⇌ E(Q), E(; residual = A - P) + Q ⇌ E(Q; residual = A - P),
             E(B; residual = A - P) + Q ⇌ E(B, Q; residual = A - P))
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
        end
    end)
    by_B = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + Q ⇌ E(Q), E(; residual = A - P) + Q ⇌ E(Q; residual = A - P))
            E(B; residual = A - P) + Q ⇌ E(B, Q; residual = A - P)
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
        end
    end)
    by_residual = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + Q ⇌ E(Q)
            (E(; residual = A - P) + Q ⇌ E(Q; residual = A - P),
             E(B; residual = A - P) + Q ⇌ E(B, Q; residual = A - P))
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
        end
    end)
    kids = EnzymeRates._expand_split_kinetic_group(m)
    @test length(kids) == 2
    @test Set(kids) == Set([by_B, by_residual])
end
```

- [ ] **Step 3: Run the file.** Focused run; both testsets 2/2. If either child set differs, print `EnzymeRates.steps.(kids)` and stop: the likely cause is the gain test (`_partition_independent_count`) counting a dead-end constant as tied, which would be a finding for Denis, not a fixture error.

- [ ] **Step 4: Commit.**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Pin split children from conformation and residual contexts"
```

---

### Task 5: Say it in the docs and the test rules

**Files:**
- Modify: `docs/src/identify/enumeration_engine.md` (the `init_mechanisms` section and "Modeling choices")
- Modify: `.claude/CLAUDE.md` ("Enumeration-engine tests")

- [ ] **Step 1: `init_mechanisms` section.** After the sentence ending `since it would split the reaction into two disconnected half-cycles.` add: "Chemistry is always written as an isomerization: the enzyme isomerizes to a product-bound form, and the release is its own step. Every move takes the isomerization steps as the chemistry steps, so a mechanism written for the derivation with chemistry folded into a release step, as on the [Ping-pong mechanisms](@ref) page, derives correctly but is not a valid parent for `expand_mechanisms`, which rejects a binding step that changes the covalent residual. A binding step may change the enzyme's conformation; a conformation change is not a chemical reaction."

- [ ] **Step 2: Modeling choices.** Add a paragraph before "**Competitive-inhibitor binding stays at equilibrium.**": "**Chemistry is the isomerization step.** The moves recognize a chemistry step by its having no ligand on it. The enumerator writes every mechanism that way, and `expand_mechanisms` rejects a parent that folds chemistry into a release; the derivation still accepts that form for hand-written textbook mechanisms."

- [ ] **Step 3: Test rules.** In `.claude/CLAUDE.md` under "### Enumeration-engine tests" change the closing sentence to "File-level test helpers are prefixed `_testhelper_` so they cannot be mistaken for package functions; a one-line closure local to a testset needs no prefix." and add rule 4 after rule 3: "4. **Write chemistry the enumerator's way.** A move fixture writes each chemistry step as an isomerization to a product-bound form followed by a release step (`E(A) <--> E(P; residual = A - P)` then `E(P; residual = A - P) ⇌ E(; residual = A - P) + P`), never folded into the release; `expand_mechanisms` rejects the folded form and the moves misread it."

- [ ] **Step 4: Build the docs** (foreground, no other Julia running, 20-minute timeout): `julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'`. Expected: clean build. Commit only the two edited files.

```bash
git add docs/src/identify/enumeration_engine.md .claude/CLAUDE.md
git commit -m "Docs and test rules: chemistry is the isomerization step"
```

---

### Task 6: Full suite

Denis's ruling after this plan: the giant `@testset "Mechanism Enumeration"` block is broken up as the last step after round 3 (and a round 4 if one is needed), in its own commit, so that the move diffs stay uncontaminated. Not part of this plan.

- [ ] **Step 1: Run the full suite detached** (command under "Running tests") and poll in bounded loops. Expected: all pass; the count is the round-2 count (24048) minus the assertions removed by Task 1 plus those added by Tasks 1, 3 and 4. Report the exact summary line. A `signal: KILL` means the OOM killer; rerun after checking `free -m` shows about 4 GB available, never by killing the VS Code language server.

- [ ] **Step 2: Report** the summary line, the commits, and the ledger.

---

## Self-review

- Spec coverage: question 1 (end-to-end conformation and residual split tests) is Task 4; question 2 and 4 (one chemistry definition, `_is_chemistry` gone, engine guarded, derivation untouched) are Tasks 2, 3, 5; the conformation comment is the second half of Task 3's test and the helper's docstring; Task 1 makes the round-2 fixtures obey the rule the docs state.
- Type consistency: `_assert_chemistry_is_iso(m::Union{Mechanism, AllostericMechanism})` is defined in Task 3 and called only there; `is_iso` is the existing `src/types.jl` accessor.
- Derivations to verify at execution: Task 1 Step 3 (four flip children), Task 1 Step 4 (six regulator children), Task 4 (one and two split children). Each is derived above from the move's code; a mismatch stops the task.
