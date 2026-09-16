# Enumeration Review Round 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close Denis's second review: chemistry steps recognized by content rather than by the iso flag (fixing the flip move on combined chemistry-release steps and stopping inhibitor-bound forms from catalyzing), context splits by conformation and residual, an `expand_mechanisms` time budget, the TestEnv recipe in `CLAUDE.md`, the docs that describe these, and a measured, formula-backed count table for the combinatorics page.

**Architecture:** One new predicate `_is_chemistry(s::Step)` in `src/mechanism_enumeration.jl` used by `_flux_carrying_groups` and by the regulator move's mirror loop. `_context_bipartitions` gains two context families over the context form: its conformation and its residual. Tests follow the round-1 rules: inline macro fixtures, exact child sets, placed after the giant testset. The combinatorics measurement is a script whose output is a committed Markdown report; the page rewrite waits for Denis's review of that report.

**Tech Stack:** Julia 1.12, EnzymeRates internals, `Test`, TestEnv for focused runs.

**Spec:** Denis's review of 2026-09-15 (second message: items 1, 3a, 3b, 4, 5, 7 and the follow-up decisions a–d) on top of `docs/superpowers/specs/2026-09-11-group-set-flips-and-context-splits-design.md`.

## Global Constraints

- 92-character lines, 4-space indentation, match surrounding style; every new test inline with the mechanism macro; move tests assert exact child sets; helpers prefixed `_testhelper_`; new testsets go after `end # top-level testset` in `test/test_mechanism_enumeration.jl` (CLAUDE.md "Enumeration-engine tests").
- `rate_equation` untouched. No numerical identifiability test in `src/`.
- Never edit an expectation to match output; print the actual result and stop.
- Focused runs via TestEnv (below); one Julia process at a time; no background monitors; the full suite needs about 4 GB free and is run detached.
- Commit after each task with the attribution trailer from the session's system reminder.

## Running tests

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_mechanism_enumeration.jl")'
```

About 4 minutes. The `UndefVarError: random_reduced_params` in `init division-freeness (bi_bi_pp)` is a known cross-file artifact. Full suite: `(nohup julia --project -e 'using Pkg; Pkg.test()' > /tmp/enzymerates_fullsuite6.log 2>&1 &)` then poll.

## Fixtures used below

The docs' ping-pong bi-bi (from `docs/src/deriving/ping_pong.md`), chemistry written as combined chemistry-release steps:

```julia
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end)
```

Its four groups are singletons: A-binding (RE), chemistry 1 (SS, releases P), B-binding (RE), chemistry 2 (SS, releases Q). The reaction with an inhibitor:

```julia
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        competitive_inhibitors: I
    end
```

---

### Task 1: Chemistry steps by content

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_flux_carrying_groups` ~line 1319)
- Test: `test/test_mechanism_enumeration.jl` (after `end # top-level testset`; the `_flux_carrying_groups` testset and the two wrapper testsets)

**Interfaces:**
- Produces: `_is_chemistry(s::Step)::Bool` — true when the step changes the enzyme's covalent state or its bound set by anything other than its own bound metabolite; false for a pure binding or release step.

- [ ] **Step 1: Write the failing tests**

Add to the `_flux_carrying_groups` testset:

```julia
    # Chemistry written as a combined chemistry-release step (no iso step at all):
    # every group lies on the catalytic cycle and is flux-carrying.
    pp = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end)
    @test all(EnzymeRates._flux_carrying_groups(pp))
```

Add to `@testset "_expand_re_to_ss (group-set flips)"`:

```julia
    @testset "_expand_re_to_ss: ping-pong with combined chemistry-release steps" begin
        m = <the docs ping-pong fixture above>
        flipA = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(; residual = A - P) + P
                E(; residual = A - P) + B ⇌ E(B; residual = A - P)
                E(B; residual = A - P) <--> E + Q
            end
        end)
        flipB = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(; residual = A - P) + P
                E(; residual = A - P) + B <--> E(B; residual = A - P)
                E(B; residual = A - P) <--> E + Q
            end
        end)
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 2
        @test Set(kids) == Set([flipA, flipB])
    end
```

- [ ] **Step 2: Run to verify they fail**

Expected: `all(_flux_carrying_groups(pp))` fails (no iso step, so no group is flux-carrying); the ping-pong flip test fails with 0 children.

- [ ] **Step 3: Implement**

Add next to `_flux_carrying_groups`:

```julia
"""
    _is_chemistry(s::Step) -> Bool

A chemistry step changes the enzyme's covalent residual, or changes its bound
set by something other than the step's own bound metabolite. A pure binding or
release step changes neither. Iso steps are chemistry; so is a combined
chemistry-release step such as `E(A) <--> E(; residual = A - P) + P`.
"""
function _is_chemistry(s::Step)
    is_iso(s) && return true
    residual(from_species(s)) == residual(to_species(s)) || return true
    bm = bound_metabolite(s)
    from_b = Metabolite[b for b in bound(from_species(s))]
    to_b = Metabolite[b for b in bound(to_species(s))]
    # remove one occurrence of the bound metabolite from whichever side holds it
    i = findfirst(==(bm), to_b)
    i === nothing || deleteat!(to_b, i)
    j = findfirst(==(bm), from_b)
    (i === nothing && j !== nothing) && deleteat!(from_b, j)
    sort(string.(name.(from_b))) != sort(string.(name.(to_b)))
end
```

Replace `is_iso(s)` with `_is_chemistry(s)` in `_flux_carrying_groups` (the `edge_is_iso` collection) and rename that local to `edge_is_chemistry`; update its docstring's "Chemistry steps are the isomerization steps" to say chemistry steps are those `_is_chemistry` recognizes. Do not touch the regulator move: Denis decided (2026-09-15) that an inhibitor-bound form may carry out a half-reaction whose ligands it does not compete with; Task 1b pins that.

- [ ] **Step 4: Run the file**

Expected: the new tests pass; every existing count unchanged (the enumerated seeds write chemistry as iso steps, which `_is_chemistry` also recognizes; the sequential regulator tests have no chemistry endpoints eligible for an inhibitor, so their counts do not move). If a pre-existing regulator-move count changes, print the mechanism and stop.

- [ ] **Step 5: Commit**

```bash
git commit -am "Recognize chemistry steps by content"
```

---

### Task 1b: Inhibitor-bound forms catalyze at most one half-reaction, never the net reaction

**Files:**
- Test: `test/test_mechanism_enumeration.jl` (new top-level testsets after the split wrapper)

Denis's decision: the competition rule (an inhibitor competes with at least one substrate and at least one product) is kept as is. Its consequence is the invariant to pin: an inhibitor-bound form can be mirrored onto a chemistry step whose ligands it does not compete with (a half-reaction in ping-pong), but no cycle of inhibitor-bound forms contains every chemistry step, because the competing substrate must bind somewhere and that form can never carry the inhibitor.

- [ ] **Step 1: Exact children of the docs ping-pong with an inhibitor**

The plain mechanism macro declares regulators with the `regulators:` label and writes an inhibitor-bound form as `E(I::Inh)` (see `_parse_plain_mechanism_body` in `src/dsl.jl` and the `E(Lactate::Inh)` forms in `test/test_identify_rate_equation.jl`). The regulator move returns children whose `reaction` carries the reaction's atoms plus the inhibitor, which a macro mechanism without atoms does not, so compare the canonical group lists: `Set(EnzymeRates.steps.(kids)) == Set(EnzymeRates.steps.(expected))`, and assert every child's reaction equals `EnzymeRates._add_competitive_inhibitor(rxn, :I)`.

```julia
@testset "_expand_add_dead_end_regulator: ping-pong, inhibitor mirrors one half-reaction" begin
    # An inhibitor competes with at least one substrate and one product. In
    # ping-pong the half-reaction whose ligands it does not compete with can still
    # run on the inhibitor-bound modified enzyme, so its chemistry step is mirrored;
    # the other half is barred by the competing substrate, so the net reaction never
    # runs with I bound. The five patterns below are: I on E; I on E and E(B; res),
    # mirroring the second half; I on E(A) and E(; res), mirroring the first half;
    # I on E(; res); I on E and E(; res).
    rxn = <the ping-pong reaction with competitive_inhibitors: I>
    m = <the docs ping-pong fixture>
    mech(block) = EnzymeRates.Mechanism(block)
    only_E = mech(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E + I ⇌ E(I::Inh)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end)
    E_and_EBres = mech(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            (E + I ⇌ E(I::Inh), E(B; residual = A - P) + I ⇌ E(B, I::Inh; residual = A - P))
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            (E(B; residual = A - P) <--> E + Q,
             E(B, I::Inh; residual = A - P) <--> E(I::Inh) + Q)
        end
    end)
    EA_and_Eres = mech(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            (E(A) + I ⇌ E(A, I::Inh), E(; residual = A - P) + I ⇌ E(I::Inh; residual = A - P))
            (E(A) <--> E(; residual = A - P) + P,
             E(A, I::Inh) <--> E(I::Inh; residual = A - P) + P)
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end)
    only_Eres = mech(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(; residual = A - P) + I ⇌ E(I::Inh; residual = A - P)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end)
    E_and_Eres = mech(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            (E + I ⇌ E(I::Inh), E(; residual = A - P) + I ⇌ E(I::Inh; residual = A - P))
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(kids) == 5
    @test Set(EnzymeRates.steps.(kids)) ==
          Set(EnzymeRates.steps.([only_E, E_and_EBres, EA_and_Eres, only_Eres, E_and_Eres]))
    @test all(c -> EnzymeRates.reaction(c) == EnzymeRates._add_competitive_inhibitor(rxn, :I), kids)
end
```

If `regulators: I` with `E(I::Inh)` does not produce the same `Species` as the move (e.g. a different regulator type or name rendering), report the two `steps` values that differ and stop; do not weaken the comparison.

- [ ] **Step 2: The invariant over the enumerated ping-pong seeds**

```julia
@testset "_expand_add_dead_end_regulator: no inhibitor-bound form runs the net reaction" begin
    # Aggregate pin over the ping-pong seed set: for every regulator child, the
    # chemistry steps mirrored onto inhibitor-bound forms are a strict subset of
    # the chemistry steps, so no cycle of inhibitor-bound forms completes the
    # reaction. Follows from competition with at least one substrate.
    rxn = <the ping-pong reaction with competitive_inhibitors: I>
    inhibitor_bound(sp) = any(b -> b isa EnzymeRates.CompetitiveInhibitor, EnzymeRates.bound(sp))
    n_children = 0; n_half = 0
    for m in EnzymeRates.init_mechanisms(rxn), c in EnzymeRates._expand_add_dead_end_regulator(m, rxn)
        chem = [s for grp in EnzymeRates.steps(c) for s in grp if EnzymeRates._is_chemistry(s)]
        mirrored = count(s -> inhibitor_bound(EnzymeRates.from_species(s)), chem)
        plain = length(chem) - mirrored
        @test mirrored < plain
        n_children += 1; mirrored > 0 && (n_half += 1)
    end
    @test n_children > 0
    @test n_half > 0          # the half-reaction mirror really occurs on the seeds
end
```

- [ ] **Step 3: Run the file; commit**

```bash
git commit -am "Pin the inhibitor invariant: at most one half-reaction with inhibitor bound"
```

---

### Task 2: Context bipartitions by conformation and residual

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_context_bipartitions`)
- Test: `test/test_mechanism_enumeration.jl` (`@testset "_context_bipartitions"` and the 102 pin)

- [ ] **Step 1: Write the failing tests**

Add to `@testset "_context_bipartitions"`:

```julia
    # Residual as a context: Q binds free E, the modified enzyme, and the
    # modified enzyme with B bound. Context B divides {E, E(;res)} from
    # {E(B;res)}; the residual divides {E} from the two modified forms.
    res = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + Q ⇌ E(Q), E(; residual = A - P) + Q ⇌ E(Q; residual = A - P),
             E(B; residual = A - P) + Q ⇌ E(B, Q; residual = A - P))
            E + A ⇌ E(A)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end)
    q_group = only(grp for grp in EnzymeRates.steps(res) if length(grp) == 3)
    rbps = EnzymeRates._context_bipartitions(q_group)
    has_res(s) = EnzymeRates.has_residual(EnzymeRates._context_form(s))
    @test length(rbps) == 2
    Er, EBr = Symbol("E_res_+A_-P"), Symbol("EB_res_+A_-P")
    @test src_forms(rbps[1][1]) == Set([:E, Er]) && src_forms(rbps[1][2]) == Set([EBr])  # by B
    @test src_forms(rbps[2][1]) == Set([:E]) && src_forms(rbps[2][2]) == Set([Er, EBr])  # residual

    # Conformation as a context: A binds E, Estar, and Estar(B).
    conf = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), Estar + A ⇌ Estar(A), Estar(B) + A ⇌ Estar(A, B))
            E(A) <--> Estar(A)
            Estar(A, B) <--> Estar(P, Q)
            Estar(P, Q) ⇌ Estar(P) + Q
            Estar(P) ⇌ E + P
        end
    end)
    a_group = only(grp for grp in EnzymeRates.steps(conf) if length(grp) == 3)
    cbps = EnzymeRates._context_bipartitions(a_group)
    @test length(cbps) == 2
    @test src_forms(cbps[1][1]) == Set([:E, :Estar]) && src_forms(cbps[1][2]) == Set([:EstarB])
    @test src_forms(cbps[2][1]) == Set([:E]) && src_forms(cbps[2][2]) == Set([:Estar, :EstarB])
```

Use the names `EnzymeRates.name` renders (check once, then delete the print). Keep the 102-children pin: measured on today's seeds, every group that spans a residual has two forms whose residual context coincides with a ligand context, so the pin must not move; if it does, print the seed and stop.

- [ ] **Step 2: Run to verify they fail**

Expected: `length(rbps) == 2` and `length(cbps) == 2` fail (today each gives 1).

- [ ] **Step 3: Implement**

In `_context_bipartitions`, build the context list from three families over `_context_form(s)`: every other bound ligand (as today), the conformation, and the residual. Each context value defines the division "steps whose context form carries it" against the rest; keep only divisions with both parts nonempty; dedupe by partition; order ligands first (by role then name), then conformations (by name), then residuals (by `string`). Rewrite the docstring: "…by binding context: another ligand already bound, the enzyme's conformation, or its covalent residual. Encodes 'the affinity for this metabolite may depend on what else is bound, on the conformation, or on the covalent state.'" A group whose forms share one conformation and one residual gets no bipartition from those families.

- [ ] **Step 4: Run the file; commit**

```bash
git commit -am "Context bipartitions by conformation and residual"
```

---

### Task 3: `expand_mechanisms` time budget

**Files:**
- Test: `test/test_mechanism_enumeration.jl` (new top-level testset after the split wrapper)

- [ ] **Step 1: Add the test**

```julia
@testset "expand_mechanisms: ter-ter random-order seed within budget" begin
    # Aggregate pin over the seed set: the seed with the most steps (55) is the
    # enumeration's worst case; measured 13 s for all seven moves.
    terter = @enzyme_reaction begin
        substrates: A[C], B[N], C[O]
        products: P[C], Q[N], R[O]
    end
    seeds = EnzymeRates.init_mechanisms(terter)
    worst = seeds[argmax([EnzymeRates.n_steps(m) for m in seeds])]
    @test EnzymeRates.n_steps(worst) == 55
    t = @elapsed kids = EnzymeRates.expand_mechanisms([worst], terter)
    @test length(kids) == 81
    @test t < 60
end
```

If the child count is not 81 after Tasks 1–2, print the count per move (`_expand_re_to_ss`, `_expand_split_kinetic_group`, `_expand_add_dead_end_regulator`, `_expand_to_allosteric`, `_expand_add_allosteric_regulator`, `_expand_change_allo_state`, `_expand_merge_regulatory_sites`) and stop; do not edit 81.

- [ ] **Step 2: Run; commit**

```bash
git commit -am "Time budget on expand_mechanisms for the worst ter-ter seed"
```

---

### Task 4: TestEnv recipe in CLAUDE.md

**Files:**
- Modify: `.claude/CLAUDE.md` "Commands" section

- [ ] **Step 1: Add the recipe** after the full-suite command:

```markdown
# Focused run of one test file (warm: skips the full suite's precompile and JIT; ~4 min for
# the enumeration file). TestEnv activates the test environment in place.
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")'
```

and, below the code block:

```markdown
Use the focused run while iterating and the full suite once before committing. Four things to know: the focused run skips Aqua and JET; a helper defined in another test file (e.g. `random_reduced_params` in `test_rate_eq_derivation.jl`) is undefined under it, so an `UndefVarError` for such a helper is an artifact of the focused run; run one Julia process at a time (the machine has 7.7 GB and no swap, and the full suite needs about 4 GB free); and never park on a background monitor — poll a log in a bounded shell loop.
```

`TestEnv` must be installed in the default environment once: `julia -e 'using Pkg; Pkg.add("TestEnv")'`; say so in the same paragraph.

- [ ] **Step 2: Commit**

```bash
git commit -am "Record the TestEnv focused-run recipe in CLAUDE.md"
```

---

### Task 5: Docs for the chemistry rule and the new contexts

**Files:**
- Modify: `docs/src/identify/enumeration_engine.md` (move 2, move 3, "Modeling choices")

- [ ] **Step 1: Move 2** — replace the sentence defining context with: "The division is always by context: the steps whose enzyme form already carries some other ligand Y, sits in a given conformation, or carries a given covalent residual, against the rest. It encodes the hypothesis that the affinity for a metabolite depends on what is bound next to it, on the enzyme's conformation, or on its covalent state; with Y a competitive inhibitor it separates a catalytic step from its inhibitor-bound mirror."

- [ ] **Step 2: Move 3** — the "Mirror steps" bullet becomes: "**Mirror steps.** If the inhibitor binds two enzyme forms that a catalytic step already connects, a mirror step is added between the two inhibitor-bound forms, so the inhibitor-bound branch stays connected to the cycle. Each mirror inherits its counterpart's kinetic group and adds no parameter. Because the inhibitor competes with at least one substrate, the form that binds that substrate never carries the inhibitor, so the inhibitor-bound branch can never complete the net reaction. In a ping-pong mechanism it can carry out the one half-reaction whose ligands the inhibitor does not compete with."

- [ ] **Step 3: Modeling choices** — replace the competitive-inhibitor paragraph with: "**Competitive-inhibitor binding stays at equilibrium.** An inhibitor bound to two forms that a catalytic step connects carries flux through its mirror step, so keeping its binding at rapid equilibrium is a modeling choice ("inhibitor binding is fast"). What competition does decide is turnover: the inhibitor competes with at least one substrate and one product, so the form that binds the competing substrate never carries it, and no inhibitor-bound branch completes the net reaction. A ping-pong mechanism is the one case with more than one chemistry step; there an inhibitor that competes with the first half-reaction's ligands can still let the modified enzyme run the second half with the inhibitor bound, which is the two-site picture, and the inhibitor must leave before the next cycle." Also extend the splits paragraph: "A group is divided only by whether another ligand is already bound, by conformation, or by covalent residual."

- [ ] **Step 4: Build the docs** (`julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'`, foreground, no other Julia running), expect a clean build; commit.

```bash
git commit -am "Docs: chemistry by content, inhibitor-bound forms never catalyze, contexts by conformation and residual"
```

---

### Task 6: Combinatorics measurement report

**Files:**
- Create: `docs/superpowers/specs/2026-09-15-combinatorics-measurement.md` (the report)
- Create: `docs/superpowers/specs/2026-09-15-combinatorics-measurement.jl` (the script that produced it, so the numbers can be regenerated)

The page rewrite is not part of this task; Denis reviews the report first.

- [ ] **Step 1: Write the script.** Reaction: bi-bi with one allosteric regulator and one competitive inhibitor:

```julia
rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    allosteric_regulators: R
    competitive_inhibitors: I
end
```

Compute and print, as Markdown tables:

1. `init_mechanisms(rxn)`: count; the number of catalytic topologies (`EnzymeRates._catalytic_topologies(rxn)`), the number of dead-end patterns (`length(EnzymeRates._competition_patterns(Set([:A,:B]), Set([:P,:Q])))`), and the count of seeds per topology (min, max), with the topology that has the fewest and the most seeds printed as step lists.
2. `seed_mechanisms(rxn, Set([:R]), Set([:I]))`: count, and how many are allosteric.
3. For each of the seven moves, over the seed set from item 2: number of children per parent — minimum, median, maximum — and the parent mechanism at the minimum and at the maximum, printed as grouped step lists (use the same printing helper as the round-1 derivation: `(step, step)` for groups, `⇌`/`<-->` for RE/SS, `E(A, B)` for forms, allosteric tags after `::`).
4. Item 3 again for the depth-1 children as parents (sample at most 200 parents deterministically: every k-th child in canonical order), so the depth dependence is visible.
5. The size of the flip closure and of the split closure from the three seeds with the most steps, and the number of distinct RE segmentations in the flip closure (as in round-1 measurement M2).

Print the RE-segment count and the independent-parameter count next to every printed mechanism.

- [ ] **Step 2: Run it** (foreground, generous timeout; the allosteric moves compile nothing, so this is enumeration time only; if it exceeds 15 minutes, reduce the depth-2 sample to 100 and say so in the report). Save stdout as the report with a short header stating the reaction, the branch commit, and the date.

- [ ] **Step 3: Commit** both files.

```bash
git add docs/superpowers/specs/2026-09-15-combinatorics-measurement.*
git commit -m "Measure per-move child counts for the combinatorics page"
```

---

### Task 7: Full suite

- [ ] Run the full suite detached and poll; expected 0 failures, 0 errors. Report the summary line verbatim. No commit.
