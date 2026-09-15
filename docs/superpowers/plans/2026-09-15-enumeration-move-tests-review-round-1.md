# Enumeration Move Tests, Review Round 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the tests of the rewritten flip and split moves reviewable and exact: every fixture written out with `@enzyme_mechanism`, every child set compared against expected mechanisms, ter-ter covered, all in `test/test_mechanism_enumeration.jl`.

**Architecture:** Test-only change plus a short rules section in `CLAUDE.md`. No `src/` change. The expected children below were derived by running the moves at branch tip `8aa9ff3` and reviewed by hand; each split child corresponds to one four-cycle of the rapid-equilibrium graph, and each flip child to one metabolite's groups (or the minimal pair once a metabolite is split).

**Tech Stack:** Julia 1.12, `Test`, EnzymeRates internals (`Mechanism`, `@enzyme_mechanism`, `@allosteric_mechanism`, `_expand_re_to_ss`, `_expand_split_kinetic_group`, `_context_bipartitions`, `_flux_carrying_groups`).

**Spec:** Denis's review of 2026-09-15 (twelve items) on top of `docs/superpowers/specs/2026-09-11-group-set-flips-and-context-splits-design.md`.

## Global Constraints

- 92-character line length (characters, not bytes), 4-space indentation; match the surrounding test style.
- Every enumeration-move test defines its mechanism inline with `@enzyme_mechanism` or `@allosteric_mechanism` so a reviewer sees the mechanism in the testset. The only exceptions are the aggregate regression pins over `init_mechanisms` (220 flips, 102 splits, 553 segment checks, 228 counter checks), which are about the seed set itself and say so in a comment.
- Every move test asserts the exact child set: `Set(children) == Set(expected)` with expected mechanisms written out, plus `length(children) == n`. `Mechanism` equality is structural and canonical, so two mechanisms written in different step orders compare equal.
- Test helpers defined in test files are prefixed `_testhelper_`.
- Never edit an expectation to match output. If a derived expectation below is wrong, report the actual children (print them) and stop.
- Run focused tests via `TestEnv.activate()` + `include`; `test/test_mechanism_enumeration.jl` takes a few minutes; run Julia in the foreground with a 20-minute timeout, never two Julia processes at once, no background monitors. A `UndefVarError: random_reduced_params` in `init division-freeness (bi_bi_pp)` under a focused run is a known cross-file artifact.
- Commit after every task with the attribution trailer from the session's system reminder.

## Running tests

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_mechanism_enumeration.jl")'
```

Full suite (Task 6 only), detached and polled:

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && (nohup julia --project -e 'using Pkg; Pkg.test()' > /tmp/enzymerates_fullsuite3.log 2>&1 &)
```

## File structure

- `test/test_mechanism_enumeration.jl` — receives everything from `test/test_flip_and_split_moves.jl`: the helpers and their testsets go into a new section `# ─── Expansion-move helpers ───` placed immediately before `@testset "_expand_re_to_ss"` (~line 1407); flip tests go inside `@testset "_expand_re_to_ss"`; split tests go inside `@testset "_expand_split_kinetic_group"`.
- `test/test_flip_and_split_moves.jl` — deleted; its `include` line in `test/runtests.jl` removed.
- `CLAUDE.md` — gains an "Enumeration tests" subsection under "Instructions for this repository".

---

### Task 1: Move the new tests into `test_mechanism_enumeration.jl`, rename the helper, fix the mirror comment

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (~line 1400, before `@testset "_expand_re_to_ss"`; inside that testset; inside `@testset "_expand_split_kinetic_group"` ~line 2018)
- Delete: `test/test_flip_and_split_moves.jl`
- Modify: `test/runtests.jl` (remove the include)

- [ ] **Step 1: Relocate, no behavior change**

Cut the whole content of `test/test_flip_and_split_moves.jl` below its `using` lines and paste it as follows:

- The constants `_bibi_rxn`, `_uni_uni_rxn`, the helpers `_flip_groups`, `_identifiable_rank`, `_closure`, and the const fixtures `_allo_uni_uni`, `_random_bibi`, together with the testsets `_independent_param_count …`, `_flux_carrying_groups`, `_re_segment_count`, `_re_segment_count_after_flip …`, `_minimal_gaining_sets`, `_context_bipartitions`, `_context_bipartitions separates …`, `_apply_bipartitions`, `_partition_independent_count …`: into a new section headed

```julia
# ─── Expansion-move helpers ──────────────────────────────────────────────
# Helpers and building blocks of the RE→SS flip move and the kinetic-group
# split move. Every fixture is written out with the mechanism macro so the
# mechanism under test is visible in the testset.
```

placed immediately before `@testset "_expand_re_to_ss" begin`.

- Every testset whose name starts with `_expand_re_to_ss:` goes to the end of `@testset "_expand_re_to_ss"` (before its closing `end`).
- Every testset whose name starts with `_expand_split_kinetic_group:` goes to the end of `@testset "_expand_split_kinetic_group"`.

Rename `_bibi_rxn` → `_testhelper_bibi_rxn`, `_uni_uni_rxn` → keep using the file's existing `uni_uni_rxn` constant (defined near line 94) and delete the duplicate, `_flip_groups` → `_testhelper_flip_groups`, `_identifiable_rank` → `_testhelper_identifiable_rank`, `_closure` → `_testhelper_closure`. `_allo_uni_uni` and `_random_bibi` are removed in Task 2; leave them for now.

Delete `test/test_flip_and_split_moves.jl` and the line `include("test_flip_and_split_moves.jl")` in `test/runtests.jl`.

- [ ] **Step 2: Fix the mirror-fixture comment in the `_flux_carrying_groups` testset**

The comment above the `mirror` fixture says "A pendant binding-only square"; it is not pendant. Replace the comment with:

```julia
    # An inhibitor bound to both E and E(S), with the mirror E(I)+S ⇌ E(I,S): the
    # inhibitor square shares the edge E→E(S) with the catalytic cycle, so the
    # cycle E→E(I)→E(I,S)→E(S)→E(P)→E passes through chemistry and every step,
    # inhibitor binding included, carries net flux. Contrast the pendant square
    # below, whose inhibitors bind only free E.
```

and make sure the pendant fixture's comment says it is the contrast: "Two inhibitors binding only free E form a binding-only square joined to the cycle at the single form E: no cycle through it contains chemistry, so its four groups are not flux-carrying."

- [ ] **Step 3: Run the enumeration file**

Expected: every relocated testset passes with the same counts as before (flip seeds 220, splits 102, segment checks 553, counter checks 228).

- [ ] **Step 4: Commit**

```bash
git add test/test_mechanism_enumeration.jl test/runtests.jl
git rm test/test_flip_and_split_moves.jl
git commit -m "Move the flip and split move tests into the enumeration test file"
```

---

### Task 2: Inline every fixture with the mechanism macro

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (the relocated testsets)

For each relocated testset that builds its mechanism from `first(EnzymeRates.init_mechanisms(...))`, `_allo_uni_uni`, or `_random_bibi`, replace the construction with an inline macro definition. The mechanisms are:

Uni-uni seed (replaces `first(init_mechanisms(uni_uni_rxn))`):

```julia
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
```

Allosteric uni-uni (replaces `_allo_uni_uni`):

```julia
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: EqualAI
            E + P ⇌ E(P)      :: EqualAI
        end
    end)
```

Random-order bi-bi with dead ends (replaces `_random_bibi`):

```julia
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(P) + A ⇌ E(A, P))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q), E(A) + P ⇌ E(A, P))
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q), E(B) + Q ⇌ E(B, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
```

- [ ] **Step 1: Replace each use**

Testsets to edit: `_independent_param_count …` (keep the loop over `MECHANISM_TEST_SPECS`, it is a spec-wide pin), `_flux_carrying_groups` (the first `m`), `_re_segment_count` (both branches), `_re_segment_count_after_flip …` (the allosteric part; the 55-seed loop stays as an aggregate pin with a comment `# Aggregate pin over the seed set; the inline fixtures below cover the shapes.`), `_expand_re_to_ss: every child raises …` (rewrite in Task 3), `_expand_re_to_ss: a group with no flux-carrying step never flips` (already inline), `_expand_re_to_ss: uni-uni flips are emitted` (inline the uni-uni), `_expand_re_to_ss: no emitted set is a superset` (rewrite in Task 3), `_expand_split_kinetic_group: random-order bi-bi reaches …` (inline `_random_bibi`), `… a rejected split is the parent's model` (inline `_random_bibi`), `… every child gains and no set contains another` (rewrite in Task 4). Delete the `_allo_uni_uni` and `_random_bibi` constants once unused.

- [ ] **Step 2: Run the enumeration file; commit**

```bash
git commit -am "Write every enumeration-move fixture inline with the mechanism macro"
```

---

### Task 3: Exact-child tests for the flip move

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` inside `@testset "_expand_re_to_ss"`

Replace the relocated testsets `_expand_re_to_ss: every child raises the RE segment count`, `… two segment-flat groups flip together`, and `… no emitted set is a superset of another` with the exact tests below. Keep `… seed child count is unchanged (220 …)`, `… a group with no flux-carrying step never flips` (extend it per (c)), and `… uni-uni flips are emitted` (replace it with (a), which subsumes it).

- [ ] **Step 1: Write the tests**

(a) Uni-uni: two children, exact.

```julia
    @testset "_expand_re_to_ss: uni-uni emits both single-group flips exactly" begin
        # Uni-uni steady-state and rapid-equilibrium laws have the same form, so
        # these two children are documented no-ops; the move still emits them
        # because their rejection has no structural proof.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
            end
        end)
        expected = [
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S <--> E(S)
                    E(S) <--> E(P)
                    E + P ⇌ E(P)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S ⇌ E(S)
                    E(S) <--> E(P)
                    E + P <--> E(P)
                end
            end),
        ]
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 2
        @test Set(kids) == Set(expected)
    end
```

(b) Random-order bi-bi, one group per metabolite: four children, exact. Write the parent and four expected children out in full (each child flips exactly one of the A, B, P, Q groups to `<-->`; the parent is the bi-bi random-order mechanism with groups `(E + A ⇌ E(A), E(B) + A ⇌ E(A, B))`, `(E + B ⇌ E(B), E(A) + B ⇌ E(A, B))`, `(E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))`, `(E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))`, `E(A, B) <--> E(P, Q)`). Assert `length(kids) == 4` and `Set(kids) == Set(expected)`.

(c) Dead-end leaf: two children, exact; the leaf group stays RE in both.

```julia
    @testset "_expand_re_to_ss: a dead-end leaf group never flips" begin
        # E(P) + S ⇌ E(P, S) is a bridge: no cycle through it contains chemistry,
        # so at steady state it carries no net flux and flipping it would add a
        # parameter the data cannot see. Only the two catalytic bindings flip.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
                E(P) + S ⇌ E(P, S)
            end
        end)
        expected = [
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S <--> E(S)
                    E(S) <--> E(P)
                    E + P ⇌ E(P)
                    E(P) + S ⇌ E(P, S)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S ⇌ E(S)
                    E(S) <--> E(P)
                    E + P <--> E(P)
                    E(P) + S ⇌ E(P, S)
                end
            end),
        ]
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 2
        @test Set(kids) == Set(expected)
    end
```

(d) Split metabolite: singles are segment-flat, the pair is emitted. Four children, exact.

```julia
    @testset "_expand_re_to_ss: a split metabolite flips as a pair, never singly" begin
        # A's two binding steps sit in separate groups. Flipping either alone
        # leaves E and E(A) joined through the other A step (E–E(B)–E(A,B)–E(A)),
        # so the segment count does not rise and no such child exists. The pair
        # cuts the square, so it is emitted; B, P, Q flip alone as before.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(B) + A ⇌ E(A, B)
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        flipA = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(B) + A <--> E(A, B)
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        flipB = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(B) + A ⇌ E(A, B)
                (E + B <--> E(B), E(A) + B <--> E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        flipP = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(B) + A ⇌ E(A, B)
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P <--> E(P), E(Q) + P <--> E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        flipQ = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(B) + A ⇌ E(A, B)
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 4
        @test Set(kids) == Set([flipA, flipB, flipP, flipQ])
        # The two single-A flips are absent: the old per-group move emitted them
        # as segment-flat no-ops.
        for single in (
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E(B) + A ⇌ E(A, B)
                    (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                    (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                    (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                    E(A, B) <--> E(P, Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E(B) + A <--> E(A, B)
                    (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                    (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                    (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                    E(A, B) <--> E(P, Q)
                end
            end))
            @test !(single in kids)
            @test EnzymeRates._re_segment_count(single) == EnzymeRates._re_segment_count(m)
        end
    end
```

(e) Ter-ter random substrate cube, ordered product release: six children, exact (one per group: A, B, C, R, P, Q). Parent:

```julia
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(C) + A ⇌ E(A, C), E(B, C) + A ⇌ E(A, B, C))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(C) + B ⇌ E(B, C), E(A, C) + B ⇌ E(A, B, C))
                (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C), E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R ⇌ E(R)
            end
        end)
```

Each expected child is the parent with exactly one group written with `<-->` (both arrows of a multi-step group). Write all six out. Assert `length(kids) == 6` and `Set(kids) == Set(expected)`.

(f) Ter-ter after one context split: segment-flat singles combine into pairs. Ten children, exact by flipped-group sets. The parent is the ter-ter split child in which A is split by B and B by A:

```julia
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(C) + A ⇌ E(A, C))
                (E(B) + A ⇌ E(A, B), E(B, C) + A ⇌ E(A, B, C))
                (E + B ⇌ E(B), E(C) + B ⇌ E(B, C))
                (E(A) + B ⇌ E(A, B), E(A, C) + B ⇌ E(A, B, C))
                (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C), E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R ⇌ E(R)
            end
        end)
```

Name the groups A1 = `(E + A, E(C) + A)`, A2 = `(E(B) + A, E(B, C) + A)`, B1 = `(E + B, E(C) + B)`, B2 = `(E(A) + B, E(A, C) + B)`, C, P, Q, R. The expected children are the parent with these group sets flipped: `{C}`, `{P}`, `{Q}`, `{R}`, `{A1, A2}`, `{B1, B2}`, `{A1, B1}`, `{A1, B2}`, `{A2, B1}`, `{A2, B2}`. Each single A or B part alone is segment-flat (E and E(A) stay joined through the other part), so no child flips exactly one of A1, A2, B1, B2. Build the expected children with `_testhelper_flip_groups(m, [indices])`, locating each group index by its step content (write a small local `group_of(step_string_match)` or, simpler, find each group by `Set(name.(from_species))` of its members, e.g. A1 is the A group whose source forms are `{:E, :EC}`), and comment each expected child with its group names. Assert `length(kids) == 10`, `Set(kids) == Set(expected)`, and that no child has exactly one of the four A/B parts steady state.

- [ ] **Step 2: Run the enumeration file**

Expected: PASS. If any exact set differs, print the offending children with their steps and stop; do not edit expectations.

- [ ] **Step 3: Commit**

```bash
git commit -am "Exact-child tests for the RE→SS flip move on uni-uni, bi-bi, and ter-ter"
```

---

### Task 4: Exact-child tests for the split move and context bipartitions

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` inside `@testset "_expand_split_kinetic_group"` and the helpers section

Replace the relocated `_expand_split_kinetic_group: every child gains and no set contains another` with the exact tests below; keep the 102-children pin, the 16-structure closure test (now inline), the rejected-split rank test, the allosteric parent test, the four-step 2+2 test, and the ter-ter budget test.

- [ ] **Step 1: Bi-bi random order without dead ends: two children, exact**

Parent as in Task 3(b). Expected: the A/B square fully separated, and the P/Q square fully separated:

```julia
        splitAB = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(B) + A ⇌ E(A, B)
                E + B ⇌ E(B)
                E(A) + B ⇌ E(A, B)
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        splitPQ = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                E + P ⇌ E(P)
                E(Q) + P ⇌ E(P, Q)
                E + Q ⇌ E(Q)
                E(P) + Q ⇌ E(P, Q)
                E(A, B) <--> E(P, Q)
            end
        end)
```

Comment: a single context split (A by B alone) is tied straight back by the A/B square, so no child splits one group; each child frees one interaction factor (independent count 5 → 6) by splitting both groups of one square. Assert `length(kids) == 2`, `Set(kids) == Set([splitAB, splitPQ])`, and `_independent_param_count` of each child is 6.

- [ ] **Step 2: Bi-bi random order with dead ends: four children, exact**

Parent as in Task 2's inline `_random_bibi`. Expected, one per four-cycle of the RE graph:

- A/P dead-end square (frees the E(A,P) interaction): A group `(E + A, E(B) + A)`, B group unchanged, P group `(E + P, E(Q) + P)`, Q group unchanged, plus singletons `E(A) + P ⇌ E(A, P)` and `E(P) + A ⇌ E(A, P)`.
- A/B square: A group `(E + A, E(P) + A)`, B group `(E + B, E(Q) + B)`, P and Q unchanged, plus singletons `E(A) + B ⇌ E(A, B)` and `E(B) + A ⇌ E(A, B)`.
- B/Q dead-end square: A unchanged, B group `(E + B, E(A) + B)`, P unchanged, Q group `(E + Q, E(P) + Q)`, plus singletons `E(B) + Q ⇌ E(B, Q)` and `E(Q) + B ⇌ E(B, Q)`.
- P/Q square: A and B unchanged, P group `(E + P, E(A) + P)`, Q group `(E + Q, E(B) + Q)`, plus singletons `E(P) + Q ⇌ E(P, Q)` and `E(Q) + P ⇌ E(P, Q)`.

Write all four out in full with the chemistry step `E(A, B) <--> E(P, Q)`. Assert `length(kids) == 4`, exact set, and independent count 6 for each.

- [ ] **Step 3: Ter-ter cube: three children, exact, and the A group's context bipartitions**

Parent as in Task 3(e). Expected children, one per substrate pair; e.g. the A/B child:

```julia
        splitAB = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(C) + A ⇌ E(A, C))
                (E(B) + A ⇌ E(A, B), E(B, C) + A ⇌ E(A, B, C))
                (E + B ⇌ E(B), E(C) + B ⇌ E(B, C))
                (E(A) + B ⇌ E(A, B), E(A, C) + B ⇌ E(A, B, C))
                (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C), E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R ⇌ E(R)
            end
        end)
```

and analogously `splitAC` (A split by C: `(E + A, E(B) + A)` and `(E(C) + A, E(B, C) + A)`; C split by A: `(E + C, E(B) + C)` and `(E(A) + C, E(A, B) + C)`; B unchanged) and `splitBC` (B split by C: `(E + B, E(A) + B)` and `(E(C) + B, E(A, C) + B)`; C split by B: `(E + C, E(A) + C)` and `(E(B) + C, E(A, B) + C)`; A unchanged). Comment: each child says "the affinity for one substrate depends on whether the other is bound", a 2 + 2 division of two four-step groups that singleton carves could never produce; independent count 7 → 8. Assert `length(kids) == 3`, exact set, count 8 each.

In the helpers section, extend `@testset "_context_bipartitions"` with the ter-ter A group: exactly two bipartitions, in this order (contexts sorted by role then name, B before C):

```julia
        a_group = only(grp for grp in EnzymeRates.steps(m)
                       if length(grp) == 4 &&
                          EnzymeRates.name(EnzymeRates.bound_metabolite(first(grp))) == :A)
        bps = EnzymeRates._context_bipartitions(a_group)
        @test length(bps) == 2
        src_forms(part) = Set(EnzymeRates.name(EnzymeRates.from_species(s)) for s in part)
        # By B: the B-free forms E, E(C) against the B-bound forms E(B), E(B,C).
        @test src_forms(bps[1][1]) == Set([:E, :EC])
        @test src_forms(bps[1][2]) == Set([:EB, :EBC])
        # By C: the C-free forms E, E(B) against the C-bound forms E(C), E(B,C).
        @test src_forms(bps[2][1]) == Set([:E, :EB])
        @test src_forms(bps[2][2]) == Set([:EC, :EBC])
```

(Check the species names `EnzymeRates.name(species)` renders as `:EC`, `:EBC`; if the renderer orders letters differently, use the rendered names, and say so in the report.)

Also add the bi-bi dead-end A group case: contexts B and P give two bipartitions: by B `{E, EP}` against `{EB}`; by P `{E, EB}` against `{EP}`.

- [ ] **Step 4: Run the enumeration file; commit**

```bash
git commit -am "Exact-child tests for the split move and context bipartitions, bi-bi and ter-ter"
```

---

### Task 5: Rules for future enumeration tests

**Files:**
- Modify: `CLAUDE.md` (under "Instructions for this repository", after the `rate_equation` perf section)
- Modify: `test/test_mechanism_enumeration.jl` (file header comment)

- [ ] **Step 1: Add the rules to CLAUDE.md**

```markdown
### Enumeration-engine tests

Tests of `init_mechanisms`, `seed_mechanisms`, `expand_mechanisms`, and the expansion moves live in `test/test_mechanism_enumeration.jl` and follow two rules:

1. **Write the mechanism in the testset.** Define every fixture inline with `@enzyme_mechanism` or `@allosteric_mechanism`, even when that repeats a mechanism used elsewhere. A fixture pulled from `init_mechanisms`, a shared constant, or a helper cannot be reviewed without leaving the testset. The only exception is an aggregate regression pin over a whole seed set (a count over all `init_mechanisms(rxn)`), which must say so in a comment.
2. **Assert the exact children.** A move test asserts `length(children) == n` and `Set(children) == Set(expected)` with every expected child written out as a mechanism. Property assertions (every child gains, no superset) are welcome in addition, never instead: a property test passes on wrong output that happens to satisfy the property.

Test helpers defined in test files are prefixed `_testhelper_` so they cannot be mistaken for package functions.
```

- [ ] **Step 2: Add a pointer to the test file header**

Below the file's two `# ABOUTME:` lines add:

```julia
# Fixture rules for this file are in CLAUDE.md, "Enumeration-engine tests":
# every mechanism inline via the macro, every move test asserts the exact
# child set.
```

- [ ] **Step 3: Commit**

```bash
git commit -am "Record the enumeration-test rules in CLAUDE.md"
```

---

### Task 6: Full suite

- [ ] **Step 1: Run the full suite** (detached, poll the log; no other Julia process running). Expected: 0 failures, 0 errors, and the total count higher than 24,376 by the number of added assertions.

- [ ] **Step 2: Report the summary line verbatim.** No commit.
