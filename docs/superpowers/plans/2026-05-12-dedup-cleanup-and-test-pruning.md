# Dedup cleanup and test pruning — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** consolidate all rate-equation dedup tests into one new top-level testset in `test/test_mechanism_enumeration.jl`, move the rate-equation canonicalizer source from `identify_rate_equation.jl` into `mechanism_enumeration.jl`, and drop the slow CSV-replay test plus the `dedup_investigation/` directory. **Pure refactor — NO behavior changes, NO API changes.**

**Architecture:** every step preserves passing tests. New tests are added against the existing code location (where the canonicalizer currently lives), code moves second, old tests are deleted last. TDD's "failing test first" ordering does NOT apply — this is a refactor, not new behavior. Tests are added as coverage migration.

**Tech Stack:** Julia 1.x, `Test.@testset`, `EnzymeRates.@enzyme_reaction` / `@enzyme_mechanism` DSL, `EnzymeMechanism(metabolites, reactions)` low-level constructor.

---

## Investigation findings baked into this plan

REPL investigation across nine topologies established:

1. **Pass-1 absorption** (user-defined kinetic-group merges) is observable
   on a 6-step random bi-uni exemplar — verified working in REPL.
2. **Pass-2 single-symbol Wegscheider absorption** does not fire on real
   or minimal mechanisms. Even the LDH 11-step Pattern-A literal — the
   canonicalizer's original motivating case — emits no
   `# Wegscheider constraints:` section. The transitive-closure code in
   `_build_kinetic_rename_map:152-167` is defensive infrastructure for a
   case that doesn't arise; not testable from minimal exemplars.
3. **Polynomial-equivalence detection** (different enzyme-form graphs,
   equivalent v polynomials) is the canonicalizer's main job and is
   verified only on the LDH Pattern-A m_a/m_b literal pair. Smaller
   topologies tried in investigation produced different canonical text
   even with matching kinetic-group partitions; the equivalence requires
   non-trivial graph overlap.

This is why the test set differs from the spec's original Source A / B / C
framing. See spec for the corrected understanding.

---

## Final testset structure

NEW top-level `@testset` block, inserted at the END of
`test/test_mechanism_enumeration.jl` AFTER the file's existing closing
`end # top-level testset` (line 4537). It is a peer to the existing
top-level testsets at lines 139, 182, 250, 487 — NOT nested inside
`@testset "Mechanism Enumeration"`.

```julia
@testset "Rate-equation canonical hash dedup"
    @testset "_factor_sort_key sort order"
    @testset "_sort_run_factors sort order"
    @testset "Hash is deterministic across repeated calls"
    @testset "Hash hex string is 16 lowercase hex chars"
    @testset "Distinct mechanisms produce distinct hashes"
    @testset "Pass-1 kinetic-group merge: User-defined section + canonical text invariant"
    @testset "LDH Pattern-A: graph-distinct mechanisms with equivalent v hash equally"
    @testset "Allosteric T-state K_i_T renamed away in canonical hash"
    @testset "rate_equation_string emits section labels"
    @testset "Hash-equivalent mechanisms share fitted_params shape"
end
```

---

## Files touched

- **Modify (new top-level block appended):** `test/test_mechanism_enumeration.jl`.
- **Modify (move canonicalizer source out, dedup-related lines deleted):** `src/identify_rate_equation.jl`.
- **Modify (canonicalizer source appended):** `src/mechanism_enumeration.jl`.
- **Modify (T-state allosteric canonicalizer assertion removed; lines 4496-4502):** `test/test_mechanism_enumeration.jl`.
- **Modify (two canonical-hash testsets at lines 492-559 deleted):** `test/test_identify_rate_equation.jl`.
- **Modify (two include lines deleted):** `test/runtests.jl`.
- **Delete (whole files):** `test/test_eq_hash_dedup.jl`, `test/test_dedup_csv_replay.jl`.
- **Delete (whole dir):** `dedup_investigation/` and every file in it (including untracked `refit_results.csv` — confirmed OK to delete).

---

## Task ordering invariant

Every commit boundary leaves `julia --project -e 'using Pkg; Pkg.test()'` passing.

1. Add new top-level testset scaffold.
2-8. Fill in sub-tests one at a time.
9. Move canonicalizer source from `identify_rate_equation.jl` → `mechanism_enumeration.jl`.
10. Delete the two old canonical-hash testsets in `test_identify_rate_equation.jl`.
11. Delete `test_eq_hash_dedup.jl` + `test_dedup_csv_replay.jl` + their `include` lines in `runtests.jl`.
12. Delete `dedup_investigation/` directory.
13. Final verification.

---

### Task 1: Add new top-level dedup testset scaffold

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (append AFTER the existing closing `end # top-level testset` line at ~4537)

- [ ] **Step 1: Locate the file's last line**

Run: `tail -3 test/test_mechanism_enumeration.jl`
Expected: ends with `end # top-level testset`. Confirm this is the close of the `@testset "Mechanism Enumeration"` block (which starts at line 487).

- [ ] **Step 2: Append a new top-level testset AFTER the existing close**

Append at the end of the file (peer placement, not nested):

```julia

# ═══════════════════════════════════════════════════════════════════════
# Rate-equation canonical hash dedup
# ═══════════════════════════════════════════════════════════════════════

@testset "Rate-equation canonical hash dedup" begin
    # Sub-tests added in subsequent tasks.
end
```

- [ ] **Step 3: Run the full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes; new empty testset reports `0 passed, 0 failed`.

- [ ] **Step 4: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Add empty top-level dedup testset scaffold

New top-level @testset block (peer to existing top-level
blocks, not nested in "Mechanism Enumeration") that will
hold all rate-equation canonical hash dedup tests.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: `_factor_sort_key` and `_sort_run_factors` pure-unit tests

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (inside the new top-level dedup block)

Both tests assert **observable ordering only**, not tuple-shape internals
(per reviewer feedback — slot-shape assertions lock implementation
details).

- [ ] **Step 1: Replace the placeholder comment with the two unit testsets**

Replace `# Sub-tests added in subsequent tasks.` with:

```julia
    @testset "_factor_sort_key sort order" begin
        # p_i atoms sort numerically by index, not lexicographically.
        @test EnzymeRates._factor_sort_key("p_1") <
              EnzymeRates._factor_sort_key("p_2")
        @test EnzymeRates._factor_sort_key("p_2") <
              EnzymeRates._factor_sort_key("p_10")
        # Non-p_i atoms (metabolite names, E_total) sort after p_i atoms.
        @test EnzymeRates._factor_sort_key("p_99") <
              EnzymeRates._factor_sort_key("E_total")
    end

    @testset "_sort_run_factors sort order" begin
        @test EnzymeRates._sort_run_factors("p_3 * p_1 * p_2") ==
              "p_1 * p_2 * p_3"
        # Exponents preserved on their atom.
        @test EnzymeRates._sort_run_factors("p_2 ^ 2 * p_1") ==
              "p_1 * p_2 ^ 2"
        # Non-p atoms sort to end.
        @test EnzymeRates._sort_run_factors("S * p_1 * p_2") ==
              "p_1 * p_2 * S"
    end
```

- [ ] **Step 2: Run the test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: 6 new `@test` invocations pass.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Test: _factor_sort_key and _sort_run_factors sort order

Observable-ordering assertions only (no tuple-shape probes
that would lock implementation details).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Hash determinism, hex shape, and distinct-mechanism tests

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Append three testsets after the `_sort_run_factors` testset**

```julia
    @testset "Hash is deterministic across repeated calls" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E_S
                E_S <--> E_P
                E + P ⇌ E_P
            end
        end
        h1 = EnzymeRates._canonical_rate_eq_hash(m)
        h2 = EnzymeRates._canonical_rate_eq_hash(m)
        @test h1 == h2
        @test h1 isa UInt64
    end

    @testset "Hash hex string is 16 lowercase hex chars" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E_S
                E_S <--> E_P
                E + P ⇌ E_P
            end
        end
        h_full, h_hex, name_map = EnzymeRates._canonical_rate_eq_hash_data(m)
        @test h_full isa UInt64
        @test length(h_hex) == 16
        @test all(c -> c in "0123456789abcdef", h_hex)
        @test name_map isa Dict{String, String}
    end

    @testset "Distinct mechanisms produce distinct hashes" begin
        # Uni-uni vs ordered bi-uni: genuinely different v polynomials.
        m_uu = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E_S
                E_S <--> E_P
                E + P ⇌ E_P
            end
        end
        m_bu = @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A ⇌ E_A
                E_A + B ⇌ E_A_B
                E_A_B <--> E_P
                E + P ⇌ E_P
            end
        end
        @test EnzymeRates._canonical_rate_eq_hash(m_uu) !=
              EnzymeRates._canonical_rate_eq_hash(m_bu)
    end
```

- [ ] **Step 2: Run the test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: 8 new `@test` invocations pass.

If the `@enzyme_mechanism` parse fails for `m_bu`, adjust step syntax
by referring to bi-uni examples in `test/test_mechanism_enumeration.jl`
around lines 924+.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Test: canonical hash determinism, hex shape, distinct mechanisms

Three minimal testsets on hand-written uni-uni and bi-uni
exemplars. Covers basic invariants previously asserted in
the deleted test_identify_rate_equation.jl 'canonical rate-
equation hash: basic' block.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: "Pass-1 kinetic-group merge" test (verified random bi-uni)

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

The exemplar is verified in REPL: a 6-step random bi-uni with
substrate-side mirror sharing (both A-binding steps share `kg=1`, both
B-binding steps share `kg=2`). Pass 1 absorbs `K2 → K1` and `K4 → K3`,
rendered as `# User defined constraints:` with `(substituted into v)`
lines.

- [ ] **Step 1: Append the testset**

```julia
    @testset "Pass-1 kinetic-group merge: User-defined section + canonical text invariant" begin
        # Random bi-uni. Both A-binding steps share kg=1 (mirror), both
        # B-binding steps share kg=2. Pass 1 absorbs K_mirror -> K_rep.
        m = EnzymeMechanism(
            ((:A, :B), (:P,), ()),
            (((:E, :A), (:E_A,), true, 1),
             ((:E_B, :A), (:E_A_B,), true, 1),
             ((:E, :B), (:E_B,), true, 2),
             ((:E_A, :B), (:E_A_B,), true, 2),
             ((:E_A_B,), (:E_P,), false, 3),
             ((:E, :P), (:E_P,), true, 4)))

        s = rate_equation_string(m)
        @test occursin("# User defined constraints:", s)
        @test occursin("(substituted into v)", s)

        # Canonical text invariant: no raw K_i / k_if / k_ir tokens
        # survive — every parameter is renamed to p_i.
        canon, _ = EnzymeRates._canonicalize_rate_eq_with_map(m)
        @test !occursin(r"\bK\d+\b", canon)
        @test !occursin(r"\bk\d+[fr]\b", canon)
    end
```

- [ ] **Step 2: Run the test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: 4 new `@test` invocations pass.

This exemplar was REPL-verified to produce:
```
# User defined constraints:
K2 = K1  (substituted into v)
K4 = K3  (substituted into v)
# Haldane constraints:
k5r = (1 / Keq) * (1 / K1) * (1 / K3) * k5f * K6
v = E_total * (...) / (...)
```
and a canonical text with no raw K_i / k_if / k_ir tokens.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Test: Pass-1 kinetic-group merge folds K into v

Random bi-uni with substrate-side mirror sharing. Renders
'# User defined constraints:' section with (substituted
into v) lines; canonical text has no raw K_i / k_if / k_ir
tokens (regression for first-appearance rename completeness).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: "LDH Pattern-A" polynomial-equivalence test (verified hash collision)

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

The two LDH literals `m_a` and `m_b` differ in their enzyme-form graphs
(different Lactate-binding paths via NADH vs NAD-bound forms) but
produce equivalent v polynomials after Pass 1. REPL-verified hash
collision: both produce `4711f0996b051276`. This is the smallest
verified Pattern-A graph-distinct-but-v-equivalent case after nine
investigation attempts; 11 steps is the minimal known case for this
property.

- [ ] **Step 1: Append the testset**

```julia
    @testset "LDH Pattern-A: graph-distinct mechanisms with equivalent v hash equally" begin
        # 11-step LDH mechanism. m_a and m_b differ in step ordering AND
        # in which intermediate forms appear (m_a has Lactate-binding
        # via E_NADH only; m_b adds a Lactate-binding via E_NAD path).
        # After Pass 1 absorption their v polynomials are equivalent.
        # 11 steps is the minimal known case for this property — smaller
        # hand-tuned topologies tried in investigation did not exhibit
        # the graph-distinct-but-v-equivalent pattern.
        m_a = EnzymeMechanism(
            ((:NADH, :Pyruvate), (:Lactate, :NAD), ()),
            (((:E, :Lactate), (:E_Lactate,), true, 1),
             ((:E, :NAD), (:E_NAD,), true, 2),
             ((:E, :NADH), (:E_NADH,), true, 3),
             ((:E, :Pyruvate), (:E_Pyruvate,), true, 4),
             ((:E_Lactate, :NAD), (:E_Lactate_NAD,), true, 2),
             ((:E_Lactate, :NADH), (:E_Lactate_NADH,), true, 3),
             ((:E_NAD, :Pyruvate), (:E_NAD_Pyruvate,), true, 4),
             ((:E_NADH, :Lactate), (:E_Lactate_NADH,), true, 1),
             ((:E_NADH, :Pyruvate), (:E_NADH_Pyruvate,), true, 4),
             ((:E_NADH_Pyruvate,), (:E_Lactate_NAD,), false, 5),
             ((:E_Pyruvate, :NAD), (:E_NAD_Pyruvate,), true, 2)))

        m_b = EnzymeMechanism(
            ((:NADH, :Pyruvate), (:Lactate, :NAD), ()),
            (((:E, :Lactate), (:E_Lactate,), true, 1),
             ((:E, :NAD), (:E_NAD,), true, 2),
             ((:E, :NADH), (:E_NADH,), true, 3),
             ((:E, :Pyruvate), (:E_Pyruvate,), true, 4),
             ((:E_Lactate, :NADH), (:E_Lactate_NADH,), true, 3),
             ((:E_NAD, :Lactate), (:E_Lactate_NAD,), true, 1),
             ((:E_NAD, :Pyruvate), (:E_NAD_Pyruvate,), true, 4),
             ((:E_NADH, :Lactate), (:E_Lactate_NADH,), true, 1),
             ((:E_NADH, :Pyruvate), (:E_NADH_Pyruvate,), true, 4),
             ((:E_NADH_Pyruvate,), (:E_Lactate_NAD,), false, 5),
             ((:E_Pyruvate, :NAD), (:E_NAD_Pyruvate,), true, 2)))

        @test EnzymeRates._canonical_rate_eq_hash(m_a) ==
              EnzymeRates._canonical_rate_eq_hash(m_b)
    end
```

- [ ] **Step 2: Run the test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: 1 new `@test` passes. REPL-verified collision:
`_canonical_rate_eq_hash(m_a) == _canonical_rate_eq_hash(m_b) == 0x4711f0996b051276`.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Test: canonicalizer collapses LDH Pattern-A duplicates

Two 11-step LDH mechanisms with graph-distinct topologies
(different Lactate-binding intermediate paths) produce
equivalent v polynomials after Pass 1 and must hash equally.
This is the polynomial-equivalence regression the canonicalizer
was built for; smaller hand-tuned topologies tried in
investigation did not exhibit the property.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Move Allosteric T-state canonicalizer assertion into dedup block

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Append the testset to the dedup block**

```julia
    @testset "Allosteric T-state K_i_T renamed away in canonical hash" begin
        # K-type allosteric uni-uni: catalytic step is :OnlyR
        # (`_t_state_dead == true`), but binding steps are :NonequalRT,
        # so K1_T and K2_T live in `den_T` of the rate equation body.
        # Canonicalizer invariant: every parameter token must rename
        # away — no raw `_T` suffixed names survive.
        m = @allosteric_mechanism begin
            substrates: S
            products: P
            site(:catalytic, 2): begin
                steps: begin
                    E_c + S ⇌ E_S    :: NonequalRT
                    E_c + P ⇌ E_P    :: NonequalRT
                    E_S <--> E_P     :: OnlyR
                end
            end
        end
        canon, _ = EnzymeRates._canonicalize_rate_eq_with_map(m)
        @test !occursin(r"\bK\d+_T\b", canon)
        @test !occursin(r"\bk\d+[fr]_T\b", canon)
    end
```

- [ ] **Step 2: Remove the canonicalizer-invariant assertions from the original allosteric block**

In `test/test_mechanism_enumeration.jl`, find the `@testset "t_state_dead with :NonequalRT: K_T in body must be in parameters(Full)"` block (currently around line 4470). Delete these lines from inside that testset (the comment block + 3 lines of canonicalizer assertions, currently at ~4497-4502):

```julia
        # Canonicalizer invariant: every parameter token in the
        # body must be renamed away. After canonicalization, no
        # raw `_T` suffixed names should survive.
        canon, _ = EnzymeRates._canonicalize_rate_eq_with_map(m)
        @test !occursin(r"\bK\d+_T\b", canon)
        @test !occursin(r"\bk\d+[fr]_T\b", canon)
```

The `@test :K1_T in params_full` and `@test :K2_T in params_full` lines (parameters-API regressions) stay in their original block.

- [ ] **Step 3: Run the test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes. Original `t_state_dead with :NonequalRT` testset now contains only parameters-API assertions.

- [ ] **Step 4: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Move allosteric T-state canonicalizer assertion into dedup block

Parameters(Full) regressions stay in their original
allosteric testset (parameters-API). Only the canonicalizer
invariant moves to the consolidated dedup block.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: `rate_equation_string` section-labels test

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

Note: `# Wegscheider constraints:` is NOT asserted here — investigation
confirmed it doesn't appear on minimal hand-synthesizable mechanisms.
The LDH Pattern-A test (Task 5) is the indirect regression for the
section-stripping pipeline as a whole.

- [ ] **Step 1: Append the testset**

```julia
    @testset "rate_equation_string emits section labels" begin
        # User-defined section emitted when a mechanism has shared
        # kinetic_groups (reuses the Pass-1 exemplar from Task 4).
        m_user = EnzymeMechanism(
            ((:A, :B), (:P,), ()),
            (((:E, :A), (:E_A,), true, 1),
             ((:E_B, :A), (:E_A_B,), true, 1),
             ((:E, :B), (:E_B,), true, 2),
             ((:E_A, :B), (:E_A_B,), true, 2),
             ((:E_A_B,), (:E_P,), false, 3),
             ((:E, :P), (:E_P,), true, 4)))
        s_user = rate_equation_string(m_user)
        @test occursin("# User defined constraints:", s_user)

        # Haldane section: any RE binding mechanism with Keq has it.
        # Investigation confirmed that # Wegscheider constraints: is
        # not emitted on minimal mechanisms; LDH Pattern-A (Task 5)
        # is the indirect regression for that section's stripping.
        m_hal = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E_S
                E_S <--> E_P
                E + P ⇌ E_P
            end
        end
        s_hal = rate_equation_string(m_hal)
        @test occursin("# Haldane constraints:", s_hal)
    end
```

- [ ] **Step 2: Run the test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: 2 new `@test` invocations pass.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Test: rate_equation_string emits User-defined and Haldane sections

# Wegscheider constraints: not asserted here — investigation
confirmed it doesn't appear on minimal mechanisms; LDH
Pattern-A (Task 5) is the indirect regression for the
section-stripping pipeline.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: "Hash-equivalent mechanisms share fitted_params shape" test

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

Reuses the LDH Pattern-A `m_a`/`m_b` pair (the only verified
hash-equivalent pair). `_fp_kind` helper goes inside a `let` block
to keep its binding scope-local (Julia's `function f(...) end` inside
a `@testset begin ... end` leaks to the file's module scope per
reviewer feedback).

- [ ] **Step 1: Append the testset**

```julia
    @testset "Hash-equivalent mechanisms share fitted_params shape" begin
        # LDH Pattern-A pair (same as Task 5). Hash-equivalent →
        # fitted_params count and kind multiset must match.
        m_a = EnzymeMechanism(
            ((:NADH, :Pyruvate), (:Lactate, :NAD), ()),
            (((:E, :Lactate), (:E_Lactate,), true, 1),
             ((:E, :NAD), (:E_NAD,), true, 2),
             ((:E, :NADH), (:E_NADH,), true, 3),
             ((:E, :Pyruvate), (:E_Pyruvate,), true, 4),
             ((:E_Lactate, :NAD), (:E_Lactate_NAD,), true, 2),
             ((:E_Lactate, :NADH), (:E_Lactate_NADH,), true, 3),
             ((:E_NAD, :Pyruvate), (:E_NAD_Pyruvate,), true, 4),
             ((:E_NADH, :Lactate), (:E_Lactate_NADH,), true, 1),
             ((:E_NADH, :Pyruvate), (:E_NADH_Pyruvate,), true, 4),
             ((:E_NADH_Pyruvate,), (:E_Lactate_NAD,), false, 5),
             ((:E_Pyruvate, :NAD), (:E_NAD_Pyruvate,), true, 2)))

        m_b = EnzymeMechanism(
            ((:NADH, :Pyruvate), (:Lactate, :NAD), ()),
            (((:E, :Lactate), (:E_Lactate,), true, 1),
             ((:E, :NAD), (:E_NAD,), true, 2),
             ((:E, :NADH), (:E_NADH,), true, 3),
             ((:E, :Pyruvate), (:E_Pyruvate,), true, 4),
             ((:E_Lactate, :NADH), (:E_Lactate_NADH,), true, 3),
             ((:E_NAD, :Lactate), (:E_Lactate_NAD,), true, 1),
             ((:E_NAD, :Pyruvate), (:E_NAD_Pyruvate,), true, 4),
             ((:E_NADH, :Lactate), (:E_Lactate_NADH,), true, 1),
             ((:E_NADH, :Pyruvate), (:E_NADH_Pyruvate,), true, 4),
             ((:E_NADH_Pyruvate,), (:E_Lactate_NAD,), false, 5),
             ((:E_Pyruvate, :NAD), (:E_NAD_Pyruvate,), true, 2)))

        # Sanity: hash-equivalent (verified in Task 5).
        @test EnzymeRates._canonical_rate_eq_hash(m_a) ==
              EnzymeRates._canonical_rate_eq_hash(m_b)

        # _fp_kind classifies a fitted_params symbol into its canonical
        # kind — invariant under rep-step renaming. `let` keeps the
        # binding scope-local; `function _fp_kind ... end` here would
        # leak the name to the file's module scope.
        let
            _fp_kind = function(s::Symbol)
                str = string(s)
                is_T = endswith(str, "_T")
                base = is_T ? str[1:end-2] : str
                kind = if startswith(base, "K") && length(base) > 1 && isdigit(base[2])
                    :K
                elseif startswith(base, "k") && length(base) > 1 && isdigit(base[2])
                    endswith(base, "f") ? :kf : endswith(base, "r") ? :kr : :other
                elseif s == :L
                    :L
                elseif startswith(str, "K_")
                    :K_reg
                else
                    :other
                end
                is_T ? Symbol(kind, :_T) : kind
            end

            fp_a = EnzymeRates.fitted_params(m_a)
            fp_b = EnzymeRates.fitted_params(m_b)
            @test length(fp_a) == length(fp_b)
            @test sort(_fp_kind.(fp_a)) == sort(_fp_kind.(fp_b))
        end
    end
```

- [ ] **Step 2: Run the test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: 3 new `@test` invocations pass.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Test: hash-equivalent LDH Pattern-A mechanisms share fitted_params shape

Uses the verified m_a/m_b literal pair. _fp_kind helper
inside a let-block to keep its binding scope-local.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 9: Move canonicalizer source from `identify_rate_equation.jl` → `mechanism_enumeration.jl`

**Files:**
- Modify: `src/identify_rate_equation.jl` (delete the canonicalizer block)
- Modify: `src/mechanism_enumeration.jl` (append the same block)

The block being moved contains five definitions:
`_canonicalize_rate_eq_with_map`, `_sort_run_factors`, `_factor_sort_key`,
`_canonical_rate_eq_hash_data`, `_canonical_rate_eq_hash`. They sit
together between the `_Stage1Failure` constructor (above) and the
`_project_cached_params` docstring (below).

**Use text-boundary matching, NOT absolute line numbers** (per reviewer
feedback — line numbers shift silently if any earlier edit lands first).

- [ ] **Step 1: Identify the block boundaries by text**

Run: `grep -n '^"""$\|^function _canonical\|^function _factor_sort_key\|^function _sort_run_factors\|^function _project_cached_params' src/identify_rate_equation.jl`

This locates the docstring openers (`"""` at column 0) and the function definitions. The canonicalizer block to move is:
- **Start**: the `"""` that opens the `_canonicalize_rate_eq_with_map` docstring (currently at line 128 — the docstring text starts with `Build the canonical text + name_map.`).
- **End**: the closing `end` of `_canonical_rate_eq_hash` (currently at line 289 — immediately followed by a blank line then `"""` opening `_project_cached_params`'s docstring).

- [ ] **Step 2: Capture the block**

Run: `awk '/^"""$/{n++} n==2{print; if(/^end$/ && captured) exit} /Build the canonical text/{captured=1} captured' src/identify_rate_equation.jl > /tmp/canonicalizer_block.jl`

Or, equivalently and more robustly, use sed bounded by text:
```bash
sed -n '/^"""$/{:a;N;/Build the canonical text/!{D;ba}};/_project_cached_params/q;p' src/identify_rate_equation.jl
```

If sed text-bounded extraction is fragile, fall back to: Read the file in the editor, locate the block boundaries by eye using the function-definition anchors, copy the block contents.

Verify `/tmp/canonicalizer_block.jl` has 5 function definitions: `_canonicalize_rate_eq_with_map`, `_sort_run_factors`, `_factor_sort_key`, `_canonical_rate_eq_hash_data`, `_canonical_rate_eq_hash`. Run:

```bash
grep -c '^function ' /tmp/canonicalizer_block.jl
```
Expected: `5`.

- [ ] **Step 3: Delete the block from `identify_rate_equation.jl`**

Use the Edit tool to remove the captured block. The Edit tool's exact-match
requirement protects against drift: find the unique start string
(`"""\nBuild the canonical text + name_map. Internal helper exposed`)
and unique end string (`function _canonical_rate_eq_hash(m::AbstractEnzymeMechanism)\n    first(_canonical_rate_eq_hash_data(m))\nend`) and delete everything from start through end (inclusive).

After deletion, the file should pick up at the `_project_cached_params`
docstring opener immediately after the deleted block.

- [ ] **Step 4: Append the block to `mechanism_enumeration.jl`**

Append at the end of `src/mechanism_enumeration.jl` (after the existing
`dedup!` closing `end`):

```

# ─── Rate-equation canonical hash ──────────────────────────────────────
```

Then append the contents of `/tmp/canonicalizer_block.jl`.

- [ ] **Step 5: Run the full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes — the dedup block's tests (which call
`EnzymeRates._canonical_rate_eq_hash` via qualified name) continue to
work from the new location.

If `UndefVarError: ANNOTATION_SUBSTITUTED` appears at load time, swap
the include order in `src/EnzymeRates.jl` (put `rate_eq_derivation.jl`
before `mechanism_enumeration.jl` — already the case currently).

- [ ] **Step 6: Commit**

```bash
git add src/identify_rate_equation.jl src/mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Move rate-equation canonicalizer to mechanism_enumeration.jl

_canonicalize_rate_eq_with_map, _canonical_rate_eq_hash_data,
_canonical_rate_eq_hash, _sort_run_factors, _factor_sort_key
moved verbatim from identify_rate_equation.jl. No behavior
change.

_project_cached_params stays in identify_rate_equation.jl
(beam-search-specific).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 10: Delete the two old canonical-hash testsets in `test_identify_rate_equation.jl`

**Files:**
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Locate the two testsets by name**

Run: `grep -n '^@testset "canonical' test/test_identify_rate_equation.jl`
Expected output:
- `@testset "canonical rate-equation hash: basic" begin` at line 492.
- `@testset "canonical hash collapses Pattern-A LDH duplicates" begin` somewhere between 492 and 559.

Verify by:
```bash
sed -n '492,494p;557,560p' test/test_identify_rate_equation.jl
```
Expected:
- Line 492: `@testset "canonical rate-equation hash: basic" begin`
- Line 559: `end` (close of `canonical hash collapses Pattern-A LDH duplicates`)
- Line 560: blank
- Line 561 or so: next testset (`@testset "_onesided_permutation_p" begin`)

- [ ] **Step 2: Delete the two testsets**

Use the Edit tool to remove the contents from the opening `@testset "canonical rate-equation hash: basic" begin` through the closing `end` of the Pattern-A testset (inclusive of both). Also remove the blank line(s) between them.

- [ ] **Step 3: Run the full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes. The deleted assertions are covered by the
new dedup block (Hash deterministic / hex / Distinct + LDH Pattern-A).

- [ ] **Step 4: Commit**

```bash
git add test/test_identify_rate_equation.jl
git commit -m "$(cat <<'EOF'
Drop canonical-hash testsets superseded by dedup block

'canonical rate-equation hash: basic' and 'canonical hash
collapses Pattern-A LDH duplicates' are covered by the
new 'Rate-equation canonical hash dedup' top-level block
in test_mechanism_enumeration.jl.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 11: Delete `test_eq_hash_dedup.jl`, `test_dedup_csv_replay.jl`, update `runtests.jl`

**Files:**
- Delete: `test/test_eq_hash_dedup.jl`
- Delete: `test/test_dedup_csv_replay.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Delete the two test files**

```bash
git rm test/test_eq_hash_dedup.jl test/test_dedup_csv_replay.jl
```

- [ ] **Step 2: Remove the two `include` lines from `runtests.jl`**

Use the Edit tool to remove these two lines from `test/runtests.jl`:

```julia
    include("test_eq_hash_dedup.jl")
    include("test_dedup_csv_replay.jl")
```

The final `runtests.jl` should match the current file minus those two lines.

- [ ] **Step 3: Run the full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes. Test-suite runtime drops noticeably (CSV
replay no longer runs).

- [ ] **Step 4: Commit**

```bash
git add test/runtests.jl
git commit -m "$(cat <<'EOF'
Drop test_eq_hash_dedup.jl and test_dedup_csv_replay.jl

CSV-replay test (~22k mechanism-type compiles, >30min runtime)
and cv_results.csv-driven cluster test are replaced by the
LDH Pattern-A literal pair (Tasks 5 and 8) plus minimal
exemplar tests in the consolidated dedup block.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 12: Delete `dedup_investigation/` directory

**Files:**
- Delete: `dedup_investigation/` (whole directory, ~22 MB; includes the currently-untracked `refit_results.csv` — user confirmed OK to delete)

- [ ] **Step 1: Confirm no source/test still references the directory**

Run: `grep -rn dedup_investigation test/ src/`
Expected: no output.

- [ ] **Step 2: Remove the directory**

```bash
git rm -r dedup_investigation/
rm -rf dedup_investigation/
```

The `git rm -r` removes tracked files; the `rm -rf` removes the
untracked `refit_results.csv` and the now-empty directory.

- [ ] **Step 3: Confirm gone**

Run: `ls dedup_investigation/ 2>&1`
Expected: `ls: cannot access 'dedup_investigation/': No such file or directory`

- [ ] **Step 4: Run the full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes.

- [ ] **Step 5: Commit**

```bash
git commit -m "$(cat <<'EOF'
Delete dedup_investigation/ artifacts

CSVs (LDH_data.csv, cv_results.csv, params_estimate_{5,6,7,8}.csv),
status docs (status_2026-05-11.md, investigation_eq_hash_duplication.md),
refit_diagnostic.jl, and untracked refit_results.csv are no
longer referenced by any test. Findings preserved in git
history.

Frees ~22 MB.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 13: Final verification

**Files:** none modified — verification only.

- [ ] **Step 1: Run the full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes. Aqua + JET clean.

- [ ] **Step 2: Verify acceptance criteria**

Run these checks:

```bash
# Canonicalizer function DEFINITIONS in mechanism_enumeration.jl
grep -n "^function _canonical\|^function _factor_sort_key\|^function _sort_run_factors\|^function _canonicalize_rate_eq" src/mechanism_enumeration.jl
# Expected: lines for _canonicalize_rate_eq_with_map, _sort_run_factors,
# _factor_sort_key, _canonical_rate_eq_hash_data, _canonical_rate_eq_hash

# Same names NOT defined (only called) in identify_rate_equation.jl
grep -n "^function _canonical\|^function _factor_sort_key\|^function _sort_run_factors\|^function _canonicalize_rate_eq" src/identify_rate_equation.jl
# Expected: no output

# Old test files gone
ls test/test_eq_hash_dedup.jl test/test_dedup_csv_replay.jl 2>&1
# Expected: 'No such file or directory' for both

# Investigation dir gone
ls dedup_investigation/ 2>&1
# Expected: 'No such file or directory'

# runtests.jl includes are clean
grep "test_eq_hash_dedup\|test_dedup_csv_replay" test/runtests.jl
# Expected: no output

# No new package exports
git diff main..HEAD -- src/EnzymeRates.jl
# Expected: no changes to the exports list
```

- [ ] **Step 3: Confirm no API surface change**

The public API (`parameters`, `fitted_params`, `rate_equation`,
`rate_equation_string`, `EnzymeMechanism`, `AllostericEnzymeMechanism`,
`identify_rate_equation`, beam-search behavior) is unchanged. The
canonicalizer is still accessible via `EnzymeRates._canonical_rate_eq_hash`
(unexported, private), now sourced from `mechanism_enumeration.jl`.

If all checks pass, this plan is complete. Task 12's commit is the
last state-changing one.
