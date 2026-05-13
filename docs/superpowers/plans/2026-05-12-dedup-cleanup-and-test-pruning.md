# Dedup cleanup and test pruning — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** consolidate all rate-equation dedup tests into one block in `test/test_mechanism_enumeration.jl`, move the rate-equation canonicalizer source from `identify_rate_equation.jl` into `mechanism_enumeration.jl`, and drop the slow CSV-replay test plus the `dedup_investigation/` directory. No behavior changes.

**Architecture:** every step preserves passing tests. Tests are added first (against the existing code location), code moves second (callers reach the canonicalizer through the `EnzymeRates._canonical_rate_eq_hash` qualified name regardless of which file hosts it), old tests are deleted last (only after new tests cover their assertions).

**Tech Stack:** Julia 1.x, `Test.@testset`, `EnzymeRates.@enzyme_reaction` / `@enzyme_mechanism` DSL, `EnzymeMechanism(metabolites, reactions)` low-level constructor.

---

## Spec deviation flagged for confirmation

The spec describes three separate tests for the rename pipeline:
1. Steps sharing a kinetic_group → only K_rep in hash (Pass 1 only)
2. Wegscheider single-symbol tie K_i = K_j → K_i absent (Pass 2 only)
3. Chained renames close transitively (Pass 1 + Pass 2 + closure)

After re-reading `_build_kinetic_rename_map` in `src/rate_eq_derivation.jl:119-170`, Pass 2 single-symbol Wegscheider absorption only fires when the thermodynamic constraint system produces a single-symbol RHS, which empirically requires Pass 1 to have already shared kinetic_groups across cycle-equivalent steps. A pure "Pass 2 alone" exemplar may not exist for minimal hand-synthesized topologies — the LDH cluster mechanisms exhibit Pass 2 firing precisely because they have Pass 1 kinetic_group sharing first.

This plan consolidates the Wegscheider and chained-rename tests into a single test: **"Kinetic-group merge + Wegscheider absorption fold into v"**, which exercises Pass 1, Pass 2, and the transitive closure simultaneously. The Pass-1-only test (kinetic_group sharing → K_rep) stays separate, since it cleanly tests Pass 1 without depending on Pass 2 firing.

If Task 4 (the combined test) cannot find a minimal topology that triggers Pass 2 absorption, the test reduces to just Pass 1 + closure on a Pass-1-only chain (which the code still handles). The implementer notes the outcome explicitly in the commit.

---

## Final testset structure

```julia
@testset "Rate-equation canonical hash dedup"
    @testset "_factor_sort_key sorts p_i atoms numerically"
    @testset "_sort_run_factors reorders mixed run by p_i"
    @testset "Hash is deterministic across repeated calls"
    @testset "Hash hex string is 16 lowercase hex chars"
    @testset "Distinct mechanisms produce distinct hashes"
    @testset "Steps sharing a kinetic_group → only representative K_i in hash"
    @testset "Kinetic-group merge + Wegscheider absorption fold into v"
    @testset "Allosteric T-state K_i_T renamed away in canonical hash"
    @testset "rate_equation_string emits constraint section labels"
    @testset "Hash-equivalent mechanisms share fitted_params shape"
end
```

10 testsets. Pure-unit testsets at the top run in microseconds. The mechanism-driven testsets at the bottom each compile one or two distinct `EnzymeMechanism` types — total dedup-block runtime target: **<30 s**.

---

## Files touched

- **Add (new testset block):** `test/test_mechanism_enumeration.jl` — new top-level `@testset "Rate-equation canonical hash dedup"` block inserted at end of file (before the file's outermost `end # top-level testset`).
- **Modify (move canonicalizer source out, dedup-related lines deleted):** `src/identify_rate_equation.jl`.
- **Modify (canonicalizer source added at bottom):** `src/mechanism_enumeration.jl`.
- **Modify (T-state allosteric assertion moved into new block; two canonical-hash testsets deleted):** `test/test_identify_rate_equation.jl`, `test/test_mechanism_enumeration.jl`.
- **Modify (two include lines deleted):** `test/runtests.jl`.
- **Delete (whole files):** `test/test_eq_hash_dedup.jl`, `test/test_dedup_csv_replay.jl`.
- **Delete (whole dir):** `dedup_investigation/` and every file in it.

---

## Task ordering invariant

Every commit boundary leaves `julia --project -e 'using Pkg; Pkg.test()'` passing. Order:

1. Add new testset block (all assertions, complete).
2. Move Allosteric T-state canonicalizer assertion into the new block.
3. Move canonicalizer source from `identify_rate_equation.jl` → `mechanism_enumeration.jl`.
4. Delete old canonical-hash testsets in `test_identify_rate_equation.jl`.
5. Delete `test_eq_hash_dedup.jl` + `test_dedup_csv_replay.jl` + their `include` lines in `runtests.jl`.
6. Delete `dedup_investigation/` directory.
7. Final verification.

---

### Task 1: Add new dedup testset scaffold (skeleton, no assertions yet)

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (insert before the closing `end # top-level testset` on the last line)

- [ ] **Step 1: Locate the closing `end` of the outermost testset**

Run: `tail -3 test/test_mechanism_enumeration.jl`
Expected: the last non-blank line is `end # top-level testset` (currently at line 4537).

- [ ] **Step 2: Insert the empty testset block above the closing `end`**

Insert this block immediately before the `end # top-level testset` line:

```julia
# ═══════════════════════════════════════════════════════════════════════
# 7. Rate-equation canonical hash dedup
# ═══════════════════════════════════════════════════════════════════════

@testset "Rate-equation canonical hash dedup" begin
    # Sub-tests added in subsequent tasks.
end
```

- [ ] **Step 3: Run the full test file to confirm parse**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes; new empty testset reports `0 passed, 0 failed` (Julia allows empty `@testset` blocks).

- [ ] **Step 4: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Add empty dedup testset scaffold

Section 7 of test_mechanism_enumeration.jl will hold all
rate-equation canonical hash dedup tests. Filled in across
the following tasks.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: `_factor_sort_key` and `_sort_run_factors` pure-unit tests

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (inside the new `@testset "Rate-equation canonical hash dedup"` block)

- [ ] **Step 1: Add the `_factor_sort_key` test**

Replace the `# Sub-tests added in subsequent tasks.` comment with:

```julia
    @testset "_factor_sort_key sorts p_i atoms numerically" begin
        # Plain p_i atoms sort by integer index, not lex.
        @test EnzymeRates._factor_sort_key("p_1") <
              EnzymeRates._factor_sort_key("p_2")
        @test EnzymeRates._factor_sort_key("p_2") <
              EnzymeRates._factor_sort_key("p_10")
        # p_i ^ k captures the exponent in the third slot.
        @test EnzymeRates._factor_sort_key("p_2 ^ 3")[3] == 3
        @test EnzymeRates._factor_sort_key("p_2")[3] == 1
        # Non-p_i atoms (e.g. metabolite names, E_total) sort
        # after p_i atoms via the leading tag `1` vs `0`.
        @test EnzymeRates._factor_sort_key("E_total")[1] == 1
        @test EnzymeRates._factor_sort_key("p_99")[1] == 0
        @test EnzymeRates._factor_sort_key("p_99") <
              EnzymeRates._factor_sort_key("E_total")
    end

    @testset "_sort_run_factors reorders mixed run by p_i" begin
        @test EnzymeRates._sort_run_factors("p_3 * p_1 * p_2") ==
              "p_1 * p_2 * p_3"
        # Exponents preserved on the atom they belong to.
        @test EnzymeRates._sort_run_factors("p_2 ^ 2 * p_1") ==
              "p_1 * p_2 ^ 2"
        # Non-p atoms sort to the end via lex-second key.
        @test EnzymeRates._sort_run_factors("S * p_1 * p_2") ==
              "p_1 * p_2 * S"
    end
```

- [ ] **Step 2: Run the new testsets in isolation**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: 8 new `@test` invocations pass. Whole suite passes.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Test: _factor_sort_key and _sort_run_factors pure-unit

8 assertions covering numeric vs lex sort orders, exponent
preservation, and non-p_i atom positioning.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Hash determinism, hex shape, and distinct-mechanism tests

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (inside the same dedup block)

- [ ] **Step 1: Append the three testsets after the existing two**

After the `_sort_run_factors` testset's closing `end`, append:

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
        # Uni-uni vs bi-bi: genuinely different v polynomials.
        m_uu = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E_S
                E_S <--> E_P
                E + P ⇌ E_P
            end
        end
        m_bb = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_A + B ⇌ E_A_B
                E_A_B <--> E_P_Q
                E_P_Q <--> E_P + Q
                E_P <--> E + P
            end
        end
        @test EnzymeRates._canonical_rate_eq_hash(m_uu) !=
              EnzymeRates._canonical_rate_eq_hash(m_bb)
    end
```

- [ ] **Step 2: Run the test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes; the three new testsets contribute 8 `@test` assertions.

If a `@enzyme_mechanism` parse fails (e.g., bi-bi step syntax doesn't match), inspect the error and adjust the step list to match the DSL grammar (refer to existing bi-bi `@enzyme_mechanism` examples around `test/test_mechanism_enumeration.jl:924` for valid step syntax).

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Test: canonical hash determinism, hex shape, distinct mechanisms

Covers basic invariants previously asserted in the deleted
test_identify_rate_equation.jl 'canonical rate-equation hash:
basic' block. Three minimal testsets on hand-written uni-uni
and bi-bi exemplars.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: "Steps sharing a kinetic_group → only K_rep in hash" test

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (inside the same dedup block)

This test asserts the Pass-1 absorption (user-defined kinetic-group merge). It uses two specs with the same partition of kinetic_groups but different step ordering and different integer labels for the groups — proving the canonical hash is invariant to those choices.

- [ ] **Step 1: Append the testset**

After the "Distinct mechanisms" testset's `end`, append:

```julia
    @testset "Steps sharing a kinetic_group → only representative K_i in hash" begin
        # Uni-uni + dead-end inhibitor I that mirrors substrate binding.
        # Mirror step shares kinetic_group with E+S binding step → K1 only
        # in v, no K_mirror symbol.
        m1 = EnzymeMechanism(
            ((:S,), (:P,), (:I,)),
            (((:E, :S), (:E_S,), true, 1),
             ((:E_S,), (:E_P,), false, 2),
             ((:E, :P), (:E_P,), true, 3),
             ((:E, :I), (:E_I,), true, 1)))   # mirror shares group 1

        # Same partition, different step ordering, different
        # kinetic_group integer labels. Canonical hash must match.
        m2 = EnzymeMechanism(
            ((:S,), (:P,), (:I,)),
            (((:E, :I), (:E_I,), true, 5),    # rep step for shared group
             ((:E, :S), (:E_S,), true, 5),
             ((:E_S,), (:E_P,), false, 8),
             ((:E, :P), (:E_P,), true, 11)))

        @test EnzymeRates._canonical_rate_eq_hash(m1) ==
              EnzymeRates._canonical_rate_eq_hash(m2)

        # Canonical text contains only renamed `p_i` tokens; no raw
        # K2, K3, K4 / k2f / k2r should survive (those would mean the
        # canonicalizer's first-appearance rename failed).
        canon, _ = EnzymeRates._canonicalize_rate_eq_with_map(m1)
        @test !occursin(r"\bK\d+\b", canon)
        @test !occursin(r"\bk\d+[fr]\b", canon)
    end
```

- [ ] **Step 2: Run the test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: testset passes — 3 `@test` assertions.

If `m1` or `m2` construction fails with "stoichiometric rank" error, the partition might violate the mechanism validator (e.g., dead-end inhibitor I that doesn't appear in regulators). Confirm the third element of the metabolites tuple is `(:I,)` (regulators).

If the canonical-hash equality assertion fails, dump both canonical strings to diff:
```julia
println(EnzymeRates._canonicalize_rate_eq_with_map(m1)[1])
println(EnzymeRates._canonicalize_rate_eq_with_map(m2)[1])
```
and adjust the partition to match.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Test: kinetic_group sharing → only K_rep in canonical hash

Two specs with same partition but different step ordering and
different kinetic_group integers must hash identically. Plus
canonical-text invariant: no raw K_i or k_if/k_ir tokens
survive (everything renamed to p_i).

Replaces the deleted Pattern-A LDH test in
test_identify_rate_equation.jl with a minimal exemplar.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: "Kinetic-group merge + Wegscheider absorption fold into v" test

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (inside the same dedup block)

This test asserts the full Pass-1 + Pass-2 + transitive-closure pipeline. The exemplar must be a mechanism whose `_build_kinetic_rename_map` produces a Pass-2 single-symbol absorption (visible as a `K_i = K_j (substituted into v)` line in the `# Wegscheider constraints:` section of `rate_equation_string`).

- [ ] **Step 1: Find an exemplar mechanism — investigation step**

The Pass-2 absorption fires when the thermodynamic constraints produce a single-symbol equality between two binding K's. Empirically this requires Pass 1 to have set up shared kinetic_groups first.

Candidate to try first — a random bi-uni with a mirror step:

```julia
m = EnzymeMechanism(
    ((:A, :B), (:P,), ()),
    (((:E, :A), (:E_A,), true, 1),
     ((:E, :B), (:E_B,), true, 2),
     ((:E_A, :B), (:E_A_B,), true, 2),    # second-binding K2 = first B-binding K2
     ((:E_B, :A), (:E_A_B,), true, 1),    # second-binding K1 = first A-binding K1
     ((:E_A_B,), (:E_P,), false, 3),
     ((:E, :P), (:E_P,), true, 4)))
```

Run in REPL: `println(rate_equation_string(m))` and check for a `(substituted into v)` line inside `# Wegscheider constraints:`. If yes, this is the exemplar.

If the Wegscheider section is empty or has only multi-symbol RHSes, try alternative topologies:
1. Bi-bi rapid-equilibrium random with one shared kinetic_group across cross-paths.
2. Uni-uni with both substrate and product dead-end inhibitors mirroring the same binding step.
3. Bi-uni ordered with iso step shared via mirror.

The acceptance criterion is: `rate_equation_string(m)` output contains a line matching `r"^\s*K\d+\s*=\s*K\d+\s+\(substituted into v\)\s*$"` inside the `# Wegscheider constraints:` section.

Record the successful exemplar's source in a code block; use it in step 2 below.

- [ ] **Step 2: Append the testset using the exemplar from step 1**

After the previous testset's `end`, append (substitute the actual successful `EnzymeMechanism(...)` literal):

```julia
    @testset "Kinetic-group merge + Wegscheider absorption fold into v" begin
        # EXEMPLAR FROM STEP 1: replace with the actual literal that
        # successfully triggers Pass-2 single-symbol Wegscheider absorption.
        m = EnzymeMechanism(
            ((:A, :B), (:P,), ()),
            (((:E, :A), (:E_A,), true, 1),
             ((:E, :B), (:E_B,), true, 2),
             ((:E_A, :B), (:E_A_B,), true, 2),
             ((:E_B, :A), (:E_A_B,), true, 1),
             ((:E_A_B,), (:E_P,), false, 3),
             ((:E, :P), (:E_P,), true, 4)))

        # Sanity: Pass-1 absorption visible — `# User defined constraints:`
        # section emitted with at least one (substituted into v) line.
        s = rate_equation_string(m)
        @test occursin("# User defined constraints:", s)
        @test occursin("(substituted into v)", s)

        # Wegscheider section also emitted. If this assertion fails,
        # the chosen exemplar doesn't trigger Pass-2 absorption — try
        # a different topology (see Task 5 Step 1 candidate list).
        @test occursin("# Wegscheider constraints:", s)

        # Canonical text invariant: every parameter token is renamed
        # away. No raw K_i, k_if, k_ir tokens survive.
        canon, _ = EnzymeRates._canonicalize_rate_eq_with_map(m)
        @test !occursin(r"\bK\d+\b", canon)
        @test !occursin(r"\bk\d+[fr]\b", canon)
    end
```

- [ ] **Step 3: Run the test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: testset passes — 5 `@test` assertions.

If the `# Wegscheider constraints:` assertion fails, the chosen exemplar from step 1 doesn't fire Pass 2. Iterate by trying the alternative topologies listed in step 1.

If after trying all four candidate topologies no Pass-2 firing is reproducible, downgrade the test: keep the Pass-1 assertions (`# User defined constraints:` and the canonical-text invariant), drop the `# Wegscheider constraints:` assertion, and note in the commit message that Pass 2 absorption couldn't be triggered by minimal hand-synthesized topologies. The canonical-text invariant still catches transitive-closure breakage on any mechanism with multiple symbol renames.

- [ ] **Step 4: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Test: Pass-1 + Pass-2 rename pipeline folds K_i into v

Asserts user-defined kinetic-group merge AND Wegscheider
single-symbol absorption both fire, and the canonicalizer's
first-appearance rename leaves no raw K_i / k_if / k_ir
tokens in the canonical text (regression for transitive-
closure correctness).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Move Allosteric T-state canonicalizer assertion into the dedup block

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

The block at `test/test_mechanism_enumeration.jl:4470-4503` (`@testset "t_state_dead with :NonequalRT: K_T in body must be in parameters(Full)"`) contains two kinds of assertions:
- `parameters(m, Full)` regressions (`:K1_T in params_full`, `:K2_T in params_full`) — these are parameters-API regressions, keep in place.
- Canonicalizer regression (`@test !occursin(r"\bK\d+_T\b", canon)` and the `k\d+[fr]_T` variant) — these are dedup regressions, move to the new block.

- [ ] **Step 1: Append the canonicalizer-only testset to the dedup block**

After the previous testset's `end`, append:

```julia
    @testset "Allosteric T-state K_i_T renamed away in canonical hash" begin
        # K-type allosteric uni-uni: catalytic step is :OnlyR
        # (`_t_state_dead == true`), but binding steps are :NonequalRT,
        # so K1_T and K2_T live in `den_T` of the rate equation body.
        # Canonicalizer invariant: every parameter token must rename
        # away. After canonicalization no raw `_T` suffixed names
        # should survive.
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

- [ ] **Step 2: Remove the canonicalizer-invariant assertions from the original block**

Edit `test/test_mechanism_enumeration.jl:4496-4502` — delete these lines from the existing `t_state_dead with :NonequalRT` testset:

```julia
        # Canonicalizer invariant: every parameter token in the
        # body must be renamed away. After canonicalization, no
        # raw `_T` suffixed names should survive.
        canon, _ = EnzymeRates._canonicalize_rate_eq_with_map(m)
        @test !occursin(r"\bK\d+_T\b", canon)
        @test !occursin(r"\bk\d+[fr]_T\b", canon)
```

The `parameters(m, Full)` assertions (`@test :K1_T in params_full`, `@test :K2_T in params_full`) stay in place. The closing `end` of the original testset comes right after them.

- [ ] **Step 3: Run the test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes. The original `t_state_dead with :NonequalRT` testset now contains only the parameters-API assertions; the new dedup block contains the canonicalizer-invariant assertions.

- [ ] **Step 4: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Move allosteric T-state canonicalizer assertion into dedup block

Keeps the parameters(m, Full) regressions in their original
allosteric block (parameters-API regression); moves only the
canonicalizer invariant (`!occursin(r\"\\bK\\d+_T\\b\", canon)`)
into the consolidated dedup testset.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: "rate_equation_string emits constraint section labels" test

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Append the testset**

After the previous testset's `end`, append:

```julia
    @testset "rate_equation_string emits constraint section labels" begin
        # User-defined: bi-uni random with mirror shares kinetic_groups,
        # so `# User defined constraints:` section is emitted.
        m_user = EnzymeMechanism(
            ((:A, :B), (:P,), ()),
            (((:E, :A), (:E_A,), true, 1),
             ((:E, :B), (:E_B,), true, 2),
             ((:E_A, :B), (:E_A_B,), true, 2),
             ((:E_B, :A), (:E_A_B,), true, 1),
             ((:E_A_B,), (:E_P,), false, 3),
             ((:E, :P), (:E_P,), true, 4)))
        s_user = rate_equation_string(m_user)
        @test occursin("# User defined constraints:", s_user)

        # Wegscheider / Haldane: any RE binding mechanism has a Haldane
        # closure involving Keq. Uni-uni 3-step is the minimal case.
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
Expected: testset passes — 2 `@test` assertions.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Test: rate_equation_string emits constraint section labels

Asserts both '# User defined constraints:' (Pass-1 kinetic-
group merge) and '# Haldane constraints:' section headers
appear in the rendered output. Replaces the deleted
section-label test in test_eq_hash_dedup.jl.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: "Hash-equivalent mechanisms share fitted_params shape" test

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Append the testset (with `_fp_kind` helper inline)**

After the previous testset's `end`, append:

```julia
    @testset "Hash-equivalent mechanisms share fitted_params shape" begin
        # `_fp_kind` classifies a fitted_params symbol into its canonical
        # kind (invariant under rep-step renaming). Same hash bucket → same
        # count and same kind multiset.
        function _fp_kind(s::Symbol)
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

        # Reuse the hash-equivalent pair from the kinetic_group sharing
        # testset.
        m1 = EnzymeMechanism(
            ((:S,), (:P,), (:I,)),
            (((:E, :S), (:E_S,), true, 1),
             ((:E_S,), (:E_P,), false, 2),
             ((:E, :P), (:E_P,), true, 3),
             ((:E, :I), (:E_I,), true, 1)))
        m2 = EnzymeMechanism(
            ((:S,), (:P,), (:I,)),
            (((:E, :I), (:E_I,), true, 5),
             ((:E, :S), (:E_S,), true, 5),
             ((:E_S,), (:E_P,), false, 8),
             ((:E, :P), (:E_P,), true, 11)))

        # Precondition: the pair is hash-equivalent. If this fails, Task 4
        # broke first.
        @test EnzymeRates._canonical_rate_eq_hash(m1) ==
              EnzymeRates._canonical_rate_eq_hash(m2)

        fp1 = EnzymeRates.fitted_params(m1)
        fp2 = EnzymeRates.fitted_params(m2)
        @test length(fp1) == length(fp2)
        @test sort(_fp_kind.(fp1)) == sort(_fp_kind.(fp2))
    end
```

- [ ] **Step 2: Run the test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: testset passes — 3 `@test` assertions.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Test: hash-equivalent mechanisms share fitted_params shape

fitted_params count and kind multiset must match within a
hash bucket. Replaces the test_dedup_csv_replay.jl hash-
equivalent shape assertion with a minimal exemplar pair.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 9: Move canonicalizer source from `identify_rate_equation.jl` → `mechanism_enumeration.jl`

**Files:**
- Modify: `src/identify_rate_equation.jl` (delete lines 128-289 — the canonicalizer block)
- Modify: `src/mechanism_enumeration.jl` (append the same block after the existing `dedup!` definition)

The block being moved (`src/identify_rate_equation.jl:128-289`) contains:
- The `_canonicalize_rate_eq_with_map` docstring and function.
- `_sort_run_factors` and `_factor_sort_key` helpers.
- `_canonical_rate_eq_hash_data` and `_canonical_rate_eq_hash`.

`_project_cached_params` (line 291+) is **not** moved.

- [ ] **Step 1: Capture the exact source to move**

Run: `sed -n '128,289p' src/identify_rate_equation.jl > /tmp/canonicalizer_block.jl`
Expected: the file `/tmp/canonicalizer_block.jl` contains the docstring starting `Build the canonical text + name_map.` through the end of `_canonical_rate_eq_hash`'s function body (final line: `end`).

Quickly verify the block boundaries:
- First line should be: `"""` (start of `_canonicalize_rate_eq_with_map` docstring)
- Last line should be: `end` (close of `_canonical_rate_eq_hash`)
- The line immediately AFTER line 289 in `src/identify_rate_equation.jl` should be the next function — `_project_cached_params` opening docstring `"""`.

- [ ] **Step 2: Delete the block from `identify_rate_equation.jl`**

Edit `src/identify_rate_equation.jl` to remove lines 128-289 (the entire canonicalizer block). The line immediately above (line 127, the closing `end` of `_Stage1Failure`) stays. The line immediately below (line 290, blank or `"""` of `_project_cached_params`) stays.

After deletion, line 128 onwards should pick up at `_project_cached_params`.

- [ ] **Step 3: Append the block to `mechanism_enumeration.jl`**

Find the end of `src/mechanism_enumeration.jl` — currently line 2270 with the closing of `dedup!`. Append (after a blank line):

```
# ─── Rate-equation canonical hash ──────────────────────────────────────
```

Then append the contents of `/tmp/canonicalizer_block.jl` (the captured block).

- [ ] **Step 4: Run the full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes — the dedup block's tests (which call `EnzymeRates._canonical_rate_eq_hash` via qualified name) continue to work from the new location.

If a `UndefVarError: ANNOTATION_SUBSTITUTED not defined` error appears, the include order is wrong. Verify `src/EnzymeRates.jl` includes `rate_eq_derivation.jl` (defines `ANNOTATION_SUBSTITUTED`) before `mechanism_enumeration.jl` (now uses it). The current order is already correct.

If a `UndefVarError: parameters not defined` or similar appears, the canonicalizer references a function defined in `identify_rate_equation.jl` — but actually it only uses `parameters` (from `types.jl`), `rate_equation_string` (from `rate_eq_derivation.jl`), `ANNOTATION_SUBSTITUTED` (from `rate_eq_derivation.jl`), and standard-lib functions. All available in `mechanism_enumeration.jl`'s scope.

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl src/mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Move rate-equation canonicalizer to mechanism_enumeration.jl

_canonicalize_rate_eq_with_map, _canonical_rate_eq_hash_data,
_canonical_rate_eq_hash, _sort_run_factors, _factor_sort_key
moved verbatim from identify_rate_equation.jl. No behavior
change — qualified-name access via EnzymeRates.<fn>
continues to work.

_project_cached_params stays in identify_rate_equation.jl
(beam-search-specific).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 10: Delete the two old canonical-hash testsets in `test_identify_rate_equation.jl`

**Files:**
- Modify: `test/test_identify_rate_equation.jl` (delete lines 492-559 — the two testsets)

- [ ] **Step 1: Verify the deletion boundaries**

Run: `sed -n '490,493p;559,562p' test/test_identify_rate_equation.jl`
Expected:
- Line 490: closing `end` of the testset PRECEDING the two to delete.
- Line 491: blank.
- Line 492: `@testset "canonical rate-equation hash: basic" begin`
- Line 559: closing `end` of the `canonical hash collapses Pattern-A LDH duplicates` testset.
- Line 560: blank.
- Line 561: `@testset "_onesided_permutation_p" begin` (next testset, stays).

- [ ] **Step 2: Delete lines 492-559**

Edit `test/test_identify_rate_equation.jl` and delete the two `@testset` blocks at lines 492-559.

- [ ] **Step 3: Run the full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes. The two deleted testsets' assertions are covered by the new dedup block testsets (determinism, hex shape, distinct, kinetic_group sharing).

- [ ] **Step 4: Commit**

```bash
git add test/test_identify_rate_equation.jl
git commit -m "$(cat <<'EOF'
Drop canonical-hash testsets superseded by dedup block

'canonical rate-equation hash: basic' and 'canonical hash
collapses Pattern-A LDH duplicates' are covered by the
minimal-exemplar tests in the new 'Rate-equation canonical
hash dedup' block of test_mechanism_enumeration.jl.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 11: Delete `test_eq_hash_dedup.jl`, `test_dedup_csv_replay.jl`, update `runtests.jl`

**Files:**
- Delete: `test/test_eq_hash_dedup.jl`
- Delete: `test/test_dedup_csv_replay.jl`
- Modify: `test/runtests.jl` (remove the two `include` lines)

- [ ] **Step 1: Delete the two test files**

```bash
git rm test/test_eq_hash_dedup.jl test/test_dedup_csv_replay.jl
```

- [ ] **Step 2: Remove the two `include` lines from `runtests.jl`**

Edit `test/runtests.jl` and delete these two lines:

```julia
    include("test_eq_hash_dedup.jl")
    include("test_dedup_csv_replay.jl")
```

The final `runtests.jl` should look like:

```julia
using Test
using EnzymeRates
using LinearAlgebra
using Random

include("mechanism_definitions_for_test_enzyme_derivation.jl")
@testset "EnzymeRates.jl" begin
    include("test_accessors.jl")
    include("test_types.jl")
    include("test_dsl.jl")
    include("test_rate_eq_derivation.jl")
    include("test_fitting.jl")
    include("test_mechanism_enumeration.jl")
    include("test_identify_rate_equation.jl")
    include("test_readme_runs.jl")
    include("test_aqua_jet.jl")
end
```

- [ ] **Step 3: Run the full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes. Test-suite runtime drops noticeably because the CSV replay no longer runs.

- [ ] **Step 4: Commit**

```bash
git add test/runtests.jl
git commit -m "$(cat <<'EOF'
Drop test_eq_hash_dedup.jl and test_dedup_csv_replay.jl

CSV-replay test (~22k mechanism-type compiles, >30min runtime)
and cv_results.csv-driven cluster test are replaced by
minimal-exemplar tests in the consolidated dedup block of
test_mechanism_enumeration.jl.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 12: Delete `dedup_investigation/` directory

**Files:**
- Delete: `dedup_investigation/` (whole directory, ~22 MB)

- [ ] **Step 1: Confirm no tests still reference the directory**

Run: `grep -rn dedup_investigation test/ src/`
Expected: no output (no file still references the directory).

- [ ] **Step 2: Remove the directory from git tracking**

```bash
git rm -r dedup_investigation/
```

Note: the directory contains one currently-untracked file (`dedup_investigation/refit_results.csv`). The `git rm -r` removes tracked files; the untracked file needs an additional `rm`:

```bash
rm -rf dedup_investigation/
```

- [ ] **Step 3: Confirm the directory is gone**

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
and refit_diagnostic.jl are no longer referenced by any test.
Findings are preserved in git history.

Frees ~22 MB.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 13: Final verification

**Files:** none modified — verification only.

- [ ] **Step 1: Run the full test suite and measure dedup block timing**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite passes. Aqua + JET clean. The new `@testset "Rate-equation canonical hash dedup"` reports ~10 sub-testsets all passing.

To measure the dedup-block time in isolation, run inside `julia --project`:

```julia
using Test, EnzymeRates
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
@time include("test/test_mechanism_enumeration.jl")
```

The full `test_mechanism_enumeration.jl` file (~4500 lines) includes the new dedup block. Expected: the dedup block itself contributes <30 s of the total time.

- [ ] **Step 2: Verify acceptance criteria from spec**

Run these checks and confirm each:

```bash
# Canonicalizer functions live in mechanism_enumeration.jl now
grep -l "_canonicalize_rate_eq_with_map\|_canonical_rate_eq_hash\b" src/
# Expected: src/mechanism_enumeration.jl  (and src/identify_rate_equation.jl
# may still reference _canonical_rate_eq_hash_data as a caller — that's fine)

# Old test files are gone
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
# Expected: no diff in the exports list
```

- [ ] **Step 3: Confirm no API surface change**

The public API (`parameters`, `fitted_params`, `rate_equation`, `rate_equation_string`, `EnzymeMechanism`, `AllostericEnzymeMechanism`, `identify_rate_equation`, beam-search behavior) is unchanged. The canonicalizer is still accessible via `EnzymeRates._canonical_rate_eq_hash` (unexported, private), now sourced from `mechanism_enumeration.jl`.

If all checks pass, this plan is complete. No final commit — Task 12's commit is the last state-changing one.
