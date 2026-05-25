# Finish Refactor — Phase 1: Fixture Rename + Legacy-Path Removal

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Migrate the 5 lumped-central-complex Segel fixtures (and the `m_manual` stress test) from opaque bare-Symbol grammar to decomposed-Species grammar by renaming, then delete the legacy DSL emission path, the 2-arg `EnzymeMechanism(metabolites, reactions)` constructor, `_mechanism_from_legacy_sig`, and the 12 dual-Sig `_is_new_sig` accessor branches — leaving one Sig encoding.

**Architecture:** The lumped central complex (`:EABEPQ`) is renamed to a single decomposed node (`E(A, B)`). A spike proved this preserves the derived `rate_equation_string` byte-for-byte (because `stoich_matrix` treats enzyme forms as opaque-by-name, and `_legacy_step_tuple_from_sig`'s bound-list-size fallback infers release direction for the fused catalytic-release step `E(A,B) <--> E(Q) + P`). Once no fixture uses opaque grammar, the DSL rejects opaque bound-form Symbols and the legacy Sig encoding is deleted; all mechanisms flow through `_sig_of`/`_mechanism_from_sig`.

**Tech Stack:** Julia 1.10+, Test.jl, Aqua.jl, JET.jl. `@generated` accessors must stay 0-alloc/<100ns on the `rate_equation` path.

**Source spec:** `docs/superpowers/specs/2026-05-22-concrete-types-refactor-continuation-design.md` §10 (finishing-phase addendum).

---

## Conventions for every task

These apply to every commit (carried from the continuation spec §8):

- **Branch:** stay on `refactor-to-concrete-types-instead-of-symbols`. This work amends PR #40. **No `--amend`** — always new commits.
- **Test integrity (spec §2 + §4):** never delete, comment, weaken, `@test_skip`, or `@test_broken` a test without a `docs/superpowers/refactor-deleted-tests.md` §2.1 log entry, and the replacement (if any) lands in the same commit. `bash scripts/check_test_integrity.sh main` must exit 0 at every commit (check exit code WITHOUT a pipe: `bash scripts/check_test_integrity.sh main; echo "EXIT=$?"`).
- **Full suite after every commit:** `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/out.log 2>&1` (run in background; ~10-12 min). Test count fluctuates ±1 (CMA-ES non-determinism in `test_identify_rate_equation.jl`); the known `mechanism recovery` flake is ignorable unless persistent.
- **Per-file iteration:** use the temp-env recipe from the handoff, or `julia --project=. -e '...'` for package-only checks. For analytical-fixture checks you need the test deps; use:
  ```bash
  julia --project=. -e '
    using Pkg; Pkg.activate(temp=true); Pkg.develop(path=".")
    Pkg.add(["Test","Aqua","JET","OptimizationBBO","OptimizationPyCMA","OrdinaryDiffEqFIRK","Tables","DataFrames","Statistics","Optimization","Random","CSV"])
    using Test, EnzymeRates
    include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
    include("test/test_rate_eq_derivation.jl")'
  ```
- **Perf gates green every commit:** `test_rate_equation_performance` (0 allocs, <100ns) and the 3 compile-budget gates in `test/test_compile_budget.jl`. Re-baseline trace-compile budgets DOWN if counts drop; never raise to mask a regression.
- **Chokepoint exclusivity:** no direct `Symbol("K…")`/`Symbol("k…")` literals outside `name`/`_param_symbol` bodies (`test/test_chokepoint.jl`).
- **No temporal-context comments** in code ("legacy", "previously", "Stage N", "will be"). Plans/specs may use them; code comments are evergreen. **Note:** when you delete the legacy path, also rename `_legacy_step_tuple_from_sig` → `_step_tuple_from_sig` (Task 8) since "legacy" no longer describes it.
- **Commit message footer:** `src delta: -X / +Y net Z, cumulative: ±W`. Compute Δ via `wc -l src/*.jl` before/after; cumulative vs main (the spec's 7,136 baseline). End commit messages with `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>`.

**Line numbers drift:** all `src/types.jl` / `src/dsl.jl` line numbers are as-of plan-writing (tip `beef41a`). Earlier tasks edit these files, so later tasks' line numbers shift. Locate code by the function signature / grep anchors given in each step, not by absolute line number.

**The rename pattern (used by Tasks 1-2, and Task 4's discovery loop):** a lumped central-complex node containing substrates `X,Y` is renamed to the single decomposed node `E(X, Y)` (substrate-side label). The binding step that *enters* it stays a binding step; the step that *leaves* it releasing a product becomes a fused catalytic-release `E(X, Y) <--> E(Z) + P` — valid because the released product need not be in the source bound list. Single-substrate forms (`EA`) → `E(A)`; product forms (`EQ`) → `E(Q)`; modified-enzyme conformations (`F`) stay bare; regulator-bound forms (`E_R1`) → `E(R1)`, (`EA_R1`) → `E(A, R1)`.

---

## File structure

| File | Tasks | Change |
|------|-------|--------|
| `test/mechanism_definitions_for_test_enzyme_derivation.jl` | 1 | Rename 4 fixtures' opaque step entries to decomposed grammar; formulas unchanged |
| `test/test_rate_eq_derivation.jl` | 2 | Rename `m_manual`'s 16 opaque step entries to decomposed grammar |
| `src/dsl.jl` | 3, 5 | Reject opaque bare-form Symbols at parse time; delete legacy emission branch + its helpers |
| `test/test_dsl.jl` | 3 | Add testset asserting opaque-form `@enzyme_mechanism` raises a clear error |
| `src/types.jl` | 6, 7, 8 | Delete `_mechanism_from_legacy_sig`, collapse `Mechanism(em)`, delete 2-arg constructor, collapse 12 `_is_new_sig` branches, delete dead helpers, rename `_legacy_step_tuple_from_sig` |
| `test/test_types.jl` | 7 | Migrate/delete the 4 `@test_throws` on the 2-arg constructor (§2.1 log) |
| `docs/superpowers/refactor-deleted-tests.md` | 7 | Append §2.1 entries |

---

## Task 1: Rename the 4 analytical Segel fixtures to decomposed grammar

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (4 `@enzyme_mechanism` blocks; the `analytical_rate_fn` formulas and all `expected_n_*` fields stay UNCHANGED)

The analytical-rate tests are the verification — they already pass for the opaque forms and must stay green after the rename (the spike proved the derived rate equation is identical). Do all four in one commit; run the rate-eq test file once at the end.

- [ ] **Step 1: Rename Segel Ordered Bi Bi (around lines 211-220)**

Replace the `steps:` block:
```julia
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> E(Q) + P
                E(Q) <--> E + Q
            end
```
(was: `E + A <--> EA` / `EA + B <--> EABEPQ` / `EABEPQ <--> EQ + P` / `EQ <--> E + Q`)

- [ ] **Step 2: Rename Segel Ping Pong Bi Bi (around lines 309-318)**

```julia
            steps: begin
                E + A <--> E(A)
                E(A) <--> F + P
                F + B <--> F(B)
                F(B) <--> E + Q
            end
```
(was: `E + A <--> EAFP` / `EAFP <--> F + P` / `F + B <--> FBEQ` / `FBEQ <--> E + Q`. `E(A)` is the substrate-side label of the `EA≡FP` central complex; the `E(A) <--> F + P` step fuses the E→F conversion with P release. `F` stays a bare conformation.)

- [ ] **Step 3: Rename Segel Ordered Ter Bi (around lines 362-367)**

```julia
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) + C <--> E(A, B, C)
                E(A, B, C) <--> E(Q) + P
                E(Q) <--> E + Q
            end
```
(was: `EA` / `EAB` / `EABCEPQ` / `EQ`.)

- [ ] **Step 4: Rename Segel Ordered Ter Ter (around lines 425-431)**

```julia
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) + C <--> E(A, B, C)
                E(A, B, C) <--> E(Q, R) + P
                E(Q, R) <--> E(R) + Q
                E(R) <--> E + R
            end
```
(was: `EA` / `EAB` / `EABCEPQR` / `EQR` / `ER`.)

- [ ] **Step 5: Run the rate-equation test file**

Run the per-file recipe (Conventions). Expected: the 4 fixtures' `@testset "Analytical Rate"`, `expected_n_states`/`expected_n_steps`/`expected_n_independent_params`/`expected_n_haldane_constraints` assertions all PASS unchanged. If any analytical-rate test fails, STOP — the rename changed the mechanism shape; re-check the substrate-side labeling against the rename pattern before proceeding.

- [ ] **Step 6: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "$(cat <<'EOF'
Migrate 4 Segel lumped-complex fixtures to decomposed grammar

Rename lumped central-complex nodes (EABEPQ, EAFP/FBEQ, EABCEPQ,
EABCEPQR) to single decomposed nodes (E(A,B), E(A), E(A,B,C), ...).
Spike-verified: derived rate_equation_string is byte-for-byte identical
to the opaque form, so analytical formulas and expected_n_* counts are
unchanged.

src delta: 0 / 0 (test-only)
Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 2: Rename `m_manual` stress test to decomposed grammar

**Files:**
- Modify: `test/test_rate_eq_derivation.jl:927-949` (the `m_manual` `@enzyme_mechanism` in `@testset "Rate equation too large error"`)

`m_manual` has no `analytical_rate_fn`; it asserts `@test_throws "polynomial terms" rate_equation_string(m_manual)`. The rename is 1:1 by symbol and preserves the 11-form/16-step structure, so the "too large" throw still triggers.

- [ ] **Step 1: Replace the `m_manual` definition**

```julia
    m_manual = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: R1
        steps: begin
            E + A <--> E(A)
            E(A) <--> F(P)
            F(P) <--> F + P
            F + B <--> F(B)
            F(B) <--> E(Q)
            E(Q) <--> E + Q
            E + R1 <--> E(R1)
            E(A) + R1 <--> E(A, R1)
            F(P) + R1 <--> F(P, R1)
            F + R1 <--> F(R1)
            F(B) + R1 <--> F(B, R1)
            E(R1) + A <--> E(A, R1)
            E(A, R1) <--> F(P, R1)
            F(P, R1) <--> F(R1) + P
            F(R1) + B <--> F(B, R1)
            F(B, R1) <--> E(R1) + Q
        end
    end
```
(Symbol map: `EA→E(A)`, `EAFP→F(P)`, `FB→F(B)`, `FBEQ→E(Q)`, `E_R1→E(R1)`, `EA_R1→E(A,R1)`, `EAFP_R1→F(P,R1)`, `F_R1→F(R1)`, `FB_R1→F(B,R1)`. The `EA <--> EAFP` and `FB <--> FBEQ` iso steps become `E(A) <--> F(P)` and `F(B) <--> E(Q)` — ligand-changing iso, the same shape as the migrated sequential bi-bi at line ~562.)

- [ ] **Step 2: Run the test**

Run the per-file recipe including `test/test_rate_eq_derivation.jl`. Expected: `@testset "Rate equation too large error"` PASSES — `@test_throws "polynomial terms"` still fires (11 forms / 16 steps preserved). If it does NOT throw, the form count dropped below the threshold; STOP and check that each renamed node is distinct (e.g., `F(P)` ≠ `E(A)`).

- [ ] **Step 3: Commit**

```bash
git add test/test_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
Migrate m_manual stress test to decomposed grammar

1:1 opaque->decomposed symbol rename preserving 11 forms / 16 steps, so
the "polynomial terms too large" throw still triggers.

src delta: 0 / 0 (test-only)
Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 3: DSL rejects opaque bound-form Symbols at parse time

**Files:**
- Modify: `src/dsl.jl` (the `_should_emit_new_grammar`/`_all_bare_terms_compatible` decision around lines 398-612, and the macro emission around lines 586-596)
- Test: `test/test_dsl.jl` (add a testset)

Today an opaque bare term makes `_all_bare_terms_compatible` return `false`, so `_should_emit_new_grammar` returns `false` and the macro emits the legacy `EnzymeMechanism(mets, rxns)` shape. After this task, an opaque bare term raises a clear error and the macro always emits the new `EnzymeMechanism(Mechanism(...))` shape.

- [ ] **Step 1: Write the failing test**

Add to `test/test_dsl.jl`, mirroring the existing rejection-test style at `@testset "@enzyme_reaction rejects bare \`regulators:\` label"` (which uses `@test_throws Exception eval(:(@enzyme_... )))`):
```julia
    @testset "opaque bound-form names are rejected" begin
        @test_throws Exception eval(:(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                ES <--> E + P
            end
        end))
    end
```
(`@enzyme_mechanism` errors are raised at macro-expansion time; `eval(:(...))` surfaces the macro's `error(...)`. Use `Exception` to match the existing pattern rather than assuming a specific wrapper type.)

- [ ] **Step 2: Run it — expect FAIL**

Per-file run of `test/test_dsl.jl`. Expected: FAIL — today the opaque form silently emits the legacy shape (no throw).

- [ ] **Step 3: Add the opaque-term rejection and always emit new grammar**

In `src/dsl.jl`, replace the emission decision in the macro body (lines ~586-596):
```julia
    rxns_expr, side_terms_per_step =
        _parse_steps_block_with_groups(steps_block, declared_mets)

    _assert_no_opaque_terms(side_terms_per_step)
    return _build_mechanism_expr(subs_list, prods_list, regs_list,
                                 role_of, side_terms_per_step)
```

Replace `_should_emit_new_grammar` / `_all_bare_terms_compatible` / `_any_call_form` (lines ~398-418, ~599-612) with the single validator:
```julia
# Raise a clear migration error if any bare-enzyme term is an opaque
# bound-form name. A bare-enzyme term `:X` is acceptable iff `:X` is a
# Call-form head seen in this steps block (`E` in `E(S)`) OR matches the
# single-cap-then-lower conformation shape (`:E`, `:Estar`, `:E_c`).
function _assert_no_opaque_terms(side_terms_per_step)
    call_heads = Set{Symbol}()
    for (_, lhs, rhs) in side_terms_per_step
        for t in (lhs..., rhs...)
            t.kind === :call && push!(call_heads, t.conformation)
        end
    end
    for (_, lhs, rhs) in side_terms_per_step
        for t in (lhs..., rhs...)
            t.kind === :bare_enzyme || continue
            (t.sym in call_heads || _is_conformation_shape(t.sym)) && continue
            error("@enzyme_mechanism: `$(t.sym)` looks like an opaque " *
                  "bound-form name; write it as decomposed call notation, " *
                  "e.g. `E(S)` or `E(A, B)`.")
        end
    end
    nothing
end
```
Keep `_is_conformation_shape` (still used). Delete the legacy `mets_expr`/`:(EnzymeMechanism($mets_expr, $rxns_expr))` lines.

- [ ] **Step 4: Run the test — expect PASS**

Per-file run of `test/test_dsl.jl`. Expected: the new testset PASSES; all existing decomposed-grammar tests still pass.

- [ ] **Step 5: Commit**

```bash
git add src/dsl.jl test/test_dsl.jl
git commit -m "$(cat <<'EOF'
DSL rejects opaque bound-form names; always emit decomposed Mechanism

Opaque bare-enzyme step entries (`:ES`, `:EABEPQ`) now raise a clear
migration error pointing at E(S)/E(A,B) call notation, instead of
falling through to the legacy EnzymeMechanism(mets, rxns) emission.

src delta: -X / +Y net Z, cumulative: ±W
Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 4: Enumerate and migrate any remaining opaque DSL consumers

**Files:**
- Modify: any test file whose `@enzyme_mechanism`/`@allosteric_mechanism` call still uses opaque grammar (surfaced by the suite)

Tasks 1-2 migrated the known opaque fixtures; Task 3 turns every other opaque DSL macro call into a hard error. Run the full suite to get the authoritative list, then migrate each with the rename pattern. (Direct `Mechanism(Species([], :ES), ...)` construction in `test/test_mechanism_enumeration.jl` is NOT affected — Task 3 only rejects at the macro parse level. Do NOT touch those.)

- [ ] **Step 1: Run the full suite, capture macro-expansion errors**

```bash
julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/out.log 2>&1; echo "EXIT=$?"
grep -nE "opaque bound-form name|looks like an opaque" /tmp/out.log
```
Expected: zero or a small number of `@enzyme_mechanism: \`X\` looks like an opaque bound-form name` errors (the test files that still use opaque grammar). If zero, skip to Step 3.

- [ ] **Step 2: Migrate each surfaced call with the rename pattern**

For each error, open the file/line, apply the rename pattern (see Conventions): lumped node `containing X,Y` → `E(X, Y)`; single-ligand → `E(A)`; product form → `E(Q)`; regulator-bound `E_R1` → `E(R1)`. Re-run that test file to confirm green. If a surfaced mechanism is a fixture with an `analytical_rate_fn`, its analytical test must stay green (same as Task 1). If a surfaced test exists ONLY to exercise the legacy/opaque grammar (it asserts the legacy shape, not a mechanism property), convert it to assert the new error (as in Task 3) — that is a same-commit replacement, not a deletion.

- [ ] **Step 3: Run the full suite + integrity check**

```bash
bash scripts/check_test_integrity.sh main; echo "EXIT=$?"
julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/out.log 2>&1; echo "EXIT=$?"; tail -20 /tmp/out.log
```
Expected: integrity EXIT=0; suite passes.

- [ ] **Step 4: Commit (only if Step 2 changed files; otherwise skip)**

```bash
git add -p   # stage only the migrated test files (run git status first)
git commit -m "$(cat <<'EOF'
Migrate remaining opaque @enzyme_mechanism calls to decomposed grammar

src delta: 0 / 0 (test-only)
Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 5: Delete the legacy DSL emission helpers

**Files:**
- Modify: `src/dsl.jl`

After Task 3, `_should_emit_new_grammar`, `_all_bare_terms_compatible`, and `_any_call_form` are deleted (Task 3 already removed their call sites). Confirm no stragglers and remove any dangling helper that referenced the legacy `mets_expr` shape.

- [ ] **Step 1: Verify no remaining references**

```bash
grep -nE "_should_emit_new_grammar|_all_bare_terms_compatible|_any_call_form" src/dsl.jl
```
Expected: no output (Task 3 removed both definitions and call sites). If any remain, delete the dangling definitions.

- [ ] **Step 2: Verify the legacy emission expr is gone**

```bash
grep -nE "EnzymeMechanism\(\\\$mets_expr|mets_expr =" src/dsl.jl
```
Expected: no output.

- [ ] **Step 3: Run per-file `test/test_dsl.jl` + commit if anything changed**

If Step 1/2 required edits:
```bash
git add src/dsl.jl
git commit -m "$(cat <<'EOF'
Remove dead legacy-emission DSL helpers

src delta: -X / +Y net Z, cumulative: ±W
Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```
(If Task 3 already removed everything, fold this verification into Task 3's review and skip the commit.)

## Task 6: Delete `_mechanism_from_legacy_sig` and collapse `Mechanism(em)`

**Files:**
- Modify: `src/types.jl:682-685` (`Mechanism(em)`), `:698-764` (`_mechanism_from_legacy_sig`)

`_mechanism_from_legacy_sig` is reachable only via `Mechanism(em)`'s `_is_new_sig(Sig)` false branch. Nothing produces a legacy-shape Sig anymore (Task 3 removed the only emitter; Task 7 removes the 2-arg constructor — verify ordering by checking callers first).

- [ ] **Step 1: Verify no producer of legacy-shape Sig remains except the 2-arg constructor**

```bash
grep -rnE "_mechanism_from_legacy_sig" src/*.jl
```
Expected: only `src/types.jl:684` (the call) + `:698` (def) + comment mentions. The 2-arg constructor (Task 7) is the last legacy-shape Sig producer; since `Mechanism(em)` dispatches on the Sig *value*, deleting the legacy decode here is safe as long as no live `EnzymeMechanism{legacy_sig}` value is constructed. Because Task 3 already routes the DSL to new-shape Sig, and Task 7 deletes the 2-arg constructor, do Task 6 and Task 7 **together** if you prefer one commit, or Task 7 first then Task 6.

- [ ] **Step 2: Collapse `Mechanism(em)` (lines 682-685)**

```julia
Mechanism(em::EnzymeMechanism{Sig}) where {Sig} = _mechanism_from_sig(Sig)
```
(was the `_is_new_sig(Sig) && return _mechanism_from_sig(Sig); _mechanism_from_legacy_sig(Sig)` two-liner.)

- [ ] **Step 3: Delete `_mechanism_from_legacy_sig` (lines 698-764) and its docstring (lines 687-697)**

Delete the entire function and the docstring block above it.

- [ ] **Step 4: Run per-file `test/test_types.jl` + `test/test_accessors.jl`**

Expected: PASS. If a `MethodError` for `_mechanism_from_legacy_sig` appears, a legacy-Sig value is still being constructed — find it via the suite and ensure Task 7 (2-arg constructor deletion) is done first.

- [ ] **Step 5: Commit** (may be combined with Task 7)

```bash
git add src/types.jl
git commit -m "$(cat <<'EOF'
Delete _mechanism_from_legacy_sig; collapse Mechanism(em) to new Sig

Nothing constructs a legacy-shape Sig after the DSL routes opaque-free
mechanisms through EnzymeMechanism(Mechanism(...)) and the 2-arg
constructor is removed.

src delta: -X / +Y net Z, cumulative: ±W
Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 7: Delete the 2-arg `EnzymeMechanism(metabolites, reactions)` constructor

**Files:**
- Modify: `src/types.jl:766-864` (docstring + constructor)
- Modify: `test/test_types.jl` (4 `@test_throws` at lines 104, 113, 355, 364)
- Modify: `docs/superpowers/refactor-deleted-tests.md` (§2.1 entries)

The 4 test_types.jl callers are `@test_throws ErrorException EnzymeMechanism(((:S,),(:P,),()), bad_rxns)` — they validate the constructor's input rejection. Once the constructor is gone they have no referent (§2.1 narrow exception). The surviving validation surface is covered by `test/test_types.jl @testset "EnzymeMechanism error cases"` (line ~391) and the `AllostericEnzymeMechanism constructor validators` (line ~482), plus the `Mechanism` constructor's own validation.

- [ ] **Step 1: Verify the validators it calls aren't shared with the new path**

```bash
grep -rnE "_validate_kinetic_groups|_validate_enzyme_connectivity" src/*.jl
```
If `_validate_enzyme_connectivity` / `_validate_kinetic_groups` are called ONLY from the 2-arg constructor body (lines 774-864), they become dead and are deleted with it. If they have other callers (e.g., the `Mechanism` constructor or `EnzymeMechanism(::Mechanism)`), KEEP them. Record which.

- [ ] **Step 2: Delete the 4 `@test_throws` tests and add §2.1 log entries**

Remove the 4 lines (104, 113, 355, 364) and any now-empty `let`/`@testset` scaffolding around them. Append to `docs/superpowers/refactor-deleted-tests.md`:
```markdown
## Phase 1 — 2-arg EnzymeMechanism constructor deletion (commit TBD-after-commit)

### test_types.jl — 4 `@test_throws` on EnzymeMechanism((mets), (rxns))
- Lines (pre-deletion): 104, 113, 355, 364
- Entity deleted: the 2-arg `EnzymeMechanism(metabolites, reactions)` constructor.
- Replacement: NONE EQUIVALENT — these asserted the constructor rejects
  malformed (mets, rxns) tuple input; the constructor no longer exists.
- Surviving validation coverage:
  - test_types.jl `@testset "EnzymeMechanism error cases"` (~line 391)
  - test_types.jl `@testset "AllostericEnzymeMechanism constructor validators"` (~line 482)
  - the `Mechanism` constructor's structural validation (exercised across
    test_mechanism_enumeration.jl and test_types.jl).
```

- [ ] **Step 3: Delete the 2-arg constructor (lines 766-864) + any now-dead validators from Step 1**

Delete the docstring (766-773) and the `function EnzymeMechanism(mets::Tuple{...}, rxns::Tuple)` body (774-864). Delete `_validate_kinetic_groups`/`_validate_enzyme_connectivity` ONLY if Step 1 found no other callers.

- [ ] **Step 4: Run full suite + integrity check**

```bash
bash scripts/check_test_integrity.sh main; echo "EXIT=$?"
julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/out.log 2>&1; echo "EXIT=$?"; tail -20 /tmp/out.log
```
Expected: integrity EXIT=0 (the §2.1 log covers the 4 deletions); suite passes.

- [ ] **Step 5: Commit, then backfill the log SHA**

```bash
git add src/types.jl test/test_types.jl docs/superpowers/refactor-deleted-tests.md
git commit -m "$(cat <<'EOF'
Delete 2-arg EnzymeMechanism(metabolites, reactions) constructor

Last producer of the legacy-shape Sig. The 4 @test_throws covering its
input validation are removed per spec §2.1 (entity gone; surviving
validation covered by EnzymeMechanism error-case + constructor-validator
testsets); logged in refactor-deleted-tests.md.

src delta: -X / +Y net Z, cumulative: ±W
Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
git rev-parse HEAD   # replace TBD-after-commit in the log, then:
git add docs/superpowers/refactor-deleted-tests.md
git commit -m "Phase 1: backfill deleted-tests log commit SHA

src delta: 0 / 0 (docs only)
Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

## Task 8: Collapse the 12 `_is_new_sig` accessor branches and delete dead helpers

**Files:**
- Modify: `src/types.jl` — accessors at lines 1346, 1359, 1372, 1390, 1479, 1493, 1502, 1509, 1520, 1531, 1554, 1606; `_is_new_sig` (1264); dead `_legacy_step_tuple` (1320) + `_species_sym` (1283); rename `_legacy_step_tuple_from_sig` (1455)

Every `EnzymeMechanism{Sig}` is now new-shape (Tasks 3, 7). Each `_is_new_sig(Sig)` is statically `true`, so the false branch is dead. Collapse each accessor to its new-shape body. `_is_new_sig`, `_legacy_step_tuple`, and `_species_sym` then have no callers. **`_legacy_step_tuple_from_sig` (1455) STAYS** — it's the new-shape `reactions` accessor's helper (line 1483); rename it to drop the misleading "legacy".

- [ ] **Step 1: Confirm `_legacy_step_tuple` (1320) and `_species_sym` (1283) are dead; `_from_sig` stays**

```bash
grep -nE "\b_legacy_step_tuple\b" src/types.jl | grep -v "_from_sig"   # def 1320 + docstring only
grep -nE "\b_species_sym\b" src/types.jl                               # def 1283 + uses inside 1320 only
grep -nE "_legacy_step_tuple_from_sig" src/types.jl                    # def 1455 + USE at 1483 -> keep
```
Expected: `_legacy_step_tuple`/`_species_sym` used only inside each other → safe to delete; `_legacy_step_tuple_from_sig` has a live caller → keep.

- [ ] **Step 2: Collapse each `_is_new_sig` accessor to its new-shape body**

For each of the 12 accessors, delete the `if _is_new_sig(Sig)` guard and the entire `else`/fallthrough legacy body, keeping only the new-shape body. Example — `reactions` (1478-1489) becomes:
```julia
@generated function reactions(::EnzymeMechanism{Sig}) where {Sig}
    tuples = Any[]
    for (g, group) in enumerate(Sig[2])
        for step_sig in group
            push!(tuples, _step_tuple_from_sig(step_sig, g))
        end
    end
    return Tuple(tuples)
end
```
`equilibrium_steps` (1492-1498) becomes `Tuple(step[4] for group in Sig[2] for step in group)`. `n_steps` (1501-1505) becomes `sum(length(group) for group in Sig[2]; init=0)`. `kinetic_groups` (1519-1525) becomes `Tuple(1:length(Sig[2]))`. Apply the same pattern to `substrates` (1345), `products` (1358), `regulators` (1371), `metabolites` (1389), `kinetic_group` (1508), `steps_in_group` (1528), `enzyme_forms` (1553), `stoich_matrix` (1605). Keep each docstring.

- [ ] **Step 3: Delete `_is_new_sig` (1261-1280), `_legacy_step_tuple` (1297-1342), `_species_sym` (1282-1283); rename `_legacy_step_tuple_from_sig` → `_step_tuple_from_sig`**

Delete the three dead definitions + their docstrings. Rename `_legacy_step_tuple_from_sig` to `_step_tuple_from_sig` at its definition (1455) and its caller (1483), and update its docstring to drop the `_legacy_step_tuple` reference (use: "Build a (lhs, rhs, is_eq, g) tuple from a Sig step; infers binding/release direction from the bound-sig lists, falling back to bound-list size when the metabolite is bound in neither species (fused catalytic-release).").

- [ ] **Step 4: Verify nothing references the deleted names**

```bash
grep -rnE "\b_is_new_sig\b|\b_legacy_step_tuple\b|\b_species_sym\b" src/*.jl
```
Expected: no output.

- [ ] **Step 5: Run full suite + integrity + perf gates**

```bash
bash scripts/check_test_integrity.sh main; echo "EXIT=$?"
julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/out.log 2>&1; echo "EXIT=$?"; tail -25 /tmp/out.log
```
Expected: integrity EXIT=0; suite passes; `test_rate_equation_performance` 0-alloc/<100ns green; compile-budget gates green; `test/test_accessors.jl` zero-alloc accessor tests green (this is the gate the accessor collapse most affects).

- [ ] **Step 6: Commit**

```bash
git add src/types.jl
git commit -m "$(cat <<'EOF'
Collapse dual-Sig accessor branches to single new-shape encoding

Every EnzymeMechanism{Sig} is now new-shape, so the 12 _is_new_sig
guards are statically true. Collapse each accessor to its new-shape
body; delete _is_new_sig, the dead _legacy_step_tuple + _species_sym;
rename _legacy_step_tuple_from_sig -> _step_tuple_from_sig.

src delta: -X / +Y net Z, cumulative: ±W
Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 9: Final verification, LOC accounting, tag

**Files:** none (verification + tag)

- [ ] **Step 1: Full green run**

```bash
bash scripts/check_test_integrity.sh main; echo "EXIT=$?"
julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/final.log 2>&1; echo "EXIT=$?"; tail -30 /tmp/final.log
```
Expected: integrity EXIT=0; all tests pass (~26,800+, ±1 CMA-ES flake); all 3 compile-budget gates green; `test_rate_equation_performance` 0-alloc/<100ns green.

- [ ] **Step 2: Confirm the legacy surface is gone**

```bash
grep -rnE "_is_new_sig|_mechanism_from_legacy_sig|_should_emit_new_grammar|EnzymeMechanism\(mets" src/*.jl
```
Expected: no output.

- [ ] **Step 3: Record code-LOC delta (Denis's metric)**

```bash
wc -l src/*.jl | tail -1
```
Note the raw delta in the final commit footer; this Phase removes the legacy decode + dual-Sig branches (expect ~250-350 raw lines down vs Task 1 start).

- [ ] **Step 4: Tag**

```bash
git tag phase1-legacy-removed
```

- [ ] **Step 5: Report to Denis**

Summarize: 5 fixtures + m_manual migrated by rename (formulas unchanged), opaque DSL grammar rejected, legacy Sig encoding + 2-arg constructor + dual-Sig branches deleted, all gates green, LOC delta. Note that Phase 2 (enumeration topology-backtracker Symbol→struct rewrite + parse-back helper removal) remains and needs its own plan.

---

## Phase 2 (separate plan, not in scope here)

Phase 2 rewrites the catalytic-topology backtracker and dead-end enumeration in `src/mechanism_enumeration.jl` to operate on `Step`/decomposed `Species` instead of opaque `Symbol` working-representation, deleting the parse-back helpers (`_parse_bound`, `_bound_mets_from_form_name`, `_dead_end_form_name`, `_is_estar_form`) and the topology atom helpers (`_atoms_dict`, `_can_pingpong`, `_subtract_atoms`). It is the high-risk part (validated only by exact topology counts: bi-bi=11, ter-ter=283, pyruvate carboxylase=312, pyruvate dehydrogenase=334) and gets its own plan written after Phase 1 lands — starting by confirming whether `init_mechanisms` currently emits opaque or decomposed Species, which sets the rewrite's shape.
```
