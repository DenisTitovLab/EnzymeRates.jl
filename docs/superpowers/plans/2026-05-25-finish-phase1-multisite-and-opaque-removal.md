# Finish Phase 1: `::Inh` Multi-Role Binding + Opaque-Fixture Removal

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:test-driven-development per task. **This work is exploratory, not mechanical — do NOT delegate whole tasks to unsupervised subagents** (a prior attempt reverted the architecture to escape a blocker). Drive directly or with tightly-scoped, monitored single-step subagents; escalate architectural blockers rather than working around them (continuation spec §11.4).

**Goal:** Express every remaining mechanism in decomposed grammar so the legacy opaque Sig path can be deleted: build a `G6P::Inh` role-tag DSL feature for metabolites that bind in multiple roles (HK), migrate the ~6 tractable opaque fixtures by rename, delete Theorell-Chance (§2.1, un-representable — confirmed by both refactor attempts), then remove the legacy Sig path and re-enable opaque rejection.

**Architecture (settled — see continuation spec §10–§11):**
- Opaque Species lose step direction in the new Sig (§11.1), so all fixtures must be *decomposed*, not routed through the new Sig as opaque.
- Multi-role metabolites: the DSL tags the inhibitor-role binding `G6P::Inh` → `CompetitiveInhibitor(:G6P)` (real chemical name preserved, so `concs.G6P` drives it). Only `::Inh` is supported (YAGNI); untagged ligands take their declared role. `name(::Species)` renders competitive-inhibitor-bound metabolites with an `inh` marker so the inhibitor form (`:E_G6Pinh`) doesn't collide with the product form (`:E_G6P`). Fitted parameter names stay positional (`K12`).

**Tech stack:** Julia 1.10+, Test/Aqua/JET. `rate_equation` 0-alloc/<100ns invariant; compile-budget gates; chokepoint param-naming.

**Base:** branch `refactor-to-concrete-types-instead-of-symbols`, tip `de41b11` (red by construction — rejection on, opaque fixtures not yet migrated; §11.1).

---

## Conventions (every task)

- **TDD**: failing test → implement → green. **No `--amend`**; stay on the branch.
- **Test integrity** (spec §2/§4): no deletion/weakening without a `docs/superpowers/refactor-deleted-tests.md` §2.1 entry in the same commit. `bash scripts/check_test_integrity.sh main; echo "EXIT=$?"` must be 0 (no pipe — the pipe masks the exit).
- **Full suite**: `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/out.log 2>&1; echo "EXIT=$?"` then read the `Test Summary` line — **`Pkg.test` exit is swallowed by the trailing `echo`; trust the summary line, not the notification's exit code** (this bit us once). Per-file iteration uses the temp-env recipe (resolve deps once):
  ```bash
  julia --project=. -e '
    using Pkg; Pkg.activate(temp=true); Pkg.develop(path=".")
    Pkg.add(["Test","Aqua","JET","OptimizationBBO","OptimizationPyCMA","OrdinaryDiffEqFIRK","Tables","DataFrames","Statistics","Optimization","Random","CSV"])
    using Test, EnzymeRates
    include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")'
  ```
- **Green-at-every-commit tension**: while opaque fixtures remain unmigrated the *full* suite is red (the file won't fully derive). Verify each step with the relevant *per-file* test, commit granularly, and treat the **end of Task 7** as the full-suite green gate. This is a documented, intentional red window on a feature branch (not pushed); do not paper over it.
- **Perf/chokepoint gates** — run `test_rate_equation_performance`, the 3 compile-budget gates, and `test_chokepoint.jl` **per-commit during the red window**, not only at the Task 7 gate. They construct their own small mechanisms and do not need the red fixtures to load, so a `_species_name_from_sig`/Step-canonicalization edit that introduces an allocation or compile-budget regression is caught at its origin commit rather than bisected back from Task 7.
- Commit footer: `src delta: -X / +Y net Z, cumulative: ±W` (`wc -l src/*.jl`; cumulative vs main 7136). End messages with `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>`.

---

## Task 0: Unblock the file (disable rejection scaffold)

Rejection is on at `de41b11`, so `mechanism_definitions` won't even load — blocking all iteration. Temporarily remove the two `_assert_no_opaque_terms(side_terms_per_step)` calls (one in `_parse_plain_mechanism_body`, one in `_parse_allosteric_mechanism_body` in `src/dsl.jl`); **keep the `_assert_no_opaque_terms` function** (re-enabled in Task 8). Also remove the two rejection testsets in `test/test_dsl.jl` (`"opaque bound-form names are rejected"`, `"@allosteric_mechanism rejects opaque bound-form names"`) — they are branch-local (not in `main`), re-added in Task 8.

- [ ] Remove the 2 calls + 2 testsets. Run `test/test_dsl.jl` per-file → loads & passes. (Full suite still red: opaque fixtures crash in derivation — expected, fixed by Tasks 5–6.)
- [ ] Commit: `Task 0: disable opaque rejection scaffold (re-enabled in Task 8)`.

## Task 1: Spike — minimal `::Inh` round-trip (design validation, TDD)

Before touching HK, validate the whole `::Inh` path on a trivial mechanism where a product also competitively inhibits.

- [ ] **Failing test** in `test/test_dsl.jl`:
```julia
@testset "::Inh role tag: product that also competitively inhibits" begin
    m = @enzyme_mechanism begin
        substrates: S
        products:   P
        catalytic_inhibitors: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E(P)
            E(P) <--> E + P
            E + P::Inh <--> E(P::Inh)
        end
    end
    forms = EnzymeRates.enzyme_forms(m)
    @test :E_P in forms        # product-bound form
    @test :E_Pinh in forms     # inhibitor-bound form — DISTINCT from :E_P
    @test :E_P != :E_Pinh
end
```
- [ ] Run → FAIL (`_call_form_term_info` errors on the `P::Inh` arg, or forms collide).
- [ ] Implement Tasks 2–4 until this passes.

## Task 2: Parse `::Inh` on call-form ligands and free metabolites

**File:** `src/dsl.jl`.

- [ ] In `_call_form_term_info` (`dsl.jl:303`), accept a ligand arg that is `Expr(:(::), name::Symbol, :Inh)` in addition to bare `Symbol`. Record the role alongside the name. Extend `_StepSideTerm` (`dsl.jl:371`) with a parallel `bound_roles::Vector{Symbol}` (`:default` or `:inh`) aligned to `bound` (keep `bound::Vector{Symbol}` holding the *real* names, e.g. `:G6P`). Reject any tag other than `:Inh` with: `error("@enzyme_mechanism: unknown role tag ::$tag on $name; only ::Inh is supported.")`.
- [ ] In the free-metabolite term parser (`_step_side_term_info` / `_parse_step_side_terms`, `dsl.jl:277-296`), accept a top-level `Expr(:(::), name::Symbol, :Inh)` free term (e.g. `... + P::Inh <--> ...`), classifying it as the inhibitor-role binder.
- [ ] Sort stability: when building the synthesized form name (Task 3), sort `bound` with its role so `E(P, P::Inh)` is deterministic.
- [ ] Verify the Task-1 parse no longer errors (form may still collide until Task 3).

## Task 3: Role-distinct form naming

**Files:** `src/types.jl` (`_species_name_from_sig` at `types.jl:1429` — **the load-bearing producer**; `name(::Species)` at `types.jl:77`), `src/dsl.jl` (synthesized name in `_call_form_term_info`).

- [ ] Decide the marker: an `inh`-suffixed leaf (`:E_Pinh`, `:E_G6P_G6Pinh`). Apply it in **all THREE** form-name producers so they agree. **`_species_name_from_sig` is primary** — it is the function the `@generated` accessors (`enzyme_forms` `types.jl:1571`, `reactions` `types.jl:1478`, `stoich_matrix`) actually call to name enzyme forms; the King-Altman derivation is keyed on its output. The other two only mirror it.
  - **`_species_name_from_sig` (`types.jl:1429`)** — the bound loop at line 1438 does `push!(parts, String(b[2]))` (bare name), ignoring `b[1]` (the kind tag). Change it to append `inh` when `b[1] === :CompetitiveInhibitor`. **Without this edit, `E(P::Inh)` and `E(P)` both render `:E_P` → the two enzyme states collapse in the King-Altman graph and stoichiometry corrupts** (the exact failure `::Inh` exists to prevent — verified: Task 1's spike asserts `:E_Pinh ∈ enzyme_forms(m)`, and `enzyme_forms` routes through this function, NOT `name(::Species)`).
  - `name(::Species)` (`types.jl:77`): when a bound metabolite `m isa CompetitiveInhibitor`, render its segment as `"$(name(m))inh"`.
  - `_call_form_term_info`'s synthesized `sym` (`src/dsl.jl`): append `inh` to an inhibitor-role ligand's segment.
- [ ] **Ripple check (broadened)**: `regulators:` ALSO maps to `role_of[r] = :CompetitiveInhibitor` (`dsl.jl:578`), so the suffix ripples to EVERY form bound by a `regulators:`/`competitive_inhibitors:` metabolite — not just the lone `:E_I` the earlier draft named (which was wrong: the "Competitive Inhibitor" regulator is `R`, so its form is `:E_R`, and there are `_R`-suffixed forms across fixtures #18, #19, #27, #29, #30, the activator fixture, etc. — ~104 `:E_*` assertions in `test/`). Enumerate affected fixtures by inspecting each mechanism's `role_of` for `:CompetitiveInhibitor` metabolites, NOT by guessing leaf names. Update those assertions to the `inh`-suffixed names IN THE SAME COMMIT (mechanical adaptation, same property, renamed form). Analytical-rate fixtures assert positional params (`K1`, `K12`), not form names, so they are unaffected.
  - **Phase 2 naming-divergence note:** the enumeration synthesizes dead-end forms via `_dead_end_form_name` (`src/mechanism_enumeration.jl`) WITHOUT this `inh` suffix. After this rename, DSL/derivation and enumeration name the same competitive-inhibitor-bound concept differently. Resolve the convention when the enumeration is rewritten (Phase 2) — flagged there, not fixed here.
- [ ] **Concentration check**: confirm `rate_equation` for the Task-1 `m` uses `concs.P` for the `E(P::Inh)` binding (real name preserved). Add to the Task-1 testset:
```julia
    concs = (S=1.0, P=2.0); pv = EnzymeRates.fitted_params(m)
    # rate_equation must reference concs.P for the inhibitor term — no concs.Pinh:
    @test_nowarn rate_equation_string(m)
    @test !occursin("Pinh", split(rate_equation_string(m), "concs")[1])  # Pinh not a concentration var
```
- [ ] Run Task-1 testset → PASS. Run `test/test_dsl.jl`, `test/test_types.jl`, `test/test_accessors.jl` per-file → green.
- [ ] Commit: `Task 1-3: ::Inh role tag — multi-role metabolite binding in decomposed DSL`.

## Task 4: kinetic-group + emission wiring for `::Inh`

**File:** `src/dsl.jl` (`_build_step_expr`, `_build_mechanism_expr`).

- [ ] Ensure `_build_step_expr` emits `EnzymeRates.CompetitiveInhibitor($(QuoteNode(name)))` for an inhibitor-role bound/free metabolite (real name), and `Substrate`/`Product` per declared role otherwise. The `role_of` default path stays for untagged ligands; `::Inh` overrides to CompetitiveInhibitor regardless of `role_of`.
- [ ] Verify a 2-ligand form `E(G6P, G6P::Inh)` emits a bound list `[Product(:G6P), CompetitiveInhibitor(:G6P)]` (or per declared role for the untagged one) and renders `:E_G6P_G6Pinh`.
- [ ] Per-file `test/test_dsl.jl` green. Commit if separate from Task 3.

## Task 5: Migrate HK to decomposed grammar with `::Inh`

**File:** `test/mechanism_definitions_for_test_enzyme_derivation.jl` (the `@allosteric_mechanism` for `name="HK"`).

Rewrite HK's `catalytic_steps` from opaque (`E_Glc`, `E_G6Pi`, `E_G6P_G6Pinh`, …) to decomposed, tagging the inhibitory-site G6P with `::Inh` and the regulatory G6P via the existing `regulatory_site`. The catalytic-site G6P stays plain `G6P` (product). Worked translation (mirror the opaque step graph exactly; reference the struct-throughout branch's HK at `origin/refactor-to-use-structs-throughout:test/...:1779` for the decomposed shape, but use `::Inh` not `G6Pinh`):
```julia
catalytic_steps: begin
    (E + Glucose ⇌ E(Glucose) :: EqualRT,
     E(ATP) + Glucose ⇌ E(Glucose, ATP) :: EqualRT,
     E(G6P::Inh) + Glucose ⇌ E(Glucose, G6P::Inh) :: EqualRT)
    (E + ATP ⇌ E(ATP) :: OnlyR,
     E(Glucose) + ATP ⇌ E(Glucose, ATP) :: OnlyR)
    E(Glucose, ATP) <--> E(G6P, ADP) :: EqualRT
    ((E(G6P, ADP) ⇌ E(ADP) + G6P) :: EqualRT,
     (E(G6P) ⇌ E + G6P) :: EqualRT,
     (E(G6P, G6P::Inh) ⇌ E(G6P::Inh) + G6P) :: EqualRT)
    ((E(G6P, ADP) ⇌ E(G6P) + ADP) :: EqualRT,
     (E(ADP) ⇌ E + ADP) :: EqualRT)
    (E + G6P::Inh ⇌ E(G6P::Inh) :: EqualRT,
     E(Glucose) + G6P::Inh ⇌ E(Glucose, G6P::Inh) :: EqualRT,
     E(G6P) + G6P::Inh ⇌ E(G6P, G6P::Inh) :: EqualRT)
end
```
Keep `allosteric_regulators: G6P::OnlyT, Pi::EqualRT`, `catalytic_inhibitors: G6P`, `regulatory_site(...): ligands: G6P, Pi`. The `analytical_rate_fn` (`hk_rate_analytical`) and all `expected_n_*` (`expected_n_states=10`, etc.) stay UNCHANGED — this is behavior-preserving.

- [ ] Apply. Run `test/test_rate_eq_derivation.jl` per-file, filter for HK: its `Analytical Rate`, `Constraints` (n_states=10), `Performance`, etc. must pass. If `Analytical Rate` fails, STOP — the decomposed step graph diverged from the opaque one; diff form-by-form against the opaque version. Do NOT touch `hk_rate_analytical`.
- [ ] Commit: `Task 5: migrate HK to decomposed grammar via ::Inh role tag`.

## Task 6: Migrate the 6 tractable opaque fixtures + delete Theorell

**File:** `test/mechanism_definitions_for_test_enzyme_derivation.jl`.

Each lumped opaque node → single decomposed node (the spike-validated rename; continuation spec §10.1 / §11.2). Formulas + `expected_n_*` unchanged. Do one fixture per commit; run its analytical test after each.

**Ping-pong fixtures are spike-verified (2026-05-25).** A standalone spike of #10 (Segel Bi Uni Uni Uni PP Ter Bi) — `EABFP → E(A,B)`, `FCEQ → F(C)`, `F` bare — derived 5 states / 5 steps / 9 independent params and matched Segel Eq. IX-228 to `rtol=1e-10` over 20 random-parameter trials via the harness's own `test_analytical_rate`/`test_constraint_counting`. The bound-list-size fallback reconstructs the fused conformation-change-plus-release step `E(A, B) <--> F + P` as a release correctly. So the 4 ping-pong fixtures (#10, #12, #13, #14) **keep their analytical formulas** — do NOT drop them or expand steps (the reference branch `refactor-to-use-structs-throughout` expanded only because it lacked this fallback). The per-fixture analytical test below is a safety net, not an expected failure point; #14 (Hexa, two conformation changes E→F→G) is the same pattern with one more conformation.

- [ ] **Segel Ordered Uni Bi** (`EAEPQ`): `E + A <--> EA` → `E + A <--> E(A)`; `EA <--> EAEPQ` (check exact opaque steps) → decomposed central node `E(A)`/`E(A,?)`; release steps → `E(A) <--> E(Q) + P` etc. Translate by mirroring the opaque step graph; verify analytical rate green.
- [ ] **RE Ordered Bi-Bi** (`EABEPQ`, RE steps): `E(A,B)` single central node; same as Segel Ordered Bi Bi rename but with `⇌` RE steps preserved.
- [ ] **Segel Bi Uni Uni Uni PP Ter Bi** (`EABFP`/`FCEQ`): rename lumped nodes to substrate-side decomposed nodes (`EABFP → E(A,B)` then iso to `F(P)`; `FCEQ → F(C)` then `E(Q)`); preserve ping-pong `F`/`E` conformations.
- [ ] **Segel Bi Uni Uni Bi PP Ter Ter** (`FCEQR`, `ER`): same pattern; `ER → E(R)`.
- [ ] **Segel Bi Bi Uni Uni PP Ter Ter** (`EABFPQ`/`FCER`/`FQ`): same pattern.
- [ ] **Segel Hexa Uni Ping Pong** (`EAFP`/`FBGQ`/`GCER`): three-conformation ping-pong (`E`/`F`/`G`); rename each lumped node to a decomposed node on the appropriate conformation.
- [ ] **Delete Theorell-Chance** (`# 6.` block): bind+release in one step is un-representable (both refactor attempts blocked it). Remove the `let ... end` block; add a §2.1 entry to `docs/superpowers/refactor-deleted-tests.md` (reason: compound bimolecular step `EA + B <--> EQ + P`; `Step` has one `bound_metabolite`; coverage adjacent via Ordered Bi Bi + Ping Pong). Commit with the log; backfill SHA.

After all: `mechanism_definitions` derives with zero opaque forms. **Run the FULL suite — this is the green gate.** Expect all green (no `Cycle N ... not proportional` errors). If any fixture still errors, it has an un-migrated opaque node or a translation error — fix that fixture (don't disable).

## Task 7: Verify no opaque forms remain; full-suite green checkpoint

- [ ] `grep -nE "<-->|⇌" test/mechanism_definitions_for_test_enzyme_derivation.jl | grep -E "[A-Z][A-Za-z0-9]*([A-Z]|_[A-Z])"` → only metabolite atoms / declared names, no opaque enzyme forms. Also grep `test_rate_eq_derivation.jl` (`m_manual` already migrated) and other test files for stray opaque `@enzyme_mechanism`/`@allosteric_mechanism` step entries; migrate any found.
- [ ] Full suite green; `check_test_integrity.sh main` EXIT=0. Commit/tag `phase1-fixtures-decomposed`.

## Task 8: Re-enable opaque rejection

> **Sequencing (roadmap ③ decision, 2026-05-25):** the legacy Sig path deletion that used to live here **moved to after the Phase 2 enumerator rewrite** — see `2026-05-25-finish-refactor-roadmap.md`. Deleting the legacy Sig before the enumerator stops emitting opaque Species risks leaving the multi-substrate `identify_rate_equation` path with neither a working legacy path nor a fixed enumerator. Re-enabling rejection (below) is independent of that — rejection gates the DSL *macros*, while the enumerator builds Species programmatically, so it stays in Phase 1.

Both macros emit decomposed and no opaque fixtures remain (Tasks 5–6), so opaque grammar can be gated.

- [ ] **Re-enable rejection**: re-add the two `_assert_no_opaque_terms(side_terms_per_step)` calls (both macros) and the two rejection testsets removed in Task 0 (strengthened to `@test_throws "opaque bound-form name"`). Per-file `test/test_dsl.jl` green.
- [ ] Full suite + integrity + perf/compile/chokepoint gates green. Commit. Tag `phase1-fixtures-decomposed-rejection-on`.

## Task 9: Phase 1 checkpoint + handoff to Phase 2

> Phase 1 ends here. The legacy Sig path is **still present** (its deletion is sequenced after the Phase 2 enumerator rewrite — roadmap ③). Do NOT assert `_is_new_sig`/`_mechanism_from_legacy_sig` are gone or mark the refactor "done" yet.

- [ ] Full `Pkg.test()` green (read the summary line, not the notification exit); all gates green; `wc -l src/*.jl` recorded for the commit footer.
- [ ] Update CLAUDE.md's DSL grammar section with the `::Inh` role-tag note.
- [ ] Report Phase 1 results to Denis (tests, LOC delta, gates), then proceed to the roadmap's Phase 2 (enumerator Symbol→struct + bi-bi exit gate), which is re-planned in its own brainstorm.

---

## Self-review notes
- **Riskiest task is 3** (role-distinct naming ripples to every competitive-inhibitor form). Task 1's spike de-risks it before HK.
- **Don't trust `Pkg.test` notification exit codes** — read the summary line.
- If Task 5 HK analytical rate fails, the decomposed step graph diverged — diff against the opaque graph form-by-form; the rename must be structure-preserving (spec §10.1). Escalate rather than editing `hk_rate_analytical`.
