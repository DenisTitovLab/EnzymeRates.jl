# Step Refactor Review Follow-ups Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Apply Denis's 2026-09-29 review of the explicit-stoichiometry `Step` branch: reorganize its tests, rename `ligand` back to `bound_metabolite`, reject mixed kinetic groups, fix a malformed fixture, rewrite serialized reproducers as DSL, unify every step parameter name as `k_X_to_Y` / `K_X_to_Y`, then (on a stacked branch) turn documentation comments into docstrings.

**Architecture:** Tasks 1–6 land on branch `step-explicit-stoichiometry` (version 0.8.0, not yet released, so the renames ship in the same breaking release). Task 7 lands on a new branch `internal-docstrings` created from the final head of Task 6.

**Tech Stack:** Julia 1.12, EnzymeRates (`src/`), `Test`, `TestEnv`, Documenter (docs doctests).

**Spec:** `docs/superpowers/specs/2026-09-28-step-explicit-stoichiometry-design.md` (the design these tasks amend; Task 6 updates its §5). Denis's decisions of 2026-09-29 are quoted in each task.

## Global Constraints

- Every mechanism the enumerator emits keeps its rate law and step/group order; Task 6 changes names only (values and meanings of every constant stay the same).
- Unified names (Denis, 2026-09-29): a side of a step is its enzyme form followed by its free metabolites, joined by `_` (a competitive-inhibitor copy is written with `inh`, as in `name(::Species)`). An SS step with canonical orientation X → Y has `k_X_to_Y` (forward) and `k_Y_to_X` (reverse). An RE step has one constant `K_X_to_Y`, the equilibrium constant of X → Y (products over reactants); a pure binding is named in the release direction, so it stays a dissociation constant (`K_ES_to_E_S` = [E][S]/[ES], the value today's `K_S_E` has); every other RE step is named in its canonical direction (`K_ES_to_EP` = [EP]/[ES], today's `Kiso_ES_to_EP`). No `kon_`, `koff_` or `Kiso_` prefixes remain. The allosteric state tag keeps its place right after the prefix (`K_A_ES_to_E_S`, `k_I_E_S_to_ES`). Regulatory-site constants (`K_<tag><ligand>reg`), `Keq`, `E_total` and `L` are unchanged.
- All `Parameter → Symbol` rendering goes through `name(p, m)`; no `Symbol("K…")`/`Symbol("k…")` literal outside a renderer (guard test in `test/test_types.jl`).
- `rate_equation` stays allocation-free and under 120 ns per call for every spec.
- No backward-compatibility shims or aliases (old names, `ligand`, the four-argument `Step`).
- CLAUDE.md: TDD; 92-character lines; match surrounding style; never remove a comment unless provably false; no temporal wording in names, comments or docs; file-level test helpers prefixed `_testhelper_`; enumeration-test rules for `test/test_mechanism_enumeration.jl`.
- One Julia process at a time (7.7 GB RAM, no swap; the VS Code language server process is fine). Focused run: `julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")'`. Full suite (~15 min) in the background: `julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])' > LOG 2>&1 &`, polled in bounded loops.
- Commits end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and `Claude-Session: https://claude.ai/code/session_014ErtWiLZFt2gypmQ95qxFM`; stage files explicitly; never skip hooks.

## Review Focus

1. A pure binding and its release written in either direction get the same `K_<bound form>_to_<free form>_<M>` name, and its value still means a dissociation constant (Task 6 Step 1).
2. A competitive-inhibitor copy of a substrate and the substrate itself bind to the same form without a name collision (`K_EAinh_to_E_Ainh` vs `K_EA_to_E_A`) (Task 6 Step 1).
3. An allosteric mechanism keeps distinct A-state and I-state names for `:NonequalAI` groups and one shared name for `:EqualAI` groups (Task 6 Step 1).
4. A DSL group that mixes RE and SS steps, or a binding with a Theorell–Chance step, is rejected at construction with a message naming the group, not at derivation (Task 3 Step 1).
5. Residual (ping-pong) forms, whose names contain `_res_+A_-P`, still render valid, unique parameter names (Task 6 Step 1).

---

### Task 1: Split `test/test_step_stoichiometry.jl` by subject

**Files:**
- Delete: `test/test_step_stoichiometry.jl`
- Modify: `test/test_types.jl` (receives the `Step` unit tests)
- Modify: `test/test_rate_eq_derivation.jl` (receives the mass-action oracle helpers, the fused/Theorell–Chance derivation cases and the reversal property test)
- Modify: `test/runtests.jl` (drop the include)

**Interfaces:**
- Produces: the helpers `_testhelper_sp`, `_testhelper_consistent_point`, `_testhelper_mass_action_rate`, `_testhelper_check_against_mass_action` and the constants `_testhelper_fused_cases`, `_testhelper_tc_cases`, each defined in exactly one file; later tasks add tests next to them.

- [ ] **Step 1: Move the tests.** Move every testset about `Step` itself — construction, `consumed`/`released`/`ligand`/`is_binding`/`is_iso`, canonical binding orientation, rejections, sort key, signature round trip, `_assert_chemistry_is_iso`, form-pair names, parameter-name collisions, the Tier 2 orientation test — to `test/test_types.jl`, next to the existing `Step`/`Species` testsets. Move the oracle helpers, the fused-step and Theorell–Chance mass-action testsets and the "reversing written steps changes nothing" testset to `test/test_rate_eq_derivation.jl`. Code moves verbatim; only placement changes. A helper used by both files goes in the file that `test/runtests.jl` includes first (`test_types.jl`), or in `test/mechanism_definitions_for_test_enzyme_derivation.jl` if a focused run of either file needs it; say which in the report. Keep the `ER` alias only if the receiving file does not already have one (reuse its existing alias otherwise, converting the moved code).
- [ ] **Step 2: Remove the file and its include** from `test/runtests.jl`.
- [ ] **Step 3: Verify.** Focused runs of `test/test_types.jl` and `test/test_rate_eq_derivation.jl`; the test counts of the moved testsets must equal those in the old file. Full suite: same total as before (54,563).
- [ ] **Step 4: Commit** ("Move the Step tests next to the types and derivation tests").

---

### Task 2: Rename `ligand` to `bound_metabolite`

**Files:** every file under `src/`, `test/`, `docs/src/` that names `ligand(` (12 uses in `src/`, about 74 in `test/` and docs; `grep -rn "ligand(" src test docs/src`). Do not touch `ligands(::RegulatorySite)` or `ligand` fields of `Kreg` (regulatory sites).

- [ ] **Step 1:** Rename the accessor and every call. Docstring: "The metabolite a pure binding step binds (it consumes exactly that metabolite, releases nothing, and `to_species` is `from_species` with it bound); `nothing` for every other step." Rename local variables named after it only where the name would now mislead.
- [ ] **Step 2:** `grep -rn "\bligand(" src test docs/src` prints nothing (regulatory `ligands(` and `p.ligand` excepted). Focused runs of `test/test_types.jl`, `test/test_rate_eq_derivation.jl`, `test/test_mechanism_enumeration.jl`; then the full suite.
- [ ] **Step 3: Commit** ("Name a pure binding's metabolite bound_metabolite").

---

### Task 3: Reject kinetic groups that mix step kinds or RE/SS flags

Denis, 2026-09-29: a Theorell–Chance step cannot share a kinetic group with a binding step. Today such a group dies at derivation with `KeyError: :k_EA_to_EAB`, and a group holding an RE and an SS step builds and derives a rate string that uses parameters it never declares.

**Files:**
- Modify: `src/types.jl` (new `_assert_uniform_groups`, called by the `Mechanism` and `AllostericMechanism` constructors)
- Test: `test/test_types.jl`

**Interfaces:**
- Produces: `_step_kind(s::Step)` and `_assert_uniform_groups(steps::Vector{Vector{Step}})`; Task 6 builds on them.

- [ ] **Step 1: Failing tests** in `test/test_types.jl`, each asserting the error message (pin text such as `"kinetic group"` and the offending step names):

```julia
@testset "kinetic groups hold one kind of step with one flag" begin
    # binding and Theorell–Chance step in one group
    err = try
        @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                (E(A) + B <--> E(A, B), E(A) + B <--> E(Q) + P)
                E(Q) <--> E + Q
            end
        end
        nothing
    catch e
        e isa LoadError ? e.error : e
    end
    @test err isa ErrorException && occursin("kinetic group", err.msg)
    # RE and SS steps in one group
    err2 = try
        @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                (E + S ⇌ E(S), E(P) + S <--> E(P, S))
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        nothing
    catch e
        e isa LoadError ? e.error : e
    end
    @test err2 isa ErrorException && occursin("kinetic group", err2.msg)
    # bindings of two different metabolites in one group
    # (write the case; expect the same rejection)
    # accepted: context-shared bindings of one metabolite, mirrored isomerizations,
    # and every mechanism in MECHANISM_TEST_SPECS (construct each; no error)
end
```

  (The `@enzyme_mechanism` macro builds the `Mechanism` when the expression runs, so the error surfaces there; adapt the `try` shape if it surfaces at expansion.)
- [ ] **Step 2: Run and see them fail** (no error raised, or the derivation-time `KeyError`).
- [ ] **Step 3: Implement.**

```julia
"""
The kind of a step for kinetic grouping: `(:binding, m)` for a pure binding of
`m`, `(:iso,)` for an isomerization, and `(:transformation, consumed, released)`
for every other step. Steps that share a kinetic group share their constants, so
they must be the same kind of reaction.
"""
function _step_kind(s::Step)
    m = bound_metabolite(s)
    m !== nothing && return (:binding, m)
    is_iso(s) && return (:iso,)
    (:transformation, consumed(s), released(s))
end

"""
Error unless every kinetic group holds steps of one kind (`_step_kind`) with one
RE/SS flag: a shared constant means the same reaction type at the same speed.
"""
function _assert_uniform_groups(steps::Vector{Vector{Step}})
    for group in steps
        s1 = first(group)
        for s in group
            _step_kind(s) == _step_kind(s1) && is_equilibrium(s) == is_equilibrium(s1) &&
                continue
            error("Mechanism: a kinetic group holds $(name(from_species(s1))) → " *
                  "$(name(to_species(s1))) and $(name(from_species(s))) → " *
                  "$(name(to_species(s))), which differ in kind or in RE/SS; the " *
                  "steps of a kinetic group share their constants, so they must be " *
                  "the same kind of step with the same flag")
        end
    end
end
```

  Call it in both constructors right after the step directions are canonicalized. Check `_step_kind` equality compares `Metabolite` values (type and name), so a substrate and its inhibitor copy are different kinds.
- [ ] **Step 4: Run.** The new tests pass; focused runs of `test/test_types.jl`, `test/test_dsl.jl`, `test/test_mechanism_enumeration.jl`, `test/test_rate_eq_derivation.jl`; then the full suite. If any existing fixture or enumerated mechanism now fails construction, stop and report it with the group that fails.
- [ ] **Step 5: Commit** ("Reject kinetic groups that mix step kinds or flags").

---

### Task 4: Correct the remaining malformed ping-pong fixture

**Files:** `test/test_mechanism_enumeration.jl` (~1025-1063, the testset asserting "5 dead-end forms → 7 variants" for `_expand_substrate_product_dead_ends`).

The fixture writes `Estar + P ⇌ Estar(A, P)` and `E(A) <--> Estar(A, P)`: A appears on the enzyme without being bound. The intended forms are `Estar(P)` (P bound to the modified enzyme), as already corrected in the other copies of this seed.

- [ ] **Step 1:** Change both lines to `Estar(P)`. Work out by hand, from `_expand_substrate_product_dead_ends`'s rules, which dead-end forms and variants the corrected topology gives; follow the enumeration-test rules (exact count and, where the testset asserts sets, the exact expected sets written out); update the testset's comment and expectations to the corrected counts. Explain the derivation of the new expectation in the report.
- [ ] **Step 2:** Run the focused `test/test_mechanism_enumeration.jl`; the corrected testset passes with the hand-derived expectation (if the function's output disagrees with your hand derivation, stop and report both).
- [ ] **Step 3: Commit** ("Bind P to the modified enzyme in the dead-end fixture").

---

### Task 5: Write the UndefVar reproducers with `@allosteric_mechanism`

**Files:** `test/reference/allosteric_undefvar_reproducers.jl` (three mechanisms stored as serialized `AllostericEnzymeMechanism{…}` type strings; included by `test/test_rate_eq_derivation.jl`).

- [ ] **Step 1:** For each of the three, write the same mechanism with `@allosteric_mechanism` (read `docs/src/deriving/mwc_allostery.md` and `test/test_allosteric_*.jl` for the syntax). Keep the three family comments. Verify each DSL mechanism equals the decoded one (`AllostericMechanism(dsl) == AllostericMechanism(decoded)`) in a scratch run before deleting the strings.
- [ ] **Step 2:** Replace the strings with the DSL forms; the file keeps its `ALLOSTERIC_UNDEFVAR_REPRODUCERS` constant and ABOUTME lines. Run the focused `test/test_rate_eq_derivation.jl`.
- [ ] **Step 3: Commit** ("Write the UndefVar reproducers in the mechanism DSL").

---

### Task 6: One naming rule for every step constant

Denis, 2026-09-29: unify and simplify the names; drop `kon`/`koff`; use the RE-direction rule in Global Constraints.

**Files:**
- Modify: `src/types.jl` (the `name(p, m)` chokepoint: `_render_binding`, `_render_iso`, the `name(::Kd|Kon|Koff|Kiso|Kfor|Krev, m)` methods; the parameter-name guard; `_assert_no_re_ss_duplicate`; the chokepoint comment listing example names)
- Modify: every test and doc that names a step constant (`grep -rn "kon_\|koff_\|Kiso_\|\bK_[A-Za-z]*_E" test docs/src README.md`, plus names such as `K_S_E`, `K_P_E` across `test/` and `docs/src/`)
- Modify: `docs/src/deriving/textbooks.md` (the jldoctest output and the "Meaning of parameter names" table), the other docs pages that name constants in prose
- Modify: `docs/superpowers/specs/2026-09-28-step-explicit-stoichiometry-design.md` §5 and the "Names" bullet of Plan A in `docs/superpowers/specs/2026-09-28-identifiable-enumeration-findings-and-plan.md` (describe the unified rule; Denis, 2026-09-29)

**Interfaces:**
- Consumes: `bound_metabolite` (Task 2), `_step_kind`/`_assert_uniform_groups` (Task 3).

- [ ] **Step 1: Failing tests** in `test/test_types.jl` (a testset "step constants are named by their reaction"):
  - Michaelis–Menten RE seed (`E + S ⇌ E(S)`, `E(S) <--> E(P)`, `E(P) ⇌ E + P`): `parameters(m, Full)` contains `:K_ES_to_E_S`, `:K_EP_to_E_P`, `:k_ES_to_EP`, `:k_EP_to_ES`.
  - The same with SS bindings: `:k_E_S_to_ES`, `:k_ES_to_E_S`, `:k_E_P_to_EP`, `:k_EP_to_E_P`.
  - Theorell–Chance (`E(A) + B <--> E(Q) + P`): `:k_EA_B_to_EQ_P`, `:k_EQ_P_to_EA_B`; an RE isomerization gives `:K_ES_to_EP`; an RE fused release `E(A,B) ⇌ E(Q) + P` gives `:K_EAB_to_EQ_P`.
  - A competitive inhibitor copy of A bound at E next to the catalytic A binding: `:K_EAinh_to_E_Ainh` and `:K_EA_to_E_A` both present.
  - A `:NonequalAI` allosteric binding group gives `:K_A_ES_to_E_S` and `:K_I_ES_to_E_S`; an `:EqualAI` one gives a single `:K_ES_to_E_S`.
  - A ping-pong residual form renders a valid name (e.g. assert the name for the B binding to `E(; residual = A - P)` equals the string you derive from the rule).
  - Writing a binding as its release, or any step backwards, yields the same names.
  - No name in `parameters(m, Full)` starts with `kon_`, `koff_` or `Kiso_` for every mechanism in `MECHANISM_TEST_SPECS`.
- [ ] **Step 2: Run and see them fail.**
- [ ] **Step 3: Implement the renderer.** Replace `_render_binding` and `_render_iso` with one side-based renderer and route every step parameter through it:

```julia
# A step side: its enzyme form followed by its free metabolites, joined by `_`;
# a competitive-inhibitor copy carries `inh`, as in `name(::Species)`.
_met_label(m::Metabolite) =
    m isa CompetitiveInhibitor ? String(name(m)) * "inh" : String(name(m))
_side_label(form::Species, mets) = join([String(name(form)); _met_label.(mets)], "_")

# `prefix` + state tag + "<a>_to_<b>": the constant of the reaction a → b.
_render_reaction(prefix::String, a::String, b::String, state::Symbol) =
    Symbol(prefix, _state_tag(state), a, "_to_", b)

_forward_sides(s::Step) = (_side_label(from_species(s), consumed(s)),
                           _side_label(to_species(s), released(s)))
```

  - `Kon`/`Kfor`: `k_` + forward sides of the group representative; `Koff`/`Krev`: `k_` + the reversed sides.
  - `Kiso`: `K_` + forward sides; `Kd` (a pure binding): `K_` + the reversed sides (release direction).
  - Update the chokepoint comment's example names.
- [ ] **Step 4: Replace the two guards with one.** With every name built from both sides of the representative step, two groups collide exactly when they contain the same reaction. Replace `_assert_unique_parameter_names` and `_assert_no_re_ss_duplicate` with one check that a reaction (the unordered pair of the two side keys, each side keyed by its form and its sorted metabolites) appears in only one step of the whole mechanism; with `_assert_uniform_groups` in place this also covers a reaction written both RE and SS. Update the error-message tests that pinned the old texts ("appears as both rapid-equilibrium and steady-state", "same parameter names") to the new message, keeping each test's mechanism and intent.
- [ ] **Step 5: Update every expectation and doc.** Tests that name constants (use the rule; do not hand-guess — print `parameters(m, Full)` for fixtures you are unsure about); `test/test_allosteric_golden.jl`'s reference file `test/reference/allosteric_golden_reference.txt` if it stores names (regenerate it only through the test's own documented procedure and review the diff: it must be a pure rename); the docs pages; the jldoctest in `docs/src/deriving/textbooks.md` (update the expected output) and its "Meaning of parameter names" table (rows for `K_ES_to_E_S`, `K_EP_to_E_P`, `k_ES_to_EP`, `k_EP_to_ES`, `k_E_S_to_ES`, `k_ES_to_E_S`, `Keq`, with the same meanings and units as today's rows); the spec §5 and the findings doc's Plan A "Names" bullet.
- [ ] **Step 6: Prove the rename changes names only.** A snapshot of the post-A, pre-rename code exists: `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/regress/snap_new_*.jsonl` (fitted names, Reduced and Full strings for all 12,556 exported mechanisms and the specs; enumeration script `regress.jl` and `enum_common.jl` alongside). Write a check script next to it that, for every mechanism: builds the old→new name map (reimplement the pre-rename renderer from `git show <pre-Task-6 commit>:src/types.jl` in the script); asserts the mapped old fitted-name set equals the new fitted-name set exactly; and evaluates the old Reduced string (parse it into a function) and the new `rate_equation` at 2 random consistent parameter/concentration points under the map, asserting relative agreement ≤ 1e-10. Report counts (mechanisms checked, failures = 0 expected). Name order and term order inside strings may differ (fitted names are sorted by name; monomial factors are ordered by symbol name).
- [ ] **Step 7: Verify.** Focused runs; the full suite; the docs doctests (`julia --project=docs -e 'using Documenter, EnzymeRates; doctest(EnzymeRates)'`); `grep -rn "kon_\|koff_\|Kiso_" src test docs/src README.md` prints nothing.
- [ ] **Step 8: Commit** ("Name every step constant by its reaction").

---

### Task 7: Documentation comments become docstrings (branch `internal-docstrings`)

Denis, 2026-09-29: several functions carry their documentation as `#` comments (for example `_canonicalize_step_directions`, `_entry_kind`); look at all the code and correct this.

**Files:** every file in `src/`.

- [ ] **Step 1: Branch.** `git checkout -b internal-docstrings` from the final head of Task 6.
- [ ] **Step 2: Convert.** For every `#` comment block that sits directly above a `function`, one-line method definition, `struct`, `const` or `macro` and documents that definition, turn it into a `"""` docstring attached to the definition, keeping its content (reflow to 92 columns; keep every true sentence). Leave `# ABOUTME:` lines, section banners (`# ─── … ───`), comments inside function bodies, and comments that describe a group of definitions rather than the one below them. A definition with both a docstring and a comment block above it gets one merged docstring. Do not change code.
- [ ] **Step 3: Verify.** `git diff --stat` shows only `src/` files; the full suite passes; the docs build (`julia --project=docs docs/make.jl`) runs without new warnings (`checkdocs = :exports` is unaffected by internal docstrings).
- [ ] **Step 4: Commit** ("Document internal functions with docstrings").
