# Required-Regulator Seeding Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Seed `identify_rate_equation`'s beam from fully-regulated mechanisms instead of fitting the partially-regulated shelf beneath them, with optional activator/inhibitor sign designation and a Δ0 competitive-site merge move.

**Architecture:** By default every declared regulator is required. A new `seed_mechanisms` grows `init_mechanisms` into the fully-required seed set by a bounded structural closure and replaces `init_mechanisms` at the beam's seed line. Sign designation, parsed from `::Activator`/`::Inhibitor` in the DSL and carried on `RegulatorMults`, is enforced by one `_filter_by_sign` helper applied in both the seed-build and `expand_mechanisms`. Part B adds a Δ0 `_expand_merge_regulatory_sites` beam move.

**Tech Stack:** Julia; the existing `Mechanism`/`AllostericMechanism` decomposed types and `_expand_*` moves in `src/mechanism_enumeration.jl`; the beam driver in `src/identify_rate_equation.jl`; the DSL in `src/dsl.jl`.

## Global Constraints

- Line length 92, 4-space indent. Match surrounding style.
- Every file starts with two `# ABOUTME:` lines.
- `rate_equation` MUST stay allocation-free and sub-120 ns; the perf gate in `test/test_rate_eq_derivation.jl` must stay green. This work does not touch `rate_equation` derivation, so the gate must remain untouched — if any change would affect it, STOP.
- The parameter-naming chokepoint guard (`test/test_types.jl`) must stay green: no stray `Symbol("K…")`/`k`/`V`/`L` literals outside a name renderer.
- Reuse the existing `_expand_*` moves; do not reimplement enumeration.
- Run `julia --project -e 'using Pkg; Pkg.test()'` green before the final commit of each Part.
- Default behavior for a reaction with NO signs and NO optional lists that declares regulators changes deliberately (seed-build now runs); a reaction with all regulators optional, or none declared, is byte-for-byte unchanged.

**Spec:** `docs/superpowers/specs/2026-07-11-required-regulator-seeding-design.md`.

---

## Part A — required-regulator seeding and sign designation

### Task 1: `RegulatorMults` carries an optional sign

**Files:**
- Modify: `src/types.jl` (`struct RegulatorMults` ~261-282)
- Test: `test/test_types.jl`

**Interfaces:**
- Produces: `RegulatorMults(regulator, allowed_multiplicities, sign::Symbol = :unspecified)`; accessor `sign(rm::RegulatorMults)::Symbol` returning `:activator`, `:inhibitor`, or `:unspecified`.

- [ ] **Step 1: Write failing tests** in `test/test_types.jl`:

```julia
@testset "RegulatorMults sign" begin
    rm0 = RegulatorMults(AllostericRegulator(:X), [2])
    @test EnzymeRates.sign(rm0) == :unspecified          # default
    rmA = RegulatorMults(AllostericRegulator(:X), [2], :activator)
    @test EnzymeRates.sign(rmA) == :activator
    @test rm0 != rmA                                      # sign participates in ==
    @test rm0 == RegulatorMults(AllostericRegulator(:X), [2])
    @test hash(rmA) != hash(rm0)
    @test_throws ErrorException RegulatorMults(AllostericRegulator(:X), [2], :bogus)
end
```

- [ ] **Step 2: Run — expect fail** (`sign` not defined / 3-arg constructor missing).

Run: `julia --project -e 'using Pkg; Pkg.test()'` scoped to the new testset, or run the file directly. Expected: FAIL.

- [ ] **Step 3: Implement.** Add a third field `sign::Symbol` with a default. Validate it is one of `(:activator, :inhibitor, :unspecified)`. Extend `==` and `hash` to include `sign`. Add `sign(rm) = rm.sign`. Keep the existing 2-arg call sites working via the default. `AllostericRegulator`/`CompetitiveInhibitor` are unchanged — the sign lives on the reaction entry, not the ligand.

```julia
struct RegulatorMults
    regulator::Regulator
    allowed_multiplicities::Vector{Int}
    sign::Symbol
    function RegulatorMults(regulator::Regulator,
                            allowed_multiplicities::Vector{Int},
                            sign::Symbol = :unspecified)
        all(m -> m ≥ 1, allowed_multiplicities) ||
            error("RegulatorMults: allowed_multiplicities must all be ≥ 1, " *
                  "got $allowed_multiplicities")
        sign in (:activator, :inhibitor, :unspecified) ||
            error("RegulatorMults: sign must be :activator, :inhibitor, or " *
                  ":unspecified, got $sign")
        new(regulator, sort(allowed_multiplicities), sign)
    end
end
sign(r::RegulatorMults) = r.sign
```
Update `==`/`hash` to fold in `sign` (add `&& a.sign == b.sign` and `hash(r.sign, …)`).

- [ ] **Step 4: Run — expect pass.**

- [ ] **Step 5: Commit** `git add src/types.jl test/test_types.jl && git commit -m "Add optional sign to RegulatorMults"`.

---

### Task 2: DSL parses `::Activator` / `::Inhibitor`

**Files:**
- Modify: `src/dsl.jl` (`_parse_regulator_entries` ~148, `_build_regulators_expr` ~223-250)
- Test: `test/test_dsl.jl` (or wherever `@enzyme_reaction` is tested — grep for `@enzyme_reaction` in `test/`)

**Interfaces:**
- Consumes: `RegulatorMults(reg, mults, sign)` from Task 1.
- Produces: `@enzyme_reaction` accepts `allosteric_regulators: X::Activator, Y::Inhibitor, Z` and stores the sign on each entry's `RegulatorMults`. Bare → `:unspecified`. `::Activator`→`:activator`, `::Inhibitor`→`:inhibitor`.

- [ ] **Step 1: Write failing tests.** Build a reaction with signed allosteric regulators; assert the stored signs; assert bare stays `:unspecified`; assert an unknown tag and a sign on a competitive inhibitor both error.

```julia
@testset "DSL regulator sign" begin
    rxn = @enzyme_reaction begin
        substrates: A[C6H12O6]
        products: B[C6H12O6]
        allosteric_regulators: X::Activator, Y::Inhibitor, Z
        oligomeric_state: 2
    end
    signs = Dict(EnzymeRates.name(EnzymeRates.regulator(rm)) => EnzymeRates.sign(rm)
                 for rm in EnzymeRates.regulators(rxn))
    @test signs[:X] == :activator
    @test signs[:Y] == :inhibitor
    @test signs[:Z] == :unspecified
    @test_throws LoadError @eval @enzyme_reaction begin
        substrates: A[C6H12O6]; products: B[C6H12O6]
        allosteric_regulators: X::Bogus
        oligomeric_state: 2
    end
    @test_throws LoadError @eval @enzyme_reaction begin
        substrates: A[C6H12O6]; products: B[C6H12O6]
        competitive_inhibitors: X::Activator
    end
end
```
(A macro-expansion error surfaces as `LoadError` wrapping the `error(...)`. Confirm the exact exception type against how existing DSL error tests are written — grep `@test_throws` in the DSL test file — and match it.)

- [ ] **Step 2: Run — expect fail.**

- [ ] **Step 3: Implement.** In `_parse_regulator_entries`, before the `Symbol`/`Expr(:call)` branches, peel a `::` annotation: an entry `X::Activator` parses as `Expr(:(::), :X, :Activator)`, and `X(1,2)::Inhibitor` as `Expr(:(::), Expr(:call,:X,1,2), :Inhibitor)`. Extract the tag (`:Activator`→`:activator`, `:Inhibitor`→`:inhibitor`, else error), recurse on the inner expression for name+mults, and thread the sign through. Only `kind === :allosteric` accepts a sign; a sign on a competitive/dead-end entry errors. Extend the parsed tuple to `(name, kind, mults, sign)` and pass `sign` into the `RegulatorMults(...)` constructor call in `_build_regulators_expr`.

- [ ] **Step 4: Run — expect pass.**

- [ ] **Step 5: Commit** `git commit -m "Parse ::Activator/::Inhibitor regulator signs in @enzyme_reaction"`.

---

### Task 3: sign helpers and `_filter_by_sign`

**Files:**
- Modify: `src/mechanism_enumeration.jl` (new helpers near the expansion moves)
- Test: `test/test_mechanism_enumeration.jl`

**Interfaces:**
- Consumes: `sign(rm)`, `regulators(rxn)`.
- Produces:
  - `_regulator_sign(reg_name::Symbol, rxn)::Symbol` — the declared sign, `:unspecified` if none.
  - `_state_respects_sign(reg_name, state::Symbol, sibling_names::Vector{Symbol}, rxn)::Bool`.
  - `_filter_by_sign(mechs::Vector, rxn)::Vector` — keeps only mechanisms whose every regulatory ligand respects its sign given its site siblings.

Sign rules (`state` is the ligand's `allo_state`; `sibling_names` are the other ligands in its site):
- `:unspecified` sign → always true.
- Never the opposite pure state: an `:activator` is never `:OnlyI`; an `:inhibitor` is never `:OnlyA`.
- `:OnlyA`/`:OnlyI` matching the sign, and `:NonequalAI`, are always allowed.
- `:EqualAI` (an antagonist) is allowed unless a sibling is a **same-sign designated** effector — an activator's EqualAI is rejected when a sibling is a designated activator (it would antagonize an activator → opposite observable), and symmetrically for inhibitor.

- [ ] **Step 1: Write failing tests** covering: unspecified passes everything; activator OnlyI rejected, OnlyA/NonequalAI accepted; activator EqualAI accepted with an inhibitor sibling, rejected with an activator sibling; `_filter_by_sign` drops a mechanism containing a sign-violating ligand and keeps a clean one. Build small `AllostericMechanism`s directly (see the existing constructor-based fixtures in the test file) for the `_filter_by_sign` cases.

- [ ] **Step 2: Run — expect fail.**

- [ ] **Step 3: Implement** the three helpers per the rules above. `_regulator_sign` scans `regulators(rxn)` for an `AllostericRegulator` of that name and returns its `sign`. `_filter_by_sign` iterates each mechanism's `regulatory_sites`, and for each ligand checks `_state_respects_sign(name(lig), state, other_names_in_site, rxn)`; a `Mechanism` (no sites) passes trivially.

- [ ] **Step 4: Run — expect pass.**

- [ ] **Step 5: Commit** `git commit -m "Add regulator-sign predicate and _filter_by_sign"`.

---

### Task 4: `seed_mechanisms` seed-build

**Files:**
- Modify: `src/mechanism_enumeration.jl` (new `seed_mechanisms`, near `init_mechanisms` ~1956)
- Test: `test/test_mechanism_enumeration.jl`

**Interfaces:**
- Consumes: `init_mechanisms`, `_expand_to_allosteric`, `_expand_add_allosteric_regulator`, `_expand_add_dead_end_regulator`, `_filter_by_sign`, `regulatory_sites`, `ligands`, `allo_states`, `cat_allo_states`, `bound_metabolite`.
- Produces: `seed_mechanisms(rxn, required_allo::Set{Symbol}, required_comp::Set{Symbol})::Vector{Union{Mechanism,AllostericMechanism}}` — every returned mechanism binds all required regulators, one required allosteric regulator per single-ligand site, cheap states only, signs respected.

**Design.** BFS closure from `init_mechanisms(rxn)`. For each node apply the three structure moves with the full `rxn`, then keep a child iff it is a **valid seed node**:
- no `:NonequalAI` tag anywhere (helper `_has_nonequalai`);
- every regulatory site is single-ligand (enforces one-site-per-regulator);
- every bound allosteric regulator ∈ `required_allo`; every bound competitive inhibitor ∈ `required_comp` (no optional regulator bound);
- passes `_filter_by_sign` (with single-ligand sites, this reduces to "designated regulator in its pure sign state").

Track visited by `hash`; enqueue valid children not seen; retain (do not store intermediates as objects beyond the queue) the valid nodes that bind **all** of `required_allo` and `required_comp` — those are the seeds. Return them deduped.

- [ ] **Step 1: Write failing tests.** Use a small reaction (uni-uni + two allosteric regulators). Assert:
  - undesignated: seed count == `skeletons × 2²` and each seed binds both regulators, each at its own single-ligand site, no `:NonequalAI`;
  - designating both signs yields exactly `skeletons` seeds, each regulator in its sign state;
  - the per-lineage floor: each seed's `length(fitted_params(compile_mechanism(seed)))` equals its base-lineage count `+ n_required + 1` (L).
  Reuse the closure/param-count helpers already proven in the scratchpad profiling if convenient, but assert on a small reaction so the numbers are hand-checkable.

- [ ] **Step 2: Run — expect fail.**

- [ ] **Step 3: Implement** `seed_mechanisms` with the BFS closure and the valid-seed-node predicate above. Keep it allocation-conscious but correctness-first; it runs serially. Add the two-line `# ABOUTME` only if creating a new file — here it is added to the existing file, so no header change.

- [ ] **Step 4: Run — expect pass.**

- [ ] **Step 5: Commit** `git commit -m "Add seed_mechanisms required-regulator seed-build"`.

---

### Task 5: beam integration + sign enforcement in `expand_mechanisms`

**Files:**
- Modify: `src/identify_rate_equation.jl` (`identify_rate_equation` ~199, `_beam_search` ~692 and its base-seed line ~721)
- Modify: `src/mechanism_enumeration.jl` (`expand_mechanisms` ~1918 — apply `_filter_by_sign`)
- Test: `test/test_identify_rate_equation.jl`, `test/test_mechanism_enumeration.jl`

**Interfaces:**
- Consumes: `seed_mechanisms`, `_filter_by_sign`, `regulators(rxn)`, `AllostericRegulator`, `CompetitiveInhibitor`.
- Produces: `identify_rate_equation(prob; optional_allosteric_regulators::Vector{Symbol}=Symbol[], optional_competitive_inhibitors::Vector{Symbol}=Symbol[], …)`.

- [ ] **Step 1: Write failing tests.**
  - `expand_mechanisms` on a mechanism whose reaction designates a regulator drops sign-violating children (e.g. no `:OnlyI` child for a designated activator), and is unchanged (same children) for a signless reaction.
  - A helper-level test of the seed branch: for a regulated reaction with empty optional lists, the computed required sets are the full declared sets; with all regulators optional, both required sets are empty. (Test the required-set computation directly if `_beam_search` is awkward to call; otherwise a small end-to-end `identify_rate_equation` that asserts the base tier is fully-regulated.)

- [ ] **Step 2: Run — expect fail.**

- [ ] **Step 3: Implement.**
  - In `expand_mechanisms`, wrap the produced `result` with `_filter_by_sign(result, rxn)` before the atom-conservation assert. No-op when no signs are declared, so signless reactions are unchanged.
  - Add `optional_allosteric_regulators`/`optional_competitive_inhibitors` keywords to `identify_rate_equation`; thread them into `_beam_search`.
  - In `_beam_search`, compute `required_allo = setdiff(declared allosteric names, optional_allosteric_regulators)` and `required_comp` likewise (declared names from `regulators(prob.reaction)` filtered by ligand type), then branch the base-seed line: `init_mechanisms` when both required sets are empty, else `seed_mechanisms(prob.reaction, required_allo, required_comp)`.

- [ ] **Step 4: Run — expect pass.** Then run the FULL suite: `julia --project -e 'using Pkg; Pkg.test()'`. Fix any regulated-reaction selection tests that shift because the default now seeds fully-regulated (update expectations or add the all-optional keyword to preserve their old behavior, per the spec's flagged default change). Keep the perf gate and chokepoint green.

- [ ] **Step 5: Commit** `git commit -m "Seed identify_rate_equation from required regulators; enforce signs in expand"`. This completes Part A.

---

## Part B — competitive-site enumeration

### Task 6: `_expand_merge_regulatory_sites`

**Files:**
- Modify: `src/mechanism_enumeration.jl` (new move near the other `_expand_*`)
- Test: `test/test_mechanism_enumeration.jl`

**Interfaces:**
- Consumes: `regulatory_sites`, `ligands`, `allo_states`, `multiplicity`, `RegulatorySite`, `AllostericMechanism`, `_filter_by_sign`, `compile_mechanism`, `fitted_params`.
- Produces: `_expand_merge_regulatory_sites(am::AllostericMechanism)::Vector{AllostericMechanism}` and a `_expand_merge_regulatory_sites(::Mechanism) = AllostericMechanism[]` no-op.

**Design.** For each unordered pair of regulatory sites, build one merged site holding both sites' ligands, and enumerate the Δ0-valid `allo_state` assignments for the merged ligands: each ligand keeps its state, OR one ligand is retagged `:EqualAI` (antagonist) — never both. Drop the all-`:EqualAI` result. Run each candidate through `_filter_by_sign` so only sign-respecting antagonist retags survive. The `AllostericMechanism` constructor canonicalizes site order, so `merge(1,2)` and `merge(2,1)` coincide; distinct merge sequences reaching the same partition dedup by `hash` at the beam's seen-set. Reuse the parent's site `multiplicity` for the merged site.

- [ ] **Step 1: Write failing tests.** Build a seed with a single-ligand `:OnlyA` site and a single-ligand `:OnlyI` site (as in the validated scratchpad probe). Assert:
  - the co-binding `{OnlyA,OnlyI}` child and each antagonist `{EqualAI,OnlyI}` / `{OnlyA,EqualAI}` child appear;
  - all children are Δ0 (`length(fitted_params(compile_mechanism(child)))` equals the parent's);
  - the parent and the merged children are pairwise `rate_equation_string`-distinct;
  - with the activator designated, `{EqualAI,OnlyI}` (activator→antagonist-of-inhibitor) survives but `{OnlyA,EqualAI}` (inhibitor→antagonist-of-activator flips the activator's role via its sibling) still respects signs — assert the surviving set matches the sign rules;
  - a `Mechanism` and a single-site `AllostericMechanism` yield `[]`.

- [ ] **Step 2: Run — expect fail.**

- [ ] **Step 3: Implement** the move per the design. Keep the tag enumeration to: both-keep, and each single-ligand→`:EqualAI`, dropping all-`:EqualAI`.

- [ ] **Step 4: Run — expect pass.**

- [ ] **Step 5: Commit** `git commit -m "Add Δ0 _expand_merge_regulatory_sites move"`.

---

### Task 7: wire the merge move into the beam

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_add_expansions_mech!` ~1931-1941)
- Test: `test/test_mechanism_enumeration.jl`

**Interfaces:**
- Consumes: `_expand_merge_regulatory_sites`.

- [ ] **Step 1: Write failing test.** `expand_mechanisms` on a two-site allosteric mechanism includes a one-site merged child; on a `Mechanism` it does not (and is otherwise unchanged).

- [ ] **Step 2: Run — expect fail.**

- [ ] **Step 3: Implement.** Add `append!(result, _expand_merge_regulatory_sites(m))` to `_add_expansions_mech!`. The `_filter_by_sign` wrapper already applied in `expand_mechanisms` (Task 5) covers the merged children too, so no extra filtering here.

- [ ] **Step 4: Run — expect pass.**

- [ ] **Step 5: Commit** `git commit -m "Apply merge move in expand_mechanisms"`.

---

### Task 8: full suite, floor-volume check, version bump

**Files:**
- Modify: `Project.toml` (version bump)
- Test: full suite

- [ ] **Step 1: Run the full suite** `julia --project -e 'using Pkg; Pkg.test()'`. All green, including the `rate_equation` perf gate and the parameter-naming chokepoint. Fix anything red.

- [ ] **Step 2: Floor-volume sanity check.** In a scratch script, build the PFK seeds (fully sign-designated → 359), run `expand_mechanisms` on a handful of them, and confirm the merge-move child count per parent is on the order of `C(k,2)` (≤ ~10 for k=5 sites) — not a blow-up. If a single parent yields an unreasonable count (hundreds), add a per-parent cap on merge children and note it. Record the observed numbers in the commit message; do not commit the scratch script.

- [ ] **Step 3: Bump the version** in `Project.toml` (patch bump).

- [ ] **Step 4: Commit** `git commit -m "Competitive-site merge move complete; vX.Y.Z"`. This completes Part B.

---

## Self-review notes

- Spec coverage: Task 1-2 → sign on the reaction; Task 3 → sign enforcement helper; Task 4 → `seed_mechanisms` (one-site, cheap, gate, floor invariant); Task 5 → the flip, two role-scoped keywords, and beam-wide sign enforcement; Task 6-7 → the Δ0 merge move (co-binding + antagonist); Task 8 → suite + volume + version. Follow-ups in the spec (V-type-optional, strict-NonequalAI, required antagonists) are intentionally not tasks.
- The `_filter_by_sign` helper is defined once (Task 3) and applied in both `seed_mechanisms` (Task 4) and `expand_mechanisms` (Task 5) — DRY.
- Default-behavior change is contained: `_filter_by_sign` is a no-op without signs, and the seed branch falls back to `init_mechanisms` when no regulator is required, so a signless/all-optional reaction is unchanged.
