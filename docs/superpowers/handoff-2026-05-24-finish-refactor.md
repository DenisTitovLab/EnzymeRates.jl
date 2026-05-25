# Handoff: finish the concrete-types refactor — opaque-helper removal + 5 fixture migration for the LOC reclaim

You are taking over the EnzymeRates concrete-types refactor at the cleanup phase. The big architectural work is done; what remains is finishing two deferred deletions that together unlock ~600 LOC of reclaim and complete the move to the new struct family. Read this prompt in full before any action.

## Repo state

- Path: `/home/denis.linux/.julia/dev/EnzymeRates`
- Branch: `refactor-to-concrete-types-instead-of-symbols` (pushed to origin)
- Tip: `5d6168d` (Stage 7e.3 — CLAUDE.md updated)
- Main branch: `main`
- Open PR: **https://github.com/DenisTitovLab/EnzymeRates.jl/pull/40** (do NOT merge yet; this work amends the PR)
- Tags applied: `stage-6beta-complete`, `stage-6-complete`, `stage-7a-complete`, `stage-7c-complete`, `stage-7d-partial-complete`, `stage-7e-complete`, `refactor-complete` (the last is a milestone tag for the current state; you'll move it forward when done)
- Address the user as **Denis** in any user-facing text
- 86 commits between main and HEAD; final src LOC: **9,385** (vs main 7,136, +2,249 lines = +31%)

## Goal

Reclaim the ~600 LOC the previous refactor session couldn't reach, completing the architectural transition to the concrete struct family. Target: **≤ 8,200 src LOC** (the original plan's success criterion).

The work falls into three coupled efforts:

1. **Delete opaque-form helpers in `src/mechanism_enumeration.jl`** (~250 LOC). Requires rewriting the topology backtracker and dead-end enumeration to operate on decomposed `Species([Substrate(:S)], :E)` form rather than opaque `Symbol` form names.

2. **Migrate the 5 remaining opaque-form test fixtures** (~350 LOC reclaim follows from this — see §B below). Requires either re-deriving the analytical rate-equation formulas for the 5-step decomposed forms OR replacing the analytical-comparison tests with synthetic-data-fit tests (the pattern `test_fitting.jl` uses).

3. **Delete the legacy DSL emission path + dual-Sig accessor branches** (~200 LOC). Becomes unblocked once #2 is done — the 5 fixtures are the last consumers of `_mechanism_from_legacy_sig` + the 2-arg `EnzymeMechanism(metabolites, reactions)` constructor + the 12 `_is_new_sig(Sig)` branches in accessors.

Doing all three completes the refactor. Each is independently shippable but #3 depends on #2.

## Authoritative docs (read in this order)

1. `.claude/CLAUDE.md` — project rules. **Rule #1**: "If you want exception to ANY rule, STOP and get explicit permission from Denis first." The CLAUDE.md has been updated to reflect the post-refactor architecture; trust it as the current state-of-truth for the codebase.
2. `docs/superpowers/specs/2026-05-22-concrete-types-refactor-continuation-design.md` — the continuation spec the previous session executed against. Useful for spec §2 (test integrity NON-NEGOTIABLE) + §4 (test-migration rule).
3. `docs/superpowers/plans/2026-05-22-concrete-types-refactor-continuation.md` — the plan the previous session executed. Tasks 7b.5, 7d.2, and the LOC target are the relevant aspirations the previous session DEFERRED.
4. `docs/superpowers/refactor-deleted-tests.md` — log of every §2.1 narrow-exception test deletion in this branch. Append your own entries here for any further deletions.
5. The PR body at https://github.com/DenisTitovLab/EnzymeRates.jl/pull/40 — honest summary of what landed and what didn't, including precise LOC accounting.
6. Auto-memory index: `/home/denis.linux/.claude/projects/-home-denis-linux--julia-dev-EnzymeRates/memory/MEMORY.md` — read the entries listed there. They capture lessons the previous session learned that you should NOT re-discover.

## What's already done

86 commits across 6 stages. The architecture goals are in place:

- One concrete struct family: `EnzymeReaction`, `Mechanism`, `AllostericMechanism`, `Step`, `Species`, `Parameter` (family).
- `EnzymeMechanism{Sig}` + `AllostericEnzymeMechanism{...}` singleton types remain for `@generated` rate-equation derivation; converted via `compile_mechanism(m::Mechanism)`.
- Chokepoint architecture: `name(p::Parameter, m::Mechanism)` (value-context) + `name(::Type{P}, idx::Int)` (type-context) + shared private `_param_symbol` formatter. Enforced by `test/test_chokepoint.jl` AST walker.
- Struct-based canonical hash (`_canonical_rate_eq_hash_data_impl_struct`); the regex pipeline is deleted.
- DSL accepts the decomposed-Species call grammar `E + S => E(S)` (parses to `Species([Substrate(:S)], :E)`).
- `EnzymeReactionLegacy` retired; `IdentifyRateEquationProblem` dispatches on the concrete `EnzymeReaction`.
- 5 zero-alloc `@generated` accessors (`substrates`/`products`/`regulators`/`reactions`/`enzyme_forms`) walking `Sig` at compile time.
- Spec types `MechanismSpec`/`StepSpec`/`AllostericMechanismSpec` renamed to internal `_RawSpec`/`_RawStep`/`_RawAllostericSpec` (leading-underscore-private).

Tests: 26,854 passing in the latest full Pkg.test() run. All 3 compile-budget gates green (`init_mechanisms` 372/750 trace-compile, `rate_equation` first-call 1.57s/2.1s, dispatch-identity check). `test_rate_equation_performance` 0-alloc/<100ns gate green. Integrity check exits 0.

## What's left — section A: opaque-helper removal (~250 LOC)

The 8 helpers and their precise caller inventory (all in `src/mechanism_enumeration.jl`):

| Helper | Callers (excluding defining function) | Role |
|---|---|---|
| `_form_name` | **0** | Synthesize opaque enzyme-form name from bound list — already unused, deletable immediately for a quick ~30 LOC win |
| `_parse_bound` | 1 — closure in `_bound_metabolites_at_forms` at line 950 | Inverse: parse opaque name back to bound-metabolite set |
| `_bound_mets_from_form_name` | 1 — line 1801 | Same as `_parse_bound` but external |
| `_dead_end_form_name` | 6 — dead-end enumeration at lines 956, 1009, 1163, 1245, 1247, 2010 | Synthesize name for "E + dead-end inhibitor bound" form |
| `_is_estar_form` | 1 — inside `_dead_end_form_name` at line 968 | Distinguish E vs Estar conformations by name prefix |
| `_atoms_dict` | 2 — topology backtracker at lines 314, 317 | Look up atom counts per metabolite |
| `_can_pingpong` | 3 — topology backtracker at lines 451, 500, 574 | Atom-balance check for ping-pong residual splits |
| `_subtract_atoms` | 4 — topology lines 278, 384, 454, 577 | Compute residual atom inventory for ping-pong |

**What's blocking removal**: the catalytic-topology backtracker (~lines 200-600 of `src/mechanism_enumeration.jl`) and the dead-end enumeration (~lines 900-2010) build opaque `Symbol` enzyme-form names internally as their working representation. The topology backtracker carries Symbols through `Vector{StepSpec-like-form}` accumulators; the dead-end enumerator carries Symbols through `bound::Dict{Symbol, Set{Symbol}}` lookups; then `_mechanism_from_raw` translates Symbols → decomposed `Species` at the boundary.

**The rewrite path** (the previous implementer called this "option b"):
- Replace the topology backtracker's working representation from `Vector{Tuple{Vector{Symbol}, Vector{Symbol}, Bool, Int}}` (lhs_syms, rhs_syms, is_eq, kg) to `Vector{Step}` with decomposed `from_species::Species`/`to_species::Species` directly.
- Replace the dead-end enumerator's `bound::Dict{Symbol, Set{Symbol}}` with reading `s.from_species.bound` / `s.to_species.bound` from existing Steps.
- `_can_pingpong` / `_subtract_atoms` work on residual atom dictionaries — they can stay as helpers OR move into the `Step` constructor's residual computation; choose whichever fits the rewritten algorithm shape.

This is a substantial algorithm rewrite (~500 LOC modified in `_apply_equivalence_grouping` and `_expand_substrate_product_dead_ends`). Plan it carefully:
- Start with the SMALLEST helper deletion (`_form_name` — 0 callers — pure win, ~5 lines).
- Then `_parse_bound` / `_bound_mets_from_form_name` consolidation (they're near-duplicates).
- Then tackle `_dead_end_form_name` + `_is_estar_form` (the dead-end enumeration cluster).
- Last: `_atoms_dict` + `_can_pingpong` + `_subtract_atoms` (the topology-backtracker cluster).

Verify after EACH helper removal that:
- Full `Pkg.test()` PASSES.
- `bash scripts/check_test_integrity.sh main` exits 0.
- The 3 compile-budget gates stay green (especially trace-compile count, which may shift).
- `test_rate_equation_performance` 0-alloc/<100ns gate stays green.
- The canonical-hash partition test in `test/test_canonical_hash_partition.jl` produces identical class counts (uni_uni=1, bi_bi=23) — this catches subtle behavior shifts in the enumeration that the unit tests miss.

## What's left — section B: 5 opaque-form fixtures

These 5 fixtures use lumped "central complex" Symbol names (`:EABEPQ`, `:EAFP`/`:FBEQ`) that fuse binding + iso into one step. Their analytical rate formulas assume the lumped step count. Splitting the central complex into explicit iso steps changes step count → analytical formula no longer matches.

| Fixture | Location | Lumped state | Steps in fixture | Analytical formula | Migration difficulty |
|---|---|---|---|---|---|
| Segel Ordered Bi Bi | `mechanism_definitions_for_test_enzyme_derivation.jl:211-252` | `:EABEPQ` | 4 | Eq. IX-87, 11-term denominator | Re-derive needed |
| Segel Ping Pong Bi Bi | `mechanism_definitions:309-352` | `:EAFP`, `:FBEQ` | 4 | Eq. IX-140, 8-term denominator | Re-derive needed |
| Segel Ordered Ter Bi | `mechanism_definitions:358-...` | `:EABCEPQ` | 5 | Eq. IX-195 | Re-derive needed |
| Segel Ordered Ter Ter | `mechanism_definitions:...` | `:EABCEPQR` | 6 | Eq. IX-200ish | Re-derive needed |
| `m_manual` stress test | `test_rate_eq_derivation.jl:927-948` | `:EAFP`, `:FBEQ` (mixed) | 16 | **None** — asserts `@test_throws "polynomial terms"` | **Easy** — no formula constraint |

### Migration approaches

Choose per fixture based on what's feasible. The previous session left them opaque because Denis authorized only mechanical migration in Stage 7b.2/7b.3; this handoff explicitly authorizes the re-derivation work.

**Option A — Re-derive analytical formula** (the cleanest, hardest):
- Use King-Altman or determinant-based derivation for the explicit-iso form.
- For Segel Ordered Bi Bi, the 5-step explicit-iso form is also documented in Segel's textbook (Chapter IX, specifically the "with explicit central complex" variant). Cite the textbook ref in the new fixture body.
- The new formula will reference 10 rate constants (k1-k5 + their reverses) instead of 8.
- ~50-100 LOC of math per fixture, NOT mechanical.

**Option B — Replace analytical-comparison test with synthetic-data-fit test** (`test_fitting.jl` pattern):
- Generate synthetic data: `rates = [rate_equation(m, c, true_params) for c in concs]`.
- Fit via CMA-ES: `FittingProblem(m, data; Keq=...)`.
- Assert: recovered params ≈ true_params within tolerance.
- The test runs ~30-60 seconds per fixture (CMA-ES is slow). The suite already accepts this — `test_fitting.jl` does exactly this for ordered_bi_bi.
- Pros: no math required. Cons: stochastic (CMA-ES); may add to the documented "CMA-ES flake" pattern.

**Option C — Change the rate-equation derivation algorithm to handle "fused-iso" steps as a first-class concept**:
- Add a `fused::Bool` flag to `Step` indicating the step combines binding + iso.
- Adapt `_dependent_param_exprs` to fold the iso transition into the binding-step rate constants.
- This would preserve the 4-step analytical form natively.
- Substantial algorithm work; probably the wrong abstraction (the chemistry doesn't actually fuse the iso; the lumping is a modeling shortcut).

**Recommendation**: Option B for `m_manual` (no formula) + Segel Ping Pong (simplest split). Option A for the three Segel Ordered fixtures (textbook references make the math tractable).

**Alternative — delete the fixtures**: spec §2.1 narrow-exception deletion is available. Loss: 5 Segel-textbook reference mechanisms whose analytical formulas validated `rate_equation_string(m)` against known-correct derivations. The canonical-hash partition test + other 30+ analytical-comparison fixtures still cover the derivation pipeline. Denis-authorized if you find Option A/B intractable.

## What's left — section C: legacy emission + dual-Sig branches

Becomes UNBLOCKED once section B is done.

Once no fixture uses opaque-form bare-Symbol step entries, you can delete:

1. **2-arg `EnzymeMechanism(metabolites, reactions)` constructor** in `src/types.jl` (lines around 749-870) — ~100 LOC.
2. **`_mechanism_from_legacy_sig`** in `src/types.jl` (lines around 698-790) — ~90 LOC.
3. **`_is_new_sig` helper + the 12 dual-Sig branches** in `EnzymeMechanism{Sig}` accessors (lines 1264, 1346, 1359, 1372, 1390, 1479, 1493, 1502, 1509, 1520, 1531, 1554, 1606 — find via `grep -n "_is_new_sig" src/types.jl`) — ~150 LOC after collapse.
4. **Narrow `_is_conformation_shape` regex** in `src/dsl.jl` — actually loosened in this branch to accept `:E_c` etc.; the narrowing direction is to reject opaque bound-form names (`:E_S`) at parse time. Once no fixture needs them, raise a clear migration error for any remaining opaque bound-form Symbol entry.

Run the integrity check + full suite after EACH deletion. After all 4: tag `refactor-final-complete`, push, and the LOC target ≤ 8,200 should be hit (the previous session's 9,385 minus ~600 reclaim = ~8,785, plus the dual-Sig collapse savings of ~150 = ~8,635 — under the target).

## Non-negotiables (will fail the refactor if violated)

These carry over from the original spec and are not negotiable:

1. **`rate_equation` perf invariant**: 0 allocations, <100 ns per call for every mechanism in `MECHANISM_TEST_SPECS`. Enforced by `test_rate_equation_performance` in `test/test_rate_eq_derivation.jl`. If a proposed change would force `rate_equation` to allocate or slow down, **STOP and ask Denis**.

2. **Test integrity** (spec §2 + continuation §4): no test deletion, weakening, `@test_skip`, or `@test_broken` without §2.1 log entry. Migration replacement testsets must land in the same commit as deletions (or migration commit first, deletion commit second). `bash scripts/check_test_integrity.sh main` must exit 0 at every commit.

3. **Compile-budget gates** stay green (`test/test_compile_budget.jl`). Current budget: trace-compile 750. Re-baseline DOWN if your changes drop the count; never raise to mask regression.

4. **No `--amend`.** Always create new commits.

5. **No temporal-context comments in code** ("Stage N", "previously", "legacy", "will be"). Documentation in plans/specs/PR descriptions is fine; code comments must be evergreen.

6. **Chokepoint exclusivity**: every `Parameter → Symbol` rendering in `src/` flows through `name(p, m)` / `name(::Type{P}, idx)` / `_param_symbol(...)`. Enforced by `test/test_chokepoint.jl` AST walker. Don't reintroduce direct `Symbol("K…")` literals.

## Workflow conventions

- **Test invocation**: use `julia --project=. -e 'using Pkg; Pkg.test()'` for the full suite (despite the original handoff's OOM warning — confirmed working through 80+ runs in the previous session). The full suite takes ~10-12 minutes.
- **Per-file testing** during iteration:
  ```bash
  julia --project=. -e '
    using Pkg; Pkg.activate(temp=true); Pkg.develop(path=".")
    Pkg.add(["Test","Aqua","JET","OptimizationBBO","OptimizationPyCMA","OrdinaryDiffEqFIRK","Tables","DataFrames","Statistics","Optimization","Random","CSV"])
    using Test, EnzymeRates
    include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
    include("test/<filename>.jl")
  '
  ```
  First run resolves deps (slow); subsequent runs in same Julia session are fast.
- **Common shell mistake**: `bash scripts/check_test_integrity.sh main 2>&1 | tail -N; echo "EXIT: $?"` prints the exit of `echo` (always 0), not the script. Check the script's exit code with `bash scripts/...; echo "EXIT=$?"` (no pipe).
- **Background long-running commands**: `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/out.log 2>&1` with `run_in_background: true` in the Bash tool. The harness notifies when complete.
- **Commit message format**: end with `src delta: -X / +Y net Z, cumulative: ±W`. Compute via `wc -l src/*.jl` before/after. Cumulative is vs main (7136 src LOC).
- **Test count fluctuates ±1** between runs (CMA-ES fits in `test_identify_rate_equation.jl > mechanism recovery` are non-deterministic). Don't over-interpret single-test count shifts. The "1 failed" in `mechanism recovery: Test Failed at test_identify_rate_equation.jl:190` is the known flake; ignore unless it persists across multiple runs.
- **Stage closeout cadence**: previous session paused at end of each Stage for Denis's go/no-go. Continue this pattern — group your work into logical commits, then check in.
- **Subagent dispatch**: use the `superpowers:subagent-driven-development` skill if the task is independent-enough. Spawn fresh subagents per task; provide full context inline (don't make them read large plan files).

## Session-specific learnings worth knowing

These are lessons the previous session paid context to discover. Use them rather than re-deriving:

1. **`_legacy_step_tuple` direction inference**: when `bound(from_species)` and `bound(to_species)` are both empty (opaque Species) AND `bound_metabolite` isn't in either, fall through to "binding direction". For decomposed Species where the released metabolite isn't bound in either (e.g., `E(S) → E + P` releasing chemical product P), compare bound-list sizes: `length(bound(from)) > length(bound(to))` → release direction. Fixed in commit `333af0e`.

2. **Step canonicalization gating**: the Step constructor canonicalizes ONLY when `is_equilibrium == true`. RE binding swaps bound-met to to-side; RE iso swaps to lex-smaller from-side. SS steps preserve user-source direction (analytical formulas reference `:kNf` as source-forward). Per CLAUDE.md "Canonical Step Form". Fixed in commits `c2f77e4` + `5ceff74`.

3. **Accessor allocation fix**: 5 EnzymeMechanism accessors (`substrates`/`products`/`regulators`/`reactions`/`enzyme_forms`) MUST be `@generated` walking `Sig` at compile time. Going through `Mechanism(em)` reconstruction allocates ~3.7 kB/call and breaks the zero-alloc gate in `test/test_accessors.jl`. Fixed in commit `bc1c592`.

4. **`_n_fit_params_estimate(::AllostericMechanism)` is unreliable**: under-counts SS `:NonequalRT` splits by 1 per split. For test assertions about parameter counts, prefer ground-truth `length(fitted_params(compile_mechanism(m)))`. Documented in memory `project-n-fit-params-estimate-undercounts.md`.

5. **`test_compile_budget.jl` warmup-reuse test was replaced** with a `which(...) === which(...)` dispatch-identity check (commit `572080a`). The pre-7d.1 warmup-reuse pattern measured benefit on per-arity specialization of `EnzymeReactionLegacy{S,P,R,N}`; non-parametric `EnzymeReaction` makes per-arity specialization impossible by construction, so the wall-clock gate was uninformative AND the subprocess OOM'd. Don't try to restore the wall-clock pattern.

6. **`_RawSpec` exists**: the previous session's Stage 7d.0 renamed `MechanismSpec`/`StepSpec`/`AllostericMechanismSpec` to `_RawSpec`/`_RawStep`/`_RawAllostericSpec` (leading-underscore-private). Functionally identical, just no longer public surface. When you read the heavy pipeline code, those names are scratch types, not spec types.

## First actions

1. Read the 6 documents listed under "Authoritative docs" (esp. CLAUDE.md, the PR body, the continuation spec §2/§4, and `refactor-deleted-tests.md`).
2. Verify repo state matches this prompt:
   ```bash
   git log --oneline | head -3   # expect 5d6168d at top
   git status                     # expect clean
   wc -l src/*.jl | tail -1       # expect 9385
   git tag | grep -E "stage-|refactor"  # expect 7 tags
   ```
3. Run `bash scripts/check_test_integrity.sh main; echo "EXIT=$?"` — expect EXIT=0.
4. Optional but recommended: run `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/baseline.out 2>&1; echo "EXIT=$?"` to confirm baseline 26,854 PASS (takes ~10 min).
5. Read the 4 memory entries in `/home/denis.linux/.claude/projects/-home-denis-linux--julia-dev-EnzymeRates/memory/`.
6. Ask Denis any clarifying questions BEFORE starting:
   - Which option (A re-derive / B synthetic-fit / C algorithm change) does he prefer for the Segel fixtures?
   - Should `m_manual` migrate or stay?
   - Tolerance for further deletion of the Segel fixtures via §2.1 narrow exception if Option A/B prove intractable?
7. Begin with section A's quick win: delete `_form_name` (0 callers). One small commit to confirm the workflow. Then move to the substantive work per section A's recommended order.

## Suggested commit pacing

Aim for ~8-15 commits total:

- Commit 1: delete `_form_name` (~5 LOC). Quick win + workflow validation.
- Commit 2: consolidate `_parse_bound` / `_bound_mets_from_form_name` (~30 LOC).
- Commit 3: rewrite dead-end enumeration to read decomposed Species; delete `_dead_end_form_name` + `_is_estar_form` (~120 LOC).
- Commit 4: rewrite topology backtracker to use Step values; delete `_atoms_dict` + `_can_pingpong` + `_subtract_atoms` (~100 LOC, biggest single commit).
- Commit 5: full suite + integrity + compile-budget verification + tag `helpers-removed`.
- Commit 6-9 (one per Segel fixture): migrate to decomposed form with re-derived analytical formula OR synthetic-data-fit test. After each: per-file test + commit.
- Commit 10: migrate `m_manual` to decomposed form (no formula concern).
- Commit 11: delete `_mechanism_from_legacy_sig` + 2-arg `EnzymeMechanism(metabolites, reactions)` constructor (~190 LOC).
- Commit 12: collapse the 12 dual-Sig accessor branches in `src/types.jl` (~150 LOC).
- Commit 13: narrow `_is_conformation_shape` regex + raise clear error for opaque bound-form Symbol step entries.
- Commit 14: final dead-code sweep + update CLAUDE.md if anything else shifted.
- Commit 15: tag `refactor-fully-complete`, push.

After all 15: amend the PR description with the new LOC numbers + the section "Future Work" can be removed.

## When you're done

Verify:
- Full `Pkg.test()` PASSES (target 26,800+ tests).
- All 3 compile-budget gates green (trace-compile may shift; re-baseline if needed).
- `test_rate_equation_performance` 0-alloc/<100ns gate green.
- `bash scripts/check_test_integrity.sh main` exits 0.
- `wc -l src/*.jl` ≤ 8,200.
- Move the `refactor-complete` tag forward to the new tip (or create a new `refactor-fully-complete` tag).
- Update the PR description with the new architecture/LOC accounting.
- Report to Denis: "ready for review/merge".

## If you get stuck

The previous session spent significant context on architectural surprises (Step canonicalization bugs, accessor allocation issues, ordering conflicts between Stage 7b and 7d). The recurring pattern was: dispatch implementer → implementer hits architectural blocker → Denis decides → fix → retry.

Recommended: when an implementer reports BLOCKED on something architectural, **don't keep retrying** the same shape of implementer prompt. Bring it back to Denis with specifics. The previous session learned that the third attempt at the same task usually means the design has issues.

The 5 fixtures + algorithm rewrite for helpers are deeper than the typical Stage commit. Budget ~1 long session OR ~2-3 shorter sessions for this work. Don't try to land all 15 commits in one pass without Denis checkpoints.

Good luck. The architecture is sound; what remains is finishing the cleanup.
