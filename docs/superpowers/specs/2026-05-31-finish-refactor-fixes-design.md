# Finish-refactor fixes — design

**Date:** 2026-05-31
**Branch:** `refactor-to-concrete-types-instead-of-symbols`
**Goal:** Close the remaining gaps on the concrete-types refactor so the branch is
PR-ready. Combines (a) the prior architecture review's confirmed P1/P2/P3 items, (b) a
second agent's five branch findings (validation, regulators, multiplicity, README,
stale comments), and (c) Denis's two requests: rationalize the three migration-era test
files into ~1:1 src↔test correspondence, and simplify the Sig machinery.

All work lands in **one PR**. The full suite currently passes (27073 pass / 1 broken /
0 fail); the one `@test_broken` is a deliberately-deferred item (Case-B `parameters(Full)`
non-injectivity, tracked by two existing specs) and is **out of scope here**.

---

## Scope decisions (settled with Denis)

| # | Decision |
|---|---|
| Multiplicity (#3) | **Enumerate** all `allowed_catalytic_multiplicities`. Update CLAUDE.md's "not enumerated" line. |
| README (#4) | **Switch the example to `allosteric_regulators: A`** so the "recovers Section 2" narrative stays true. |
| Validation (#1) | **Restore all of `main`'s checks** on the struct constructors, incl. atom mass-balance. |
| Atom-balance × `@enzyme_mechanism` placeholders | **Balanced placeholders.** Keep ALL checks (incl. balance) on the `EnzymeReaction` ctor; fix `@enzyme_mechanism`'s `_build_mechanism_expr` to emit *balanced* placeholder atoms (each substrate gets `:C => n_prods`, each product gets `:C => n_subs`, so per-element totals are equal). Discovered during A2: the DSL assigns `:C => 1` to every substrate/product, so multi-product mechanisms (uni-bi etc.) fail the restored balance check. Placeholder atom *values* are never read for hand-built `@enzyme_mechanism` (the enumerator's `_atoms_dict` topology path only runs on real `@enzyme_reaction` reactions), but they ARE encoded into `Sig` type identity — so the full suite (esp. compile-budget goldens) must stay green after the change. See [[project-atom-balance-vs-enzyme-mechanism-placeholders]]. |
| Sig (Q2) | **Keep the explicit per-type methods** (simple > clever; decode is compile-time-free). Take only the free wins: collapse the 4 metabolite encoders into 1, and fix the `metabolites()` hand-decode bypass. No full reflection rewrite, no Sig layout change. |
| PR scope | All four clusters (A bugs, B dedup/cleanup, C test-reorg, D doc sweep) in this PR. |
| Wegscheider map (B10) | **Skip** the move — it has callers in both files, so relocating only shifts the back-edge. |
| Public surface (P1.4) | **Soften CLAUDE.md wording** ("accessed as `EnzymeRates.X`"); do not expand the 18-name export set. |

---

## Cluster A — Bug fixes (PR blockers)

### A1 + A2 — Name-identity validation on `EnzymeReaction` (one root cause)

`origin/main`'s `EnzymeReaction(subs, prods, regs)` validated mandatory atoms, positive
atom counts, non-empty substrate/product sets, duplicate names, and atom mass-balance
(`main:src/types.jl:21-75`). The refactor's concrete `EnzymeReaction` ctor
(`types.jl:271-282`) validates only `multiplicities ≥ 1`; `ReactantAtoms`
(`types.jl:224-231`) validates nothing. Restore parity **and** add cross-category
name-uniqueness, which also fixes finding #2's regulator/substrate name collision.

**`ReactantAtoms` inner constructor — add:**
- atoms tuple non-empty (atoms are mandatory).
- each atom count is a positive `Integer` and not `Bool` (`c isa Integer && !(c isa Bool)
  && c > 0`).

**`EnzymeReaction` inner constructor — add (keeping the existing `multiplicities ≥ 1`):**
- substrates non-empty; products non-empty.
- no duplicate substrate names; no duplicate product names; no duplicate regulator names.
- **atom mass-balance:** per-element sum over substrate `ReactantAtoms` equals the sum
  over product `ReactantAtoms`, element by element. Reuse the existing atom payload on
  `ReactantAtoms` (no new `_sum_atoms` helper on `Tuple` needed — sum directly over the
  `atoms` vectors of `substrates(r)` / `products(r)`).
- ~~no name collision across categories~~ — **DROPPED during implementation.** This check
  was proposed here as a structural fix for finding #2, but it directly contradicts the
  documented `::Inh` role-tag feature: a metabolite may legitimately bind as a
  `CompetitiveInhibitor` under its real name (e.g. G6P is both a product AND a competitive
  inhibitor in the Hexokinase fixture, driven by one `concs.G6P`). `main` never had this
  check (it validated only *within*-category duplicates). So `EnzymeReaction` restores
  exactly `main`'s four checks: non-empty subs/prods, intra-category duplicate names, atom
  balance — NOT cross-category uniqueness.

Finding #2's regulator-duplication is fixed entirely by A3 (idempotent
`_add_competitive_inhibitor`), not by a cross-category check. **A2 and A3 are
interdependent and must land in one commit:** restoring the `EnzymeReaction` duplicate-
regulator check (A2) exposes that `expand_mechanisms` → `_add_competitive_inhibitor`
re-adds an already-declared dead-end regulator, which the hardened ctor then rejects. A3's
idempotency guard is what keeps enumeration working. (On this branch `EnzymeReaction(...)`
is an internal reconstruction chokepoint — called by Sig round-trip,
`_drop_unbound_regulators` on every `compile_mechanism`, and the enumeration moves — unlike
`main` where it was user-facing only. That is why ctor validation has broad reach and why
two deliberately-unbalanced toy test fixtures had to be balanced.)

**`_add_competitive_inhibitor` (`mechanism_enumeration.jl:1230`) — make idempotent:** if
`reg_name` already names a regulator on `rxn`, return `rxn` unchanged rather than pushing a
duplicate `RegulatorMults` (which the hardened `EnzymeReaction` ctor would now reject). This
fixes finding #2's "dead-end :I expansion produced two :I regulators."

**Implementation note (TDD):** add each check test-first, then run the suite and repair any
fixtures that relied on the absence of validation. **If a fixture is a deliberately
unbalanced toy reaction, STOP and ask Denis** rather than weakening the balance check.

### A3 — Enumerate catalytic multiplicity (#3)

`_expand_to_allosteric` (`mechanism_enumeration.jl:1445-1462`) currently does
`cn = only(allowed_catalytic_multiplicities(rxn))`, throwing `ArgumentError` for any
multi-valued list. Change to iterate:

```julia
for cn in allowed_catalytic_multiplicities(rxn)
    # build the base :EqualAI variant + one :OnlyA-per-group variant at this cn
end
```

Single-valued reactions (every current fixture, incl. `oligomeric_state: N` → `[N]`)
produce byte-identical output. Only genuinely multi-valued reactions fan out into more
allosteric candidates. Update CLAUDE.md's "oligomeric_state … not enumerated" line to
reflect that `allowed_catalytic_multiplicities` IS enumerated.

### A4 — README allosteric-recovery truthfulness (#4)

`_expand_add_allosteric_regulator` (`mechanism_enumeration.jl:1513-1518`) only promotes
declared `AllostericRegulator` entries, so a `competitive_inhibitors: A` reaction yields
zero allosteric-regulator expansions. README §"Recover the mechanism" (`README.md:138-150`)
claims the search recovers the Section-2 MWC model from `competitive_inhibitors: A` — which
it cannot.

**Fix:** change the README reaction to `allosteric_regulators: A` (valid DSL) so the
"recovers the mechanism we used to generate the data" narrative is true. Adjust the
surrounding prose accordingly (drop the "enumerates dead-end-inhibitor and allosteric
variants" framing; the example now declares allostery directly). `test_readme_runs.jl`
must still pass.

---

## Cluster B — Dedup / cleanup (behavior-preserving)

| ID | Where | Change |
|----|-------|--------|
| **B1** | `types.jl:1055-1075` | Rewrite `@generated metabolites(::EnzymeMechanism{Sig})` to `m = Mechanism(EnzymeMechanism{Sig}())` and derive substrate/product/regulator names from the concrete reaction. Removes the only `Sig[…]` hand-decode in `src/`. Stays `@generated` (hot-path type-stability preserved). |
| **B2** | `types.jl:387-391` + `449-453` | Extract `_canonicalize_iso_groups(reaction, groups::Vector{Vector{Step}})`; call from both `Mechanism` and `AllostericMechanism` constructors. Depends only on reaction-derived sets + flattened binding steps — mechanical. |
| **B3** | `rate_eq_derivation.jl:1016` → `types.jl`; `types.jl:1376/1398/1433`; `rate_eq_derivation.jl:1049/1088` | Move `_emit_cat_params_for_rep` to `types.jl` beside the Parameter family. Route the three `types.jl` walkers (`_onlyA_parameters_for_sym`, `_all_params_for_sym`, `_enumerate_parameters_full`) through it. Give the two rename-pair builders (`_I_rename_parameters`, `_A_rename_parameters`) a shared small variant (zip `:A`/`:I`). Make the helper's "centralizing" docstring true. |
| **B4** | `types.jl:556-559` | Collapse the 4 single-name metabolite encoders into one: `_to_sig(m::Metabolite) = (nameof(typeof(m)), name(m))`. No Sig layout change (the tag Symbol is identical to today). The free half of the Sig decision. |
| **B5** | `dsl.jl:563-580` + `1171-1180` | Extract `_reject_opaque_bound_forms(side_terms_per_step, macro_name)`; call from both macros. Fix the `@allosteric_mechanism` path's error (`dsl.jl:1175`) to name `@allosteric_mechanism`, not `@enzyme_mechanism`. |
| **B6** | `test_chokepoint.jl:50` | Broaden the renderer-body classifier to match all 10 Parameter subtype names (`Kd|Kiso|Kon|Koff|Kfor|Krev|Kreg|Keq|Etot|Lallo`) or test `<:Parameter` structurally, so `name(::Etot)`/`name(::Lallo)` are recognized. Applied as part of the C1 move. |
| **B7** | `test_rate_eq_derivation.jl` perf gate | Replace mean-over-discard-loop with: sink each result into an escaping `Ref`/accumulator and `@test` it (defeats DCE), and take a **minimum** over several `@elapsed` batches instead of a mean. Keep the `allocs == 0` half. Verify the loop isn't elided with `@code_typed`. |
| **B8** | `mechanism_enumeration.jl:1765-2009` → `identify_rate_equation.jl` | Move the ~245-line canonical-rate-eq-hash block to its sole consumer. Update the source file ABOUTME lines. Drives C3's test destination. |
| **B9** | `types.jl:417/770` + `111/788` | Singleton `AllostericEnzymeMechanism` ctor references `_VALID_CAT_ALLO_STATES` instead of an inline tuple literal. Extract `_VALID_REG_ALLO_STATES` const, referenced by both `RegulatorySite` ctor and the singleton ctor. **Keep** the singleton-side validation (it checks type-param-encoded data the concrete ctor never sees). |

**Explicitly skipped:** B10 (`_build_wegscheider_rename_map` relocation) — callers in both
files mean a move only shifts the cross-file back-edge, so it's churn without benefit.

---

## Cluster C — Test reorganization (~1:1 src↔test)

All three files guard **live invariants** (not migration scaffolds); deleting any would
reduce coverage (forbidden). Keep the content, move it into the natural home file, and drop
the standalone `include` from `runtests.jl`.

| Source file | Destination | Rationale |
|---|---|---|
| **C1** `test_chokepoint.jl` | append to `test_types.jl` | The chokepoint is `name(p, m)` in `types.jl`. Apply B6's classifier broadening during the move. |
| **C2** `test_dep_set_invariance.jl` | append to `test_rate_eq_derivation.jl` | Guards the dependent-parameter machinery (`_dependent_param_exprs`), which lives in the rate-eq / thermo layer covered by this file. |
| **C3** `test_canonical_hash_partition.jl` | append to `test_identify_rate_equation.jl` | Guards `_canonical_rate_eq_hash` / dedup determinism, whose sole consumer is `identify_rate_equation.jl` — and B8 moves the code there too, so test and code stay co-located. |

`runtests.jl` loses three `include` lines; the moved testsets run within their host file's
existing scope. Verify the host files already `using` everything the moved testsets need
(or add the imports).

---

## Cluster D — Stale comment / doc sweep (#5)

- Remove the 6 temporal labels in tests: `test_dsl.jl:79` ("new grammar"),
  `test_types.jl:2, 99, 296, 299, 834` ("new design" / "new concrete"). Rename to
  evergreen descriptions of what each testset covers.
- Fix stale src comments flagged by the agent — **verify exact current locations during
  implementation** (the agent's line numbers may be from a stale tree):
  binding-step canonicalization comments near `types.jl:129` and `:319`; the indexed
  `k1f/k1r` doc near `types.jl:864` (the structural-naming refactor replaced indexed
  names).
- CLAUDE.md: update the multiplicity line (A3 — now enumerated); **soften** the
  "public mechanism-construction surface" wording for
  `Mechanism`/`AllostericMechanism`/`init_mechanisms` to say they are accessed as
  `EnzymeRates.X` (the 18-name export set is correct and unchanged); review the guidance
  near the line the agent cited (`:245`).

---

## Sig machinery — what we are NOT doing, and why

The full reflection rewrite (generic `@generated` fieldwalk encoder + tag-dispatched decode
table) would cut the ~93-line Sig block to ~40 (~55% less). We are **not** doing it:

- It is cleverer, not simpler (CLAUDE.md dispreference). The 19 current methods are boring,
  auditable one-liners.
- The decode runs once at compile time inside `@generated`, so its line-count costs nothing
  at runtime — the rewrite buys source lines, not performance.
- It would change the Sig layout (adds explicit type tags), re-baselining Sig-shape goldens
  and degrading isbits-safety from "auditable per method" to "trusted recursion."

We take only the two zero-risk wins already captured above: **B4** (4 metabolite encoders →
1, no layout change) and **B1** (fix the `metabolites()` hand-decode bypass).

---

## Implementation ordering

The clusters are sequenced because one PR means later steps depend on earlier ones:

1. **Cluster A** first (validation + bug fixes), TDD. Restoring balance/duplicate checks may
   surface fixture breakage — fix fixtures or stop-and-ask per A1's note. Run full suite.
2. **Cluster B** behavior-preserving refactors (B1–B9). Run full suite.
3. **Cluster C** test moves — after B6 (classifier) and B8 (code relocation), since
   destinations depend on them. Run full suite.
4. **Cluster D** comment/doc sweep — last (or interleaved where convenient).
5. Full `Pkg.test()` green after each cluster; final green before PR.

**Memory note:** a second agent's full-suite run was OOM-killed (~2.1 GB RSS) at the
compile-budget section; the reviewer's run completed in 9m34s and passed. This looks like
environment memory pressure, not a code defect. Watch RSS during runs; if it recurs,
investigate before assuming a regression.

---

## Testing strategy

- **A1/A2:** direct-constructor tests for each restored check (empty atoms, non-positive /
  Bool counts, empty subs/prods, duplicate names within and across categories, atom
  imbalance), plus an idempotency test for `_add_competitive_inhibitor` on an already-declared
  regulator. These belong in `test_types.jl` (struct ctors) and `test_mechanism_enumeration.jl`
  (the helper).
- **A3:** an enumeration test with a multi-valued `allowed_catalytic_multiplicities` reaction
  asserting candidates appear at each multiplicity; confirm single-valued output is unchanged.
- **A4:** `test_readme_runs.jl` continues to pass with the edited example.
- **B1–B9:** the existing suite is the regression guard (these are behavior-preserving). The
  `rate_equation` perf gate (B7) and chokepoint guard (B6) must still pass; B2/B3 must not
  change `parameters` / `rate_equation_string` golden outputs.
- **C:** the moved testsets run identically inside their host files; partition/invariance
  goldens unchanged.
- Full `Pkg.test()` is the final gate.

## Risks

- **Fixture breakage from restored validation** (A1) — most likely friction point; handled by
  TDD + stop-and-ask on legitimately-unbalanced fixtures.
- **Enumeration count drift** (A3) — only if a test uses a multi-valued multiplicity reaction;
  none currently do, so existing counts should hold.
- **Compile-budget / RSS** — watch memory during full runs (see note above).
- **`@generated` correctness** (B1, B3 move) — verify `metabolites` / parameter outputs are
  byte-identical before and after; these feed the fitting hot path.
