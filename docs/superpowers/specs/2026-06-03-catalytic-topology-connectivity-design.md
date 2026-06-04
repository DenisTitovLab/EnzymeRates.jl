# Catalytic-topology connectivity fix — design

**Date:** 2026-06-03
**Branch:** `change-enumeration-strategy-a-bit`
**Goal:** Fix a structural bug in `_catalytic_topologies` that emits enzyme-form graphs
with **dangling single-metabolite forms** (two forms that differ by exactly one bound
metabolite but have no binding edge between them). Land the combiner fix **and** an
enumeration-test overhaul (a connectivity invariant on every produced mechanism, plus
specific-mechanism assertions replacing bare count checks) in **one PR**. The
beam-search / dedup / base-tier redesign is explicitly deferred and will be re-premised on
the corrected enumeration numbers this fix produces.

---

## Background — the bug (verified)

`_catalytic_topologies` produces topologies whose enzyme-form graph is not fully
connected: a form is created (e.g. `E -(Pyruvate)-> EPyruvate`) but its one-metabolite
neighbour toward the catalytic complex is missing (no `EPyruvate -(NADH)-> ENADHPyruvate`),
leaving `EPyruvate` as a single-metabolite dead-end. **8 of 11 LDH bi-bi topologies are
malformed this way.**

Verbatim example (topology 1, LDH `NADH+Pyruvate ⇌ Lactate+NAD`):

```
forms: E, ENADH, EPyruvate, ENADHPyruvate, ENAD, ELactate, ELactateNAD
  E       -(NADH)->     ENADH
  E       -(Pyruvate)-> EPyruvate         ← created…
  ENADH   -(Pyruvate)-> ENADHPyruvate
  E       -(Lactate)-> ELactate
  E       -(NAD)->     ENAD               ← created…
  ELactate-(NAD)->     ELactateNAD
  ENADHPyruvate -(ISO,SS)-> ELactateNAD
  MISSING: EPyruvate→ENADHPyruvate, ENAD→ELactateNAD
```

**Consequences.** The malformed sequential topologies all carry the *same* full 7-form set
(free E binds all four metabolites), differing only in which second-binding edges exist.
Under rapid equilibrium a form's concentration depends only on the form set, so they
collapse to the **same rate equation** — manufacturing redundant work and confusing output:

- LDH `init_mechanisms` = **77 structures → only 32 distinct rate equations** (45 redundant).
- ter-ter `init_mechanisms` = **74,995 structures**.

**Provenance (git bisect — NOT the recent refactors).** `init_mechanisms` counts are
byte-identical from **#29 (`b3038e4`, "Add biochemical constraints to catalytic topology
enumeration", 2026-04-08)** through the current concrete-type tree:

| commit | LDH | uni-bi | ter-ter |
|---|---|---|---|
| #28 `7654034` (competition patterns) | 55 | 3 | hangs/OOM |
| **#29 `b3038e4` (weak-ordering combiner)** | **77** | 3 | **74,995** |
| #31 `95adc10` | 77 | 3 | 74,995 |
| #38 `7f13907` (just before refactor) | 77 | 3 | 74,995 |
| #40 `80367d0` / present | 77 | 3 | 74,995 |

The dangling-form behaviour was introduced with #29's weak-ordering combiner. The
struct/parameter/concrete-type refactors changed nothing here.

**Why the count tests never caught it.** `@testset "_catalytic_topologies"` asserts only
`length(topos) == N`. The connectivity invariant was never asserted, and **the bug does not
change cardinality** — the correct generator yields the *same* 11 for bi-bi (3 substrate
orderings × 3 product orderings + 2 ping-pong). So `== 11` passes for both the buggy and the
corrected generator. The cardinality tests are change-detectors, not correctness checks.

---

## Root cause (located)

- The backtracker (`mechanism_enumeration.jl:197-504`) builds complete, correct catalytic
  **paths** — each a full route `E → bind substrates → isomerize → release products → E`,
  carrying full binding/release **history** and correctly handling atoms, isomerization and
  ping-pong. **These paths are correct.**
- The bug is in the **combiner** that reassembles paths into topologies.
  `_steps_for_ordering` (`mechanism_enumeration.jl:570-595`) selects *individual* binding
  steps, including a step `src -(m)-> dst` iff `src`'s bound metabolites are all in the
  cumulative weak-ordering levels seen so far — an **upper bound only**. Free `E` (empty
  bound-set) is trivially "accessible", so the combiner emits `E -(m)-> ·` for *every*
  metabolite `m`, including metabolites the ordering places later, while omitting the
  matching second-binding edge. That is the dangling form.
- **A per-step lower-bound check does not work.** Requiring `src` to also contain all
  earlier-level metabolites breaks **ping-pong**: after the covalent step a consumed
  substrate leaves the bound-set, so a legitimate `Estar -(m)-> EstarM` edge is, by
  bound-set alone, indistinguishable from an illegitimate `E -(m)-> EM`. Measured: that fix
  raised LDH 77→181, ter-ter 74,995→320,245, and left violations on the product and
  ping-pong sides. The disambiguating history lives only in the whole paths.

---

## The invariant (exact)

> For any two enzyme forms `A`, `B` present in a mechanism that are **identical in
> conformation and residual** and whose bound-metabolite sets satisfy
> `bound(B) = bound(A) ∪ {m}` for **exactly one** metabolite `m` (they share every bound
> metabolite except one), there MUST be a binding step connecting them (canonical direction
> `A + m → B`).

Equivalently: no form may have a one-metabolite neighbour present in the mechanism without
an edge to it. This forbids single-metabolite dead-ends and partially-connected binding
"squares". It is scoped *within* a conformation+residual, so isomerization / ping-pong
conformational changes are never flagged.

---

## The fix — union of whole consistent paths

Replace step-cherry-picking with **whole-path union**:

> `topology(sub_weak_order, prod_weak_order)` = ⋃ { `path` in the iso-group :
> `substrate_consumption_order(path)` linearizes `sub_weak_order` **AND**
> `product_release_order(path)` linearizes `prod_weak_order` }.

**Why it's correct.** Each path is a complete route, so no form in it dangles. The union of
paths that all linearize *one* weak ordering stays fully connected: the only way to get a
dangling square (`{x}` and `{x,m}` present, no `{x}→{x,m}` edge) is to union paths of
*contradictory* binding orders, and contradictory orders can never both linearize a single
weak ordering — so they are never unioned. Verified by hand for strict orders, ties
(same-level metabolites), ≥3 substrates, and ping-pong. Ping-pong is handled for free:
orderings inconsistent with the iso pattern have *no* matching paths → empty topology →
skipped; the consumption order is read from the path, never inferred from a bound-set.

**Reuse.** The backtracker, `unique_paths`, the iso-pattern grouping, and `_weak_orderings`
are all unchanged and correct. Only the assembly (`mechanism_enumeration.jl:~597-672`)
changes.

**New helpers (small, pure):**
- `_binding_order(path)::Vector{Symbol}` — substrate `bound_metabolite`s in route (path)
  order, binding steps only.
- `_release_order(path)::Vector{Symbol}` — product `bound_metabolite`s in route order.
- `_linearizes(order::Vector{Symbol}, weak_order::Vector{Vector{Symbol}})::Bool` — `true`
  iff every metabolite of `weak_order` appears exactly once in `order` and the level index
  along `order` is non-decreasing (earlier levels strictly before later; any order within a
  level).

**Delete** `_steps_for_ordering` and the `sub_keys`/`prod_keys`/`union` block it feeds.

---

## Count & downstream impact (re-derive, never assume)

- **bi-bi topologies are expected to stay 11**; the 8 malformed ones become 8 clean ones
  with *distinct* form-sets (NADH-first genuinely lacks `EPyruvate`). Verify, don't assume.
- uni-uni (1), uni-bi (3), ter-ter (283), bi-bi ping-pong (10), and the 51-case: **re-derive
  each** and update the test to the corrected value. Several may be unchanged (the bug is
  content, not cardinality), but each must be re-derived.
- `init_mechanisms` and `_expand_substrate_product_dead_ends` counts **will change** — clean
  topologies have different form-sets, hence different dead-end opportunities. LDH init
  77 → (re-derive); ter-ter 74,995 → (re-derive, expected to drop sharply). These numbers
  feed the deferred beam-search redesign.
- **Compile-budget gate** (`test/test_compile_budget.jl`, ~750 trace-compiles for a bi-bi
  init): if the bi-bi init count changes, this golden may shift. Re-baseline with a one-line
  justification if it moves; do not silently bump it.

---

## Test plan (TDD)

1. **Predicate.** Add a test-suite helper `_connectivity_violations(steps) ->
   Vector{Tuple{Symbol,Symbol}}` implementing the invariant above (compare forms by
   `(conformation, residual, sorted bound-metabolite names)`; flag any pair whose bound-sets
   differ by exactly one metabolite, same conformation+residual, with no connecting step).
2. **RED — invariant everywhere.** Assert `isempty(_connectivity_violations(...))` for every
   produced mechanism, added to **all** testsets covering the three layers Denis named:
   - `@testset "_catalytic_topologies"` — every topology, every reaction case.
   - `@testset "init_mechanisms"` — every produced mechanism.
   - the `_expand_substrate_product_dead_ends` testset — every produced mechanism.
   These fail on current output (8/11 bi-bi topologies, 56/77 LDH init, …).
3. **RED — specific mechanisms.** For uni-uni, uni-bi, bi-bi, a bi-bi ping-pong, and a small
   ter case, assert the **exact expected topology form-graphs** (sets of `(from, metabolite,
   to)` edges), hand-derived, replacing bare `length ==` checks. Ping-pong is the
   highest-risk case and its specific assertion is the key guard.
4. **GREEN.** Implement the path-union combiner; delete `_steps_for_ordering`.
5. **Counts as secondary check.** Update the cardinality assertions to the re-derived
   correct values.
6. `julia --project -e 'using Pkg; Pkg.test()'` green.

---

## Scope

**In scope (one PR):**
- Path-union combiner in `_catalytic_topologies` + delete `_steps_for_ordering`; add the
  three small helpers.
- Connectivity-invariant predicate + assertions in the three enumeration testsets.
- Specific-mechanism assertions for the small reactions (replacing/augmenting count checks).
- Re-derived count assertions; CLAUDE.md "Verified topology counts" updated if any move.

**Out of scope (deferred):**
- Beam-search / dedup / base-tier "fit-all + iteration naming" redesign
  (`.claude/2026-06-03-base-tier-fit-all-and-iteration-naming.md`) — revisited after, on the
  corrected numbers.
- Removing the canonical eq-hash / fit-reuse subsystem — **proven load-bearing**: clean
  topologies still exhibit binding-order rapid-equilibrium degeneracy (different orderings,
  same rate law), so equation-level fit reuse stays. KEEP.
- `rate_equation` derivation — untouched.

**Preserved invariants:**
- `rate_equation` allocation-free and <100 ns/call (perf gate green — the perf fixtures are
  hand-defined `MECHANISM_TEST_SPECS`, not enumerated, so they are unaffected).
- Export count = 18; `rate_equation` derivation untouched.
- NOTE: the generated *mechanism set* changes (malformed topologies replaced by clean ones
  with distinct form-sets), so `Sig` types and trace-compile counts shift. Expected, not a
  regression — see the compile-budget re-baseline note above.

---

## Acceptance criteria

- [ ] `_connectivity_violations` is empty for every topology, init mechanism, and dead-end
      mechanism across all test reactions (uni-uni, uni-bi, bi-bi, bi-bi ping-pong, ter-ter).
- [ ] Specific-mechanism assertions pass against hand-derived clean topologies for the small
      reactions.
- [ ] bi-bi topology count verified (= 11 expected); all other enumeration counts re-derived
      and asserted at their corrected values.
- [ ] Full suite green; perf and (re-baselined if needed) compile-budget gates green.

## Risks

- **Ping-pong** is where the first naive fix broke. The path-union design handles it by
  construction (history from paths), but the ping-pong specific-mechanism test + the
  invariant assertion on the ping-pong reaction are the guards.
- **Hand-deriving correct topologies** for the specific-mechanism tests is error-prone;
  cross-check each against the connectivity invariant, the re-derived count, and manual
  review before locking it in.
- **Larger reactions:** re-derive ter-ter and confirm enumeration wall-clock does not regress
  (path-union iterates weak-orderings × paths, same order of work as today).
