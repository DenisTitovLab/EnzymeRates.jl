# Canonical kinetic-group partition at mechanism construction

## Problem

The LDH search fits the same rate equation many times under different constant
names. On the 2026-07-03 HPC run, the search performed **68,085 distinct-`eq_hash`
fits, but only 23,080 are distinct rate functions — 45,005 (66%) are redundant
re-fits.** Fits are memoized by `eq_hash` (`identify_rate_equation.jl:638`), so the
union count is the true fit count; canonical dedup cuts fitting roughly threefold.
The rate grows with parameter count: 10% at iteration 1, 63–65% in the high-param
allosteric iterations.

`eq_hash` misses these because `_rate_eq_dedup_key` (`identify_rate_equation.jl:311`)
hashes the rendered equation string, and the same rate law renders under different
constant names.

### Root cause: redundant kinetic-group partitions

The dupes are the **same mechanism graph** — identical Species, Steps, and RE/SS
assignment. They differ only in how a set of binding constants is made equal: some
mechanisms place those steps in one kinetic group (a shared constant), others split
them into separate groups whose constants a Wegscheider tie then forces equal. Both
routes produce identical constant values and the identical rate law, but a different
`mechanism_type` and a different rendered string, so `eq_hash` treats them as
distinct. (Across 11,447 dup classes, no two members share a `mechanism_type`,
confirming the split is a real encoding difference, not nondeterministic rendering.)

Of 11,447 duplicate classes:

| | Share | Cause | Reachable at construction? |
|---|---|---|---|
| Same graph, different partition | **84.8%** | redundant kinetic-group split that Wegscheider re-merges | **Yes** |
| Different graphs, same rate function | 15.2% | cross-metabolite algebraic coincidence (`kon_NAD`↔`kon_NADH`) | No — needs an algebraic key |

Confirmed pair `3a2788df` / `f5f7e53b` (n_params 6): identical 13-step set, 8 vs 9
kinetic groups. One puts "Lactate binds E·NADH" in the shared Lactate group (one
`K_Lactate`); the other splits it into its own group (`K_Lactate_ENADH`), which a
thermodynamic box then forces equal to `K_Lactate_E`. The split is redundant: it
spawns a separate constant that Wegscheider immediately re-ties.

Measured coverage of a binding-K merge: of the 9,710 same-graph dup classes,
**9,707 (100.0%)** differ only in how RE-binding steps are grouped. Merging kinetic
groups whose binding-K representatives are Wegscheider-tied therefore collapses the
full same-graph population. The 3 remaining classes (0.03%) involve a steady-state
step and fall to the deferred residual.

## Goal

Collapse the same-graph dupes at construction, before derivation or fitting, by
making the kinetic-group partition canonical. Three wins:

- **Fitting:** removes ~93% of the fitting waste (the 84.8% of classes, 74% of dup
  members).
- **Compiled code and memory:** the dedup cuts distinct `@generated` mechanism
  instantiations ~2.6× (68,085 → ~26,000). RAM growth tracks that count
  (`docs/superpowers/specs/2026-07-03-ldh-memory-findings.md`).
- **Code:** the single-symbol Wegscheider collapse moves from the derivation to the
  constructor, so the derivation-time Pass-2 absorption is removed — one fewer code
  path, and a smaller generated expression per mechanism.

## Design

Add **partition canonicalization** to the `Mechanism` constructor, after the
existing direction/order canonicalization:

1. Build a provisional mechanism carrying the caller's partition, already
   canonical in step direction and step/group order.
2. Compute the single-symbol binding-K Wegscheider ties on that mechanism — the
   relation `_build_wegscheider_rename_map` already finds.
3. Union-find the kinetic groups: union two groups when their representative
   binding-K's are tied. Merge each connected component into one group (its steps
   share one kinetic parameter).
4. Re-run step/group ordering on the merged partition and store.

Two mechanisms with the identical step-set and RE/SS that differ only by a
redundant split now build to the **same** `Mechanism`. `Base.hash`/`==` already key
on `(reaction, steps)`, so the existing `unique!` in `_beam_search`
(`identify_rate_equation.jl:644`) collapses them, and the fit memo never fits the
duplicate. `AllostericMechanism` gets the same treatment through the `_state_*` tie
analogs.

### Why merging preserves the model

A split whose constant is Wegscheider-forced equal to another group's constant has
the same rate law **and** the same parameter count as the merged form — that is why
these mechanisms already share a rate function. Merging converts a thermodynamic
constraint into a kinetic-group identity, which also simplifies the derivation: the
merged mechanism carries fewer Wegscheider constraints.

### Feasibility

The tie machinery already operates on the plain struct.
`_dependent_param_exprs_kernel` takes `mech::Mechanism`
(`thermodynamic_constr_for_rate_eq_derivation.jl:305`);
`_build_wegscheider_rename_map` only reconstructs a `Mechanism` from the lifted type
as a convenience (`rate_eq_derivation.jl:118`). Refactor it to accept a `Mechanism`
so the constructor can call it on the provisional mechanism.

**Recursion:** computing ties needs a `Mechanism`, and the constructor is building
one. Build the provisional mechanism through a private path that runs direction and
order canonicalization but skips the partition merge, so the tie computation cannot
re-enter the merging constructor.

## Scope decisions

- **Binding-K single-symbol ties only.** Measured sufficient (100.0% of same-graph
  dupes). Reuses proven machinery; no new tie families, no sign-flip exclusion to
  navigate.
- **Run the tie analysis unconditionally in the constructor.** Correctness first;
  profile enumeration afterward and add a structural pre-check gate only if it
  proves a bottleneck.
- **Defer the 15.2% cross-graph residual.** Genuinely different graphs that coincide
  in rate function need a canonical algebraic key; out of scope here.

## Replace the derivation-time Pass-2 absorption

The constructor merge subsumes `_build_wegscheider_rename_map` Pass-2. Both collapse
the same relation — single-symbol binding-K↔binding-K ties
(`rate_eq_derivation.jl:132-133`) — the constructor structurally (one kinetic group),
Pass-2 by renaming in the polynomial. Once the graph arrives pre-merged, Pass-2 has
nothing to absorb, so remove it: the equality now lives in the mechanism structure,
not a derivation-time rename.

Remove it on evidence, not faith. `project_dedup_pass2_dead_code` records Pass-2 as
load-bearing for the non-competitive-inhibitor and non-essential-activator families,
which produce single-symbol binding-K ties the constructor merge should now cover.
Sequence: (1) implement the constructor merge; (2) assert Pass-2 finds no tie to
absorb across the full suite and those two families specifically; (3) delete Pass-2
only once that assertion is green. The general multi-symbol Wegscheider / Haldane
elimination is a separate path and stays in the derivation — only the single-symbol
binding-K case moves to the constructor.

## Blast radius and tests

The canonical mechanism becomes the merged form, so these move and must be
re-baselined:

- `fitted_params` and the reduced-equation string for affected mechanisms. The
  fitted model is unchanged (same rate law, same parameter count); only the rendered
  form and parameter names change.
- Distinct-mechanism counts drop. The bi-bi partition-count assertion
  (`test/test_identify_rate_equation.jl:928-951`, "exactly 55 classes") should fall;
  confirm the new count and update the comment to record the merge as the cause.
- Golden references: the allosteric golden, flat-string and Expr-shape regression
  tests in `test/test_rate_eq_derivation.jl`, and any `fitted_params` golden.

Unaffected:

- The `rate_equation` 0-alloc / sub-100 ns runtime contract. The merge is
  compile-time construction; the generated numeric body is untouched.
- The parameter-naming chokepoint AST guard (`test/test_types.jl:1577-1644`). Names
  still flow through `name(p, m)`.

The Canonical Step Form guard **extends** here: partition canonicalization joins the
existing direction and order canonicalization as a third construction-time
invariant. Its tests move accordingly.

## Validation

Unit tests on the confirmed pairs are the acceptance bar.

- **Mechanism invariant (TDD, acceptance):** construct the confirmed split and merged
  pair (`3a2788df` / `f5f7e53b`) and a few synthetic same-graph pairs; assert they
  build to an equal `Mechanism` (`==`) and share an `eq_hash`. Write these first,
  watch them fail, implement, watch them pass.
- **Pass-2 removal guard:** before deleting Pass-2, assert it finds no single-symbol
  binding-K tie to absorb across the full suite and the non-competitive-inhibitor and
  non-essential-activator families.
- **Optional aggregate cross-check:** re-run the LDH pipeline and the independent
  fingerprint canonicalizer; the distinct-`eq_hash` count should fall from 68,085
  toward the ~26,000 cross-graph floor. Available if wanted, not required for
  acceptance.

## Out of scope

- The algebraic dedup key for the cross-graph 15.2%.
- Canonical representative naming and allosteric factor ordering (old "option 1").
  These were symptoms of the redundant splits; canonicalizing the partition removes
  most of the ties that drove them. Revisit only if a measured residual justifies it.
