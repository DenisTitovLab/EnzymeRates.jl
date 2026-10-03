# Sub-project B: over-engineering review (read-only)

Branch `step-explicit-stoichiometry`, head b340822, compared with a2a02b1. B's net source
addition is about 525 lines (`mechanism_enumeration.jl` +628/−113,
`thermodynamic_constr_for_rate_eq_derivation.jl` +27/−17). Roughly half of it is docstrings.
All line numbers below refer to `src/mechanism_enumeration.jl` at b340822 unless a file is
named. Test line numbers refer to `test/test_mechanism_enumeration.jl`.

Timing context, from `.superpowers/sdd/2026-09-30-exact-filters/gaugefix-report.md` §5: R6
levels 0–2 take 8.85 s against 7.50 s at a2a02b1, so +18% against a 20% limit. Building
dead-end children before running the gauge cost 9.22 s, so the Mechanism shortcut in
`_dead_end_child` is worth about 0.37 s, or 4%. The split costs 1.18 ms per parent against
1.02 ms before B.

---

## 1. Merge the four "every complex has a twin" predicates into one

**Where.** `_copy_twin_test` (1979–2011, two methods), `_twin_only` (2013–2025), `_all_twin`
(2027–2037), and their callers: the split (1727–1728), `_redundant_copy_groups` (2143, 2156,
2162), and the dead-end kernel (2361–2362, 2408–2409, 2444–2446).

**Today.** The same test, "every complex of this copy group has a productive twin", is
written four times.
- `_productive_twin(groups)` answers it for one site on one state's graph.
- `_all_twin(groups, g, twin)` applies it to group g.
- `_twin_only(twin, m, part, g)` applies it to a part and reads g's tag only to pass it on.
- `_copy_twin_test(m)` applies it per site in every state where the copy binds. Its allosteric
  method builds `_state_mechanism(am, :I)` and special-cases a free enzyme that the inactive
  graph lacks.

The dead-end kernel computes `_copy_twin_test(m)` and also `_productive_twin(steps(m))`. For a
`Mechanism` that computes the same twin function twice (handoff cleanup 1). For an allosteric
parent, the second twin function goes unused, because the allosteric `_dead_end_child` ignores
its last argument (`_`).

**Change.**
- Keep one predicate on one state's graph:
  `_all_twin(part::Vector{Step}, twin) = (L = bound_metabolite(first(part))) isa
  CompetitiveInhibitor && all(s -> twin(from_species(s), L) !== nothing, part)`.
- Make `_gauge_rescaling` start with `_all_twin(groups[g], twin) || return nothing`. Every
  `_all_twin(...) && _gauge_rescaling(...) !== nothing` pair then collapses to the second
  half. Both outcomes mean "not redundant", so nothing is lost.
- Split: `twin = _productive_twin(groups)` and `redundant_part(g, bp) = any(part ->
  _all_twin(part, twin), bp) && begin ... end`. The `begin` block, which builds the child and
  calls `_redundant_copy_groups`, is unchanged.
- Dead-end: drop the `all_twin` argument. The `Mechanism` method becomes
  `_gauge_rescaling(groups, length(groups), twin, identity) !== nothing && return nothing`.
  The allosteric method becomes `_all_twin(groups[end], twin) &&
  !isempty(_redundant_copy_groups(child))`. Both methods now use `twin`, and the kernel computes
  `_productive_twin(steps(m))` once.
- Delete `_copy_twin_test` and `_twin_only`.

**Why the outputs stay identical.**
- *Mechanism.* `_copy_twin_test(m::Mechanism)` equals `_productive_twin(...) !== nothing`.
  The same function is applied, so the verdicts are the same.
- *Allosteric.* The active-state twin test is the active half of `_copy_twin_test`. It is a
  necessary condition for redundancy, so as a guard it can only send more candidates to the
  full test. The full per-state test runs on the built child (`_redundant_copy_groups`), and
  it agrees with `_copy_twin_test`'s inactive-state verdict in every branch:
  - Site in the parent's inactive graph: the copy forms are not productive and mirrors join
    only copy complexes, so segments and offsets of productive forms are unchanged.
  - Free enzyme absent from the inactive graph: `_reachable_from_free`
    (`rate_eq_derivation.jl:1054`) treats E as a root once the `:EqualAI` copy step exists. E
    and E(L\*) are kept, and E(L\*) has no productive twin, so the copy is not all-twin. This
    matches `_copy_twin_test` returning false.
  - Any other site the inactive graph lacks: the copy step is stranded and pruned, so it binds
    nothing there. This matches `_copy_twin_test` returning true.
  - Split: the child is a regrouping of the parent with tags inherited, so the parent's and
    the child's inactive graphs coincide.
- *One deliberate difference.* Today the allosteric dead-end path checks every group of the
  child, while the `Mechanism` path checks only the new group. The proposal keeps that. It can
  differ only if a placement makes an older copy redundant, which B's own argument excludes. If
  it did happen, today's code would emit a child that `_assert_emission_rules` rejects at its
  next expansion, so the stricter check is safer.

**Cost.** R6 does not change: the `Mechanism` work is identical, with one twin construction
fewer per parent. An allosteric candidate that is all-twin in the active state but not in the
inactive one now pays one `_redundant_copy_groups(child)` call. That child was already built.

**Tests that pin it.**
- 6998 `_productive_twin`; 7224–7613 `_redundant_copy_groups`.
- 3842, 3897, 3940: allosteric dead-end and relaxation, exact children plus rank. These pin
  the free-enzyme behaviour at move level.
- 2536, 2575, 2619: split.
- 9700 and 9716: population pins.
- To rewrite:
  - 7194–7221 call `_copy_twin_test` directly. Replace them with `_redundant_copy_groups` on
    built children. 3842 and 3897 already cover the `:EqualAI` cases at move level, and the
    `:OnlyA`-copy case needs one new fixture with a rank pin.
  - 7275 and 7318 call `_twin_only`. Use `_all_twin(steps(case3)[g3], _productive_twin(...))`.
  - `docs/src/developer.md:104` names `_copy_twin_test`.
- Dropping the unit tests of a deleted helper needs Denis's agreement (CLAUDE.md coverage
  rule).

**Size and risk.** About −60 source lines. Risk: low to medium. The allosteric equivalence
rests on the pruning argument above, and 3842 and 3897 pin it.

## 2. Deduplicate the docstrings, and fix the stale and inaccurate ones

The same argument is written in three to five places. Keep each in one home:

| Argument | Keep at | Delete or shorten at |
|---|---|---|
| A complex that duplicates only a copy's complex is no twin; the sites that pin either copy keep both constants separable | `_productive_twin` 1954–1957 | 183–185, 1712–1715, 2122–2125, 2249–2252 |
| All-twin with a failing gauge is kept; the gauge is the proof | `_redundant_copy_groups` 2117–2122 | 182–184, 2245–2249 |
| A placement never makes an older copy redundant, so only the new group is tested | `_dead_end_child` 2292–2296 | 189–191, 2252–2257 |
| Theorem 2's finite merge, K_h + K\* | `_gauge_rescaling` 2060–2062 | 2114–2117, 2243–2245 |
| Twin definition: same RE segment, same offsets, proportional weights | `_productive_twin` 1940–1953 | 179–182, 2111–2112, 2240–2242 |
| A zero-flux SS group's two constants enter only as their ratio | `_flux_carrying_groups` 1495–1500 | 179–180, 1290–1294, 1701–1705, 1891–1893 |
| Competition is decided per site | `_expand_add_dead_end_regulator` 2226–2235 | `_forms_with_binding_step_native` 2177–2184 (reduce to its contract) |

Also:
- **Remove evidence and opaque labels.**
  - 1952 "(four such copies at depth 2 of R6 have full rank)": a measurement, which belongs in
    the record.
  - 2122 "(Case 3's is)": a spec case label, opaque in source.
  - 1713–1715: narration of the tolerated two-copies class, which belongs in the spec.
- **Handoff cleanup 3.** `_gauge_rescaling` 2064–2066 calls the finite ρ = K_h/(K_h + K\*)
  "the" shared factor, while 2047–2053 describe a first-order rescaling. Say instead: "share
  one symbol; Theorem 2's finite merge gives ρ = K_h/(K_h + K\*)".
- **Handoff cleanup 4.** 2238–2239 says the dead-end check runs "`_redundant_copy_groups`,
  through `_dead_end_child`". A `Mechanism` takes the `_gauge_rescaling` path. This fix falls
  out of the deduplication.
- **Stale text, disagreeing with the code.** The `expand_mechanisms` docstring (2975–2977)
  says each parent's "every competitive-inhibitor group creates a new complex". That is the
  withdrawn plain copy rule; under option 3 an all-twin group whose gauge fails passes.
- `_assert_emission_rules` (177–196) shrinks to about 8 lines: the two rules, why a parent
  must obey them (the flip tests only flipped groups), and which moves filter. Its code is 15
  lines.

**Size and risk.** About −50 to −70 lines. Risk: none, since only docs change. Pinned by
nothing.

## 3. Merge the dead-end wrappers into the kernel

**Where.** The two `_expand_add_dead_end_regulator` methods (2265–2285), the
`_expand_add_dead_end_regulator_native` docstring and signature (2312–2323), and the
two-signature header (2215–2219).

**Why it is possible.** B made the kernel return typed children, so the wrappers now only pass
`Set{Symbol}()` or the allosteric ligand names. Grep finds no other caller of `_native`.

**Change.** Write one method `_expand_add_dead_end_regulator(m::Union{Mechanism,
AllostericMechanism}, rxn; exclude_regs)` with `additional_excluded = m isa AllostericMechanism
? Set(name(l) for site in regulatory_sites(m) for l in ligands(site)) : Set{Symbol}()`. This is
a pure refactor, so the outputs are identical.

**Tests.** Every dead-end testset (2669 onward), and seed tests 10330 and 10477.

**Size and risk.** About −18 lines. Risk: low.

## 4. The split's unit pre-filter is not a pure optimization: the spec's claim fails

**Where.** `redundant_part` (1728–1732, 1736) and the docstring at 1708–1712. Spec §4 says:
"Every superset of such a unit keeps the part, and splitting other groups only relaxes the
part's gauge."

**The counterexample is in the test file.** Take H1 (7303–7320, 7 fitted, rank 7, kept) with
its two A groups merged into one, {E + A, E(B) + A}. Both twins, E(A) and E(A, B), are then
formed by one A group and share a class.
- B group 2 reads (ρ,ρ) → (nothing, nothing) beside (1,1) → (nothing, nothing), which is
  consistent.
- Every other group is consistent too, so the merged mechanism's copy is redundant.

Splitting the A group (a non-copy group) separates that class, and B group 2 becomes
(ρ₁, ρ₂) beside (1,1), so the gauge fails. The class rule is anti-monotone under a split of
the group that forms the twins.

**Effect.** Removing the pre-filter can only add children:
- Sets of units that do not contain a pre-filtered unit U are generated and tested identically
  by `_minimal_gaining_sets`. Their candidates and their failed subsets never involve U, and
  inserting U renumbers the units monotonically, so the order is preserved.
- If {U} gains alone, as in test 2536, the post-filter rejects the child for the same reason,
  so nothing changes.
- A difference arises only when {U} fails alone and some gaining superset's partner split
  breaks the part's gauge. Such a child has no redundant group: it satisfies rule 2 and today
  it is lost.

My attempted construction (H1's merged sibling plus A\* at E(Q)) did not produce a difference,
because a Q mirror blocks the gauge until Q is split. The case may not arise in R1–R6.

**Proposal.**
1. Delete the pre-filter on a scratch branch and rerun B's regression harness at R6 levels 2
   and 3 and on the allosteric sample.
2. If nothing changes, delete it: −8 code lines plus its docstring sentence. Item 1's split
   guard then goes with it.
3. If children are added, they are rule-conforming. Denis decides, and spec §4's sentence
   needs correcting either way.

**Cost.** Unmeasured. Apriori would extend supersets of U, while the pre-filter's own child
construction per twin-only unit would disappear.

**Tests.** 2536 passes either way. Pin 9716 (31,730) would catch any R6 level-2 change.

**Risk.** Medium, because the behaviour may change.

## 5. Keep: four measured or proved optimizations

These answer "pre-filters that are pure optimizations of a post-filter".

- **Flip all-steady-state unit filter (1308).** It is pure: by the lift lemma, a set that
  contains a group with no flux on the all-SS graph never gains, so outputs are identical
  without it. It is also essential, because without it Apriori extends every hopeless set.
- **Split `reverted` / bottomless guard (1743–1746).** It is pure, and it skips a
  `_re_segment_extras` call for each candidate without a reverted unit.
- **Pre-construction gauge in the `Mechanism` `_dead_end_child`.** One path (build the child,
  then `_redundant_copy_groups`) would give identical output. It costs 0.37 s on R6 levels 0–2
  (measured), which takes B from +18% to about +23%, over the limit.
- **Possible headroom.** Spec §3 suggests sharing `_re_segment_extras` between
  `_bottomless_re_segment` and `_flux_carrying_groups` in the flip's `gains`. That is
  performance only, but it could buy room for unifying the dead-end path later.

## 6. `_redundant_copy_groups(am)`: two dead clauses and a shorter A/I form

**Where.** 2161–2166.

**Dead clauses.**
- In `binds_inactive = tags[g] !== :OnlyA && !isempty(inactive[g])`, the first clause is
  dead. `_state_allo_mechanism` drops `:OnlyA` groups, and `_assert_each_reaction_once` stops
  any other group from holding the same reaction, so an `:OnlyA` group's `inactive[g]` is
  always empty.
- In the "binds nothing" `Dict`, `h != g` is dead for the same reason: `inactive[g]` is empty.

**Change.** With item 1:
`ra = _gauge_rescaling(active, …, label(:A))`; `ra === nothing && return false`;
`ri = isempty(inactive[g]) ? Dict(h => (nothing, nothing) for h in eachindex(inactive) if
!isempty(inactive[h])) : _gauge_rescaling(inactive, …, label(:I))`; then the existing
`:EqualAI` comparison. A general per-state loop would not be shorter than this two-line form.
Merging the `Mechanism` and allosteric methods is not worth it: one branch gets added, and the
`Mechanism` method is the hot path.

**Tests.** 7538/7539, 7569, 7613, 3842–4016.

**Size and risk.** About −5 lines. Risk: low.

## 7. One helper for the RE-segment index preamble

**Where.** Three functions open with the same six lines: `species, segments, extras =
_re_segment_extras(groups)`; `idx = Dict(...)`; `seg = zeros(...)`; a fill loop. They are
`_flux_carrying_steps` (1427–1432), `_productive_twin` (1960–1965) and the pre-B
`_hyperbolic_catalysis` (1549–1554).

**Change.** Add a wrapper in `mechanism_enumeration.jl` that returns `index` and `segment_of`
as well. Do not change `_re_segment_extras` itself: `_bottomless_re_segment` calls it on every
`Mechanism` construction. This is a pure refactor.

**Tests.** 6778, 6998, 7618, 7902.

**Size and risk.** About −10 lines. Risk: low.

## 8. One "rebuild a step with a new RE/SS flag" helper

**Where.** Four sites write `Step(from_species(s), to_species(s), consumed(s), released(s),
flag)`: 852 (init), 1343 (`_flip_group_to_ss`), 1400 (`_all_steady_state`) and 1899
(`_revert_zero_flux_parts`).

**Change.** Add `_with_equilibrium(s, flag)`. Then `_flip_group_to_ss` (1338–1352) becomes a
one-line comprehension, and `_all_steady_state` becomes `[ _with_equilibrium.(g, false) for g
in groups ]`.

**Tests.** 2227, 8000–8040, 6790–6967, 8285.

**Size and risk.** About −12 lines. Risk: low.

## 9. Inline two single-use helpers

- `_copy_complex` (1933–1935) has one caller, at 2415. It was shared once; deferred minor
  Task 3 records why.
- `_composition` (1922–1931) is used only inside `_productive_twin`. Make it a local closure
  there, and rewrite test line 7075.

**Size and risk.** About −7 lines. Risk: low.

## 10. Tidy `_split_gain_test`

**Where.** 1766–1792.

**Change.** Write one loop over both parts: set the kind for each step, and the new id for
part 2. `base` can use the counter's default kinds. Optionally,
`_revert_zero_flux_parts` can return `(parts, reverted)` instead of relying on `parts !== bp`
identity.

**Tests.** 8348, and the split exact-children tests.

**Size and risk.** About −4 lines. Risk: low.

## 11. `_block_balanced`: one forest BFS (needs Denis's permission, since it is a rewrite)

**Where.** 1426–1508.

**Why it works.** The tree path between two vertices of one biconnected block stays inside
that block. So one BFS over the whole segment graph gives potentials, and a block is
unbalanced iff it holds an edge whose weight is not the potential difference of its ends.

**What it replaces.** `_block_balanced` (24 lines) and the `members` dictionary. The
mathematics is identical.

**Tests.** 6778, plus the population pins.

**Size and risk.** About −12 lines. Risk: low to medium.

## 12. `_redundant_copy_groups(am)`: build the inactive view without a Mechanism round trip

**Where.** 2151–2153.

**Today.** The code builds `_state_mechanism(am, :I)`, which constructs a `Mechanism`, and
then matches steps by `minmax(_forward_sides(s)...)` string labels.

**Change.** Call `_reachable_from_free` on the non-`:OnlyA` groups, and keep a step iff both of
its endpoints are reachable. This is exactly `_state_allo_mechanism`'s rule, so the result is
identical. It touches a helper next to the derivation code, so it is low priority.

**Size and risk.** About −2 lines, and faster on allosteric mechanisms. Risk: low.

## 13. Test only: pin the `ss_mirror` rank (handoff cleanup 2)

**Where.** 7417.

**Change.** Add `@test fitted(ss_mirror) == 8 && _testhelper_identifiable_rank(ss_mirror) ==
6`. Spec §2's table gives 8/6.

## 14. Shorten the `seed_mechanisms` error and docstring

**Where.** 3067–3071 and 3104–3112.

**Change.** The docstring sentence repeats the error message, and the nine-line message
restates the copy rule. Keep the regulator names, one sentence on why, and the keywords.

**Tests.** 10477 asserts, among other phrases, the exact wording "a uni-uni seed with one
conformation has no such site". Shortening the message means editing that assertion with it.

**Size and risk.** About −5 lines. Risk: low.

---

## Places where the code and the spec or docs disagree

1. The split pre-filter's monotonicity claim, in spec §4 and implicitly in docstring 1708–1712
   (item 4).
2. The stale "creates a new complex" in the `expand_mechanisms` docstring at 2975–2977
   (item 2).
3. The `Mechanism` `_gauge_rescaling` path versus the docstring's `_redundant_copy_groups`
   wording at 2238–2239 (item 2; handoff cleanup 4).
4. Scope of the dead-end check: the allosteric path checks every group of the child, the
   `Mechanism` path only the new group. The docstrings state this. It is harmless given B's
   argument, and the allosteric path is the safer of the two.
5. `docs/src/developer.md:104` will name a deleted function if item 1 lands.

## Estimated total

About 150–200 source lines, of B's roughly 525 net added. About 70 of these are docstring
lines. Items 4 and 11 (about 20 lines) depend on a measurement and on Denis's permission.
