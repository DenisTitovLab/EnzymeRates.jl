# Multi-`:OnlyA` catalytic states via a Δ0 promote move

## Context

`identify_rate_equation` enumerates biochemically valid mechanisms and fits each. For an
allosteric (MWC) mechanism, every catalytic kinetic group carries one conformational-selectivity
tag — `:EqualAI` (binds both conformations equally), `:OnlyA` (binds only the active
conformation), or `:NonequalAI` (binds both, with different constants). The enumeration cannot
currently produce a mechanism where **two or more distinct catalytic groups are `:OnlyA`** — for
example, an ordered bi-substrate enzyme where both substrates bind only the active conformation.
Such mechanisms are biochemically valid, sit at the same parameter count as their single-`:OnlyA`
siblings, and derive to distinct rate equations, so their absence is a real gap in model
selection.

### The gap, confirmed three ways

**Structural.** The only move that produces a catalytic `:OnlyA` tag is `_expand_to_allosteric`
(`src/mechanism_enumeration.jl`), which fires only on a non-allosteric `Mechanism` and sets
*exactly one* group to `:OnlyA` per variant. Once a mechanism is an `AllostericMechanism`,
`_expand_to_allosteric` is a no-op, and no other move ever adds a second catalytic `:OnlyA`:
`_expand_change_allo_state` only relaxes toward `:NonequalAI`; the regulator and merge moves touch
regulatory sites, not catalytic groups; `_expand_split_kinetic_group` can duplicate a group's tag
but only for the *same* metabolite (and canonicalization usually merges the split back).

**Empirical reachability.** The full breadth-first closure of the expansion moves for a
bi-substrate reaction with one allosteric regulator reaches 26,362 allosteric mechanisms (the
frontier drains — this is the complete reachable set). The maximum number of distinct-metabolite
catalytic `:OnlyA` groups across all of them is **1**; zero mechanisms carry two or more.

**Distinguishability.** Hand-building an ordered bi-uni mechanism and deriving it in both the
1-`:OnlyA` (substrate A only-active) and 2-`:OnlyA` (A and B only-active) assignments yields **the
same five fitted parameters** but structurally different rate equations: with B set to `:OnlyA`,
all B-dependence drops out of the inactive-state (`L`) term of the denominator. The derivation
engine handles multi-`:OnlyA` correctly — the omission is purely in enumeration.

## Goal

Make every catalytic conformational-selectivity assignment reachable — specifically, any subset of
catalytic groups tagged `:OnlyA` — in both the seeded (required-regulator) and general searches,
without a combinatorial blow-up of the seed set or the fit count.

## Non-goals

- **Precalculating multi-`:OnlyA` into the seeds** (folding subset enumeration into
  `_expand_to_allosteric`). This front-loads a `2^k` fan-out into the unpruned base tier and
  multiplies the required-regulator seed set — the cost that required-regulator seeding was built
  to avoid. The cost trade-off is analyzed below and deferred, not adopted.
- **The all-`:EqualAI` catalytic mechanism.** It is genuinely degenerate: if every catalytic group
  (including the catalytic/kcat step) binds both conformations equally, the active and inactive
  states are catalytically identical, `L` folds into `kcat`, and no regulator has an observable
  effect. `_expand_to_allosteric`'s invariant of forcing ≥1 catalytic `:OnlyA` is the correct
  non-degeneracy guard and stays as is.
- **Promoting regulatory ligands.** A regulatory ligand reaches `:OnlyA` directly through
  `_expand_add_allosteric_regulator` (including at a shared site), so promoting a regulatory
  `:EqualAI` ligand would be redundant. The new move touches catalytic groups only.
- Changing fitting, model selection, or the derivation engine.

## The move

Add one expansion move, `_expand_promote_catalytic_to_onlya`, alongside the other
`AllostericMechanism` moves in `src/mechanism_enumeration.jl`:

```julia
"""
    _expand_promote_catalytic_to_onlya(am::AllostericMechanism)
        → Vector{AllostericMechanism}

Δ0 catalytic-state move. For each catalytic kinetic group tagged `:EqualAI`, emit one
variant with that group set to `:OnlyA` — binding (K-type) and iso/catalytic (V-type)
groups alike. The catalytic steps, multiplicity, regulatory sites, and every other tag
pass through unchanged. No-op on a non-allosteric `Mechanism`.
"""
function _expand_promote_catalytic_to_onlya(am::AllostericMechanism)
    results = AllostericMechanism[]
    for g in 1:length(cat_allo_states(am))
        cat_allo_states(am)[g] == :EqualAI || continue
        new_states = copy(cat_allo_states(am))
        new_states[g] = :OnlyA
        push!(results, _with_cat_allo_states(am, new_states))
    end
    results
end

_expand_promote_catalytic_to_onlya(::Mechanism) = AllostericMechanism[]
```

It mirrors the catalytic loop of `_expand_change_allo_state` and reuses the existing
`_with_cat_allo_states` helper (`src/types.jl:701`), which rebuilds the mechanism through the
`AllostericMechanism` constructor and so re-canonicalizes step and tag order.

**Wiring.** Append the move to `_add_expansions_mech!` next to `_expand_change_allo_state`, so it
runs in the general beam (`expand_mechanisms`). It is deliberately **not** added to
`_seed_children`: the seed set stays minimal (one catalytic `:OnlyA` per seed), and the beam
elaborates multi-`:OnlyA` at the same parameter count. Update the `expand_mechanisms` docstring's
move list to include it.

## Why this is correct and safe

- **Covers all catalytic steps.** The loop promotes every `:EqualAI` catalytic group regardless of
  whether its representative step is a binding step or an isomerization. On an already-allosteric
  input, `L` is already revealed by whatever made the mechanism allosteric, so promoting an
  iso/kcat group is observable with no regulator pairing — the reason the lift move
  (`_expand_to_allosteric`) needs the pairing does not apply here.
- **Never degenerate.** The input `AllostericMechanism` already carries ≥1 non-`:EqualAI` catalytic
  group. Promoting only adds conformational selectivity, so the non-degeneracy invariant is
  preserved; the move can never produce an all-`:EqualAI` mechanism.
- **Typically Δ0, and the beam does not depend on it being exactly Δ0.** `:EqualAI` and `:OnlyA`
  each cost one dissociation constant, and `L` already exists, so the derived parameter count is
  normally unchanged (confirmed: 1-`:OnlyA` and 2-`:OnlyA` bi-uni both derive to five fitted
  parameters). If a specific promotion ever shifts identifiability and changes the derived count,
  `_process_batch` buckets each child by its *actual* derived parameter count, so the child is
  still placed and fit correctly. "Δ0" describes the typical cost, not a load-bearing invariant.
- **Reachability and lattice coverage.** `_expand_to_allosteric` seeds one `:OnlyA` group; repeated
  promotion reaches every `:OnlyA` subset; `_expand_change_allo_state` relaxes any group to
  `:NonequalAI`. Together they cover the entire catalytic state lattice except all-`:EqualAI`
  (correctly excluded).
- **Termination.** The move is monotone — it only ever turns `:EqualAI` into `:OnlyA`, never the
  reverse — and bounded by the number of catalytic groups. The `_process_batch` structural
  seen-set (PR #66) is the backstop, exactly as for the other Δ0 moves.
- **Dedup.** Because the constructor canonicalizes tag order, promoting groups `i` then `j` yields
  the same structure as `j` then `i`; duplicates collapse by `hash`.
- **Filters.** Catalytic groups carry no regulator type, so `_filter_by_reg_type` passes trivially;
  the move changes no steps or species, so `_assert_atom_conserving` is unaffected.

## Cost trade-off (recorded, deferred)

This move puts the `:OnlyA` subset combinatorics *inside* the beam, where it is pruned like every
other structural axis, instead of *precalculating* it into the seeds. The two strategies were
compared by number of equations fit:

- **Precalc** pays `2^k` (all `:OnlyA` subsets, `k` = catalytic groups) once at the base tier,
  which is **unpruned** — paid even for subsets the data rejects — and decouples the `:OnlyA` axis
  from the downstream search.
- **The Δ0 move** re-explores the `:OnlyA` lattice at each downstream skeleton the beam reaches
  (RE/SS choice, splits, optional inhibitors/regulators, site merges, `:NonequalAI`
  relaxations). Because `multi-OnlyA + skeleton-A` and `multi-OnlyA + skeleton-B` are distinct
  equations, the seen-set does not dedup them, so the `:OnlyA` cost couples to the number of
  reachable skeletons `D`.

The crossover is not `D` vs `2^k`. Both sides are pruned, so what matters is the per-skeleton
*extra* work `δ` the move does beyond what precalc's base pruning already discarded — the immediate
bad supersets it fits-then-prunes at each skeleton, `δ ≈ k` when `:OnlyA` is discriminable, growing
only when it is not. The comparison is therefore **`D × δ` vs `2^k`**, and three factors bias it
toward the Δ0 move: `δ` is small in the common discriminable case; the move fires only where
single-`:OnlyA` survived (conditional, not the unconditional `2^k`); and it re-decides `:OnlyA` in
each real context, catching context-dependent selectivity (e.g. `:OnlyA` × RE/SS) that base-context
precalc can miss entirely. This move axis is no more explosive than the `:NonequalAI` or
site-merge axes already explored in-beam, and singling `:OnlyA` out for precalc would make it the
lone special-cased axis.

The genuinely precalc-favoring corner — a rich skeleton space *and* `:OnlyA` that stays expensive
to discriminate downstream — is narrow and empirical. **If the fit count turns out huge in
practice, revisit the trade-off then**; the targeted fix at that point is a local guard (skip the
promote past some skeleton complexity, or cache per-topology which `:OnlyA` patterns never
survive), not moving `2^k` mechanisms into the unpruned base tier.

## Testing (TDD)

Unit tests for the move (in `test/test_mechanism_enumeration.jl`, next to the
`_expand_change_allo_state` tests):

- From an `AllostericMechanism` with `cat_allo_states == [:OnlyA, :EqualAI, :EqualAI]`, the move
  emits exactly two children — `[:OnlyA, :OnlyA, :EqualAI]` and `[:OnlyA, :EqualAI, :OnlyA]` — and
  nothing else.
- No-op when no catalytic group is `:EqualAI` (all-`:OnlyA` or already-relaxed input).
- No-op on a `Mechanism`.
- Promoting an **iso/catalytic** group (not just a binding group) is included among the children,
  confirming V-type coverage.
- Order independence: promoting `i` then `j` and `j` then `i` produce `==` mechanisms.

Behavioral tests:

- **Δ0:** a promoted child's `fitted_params` count equals its parent's (bi-uni 1-`:OnlyA` →
  2-`:OnlyA`, both five).
- **Distinctness:** the promoted child's `rate_equation_string` differs from the parent's.
- **Reachability (integration):** the breadth-first closure for a bi-substrate + regulator reaction
  now reaches a mechanism with ≥2 distinct-metabolite catalytic `:OnlyA` groups (the probe that
  currently returns 0), and the enumeration still terminates.

Contract guard:

- Add a multi-`:OnlyA` mechanism (e.g. an ordered bi-substrate with two catalytic `:OnlyA` groups)
  to `MECHANISM_TEST_SPECS` so the golden derivation and the `rate_equation` performance contract
  (`allocs == 0`, `t < 120e-9`, `test/test_rate_eq_derivation.jl`) cover the newly reachable
  mechanisms.

## Files touched

- `src/mechanism_enumeration.jl` — new move + one line in `_add_expansions_mech!` + docstring update.
- `test/test_mechanism_enumeration.jl` — move unit tests + reachability integration test.
- `test/mechanism_definitions_for_test_enzyme_derivation.jl` — one `MECHANISM_TEST_SPECS` entry.
- `test/test_rate_eq_derivation.jl` / golden — coverage follows from the new spec entry.
