# Skip redundant regulatory-site merges

**Date:** 2026-07-18
**Status:** Approved, ready to plan
**Branch:** `skip-redundant-regulatory-merges`

## Problem

The `_expand_merge_regulatory_sites` enumeration move produces most of the
children in an allosteric search — 64% of them in the PFKP HPC run
(`docs/hpc_results/pfkp_hpc_results/2026_07_18_results`). Almost all of those
children never improve the beam: 9 of 6072 were the best at their parameter
count, and none were expanded further.

A large share of them are also statistically indistinguishable from one another.
At 11 parameters, 5459 of 6072 fitted equations sat within `1e-4` loss of a
neighbor; at 12 parameters, 3192 of 3333. The search fits thousands of distinct
equations that the data cannot tell apart.

The cause is a specific redundancy. When the move merges an activator that binds
only the active conformation (`:OnlyA`) with an inhibitor that binds only the
inactive conformation (`:OnlyI`), the merged single-site mechanism derives to the
**same rate equation** as leaving the two regulators on separate sites. An
`:OnlyA` ligand contributes only to the active-state partition function and an
`:OnlyI` ligand only to the inactive-state one, so co-locating them on one site
changes nothing: separate sites give `(1 + a/K_a)^n` for the active branch and
`(1 + i/K_i)^n` for the inactive branch, and so does the merged site.

The two forms are mathematically equal but render to different parameter
structure, so they hash to different `eq_hash` keys. The fit memo keys on
`eq_hash`, so it never deduplicates them: the search fits the same equation twice.

## Verification

`scratchpad/verify_merge.jl` derived both forms for a PFKP catalytic base and
compared them numerically over 50 random parameter and concentration points:

| merge | states | merged vs separate | keep? |
| --- | --- | --- | --- |
| activator + inhibitor | `:OnlyA` + `:OnlyI` | identical (max rel. diff `2.6e-16`) | **skip** |
| two activators | `:OnlyA` + `:OnlyA` | different (rel. diff `0.07`) | keep |
| antagonist | one ligand retagged `:EqualAI` | different | keep |

The dedup keys of the identical pair differ, confirming the memo misses them.

## Design

Skip only the redundant merge. A merge of sites `i` and `j` is redundant exactly
when the two sites act on disjoint conformations — one site's ligands all bind
the active state and the other's all bind the inactive state.

A ligand's active conformations follow from its allosteric state:

- `:OnlyA` → active state only
- `:OnlyI` → inactive state only
- `:EqualAI`, `:NonequalAI` → both states

So disjoint coverage arises only between an all-`:OnlyA` site and an all-`:OnlyI`
site. Any `:EqualAI` or `:NonequalAI` ligand covers both states and makes the
merge non-redundant, matching the verified cases above.

### Changes (`src/mechanism_enumeration.jl`)

1. A helper returning the conformations a site's ligands act on:

   ```julia
   function _site_active_states(site::RegulatorySite)
       active = Set{Symbol}()
       for st in allo_states(site)
           st in (:OnlyA, :EqualAI, :NonequalAI) && push!(active, :active)
           st in (:OnlyI, :EqualAI, :NonequalAI) && push!(active, :inactive)
       end
       active
   end
   ```

2. `_merged_site_state_assignments` gains a `drop_all_keep` keyword that omits the
   all-keep assignment while still returning the `:EqualAI`-retagged (antagonist)
   variants:

   ```julia
   function _merged_site_state_assignments(base_states::Vector{Symbol};
                                           drop_all_keep::Bool=false)
       assignments = Vector{Symbol}[]
       drop_all_keep || push!(assignments, copy(base_states))
       for i in eachindex(base_states)
           base_states[i] == :EqualAI && continue
           retagged = copy(base_states)
           retagged[i] = :EqualAI
           all(==(:EqualAI), retagged) && continue
           push!(assignments, retagged)
       end
       assignments
   end
   ```

3. `_expand_merge_regulatory_sites` drops the all-keep assignment for a redundant
   pair:

   ```julia
   redundant = isempty(intersect(_site_active_states(sites[i]),
                                  _site_active_states(sites[j])))
   for states in _merged_site_state_assignments(base_states; drop_all_keep=redundant)
   ```

## Why this loses nothing

The dropped all-keep child carries the two regulators on one site with their
current states. Its equation equals that of the mechanism being expanded, which
carries the same regulators on separate sites and is already in the search. The
skip removes a re-render of an equation the search already holds; it removes no
reachable mechanism.

The antagonist retags (`:EqualAI` on one ligand) survive, so the search still
reaches every genuinely distinct merged form.

## Tests (TDD)

1. **Move behavior (the red-green driver)** — `_expand_merge_regulatory_sites` on
   a mechanism with an all-`:OnlyA` site and an all-`:OnlyI` site emits the
   antagonist retags but not the pure all-keep merge. Before the fix the move
   emits that merge, so the test fails; after the fix it passes. The same test
   asserts that same-direction (`:OnlyA` + `:OnlyA`) and any pair carrying an
   `:EqualAI`/`:NonequalAI` ligand still produce the all-keep merge.
2. **Derivation equivalence (premise guard)** — an `:OnlyA` + `:OnlyI` mechanism
   merged onto one site derives to the same rate equation as the separate-site
   form, numerically over random parameter and concentration points. This holds
   before and after the fix; it guards the claim that the skipped child is truly
   redundant, so the skip loses nothing.
3. **`_site_active_states`** — unit coverage of each state's conformation set.
4. **Regression** — the full suite stays green, including the enumeration and
   `rate_equation` performance tests.

## Scope and non-goals

- **In scope:** `_expand_merge_regulatory_sites` only.
- **Out of scope, noted for follow-up:** `_expand_add_allosteric_regulator` can
  place an `:OnlyA` and an `:OnlyI` ligand on one site through a different path
  and may carry the same redundancy. Measure and address separately.
- **Parked:** idea 7's "mirror inhibitor-containing/lacking steps" verification
  awaits a precise definition of "mirror." The primary invariant (inhibitor
  binding steps are always rapid-equilibrium) is already verified: 0 violations
  across 610,946 inhibitor steps in the LDH run.
- **Dropped:** idea 3 (strict-sequential parameter advancing). A faithful
  simulation over the LDH run's parent→child graph avoided 0 fits — the beam
  already advances one parameter count at a time in its parent selection, so
  there is no waste to remove.
