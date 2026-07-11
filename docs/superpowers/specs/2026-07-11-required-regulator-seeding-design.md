# Required-regulator seeding, sign designation, and competitive-site enumeration

## Context

`identify_rate_equation` seeds its beam search with `init_mechanisms`: non-allosteric,
minimum-parameter `Mechanism`s that bind no regulator. The beam then climbs one parameter
count at a time, fitting every surviving mechanism at each level. Reaching a fully-regulated
mechanism therefore means fitting the entire partially-regulated shelf beneath it first.

For phosphofructokinase this is ruinous:

```julia
rxn = @enzyme_reaction begin
    substrates: F6P[C6H13O9P], ATP[C10H16N5O13P3]
    products: F16BP[C6H14O12P2], ADP[C10H15N5O10P2]
    allosteric_regulators: ATP, ADP, Phosphate, F26BP, Citrate
    oligomeric_state: 4
end
```

`init_mechanisms` returns 69 seeds at minimum parameter count — 5 for most, 6 for a few —
binding zero regulators. A fully-regulated mechanism — one dissociation constant per metabolite
and per regulator, plus `kcat` and `L` — sits at 11 parameters, six expansion generations deep. A breadth-first
expansion reaches over 300,000 distinct mechanisms within three generations, and **95.8% of
the reachable space lies below 11 parameters**. The search fits on the order of 100,000
distinct sub-11-parameter equations, every one missing a regulator the experimenter already
knows matters, before it reaches the first useful mechanism. On 1000 cores this still takes
hours.

The experimenter, however, knows the answer's shape in advance: PFK needs all five effectors.
Fitting any equation with fewer bindings is waste.

This design lets the user declare that knowledge and start the search where it should start.

### Relationship to #66

PR #66 gave the beam a structural `seen` set in `_process_batch`: each structurally-distinct
mechanism is processed at most once, so the frontier drains and the search terminates at a
bounded parameter count. This design builds directly on that guarantee — every new move below
is Δ ≥ 0 or Δ = 0 with the seen-set as its termination backstop. No termination machinery is
added here.

## Goals

1. Default the search to requiring every declared regulator, so it starts from fully-regulated
   mechanisms and never fits a partially-regulated one. Let the user move specific regulators to
   optional, per role, when a regulator should be explored rather than assumed.
2. Let the user declare each allosteric regulator's sign (activator or inhibitor), collapsing
   the state combinatorics and expressing real prior knowledge.
3. Explore competitive site-sharing and the activator/antagonist ambiguity, which a
   single-site seed set cannot express, through a fit-gated beam move.

## Non-goals

- Changing fitting or model selection. This design changes only which mechanisms seed the beam.
  It does deliberately change the *default* seed for a reaction that declares regulators (see
  "User surface"), but the previous behavior remains reachable.
- Estimating `Keq` from data (always user-supplied).
- Site-sharing among *required* regulators inside the seed-build. The seed-build keeps one
  site per required regulator; the beam merge move (Part B) reaches shared sites.

## Overview

The work splits into two composable parts.

**Part A — required-regulator seeding and sign designation.** By default every declared
regulator is required: a new `seed_mechanisms` function replaces `init_mechanisms` as the beam's
seed, binding all of them, one site each, in cheap allosteric states. Two keywords —
`optional_allosteric_regulators` and `optional_competitive_inhibitors` — move named regulators
back to optional, so the beam adds them as refinements rather than forcing them into every seed.
Optional DSL sign annotations pin each regulator's state. Part A ships the efficiency win on its
own.

**Part B — competitive-site enumeration.** A new Δ0 beam move, `_expand_merge_regulatory_sites`,
merges two regulatory sites and enumerates the sign-respecting state assignments of the merged
ligands. This reaches two-effector co-binding and the activator/antagonist forms that data
struggles to separate.

Part B builds on Part A but lands independently.

---

## Part A — required-regulator seeding and sign designation

### User surface

**Default: every declared regulator is required.** A reaction that declares
`allosteric_regulators` or `competitive_inhibitors` seeds, by default, only mechanisms that bind
all of them — the common case, where the experimenter already knows the effectors matter. This
deliberately changes the default seed for a regulated reaction; the previous behavior (seeding
from non-regulated `init_mechanisms`) is recovered by marking every regulator optional. Existing
selection tests on regulated reactions therefore change and need review.

**Opt-out keywords, one per role.**
`identify_rate_equation(prob; optional_allosteric_regulators::Vector{Symbol} = Symbol[],
optional_competitive_inhibitors::Vector{Symbol} = Symbol[], …)`. A name in a list moves that
regulator to optional: the beam may add it as a refinement, but no seed is forced to bind it.

The keywords are split by role because a single metabolite may be declared as *both* an
allosteric regulator and a competitive inhibitor — two distinct bindings, for example a
substrate that occupies a regulatory site and also competitively blocks the catalytic site. One
flat list could not say "explore its competitive binding but require its allosteric binding," so
the opt-out is separate per role. The required sets are the declared regulators of each role
minus that role's optional list. When both required sets are empty — every regulator optional,
or none declared — the seed-build is skipped and the base tier is today's `init_mechanisms`.

**DSL sign designation.** Each `allosteric_regulators` entry accepts an optional
`::Activator` or `::Inhibitor` annotation:

```julia
allosteric_regulators: ATP::Inhibitor, ADP::Activator, Phosphate, F26BP::Activator, Citrate::Inhibitor
```

A bare name stays undesignated. The annotation declares the regulator's **observable sign** —
whether it raises or lowers activity — not its molecular mechanism. `Activator` maps to the
`:OnlyA` seed state (binds the active conformation); `Inhibitor` maps to `:OnlyI`.

### Parameter floor

The floor is structural, not a magic number. A seed qualifies when it binds every required
regulator. The parameter count then follows per lineage:

```
floor = base_seed_params + n_required + (1 for L, when a required regulator is allosteric)
```

Each regulator — allosteric or competitive — adds one dissociation constant; a required
allosteric regulator additionally lifts the mechanism to two conformations, adding `L`. For
PFK's five allosteric regulators: `5 + 5 + 1 = 11`. A seed from a 6-parameter init lineage
lands at 12. The seed-build asserts this relationship as an invariant — a "regulated" seed that
misses the expected count signals a bug, such as an accidental `:NonequalAI` doubling.

### `seed_mechanisms(rxn, required_allosteric, required_competitive)`

A new function in `mechanism_enumeration.jl`, taking the two required-regulator sets. It grows
`init_mechanisms` into the fully-required seed set by a breadth-first closure under the structure
moves alone:

- `_expand_to_allosteric` — the only lift from `Mechanism` to `AllostericMechanism`.
- `_expand_add_allosteric_regulator` — restricted to **new sites only**.
- `_expand_add_dead_end_regulator` — for required competitive inhibitors.

Three constraints bound the closure and shape the seeds:

1. **Cheap states only.** Drop any child carrying a `:NonequalAI` tag. The beam reaches
   `:NonequalAI` later through `change_allo_state`.
2. **One site per required regulator.** Keep only new-site additions from
   `_expand_add_allosteric_regulator` (a child with one more regulatory site than its parent).
   This is the constraint that makes the closure tractable — see "Why one site per regulator."
3. **Sign designation.** For a designated regulator, keep only the child in its declared state
   (`:OnlyA` for an activator, `:OnlyI` for an inhibitor). For an undesignated regulator, keep
   both.

The detail moves — `re_to_ss`, `split_kinetic_group`, `change_allo_state` — never run in the
seed-build. The beam applies them later.

The closure keeps a `Set{UInt64}` of `hash(mech)` for termination and to skip repeated work,
mirroring the beam's seen-set. It retains only the seeds — mechanisms that bind all required
regulators — as objects; intermediate nodes contribute their hash and are discarded. The
seed-build is therefore memory-light and runs serially in tens of seconds. Fitting, the
expensive step, stays on the existing `pmap`.

Sub-floor intermediates that a move overshoots (a `+2` lift landing directly at an
all-required-bound child) are kept, so no reachable seed is lost.

### Beam integration

The seed line in `_beam_search` computes the required sets and branches:

```julia
required_allo = setdiff(declared_allosteric(prob.reaction), optional_allosteric_regulators)
required_comp = setdiff(declared_competitive(prob.reaction), optional_competitive_inhibitors)
base = (isempty(required_allo) && isempty(required_comp)) ?
    unique!(collect(init_mechanisms(prob.reaction))) :
    unique!(collect(seed_mechanisms(prob.reaction, required_allo, required_comp)))
```

`identify_rate_equation` threads the two optional keywords into `_beam_search`. Everything
downstream — the seen-set, `_process_batch`, the advancing sweep, the `pmap` fitting — is
untouched. The beam expands seeds with the full move set: it adds optional regulators, applies
the detail moves, and relaxes states to `:NonequalAI`, all as refinements above the floor.

### Sign designation is enforced everywhere

A designated sign constrains the regulator wherever a move chooses its state, in the beam as
well as the seed-build. A designated inhibitor is `:OnlyI` in seeds and may relax to
`:NonequalAI` (the general two-constant model, which parsimony rejects unless the data demands
it), but the beam never emits it as `:OnlyA` or as an EqualAI antagonist of an activator.
Because `change_allo_state` only moves toward `:NonequalAI`, the observable sign can never flip.
An undesignated regulator keeps the full state range.

### Why one site per regulator

`_expand_add_allosteric_regulator` places each new regulator at a new site *or any existing
site*. Across five regulators this generates the full set-partition (Bell-number) space of
site groupings, multiplied by the state combinatorics. Profiled on PFK, the unconstrained
closure exceeds 2,000,000 nodes without converging and exhausts memory.

Restricting the seed-build to new sites collapses this. Profiled on PFK:

| | Unconstrained | One site per regulator |
|---|---|---|
| Closure size | > 2,000,000, non-convergent | 87,223 nodes, ~42 s |
| Fully-regulated seeds | never reached | 11,488 |

The 11,488 undesignated seeds factor exactly as `359 × 2⁵` — 359 catalytic-allostery
skeletons (5–6 per topology across the 69 topologies, each a distinct choice of which single
group becomes the allosteric group) times 2⁵ regulator-state assignments. Every skeleton
yields exactly 32 seeds; every seed carries exactly five one-ligand sites, one per regulator.
The count is fully explained by the combinatorics, with no duplicate inflation.

Sign designation collapses the 2⁵ factor. Designating all five signs yields **359 seeds**;
designating three yields `359 × 2² = 1,436`. Against today's 100,000-plus sub-11 fits, either
is a decisive reduction, and every seed is useful.

Competitive site-sharing among required regulators is not lost — the beam merge move (Part B)
reaches it, fit-gated, so only promising seeds spend effort on it.

---

## Part B — competitive-site enumeration (`_expand_merge_regulatory_sites`)

### The move

For an `AllostericMechanism` with two or more regulatory sites, merge a pair of sites into one
shared site and enumerate the Δ0-valid state assignments of the merged ligands, filtered by
each ligand's declared sign. Each ligand keeps its single dissociation constant, so every
child holds the same parameter count as the parent — the move is Δ0.

A shared site models mutually-exclusive (competitive) binding; separate sites model independent
binding. Merging therefore explores a genuinely different hypothesis at no parameter cost.

### Sign-respecting state assignments

For a merged {activator, inhibitor} pair the move emits three children:

| Merged states | Meaning |
|---|---|
| `{OnlyA, OnlyI}` | two effectors co-bind, competing for the shared site |
| `{EqualAI, OnlyI}` | the activator becomes a pure antagonist of the inhibitor (↑ activity) |
| `{OnlyA, EqualAI}` | the inhibitor becomes a pure antagonist of the activator (↓ activity) |

`{EqualAI, EqualAI}` is degenerate — an all-EqualAI site has no allosteric effect — and is
dropped. The retag preserves each ligand's observable sign: an activator can become an
antagonist of an inhibitor (still raising activity) but never an antagonist of an activator.
Same-sign pairs therefore emit only the co-binding child; a retag there would flip a sign.
Undesignated ligands get every Δ0-valid assignment.

This is why the sign designation is the observable sign, not the molecular mechanism: each
sign spans a direct form and an antagonist form that data cannot easily separate, and the
search offers both.

### Δ0 and distinctness are confirmed

The move rests on two claims, both verified on a PFK seed. Building the unmerged form (activator
and inhibitor at separate sites) and the two merged forms and compiling each:

```
unmerged   {OnlyA | OnlyI}  separate sites   np = 11
plain-merge {OnlyA, OnlyI}   one site         np = 11
antagonist {EqualAI, OnlyI}  one site         np = 11
```

All three carry 11 parameters (Δ0 confirmed), and all three are pairwise rate-equation-distinct
(distinct `eq_hash`). The antagonist form is not a renaming variant of the direct form, so
enumerating both is meaningful and cross-validation can, given enough data, prefer one.

### Bounding the Δ0 volume

The merge move reintroduces the site-partition combinatorics that the seed-build excludes, now
at the floor parameter level. Three mechanisms bound it:

1. **The seen-set** (#66) guarantees termination: each merged structure is processed once.
2. **Fit-gating**: the move runs only on parents the beam selects, so only promising mechanisms
   spawn site-sharing variants.
3. **Merge canonicalization**: merge order must not matter. `merge(1,2)` followed by a merge
   with site 3 must canonicalize to the same three-way shared site as any other order, or Δ0
   variants multiply. The canonicalization reuses the mechanism constructor's existing
   site-ordering.

The residual risk is volume at the floor level when the beam retains many seeds, each spawning
on the order of Bell(k) merge children. The plan includes a volume measurement on PFK and, if
needed, a per-parent cap on merge children. This tuning does not affect Part A's win, which
skips the sub-floor shelf regardless.

### Placement

`_expand_merge_regulatory_sites` joins `_add_expansions_mech!` as a seventh move, dispatched on
`AllostericMechanism`; a `Mechanism` overload returns empty. It is a beam move only and never
runs in the seed-build.

---

## Reachability

Restricting the seed-build to cheap states and one site per regulator loses no reachable
mechanism above the floor:

- **`:NonequalAI` states** — reached from a cheap-state seed by `change_allo_state` (+1). Every
  `:NonequalAI` variant of a seed is a beam refinement.
- **Steady-state groups and splits** — reached by `re_to_ss` and `split_kinetic_group`, the
  detail moves the beam applies.
- **Optional regulators** — added by the beam's `add_allosteric_regulator`, at new or existing
  sites, as effectors or (undesignated) EqualAI antagonists of seeded effectors.
- **Competitive site-sharing among required regulators** — reached by the Part B merge move.

Because every move is Δ ≥ 0, no descendant of a floor seed falls below the floor. The
sub-floor space is traversed once, in the seed-build, and never fit.

## Testing

**Seed-build (`seed_mechanisms`).**
- Seed count on a small reaction matches the `skeletons × states` combinatorics.
- Every seed binds all required regulators; every required regulator sits at its own
  single-ligand site.
- A designated sign collapses that regulator to one state; undesignated seeds carry both.
- No seed carries a `:NonequalAI` tag or a detail-move artifact (steady-state group, split).
- The per-lineage floor invariant holds (`base + n_required + L`).
- A dual-role metabolite (declared as both an allosteric regulator and a competitive inhibitor):
  moving one role to its optional list keeps the other role required and bound in every seed.
- With every declared regulator marked optional (or none declared), `seed_mechanisms` is not
  called and the base tier equals today's `init_mechanisms` output.

**DSL.**
- `::Activator` / `::Inhibitor` parse and attach the sign; bare names stay undesignated.
- An unknown tag (`::Foo`) errors with a clear message.
- Multiplicities and signs combine (`ATP(1,2)::Inhibitor`).

**Merge move (`_expand_merge_regulatory_sites`).**
- A merged pair holds the same parameter count as its unmerged parent (Δ0).
- The unmerged, co-binding, and antagonist forms are pairwise rate-equation-distinct.
- Sign filtering: a designated activator never appears as an antagonist of an activator.
- Merge order canonicalizes: distinct merge sequences reaching the same partition produce one
  structure.
- A `Mechanism` and a single-site `AllostericMechanism` yield no children.

**End-to-end.**
- `identify_rate_equation` on a small regulated reaction (default, no optional lists) fits only
  fully-regulated mechanisms and terminates.
- The same reaction with every regulator marked optional selects the same model as today.

## Follow-ups

Recorded, not built here:

- **V-type catalysis driven by an optional regulator.** `_expand_to_allosteric` fires only on a
  non-allosteric `Mechanism`, so a catalytic-`OnlyA` (V-type) mechanism paired with an *optional*
  regulator is unreachable from an allosteric seed. Niche.
- **Strict designation.** An option to forbid `:NonequalAI` on a designated regulator, for the
  maximal filter, at the cost of the general two-constant model.
- **Required antagonists.** Seeding a required regulator as an EqualAI antagonist, or two
  required effectors sharing a site, needs site-sharing in the seed-build and a re-tag
  capability no current move provides.

## Sequencing

Land Part A first: the two optional-regulator keywords, the DSL sign annotation, and
`seed_mechanisms`, each
with its tests. Part A ships the efficiency win and is behavior-preserving when no regulators
are required. Then land Part B: the merge move, its canonicalization, and the floor-volume
measurement. Bump the package version after Part A and again after Part B.
