# Developer / Architecture

This page is for maintainers and contributors. It documents the internal
architecture of EnzymeRates.jl; the public API lives on the
[API Reference](@ref) page.

The package keeps two parallel mechanism representations, split along a
deliberate performance boundary. **Singleton types** drive the compile-time
rate-equation derivation; **concrete types** drive the mechanism enumeration.
The two halves and the optimizer layer are described below.

## Derivation architecture

A mechanism's rate equation is derived once, at compile time, by the
King–Altman/Cha method in `@generated` functions (`src/rate_eq_derivation.jl`).
To make that possible, a mechanism is encoded as a **singleton type** — the
exported `EnzymeMechanism{Sig}` or `AllostericEnzymeMechanism{CatalyticMech,
CatSites, RegSites}` (`src/types.jl`) — that carries all of the mechanism's
structure as Julia type parameters, canonicalized so that the order or direction
in which the steps were written does not change the resulting type. The compiler
specializes the derivation per type, so the symbolic King–Altman work runs
during precompilation rather than at every call.

This is the most important architectural decision in the package, and it is
deliberate. Moving the derivation to compile time leaves
`rate_equation(m, conc, params)` as a flat numeric expression that must be
**allocation-free and sub-120 ns per call**, enforced by
`test_rate_equation_performance` (`test/test_rate_eq_derivation.jl`, asserting
`allocs == 0` and `t < 120e-9` for every fixture mechanism). That speed is the
binding constraint on the whole package: the fitter is a multi-start, global,
gradient-free optimizer that evaluates `rate_equation` millions of times per
fit, and a single rate equation can take minutes to fit, so any per-call
allocation or microsecond-scale overhead would make fitting — and therefore
`identify_rate_equation`, which fits thousands of candidates — impractical.
`loss!` is held to the same standard: `FittingProblem` pre-allocates its
`log_ratios_buffer` once and `loss!` reuses it, so the inner optimization loop
allocates nothing. The cost is the flip side of the benefit: each unique
singleton type triggers a full symbolic derivation, and a mechanism with many
enzyme forms or steps can be slow to compile, exhaust memory, or `StackOverflow`.

`compile_mechanism` (internal) is the boundary between the two representations:
`compile_mechanism(m::Mechanism) = EnzymeMechanism(m)` and
`compile_mechanism(am::AllostericMechanism) = AllostericEnzymeMechanism(am)`. The
lift `EnzymeMechanism(m::Mechanism)` first drops regulators declared on the
reaction but bound by no step, so they add no parameter and the reaction of
`Mechanism(em)`, which lifts back, does not list them.

Only the `@generated` methods (`rate_equation`, `parameters`, `fitted_params`,
`metabolites`, `_kcat_forward`) do work per singleton type: `rate_equation` and
`_kcat_forward` run the derivation, `fitted_params` and the reduced `parameters` run
its constraint solve, and `metabolites` and the full `parameters` read names off the
lifted mechanism. The lifts `Mechanism(em)` and `AllostericMechanism(aem)` read the
type parameters at run
time, and the derivation helpers that take a singleton or its type take it
`@nospecialize`, so they compile once for all mechanisms. A new mechanism type
then pays for its generated bodies and the derivation they run, plus a few thin
methods that still specialize on it: the forwarders that supply the default mode
(such as `rate_equation(m, concs, params)` and `rate_equation_string(m)`), `show`,
`catalytic_mechanism` and `catalytic_multiplicity`. `FittingProblem` and `loss!`
specialize on the type on purpose: `loss!` then reads `fitted_params` and
`metabolites` as constants and calls the generated `rate_equation` without
dispatch.

## Enumeration engine architecture

Mechanism enumeration uses the **concrete types** `Mechanism` and
`AllostericMechanism` (`src/types.jl`) — internal, not exported — built directly
from `Step` and `Species` values. `Mechanism` has two semantic fields:
`reaction::EnzymeReaction` and `steps::Vector{Vector{Step}}` — kinetic groups,
one inner vector per group holding the steps that share that group's parameters.
It also carries `naming`, a lazily filled cache of data derived from the steps (the
free-enzyme forms and each group's naming representative and rendered sides), which
takes no part in `==`, `hash`, the compiled type or display; `AllostericMechanism`
carries the same cache.
`Step` has `from_species`, `to_species`, `consumed`, `released`, and
`is_equilibrium`: going from `from_species` to `to_species`, a step takes up the
metabolites in `consumed` from solution and gives off those in `released`. Steps
are classified by these lists. A binding takes up exactly one metabolite and gives
off none (`bound_metabolite`, `is_binding`). It is plain when `to_species` is
`from_species` with that metabolite added, the residual unchanged and the
conformation free to change (`_binds_ligand`), and fused otherwise, as when the
last substrate binds straight into the product-bound form
(`E(A) + B → E(P, Q)`). An isomerization has both lists empty (`is_iso`); a
Theorell–Chance step takes up one metabolite and gives off another. `_is_chemistry`
is true for every step but a plain binding; the allosteric moves,
`_onlya_haldane_violation`, the steady-state pivot tie-break in
`_assemble_constraints` and `_fused_substrate_binding` read it to tell catalysis from
binding. The derivation reads every binding alike, so a fused binding's rapid-equilibrium
constant is a dissociation constant and its steady-state pair a binding and a
release rate constant. Every step's constants are named by its two sides — each
side's enzyme form followed by its free metabolites (`K_ES_to_E_S` for `E + S ⇌ E(S)`,
`K_EPQ_to_EA_B` for `E(A) + B ⇌ E(P, Q)`, `k_EA_B_to_EQ_P` for
`E(A) + B <--> E(Q) + P`). Like the singleton types, these are canonicalized so
that the order or direction in which steps are written does not change the
resulting mechanism. The `Step` constructor stores a step that takes up nothing
and gives off one metabolite as the binding it reverses, so every binding is
stored with its metabolite consumed; the `Mechanism` and `AllostericMechanism`
constructors orient every other step (`_canonical_step_direction`) and sort steps
and groups. The steps of one group that metabolite progression leaves tied, such as
a conformational change and its mirror at inhibitor-bound forms, turn as one: they
run between the same two conformations in the same order (`_orient_tied_steps`), so
the group's shared constants describe one physical direction. A group whose tied
steps change different pairs of conformations is rejected. The `RegulatorySite`
constructor sorts a site's ligands by name, and the `AllostericMechanism`
constructor sorts the sites. The two mechanism constructors
also enforce the kinetic-group rules. A group's steps take up and
give off the same metabolites (`_step_kind`, the pair `(consumed, released)`) and
carry one RE/SS flag (`_assert_uniform_groups`), so a fused and a plain binding
of one metabolite may share a group and a Theorell–Chance step cannot share one
with a binding. A reaction — the pair of a step's two sides — appears in one step
of one group (`_assert_each_reaction_once`). No two distinct forms render one name
(`_assert_distinct_form_names`): a form's name joins its conformation and bound
metabolites without a separator, so E with NAD and P bound and E with NADP bound
would both be `ENADP` and share their constants' names.

These are ordinary value types to avoid excessive precompilation costs. The enumeration builds,
expands, and deduplicates many thousands of candidate mechanisms (see
[The enumeration engine](@ref)). If each candidate were a singleton parametric
type, merely constructing it would trigger compiler specialization — the same
precompilation cost the derivation pays — multiplied across thousands of
mechanisms, which would be prohibitive. Carrying mechanism structure as runtime
values instead means enumeration and deduplication cost no compilation at all.
Only the candidates the search actually fits are lifted to singleton types, one
at a time, through `compile_mechanism`.

The moves never emit a child that is provably a reparameterization
of its parent. The split move divides a group by binding context and accepts a
set of splits only if the independent-parameter count rises; the count comes from
the thermodynamic constraint solve, and for a `Mechanism` it is evaluated without
building the child, from a cycle basis computed once per parent
(`_partition_independent_count`). The RE→SS move flips whole groups and accepts a
set only if the rapid-equilibrium segment count rises and every flipped group
carries net flux in the child. In a `Mechanism` it never offers as a unit a flank
of a qualifying chain whose isomerization is steady state (`_chain_flank_groups`):
by the chain lemma that flip adds parameters and no rate law. An
`AllostericMechanism` keeps those flips, since over more than one catalytic
subunit a flank's flip can raise the rank. Flux is decided on the graph of
rapid-equilibrium segments (`_flux_carrying_groups`): each steady-state step is an
edge weighted by
its net uptake of substrates minus products, segment offsets included, and a step
carries flux exactly when its biconnected block holds a cycle of nonzero weight. A
steady-state group with no such step exposes only the ratio of its constants, so a
flip set that leaves one counts as failed and is extended, and a split part with
none is emitted at rapid equilibrium. A dead-end copy of a substrate or product
is redundant, and never emitted, when in every conformation where it binds each of
its complexes has a productive twin (a form bound to no competitive inhibitor, in
the complex's rapid-equilibrium segment with its offsets, so the two weights are
proportional: `_productive_twin`, tested over a group by `_all_twin`; every
conformation holds its free enzyme, `_redundant_copy_groups`), and the dwell gauge,
which rescales each twin by
its own factor and all of the copy's complexes by one factor, is consistent with
every kinetic group's shared constants, the copy's own included, an `:EqualAI`
group taking one rescaling in both conformations (`_gauge_rescaling`, bookkeeping
over the steps); a copy whose complexes all have twins but whose gauge fails is
kept, since only the gauge proves the copy's constant invisible. The
dead-end move skips a redundant placement, the split drops every child that holds
a redundant copy group (a split can complete a gauge, and splitting a group that
forms a copy's twins can break one, so only the whole child decides), and
`_expand_change_allo_state` and the parent check reject a redundant group
(`_redundant_copy_groups`). A complex
that duplicates only a copy's complex never rejects, since the sites that pin
either copy may keep both constants separable.
`expand_mechanisms` asserts both rules on every parent (`_assert_emission_rules`).
Both refinement moves share one minimal-set search (`_minimal_gaining_sets`).
Duplicate equations that survive these proofs are collapsed at compile time by
`eq_hash`. A numerical identifiability rank exists only in the test suite, as an
oracle for the proofs; nothing in `src/` estimates identifiability numerically.

`init_mechanisms` appends to the seeds their merged and Theorell–Chance variants
(`_seed_variants`, plain `Mechanism`s only). `_merge_isomerization` removes a
chemistry isomerization and moves every other step at its substrate side onto its
product side, so the last substrate's binding becomes fused; `_eliminate_form`
replaces a form whose only two steps are bindings into it by one Theorell–Chance
step in a group of its own. Every step of a base starts at rapid equilibrium. The
variants are the inclusion-minimal sets of the base's flux-carrying groups
(`_minimal_flips`) whose flip to steady state leaves no rapid-equilibrium
turnover cycle (`_re_turnover_cycle`: each RE step weighted by its uptake of
substrates minus products, a turnover cycle being a cycle of nonzero weight),
keeps a maximal rate both ways (`_has_vmax`, condition V), keeps chemistry out of
equilibrium with both sides (`_chemistry_equilibrates_both_sides`, condition C),
leaves no bottomless segment and keeps every steady-state group flux-carrying.
`_seed_candidate_screen` runs the first three tests on index arrays with a
weighted union-find, so the other two run only on a candidate that passes them. A
merged variant whose every merged complex has two steady-state steps, each alone
in its group, has the family of the unmerged chain with rapid-equilibrium flanks
at the same count and is skipped. `_degenerate`, the failure of V or C on the
active-state graph, marks the starting mechanisms that `_base_tier`
(`src/identify_rate_equation.jl`) replaces by their flip children that are not
degenerate.

Conformational mechanism types declare `_requires_hyperbolic_catalysis` (true for
an `AllostericMechanism` whose catalytic multiplicity is above 1), and the moves
that can give such a mechanism a non-hyperbolic catalytic scheme consult
`_hyperbolic_catalysis`: `_expand_to_allosteric` promotes a non-hyperbolic
parent at multiplicity 1 only, and `_expand_re_to_ss` filters its emitted
children. The predicate scores
the rapid-equilibrium segment graph of the catalytic scheme, leaving out every
step that touches a form carrying a declared inhibitor: for each substrate and
product, a directed steady-state edge scores one if the step binds it in that
direction plus the count the edge's source form carries beyond its segment's
bottom form (`_re_segment_extras`, shared with `_bottomless_re_segment`), and a
segment scores the most any of its forms carries. The equation's degree in the
metabolite is the best score over root segments and spanning arborescences toward
the root; the predicate decides "at most 1" with reachability checks
(`_all_reach`) rather than a max-arborescence solve, since a weight-2 pattern is
one segment or edge scoring 2, a scoring root plus a compatible scoring edge, or
two compatible scoring edges. The exactness of the predicate against the derived
denominator is pinned by a test over the uni-bi enumeration and the
ping-pong-capable bi-bi enumeration, which contains every bi-bi mechanism. A future
conformational type (KNF, mnemonic, slow isomerization) adds one
`_requires_hyperbolic_catalysis` method and calls the predicate in its promotion
move.

## Optimization algorithm architecture

Fitting depends only on Optimization.jl; the package ships no solver of its own.
`fit_rate_equation` (`src/fitting.jl`) wraps `loss!` into an
`Optimization.OptimizationFunction`, builds an `OptimizationProblem`, and calls
`Optimization.solve` with whatever optimizer the caller passes. This gives the
package access to the global, gradient-free optimizers that non-convex
rate-equation fitting needs — CMA-ES (`OptimizationCMAEvolutionStrategy`) and
BBO differential evolution (`OptimizationBBO`) are the tested choices — without
taking on a solver dependency or locking users into one algorithm. Optimization.jl
common options (`maxtime`, `maxiters`, …) are named keyword arguments, and
solver-specific options pass through a `solver_kwargs` named tuple; see
[Loss & optimizers](@ref) for the user-facing view.
