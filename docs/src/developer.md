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
reaction but bound by no step, so they neither appear in `regulators` nor add a
parameter; `Mechanism(em)` lifts back.

## Enumeration engine architecture

Mechanism enumeration uses the **concrete types** `Mechanism` and
`AllostericMechanism` (`src/types.jl`) — internal, not exported — built directly
from `Step` and `Species` values. `Mechanism` has two fields:
`reaction::EnzymeReaction` and `steps::Vector{Vector{Step}}` — kinetic groups,
one inner vector per group holding the steps that share that group's parameters.
`Step` has `from_species`, `to_species`, `consumed`, `released`, and
`is_equilibrium`: going from `from_species` to `to_species`, a step takes up the
metabolites in `consumed` from solution and gives off those in `released`. A pure
binding consumes one metabolite that `to_species` then carries
(`bound_metabolite`); an isomerization has both lists empty (`is_iso`); every
other step, such as fused chemistry and release, a Theorell–Chance step, or
several metabolites on one side, is a transformation. Every step's constants are
named by its two sides — each side's enzyme form followed by its free
metabolites (`K_ES_to_E_S` for `E + S ⇌ E(S)`, `k_EA_B_to_EQ_P` for
`E(A) + B <--> E(Q) + P`). Like the singleton types, these are canonicalized so
that the order or direction in which steps are written does not change the
resulting mechanism. The `Step` constructor stores a pure binding with its
metabolite consumed; the `Mechanism` and `AllostericMechanism` constructors
orient every other step (`_canonical_step_direction`) and sort steps and groups.
They also enforce the kinetic-group rules. A group holds one kind of step
(`_step_kind`: bindings of one metabolite, isomerizations, or transformations
that take up and give off the same metabolites) with one RE/SS flag
(`_assert_uniform_groups`), so a Theorell–Chance step cannot share a group with
a binding. A reaction — the pair of a step's two sides — appears in one step of
one group (`_assert_each_reaction_once`).

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
carries net flux in the child. Flux is decided on the graph of rapid-equilibrium
segments (`_flux_carrying_groups`): each steady-state step is an edge weighted by
its net uptake of substrates minus products, segment offsets included, and a step
carries flux exactly when its biconnected block holds a cycle of nonzero weight. A
steady-state group with no such step exposes only the ratio of its constants, so a
flip set that leaves one counts as failed and is extended, and a split part with
none is emitted at rapid equilibrium. A dead-end copy of a substrate or product
must create a complex that no productive form (one bound to no competitive
inhibitor) duplicates, by composition or by segment and offsets
(`_twin_site_test`), in every conformation where the copy binds; every
conformation holds its free enzyme (`_copy_twin_test`). The dead-end move skips a
pattern whose sites are all twins, the split never isolates twin-only copy sites,
and `_expand_change_allo_state` and the parent check reject a group with only twin
sites (`_duplicate_copy_groups`). A copy whose complexes all duplicate productive
forms has a dwell gauge that absorbs its constant; a complex that duplicates only a
copy's complex never rejects, since the sites that pin either copy may keep both
constants separable.
`expand_mechanisms` asserts both rules on every parent (`_assert_emission_rules`).
Both refinement moves share one minimal-set search (`_minimal_gaining_sets`).
Duplicate equations that survive these proofs are collapsed at compile time by
`eq_hash`. A numerical identifiability rank exists only in the test suite, as an
oracle for the proofs; nothing in `src/` estimates identifiability numerically.

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
denominator is pinned by a test over the uni-bi, bi-bi, and ping-pong
enumerations. A future conformational type (KNF, mnemonic, slow isomerization)
adds one `_requires_hyperbolic_catalysis` method and calls the predicate in its
promotion move.

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
