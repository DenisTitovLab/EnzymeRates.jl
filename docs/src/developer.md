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
To make that possible, a mechanism is encoded as a **singleton type** —
`EnzymeMechanism{Sig}` or `AllostericEnzymeMechanism{CatalyticMech, CatSites,
RegSites}` (`src/types.jl`) — that carries all of the mechanism's structure as
Julia type parameters. The compiler specializes the derivation per type, so the
symbolic King–Altman work runs during precompilation rather than at every call.

This is the most important architectural decision in the package, and it is
deliberate. Moving the derivation to compile time leaves
`rate_equation(m, conc, params)` as a flat numeric expression that must be
**allocation-free and sub-100 ns per call**, enforced by
`test_rate_equation_performance` (`test/test_rate_eq_derivation.jl`, asserting
`allocs == 0` and `t < 100e-9` for every fixture mechanism). That speed is the
binding constraint on the whole package: the fitter is a multi-start, global,
gradient-free optimizer that evaluates `rate_equation` millions of times per
fit, and a single rate equation can take minutes to fit, so any per-call
allocation or microsecond-scale overhead would make fitting — and therefore
`identify_rate_equation`, which fits thousands of candidates — impractical.
`loss!` is held to the same standard: `FittingProblem` pre-allocates its
`log_ratios_buffer` once and `loss!` reuses it, so the inner optimization loop
allocates nothing.

The cost of compile-time derivation is the flip side of its benefit. Each
unique singleton type triggers a full symbolic derivation, and a mechanism with
many enzyme forms or steps can be slow to compile, exhaust memory, or
`StackOverflow`. `identify_rate_equation` caps expansion via `max_param_count`,
which bounds search depth but not per-mechanism compile cost, so one very large
mechanism can be slow even at a low cap.

`Sig` is purely structural: two mechanisms that differ only in the order their
steps were written collapse to the same `EnzymeMechanism` type, because the
`Mechanism` constructor canonicalizes step and group order before `_sig_of`
runs (see *Enumeration engine architecture* below). `AllostericEnzymeMechanism`
uses three separate type parameters rather than one value-tuple `Sig` because
its first slot is a `DataType` (a concrete `EnzymeMechanism` subtype), and Julia
rejects a `DataType` in the value-tuple position of a type parameter.

`compile_mechanism` (internal) is the boundary between the two representations:
`compile_mechanism(m::Mechanism) = EnzymeMechanism(m)` and
`compile_mechanism(am::AllostericMechanism) = AllostericEnzymeMechanism(am)`. The
lift `EnzymeMechanism(m::Mechanism)` first drops regulators declared on the
reaction but bound by no step, so they neither appear in `regulators` nor add a
parameter; `Mechanism(em)` lifts back.

## Enumeration engine architecture

Mechanism enumeration uses the **concrete types** `Mechanism` and
`AllostericMechanism` (`src/types.jl`), built directly from `Step` and `Species`
values. `Mechanism` has two fields: `reaction::EnzymeReaction` and
`steps::Vector{Vector{Step}}` — kinetic groups, one inner vector per group
holding the steps that share that group's parameters. `Step` has
`from_species`, `to_species`, `bound_metabolite`, and `is_equilibrium`.

These are ordinary value types, and that is the point. The enumeration builds,
expands, and deduplicates many thousands of candidate mechanisms (see
[The enumeration engine](@ref)). If each candidate were a singleton parametric
type, merely constructing it would trigger compiler specialization — the same
precompilation cost the derivation pays — multiplied across thousands of
mechanisms, which would be prohibitive. Carrying mechanism structure as runtime
values instead means enumeration and deduplication cost no compilation at all.
Only the candidates the search actually fits are lifted to singleton types, one
at a time, through `compile_mechanism`.

Canonicalization makes the concrete representation well-behaved, and it is
load-bearing rather than cosmetic. The `Step` and `Mechanism` /
`AllostericMechanism` constructors put every mechanism into a canonical form:
binding steps store the bound metabolite on `to_species` (so a product-release
step and a product-binding step share one representation); iso steps are
oriented to the physical-forward direction by `_canonical_iso_direction`
(atom balance first, then binding-graph context, then a lexical fallback); and
steps and groups are sorted by a canonical key. Two consequences follow. First,
deduplication is a pure `==` / `hash` comparison with no mutation, because
construction has already normalized order. Second, the Haldane/Wegscheider
reduction picks dependent parameters by step order, so without canonical order
the reduced rate equation and `fitted_params` would differ between two
mechanically identical mechanisms written in different orders.

The pipeline — `init_mechanisms`, `expand_mechanisms`, and the `unique!`-based
deduplication — runs end to end on these concrete structs; there is no separate
intermediate representation.

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

`FittingProblem` carries `scale_k_to_kcat::Union{Real,Nothing}`: a `Real`
selects relative data (per-group-centered loss plus rescaling of the SS rate
constants to that kcat), and `nothing` selects absolute per-enzyme turnover
(uncentered loss, no rescale). `fit_rate_equation` returns
`(params, loss, retcode)`, where `retcode` is the best restart's
`Symbol(sol.retcode)`.
