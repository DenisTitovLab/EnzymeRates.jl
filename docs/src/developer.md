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
`Step` has `from_species`, `to_species`, `bound_metabolite`, and
`is_equilibrium`. Like the singleton types, these are canonicalized so that the
order or direction in which steps are written does not change the resulting
mechanism.

These are ordinary value types to avoid excessive precompilation costs. The enumeration builds,
expands, and deduplicates many thousands of candidate mechanisms (see
[The enumeration engine](@ref)). If each candidate were a singleton parametric
type, merely constructing it would trigger compiler specialization — the same
precompilation cost the derivation pays — multiplied across thousands of
mechanisms, which would be prohibitive. Carrying mechanism structure as runtime
values instead means enumeration and deduplication cost no compilation at all.
Only the candidates the search actually fits are lifted to singleton types, one
at a time, through `compile_mechanism`.

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
