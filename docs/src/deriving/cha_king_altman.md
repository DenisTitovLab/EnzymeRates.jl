# The Cha / King–Altman algorithm

EnzymeRates derives rate equations symbolically with the King–Altman schematic
method [KingAltman1956](@cite), extended by Cha's rapid-equilibrium treatment
[Cha1968](@cite). What makes this practical is Julia's metaprogramming: the
whole derivation runs **at compile time** inside `@generated` functions. The
first time [`rate_equation`](@ref) is evaluated for a given mechanism, the
package performs the full symbolic derivation once and writes the resulting
closed-form expression directly into the compiled method. Every later call is
just that expression — a plain arithmetic formula with no symbolic machinery
left behind, no heap allocations, and a runtime on the order of **100 ns** on a
modern processor. The fitter can therefore evaluate the rate equation millions
of times per cross-validation fold and never pay the symbolic cost again.

## Overview of the algorithm

The mechanism graph is partitioned into **rapid-equilibrium (RE) segments**
— sets of enzyme forms linked by RE steps — connected to each other by SS
steps.

**Within each RE segment**, the first form in the mechanism's form order is the
**segment root**, so free E, the first form, roots its own segment. Every other form in
the segment gets an alpha factor, its abundance relative to the root. The law in lowest
terms does not depend on which form roots a segment. The readable
`1 + S / K_ES_to_E_S + P / K_EP_to_E_P` denominator comes from the printer, which divides
the law by the denominator's constant term ([Normalization](@ref)), here free E's weight.


**Across segments**, the SS steps form an inter-segment rate matrix. The denominator
sums, over segments, each segment's sigma (the sum of its forms' alpha factors) times its
cofactor of that matrix, the sum of the spanning trees of the segment graph directed
toward it. The numerator sums each SS step's forward flux times the trees toward its
source segment, less its reverse flux times the trees toward its target, weighted by the
copies of the first substrate the step moves.


## Keeping the equation finite at zero concentration

The derived rate equation is **division-free in the concentrations**: no
metabolite concentration ever appears in a denominator, so the equation contains
no `1/[S]`-style terms anywhere.
The printed law keeps that property: no printed term divides by a concentration
([Normalization](@ref)).

This matters for fitting real data. Datasets often include measurements where a
substrate or product concentration is exactly zero — a saturation curve that
starts at zero substrate, or an initial-rate measurement taken before any
product has formed. An equation with a `1/[S]` term would diverge at those
points and the fit would fail. Because the rate equation has no such terms, it
stays finite at every concentration, including zero, and those data points are
fit like any other.



