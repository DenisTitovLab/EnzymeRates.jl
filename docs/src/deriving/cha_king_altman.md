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

**Within each RE segment**, the package picks the form with the fewest bound
metabolites as the **segment root** (the free-enzyme reference form), with
ties broken toward no covalent residual and then by name.
Every other form in the segment gets an alpha factor — its relative abundance
expressed as a ratio to the root.
Referencing to the free-enzyme form is what produces the readable
`1 + S / K_S_E + P / K_P_E` denominator in the final equation.


**Across segments**, the SS steps form an inter-segment rate matrix.
The Cha numerator and denominator come from cofactor determinants of that
matrix, with each diagonal entry weighted by the sigma sum of alpha factors
for its segment.


## Keeping the equation finite at zero concentration

The derived rate equation is **division-free in the concentrations**: no
metabolite concentration ever appears in a denominator, so the equation contains
no `1/[S]`-style terms anywhere.

This matters for fitting real data. Datasets often include measurements where a
substrate or product concentration is exactly zero — a saturation curve that
starts at zero substrate, or an initial-rate measurement taken before any
product has formed. An equation with a `1/[S]` term would diverge at those
points and the fit would fail. Because the rate equation has no such terms, it
stays finite at every concentration, including zero, and those data points are
fit like any other.



