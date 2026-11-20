# Roadmap

The current scope of EnzymeRates.jl covers the three pillars described in
[Getting Started](@ref): deriving rate equations, fitting them to data, and
identifying the best mechanism by cross-validation. The directions below are
under consideration for future development.

## Support more mechanisms

EnzymeRates derives and identifies a broad space of mechanisms, but several
families are not yet covered.

**KNF (sequential) allostery.** The current MWC model treats all subunits as
switching conformation in concert. A Koshland–Némethy–Filmer (sequential) model
would let subunits change conformation individually, capturing a complementary
class of cooperative behavior.

**Iso mechanisms.** In an iso mechanism the enzyme takes on a different
conformation as it converts from the substrate-bound to the product-bound
species, and that conformation relaxes back to the original only through a
separate, slow steady-state isomerization step ([Iso mechanisms](@ref) works
through an example). The derivation already handles these, but the enumeration
does not yet generate iso mechanisms, so `identify_rate_equation` cannot reach
them.

**Allosteric ping-pong that stops both half-reactions.** An allosteric ping-pong may stop
both half-reactions in the inactive conformation when its covalent intermediate can lie
higher in energy there: each half's step into it vanishes while its step out keeps its
rate, so the intermediate empties ([How an excluded step is switched off](@ref)). Where
the intermediate empties only by binding the second substrate, no choice of energies gives
the limit, and `@allosteric_mechanism` rejects the mechanism. Supporting it needs the conformations to
interconvert at the covalent intermediate as well, which will come with iso mechanisms.

## Support more measurement types

Fitting and identification currently work from reaction-rate (turnover) data. A
planned extension would support binding and affinity observables as well: the
fractional binding of one or more metabolites to the enzyme, an apparent `Km`,
or a metabolite dissociation constant `Kd`. A mechanism could then be derived
against and constrained by binding measurements, not rate alone, and report
these quantities directly.

## Parameter estimation and identifiability

Beyond selecting the best mechanism, a planned analysis would characterize its
parameters. Bootstrap resampling of the data would estimate a distribution — and
hence a confidence interval — for each fitted parameter rather than a single
point value. An accompanying identifiability analysis would flag parameters,
individually or in combination, that the data cannot constrain, giving a clear
picture of what the experiment actually resolves.

## Outlier dataset identification

The cross-validation framework already computes per-fold losses. A planned
extension would use those fold losses to flag measurement groups that the
selected mechanism fits poorly — surfacing suspect datasets or experimental
conditions that violate the assumed mechanism.
