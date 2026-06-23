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

**V-type allosteric mechanisms.** A purely V-type allosteric mechanism — one
where only the catalytic step differs between conformations, so the inactive
state binds substrate but cannot turn it over — is not currently reachable by
the search: the intermediate that introduces the conformational equilibrium `L`
is non-identifiable without a regulator (see [The enumeration engine](@ref)). A
planned `+2` move would add the `:OnlyA` catalytic step and a regulator together,
so `L` becomes identifiable and the V-type mechanism becomes reachable.

**Iso mechanisms.** In an iso mechanism the enzyme takes on a different
conformation as it converts from the substrate-bound to the product-bound
species, and that conformation relaxes back to the original only through a
separate, slow steady-state isomerization step ([Iso mechanisms](@ref) works
through an example). The derivation already handles these, but the enumeration
does not yet generate every iso family, so `identify_rate_equation` cannot reach
all of them.

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

## Plotting

There is no built-in visualization. A planned plotting module would overlay
fitted rate equations on data, making it straightforward to inspect fit quality
and compare the best and runner-up mechanisms visually.

## Outlier dataset identification

The cross-validation framework already computes per-fold losses. A planned
extension would use those fold losses to flag measurement groups that the
selected mechanism fits poorly — surfacing suspect datasets or experimental
conditions that violate the assumed mechanism.
