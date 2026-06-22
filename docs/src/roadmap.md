# Roadmap

The current scope of EnzymeRates.jl covers the three pillars described in
[Getting Started](@ref): deriving rate equations, fitting them to data, and
identifying the best mechanism by cross-validation. The directions below are
under consideration for future development. They are not commitments and carry
no planned schedule.

## KNF allostery model

The existing MWC (concerted, A/I) model treats all subunits as switching
between two states in concert. A Koshland–Némethy–Filmer (sequential) model
would allow subunits to change conformation individually, capturing a broader
class of cooperative behavior.

## V-type allosteric mechanisms

The enumeration cannot currently reach a purely V-type allosteric mechanism —
one where only the catalytic step differs between conformations, so the inactive
state binds substrate but cannot turn it over. The move that makes a mechanism
allosteric can set the catalytic step to `:OnlyA`, but without a regulator the
conformational equilibrium `L` only rescales the rate and is non-identifiable,
so the search prunes that intermediate before a regulator is added (see
[The enumeration engine](@ref)). A planned `+2` move would add an `:OnlyA`
catalytic step together with a regulator in one step, so the regulator makes `L`
identifiable and the V-type mechanism becomes reachable.

## Parameter identifiability

The search selects the best model from data but does not yet report which
parameters the data actually constrain. A planned identifiability analysis
would flag individually unidentifiable parameters and combinations of parameters
that the data cannot separate, giving practitioners a clearer picture of what
the experiment resolves.

## Iso mechanisms

The enumeration engine generates catalytic topologies with isomerization steps,
but certain iso mechanism families — in particular iso mechanisms with multiple
independent isomerization sequences — are not yet fully enumerated. Broader
support would extend the reach of `identify_rate_equation` to enzymes where
isomerization is known to be rate-limiting.

## Plotting

There is no built-in visualization. A planned plotting module would overlay
fitted rate equations on data, making it straightforward to inspect fit quality
and compare the best and runner-up mechanisms visually.

## Outlier dataset identification

The cross-validation framework already computes per-fold losses. A planned
extension would use those fold losses to flag measurement groups that the
selected mechanism fits poorly — surfacing suspect datasets or experimental
conditions that violate the assumed mechanism.
