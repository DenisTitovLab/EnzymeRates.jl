module EnzymeRates

# Types
export AbstractEnzymeReaction, EnzymeReaction
export AbstractEnzymeMechanism, EnzymeMechanism

# DSL
export @enzyme_reaction, @enzyme_mechanism

# Core API
export rate_equation, rate_equation_string, parameters, metabolites

# Identifiability
export is_identifiable, structural_identifiability_deficit

include("types.jl")
include("dsl.jl")
include("sym_poly_for_rate_eq_derivation.jl")
include("rate_eq_derivation.jl")
include("rate_eq_rewriting.jl")

end
