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
include("symbolic_poly.jl")
include("rate_equation_constraints.jl")
include("rate_equation_derivation.jl")

end
