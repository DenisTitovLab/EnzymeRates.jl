module EnzymeRates

export AbstractEnzymeReaction, EnzymeReaction
export AbstractEnzymeMechanism, EnzymeMechanism
export @enzyme_reaction, @mechanism
export substrates, products, regulators
export enzyme_forms, metabolites, n_states, graph, stoich_matrix, reactions, n_steps

# Mode types for rate equations
export RateEquationMode, RawMode, HaldaneWegscheiderMode, IdentifiableHaldaneWegscheiderMode
export Raw, HaldaneWegscheider, IdentifiableHaldaneWegscheider

# Core API (mode-dispatched)
export parameters, rate_equation, rate_equation_string

# Legacy/helper functions
export all_parameters, independent_parameters, dependent_parameters
export constraint_strings
export is_identifiable, non_identifiable_directions, identifiable_combinations

include("types.jl")
include("expr_utils.jl")
include("dsl.jl")
include("accessors.jl")
include("haldane_wegscheider_conditions.jl")
include("qssa.jl")
include("identifiability.jl")

end
