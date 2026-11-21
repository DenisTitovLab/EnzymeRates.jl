# ABOUTME: Top-level module; defines the exported names and include order.
# ABOUTME: Exports the 15 names promised to users and wires together src/ source files.
"""Finds the best mechanism for a reaction from measured rates."""
module EnzymeRates

# Types
export EnzymeReaction, Mechanism
export FittingProblem
export IdentifyRateEquationProblem, IdentifyRateEquationResults

# DSL
export @enzyme_reaction, @enzyme_mechanism

# Core API
export compile_rate_equation, rate_equation_string, parameters,
    dependent_parameter_values, metabolites

# Fitting & model selection
export fit_rate_equation, identify_rate_equation

# Plotting
export plot_rate_equation_on_data


using Dates
using Distributed

include("Data/Data.jl")
using .Data: Data
include("Mechanisms/Mechanisms.jl")
using .Mechanisms: Mechanisms, EnzymeReaction, Mechanism, metabolites
include("Derivation/Derivation.jl")
using .Derivation: Derivation, compile_rate_equation, rate_equation_string, parameters,
    dependent_parameter_values
include("DSL/DSL.jl")
using .DSL: DSL, @enzyme_reaction, @enzyme_mechanism
include("Fitting/Fitting.jl")
using .Fitting: Fitting, FittingProblem, fit_rate_equation, test_rate_equation
include("Enumeration/Enumeration.jl")
include("identify_rate_equation.jl")
include("Plotting/Plotting.jl")
using .Plotting: Plotting, plot_rate_equation_on_data

end
