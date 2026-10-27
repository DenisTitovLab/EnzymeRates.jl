# ABOUTME: The Fitting module: FittingProblem, the loss, multi-start fits with Optimization.jl
# ABOUTME: and rate constants rescaled to reference turnovers.
module Fitting

using ..Mechanisms: Mechanisms
using ..Derivation: Derivation
using FunctionWrappers: FunctionWrapper
using Optimization: Optimization
using Tables: Tables

const API = (:FittingProblem, :_rate_table, :_reference, :fit_rate_equation, :loss!,
    :rescale_parameter_values)

include("fitting.jl")

end
