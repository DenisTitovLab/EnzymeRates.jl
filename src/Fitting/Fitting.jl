# ABOUTME: The Fitting module: FittingProblem, the loss, multi-start fits with
# ABOUTME: Optimization.jl and rate constants rescaled to reference turnovers.
module Fitting

using ..Mechanisms: Mechanisms
using ..Derivation: Derivation
using ..Data: Data
using FunctionWrappers: FunctionWrapper
using Optimization: Optimization

const API = (:FittingProblem, :cross_validation_loss, :fit_rate_equation,
    :rescale_parameter_values)

include("fitting.jl")

end
