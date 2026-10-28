# ABOUTME: The Fitting module: FittingProblem, the loss, multi-start fits with
# ABOUTME: Optimization.jl and rate constants scaled to a reference scale.
module Fitting

using ..Mechanisms: Mechanisms
using ..Derivation: Derivation
using ..Data: Data
using FunctionWrappers: FunctionWrapper
using Optimization: Optimization

const API = (:FittingProblem, :fit_rate_equation, :test_rate_equation)

include("fitting.jl")

end
