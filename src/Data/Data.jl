# ABOUTME: The Data module: RateData, the validated rate table with its Keq and reference
# ABOUTME: anchors, held by FittingProblem and IdentifyRateEquationProblem.
module Data

using Tables: Tables

const API = (:RateData, :check_anchors)

include("rate_data.jl")

end
