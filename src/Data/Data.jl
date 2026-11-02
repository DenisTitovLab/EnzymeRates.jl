# ABOUTME: The Data module: RateData, the validated rate table with its Keq and reference
# ABOUTME: anchors, held by FittingProblem and IdentifyRateEquationProblem.
"""Checks and holds a table of measured rates."""
module Data

using Tables: Tables

const API = (:RateData,)

include("rate_data.jl")

end
