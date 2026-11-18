# ABOUTME: The Plotting module: plot_rate_equation_on_data and the panels it prepares for
# ABOUTME: the Makie extension to draw.
"""Draws a rate equation over measured rates."""
module Plotting

using ..Mechanisms: Mechanisms
using ..Derivation: Derivation
using ..Data: Data
using Statistics: mean
using Tables: Tables

public plot_rate_equation_on_data

include("plot_rate_equation_on_data.jl")

end
