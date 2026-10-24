# ABOUTME: Adds OptimizationCMAEvolutionStrategy to the gcprobe 1.13 environment (offline), for
# ABOUTME: the identify workload's optimizer.
using Pkg; Pkg.activate(joinpath(@__DIR__, "env113")); Pkg.add("OptimizationCMAEvolutionStrategy")
using OptimizationCMAEvolutionStrategy; println("ok")
