# ABOUTME: Builds this study's environment (census mkenv recipe plus the fit and CSV packages):
# ABOUTME: develops the graph_mechanism_wt worktree by path. Usage: julia mkenv.jl <worktree> <envdir>
using Pkg
Pkg.offline(true)
Pkg.activate(ARGS[2])
Pkg.develop(path = ARGS[1])
Pkg.add(["LinearAlgebra", "Random", "SHA", "Test", "TOML", "Serialization", "Printf",
         "Statistics", "OptimizationCMAEvolutionStrategy", "CSV", "DataFrames"])
Pkg.instantiate()
Pkg.precompile()
using EnzymeRates
println("ok ", pkgversion(EnzymeRates), " from ", pathof(EnzymeRates))
