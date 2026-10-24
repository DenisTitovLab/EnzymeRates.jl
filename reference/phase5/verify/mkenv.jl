# ABOUTME: Creates the verifier's environment, developing graph_mechanism_wt by path (census mkenv
# ABOUTME: recipe plus Printf). Usage: julia mkenv.jl <worktree> <envdir>
using Pkg
Pkg.offline(true)
Pkg.activate(ARGS[2])
Pkg.develop(path = ARGS[1])
Pkg.add(["LinearAlgebra", "Random", "SHA", "Test", "TOML", "Printf"])
Pkg.instantiate()
Pkg.precompile()
using EnzymeRates
println("ok ", pkgversion(EnzymeRates), " from ", pathof(EnzymeRates))
