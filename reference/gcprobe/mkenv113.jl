# ABOUTME: Creates a Julia 1.13 environment in gcprobe/env113 that develops the joint-solve
# ABOUTME: worktree by path, so 1.13 probes never rewrite the worktree's own manifest.
using Pkg
Pkg.activate(joinpath(@__DIR__, "env113"))
Pkg.develop(path = ARGS[1])
Pkg.instantiate()
using EnzymeRates
println("ok ", pkgversion(EnzymeRates))
