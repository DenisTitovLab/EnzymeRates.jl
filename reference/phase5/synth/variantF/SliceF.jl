# ABOUTME: Law-type variant F of the slice as module SliceF: this directory's patched files
# ABOUTME: (made by ../mkF.py) in place of the slice's originals.
module SliceF
for file in ("types.jl", "names.jl", "poly.jl", "convert.jl", "thermo.jl", "api.jl", "engine.jl", "normalize.jl", "print.jl", "allosteric.jl", "admission.jl", "print_allo.jl", "codegen.jl", "enum.jl", "key.jl")
    local p = joinpath(@__DIR__, file)
    include(isfile(p) ? p : joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/src", file))
end
end # module
