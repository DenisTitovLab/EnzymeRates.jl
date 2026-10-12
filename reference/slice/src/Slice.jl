# ABOUTME: Vertical-slice prototype of the graph-based rate-law derivation (module Slice):
# ABOUTME: loads the core files, then each later component file present, in pipeline order.
module Slice

for file in ("types.jl", "names.jl", "poly.jl", "convert.jl")
    include(file)
end
for file in ("thermo.jl", "api.jl", "engine.jl", "normalize.jl", "print.jl", "kcat.jl",
             "allosteric.jl", "admission.jl", "print_allo.jl", "codegen.jl", "enum.jl",
             "key.jl", "pipeline.jl")
    isfile(joinpath(@__DIR__, file)) && include(file)
end

end # module
