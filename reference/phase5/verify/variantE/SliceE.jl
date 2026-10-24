# ABOUTME: Law-type variant E for the verifier: C's Int-only engine and allosteric law with Poly a
# ABOUTME: mutable struct (variantE/poly.jl); other files from the stock slice.
module SliceE
const SRC = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/src"
const FROM = Dict("poly.jl" => "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/verify/variantE", "engine.jl" => "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/law/variants/C", "allosteric.jl" => "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/law/variants/C")
for file in ("types.jl", "names.jl", "poly.jl", "convert.jl", "thermo.jl", "api.jl", "engine.jl",
             "normalize.jl", "print.jl", "allosteric.jl", "admission.jl", "print_allo.jl",
             "codegen.jl", "enum.jl", "key.jl")
    include(joinpath(get(FROM, file, SRC), file))
end
end # module
