# ABOUTME: Law-type variant D for the verifier: B's Dict-alias Poly (and its key, codegen, print
# ABOUTME: reads) with C's Int-only engine and allosteric law; other files from the stock slice.
module SliceD
const SRC = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/src"
const FROM = Dict("poly.jl" => "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/law/variants/B", "key.jl" => "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/law/variants/B", "codegen.jl" => "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/law/variants/B",
                  "print.jl" => "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/law/variants/B", "engine.jl" => "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/law/variants/C", "allosteric.jl" => "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/law/variants/C")
for file in ("types.jl", "names.jl", "poly.jl", "convert.jl", "thermo.jl", "api.jl", "engine.jl",
             "normalize.jl", "print.jl", "allosteric.jl", "admission.jl", "print_allo.jl",
             "codegen.jl", "enum.jl", "key.jl")
    include(joinpath(get(FROM, file, SRC), file))
end
end # module
