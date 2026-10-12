# ABOUTME: Loads the full Slice.jl (every stage file present) to check that enum.jl loads beside
# ABOUTME: the other stages' files without overwriting their methods; prints the files present.
const SRC = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/src"
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common/Population.jl")
include(joinpath(SRC, "Slice.jl"))
println("files: ", filter(f -> endswith(f, ".jl"), readdir(SRC)))
println("flip_children defined: ", isdefined(Slice, :flip_children))
