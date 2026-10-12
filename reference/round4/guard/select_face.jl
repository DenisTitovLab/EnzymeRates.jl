# ABOUTME: Writes the keys of accepted records whose law vanishes (0/0) on a face, for the
# ABOUTME: package-side evaluation of rate_equation on that face.
using Serialization
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"
out = Tuple{String, String, Symbol}[]
for mode in ("full", "hand")
    for s in deserialize(joinpath(R, "data", "slice_$mode.jls"))
        haskey(s, :law) && haskey(s.law, :reduced) || continue
        t = s.law[:reduced]
        t.dfwd && push!(out, (mode, s.key, :fwd))
        t.drev && push!(out, (mode, s.key, :rev))
    end
end
serialize(joinpath(R, "data", "face_keys.jls"), out)
println(length(out), " faces")
