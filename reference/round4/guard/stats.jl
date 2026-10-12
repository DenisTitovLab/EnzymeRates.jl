# ABOUTME: Size of the face discrepancies where the RE limit depends on the speed pattern:
# ABOUTME: spread between two patterns and the law's distance from the oracle, per record.
using Serialization, Printf, Statistics
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"
sl = deserialize(joinpath(R, "data", "slice_full.jls"))
pat, dev, mis = Float64[], Float64[], Float64[]
nz = 0
for s in sl
    haskey(s, :num) && s.num isa Dict && !haskey(s.num, :err) || continue
    for f in (:fwd, :rev)
        n = s.num[f]
        if n.verdict === :pattern
            o2, o3 = n.o[2], n.o[3]
            push!(pat, abs(o2 - o3) / max(abs(o2), abs(o3)))
            push!(dev, max(abs(n.vlaw - o2), abs(n.vlaw - o3)) / max(abs(o2), abs(o3), abs(n.vlaw)))
        elseif n.verdict === :mismatch
            n.o[2] == 0 && (global nz += 1)
            push!(mis, abs(n.vlaw - n.o[2]) / max(abs(n.vlaw), abs(n.o[2])))
        end
    end
end
qs(v) = @sprintf("n=%d min %.1e  median %.1e  90%% %.1e  max %.1e", length(v), minimum(v), median(v), quantile(v, 0.9), maximum(v))
println("pattern-dependent faces: spread between two speed patterns |o1-o2|/max: ", qs(pat))
println("                         law's largest relative distance from the two:  ", qs(dev))
println("mismatch faces (unique limit ≠ law): |law-o|/max: ", qs(mis), "; oracle exactly 0 in $nz")
