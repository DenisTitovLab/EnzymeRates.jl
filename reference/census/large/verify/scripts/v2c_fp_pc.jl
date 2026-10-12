# ABOUTME: The fingerprint-only matches of pyruvate carboxylase tested for a shared law by fitting
# ABOUTME: both ways (the part of v2b_dz.jl that stopped on a non-finite fit Jacobian).
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "vfit.jl"))
name = :pyruvate_carboxylase_rxn
data = load_run(name)
byid = Dict(r[:id] => r for recs in data.rounds for r in recs)
rng = Xoshiro(20261011 + 16)
sp = deserialize(joinpath(VER, "data", "v2_splits_$(name).jls"))
sample = sp[sort!(randperm(rng, length(sp))[1:min(25, length(sp))])]
tally = Dict{String, Int}()
for (id, others) in sample
    r, o = byid[id], byid[others[1]]
    a, b = both_ways(r[:smech], o[:smech], rng)
    v = a < 1e-8 && b < 1e-8 ? "same law" : "different"
    tally[v] = get(tally, v, 0) + 1
    @printf("   fp-only %s vs %s (slice %d/%d vs %d/%d): fits %.1e %.1e -> %s\n", id, o[:id],
            r[:rank], r[:n_slice], o[:rank], o[:n_slice], a, b, v)
end
println("(b) $name fingerprint-only matches: $(length(sp)); sampled $(length(sample)): ", tally)
