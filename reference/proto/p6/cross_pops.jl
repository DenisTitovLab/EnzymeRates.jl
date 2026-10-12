# ABOUTME: Counts the P6 equality classes that span populations, by population pair.
# ABOUTME: Such classes are the same graph on different (or identical) reactions.
include(joinpath(@__DIR__, "harness.jl"))
ms = _testhelper_load_all()
pop(id) = split(id, "/")[1]
cls = Dict{GM.GMech, Vector{Int}}()
for (i, m) in enumerate(ms); push!(get!(cls, _testhelper_gm(_testhelper_toraw(m)), Int[]), i); end
pairs = Dict{String, Int}()
for c in values(cls)
    ps = sort(unique(pop.(getfield.(ms[c], :id))))
    length(ps) > 1 && (k = join(ps, " + "); pairs[k] = get(pairs, k, 0) + 1)
end
for (k, v) in sort!(collect(pairs); by = last, rev = true); println(rpad(k, 40), v); end
