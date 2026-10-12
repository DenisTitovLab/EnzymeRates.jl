# ABOUTME: Lists why the enumerator children without a law failed to derive (error messages).
include(joinpath(@__DIR__, "lib.jl"))
parents = [c.mech for set in [PLAIN_ENUM; ALLO_ENUM]
           for c in deserialize(joinpath(SP, "slice/core/data/$set.jls"))]
kids(m) = unique!([S.flip_children(m); S.split_children(m)])
lvl1 = unique!(reduce(vcat, (kids(m) for m in parents); init = S.Mechanism[]))
sample = lvl1[randperm(Xoshiro(11), length(lvl1))[1:600]]
lvl2 = unique!(reduce(vcat, (kids(c) for c in sample); init = S.Mechanism[]))
test = unique!([lvl1[randperm(Xoshiro(13), length(lvl1))[1:2000]]; lvl2])
cnt = Dict{String, Int}()
for c in test
    lt = law_test(c)
    haskey(lt, :law) && continue
    k = first(string(get(lt, :reduced_err, get(lt, :bare_err, "?"))), 90)
    cnt[k] = get(cnt, k, 0) + 1
end
foreach(println, cnt)
