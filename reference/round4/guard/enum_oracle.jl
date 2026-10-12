# ABOUTME: The RE-limit oracle on unflagged enumerator children (a level-1 sample and every
# ABOUTME: level-2 child): measures what the face check and the law-based test miss there.
include(joinpath(@__DIR__, "lib.jl"))
const T0 = time()
parents = [c.mech for set in [PLAIN_ENUM; ALLO_ENUM]
           for c in deserialize(joinpath(SP, "slice/core/data/$set.jls"))]
kids(m) = unique!([S.flip_children(m); S.split_children(m)])
lvl1 = unique!(reduce(vcat, (kids(m) for m in parents); init = S.Mechanism[]))
rng = Xoshiro(11)
sample = lvl1[randperm(rng, length(lvl1))[1:600]]
lvl2 = unique!(reduce(vcat, (kids(c) for c in sample); init = S.Mechanism[]))
test = unique!([lvl1[randperm(Xoshiro(13), length(lvl1))[1:2000]]; lvl2])
println("children tested: ", length(test))
cnt = Dict{String, Int}()
rng = Xoshiro(14)
for (i, c) in enumerate(test)
    fc = something(face_check(c, c.products), face_check(c, c.substrates), "")
    lt = law_test(c)
    k = if !haskey(lt, :law)
        "no law"
    else
        num = numeric_checks(c, lt[:law], rng)
        vs = [num[f].verdict for f in (:pos, :fwd, :rev)]
        o = all(==(:agree), vs) ? "oracle agree" : "oracle " * join(string.(vs), "/")
        l = (t = lt[:reduced]; t.b ? "N≡0" : (t.dfwd || t.drev) ? "D=0 face" : "law ok")
        "$(isempty(fc) ? "face pass" : "face flag") | $l | $o | $(S.admit(c) === nothing ? "admitted" : "not admitted")"
    end
    cnt[k] = get(cnt, k, 0) + 1
    i % 1000 == 0 && @printf("%d %.0f s\n", i, time() - T0)
end
foreach(p -> println(p), sort!(collect(cnt); by = last, rev = true))
@printf("done %.0f s\n", time() - T0)
