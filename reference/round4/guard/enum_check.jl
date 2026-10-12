# ABOUTME: Applies the guard candidates to the enumerator's own children (RE→SS flips and
# ABOUTME: kinetic-group splits of every enumerated parent, and a second level from a sample).
include(joinpath(@__DIR__, "lib.jl"))
const T0 = time()
parents = [c.mech for set in [PLAIN_ENUM; ALLO_ENUM]
           for c in deserialize(joinpath(SP, "slice/core/data/$set.jls"))]
println("parents: ", length(parents))
const PARENT = Dict{S.Mechanism, Tuple{S.Mechanism, Symbol}}()
function children(m)
    out = S.Mechanism[]
    for (mv, f) in ((:flip, S.flip_children), (:split, S.split_children))
        try
            cs = f(m)
            foreach(c -> get!(PARENT, c, (m, mv)), cs)
            append!(out, cs)
        catch e
            println("move error: ", errmsg(e))
        end
    end
    unique!(out)
end
function verdicts(c)
    fwd = face_check(c, c.products)
    rev = face_check(c, c.substrates)
    if fwd === nothing && rev === nothing && S.is_allosteric(c) && any(S.inactive_mask(c))
        im = S.inactive_mask(c)
        fwd, rev = face_check(c, c.products, im), face_check(c, c.substrates, im)
    end
    (bottomless = S.bottomless(c) !== nothing, face = something(fwd, rev, ""),
     re = re_turnover(c) !== nothing,
     admitted = S.admit(c) === nothing)
end
lvl1 = S.Mechanism[]
for (i, m) in enumerate(parents)
    append!(lvl1, children(m))
    i % 1000 == 0 && @printf("%d parents, %d children, %.0f s\n", i, length(lvl1), time() - T0)
end
unique!(lvl1)
rng = Xoshiro(11)
sample = lvl1[randperm(rng, length(lvl1))[1:min(600, length(lvl1))]]
lvl2 = unique!(reduce(vcat, (children(c) for c in sample); init = S.Mechanism[]))
@printf("level 1: %d distinct children; level 2 (from 600 level-1): %d (%.0f s)\n",
        length(lvl1), length(lvl2), time() - T0)
flagged = Tuple{String, S.Mechanism, Any}[]
for (name, cs) in (("level1", lvl1), ("level2", lvl2))
    cnt = Dict{String, Int}()
    for c in cs
        v = verdicts(c)
        k = v.bottomless ? "bottomless" : v.re ? "RE turnover" :
            !isempty(v.face) ? "face: " * split(v.face, " ")[1] : "pass"
        k = k * (v.admitted ? "" : " (not admitted)")
        cnt[k] = get(cnt, k, 0) + 1
        k == "pass" || k == "pass (not admitted)" || push!(flagged, (name, c, v))
    end
    println(name, ": ", sort!(collect(cnt)))
end
println("\nflagged children: ", length(flagged))
for (name, c, v) in flagged[1:min(end, 6)]
    println("  [$name] ", v.face, v.bottomless ? " bottomless" : "", v.re ? " RE turnover" : "")
    foreach(e -> println("      ", S.step_text(c, e), " ::g", c.group[findfirst(==(e), c.edges)],
                         S.is_allosteric(c) ? " " * string(c.tags[c.group[findfirst(==(e), c.edges)]]) : ""), c.edges)
end
serialize(joinpath(R, "data", "enum_flagged.jls"),
          [(level = name, parent = spec(PARENT[c][1]), move = PARENT[c][2], child = spec(c), face = v.face)
           for (name, c, v) in flagged])
# the oracle on flagged children
rng = Xoshiro(12)
oc = Dict{String, Int}()
for (name, c, v) in flagged
    lt = law_test(c)
    haskey(lt, :law) || (oc["no law"] = get(oc, "no law", 0) + 1; continue)
    num = numeric_checks(c, lt[:law], rng)
    k = join(string.([num[f].verdict for f in (:pos, :fwd, :rev)]), "/") * " law " * string(lt[:reduced])
    oc[k] = get(oc, k, 0) + 1
end
println("oracle on flagged: ", oc)
@printf("done %.0f s\n", time() - T0)
