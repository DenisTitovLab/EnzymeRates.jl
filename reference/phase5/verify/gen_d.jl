# ABOUTME: Verifier run D (package): rebuilds the phantom population of gen_a.jl and lists every
# ABOUTME: parent edge of chosen gids, so a phantom's lineage is not limited to its first parent.
using EnzymeRates, Random, Printf
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SP, "phase45/research5/verify/pkgside.jl"))
const WATCH = parse.(Int, ARGS)
rxn2 = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end
mechs, gid, info, edges = Any[], Dict{Any, Int}(), Any[], Tuple{Int, Int, Symbol}[]
function reg!(m, round, par, mv)
    g = get(gid, m, 0)
    g == 0 || return g
    push!(mechs, m); g = gid[m] = length(mechs)
    cx = cxA(m); np = ER._independent_param_count(m)
    push!(info, (; round, ok = cx <= 337 && np <= 20))
    g
end
okg(r) = [g for (g, x) in enumerate(info) if x.round == r && x.ok]
seeds = unique!(collect(ER.init_mechanisms(rxn2)))
base, _ = ER._base_tier(seeds, rxn2)
foreach(m -> reg!(m, 0, 0, :seed), base)
for p in okg(0), (c, mv) in children_labeled(mechs[p], rxn2)
    push!(edges, (p, reg!(c, 1, p, mv), mv))
end
pool1 = okg(1)
theirs = pool1[sort!(randperm(Xoshiro(60), length(pool1))[1:1000])]
rest = setdiff(pool1, theirs)
mine = rest[sort!(randperm(Xoshiro(424242), length(rest))[1:300])]
for p in mine, (c, mv) in children_labeled(mechs[p], rxn2)
    push!(edges, (p, reg!(c, 2, p, mv), mv))
end
println("records ", length(mechs), ", edges ", length(edges))
first_parent(g) = (i = findfirst(e -> e[2] == g, edges); i === nothing ? 0 : edges[i][1])
for g in WATCH
    ps = [(e[1], e[3]) for e in edges if e[2] == g]
    println("gid $g: parents ", ps, "; their first-parent chains ",
            [(p, first_parent(p)) for (p, _) in ps])
end
