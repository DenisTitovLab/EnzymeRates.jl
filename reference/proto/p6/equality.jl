# ABOUTME: P6 equality checks on every population: the == partition against the package's,
# ABOUTME: equivalent variants (equal, same hash) and single mutations (unequal, new hash).
include(joinpath(@__DIR__, "harness.jl"))
using Printf

ms = _testhelper_load_all()
pop(id) = split(id, "/")[1]
gms = [_testhelper_gm(_testhelper_toraw(m)) for m in ms]
pkeys = [_testhelper_pkgkey(m) for m in ms]

# ── Partition ────────────────────────────────────────────────────────────────────
mine = Dict{GM.GMech, Vector{Int}}()
for (i, g) in enumerate(gms); push!(get!(mine, g, Int[]), i); end
pk = Dict{String, Vector{Int}}()
for (i, k) in enumerate(pkeys); push!(get!(pk, k, Int[]), i); end
false_merges = [c for c in values(mine) if length(unique(pkeys[c])) > 1]
false_splits = [c for c in values(pk) if length(unique(gms[c])) > 1]
cross = [c for c in values(mine) if length(unique(pop.(getfield.(ms[c], :id)))) > 1]
within = [c for c in values(pk) if length(c) > length(unique(pop.(getfield.(ms[c], :id))))]
@printf("mechanisms %d: %d GM classes, %d package classes (reaction ignored)\n",
        length(ms), length(mine), length(pk))
@printf("false merges %d, false splits %d\n", length(false_merges), length(false_splits))
@printf("classes spanning populations (same graph, different reaction): %d, %s\n", length(cross),
        join([join(getfield.(ms[c], :id), "=") for c in first(cross, 6)], "; "))
@printf("package-equal duplicates inside one population: %d\n", length(within))
@printf("distinct hash values over the %d classes: %d\n", length(mine),
        length(unique(hash.(collect(keys(mine))))))
phys = Dict{GM.GMech, Vector{Int}}()
for (i, g) in enumerate(gms); push!(get!(phys, GM.physical(g), Int[]), i); end
@printf("orientation rule 2 (physical-forward storage) gives the same partition: %s\n",
        Set(Set.(values(phys))) == Set(Set.(values(mine))))

# ── Equivalent variants ──────────────────────────────────────────────────────────
nvar = 0; neq_fail = String[]; nhash_fail = String[]
for (i, m) in enumerate(ms)
    n = length(m.steps)
    for v in 1:6
        rng = Xoshiro(1000i + v)
        flip = v == 1 ? trues(n) : v == 2 ? falses(n) : bitrand(rng, n)
        gv = _testhelper_gm(_testhelper_toraw(m; rng = v == 1 ? nothing : rng, flip))
        global nvar += 1
        gv == gms[i] || push!(neq_fail, "$(m.id)#$v")
        hash(gv) == hash(gms[i]) || push!(nhash_fail, "$(m.id)#$v")
    end
end
@printf("variants %d (all steps reversed; permuted; permuted + random reversals x4): unequal %d, hash differs %d %s\n",
        nvar, length(neq_fail), length(nhash_fail), join(first(neq_fail, 5), ", "))

# ── Single mutations of a shuffled input ─────────────────────────────────────────
function mutate(raw, kind, m)
    subs, prods, steps, group, kw = deepcopy(raw)
    tags, mult, sites = collect(kw.tags), kw.multiplicity, collect(kw.sites)
    G = maximum(group)
    relabel_out!(l) = (for j in eachindex(group); group[j] > l && (group[j] -= 1); end;
                       isempty(tags) || deleteat!(tags, l))
    if kind === :re                       # flip RE/SS of one whole group
        for j in eachindex(steps)
            group[j] == 1 || continue
            s = steps[j]; steps[j] = GM.RawStep(s.from, s.to, s.consumed, s.released, !s.re)
        end
    elseif kind === :split                # one step of a multi-step group leaves it
        l = findfirst(l -> count(==(l), group) > 1, 1:G)
        l === nothing && return nothing
        group[findfirst(==(l), group)] = G + 1
        isempty(tags) || push!(tags, tags[l])
    elseif kind === :merge                # two groups of one kind become one
        kindof(s) = (s.re, Set([(sort(s.consumed), sort(s.released)), (sort(s.released), sort(s.consumed))]))
        pair = nothing
        for l1 in 1:G, l2 in l1+1:G
            k1 = [kindof(steps[j]) for j in eachindex(steps) if group[j] == l1]
            k2 = [kindof(steps[j]) for j in eachindex(steps) if group[j] == l2]
            length(unique([k1; k2])) == 1 && (pair = (l1, l2); break)
        end
        pair === nothing && return nothing
        for j in eachindex(group); group[j] == pair[2] && (group[j] = pair[1]); end
        relabel_out!(pair[2])
    elseif kind === :tag
        isempty(tags) && return nothing
        tags[1] = GM.TAGS[mod1(findfirst(==(tags[1]), GM.TAGS) + 1, 3)]
    elseif kind === :mult
        isempty(tags) && return nothing
        mult += 1
    elseif kind === :site
        isempty(sites) && return nothing
        sites[1] = GM.Site(sites[1].ligands, sites[1].states, sites[1].multiplicity + 1)
    elseif kind === :inh                  # a reactant binds as an inhibitor copy throughout
        x = m.substrates[1]
        sw(v) = [y == GM.Met(x, false) ? GM.Met(x, true) : y for y in v]
        swf(f) = GM.RawForm(f.conf, sw(f.bound), f.added, f.subtracted)
        any(y -> y == GM.Met(x, false), [s.from.bound for s in steps] |> Iterators.flatten) ||
            return nothing
        steps = [GM.RawStep(swf(s.from), swf(s.to), sw(s.consumed), sw(s.released), s.re) for s in steps]
    elseif kind === :conf                 # rename a conformation
        rn(f) = GM.RawForm(f.conf === :E ? :Z : f.conf, f.bound, f.added, f.subtracted)
        steps = [GM.RawStep(rn(s.from), rn(s.to), s.consumed, s.released, s.re) for s in steps]
    elseif kind === :roles                # swap a substrate and a product
        subs[1], prods[1] = prods[1], subs[1]
    elseif kind === :drop                 # drop one step
        l = group[end]
        pop!(steps); pop!(group)
        l in group || relabel_out!(l)
        isempty(steps) && return nothing
    end
    (subs, prods, steps, group, (; tags, multiplicity = mult, sites))
end
kinds = (:re, :split, :merge, :tag, :mult, :site, :inh, :conf, :roles, :drop)
stats = Dict(k => [0, 0, 0, 0] for k in kinds)    # applied, rejected, equal, same hash
for (i, m) in enumerate(ms)
    raw = _testhelper_toraw(m; rng = Xoshiro(77i), flip = bitrand(Xoshiro(78i), length(m.steps)))
    for k in kinds
        r = mutate(raw, k, m)
        r === nothing && continue
        st = stats[k]; st[1] += 1
        gmut = try _testhelper_gm(r) catch; st[2] += 1; continue end
        gmut == gms[i] && (st[3] += 1)
        hash(gmut) == hash(gms[i]) && (st[4] += 1)
    end
end
println("mutation    applied  rejected-by-ctor  equal(BAD)  same-hash(BAD)")
for k in kinds
    st = stats[k]
    @printf("%-10s %8d %17d %11d %15d\n", k, st...)
end
