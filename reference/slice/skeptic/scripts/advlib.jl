# ABOUTME: Builds the skeptic's adversarial mechanisms as Pop records (one spec, both routes) and
# ABOUTME: the ε-limit oracle of decision 6 (the :OnlyA steps vanishing along a consistent path).
include(joinpath(@__DIR__, "spotlib.jl"))

"""A modifier copy (`A::Mod`) or a modifier-only name in a step or form."""
struct Mod
    name::Symbol
end

"""A form spec: conformation, bound metabolites (`Symbol` or `Mod`), residual."""
F(conf, bound...; add = Symbol[], sub = Symbol[]) = (conf, Any[bound...], add, sub)

"""A step spec in its physical-forward direction: `f` takes up `c` and gives off `r` to
become `t`; group `g`."""
st(f, t; c = Any[], r = Any[], re = false, g = 1) =
    (f = f, t = t, c = Any[c...], r = Any[r...], re = re, g = g)

"""The Pop record of a spec: species, steps in the written direction, groups by label."""
function build(id, subs, prods, steps; tags = Symbol[], mult = 1, sites = [])
    subs, prods = sort(subs), sort(prods)
    role(x) = x isa Mod || !(x in subs || x in prods) ? :Inh : x in subs ? :S : :P
    nm(x) = x isa Mod ? x.name : x
    pmet(x) = Pop.PMet(nm(x), role(x))
    key(f) = (f[1], sort([(nm(x), role(x)) for x in f[2]]), sort(f[3]), sort(f[4]))
    keys_ = unique([key(x) for s in steps for x in (s.f, s.t)])
    idx = Dict(k => i for (i, k) in enumerate(keys_))
    species = [Pop.PSpecies([Pop.PMet(n, r) for (n, r) in k[2]], k[1], k[3], k[4],
                            Symbol(k[1], join(string.(first.(k[2]))))) for k in keys_]
    psteps = [Pop.PStep(idx[key(s.f)], idx[key(s.t)], pmet.(s.c), pmet.(s.r), s.re)
              for s in steps]
    G = maximum(s.g for s in steps)
    groups = [findall(s -> s.g == g, steps) for g in 1:G]
    inh = sort(unique([nm(x) for s in steps for x in [s.c; s.r; s.f[2]; s.t[2]]
                       if role(x) === :Inh && !(nm(x) in subs || nm(x) in prods)]))
    regs = sort(unique([l for s in sites for l in s.ligands]))
    Pop.PMech(; id, source = id, substrates = subs, products = prods, inhibitors = inh,
              regulators = regs, metabolites = [subs; prods; sort(unique([inh; regs]))],
              species, steps = psteps, groups, allosteric = !isempty(tags),
              cat_multiplicity = mult, group_tags = tags,
              sites = [Pop.PSite(s.ligands, s.multiplicity, s.states) for s in sites])
end

"""
    eps_limit(p, tags, kA, kI, L, concs, ε) -> BigFloat

Decision 6's limit for multiplicity 1 without sites, as one network of both conformations
joined at the flip species by `[I]/[A] = L`: the inactive `:OnlyA` steps move along the
potential of the components that the other steps join (each `:OnlyA` plain binding rising
toward its bound side), so a binding's affinity vanishes like `ε^Δφ` and every other `:OnlyA`
step is also slowed by `ε^(1 + |Δφ|)`; thermodynamics holds for every `ε`.
"""
function eps_limit(p, tags, kA, kI, L, concs, ε)
    g = Pop.group_of_step(p)
    onlyA(j) = tags[g[j]] === :OnlyA
    n = length(p.species)
    comp = collect(1:n)
    top(i) = (while comp[i] != i; i = comp[i]; end; i)
    for (j, s) in enumerate(p.steps)
        onlyA(j) || (comp[top(s.from)] = top(s.to))
    end
    comp = top.(1:n)
    plain(s) = Pop.is_binding(s) && p.species[s.from].res_added == p.species[s.to].res_added &&
        p.species[s.from].res_subtracted == p.species[s.to].res_subtracted &&
        sort([(x.name, x.role) for x in [p.species[s.from].bound; s.consumed]]) ==
        sort([(x.name, x.role) for x in p.species[s.to].bound])
    arcs = [(comp[s.from], comp[s.to]) for (j, s) in enumerate(p.steps) if onlyA(j) && plain(s)]
    phi = zeros(Int, n)                       # longest-path rank along the arcs
    for _ in 1:n, (u, v) in arcs
        phi[v] = max(phi[v], phi[u] + 1)
    end
    φ(i) = phi[comp[i]]
    kIe, reI = Vector{BigFloat}[], Bool[]
    for (j, s) in enumerate(p.steps)
        if !onlyA(j)
            push!(kIe, kI[j]); push!(reI, s.re)
            continue
        end
        ka = s.re ? (Pop.is_binding(s) ? 1 / kI[j][1] : kI[j][1]) : kI[j][1] / kI[j][2]
        dphi = φ(s.to) - φ(s.from)
        kanew = ka * ε^dphi
        if plain(s)
            if s.re
                push!(kIe, [1 / kanew]); push!(reI, true)
            else
                push!(kIe, [kI[j][2] * kanew, kI[j][2]]); push!(reI, false)
            end
        else
            kr = (s.re ? big(1) : kI[j][2]) * ε^(1 + abs(dphi))
            push!(kIe, [kr * kanew, kr]); push!(reI, false)
        end
    end
    flips = [i for (i, sp) in enumerate(p.species) if Pop.is_free(sp)]
    flip = flips[argmin([String(p.species[i].conformation) for i in flips])]
    Λ = big(10)^60
    M = zeros(BigFloat, 2n, 2n)
    x1 = p.substrates[1]
    fl = Tuple{Int, Int, BigFloat, BigFloat, Int}[]
    for (off, k, re) in ((0, kA, [s.re for s in p.steps]), (n, kIe, reI)),
        (j, s) in enumerate(p.steps)
        a, b = s.from + off, s.to + off
        f, r = re[j] ? (Pop.is_binding(s) ? (Λ, Λ * k[j][1]) : (Λ, Λ / k[j][1])) :
               (k[j][1], k[j][2])
        rf = f * prod((BigFloat(concs[x.name]) for x in s.consumed); init = big(1))
        rr = r * prod((BigFloat(concs[x.name]) for x in s.released); init = big(1))
        M[b, a] += rf; M[a, a] -= rf; M[a, b] += rr; M[b, b] -= rr
        w = count(x -> x.name === x1, s.consumed) - count(x -> x.name === x1, s.released)
        push!(fl, (a, b, rf, rr, w))
    end
    a, b = flip, flip + n
    M[b, a] += Λ; M[a, a] -= Λ; M[a, b] += Λ / L; M[b, b] -= Λ / L
    M[2n, :] .= 1
    rhs = zeros(BigFloat, 2n); rhs[2n] = 1
    x = M \ rhs
    sum((w * (rf * x[a] - rr * x[b]) for (a, b, rf, rr, w) in fl if w != 0); init = big(0))
end
