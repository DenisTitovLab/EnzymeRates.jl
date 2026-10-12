# ABOUTME: Mutants of enumerated mechanisms for the guard study: one or two kinetic groups'
# ABOUTME: RE/SS flags flipped, the all-RE version, one catalytic step deleted; with provenance.

const KINDS = (:to_re, :to_ss, :flip2, :all_re, :del_cat)

"""Mutant `kind` of `m` as `(desc, spec, mech or nothing, ctor_err)`, or `nothing` when the
kind does not apply (no SS group to flip, already all RE, ...). Group flips keep each group
uniform (its members share their flag)."""
function mutate(m, kind, rng)
    (subs, prods, steps0, group0), kw = S.raw(m)
    steps, group, tags = copy(steps0), copy(group0), copy(kw.tags)
    G = maximum(group)
    flag(g) = steps[findfirst(==(g), group)].re
    function setre!(g, re)
        for j in findall(==(g), group)
            s = steps[j]
            steps[j] = S.RawStep(s.from, s.to, s.consumed, s.released, re)
        end
    end
    gtext(g) = join((S.step_text(m, m.edges[k]) for k in findall(==(g), m.group)), " | ")
    desc = ""
    if kind === :to_re || kind === :to_ss
        gs = [g for g in 1:G if flag(g) == (kind === :to_ss)]
        isempty(gs) && return nothing
        g = rand(rng, gs)
        desc = "group $g -> $(kind === :to_re ? "RE" : "SS"): " * gtext(g)
        setre!(g, kind === :to_re)
    elseif kind === :flip2
        G < 2 && return nothing
        g1, g2 = randperm(rng, G)[1:2]
        desc = "groups $g1, $g2 toggled: " * gtext(g1) * " || " * gtext(g2)
        f1, f2 = flag(g1), flag(g2)
        setre!(g1, !f1)
        setre!(g2, !f2)
    elseif kind === :all_re
        all(flag, 1:G) && return nothing
        foreach(g -> setre!(g, true), 1:G)
        desc = "all RE"
    elseif kind === :del_cat
        core = findall(S.catalytic_core(m))
        isempty(core) && return nothing
        k = rand(rng, core)
        desc = "deleted edge $k: " * S.step_text(m, m.edges[k])
        g = group[k]
        deleteat!(steps, k)
        deleteat!(group, k)
        if !(g in group)
            group = [x > g ? x - 1 : x for x in group]
            isempty(tags) || deleteat!(tags, g)
        end
        isempty(steps) && return nothing
    end
    sp = spec(subs, prods, steps, group, tags, kw.multiplicity, kw.sites)
    mech, err = try
        S.Mechanism(subs, prods, steps, group; tags, multiplicity = kw.multiplicity,
                    sites = kw.sites), nothing
    catch e
        nothing, errmsg(e)
    end
    (desc = desc, spec = sp, mech = mech, ctor_err = err)
end

"""Mutants of `nper` random parents from each of the plain and the allosteric enumerated
sets, every kind per parent, each distinct mechanism once (later copies add provenance)."""
function mutant_records(nper, seed)
    rng = Xoshiro(seed)
    recs = NamedTuple[]
    seen = Dict{Any, Int}()
    for sets in (PLAIN_ENUM, ALLO_ENUM)
        pool = [(set, c) for set in sets
                for c in deserialize(joinpath(SP, "slice/core/data/$set.jls"))]
        picks = randperm(rng, length(pool))[1:min(nper, length(pool))]
        for i in picks
            set, c = pool[i]
            for kind in KINDS
                mu = mutate(c.mech, kind, rng)
                mu === nothing && continue
                prov = "$set/$(c.index) $kind: $(mu.desc)"
                id = mu.mech === nothing ? mu.spec : mu.mech
                j = get(seen, id, 0)
                if j > 0
                    r = recs[j]
                    recs[j] = merge(r, (prov = r.prov * " ;; " * prov,))
                    continue
                end
                push!(recs, (key = "mut/$set/$(c.index)/$kind", set = "mut_$kind",
                             prov = prov, mech = mu.mech, spec = mu.spec,
                             ctor_err = mu.ctor_err))
                seen[id] = length(recs)
            end
        end
    end
    recs
end
