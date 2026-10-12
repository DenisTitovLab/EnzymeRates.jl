# ABOUTME: The skeptic's adversarial-spec builders (Mod, F, st, build) copied from
# ABOUTME: slice/skeptic/scripts/advlib.jl without its oracle includes, then its 40 cases.

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

include(joinpath(SP, "slice/skeptic/scripts/adv_cases.jl"))

"""The skeptic's adversarial cases as records (a spec even when the slice rejects one)."""
function adv_records()
    recs = NamedTuple[]
    for (i, cs) in enumerate(CASES)
        p = try
            cs.thunk()
        catch e
            println("adv $(cs.name): build error ", errmsg(e))
            continue
        end
        mech, err = try
            S.Converted(p, i).mech, nothing
        catch e
            nothing, errmsg(e)
        end
        sp = try
            spec(S.pop_raw(p))
        catch e
            println("adv $(cs.name): no spec ", errmsg(e))
            continue
        end
        push!(recs, (key = "skeptic_adv/$(cs.name)", set = "skeptic_adv",
                     prov = "adv:$(cs.name) expect=$(cs.expect)", mech = mech, spec = sp,
                     ctor_err = err))
    end
    recs
end
