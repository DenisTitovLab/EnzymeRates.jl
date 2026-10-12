# ABOUTME: Test plumbing for probe P6: Pop.PMech records as raw GM input (with equivalent and
# ABOUTME: mutated variants), the package's structural key, and a batch timer; stdlib only.
isdefined(Main, :Pop) || include(joinpath(@__DIR__, "..", "common", "Population.jl"))
isdefined(Main, :GM) || include(joinpath(@__DIR__, "GraphMech.jl"))
using Random

_testhelper_met(x::Pop.PMet) =
    x.role in (:S, :P, :Inh) ? GM.Met(x.name, x.role === :Inh) : error("role $(x.role) in a step")

"""Raw GM input of `m`: `(substrates, products, steps, group, kw)`. With `rng`, the steps,
group labels, sites, ligands and every list come in a random order; `flip[j]` writes step
`j` in reverse. Every step gets its own form objects, as a DSL would write them."""
function _testhelper_toraw(m::Pop.PMech; rng = nothing, flip = falses(length(m.steps)))
    sh(v) = rng === nothing ? copy(v) : shuffle(rng, v)
    perm(n) = rng === nothing ? collect(1:n) : randperm(rng, n)
    form(i) = (sp = m.species[i]; GM.RawForm(sp.conformation, sh(_testhelper_met.(sp.bound)),
                                             sh(sp.res_added), sh(sp.res_subtracted)))
    gof = Pop.group_of_step(m)
    label = perm(length(m.groups))
    steps, group = GM.RawStep[], Int[]
    for j in perm(length(m.steps))
        s = m.steps[j]
        f, t, c, r = form(s.from), form(s.to), sh(_testhelper_met.(s.consumed)), sh(_testhelper_met.(s.released))
        push!(steps, flip[j] ? GM.RawStep(t, f, r, c, s.re) : GM.RawStep(f, t, c, r, s.re))
        push!(group, label[gof[j]])
    end
    tags = Symbol[]
    if m.allosteric
        tags = Vector{Symbol}(undef, length(m.groups))
        for g in eachindex(m.groups); tags[label[g]] = m.group_tags[g]; end
    end
    sites = GM.Site[]
    for s in m.sites
        p = perm(length(s.ligands))
        push!(sites, GM.Site(s.ligands[p], s.states[p], s.multiplicity))
    end
    rng === nothing || shuffle!(rng, sites)
    (sh(m.substrates), sh(m.products), steps, group,
     (; tags, multiplicity = m.cat_multiplicity, sites))
end

_testhelper_build(raw) = GM.canonicalize(raw[1], raw[2], raw[3], raw[4]; raw[5]...)
_testhelper_gm(raw) = first(_testhelper_build(raw))

"""The package's canonical structure of `m` as text (forms by rendered name, steps in
stored order with their group, tags, sites): equal exactly when the package's mechanisms are
equal, the reaction aside."""
_testhelper_pkgkey(m::Pop.PMech) = join([join(m.substrates, ","), join(m.products, ","),
    string(m.allosteric, m.cat_multiplicity), join(m.group_tags, ","),
    [string(Pop.step_text(m, s), " @", g) for (s, g) in zip(m.steps, Pop.group_of_step(m))]...,
    [string(s.ligands, s.states, s.multiplicity) for s in m.sites]...], " | ")

"""Median over `reps` batches of the per-call time in ns of `f()`, run `k` times a batch."""
function _testhelper_time_ns(f, k::Int; reps = 11)
    Base.donotdelete(f())
    ts = Vector{Float64}(undef, reps)
    for r in 1:reps
        t0 = time_ns()
        for _ in 1:k; Base.donotdelete(f()); end
        ts[r] = (time_ns() - t0) / k
    end
    sort!(ts)[(reps + 1) ÷ 2]
end

_testhelper_load_all() =
    [m for pop in ("specs", "handwritten", "init_uniuni", "init_unibi", "init_bibi",
                   "init_bibi_pp", "expand_bibi", "allo_uniuni", "allo_bibi")
     for m in Pop.load_population(pop)]
