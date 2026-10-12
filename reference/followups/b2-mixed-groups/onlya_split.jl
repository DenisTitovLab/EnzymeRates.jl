# ABOUTME: Splits the :OnlyA mixed groups of the allosteric populations, and the P5a (i) failures,
# ABOUTME: per mechanism by RE/SS flag and by whether a plain member is a dead-end binding.
include(joinpath(@__DIR__, "..", "..", "proto", "common", "Population.jl"))
include(joinpath(@__DIR__, "mixed_lib.jl"))
using Serialization

const P5A = joinpath(@__DIR__, "..", "..", "proto", "p5a", "data")

"""Per-mechanism profile of its :OnlyA mixed groups: (any RE, any SS, any with a dead-end
plain member, all plain members catalytic)."""
function profile(m)
    mg = [g for g in mixed_groups(m) if m.group_tags[g] === :OnlyA]
    isempty(mg) && return nothing
    re = [m.steps[first(m.groups[g])].re for g in mg]
    de = [any(j -> step_class(m, m.steps[j]) === :plain && deadend(m.species[m.steps[j].to]),
              m.groups[g]) for g in mg]
    (re = any(re), ss = any(!, re), dead = any(de), nmg = length(mg),
     re_groups = count(re), ss_groups = count(!, re))
end

report(label, ps) = println(label, ": mechanisms ", length(ps),
    "; with an RE :OnlyA mixed group ", count(p -> p.re, ps),
    ", with an SS one ", count(p -> p.ss, ps), ", both ", count(p -> p.re && p.ss, ps),
    "; groups RE ", sum(p -> p.re_groups, ps; init = 0), ", SS ", sum(p -> p.ss_groups, ps; init = 0),
    "; with a dead-end plain member ", count(p -> p.dead, ps),
    ", catalytic plain members only ", count(p -> !p.dead, ps))

for name in ("expand_bibi", "allo_bibi")
    pop = Pop.load_population(name)
    ps = filter(!isnothing, [profile(m) for m in pop if m.allosteric])
    report("$name, all :OnlyA mixed", ps)
    byid = Dict(m.id => m for m in pop)
    recs = filter(r -> r[:conv] === :i && !haskey(r, :error), deserialize(joinpath(P5A, "pop_$name.jls")))
    good = [r[:id] for r in recs
            if haskey(r, :cone_v) && r[:cone_spread] < 1e-20 && r[:cone_vs_del] < 1e-20]
    bad = setdiff([r[:id] for r in recs], good)
    report("$name, P5a (i) failures", filter(!isnothing, [profile(byid[id]) for id in bad]))
    report("$name, P5a (i) unique = deletion", filter(!isnothing, [profile(byid[id]) for id in good]))
end
