# ABOUTME: Re-reads the P5a convention-(i) population records and classifies the mixed :OnlyA
# ABOUTME: groups of the failing mechanisms (member kinds, RE/SS flag); prints example mechanisms.
include(joinpath(@__DIR__, "..", "..", "proto", "common", "Population.jl"))
include(joinpath(@__DIR__, "mixed_lib.jl"))
using Serialization, Printf

const P5A = joinpath(@__DIR__, "..", "..", "proto", "p5a", "data")

for name in ("expand_bibi", "allo_bibi")
    pop = Dict(m.id => m for m in Pop.load_population(name))
    recs = filter(r -> r[:conv] === :i && !haskey(r, :error),
                  deserialize(joinpath(P5A, "pop_$name.jls")))
    fail = [r for r in recs
            if !(haskey(r, :cone_v) && r[:cone_spread] < 1e-20 && r[:cone_vs_del] < 1e-20)]
    with_mixed = 0; flags = Dict{Bool, Int}(); members = Dict{Symbol, Int}()
    plain_kind = Dict{Symbol, Int}(); tc_elsewhere = 0; pingpong = 0
    shapes = Dict{String, Int}()
    for r in fail
        m = pop[r[:id]]
        mg = [g for g in mixed_groups(m) if m.group_tags[g] === :OnlyA]
        :pingpong in m.family && (pingpong += 1)
        isempty(mg) && continue
        with_mixed += 1
        any(s -> step_class(m, s) === :tc, m.steps) && (tc_elsewhere += 1)
        for g in mg
            re = m.steps[first(m.groups[g])].re
            flags[re] = get(flags, re, 0) + 1
            for j in m.groups[g]
                c = step_class(m, m.steps[j])
                if c === :plain
                    k = deadend(m.species[m.steps[j].to]) ? :deadend : :catalytic
                    plain_kind[k] = get(plain_kind, k, 0) + 1
                else
                    members[c] = get(members, c, 0) + 1
                end
            end
            sig = (re ? "RE " : "SS ") * group_text(m, g)
            shapes[sig] = get(shapes, sig, 0) + 1
        end
    end
    @printf("== %s, convention (i): %d records, %d not (unique and = deletion)\n",
            name, length(recs), length(fail))
    @printf("   failing with an :OnlyA mixed group: %d (ping-pong among all failures: %d)\n",
            with_mixed, pingpong)
    println("   their mixed groups by flag (true = RE): ", sort!(collect(flags)))
    println("   chemistry members by class: ", sort!(collect(members)),
            "; plain members: ", sort!(collect(plain_kind)))
    println("   failing mixed-group mechanisms holding a Theorell-Chance step anywhere: ",
            tc_elsewhere)
    for (sig, c) in sort!(collect(shapes); by = x -> (-x[2], x[1]))
        @printf("   %4d  %s\n", c, sig)
    end
end

println("\n== example mechanisms")
for (name, id) in (("init_bibi", "init_bibi/56"), ("init_bibi", "init_bibi/57"),
                   ("allo_bibi", "allo_bibi/826"), ("expand_bibi", "expand_bibi/439"))
    m = only(filter(x -> x.id == id, Pop.load_population(name)))
    print_mech(m)
end
