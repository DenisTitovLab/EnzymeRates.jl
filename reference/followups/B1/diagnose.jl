# ABOUTME: Follow-up B1: prints chosen mechanisms in DSL style with tags, door types, components, the
# ABOUTME: structural classes under D/iu and P5a's verdicts for every convention. Usage: julia diagnose.jl id...
include(joinpath(@__DIR__, "b1_lib.jl"))
ids = ARGS
bypop = Dict{String, Vector{String}}()
for id in ids; push!(get!(bypop, split(id, "/")[1], String[]), id); end
for (name, want) in bypop
    pop = Dict(m.id => m for m in load_allo(name))
    recs = load_records(name)
    for id in want
        m = pop[id]
        println("\n", "="^100)
        for rd in (:D, :iu)
            a = analyze(m, rd)
            println("  reading $rd: class=", a.cls, " islands=", a.nisl, " signs=", a.isl_sign, " affinity doors=",
                    a.nadoor, " kinetic doors=", a.nkdoor, " (self-loops ", a.nkself, ", forced ", a.forced,
                    ", cluster-graph edges ", a.nkedge, ", bridges ", a.bridges, ")")
        end
        print_mech(stdout, m, :D)
        for conv in (:i, :ii, :iu)
            r = get(recs, (id, conv), nothing)
            r === nothing && continue
            @printf("     P5a (%s): %-8s cone_dim=%s n=%s spread=%.1e vs_del=%.1e seq=%s\n", conv, verdict(r),
                    get(r, :cone_dim, "-"), get(r, :cone_n, "-"), get(r, :cone_spread, NaN),
                    get(r, :cone_vs_del, NaN), get(r, :cone_seq, "-"))
        end
    end
end
