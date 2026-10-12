# ABOUTME: Counts laws whose kcat is identically zero (every candidate's numerator vanishes) in a
# ABOUTME: random sample of every population, and prints the steps of one ter-ter example.
include(joinpath(@__DIR__, "indep.jl"))

function main()
    rng = Xoshiro(17)
    rows = [x for f in readdir(IDATA_IN) if startswith(f, "pipe_")
            for x in deserialize(joinpath(IDATA_IN, f)) if x.ok && !(x.kcat isa String)]
    bypop = Dict{String, Vector{Any}}()
    foreach(x -> push!(get!(bypop, startswith(x.set, "init_terter") ? "init_terter" : x.set,
                            Any[]), x), rows)
    example = nothing
    for (pop, xs) in sort!(collect(bypop); by = first)
        pick = pop == "init_terter" ? xs[randperm(rng, length(xs))[1:3000]] : xs
        nzero = 0
        for set in unique(x.set for x in pick)
            idx = [x.index for x in pick if x.set == set]
            mechs = if startswith(set, "init_terter") || set == "expand_terter"
                ps = deserialize(joinpath(IDATA_IN, "pop_$set.jls"))
                [S.Converted(ps[i], i).mech for i in idx]
            else
                recs = Dict(r.c.index => r for r in skeptic_records(set))
                [recs[i].c.mech for i in idx]
            end
            GC.gc()
            for m in mechs
                _, terms = S.kcat_candidates(m)
                if all(t -> iszero(t[1]), terms)
                    nzero += 1
                    example === nothing && pop == "init_terter" && (example = m)
                end
            end
        end
        nzero > 0 && @printf("%-16s laws checked %6d, kcat identically 0: %d\n", pop,
                             length(pick), nzero)
        flush(stdout)
    end
    if example !== nothing
        fwd = S.orient(example)
        println("example (ter-ter):")
        for k in eachindex(example.edges)
            println("  ", S.step_text(example, S.oriented(example, k, fwd)), "   [g",
                    example.group[k], "]")
        end
        println("  catalytic core: ", Int.(S.catalytic_core(example)))
    end
end

main()
