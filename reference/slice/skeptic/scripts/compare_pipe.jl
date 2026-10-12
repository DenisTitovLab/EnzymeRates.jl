# ABOUTME: Compares the skeptic's clean rerun of the pipeline with the integration's saved rows,
# ABOUTME: record by record, and scans every K1 class for members with different metabolites.
include(joinpath(@__DIR__, "indep.jl"))

const MINE = joinpath(SCR, "slice", "skeptic", "rerun", "data")

function main()
    labels = [replace(f, "pipe_" => "", ".jls" => "") for f in readdir(MINE)
              if startswith(f, "pipe_")]
    fields = (:ok, :msg, :key, :order, :nfit, :rref, :ncols, :trees, :complexity, :nQ, :nD,
              :rules, :ops, :printed, :coefficients, :kcat, :halves, :big)
    for l in sort(labels)
        a = Dict((x.set, x.index) => x for x in deserialize(joinpath(IDATA_IN, "pipe_$l.jls")))
        b = Dict((x.set, x.index) => x for x in deserialize(joinpath(MINE, "pipe_$l.jls")))
        diff = Dict{Symbol, Int}()
        for (k, x) in a
            y = get(b, k, nothing)
            y === nothing && (diff[:missing] = get(diff, :missing, 0) + 1; continue)
            for f in fields
                hasproperty(x, f) || continue
                getproperty(x, f) == getproperty(y, f) || (diff[f] = get(diff, f, 0) + 1)
            end
        end
        @printf("%-14s records %6d (mine %6d): %s\n", l, length(a), length(b),
                isempty(diff) ? "identical in every compared field" :
                join(("$k differs in $v" for (k, v) in diff), ", "))
    end
    # K1 classes whose members do not read the same metabolites (merges across names)
    rows = [x for f in readdir(IDATA_IN) if startswith(f, "pipe_")
            for x in deserialize(joinpath(IDATA_IN, f)) if x.ok]
    byset = Dict{String, Vector{Any}}()
    for x in rows
        push!(get!(byset, x.set, Any[]), x)
    end
    println("== K1 classes (within each set) whose members read different metabolite names")
    for (set, xs) in sort!(collect(byset); by = first)
        startswith(set, "init_terter") && continue
        cls = Dict{String, Vector{Int}}()
        foreach(x -> push!(get!(cls, x.key, Int[]), x.index), xs)
        multi = [v for v in values(cls) if length(v) > 1]
        isempty(multi) && continue
        recs = Dict(r.c.index => r for r in skeptic_records(set))
        bad = 0
        example = ""
        for v in multi
            names = unique([S.metabolites(recs[i].c.mech) for i in v])
            length(names) > 1 || continue
            bad += 1
            isempty(example) && (example = join(string.(names), " vs "))
        end
        @printf("  %-18s classes with 2+ members %4d; reading different names %d %s\n", set,
                length(multi), bad, example)
    end
end

main()
