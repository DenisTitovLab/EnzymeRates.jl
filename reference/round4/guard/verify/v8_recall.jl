# ABOUTME: Recall and precision of the probe's structural face check against the independent
# ABOUTME: oracle (package constants) on every plain level-1 enumerator child and a level-2 sample.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard/lib.jl")
include(joinpath(@__DIR__, "vpkglib.jl"))
setprecision(BigFloat, 512)

# to_pkg and face_gaps as in v5_census.jl
include(joinpath(@__DIR__, "v5_parts.jl"))

const PLAIN = ["init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi", "p1_enum_inh"]

function tab(label, ms, rng; ndraw = 2)
    cnt = Dict{Tuple{String, String}, Int}()
    exs = Dict{Tuple{String, String}, Vector{Any}}()
    for m in ms
        S.is_allosteric(m) && continue
        probe = try
            f = something(face_check(m, m.products), face_check(m, m.substrates), "")
            isempty(f) ? "pass" : "flag: " * split(f, " ")[1]
        catch e
            "face_check error"
        end
        oracle = try
            pm = to_pkg(m)
            gaps = Float64[]; nonu = false
            for _ in 1:ndraw
                g, subs, prods = face_gaps(pm, rng)
                for z in (prods, subs)
                    g[z] === :nonunique ? (nonu = true) : push!(gaps, Float64(g[z]))
                end
            end
            nonu ? "nonunique" : maximum(gaps) > 1e-6 ? "gap" : "no gap"
        catch e
            "error: " * first(errline(e), 40)
        end
        k = (probe, oracle)
        cnt[k] = get(cnt, k, 0) + 1
        length(get!(exs, k, Any[])) < 2 && push!(exs[k], m)
    end
    println("\n  $label: (probe face check, oracle) counts")
    for (k, v) in sort!(collect(cnt))
        println("    ", rpad(string(k), 50), v)
    end
    for k in keys(exs)
        (startswith(k[1], "flag") && k[2] == "no gap" || k[1] == "pass" && k[2] in ("gap", "nonunique")) || continue
        for m in exs[k]
            println("    disagreement $k: ", join([S.step_text(m, e) for e in m.edges], "; "))
        end
    end
end

function main()
    rng = Xoshiro(31)
    parents = [c.mech for set in PLAIN for c in deserialize(joinpath(SP, "slice/core/data/$set.jls"))]
    kids = S.Mechanism[]
    for m in parents, f in (S.flip_children, S.split_children)
        try append!(kids, f(m)) catch end
    end
    unique!(kids)
    plainkids = filter(!S.is_allosteric, kids)
    tab("plain level-1 children ($(length(plainkids)))", plainkids, rng)
    sample = plainkids[randperm(rng, length(plainkids))[1:800]]
    kids2 = S.Mechanism[]
    for m in sample, f in (S.flip_children, S.split_children)
        try append!(kids2, f(m)) catch end
    end
    unique!(kids2)
    tab("level-2 children of 800 random plain level-1 children ($(length(kids2)))", kids2, rng)
end

main()
