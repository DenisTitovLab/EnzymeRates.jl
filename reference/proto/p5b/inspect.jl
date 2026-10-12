# ABOUTME: Inspects P5b results: artifact matches in the package strings, rank deficits where counts differ.
# ABOUTME: Stdlib only.
using Serialization
const ARTIFACT = r"(?<![\w.])1 \^ \d|\^ 1(?![\d.])|(?<![\w.^])1 \* |\* 1(?![\d.\w])|L \* 0"
for name in ("specs", "allo_uniuni", "allo_bibi", "expand_bibi")
    rs = deserialize(joinpath(@__DIR__, "data", "results_$(name).jls"))
    kinds = Dict{String, Int}()
    for r in rs
        include_string(Main, "1") # no-op
    end
    pop = Dict(r[:id] => r for r in rs)
    println("== $name")
    # artifact kinds in the package strings (context of each match)
    include(joinpath(@__DIR__, "..", "common", "Population.jl"))
    P = Dict(m.id => m for m in Pop.load_population(name))
    for r in rs
        s = P[r[:id]].rate_equation_string
        for mt in eachmatch(ARTIFACT, s)
            k = mt.match
            kinds[k] = get(kinds, k, 0) + 1
        end
    end
    println("   package artifact matches: ", kinds)
    ex = findfirst(r -> r[:artifacts_pkg] > 0, rs)
    if ex !== nothing
        s = P[rs[ex][:id]].rate_equation_string
        mt = match(ARTIFACT, s)
        println("   e.g. ", rs[ex][:id], ": …", s[max(1, mt.offset - 120):min(end, mt.offset + 40)], "…")
    end
    dif = filter(r -> r[:count] != r[:pkg_count] && haskey(r, :rank_mine), rs)
    tab = Dict{Tuple{Int, Int, Int}, Int}()
    for r in dif
        k = (r[:pkg_count] - r[:count], r[:count] - r[:rank_mine][1], r[:pkg_count] - r[:rank_pkg][1])
        tab[k] = get(tab, k, 0) + 1
    end
    println("   count-differing (pkg-mine, mine deficit, pkg deficit) => n: ", sort!(collect(tab)))
    same = filter(r -> r[:count] == r[:pkg_count] && haskey(r, :rank_mine), rs)
    println("   equal counts: rank mine == rank pkg in ", count(r -> r[:rank_mine][1] == r[:rank_pkg][1], same),
            " of ", length(same), "; deficits ", sort!(collect(Dict(d => count(r -> r[:count] - r[:rank_mine][1] == d, same)
                                                                  for d in unique(r[:count] - r[:rank_mine][1] for r in same)))))
end
