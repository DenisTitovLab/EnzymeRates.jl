# ABOUTME: Compares the package's allosteric parent counts (less L and the regulator constants)
# ABOUTME: with the slice's solve count, and tallies the offsets.
include(joinpath(@__DIR__, "load.jl"))
for set in ("expand_bibi", "allo_bibi")
    rows = filter(r -> r.allosteric, deserialize(joinpath(S, "slice/enum/data/cmp_$set.jls")))
    off = Dict{Int, Int}()
    for r in rows
        nreg = sum(length(ns) for (_, ns) in Slice.regulator_names(r.mech); init = 0)
        d = r.pkg_count - 1 - nreg - r.count
        off[d] = get(off, d, 0) + 1
    end
    println(set, ": package count minus (L, regulator constants) minus slice count => parents ",
            sort!(collect(off)))
end
