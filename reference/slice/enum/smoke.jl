# ABOUTME: Smoke test of enum.jl: runs the three moves on a few init_bibi_pp parents and
# ABOUTME: compares their children with the package's.
include(joinpath(@__DIR__, "load.jl"))
d = deserialize(joinpath(S, "slice/enum/data/pkg_init_bibi_pp.jls"))
canon(p) = Slice.Converted(p, 0).mech
for r in d.records[1:20]
    m = canon(r.parent)
    for (mv, f) in ((:flip, Slice.flip_children), (:split, Slice.split_children),
                    (:dead, x -> Slice.dead_end_children(x, [:I, :A, :Q])))
        mine = f(m)
        pkg = canon.(getfield(r, mv).children)
        println(r.index, " ", mv, ": mine ", length(mine), " pkg ", length(pkg), " same set ",
                Set(mine) == Set(pkg), " err '", getfield(r, mv).error, "'")
    end
end
