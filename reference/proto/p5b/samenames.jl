# ABOUTME: Counts same-model mechanisms whose fitted name lists are identical to the package's (P5b).
# ABOUTME: Stdlib only.
using Serialization
for p in ("specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi")
    rs = deserialize(joinpath(@__DIR__, "data", "results_$(p).jls"))
    eq = filter(r -> r[:count] == r[:pkg_count], rs)
    println(p, ": same model ", length(eq), ", identical fitted lists ", count(r -> collect(r[:fitted]) == collect(r[:pkg_fitted]), eq),
            ", identical as sets ", count(r -> Set(r[:fitted]) == Set(r[:pkg_fitted]), eq))
end
