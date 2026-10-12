# ABOUTME: Compares the exact block-level speed-free predicates with the numeric speed-free lists that the
# ABOUTME: re-run validation stored for every plain population mechanism (thermo and bare constants).
include(joinpath(@__DIR__, "SkLib.jl"))
using Serialization
rows = Dict(r[:id] => r for r in deserialize(joinpath(@__DIR__, "rerun", "results.jls")))
nt = nb = n = 0
for p in ["specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi", "handwritten"]
    for m in Pop.load_population(p)
        m.allosteric && continue
        r = rows[m.id]; global n += 1
        th, ba = Sk.speed_free_exact(m)
        global nt += th == r[:speed_free_thermo]; global nb += ba == r[:speed_free_bare]
    end
end
println("exact block predicate = numeric thermo in $nt of $n; exact bridge predicate = numeric bare in $nb of $n")
