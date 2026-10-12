# ABOUTME: Smoke test of Allo.derive on the allosteric specs: fitted lists, folds and printed laws.
# ABOUTME: Probe P5b, stdlib only.
include(joinpath(@__DIR__, "Allo.jl"))
pop = filter(m -> m.allosteric, Pop.load_population(ARGS[1]))
sel = length(ARGS) > 1 ? parse.(Int, ARGS[2:end]) : eachindex(pop)
for m in pop[sel]
    md = try Allo.derive(m) catch e; println("== ", m.id, " ERROR ", sprint(showerror, e)); continue end
    println("== ", m.id, " ", m.source)
    println("  keepI matches package: ", md.keepI == m.inactive_kept, "; rows ", md.solve_stats.rows, ", cols ", length(md.cols))
    println("  mine (", length(md.names.fitted), "): ", md.names.fitted)
    println("  pkg  (", length(m.fitted_params), "): ", m.fitted_params)
    println("  folds: ", md.names.folds, " kinds ", md.names.tie_kind)
    println("  deps: ", md.names.deps)
    println("  norm case: ", Allo.norm_case(md))
    println(Allo.law_text(md))
end
