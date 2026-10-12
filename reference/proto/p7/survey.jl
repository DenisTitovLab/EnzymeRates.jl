# ABOUTME: Survey of the P7 candidate mechanisms: denominator sizes, costs and families.
# ABOUTME: Prints the specs and the 30 largest-denominator mechanisms of expand_bibi and init_bibi.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
using Printf
for pop in ("specs", "init_bibi", "expand_bibi")
    ms = Pop.load_population(pop)
    nd(m) = m.allosteric ? sum(length(p[2]) for p in (m.state_polys_A, m.state_polys_I) if !isempty(p); init=0) : length(m.pkg_den)
    nn(m) = m.allosteric ? sum(length(p[1]) for p in (m.state_polys_A, m.state_polys_I) if !isempty(p); init=0) : length(m.pkg_num)
    println("== $pop: $(length(ms)) mechanisms, allosteric $(count(m -> m.allosteric, ms))")
    sel = pop == "specs" ? ms : sort(ms; by = m -> -nd(m))[1:40]
    for m in sel
        @printf("%-16s allo=%d den=%4d num=%4d cx=%5d eval=%6.1fns alloc=%d first=%6.3fs nfit=%2d ndep=%2d fam=%s %s\n",
            m.id, m.allosteric, nd(m), nn(m), m.eq_complexity, m.eval_ns, m.eval_allocs, m.first_call_s,
            length(m.fitted_params), length(m.dependents), join(string.(m.family), ","), first(m.source, 50))
    end
end
