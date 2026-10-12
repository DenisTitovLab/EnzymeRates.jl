# ABOUTME: The slowest laws of a saved ter-ter pipeline run: stage times, sizes and complexity,
# ABOUTME: then each stage of the slowest re-run with its own timing; usage: peek_slow.jl label set
include(joinpath(@__DIR__, "lib.jl"))
rows = filter(x -> x.ok, deserialize(joinpath(IDATA, "pipe_$(ARGS[1]).jls")))
st(x) = x.t_solve + x.t_law + x.t_key + x.t_print + x.t_code + x.t_kcat
for x in first(sort(rows; by = x -> -st(x)), 8)
    @printf("%-20s sum %7.0f ms: law %.0f key %.0f print %.0f code %.0f kcat %.0f | nD %d nQ %d ops %d trees %s complexity %d forms %d edges %d\n",
            x.id, 1e3st(x), 1e3x.t_law, 1e3x.t_key, 1e3x.t_print, 1e3x.t_code, 1e3x.t_kcat,
            x.nD, x.nQ, x.ops, x.trees, x.complexity, x.nforms, x.nedges)
end
@printf("laws with V×τ above 5000 (the package's guard would refuse them): %d of %d; largest D %d\n",
        count(x -> x.complexity > 5000, rows), length(rows), maximum(x -> x.nD, rows))
x = first(sort(rows; by = x -> -st(x)))
r = only(filter(r -> r.c.index == x.index, records_of(ARGS[2])))
m = r.c.mech
sol = S.solve(m)
t = @elapsed L = S.rate_law(m, sol)
@printf("re-run %s: law %.2f s\n", x.id, t)
t = @elapsed c = S.canonical(L); @printf("  key %.2f s\n", t)
t = @elapsed S.rate_equation_string(m, L, sol); @printf("  printed %.2f s\n", t)
t = @elapsed S.coefficient_form(m, L); @printf("  coefficients %.2f s\n", t)
t = @elapsed S.normalizer(m, L); @printf("  normalizer %.2f s\n", t)
t = @elapsed S.body(S.program(m, sol, L)); @printf("  code %.2f s\n", t)
t = @elapsed try S.kcat_candidates(m, sol) catch end; @printf("  kcat %.2f s\n", t)
