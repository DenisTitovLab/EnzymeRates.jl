# ABOUTME: Compares which allosteric records have an inactive state that turns over, by my
# ABOUTME: inactive law and by the package's inactive state polynomials.
include(joinpath(@__DIR__, "lib.jl"))
recs = [r for r in load_records() if S.is_allosteric(r.c.mech) && !isempty(r.p.state_polys_I)]
agree, n, ex = 0, 0, String[]
for r in recs
    m = r.c.mech
    sol = try S.solve(m) catch; continue end
    mask = S.inactive_mask(m)
    mine = any(mask) && !iszero(S.law(m, sol, mask; inactive = true).Q)
    pkg = !isempty(r.p.state_polys_I[1])
    global n += 1
    mine == pkg ? (global agree += 1) : push!(ex, "$(r.set) $(r.p.id): mine $mine package $pkg")
end
println("inactive state turns over: same verdict as the package in $agree of $n")
foreach(println, first(ex, 10))
println("differences by set: ", Dict(s => count(x -> startswith(x, s * " "), ex) for s in unique(first.(split.(ex, " ")))))
