# ABOUTME: Solve-stage check (d): dependents with rational exponents on the skeptics' cases
# ABOUTME: (half_keq, half_exponent, A10g), against the package's, and a census of all records.
include(joinpath(@__DIR__, "lib.jl"))

records = [r for r in load_records() if (try S.solve(r.c.mech); true catch; false end)]
fmt(e) = isinteger(e) ? string(numerator(e)) : string(numerator(e), "/", denominator(e))
function show_case(r)
    m = r.c.mech
    sol = S.solve(m)
    names = S.column_names(m, sol)
    println("== ", r.set, " ", r.p.id, "  ", r.p.source)
    fwd = S.orient(m)
    foreach(k -> println("   ", S.step_text(m, S.oriented(m, k, fwd)), "   [g", m.group[k], "]"),
            eachindex(m.edges))
    for d in sol.deps
        d.label === :kinetic_group && continue
        rhs = join(["$(names[i])^$(fmt(e))" for (i, e) in d.terms], " * ")
        println("   mine: ", names[d.col], " = Keq^", fmt(d.keq), " * ", rhs, "  [", d.label, "]")
    end
    for (d, mono) in r.p.dependent_monomials
        println("   pkg:  ", d, " = ", join(["$(x)^$(fmt(e))" for (x, e) in mono], " * "))
    end
    println("   fitted mine ", length(sol.indep), ", package ", length(r.p.fitted_params))
end
for r in records
    any(occursin(r.p.source), ("half_keq", "half_exponent", "A10g")) && show_case(r)
end

# Census: every record whose dependents carry a non-integer exponent.
nonint, keqfrac, maxden = String[], 0, 1
for r in records
    sol = S.solve(r.c.mech)
    es = [e for d in sol.deps for e in [last.(d.terms); d.keq]]
    any(!isinteger, es) || continue
    push!(nonint, "$(r.set) $(r.p.id) $(r.p.source)")
    any(d -> !isinteger(d.keq), sol.deps) && (global keqfrac += 1)
    global maxden = max(maxden, maximum(denominator, es))
end
@printf("census: %d of %d records have a non-integer exponent (%d with Keq^(p/q)); \
        largest denominator %d\n", length(nonint), length(records), keqfrac, maxden)
foreach(x -> println("    ", x), nonint)
