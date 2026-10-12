# ABOUTME: What the slice derives for an all-rapid-equilibrium mechanism and for a mechanism
# ABOUTME: whose declared substrate never takes part: the law, printed form and values.
include(joinpath(@__DIR__, "..", "advlib.jl"))
E(b...; kw...) = F(:E, b...; kw...)
p = build("all_re", [:S], [:P], [
    st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(:P); re = true, g = 2),
    st(E(:P), E(); r = [:P], re = true, g = 3)])
m = S.Converted(p, 0).mech
d = S.derive(m)
println("fitted ", S.fitted_params(m))
println("Q terms ", S.nterms(d.law.Q), "; D terms ", S.nterms(d.law.D))
println(d.printed)
v = S.rate_equation(m, (S = 2.0, P = 0.5), (; (x => 1.5 for x in S.parameters(m))...))
println("rate at S = 2, P = 0.5 with every constant 1.5: ", v)
println("kcat: ", try S.kcat(m, (; (x => 1.5 for x in S.parameters(m))...)) catch e; sprint(showerror, e) end)
q = build("missing", [:A, :B], [:P], [
    st(E(), E(:A); c = [:A], re = true, g = 1), st(E(:A), E(); r = [:P], g = 2)])
mq = S.Converted(q, 0).mech
println("missing substrate: constructor accepted; check_boundary: ",
        try S.check_boundary(mq); "passes" catch e; sprint(showerror, e) end)
println("derive: ", try S.derive(mq); "derived" catch e; first(sprint(showerror, e), 200) end)
