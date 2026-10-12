# ABOUTME: Skeptic P5a: shows why the package's rate is NaN for AD2c: its inactive-state polynomials
# ABOUTME: (num, den, free weight) over a pruned graph with two disconnected components.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
using Serialization
m = only(filter(x -> x.id == "AD2c_two_free_disconnected", deserialize(joinpath(@__DIR__, "data", "adv_cases.jls"))))
println("inactive_kept: ", m.inactive_kept)
for (lab, p) in zip(("num_I", "den_I", "dfree_I"), m.state_polys_I)
    println(lab, " = ", isempty(p) ? "0 (empty polynomial)" : p)
end
for (lab, p) in zip(("num_A", "den_A", "dfree_A"), m.state_polys_A)
    println(lab, " terms: ", length(p))
end
println("rate_equation_string: ", first(m.rate_equation_string, 400))
println("v samples: ", [s.v for s in m.samples])
