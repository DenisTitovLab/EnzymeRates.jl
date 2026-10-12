# ABOUTME: Prints the package-side fields of the slowest p1_enum_inh records (timings, sizes).
include(joinpath(@__DIR__, "..", "..", "proto", "common", "Population.jl"))
using Serialization
pops = deserialize(joinpath(@__DIR__, "..", "..", "proto", "skeptic_P1-engine", "data", "enum_inh.jls"))
for id in ("bi_bi_with_inh/1013", "bi_bi_with_inh/959", "bi_bi_with_inh/584", "bi_bi_with_inh/977")
    p = only(filter(x -> x.id == id, pops))
    println(id, ": ok ", p.ok, " eval_ns ", p.eval_ns, " first_call_s ", p.first_call_s,
            " eq_complexity ", p.eq_complexity, " forms ", length(p.species), " steps ",
            length(p.steps), " re ", count(s -> s.re, p.steps), " fitted ", length(p.fitted_params),
            " inhibitors ", p.inhibitors, " pkg_den ", length(p.pkg_den))
end
