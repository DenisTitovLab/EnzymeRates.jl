# ABOUTME: Lists, for the same-model pairs the LM estimate found in expand_bibi, whether each
# ABOUTME: mechanism is allosteric and how many fitted parameters it has.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
ms = Pop.load_population("expand_bibi"); byid = Dict(m.id => m for m in ms)
for (a, b) in [(304, 305), (303, 304), (556, 614), (22, 23), (201, 202), (181, 182), (42, 43), (139, 140)]
    ma, mb = byid["expand_bibi/$a"], byid["expand_bibi/$b"]
    println(a, " ~ ", b, ": allosteric ", ma.allosteric, "/", mb.allosteric, ", fitted ",
            length(ma.fitted_params), "/", length(mb.fitted_params), ", RE/SS steps ",
            count(s -> s.re, ma.steps), "+", count(s -> !s.re, ma.steps), " vs ",
            count(s -> s.re, mb.steps), "+", count(s -> !s.re, mb.steps))
end
