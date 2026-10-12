# ABOUTME: Skeptic P5a: counts package samples per allosteric :OnlyA mechanism in each population,
# ABOUTME: to check how many mechanisms V4 actually compared (an empty sample list compares nothing).
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
for name in ("specs", "handwritten", "allo_uniuni", "expand_bibi", "allo_bibi")
    ms = filter(m -> m.allosteric && :OnlyA in m.group_tags, Pop.load_population(name))
    ns = [length(m.samples) for m in ms]
    nnan = sum(count(s -> isnan(s.v), m.samples) for m in ms)
    println(rpad(name, 12), " :OnlyA mechs ", length(ms), "  samples per mech: ", sort(unique(ns)),
            "  with 0 samples: ", count(==(0), ns), "  NaN v samples: ", nnan,
            "  ok=false: ", count(m -> !m.ok, ms))
end
