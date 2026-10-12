# ABOUTME: Follow-up B1: counts mechanisms whose blocked-step door typing differs between R_D (per step) and
# ABOUTME: P5a's group-tied (ii), and between R_iu and (i); stdlib only.
include(joinpath(@__DIR__, "b1_lib.jl"))
for name in POPS
    pop = filter(m -> :OnlyA in m.group_tags, load_allo(name))
    d1 = count(m -> (b = blocked_flags(m); affinity_flags(m, :D)[b] != affinity_flags(m, :ii)[b]), pop)
    d2 = count(m -> (b = blocked_flags(m); affinity_flags(m, :iu)[b] != affinity_flags(m, :i)[b]), pop)
    println(rpad(name, 16), " D≠ii typing: ", d1, "   iu≠i typing: ", d2, "   of ", length(pop))
end
