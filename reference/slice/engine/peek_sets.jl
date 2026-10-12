# ABOUTME: Lists each converted set's plain and allosteric counts and the ids of the
# ABOUTME: hand-built and adversarial records, to plan the engine-stage checks.
include(joinpath(@__DIR__, "..", "solve", "lib.jl"))
recs = load_records()
for set in unique(r.set for r in recs)
    rs = filter(r -> r.set == set, recs)
    np = count(r -> !S.is_allosteric(r.c.mech), rs)
    println(rpad(set, 20), " records ", length(rs), " plain ", np)
end
for r in recs
    r.set in ("p1_adv_pkg", "p1_adv_pkg2", "p1_direct", "p2_adv_pkg", "p2_adv_pkg2", "p2_hand",
              "p2_fuzz_dump", "p2_fuzz_dump_flip") || continue
    println(r.set, " ", r.p.id, " | ", r.p.source, " | ok=", r.p.ok, " allo=", S.is_allosteric(r.c.mech))
end
