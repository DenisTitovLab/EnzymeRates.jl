# ABOUTME: Follow-up B1: lists ids of one structural class under a reading with P5a's verdict.
# ABOUTME: Usage: julia list_class.jl pop reading class conv
include(joinpath(@__DIR__, "b1_lib.jl"))
name, rd, cls, conv = ARGS[1], Symbol(ARGS[2]), ARGS[3], Symbol(ARGS[4])
recs = load_records(name)
for m in load_allo(name)
    a = analyze(m, rd)
    a.cls == cls && println(m.id, "  ", verdict(get(recs, (m.id, conv), nothing)), "  kdoors=", a.nkdoor, " adoors=", a.nadoor, " isl=", a.nisl, " signs=", a.isl_sign)
end
