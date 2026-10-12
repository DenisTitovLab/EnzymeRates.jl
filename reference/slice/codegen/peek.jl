# ABOUTME: Prints the program structure (locals, term counts, ops) and optionally the body of
# ABOUTME: the laws named in ARGS as set/id pairs, e.g. expand_bibi expand_bibi/626.
include(joinpath(@__DIR__, "lib.jl"))
const SHOW = get(ENV, "SHOW_BODY", "0") == "1"
for (set, id) in Iterators.partition(ARGS, 2)
    r = only(filter(x -> x.p.id == id, load_records((set,))))
    m = r.c.mech
    P = S.program(m)
    b = S.body(P)
    println("== $id: allosteric $(S.is_allosteric(m)), n $(m.multiplicity), ",
            "locals $(length(P.defs)), ops $(countops(b))")
    for (v, q) in P.defs
        println("  local $(P.vars.names[v]): $(S.nterms(q)) terms, ",
                "$(sum(length, keys(q.terms))) factors")
    end
    println("  num $(S.nterms(P.num)) terms, den $(S.nterms(P.den)) terms")
    SHOW && println(b)
end
