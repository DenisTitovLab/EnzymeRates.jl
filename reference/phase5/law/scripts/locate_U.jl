# ABOUTME: Finds which slice function holds variant U's hoisted boxes: calls each derivation stage once,
# ABOUTME: printing a marker to stderr before each first call, so its compile remarks fall in its segment.
include(joinpath(@__DIR__, "lib.jl"))
const S = load_variant(get(ENV, "VARIANT", "U"))
rows = census_rows()
pick(set) = first(r for r in rows if r.set == set && S.is_allosteric(slice_mech(S, r.net)))
ma = slice_mech(S, first(r for r in rows if group(r.set) == "pingpong" &&
                         S.is_allosteric(slice_mech(S, r.net))).net)
mp = slice_mech(S, first(r for r in rows if group(r.set) == "bibi").net)
mark(s) = (println(stderr, "=== VARIANT ", s); flush(stderr); println("=== VARIANT ", s); flush(stdout))
for (label, m) in (("plain", mp), ("allo", ma))
    mark("$label:solve");   sol = S.solve(m)
    mark("$label:law");     Lw = S.law(m, sol)
    mark("$label:rate_law"); L = S.rate_law(m, sol)
    mark("$label:print");   S.rate_equation_string(m, L, sol)
    mark("$label:coef");    S.coefficient_form(m, L)
    mark("$label:key");     S.canonical(L)
    mark("$label:program"); p = S.program(m, sol, L)
    mark("$label:body");    S.body(p)
    mark("$label:class5");  class5(m, L)
    mark("$label:compile"); em = S.compile(m)
    mark("$label:generated"); S.rate_equation(em, NamedTuple{Tuple(S.metabolites(m))}(Tuple(1.0 for _ in S.metabolites(m))),
                                              NamedTuple{Tuple(S.parameters(m))}(Tuple(1.0 for _ in S.parameters(m))))
end
mark("done")
