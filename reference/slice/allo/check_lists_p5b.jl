# ABOUTME: Attributes the fitted-list moves of check (e): P5b's fitted lists, mine and the package's,
# ABOUTME: compared by kinetic class on the 1,844 equal-length mechanisms of the 1,993.
include(joinpath(@__DIR__, "lib.jl"))

"""Each name of `names_by_step` (Pop step => its constants' names, `rev` the steps whose names
run against the Pop step) as the set of kinetic classes of the steps carrying it."""
function classes_of(r, names_by_step, inactive; rev = falses(length(r.p.steps)))
    m, c = r.c.mech, r.c
    fwd = S.orient(m)
    out = Dict{Symbol, Set{Any}}()
    for (j, names) in enumerate(names_by_step), (q, n) in enumerate(names)
        k = c.edge_of[j]
        slot = first(pop_slot(c, fwd, j, q; reversed = rev[j]))
        col = S.Col(k, slot, inactive && m.tags[m.group[k]] !== :EqualAI)
        push!(get!(Set{Any}, out, n), S.class(m, col))
    end
    out
end

p5b = Dict(x[:id] => x for f in readdir(joinpath(PROTO, "p5b", "data"); join = true)
           for x in deserialize(f) if x isa AbstractDict && haskey(x, :fitted))
records = [r for r in load_records(("specs", "handwritten", "allo_uniuni", "allo_bibi",
                                    "expand_bibi")) if S.is_allosteric(r.c.mech)]
tab = Dict{Any, Int}()
for r in records
    m, p = r.c.mech, r.p
    haskey(p5b, p.id) || continue
    sol = S.solve(m)
    fitted = S.fitted_params(m)
    length(fitted) == length(p.fitted_params) || continue
    rev = isempty(p.inactive_reversed) ? falses(length(p.steps)) : p.inactive_reversed
    named = [classes_of(r, p.step_names_A, false), classes_of(r, p.step_names_I, true),
             classes_of(r, p.step_names_I_pkg, true; rev)]
    cls(n, which) = reduce(union, (get(named[w], n, Set()) for w in which); init = Set())
    keep(n) = n !== :L && !endswith(String(n), "reg")
    pkg = reduce(union, (cls(n, (1, 3)) for n in p.fitted_params if keep(n)); init = Set())
    q = reduce(union, (cls(n, (1, 2)) for n in p5b[p.id][:fitted] if keep(n)); init = Set())
    mine = Set(S.class(m, sol.cols[i]) for i in sol.indep)
    key = (mine == pkg ? "mine = package" : "mine ≠ package",
           q == pkg ? "P5b = package" : "P5b ≠ package", mine == q ? "mine = P5b" : "mine ≠ P5b")
    tab[key] = get(tab, key, 0) + 1
end
foreach(((k, v),) -> println(k, " => ", v), sort!(collect(tab); by = string))
