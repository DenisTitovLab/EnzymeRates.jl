# ABOUTME: Explains Haldane classes (my solve) that are not a chemistry reverse constant: is every
# ABOUTME: chemistry reverse class tied to a free-enzyme binding member (and so kept independent)?
include(joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common", "Thermo.jl"))
include(joinpath(@__DIR__, "census_helpers.jl"))
const T = Thermo
cnt = Dict{String, Int}(); ex = Dict{String, String}()
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    sol = T.solve(m)
    hal = [my_slot(sol, d.col) for d in sol.deps if d.label === :haldane]
    (isempty(hal) || all(cl -> textbook_class(sol, cl), hal)) && continue
    gos = Pop.group_of_step(m)
    # every chemistry step's against constant (or RE K): its class (group, stored slot)
    chemcls = Set{Tuple{Int, Int}}()
    for (i, c) in enumerate(sol.cols)
        s = sol.cs[c.j]
        s.kind >= 2 && (s.re || c.prio[4] == 1) && push!(chemcls, my_slot(sol, i))
    end
    indep_cls = Set(my_slot(sol, i) for i in sol.indep)
    # a class is pinned to the free-enzyme binding priority when it holds a plain binding of a free form
    freebind(cl) = any(j -> sol.cs[j].kind == 1 && sol.cs[j].occ == 0, m.groups[cl[1]])
    key = isempty(chemcls) ? "no chemistry constant" :
          all(cl -> cl in indep_cls && freebind(cl), chemcls) ? "every chemistry reverse class holds a free-enzyme binding (kept)" :
          all(cl -> freebind(cl) || !(cl in indep_cls), chemcls) ? "chemistry reverse classes are free-binding or already dependent" :
          all(cl -> any(j -> sol.cs[j].kind == 1, m.groups[cl[1]]) || !(cl in indep_cls), chemcls) ?
              "every kept chemistry reverse class holds a plain binding member (class takes the binding's priority)" :
          "other"
    cnt[key] = get(cnt, key, 0) + 1
    haskey(ex, key) || (ex[key] = m.id)
end
for (k, v) in sort(collect(cnt)); println(v, "  ", k, "   e.g. ", ex[k]); end
