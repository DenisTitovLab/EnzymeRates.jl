# ABOUTME: Allosteric checks (d), (e), (f): no phantom names or fold hazards in the MWC laws, fitted
# ABOUTME: lists against the package's by kinetic class, and admission and derivation times.
include(joinpath(@__DIR__, "lib.jl"))

const POP1993 = ("specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi")
const POP848 = ("allo_uniuni_inh", "allo_unibi", "allo_bibi_pp", "allo_terbi")

"""Each package fitted name of record `r` as the set of kinetic classes (`Slice.class`) of
the steps that carry it, active and inactive (the package's own pruned inactive names)."""
function pkg_classes(r)
    m, c, p = r.c.mech, r.c, r.p
    fwd = S.orient(m)
    out = Dict{Symbol, Set{Any}}()
    for (j, s) in enumerate(p.steps), st in (:A, :I)
        names = st === :A ? p.step_names_A[j] : p.step_names_I_pkg[j]
        rev = st === :I && !isempty(p.inactive_reversed) && p.inactive_reversed[j]
        k = c.edge_of[j]
        for (q, n) in enumerate(names)
            slot = first(pop_slot(c, fwd, j, q; reversed = rev))
            col = S.Col(k, slot, st === :I && m.tags[m.group[k]] !== :EqualAI)
            push!(get!(Set{Any}, out, n), S.class(m, col))
        end
    end
    out
end

records = [r for r in load_records() if S.is_allosteric(r.c.mech)]
phantom, unused, hazard, nlaw = 0, 0, 0, 0
lists = Dict{Any, Int}()
ta, tl = Float64[], Float64[]
foreach(r -> (S.admit(r.c.mech); S.allo_law(r.c.mech)), records[1:20])     # compile
for r in records
    m, p = r.c.mech, r.p
    push!(ta, minimum(@elapsed(S.admit(m)) for _ in 1:3))
    push!(tl, minimum(@elapsed(S.allo_law(m, S.solve(m))) for _ in 1:3))
    A = S.allo_law(m, S.solve(m))
    sol = S.solve(m)
    fitted = S.fitted_params(m)
    # (d) every constant the law reads is a parameter, and every fitted one is read
    used = Set{Symbol}()
    polys = [A.H, A.active.Q, A.active.D, A.fA, A.fI, (s[i] for s in A.sites for i in 1:2)...]
    A.inactive === nothing || push!(polys, A.inactive.Q, A.inactive.D)
    foreach(q -> foreach(mo -> foreach(((v, _),) -> push!(used, A.vars.names[v]), mo),
                         keys(q.terms)), polys)
    consts = setdiff(used, S.metabolites(m))
    global phantom += length(setdiff(consts, Set([fitted...; :Keq; :L])))
    nu = setdiff(Set(fitted), used, (:L,))
    isempty(nu) || println("    fitted but unread: ", r.set, " ", p.id, " ", collect(nu))
    global unused += length(nu)
    names = S.column_names(m, sol)
    global hazard += !allunique(names[sol.indep])
    global nlaw += 1
    # (e) fitted lists against the package's, by kinetic class
    r.set in (POP1993..., POP848...) || continue
    d = length(fitted) - length(p.fitted_params)
    regs = Set(x for (_, ns) in S.regulator_names(m) for x in ns)
    mine = Set(S.class(m, sol.cols[i]) for i in sol.indep)
    pc = pkg_classes(r)
    pkg = Set(); missing_ = 0
    for n in p.fitted_params
        (n === :L || endswith(String(n), "reg")) && continue
        haskey(pc, n) ? union!(pkg, pc[n]) : (missing_ += 1)
    end
    key = (r.set in POP1993 ? "1993" : "848", d, d == 0 ? (mine == pkg ? "same classes" :
           "different classes") : "-", missing_ > 0 ? "package names without steps" : "")
    lists[key] = get(lists, key, 0) + 1
end
println("(d) laws: $nlaw; constants read by a law that are not parameters: $phantom; ",
        "fitted constants no law polynomial reads (L excluded): $unused; laws whose ",
        "independent columns share a name: $hazard")
println("(e) fitted lists (population, mine - package length, classes, note) => mechanisms:")
foreach(((k, v),) -> println("    ", k, " => ", v), sort!(collect(lists); by = string))
q(t, x) = 1e6 * quantile(t, x)
big3 = sortperm(tl; rev = true)[1:3]
println("slowest derivations: ", [(records[i].p.id, round(1e3 * tl[i]; digits = 2)) for i in big3], " ms")
@printf("(f) admit: median %.1f µs, p90 %.1f, max %.1f; solve + allo_law: median %.1f µs, p90 %.1f, max %.0f (%d mechanisms)\n",
        q(ta, 0.5), q(ta, 0.9), 1e6 * maximum(ta), q(tl, 0.5), q(tl, 0.9), 1e6 * maximum(tl),
        length(tl))
