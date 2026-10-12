# ABOUTME: Exploration for the MWC printer: how often the cross weights are 1, how each state's
# ABOUTME: normalizer relates to its flip weight, and whether the normalized cross weights agree.
include(joinpath(@__DIR__, "lib.jl"))

"""The flip weight of law `L` divided by its own normalizer, as `printed` groups."""
flipgroups(m, L, f) = S.printed(m, S.Law(L.vars, L.H, f, L.D, L.sigma, L.trees, L.shift)).q
isone_groups(g) = length(g) == 1 && isempty(g[1][1]) && g[1][2] == (1 // 1, S.Mono(), nothing, false)

records = [r for r in load_records() if S.is_allosteric(r.c.mech) && S.admit(r.c.mech) === nothing]
tab = Dict{Any, Int}()
bump!(k) = (tab[k] = get(tab, k, 0) + 1)
for r in records
    m = r.c.mech
    A = S.allo_law(m)
    raw_same = isone(S.nterms(A.fA)) && isone(S.nterms(A.fI)) &&
               A.fA == one(typeof(A.fA)) && A.fI == one(typeof(A.fI))
    bump!(("cross weights both 1", raw_same))
    bump!(("cross weights monomials", S.nterms(A.fA) == 1 && S.nterms(A.fI) == 1))
    for (st, L) in (("A", A.active), ("I", A.inactive))
        L === nothing && (bump!(("inactive state empty",)); continue)
        rule, _ = S.normalizer(m, L)
        P = S.printed(m, L)
        f = L.trees[1] * L.shift
        g = flipgroups(m, L, f)
        bump!(("state", st, rule, P.N === nothing ? "N monomial" : "N apart",
               isone_groups(g) ? "flip weight = N" : "flip weight ≠ N"))
    end
    if A.inactive !== nothing
        fa = flipgroups(m, A.active, A.active.trees[1] * A.active.shift)
        fi = flipgroups(m, A.inactive, A.inactive.trees[1] * A.inactive.shift)
        bump!(("normalized flip weights equal (as groups)", fa == fi,
               "N apart: ", S.printed(m, A.active).N !== nothing, S.printed(m, A.inactive).N !== nothing))
    end
end
foreach(((k, v),) -> println(k, " => ", v), sort!(collect(tab); by = string))
