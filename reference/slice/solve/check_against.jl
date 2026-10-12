# ABOUTME: Compares the solve's `against` rule (against the net reaction) with the literal
# ABOUTME: "slot 2" rule: forward Haldane dependents, dependent sets, and label statistics.
include(joinpath(@__DIR__, "lib.jl"))

const SHARED = ("specs", "handwritten", "init_uniuni", "init_unibi", "init_bibi",
                "init_bibi_pp", "expand_bibi")
records = [r for r in load_records(SHARED) if !S.is_allosteric(r.c.mech)]

"""Whether column `c` drives with the net reaction (a forward rate constant), by the solve's
rule: the reverse of `against`."""
function forward(m, fwd, c)
    o = S.oriented(m, c.edge, fwd)
    o.re && return false
    takes = S.is_binding(o) && any(i -> S.count_at(o.cons, i) > 0 && !m.mets[i].modifier &&
                                        m.mets[i].name in m.products, eachindex(m.mets))
    (c.slot == 2) == takes
end

"""The Haldane-labelled dependent columns (as `Col`s) when the columns are eliminated in the
order `cols2` instead of the solve's."""
function haldane_cols(m, cols, A, b, cols2)
    perm = [findfirst(==(c), cols) for c in cols2]
    M = Rational{BigInt}.([A[:, perm] b])
    n = length(cols2)
    piv = S.rref!(M, n)
    rep = Dict(S.class(m, c) => i for (i, c) in enumerate(cols2))
    [cols2[p] for (r, p) in enumerate(piv) if !iszero(M[r, n + 1]) && rep[S.class(m, cols2[p])] == p]
end

fwd_mine, fwd_lit, same_set, two_haldane, labels = 0, 0, 0, 0, Dict{Symbol, Int}()
for r in records
    m = r.c.mech
    fwd = S.orient(m)
    cols, A, b = S.constraint_system(m)
    sol = S.solve(m)
    foreach(d -> labels[d.label] = get(labels, d.label, 0) + 1, sol.deps)
    mine = [sol.cols[d.col] for d in sol.deps if d.label === :haldane]
    count(d -> d.label === :haldane, sol.deps) >= 2 && (global two_haldane += 1)
    # the literal rule: against = slot 2, the rest of the key unchanged
    keys = map(cols) do c
        o = S.oriented(m, c.edge, fwd)
        kind = o.cons == 0 && o.rel == 0 ? 3 : S.is_plain_binding(m, o) ? 1 : 2
        occupied = m.forms[o.a].bound != 0 && m.forms[o.b].bound != 0
        (kind, occupied, !o.re, c.slot == 2, o.a, o.b, o.cons, o.rel)
    end
    lit = haldane_cols(m, cols, A, b, cols[sortperm(keys; rev = true)])
    any(c -> forward(m, fwd, c), mine) && (global fwd_mine += 1)
    any(c -> forward(m, fwd, c), lit) && (global fwd_lit += 1)
    Set(mine) == Set(lit) && (global same_set += 1)
end
@printf("%d plain mechanisms of the shared populations\n", length(records))
@printf("a Haldane dependent is a forward rate constant: solve's rule %d, literal slot-2 rule \
        %d; the two rules pick the same Haldane dependents in %d\n", fwd_mine, fwd_lit,
        same_set)
@printf("mechanisms with two or more Haldane dependents: %d\n", two_haldane)
println("dependent labels over all of them: ", sort!(collect(labels)))
