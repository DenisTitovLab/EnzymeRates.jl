# ABOUTME: kcat, the saturating turnover: the law of the mechanism without its dead-end
# ABOUTME: branches, products at zero, substrates growing, the maximum over modifier sets.

"""
    catalytic_core(m) -> BitVector

The edges outside `m`'s dead-end branches: biconnected blocks that hold no cycle turning the
reaction over (`unbalanced_blocks`) are removed while one meets the rest at most at one
form.
"""
function catalytic_core(m::Mechanism)
    sg = [x.name in m.substrates ? 1 : x.name in m.products ? -1 : 0 for x in m.mets]
    up(p) = sum(sg[i] * count_at(p, i) for i in eachindex(m.mets))
    block, turns = unbalanced_blocks(length(m.forms),
                                     [(Int(e.a), Int(e.b)) for e in m.edges],
                                     [up(e.cons) - up(e.rel) for e in m.edges])
    keep = trues(length(m.edges))
    while true
        at = [unique(block[k] for k in findall(keep) if i in (m.edges[k].a, m.edges[k].b))
              for i in eachindex(m.forms)]
        cuts = [count(s -> length(s) > 1 && b in s, at) for b in 1:maximum(block)]
        leaf = [k for k in findall(keep) if !(block[k] in turns) && cuts[block[k]] <= 1]
        isempty(leaf) && return keep
        keep[leaf] .= false
    end
end

"""The error kcat raises when there is no kcat: `what` names the problem."""
kcat_error(what) = error("kcat: $what; fit in absolute mode (scale_k_to_kcat = nothing)")

"""
    face_ratio(n, d, growing) -> (numerator, denominator)

kcat's candidate from the leading faces `n` and `d` (concentration monomial => coefficient)
of the numerator and the denominator: their one ratio, or an error when the ratio differs
across the face, so that the limit depends on how the `growing` names grow.
"""
function face_ratio(n, d, growing)
    f = first(sort!(collect(keys(d))))
    nf(mo) = get(n, mo, zero(typeof(d[f])))
    keys(n) ⊆ keys(d) && all(mo -> nf(mo) * d[f] == nf(f) * d[mo], keys(d)) ||
        kcat_error("the saturating rate depends on how $(join(growing, ", ")) grow " *
                   "relative to each other")
    nf(f), d[f]
end

"""
    kcat_terms(m, sol = solve(m), mask = catalytic_core(m)) -> (vars, terms)

kcat's candidates as `(numerator, denominator)` polynomials over `vars` in the constants.
On the law of the edges in `mask`, for each set of modifiers (names that bind only as
modifiers) that one of its terms holds: the limit of `v/E_total` with the products and the
other modifiers at zero as the substrates grow, then that set's modifiers. Each stage grows
its concentrations together, so the limit is the ratio of the numerator and denominator
coefficients on the denominator's leading face (`face_ratio`). The maximum over these sets
is the maximum over every corner, each modifier at zero or saturating: a corner's leading
face holds a term whose own set gives the same limit. An error names the problem when the
rate grows without bound or its limit depends on how the concentrations grow.
"""
function kcat_terms(m::Mechanism, sol::Solution = solve(m), mask = catalytic_core(m))
    L, nc = law(m, sol, mask), length(metabolites(m))
    names = L.vars.names
    subs = [i for i in 1:nc if names[i] in m.substrates]
    N = L.Q * conc_mono(L.vars, m.substrates)
    held(mo) = [v for (v, _) in mo if v <= nc && !is_reactant(m, names[v])]
    terms = map(sort!(unique(held(mo) for p in (N, L.D) for mo in keys(p.terms)))) do sat
        off = [i for i in 1:nc if !(i in subs || i in sat)]
        kept(p) = split_conc(Poly(filter(t -> !any(v -> v.first in off, t.first), p.terms)),
                             L.vars)
        n, d = kept(N), kept(L.D)
        isempty(d) && kcat_error("the rate is undefined (0/0) with the products at zero")
        deg(mo) = (sum(e for (v, e) in mo if v in subs; init = zero(Exp)),
                   sum(e for (v, e) in mo if v in sat; init = zero(Exp)))
        top = maximum(deg, keys(d))
        for mo in keys(n)
            deg(mo) > top || continue
            kcat_error("the rate grows without bound with " *
                       join(names[deg(mo)[1] > top[1] ? subs : sat], ", ") *
                       ", so the mechanism never saturates")
        end
        face(x) = filter(t -> deg(t.first) == top, x)
        face_ratio(face(n), face(d), names[[subs; sat]])
    end
    L.vars, terms
end
