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

"""
    kcat_terms(m, sol = solve(m), mask = catalytic_core(m)) -> (vars, terms)

kcat's candidates as `(numerator, denominator)` polynomials over `vars` in the constants.
On the law of the edges in `mask`, for each set of modifiers (names that bind only as
modifiers) that one of its terms holds: the limit of `v/E_total` with the products and the
other modifiers at zero as the substrates grow, then that set's modifiers. Each stage grows
its concentrations together, so the limit is the ratio of the numerator and denominator
coefficients on the denominator's leading face. The maximum over these sets is the maximum
over every corner, each modifier at zero or saturating: a corner's leading face holds a term
whose own set gives the same limit. An error names the problem when the rate grows without
bound or its limit depends on how the concentrations grow relative to each other.
"""
function kcat_terms(m::Mechanism, sol::Solution = solve(m), mask = catalytic_core(m))
    L, nc = law(m, sol, mask), length(metabolites(m))
    names, absolute = L.vars.names, "; fit in absolute mode (scale_k_to_kcat = nothing)"
    subs = [i for i in 1:nc if names[i] in m.substrates]
    N = L.Q * conc_mono(L.vars, m.substrates)
    held(mo) = [v for (v, _) in mo if v <= nc && !is_reactant(m, names[v])]
    terms = map(sort!(unique(held(mo) for p in (N, L.D) for mo in keys(p.terms)))) do sat
        off = [i for i in 1:nc if !(i in subs || i in sat)]
        kept(p) = split_conc(Poly(filter(t -> !any(v -> v.first in off, t.first), p.terms)),
                             L.vars)
        n, d = kept(N), kept(L.D)
        isempty(d) && error("kcat: the rate is undefined (0/0) with the products at zero" *
                            absolute)
        deg(mo) = (sum(e for (v, e) in mo if v in subs; init = zero(Exp)),
                   sum(e for (v, e) in mo if v in sat; init = zero(Exp)))
        top = maximum(deg, keys(d))
        face = sort!([mo for mo in keys(d) if deg(mo) == top])
        for mo in keys(n)
            deg(mo) > top || continue
            error("kcat: the rate grows without bound with " *
                  join(names[deg(mo)[1] > top[1] ? subs : sat], ", ") *
                  ", so the mechanism never saturates" * absolute)
        end
        nf(mo) = get(n, mo, zero(typeof(d[mo])))
        f = first(face)
        all(mo -> nf(mo) * d[f] == nf(f) * d[mo], face) &&
            all(mo -> deg(mo) < top || haskey(d, mo), keys(n)) ||
            error("kcat: the saturating rate depends on how " *
                  join(names[[subs; sat]], ", ") * " grow relative to each other" *
                  absolute)
        (nf(f), d[f])
    end
    L.vars, terms
end

"""kcat at the parameter values `params` (by name, as `parameters(m)`): the largest
candidate of `kcat_terms`."""
function kcat(m::Mechanism, params::NamedTuple)
    vars, terms = kcat_terms(m)
    T = typeof(float(params.Keq))
    x = T[vars.conc[i] ? one(T) : params[s] for (i, s) in enumerate(vars.names)]
    maximum(evaluate(n, x) / evaluate(d, x) for (n, d) in terms)
end
