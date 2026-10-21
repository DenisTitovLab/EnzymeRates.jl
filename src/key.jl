# ABOUTME: The dedup key: a rate law brought to lowest terms, its fitted constants labelled
# ABOUTME: canonically by colour refinement and individualization, and that text's SHA-256.

"""
A rate law's canonical form (spec section 10): `order`, its fitted constants in canonical
order, the constants no term holds last by name (`L`, a fixed variable, is not listed);
`text`, the canonical text; `key`, the text's SHA-256 in 64 hex digits; `leaves`, the
labellings the search compared; `capped`, whether the search stopped at `_KEY_LEAF_CAP`,
in which case `text` is the mechanism's own (`"mechanism|"` and its id) and `order` is the
variable order. Below the cap, two laws get one key exactly when one is the other with its
fitted constants renamed, each brought to lowest terms unit by unit; concentrations count
by metabolite name.
"""
struct _Canonical
    order::Vector{Symbol}
    text::String
    key::String
    leaves::Int
    capped::Bool
end

"""The most labellings the key's search compares before the law takes its mechanism's
key. Enumerated laws need one; the cap bounds the search on laws with many symmetric
constants."""
const _KEY_LEAF_CAP = 100

"""
    _key_parts(L) -> (components, units, header)

The polynomials of law `L` in role order, the units of components that rescale together
without changing the rate, and a header naming the law's shape: `H`, then `(Q, D)` for a
plain law; for an MWC law `H`, each conformation's `(Q, D, free-E weight)` and each
regulatory site's binding polynomials in both conformations, the header naming the
catalytic and the site multiplicities.
"""
_key_parts(L::_Law) = (_Poly[L.H, L.Q, L.D], [[1], [2, 3]], "plain")

function _key_parts(A::_AlloLaw)
    comps = _Poly[A.H, A.active.Q, A.active.D, A.fA, A.inactive.Q, A.inactive.D, A.fI]
    units = [[1], [2, 3, 4], [5, 6, 7]]
    for (WA, WI, _) in A.sites
        push!(comps, WA, WI)
        push!(units, [length(comps) - 1, length(comps)])
    end
    comps, units, string("mwc|n=", A.n, "|sites=", join((s[3] for s in A.sites), ","))
end

"""The rank of each element of `xs` among its distinct values, 1 for the smallest."""
_ranks(xs::AbstractVector) = (u = sort!(unique(xs)); [searchsortedfirst(u, x) for x in xs])

"""Colour refinement of the constants and the terms from the colours `pcol` and `tcol`,
over the incidences `inc[t] = [(constant, exponent)]` and `occ[p] = [(term, exponent)]`:
the coarsest stable ordered partition below the given one."""
function _refine(pcol::Vector{Int}, tcol::Vector{Int},
                 inc::Vector{Vector{Tuple{Int, _Exp}}},
                 occ::Vector{Vector{Tuple{Int, _Exp}}})
    pcol, tcol = _ranks(pcol), _ranks(tcol)
    while true
        t2 = _ranks([(tcol[t], sort!([(pcol[p], e) for (p, e) in inc[t]]))
                     for t in eachindex(inc)])
        p2 = _ranks([(pcol[p], sort!([(t2[t], e) for (t, e) in occ[p]]))
                     for p in eachindex(occ)])
        maximum(p2; init = 0) == maximum(pcol; init = 0) &&
            maximum(t2; init = 0) == maximum(tcol; init = 0) && return p2, t2
        pcol, tcol = p2, t2
    end
end

"""The terms as sorted lines of text, constant `p` written `x<pcol[p]>`."""
function _encode(pcol::Vector{Int}, inc::Vector{Vector{Tuple{Int, _Exp}}},
                 labels::Vector{String})
    factors(t) = join(("x$(pcol[p])^$e" for (p, e) in sort(inc[t]; by = x -> pcol[x[1]])),
                      "*")
    join(sort!([string(labels[t], "|", factors(t)) for t in eachindex(inc)]), "\n")
end

"""The key search's state: the smallest text found, the constants' order in it, and the
labellings compared."""
mutable struct _KeySearch
    text::String
    order::Vector{Int}
    leaves::Int
end

"""Refine the colours `pcol`, `tcol`; at a labelling, keep it in `st` when its text is the
smallest so far; else individualize each constant of the first tied class in turn and
search on, until `_KEY_LEAF_CAP` labellings."""
function _search!(st::_KeySearch, pcol::Vector{Int}, tcol::Vector{Int}, inc, occ,
                  labels::Vector{String})
    st.leaves >= _KEY_LEAF_CAP && return st
    pcol, tcol = _refine(pcol, tcol, inc, occ)
    if allunique(pcol)
        st.leaves += 1
        text = _encode(pcol, inc, labels)
        if st.leaves == 1 || text < st.text
            st.text, st.order = text, sortperm(pcol)
        end
        return st
    end
    c = minimum(x for x in pcol if count(==(x), pcol) > 1)
    for p in findall(==(c), pcol)
        q = 2 .* pcol
        q[p] -= 1
        _search!(st, q, tcol, inc, occ, labels)
    end
    st
end

"""
    _canonical(L, mechanism_id) -> _Canonical

The canonical form of the rate law `L` (a `_Law` or an `_AlloLaw`). Each unit of
`_key_parts` is multiplied by the monomial that brings every variable to lowest terms
(`_lowest_shift`) and divided by its integer content. A term is labelled by its
component, coefficient and fixed factors (concentrations as `c_<name>`, `Keq`, `L`),
sorted. The fitted constants, every other variable a term holds, are coloured by
refinement on the term–constant incidence; while colours tie, each constant of the first
tied class is individualized in turn, keeping the smallest text (`_search!`). At
`_KEY_LEAF_CAP` labellings the law takes the key of `"mechanism|" * mechanism_id`.
"""
function _canonical(L::Union{_Law, _AlloLaw}, mechanism_id::AbstractString)
    comps, units, header = _key_parts(L)
    vars = L.vars
    fixed(v) = vars.conc[v] || vars.names[v] === :Keq || vars.names[v] === :L
    terms = Tuple{Int, _Mono, Int}[]
    for u in units
        s = _lowest_shift(comps[u])
        g = foldl(gcd, (c for i in u for c in values(comps[i].terms)); init = 0)
        for i in u, (mo, c) in comps[i].terms
            push!(terms, (i, _mmul(mo, s), div(c, g)))
        end
    end
    params = [v for v in eachindex(vars.names) if !fixed(v)]
    held = Set(v for (_, mo, _) in terms for (v, _) in mo)
    present = filter(in(held), params)
    at = Dict(v => j for (j, v) in enumerate(present))
    label(v) = vars.conc[v] ? "c_" * String(vars.names[v]) : String(vars.names[v])
    labels = [string(i, "|", c, "|",
                     join(sort!([string(label(v), "^", e) for (v, e) in mo if fixed(v)]),
                          "*"))
              for (i, mo, c) in terms]
    inc = [Tuple{Int, _Exp}[(at[v], e) for (v, e) in mo if !fixed(v)]
           for (_, mo, _) in terms]
    occ = [Tuple{Int, _Exp}[] for _ in present]
    for (t, ps) in enumerate(inc), (p, e) in ps
        push!(occ[p], (t, e))
    end
    st = _search!(_KeySearch("", Int[], 0), ones(Int, length(present)), _ranks(labels),
                  inc, occ, labels)
    if st.leaves >= _KEY_LEAF_CAP
        text = "mechanism|" * mechanism_id
        return _Canonical(vars.names[params], text, bytes2hex(sha256(text)), st.leaves,
                          true)
    end
    absent = sort!(vars.names[setdiff(params, present)])
    text = string(header, "|P=", length(params), "|absent=", length(absent), "\n", st.text)
    _Canonical([vars.names[present[st.order]]; absent], text, bytes2hex(sha256(text)),
               st.leaves, false)
end
