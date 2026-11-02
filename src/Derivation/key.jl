# ABOUTME: The canonical key of a rate law, by colour refinement and individualization,
# ABOUTME: and the rescaling certificate that flags a one-subunit L phantom.

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

The polynomials of law `L` in role order, the units (`H`, fixed by the reaction, alone;
then each group of components that rescale together without changing the rate), and a
header naming the law's shape: `H`, then `(Q, D)` for a plain law; for an MWC law `H`, each
conformation's `(Q, D, free-E weight)` and each regulatory site's binding polynomials in
both conformations, the header naming the catalytic and the site multiplicities.
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
        # Each comprehension captures `let` rebindings of the colours the loop reassigns,
        # so none is boxed and every vector `_ranks` sorts has a concrete element type.
        t2 = let pcol = pcol, tcol = tcol
            _ranks([(tcol[t], sort!([(pcol[p], e) for (p, e) in inc[t]]))
                    for t in eachindex(inc)])
        end
        p2 = let pcol = pcol
            _ranks([(pcol[p], sort!([(t2[t], e) for (t, e) in occ[p]]))
                    for p in eachindex(occ)])
        end
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
`mechanism_id` may be the mechanism itself, whose id (the first 16 hex digits of the
SHA-256 of its text, `show`) is then taken only at the cap.
"""
function _canonical(L::Union{_Law, _AlloLaw},
                    mechanism_id::Union{AbstractString, _AnyMech})
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
        text = "mechanism|" * (mechanism_id isa AbstractString ? mechanism_id :
                               bytes2hex(sha256(sprint(show, mechanism_id)))[1:16])
        return _Canonical(vars.names[params], text, bytes2hex(sha256(text)), st.leaves,
                          true)
    end
    absent = sort!(vars.names[setdiff(params, present)])
    text = string(header, "|P=", length(params), "|absent=", length(absent), "\n", st.text)
    _Canonical([vars.names[present[st.order]]; absent], text, bytes2hex(sha256(text)),
               st.leaves, false)
end

# ─── One-subunit L phantoms ─────────────────────────────────────────────────

"""
    _phantom_weights(A) -> w

The rescaling certificate of the MWC law `A` (spec section 13, "One-subunit L phantoms").
At one catalytic subunit, with an inactive conformation that turns nothing over
(`QI = 0`), the law is `v = N/(Y + L·X)` with `N = H·fI·QA·WA`, `Y = fI·DA·WA` and
`X = fA·DI·WI`. When every term of `X` is a term of `Y` with the same coefficient,
`Y + L·X = (Y − X) + (1 + L)·X`; weights `w` over the fitted constants other than `L`, and
a shift `s`, with `e·w = s` on the constant exponents `e` of every term of `N` and of
`Y − X` and `e·w = s + 1` on those of every term of `X`, then give
`v(θ, L) = v₀(θ ⊙ (1 + L)^w)` for all positive `θ` and `L`, `v₀ = N/Y` being the law at
`L = 0`: the two laws have one family. Returns `w` in `A.vars` order (the independents,
then the regulator constants), or an empty vector when no weights exist, when `A` has
another shape, or when the arithmetic leaves `Int64`, so that such a variant is fitted.
"""
function _phantom_weights(A::_AlloLaw)
    none = Rational{Int}[]
    A.n == 1 && iszero(A.inactive.Q) || return none
    try
        WA = foldl(*, (w for (w, _, k) in A.sites for _ in 1:k); init = one(_Poly))
        WI = foldl(*, (w for (_, w, k) in A.sites for _ in 1:k); init = one(_Poly))
        N = A.H * A.fI * A.active.Q * WA
        Y = A.fI * A.active.D * WA
        X = A.fA * A.inactive.D * WI
        all(((m, c),) -> get(Y.terms, m, 0) == c, X.terms) || return none
        vars = A.vars
        θ = [v for v in eachindex(vars.names)
             if !vars.conc[v] && vars.names[v] !== :Keq && vars.names[v] !== :L]
        at = Dict(v => j for (j, v) in enumerate(θ))
        n = length(θ)
        # One row per distinct exponent vector: [e; -1] · [w; s] = rhs.
        rows = Dict{Vector{Rational{Int}}, Rational{Int}}()
        function row!(p::_Poly, rhs::Rational{Int})
            for m in keys(p.terms)
                e = zeros(Rational{Int}, n + 1)
                for (v, x) in m
                    j = get(at, v, 0)
                    j == 0 || (e[j] += x)
                end
                e[end] = -1
                get!(rows, e, rhs) == rhs || return false
            end
            true
        end
        row!(N, 0 // 1) && row!(Y - X, 0 // 1) && row!(X, 1 // 1) || return none
        M = zeros(Rational{Int}, length(rows), n + 2)
        for (i, (e, rhs)) in enumerate(rows)
            M[i, 1:(n + 1)] = e
            M[i, n + 2] = rhs
        end
        pivots = _rref!(M, n + 1)
        all(iszero, view(M, (length(pivots) + 1):size(M, 1), n + 2)) || return none
        w = zeros(Rational{Int}, n + 1)
        for (r, c) in enumerate(pivots)
            w[c] = M[r, n + 2]
        end
        w[1:n]
    catch err
        err isa OverflowError || rethrow()
        none
    end
end

"""
    _is_phantom(A) -> Bool

Whether a mechanism whose MWC law is `A` is a one-subunit L phantom (spec section 13): `A`
has a rescaling certificate (`_phantom_weights`), which requires one catalytic subunit and
an inactive conformation that turns nothing over, so the mechanism's family is the family
of its law at `L = 0`.
"""
_is_phantom(A::_AlloLaw) = !isempty(_phantom_weights(A))

# ─── The key a mechanism is fitted under ───────────────────────────────────

"""
    rate_equation_dedup_key(m) -> (; key, parameters, fit_key, fit_parameters, capped)

The rate-law key of the well-posed mechanism `m` and the key of the fit it reads, from its
law (through the memo, `_Memo`). `key` and `parameters` are its law's canonical key and its
fitted constants in canonical order (`_canonical`); `capped` says that the labelling
reached the leaf cap, so the law took the key of `m`'s text id (`_canonical`).
`fit_key` and `fit_parameters` name the law whose fit `m` reads: its own, except for an L
phantom (`_is_phantom`), which reads its active law, its law at `L = 0`, when that law
holds every fitted constant of `m` but `L`. Two mechanisms with one key share one fit.
"""
function rate_equation_dedup_key(m::_AnyMech)
    check_valid_rate_equation(m)
    _with_memo_law(m) do m, s, L
        own = _canonical(L, m)
        fit = own
        # An L phantom reads the fit of its law at L = 0, its active law, when that law
        # holds every fitted constant but L. A regulatory site whose ligands are all
        # :EqualAI has one factor in both conformations, so the certificate admits it,
        # yet the active law lacks its constants: such a phantom reads its own fit, so a
        # row never reads names its fit lacks.
        if m isa Mechanisms.AllostericMechanism && _is_phantom(L)
            active = _canonical(L.active, m)
            Set(active.order) == Set(s.fitted[1:(end - 1)]) && (fit = active)
        end
        (; key = own.key, parameters = own.order, fit_key = fit.key,
         fit_parameters = fit.order, capped = own.capped)
    end
end
