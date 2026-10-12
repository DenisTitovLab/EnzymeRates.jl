# ABOUTME: The dedup key (K1, decision 10): a rate law brought to lowest terms, its fitted
# ABOUTME: constants labelled by colour refinement, and the SHA-256 of that canonical text.

using SHA

"""
A rate law's canonical form: its fitted constants in canonical order (any that no term holds
last, by name), its canonical text, and the key, the text's SHA-256. Two laws get one key
when one is the other with its fitted constants renamed, each in lowest terms (a common
factor that is no monomial stays). Concentrations count by their position in
`metabolites(m)`, so renaming metabolites without changing their order keeps the key.
"""
struct Canonical
    order::Vector{Symbol}
    text::String
    key::String
end

"""
    key_parts(L) -> (components, units, header)

The polynomials of law `L` in role order, the units of components that rescale together
without changing the rate, and a header for the law's shape: `H` alone and `(Q, D)` for a
plain law; for an MWC law `H`, each conformation's `(Q, D, flip weight)` and each regulatory
site's polynomials in both conformations, the header naming the multiplicities.
"""
key_parts(L::Law) = [L.H, L.Q, L.D], [[1], [2, 3]], "plain"

function key_parts(A::AlloLaw)
    comps = [A.H, A.active.Q, A.active.D, A.fA, inactive_terms(A)..., A.fI]
    units = [[1], [2, 3, 4], [5, 6, 7]]
    for (WA, WI, _) in A.sites
        push!(comps, WA, WI)
        push!(units, [length(comps) - 1, length(comps)])
    end
    comps, units, string("mwc|n=", A.n, "|sites=", join((s[3] for s in A.sites), ","))
end

"""The rank of each element of `xs` among its distinct values (1 for the smallest)."""
ranks(xs) = (u = sort!(unique(xs)); [searchsortedfirst(u, x) for x in xs])

"""Colour refinement of constants and terms, from the colours `pcol` and `tcol`, over the
incidences `inc[t] = [(constant, exponent)]` and `occ[p] = [(term, exponent)]`: the coarsest
stable ordered partition below the given one."""
function refine(pcol, tcol, inc, occ)
    pcol, tcol = ranks(pcol), ranks(tcol)
    while true
        t2 = ranks([(tcol[t], sort!([(pcol[p], e) for (p, e) in inc[t]]))
                    for t in eachindex(inc)])
        p2 = ranks([(pcol[p], sort!([(t2[t], e) for (t, e) in occ[p]]))
                    for p in eachindex(occ)])
        maximum(p2; init = 0) == maximum(pcol; init = 0) &&
            maximum(t2; init = 0) == maximum(tcol; init = 0) && return p2, t2
        pcol, tcol = p2, t2
    end
end

"""The terms as sorted lines of text, constant `p` written `x<pcol[p]>`."""
encode(pcol, inc, labels) =
    join(sort!([string(labels[t], "|", join(("x$(pcol[p])^$e" for (p, e) in
                                             sort(inc[t]; by = x -> pcol[x[1]])), "*"))
                for t in eachindex(inc)]), "\n")

"""
    canonical(L) -> Canonical

The canonical form of the rate law `L` (a `Law` or an `AlloLaw`). Each unit of `key_parts`
is multiplied by the monomial that brings every variable to lowest terms (`lowest_shift`)
and divided by its integer content. A term is labelled by its component, coefficient and
fixed factor (concentrations by position, `Keq`, `L`). The fitted constants, every other
variable a term holds, are coloured by refinement on the term–constant incidence; while
colours tie, each constant of the first tied class is individualized in turn, keeping the
smallest text.
"""
function canonical(L)
    comps, units, header = key_parts(L)
    vars = L.vars
    fixed(v) = vars.conc[v] || vars.names[v] in (:Keq, :L)
    terms = Tuple{Int, Mono, BigInt}[]
    for u in units
        s = lowest_shift(comps[u])
        g = foldl(gcd, (c for i in u for c in values(comps[i].terms)); init = big(0))
        append!(terms, ((i, mmul(mo, s), div(big(c), g)) for i in u
                        for (mo, c) in comps[i].terms))
    end
    params = [v for v in eachindex(vars.names) if !fixed(v)]
    held = Set(v for (_, mo, _) in terms for (v, _) in mo)
    present = filter(in(held), params)
    at = Dict(v => j for (j, v) in enumerate(present))
    name(v) = vars.conc[v] ? "c$v" : string(vars.names[v])
    labels = [string(i, "|", c, "|", join(("$(name(v))^$e" for (v, e) in mo if fixed(v)),
                                          "*")) for (i, mo, c) in terms]
    inc = [[(at[v], e) for (v, e) in mo if !fixed(v)] for (_, mo, _) in terms]
    occ = [Tuple{Int, Exp}[] for _ in present]
    for (t, ps) in enumerate(inc), (p, e) in ps
        push!(occ[p], (t, e))
    end
    best, order, leaves = nothing, Int[], 0
    function search(pcol, tcol)
        pcol, tcol = refine(pcol, tcol, inc, occ)
        if allunique(pcol)
            (leaves += 1) > 100_000 && error("canonical labelling: over 100000 leaves")
            text = encode(pcol, inc, labels)
            (best === nothing || text < best) && ((best, order) = (text, sortperm(pcol)))
            return
        end
        c = minimum(x for x in pcol if count(==(x), pcol) > 1)
        for p in findall(==(c), pcol)
            q = 2 .* pcol
            q[p] -= 1
            search(q, tcol)
        end
    end
    search(ones(Int, length(present)), ranks(labels))
    absent = sort!(vars.names[setdiff(params, present)])
    text = string(header, "|P=", length(params), "|absent=", length(absent), "\n", best)
    Canonical([vars.names[present[order]]; absent], text, bytes2hex(sha256(text)))
end
