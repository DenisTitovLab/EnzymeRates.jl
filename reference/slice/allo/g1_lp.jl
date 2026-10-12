# ABOUTME: G1 two ways for the allosteric-stage checks: an exact rational LP (phase-1 simplex,
# ABOUTME: Bland's rule) on the joint system's rows, and the graph form; `S` is the Slice module.
using LinearAlgebra

const QB = Rational{BigInt}

"""Whether `{x ≥ 0 : M x = b}` is nonempty: phase 1 of the simplex method, Bland's rule,
exact rational arithmetic."""
function feasible(M::Matrix{QB}, b::Vector{QB})
    m, n = size(M)
    flip = [x < 0 ? -1 : 1 for x in b]
    T = hcat(flip .* M, Matrix{QB}(I, m, m), flip .* b)
    basis = collect(n + 1:n + m)
    cost = [zeros(QB, n); ones(QB, m)]
    while true
        d(j) = cost[j] - sum(cost[basis[i]] * T[i, j] for i in 1:m)
        j = findfirst(j -> d(j) < 0, 1:n + m)
        j === nothing && break
        rows = [i for i in 1:m if T[i, j] > 0]
        r = minimum(i -> T[i, end] / T[i, j], rows)
        i = argmin(i -> basis[i], [i for i in rows if T[i, end] / T[i, j] == r])
        T[i, :] ./= T[i, j]
        for k in 1:m
            k == i || iszero(T[k, j]) || (T[k, :] .-= T[k, j] .* T[i, :])
        end
        basis[i] = j
    end
    iszero(sum(cost[basis[i]] * T[i, end] for i in 1:m))
end

"""G1 as an LP on `constraint_system(m)`: every column fixed except the inactive `:OnlyA` ones;
a plain binding's K → ∞ (RE) or kon → 0 with koff fixed (SS); other steps: both SS rates → 0,
an RE K free. Feasible when some direction keeps every row."""
function g1_lp(m)
    cols, A, _ = S.constraint_system(m)
    fwd = S.orient(m)
    P, s, F = Int[], Int[], Int[]
    for (i, c) in enumerate(cols)
        c.inactive && m.tags[m.group[c.edge]] === :OnlyA || continue
        o = S.oriented(m, c.edge, fwd)
        bind = S.is_plain_binding(m, o)
        if o.re
            bind ? (push!(P, i); push!(s, -1)) : push!(F, i)
        elseif !bind || c.slot == 1
            push!(P, i); push!(s, 1)
        end
    end
    isempty(P) && return true
    AP, AF = QB.(A[:, P]) .* permutedims(s), QB.(A[:, F])
    feasible([AP AF -AF], -AP * ones(QB, length(P)))
end

"""G1 on the graph: no `:OnlyA` plain binding inside a component of the edges outside `:OnlyA`
groups, and no directed cycle of such bindings between components."""
function g1_graph(m)
    kept = BitVector(m.tags[g] !== :OnlyA for g in m.group)
    comp = collect(eachindex(m.forms))
    top(i) = (while comp[i] != i; i = comp[i]; end; i)
    foreach(k -> kept[k] && (comp[top(m.edges[k].a)] = top(m.edges[k].b)), eachindex(kept))
    comp .= top.(eachindex(comp))
    fwd = S.orient(m)
    arcs = [(comp[o.a], comp[o.b]) for o in (S.oriented(m, k, fwd) for k in findall(.!kept))
            if S.is_plain_binding(m, o)]
    any(a -> a[1] == a[2], arcs) && return false
    reach(u, seen = Set{Int}()) = (for (x, y) in arcs
                                       x == u && !(y in seen) && (push!(seen, y); reach(y, seen))
                                   end; seen)
    !any(((u, v),) -> u in reach(v), arcs)
end
