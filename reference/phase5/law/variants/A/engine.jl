# ABOUTME: The rate law of a mechanism: rapid-equilibrium segments contracted, spanning
# ABOUTME: trees of the segment graph, the edge-sum numerator, the solve and H divided out.

using LinearAlgebra

"""The most spanning-tree products a law may expand: segments times the spanning trees of
the segment graph, parallel steady-state steps counted apart (`complexity`)."""
const MAX_TERMS = 5000

"""
A rate law `v = E_total·H·Q/D` over `vars`: the concentrations (`metabolites(m)` order),
the independents of the solve (parameter order, both conformations) and `:Keq`.
`H = ΠS − ΠP/Keq` and `Q` has positive coefficients. `D = shift·Σ_g sigma[g]·trees[g]`:
`sigma[g]` sums the weights of RE segment `g`'s forms, `trees[g]` is its spanning-tree sum,
and `shift`, the concentration monomial that brings `Q` and `D` to lowest terms, is already
in `Q` and `D`.
"""
struct Law{T <: Integer}
    vars::Vars
    H::Poly{T}
    Q::Poly{T}
    D::Poly{T}
    sigma::Vector{Poly{T}}
    trees::Vector{Poly{T}}
    shift::Mono
end

"""A variable table holding the concentrations of `metabolites(m)`, in that order."""
function conc_vars(m::Mechanism)
    t = Vars()
    foreach(x -> var!(t, x; conc = true), metabolites(m))
    t
end

"""The product of the concentrations of `names` in `t`."""
conc_mono(t::Vars, names) = mono((t.index[x] => 1 for x in names)...)

"""
    segment_nodes(m, mask) -> (node, G, parent, order)

The segment graph of the edges in `mask`: `re_forest`'s `parent` and `order`, and each
form's node, its segment numbered over the forms the edges touch (0 for any other form).
"""
function segment_nodes(m::Mechanism, mask)
    segment, parent, order = re_forest(m, mask)
    touched, used = falses(length(m.forms)), falses(length(m.forms))
    foreach(e -> touched[e.a] = touched[e.b] = true, m.edges[mask])
    used[segment[touched]] .= true
    node = [touched[i] ? count(view(used, 1:segment[i])) : 0 for i in eachindex(segment)]
    node, maximum(node; init = 0), parent, order
end

"""
    spanning_trees(n, pairs) -> trees

The spanning trees of the graph on the nodes `1:n` with the edges `pairs`, each as a list of
edge positions (each edge included or excluded in turn, union-find against cycles).
"""
function spanning_trees(n, pairs)
    trees = Vector{Int}[]
    function grow(i, chosen, root)
        if length(chosen) == n - 1
            push!(trees, copy(chosen))
            return
        end
        length(pairs) - i + 1 < n - 1 - length(chosen) && return
        top(j) = (while root[j] != j; j = root[j]; end; j)
        a, b = top(pairs[i][1]), top(pairs[i][2])
        if a != b
            r = copy(root)
            r[a] = b
            push!(chosen, i)
            grow(i + 1, chosen, r)
            pop!(chosen)
        end
        grow(i + 1, chosen, root)
    end
    grow(1, Int[], collect(1:n))
    trees
end

"""The product of the weights `W[(child, parent)]` of the tree `t` (positions in `pairs`)
oriented toward node `g`."""
function toward(t, g, pairs, W::Dict{Tuple{Int, Int}, Poly{T}}) where {T}
    p, seen, todo = one(Poly{T}), Set(g), [g]
    while !isempty(todo)
        q = pop!(todo)
        for j in t
            u, v = pairs[j]
            c = u == q && !(v in seen) ? v : v == q && !(u in seen) ? u : 0
            c == 0 && continue
            p *= W[(c, q)]
            push!(seen, c)
            push!(todo, c)
        end
    end
    p
end

"""
    bare_law(m, mask = all edges, T = Int) -> (; vars, cols, num, den, sigma, trees)

The steady-state law of the edges in `mask` (at least one) with one variable per bare
constant: the concentrations (`conc_vars`), then variable `nc + j` for slot `cols[j][2]` of
edge `cols[j][1]`. `num/den` is the net consumption of `m.substrates[1]` per enzyme. The
rapid-equilibrium segments (`re_forest`: rooted at their lowest form, trees by content)
contract to nodes; a form weighs its root's weight times each tree edge's equilibrium
ratio. `den = Σ_g sigma[g]·trees[g]`: `sigma[g]` sums segment `g`'s weights and `trees[g]`
its spanning trees of the segment graph oriented toward it, each tree edge weighing the
steady-state flux out of its source. The numerator sums each steady-state edge's net flux
times the copies of `m.substrates[1]` it moves. An error when `complexity` exceeds
`MAX_TERMS`.
"""
function bare_law(m::Mechanism, mask = trues(length(m.edges)), ::Type{T} = Int) where {T}
    vars, fwd = conc_vars(m), orient(m)
    cols, at = Tuple{Int, Int}[], Dict{Tuple{Int, Int}, Int}()
    for k in findall(mask), s in 1:(m.edges[k].re ? 1 : 2)
        push!(cols, (k, s))
        at[(k, s)] = var!(vars, Symbol("c", k, "_", s))
    end
    con(k, s, e = 1) = mono(at[(k, s)] => e)
    cm(p) = mono((vars.index[x.name] => count_at(p, i) for (i, x) in enumerate(m.mets))...)
    node, G, parent, order = segment_nodes(m, mask)
    ratio(k, e = m.edges[k]) = mmul(con(k, 1, is_binding(e) || !fwd[k] ? -1 : 1),
                                    mmul(cm(e.cons), mpow(cm(e.rel), -1)))
    alpha = fill(Mono(), length(m.forms))
    for i in order
        k = parent[i]
        k == 0 && continue
        e = m.edges[abs(k)]
        alpha[i] = mmul(alpha[k > 0 ? e.a : e.b], mpow(ratio(abs(k)), sign(k)))
    end
    pairs, W = Tuple{Int, Int}[], Dict{Tuple{Int, Int}, Poly{T}}()
    flux, x = Tuple{Int, Int, Int, Mono, Mono}[], vars.index[m.substrates[1]]
    for k in findall(mask)
        e = m.edges[k]
        e.re && continue
        f = mmul(mmul(con(k, fwd[k] ? 1 : 2), cm(e.cons)), alpha[e.a])
        r = mmul(mmul(con(k, fwd[k] ? 2 : 1), cm(e.rel)), alpha[e.b])
        a, b = node[e.a], node[e.b]
        push!(flux, (Int(mexp(mmul(f, mpow(r, -1)), x)), a, b, f, r))
        a == b && continue
        minmax(a, b) in pairs || push!(pairs, minmax(a, b))
        addto!(get!(Poly{T}, W, (a, b)), f, 1)
        addto!(get!(Poly{T}, W, (b, a)), r, 1)
    end
    complexity(m, mask) > MAX_TERMS &&
        error("the rate law of this mechanism expands more than $MAX_TERMS " *
              "spanning-tree products; laws this large take very long to compile and are " *
              "unlikely to be useful for fitting")
    trees = spanning_trees(G, pairs)
    G > 1 && isempty(trees) &&
        error("the steps do not join every rapid-equilibrium segment")
    D = [sum((toward(t, g, pairs, W) for t in trees); init = zero(Poly{T})) for g in 1:G]
    sigma = [zero(Poly{T}) for _ in 1:G]
    foreach(i -> node[i] > 0 && addto!(sigma[node[i]], alpha[i], 1), eachindex(m.forms))
    num = zero(Poly{T})
    for (w, a, b, f, r) in flux
        w == 0 || (num += scale(D[a] * f - D[b] * r, w))
    end
    (; vars, cols, num, den = sum(sigma .* D; init = zero(Poly{T})), sigma, trees = D)
end

"""
    law(m, sol = solve(m), mask = all edges; inactive = false) -> Law

The rate law of the edges in `mask`, in the inactive conformation's constants when
`inactive`, over the independents of `sol`: the bare law with each bare constant replaced by
its column's value (an independent, or a dependent's power product of independents and
`Keq`), the numerator divided exactly by `H`, then brought to lowest terms. Computed in
checked `Int`, retried in `BigInt` on overflow.
"""
function law(m::Mechanism, sol::Solution = solve(m), mask = trues(length(m.edges));
             inactive::Bool = false)
    try
        widen_law(law(Int, m, sol, mask, inactive))
    catch err
        err isa OverflowError || rethrow()
        law(BigInt, m, sol, mask, inactive)
    end
end

"""`L` with every coefficient a `BigInt`, the one law type callers see."""
widen_law(L::Law) = Law(L.vars, convert(Poly{BigInt}, L.H), convert(Poly{BigInt}, L.Q),
                        convert(Poly{BigInt}, L.D), [convert(Poly{BigInt}, p) for p in L.sigma],
                        [convert(Poly{BigInt}, p) for p in L.trees], L.shift)

function law(::Type{T}, m::Mechanism, sol::Solution, mask, inactive::Bool) where {T}
    b, vars, names = bare_law(m, mask, T), conc_vars(m), column_names(m, sol)
    value = Dict(i => mono(var!(vars, names[i]) => 1) for i in sol.indep)
    keq, nc = var!(vars, :Keq), length(metabolites(m))
    for d in sol.deps
        value[d.col] = foldl(mmul, (mpow(value[f], e) for (f, e) in d.terms);
                             init = mono(keq => d.keq))
    end
    sub(p) = subst(p, v -> v <= nc ? nothing :
                           value[sol.index[Col(b.cols[v - nc]..., inactive)]])
    H = term(T, conc_mono(vars, m.substrates)) -
        term(T, mmul(conc_mono(vars, m.products), mono(keq => -1)))
    Q, D = divexact(sub(b.num), H), sub(b.den)
    shift = conc_shift((Q, D), vars)
    Law(vars, H, Q * shift, D * shift, map(sub, b.sigma), map(sub, b.trees), shift)
end

"""
    complexity(m, mask = all edges) -> BigInt

The size of the law of the edges in `mask` before kinetic grouping, today's
`_eq_complexity`: the number of segments times the spanning trees of the segment graph
counted with parallel steady-state edges (Kirchhoff's matrix-tree theorem).
"""
function complexity(m::Mechanism, mask = trues(length(m.edges)))
    node, G = segment_nodes(m, mask)
    G <= 1 && return big(G)
    lap = zeros(BigInt, G, G)
    for e in m.edges[mask]
        a, b = node[e.a], node[e.b]
        (e.re || a == b) && continue
        lap[a, a] += 1
        lap[b, b] += 1
        lap[a, b] -= 1
        lap[b, a] -= 1
    end
    G * det(lap[2:end, 2:end])
end
