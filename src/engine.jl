# ABOUTME: Rate laws on edge masks: one conformation's law (segments, spanning trees, H
# ABOUTME: divided out), the MWC law on both conformations, V×τ and the summed term cap.

"""
The rate law of one conformation, `v = E_total·H·Q/D`, over `vars`: the concentrations,
then the solve's independent constants (`_Solve.independents` order), then `:Keq`.
`H = ΠS − ΠP/Keq`, and `Q` has positive coefficients. `D = shift·Σ_g sigma[g]·trees[g]`:
`sigma[g]` sums the weights of rapid-equilibrium segment `g`'s forms, `trees[g]` is its
spanning-tree sum, and `shift`, the concentration monomial that brings `Q` and `D` to
lowest terms, is already in `Q` and `D`. Two laws are equal when their fields are.
"""
mutable struct _Law
    const vars::_Vars
    const H::_Poly
    const Q::_Poly
    const D::_Poly
    const sigma::Vector{_Poly}
    const trees::Vector{_Poly}
    const shift::_Mono
end
Base.:(==)(a::_Law, b::_Law) =
    a.vars == b.vars && a.H == b.H && a.Q == b.Q && a.D == b.D && a.sigma == b.sigma &&
    a.trees == b.trees && a.shift == b.shift
Base.hash(L::_Law, h::UInt) = hash((L.vars, L.H, L.Q, L.D), hash(:_Law, h))

"""
    _spanning_trees(n, pairs) -> trees

The spanning trees of the graph on the nodes `1:n` with the edges `pairs`, each as a list
of positions in `pairs`: each edge is included or left out in turn, and a union-find
rejects a cycle. A disconnected graph has none; a graph of at most one node has one, with
no edge.
"""
function _spanning_trees(n::Int, pairs::Vector{Tuple{Int, Int}})
    n <= 1 && return [Int[]]
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
function _toward(t::Vector{Int}, g::Int, pairs::Vector{Tuple{Int, Int}},
                 W::Dict{Tuple{Int, Int}, _Poly})
    p, seen, todo = one(_Poly), Set([g]), [g]
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
    _bare_law(m, mask, concentrations) -> (; cols, num, den, sigma, trees)

The steady-state law of the edges of `m` that `mask` keeps (at least one), with one
variable per bare constant: the concentrations `concentrations`, then variable
`length(concentrations) + j` for slot `cols[j][2]` of edge `cols[j][1]`, in edge order.
`num/den` is the net consumption of `substrates(m)[1]` per enzyme. Each rapid-equilibrium
segment (`_segments`) is a node, rooted at its first form, so free E (form 1) roots its
own; a form weighs its root's weight, 1, times the equilibrium ratio of each edge on its
path from the root. The forms join the tree in `_segments`' order, each by the first kept
rapid-equilibrium edge in edge (content) order from a form placed before it.
`den = Σ_g sigma[g]·trees[g]`: `sigma[g]` sums segment `g`'s weights and `trees[g]` its
spanning trees of the segment graph oriented toward it, each tree edge weighing the
steady-state flux out of its source. The numerator sums each steady-state edge's net flux
times the copies of the first substrate it moves, its forms' weights included.

A rapid-equilibrium edge off the tree must close a cycle with no net uptake of any
metabolite. Class 2 of the guard rejects every mechanism where one does not, so this is
an internal error.
"""
function _bare_law(m::Mechanism, mask::BitVector, concentrations::Vector{Symbol})
    vars = _Vars()
    foreach(x -> _var!(vars, x; conc = true), concentrations)
    cols = Tuple{Int, Int}[]
    at = Dict{Tuple{Int, Int}, Int}()
    for (k, e) in enumerate(m.edges), slot in 1:(e.rapid_equilibrium ? 1 : 2)
        mask[k] || continue
        push!(cols, (k, slot))
        at[(k, slot)] = _var!(vars, Symbol("c", k, "_", slot))
    end
    conc(c::Counts) = _monomial((vars.index[m.table[i].name] => copies(c, i)
                                 for i in eachindex(m.table))...)
    bare(k, slot) = _monomial(at[(k, slot)] => 1)
    # The equilibrium ratio w(to)/w(from) of rapid-equilibrium edge k: [M]/K for a binding
    # (K a dissociation constant), K·Π[consumed]/Π[released] for any other edge.
    ratio(k) = _mmul(_mpow(bare(k, 1), _is_binding(m.edges[k]) ? -1 : 1),
                     _mmul(conc(m.edges[k].consumed), _mpow(conc(m.edges[k].released), -1)))
    s = _segments(m, mask, _flags(m))
    re = [k for k in eachindex(m.edges) if mask[k] && m.edges[k].rapid_equilibrium]
    alpha = fill(_Mono(), length(m.forms))
    placed = falses(length(m.forms))
    for segment in s.segments
        placed[segment[1]] = true
        for v in segment[2:end]
            k = first(k for k in re if (m.edges[k].to == v && placed[m.edges[k].from]) ||
                                       (m.edges[k].from == v && placed[m.edges[k].to]))
            e = m.edges[k]
            alpha[v] = e.to == v ? _mmul(alpha[e.from], ratio(k)) :
                                   _mmul(alpha[e.to], _mpow(ratio(k), -1))
            placed[v] = true
        end
    end
    for k in re
        e = m.edges[k]
        any(q -> vars.conc[q.first],
            _mmul(_mmul(alpha[e.from], ratio(k)), _mpow(alpha[e.to], -1))) && error(
            "internal error: a rapid-equilibrium cycle of this mechanism turns the " *
            "reaction over (class 2)")
    end
    pairs, W = Tuple{Int, Int}[], Dict{Tuple{Int, Int}, _Poly}()
    flux = Tuple{Int, Int, Int, _Mono, _Mono}[]
    x = vars.index[first(substrates(m))]
    for (k, e) in enumerate(m.edges)
        mask[k] && !e.rapid_equilibrium || continue
        f = _mmul(_mmul(bare(k, 1), conc(e.consumed)), alpha[e.from])
        r = _mmul(_mmul(bare(k, 2), conc(e.released)), alpha[e.to])
        a, b = s.segment_of[e.from], s.segment_of[e.to]
        push!(flux, (Int(_mexp(_mmul(f, _mpow(r, -1)), x)), a, b, f, r))
        a == b && continue
        minmax(a, b) in pairs || push!(pairs, minmax(a, b))
        _addto!(get!(_Poly, W, (a, b)), f, 1)
        _addto!(get!(_Poly, W, (b, a)), r, 1)
    end
    G = length(s.segments)
    spanning = _spanning_trees(G, pairs)
    isempty(spanning) && error("internal error: the steady-state steps of this mechanism " *
                               "do not join its rapid-equilibrium segments")
    trees = [sum((_toward(t, g, pairs, W) for t in spanning); init = zero(_Poly))
             for g in 1:G]
    sigma = [zero(_Poly) for _ in 1:G]
    for (g, segment) in enumerate(s.segments), f in segment
        _addto!(sigma[g], alpha[f], 1)
    end
    num = zero(_Poly)
    for (w, a, b, f, r) in flux
        w == 0 || (num += _scale(trees[a] * f - trees[b] * r, w))
    end
    (; cols, num, den = sum(sigma .* trees; init = zero(_Poly)), sigma, trees)
end

"""
    _column_values(vars, s, wrap = identity) -> Dict{Int, _Mono}

Each column of the solve `s` as a monomial over `vars`: an independent's variable, or a
dependent's power product of `Keq` and independents passed through `wrap`, in
`s.dependents` order (so a `wrap` that defines locals numbers them by column).
"""
function _column_values(vars::_Vars, s::_Solve, wrap = identity)
    value = Dict{Int, _Mono}(c => _monomial(vars.index[s.names[c]] => 1)
                             for c in s.independents)
    for d in s.dependents
        r = _monomial(vars.index[:Keq] => d.keq,
                      (vars.index[s.names[f]] => e for (f, e) in d.factors)...)
        value[d.col] = wrap(r)
    end
    value
end

"""
    _law(m, s, mask; inactive = false, concentrations = metabolites(m), shown = m) -> _Law

The rate law of the edges of `m` that `mask` keeps, in the inactive conformation's
constants when `inactive`, over the independents of the solve `s`: the bare law
(`_bare_law`) with each bare constant replaced by its column's value (an independent, or a
dependent's power product of independents and `Keq`), the numerator divided exactly by `H`,
then both brought to lowest terms in the concentrations. Every law derived from one solve
numbers its variables alike: `concentrations`, then `s.independents`, then `Keq`.

Runs `_assert_derivable` first. Checked `Int64` arithmetic throughout: under the term cap
no coefficient leaves it, so an `OverflowError` becomes an internal error that names the
mechanism `shown`, the allosteric mechanism for an MWC conformation. Class 5
(`_assert_support`) is left to the callers that use the law.
"""
function _law(m::Mechanism, s::_Solve, mask::BitVector; inactive::Bool = false,
              concentrations::Vector{Symbol} = metabolites(m), shown::_AnyMech = m)
    _assert_derivable(m, mask)
    try
        b = _bare_law(m, mask, concentrations)
        vars = _Vars()
        foreach(x -> _var!(vars, x; conc = true), concentrations)
        foreach(i -> _var!(vars, s.names[i]), s.independents)
        keq = _var!(vars, :Keq)
        value = _column_values(vars, s)
        nc = length(concentrations)
        sub(p) = _subst(v -> v <= nc ? nothing :
                             value[s.column[_Column(b.cols[v - nc]..., inactive)]], p)
        side(xs) = _monomial((vars.index[x] => 1 for x in xs)...)
        H = _term(side(substrates(m))) -
            _term(_mmul(side(products(m)), _monomial(keq => -1)))
        Q, D = _divexact(sub(b.num), H), sub(b.den)
        shift = _conc_shift((Q, D), vars)
        _Law(vars, H, Q * shift, D * shift, map(sub, b.sigma), map(sub, b.trees), shift)
    catch err
        err isa OverflowError || rethrow()
        error("internal error: the rate law of this mechanism leaves Int64 arithmetic, " *
              "which the term cap rules out:\n" * sprint(show, shown))
    end
end

"""The law of the edges of `m` over the independents of `s`, every edge kept."""
_law(m::Mechanism, s::_Solve) = _law(m, s, trues(length(m.edges)))

# ─── The MWC law ─────────────────────────────────────────────────────────

"""
An MWC rate law over `vars`: `active.vars` (the concentrations of `metabolites(am)`, the
independents of both conformations and `:Keq`), then the regulator constants (in
`_Solve.fitted` order) and `:L`. With `QA, DA` from `active` and `QI, DI` from `inactive`
(0 and 1 for a conformation that is free E alone, `_flip_only_law`), and
`WA = Π_s sites[s][1]^sites[s][3]`, `WI = Π_s sites[s][2]^sites[s][3]`:

    v = E_total·H·(fI^n·QA·DA^(n−1)·WA + L·fA^n·QI·DI^(n−1)·WI) /
                  (fI^n·DA^n·WA + L·fA^n·DI^n·WI)

`fA` and `fI` are free E's weights in the two laws over their common factor
(`_cross_weights`); `n` is the catalytic multiplicity.

The law is per active site: `E_total` is the active-site concentration, so `n` enters
only as an exponent (of `DA`, `DI` and the free-E weights), never as a leading factor, and
each regulatory site's factor enters numerator and denominator at the site's own
multiplicity, whatever `n` is.
"""
mutable struct _AlloLaw
    const vars::_Vars
    const H::_Poly
    const active::_Law
    const inactive::_Law
    const fA::_Poly
    const fI::_Poly
    const sites::Vector{Tuple{_Poly, _Poly, Int}}
    const n::Int
end
Base.:(==)(a::_AlloLaw, b::_AlloLaw) =
    a.vars == b.vars && a.H == b.H && a.active == b.active && a.inactive == b.inactive &&
    a.fA == b.fA && a.fI == b.fI && a.sites == b.sites && a.n == b.n
Base.hash(A::_AlloLaw, h::UInt) =
    hash((A.vars, A.active, A.inactive, A.n), hash(:_AlloLaw, h))

"""The law of a conformation that is free E alone, over the variables of `L`, the other
conformation's law: it turns nothing over (`Q = 0`), weighs 1 (`D = 1`) and has no
segment."""
_flip_only_law(L::_Law) =
    _Law(L.vars, L.H, zero(_Poly), one(_Poly), _Poly[], [one(_Poly)], _Mono())

"""Whether `L` is the law of a conformation that is free E alone (`_flip_only_law`)."""
_flip_only(L::_Law) = isempty(L.sigma)

"""
    _cross_weights(fA, fI) -> (fA′, fI′)

The free-E weights `fA` and `fI` of the two conformations over their common factor: each
divided by the greatest common divisor of their integer contents and by the monomial of
the lower exponent of each variable in their monomial contents (`_content`). When their
primitive parts are equal those cancel too, leaving monomials, 1 for equal weights.
"""
function _cross_weights(fA::_Poly, fI::_Poly)
    (cA, uA, rA), (cI, uI, rI) = _content(fA), _content(fI)
    same = rA == rI
    g = gcd(cA, cI)
    s = _lowest_shift((_term(uA), _term(uI)))
    weight(c, mu, r) = (p = _term(_mmul(mu, s), div(c, g)); same ? p : p * r)
    weight(cA, uA, rA), weight(cI, uI, rI)
end

"""
The binding polynomial of regulatory site `site` of `am` in the active (`inactive =
false`) or inactive conformation, over `vars`: 1 plus each ligand's concentration over its
dissociation constant (`_reg_K`), an `:OnlyI` ligand left out of the active conformation
and an `:OnlyA` one out of the inactive.
"""
function _site_weight(am::AllostericMechanism, vars::_Vars, site::RegulatorySite,
                      inactive::Bool)
    W = one(_Poly)
    for (lig, tag) in zip(ligands(site), allo_states(site))
        K = _reg_K(am, lig, tag, inactive)
        K === nothing && continue
        W = W + _term(_monomial(vars.index[lig] => 1, vars.index[K] => -1))
    end
    W
end

"""
    _law(am, s = _solve(am)) -> _AlloLaw

The MWC law of `am` over the independents of `s`: the active conformation's law on
`_active_mask(am)` and the inactive one's on `_inactive_mask(am)`, `_flip_only_law` for a
conformation whose mask is empty, both over the concentrations of `metabolites(am)`, so
their polynomials multiply directly; free E's weight in each over their common factor; and
each regulatory site's binding polynomials. Free E roots segment 1, so its weight is that
segment's tree sum times the law's shift. A law over the summed cap is refused before any
polynomial product (`_assert_derivable`). The guard rejects an empty active mask
(`_NO_ACTIVE_TURNOVER`), but the oracles read its law.
"""
function _law(am::AllostericMechanism, s::_Solve)
    _assert_derivable(am)
    m, amask, imask = catalytic(am), _active_mask(am), _inactive_mask(am)
    concentrations = metabolites(am)
    LA = any(amask) ? _law(m, s, amask; concentrations, shown = am) :
         _flip_only_law(_law(m, s, imask; inactive = true, concentrations, shown = am))
    LI = any(imask) ? _law(m, s, imask; inactive = true, concentrations, shown = am) :
         _flip_only_law(LA)
    vars = deepcopy(LA.vars)
    for x in view(s.fitted, (length(s.independents) + 1):(length(s.fitted) - 1))
        _var!(vars, x)
    end
    _var!(vars, :L)
    fA, fI = _cross_weights(LA.trees[1] * LA.shift, LI.trees[1] * LI.shift)
    sites = Tuple{_Poly, _Poly, Int}[(_site_weight(am, vars, site, false),
                                      _site_weight(am, vars, site, true),
                                      multiplicity(site))
                                     for site in regulatory_sites(am)]
    _AlloLaw(vars, LA.H, LA, LI, fA, fI, sites, catalytic_multiplicity(am))
end
_law(m::_AnyMech) = _law(m, _solve(m))

"""
    _with_law(f, m, s = _solve(m)) -> f(m, s, L)

`f` called on `m`, its solve `s` and its law `L = _law(m, s)` once class 5 passes on that
law (`_assert_support`). One method per mechanism type, so `L` is concrete in each and `f`
receives a matching pair.
"""
function _with_law(f, m::Mechanism, s::_Solve = _solve(m))
    L = _law(m, s)
    _assert_support(m, L)
    f(m, s, L)
end
function _with_law(f, am::AllostericMechanism, s::_Solve = _solve(am))
    A = _law(am, s)
    _assert_support(am, A)
    f(am, s, A)
end

"""
The spanning trees of the multigraph whose Laplacian is `lap`: the determinant of `lap`
without its first row and column (Kirchhoff's matrix-tree theorem), by fraction-free
(Bareiss) elimination in checked `Int`. After an `OverflowError` it takes the `BigInt`
determinant, saturated at `typemax(Int)`, which is above every cap. A graph of at most one
node has one tree; a disconnected graph has none. A zero pivot means no tree: the reduced
Laplacian is positive semidefinite, so a zero pivot is a zero leading principal minor and
the determinant is 0; no row swap is needed.
"""
function _tree_count(lap::Matrix{Int})
    n = size(lap, 1) - 1
    n <= 0 && return 1
    M = lap[2:end, 2:end]
    prev = 1
    try
        for k in 1:(n - 1)
            M[k, k] == 0 && return 0
            for i in (k + 1):n, j in (k + 1):n
                M[i, j] = div(Base.checked_sub(Base.checked_mul(M[i, j], M[k, k]),
                                               Base.checked_mul(M[i, k], M[k, j])), prev)
            end
            prev = M[k, k]
        end
        M[n, n]
    catch err
        err isa OverflowError || rethrow()
        Int(min(det(BigInt.(lap[2:end, 2:end])), big(typemax(Int))))
    end
end

"""
Estimate the term count of the law of the edges of `m` that `mask` keeps, V×τ, from its
segment graph, without deriving the law: `V` is the number of rapid-equilibrium segments
(`_segments`), `τ` the spanning trees of the segment graph whose edges are the kept
steady-state edges between segments, parallel ones counted apart (`_tree_count`). V×τ is
the number of spanning-tree products in the law's denominator before kinetic grouping, so
it bounds fit cost. It guards the derivation (`_assert_derivable` aborts when V×τ exceeds
`MAX_RATE_EQUATION_TERMS`, before any polynomial is built) and backs the
`eq_complexity_filter` keyword of [`identify_rate_equation`](@ref). Saturates at
`typemax(Int)`.
"""
function _eq_complexity(m::Mechanism, mask::BitVector)
    s = _segments(m, mask, _flags(m))
    G = length(s.segments)
    G <= 1 && return G
    lap = zeros(Int, G, G)
    for (k, e) in enumerate(m.edges)
        mask[k] && !e.rapid_equilibrium || continue
        g1, g2 = s.segment_of[e.from], s.segment_of[e.to]
        g1 == g2 && continue
        lap[g1, g1] += 1; lap[g2, g2] += 1; lap[g1, g2] -= 1; lap[g2, g1] -= 1
    end
    τ = _tree_count(lap)
    τ > typemax(Int) ÷ G ? typemax(Int) : G * τ
end

"""
    _eq_complexity(m) -> Int

The size of `m`'s rate law in spanning-tree products, V×τ (`_eq_complexity(m, mask)`):
over every edge of a plain mechanism; for an allosteric one, the active law's count plus
the inactive law's, a conformation that is free E alone counting 0, saturated at
`typemax(Int)`. The cap (`_assert_derivable`) and identify's `eq_complexity_filter` read
it.
"""
_eq_complexity(m::Mechanism) = _eq_complexity(m, trues(length(m.edges)))
function _eq_complexity(am::AllostericMechanism)
    m, imask = catalytic(am), _inactive_mask(am)
    a = _eq_complexity(m, _active_mask(am))
    i = any(imask) ? _eq_complexity(m, imask) : 0
    a > typemax(Int) - i ? typemax(Int) : a + i
end

"""
The term cap: the most spanning-tree products (V×τ, `_eq_complexity`) a derived rate law
may have, both conformations of an MWC law together. A larger law takes too long to
compile and is unlikely to be useful for fitting.
"""
const MAX_RATE_EQUATION_TERMS = 5000

"""
The error `_assert_derivable` raises for a mechanism whose rate equation would have more
than `MAX_RATE_EQUATION_TERMS` terms.
"""
struct RateEquationTooLargeError <: Exception end

Base.showerror(io::IO, ::RateEquationTooLargeError) = print(io,
    "Rate equation for this mechanism has more than $MAX_RATE_EQUATION_TERMS " *
    "polynomial terms (limit: $MAX_RATE_EQUATION_TERMS). Equations this large " *
    "take a very long time to compile and are unlikely to be practically " *
    "useful for parameter fitting.")

"""
    _assert_derivable(m, mask)
    _assert_derivable(am)

Abort a derivation whose term count (V×τ, `_eq_complexity` of the same arguments) exceeds
`MAX_RATE_EQUATION_TERMS`: the law of the edges of `m` that `mask` keeps, or both
conformations of `am`'s MWC law together. A quick count on the segment graph that throws
a `RateEquationTooLargeError` before any polynomial is built. Every derivation entry calls
it (`_law`).
"""
function _assert_derivable(args...)
    _eq_complexity(args...) > MAX_RATE_EQUATION_TERMS && throw(RateEquationTooLargeError())
    return nothing
end
