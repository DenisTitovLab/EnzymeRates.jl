# ABOUTME: A mechanism's rate law: sparse polynomials, the law on edge masks and the MWC
# ABOUTME: law, its code (compile_rate_equation), its text and the unit rule it prints by.

# ─── Sparse polynomials ─────────────────────────────────────────────────

"""
A variable table: variable `i` is `names[i]`, a concentration when `conc[i]`; `index` maps
each name to its variable. Two tables are equal when they name the same variables in the
same order.
"""
struct _Vars
    names::Vector{Symbol}
    conc::BitVector
    index::Dict{Symbol, Int}
end
_Vars() = _Vars(Symbol[], falses(0), Dict{Symbol, Int}())
Base.:(==)(a::_Vars, b::_Vars) = a.names == b.names && a.conc == b.conc
Base.hash(t::_Vars, h::UInt) = hash(t.conc, hash(t.names, hash(:_Vars, h)))
Base.copy(t::_Vars) = _Vars(copy(t.names), copy(t.conc), copy(t.index))

"""The index of variable `s` in `t`, added when new; `conc` marks a concentration. Errors
when `s` is already a variable of the other kind."""
function _var!(t::_Vars, s::Symbol; conc::Bool = false)
    i = get(t.index, s, 0)
    i == 0 || return t.conc[i] == conc ? i :
        error("internal error: $s is a concentration in one use and a constant in another")
    push!(t.names, s)
    push!(t.conc, conc)
    t.index[s] = length(t.names)
end

"""An exponent: an integer for a concentration, a rational for a constant."""
const _Exp = Rational{Int}

"""A monomial: `variable => exponent` pairs sorted by variable, no zero exponent. It is
never mutated once built, since polynomials share their monomials as `Dict` keys."""
const _Mono = Vector{Pair{Int, _Exp}}

"""The product of two monomials (a merge of the sorted pairs)."""
function _mmul(a::_Mono, b::_Mono)
    r = _Mono()
    i, j = 1, 1
    while i <= length(a) || j <= length(b)
        if j > length(b) || (i <= length(a) && a[i].first < b[j].first)
            push!(r, a[i])
            i += 1
        elseif i > length(a) || b[j].first < a[i].first
            push!(r, b[j])
            j += 1
        else
            e = a[i].second + b[j].second
            iszero(e) || push!(r, a[i].first => e)
            i += 1
            j += 1
        end
    end
    r
end

_mpow(a::_Mono, e::Union{Integer, _Exp}) =
    iszero(e) ? _Mono() : _Mono([v => x * e for (v, x) in a])

"""The monomial of the `variable => exponent` pairs `ps`; a variable may repeat."""
function _monomial(ps::Pair{Int}...)
    # A plain loop: a reduction over the tuple would build iterator types from the tuple's
    # runtime length and compile again for each new number of pairs.
    m = _Mono()
    for (v, e) in ps
        iszero(e) || (m = _mmul(m, _Mono([v => _Exp(e)])))
    end
    m
end

"""The exponent of variable `v` in `a`."""
function _mexp(a::_Mono, v::Int)
    i = findfirst(p -> p.first == v, a)
    i === nothing ? zero(_Exp) : a[i].second
end

"""
A polynomial: monomial => nonzero coefficient. Every sum and product is checked, so a
coefficient that leaves `Int64` throws `OverflowError` instead of wrapping. It is mutable,
with a `const` field, so a loop that merges it with another type into a `Union` holds a
reference, never an inline box the compiler may hoist. Two polynomials are equal when
their terms are.
"""
mutable struct _Poly
    const terms::Dict{_Mono, Int}
end
_Poly() = _Poly(Dict{_Mono, Int}())

"""Add `c·m` into `p`, dropping a zero sum; returns `p`."""
function _addto!(p::_Poly, m::_Mono, c::Integer)
    v = Base.checked_add(get(p.terms, m, 0), Int(c))
    iszero(v) ? delete!(p.terms, m) : (p.terms[m] = v)
    p
end

"""The one-term polynomial `c·m`."""
_term(m::_Mono, c::Integer = 1) = _addto!(_Poly(), m, c)

Base.one(::Type{_Poly}) = _term(_Mono())
Base.zero(::Type{_Poly}) = _Poly()
Base.iszero(p::_Poly) = isempty(p.terms)
Base.:(==)(p::_Poly, q::_Poly) = p.terms == q.terms
Base.hash(p::_Poly, h::UInt) = hash(p.terms, hash(:_Poly, h))

_nterms(p::_Poly) = length(p.terms)

function Base.:+(p::_Poly, q::_Poly)
    r = _Poly(copy(p.terms))
    for (m, c) in q.terms
        _addto!(r, m, c)
    end
    r
end
Base.:-(p::_Poly) = _scale(p, -1)
Base.:-(p::_Poly, q::_Poly) = p + (-q)

function Base.:*(p::_Poly, q::_Poly)
    r = _Poly()
    for (m1, c1) in p.terms, (m2, c2) in q.terms
        _addto!(r, _mmul(m1, m2), Base.checked_mul(c1, c2))
    end
    r
end
# Multiplying by a monomial maps distinct monomials to distinct ones, so no terms merge.
Base.:*(p::_Poly, m::_Mono) =
    _Poly(Dict{_Mono, Int}(_mmul(k, m) => c for (k, c) in p.terms))

_scale(p::_Poly, c::Integer) = iszero(c) ? _Poly() :
    _Poly(Dict{_Mono, Int}(m => Base.checked_mul(x, Int(c)) for (m, x) in p.terms))

"""
`p` with each variable `v` for which `f(v)` returns a monomial (a power product) replaced
by it, `v^e` by `f(v)^e`; `f(v) === nothing` keeps `v`.
"""
function _subst(f, p::_Poly)
    r = _Poly()
    for (m, c) in p.terms
        out = _Mono()
        for (v, e) in m
            s = f(v)
            out = _mmul(out, s === nothing ? _Mono([v => e]) : _mpow(s, e))
        end
        _addto!(r, out, c)
    end
    r
end

"""
    _divexact(p, d) -> p / d

The exact quotient of `p` by a two-term `d` (the Haldane factor); an `ArgumentError` when
`d` does not divide `p`. With `d = c_u·u + c_w·w` and `t = w/u`, the terms of `p/u` fall
into chains `key·t^k` over integer `k`; dividing a chain by `c_u + c_w·t` is the recurrence
`b_k = (a_k − c_w·b_{k−1}) / c_u`, which must end at zero.
"""
function _divexact(p::_Poly, d::_Poly)
    _nterms(d) == 2 || throw(ArgumentError("_divexact divides by a two-term polynomial"))
    (u, cu), (w, cw) = collect(d.terms)
    t = _mmul(w, _mpow(u, -1))
    v, et = first(t)
    chains = Dict{_Mono, Vector{Pair{Int, Int}}}()
    for (m, c) in p.terms
        m = _mmul(m, _mpow(u, -1))
        k = floor(Int, _mexp(m, v) / et)
        push!(get!(chains, _mmul(m, _mpow(t, -k)), Pair{Int, Int}[]), k => c)
    end
    q = _Poly()
    for (key, chain) in chains
        sort!(chain; by = first)
        b, j = 0, 1
        for k in chain[1].first:chain[end].first
            a = 0
            if chain[j].first == k
                a = chain[j].second
                j += 1
            end
            r = Base.checked_sub(a, Base.checked_mul(cw, b))
            rem(r, cu) == 0 ||
                throw(ArgumentError("_divexact: the quotient is not integral"))
            b = div(r, cu)
            k == chain[end].first || iszero(b) || _addto!(q, _mmul(key, _mpow(t, k)), b)
        end
        iszero(b) || throw(ArgumentError("_divexact: the divisor does not divide"))
    end
    q
end

"""
    _lowest_shift(ps, counted = every variable) -> _Mono

The monomial in the `counted` variables that, multiplied into every polynomial of `ps`,
brings each such variable's lowest exponent over all their terms to zero, a term without
the variable counting as exponent 0: the shift to lowest terms.
"""
function _lowest_shift(ps, counted = _ -> true)
    low, seen = Dict{Int, _Exp}(), Dict{Int, Int}()
    n = sum(p -> length(p.terms), ps; init = 0)
    for p in ps, m in keys(p.terms)
        for (v, e) in m
            counted(v) || continue
            low[v] = min(get(low, v, e), e)
            seen[v] = get(seen, v, 0) + 1
        end
    end
    _monomial((v => -(seen[v] < n ? min(e, 0) : e) for (v, e) in low)...)
end

"""The shift to lowest terms in the concentrations of `t` alone (`_lowest_shift`)."""
_conc_shift(ps, t::_Vars) = _lowest_shift(ps, v -> t.conc[v])

"""`p` as `concentration monomial => its coefficient`, a polynomial in the constants."""
function _split_conc(p::_Poly, t::_Vars)
    out = Dict{_Mono, _Poly}()
    for (m, c) in p.terms
        _addto!(get!(_Poly, out, filter(q -> t.conc[q.first], m)),
                filter(q -> !t.conc[q.first], m), c)
    end
    out
end

"""
The nonzero polynomial `p` as `c·mu·rest`: its integer content `c` (the coefficients'
gcd, positive, or a one-term polynomial's own coefficient), its monomial content `mu`
(each variable's lowest exponent over the terms, a term without the variable counting as
exponent 0) and the rest, `one(_Poly)` for a one-term polynomial.
"""
function _content(p::_Poly)
    c, s = foldl(gcd, values(p.terms)), _lowest_shift((p,))
    c, _mpow(s, -1), _Poly(Dict{_Mono, Int}(_mmul(mo, s) => div(x, c)
                                            for (mo, x) in p.terms))
end

# ─── One conformation's law ─────────────────────────────────────────────

"""
The rate law of one conformation, `v = H·Q/D` per unit enzyme, over `vars`: the
concentrations, then the solve's independent constants (`_Solve.independents` order), then
`:Keq`.
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
of positions in `pairs`. A disconnected graph has none; a graph of at most one node has
one, with no edge.
"""
function _spanning_trees(n::Int, pairs::Vector{Tuple{Int, Int}})
    n <= 1 && return [Int[]]
    trees = Vector{Int}[]
    _grow_spanning_trees!(trees, n, pairs, 1, Int[], collect(1:n))
    trees
end

"""Add to `trees` each spanning tree that extends `chosen` with edges of `pairs` from
position `i` on, `root` being the union-find parents under `chosen`."""
function _grow_spanning_trees!(trees::Vector{Vector{Int}}, n::Int,
                               pairs::Vector{Tuple{Int, Int}}, i::Int,
                               chosen::Vector{Int}, root::Vector{Int})
    if length(chosen) == n - 1
        push!(trees, copy(chosen))
        return
    end
    length(pairs) - i + 1 < n - 1 - length(chosen) && return
    a, b = _union_find_root(root, pairs[i][1]), _union_find_root(root, pairs[i][2])
    if a != b
        r = copy(root)
        r[a] = b
        push!(chosen, i)
        _grow_spanning_trees!(trees, n, pairs, i + 1, chosen, r)
        pop!(chosen)
    end
    _grow_spanning_trees!(trees, n, pairs, i + 1, chosen, root)
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
`num/den` is the net consumption of `m.substrates[1]` per enzyme. Each rapid-equilibrium
segment (`Mechanisms.rapid_equilibrium_segments`) is a node, rooted at its first form, so
free E (form 1) roots its own; a form weighs its root's weight, 1, times the equilibrium
ratio of each edge on its path from the root. The forms join the tree in the order of
`Mechanisms.rapid_equilibrium_segments`, each by the first kept
rapid-equilibrium edge in step order from a form placed before it.
`den = Σ_g sigma[g]·trees[g]`: `sigma[g]` sums segment `g`'s weights and `trees[g]` its
spanning trees of the segment graph oriented toward it, each tree edge weighing the
steady-state flux out of its source. The numerator sums each steady-state edge's net flux
times the copies of the first substrate it moves, its forms' weights included.

A rapid-equilibrium edge off the tree must close a cycle with no net uptake of any
metabolite. Class 2 of the guard rejects every mechanism where one does not, so this is
an internal error.
"""
function _bare_law(m::Mechanisms.CatalyticScheme, mask::BitVector,
                   concentrations::Vector{Symbol})
    vars = _Vars()
    foreach(x -> _var!(vars, x; conc = true), concentrations)
    cols = Tuple{Int, Int}[]
    at = Dict{Tuple{Int, Int}, Int}()
    for k in eachindex(m.steps), slot in 1:(m.steps[k].rapid_equilibrium ? 1 : 2)
        mask[k] || continue
        push!(cols, (k, slot))
        at[(k, slot)] = _var!(vars, Symbol("c", k, "_", slot))
    end
    conc(held::Vector{Int}) =
        _monomial((vars.index[m.metabolites[i]] => 1 for i in held)...)
    bare(k, slot) = _monomial(at[(k, slot)] => 1)
    # The equilibrium ratio w(to)/w(from) of rapid-equilibrium edge k:
    # K·Π[taken up]/Π[given off] for a constant K named after the stored direction, the
    # same with 1/K for one named after the reverse ([M]/K for a binding stored as its
    # uptake, K a dissociation constant).
    ratio(k) = _mmul(_mpow(bare(k, 1), _constant_named_in_reverse(m.steps[k]) ? -1 : 1),
                     _mmul(conc(m.steps[k].taken_up),
                           _mpow(conc(m.steps[k].given_off), -1)))
    s = Mechanisms.rapid_equilibrium_segments(m, mask)
    re = [k for k in eachindex(m.steps)
          if mask[k] && m.steps[k].rapid_equilibrium]
    nforms = length(m.enzyme_forms)
    alpha = fill(_Mono(), nforms)
    placed = falses(nforms)
    # Whether edge k joins form v to a placed form.
    joins_placed(k, v) = (a = m.steps[k].from; b = m.steps[k].to;
                          (b == v && placed[a]) || (a == v && placed[b]))
    for segment in s.segments
        placed[segment[1]] = true
        for v in segment[2:end]
            k = first(k for k in re if joins_placed(k, v))
            a, b = m.steps[k].from, m.steps[k].to
            alpha[v] = b == v ? _mmul(alpha[a], ratio(k)) :
                                _mmul(alpha[b], _mpow(ratio(k), -1))
            placed[v] = true
        end
    end
    for k in re
        a, b = m.steps[k].from, m.steps[k].to
        any(q -> vars.conc[q.first],
            _mmul(_mmul(alpha[a], ratio(k)), _mpow(alpha[b], -1))) && error(
            "internal error: a rapid-equilibrium cycle of this mechanism turns the " *
            "reaction over (class 2)")
    end
    pairs, W = Tuple{Int, Int}[], Dict{Tuple{Int, Int}, _Poly}()
    flux = Tuple{Int, Int, Int, _Mono, _Mono}[]
    x = vars.index[m.metabolites[first(m.substrates)]]
    for k in eachindex(m.steps)
        mask[k] && !m.steps[k].rapid_equilibrium || continue
        from, to = m.steps[k].from, m.steps[k].to
        f = _mmul(_mmul(bare(k, 1), conc(m.steps[k].taken_up)), alpha[from])
        r = _mmul(_mmul(bare(k, 2), conc(m.steps[k].given_off)), alpha[to])
        a, b = s.segment_of[from], s.segment_of[to]
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
    _law(am, s, mask; inactive = false) -> _Law

The rate law of the edges of `am`'s catalytic scheme that `mask` keeps, over the
concentrations of `metabolites(am)`, in the inactive conformation's constants when
`inactive`, over the independents of the solve `s`: the bare law (`_bare_law`) with each
bare constant replaced by its column's value (an independent, or a dependent's power
product of independents and `Keq`), the numerator divided exactly by `H`, then both
brought to lowest terms in the concentrations. Every law derived from one solve numbers
its variables alike: the concentrations, then `s.independents`, then `Keq`.

Checked `Int64` arithmetic throughout: under the term cap, which `_law(am, s)` checks, no
coefficient leaves it, so an `OverflowError` becomes an internal error that names `am`.
"""
function _law(am::Mechanisms.Mechanism, s::_Solve, mask::BitVector; inactive::Bool = false)
    m, concentrations = am.catalytic, Mechanisms.metabolites(am)
    try
        b = _bare_law(m, mask, concentrations)
        vars = _Vars()
        foreach(x -> _var!(vars, x; conc = true), concentrations)
        foreach(i -> _var!(vars, s.names[i]), s.independents)
        keq = _var!(vars, :Keq)
        value = _column_values(vars, s)
        nc = length(concentrations)
        sub(p) = _subst(v -> v <= nc ? nothing :
                             value[s.column[_StepParameter(b.cols[v - nc]..., inactive)]],
                      p)
        side(xs) = _monomial((vars.index[m.metabolites[x]] => 1 for x in xs)...)
        H = _term(side(m.substrates)) -
            _term(_mmul(side(m.products), _monomial(keq => -1)))
        Q, D = _divexact(sub(b.num), H), sub(b.den)
        shift = _conc_shift((Q, D), vars)
        _Law(vars, H, Q * shift, D * shift, map(sub, b.sigma), map(sub, b.trees), shift)
    catch err
        err isa OverflowError || rethrow()
        error("internal error: the rate law of this mechanism leaves Int64 arithmetic, " *
              "which the term cap rules out:\n" * sprint(show, am))
    end
end

# ─── The MWC law ────────────────────────────────────────────────────────

"""
An MWC rate law over `vars`: `active.vars` (the concentrations of `metabolites(am)`, the
independents of both conformations and `:Keq`), then the regulator constants (in
`_Solve.fitted` order) and `:L`. With `QA, DA` from `active` and `QI, DI` from `inactive`
(0 and 1 for an inactive conformation that is free E alone, `_flip_only_law`), and
`WA = Π_s sites[s][1]^sites[s][3]`, `WI = Π_s sites[s][2]^sites[s][3]`:

    v = H·(fI^n·QA·DA^(n−1)·WA + L·fA^n·QI·DI^(n−1)·WI) /
          (fI^n·DA^n·WA + L·fA^n·DI^n·WI)

`fA` and `fI` are free E's weights in the two laws over their common factor
(`_cross_weights`); `n` is the oligomeric state.

The law is per active site: `v` is the rate per catalytic site, so `n` enters
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
dissociation constant (`_regulator_constant`), an `:OnlyI` ligand left out of the active
conformation and an `:OnlyA` one out of the inactive.
"""
function _site_weight(am::Mechanisms.Mechanism, vars::_Vars,
                      site::Mechanisms.RegulatorySite, inactive::Bool)
    W = one(_Poly)
    for (lig, tag) in zip(site.ligands, site.allosteric_states)
        K = _regulator_constant(am, lig, tag, inactive)
        K === nothing && continue
        W = W + _term(_monomial(vars.index[lig] => 1, vars.index[K] => -1))
    end
    W
end

"""
    _law(m, s) -> Union{_Law, _AlloLaw}

The law of mechanism `m` over the independents of its solve `s`: a plain mechanism's law
on every step, a `_Law`; an allosteric mechanism's MWC law, an `_AlloLaw`: each
conformation's law on its steps (`steps_running_in_allosteric_conformations`),
`_flip_only_law` for an inactive conformation that runs none, both over the concentrations
of `metabolites(m)`, so their polynomials multiply directly; free E's weight in each over
their common factor; and each regulatory site's binding polynomials. Free E roots segment 1,
so its weight is that segment's tree sum times the law's shift. A law over the term cap
(`rate_equation_term_count`, `MAX_RATE_EQUATION_TERMS`) throws a `RateEquationTooLargeError`
before any polynomial is built. `m` passes the guard, which rejects an active conformation
that runs no step (`_NO_ACTIVE_TURNOVER`).
"""
function _law(am::Mechanisms.Mechanism, s::_Solve)
    rate_equation_term_count(am) > MAX_RATE_EQUATION_TERMS &&
        throw(RateEquationTooLargeError())
    _is_allosteric(am) || return _law(am, s, trues(length(am.catalytic.steps)))
    runs = Mechanisms.steps_running_in_allosteric_conformations(am)
    LA = _law(am, s, runs.active)
    LI = any(runs.inactive) ? _law(am, s, runs.inactive; inactive = true) :
         _flip_only_law(LA)
    vars = copy(LA.vars)
    for x in view(s.fitted, (length(s.independents) + 1):(length(s.fitted) - 1))
        _var!(vars, x)
    end
    _var!(vars, :L)
    fA, fI = _cross_weights(LA.trees[1] * LA.shift, LI.trees[1] * LI.shift)
    sites = Tuple{_Poly, _Poly, Int}[(_site_weight(am, vars, site, false),
                                      _site_weight(am, vars, site, true),
                                      site.multiplicity)
                                     for site in am.allosteric.regulatory_sites]
    _AlloLaw(vars, LA.H, LA, LI, fA, fI, sites, am.allosteric.oligomeric_state)
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
    _term_count(m, mask) -> Int

The term count V×τ (`rate_equation_term_count`) of the law of the edges of catalytic
scheme `m` that `mask` keeps, read off its segment graph
(`Mechanisms.rapid_equilibrium_segments`, `_tree_count`). Saturates at `typemax(Int)`.
"""
function _term_count(m::Mechanisms.CatalyticScheme, mask::BitVector)
    s = Mechanisms.rapid_equilibrium_segments(m, mask)
    G = length(s.segments)
    G <= 1 && return G
    lap = zeros(Int, G, G)
    for k in s.steps
        a, b = m.steps[k].from, m.steps[k].to
        g1, g2 = s.segment_of[a], s.segment_of[b]
        lap[g1, g1] += 1; lap[g2, g2] += 1; lap[g1, g2] -= 1; lap[g2, g1] -= 1
    end
    τ = _tree_count(lap)
    τ > typemax(Int) ÷ G ? typemax(Int) : G * τ
end

"""
    rate_equation_term_count(m) -> Int

The size of `m`'s rate law in spanning-tree products, V×τ, read off its graph without
deriving the law. `V` is the number of rapid-equilibrium segments, and `τ` the number of
spanning trees of the segment graph, whose edges are the steady-state steps between
segments, parallel ones counted apart. V×τ is the number of spanning-tree products in the
law's denominator before kinetic grouping, so it bounds fit cost. A plain mechanism's
count reads every step; an allosteric one's is the active law's count plus the inactive
law's, a conformation that is free E alone counting 0. Saturates at `typemax(Int)`. The
term cap, `MAX_RATE_EQUATION_TERMS`, bounds this count.
"""
function rate_equation_term_count(am::Mechanisms.Mechanism)
    m = am.catalytic
    (; active, inactive) = Mechanisms.steps_running_in_allosteric_conformations(am)
    a = _term_count(m, active)
    _is_allosteric(am) || return a
    i = any(inactive) ? _term_count(m, inactive) : 0
    a > typemax(Int) - i ? typemax(Int) : a + i
end

"""
The term cap: the most spanning-tree products (V×τ, `rate_equation_term_count`) a derived
rate law may have, both conformations of an MWC law together. A larger law takes too long to
compile and is unlikely to be useful for fitting.
"""
const MAX_RATE_EQUATION_TERMS = 5000

"""
The error `_law` raises for a mechanism whose rate equation would have more than
`MAX_RATE_EQUATION_TERMS` terms.
"""
struct RateEquationTooLargeError <: Exception end

Base.showerror(io::IO, ::RateEquationTooLargeError) = print(io,
    "Rate equation for this mechanism has more than $MAX_RATE_EQUATION_TERMS " *
    "polynomial terms (limit: $MAX_RATE_EQUATION_TERMS). Equations this large " *
    "take a very long time to compile and are unlikely to be practically " *
    "useful for parameter fitting.")

# ─── A law as straight-line code ────────────────────────────────────────

"""
A rate law as straight-line code, `v = num/den`, the rate per unit enzyme. `vars` lists
the law's inputs (concentrations, constants and `Keq`), then its locals; `defs[i] = v => p`
defines local `v` as the polynomial `p` over the inputs and earlier locals. Every term of
`num` and `den` holds each concentration to a nonnegative power.
"""
struct _Program
    vars::_Vars
    defs::Vector{Pair{Int, _Poly}}
    num::_Poly
    den::_Poly
end

"""The local of `vars` that `defs` defines as `q`, defined anew unless one already is."""
function _local!(vars::_Vars, defs::Vector{Pair{Int, _Poly}}, q::_Poly)
    i = findfirst(d -> d.second == q, defs)
    i === nothing || return defs[i].first
    push!(defs, _var!(vars, Symbol("#", length(defs) + 1)) => q)
    defs[end].first
end

"""`p` as one term: `p` itself when it has at most one, else the monomial of its terms'
lowest concentration exponents times a local holding the rest, which is then in lowest
terms and so finite at zero concentrations."""
function _factor!(vars::_Vars, defs::Vector{Pair{Int, _Poly}}, p::_Poly)
    _nterms(p) <= 1 && return p
    s = _conc_shift((p,), vars)
    _term(_mmul(_monomial(_local!(vars, defs, p * s) => 1), _mpow(s, -1)))
end

"""
    _state!(vars, defs, m, s, L, mask, inactive) -> (Q, D)

The terms over `vars` of the law `L` of the edges `mask` of `m`, in the inactive
conformation's constants when `inactive`, with locals added to `defs`. Each bare constant
reads its column's value (`_column_values`), a dependent's from a local unless it is one
symbol; those locals are defined by column, so a body's numbering never follows a hash.
Each segment's weight sum `σ_g` and tree sum `T_g`, and `Q`, get a local (`_factor!`).
Returns `Q` as one term and `D = shift·Σ_g σ_g·T_g`.
"""
function _state!(vars::_Vars, defs::Vector{Pair{Int, _Poly}}, m::Mechanisms.CatalyticScheme,
                 s::_Solve, L::_Law, mask::BitVector, inactive::Bool)
    value = _column_values(vars, s, r -> length(r) <= 1 ? r :
                                         _monomial(_local!(vars, defs, _term(r)) => 1))
    b = _bare_law(m, mask, L.vars.names[L.vars.conc])
    nc = count(L.vars.conc)
    bare(v) = v <= nc ? nothing :
              value[s.column[_StepParameter(b.cols[v - nc]..., inactive)]]
    shift = _term(L.shift)
    D = zero(_Poly)
    for (σ, T) in zip(b.sigma, b.trees)
        D = D + shift * _factor!(vars, defs, _subst(bare, σ)) *
                _factor!(vars, defs, _subst(bare, T))
    end
    _factor!(vars, defs, L.Q), D
end

"""The program of a plain law: `num = H·Q`, `den = D` (`_state!`)."""
function _program(m::Mechanisms.Mechanism, s::_Solve, L::_Law)
    vars, defs = copy(L.vars), Pair{Int, _Poly}[]
    Q, D = _state!(vars, defs, m.catalytic, s, L, trues(length(m.catalytic.steps)), false)
    _Program(vars, defs, _factor!(vars, defs, L.H) * Q, D)
end

"""
The program of an MWC law: each conformation's terms (`_state!`; 0 and 1 for an inactive
conformation that is free E alone), locals for the free-E weights `fA`, `fI`, the
conformation terms `xA = fI·DA` and `xI = fA·DI` and each site's two binding polynomials;
`num = H·(fI·QA·xA^(n−1)·WA + L·fA·QI·xI^(n−1)·WI)`, `den = xA^n·WA + L·xI^n·WI`.
"""
function _program(am::Mechanisms.Mechanism, s::_Solve, A::_AlloLaw)
    m, n = am.catalytic, A.n
    vars, defs = copy(A.vars), Pair{Int, _Poly}[]
    runs = Mechanisms.steps_running_in_allosteric_conformations(am)
    QA, DA = _state!(vars, defs, m, s, A.active, runs.active, false)
    QI, DI = _flip_only(A.inactive) ? (zero(_Poly), one(_Poly)) :
             _state!(vars, defs, m, s, A.inactive, runs.inactive, true)
    f(p) = _factor!(vars, defs, p)
    pw(p, k) = foldl(*, (p for _ in 1:k); init = one(_Poly))
    fA, fI = f(A.fA), f(A.fI)
    xA, xI = f(fI * DA), f(fA * DI)
    WA = foldl(*, (pw(f(w), k) for (w, _, k) in A.sites); init = one(_Poly))
    WI = foldl(*, (pw(f(w), k) for (_, w, k) in A.sites); init = one(_Poly))
    L = _term(_monomial(vars.index[:L] => 1))
    num = f(A.H) * (fI * QA * pw(xA, n - 1) * WA + L * fA * QI * pw(xI, n - 1) * WI)
    _Program(vars, defs, num, pw(xA, n) * WA + L * pw(xI, n) * WI)
end

"""The local holding variable `v` to the unit power `d` (`1`, `-1`, `1//q` or `-1//q`):
`x<v>`, its reciprocal `x<v>m1`, its root `x<v>p<q>` or that root's reciprocal
`x<v>m<q>`."""
_loc(v::Int, d::_Exp) =
    d == 1 ? Symbol(:x, v) : Symbol(:x, v, d < 0 ? :m : :p, denominator(d))

"""`x` to the unit power `d`: a reciprocal (`inv`), a square root, or for `q > 2` the root
`x^(one(x)/q)`, computed in `x`'s own precision."""
function _power(x::Symbol, d::_Exp)
    q = denominator(d)
    base = q == 1 ? x : q == 2 ? :(sqrt($x)) : :($x^(one($x) / $q))
    d < 0 ? :(inv($base)) : base
end

"""Greedy Horner form of the terms `ts`, each `(coefficient, local => power)`: take out
the local most terms hold (ties by name), then recurse on the terms with it and on the
rest, with binary `+` and `*` only."""
function _horner(ts::Vector{Tuple{Int, Dict{Symbol, Int}}})
    c0 = sum((c for (c, a) in ts if isempty(a)); init = 0)
    rest = Tuple{Int, Dict{Symbol, Int}}[(c, a) for (c, a) in ts if !isempty(a)]
    isempty(rest) && return c0
    held = Dict{Symbol, Int}()
    for (_, a) in rest, x in keys(a)
        held[x] = get(held, x, 0) + 1
    end
    best = first(sort!(collect(keys(held)); by = x -> (-held[x], String(x))))
    lower(a) = (b = copy(a); b[best] == 1 ? delete!(b, best) : (b[best] -= 1); b)
    inner = _horner(Tuple{Int, Dict{Symbol, Int}}[(c, lower(a))
                                                  for (c, a) in rest if haskey(a, best)])
    ex = inner == 1 ? best : :($best * $inner)
    others = Tuple{Int, Dict{Symbol, Int}}[(c, a) for (c, a) in rest if !haskey(a, best)]
    iszero(c0) || push!(others, (c0, Dict{Symbol, Int}()))
    isempty(others) ? ex : :($ex + $(_horner(others)))
end

"""`p` as code, `atoms(mono)` giving each monomial's locals and their powers: the Horner
form of its positive terms minus that of its negative ones."""
function _emit(p::_Poly, atoms)
    ts = Tuple{Int, Dict{Symbol, Int}}[(c, atoms(mo))
                                       for (mo, c) in sort!(collect(p.terms); by = first)]
    pos = Tuple{Int, Dict{Symbol, Int}}[t for t in ts if t[1] > 0]
    neg = Tuple{Int, Dict{Symbol, Int}}[(-c, a) for (c, a) in ts if c < 0]
    isempty(neg) ? _horner(pos) : isempty(pos) ? :(-$(_horner(neg))) :
        :($(_horner(pos)) - $(_horner(neg)))
end

"""
    _body(p) -> Expr

The body of a compiled rate law's call for the program `p`: each variable it reads, in
variable order, bound to a positional local, an input (a concentration or a constant read
by name, `getfield(concs, :X)`, `getfield(params, :k)`, or the argument `Keq`) or a local
in greedy Horner form, followed by the reciprocals and roots the body reads (`_power`);
then `num/den`. `x^e` reads `x` to the integer part of `e` and `x^(±1/q)` to the rest
(`x^(3/2) = x·sqrt(x)`).
"""
function _body(p::_Program)
    used = Set{Tuple{Int, _Exp}}()
    function atoms(mo::_Mono)
        a = Dict{Symbol, Int}()
        for (v, e) in mo
            n = trunc(Int, e)
            r = e - n
            for (d, k) in ((sign(n) // 1, abs(n)),
                           (sign(r) // denominator(r), numerator(abs(r))))
                k == 0 && continue
                push!(used, (v, d))
                a[_loc(v, d)] = k
            end
        end
        a
    end
    def = Dict{Int, Any}(v => _emit(q, atoms) for (v, q) in p.defs)
    num, den = _emit(p.num, atoms), _emit(p.den, atoms)
    lines = Expr[]
    for v in sort!(unique!([first(u) for u in used]))
        name = p.vars.names[v]
        x = get(def, v) do
            name === :Keq ? :Keq :
                :(getfield($(p.vars.conc[v] ? :concs : :params), $(QuoteNode(name))))
        end
        push!(lines, :($(_loc(v, 1 // 1)) = $x))
        for d in sort!([d for (w, d) in used if w == v && d != 1])
            push!(lines, :($(_loc(v, d)) = $(_power(_loc(v, 1 // 1), d))))
        end
    end
    Expr(:block, lines..., :(return $num / $den))
end

"""The body of the generated call of `m`'s compiled rate law: the program's code
(`_body`) of its law (`_law`)."""
function _rate_body(m::Mechanisms.Mechanism)
    s = _solve(m)
    _body(_program(m, s, _law(m, s)))
end

# ─── compile_rate_equation ──────────────────────────────────────────────
#
# A mechanism's Sig, the type parameter of its compiled rate law, is the mechanism as nested
# tuples of Symbols, Ints and Bools, its fields in declaration order:
#   (metabolites, enzyme_forms, substrates, products,
#    ((from, to, taken_up, given_off, rapid_equilibrium, kinetic_group), ...),
#    allosteric_states, oligomeric_state, ((ligands, multiplicity, allosteric_states), ...))
# each vector a tuple of its entries.

"""`m` as nested tuples, its Sig (layout above): the only type parameter of its compiled
rate law. Equal mechanisms have equal Sigs, since a mechanism is stored in canonical
form."""
function _sig_of(m::Mechanisms.Mechanism)
    c, a = m.catalytic, m.allosteric
    # Each tuple is a vector splatted into the builtin `tuple`, and the vectors of steps
    # and sites have their element type given, so the encoder dispatches nothing and
    # compiles once, whatever the Sig's type.
    steps = Tuple{Int, Int, Tuple{Vararg{Int}}, Tuple{Vararg{Int}}, Bool, Int}[
        (e.from, e.to, (e.taken_up...,), (e.given_off...,), e.rapid_equilibrium,
         e.kinetic_group) for e in c.steps]
    sites = Tuple{Tuple{Vararg{Symbol}}, Int, Tuple{Vararg{Symbol}}}[
        ((s.ligands...,), s.multiplicity, (s.allosteric_states...,))
        for s in a.regulatory_sites]
    ((c.metabolites...,), (c.enzyme_forms...,), (c.substrates...,), (c.products...,),
     (steps...,), (a.allosteric_states...,), a.oligomeric_state, (sites...,))
end

"""The mechanism whose Sig is `sig`, built by the constructors, which check it. It takes
the Sig unspecialized, so it compiles once, whatever the mechanism."""
function _mechanism_from_sig(@nospecialize(sig::Tuple))
    metabolites, forms, substrates, products, steps, states, oligomeric_state, sites = sig
    catalytic = Mechanisms.CatalyticScheme(
        Symbol[metabolites...], Symbol[forms...], Int[substrates...], Int[products...],
        Mechanisms.Step[Mechanisms.Step(s[1], s[2], Int[s[3]...], Int[s[4]...], s[5], s[6])
                        for s in steps])
    regulatory_sites = Mechanisms.RegulatorySite[
        Mechanisms.RegulatorySite(Symbol[s[1]...], s[2], Symbol[s[3]...]) for s in sites]
    Mechanisms.Mechanism(catalytic,
                         Mechanisms.Allostery(Symbol[states...], oligomeric_state,
                                              regulatory_sites))
end

"""The rate law of the mechanism whose Sig (`_sig_of`) is `Sig`, as a function
(`compile_rate_equation`). Its only parameter is the Sig, so equal mechanisms share one
instance and its call is generated once per mechanism."""
struct _CompiledRateEquation{Sig} end

"""The mechanism whose compiled rate law has the type `R`, decoded from its Sig
(`_mechanism_from_sig`) and checked as `check_valid_rate_equation` checks it, its guard
verdict taken afresh rather than through the memo (`_Memo`): what the law's generated
methods derive from. It reads the Sig off `R` at run time, so it and the decoder compile
once, whatever the mechanism."""
function _derivation_input(@nospecialize(R::Type{<:_CompiledRateEquation}))
    m = _mechanism_from_sig(R.parameters[1])
    _check_verdict(m, _guard(m))
    m
end

"""
    compile_rate_equation(m) -> rate

The rate law of mechanism `m` as a function: `rate(concs, params, Keq)` returns the net
reaction rate per unit enzyme at the metabolite concentrations `concs` and the constants
`params` (both `NamedTuple`s), with the equilibrium constant `Keq`. Multiply by the enzyme
concentration for an absolute rate. The law is derived from the mechanism's elementary
steps by the King–Altman/Cha method, so one call gives the turnover of the whole
mechanism, not the rate of one step.

`params` holds the constants a fit estimates, [`parameters`](@ref)`(rate)`, the same names
as `parameters(m)`, and `concs` the concentrations of
[`metabolites`](@ref EnzymeRates.Mechanisms.metabolites)`(rate)`; both are tuples of
Symbols. The law computes each dependent constant from `params` and `Keq` by
its Haldane or Wegscheider relation. Both `NamedTuple`s are read by name, so their fields
may come in any order, and fields the law does not read are ignored.

The first call derives the law and compiles it to straight-line arithmetic, once per
mechanism and session; every later call is fast and allocates nothing, so compile once,
outside a loop. Equal mechanisms give the same `rate`. Calling `rate` errors when the rate
law of `m` is not defined (`check_valid_rate_equation`), and with a
`RateEquationTooLargeError` when its term count (`rate_equation_term_count`) exceeds the
cap of 5,000 (`MAX_RATE_EQUATION_TERMS`). [`rate_equation_string`](@ref)
prints the law.

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S ⇌ ES
               ES <--> EP
               EP ⇌ E + P
           end
       end;

julia> rate = compile_rate_equation(m)
compiled rate equation

julia> parameters(rate)
(:K_EP_to_E_P, :K_ES_to_E_S, :k_ES_to_EP)

julia> params = (K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0);

julia> rate((S = 1.0, P = 0.5), params, 4.0)
0.65625
```
"""
compile_rate_equation(m::Mechanisms.Mechanism) = _CompiledRateEquation{_sig_of(m)}()

# In each generator, the argument `compiled_law` is the law's type. Its name holds `_`, as
# no metabolite name does, and the body reads every input by name, so a metabolite named
# like it reads as a concentration.
@generated function (compiled_law::_CompiledRateEquation)(concs::NamedTuple,
                                                         params::NamedTuple, Keq::Real)
    # The generator decodes the mechanism and checks its rate law once
    # (`_derivation_input`).
    _rate_body(_derivation_input(compiled_law))
end

# `Fitting.loss!` builds its `NamedTuple`s from these names, so they must be compile-time
# constants: each is generated from the Sig, once per mechanism. The metabolites need no
# check of the law.
@generated parameters(compiled_law::_CompiledRateEquation) =
    Tuple(_solve(_derivation_input(compiled_law)).fitted)
@generated Mechanisms.metabolites(compiled_law::_CompiledRateEquation) =
    Tuple(Mechanisms.metabolites(_mechanism_from_sig(compiled_law.parameters[1])))

# A compiled law prints one line whatever its mechanism; `@nospecializeinfer` keeps a caller
# compiled for one law, such as `sprint`, from inferring a show for that law.
Base.@nospecializeinfer Base.show(io::IO, @nospecialize(rate::_CompiledRateEquation)) =
    print(io, "compiled rate equation")

# ─── Text of powers, terms and sums ─────────────────────────────────────

"""`x^e` as text for a positive exponent `e`: `x`, `x^2`, `sqrt(x)`, `x^(3/2)`."""
_power_text(x, e::_Exp) = e == 1 ? string(x) : e == 1 // 2 ? "sqrt($x)" :
    isinteger(e) ? "$x^$(numerator(e))" : "$x^($(numerator(e))/$(denominator(e)))"

"""Factors as a ratio, `a * b / (c * d)`; `1` when no factor is above the bar."""
_ratio_text(up::Vector{String}, dn::Vector{String}) =
    string(isempty(up) ? "1" : join(up, " * "),
           isempty(dn) ? "" :
           " / " * (length(dn) == 1 ? only(dn) : "(" * join(dn, " * ") * ")"))

"""
    _term_text(c, mo, names; rest = one(_Poly), over = false) -> String

The positive term `c·mo` as a ratio of products over the variable names `names`, times the
polynomial `rest` and over the normalizer's printed factor `N` when `over`:
`2 * A * k1 * (k2 + k3) / (k4 * N)`, and `1` for the empty term.
"""
function _term_text(c, mo::_Mono, names::Vector{Symbol}; rest::_Poly = one(_Poly),
                    over::Bool = false)
    up = String[_power_text(names[v], e) for (v, e) in mo if e > 0]
    dn = String[_power_text(names[v], -e) for (v, e) in mo if e < 0]
    numerator(c) == 1 || pushfirst!(up, string(numerator(c)))
    denominator(c) == 1 || pushfirst!(dn, string(denominator(c)))
    rest == one(_Poly) || push!(up, "(" * _sum_text(rest, names) * ")")
    over && push!(dn, "N")
    _ratio_text(up, dn)
end

"""The nonzero polynomial `p` over the variable names `names` as a sum, `a + 2 * b - c`,
its terms in monomial order."""
function _sum_text(p::_Poly, names::Vector{Symbol})
    s = join(string(c < 0 ? " - " : " + ", _term_text(abs(c), mo, names))
             for (mo, c) in sort!(collect(p.terms); by = first))
    (startswith(s, " + ") ? "" : "-") * chop(s; head = 3, tail = 0)
end

const _SUPERSCRIPTS = Dict(zip("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))

"""
A coefficient's name for the concentration monomial `mo` over the names `names`: `0` for
none, else each concentration in variable order with its exponent as a superscript, joined
by `_`: `A²_B`.
"""
function _conc_label(mo::_Mono, names::Vector{Symbol})
    isempty(mo) && return "0"
    sup(e) = e == 1 ? "" : String([_SUPERSCRIPTS[d] for d in string(numerator(e))])
    join((string(names[v], sup(e)) for (v, e) in mo), "_")
end

# ─── A law's polynomials over its normalizer ────────────────────────────

"""
One printed term of a law's polynomial: the concentration monomial `conc` times its
coefficient `c·mo·rest`, over the normalizer's printed factor `N` when `over`. `c` and `mo`
have the normalizer's integer and monomial contents divided out; `rest` is the
coefficient's polynomial factor, `one(_Poly)` when there is none.
"""
struct _Term
    conc::_Mono
    c::Rational{Int}
    mo::_Mono
    rest::_Poly
    over::Bool
end

"""Multiply the factors printed apart, `apart` (polynomial and exponent), by `p^e`: an
equal polynomial adds its exponent, and a factor whose exponent sums to 0 goes."""
function _apart!(apart::Vector{Tuple{_Poly, _Exp}}, p::_Poly, e::_Exp)
    i = findfirst(x -> x[1] == p, apart)
    i === nothing ? push!(apart, (p, e)) : (apart[i] = (p, apart[i][2] + e))
    filter!(x -> !iszero(x[2]), apart)
end

"""
    _printed(L, nz, polys, subs) -> (; N, terms)

The polynomials `polys` of the conformation law `L` divided by its normalizer `nz`
(`_normalizer`, `subs` the substrates' variables). The normalizer's integer contents to
integer powers and monomial contents divide every coefficient; the other factors, an
integer content to a fractional power included, print apart: `N` is their text, empty
when there is none. `terms[j]` holds a `_Term` for each concentration monomial of
`polys[j]`, in the unit rule's order (`_term_key`). A coefficient whose polynomial factor
is `N`'s one factor drops it and is not over `N`.
"""
function _printed(L::_Law, nz::Vector{Pair{_Mono, _Exp}}, polys::Vector{_Poly},
                  subs::Vector{Int})
    D, names = _split_conc(L.D, L.vars), L.vars.names
    k, mu, apart = 1 // 1, _Mono(), Tuple{_Poly, _Exp}[]
    for (mo, e) in nz
        c, u, rest = _content(D[mo])
        if isinteger(e)
            k *= (c // 1)^numerator(e)
        elseif c != 1
            _apart!(apart, _term(_Mono(), c), e)
        end
        mu = _mmul(mu, _mpow(u, e))
        rest == one(_Poly) || _apart!(apart, rest, e)
    end
    paren(p) = "(" * _sum_text(p, names) * ")"
    N = isempty(apart) ? "" :
        length(apart) == 1 && apart[1][2] == 1 ? _sum_text(apart[1][1], names) :
        _ratio_text(String[_power_text(paren(p), e) for (p, e) in apart if e > 0],
                    String[_power_text(paren(p), -e) for (p, e) in apart if e < 0])
    whole = length(apart) == 1 && apart[1][2] == 1 ? apart[1][1] : zero(_Poly)
    # `term` captures only names assigned once, so none is boxed (the loop reassigns `k`).
    scale, inv_mu = k, _mpow(mu, -1)
    function term(conc, q)
        c, u, rest = _content(q)
        w = rest == whole
        _Term(conc, (c // 1) / scale, _mmul(u, inv_mu), w ? one(_Poly) : rest,
              !isempty(apart) && !w)
    end
    terms = Vector{_Term}[
        _Term[term(conc, q) for (conc, q) in sort!(collect(_split_conc(p, L.vars));
                                                   by = x -> _term_key(x[1], L.vars, subs))]
        for p in polys]
    (; N, terms)
end

"""The terms `ts` (`_printed`) as a sum, each concentration monomial times its
coefficient."""
_group_sum(ts::Vector{_Term}, names::Vector{Symbol}) =
    join((_term_text(t.c, _mmul(t.conc, t.mo), names; rest = t.rest, over = t.over)
          for t in ts), " + ")

"""
    _coefficient_sum(prefix, ts, names) -> (sum, lines)

The terms `ts` (`_printed`) as a sum of named coefficients, `prefix_A_B * A * B`, with
`prefix_0` alone for the constant monomial and `1` when its coefficient is 1, and one
line per named coefficient giving its value in the fitted constants.
"""
function _coefficient_sum(prefix::String, ts::Vector{_Term}, names::Vector{Symbol})
    unit(t) = isempty(t.conc) && t.c == 1 && isempty(t.mo) && t.rest == one(_Poly) &&
              !t.over
    label(t) = string(prefix, "_", _conc_label(t.conc, names))
    sum = join((unit(t) ? "1" : isempty(t.conc) ? label(t) :
                string(label(t), " * ", _term_text(1, t.conc, names)) for t in ts),
               " + ")
    lines = String[string(label(t), " = ",
                          _term_text(t.c, t.mo, names; rest = t.rest, over = t.over))
                   for t in ts if !unit(t)]
    sum, lines
end

"""The `v` line, `v = (Q) * (H) / (D)`: `Q` bare when it is one name."""
function _v_line(q::AbstractString, h::AbstractString, d::AbstractString)
    lead = occursin(r"^[^ ()*/+^-]+$", q) ? q * " * " : "(" * q * ") * "
    "v = $(lead)($h) / ($d)"
end

const _COEFFICIENTS_COMMENT =
    "# The coefficients are functions of the fitted constants, not parameters:"

# ─── Reference lines and header ─────────────────────────────────────────

"""
The reference lines of the solve `s` over the variables `vars`: a `# Wegscheider
constraints:` section, then a `# Haldane constraints:` section, each dependent constant by
name as a ratio of products of `Keq` and the independents with exact exponents (`sqrt`,
`^(3/2)`), `= 1` for a constant that thermodynamics pins to 1. A section without a line is
left out; a kinetic group's other members equal its kept constant and print no line.
"""
function _reference_lines(s::_Solve, vars::_Vars)
    lines, value = String[], _column_values(vars, s)
    for (label, title) in ((:wegscheider, "Wegscheider"), (:haldane, "Haldane"))
        ds = sort!(filter(d -> d.label === label, s.dependents); by = d -> s.names[d.col])
        isempty(ds) || push!(lines, "# $title constraints:")
        for d in ds
            push!(lines, "$(s.names[d.col]) = $(_term_text(1, value[d.col], vars.names))")
        end
    end
    lines
end

"""The `concs` destructuring line of `m`'s printed law."""
_concs_line(m::Mechanisms.Mechanism) =
    "(; $(join(Mechanisms.metabolites(m), ", "))) = concs"

"""The `params` and `concs` destructuring lines of `m`'s printed law over the solve `s`."""
_header_lines(m::Mechanisms.Mechanism, s::_Solve) =
    String["(; $(join(s.fitted, ", "))) = params", _concs_line(m)]

# ─── Plain and MWC laws ─────────────────────────────────────────────────

"""
    _law_text(m, s, L; coefficients = false, coefficient_values = true) -> String

The printed law `L` of `m` over the solve `s`, laid out line by line as
`rate_equation_string` describes, in the order the lines evaluate.
"""
function _law_text(m::Mechanisms.Mechanism, s::_Solve, L::_Law; coefficients::Bool = false,
                   coefficient_values::Bool = true)
    subs = _substrate_vars(m, L.vars)
    P = _printed(L, _normalizer(L.D, L.vars, subs), [L.Q, L.D], subs)
    names, (q, d) = L.vars.names, P.terms
    h = _sum_text(L.H, names)
    norms = isempty(P.N) ? String[] : ["N = $(P.N)"]
    if coefficients
        (qs, lq) = _coefficient_sum("n", q, names)
        (ds, ld) = _coefficient_sum("d", d, names)
        coefficient_values || return join([_concs_line(m), _v_line(qs, h, ds)], "\n")
        return join([_header_lines(m, s); _COEFFICIENTS_COMMENT; norms; lq; ld;
                     _v_line(qs, h, ds)], "\n")
    end
    join([_header_lines(m, s); _reference_lines(s, L.vars); norms;
          _v_line(_group_sum(q, names), h, _group_sum(d, names))], "\n")
end

"""
The conformation law `L` of an MWC law, tagged `tag` (`"A"` or `"I"`), printed over its own
normalizer: the terms (`_printed`) of its `Q`, its `D` and its free-E weight (the tree sum
of free E's segment times the law's shift), each over the normalizer's printed factor,
which reads as the name `N_<tag>` appended to the variable names.
"""
function _conformation(L::_Law, tag::String, subs::Vector{Int})
    names = [L.vars.names; Symbol("N_", tag)]
    n_var = _monomial(length(names) => -1)
    over(ts) = _Term[_Term(t.conc, t.c, t.over ? _mmul(t.mo, n_var) : t.mo, t.rest, false)
                     for t in ts]
    P = _printed(L, _normalizer(L.D, L.vars, subs), [L.Q, L.D, L.trees[1] * L.shift], subs)
    (; tag, names, N = P.N, q = over(P.terms[1]), d = over(P.terms[2]),
     f = over(P.terms[3]))
end

"""Push `W_<tag> = (1 + …)^k * …`, the product of the regulatory sites' binding
polynomials `ws` (polynomial and multiplicity) in one conformation, onto `defs`; a site
whose polynomial is 1 is left out, and so is the line when every site is."""
function _site_line!(defs::Vector{String}, tag::String, ws::Vector{Tuple{_Poly, Int}},
                     names::Vector{Symbol})
    parts = [string("(", _sum_text(w, names), ")", k == 1 ? "" : "^$k")
             for (w, k) in ws if _nterms(w) > 1]
    isempty(parts) || push!(defs, "W_$tag = " * join(parts, " * "))
    defs
end

"""
The printed MWC law `A` of `am` over the solve `s`, as `rate_equation_string` describes,
with these definitions:
- each conformation `c`, `A` and then `I` unless it is free E alone, defines `Q_c` (unless
  it turns nothing over) and `D_c` over its own normalizer `N_c`, which prints as a line
  when it has a factor that is not a single term, or an integer to a fractional power;
- `F_A` and `F_I`, the free-E weights over the normalizers, print only when they differ,
  and each only when it is not 1;
- `W_A` and `W_I` multiply the regulatory sites' factors.
The `v` line is `H·(F_I^n·Q_A·D_A^(n−1)·W_A + L·F_A^n·Q_I·D_I^(n−1)·W_I) /
(F_I^n·D_A^n·W_A + L·F_A^n·D_I^n·W_I)` without its factors of 1. Both forms open with the
header and end with `v`, in the order they evaluate. In full, the reference lines and the
definitions come between. With `coefficients`, the coefficient comment, the normalizers,
one line per coefficient, named by conformation (`nA_0`, `dA_S`, `fI_B`), and the
definitions come between: `Q_c`, `D_c` and `F_c` read in those coefficients, `W_A` and
`W_I` in the fitted constants. Without `coefficient_values`, the `concs` line, the
definitions and `v`: the lines that name coefficients without defining them.
"""
function _law_text(am::Mechanisms.Mechanism, s::_Solve, A::_AlloLaw;
                   coefficients::Bool = false, coefficient_values::Bool = true)
    subs = _substrate_vars(am, A.vars)
    confs = [_conformation(A.active, "A", subs)]
    _flip_only(A.inactive) || push!(confs, _conformation(A.inactive, "I", subs))
    defs, norms, coefs = String[], String[], String[]
    function define!(name, ts, c)
        if coefficients
            prefix = string(name[1] == 'Q' ? 'n' : lowercase(name[1]), c.tag)
            text, lines = _coefficient_sum(prefix, ts, c.names)
            push!(defs, "$name = $text")
            append!(coefs, lines)
        else
            push!(defs, "$name = $(_group_sum(ts, c.names))")
        end
    end
    for c in confs
        isempty(c.N) || push!(coefficients ? norms : defs, "N_$(c.tag) = $(c.N)")
        isempty(c.q) || define!("Q_$(c.tag)", c.q, c)
        define!("D_$(c.tag)", c.d, c)
    end
    flips = [_group_sum(c.f, c.names) for c in confs]
    if flips[1] != (length(flips) == 2 ? flips[2] : "1")
        for (c, t) in zip(confs, flips)
            t == "1" || define!("F_$(c.tag)", c.f, c)
        end
    end
    _site_line!(defs, "A", Tuple{_Poly, Int}[(w, k) for (w, _, k) in A.sites], A.vars.names)
    _site_line!(defs, "I", Tuple{_Poly, Int}[(w, k) for (_, w, k) in A.sites], A.vars.names)
    has(x) = any(l -> startswith(l, x * " = "), defs)
    pw(x, e) = !has(x) || e == 0 ? "" : e == 1 ? x : "$x^$e"
    product(fs...) = (f = filter(!isempty, String[fs...]);
                      isempty(f) ? "1" : join(f, " * "))
    branch(c, o, e) = product(pw("F_$o", A.n), pw("Q_$c", e), pw("D_$c", A.n - e),
                              pw("W_$c", 1))
    num = String[]
    has("Q_A") && push!(num, branch("A", "I", 1))
    has("Q_I") && push!(num, "L * " * branch("I", "A", 1))
    rest = branch("I", "A", 0)
    den = branch("A", "I", 0) * " + L" * (rest == "1" ? "" : " * " * rest)
    v = _v_line(join(num, " + "), _sum_text(A.H, A.vars.names), den)
    coefficient_values || return join([_concs_line(am); defs; v], "\n")
    head = _header_lines(am, s)
    coefficients ? join([head; _COEFFICIENTS_COMMENT; norms; coefs; defs; v], "\n") :
                   join([head; _reference_lines(s, A.vars); defs; v], "\n")
end

# ─── rate_equation_string ───────────────────────────────────────────────

"""
    rate_equation_string(m; coefficients = false, coefficient_values = true) -> String

Return the rate law of mechanism `m` as a multi-line `String` (it returns the text; it does
not print).

The text runs as Julia code once `params`, `concs` and `Keq` are bound, as they are in
a call of [`compile_rate_equation`](@ref)`(m)`. It destructures
[`parameters`](@ref)`(m)` from `params` and the metabolites from `concs`, and reads the
equilibrium constant `Keq` as an input of its own. It lists each dependent constant under
`# Wegscheider constraints:` or `# Haldane constraints:` as a reference line in `Keq` and
the independent constants. It ends with `v = (Q) * (H) / (D)`, the rate per unit enzyme,
in the independent constants. `H = Π[S] − Π[P]/Keq` vanishes at equilibrium; `Q` and `D`
have positive terms. Numerator and denominator are divided by a normalizer that makes
every term of `D` dimensionless, when one exists, chosen by the rules of
[Normalization](@ref); a law with none prints as derived, and its terms keep their units.
A factor of the normalizer that is not a single term, or an integer to a fractional power,
prints once as `N = …`, after the reference lines.
An allosteric mechanism prints each conformation's numerator `Q_A`, `Q_I`, denominator
`D_A`, `D_I`, normalizer `N_A`, `N_I`, free-enzyme weight `F_A`, `F_I` (when they differ)
and regulatory-site factor `W_A`, `W_I`, then the MWC combination.

With `coefficients = true` the law reads in named coefficients,
`v = n_0 * (H) / (1 + d_S * S + …)`, each named by the concentrations it
multiplies (`d_A²_B` for `A^2 * B`). That text runs as Julia code too. It destructures
`params` and `concs`, then prints the coefficient comment, `N` (`N_A`, `N_I`) when the
normalizer prints, and each coefficient's value in the fitted constants, one line each.
For an allosteric law it prints `Q_c`, `D_c` and `F_c` in those coefficients and the site
factors `W_c`. It ends with the `v` line; the reference lines print only in the full law.

`coefficient_values = false` (with `coefficients = true`; `ArgumentError` otherwise) leaves
out everything that defines a coefficient: the `params` line, the comment, `N` and each
coefficient's value. What remains is the `concs` line, for an allosteric law the `Q_c`,
`D_c`, `F_c` and `W_c` lines, and the `v` line. This text names coefficients it does not
define, so it does not run; it is the readable form the search CSVs hold.

Errors when the rate law of `m` is not defined (`check_valid_rate_equation`), and with a
`RateEquationTooLargeError` when its term count (`rate_equation_term_count`) exceeds the
cap of 5,000 (`MAX_RATE_EQUATION_TERMS`).

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S ⇌ ES
               ES <--> EP
               EP ⇌ E + P
           end
       end;

julia> print(rate_equation_string(m))
(; K_EP_to_E_P, K_ES_to_E_S, k_ES_to_EP) = params
(; S, P) = concs
# Haldane constraints:
k_EP_to_ES = K_EP_to_E_P * k_ES_to_EP / (K_ES_to_E_S * Keq)
v = (k_ES_to_EP / K_ES_to_E_S) * (S - P / Keq) / (1 + S / K_ES_to_E_S + P / K_EP_to_E_P)

julia> print(rate_equation_string(m; coefficients = true))
(; K_EP_to_E_P, K_ES_to_E_S, k_ES_to_EP) = params
(; S, P) = concs
# The coefficients are functions of the fitted constants, not parameters:
n_0 = k_ES_to_EP / K_ES_to_E_S
d_S = 1 / K_ES_to_E_S
d_P = 1 / K_EP_to_E_P
v = n_0 * (S - P / Keq) / (1 + d_S * S + d_P * P)

julia> print(rate_equation_string(m; coefficients = true, coefficient_values = false))
(; S, P) = concs
v = n_0 * (S - P / Keq) / (1 + d_S * S + d_P * P)
```
"""
function rate_equation_string(m::Mechanisms.Mechanism; coefficients::Bool = false,
                              coefficient_values::Bool = true)
    coefficient_values || coefficients ||
        throw(ArgumentError("coefficient_values = false needs coefficients = true"))
    check_valid_rate_equation(m)
    _with_memo_law((m, s, L) -> _law_text(m, s, L; coefficients, coefficient_values), m)
end

# ─── The unit rule ──────────────────────────────────────────────────────

"""The variables of `vars` that hold `m`'s substrates (the unit rule's `substrates`)."""
_substrate_vars(m::Mechanisms.Mechanism, vars::_Vars) =
    Int[vars.index[m.catalytic.metabolites[x]] for x in m.catalytic.substrates]

"""
    _term_key(mo, vars, substrates) -> key

The unit rule's order on the concentration monomials of a denominator (`_normalizer`):
total degree first, then the number of distinct concentrations that are not substrates
(`substrates` holds the substrates' variables), then the exponent vector over the
concentrations of `vars`, in `vars` order, with an earlier name's larger exponent first,
so `A²` comes before `A·B`.
"""
function _term_key(mo::_Mono, vars::_Vars, substrates::Vector{Int})
    conc = findall(vars.conc)
    ex = _Exp[_mexp(mo, v) for v in conc]
    (sum(ex), count(i -> !iszero(ex[i]) && !(conc[i] in substrates), eachindex(conc)), -ex)
end

"""
    _normalizer(D, vars, substrates) -> Vector{Pair{_Mono, _Exp}}

The normalizer of the denominator `D` over `vars` by the unit rule (the derivation
tutorial's "Normalization" states it for users), as the pairs `m => λ_m` of
`N = Π_m c(m)^λ_m`, `c(m)` being the coefficient of the concentration monomial `m` in `D`;
no pairs for `N = 1`. Every term of `D` has one unit, each concentration counting as a
dimension of its own, and every term divided by `N` is dimensionless unless the last rule
applies. The first rule that applies gives `N`:

1. The constant term: `D` has a term free of concentrations, and `N` is its coefficient.
2. Cleland's terms: with `ΠS` the product of the `k` substrates, `D` has the terms
   `ΠS/X` for every substrate `X` and `ΠS`, and `N = Π_X c(ΠS/X) / c(ΠS)^(k - 1)`,
   Cleland's normalizer: each substrate term then reads as its concentrations over their
   Michaelis constants. One substrate never takes it, since `ΠS/A` is the absent
   constant term.
3. A graded basis: taking the concentration monomials `m` of `D` in `_term_key` order,
   the first ones whose vectors `[1; a_m]` (`a_m` the exponent vector) are affinely
   independent form a basis; when weights on the basis satisfy `Σ λ_m = 1` and
   `Σ λ_m·a_m = 0`, `N = Π c(m)^λ_m`. The weights are exact and may be fractional.
4. None: when no such weights exist, no product of the coefficients has `D`'s unit: some
   weighting of the concentrations gives every term of `D` the same nonzero degree. In
   `E + A <--> EA`, `EA + B ⇌ EQ + P`, `EQ <--> E + Q`, for example, the
   rapid-equilibrium Theorell–Chance step makes every term carry exactly one of `B` and
   `P`. `N = 1`: the law then prints as derived, and its terms keep their units.

Exact `Rational{Int}` arithmetic; an `OverflowError` escapes.
"""
function _normalizer(D::_Poly, vars::_Vars, substrates::Vector{Int})
    Dc = _split_conc(D, vars)
    haskey(Dc, _Mono()) && return Pair{_Mono, _Exp}[_Mono() => 1]
    S = _monomial((v => 1 for v in substrates)...)
    cleland = Pair{_Mono, _Exp}[_mmul(S, _monomial(v => -1)) => 1 for v in substrates]
    push!(cleland, S => 1 - length(substrates))
    all(p -> haskey(Dc, p.first), cleland) && return cleland
    conc = findall(vars.conc)
    monos = sort!(collect(keys(Dc)); by = mo -> _term_key(mo, vars, substrates))
    column(mo) = _Exp[1; [_mexp(mo, v) for v in conc]]
    basis = _Mono[]
    for mo in monos
        candidate = reduce(hcat, [column(b) for b in [basis; [mo]]])
        length(_rref!(candidate, length(basis) + 1)) > length(basis) && push!(basis, mo)
    end
    k = length(basis)
    M = [reduce(hcat, [column(b) for b in basis]) _Exp[1; zeros(Int, length(conc))]]
    _rref!(M, k)
    all(iszero, view(M, (k + 1):size(M, 1), k + 1)) || return Pair{_Mono, _Exp}[]
    Pair{_Mono, _Exp}[b => M[i, k + 1] for (i, b) in enumerate(basis)
                      if !iszero(M[i, k + 1])]
end
