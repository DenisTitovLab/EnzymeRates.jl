# ABOUTME: The unit rule of a rate law: the normalizer that makes every denominator term
# ABOUTME: dimensionless, by the constant term, Cleland's terms, a graded basis or a ratio.

"""
The normalizer of one conformation's denominator (`_normalizer`): `N = Π_m c(m)^λ_m ·
x^shift` over the pairs `m => λ_m` of `lambda`, `c(m)` being the coefficient of the
concentration monomial `m` in the denominator. `rule` names the case that chose it:
`:constant`, the constant term; `:cleland`, Cleland's terms ΠS/X for each substrate X and
ΠS; `:graded`, weights on the first affinely independent terms; `:ratio`, the first term,
whose concentration monomial `shift` stays in `N`. `shift` is empty under the other three
rules.
"""
struct _Norm
    rule::Symbol
    lambda::Vector{Pair{_Mono, _Exp}}
    shift::_Mono
end

"""The variables of `vars` that hold `m`'s substrates (the unit rule's `substrates`)."""
_substrate_vars(m::Union{Mechanism, AllostericMechanism}, vars::_Vars) =
    Int[vars.index[x] for x in substrates(m)]

"""
    _term_key(mo, vars, substrates) -> key

The unit rule's order on the concentration monomials of a denominator (spec section 9):
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
    _normalizer(D, vars, substrates) -> _Norm

The normalizer of the denominator `D` over `vars` (spec section 9, "Normalization: the
unit rule"). Every term of `D` has one unit, and every term divided by `N` is
dimensionless. The first rule that applies gives `N`:

1. `:constant`: `D` has a term free of concentrations, and `N` is its coefficient.
2. `:cleland`: with `ΠS` the product of the `k` distinct substrates, `D` has the terms
   `ΠS/X` for every substrate `X` and `ΠS`, and `N = Π_X c(ΠS/X) / c(ΠS)^(k - 1)`,
   Cleland's normalizer: the substrate terms then read as products of Michaelis
   constants. One substrate never takes it, since `ΠS/A` is the absent constant term.
3. `:graded`: taking the concentration monomials `m` of `D` in `_term_key` order, the
   first ones whose vectors `[1; a_m]` (`a_m` the exponent vector) are affinely
   independent form a basis; when weights on the basis satisfy `Σ λ_m = 1` and
   `Σ λ_m·a_m = 0`, `N = Π c(m)^λ_m`. The weights are exact and may be fractional.
4. `:ratio`: otherwise the law reads some concentrations only through their ratio, and
   `N` is the first term in that order, its concentration monomial included.

Exact `Rational{Int}` arithmetic; an `OverflowError` escapes.
"""
function _normalizer(D::_Poly, vars::_Vars, substrates::Vector{Int})
    Dc = _split_conc(D, vars)
    haskey(Dc, _Mono()) && return _Norm(:constant, Pair{_Mono, _Exp}[_Mono() => 1], _Mono())
    subs = unique(substrates)
    S = _monomial((v => 1 for v in subs)...)
    cleland = Pair{_Mono, _Exp}[_mmul(S, _monomial(v => -1)) => 1 for v in subs]
    push!(cleland, S => 1 - length(subs))
    all(p -> haskey(Dc, p.first), cleland) && return _Norm(:cleland, cleland, _Mono())
    conc = findall(vars.conc)
    monos = sort!(collect(keys(Dc)); by = mo -> _term_key(mo, vars, substrates))
    column(mo) = _Exp[1; [_mexp(mo, v) for v in conc]]
    basis = _Mono[]
    for mo in monos
        M = reduce(hcat, [column(b) for b in [basis; [mo]]])
        length(_rref!(M, length(basis) + 1)) > length(basis) && push!(basis, mo)
    end
    k = length(basis)
    M = [reduce(hcat, [column(b) for b in basis]) _Exp[1; zeros(Int, length(conc))]]
    _rref!(M, k)
    all(iszero, view(M, (k + 1):size(M, 1), k + 1)) ||
        return _Norm(:ratio, Pair{_Mono, _Exp}[first(monos) => 1], first(monos))
    _Norm(:graded, Pair{_Mono, _Exp}[b => M[i, k + 1] for (i, b) in enumerate(basis)
                                     if !iszero(M[i, k + 1])], _Mono())
end
