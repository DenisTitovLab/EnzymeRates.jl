# ABOUTME: Code generation: the generated rate_equation emits the derived law as
# ABOUTME: straight-line locals in greedy Horner form, reading its inputs by name.

"""
A rate law as straight-line code, `v = E_total·num/den`. `vars` lists the law's inputs
(concentrations and parameters, read by name), then its locals; `defs[i] = v => p`
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
conformation's constants when `inactive`, with locals added to `defs`: each bare constant
reads its column's value (`_column_values`), a dependent's from a local unless it is one
symbol, the locals defined by column so a body's numbering never follows a hash; each
segment's weight sum `σ_g` and tree sum `T_g` and `Q` get a local (`_factor!`). Returns `Q`
as one term and `D = shift·Σ_g σ_g·T_g`.
"""
function _state!(vars::_Vars, defs::Vector{Pair{Int, _Poly}}, m::Mechanisms.Mechanism,
                 s::_Solve, L::_Law, mask::BitVector, inactive::Bool)
    value = _column_values(vars, s, r -> length(r) <= 1 ? r :
                                         _monomial(_local!(vars, defs, _term(r)) => 1))
    b = _bare_law(m, mask, L.vars.names[L.vars.conc])
    nc = count(L.vars.conc)
    bare(v) = v <= nc ? nothing :
              value[s.column[_Column((b.cols[v - nc]..., inactive))]]
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
    vars, defs = deepcopy(L.vars), Pair{Int, _Poly}[]
    Q, D = _state!(vars, defs, m, s, L, trues(length(Mechanisms.steps(m))), false)
    _Program(vars, defs, _factor!(vars, defs, L.H) * Q, D)
end

"""
The program of an MWC law: each conformation's terms (`_state!`; 0 and 1 for an inactive
conformation that is free E alone), locals for the free-E weights `fA`, `fI`, the
conformation terms `xA = fI·DA` and `xI = fA·DI` and each site's two binding polynomials;
`num = H·(fI·QA·xA^(n−1)·WA + L·fA·QI·xI^(n−1)·WI)`, `den = xA^n·WA + L·xI^n·WI`.
"""
function _program(am::Mechanisms.AllostericMechanism, s::_Solve, A::_AlloLaw)
    m, n = Mechanisms.catalytic(am), A.n
    vars, defs = deepcopy(A.vars), Pair{Int, _Poly}[]
    QA, DA = _state!(vars, defs, m, s, A.active, Mechanisms.conformation_steps(am, :active),
                     false)
    QI, DI = _flip_only(A.inactive) ? (zero(_Poly), one(_Poly)) :
             _state!(vars, defs, m, s, A.inactive,
                     Mechanisms.conformation_steps(am, :inactive), true)
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

The body of `rate_equation` for the program `p`: each variable it reads, in variable
order, bound to a positional local, an input read by name (`getfield(concs, :X)`,
`getfield(params, :k)`) or a local in greedy Horner form, followed by the reciprocals and
roots the body reads (`_power`); then `E_total·num/den`. `x^e` reads `x` to the integer
part of `e` and `x^(±1/q)` to the rest (`x^(3/2) = x·sqrt(x)`).
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
        x = get(def, v) do
            :(getfield($(p.vars.conc[v] ? :concs : :params), $(QuoteNode(p.vars.names[v]))))
        end
        push!(lines, :($(_loc(v, 1 // 1)) = $x))
        for d in sort!([d for (w, d) in used if w == v && d != 1])
            push!(lines, :($(_loc(v, d)) = $(_power(_loc(v, 1 // 1), d))))
        end
    end
    Expr(:block, lines..., :(return getfield(params, :E_total) * $num / $den))
end

"""The body of the generated `rate_equation` of `m`: the program's code (`_body`) of its
law (`_with_law`)."""
_rate_body(m::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}) =
    _with_law((m, s, L) -> _body(_program(m, s, L)), m)

# ─── rate_equation ─────────────────────────────────────────────────────────────

"""
    rate_equation(m, concs, params)

Return the net reaction rate of mechanism `m` at the metabolite concentrations
`concs` and parameters `params` (both `NamedTuple`s). The rate equation itself
is derived from the mechanism's elementary steps by the King–Altman/Cha method:
rapid-equilibrium steps collapse to binding constants, steady-state steps are
assembled into the King–Altman determinant, and the two combine into the
quasi-steady-state flux through the whole mechanism — so one call returns the
overall turnover, not the rate of a single step.

`params` holds the independent constants (`parameters(m)`: the fitted constants, `Keq`
and `E_total`); the rate law computes each dependent constant from them by its
Haldane or Wegscheider relation. Both are read by name, so their fields may come in any
order, and fields the law does not read are ignored.

The derivation runs once, at compile time: the body is generated as straight-line
arithmetic with no allocations, loops, or matrix operations. Use
[`rate_equation_string`](@ref) to inspect the derived equation symbolically.
"""
@generated function rate_equation(m::M_type, concs::NamedTuple,
        params::NamedTuple) where {M_type <: Mechanisms.AbstractEnzymeMechanism}
    # The generator decodes the mechanism and checks its rate law once
    # (`_derivation_input`).
    _rate_body(_derivation_input(M_type))
end

# Concrete-mechanism convenience: lift to the singleton (allocates, not the hot path).
rate_equation(m::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}, concs,
              params) =
    rate_equation(Mechanisms.compile_mechanism(m), concs, params)
