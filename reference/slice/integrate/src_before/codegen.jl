# ABOUTME: Code generation: rate_equation of a compiled mechanism, generated from its Sig;
# ABOUTME: the law is derived in the generator and emitted as locals by greedy Horner.

"""
A rate law as straight-line code, `v = E_total·num/den`. `vars` lists the law's inputs
(concentrations and parameters, read by name), then the locals `#1`, `#2`, …; `defs[i] =
v => p` defines local `v` as the polynomial `p` over the inputs and earlier locals. Every
term holds each concentration to a nonnegative power.
"""
struct Program{T <: Integer}
    vars::Vars
    defs::Vector{Pair{Int, Poly{T}}}
    num::Poly{T}
    den::Poly{T}
end

"""The local of `vars` that `defs` defines as `q`, defined anew (as `#k`) unless one is."""
function local!(vars::Vars, defs, q)
    i = findfirst(d -> d.second == q, defs)
    i === nothing || return defs[i].first
    push!(defs, var!(vars, Symbol("#", length(defs) + 1)) => q)
    defs[end].first
end

"""`p` as one term: `p` itself when it has at most one, else the monomial of its terms'
lowest concentration exponents times a local holding the rest (in lowest terms, so finite at
zero concentrations)."""
function factor!(vars::Vars, defs, p::Poly{T}) where {T}
    nterms(p) <= 1 && return p
    s = conc_shift((p,), vars)
    term(T, mmul(mono(local!(vars, defs, p * s) => 1), mpow(s, -1)))
end

"""
    state!(vars, defs, m, sol, L, mask, inactive) -> (Q, D)

The terms of the law `L = law(m, sol, mask; inactive)` over `vars`, with locals added to
`defs`: one for each power product of independents and `Keq` that a dependent constant
equals (a dependent equal to one symbol reads it), one for `Q`, and one for each segment's
weight sum `σ_g` and tree sum `T_g` over those constants. Returns `Q` as one term and
`D = shift·Σ_g σ_g·T_g`.
"""
function state!(vars, defs, m::Mechanism, sol::Solution, L::Law{T}, mask,
                inactive) where {T}
    b, names, nc = bare_law(m, mask, T), column_names(m, sol), length(metabolites(m))
    var(i) = vars.index[names[i]]
    deps, memo = Dict(d.col => d for d in sol.deps), Dict{Int, Mono}()
    function constant(c)
        haskey(deps, c) || return mono(var(c) => 1)
        d = deps[c]
        r = mono((var(f) => e for (f, e) in d.terms)..., vars.index[:Keq] => d.keq)
        length(r) <= 1 ? r : mono(local!(vars, defs, term(T, r)) => 1)
    end
    col(v) = get!(() -> constant(sol.index[Col(b.cols[v - nc]..., inactive)]), memo, v)
    sub(p) = subst(p, v -> v <= nc ? nothing : col(v))
    f(p) = factor!(vars, defs, p)
    D = sum(term(T, L.shift) * f(sub(s)) * f(sub(t)) for (s, t) in zip(b.sigma, b.trees))
    f(L.Q), D
end

"""The program of a plain law: `num = H·Q`, `den = D` (`state!`)."""
function program(m::Mechanism, sol::Solution, L::Law{T}) where {T}
    vars, defs = deepcopy(L.vars), Pair{Int, Poly{T}}[]
    Q, D = state!(vars, defs, m, sol, L, trues(length(m.edges)), false)
    Program(vars, defs, factor!(vars, defs, L.H) * Q, D)
end

"""
The program of an allosteric law: each state's terms (`state!`), locals for the cross
weights `fA`, `fI`, the state terms `xA = fI·DA` and `xI = fA·DI` and each site's two bases;
`num = H·(fI·QA·xA^(n-1)·WA + L·fA·QI·xI^(n-1)·WI)`, `den = xA^n·WA + L·xI^n·WI`.
"""
function program(m::Mechanism, sol::Solution, A::AlloLaw{T}) where {T}
    vars, defs, n = deepcopy(A.vars), Pair{Int, Poly{T}}[], A.n
    QA, DA = state!(vars, defs, m, sol, A.active, trues(length(m.edges)), false)
    QI, DI = A.inactive === nothing ? (zero(Poly{T}), one(Poly{T})) :
             state!(vars, defs, m, sol, A.inactive, inactive_mask(m), true)
    f(p) = factor!(vars, defs, p)
    pw(p, k) = prod(fill(p, k); init = one(Poly{T}))
    fA, fI = f(A.fA), f(A.fI)
    xA, xI = f(fI * DA), f(fA * DI)
    WA, WI = (prod((pw(f(s[i]), s[3]) for s in A.sites); init = one(Poly{T})) for i in 1:2)
    L = term(T, mono(vars.index[:L] => 1))
    num = f(A.H) * (fI * QA * pw(xA, n - 1) * WA + L * fA * QI * pw(xI, n - 1) * WI)
    Program(vars, defs, num, pw(xA, n) * WA + L * pw(xI, n) * WI)
end

"""The local holding variable `v` to the unit power `d` (`1`, `-1` or `±1/q`)."""
loc(v, d) = d == 1 ? Symbol(:x, v) : Symbol(:x, v, d < 0 ? :m : :p, denominator(d))

"""`x` to the unit power `d`: a reciprocal, `sqrt` for a half, else a float power."""
function power(x, d)
    q = denominator(d)
    base = q == 1 ? x : q == 2 ? :(sqrt($x)) : :($x^$(1 / q))
    d < 0 ? :(inv($base)) : base
end

"""A coefficient as a literal: an `Int`, or a `Float64` beyond `Int`'s range."""
literal(c) = typemin(Int) <= c <= typemax(Int) ? Int(c) : Float64(c)

"""Greedy Horner form of the terms `ts`, each `(coefficient, local => power)`: take out the
local most terms hold (ties by name), then recurse on the terms with it and on the rest,
with binary `+` and `*` only."""
function horner(ts)
    c0 = sum((c for (c, a) in ts if isempty(a)); init = 0)
    rest = [(c, a) for (c, a) in ts if !isempty(a)]
    isempty(rest) && return literal(c0)
    held = Dict{Symbol, Int}()
    foreach(((_, a),) -> foreach(s -> held[s] = get(held, s, 0) + 1, keys(a)), rest)
    best = first(sort!(collect(keys(held)); by = s -> (-held[s], string(s))))
    lower(a) = (b = copy(a); b[best] == 1 ? delete!(b, best) : (b[best] -= 1); b)
    inner = horner([(c, lower(a)) for (c, a) in rest if haskey(a, best)])
    ex = inner == 1 ? best : :($best * $inner)
    others = [(c, a) for (c, a) in rest if !haskey(a, best)]
    iszero(c0) || push!(others, (c0, Dict{Symbol, Int}()))
    isempty(others) ? ex : :($ex + $(horner(others)))
end

"""`p` as code, `atoms(mono)` giving each monomial's locals: the Horner form of its positive
terms minus that of its negative ones."""
function emit(p::Poly, atoms)
    ts = [(c, atoms(mo)) for (mo, c) in sort!(collect(p.terms); by = first)]
    pos, neg = [t for t in ts if t[1] > 0], [(-c, a) for (c, a) in ts if c < 0]
    isempty(neg) ? horner(pos) : isempty(pos) ? :(-$(horner(neg))) :
        :($(horner(pos)) - $(horner(neg)))
end

"""
    body(p) -> Expr

The body of `rate_equation` for program `p`: each variable it reads, in order, loaded by
name into a positional local (an input) or computed by greedy Horner (a local), followed by
its reciprocal and root if read; then `E_total·num/den`. `x^e` reads `x` to the integer
part of `e` and `x^(±1/q)` to the rest (`x^(3/2) = x·sqrt(x)`).
"""
function body(p::Program)
    used = Set{Tuple{Int, Exp}}()
    function atoms(mo)
        a = Dict{Symbol, Int}()
        for (v, e) in mo
            n = trunc(Int, e)
            r = e - n
            for (d, k) in ((sign(n) // 1, abs(n)),
                           (sign(r) // denominator(r), numerator(abs(r))))
                k == 0 && continue
                push!(used, (v, d))
                a[loc(v, d)] = k
            end
        end
        a
    end
    def = Dict(v => emit(q, atoms) for (v, q) in p.defs)
    num, den = emit(p.num, atoms), emit(p.den, atoms)
    lines = Expr[]
    for v in sort!(unique(first.(used)))
        x = get(def, v) do
            :(getfield($(p.vars.conc[v] ? :concs : :params), $(QuoteNode(p.vars.names[v]))))
        end
        push!(lines, :($(loc(v, 1)) = $x))
        for d in sort!([d for (w, d) in used if w == v && d != 1])
            push!(lines, :($(loc(v, d)) = $(power(loc(v, 1), d))))
        end
    end
    Expr(:block, lines..., :(return getfield(params, :E_total) * $num / $den))
end

"""The program of `m`'s rate law, plain or allosteric."""
program(m::Mechanism) =
    (sol = solve(m); program(m, sol, is_allosteric(m) ? allo_law(m, sol) : law(m, sol)))

"""
    rate_equation(em, concs, params)

The rate of `em` at the concentrations `concs` and the parameters `params` (`NamedTuple`s
read by name; `params` holds `parameters(em)`). The generator decodes the Sig, derives the
law and emits it as straight-line code: allocation-free, no loop, no division by a
polynomial.
"""
@generated function rate_equation(::EnzymeMechanism{S}, concs::NamedTuple,
                                  params::NamedTuple) where {S}
    body(program(decode(S...)))
end

rate_equation(m::Mechanism, concs, params) = rate_equation(compile(m), concs, params)

"""`fitted_params`, `parameters` and `metabolites` of a compiled mechanism as constant
tuples, so that the `NamedTuple`s a fitter builds from them are type-stable."""
@generated fitted_params(::EnzymeMechanism{S}) where {S} = fitted_params(decode(S...))
@generated parameters(::EnzymeMechanism{S}) where {S} = parameters(decode(S...))
@generated metabolites(::EnzymeMechanism{S}) where {S} = Tuple(metabolites(decode(S...)))
