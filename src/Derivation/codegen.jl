# ABOUTME: Code generation: compile_rate_equation's generated call emits the derived law
# ABOUTME: as straight-line locals in greedy Horner form, reading its inputs by name.

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
    Q, D = _state!(vars, defs, m, s, L, trues(length(m.steps)), false)
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
    vars, defs = copy(A.vars), Pair{Int, _Poly}[]
    QA, DA = _state!(vars, defs, m, s, A.active,
                     Mechanisms.steps_in_active_conformation(am), false)
    QI, DI = _flip_only(A.inactive) ? (zero(_Poly), one(_Poly)) :
             _state!(vars, defs, m, s, A.inactive,
                     Mechanisms.steps_in_inactive_conformation(am), true)
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
(`_body`) of its law (`_with_law`)."""
_rate_body(m::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}) =
    _with_law((m, s, L) -> _body(_program(m, s, L)), m)

# ─── compile_rate_equation ─────────────────────────────────────────────────────
#
# A mechanism's Sig, the type parameter of its compiled rate law, every entry an `Int`:
#   names       their number, then each name as its number of characters and each
#               character's code point: the sorted distinct names of its metabolites, its
#               conformation labels and its allosteric states. After them a name is its
#               index among them, a list is its length followed by its entries, and a
#               metabolite is its name followed by its modifier flag.
#   catalytic   its substrates and its products; the number of forms, then each form
#               (label, bound, added, subtracted); the number of steps, then each step
#               (from, to, taken up, given off, flag, group).
#   allostery   one allosteric state per kinetic group, none for a plain mechanism; the
#               oligomeric state; the number of sites, then each site (ligands,
#               multiplicity, states).

"""The allosteric states of `m`, one per kinetic group, its oligomeric state and its
regulatory sites: none, 1 and none for a plain mechanism."""
_allostery(::Mechanisms.Mechanism) = (Symbol[], 1, Mechanisms.RegulatorySite[])
_allostery(am::Mechanisms.AllostericMechanism) =
    (am.allosteric_states, am.oligomeric_state, am.regulatory_sites)

"""`m` as one tuple of integers, its Sig (layout above): the only type parameter of its
compiled rate law. Equal mechanisms have equal Sigs, since a mechanism is stored in
canonical form."""
function _sig_of(m::_AnyMech)
    cm = Mechanisms.catalytic(m)
    states, oligomeric_state, sites = _allostery(m)
    names = sort!(unique!([Mechanisms.metabolites(m); [n.conformation for n in cm.nodes];
                           states; [x for s in sites for x in s.allosteric_states]]))
    out = Int[length(names)]
    for x in names
        text = String(x)
        push!(out, length(text))
        foreach(c -> push!(out, Int(c)), text)
    end
    index(x::Symbol) = searchsortedfirst(names, x)
    function symbols!(list::Vector{Symbol})
        push!(out, length(list))
        foreach(x -> push!(out, index(x)), list)
    end
    function metabolites!(list::Vector{Mechanisms.Metabolite})
        push!(out, length(list))
        foreach(x -> push!(out, index(x.name), x.modifier), list)
    end
    symbols!(cm.substrates)
    symbols!(cm.products)
    push!(out, length(cm.nodes))
    for n in cm.nodes
        push!(out, index(n.conformation))
        metabolites!(n.bound)
        symbols!(n.added)
        symbols!(n.subtracted)
    end
    push!(out, length(cm.steps))
    for e in cm.steps
        push!(out, e.from, e.to)
        metabolites!(e.taken_up)
        metabolites!(e.given_off)
        push!(out, e.rapid_equilibrium, e.group)
    end
    symbols!(states)
    push!(out, oligomeric_state, length(sites))
    for s in sites
        symbols!(s.ligands)
        push!(out, s.multiplicity)
        symbols!(s.allosteric_states)
    end
    Tuple(out)
end

"""A Sig read from its start: a struct, not captured variables, so nothing is boxed."""
mutable struct _SigReader
    sig::Vector{Int}
    at::Int
end
_next!(r::_SigReader) = r.sig[r.at += 1]

"""The mechanism whose Sig is `sig`, splatted, built by the constructors, which check it:
a plain mechanism when the Sig holds no allosteric state. The varargs method compiles
once, whatever the Sig length, as long as callers splat a Sig whose length their code
leaves open."""
function _mechanism_from_sig(sig::Int...)
    r = _SigReader(collect(sig), 0)
    names = [Symbol(String([Char(_next!(r)) for _ in 1:_next!(r)])) for _ in 1:_next!(r)]
    name() = names[_next!(r)]
    symbols() = Symbol[name() for _ in 1:_next!(r)]
    metabolites() = Mechanisms.Metabolite[Mechanisms.Metabolite(name(), _next!(r) == 1)
                                          for _ in 1:_next!(r)]
    substrates = symbols()
    products = symbols()
    forms = [Mechanisms.EnzymeForm(name(), metabolites(), symbols(), symbols())
             for _ in 1:_next!(r)]
    steps = [Mechanisms.Step(_next!(r), _next!(r), metabolites(), metabolites(),
                             _next!(r) == 1, _next!(r))
             for _ in 1:_next!(r)]
    m = Mechanisms.Mechanism(substrates, products, forms, steps)
    states = symbols()
    oligomeric_state = _next!(r)
    sites = [Mechanisms.RegulatorySite(symbols(), _next!(r), symbols())
             for _ in 1:_next!(r)]
    isempty(states) ? m :
        Mechanisms.AllostericMechanism(m, states, oligomeric_state, sites)
end

"""The rate law of the mechanism whose Sig (`_sig_of`) is `Sig`, as a function
(`compile_rate_equation`). Its only parameter is the Sig, so equal mechanisms share one
instance and its call is generated once per mechanism."""
struct _CompiledRateEquation{Sig} end

"""The mechanism whose compiled rate law has the type `R`, decoded from its Sig
(`_mechanism_from_sig`) and checked as `check_valid_rate_equation` checks it, without the
memo (`_well_posed`): what the law's generated methods derive from. It reads the Sig off
`R` at run time, so it and the decoder compile once, whatever the mechanism."""
_derivation_input(@nospecialize(R::Type{<:_CompiledRateEquation})) =
    _well_posed(_mechanism_from_sig(R.parameters[1]...))

"""
    compile_rate_equation(m) -> rate

The rate law of mechanism `m` as a function: `rate(concs, params, Keq)` returns the net
reaction rate per unit enzyme at the metabolite concentrations `concs` and the constants
`params` (both `NamedTuple`s), with the equilibrium constant `Keq`. Multiply by the enzyme
concentration for an absolute rate. The law is derived from the mechanism's elementary
steps by the King–Altman/Cha method, so one call gives the turnover of the whole
mechanism, not the rate of one step.

`params` holds the constants a fit estimates, [`parameters`](@ref)`(rate)`, the same names
as `parameters(m)`, and `concs` the concentrations of [`metabolites`](@ref)`(rate)`; both
are tuples of Symbols. The law computes each dependent constant from `params` and `Keq` by
its Haldane or Wegscheider relation. Both `NamedTuple`s are read by name, so their fields
may come in any order, and fields the law does not read are ignored.

The first call derives the law and compiles it to straight-line arithmetic, once per
mechanism and session; every later call is fast and allocates nothing, so compile once,
outside a loop. Equal mechanisms give the same `rate`. Calling `rate` errors when the rate
law of `m` is not defined (`check_valid_rate_equation`). [`rate_equation_string`](@ref)
prints the law.

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S ⇌ E(S)
               E(S) <--> E(P)
               E(P) ⇌ E + P
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
compile_rate_equation(m::_AnyMech) = _CompiledRateEquation{_sig_of(m)}()

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
    Tuple(Mechanisms.metabolites(_mechanism_from_sig(compiled_law.parameters[1]...)))

# A compiled law prints one line whatever its mechanism; `@nospecializeinfer` keeps a caller
# compiled for one law, such as `sprint`, from inferring a show for that law.
Base.@nospecializeinfer Base.show(io::IO, @nospecialize(rate::_CompiledRateEquation)) =
    print(io, "compiled rate equation")
