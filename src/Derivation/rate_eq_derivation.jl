# ABOUTME: The parameter API: parameters, the constants a fit estimates, parameter_units,
# ABOUTME: their units, and dependent_parameter_values, the constants the law computes.

# ─── Parameters API ─────────────────────────────────────────

# Every mechanism is a singleton type of its own, so a method specialized on one
# compiles again for each mechanism. The derivation helpers that lift a singleton,
# or its type, to its concrete form therefore take it `@nospecialize`.

"""
    parameters(m) -> Tuple{Vararg{Symbol}}
    parameters(rate) -> Tuple{Vararg{Symbol}}

The names of the constants of the rate law of `m` that a fit estimates, as a tuple of
Symbols: the independent constants in name order, then an allosteric mechanism's regulator
constants and `L`. These are the names `params` holds in a call of
[`compile_rate_equation`](@ref)`(m)`, whose result `rate` gives the same names as a
compile-time constant; the equilibrium constant `Keq` is not among them, since it is
given, never fitted. The rate law computes every other constant from these and `Keq` (see
[`dependent_parameter_values`](@ref)). Errors when the rate law of `m` is not defined
(`check_valid_rate_equation`).
"""
@generated function parameters(::M) where {M <: Mechanisms.AbstractEnzymeMechanism}
    Tuple(_solve(M).fitted)
end

# The concrete Mechanism / AllostericMechanism reads the joint solve of the checked
# mechanism through the memo (`_memo_solve`), compiling nothing. Callers that hold a
# `Mechanism` (e.g. the `identify_rate_equation` pipeline) use it without an explicit
# lift.
function parameters(m::_AnyMech)
    check_valid_rate_equation(m)
    Tuple(_memo_solve(m).fitted)
end

"""
    parameter_units(m) -> (; rate::Vector{Int}, concentrations::Matrix{Int})

The units of the constants of [`parameters`](@ref)`(m)`, as integer powers: constant
`parameters(m)[i]` has the units of the rate per unit enzyme to the power `rate[i]` times
the concentration of each metabolite `metabolites(m)[j]` to the power
`concentrations[i, j]`. A steady-state step's rate constant has `rate` 1 and −1 for each
metabolite the step takes up in its direction. Every other constant has `rate` 0: a
rapid-equilibrium binding's dissociation constant +1 for the metabolite that binds, any
other rapid-equilibrium constant +1 for each metabolite given off and −1 for each taken
up, a regulator's constant +1 for its ligand, and `L` none. A reactant's modifier copy
counts as the reactant. Errors when the rate law of `m` is not defined
(`check_valid_rate_equation`).
"""
function parameter_units(m::_AnyMech)
    check_valid_rate_equation(m)
    _parameter_units(m, _memo_solve(m).info)
end
function parameter_units(@nospecialize(m::Mechanisms.AbstractEnzymeMechanism))
    cm = _derivation_input(typeof(m))
    _parameter_units(cm, _solve(cm).info)
end

"""`parameter_units` of the checked mechanism `m`, whose fitted constants are `info`."""
function _parameter_units(m::_AnyMech, info::Vector{_Parameter})
    column = Dict(x => j for (j, x) in enumerate(Mechanisms.metabolites(m)))
    units = (; rate = zeros(Int, length(info)),
             concentrations = zeros(Int, length(info), length(column)))
    for (i, p) in enumerate(info)
        _write_units!(units, i, p, m, column)
    end
    units
end

"""Writes the units of `p`, the `i`-th fitted constant of `m`, into row `i` of `units`
(`parameter_units`), whose rows start at zero; `column` maps each metabolite of `m` to its
column."""
function _write_units!(units, i::Int, p::_StepParameter, m::_AnyMech,
                       column::Dict{Symbol, Int})
    e = Mechanisms.catalytic(m).steps[p.step]
    add!(xs, power) = foreach(x -> units.concentrations[i, column[x.name]] += power, xs)
    if e.rapid_equilibrium
        # An equilibrium constant is what the direction it is named after gives off over
        # what it takes up: a binding's dissociation constant has what binds.
        named = _constant_named_in_reverse(e) ? -1 : 1
        add!(e.given_off, named)
        add!(e.taken_up, -named)
    else
        units.rate[i] = 1
        add!(p.slot == 1 ? e.taken_up : e.given_off, -1)
    end
    nothing
end
_write_units!(units, i::Int, p::_RegulatorParameter, m::_AnyMech,
              column::Dict{Symbol, Int}) =
    (units.concentrations[i, column[p.ligand]] += 1; nothing)
_write_units!(units, i::Int, ::_AllostericParameter, m::_AnyMech,
              column::Dict{Symbol, Int}) = nothing

"""The dependent constants of the joint solve `s` in name order, less the kinetic-group
tie pins, which are not constants of the rate law."""
_dependents(s) =
    sort!(filter(d -> d.label !== :kinetic_group, s.dependents); by = d -> s.names[d.col])

"""
    dependent_parameter_values(m, params::NamedTuple, Keq::Real) -> NamedTuple

The values of the dependent constants of `m`'s rate law, by name in name order, each
computed from the fitted constants in `params` and the equilibrium constant `Keq` by its
Haldane or Wegscheider relation (the relations [`rate_equation_string`](@ref) prints).
`params` holds the values of [`parameters`](@ref)`(m)` by name; other names are ignored. An
error names every constant of `parameters(m)` that `params` lacks.
`merge(params, dependent_parameter_values(m, params, Keq))` is every constant.

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

julia> p = (K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0);

julia> dependent_parameter_values(m, p, 4.0)
(k_EP_to_ES = 0.375,)
```
"""
function dependent_parameter_values(
        m::Union{Mechanisms.AbstractEnzymeMechanism, Mechanisms.Mechanism,
                 Mechanisms.AllostericMechanism},
        params::NamedTuple, Keq::Real)
    s = _solve(m isa Mechanisms.AbstractEnzymeMechanism ? typeof(m) : _well_posed(m))
    absent = Symbol[k for k in s.fitted if !haskey(params, k)]
    isempty(absent) || error("dependent_parameter_values: params must hold every name " *
                             "of parameters(m); it lacks " * join(absent, ", "))
    # A negative integer exponent inverts the positive power; a fractional one is taken
    # in the value's own precision, so a BigFloat value keeps its full precision.
    pow(x, e) = !isone(denominator(e)) ?
                    x^(oftype(float(x), numerator(e)) / denominator(e)) :
                numerator(e) >= 0 ? x^numerator(e) : inv(x^-numerator(e))
    value(d) = prod((pow(params[s.names[f]], e) for (f, e) in d.factors);
                    init = pow(Keq, d.keq))
    deps = _dependents(s)
    NamedTuple{Tuple(s.names[d.col] for d in deps)}(Tuple(value(d) for d in deps))
end
