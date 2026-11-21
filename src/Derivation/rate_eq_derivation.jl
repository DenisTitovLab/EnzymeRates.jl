# ABOUTME: The parameter API: parameters, the constants a fit estimates, parameter_info
# ABOUTME: for Fitting, and dependent_parameter_values, the constants the law computes.

# ─── Parameters API ─────────────────────────────────────────

# Every mechanism is a singleton type of its own, so a method specialized on one
# compiles again for each mechanism. The derivation helpers that lift a singleton,
# or its type, to its concrete form therefore take it `@nospecialize`.

"""
    parameters(m) -> Tuple{Vararg{Symbol}}

The names of the constants of the rate law of `m` that a fit estimates, as a tuple of
Symbols: the independent constants in name order, then an allosteric mechanism's regulator
constants and `L`. These are the names `params` holds in [`rate_equation`](@ref); the
equilibrium constant `Keq` is not among them, since it is given, never fitted. The rate law
computes every other constant from these and `Keq` (see
[`dependent_parameter_values`](@ref)). Errors when the rate law of `m` is not defined
(`check_valid_rate_equation`).
"""
@generated function parameters(::M) where {M <: Mechanisms.AbstractEnzymeMechanism}
    Tuple(_solve(M).fitted)
end

# The concrete Mechanism / AllostericMechanism reads the joint solve of the checked
# mechanism through the memo (`_memo_solve`), compiling nothing. Callers that hold a
# `Mechanism` (e.g. the `identify_rate_equation` pipeline) use it without an explicit
# lift. Pass an `EnzymeMechanism` / `AllostericEnzymeMechanism` for the allocation-free
# `rate_equation` hot path.
function parameters(m::_AnyMech)
    check_valid_rate_equation(m)
    Tuple(_memo_solve(m).fitted)
end

"""
    parameter_info(m) -> Vector{Parameter}

One `Parameter` per fitted parameter of the concrete mechanism `m`, in `parameters(m)`
order: each independent column, then the regulator parameters, then
`AllostericParameter()`. `name(p, m)` gives each one's name. Errors when the rate law of
`m` is not defined (`check_valid_rate_equation`).
"""
function parameter_info(m::_AnyMech)
    check_valid_rate_equation(m)
    copy(_memo_solve(m).info)
end

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
