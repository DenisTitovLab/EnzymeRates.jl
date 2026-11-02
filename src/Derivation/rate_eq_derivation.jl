# ABOUTME: The parameter API (parameters and its kinds, dependent_parameter_values) and the
# ABOUTME: constant names of each conformation's law: edge names, rate-constant names.

# ─── Parameters API ─────────────────────────────────────────

# Every mechanism is a singleton type of its own, so a method specialized on one
# compiles again for each mechanism. The derivation helpers that lift a singleton,
# or its type, to its concrete form therefore take it `@nospecialize`.

"""The parameters `parameters(m)` lists: the independent constants and `Keq`, `E_total`."""
struct IndependentParameters end

"""The independent constants of a rate law, the ones a fit estimates."""
struct FittedParameters end

"""Every constant of a rate law: the independent ones, then the dependent ones."""
struct AllParameters end

"""The steady-state rate constants of a rate law, those a common factor scales."""
struct RateConstants end

const _ParameterKind = Union{IndependentParameters, FittedParameters, AllParameters}

"""
    parameters(m) -> Tuple{Vararg{Symbol}}
    parameters(m, kind) -> Tuple{Vararg{Symbol}}

The names of the parameters of the rate law of `m`, as a tuple of Symbols. `kind` picks
which:

  - `IndependentParameters()` (the default): the independent constants, then `Keq` and
    `E_total`. These are the parameters `rate_equation` takes.
  - `FittedParameters()`: the independent constants alone (no `Keq`, no `E_total`).
  - `AllParameters()`: the `IndependentParameters()` names, then the dependent constants
    in name order. The rate law computes each dependent constant from the independent
    ones by its Haldane or Wegscheider relation (see [`dependent_parameter_values`](@ref)).
"""
parameters(m) = parameters(m, IndependentParameters())

"""The dependent constants of the joint solve `s` in name order, less the kinetic-group
tie pins, which are not constants of the rate law."""
_dependents(s) =
    sort!(filter(d -> d.label !== :kinetic_group, s.dependents); by = d -> s.names[d.col])

"""The names of the `kind` parameters in the joint solve `s`."""
_parameter_names(s, ::FittedParameters) = Tuple(s.fitted)
_parameter_names(s, ::IndependentParameters) = (s.fitted..., :Keq, :E_total)
_parameter_names(s, ::AllParameters) =
    (_parameter_names(s, IndependentParameters())...,
     (s.names[d.col] for d in _dependents(s))...)

@generated function parameters(::M, ::K) where {M <: Mechanisms.AbstractEnzymeMechanism,
                                                K <: _ParameterKind}
    _parameter_names(_solve(M), K())
end

# The concrete Mechanism / AllostericMechanism reads the joint solve of the checked
# mechanism through the memo (`_memo_solve`), compiling nothing. Callers that hold a
# `Mechanism` (e.g. the `identify_rate_equation` pipeline) use it without an explicit
# lift. Pass an `EnzymeMechanism` / `AllostericEnzymeMechanism` for the allocation-free
# `rate_equation` hot path.
function parameters(m::_AnyMech, kind::_ParameterKind)
    check_valid_rate_equation(m)
    _parameter_names(_memo_solve(m), kind)
end

"""
    dependent_parameter_values(m, params::NamedTuple) -> NamedTuple

The values of the dependent constants of `m`'s rate law, in name order (the ones
[`parameters`](@ref)`(m, AllParameters())` lists after the independent ones), each computed
from the independent constants and `Keq` in `params` by its Haldane or Wegscheider relation
(the relations [`rate_equation_string`](@ref) prints). `params` holds the values of
[`parameters`](@ref)`(m)` by name; `E_total` may be left out, and a dependent constant
`params` holds is ignored. An error names every independent constant or `Keq` that `params`
lacks. `merge(params, dependent_parameter_values(m, params))` is every constant.

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

julia> p = (K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0, Keq = 4.0);

julia> dependent_parameter_values(m, p)
(k_EP_to_ES = 0.375,)
```
"""
function dependent_parameter_values(
        m::Union{Mechanisms.AbstractEnzymeMechanism, Mechanisms.Mechanism,
                 Mechanisms.AllostericMechanism},
        params::NamedTuple)
    s = _solve(m isa Mechanisms.AbstractEnzymeMechanism ? typeof(m) : _well_posed(m))
    absent = [k for k in (s.fitted..., :Keq) if !haskey(params, k)]
    isempty(absent) || error("dependent_parameter_values: params must hold every name " *
                             "of parameters(m) except E_total; it lacks " *
                             join(absent, ", "))
    # A negative integer exponent inverts the positive power; a fractional one is taken
    # in the value's own precision, so a BigFloat value keeps its full precision.
    pow(x, e) = !isone(denominator(e)) ?
                    x^(oftype(float(x), numerator(e)) / denominator(e)) :
                numerator(e) >= 0 ? x^numerator(e) : inv(x^-numerator(e))
    value(d) = prod((pow(params[s.names[f]], e) for (f, e) in d.factors);
                    init = pow(params.Keq, d.keq))
    deps = _dependents(s)
    NamedTuple{Tuple(s.names[d.col] for d in deps)}(Tuple(value(d) for d in deps))
end

# ─── Edge masks and constant names ──────────────────────────

"""
The parameter Symbols of each edge's constants in conformation `state` of `m` — `[K]` at
rapid equilibrium, `[kf, kr]` at steady state — as the joint solve names its columns
(`constant_columns`): `names[k]` for edge `k` of the state's law (`active_steps`,
`inactive_steps`; every edge for `:None`), empty for every other edge. `state` is
`:None` for a plain mechanism and `:A` or `:I` for an allosteric one; an inactive
constant reads its own column when its group has one, so an `:EqualAI` group's
constants take their one shared name in both.
"""
function _edge_names(m::_AnyMech, state::Symbol, cols = constant_columns(m))
    cm = Mechanisms.catalytic(m)
    mask = _conformation_steps(m, state === :I ? 2 : 1)
    index = Dict(c => i for (i, c) in enumerate(cols.columns))
    column(k, slot) = get(() -> index[_Column((k, slot, false))], index,
                          _Column((k, slot, state === :I)))
    [mask[k] ? Symbol[cols.names[column(k, slot)]
                      for slot in 1:(cm.steps[k].rapid_equilibrium ? 1 : 2)] :
               Symbol[]
     for k in eachindex(cm.steps)]
end

# ─── Steady-state rate-constant names ──────────────────────────

"""
The names of the steady-state rate constants of `m`'s law (`_edge_names`), as a tuple: for
an `AllostericMechanism`, those of both conformations, the inactive one's `I_` names
included. They are the constants a common factor scales, leaving the rapid-equilibrium
constants, `Keq`, `E_total`, `L` and the regulator constants alone.
"""
function _rate_constant_names(m::_AnyMech)
    cm = Mechanisms.catalytic(m)
    states = m isa Mechanisms.Mechanism ? (:None,) : (:A, :I)
    cols = constant_columns(m)
    Tuple(unique(n for st in states for (k, ns) in enumerate(_edge_names(m, st, cols))
                 if !cm.steps[k].rapid_equilibrium for n in ns))
end

@generated function parameters(::M, ::RateConstants) where {
        M <: Mechanisms.AbstractEnzymeMechanism}
    _rate_constant_names(_derivation_input(M))
end

parameters(m::_AnyMech, ::RateConstants) = _rate_constant_names(_well_posed(m))
