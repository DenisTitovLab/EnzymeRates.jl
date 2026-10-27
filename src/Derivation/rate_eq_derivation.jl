# ABOUTME: The parameter API (parameters, fitted_params, all_parameters) and the constant
# ABOUTME: names of each conformation's law: edge masks and edge names.

# ─── Parameters API ─────────────────────────────────────────

# Every mechanism is a singleton type of its own, so a method specialized on one
# compiles again for each mechanism. The derivation helpers that lift a singleton,
# or its type, to its concrete form therefore take it `@nospecialize`.

"""
    parameters(m)

The names of the parameters the rate law of `m` takes, as a tuple of Symbols: its
independent constants (`fitted_params`), then `Keq` and `E_total`. The rate law computes
every other constant from these (see [`all_parameters`](@ref)).
"""
@generated function parameters(::M) where {M <: Mechanisms.AbstractEnzymeMechanism}
    (_solve(M).fitted..., :Keq, :E_total)
end

"""Independent rate constant names for fitting (excludes Keq, E_total)."""
@generated function fitted_params(::M) where {M <: Mechanisms.AbstractEnzymeMechanism}
    Tuple(_solve(M).fitted)
end

# Convenience methods for the concrete Mechanism / AllostericMechanism
# working representation: lift to the compiled singleton via
# `compile_mechanism`, then dispatch to the singleton method. Callers that
# hold a `Mechanism` (e.g. the `identify_rate_equation` pipeline) use these
# without an explicit lift. `compile_mechanism` does Sig conversion (it
# allocates) — it is NOT the rate_equation hot path. Pass an
# `EnzymeMechanism` / `AllostericEnzymeMechanism` for the allocation-free
# `rate_equation` hot path.
parameters(m::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}) =
    parameters(Mechanisms.compile_mechanism(m))
fitted_params(m::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}) =
    fitted_params(Mechanisms.compile_mechanism(m))

"""
    all_parameters(m, params) -> NamedTuple

Return `params` followed by the dependent constants of `m`'s rate law that it lacks, in
name order, each computed from the independent constants and `Keq` in `params` by its
Haldane or Wegscheider relation (the relations [`rate_equation_string`](@ref) prints).
`params` holds the values of [`parameters`](@ref)`(m)` by name; `E_total` may be left
out. A dependent constant already in `params` is recomputed in place. An error names
every independent constant or `Keq` that `params` lacks.

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

julia> all_parameters(m, (K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0,
                          Keq = 4.0))
(K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0, Keq = 4.0, k_EP_to_ES = 0.375)
```
"""
function all_parameters(m::Union{Mechanisms.AbstractEnzymeMechanism, Mechanisms.Mechanism,
                                  Mechanisms.AllostericMechanism},
                        params::NamedTuple)
    s = _solve(m isa Mechanisms.AbstractEnzymeMechanism ? typeof(m) : _well_posed(m))
    absent = [k for k in (s.fitted..., :Keq) if !haskey(params, k)]
    isempty(absent) || error("all_parameters: params must hold every name of " *
                             "parameters(m) except E_total; it lacks " *
                             join(absent, ", "))
    # A negative integer exponent inverts the positive power; a fractional one is taken
    # in the value's own precision, so a BigFloat value keeps its full precision.
    pow(x, e) = !isone(denominator(e)) ?
                    x^(oftype(float(x), numerator(e)) / denominator(e)) :
                numerator(e) >= 0 ? x^numerator(e) : inv(x^-numerator(e))
    value(d) = prod((pow(params[s.names[f]], e) for (f, e) in d.factors);
                    init = pow(params.Keq, d.keq))
    deps = sort!(filter(d -> d.label !== :kinetic_group, s.dependents);
                 by = d -> s.names[d.col])
    merge(params, NamedTuple{Tuple(s.names[d.col] for d in deps)}(
        Tuple(value(d) for d in deps)))
end

# ─── Edge masks and constant names ──────────────────────────

"""The edges of `m`'s law in conformation `state`: every edge for `:None`; for `:A` the
active conformation's (`_active_mask`), for `:I` the inactive one's
(`_inactive_mask`)."""
_state_mask(m::Mechanisms.Mechanism, ::Symbol) = trues(length(m.edges))
_state_mask(am::Mechanisms.AllostericMechanism, state::Symbol) =
    state === :I ? Mechanisms._inactive_mask(am) : Mechanisms._active_mask(am)

"""
The parameter Symbols of each edge's constants in conformation `state` of `m` — `[K]` at
rapid equilibrium, `[kf, kr]` at steady state — through the naming chokepoint
(`name(p, m, kept)`, `kept = _kept_members(m)`): `names[k]` for edge `k` of the state's
law (`_state_mask`), empty for every other edge. `state` is `:None` for a plain
mechanism and `:A` or `:I` for an allosteric one, an `:EqualAI` group's constants taking
its one shared name in both.
"""
function _edge_names(m::Mechanisms._AnyMech, state::Symbol,
                     kept = Mechanisms._kept_members(m))
    cm = Mechanisms.catalytic(m)
    mask = _state_mask(m, state)
    tag(g) = m isa Mechanisms.AllostericMechanism &&
             Mechanisms.cat_allo_state(m, g) === :EqualAI ? :EqualAI : state
    [mask[k] ? Symbol[Mechanisms.name(p, m, kept)
                      for p in Mechanisms._edge_constants(cm, k, tag(cm.group[k]))] :
               Symbol[]
     for k in eachindex(cm.edges)]
end

# ─── Steady-state rate-constant names ──────────────────────────

"""
The names of the steady-state rate constants of `em`'s law (`_edge_names`): for an
`AllostericEnzymeMechanism`, those of both conformations, the inactive one's `I_` names
included. `_scale_rate_constants` scales exactly these, leaving the rapid-equilibrium
constants, `Keq`, `E_total`, `L` and the regulator constants alone.
"""
function _ss_rate_constant_names(@nospecialize(em::Mechanisms.AbstractEnzymeMechanism))
    m = _derivation_input(typeof(em))
    edges = Mechanisms.catalytic(m).edges
    states = m isa Mechanisms.Mechanism ? (:None,) : (:A, :I)
    kept = Mechanisms._kept_members(m)
    Set{Symbol}(n for st in states for (k, ns) in enumerate(_edge_names(m, st, kept))
                if !edges[k].rapid_equilibrium for n in ns)
end
