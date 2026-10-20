# ABOUTME: King-Altman/Cha rate equation derivation via @generated functions.
# ABOUTME: Parameters API and MWC allosteric rate assembly.

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
@generated function parameters(::M) where {M <: AbstractEnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    (indep..., :Keq, :E_total)
end

"""Independent rate constant names for fitting (excludes Keq, E_total)."""
@generated function fitted_params(::M) where {M <: AbstractEnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    indep
end

# Convenience methods for the concrete Mechanism / AllostericMechanism
# working representation: lift to the compiled singleton via
# `compile_mechanism`, then dispatch to the singleton method. Callers that
# hold a `Mechanism` (e.g. the `identify_rate_equation` pipeline) use these
# without an explicit lift. `compile_mechanism` does Sig conversion (it
# allocates) — it is NOT the rate_equation hot path. Pass an
# `EnzymeMechanism` / `AllostericEnzymeMechanism` for the allocation-free
# `rate_equation` hot path.
parameters(m::Union{Mechanism, AllostericMechanism}) = parameters(compile_mechanism(m))
fitted_params(m::Union{Mechanism, AllostericMechanism}) =
    fitted_params(compile_mechanism(m))

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
function all_parameters(m::Union{AbstractEnzymeMechanism, Mechanism, AllostericMechanism},
                        params::NamedTuple)
    s = _solve(m isa AbstractEnzymeMechanism ? typeof(m) : _well_posed(m))
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

# ─── RE Group Helpers ───────────────────────────────────────

"""
Free enzyme of a rapid-equilibrium segment `segment` (forms of `m`): the form with the
fewest bound metabolites, tie-broken toward no covalent residual, then toward the first
in the canonical form order (the smaller index). Referencing each segment's alphas to
this form makes the derivation independent of the step order and yields the readable
`1 + [S]/K + …` form.
"""
_segment_root(m::Mechanism, segment::Vector{Int}) =
    argmin(f -> (total(m.forms[f].bound), _has_residual(m.forms[f]), f), segment)

"""
    _re_weight_ratio(m, e, K) -> POLY

w(to)/w(from) of rapid-equilibrium edge `e` of `m` as a monomial: [M]/K for a binding of
M, plain or fused (K a dissociation constant), and K·Π[consumed]/Π[released] for every
other edge (K the equilibrium constant of the stored direction, products over
reactants: [to]·Π[released] / ([from]·Π[consumed]), as its name `K_<from>_to_<to>`
says). A metabolite counts by name, so a modifier copy reads its reactant's
concentration.
"""
function _re_weight_ratio(m::Mechanism, e::Edge, K::Symbol)
    d = Dict{Symbol, Int}(K => _is_binding(e) ? -1 : 1)
    for i in entries(e.consumed); x = m.table[i].name; d[x] = get(d, x, 0) + 1; end
    for i in entries(e.released); x = m.table[i].name; d[x] = get(d, x, 0) - 1; end
    filter!(p -> p.second != 0, d)
    POLY(sort!(MONO(collect(d)); by = first) => 1)
end

"""
Compute alpha factors (relative concentrations within RE segments) as POLY values, one
per form of `m`. Reads the edges `flat` of a sub-graph of `m` in that order; each RE
edge's weight ratio comes from `_re_weight_ratio`, entering its target from its source
or the reverse. `segments` are the sub-graph's RE segments (`_segments`) and `names[k]` the
parameter Symbols of edge `k` (the K first). The guard's class 2 rejects every mechanism
whose RE edges close a cycle that turns the reaction over (`_re_turnover_transitions`), so
the weights reproduce each RE edge's own ratio.
"""
function _compute_alpha(m::Mechanism, flat::Vector{Int},
                        segments::Vector{Vector{Int}}, names::Vector{Vector{Symbol}})
    alpha = Vector{POLY}(fill(poly_one(), length(m.forms)))

    for segment in segments
        length(segment) == 1 && continue
        root = _segment_root(m, segment)
        visited = Set{Int}([root])
        queue = [root]
        while !isempty(queue)
            cur = popfirst!(queue)
            for k in flat
                e = m.edges[k]
                e.rapid_equilibrium || continue
                K = names[k][1]
                if e.from == cur && e.to ∉ visited
                    alpha[e.to] = poly_mul(alpha[cur], _re_weight_ratio(m, e, K))
                    push!(visited, e.to); push!(queue, e.to)
                elseif e.to == cur && e.from ∉ visited
                    alpha[e.from] = poly_mul(alpha[cur],
                                             _invert_monomial(_re_weight_ratio(m, e, K)))
                    push!(visited, e.from); push!(queue, e.from)
                end
            end
        end
    end
    alpha
end

"""Build rate poly for one SS step direction in Laurent (fractional) form: the rate
constant `k_poly` times the metabolites named `mets` it takes up times the weight of its
source form `i_form`."""
function _ss_contrib(k_poly, mets::Vector{Symbol}, i_form, alpha)
    r = isempty(mets) ? k_poly : poly_mul(k_poly, reduce(poly_mul, poly_sym.(mets)))
    poly_mul(r, alpha[i_form])
end

# ─── Raw Rate Equation Derivation (Unified Cha / King-Altman) ───

"""
    _raw_symbolic_rate_polys(m, mask, names) -> (num, den, d_free)

Build raw numerator and denominator POLYs for the rate equation of the edges of `m`'s
catalytic graph that `mask` keeps: all of them, or one conformation's
(`_state_mask`). `names[k]` holds the parameter Symbols of edge `k`'s constants (`[K]`
or `[kf, kr]`), which are the leaves of `num`/`den`. The law covers only the forms the
kept edges touch (`_present_forms`). Aborts, before the King–Altman expansion, a graph
whose denominator would be too large to derive (`_assert_derivable`).

Also returns `d_free`, the weight in the returned denominator of free E, form 1, where the
conformations interconvert: the spanning-tree weight `D[g_free]` of its segment times the
concentration monomial that brings `num`/`den` to lowest terms. Free E roots its segment,
so its own `alpha` is 1.
"""
function _raw_symbolic_rate_polys(m::_AnyMech, mask::BitVector,
                                  names::Vector{Vector{Symbol}})
    cm = catalytic(m)
    _assert_derivable(cm, mask)
    segments = _segments(cm, mask, _flags(cm))
    # A fully-inert conformation (every binding pruned) has no form; it exists only as
    # free enzyme — no flux, partition 1, D[g_free] 1.
    isempty(segments.segments) && return poly_zero(), poly_one(), poly_one()
    flat = [k for k in _flat_edges(cm) if mask[k]]
    # The concentration symbols: the same set `_poly_to_expr` treats as concentrations,
    # and the only symbols `_reduce_conc_lowest_terms` may shift. No parameter symbol is
    # in it, so the concentration-GCD never drops a fitted parameter.
    conc_set = Set{Symbol}(metabolites(m))
    alpha = _compute_alpha(cm, flat, segments.segments, names)
    seg = segments.segment_of
    G = length(segments.segments)

    # The Laplacian of the segment graph: each SS edge joining two segments adds its
    # forward flux out of its source segment and its reverse flux out of its target. An
    # edge within one segment enters only the numerator. `ss` keeps every SS edge with
    # its two forms and fluxes for the numerator.
    L = [poly_zero() for _ in 1:G, _ in 1:G]
    ss = Tuple{Edge, Int, Int, POLY, POLY}[]
    for k in flat
        e = cm.edges[k]
        e.rapid_equilibrium && continue
        a, b = e.from, e.to
        fwd = _ss_contrib(poly_sym(names[k][1]), _entry_names(cm, e.consumed), a, alpha)
        rev = _ss_contrib(poly_sym(names[k][2]), _entry_names(cm, e.released), b, alpha)
        push!(ss, (e, a, b, fwd, rev))
        g1, g2 = seg[a], seg[b]
        g1 == g2 && continue
        L[g1, g1] = poly_add(L[g1, g1], fwd); L[g1, g2] = poly_sub(L[g1, g2], fwd)
        L[g2, g2] = poly_add(L[g2, g2], rev); L[g2, g1] = poly_sub(L[g2, g1], rev)
    end
    # D[g], the spanning-tree weight of segment g: L without g's row and column.
    D = [(o = setdiff(1:G, g); sym_det(L[o, o])) for g in 1:G]

    # The conformations interconvert at free E, form 1, which every state's edges reach
    # (`_state_mask`). It roots its segment (`_segment_root` ranks it first), so its own
    # alpha is 1 and its weight is its segment's tree weight.
    d_free = D[seg[1]]

    den = poly_zero()
    for g in 1:G
        sigma = reduce(poly_add, (alpha[f] for f in segments.segments[g]); init=poly_zero())
        den = poly_add(den, poly_mul(sigma, D[g]))
    end

    # Numerator of the rate: v·den summed over steady-state edges. With u(f) the
    # first substrate's exponent in form f's RE weight, each SS edge e contributes
    # ω_e·(forward − reverse flux) with ω_e = (copies of the first substrate e
    # consumes − copies it releases) + u(from_e) − u(to_e). Flux conservation at every
    # form makes the sum equal the net consumption of the first substrate for every
    # parameter value, so it needs no choice of reaction cut and does not depend on
    # the direction an edge is written in. Edges with ω_e = 0 contribute nothing.
    x = first(substrates(m))
    expo(p) = (mono = only(keys(p));
               k = findfirst(q -> q.first == x, mono);
               k === nothing ? 0 : mono[k].second)
    num = poly_zero()
    for (e, a, b, fwd, rev) in ss
        ω = count(==(x), _entry_names(cm, e.consumed)) -
            count(==(x), _entry_names(cm, e.released)) + expo(alpha[a]) - expo(alpha[b])
        ω == 0 && continue
        term = poly_sub(poly_mul(fwd, D[seg[a]]), poly_mul(rev, D[seg[b]]))
        num = poly_add(num, POLY(key => ω * c for (key, c) in term))
    end

    _reduce_conc_lowest_terms(num, den, d_free, conc_set)
end

_raw_symbolic_rate_polys(m::Mechanism) =
    _raw_symbolic_rate_polys(m, trues(length(m.edges)), _edge_names(m, :None))
_raw_symbolic_rate_polys(@nospecialize(M::Type{<:EnzymeMechanism})) =
    _raw_symbolic_rate_polys(Mechanism(M()))

"""The edges of `m`'s law in conformation `state`: every edge for `:None`; for `:A` the
active conformation's (`_active_mask`), for `:I` the inactive one's
(`_inactive_mask`)."""
_state_mask(m::Mechanism, ::Symbol) = trues(length(m.edges))
_state_mask(am::AllostericMechanism, state::Symbol) =
    state === :I ? _inactive_mask(am) : _active_mask(am)

"""
The parameter Symbols of each edge's constants in conformation `state` of `m` — `[K]` at
rapid equilibrium, `[kf, kr]` at steady state — through the naming chokepoint
(`name(p, m, kept)`, `kept = _kept_members(m)`): `names[k]` for edge `k` of the state's
law (`_state_mask`), empty for every other edge. `state` is `:None` for a plain
mechanism and `:A` or `:I` for an allosteric one, an `:EqualAI` group's constants taking
its one shared name in both.
"""
function _edge_names(m::_AnyMech, state::Symbol, kept = _kept_members(m))
    cm = catalytic(m)
    mask = _state_mask(m, state)
    tag(g) = m isa AllostericMechanism && cat_allo_state(m, g) === :EqualAI ?
             :EqualAI : state
    [mask[k] ? Symbol[name(p, m, kept) for p in _edge_constants(cm, k, tag(cm.group[k]))] :
               Symbol[]
     for k in eachindex(cm.edges)]
end

"""
Estimate a mechanism's King–Altman denominator term count, V×τ, from its
catalytic segment graph WITHOUT deriving the rate equation.

`V` is the number of RE-connected segments (`_segments`) of the edges `mask` keeps; `τ` is
the spanning-tree count of the segment graph, whose edges are the SS edges: by the
Matrix–Tree theorem, the
determinant of its Laplacian with the first row and column dropped, computed exactly (`det`
of a `BigInt` matrix is fraction-free Bareiss; a disconnected graph gives 0). `V×τ` equals
the number of spanning-tree products in the denominator — the products the compiled equation
evaluates on every call — so it bounds fit cost. It guards the derivation itself
(`_assert_derivable` aborts when V×τ exceeds `MAX_RATE_EQUATION_TERMS`, before the O(G!)
symbolic cofactor expansion) and backs the `eq_complexity_filter` keyword of
[`identify_rate_equation`](@ref), which skips over-complex mechanisms before fitting. The
method for a mechanism counts every edge of its catalytic graph, which for an allosteric
mechanism with `:OnlyI` groups holds more than its active conformation's law.
"""
function _eq_complexity(m::Mechanism, mask::BitVector)
    s = _segments(m, mask, _flags(m))
    G = length(s.segments)
    G <= 1 && return big(G)                          # single RE segment ⇒ τ = 1
    lap = zeros(BigInt, G, G)
    for (k, e) in enumerate(m.edges)
        # SS edges are the segment-graph edges
        mask[k] && !e.rapid_equilibrium || continue
        g1, g2 = s.segment_of[e.from], s.segment_of[e.to]
        g1 == g2 && continue
        lap[g1, g1] += 1; lap[g2, g2] += 1; lap[g1, g2] -= 1; lap[g2, g1] -= 1
    end
    G * det(lap[2:end, 2:end])
end
_eq_complexity(m::_AnyMech) =
    _eq_complexity(catalytic(m), trues(length(catalytic(m).edges)))

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
Abort the rate-equation derivation for a mechanism whose denominator term count
(V×τ) exceeds `MAX_RATE_EQUATION_TERMS`. Called at each num/den derivation entry
point, before the O(G!) symbolic cofactor expansion in `sym_det` — a numeric
O(G³) check that throws a `RateEquationTooLargeError` quickly instead of hanging on
the factorial expansion.
"""
function _assert_derivable(m::Mechanism, mask::BitVector)
    _eq_complexity(m, mask) > MAX_RATE_EQUATION_TERMS && throw(RateEquationTooLargeError())
    return nothing
end

# ─── Expr generation from POLY ──────────────────────────────

"""The numerator and denominator Exprs of the rate equation of `m`, from its raw
symbolic rate polys: flat sums for the printed equation, or in Horner form
(`_poly_to_horner`) for the generated code when `horner` is true. Errors when the
denominator vanishes on a face of the guard (`_assert_support`)."""
function _num_den_exprs(m::Mechanism; horner = false)
    names = _edge_names(m, :None)
    num, den, _ = _raw_symbolic_rate_polys(m, trues(length(m.edges)), names)
    _assert_support(m, _law_branches(den))
    horner && return _poly_to_horner(num), _poly_to_horner(den)
    param_syms = Set{Symbol}(n for ns in names for n in ns)
    _poly_to_expr(num, param_syms), _poly_to_expr(den, param_syms)
end

"""Build destructuring Expr: (; a, b, c) = source"""
function _destructuring_expr(syms, source::Symbol)
    Expr(:(=), Expr(:tuple, Expr(:parameters, syms...)), source)
end

"""
The assignments `sym = rhs` of the dependent parameters `dep`, in name order.
The joint solve is a full Gauss–Jordan elimination (`_reduce_joint_system`), so every
right-hand side reads only `Keq` and independent parameters, and the assignments may run
in any order.
"""
_dep_assignments(dep) =
    [Expr(:(=), sym, rhs) for (sym, rhs) in sort!(collect(dep); by = first)]

"""
Error if a metabolite of `m` has the name of a parameter in `param_syms` or `dep`. The
rate equation binds the parameters, the concentrations and the dependent parameters as
local names of one scope, so a shared name would let one overwrite another. Every
parameter name holds `_` or is `Keq` or `L`, and the reaction and mechanism constructors
reject those metabolite names (`_check_metabolite_names`), so this check cannot fire for
a mechanism the constructors build.
"""
function _assert_no_name_clash(m::_AnyMech, param_syms, dep)
    for x in metabolites(m)
        (x in param_syms || haskey(dep, x)) && error(
            "Metabolite $x has the name of a parameter of this mechanism's rate " *
            "equation; rename the metabolite")
    end
end

"""
The body of the generated `rate_equation` of `m`, the mechanism its compiled type
decodes to: destructure `param_syms` from
`params` and the metabolites from `concs`, assign the dependent parameters `dep`, and
return `E_total * (num) / (den)` with `num` and `den` from `_num_den_exprs` in Horner
form. It first runs `_assert_no_name_clash`, which cannot fire for a mechanism the
constructors build.
"""
function _rate_body(m::_AnyMech, param_syms, dep)
    _assert_no_name_clash(m, param_syms, dep)
    num, den = _num_den_exprs(m; horner = true)
    Expr(:block,
        _destructuring_expr(param_syms, :params),
        _destructuring_expr(metabolites(m), :concs),
        _dep_assignments(dep)...,
        :(return E_total * ($num) / ($den)))
end

# ─── rate_equation ─────────────────────────────────────────────

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
Haldane or Wegscheider relation.

The derivation runs once, at compile time: the body is generated as a single
arithmetic expression with no allocations, loops, or matrix operations. Use
[`rate_equation_string`](@ref) to inspect the derived equation symbolically.
"""
@generated function rate_equation(m::M_type, concs::NamedTuple,
        params::NamedTuple) where {M_type <: AbstractEnzymeMechanism}
    # The body binds each metabolite's name as a local, so the type parameter's name
    # holds `_`, which no metabolite name does. The generator decodes the mechanism and
    # checks its rate law once (`_derivation_input`).
    m_decoded = _derivation_input(M_type)
    dep, indep = _dependent_param_exprs(m_decoded)
    _rate_body(m_decoded, (indep..., :Keq, :E_total), dep)
end

# Concrete-mechanism convenience: lift to the singleton (see the note on
# the parameters/fitted_params methods above — allocates, not the hot path).
rate_equation(m::Union{Mechanism, AllostericMechanism}, concs, params) =
    rate_equation(compile_mechanism(m), concs, params)

# ─── String Representation ────────────────────────────────────

"""
    rate_equation_string(m) -> String

Return the symbolic rate equation for mechanism `m` as a multi-line
`String` (it returns the text — it does not print). Pass a concrete `Mechanism` /
`AllostericMechanism` or its compiled [`EnzymeMechanism`](@ref) singleton.

The string spells out how [`rate_equation`](@ref)
evaluates: a `(; …) = params` destructure line, a `(; …) = concs`
destructure line, then the `v = E_total * (num) / (den)` line. Dependent rate
constants are listed first under `# Wegscheider constraints:` and
`# Haldane constraints:` headers — the thermodynamic identities that eliminate
parameters — and only the independent set appears in the `params` destructure.

Every name prints as a Julia identifier, the residual forms of ping-pong
mechanisms included (`k_EA_to_EP_res_plusA_minusP`), so with `params` and
`concs` bound the string runs as Julia code.

Use `print` on the result to see the multi-line layout without escaped
newlines.

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

julia> print(rate_equation_string(m))
(; K_EP_to_E_P, K_ES_to_E_S, k_ES_to_EP, Keq, E_total) = params
(; S, P) = concs
# Haldane constraints:
k_EP_to_ES = (1 / Keq) * K_EP_to_E_P * (1 / K_ES_to_E_S) * k_ES_to_EP
v = E_total * (k_ES_to_EP * S / K_ES_to_E_S - k_EP_to_ES * P / K_EP_to_E_P) / (1 + P / K_EP_to_E_P + S / K_ES_to_E_S)
```
"""
function rate_equation_string end

# Concrete-mechanism convenience: lift to the singleton.
rate_equation_string(m::Union{Mechanism, AllostericMechanism}) =
    rate_equation_string(compile_mechanism(m))

"""
The `rate_equation_string` text of `m`: the `params` destructure of `param_syms`, the
`concs` destructure, one line `sym = rhs` per dependent parameter of `dep`, and the
`v = E_total * (num) / (den)` line from `_num_den_exprs`. A dependent whose right-hand
side mentions `Keq` goes under `# Haldane constraints:`, every other one under
`# Wegscheider constraints:`. Each section is sorted by name, which is load-bearing: the
eq_hash dedup of rate-equivalent mechanisms compares these strings, so their line order
must not depend on the order the solve emits its dependents in. It first runs
`_assert_no_name_clash`, which cannot fire for a mechanism the constructors build.
"""
function _equation_text(m::_AnyMech, param_syms, dep)
    _assert_no_name_clash(m, param_syms, dep)
    weg, hal = String[], String[]
    for (sym, rhs) in sort!(collect(dep); by = first)
        push!(_mentions(rhs, :Keq) ? hal : weg, "$sym = $(_expr_to_string(rhs))")
    end
    num, den = _num_den_exprs(m)
    lines = ["(; $(join(param_syms, ", "))) = params",
             "(; $(join(metabolites(m), ", "))) = concs"]
    isempty(weg) || push!(lines, "# Wegscheider constraints:", weg...)
    isempty(hal) || push!(lines, "# Haldane constraints:", hal...)
    push!(lines, "v = E_total * ($(_expr_to_string(num))) / ($(_expr_to_string(den)))")
    join(lines, "\n")
end

function rate_equation_string(@nospecialize(m::AbstractEnzymeMechanism))
    mech = _derivation_input(typeof(m))
    dep, indep = _dependent_param_exprs(mech)
    _equation_text(mech, (indep..., :Keq, :E_total), dep)
end

# ─── Steady-state rate-constant names ──────────────────────────

"""
The names of the steady-state rate constants of `em`'s law (`_edge_names`): for an
`AllostericEnzymeMechanism`, those of both conformations, the inactive one's `I_` names
included. `_scale_rate_constants` scales exactly these, leaving the rapid-equilibrium
constants, `Keq`, `E_total`, `L` and the regulator constants alone.
"""
function _ss_rate_constant_names(@nospecialize(em::AbstractEnzymeMechanism))
    m = _derivation_input(typeof(em))
    edges = catalytic(m).edges
    states = m isa Mechanism ? (:None,) : (:A, :I)
    kept = _kept_members(m)
    Set{Symbol}(n for st in states for (k, ns) in enumerate(_edge_names(m, st, kept))
                if !edges[k].rapid_equilibrium for n in ns)
end

# ═══════════════════════════════════════════════════════════════════
# AllostericEnzymeMechanism rate equations (MWC)
#
# Per-active-site (per-protomer) normalization: `E_total` is the active-
# site concentration, so the leading `cat_n` multiplicity coefficient is
# absorbed into E_total and does not appear in the rate. The
# `Q_cat_c^(cat_n - 1)` / `Q_cat_c^cat_n` binding-statistics powers stay.
#
# MWC rate formula (per conformation c, summed over conformations).
# Let cat_n = catalytic_multiplicity(m):
#   num = sum_c( L_c * N_cat_c * Q_cat_c^(cat_n - 1)
#             * prod(Q_reg_i_c^n_reg_i for all regulatory sites i) )
#   den = sum_c( L_c * Q_cat_c^cat_n * prod(Q_reg_i_c^n_reg_i for all regulatory sites i) )
#   v = E_total * num / den
#
# Regulatory sites contribute to BOTH numerator and denominator at their
# multiplicity, regardless of whether n_reg_i matches cat_n.
# ═══════════════════════════════════════════════════════════════════

"""
    _state_rate_polys(am, state) -> (num, den, d_free)

The rate polynomials of `am`'s catalytic graph in conformation `state` (`:A` or `:I`),
in that state's names (`_edge_names`, `kept` as `_kept_members(am)`): each conformation on
its own pruned graph (`_state_mask`: the active one without its `:OnlyI` groups, the
inactive one without its `:OnlyA` groups, each without the forms that strands), each edge
in `am`'s stored direction. King–Altman on the pruned graph derives the
inactive conformation's binding partition, and for a dead cycle its steady-state fluxes
cancel, so `num` is 0. Every derivation entry first runs the boundary check
(`_check_well_posed`), whose class 3 rejects, among other faults, a face where a segment
of either conformation's law network has no bottom form.
"""
_state_rate_polys(am::AllostericMechanism, state::Symbol, kept = _kept_members(am)) =
    _raw_symbolic_rate_polys(am, _state_mask(am, state), _edge_names(am, state, kept))

"""
The regulator constants of `am` listed with conformation `state`, site by site: for `:A`
every ligand the active conformation binds, an `:EqualAI` ligand under its one bare name
(`K_<ligand>reg`, state `:EqualAI`) and an `:OnlyA` or `:NonequalAI` one under its
`K_A_` name; for `:I` the `K_I_` constant of every `:OnlyI` or `:NonequalAI` ligand.
"""
_kreg_params(am::AllostericMechanism, state::Symbol) =
    state === :A ?
    Kreg[Kreg(lig, tag === :EqualAI ? :EqualAI : :A)
         for site in regulatory_sites(am)
         for (lig, tag) in zip(ligands(site), allo_states(site)) if tag !== :OnlyI] :
    Kreg[Kreg(lig, :I) for site in regulatory_sites(am)
         for (lig, tag) in zip(ligands(site), allo_states(site))
         if tag in (:OnlyI, :NonequalAI)]

# ─── Rate body building helpers ───────────────────────────────────

"""
The dissociation-constant name of the ligand named `lig` (allosteric tag `tag`) of `am`
in the active (`inactive = false`) or inactive conformation, or `nothing` when the
ligand is absent from that conformation (`:OnlyI` in the active state, `:OnlyA` in the
inactive one). An `:EqualAI` ligand has one bare name in both conformations
(`K_<ligand>reg`); an `:OnlyA`, `:OnlyI` or `:NonequalAI` ligand has its conformation's
`K_A_` or `K_I_` name. Renders the name via the `name(::Kreg, am)` chokepoint.
"""
function _reg_K(am::AllostericMechanism, lig::Symbol, tag::Symbol, inactive::Bool)
    tag === (inactive ? :OnlyA : :OnlyI) && return nothing
    name(Kreg(lig, tag === :EqualAI ? :EqualAI : inactive ? :I : :A), am)
end

"""Build the regulatory site partition function expression: 1 + lig/K_lig_reg_i + ...
over the ligands of `site` present in the given conformation (`_reg_K`)."""
function _reg_site_expr(am::AllostericMechanism, site, inactive::Bool)
    terms = Any[1]
    for (lig, tag) in zip(ligands(site), allo_states(site))
        K = _reg_K(am, lig, tag, inactive)
        K === nothing || push!(terms, :($lig / $K))
    end
    _nest_binary(:+, terms)
end

"""Raise an expression to a positive integer power (returns `expr` itself for n=1)."""
_power_expr(expr, n::Int) = n == 1 ? expr : :(($expr)^$n)

"""MWC active + L·inactive state-combine `A + L * B` of the retained num/den sums
(`_num_den_exprs`)."""
_mwc_combine(a, b) = :($a + L * $b)

"""Cross-weight an MWC state term by the OTHER conformation's free-enzyme weight
`D_other^n` (`n = catalytic_multiplicity`). Restores a common free-enzyme basis
when the two conformations' free-enzyme weights differ and cannot be divided out (a
metabolite-bearing or multi-term `D`). A no-op when `d_other_expr == 1`."""
_mwc_cross_weight(term, d_other_expr, n) =
    d_other_expr == 1 ? term : :($(_power_expr(d_other_expr, n)) * $term)

"""Inverse of a single-term (monomial) `POLY`: negate every exponent and invert
the coefficient. Errors unless `p` is exactly one term."""
function _invert_monomial(p::POLY)
    length(p) == 1 || error("_invert_monomial: not a monomial: $p")
    mono, coef = first(p)
    POLY(_mono((s => -e for (s, e) in mono)...) => inv(coef))
end

"""True when `p` is a single term whose monomial names no symbol in `mets`.
A metabolite-free monomial free-enzyme weight can be divided out of a `POLY` as a
Laurent factor; a metabolite-bearing or multi-term `D` cannot (division would
put a concentration or a rational in a denominator), so it is cross-weighted."""
function _is_metabolite_free_monomial(p::POLY, mets::Set{Symbol})
    length(p) == 1 || return false
    mono, _ = first(p)
    !any(s in mets for (s, _) in mono)
end

"""
    _mwc_state_polys(am, mets, kept) -> (num_A, den_A, num_I, den_I, d_A, d_I)

The catalytic rate polys of both conformations of `am` (`_state_rate_polys`, with `kept`
as `_kept_members(am)`) under the
formulation-1 per-state free-enzyme normalization that `rate_equation`
(`_num_den_exprs`) uses. Each state's polys come from its conformation's pruned graph
(`_state_mask`: `_active_mask` drops the `:OnlyI` groups and `_inactive_mask` the
`:OnlyA` ones, each with every form they disconnect from free E). King–Altman on that
subgraph derives the inactive conformation's binding partition, and for a dead cycle the
pruned graph's steady-state fluxes cancel exactly, so `num_I` is 0.

The free-enzyme weights `d_A` and `d_I` combine three ways, all rendering the same
value (`mets` are the catalytic metabolites):
- `d_A == d_I`: the polys stay raw (identical conformations; the factor cancels);
- both metabolite-free monomials: each state's polys are divided by its own weight
  (the clean standard-MWC form);
- otherwise: the polys stay raw, and each state is cross-weighted by the other state's
  weight (`d^n` on the rate law's terms).
Only the third case returns the weights; the first two return `poly_one()` for both.
"""
function _mwc_state_polys(am::AllostericMechanism, mets::Set{Symbol},
                          kept::Dict{Tuple{Int, Int, Bool}, Int})
    num_A, den_A, d_A = _state_rate_polys(am, :A, kept)
    num_I, den_I, d_I = _state_rate_polys(am, :I, kept)
    if d_A == d_I
        d_A = d_I = poly_one()
    elseif _is_metabolite_free_monomial(d_A, mets) &&
           _is_metabolite_free_monomial(d_I, mets)
        inv_A, inv_I = _invert_monomial(d_A), _invert_monomial(d_I)
        num_A, den_A = poly_mul(num_A, inv_A), poly_mul(den_A, inv_A)
        num_I, den_I = poly_mul(num_I, inv_I), poly_mul(den_I, inv_I)
        d_A = d_I = poly_one()
    end
    num_A, den_A, num_I, den_I, d_A, d_I
end

"""
Assemble the MWC numerator and denominator Exprs of `am`, the mechanism its compiled type
decodes to.
Returns `(full_num, full_den)`, written with flat sums for the printed equation, or in
Horner form (`_poly_to_horner`) for the generated code when `horner` is true.
Per-active-site normalization: the
numerator carries no leading `catalytic_multiplicity` factor; only the
`Q_cat^(CatN-1)` / `Q_cat^CatN` binding-statistics powers remain. Errors when the
denominator vanishes on a face of the guard (`_assert_support`).
"""
function _num_den_exprs(am::AllostericMechanism; horner = false)
    CatN = catalytic_multiplicity(am)
    RS = regulatory_sites(am)
    kept = _kept_members(am)

    # The active conformation's catalytic names drive `_poly_to_expr`'s
    # param/metabolite ordering split; the I-poly's `I_` names sort as non-params.
    cat_params = Set{Symbol}(n for ns in _edge_names(am, :A, kept) for n in ns)
    to_expr(p) = horner ? _poly_to_horner(p) : _poly_to_expr(p, cat_params)
    num_A_poly, den_A_poly, num_i_poly, den_i_poly, d_A, d_I =
        _mwc_state_polys(am, Set{Symbol}(metabolites(am.catalytic)), kept)
    # Each conformation's denominator times the other's free-enzyme weight
    # (`_law_branches`): every term is positive, so no face can cancel them.
    _assert_support(am, _law_branches(den_A_poly, den_i_poly, d_A, d_I))
    # A free-enzyme weight left for cross-weighting multiplies the other state's terms;
    # a normalized one renders as 1, which `_mwc_cross_weight` skips.
    D_A_expr = to_expr(d_A)
    D_I_expr = to_expr(d_I)

    N_A = to_expr(num_A_poly)
    Q_A = to_expr(den_A_poly)
    N_I = to_expr(num_i_poly)
    Q_I = to_expr(den_i_poly)

    # Each regulatory site's factor at its multiplicity, per conformation.
    reg_A = Any[_power_expr(_reg_site_expr(am, s, false), multiplicity(s)) for s in RS]
    reg_I = Any[_power_expr(_reg_site_expr(am, s, true), multiplicity(s)) for s in RS]
    # Numerator: N × Q_cat^(CatN-1) × the reg-site factors; denominator: Q_cat^CatN × them.
    num_term(N, Q, reg) =
        _nest_binary(:*, Any[N, (CatN > 1 ? (_power_expr(Q, CatN - 1),) : ())..., reg...])
    den_term(Q, reg) = _nest_binary(:*, Any[_power_expr(Q, CatN), reg...])

    num_A = _mwc_cross_weight(num_term(N_A, Q_A, reg_A), D_I_expr, CatN)
    den_A = _mwc_cross_weight(den_term(Q_A, reg_A), D_I_expr, CatN)
    den_I = _mwc_cross_weight(den_term(Q_I, reg_I), D_A_expr, CatN)
    full_den = _mwc_combine(den_A, den_I)

    if isempty(num_i_poly)
        # Native I-state numerator is zero: the steady-state fluxes of a broken cycle's
        # reachable-form-pruned I-graph cancel exactly, so the I-state cycle is dead and
        # the L*num_I term is dropped entirely (skip dead numerator branch). Q_I still
        # contributes to denominator as enzyme mass. A live redundant-path `:OnlyA`
        # mechanism (num_I ≠ 0) keeps the term.
        num_A, full_den
    else
        num_I = _mwc_cross_weight(num_term(N_I, Q_I, reg_I), D_A_expr, CatN)
        _mwc_combine(num_A, num_I), full_den
    end
end
