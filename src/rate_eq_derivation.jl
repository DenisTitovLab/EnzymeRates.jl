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
    dep, indep = _dependent_param_exprs(m isa AbstractEnzymeMechanism ? _concrete(m) : m)
    absent = [k for k in (indep..., :Keq) if !haskey(params, k)]
    isempty(absent) || error("all_parameters: params must hold every name of " *
                             "parameters(m) except E_total; it lacks " *
                             join(absent, ", "))
    names = sort!(collect(keys(dep)))
    merge(params, NamedTuple{Tuple(names)}(
        Tuple(_evaluate_dependent(dep[k], params) for k in names)))
end

"""The value at `params` of a dependent constant's expression (`build_power_expr`): a
power product of `Keq` and independent constants, or the number it is pinned to."""
_evaluate_dependent(x::Symbol, params) = params[x]
_evaluate_dependent(x::Real, params) = x
function _evaluate_dependent(e::Expr, params)
    op, args = e.args[1], [_evaluate_dependent(a, params) for a in e.args[2:end]]
    op === :* && return prod(args)
    op === :/ && return args[1] / args[2]
    op === :^ && return args[1]^args[2]
    error("_evaluate_dependent: no rule for the operator $op")
end

# ─── RE Group Helpers ───────────────────────────────────────

"""
Concentration symbols (substrates ∪ products ∪ regulators) that may appear in
the rate-equation polynomials — the same set `_poly_to_expr` treats as
concentrations and the only symbols `_reduce_conc_lowest_terms` is allowed to
shift. Parameter symbols are never in this set, so the concentration-GCD can
never drop a fitted parameter.
"""
_concentration_symbols(rxn::EnzymeReaction) = Set(_metabolite_names(rxn))

"""
Free enzyme of a rapid-equilibrium segment: the form with the fewest bound
metabolites, tie-broken toward no covalent residual, then a deterministic
name. Referencing each segment's alphas to this form makes the derivation
independent of the step order and yields the readable `1 + [S]/K + …` form.
"""
function _segment_root(segment, species)
    argmin(i -> (length(bound(species[i])),
                 has_residual(species[i]) ? 1 : 0,
                 string(name(species[i]))), segment)
end

"""
    _re_weight_ratio(s, K) -> POLY

w(to)/w(from) of rapid-equilibrium step `s` as a monomial: [M]/K for a binding of M,
plain or fused (K a dissociation constant), and K·Π[consumed]/Π[released] for every
other step (K the equilibrium constant of the stored direction, products over
reactants: [to]·Π[released] / ([from]·Π[consumed]), as its name `K_<from>_to_<to>`
says).
"""
function _re_weight_ratio(s::Step, K::Symbol)
    d = Dict{Symbol, Int}(K => is_binding(s) ? -1 : 1)
    for m in consumed(s); d[name(m)] = get(d, name(m), 0) + 1; end
    for m in released(s); d[name(m)] = get(d, name(m), 0) - 1; end
    filter!(p -> p.second != 0, d)
    POLY(sort!(MONO(collect(d)); by = first) => 1)
end

"""
Compute alpha factors (relative concentrations within RE segments) as POLY
values. Iterates the flat steps `flat` of a mechanism of reaction `rxn` directly; each
RE step's weight ratio comes from `_re_weight_ratio`, entering `to_species` from
`from_species` or the reverse. `idx[sp]` is the position of form `sp` in `species`;
`step_to_K[j]` is the parameter Symbol for the RE step at flat position `j`, and
`conc_set` the concentration symbols. Raises when the RE steps close a catalytic cycle,
which gives the mechanism no finite rate.
"""
function _compute_alpha(rxn, flat, conc_set, species, idx, segments, step_to_K)
    N = length(species)
    alpha = Vector{POLY}(fill(poly_one(), N))

    for segment in segments
        length(segment) == 1 && continue
        root = _segment_root(segment, species)
        visited = Set{Int}([root])
        queue = [root]
        while !isempty(queue)
            cur = popfirst!(queue)
            for (j, (s, _)) in enumerate(flat)
                is_equilibrium(s) || continue
                i_f = idx[from_species(s)]
                j_f = idx[to_species(s)]
                K = step_to_K[j]
                if i_f == cur && j_f ∉ visited
                    alpha[j_f] = poly_mul(alpha[cur], _re_weight_ratio(s, K))
                    push!(visited, j_f); push!(queue, j_f)
                elseif j_f == cur && i_f ∉ visited
                    alpha[i_f] = poly_mul(alpha[cur],
                                          _invert_monomial(_re_weight_ratio(s, K)))
                    push!(visited, i_f); push!(queue, i_f)
                end
            end
        end
    end

    # Every RE step must reproduce its own concentration ratio from the weights; a
    # mismatch means a cycle of RE steps performs turnover.
    conc(p) = Dict(k => v for (k, v) in only(keys(p)) if k in conc_set)
    for (j, (s, _)) in enumerate(flat)
        is_equilibrium(s) || continue
        a, b = idx[from_species(s)], idx[to_species(s)]
        ratio = poly_mul(alpha[b], _invert_monomial(alpha[a]))
        conc(ratio) == conc(_re_weight_ratio(s, step_to_K[j])) || error(
            "rate_equation: the rapid-equilibrium steps close a catalytic cycle " *
            "(through $(_form_label(from_species(s), rxn)) ⇌ " *
            "$(_form_label(to_species(s), rxn))), so the " *
            "mechanism has no finite rate. Make one step of the cycle steady-state.")
    end
    alpha
end

"""Build rate poly for one SS step direction in Laurent (fractional) form: the rate
constant `k_poly` times the metabolites `mets` it takes up times the weight of its source
form `i_form`."""
function _ss_contrib(k_poly, mets, i_form, alpha)
    r = isempty(mets) ? k_poly : poly_mul(k_poly, reduce(poly_mul, poly_sym.(name.(mets))))
    poly_mul(r, alpha[i_form])
end

# ─── Raw Rate Equation Derivation (Unified Cha / King-Altman) ───

"""
    _raw_symbolic_rate_polys(rxn, groups, names) -> (num, den, d_free)

Build raw numerator and denominator POLYs for the rate equation of the kinetic groups
`groups` of a mechanism of reaction `rxn`. `names[j]` holds the parameter Symbols of
flat step `j`'s constants (`[K]` or `[kf, kr]`), which are the leaves of `num`/`den`.
Aborts, before the King–Altman expansion, a graph whose denominator would be too large
to derive (`_assert_derivable`).

Also returns `d_free`, the weight in the returned denominator of the free
resting enzyme (the form with empty `bound` and empty `residual`): the
spanning-tree weight `D[g_free]` of its segment times the concentration monomial
that brings `num`/`den` to lowest terms. Free E roots its segment, so its own
`alpha` is 1.
"""
function _raw_symbolic_rate_polys(rxn::EnzymeReaction, groups::Vector{Vector{Step}},
                                  names::Vector{Vector{Symbol}})
    _assert_derivable(groups)
    species, segments, _, idx, seg = _re_segment_extras(groups)
    # A fully-inert conformation (every binding pruned) has no enumerated form; it
    # exists only as free enzyme — no flux, partition 1, D[g_free] 1.
    isempty(species) && return poly_zero(), poly_one(), poly_one()
    flat = _flat_steps(groups)
    conc_set = _concentration_symbols(rxn)
    step_to_K = Dict{Int, Symbol}(
        i => names[i][1] for i in eachindex(flat) if is_equilibrium(flat[i][1]))
    alpha = _compute_alpha(rxn, flat, conc_set, species, idx, segments, step_to_K)
    G = length(segments)

    # The Laplacian of the segment graph: each SS step joining two segments adds its
    # forward flux out of its source segment and its reverse flux out of its target. A
    # step within one segment enters only the numerator. `ss` keeps every SS step with
    # its two forms and fluxes for the numerator.
    L = [poly_zero() for _ in 1:G, _ in 1:G]
    ss = Tuple{Step, Int, Int, POLY, POLY}[]
    for (j, (s, _)) in enumerate(flat)
        is_equilibrium(s) && continue
        a, b = idx[from_species(s)], idx[to_species(s)]
        fwd = _ss_contrib(poly_sym(names[j][1]), consumed(s), a, alpha)
        rev = _ss_contrib(poly_sym(names[j][2]), released(s), b, alpha)
        push!(ss, (s, a, b, fwd, rev))
        g1, g2 = seg[a], seg[b]
        g1 == g2 && continue
        L[g1, g1] = poly_add(L[g1, g1], fwd); L[g1, g2] = poly_sub(L[g1, g2], fwd)
        L[g2, g2] = poly_add(L[g2, g2], rev); L[g2, g1] = poly_sub(L[g2, g1], rev)
    end
    # D[g], the spanning-tree weight of segment g: L without g's row and column.
    D = [(o = setdiff(1:G, g); sym_det(L[o, o])) for g in 1:G]

    # A graph with no unbound, residual-free form (e.g. a hand-written mechanism whose
    # every form carries a residual) has no free enzyme to weight, so `d_free` is 1.
    i_free = findfirst(f -> isempty(bound(f)) && isempty(residual(f)), species)
    d_free = i_free === nothing ? poly_one() : D[seg[i_free]]

    den = poly_zero()
    for g in 1:G
        sigma = reduce(poly_add, (alpha[i] for i in segments[g]); init=poly_zero())
        den = poly_add(den, poly_mul(sigma, D[g]))
    end

    # Numerator of the rate: v·den summed over steady-state steps. With u(f) the
    # first substrate's exponent in form f's RE weight, each SS step e contributes
    # ω_e·(forward − reverse flux) with ω_e = (copies of the first substrate e
    # consumes − copies it releases) + u(from_e) − u(to_e). Flux conservation at every
    # form makes the sum equal the net consumption of the first substrate for every
    # parameter value, so it needs no choice of reaction cut and does not depend on
    # the direction a step is written in. Steps with ω_e = 0 contribute nothing.
    x = name(first(substrates(rxn)))
    expo(p) = (mono = only(keys(p));
               k = findfirst(q -> q.first == x, mono);
               k === nothing ? 0 : mono[k].second)
    num = poly_zero()
    for (s, a, b, fwd, rev) in ss
        ω = count(m -> name(m) == x, consumed(s)) - count(m -> name(m) == x, released(s)) +
            expo(alpha[a]) - expo(alpha[b])
        ω == 0 && continue
        term = poly_sub(poly_mul(fwd, D[seg[a]]), poly_mul(rev, D[seg[b]]))
        num = poly_add(num, POLY(key => ω * c for (key, c) in term))
    end

    _reduce_conc_lowest_terms(num, den, d_free, conc_set)
end

_raw_symbolic_rate_polys(m::Mechanism) =
    _raw_symbolic_rate_polys(reaction(m), steps(m), _step_names(m, :None))
_raw_symbolic_rate_polys(@nospecialize(M::Type{<:EnzymeMechanism})) =
    _raw_symbolic_rate_polys(Mechanism(M()))

"""The kinetic groups of `m`'s law in conformation `state`: every group for `:None` and
`:A`; for `:I`, the inactive conformation's pruned graph (`_inactive_groups`)."""
_state_groups(m::Mechanism, ::Symbol) = steps(m)
_state_groups(am::AllostericMechanism, state::Symbol) =
    state === :I ? _inactive_groups(am) : steps(am)

"""
The parameter Symbols of each flat step's constants in conformation `state` of `m` —
`[K]` at rapid equilibrium, `[kf, kr]` at steady state — through the naming chokepoint
(`name(p, m)`), over the steps of `_state_groups(m, state)`: `:None` for a plain
mechanism; `:A` or `:I` for an allosteric one, an `:EqualAI` group's constants taking its
one shared name in both.
"""
function _step_names(m::_AnyMech, state::Symbol)
    tag(g) = m isa AllostericMechanism && cat_allo_state(m, g) === :EqualAI ?
             :EqualAI : state
    [Symbol[name(p, m) for p in _step_constants(s, tag(g))]
     for (s, g) in _flat_steps(_state_groups(m, state))]
end

"""
Estimate a mechanism's King–Altman denominator term count, V×τ, from its
catalytic segment graph WITHOUT deriving the rate equation.

`V` is the number of RE-connected segments (`_re_segment_extras`); `τ` is the spanning-tree
count of the segment graph, whose edges are the SS steps: by the Matrix–Tree theorem, the
determinant of its Laplacian with the first row and column dropped, computed exactly (`det`
of a `BigInt` matrix is fraction-free Bareiss; a disconnected graph gives 0). `V×τ` equals
the number of spanning-tree products in the denominator — the products the compiled equation
evaluates on every call — so it bounds fit cost. It guards the derivation itself
(`_assert_derivable` aborts when V×τ exceeds `MAX_RATE_EQUATION_TERMS`, before the O(G!)
symbolic cofactor expansion) and backs the `eq_complexity_filter` keyword of
[`identify_rate_equation`](@ref), which skips over-complex mechanisms before fitting. The
method for a mechanism counts an allosteric mechanism on `steps(am)`, its active
conformation.
"""
function _eq_complexity(groups::Vector{Vector{Step}})
    _, segments, _, idx, seg = _re_segment_extras(groups)
    G = length(segments)
    G <= 1 && return big(G)                          # single RE segment ⇒ τ = 1
    lap = zeros(BigInt, G, G)
    for group in groups, s in group
        is_equilibrium(s) && continue                # SS steps are the segment-graph edges
        g1, g2 = seg[idx[from_species(s)]], seg[idx[to_species(s)]]
        g1 == g2 && continue
        lap[g1, g1] += 1; lap[g2, g2] += 1; lap[g1, g2] -= 1; lap[g2, g1] -= 1
    end
    G * det(lap[2:end, 2:end])
end
_eq_complexity(m::_AnyMech) = _eq_complexity(steps(m))

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
function _assert_derivable(groups::Vector{Vector{Step}})
    _eq_complexity(groups) > MAX_RATE_EQUATION_TERMS && throw(RateEquationTooLargeError())
    return nothing
end

# ─── Expr generation from POLY ──────────────────────────────

"""The numerator and denominator Exprs of the rate equation of `M`, from its raw
symbolic rate polys: flat sums for the printed equation, or in Horner form
(`_poly_to_horner`) for the generated code when `horner` is true."""
function _num_den_exprs(@nospecialize(M::Type{<:EnzymeMechanism}); horner = false)
    m = Mechanism(M())
    names = _step_names(m, :None)
    num, den, _ = _raw_symbolic_rate_polys(reaction(m), steps(m), names)
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
Error if a metabolite of `M` has the name of a parameter in `param_syms` or `dep`. The
rate equation binds the parameters, the concentrations and the dependent parameters as
local names of one scope, so a shared name would let one overwrite another. Every
parameter name holds `_` or is `Keq` or `L`, and the reaction and mechanism constructors
reject those metabolite names (`_check_metabolite_names`), so this check cannot fire for
a mechanism the constructors build.
"""
function _assert_no_name_clash(@nospecialize(M::Type{<:AbstractEnzymeMechanism}),
                               param_syms, dep)
    for m in metabolites(M())
        (m in param_syms || haskey(dep, m)) && error(
            "Metabolite $m has the name of a parameter of this mechanism's rate " *
            "equation; rename the metabolite")
    end
end

"""
The body of the generated `rate_equation` of `M`: destructure `param_syms` from
`params` and the metabolites from `concs`, assign the dependent parameters `dep`, and
return `E_total * (num) / (den)` with `num` and `den` from `_num_den_exprs` in Horner
form. It first runs `_assert_no_name_clash`, which cannot fire for a mechanism the
constructors build.
"""
function _rate_body(@nospecialize(M::Type{<:AbstractEnzymeMechanism}), param_syms, dep)
    _assert_no_name_clash(M, param_syms, dep)
    num, den = _num_den_exprs(M; horner = true)
    Expr(:block,
        _destructuring_expr(param_syms, :params),
        _destructuring_expr(metabolites(M()), :concs),
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
    # holds `_`, which no metabolite name does.
    dep, indep = _dependent_param_exprs(M_type)
    _rate_body(M_type, (indep..., :Keq, :E_total), dep)
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
The `rate_equation_string` text of `M`: the `params` destructure of `param_syms`, the
`concs` destructure, one line `sym = rhs` per dependent parameter of `dep`, and the
`v = E_total * (num) / (den)` line from `_num_den_exprs`. A dependent whose right-hand
side mentions `Keq` goes under `# Haldane constraints:`, every other one under
`# Wegscheider constraints:`. Each section is sorted by name, which is load-bearing: the
eq_hash dedup of rate-equivalent mechanisms compares these strings, so their line order
must not depend on the order the solve emits its dependents in. It first runs
`_assert_no_name_clash`, which cannot fire for a mechanism the constructors build.
"""
function _equation_text(@nospecialize(M::Type{<:AbstractEnzymeMechanism}), param_syms,
                        dep)
    _assert_no_name_clash(M, param_syms, dep)
    weg, hal = String[], String[]
    for (sym, rhs) in sort!(collect(dep); by = first)
        push!(_mentions(rhs, :Keq) ? hal : weg, "$sym = $(_expr_to_string(rhs))")
    end
    num, den = _num_den_exprs(M)
    lines = ["(; $(join(param_syms, ", "))) = params",
             "(; $(join(metabolites(M()), ", "))) = concs"]
    isempty(weg) || push!(lines, "# Wegscheider constraints:", weg...)
    isempty(hal) || push!(lines, "# Haldane constraints:", hal...)
    push!(lines, "v = E_total * ($(_expr_to_string(num))) / ($(_expr_to_string(den)))")
    join(lines, "\n")
end

function rate_equation_string(@nospecialize(m::AbstractEnzymeMechanism))
    dep, indep = _dependent_param_exprs(typeof(m))
    _equation_text(typeof(m), (indep..., :Keq, :E_total), dep)
end

# ─── Steady-state rate-constant names ──────────────────────────

"""
The names of the steady-state rate constants of `em`'s law (`_step_names`): for an
`AllostericEnzymeMechanism`, those of both conformations, the inactive one's `I_` names
included. `_scale_rate_constants` scales exactly these, leaving the rapid-equilibrium
constants, `Keq`, `E_total`, `L` and the regulator constants alone.
"""
function _ss_rate_constant_names(@nospecialize(em::AbstractEnzymeMechanism))
    m = _concrete(em)
    states = m isa Mechanism ? (:None,) : (:A, :I)
    Set{Symbol}(n for st in states
                for ((s, _), ns) in zip(_flat_steps(_state_groups(m, st)),
                                        _step_names(m, st))
                if !is_equilibrium(s) for n in ns)
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
Names of enzyme forms in the connected component of the free enzyme over ALL
steps of `groups` (rapid-equilibrium and steady-state alike). Every form
carrying neither a bound metabolite nor a residual seeds the search: such a form
interconverts between conformations under formulation 1, and so is a root the
inactive conformation can be entered through. Enumeration yields exactly one
such form; the hand-written DSL admits a second one in rapid equilibrium with the
first (`_assert_one_flip_segment`), and every match then seeds. A ping-pong covalent intermediate carries a residual and is
therefore not a root: a component the free enzyme cannot reach holds no inactive
mass, and leaving it in place would strand the free-enzyme spanning tree
(`D[g_free] = 0`).
"""
function _reachable_from_free(groups)
    reach = Set{Symbol}(name(f) for grp in groups for s in grp
                        for f in (from_species(s), to_species(s))
                        if isempty(bound(f)) && isempty(residual(f)))
    changed = true
    while changed
        changed = false
        for grp in groups, s in grp
            fn, tn = name(from_species(s)), name(to_species(s))
            if fn in reach && tn ∉ reach
                push!(reach, tn); changed = true
            elseif tn in reach && fn ∉ reach
                push!(reach, fn); changed = true
            end
        end
    end
    reach
end

"""
The inactive conformation's step graph of `am`, aligned with `steps(am)`: an `:OnlyA`
group is empty, and every other group keeps the steps whose two forms are reachable
from the free enzyme over the non-`:OnlyA` groups (`_reachable_from_free`). The steps
are `am`'s own, in their stored direction, so the inactive law reads `am`'s orientation
and names. The derivation (`_state_rate_polys`), the joint solve (`_inactive_law`) and
the copy rule (`_redundant_copy_groups`) read the inactive state from it.
"""
function _inactive_groups(am::AllostericMechanism)
    tags = cat_allo_states(am)
    reach = _reachable_from_free(
        [group for (g, group) in enumerate(steps(am)) if tags[g] !== :OnlyA])
    [tags[g] === :OnlyA ? Step[] :
     Step[s for s in group if name(from_species(s)) in reach &&
                              name(to_species(s)) in reach]
     for (g, group) in enumerate(steps(am))]
end

"""
    _state_rate_polys(am, state) -> (num, den, d_free)

The rate polynomials of `am`'s catalytic graph in conformation `state` (`:A` or `:I`),
in that state's names (`_step_names`): the active conformation on the full graph, the
inactive one on its pruned graph (`_inactive_groups`), each step in `am`'s stored
direction. King–Altman on the pruned graph derives the inactive conformation's binding
partition, and for a dead cycle its steady-state fluxes cancel, so `num` is 0. The
inactive graph is checked for a rapid-equilibrium segment with no bottom form
(`_assert_bottomed_re_segments`), as the mechanism constructors check the full graph.
"""
function _state_rate_polys(am::AllostericMechanism, state::Symbol)
    groups = _state_groups(am, state)
    state === :I && _assert_bottomed_re_segments(reaction(am), groups)
    _raw_symbolic_rate_polys(reaction(am), groups, _step_names(am, state))
end

"""
The regulator constants of `am` listed with conformation `state`, site by site: for `:A`
every ligand the active conformation binds, an `:EqualAI` ligand under its one bare name
(`K_<ligand>reg`, state `:EqualAI`) and an `:OnlyA` or `:NonequalAI` one under its
`K_A_` name; for `:I` the `K_I_` constant of every `:OnlyI` or `:NonequalAI` ligand.
"""
_kreg_params(am::AllostericMechanism, state::Symbol) =
    state === :A ?
    Kreg[Kreg(site, lig, tag === :EqualAI ? :EqualAI : :A)
         for site in regulatory_sites(am)
         for (lig, tag) in zip(ligands(site), allo_states(site)) if tag !== :OnlyI] :
    Kreg[Kreg(site, lig, :I) for site in regulatory_sites(am)
         for (lig, tag) in zip(ligands(site), allo_states(site))
         if tag in (:OnlyI, :NonequalAI)]

# ─── Rate body building helpers ───────────────────────────────────

"""
The dissociation-constant name of ligand `lig` (allosteric tag `tag`) at regulatory
`site` of `am` in the active (`inactive = false`) or inactive conformation, or
`nothing` when the ligand is absent from that conformation (`:OnlyI` in the active
state, `:OnlyA` in the inactive one). An `:EqualAI` ligand has one bare name in both
conformations (`K_<ligand>reg`); an `:OnlyA`, `:OnlyI` or `:NonequalAI` ligand has its
conformation's `K_A_` or `K_I_` name. Renders the name via the `name(::Kreg, am)`
chokepoint.
"""
function _reg_K(am::AllostericMechanism, site, lig, tag::Symbol, inactive::Bool)
    tag === (inactive ? :OnlyA : :OnlyI) && return nothing
    name(Kreg(site, lig, tag === :EqualAI ? :EqualAI : inactive ? :I : :A), am)
end

"""Build the regulatory site partition function expression: 1 + lig/K_lig_reg_i + ...
over the ligands of `site` present in the given conformation (`_reg_K`)."""
function _reg_site_expr(am::AllostericMechanism, site, inactive::Bool)
    terms = Any[1]
    for (lig, tag) in zip(ligands(site), allo_states(site))
        K = _reg_K(am, site, lig, tag, inactive)
        K === nothing || push!(terms, :($(name(lig)) / $K))
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
    _mwc_state_polys(am, mets) -> (num_A, den_A, num_I, den_I, d_A, d_I)

The catalytic rate polys of both conformations of `am` (`_state_rate_polys`) under the
formulation-1 per-state free-enzyme normalization that `rate_equation`
(`_num_den_exprs`) uses. The I-state polys come from the inactive conformation's pruned
graph (`_inactive_groups` drops `:OnlyA` groups and every form they disconnect from free
E). King–Altman on that subgraph derives the inactive conformation's binding partition,
and for a dead cycle the pruned graph's steady-state fluxes cancel exactly, so `num_I`
is 0.

The free-enzyme weights `d_A` and `d_I` combine three ways, all rendering the same
value (`mets` are the catalytic metabolites):
- `d_A == d_I`: the polys stay raw (identical conformations; the factor cancels);
- both metabolite-free monomials: each state's polys are divided by its own weight
  (the clean standard-MWC form);
- otherwise: the polys stay raw, and each state is cross-weighted by the other state's
  weight (`d^n` on the rate law's terms).
Only the third case returns the weights; the first two return `poly_one()` for both.
"""
function _mwc_state_polys(am::AllostericMechanism, mets::Set{Symbol})
    num_A, den_A, d_A = _state_rate_polys(am, :A)
    num_I, den_I, d_I = _state_rate_polys(am, :I)
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
Assemble the MWC numerator and denominator Exprs.
Returns `(full_num, full_den)`, written with flat sums for the printed equation, or in
Horner form (`_poly_to_horner`) for the generated code when `horner` is true.
Per-active-site normalization: the
numerator carries no leading `catalytic_multiplicity` factor; only the
`Q_cat^(CatN-1)` / `Q_cat^CatN` binding-statistics powers remain.
"""
function _num_den_exprs(@nospecialize(M_type::Type{<:AllostericEnzymeMechanism});
                        horner = false)
    m = M_type()
    am = AllostericMechanism(m)
    CatN = catalytic_multiplicity(m)
    RS = regulatory_sites(am)

    # The active conformation's catalytic names drive `_poly_to_expr`'s
    # param/metabolite ordering split; the I-poly's `I_` names sort as non-params.
    cat_params = Set{Symbol}(n for ns in _step_names(am, :A) for n in ns)
    to_expr(p) = horner ? _poly_to_horner(p) : _poly_to_expr(p, cat_params)
    num_A_poly, den_A_poly, num_i_poly, den_i_poly, d_A, d_I =
        _mwc_state_polys(am, Set{Symbol}(metabolites(catalytic_mechanism(m))))
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
