# ABOUTME: King-Altman/Cha rate equation derivation via @generated functions.
# ABOUTME: Parameters API, kcat computation, and MWC allosteric rate assembly.

# ─── Parameters API ─────────────────────────────────────────

# Every mechanism is a singleton type of its own, so a method specialized on one
# compiles again for each mechanism. The derivation helpers that lift a singleton,
# or its type, to its concrete form therefore take it `@nospecialize`.

"""
Suffix appended to single-symbol equality lines whose LHS got folded
into the kinetic-group rename map. Both display sites (User defined
kinetic-group merges and absorbed single-symbol Wegscheider ties) emit
this exact string so the rate-equation dedup key in
`identify_rate_equation.jl` can strip these provenance lines.
"""
const ANNOTATION_SUBSTITUTED = "  (substituted into v)"

"""
    parameters(m::EnzymeMechanism, [mode])
    parameters(m::AllostericEnzymeMechanism, [mode])

Return the parameter names required for the given mode as a tuple of Symbols.

# Modes
- `Reduced` (default): independent k's + Keq + E_total. The set of
  symbols the user supplies to evaluate the Haldane-reduced rate
  equation. Returned for both `EnzymeMechanism` and
  `AllostericEnzymeMechanism`.
- `Full`: all raw rate-constant symbols + `E_total`, that is "all 2N k's +
  `E_total`." Full mode is for plain `EnzymeMechanism`s only; an
  `AllostericEnzymeMechanism` has no Full rate equation and no Full parameter list.
"""
function parameters end

parameters(m::Union{AbstractEnzymeMechanism, Mechanism, AllostericMechanism}) =
    parameters(m, Reduced)

@generated function parameters(
    ::EnzymeMechanism{Sig}, ::FullMode,
) where {Sig}
    mech = Mechanism(EnzymeMechanism{Sig}())
    params = _enumerate_parameters_full(mech)
    names = Tuple(name(p, mech) for p in params)
    Tuple((names..., :E_total))
end

@generated function parameters(::M, ::ReducedMode) where {M <: AbstractEnzymeMechanism}
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
# `EnzymeMechanism` / `AllostericEnzymeMechanism` for the 0-alloc/<100ns
# `rate_equation` guarantee.
parameters(m::Union{Mechanism, AllostericMechanism},
           mode::AbstractRateEquationMode) =
    parameters(compile_mechanism(m), mode)
fitted_params(m::Union{Mechanism, AllostericMechanism}) =
    fitted_params(compile_mechanism(m))

"""
Build a renaming map for single-symbol Wegscheider RE ties between two
binding K's. Calls `_dependent_param_exprs_kernel` to discover
binding-K-to-binding-K Wegscheider closures of the form `K_a = K_b`
(RHS is a bare Symbol). Both sides must be binding K's (RE step with
metabolite on LHS) — absorbing a binding-K-to-iso-K tie would produce
inconsistent sign-flips when the kernel runs with the full rename,
since the binding-K column is sign-flipped (Kd convention) but the
iso-K column is not.

The rename means the polynomial in `v` uses the representative symbol
directly, so Source-C duplicates (split kinetic groups that
Wegscheider ties back together) collapse at hash time.
"""
function _build_wegscheider_rename_map(mech::Mechanism)
    rename = Dict{Symbol, Symbol}()
    step_params = _step_parameters(mech)
    # binding-K set: value-context rep name of each RE binding step. Walk
    # Mechanism.steps directly — an RE step that is a binding (`is_binding`,
    # plain or fused) is a binding step; step_params is indexed in the same
    # flat order.
    binding_set = Set{Symbol}()
    for (idx, (s, _)) in enumerate(_flat_steps(mech))
        is_equilibrium(s) && is_binding(s) || continue
        push!(binding_set, name(step_params[idx][1], mech))
    end
    # Pass 2: single-symbol Wegscheider RE ties between two binding K's.
    dep_raw, _ = _dependent_param_exprs_kernel(mech, rename)
    for (lhs, rhs) in dep_raw
        rhs isa Symbol || continue
        lhs in binding_set && rhs in binding_set || continue
        target = get(rename, rhs, rhs)
        rename[lhs] = target
        for k in collect(keys(rename))
            rename[k] == lhs && (rename[k] = target)
        end
    end
    rename
end

_build_wegscheider_rename_map(@nospecialize(M::Type{<:EnzymeMechanism})) =
    _build_wegscheider_rename_map(Mechanism(M()))
_build_wegscheider_rename_map(@nospecialize(m::EnzymeMechanism)) =
    _build_wegscheider_rename_map(typeof(m))

# ─── RE Group Helpers ───────────────────────────────────────

"""
Concentration symbols (substrates ∪ products ∪ regulators) that may appear in
the rate-equation polynomials — the same set `_poly_to_expr` treats as
concentrations and the only symbols `_reduce_conc_lowest_terms` is allowed to
shift. Parameter symbols are never in this set, so the concentration-GCD can
never drop a fitted parameter.
"""
_concentration_symbols(mech::Mechanism) = Set(_metabolite_names(reaction(mech)))

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
values. Iterates `mech.steps` directly; each RE step's weight ratio comes
from `_re_weight_ratio`, entering `to_species` from `from_species` or the
reverse. `idx[sp]` is the position of form `sp` in `species`; `step_to_K[j]` is
the parameter Symbol for the RE step at flat position `j` (rep-renamed via the
`name(p, m)` chokepoint). Raises when the RE steps close a catalytic cycle,
which gives the mechanism no finite rate.
"""
function _compute_alpha(mech::Mechanism, species, idx, segments, step_to_K)
    N = length(species)
    alpha = Vector{POLY}(fill(poly_one(), N))
    flat = _flat_steps(mech)

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
    conc_set = _concentration_symbols(mech)
    conc(p) = Dict(k => v for (k, v) in only(keys(p)) if k in conc_set)
    for (j, (s, _)) in enumerate(flat)
        is_equilibrium(s) || continue
        a, b = idx[from_species(s)], idx[to_species(s)]
        ratio = poly_mul(alpha[b], _invert_monomial(alpha[a]))
        conc(ratio) == conc(_re_weight_ratio(s, step_to_K[j])) || error(
            "rate_equation: the rapid-equilibrium steps close a catalytic cycle " *
            "(through $(name(from_species(s))) ⇌ $(name(to_species(s)))), so the " *
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
Build raw numerator and denominator POLYs for the rate equation by
walking the lifted `Mechanism`. Parameter Symbols on the leaves of
`num`/`den` are produced via the `name(p, mech)` chokepoint (which
collapses kinetic-group members to their rep's name). `rename_map` then
applies any single-symbol Wegscheider ties as a post-pass.

Also returns `d_free`, the weight in the returned denominator of the free
resting enzyme (the form with empty `bound` and empty `residual`): the
spanning-tree weight `D[g_free]` of its segment times the concentration monomial
that brings `num`/`den` to lowest terms. Free E roots its segment, so its own
`alpha` is 1.
"""
function _raw_symbolic_rate_polys(mech::Mechanism, step_params, rename_map)
    species, segments, _, idx, seg = _re_segment_extras(steps(mech))
    # A fully-inert conformation (every binding pruned) has no enumerated form; it
    # exists only as free enzyme — no flux, partition 1, D[g_free] 1.
    isempty(species) && return poly_zero(), poly_one(), poly_one()
    flat = _flat_steps(mech)
    step_to_K = Dict{Int, Symbol}(
        i => name(step_params[i][1], mech)
        for i in eachindex(flat) if is_equilibrium(flat[i][1]))
    alpha = _compute_alpha(mech, species, idx, segments, step_to_K)
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
        fwd = _ss_contrib(poly_sym(name(step_params[j][1], mech)), consumed(s), a, alpha)
        rev = _ss_contrib(poly_sym(name(step_params[j][2], mech)), released(s), b, alpha)
        push!(ss, (s, a, b, fwd, rev))
        g1, g2 = seg[a], seg[b]
        g1 == g2 && continue
        L[g1, g1] = poly_add(L[g1, g1], fwd); L[g1, g2] = poly_sub(L[g1, g2], fwd)
        L[g2, g2] = poly_add(L[g2, g2], rev); L[g2, g1] = poly_sub(L[g2, g1], rev)
    end
    # D[g], the spanning-tree weight of segment g: L without g's row and column.
    D = [(o = setdiff(1:G, g); sym_det(L[o, o])) for g in 1:G]

    # A fully-inert conformation (every binding pruned, e.g. all-`:OnlyA` in the
    # inactive state) has no reactions and so no enumerated form; its free enzyme
    # spans the whole (empty) graph, so `D[g_free] = 1`.
    i_free = findfirst(f -> isempty(bound(f)) && isempty(residual(f)), species)
    d_free = i_free === nothing ? poly_one() :
             _rename_symbols(D[seg[i_free]], rename_map)

    den = poly_zero()
    for g in 1:G
        sigma = reduce(poly_add, (alpha[i] for i in segments[g]); init=poly_zero())
        csigma = _rename_symbols(sigma, rename_map)
        den = poly_add(den, poly_mul(csigma, D[g]))
    end

    # Numerator of the rate: v·den summed over steady-state steps. With u(f) the
    # first substrate's exponent in form f's RE weight, each SS step e contributes
    # ω_e·(forward − reverse flux) with ω_e = (copies of the first substrate e
    # consumes − copies it releases) + u(from_e) − u(to_e). Flux conservation at every
    # form makes the sum equal the net consumption of the first substrate for every
    # parameter value, so it needs no choice of reaction cut and does not depend on
    # the direction a step is written in. Steps with ω_e = 0 contribute nothing.
    x = name(first(substrates(reaction(mech))))
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

    num = _rename_symbols(num, rename_map)
    den = _rename_symbols(den, rename_map)
    conc_set = _concentration_symbols(mech)
    _reduce_conc_lowest_terms(num, den, d_free, conc_set)
end

function _raw_symbolic_rate_polys(@nospecialize(M::Type{<:EnzymeMechanism}))
    mech = Mechanism(M())
    _assert_derivable(mech)
    step_params = _step_parameters(mech)
    rename_map = _build_wegscheider_rename_map(M)
    _raw_symbolic_rate_polys(mech, step_params, rename_map)
end

"""
Estimate a mechanism's King–Altman denominator term count, V×τ, from its
catalytic segment graph WITHOUT deriving the rate equation.

`V` is the number of RE-connected segments (`_re_segment_extras`); `τ` is the
spanning-tree count of the segment graph, whose edges are the SS steps: by the
Matrix–Tree theorem, the determinant of its Laplacian with the first row and column
dropped, computed exactly (`det` of a `BigInt` matrix is fraction-free Bareiss; a
disconnected graph gives 0). `V×τ` equals the number of spanning-tree
products in the denominator — the products the compiled equation evaluates on
every call — so it bounds fit cost. It guards the derivation itself
(`_assert_derivable` aborts when V×τ exceeds `MAX_RATE_EQUATION_TERMS`, before the
O(G!) symbolic cofactor expansion) and backs the `eq_complexity_filter` keyword
of [`identify_rate_equation`](@ref), which skips over-complex mechanisms before
fitting.
"""
function _eq_complexity(mech::Mechanism)
    _, segments, _, idx, seg = _re_segment_extras(steps(mech))
    G = length(segments)
    G <= 1 && return big(G)                          # single RE segment ⇒ τ = 1
    lap = zeros(BigInt, G, G)
    for group in steps(mech), s in group
        is_equilibrium(s) && continue                # SS steps are the segment-graph edges
        g1, g2 = seg[idx[from_species(s)]], seg[idx[to_species(s)]]
        g1 == g2 && continue
        lap[g1, g1] += 1; lap[g2, g2] += 1; lap[g1, g2] -= 1; lap[g2, g1] -= 1
    end
    G * det(lap[2:end, 2:end])
end
_eq_complexity(m::AllostericMechanism) = _eq_complexity(_state_mechanism(m, :A))

"""
Abort the rate-equation derivation for a mechanism whose denominator term count
(V×τ) exceeds `MAX_RATE_EQUATION_TERMS`. Called at each num/den derivation entry
point, before the O(G!) symbolic cofactor expansion in `sym_det` — a numeric
O(G³) check that errors quickly instead of hanging on the factorial expansion.
"""
function _assert_derivable(mech::Mechanism)
    _eq_complexity(mech) > MAX_RATE_EQUATION_TERMS && error(
        "Rate equation for this mechanism has more than $MAX_RATE_EQUATION_TERMS " *
        "polynomial terms (limit: $MAX_RATE_EQUATION_TERMS). Equations this large " *
        "take a very long time to compile and are unlikely to be practically " *
        "useful for parameter fitting.")
    return nothing
end

# ─── Expr generation from POLY ──────────────────────────────

"""
Compute the raw rate expression (bare symbols) and sorted parameter/concentration symbols.
Returns `(expr, all_params, sorted_concs)`.
"""
function _raw_rate_expr_and_symbols(@nospecialize(M::Type{<:EnzymeMechanism}))
    num, den, _ = _raw_symbolic_rate_polys(M)
    m = M()
    param_syms = Set{Symbol}(_raw_param_symbols(m))
    num_expr = _poly_to_expr(num, param_syms)
    den_expr = _poly_to_expr(den, param_syms)
    expr = :(E_total * ($num_expr) / ($den_expr))
    all_params = _sorted_raw_param_symbols(M)
    return expr, all_params, metabolites(m)
end

# ─── Mode-dispatched rate_equation ────────────────────────────

"""
    rate_equation(m, concs, params, [mode])

Return the net reaction rate of mechanism `m` at the metabolite concentrations
`concs` and parameters `params` (both `NamedTuple`s). The rate equation itself
is derived from the mechanism's elementary steps by the King–Altman/Cha method:
rapid-equilibrium steps collapse to binding constants, steady-state steps are
assembled into the King–Altman determinant, and the two combine into the
quasi-steady-state flux through the whole mechanism — so one call returns the
overall turnover, not the rate of a single step.

`mode` selects how parameters enter the equation. The default [`Reduced`](@ref)
applies the Haldane/Wegscheider reduction, deriving the dependent rate constants
from `Keq` and the independent parameters; [`Full`](@ref) instead takes every
rate constant as independent.

The derivation runs once, at compile time: the body is generated as a single
arithmetic expression with no allocations, loops, or matrix operations. Use
[`rate_equation_string`](@ref) to inspect the derived equation symbolically.
"""
function rate_equation end

rate_equation(m::Union{AbstractEnzymeMechanism, Mechanism, AllostericMechanism},
              concs, params) = rate_equation(m, concs, params, Reduced)

# Concrete-mechanism convenience: lift to the singleton (see the note on
# the parameters/fitted_params methods above — allocates, not the hot path).
rate_equation(m::Union{Mechanism, AllostericMechanism}, concs, params,
              mode::AbstractRateEquationMode) =
    rate_equation(compile_mechanism(m), concs, params, mode)

@generated function rate_equation(
    m::M, concs::NamedTuple, params::NamedTuple, ::FullMode,
) where {M <: EnzymeMechanism}
    _build_rate_body(M, FullMode)
end

@generated function rate_equation(
    m::M, concs::NamedTuple, params::NamedTuple,
    ::ReducedMode,
) where {M <: EnzymeMechanism}
    _build_rate_body(M, ReducedMode)
end

# ─── String Representation ────────────────────────────────────

"""
    rate_equation_string(m, [mode]) -> String

Return the symbolic rate equation for mechanism `m` as a multi-line
`String` (it returns the text — it does not print). `mode` is `Reduced`
(default) or `Full`; pass a concrete `Mechanism` / `AllostericMechanism`
or its compiled [`EnzymeMechanism`](@ref) singleton.

The string is a runnable transcript of how [`rate_equation`](@ref)
evaluates: a `(; …) = params` destructure line, a `(; …) = concs`
destructure line, then the `v = E_total * (num) / (den)` line. In
`Reduced` mode, dependent rate constants are listed first under
`# Wegscheider constraints:` and `# Haldane constraints:` headers — the
thermodynamic identities that eliminate parameters — and only the
independent set appears in the `params` destructure. In `Full` mode every
rate constant is independent, so there is no constraint section. `Full`
mode is defined for `EnzymeMechanism` only; an `AllostericEnzymeMechanism`
supports `Reduced` mode only.

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

rate_equation_string(m::Union{AbstractEnzymeMechanism, Mechanism, AllostericMechanism}) =
    rate_equation_string(m, Reduced)

# Concrete-mechanism convenience: lift to the singleton.
rate_equation_string(m::Union{Mechanism, AllostericMechanism},
                     mode::AbstractRateEquationMode) =
    rate_equation_string(compile_mechanism(m), mode)

"""Build the `v = E_total * (num) / (den)` line from the raw symbolic rate polys."""
function _rate_v_line(@nospecialize(M::Type{<:EnzymeMechanism}))
    num, den, _ = _raw_symbolic_rate_polys(M)
    m = M()
    ps = Set{Symbol}(_raw_param_symbols(m))
    "v = E_total * ($(_expr_to_string(_poly_to_expr(num, ps)))) / " *
        "($(_expr_to_string(_poly_to_expr(den, ps))))"
end

function rate_equation_string(@nospecialize(em::EnzymeMechanism), ::FullMode)
    mech = Mechanism(em)
    param_names = Symbol[name(p, mech) for p in _enumerate_parameters_full(mech)]
    lines = ["(; $(join((param_names..., :E_total), ", "))) = params",
             "(; $(join(metabolites(em), ", "))) = concs"]
    push!(lines, _rate_v_line(typeof(em)))
    join(lines, "\n")
end

"""
Render each dep-map entry as a constraint line and append it to `weg_lines` or
`hal_lines`, split by whether its RHS references `Keq` (Haldane) or not
(Wegscheider). Single-symbol RHSes get the substituted-into-v annotation;
multi-symbol RHSes get runtime assignment in `_build_rate_body` (no annotation).
Entries are visited in lexicographic LHS order.
"""
function _partition_constraint_lines!(weg_lines, hal_lines, dep)
    for (sym, expr) in sort(collect(dep); by=p -> string(p[1]))
        is_haldane = _mentions(expr, :Keq)
        suffix = expr isa Symbol ? ANNOTATION_SUBSTITUTED : ""
        push!(is_haldane ? hal_lines : weg_lines, "$sym = $(string(expr))$suffix")
    end
end

"""Append the `# Wegscheider constraints:` and `# Haldane constraints:` sections
(each skipped when empty) to `lines`."""
function _append_constraint_sections!(lines, weg_lines, hal_lines)
    isempty(weg_lines)  ||
        (push!(lines, "# Wegscheider constraints:");  append!(lines, weg_lines))
    isempty(hal_lines)  ||
        (push!(lines, "# Haldane constraints:");      append!(lines, hal_lines))
end

function rate_equation_string(@nospecialize(m::EnzymeMechanism), ::ReducedMode)
    M = typeof(m)
    _, indep = _dependent_param_exprs(M)

    dep_raw, _ = _dependent_param_exprs_kernel(M, Dict{Symbol, Symbol}())
    weg_lines, hal_lines = String[], String[]
    _partition_constraint_lines!(weg_lines, hal_lines, dep_raw)

    lines = ["(; $(join((indep..., :Keq, :E_total), ", "))) = params",
             "(; $(join(metabolites(m), ", "))) = concs"]
    _append_constraint_sections!(lines, weg_lines, hal_lines)
    push!(lines, _rate_v_line(M))
    join(lines, "\n")
end

# ─── kcat Computation Helpers ──────────────────────────────────

"""
Set of Symbol names for SS rate-constant parameters (Kon, Koff, Kfor,
Krev) of `em`. For `AllostericEnzymeMechanism`, also includes the I-state
names (`I_` tag after the prefix, e.g. `k_I_ES_to_EP`) of every SS rate
constant that lives in the inactive state polynomial. Routes Symbol production through the
`name(p, m)` chokepoint via Parameter-subtype dispatch. Used by
`rescale_parameter_values` to scale only SS k's without touching RE
Kd's, Keq, E_total, L, or regulatory K's.
"""
function _ss_rate_constant_names(@nospecialize(em::AbstractEnzymeMechanism))
    m = _concrete(em)
    params = m isa Mechanism ? _enumerate_parameters_full(m) :
        [_cat_params(m, :A); _cat_params(m, :I)]
    Set{Symbol}(name(p, m) for p in params if p isa Union{Kon, Koff, Kfor, Krev})
end

"""Group `num` and `den` POLYs by metabolite monomial pattern, using
`k_param_names` to classify each symbol. Returns `(num_groups,
den_groups)` where each value is a POLY of k-monomials sharing the
same met-monomial. Reverse (negative-coefficient) terms are dropped
from the numerator. Used by `_kcat_forward` to compare saturating
metabolite patterns across A/I states."""
function _kcat_groups_from_polys(num::POLY, den::POLY,
                                  k_param_names::Set{Symbol})
    function split_mono(mono::MONO)
        k_mono = MONO()
        met_mono = MONO()
        for (s, e) in mono
            if s in k_param_names || s == :Keq
                push!(k_mono, s => e)
            elseif s != :E_total
                push!(met_mono, s => e)
            end
        end
        sort!(k_mono; by=first), sort!(met_mono; by=first)
    end

    num_groups = Dict{MONO, POLY}()
    for (mono, coeff) in num
        coeff > 0 || continue
        k_part, met_part = split_mono(mono)
        p = get!(num_groups, met_part, POLY())
        p[k_part] = get(p, k_part, Rational{Int}(0)) + coeff
    end
    den_groups = Dict{MONO, POLY}()
    for (mono, coeff) in den
        k_part, met_part = split_mono(mono)
        p = get!(den_groups, met_part, POLY())
        p[k_part] = get(p, k_part, Rational{Int}(0)) + coeff
    end
    num_groups, den_groups
end

"""
    _kcat_forward(m::EnzymeMechanism, params) → Float64

Compute kcat (forward) analytically from the polynomial structure.
kcat is the maximum rate at saturating substrates, zero products,
and E_total=1. For mechanisms with multiple catalytic paths
(e.g., non-essential activator), returns the max over all paths.

Groups the rate-equation numerator (forward terms only, positive coefficients)
and denominator by metabolite pattern. For each matching pair, builds k-only
(num_k_expr / den_k_expr) candidates and emits `max(nk/dk, ...)`. Equilibrium
constants cancel between num and den at matching metabolite levels.

Multiple candidates arise for mechanisms with alternative catalytic pathways
(e.g., non-essential activator with/without activator bound).
"""
@generated function _kcat_forward(
    ::M, params::NamedTuple,
) where {M <: EnzymeMechanism}
    num, den, _ = _raw_symbolic_rate_polys(M)
    k_param_names = Set{Symbol}(_raw_param_symbols(M()))
    num_groups, den_groups = _kcat_groups_from_polys(num, den, k_param_names)

    # Build kcat candidates: for each forward numerator metabolite group
    # with a matching denominator group, create (num_k_expr, den_k_expr).
    # kcat is evaluated at products = 0, so product-containing monomials are
    # outside its domain — King–Altman net-flux cross-terms like A·B·P yield
    # spurious candidates that can win the max. Keep substrate-only patterns.
    prod_syms = Set{Symbol}(name(p) for p in products(reaction(Mechanism(M()))))
    components = Tuple{Any, Any}[]
    for (met_key, num_k) in sort!(collect(num_groups); by=first)
        den_k = get(den_groups, met_key, nothing)
        den_k === nothing && continue
        any(first(s) in prod_syms for s in met_key) && continue
        num_expr = _poly_to_expr(num_k)
        den_expr = _poly_to_expr(den_k)
        push!(components, (num_expr, den_expr))
    end

    dep_exprs, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq)
    assignments = [Expr(:(=), sym, dep_exprs[sym])
                   for (sym, _) in sort(collect(dep_exprs); by=first)]
    candidates = [:($nk / $dk) for (nk, dk) in components]
    result = length(candidates) == 1 ?
        candidates[1] : Expr(:call, :max, candidates...)
    Expr(:block,
        _destructuring_expr(hw_params, :params),
        assignments...,
        :(return $result))
end

"""
    _kcat_forward(m::AllostericEnzymeMechanism, params) → Float64

Returns the peak achievable forward turnover: `max` over saturating
substrate patterns and regulator corners (each regulator 0 or saturating)
at products = 0, E_total = 1. Equals the numerical grid-peak forward rate.
Per active site (protomer) — `E_total` is the active-site concentration,
so this carries no `catalytic_multiplicity` factor.
"""
@generated function _kcat_forward(
    ::AllostericEnzymeMechanism{CM,CS,RS},
    params::NamedTuple,
) where {CM,CS,RS}
    M_type = AllostericEnzymeMechanism{CM,CS,RS}
    aem = M_type()
    am  = AllostericMechanism(aem)
    CatN = catalytic_multiplicity(am)

    # Build A-state and I-state polynomials natively per conformation
    # (`_state_rate_polys`), so the saturating metabolite pattern can be matched
    # across conformations. The I-run derives on the reachable-form-pruned graph:
    # a pattern that only exists via an `:OnlyA` step drops out of I-state, and a
    # dead I-cycle yields `num_I = poly_zero()` (empty groups), so its saturating
    # contribution vanishes.
    num_A_poly, den_A_poly, d_free_A = _state_rate_polys(am, :A)
    num_I_poly, den_I_poly, d_free_I = _state_rate_polys(am, :I)
    cat_mets = Set{Symbol}(metabolites(CM()))

    # Same per-state free-enzyme normalization as `_allosteric_num_den_exprs`,
    # applied at the POLY level (this function groups saturating metabolite
    # patterns directly off `num`/`den` polys, not Exprs). The normalization is
    # a common factor of the saturating-limit ratio, so it leaves kcat's value
    # unchanged; matching the same branch as `rate_equation` keeps the
    # saturating-pattern grouping below consistent with it.
    if d_free_A == d_free_I
        # raw — the free-enzyme weight is common to both states and cancels; leave the polys
        # as captured
    elseif _is_metabolite_free_monomial(d_free_A, cat_mets) &&
           _is_metabolite_free_monomial(d_free_I, cat_mets)
        inv_A = _invert_monomial(d_free_A); inv_I = _invert_monomial(d_free_I)
        num_A_poly = poly_mul(num_A_poly, inv_A); den_A_poly = poly_mul(den_A_poly, inv_A)
        num_I_poly = poly_mul(num_I_poly, inv_I); den_I_poly = poly_mul(den_I_poly, inv_I)
    else
        num_A_poly = poly_mul(num_A_poly, d_free_I)
        den_A_poly = poly_mul(den_A_poly, d_free_I)
        num_I_poly = poly_mul(num_I_poly, d_free_A)
        den_I_poly = poly_mul(den_I_poly, d_free_A)
    end

    # Catalytic param-name sets for the metabolite/k split. The A-set is the
    # A-state tagged column set plus any non-metabolite symbol the fold above
    # introduced into the A-polys (e.g. an I-state param pulled in by
    # cross-weighting); the I-set adds the I-polynomials' own params
    # (`:I` mirrors plus the native `:NonequalAI` I-names), which are exactly the
    # non-metabolite symbols the I-polys reference.
    a_param_names = union(
        Set(_param_columns(_state_mechanism(am, :A), _state_step_params(am, :A))),
        setdiff(union(_poly_param_syms(num_A_poly), _poly_param_syms(den_A_poly)),
                cat_mets))
    i_param_names = union(a_param_names,
        setdiff(union(_poly_param_syms(num_I_poly), _poly_param_syms(den_I_poly)),
                cat_mets))
    num_A_groups, den_A_groups =
        _kcat_groups_from_polys(num_A_poly, den_A_poly, a_param_names)
    num_I_groups, den_I_groups =
        _kcat_groups_from_polys(num_I_poly, den_I_poly, i_param_names)

    # kcat = peak forward turnover at saturation: max over saturating patterns
    # (met_key) and regulator corners. Only substrate-saturating patterns are
    # valid at products=0; product-containing patterns are excluded.
    prod_syms = Set{Symbol}(name(p) for p in products(am.reaction))
    a_keys = sort!([k for k in keys(num_A_groups)
                    if haskey(den_A_groups, k) && !any(first(s) in prod_syms for s in k)])
    isempty(a_keys) &&
        error("_kcat_forward: AllostericEnzymeMechanism produced no kcat " *
              "components — saturating-substrate pattern not found in numerator")

    a_assignments, i_assignments_ = _build_dep_assignments(M_type)
    # Keep inactive-state assignments unconditionally: B_I references them, and
    # every I-state dependent the combined solve emits is expressed purely in
    # already-solved columns (independent params or other dependents), so
    # nothing is left undefined.
    i_assignments = i_assignments_
    _, indep = _dependent_param_exprs(M_type)
    hw_params = (indep..., :Keq)
    i_state_dead = _i_state_num_zero(am)

    # Regulator-corner setup (independent of the saturating pattern).
    all_ligs = AllostericRegulator[]
    for site in am.regulatory_sites
        for lig in site.ligands
            lig in all_ligs || push!(all_ligs, lig)
        end
    end
    n_ligs = length(all_ligs)
    lig_idx = Dict(lig => i - 1 for (i, lig) in enumerate(all_ligs))

    kcat_exprs = Any[]
    for met_key in a_keys
        num_k_A_expr = _poly_to_expr(num_A_groups[met_key])
        den_k_A_expr = _poly_to_expr(den_A_groups[met_key])
        num_I_p = get(num_I_groups, met_key, nothing)
        den_I_p = get(den_I_groups, met_key, nothing)
        num_k_I_expr = num_I_p === nothing ? 0 : _poly_to_expr(num_I_p)
        den_k_I_expr = den_I_p === nothing ? 0 : _poly_to_expr(den_I_p)
        i_pattern_dead = den_I_p === nothing

        A_A, B_A = _mwc_power_pair(num_k_A_expr, den_k_A_expr, CatN)
        if i_pattern_dead
            A_I = 0
            B_I = 0
        else
            A_I_live, B_I = _mwc_power_pair(num_k_I_expr, den_k_I_expr, CatN)
            A_I = i_state_dead ? 0 : A_I_live
        end

        if isempty(RS)
            push!(kcat_exprs,
                  :($(_mwc_combine(A_A, A_I)) / $(_mwc_combine(B_A, B_I))))
        else
            for mask in 0:(2^n_ligs - 1)
                W_A_factors = Any[]
                W_I_factors = Any[]
                for (site_idx, site) in enumerate(am.regulatory_sites)
                    n_reg = site.multiplicity
                    sat_terms_A = Any[]
                    sat_terms_I = Any[]
                    for (lig, tag) in zip(site.ligands, site.allo_states)
                        if (mask >> lig_idx[lig]) & 1 == 1
                            if tag !== :OnlyI
                                K_A_sym = name(Kreg(site, lig, :A), am)
                                push!(sat_terms_A, :(inv($K_A_sym)))
                            end
                            if tag !== :OnlyA
                                # `:EqualAI` ligands share the A-state symbol; the
                                # body emits an `:EqualAI` ligand's I-state slot
                                # via the A-state name (no `_T` rename).
                                K_I_state = tag === :EqualAI ? :A : :I
                                K_I_sym = name(Kreg(site, lig, K_I_state), am)
                                push!(sat_terms_I, :(inv($K_I_sym)))
                            end
                        end
                    end
                    if !isempty(sat_terms_A)
                        q_A = length(sat_terms_A) == 1 ?
                            sat_terms_A[1] : _nest_binary(:+, sat_terms_A)
                        push!(W_A_factors, _power_expr(q_A, n_reg))
                    end
                    if !isempty(sat_terms_I)
                        q_I = length(sat_terms_I) == 1 ?
                            sat_terms_I[1] : _nest_binary(:+, sat_terms_I)
                        push!(W_I_factors, _power_expr(q_I, n_reg))
                    end
                end
                if isempty(W_A_factors) && isempty(W_I_factors)
                    kcat_expr =
                        :($(_mwc_combine(A_A, A_I)) / $(_mwc_combine(B_A, B_I)))
                else
                    W_A = isempty(W_A_factors) ? 1 :
                        length(W_A_factors) == 1 ? W_A_factors[1] :
                        _nest_binary(:*, W_A_factors)
                    W_I = isempty(W_I_factors) ? 1 :
                        length(W_I_factors) == 1 ? W_I_factors[1] :
                        _nest_binary(:*, W_I_factors)
                    kcat_expr = :(($(A_A) * $(W_A) + L * $(A_I) * $(W_I)) /
                        ($(B_A) * $(W_A) + L * $(B_I) * $(W_I)))
                end
                push!(kcat_exprs, kcat_expr)
            end
        end
    end

    result = length(kcat_exprs) == 1 ? kcat_exprs[1] :
        Expr(:call, :max, kcat_exprs...)
    return Expr(:block,
        _destructuring_expr(hw_params, :params),
        a_assignments...,
        i_assignments...,
        :(return $result))
end

# ─── Public API: rescale_parameter_values ──────────────────────────

"""
    rescale_parameter_values(m, params::NamedTuple; scale_k_to_kcat=1.0)

Rescale SS rate constants so that `_kcat_forward(m, result) ≈ scale_k_to_kcat`.
Non-SS parameters (K's, Keq, E_total, L, regulatory K's) are unchanged.
"""
function rescale_parameter_values(
    @nospecialize(m::AbstractEnzymeMechanism), params::NamedTuple; scale_k_to_kcat=1.0,
)
    kcat_current = _kcat_forward(m, params)
    scale = scale_k_to_kcat / kcat_current
    ss_names = _ss_rate_constant_names(m)
    NamedTuple{keys(params)}(Tuple(
        k in ss_names ? v * scale : v
        for (k, v) in zip(keys(params), values(params))
    ))
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
The native I-state catalytic numerator polynomial is empty (zero): the
steady-state fluxes of a broken cycle's reachable-form-pruned I-graph cancel
exactly, so the numerator `_raw_symbolic_rate_polys` builds is `poly_zero()`
natively — no forced zero. It decides, at both consumer sites
(`_allosteric_num_den_exprs` and `_kcat_forward`), whether the `L·num_I` term is
emitted and whether `kcat` carries the I-state term. A live redundant-path `:OnlyA`
mechanism (num_I ≠ 0) is thereby handled consistently.
"""
_i_state_num_zero(am::AllostericMechanism) =
    isempty(first(_state_rate_polys(am, :I)))

"""
Names of enzyme forms in the connected component of the free enzyme over ALL
steps of `groups` (rapid-equilibrium and steady-state alike). Every form
carrying neither a bound metabolite nor a residual seeds the search: such a form
interconverts between conformations under formulation 1, and so is a root the
inactive conformation can be entered through. Enumeration yields exactly one
such form; the hand-written DSL admits a second conformation name, and every
match then seeds. A ping-pong covalent intermediate carries a residual and is
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
from the free enzyme over the non-`:OnlyA` groups (`_reachable_from_free`). The
derivation (`_state_allo_mechanism`) and the copy rule (`_redundant_copy_groups`) read
the inactive state from it.
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
The `AllostericMechanism` for `am` in conformational `state`: `am` itself for
`:A`; for `:I`, a fresh `AllostericMechanism` from the nonempty groups of
`_inactive_groups(am)`, with `:OnlyA` catalytic groups
dropped AND every enzyme form disconnected from free E by that drop pruned at
the step level. After removing the `:OnlyA` groups, a form is kept iff it lies
in the connected component of the free-enzyme root — a form carrying neither a
bound metabolite nor a residual, which under formulation 1 is what interconverts
between conformations; enumeration yields exactly one such form, while the
hand-written DSL admits a second conformation name — over ALL remaining steps
(rapid-equilibrium and steady-state alike); a step is kept iff both its endpoints
are kept, so a kinetic group with all its steps dropped disappears. Forms whose only route
back to free E ran through an `:OnlyA` group become disconnected and drop out;
forms still reachable through the surviving steps — including a substrate
complex repopulated by reverse catalysis — are retained. A ping-pong covalent
intermediate carries a residual and so is never itself a root. A mechanism
with no `:OnlyA` catalytic group keeps its whole graph and re-derives its full
native I-state.

Routing both the derivation mechanism and the step_params through this one
struct keeps steps and allo-state tags aligned: the `AllostericMechanism`
constructor canonicalizes `cat_steps` and applies the SAME permutation to
`cat_allo_states`, so pruning-induced reorder/iso-flips can't desync the two.
"""
function _state_allo_mechanism(am::AllostericMechanism, state::Symbol)
    state === :I || return am
    groups = _inactive_groups(am)
    keep = findall(!isempty, groups)
    AllostericMechanism(reaction(am), groups[keep], cat_allo_states(am)[keep],
                        catalytic_multiplicity(am), regulatory_sites(am))
end

"""
State-tagged `step_params` for `am`'s catalytic mechanism in conformational
`state` (`:A` or `:I`). Same shape as `_step_parameters(::Mechanism)` — a
per-flat-step vector of `Parameter`s — but each catalytic group's `Parameter`s
carry the group's state tag: a `:NonequalAI`/`:OnlyA` group is tagged with
`state`; an `:EqualAI` group is tagged `:EqualAI` (so `name(p, am)` renders the
shared bare Symbol in both states). For `state == :I`, `:OnlyA` groups are
already pruned from `_state_allo_mechanism(am, :I)`, matching the broken-cycle
graph from `_state_mechanism(am, :I)`.

Walks `_state_allo_mechanism`'s already-canonical, aligned steps/states so the
per-flat-step order matches `_flat_steps(_state_mechanism(am, state))`. Each
`Parameter` is anchored on its own step (not the rep), exactly as
`_step_parameters` does, so arity follows the step's RE/SS type while
`name(p, am)` collapses to the rep's structural Symbol via the chokepoint.
"""
function _state_step_params(am::AllostericMechanism, state::Symbol)
    sam = _state_allo_mechanism(am, state)
    out = Vector{Vector{Parameter}}()
    for (g, group) in enumerate(steps(sam))
        tag = cat_allo_state(sam, g) === :EqualAI ? :EqualAI : state
        for s in group
            push!(out, is_equilibrium(s) ?
                Parameter[is_binding(s) ? Kd(s, tag) : Kiso(s, tag)] :
                Parameter[is_binding(s) ? Kon(s, tag)  : Kfor(s, tag),
                          is_binding(s) ? Koff(s, tag) : Krev(s, tag)])
        end
    end
    out
end

"""
The catalytic `Mechanism` for `am` in conformational `state`. `:A` keeps the
full catalytic graph; `:I` keeps only the reachable-form subgraph (`:OnlyA`
groups and the forms they disconnect from free E pruned, via
`_state_allo_mechanism`) so King–Altman re-derives the broken-cycle I-state law
natively.
"""
function _state_mechanism(am::AllostericMechanism, state::Symbol)
    sam = _state_allo_mechanism(am, state)
    Mechanism(reaction(sam), steps(sam))
end

"""
Derive `(num_poly, den_poly, d_free_poly)` for `am`'s catalytic mechanism in
conformational `state`, natively in that state's parameter names. Runs the
shared King–Altman engine on the state-tagged `step_params` and state graph, so
no post-hoc rename is needed (`:EqualAI` groups render the shared bare Symbol
automatically). The `:I` polynomials reference each `:NonequalAI` group's
native `K_I_…`/`k_I_…` symbol; a forbidden split's `K_I_…` is defined by the
combined constraint solve's dependent assignment (`_build_dep_assignments`).
`d_free_poly` is free E's weight in that state's denominator (see
`_raw_symbolic_rate_polys`).
"""
function _state_rate_polys(am::AllostericMechanism, state::Symbol)
    cm = _state_mechanism(am, state)
    _assert_derivable(cm)
    sp = _state_step_params(am, state)
    @assert length(sp) == length(_flat_steps(cm)) "state step_params/steps misaligned"
    _raw_symbolic_rate_polys(cm, sp, _state_wegscheider_rename_map(am, state))
end

"""
State-tagged Wegscheider rename map for `am`'s catalytic sub-mechanism in
conformational `state` — the state-aware analog of `_build_wegscheider_rename_map`
(which runs the kernel with `:None` step_params and so cannot see the `:A`/`:I`
tags). Discovers single-symbol RE binding-K Wegscheider ties (`K_a = K_b`, both
binding K's) under the state-tagged `step_params` and folds each
absorbed symbol into its target. Empty for all current specs (catalysis is
steady-state, so no fully-RE catalytic box), but a fully-RE catalytic core would
now collapse its tie natively — the same way the non-allosteric path does.
"""
function _state_wegscheider_rename_map(am::AllostericMechanism, state::Symbol)
    cm = _state_mechanism(am, state)
    sp = _state_step_params(am, state)
    rename = Dict{Symbol, Symbol}()
    # binding-K set: value-context rep name of each RE binding step.
    binding_set = Set{Symbol}()
    for (idx, (s, _)) in enumerate(_flat_steps(cm))
        is_equilibrium(s) && is_binding(s) || continue
        push!(binding_set, name(sp[idx][1], cm))
    end
    # Single-symbol Wegscheider RE ties between two binding K's.
    dep_raw, _ = _dependent_param_exprs_kernel(cm, rename; step_params = sp)
    for (lhs, rhs) in dep_raw
        rhs isa Symbol || continue
        lhs in binding_set && rhs in binding_set || continue
        target = get(rename, rhs, rhs)
        rename[lhs] = target
        for k in collect(keys(rename))
            rename[k] == lhs && (rename[k] = target)
        end
    end
    rename
end

"""
Catalytic `Parameter`s of `am` in conformation `state` (`:A` or `:I`): the
constants of each kinetic group's rep step (`Kd`/`Kiso`/`Kon`+`Koff`/`Kfor`+`Krev`),
in group order. In `:A` an `:EqualAI` group takes the `:EqualAI` tag, because its
symbol is shared with the I-state (the chokepoint `name(p, m)` renders both to the
same `Symbol`), and every other group takes `:A`. In `:I` the `:OnlyA` groups are
skipped and every other group takes `:I`. Synthesized-dep I-mirrors (deps whose RHS
references a `:NonequalAI` symbol) belong to dep-parameter machinery and are emitted
Symbol-level by the dep-assignment builder.
"""
function _cat_params(am::AllostericMechanism, state::Symbol)
    out = Parameter[]
    for (g, rep) in enumerate(_group_reps(am))
        tag = cat_allo_state(am, g)
        state === :I && tag === :OnlyA && continue
        st = state === :A && tag === :EqualAI ? :EqualAI : state
        append!(out, _step_constants(rep, st))
    end
    out
end

"""Regulator-site `Kreg`s of `am` in conformation `state` (`:A` or `:I`), site by site:
a ligand absent from that conformation (`:OnlyI` in `:A`, `:OnlyA` in `:I`) has none."""
_kreg_params(am::AllostericMechanism, state::Symbol) =
    Kreg[Kreg(site, lig, state) for site in regulatory_sites(am)
         for (lig, tag) in zip(ligands(site), allo_states(site))
         if tag !== (state === :A ? :OnlyI : :OnlyA)]

# ─── Dependent parameter expressions ─────────────────────────────

"""
    _combined_state_dependent_exprs(am::AllostericMechanism)

Stack the A-state and I-state thermodynamic constraint systems
(`_assemble_constraints`, tagged `is_i_state`) over their combined column
space and solve once with the shared solver (`_solve_dependent_set`). The
`(is_i_state, type)` pivot priority alone gives the I-above-A collapse
direction — an I-state column always outranks its A-state counterpart, so a
cross-state affinity tie (e.g. a live-forbidden `:NonequalAI` split) is
expressed onto the free A-side directly, with no post-hoc merge step; within
a state the catalytic pivot order (`_step_priority`) is unchanged. Returns
`(dep_exprs, indep_params)` over the union of both states' catalytic columns.
"""
function _combined_state_dependent_exprs(am::AllostericMechanism)
    function state_system(state)
        cm = _state_mechanism(am, state)
        sp = _state_step_params(am, state)
        _assemble_constraints(cm, _state_wegscheider_rename_map(am, state);
                              step_params = sp, is_i_state = (state === :I))
    end
    A_A, rhs_A, cols_A, pri_A = state_system(:A)
    A_I, rhs_I, cols_I, pri_I = state_system(:I)

    # Union columns: A-state first, then any I-only column (a shared `:EqualAI` group
    # carries the same bare Symbol in both states and coincides — it keeps its A tag).
    columns = copy(cols_A)
    col_index = Dict(c => i for (i, c) in enumerate(columns))
    priority = copy(pri_A)
    for (j, c) in enumerate(cols_I)
        haskey(col_index, c) && continue
        push!(columns, c)
        col_index[c] = length(columns)
        push!(priority, pri_I[j])
    end

    # Stack the two per-state constraint blocks over the combined column space and
    # solve once. Cross-state ties emerge as `I-row − A-row` (the `log Keq` cancels).
    A = zeros(Rational{BigInt}, size(A_A, 1) + size(A_I, 1), length(columns))
    for i in axes(A_A, 1), (j, c) in enumerate(cols_A)
        A_A[i, j] == 0 || (A[i, col_index[c]] = A_A[i, j])
    end
    off = size(A_A, 1)
    for i in axes(A_I, 1), (j, c) in enumerate(cols_I)
        A_I[i, j] == 0 || (A[off + i, col_index[c]] = A_I[i, j])
    end
    rhs = vcat(rhs_A, rhs_I)
    return _solve_dependent_set(A, rhs, columns, priority)
end

"""
    _dependent_param_exprs(am::AllostericMechanism)

Return `(dep_exprs, indep_params)` for an allosteric mechanism from the
single combined constraint solve (`_combined_state_dependent_exprs`), which
stacks the A-state and I-state constraint rows and solves them once — a
cross-state tie (e.g. a live-forbidden `:NonequalAI` split) falls out of the
stacked system directly, with no post-hoc merge step.

Regulator-site affinities complete no catalytic thermodynamic cycle, so they
are independent on top of the combined solve — except an `:EqualAI`
regulator, whose I-name mirrors its shared A-name (`K_I_reg = K_A_reg`, added
to `dep`). `L` (the conformational constant) is always independent. Any
symbol the combined solve already made dependent is dropped from `indep`. The
`Type{<:AbstractEnzymeMechanism}` method lifts with `_concrete` and delegates here.
"""
function _dependent_param_exprs(am::AllostericMechanism)
    dep, indep = _combined_state_dependent_exprs(am)

    # A per-state Wegscheider rename folds a single-symbol binding-K tie onto one
    # representative; the folded symbol enters the combined solve as a zero-column and
    # lands in `indep` as a fittable dummy. Drop it — unless a retained polynomial
    # still references it (a shared `:EqualAI` symbol can be folded in one state yet
    # used by the other state's polynomial, which keeps its own rename). The
    # non-allosteric analog needs no reference guard because it has a single,
    # consistent rename. No-op when both state renames are empty (all current specs).
    rename = merge(_state_wegscheider_rename_map(am, :A),
                   _state_wegscheider_rename_map(am, :I))
    if !isempty(rename)
        refs = Set{Symbol}()
        for st in (:A, :I)
            num, den, _ = _state_rate_polys(am, st)
            union!(refs, _poly_param_syms(den))
            isempty(num) || union!(refs, _poly_param_syms(num))
        end
        indep = Tuple(p for p in indep if get(rename, p, p) == p || p in refs)
    end

    # Regulator-site affinities complete no catalytic thermodynamic cycle, so they
    # are independent — except an `:EqualAI` regulator, whose I-name mirrors its
    # shared A-name: the mirror joins `dep`, which drops it from the independent list.
    # `L` (the conformational constant) is always independent.
    for site in regulatory_sites(am), (lig, tag) in zip(ligands(site), allo_states(site))
        tag === :EqualAI || continue
        dep[name(Kreg(site, lig, :I), am)] = name(Kreg(site, lig, :A), am)
    end
    reg_params = Symbol[name(p, am) for p in [_kreg_params(am, :A); _kreg_params(am, :I)]]
    return dep, Tuple(p for p in (indep..., reg_params..., :L) if p ∉ keys(dep))
end

_dependent_param_exprs(@nospecialize(M::Type{<:AbstractEnzymeMechanism})) =
    _dependent_param_exprs(_concrete(M()))

# `parameters` and `fitted_params` for `AllostericEnzymeMechanism` are the
# `AbstractEnzymeMechanism` methods at the top of this file.

# ─── Rate body building helpers ───────────────────────────────────

"""Build the regulatory site partition function expression: 1 + lig/K_lig_reg_i + ...
Skips ligands absent from the given conformation (`:OnlyA` in I-state, `:OnlyI`
in A-state). Uses the A-state K symbol when the ligand tag is `:EqualAI`.
Renders K-names via the `name(::Kreg, am)` chokepoint."""
function _reg_site_expr(am::AllostericMechanism, site_idx::Int, inactive::Bool)
    site = regulatory_sites(am)[site_idx]
    terms = Any[1]
    for (lig, tag) in zip(ligands(site), allo_states(site))
        if inactive
            tag === :OnlyA && continue
        else
            tag === :OnlyI && continue
        end
        # `:EqualAI` ligands share the A-state symbol in both conformations;
        # `:NonequalAI` / `:OnlyI` ligands carry a distinct I-state K name.
        state = (inactive && tag in (:NonequalAI, :OnlyI)) ? :I : :A
        K_sym = name(Kreg(site, lig, state), am)
        push!(terms, :($(name(lig)) / $K_sym))
    end
    _nest_binary(:+, terms)
end

"""Raise an expression to an integer power (returns 1 for n=0, expr for n=1)."""
function _power_expr(expr, n::Int)
    n == 0 && return 1
    n == 1 && return expr
    :(($expr)^$n)
end

"""MWC active + L·inactive state-combine `A + L * B`. Shared by `_kcat_forward`
(numerator/denominator halves of the state ratio) and `_allosteric_num_den_exprs`
(the retained num/den sums)."""
_mwc_combine(a, b) = :($a + L * $b)

"""MWC binding-statistics power-pair `(X * Y^(n-1), Y^n)` for a saturating pattern:
the numerator carries one fewer denominator power than the denominator, where
`n = catalytic_multiplicity`. Used by `_kcat_forward` per conformation."""
_mwc_power_pair(x, y, n) =
    (n == 1 ? x : :($x * $y^$(n - 1)), :($y^$n))

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
The set of Symbols that belong to the I-block when `_build_dep_assignments`
splits the combined solve's dependents for emission: catalytic columns that
appear only in the I-state's parameter set (`cols_I \\ cols_A` — a shared
`:EqualAI` catalytic Symbol coincides with its A-state column, so it stays in
the A-block), plus every non-`:OnlyA` regulator's I-name.
"""
function _i_state_symbol_set(am::AllostericMechanism)
    cols_A = _param_columns(_state_mechanism(am, :A), _state_step_params(am, :A))
    cols_I = _param_columns(_state_mechanism(am, :I), _state_step_params(am, :I))
    syms = Set{Symbol}(setdiff(cols_I, cols_A))
    union!(syms, (name(p, am) for p in _kreg_params(am, :I)))
end

"""
Build active-state and inactive-state dep-param assignment Exprs from the
single combined solve (`_dependent_param_exprs`). Returns `(a_assignments::
Vector{Expr}, i_assignments::Vector{Expr})`. Shared by
`_build_allosteric_rate_body` and `rate_equation_string`.

Splits `dep` by `_i_state_symbol_set`: a dependent whose LHS is an I-only
catalytic column or a non-`:OnlyA` regulator I-name (including an `:EqualAI`
reg mirror, `K_I_reg = K_A_reg`) is emitted in the I-block; everything else
(including a shared `:EqualAI` catalytic dependent, whose bare Symbol
coincides with its A-state column) is emitted in the A-block. All Symbols
route through the `name(p, am)` chokepoint (native derivation + `Kreg`).
"""
function _build_dep_assignments(
    @nospecialize(M_type::Type{<:AllostericEnzymeMechanism}),
)
    am = AllostericMechanism(M_type())
    dep, _ = _dependent_param_exprs(am)
    # A-block first so an I-block `:EqualAI` regulator mirror (`K_I_reg = K_A_reg`)
    # finds its A-name defined. The combined solve expresses every dependent purely
    # in independent columns, so no dependent reads another — order within a block
    # is free.
    i_syms = _i_state_symbol_set(am)
    a_assignments = Expr[]
    i_assignments = Expr[]
    for (sym, rhs) in sort(collect(dep); by = first)
        push!(sym in i_syms ? i_assignments : a_assignments, Expr(:(=), sym, rhs))
    end
    return a_assignments, i_assignments
end

"""
Assemble the MWC numerator and denominator Exprs.
Returns `(full_num, full_den)`. Per-active-site normalization: the
numerator carries no leading `catalytic_multiplicity` factor; only the
`Q_cat^(CatN-1)` / `Q_cat^CatN` binding-statistics powers remain.
"""
function _allosteric_num_den_exprs(@nospecialize(M_type::Type{<:AllostericEnzymeMechanism}))
    m = M_type()
    am = AllostericMechanism(m)
    CM = typeof(catalytic_mechanism(m))
    CatN = catalytic_multiplicity(m)
    RS = regulatory_sites(am)

    num_A_poly, den_A_poly, d_free_A = _state_rate_polys(am, :A)
    # A-state catalytic param symbols (the tagged column set) drive `_poly_to_expr`'s
    # param/metabolite ordering split; the I-poly's `:I` symbols sort as non-params.
    cat_params = Set(_param_columns(_state_mechanism(am, :A), _state_step_params(am, :A)))
    cat_mets = Set{Symbol}(metabolites(CM()))

    # I-state catalytic polys, always re-derived natively on the reachable-form
    # subgraph (`_state_allo_mechanism(am, :I)` drops `:OnlyA` groups and every
    # form they disconnect from free E). Reachable-subgraph King–Altman gives the
    # same binding partition monomial-zeroing produced, and for a dead cycle the
    # pruned graph's steady-state fluxes cancel exactly, so the numerator
    # `_raw_symbolic_rate_polys` builds is 0 natively — no forced zero needed.
    num_i_poly, den_i_poly, d_free_I = _state_rate_polys(am, :I)

    # Formulation-1 per-state free-enzyme normalization. Render the same value
    # three ways by how the two free-enzyme weights combine:
    #   D_A == D_I               → raw (identical conformations; the factor cancels)
    #   both metabolite-free monomials → divide Q/D (clean standard-MWC form)
    #   otherwise                → cross-weight by the other state's D^n (polynomial)
    D_A_expr = 1
    D_I_expr = 1
    if d_free_A == d_free_I
        # raw combine — leave the polynomials and D exprs as identities
    elseif _is_metabolite_free_monomial(d_free_A, cat_mets) &&
           _is_metabolite_free_monomial(d_free_I, cat_mets)
        inv_A = _invert_monomial(d_free_A)
        inv_I = _invert_monomial(d_free_I)
        num_A_poly = poly_mul(num_A_poly, inv_A); den_A_poly = poly_mul(den_A_poly, inv_A)
        num_i_poly = poly_mul(num_i_poly, inv_I); den_i_poly = poly_mul(den_i_poly, inv_I)
    else
        D_A_expr = _poly_to_expr(d_free_A, cat_params)
        D_I_expr = _poly_to_expr(d_free_I, cat_params)
    end

    N_A = _poly_to_expr(num_A_poly, cat_params)
    Q_A = _poly_to_expr(den_A_poly, cat_params)
    N_I = _poly_to_expr(num_i_poly, cat_params)
    Q_I = _poly_to_expr(den_i_poly, cat_params)

    reg_Q_A = Any[_reg_site_expr(am, i, false) for i in eachindex(RS)]
    reg_Q_I = Any[_reg_site_expr(am, i, true) for i in eachindex(RS)]

    # Numerator: N × Q_cat^(CatN-1) × all reg-site factors at multiplicity.
    function make_num_term(N, Q, reg_Qs)
        factors = Any[N]
        CatN > 1 && push!(factors, _power_expr(Q, CatN - 1))
        for i in eachindex(RS)
            push!(factors, _power_expr(reg_Qs[i], multiplicity(RS[i])))
        end
        _nest_binary(:*, factors)
    end

    # Denominator: Q_cat^CatN × all reg-site factors at multiplicity.
    function make_den_term(Q, reg_Qs)
        factors = Any[_power_expr(Q, CatN)]
        for i in eachindex(RS)
            push!(factors, _power_expr(reg_Qs[i], multiplicity(RS[i])))
        end
        _nest_binary(:*, factors)
    end

    num_A = _mwc_cross_weight(make_num_term(N_A, Q_A, reg_Q_A), D_I_expr, CatN)
    den_A = _mwc_cross_weight(make_den_term(Q_A, reg_Q_A), D_I_expr, CatN)
    den_I = _mwc_cross_weight(make_den_term(Q_I, reg_Q_I), D_A_expr, CatN)
    full_den = _mwc_combine(den_A, den_I)

    if isempty(num_i_poly)
        # Native I-state numerator is zero (`_i_state_num_zero`): the I-state
        # cycle is dead, so drop the L*num_I term entirely (skip dead numerator
        # branch). Q_I still contributes to denominator as enzyme mass.
        num_A, full_den
    else
        num_I = _mwc_cross_weight(make_num_term(N_I, Q_I, reg_Q_I), D_A_expr, CatN)
        _mwc_combine(num_A, num_I), full_den
    end
end

"""Build the MWC rate equation body as an Expr block."""
function _build_allosteric_rate_body(
    @nospecialize(M_type::Type{<:AllostericEnzymeMechanism}),
)
    full_num, full_den = _allosteric_num_den_exprs(M_type)
    rate_expr = :(E_total * ($full_num) / ($full_den))

    a_assignments, i_assignments_ = _build_dep_assignments(M_type)
    # Keep inactive-state assignments unconditionally: the retained Q_I
    # (`L * den_I`) references them, and every I-state dependent the combined
    # solve emits is expressed purely in already-solved columns (independent
    # params or other dependents), so nothing is left undefined.
    i_assignments = i_assignments_

    _, indep = _dependent_param_exprs(M_type)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(M_type())

    Expr(:block,
        _destructuring_expr(hw_params, :params),
        _destructuring_expr(mets, :concs),
        a_assignments...,
        i_assignments...,
        :(return $rate_expr))
end

# ─── Rate equation dispatch ───────────────────────────────────────

@generated function rate_equation(
    ::AllostericEnzymeMechanism{CM,CS,RS},
    concs::NamedTuple, params::NamedTuple, ::ReducedMode,
) where {CM,CS,RS}
    _build_allosteric_rate_body(AllostericEnzymeMechanism{CM,CS,RS})
end

# ─── String representation ────────────────────────────────────────

function rate_equation_string(
    @nospecialize(m::AllostericEnzymeMechanism), ::ReducedMode,
)
    M = typeof(m)
    _, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(m)

    # Every dependent assignment comes from the single combined solve — the same set
    # the compiled body assigns — split into Wegscheider/Haldane by Keq-reference.
    a_assignments, i_assignments = _build_dep_assignments(M)
    weg_lines, hal_lines = String[], String[]
    for a in (a_assignments..., i_assignments...)
        sym = a.args[1]
        expr = a.args[2]
        is_haldane = _mentions(expr, :Keq)
        line = "$sym = $(_expr_to_string(expr))"
        push!(is_haldane ? hal_lines : weg_lines, line)
    end

    # Sort each section lexicographically — load-bearing for eq_hash
    # dedup of allosteric Source-C clusters, since inactive-state lines are
    # appended in iteration order rather than the lexicographic active-state
    # order.
    sort!(weg_lines)
    sort!(hal_lines)

    full_num, full_den = _allosteric_num_den_exprs(M)
    v_line = "v = E_total * ($(_expr_to_string(full_num))) / ($(_expr_to_string(full_den)))"

    lines = ["(; $(join(hw_params, ", "))) = params",
             "(; $(join(mets, ", "))) = concs"]
    _append_constraint_sections!(lines, weg_lines, hal_lines)
    push!(lines, v_line)
    join(lines, "\n")
end
