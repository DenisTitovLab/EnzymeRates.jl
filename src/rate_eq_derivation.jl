# ABOUTME: King-Altman/Cha rate equation derivation via @generated functions.
# ABOUTME: Parameters API, kcat computation, and MWC allosteric rate assembly.

# ─── Parameters API ─────────────────────────────────────────

# Every mechanism is a singleton type of its own, so a method specialized on one
# compiles again for each mechanism. The derivation helpers that lift a singleton,
# or its type, to its concrete form therefore take it `@nospecialize`.

"""
Suffix the `Reduced` `rate_equation_string` of an `EnzymeMechanism` appends to the line
of each single-symbol Wegscheider tie that the rename (`_build_wegscheider_rename_map`)
folds into `v`. The rate-equation dedup key in `identify_rate_equation.jl` strips these
provenance lines.
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

@generated parameters(::EnzymeMechanism{Sig}, ::FullMode) where {Sig} =
    (_raw_param_symbols(EnzymeMechanism{Sig}())..., :E_total)

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
binding K's of `mech` under `step_params`. Solves the constraint system with
no rename to discover binding-K-to-binding-K Wegscheider closures of the form
`K_a = K_b` (RHS is a bare Symbol). Both sides must be binding K's (RE step with
metabolite on LHS) — absorbing a binding-K-to-iso-K tie would produce
inconsistent sign-flips when the solve runs with the full rename,
since the binding-K column is sign-flipped (Kd convention) but the
iso-K column is not. No tie chains onto another: the solve is a full
Gauss–Jordan elimination, so a dependent's right-hand side holds only
independent columns, and no target is itself renamed.

`step_params` defaults to the mechanism's own `:None`-state constants. The
allosteric derivation passes each conformation's state-tagged ones
(`_state_parts`), so the map holds that state's `K_A_…`/`K_I_…` names. For
example, a tie arises from an RE binding square whose two bindings of one ligand
share a group while the other ligand's two bindings do not: the square then
equates the other ligand's two K's.

The rename means the polynomial in `v` uses the representative symbol
directly, so Source-C duplicates (split kinetic groups that
Wegscheider ties back together) collapse at hash time.
"""
function _build_wegscheider_rename_map(mech::Mechanism;
                                       step_params = _step_parameters(mech))
    # binding-K set: value-context rep name of each RE binding step. Walk
    # Mechanism.steps directly — an RE step that is a binding (`is_binding`,
    # plain or fused) is a binding step; step_params is indexed in the same
    # flat order.
    binding_set = Set{Symbol}(name(step_params[j][1], mech)
                              for (j, (s, _)) in enumerate(_flat_steps(mech))
                              if is_equilibrium(s) && is_binding(s))
    # Pass 2: single-symbol Wegscheider RE ties between two binding K's.
    dep_raw, _ = _solve_dependent_set(
        _assemble_constraints(mech, Dict{Symbol, Symbol}(); step_params)...)
    Dict{Symbol, Symbol}(lhs => rhs for (lhs, rhs) in dep_raw
                         if rhs isa Symbol && lhs in binding_set && rhs in binding_set)
end

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
`name(p, m)` chokepoint, then through any Wegscheider tie). Raises when the RE
steps close a catalytic cycle, which gives the mechanism no finite rate.
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
collapses kinetic-group members to their rep's name). `rename_map` applies any
single-symbol Wegscheider ties at the leaves: each RE step's constant is renamed
as it enters `step_to_K`, the only route by which an RE binding K (every rename key
and target is one) reaches the polynomials. `step_params` defaults
to the mechanism's own `:None`-state constants and `rename_map` to the Wegscheider
rename under them; the allosteric derivation passes each conformation's
state-tagged constants (`_state_parts`). Aborts, before the King–Altman expansion,
a mechanism whose denominator would be too large to derive (`_assert_derivable`).

Also returns `d_free`, the weight in the returned denominator of the free
resting enzyme (the form with empty `bound` and empty `residual`): the
spanning-tree weight `D[g_free]` of its segment times the concentration monomial
that brings `num`/`den` to lowest terms. Free E roots its segment, so its own
`alpha` is 1.
"""
function _raw_symbolic_rate_polys(
    mech::Mechanism,
    step_params = _step_parameters(mech),
    rename_map = _build_wegscheider_rename_map(mech; step_params),
)
    _assert_derivable(mech)
    species, segments, _, idx, seg = _re_segment_extras(steps(mech))
    # A fully-inert conformation (every binding pruned) has no enumerated form; it
    # exists only as free enzyme — no flux, partition 1, D[g_free] 1.
    isempty(species) && return poly_zero(), poly_one(), poly_one()
    flat = _flat_steps(mech)
    step_to_K = Dict{Int, Symbol}(
        i => (K = name(step_params[i][1], mech); get(rename_map, K, K))
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

    conc_set = _concentration_symbols(mech)
    _reduce_conc_lowest_terms(num, den, d_free, conc_set)
end

_raw_symbolic_rate_polys(@nospecialize(M::Type{<:EnzymeMechanism})) =
    _raw_symbolic_rate_polys(Mechanism(M()))

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

"""The numerator and denominator Exprs of the rate equation of `M`, from its raw
symbolic rate polys."""
function _num_den_exprs(@nospecialize(M::Type{<:EnzymeMechanism}))
    num, den, _ = _raw_symbolic_rate_polys(M)
    param_syms = Set{Symbol}(_raw_param_symbols(M()))
    _poly_to_expr(num, param_syms), _poly_to_expr(den, param_syms)
end

"""Build destructuring Expr: (; a, b, c) = source"""
function _destructuring_expr(syms, source::Symbol)
    Expr(:(=), Expr(:tuple, Expr(:parameters, syms...)), source)
end

"""
The assignments `sym = rhs` of the dependent parameters `dep`, in name order.
`_solve_dependent_set` is a full Gauss–Jordan elimination, so every right-hand side
reads only `Keq` and independent parameters — as does the `:EqualAI` regulator mirror
`K_I_reg = K_A_reg` that the allosteric `_dependent_param_exprs` adds — and the
assignments may run in any order.
"""
_dep_assignments(dep) =
    [Expr(:(=), sym, rhs) for (sym, rhs) in sort!(collect(dep); by = first)]

"""
The body of the generated `rate_equation` of `M`: destructure `param_syms` from
`params` and the metabolites from `concs`, assign the dependent parameters `dep`, and
return `E_total * (num) / (den)` with `num` and `den` from `_num_den_exprs`.
"""
function _rate_body(@nospecialize(M::Type{<:AbstractEnzymeMechanism}), param_syms, dep)
    num, den = _num_den_exprs(M)
    Expr(:block,
        _destructuring_expr(param_syms, :params),
        _destructuring_expr(metabolites(M()), :concs),
        _dep_assignments(dep)...,
        :(return E_total * ($num) / ($den)))
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
    _rate_body(M, (_raw_param_symbols(M())..., :E_total), Dict())
end

@generated function rate_equation(
    m::M, concs::NamedTuple, params::NamedTuple,
    ::ReducedMode,
) where {M <: AbstractEnzymeMechanism}
    dep, indep = _dependent_param_exprs(M)
    _rate_body(M, (indep..., :Keq, :E_total), dep)
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

"""
The `rate_equation_string` text of `M`: the `params` destructure of `param_syms`, the
`concs` destructure, one line `sym = rhs` per dependent parameter of `dep` with `rhs`
rendered by `rhs_string`, and the `v = E_total * (num) / (den)` line from
`_num_den_exprs`. A dependent whose right-hand side mentions `Keq` goes under
`# Haldane constraints:`, every other one under `# Wegscheider constraints:`; the line
of a dependent in `substituted` ends in `ANNOTATION_SUBSTITUTED`. Each section is sorted
by name, which is load-bearing: the eq_hash dedup of rate-equivalent mechanisms compares
these strings, so their line order must not depend on the order the solve emits its
dependents in.
"""
function _equation_text(@nospecialize(M::Type{<:AbstractEnzymeMechanism}), param_syms,
                        dep; substituted = (), rhs_string = _expr_to_string)
    weg, hal = String[], String[]
    for (sym, rhs) in sort!(collect(dep); by = first)
        suffix = sym in substituted ? ANNOTATION_SUBSTITUTED : ""
        push!(_mentions(rhs, :Keq) ? hal : weg, "$sym = $(rhs_string(rhs))$suffix")
    end
    num, den = _num_den_exprs(M)
    lines = ["(; $(join(param_syms, ", "))) = params",
             "(; $(join(metabolites(M()), ", "))) = concs"]
    isempty(weg) || push!(lines, "# Wegscheider constraints:", weg...)
    isempty(hal) || push!(lines, "# Haldane constraints:", hal...)
    push!(lines, "v = E_total * ($(_expr_to_string(num))) / ($(_expr_to_string(den)))")
    join(lines, "\n")
end

rate_equation_string(@nospecialize(m::EnzymeMechanism), ::FullMode) =
    _equation_text(typeof(m), (_raw_param_symbols(m)..., :E_total), Dict())

# The display solves with no Wegscheider rename, so an absorbed single-symbol tie stays
# visible under `# Wegscheider constraints:`, marked as substituted into v. Its
# right-hand sides print with Base `string` because the eq_hash of a plain mechanism is
# computed over this text. Base `string` writes a symbol that is not an identifier (a
# ping-pong residual name such as `K_EB_res_+A_-P_to_EQ`) inside a product as `var"…"`.
function rate_equation_string(@nospecialize(m::EnzymeMechanism), ::ReducedMode)
    _, indep = _dependent_param_exprs(typeof(m))
    mech = Mechanism(m)
    dep, _ = _solve_dependent_set(_assemble_constraints(mech, Dict{Symbol, Symbol}())...)
    _equation_text(typeof(m), (indep..., :Keq, :E_total), dep;
                   substituted = keys(_build_wegscheider_rename_map(mech)),
                   rhs_string = string)
end

function rate_equation_string(@nospecialize(m::AllostericEnzymeMechanism), ::ReducedMode)
    dep, indep = _dependent_param_exprs(typeof(m))
    _equation_text(typeof(m), (indep..., :Keq, :E_total), dep)
end

# ─── kcat Computation Helpers ──────────────────────────────────

"""
Set of Symbol names for SS rate-constant parameters (`Kfor`, `Krev`) of
`em`. For `AllostericEnzymeMechanism`, also includes the I-state
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
    Set{Symbol}(name(p, m) for p in params if p isa Union{Kfor, Krev})
end

"""Group `num` and `den` POLYs by metabolite monomial pattern: each monomial splits
into its part over the metabolites `mets` (the group key) and its parameter part.
Returns `(num_groups, den_groups)` where each value is a POLY of parameter monomials
sharing the same metabolite monomial. Reverse (negative-coefficient) terms are dropped
from the numerator. Used by `_kcat_forward` to compare saturating metabolite patterns
across A/I states."""
function _kcat_groups_from_polys(num::POLY, den::POLY, mets::Set{Symbol})
    function group(p, forward_only)
        groups = Dict{MONO, POLY}()
        for (mono, coeff) in p
            forward_only && coeff <= 0 && continue
            g = get!(groups, filter(q -> q.first in mets, mono), POLY())
            k = filter(q -> q.first ∉ mets, mono)
            g[k] = get(g, k, 0) + coeff
        end
        groups
    end
    group(num, true), group(den, false)
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
    num_groups, den_groups = _kcat_groups_from_polys(num, den, Set(metabolites(M())))
    prods = Set(name.(products(reaction(Mechanism(M())))))
    candidates = [:($(_poly_to_expr(num_groups[k])) / $(_poly_to_expr(den_groups[k])))
                  for k in _kcat_keys(num_groups, den_groups, prods)]
    dep, indep = _dependent_param_exprs(M)
    Expr(:block,
        _destructuring_expr((indep..., :Keq), :params),
        _dep_assignments(dep)...,
        :(return $(_max_expr(candidates))))
end

"""
The saturating metabolite patterns `_kcat_forward` takes its candidates from: the keys
of `num_groups` that are also keys of `den_groups` and name no product of `prods`, in
sorted order. kcat is evaluated at products = 0, so product-containing monomials are
outside its domain — King–Altman net-flux cross-terms like A·B·P yield spurious
candidates that can win the max. Errors when no pattern remains, as when catalysis runs
only on a product-bound form.
"""
function _kcat_keys(num_groups, den_groups, prods)
    ks = sort!([k for k in keys(num_groups)
                if haskey(den_groups, k) && !any(first(s) in prods for s in k)])
    isempty(ks) && error("_kcat_forward: no kcat components — no product-free " *
        "saturating-substrate pattern appears in both the rate numerator and denominator")
    ks
end

"""The `max` of the kcat candidate Exprs `cands`, or the lone candidate itself."""
_max_expr(cands) = length(cands) == 1 ? only(cands) : Expr(:call, :max, cands...)

"""
    _kcat_forward(m::AllostericEnzymeMechanism, params) → Float64

Returns the peak achievable forward turnover: `max` over saturating
substrate patterns and regulator corners (each regulator 0 or saturating)
at products = 0, E_total = 1. Equals the numerical grid-peak forward rate. At a
corner whose saturating regulators bind one conformation at a higher total site
multiplicity than the other (an `:OnlyA` or `:OnlyI` ligand), that conformation holds
all the enzyme, so the corner is its own turnover.
Per active site (protomer) — `E_total` is the active-site concentration,
so this carries no `catalytic_multiplicity` factor.
"""
@generated function _kcat_forward(
    ::AllostericEnzymeMechanism{CM,CS,RS},
    params::NamedTuple,
) where {CM,CS,RS}
    am = AllostericMechanism(AllostericEnzymeMechanism{CM,CS,RS}())
    CatN = catalytic_multiplicity(am)

    # Both conformations' polys under the free-enzyme normalization `rate_equation` uses
    # (`_mwc_state_polys`), each multiplied by the other conformation's free-enzyme weight
    # left for cross-weighting, so the saturating metabolite pattern can be matched
    # across conformations. A pattern that only exists via an `:OnlyA` step is absent
    # from the I-state groups, and a dead I-cycle leaves no I-state numerator groups.
    cat_mets = Set{Symbol}(metabolites(CM()))
    num_A, den_A, num_I, den_I, d_A, d_I = _mwc_state_polys(am, cat_mets)
    num_A_groups, den_A_groups =
        _kcat_groups_from_polys(poly_mul(num_A, d_I), poly_mul(den_A, d_I), cat_mets)
    num_I_groups, den_I_groups =
        _kcat_groups_from_polys(poly_mul(num_I, d_A), poly_mul(den_I, d_A), cat_mets)

    # kcat = peak forward turnover at saturation: max over saturating patterns
    # (met_key) and regulator corners.
    a_keys = _kcat_keys(num_A_groups, den_A_groups, Set(name.(products(reaction(am)))))

    dep, indep = _dependent_param_exprs(am)

    # Regulator corners: each ligand 0 or saturating, independently of the pattern.
    all_ligs = allosteric_regulators(am)
    lig_idx = Dict(lig => i - 1 for (i, lig) in enumerate(all_ligs))

    kcat_exprs = Any[]
    for met_key in a_keys
        # A pattern missing from the I-state numerator, or a dead I-cycle, gives A_I = 0.
        A_A, B_A = _mwc_power_pair(_poly_to_expr(num_A_groups[met_key]),
                                   _poly_to_expr(den_A_groups[met_key]), CatN)
        A_I, B_I = haskey(den_I_groups, met_key) ?
            _mwc_power_pair(_poly_to_expr(get(num_I_groups, met_key, poly_zero())),
                            _poly_to_expr(den_I_groups[met_key]), CatN) : (0, 0)
        for mask in 0:(2^length(all_ligs) - 1)
            W_A_factors = Any[]
            W_I_factors = Any[]
            deg_A = deg_I = 0
            for site in regulatory_sites(am)
                sat_terms_A = Any[]
                sat_terms_I = Any[]
                for (lig, tag) in zip(ligands(site), allo_states(site))
                    (mask >> lig_idx[lig]) & 1 == 1 || continue
                    for (terms, inactive) in ((sat_terms_A, false), (sat_terms_I, true))
                        K = _reg_K(am, site, lig, tag, inactive)
                        K === nothing || push!(terms, :(inv($K)))
                    end
                end
                n = multiplicity(site)
                if !isempty(sat_terms_A)
                    push!(W_A_factors, _power_expr(_nest_binary(:+, sat_terms_A), n))
                    deg_A += n
                end
                if !isempty(sat_terms_I)
                    push!(W_I_factors, _power_expr(_nest_binary(:+, sat_terms_I), n))
                    deg_I += n
                end
            end
            # The saturating ligands grow a conformation's regulator factor as the
            # concentration to the power `deg`. With unequal powers the conformation with
            # the larger one holds all the enzyme in the limit, so the corner is that
            # conformation's own turnover; it is skipped when that conformation lacks the
            # pattern (0/0), which is safe: taking the substrates to saturation first then
            # gives A_A/B_A, which the corner with no regulator bound already supplies.
            # Equal powers cancel and leave the weighted combination.
            if deg_A > deg_I
                push!(kcat_exprs, :($A_A / $B_A))
            elseif deg_I > deg_A
                B_I == 0 || push!(kcat_exprs, :($A_I / $B_I))
            elseif deg_A == 0
                push!(kcat_exprs, :($(_mwc_combine(A_A, A_I)) / $(_mwc_combine(B_A, B_I))))
            else
                W_A = _nest_binary(:*, W_A_factors)
                W_I = _nest_binary(:*, W_I_factors)
                push!(kcat_exprs, :(($(A_A) * $(W_A) + L * $(A_I) * $(W_I)) /
                    ($(B_A) * $(W_A) + L * $(B_I) * $(W_I))))
            end
        end
    end

    return Expr(:block,
        _destructuring_expr((indep..., :Keq), :params),
        _dep_assignments(dep)...,
        :(return $(_max_expr(kcat_exprs))))
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
    _state_parts(am, state) -> (cm, sp)

The catalytic `Mechanism` `cm` of `am` in conformational `state` (`:A` or `:I`) and
its state-tagged step constants `sp`, both read from one
`_state_allo_mechanism(am, state)`. `cm` is `_state_mechanism(am, state)`. `sp` has
the shape of `_step_parameters(cm)` — a per-flat-step vector of `Parameter`s — but
each catalytic group's `Parameter`s carry the group's state tag: a
`:NonequalAI`/`:OnlyA` group is tagged with `state`; an `:EqualAI` group is tagged
`:EqualAI` (so `name(p, am)` renders the shared bare Symbol in both states). For
`state == :I`, `:OnlyA` groups are already pruned from
`_state_allo_mechanism(am, :I)`, matching the broken-cycle graph.

Walks `_state_allo_mechanism`'s already-canonical, aligned steps/states so the
per-flat-step order matches `_flat_steps(cm)`. Each `Parameter` is anchored on its
own step (not the rep), exactly as `_step_parameters` does, so arity follows the
step's RE/SS type while `name(p, am)` collapses to the rep's structural Symbol via
the chokepoint.
"""
function _state_parts(am::AllostericMechanism, state::Symbol)
    sam = _state_allo_mechanism(am, state)
    cm = Mechanism(reaction(sam), steps(sam))
    sp = Vector{Parameter}[
        _step_constants(s, cat_allo_state(sam, g) === :EqualAI ? :EqualAI : state)
        for (g, group) in enumerate(steps(sam)) for s in group]
    @assert length(sp) == length(_flat_steps(cm)) "state step_params/steps misaligned"
    cm, sp
end

"""
The catalytic `Mechanism` for `am` in conformational `state`. `:A` keeps the
full catalytic graph; `:I` keeps only the reachable-form subgraph (`:OnlyA`
groups and the forms they disconnect from free E pruned, via
`_state_allo_mechanism`) so King–Altman re-derives the broken-cycle I-state law
natively.
"""
_state_mechanism(am::AllostericMechanism, state::Symbol) = first(_state_parts(am, state))

"""
Derive `(num_poly, den_poly, d_free_poly)` for `am`'s catalytic mechanism in
conformational `state`, natively in that state's parameter names. Runs the
shared King–Altman engine on the state-tagged `step_params` and state graph, so
no post-hoc rename is needed (`:EqualAI` groups render the shared bare Symbol
automatically). The `:I` polynomials reference each `:NonequalAI` group's
native `K_I_…`/`k_I_…` symbol; a forbidden split's `K_I_…` is defined by the
combined constraint solve's dependent assignment (`_dep_assignments`).
`d_free_poly` is free E's weight in that state's denominator (see
`_raw_symbolic_rate_polys`).
"""
_state_rate_polys(am::AllostericMechanism, state::Symbol) =
    _raw_symbolic_rate_polys(_state_parts(am, state)...)

"""
Catalytic `Parameter`s of `am` in conformation `state` (`:A` or `:I`): the
constants of each kinetic group's rep step (`Krapid` or `Kfor`+`Krev`),
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
    _dependent_param_exprs(am::AllostericMechanism)

Return `(dep_exprs, indep_params)` for an allosteric mechanism from one combined
constraint solve. Stacks the A-state and I-state thermodynamic constraint systems
(`_assemble_constraints` under each state's Wegscheider rename, tagged
`is_i_state`) over their combined column space and solves once with the shared
solver (`_solve_dependent_set`). The `(is_i_state, type)` pivot priority alone
gives the I-above-A collapse direction — an I-state column always outranks its
A-state counterpart, so a cross-state affinity tie (e.g. a live-forbidden
`:NonequalAI` split) is expressed onto the free A-side directly, with no post-hoc
merge step; within a state the catalytic pivot order (`_step_priority`) is
unchanged.

Regulator-site affinities complete no catalytic thermodynamic cycle, so they
are independent on top of the combined solve — except an `:EqualAI`
regulator, whose I-name mirrors its shared A-name (`K_I_reg = K_A_reg`, added
to `dep`). `L` (the conformational constant) is always independent. The
independents list the solve's catalytic ones by name, then the active-state
regulator constants site by site, then the inactive-state ones site by site, then
`L`. Any symbol the combined solve already made dependent is dropped from `indep`.
The `Type{<:AbstractEnzymeMechanism}` method lifts with `_concrete` and delegates
here.
"""
function _dependent_param_exprs(am::AllostericMechanism)
    function state_system(state)
        cm, sp = _state_parts(am, state)
        state_rename = _build_wegscheider_rename_map(cm; step_params = sp)
        state_rename, _assemble_constraints(cm, state_rename; step_params = sp,
                                            is_i_state = state === :I)
    end
    rename_A, (A_A, rhs_A, cols_A, pri_A) = state_system(:A)
    rename_I, (A_I, rhs_I, cols_I, pri_I) = state_system(:I)

    # Union columns: A-state first, then any I-only column (a shared `:EqualAI` group
    # carries the same bare Symbol in both states and coincides — it keeps its A tag).
    columns = unique([cols_A; cols_I])
    col = Dict(c => i for (i, c) in enumerate(columns))
    priority = merge(Dict(zip(cols_I, pri_I)), Dict(zip(cols_A, pri_A)))

    # Stack the two per-state constraint blocks over the combined column space and
    # solve once. Cross-state ties emerge as `I-row − A-row` (the `log Keq` cancels).
    A = zeros(Rational{BigInt}, size(A_A, 1) + size(A_I, 1), length(columns))
    A[axes(A_A, 1), [col[c] for c in cols_A]] = A_A
    A[size(A_A, 1) .+ axes(A_I, 1), [col[c] for c in cols_I]] = A_I
    dep, indep = _solve_dependent_set(A, [rhs_A; rhs_I], columns,
                                      [priority[c] for c in columns])

    # A per-state Wegscheider rename folds a single-symbol binding-K tie onto one
    # representative; the folded symbol enters the combined solve as a zero-column and
    # lands in `indep` as a fittable dummy. Drop it — unless a retained polynomial
    # still references it (a shared `:EqualAI` symbol can be folded in one state yet
    # used by the other state's polynomial, which keeps its own rename). The
    # non-allosteric analog needs no reference guard because it has a single,
    # consistent rename. No-op when both state renames are empty (all current specs).
    rename = merge(rename_A, rename_I)
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

"""
The dissociation-constant name of ligand `lig` (allosteric tag `tag`) at regulatory
`site` of `am` in the active (`inactive = false`) or inactive conformation, or
`nothing` when the ligand is absent from that conformation (`:OnlyI` in the active
state, `:OnlyA` in the inactive one). `:EqualAI` ligands share the A-state symbol in
both conformations (no I-state rename); `:NonequalAI` / `:OnlyI` ligands carry a
distinct I-state name. Renders the name via the `name(::Kreg, am)` chokepoint.
"""
function _reg_K(am::AllostericMechanism, site, lig, tag::Symbol, inactive::Bool)
    tag === (inactive ? :OnlyA : :OnlyI) && return nothing
    name(Kreg(site, lig, inactive && tag !== :EqualAI ? :I : :A), am)
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

"""MWC active + L·inactive state-combine `A + L * B`. Shared by `_kcat_forward`
(numerator/denominator halves of the state ratio) and `_num_den_exprs`
(the retained num/den sums)."""
_mwc_combine(a, b) = :($a + L * $b)

"""MWC binding-statistics power-pair `(X * Y^(n-1), Y^n)` for a saturating pattern:
the numerator carries one fewer denominator power than the denominator, where
`n = catalytic_multiplicity`. A zero numerator `X` stays the literal 0. Used by
`_kcat_forward` per conformation."""
_mwc_power_pair(x, y, n) =
    (n == 1 || x == 0 ? x : :($x * $y^$(n - 1)), :($y^$n))

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
(`_num_den_exprs`) and `_kcat_forward` share. The I-state polys are always re-derived
natively on the reachable-form subgraph (`_state_allo_mechanism(am, :I)` drops `:OnlyA`
groups and every form they disconnect from free E). King–Altman on that subgraph
derives the inactive conformation's binding partition, and for a dead cycle the pruned
graph's steady-state fluxes cancel exactly, so `num_I` is 0.

The free-enzyme weights `d_A` and `d_I` combine three ways, all rendering the same
value (`mets` are the catalytic metabolites):
- `d_A == d_I`: the polys stay raw (identical conformations; the factor cancels);
- both metabolite-free monomials: each state's polys are divided by its own weight
  (the clean standard-MWC form);
- otherwise: the polys stay raw, and each state is cross-weighted by the other state's
  weight (`d^n` on the rate law's terms).
Only the third case returns the weights; the first two return `poly_one()` for both.
The normalization is a common factor of kcat's saturating-limit ratio, so it leaves
kcat's value unchanged; sharing it keeps kcat's saturating-pattern grouping consistent
with `rate_equation`.
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
Returns `(full_num, full_den)`. Per-active-site normalization: the
numerator carries no leading `catalytic_multiplicity` factor; only the
`Q_cat^(CatN-1)` / `Q_cat^CatN` binding-statistics powers remain.
"""
function _num_den_exprs(@nospecialize(M_type::Type{<:AllostericEnzymeMechanism}))
    m = M_type()
    am = AllostericMechanism(m)
    CatN = catalytic_multiplicity(m)
    RS = regulatory_sites(am)

    # A-state catalytic param symbols (the tagged column set) drive `_poly_to_expr`'s
    # param/metabolite ordering split; the I-poly's `:I` symbols sort as non-params.
    cat_params = Set(_param_columns(_state_parts(am, :A)...))
    num_A_poly, den_A_poly, num_i_poly, den_i_poly, d_A, d_I =
        _mwc_state_polys(am, Set{Symbol}(metabolites(catalytic_mechanism(m))))
    # A free-enzyme weight left for cross-weighting multiplies the other state's terms;
    # a normalized one renders as 1, which `_mwc_cross_weight` skips.
    D_A_expr = _poly_to_expr(d_A, cat_params)
    D_I_expr = _poly_to_expr(d_I, cat_params)

    N_A = _poly_to_expr(num_A_poly, cat_params)
    Q_A = _poly_to_expr(den_A_poly, cat_params)
    N_I = _poly_to_expr(num_i_poly, cat_params)
    Q_I = _poly_to_expr(den_i_poly, cat_params)

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
