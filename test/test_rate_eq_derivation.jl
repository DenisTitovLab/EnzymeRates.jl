# Tests for enzyme rate equation derivation
# Validates structure, constraints,
# and correctness of derived rate equations

using OrdinaryDiffEqFIRK

# ── Helper functions ────────────────────────────────────────────────────────

# ── Reference QSSA implementation ───────────────────────────────────────────

"""
Independent reference: compute QSSA rate using Laplacian cofactor method.
Works directly with EnzymeMechanism type parameters.
"""
function reference_qssa(
    m::EnzymeMechanism,
    params::NamedTuple,
    concs::NamedTuple,
)
    # Walk reactions via accessor so this helper works for both legacy
    # `(mets, rxns)` Sig shapes and the new `(reaction_sig, steps_sig)`
    # shape produced by `EnzymeMechanism(::Mechanism)`.
    Reactions = EnzymeRates.reactions(m)
    enz_names = EnzymeRates.enzyme_forms(m)
    n = length(enz_names)
    name_to_idx = Dict(nm => i for (i, nm) in enumerate(enz_names))
    enz_set = Set(enz_names)

    ref_name, nu_ref = _reference_metabolite(m)

    # Build rate matrix R[i,j] = pseudo-first-order rate from i to j
    R = zeros(n, n)
    for (step_idx, (lhs, rhs)) in enumerate(Reactions)
        e_lhs = first(s for s in lhs if s in enz_set)
        e_rhs = first(s for s in rhs if s in enz_set)
        i = name_to_idx[e_lhs]
        j = name_to_idx[e_rhs]

        m_lhs = [s for s in lhs if s ∉ enz_set]
        m_rhs = [s for s in rhs if s ∉ enz_set]

        kf = params[Symbol("k$(step_idx)f")]
        kr = params[Symbol("k$(step_idx)r")]

        rf = isempty(m_lhs) ? kf : kf * concs[m_lhs[1]]
        R[i, j] += rf

        rr = isempty(m_rhs) ? kr : kr * concs[m_rhs[1]]
        R[j, i] += rr
    end

    # Build Laplacian
    L = zeros(n, n)
    for i in 1:n
        for j in 1:n
            if i != j
                L[i, j] = -R[i, j]
                L[i, i] += R[i, j]
            end
        end
    end

    # Cofactors
    D = zeros(n)
    for i in 1:n
        rows = [r for r in 1:n if r != i]
        cols = [c for c in 1:n if c != i]
        D[i] = det(L[rows, cols])
    end

    D_total = sum(D)
    E_conc = D ./ D_total .* params.E_total

    # Compute net consumption of reference substrate
    v = 0.0
    for (step_idx, (lhs, rhs)) in enumerate(Reactions)
        e_lhs = first(s for s in lhs if s in enz_set)
        e_rhs = first(s for s in rhs if s in enz_set)
        i = name_to_idx[e_lhs]
        j = name_to_idx[e_rhs]

        m_lhs = [s for s in lhs if s ∉ enz_set]
        m_rhs = [s for s in rhs if s ∉ enz_set]

        kf = params[Symbol("k$(step_idx)f")]
        kr = params[Symbol("k$(step_idx)r")]

        rf = kf * (isempty(m_lhs) ? 1.0 : concs[m_lhs[1]])
        rr = kr * (isempty(m_rhs) ? 1.0 : concs[m_rhs[1]])

        flux = rf * E_conc[i] - rr * E_conc[j]
        if !isempty(m_lhs) && m_lhs[1] == ref_name
            v += flux
        elseif !isempty(m_rhs) && m_rhs[1] == ref_name
            v -= flux
        end
    end

    return v / abs(nu_ref)
end

# ── Parameter generation helpers ────────────────────────────────────────────

"""
Re-key a structural-named parameter NamedTuple to the **per-step** positional
names (`K1`, `k1f`, `k1r`, …) that hand-derived analytical oracles destructure.
Walks every step in source order; emits one positional name per RE step or two
per SS step keyed on the step's source position; each value is looked up by the
rep parameter's structural name in `nt` (so group members share the rep's
value, matching how the package's destructure-by-name works today).
Permanent test utility — oracles are inherently positional and source-indexed.
"""
function positional_params(m, nt::NamedTuple)
    mech = m isa EnzymeRates.Mechanism             ? m :
           m isa EnzymeRates.AllostericMechanism   ? m :
           m isa EnzymeRates.AllostericEnzymeMechanism ?
               EnzymeRates.AllostericMechanism(m) :
           EnzymeRates.Mechanism(m)
    names = Symbol[]
    vals  = Any[]
    for group in EnzymeRates.steps(mech)
        rep = first(group)
        if EnzymeRates.is_equilibrium(rep)
            rep_name = EnzymeRates.name(EnzymeRates.Kd(rep, :None), mech)
            for s in group
                push!(names, Symbol("K", EnzymeRates.source_idx(s)))
                push!(vals,  nt[rep_name])
            end
        else
            fwd_name = EnzymeRates.name(EnzymeRates.Kfor(rep, :None), mech)
            rev_name = EnzymeRates.name(EnzymeRates.Krev(rep, :None), mech)
            for s in group
                idx = EnzymeRates.source_idx(s)
                push!(names, Symbol("k", idx, "f")); push!(vals, nt[fwd_name])
                push!(names, Symbol("k", idx, "r")); push!(vals, nt[rev_name])
            end
        end
    end
    # Pass through any remaining params not covered by the step loop (T-state
    # params, Lallo, Kreg, Keq, E_total, etc.) under their existing names.
    emitted = Set(names)
    for k in keys(nt)
        k ∈ emitted && continue
        push!(names, k)
        push!(vals, nt[k])
    end
    NamedTuple{Tuple(names)}(Tuple(vals))
end

"""Generate random reduced (fitted) params + Keq + E_total for a mechanism."""
function random_reduced_params(m; rng=Random.default_rng())
    fp = EnzymeRates.fitted_params(m)
    vals = Tuple(0.1 + 9.9 * rand(rng) for _ in fp)
    Keq_val = 0.1 + 9.9 * rand(rng)
    E_total_val = 0.1 + 9.9 * rand(rng)
    keys_out = (fp..., :Keq, :E_total)
    vals_out = (vals..., Keq_val, E_total_val)
    NamedTuple{keys_out}(vals_out)
end

"""Check if mechanism has any rapid-equilibrium steps."""
_has_re_steps(m) = any(EnzymeRates.equilibrium_steps(m))

"""
Test that `rate_equation` is non-allocating and fast for the given mechanism.
Must be a standalone function to avoid @testset closure boxing.
"""
function test_rate_equation_performance(m, params, concs)
    rate_equation(m, concs, params) # warmup/compile
    allocs = @allocated rate_equation(m, concs, params)
    t = @elapsed for _ in 1:10_000; rate_equation(m, concs, params); end
    return allocs, t / 10_000
end

"""
Get independent parameter symbols from mechanism using internal API.
"""
function _get_independent_params(m)
    _, indep = EnzymeRates._dependent_param_exprs(typeof(m))
    return indep
end

"""
Get dependent parameter expressions from mechanism using internal API.
Returns vector of (symbol, expression_string) pairs.
"""
function _get_dependent_params(m)
    dep_exprs, _ = EnzymeRates._dependent_param_exprs(typeof(m))
    pairs = Tuple{Symbol, String}[]
    for (sym, expr) in sort(collect(dep_exprs); by=first)
        push!(pairs, (sym, string(expr)))
    end
    return pairs
end

"""
Build a NamedTuple with only independent params + Keq + E_total,
given all_params (with all k's + E_total) and a Keq value.
"""
function make_independent_params(m, all_params, Keq)
    indep = _get_independent_params(m)
    keys_out = (indep..., :Keq, :E_total)
    vals_out = Tuple(
        k == :Keq ? Keq :
        k == :E_total ? all_params.E_total :
        all_params[k]
        for k in keys_out
    )
    return NamedTuple{keys_out}(vals_out)
end

"""
Compute all k values from independent params + Keq by
evaluating dependent parameter expressions.
Returns a NamedTuple with all k's + E_total.
"""
function compute_all_params(m, new_params)
    eq = EnzymeRates.equilibrium_steps(m)
    ns = EnzymeRates.n_steps(m)
    dep = _get_dependent_params(m)
    # Build all raw param keys: K_i for RE steps, k_if/k_ir for SS steps
    all_keys = Symbol[]
    for i in 1:ns
        if eq[i]
            push!(all_keys, Symbol("K$i"))
        else
            push!(all_keys, Symbol("k$(i)f"))
            push!(all_keys, Symbol("k$(i)r"))
        end
    end
    push!(all_keys, :E_total)
    # Evaluate dependent expressions using new_params as the "params" namespace
    dep_dict = Dict{Symbol, Float64}()
    for (sym, expr_str) in dep
        val = _eval_dep_expr(expr_str, new_params)
        dep_dict[sym] = val
    end
    # Resolve kinetic-group aliases (e.g., K2 → K1 when steps share group)
    rename = EnzymeRates._build_kinetic_rename_map(m)
    for (alias, rep) in rename
        if haskey(dep_dict, rep)
            dep_dict[alias] = dep_dict[rep]
        elseif haskey(new_params, rep)
            dep_dict[alias] = Float64(new_params[rep])
        end
    end
    all_vals = Float64[haskey(dep_dict, k) ? dep_dict[k] :
                       haskey(new_params, k) ? Float64(new_params[k]) :
                       error("Missing parameter $k")
                       for k in all_keys]
    return NamedTuple{Tuple(all_keys)}(Tuple(all_vals))
end

"""
AllostericEnzymeMechanism version of compute_all_params.
Returns all independent + dependent (Haldane-derived) params + Keq + E_total.
"""
function compute_all_params(m::EnzymeRates.AllostericEnzymeMechanism, new_params)
    indep = _get_independent_params(m)
    dep = _get_dependent_params(m)
    dep_dict = Dict{Symbol, Float64}()
    for (sym, expr_str) in dep
        dep_dict[sym] = _eval_dep_expr(expr_str, new_params)
    end
    all_keys = (indep..., Tuple(keys(dep_dict))..., :Keq, :E_total)
    all_vals = Tuple(Float64[
        haskey(dep_dict, k) ? dep_dict[k] :
        k == :Keq ? Float64(new_params[:Keq]) :
        k == :E_total ? Float64(new_params[:E_total]) :
        Float64(new_params[k])
        for k in all_keys
    ])
    NamedTuple{all_keys}(all_vals)
end

function _eval_dep_expr(expr_str::String, params::NamedTuple)
    # Build let bindings from params - bind each key directly so bare names work
    bindings = ["$(k) = $(params[k])" for k in keys(params)]
    code = "let $(join(bindings, ", "))\n  $expr_str\nend"
    return Float64(eval(Meta.parse(code)))
end

"""
Generate random independent params + Keq + E_total for testing.
Also returns all_params (old-style with all k's + E_total) for reference comparison.
"""
function random_independent_params_concs(
    m, met_names::Vector{Symbol}; rng=Random.default_rng()
)
    indep = _get_independent_params(m)
    # Generate random values for independent params + Keq + E_total
    param_keys = (indep..., :Keq, :E_total)
    param_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in param_keys)
    new_params = NamedTuple{param_keys}(param_vals)

    conc_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in met_names)
    concs = NamedTuple{Tuple(met_names)}(conc_vals)

    all_params = compute_all_params(m, new_params)
    return new_params, concs, all_params
end

"""
Convert raw params (K_g / k_g_f / k_g_r for each kinetic group's
representative step) to ODE params (large k_if/k_ir for all steps).
Steps sharing a kinetic_group share the same K (or k_f/k_r) param.

For binding RE steps (metabolite on LHS, canonical form):
    K = Kd = kr/kf, so k_if = 1e6, k_ir = 1e6 * K.
For RE isomerization steps (no metabolite, enzyme-only):
    K = Ka = kf/kr, so k_if = 1e6 * K, k_ir = 1e6.
"""
function raw_to_ode_params(m, raw_params)
    eq = EnzymeRates.equilibrium_steps(m)
    ns = EnzymeRates.n_steps(m)
    rxns = EnzymeRates.reactions(m)
    enz_set = Set(EnzymeRates.enzyme_forms(m))
    # A canonical RE binding step has a metabolite on LHS (canonical form
    # invariant: all RE binding steps are written `E + S ⇌ ES`).
    is_binding_step = Bool[
        eq[i] && any(s ∉ enz_set for s in rxns[i][1])
        for i in 1:ns
    ]
    param_keys = Symbol[]
    param_vals = Float64[]
    for i in 1:ns
        g = EnzymeRates.kinetic_group(m, i)
        rep = first(EnzymeRates.steps_in_group(m, g))
        push!(param_keys, Symbol("k$(i)f"))
        push!(param_keys, Symbol("k$(i)r"))
        if eq[i]
            K = Float64(raw_params[Symbol("K$rep")])
            if is_binding_step[i]
                # Binding step (metabolite on LHS): K = Kd = kr/kf
                push!(param_vals, 1e6)
                push!(param_vals, 1e6 * K)
            else
                # RE isomerization (no metabolite): K = Ka = kf/kr
                push!(param_vals, 1e6 * K)
                push!(param_vals, 1e6)
            end
        else
            push!(param_vals, Float64(raw_params[Symbol("k$(rep)f")]))
            push!(param_vals, Float64(raw_params[Symbol("k$(rep)r")]))
        end
    end
    push!(param_keys, :E_total)
    push!(param_vals, Float64(raw_params[:E_total]))
    return NamedTuple{Tuple(param_keys)}(Tuple(param_vals))
end

function _reference_metabolite(m)
    subs = EnzymeRates.substrates(m)
    isempty(subs) && error("No substrate found in mechanism")
    name = subs[1]
    coeff = -count(==(name), subs)
    return name, coeff
end

# ── ODE steady-state helpers ────────────────────────────────────────────────

function build_ode_rhs(
    m::EnzymeMechanism,
    params, concs,
)
    # Walk reactions via accessor so this helper works for both legacy
    # and new Sig shapes.
    Reactions = EnzymeRates.reactions(m)
    enz_names = EnzymeRates.enzyme_forms(m)
    name_to_idx = Dict(nm => i for (i, nm) in enumerate(enz_names))
    enz_set = Set(enz_names)

    step_data = []
    for (step_idx, (lhs, rhs)) in enumerate(Reactions)
        e_lhs = first(s for s in lhs if s in enz_set)
        e_rhs = first(s for s in rhs if s in enz_set)
        i = name_to_idx[e_lhs]
        j = name_to_idx[e_rhs]

        m_lhs = [s for s in lhs if s ∉ enz_set]
        m_rhs = [s for s in rhs if s ∉ enz_set]

        kf = Float64(params[Symbol("k$(step_idx)f")])
        kr = Float64(params[Symbol("k$(step_idx)r")])

        rf = isempty(m_lhs) ? kf : kf * concs[m_lhs[1]]
        rr = isempty(m_rhs) ? kr : kr * concs[m_rhs[1]]

        push!(step_data, (i, j, rf, rr))
    end

    n = length(enz_names)
    function rhs!(du, u, p, t)
        fill!(du, 0.0)
        for (i, j, rf, rr) in step_data
            flux = rf * u[i] - rr * u[j]
            du[i] -= flux
            du[j] += flux
        end
    end
    return rhs!
end

function ode_steady_state_flux(
    m::EnzymeMechanism,
    params, concs,
)
    Reactions = EnzymeRates.reactions(m)
    E_total = params.E_total
    enz_names = EnzymeRates.enzyme_forms(m)
    n = length(enz_names)
    enz_set = Set(enz_names)
    ref_name, nu_ref = _reference_metabolite(m)

    u0 = zeros(n)
    u0[1] = E_total

    rhs! = build_ode_rhs(m, params, concs)
    prob = ODEProblem(rhs!, u0, (0.0, 1e6))
    sol = solve(prob, RadauIIA9(); abstol=1e-12, reltol=1e-12)
    u_ss = sol.u[end]

    name_to_idx = Dict(nm => i for (i, nm) in enumerate(enz_names))
    v = 0.0
    for (step_idx, (lhs, rhs)) in enumerate(Reactions)
        e_lhs = first(s for s in lhs if s in enz_set)
        e_rhs = first(s for s in rhs if s in enz_set)
        i = name_to_idx[e_lhs]
        j = name_to_idx[e_rhs]

        m_lhs = [s for s in lhs if s ∉ enz_set]
        m_rhs = [s for s in rhs if s ∉ enz_set]

        kf = Float64(params[Symbol("k$(step_idx)f")])
        kr = Float64(params[Symbol("k$(step_idx)r")])
        rf = isempty(m_lhs) ? kf : kf * concs[m_lhs[1]]
        rr = isempty(m_rhs) ? kr : kr * concs[m_rhs[1]]

        flux = rf * u_ss[i] - rr * u_ss[j]
        if !isempty(m_lhs) && m_lhs[1] == ref_name
            v += flux
        elseif !isempty(m_rhs) && m_rhs[1] == ref_name
            v -= flux
        end
    end

    return v / abs(nu_ref)
end

# ── Rate equation string evaluation helper ──────────────────────────────────

function _eval_rate_string(s, params, concs)
    # Filter out destructuring lines ("= params", "= concs"), keep constraint + v = lines
    lines = split(s, "\n")
    eval_lines = String[]
    for line in lines
        stripped = strip(line)
        isempty(stripped) && continue
        endswith(stripped, "= params") && continue
        endswith(stripped, "= concs") && continue
        push!(eval_lines, stripped)
    end
    # The last line should be "v = ..." — extract the expression after "v = "
    code_body = join(eval_lines, "\n")
    eq_line = last(split(code_body, "v = "; limit=2))
    bindings = vcat(
        ["$k = $(params[k])" for k in keys(params)],
        ["$k = $(concs[k])" for k in keys(concs)],
    )
    code = "let $(join(bindings, ", "))\n  $eq_line\nend"
    eval(Meta.parse(code))
end

# ── Modular test functions for MechanismTestSpec ────────────────────────────

function test_structure(spec::MechanismTestSpec)
    m = spec.mechanism
    @testset "Structure" begin
        @test EnzymeRates.n_states(m) == spec.expected_n_states
        @test EnzymeRates.n_steps(m) == spec.expected_n_steps
        @test length(metabolites(m)) == spec.expected_n_metabolites
    end
end

"""Classify a dep expression as Haldane (RHS references Keq), Mirror
(RHS is a single Symbol), or Wegscheider (RHS Expr without Keq)."""
function _classify_dep_expr(expr)
    if expr isa Symbol
        return :mirror
    elseif EnzymeRates._expr_references_any(expr, Set([:Keq]))
        return :haldane
    else
        return :wegscheider
    end
end

function test_constraint_counting(spec::MechanismTestSpec)
    m = spec.mechanism
    @testset "Constraints" begin
        dep_exprs, indep = EnzymeRates._dependent_param_exprs(typeof(m))
        n_haldane = 0
        n_mirror = 0
        n_wegscheider = 0
        for (_, expr) in dep_exprs
            cat = _classify_dep_expr(expr)
            if cat == :haldane
                n_haldane += 1
            elseif cat == :mirror
                n_mirror += 1
            else
                n_wegscheider += 1
            end
        end
        @test n_haldane == spec.expected_n_haldane_constraints
        @test n_mirror == spec.expected_n_mirror_constraints
        @test n_wegscheider == spec.expected_n_wegscheider_constraints
        @test length(indep) == spec.expected_n_independent_params
    end
end

function test_reference_qssa(spec::MechanismTestSpec; n_trials=20, seed=42)
    m = spec.mechanism
    # Reference QSSA only works for all-SS mechanisms
    _has_re_steps(m) && return
    met_names = spec.metabolite_names
    @testset "Reference QSSA" begin
        rng = Random.MersenneTwister(seed)
        @test all(1:n_trials) do _
            new_params, concs, all_params =
                random_independent_params_concs(
                    m, met_names; rng=rng)
            isapprox(
                rate_equation(m, concs, new_params),
                reference_qssa(m, all_params, concs);
                rtol=spec.reference_rtol)
        end
    end
end

function test_analytical_rate(spec::MechanismTestSpec; n_trials=20, seed=1001)
    # Skip if no analytical rate function provided
    spec.analytical_rate_fn === nothing && return

    m = spec.mechanism
    met_names = spec.metabolite_names
    @testset "Analytical Rate" begin
        rng = Random.MersenneTwister(seed)
        @test all(1:n_trials) do _
            new_params, concs, all_params =
                random_independent_params_concs(
                    m, met_names; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(positional_params(m, all_params), (Et=Et,))
            p_pkg = merge(new_params, (E_total=Et,))
            isapprox(
                rate_equation(m, concs, p_pkg),
                spec.analytical_rate_fn(p, concs);
                rtol=1e-10)
        end
    end
end

function test_haldane_equilibrium(spec::MechanismTestSpec; seed=42)
    m = spec.mechanism
    met_names = spec.metabolite_names
    @testset "Haldane Equilibrium" begin
        rng = Random.MersenneTwister(seed)
        new_params, _, _ = random_independent_params_concs(m, met_names; rng=rng)
        Keq = new_params.Keq
        n_prods = length(EnzymeRates.products(m))
        # Build equilibrium concentrations: prod(P_i) / prod(S_i) = Keq
        # Set all substrates to 1.0, distribute Keq^(1/n_prods) across products
        sub_names = collect(EnzymeRates.substrates(m))
        prod_names = collect(EnzymeRates.products(m))
        eq_vals = Dict{Symbol,Float64}()
        for s in sub_names; eq_vals[s] = 1.0; end
        p_each = Keq^(1.0 / n_prods)
        for p in prod_names; eq_vals[p] = p_each; end
        # Regulators don't affect equilibrium — set to arbitrary value
        for s in met_names
            haskey(eq_vals, s) || (eq_vals[s] = 2.5)
        end
        eq_concs = NamedTuple{Tuple(met_names)}(Tuple(eq_vals[s] for s in met_names))
        v_eq = rate_equation(m, eq_concs, new_params)
        @test abs(v_eq) < 1e-10
    end
end

function test_performance(spec::MechanismTestSpec; seed=42)
    m = spec.mechanism
    met_names = spec.metabolite_names
    @testset "Performance" begin
        rng = Random.MersenneTwister(seed)
        params, concs, _ = random_independent_params_concs(m, met_names; rng=rng)
        allocs, t = test_rate_equation_performance(m, params, concs)
        @test allocs == 0
        @test t < 100e-9
    end
end

function test_ode_steadystate(spec::MechanismTestSpec; n_trials=10, seed=42)
    m = spec.mechanism
    met_names = spec.metabolite_names
    @testset "ODE Steady-State" begin
        rng = Random.MersenneTwister(seed)
        has_re = _has_re_steps(m)
        @test all(1:n_trials) do _
            new_params, concs, all_params =
                random_independent_params_concs(
                    m, met_names; rng=rng)
            # Convert K_i to large k_if/k_ir for ODE
            ode_params = has_re ?
                raw_to_ode_params(m, all_params) :
                all_params
            v_ode = ode_steady_state_flux(m, ode_params, concs)
            v_ka = rate_equation(m, concs, new_params)
            # Use looser tolerance for RE mechanisms (large rate approximation)
            rtol = has_re ? 1e-3 : spec.ode_rtol
            isapprox(v_ode, v_ka; rtol=rtol)
        end
    end
end

function test_rate_equation_string(spec::MechanismTestSpec)
    m = spec.mechanism
    met_names = spec.metabolite_names
    @testset "Rate Equation String" begin
        s = rate_equation_string(m)
        @test !occursin("params.", s)   # No field-access prefixes
        @test !occursin("concs.", s)    # No field-access prefixes
        @test !occursin("+ -", s)       # No malformed signs
        @test !occursin("- -", s)       # No malformed signs
        @test occursin("v = E_total * (", s)  # Proper format
        @test occursin(") / (", s)      # Has denominator
        @test occursin("= params", s)   # Has params destructuring
        @test occursin("= concs", s)    # Has concs destructuring

        # Numerical equivalence test
        rng = Random.MersenneTwister(9000 + hash(spec.name) % 1000)
        @test all(1:10) do _
            new_params, concs, all_params =
                random_independent_params_concs(
                    m, met_names; rng=rng)
            isapprox(
                rate_equation(m, concs, new_params),
                _eval_rate_string(s, all_params, concs);
                rtol=1e-10)
        end
    end
end

"""
Extract numerator and denominator strings from rate equation v-line.
Handles nested parentheses via depth counting.
"""
function _extract_num_denom(v_line::AbstractString)
    marker = "E_total * ("
    start = findfirst(marker, v_line)
    start === nothing && return nothing, nothing
    pos = last(start) + 1
    depth = 1
    while depth > 0 && pos <= length(v_line)
        c = v_line[pos]
        depth += (c == '(') - (c == ')')
        pos += 1
    end
    num_str = v_line[last(start)+1:pos-2]
    rest = v_line[pos:end]
    div_start = findfirst(" / (", rest)
    div_start === nothing && return String(num_str), nothing
    denom_begin = pos + last(div_start)
    denom_str = v_line[denom_begin:end-1]
    return String(num_str), String(denom_str)
end

function test_factored_form(spec::MechanismTestSpec)
    has_num = spec.expected_factored_num !== nothing
    has_denom = spec.expected_factored_denom !== nothing
    (has_num || has_denom) || return
    m = spec.mechanism
    @testset "Factored Form" begin
        s = rate_equation_string(m)
        v_line = last(split(s, "\n"))
        num_str, denom_str = _extract_num_denom(v_line)
        @test num_str !== nothing
        @test denom_str !== nothing
        if num_str !== nothing && denom_str !== nothing
            has_num && @test num_str == spec.expected_factored_num
            has_denom && @test denom_str == spec.expected_factored_denom
        end
    end
end

function test_analytical_kcat(spec::MechanismTestSpec; seed=42)
    spec.analytical_kcat_fn === nothing && return
    m = spec.mechanism
    @testset "Analytical kcat" begin
        rng = Random.MersenneTwister(seed)
        params = random_reduced_params(m; rng)
        kcat = EnzymeRates._kcat_forward(m, params)
        @test kcat ≈ spec.analytical_kcat_fn(params) rtol=1e-10
    end
end

function test_kcat_rescaling(spec::MechanismTestSpec; seed=100)
    m = spec.mechanism
    @testset "kcat rescaling" begin
        rng = Random.MersenneTwister(seed)
        params = random_reduced_params(m; rng)

        # kcat should be positive
        kcat_orig = EnzymeRates._kcat_forward(m, params)
        @test kcat_orig > 0

        # Scale invariance: scaling SS k's by α scales kcat by α
        ss_names = EnzymeRates._ss_rate_constant_names(m)
        α = 0.1 + 9.9 * rand(rng)
        scaled_params = NamedTuple{keys(params)}(Tuple(
            k in ss_names ? v * α : v
            for (k, v) in zip(keys(params), values(params))
        ))
        @test EnzymeRates._kcat_forward(m, scaled_params) ≈
            α * kcat_orig rtol=1e-10

        # Rescale so kcat = 1
        norm = rescale_parameter_values(m, params)
        kcat_norm = EnzymeRates._kcat_forward(m, norm)
        @test kcat_norm ≈ 1.0 rtol=1e-10

        # K values (non-SS params) unchanged
        for k in keys(params)
            if !(k in ss_names)
                @test norm[k] == params[k]
            end
        end

        # Custom kcat target
        kcat_target = 0.1 + 9.9 * rand(rng)
        norm_custom = rescale_parameter_values(m, params; kcat=kcat_target)
        @test EnzymeRates._kcat_forward(m, norm_custom) ≈
            kcat_target rtol=1e-10

        # Rate proportionality: v_norm / v_orig = 1 / kcat_orig
        met_names = metabolites(m)
        conc_vals = Tuple(0.5 + rand(rng) for _ in met_names)
        concs = NamedTuple{Tuple(met_names)}(conc_vals)
        v_orig = rate_equation(m, concs, params)
        v_norm = rate_equation(m, concs, norm)
        @test v_norm / v_orig ≈ 1.0 / kcat_orig rtol=1e-8

        # V ≈ 1 at saturating substrates, products=0
        sub_names = collect(EnzymeRates.substrates(m))
        prod_names = collect(EnzymeRates.products(m))
        reg_names = collect(EnzymeRates.regulators(m))
        n_reg = length(reg_names)

        norm_e1 = merge(norm, (E_total=1.0,))

        BIG = 1e6
        max_rate = 0.0
        for mask in 0:(2^n_reg - 1)
            conc_dict = Dict{Symbol,Float64}()
            for s in sub_names; conc_dict[s] = BIG; end
            for p in prod_names; conc_dict[p] = 0.0; end
            for (i, r) in enumerate(reg_names)
                conc_dict[r] = ((mask >> (i - 1)) & 1) == 1 ? BIG : 0.0
            end
            concs = NamedTuple{Tuple(met_names)}(
                Tuple(conc_dict[n] for n in met_names))
            v = rate_equation(m, concs, norm_e1)
            max_rate = max(max_rate, v)
        end
        @test max_rate ≈ 1.0 rtol=1e-3
    end
end

"""
Run all tests for a mechanism specification.
Organizes tests by mechanism: all tests for one mechanism together.
"""
function run_all_tests(spec::MechanismTestSpec)
    @testset "$(spec.name)" begin
        test_structure(spec)
        test_constraint_counting(spec)
        test_reference_qssa(spec)
        test_analytical_rate(spec)      # Only runs if analytical_rate_fn provided
        test_haldane_equilibrium(spec)
        test_performance(spec)
        test_rate_equation_string(spec)
        test_factored_form(spec)        # Only runs if expected strings provided
        spec.run_ode_test && test_ode_steadystate(spec)
        test_analytical_kcat(spec)      # Only runs if analytical_kcat_fn provided
        test_kcat_rescaling(spec)
    end
end

# ── Main test loop ──────────────────────────────────────────────────────────

@testset "Enzyme Derivation Tests" begin
    for spec in MECHANISM_TEST_SPECS
        run_all_tests(spec)
    end
end

# ── Standalone kcat tests ──────────────────────────────────────────────────────

@testset "rate_equation polynomial body uses 2-arg +/* calls" begin
    # The fitter calls rate_equation millions of times per CV fold. The
    # polynomial body emitted by _poly_to_expr (via _nest_binary) MUST
    # have exactly 2 operands per +/* call so LLVM inlines the binary
    # Float64 path; n-ary varargs above ~30 terms boxes the argument
    # tuple and turns 100ns/0B into 1µs/2KB per call. See
    # docs/superpowers/specs/2026-05-11-rate-eq-emission-perf-fix-design.md.
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "Random-order Bi-Bi")
    rate_expr, _, _ = EnzymeRates._raw_rate_expr_and_symbols(
        typeof(spec.mechanism))
    bad = Expr[]
    function walk!(e)
        if e isa Expr
            if e.head == :call && !isempty(e.args) &&
               e.args[1] isa Symbol && e.args[1] in (:+, :*) &&
               length(e.args) != 3
                push!(bad, e)
            end
            start = e.head == :call ? 2 : 1
            for i in start:length(e.args)
                walk!(e.args[i])
            end
        end
    end
    walk!(rate_expr)
    @test isempty(bad)
end

@testset "rate_equation_string prints flat +/* sums (no nested parens)" begin
    # Orthogonal guard (NOT a perf-fix backstop): _expr_to_string is
    # precedence-aware and flattens nested +/* nodes
    # transparently: balanced +(+(a,b), +(c,d)) prints as "a + b + c + d"
    # with no added parens. If that flattening regresses, the printed
    # output would gain nested parens but still evaluate correctly —
    # easy to miss without an explicit guard. Witness count for
    # Random-order Bi-Bi is 5 in Full mode (params destructure, concs
    # destructure, num wrap, den wrap, one negatives wrap inside num).
    # Full mode is the right witness here: Reduced mode adds dep-expr
    # 1/(…) divisor assignments and lands at ~13 parens, which is
    # structural noise that would mask a flattening regression.
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "Random-order Bi-Bi")
    s = rate_equation_string(spec.mechanism, EnzymeRates.Full)
    @test count(==('('), s) <= 6
end

@testset "_ss_rate_constant_names" begin
    # SS-only Uni-Uni: 2 SS binding steps → k1f, k1r, k2f, k2r are SS.
    uni_uni = only(s for s in MECHANISM_TEST_SPECS
                   if s.name == "Uni-Uni").mechanism
    names = EnzymeRates._ss_rate_constant_names(uni_uni)
    for sym in (:k1f, :k1r, :k2f, :k2r)
        @test sym in names
    end

    # Mixed RE/SS: RE binding (K1) + SS catalysis (k2f, k2r). Only the
    # SS k's are returned; RE binding K is excluded.
    re_uu = only(s for s in MECHANISM_TEST_SPECS
                 if s.name == "RE Uni-Uni").mechanism
    re_uu_names = EnzymeRates._ss_rate_constant_names(re_uu)
    @test :k2f in re_uu_names && :k2r in re_uu_names
    for sym in (:K1, :Keq, :L, :E_total)
        @test !(sym in re_uu_names)
    end

    # Allosteric: T-state versions of every SS rate constant are also
    # included so `rescale_parameter_values` scales them in tandem.
    mwc = only(s for s in MECHANISM_TEST_SPECS
               if s.name == "MWC Dimer [AllostericEnzymeMechanism]").mechanism
    mwc_names = EnzymeRates._ss_rate_constant_names(mwc)
    for sym in (:k3f, :k3r, :k3f_T, :k3r_T)
        @test sym in mwc_names
    end
    for sym in (:K1, :K2, :K1_T, :K2_T, :Keq, :L, :E_total)
        @test !(sym in mwc_names)
    end
end

# ── Degenerate constraint handling ────────────────────────────────────────────

@testset "Degenerate constraint handling" begin
    # ── Unit tests: build_power_expr return types ─────────────────

    @testset "build_power_expr return types" begin
        bpe = EnzymeRates.build_power_expr
        R = Rational{BigInt}
        # Single symbol factor (exp=1): returns bare Symbol
        @test bpe(R(0), [(:k1f, R(1))]) === :k1f

        # Single inverse factor: returns Expr
        @test bpe(R(0), [(:k1f, R(-1))]) isa Expr

        # Keq-only: returns Symbol
        @test bpe(R(1), Tuple{Symbol, R}[]) === :Keq

        # Multiple factors: returns Expr
        @test bpe(R(0), [(:k1f, R(1)), (:k2f, R(1))]) isa Expr

        # No factors and zero Keq: returns Int literal 1
        @test bpe(R(0), Tuple{Symbol, R}[]) === 1

        # All return types are valid AST nodes
        for r in [
            bpe(R(0), [(:k1f, R(1))]),
            bpe(R(0), [(:k1f, R(-1))]),
            bpe(R(1), Tuple{Symbol, R}[]),
            bpe(R(0), [(:k1f, R(1)), (:k2f, R(1))]),
            bpe(R(1), [(:k1f, R(1))]),
        ]
            @test r isa Union{Int, Symbol, Expr}
        end
    end

end

# ── Large equation compilation regression test ────────────────────────────

@testset "Large equation compilation (<20s)" begin
    # Use the manually-defined large mechanism (11 forms, 16 steps)
    # from the "Rate equation too large error" test below, but with
    # all steps as RE to keep it compilable.
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        competitive_inhibitors: R1
    end
    # Add dead-end regulator R1 to the largest topology
    n_forms_of(mech) = length(Set(
        sp for grp in EnzymeRates.steps(mech) for st in grp
        for sp in (EnzymeRates.from_species(st),
                   EnzymeRates.to_species(st))))
    topos = EnzymeRates.init_mechanisms(rxn)
    sort!(topos; by=n_forms_of, rev=true)
    variants = EnzymeRates._expand_add_dead_end_regulator(topos[1], rxn)
    # Find largest compilable mechanism
    sort!(variants; by=n_forms_of, rev=true)
    m = nothing
    for s in variants
        try
            m = EnzymeMechanism(s)
            parameters(m)
            break
        catch
            m = nothing
        end
    end
    @test m !== nothing

    metabs = metabolites(m)
    params_tup = parameters(m)
    concs = NamedTuple{metabs}(ones(length(metabs)))
    pvals = NamedTuple{params_tup}(ones(length(params_tup)))

    t_compile = @elapsed begin
        v = rate_equation(m, concs, pvals)
        s = rate_equation_string(m)
    end
    @test v isa Float64
    @test s isa String
    @test t_compile < 20.0
end

@testset "Rate equation too large error" begin
    # Manually defined mechanism (11 forms, 16 steps, ~29k terms)
    # triggers the post-hoc check in _raw_symbolic_rate_polys.
    m_manual = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: R1
        steps: begin
            E + A <--> E(A)
            E(A) <--> F(P)
            F(P) <--> F + P
            F + B <--> F(B)
            F(B) <--> E(Q)
            E(Q) <--> E + Q
            E + R1 <--> E(R1)
            E(A) + R1 <--> E(A, R1)
            F(P) + R1 <--> F(P, R1)
            F + R1 <--> F(R1)
            F(B) + R1 <--> F(B, R1)
            E(R1) + A <--> E(A, R1)
            E(A, R1) <--> F(P, R1)
            F(P, R1) <--> F(R1) + P
            F(R1) + B <--> F(B, R1)
            F(B, R1) <--> E(R1) + Q
        end
    end
    @test_throws "polynomial terms" rate_equation_string(m_manual)

    # A directly-constructed random Bi-Bi with an R1 dead-end that binds
    # the free enzyme AND every catalytic form, so each `E(X, R1)` form is
    # reachable two ways (`E(R1)+X` and `E(X)+R1`). Those cycles multiply
    # the King-Altman spanning-tree count past MAX_RATE_EQUATION_TERMS, so
    # the all-SS derivation aborts inside `sym_det`. Built directly rather
    # than enumerated so the guard is exercised deterministically.
    m_cyclic = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: R1
        steps: begin
            E + A <--> E(A)
            E + B <--> E(B)
            E(A) + B <--> E(A, B)
            E(B) + A <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(P) + Q
            E(P, Q) <--> E(Q) + P
            E + P <--> E(P)
            E + Q <--> E(Q)
            E + R1 <--> E(R1)
            E(R1) + A <--> E(A, R1)
            E(A) + R1 <--> E(A, R1)
            E(R1) + B <--> E(B, R1)
            E(B) + R1 <--> E(B, R1)
            E(R1) + P <--> E(P, R1)
            E(P) + R1 <--> E(P, R1)
            E(R1) + Q <--> E(Q, R1)
            E(Q) + R1 <--> E(Q, R1)
        end
    end
    @test_throws "polynomial terms" rate_equation_string(m_cyclic)
end


# ── Single-feature edge cases ─────────────────────────────────────────────
@testset "Allosteric edge cases" begin
    # OnlyR substrate: S binds only in R-state (R-state-active convention).
    # T-state cycle is dead, so all forward catalysis happens through R.
    # As K1 → ∞ (weaker R-state binding), rate vanishes.
    onlyR_sub = @allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)     :: OnlyR
            E(S) <--> E(P)   :: EqualRT
            E(P) ⇌ E + P     :: EqualRT
        end
    end
    concs = (S=1.0, P=0.001)
    base_params = (k2f=10.0, K3=0.5, L=10.0, Keq=1000.0, E_total=1.0)
    rate_strong = rate_equation(onlyR_sub, concs, merge(base_params, (K1=0.01,)))
    rate_weak   = rate_equation(onlyR_sub, concs, merge(base_params, (K1=1e6,)))
    @test rate_strong > 1.0
    @test rate_weak < 1e-3
    @test rate_weak / rate_strong < 1e-5

    # V-type only: all bindings :EqualRT but catalysis :OnlyR. T-state binds
    # substrate normally but cannot catalyze (k_T = 0), so N_T = 0. As L → ∞
    # (T-state dominant) the rate vanishes.
    vtype = @allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)     :: EqualRT
            E(S) <--> E(P)   :: OnlyR
            E(P) ⇌ E + P     :: EqualRT
        end
    end
    vparams = (K1=0.1, k2f=10.0, K3=0.5, Keq=1000.0, E_total=1.0)
    rate_R = rate_equation(vtype, concs, merge(vparams, (L=0.0,)))
    rate_T = rate_equation(vtype, concs, merge(vparams, (L=1e10,)))
    @test rate_R > 1.0
    @test rate_T < 1e-6
    # T-state numerator branch is elided when t_state_dead (any :OnlyR catalytic group);
    # rate is E_total · catN · num_R / (Q_R^catN + L · Q_T^catN). At large L, the T-state
    # enzyme mass dominates the denominator → rate ∝ 1/(1+L).
    @test rate_T * 1e10 < 100.0    # bounded as L grows

    # :OnlyT on a substrate-binding catalytic group → constructor error
    # (R-state convention: relabel so the active state is R, i.e. use :OnlyR).
    @test_throws Exception eval(:(@allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)     :: OnlyT
            E(S) <--> E(P)   :: EqualRT
            E(P) ⇌ E + P     :: EqualRT
        end
    end))

    # :OnlyT on a product-binding catalytic group → constructor error
    @test_throws Exception eval(:(@allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)     :: EqualRT
            E(S) <--> E(P)   :: EqualRT
            E(P) ⇌ E + P     :: OnlyT
        end
    end))

    # :OnlyT on the catalysis SS step → constructor error
    @test_throws Exception eval(:(@allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)     :: EqualRT
            E(S) <--> E(P)   :: OnlyT
            E(P) ⇌ E + P     :: EqualRT
        end
    end))


    # Single-ligand :EqualRT reg site cancels identically → constructor error
    @test_throws Exception eval(:(@allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: I::EqualRT
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)     :: EqualRT
            E(S) <--> E(P)   :: EqualRT
            E(P) ⇌ E + P     :: EqualRT
        end
    end))

    # Stoichiometric infeasibility: Q listed as product but never appears
    # in any reaction step → invariant check rejects.
    rxn_no_q = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C], Q[N]
    end
    s_q1 = EnzymeRates.Step(EnzymeRates.Species(EnzymeRates.Metabolite[], :E),
                            EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                            EnzymeRates.Substrate(:S), true; source_idx = 1)
    s_q2 = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                            EnzymeRates.Species([EnzymeRates.Product(:P)], :E_P),
                            nothing, false; source_idx = 2)
    s_q3 = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Product(:P)], :E_P),
                            EnzymeRates.Species(EnzymeRates.Metabolite[], :E),
                            EnzymeRates.Product(:P), true; source_idx = 3)
    m_no_q = EnzymeRates.Mechanism(rxn_no_q, [[s_q1], [s_q2], [s_q3]])
    @test_throws ErrorException EnzymeRates._assert_mechanism_invariants(m_no_q)

    # Same-kinetics group across different metabolites: group 1 contains both
    # an S-binding and an A-binding step → invariant check rejects.
    rxn_sa = @enzyme_reaction begin
        substrates: S[C], A[N]
        products:   P[CN]
    end
    g_sa1 = EnzymeRates.Step(EnzymeRates.Species(EnzymeRates.Metabolite[], :E),
                             EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                             EnzymeRates.Substrate(:S), true; source_idx = 1)
    g_sa2 = EnzymeRates.Step(EnzymeRates.Species(EnzymeRates.Metabolite[], :E),
                             EnzymeRates.Species([EnzymeRates.Substrate(:A)], :E_A),
                             EnzymeRates.Substrate(:A), true; source_idx = 2)
    g_sa3 = EnzymeRates.Step(
        EnzymeRates.Species([EnzymeRates.Substrate(:S), EnzymeRates.Substrate(:A)], :E_S_A),
        EnzymeRates.Species([EnzymeRates.Product(:P)], :E_P),
        nothing, false; source_idx = 3)
    m_sa = EnzymeRates.Mechanism(rxn_sa, [[g_sa1, g_sa2], [g_sa3]])
    @test_throws ErrorException EnzymeRates._assert_mechanism_invariants(m_sa)

    # Regression: T-state binding K's must be in Kd convention even when
    # `:OnlyR` and `:NonequalRT` catalytic groups coexist. Without the fix,
    # the flat-poly path in _allosteric_num_den_exprs renders T-state K's
    # as `K_T * met` (Ka) instead of `met / K_T` (Kd), silently producing
    # wrong rates whenever a mechanism mixes these two tags. Regression
    # for src/rate_eq_derivation.jl:1395-1396.
    cm_mix = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    m_mix = EnzymeRates.AllostericEnzymeMechanism(
        cm_mix,
        (2, (:NonequalRT, :OnlyR, :NonequalRT)),
        (((:I,), 2, (:OnlyT,)),),
    )
    p_mix = (K1=0.1, k2f=10.0, K3=0.5,
             K1_T=10.0, K3_T=10.0,
             K_I_T_reg1=1.0, L=1.0, Keq=1000.0, E_total=1.0)
    rate_mix = rate_equation(m_mix, (S=10.0, P=0.0, I=0.0), p_mix)
    # With Kd convention (correct): rate ≈ 19.79 (R-state catalysis dominates).
    # With Ka convention (bug): rate ≈ 9.9 — half the correct value.
    @test isapprox(rate_mix, 19.79; rtol=0.05)

    # Sanity: rate_equation_string emits Kd form for T-state K's.
    @test occursin("S / K1_T", rate_equation_string(m_mix))
    @test occursin("P / K3_T", rate_equation_string(m_mix))

    # Regression: :NonequalRT substrate + :EqualRT catalysis must produce
    # zero rate at chemical equilibrium. The framework derives a T-state
    # Haldane (k2r_T) from the :EqualRT k2f because the dep expression
    # for k2r references :NonequalRT K1, so the dep-assignment builder
    # synthesizes a T-name for k2r and substitutes it into N_T.
    cm_mixed = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end
    m_mixed = EnzymeRates.AllostericEnzymeMechanism(
        cm_mixed,
        (2, (:NonequalRT, :EqualRT, :EqualRT)),
        (((:I,), 2, (:NonequalRT,)),),
    )
    Keq_val = 5.0
    p_eq = (K1=0.3, k2f=8.0, K3=0.7,
            K1_T=2.5,
            K_I_reg1=1.0, K_I_T_reg1=4.0,
            L=2.0, Keq=Keq_val, E_total=1.0)
    # At chemical equilibrium: P = Keq · S
    S_eq = 1.5
    P_eq = Keq_val * S_eq
    rate_eq = rate_equation(m_mixed, (S=S_eq, P=P_eq, I=0.5), p_eq)
    @test isapprox(rate_eq, 0.0; atol=1e-10)

    # Empty ligand list at reg site → constructor error
    cm_simple = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
        cm_simple,
        (2, (:NonequalRT, :EqualRT, :EqualRT)),
        (((), 2, ()),),  # empty ligand tuple
    )
end

@testset "rate_equation_string allosteric byte-identical fixture" begin
    m_allo = @allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: R::NonequalRT
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)     :: NonequalRT
            E + S ⇌ E(S)     :: NonequalRT
            E(S) <--> E(P)   :: NonequalRT
        end
    end
    actual = rate_equation_string(m_allo)
    expected = raw"""(; K1, K2, k3f, K1_T, K2_T, k3f_T, K_R_reg1, K_R_T_reg1, L, Keq, E_total) = params
(; S, P, R) = concs
# Haldane constraints:
k3r = (1 / Keq) * K1 * (1 / K2) * k3f
k3r_T = (1 / Keq) * K1_T * (1 / K2_T) * k3f_T
v = E_total * (2 * ((k3f * S / K2 - k3r * P / K1) * (1 + P / K1 + S / K2) * (1 + R / K_R_reg1) ^ 2 + L * (S * k3f_T / K2_T - P * k3r_T / K1_T) * (1 + P / K1_T + S / K2_T) * (1 + R / K_R_T_reg1) ^ 2)) / ((1 + P / K1 + S / K2) ^ 2 * (1 + R / K_R_reg1) ^ 2 + L * (1 + P / K1_T + S / K2_T) ^ 2 * (1 + R / K_R_T_reg1) ^ 2)"""
    @test actual == expected
end

@testset "Parameter-struct allosteric helpers" begin
    cm = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end

    @testset "_onlyR_parameters" begin
        # One :OnlyR catalytic group (RE binding) → single Kd(_, :R).
        aem = EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:OnlyR, :EqualRT, :EqualRT)),
            (((:I,), 1, (:OnlyT,)),),
        )
        am = EnzymeRates.AllostericMechanism(aem)
        params = EnzymeRates._onlyR_parameters(am)
        rendered = Set(EnzymeRates.name(p, am) for p in params)
        # RE binding rep step 1 → single Kd named :K1.
        @test rendered == Set([:K1])

        # SS iso group with :OnlyR → Kfor + Krev pair.
        aem_ss = EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:EqualRT, :OnlyR, :EqualRT)),
            (((:I,), 1, (:OnlyT,)),),
        )
        am_ss = EnzymeRates.AllostericMechanism(aem_ss)
        ps = EnzymeRates._onlyR_parameters(am_ss)
        @test length(ps) == 2
        rendered_ss = Set(EnzymeRates.name(p, am_ss) for p in ps)
        # SS iso rep step 2 → :k2f + :k2r pair.
        @test rendered_ss == Set([:k2f, :k2r])

        # No :OnlyR groups → empty.
        aem_none = EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:EqualRT, :EqualRT, :EqualRT)),
            (((:I,), 1, (:OnlyT,)),),
        )
        @test isempty(EnzymeRates._onlyR_parameters(
            EnzymeRates.AllostericMechanism(aem_none)))
    end

    @testset "_T_rename_parameters" begin
        # Mix of :NonequalRT (rename), :EqualRT (skip), :OnlyR (skip).
        aem = EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalRT, :EqualRT, :NonequalRT)),
            (((:I,), 1, (:OnlyT,)),),
        )
        am = EnzymeRates.AllostericMechanism(aem)
        rename = EnzymeRates._T_rename_parameters(am)

        rendered = Dict{Symbol, Symbol}()
        for (r_p, t_p) in rename
            rendered[EnzymeRates.name(r_p, am)] =
                EnzymeRates.name(t_p, am)
        end
        # Groups 1 (rep=1) and 3 (rep=3) are :NonequalRT RE bindings →
        # Kd(rep, :R) → Kd(rep, :T) for each. Group 2 (rep=2) is :EqualRT
        # and is skipped.
        @test rendered == Dict(:K1 => :K1_T, :K3 => :K3_T)

        # Empty case: no :NonequalRT groups → empty rename.
        aem_none = EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:OnlyR, :EqualRT, :EqualRT)),
            (((:I,), 1, (:OnlyT,)),),
        )
        @test isempty(EnzymeRates._T_rename_parameters(
            EnzymeRates.AllostericMechanism(aem_none)))
    end

    @testset "_all_t_state_parameters" begin
        # :NonequalRT cat group + :NonequalRT reg ligand → both contribute.
        aem = EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalRT, :EqualRT, :NonequalRT)),
            (((:R,), 1, (:NonequalRT,)),),
        )
        am = EnzymeRates.AllostericMechanism(aem)
        params = EnzymeRates._all_t_state_parameters(am)
        rendered = [EnzymeRates.name(p, am) for p in params]

        # Catalytic side: every non-:OnlyR group contributes a T-state
        # parameter (Kd/Kiso for RE, Kfor/Krev for SS).
        @test :K1_T in rendered
        @test :k2f_T in rendered
        @test :k2r_T in rendered
        @test :K3_T in rendered
        # Regulator side: :R is :NonequalRT → K_R_T_reg1 appears.
        @test :K_R_T_reg1 in rendered

        # :OnlyR cat group + :OnlyR reg ligand are skipped.
        aem_skip = EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:OnlyR, :NonequalRT, :EqualRT)),
            (((:R,), 1, (:OnlyR,)),),
        )
        am_skip = EnzymeRates.AllostericMechanism(aem_skip)
        rendered_skip = [EnzymeRates.name(p, am_skip)
                         for p in EnzymeRates._all_t_state_parameters(am_skip)]
        @test :K1_T ∉ rendered_skip          # :OnlyR cat group skipped
        @test :k2f_T in rendered_skip        # :NonequalRT SS iso emits both
        @test :k2r_T in rendered_skip
        @test :K_R_T_reg1 ∉ rendered_skip    # :OnlyR reg ligand skipped
    end
end
