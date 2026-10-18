# ABOUTME: Tests for rate-equation derivation (King-Altman/Cha, mass-action oracle),
# ABOUTME: the parameters API, peak rescaling, rate-equation strings, and the perf gate.

using OrdinaryDiffEqFIRK

include("reference/allosteric_undefvar_reproducers.jl")

const ER = EnzymeRates

# ── Helper functions ────────────────────────────────────────────────────────

# ── Reference QSSA implementation ───────────────────────────────────────────

"""
Independent reference: compute QSSA rate using Laplacian cofactor method.
Works directly with EnzymeMechanism type parameters.
"""
function _testhelper_reference_qssa(
    m::EnzymeMechanism,
    params::NamedTuple,
    concs::NamedTuple,
)
    # Walk the steps of the lifted mechanism, so this helper follows the compiled
    # mechanism representation.
    flat = _testhelper_flat_steps(m)
    enz_names = _testhelper_enzyme_forms(flat)
    n = length(enz_names)
    name_to_idx = Dict(nm => i for (i, nm) in enumerate(enz_names))

    ref_name, nu_ref = _testhelper_reference_metabolite(m)

    # Build rate matrix R[i,j] = pseudo-first-order rate from i to j
    R = zeros(n, n)
    for (step_idx, step) in enumerate(flat)
        i = name_to_idx[EnzymeRates.name(EnzymeRates.from_species(step))]
        j = name_to_idx[EnzymeRates.name(EnzymeRates.to_species(step))]

        m_lhs = EnzymeRates.name.(EnzymeRates.consumed(step))
        m_rhs = EnzymeRates.name.(EnzymeRates.released(step))

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
    for (step_idx, step) in enumerate(flat)
        i = name_to_idx[EnzymeRates.name(EnzymeRates.from_species(step))]
        j = name_to_idx[EnzymeRates.name(EnzymeRates.to_species(step))]

        m_lhs = EnzymeRates.name.(EnzymeRates.consumed(step))
        m_rhs = EnzymeRates.name.(EnzymeRates.released(step))

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
As-written index of canonical regulatory site `site` (at canonical position
`pos`) among `written_sites`, the sites of a raw value in written order: it
matches by ligand-name set so positional Kreg names (`:K_<lig>_reg{site}`) keep
the oracle's written site numbering after canonical site reordering.
"""
function _testhelper_positional_site_idx(site, pos, written_sites)
    want = Set(EnzymeRates.ligands(site))
    j = findfirst(ws -> Set(ws.ligands) == want, written_sites)
    j === nothing ? pos : j
end

"""
The constants of the plain mechanism `m` (values by name in `all_params`, from
`all_parameters`) under the positional names the QSSA and ODE cross-checks read: stored
flat step `i` gives `K{i}` at rapid equilibrium and `k{i}f`, `k{i}r` at steady state, in
its stored direction. Every entry of `all_params` passes through.
"""
function _testhelper_stored_params(m, all_params::NamedTuple)
    mech = _testhelper_lift(m)
    out = Pair{Symbol, Any}[]
    for (i, (s, _)) in enumerate(EnzymeRates._flat_steps(EnzymeRates.steps(mech)))
        if EnzymeRates.is_equilibrium(s)
            push!(out, Symbol("K", i) =>
                all_params[EnzymeRates.name(EnzymeRates.Krapid(s, :None), mech)])
        else
            push!(out,
                Symbol("k", i, "f") =>
                    all_params[EnzymeRates.name(EnzymeRates.Kfor(s, :None), mech)],
                Symbol("k", i, "r") =>
                    all_params[EnzymeRates.name(EnzymeRates.Krev(s, :None), mech)])
        end
    end
    merge((; out...), all_params)
end

"""
The constants of `m` (values by name in `all_params`, from `all_parameters`) under the
positional names a textbook rate law reads, by role. Written step `i` (in the order of
`raw`'s written steps) takes the constants of the stored step it matches, in the
written direction: `k{i}f` and `k{i}r` at steady state, and at rapid equilibrium `K{i}`,
a binding's dissociation constant or, for any other step, the equilibrium constant of the
written direction. The inactive conformation of a `:NonequalAI` group adds `k{i}f_T`,
`k{i}r_T` or `K{i}_T`. The constructor stores every binding with its metabolite taken up,
so a written binding's direction is read from its metabolite's role: a product's binding
is written as its release. Every spec follows that convention; a binding written against
it (`E + P <--> E(P)`) would get its constants swapped. No spec in `MECHANISM_TEST_SPECS`
writes a non-binding step against its stored direction, so the suite never inverts `K`
or swaps `k{i}f` and `k{i}r` for such a step. A ligand of regulatory site `j`
(the written order of `raw`'s sites) gives `K_<ligand>_reg{j}` and, when its inactive
conformation has a constant of its own, `K_<ligand>_T_reg{j}`. `raw` is the mechanism as
written, the DSL walker's value (`_testhelper_MechanismTestSpec.raw`). Every entry of
`all_params` passes through.
"""
function _testhelper_written_params(m, all_params::NamedTuple; raw)
    mech = _testhelper_lift(m)
    allo = mech isa EnzymeRates.AllostericMechanism
    stored = EnzymeRates._flat_steps(EnzymeRates.steps(mech))
    written = EnzymeRates._step_groups(allo ? raw.catalytic : raw)
    out = Pair{Symbol, Any}[]
    for (i, w) in enumerate(Iterators.flatten(written))
        j = findfirst(x -> first(x) == w || first(x) == EnzymeRates._reversed(w), stored)
        j === nothing && error("written step $i matches no stored step")
        s, g = stored[j]
        forward = EnzymeRates.is_binding(w) ?
                  !(EnzymeRates.bound_metabolite(w) isa EnzymeRates.Product) : s == w
        tag = allo ? EnzymeRates.cat_allo_state(mech, g) : :None
        states = tag === :NonequalAI ? ((:A, ""), (:I, "_T")) :
                 ((tag === :OnlyA ? :A : tag, ""),)
        for (state, suffix) in states
            if EnzymeRates.is_equilibrium(s)
                key = EnzymeRates.name(EnzymeRates.Krapid(s, state), mech)
                state === :I && !haskey(all_params, key) && continue
                K = all_params[key]
                push!(out, Symbol("K", i, suffix) =>
                    EnzymeRates.is_binding(s) || forward ? K : 1 / K)
            else
                kf_key = EnzymeRates.name(EnzymeRates.Kfor(s, state), mech)
                state === :I && !haskey(all_params, kf_key) && continue
                kf = all_params[kf_key]
                kr = all_params[EnzymeRates.name(EnzymeRates.Krev(s, state), mech)]
                f, r = forward ? (kf, kr) : (kr, kf)
                push!(out, Symbol("k", i, "f", suffix) => f,
                      Symbol("k", i, "r", suffix) => r)
            end
        end
    end
    if allo
        for (pos, site) in enumerate(EnzymeRates.regulatory_sites(mech))
            j = _testhelper_positional_site_idx(site, pos, raw.sites)
            for (lig, tag) in zip(EnzymeRates.ligands(site), EnzymeRates.allo_states(site))
                for (inactive, suffix) in ((false, ""), (true, "_T"))
                    tag === :EqualAI && inactive && continue
                    K = EnzymeRates._reg_K(mech, lig, tag, inactive)
                    K === nothing ||
                        push!(out, Symbol("K_", lig, suffix, "_reg", j) => all_params[K])
                end
            end
        end
    end
    merge((; out...), all_params)
end

"""Check if mechanism has any rapid-equilibrium steps."""
_testhelper_has_re_steps(m) = any(EnzymeRates.is_equilibrium, _testhelper_flat_steps(m))

"""
Test that `rate_equation` is non-allocating and fast for the given mechanism.
Must be a standalone function to avoid @testset closure boxing.
"""
function _testhelper_test_rate_equation_performance(m, params, concs)
    rate_equation(m, concs, params) # warmup/compile
    allocs = @allocated rate_equation(m, concs, params)
    # Minimum over several batches defeats the GC/scheduling inflation a
    # single mean suffers; accumulating into `acc` (and observing it via the
    # finite-result check) prevents the optimizer from eliding the calls.
    best = Inf
    acc = 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += rate_equation(m, concs, params)
        end
        best = min(best, t / 10_000)
    end
    isfinite(acc) || error("rate_equation produced a non-finite result")
    return allocs, best
end

"""
Generate random independent params + Keq + E_total for testing.
Also returns all_params (the full set of k's + E_total) for reference comparison.
"""
function _testhelper_random_independent_params_concs(
    m, met_names::Vector{Symbol}; rng=Random.default_rng()
)
    # Generate random values for independent params + Keq + E_total
    param_keys = (EnzymeRates.fitted_params(m)..., :Keq, :E_total)
    param_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in param_keys)
    new_params = NamedTuple{param_keys}(param_vals)

    conc_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in met_names)
    concs = NamedTuple{Tuple(met_names)}(conc_vals)

    all_params = all_parameters(m, new_params)
    return new_params, concs, all_params
end

"""
Convert structural-named params to ODE params (large k_if/k_ir for all steps).
Steps sharing a kinetic_group share the same K (or k_f/k_r) param.

For binding RE steps (metabolite on LHS, canonical form):
    K = Kd = kr/kf, so k_if = 1e6, k_ir = 1e6 * K.
For RE isomerization steps (no metabolite, enzyme-only):
    K = Ka = kf/kr, so k_if = 1e6 * K, k_ir = 1e6.
"""
function _testhelper_raw_to_ode_params(m, raw_params)
    mech = m isa EnzymeRates.Mechanism ? m : EnzymeRates.Mechanism(m)
    flat = EnzymeRates._flat_steps(mech)
    # A canonical RE binding step has a metabolite on LHS (canonical form
    # invariant: all RE binding steps are written `E + S ⇌ ES`).
    is_binding_step = Bool[
        EnzymeRates.is_equilibrium(s) && !isempty(EnzymeRates.consumed(s))
        for (s, _) in flat
    ]
    _testhelper_lookup(k) = Float64(raw_params[k])
    param_keys = Symbol[]
    param_vals = Float64[]
    for (i, (step, _)) in enumerate(flat)
        push!(param_keys, Symbol("k$(i)f"))
        push!(param_keys, Symbol("k$(i)r"))
        if EnzymeRates.is_equilibrium(step)
            # Look up the step's structural K key, which its kinetic group shares
            K = _testhelper_lookup(
                EnzymeRates.name(EnzymeRates.Krapid(step, :None), mech))
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
            # SS step: its forward and reverse rate constants
            fwd_key = EnzymeRates.name(EnzymeRates.Kfor(step, :None), mech)
            rev_key = EnzymeRates.name(EnzymeRates.Krev(step, :None), mech)
            push!(param_vals, _testhelper_lookup(fwd_key))
            push!(param_vals, _testhelper_lookup(rev_key))
        end
    end
    push!(param_keys, :E_total)
    push!(param_vals, Float64(raw_params[:E_total]))
    return NamedTuple{Tuple(param_keys)}(Tuple(param_vals))
end

function _testhelper_reference_metabolite(m)
    subs = _testhelper_substrates(m)
    isempty(subs) && error("No substrate found in mechanism")
    name = subs[1]
    coeff = -count(==(name), subs)
    return name, coeff
end

# ── ODE steady-state helpers ────────────────────────────────────────────────

function _testhelper_build_ode_rhs(
    m::EnzymeMechanism,
    params, concs,
)
    # Walk the steps of the lifted mechanism, so this helper follows the compiled
    # mechanism representation.
    flat = _testhelper_flat_steps(m)
    enz_names = _testhelper_enzyme_forms(flat)
    name_to_idx = Dict(nm => i for (i, nm) in enumerate(enz_names))

    step_data = []
    for (step_idx, step) in enumerate(flat)
        i = name_to_idx[EnzymeRates.name(EnzymeRates.from_species(step))]
        j = name_to_idx[EnzymeRates.name(EnzymeRates.to_species(step))]

        m_lhs = EnzymeRates.name.(EnzymeRates.consumed(step))
        m_rhs = EnzymeRates.name.(EnzymeRates.released(step))

        kf = Float64(params[Symbol("k$(step_idx)f")])
        kr = Float64(params[Symbol("k$(step_idx)r")])

        rf = isempty(m_lhs) ? kf : kf * concs[m_lhs[1]]
        rr = isempty(m_rhs) ? kr : kr * concs[m_rhs[1]]

        push!(step_data, (i, j, rf, rr))
    end

    n = length(enz_names)
    function _testhelper_rhs!(du, u, p, t)
        fill!(du, 0.0)
        for (i, j, rf, rr) in step_data
            flux = rf * u[i] - rr * u[j]
            du[i] -= flux
            du[j] += flux
        end
    end
    return _testhelper_rhs!
end

function _testhelper_ode_steady_state_flux(
    m::EnzymeMechanism,
    params, concs,
)
    flat = _testhelper_flat_steps(m)
    E_total = params.E_total
    enz_names = _testhelper_enzyme_forms(flat)
    n = length(enz_names)
    ref_name, nu_ref = _testhelper_reference_metabolite(m)

    u0 = zeros(n)
    u0[1] = E_total

    rhs! = _testhelper_build_ode_rhs(m, params, concs)
    prob = ODEProblem(rhs!, u0, (0.0, 1e6))
    # 1e-13 (vs 1e-12) keeps the stiff multi-step solves converged to the
    # exact rate equation regardless of catalytic-step storage order; at 1e-12
    # a ter-bi corner drifts to ~1e-6, exceeding the cross-check rtol.
    sol = solve(prob, RadauIIA9(); abstol=1e-13, reltol=1e-13)
    u_ss = sol.u[end]

    name_to_idx = Dict(nm => i for (i, nm) in enumerate(enz_names))
    v = 0.0
    for (step_idx, step) in enumerate(flat)
        i = name_to_idx[EnzymeRates.name(EnzymeRates.from_species(step))]
        j = name_to_idx[EnzymeRates.name(EnzymeRates.to_species(step))]

        m_lhs = EnzymeRates.name.(EnzymeRates.consumed(step))
        m_rhs = EnzymeRates.name.(EnzymeRates.released(step))

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

# Runs the rendered code as written: its destructuring lines read `params` and
# `concs`, its constraint lines define the dependent parameters, and its `v` line
# is the value.
function _testhelper_eval_rate_string(s, params, concs)
    code = "let params = $params, concs = $concs\n" *
           s * "\nend"
    eval(Meta.parse(code))
end

# ── Modular test functions for _testhelper_MechanismTestSpec ────────────────

function _testhelper_test_structure(spec::_testhelper_MechanismTestSpec)
    m = spec.mechanism
    @testset "Structure" begin
        flat = _testhelper_flat_steps(m)
        @test length(_testhelper_enzyme_forms(flat)) == spec.expected_n_states
        @test length(flat) == spec.expected_n_steps
        @test metabolites(m) == Tuple(spec.metabolite_names)
        # Structural parameter names must be injective on the load-bearing
        # paths the fitter consumes (parameters + fitted_params). A collision
        # would silently shorten the destructured params NamedTuple.
        @test allunique(EnzymeRates.parameters(m))
        @test allunique(EnzymeRates.fitted_params(m))
    end
end

"""Classify a dep expression as Haldane (RHS references Keq), Mirror
(RHS is a single Symbol), or Wegscheider (RHS Expr without Keq)."""
function _testhelper_classify_dep_expr(expr)
    if expr isa Symbol
        return :mirror
    elseif EnzymeRates._mentions(expr, :Keq)
        return :haldane
    else
        return :wegscheider
    end
end

# Parameter-name symbols referenced on the RHS of a dependent-param expression
# `x`, skipping the operator/function-name slot of a `:call` Expr (so `:+`,
# `:*`, `:/`, `:^`, `:sqrt`, ... never appear as candidate parameter names).
_testhelper_rhs_syms(x) = x isa Symbol ? [x] :
    x isa Expr ? reduce(
        vcat,
        map(_testhelper_rhs_syms, x.head === :call ? x.args[2:end] : x.args);
        init=Symbol[]) :
    Symbol[]

"""
The dependent-parameter graph `dep` (with independent parameters `indep`) is sound
iff (a) no symbol is both independent (fitted) and dependent, (b) the dep→dep
dependency graph is acyclic, and (c) every dependent expression's RHS symbol
is a fitted param, another dependent param, or `Keq`.
"""
function _testhelper_dep_graph_is_sound(dep, indep)
    indepset = Set(indep)
    isempty(intersect(indepset, keys(dep))) || return false

    depset = Set(keys(dep))
    edges = Dict(k => Symbol[s for s in _testhelper_rhs_syms(v) if s in depset]
                 for (k, v) in dep)
    state = Dict{Symbol, Int}()  # 0 = unseen, 1 = on-stack, 2 = done
    ok = true
    function _testhelper_dfs(n)
        get(state, n, 0) == 2 && return
        if get(state, n, 0) == 1
            ok = false
            return
        end
        state[n] = 1
        for m in get(edges, n, Symbol[])
            _testhelper_dfs(m)
            ok || return
        end
        state[n] = 2
    end
    for k in keys(dep)
        _testhelper_dfs(k)
        ok || break
    end
    ok || return false

    known = union(indepset, depset, Set([:Keq]))
    for (_, v) in dep, s in _testhelper_rhs_syms(v)
        s in known || return false
    end
    true
end

function _testhelper_test_constraint_counting(spec::_testhelper_MechanismTestSpec)
    m = spec.mechanism
    @testset "Constraints" begin
        dep_exprs, indep = EnzymeRates._dependent_param_exprs(typeof(m))
        # A dependent (Haldane/Wegscheider-derived) parameter must never appear in
        # the independent set returned as fitted_params. The LDH i-state specs have a
        # shared :EqualAI catalytic reverse rate that is Haldane-dependent in the
        # A-state and that the dead I-state's law reads too: one class in both
        # conformations, so one dependent.
        @test _testhelper_dep_graph_is_sound(dep_exprs, indep)
        n_haldane = 0
        n_mirror = 0
        n_wegscheider = 0
        for (_, expr) in dep_exprs
            cat = _testhelper_classify_dep_expr(expr)
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
        @test length(indep) == _testhelper_declared_dimension(m)
    end
end

function _testhelper_test_reference_qssa(
    spec::_testhelper_MechanismTestSpec; n_trials=20, seed=42
)
    m = spec.mechanism
    # Reference QSSA only works for all-SS mechanisms
    _testhelper_has_re_steps(m) && return
    met_names = spec.metabolite_names
    @testset "Reference QSSA" begin
        rng = Random.MersenneTwister(seed)
        @test all(1:n_trials) do _
            new_params, concs, all_params =
                _testhelper_random_independent_params_concs(
                    m, met_names; rng=rng)
            # _testhelper_reference_qssa reads stored-step k$(i)f/k$(i)r keys
            pos_params = _testhelper_stored_params(m, all_params)
            isapprox(
                rate_equation(m, concs, new_params),
                _testhelper_reference_qssa(m, pos_params, concs);
                rtol=spec.reference_rtol)
        end
    end
end

function _testhelper_test_analytical_rate(
    spec::_testhelper_MechanismTestSpec; n_trials=20, seed=1001
)
    # Skip if no analytical rate function provided
    spec.analytical_rate_fn === nothing && return

    m = spec.mechanism
    met_names = spec.metabolite_names
    @testset "Analytical Rate" begin
        rng = Random.MersenneTwister(seed)
        @test all(1:n_trials) do _
            new_params, concs, all_params =
                _testhelper_random_independent_params_concs(
                    m, met_names; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(_testhelper_written_params(m, all_params; raw=spec.raw),
                      (Et=Et,))
            p_pkg = merge(new_params, (E_total=Et,))
            isapprox(
                rate_equation(m, concs, p_pkg),
                spec.analytical_rate_fn(p, concs);
                rtol=1e-10)
        end
    end
end

function _testhelper_test_haldane_equilibrium(spec::_testhelper_MechanismTestSpec; seed=42)
    m = spec.mechanism
    met_names = spec.metabolite_names
    @testset "Haldane Equilibrium" begin
        rng = Random.MersenneTwister(seed)
        new_params, _, _ = _testhelper_random_independent_params_concs(
            m, met_names; rng=rng)
        Keq = new_params.Keq
        n_prods = length(_testhelper_products(m))
        # Build equilibrium concentrations: prod(P_i) / prod(S_i) = Keq
        # Set all substrates to 1.0, distribute Keq^(1/n_prods) across products
        sub_names = _testhelper_substrates(m)
        prod_names = _testhelper_products(m)
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

function _testhelper_test_performance(spec::_testhelper_MechanismTestSpec; seed=42)
    m = spec.mechanism
    met_names = spec.metabolite_names
    @testset "Performance" begin
        rng = Random.MersenneTwister(seed)
        params, concs, _ = _testhelper_random_independent_params_concs(
            m, met_names; rng=rng)
        allocs, t = _testhelper_test_rate_equation_performance(m, params, concs)
        @test allocs == 0
        # Locally 120 ns leaves about 4× headroom over the slowest fixture laws, PFK-1 and
        # "MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites" at about 30 ns. Shared CI
        # runners occasionally run those two 4-7× slower and miss it; rerun such a job
        # rather than loosening the bound. The 0-alloc bound above is the strict part.
        @test t < 120e-9
    end
end

function _testhelper_test_ode_steadystate(
    spec::_testhelper_MechanismTestSpec; n_trials=10, seed=42
)
    m = spec.mechanism
    met_names = spec.metabolite_names
    @testset "ODE Steady-State" begin
        rng = Random.MersenneTwister(seed)
        has_re = _testhelper_has_re_steps(m)
        @test all(1:n_trials) do _
            new_params, concs, all_params =
                _testhelper_random_independent_params_concs(
                    m, met_names; rng=rng)
            # Convert structural params to positional k_if/k_ir for ODE
            pos_params = _testhelper_stored_params(m, all_params)
            ode_params = has_re ?
                _testhelper_raw_to_ode_params(m, all_params) :
                pos_params
            v_ode = _testhelper_ode_steady_state_flux(m, ode_params, concs)
            v_ka = rate_equation(m, concs, new_params)
            # Use looser tolerance for RE mechanisms (large rate approximation)
            rtol = has_re ? 1e-3 : spec.ode_rtol
            isapprox(v_ode, v_ka; rtol=rtol)
        end
    end
end

function _testhelper_test_rate_equation_string(spec::_testhelper_MechanismTestSpec)
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
            new_params, concs, _ =
                _testhelper_random_independent_params_concs(
                    m, met_names; rng=rng)
            isapprox(
                rate_equation(m, concs, new_params),
                _testhelper_eval_rate_string(s, new_params, concs);
                rtol=1e-10)
        end

        has_num = spec.expected_factored_num !== nothing
        has_denom = spec.expected_factored_denom !== nothing
        if has_num || has_denom
            num_str, denom_str = _testhelper_extract_num_denom(last(split(s, "\n")))
            @test num_str !== nothing
            @test denom_str !== nothing
            has_num && @test num_str == spec.expected_factored_num
            has_denom && @test denom_str == spec.expected_factored_denom
        end
    end
end

"""
Extract numerator and denominator strings from rate equation v-line.
Handles nested parentheses via depth counting.
"""
function _testhelper_extract_num_denom(v_line::AbstractString)
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

"""
Test that the peak turnover of `spec`'s mechanism at `params` (fitted constants, `Keq`
and `E_total`; seed 42's draw unless given) equals `spec.analytical_kcat_fn`'s value.
"""
function _testhelper_test_analytical_kcat(spec::_testhelper_MechanismTestSpec;
                                          params=nothing)
    spec.analytical_kcat_fn === nothing && return
    m = spec.mechanism
    if params === nothing
        params = _testhelper_random_independent_params(m; rng=Random.MersenneTwister(42))
    end
    @testset "Peak turnover approaches the analytical kcat" begin
        # Every fitted constant and Keq lies in 0.1–10, so over a range reaching 1e12 the
        # peak sits at saturating substrates and equals the spec's analytical value to
        # within K/1e12: the kcat at zero products or, where a product raises the
        # turnover, the supremum over product levels.
        peak = EnzymeRates._peak_turnover(m, params, (1e-12, 1e12)).value
        # The textbook formula may read a dependent rate constant, so the bridge
        # reads every constant through all_parameters.
        p = merge(_testhelper_written_params(m, all_parameters(m, params); raw=spec.raw),
                  (Et=params.E_total,))
        @test peak ≈ spec.analytical_kcat_fn(p) rtol=1e-9
    end
end

function _testhelper_test_peak_rescaling(spec::_testhelper_MechanismTestSpec; seed=100)
    m = spec.mechanism
    @testset "Peak rescaling" begin
        rng = Random.MersenneTwister(seed)
        params = _testhelper_random_independent_params(m; rng)
        range = (1e-12, 1e3)
        peak = EnzymeRates._peak_turnover(m, params, range).value
        @test peak > 0

        # Scaling every steady-state rate constant by α scales the peak by α.
        ss_names = EnzymeRates._ss_rate_constant_names(m)
        α = 0.1 + 9.9 * rand(rng)
        scaled_params = NamedTuple{keys(params)}(Tuple(
            k in ss_names ? v * α : v
            for (k, v) in zip(keys(params), values(params))
        ))
        @test EnzymeRates._peak_turnover(m, scaled_params, range).value ≈
            α * peak rtol=1e-9

        # The peak scale makes the peak 1 and leaves every other constant alone.
        norm, spread = rescale_parameter_values(m, params)
        @test spread === nothing
        @test EnzymeRates._peak_turnover(m, norm, range).value ≈ 1.0 rtol=1e-9
        for k in keys(params)
            if !(k in ss_names)
                @test norm[k] == params[k]
            end
        end

        # Every rate changes by the same factor.
        met_names = metabolites(m)
        concs = NamedTuple{Tuple(met_names)}(Tuple(0.5 + rand(rng) for _ in met_names))
        v_orig = rate_equation(m, concs, params)
        @test rate_equation(m, concs, norm) / v_orig ≈ 1.0 / peak rtol=1e-8

        # The peak bounds the rescaled turnover at every corner of the range.
        norm_turnover = merge(norm, (E_total=1.0,))
        corners = Iterators.product(ntuple(_ -> range, length(met_names))...)
        @test all(c -> rate_equation(m, NamedTuple{met_names}(c), norm_turnover) <=
                       1 + 1e-9, corners)

        # An anchor at those concentrations sets the turnover there.
        target = sign(v_orig) * (0.1 + 9.9 * rand(rng))
        anchored, spread_one = rescale_parameter_values(m, params, [concs => target])
        @test spread_one == 0.0
        # The anchors' turnover is per enzyme, so the scale leaves E_total as passed in.
        @test anchored.E_total == params.E_total
        @test rate_equation(m, concs, merge(anchored, (E_total=1.0,))) ≈ target rtol=1e-9
    end
end

"""
Assert `rate_equation` stays finite when any single metabolite concentration
is zero (real kinetic data routinely has zeros). Zeroing a substrate or product
must still leave a nonzero net rate — the opposite-direction metabolites keep
driving flux — so a `0.0` there normally signals a `1/conc` term that blew the
denominator to `Inf`. A regulator zero (essential activator) can legitimately
zero the rate, so only finiteness is required there. Zeroing a metabolite that
sits in the free-enzyme weight `D[g_free]` also keeps a finite rate with a
nonzero reverse flux; the mass-action gates in `allosteric_ground_truth.jl`
check that case against the ground truth.
"""
function _testhelper_test_zero_metabolite_finite(spec::_testhelper_MechanismTestSpec)
    m = spec.mechanism
    @testset "Zero-metabolite finiteness" begin
        rng = Random.MersenneTwister(777 + hash(spec.name) % 1000)
        mets = collect(metabolites(m))
        sub_prod = Set{Symbol}(_testhelper_substrates(m))
        union!(sub_prod, _testhelper_products(m))
        params = _testhelper_random_independent_params(m; rng)
        for zeroed in mets
            cvals = Tuple(n == zeroed ? 0.0 : 0.5 + rand(rng) for n in mets)
            concs = NamedTuple{Tuple(mets)}(cvals)
            v = rate_equation(m, concs, params)
            @test isfinite(v)
            zeroed in sub_prod && @test v != 0.0
        end
    end
end

"""The spec's rate law against the exact steady state of its mass-action network
(`_testhelper_check_law_against_oracle`); an allosteric spec runs with one catalytic
subunit, where the MWC law is the exact steady state."""
function _testhelper_test_mass_action_oracle(spec::_testhelper_MechanismTestSpec)
    @testset "Mass-action oracle" begin
        _testhelper_check_law_against_oracle(spec.mechanism)
    end
end

"""
Run all tests for a mechanism specification.
Organizes tests by mechanism: all tests for one mechanism together.
"""
function _testhelper_run_all_tests(spec::_testhelper_MechanismTestSpec)
    @testset "$(spec.name)" begin
        _testhelper_test_structure(spec)
        _testhelper_test_constraint_counting(spec)
        _testhelper_test_mass_action_oracle(spec)
        _testhelper_test_reference_qssa(spec)
        # Only runs if analytical_rate_fn provided
        _testhelper_test_analytical_rate(spec)
        _testhelper_test_haldane_equilibrium(spec)
        _testhelper_test_performance(spec)
        _testhelper_test_rate_equation_string(spec)
        _testhelper_test_zero_metabolite_finite(spec)
        spec.run_ode_test && _testhelper_test_ode_steadystate(spec)
        # Only runs if analytical_kcat_fn provided
        _testhelper_test_analytical_kcat(spec)
        _testhelper_test_peak_rescaling(spec)
    end
end

# ── Main test loop ──────────────────────────────────────────────────────────

@testset "Enzyme Derivation Tests" begin
    # A truncated table would let every per-spec gate pass vacuously.
    @test length(MECHANISM_TEST_SPECS) == 42
    for spec in MECHANISM_TEST_SPECS
        _testhelper_run_all_tests(spec)
    end
end

# Seed 42 puts m_OnlyA_prod's peak turnover at zero products on some Julia versions and
# above it on others (MersenneTwister's stream differs), so fixed constants test the
# peak at a nonzero product on every version.
@testset "m_OnlyA_prod's peak turnover at a nonzero product" begin
    spec = only(s for s in MECHANISM_TEST_SPECS if s.name == "m_OnlyA_prod")
    # K1 = K3 = k2f = 1 and Keq = 10, so k2r = k2f·K3/(Keq·K1) = 0.1 and
    # k2f·(L − 1) = 8 exceeds k2r·(1 + L) = 1: the product raises the turnover.
    params = (K_A_EP_to_E_P=1.0, K_ES_to_E_S=1.0, k_A_ES_to_EP=1.0, L=9.0, Keq=10.0,
              E_total=1.0)
    # Above the zero-product kcat k2f/(1 + L) = 0.1 (the peak is 0.140).
    @test EnzymeRates._peak_turnover(spec.mechanism, params, (1e-12, 1e12)).value > 0.13
    _testhelper_test_analytical_kcat(spec; params)
end

# A fused and a plain binding of B share one rapid-equilibrium kinetic group.
const _testhelper_fused_and_plain_binding =
    @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
        E + A <--> E(A); (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))
        E(Q) + P <--> E(P, Q); E + Q ⇌ E(Q) end end)

# Mechanisms with fused steps, shared by the mass-action and orientation testsets.
const _testhelper_fused_cases = [
        # uni-uni: fused release / fused binding, both orientations, both flags
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> E(S); E(S) <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S ⇌ E(S); E(S) <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> E(S); E(S) ⇌ E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> E(S); E + P <--> E(S) end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> E(P); E(P) <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E(P) <--> E + S; E + P <--> E(P) end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E(P) ⇌ E + S; E(P) <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E(P) <--> E + S; E(P) ⇌ E + P end end),
        # ordered bi-bi, one central complex (Segel IX-87) and mixes
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> E(A); E(A) + B <--> E(A, B); E(A, B) <--> E(Q) + P
            E(Q) <--> E + Q end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> E(A); E(A) + B <--> E(A, B); E(A, B) ⇌ E(Q) + P
            E(Q) <--> E + Q end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A ⇌ E(A); E(P, Q) <--> E(A) + B; E(P, Q) <--> E(Q) + P
            E(Q) ⇌ E + Q end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A ⇌ E(A); E(P, Q) <--> E(A) + B; E(P, Q) ⇌ E(Q) + P
            E(Q) ⇌ E + Q end end),
        # ping-pong with an RE fused release
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> E(A); E(A) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B <--> E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q end end),
        # merged random bi-bi with RE bindings and releases, and a reversed fused
        # binding inside an all-SS merged random bi-bi
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A ⇌ E(A); E + B ⇌ E(B); E(A) + B ⇌ E(A, B); E(B) + A ⇌ E(A, B)
            E(A, B) <--> E(Q) + P; E(A, B) <--> E(P) + Q; E(Q) ⇌ E + Q; E(P) ⇌ E + P
        end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> E(A); E + B <--> E(B); E(A) + B <--> E(P, Q)
            E(P, Q) <--> E(B) + A; E(P, Q) <--> E(Q) + P; E(P, Q) <--> E(P) + Q
            E(Q) <--> E + Q; E(P) <--> E + P end end),
        # dead ends on a merged complex
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; modifiers: I
            steps: begin
            E + A ⇌ E(A); E(A) + B ⇌ E(A, B); E(A, B) <--> E(Q) + P; E(Q) ⇌ E + Q
            E + I ⇌ E(I); E(Q) + I ⇌ E(I, Q) end end),
        # ping-pong: parallel product-release route via a mixed substrate/product complex
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; modifiers: I
            steps: begin
            E + A ⇌ E(A); E(I) + A ⇌ E(A, I); E(Q) + A ⇌ E(A, Q)
            E + I ⇌ E(I); E(A) + I ⇌ E(A, I)
            E(P; residual = A - P) + I ⇌ E(I, P; residual = A - P)
            E(; residual = A - P) + I ⇌ E(I; residual = A - P)
            E + Q <--> E(Q); E(A) + Q <--> E(A, Q)
            E(A) <--> E(P; residual = A - P)
            E(A, I) <--> E(I, P; residual = A - P)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            E(I; residual = A - P) + P ⇌ E(I, P; residual = A - P)
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E(P; residual = A - P) + B ⇌ E(B, P; residual = A - P)
            E(; residual = A - P) + B ⇌ E(B; residual = A - P) end end),
        # both substrates bind in one rapid-equilibrium step
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A + B ⇌ E(A, B); E(A, B) <--> E(P, Q); E(P, Q) <--> E(Q) + P
            E(Q) <--> E + Q end end),
        # the first substrate also binds as a dead-end modifier copy, at rapid
        # equilibrium on E and at steady state on E(Q); both the derivation and
        # the oracle count the first substrate by name
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; modifiers: A
            steps: begin
            E + A ⇌ E(A); E(A) + B <--> E(A, B); E(A, B) <--> E(Q) + P; E(Q) <--> E + Q
            E + A::Mod ⇌ E(A::Mod); E(Q) + A::Mod <--> E(A::Mod, Q) end end),
        _testhelper_fused_and_plain_binding,
]

@testset "the free-enzyme weight is free E's weight in the reduced denominator" begin
    # Ping-pong with a rapid-equilibrium second chemistry step: every form lies in one
    # rapid-equilibrium segment rooted at free E, where E(; residual = A - P) sits at
    # Q·K/B relative to E. Reducing the polynomials to lowest terms in the
    # concentrations multiplies them by B, so free E weighs B in the reduced
    # denominator; the MWC combination divides each conformation by this weight.
    em = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E + Q ⇌ E(Q)
        end
    end
    _, _, d_free = EnzymeRates._raw_symbolic_rate_polys(typeof(em))
    @test d_free == EnzymeRates.poly_sym(:B)
end

@testset "a fused release keeps its catalytic rate constant fitted" begin
    # E(S) <--> E + P gives off P while turning S into P. It is stored as the binding it
    # reverses, E + P → E(S), whose forward constant k_E_P_to_ES runs against the
    # reaction. The Haldane relation makes that constant dependent, as k₋₂ of the
    # textbook mechanism, and the catalytic constant k_ES_to_E_P stays fitted.
    ss = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end
    @test Set(EnzymeRates.fitted_params(ss)) ==
          Set([:k_E_S_to_ES, :k_ES_to_E_S, :k_ES_to_E_P])
    re = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
        end
    end
    @test Set(EnzymeRates.fitted_params(re)) == Set([:K_ES_to_E_S, :k_ES_to_E_P])
end

@testset "fused steps derive the mass-action rate" begin
    for em in _testhelper_fused_cases
        _testhelper_check_law_against_oracle(em)
    end
    # The B group holds the fused E(A) + B → E(P, Q) and the dead-end E(Q) + B ⇌ E(B, Q).
    # Fitted: A (2) + B (1) + P (2) + Q (1) − 1 Haldane = 5.
    @test length(ER.fitted_params(_testhelper_fused_and_plain_binding)) == 5
end

# Theorell–Chance and two-metabolite mechanisms, shared by the mass-action and
# orientation testsets.
const _testhelper_tc_cases = (
    tc_ss = @enzyme_mechanism(begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end),
    tc_re_outer = @enzyme_mechanism(begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(Q) + P
            E(Q) ⇌ E + Q
        end
    end),
    tc_re_step = @enzyme_mechanism(begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(Q) + P
            E(Q) <--> E + Q
        end
    end),
    two_in = @enzyme_mechanism(begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A + B <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E + P + Q
        end
    end),
)

@testset "Theorell–Chance and two-metabolite steps" begin
    (; tc_ss, tc_re_outer, tc_re_step, two_in) = _testhelper_tc_cases
    @test length(ER.fitted_params(tc_ss)) == 5
    @test length(ER.fitted_params(tc_re_outer)) == 3
    for em in (tc_ss, tc_re_outer, tc_re_step, two_in)
        _testhelper_check_law_against_oracle(em)
    end
    # canonical orientation: the substrate side is `from`
    s = only(s for g in ER.steps(ER.Mechanism(tc_ss)) for s in g
             if !ER.is_iso(s) && !ER.is_binding(s))
    @test ER.name(ER.from_species(s)) == :EA &&
          ER.consumed(s) == ER.Metabolite[ER.Substrate(:B)]
    @test :k_EA_B_to_EQ_P in _testhelper_step_constant_names(tc_ss)
    @test :k_EQ_P_to_EA_B in _testhelper_step_constant_names(tc_ss)
end

@testset "reversing written steps changes nothing" begin
    _testhelper_rev(s) = ER.Step(ER.to_species(s), ER.from_species(s), ER.released(s),
                                 ER.consumed(s), ER.is_equilibrium(s))
    rng = Random.MersenneTwister(7)
    # Every step reversed, a random subset, and every group's steps in reverse order.
    _testhelper_variants(groups) = (
        [[_testhelper_rev(s) for s in g] for g in groups],
        [[rand(rng, Bool) ? _testhelper_rev(s) : s for s in g] for g in groups],
        [reverse(g) for g in groups])
    # The Theorell–Chance step of `tc_ss` written backwards.
    tc_backwards = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(Q) + P <--> E(A) + B
            E(Q) <--> E + Q
        end
    end
    @test ER.Mechanism(tc_backwards) == ER.Mechanism(_testhelper_tc_cases.tc_ss)
    @test ER.rate_equation_string(tc_backwards) ==
          ER.rate_equation_string(_testhelper_tc_cases.tc_ss)
    ems = Any[spec.mechanism for spec in MECHANISM_TEST_SPECS
              if spec.mechanism isa ER.EnzymeMechanism]
    append!(ems, _testhelper_fused_cases, values(_testhelper_tc_cases), [tc_backwards])
    for em in ems
        m = ER.Mechanism(em)
        for flipped in _testhelper_variants(ER.steps(m))
            m2 = ER.Mechanism(ER.reaction(m), flipped)
            @test m2 == m
            @test ER.rate_equation_string(ER.compile_mechanism(m2)) ==
                  ER.rate_equation_string(em)
        end
    end
    for spec in MECHANISM_TEST_SPECS
        spec.mechanism isa ER.AllostericEnzymeMechanism || continue
        am = ER.AllostericMechanism(spec.mechanism)
        for flipped in _testhelper_variants(ER.steps(am))
            am2 = ER.AllostericMechanism(ER.reaction(am), flipped, ER.cat_allo_states(am),
                                         ER.catalytic_multiplicity(am),
                                         ER.regulatory_sites(am))
            @test am2 == am
            @test ER.rate_equation_string(ER.compile_mechanism(am2)) ==
                  ER.rate_equation_string(spec.mechanism)
        end
    end
end

# Canonicalizing catalytic-step order in the Mechanism constructor must keep
# each catalytic kinetic group's allosteric-state tag bound to the SAME
# physical step. HK declares BOTH ATP binding and the catalytic conversion as
# :OnlyA; after construction the two :OnlyA tags must still sit on those two
# physical steps (not on some other group canonical reordering moved into their
# source slots): exactly one :OnlyA group is the ATP-binding step, and the other
# is the metabolite-free catalytic conversion.
@testset "Allosteric cat_allo_state stays bound to its step" begin
    hk = only(s for s in MECHANISM_TEST_SPECS if s.name == "HK")
    am = EnzymeRates.AllostericMechanism(hk.mechanism)
    onlyA_groups = [g for g in eachindex(EnzymeRates.steps(am))
                    if EnzymeRates.cat_allo_state(am, g) === :OnlyA]
    @test length(onlyA_groups) == 2
    onlyA_bms = [EnzymeRates.bound_metabolite(first(EnzymeRates.steps(am)[g]))
                 for g in onlyA_groups]
    @test count(bm -> bm !== nothing && EnzymeRates.name(bm) === :ATP,
                onlyA_bms) == 1
    @test count(isnothing, onlyA_bms) == 1
end

# ── §5a regression: every inactive-state parameter must be defined ──────────
# Allosteric mechanisms from an LDH `identify_rate_equation` run crashed with
# `UndefVarError` on undefined I-state parameters: the I-state polynomials
# referenced names the dep-assignment and destructuring machinery never emitted.
# The undefined constants were an I-state lactate dissociation constant and an
# I-state reverse chemistry rate in dead inactive states, and an I-state NAD
# binding rate in a live one. Those mechanisms are not kept here. The regression
# runs the three `LDH i-state …` specs of MECHANISM_TEST_SPECS, which exposed the
# Bug-2 fitted_params leak (defined in
# mechanism_definitions_for_test_enzyme_derivation.jl), and checks that each rate
# equation defines every name it references.
const _LDH_ISTATE_MECHS = [spec.mechanism for spec in MECHANISM_TEST_SPECS
                           if startswith(spec.name, "LDH i-state")]

# Parameter Symbols referenced on an assignment/`v` RHS but never defined
# (destructured from `params`/`concs` or assigned as an LHS). Empty ⟺ the
# rendered rate equation is closed: every referenced name has a definition.
function _testhelper_undefined_rhs_symbols(s::AbstractString)
    ident = r"[A-Za-z_][A-Za-z0-9_]*"
    defined = Set{Symbol}()
    referenced = Set{Symbol}()
    for line in split(s, '\n')
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, "#")) && continue
        if occursin("= params", line) || occursin("= concs", line)
            for mt in eachmatch(ident, split(line, '=')[1])
                push!(defined, Symbol(mt.match))
            end
        else
            lhs, rhs = split(line, '='; limit=2)
            push!(defined, Symbol(strip(lhs)))
            for mt in eachmatch(ident, rhs)
                push!(referenced, Symbol(mt.match))
            end
        end
    end
    setdiff(referenced, defined)
end

@testset "§5a I-state parameters are all defined (regression)" begin
    # A renamed spec would silently empty the name filter.
    @test length(_LDH_ISTATE_MECHS) == 3
    for em in _LDH_ISTATE_MECHS
        pnames = EnzymeRates.fitted_params(em)
        mets = EnzymeRates.metabolites(em)
        prods = Set(_testhelper_products(em))
        params = merge(NamedTuple{pnames}(ntuple(_ -> 1.3, length(pnames))),
                       (Keq = 20000.0, E_total = 1.0))
        # No UndefVarError: the @generated body compiles and evaluates finite at
        # products = 0 (initial-rate conditions).
        concs0 = NamedTuple{mets}(
            ntuple(i -> mets[i] in prods ? 0.0 : 1.5, length(mets)))
        @test isfinite(EnzymeRates.rate_equation(em, concs0, params))
        # DEFINED ⊇ REFERENCED on the rendered transcript.
        @test isempty(
            _testhelper_undefined_rhs_symbols(EnzymeRates.rate_equation_string(em)))
    end
end

@testset "allosteric dependent-param graph is sound" begin
    for T in ALLOSTERIC_UNDEFVAR_REPRODUCERS
        @test _testhelper_dep_graph_is_sound(EnzymeRates._dependent_param_exprs(T)...)
    end
end

@testset "allosteric reproducers: detailed balance" begin
    # The dep/indep graph is not always sufficient: a dep expression can be
    # structurally sound (acyclic, every RHS symbol defined) while the
    # code-generated assignment order still leaves a forward reference
    # unresolved, producing a runtime UndefVarError. This is the
    # user-visible symptom for every reproducer, independent of which
    # structural check above flags it. Being callable is necessary but not
    # sufficient: the equation must also satisfy `v = 0` at `Q = Keq` for
    # arbitrary parameter values. The single combined constraint solve
    # in `_dependent_param_exprs` ties every cross-state affinity split
    # directly, so all three reproducers (D1's `:EqualAI`-shared `koff` merge, D2's
    # `:NonequalAI` split, D3's steady-state speed) are detailed-balance-correct.
    for T in ALLOSTERIC_UNDEFVAR_REPRODUCERS
        m = T()
        pn = collect(EnzymeRates.fitted_params(m))
        subs = _testhelper_substrates(m)
        prods = _testhelper_products(m)
        mets = collect(EnzymeRates.metabolites(m))
        Keq = 20000.0
        p_each = Keq^(1 / length(prods))
        rng = Random.MersenneTwister(1)
        maxv = 0.0
        for _ in 1:100
            pv = Dict(p => exp(2 * randn(rng)) for p in pn)
            params = NamedTuple{(pn..., :Keq, :E_total)}(
                (Tuple(pv[p] for p in pn)..., Keq, 1.0))
            cv = Dict{Symbol,Float64}()
            for s in subs; cv[s] = 1.0; end
            for p in prods; cv[p] = p_each; end
            for s in mets; haskey(cv, s) || (cv[s] = 2.5); end
            concs = NamedTuple{Tuple(mets)}(Tuple(cv[s] for s in mets))
            maxv = max(maxv, abs(rate_equation(m, concs, params)))
        end
        @test maxv < 1e-8
    end
end

# ── Standalone tests ───────────────────────────────────────────────────────────

"""The `return E_total * (num) / (den)` expression of the generated `rate_equation` body of
`m`."""
function _testhelper_generated_rate_return(m)
    dep, indep = ER._dependent_param_exprs(typeof(m))
    last(ER._rate_body(typeof(m), (indep..., :Keq, :E_total), dep).args)
end

@testset "rate_equation polynomial body uses 2-arg +/* calls" begin
    # The fitter calls rate_equation millions of times per CV fold. The
    # polynomial body of the generated code MUST
    # have exactly 2 operands per +/* call so LLVM inlines the binary
    # Float64 path; n-ary varargs above ~30 terms boxes the argument
    # tuple and turns 100ns/0B into 1µs/2KB per call.
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "Random-order Bi-Bi")
    bad = Expr[]
    function _testhelper_walk!(e)
        if e isa Expr
            if e.head == :call && !isempty(e.args) &&
               e.args[1] isa Symbol && e.args[1] in (:+, :*) &&
               length(e.args) != 3
                push!(bad, e)
            end
            start = e.head == :call ? 2 : 1
            for i in start:length(e.args)
                _testhelper_walk!(e.args[i])
            end
        end
    end
    _testhelper_walk!(_testhelper_generated_rate_return(spec.mechanism))
    @test isempty(bad)
end

# The three largest rate laws among the children of one `expand_mechanisms` round over
# 60 bi-bi seeds: random-order merged complexes with dead-end bindings, every step steady
# state, 396 denominator terms each.
const _testhelper_largest_enumerated_mechanisms = [
    @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
        (E + A <--> E(A), E(B) + A <--> E(P, Q))
        (E + B <--> E(B), E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q))
        (E + P <--> E(P), E(Q) + P <--> E(P, Q))
        (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q)) end end),
    @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
        (E + A <--> E(A), E(B) + A <--> E(P, Q), E(Q) + A <--> E(A, Q))
        (E + B <--> E(B), E(A) + B <--> E(P, Q))
        (E + P <--> E(P), E(Q) + P <--> E(P, Q))
        (E + Q <--> E(Q), E(A) + Q <--> E(A, Q), E(P) + Q <--> E(P, Q)) end end),
    @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
        (E + A <--> E(A), E(B) + A <--> E(P, Q), E(P) + A <--> E(A, P))
        (E + B <--> E(B), E(A) + B <--> E(P, Q))
        (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
        (E + Q <--> E(Q), E(P) + Q <--> E(P, Q)) end end),
]

# Three random-order seeds of A + B + C ⇌ P + Q, larger than the bi-bi laws above: 589,
# 799 and 858 denominator terms.
const _testhelper_large_ter_bi_mechanisms = [
    @enzyme_mechanism(begin substrates: A, B, C; products: P, Q; steps: begin
        (E + A <--> E(A), E(B) + A <--> E(A, B), E(B, C) + A <--> E(P, Q),
         E(C) + A <--> E(A, C), E(P) + A <--> E(A, P))
        (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(A, C) + B ⇌ E(P, Q), E(C) + B ⇌ E(B, C),
         E(C, Q) + B ⇌ E(B, C, Q), E(Q) + B ⇌ E(B, Q))
        (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(A, B) + C ⇌ E(P, Q), E(B) + C ⇌ E(B, C),
         E(B, Q) + C ⇌ E(B, C, Q), E(Q) + C ⇌ E(C, Q))
        (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
        (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(B, C) + Q <--> E(B, C, Q),
         E(C) + Q <--> E(C, Q), E(P) + Q <--> E(P, Q)) end end),
    @enzyme_mechanism(begin substrates: A, B, C; products: P, Q; steps: begin
        (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(B, C) + A ⇌ E(P, Q), E(C) + A ⇌ E(A, C),
         E(P) + A ⇌ E(A, P))
        (E + B <--> E(B), E(A) + B <--> E(A, B), E(A, C) + B <--> E(P, Q),
         E(C) + B <--> E(B, C), E(C, Q) + B <--> E(B, C, Q), E(Q) + B <--> E(B, Q))
        (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(A, B) + C ⇌ E(P, Q), E(B) + C ⇌ E(B, C),
         E(B, Q) + C ⇌ E(B, C, Q), E(Q) + C ⇌ E(C, Q))
        (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
        (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(B, C) + Q <--> E(B, C, Q),
         E(C) + Q <--> E(C, Q), E(P) + Q <--> E(P, Q)) end end),
    @enzyme_mechanism(begin substrates: A, B, C; products: P, Q; steps: begin
        (E + A <--> E(A), E(C) + A <--> E(Q; residual = A + C - Q),
         E(C, P) + A <--> E(A, C, P), E(P) + A <--> E(A, P))
        (E + C <--> E(C), E(A) + C <--> E(Q; residual = A + C - Q),
         E(A, P) + C <--> E(A, C, P), E(P) + C <--> E(C, P))
        (E + P <--> E(P), E(A) + P <--> E(A, P), E(C) + P <--> E(C, P),
         E(Q; residual = A + C - Q) + P <--> E(A, C, P))
        (E(P) + Q ⇌ E(B, Q; residual = A + C - Q),
         E(; residual = A + C - Q) + Q ⇌ E(Q; residual = A + C - Q))
        (E(Q; residual = A + C - Q) + B ⇌ E(B, Q; residual = A + C - Q),
         E(; residual = A + C - Q) + B ⇌ E(P)) end end),
]

# Seconds per call allowed per denominator term: each of the six laws above must run in
# under the larger of 120 ns and this constant times its denominator terms. It covers
# the worst CI readings on each runner, in ns per term, ter-bi / bi-bi: ubuntu x64
# 0.305 / 0.313, windows x64 0.346 / 0.459, macOS aarch64 0.248 / 0.732.
const _testhelper_seconds_per_denominator_term = 0.8e-9

"""The number of terms in the denominator of `m`'s derived rate law."""
_testhelper_denominator_terms(m) =
    length(EnzymeRates._raw_symbolic_rate_polys(typeof(m))[2])

@testset "rate_equation time per call grows at most in proportion to the law's size" begin
    @test _testhelper_denominator_terms.(_testhelper_largest_enumerated_mechanisms) ==
          [396, 396, 396]
    @test _testhelper_denominator_terms.(_testhelper_large_ter_bi_mechanisms) ==
          [589, 799, 858]
    for m in vcat(_testhelper_largest_enumerated_mechanisms,
                  _testhelper_large_ter_bi_mechanisms)
        params, concs, _ = _testhelper_random_independent_params_concs(
            m, collect(metabolites(m)); rng = Random.MersenneTwister(42))
        allocs, t = _testhelper_test_rate_equation_performance(m, params, concs)
        n = _testhelper_denominator_terms(m)
        @info "rate_equation: $(round(1e9 * t; digits = 1)) ns per call, $n denominator " *
              "terms, $(round(1e9 * t / n; digits = 3)) ns per term"
        @test allocs == 0
        @test t < max(120e-9, _testhelper_seconds_per_denominator_term * n)
    end
end

"""The summands of the sum `e`, through `+`, `-` and unary minus."""
function _testhelper_summands(e)
    e isa Expr && e.head == :call && e.args[1] in (:+, :-) || return Any[e]
    reduce(vcat, (_testhelper_summands(a) for a in e.args[2:end]))
end

"""The factors of the product `e`, through `*` and the numerator of `/`; a power counts
as its base."""
function _testhelper_factors(e)
    e isa Expr && e.head == :call || return Any[e]
    op = e.args[1]
    op == :* && return reduce(vcat, (_testhelper_factors(a) for a in e.args[2:end]))
    op == :/ && return _testhelper_factors(e.args[2])
    op == :^ && return Any[e.args[2]]
    op == :- && length(e.args) == 2 && return _testhelper_factors(e.args[2])
    Any[e]
end

"""The sums in `e` in which some metabolite of `mets` divides more than one summand."""
function _testhelper_unfactored_sums(e, mets)
    e isa Expr || return Any[]
    found = reduce(vcat, (_testhelper_unfactored_sums(a, mets) for a in e.args);
                   init = Any[])
    e.head == :call && e.args[1] in (:+, :-) || return found
    shared(x) = count(s -> x in _testhelper_factors(s), _testhelper_summands(e)) > 1
    any(shared, mets) ? push!(found, e) : found
end

@testset "the generated rate_equation body is in Horner form" begin
    # Horner form multiplies by each factor once per branch instead of once per term,
    # which keeps the 396-term laws above at 54-62 ns per call locally. Greedy Horner
    # takes a factor out of every summand it divides, so no concentration divides two
    # summands of one sum.
    for m in (_testhelper_largest_enumerated_mechanisms[2],
              only(s for s in MECHANISM_TEST_SPECS
                   if s.name == "MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites"
                  ).mechanism)
        @test isempty(_testhelper_unfactored_sums(
            _testhelper_generated_rate_return(m), metabolites(m)))
    end
end

@testset "rate_equation_string prints every tie as a constraint line" begin
    # A binds E and E(B) under one constant while B's two bindings keep their own, so
    # the RE binding square ties B's two dissociation constants. The occupied binding's
    # constant is the dependent one, and v reads it.
    folded = @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
            E + B ⇌ E(B)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    @test occursin("# Wegscheider constraints:\nK_EAB_to_EA_B = K_EB_to_E_B\n",
                   rate_equation_string(folded))
    _testhelper_check_law_against_oracle(folded)
    # S binds E and Estar under one constant, so the square of the two isomerizations
    # ties their equilibrium constants; the isomerization between occupied forms is the
    # dependent one.
    iso_tie = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E + S ⇌ E(S), Estar + S ⇌ Estar(S))
            E ⇌ Estar
            E(S) ⇌ Estar(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    @test occursin("\nK_ES_to_EstarS = K_E_to_Estar\n", rate_equation_string(iso_tie))
    _testhelper_check_law_against_oracle(iso_tie)
    # I binding E and I binding Estar share one group and run opposite ways around the
    # cycle, so the Haldane relation fixes the Theorell–Chance constant at Keq; v reads it.
    keq_tie = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ Estar + P
            (E + I <--> E(I), Estar + I <--> E(I))
        end
    end
    @test occursin("\nK_E_S_to_Estar_P = Keq\n", rate_equation_string(keq_tie))
    _testhelper_check_law_against_oracle(keq_tie)
end

@testset "rate_equation_string prints residual names that run as written" begin
    # Every step of this ping-pong is at steady state, so the Haldane product holds the
    # rate constants of the residual-bearing forms. Their names are Julia identifiers, so
    # the printed law, run as written, equals the compiled law, which the oracle checks.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B <--> E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) <--> E + Q
        end
    end
    s = rate_equation_string(m)
    @test occursin(
        r"\nk_EQ_to_EB_res_plusA_minusP = .* \* k_EA_to_EP_res_plusA_minusP \* ", s)
    @test !occursin("var\"", s)
    params = _testhelper_random_independent_params(m; rng = Random.MersenneTwister(6))
    concs = (A = 0.7, B = 1.3, P = 0.4, Q = 0.9)
    @test isapprox(_testhelper_eval_rate_string(s, params, concs),
                   rate_equation(m, concs, params); rtol = 1e-12)
    _testhelper_check_law_against_oracle(m)
end

@testset "rate_equation_string runs a half-integer dependent as written" begin
    # The conformational square ties two isomerizations of one group, so the
    # Haldane relation fixes a rate constant through square roots. The printed law,
    # run as written, equals the compiled law, which the oracle checks.
    half = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E(S) <--> F(S), E <--> F)
            F(S) <--> E + P
            F + S <--> E(S)
        end
    end
    s = rate_equation_string(half)
    @test occursin("0.5", s)
    rng = Random.MersenneTwister(4)
    point = setprecision(() -> _testhelper_consistent_point(half, rng), BigFloat, 256)
    fitted = EnzymeRates._dependent_param_exprs(typeof(half))[2]
    params = NamedTuple{(fitted..., :Keq, :E_total)}(
        (Float64.(getindex.(Ref(point.vals), fitted))..., Float64(point.Keq), 1.0))
    concs = (S = 0.7, P = 1.3)
    @test isapprox(_testhelper_eval_rate_string(s, params, concs),
                   rate_equation(half, concs, params); rtol = 1e-12)
end

@testset "rate_equation_string prints flat +/* sums (no nested parens)" begin
    # Orthogonal guard (NOT a perf-fix backstop): _expr_to_string is
    # precedence-aware and flattens nested +/* nodes
    # transparently: balanced +(+(a,b), +(c,d)) prints as "a + b + c + d"
    # with no added parens. If that flattening regresses, the printed
    # output would gain nested parens but still evaluate correctly —
    # easy to miss without an explicit guard. The `v =` line of Random-order
    # Bi-Bi holds three: the numerator, the denominator and the negative terms
    # inside the numerator.
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "Random-order Bi-Bi")
    s = last(split(rate_equation_string(spec.mechanism), '\n'))
    @test startswith(s, "v = ")
    @test count(==('('), s) <= 3
end

@testset "_ss_rate_constant_names" begin
    # SS-only Uni-Uni: 2 SS binding steps → kon/koff names are SS.
    uni_uni = only(s for s in MECHANISM_TEST_SPECS
                   if s.name == "Uni-Uni").mechanism
    names = EnzymeRates._ss_rate_constant_names(uni_uni)
    for sym in (:k_E_S_to_ES, :k_ES_to_E_S, :k_ES_to_E_P, :k_E_P_to_ES)
        @test sym in names
    end

    # Mixed RE/SS: RE binding (K_EA_to_E_A) + SS catalysis. Only the
    # SS k's are returned; RE binding K is excluded.
    re_uu = only(s for s in MECHANISM_TEST_SPECS
                 if s.name == "RE Uni-Uni").mechanism
    re_uu_names = EnzymeRates._ss_rate_constant_names(re_uu)
    @test :k_EA_to_E_P in re_uu_names && :k_E_P_to_EA in re_uu_names
    for sym in (:K_EA_to_E_A, :Keq, :L, :E_total)
        @test !(sym in re_uu_names)
    end

    # Allosteric: I-state versions of every SS rate constant are also
    # included so `rescale_parameter_values` scales them in tandem.
    mwc = only(s for s in MECHANISM_TEST_SPECS
               if s.name == "MWC Dimer [AllostericEnzymeMechanism]").mechanism
    mwc_names = EnzymeRates._ss_rate_constant_names(mwc)
    for sym in (:k_A_ES_to_EP, :k_A_EP_to_ES, :k_I_ES_to_EP, :k_I_EP_to_ES)
        @test sym in mwc_names
    end
    for sym in (:K_A_ES_to_E_S, :K_A_EP_to_E_P, :K_I_ES_to_E_S, :K_I_EP_to_E_P,
                :Keq, :L, :E_total)
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

        # No factors and zero Keq: returns the Float64 literal 1.0
        @test bpe(R(0), Tuple{Symbol, R}[]) === 1.0

        # Keq times a factor is a valid AST node
        @test bpe(R(1), [(:k1f, R(1))]) isa Union{Float64, Symbol, Expr}
    end

    @testset "a constant thermodynamics pins to 1 is the dependent constant 1.0" begin
        # I binds E and F under one dissociation constant, so the cycle E → E(I) → F → E
        # forces K_E_to_F = 1 and leaves no other constant in its row.
        em = @enzyme_mechanism begin
            substrates: A
            products: P
            modifiers: I
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(P)
                E(P) <--> E + P
                E ⇌ F
                (E + I ⇌ E(I), F + I ⇌ E(I))
            end
        end
        @test occursin("# Wegscheider constraints:\nK_E_to_F = 1.0\n",
                       rate_equation_string(em))
        @test :K_E_to_F ∉ EnzymeRates.fitted_params(em)
        _testhelper_check_law_against_oracle(em)
        # The generated rate equation assigns the literal 1.0 and stays allocation-free.
        syms = parameters(em)
        params = NamedTuple{syms}(ntuple(i -> 0.5 + i / 10, length(syms)))
        allocs, _ = _testhelper_test_rate_equation_performance(
            em, params, (A = 1.0, P = 0.5, I = 0.2))
        @test allocs == 0
    end

end

# ── Large equation compilation regression test ────────────────────────────

@testset "Large equation compilation (<20s)" begin
    # A ping-pong bi-bi with a dead-end inhibitor R1: 13 forms and 13 steps, every
    # step at rapid equilibrium except the A → P isomerization.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: R1
        steps: begin
            (E + A ⇌ E(A),
             E(P; residual = A - P) + A ⇌ E(A, P; residual = A - P))
            (E + Q ⇌ E(Q),
             E(B; residual = A - P) + Q ⇌ E(B, Q; residual = A - P))
            (E + R1 ⇌ E(R1),
             E(B; residual = A - P) + R1 ⇌ E(B, R1; residual = A - P),
             E(P; residual = A - P) + R1 ⇌ E(P, R1; residual = A - P))
            (E(A) + P ⇌ E(A, P),
             E(; residual = A - P) + P ⇌ E(P; residual = A - P))
            E(A) <--> E(P; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            (E(Q) + B ⇌ E(B, Q),
             E(; residual = A - P) + B ⇌ E(B; residual = A - P))
        end
    end

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

@testset "All-RE catalytic cycle raises" begin
    # Binding stage mixed (S2 SS, S1 RE), release stage mixed (P1 SS, P2 RE),
    # chemistry RE ⇒ a complete all-RE catalytic cycle exists ⇒ no finite rate.
    m_allre = @enzyme_mechanism begin
        substrates: S1, S2
        products: P1, P2
        steps: begin
            E + S1 ⇌ E(S1)
            E + S2 ⇌ E(S2)
            E(S1) + S2 <--> E(S1, S2)
            E(S2) + S1 ⇌ E(S1, S2)
            E(S1, S2) ⇌ E(P1, P2)
            E(P1, P2) ⇌ E(P1) + P2
            E(P1, P2) <--> E(P2) + P1
            E(P1) ⇌ E + P1
            E(P2) ⇌ E + P2
        end
    end
    # The activator route E(R) + S ⇌ E(S, R) ⇌ E(P, R) ⇌ E(R) + P is an all-RE
    # catalytic cycle, so the mechanism has no finite rate.
    m_sib = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: R
        steps: begin
            E + S ⇌ E(S)
            E(S) ⇌ E(P)
            E + P <--> E(P)
            E(R) + S ⇌ E(S, R)
            E(S, R) ⇌ E(P, R)
            E(R) + P ⇌ E(P, R)
            (E + R ⇌ E(R),
             E(S) + R ⇌ E(S, R),
             E(P) + R ⇌ E(P, R))
        end
    end
    for m in (m_allre, m_sib)
        err = try
            rate_equation_string(m); nothing
        catch e; e end
        @test err isa ErrorException
        @test occursin("no finite rate", err.msg)
    end
end

@testset "_eq_complexity (V×τ term-count estimate)" begin
    raw(m) = m isa EnzymeRates.Mechanism ? m : EnzymeRates.Mechanism(m)
    # V×τ = the King–Altman denominator term count, computed pre-derivation from
    # the catalytic segment graph. Values are the exact denominator product counts.
    ordered_bibi = @enzyme_mechanism begin
        substrates: S1, S2
        products: P1, P2
        steps: begin
            E + S1 <--> E(S1)
            E(S1) + S2 <--> E(S1, S2)
            E(S1, S2) <--> E(P1, P2)
            E(P1, P2) <--> E(P1) + P2
            E(P1) <--> E + P1
        end
    end
    @test EnzymeRates._eq_complexity(raw(ordered_bibi)) == 25

    random_bibi = @enzyme_mechanism begin
        substrates: S1, S2
        products: P1, P2
        steps: begin
            E + S1 <--> E(S1)
            E + S2 <--> E(S2)
            E(S1) + S2 <--> E(S1, S2)
            E(S2) + S1 <--> E(S1, S2)
            E(S1, S2) <--> E(P1, P2)
            E(P1, P2) <--> E(P1) + P2
            E(P1, P2) <--> E(P2) + P1
            E(P1) <--> E + P1
            E(P2) <--> E + P2
        end
    end
    @test EnzymeRates._eq_complexity(raw(random_bibi)) == 336

    ordered_terter = @enzyme_mechanism begin
        substrates: S1, S2, S3
        products: P1, P2, P3
        steps: begin
            E + S1 <--> E(S1)
            E(S1) + S2 <--> E(S1, S2)
            E(S1, S2) + S3 <--> E(S1, S2, S3)
            E(S1, S2, S3) <--> E(P1, P2, P3)
            E(P1, P2, P3) <--> E(P1, P2) + P3
            E(P1, P2) <--> E(P1) + P2
            E(P1) <--> E + P1
        end
    end
    @test EnzymeRates._eq_complexity(raw(ordered_terter)) == 49

    # Allosteric: computed on the catalytic core, so it equals the catalytic V×τ
    # (an all-SS uni-uni: G = 3 segments, τ = 3 spanning trees ⇒ 9).
    allo = @allosteric_mechanism begin
        substrates: F6P
        products: F16BP
        catalytic_multiplicity: 2
        allosteric_regulators: I
        steps: begin
            E + F6P <--> E(F6P)      :: EqualAI
            E(F6P) <--> E(F16BP)     :: EqualAI
            E(F16BP) <--> E + F16BP  :: EqualAI
        end
        regulatory_site(2): I::OnlyI
    end
    @test EnzymeRates._eq_complexity(EnzymeRates.AllostericMechanism(allo)) == 9
end

@testset "Rate equation too large error" begin
    # Manually defined mechanism (11 forms, 16 steps; V×τ ≈ 29k denominator
    # terms) exceeds MAX_RATE_EQUATION_TERMS, so the upfront V×τ check
    # (_assert_derivable) aborts the derivation before it builds the polynomial.
    m_manual = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: R1
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
end

# The generated rate_equation reads its mechanism type and its arguments under names of
# its own; a metabolite may take the names of the type parameter and the first argument.
@testset "a metabolite named M or m reads as a concentration" begin
    upper = @enzyme_mechanism begin
        substrates: M
        products:   P
        steps: begin
            E + M ⇌ E(M)
            E(M) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    _testhelper_check_law_against_oracle(upper)
    lower = @enzyme_mechanism begin
        substrates: m
        products:   P
        steps: begin
            E + m ⇌ E(m)
            E(m) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    _testhelper_check_law_against_oracle(lower)
end

# The rate equation binds its constants and the concentrations as local names of one
# scope, so no metabolite takes a constant's name: a reserved name, or one holding `_`,
# is rejected when a mechanism macro's code runs.
@testset "a metabolite named after a constant is rejected" begin
    # An allosteric regulator named L would take the conformational constant's name.
    @test_throws(ErrorException(
        "AllostericMechanism: L is a reserved name; E_total, Keq, L, N, params, concs " *
        "and v are kept for a rate equation's constants and local names"),
        @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: L
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)     :: EqualAI
                E(S) <--> E(P)   :: OnlyA
                E(P) ⇌ E + P     :: EqualAI
            end
            regulatory_site(2): L::OnlyI
        end)
    # An inhibitor named after the reverse chemistry constant holds `_`.
    @test_throws(ErrorException(
        "EnzymeReaction: the metabolite name k_EP_to_ES holds `_`; constant names join " *
        "metabolite names with `_`, so a metabolite name holds none"),
        @enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: k_EP_to_ES
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                E + k_EP_to_ES ⇌ E(k_EP_to_ES)
            end
        end)
end

# Every constant of a rate equation needs a name of its own. An :EqualAI ligand's
# constant is bare, K_<ligand>reg, while an :OnlyA or :NonequalAI ligand's names its
# conformation, K_A_<ligand>reg or K_I_<ligand>reg, so a ligand named A_R or I_X would
# take the constant of a ligand R or X. A ligand name holds no `_`.
@testset "a ligand name holding `_` is rejected" begin
    @test_throws(ErrorException(
        "AllostericMechanism: the metabolite name A_R holds `_`; constant names join " *
        "metabolite names with `_`, so a metabolite name holds none"),
        @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            allosteric_regulators: R, A_R
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: NonequalAI
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): R::OnlyA
            regulatory_site(2): A_R::EqualAI
        end)
    @test_throws(ErrorException(
        "AllostericMechanism: the metabolite name I_X holds `_`; constant names join " *
        "metabolite names with `_`, so a metabolite name holds none"),
        @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            allosteric_regulators: X, I_X
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: NonequalAI
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): X::NonequalAI, I_X::EqualAI
        end)

    # Renamed, the :EqualAI ligand has a constant of its own.
    renamed = @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        allosteric_regulators: R, T
        steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: NonequalAI
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::OnlyA
        regulatory_site(2): T::EqualAI
    end
    @test count(==(:K_A_Rreg), parameters(renamed)) == 1
    @test :K_Treg in parameters(renamed)
end

# Every constant of a rate equation needs a name of its own. An :EqualAI step's constants
# are bare, while a :NonequalAI step's name their conformation, k_A_… for the active one
# and k_I_… for the inactive one. A conformation labelled A or I renders the same letter
# at the start of a step side, so an :EqualAI step can take a :NonequalAI step's name.
@testset "constants that share a name raise" begin
    # The uptake constant of the :EqualAI step A(Eb) <--> A + Eb is k_A_Eb_to_AEb, the
    # name of the active constant of the :NonequalAI step E(b) <--> A(Eb).
    conformation_a = @allosteric_mechanism begin
        substrates: b
        products:   Eb
        steps: begin
            E + b ⇌ E(b)          :: EqualAI
            E(b) <--> A(Eb)       :: NonequalAI
            A(Eb) <--> A + Eb     :: EqualAI
            A ⇌ E                 :: EqualAI
        end
    end
    clash_a = ErrorException("Constants of this mechanism's rate equation share a name " *
                             "(k_A_Eb_to_AEb); rename a metabolite or a conformation " *
                             "label so that each constant has a name of its own")
    params = (K_Eb_to_E_b = 1.0, K_A_to_E = 2.0, k_A_Eb_to_AEb = 3.0, k_AEb_to_A_Eb = 4.0,
              k_I_Eb_to_AEb = 0.5, L = 10.0, Keq = 4.0, E_total = 1.0)
    @test_throws clash_a parameters(conformation_a)
    @test_throws clash_a EnzymeRates.fitted_params(conformation_a)
    @test_throws clash_a rate_equation(conformation_a, (b = 1.0, Eb = 0.1), params)
    @test_throws clash_a rate_equation_string(conformation_a)
    @test_throws clash_a all_parameters(conformation_a, params)
    @test_throws clash_a EnzymeRates._independent_param_count(
        EnzymeRates.AllostericMechanism(conformation_a))

    # Labelled I, the :EqualAI uptake constant k_I_Eb_to_IEb is the name of the
    # :NonequalAI step's inactive constant.
    conformation_i = @allosteric_mechanism begin
        substrates: b
        products:   Eb
        steps: begin
            E + b ⇌ E(b)          :: EqualAI
            E(b) <--> I(Eb)       :: NonequalAI
            I(Eb) <--> I + Eb     :: EqualAI
            I ⇌ E                 :: EqualAI
        end
    end
    clash_i = ErrorException("Constants of this mechanism's rate equation share a name " *
                             "(k_I_Eb_to_IEb); rename a metabolite or a conformation " *
                             "label so that each constant has a name of its own")
    @test_throws clash_i parameters(conformation_i)
    @test_throws clash_i rate_equation_string(conformation_i)

    # Labelled F, each constant has a name of its own.
    renamed = @allosteric_mechanism begin
        substrates: b
        products:   Eb
        steps: begin
            E + b ⇌ E(b)          :: EqualAI
            E(b) <--> F(Eb)       :: NonequalAI
            F(Eb) <--> F + Eb     :: EqualAI
            F ⇌ E                 :: EqualAI
        end
    end
    @test parameters(renamed) == (:K_Eb_to_E_b, :K_F_to_E, :k_A_Eb_to_FEb, :k_FEb_to_F_Eb,
                                  :k_F_Eb_to_FEb, :k_I_Eb_to_FEb, :L, :Keq, :E_total)
end

"""Every name the expression `ex` binds or reads: each symbol except a called
function's."""
function _testhelper_names_read(ex)
    ex isa Symbol && return Set([ex])
    ex isa Expr || return Set{Symbol}()
    union(Set{Symbol}(),
          (_testhelper_names_read(a) for a in (ex.head === :call ? ex.args[2:end] :
                                                ex.args))...)
end

# The names no metabolite takes: E_total, Keq, L, N, params, concs, v and every name
# holding `_` (test_types.jl checks that the reaction rejects each). Besides the
# metabolites, the generated rate_equation body and the printed law bind or read only
# such names, so no concentration can overwrite a constant or a local.
@testset "the reserved names cover every non-metabolite name of a law" begin
    reserved = (:E_total, :Keq, :L, :N, :params, :concs, :v)
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        dep, indep = EnzymeRates._dependent_param_exprs(typeof(m))
        body = EnzymeRates._rate_body(typeof(m), (indep..., :Keq, :E_total), dep)
        printed = Meta.parseall(rate_equation_string(m))
        others = setdiff(_testhelper_names_read(body) ∪ _testhelper_names_read(printed),
                         metabolites(m))
        @test all(n -> occursin('_', String(n)) || n in reserved, others)
        # The law destructures params and concs and reads Keq; the printed one assigns v.
        @test (:params, :concs, :Keq) ⊆ others
        @test :v in _testhelper_names_read(printed)
    end
end


# ── Single-feature edge cases ─────────────────────────────────────────────
@testset "Allosteric edge cases" begin
    # OnlyA substrate + OnlyA catalysis: S binds only in the R-state and only
    # the R-state catalyzes (k_T = 0). The T-state therefore carries no flux
    # (N_T = 0) and never populates E(S), so no reverse-catalysis weight leaks
    # into Q_I. As K_A_ES_to_E_S → ∞ (weaker R-state binding), rate vanishes.
    onlyR_sub = @allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)     :: OnlyA
            E(S) <--> E(P)   :: OnlyA
            E(P) ⇌ E + P     :: EqualAI
        end
    end
    concs = (S=1.0, P=0.001)
    base_params = (k_A_ES_to_EP=10.0, K_EP_to_E_P=0.5, L=10.0, Keq=1000.0, E_total=1.0)
    rate_strong = rate_equation(onlyR_sub, concs, merge(base_params, (K_A_ES_to_E_S=0.01,)))
    rate_weak   = rate_equation(onlyR_sub, concs, merge(base_params, (K_A_ES_to_E_S=1e6,)))
    @test rate_strong > 1.0
    @test rate_weak < 1e-3
    @test rate_weak / rate_strong < 1e-5

    # V-type only: all bindings :EqualAI but catalysis :OnlyA. T-state binds
    # substrate normally but cannot catalyze (k_T = 0), so N_T = 0. As L → ∞
    # (T-state dominant) the rate vanishes.
    vtype = @allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)     :: EqualAI
            E(S) <--> E(P)   :: OnlyA
            E(P) ⇌ E + P     :: EqualAI
        end
    end
    vparams = (K_ES_to_E_S=0.1, k_A_ES_to_EP=10.0, K_EP_to_E_P=0.5, Keq=1000.0, E_total=1.0)
    rate_R = rate_equation(vtype, concs, merge(vparams, (L=0.0,)))
    rate_T = rate_equation(vtype, concs, merge(vparams, (L=1e10,)))
    @test rate_R > 1.0
    @test rate_T < 1e-6
    # T-state numerator branch is elided when the T-state cycle is dead (:OnlyA catalysis);
    # rate is E_total · num_R / (Q_R^catN + L · Q_T^catN). At large L, the T-state
    # enzyme mass dominates the denominator → rate ∝ 1/(1+L).
    @test rate_T * 1e10 < 100.0    # bounded as L grows

    # Regression: T-state binding K's must be in Kd convention even when
    # `:OnlyA` and `:NonequalAI` catalytic groups coexist. Without the fix,
    # the flat-poly path in _num_den_exprs renders T-state K's
    # as `K_T * met` (Ka) instead of `met / K_T` (Kd), silently producing
    # wrong rates whenever a mechanism mixes these two tags. Regression
    # for src/rate_eq_derivation.jl:1395-1396.
    # S binding :NonequalAI, catalysis :OnlyA, P binding :NonequalAI.
    m_mix = @allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: I
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: OnlyA
            E(P) ⇌ E + P      :: NonequalAI
        end
        regulatory_site(2): I::OnlyI
    end
    p_mix = (K_A_ES_to_E_S=0.1, k_A_ES_to_EP=10.0, K_A_EP_to_E_P=0.5,
             K_I_ES_to_E_S=10.0, K_I_EP_to_E_P=10.0,
             K_I_Ireg=1.0, L=1.0, Keq=1000.0, E_total=1.0)
    rate_mix = rate_equation(m_mix, (S=10.0, P=0.0, I=0.0), p_mix)
    # With Kd convention (correct): rate ≈ 9.90 (R-state catalysis dominates),
    # per active site (E_total = active-site concentration, catN = 2).
    # With Ka convention (bug): rate ≈ 4.95 — half the correct value.
    @test isapprox(rate_mix, 9.90; rtol=0.05)

    # Sanity: rate_equation_string emits Kd form for T-state K's.
    @test occursin("S / K_I_ES_to_E_S", rate_equation_string(m_mix))
    @test occursin("P / K_I_EP_to_E_P", rate_equation_string(m_mix))

    # Wegscheider-cycle EqualAI×NonequalAI: the Random-order Bi-Bi mechanism has a
    # genuine independent Wegscheider cycle. Group 2 (the steady-state B-binding)
    # is :NonequalAI while its box partners are :EqualAI, so its affinity is
    # forbidden and collapses
    # (k_I_EB_to_E_B = k_A_EB_to_E_B·k_I_E_B_to_EB/k_A_E_B_to_EB) while
    # its speed stays free. Asserts the mechanism-agnostic invariant: zero net rate
    # at chemical equilibrium. (Over-parametrized; the enumerator will skip such
    # degenerate configs in a follow-up PR.)
    # B binding (step 2) :NonequalAI, rest :EqualAI.
    m_ro = @allosteric_mechanism begin
        substrates: A, B
        products:   P, Q
        allosteric_regulators: I
        catalytic_multiplicity: 2
        steps: begin
            E + A <--> E(A)            :: EqualAI
            E + B <--> E(B)            :: NonequalAI
            E(A) + B <--> E(A, B)      :: EqualAI
            E(B) + A <--> E(A, B)      :: EqualAI
            E(A, B) <--> E(P, Q)       :: EqualAI
            E(P, Q) <--> E(Q) + P      :: EqualAI
            E(Q) <--> E + Q            :: EqualAI
        end
        regulatory_site(2): I::NonequalAI
    end
    # Overall A + B ⇌ P + Q, so Keq = P·Q/(A·B).
    Keq_ro = 4.0
    A_eq, B_eq = 1.5, 2.0
    P_eq = 3.0
    Q_eq = Keq_ro * A_eq * B_eq / P_eq
    # Which koff is Wegscheider-dependent follows the solve's priority. The values
    # come from the mechanism's own fitted_params, so the Haldane property (zero net
    # rate at equilibrium) is verified whatever the priority picks.
    fp_ro = EnzymeRates.fitted_params(m_ro)
    rng_ro = Random.MersenneTwister(123)
    p_ro = NamedTuple{(fp_ro..., :Keq, :E_total)}(
        (ntuple(_ -> 0.5 + rand(rng_ro), length(fp_ro))..., Keq_ro, 1.0))
    @test isapprox(
        rate_equation(m_ro, (A=A_eq, B=B_eq, P=P_eq, Q=Q_eq, I=0.5), p_ro), 0.0;
        atol=1e-9)
end


@testset "rate_equation_string allosteric byte-identical fixture" begin
    m_allo = @allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + P ⇌ E(P)     :: NonequalAI
            E + S ⇌ E(S)     :: NonequalAI
            E(S) <--> E(P)   :: NonequalAI
        end
        regulatory_site(2): R::NonequalAI
    end
    actual = rate_equation_string(m_allo)
    expected = raw"""(; K_A_EP_to_E_P, K_A_ES_to_E_S, K_I_EP_to_E_P, K_I_ES_to_E_S, k_A_ES_to_EP, k_I_ES_to_EP, K_A_Rreg, K_I_Rreg, L, Keq, E_total) = params
(; S, P, R) = concs
# Haldane constraints:
k_A_EP_to_ES = (1 / Keq) * K_A_EP_to_E_P * (1 / K_A_ES_to_E_S) * k_A_ES_to_EP
k_I_EP_to_ES = (1 / Keq) * K_I_EP_to_E_P * (1 / K_I_ES_to_E_S) * k_I_ES_to_EP
v = E_total * ((k_A_ES_to_EP * S / K_A_ES_to_E_S - k_A_EP_to_ES * P / K_A_EP_to_E_P) * (1 + P / K_A_EP_to_E_P + S / K_A_ES_to_E_S) * (1 + R / K_A_Rreg) ^ 2 + L * (S * k_I_ES_to_EP / K_I_ES_to_E_S - P * k_I_EP_to_ES / K_I_EP_to_E_P) * (1 + P / K_I_EP_to_E_P + S / K_I_ES_to_E_S) * (1 + R / K_I_Rreg) ^ 2) / ((1 + P / K_A_EP_to_E_P + S / K_A_ES_to_E_S) ^ 2 * (1 + R / K_A_Rreg) ^ 2 + L * (1 + P / K_I_EP_to_E_P + S / K_I_ES_to_E_S) ^ 2 * (1 + R / K_I_Rreg) ^ 2)"""
    @test actual == expected
end

@testset "Parameter-struct allosteric helpers" begin
    @testset "_step_names and _kreg_params" begin
        names_of(am, params) = [EnzymeRates.name(p, am) for p in params]
        flat_names(am, state) = reduce(vcat, EnzymeRates._step_names(am, state))
        # :NonequalAI cat groups + :NonequalAI reg ligand → both conformations named.
        aem = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: R
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: NonequalAI
            end
            regulatory_site(1): R::NonequalAI
        end
        am = EnzymeRates.AllostericMechanism(aem)
        # Catalytic side: every step of the inactive law takes its inactive name, the
        # :EqualAI iso its one shared name.
        @test Set(flat_names(am, :I)) ==
              Set([:K_I_ES_to_E_S, :k_ES_to_EP, :k_EP_to_ES, :K_I_EP_to_E_P])
        @test Set(flat_names(am, :A)) ==
              Set([:K_A_ES_to_E_S, :k_ES_to_EP, :k_EP_to_ES, :K_A_EP_to_E_P])
        @test names_of(am, EnzymeRates._kreg_params(am, :A)) == [:K_A_Rreg]
        @test names_of(am, EnzymeRates._kreg_params(am, :I)) == [:K_I_Rreg]

        # Balanced :OnlyA bindings (both S and P) keep the :NonequalAI iso
        # Haldane-valid. The inactive conformation keeps no form but E, so its law has
        # no step and the :OnlyA reg ligand no inactive constant.
        aem_skip = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: R
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: OnlyA
                E(S) <--> E(P)    :: NonequalAI
                E(P) ⇌ E + P      :: OnlyA
            end
            regulatory_site(1): R::OnlyA
        end
        am_skip = EnzymeRates.AllostericMechanism(aem_skip)
        @test isempty(EnzymeRates._step_names(am_skip, :I))
        @test isempty(EnzymeRates._kreg_params(am_skip, :I))
        @test names_of(am_skip, EnzymeRates._kreg_params(am_skip, :A)) == [:K_A_Rreg]

        # An :OnlyI reg ligand has an inactive constant and no active one; an :EqualAI
        # ligand one bare constant, listed with the active ones.
        am_onlyi = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: R
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: NonequalAI
            end
            regulatory_site(1): R::OnlyI
        end)
        @test isempty(EnzymeRates._kreg_params(am_onlyi, :A))
        @test names_of(am_onlyi, EnzymeRates._kreg_params(am_onlyi, :I)) == [:K_I_Rreg]
        am_equal = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: R
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: NonequalAI
            end
            regulatory_site(1): R::EqualAI
        end)
        @test names_of(am_equal, EnzymeRates._kreg_params(am_equal, :A)) == [:K_Rreg]
        @test isempty(EnzymeRates._kreg_params(am_equal, :I))
    end
end

@testset "structural names + synth-dep routing: allosteric NonequalAI" begin
    # NonequalAI substrate binding (PEP), EqualAI catalysis.
    # k_cat_rev (Haldane dep) references the NonequalAI PEP binding K,
    # so an I-state dep name is produced. It must use the structural
    # mid-name I_ token, not string(active) * "_T"; a rate polynomial with
    # structural I-names beside _T-suffixed independent names would make
    # rate_equation error with a KeyError.
    m = @allosteric_mechanism begin
        substrates: PEP, ADP
        products:   Pyruvate, ATP
        allosteric_regulators: ATP, F16BP
        catalytic_multiplicity: 4
        steps: begin
            (E + PEP ⇌ E(PEP),
             E(ADP) + PEP ⇌ E(PEP, ADP))                          :: NonequalAI
            (E + ADP ⇌ E(ADP),
             E(PEP) + ADP ⇌ E(PEP, ADP))                          :: EqualAI
            E(PEP, ADP) <--> E(Pyruvate, ATP)                     :: EqualAI
            (E(Pyruvate, ATP) ⇌ E(ATP) + Pyruvate,
             E(Pyruvate) ⇌ E + Pyruvate)                          :: EqualAI
            (E(Pyruvate, ATP) ⇌ E(Pyruvate) + ATP,
             E(ATP) ⇌ E + ATP)                                    :: EqualAI
        end
        regulatory_site(2): ATP::OnlyI
        regulatory_site(4): F16BP::OnlyA
    end
    ps_str = String.(collect(EnzymeRates.parameters(m)))
    # I-state param names use mid-name I_ token, not _T suffix.
    @test any(s -> startswith(s, "K_I_"), ps_str)
    @test !any(s -> endswith(s, "_T"), ps_str)
    # rate_equation must not error (intermediate state KeyErrors here).
    rng = Random.MersenneTwister(9999)
    met_names = [:PEP, :ADP, :Pyruvate, :ATP, :F16BP]
    concs_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in met_names)
    concs = NamedTuple{Tuple(met_names)}(concs_vals)
    indep = EnzymeRates.fitted_params(m)
    param_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in indep)
    params = NamedTuple{indep}(param_vals)
    params = merge(params, (Keq=1.0, E_total=1.0))
    v = rate_equation(m, concs, params)
    @test isfinite(v)
end

@testset "an allosteric Wegscheider tie leaves the grouped mechanism's model" begin
    # A binds E and E(B) under one shared K while B's two bindings keep their own, so
    # the RE binding square ties B's two K's (K_EAB_to_EA_B = K_EB_to_E_B) in each
    # conformation. The tied constants are dependents, so the mechanism fits the same
    # constants as the one that groups B's two bindings outright, and its law takes the
    # same values.
    tie = @allosteric_mechanism begin
        substrates: A, B
        products: P
        catalytic_multiplicity: 2
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B)) :: NonequalAI
            E + B ⇌ E(B)                       :: NonequalAI
            E(A) + B ⇌ E(A, B)                 :: NonequalAI
            E(A, B) <--> E(P)                  :: NonequalAI
            E(P) ⇌ E + P                       :: NonequalAI
        end
    end
    grouped = @allosteric_mechanism begin
        substrates: A, B
        products: P
        catalytic_multiplicity: 2
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B)) :: NonequalAI
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B)) :: NonequalAI
            E(A, B) <--> E(P)                  :: NonequalAI
            E(P) ⇌ E + P                       :: NonequalAI
        end
    end
    @test EnzymeRates.fitted_params(tie) == EnzymeRates.fitted_params(grouped) ==
        (:K_A_EA_to_E_A, :K_A_EB_to_E_B, :K_A_EP_to_E_P, :K_I_EA_to_E_A,
         :K_I_EB_to_E_B, :K_I_EP_to_E_P, :k_A_EAB_to_EP, :k_I_EAB_to_EP, :L)
    rng = Random.MersenneTwister(4)
    for _ in 1:5
        params = _testhelper_random_independent_params(tie; rng)
        concs = (A = 0.5 + rand(rng), B = 0.5 + rand(rng), P = 0.5 + rand(rng))
        @test rate_equation(tie, concs, params) ≈
              rate_equation(grouped, concs, params) rtol = 1e-12
    end
end

@testset "Fix A: dead-inactive-state allosteric body defines all I-state symbols" begin
    # Random-order allosteric bi-bi with an :OnlyA catalytic step → dead inactive
    # state. Verified pre-fix to crash with an UndefVarError on the I-state
    # release rate of A from E(A).
    m = @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 2
        steps: begin
            E + A <--> E(A)        :: NonequalAI
            E + B <--> E(B)        :: NonequalAI
            E(A) + B <--> E(A, B)  :: NonequalAI
            E(B) + A <--> E(A, B)  :: NonequalAI
            E(A, B) <--> E(P, Q)   :: OnlyA
            E(P, Q) <--> E(P) + Q  :: NonequalAI
            E(P, Q) <--> E(Q) + P  :: NonequalAI
            E(P) <--> E + P        :: NonequalAI
            E(Q) <--> E + Q        :: NonequalAI
        end
    end
    pn = EnzymeRates.fitted_params(m)
    params = merge(NamedTuple{pn}(ntuple(_ -> 1.0, length(pn))),
                   (Keq = 1.0, E_total = 1.0))
    concs = (A = 1.0, B = 1.0, P = 1.0, Q = 1.0)
    @test isfinite(rate_equation(m, concs, params))
end

@testset "the peak of an allosteric random-order bi-bi reaches the grid maximum" begin
    # Random-order allosteric bi-bi, live I-state.
    m = @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 2
        steps: begin
            E + A <--> E(A)        :: NonequalAI
            E + B <--> E(B)        :: NonequalAI
            E(A) + B <--> E(A, B)  :: NonequalAI
            E(B) + A <--> E(A, B)  :: NonequalAI
            E(A, B) <--> E(P, Q)   :: NonequalAI
            E(P, Q) <--> E(P) + Q  :: NonequalAI
            E(P, Q) <--> E(Q) + P  :: NonequalAI
            E(P) <--> E + P        :: NonequalAI
            E(Q) <--> E + Q        :: NonequalAI
        end
    end
    rng = Random.MersenneTwister(1)
    pn = EnzymeRates.fitted_params(m)
    pv = NamedTuple{pn}(Tuple(0.2 + 2 * rand(rng) for _ in pn))
    peak = EnzymeRates._peak_turnover(m, merge(pv, (Keq = 1.0,)), (1e-12, 1e9))
    @test isfinite(peak.value)
    fp = merge(pv, (Keq = 1.0, E_total = 1.0))
    vgrid = maximum(rate_equation(m, (A = x, B = y, P = p, Q = q), fp)
                    for x in 10.0 .^ (-12:9), y in 10.0 .^ (-12:9),
                        p in 10.0 .^ (-12:3:9), q in 10.0 .^ (-12:3:9))
    @test peak.value >= vgrid * (1 - 1e-12)
end

@testset "the peak turnover covers every regulator corner at saturating substrate" begin
    # A regulator that binds one conformation only (:OnlyA or :OnlyI) puts all the enzyme
    # in that conformation when it saturates, so the rate tends to that conformation's own
    # saturating turnover; with the regulator absent both conformations contribute. The
    # peak turnover is at least every regulator corner's turnover at saturating S.
    # :NonequalAI catalysis: both conformations turn over.
    only_a = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: X
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: NonequalAI
            E(P) ⇌ E + P      :: NonequalAI
        end
        regulatory_site(2): X::OnlyA
    end
    only_i = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: X
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: NonequalAI
            E(P) ⇌ E + P      :: NonequalAI
        end
        regulatory_site(2): X::OnlyI
    end
    # :OnlyA catalysis with :EqualAI binding: the inactive conformation binds S but does
    # not turn over, so saturating X drives the rate to 0.
    only_i_dead = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: X
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: OnlyA
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): X::OnlyI
    end
    # Two sites: with X and Y saturating, the inactive conformation's regulator power is
    # 2 + 1 against the active one's 1.
    two_site = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: X, Y
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: OnlyA
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): X::OnlyI
        regulatory_site(1): Y::NonequalAI
    end
    BIG = 1e9
    for m in (only_a, only_i, only_i_dead, two_site)
        pn = EnzymeRates.fitted_params(m)
        regs = Tuple(x for x in metabolites(m) if x ∉ (:S, :P))
        for seed in 1:5
            rng = Random.MersenneTwister(seed)
            p = merge(NamedTuple{pn}(Tuple(0.2 + 5 * rand(rng) for _ in pn)),
                      (Keq = 2.0, E_total = 1.0))
            # Corners of the searched range: each regulator at a bound, saturating S,
            # the product at the lower bound.
            vmax = maximum(0:(2^length(regs) - 1)) do mask
                r = Tuple((mask >> (i - 1)) & 1 == 1 ? BIG : 1e-12
                          for i in eachindex(regs))
                rate_equation(m, NamedTuple{(:S, :P, regs...)}((BIG, 1e-12, r...)), p)
            end
            @test EnzymeRates._peak_turnover(m, p, (1e-12, BIG)).value >=
                  vmax * (1 - 1e-12)
        end
    end
end

@testset "the peak turnover of a random-order bi-bi reaches the grid maximum" begin
    # 40 parameter draws; the grid spans every metabolite, products included.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E + B <--> E(B)
            E(A) + B <--> E(A, B)
            E(B) + A <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(P) + Q
            E(P, Q) <--> E(Q) + P
            E(P) <--> E + P
            E(Q) <--> E + Q
        end
    end
    pn = EnzymeRates.fitted_params(m)
    for seed in 1:40
        rng = Random.MersenneTwister(seed)
        pv = NamedTuple{pn}(Tuple(0.2 + 5 * rand(rng) for _ in pn))
        peak = EnzymeRates._peak_turnover(m, merge(pv, (Keq = 1.0,)), (1e-12, 1e3))
        fp = merge(pv, (Keq = 1.0, E_total = 1.0))
        vgrid = maximum(rate_equation(m, (A = x, B = y, P = p, Q = q), fp)
                        for x in 10.0 .^ (-12:3), y in 10.0 .^ (-12:3),
                            p in 10.0 .^ (-12:3:3), q in 10.0 .^ (-12:3:3))
        @test peak.value >= vgrid * (1 - 1e-12)
    end
end

@testset "the peak of a law that needs bound product lies at a positive product level" begin
    # Catalysis runs only on the product-bound form E(P), so every forward numerator
    # term carries P and the turnover vanishes at zero product.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E(P)
            E(P) + S ⇌ E(P, S)
            E(P, S) <--> E(P, P)
            E(P, P) ⇌ E(P) + P
        end
    end
    pn = EnzymeRates.fitted_params(m)
    params = merge(NamedTuple{pn}(ntuple(_ -> 1.5, length(pn))), (Keq = 2.0,))
    peak = EnzymeRates._peak_turnover(m, params, (1e-12, 1e3))
    @test peak.value > 0
    @test 1e-3 < peak.concs.P < 1e3
end

@testset "rendering helpers" begin
    k = ER.POLY(ER._mono(:k_ES_to_EP => 1) => 1)
    @test ER._invert_monomial(k) == ER.POLY(ER._mono(:k_ES_to_EP => -1) => 1)
    @test ER._invert_monomial(ER.poly_one()) == ER.poly_one()
    @test ER._is_metabolite_free_monomial(k, Set([:S, :P]))
    kb = ER.poly_mul(k, ER.POLY(ER._mono(:B => 1) => 1))          # k * B (has metabolite)
    @test !ER._is_metabolite_free_monomial(kb, Set([:B]))
    twoterm = ER.poly_add(k, ER.poly_one())                       # k + 1 (not a monomial)
    @test !ER._is_metabolite_free_monomial(twoterm, Set([:S]))
    @test_throws ErrorException ER._invert_monomial(twoterm)      # non-monomial errors
    @test ER._mwc_cross_weight(:foo, 1, 2) == :foo               # no-op when D==1
    @test ER._mwc_cross_weight(:foo, :D, 1) == :(D * foo)
end
