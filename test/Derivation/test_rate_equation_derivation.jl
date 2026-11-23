# ABOUTME: Tests of src/Derivation/rate_equation_derivation.jl: every spec's law and its
# ABOUTME: speed gate, the engine, polynomials, the Sig, the compiled and printed laws.

# ── Every spec's law, the compiled law and its speed ──────────────────────────────────

using OrdinaryDiffEqFIRK

# ── Helper functions ────────────────────────────────────────────────────────

# ── Reference QSSA implementation ───────────────────────────────────────────

"""
Independent reference: the QSSA rate per unit enzyme by the Laplacian cofactor method,
walking the steps of the plain mechanism `m`.
"""
function _testhelper_reference_qssa(
    m::EnzymeRates.Mechanism,
    params::NamedTuple,
    concs::NamedTuple,
)
    g = _testhelper_graph(m)
    n = length(g.forms)

    ref_name, nu_ref = _testhelper_reference_metabolite(m)

    # Build rate matrix R[i,j] = pseudo-first-order rate from i to j
    R = zeros(n, n)
    for (step_idx, step) in enumerate(g.steps)
        i, j = step.from, step.to

        m_lhs = step.consumed
        m_rhs = step.released

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

    E_conc = D ./ sum(D)

    # Compute net consumption of reference substrate
    v = 0.0
    for (step_idx, step) in enumerate(g.steps)
        i, j = step.from, step.to

        m_lhs = step.consumed
        m_rhs = step.released

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
matches by ligand-name set so positional regulator names (`:K_<lig>_reg{site}`) keep
the oracle's written site numbering after canonical site reordering.
"""
function _testhelper_positional_site_idx(site, pos, written_sites)
    want = Set(site.ligands)
    j = findfirst(ws -> Set(ws.ligands) == want, written_sites)
    j === nothing ? pos : j
end

"""
The constants of the plain mechanism `m` (values by name in `all_params`, `params` merged
with `dependent_parameter_values(m, params, Keq)`) under the positional names the QSSA and
ODE cross-checks read: stored flat step `i` gives `K{i}` at rapid equilibrium and `k{i}f`,
`k{i}r` at steady state, in its stored direction. Every entry of `all_params` passes
through.
"""
function _testhelper_stored_params(m, all_params::NamedTuple)
    out = Pair{Symbol, Any}[]
    for (i, s) in enumerate(_testhelper_graph(m).steps)
        if s.rapid_equilibrium
            push!(out, Symbol("K", i) => all_params[only(s.names)])
        else
            push!(out, Symbol("k", i, "f") => all_params[s.names[1]],
                  Symbol("k", i, "r") => all_params[s.names[2]])
        end
    end
    merge((; out...), all_params)
end

"""
The constants of `m` (values by name in `all_params`, `params` merged with
`dependent_parameter_values(m, params, Keq)`) under the positional names a textbook rate
law reads, by role. Written step `i` (in the order of `raw`'s written steps) takes the
constants of the step it became (`_testhelper_stored_steps`), in the written direction:
`k{i}f` and `k{i}r` at steady state, and at rapid equilibrium `K{i}`, a binding's
dissociation constant or, for any other step, the equilibrium constant of the written
direction. The inactive law's kept steps in a `:NonequalAI` group add `k{i}f_T`,
`k{i}r_T` or `K{i}_T`. No spec in `MECHANISM_TEST_SPECS` writes a non-binding step against
its stored direction, so the suite never inverts `K` or swaps `k{i}f` and `k{i}r` for such a
step. A ligand of regulatory site `j` (the written order of `raw`'s sites) gives
`K_<ligand>_reg{j}` and, when its inactive conformation has a constant of its own,
`K_<ligand>_T_reg{j}`. `raw` is the mechanism as written, the DSL walker's value
(`_testhelper_MechanismTestSpec.raw`); errors unless its catalytic steps canonicalize to
`m`'s catalytic graph. Every entry of `all_params` passes through.
"""
function _testhelper_written_params(m, all_params::NamedTuple; raw)
    allo = _testhelper_is_allosteric(m)
    g = _testhelper_graph(m)
    _testhelper_build(raw).catalytic == m.catalytic ||
        error("raw writes another mechanism than m")
    c = _testhelper_stored_steps(raw, m.catalytic)
    # The test reader's steps run group by group, each group in edge order.
    flat = EnzymeRates.Derivation._flat_edges(m.catalytic)
    out = Pair{Symbol, Any}[]
    for i in eachindex(c.edge_of)
        s = g.steps[findfirst(==(c.edge_of[i]), flat)]
        forward = !c.flipped[i]
        tag = allo ? g.tags[s.group] : :None
        states = tag === :NonequalAI ? ((s.names, ""), (s.inactive_names, "_T")) :
                 ((s.names, ""),)
        # A step the inactive law does not keep has no inactive constants.
        for (names, suffix) in states
            isempty(names) && suffix == "_T" && continue
            if s.rapid_equilibrium
                K = all_params[only(names)]
                push!(out, Symbol("K", i, suffix) =>
                    _testhelper_is_binding(s) || forward ? K : 1 / K)
            else
                kf_key, kr_key = names
                kf = all_params[kf_key]
                kr = all_params[kr_key]
                f, r = forward ? (kf, kr) : (kr, kf)
                push!(out, Symbol("k", i, "f", suffix) => f,
                      Symbol("k", i, "r", suffix) => r)
            end
        end
    end
    if allo
        for (pos, site) in enumerate(m.allosteric.regulatory_sites)
            j = _testhelper_positional_site_idx(site, pos, raw.regulatory_sites)
            for (lig, tag) in zip(site.ligands,
                                  site.allosteric_states)
                for (inactive, suffix) in ((false, ""), (true, "_T"))
                    tag === :EqualAI && inactive && continue
                    K = EnzymeRates.Derivation._regulator_constant(m, lig, tag, inactive)
                    K === nothing ||
                        push!(out, Symbol("K_", lig, suffix, "_reg", j) => all_params[K])
                end
            end
        end
    end
    merge((; out...), all_params)
end

"""Check if mechanism has any rapid-equilibrium steps."""
_testhelper_has_re_steps(m) = any(s -> s.rapid_equilibrium, _testhelper_graph(m).steps)

"""
Convert structural-named params to ODE params (large k_if/k_ir for all steps).
Steps sharing a kinetic_group share the same K (or k_f/k_r) param.

For binding RE steps (metabolite on LHS, canonical form):
    K = Kd = kr/kf, so k_if = 1e6, k_ir = 1e6 * K.
For RE isomerization steps (no metabolite, enzyme-only):
    K = Ka = kf/kr, so k_if = 1e6 * K, k_ir = 1e6.
"""
function _testhelper_raw_to_ode_params(m, raw_params)
    flat = _testhelper_graph(m).steps
    # A canonical RE binding step has a metabolite on LHS (canonical form
    # invariant: all RE binding steps are written `E + S ⇌ ES`).
    is_binding_step = Bool[s.rapid_equilibrium && !isempty(s.consumed) for s in flat]
    _testhelper_lookup(k) = Float64(raw_params[k])
    param_keys = Symbol[]
    param_vals = Float64[]
    for (i, step) in enumerate(flat)
        push!(param_keys, Symbol("k$(i)f"))
        push!(param_keys, Symbol("k$(i)r"))
        if step.rapid_equilibrium
            # Look up the step's structural K key, which its kinetic group shares
            K = _testhelper_lookup(only(step.names))
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
            fwd_key, rev_key = step.names
            push!(param_vals, _testhelper_lookup(fwd_key))
            push!(param_vals, _testhelper_lookup(rev_key))
        end
    end
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
    m::EnzymeRates.Mechanism,
    params, concs,
)
    step_data = []
    for (step_idx, step) in enumerate(_testhelper_graph(m).steps)
        i, j = step.from, step.to

        m_lhs = step.consumed
        m_rhs = step.released

        kf = Float64(params[Symbol("k$(step_idx)f")])
        kr = Float64(params[Symbol("k$(step_idx)r")])

        rf = isempty(m_lhs) ? kf : kf * concs[m_lhs[1]]
        rr = isempty(m_rhs) ? kr : kr * concs[m_rhs[1]]

        push!(step_data, (i, j, rf, rr))
    end

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
    m::EnzymeRates.Mechanism,
    params, concs,
)
    g = _testhelper_graph(m)
    n = length(g.forms)
    ref_name, nu_ref = _testhelper_reference_metabolite(m)

    # One unit of enzyme, all free at the start, so the flux is the rate per unit enzyme.
    u0 = zeros(n)
    u0[1] = 1.0

    rhs! = _testhelper_build_ode_rhs(m, params, concs)
    prob = ODEProblem(rhs!, u0, (0.0, 1e6))
    # 1e-13 (vs 1e-12) keeps the stiff multi-step solves converged to the
    # exact rate equation regardless of catalytic-step storage order; at 1e-12
    # a ter-bi corner drifts to ~1e-6, exceeding the cross-check rtol.
    sol = solve(prob, RadauIIA9(); abstol=1e-13, reltol=1e-13)
    u_ss = sol.u[end]

    v = 0.0
    for (step_idx, step) in enumerate(g.steps)
        i, j = step.from, step.to

        m_lhs = step.consumed
        m_rhs = step.released

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

# ── Modular test functions for _testhelper_MechanismTestSpec ────────────────

function _testhelper_test_structure(spec::_testhelper_MechanismTestSpec)
    m = spec.mechanism
    @testset "Structure" begin
        g = _testhelper_graph(m)
        @test length(g.forms) == spec.expected_n_states
        @test length(g.steps) == spec.expected_n_steps
        @test metabolites(m) == spec.metabolite_names
        # Structural parameter names must be injective on the load-bearing
        # paths the fitter consumes (parameters + the fitted parameters). A collision
        # would silently shorten the destructured params NamedTuple.
        @test allunique(parameters(m))
    end
end

"""Whether the expression `ex` mentions the symbol `s`."""
_testhelper_mentions(ex, s::Symbol) =
    ex === s || (ex isa Expr && any(a -> _testhelper_mentions(a, s), ex.args))

"""Classify a dep expression as Haldane (RHS references Keq), Mirror
(RHS is a single Symbol), or Wegscheider (RHS Expr without Keq)."""
function _testhelper_classify_dep_expr(expr)
    if expr isa Symbol
        return :mirror
    elseif _testhelper_mentions(expr, :Keq)
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
        dep_exprs, indep = _testhelper_dependent_param_exprs(m)
        # A dependent (Haldane/Wegscheider-derived) parameter must never appear in
        # the independent set returned as the fitted parameters. The LDH i-state specs
        # have a shared :EqualAI catalytic reverse rate that is Haldane-dependent in the
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
        rate = compile_rate_equation(m)
        @test all(1:n_trials) do _
            params, Keq, concs, all_params =
                _testhelper_random_independent_params_concs(
                    m, met_names; rng=rng)
            # _testhelper_reference_qssa reads stored-step k$(i)f/k$(i)r keys
            pos_params = _testhelper_stored_params(m, all_params)
            isapprox(
                rate(concs, params, Keq),
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
        rate = compile_rate_equation(m)
        @test all(1:n_trials) do _
            params, Keq, concs, all_params =
                _testhelper_random_independent_params_concs(
                    m, met_names; rng=rng)
            # The textbook law multiplies by the total enzyme, here one unit of it.
            p = merge(_testhelper_written_params(m, all_params; raw=spec.raw),
                      (Keq=Keq, Et=1.0))
            isapprox(
                rate(concs, params, Keq),
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
        params, Keq, _, _ = _testhelper_random_independent_params_concs(
            m, met_names; rng=rng)
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
        v_eq = compile_rate_equation(m)(eq_concs, params, Keq)
        @test abs(v_eq) < 1e-10
    end
end

"""The number of terms in the denominator of `m`'s derived rate law after the solve is
substituted: an MWC law counts both conformations, an inactive conformation that is free E
alone counting 1."""
function _testhelper_denominator_terms(m)
    L = EnzymeRates.Derivation._memo_law(m)
    L isa EnzymeRates.Derivation._Law && return EnzymeRates.Derivation._nterms(L.D)
    EnzymeRates.Derivation._nterms(L.active.D) +
        EnzymeRates.Derivation._nterms(L.inactive.D)
end

# Small laws take 1-8 ns on an aarch64 development machine, and 5-38 ns on CI under coverage
# counters, varying by runner, with single readings up to 56 ns. Their intercept is
# therefore wide and catches only a large slowdown; the 0.5 ns per term binds on large laws.
# A job that misses the bound is rerun rather than the bound loosened.
_testhelper_rate_equation_time_limit(n_terms) = 50e-9 + 0.5e-9 * n_terms

# The times come from the speed child, which runs without CI's coverage counters and bounds
# checks (`_testhelper_gate_readings`). The allocations are counted in this process, so CI
# checks them under both.
function _testhelper_test_performance(spec::_testhelper_MechanismTestSpec)
    m = spec.mechanism
    met_names = spec.metabolite_names
    @testset "Performance" begin
        rng = Random.MersenneTwister(42)
        params, Keq, concs, _ = _testhelper_random_independent_params_concs(
            m, met_names; rng=rng)
        allocs, _ = _testhelper_test_rate_equation_performance(
            compile_rate_equation(m), params, Keq, concs)
        @test allocs == 0
        n = _testhelper_denominator_terms(m)
        (; law, per_row, allowed) = _testhelper_gate_readings(m)
        @test law < _testhelper_rate_equation_time_limit(n)
        @test _testhelper_loss_allocations(_testhelper_loss_problem(m)) == 0
        @test _testhelper_loss_allocations(
            _testhelper_loss_problem(m; absolute = true)) == 0
        @info "$(spec.name): rate law $(round(1e9law; digits = 1)) ns, $n terms, " *
              "loss! $(round(1e9per_row; digits = 1)) ns per row, " *
              "allowed $(round(1e9allowed; digits = 1))"
        @test per_row <= allowed
        rate = compile_rate_equation(m)
        @inferred metabolites(rate)
        @inferred parameters(rate)
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
        rate = compile_rate_equation(m)
        @test all(1:n_trials) do _
            params, Keq, concs, all_params =
                _testhelper_random_independent_params_concs(
                    m, met_names; rng=rng)
            # Convert structural params to positional k_if/k_ir for ODE
            pos_params = _testhelper_stored_params(m, all_params)
            ode_params = has_re ?
                _testhelper_raw_to_ode_params(m, all_params) :
                pos_params
            v_ode = _testhelper_ode_steady_state_flux(m, ode_params, concs)
            v_ka = rate(concs, params, Keq)
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
        @test startswith(last(split(s, "\n")), "v = ")  # Proper format
        @test occursin(") / (", s)      # Has denominator
        @test occursin("= params", s)   # Has params destructuring
        @test occursin("= concs", s)    # Has concs destructuring
        @test isempty(_testhelper_print_artifacts(s))

        # Numerical equivalence test
        rng = Random.MersenneTwister(9000 + hash(spec.name) % 1000)
        rate = compile_rate_equation(m)
        @test all(1:10) do _
            params, Keq, concs, _ =
                _testhelper_random_independent_params_concs(
                    m, met_names; rng=rng)
            isapprox(
                rate(concs, params, Keq),
                _testhelper_eval_rate_string(s, params, concs, Keq);
                rtol=1e-10)
        end

        # The eleven pinned specs turn S into P, so H reads S - P / Keq.
        if spec.expected_factored_num !== nothing ||
           spec.expected_factored_denom !== nothing
            q, h, d = _testhelper_extract_factors(last(split(s, "\n")))
            @test h == "S - P / Keq"
            spec.expected_factored_num === nothing || @test q == spec.expected_factored_num
            spec.expected_factored_denom === nothing ||
                @test d == spec.expected_factored_denom
        end
    end
end

"""
The three factors of a printed `v` line, `v = (Q) * (H) / (D)`, as written
inside their parentheses: `Q` (a bare name as it is, `"1"` when the line leaves it out),
`H` and `D`. Each factor is a parenthesized group or a bare name, so the split runs at the
` * (` and ` / (` outside every parenthesis.
"""
function _testhelper_extract_factors(v_line::AbstractString)
    s = String(chopprefix(v_line, "v = "))
    parts, depth, from, i = String[], 0, firstindex(s), firstindex(s)
    while i <= lastindex(s)
        depth += (s[i] == '(') - (s[i] == ')')
        if depth == 0 && (startswith(SubString(s, i), " * (") ||
                          startswith(SubString(s, i), " / ("))
            push!(parts, s[from:prevind(s, i)])
            i += 3
            from = i
            continue
        end
        i = nextind(s, i)
    end
    push!(parts, s[from:end])
    inner(t) = startswith(t, "(") && endswith(t, ")") ? t[2:prevind(t, lastindex(t))] : t
    length(parts) == 2 ? ("1", inner(parts[1]), inner(parts[2])) :
                         (inner(parts[1]), inner(parts[2]), inner(parts[3]))
end

function _testhelper_test_rescaling(spec::_testhelper_MechanismTestSpec; seed=100)
    m = spec.mechanism
    @testset "Rescaling to a measured turnover" begin
        rng = Random.MersenneTwister(seed)
        params, Keq = _testhelper_random_independent_params(m; rng)
        met_names = metabolites(m)
        concs = NamedTuple{Tuple(met_names)}(Tuple(0.5 + rand(rng) for _ in met_names))
        rate = compile_rate_equation(m)
        v_orig = rate(concs, params, Keq)
        @test isfinite(v_orig) && v_orig != 0

        # Scaling every steady-state rate constant by α scales every rate by α.
        ss_names = _testhelper_steady_state_names(m)
        α = 0.1 + 9.9 * rand(rng)
        scaled_params = NamedTuple{keys(params)}(Tuple(
            k in ss_names ? v * α : v
            for (k, v) in zip(keys(params), values(params))
        ))
        @test rate(concs, scaled_params, Keq) ≈ α * v_orig rtol=1e-9

        # An anchor at those concentrations sets the turnover there and leaves every
        # constant but the steady-state rate constants alone.
        target = sign(v_orig) * (0.1 + 9.9 * rand(rng))
        anchored = EnzymeRates.Fitting._scale_rate_constants(m, params, target / v_orig)
        for k in keys(params)
            if !(k in ss_names)
                @test anchored[k] == params[k]
            end
        end
        @test rate(concs, anchored, Keq) ≈ target rtol=1e-9

        # Every rate changes by the same factor.
        other = NamedTuple{Tuple(met_names)}(Tuple(0.5 + rand(rng) for _ in met_names))
        @test rate(other, anchored, Keq) / rate(other, params, Keq) ≈ target / v_orig rtol =
              1e-8
    end
end

"""
Assert the compiled rate law stays finite when any single metabolite concentration
is zero (real kinetic data routinely has zeros). Zeroing a substrate or product
must still leave a nonzero net rate — the opposite-direction metabolites keep
driving flux — so a `0.0` there normally signals a `1/conc` term that blew the
denominator to `Inf`. A regulator zero (essential activator) can legitimately
zero the rate, so only finiteness is required there. Zeroing a metabolite that
sits in the free-enzyme weight `D[g_free]` also keeps a finite rate with a
nonzero reverse flux; the mass-action gates in `test_rate_equation_derivation.jl`
check that case against the ground truth.
"""
function _testhelper_test_zero_metabolite_finite(spec::_testhelper_MechanismTestSpec)
    m = spec.mechanism
    @testset "Zero-metabolite finiteness" begin
        rng = Random.MersenneTwister(777 + hash(spec.name) % 1000)
        mets = collect(metabolites(m))
        sub_prod = Set{Symbol}(_testhelper_substrates(m))
        union!(sub_prod, _testhelper_products(m))
        params, Keq = _testhelper_random_independent_params(m; rng)
        rate = compile_rate_equation(m)
        for zeroed in mets
            cvals = Tuple(n == zeroed ? 0.0 : 0.5 + rand(rng) for n in mets)
            concs = NamedTuple{Tuple(mets)}(cvals)
            v = rate(concs, params, Keq)
            @test isfinite(v)
            zeroed in sub_prod && @test v != 0.0
        end
    end
end

# The plain specs whose compiled body also reads BigFloat inputs in the mass-action check,
# beside every allosteric spec: an isomerization, a ping-pong, a ter-bi and a law with
# both a Haldane and a Wegscheider dependent.
const _testhelper_bigfloat_body_specs = ("Segel Iso Uni Uni", "Segel Ping Pong Bi Bi",
                                         "Segel Ordered Ter Bi", "Random-order Bi-Bi")

"""The spec's rate law against the exact steady state of its mass-action network
(`_testhelper_check_law_against_oracle`); an allosteric spec runs with one catalytic
subunit, where the MWC law is the exact steady state. The allosteric specs and
`_testhelper_bigfloat_body_specs` also run 256-bit BigFloat inputs through the compiled
body, to 1e-20."""
function _testhelper_test_mass_action_oracle(spec::_testhelper_MechanismTestSpec)
    @testset "Mass-action oracle" begin
        bigfloat = _testhelper_is_allosteric(spec.mechanism) ||
                   spec.name in _testhelper_bigfloat_body_specs
        _testhelper_check_law_against_oracle(spec.mechanism; bigfloat)
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
        _testhelper_test_rescaling(spec)
    end
end

"""The allocations of `parameters(rate)` of a compiled rate law `rate`, once compiled,
measured inside a function."""
function _testhelper_parameters_allocations(rate)
    parameters(rate)
    @allocated parameters(rate)
end

@testset "Keq is an argument" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    names = parameters(m)
    @test :Keq ∉ names && :E_total ∉ names
    params = NamedTuple{names}(Tuple(rand(MersenneTwister(1), length(names)) .+ 0.5))
    concs = (S = 0.7, P = 0.3)
    rate = compile_rate_equation(m)
    v = rate(concs, params, 2.0)
    # The law is the textbook rapid-equilibrium uni-uni law per unit enzyme, its reverse
    # rate constant set by the Haldane relation
    # Keq = k_ES_to_EP * K_EP_to_E_P / (k_EP_to_ES * K_ES_to_E_S).
    (; K_ES_to_E_S, K_EP_to_E_P, k_ES_to_EP) = params
    k_EP_to_ES = k_ES_to_EP * K_EP_to_E_P / (K_ES_to_E_S * 2.0)
    @test v ≈ (k_ES_to_EP * concs.S / K_ES_to_E_S - k_EP_to_ES * concs.P / K_EP_to_E_P) /
              (1 + concs.S / K_ES_to_E_S + concs.P / K_EP_to_E_P) rtol = 1e-12
    # Extra names in params are ignored.
    @test rate(concs, merge(params, (Keq = 99.0, E_total = 99.0)), 2.0) == v
    @test rate((S = 0.3 / 2.0, P = 0.3), params, 2.0) ≈ 0 atol = 1e-12
    dep = dependent_parameter_values(m, params, 2.0)
    @test dep isa NamedTuple && isdisjoint(keys(dep), names)
    @test dep.k_EP_to_ES ≈ k_EP_to_ES rtol = 1e-12
    @test !occursin("E_total", rate_equation_string(m))
    @test occursin("Keq", rate_equation_string(m))
end

@testset "compile_rate_equation" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    names = parameters(m)
    params = NamedTuple{names}(Tuple(rand(MersenneTwister(1), length(names)) .+ 0.5))
    concs = (S = 0.7, P = 0.3)
    rate = compile_rate_equation(m)
    @test parameters(rate) == parameters(m)
    @test metabolites(rate) == Tuple(metabolites(m))
    # Measured behind a function: before Julia 1.12, `@allocated` in this body counts the
    # 48 bytes of dispatching the call on `rate`, whose type the body does not know.
    allocs, _ = _testhelper_test_rate_equation_performance(rate, params, 2.0, concs)
    @test allocs == 0
    @test compile_rate_equation(m) === rate
    # It prints one line, which names what it is, not its type parameter.
    @test sprint(show, rate) == "compiled rate equation"
    @test repr(MIME"text/plain"(), rate) == "compiled rate equation"

    # An allosteric mechanism compiles alike; its law equals the mass-action steady state
    # of one protomer (`_testhelper_check_law_against_oracle`).
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES     :: EqualAI
            ES <--> EP   :: OnlyA
            EP ⇌ E + P     :: EqualAI
        end
        regulatory_site(2): X::OnlyA
    end
    arate = compile_rate_equation(am)
    @test parameters(arate) == parameters(am)
    @test metabolites(arate) == Tuple(metabolites(am)) == (:S, :P, :X)
    @test compile_rate_equation(am) === arate && arate !== rate
    anames = parameters(am)
    aparams = NamedTuple{anames}(Tuple(rand(MersenneTwister(2), length(anames)) .+ 0.5))
    aconcs = (S = 0.7, P = 0.3, X = 0.2)
    @test isfinite(arate(aconcs, aparams, 2.0))
    aallocs, _ = _testhelper_test_rate_equation_performance(arate, aparams, 2.0, aconcs)
    @test aallocs == 0
    _testhelper_check_law_against_oracle(am)

    # A fit evaluates the law `compile_rate_equation` gives its mechanism.
    table = (group = ["g", "g"], Rate = [1.0, 2.0], S = [1.0, 2.0], P = [0.1, 0.1])
    @test FittingProblem(m, table, 2.0).rate === rate
    @test FittingProblem(am, (; table..., X = [0.5, 0.5]), 2.0).rate === arate
end

@testset "parameters agree across representations" begin
    # A mechanism and its compiled rate law, whose names are compile-time constants.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        rate = compile_rate_equation(m)
        @test parameters(m) == parameters(rate)
        params, Keq = _testhelper_random_independent_params(m)
        dependents = keys(dependent_parameter_values(m, params, Keq))
        @test issorted(dependents) && isdisjoint(dependents, parameters(m))
        @test _testhelper_parameters_allocations(rate) == 0
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

# A fused and a plain binding of B share one rapid-equilibrium kinetic group.
const _testhelper_fused_and_plain_binding =
    @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
        E + A <--> EA; (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)
        EQ + P <--> EPQ; E + Q ⇌ EQ end end)

# Mechanisms with fused steps, shared by the mass-action and orientation testsets.
const _testhelper_fused_cases = [
        # uni-uni: fused release / fused binding, both orientations, both flags
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> ES; ES <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S ⇌ ES; ES <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> ES; ES ⇌ E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> ES; E + P <--> ES end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> EP; EP <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            EP <--> E + S; E + P <--> EP end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            EP ⇌ E + S; EP <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            EP <--> E + S; EP ⇌ E + P end end),
        # ordered bi-bi, one central complex (Segel IX-87) and mixes
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> EA; EA + B <--> EAB; EAB <--> EQ + P
            EQ <--> E + Q end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> EA; EA + B <--> EAB; EAB ⇌ EQ + P
            EQ <--> E + Q end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A ⇌ EA; EPQ <--> EA + B; EPQ <--> EQ + P
            EQ ⇌ E + Q end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A ⇌ EA; EPQ <--> EA + B; EPQ ⇌ EQ + P
            EQ ⇌ E + Q end end),
        # ping-pong with an RE fused release
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> EA; EA ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP
            EB_res_plusA_minusP <--> E + Q end end),
        # merged random bi-bi with RE bindings and releases, and a reversed fused
        # binding inside an all-SS merged random bi-bi
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A ⇌ EA; E + B ⇌ EB; EA + B ⇌ EAB; EB + A ⇌ EAB
            EAB <--> EQ + P; EAB <--> EP + Q; EQ ⇌ E + Q; EP ⇌ E + P
        end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> EA; E + B <--> EB; EA + B <--> EPQ
            EPQ <--> EB + A; EPQ <--> EQ + P; EPQ <--> EP + Q
            EQ <--> E + Q; EP <--> E + P end end),
        # dead ends on a merged complex
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; modifiers: I
            steps: begin
            E + A ⇌ EA; EA + B ⇌ EAB; EAB <--> EQ + P; EQ ⇌ E + Q
            E + I ⇌ EI; EQ + I ⇌ EIQ end end),
        # ping-pong: parallel product-release route via a mixed substrate/product complex
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; modifiers: I
            steps: begin
            E + A ⇌ EA; EI + A ⇌ EAI; EQ + A ⇌ EAQ
            E + I ⇌ EI; EA + I ⇌ EAI
            EP_res_plusA_minusP + I ⇌ EIP_res_plusA_minusP
            E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP
            E + Q <--> EQ; EA + Q <--> EAQ
            EA <--> EP_res_plusA_minusP
            EAI <--> EIP_res_plusA_minusP
            EB_res_plusA_minusP + P ⇌ EBP_res_plusA_minusP
            EI_res_plusA_minusP + P ⇌ EIP_res_plusA_minusP
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            EP_res_plusA_minusP + B ⇌ EBP_res_plusA_minusP
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP end end),
        # both substrates bind in one rapid-equilibrium step
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A + B ⇌ EAB; EAB <--> EPQ; EPQ <--> EQ + P
            EQ <--> E + Q end end),
        # the first substrate also binds as a dead-end modifier copy, at rapid
        # equilibrium on E and at steady state on EQ; both the derivation and
        # the oracle count the first substrate by name
        @enzyme_mechanism(begin substrates: A, B; products: P, Q
            steps: begin
            E + A ⇌ EA; EA + B <--> EAB; EAB <--> EQ + P; EQ <--> E + Q
            E + A ⇌ EAmod; EQ + A <--> EAmodQ end end),
        _testhelper_fused_and_plain_binding,
]

@testset "each rapid-equilibrium segment is rooted at its first form" begin
    # The engine weighs each form of a rapid-equilibrium segment against the segment's
    # first form in the canonical form order (`rapid_equilibrium_segments`), so free E,
    # form 1, roots its own segment. EAb comes before EAmod, where A binds without
    # reacting, there, by name. The law in lowest terms does not depend on the root.
    m = @enzyme_mechanism begin
        substrates: A, Ab
        products: P
        steps: begin
            E + A ⇌ EA
            EA + Ab ⇌ EAAb
            EAAb <--> EP
            EP ⇌ E + P
            E + A ⇌ EAmod
            E + Ab ⇌ EAb
            EAb + A ⇌ EAAb
            EA <--> E_res_plusA_minusP + P
        end
    end
    form(name) = only(findall(==(name), m.catalytic.enzyme_forms))
    a_copy, ab, e = form(:EAmod), form(:EAb), form(:E)
    @test e == 1 && ab < a_copy
    s = ER.Mechanisms.rapid_equilibrium_segments(m.catalytic,
                                                 trues(length(m.catalytic.steps)))
    @test all(seg -> seg[1] == minimum(seg), s.segments)
    @test s.segments[s.segment_of[e]][1] == e
    @test s.segment_of[a_copy] == s.segment_of[ab] == s.segment_of[e]
    covalent = form(:E_res_plusA_minusP)
    @test s.segments[s.segment_of[covalent]] == [covalent]
    _testhelper_check_law_against_oracle(m)
end

@testset "the free-enzyme weight is free E's weight in the reduced denominator" begin
    # Ping-pong with a rapid-equilibrium second chemistry step: every form lies in one
    # rapid-equilibrium segment rooted at free E, where E_res_plusA_minusP sits at
    # Q·K/B relative to E. Reducing the polynomials to lowest terms in the
    # concentrations multiplies them by B, so free E weighs B in the reduced
    # denominator; the MWC law reads this weight as free E's (`_cross_weights`).
    m = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            E + Q ⇌ EQ
        end
    end
    s = EnzymeRates.Derivation._solve(m)
    L = EnzymeRates.Derivation._law(m, s)
    @test L.trees[1] * L.shift == EnzymeRates.Derivation._term(
        EnzymeRates.Derivation._monomial(L.vars.index[:B] => 1))
    # Free E's weight over the denominator is its steady-state share of the enzyme.
    point, concs, params = _testhelper_law_inputs(m, s.fitted, Random.MersenneTwister(8))
    env = merge(params, concs)
    share = setprecision(BigFloat, 256) do
        _testhelper_poly_value(L.trees[1] * L.shift, L.vars, env) /
            _testhelper_poly_value(L.D, L.vars, env)
    end
    @test abs(share - _testhelper_steady_state_rate(m, point, concs; occupancy = 1)) <=
          big(10)^-20 * share
end

@testset "fused steps derive the mass-action rate" begin
    for m in _testhelper_fused_cases
        _testhelper_check_law_against_oracle(m)
    end
    # The B group holds the fused EA + B → EPQ and the dead-end EQ + B ⇌ EBQ.
    # Fitted: A (2) + B (1) + P (2) + Q (1) − 1 Haldane = 5.
    @test length(parameters(_testhelper_fused_and_plain_binding)) == 5
end

# Theorell–Chance and two-metabolite mechanisms, shared by the mass-action and
# orientation testsets.
const _testhelper_tc_cases = (
    tc_ss = @enzyme_mechanism(begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EQ + P
            EQ <--> E + Q
        end
    end),
    tc_re_outer = @enzyme_mechanism(begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EQ + P
            EQ ⇌ E + Q
        end
    end),
    tc_re_step = @enzyme_mechanism(begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B ⇌ EQ + P
            EQ <--> E + Q
        end
    end),
    two_in = @enzyme_mechanism(begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A + B <--> EAB
            EAB <--> EPQ
            EPQ <--> E + P + Q
        end
    end),
)

@testset "Theorell–Chance and two-metabolite steps" begin
    (; tc_ss, tc_re_outer, tc_re_step, two_in) = _testhelper_tc_cases
    @test length(parameters(tc_ss)) == 5
    @test length(parameters(tc_re_outer)) == 3
    for m in (tc_ss, tc_re_outer, tc_re_step, two_in)
        _testhelper_check_law_against_oracle(m)
    end
    # canonical orientation: the substrate side is `from`
    e = only(e for e in tc_ss.catalytic.steps
             if !isempty(e.taken_up) && !isempty(e.given_off))
    @test tc_ss.catalytic.enzyme_forms[e.from] === :EA &&
          tc_ss.catalytic.metabolites[e.taken_up] == [:B]
    @test :k_EA_B_to_EQ_P in _testhelper_step_constant_names(tc_ss)
    @test :k_EQ_P_to_EA_B in _testhelper_step_constant_names(tc_ss)
end

@testset "reversing written steps changes nothing" begin
    rng = Random.MersenneTwister(7)
    # The written form `raw` (`_testhelper_written`) with every step reversed, a random
    # subset, and every group's steps in reverse order.
    function _testhelper_variants(raw)
        group = [e.kinetic_group for e in raw.steps]
        within = reduce(vcat, [reverse(findall(==(l), group)) for l in unique(group)])
        ((; raw..., steps = _testhelper_turned.(raw.steps)),
         (; raw..., steps = [rand(rng, Bool) ? _testhelper_turned(e) : e
                             for e in raw.steps]),
         (; raw..., steps = raw.steps[within]))
    end
    # The Theorell–Chance step of `tc_ss` written backwards.
    tc_backwards = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EQ + P <--> EA + B
            EQ <--> E + Q
        end
    end
    @test tc_backwards == _testhelper_tc_cases.tc_ss
    @test ER.rate_equation_string(tc_backwards) ==
          ER.rate_equation_string(_testhelper_tc_cases.tc_ss)
    mechanisms = Any[spec.mechanism for spec in MECHANISM_TEST_SPECS
              if !_testhelper_is_allosteric(spec.mechanism)]
    append!(mechanisms, _testhelper_fused_cases, values(_testhelper_tc_cases),
            [tc_backwards])
    for m in mechanisms
        for flipped in _testhelper_variants(_testhelper_written(m))
            m2 = _testhelper_build(flipped)
            @test m2 == m
            @test ER.rate_equation_string(m2) ==
                  ER.rate_equation_string(m)
        end
    end
    for spec in MECHANISM_TEST_SPECS
        _testhelper_is_allosteric(spec.mechanism) || continue
        am = spec.mechanism
        for flipped in _testhelper_variants(_testhelper_written(am))
            am2 = _testhelper_build(flipped)
            @test am2 == am
            @test ER.rate_equation_string(am2) ==
                  ER.rate_equation_string(spec.mechanism)
        end
    end
end

# ── Regression: every inactive-state parameter must be defined ──────────────
# Allosteric mechanisms from an LDH `identify_rate_equation` run crashed with
# `UndefVarError` on undefined I-state parameters: the I-state polynomials
# referenced names the dep-assignment and destructuring machinery never emitted.
# The undefined constants were an I-state lactate dissociation constant and an
# I-state reverse chemistry rate in dead inactive states, and an I-state NAD
# binding rate in a live one. Those mechanisms are not kept here. The regression
# runs the three `LDH i-state …` specs of MECHANISM_TEST_SPECS (defined in
# testhelpers.jl), whose Haldane-dependent
# reverse rate must stay out of the fitted set, and checks that each rate
# equation defines every name it references.
const _LDH_ISTATE_MECHS = [spec.mechanism for spec in MECHANISM_TEST_SPECS
                           if startswith(spec.name, "LDH i-state")]

# Parameter Symbols referenced on an assignment/`v` RHS but never defined
# (destructured from `params`/`concs`, assigned as an LHS, or `Keq`, an input of its own
# that is bound before the law runs). Empty ⟺ the rendered rate equation is closed:
# every referenced name has a definition.
function _testhelper_undefined_rhs_symbols(s::AbstractString)
    # A name may hold superscripts (a coefficient's); a called function (`sqrt`) is no name.
    ident = r"[A-Za-z_][A-Za-z0-9_⁰¹²³⁴⁵⁶⁷⁸⁹]*(?![A-Za-z0-9_⁰¹²³⁴⁵⁶⁷⁸⁹(])"
    defined = Set{Symbol}([:Keq])
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
    for m in _LDH_ISTATE_MECHS
        pnames = parameters(m)
        mets = Tuple(EnzymeRates.metabolites(m))
        prods = Set(_testhelper_products(m))
        params = NamedTuple{pnames}(ntuple(_ -> 1.3, length(pnames)))
        # No UndefVarError: the @generated body compiles and evaluates finite at
        # products = 0 (initial-rate conditions).
        concs0 = NamedTuple{mets}(
            ntuple(i -> mets[i] in prods ? 0.0 : 1.5, length(mets)))
        @test isfinite(compile_rate_equation(m)(concs0, params, 20000.0))
        # DEFINED ⊇ REFERENCED on the rendered transcript.
        @test isempty(
            _testhelper_undefined_rhs_symbols(EnzymeRates.rate_equation_string(m)))
    end
end

@testset "allosteric dependent-param graph is sound" begin
    for m in ALLOSTERIC_UNDEFVAR_REPRODUCERS
        @test _testhelper_dep_graph_is_sound(_testhelper_dependent_param_exprs(m)...)
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
    # in `_solve` ties every cross-state affinity split
    # directly, so all three reproducers (D1's `:EqualAI`-shared `koff` merge, D2's
    # `:NonequalAI` split, D3's steady-state speed) are detailed-balance-correct.
    for m in ALLOSTERIC_UNDEFVAR_REPRODUCERS
        pn = collect(parameters(m))
        subs = _testhelper_substrates(m)
        prods = _testhelper_products(m)
        mets = collect(EnzymeRates.metabolites(m))
        Keq = 20000.0
        p_each = Keq^(1 / length(prods))
        rng = Random.MersenneTwister(1)
        maxv = 0.0
        rate = compile_rate_equation(m)
        for _ in 1:100
            pv = Dict(p => exp(2 * randn(rng)) for p in pn)
            params = NamedTuple{Tuple(pn)}(Tuple(pv[p] for p in pn))
            cv = Dict{Symbol,Float64}()
            for s in subs; cv[s] = 1.0; end
            for p in prods; cv[p] = p_each; end
            for s in mets; haskey(cv, s) || (cv[s] = 2.5); end
            concs = NamedTuple{Tuple(mets)}(Tuple(cv[s] for s in mets))
            maxv = max(maxv, abs(rate(concs, params, Keq)))
        end
        @test maxv < 1e-8
    end
end

# ── Standalone tests ───────────────────────────────────────────────────────────

"""The body of the generated call of `m`'s compiled rate law, as the generator builds
it."""
_testhelper_rate_body(m) = ER.Derivation._rate_body(
    ER.Derivation._derivation_input(typeof(compile_rate_equation(m))))

"""Each local the generated rate-law body of `m` binds, in order, with its right-hand
side."""
_testhelper_body_locals(m) =
    Pair{Symbol, Any}[ex.args[1] => ex.args[2] for ex in _testhelper_rate_body(m).args
                      if ex isa Expr && ex.head === :(=) && ex.args[1] isa Symbol]

"""The expression the generated rate-law body of `m` returns."""
_testhelper_generated_rate_return(m) =
    only(ex.args[1] for ex in _testhelper_rate_body(m).args
         if ex isa Expr && ex.head === :return)

@testset "the rate law's polynomial body uses 2-arg +/* calls" begin
    # The fitter calls the rate law millions of times per CV fold. The
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
    # Every local the body binds, and what it returns, for a plain law and an MWC law.
    mwc = only(s for s in MECHANISM_TEST_SPECS
               if s.name == "MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites")
    for m in (spec.mechanism, mwc.mechanism)
        foreach(_testhelper_walk!, last.(_testhelper_body_locals(m)))
        _testhelper_walk!(_testhelper_generated_rate_return(m))
    end
    @test isempty(bad)
end

@testset "the rate law's time per call grows at most in proportion to its size" begin
    @test _testhelper_denominator_terms.(_testhelper_largest_enumerated_mechanisms) ==
          [396, 396, 396]
    @test _testhelper_denominator_terms.(_testhelper_large_ter_bi_mechanisms) ==
          [589, 799, 858]
    for (i, m) in enumerate(vcat(_testhelper_largest_enumerated_mechanisms,
                                 _testhelper_large_ter_bi_mechanisms))
        params, Keq, concs, _ = _testhelper_random_independent_params_concs(
            m, collect(metabolites(m)); rng = Random.MersenneTwister(42))
        allocs, _ = _testhelper_test_rate_equation_performance(
            compile_rate_equation(m), params, Keq, concs)
        n = _testhelper_denominator_terms(m)
        (; law, per_row, allowed) = _testhelper_gate_readings(m)
        @test allocs == 0
        @test law < _testhelper_rate_equation_time_limit(n)
        @test _testhelper_loss_allocations(_testhelper_loss_problem(m)) == 0
        @test _testhelper_loss_allocations(
            _testhelper_loss_problem(m; absolute = true)) == 0
        @info "large law $i: rate law $(round(1e9law; digits = 1)) ns, $n terms, " *
              "loss! $(round(1e9per_row; digits = 1)) ns per row, " *
              "allowed $(round(1e9allowed; digits = 1))"
        @test per_row <= allowed
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

@testset "the generated rate-law body is in Horner form" begin
    # Horner form multiplies by each factor once per branch instead of once per term,
    # which keeps the 396-term gate laws at tens of ns per call. Greedy Horner takes a
    # factor out of every summand it divides, so no local divides two summands of one sum,
    # in any local's definition or in what the body returns.
    for m in (_testhelper_largest_enumerated_mechanisms[2],
              only(s for s in MECHANISM_TEST_SPECS
                   if s.name == "MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites"
                  ).mechanism)
        locals = _testhelper_body_locals(m)
        names = first.(locals)
        for ex in [last.(locals); _testhelper_generated_rate_return(m)]
            @test isempty(_testhelper_unfactored_sums(ex, names))
        end
    end
end

@testset "the rate law reads concs and params by name in any order" begin
    # The body reads each input by name (`getfield`), so reordered NamedTuples, a data row
    # with Rate and group fields, and `params` merged with
    # `dependent_parameter_values(m, params, Keq)` give the same value without allocating.
    # A metabolite named x1 reads as a concentration beside the body's positional local x1.
    m = only(s for s in MECHANISM_TEST_SPECS if s.name == "Random-order Bi-Bi").mechanism
    params, Keq, concs, all_params = _testhelper_random_independent_params_concs(
        m, collect(metabolites(m)); rng = Random.MersenneTwister(91))
    rate = compile_rate_equation(m)
    v = rate(concs, params, Keq)
    flipped(nt) = NamedTuple{reverse(keys(nt))}(reverse(values(nt)))
    row = merge((Rate = 1.0, group = 2), flipped(concs))
    for (c, p) in ((flipped(concs), params), (concs, flipped(params)),
                   (row, all_params), (row, flipped(all_params)))
        @test rate(c, p, Keq) == v
        allocs, _ = _testhelper_test_rate_equation_performance(rate, p, Keq, c)
        @test allocs == 0
    end
    x1 = @enzyme_mechanism begin
        substrates: x1
        products: P
        steps: begin
            E + x1 ⇌ Ex1
            Ex1 <--> EP
            EP ⇌ E + P
        end
    end
    _testhelper_check_law_against_oracle(x1)
    params, Keq, concs, _ = _testhelper_random_independent_params_concs(
        x1, collect(metabolites(x1)); rng = Random.MersenneTwister(92))
    allocs, _ = _testhelper_test_rate_equation_performance(compile_rate_equation(x1),
                                                           params, Keq, concs)
    @test allocs == 0
end

@testset "a cube-root dependent emits x^(one(x)/3) and matches the oracle in BigFloat" begin
    # One isomerization group, crossed three times by the only cycle: its reverse rate is
    # the cube root of a power product (Keq^(-1/3)·…). The body computes each root once as
    # x^(one(x)/3), at BigFloat precision, where a Float64 exponent 1/3 would err by 1e-17.
    m = @enzyme_mechanism begin
        substrates: S, B
        products: P
        steps: begin
            (E <--> F, ES <--> FS, EP <--> FP)
            F + S <--> ES
            FS + B <--> EP
            FP <--> E + P
        end
    end
    @test any(d -> denominator(d.keq) == 3, EnzymeRates.Derivation._solve(m).dependents)
    roots = [(name, rhs) for (name, rhs) in _testhelper_body_locals(m)
             if occursin(r"^x\d+[mp]3$", String(name))]
    @test !isempty(roots)
    for (name, rhs) in roots
        x = Symbol(match(r"^(x\d+)", String(name)).captures[1])
        root = :($x^(one($x) / 3))
        @test rhs == (endswith(String(name), "p3") ? root : :(inv($root)))
    end
    _testhelper_check_law_against_oracle(m; bigfloat = true, seed = 93)
end

@testset "a half-integer dependent emits sqrt and inv(sqrt)" begin
    # The conformational square of "rate_equation_string runs a half-integer dependent as
    # written": the Haldane relation fixes a rate constant through square roots, so the
    # body computes each root once as sqrt(x) and each reciprocal root as inv(sqrt(x)).
    half = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (ES <--> FS, E <--> F)
            FS <--> E + P
            F + S <--> ES
        end
    end
    roots = [(name, rhs) for (name, rhs) in _testhelper_body_locals(half)
             if occursin(r"^x\d+[mp]2$", String(name))]
    @test Set(String(name)[(end - 1):end] for (name, _) in roots) == Set(["m2", "p2"])
    for (name, rhs) in roots
        x = Symbol(match(r"^(x\d+)", String(name)).captures[1])
        @test rhs == (endswith(String(name), "p2") ? :(sqrt($x)) : :(inv(sqrt($x))))
    end
end

@testset "dependent locals are numbered by solve column" begin
    # The program defines a local for each dependent that is a power product of two or more
    # symbols before any other local, in the solve's column order, so a body's numbering
    # never follows a hash; a dependent equal to one symbol takes no local. The random-order
    # bi-bi has a Haldane dependent and a Wegscheider one, each a product of several
    # constants.
    m = only(s for s in MECHANISM_TEST_SPECS
             if s.name == "Random-order Bi-Bi").mechanism
    s = EnzymeRates.Derivation._solve(m)
    L = EnzymeRates.Derivation._law(m, s)
    p = EnzymeRates.Derivation._program(m, s, L)
    @test issorted([d.col for d in s.dependents])
    function _testhelper_power(d, s, L)
        EnzymeRates.Derivation._monomial(
            L.vars.index[:Keq] => d.keq,
            (L.vars.index[s.names[f]] => e for (f, e) in d.factors)...)
    end
    products = unique!([_testhelper_power(d, s, L) for d in s.dependents
                        if length(_testhelper_power(d, s, L)) >= 2])
    @test length(s.dependents) == length(products) == 2
    @test [q for (_, q) in p.defs[1:length(products)]] ==
          EnzymeRates.Derivation._term.(products)
    # The labels fixture of "_solve labels kinetic-group, Haldane and Wegscheider
    # dependents" (test_thermodynamic_constraints.jl): a kinetic-group dependent and a
    # Wegscheider one that each equal one constant, read directly, and a Haldane
    # product, the program's first local. Every later local holds a sum of two or more
    # terms.
    folded = @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB)
            E + B ⇌ EB
            EA + B ⇌ EAB
            EAB <--> EP
            EP ⇌ E + P
        end
    end
    sf = EnzymeRates.Derivation._solve(folded)
    Lf = EnzymeRates.Derivation._law(folded, sf)
    pf = EnzymeRates.Derivation._program(folded, sf, Lf)
    @test [d.label for d in sf.dependents] == [:haldane, :kinetic_group, :wegscheider]
    haldane = only(d for d in sf.dependents if d.label === :haldane)
    @test first(pf.defs).second ==
          EnzymeRates.Derivation._term(_testhelper_power(haldane, sf, Lf))
    @test count(d -> EnzymeRates.Derivation._nterms(d.second) == 1, pf.defs) == 1
end

@testset "the generated body is finite at zero concentrations" begin
    # Every term of the program holds each concentration to a nonnegative power, so with
    # one concentration at zero the body returns the law's value there, which the law in
    # lowest terms gives exactly; "Zero-metabolite finiteness" checks finiteness alone.
    for spec in MECHANISM_TEST_SPECS
        @testset "$(spec.name)" begin
            m = spec.mechanism
            L = EnzymeRates.Derivation._memo_law(m)
            rng = Random.MersenneTwister(778 + hash(spec.name) % 1000)
            params, Keq = _testhelper_random_independent_params(m; rng)
            mets = collect(metabolites(m))
            rate = compile_rate_equation(m)
            for zeroed in mets
                concs = NamedTuple{Tuple(mets)}(
                    Tuple(n == zeroed ? 0.0 : 0.5 + rand(rng) for n in mets))
                v = rate(concs, params, Keq)
                env = Dict{Symbol, BigFloat}(
                    k => BigFloat(x) for (k, x) in pairs(merge(params, concs, (; Keq))))
                limit = _testhelper_evaluate_law(L, env)
                @test isfinite(v) && isapprox(v, Float64(limit); rtol = 1e-8)
            end
        end
    end
end

@testset "256-bit inputs flow through a four-subunit MWC body" begin
    # The mass-action check reads the allosteric specs' bodies in BigFloat at one
    # catalytic subunit only. Here the body with its D^(n−1), D^n and site powers at
    # n = 4 reads inputs that use all 256 bits and equals the law there to 1e-20, which a
    # Float64 step anywhere in it would miss.
    for name in ("PFK-1", "MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites")
        @testset "$name" begin
            m = only(s for s in MECHANISM_TEST_SPECS if s.name == name).mechanism
            L = EnzymeRates.Derivation._memo_law(m)
            @test L.n == 4
            rng = Random.MersenneTwister(779 + hash(name) % 1000)
            setprecision(BigFloat, 256) do
                params, Keq = _testhelper_random_independent_params(m; rng)
                pbig = map(x -> x * (1 + rand(rng, BigFloat)), params)
                keqbig = Keq * (1 + rand(rng, BigFloat))
                mets = metabolites(m)
                cbig = NamedTuple{Tuple(mets)}(Tuple(big(0.5) + rand(rng, BigFloat)
                                                     for _ in mets))
                env = Dict{Symbol, BigFloat}(pairs(merge(pbig, cbig, (Keq = keqbig,))))
                law = _testhelper_evaluate_law(L, env)
                @test abs(compile_rate_equation(m)(cbig, pbig, keqbig) - law) <=
                      big(10)^-20 * abs(law)
            end
        end
    end
end

"""Whether each reference line of `m`'s printed law, read back in 256-bit BigFloat
(`_testhelper_eval_rate_lines`), equals the dependent constant
`dependent_parameter_values(m, params, Keq)` computes, to 1e-12 relative, at a random
point; false when the law prints no reference line."""
function _testhelper_reference_lines_match(m)
    rng = Random.MersenneTwister(4)
    p, Keq = _testhelper_random_independent_params(m; rng)
    c = NamedTuple{Tuple(metabolites(m))}(Tuple(0.5 + rand(rng) for _ in metabolites(m)))
    env = _testhelper_eval_rate_lines(rate_equation_string(m), p, c, Keq)
    ap = dependent_parameter_values(m, p, Keq)
    dependents = collect(keys(ap))
    !isempty(dependents) &&
        all(k -> abs(env[k] - ap[k]) <= 1e-12 * abs(ap[k]), dependents)
end

@testset "rate_equation_string prints every tie as a constraint line" begin
    # A binds E and EB under one constant while B's two bindings keep their own, so
    # the RE binding square ties B's two dissociation constants. The occupied binding's
    # constant is the dependent one, and v reads it.
    folded = @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB)
            E + B ⇌ EB
            EA + B ⇌ EAB
            EAB <--> EP
            EP ⇌ E + P
        end
    end
    @test occursin("# Wegscheider constraints:\nK_EAB_to_EA_B = K_EB_to_E_B\n",
                   rate_equation_string(folded))
    _testhelper_check_law_against_oracle(folded)
    @test _testhelper_reference_lines_match(folded)
    # S binds E and Estar under one constant, so the square of the two isomerizations
    # ties their equilibrium constants; the isomerization neither of whose forms is free E
    # is the dependent one. Both run the way the reaction turns over, E → Estar and
    # EstarS → ES, so round the square the bound one's constant is the inverse.
    iso_tie = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E + S ⇌ ES, Estar + S ⇌ EstarS)
            E ⇌ Estar
            ES ⇌ EstarS
            ES <--> EP
            EP ⇌ E + P
        end
    end
    @test occursin("\nK_EstarS_to_ES = 1 / K_E_to_Estar\n", rate_equation_string(iso_tie))
    _testhelper_check_law_against_oracle(iso_tie)
    @test _testhelper_reference_lines_match(iso_tie)
    # I binding E and I binding Estar share one group and run opposite ways around the
    # cycle, so the Haldane relation fixes the Theorell–Chance constant at Keq; v reads it.
    keq_tie = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ Estar + P
            (E + I <--> EI, Estar + I <--> EI)
        end
    end
    @test occursin("\nK_E_S_to_Estar_P = Keq\n", rate_equation_string(keq_tie))
    _testhelper_check_law_against_oracle(keq_tie)
    @test _testhelper_reference_lines_match(keq_tie)
end

@testset "rate_equation_string prints residual names that run as written" begin
    # Every step of this ping-pong is at steady state, so the Haldane product holds the
    # rate constants of the residual-bearing forms. Their names are Julia identifiers, so
    # the printed law, run as written, equals the compiled law, which the oracle checks.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP <--> E_res_plusA_minusP + P
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
            EQ <--> E + Q
        end
    end
    s = rate_equation_string(m)
    @test occursin(
        r"\nk_EQ_to_EB_res_plusA_minusP = k_EA_to_EP_res_plusA_minusP \* ", s)
    @test _testhelper_reference_lines_match(m)
    @test !occursin("var\"", s)
    params, Keq = _testhelper_random_independent_params(m; rng = Random.MersenneTwister(6))
    concs = (A = 0.7, B = 1.3, P = 0.4, Q = 0.9)
    @test isapprox(_testhelper_eval_rate_string(s, params, concs, Keq),
                   compile_rate_equation(m)(concs, params, Keq); rtol = 1e-12)
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
            (ES <--> FS, E <--> F)
            FS <--> E + P
            F + S <--> ES
        end
    end
    s = rate_equation_string(half)
    @test occursin("\nk_F_to_E = k_E_to_F * sqrt(k_FS_to_E_P) * sqrt(k_F_S_to_ES) / " *
                   "(sqrt(k_ES_to_F_S) * sqrt(k_E_P_to_FS) * sqrt(Keq))\n", s)
    @test _testhelper_reference_lines_match(half)
    rng = Random.MersenneTwister(4)
    point = setprecision(() -> _testhelper_consistent_point(half, rng), BigFloat, 256)
    fitted = _testhelper_dependent_param_exprs(half)[2]
    params = NamedTuple{Tuple(fitted)}(Float64.(getindex.(Ref(point.vals), fitted)))
    Keq = Float64(point.Keq)
    concs = (S = 0.7, P = 1.3)
    @test isapprox(_testhelper_eval_rate_string(s, params, concs, Keq),
                   compile_rate_equation(half)(concs, params, Keq); rtol = 1e-12)
end

"""The deepest parenthesis level of the printed text `t` that holds a `+` or a `-`; -1
when it holds neither."""
function _testhelper_sum_depth(t::AbstractString)
    depth, worst = 0, -1
    for c in t
        depth += (c == '(') - (c == ')')
        c in ('+', '-') && (worst = max(worst, depth))
    end
    worst
end

@testset "rate_equation_string prints Q, H and D as flat sums" begin
    # Each of Q, H and D prints as one sum of terms, one per concentration monomial. A
    # coefficient that is a polynomial prints in parentheses inside its term, and none
    # holds another sum, so no sum sits deeper than one parenthesis. Random-order Bi-Bi is
    # steady state throughout, so its coefficients are polynomials; H holds none.
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "Random-order Bi-Bi")
    s = last(split(rate_equation_string(spec.mechanism), '\n'))
    @test startswith(s, "v = (")
    q, h, d = _testhelper_extract_factors(s)
    @test (_testhelper_sum_depth(q), _testhelper_sum_depth(h), _testhelper_sum_depth(d)) ==
          (1, 0, 1)
    @test h == "A * B - P * Q / Keq"
    @test _testhelper_printed_law_error(spec.mechanism) <= big(10)^-20
end

# ── The derivation memo ───────────────────────────────────────────────────────────────

@testset "per-mechanism questions answer as on a fresh mechanism, in any order" begin
    # Identify asks one concrete mechanism several questions in a row. Each answer (a
    # value, or the error the call throws) equals the answer on a fresh copy of the
    # mechanism asked alone: in order, in reverse, and with two mechanisms alternated.
    # A mechanism the guard rejects (a bottomless rapid-equilibrium segment) is among
    # them.
    D = EnzymeRates.Derivation
    rejected = EnzymeRates.DSL._construct(EnzymeRates.DSL.raw_mechanism(:(begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            EPQ ⇌ EQ + P
            EPQ ⇌ EP + Q
            EP <--> E + P
            EQ <--> E + Q
        end
    end)))
    mechs = [[s.mechanism for s in MECHANISM_TEST_SPECS]; rejected]
    questions = (m -> D.check_valid_rate_equation(Bool, m),
                 D.check_valid_rate_equation,
                 m -> parameters(m), rate_equation_string,
                 m -> rate_equation_string(m; coefficients = true),
                 D.rate_equation_dedup_key)
    answer(q, m) = try q(m) catch e; e end
    expected = [[answer(q, deepcopy(m)) for q in questions] for m in mechs]
    @test expected[end][1] === false && expected[end][2] isa ErrorException
    for (m, want) in zip(mechs, expected)
        @test [answer(q, m) for q in questions] == want
        again = deepcopy(m)
        @test reverse([answer(q, again) for q in reverse(questions)]) == want
    end
    for i in 1:(length(mechs) - 1)
        a, b = deepcopy(mechs[i]), deepcopy(mechs[i + 1])
        got = [(answer(q, a), answer(q, b)) for q in questions]
        @test first.(got) == expected[i] && last.(got) == expected[i + 1]
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
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP)
            (E + Q ⇌ EQ,
             EB_res_plusA_minusP + Q ⇌ EBQ_res_plusA_minusP)
            (E + R1 ⇌ ER1,
             EB_res_plusA_minusP + R1 ⇌ EBR1_res_plusA_minusP,
             EP_res_plusA_minusP + R1 ⇌ EPR1_res_plusA_minusP)
            (EA + P ⇌ EAP,
             E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP)
            EA <--> EP_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            (EQ + B ⇌ EBQ,
             E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP)
        end
    end

    metabs = Tuple(metabolites(m))
    params_tup = parameters(m)
    concs = NamedTuple{metabs}(ones(length(metabs)))
    pvals = NamedTuple{params_tup}(ones(length(params_tup)))

    t_compile = @elapsed begin
        v = compile_rate_equation(m)(concs, pvals, 1.0)
        s = rate_equation_string(m)
    end
    @test v isa Float64
    @test s isa String
    @test t_compile < 20.0
end

@testset "rate_equation_term_count (V×τ term-count estimate)" begin
    # V×τ = the King–Altman denominator term count, computed pre-derivation from
    # the catalytic segment graph. Values are the exact denominator product counts.
    ordered_bibi = @enzyme_mechanism begin
        substrates: S1, S2
        products: P1, P2
        steps: begin
            E + S1 <--> ES1
            ES1 + S2 <--> ES1S2
            ES1S2 <--> EP1P2
            EP1P2 <--> EP1 + P2
            EP1 <--> E + P1
        end
    end
    @test EnzymeRates.Derivation.rate_equation_term_count(ordered_bibi) == 25

    random_bibi = @enzyme_mechanism begin
        substrates: S1, S2
        products: P1, P2
        steps: begin
            E + S1 <--> ES1
            E + S2 <--> ES2
            ES1 + S2 <--> ES1S2
            ES2 + S1 <--> ES1S2
            ES1S2 <--> EP1P2
            EP1P2 <--> EP1 + P2
            EP1P2 <--> EP2 + P1
            EP1 <--> E + P1
            EP2 <--> E + P2
        end
    end
    @test EnzymeRates.Derivation.rate_equation_term_count(random_bibi) == 336

    ordered_terter = @enzyme_mechanism begin
        substrates: S1, S2, S3
        products: P1, P2, P3
        steps: begin
            E + S1 <--> ES1
            ES1 + S2 <--> ES1S2
            ES1S2 + S3 <--> ES1S2S3
            ES1S2S3 <--> EP1P2P3
            EP1P2P3 <--> EP1P2 + P3
            EP1P2 <--> EP1 + P2
            EP1 <--> E + P1
        end
    end
    @test EnzymeRates.Derivation.rate_equation_term_count(ordered_terter) == 49

    # Allosteric: both conformations count. An all-SS uni-uni has G = 3 segments and
    # τ = 3 spanning trees, 9, and its inactive conformation keeps every step: 9 more.
    allo = @enzyme_mechanism begin
        substrates: F6P
        products: F16BP
        oligomeric_state: 2
        steps: begin
            E + F6P <--> EF6P      :: EqualAI
            EF6P <--> EF16BP     :: EqualAI
            EF16BP <--> E + F16BP  :: EqualAI
        end
        regulatory_site(2): I::OnlyI
    end
    runs = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(allo)
    @test _testhelper_laplacian_tree_count(allo.catalytic, runs.active) == 9
    @test _testhelper_laplacian_tree_count(allo.catalytic, runs.inactive) == 9
    @test EnzymeRates.Derivation.rate_equation_term_count(allo) == 18
end

@testset "Rate equation too large error" begin
    # Manually defined mechanism (11 forms, 16 steps; V×τ ≈ 29k denominator
    # terms) exceeds MAX_RATE_EQUATION_TERMS, so the upfront V×τ check in `_law`
    # aborts the derivation before it builds the polynomial.
    m_manual = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: R1
        steps: begin
            E + A <--> EA
            EA <--> FP
            FP <--> F + P
            F + B <--> FB
            FB <--> EQ
            EQ <--> E + Q
            E + R1 <--> ER1
            EA + R1 <--> EAR1
            FP + R1 <--> FPR1
            F + R1 <--> FR1
            FB + R1 <--> FBR1
            ER1 + A <--> EAR1
            EAR1 <--> FPR1
            FPR1 <--> FR1 + P
            FR1 + B <--> FBR1
            FBR1 <--> ER1 + Q
        end
    end
    @test_throws "polynomial terms" rate_equation_string(m_manual)
    @testset "a user catches the error by its public name" begin
        # Only exported and public names, written as a user outside the package writes them.
        TooLarge = EnzymeRates.Derivation.RateEquationTooLargeError
        @test Base.ispublic(EnzymeRates.Derivation, :RateEquationTooLargeError)
        @test_throws TooLarge rate_equation_string(m_manual)
        @test_throws TooLarge EnzymeRates.Derivation.rate_equation_dedup_key(m_manual)
        ones_for(names) = NamedTuple{Tuple(names)}(ntuple(_ -> 1.0, length(names)))
        rate = compile_rate_equation(m_manual)
        @test_throws TooLarge rate(ones_for(metabolites(m_manual)),
                                   ones_for(parameters(m_manual)), 2.0)
    end
end

# The generated call of a compiled rate law names its type parameter M_type, which holds
# `_` as no metabolite name does; metabolites named M and m read as concentrations.
@testset "a metabolite named M or m reads as a concentration" begin
    upper = @enzyme_mechanism begin
        substrates: M
        products:   P
        steps: begin
            E + M ⇌ EM
            EM <--> EP
            EP ⇌ E + P
        end
    end
    _testhelper_check_law_against_oracle(upper)
    lower = @enzyme_mechanism begin
        substrates: m
        products:   P
        steps: begin
            E + m ⇌ Em
            Em <--> EP
            EP ⇌ E + P
        end
    end
    _testhelper_check_law_against_oracle(lower)
end

# The names no metabolite takes: Keq, L, N, params, concs, v and every name
# holding `_` (test_enzyme_reaction_type.jl checks that the reaction rejects each). Besides
# the metabolites, the printed law binds or reads only such names; the generated body binds
# only positional locals (the next testset), so no concentration can overwrite a constant
# or a local.
@testset "the reserved names cover every non-metabolite name of a law" begin
    reserved = (:Keq, :L, :N, :params, :concs, :v)
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        printed = Meta.parseall(rate_equation_string(m))
        coefficients = Meta.parseall(rate_equation_string(m; coefficients = true))
        for p in (printed, coefficients)
            others = setdiff(_testhelper_names_read(p), metabolites(m))
            # N, N_A, Q_A, D_A, F_A, W_A, n_0, d_A²_B, nA_0 and the rest hold `_` or are N.
            @test all(n -> occursin('_', String(n)) || n in reserved, others)
            @test :v in _testhelper_names_read(p)
        end
        # The full law destructures params and concs and reads Keq.
        @test (:params, :concs, :Keq) ⊆ setdiff(_testhelper_names_read(printed),
                                               metabolites(m))
    end
end

@testset "the generated body binds positional locals and reads inputs by getfield" begin
    # Each local is x<i>, its reciprocal x<i>m<q> or its root x<i>p<q>, and each input is
    # read as getfield(concs, :X) of a metabolite, getfield(params, :k) of a parameter, or
    # the argument Keq, so no metabolite name can shadow a constant or a local, and Keq is
    # never read from params.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        body = _testhelper_rate_body(m)
        locals = _testhelper_body_locals(m)
        @test all(p -> occursin(r"^x\d+([mp]\d+)?$", String(p.first)), locals)
        @test _testhelper_names_read(body) ⊆
              union(Set(first.(locals)), [:concs, :params, :Keq])
        @test :Keq in _testhelper_names_read(body)
        reads = Expr[]
        function _testhelper_getfields!(e)
            e isa Expr || return
            e.head === :call && e.args[1] === :getfield && push!(reads, e)
            foreach(_testhelper_getfields!, e.args)
        end
        _testhelper_getfields!(body)
        @test all(r -> r.args[2] in (:concs, :params) && r.args[3] isa QuoteNode, reads)
        @test Set(r.args[3].value for r in reads if r.args[2] === :concs) ⊆
              Set(metabolites(m))
        @test Set(r.args[3].value for r in reads if r.args[2] === :params) ⊆
              Set(parameters(m))
    end
end


# ── Single-feature edge cases ─────────────────────────────────────────────
@testset "Allosteric edge cases" begin
    # OnlyA substrate + OnlyA catalysis: S binds only in the R-state and only
    # the R-state catalyzes (k_T = 0). The T-state therefore carries no flux
    # (N_T = 0) and never populates ES, so no reverse-catalysis weight leaks
    # into Q_I. As K_A_ES_to_E_S → ∞ (weaker R-state binding), rate vanishes.
    onlyR_sub = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES     :: OnlyA
            ES <--> EP   :: OnlyA
            EP ⇌ E + P     :: EqualAI
        end
    end
    concs = (S=1.0, P=0.001)
    base_params = (k_A_ES_to_EP=10.0, K_EP_to_E_P=0.5, L=10.0)
    onlyR_rate = compile_rate_equation(onlyR_sub)
    rate_strong = onlyR_rate(concs, merge(base_params, (K_A_ES_to_E_S=0.01,)), 1000.0)
    rate_weak   = onlyR_rate(concs, merge(base_params, (K_A_ES_to_E_S=1e6,)), 1000.0)
    @test rate_strong > 1.0
    @test rate_weak < 1e-3
    @test rate_weak / rate_strong < 1e-5

    # V-type only: all bindings :EqualAI but catalysis :OnlyA. T-state binds
    # substrate normally but cannot catalyze (k_T = 0), so N_T = 0. As L → ∞
    # (T-state dominant) the rate vanishes.
    vtype = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES     :: EqualAI
            ES <--> EP   :: OnlyA
            EP ⇌ E + P     :: EqualAI
        end
    end
    vparams = (K_ES_to_E_S=0.1, k_A_ES_to_EP=10.0, K_EP_to_E_P=0.5)
    vtype_rate = compile_rate_equation(vtype)
    rate_R = vtype_rate(concs, merge(vparams, (L=0.0,)), 1000.0)
    rate_T = vtype_rate(concs, merge(vparams, (L=1e10,)), 1000.0)
    @test rate_R > 1.0
    @test rate_T < 1e-6
    # T-state numerator branch is elided when the T-state cycle is dead (:OnlyA catalysis);
    # rate is num_R / (Q_R^catN + L · Q_T^catN). At large L, the T-state
    # enzyme mass dominates the denominator → rate ∝ 1/(1+L).
    @test rate_T * 1e10 < 100.0    # bounded as L grows

    # The inactive conformation's binding constants are dissociation constants (`S / K`)
    # also when `:OnlyA` and `:NonequalAI` catalytic groups coexist; read as association
    # constants (`K * S`) they would halve the rate below.
    # S binding :NonequalAI, catalysis :OnlyA, P binding :NonequalAI.
    m_mix = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: OnlyA
            EP ⇌ E + P      :: NonequalAI
        end
        regulatory_site(2): I::OnlyI
    end
    p_mix = (K_A_ES_to_E_S=0.1, k_A_ES_to_EP=10.0, K_A_EP_to_E_P=0.5,
             K_I_ES_to_E_S=10.0, K_I_EP_to_E_P=10.0,
             K_I_Ireg=1.0, L=1.0)
    rate_mix = compile_rate_equation(m_mix)((S=10.0, P=0.0, I=0.0), p_mix, 1000.0)
    # With Kd convention (correct): rate ≈ 9.90 (R-state catalysis dominates),
    # per active site (catN = 2).
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
    # at chemical equilibrium. (Over-parametrized; the enumerator emits such
    # degenerate configs.)
    # B binding (step 2) :NonequalAI, rest :EqualAI.
    m_ro = @enzyme_mechanism begin
        substrates: A, B
        products:   P, Q
        oligomeric_state: 2
        steps: begin
            E + A <--> EA            :: EqualAI
            E + B <--> EB            :: NonequalAI
            EA + B <--> EAB      :: EqualAI
            EB + A <--> EAB      :: EqualAI
            EAB <--> EPQ       :: EqualAI
            EPQ <--> EQ + P      :: EqualAI
            EQ <--> E + Q            :: EqualAI
        end
        regulatory_site(2): I::NonequalAI
    end
    # Overall A + B ⇌ P + Q, so Keq = P·Q/(A·B).
    Keq_ro = 4.0
    A_eq, B_eq = 1.5, 2.0
    P_eq = 3.0
    Q_eq = Keq_ro * A_eq * B_eq / P_eq
    # Which koff is Wegscheider-dependent follows the solve's priority. The values
    # come from the mechanism's own fitted parameters, so the Haldane property (zero net
    # rate at equilibrium) is verified whatever the priority picks.
    fp_ro = parameters(m_ro)
    rng_ro = Random.MersenneTwister(123)
    p_ro = NamedTuple{fp_ro}(ntuple(_ -> 0.5 + rand(rng_ro), length(fp_ro)))
    @test isapprox(
        compile_rate_equation(m_ro)((A=A_eq, B=B_eq, P=P_eq, Q=Q_eq, I=0.5), p_ro, Keq_ro),
        0.0; atol=1e-9)
end


@testset "rate_equation_string allosteric byte-identical fixture" begin
    m_allo = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP     :: NonequalAI
            E + S ⇌ ES     :: NonequalAI
            ES <--> EP   :: NonequalAI
        end
        regulatory_site(2): R::NonequalAI
    end
    actual = rate_equation_string(m_allo)
    expected = raw"""(; K_A_EP_to_E_P, K_A_ES_to_E_S, K_I_EP_to_E_P, K_I_ES_to_E_S, k_A_ES_to_EP, k_I_ES_to_EP, K_A_Rreg, K_I_Rreg, L) = params
(; S, P, R) = concs
# Haldane constraints:
k_A_EP_to_ES = K_A_EP_to_E_P * k_A_ES_to_EP / (K_A_ES_to_E_S * Keq)
k_I_EP_to_ES = K_I_EP_to_E_P * k_I_ES_to_EP / (K_I_ES_to_E_S * Keq)
Q_A = k_A_ES_to_EP / K_A_ES_to_E_S
D_A = 1 + S / K_A_ES_to_E_S + P / K_A_EP_to_E_P
Q_I = k_I_ES_to_EP / K_I_ES_to_E_S
D_I = 1 + S / K_I_ES_to_E_S + P / K_I_EP_to_E_P
W_A = (1 + R / K_A_Rreg)^2
W_I = (1 + R / K_I_Rreg)^2
v = (Q_A * D_A * W_A + L * Q_I * D_I * W_I) * (S - P / Keq) / (D_A^2 * W_A + L * D_I^2 * W_I)"""
    @test actual == expected
    expected_coefficients = raw"""(; K_A_EP_to_E_P, K_A_ES_to_E_S, K_I_EP_to_E_P, K_I_ES_to_E_S, k_A_ES_to_EP, k_I_ES_to_EP, K_A_Rreg, K_I_Rreg, L) = params
(; S, P, R) = concs
# The coefficients are functions of the fitted constants, not parameters:
nA_0 = k_A_ES_to_EP / K_A_ES_to_E_S
dA_S = 1 / K_A_ES_to_E_S
dA_P = 1 / K_A_EP_to_E_P
nI_0 = k_I_ES_to_EP / K_I_ES_to_E_S
dI_S = 1 / K_I_ES_to_E_S
dI_P = 1 / K_I_EP_to_E_P
Q_A = nA_0
D_A = 1 + dA_S * S + dA_P * P
Q_I = nI_0
D_I = 1 + dI_S * S + dI_P * P
W_A = (1 + R / K_A_Rreg)^2
W_I = (1 + R / K_I_Rreg)^2
v = (Q_A * D_A * W_A + L * Q_I * D_I * W_I) * (S - P / Keq) / (D_A^2 * W_A + L * D_I^2 * W_I)"""
    @test rate_equation_string(m_allo; coefficients = true) == expected_coefficients
    # Both texts equal the derived law (value check of this pin).
    @test _testhelper_printed_law_error(m_allo) <= big(10)^-20
end

@testset "Fix A: dead-inactive-state allosteric body defines all I-state symbols" begin
    # Random-order allosteric bi-bi with an :OnlyA catalytic step → dead inactive
    # state. Verified pre-fix to crash with an UndefVarError on the I-state
    # release rate of A from EA.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A <--> EA        :: NonequalAI
            E + B <--> EB        :: NonequalAI
            EA + B <--> EAB  :: NonequalAI
            EB + A <--> EAB  :: NonequalAI
            EAB <--> EPQ   :: OnlyA
            EPQ <--> EP + Q  :: NonequalAI
            EPQ <--> EQ + P  :: NonequalAI
            EP <--> E + P        :: NonequalAI
            EQ <--> E + Q        :: NonequalAI
        end
    end
    pn = parameters(m)
    params = NamedTuple{pn}(ntuple(_ -> 1.0, length(pn)))
    concs = (A = 1.0, B = 1.0, P = 1.0, Q = 1.0)
    @test isfinite(compile_rate_equation(m)(concs, params, 1.0))
end

# ── The Sig ───────────────────────────────────────────────────────────────────────────

"""Whether every leaf of the nested tuple `x`, a Sig or a part of one, is a Symbol, an Int
or a Bool."""
_testhelper_sig_leaves(x) =
    x isa Tuple ? all(_testhelper_sig_leaves, x) : x isa Union{Symbol, Int, Bool}

@testset "a mechanism's Sig decodes to the mechanism" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        sig = ER.Derivation._sig_of(x)
        @test sig isa Tuple && _testhelper_sig_leaves(sig)
        @test _testhelper_decode(sig) == x
        rate = compile_rate_equation(x)
        @test typeof(rate).parameters[1] === sig
        @test _testhelper_decode(typeof(rate).parameters[1]) == x
    end
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ EP
            E + S ⇌ ES
            ES <--> EP
        end
    end
    @test _testhelper_decode(ER.Derivation._sig_of(m)) == m
    am = ER.Mechanism(
        m.catalytic, ER.Mechanisms.Allostery([:EqualAI, :NonequalAI, :OnlyA], 4,
        [ER.Mechanisms.RegulatorySite([:X, :R], 2, [:OnlyA, :OnlyI])]))
    @test _testhelper_decode(ER.Derivation._sig_of(am)) == am
    # A name outside ASCII, sorting last among 31 substrates.
    names = [[Symbol("S", i) for i in 1:30]; :αKG; :P]
    wide = _testhelper_plain_mechanism(
        names, [:E; [Symbol("E", x) for x in names]], collect(1:31), [32],
        [ER.Mechanisms.Step(1, i + 1, [i], Int[], true, i) for i in 1:32])
    @test wide.catalytic.metabolites[31] === :αKG
    @test _testhelper_decode(ER.Derivation._sig_of(wide)) == wide
end

"""Uni-uni with `n` modifiers, `I1` to `In`, each binding free E alone at rapid
equilibrium: a Sig whose type no other test decodes, its length set by `n`."""
function _testhelper_n_modifiers(n)
    mods = [Symbol("I", i) for i in 1:n]
    forms = [:E, :EP, :ES, [Symbol("E", x) for x in mods]...]
    steps = [ER.Mechanisms.Step(1, 3, [1], Int[], true, 1),
             ER.Mechanisms.Step(3, 2, Int[], Int[], false, 2),
             ER.Mechanisms.Step(1, 2, [2], Int[], true, 3),
             [ER.Mechanisms.Step(1, 3 + i, [2 + i], Int[], true, 3 + i) for i in 1:n]...]
    _testhelper_plain_mechanism([:S; :P; mods], forms, [1], [2], steps)
end

@testset "the Sig decoder compiles once, whatever the mechanism" begin
    # A compiled law's generated methods decode its Sig without specializing on its type,
    # so the decoder compiles once: after one decode, decoding the Sigs of three other
    # mechanisms, each of a Sig type of its own, adds no method instance.
    instances(f) = count(_ -> true, Base.specializations(only(methods(f))))
    _testhelper_allosteric_n_modifiers(n) =
        ER.Mechanism(_testhelper_n_modifiers(n).catalytic,
                     ER.Mechanisms.Allostery(fill(:EqualAI, n + 3), 2,
                                             ER.Mechanisms.RegulatorySite[]))
    decoder = ER.Derivation._mechanism_from_sig
    for build in (_testhelper_n_modifiers, _testhelper_allosteric_n_modifiers)
        counts = 41:44
        sigs = [typeof(compile_rate_equation(build(n))).parameters[1] for n in counts]
        @test allunique(typeof.(sigs))
        @test _testhelper_decode(sigs[1]) == build(first(counts))
        before = instances(decoder)
        foreach(_testhelper_decode, sigs[2:4])
        @test instances(decoder) == before
    end
end

@testset "_sig_of / _mechanism_from_sig roundtrip" begin
    # One fixture built from the walker's value: two substrates, two products, a
    # modifier-only name bound by a step, a covalent residual, a Theorell–Chance step,
    # and steps at rapid equilibrium and at steady state. It holds no reactant's
    # modifier copy.
    m = ER.DSL._construct(ER.DSL.raw_mechanism(:(begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA <--> E_res_plusA_minusP + P
            E_res_plusA_minusP + B <--> EQ
            EA + B <--> EQ + P
            E + Q ⇌ EQ
            E + I ⇌ EI
        end
    end)))

    sig = ER.Derivation._sig_of(m)
    @test sig isa Tuple && _testhelper_sig_leaves(sig)
    # Its leaves are Symbols for the names, Ints for indices and the oligomeric state, and
    # Bools for the rapid-equilibrium flags.
    @test sig[1] === (:A, :B, :P, :Q, :I) && sig[7] === 1
    @test all(s -> s[5] isa Bool, sig[5])
    @test ER.Derivation._mechanism_from_sig(sig) == m   # roundtrip

    # The Sig is the compiled rate law's only parameter.
    rate = compile_rate_equation(m)
    @test typeof(rate) === ER.Derivation._CompiledRateEquation{sig}

    # Roundtrip through the type-parameter form preserves everything.
    @test ER.Derivation._mechanism_from_sig(typeof(rate).parameters[1]) == m

    # A regulator that no step binds is no part of a mechanism, so the init mechanism
    # of a reaction that declares one compiles to the law of the reaction's without
    # it.
    without_u = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    with_u = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: U
    end
    @test compile_rate_equation(
        only(_testhelper_mechanisms(ER.Enumeration.init_mechanisms(with_u)))) ===
          compile_rate_equation(
              only(_testhelper_mechanisms(ER.Enumeration.init_mechanisms(without_u))))
end

# ── Polynomials ───────────────────────────────────────────────────────────────────────

@testset "monomials multiply, raise and read exponents exactly" begin
    P = EnzymeRates.Derivation
    @test P._monomial(2 => 1, 1 => 1//2, 2 => 2) == [1 => 1//2, 2 => 3//1]
    @test P._monomial(1 => 1, 1 => -1) == P._Mono()
    @test P._monomial(1 => 0, 3 => 2) == [3 => 2//1]
    @test P._mmul(P._monomial(1 => 1//2), P._monomial(1 => -1//2, 3 => 2)) == [3 => 2//1]
    @test P._mpow(P._monomial(1 => 1//2, 2 => -1), 3) == [1 => 3//2, 2 => -3//1]
    @test P._mpow(P._monomial(1 => 2), 1//2) == [1 => 1//1]
    @test P._mpow(P._monomial(1 => 2), 0) == P._Mono()
    @test P._mexp(P._monomial(1 => 1//2), 1) == 1//2
    @test P._mexp(P._monomial(1 => 1//2), 2) == 0
end

# A law splats a mechanism-dependent number of pairs into `_monomial`; if a new count
# compiled code, a fresh mechanism would compile during identify's PASS1, and which code
# it compiled would depend on compile order.
@testset "a monomial of a new number of pairs compiles nothing" begin
    P = EnzymeRates.Derivation
    if isdefined(Core, :methodtable)
        function _testhelper_instance_count()
            meths = Method[]
            Base.visit(m -> push!(meths, m), Core.methodtable)
            sum(m -> count(_ -> true, Base.specializations(m)), meths)
        end
        monomial_of(n) = P._monomial((v => 1 for v in 1:n)...)
        expected = [v => 1//1 for v in 1:29]
        monomial_of(1) == monomial_of(2) == monomial_of(3) == expected
        n0 = _testhelper_instance_count()
        m = monomial_of(29)
        n1 = _testhelper_instance_count()
        @test n1 == n0
        @test m == expected
    else
        @test_skip "counting method instances needs the global method table of Julia 1.12"
    end
end

@testset "polynomials add, subtract and multiply exactly" begin
    P = EnzymeRates.Derivation
    x, y = P._term(P._monomial(1 => 1)), P._term(P._monomial(2 => 1))
    p = x + y
    @test (p * p).terms == Dict(P._monomial(1 => 2) => 1, P._monomial(1 => 1, 2 => 1) => 2,
                                P._monomial(2 => 2) => 1)
    @test P._nterms(p * p) == 3
    @test (x - y) * (x + y) == x * x - y * y
    @test p - y == x && -x == P._scale(x, -1)
    @test iszero(p - p) && iszero(zero(P._Poly)) && !iszero(one(P._Poly))
    @test p * P._monomial(1 => -1) == one(P._Poly) + P._term(P._monomial(1 => -1, 2 => 1))
    @test P._scale(p, 0) == zero(P._Poly) && P._term(P._Mono(), 0) == zero(P._Poly)
    @test P._addto!(P._term(P._monomial(1 => 1), 2), P._monomial(1 => 1), -2) ==
          zero(P._Poly)
    # Equal by terms, hashed by terms, and distinct objects.
    @test x + y == y + x && hash(x + y) == hash(y + x) && x + y !== y + x
    @test ismutabletype(P._Poly) && isconst(P._Poly, :terms)
    @test fieldtypes(P._Poly) == (Dict{P._Mono, Int},)
end

@testset "an Int64 overflow throws OverflowError" begin
    P = EnzymeRates.Derivation
    m = P._monomial(1 => 1)
    @test_throws OverflowError P._term(m, typemax(Int)) + P._term(m, 1)
    @test_throws OverflowError P._term(m, typemax(Int) ÷ 2 + 1) * P._term(P._Mono(), 2)
    @test_throws OverflowError P._scale(P._term(m, typemax(Int)), 2)
end

@testset "_divexact divides out the Haldane factor" begin
    P = EnzymeRates.Derivation
    # Variables: 1 = A, 2 = P, 3 = Keq.
    mono(ps...) = P._term(P._monomial(ps...))
    H = mono(1 => 1) - mono(2 => 1, 3 => -1)
    q = mono(1 => 1, 3 => 1//2) + P._scale(mono(2 => 1, 3 => -1//2), 2) + one(P._Poly)
    @test P._divexact(H * q, H) == q
    # A substrate with two copies: H = A² − P/Keq.
    H2 = mono(1 => 2) - mono(2 => 1, 3 => -1)
    q2 = mono(1 => 1) + mono(2 => 1) + mono(1 => 1, 2 => 1, 3 => 1//2)
    @test P._divexact(H2 * q2, H2) == q2
    # A chain with a missing power: H·(A³ + P³/Keq³) has no A²·P²/Keq² term.
    qg = mono(1 => 3) + mono(2 => 3, 3 => -3)
    @test P._divexact(H * qg, H) == qg
    @test iszero(P._divexact(zero(P._Poly), H))
    @test_throws ArgumentError P._divexact(H * q + one(P._Poly), H)
    @test_throws ArgumentError P._divexact(q, q)
end

@testset "_subst replaces variables by power products" begin
    P = EnzymeRates.Derivation
    mono(ps...) = P._term(P._monomial(ps...))
    p = mono(1 => 1, 2 => 2) + P._scale(mono(3 => 1), 3)
    f(v) = v == 2 ? P._monomial(4 => 1//2, 5 => -1) : nothing
    @test P._subst(f, p) == mono(1 => 1, 4 => 1, 5 => -2) + P._scale(mono(3 => 1), 3)
    # Terms that substitution makes alike merge, and cancel to zero.
    g(v) = v == 2 ? P._monomial(1 => 1) : nothing
    @test P._subst(g, mono(1 => 1) + mono(2 => 1)) == P._scale(mono(1 => 1), 2)
    @test iszero(P._subst(g, mono(1 => 1) - mono(2 => 1)))
end

@testset "the shifts bring polynomials to lowest terms in the concentrations" begin
    P = EnzymeRates.Derivation
    t = P._Vars()
    a, b, k = P._var!(t, :A; conc = true), P._var!(t, :B; conc = true), P._var!(t, :k)
    @test (a, b, k) == (1, 2, 3) && P._var!(t, :A; conc = true) == 1
    @test t.names == [:A, :B, :k] && t.conc == [true, true, false]
    @test_throws ErrorException P._var!(t, :k; conc = true)
    @test_throws ErrorException("internal error: A is a concentration in one use and a " *
                                "constant in another") P._var!(t, :A)
    @test t == P._Vars([:A, :B, :k], BitVector([1, 1, 0]), Dict(:A => 1, :B => 2, :k => 3))
    mono(ps...) = P._term(P._monomial(ps...))
    Q = mono(a => 1, b => 1, k => 1) + mono(a => 2, k => -1)
    D = mono(a => 1) + mono(a => 1, b => 1)
    # A is in every term, lowest 1; B is missing from some, lowest 0; k reaches −1.
    @test P._lowest_shift((Q, D)) == P._monomial(a => -1, k => 1)
    @test P._conc_shift((Q, D), t) == P._monomial(a => -1)
    D2 = mono(a => 1, k => 1) + mono(a => 1, b => 1) + P._scale(mono(k => 1//2), 2) +
         mono(b => 1, k => 1)
    @test P._split_conc(D2, t) == Dict(P._monomial(a => 1) => mono(k => 1),
                                      P._monomial(a => 1, b => 1) => one(P._Poly),
                                      P._Mono() => P._scale(mono(k => 1//2), 2),
                                      P._monomial(b => 1) => mono(k => 1))
end

@testset "_content splits off the integer and monomial content" begin
    P = EnzymeRates.Derivation
    mono(ps...) = P._term(P._monomial(ps...))
    p = P._scale(mono(1 => 2, 3 => 1), 6) + P._scale(mono(1 => 1, 3 => 1//2), 4)
    @test P._content(p) ==
          (2, P._monomial(1 => 1, 3 => 1//2),
           P._scale(mono(1 => 1, 3 => 1//2), 3) + P._scale(one(P._Poly), 2))
    @test P._content(P._scale(mono(1 => 1), -2) + P._scale(one(P._Poly), -4)) ==
          (2, P._Mono(), -mono(1 => 1) - P._scale(one(P._Poly), 2))
    @test P._content(P._scale(mono(1 => -1, 2 => 1), 5)) ==
          (5, P._monomial(1 => -1, 2 => 1), one(P._Poly))
end

# ── The engine ────────────────────────────────────────────────────────────────────────

@testset "_spanning_trees enumerates each spanning tree once" begin
    trees(n, pairs) = Set(Set.(EnzymeRates.Derivation._spanning_trees(n, pairs)))
    @test trees(3, [(1, 2), (2, 3), (1, 3)]) == Set([Set([1, 2]), Set([1, 3]), Set([2, 3])])
    @test length(
        EnzymeRates.Derivation._spanning_trees(4, [(1, 2), (2, 3), (3, 4), (1, 4)])) == 4
    @test EnzymeRates.Derivation._spanning_trees(2, [(1, 2), (1, 2)]) == [[1], [2]]
    @test isempty(EnzymeRates.Derivation._spanning_trees(3, [(1, 2)]))
    @test EnzymeRates.Derivation._spanning_trees(1, Tuple{Int, Int}[]) == [Int[]]
end

@testset "_tree_count is the BigInt determinant, saturated at typemax(Int)" begin
    complete(n) = [i == j ? n - 1 : -1 for i in 1:n, j in 1:n]
    @test EnzymeRates.Derivation._tree_count(complete(5)) == 125
    @test EnzymeRates.Derivation._tree_count(complete(20)) == typemax(Int)
    @test EnzymeRates.Derivation._tree_count(zeros(Int, 1, 1)) == 1
    @test EnzymeRates.Derivation._tree_count([1 -1 0; -1 1 0; 0 0 0]) == 0
    @test EnzymeRates.Derivation._tree_count([1 0 -1; 0 0 0; -1 0 1]) == 0
    rng = Random.MersenneTwister(31)
    for _ in 1:50
        n = rand(rng, 2:14)
        lap = zeros(Int, n, n)
        for a in 1:n, b in (a + 1):n
            rand(rng) < 0.5 || continue
            k = rand(rng, 1:3)
            lap[a, a] += k; lap[b, b] += k; lap[a, b] -= k; lap[b, a] -= k
        end
        exact = det(BigInt.(lap[2:end, 2:end]))
        @test EnzymeRates.Derivation._tree_count(lap) ==
              (exact > typemax(Int) ? typemax(Int) : Int(exact))
    end
end

@testset "rate_equation_term_count is segments times spanning trees" begin
    n_plain = 0
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        _testhelper_is_allosteric(m) && continue
        n_plain += 1
        @test EnzymeRates.Derivation.rate_equation_term_count(m) ===
              Int(_testhelper_laplacian_tree_count(m, trues(length(m.catalytic.steps))))
    end
    @test n_plain == 31
    random_bibi = @enzyme_mechanism begin
        substrates: S1, S2
        products: P1, P2
        steps: begin
            E + S1 <--> ES1
            E + S2 <--> ES2
            ES1 + S2 <--> ES1S2
            ES2 + S1 <--> ES1S2
            ES1S2 <--> EP1P2
            EP1P2 <--> EP1 + P2
            EP1P2 <--> EP2 + P1
            EP1 <--> E + P1
            EP2 <--> E + P2
        end
    end
    @test EnzymeRates.Derivation.rate_equation_term_count(random_bibi) === 336
    one_segment = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
        end
    end
    @test EnzymeRates.Derivation.rate_equation_term_count(one_segment) === 1
end

@testset "the plain law equals the exact steady state on every plain spec" begin
    for (i, spec) in enumerate(MECHANISM_TEST_SPECS)
        m = spec.mechanism
        _testhelper_is_allosteric(m) && continue
        @testset "$(spec.name)" begin
            L = EnzymeRates.Derivation._memo_law(m)
            rng = Random.MersenneTwister(500 + i)
            for _ in 1:3
                point, concs, params =
                    _testhelper_law_inputs(m, EnzymeRates.Derivation._solve(m).fitted, rng)
                oracle = _testhelper_steady_state_rate(m, point, concs)
                @test abs(_testhelper_evaluate_law(L, merge(params, concs)) - oracle) <=
                      big(10)^-20 * abs(oracle)
            end
        end
    end
end

@testset "a plain law is H times a positive polynomial in lowest terms" begin
    rng = Random.MersenneTwister(7)
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        _testhelper_is_allosteric(m) && continue
        @testset "$(spec.name)" begin
            L = EnzymeRates.Derivation._memo_law(m)
            @test L.vars.names == [EnzymeRates.metabolites(m);
                                   EnzymeRates.Derivation._solve(m).fitted; :Keq]
            env = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1)
                                         for x in L.vars.names)
            H = prod(env[x] for x in _testhelper_substrates(m)) -
                prod(env[x] for x in _testhelper_products(m)) / env[:Keq]
            @test isapprox(_testhelper_poly_value(L.H, L.vars, env), H;
                           rtol = big(10)^-60)
            @test !iszero(L.Q) && all(>(0), values(L.Q.terms))
            @test all(>(0), values(L.D.terms))
            @test isempty(EnzymeRates.Derivation._conc_shift((L.Q, L.D), L.vars))
            @test L.D == sum(L.sigma .* L.trees;
                             init = zero(EnzymeRates.Derivation._Poly)) * L.shift
        end
    end
end

@testset "a conformation's law on its mask is the law of the steps it keeps" begin
    # An inhibitor that binds only the inactive conformation (:OnlyI). The active mask
    # drops its binding and the form EI, keeping the Michaelis–Menten steps, whose
    # constants both conformations share (:EqualAI); the inactive mask keeps all three
    # steps, the binding under its inactive constant. Both laws number their variables
    # alike: the mechanism's concentrations, the solve's independents, then Keq.
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + I ⇌ EI            :: OnlyI
        end
    end
    mm = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
        end
    end
    competitive = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            E + I ⇌ EI
        end
    end
    s = EnzymeRates.Derivation._solve(am)
    mets = EnzymeRates.metabolites(am)
    runs = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(am)
    active = EnzymeRates.Derivation._law(am, s, runs.active)
    inactive = EnzymeRates.Derivation._law(am, s, runs.inactive; inactive = true)
    @test active.vars == inactive.vars
    @test active.vars.names == [mets; s.names[s.independents]; :Keq]
    rng = Random.MersenneTwister(71)
    for _ in 1:3
        env = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1)
                                     for x in active.vars.names)
        env[:K_EI_to_E_I] = env[:K_I_EI_to_E_I]
        @test isapprox(_testhelper_evaluate_law(active, env),
                       _testhelper_evaluate_law(EnzymeRates.Derivation._memo_law(mm), env);
                       rtol = big(10)^-60)
        @test isapprox(
            _testhelper_evaluate_law(inactive, env),
            _testhelper_evaluate_law(EnzymeRates.Derivation._memo_law(competitive), env);
            rtol = big(10)^-60)
    end
end

@testset "rewriting or regrouping a mechanism changes no plain law" begin
    # Written backwards: every step the other way round, in the reverse order, which
    # canonical form gives back.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        _testhelper_is_allosteric(m) && continue
        c = m.catalytic
        back = _testhelper_plain_mechanism(c.metabolites, c.enzyme_forms, c.substrates,
                                           c.products,
                                           reverse(_testhelper_turned.(c.steps)))
        @test EnzymeRates.Derivation._memo_law(back) == EnzymeRates.Derivation._memo_law(m)
    end
    # A same-model regrouping: the binding square ties B's two dissociation constants,
    # so tying them as one group as well leaves the model and its fitted names alone.
    folded = @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB)
            E + B ⇌ EB
            EA + B ⇌ EAB
            EAB <--> EP
            EP ⇌ E + P
        end
    end
    tied = @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB)
            (E + B ⇌ EB, EA + B ⇌ EAB)
            EAB <--> EP
            EP ⇌ E + P
        end
    end
    @test folded != tied
    @test EnzymeRates.Derivation._solve(folded).fitted ==
          EnzymeRates.Derivation._solve(tied).fitted
    @test EnzymeRates.Derivation._memo_law(folded) == EnzymeRates.Derivation._memo_law(tied)
end

@testset "plain laws satisfy detailed balance" begin
    for (i, spec) in enumerate(MECHANISM_TEST_SPECS)
        m = spec.mechanism
        _testhelper_is_allosteric(m) && continue
        L = EnzymeRates.Derivation._memo_law(m)
        rng = Random.MersenneTwister(600 + i)
        env = Dict{Symbol, BigFloat}(x => exp(4 * big(rand(rng)) - 2)
                                     for x in L.vars.names)
        p1 = first(_testhelper_products(m))
        env[p1] = env[:Keq] * prod(env[x] for x in _testhelper_substrates(m)) /
                  prod((env[x] for x in _testhelper_products(m) if x != p1);
                       init = big(1))
        off = copy(env)
        off[p1] /= 10
        @test abs(_testhelper_evaluate_law(L, env)) <=
              big(10)^-30 * abs(_testhelper_evaluate_law(L, off))
    end
end

@testset "the engine refuses a rapid-equilibrium cycle that turns the reaction over" begin
    # E + S ⇌ ES ⇌ EP ⇌ E + P closes a rapid-equilibrium cycle that turns S into P
    # beside a steady-state route; the macro rejects it (class 2), and the engine,
    # reading it unchecked, refuses to derive a law for it.
    beside = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES ⇌ EP
            E + P ⇌ EP
            ES <--> FP
            FP <--> E + P
        end
    end
    msg = "internal error: a rapid-equilibrium cycle of this mechanism turns the " *
          "reaction over (class 2)"
    @test_throws ErrorException(msg) EnzymeRates.Derivation._memo_law(beside)
    err = @test_throws ErrorException @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES ⇌ EP
            E + P ⇌ EP
            ES <--> FP
            FP <--> E + P
        end
    end
    @test startswith(err.value.msg,
                     "Mechanism: the rapid-equilibrium steps E + P ⇌ EP, " *
                     "E + S ⇌ ES and ES ⇌ EP close a cycle")
end

@testset "_law refuses a mask whose steady-state steps leave its segments apart" begin
    # The mask keeps two steady-state steps that share no form, so the segment graph has
    # no spanning tree: V×τ is 0, which the cap lets through, and the engine refuses.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EAB
            EAB <--> EPQ
            EPQ <--> EQ + P
            EQ <--> E + Q
        end
    end
    mask = BitVector([EnzymeRates.Mechanisms.step_text(m.catalytic, k) in
                      ("E + A <--> EA", "EQ + P <--> EPQ")
                      for k in eachindex(m.catalytic.steps)])
    @test EnzymeRates.Derivation._term_count(m.catalytic, mask) == 0
    msg = "internal error: the steady-state steps of this mechanism " *
          "do not join its rapid-equilibrium segments"
    @test_throws ErrorException(msg) EnzymeRates.Derivation._law(
        m, EnzymeRates.Derivation._solve(m), mask)
end

@testset "_law refuses a law over the term cap" begin
    # The mechanism of "Rate equation too large error" (test_rate_equation_derivation.jl):
    # 11 forms, 16 steps, V×τ about 29,000.
    m_manual = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: R1
        steps: begin
            E + A <--> EA
            EA <--> FP
            FP <--> F + P
            F + B <--> FB
            FB <--> EQ
            EQ <--> E + Q
            E + R1 <--> ER1
            EA + R1 <--> EAR1
            FP + R1 <--> FPR1
            F + R1 <--> FR1
            FB + R1 <--> FBR1
            ER1 + A <--> EAR1
            EAR1 <--> FPR1
            FPR1 <--> FR1 + P
            FR1 + B <--> FBR1
            FBR1 <--> ER1 + Q
        end
    end
    @test EnzymeRates.Derivation.rate_equation_term_count(m_manual) >
          EnzymeRates.Derivation.MAX_RATE_EQUATION_TERMS
    err = @test_throws(EnzymeRates.Derivation.RateEquationTooLargeError,
                       EnzymeRates.Derivation._memo_law(m_manual))
    @test sprint(showerror, err.value) ==
          "Rate equation for this mechanism has more than 5000 polynomial terms (limit: " *
          "5000). Equations this large take a very long time to compile and are " *
          "unlikely to be practically useful for parameter fitting."
end

@testset "an Int64 overflow in a law is an internal error naming the mechanism" begin
    # Under the term cap no coefficient or exponent leaves Int64. A doctored solve whose
    # first dependent reads one constant twice, each time to the power typemax(Int),
    # forces an exponent past it.
    function _testhelper_overflowing(s)
        d = first(s.dependents)
        f = first(d.factors).first
        big_power = EnzymeRates.Derivation._Dependent(
            d.col, d.keq, [f => typemax(Int) // 1, f => typemax(Int) // 1], d.label)
        EnzymeRates.Derivation._Solve(s.columns, s.names, s.column, s.independents,
                                      [big_power; s.dependents[2:end]], s.info, s.fitted)
    end
    mm = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            ES <--> E + P
        end
    end
    bad = _testhelper_overflowing(EnzymeRates.Derivation._solve(mm))
    err = @test_throws ErrorException EnzymeRates.Derivation._law(mm, bad)
    @test startswith(err.value.msg, "internal error: the rate law of this mechanism " *
                                    "leaves Int64 arithmetic")
    @test endswith(err.value.msg, sprint(show, mm))
    # An MWC law's message names the allosteric mechanism, not its catalytic graph.
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES       :: EqualAI
            ES <--> EP     :: NonequalAI
            E + P ⇌ EP       :: EqualAI
        end
    end
    bad = _testhelper_overflowing(EnzymeRates.Derivation._solve(am))
    err = @test_throws ErrorException EnzymeRates.Derivation._law(am, bad)
    @test startswith(err.value.msg, "internal error: the rate law of this mechanism " *
                                    "leaves Int64 arithmetic")
    @test endswith(err.value.msg, sprint(show, am))
end

# ── The MWC law ─────────────────────────────────────────────────────────────

@testset "the MWC law equals the one-subunit oracle on every allosteric spec" begin
    # At one catalytic subunit the MWC law is the exact steady state of one protomer that
    # flips at free E: 3 points per spec, 1e-20 relative in 256-bit BigFloat. Each
    # conformation's Q has positive coefficients, the law vanishes at equilibrium in each
    # conformation and combined, and the fitted names count the declared dimension. A
    # second derivation gives an equal law with an equal hash.
    specs = [s for s in MECHANISM_TEST_SPECS
             if _testhelper_is_allosteric(s.mechanism)]
    @test length(specs) == 11
    for (i, spec) in enumerate(specs)
        @testset "$(spec.name)" begin
            am = _testhelper_single_subunit(spec.mechanism)
            s = EnzymeRates.Derivation._solve(am)
            A = EnzymeRates.Derivation._law(am, s)
            rng = Random.MersenneTwister(i)
            for _ in 1:3
                point, concs, params = _testhelper_law_inputs(am, s.fitted, rng)
                oracle = _testhelper_steady_state_rate(am, point, concs)
                law = _testhelper_evaluate_law(A, merge(params, concs))
                @test abs(law - oracle) <= big(10)^-20 * abs(oracle)
            end
            @test all(>(0), values(A.active.Q.terms))
            @test all(>(0), values(A.inactive.Q.terms))
            @test all(<(big(10)^-30),
                      _testhelper_detailed_balance(am, Random.MersenneTwister(100 + i)))
            @test length(s.fitted) == _testhelper_declared_dimension(am)
            again = EnzymeRates.Derivation._law(am, s)
            @test again == A && hash(again) == hash(A)
        end
    end
end

@testset "the MWC law equals the n-protomer oracle at 1–3 protomers" begin
    # The bi-uni of the :NonequalAI ground-truth gate: A binds at steady state, B and P at
    # rapid equilibrium, and the chemistry differs between the conformations, so the
    # inactive numerator is live. The oracle solves the concerted network of n protomers,
    # each moving through its conformation's forms, all flipping together at free E.
    mechs = (
        @enzyme_mechanism(begin
            substrates: A, B
            products: P
            oligomeric_state: 1
            steps: begin
                E + A <--> EA        :: EqualAI
                EA + B ⇌ EAB     :: EqualAI
                EAB <--> EP      :: NonequalAI
                E + P ⇌ EP           :: EqualAI
            end
        end),
        @enzyme_mechanism(begin
            substrates: A, B
            products: P
            oligomeric_state: 2
            steps: begin
                E + A <--> EA        :: EqualAI
                EA + B ⇌ EAB     :: EqualAI
                EAB <--> EP      :: NonequalAI
                E + P ⇌ EP           :: EqualAI
            end
        end),
        @enzyme_mechanism(begin
            substrates: A, B
            products: P
            oligomeric_state: 3
            steps: begin
                E + A <--> EA        :: EqualAI
                EA + B ⇌ EAB     :: EqualAI
                EAB <--> EP      :: NonequalAI
                E + P ⇌ EP           :: EqualAI
            end
        end))
    for (n, am) in enumerate(mechs)
        s = EnzymeRates.Derivation._solve(am)
        A = EnzymeRates.Derivation._law(am, s)
        @test A.n == n
        rng = Random.MersenneTwister(30 + n)
        for _ in 1:2
            point, concs, params = _testhelper_law_inputs(am, s.fitted, rng)
            oracle = setprecision(() -> _testhelper_protomer_rate(am, point, concs; n),
                                  BigFloat, 256)
            law = _testhelper_evaluate_law(A, merge(params, concs))
            @test abs(law - oracle) <= big(10)^-20 * abs(oracle)
        end
    end
end

@testset "an :OnlyI law is the exact limit of its vanishing constants" begin
    # The twelve mechanisms whose :OnlyI groups bind modifiers alone and whose laws the
    # package derives, read without the guard (all but the essential activator, whose
    # active conformation runs no step): each law equals the rapid-equilibrium limit of
    # the oracle's own deletion network (`_testhelper_limit_law`) at 3 points, and the
    # active conformation's denominator holds no constant of an :OnlyI group, since the
    # conformation's steps (`steps_running_in_allosteric_conformations`) drop them.
    for (i, (label, am)) in enumerate(pairs(_testhelper_onlyi_fixtures()))
        label === :essential_activator && continue
        @testset "$label" begin
            sub = _testhelper_single_subunit(am)
            s = EnzymeRates.Derivation._solve(sub)
            A = EnzymeRates.Derivation._law(sub, s)
            rng = Random.MersenneTwister(400 + i)
            for _ in 1:3
                point, concs, params = _testhelper_law_inputs(sub, s.fitted, rng)
                exact = setprecision(BigFloat, 256) do
                    _testhelper_limit_law(sub, point, concs, Set{Symbol}(); rng)
                end
                law = _testhelper_evaluate_law(A, merge(params, concs))
                @test abs(law - exact) <= big(10)^-20 * abs(exact)
            end
            cm = sub.catalytic
            onlyi =
                Set(x for (k, ns) in enumerate(_testhelper_edge_names(sub, :I))
                    if sub.allosteric.allosteric_states[cm.steps[k].kinetic_group] ===
                       :OnlyI
                    for x in ns)
            held = Set(A.vars.names[v] for mo in keys(A.active.D.terms) for (v, _) in mo)
            @test !isempty(onlyi)
            @test isempty(intersect(onlyi, held))
        end
    end
end

@testset "a flip-only inactive conformation is Q = 0, D = 1" begin
    # Both bindings :OnlyA: the inactive conformation's free E binds nothing, so that
    # conformation is free E alone. It turns nothing over and weighs 1, so the law is
    # H·QA/(DA + L), and it adds nothing to the term count.
    inert = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: OnlyA
        end
    end
    @test !any(EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(
        inert).inactive)
    s = EnzymeRates.Derivation._solve(inert)
    A = EnzymeRates.Derivation._law(inert, s)
    @test EnzymeRates.Derivation._flip_only(A.inactive) &&
          !EnzymeRates.Derivation._flip_only(A.active)
    @test iszero(A.inactive.Q) && A.inactive.D == one(EnzymeRates.Derivation._Poly)
    @test A.inactive.trees == [one(EnzymeRates.Derivation._Poly)] &&
          isempty(A.inactive.shift)
    @test A.fA == one(EnzymeRates.Derivation._Poly) &&
          A.fI == one(EnzymeRates.Derivation._Poly)
    @test EnzymeRates.Derivation.rate_equation_term_count(inert) == 1
    rng = Random.MersenneTwister(77)
    for _ in 1:3
        point, concs, params = _testhelper_law_inputs(inert, s.fitted, rng)
        oracle = _testhelper_steady_state_rate(inert, point, concs)
        law = _testhelper_evaluate_law(A, merge(params, concs))
        @test abs(law - oracle) <= big(10)^-20 * abs(oracle)
    end
end

@testset "_cross_weights divides out the common integer and monomial content" begin
    vars = EnzymeRates.Derivation._Vars()
    B = EnzymeRates.Derivation._var!(vars, :B; conc = true)
    k = EnzymeRates.Derivation._var!(vars, :k_ES_to_EP)
    q = EnzymeRates.Derivation._var!(vars, :K_EQ_to_E_Q)
    t(c, ps...) = EnzymeRates.Derivation._term(EnzymeRates.Derivation._monomial(ps...), c)
    # Two monomials, 2·k·B and 4·B², share the content 2·B.
    @test EnzymeRates.Derivation._cross_weights(t(2, k => 1, B => 1), t(4, B => 2)) ==
          (t(1, k => 1), t(2, B => 1))
    # Equal weights both become 1, whatever their polynomial part.
    w = t(1, k => 1) + t(1, B => 1)
    @test EnzymeRates.Derivation._cross_weights(w, w) ==
          (one(EnzymeRates.Derivation._Poly), one(EnzymeRates.Derivation._Poly))
    # A two-term weight keeps its primitive part; only the shared B comes out.
    @test EnzymeRates.Derivation._cross_weights(t(1, q => 1, B => 1) * w, t(1, B => 1)) ==
          (t(1, q => 1) * w, one(EnzymeRates.Derivation._Poly))
end

@testset "the cap counts both conformations" begin
    # A random-order bi-bi with the abortive complexes EAQ and EBP, every step
    # steady state: 9 forms and 480 spanning trees, 4,320 products. With its chemistry
    # :NonequalAI at one catalytic subunit both conformations keep every step, so the law
    # counts 4,320 + 4,320 = 8,640: over the cap of 5,000, though each conformation is
    # under it. The engine refuses it before any polynomial product, with
    # RateEquationTooLargeError's message.
    am = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA          :: EqualAI
            E + B <--> EB          :: EqualAI
            EA + B <--> EAB    :: EqualAI
            EB + A <--> EAB    :: EqualAI
            EAB <--> EPQ     :: NonequalAI
            EPQ <--> EP + Q    :: EqualAI
            EPQ <--> EQ + P    :: EqualAI
            EP <--> E + P          :: EqualAI
            EQ <--> E + Q          :: EqualAI
            EA + Q <--> EAQ    :: EqualAI
            EQ + A <--> EAQ    :: EqualAI
            EB + P <--> EBP    :: EqualAI
            EP + B <--> EBP    :: EqualAI
        end
    end
    cm = am.catalytic
    runs = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(am)
    @test _testhelper_laplacian_tree_count(_testhelper_plain_of(am), runs.active) == 4320
    @test _testhelper_laplacian_tree_count(_testhelper_plain_of(am), runs.inactive) == 4320
    @test EnzymeRates.Derivation._term_count(cm, runs.active) <=
          EnzymeRates.Derivation.MAX_RATE_EQUATION_TERMS
    @test EnzymeRates.Derivation.rate_equation_term_count(am) == 8640
    @test_throws(EnzymeRates.Derivation.RateEquationTooLargeError,
                 EnzymeRates.Derivation._memo_law(am))
    @test_throws "has more than 5000 polynomial terms" EnzymeRates.Derivation._memo_law(am)
    @test_throws EnzymeRates.Derivation.RateEquationTooLargeError rate_equation_string(am)
end

# ── MWC laws against mass-action ground truth ─────────────────────────────────────────

# ── The general solver ──────────────────────────────────────────────────────
# Build the pseudo-first-order rate matrix from a directed edge list, impose
# mass conservation Σc = Etot by replacing the first row, solve the linear
# steady state, and return the net catalytic flux Σ(kf·c[reactant] − kr·c[product])
# over the catalytic (product-forming) edges — the reaction velocity.
#
# `species`  :: Vector{Symbol}   (species[1] carries the conservation row)
# `edges`    :: Vector of (from::Symbol, to::Symbol, rate::Float64), directed,
#               rate = pseudo-first-order
# `cat_edges`:: Vector of (reactant::Symbol, product::Symbol, kf::Float64, kr::Float64)
"Net catalytic flux at steady state for an explicit two-conformation network."
function _testhelper_mwc_ground_truth_flux(species, edges, cat_edges, Etot)
    n = length(species); idx = Dict(s => i for (i, s) in enumerate(species))
    M = zeros(n, n)
    for (a, b, r) in edges
        M[idx[b], idx[a]] += r
        M[idx[a], idx[a]] -= r
    end
    A = copy(M); A[1, :] .= 1.0
    rhs = zeros(n); rhs[1] = Etot
    c = A \ rhs
    sum(kf * c[idx[r]] - kr * c[idx[p]] for (r, p, kf, kr) in cat_edges)
end

# ── Uni-uni :OnlyA network (the exemplar the fix must get right) ─────────────
# S binds :OnlyA (only the active conformation) and catalysis is :OnlyA (inactive
# rate k_I = 0), P :EqualAI. Forms: E_A, ES_A, EP_A (active) and E_I, EP_I
# (inactive). The inactive conformation has no catalytic edge — that is what
# :OnlyA catalysis means, and with no inactive S binding it never reaches ES.
# Fast RE bindings and flips use FAST; catalysis is O(1). Detailed-balance flip
# ratio [X_I]/[X_A] = L·∏(K_A_i/K_I_i); every present flip here carries an
# :EqualAI ligand (or none), so the ratio is L.
function _testhelper_uni_onlyA_flux(KA, KP, k; L, Keq, S, P, FAST=1e7)
    kr = k * KP / (Keq * KA)
    species = [:E_A, :ES_A, :EP_A, :E_I, :EP_I]
    edges = [
        (:E_A, :ES_A, FAST * S / KA), (:ES_A, :E_A, FAST),   # active S binding (RE)
        (:E_A, :EP_A, FAST * P / KP), (:EP_A, :E_A, FAST),   # active P binding (RE)
        (:E_I, :EP_I, FAST * P / KP), (:EP_I, :E_I, FAST),   # inactive P binding (EqualAI, RE)
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),          # free-enzyme flip, ratio L
        (:EP_A, :EP_I, FAST * L), (:EP_I, :EP_A, FAST),      # EP flip, ratio L
        (:ES_A, :EP_A, k), (:EP_A, :ES_A, kr),              # active catalysis (SS)
    ]
    cat_edges = [(:ES_A, :EP_A, k, kr)]
    _testhelper_mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# ── Multi-:OnlyA bi-uni network (both substrates bind the active state only) ──
# A + B ⇌ P. Both A and B bind :OnlyA (active conformation only) and catalysis is
# :OnlyA (inactive rate k_I = 0), P :EqualAI. Active forms E_A, EA_A, EAB_A, EP_A.
# Inactive forms E_I, EP_I, RE-connected via P binding — a single segment with no
# catalytic edge (that is what :OnlyA catalysis means; with A,B :OnlyA the inactive
# state never reaches EA or EAB), so free E weighs 1 in the inactive conformation.
# Flip ratio [X_I]/[X_A] = L·∏(K_A_i/K_I_i): free enzyme and EP (:EqualAI or bare)
# flip with ratio L; an :OnlyA-ligand-bearing state has ratio 0 (no flip), so EA
# and EAB never flip.
function _testhelper_multi_onlyA_flux(KA, KB, KP, k; L, Keq, A, B, P, FAST=1e7)
    kr = k * KP / (Keq * KA * KB)
    species = [:E_A, :EA_A, :EAB_A, :EP_A, :E_I, :EP_I]
    edges = [
        (:E_A, :EA_A, FAST * A / KA), (:EA_A, :E_A, FAST),    # active A binding (RE)
        (:EA_A, :EAB_A, FAST * B / KB), (:EAB_A, :EA_A, FAST),# active B binding (RE)
        (:E_A, :EP_A, FAST * P / KP), (:EP_A, :E_A, FAST),    # active P binding (RE)
        (:E_I, :EP_I, FAST * P / KP), (:EP_I, :E_I, FAST),    # inactive P binding (EqualAI, RE)
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),           # free-enzyme flip, ratio L
        (:EP_A, :EP_I, FAST * L), (:EP_I, :EP_A, FAST),       # EP flip, ratio L
        (:EAB_A, :EP_A, k), (:EP_A, :EAB_A, kr),             # active catalysis (SS)
    ]
    cat_edges = [(:EAB_A, :EP_A, k, kr)]
    _testhelper_mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# ── The gate: :OnlyA MWC derivation matches mass-action ground truth ─────────
# S binds :OnlyA and catalysis is :OnlyA (k_I = 0), so the inactive conformation
# runs no catalysis: each conformation's law is one rapid-equilibrium segment and free E
# weighs 1 in both (fA = fI = 1), so the law is H·Q_A/(D_A + L·D_I). The gate checks the
# derived rate against the independent mass-action ground truth.
@testset "OnlyA MWC derivation matches mass-action ground truth" begin
    onlyA = @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES     :: OnlyA
            ES <--> EP   :: OnlyA
            E + P ⇌ EP     :: EqualAI
        end
    end
    # (:K_A_ES_to_E_S, :K_EP_to_E_P, :k_A_ES_to_EP, :L)
    fp = parameters(onlyA)
    onlyA_rate = ER.compile_rate_equation(onlyA)
    @test fp == (:K_A_ES_to_E_S, :K_EP_to_E_P, :k_A_ES_to_EP, :L)

    rng = MersenneTwister(20260713)
    for _ in 1:5
        KA = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng); k = 0.5 + 2rand(rng)
        L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        S = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map parameters(m) -> ground-truth params:
        #   K_A_ES_to_E_S=KA, K_EP_to_E_P=KP, k_A_ES_to_EP=k.
        prm = NamedTuple{fp}((KA, KP, k, L))
        v_code = real(onlyA_rate((S=S, P=P), prm, Keq))
        v_gt = _testhelper_uni_onlyA_flux(KA, KP, k, L=L, Keq=Keq, S=S, P=P)
        @test isapprox(v_code, v_gt; rtol=1e-4)
        # Far from equilibrium the law runs forward.
        @test real(onlyA_rate((S=10.0, P=1e-3), prm, Keq)) > 0
        # self-validation: L = 0 → inactive conformation unpopulated → the
        # single-conformation (non-allosteric) uni-uni rate.
        kr = k * KP / (Keq * KA)
        nonallo = (k * S / KA - kr * P / KP) / (1 + S / KA + P / KP)
        @test isapprox(
            _testhelper_uni_onlyA_flux(KA, KP, k, L=0.0, Keq=Keq, S=S, P=P), nonallo;
            rtol=1e-4)
    end

    # The inactive state neither binds S nor catalyzes, so its free-enzyme graph is a
    # single segment and both conformations weigh free E 1. The cross-weight regime, where
    # free E's weight differs between the conformations (fA ≠ fI), needs a steady-state
    # :OnlyA binding and is covered by the metabolite-bearing-D :OnlyA gate below.
    A = ER.Derivation._memo_law(onlyA)
    # single active segment
    @test A.active.trees[1] * A.active.shift == one(ER.Derivation._Poly)
    # one inactive segment
    @test A.inactive.trees[1] * A.inactive.shift == one(ER.Derivation._Poly)
    @test A.fA == A.fI == one(ER.Derivation._Poly)

    # Scaling the rate constants to a measured turnover stays consistent with the law.
    prm = NamedTuple{fp}((1.3, 0.9, 2.1, 0.7))
    Keq = 3.0
    v0 = onlyA_rate((S = 2.0, P = 0.5), prm, Keq)
    rescaled = EnzymeRates.Fitting._scale_rate_constants(onlyA, prm, 1.5 / v0)
    @test isapprox(real(onlyA_rate((S = 2.0, P = 0.5), rescaled, Keq)), 1.5;
                   rtol=1e-9)
end

# ── The gate: multi-:OnlyA bi-uni derivation matches mass-action ground truth ─
# Both A and B bind :OnlyA and catalysis is :OnlyA (k_I = 0), so the inactive
# conformation runs no catalysis: each conformation's law is one rapid-equilibrium
# segment and free E weighs 1 in both (fA = fI = 1), so the law is H·Q_A/(D_A + L·D_I).
# The gate checks the derived rate against the independent mass-action ground truth.
@testset "multi-OnlyA MWC derivation matches mass-action ground truth" begin
    multiA = @enzyme_mechanism begin
        substrates: A, B ; products: P ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA          :: OnlyA
            EA + B ⇌ EAB    :: OnlyA
            EAB <--> EP     :: OnlyA
            E + P ⇌ EP          :: EqualAI
        end
    end
    # (:K_A_EAB_to_EA_B, :K_A_EA_to_E_A, :K_EP_to_E_P, :k_A_EAB_to_EP, :L)
    fp = parameters(multiA)
    multiA_rate = ER.compile_rate_equation(multiA)
    @test fp == (:K_A_EAB_to_EA_B, :K_A_EA_to_E_A, :K_EP_to_E_P, :k_A_EAB_to_EP, :L)

    rng = MersenneTwister(20260713)
    for _ in 1:5
        KA = 0.5 + 2rand(rng); KB = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng)
        k = 0.5 + 2rand(rng); L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        A = 0.5 + 2rand(rng); B = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map parameters(m) -> ground-truth params:
        #   K_A_EA_to_E_A=KA, K_EP_to_E_P=KP, K_A_EAB_to_EA_B=KB, k_A_EAB_to_EP=k.
        prm = NamedTuple{fp}((KB, KA, KP, k, L))
        v_code = real(multiA_rate((A=A, B=B, P=P), prm, Keq))
        v_gt = _testhelper_multi_onlyA_flux(KA, KB, KP, k, L=L, Keq=Keq, A=A, B=B, P=P)
        @test isapprox(v_code, v_gt; rtol=1e-4)
        @test real(multiA_rate((A=10.0, B=10.0, P=1e-3), prm, Keq)) > 0
        # self-validation: L = 0 → inactive conformation unpopulated → the
        # single-conformation (non-allosteric) bi-uni rate.
        kr = k * KP / (Keq * KA * KB)
        nonallo = (k*A*B/(KA*KB) - kr*P/KP) / (1 + A/KA + A*B/(KA*KB) + P/KP)
        @test isapprox(
            _testhelper_multi_onlyA_flux(KA, KB, KP, k, L=0.0, Keq=Keq, A=A, B=B, P=P),
            nonallo; rtol=1e-4)
    end
end

# ── Metabolite-bearing-D ordered bi-uni :OnlyA network (the LDH i-state case) ──
# S + B ⇌ P. S binds :OnlyA via a STEADY-STATE step (E + S <--> ES); B binds
# :EqualAI at rapid equilibrium on ES; catalysis EBS <--> EP is
# steady-state :OnlyA (inactive rate k_I = 0); P binds :EqualAI at rapid
# equilibrium. Because B binds rapidly on the steady-state catalytic path, free E's
# weight in the active conformation CARRIES the metabolite B: koff + k·B/K_B. The
# inactive conformation runs no catalysis and never binds S or B (both :OnlyA), so its
# graph is a single rapid-equilibrium segment {E_I, EP_I} and free E weighs 1 there.
# Free E's weight differs between the conformations (fA ≠ fI), a metabolite-bearing
# active weight against a bare inactive one: the cross-weight regime, which has no
# other ground truth. Flip ratio [X_I]/[X_A] = L·∏(K_A/K_I): free enzyme and
# EP flip (ratio L); an S-bearing state has ratio 0 (S :OnlyA), so ES and
# EBS never flip. Reverse catalysis kr from the Haldane relation.
function _testhelper_metab_dfree_onlyA_flux(kon, koff, KB, KP, k; L, Keq, S, B, P, FAST=1e7)
    kr = k * kon * KP / (koff * KB * Keq)
    species = [:E_A, :ES_A, :ESB_A, :EP_A, :E_I, :EP_I]
    edges = [
        (:E_A, :ES_A, kon * S), (:ES_A, :E_A, koff),          # S binding (SS, active only)
        (:ES_A, :ESB_A, FAST * B / KB), (:ESB_A, :ES_A, FAST),# B binding (RE, active)
        (:E_A, :EP_A, FAST * P / KP), (:EP_A, :E_A, FAST),    # P binding (RE, active)
        (:E_I, :EP_I, FAST * P / KP), (:EP_I, :E_I, FAST),    # P binding (RE, inactive)
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),           # free-enzyme flip, ratio L
        (:EP_A, :EP_I, FAST * L), (:EP_I, :EP_A, FAST),       # EP flip, ratio L
        (:ESB_A, :EP_A, k), (:EP_A, :ESB_A, kr),            # active catalysis (SS)
    ]
    cat_edges = [(:ESB_A, :EP_A, k, kr)]
    _testhelper_mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# ── Single-conformation (non-allosteric) reference for the ordered bi-uni ─────
# The active mechanism alone (E, ES, EBS, EP); no inactive conformation,
# no flips. The L → 0 limit of the :OnlyA flux must equal this rate.
function _testhelper_metab_dfree_base_flux(kon, koff, KB, KP, k, Keq, S, B, P; FAST=1e7)
    kr = k * kon * KP / (koff * KB * Keq)
    species = [:E, :ES, :ESB, :EP]
    edges = [
        (:E, :ES, kon * S), (:ES, :E, koff),
        (:ES, :ESB, FAST * B / KB), (:ESB, :ES, FAST),
        (:E, :EP, FAST * P / KP), (:EP, :E, FAST),
        (:ESB, :EP, k), (:EP, :ESB, kr),
    ]
    _testhelper_mwc_ground_truth_flux(species, edges, [(:ESB, :EP, k, kr)], 1.0)
end

# ── The gate: metabolite-bearing-D :OnlyA derivation matches mass-action GT ────
# The LDH i-state case. S binds :OnlyA on the steady-state catalytic path
# and B binds rapidly there too, so free E's weight carries the metabolite B and differs
# between the conformations (fA ≠ fI). The cross weights fA, fI supply the common
# free-enzyme basis this regime needs; dropping them gets it wrong. Regression guard.
@testset "metabolite-bearing-D MWC derivation matches mass-action ground truth" begin
    metabD = @enzyme_mechanism begin
        substrates: S, B ; products: P ; oligomeric_state: 1
        steps: begin
            E + S <--> ES      :: OnlyA
            ES + B ⇌ EBS   :: EqualAI
            EBS <--> EP    :: OnlyA
            E + P ⇌ EP         :: EqualAI
        end
    end
    # (:K_EBS_to_ES_B,:K_EP_to_E_P,:k_A_EBS_to_EP,:k_A_ES_to_E_S,:k_A_E_S_to_ES,:L)
    fp = parameters(metabD)
    metabD_rate = ER.compile_rate_equation(metabD)
    @test fp == (:K_EBS_to_ES_B, :K_EP_to_E_P, :k_A_EBS_to_EP, :k_A_ES_to_E_S,
                 :k_A_E_S_to_ES, :L)

    rng = MersenneTwister(20260713)
    for _ in 1:5
        kon = 0.5 + 2rand(rng); koff = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng)
        KB = 0.5 + 2rand(rng); k = 0.5 + 2rand(rng); L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        S = 0.5 + 2rand(rng); B = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map parameters(m) -> ground-truth params:
        #   k_A_E_S_to_ES=kon, k_A_ES_to_E_S=koff, K_EBS_to_ES_B=KB, K_EP_to_E_P=KP,
        #   k_A_EBS_to_EP=k.
        prm = NamedTuple{fp}((KB, KP, k, koff, kon, L))
        v_code = real(metabD_rate((S=S, B=B, P=P), prm, Keq))
        v_gt = _testhelper_metab_dfree_onlyA_flux(
            kon, koff, KB, KP, k, L=L, Keq=Keq, S=S, B=B, P=P)
        @test isapprox(v_code, v_gt; rtol=1e-4)
        @test real(metabD_rate((S=10.0, B=10.0, P=1e-3), prm, Keq)) > 0
        # self-validation: L = 0 → inactive conformation unpopulated → the
        # single-conformation (non-allosteric) ordered bi-uni rate.
        @test isapprox(
            _testhelper_metab_dfree_onlyA_flux(
                kon, koff, KB, KP, k, L=0.0, Keq=Keq, S=S, B=B, P=P),
            _testhelper_metab_dfree_base_flux(
                kon, koff, KB, KP, k, Keq, S, B, P); rtol=1e-4)
    end

    # B = 0 zeroes the metabolite B in free E's weight. The active steady-state S
    # binding leaves a bare koff term in free E's weight, so a reverse path from
    # product back to free enzyme stays open and the reverse flux survives.
    kon, koff, KP, KB, k, L, Keq, S, P = 1.7, 1.1, 0.9, 0.8, 2.1, 0.7, 3.0, 1.1, 0.6
    prm = NamedTuple{fp}((KB, KP, k, koff, kon, L))
    v0 = real(metabD_rate((S=S, B=0.0, P=P), prm, Keq))
    @test isapprox(v0,
        _testhelper_metab_dfree_onlyA_flux(
            kon, koff, KB, KP, k; L=L, Keq=Keq, S=S, B=0.0, P=P);
        rtol=1e-4)
    @test v0 < 0    # reverse flux survives at B=0 (invalid :EqualAI-cat trap→0 was fake)
    @test real(metabD_rate((S=10.0, B=10.0, P=1e-3), prm, Keq)) > 0
end

# ── Free-flip-only ordered bi-uni :NonequalAI-catalysis reference ────────────
# A + B ⇌ P. A binds :EqualAI via a STEADY-STATE step (E + A <--> EA, kon·A /
# koff); B binds :EqualAI at rapid equilibrium on EA; catalysis EAB <--> EP
# is STEADY-STATE :NonequalAI (rate k_A in the active conformation, k_I in the
# inactive); P binds :EqualAI at rapid equilibrium. Every ligand is :EqualAI
# (K_I = K_A).
#
# Formulation-1 reference: only the free enzyme flips conformation. Each
# conformation runs its own catalytic cycle, coupled only through the shared
# free-enzyme pool — the model this package derives (commit-when-free). Reverse
# catalysis kr from the Haldane relation. `_testhelper_biuni_mwc_oligomer_flux(1, …;
# freeflip=false)` is the per-form-flip model with the same forms and edges, which
# also flips every catalytic intermediate.
function _testhelper_biuni_nonequalAI_freeflip_flux(
        kon, koff, KB, KP; k_A, k_I, L, Keq, A, B, P, FAST=1e7)
    krA = k_A * kon * KP / (koff * KB * Keq); krI = k_I * kon * KP / (koff * KB * Keq)
    species = [:E_A, :EA_A, :EAB_A, :EP_A, :E_I, :EA_I, :EAB_I, :EP_I]
    edges = [
        (:E_A, :EA_A, kon*A), (:EA_A, :E_A, koff), (:E_I, :EA_I, kon*A), (:EA_I, :E_I, koff),
        (:EA_A, :EAB_A, FAST*B/KB), (:EAB_A, :EA_A, FAST),
        (:EA_I, :EAB_I, FAST*B/KB), (:EAB_I, :EA_I, FAST),
        (:E_A, :EP_A, FAST*P/KP), (:EP_A, :E_A, FAST),
        (:E_I, :EP_I, FAST*P/KP), (:EP_I, :E_I, FAST),
        (:E_A, :E_I, FAST*L), (:E_I, :E_A, FAST),                  # only free enzyme flips
        (:EAB_A, :EP_A, k_A), (:EP_A, :EAB_A, krA),
        (:EAB_I, :EP_I, k_I), (:EP_I, :EAB_I, krI),
    ]
    cat_edges = [(:EAB_A, :EP_A, k_A, krA), (:EAB_I, :EP_I, k_I, krI)]
    _testhelper_mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# ── :NonequalAI-catalysis harness self-validation ────────────────────────────
# First ground the hand-written mass-action reference `_testhelper_metab_dfree_base_flux`
# (the single-conformation ordered bi-uni, reused as the non-allosteric rate at
# k = k_A) against the ODE-validated non-allosteric rate law — an
# independent, separately ODE-cross-checked code path. This closes the trust loop
# for the whole `_testhelper_mwc_ground_truth_flux` harness. Then confirm the free-flip-only
# reference degenerates correctly: (a) L = 0 (inactive unpopulated) → the base rate
# at k = k_A; (b) k_I = k_A (conformations identical) → the base rate, independent
# of L.
@testset ":NonequalAI-catalysis ground-truth harness self-validation" begin
    base = @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A <--> EA
            EA + B ⇌ EAB
            EAB <--> EP
            E + P ⇌ EP
        end
    end
    bfp = parameters(base)
    base_rate_law = ER.compile_rate_equation(base)

    rng = MersenneTwister(20260713)
    for _ in 1:5
        kon = 0.5 + 2rand(rng); koff = 0.5 + 2rand(rng)
        KB = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng)
        kA = 0.5 + 2rand(rng); kI = 0.5 + 2rand(rng); Keq = 2.0 + 2rand(rng)
        A = 0.5 + 2rand(rng); B = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        L = 0.5 + rand(rng)

        base_rate = _testhelper_metab_dfree_base_flux(kon, koff, KB, KP, kA, Keq, A, B, P)

        # `_testhelper_metab_dfree_base_flux` vs the ODE-validated non-allosteric
        # rate law.
        bd = Dict(:k_E_A_to_EA=>kon, :k_EA_to_E_A=>koff, :K_EAB_to_EA_B=>KB,
                  :K_EP_to_E_P=>KP, :k_EAB_to_EP=>kA)
        bprm = NamedTuple{bfp}(Tuple(bd[s] for s in bfp))
        @test isapprox(base_rate,
            real(base_rate_law((A=A, B=B, P=P), bprm, Keq)); rtol=1e-4)

        # (a) L = 0 : inactive conformation unpopulated → base rate at k_A.
        f0 = _testhelper_biuni_nonequalAI_freeflip_flux(kon, koff, KB, KP;
            k_A=kA, k_I=kI, L=0.0, Keq=Keq, A=A, B=B, P=P)
        @test isapprox(f0, base_rate; rtol=1e-4)

        # (b) k_I = k_A : conformations identical → base rate, independent of L.
        fe = _testhelper_biuni_nonequalAI_freeflip_flux(kon, koff, KB, KP;
            k_A=kA, k_I=kA, L=L, Keq=Keq, A=A, B=B, P=P)
        fe5 = _testhelper_biuni_nonequalAI_freeflip_flux(kon, koff, KB, KP;
            k_A=kA, k_I=kA, L=5.0, Keq=Keq, A=A, B=B, P=P)
        @test isapprox(fe, base_rate; rtol=1e-4)
        @test isapprox(fe, fe5; rtol=1e-4)
    end
end

# ── The gate: allosteric ping-pong (single free form, formulation-1 normalization) ─
# A ping-pong mechanism has ONE free enzyme form — the covalent intermediate
# `E_res_plusA_minusP` carries a residual, so it is not free — hence free E's
# weight is well-defined and the derivation normalizes it like any other mechanism (no
# special case, no guard). This gate confirms the derivation handles the
# four-segment residual-bearing King–Altman and that the normalization collapses
# correctly for identical conformations: an all-:EqualAI ping-pong must equal the
# non-allosteric ping-pong and be independent of L.
@testset "allosteric ping-pong self-consistency" begin
    nonallo = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q
        steps: begin
            E + A <--> EA
            EA <--> E_res_plusA_minusP + P
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP
            EB_res_plusA_minusP <--> E + Q
        end
    end
    alloEq = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: EqualAI
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: EqualAI
        end
    end
    fpn = parameters(nonallo)
    fpa = parameters(alloEq)
    nonallo_rate = ER.compile_rate_equation(nonallo)
    alloEq_rate = ER.compile_rate_equation(alloEq)
    rng = MersenneTwister(20260714)
    for _ in 1:5
        vals = Dict(s => 0.5 + 2rand(rng) for s in fpn)
        Keq = 2.0 + 2rand(rng)
        concs = (A=0.5+2rand(rng), B=0.5+2rand(rng), P=0.5+2rand(rng), Q=0.5+2rand(rng))
        pn = NamedTuple{fpn}(Tuple(vals[s] for s in fpn))
        vn = real(nonallo_rate(concs, pn, Keq))
        # all-:EqualAI is L-independent and equals non-allosteric
        for L in (0.7, 5.0)
            pa = NamedTuple{fpa}(Tuple(s === :L ? L : vals[s] for s in fpa))
            @test isapprox(real(alloEq_rate(concs, pa, Keq)), vn; rtol=1e-6)
        end
    end
end

# ── Two-conformation ping-pong bi-bi network (formulation 1) ─────────────────
# Ping-pong has TWO empty-bound forms — free E and the covalent intermediate F
# (`E_res_plusA_minusP`). Only free E flips; F does not. The E_A<->E_I edge is
# therefore the ONLY cut between the two conformation subnetworks, so at steady
# state it carries zero net flux and E_I/E_A sits at exactly L for ANY FAST — the
# fast-flip limit is exact for this topology rather than a large-FAST limit. Each
# conformation turns its own four-form cycle (E, EA, F, FB), coupled only through
# the shared free-enzyme pool.
#
# Thermodynamics pins ONE combination per conformation, not one per half-reaction.
# Detailed balance around the closed cycle E → EA → F → FB → E requires
#   (k1f·A/k1r)·(k2f/(k2r·P))·(k3f·B/k3r)·(k4f/(k4r·Q)) = 1
# at the equilibrium ratio P·Q/(A·B) = Keq, i.e. the overall Haldane relation
#   k1f·k2f·k3f·k4f = Keq · k1r·k2r·k3r·k4r,
# matching the numerator of Segel Eq. IX-140 (k1f·k2f·k3f·k4f·A·B −
# k1r·k2r·k3r·k4r·P·Q). The two half-reactions' equilibrium constants are NOT
# separately fixed — only their product is — so exactly one reverse constant per
# conformation is dependent. `k1r` (shared, :EqualAI) closes the active cycle;
# `k2r_I` then closes the inactive one against that same `k1r`.
function _testhelper_pingpong_nonequalAI_freeflip_flux(k1f, k3f, k3r;
        k2f_A, k2r_A, k4f_A, k4r_A, k2f_I, k4f_I, k4r_I,
        L, Keq, A, B, P, Q, FAST=1e7)
    k1r   = k1f * k2f_A * k3f * k4f_A / (Keq * k2r_A * k3r * k4r_A)
    k2r_I = k1f * k2f_I * k3f * k4f_I / (Keq * k1r * k3r * k4r_I)
    species = [:E_A, :EA_A, :F_A, :FB_A, :E_I, :EA_I, :F_I, :FB_I]
    edges = Tuple{Symbol,Symbol,Float64}[]
    cat_edges = Tuple{Symbol,Symbol,Float64,Float64}[]
    for (c, k2f, k2r, k4f, k4r) in ((:A, k2f_A, k2r_A, k4f_A, k4r_A),
                                    (:I, k2f_I, k2r_I, k4f_I, k4r_I))
        e, ea = Symbol(:E_, c), Symbol(:EA_, c)
        f, fb = Symbol(:F_, c), Symbol(:FB_, c)
        append!(edges, [
            (e, ea, k1f * A), (ea, e, k1r),               # E + A ⇌ EA
            (ea, f, k2f), (f, ea, k2r * P),               # EA ⇌ F + P
            (f, fb, k3f * B), (fb, f, k3r),               # F + B ⇌ FB
            (fb, e, k4f), (e, fb, k4r * Q),               # FB ⇌ E + Q
        ])
        push!(cat_edges, (ea, f, k2f, k2r * P))           # net flux across the P cut
    end
    push!(edges, (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST))   # only free enzyme flips
    _testhelper_mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# `_testhelper_rate_ping_pong_bi_bi` in
# `test/testhelpers.jl` transcribes this same Segel
# formula, and both transcriptions are live. Keep them independent rather than sharing one:
# a shared transcription error would green this gate and that one at once, whereas two
# independent transcriptions cross-check each other. Sharing would also couple this gate to
# the MECHANISM_TEST_SPECS fixture, whose copy takes
# `(params::NamedTuple, concs::NamedTuple)` with an `Etotal` rather than the 12 positional
# scalars this one takes.
"Segel Eq. IX-140 ping-pong bi-bi rate: E + A ⇌ EA ⇌ F + P; F + B ⇌ FB ⇌ E + Q."
function _testhelper_segel_pingpong_flux(k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, A, B, P, Q)
    num = k1f*k2f*k3f*k4f*A*B - k1r*k2r*k3r*k4r*P*Q
    den = k1f*k2f*(k3r+k4f)*A + k3f*k4f*(k1r+k2f)*B +
          k1r*k2r*(k3r+k4f)*P + k3r*k4r*(k1r+k2f)*Q +
          k1f*k3f*(k2f+k4f)*A*B + k1f*k2r*(k3r+k4f)*A*P +
          k3f*k4r*(k1r+k2f)*B*Q + k2r*k4r*(k1r+k3r)*P*Q
    num / den
end

# ── Ping-pong ground-truth harness self-validation ───────────────────────────
# (d) is the load-bearing one: FAST-invariance is exact physics for this topology
# (the free-enzyme flip is the only cut, so it carries zero net flux), and it is
# what licenses the tight gate below. (a) anchors the network against an
# independent closed form rather than another network solve.
@testset "ping-pong free-flip ground-truth harness self-validation" begin
    rng = MersenneTwister(20260716)
    for _ in 1:4
        k1f = 0.5+2rand(rng); k3f = 0.5+2rand(rng); k3r = 0.5+2rand(rng)
        k2f_A = 0.5+2rand(rng); k2r_A = 0.5+2rand(rng)
        k4f_A = 0.5+2rand(rng); k4r_A = 0.5+2rand(rng)
        k2f_I = 0.5+2rand(rng); k4f_I = 0.5+2rand(rng); k4r_I = 0.5+2rand(rng)
        L = 0.5+rand(rng); Keq = 2.0+2rand(rng)
        A = 0.5+2rand(rng); B = 0.5+2rand(rng); P = 0.5+2rand(rng); Q = 0.5+2rand(rng)
        act = (k2f_A=k2f_A, k2r_A=k2r_A, k4f_A=k4f_A, k4r_A=k4r_A)
        ina = (k2f_I=k2f_I, k4f_I=k4f_I, k4r_I=k4r_I)

        # (a) L = 0 : inactive unpopulated → the active-only ping-pong, checked
        #     against the Segel closed form. `k1r` is the dependent reverse
        #     constant the Haldane fixes; Segel's k1r is that same constant.
        k1r = k1f*k2f_A*k3f*k4f_A / (Keq*k2r_A*k3r*k4r_A)
        f0 = _testhelper_pingpong_nonequalAI_freeflip_flux(k1f, k3f, k3r; act..., ina...,
            L=0.0, Keq=Keq, A=A, B=B, P=P, Q=Q)
        @test isapprox(f0, _testhelper_segel_pingpong_flux(k1f, k1r, k2f_A, k2r_A,
            k3f, k3r, k4f_A, k4r_A, A, B, P, Q); rtol=1e-9)

        # (b) identical conformations → the active-only rate, independent of L.
        same = (k2f_I=k2f_A, k4f_I=k4f_A, k4r_I=k4r_A)
        fe = _testhelper_pingpong_nonequalAI_freeflip_flux(k1f, k3f, k3r; act..., same...,
            L=L, Keq=Keq, A=A, B=B, P=P, Q=Q)
        fe5 = _testhelper_pingpong_nonequalAI_freeflip_flux(k1f, k3f, k3r; act..., same...,
            L=5.0, Keq=Keq, A=A, B=B, P=P, Q=Q)
        @test isapprox(fe, f0; rtol=1e-9)
        @test isapprox(fe, fe5; rtol=1e-9)

        # (c) v = 0 at the equilibrium metabolite ratio P·Q/(A·B) = Keq. Both
        #     conformations are live and unequal, so this pins both Haldanes.
        Qeq = Keq * A * B / P
        @test abs(_testhelper_pingpong_nonequalAI_freeflip_flux(k1f, k3f, k3r; act...,
            ina..., L=L, Keq=Keq, A=A, B=B, P=P, Q=Qeq)) < 1e-9

        # (d) FAST-invariance: the free-enzyme flip is the only cut between the
        #     conformation subnetworks, so it carries zero net flux and E_I/E_A is
        #     exactly L for any FAST. The fast-flip limit is therefore exact here.
        v_live = _testhelper_pingpong_nonequalAI_freeflip_flux(
            k1f, k3f, k3r; act..., ina..., L=L, Keq=Keq, A=A, B=B, P=P, Q=Q)
        for fast in (1e2, 1e12)
            @test isapprox(v_live, _testhelper_pingpong_nonequalAI_freeflip_flux(k1f, k3f,
                k3r; act..., ina..., L=L, Keq=Keq, A=A, B=B, P=P, Q=Q, FAST=fast);
                rtol=1e-9)
        end

        # (e) a live inactive conformation must move the flux, or the gate below
        #     would pass without ever exercising the cross term.
        @test !isapprox(v_live, f0; rtol=1e-3)
    end
end

# ── The gate: :NonequalAI ping-pong matches mass-action ground truth ──────────
# The two empty-bound forms (free E and the covalent F) are what a cross-weighted
# combine must get right: only E flips, so the derivation must normalize the two
# conformations on the free-enzyme segment alone and leave F out of the flip. Both
# catalytic steps are :NonequalAI, so both conformations turn a productive cycle
# with different rate constants (D_A ≠ D_I) and the cross term is live. The gate
# goes red if the derivation flips F, drops the cross weights, or mis-renders the
# normalization. Only the law without cross weights has a measured margin: against a
# per-form-flip (formulation-2) oracle it deviates 0.95%-93% from the derivation,
# orders of magnitude above the 1e-10 tolerance below.
#
# Because the free-enzyme flip is the only cut, the combine is algebraically exact
# here rather than a large-FAST limit, so this gate runs far tighter than the 1e-4
# gates above. Measured: the derivation tracks the oracle to 2.7e-13 worst case
# over 300 random draws, and to 2.2e-15 against the same oracle solved in
# BigFloat — the derivation is exact to a few ulp, and the Float64 residual is the
# oracle's own linear solve, whose FAST-invariance (exact physics) itself only
# holds to 1.3e-13. rtol 1e-10 clears that noise floor with three orders to spare
# and is still seven orders tighter than any real error mode.
@testset ":NonequalAI ping-pong MWC derivation matches mass-action ground truth" begin
    allo = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: NonequalAI
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: NonequalAI
        end
    end
    fp = parameters(allo)
    allo_rate = ER.compile_rate_equation(allo)
    @test fp == (:k_A_EA_to_E_res_plusA_minusP_P, :k_A_EB_res_plusA_minusP_to_E_Q,
                 :k_A_E_Q_to_EB_res_plusA_minusP, :k_A_E_res_plusA_minusP_P_to_EA,
                 :k_EA_to_E_A, :k_E_A_to_EA,
                 :k_E_res_plusA_minusP_B_to_EB_res_plusA_minusP,
                 :k_I_EA_to_E_res_plusA_minusP_P, :k_I_EB_res_plusA_minusP_to_E_Q,
                 :k_I_E_Q_to_EB_res_plusA_minusP, :L)

    rng = MersenneTwister(20260716)
    for _ in 1:6
        k1f = 0.5+2rand(rng); k3f = 0.5+2rand(rng); k3r = 0.5+2rand(rng)
        k2f_A = 0.5+2rand(rng); k2r_A = 0.5+2rand(rng)
        k4f_A = 0.5+2rand(rng); k4r_A = 0.5+2rand(rng)
        k2f_I = 0.5+2rand(rng); k4f_I = 0.5+2rand(rng); k4r_I = 0.5+2rand(rng)
        L = 0.5+rand(rng); Keq = 2.0+2rand(rng)
        A = 0.5+2rand(rng); B = 0.5+2rand(rng); P = 0.5+2rand(rng); Q = 0.5+2rand(rng)
        # Map parameters(m) -> ground-truth params. A step that releases a
        # product while changing the residual is stored as the binding it
        # reverses, F + P → EA, however the step is written; its rate constants
        # are named by its two sides, so `k_…_EA_to_E_res_…_P` is the
        # product-releasing rate and `k_…_E_res_…_P_to_EA` the product-rebinding
        # one — the binding steps read the usual way round.
        #   k_E_A_to_EA=k1f                                  (E + A ⇌ EA, shared)
        #   k_E_res_plusA_minusP_B_to_EB_res_plusA_minusP=k3f,
        #   k_EB_res_plusA_minusP_to_E_res_plusA_minusP_B=k3r  (F + B ⇌ FB, shared)
        #   k_A_EA_to_E_res_plusA_minusP_P=k2f_A,
        #   k_A_E_res_plusA_minusP_P_to_EA=k2r_A             (EA ⇌ F + P, active)
        #   k_I_EA_to_E_res_plusA_minusP_P=k2f_I,
        #   k_I_E_res_plusA_minusP_P_to_EA=k2r_I             (EA ⇌ F + P, inactive)
        #   k_A_EB_res_plusA_minusP_to_E_Q=k4f_A,
        #   k_A_E_Q_to_EB_res_plusA_minusP=k4r_A             (FB ⇌ E + Q, active)
        #   k_I_EB_res_plusA_minusP_to_E_Q=k4f_I,
        #   k_I_E_Q_to_EB_res_plusA_minusP=k4r_I             (FB ⇌ E + Q, inactive)
        # The release of B, `k_EB_res_plusA_minusP_to_E_res_plusA_minusP_B`, and the
        # inactive P rebinding rate, `k_I_E_res_plusA_minusP_P_to_EA`, are absent from
        # parameters(m): each conformation's Haldane relation makes one constant that
        # runs against the reaction dependent. The oracle takes k2r_A and k4r_I and
        # derives k1r and k2r_I by the Haldane relations, so both are computed here the
        # same way; every constant below satisfies both relations, and the derivation
        # rederives k3r and k2r_I from the fitted ones.
        k1r = k1f * k2f_A * k3f * k4f_A / (Keq * k2r_A * k3r * k4r_A)
        k2r_I = k1f * k2f_I * k3f * k4f_I / (Keq * k1r * k3r * k4r_I)
        d = Dict(:k_E_A_to_EA => k1f,
                 :k_EA_to_E_A => k1r,
                 :k_E_res_plusA_minusP_B_to_EB_res_plusA_minusP => k3f,
                 :k_EB_res_plusA_minusP_to_E_res_plusA_minusP_B => k3r,
                 :k_A_EA_to_E_res_plusA_minusP_P => k2f_A,
                 :k_A_E_res_plusA_minusP_P_to_EA => k2r_A,
                 :k_A_EB_res_plusA_minusP_to_E_Q => k4f_A,
                 :k_A_E_Q_to_EB_res_plusA_minusP => k4r_A,
                 :k_I_EA_to_E_res_plusA_minusP_P => k2f_I,
                 :k_I_E_res_plusA_minusP_P_to_EA => k2r_I,
                 :k_I_EB_res_plusA_minusP_to_E_Q => k4f_I,
                 :k_I_E_Q_to_EB_res_plusA_minusP => k4r_I,
                 :L => L)
        prm = NamedTuple{fp}(Tuple(d[s] for s in fp))
        v_code = real(allo_rate((A=A, B=B, P=P, Q=Q), prm, Keq))
        v_gt = _testhelper_pingpong_nonequalAI_freeflip_flux(k1f, k3f, k3r;
            k2f_A=k2f_A, k2r_A=k2r_A, k4f_A=k4f_A, k4r_A=k4r_A,
            k2f_I=k2f_I, k4f_I=k4f_I, k4r_I=k4r_I,
            L=L, Keq=Keq, A=A, B=B, P=P, Q=Q)
        @test isapprox(v_code, v_gt; rtol=1e-10)
    end
end

# ── Ping-pong network with rapid-equilibrium chemistry, at zero products ──────
# The undecorated ping-pong: E + A ⇌ EA → FP ⇌ F + P, F + B ⇌ FB ⇌ EQ ⇌ E + Q,
# with both chemistry steps :OnlyA and every binding :EqualAI. Only EA → FP is
# steady state; the second chemistry FB ⇌ EQ is rapid equilibrium, so free E and
# the covalent F lie in one active rapid-equilibrium segment. Rapid steps use FAST.
# Only free E flips (formulation 1, E_I/E_A = L). The inactive conformation runs no
# chemistry, so it holds E_I, EA_I and EQ_I; its covalent forms, which free E cannot
# reach, hold no mass. At P = Q = 0 the P and Q releases are one-way and the
# chemistry's reverse never fires, so every form past EA → FP drains back to E,
# and the rate is k·(A/KA)/((1 + A/KA)(1 + L)), free of B.
function _testhelper_pingpong_re_chemistry_flux(; KA, KB, KP, KQ, K2, k, L, A, B, FAST=1e9)
    species = [:E_A, :EA_A, :FP_A, :F_A, :FB_A, :EQ_A, :E_I, :EA_I, :EQ_I]
    edges = [
        (:E_A, :EA_A, FAST * A), (:EA_A, :E_A, FAST * KA),   # E + A ⇌ EA (RE)
        (:EA_A, :FP_A, k),                                   # EA → FP (SS)
        (:FP_A, :F_A, FAST * KP),                            # FP → F + P (RE, P = 0)
        (:F_A, :FB_A, FAST * B), (:FB_A, :F_A, FAST * KB),   # F + B ⇌ FB (RE)
        (:FB_A, :EQ_A, FAST * K2), (:EQ_A, :FB_A, FAST),     # FB ⇌ EQ (RE)
        (:EQ_A, :E_A, FAST * KQ),                            # EQ → E + Q (RE, Q = 0)
        (:E_I, :EA_I, FAST * A), (:EA_I, :E_I, FAST * KA),   # inactive A binding
        (:EQ_I, :E_I, FAST * KQ),                            # inactive Q release
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),          # free-enzyme flip, ratio L
    ]
    _testhelper_mwc_ground_truth_flux(species, edges, [(:EA_A, :FP_A, k, 0.0)], 1.0)
end

# ── Ping-pong rapid-equilibrium-chemistry harness self-validation ─────────────
# The network against its closed form, live (L = 3) and unpopulated (L = 0), across
# four decades of B: the B-binding rate FAST·B stays far above k, so the
# rapid-equilibrium limit holds to about k/(FAST·B).
@testset "ping-pong RE-chemistry ground-truth harness self-validation" begin
    p = (KA=0.7, KB=1.3, KP=0.9, KQ=1.1, K2=2.0, k=1.7, A=1.5)
    closed(L) = p.k * (p.A / p.KA) / ((1 + p.A / p.KA) * (1 + L))
    for L in (0.0, 3.0), B in (1e-2, 1.0, 1e2)
        @test isapprox(
            _testhelper_pingpong_re_chemistry_flux(; p..., L=L, B=B), closed(L); rtol=1e-5)
    end
end

# ── The gate: a rapid-equilibrium segment holding free E and F ────────────────
# The two ping-pong value gates above write every step steady state. This gate checks,
# against the exact two-conformation steady state, a mechanism whose active
# rapid-equilibrium segment holds both free E and F. Clearing that segment's polynomial
# gives free E the weight B (F sits at Q/B relative to E), so the inactive term must
# carry the same factor beside L.
@testset "ping-pong MWC derivation with free E and F in one RE segment" begin
    # Both half-reactions :OnlyA: in the inactive conformation the covalent forms between
    # the halves join free E only through the tagged chemistry. The guard admits the
    # mechanism under a reading in which EA <--> EP_res_plusA_minusP loses its stored
    # direction and EB_res_plusA_minusP ⇌ EQ its reverse, both into the covalent forms,
    # and both oracles find the law the limit on every face. Its active rapid-equilibrium
    # segment holds free E and F = E_res_plusA_minusP, and the law equals the exact
    # two-conformation steady state. With the second half alone :OnlyA, the inactive
    # conformation's other steps join the two sides of EB_res_plusA_minusP ⇌ EQ, a
    # rapid-equilibrium step no reading switches off. (Half 1 alone fails class 3 at zero
    # substrates and at B = 0.)
    both = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                          :: EqualAI
            EA <--> EP_res_plusA_minusP                      :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP    :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                         :: OnlyA
            E + Q ⇌ EQ                                          :: EqualAI
        end
    end
    @test all(_testhelper_admission_agreement(both;
                                              rng = MersenneTwister(663), corners = true))
    _testhelper_check_law_against_limit(both)
    _testhelper_check_law_against_oracle(both)
    @test parameters(both) ==
          (:K_A_EB_res_plusA_minusP_to_EQ, :K_EA_to_E_A,
           :K_EB_res_plusA_minusP_to_E_res_plusA_minusP_B,
           :K_EP_res_plusA_minusP_to_E_res_plusA_minusP_P, :K_EQ_to_E_Q,
           :k_A_EA_to_EP_res_plusA_minusP, :L)
    half2 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                          :: EqualAI
            EA <--> EP_res_plusA_minusP                      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP    :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                         :: OnlyA
            E + Q ⇌ EQ                                          :: EqualAI
        end
    end
    @test ER.Derivation._guard(half2) ==
          ER.Derivation._Verdict(ER.Derivation._FAST_STEP_PINNED,
                                 ER.Derivation._LIMIT_UNTESTED, 0, 2,
                                 _testhelper_edge(half2, "EB_res_plusA_minusP ⇌ EQ"))
end

# ── Fully-inert inactive network (both substrate and product bind :OnlyA) ─────
# S and P both bind :OnlyA, so the inactive conformation binds nothing — it is a
# free-enzyme reservoir of mass L, coupled to the active cycle only by the E flip.
# Its partition is 1, so the denominator must carry L·1 (not 0). Fast RE bindings
# and the flip use FAST; catalysis is O(1). Reverse catalysis from the Haldane
# relation. At L = 0 the reservoir is unpopulated → the non-allosteric rate.
function _testhelper_inert_inactive_flux(; KS, KP, k, L, Keq, S, P, FAST=1e7)
    kr = k * KP / (Keq * KS)
    species = [:E_A, :ES_A, :EP_A, :E_I]
    edges = [
        (:E_A, :ES_A, FAST * S / KS), (:ES_A, :E_A, FAST),   # active S binding (RE)
        (:E_A, :EP_A, FAST * P / KP), (:EP_A, :E_A, FAST),   # active P binding (RE)
        (:ES_A, :EP_A, k), (:EP_A, :ES_A, kr),              # active catalysis (SS)
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),          # free-enzyme flip, ratio L
    ]                                                         # E_I inert: no other edges
    _testhelper_mwc_ground_truth_flux(species, edges, [(:ES_A, :EP_A, k, kr)], 1.0)
end

# ── The gate: a fully-inert inactive contributes L to the denominator ─────────
# `den = D_A^n + L`: the inactive conformation is free E alone (`_flip_only_law`:
# Q = 0, D = 1).
@testset "fully-inert inactive contributes L to the denominator" begin
    inert = @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES   :: OnlyA
            ES <--> EP :: EqualAI
            E + P ⇌ EP   :: OnlyA
        end
    end
    # (:K_A_EP_to_E_P, :K_A_ES_to_E_S, :k_ES_to_EP, :L)
    fp = parameters(inert)
    inert_rate = ER.compile_rate_equation(inert)
    @test fp == (:K_A_EP_to_E_P, :K_A_ES_to_E_S, :k_ES_to_EP, :L)

    rng = MersenneTwister(20260714)
    for _ in 1:5
        KS = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng); k = 0.5 + 2rand(rng)
        L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        S = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map parameters(m) -> ground-truth params:
        #   K_A_ES_to_E_S=KS, K_A_EP_to_E_P=KP, k_ES_to_EP=k.
        prm = NamedTuple{fp}((KP, KS, k, L))
        v_code = real(inert_rate((S=S, P=P), prm, Keq))
        @test isapprox(v_code,
            _testhelper_inert_inactive_flux(;
                KS=KS, KP=KP, k=k, L=L, Keq=Keq, S=S, P=P); rtol=1e-4)
        # self-validation: L = 0 → reservoir unpopulated → non-allosteric active rate.
        nonallo = (k * S / KS - (k * KP / (Keq * KS)) * P / KP) / (1 + S / KS + P / KP)
        @test isapprox(
            _testhelper_inert_inactive_flux(; KS=KS, KP=KP, k=k, L=0.0, Keq=Keq, S=S, P=P),
            nonallo; rtol=1e-4)
    end
end

# ── The gate: a ping-pong :OnlyA I-state must keep a reachable free-enzyme root ──
# A covalent intermediate carries no bound metabolite but does carry a residual.
# Seeding I-state reachability from it makes the pruned inactive graph retain a
# covalent island that free E cannot reach, so no spanning tree rooted at free E
# exists and free E's weight is 0 — the normalization then divides by zero and
# the rate law is NaN at every concentration. Under formulation 1 only free
# enzyme flips, so a component free E cannot reach holds no inactive mass and
# must be stranded. With the F6P binding and the first half-reaction :OnlyA, the covalent
# forms E(; residual) and E(F16BP; residual) join free E in the inactive conformation only
# through tagged chemistry: under each consistent reading the F6P binding runs one way
# into E(; residual) at zero substrates, which fills faster than it empties (`_FILLED`),
# the ε-path oracle finding no limit there, and the macro rejects it. With the F6P
# binding and the second half :OnlyA, the rapid-equilibrium E(F6P; residual) ⇌ EADP
# bridges two groups of the inactive conformation (`_BRIDGE`), the ε-path oracle finding
# no limit under its one consistent reading. With that step steady-state, the inactive
# conformation keeps its free-enzyme root and every face passes; the derivation must
# handle it.
@testset "ping-pong :OnlyA I-state keeps a reachable free-enzyme root" begin
    err1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        oligomeric_state: 1
        steps: begin
            E + ATP ⇌ EATP                                                      :: EqualAI
            EATP <--> EF16BP_res_plusATP_minusF16BP                          :: OnlyA
            E_res_plusATP_minusF16BP + F16BP ⇌ EF16BP_res_plusATP_minusF16BP:: EqualAI
            E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP    :: OnlyA
            EF6P_res_plusATP_minusF16BP ⇌ EADP                               :: EqualAI
            E + ADP ⇌ EADP                                                      :: EqualAI
        end
    end
    v = ER.Derivation._guard(err1)
    @test (v.fault, v.face, v.conformation) == (ER.Derivation._FILLED, 3, 2)
    readings = _testhelper_consistent_readings(err1)
    @test length(readings) == 2
    zero_substrates = _testhelper_faces(err1)[3]
    for reading in readings
        @test _testhelper_epath_verdict(err1; rng = MersenneTwister(761),
                                        face = zero_substrates, corners = true, reading) ===
              :direction_dependent
    end
    @test_throws ErrorException(
        ER.Derivation._guard_message(err1)) @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        oligomeric_state: 1
        steps: begin
            E + ATP ⇌ EATP                                                      :: EqualAI
            EATP <--> EF16BP_res_plusATP_minusF16BP                          :: OnlyA
            E_res_plusATP_minusF16BP + F16BP ⇌ EF16BP_res_plusATP_minusF16BP:: EqualAI
            E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP    :: OnlyA
            EF6P_res_plusATP_minusF16BP ⇌ EADP                               :: EqualAI
            E + ADP ⇌ EADP                                                      :: EqualAI
        end
    end
    half2 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        oligomeric_state: 1
        steps: begin
            E + ATP ⇌ EATP                                                      :: EqualAI
            EATP <--> EF16BP_res_plusATP_minusF16BP                          :: EqualAI
            E_res_plusATP_minusF16BP + F16BP ⇌ EF16BP_res_plusATP_minusF16BP:: EqualAI
            E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP    :: OnlyA
            EF6P_res_plusATP_minusF16BP ⇌ EADP                               :: OnlyA
            E + ADP ⇌ EADP                                                      :: EqualAI
        end
    end
    @test ER.Derivation._guard(half2).fault === ER.Derivation._BRIDGE
    readings = _testhelper_consistent_readings(half2)
    @test length(readings) == 1
    @test any(face -> _testhelper_epath_verdict(half2; rng = MersenneTwister(762), face,
                                                reading = only(readings)) !== :converges,
              _testhelper_faces(half2))
    half2_ss = @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        oligomeric_state: 1
        steps: begin
            E + ATP ⇌ EATP                                                      :: EqualAI
            EATP <--> EF16BP_res_plusATP_minusF16BP                          :: EqualAI
            E_res_plusATP_minusF16BP + F16BP ⇌ EF16BP_res_plusATP_minusF16BP:: EqualAI
            E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP    :: OnlyA
            EF6P_res_plusATP_minusF16BP <--> EADP                            :: OnlyA
            E + ADP ⇌ EADP                                                      :: EqualAI
        end
    end
    _testhelper_check_law_against_oracle(half2_ss)

    inactive = ER.Derivation._memo_law(half2_ss).inactive
    @test !iszero(inactive.trees[1] * inactive.shift)

    fp = parameters(half2_ss)
    half2_ss_rate = ER.compile_rate_equation(half2_ss)
    prm = NamedTuple{fp}(Tuple(1.3 for _ in fp))
    Keq = 3.0
    concs = (ATP = 1.1, F6P = 0.7, ADP = 0.6, F16BP = 0.9)
    @test isfinite(real(half2_ss_rate(concs, prm, Keq)))
    # With both products at zero the inactive conformation holds its enzyme in the
    # covalent forms: half 2 cannot run there, and half 1 runs back only with F16BP. The
    # rate is then 0, as the mass-action limit is.
    @test real(half2_ss_rate(merge(concs, (ADP = 0.0, F16BP = 0.0)), prm, Keq)) == 0
    # Far from equilibrium, with a trace of each product, the law runs forward.
    forward = (ATP = 10.0, F6P = 10.0, ADP = 1e-3, F16BP = 1e-3)
    @test real(half2_ss_rate(forward, prm, Keq)) > 0
end

# ── n-protomer concerted-MWC oracle (formulation 1) ─────────────────────────
# Concerted: all protomers share one conformation. Within a conformation the
# protomers are independent, so the joint occupancy state is a tuple. Only the
# FULLY-unliganded oligomer flips (`freeflip=true`) — the n-protomer extension
# of the free-flip-only model this package derives. `freeflip=false` flips every
# joint state (the classic per-form-flip MWC, formulation 2) and is kept ONLY as a
# self-validation reference: it must NOT match the derivation. At nprot = 1 the two
# differ by ~0.1–3% whenever k_A ≠ k_I, because per-form flipping routes turnover
# through the faster conformation.
const OCC = (:E, :EA, :EAB, :EP)

function _testhelper_biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP; k_A, k_I, L, Keq,
                                 A, B, P, FAST=1e7, freeflip=true)
    krA = k_A * kon * KP / (koff * KB * Keq)
    krI = k_I * kon * KP / (koff * KB * Keq)
    _testhelper_prot_edges(kX, krX) = [
        (:E, :EA, kon * A), (:EA, :E, koff),
        (:EA, :EAB, FAST * B / KB), (:EAB, :EA, FAST),
        (:E, :EP, FAST * P / KP), (:EP, :E, FAST),
        (:EAB, :EP, kX), (:EP, :EAB, krX),
    ]
    tbl = Dict(:A => _testhelper_prot_edges(k_A, krA),
               :I => _testhelper_prot_edges(k_I, krI))
    catrate = Dict(:A => (k_A, krA), :I => (k_I, krI))

    occs = collect(Iterators.product(ntuple(_ -> OCC, nprot)...))
    sp(conf, o) = Symbol(conf, "_", join(o, "_"))
    species = Symbol[sp(conf, o) for conf in (:A, :I) for o in occs]
    setidx(o, i, v) = ntuple(j -> j == i ? v : o[j], length(o))

    edges = Tuple{Symbol,Symbol,Float64}[]
    cat_edges = Tuple{Symbol,Symbol,Float64,Float64}[]
    for conf in (:A, :I), o in occs, i in 1:nprot
        for (f, t, r) in tbl[conf]
            o[i] == f || continue
            push!(edges, (sp(conf, o), sp(conf, setidx(o, i, t)), r))
        end
        if o[i] == :EAB
            kf, kr = catrate[conf]
            push!(cat_edges, (sp(conf, o), sp(conf, setidx(o, i, :EP)), kf, kr))
        end
    end
    empty_o = ntuple(_ -> :E, nprot)
    if freeflip
        push!(edges, (sp(:A, empty_o), sp(:I, empty_o), FAST * L))
        push!(edges, (sp(:I, empty_o), sp(:A, empty_o), FAST))
    else
        for o in occs
            push!(edges, (sp(:A, o), sp(:I, o), FAST * L))
            push!(edges, (sp(:I, o), sp(:A, o), FAST))
        end
    end
    _testhelper_mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# Self-validation. Check (c) is the load-bearing one: it pins the oracle to
# formulation 1. Checks (a), (b) and (d) pass for BOTH formulations and so cannot
# distinguish them on their own; they run for both `freeflip` settings.
@testset "concerted-MWC oligomer oracle self-validation" begin
    rng = MersenneTwister(11)
    for nprot in (1, 2, 3), freeflip in (true, false), _ in 1:3
        kon = 0.5+2rand(rng); koff = 0.5+2rand(rng)
        KB = 0.5+2rand(rng); KP = 0.5+2rand(rng)
        kA = 0.5+2rand(rng); kI = 0.5+2rand(rng); Keq = 2.0+2rand(rng)
        A = 0.5+2rand(rng); B = 0.5+2rand(rng); P = 0.5+2rand(rng)
        L = 0.5+rand(rng)
        base = _testhelper_metab_dfree_base_flux(kon, koff, KB, KP, kA, Keq, A, B, P)

        # (a) L = 0 : inactive unpopulated -> nprot x the single-protomer rate.
        @test isapprox(_testhelper_biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=0.0, Keq=Keq, A=A, B=B, P=P, freeflip),
            nprot * base; rtol=1e-4)

        # (b) k_I = k_A : conformations identical -> L-independent.
        f1 = _testhelper_biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kA, L=L, Keq=Keq, A=A, B=B, P=P, freeflip)
        f5 = _testhelper_biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kA, L=5.0, Keq=Keq, A=A, B=B, P=P, freeflip)
        @test isapprox(f1, nprot * base; rtol=1e-4)
        @test isapprox(f1, f5; rtol=1e-4)

        # (d) v = 0 at the equilibrium metabolite ratio.
        Peq = Keq * A * B
        @test abs(_testhelper_biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=L, Keq=Keq, A=A, B=B, P=Peq, freeflip)) < 1e-6
    end

    # (e) a live inactive conformation must move the flux, or the gate below
    #     would pass without ever exercising the cross term.
    v_live = _testhelper_biuni_mwc_oligomer_flux(2, 1.7, 1.1, 0.8, 0.9;
            k_A=2.5, k_I=0.4, L=0.7, Keq=3.0, A=1.1, B=0.5, P=0.6)
    v_dead = _testhelper_biuni_mwc_oligomer_flux(2, 1.7, 1.1, 0.8, 0.9;
            k_A=2.5, k_I=0.0, L=0.7, Keq=3.0, A=1.1, B=0.5, P=0.6)
    @test !isapprox(v_live, v_dead; rtol=1e-3)

    # (c) THE formulation-1 pin: at nprot = 1 the oracle must reproduce the
    #     established free-flip reference, and must NOT equal the per-form-flip
    #     model. Without this, a formulation-2 oracle would pass (a), (b) and (d)
    #     and then disagree with the derivation by 0.1-3% for live :NonequalAI —
    #     a real number that is not a bug.
    args = (1.7, 1.1, 0.8, 0.9)
    kw = (k_A=2.5, k_I=0.4, L=0.7, Keq=3.0, A=1.1, B=0.5, P=0.6)
    @test isapprox(_testhelper_biuni_mwc_oligomer_flux(1, args...; kw...),
                   _testhelper_biuni_nonequalAI_freeflip_flux(args...; kw...); rtol=1e-4)
    @test !isapprox(_testhelper_biuni_mwc_oligomer_flux(1, args...; freeflip=false, kw...),
                    _testhelper_biuni_nonequalAI_freeflip_flux(args...; kw...); rtol=1e-4)
end

# ── The gate: :NonequalAI catalysis for 1, 2 and 3 protomers ──────────────────
# The productive-inactive case. Catalysis is :NonequalAI (k_A ≠ k_I) with shared
# graph topology, so D_A and D_I differ only by the catalytic rate constant. Under
# formulation 1 (the enzyme commits to one conformation per catalytic cycle) the
# two conformations normalize per-state, so the derivation cross-weights them. The
# gate goes red if the derivation drops the cross weights (which
# matches a per-form-flip model, not formulation 1) or mis-renders the
# normalization. At nprot = 1 the oligomer oracle is the free-flip-only reference
# `_testhelper_biuni_nonequalAI_freeflip_flux` (pinned in the oligomer self-validation).
#
# The ^n combine also needs a LIVE inactive numerator: `:OnlyA` always yields a dead
# inactive cycle (the guard forces an `:OnlyA` catalytic tag alongside an `:OnlyA`
# binding), so the numerator's inactive term `L*F_A^n*Q_I*D_I^(n-1)*W_I` is live only for
# `:NonequalAI`. This is the only gate that exercises it, and the only mass-action
# gate at n >= 2 for any family.
@testset ":NonequalAI MWC derivation matches mass-action ground truth, 1-3 protomers" begin
    rng = MersenneTwister(20260716)
    for nprot in (1, 2, 3)
        allo = nprot == 1 ?
            @enzyme_mechanism(begin
                substrates: A, B ; products: P ; oligomeric_state: 1
                steps: begin
                    E + A <--> EA        :: EqualAI
                    EA + B ⇌ EAB     :: EqualAI
                    EAB <--> EP      :: NonequalAI
                    E + P ⇌ EP           :: EqualAI
                end
            end) : nprot == 2 ?
            @enzyme_mechanism(begin
                substrates: A, B ; products: P ; oligomeric_state: 2
                steps: begin
                    E + A <--> EA        :: EqualAI
                    EA + B ⇌ EAB     :: EqualAI
                    EAB <--> EP      :: NonequalAI
                    E + P ⇌ EP           :: EqualAI
                end
            end) :
            @enzyme_mechanism(begin
                substrates: A, B ; products: P ; oligomeric_state: 3
                steps: begin
                    E + A <--> EA        :: EqualAI
                    EA + B ⇌ EAB     :: EqualAI
                    EAB <--> EP      :: NonequalAI
                    E + P ⇌ EP           :: EqualAI
                end
            end)
        fp = parameters(allo)
        allo_rate = ER.compile_rate_equation(allo)
        @test fp == (:K_EAB_to_EA_B, :K_EP_to_E_P, :k_A_EAB_to_EP, :k_EA_to_E_A,
                     :k_E_A_to_EA, :k_I_EAB_to_EP, :L)
        # Map parameters(m) -> ground-truth params:
        #   k_E_A_to_EA=kon, k_EA_to_E_A=koff, K_EAB_to_EA_B=KB, K_EP_to_E_P=KP,
        #   k_A_EAB_to_EP=k_A, k_I_EAB_to_EP=k_I.
        function _testhelper_params_for(kon, koff, KP, KB, kA, kI, L)
            d = Dict(:k_E_A_to_EA=>kon, :k_EA_to_E_A=>koff, :K_EP_to_E_P=>KP,
                     :K_EAB_to_EA_B=>KB, :k_A_EAB_to_EP=>kA, :k_I_EAB_to_EP=>kI, :L=>L)
            NamedTuple{fp}(Tuple(d[s] for s in fp))
        end
        for _ in 1:6
            kon = 0.5+2rand(rng); koff = 0.5+2rand(rng)
            KP = 0.5+2rand(rng); KB = 0.5+2rand(rng)
            kA = 0.5+2rand(rng); kI = 0.5+2rand(rng)
            L = 0.5+rand(rng); Keq = 2.0+2rand(rng)
            A = 0.5+2rand(rng); B = 0.5+2rand(rng); P = 0.5+2rand(rng)
            prm = _testhelper_params_for(kon, koff, KP, KB, kA, kI, L)
            # The rate law is per active site; the oracle is per oligomer.
            v_code = nprot * real(allo_rate((A=A, B=B, P=P), prm, Keq))
            v_gt = _testhelper_biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=L, Keq=Keq, A=A, B=B, P=P)
            @test isapprox(v_code, v_gt; rtol=1e-4)
        end

        if nprot == 1
            # B = 0 zeroes the metabolite B in free E's weight. The active
            # steady-state A binding leaves a bare koff term in free E's weight, so a
            # reverse path from product back to free enzyme stays open and the reverse
            # flux survives.
            kon, koff, KP, KB = 1.7, 1.1, 0.9, 0.8
            kA, kI, L, Keq, A, P = 2.5, 0.4, 0.7, 3.0, 1.1, 0.9
            vN = real(allo_rate((A=A, B=0.0, P=P),
                _testhelper_params_for(kon, koff, KP, KB, kA, kI, L), Keq))
            @test isapprox(vN, _testhelper_biuni_nonequalAI_freeflip_flux(kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=L, Keq=Keq, A=A, B=0.0, P=P); rtol=1e-4)
            @test abs(vN) > 1e-3
        end
    end
end

# ── The printed law ───────────────────────────────────────────────────────────────────

# The unit rule's normalizer of the plain mechanism `m`'s law, the law and the concrete
# mechanism.
function _testhelper_normalizer(m)
    L = EnzymeRates.Derivation._memo_law(m)
    EnzymeRates.Derivation._normalizer(
        L.D, L.vars, EnzymeRates.Derivation._substrate_vars(m, L.vars)), L, m
end

# The concentration monomial of law `L` with the exponents `x => e`.
_testhelper_conc_mono(L, ps::Pair{Symbol}...) =
    EnzymeRates.Derivation._monomial((L.vars.index[x] => e for (x, e) in ps)...)

# The right-hand side of the `N = …` line of mechanism `m`'s printed law.
_testhelper_n_text(m) =
    chopprefix(only(filter(startswith("N = "), split(rate_equation_string(m), '\n'))),
               "N = ")

# Both printed forms of `m`, pasted into a fresh module where `params`, `concs` and `Keq`
# are bound at random values drawn from `rng`, run top to bottom and give the compiled
# law's value. `at_zero`, when a metabolite name, sets that concentration to 0.
function _testhelper_check_pasted_forms(m, rng; at_zero = nothing)
    params, Keq = _testhelper_random_independent_params(m; rng)
    concs = NamedTuple{Tuple(metabolites(m))}(
        Tuple(x === at_zero ? 0.0 : 0.1 + 9.9 * rand(rng) for x in metabolites(m)))
    rate = compile_rate_equation(m)
    for coefficients in (false, true)
        text = rate_equation_string(m; coefficients)
        lines = split(text, '\n')
        @test endswith(lines[1], ") = params")
        @test endswith(lines[2], ") = concs")
        @test startswith(last(lines), "v = ")
        if coefficients
            @test !any(in(("# Haldane constraints:", "# Wegscheider constraints:")),
                       lines)
        end
        mod = Module()
        Core.eval(mod, :(params = $params))
        Core.eval(mod, :(concs = $concs))
        Core.eval(mod, :(Keq = $Keq))
        @test isapprox(Base.include_string(mod, text), rate(concs, params, Keq);
                       rtol = 1e-10)
    end
end

# ── The unit rule ────────────────────────────────────────────────────────────

@testset "Michaelis–Menten takes its constant term" begin
    # Rule 1: the denominator has a term free of concentrations, the free enzyme's, and
    # N is its coefficient, so D sums the forms' abundances relative to free E.
    mm = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            ES <--> E + P
        end
    end
    nz, L, mech = _testhelper_normalizer(mm)
    @test nz == [EnzymeRates.Derivation._Mono() => 1 // 1]
    _testhelper_check_units(mech, L)
end

@testset "Segel Ping Pong Bi Bi takes Cleland's normalizer" begin
    # Rule 2: no constant term, and the denominator holds ΠS/A = B, ΠS/B = A and ΠS = A·B,
    # so N is Cleland's, c(B)·c(A)/c(A·B): the weights B:1, A:1, A·B:−1. The graded basis
    # A, B, P, Q, A·B puts the same weights on the same terms, so N has its value.
    pp = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA <--> F + P
            F + B <--> FB
            FB <--> E + Q
        end
    end
    nz, L, mech = _testhelper_normalizer(pp)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz == [m(:B => 1) => 1 // 1, m(:A => 1) => 1 // 1,
                 m(:A => 1, :B => 1) => -1 // 1]
    @test _testhelper_n_text(pp) ==
          "(k_EA_to_E_A + k_EA_to_F_P) * (k_EA_to_E_A * k_E_Q_to_FB * k_F_P_to_EA * Keq " *
          "+ k_EA_to_F_P * k_E_A_to_EA * k_F_B_to_FB) / (k_EA_to_F_P + k_FB_to_E_Q)"
    _testhelper_check_units(mech, L)
end

@testset "Segel Bi Uni Uni Uni PP Ter Bi takes Cleland's normalizer" begin
    # The denominator holds ΠS/X for each substrate, B·C, A·C and A·B, and ΠS = A·B·C, so
    # N is Cleland's, c(B·C)·c(A·C)·c(A·B)/c(A·B·C)², before the graded basis, which takes
    # the degree-1 term C first (C:2, A·B:1, A·C:−1, B·C:−1) and prints the product of two
    # sums. Cleland's N is that product over a third sum squared.
    terbi = @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EAB
            EAB <--> F + P
            F + C <--> FC
            FC <--> E + Q
        end
    end
    nz, L, mech = _testhelper_normalizer(terbi)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz == [m(:B => 1, :C => 1) => 1 // 1, m(:A => 1, :C => 1) => 1 // 1,
                 m(:A => 1, :B => 1) => 1 // 1,
                 m(:A => 1, :B => 1, :C => 1) => -2 // 1]
    _testhelper_check_units(mech, L)
end

@testset "the ping-pong power case takes Cleland's normalizer" begin
    # Census ping-pong mechanism 6e46f55b0b68e5d9. Its denominator holds B, A and A·B, so
    # N is Cleland's, with the weights B:1, A:1, A·B:−1, before the graded basis, which
    # meets A² before A·B and puts its weights on A and A².
    power = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EP_res_plusA_minusP, EQ + A <--> EAQ)
            (E + Q <--> EQ, EP_res_plusA_minusP + Q <--> EAQ)
            (E_res_plusA_minusP + B <--> EQ,
             EP_res_plusA_minusP + B <--> EBP_res_plusA_minusP)
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EQ + P ⇌ EBP_res_plusA_minusP)
        end
    end
    nz, L, mech = _testhelper_normalizer(power)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz == [m(:B => 1) => 1 // 1, m(:A => 1) => 1 // 1,
                 m(:A => 1, :B => 1) => -1 // 1]
    _testhelper_check_units(mech, L)
end

@testset "_term_key puts A² before A·B" begin
    # A² and A·B are both degree 2 with no non-substrate; the exponent vector puts the
    # earlier name's larger exponent first, so the graded basis meets A² first.
    vars = EnzymeRates.Derivation._Vars()
    a = EnzymeRates.Derivation._var!(vars, :A; conc = true)
    b = EnzymeRates.Derivation._var!(vars, :B; conc = true)
    ab, a2 = EnzymeRates.Derivation._monomial(a => 1, b => 1),
        EnzymeRates.Derivation._monomial(a => 2)
    @test sort([ab, a2]; by = mo -> EnzymeRates.Derivation._term_key(mo, vars, [a, b])) ==
          [a2, ab]
end

@testset "F1 takes no normalizer" begin
    # Rule 4: two half-reactions between the free forms E and F, each A → P. Every term
    # holds A or P, and no weights on A and P solve the units equation, so no product of
    # the coefficients has the denominator's unit, and the law prints as derived.
    f1 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
        end
    end
    nz, L, mech = _testhelper_normalizer(f1)
    @test isempty(nz)
    _testhelper_check_units(mech, L)
end

@testset "F2 takes fractional weights" begin
    # F1 with F binding two A at rapid equilibrium: the terms A, P, A³ and A²·P. The basis
    # A, P, A³ solves the units equation with A:3/2, A³:−1/2.
    f2 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A + A ⇌ FAA
        end
    end
    nz, L, mech = _testhelper_normalizer(f2)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz == [m(:A => 1) => 3 // 2, m(:A => 3) => -1 // 2]
    _testhelper_check_units(mech, L)
end

@testset "F3, a rapid-equilibrium Theorell–Chance step, takes no normalizer" begin
    # Every term holds B or P (the rapid-equilibrium step EA + B ⇌ EQ + P ties them),
    # and the units equation has no solution on the basis B, P, A·B, B·Q, so no product of
    # the coefficients has the denominator's unit, each concentration counting as a
    # dimension of its own, and the law prints as derived.
    f3 = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B ⇌ EQ + P
            EQ <--> E + Q
        end
    end
    nz, L, mech = _testhelper_normalizer(f3)
    @test isempty(nz)
    _testhelper_check_units(mech, L)
end

@testset "a covalent ping-pong without an A·B term takes no normalizer" begin
    # Rule 4: the denominator holds A, B, P and Q, so ΠS/B = A and ΠS/A = B are there but
    # ΠS = A·B is not, and Cleland's rule does not apply. No weights on four degree-1
    # terms solve the units equation, so no product of the coefficients has the
    # denominator's unit, and the law prints as derived.
    covalent = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E_res_plusA_minusP + P
            E_res_plusA_minusP + B <--> E + Q
        end
    end
    nz, L, mech = _testhelper_normalizer(covalent)
    @test isempty(nz)
    _testhelper_check_units(mech, L)
    @test _testhelper_printed_law_error(covalent) <= big(10)^-20
end

@testset "F4 takes integer weights A:2, A²:−1" begin
    # F1 with F binding A twice in two rapid-equilibrium steps: the terms A, P, A², A·P,
    # A³ and A²·P. The basis A, P, A² solves the units equation with A:2, A²:−1.
    f4 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A ⇌ FA
            FA + A ⇌ FAA
        end
    end
    nz, L, mech = _testhelper_normalizer(f4)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz == [m(:A => 1) => 2 // 1, m(:A => 2) => -1 // 1]
    _testhelper_check_units(mech, L)
end

@testset "each MWC conformation takes its own rule" begin
    # An MWC law normalizes each conformation's polynomials by that conformation's rule,
    # dividing its Q, D and free-E weight alike, so the MWC combination does not change. A
    # two-subunit ping-pong whose chemistry differs between the conformations takes
    # Cleland's normalizer in both, each on its own terms. F3 with every step :OnlyA, beside
    # a steady-state dead-end inhibitor I that binds free E in both conformations, takes no
    # normalizer in the active conformation and the constant term in the inactive one, which
    # keeps free E and EI. With A's binding and the Theorell–Chance step alone :OnlyA, the
    # inactive conformation would keep free E and EQ, but the rapid-equilibrium
    # Theorell–Chance step, which then lies between two of its groups, fills EA faster than
    # EA empties at zero substrates, and the guard refuses it; under none of its consistent
    # readings does the ε-path oracle find a limit.
    filled = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA            :: OnlyA
            EA + B ⇌ EQ + P        :: OnlyA
            EQ <--> E + Q            :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(filled)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._FILLED, 3, 2)
    readings = _testhelper_consistent_readings(filled)
    @test length(readings) == 2
    rng = Random.MersenneTwister(31)
    for reading in readings
        @test any(face -> _testhelper_epath_verdict(filled; rng, face, reading) !==
                          :converges, _testhelper_faces(filled))
    end
    pp2 = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A <--> EA       :: EqualAI
            EA <--> F + P       :: NonequalAI
            F + B <--> FB       :: EqualAI
            FB <--> E + Q       :: EqualAI
        end
    end
    tc = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A <--> EA            :: OnlyA
            EA + B ⇌ EQ + P        :: OnlyA
            EQ <--> E + Q            :: OnlyA
            E + I <--> EI            :: EqualAI
        end
    end
    for (i, am) in enumerate((pp2, tc))
        s = EnzymeRates.Derivation._solve(am)
        A = EnzymeRates.Derivation._law(am, s)
        subs = EnzymeRates.Derivation._substrate_vars(am, A.vars)
        nA = EnzymeRates.Derivation._normalizer(A.active.D, A.active.vars, subs)
        nI = EnzymeRates.Derivation._normalizer(A.inactive.D, A.inactive.vars, subs)
        m(ps...) = _testhelper_conc_mono(A.active, ps...)
        if i == 1
            cleland = [m(:B => 1) => 1 // 1, m(:A => 1) => 1 // 1,
                       m(:A => 1, :B => 1) => -1 // 1]
            @test nA == nI == cleland
        else
            @test isempty(nA)
            @test nI == [EnzymeRates.Derivation._Mono() => 1 // 1]
        end
        _testhelper_check_units(am, A)
        # The normalized polynomials in the MWC combination give the law. That holds for
        # any nonzero N_A and N_I, so it pins the invariance; the rule and weight pins
        # above test the normalizer itself.
        rng = Random.MersenneTwister(32 + i)
        for _ in 1:3
            _, concs, params = _testhelper_law_inputs(am, s.fitted, rng)
            env = merge(params, concs)
            value(p) = _testhelper_poly_value(p, A.vars, env)
            NA = _testhelper_norm_value(nA, A.active.D, A.active.vars, env)
            NI = _testhelper_norm_value(nI, A.inactive.D, A.inactive.vars, env)
            n, l = A.n, env[:L]
            qa, da, fa = value(A.active.Q) / NA, value(A.active.D) / NA, value(A.fA) / NA
            qi, di, fi = (value(A.inactive.Q) / NI, value(A.inactive.D) / NI,
                          value(A.fI) / NI)
            num = fi^n * qa * da^(n - 1) + l * fa^n * qi * di^(n - 1)
            normalized = value(A.H) * num / (fi^n * da^n + l * fa^n * di^n)
            law = _testhelper_evaluate_law(A, env)
            @test abs(normalized - law) <= big(10)^-20 * abs(law)
        end
    end
end

@testset "every spec's law has one unit per denominator and a rate per time" begin
    # Units on the 42 specs, each at one catalytic subunit.
    for spec in MECHANISM_TEST_SPECS
        @testset "$(spec.name)" begin
            am = _testhelper_single_subunit(spec.mechanism)
            _testhelper_check_units(
                am, EnzymeRates.Derivation._law(am, EnzymeRates.Derivation._solve(am)))
        end
    end
end

@testset "the printed law and its coefficient form equal the law on every spec" begin
    # Both forms, read back in 256-bit BigFloat, equal the derived law to 1e-20 at a
    # consistent point of every spec, and neither form holds a printing artifact.
    for spec in MECHANISM_TEST_SPECS
        @testset "$(spec.name)" begin
            m = spec.mechanism
            full = rate_equation_string(m)
            @test isempty(_testhelper_print_artifacts(full))
            @test isempty(_testhelper_print_artifacts(
                rate_equation_string(m; coefficients = true)))
            @test _testhelper_printed_law_error(m; points = 1) <= big(10)^-20
        end
    end
end

@testset "both printed forms paste and evaluate as Julia" begin
    # Each form, pasted into a fresh module where `params` and `concs` are bound, runs top
    # to bottom: it opens with their destructuring lines and ends with the `v` line, whose
    # value is the compiled law's. The coefficient form leaves the reference lines to the
    # full law. The cases are every spec and a two-subunit ping-pong with a regulatory
    # site, whose coefficient form defines N_A and N_I before the coefficients they divide
    # and W_A and W_I before `v`; the coefficient-name testset below checks F1-F4 and its
    # other hand-written laws the same way.
    regulated_pp = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A <--> EA       :: EqualAI
            EA <--> F + P       :: NonequalAI
            F + B <--> FB       :: EqualAI
            FB <--> E + Q       :: EqualAI
        end
        regulatory_site(2): R::NonequalAI
    end
    coef = split(rate_equation_string(regulated_pp; coefficients = true), '\n')
    @test [l[1:3] for l in coef if occursin(r"^[QDFWN]_[AI] = ", l)] ==
          ["N_A", "N_I", "Q_A", "D_A", "Q_I", "D_I", "F_A", "F_I", "W_A", "W_I"]
    cases = [[(spec.name, spec.mechanism) for spec in MECHANISM_TEST_SPECS];
             [("regulated ping-pong", regulated_pp)]]
    rng = Random.MersenneTwister(42)
    for (name, m) in cases
        @testset "$name" begin
            _testhelper_check_pasted_forms(m, rng)
        end
    end
end

@testset "coefficient names read powers, and N prints the normalizer once" begin
    # F1's denominator holds only A and P, each to the first power, so it takes no
    # normalizer: its coefficients print as derived, and no N line prints.
    f1 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
        end
    end
    c1 = split(rate_equation_string(f1; coefficients = true), '\n')
    @test c1[1:2] == ["(; k_E_A_to_F_P, k_F_A_to_E_P, k_F_P_to_E_A) = params",
                      "(; A, P) = concs"]
    @test last(c1) == "v = (n_A * A + n_P * P) * (A - P / Keq) / " *
                      "(d_A * A + d_P * P)"
    @test c1[3] == "# The coefficients are functions of the fitted constants, not " *
                   "parameters:"
    @test c1[4:end - 1] == [
        "n_A = 2 * k_E_A_to_F_P * k_F_A_to_E_P",
        "n_P = 2 * k_E_A_to_F_P * k_F_A_to_E_P / Keq",
        "d_A = (k_E_A_to_F_P + k_F_A_to_E_P)",
        "d_P = (k_E_A_to_F_P * k_F_A_to_E_P + k_F_P_to_E_A^2 * Keq^2) / " *
        "(k_F_P_to_E_A * Keq^2)"]
    # F3's rapid-equilibrium Theorell–Chance step makes every denominator term carry B or
    # P, so it takes no normalizer: its coefficients print as derived, and no N line
    # prints.
    f3 = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B ⇌ EQ + P
            EQ <--> E + Q
        end
    end
    c3 = split(rate_equation_string(f3; coefficients = true), '\n')
    @test last(c3) == "v = n_0 * (A * B - P * Q / Keq) / (d_B * B + d_P * P + " *
                      "d_A_B * A * B + d_A_P * A * P + d_B_Q * B * Q + d_P_Q * P * Q)"
    @test "n_0 = k_EA_to_E_A * k_E_Q_to_EQ * Keq" in c3
    @test "d_B = k_EA_to_E_A * k_E_Q_to_EQ * Keq / k_E_A_to_EA" in c3
    @test !any(startswith("N = "), c3)
    # F2's weights are fractional: N is a power 3/2, and its coefficients hold roots.
    f2 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A + A ⇌ FAA
        end
    end
    @test occursin("\nN = (k_E_A_to_F_P + k_F_A_to_E_P)^(3/2)\n", rate_equation_string(f2))
    @test "n_A = 2 * sqrt(K_F_A_A_to_FAA) * k_E_A_to_F_P^(3/2) * k_F_A_to_E_P / N" in
          split(rate_equation_string(f2; coefficients = true), '\n')
    # F2 with its two steps in one group: the weight 3/2 falls on the coefficient of A,
    # whose integer content is 2. An integer to an integer power divides the coefficients;
    # to the power 3/2 it prints apart, as N.
    f2g = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            (E + A <--> F + P, F + A <--> E + P)
            F + A + A ⇌ FAA
        end
    end
    @test occursin("\nN = (2)^(3/2)\n", rate_equation_string(f2g))
    # The random-order bi-bi's denominator holds A^2 * B, named by its superscript.
    spec_mechanism(name) = only(s for s in MECHANISM_TEST_SPECS if s.name == name).mechanism
    ro = spec_mechanism("Random-order Bi-Bi")
    @test occursin(" + d_A²_B * A^2 * B + ", rate_equation_string(ro; coefficients = true))
    # The unit-rule testsets above pin the weights of F4, Segel's ter bi and the ping-pong
    # power case; N prints the normalizer's polynomial factor alone, its integer
    # and monomial content moved into the coefficients. F4's is a sum squared. Segel's,
    # Cleland's normalizer, is the product of two sums over a third sum squared. The power
    # case's, Cleland's too, is one sum. Segel's two ter-ter ping-pongs and the hexa uni
    # ping-pong take Cleland's normalizer too: the product of two sums, or of three for the
    # hexa uni, over a sum of three products squared.
    f4 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A ⇌ FA
            FA + A ⇌ FAA
        end
    end
    @test occursin("\nN = (k_E_A_to_F_P + k_F_A_to_E_P)^2\n", rate_equation_string(f4))
    terbi = spec_mechanism("Segel Bi Uni Uni Uni PP Ter Bi")
    bi_uni_uni_bi = spec_mechanism("Segel Bi Uni Uni Bi PP Ter Ter")
    bi_bi_uni_uni = spec_mechanism("Segel Bi Bi Uni Uni PP Ter Ter")
    hexa_uni = spec_mechanism("Segel Hexa Uni Ping Pong")
    power = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EP_res_plusA_minusP, EQ + A <--> EAQ)
            (E + Q <--> EQ, EP_res_plusA_minusP + Q <--> EAQ)
            (E_res_plusA_minusP + B <--> EQ,
             EP_res_plusA_minusP + B <--> EBP_res_plusA_minusP)
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EQ + P ⇌ EBP_res_plusA_minusP)
        end
    end
    @test _testhelper_n_text(terbi) ==
          "(k_EAB_to_EA_B + k_EAB_to_F_P) * (k_EAB_to_EA_B * k_EA_to_E_A * k_E_Q_to_FC * " *
          "k_F_P_to_EAB * Keq + k_EAB_to_F_P * k_EA_B_to_EAB * k_E_A_to_EA * " *
          "k_F_C_to_FC) / (k_EAB_to_F_P + k_FC_to_E_Q)^2"
    @test _testhelper_n_text(power) ==
          "K_EBP_res_plusA_minusP_to_EQ_P * " *
          "k_EP_res_plusA_minusP_B_to_EBP_res_plusA_minusP * " *
          "k_EP_res_plusA_minusP_to_E_A * k_E_A_to_EP_res_plusA_minusP + " *
          "K_EBP_res_plusA_minusP_to_EQ_P * " *
          "k_EP_res_plusA_minusP_B_to_EBP_res_plusA_minusP * k_EQ_to_E_Q * " *
          "k_E_A_to_EP_res_plusA_minusP + k_EP_res_plusA_minusP_to_E_A * k_EQ_to_E_Q * " *
          "k_E_Q_to_EQ * Keq + k_EP_res_plusA_minusP_to_E_A^2 * k_E_Q_to_EQ * Keq"
    @test _testhelper_n_text(bi_uni_uni_bi) ==
          "(k_EAB_to_EA_B + k_EAB_to_F_P) * (k_EAB_to_EA_B * k_EA_to_E_A * k_ER_Q_to_FC " *
          "* k_E_R_to_ER * k_F_P_to_EAB * Keq + k_EAB_to_F_P * k_EA_B_to_EAB * " *
          "k_ER_to_E_R * k_E_A_to_EA * k_F_C_to_FC) / (k_EAB_to_F_P * k_ER_to_E_R + " *
          "k_EAB_to_F_P * k_FC_to_ER_Q + k_ER_to_E_R * k_FC_to_ER_Q)^2"
    @test _testhelper_n_text(bi_bi_uni_uni) ==
          "(k_EAB_to_EA_B + k_EAB_to_FQ_P) * (k_FC_to_E_R + k_FC_to_F_C) / " *
          "(k_EAB_to_FQ_P * k_FC_to_E_R + k_EAB_to_FQ_P * k_FQ_to_F_Q + " *
          "k_FC_to_E_R * k_FQ_to_F_Q)^2"
    @test _testhelper_n_text(hexa_uni) ==
          "(k_EA_to_E_A + k_EA_to_F_P) * (k_FB_to_F_B + k_FB_to_G_Q) * (k_EA_to_E_A * " *
          "k_E_R_to_GC * k_FB_to_F_B * k_F_P_to_EA * k_G_Q_to_FB * Keq + k_EA_to_F_P * " *
          "k_E_A_to_EA * k_FB_to_G_Q * k_F_B_to_FB * k_G_C_to_GC) / (k_EA_to_F_P * " *
          "k_FB_to_G_Q + k_EA_to_F_P * k_GC_to_E_R + k_FB_to_G_Q * k_GC_to_E_R)^2"
    rng = Random.MersenneTwister(51)
    for m in (f1, f2, f2g, f3, f4, terbi, power, bi_uni_uni_bi, bi_bi_uni_uni, hexa_uni)
        @test _testhelper_printed_law_error(m) <= big(10)^-20
        @test isempty(_testhelper_print_artifacts(rate_equation_string(m))) &&
              isempty(_testhelper_print_artifacts(
                  rate_equation_string(m; coefficients = true)))
        _testhelper_check_pasted_forms(m, rng)
    end
end

@testset "a law with no normalizer prints finite where a concentration is zero" begin
    # The printed law divides by no concentration, so both forms equal the compiled law at
    # a zero of each one. F1, F3 and the covalent ping-pong take no normalizer, and
    # neither do the active conformation of `tc`, F3 with every step :OnlyA beside a
    # steady-state dead-end inhibitor I, and the inactive one of `pp`; the rule pin keeps
    # `pp` an MWC law whose other conformation has a normalizer. `pp`'s inactive and `tc`'s
    # active conformation print their polynomials as derived, with no `N` line, and `tc`'s
    # inactive one, free E and EI, prints its constant term's D_I line.
    f1 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
        end
    end
    f3 = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B ⇌ EQ + P
            EQ <--> E + Q
        end
    end
    covalent = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E_res_plusA_minusP + P
            E_res_plusA_minusP + B <--> E + Q
        end
    end
    tc = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A <--> EA            :: OnlyA
            EA + B ⇌ EQ + P        :: OnlyA
            EQ <--> E + Q            :: OnlyA
            E + I <--> EI            :: EqualAI
        end
    end
    pp = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A <--> E_res_plusA_minusP + P :: EqualAI
            E + Q <--> EQ :: EqualAI
            E_res_plusA_minusP + B <--> EQ :: OnlyA
        end
        regulatory_site(2): R::NonequalAI
    end
    law = EnzymeRates.Derivation._law(pp, EnzymeRates.Derivation._solve(pp))
    subs = EnzymeRates.Derivation._substrate_vars(pp, law.vars)
    mono(ps...) = _testhelper_conc_mono(law.active, ps...)
    # Cleland's terms normalize the active law; no product of coefficients has the
    # inactive denominator's unit.
    @test EnzymeRates.Derivation._normalizer(law.active.D, law.active.vars, subs) ==
          [mono(:B => 1) => 1 // 1, mono(:A => 1) => 1 // 1,
           mono(:A => 1, :B => 1) => -1 // 1]
    @test isempty(EnzymeRates.Derivation._normalizer(law.inactive.D, law.inactive.vars,
                                                     subs))
    full = split(rate_equation_string(pp), '\n')
    coef = split(rate_equation_string(pp; coefficients = true), '\n')
    heads(t) = [l[1:3] for l in t if occursin(r"^[QDFWN]_[AI] = ", l)]
    @test heads(full) == ["N_A", "Q_A", "D_A", "D_I", "F_A", "F_I", "W_A", "W_I"]
    @test heads(coef) == heads(full)
    @test "D_I = A * k_EQ_to_E_Q * k_E_A_to_E_res_plusA_minusP_P + " *
          "P * k_A_E_res_plusA_minusP_B_to_EQ * k_EQ_to_E_Q^2 * " *
          "k_E_A_to_E_res_plusA_minusP_P / " *
          "(k_A_EQ_to_E_res_plusA_minusP_B * k_E_Q_to_EQ * Keq) + " *
          "P * Q * k_A_E_res_plusA_minusP_B_to_EQ * k_EQ_to_E_Q * " *
          "k_E_A_to_E_res_plusA_minusP_P / (k_A_EQ_to_E_res_plusA_minusP_B * Keq)" in full
    @test "F_I = P * k_A_E_res_plusA_minusP_B_to_EQ * k_EQ_to_E_Q^2 * " *
          "k_E_A_to_E_res_plusA_minusP_P / " *
          "(k_A_EQ_to_E_res_plusA_minusP_B * k_E_Q_to_EQ * Keq)" in full
    @test ["dI_A = k_EQ_to_E_Q * k_E_A_to_E_res_plusA_minusP_P",
           "dI_P = k_A_E_res_plusA_minusP_B_to_EQ * k_EQ_to_E_Q^2 * " *
           "k_E_A_to_E_res_plusA_minusP_P / " *
           "(k_A_EQ_to_E_res_plusA_minusP_B * k_E_Q_to_EQ * Keq)",
           "dI_P_Q = k_A_E_res_plusA_minusP_B_to_EQ * k_EQ_to_E_Q * " *
           "k_E_A_to_E_res_plusA_minusP_P / (k_A_EQ_to_E_res_plusA_minusP_B * Keq)",
           "fI_P = k_A_E_res_plusA_minusP_B_to_EQ * k_EQ_to_E_Q^2 * " *
           "k_E_A_to_E_res_plusA_minusP_P / " *
           "(k_A_EQ_to_E_res_plusA_minusP_B * k_E_Q_to_EQ * Keq)"] ⊆ coef
    tc_full = split(rate_equation_string(tc), '\n')
    tc_coef = split(rate_equation_string(tc; coefficients = true), '\n')
    @test heads(tc_full) == ["Q_A", "D_A", "D_I", "F_A"]
    @test heads(tc_coef) == heads(tc_full)
    @test "D_I = 1 + I * k_E_I_to_EI / k_EI_to_E_I" in tc_full
    @test ["dA_P = k_A_EA_to_E_A * k_EI_to_E_I", "dA_P_Q = k_A_E_Q_to_EQ * k_EI_to_E_I",
           "dI_I = k_E_I_to_EI / k_EI_to_E_I"] ⊆ tc_coef
    rng = Random.MersenneTwister(41)
    for m in (f1, f3, covalent, tc, pp), x in metabolites(m), _ in 1:2
        _testhelper_check_pasted_forms(m, rng; at_zero = x)
    end
end

@testset "an :OnlyI inhibitor prints a law equal to its limit" begin
    # I binds only the inactive conformation, at rapid equilibrium and at
    # steady state. Both printed forms equal the rapid-equilibrium limit of the oracle's
    # own network (`_testhelper_limit_law`) to 1e-20 at two points; the active
    # conformation's denominator holds no constant of I's binding, the inactive one's does.
    competitive_RE = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + I ⇌ EI            :: OnlyI
        end
    end
    competitive_SS = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: NonequalAI
            E + I <--> EI         :: OnlyI
        end
    end
    cases = ((competitive_RE, "D_I = 1 + S / K_ES_to_E_S + I / K_I_EI_to_E_I",
              "dI_I = 1 / K_I_EI_to_E_I"),
             (competitive_SS,
              "D_I = 1 + S / K_ES_to_E_S + I * k_I_E_I_to_EI / k_I_EI_to_E_I",
              "dI_I = k_I_E_I_to_EI / k_I_EI_to_E_I"))
    for (i, (m, d_i, d_i_coefficient)) in enumerate(cases)
        full = rate_equation_string(m)
        coef = rate_equation_string(m; coefficients = true)
        @test "D_A = 1 + S / K_ES_to_E_S" in split(full, '\n')
        @test d_i in split(full, '\n')
        @test last(split(full, '\n')) ==
              "v = (Q_A + L * Q_I) * (S - P / Keq) / (D_A + L * D_I)"
        @test "D_A = 1 + dA_S * S" in split(coef, '\n')
        @test d_i_coefficient in split(coef, '\n')
        cm = m.catalytic
        onlyi =
            [String(x) for (k, ns) in enumerate(_testhelper_edge_names(m, :I))
             if m.allosteric.allosteric_states[cm.steps[k].kinetic_group] === :OnlyI
             for x in ns]
        d_a = only(l for l in split(full, '\n') if startswith(l, "D_A = "))
        @test !isempty(onlyi)
        @test !any(x -> occursin(x, d_a), onlyi)
        @test any(x -> occursin(x, d_i), onlyi)
        s = EnzymeRates.Derivation._solve(m)
        rng = Random.MersenneTwister(500 + i)
        for _ in 1:2
            point, concs, params = _testhelper_law_inputs(m, s.fitted, rng)
            exact = setprecision(BigFloat, 256) do
                _testhelper_limit_law(m, point, concs, Set{Symbol}(); rng)
            end
            for text in (full, coef)
                v = _testhelper_eval_rate_string(text, params, concs, params[:Keq])
                @test abs(v - exact) <= big(10)^-20 * abs(exact)
            end
        end
        @test _testhelper_printed_law_error(m) <= big(10)^-20
    end
end

@testset "the MWC lines name each conformation's polynomials" begin
    # A two-site homodimer prints Q, D and W for each conformation and no F: both
    # free-enzyme weights are 1. test_rate_equation_derivation.jl pins its text.
    two_site = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP     :: NonequalAI
            E + S ⇌ ES     :: NonequalAI
            ES <--> EP   :: NonequalAI
        end
        regulatory_site(2): R::NonequalAI
    end
    # A ping-pong whose second half runs in the active conformation only, its last step
    # steady-state: neither conformation's denominator has a constant term, so each prints
    # its own normalizer, and free E weighs F6P and F16BP in one and F16BP alone in the
    # other, so F_A and F_I print. With that last step at rapid equilibrium it bridges two
    # groups of the inactive conformation, and the guard refuses it.
    bridged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        oligomeric_state: 1
        steps: begin
            E + ATP ⇌ EATP                                                      :: EqualAI
            EATP <--> EF16BP_res_plusATP_minusF16BP                          :: EqualAI
            E_res_plusATP_minusF16BP + F16BP ⇌ EF16BP_res_plusATP_minusF16BP:: EqualAI
            E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP    :: OnlyA
            EF6P_res_plusATP_minusF16BP ⇌ EADP                               :: OnlyA
            E + ADP ⇌ EADP                                                      :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(bridged).fault === EnzymeRates.Derivation._BRIDGE
    half2 = @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        oligomeric_state: 1
        steps: begin
            E + ATP ⇌ EATP                                                      :: EqualAI
            EATP <--> EF16BP_res_plusATP_minusF16BP                          :: EqualAI
            E_res_plusATP_minusF16BP + F16BP ⇌ EF16BP_res_plusATP_minusF16BP:: EqualAI
            E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP    :: OnlyA
            EF6P_res_plusATP_minusF16BP <--> EADP                            :: OnlyA
            E + ADP ⇌ EADP                                                      :: EqualAI
        end
    end
    # Both bindings :OnlyA: the inactive conformation is free E alone, so the
    # denominator's inactive term is L with no factor.
    inert = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: OnlyA
        end
    end
    heads(t) = [l[1:3] for l in split(t, '\n') if occursin(r"^[QDFWN]_[AI] = ", l)]
    @test heads(rate_equation_string(two_site)) ==
          ["Q_A", "D_A", "Q_I", "D_I", "W_A", "W_I"]
    t2 = rate_equation_string(half2)
    @test heads(t2) == ["N_A", "Q_A", "D_A", "N_I", "D_I", "F_A", "F_I"]
    @test "F_A = F6P / (K_A_EF6P_res_plusATP_minusF16BP_to_E_res_plusATP_minusF16BP_F6P " *
          "* k_EATP_to_EF16BP_res_plusATP_minusF16BP * N_A) + F16BP * K_EADP_to_E_ADP / " *
          "(K_A_EF6P_res_plusATP_minusF16BP_to_E_res_plusATP_minusF16BP_F6P * " *
          "K_EATP_to_E_ATP * k_A_EADP_to_EF6P_res_plusATP_minusF16BP * Keq * N_A)" in
          split(t2, '\n')
    @test "F_I = F16BP / " *
          "(K_A_EF6P_res_plusATP_minusF16BP_to_E_res_plusATP_minusF16BP_F6P * " *
          "K_EATP_to_E_ATP * " *
          "K_EF16BP_res_plusATP_minusF16BP_to_E_res_plusATP_minusF16BP_F16BP * " *
          "k_A_EADP_to_EF6P_res_plusATP_minusF16BP * Keq * N_I)" in split(t2, '\n')
    @test last(split(t2, '\n')) == "v = (F_I * Q_A) * (ATP * F6P - ADP * " *
                                    "F16BP / Keq) / (F_I * D_A + L * F_A * D_I)"
    ti = rate_equation_string(inert)
    _testhelper_short_heads(m) = [l[1:3] for l in split(_testhelper_short_law(m), '\n')
                                  if occursin(r"^[QDFWN]_[AI] = ", l)]
    @test _testhelper_short_heads(two_site) == ["Q_A", "D_A", "Q_I", "D_I", "W_A", "W_I"]
    @test _testhelper_short_heads(half2) == ["Q_A", "D_A", "D_I", "F_A", "F_I"]
    @test heads(ti) == ["Q_A", "D_A"]
    @test last(split(ti, '\n')) == "v = (Q_A * D_A) * (S - P / Keq) / (D_A^2 + L)"
    for m in (two_site, half2, inert)
        @test _testhelper_printed_law_error(m) <= big(10)^-20
        @test isempty(_testhelper_print_artifacts(rate_equation_string(m))) &&
              isempty(_testhelper_print_artifacts(
                  rate_equation_string(m; coefficients = true)))
    end
end

@testset "the printed-law reader runs a law only as written" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    s = rate_equation_string(m)
    p, Keq = _testhelper_random_independent_params(m; rng = Random.MersenneTwister(1))
    c = (S = 0.4, P = 0.9)
    @test isapprox(_testhelper_eval_rate_string(s, p, c, Keq),
                   compile_rate_equation(m)(c, p, Keq); rtol = 1e-12)
    # A name the params or concs line leaves out, and a line read before the lines that
    # bind its names, fail as they do in Julia.
    for t in (replace(s, "K_ES_to_E_S, " => ""; count = 1),
              replace(s, "(; S, P) = concs" => "(; S) = concs"),
              join(reverse(split(s, '\n')), '\n'))
        @test_throws "before it binds it" _testhelper_eval_rate_string(t, p, c, Keq)
    end
end

@testset "the short law names the coefficients and leaves their values out" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    @test _testhelper_short_law(m) == "(; S, P) = concs\n" *
        "v = n_0 * (S - P / Keq) / (1 + d_S * S + d_P * P)"
    @test rate_equation_string(m; coefficients = true,
                               coefficient_values = false) == _testhelper_short_law(m)
    @test_throws ArgumentError rate_equation_string(m; coefficient_values = false)
    @test_throws ArgumentError rate_equation_string(m;
                                                    coefficient_values = false)

    # The MWC short law keeps the conformations' definitions in their coefficients.
    inert = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: OnlyA
        end
    end
    @test _testhelper_short_law(inert) == "(; S, P) = concs\n" *
        "Q_A = nA_0\n" *
        "D_A = 1 + dA_S * S + dA_P * P\n" *
        "v = (Q_A * D_A) * (S - P / Keq) / (D_A^2 + L)"

    # Two conformations that differ, with a regulatory site: each prints its Q and D, and
    # each site factor is the binding polynomial (1 + R / K)^2 of a two-site homodimer.
    regulated = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP     :: NonequalAI
            E + S ⇌ ES     :: NonequalAI
            ES <--> EP   :: NonequalAI
        end
        regulatory_site(2): R::NonequalAI
    end
    @test _testhelper_short_law(regulated) == "(; S, P, R) = concs\n" *
        "Q_A = nA_0\n" *
        "D_A = 1 + dA_S * S + dA_P * P\n" *
        "Q_I = nI_0\n" *
        "D_I = 1 + dI_S * S + dI_P * P\n" *
        "W_A = (1 + R / K_A_Rreg)^2\n" *
        "W_I = (1 + R / K_I_Rreg)^2\n" *
        "v = (Q_A * D_A * W_A + L * Q_I * D_I * W_I) * (S - P / Keq) / " *
        "(D_A^2 * W_A + L * D_I^2 * W_I)"

    # An aggregate check over every spec and both large-law lists: the short form is the
    # coefficient form minus its `params` line, its comment and the value lines (`n…`,
    # `d…`, `f…`, `N…` on the left), in order.
    # A coefficient's value line: `n_0 = `, `dA_S = `, `fI_B = `, `N = `, `N_A = `.
    value_line = r"^(?:[ndf][AI]?(?:_\S*)?|N(?:_[AI])?) = "
    for m in [[s.mechanism for s in MECHANISM_TEST_SPECS];
              _testhelper_largest_enumerated_mechanisms;
              _testhelper_large_ter_bi_mechanisms]
        full = split(rate_equation_string(m; coefficients = true), '\n')
        expected = [l for l in full[2:end]
                    if !startswith(l, "# ") && !occursin(value_line, l)]
        short = split(_testhelper_short_law(m), '\n')
        @test short == expected
        @test startswith(short[1], "(; ") && endswith(short[1], ") = concs")
        @test startswith(last(short), "v = ")
    end
end

# ── The allosteric golden reference ───────────────────────────────────────────────────

const _ALLO_GOLDEN_PATH =
    joinpath(@__DIR__, "..", "reference", "allosteric_golden_reference.txt")

"""Canonical serialization of every allosteric spec's derivation output."""
function _testhelper_allosteric_golden_lines()
    lines = String[]
    for spec in MECHANISM_TEST_SPECS
        _testhelper_is_allosteric(spec.mechanism) || continue
        m = spec.mechanism
        push!(lines, "### " * spec.name)
        eq_string = EnzymeRates.rate_equation_string(m)
        push!(lines, "RATE_EQUATION_STRING " * replace(eq_string, "\n" => "\\n"))
        push!(lines, "RATE_EQUATION_COEFFICIENTS " *
                     replace(EnzymeRates.rate_equation_string(m; coefficients = true),
                             "\n" => "\\n"))
    end
    lines
end

@testset "allosteric golden reference (D1)" begin
    # The per-spec mass-action oracle of test_rate_equation_derivation.jl
    # (`_testhelper_test_mass_action_oracle`) checks each spec's law, with one catalytic
    # subunit, against the exact steady state of the spec's mass-action network; the
    # per-spec string check (`_testhelper_test_rate_equation_string`) ties the printed
    # string to the compiled law, and test_rate_equation_derivation.jl ties both printed
    # forms to the derived law.
    @test isfile(_ALLO_GOLDEN_PATH)
    current = _testhelper_allosteric_golden_lines()
    reference = readlines(_ALLO_GOLDEN_PATH)
    @test length(current) == length(reference)
    for (c, r) in zip(current, reference)
        @test c == r
    end
    # The fitted parameters list the catalytic constants by name, so rate-equivalent
    # mechanisms render one params line, then the regulator constants, then L.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        _testhelper_is_allosteric(m) || continue
        fp = parameters(m)
        reg = Tuple(K for st in (:A, :I) for K in _testhelper_regulator_constants(m, st))
        cat = fp[1:end - 1 - count(in(fp), reg)]
        @test issorted(collect(cat); by = string)
        @test fp == (cat..., filter(in(fp), reg)..., :L)
    end
end
