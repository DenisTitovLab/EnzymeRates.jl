# ABOUTME: Tests for rate-equation derivation (King-Altman/Cha, mass-action oracle),
# ABOUTME: the parameters API, rescaling to turnovers, rate-equation strings, perf gate.

using OrdinaryDiffEqFIRK

include("../reference/allosteric_undefvar_reproducers.jl")

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

    D_total = sum(D)
    E_conc = D ./ D_total .* params.E_total

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
matches by ligand-name set so positional Kreg names (`:K_<lig>_reg{site}`) keep
the oracle's written site numbering after canonical site reordering.
"""
function _testhelper_positional_site_idx(site, pos, written_sites)
    want = Set(site.ligands)
    j = findfirst(ws -> Set(ws.ligands) == want, written_sites)
    j === nothing ? pos : j
end

"""
The constants of the plain mechanism `m` (values by name in `all_params`, `params` merged
with `dependent_parameter_values(m, params)`) under the positional names the QSSA and ODE
cross-checks read: stored flat step `i` gives `K{i}` at rapid equilibrium and `k{i}f`,
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
`dependent_parameter_values(m, params)`) under the positional names a textbook rate law
reads, by role. Written step `i` (in the order of `raw`'s written steps) takes the
constants of the edge it became (`_canonicalize`), in the written direction: `k{i}f` and `k{i}r` at steady state, and at rapid equilibrium `K{i}`, a
binding's dissociation constant or, for any other step, the equilibrium constant of the
written direction. The inactive law's kept steps in a `:NonequalAI` group add `k{i}f_T`,
`k{i}r_T` or `K{i}_T`. No spec in `MECHANISM_TEST_SPECS` writes a non-binding step against
its stored direction, so the suite never inverts `K` or swaps `k{i}f` and `k{i}r` for such a
step. A ligand of regulatory site `j` (the written order of `raw`'s sites) gives
`K_<ligand>_reg{j}` and, when its inactive conformation has a constant of its own,
`K_<ligand>_T_reg{j}`. `raw` is the mechanism as written, the DSL walker's value
(`_testhelper_MechanismTestSpec.raw`); errors unless its catalytic steps canonicalize to
`m`'s catalytic graph. Every entry of `all_params` passes through.
"""
function _testhelper_written_params(m, all_params::NamedTuple; raw)
    mech = _testhelper_lift(m)
    allo = mech isa EnzymeRates.AllostericMechanism
    g = _testhelper_graph(mech)
    c = _testhelper_canonicalize(raw)
    c.mechanism == EnzymeRates.Mechanisms.catalytic(mech) ||
        error("raw writes another mechanism than m")
    # The test reader's steps run group by group, each group in edge order.
    flat = reduce(
        vcat, EnzymeRates.Mechanisms.catalytic(mech).groups)
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
        for (pos, site) in enumerate(mech.regulatory_sites)
            j = _testhelper_positional_site_idx(site, pos, raw.regulatory_sites)
            for (lig, tag) in zip(site.ligands,
                                  site.mwc_tags)
                for (inactive, suffix) in ((false, ""), (true, "_T"))
                    tag === :EqualAI && inactive && continue
                    K = EnzymeRates.Derivation.regulator_constant(mech, lig, tag, inactive)
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

"""The call the gate times: `rate_equation` behind a function boundary the compiler may not
inline, so it cannot hoist the call out of the timing loop."""
@noinline _testhelper_rate_barrier(m, concs, params) = rate_equation(m, concs, params)

"""
Test that `rate_equation` is non-allocating and fast for the given mechanism.
Must be a standalone function to avoid @testset closure boxing. Each call goes through
`_testhelper_rate_barrier`, so every one of the 10,000 calls per batch runs.
"""
function _testhelper_test_rate_equation_performance(m, params, concs)
    _testhelper_rate_barrier(m, concs, params) # warmup/compile
    allocs = @allocated _testhelper_rate_barrier(m, concs, params)
    # Minimum over several batches defeats the GC/scheduling inflation a
    # single mean suffers; accumulating into `acc` (and observing it via the
    # finite-result check) prevents the optimizer from eliding the calls.
    best = Inf
    acc = 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += _testhelper_rate_barrier(m, concs, params)
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
    param_keys = (parameters(m, FittedParameters())..., :Keq, :E_total)
    param_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in param_keys)
    new_params = NamedTuple{param_keys}(param_vals)

    conc_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in met_names)
    concs = NamedTuple{Tuple(met_names)}(conc_vals)

    all_params = merge(new_params, dependent_parameter_values(m, new_params))
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
    m::EnzymeMechanism,
    params, concs,
)
    g = _testhelper_graph(m)
    E_total = params.E_total
    n = length(g.forms)
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
        @test metabolites(m) == Tuple(spec.metabolite_names)
        # Structural parameter names must be injective on the load-bearing
        # paths the fitter consumes (parameters + the fitted parameters). A collision
        # would silently shorten the destructured params NamedTuple.
        @test allunique(EnzymeRates.parameters(m))
        @test allunique(parameters(m, FittedParameters()))
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
        dep_exprs, indep = _testhelper_dependent_param_exprs(typeof(m))
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
        # Locally 120 ns leaves about 4× headroom over the slowest fixture law, "LDH
        # i-state NonequalAI 6-group, steady-state NAD group" at about 26 ns, and about 9×
        # over the next, "Random-order Bi-Bi" at about 13 ns. Shared CI runners can run a
        # law 4-7× slower, so a job that misses the bound is rerun rather than the bound
        # loosened. The 0-alloc bound above is the strict part.
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
        @test startswith(last(split(s, "\n")), "v = E_total * ")  # Proper format
        @test occursin(") / (", s)      # Has denominator
        @test occursin("= params", s)   # Has params destructuring
        @test occursin("= concs", s)    # Has concs destructuring
        @test isempty(_testhelper_print_artifacts(s))

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
The three factors of a printed `v` line, `v = E_total * (Q) * (H) / (D)`, as written
inside their parentheses: `Q` (a bare name as it is, `"1"` when the line leaves it out),
`H` and `D`. Each factor is a parenthesized group or a bare name, so the split runs at the
` * (` and ` / (` outside every parenthesis.
"""
function _testhelper_extract_factors(v_line::AbstractString)
    s = String(chopprefix(v_line, "v = E_total * "))
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
        params = _testhelper_random_independent_params(m; rng)
        met_names = metabolites(m)
        concs = NamedTuple{Tuple(met_names)}(Tuple(0.5 + rand(rng) for _ in met_names))
        v_orig = rate_equation(m, concs, params)
        @test isfinite(v_orig) && v_orig != 0

        # Scaling every steady-state rate constant by α scales every rate by α.
        ss_names = parameters(m, EnzymeRates.Derivation.RateConstants())
        α = 0.1 + 9.9 * rand(rng)
        scaled_params = NamedTuple{keys(params)}(Tuple(
            k in ss_names ? v * α : v
            for (k, v) in zip(keys(params), values(params))
        ))
        @test rate_equation(m, concs, scaled_params) ≈ α * v_orig rtol=1e-9

        # An anchor at those concentrations sets the turnover there and leaves every
        # constant but the steady-state rate constants alone.
        target = sign(v_orig) * (0.1 + 9.9 * rand(rng))
        v_unit = rate_equation(m, concs, merge(params, (E_total = 1.0,)))
        anchored = _testhelper_scale_rate_constants(m, params, target / v_unit)
        for k in keys(params)
            if !(k in ss_names)
                @test anchored[k] == params[k]
            end
        end
        # The anchors' turnover is per enzyme, so the scale leaves E_total as passed in.
        @test anchored.E_total == params.E_total
        @test rate_equation(m, concs, merge(anchored, (E_total=1.0,))) ≈ target rtol=1e-9

        # Every rate changes by the same factor.
        other = NamedTuple{Tuple(met_names)}(Tuple(0.5 + rand(rng) for _ in met_names))
        @test rate_equation(m, other, anchored) / rate_equation(m, other, params) ≈
              target / rate_equation(m, concs, merge(params, (E_total=1.0,))) rtol=1e-8
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
        bigfloat = spec.mechanism isa EnzymeRates.AllostericEnzymeMechanism ||
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

"""The allocations of `parameters(m, kind)` once compiled, measured inside a function."""
function _testhelper_parameters_allocations(m, kind)
    parameters(m, kind)
    @allocated parameters(m, kind)
end

@testset "parameters kinds agree across representations" begin
    for spec in MECHANISM_TEST_SPECS
        em = spec.mechanism
        m = _testhelper_lift(em)
        for kind in (IndependentParameters(), FittedParameters(), AllParameters())
            @test parameters(m, kind) == parameters(em, kind)
        end
        @test parameters(em) == parameters(em, IndependentParameters())
        @test parameters(m) == parameters(em)
        @test parameters(em, IndependentParameters()) ==
              (parameters(em, FittedParameters())..., :Keq, :E_total)
        @test parameters(em, AllParameters())[1:length(parameters(em))] == parameters(em)
        @test issorted(parameters(em, AllParameters())[length(parameters(em)) + 1:end])
        @test _testhelper_parameters_allocations(em, FittedParameters()) == 0
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

@testset "each rapid-equilibrium segment is rooted at its first form" begin
    # The engine weighs each form of a rapid-equilibrium segment against the segment's
    # first form in the canonical form order (`rapid_equilibrium_segments`), so free E,
    # form 1, roots its own segment. A's modifier copy comes before Ab there: the table
    # sorts A's copy just after A, before Ab. The law in lowest terms does not depend on the root.
    em = @enzyme_mechanism begin
        substrates: A, Ab
        products: P
        modifiers: A
        steps: begin
            E + A ⇌ E(A)
            E(A) + Ab ⇌ E(A, Ab)
            E(A, Ab) <--> E(P)
            E(P) ⇌ E + P
            E + A::Mod ⇌ E(A::Mod)
            E + Ab ⇌ E(Ab)
            E(Ab) + A ⇌ E(A, Ab)
            E(A) <--> E(; residual = A - P) + P
        end
    end
    m = ER.Mechanism(em)
    form(text) =
        only(f for f in eachindex(m.nodes) if ER.Mechanisms._form_name(m, f) == text)
    a_copy, ab, e = form("EAmod"), form("EAb"), form("E")
    @test e == 1 && a_copy < ab
    s = ER.Mechanisms.rapid_equilibrium_segments(m, trues(length(m.steps)))
    @test all(seg -> seg[1] == minimum(seg), s.segments)
    @test s.segments[s.segment_of[e]][1] == e
    @test s.segment_of[a_copy] == s.segment_of[ab] == s.segment_of[e]
    covalent = form("E_res_plusA_minusP")
    @test s.segments[s.segment_of[covalent]] == [covalent]
    _testhelper_check_law_against_oracle(em)
end

@testset "the free-enzyme weight is free E's weight in the reduced denominator" begin
    # Ping-pong with a rapid-equilibrium second chemistry step: every form lies in one
    # rapid-equilibrium segment rooted at free E, where E(; residual = A - P) sits at
    # Q·K/B relative to E. Reducing the polynomials to lowest terms in the
    # concentrations multiplies them by B, so free E weighs B in the reduced
    # denominator; the MWC law reads this weight as free E's (`_cross_weights`).
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
    m = EnzymeRates.Mechanism(em)
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
    @test Set(parameters(ss, FittedParameters())) ==
          Set([:k_E_S_to_ES, :k_ES_to_E_S, :k_ES_to_E_P])
    re = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
        end
    end
    @test Set(parameters(re, FittedParameters())) == Set([:K_ES_to_E_S, :k_ES_to_E_P])
end

@testset "fused steps derive the mass-action rate" begin
    for em in _testhelper_fused_cases
        _testhelper_check_law_against_oracle(em)
    end
    # The B group holds the fused E(A) + B → E(P, Q) and the dead-end E(Q) + B ⇌ E(B, Q).
    # Fitted: A (2) + B (1) + P (2) + Q (1) − 1 Haldane = 5.
    @test length(parameters(_testhelper_fused_and_plain_binding, FittedParameters())) == 5
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
    @test length(parameters(tc_ss, FittedParameters())) == 5
    @test length(parameters(tc_re_outer, FittedParameters())) == 3
    for em in (tc_ss, tc_re_outer, tc_re_step, two_in)
        _testhelper_check_law_against_oracle(em)
    end
    # canonical orientation: the substrate side is `from`
    m = ER.Mechanism(tc_ss)
    e = only(e for e in m.steps
             if !ER.Mechanisms._is_iso(e) && !ER.Mechanisms._is_binding(e))
    @test ER.Mechanisms._form_name(m, e.from) == "EA" &&
          e.taken_up == [ER.Mechanisms.Met(:B, false)]
    @test :k_EA_B_to_EQ_P in _testhelper_step_constant_names(tc_ss)
    @test :k_EQ_P_to_EA_B in _testhelper_step_constant_names(tc_ss)
end

@testset "reversing written steps changes nothing" begin
    rng = Random.MersenneTwister(7)
    # The written form `raw` (`_testhelper_written`) with every step reversed, a random
    # subset, and every group's steps in reverse order.
    function _testhelper_variants(raw)
        group = [e.group for e in raw.steps]
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
        for flipped in _testhelper_variants(_testhelper_written(m))
            m2 = _testhelper_build(flipped)
            @test m2 == m
            @test ER.rate_equation_string(ER.compile_mechanism(m2)) ==
                  ER.rate_equation_string(em)
        end
    end
    for spec in MECHANISM_TEST_SPECS
        spec.mechanism isa ER.AllostericEnzymeMechanism || continue
        am = ER.AllostericMechanism(spec.mechanism)
        for flipped in _testhelper_variants(_testhelper_written(am))
            am2 = _testhelper_build(flipped)
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
@testset "Allosteric group tag stays bound to its step" begin
    hk = only(s for s in MECHANISM_TEST_SPECS if s.name == "HK")
    am = EnzymeRates.AllostericMechanism(hk.mechanism)
    cm = am.catalytic
    members = cm.groups
    onlyA_groups = [g for g in eachindex(members)
                    if am.mwc_tags[g] === :OnlyA]
    @test length(onlyA_groups) == 2
    onlyA_firsts = [cm.steps[first(members[g])] for g in onlyA_groups]
    @test count(e -> EnzymeRates.Mechanisms._is_binding(e) &&
                     only(e.taken_up).name === :ATP,
                onlyA_firsts) == 1
    @test count(!EnzymeRates.Mechanisms._is_binding, onlyA_firsts) == 1
end

# ── §5a regression: every inactive-state parameter must be defined ──────────
# Allosteric mechanisms from an LDH `identify_rate_equation` run crashed with
# `UndefVarError` on undefined I-state parameters: the I-state polynomials
# referenced names the dep-assignment and destructuring machinery never emitted.
# The undefined constants were an I-state lactate dissociation constant and an
# I-state reverse chemistry rate in dead inactive states, and an I-state NAD
# binding rate in a live one. Those mechanisms are not kept here. The regression
# runs the three `LDH i-state …` specs of MECHANISM_TEST_SPECS (defined in
# mechanism_definitions_for_test_enzyme_derivation.jl), whose Haldane-dependent
# reverse rate must stay out of the fitted set, and checks that each rate
# equation defines every name it references.
const _LDH_ISTATE_MECHS = [spec.mechanism for spec in MECHANISM_TEST_SPECS
                           if startswith(spec.name, "LDH i-state")]

# Parameter Symbols referenced on an assignment/`v` RHS but never defined
# (destructured from `params`/`concs` or assigned as an LHS). Empty ⟺ the
# rendered rate equation is closed: every referenced name has a definition.
function _testhelper_undefined_rhs_symbols(s::AbstractString)
    # A name may hold superscripts (a coefficient's); a called function (`sqrt`) is no name.
    ident = r"[A-Za-z_][A-Za-z0-9_⁰¹²³⁴⁵⁶⁷⁸⁹]*(?![A-Za-z0-9_⁰¹²³⁴⁵⁶⁷⁸⁹(])"
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
        pnames = parameters(em, FittedParameters())
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
        @test _testhelper_dep_graph_is_sound(_testhelper_dependent_param_exprs(T)...)
    end
end

@testset "the third reproducer with a rapid-equilibrium Lactate release traps E(NAD)" begin
    # Written with E(NAD) + Lactate ⇌ E(Lactate, NAD) at rapid equilibrium beside the
    # steady-state NAD binding, the mechanism fills E(NAD) by the fast Lactate release at
    # zero products, and only steady-state steps empty it: a trap (`_TRAP`), which the
    # RE-speed oracle sees as a pattern, as it does the split at zero Lactate. The
    # reproducer takes the fix the message prints, that step steady-state, which passes
    # every face and keeps the steady-state speed k_I_E_NAD_to_ENAD dependent.
    NO = EnzymeRates.Derivation._NO_FAULT
    trap = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: NADH, Pyruvate
        products: Lactate, NAD
        catalytic_multiplicity: 4
        steps: begin
            E + Lactate ⇌ E(Lactate) :: EqualAI
            E + NAD <--> E(NAD) :: NonequalAI
            E + NADH ⇌ E(NADH) :: OnlyA
            (E + Pyruvate <--> E(Pyruvate),
             E(NADH) + Pyruvate <--> E(NADH, Pyruvate)) :: EqualAI
            E(Lactate) + NAD ⇌ E(Lactate, NAD) :: EqualAI
            (E(Lactate) + NADH ⇌ E(Lactate, NADH),
             E(Pyruvate) + NADH ⇌ E(NADH, Pyruvate)) :: OnlyA
            E(NAD) + Lactate ⇌ E(Lactate, NAD) :: EqualAI
            E(NAD) + Pyruvate <--> E(NAD, Pyruvate) :: EqualAI
            E(NADH) + Lactate ⇌ E(Lactate, NADH) :: EqualAI
            E(NADH, Pyruvate) <--> E(Lactate, NAD) :: OnlyA
            E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate) :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(trap)
    @test (v.fault, v.limit, v.face) ==
          (EnzymeRates.Derivation._TRAP, EnzymeRates.Derivation._LIMIT_SPEED, 2)
    @test EnzymeRates.Mechanisms.step_text(trap.catalytic, v.edge) ==
          "E(NAD) + Lactate ⇌ E(Lactate, NAD)"
    @test _testhelper_admission_faults(trap) ==
          [NO, EnzymeRates.Derivation._TRAP, NO, NO, NO, EnzymeRates.Derivation._SPLIT, NO]
    @test _testhelper_respeed_verdicts(trap; rng = Random.MersenneTwister(1),
                                       law = :limit) ==
          [:agree, :pattern, :agree, :agree, :agree, :pattern, :agree]
    msg = EnzymeRates.Derivation._guard_message(trap, v)
    @test_throws ErrorException(msg) @allosteric_mechanism begin
        substrates: NADH, Pyruvate
        products: Lactate, NAD
        catalytic_multiplicity: 4
        steps: begin
            E + Lactate ⇌ E(Lactate) :: EqualAI
            E + NAD <--> E(NAD) :: NonequalAI
            E + NADH ⇌ E(NADH) :: OnlyA
            (E + Pyruvate <--> E(Pyruvate),
             E(NADH) + Pyruvate <--> E(NADH, Pyruvate)) :: EqualAI
            E(Lactate) + NAD ⇌ E(Lactate, NAD) :: EqualAI
            (E(Lactate) + NADH ⇌ E(Lactate, NADH),
             E(Pyruvate) + NADH ⇌ E(NADH, Pyruvate)) :: OnlyA
            E(NAD) + Lactate ⇌ E(Lactate, NAD) :: EqualAI
            E(NAD) + Pyruvate <--> E(NAD, Pyruvate) :: EqualAI
            E(NADH) + Lactate ⇌ E(Lactate, NADH) :: EqualAI
            E(NADH, Pyruvate) <--> E(Lactate, NAD) :: OnlyA
            E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate) :: EqualAI
        end
    end
    fixed = EnzymeRates.AllostericMechanism(ALLOSTERIC_UNDEFVAR_REPRODUCERS[3]())
    @test _testhelper_admission_faults(fixed) == fill(NO, 7)
    @test _testhelper_respeed_verdicts(fixed; rng = Random.MersenneTwister(2)) ==
          fill(:agree, 7)
    @test :k_I_E_NAD_to_ENAD in keys(_testhelper_dependent_param_exprs(fixed)[1])
    _testhelper_check_law_against_oracle(ALLOSTERIC_UNDEFVAR_REPRODUCERS[3]())
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
    for T in ALLOSTERIC_UNDEFVAR_REPRODUCERS
        m = T()
        pn = collect(parameters(m, FittedParameters()))
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

"""The body of the generated `rate_equation` of `m`, as the generator builds it."""
_testhelper_rate_body(m) =
    ER.Derivation._rate_body(ER.Derivation._derivation_input(typeof(m)))

"""Each local the generated `rate_equation` body of `m` binds, in order, with its
right-hand side."""
_testhelper_body_locals(m) =
    Pair{Symbol, Any}[ex.args[1] => ex.args[2] for ex in _testhelper_rate_body(m).args
                      if ex isa Expr && ex.head === :(=) && ex.args[1] isa Symbol]

"""The expression the generated `rate_equation` body of `m` returns."""
_testhelper_generated_rate_return(m) =
    only(ex.args[1] for ex in _testhelper_rate_body(m).args
         if ex isa Expr && ex.head === :return)

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
    # Every local the body binds, and what it returns, for a plain law and an MWC law.
    mwc = only(s for s in MECHANISM_TEST_SPECS
               if s.name == "MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites")
    for m in (spec.mechanism, mwc.mechanism)
        foreach(_testhelper_walk!, last.(_testhelper_body_locals(m)))
        _testhelper_walk!(_testhelper_generated_rate_return(m))
    end
    @test isempty(bad)
end

# Seconds per call allowed per denominator term: each of the six gate laws (the fixtures
# file) must run in under the larger of 120 ns and this constant times its denominator
# terms. It is twice the worst CI reading over two runs, rounded up to one significant
# figure; the worst readings on each runner, in ns per term, ter-bi / bi-bi: ubuntu x64
# 0.223 / 0.299, windows x64 0.352 / 0.297, macOS aarch64 0.143 / 0.198.
const _testhelper_seconds_per_denominator_term = 0.8e-9

"""The number of terms in the denominator of `m`'s derived rate law after the solve is
substituted: an MWC law counts both conformations, an inactive conformation that is free E
alone counting 1."""
function _testhelper_denominator_terms(m)
    L = EnzymeRates.Derivation._law(_testhelper_lift(m))
    L isa EnzymeRates.Derivation._Law && return EnzymeRates.Derivation._nterms(L.D)
    EnzymeRates.Derivation._nterms(L.active.D) +
        EnzymeRates.Derivation._nterms(L.inactive.D)
end

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

@testset "rate_equation reads concs and params by name in any order" begin
    # The body reads each input by name (`getfield`), so reordered NamedTuples, a data row
    # with Rate and group fields, and `params` merged with
    # `dependent_parameter_values(m, params)` give the same value without allocating. A metabolite named x1 reads as a
    # concentration beside the body's positional local x1.
    m = only(s for s in MECHANISM_TEST_SPECS if s.name == "Random-order Bi-Bi").mechanism
    params, concs, all_params = _testhelper_random_independent_params_concs(
        m, collect(metabolites(m)); rng = Random.MersenneTwister(91))
    v = rate_equation(m, concs, params)
    flipped(nt) = NamedTuple{reverse(keys(nt))}(reverse(values(nt)))
    row = merge((Rate = 1.0, group = 2), flipped(concs))
    for (c, p) in ((flipped(concs), params), (concs, flipped(params)),
                   (row, all_params), (row, flipped(all_params)))
        @test rate_equation(m, c, p) == v
        allocs, _ = _testhelper_test_rate_equation_performance(m, p, c)
        @test allocs == 0
    end
    x1 = @enzyme_mechanism begin
        substrates: x1
        products: P
        steps: begin
            E + x1 ⇌ E(x1)
            E(x1) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    _testhelper_check_law_against_oracle(x1)
    params, concs, _ = _testhelper_random_independent_params_concs(
        x1, collect(metabolites(x1)); rng = Random.MersenneTwister(92))
    allocs, _ = _testhelper_test_rate_equation_performance(x1, params, concs)
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
            (E <--> F, E(S) <--> F(S), E(P) <--> F(P))
            F + S <--> E(S)
            F(S) + B <--> E(P)
            F(P) <--> E + P
        end
    end
    mech = _testhelper_lift(m)
    @test any(d -> denominator(d.keq) == 3, EnzymeRates.Derivation._solve(mech).dependents)
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
            (E(S) <--> F(S), E <--> F)
            F(S) <--> E + P
            F + S <--> E(S)
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
    # The program defines a local for each dependent that is a power product of two or
    # more symbols before any other local, in the solve's column order, so a body's
    # numbering never follows a hash (spec section 2, criterion 3); a dependent equal to
    # one symbol takes no local. The random-order bi-bi has a Haldane dependent and a
    # Wegscheider one, each a product of several constants.
    m = _testhelper_lift(only(s for s in MECHANISM_TEST_SPECS
                              if s.name == "Random-order Bi-Bi").mechanism)
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
    # dependents" (test_joint_solve.jl): a kinetic-group dependent and a Wegscheider one
    # that each equal one constant, read directly, and a Haldane product, the program's
    # first local. Every later local holds a sum of two or more terms.
    folded = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
            E + B ⇌ E(B)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
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
            L = EnzymeRates.Derivation._law(_testhelper_lift(m))
            rng = Random.MersenneTwister(778 + hash(spec.name) % 1000)
            params = _testhelper_random_independent_params(m; rng)
            mets = collect(metabolites(m))
            for zeroed in mets
                concs = NamedTuple{Tuple(mets)}(
                    Tuple(n == zeroed ? 0.0 : 0.5 + rand(rng) for n in mets))
                v = rate_equation(m, concs, params)
                env = Dict{Symbol, BigFloat}(k => BigFloat(x)
                                             for (k, x) in pairs(merge(params, concs)))
                limit = params.E_total * _testhelper_evaluate_law(L, env)
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
            L = EnzymeRates.Derivation._law(_testhelper_lift(m))
            @test L.n == 4
            rng = Random.MersenneTwister(779 + hash(name) % 1000)
            setprecision(BigFloat, 256) do
                params = _testhelper_random_independent_params(m; rng)
                pbig = map(x -> x * (1 + rand(rng, BigFloat)), params)
                mets = metabolites(m)
                cbig = NamedTuple{Tuple(mets)}(Tuple(big(0.5) + rand(rng, BigFloat)
                                                     for _ in mets))
                env = Dict{Symbol, BigFloat}(pairs(merge(pbig, cbig)))
                law = pbig.E_total * _testhelper_evaluate_law(L, env)
                @test abs(rate_equation(m, cbig, pbig) - law) <= big(10)^-20 * abs(law)
            end
        end
    end
end

"""Whether each reference line of `m`'s printed law, read back in 256-bit BigFloat
(`_testhelper_eval_rate_lines`), equals the dependent constant `dependent_parameter_values(m, params)`
computes, to 1e-12 relative, at a random point; false when the law prints no reference
line."""
function _testhelper_reference_lines_match(m)
    rng = Random.MersenneTwister(4)
    p = _testhelper_random_independent_params(m; rng)
    c = NamedTuple{metabolites(m)}(Tuple(0.5 + rand(rng) for _ in metabolites(m)))
    env = _testhelper_eval_rate_lines(rate_equation_string(m), p, c)
    ap = dependent_parameter_values(m, p)
    dependents = collect(keys(ap))
    !isempty(dependents) &&
        all(k -> abs(env[k] - ap[k]) <= 1e-12 * abs(ap[k]), dependents)
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
    @test _testhelper_reference_lines_match(folded)
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
    @test _testhelper_reference_lines_match(iso_tie)
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
        r"\nk_EQ_to_EB_res_plusA_minusP = k_EA_to_EP_res_plusA_minusP \* ", s)
    @test _testhelper_reference_lines_match(m)
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
    @test occursin("\nk_F_to_E = k_E_to_F * sqrt(k_FS_to_E_P) * sqrt(k_F_S_to_ES) / " *
                   "(sqrt(k_ES_to_F_S) * sqrt(k_E_P_to_FS) * sqrt(Keq))\n", s)
    @test _testhelper_reference_lines_match(half)
    rng = Random.MersenneTwister(4)
    point = setprecision(() -> _testhelper_consistent_point(half, rng), BigFloat, 256)
    fitted = _testhelper_dependent_param_exprs(typeof(half))[2]
    params = NamedTuple{(fitted..., :Keq, :E_total)}(
        (Float64.(getindex.(Ref(point.vals), fitted))..., Float64(point.Keq), 1.0))
    concs = (S = 0.7, P = 1.3)
    @test isapprox(_testhelper_eval_rate_string(s, params, concs),
                   rate_equation(half, concs, params); rtol = 1e-12)
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
    @test startswith(s, "v = E_total * (")
    q, h, d = _testhelper_extract_factors(s)
    @test (_testhelper_sum_depth(q), _testhelper_sum_depth(h), _testhelper_sum_depth(d)) ==
          (1, 0, 1)
    @test h == "A * B - P * Q / Keq"
    @test _testhelper_printed_law_error(spec.mechanism) <= big(10)^-20
end

@testset "RateConstants lists what Fitting scales" begin
    # Every steady-state rate constant of both conformations takes the factor; scaling
    # exactly these multiplies the rate by it, on concrete and compiled mechanisms alike.
    RC = EnzymeRates.Derivation.RateConstants()
    rng = Random.MersenneTwister(5)
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        names = parameters(m, RC)
        @test names isa Tuple{Vararg{Symbol}}
        @test allunique(names) && !isempty(names)
        @test parameters(_testhelper_lift(m), RC) === names
        @test issubset(names, parameters(m, AllParameters()))
        @test isdisjoint(names, (:Keq, :E_total, :L))
        params = _testhelper_random_independent_params(m; rng)
        met_names = metabolites(m)
        concs = NamedTuple{Tuple(met_names)}(Tuple(0.5 + rand(rng) for _ in met_names))
        v = rate_equation(m, concs, params)
        scaled = _testhelper_scale_rate_constants(m, params, 3.0)
        @test rate_equation(m, concs, scaled) ≈ 3.0 * v rtol = 1e-9
    end
    # The exact names, in the order `RateConstants()` returns them: the set is the
    # steady-state constants of every conformation's law (the inactive one's `I_` names
    # included for MWC Dimer), and rapid-equilibrium constants are absent (RE Uni-Uni).
    expected = Dict(
        "Uni-Uni" => (:k_E_P_to_ES, :k_ES_to_E_P, :k_E_S_to_ES, :k_ES_to_E_S),
        "RE Uni-Uni" => (:k_E_P_to_EA, :k_EA_to_E_P),
        "MWC Dimer [AllostericEnzymeMechanism]" =>
            (:k_A_ES_to_EP, :k_A_EP_to_ES, :k_I_ES_to_EP, :k_I_EP_to_ES),
        "m_OnlyA_prod" => (:k_A_ES_to_EP, :k_A_EP_to_ES),
        "LDH i-state NonequalAI 5-group" =>
            (:k_A_E_NADH_to_ENADH, :k_A_ENADH_to_E_NADH,
             :k_A_ENADHPyruvate_to_ELactateNAD, :k_A_ELactateNAD_to_ENADHPyruvate))
    for (name, names) in expected
        m = only(s for s in MECHANISM_TEST_SPECS if s.name == name).mechanism
        @test parameters(m, RC) === names
        @test Set(parameters(_testhelper_lift(m), RC)) == Set(names)
    end
end

@testset "scaling the rate constants multiplies the steady-state ones only" begin
    for name in ("RE Uni-Uni", "MWC Dimer [AllostericEnzymeMechanism]")
        m = only(s for s in MECHANISM_TEST_SPECS if s.name == name).mechanism
        ks = (parameters(m)...,)
        params = NamedTuple{ks}(Tuple(Float64.(eachindex(ks))))
        scaled = _testhelper_scale_rate_constants(m, params, 10.0)
        ss = parameters(m, EnzymeRates.Derivation.RateConstants())
        @test keys(scaled) == ks
        @test !isempty(ss) && any(k -> !(k in ss), ks)
        @test all(scaled[k] == (k in ss ? 10 * params[k] : params[k]) for k in ks)
    end
end

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
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            E(P, Q) ⇌ E(Q) + P
            E(P, Q) ⇌ E(P) + Q
            E(P) <--> E + P
            E(Q) <--> E + Q
        end
    end), false))
    mechs = [[_testhelper_lift(s.mechanism) for s in MECHANISM_TEST_SPECS]; rejected]
    questions = (m -> D.check_valid_rate_equation(m; throw = false),
                 D.check_valid_rate_equation,
                 m -> parameters(m, FittedParameters()), rate_equation_string,
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

# ── Degenerate constraint handling ────────────────────────────────────────────

@testset "Degenerate constraint handling" begin
    # ── Unit tests: a rendered dependent's types ─────────────────

    @testset "a dependent renders as a Symbol, an Expr or 1.0 with exact exponents" begin
        render = _testhelper_power_expr
        R = Rational{BigInt}
        # Single symbol factor (exp=1): returns bare Symbol
        @test render(R(0), [(:k1f, R(1))]) === :k1f

        # Single inverse factor: returns Expr
        @test render(R(0), [(:k1f, R(-1))]) isa Expr

        # Keq-only: returns Symbol
        @test render(R(1), Tuple{Symbol, R}[]) === :Keq

        # Multiple factors: returns Expr
        @test render(R(0), [(:k1f, R(1)), (:k2f, R(1))]) isa Expr

        # No factors and zero Keq: returns the Float64 literal 1.0
        @test render(R(0), Tuple{Symbol, R}[]) === 1.0

        # Keq times a factor is a valid AST node
        @test render(R(1), [(:k1f, R(1))]) isa Union{Float64, Symbol, Expr}

        # A fractional exponent is an exact Rational{Int} literal, so the expression
        # evaluates exactly in BigFloat.
        @test render(R(-1//2), [(:k1f, R(1//2))]) ==
              Expr(:call, :*, Expr(:call, :^, :Keq, -1//2), Expr(:call, :^, :k1f, 1//2))
        @test render(R(-1//2), [(:k1f, R(1//2))]).args[2].args[3] isa Rational{Int}
        half = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                (E(S) <--> F(S), E <--> F)
                F(S) <--> E + P
                F + S <--> E(S)
            end
        end
        dep, fitted = _testhelper_dependent_param_exprs(typeof(half))
        setprecision(BigFloat, 256) do
            p = NamedTuple{(fitted..., :Keq)}(
                Tuple(big(1) + big(i) / 7 for i in 1:(length(fitted) + 1)))
            env = Dict{Symbol, BigFloat}(k => v for (k, v) in pairs(p))
            @test abs(_testhelper_printed_value(dep[:k_F_to_E], env) -
                      dependent_parameter_values(half, p).k_F_to_E) <=
                  big(10)^-60 * dependent_parameter_values(half, p).k_F_to_E
        end
    end

    @testset "a constant thermodynamics pins to 1 prints as 1" begin
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
        @test occursin("# Wegscheider constraints:\nK_E_to_F = 1\n",
                       rate_equation_string(em))
        @test :K_E_to_F ∉ parameters(em, FittedParameters())
        _testhelper_check_law_against_oracle(em)
        # The generated rate equation folds the pinned constant into its coefficients
        # and stays allocation-free.
        syms = parameters(em)
        params = NamedTuple{syms}(ntuple(i -> 0.5 + i / 10, length(syms)))
        @test dependent_parameter_values(em, params).K_E_to_F == 1
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
    # A cycle of rapid-equilibrium steps that turns the reaction over gives the mechanism
    # no finite rate. The guard's class 2 finds it, and the macro throws when its expanded
    # code runs, naming the cycle's steps.
    msg = "close a cycle that turns the reaction over, so the mechanism has no finite rate"
    # Binding stage mixed (S2 SS, S1 RE), release stage mixed (P1 SS, P2 RE), chemistry
    # RE: a complete all-RE catalytic cycle.
    @test_throws msg @enzyme_mechanism begin
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
    # The activator route E(R) + S ⇌ E(S, R) ⇌ E(P, R) ⇌ E(R) + P is an all-RE catalytic
    # cycle.
    @test_throws msg @enzyme_mechanism begin
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
end

@testset "rate_equation_complexity (V×τ term-count estimate)" begin
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
    @test EnzymeRates.Derivation.rate_equation_complexity(raw(ordered_bibi)) == 25

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
    @test EnzymeRates.Derivation.rate_equation_complexity(raw(random_bibi)) == 336

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
    @test EnzymeRates.Derivation.rate_equation_complexity(raw(ordered_terter)) == 49

    # Allosteric: both conformations count. An all-SS uni-uni has G = 3 segments and
    # τ = 3 spanning trees, 9, and its inactive conformation keeps every step: 9 more.
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
    am = EnzymeRates.AllostericMechanism(allo)
    @test _testhelper_laplacian_tree_count(
        am.catalytic, am.active_steps) == 9
    @test _testhelper_laplacian_tree_count(
        am.catalytic, am.inactive_steps) == 9
    @test EnzymeRates.Derivation.rate_equation_complexity(am) == 18
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

# The generated rate_equation names its type parameter M_type, which holds `_` as no
# metabolite name does, and its first argument m; metabolites named M and m read as
# concentrations.
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
        "Mechanism: the metabolite name k_EP_to_ES holds `_`; constant names join " *
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
# and k_I_… for the inactive one. A conformation labelled A or I would render the same
# letter at the start of a step side, so an :EqualAI step could take a :NonequalAI step's
# name; the macros and the constructors reject those two labels, and the derivation
# checks the names of a graph built without them. A metabolite named minusQ imitates a
# residual, so the derivation also rejects a mechanism the macros build.

"""The block of an allosteric mechanism whose second conformation is labelled `label`:
`E(b) <--> label(Eb)` at :NonequalAI, then `label(Eb) <--> label + Eb` and `label ⇌ E`
at :EqualAI."""
_testhelper_relabelled_block(label) = :(begin
    substrates: b
    products:   Eb
    steps: begin
        E + b ⇌ E(b)             :: EqualAI
        E(b) <--> $label(Eb)     :: NonequalAI
        $label(Eb) <--> $label + Eb :: EqualAI
        $label ⇌ E               :: EqualAI
    end
end)

"""The raw value of `_testhelper_relabelled_block(label)`, written field by field so that
no label check reads it."""
function _testhelper_relabelled_raw(label)
    M = EnzymeRates.Mechanisms
    b, eb, none = M.Met(:b, false), M.Met(:Eb, false), M.Met[]
    node(c, bound) = M.Node(c, bound, Symbol[], Symbol[])
    (; substrates = [:b], products = [:Eb], modifiers = Symbol[],
     nodes = [node(:E, none), node(:E, [b]), node(:E, [b]), node(label, [eb]),
              node(label, [eb]), node(label, none), node(label, none), node(:E, none)],
     steps = [M.Step(1, 2, [b], none, true, 1), M.Step(3, 4, none, none, false, 2),
              M.Step(5, 6, none, [eb], false, 3), M.Step(7, 8, none, none, true, 4)],
     mwc_tags = [:EqualAI, :NonequalAI, :EqualAI, :EqualAI], catalytic_multiplicity = 1,
     regulatory_sites = EnzymeRates.DSL._WrittenSite[])
end

"""The allosteric graph `raw` writes, in canonical form and composed with its tags,
multiplicity and sites, with no representation check: the constructors' label check
would reject it."""
function _testhelper_unchecked(raw::NamedTuple)
    c = _testhelper_canonicalize(raw)
    tags = similar(raw.mwc_tags)
    for (l, t) in enumerate(raw.mwc_tags)
        tags[c.group_of[l]] = t
    end
    EnzymeRates.AllostericMechanism(c.mechanism, tags, raw.catalytic_multiplicity,
        [EnzymeRates.Mechanisms.RegulatorySite(s.ligands, s.multiplicity, s.mwc_tags)
         for s in raw.regulatory_sites])
end

@testset "constants that would share a name are rejected" begin
    # The walker reads the block the field-by-field raw value writes.
    @test EnzymeRates.DSL.raw_mechanism(_testhelper_relabelled_block(:F), true) ==
          _testhelper_relabelled_raw(:F)
    for label in (:A, :I)
        reason = "`$label` cannot be a conformation label: constant names tag a " *
                 "conformation with A_ or I_ (`k_A_ES_to_EP`), so the labels A and I " *
                 "could give two constants one name; choose another label."
        @test_throws("@allosteric_mechanism: " * reason,
                     eval(:(@allosteric_mechanism $(_testhelper_relabelled_block(label)))))
        @test_throws("@enzyme_mechanism: " * reason, eval(:(@enzyme_mechanism begin
            substrates: b
            products:   Eb
            steps: begin
                E + b ⇌ E(b)
                E(b) <--> $label(Eb)
                $label(Eb) <--> $label + Eb
                $label ⇌ E
            end
        end)))
        @test_throws(ErrorException("Mechanism: " * reason),
                     EnzymeRates.DSL._construct(_testhelper_relabelled_raw(label)))
    end

    # The derivation checks the names of a graph built without the label check. The
    # uptake constant of the :EqualAI step A(Eb) <--> A + Eb is k_A_Eb_to_AEb, the name
    # of the active constant of the :NonequalAI step E(b) <--> A(Eb).
    am_a = _testhelper_unchecked(_testhelper_relabelled_raw(:A))
    conformation_a = EnzymeRates.compile_mechanism(am_a)
    clash_a = ErrorException("Constants of this mechanism's rate equation share a name " *
                             "(k_A_Eb_to_AEb); rename a metabolite or a conformation " *
                             "label so that each constant has a name of its own")
    params = (K_Eb_to_E_b = 1.0, K_A_to_E = 2.0, k_A_Eb_to_AEb = 3.0, k_AEb_to_A_Eb = 4.0,
              k_I_Eb_to_AEb = 0.5, L = 10.0, Keq = 4.0, E_total = 1.0)
    @test_throws clash_a parameters(conformation_a)
    @test_throws clash_a parameters(conformation_a, FittedParameters())
    @test_throws clash_a rate_equation(conformation_a, (b = 1.0, Eb = 0.1), params)
    @test_throws clash_a rate_equation_string(conformation_a)
    @test_throws clash_a dependent_parameter_values(conformation_a, params)
    @test_throws clash_a EnzymeRates.Derivation.parameter_count(am_a)

    # Labelled I, the :EqualAI uptake constant k_I_Eb_to_IEb is the name of the
    # :NonequalAI step's inactive constant.
    conformation_i =
        EnzymeRates.compile_mechanism(_testhelper_unchecked(_testhelper_relabelled_raw(:I)))
    clash_i = ErrorException("Constants of this mechanism's rate equation share a name " *
                             "(k_I_Eb_to_IEb); rename a metabolite or a conformation " *
                             "label so that each constant has a name of its own")
    @test_throws clash_i parameters(conformation_i)
    @test_throws clash_i rate_equation_string(conformation_i)

    # Labelled F, each constant has a name of its own.
    renamed = eval(:(@allosteric_mechanism $(_testhelper_relabelled_block(:F))))
    @test parameters(renamed) == (:K_Eb_to_E_b, :K_F_to_E, :k_A_Eb_to_FEb, :k_FEb_to_F_Eb,
                                  :k_F_Eb_to_FEb, :k_I_Eb_to_FEb, :L, :Keq, :E_total)
    @test EnzymeRates.AllostericMechanism(renamed) ==
          _testhelper_unchecked(_testhelper_relabelled_raw(:F))

    # A modifier named minusQ imitates a residual: the side E(; residual = A - P) + minusQ
    # and the form E(; residual = A - P - Q) both render E_res_plusA_minusP_minusQ, so each
    # constant of the minusQ binding takes the name of a constant of the step from
    # E(minusQ; residual = A - P) to that form.
    residual_clash = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: minusQ
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
            E(; residual = A - P) + minusQ <--> E(minusQ; residual = A - P)
            E(minusQ; residual = A - P) <--> E(; residual = A - P - Q)
        end
    end
    @test_throws(ErrorException(
        "Constants of this mechanism's rate equation share a name " *
        "(k_EminusQ_res_plusA_minusP_to_E_res_plusA_minusP_minusQ, " *
        "k_E_res_plusA_minusP_minusQ_to_EminusQ_res_plusA_minusP); rename a metabolite " *
        "or a conformation label so that each constant has a name of its own"),
        parameters(residual_clash, FittedParameters()))
end

# The names no metabolite takes: E_total, Keq, L, N, params, concs, v and every name
# holding `_` (test_types.jl checks that the reaction rejects each). Besides the
# metabolites, the printed law binds or reads only such names; the generated body binds
# only positional locals (the next testset), so no concentration can overwrite a constant
# or a local.
@testset "the reserved names cover every non-metabolite name of a law" begin
    reserved = (:E_total, :Keq, :L, :N, :params, :concs, :v)
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
    # read as getfield(concs, :X) of a metabolite or getfield(params, :k) of a parameter,
    # so no metabolite name can shadow a constant or a local.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        body = _testhelper_rate_body(m)
        locals = _testhelper_body_locals(m)
        @test all(p -> occursin(r"^x\d+([mp]\d+)?$", String(p.first)), locals)
        @test _testhelper_names_read(body) ⊆ union(Set(first.(locals)), [:concs, :params])
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

    # The inactive conformation's binding constants are dissociation constants (`S / K`)
    # also when `:OnlyA` and `:NonequalAI` catalytic groups coexist; read as association
    # constants (`K * S`) they would halve the rate below.
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
    # come from the mechanism's own fitted parameters, so the Haldane property (zero net
    # rate at equilibrium) is verified whatever the priority picks.
    fp_ro = parameters(m_ro, FittedParameters())
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
k_A_EP_to_ES = K_A_EP_to_E_P * k_A_ES_to_EP / (K_A_ES_to_E_S * Keq)
k_I_EP_to_ES = K_I_EP_to_E_P * k_I_ES_to_EP / (K_I_ES_to_E_S * Keq)
Q_A = k_A_ES_to_EP / K_A_ES_to_E_S
D_A = 1 + S / K_A_ES_to_E_S + P / K_A_EP_to_E_P
Q_I = k_I_ES_to_EP / K_I_ES_to_E_S
D_I = 1 + S / K_I_ES_to_E_S + P / K_I_EP_to_E_P
W_A = (1 + R / K_A_Rreg)^2
W_I = (1 + R / K_I_Rreg)^2
v = E_total * (Q_A * D_A * W_A + L * Q_I * D_I * W_I) * (S - P / Keq) / (D_A^2 * W_A + L * D_I^2 * W_I)"""
    @test actual == expected
    expected_coefficients = raw"""(; K_A_EP_to_E_P, K_A_ES_to_E_S, K_I_EP_to_E_P, K_I_ES_to_E_S, k_A_ES_to_EP, k_I_ES_to_EP, K_A_Rreg, K_I_Rreg, L, Keq, E_total) = params
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
v = E_total * (Q_A * D_A * W_A + L * Q_I * D_I * W_I) * (S - P / Keq) / (D_A^2 * W_A + L * D_I^2 * W_I)"""
    @test rate_equation_string(m_allo; coefficients = true) == expected_coefficients
    # Both texts equal the derived law (value check of this pin).
    @test _testhelper_printed_law_error(m_allo) <= big(10)^-20
end

@testset "Parameter-struct allosteric helpers" begin
    @testset "_edge_names and the regulator constants" begin
        flat_names(am, state) = reduce(vcat, EnzymeRates.Derivation._edge_names(am, state))
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
        @test _testhelper_regulator_constants(am, :A) == [:K_A_Rreg]
        @test _testhelper_regulator_constants(am, :I) == [:K_I_Rreg]
        @test EnzymeRates.Derivation._regulator_names(am) == [:K_A_Rreg, :K_I_Rreg]

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
        @test all(isempty, EnzymeRates.Derivation._edge_names(am_skip, :I))
        @test isempty(_testhelper_regulator_constants(am_skip, :I))
        @test EnzymeRates.Derivation._regulator_names(am_skip) == [:K_A_Rreg]
        @test _testhelper_regulator_constants(am_skip, :A) ==
              [:K_A_Rreg]

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
        @test isempty(_testhelper_regulator_constants(am_onlyi, :A))
        @test EnzymeRates.Derivation._regulator_names(am_onlyi) == [:K_I_Rreg]
        @test _testhelper_regulator_constants(am_onlyi, :I) ==
              [:K_I_Rreg]
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
        @test _testhelper_regulator_constants(am_equal, :A) ==
              [:K_Rreg]
        @test isempty(_testhelper_regulator_constants(am_equal, :I))
        @test EnzymeRates.Derivation._regulator_names(am_equal) == [:K_Rreg]
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
    indep = parameters(m, FittedParameters())
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
    @test parameters(tie, FittedParameters()) == parameters(grouped, FittedParameters()) ==
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
    pn = parameters(m, FittedParameters())
    params = merge(NamedTuple{pn}(ntuple(_ -> 1.0, length(pn))),
                   (Keq = 1.0, E_total = 1.0))
    concs = (A = 1.0, B = 1.0, P = 1.0, Q = 1.0)
    @test isfinite(rate_equation(m, concs, params))
end
