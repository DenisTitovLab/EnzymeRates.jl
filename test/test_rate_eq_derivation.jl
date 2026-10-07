# ABOUTME: Tests for rate-equation derivation (King-Altman/Cha, mass-action oracle),
# ABOUTME: the parameters API, kcat/rescaling, rate-equation strings, and the perf gate.

using OrdinaryDiffEqFIRK

include("reference/allosteric_undefvar_reproducers.jl")

const ER = EnzymeRates

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
    # Walk the steps of the lifted mechanism, so this helper follows the compiled
    # mechanism representation.
    flat = _testhelper_flat_steps(m)
    enz_names = _testhelper_enzyme_forms(flat)
    n = length(enz_names)
    name_to_idx = Dict(nm => i for (i, nm) in enumerate(enz_names))

    ref_name, nu_ref = _reference_metabolite(m)

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

# ── Brute-force mass-action oracle (no Cha segments) ────────────────────────
# A consistent parameter point: free energies g for every enzyme form and
# metabolite make each step's association constant
#   Ka = exp(g_from + Σ g_consumed − g_to − Σ g_released),
# so every cycle multiplies to Keq^n by construction. The steps of a kinetic
# group share their constants: the form free energies are projected onto
# g_from − g_to = g_from₁ − g_to₁ for every step after its group's first (the
# metabolites cancel, since a group's steps take up and give off the same ones),
# and a steady-state group draws one forward rate.
function _testhelper_consistent_point(em, rng)
    mech = ER.Mechanism(em)
    flat = ER._flat_steps(mech)
    forms = unique(vcat([[ER.from_species(s), ER.to_species(s)] for (s, _) in flat]...))
    g = Dict{Any, BigFloat}(f => 4 * big(rand(rng)) - 2 for f in forms)
    basis = Vector{BigFloat}[]
    for (s, gi) in flat
        s1 = first(ER.steps(mech)[gi])
        s == s1 && continue
        v = BigFloat[(f == ER.from_species(s)) - (f == ER.to_species(s)) -
                     (f == ER.from_species(s1)) + (f == ER.to_species(s1)) for f in forms]
        for b in basis; v -= (b' * v) * b; end
        norm_v = sqrt(v' * v)
        norm_v > big(10)^-30 && push!(basis, v / norm_v)
    end
    x = BigFloat[g[f] for f in forms]
    for b in basis; x -= (b' * x) * b; end
    for (f, xf) in zip(forms, x); g[f] = xf; end
    rxn = ER.reaction(mech)
    mets = vcat([ER.name(x) for x in ER.substrates(rxn)],
                [ER.name(x) for x in ER.products(rxn)],
                [ER.name(ER.regulator(r)) for r in ER.regulators(rxn)])
    for x in mets; g[x] = 4 * big(rand(rng)) - 2; end
    Keq = exp(sum(g[ER.name(x)] for x in ER.substrates(rxn)) -
              sum(g[ER.name(x)] for x in ER.products(rxn)))
    steps = []
    vals = Dict{Symbol, BigFloat}()
    sp = ER._step_parameters(mech)
    kf_of_group = Dict{Int, BigFloat}()
    for (idx, (s, gi)) in enumerate(flat)
        Ka = exp(g[ER.from_species(s)] +
                 sum((g[ER.name(m)] for m in ER.consumed(s)); init = big(0)) -
                 g[ER.to_species(s)] -
                 sum((g[ER.name(m)] for m in ER.released(s)); init = big(0)))
        if ER.is_equilibrium(s)
            p = sp[idx][1]
            vals[ER.name(p, mech)] = p isa ER.Kd ? 1 / Ka : Ka
            push!(steps, (s, Ka, nothing))
        else
            kf = get!(() -> exp(4 * big(rand(rng)) - 2), kf_of_group, gi)
            vals[ER.name(sp[idx][1], mech)] = kf
            vals[ER.name(sp[idx][2], mech)] = kf / Ka
            push!(steps, (s, Ka, kf))
        end
    end
    (; steps, vals, Keq)
end

# v per unit enzyme from the full linear steady state over every form; RE steps
# run at rate scale Λ. v = net consumption of the first substrate.
function _testhelper_mass_action_rate(em, point, concs::Dict{Symbol, BigFloat};
                                      Λ = big(10)^60)
    setprecision(BigFloat, 512) do
        mech = ER.Mechanism(em)
        forms = unique(vcat([[ER.from_species(s), ER.to_species(s)]
                             for (s, _, _) in point.steps]...))
        idx = Dict(f => i for (i, f) in enumerate(forms))
        n = length(forms)
        M = zeros(BigFloat, n, n)
        rates = []
        for (s, Ka, kf) in point.steps
            f_ = kf === nothing ? Λ : kf
            r_ = kf === nothing ? Λ / Ka : kf / Ka
            rf = f_ * prod((concs[ER.name(m)] for m in ER.consumed(s)); init = big(1))
            rr = r_ * prod((concs[ER.name(m)] for m in ER.released(s)); init = big(1))
            a, b = idx[ER.from_species(s)], idx[ER.to_species(s)]
            M[b, a] += rf; M[a, a] -= rf; M[a, b] += rr; M[b, b] -= rr
            push!(rates, (s, a, b, rf, rr))
        end
        M[n, :] .= 1
        rhs = zeros(BigFloat, n); rhs[n] = 1
        x = M \ rhs
        first_sub = ER.name(first(ER.substrates(ER.reaction(mech))))
        sum(((count(m -> ER.name(m) == first_sub, ER.consumed(s)) -
              count(m -> ER.name(m) == first_sub, ER.released(s))) * (rf * x[a] - rr * x[b])
             for (s, a, b, rf, rr) in rates); init = big(0))
    end
end

function _testhelper_check_against_mass_action(em; n_points = 3, seed = 1)
    rng = Random.MersenneTwister(seed)
    fitted = ER.fitted_params(em)
    mets = ER.metabolites(em)
    for _ in 1:n_points
        point = _testhelper_consistent_point(em, rng)
        concs = Dict{Symbol, BigFloat}(x => exp(3 * big(rand(rng)) - 1.5) for x in mets)
        params = NamedTuple{(fitted..., :Keq, :E_total)}(
            (Float64.(getindex.(Ref(point.vals), fitted))..., Float64(point.Keq), 1.0))
        concs_nt = NamedTuple{Tuple(mets)}(Tuple(Float64(concs[x]) for x in mets))
        v_pkg = ER.rate_equation(em, concs_nt, params)
        v_ref = _testhelper_mass_action_rate(em, point, concs)
        @test isapprox(v_pkg, Float64(v_ref); rtol = 1e-8)
    end
end

# ── Parameter generation helpers ────────────────────────────────────────────

"""
Per-step flat index for positional oracle naming. `source_steps === nothing`
keys each canonical step on its canonical stored flat position (group-major) —
self-consistent for the ODE/QSSA cross-checks. When the as-written source-order
step groups are supplied, each canonical step instead keys on the SOURCE flat
index of its structurally-matching as-written step, so positional textbook
oracles (numbered k1,k2,… in source order) line up after the constructor
canonicalizes step order.
"""
function _positional_flat_idx(mech, source_steps)
    if source_steps === nothing
        flat_idx = Vector{Vector{Int}}()
        pos = 0
        for group in EnzymeRates.steps(mech)
            idxs = Int[]
            for _ in group
                pos += 1
                push!(idxs, pos)
            end
            push!(flat_idx, idxs)
        end
        return flat_idx
    end
    # Match canonical steps to as-written steps by structural `Step ==`. The
    # source groups are direction-canonicalized so their stored direction matches
    # the mechanism's; group/within-group order is preserved, so flat position in
    # `src_flat` IS the as-written step index the oracle numbers k1,k2,….
    src_canon = EnzymeRates._canonicalize_step_directions(
        EnzymeRates.reaction(mech), source_steps)
    src_flat = EnzymeRates.Step[s for g in src_canon for s in g]
    used = falses(length(src_flat))
    flat_idx = Vector{Vector{Int}}()
    for group in EnzymeRates.steps(mech)
        idxs = Int[]
        for s in group
            j = findfirst(k -> !used[k] && src_flat[k] == s, eachindex(src_flat))
            j === nothing && error(
                "positional_params bridge: canonical step $s has no as-written match")
            used[j] = true
            push!(idxs, j)
        end
        push!(flat_idx, idxs)
    end
    flat_idx
end

"""
As-written index of canonical regulatory site `site` (at canonical position
`pos`). `source_reg_sites === nothing` returns the canonical position; otherwise
matches by ligand-name set so positional Kreg names (`:K_<lig>_reg{site}`) keep
the oracle's source site numbering after canonical site reordering.
"""
function _positional_site_idx(site, pos, source_reg_sites)
    source_reg_sites === nothing && return pos
    want = Set(EnzymeRates.name(l) for l in EnzymeRates.ligands(site))
    j = findfirst(source_reg_sites) do ss
        Set(EnzymeRates.name(l) for l in EnzymeRates.ligands(ss)) == want
    end
    j === nothing ? pos : j
end

"""
Re-key a structural-named parameter NamedTuple to the **per-step** positional
names (`K1`, `k1f`, `k1r`, …) that hand-derived analytical oracles destructure.
Walks every kinetic group; for allosteric mechanisms also emits the T-state
positional variants (`K1_T`, `k1f_T`, …) for :NonequalAI groups. Kreg
structural keys (`:K_A_<lig>reg`, `:K_I_<lig>reg`, `:K_<lig>reg`) are mapped
to positional keys (`:K_<lig>_reg{site}`, `:K_<lig>_T_reg{site}`). Any
remaining keys (L, Keq, E_total, etc.) are passed through unchanged.
Permanent test utility — oracles are inherently positional and index by
flat-iteration order. `source_steps`/`source_reg_sites` (the as-written orders)
bridge the oracle's source numbering to canonical stored order; omitting them
keeps the canonical numbering the ODE/QSSA cross-checks rely on.
"""
function positional_params(m, nt::NamedTuple;
                           source_steps=nothing, source_reg_sites=nothing)
    mech = m isa EnzymeRates.Mechanism             ? m :
           m isa EnzymeRates.AllostericMechanism   ? m :
           m isa EnzymeRates.AllostericEnzymeMechanism ?
               EnzymeRates.AllostericMechanism(m) :
           EnzymeRates.Mechanism(m)
    is_allo = mech isa EnzymeRates.AllostericMechanism
    names = Symbol[]
    vals  = Any[]

    flat_idx = _positional_flat_idx(mech, source_steps)

    fes = EnzymeRates._free_enz_set(mech)
    for (g, group) in enumerate(EnzymeRates.steps(mech))
        gidx = flat_idx[g]
        rep = EnzymeRates._group_rep(group, fes)
        cat_st = is_allo ? EnzymeRates.cat_allo_state(mech, g) : :None
        # Determine which active-branch state token to pass to name()
        act_st = (cat_st === :EqualAI || cat_st === :None) ? cat_st : :A
        # Inactive branch exists for :NonequalAI groups only
        has_inactive = is_allo && cat_st === :NonequalAI

        if EnzymeRates.is_equilibrium(rep)
            # RE step: binding → Kd; iso → Kiso
            act_key = if EnzymeRates.is_binding(rep)
                EnzymeRates.name(EnzymeRates.Kd(rep, act_st), mech)
            else
                EnzymeRates.name(EnzymeRates.Kiso(rep, act_st), mech)
            end
            if haskey(nt, act_key)
                for idx in gidx
                    push!(names, Symbol("K", idx))
                    push!(vals,  nt[act_key])
                end
            end
            if has_inactive
                ina_key = if EnzymeRates.is_binding(rep)
                    EnzymeRates.name(EnzymeRates.Kd(rep, :I), mech)
                else
                    EnzymeRates.name(EnzymeRates.Kiso(rep, :I), mech)
                end
                if haskey(nt, ina_key)
                    for idx in gidx
                        push!(names, Symbol("K", idx, "_T"))
                        push!(vals,  nt[ina_key])
                    end
                end
            end
        else
            # SS step: binding → Kon/Koff; iso → Kfor/Krev
            act_fwd, act_rev = if EnzymeRates.is_binding(rep)
                EnzymeRates.name(EnzymeRates.Kon(rep, act_st), mech),
                EnzymeRates.name(EnzymeRates.Koff(rep, act_st), mech)
            else
                EnzymeRates.name(EnzymeRates.Kfor(rep, act_st), mech),
                EnzymeRates.name(EnzymeRates.Krev(rep, act_st), mech)
            end
            for idx in gidx
                if haskey(nt, act_fwd)
                    push!(names, Symbol("k", idx, "f")); push!(vals, nt[act_fwd])
                end
                if haskey(nt, act_rev)
                    push!(names, Symbol("k", idx, "r")); push!(vals, nt[act_rev])
                end
            end
            if has_inactive
                ina_fwd, ina_rev = if EnzymeRates.is_binding(rep)
                    EnzymeRates.name(EnzymeRates.Kon(rep, :I), mech),
                    EnzymeRates.name(EnzymeRates.Koff(rep, :I), mech)
                else
                    EnzymeRates.name(EnzymeRates.Kfor(rep, :I), mech),
                    EnzymeRates.name(EnzymeRates.Krev(rep, :I), mech)
                end
                for idx in gidx
                    if haskey(nt, ina_fwd)
                        push!(names, Symbol("k", idx, "f_T")); push!(vals, nt[ina_fwd])
                    end
                    if haskey(nt, ina_rev)
                        push!(names, Symbol("k", idx, "r_T")); push!(vals, nt[ina_rev])
                    end
                end
            end
        end
    end

    # Kreg: emit positional :K_<lig>_reg{site} / :K_<lig>_T_reg{site} keys.
    if is_allo
        for (pos, site) in enumerate(EnzymeRates.regulatory_sites(mech))
            site_idx = _positional_site_idx(site, pos, source_reg_sites)
            for (lig, tag) in zip(EnzymeRates.ligands(site),
                                  EnzymeRates.allo_states(site))
                lig_str = String(EnzymeRates.name(lig))
                act_pos = Symbol("K_", lig_str, "_reg", site_idx)
                ina_pos = Symbol("K_", lig_str, "_T_reg", site_idx)
                if tag === :EqualAI
                    # EqualAI emits K_A_<lig>reg as the independent name (with
                    # K_I_<lig>reg as a Haldane-derived dep equal to it). We look
                    # up whichever of the two is present in nt.
                    struct_key_a = EnzymeRates.name(
                        EnzymeRates.Kreg(site, lig, :A), mech)
                    struct_key_i = EnzymeRates.name(
                        EnzymeRates.Kreg(site, lig, :I), mech)
                    struct_key = haskey(nt, struct_key_a) ? struct_key_a : struct_key_i
                    if haskey(nt, struct_key)
                        push!(names, act_pos); push!(vals, nt[struct_key])
                    end
                elseif tag === :OnlyA
                    struct_key = EnzymeRates.name(
                        EnzymeRates.Kreg(site, lig, :A), mech)
                    if haskey(nt, struct_key)
                        push!(names, act_pos); push!(vals, nt[struct_key])
                    end
                elseif tag === :OnlyI
                    struct_key = EnzymeRates.name(
                        EnzymeRates.Kreg(site, lig, :I), mech)
                    if haskey(nt, struct_key)
                        push!(names, ina_pos); push!(vals, nt[struct_key])
                    end
                else  # :NonequalAI
                    act_key = EnzymeRates.name(
                        EnzymeRates.Kreg(site, lig, :A), mech)
                    ina_key = EnzymeRates.name(
                        EnzymeRates.Kreg(site, lig, :I), mech)
                    if haskey(nt, act_key)
                        push!(names, act_pos); push!(vals, nt[act_key])
                    end
                    if haskey(nt, ina_key)
                        push!(names, ina_pos); push!(vals, nt[ina_key])
                    end
                end
            end
        end
    end

    # Pass through any remaining keys (L, Keq, E_total, etc.) unchanged.
    emitted = Set(names)
    for k in keys(nt)
        k ∈ emitted && continue
        push!(names, k)
        push!(vals, nt[k])
    end
    NamedTuple{Tuple(names)}(Tuple(vals))
end

"""
Positional params for the **hand-written analytical oracles**, which fix
`k{idx}f` as the chemically-forward (substrate→product) rate of source step
`idx`. A step that only gives off one product is stored as the binding it reverses
(`E + P → EP` for a plain release, `E + P → ES` for a fused one), so the
package's stored-forward rate (the one `positional_params` puts on `k{idx}f`)
is actually the chemical REVERSE (binding) of the oracle's forward (release
`EP → E + P`). Swap the `k{idx}f`/`k{idx}r` values for those steps so the
oracle's forward keeps its release meaning. (The QSSA / ODE oracles read the
canonical stored direction directly and need the un-swapped
`positional_params`.)
"""
function analytical_oracle_params(m, nt::NamedTuple;
                                  source_steps=nothing, source_reg_sites=nothing)
    mech = m isa EnzymeRates.AllostericEnzymeMechanism ?
               EnzymeRates.AllostericMechanism(m) :
           m isa EnzymeRates.EnzymeMechanism ? EnzymeRates.Mechanism(m) : m
    pos = positional_params(m, nt; source_steps=source_steps,
                            source_reg_sites=source_reg_sites)
    # swap_idxs must live in the SAME index space as `pos`'s positional keys —
    # the as-written source index when bridged, canonical position otherwise.
    flat_idx = _positional_flat_idx(mech, source_steps)
    swap_idxs = Set{Int}()
    for (g, group) in enumerate(EnzymeRates.steps(mech))
        for (within, s) in enumerate(group)
            if EnzymeRates.bound_metabolite(s) isa EnzymeRates.Product
                push!(swap_idxs, flat_idx[g][within])
            end
        end
    end
    isempty(swap_idxs) && return pos
    names = Symbol[]; vals = Any[]
    for (k, v) in pairs(pos)
        ks = String(k)
        swapped = k
        for i in swap_idxs
            if ks == "k$(i)f";      swapped = Symbol("k", i, "r"); break
            elseif ks == "k$(i)r";  swapped = Symbol("k", i, "f"); break
            elseif ks == "k$(i)f_T"; swapped = Symbol("k", i, "r_T"); break
            elseif ks == "k$(i)r_T"; swapped = Symbol("k", i, "f_T"); break
            end
        end
        push!(names, swapped); push!(vals, v)
    end
    NamedTuple{Tuple(names)}(Tuple(vals))
end

"""Check if mechanism has any rapid-equilibrium steps."""
_has_re_steps(m) = any(EnzymeRates.is_equilibrium, _testhelper_flat_steps(m))

"""
Test that `rate_equation` is non-allocating and fast for the given mechanism.
Must be a standalone function to avoid @testset closure boxing.
"""
function test_rate_equation_performance(m, params, concs)
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
Compute all structural-named params (independent + Haldane-derived dependents)
plus Keq + E_total for a mechanism. Returns a NamedTuple with structural keys
(e.g. :K_ES_to_E_S, :k_ES_to_EP) that positional_params can remap to oracle-style
positional names.
"""
function compute_all_params(m, new_params)
    dep, _ = EnzymeRates._dependent_param_exprs(typeof(m))
    dep_vals = (k => Float64(_eval_dep_expr(e, new_params)) for (k, e) in dep)
    merge(new_params, (; dep_vals...))
end

"""
Evaluate a dependent-parameter expression: a Symbol is looked up in `params`, a
number is itself, and a call applies the named `Base` function to its evaluated
arguments.
"""
_eval_dep_expr(x::Symbol, params) = params[x]
_eval_dep_expr(x::Real, params) = x
_eval_dep_expr(e::Expr, params) =
    getfield(Base, e.args[1])((_eval_dep_expr(a, params) for a in e.args[2:end])...)

"""
Generate random independent params + Keq + E_total for testing.
Also returns all_params (the full set of k's + E_total) for reference comparison.
"""
function random_independent_params_concs(
    m, met_names::Vector{Symbol}; rng=Random.default_rng()
)
    # Generate random values for independent params + Keq + E_total
    param_keys = (EnzymeRates.fitted_params(m)..., :Keq, :E_total)
    param_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in param_keys)
    new_params = NamedTuple{param_keys}(param_vals)

    conc_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in met_names)
    concs = NamedTuple{Tuple(met_names)}(conc_vals)

    all_params = compute_all_params(m, new_params)
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
function raw_to_ode_params(m, raw_params)
    mech = m isa EnzymeRates.Mechanism ? m : EnzymeRates.Mechanism(m)
    flat = EnzymeRates._flat_steps(mech)
    # Kinetic-group rename map: maps Wegscheider-equivalent RE binding K names
    # (e.g. K_ERinhS_to_ERinh_S → K_ES_to_E_S) so the param lookup succeeds even when the
    # mechanism has sharing via Wegscheider constraints.
    rename = EnzymeRates._build_wegscheider_rename_map(mech)
    # A canonical RE binding step has a metabolite on LHS (canonical form
    # invariant: all RE binding steps are written `E + S ⇌ ES`).
    is_binding_step = Bool[
        EnzymeRates.is_equilibrium(s) && !isempty(EnzymeRates.consumed(s))
        for (s, _) in flat
    ]
    # Resolve a structural param key through the rename map if not present
    _lookup(k) = haskey(raw_params, k) ? Float64(raw_params[k]) :
                 haskey(rename, k) ? Float64(raw_params[rename[k]]) :
                 error("raw_to_ode_params: missing param $k")
    param_keys = Symbol[]
    param_vals = Float64[]
    for (i, (step, g)) in enumerate(flat)
        rep_step = first(EnzymeRates.steps(mech)[g])
        push!(param_keys, Symbol("k$(i)f"))
        push!(param_keys, Symbol("k$(i)r"))
        if EnzymeRates.is_equilibrium(step)
            # Look up structural K key for the rep step (Kd for binding, Kiso for iso)
            if is_binding_step[i]
                K_key = EnzymeRates.name(EnzymeRates.Kd(rep_step, :None), mech)
                K = _lookup(K_key)
                # Binding step (metabolite on LHS): K = Kd = kr/kf
                push!(param_vals, 1e6)
                push!(param_vals, 1e6 * K)
            else
                K_key = EnzymeRates.name(EnzymeRates.Kiso(rep_step, :None), mech)
                K = _lookup(K_key)
                # RE isomerization (no metabolite): K = Ka = kf/kr
                push!(param_vals, 1e6 * K)
                push!(param_vals, 1e6)
            end
        else
            # SS step: binding → Kon/Koff; iso → Kfor/Krev
            fwd_key, rev_key = if EnzymeRates.is_binding(rep_step)
                EnzymeRates.name(EnzymeRates.Kon(rep_step, :None), mech),
                EnzymeRates.name(EnzymeRates.Koff(rep_step, :None), mech)
            else
                EnzymeRates.name(EnzymeRates.Kfor(rep_step, :None), mech),
                EnzymeRates.name(EnzymeRates.Krev(rep_step, :None), mech)
            end
            push!(param_vals, _lookup(fwd_key))
            push!(param_vals, _lookup(rev_key))
        end
    end
    push!(param_keys, :E_total)
    push!(param_vals, Float64(raw_params[:E_total]))
    return NamedTuple{Tuple(param_keys)}(Tuple(param_vals))
end

function _reference_metabolite(m)
    subs = _testhelper_substrates(m)
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
    flat = _testhelper_flat_steps(m)
    E_total = params.E_total
    enz_names = _testhelper_enzyme_forms(flat)
    n = length(enz_names)
    ref_name, nu_ref = _reference_metabolite(m)

    u0 = zeros(n)
    u0[1] = E_total

    rhs! = build_ode_rhs(m, params, concs)
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
function _eval_rate_string(s, params, concs)
    code = "let params = $params, concs = $concs\n" *
           replace(s, EnzymeRates.ANNOTATION_SUBSTITUTED => "") * "\nend"
    eval(Meta.parse(code))
end

# ── Modular test functions for MechanismTestSpec ────────────────────────────

function test_structure(spec::MechanismTestSpec)
    m = spec.mechanism
    @testset "Structure" begin
        flat = _testhelper_flat_steps(m)
        @test length(_testhelper_enzyme_forms(flat)) == spec.expected_n_states
        @test length(flat) == spec.expected_n_steps
        @test metabolites(m) == Tuple(spec.metabolite_names)
        # Structural parameter names must be injective on the load-bearing
        # paths the fitter consumes (Reduced + fitted_params). A collision
        # would silently shorten the destructured params NamedTuple.
        @test allunique(EnzymeRates.parameters(m))           # Reduced (default)
        @test allunique(EnzymeRates.fitted_params(m))
    end
end

"""Classify a dep expression as Haldane (RHS references Keq), Mirror
(RHS is a single Symbol), or Wegscheider (RHS Expr without Keq)."""
function _classify_dep_expr(expr)
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
rhs_syms(x) = x isa Symbol ? [x] :
    x isa Expr ? reduce(vcat, map(rhs_syms, x.head === :call ? x.args[2:end] : x.args);
                         init=Symbol[]) :
    Symbol[]

"""
The dependent-parameter graph `dep` (with independent parameters `indep`) is sound
iff (a) no symbol is both independent (fitted) and dependent, (b) the dep→dep
dependency graph is acyclic, and (c) every dependent expression's RHS symbol
is a fitted param, another dependent param, or `Keq`.
"""
function _dep_graph_is_sound(dep, indep)
    indepset = Set(indep)
    isempty(intersect(indepset, keys(dep))) || return false

    depset = Set(keys(dep))
    edges = Dict(k => Symbol[s for s in rhs_syms(v) if s in depset] for (k, v) in dep)
    state = Dict{Symbol, Int}()  # 0 = unseen, 1 = on-stack, 2 = done
    ok = true
    function dfs(n)
        get(state, n, 0) == 2 && return
        if get(state, n, 0) == 1
            ok = false
            return
        end
        state[n] = 1
        for m in get(edges, n, Symbol[])
            dfs(m)
            ok || return
        end
        state[n] = 2
    end
    for k in keys(dep)
        dfs(k)
        ok || break
    end
    ok || return false

    known = union(indepset, depset, Set([:Keq]))
    for (_, v) in dep, s in rhs_syms(v)
        s in known || return false
    end
    true
end

function test_constraint_counting(spec::MechanismTestSpec)
    m = spec.mechanism
    @testset "Constraints" begin
        dep_exprs, indep = EnzymeRates._dependent_param_exprs(typeof(m))
        # A dependent (Haldane/Wegscheider-derived) parameter must never appear in
        # the independent set returned as fitted_params. The LDH i-state specs have a
        # shared :EqualAI catalytic reverse rate that is Haldane-dependent in the
        # A-state while the dead I-state references it unpinned — the exact leak the
        # uniform dep-filter closes.
        @test _dep_graph_is_sound(dep_exprs, indep)
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
            # reference_qssa uses positional step-indexed k$(i)f/k$(i)r keys
            pos_params = positional_params(m, all_params)
            isapprox(
                rate_equation(m, concs, new_params),
                reference_qssa(m, pos_params, concs);
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
            p = merge(analytical_oracle_params(
                          m, all_params;
                          source_steps=spec.source_steps,
                          source_reg_sites=spec.source_reg_sites),
                      (Et=Et,))
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

function test_performance(spec::MechanismTestSpec; seed=42)
    m = spec.mechanism
    met_names = spec.metabolite_names
    @testset "Performance" begin
        rng = Random.MersenneTwister(seed)
        params, concs, _ = random_independent_params_concs(m, met_names; rng=rng)
        allocs, t = test_rate_equation_performance(m, params, concs)
        @test allocs == 0
        # 120ns, not the ~tens-of-ns real per-call cost: shared CI runners'
        # best-case timing runs slower than a dedicated box, so 120ns keeps
        # margin without flaking. The 0-alloc bound above is the strict part.
        @test t < 120e-9
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
            # Convert structural params to positional k_if/k_ir for ODE
            pos_params = positional_params(m, all_params)
            ode_params = has_re ?
                raw_to_ode_params(m, all_params) :
                pos_params
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
            new_params, concs, _ =
                random_independent_params_concs(
                    m, met_names; rng=rng)
            isapprox(
                rate_equation(m, concs, new_params),
                _eval_rate_string(s, new_params, concs);
                rtol=1e-10)
        end

        has_num = spec.expected_factored_num !== nothing
        has_denom = spec.expected_factored_denom !== nothing
        if has_num || has_denom
            num_str, denom_str = _extract_num_denom(last(split(s, "\n")))
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

function test_analytical_kcat(spec::MechanismTestSpec; seed=42)
    spec.analytical_kcat_fn === nothing && return
    m = spec.mechanism
    @testset "Analytical kcat" begin
        rng = Random.MersenneTwister(seed)
        params = random_reduced_params(m; rng)
        kcat = EnzymeRates._kcat_forward(m, params)
        # The oracle's positional formula may reference a forward rate that is
        # Haldane-DEPENDENT under the canonical step order (absent from the
        # reduced params), so bridge the FULL param set (dependent values
        # included) to positional names.
        p = merge(analytical_oracle_params(
                      m, compute_all_params(m, params);
                      source_steps=spec.source_steps,
                      source_reg_sites=spec.source_reg_sites),
                  (Et=params.E_total,))
        @test kcat ≈ spec.analytical_kcat_fn(p) rtol=1e-10
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
        norm_custom = rescale_parameter_values(m, params; scale_k_to_kcat=kcat_target)
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
        sub_names = _testhelper_substrates(m)
        prod_names = _testhelper_products(m)
        reg_names = _testhelper_regulators(m)
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
function test_zero_metabolite_finite(spec::MechanismTestSpec)
    m = spec.mechanism
    @testset "Zero-metabolite finiteness" begin
        rng = Random.MersenneTwister(777 + hash(spec.name) % 1000)
        mets = collect(metabolites(m))
        sub_prod = Set{Symbol}(_testhelper_substrates(m))
        union!(sub_prod, _testhelper_products(m))
        params = random_reduced_params(m; rng)
        for zeroed in mets
            cvals = Tuple(n == zeroed ? 0.0 : 0.5 + rand(rng) for n in mets)
            concs = NamedTuple{Tuple(mets)}(cvals)
            v = rate_equation(m, concs, params)
            @test isfinite(v)
            zeroed in sub_prod && @test v != 0.0
        end
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
        test_zero_metabolite_finite(spec)
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
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; regulators: I
            steps: begin
            E + A ⇌ E(A); E(A) + B ⇌ E(A, B); E(A, B) <--> E(Q) + P; E(Q) ⇌ E + Q
            E + I ⇌ E(I); E(Q) + I ⇌ E(I, Q) end end),
        # ping-pong: parallel product-release route via a mixed substrate/product complex
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; regulators: I
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
        # the first substrate also binds as a dead-end inhibitor copy, at rapid
        # equilibrium on E and at steady state on E(Q); both the derivation and
        # the oracle count the first substrate by name
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; regulators: A
            steps: begin
            E + A ⇌ E(A); E(A) + B <--> E(A, B); E(A, B) <--> E(Q) + P; E(Q) <--> E + Q
            E + A::Inh ⇌ E(A::Inh); E(Q) + A::Inh <--> E(A::Inh, Q) end end),
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
        _testhelper_check_against_mass_action(em)
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
        _testhelper_check_against_mass_action(em)
    end
    # canonical orientation: the substrate side is `from`
    s = only(s for g in ER.steps(ER.Mechanism(tc_ss)) for s in g
             if !ER.is_iso(s) && !ER.is_binding(s))
    @test ER.name(ER.from_species(s)) == :EA &&
          ER.consumed(s) == ER.Metabolite[ER.Substrate(:B)]
    @test :k_EA_B_to_EQ_P in ER.parameters(tc_ss, ER.Full)
    @test :k_EQ_P_to_EA_B in ER.parameters(tc_ss, ER.Full)
end

@testset "reversing written steps changes nothing" begin
    rev(s) = ER.Step(ER.to_species(s), ER.from_species(s), ER.released(s), ER.consumed(s),
                     ER.is_equilibrium(s))
    rng = Random.MersenneTwister(7)
    # Every step reversed, a random subset, and every group's steps in reverse order.
    variants(groups) = ([[rev(s) for s in g] for g in groups],
                        [[rand(rng, Bool) ? rev(s) : s for s in g] for g in groups],
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
        for flipped in variants(ER.steps(m))
            m2 = ER.Mechanism(ER.reaction(m), flipped)
            @test m2 == m
            @test ER.rate_equation_string(ER.compile_mechanism(m2)) ==
                  ER.rate_equation_string(em)
        end
    end
    for spec in MECHANISM_TEST_SPECS
        spec.mechanism isa ER.AllostericEnzymeMechanism || continue
        am = ER.AllostericMechanism(spec.mechanism)
        for flipped in variants(ER.steps(am))
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
    onlyA_groups = [g for g in EnzymeRates.kinetic_groups(am)
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
function _undefined_rhs_symbols(s::AbstractString)
    ident = r"[A-Za-z_][A-Za-z0-9_]*"
    defined = Set{Symbol}()
    referenced = Set{Symbol}()
    for raw in split(s, '\n')
        line = replace(raw, EnzymeRates.ANNOTATION_SUBSTITUTED => "")
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
    for em in _LDH_ISTATE_MECHS
        pnames = EnzymeRates.fitted_params(em)
        mets = EnzymeRates.metabolites(em)
        prods = Set(_testhelper_products(em))
        params = merge(NamedTuple{pnames}(ntuple(_ -> 1.3, length(pnames))),
                       (Keq = 20000.0, E_total = 1.0))
        # No UndefVarError: the @generated body compiles and evaluates finite at
        # products = 0 (the kcat evaluation domain).
        concs0 = NamedTuple{mets}(
            ntuple(i -> mets[i] in prods ? 0.0 : 1.5, length(mets)))
        @test isfinite(EnzymeRates.rate_equation(em, concs0, params))
        # DEFINED ⊇ REFERENCED on the rendered transcript.
        @test isempty(_undefined_rhs_symbols(EnzymeRates.rate_equation_string(em)))
    end
end

@testset "allosteric dependent-param graph is sound" begin
    for T in ALLOSTERIC_UNDEFVAR_REPRODUCERS
        @test _dep_graph_is_sound(EnzymeRates._dependent_param_exprs(T)...)
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

# ── Standalone kcat tests ──────────────────────────────────────────────────────

@testset "rate_equation polynomial body uses 2-arg +/* calls" begin
    # The fitter calls rate_equation millions of times per CV fold. The
    # polynomial body emitted by _poly_to_expr (via _nest_binary) MUST
    # have exactly 2 operands per +/* call so LLVM inlines the binary
    # Float64 path; n-ary varargs above ~30 terms boxes the argument
    # tuple and turns 100ns/0B into 1µs/2KB per call.
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

        # No factors and zero Keq: returns Int literal 1
        @test bpe(R(0), Tuple{Symbol, R}[]) === 1

        # Keq times a factor is a valid AST node
        @test bpe(R(1), [(:k1f, R(1))]) isa Union{Int, Symbol, Expr}
    end

end

@testset "a reaction with no net change keeps its binding-square constraint" begin
    # Every product shares a substrate's name, so the net reaction changes no
    # metabolite. The random-order binding square closes one Wegscheider cycle
    # (log Keq exponent 0), and the :OnlyA check reads that cycle.
    m = @enzyme_mechanism begin
        substrates: A, B
        products:   A, B
        steps: begin
            E + A ⇌ E(A)
            E + B ⇌ E(B)
            E(A) + B ⇌ E(A, B)
            E(B) + A ⇌ E(A, B)
        end
    end
    C, rhs = EnzymeRates._thermodynamic_constraints(EnzymeRates.Mechanism(m))
    @test C == [1 -1 1 -1]
    @test rhs == [0]
    @test EnzymeRates._dependent_param_exprs(EnzymeRates.Mechanism(m))[2] ==
          (:K_EAB_to_EB_A, :K_EA_to_E_A, :K_EB_to_E_B)

    # A paired :OnlyA binding balances the cycle; a lone one leaves it unsatisfiable.
    @test @allosteric_mechanism(begin
        substrates: A, B
        products:   A, B
        catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)          :: OnlyA
            E + B ⇌ E(B)          :: EqualAI
            E(A) + B ⇌ E(A, B)    :: EqualAI
            E(B) + A ⇌ E(A, B)    :: OnlyA
        end
    end) isa AllostericEnzymeMechanism
    lone_onlya = ErrorException(
        "AllostericMechanism: an :OnlyA binding (K_EA_to_E_A) leaves a thermodynamic " *
        "(Haldane/Wegscheider) cycle unsatisfiable: the inactive conformation cannot " *
        "close that cycle at finite nonzero affinity. Tag the cycle's chemical step " *
        ":OnlyA, or tag an opposing binding :OnlyA so the affinities diverge together.")
    @test_throws lone_onlya @allosteric_mechanism begin
        substrates: A, B
        products:   A, B
        catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)          :: OnlyA
            E + B ⇌ E(B)          :: EqualAI
            E(A) + B ⇌ E(A, B)    :: EqualAI
            E(B) + A ⇌ E(A, B)    :: EqualAI
        end
    end
end

# ── Large equation compilation regression test ────────────────────────────

@testset "Large equation compilation (<20s)" begin
    # A ping-pong bi-bi with a dead-end inhibitor R1: 13 forms and 13 steps, every
    # step at rapid equilibrium except the A → P isomerization.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: R1
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
        regulators: R
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
        allosteric_regulators: I::OnlyI
        catalytic_steps: begin
            E + F6P <--> E(F6P)      :: EqualAI
            E(F6P) <--> E(F16BP)     :: EqualAI
            E(F16BP) <--> E + F16BP  :: EqualAI
        end
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
        catalytic_steps: begin
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
        catalytic_steps: begin
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
    # T-state numerator branch is elided when t_state_dead (any :OnlyA catalytic group);
    # rate is E_total · num_R / (Q_R^catN + L · Q_T^catN). At large L, the T-state
    # enzyme mass dominates the denominator → rate ∝ 1/(1+L).
    @test rate_T * 1e10 < 100.0    # bounded as L grows

    # Regression: T-state binding K's must be in Kd convention even when
    # `:OnlyA` and `:NonequalAI` catalytic groups coexist. Without the fix,
    # the flat-poly path in _allosteric_num_den_exprs renders T-state K's
    # as `K_T * met` (Ka) instead of `met / K_T` (Kd), silently producing
    # wrong rates whenever a mechanism mixes these two tags. Regression
    # for src/rate_eq_derivation.jl:1395-1396.
    cm_mix, src_mix = @enzyme_mechanism_src begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    # Bind allosteric states to the steps AS WRITTEN: S binding :NonequalAI,
    # catalysis :OnlyA, P binding :NonequalAI.
    m_mix = allo_from_source(
        (cm_mix, src_mix), (2, (:NonequalAI, :OnlyA, :NonequalAI)),
        (((:I,), 2, (:OnlyI,)),))
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
    cm_ro, src_ro = @enzyme_mechanism_src begin
        substrates: A, B
        products:   P, Q
        steps: begin
            E + A <--> E(A)
            E + B <--> E(B)
            E(A) + B <--> E(A, B)
            E(B) + A <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end
    # B binding (step 2) :NonequalAI, rest :EqualAI — bound to the steps AS WRITTEN.
    m_ro = allo_from_source(
        (cm_ro, src_ro),
        (2, (:EqualAI, :NonequalAI, :EqualAI, :EqualAI,
             :EqualAI, :EqualAI, :EqualAI)),
        (((:I,), 2, (:NonequalAI,)),))
    # Overall A + B ⇌ P + Q, so Keq = P·Q/(A·B).
    Keq_ro = 4.0
    A_eq, B_eq = 1.5, 2.0
    P_eq = 3.0
    Q_eq = Keq_ro * A_eq * B_eq / P_eq
    # This degenerate over-parametrized Wegscheider cycle has an ambiguous
    # independent-parameter basis (which koff is Wegscheider-dependent can
    # shift). Generate the param values from the mechanism's actual
    # fitted_params so the Haldane property (zero net rate at equilibrium)
    # is verified regardless of the basis the derivation selects.
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
        allosteric_regulators: R::NonequalAI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)     :: NonequalAI
            E + S ⇌ E(S)     :: NonequalAI
            E(S) <--> E(P)   :: NonequalAI
        end
    end
    actual = rate_equation_string(m_allo)
    expected = raw"""(; K_A_EP_to_E_P, K_A_ES_to_E_S, k_A_ES_to_EP, K_I_EP_to_E_P, K_I_ES_to_E_S, k_I_ES_to_EP, K_A_Rreg, K_I_Rreg, L, Keq, E_total) = params
(; S, P, R) = concs
# Haldane constraints:
k_A_EP_to_ES = (1 / Keq) * K_A_EP_to_E_P * (1 / K_A_ES_to_E_S) * k_A_ES_to_EP
k_I_EP_to_ES = (1 / Keq) * K_I_EP_to_E_P * (1 / K_I_ES_to_E_S) * k_I_ES_to_EP
v = E_total * ((k_A_ES_to_EP * S / K_A_ES_to_E_S - k_A_EP_to_ES * P / K_A_EP_to_E_P) * (1 + P / K_A_EP_to_E_P + S / K_A_ES_to_E_S) * (1 + R / K_A_Rreg) ^ 2 + L * (S * k_I_ES_to_EP / K_I_ES_to_E_S - P * k_I_EP_to_ES / K_I_EP_to_E_P) * (1 + P / K_I_EP_to_E_P + S / K_I_ES_to_E_S) * (1 + R / K_I_Rreg) ^ 2) / ((1 + P / K_A_EP_to_E_P + S / K_A_ES_to_E_S) ^ 2 * (1 + R / K_A_Rreg) ^ 2 + L * (1 + P / K_I_EP_to_E_P + S / K_I_ES_to_E_S) ^ 2 * (1 + R / K_I_Rreg) ^ 2)"""
    @test actual == expected
end

@testset "Parameter-struct allosteric helpers" begin
    cm_src = @enzyme_mechanism_src begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end

    @testset "_cat_params and _kreg_params" begin
        names_of(am, params) = [EnzymeRates.name(p, am) for p in params]
        i_params(am) = [EnzymeRates._cat_params(am, :I); EnzymeRates._kreg_params(am, :I)]
        # :NonequalAI cat group + :NonequalAI reg ligand → both contribute.
        aem = allo_from_source(
            cm_src, (2, (:NonequalAI, :EqualAI, :NonequalAI)),
            (((:R,), 1, (:NonequalAI,)),),
        )
        am = EnzymeRates.AllostericMechanism(aem)
        rendered = names_of(am, i_params(am))

        # Catalytic side: every non-:OnlyA group contributes an I-state
        # parameter (Kd for RE binding, Kfor/Krev for SS).
        @test :K_I_ES_to_E_S in rendered
        @test :k_I_ES_to_EP in rendered
        @test :k_I_EP_to_ES in rendered
        @test :K_I_EP_to_E_P in rendered
        # Regulator side: :R is :NonequalAI → K_I_Rreg appears.
        @test :K_I_Rreg in rendered
        # Active state: the :EqualAI iso takes its shared untagged names; every other
        # group and the regulator take A-state names.
        @test names_of(am, EnzymeRates._cat_params(am, :A)) ==
              [:K_A_EP_to_E_P, :K_A_ES_to_E_S, :k_ES_to_EP, :k_EP_to_ES]
        @test names_of(am, EnzymeRates._kreg_params(am, :A)) == [:K_A_Rreg]

        # Balanced :OnlyA bindings (both S and P) keep the :NonequalAI iso
        # Haldane-valid — a one-sided :OnlyA binding under a non-:OnlyA iso is
        # rejected at construction. Both :OnlyA cat groups and the :OnlyA reg
        # ligand are skipped; the :NonequalAI iso still emits both k params.
        aem_skip = allo_from_source(
            cm_src, (2, (:OnlyA, :NonequalAI, :OnlyA)),
            (((:R,), 1, (:OnlyA,)),),
        )
        am_skip = EnzymeRates.AllostericMechanism(aem_skip)
        rendered_skip = names_of(am_skip, i_params(am_skip))
        @test :K_I_ES_to_E_S ∉ rendered_skip    # :OnlyA cat group skipped
        @test :k_I_ES_to_EP in rendered_skip    # :NonequalAI SS iso emits both
        @test :k_I_EP_to_ES in rendered_skip
        @test :K_I_Rreg ∉ rendered_skip    # :OnlyA reg ligand skipped
        @test names_of(am_skip, EnzymeRates._kreg_params(am_skip, :A)) == [:K_A_Rreg]

        # An :OnlyI reg ligand has an I-state Kreg and no A-state one.
        am_onlyi = EnzymeRates.AllostericMechanism(allo_from_source(
            cm_src, (2, (:NonequalAI, :EqualAI, :NonequalAI)),
            (((:R,), 1, (:OnlyI,)),),
        ))
        @test isempty(EnzymeRates._kreg_params(am_onlyi, :A))
        @test names_of(am_onlyi, EnzymeRates._kreg_params(am_onlyi, :I)) == [:K_I_Rreg]
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
        allosteric_regulators: ATP::OnlyI, F16BP::OnlyA
        catalytic_multiplicity: 4
        catalytic_steps: begin
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
        regulatory_site(multiplicity = 2): begin
            ligands: ATP
        end
        regulatory_site(multiplicity = 4): begin
            ligands: F16BP
        end
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

@testset "allosteric single-symbol Wegscheider tie is absorbed in both states" begin
    # A binds E and E(B) under one shared K while B's two bindings keep their own, so
    # the RE binding square ties B's two K's (K_EAB_to_EA_B = K_EB_to_E_B) in each
    # conformation. Each state's Wegscheider rename folds its tie onto K_EB_to_E_B,
    # and the folded symbols leave fitted_params, so the mechanism derives exactly as
    # the one that groups B's two bindings outright.
    tie = @allosteric_mechanism begin
        substrates: A, B
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
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
        catalytic_steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B)) :: NonequalAI
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B)) :: NonequalAI
            E(A, B) <--> E(P)                  :: NonequalAI
            E(P) ⇌ E + P                       :: NonequalAI
        end
    end
    @test EnzymeRates.fitted_params(tie) ==
        (:K_A_EA_to_E_A, :K_A_EB_to_E_B, :K_A_EP_to_E_P, :k_A_EAB_to_EP,
         :K_I_EA_to_E_A, :K_I_EB_to_E_B, :K_I_EP_to_E_P, :k_I_EAB_to_EP, :L)
    @test rate_equation_string(tie) == rate_equation_string(grouped)
end

@testset "Fix A: dead-inactive-state allosteric body defines all I-state symbols" begin
    # Random-order allosteric bi-bi with an :OnlyA catalytic step → dead inactive
    # state. Verified pre-fix to crash with an UndefVarError on the I-state
    # release rate of A from E(A).
    m = @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 2
        catalytic_steps: begin
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
    @test isfinite(rate_equation(m, concs, params, Reduced))
end

@testset "Fix B: _kcat_forward handles multiple saturating patterns" begin
    # Random-order allosteric bi-bi, live I-state. Verified pre-fix to raise
    # "multiple saturating-substrate kcat components (9 found)".
    m = @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 2
        catalytic_steps: begin
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
    kc = EnzymeRates._kcat_forward(m, merge(pv, (Keq = 1.0,)))
    @test isfinite(kc)
    # Peak-productive-turnover contract: equals the numerical peak forward rate.
    fp = merge(pv, (Keq = 1.0, E_total = 1.0))
    vmax = maximum(rate_equation(m, (A = x, B = y, P = 0.0, Q = 0.0), fp, Reduced)
                   for x in 10.0 .^ (0:1:9), y in 10.0 .^ (0:1:9))
    @test kc ≈ vmax rtol = 1e-3
end

@testset "Fix B: non-allosteric random-order bi-bi kcat = peak (contract guard)" begin
    # Sweeps 40 parameter draws; product-containing King–Altman cross-terms
    # (e.g. A·B·P monomials) are spurious candidates at products=0 and can
    # inflate _kcat_forward by >10× before the product filter is applied.
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
        kc = EnzymeRates._kcat_forward(m, merge(pv, (Keq = 1.0,)))
        fp = merge(pv, (Keq = 1.0, E_total = 1.0))
        vmax = maximum(rate_equation(m, (A = x, B = y, P = 0.0, Q = 0.0), fp, Reduced)
                       for x in 10.0 .^ (0:1:9), y in 10.0 .^ (0:1:9))
        @test kc ≈ vmax rtol = 1e-3
    end
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
