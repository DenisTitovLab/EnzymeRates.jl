using LinearAlgebra
using Random

"""
Independent reference: compute QSSA rate using Laplacian cofactor method.
Works directly with EnzymeMechanism type parameters.
"""
function reference_qssa(m::EnzymeMechanism{Species, Reactions}, params::NamedTuple, concs::NamedTuple) where {Species, Reactions}
    enzs = enzyme_forms(m)
    n = length(enzs)
    enz_names = Tuple(e[1] for e in enzs)
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

"""
    _unicyclic_denominator(f, r)

QSSA denominator for a unicyclic network with forward rates `f[i]`
and reverse rates `r[i]` around the cycle (1-indexed; cyclic).
"""
function _unicyclic_denominator(f::AbstractVector, r::AbstractVector)
    n = length(f)
    @assert length(r) == n
    T = promote_type(eltype(f), eltype(r))
    D = zero(T)

    for i in 1:n
        Ti = zero(T)
        for m in 0:n-1
            pr = one(T)
            for j in 0:m-1
                pr *= r[mod1(i + j, n)]
            end
            pf = one(T)
            for j in (m+1):(n-1)
                pf *= f[mod1(i + j, n)]
            end
            Ti += pr * pf
        end
        D += Ti
    end

    return D
end

"""
    _unicyclic_flux(f, r, Et) -> v

Closed-form QSSA flux for a unicyclic enzyme-state cycle.
v = Et*(prod(f)-prod(r)) / sum(tree_weights)
"""
function _unicyclic_flux(f::AbstractVector, r::AbstractVector, Et)
    pf = prod(f)
    pr = prod(r)
    D  = _unicyclic_denominator(f, r)
    return Et * (pf - pr) / D
end


function random_params_concs(m, met_names::Vector{Symbol}; rng=Random.default_rng())
    ns = n_steps(m)
    param_keys = Symbol[]
    param_vals = Float64[]
    for i in 1:ns
        push!(param_keys, Symbol("k$(i)f"))
        push!(param_vals, 0.1 + 9.9 * rand(rng))
        push!(param_keys, Symbol("k$(i)r"))
        push!(param_vals, 0.1 + 9.9 * rand(rng))
    end
    push!(param_keys, :E_total)
    push!(param_vals, 0.1 + 9.9 * rand(rng))
    params = NamedTuple{Tuple(param_keys)}(Tuple(param_vals))

    conc_vals = [0.1 + 9.9 * rand(rng) for _ in met_names]
    concs = NamedTuple{Tuple(met_names)}(Tuple(conc_vals))

    return params, concs
end

"""
Test that `rate_equation` is non-allocating and fast for the given mechanism.
Must be a standalone function to avoid @testset closure boxing.
"""
function test_rate_equation_performance(m, params, concs)
    rate_equation(m, params, concs) # warmup/compile
    allocs = @allocated rate_equation(m, params, concs)
    t = @elapsed for _ in 1:10_000; rate_equation(m, params, concs); end
    return allocs, t / 10_000
end

"""
Build a NamedTuple with only independent params + Keq + E_total,
given all_params (with all k's + E_total) and a Keq value.
"""
function make_independent_params(m, all_params, Keq)
    indep = independent_parameters(m)
    keys_out = (indep..., :Keq, :E_total)
    vals_out = Tuple(k == :Keq ? Keq : k == :E_total ? all_params.E_total : all_params[k] for k in keys_out)
    return NamedTuple{keys_out}(vals_out)
end

"""
Compute all k values from independent params + Keq by evaluating dependent parameter expressions.
Returns a NamedTuple with all k's + E_total.
"""
function compute_all_params(m, new_params)
    ns = n_steps(m)
    dep = dependent_parameters(m)
    # Build all k values: start with independent ones from new_params, add dependent computed ones
    all_keys = Symbol[]
    all_vals = Float64[]
    for i in 1:ns
        push!(all_keys, Symbol("k$(i)f"))
        push!(all_keys, Symbol("k$(i)r"))
    end
    push!(all_keys, :E_total)
    # Evaluate dependent expressions using new_params as the "params" namespace
    dep_dict = Dict{Symbol, Float64}()
    for (sym, expr_str) in dep
        # Replace params.X with actual values
        val = _eval_dep_expr(expr_str, new_params)
        dep_dict[sym] = val
    end
    all_vals = Float64[haskey(dep_dict, k) ? dep_dict[k] :
                       haskey(new_params, k) ? Float64(new_params[k]) :
                       error("Missing parameter $k")
                       for k in all_keys]
    return NamedTuple{Tuple(all_keys)}(Tuple(all_vals))
end

function _eval_dep_expr(expr_str::String, params::NamedTuple)
    # Build let bindings from params
    bindings = ["$(k) = $(params[k])" for k in keys(params)]
    code = "let params = (; $(join(bindings, ", ")))\n  $expr_str\nend"
    return Float64(eval(Meta.parse(code)))
end

"""
Generate random independent params + Keq + E_total for testing.
Also returns all_params (old-style with all k's + E_total) for reference comparison.
"""
function random_independent_params_concs(m, met_names::Vector{Symbol}; rng=Random.default_rng())
    indep = independent_parameters(m)
    # Generate random values for independent params + Keq + E_total
    param_keys = (indep..., :Keq, :E_total)
    param_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in param_keys)
    new_params = NamedTuple{param_keys}(param_vals)

    conc_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in met_names)
    concs = NamedTuple{Tuple(met_names)}(conc_vals)

    all_params = compute_all_params(m, new_params)
    return new_params, concs, all_params
end

function _reference_metabolite(m)
    subs = substrates(m)
    isempty(subs) && error("No substrate found in mechanism")
    name = subs[1][1]
    coeff = -count(s -> s[1] == name, subs)
    return name, coeff
end

# ── Mechanism constructors ──────────────────────────────────────────────────
# Each returns (mechanism, met_names::Vector{Symbol})

function make_uni_uni()
    species = (
        ((:S, ((:C, 1),)),), ((:P, ((:C, 1),)),), (),
        ((:E, ()), (:ES, ((:C, 1),))),
    )
    rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)))
    EnzymeMechanism(species, rxns), [:S, :P]
end

function make_seq_unibi()
    species = (
        ((:S1, ((:C, 1), (:H, 1))),),
        ((:P1, ((:C, 1),)), (:P2, ((:H, 1),))),
        (),
        ((:E, ()), (:ES1, ((:C, 1), (:H, 1))), (:EP1P2, ((:C, 1), (:H, 1))), (:EP2, ((:H, 1),))),
    )
    rxns = (
        ((:E, :S1), (:ES1,)), ((:ES1,), (:EP1P2,)),
        ((:EP1P2,), (:EP2, :P1)), ((:EP2,), (:E, :P2)),
    )
    EnzymeMechanism(species, rxns), [:S1, :P1, :P2]
end

function make_seq_biuni()
    species = (
        ((:S1, ((:C, 1),)), (:S2, ((:H, 1),))),
        ((:P1, ((:C, 1), (:H, 1))),),
        (),
        ((:E, ()), (:ES1, ((:C, 1),)), (:ES1S2, ((:C, 1), (:H, 1))), (:EP1, ((:C, 1), (:H, 1)))),
    )
    rxns = (
        ((:E, :S1), (:ES1,)), ((:ES1, :S2), (:ES1S2,)),
        ((:ES1S2,), (:EP1,)), ((:EP1,), (:E, :P1)),
    )
    EnzymeMechanism(species, rxns), [:S1, :S2, :P1]
end

function make_seq_bibi()
    species = (
        ((:S1, ((:C, 1),)), (:S2, ((:H, 1),))),
        ((:P1, ((:C, 1),)), (:P2, ((:H, 1),))),
        (),
        ((:E, ()), (:ES1, ((:C, 1),)), (:ES1S2, ((:C, 1), (:H, 1))),
         (:EP1P2, ((:C, 1), (:H, 1))), (:EP2, ((:H, 1),))),
    )
    rxns = (
        ((:E, :S1), (:ES1,)), ((:ES1, :S2), (:ES1S2,)),
        ((:ES1S2,), (:EP1P2,)), ((:EP1P2,), (:EP2, :P1)),
        ((:EP2,), (:E, :P2)),
    )
    EnzymeMechanism(species, rxns), [:S1, :S2, :P1, :P2]
end

function make_seq_biter()
    species = (
        ((:S1, ((:C, 1),)), (:S2, ((:H, 1), (:N, 1)))),
        ((:P1, ((:C, 1),)), (:P2, ((:H, 1),)), (:P3, ((:N, 1),))),
        (),
        ((:E, ()), (:ES1, ((:C, 1),)), (:ES1S2, ((:C, 1), (:H, 1), (:N, 1))),
         (:EP1P2P3, ((:C, 1), (:H, 1), (:N, 1))), (:EP2P3, ((:H, 1), (:N, 1))), (:EP3, ((:N, 1),))),
    )
    rxns = (
        ((:E, :S1), (:ES1,)), ((:ES1, :S2), (:ES1S2,)),
        ((:ES1S2,), (:EP1P2P3,)), ((:EP1P2P3,), (:EP2P3, :P1)),
        ((:EP2P3,), (:EP3, :P2)), ((:EP3,), (:E, :P3)),
    )
    EnzymeMechanism(species, rxns), [:S1, :S2, :P1, :P2, :P3]
end

function make_seq_terbi()
    species = (
        ((:S1, ((:C, 1),)), (:S2, ((:H, 1),)), (:S3, ((:N, 1),))),
        ((:P1, ((:C, 1), (:H, 1))), (:P2, ((:N, 1),))),
        (),
        ((:E, ()), (:ES1, ((:C, 1),)), (:ES1S2, ((:C, 1), (:H, 1))),
         (:ES1S2S3, ((:C, 1), (:H, 1), (:N, 1))), (:EP1P2, ((:C, 1), (:H, 1), (:N, 1))), (:EP2, ((:N, 1),))),
    )
    rxns = (
        ((:E, :S1), (:ES1,)), ((:ES1, :S2), (:ES1S2,)),
        ((:ES1S2, :S3), (:ES1S2S3,)), ((:ES1S2S3,), (:EP1P2,)),
        ((:EP1P2,), (:EP2, :P1)), ((:EP2,), (:E, :P2)),
    )
    EnzymeMechanism(species, rxns), [:S1, :S2, :S3, :P1, :P2]
end

function make_seq_terter()
    species = (
        ((:S1, ((:C, 1),)), (:S2, ((:H, 1),)), (:S3, ((:N, 1),))),
        ((:P1, ((:C, 1),)), (:P2, ((:H, 1),)), (:P3, ((:N, 1),))),
        (),
        ((:E, ()), (:ES1, ((:C, 1),)), (:ES1S2, ((:C, 1), (:H, 1))),
         (:ES1S2S3, ((:C, 1), (:H, 1), (:N, 1))), (:EP1P2P3, ((:C, 1), (:H, 1), (:N, 1))),
         (:EP2P3, ((:H, 1), (:N, 1))), (:EP3, ((:N, 1),))),
    )
    rxns = (
        ((:E, :S1), (:ES1,)), ((:ES1, :S2), (:ES1S2,)),
        ((:ES1S2, :S3), (:ES1S2S3,)), ((:ES1S2S3,), (:EP1P2P3,)),
        ((:EP1P2P3,), (:EP2P3, :P1)), ((:EP2P3,), (:EP3, :P2)),
        ((:EP3,), (:E, :P3)),
    )
    EnzymeMechanism(species, rxns), [:S1, :S2, :S3, :P1, :P2, :P3]
end

function make_pingpong_bibi()
    species = (
        ((:A, ((:C, 2), (:N, 1))), (:B, ((:C, 3),))),
        ((:P, ((:C, 2),)), (:Q, ((:C, 3), (:N, 1)))),
        (),
        ((:E, ()), (:EA, ((:C, 2), (:N, 1))), (:FP, ((:C, 2), (:N, 1))),
         (:F, ((:N, 1),)), (:FB, ((:C, 3), (:N, 1))), (:EQ, ((:C, 3), (:N, 1)))),
    )
    rxns = (
        ((:E, :A), (:EA,)), ((:EA,), (:FP,)), ((:FP,), (:F, :P)),
        ((:F, :B), (:FB,)), ((:FB,), (:EQ,)), ((:EQ,), (:E, :Q)),
    )
    EnzymeMechanism(species, rxns), [:A, :P, :B, :Q]
end

function make_random_bibi()
    species = (
        ((:A, ((:C, 1),)), (:B, ((:N, 1),))),
        ((:P, ((:C, 1),)), (:Q, ((:N, 1),))),
        (),
        ((:E, ()), (:EA, ((:C, 1),)), (:EB, ((:N, 1),)),
         (:EAB, ((:C, 1), (:N, 1))), (:EPQ, ((:C, 1), (:N, 1))), (:EQ, ((:N, 1),))),
    )
    rxns = (
        ((:E, :A), (:EA,)), ((:E, :B), (:EB,)),
        ((:EA, :B), (:EAB,)), ((:EB, :A), (:EAB,)),
        ((:EAB,), (:EPQ,)), ((:EPQ,), (:EQ, :P)), ((:EQ,), (:E, :Q)),
    )
    EnzymeMechanism(species, rxns), [:A, :B, :P, :Q]
end

function make_doubly_branched()
    species = (
        ((:S, ((:C, 1),)),), ((:P, ((:C, 1),)),), (),
        ((:E, ()), (:EA, ((:C, 1),)), (:EB, ((:C, 1),)), (:EC, ((:C, 1),))),
    )
    rxns = (
        ((:E, :S), (:EA,)), ((:EA,), (:EB,)), ((:EA,), (:EC,)),
        ((:EB,), (:E, :P)), ((:EC,), (:E, :P)),
    )
    EnzymeMechanism(species, rxns), [:S, :P]
end

# ── Shared parametric test helpers ──────────────────────────────────────────

"""Run reference_qssa comparison for n_trials random param sets."""
function test_reference_comparison(m, met_names; n=10, seed=42, rtol=1e-8)
    rng = Random.MersenneTwister(seed)
    for _ in 1:n
        new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
        @test rate_equation(m, new_params, concs) ≈ reference_qssa(m, all_params, concs) rtol=rtol
    end
end

"""Test zero-allocation performance."""
function test_performance(m, met_names; seed=42)
    rng = Random.MersenneTwister(seed)
    params, concs, _ = random_independent_params_concs(m, met_names; rng=rng)
    allocs, t = test_rate_equation_performance(m, params, concs)
    @test allocs == 0
    @test t < 100e-9
end

"""Test Haldane: rate=0 at equilibrium for a mechanism with given substrate/product structure."""
function test_haldane_equilibrium(m, met_names; seed=42)
    rng = Random.MersenneTwister(seed)
    new_params, _, _ = random_independent_params_concs(m, met_names; rng=rng)
    Keq = new_params.Keq
    n_subs = length(substrates(m))
    n_prods = length(products(m))
    # Build equilibrium concentrations: prod(P_i) / prod(S_i) = Keq
    # Set all substrates to 1.0, distribute Keq^(1/n_prods) across products
    sub_names = [s[1] for s in substrates(m)]
    prod_names = [p[1] for p in products(m)]
    eq_vals = Dict{Symbol,Float64}()
    for s in sub_names; eq_vals[s] = 1.0; end
    p_each = Keq^(1.0 / n_prods)
    for p in prod_names; eq_vals[p] = p_each; end
    eq_concs = NamedTuple{Tuple(met_names)}(Tuple(eq_vals[s] for s in met_names))
    v_eq = rate_equation(m, new_params, eq_concs)
    @test abs(v_eq) < 1e-10
end
