# ABOUTME: n=1 two-conformation mass-action ground truth for allosteric MWC rate equations.
# ABOUTME: Small explicit networks solved as a linear steady state; the acceptance gate for the normalization fix.
using Test, EnzymeRates, LinearAlgebra, Random
const ER = EnzymeRates

# ── The general solver ──────────────────────────────────────────────────────
# Build the pseudo-first-order rate matrix from a directed edge list, impose
# mass conservation Σc = E_total by replacing the first row, solve the linear
# steady state, and return the net catalytic flux Σ(kf·c[reactant] − kr·c[product])
# over the catalytic (product-forming) edges — the reaction velocity.
#
# `species`  :: Vector{Symbol}                            (species[1] carries the conservation row)
# `edges`    :: Vector of (from::Symbol, to::Symbol, rate::Float64)   directed, rate = pseudo-first-order
# `cat_edges`:: Vector of (reactant::Symbol, product::Symbol, kf::Float64, kr::Float64)
"Net catalytic flux at steady state for an explicit two-conformation network."
function mwc_ground_truth_flux(species, edges, cat_edges, Etot)
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
# S binds :OnlyA (only the active conformation), catalysis :EqualAI, P :EqualAI.
# Forms: E_A, ES_A, EP_A (active) and E_I, EP_I (inactive) plus the dead-end
# ES_I, reached only by inactive catalysis (no inactive S binding, no ES flip).
# Fast RE bindings and flips use FAST; catalysis is O(1). Detailed-balance flip
# ratio [X_I]/[X_A] = L·∏(K_A_i/K_I_i); every present flip here carries an
# :EqualAI ligand (or none), so the ratio is L.
function uni_onlyA_flux(KA, KP, k; L, Keq, S, P, FAST=1e7)
    kr = k * KP / (Keq * KA)
    species = [:E_A, :ES_A, :EP_A, :E_I, :EP_I, :ES_I]
    edges = [
        (:E_A, :ES_A, FAST * S / KA), (:ES_A, :E_A, FAST),   # active S binding (RE)
        (:E_A, :EP_A, FAST * P / KP), (:EP_A, :E_A, FAST),   # active P binding (RE)
        (:E_I, :EP_I, FAST * P / KP), (:EP_I, :E_I, FAST),   # inactive P binding (EqualAI, RE)
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),          # free-enzyme flip, ratio L
        (:EP_A, :EP_I, FAST * L), (:EP_I, :EP_A, FAST),      # EP flip, ratio L
        (:ES_A, :EP_A, k), (:EP_A, :ES_A, kr),              # active catalysis (SS)
        (:ES_I, :EP_I, k), (:EP_I, :ES_I, kr),              # inactive catalysis (SS, dead-end)
    ]
    cat_edges = [(:ES_A, :EP_A, k, kr), (:ES_I, :EP_I, k, kr)]
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# ── Uni-uni all-:EqualAI network (both conformations identical) ──────────────
# Same topology, but S also binds the inactive conformation with K_I = K_A, and
# the ES flip is present. With K_I = K_A everywhere the flip ratio is L on every
# edge, so the two conformations are indistinguishable and the total catalytic
# flux is the non-allosteric rate, independent of L.
function uni_equalAI_flux(KA, KP, k, L, Keq, S, P; FAST=1e7)
    kr = k * KP / (Keq * KA)
    species = [:E_A, :ES_A, :EP_A, :E_I, :ES_I, :EP_I]
    edges = [
        (:E_A, :ES_A, FAST * S / KA), (:ES_A, :E_A, FAST),   # active S binding
        (:E_A, :EP_A, FAST * P / KP), (:EP_A, :E_A, FAST),   # active P binding
        (:E_I, :ES_I, FAST * S / KA), (:ES_I, :E_I, FAST),   # inactive S binding (K_I = K_A)
        (:E_I, :EP_I, FAST * P / KP), (:EP_I, :E_I, FAST),   # inactive P binding (K_I = K_A)
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),          # free-enzyme flip, ratio L
        (:ES_A, :ES_I, FAST * L), (:ES_I, :ES_A, FAST),      # ES flip, ratio L
        (:EP_A, :EP_I, FAST * L), (:EP_I, :EP_A, FAST),      # EP flip, ratio L
        (:ES_A, :EP_A, k), (:EP_A, :ES_A, kr),              # active catalysis
        (:ES_I, :EP_I, k), (:EP_I, :ES_I, kr),              # inactive catalysis
    ]
    cat_edges = [(:ES_A, :EP_A, k, kr), (:ES_I, :EP_I, k, kr)]
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# ── Self-validation: the harness physics, independent of any derivation code ─
@testset "ground-truth harness self-validation" begin
    KA, KP, k, L, Keq, S, P = 1.3, 0.9, 2.1, 0.7, 3.0, 1.1, 0.6
    kr = k * KP / (Keq * KA)
    nonallo = (k * S / KA - kr * P / KP) / (1 + S / KA + P / KP)

    # (a) L = 0 : the inactive conformation is unpopulated, so the :OnlyA flux
    #     reduces to the single-conformation (non-allosteric) uni-uni rate.
    f0 = uni_onlyA_flux(KA, KP, k, L=0.0, Keq=Keq, S=S, P=P)
    @test isapprox(f0, nonallo; rtol=1e-4)

    # (b) all-:EqualAI : both conformations are identical, so the flux is the
    #     base rate and is independent of L.
    fa = uni_equalAI_flux(KA, KP, k, L, Keq, S, P)
    @test isapprox(fa, nonallo; rtol=1e-4)
    @test isapprox(fa, uni_equalAI_flux(KA, KP, k, 5.0, Keq, S, P); rtol=1e-4)
end

# ── The gate: current :OnlyA derivation vs ground truth ──────────────────────
# RED until Task 3. The current derivation combines conformations as Q_A + L·Q_I
# and leaks a bare catalytic rate constant into the L-term when the inactive
# graph fragments; this @test MUST fail against the current code. Task 3's
# free-enzyme cross-weighting fix turns it green (then this stays as the
# regression guard).
@testset "current OnlyA derivation vs ground truth (RED until Task 3)" begin
    onlyA = @allosteric_mechanism begin
        substrates: S ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + S ⇌ E(S)     :: OnlyA
            E(S) <--> E(P)   :: EqualAI
            E + P ⇌ E(P)     :: EqualAI
        end
    end
    fp = ER.fitted_params(onlyA)          # (:K_P_E, :K_A_S_E, :k_ES_to_EP, :L)
    @test fp == (:K_P_E, :K_A_S_E, :k_ES_to_EP, :L)

    rng = MersenneTwister(20260713)
    for _ in 1:5
        KA = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng); k = 0.5 + 2rand(rng)
        L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        S = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map fitted_params -> ground-truth params: K_A_S_E=KA, K_P_E=KP, k_ES_to_EP=k.
        prm = NamedTuple{(fp..., :Keq, :E_total)}((KP, KA, k, L, Keq, 1.0))
        v_code = real(ER.rate_equation(onlyA, (S=S, P=P), prm))
        v_gt = uni_onlyA_flux(KA, KP, k, L=L, Keq=Keq, S=S, P=P)
        @test isapprox(v_code, v_gt; rtol=1e-4)
    end
end
