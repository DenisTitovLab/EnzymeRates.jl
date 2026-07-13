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

# ── Multi-:OnlyA bi-uni network (the fragmenting case the fix must get right) ──
# A + B ⇌ P. Both A and B bind :OnlyA (active conformation only); catalysis
# :EqualAI, P :EqualAI. Active forms E_A, EA_A, EAB_A, EP_A. Inactive forms
# E_I, EP_I (RE-connected via P binding) and the dead-end EAB_I, reached only by
# reverse catalysis from EP_I (A,B :OnlyA → no inactive A/B binding, no EAB flip).
# The inactive graph fragments: {E_I, EP_I} is one rapid-equilibrium segment and
# {EAB_I} its own segment behind the steady-state catalytic edge, so the inactive
# free-enzyme spanning-tree weight D_I = k while D_A = 1 — the cross-weighting the
# fix supplies. Flip ratio [X_I]/[X_A] = L·∏(K_A_i/K_I_i): free enzyme and EP
# (:EqualAI or bare) flip with ratio L; an :OnlyA-ligand-bearing state has ratio 0
# (no flip), so EA and EAB never flip.
function multi_onlyA_flux(KA, KB, KP, k; L, Keq, A, B, P, FAST=1e7)
    kr = k * KP / (Keq * KA * KB)
    species = [:E_A, :EA_A, :EAB_A, :EP_A, :E_I, :EP_I, :EAB_I]
    edges = [
        (:E_A, :EA_A, FAST * A / KA), (:EA_A, :E_A, FAST),    # active A binding (RE)
        (:EA_A, :EAB_A, FAST * B / KB), (:EAB_A, :EA_A, FAST),# active B binding (RE)
        (:E_A, :EP_A, FAST * P / KP), (:EP_A, :E_A, FAST),    # active P binding (RE)
        (:E_I, :EP_I, FAST * P / KP), (:EP_I, :E_I, FAST),    # inactive P binding (EqualAI, RE)
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),           # free-enzyme flip, ratio L
        (:EP_A, :EP_I, FAST * L), (:EP_I, :EP_A, FAST),       # EP flip, ratio L
        (:EAB_A, :EP_A, k), (:EP_A, :EAB_A, kr),             # active catalysis (SS)
        (:EAB_I, :EP_I, k), (:EP_I, :EAB_I, kr),             # inactive catalysis (SS, dead-end)
    ]
    cat_edges = [(:EAB_A, :EP_A, k, kr), (:EAB_I, :EP_I, k, kr)]
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# ── Multi-:EqualAI bi-uni network (both conformations identical) ──────────────
# Same A + B ⇌ P topology, but A and B also bind the inactive conformation with
# K_I = K_A, and every form flips (ratio L). With K_I = K_A everywhere the two
# conformations are indistinguishable, so the flux is the non-allosteric bi-uni
# rate, independent of L.
function multi_equalAI_flux(KA, KB, KP, k, L, Keq, A, B, P; FAST=1e7)
    kr = k * KP / (Keq * KA * KB)
    species = [:E_A, :EA_A, :EAB_A, :EP_A, :E_I, :EA_I, :EAB_I, :EP_I]
    edges = [
        (:E_A, :EA_A, FAST * A / KA), (:EA_A, :E_A, FAST),    # active A binding
        (:EA_A, :EAB_A, FAST * B / KB), (:EAB_A, :EA_A, FAST),# active B binding
        (:E_A, :EP_A, FAST * P / KP), (:EP_A, :E_A, FAST),    # active P binding
        (:E_I, :EA_I, FAST * A / KA), (:EA_I, :E_I, FAST),    # inactive A binding (K_I = K_A)
        (:EA_I, :EAB_I, FAST * B / KB), (:EAB_I, :EA_I, FAST),# inactive B binding (K_I = K_A)
        (:E_I, :EP_I, FAST * P / KP), (:EP_I, :E_I, FAST),    # inactive P binding (K_I = K_A)
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),           # free-enzyme flip, ratio L
        (:EA_A, :EA_I, FAST * L), (:EA_I, :EA_A, FAST),       # EA flip, ratio L
        (:EAB_A, :EAB_I, FAST * L), (:EAB_I, :EAB_A, FAST),   # EAB flip, ratio L
        (:EP_A, :EP_I, FAST * L), (:EP_I, :EP_A, FAST),       # EP flip, ratio L
        (:EAB_A, :EP_A, k), (:EP_A, :EAB_A, kr),             # active catalysis
        (:EAB_I, :EP_I, k), (:EP_I, :EAB_I, kr),             # inactive catalysis
    ]
    cat_edges = [(:EAB_A, :EP_A, k, kr), (:EAB_I, :EP_I, k, kr)]
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

# ── Multi-:OnlyA harness self-validation (same physics, bi-uni topology) ─────
@testset "multi-OnlyA ground-truth harness self-validation" begin
    KA, KB, KP, k, L, Keq, A, B, P = 1.3, 0.8, 0.9, 2.1, 0.7, 3.0, 1.1, 0.5, 0.6
    kr = k * KP / (Keq * KA * KB)
    nonallo = (k*A*B/(KA*KB) - kr*P/KP) / (1 + A/KA + A*B/(KA*KB) + P/KP)

    # (a) L = 0 : the inactive conformation is unpopulated, so the multi-:OnlyA
    #     flux reduces to the single-conformation (non-allosteric) bi-uni rate.
    f0 = multi_onlyA_flux(KA, KB, KP, k, L=0.0, Keq=Keq, A=A, B=B, P=P)
    @test isapprox(f0, nonallo; rtol=1e-4)

    # (b) all-:EqualAI : both conformations are identical, so the flux is the
    #     base bi-uni rate and is independent of L.
    fa = multi_equalAI_flux(KA, KB, KP, k, L, Keq, A, B, P)
    @test isapprox(fa, nonallo; rtol=1e-4)
    @test isapprox(fa, multi_equalAI_flux(KA, KB, KP, k, 5.0, Keq, A, B, P); rtol=1e-4)
end

# ── The gate: :OnlyA MWC derivation matches mass-action ground truth ─────────
# The naive combination Q_A + L·Q_I leaks a bare catalytic rate constant into
# the L-term when the inactive graph fragments; the free-enzyme cross-weighting
# in `_allosteric_num_den_exprs` removes the leak. This is the regression guard.
@testset "OnlyA MWC derivation matches mass-action ground truth" begin
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

# ── The gate: multi-:OnlyA bi-uni derivation matches mass-action ground truth ─
# Both A and B bind :OnlyA, so the inactive graph fragments (D_I = k, D_A = 1)
# and the naive Q_A + L·Q_I combination would leak the bare k_EAB_to_EP constant
# into the L-term. The free-enzyme cross-weighting removes the leak; regression
# guard.
@testset "multi-OnlyA MWC derivation matches mass-action ground truth" begin
    multiA = @allosteric_mechanism begin
        substrates: A, B ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)          :: OnlyA
            E(A) + B ⇌ E(A, B)    :: OnlyA
            E(A, B) <--> E(P)     :: EqualAI
            E + P ⇌ E(P)          :: EqualAI
        end
    end
    fp = ER.fitted_params(multiA)   # (:K_A_A_E, :K_P_E, :K_A_B_EA, :k_EAB_to_EP, :L)
    @test fp == (:K_A_A_E, :K_P_E, :K_A_B_EA, :k_EAB_to_EP, :L)

    rng = MersenneTwister(20260713)
    for _ in 1:5
        KA = 0.5 + 2rand(rng); KB = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng)
        k = 0.5 + 2rand(rng); L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        A = 0.5 + 2rand(rng); B = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map fitted_params -> ground-truth params:
        #   K_A_A_E=KA, K_P_E=KP, K_A_B_EA=KB, k_EAB_to_EP=k.
        prm = NamedTuple{(fp..., :Keq, :E_total)}((KA, KP, KB, k, L, Keq, 1.0))
        v_code = real(ER.rate_equation(multiA, (A=A, B=B, P=P), prm))
        v_gt = multi_onlyA_flux(KA, KB, KP, k, L=L, Keq=Keq, A=A, B=B, P=P)
        @test isapprox(v_code, v_gt; rtol=1e-4)
    end
end
