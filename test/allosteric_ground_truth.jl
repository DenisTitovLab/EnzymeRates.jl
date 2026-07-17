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
# S binds :OnlyA (only the active conformation) and catalysis is :OnlyA (inactive
# rate k_I = 0), P :EqualAI. Forms: E_A, ES_A, EP_A (active) and E_I, EP_I
# (inactive). The inactive conformation has no catalytic edge — that is what
# :OnlyA catalysis means, and with no inactive S binding it never reaches ES.
# Fast RE bindings and flips use FAST; catalysis is O(1). Detailed-balance flip
# ratio [X_I]/[X_A] = L·∏(K_A_i/K_I_i); every present flip here carries an
# :EqualAI ligand (or none), so the ratio is L.
function uni_onlyA_flux(KA, KP, k; L, Keq, S, P, FAST=1e7)
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

# ── Multi-:OnlyA bi-uni network (both substrates bind the active state only) ──
# A + B ⇌ P. Both A and B bind :OnlyA (active conformation only) and catalysis is
# :OnlyA (inactive rate k_I = 0), P :EqualAI. Active forms E_A, EA_A, EAB_A, EP_A.
# Inactive forms E_I, EP_I, RE-connected via P binding — a single segment with no
# catalytic edge (that is what :OnlyA catalysis means; with A,B :OnlyA the inactive
# state never reaches EA or EAB), so the inactive free-enzyme weight is D_I = 1.
# Flip ratio [X_I]/[X_A] = L·∏(K_A_i/K_I_i): free enzyme and EP (:EqualAI or bare)
# flip with ratio L; an :OnlyA-ligand-bearing state has ratio 0 (no flip), so EA
# and EAB never flip.
function multi_onlyA_flux(KA, KB, KP, k; L, Keq, A, B, P, FAST=1e7)
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
# S binds :OnlyA and catalysis is :OnlyA (k_I = 0), so the inactive conformation
# runs no catalysis: its free-enzyme graph is a single RE segment (d_free_I = 1)
# and the derivation takes the raw Q_A + L·Q_I normalization. The gate checks the
# derived rate against the independent mass-action ground truth.
@testset "OnlyA MWC derivation matches mass-action ground truth" begin
    onlyA = @allosteric_mechanism begin
        substrates: S ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + S ⇌ E(S)     :: OnlyA
            E(S) <--> E(P)   :: OnlyA
            E + P ⇌ E(P)     :: EqualAI
        end
    end
    fp = ER.fitted_params(onlyA)          # (:K_P_E, :K_A_S_E, :k_A_ES_to_EP, :L)
    @test fp == (:K_P_E, :K_A_S_E, :k_A_ES_to_EP, :L)

    rng = MersenneTwister(20260713)
    for _ in 1:5
        KA = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng); k = 0.5 + 2rand(rng)
        L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        S = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map fitted_params -> ground-truth params: K_A_S_E=KA, K_P_E=KP, k_A_ES_to_EP=k.
        prm = NamedTuple{(fp..., :Keq, :E_total)}((KP, KA, k, L, Keq, 1.0))
        v_code = real(ER.rate_equation(onlyA, (S=S, P=P), prm))
        v_gt = uni_onlyA_flux(KA, KP, k, L=L, Keq=Keq, S=S, P=P)
        @test isapprox(v_code, v_gt; rtol=1e-4)
        @test isfinite(ER._kcat_forward(onlyA, prm))
    end
end

# ── The gate: multi-:OnlyA bi-uni derivation matches mass-action ground truth ─
# Both A and B bind :OnlyA and catalysis is :OnlyA (k_I = 0), so the inactive
# conformation runs no catalysis: its free-enzyme graph is a single RE segment
# (D_I = 1) and the derivation takes the raw Q_A + L·Q_I normalization. The gate
# checks the derived rate against the independent mass-action ground truth.
@testset "multi-OnlyA MWC derivation matches mass-action ground truth" begin
    multiA = @allosteric_mechanism begin
        substrates: A, B ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)          :: OnlyA
            E(A) + B ⇌ E(A, B)    :: OnlyA
            E(A, B) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)          :: EqualAI
        end
    end
    fp = ER.fitted_params(multiA)   # (:K_A_A_E, :K_P_E, :K_A_B_EA, :k_A_EAB_to_EP, :L)
    @test fp == (:K_A_A_E, :K_P_E, :K_A_B_EA, :k_A_EAB_to_EP, :L)

    rng = MersenneTwister(20260713)
    for _ in 1:5
        KA = 0.5 + 2rand(rng); KB = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng)
        k = 0.5 + 2rand(rng); L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        A = 0.5 + 2rand(rng); B = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map fitted_params -> ground-truth params:
        #   K_A_A_E=KA, K_P_E=KP, K_A_B_EA=KB, k_A_EAB_to_EP=k.
        prm = NamedTuple{(fp..., :Keq, :E_total)}((KA, KP, KB, k, L, Keq, 1.0))
        v_code = real(ER.rate_equation(multiA, (A=A, B=B, P=P), prm))
        v_gt = multi_onlyA_flux(KA, KB, KP, k, L=L, Keq=Keq, A=A, B=B, P=P)
        @test isapprox(v_code, v_gt; rtol=1e-4)
        @test isfinite(ER._kcat_forward(multiA, prm))
    end
end

# ── Metabolite-bearing-D ordered bi-uni :OnlyA network (the LDH i-state case) ──
# S + B ⇌ P. S binds :OnlyA via a STEADY-STATE step (E + S <--> E(S)); B binds
# :EqualAI at rapid equilibrium on E(S); catalysis E(S,B) <--> E(P) is
# steady-state :OnlyA (inactive rate k_I = 0); P binds :EqualAI at rapid
# equilibrium. Because B binds rapidly on the steady-state catalytic path, the
# active free-enzyme spanning-tree weight D[g_free] CARRIES the metabolite B:
# D_A = koff_S + k·B/K_B. The inactive conformation runs no catalysis and never
# binds S or B (both :OnlyA), so its graph is a single rapid-equilibrium segment
# {E_I, E(P)_I} and D_I = 1. D_A ≠ D_I — a metabolite-bearing active weight against
# a bare inactive one — the cross-weight regime the fix re-baselines and which has
# no other ground truth. Flip ratio [X_I]/[X_A] = L·∏(K_A/K_I): free enzyme and
# E(P) flip (ratio L); an S-bearing state has ratio 0 (S :OnlyA), so E(S) and
# E(S,B) never flip. Reverse catalysis kr from the Haldane relation.
function metab_dfree_onlyA_flux(kon, koff, KB, KP, k; L, Keq, S, B, P, FAST=1e7)
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
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# ── Metabolite-bearing-D all-:EqualAI network (both conformations identical) ───
# Same topology, but S also binds the inactive conformation (SS, K_I = K_A) and
# every form flips (ratio L). With K_I = K_A everywhere the two conformations are
# indistinguishable, so the flux is the non-allosteric ordered bi-uni rate,
# independent of L.
function metab_dfree_equalAI_flux(kon, koff, KB, KP, k, L, Keq, S, B, P; FAST=1e7)
    kr = k * kon * KP / (koff * KB * Keq)
    species = [:E_A, :ES_A, :ESB_A, :EP_A, :E_I, :ES_I, :ESB_I, :EP_I]
    edges = [
        (:E_A, :ES_A, kon * S), (:ES_A, :E_A, koff),
        (:E_I, :ES_I, kon * S), (:ES_I, :E_I, koff),          # inactive S binding (K_I = K_A)
        (:ES_A, :ESB_A, FAST * B / KB), (:ESB_A, :ES_A, FAST),
        (:ES_I, :ESB_I, FAST * B / KB), (:ESB_I, :ES_I, FAST),
        (:E_A, :EP_A, FAST * P / KP), (:EP_A, :E_A, FAST),
        (:E_I, :EP_I, FAST * P / KP), (:EP_I, :E_I, FAST),
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),
        (:ES_A, :ES_I, FAST * L), (:ES_I, :ES_A, FAST),       # ES flip, ratio L
        (:ESB_A, :ESB_I, FAST * L), (:ESB_I, :ESB_A, FAST),   # ESB flip, ratio L
        (:EP_A, :EP_I, FAST * L), (:EP_I, :EP_A, FAST),
        (:ESB_A, :EP_A, k), (:EP_A, :ESB_A, kr),
        (:ESB_I, :EP_I, k), (:EP_I, :ESB_I, kr),
    ]
    cat_edges = [(:ESB_A, :EP_A, k, kr), (:ESB_I, :EP_I, k, kr)]
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# ── Single-conformation (non-allosteric) reference for the ordered bi-uni ─────
# The active mechanism alone (E, E(S), E(S,B), E(P)); no inactive conformation,
# no flips. The L → 0 limit and the all-:EqualAI flux must both equal this rate.
function metab_dfree_base_flux(kon, koff, KB, KP, k, Keq, S, B, P; FAST=1e7)
    kr = k * kon * KP / (koff * KB * Keq)
    species = [:E, :ES, :ESB, :EP]
    edges = [
        (:E, :ES, kon * S), (:ES, :E, koff),
        (:ES, :ESB, FAST * B / KB), (:ESB, :ES, FAST),
        (:E, :EP, FAST * P / KP), (:EP, :E, FAST),
        (:ESB, :EP, k), (:EP, :ESB, kr),
    ]
    mwc_ground_truth_flux(species, edges, [(:ESB, :EP, k, kr)], 1.0)
end

# ── Metabolite-bearing-D harness self-validation ─────────────────────────────
@testset "metabolite-bearing-D ground-truth harness self-validation" begin
    kon, koff, KB, KP, k, Keq, S, B, P = 1.7, 1.1, 0.8, 0.9, 2.1, 3.0, 1.1, 0.5, 0.6
    base = metab_dfree_base_flux(kon, koff, KB, KP, k, Keq, S, B, P)

    # (a) L = 0 : the inactive conformation is unpopulated, so the :OnlyA flux
    #     reduces to the single-conformation (non-allosteric) ordered bi-uni rate.
    f0 = metab_dfree_onlyA_flux(kon, koff, KB, KP, k, L=0.0, Keq=Keq, S=S, B=B, P=P)
    @test isapprox(f0, base; rtol=1e-4)

    # (b) all-:EqualAI : both conformations are identical, so the flux is the base
    #     rate and is independent of L.
    fa = metab_dfree_equalAI_flux(kon, koff, KB, KP, k, 0.7, Keq, S, B, P)
    @test isapprox(fa, base; rtol=1e-4)
    @test isapprox(fa, metab_dfree_equalAI_flux(kon, koff, KB, KP, k, 5.0, Keq, S, B, P); rtol=1e-4)
end

# ── The gate: metabolite-bearing-D :OnlyA derivation matches mass-action GT ────
# The LDH i-state re-baseline. S binds :OnlyA on the steady-state catalytic path
# and B binds rapidly there too, so the free-enzyme spanning-tree weight D[g_free]
# carries the metabolite B and D_A ≠ D_I. The cross-weighting `den = D_I·Q_A +
# L·D_A·Q_I` supplies the common free-enzyme basis this regime needs; the naive
# Q_A + L·Q_I combination gets it wrong. Regression guard.
@testset "metabolite-bearing-D MWC derivation matches mass-action ground truth" begin
    metabD = @allosteric_mechanism begin
        substrates: S, B ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + S <--> E(S)      :: OnlyA
            E(S) + B ⇌ E(S, B)   :: EqualAI
            E(S, B) <--> E(P)    :: OnlyA
            E + P ⇌ E(P)         :: EqualAI
        end
    end
    fp = ER.fitted_params(metabD)  # (:K_P_E,:kon_A_S_E,:koff_A_S_E,:k_A_EBS_to_EP,:K_B_ES,:L)
    @test fp == (:K_P_E, :kon_A_S_E, :koff_A_S_E, :k_A_EBS_to_EP, :K_B_ES, :L)

    rng = MersenneTwister(20260713)
    for _ in 1:5
        kon = 0.5 + 2rand(rng); koff = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng)
        KB = 0.5 + 2rand(rng); k = 0.5 + 2rand(rng); L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        S = 0.5 + 2rand(rng); B = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map fitted_params -> ground-truth params:
        #   kon_A_S_E=kon, koff_A_S_E=koff, K_B_ES=KB, K_P_E=KP, k_A_EBS_to_EP=k.
        prm = NamedTuple{(fp..., :Keq, :E_total)}((KP, kon, koff, k, KB, L, Keq, 1.0))
        v_code = real(ER.rate_equation(metabD, (S=S, B=B, P=P), prm))
        v_gt = metab_dfree_onlyA_flux(kon, koff, KB, KP, k, L=L, Keq=Keq, S=S, B=B, P=P)
        @test isapprox(v_code, v_gt; rtol=1e-4)
        @test isfinite(ER._kcat_forward(metabD, prm))
    end
end

# ── Ordered bi-uni :NonequalAI-catalysis network (the per-form-flip reference) ─
# A + B ⇌ P. A binds :EqualAI via a STEADY-STATE step (E + A <--> E(A), kon·A /
# koff); B binds :EqualAI at rapid equilibrium on E(A); catalysis E(A,B) <--> E(P)
# is STEADY-STATE :NonequalAI (rate k_A in the active conformation, k_I in the
# inactive); P binds :EqualAI at rapid equilibrium. Every ligand is :EqualAI
# (K_I = K_A), so every form — including the catalytic intermediates — flips with
# ratio [X_I]/[X_A] = L.
#
# This is the PER-FORM-FLIP model (the enzyme may change conformation mid-cycle),
# retained here only to self-validate the harness. It is NOT the model the package
# derives: `biuni_nonequalAI_freeflip_flux` (only the free enzyme flips) is the
# formulation-1 reference the derivation gate checks. The two differ by ~0.1–3%
# whenever k_A ≠ k_I, because per-form flipping routes turnover through the faster
# conformation. Reverse catalysis kr from the Haldane relation.
function biuni_nonequalAI_flux(kon, koff, KB, KP; k_A, k_I, L, Keq, A, B, P, FAST=1e7)
    krA = k_A * kon * KP / (koff * KB * Keq)
    krI = k_I * kon * KP / (koff * KB * Keq)
    species = [:E_A, :EA_A, :EAB_A, :EP_A, :E_I, :EA_I, :EAB_I, :EP_I]
    edges = [
        (:E_A, :EA_A, kon * A), (:EA_A, :E_A, koff),          # active A binding (SS)
        (:E_I, :EA_I, kon * A), (:EA_I, :E_I, koff),          # inactive A binding (SS, K_I = K_A)
        (:EA_A, :EAB_A, FAST * B / KB), (:EAB_A, :EA_A, FAST),# active B binding (RE)
        (:EA_I, :EAB_I, FAST * B / KB), (:EAB_I, :EA_I, FAST),# inactive B binding (RE)
        (:E_A, :EP_A, FAST * P / KP), (:EP_A, :E_A, FAST),    # active P binding (RE)
        (:E_I, :EP_I, FAST * P / KP), (:EP_I, :E_I, FAST),    # inactive P binding (RE)
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),           # free-enzyme flip, ratio L
        (:EA_A, :EA_I, FAST * L), (:EA_I, :EA_A, FAST),       # EA flip, ratio L
        (:EAB_A, :EAB_I, FAST * L), (:EAB_I, :EAB_A, FAST),   # EAB flip, ratio L
        (:EP_A, :EP_I, FAST * L), (:EP_I, :EP_A, FAST),       # EP flip, ratio L
        (:EAB_A, :EP_A, k_A), (:EP_A, :EAB_A, krA),         # active catalysis (SS, k_A)
        (:EAB_I, :EP_I, k_I), (:EP_I, :EAB_I, krI),         # inactive catalysis (SS, k_I)
    ]
    cat_edges = [(:EAB_A, :EP_A, k_A, krA), (:EAB_I, :EP_I, k_I, krI)]
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# ── :NonequalAI-catalysis harness self-validation ────────────────────────────
# First ground the hand-written mass-action reference `metab_dfree_base_flux`
# (the single-conformation ordered bi-uni, reused as the non-allosteric rate at
# k = k_A) against the ODE-validated non-allosteric `rate_equation` — an
# independent, separately ODE-cross-checked code path. This closes the trust loop
# for the whole `mwc_ground_truth_flux` harness. Then confirm the two-conformation
# `biuni_nonequalAI_flux` degenerates correctly: (a) L = 0 (inactive unpopulated)
# → the base rate at k = k_A; (b) k_I = k_A (conformations identical) → the base
# rate, independent of L.
@testset ":NonequalAI-catalysis ground-truth harness self-validation" begin
    base = @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P)
            E + P ⇌ E(P)
        end
    end
    bfp = ER.fitted_params(base)

    rng = MersenneTwister(20260713)
    for _ in 1:5
        kon = 0.5 + 2rand(rng); koff = 0.5 + 2rand(rng)
        KB = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng)
        kA = 0.5 + 2rand(rng); kI = 0.5 + 2rand(rng); Keq = 2.0 + 2rand(rng)
        A = 0.5 + 2rand(rng); B = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        L = 0.5 + rand(rng)

        base_rate = metab_dfree_base_flux(kon, koff, KB, KP, kA, Keq, A, B, P)

        # `metab_dfree_base_flux` vs the ODE-validated non-allosteric rate_equation.
        bd = Dict(:kon_A_E=>kon, :koff_A_E=>koff, :K_B_EA=>KB, :K_P_E=>KP,
                  :k_EAB_to_EP=>kA)
        bprm = NamedTuple{(bfp..., :Keq, :E_total)}(((bd[s] for s in bfp)..., Keq, 1.0))
        @test isapprox(base_rate,
            real(ER.rate_equation(base, (A=A, B=B, P=P), bprm)); rtol=1e-4)

        # (a) L = 0 : inactive conformation unpopulated → base rate at k_A.
        f0 = biuni_nonequalAI_flux(kon, koff, KB, KP;
            k_A=kA, k_I=kI, L=0.0, Keq=Keq, A=A, B=B, P=P)
        @test isapprox(f0, base_rate; rtol=1e-4)

        # (b) k_I = k_A : conformations identical → base rate, independent of L.
        fe = biuni_nonequalAI_flux(kon, koff, KB, KP;
            k_A=kA, k_I=kA, L=L, Keq=Keq, A=A, B=B, P=P)
        fe5 = biuni_nonequalAI_flux(kon, koff, KB, KP;
            k_A=kA, k_I=kA, L=5.0, Keq=Keq, A=A, B=B, P=P)
        @test isapprox(fe, base_rate; rtol=1e-4)
        @test isapprox(fe, fe5; rtol=1e-4)
    end
end

# ── Free-flip-only ordered bi-uni :NonequalAI-catalysis reference ────────────
# Formulation-1 reference: only the free enzyme flips conformation. Each
# conformation runs its own catalytic cycle, coupled only through the shared
# free-enzyme pool. Same forms/edges as `biuni_nonequalAI_flux` but with only
# the E_A<->E_I flip — the model this package derives (commit-when-free).
function biuni_nonequalAI_freeflip_flux(kon, koff, KB, KP; k_A, k_I, L, Keq, A, B, P, FAST=1e7)
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
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# ── Free-flip-only :NonequalAI-catalysis harness self-validation ────────────
# Same physics checks as the topology-sharing harness above, applied to the
# free-flip-only reference: (a) L = 0 (inactive unpopulated) → the base rate at
# k = k_A; (b) k_I = k_A (conformations identical) → the base rate, independent
# of L.
@testset ":NonequalAI-catalysis free-flip-only ground-truth harness self-validation" begin
    rng = MersenneTwister(20260713)
    for _ in 1:5
        kon = 0.5 + 2rand(rng); koff = 0.5 + 2rand(rng)
        KB = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng)
        kA = 0.5 + 2rand(rng); kI = 0.5 + 2rand(rng); Keq = 2.0 + 2rand(rng)
        A = 0.5 + 2rand(rng); B = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        L = 0.5 + rand(rng)

        base_rate = metab_dfree_base_flux(kon, koff, KB, KP, kA, Keq, A, B, P)

        # (a) L = 0 : inactive conformation unpopulated → base rate at k_A.
        f0 = biuni_nonequalAI_freeflip_flux(kon, koff, KB, KP;
            k_A=kA, k_I=kI, L=0.0, Keq=Keq, A=A, B=B, P=P)
        @test isapprox(f0, base_rate; rtol=1e-4)

        # (b) k_I = k_A : conformations identical → base rate, independent of L.
        fe = biuni_nonequalAI_freeflip_flux(kon, koff, KB, KP;
            k_A=kA, k_I=kA, L=L, Keq=Keq, A=A, B=B, P=P)
        fe5 = biuni_nonequalAI_freeflip_flux(kon, koff, KB, KP;
            k_A=kA, k_I=kA, L=5.0, Keq=Keq, A=A, B=B, P=P)
        @test isapprox(fe, base_rate; rtol=1e-4)
        @test isapprox(fe, fe5; rtol=1e-4)
    end

    # (c) k_I ≠ k_A, L > 0 : the free-flip (formulation-1) reference must DIFFER
    #     from the per-form-flip one. That ~0.1–3% gap is the whole point of the
    #     fix, and it guards the :NonequalAI derivation gate against being pointed
    #     back at the per-form-flip reference (which the raw Q_A + L·Q_I combine
    #     matches). Fixed, clearly-unequal rate constants keep the gap unambiguous.
    v_free = biuni_nonequalAI_freeflip_flux(1.7, 1.1, 0.8, 0.9;
        k_A=2.5, k_I=0.4, L=0.7, Keq=3.0, A=1.1, B=0.5, P=0.6)
    v_perform = biuni_nonequalAI_flux(1.7, 1.1, 0.8, 0.9;
        k_A=2.5, k_I=0.4, L=0.7, Keq=3.0, A=1.1, B=0.5, P=0.6)
    @test !isapprox(v_free, v_perform; rtol=1e-4)
end

# ── The gate: :NonequalAI-catalysis derivation matches mass-action GT ─────────
# The productive-inactive case. Catalysis is :NonequalAI (k_A ≠ k_I) with shared
# graph topology, so D_A and D_I differ only by the catalytic rate constant. Under
# formulation 1 (the enzyme commits to one conformation per catalytic cycle) the
# two conformations normalize per-state, so the derivation cross-weights them and
# this gate checks against the free-flip-only reference `biuni_nonequalAI_freeflip_flux`.
# The gate goes red if the derivation reverts to the raw Q_A + L·Q_I combine (which
# matches a per-form-flip model, not formulation 1) or mis-renders the normalization.
@testset ":NonequalAI-catalysis MWC derivation matches mass-action ground truth" begin
    allo = @allosteric_mechanism begin
        substrates: A, B ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A)        :: EqualAI
            E(A) + B ⇌ E(A, B)     :: EqualAI
            E(A, B) <--> E(P)      :: NonequalAI
            E + P ⇌ E(P)           :: EqualAI
        end
    end
    fp = ER.fitted_params(allo)
    @test fp == (:kon_A_E, :koff_A_E, :K_P_E, :K_B_EA,
                 :k_A_EAB_to_EP, :k_I_EAB_to_EP, :L)

    rng = MersenneTwister(20260713)
    for _ in 1:6
        kon = 0.5 + 2rand(rng); koff = 0.5 + 2rand(rng)
        KP = 0.5 + 2rand(rng); KB = 0.5 + 2rand(rng)
        kA = 0.5 + 2rand(rng); kI = 0.5 + 2rand(rng)
        L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        A = 0.5 + 2rand(rng); B = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map fitted_params -> ground-truth params:
        #   kon_A_E=kon, koff_A_E=koff, K_B_EA=KB, K_P_E=KP,
        #   k_A_EAB_to_EP=k_A, k_I_EAB_to_EP=k_I.
        d = Dict(:kon_A_E=>kon, :koff_A_E=>koff, :K_P_E=>KP, :K_B_EA=>KB,
                 :k_A_EAB_to_EP=>kA, :k_I_EAB_to_EP=>kI, :L=>L)
        prm = NamedTuple{(fp..., :Keq, :E_total)}(((d[s] for s in fp)..., Keq, 1.0))
        v_code = real(ER.rate_equation(allo, (A=A, B=B, P=P), prm))
        v_gt = biuni_nonequalAI_freeflip_flux(kon, koff, KB, KP;
            k_A=kA, k_I=kI, L=L, Keq=Keq, A=A, B=B, P=P)
        @test isapprox(v_code, v_gt; rtol=1e-4)
    end
end

# ── The gate: a metabolite inside D[g_free] at zero concentration ─────────────
# When D[g_free] carries a metabolite, the derivation must handle that metabolite
# going to zero correctly. Both mechanisms here keep a finite free-enzyme weight
# at B = 0: the active steady-state S/A-binding leaves a bare koff term in
# D[g_free], so a reverse path from product back to free enzyme stays open and the
# reverse flux is nonzero. The :OnlyA metabolite-D mechanism and the productive
# :NonequalAI mechanism both retain surviving reverse flux at B = 0. Both are
# checked against the ground truth.
@testset "metabolite in D[g_free] at zero concentration" begin
    metabD = @allosteric_mechanism begin
        substrates: S, B ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + S <--> E(S) :: OnlyA ; E(S) + B ⇌ E(S, B) :: EqualAI
            E(S, B) <--> E(P) :: OnlyA ; E + P ⇌ E(P) :: EqualAI
        end
    end
    fp = ER.fitted_params(metabD)   # (:K_P_E,:kon_A_S_E,:koff_A_S_E,:k_A_EBS_to_EP,:K_B_ES,:L)
    kon, koff, KP, KB, k, L, Keq, S, P = 1.7, 1.1, 0.9, 0.8, 2.1, 0.7, 3.0, 1.1, 0.6
    prm = NamedTuple{(fp..., :Keq, :E_total)}((KP, kon, koff, k, KB, L, Keq, 1.0))
    v0 = real(ER.rate_equation(metabD, (S=S, B=0.0, P=P), prm))
    @test isapprox(v0,
        metab_dfree_onlyA_flux(kon, koff, KB, KP, k; L=L, Keq=Keq, S=S, B=0.0, P=P); rtol=1e-4)
    @test v0 < 0    # reverse flux survives at B=0 (invalid :EqualAI-cat trap→0 was fake)
    @test isfinite(ER._kcat_forward(metabD, prm))

    allo = @allosteric_mechanism begin
        substrates: A, B ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A) :: EqualAI ; E(A) + B ⇌ E(A, B) :: EqualAI
            E(A, B) <--> E(P) :: NonequalAI ; E + P ⇌ E(P) :: EqualAI
        end
    end
    fpn = ER.fitted_params(allo)
    kon2, koff2, KP2, KB2 = 1.7, 1.1, 0.9, 0.8
    kA, kI, L2, Keq2, A2, P2 = 2.5, 0.4, 0.7, 3.0, 1.1, 0.9
    d = Dict(:kon_A_E=>kon2, :koff_A_E=>koff2, :K_P_E=>KP2, :K_B_EA=>KB2,
             :k_A_EAB_to_EP=>kA, :k_I_EAB_to_EP=>kI, :L=>L2)
    prm2 = NamedTuple{(fpn..., :Keq, :E_total)}(((d[s] for s in fpn)..., Keq2, 1.0))
    vN = real(ER.rate_equation(allo, (A=A2, B=0.0, P=P2), prm2))
    @test isapprox(vN, biuni_nonequalAI_freeflip_flux(kon2, koff2, KB2, KP2;
        k_A=kA, k_I=kI, L=L2, Keq=Keq2, A=A2, B=0.0, P=P2); rtol=1e-4)
    @test abs(vN) > 1e-3                                     # reverse flux survives at B = 0
end

# ── The gate: allosteric ping-pong (single free form, formulation-1 normalization) ─
# A ping-pong mechanism has ONE free enzyme form — the covalent intermediate
# `E(; residual = A - P)` carries a residual, so it is not free — hence D[g_free]
# is well-defined and the derivation normalizes it like any other mechanism (no
# special case, no guard). This gate confirms the derivation handles the
# four-segment residual-bearing King–Altman and that the normalization collapses
# correctly for identical conformations: an all-:EqualAI ping-pong must equal the
# non-allosteric ping-pong and be independent of L.
@testset "allosteric ping-pong self-consistency" begin
    nonallo = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B <--> E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end
    alloEq = @allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A)                                       :: EqualAI
            E(A) <--> E(; residual = A - P) + P                   :: EqualAI
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
            E(B; residual = A - P) <--> E + Q                     :: EqualAI
        end
    end
    fpn = ER.fitted_params(nonallo); fpa = ER.fitted_params(alloEq)
    rng = MersenneTwister(20260714)
    for _ in 1:5
        vals = Dict(s => 0.5 + 2rand(rng) for s in fpn)
        Keq = 2.0 + 2rand(rng)
        concs = (A=0.5+2rand(rng), B=0.5+2rand(rng), P=0.5+2rand(rng), Q=0.5+2rand(rng))
        pn = NamedTuple{(fpn..., :Keq, :E_total)}(((vals[s] for s in fpn)..., Keq, 1.0))
        vn = real(ER.rate_equation(nonallo, concs, pn))
        for L in (0.7, 5.0)          # all-:EqualAI is L-independent and equals non-allosteric
            pa = NamedTuple{(fpa..., :Keq, :E_total)}(
                ((s === :L ? L : vals[s] for s in fpa)..., Keq, 1.0))
            @test isapprox(real(ER.rate_equation(alloEq, concs, pa)), vn; rtol=1e-6)
        end
    end
end

# ── Two-conformation ping-pong bi-bi network (formulation 1) ─────────────────
# Ping-pong has TWO empty-bound forms — free E and the covalent intermediate F
# (`E(; residual = A - P)`). Only free E flips; F does not. The E_A<->E_I edge is
# therefore the ONLY cut between the two conformation subnetworks, so at steady
# state it carries zero net flux and E_I/E_A sits at exactly L for ANY FAST — the
# fast-flip limit is exact for this topology rather than a large-FAST limit. Each
# conformation turns its own four-form cycle (E, EA, F, FB), coupled only through
# the shared free-enzyme pool.
#
# Thermodynamics pins ONE combination per conformation, not one per half-reaction.
# Detailed balance around the closed cycle E → EA → F → FB → E requires
#   (kon_A·A/koff_A)·(k_P/(koff_P·P))·(kon_B·B/koff_B)·(k_Q/(koff_Q·Q)) = 1
# at the equilibrium ratio P·Q/(A·B) = Keq, i.e. the overall Haldane relation
#   kon_A·k_P·kon_B·k_Q = Keq · koff_A·koff_P·koff_B·koff_Q,
# matching the numerator of Segel Eq. IX-140 (k1f·k2f·k3f·k4f·A·B −
# k1r·k2r·k3r·k4r·P·Q). The two half-reactions' equilibrium constants are NOT
# separately fixed — only their product is — so exactly one reverse constant per
# conformation is dependent. `koff_A` (shared, :EqualAI) closes the active cycle;
# `koff_P_I` then closes the inactive one against that same `koff_A`.
function pingpong_nonequalAI_freeflip_flux(kon_A, kon_B, koff_B;
        k_P_A, koff_P_A, k_Q_A, koff_Q_A, k_P_I, k_Q_I, koff_Q_I,
        L, Keq, A, B, P, Q, FAST=1e7)
    koff_A   = kon_A * k_P_A * kon_B * k_Q_A / (Keq * koff_P_A * koff_B * koff_Q_A)
    koff_P_I = kon_A * k_P_I * kon_B * k_Q_I / (Keq * koff_A * koff_B * koff_Q_I)
    species = [:E_A, :EA_A, :F_A, :FB_A, :E_I, :EA_I, :F_I, :FB_I]
    edges = Tuple{Symbol,Symbol,Float64}[]
    cat_edges = Tuple{Symbol,Symbol,Float64,Float64}[]
    for (c, k_P, koff_P, k_Q, koff_Q) in ((:A, k_P_A, koff_P_A, k_Q_A, koff_Q_A),
                                          (:I, k_P_I, koff_P_I, k_Q_I, koff_Q_I))
        e, ea = Symbol(:E_, c), Symbol(:EA_, c)
        f, fb = Symbol(:F_, c), Symbol(:FB_, c)
        append!(edges, [
            (e, ea, kon_A * A), (ea, e, koff_A),          # E + A ⇌ EA
            (ea, f, k_P), (f, ea, koff_P * P),            # EA ⇌ F + P
            (f, fb, kon_B * B), (fb, f, koff_B),          # F + B ⇌ FB
            (fb, e, k_Q), (e, fb, koff_Q * Q),            # FB ⇌ E + Q
        ])
        push!(cat_edges, (ea, f, k_P, koff_P * P))        # net flux across the P cut
    end
    push!(edges, (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST))   # only free enzyme flips
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# `test/mechanism_definitions_for_test_enzyme_derivation.jl:328` transcribes this
# same Segel formula, and both transcriptions are live. Keep them independent
# rather than sharing one: a shared transcription error would green this gate and
# that one at once, whereas two independent transcriptions cross-check each other.
# Sharing would also couple this gate to the MECHANISM_TEST_SPECS fixture, whose
# copy takes `(params::NamedTuple, concs::NamedTuple)` with an `Etotal` rather
# than the 12 positional scalars this one takes.
"Segel Eq. IX-140 ping-pong bi-bi rate: E + A ⇌ EA ⇌ F + P; F + B ⇌ FB ⇌ E + Q."
function segel_pingpong_flux(k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, A, B, P, Q)
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
        kon_A = 0.5+2rand(rng); kon_B = 0.5+2rand(rng); koff_B = 0.5+2rand(rng)
        k_P_A = 0.5+2rand(rng); koff_P_A = 0.5+2rand(rng)
        k_Q_A = 0.5+2rand(rng); koff_Q_A = 0.5+2rand(rng)
        k_P_I = 0.5+2rand(rng); k_Q_I = 0.5+2rand(rng); koff_Q_I = 0.5+2rand(rng)
        L = 0.5+rand(rng); Keq = 2.0+2rand(rng)
        A = 0.5+2rand(rng); B = 0.5+2rand(rng); P = 0.5+2rand(rng); Q = 0.5+2rand(rng)
        act = (k_P_A=k_P_A, koff_P_A=koff_P_A, k_Q_A=k_Q_A, koff_Q_A=koff_Q_A)
        ina = (k_P_I=k_P_I, k_Q_I=k_Q_I, koff_Q_I=koff_Q_I)

        # (a) L = 0 : inactive unpopulated → the active-only ping-pong, checked
        #     against the Segel closed form. `koff_A` is the dependent reverse
        #     constant the Haldane fixes; Segel's k1r is that same constant.
        koff_A = kon_A*k_P_A*kon_B*k_Q_A / (Keq*koff_P_A*koff_B*koff_Q_A)
        f0 = pingpong_nonequalAI_freeflip_flux(kon_A, kon_B, koff_B; act..., ina...,
            L=0.0, Keq=Keq, A=A, B=B, P=P, Q=Q)
        @test isapprox(f0, segel_pingpong_flux(kon_A, koff_A, k_P_A, koff_P_A,
            kon_B, koff_B, k_Q_A, koff_Q_A, A, B, P, Q); rtol=1e-9)

        # (b) identical conformations → the active-only rate, independent of L.
        same = (k_P_I=k_P_A, k_Q_I=k_Q_A, koff_Q_I=koff_Q_A)
        fe = pingpong_nonequalAI_freeflip_flux(kon_A, kon_B, koff_B; act..., same...,
            L=L, Keq=Keq, A=A, B=B, P=P, Q=Q)
        fe5 = pingpong_nonequalAI_freeflip_flux(kon_A, kon_B, koff_B; act..., same...,
            L=5.0, Keq=Keq, A=A, B=B, P=P, Q=Q)
        @test isapprox(fe, f0; rtol=1e-9)
        @test isapprox(fe, fe5; rtol=1e-9)

        # (c) v = 0 at the equilibrium metabolite ratio P·Q/(A·B) = Keq. Both
        #     conformations are live and unequal, so this pins both Haldanes.
        Qeq = Keq * A * B / P
        @test abs(pingpong_nonequalAI_freeflip_flux(kon_A, kon_B, koff_B; act...,
            ina..., L=L, Keq=Keq, A=A, B=B, P=P, Q=Qeq)) < 1e-9

        # (d) FAST-invariance: the free-enzyme flip is the only cut between the
        #     conformation subnetworks, so it carries zero net flux and E_I/E_A is
        #     exactly L for any FAST. The fast-flip limit is therefore exact here.
        v_live = pingpong_nonequalAI_freeflip_flux(kon_A, kon_B, koff_B; act..., ina...,
            L=L, Keq=Keq, A=A, B=B, P=P, Q=Q)
        for fast in (1e2, 1e12)
            @test isapprox(v_live, pingpong_nonequalAI_freeflip_flux(kon_A, kon_B,
                koff_B; act..., ina..., L=L, Keq=Keq, A=A, B=B, P=P, Q=Q, FAST=fast);
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
# goes red if the derivation flips F, reverts to the raw Q_A + L·Q_I combine, or
# mis-renders the normalization. Only the raw-combine mode has a measured margin:
# against a per-form-flip (formulation-2) oracle it deviates 0.95%-93% from the
# derivation, orders of magnitude above the 1e-10 tolerance below.
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
    allo = @allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A)                                       :: EqualAI
            E(A) <--> E(; residual = A - P) + P                   :: NonequalAI
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
            E(B; residual = A - P) <--> E + Q                     :: NonequalAI
        end
    end
    fp = ER.fitted_params(allo)
    @test fp == (:kon_A_E, :kon_A_P_EA, :koff_A_P_EA,
                 Symbol("kon_A_Q_EB_res_+A_-P"), Symbol("koff_A_Q_EB_res_+A_-P"),
                 Symbol("kon_B_E_res_+A_-P"), Symbol("koff_B_E_res_+A_-P"),
                 :kon_I_P_EA, Symbol("kon_I_Q_EB_res_+A_-P"),
                 Symbol("koff_I_Q_EB_res_+A_-P"), :L)

    rng = MersenneTwister(20260716)
    for _ in 1:6
        kon_A = 0.5+2rand(rng); kon_B = 0.5+2rand(rng); koff_B = 0.5+2rand(rng)
        k_P_A = 0.5+2rand(rng); koff_P_A = 0.5+2rand(rng)
        k_Q_A = 0.5+2rand(rng); koff_Q_A = 0.5+2rand(rng)
        k_P_I = 0.5+2rand(rng); k_Q_I = 0.5+2rand(rng); koff_Q_I = 0.5+2rand(rng)
        L = 0.5+rand(rng); Keq = 2.0+2rand(rng)
        A = 0.5+2rand(rng); B = 0.5+2rand(rng); P = 0.5+2rand(rng); Q = 0.5+2rand(rng)
        # Map fitted_params -> ground-truth params. On a release step the
        # canonical direction runs E(A) → F + P, so `kon_…` is the forward
        # (product-releasing) rate and `koff_…` the reverse (product-rebinding)
        # one — the binding steps read the usual way round.
        #   kon_A_E=kon_A                                (E + A ⇌ EA, shared)
        #   kon_B_E_res_+A_-P=kon_B, koff_B_E_res_+A_-P=koff_B  (F + B ⇌ FB, shared)
        #   kon_A_P_EA=k_P_A, koff_A_P_EA=koff_P_A       (EA ⇌ F + P, active)
        #   kon_I_P_EA=k_P_I                             (EA ⇌ F + P, inactive)
        #   kon_A_Q_EB_res_+A_-P=k_Q_A,
        #   koff_A_Q_EB_res_+A_-P=koff_Q_A               (FB ⇌ E + Q, active)
        #   kon_I_Q_EB_res_+A_-P=k_Q_I,
        #   koff_I_Q_EB_res_+A_-P=koff_Q_I               (FB ⇌ E + Q, inactive)
        # `koff_A_E` and `koff_I_P_EA` are absent from fitted_params: each
        # conformation's Haldane makes one reverse constant dependent, and the
        # oracle derives exactly those two.
        d = Dict(:kon_A_E => kon_A,
                 Symbol("kon_B_E_res_+A_-P") => kon_B,
                 Symbol("koff_B_E_res_+A_-P") => koff_B,
                 :kon_A_P_EA => k_P_A, :koff_A_P_EA => koff_P_A,
                 Symbol("kon_A_Q_EB_res_+A_-P") => k_Q_A,
                 Symbol("koff_A_Q_EB_res_+A_-P") => koff_Q_A,
                 :kon_I_P_EA => k_P_I,
                 Symbol("kon_I_Q_EB_res_+A_-P") => k_Q_I,
                 Symbol("koff_I_Q_EB_res_+A_-P") => koff_Q_I,
                 :L => L)
        prm = NamedTuple{(fp..., :Keq, :E_total)}(((d[s] for s in fp)..., Keq, 1.0))
        v_code = real(ER.rate_equation(allo, (A=A, B=B, P=P, Q=Q), prm))
        v_gt = pingpong_nonequalAI_freeflip_flux(kon_A, kon_B, koff_B;
            k_P_A=k_P_A, koff_P_A=koff_P_A, k_Q_A=k_Q_A, koff_Q_A=koff_Q_A,
            k_P_I=k_P_I, k_Q_I=k_Q_I, koff_Q_I=koff_Q_I,
            L=L, Keq=Keq, A=A, B=B, P=P, Q=Q)
        @test isapprox(v_code, v_gt; rtol=1e-10)
    end
end

# ── Fully-inert inactive network (both substrate and product bind :OnlyA) ─────
# S and P both bind :OnlyA, so the inactive conformation binds nothing — it is a
# free-enzyme reservoir of mass L, coupled to the active cycle only by the E flip.
# Its partition is 1, so the denominator must carry L·1 (not 0). Fast RE bindings
# and the flip use FAST; catalysis is O(1). Reverse catalysis from the Haldane
# relation. At L = 0 the reservoir is unpopulated → the non-allosteric rate.
function inert_inactive_flux(; KS, KP, k, L, Keq, S, P, FAST=1e7)
    kr = k * KP / (Keq * KS)
    species = [:E_A, :ES_A, :EP_A, :E_I]
    edges = [
        (:E_A, :ES_A, FAST * S / KS), (:ES_A, :E_A, FAST),   # active S binding (RE)
        (:E_A, :EP_A, FAST * P / KP), (:EP_A, :E_A, FAST),   # active P binding (RE)
        (:ES_A, :EP_A, k), (:EP_A, :ES_A, kr),              # active catalysis (SS)
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),          # free-enzyme flip, ratio L
    ]                                                         # E_I inert: no other edges
    mwc_ground_truth_flux(species, edges, [(:ES_A, :EP_A, k, kr)], 1.0)
end

# ── The gate: a fully-inert inactive contributes L to the denominator ─────────
# `den = Q_A^cat_n + L·1^cat_n`, not `Q_A^cat_n + L·0`. The inactive-state graph
# is step-less (every binding pruned), and the derivation returns `Q_I = 1` for it.
@testset "fully-inert inactive contributes L to the denominator" begin
    inert = @allosteric_mechanism begin
        substrates: S ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + S ⇌ E(S)   :: OnlyA
            E(S) <--> E(P) :: EqualAI
            E + P ⇌ E(P)   :: OnlyA
        end
    end
    fp = ER.fitted_params(inert)      # (:K_A_P_E, :K_A_S_E, :k_ES_to_EP, :L)
    @test fp == (:K_A_P_E, :K_A_S_E, :k_ES_to_EP, :L)

    rng = MersenneTwister(20260714)
    for _ in 1:5
        KS = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng); k = 0.5 + 2rand(rng)
        L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        S = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map fitted_params -> ground-truth params: K_A_S_E=KS, K_A_P_E=KP, k_ES_to_EP=k.
        prm = NamedTuple{(fp..., :Keq, :E_total)}((KP, KS, k, L, Keq, 1.0))
        v_code = real(ER.rate_equation(inert, (S=S, P=P), prm))
        @test isapprox(v_code,
            inert_inactive_flux(; KS=KS, KP=KP, k=k, L=L, Keq=Keq, S=S, P=P); rtol=1e-4)
        # self-validation: L = 0 → reservoir unpopulated → non-allosteric active rate.
        nonallo = (k * S / KS - (k * KP / (Keq * KS)) * P / KP) / (1 + S / KS + P / KP)
        @test isapprox(inert_inactive_flux(; KS=KS, KP=KP, k=k, L=0.0, Keq=Keq, S=S, P=P),
                       nonallo; rtol=1e-4)
    end
end

# ── The gate: a ping-pong :OnlyA I-state must keep a reachable free-enzyme root ──
# A covalent intermediate carries no bound metabolite but does carry a residual.
# Seeding I-state reachability from it makes the pruned inactive graph retain a
# covalent island that free E cannot reach, so no spanning tree rooted at free E
# exists and D[g_free] = 0 — the normalization then divides by zero and
# `rate_equation` is NaN at every concentration. Under formulation 1 only free
# enzyme flips, so a component free E cannot reach holds no inactive mass and
# must be stranded. The mechanism below is accepted by the `:OnlyA`
# thermodynamic guard — it is valid, and the derivation must handle it.
@testset "ping-pong :OnlyA I-state keeps a reachable free-enzyme root" begin
    err1 = @allosteric_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        catalytic_multiplicity: 1
        catalytic_steps: begin
            E + ATP ⇌ E(ATP)                                                      :: EqualAI
            E(ATP) <--> E(F16BP; residual = ATP - F16BP)                          :: OnlyA
            E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP):: EqualAI
            E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)    :: OnlyA
            E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                               :: EqualAI
            E + ADP ⇌ E(ADP)                                                      :: EqualAI
        end
    end
    am = ER.AllostericMechanism(err1)
    @test ER._onlya_haldane_violation(ER.reaction(am), ER.steps(am),
                                      ER.cat_allo_states(am)) === nothing

    _, _, d_free_I = ER._state_rate_polys(am, :I)
    @test !isempty(d_free_I)

    fp = ER.fitted_params(err1)
    prm = NamedTuple{(fp..., :Keq, :E_total)}(((1.3 for _ in fp)..., 3.0, 1.0))
    concs = (ATP = 1.1, F6P = 0.7, ADP = 0.6, F16BP = 0.9)
    @test isfinite(real(ER.rate_equation(err1, concs, prm)))
    @test isfinite(ER._kcat_forward(err1, prm))
end

# ── n-protomer concerted-MWC oracle (formulation 1) ─────────────────────────
# Concerted: all protomers share one conformation. Within a conformation the
# protomers are independent, so the joint occupancy state is a tuple. Only the
# FULLY-unliganded oligomer flips (`freeflip=true`) — the n-protomer extension
# of the free-flip-only model this package derives. `freeflip=false` flips every
# joint state (the classic per-form-flip MWC, formulation 2) and is kept ONLY as
# a discriminator: it must NOT match the derivation.
const OCC = (:E, :EA, :EAB, :EP)

function biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP; k_A, k_I, L, Keq,
                                 A, B, P, FAST=1e7, freeflip=true)
    krA = k_A * kon * KP / (koff * KB * Keq)
    krI = k_I * kon * KP / (koff * KB * Keq)
    prot_edges(kX, krX) = [
        (:E, :EA, kon * A), (:EA, :E, koff),
        (:EA, :EAB, FAST * B / KB), (:EAB, :EA, FAST),
        (:E, :EP, FAST * P / KP), (:EP, :E, FAST),
        (:EAB, :EP, kX), (:EP, :EAB, krX),
    ]
    tbl = Dict(:A => prot_edges(k_A, krA), :I => prot_edges(k_I, krI))
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
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# Self-validation. Check (c) is the load-bearing one: it pins the oracle to
# formulation 1. Checks (a)/(b) pass for BOTH formulations and so cannot
# distinguish them on their own.
@testset "concerted-MWC oligomer oracle self-validation" begin
    rng = MersenneTwister(11)
    for nprot in (1, 2, 3), _ in 1:3
        kon = 0.5+2rand(rng); koff = 0.5+2rand(rng)
        KB = 0.5+2rand(rng); KP = 0.5+2rand(rng)
        kA = 0.5+2rand(rng); kI = 0.5+2rand(rng); Keq = 2.0+2rand(rng)
        A = 0.5+2rand(rng); B = 0.5+2rand(rng); P = 0.5+2rand(rng)
        L = 0.5+rand(rng)
        base = metab_dfree_base_flux(kon, koff, KB, KP, kA, Keq, A, B, P)

        # (a) L = 0 : inactive unpopulated -> nprot x the single-protomer rate.
        @test isapprox(biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=0.0, Keq=Keq, A=A, B=B, P=P),
            nprot * base; rtol=1e-4)

        # (b) k_I = k_A : conformations identical -> L-independent.
        f1 = biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kA, L=L, Keq=Keq, A=A, B=B, P=P)
        f5 = biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kA, L=5.0, Keq=Keq, A=A, B=B, P=P)
        @test isapprox(f1, nprot * base; rtol=1e-4)
        @test isapprox(f1, f5; rtol=1e-4)

        # (d) v = 0 at the equilibrium metabolite ratio.
        Peq = Keq * A * B
        @test abs(biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=L, Keq=Keq, A=A, B=B, P=Peq)) < 1e-6
    end

    # (e) a live inactive conformation must move the flux, or the gate below
    #     would pass without ever exercising the cross term.
    v_live = biuni_mwc_oligomer_flux(2, 1.7, 1.1, 0.8, 0.9;
            k_A=2.5, k_I=0.4, L=0.7, Keq=3.0, A=1.1, B=0.5, P=0.6)
    v_dead = biuni_mwc_oligomer_flux(2, 1.7, 1.1, 0.8, 0.9;
            k_A=2.5, k_I=0.0, L=0.7, Keq=3.0, A=1.1, B=0.5, P=0.6)
    @test !isapprox(v_live, v_dead; rtol=1e-3)

    # (c) THE formulation-1 pin: at nprot = 1 the oracle must reproduce the
    #     established free-flip reference, and must NOT equal the per-form-flip
    #     model. Without this, a formulation-2 oracle would pass (a), (b) and (d)
    #     and then disagree with the derivation by 0.1-3% for live :NonequalAI —
    #     a real number that is not a bug.
    args = (1.7, 1.1, 0.8, 0.9)
    kw = (k_A=2.5, k_I=0.4, L=0.7, Keq=3.0, A=1.1, B=0.5, P=0.6)
    @test isapprox(biuni_mwc_oligomer_flux(1, args...; kw...),
                   biuni_nonequalAI_freeflip_flux(args...; kw...); rtol=1e-4)
    @test !isapprox(biuni_mwc_oligomer_flux(1, args...; freeflip=false, kw...),
                    biuni_nonequalAI_freeflip_flux(args...; kw...); rtol=1e-4)
end

# ── The gate: the ^n combine with a LIVE inactive numerator ─────────────────
# `:OnlyA` always yields a dead inactive cycle (the guard forces an `:OnlyA`
# catalytic tag alongside an `:OnlyA` binding), so the numerator cross term
# `L*N_I*D_I^(n-1)` is live only for `:NonequalAI`. This is the only gate that
# exercises it, and the only mass-action gate at n >= 2 for any family.
@testset ":NonequalAI ^n cross term matches multi-protomer ground truth" begin
    rng = MersenneTwister(20260716)
    for nprot in (2, 3)
        allo = nprot == 2 ?
            @allosteric_mechanism(begin
                substrates: A, B ; products: P ; catalytic_multiplicity: 2
                catalytic_steps: begin
                    E + A <--> E(A)        :: EqualAI
                    E(A) + B ⇌ E(A, B)     :: EqualAI
                    E(A, B) <--> E(P)      :: NonequalAI
                    E + P ⇌ E(P)           :: EqualAI
                end
            end) :
            @allosteric_mechanism(begin
                substrates: A, B ; products: P ; catalytic_multiplicity: 3
                catalytic_steps: begin
                    E + A <--> E(A)        :: EqualAI
                    E(A) + B ⇌ E(A, B)     :: EqualAI
                    E(A, B) <--> E(P)      :: NonequalAI
                    E + P ⇌ E(P)           :: EqualAI
                end
            end)
        fp = ER.fitted_params(allo)
        @test fp == (:kon_A_E, :koff_A_E, :K_P_E, :K_B_EA,
                     :k_A_EAB_to_EP, :k_I_EAB_to_EP, :L)
        for _ in 1:6
            kon = 0.5+2rand(rng); koff = 0.5+2rand(rng)
            KP = 0.5+2rand(rng); KB = 0.5+2rand(rng)
            kA = 0.5+2rand(rng); kI = 0.5+2rand(rng)
            L = 0.5+rand(rng); Keq = 2.0+2rand(rng)
            A = 0.5+2rand(rng); B = 0.5+2rand(rng); P = 0.5+2rand(rng)
            d = Dict(:kon_A_E=>kon, :koff_A_E=>koff, :K_P_E=>KP, :K_B_EA=>KB,
                     :k_A_EAB_to_EP=>kA, :k_I_EAB_to_EP=>kI, :L=>L)
            prm = NamedTuple{(fp..., :Keq, :E_total)}(((d[s] for s in fp)..., Keq, 1.0))
            # `rate_equation` is per active site; the oracle is per oligomer.
            v_code = nprot * real(ER.rate_equation(allo, (A=A, B=B, P=P), prm))
            v_gt = biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=L, Keq=Keq, A=A, B=B, P=P)
            @test isapprox(v_code, v_gt; rtol=1e-4)
        end
    end
end
