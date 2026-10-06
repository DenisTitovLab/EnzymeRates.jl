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
    fp = ER.fitted_params(onlyA)  # (:K_EP_to_E_P, :K_A_ES_to_E_S, :k_A_ES_to_EP, :L)
    @test fp == (:K_EP_to_E_P, :K_A_ES_to_E_S, :k_A_ES_to_EP, :L)

    rng = MersenneTwister(20260713)
    for _ in 1:5
        KA = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng); k = 0.5 + 2rand(rng)
        L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        S = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map fitted_params -> ground-truth params:
        #   K_A_ES_to_E_S=KA, K_EP_to_E_P=KP, k_A_ES_to_EP=k.
        prm = NamedTuple{(fp..., :Keq, :E_total)}((KP, KA, k, L, Keq, 1.0))
        v_code = real(ER.rate_equation(onlyA, (S=S, P=P), prm))
        v_gt = uni_onlyA_flux(KA, KP, k, L=L, Keq=Keq, S=S, P=P)
        @test isapprox(v_code, v_gt; rtol=1e-4)
        @test isfinite(ER._kcat_forward(onlyA, prm))
        # self-validation: L = 0 → inactive conformation unpopulated → the
        # single-conformation (non-allosteric) uni-uni rate.
        kr = k * KP / (Keq * KA)
        nonallo = (k * S / KA - kr * P / KP) / (1 + S / KA + P / KP)
        @test isapprox(uni_onlyA_flux(KA, KP, k, L=0.0, Keq=Keq, S=S, P=P), nonallo;
                       rtol=1e-4)
    end

    # The inactive state neither binds S nor catalyzes, so its free-enzyme graph is a
    # single segment and both states surface D = 1. The genuinely fragmenting
    # cross-weight regime (D_A ≠ D_I, a metabolite-bearing inactive weight) needs a
    # steady-state :OnlyA binding and is covered by the metabolite-bearing-D :OnlyA
    # gate below.
    am = ER.AllostericMechanism(onlyA)
    _, _, dA = ER._state_rate_polys(am, :A)
    _, _, dI = ER._state_rate_polys(am, :I)
    @test dA == ER.poly_one()                       # single active segment
    @test dI == ER.poly_one()                       # single inactive segment

    # kcat stays consistent with rate_equation under normalization.
    prm = NamedTuple{(fp..., :Keq, :E_total)}((0.9, 1.3, 2.1, 0.7, 3.0, 1.0))
    rescaled = ER.rescale_parameter_values(onlyA, prm; scale_k_to_kcat=5.0)
    @test isapprox(ER._kcat_forward(onlyA, rescaled), 5.0; rtol=1e-6)
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
    # (:K_A_EA_to_E_A, :K_EP_to_E_P, :K_A_EAB_to_EA_B, :k_A_EAB_to_EP, :L)
    fp = ER.fitted_params(multiA)
    @test fp == (:K_A_EA_to_E_A, :K_EP_to_E_P, :K_A_EAB_to_EA_B, :k_A_EAB_to_EP, :L)

    rng = MersenneTwister(20260713)
    for _ in 1:5
        KA = 0.5 + 2rand(rng); KB = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng)
        k = 0.5 + 2rand(rng); L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        A = 0.5 + 2rand(rng); B = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map fitted_params -> ground-truth params:
        #   K_A_EA_to_E_A=KA, K_EP_to_E_P=KP, K_A_EAB_to_EA_B=KB, k_A_EAB_to_EP=k.
        prm = NamedTuple{(fp..., :Keq, :E_total)}((KA, KP, KB, k, L, Keq, 1.0))
        v_code = real(ER.rate_equation(multiA, (A=A, B=B, P=P), prm))
        v_gt = multi_onlyA_flux(KA, KB, KP, k, L=L, Keq=Keq, A=A, B=B, P=P)
        @test isapprox(v_code, v_gt; rtol=1e-4)
        @test isfinite(ER._kcat_forward(multiA, prm))
        # self-validation: L = 0 → inactive conformation unpopulated → the
        # single-conformation (non-allosteric) bi-uni rate.
        kr = k * KP / (Keq * KA * KB)
        nonallo = (k*A*B/(KA*KB) - kr*P/KP) / (1 + A/KA + A*B/(KA*KB) + P/KP)
        @test isapprox(multi_onlyA_flux(KA, KB, KP, k, L=0.0, Keq=Keq, A=A, B=B, P=P),
                       nonallo; rtol=1e-4)
    end
end

# ── Metabolite-bearing-D ordered bi-uni :OnlyA network (the LDH i-state case) ──
# S + B ⇌ P. S binds :OnlyA via a STEADY-STATE step (E + S <--> E(S)); B binds
# :EqualAI at rapid equilibrium on E(S); catalysis E(S,B) <--> E(P) is
# steady-state :OnlyA (inactive rate k_I = 0); P binds :EqualAI at rapid
# equilibrium. Because B binds rapidly on the steady-state catalytic path, the
# active free-enzyme spanning-tree weight D[g_free] CARRIES the metabolite B:
# D_A = koff + k·B/K_B. The inactive conformation runs no catalysis and never
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

# ── Single-conformation (non-allosteric) reference for the ordered bi-uni ─────
# The active mechanism alone (E, E(S), E(S,B), E(P)); no inactive conformation,
# no flips. The L → 0 limit of the :OnlyA flux must equal this rate.
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
    # (:K_EP_to_E_P,:k_A_E_S_to_ES,:k_A_ES_to_E_S,:k_A_EBS_to_EP,:K_EBS_to_ES_B,:L)
    fp = ER.fitted_params(metabD)
    @test fp == (:K_EP_to_E_P, :k_A_E_S_to_ES, :k_A_ES_to_E_S, :k_A_EBS_to_EP,
                 :K_EBS_to_ES_B, :L)

    rng = MersenneTwister(20260713)
    for _ in 1:5
        kon = 0.5 + 2rand(rng); koff = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng)
        KB = 0.5 + 2rand(rng); k = 0.5 + 2rand(rng); L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        S = 0.5 + 2rand(rng); B = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map fitted_params -> ground-truth params:
        #   k_A_E_S_to_ES=kon, k_A_ES_to_E_S=koff, K_EBS_to_ES_B=KB, K_EP_to_E_P=KP,
        #   k_A_EBS_to_EP=k.
        prm = NamedTuple{(fp..., :Keq, :E_total)}((KP, kon, koff, k, KB, L, Keq, 1.0))
        v_code = real(ER.rate_equation(metabD, (S=S, B=B, P=P), prm))
        v_gt = metab_dfree_onlyA_flux(kon, koff, KB, KP, k, L=L, Keq=Keq, S=S, B=B, P=P)
        @test isapprox(v_code, v_gt; rtol=1e-4)
        @test isfinite(ER._kcat_forward(metabD, prm))
        # self-validation: L = 0 → inactive conformation unpopulated → the
        # single-conformation (non-allosteric) ordered bi-uni rate.
        @test isapprox(
            metab_dfree_onlyA_flux(kon, koff, KB, KP, k, L=0.0, Keq=Keq, S=S, B=B, P=P),
            metab_dfree_base_flux(kon, koff, KB, KP, k, Keq, S, B, P); rtol=1e-4)
    end

    # B = 0 puts the metabolite B in D[g_free] at zero concentration. The active
    # steady-state S binding leaves a bare koff term in D[g_free], so a reverse path
    # from product back to free enzyme stays open and the reverse flux survives.
    kon, koff, KP, KB, k, L, Keq, S, P = 1.7, 1.1, 0.9, 0.8, 2.1, 0.7, 3.0, 1.1, 0.6
    prm = NamedTuple{(fp..., :Keq, :E_total)}((KP, kon, koff, k, KB, L, Keq, 1.0))
    v0 = real(ER.rate_equation(metabD, (S=S, B=0.0, P=P), prm))
    @test isapprox(v0,
        metab_dfree_onlyA_flux(kon, koff, KB, KP, k; L=L, Keq=Keq, S=S, B=0.0, P=P);
        rtol=1e-4)
    @test v0 < 0    # reverse flux survives at B=0 (invalid :EqualAI-cat trap→0 was fake)
    @test isfinite(ER._kcat_forward(metabD, prm))
end

# ── Free-flip-only ordered bi-uni :NonequalAI-catalysis reference ────────────
# A + B ⇌ P. A binds :EqualAI via a STEADY-STATE step (E + A <--> E(A), kon·A /
# koff); B binds :EqualAI at rapid equilibrium on E(A); catalysis E(A,B) <--> E(P)
# is STEADY-STATE :NonequalAI (rate k_A in the active conformation, k_I in the
# inactive); P binds :EqualAI at rapid equilibrium. Every ligand is :EqualAI
# (K_I = K_A).
#
# Formulation-1 reference: only the free enzyme flips conformation. Each
# conformation runs its own catalytic cycle, coupled only through the shared
# free-enzyme pool — the model this package derives (commit-when-free). Reverse
# catalysis kr from the Haldane relation. `biuni_mwc_oligomer_flux(1, …;
# freeflip=false)` is the per-form-flip model with the same forms and edges, which
# also flips every catalytic intermediate.
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

# ── :NonequalAI-catalysis harness self-validation ────────────────────────────
# First ground the hand-written mass-action reference `metab_dfree_base_flux`
# (the single-conformation ordered bi-uni, reused as the non-allosteric rate at
# k = k_A) against the ODE-validated non-allosteric `rate_equation` — an
# independent, separately ODE-cross-checked code path. This closes the trust loop
# for the whole `mwc_ground_truth_flux` harness. Then confirm the free-flip-only
# reference degenerates correctly: (a) L = 0 (inactive unpopulated) → the base rate
# at k = k_A; (b) k_I = k_A (conformations identical) → the base rate, independent
# of L.
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
        bd = Dict(:k_E_A_to_EA=>kon, :k_EA_to_E_A=>koff, :K_EAB_to_EA_B=>KB,
                  :K_EP_to_E_P=>KP, :k_EAB_to_EP=>kA)
        bprm = NamedTuple{(bfp..., :Keq, :E_total)}(((bd[s] for s in bfp)..., Keq, 1.0))
        @test isapprox(base_rate,
            real(ER.rate_equation(base, (A=A, B=B, P=P), bprm)); rtol=1e-4)

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
#   (k1f·A/k1r)·(k2f/(k2r·P))·(k3f·B/k3r)·(k4f/(k4r·Q)) = 1
# at the equilibrium ratio P·Q/(A·B) = Keq, i.e. the overall Haldane relation
#   k1f·k2f·k3f·k4f = Keq · k1r·k2r·k3r·k4r,
# matching the numerator of Segel Eq. IX-140 (k1f·k2f·k3f·k4f·A·B −
# k1r·k2r·k3r·k4r·P·Q). The two half-reactions' equilibrium constants are NOT
# separately fixed — only their product is — so exactly one reverse constant per
# conformation is dependent. `k1r` (shared, :EqualAI) closes the active cycle;
# `k2r_I` then closes the inactive one against that same `k1r`.
function pingpong_nonequalAI_freeflip_flux(k1f, k3f, k3r;
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
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# `rate_ping_pong_bi_bi` in `test/mechanism_definitions_for_test_enzyme_derivation.jl`
# transcribes this same Segel formula, and both transcriptions are live. Keep them
# independent rather than sharing one: a shared transcription error would green this
# gate and that one at once, whereas two independent transcriptions cross-check each
# other. Sharing would also couple this gate to the MECHANISM_TEST_SPECS fixture,
# whose copy takes `(params::NamedTuple, concs::NamedTuple)` with an `Etotal` rather
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
        f0 = pingpong_nonequalAI_freeflip_flux(k1f, k3f, k3r; act..., ina...,
            L=0.0, Keq=Keq, A=A, B=B, P=P, Q=Q)
        @test isapprox(f0, segel_pingpong_flux(k1f, k1r, k2f_A, k2r_A,
            k3f, k3r, k4f_A, k4r_A, A, B, P, Q); rtol=1e-9)

        # (b) identical conformations → the active-only rate, independent of L.
        same = (k2f_I=k2f_A, k4f_I=k4f_A, k4r_I=k4r_A)
        fe = pingpong_nonequalAI_freeflip_flux(k1f, k3f, k3r; act..., same...,
            L=L, Keq=Keq, A=A, B=B, P=P, Q=Q)
        fe5 = pingpong_nonequalAI_freeflip_flux(k1f, k3f, k3r; act..., same...,
            L=5.0, Keq=Keq, A=A, B=B, P=P, Q=Q)
        @test isapprox(fe, f0; rtol=1e-9)
        @test isapprox(fe, fe5; rtol=1e-9)

        # (c) v = 0 at the equilibrium metabolite ratio P·Q/(A·B) = Keq. Both
        #     conformations are live and unequal, so this pins both Haldanes.
        Qeq = Keq * A * B / P
        @test abs(pingpong_nonequalAI_freeflip_flux(k1f, k3f, k3r; act...,
            ina..., L=L, Keq=Keq, A=A, B=B, P=P, Q=Qeq)) < 1e-9

        # (d) FAST-invariance: the free-enzyme flip is the only cut between the
        #     conformation subnetworks, so it carries zero net flux and E_I/E_A is
        #     exactly L for any FAST. The fast-flip limit is therefore exact here.
        v_live = pingpong_nonequalAI_freeflip_flux(k1f, k3f, k3r; act..., ina...,
            L=L, Keq=Keq, A=A, B=B, P=P, Q=Q)
        for fast in (1e2, 1e12)
            @test isapprox(v_live, pingpong_nonequalAI_freeflip_flux(k1f, k3f,
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
    @test fp == (:k_E_A_to_EA, Symbol("k_A_E_Q_to_EB_res_+A_-P"),
                 Symbol("k_A_EB_res_+A_-P_to_E_Q"), Symbol("k_A_E_res_+A_-P_P_to_EA"),
                 Symbol("k_A_EA_to_E_res_+A_-P_P"),
                 Symbol("k_E_res_+A_-P_B_to_EB_res_+A_-P"),
                 Symbol("k_EB_res_+A_-P_to_E_res_+A_-P_B"),
                 Symbol("k_I_EB_res_+A_-P_to_E_Q"), Symbol("k_I_E_res_+A_-P_P_to_EA"),
                 Symbol("k_I_EA_to_E_res_+A_-P_P"), :L)

    rng = MersenneTwister(20260716)
    for _ in 1:6
        k1f = 0.5+2rand(rng); k3f = 0.5+2rand(rng); k3r = 0.5+2rand(rng)
        k2f_A = 0.5+2rand(rng); k2r_A = 0.5+2rand(rng)
        k4f_A = 0.5+2rand(rng); k4r_A = 0.5+2rand(rng)
        k2f_I = 0.5+2rand(rng); k4f_I = 0.5+2rand(rng); k4r_I = 0.5+2rand(rng)
        L = 0.5+rand(rng); Keq = 2.0+2rand(rng)
        A = 0.5+2rand(rng); B = 0.5+2rand(rng); P = 0.5+2rand(rng); Q = 0.5+2rand(rng)
        # Map fitted_params -> ground-truth params. A step that releases a
        # product while changing the residual is stored as the binding it
        # reverses, F + P → E(A), however the step is written; its rate constants
        # are named by its two sides, so `k_…_EA_to_E_res_…_P` is the
        # product-releasing rate and `k_…_E_res_…_P_to_EA` the product-rebinding
        # one — the binding steps read the usual way round.
        #   k_E_A_to_EA=k1f                          (E + A ⇌ EA, shared)
        #   k_E_res_+A_-P_B_to_EB_res_+A_-P=k3f,
        #   k_EB_res_+A_-P_to_E_res_+A_-P_B=k3r      (F + B ⇌ FB, shared)
        #   k_A_EA_to_E_res_+A_-P_P=k2f_A,
        #   k_A_E_res_+A_-P_P_to_EA=k2r_A            (EA ⇌ F + P, active)
        #   k_I_EA_to_E_res_+A_-P_P=k2f_I,
        #   k_I_E_res_+A_-P_P_to_EA=k2r_I            (EA ⇌ F + P, inactive)
        #   k_A_EB_res_+A_-P_to_E_Q=k4f_A,
        #   k_A_E_Q_to_EB_res_+A_-P=k4r_A            (FB ⇌ E + Q, active)
        #   k_I_EB_res_+A_-P_to_E_Q=k4f_I,
        #   k_I_E_Q_to_EB_res_+A_-P=k4r_I            (FB ⇌ E + Q, inactive)
        # `k_EA_to_E_A` and `k_I_E_Q_to_EB_res_+A_-P` are absent from
        # fitted_params: each conformation's Haldane makes one constant dependent.
        # The oracle takes k4r_I and derives k2r_I by the inactive Haldane, so
        # k2r_I is computed here the same way, and the derivation rederives k4r_I.
        k1r = k1f * k2f_A * k3f * k4f_A / (Keq * k2r_A * k3r * k4r_A)
        k2r_I = k1f * k2f_I * k3f * k4f_I / (Keq * k1r * k3r * k4r_I)
        d = Dict(:k_E_A_to_EA => k1f,
                 Symbol("k_E_res_+A_-P_B_to_EB_res_+A_-P") => k3f,
                 Symbol("k_EB_res_+A_-P_to_E_res_+A_-P_B") => k3r,
                 Symbol("k_A_EA_to_E_res_+A_-P_P") => k2f_A,
                 Symbol("k_A_E_res_+A_-P_P_to_EA") => k2r_A,
                 Symbol("k_A_EB_res_+A_-P_to_E_Q") => k4f_A,
                 Symbol("k_A_E_Q_to_EB_res_+A_-P") => k4r_A,
                 Symbol("k_I_EA_to_E_res_+A_-P_P") => k2f_I,
                 Symbol("k_I_E_res_+A_-P_P_to_EA") => k2r_I,
                 Symbol("k_I_EB_res_+A_-P_to_E_Q") => k4f_I,
                 Symbol("k_I_E_Q_to_EB_res_+A_-P") => k4r_I,
                 :L => L)
        prm = NamedTuple{(fp..., :Keq, :E_total)}(((d[s] for s in fp)..., Keq, 1.0))
        v_code = real(ER.rate_equation(allo, (A=A, B=B, P=P, Q=Q), prm))
        v_gt = pingpong_nonequalAI_freeflip_flux(k1f, k3f, k3r;
            k2f_A=k2f_A, k2r_A=k2r_A, k4f_A=k4f_A, k4r_A=k4r_A,
            k2f_I=k2f_I, k4f_I=k4f_I, k4r_I=k4r_I,
            L=L, Keq=Keq, A=A, B=B, P=P, Q=Q)
        @test isapprox(v_code, v_gt; rtol=1e-10)
    end
end

# ── Ping-pong network with rapid-equilibrium chemistry, at zero products ──────
# The undecorated ping-pong: E + A ⇌ EA → F(P) ⇌ F + P, F + B ⇌ FB ⇌ EQ ⇌ E + Q,
# with both chemistry steps :OnlyA and every binding :EqualAI. Only EA → F(P) is
# steady state; the second chemistry FB ⇌ EQ is rapid equilibrium, so free E and
# the covalent F lie in one active rapid-equilibrium segment. Rapid steps use FAST.
# Only free E flips (formulation 1, E_I/E_A = L). The inactive conformation runs no
# chemistry, so it holds E_I, EA_I and EQ_I; its covalent forms, which free E cannot
# reach, hold no mass. At P = Q = 0 the P and Q releases are one-way and the
# chemistry's reverse never fires, so every form past EA → F(P) drains back to E,
# and the rate is k·(A/KA)/((1 + A/KA)(1 + L)), free of B.
function pingpong_re_chemistry_flux(; KA, KB, KP, KQ, K2, k, L, A, B, FAST=1e9)
    species = [:E_A, :EA_A, :FP_A, :F_A, :FB_A, :EQ_A, :E_I, :EA_I, :EQ_I]
    edges = [
        (:E_A, :EA_A, FAST * A), (:EA_A, :E_A, FAST * KA),   # E + A ⇌ EA (RE)
        (:EA_A, :FP_A, k),                                   # EA → F(P) (SS)
        (:FP_A, :F_A, FAST * KP),                            # F(P) → F + P (RE, P = 0)
        (:F_A, :FB_A, FAST * B), (:FB_A, :F_A, FAST * KB),   # F + B ⇌ FB (RE)
        (:FB_A, :EQ_A, FAST * K2), (:EQ_A, :FB_A, FAST),     # FB ⇌ EQ (RE)
        (:EQ_A, :E_A, FAST * KQ),                            # EQ → E + Q (RE, Q = 0)
        (:E_I, :EA_I, FAST * A), (:EA_I, :E_I, FAST * KA),   # inactive A binding
        (:EQ_I, :E_I, FAST * KQ),                            # inactive Q release
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),          # free-enzyme flip, ratio L
    ]
    mwc_ground_truth_flux(species, edges, [(:EA_A, :FP_A, k, 0.0)], 1.0)
end

# ── Ping-pong rapid-equilibrium-chemistry harness self-validation ─────────────
# The network against its closed form, live (L = 3) and unpopulated (L = 0), across
# four decades of B: the B-binding rate FAST·B stays far above k, so the
# rapid-equilibrium limit holds to about k/(FAST·B).
@testset "ping-pong RE-chemistry ground-truth harness self-validation" begin
    p = (KA=0.7, KB=1.3, KP=0.9, KQ=1.1, K2=2.0, k=1.7, A=1.5)
    closed(L) = p.k * (p.A / p.KA) / ((1 + p.A / p.KA) * (1 + L))
    for L in (0.0, 3.0), B in (1e-2, 1.0, 1e2)
        @test isapprox(pingpong_re_chemistry_flux(; p..., L=L, B=B), closed(L); rtol=1e-5)
    end
end

# ── The gate: a rapid-equilibrium segment holding free E and F ────────────────
# The two ping-pong value gates above write every step steady state, and the two
# ping-pong :OnlyA gates below, whose active segments do hold free E and F, check
# finiteness, the equilibrium ratio and the L = 0 limit only. This gate compares an
# L > 0 rate against ground truth for a mechanism whose active rapid-equilibrium
# segment holds both free E and F. Clearing that segment's polynomial gives free E the
# weight B (F sits at Q/B relative to E), so the inactive term must carry the same
# factor beside L. The ground truth is independent of B (0.28977 at these parameters).
@testset "ping-pong MWC derivation with free E and F in one RE segment" begin
    allo = @allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)                                          :: EqualAI
            E(A) <--> E(P; residual = A - P)                      :: OnlyA
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
            E(B; residual = A - P) ⇌ E(Q)                         :: OnlyA
            E + Q ⇌ E(Q)                                          :: EqualAI
        end
    end
    fp = ER.fitted_params(allo)
    @test fp == (:K_EA_to_E_A, :K_EQ_to_E_Q, Symbol("k_A_EA_to_EP_res_+A_-P"),
                 Symbol("K_A_EB_res_+A_-P_to_EQ"),
                 Symbol("K_EB_res_+A_-P_to_E_res_+A_-P_B"),
                 Symbol("K_EP_res_+A_-P_to_E_res_+A_-P_P"), :L)
    p = (KA=0.7, KB=1.3, KP=0.9, KQ=1.1, K2=2.0, k=1.7, L=3.0, A=1.5)
    # Map fitted_params -> ground-truth params. Each K is the ratio of its to-side to
    # its from-side, so the binding Ks are dissociation constants:
    #   K_EA_to_E_A=KA, K_EQ_to_E_Q=KQ, K_EB_res_…_to_E_res_…_B=KB,
    #   K_EP_res_…_to_E_res_…_P=KP, K_A_EB_res_…_to_EQ=K2 ([EQ]/[FB]),
    #   k_A_EA_to_EP_res_…=k. At P = Q = 0 only KA, k and L enter the rate.
    d = Dict(:K_EA_to_E_A => p.KA, :K_EQ_to_E_Q => p.KQ,
             Symbol("k_A_EA_to_EP_res_+A_-P") => p.k,
             Symbol("K_A_EB_res_+A_-P_to_EQ") => p.K2,
             Symbol("K_EB_res_+A_-P_to_E_res_+A_-P_B") => p.KB,
             Symbol("K_EP_res_+A_-P_to_E_res_+A_-P_P") => p.KP, :L => p.L)
    prm = NamedTuple{(fp..., :Keq, :E_total)}(((d[s] for s in fp)..., 2.0, 1.0))
    v_code(B) = real(ER.rate_equation(allo, (A=p.A, B=B, P=0.0, Q=0.0), prm))
    v_gt(B) = pingpong_re_chemistry_flux(; p..., B=B)
    for B in (1e-3, 0.1, 1.0, 10.0, 1e3)
        @test isapprox(v_code(B), v_gt(B); rtol=1e-5)
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
    fp = ER.fitted_params(inert)      # (:K_A_EP_to_E_P, :K_A_ES_to_E_S, :k_ES_to_EP, :L)
    @test fp == (:K_A_EP_to_E_P, :K_A_ES_to_E_S, :k_ES_to_EP, :L)

    rng = MersenneTwister(20260714)
    for _ in 1:5
        KS = 0.5 + 2rand(rng); KP = 0.5 + 2rand(rng); k = 0.5 + 2rand(rng)
        L = 0.5 + rand(rng); Keq = 2.0 + 2rand(rng)
        S = 0.5 + 2rand(rng); P = 0.5 + 2rand(rng)
        # Map fitted_params -> ground-truth params:
        #   K_A_ES_to_E_S=KS, K_A_EP_to_E_P=KP, k_ES_to_EP=k.
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
    # The constructor errors on any `_onlya_haldane_violation`, so building `am`
    # proves the guard accepts it.
    am = ER.AllostericMechanism(err1)

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
# joint state (the classic per-form-flip MWC, formulation 2) and is kept ONLY as a
# self-validation reference: it must NOT match the derivation. At nprot = 1 the two
# differ by ~0.1–3% whenever k_A ≠ k_I, because per-form flipping routes turnover
# through the faster conformation.
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
        base = metab_dfree_base_flux(kon, koff, KB, KP, kA, Keq, A, B, P)

        # (a) L = 0 : inactive unpopulated -> nprot x the single-protomer rate.
        @test isapprox(biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=0.0, Keq=Keq, A=A, B=B, P=P, freeflip),
            nprot * base; rtol=1e-4)

        # (b) k_I = k_A : conformations identical -> L-independent.
        f1 = biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kA, L=L, Keq=Keq, A=A, B=B, P=P, freeflip)
        f5 = biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kA, L=5.0, Keq=Keq, A=A, B=B, P=P, freeflip)
        @test isapprox(f1, nprot * base; rtol=1e-4)
        @test isapprox(f1, f5; rtol=1e-4)

        # (d) v = 0 at the equilibrium metabolite ratio.
        Peq = Keq * A * B
        @test abs(biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=L, Keq=Keq, A=A, B=B, P=Peq, freeflip)) < 1e-6
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

# ── The gate: :NonequalAI catalysis for 1, 2 and 3 protomers ──────────────────
# The productive-inactive case. Catalysis is :NonequalAI (k_A ≠ k_I) with shared
# graph topology, so D_A and D_I differ only by the catalytic rate constant. Under
# formulation 1 (the enzyme commits to one conformation per catalytic cycle) the
# two conformations normalize per-state, so the derivation cross-weights them. The
# gate goes red if the derivation reverts to the raw Q_A + L·Q_I combine (which
# matches a per-form-flip model, not formulation 1) or mis-renders the
# normalization. At nprot = 1 the oligomer oracle is the free-flip-only reference
# `biuni_nonequalAI_freeflip_flux` (pinned in the oligomer self-validation).
#
# The ^n combine also needs a LIVE inactive numerator: `:OnlyA` always yields a dead
# inactive cycle (the guard forces an `:OnlyA` catalytic tag alongside an `:OnlyA`
# binding), so the numerator cross term `L*N_I*D_I^(n-1)` is live only for
# `:NonequalAI`. This is the only gate that exercises it, and the only mass-action
# gate at n >= 2 for any family.
@testset ":NonequalAI MWC derivation matches mass-action ground truth, 1-3 protomers" begin
    rng = MersenneTwister(20260716)
    for nprot in (1, 2, 3)
        allo = nprot == 1 ?
            @allosteric_mechanism(begin
                substrates: A, B ; products: P ; catalytic_multiplicity: 1
                catalytic_steps: begin
                    E + A <--> E(A)        :: EqualAI
                    E(A) + B ⇌ E(A, B)     :: EqualAI
                    E(A, B) <--> E(P)      :: NonequalAI
                    E + P ⇌ E(P)           :: EqualAI
                end
            end) : nprot == 2 ?
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
        @test fp == (:k_E_A_to_EA, :k_EA_to_E_A, :K_EP_to_E_P, :K_EAB_to_EA_B,
                     :k_A_EAB_to_EP, :k_I_EAB_to_EP, :L)
        # Map fitted_params -> ground-truth params:
        #   k_E_A_to_EA=kon, k_EA_to_E_A=koff, K_EAB_to_EA_B=KB, K_EP_to_E_P=KP,
        #   k_A_EAB_to_EP=k_A, k_I_EAB_to_EP=k_I.
        function params_for(kon, koff, KP, KB, kA, kI, L, Keq)
            d = Dict(:k_E_A_to_EA=>kon, :k_EA_to_E_A=>koff, :K_EP_to_E_P=>KP,
                     :K_EAB_to_EA_B=>KB, :k_A_EAB_to_EP=>kA, :k_I_EAB_to_EP=>kI, :L=>L)
            NamedTuple{(fp..., :Keq, :E_total)}(((d[s] for s in fp)..., Keq, 1.0))
        end
        for _ in 1:6
            kon = 0.5+2rand(rng); koff = 0.5+2rand(rng)
            KP = 0.5+2rand(rng); KB = 0.5+2rand(rng)
            kA = 0.5+2rand(rng); kI = 0.5+2rand(rng)
            L = 0.5+rand(rng); Keq = 2.0+2rand(rng)
            A = 0.5+2rand(rng); B = 0.5+2rand(rng); P = 0.5+2rand(rng)
            prm = params_for(kon, koff, KP, KB, kA, kI, L, Keq)
            # `rate_equation` is per active site; the oracle is per oligomer.
            v_code = nprot * real(ER.rate_equation(allo, (A=A, B=B, P=P), prm))
            v_gt = biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=L, Keq=Keq, A=A, B=B, P=P)
            @test isapprox(v_code, v_gt; rtol=1e-4)
        end

        if nprot == 1
            # B = 0 puts the metabolite B in D[g_free] at zero concentration. The
            # active steady-state A binding leaves a bare koff term in D[g_free], so a
            # reverse path from product back to free enzyme stays open and the reverse
            # flux survives.
            kon, koff, KP, KB = 1.7, 1.1, 0.9, 0.8
            kA, kI, L, Keq, A, P = 2.5, 0.4, 0.7, 3.0, 1.1, 0.9
            vN = real(ER.rate_equation(allo, (A=A, B=0.0, P=P),
                                       params_for(kon, koff, KP, KB, kA, kI, L, Keq)))
            @test isapprox(vN, biuni_nonequalAI_freeflip_flux(kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=L, Keq=Keq, A=A, B=0.0, P=P); rtol=1e-4)
            @test abs(vN) > 1e-3
        end
    end
end

# ── The gate: dead-inactive ping-pong :OnlyA (all chemical steps :OnlyA) ──────
# A ping-pong bi-bi where the F6P substrate binding is :OnlyA and BOTH chemical
# (isomerisation) steps are :OnlyA. That combination makes the inactive
# conformation catalytically dead: dropping the :OnlyA iso groups strands the
# covalent-intermediate branch, so the inactive graph collapses to the free
# enzyme plus its rapid-equilibrium ATP and ADP bindings. Its free-enzyme weight
# is therefore the bare constant d_free_I = 1 and it contributes only L·Q_I (no
# flux) to the denominator. The multi-chemical-step :OnlyA form the enumeration
# must feed: the derivation admits it, keeps a finite rate as a product goes to
# zero (no covalent sink), and yields a finite positive kcat.
@testset "dead-inactive ping-pong :OnlyA derivation" begin
    dead = @allosteric_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        catalytic_multiplicity: 1
        catalytic_steps: begin
            E + ATP ⇌ E(ATP)                                                      :: EqualAI
            E(ATP) <--> E(F16BP; residual = ATP - F16BP)                          :: OnlyA
            E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP):: EqualAI
            E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)    :: OnlyA
            E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                               :: OnlyA
            E + ADP ⇌ E(ADP)                                                      :: EqualAI
        end
    end
    afp = ER.fitted_params(dead)
    @test afp == (:K_EADP_to_E_ADP, :K_EATP_to_E_ATP,
                  Symbol("k_A_EATP_to_EF16BP_res_+ATP_-F16BP"),
                  Symbol("K_A_EF6P_res_+ATP_-F16BP_to_EADP"),
                  Symbol("K_EF16BP_res_+ATP_-F16BP_to_E_res_+ATP_-F16BP_F16BP"),
                  Symbol("K_A_EF6P_res_+ATP_-F16BP_to_E_res_+ATP_-F16BP_F6P"), :L)

    # The constructor errors on any `_onlya_haldane_violation`, so building `am`
    # proves the guard accepts it.
    am = ER.AllostericMechanism(dead)
    _, _, d_free_I = ER._state_rate_polys(am, :I)
    @test ER._poly_to_expr(d_free_I, Set{Symbol}(), Set{Symbol}()) == 1

    # Non-allosteric twin: the SAME six steps with no allosteric tags. At L = 0 the
    # inactive conformation is unpopulated, so the allosteric rate must reduce to
    # this twin's rate — an independent re-derivation through the non-allosteric
    # King–Altman path. Map allo params -> twin params: the three :EqualAI bindings
    # (K_EADP_to_E_ADP, K_EATP_to_E_ATP, K_EF16BP_res_…) share names; the three
    # :OnlyA A-tagged params drop the "A_" tag (k_A_EATP…→k_EATP…, and the two
    # K_A_EF6P…→K_EF6P…).
    nonallo = @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        steps: begin
            E + ATP ⇌ E(ATP)
            E(ATP) <--> E(F16BP; residual = ATP - F16BP)
            E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP)
            E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)
            E(F6P; residual = ATP - F16BP) ⇌ E(ADP)
            E + ADP ⇌ E(ADP)
        end
    end
    nfp = ER.fitted_params(nonallo)
    allo_to_twin = Dict(
        Symbol("k_A_EATP_to_EF16BP_res_+ATP_-F16BP") =>
            Symbol("k_EATP_to_EF16BP_res_+ATP_-F16BP"),
        Symbol("K_A_EF6P_res_+ATP_-F16BP_to_EADP") =>
            Symbol("K_EF6P_res_+ATP_-F16BP_to_EADP"),
        Symbol("K_A_EF6P_res_+ATP_-F16BP_to_E_res_+ATP_-F16BP_F6P") =>
            Symbol("K_EF6P_res_+ATP_-F16BP_to_E_res_+ATP_-F16BP_F6P"))
    twin_of(s) = get(allo_to_twin, s, s)

    rng = MersenneTwister(20260717)
    for _ in 1:6
        vals = Dict(s => 0.5 + 2rand(rng) for s in nfp)
        Keq = 2.0 + 2rand(rng)
        ATP = 0.5+2rand(rng); F6P = 0.5+2rand(rng)
        ADP = 0.5+2rand(rng); F16BP = 0.5+2rand(rng)
        concs = (ATP=ATP, F6P=F6P, ADP=ADP, F16BP=F16BP)
        prm = NamedTuple{(afp..., :Keq, :E_total)}(
            ((s === :L ? 0.7 : vals[twin_of(s)] for s in afp)..., Keq, 1.0))

        # rate finite at a normal point AND as F16BP -> 0 (the no-sink check: a
        # dead inactive conformation must not strand enzyme in a covalent form).
        @test isfinite(real(ER.rate_equation(dead, concs, prm)))
        @test isfinite(real(ER.rate_equation(dead,
            (ATP=ATP, F6P=F6P, ADP=ADP, F16BP=1e-8), prm)))
        kcat = ER._kcat_forward(dead, prm)
        @test isfinite(kcat) && kcat > 0

        # v = 0 at the equilibrium metabolite ratio ADP·F16BP = Keq·ATP·F6P.
        @test isapprox(real(ER.rate_equation(dead,
            (ATP=ATP, F6P=F6P, ADP=Keq*ATP*F6P/F16BP, F16BP=F16BP), prm)),
            0.0; atol=1e-8)

        # L = 0 cross-check: the allosteric rate reduces to the non-allosteric twin.
        pn = NamedTuple{(nfp..., :Keq, :E_total)}(((vals[s] for s in nfp)..., Keq, 1.0))
        p0 = NamedTuple{(afp..., :Keq, :E_total)}(
            ((s === :L ? 0.0 : vals[twin_of(s)] for s in afp)..., Keq, 1.0))
        @test isapprox(real(ER.rate_equation(dead, concs, p0)),
                       real(ER.rate_equation(nonallo, concs, pn)); rtol=1e-4)
    end
end
