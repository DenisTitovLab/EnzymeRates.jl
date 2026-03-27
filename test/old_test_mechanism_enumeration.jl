# ABOUTME: Tests for the staged mechanism enumeration pipeline
# ABOUTME: Organized by stage with hand-calculated expected values

using Random

# ── Helper: convert EnzymeMechanism → MechanismSpec ──────────

"""
Convert a compiled EnzymeMechanism back to a MechanismSpec.
Reconstructs StepSpec entries from the mechanism's reactions.
"""
function mechanism_spec_from_mechanism(
    m::EnzymeMechanism,
    @nospecialize(rxn::EnzymeReaction))
    rxns = EnzymeRates.reactions(m)
    eq_steps_tuple = EnzymeRates.equilibrium_steps(m)
    pc = EnzymeRates.param_constraints(m)

    steps = EnzymeRates.StepSpec[]
    for (i, (lhs, rhs)) in enumerate(rxns)
        push!(steps, EnzymeRates.StepSpec(
            collect(lhs), collect(rhs),
            eq_steps_tuple[i]))
    end

    constraints = [
        (t, c, [(s, sc) for (s, sc) in f])
        for (t, c, f) in pc
    ]

    EnzymeRates.MechanismSpec(
        rxn, steps, constraints,
        length(parameters(m)))
end

# ── Reaction definitions ─────────────────────────────────────

# --- No regulators ---

const uni_uni = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end

const uni_bi = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
end

const bi_bi = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

const bi_bi_ping_pong = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end

const ter_ter = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
end

# --- Dead-end inhibitor variants ---

const uni_uni_dead_end_I = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    dead_end_inhibitors: I
end

const uni_bi_dead_end_I = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    dead_end_inhibitors: I
end

const bi_bi_dead_end_I = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    dead_end_inhibitors: I
end

const bi_bi_ping_pong_dead_end_I = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    dead_end_inhibitors: I
end

const uni_uni_dead_end_I_J = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    dead_end_inhibitors: I, J
end

# --- Allosteric regulator variants ---

const uni_uni_allosteric_R = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R
end

const uni_bi_allosteric_R = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    allosteric_regulators: R
end

const bi_bi_ping_pong_allosteric_R = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
end

const uni_bi_allosteric_R_cn2 = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    allosteric_regulators: R
end

# --- Mixed ---

const bi_bi_dead_end_I_allosteric_R = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    dead_end_inhibitors: I
    allosteric_regulators: R
end

const bi_bi_allosteric_R1_R2 = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    allosteric_regulators: R1, R2
end

# --- Unknown role (for end-to-end partitioning) ---

const uni_uni_reg_unknown = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    regulators: I
end

const uni_bi_reg_unknown = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    regulators: I
end

const bi_bi_ping_pong_reg_unknown = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    regulators: I
end

# ── Helper: run pipeline stage by stage (all partitions) ─────

function _run_full_pipeline_stages(rxn; catalytic_n::Int=0,
                                   max_re_groups::Int=7)
    catalytic = EnzymeRates._catalytic_topologies(rxn)

    roles = EnzymeRates.regulator_roles(rxn)
    fixed_de = Symbol[r[1] for r in roles
                       if r[2] == :dead_end]
    fixed_allo = Symbol[r[1] for r in roles
                         if r[2] == :allosteric]
    unknown = Symbol[r[1] for r in roles
                      if r[2] == :unknown]
    n_unknown = length(unknown)

    n_ress = 0; n_de = 0; n_eq = 0; n_dd = 0
    n_allo = 0; n_tr = 0; n_allo_dd = 0

    for reg_mask in 0:(1 << n_unknown) - 1
        de = Symbol[fixed_de;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i - 1)) & 1 == 0]]
        al = Symbol[fixed_allo;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i - 1)) & 1 == 1]]

        ress = EnzymeRates._expand_ress_variants(
            catalytic, rxn; max_re_groups)
        n_ress += length(ress)

        with_de = EnzymeRates._expand_dead_end(
            ress, rxn; dead_end_regs=de,
            include_substrate_product=true)
        n_de += length(with_de)

        with_eq =
            EnzymeRates._expand_equivalence_constraints(
                with_de, rxn)
        n_eq += length(with_eq)

        dd = EnzymeRates._deduplicate(with_eq, rxn)
        n_dd += length(dd)

        if !isempty(al)
            cn = catalytic_n > 0 ? catalytic_n : 1
            allo = EnzymeRates._expand_allosteric(
                dd, rxn; catalytic_n=cn,
                allosteric_regs=al)
            n_allo += length(allo)
            allo = EnzymeRates._expand_tr_equivalence(
                allo, rxn)
            n_tr += length(allo)
            allo = EnzymeRates._deduplicate_allosteric(
                allo, rxn)
            n_allo_dd += length(allo)
        end
    end

    (catalytic=length(catalytic), ress=n_ress,
     dead_end=n_de, equivalence=n_eq, dedup=n_dd,
     allosteric=n_allo, tr_equiv=n_tr,
     allosteric_dedup=n_allo_dd)
end

@testset "Mechanism Enumeration Pipeline" begin

@testset "Stage 1: Catalytic topologies" begin

    @testset "Uni-Uni" begin
        topos = EnzymeRates._catalytic_topologies(uni_uni)
        @test length(topos) == 1

        # E ⇌ ES, E ⇌ EP, ES <--> EP
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end

        @test compile_mechanism(topos[1]) === m_uu

        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps
            ) == 1
        end
    end

    @testset "Uni-Bi" begin
        topos = EnzymeRates._catalytic_topologies(uni_bi)
        @test length(topos) == 3

        # Topo 1: ordered release P-first
        # Path: E→ES→EPQ→EQ→E
        m_ub1 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P_Q[AB],
                    E_Q[B], E_S[AB]
            end
            steps: begin
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end

        # Topo 2: random release (both P and Q paths)
        m_ub2 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P[A], E_P_Q[AB],
                    E_Q[B], E_S[AB]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end

        # Topo 3: ordered release Q-first
        # Path: E→ES→EPQ→EP→E
        m_ub3 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P[A],
                    E_P_Q[AB], E_S[AB]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end

        # Round-trip: each hand-defined mechanism matches
        # a compiled topology
        defined = [m_ub1, m_ub2, m_ub3]
        for (i, m) in enumerate(defined)
            compiled = compile_mechanism(topos[i])
            @test compiled === m
        end

        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps
            ) == 1
        end
    end

    @testset "Bi-Bi" begin
        topos = EnzymeRates._catalytic_topologies(bi_bi)
        @test length(topos) == 9

        # Topo 1: sequential bind A-first,
        #         sequential release Q-first
        # Path: E→EA→EAB→EPQ→EQ→E
        m_bb1 = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        @test compile_mechanism(topos[1]) === m_bb1

        # Topo 2: sequential bind B-first,
        #         sequential release P-first
        # Path: E→EB→EAB→EPQ→EP→E
        m_bb2 = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A_B[CN],
                    E_B[N], E_P[C],
                    E_P_Q[CN]
            end
            steps: begin
                [E_B, A] ⇌ [E_A_B]
                [E, B] ⇌ [E_B]
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        @test compile_mechanism(topos[2]) === m_bb2

        # Topo 6: sequential bind A-first,
        #         sequential release P-first
        # Path: E→EA→EAB→EPQ→EP→E
        m_bb6 = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P[C],
                    E_P_Q[CN]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        @test compile_mechanism(topos[6]) === m_bb6

        # Topo 9: random bind, random release (fully
        # random on both sides)
        m_bb9 = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_B, A] ⇌ [E_A_B]
                [E, B] ⇌ [E_B]
                [E_A, B] ⇌ [E_A_B]
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        @test compile_mechanism(topos[9]) === m_bb9

        # Verify structural properties hold for all
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps
            ) == 1
        end

        # Verify 9 topologies decompose into known
        # categories by form count:
        # 5 forms: sequential/sequential (4 topos)
        # 6 forms: random-one-side (4 topos)
        # 7 forms: fully random (1 topo)
        form_counts = [
            length(
                EnzymeRates.enzyme_forms(
                    compile_mechanism(t))
            ) for t in topos
        ]
        @test count(==(5), form_counts) == 4
        @test count(==(6), form_counts) == 4
        @test count(==(7), form_counts) == 1
    end

    @testset "Bi-Bi Ping-Pong" begin
        topos =
            EnzymeRates._catalytic_topologies(
                bi_bi_ping_pong)
        @test length(topos) == 10

        # Topo 6: classic ping-pong with Estar
        # E→EA→Estar(+P)→Estar_B→E(+Q)
        m_pp = @enzyme_mechanism begin
            species: begin
                substrates: A[CX], B[N]
                products: P[C], Q[NX]
                enzymes: E,
                    E_A[CX], E_Q[NX],
                    Estar[X], Estar_A_P[CX],
                    Estar_B[NX]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [Estar, B] ⇌ [Estar_B]
                [E, Q] ⇌ [E_Q]
                [Estar, P] ⇌ [Estar_A_P]
                [E_A] <--> [Estar_A_P]
                [Estar_B] ⇌ [E_Q]
            end
        end
        @test compile_mechanism(topos[6]) === m_pp

        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps
            ) == 1
        end
    end

    @testset "Ter-Ter" begin
        topos = EnzymeRates._catalytic_topologies(
            ter_ter)
        # 3969 = 63 × 63 = (2^(3!) - 1)². Each side
        # (binding/release) has 3!=6 permutation paths
        # through Boolean lattice B_3; all 2^6-1=63
        # non-empty path subsets produce distinct edge
        # sets; sides are independent.
        @test length(topos) == 3969
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps
            ) == 1
        end
    end

    @testset "Ter-Bi" begin
        ter_bi = @enzyme_reaction begin
            substrates: A[C], B[N], D[X]
            products: P[CN], Q[X]
        end
        topos = EnzymeRates._catalytic_topologies(
            ter_bi)
        # 204 = 189 sequential + 15 ping-pong.
        # Sequential: (2^(3!) - 1) × (2^(2!) - 1)
        #   = 63 × 3 = 189.
        # Ping-pong: D[X]→Q[X] can isomerize
        # independently, creating 15 Estar topologies.
        @test length(topos) == 204
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps
            ) == 1
        end
    end

end

@testset "Stage 2: RE/SS expansion" begin
    @testset "Uni-Uni: 2^3 - 1 = 7 variants" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni)
        result = EnzymeRates._expand_ress_variants(
            [topo], uni_uni)
        @test length(result) == 7

        for s in result
            @test any(
                !s2.is_equilibrium for s2 in s.steps
            )
        end
    end

    @testset "Uni-Bi" begin
        # Ordered Q-first release: E→ES→EPQ→EQ→E
        m_ub1 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P_Q[AB],
                    E_Q[B], E_S[AB]
            end
            steps: begin
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_ub1, uni_bi)
        n_steps = length(topo.steps)
        result = EnzymeRates._expand_ress_variants(
            [topo], uni_bi)
        @test length(result) == 2^n_steps - 1

        for s in result
            @test any(
                !s2.is_equilibrium for s2 in s.steps
            )
        end
    end

    @testset "Bi-Bi" begin
        # Sequential A-first bind, Q-first release
        m_bb_seq = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_seq, bi_bi)
        n_steps = length(topo.steps)
        result = EnzymeRates._expand_ress_variants(
            [topo], bi_bi)
        @test length(result) == 2^n_steps - 1

        for s in result
            @test any(
                !s2.is_equilibrium for s2 in s.steps
            )
        end
    end

    @testset "max_re_groups filtering" begin
        # Fully random Bi-Bi
        m_bb_random = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_B, A] ⇌ [E_A_B]
                [E, B] ⇌ [E_B]
                [E_A, B] ⇌ [E_A_B]
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_random, bi_bi)
        result_default = EnzymeRates._expand_ress_variants(
            [topo], bi_bi)
        @test length(result_default) == 511
        result_mrg2 = EnzymeRates._expand_ress_variants(
            [topo], bi_bi; max_re_groups=2)
        @test length(result_mrg2) == 241
        result_mrg3 = EnzymeRates._expand_ress_variants(
            [topo], bi_bi; max_re_groups=3)
        @test length(result_mrg3) == 379
        result_mrg1 = EnzymeRates._expand_ress_variants(
            [topo], bi_bi; max_re_groups=1)
        # 89 = 9 single-SS + 32 two-SS + 48 three-SS
        # where remaining RE steps form 1 connected
        # component in each case
        @test length(result_mrg1) == 89
    end
end

@testset "Stage 3a: Substrate/product dead-end expansion" begin

    @testset "Uni-Uni: passthrough (no off-cycle forms)" begin
        # Uni-Uni: 3 forms (E, E_S, E_P), all on-cycle.
        # E_S has all subs, E_P has all prods, and E
        # can only bind S→E_S or P→E_P (both catalytic).
        # → 0 dead-end forms, passthrough.
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni)
        ress = EnzymeRates._expand_ress_variants(
            [topo], uni_uni)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                ress, uni_uni)
        @test length(result) == length(ress)
    end

    @testset "Bi-Bi random: 4 mixed dead-end forms" begin
        # Fully-random Bi-Bi (7 forms):
        # E, E_A, E_B, E_A_B, E_P, E_Q, E_P_Q
        # Dead-ends require mixed substrate+product.
        # Eligible forms (not all-subs or all-prods):
        #   E: all 4 mets → catalytic (no dead-ends)
        #   E_A: +P→E_A_P(mixed✓), +Q→E_A_Q(mixed✓)
        #   E_B: +P→E_B_P(mixed✓), +Q→E_B_Q(mixed✓)
        #   E_P: +A→E_A_P(mixed✓), +B→E_B_P(mixed✓)
        #   E_Q: +A→E_A_Q(mixed✓), +B→E_B_Q(mixed✓)
        # 4 unique dead-end forms → 2^4 = 16 variants
        m_bb_random = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_B, A] ⇌ [E_A_B]
                [E, B] ⇌ [E_B]
                [E_A, B] ⇌ [E_A_B]
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_random, bi_bi)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [topo], bi_bi)
        @test length(result) == 16
    end

    @testset "Uni-Bi ordered: passthrough (no mixed forms)" begin
        # Uni-Bi ordered Q-first release:
        # E, E_S, E_P_Q, E_Q (4 forms)
        # E+P→E_P would be single-product → rejected
        #   (mixed substrate+product required)
        # E_Q+S→E_S_Q has all subs (S only) → rejected
        # → 0 dead-end forms, passthrough.
        m_ub1 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P_Q[AB],
                    E_Q[B], E_S[AB]
            end
            steps: begin
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_ub1, uni_bi)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [topo], uni_bi)
        @test length(result) == 1
    end

    @testset "Bi-Bi Ping-Pong: dead-ends with Estar" begin
        # Classic ping-pong: E→EA→Estar(+P)→EstarB→E(+Q)
        # Forms: E, E_A, Estar, Estar_A_P, Estar_B, E_Q
        # Mixed dead-end opportunities from E_A:
        #   +P→E_A_P(mixed), +Q→E_A_Q(mixed)
        # From E_Q (has product Q):
        #   +B→E_B_Q(mixed)
        # 3 dead-end forms → 2^3 = 8 variants
        m_pp = @enzyme_mechanism begin
            species: begin
                substrates: A[CX], B[N]
                products: P[C], Q[NX]
                enzymes: E,
                    E_A[CX], E_Q[NX],
                    Estar[X], Estar_A_P[CX],
                    Estar_B[NX]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [Estar, B] ⇌ [Estar_B]
                [E, Q] ⇌ [E_Q]
                [Estar, P] ⇌ [Estar_A_P]
                [E_A] <--> [Estar_A_P]
                [Estar_B] ⇌ [E_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_pp, bi_bi_ping_pong)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [topo], bi_bi_ping_pong)
        @test length(result) == 8
    end
end

@testset "Stage 3b: Regulator dead-end expansion" begin
    @testset "No regulators: passthrough" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni)
        result = EnzymeRates._expand_dead_end(
            [topo], uni_uni; dead_end_regs=Symbol[])
        @test length(result) == 1
    end

    @testset "Uni-Uni + I: eligible forms" begin
        # Catalytic topology same as base Uni-Uni
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni_dead_end_I)
        result = EnzymeRates._expand_dead_end(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        # Only E eligible (ES fully occupied with all
        # substrates, EP fully occupied with all
        # products) → (2^1)^1 = 2 variants
        @test length(result) == 2
    end

    @testset "Uni-Bi + I" begin
        # Ordered Q-first release
        m_ub1 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P_Q[AB],
                    E_Q[B], E_S[AB]
            end
            steps: begin
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_ub1, uni_bi_dead_end_I)
        result = EnzymeRates._expand_dead_end(
            [topo], uni_bi_dead_end_I;
            dead_end_regs=[:I])
        # Q-first release: E, ES, EQ, EPQ (4 forms).
        # ES has all subs, EPQ has all prods
        # → eligible = {E, EQ} = 2 forms
        # Sub/prod dead-ends: 0 (E_P is single-product,
        #   E_S_Q has all subs — both rejected)
        # Regulator dead-ends: E_I__reg1, E_I__reg1_Q
        # Combined: 2 dead-end forms → 2^2 = 4
        @test length(result) == 4
    end

    @testset "Bi-Bi + I" begin
        # Sequential A-first bind, Q-first release
        m_bb_seq = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_seq, bi_bi_dead_end_I)
        result = EnzymeRates._expand_dead_end(
            [topo], bi_bi_dead_end_I;
            dead_end_regs=[:I])
        # Sequential (5 forms: E, EA, EAB,
        # EQ, EPQ). EAB has all subs, EPQ all prods
        # → eligible = {E, EA, EQ} = 3 forms
        # Sub/prod dead-ends: E_A_P, E_A_Q, E_B_Q
        #   (3 mixed forms; E_B and E_P rejected as
        #   single-type, not mixed sub+prod)
        # Regulator dead-ends: E_I__reg1,
        #   E_A_I__reg1, E_I__reg1_Q (3 forms)
        # Combined: 6 dead-end forms → 2^6 = 64
        @test length(result) == 64
    end

    @testset "2 inhibitors: Uni-Uni + I, J" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni_dead_end_I_J)
        result = EnzymeRates._expand_dead_end(
            [topo], uni_uni_dead_end_I_J;
            dead_end_regs=[:I, :J])
        # Only E eligible → (2^2)^1 = 4
        @test length(result) == 4
    end

    @testset "Allosteric-only: passthrough" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni_allosteric_R)
        result = EnzymeRates._expand_dead_end(
            [topo], uni_uni_allosteric_R;
            dead_end_regs=Symbol[])
        @test length(result) == 1
    end

    @testset "Properties" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni_dead_end_I)
        result = EnzymeRates._expand_dead_end(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        for s in result
            if length(s.steps) > length(topo.steps)
                @test s.param_count > topo.param_count
            else
                @test s.param_count == topo.param_count
            end
        end
    end
end

@testset "Stage 4: Equivalence constraints" begin
    @testset "No equiv groups: passthrough" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni)
        result =
            EnzymeRates._expand_equivalence_constraints(
                [topo], uni_uni)
        @test length(result) == 1
    end

    @testset "1 equiv group: 2 variants" begin
        # Bi-Bi sequential: dead-end P-binding at EA
        # creates EA+P→EAP step. Combined with
        # catalytic EQ+P→EPQ, P appears in 2 RE steps
        # → 1 equiv group → 2 variants.
        m_bb_seq = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_seq, bi_bi)
        de = EnzymeRates._expand_dead_end(
            [topo], bi_bi; dead_end_regs=Symbol[])
        # Variant with E_A_P dead-end: P appears in 2
        # RE steps (EQ+P→EPQ, EA+P→EAP) → 1 equiv
        # group → 2 variants
        target = Set([:E, :E_A, :E_A_B, :E_Q,
            :E_P_Q, :E_A_P])
        s = first(
            x for x in de
            if Set(EnzymeRates.all_form_names(x)) ==
                target)
        eq =
            EnzymeRates._expand_equivalence_constraints(
                [s], bi_bi)
        @test length(eq) == 2
    end

    @testset "Multiple equiv groups: multiplicative" begin
        # Bi-Bi sequential: dead-end form E_A_Q binds
        # both A (at E_Q) and Q (at E_A), creating
        # mirror steps. A appears in 2 RE steps and Q
        # in 2 RE steps → 2 equiv groups → 4 variants.
        m_bb_seq = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_seq, bi_bi)
        de = EnzymeRates._expand_dead_end(
            [topo], bi_bi; dead_end_regs=Symbol[])
        # Variant with E_A_Q dead-end: A in 2 RE steps
        # (E+A→EA, EQ+A→EAQ) and Q in 2 RE steps
        # (E+Q→EQ, EA+Q→EAQ) → 2 equiv groups
        target = Set([:E, :E_A, :E_A_B, :E_Q,
            :E_P_Q, :E_A_Q])
        s = first(
            x for x in de
            if Set(EnzymeRates.all_form_names(x)) ==
                target)
        eq = EnzymeRates._expand_equivalence_constraints(
            [s], bi_bi)
        @test length(eq) == 4
    end

    @testset "Substrate/regulator same metabolite" begin
        # Dummy names (:I__reg1) prevent grouping
        # regulator binding steps with catalytic steps
        # that bind a different real metabolite. With
        # Uni-Uni + I, the dead-end spec has :S and
        # :I__reg1 as distinct metabolites → no equiv
        # group → only 1 variant (passthrough).
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni_dead_end_I)
        de = EnzymeRates._expand_dead_end(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        s = de[argmax(
            length(x.steps) for x in de)]
        eq =
            EnzymeRates._expand_equivalence_constraints(
                [s], uni_uni_dead_end_I)
        # Only 1 variant: no equiv groups because
        # each metabolite (S, P, I__reg1) appears in
        # exactly 1 binding step
        @test length(eq) == 1
        @test isempty(eq[1].param_constraints)
    end

    @testset "Substrate/product dead-end equiv" begin
        # Equivalence constraints apply to dead-end
        # complexes from substrate/product binding.
        # Bi-Bi sequential: EA+B→EAB (catalytic).
        # If E also has dead-end B-binding (E+B→E_B),
        # both steps bind B → equiv group.
        m_bb_seq = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_seq, bi_bi)
        de = EnzymeRates._expand_dead_end(
            [topo], bi_bi; dead_end_regs=Symbol[])
        # Find specs with substrate/product dead-ends
        has_de = filter(
            x -> length(x.steps) > length(topo.steps),
            de)
        @test !isempty(has_de)
        s = has_de[argmax(
            length(x.steps) for x in has_de)]
        eq =
            EnzymeRates._expand_equivalence_constraints(
                [s], bi_bi)
        # Should have equivalence variants
        @test length(eq) >= 2
        constrained = filter(
            x -> !isempty(x.param_constraints),
            eq)
        @test !isempty(constrained)
    end

    @testset "Properties" begin
        # Sequential A-first bind, Q-first release
        m_bb_seq = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_seq, bi_bi_dead_end_I)
        de = EnzymeRates._expand_dead_end(
            [topo], bi_bi_dead_end_I;
            dead_end_regs=[:I])
        eq =
            EnzymeRates._expand_equivalence_constraints(
                de, bi_bi_dead_end_I)
        @test length(eq) >= length(de)

        unconstrained = filter(
            s -> isempty(s.param_constraints), eq)
        constrained = filter(
            s -> !isempty(s.param_constraints), eq)
        @test !isempty(constrained)
        @test !isempty(unconstrained)
        # Every constrained spec has fewer params than
        # at least one unconstrained spec with same steps
        for c in constrained
            matching = filter(
                u -> length(u.steps) == length(c.steps),
                unconstrained)
            @test !isempty(matching)
            @test c.param_count <
                maximum(
                    u.param_count for u in matching)
        end
    end
end

@testset "Stage 5: Deduplication" begin
    @testset "Uni-Uni dedup" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni)
        ress = EnzymeRates._expand_ress_variants(
            [topo], uni_uni)
        deduped = EnzymeRates._deduplicate(ress, uni_uni)
        # 3 distinct concentration fingerprints survive:
        #   {1,[S],[P]}: both binding steps RE — classic
        #     reversible Michaelis-Menten denominator
        #   {1,[S]}: only S-binding RE — denominator lacks
        #     product term
        #   {1,[P]}: only P-binding RE — denominator lacks
        #     substrate term
        # These are genuinely different rate equation forms.
        @test length(deduped) == 3
        for d in deduped
            @test d.param_count <=
                minimum(s.param_count for s in ress)
        end
    end

    @testset "Duplicate removal" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni)
        result = EnzymeRates._deduplicate(
            [topo, deepcopy(topo)], uni_uni)
        @test length(result) == 1
    end

    @testset "Bi-Bi ordered: single-SS fingerprints" begin
        # For an ordered Bi-Bi topology, create single-SS
        # variants by toggling each RE step to SS one at a
        # time. The step-based pipeline correctly toggles
        # all steps including isomerization.
        m_bb_seq = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_seq, bi_bi)
        ress = EnzymeRates._expand_ress_variants(
            [topo], bi_bi)
        # Filter to single-SS variants: exactly one SS step
        single_ss = filter(ress) do spec
            count(!s.is_equilibrium for s in spec.steps) == 1
        end
        # All 5 steps toggleable (including isomerization)
        @test length(single_ss) == 5
        deduped = EnzymeRates._deduplicate(
            single_ss, bi_bi)
        # All 5 single-SS variants have distinct fingerprints
        @test length(deduped) == 5
    end

    @testset "Keeps lower param_count" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni)
        topo2 = EnzymeRates.MechanismSpec(
            topo.reaction, topo.steps,
            topo.param_constraints,
            topo.param_count + 5)
        result = EnzymeRates._deduplicate(
            [topo, topo2], uni_uni)
        @test length(result) == 1
        @test result[1].param_count == topo.param_count
    end

    @testset "Bi-Bi random: single-SS dedup" begin
        # Random BiBi has symmetric binding paths.
        # Some single-SS variants produce identical
        # concentration fingerprints → dedup removes them.
        m_bb_random = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_B, A] ⇌ [E_A_B]
                [E, B] ⇌ [E_B]
                [E_A, B] ⇌ [E_A_B]
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_random, bi_bi)
        ress = EnzymeRates._expand_ress_variants(
            [topo], bi_bi)
        single_ss = filter(ress) do spec
            count(!s.is_equilibrium
                for s in spec.steps) == 1
        end
        @test length(single_ss) == 9
        deduped = EnzymeRates._deduplicate(
            single_ss, bi_bi)
        # Symmetric paths collapse: 9 → 5
        @test length(deduped) == 5
    end

    @testset "Uni-Uni + I: regulator dedup" begin
        # Dead-end regulator creates additional forms.
        # Some RE/SS variants produce identical
        # fingerprints → dedup removes them.
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni_dead_end_I)
        ress = EnzymeRates._expand_ress_variants(
            [topo], uni_uni_dead_end_I)
        de = EnzymeRates._expand_dead_end(
            ress, uni_uni_dead_end_I;
            dead_end_regs=[:I],
            include_substrate_product=true)
        eq = EnzymeRates._expand_equivalence_constraints(
            de, uni_uni_dead_end_I)
        deduped = EnzymeRates._deduplicate(
            eq, uni_uni_dead_end_I)
        @test length(eq) == 14
        @test length(deduped) == 6
    end

    @testset "Bi-Bi random: full dedup" begin
        m_bb_random = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_B, A] ⇌ [E_A_B]
                [E, B] ⇌ [E_B]
                [E_A, B] ⇌ [E_A_B]
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_random, bi_bi)
        ress = EnzymeRates._expand_ress_variants(
            [topo], bi_bi)
        @test length(ress) == 511
        deduped = EnzymeRates._deduplicate(
            ress, bi_bi)
        @test length(deduped) == 146
    end
end

@testset "Stage 6: Allosteric expansion" begin
    @testset "No allosteric regs: passthrough" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni)
        dd = EnzymeRates._deduplicate([topo], uni_uni)
        result = EnzymeRates._expand_allosteric(
            dd, uni_uni; catalytic_n=1,
            allosteric_regs=Symbol[])
        @test length(result) == 1
    end

    @testset "1 reg, catalytic_n=1" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni_allosteric_R)
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        result = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R; catalytic_n=1,
            allosteric_regs=[:R])
        @test length(result) == 1  # m=1 only
    end

    @testset "1 reg, catalytic_n=2" begin
        # Ordered Q-first release
        m_ub1 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P_Q[AB],
                    E_Q[B], E_S[AB]
            end
            steps: begin
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_ub1, uni_bi_allosteric_R)
        dd = EnzymeRates._deduplicate(
            [topo], uni_bi_allosteric_R)
        result = EnzymeRates._expand_allosteric(
            dd, uni_bi_allosteric_R; catalytic_n=2,
            allosteric_regs=[:R])
        @test length(result) == 2  # m=1, m=2
    end

    @testset "1 reg, catalytic_n=3" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni_allosteric_R)
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        result = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R; catalytic_n=3,
            allosteric_regs=[:R])
        @test length(result) == 3  # m=1,2,3
    end

    @testset "2 regs, catalytic_n=1" begin
        # Sequential A-first bind, Q-first release
        m_bb_seq = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_seq, bi_bi_allosteric_R1_R2)
        dd = EnzymeRates._deduplicate(
            [topo], bi_bi_allosteric_R1_R2)
        result = EnzymeRates._expand_allosteric(
            dd, bi_bi_allosteric_R1_R2; catalytic_n=1,
            allosteric_regs=[:R1, :R2])
        # 2 set partitions: {R1},{R2} and {R1,R2}
        # Each with m=1 only → 2 variants
        @test length(result) == 2
    end

    @testset "2 regs, catalytic_n=2" begin
        # Sequential A-first bind, Q-first release
        m_bb_seq = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_seq, bi_bi_allosteric_R1_R2)
        dd = EnzymeRates._deduplicate(
            [topo], bi_bi_allosteric_R1_R2)
        result = EnzymeRates._expand_allosteric(
            dd, bi_bi_allosteric_R1_R2; catalytic_n=2,
            allosteric_regs=[:R1, :R2])
        # Separate sites: 2 sites × 2 mults each = 4
        # Same site: 2 mults = 2. Total = 6
        @test length(result) == 6
    end

    @testset "Dead-end edges: passthrough" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni_dead_end_I)
        de = EnzymeRates._expand_dead_end(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I],
            include_substrate_product=false)
        dd = EnzymeRates._deduplicate(
            de, uni_uni_dead_end_I)
        result = EnzymeRates._expand_allosteric(
            dd, uni_uni_dead_end_I; catalytic_n=1,
            allosteric_regs=Symbol[])
        @test length(result) == length(dd)
    end

    @testset "Dead-end I + allosteric R" begin
        # Sequential A-first bind, Q-first release
        m_bb_seq = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_seq, bi_bi_dead_end_I_allosteric_R)
        de = EnzymeRates._expand_dead_end(
            [topo], bi_bi_dead_end_I_allosteric_R;
            dead_end_regs=[:I],
            include_substrate_product=false)
        dd = EnzymeRates._deduplicate(
            de, bi_bi_dead_end_I_allosteric_R)
        result = EnzymeRates._expand_allosteric(
            dd, bi_bi_dead_end_I_allosteric_R;
            catalytic_n=1, allosteric_regs=[:R])
        # Each deduplicated spec produces 1 allosteric
        # variant (m=1) — dead-end edges pass through
        @test length(result) == length(dd)
    end

    @testset "Properties" begin
        # Ordered Q-first release
        m_ub1 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P_Q[AB],
                    E_Q[B], E_S[AB]
            end
            steps: begin
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_ub1, uni_bi_allosteric_R)
        dd = EnzymeRates._deduplicate(
            [topo], uni_bi_allosteric_R)
        result = EnzymeRates._expand_allosteric(
            dd, uni_bi_allosteric_R;
            catalytic_n=2, allosteric_regs=[:R])
        for s in result
            @test s.catalytic_n == 2
            @test !isempty(s.allosteric_reg_sites)
        end
    end
end

@testset "Stage 7: TR equivalence" begin
    @testset "Uni-Uni + R: 2^4 = 16 variants" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni_allosteric_R)
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        allo = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R; catalytic_n=1,
            allosteric_regs=[:R])
        tr = EnzymeRates._expand_tr_equivalence(
            allo, uni_uni_allosteric_R)
        # 3 metabolites (4 modes each) + 1 SS step (3 modes) = 4^3 * 3^1
        @test length(tr) == 192
    end

    @testset "Uni-Bi + R: 2^5 = 32 variants" begin
        # Ordered Q-first release
        m_ub1 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P_Q[AB],
                    E_Q[B], E_S[AB]
            end
            steps: begin
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_ub1, uni_bi_allosteric_R)
        dd = EnzymeRates._deduplicate(
            [topo], uni_bi_allosteric_R)
        allo = EnzymeRates._expand_allosteric(
            dd, uni_bi_allosteric_R; catalytic_n=1,
            allosteric_regs=[:R])
        tr = EnzymeRates._expand_tr_equivalence(
            allo, uni_bi_allosteric_R)
        # 4 metabolites (4 modes each) + 1 SS step (3 modes) = 4^4 * 3^1
        @test length(tr) == 768
    end

    @testset "Bi-Bi + R1, R2: 2^7 = 128 per spec" begin
        # Sequential A-first bind, Q-first release
        m_bb_seq = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_seq, bi_bi_allosteric_R1_R2)
        dd = EnzymeRates._deduplicate(
            [topo], bi_bi_allosteric_R1_R2)
        allo = EnzymeRates._expand_allosteric(
            dd, bi_bi_allosteric_R1_R2;
            catalytic_n=1, allosteric_regs=[:R1, :R2])
        # Use first allosteric spec only
        tr = EnzymeRates._expand_tr_equivalence(
            [allo[1]], bi_bi_allosteric_R1_R2)
        # 6 metabolites (4 modes each) + 1 SS step (3 modes) = 4^6 * 3^1
        @test length(tr) == 12288
    end

    @testset "Properties" begin
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni_allosteric_R)
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        allo = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R;
            catalytic_n=1, allosteric_regs=[:R])
        tr = EnzymeRates._expand_tr_equivalence(
            allo, uni_uni_allosteric_R)

        no_tr = filter(
            s -> isempty(s.tr_equiv_metabolites), tr)
        with_tr = filter(
            s -> !isempty(s.tr_equiv_metabolites), tr)
        @test !isempty(no_tr)
        @test !isempty(with_tr)

        # TR equiv reduces parameter count
        m_no = compile_mechanism(no_tr[1])
        m_with = compile_mechanism(with_tr[end])
        @test length(parameters(m_with)) <
            length(parameters(m_no))
    end
end

@testset "Stage 8: Allosteric deduplication" begin
    @testset "T/R mirrors dedup" begin
        # Use Bi-Bi + R1, R2 (6 metabolites) — with even
        # count, complementary subsets of equal size exist
        # (true T↔R mirrors with same param count).
        m_bb_seq = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_seq, bi_bi_allosteric_R1_R2)
        dd = EnzymeRates._deduplicate(
            [topo], bi_bi_allosteric_R1_R2)
        allo = EnzymeRates._expand_allosteric(
            dd, bi_bi_allosteric_R1_R2;
            catalytic_n=1, allosteric_regs=[:R1, :R2])
        # Use first allosteric spec
        tr = EnzymeRates._expand_tr_equivalence(
            [allo[1]], bi_bi_allosteric_R1_R2)
        @test length(tr) == 12288  # 4^6 * 3^1 (6 mets + 1 SS step)
        deduped = EnzymeRates._deduplicate_allosteric(
            tr, bi_bi_allosteric_R1_R2)
        @test length(deduped) < length(tr)
    end

    @testset "Uni-Uni + R: mirrors removed" begin
        # 3 metabolites (4 modes each) + 1 SS step (3 modes):
        # 4^3 * 3^1 = 192 variants.
        # Mirror swaps r_only↔t_only (mets) and r_only↔both (steps).
        # By Burnside: orbits = (4^n + 2^n)/2 per item type.
        # Mets (n=3): (64+8)/2 = 36. Steps (n=1): (3+1)/2 = 2.
        # Total after dedup: 36 * 2 = 72.
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni_allosteric_R)
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        allo = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R;
            catalytic_n=1, allosteric_regs=[:R])
        tr = EnzymeRates._expand_tr_equivalence(
            allo, uni_uni_allosteric_R)
        deduped = EnzymeRates._deduplicate_allosteric(
            tr, uni_uni_allosteric_R)
        @test length(deduped) == 72
    end

    @testset "Keeps lower param_count on mirror" begin
        # Sequential A-first bind, Q-first release
        m_bb_seq = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C],
                    E_A_B[CN], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        topo = mechanism_spec_from_mechanism(
            m_bb_seq, bi_bi_allosteric_R1_R2)
        dd = EnzymeRates._deduplicate(
            [topo], bi_bi_allosteric_R1_R2)
        allo = EnzymeRates._expand_allosteric(
            dd, bi_bi_allosteric_R1_R2;
            catalytic_n=1, allosteric_regs=[:R1, :R2])
        tr = EnzymeRates._expand_tr_equivalence(
            [allo[1]], bi_bi_allosteric_R1_R2)
        deduped = EnzymeRates._deduplicate_allosteric(
            tr, bi_bi_allosteric_R1_R2)
        @test length(deduped) < length(tr)
        for s in deduped
            @test s isa EnzymeRates.AllostericMechanismSpec
        end
    end

    @testset "Different base mechanisms survive" begin
        # Ordered Q-first release: E→ES→EPQ→EQ→E
        m_ub1 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P_Q[AB],
                    E_Q[B], E_S[AB]
            end
            steps: begin
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end
        # Random release: E→ES→EPQ→{EP,EQ}→E
        m_ub2 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P[A], E_P_Q[AB],
                    E_Q[B], E_S[AB]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end
        t1 = mechanism_spec_from_mechanism(
            m_ub1, uni_bi_allosteric_R)
        t2 = mechanism_spec_from_mechanism(
            m_ub2, uni_bi_allosteric_R)
        dd = EnzymeRates._deduplicate(
            [t1, t2], uni_bi_allosteric_R)
        @test length(dd) >= 2
        allo = EnzymeRates._expand_allosteric(
            dd, uni_bi_allosteric_R;
            catalytic_n=1, allosteric_regs=[:R])
        tr = EnzymeRates._expand_tr_equivalence(
            allo, uni_bi_allosteric_R)
        deduped = EnzymeRates._deduplicate_allosteric(
            tr, uni_bi_allosteric_R)
        @test length(deduped) >= 2
    end
end

@testset "Cross-stage properties" begin
    @testset "RE partition bounds" begin
        # Ordered Q-first release
        m_ub1 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P_Q[AB],
                    E_Q[B], E_S[AB]
            end
            steps: begin
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end
        spec = mechanism_spec_from_mechanism(
            m_ub1, uni_bi)
        ress = EnzymeRates._expand_ress_variants(
            [spec], uni_bi)
        for s in ress
            partition =
                EnzymeRates._compute_re_partition_from_steps(
                    s.steps)
            @test 1 <= length(partition) <= 7
        end
    end

    @testset "Stage monotonicity" begin
        for rxn in [uni_bi_reg_unknown]
            counts = _run_full_pipeline_stages(rxn)
            @test counts.dead_end >= counts.ress
            @test counts.equivalence >= counts.dead_end
            @test counts.dedup <= counts.equivalence
        end
    end

    @testset "Regulator roles affect partitioning" begin
        c_unk = _run_full_pipeline_stages(
            uni_uni_reg_unknown)
        c_de = _run_full_pipeline_stages(
            uni_uni_dead_end_I)
        c_al = _run_full_pipeline_stages(
            uni_uni_allosteric_R)
        @test c_unk.catalytic == c_de.catalytic
        @test c_unk.catalytic == c_al.catalytic
        @test c_de.dedup >= c_al.dedup
    end

    @testset "param_count accuracy per stage" begin
        rng = Random.MersenneTwister(99)

        # Stage 1: catalytic topologies
        catalytic = EnzymeRates._catalytic_topologies(
            uni_bi)
        for s in catalytic
            m = compile_mechanism(s)
            @test s.param_count ==
                length(parameters(m))
        end

        # Stage 2: RE/SS expansion
        ress = EnzymeRates._expand_ress_variants(
            catalytic, uni_bi)
        for s in ress[1:min(10, length(ress))]
            m = compile_mechanism(s)
            @test s.param_count ==
                length(parameters(m))
        end

        # Stage 4: equivalence constraints
        eq =
            EnzymeRates._expand_equivalence_constraints(
                ress, uni_bi)
        for s in eq[1:min(10, length(eq))]
            m = compile_mechanism(s)
            @test s.param_count ==
                length(parameters(m))
        end

        # With regulators: dead-end + equivalence
        cat_r = EnzymeRates._catalytic_topologies(
            uni_bi_dead_end_I)
        ress_r = EnzymeRates._expand_ress_variants(
            cat_r, uni_bi_dead_end_I)
        de_r =
            EnzymeRates._expand_dead_end(
                ress_r, uni_bi_dead_end_I;
                dead_end_regs=[:I],
                include_substrate_product=true)
        sample_de = de_r[randperm(
            rng, length(de_r))[
            1:min(10, length(de_r))]]
        for s in sample_de
            m = compile_mechanism(s)
            @test s.param_count ==
                length(parameters(m))
        end
    end

    @testset "compile_mechanism round-trip" begin
        for rxn in [uni_uni, uni_bi]
            all_specs = collect(
                EnzymeRates.old_enumerate_mechanisms(rxn))
            cat_specs = filter(
                s -> s isa EnzymeRates.MechanismSpec, all_specs)
            for s in first(cat_specs, 3)
                m = compile_mechanism(s)
                @test m isa EnzymeMechanism
                s2 = mechanism_spec_from_mechanism(
                    m, rxn)
                @test s2.steps == s.steps
                @test s2.param_count >= s.param_count
            end
        end

        @testset "allosteric compilation" begin
            all_specs = collect(
                EnzymeRates.old_enumerate_mechanisms(
                    uni_uni_allosteric_R))
            allo_specs = filter(
                s -> s isa EnzymeRates.AllostericMechanismSpec,
                all_specs)
            for s in first(allo_specs, 2)
                m = compile_mechanism(s)
                @test m isa AllostericEnzymeMechanism
                @test length(metabolites(m)) > 0
                @test length(parameters(m)) > 0
            end
        end
    end
end

@testset "End-to-end pipeline" begin
    @testset "Uni-Uni, no regs" begin
        result = collect(
            EnzymeRates.old_enumerate_mechanisms(uni_uni))
        @test length(result) == 3
    end

    @testset "Uni-Bi, no regs" begin
        result = collect(
            EnzymeRates.old_enumerate_mechanisms(uni_bi))
        @test length(result) == 56
    end

    @testset "Bi-Bi, no regs" begin
        stats = @timed collect(
            EnzymeRates.old_enumerate_mechanisms(bi_bi))
        @test length(stats.value) == 63762
        # Performance: ~3s / 5GB baseline
        @test stats.bytes < 15 * 1024^3
        @test stats.time < 30
    end
    @testset "Bi-Bi Ping-Pong, no regs" begin
        stats = @timed collect(
            EnzymeRates.old_enumerate_mechanisms(
                bi_bi_ping_pong))
        @test length(stats.value) == 64276
        # Performance: ~7s / 6GB baseline
        @test stats.bytes < 15 * 1024^3
        @test stats.time < 30
    end

    @testset "Uni-Uni + 1 unknown reg" begin
        result = collect(
            EnzymeRates.old_enumerate_mechanisms(
                uni_uni_reg_unknown))
        @test length(result) == 101
    end

    @testset "Uni-Bi + 1 unknown reg" begin
        result = collect(
            EnzymeRates.old_enumerate_mechanisms(
                uni_bi_reg_unknown))
        @test length(result) == 3794
    end
end

@testset "param_count accuracy" begin
    @testset "All Uni-Bi specs" begin
        all_specs = collect(
            EnzymeRates.old_enumerate_mechanisms(uni_bi))
        @test length(all_specs) == 56
        n_match = count(all_specs) do s
            m = compile_mechanism(s)
            s.param_count == length(parameters(m))
        end
        @test n_match == 56
    end

    # AllostericMechanismSpec.base.param_count only covers
    # catalytic parameters, not total allosteric params.
    # Accuracy check is limited to verifying compilation
    # succeeds and produces parameters.
    @testset "Allosteric specs (compilation)" begin
        all_specs = collect(
            EnzymeRates.old_enumerate_mechanisms(
                uni_uni_allosteric_R))
        allo_specs = filter(
            s -> s isa EnzymeRates.AllostericMechanismSpec,
            all_specs)
        @test !isempty(allo_specs)
        for s in first(allo_specs, 3)
            m = compile_mechanism(s)
            @test length(parameters(m)) > 0
        end
    end

    @testset "Sampled Bi-Bi specs" begin
        rng_bb = Random.MersenneTwister(42)
        all_specs = collect(
            EnzymeRates.old_enumerate_mechanisms(bi_bi))
        @test length(all_specs) == 63762
        sample = all_specs[randperm(
            rng_bb, length(all_specs))[1:min(50, end)]]
        n_match = count(sample) do s
            s isa EnzymeRates.MechanismSpec || return true
            m = compile_mechanism(s)
            s.param_count == length(parameters(m))
        end
        @test n_match == length(sample)
    end
end

end # outer testset
