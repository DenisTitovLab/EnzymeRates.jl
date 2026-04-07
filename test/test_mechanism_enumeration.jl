# ABOUTME: Tests for mechanism enumeration pipeline
# ABOUTME: Unit tests per move + integration tests per reaction

using EnzymeRates: StepSpec, MechanismSpec, AllostericMechanismSpec,
    ParamConstraint, AbstractMechanismSpec

# Helper: convert EnzymeMechanism → MechanismSpec
function mechanism_spec_from_mechanism(
    m::EnzymeMechanism,
    @nospecialize(rxn::EnzymeReaction))
    rxns = EnzymeRates.reactions(m)
    eq_steps_tuple = EnzymeRates.equilibrium_steps(m)
    pc = EnzymeRates.param_constraints(m)
    steps = StepSpec[]
    for (i, (lhs, rhs)) in enumerate(rxns)
        push!(steps, StepSpec(
            collect(lhs), collect(rhs),
            eq_steps_tuple[i]))
    end
    constraints = [
        (t, c, [(s, sc) for (s, sc) in f])
        for (t, c, f) in pc
    ]
    MechanismSpec(
        rxn, steps, constraints,
        length(parameters(m)))
end

const uni_uni_rxn = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end

const uni_bi_rxn = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
end

const bi_bi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

const bi_bi_pp_rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end

const uni_uni_with_reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    dead_end_inhibitors: I
end

const uni_uni_allo = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    oligomeric_state: 2
end

const uni_uni_allo_reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R
    oligomeric_state: 2
end

const uni_uni_allo_2reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R1, R2
    oligomeric_state: 2
end

const ter_ter_rxn = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
end

const ter_bi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[CN], Q[X]
end

# Pyruvate carboxylase: Pyr + HCO3 + ATP = OAA + ADP + Pi
# Mechanism: ATP+HCO3 → ADP+Pi+CO2_residual,
#            then Pyr+CO2 → OAA
const pyruvate_carboxylase_rxn = @enzyme_reaction begin
    substrates: Pyr[C3H3O3], HCO3[HCO3], ATP[C10H16N5O13P3]
    products: OAA[C4H3O5], ADP[C10H15N5O10P2], Pi[H2PO4]
end

# Pyruvate dehydrogenase: Pyr + NAD + CoA = AcCoA + NADH + CO2
# Mechanism: Pyr → CO2+residual, CoA+residual → AcCoA+residual,
#            NAD+residual → NADH
const pyruvate_dehydrogenase_rxn = @enzyme_reaction begin
    substrates: Pyr[C3H3O3], NAD[C21H28N7O14P2], CoA[C21H36N7O16P3S]
    products: AcCoA[C23H38N7O17P3S], NADH[C21H29N7O14P2], CO2[CO2]
end

"""Collect all mechanisms by running the full enumeration loop."""
function enumerate_all(
    @nospecialize(reaction::EnzymeReaction);
    max_params::Int=20)
    cache = Dict{Int, Vector{EnzymeRates.AbstractMechanismSpec}}()

    init_specs = EnzymeRates.init_mechanisms(reaction)
    min_pc = init_specs[1].param_count
    cache[min_pc] = EnzymeRates.AbstractMechanismSpec[init_specs...]
    EnzymeRates.dedup!(cache)

    results = Dict{Int, Vector{EnzymeRates.AbstractMechanismSpec}}()

    for pc in min_pc:max_params
        level = pop!(cache, pc, EnzymeRates.AbstractMechanismSpec[])
        isempty(level) && (isempty(cache) ? break : continue)

        results[pc] = level

        new_specs = EnzymeRates.expand_mechanisms(level, reaction)
        for (target_pc, specs) in new_specs
            target_pc > max_params && continue
            append!(get!(cache, target_pc,
                EnzymeRates.AbstractMechanismSpec[]), specs)
        end
        EnzymeRates.dedup!(cache)
    end
    results
end

@testset "Mechanism Enumeration" begin

@testset "test reaction atom balance" begin
    for rxn in [pyruvate_carboxylase_rxn,
                pyruvate_dehydrogenase_rxn]
        sub_atoms = Dict{Symbol,Int}()
        for (_, atoms) in EnzymeRates.substrates(rxn)
            for (a, c) in atoms
                sub_atoms[a] = get(sub_atoms, a, 0) + c
            end
        end
        prod_atoms = Dict{Symbol,Int}()
        for (_, atoms) in EnzymeRates.products(rxn)
            for (a, c) in atoms
                prod_atoms[a] = get(prod_atoms, a, 0) + c
            end
        end
        @test sub_atoms == prod_atoms
    end
end

@testset "AllostericEnzymeMechanism TR equivalence" begin
    base_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    base_steps = [
        StepSpec([:E, :S], [:E_S], true),
        StepSpec([:E, :P], [:E_P], true),
        StepSpec([:E_S], [:E_P], false),
    ]
    base_spec = MechanismSpec(
        base_rxn, base_steps, ParamConstraint[], 3)

    # S is TR-equivalent (K_T_S = K_R_S),
    # P is R-only (absent from T-state)
    allo_spec = AllostericMechanismSpec(
        base_spec, 2,
        Vector{Symbol}[], Int[],
        [:S], Int[],
        Symbol[], [:P], Int[],
        base_spec.param_count + 1)

    m_compiled = AllostericEnzymeMechanism(allo_spec)
    cat_sites = typeof(m_compiled).parameters[3]
    @test :S in cat_sites[3]   # tr_equiv_mets
    @test :P in cat_sites[6]   # t_only_mets
end

@testset "Catalytic topologies" begin

    @testset "Uni-Uni" begin
        topos = EnzymeRates._catalytic_topologies(uni_uni_rxn)
        @test length(topos) == 1

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
        spec_rt = mechanism_spec_from_mechanism(
            m_uu, uni_uni_rxn)
        @test EnzymeMechanism(spec_rt) === m_uu
        @test EnzymeMechanism(topos[1]) === m_uu
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps) == 1
        end
    end

    @testset "Uni-Bi" begin
        topos = EnzymeRates._catalytic_topologies(uni_bi_rxn)
        @test length(topos) == 3
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps) == 1
            m = EnzymeMechanism(t)
            @test m isa EnzymeMechanism
        end
    end

    @testset "Bi-Bi" begin
        topos = EnzymeRates._catalytic_topologies(bi_bi_rxn)
        @test length(topos) == 9
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps) == 1
        end
    end

    @testset "Bi-Bi Ping-Pong" begin
        topos = EnzymeRates._catalytic_topologies(
            bi_bi_pp_rxn)
        @test length(topos) == 10
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps) == 1
        end
    end

    @testset "Ter-Ter" begin
        topos = EnzymeRates._catalytic_topologies(
            ter_ter_rxn)
        # 169 = 13 × 13. Weak orderings of 3 items = 13
        # (Fubini number F(3)). Substrate and product
        # orderings are independent.
        @test length(topos) == 169
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps) == 1
        end
    end

    @testset "Ter-Bi" begin
        topos = EnzymeRates._catalytic_topologies(
            ter_bi_rxn)
        # 45 = 39 sequential + 6 ping-pong
        # Sequential: F(3) × F(2) = 13 × 3 = 39
        # Ping-pong: 6 topologies from D[X]→Q[X] iso group
        @test length(topos) == 45
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps) == 1
        end
    end

    @testset "weak-ordering combining" begin
        # For bi-bi: 2 subs × 2 prods, all sequential
        # Weak orderings of 2 items = 3 (F(2))
        # Total: 3 × 3 = 9 topologies
        bi_bi_rxn_test = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
        end
        topos = EnzymeRates._catalytic_topologies(
            bi_bi_rxn_test)
        @test length(topos) == 9

        # For ter-ter: 3 subs × 3 prods, all sequential
        # Weak orderings of 3 items = 13 (F(3))
        # Total: 13 × 13 = 169 topologies
        topos_tt = EnzymeRates._catalytic_topologies(
            ter_ter_rxn)
        @test length(topos_tt) == 169
    end
end

@testset "init_mechanisms" begin

    @testset "Param count invariant" begin
        for (rxn, n_s, n_p) in [
            (uni_uni_rxn, 1, 1),
            (uni_bi_rxn, 1, 2),
            (bi_bi_rxn, 2, 2),
            (bi_bi_pp_rxn, 2, 2),
        ]
            specs = EnzymeRates.init_mechanisms(rxn)
            expected_pc = n_s + n_p + 3
            for s in specs
                @test s.param_count == expected_pc
            end
        end
    end

    @testset "All have exactly 1 SS step" begin
        for rxn in [uni_uni_rxn, uni_bi_rxn,
                    bi_bi_rxn, bi_bi_pp_rxn]
            specs = EnzymeRates.init_mechanisms(rxn)
            for s in specs
                @test count(
                    !st.is_equilibrium
                    for st in s.steps) == 1
            end
        end
    end

    @testset "Uni-Uni: no dead-end forms" begin
        specs = EnzymeRates.init_mechanisms(
            uni_uni_rxn)
        @test length(specs) == 1
    end

    @testset "Dead-end substrate/product expansion" begin

        @testset "Uni-Uni: no dead-end forms" begin
            # 3 forms: E, E_S[C], E_P[C]. E_S has all subs,
            # E_P has all prods. No mixed dead-end possible.
            # → 0 dead-end forms, 1 variant (bare topology)
            m = @enzyme_mechanism begin
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
            spec = mechanism_spec_from_mechanism(m, uni_uni_rxn)
            @test EnzymeMechanism(spec) === m
            result =
                EnzymeRates._expand_substrate_product_dead_ends(
                    [spec], uni_uni_rxn)
            @test length(result) == 1
        end

        @testset "Bi-Bi random: 4 dead-end forms" begin
            # 7 forms: E, E_A, E_B, E_A_B, E_P, E_Q, E_P_Q
            # Eligible dead-end forms (mixed sub+prod binding):
            #   E_A: +P→E_A_P(mixed✓), +Q→E_A_Q(mixed✓)
            #   E_B: +P→E_B_P(mixed✓), +Q→E_B_Q(mixed✓)
            #   E_P: +A→E_A_P(same), +B→E_B_P(same)
            #   E_Q: +A→E_A_Q(same), +B→E_B_Q(same)
            # 4 unique dead-end forms → 2^4 = 16 variants
            m = @enzyme_mechanism begin
                species: begin
                    substrates: A[C], B[N]
                    products: P[C], Q[N]
                    enzymes: E, E_A[C], E_A_B[CN],
                        E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
            spec = mechanism_spec_from_mechanism(m, bi_bi_rxn)
            @test EnzymeMechanism(spec) === m
            result =
                EnzymeRates._expand_substrate_product_dead_ends(
                    [spec], bi_bi_rxn)
            # 4 unique dead-end forms, 7 competition patterns,
            # all 7 produce distinct dead-end sets → 7 variants
            @test length(result) == 7
        end

        @testset "Uni-Bi ordered: no dead-end forms" begin
            # 4 forms: E, E_S, E_P_Q, E_Q
            # E+P→E_P: single-product → rejected (need mixed)
            # E_Q+S→E_S_Q: has all subs → rejected
            # → 0 dead-end forms, 1 variant
            m = @enzyme_mechanism begin
                species: begin
                    substrates: S[AB]
                    products: P[A], Q[B]
                    enzymes: E, E_P_Q[AB], E_Q[B], E_S[AB]
                end
                steps: begin
                    [E, Q] ⇌ [E_Q]
                    [E_Q, P] ⇌ [E_P_Q]
                    [E, S] ⇌ [E_S]
                    [E_S] <--> [E_P_Q]
                end
            end
            spec = mechanism_spec_from_mechanism(m, uni_bi_rxn)
            @test EnzymeMechanism(spec) === m
            result =
                EnzymeRates._expand_substrate_product_dead_ends(
                    [spec], uni_bi_rxn)
            @test length(result) == 1
        end

        @testset "Bi-Bi Ping-Pong: 3 dead-end forms" begin
            # Forms: E, E_A, Estar, Estar_A_P, Estar_B, E_Q
            # E_A: +P→E_A_P(mixed✓), +Q→E_A_Q(mixed✓)
            # E_Q: +B→E_B_Q(mixed✓)
            # 3 dead-end forms → 2^3 = 8 variants
            m = @enzyme_mechanism begin
                species: begin
                    substrates: A[CX], B[N]
                    products: P[C], Q[NX]
                    enzymes: E, E_A[CX], E_Q[NX],
                        Estar[X], Estar_A_P[CX], Estar_B[NX]
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
            spec = mechanism_spec_from_mechanism(
                m, bi_bi_pp_rxn)
            @test EnzymeMechanism(spec) === m
            result =
                EnzymeRates._expand_substrate_product_dead_ends(
                    [spec], bi_bi_pp_rxn)
            # 3 dead-end forms, 7 competition patterns,
            # 5 unique dead-end sets after dedup → 5 variants
            @test length(result) == 5
        end
    end

    @testset "Competition patterns" begin
        # Uni-uni: 1×1, only 1 pattern (single edge)
        pats_11 = EnzymeRates._competition_patterns(
            Set([:S]), Set([:P]))
        @test length(pats_11) == 1
        @test pats_11[1] == Set([(:S, :P)])

        # Uni-bi: 1×2, S competes with both P and Q
        pats_12 = EnzymeRates._competition_patterns(
            Set([:S]), Set([:P, :Q]))
        @test length(pats_12) == 1
        @test pats_12[1] == Set([(:S, :P), (:S, :Q)])

        # Bi-uni: symmetric
        pats_21 = EnzymeRates._competition_patterns(
            Set([:A, :B]), Set([:P]))
        @test length(pats_21) == 1
        @test pats_21[1] == Set([(:A, :P), (:B, :P)])

        # Bi-bi: 7 patterns
        pats_22 = EnzymeRates._competition_patterns(
            Set([:A, :B]), Set([:P, :Q]))
        @test length(pats_22) == 7
        # Every pattern covers all vertices
        for pat in pats_22
            for s in [:A, :B]
                @test any(
                    p -> (s, p) in pat, [:P, :Q])
            end
            for p in [:P, :Q]
                @test any(
                    s -> (s, p) in pat, [:A, :B])
            end
        end
        # Invalid: {A↔P, B↔P} leaves Q uncovered
        @test Set([(:A, :P), (:B, :P)]) ∉ pats_22

        # Ter-ter: 265 patterns
        pats_33 = EnzymeRates._competition_patterns(
            Set([:A, :B, :C]),
            Set([:P, :Q, :R]))
        @test length(pats_33) == 265
        for pat in pats_33
            for s in [:A, :B, :C]
                @test any(
                    p -> (s, p) in pat,
                    [:P, :Q, :R])
            end
            for p in [:P, :Q, :R]
                @test any(
                    s -> (s, p) in pat,
                    [:A, :B, :C])
            end
        end
    end

    @testset "Dead-end filtering by competition" begin

        # Shared bi-bi random mechanism for multiple tests
        m_bb = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        spec_bb = mechanism_spec_from_mechanism(
            m_bb, bi_bi_rxn)

        @testset "Bi-bi random: 7 variants (was 16)" begin
            # 4 dead-end forms × 7 competition patterns.
            # Each pattern yields a distinct dead-end set:
            #   {A↔P,B↔Q}: {E_A_Q,E_B_P}
            #   {A↔Q,B↔P}: {E_A_P,E_B_Q}
            #   {A↔P,A↔Q,B↔P}: {E_B_Q}
            #   {A↔P,A↔Q,B↔Q}: {E_B_P}
            #   {A↔P,B↔P,B↔Q}: {E_A_Q}
            #   {A↔Q,B↔P,B↔Q}: {E_A_P}
            #   {A↔P,A↔Q,B↔P,B↔Q}: {} (no dead-ends)
            # All 7 sets are distinct → 7 variants after dedup
            result =
                EnzymeRates._expand_substrate_product_dead_ends(
                    [spec_bb], bi_bi_rxn)
            @test length(result) == 7
        end

        @testset "Bi-bi random: complete competition → bare topology" begin
            result =
                EnzymeRates._expand_substrate_product_dead_ends(
                    [spec_bb], bi_bi_rxn)
            # Complete pattern {A↔P,A↔Q,B↔P,B↔Q} forbids
            # all dead-end forms → 1 variant has no dead-end
            # steps (same step count as original)
            bare = filter(
                r -> length(r.steps) == length(spec_bb.steps),
                result)
            @test length(bare) == 1
        end

        @testset "Bi-bi random: diagonal has exactly 2 dead-end forms" begin
            result =
                EnzymeRates._expand_substrate_product_dead_ends(
                    [spec_bb], bi_bi_rxn)
            # Diagonal patterns {A↔P,B↔Q} and {A↔Q,B↔P}
            # each allow exactly 2 dead-end forms.
            two_de = filter(result) do r
                de_forms = setdiff(
                    EnzymeRates.all_form_names(r),
                    EnzymeRates.all_form_names(spec_bb))
                length(de_forms) == 2
            end
            @test length(two_de) == 2  # diagonal + anti-diagonal
        end

        @testset "Ter-ter per-topology (OOM on full init)" begin
            # Test that competition filtering works
            # on representative ter-ter topologies.
            topos = EnzymeRates._catalytic_topologies(
                ter_ter_rxn)
            @test length(topos) == 169
            # Test first (random, most forms) and last topology
            for topo in [topos[1], topos[end]]
                result =
                    EnzymeRates._expand_substrate_product_dead_ends(
                        [topo], ter_ter_rxn)
                # Competition patterns reduce 2^27 to
                # ≤265 variants per topology
                @test length(result) > 0
                @test length(result) <= 265
                for spec in result
                    # Verify param_count set correctly
                    @test spec.param_count >=
                        length(EnzymeRates.substrates(
                            ter_ter_rxn)) +
                        length(EnzymeRates.products(
                            ter_ter_rxn)) + 3
                end
            end
        end

        @testset "Ter-ter diagonal: 12 of 27 allowed" begin
            # Random ter-ter topology has 27 possible
            # dead-end forms. With diagonal competition
            # {A↔P, B↔Q, D↔R}:
            #   1S+1P: 6 allowed, 3 forbidden (the 3
            #     diagonal pairs A-P, B-Q, D-R)
            #   2S+1P: 3 allowed (E_A_B_R, E_A_D_Q,
            #     E_B_D_P), 6 forbidden
            #   1S+2P: 3 allowed (E_A_Q_R, E_B_P_R,
            #     E_D_P_Q), 6 forbidden
            #   Total: 12 allowed out of 27
            topos = EnzymeRates._catalytic_topologies(
                ter_ter_rxn)
            random_topo = topos[1]  # most forms
            bound =
                EnzymeRates._bound_metabolites_at_forms(
                    random_topo, ter_ter_rxn)
            sub_names = Set([:A, :B, :D])
            prod_names = Set([:P, :Q, :R])
            cat_forms =
                EnzymeRates.all_form_names(random_topo)
            _sp_de_opps =
                EnzymeRates._substrate_product_dead_end_opportunities
            de_opps = _sp_de_opps(
                bound, cat_forms,
                sub_names, prod_names)
            # Group dead-end forms
            de_forms = Dict{Symbol,
                Vector{Tuple{Symbol, Symbol}}}()
            for (f, m) in de_opps
                de_name =
                    EnzymeRates._dead_end_form_name(
                        bound[f], m)
                push!(get!(de_forms, de_name,
                    Tuple{Symbol, Symbol}[]), (f, m))
            end
            de_form_names =
                sort(collect(keys(de_forms)))
            @test length(de_form_names) == 27

            # Build de_bound mapping
            de_bound = Dict{Symbol, Set{Symbol}}()
            for de_name in de_form_names
                f, m = first(de_forms[de_name])
                de_bound[de_name] = union(
                    bound[f], Set([m]))
            end

            # Apply diagonal competition filter
            diagonal =
                Set([(:A, :P), (:B, :Q), (:D, :R)])
            allowed = Symbol[]
            for de_name in de_form_names
                mets = de_bound[de_name]
                de_subs = intersect(mets, sub_names)
                de_prods =
                    intersect(mets, prod_names)
                has_conflict = any(
                    (s, p) in diagonal
                    for s in de_subs
                    for p in de_prods)
                has_conflict ||
                    push!(allowed, de_name)
            end
            @test length(allowed) == 12

            # Verify specific allowed forms
            @test :E_A_Q in allowed   # 1S+1P
            @test :E_A_R in allowed
            @test :E_B_P in allowed
            @test :E_B_R in allowed
            @test :E_D_P in allowed
            @test :E_D_Q in allowed
            @test :E_A_B_R in allowed # 2S+1P
            @test :E_A_D_Q in allowed
            @test :E_B_D_P in allowed
            @test :E_A_Q_R in allowed # 1S+2P
            @test :E_B_P_R in allowed
            @test :E_D_P_Q in allowed

            # Verify specific forbidden forms
            @test :E_A_P ∉ allowed    # A↔P diagonal
            @test :E_B_Q ∉ allowed    # B↔Q diagonal
            @test :E_D_R ∉ allowed    # D↔R diagonal
        end

        @testset "Round-trip: competition-filtered specs compile" begin
            for rxn in [uni_uni_rxn, bi_bi_rxn, bi_bi_pp_rxn]
                specs =
                    EnzymeRates.init_mechanisms(rxn)
                # Test first 5 specs (compilation can be slow)
                for spec in first(specs, 5)
                    m = EnzymeMechanism(spec)
                    @test m isa EnzymeMechanism
                    @test length(parameters(m)) <=
                        spec.param_count
                end
            end
        end
    end
end

@testset "RE→SS conversion" begin
    @testset "Multiple RE steps" begin
        m = @enzyme_mechanism begin
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
        spec = mechanism_spec_from_mechanism(m, uni_uni_rxn)
        @test EnzymeMechanism(spec) === m
        result = EnzymeRates._expand_re_to_ss(spec)
        @test length(result) == 2
        for r in result
            @test r.param_count == spec.param_count + 1
        end
    end

    @testset "All SS → yields nothing" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] <--> [E_P]
                [E, S] <--> [E_S]
                [E_S] <--> [E_P]
            end
        end
        spec = mechanism_spec_from_mechanism(m, uni_uni_rxn)
        @test EnzymeMechanism(spec) === m
        result = EnzymeRates._expand_re_to_ss(spec)
        @test isempty(result)
    end

    @testset "All-SS with constrained dead-end RE → nothing" begin
        # Uni-uni where all catalytic steps are SS.
        # Dead-end inhibitor I binds to 2 forms (E and E_P),
        # creating 2 RE binding steps with K constrained equal.
        # These constrained RE steps should be skipped.
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] <--> [E_P]
                [E, S] <--> [E_S]
                [E_S] <--> [E_P]
            end
        end
        spec = mechanism_spec_from_mechanism(m, uni_uni_rxn)
        @test EnzymeMechanism(spec) === m
        # Add dead-end inhibitor I binding to 2 forms
        rxn_i = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
        end
        de_specs = EnzymeRates._expand_add_dead_end_regulator(
            spec, rxn_i)
        # Find a spec where I binds to 2+ forms (creating
        # constrained RE binding steps)
        multi_form = filter(de_specs) do s
            n = count(s.steps) do st
                any(contains(string(sym), "I__reg")
                    for sym in Iterators.flatten(
                        (st.reactants, st.products)))
            end
            n >= 2
        end
        if !isempty(multi_form)
            spec_de = first(multi_form)
            # The dead-end RE binding steps should be
            # constrained → _expand_re_to_ss yields nothing
            # (all catalytic steps are already SS, and the
            # only RE steps are constrained dead-end bindings)
            result = EnzymeRates._expand_re_to_ss(spec_de)
            @test isempty(result)
        end
    end

    @testset "Bi-bi: exact RE→SS count" begin
        # Bi-bi random: 9 steps. In init_mechanisms form:
        # 1 SS step (isomerization), 8 RE binding steps.
        # With max constraints (K_A, K_B, K_P, K_Q each
        # constrained across 2 forms): 4 constraint groups,
        # 4 leaders + 4 followers → all 8 RE steps constrained.
        # After removing one constraint: 6 constrained + 2 free.
        # The 2 freed RE steps are eligible for RE→SS.
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        spec = mechanism_spec_from_mechanism(m, bi_bi_rxn)
        @test EnzymeMechanism(spec) === m
        # Add max constraints
        spec_c = EnzymeRates.MechanismSpec(
            spec.reaction, spec.steps,
            EnzymeRates._max_equivalence_constraints(spec),
            spec.param_count)
        # With all K's constrained, no RE→SS possible
        @test isempty(EnzymeRates._expand_re_to_ss(spec_c))
        # Remove one constraint to free 2 RE steps
        unconstrained = first(
            EnzymeRates._expand_remove_constraint(spec_c))
        result = EnzymeRates._expand_re_to_ss(unconstrained)
        # Should have results (freed RE steps now eligible)
        @test !isempty(result)
    end
end

@testset "Remove equivalence constraint" begin
    @testset "Multiple constraints: exact count" begin
        # Bi-bi random with max constraints:
        # K_A constrained across 2 forms, K_B same,
        # K_P same, K_Q same = 4 constraint groups.
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        spec = mechanism_spec_from_mechanism(m, bi_bi_rxn)
        @test EnzymeMechanism(spec) === m
        spec_c = EnzymeRates.MechanismSpec(
            spec.reaction, spec.steps,
            EnzymeRates._max_equivalence_constraints(spec),
            spec.param_count)
        n_constraints = length(spec_c.param_constraints)
        @test n_constraints == 4  # K_A, K_B, K_P, K_Q
        result = EnzymeRates._expand_remove_constraint(spec_c)
        # 4 RE constraint groups → 4 results, each at +1
        @test length(result) == n_constraints
        for r in result
            @test length(r.param_constraints) ==
                n_constraints - 1
            @test r.param_count == spec_c.param_count + 1
        end
    end

    @testset "No constraints → yields nothing" begin
        m = @enzyme_mechanism begin
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
        spec = mechanism_spec_from_mechanism(m, uni_uni_rxn)
        @test EnzymeMechanism(spec) === m
        @test isempty(spec.param_constraints)
        result = EnzymeRates._expand_remove_constraint(spec)
        @test isempty(result)
    end

    @testset "SS constraints removed as kf/kr pairs" begin
        # Build a spec with SS constraints: take a bi-bi
        # mechanism, convert constrained steps to SS, then
        # rebuild constraints. This produces kf/kr pairs
        # that must be removed together.
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        spec = mechanism_spec_from_mechanism(m, bi_bi_rxn)
        @test EnzymeMechanism(spec) === m
        max_c = EnzymeRates._max_equivalence_constraints(spec)
        c_idxs = EnzymeRates._constrained_step_indices(max_c)
        new_steps = [
            EnzymeRates.StepSpec(
                s.reactants, s.products,
                (i in c_idxs) ? false : s.is_equilibrium)
            for (i, s) in enumerate(spec.steps)]
        ss_spec = EnzymeRates.MechanismSpec(
            spec.reaction, new_steps,
            ParamConstraint[], spec.param_count)
        # Rebuild constraints; the SS group gets kf/kr pairs
        new_constraints =
            EnzymeRates._max_equivalence_constraints(ss_spec)
        has_ss_pair = any(
            endswith(string(c[1]), "f")
            for c in new_constraints)
        if has_ss_pair
            spec_with_ss = EnzymeRates.MechanismSpec(
                ss_spec.reaction, ss_spec.steps,
                new_constraints, spec.param_count)
            result = EnzymeRates._expand_remove_constraint(
                spec_with_ss)
            # Each result must not have an orphaned kf
            # without its matching kr (or vice versa)
            for r in result
                for c in r.param_constraints
                    s = string(c[1])
                    if endswith(s, "f")
                        kr_sym = Symbol(
                            s[1:end-1] * "r")
                        @test any(
                            c2[1] == kr_sym
                            for c2 in r.param_constraints)
                    elseif endswith(s, "r")
                        kf_sym = Symbol(
                            s[1:end-1] * "f")
                        @test any(
                            c2[1] == kf_sym
                            for c2 in r.param_constraints)
                    end
                end
            end
            # param_count delta: +1 (RE) or +2 (SS pair)
            for r in result
                delta = r.param_count -
                    spec_with_ss.param_count
                @test delta == 1 || delta == 2
            end
        end
    end
end

@testset "Inhibitor competition patterns" begin
    # Uni-uni, no existing inhibitors
    pats = EnzymeRates._inhibitor_competition_patterns(
        Set([:S]), Set([:P]), Symbol[])
    @test length(pats) == 1
    @test pats[1] == (Set([:S]), Set([:P]), Set{Symbol}())

    # Bi-bi, no existing inhibitors: 3×3 = 9
    pats_bb = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B]), Set([:P, :Q]), Symbol[])
    @test length(pats_bb) == 9

    # Ter-ter, no existing inhibitors: 7×7 = 49
    pats_tt = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B, :C]), Set([:P, :Q, :R]), Symbol[])
    @test length(pats_tt) == 49

    # Bi-bi, 1 existing inhibitor: 9 × 2 = 18
    pats_1i = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B]), Set([:P, :Q]), [:I1__reg])
    @test length(pats_1i) == 18

    # Bi-bi, 2 existing inhibitors: 9 × 4 = 36
    pats_2i = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B]), Set([:P, :Q]),
        [:I1__reg, :I2__reg])
    @test length(pats_2i) == 36
end

@testset "Forms with binding step" begin
    # Uni-uni: S binds to E, P binds to E
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
    spec_uu = mechanism_spec_from_mechanism(
        m_uu, uni_uni_rxn)
    @test EnzymeRates._forms_with_binding_step(
        spec_uu.steps, :S) == Set([:E])
    @test EnzymeRates._forms_with_binding_step(
        spec_uu.steps, :P) == Set([:E])

    # Bi-bi random: B binds to E and E_A
    m_bb = @enzyme_mechanism begin
        species: begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            enzymes: E, E_A[C], E_A_B[CN],
                E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
    spec_bb = mechanism_spec_from_mechanism(
        m_bb, bi_bi_rxn)
    @test EnzymeRates._forms_with_binding_step(
        spec_bb.steps, :B) == Set([:E, :E_A])
    @test EnzymeRates._forms_with_binding_step(
        spec_bb.steps, :A) == Set([:E, :E_B])
    @test EnzymeRates._forms_with_binding_step(
        spec_bb.steps, :P) == Set([:E, :E_Q])
    @test EnzymeRates._forms_with_binding_step(
        spec_bb.steps, :Q) == Set([:E, :E_P])
end

@testset "Add dead-end regulator" begin
    # Helper: build MechanismSpec from @enzyme_mechanism
    # with a different reaction (e.g. one with regulators)
    function _spec_with_rxn(m, rxn_base, rxn_target)
        spec = mechanism_spec_from_mechanism(
            m, rxn_base)
        @test EnzymeMechanism(spec) === m
        MechanismSpec(rxn_target, spec.steps,
            spec.param_constraints, spec.param_count)
    end

    @testset "Uni-uni + new regulator" begin
        m = @enzyme_mechanism begin
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
        spec = _spec_with_rxn(
            m, uni_uni_rxn, uni_uni_with_reg)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_with_reg)
        # Uni-uni: sub_names={S}, prod_names={P}
        # E: bound={} → eligible
        # E_S: bound={S} → all subs → NOT eligible
        # E_P: bound={P} → all prods → NOT eligible
        # 1 eligible form → 2^1 - 1 = 1 variant
        @test length(result) == 1
        for r in result
            @test r.param_count == spec.param_count + 1
        end
    end

    @testset "No regulators → yields nothing" begin
        m = @enzyme_mechanism begin
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
        spec = mechanism_spec_from_mechanism(
            m, uni_uni_rxn)
        @test EnzymeMechanism(spec) === m
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_rxn)
        @test isempty(result)
    end

    @testset "All results compile" begin
        m = @enzyme_mechanism begin
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
        spec = _spec_with_rxn(
            m, uni_uni_rxn, uni_uni_with_reg)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_with_reg)
        for r in result
            compiled = EnzymeMechanism(r)
            @test compiled isa EnzymeMechanism
        end
    end

    @testset "exclude_regs works" begin
        m = @enzyme_mechanism begin
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
        spec = _spec_with_rxn(
            m, uni_uni_rxn, uni_uni_with_reg)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_with_reg;
            exclude_regs=Set([:I]))
        @test isempty(result)
    end

    @testset "Mirror steps created" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        bi_bi_with_reg = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        spec = _spec_with_rxn(
            m, bi_bi_rxn, bi_bi_with_reg)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, bi_bi_with_reg)
        # Find a variant with multiple eligible forms
        multi_form = filter(
            r -> length(r.steps) >
                length(spec.steps) + 2,
            result)
        if !isempty(multi_form)
            r = first(multi_form)
            @test length(r.steps) >
                length(spec.steps) + 2
        end
    end

    @testset "Equivalence constraints on binding K's" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        bi_bi_with_reg = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        spec = _spec_with_rxn(
            m, bi_bi_rxn, bi_bi_with_reg)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, bi_bi_with_reg)
        for r in result
            n_reg_binding = count(
                s -> length(s.reactants) == 2 &&
                    contains(string(s.reactants[2]),
                        "__reg"),
                r.steps)
            new_constraints = setdiff(
                r.param_constraints,
                spec.param_constraints)
            if n_reg_binding >= 2
                @test length(new_constraints) ==
                    n_reg_binding - 1
            else
                @test isempty(new_constraints)
            end
        end
    end

    @testset "Two regulators: both bound" begin
        m = @enzyme_mechanism begin
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
        rxn_ij = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I, J
        end
        spec = _spec_with_rxn(m, uni_uni_rxn, rxn_ij)
        # Add I
        i_specs =
            EnzymeRates._expand_add_dead_end_regulator(
                spec, rxn_ij)
        with_i = first(i_specs)
        # Add J to mechanism that already has I
        j_specs =
            EnzymeRates._expand_add_dead_end_regulator(
                with_i, rxn_ij)
        @test !isempty(j_specs)
        for s in j_specs
            has_j = any(
                contains(string(sym), "J__reg")
                for st in s.steps
                for sym in Iterators.flatten(
                    (st.reactants, st.products)))
            @test has_j
        end
    end

    @testset "Bi-bi: exact dead-end regulator count" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        rxn_i = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        spec = _spec_with_rxn(m, bi_bi_rxn, rxn_i)
        result =
            EnzymeRates._expand_add_dead_end_regulator(
                spec, rxn_i)
        # E_A_B has all subs → ineligible
        # E_P_Q has all prods → ineligible
        # 5 eligible forms, 9 inhibitor patterns
        # (3 sub subsets × 3 prod subsets) → 9
        @test length(result) == 9
    end

    @testset "Dead-end reg on allosteric spec" begin
        m = @enzyme_mechanism begin
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
        # Allosteric-only regulator → no dead-end
        rxn_allo_only = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        spec = _spec_with_rxn(
            m, uni_uni_rxn, rxn_allo_only)
        allo = first(EnzymeRates._expand_to_allosteric(
            spec, rxn_allo_only))
        result =
            EnzymeRates._expand_add_dead_end_regulator(
                allo, rxn_allo_only)
        @test isempty(result)

        # Mixed: I is dead-end, R is allosteric
        rxn_mixed = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
            allosteric_regulators: R
            oligomeric_state: 2
        end
        spec_m = _spec_with_rxn(
            m, uni_uni_rxn, rxn_mixed)
        allo_m = first(EnzymeRates._expand_to_allosteric(
            spec_m, rxn_mixed))
        result_m =
            EnzymeRates._expand_add_dead_end_regulator(
                allo_m, rxn_mixed)
        @test !isempty(result_m)
        for r in result_m
            @test r isa AllostericMechanismSpec
        end
    end

    @testset "Ping-pong: exact dead-end count" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[CX], B[N]
                products: P[C], Q[NX]
                enzymes: E, E_A[CX], E_Q[NX],
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
        rxn_pp_i = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products: P[C], Q[NX]
            dead_end_inhibitors: I
        end
        spec = _spec_with_rxn(
            m, bi_bi_pp_rxn, rxn_pp_i)
        result =
            EnzymeRates._expand_add_dead_end_regulator(
                spec, rxn_pp_i)
        # Estar_A_P has all prods → ineligible
        # Estar_B has all subs → ineligible
        # 4 eligible, 9 patterns → 3 unique form sets
        @test length(result) == 3
    end

    @testset "Topology-aware: sequential bi-bi" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_P_Q[CN], E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_A, B] ⇌ [E_A_B]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        rxn_seq = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        spec = _spec_with_rxn(m, bi_bi_rxn, rxn_seq)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, rxn_seq)
        # Sequential bi-bi: binding step sources:
        #   A: from E    B: from E_A
        #   P: from E_Q  Q: from E
        # Eligible forms: E, E_A, E_Q
        # 9 inhibitor patterns → 4 unique form sets:
        #   {E,E_Q}: ({A},{P}), ({A,B},{P})
        #   {E,E_A}: ({B},{Q}), ({B},{P,Q})
        #   {E_A,E_Q}: ({B},{P})
        #   {E}: ({A},{Q}), ({A},{P,Q}),
        #         ({A,B},{Q}), ({A,B},{P,Q})
        @test length(result) == 4
        for r in result
            @test r.param_count == spec.param_count + 1
        end
    end

    @testset "Bi-bi random: 9 inhibitor variants" begin
        m = @enzyme_mechanism begin
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
        bb_with_reg = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        spec = _spec_with_rxn(
            m, bi_bi_rxn, bb_with_reg)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, bb_with_reg)
        @test length(result) == 9
    end

    @testset "Two inhibitors: compete vs not" begin
        # Use bi-bi random: I1 binds to multiple forms,
        # creating mirror steps. When adding I2, compete
        # vs not-compete with I1 produces different forms.
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        rxn_2i = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I1, I2
        end
        spec = _spec_with_rxn(m, bi_bi_rxn, rxn_2i)
        # Adds both eligible regs (I1 + I2): 9 each
        result1 =
            EnzymeRates._expand_add_dead_end_regulator(
                spec, rxn_2i)
        @test length(result1) == 18
        # Pick variant where I1 binds to multiple forms
        multi = filter(result1) do r
            inh_forms = filter(
                f -> contains(string(f), "I1__reg"),
                collect(
                    EnzymeRates.all_form_names(r)))
            length(inh_forms) >= 2
        end
        @test !isempty(multi)
        spec_i1 = first(multi)
        # Add I2 to spec that already has I1
        result2 =
            EnzymeRates._expand_add_dead_end_regulator(
                spec_i1, rxn_2i)
        # With I1 present, competition patterns produce
        # more I2 variants than the base 9.
        @test length(result2) == 17
        # Not-competing variant: I2 coexists with I1
        has_coexist = any(result2) do r
            any(
                f -> contains(string(f), "I1__reg") &&
                     contains(string(f), "I2__reg"),
                collect(
                    EnzymeRates.all_form_names(r)))
        end
        @test has_coexist
        # Competing variant: I2 forms but no coexistence
        has_compete = any(result2) do r
            forms = collect(
                EnzymeRates.all_form_names(r))
            has_i2 = any(
                f -> contains(string(f), "I2__reg"),
                forms)
            no_coexist = !any(
                f -> contains(string(f), "I1__reg") &&
                     contains(string(f), "I2__reg"),
                forms)
            has_i2 && no_coexist
        end
        @test has_compete  # competing variant exists
    end
end

@testset "Regulator dummy naming stability" begin
    m = @enzyme_mechanism begin
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
    rxn2 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: I, J
    end
    base = mechanism_spec_from_mechanism(
        m, uni_uni_rxn)
    @test EnzymeMechanism(base) === m
    spec = MechanismSpec(rxn2, base.steps,
        base.param_constraints, base.param_count)
    # Add I first
    i_specs = EnzymeRates._expand_add_dead_end_regulator(
        spec, rxn2)
    with_i = first(filter(i_specs) do s
        any(contains(string(sym), "I__reg")
            for st in s.steps
            for sym in Iterators.flatten(
                (st.reactants, st.products)))
    end)
    # Now add J to a mechanism that already has I
    j_specs = EnzymeRates._expand_add_dead_end_regulator(
        with_i, rxn2)
    # J should use J__reg (no numeric suffix)
    for s in j_specs
        j_syms = [sym for st in s.steps
            for sym in Iterators.flatten(
                (st.reactants, st.products))
            if contains(string(sym), "J__reg")]
        for sym in j_syms
            @test !contains(string(sym), r"__reg\d")
        end
    end
end

@testset "Allosteric conversion" begin
    @testset "Uni-uni: K-type + V-type" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        # K-type: 1×1 = 1, V-type: 1. Total = 2
        @test length(result) == 2
        for r in result
            @test r isa AllostericMechanismSpec
            @test r.catalytic_n == 2
        end
    end

    @testset "Bi-bi: all substrate+product combos" begin
        bi_bi_allo = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        specs = EnzymeRates.init_mechanisms(bi_bi_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, bi_bi_allo)
        # K-type: 3×3 = 9, V-type: 1. Total = 10
        @test length(result) == 10
    end

    @testset "All are +1 param" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        for r in result
            @test r.param_count == spec.param_count + 1
        end
    end

    @testset "K-type: cat steps stay tr_equiv" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        k_type = filter(
            r -> !isempty(r.r_only_metabolites), result)
        @test !isempty(k_type)
        for r in k_type
            @test isempty(r.r_only_cat_steps)
            @test isempty(r.t_only_metabolites)
        end
    end

    @testset "V-type: all metabolites tr_equiv" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        v_type = filter(
            r -> !isempty(r.r_only_cat_steps) &&
                 isempty(r.r_only_metabolites),
            result)
        @test length(v_type) == 1
        sub_names = [s[1] for s in EnzymeRates.substrates(
            uni_uni_allo)]
        prod_names = [p[1] for p in EnzymeRates.products(
            uni_uni_allo)]
        for r in v_type
            for m in Symbol[sub_names; prod_names]
                @test m in r.tr_equiv_metabolites
            end
        end
    end

    @testset "All compile" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        for r in result
            m = AllostericEnzymeMechanism(r)
            @test m isa AllostericEnzymeMechanism
        end
    end

    @testset "Already allosteric → empty" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo = first(EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo))
        @test isempty(EnzymeRates._expand_to_allosteric(
            allo, uni_uni_allo))
    end

    @testset "oligomeric_state from reaction" begin
        rxn4 = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            oligomeric_state: 4
        end
        specs = EnzymeRates.init_mechanisms(rxn4)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, rxn4)
        for r in result
            @test r.catalytic_n == 4
        end
    end
end

@testset "Add allosteric regulator" begin
    @testset "Add regulator to allosteric spec" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo_reg)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo_reg)
        @test !isempty(allo_specs)
        allo = first(allo_specs)
        result = EnzymeRates._expand_add_allosteric_regulator(
            allo, uni_uni_allo_reg)
        # R not yet added: 3 flavors × 1 site option
        # (new site only, no existing reg sites) = 3
        @test length(result) == 3
    end

    @testset "Non-allosteric → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo_reg)
        spec = first(specs)
        result = EnzymeRates._expand_add_allosteric_regulator(
            spec, uni_uni_allo_reg)
        @test isempty(result)
    end

    @testset "Second regulator with site options" begin
        specs = EnzymeRates.init_mechanisms(
            uni_uni_allo_2reg)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo_2reg)
        allo = first(allo_specs)
        # Add R1 first
        r1_added = EnzymeRates._expand_add_allosteric_regulator(
            allo, uni_uni_allo_2reg)
        @test !isempty(r1_added)
        # Now add R2 to one with R1
        with_r1 = first(r1_added)
        r2_added = EnzymeRates._expand_add_allosteric_regulator(
            with_r1, uni_uni_allo_2reg)
        # R2: 3 flavors × 2 site options
        # (new site + R1's site) = 6
        @test length(r2_added) == 6
    end
end

@testset "Remove TR equivalence" begin
    @testset "Remove metabolite TR equiv" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        allo = first(allo_specs)
        n_tr = length(allo.tr_equiv_metabolites) +
               length(allo.tr_equiv_cat_steps)
        result = EnzymeRates._expand_remove_tr_equiv(
            allo, uni_uni_allo)
        @test length(result) == n_tr
    end

    @testset "No TR equivs left → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        allo = first(allo_specs)
        fully_relaxed = allo
        while true
            r = EnzymeRates._expand_remove_tr_equiv(
                fully_relaxed, uni_uni_allo)
            isempty(r) && break
            fully_relaxed = first(r)
        end
        @test isempty(
            EnzymeRates._expand_remove_tr_equiv(
                fully_relaxed, uni_uni_allo))
    end

    @testset "MechanismSpec → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_remove_tr_equiv(
            spec, uni_uni_allo)
        @test isempty(result)
    end

    @testset "V-type can remove r_only_cat_step" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        v_type = first(filter(
            r -> !isempty(r.r_only_cat_steps), allo_specs))
        @test !isempty(v_type.r_only_cat_steps)
        result = EnzymeRates._expand_remove_tr_equiv(
            v_type, uni_uni_allo)
        step_removals = filter(result) do r
            length(r.r_only_cat_steps) <
                length(v_type.r_only_cat_steps)
        end
        @test !isempty(step_removals)
        for r in step_removals
            @test r.param_count == v_type.param_count + 1
        end
    end

    @testset "Blocked when metabolites are r_only" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        k_type = first(filter(
            r -> !isempty(r.r_only_metabolites), allo_specs))
        mixed = AllostericMechanismSpec(
            k_type.base, k_type.catalytic_n,
            deepcopy(k_type.allosteric_reg_sites),
            copy(k_type.allosteric_multiplicities),
            copy(k_type.tr_equiv_metabolites),
            copy(k_type.tr_equiv_cat_steps),
            copy(k_type.r_only_metabolites),
            copy(k_type.t_only_metabolites),
            [1],  # r_only_cat_steps
            k_type.param_count)
        result = EnzymeRates._expand_remove_tr_equiv(
            mixed, uni_uni_allo)
        step_removals = filter(result) do r
            length(r.r_only_cat_steps) <
                length(mixed.r_only_cat_steps)
        end
        @test isempty(step_removals)
    end

    @testset "TR equiv removal delta for allosteric regulators" begin
        rxn_r = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        specs = EnzymeRates.init_mechanisms(rxn_r)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(spec, rxn_r)
        allo = first(allo_specs)
        reg_specs = EnzymeRates._expand_add_allosteric_regulator(allo, rxn_r)
        tr_spec = first(filter(r -> :R in r.tr_equiv_metabolites, reg_specs))
        pc_before = tr_spec.param_count
        result = EnzymeRates._expand_remove_tr_equiv(tr_spec, rxn_r)
        r_removal = filter(result) do r
            :R ∉ r.tr_equiv_metabolites &&
            :R ∉ r.r_only_metabolites &&
            :R ∉ r.t_only_metabolites
        end
        @test !isempty(r_removal)
        for r in r_removal
            @test r.param_count == pc_before + 1
        end
    end

    @testset "TR equiv removal delta skips constrained follower steps" begin
        # bi-bi random has 2 binding steps for A with K_follower = K_leader
        # constraint. Removing TR equiv for A should add +1 (one K_A_T),
        # not +2 (which would count both binding steps independently).
        const_bi_bi = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        specs = EnzymeRates.init_mechanisms(const_bi_bi)
        # Find a spec with 2 A-binding steps (constrained equal)
        two_a = filter(specs) do s
            length(filter(st -> EnzymeRates.step_metabolite(st) === :A,
                s.steps)) >= 2
        end
        @test !isempty(two_a)
        spec = first(two_a)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, const_bi_bi)
        # Find an allosteric spec with A in tr_equiv_metabolites
        tr_specs = filter(s -> :A in s.tr_equiv_metabolites, allo_specs)
        @test !isempty(tr_specs)
        tr_spec = first(tr_specs)
        result = EnzymeRates._expand_remove_tr_equiv(
            tr_spec, const_bi_bi)
        a_removals = filter(result) do r
            :A ∉ r.tr_equiv_metabolites &&
            :A ∉ r.r_only_metabolites &&
            :A ∉ r.t_only_metabolites
        end
        @test !isempty(a_removals)
        for r in a_removals
            m = EnzymeRates.compile_mechanism(r)
            @test length(parameters(m)) == r.param_count
        end
    end
end

@testset "Dedup" begin
    @testset "Same mechanism, different step order" begin
        spec1 = MechanismSpec(
            uni_uni_rxn,
            [StepSpec([:E, :S], [:E_S], true),
             StepSpec([:E, :P], [:E_P], true),
             StepSpec([:E_S], [:E_P], false)],
            ParamConstraint[], 5)
        spec2 = MechanismSpec(
            uni_uni_rxn,
            [StepSpec([:E, :P], [:E_P], true),
             StepSpec([:E_S], [:E_P], false),
             StepSpec([:E, :S], [:E_S], true)],
            ParamConstraint[], 5)
        cache = Dict(5 => AbstractMechanismSpec[spec1, spec2])
        EnzymeRates.dedup!(cache)
        @test length(cache[5]) == 1
    end

    @testset "Different mechanisms preserved" begin
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        pc = first(specs).param_count
        cache = Dict(pc => AbstractMechanismSpec[specs...])
        EnzymeRates.dedup!(cache)
        @test length(cache[pc]) >= 1
        @test length(cache[pc]) <= length(specs)
    end

    @testset "Idempotent" begin
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        pc = first(specs).param_count
        cache = Dict(pc => AbstractMechanismSpec[specs...])
        EnzymeRates.dedup!(cache)
        n1 = length(cache[pc])
        EnzymeRates.dedup!(cache)
        @test length(cache[pc]) == n1
    end

    @testset "Allosteric dedup: site order" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo = first(EnzymeRates._expand_to_allosteric(spec, uni_uni_allo))
        spec_ab = AllostericMechanismSpec(
            allo.base, allo.catalytic_n,
            [[:A], [:B]], [2, 2],
            copy(allo.tr_equiv_metabolites),
            copy(allo.tr_equiv_cat_steps),
            copy(allo.r_only_metabolites),
            copy(allo.t_only_metabolites),
            copy(allo.r_only_cat_steps),
            allo.param_count + 2)
        spec_ba = AllostericMechanismSpec(
            allo.base, allo.catalytic_n,
            [[:B], [:A]], [2, 2],
            copy(allo.tr_equiv_metabolites),
            copy(allo.tr_equiv_cat_steps),
            copy(allo.r_only_metabolites),
            copy(allo.t_only_metabolites),
            copy(allo.r_only_cat_steps),
            allo.param_count + 2)
        pc = spec_ab.param_count
        cache = Dict(pc => EnzymeRates.AbstractMechanismSpec[spec_ab, spec_ba])
        EnzymeRates.dedup!(cache)
        @test length(cache[pc]) == 1
    end
end

@testset "expand_mechanisms" begin
    @testset "Returns dict keyed by param count" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        result = EnzymeRates.expand_mechanisms(
            specs, uni_uni_rxn)
        @test result isa Dict{Int,
            Vector{AbstractMechanismSpec}}
        base_pc = first(specs).param_count
        @test haskey(result, base_pc + 1)
    end

    @testset "Allosteric expansion included" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        result = EnzymeRates.expand_mechanisms(
            specs, uni_uni_allo)
        base_pc = first(specs).param_count
        has_allo = any(
            any(s isa AllostericMechanismSpec
                for s in ss)
            for (_, ss) in result)
        @test has_allo
    end

    @testset "No self-expansion to same param count" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        base_pc = first(specs).param_count
        result = EnzymeRates.expand_mechanisms(
            specs, uni_uni_rxn)
        # All results should have param_count > base
        for (pc, _) in result
            @test pc > base_pc
        end
    end

    @testset "Allosteric rewrap preserves structure" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        allo = first(allo_specs)
        result = EnzymeRates.expand_mechanisms(
            [allo], uni_uni_allo)
        # Should have expansions from base moves (RE→SS)
        # rewrapped as AllostericMechanismSpec
        has_rewrapped = any(
            any(s isa AllostericMechanismSpec
                for s in ss)
            for (_, ss) in result)
        @test has_rewrapped
    end

    @testset "Dead-end excludes allosteric regs" begin
        specs = EnzymeRates.init_mechanisms(
            uni_uni_allo_reg)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo_reg)
        allo = first(allo_specs)
        # Add R as allosteric regulator
        with_reg = first(
            EnzymeRates._expand_add_allosteric_regulator(
                allo, uni_uni_allo_reg))
        result = EnzymeRates.expand_mechanisms(
            [with_reg], uni_uni_allo_reg)
        # R should NOT appear as dead-end in any expansion
        for (_, ss) in result
            for s in ss
                base = s isa AllostericMechanismSpec ?
                    s.base : s
                for step in base.steps
                    for sym in Iterators.flatten(
                            (step.reactants, step.products))
                        @test !contains(
                            string(sym), "R__reg")
                    end
                end
            end
        end
    end
end

@testset "Integration" begin
    @testset "Uni-uni full enumeration" begin
        results = enumerate_all(uni_uni_rxn; max_params=8)
        @test !isempty(results)
        pcs = sort(collect(keys(results)))
        @test issorted(pcs)
        # Every mechanism compiles
        allo_count = 0
        for (pc, specs) in results
            for spec in specs
                if spec isa EnzymeRates.MechanismSpec
                    m = EnzymeMechanism(spec)
                    @test length(parameters(m)) <= pc
                elseif spec isa EnzymeRates.AllostericMechanismSpec
                    allo_count += 1
                    allo_count > 3 && continue
                    m = AllostericEnzymeMechanism(spec)
                    @test length(parameters(m)) <= spec.param_count
                end
            end
        end
    end

    @testset "Bi-bi full enumeration" begin
        results = enumerate_all(bi_bi_rxn; max_params=10)
        @test !isempty(results)
        allo_count = 0
        for (pc, specs) in results
            for spec in specs
                if spec isa EnzymeRates.MechanismSpec
                    m = EnzymeMechanism(spec)
                    @test length(parameters(m)) <= pc
                elseif spec isa EnzymeRates.AllostericMechanismSpec
                    allo_count += 1
                    allo_count > 3 && continue
                    m = AllostericEnzymeMechanism(spec)
                    @test length(parameters(m)) <= spec.param_count
                end
            end
        end
    end

    @testset "With allosteric regulators" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        results = enumerate_all(rxn; max_params=10)
        has_allo = any(
            any(s isa EnzymeRates.AllostericMechanismSpec for s in specs)
            for (_, specs) in results)
        @test has_allo
        # Every mechanism compiles
        allo_count = 0
        for (pc, specs) in results
            for spec in specs
                if spec isa EnzymeRates.MechanismSpec
                    m = EnzymeMechanism(spec)
                    @test length(parameters(m)) <= pc
                elseif spec isa EnzymeRates.AllostericMechanismSpec
                    allo_count += 1
                    allo_count > 3 && continue
                    m = AllostericEnzymeMechanism(spec)
                    @test length(parameters(m)) <= spec.param_count
                end
            end
        end
    end

    @testset "With dead-end regulator" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
        end
        results = enumerate_all(rxn; max_params=8)
        @test !isempty(results)
        # Should have more mechanisms than plain uni-uni
        plain = enumerate_all(uni_uni_rxn; max_params=8)
        total_with_reg = sum(length(v) for v in values(results))
        total_plain = sum(length(v) for v in values(plain))
        @test total_with_reg > total_plain
    end

    @testset "Multiple levels populated" begin
        results = enumerate_all(uni_uni_rxn; max_params=8)
        @test length(results) >= 2  # At least 2 param count levels
    end
end

@testset "r_only params excluded from parameter list" begin
    specs = EnzymeRates.init_mechanisms(uni_uni_allo)
    spec = first(specs)
    allo_specs = EnzymeRates._expand_to_allosteric(
        spec, uni_uni_allo)

    @testset "K-type: no K_T params for r_only metabolites" begin
        k_type = first(filter(
            r -> !isempty(r.r_only_metabolites), allo_specs))
        m = AllostericEnzymeMechanism(k_type)
        params = parameters(m)
        @test length(params) == k_type.param_count
        t_params = filter(
            p -> endswith(string(p), "_T"), params)
        @test isempty(t_params)
    end

    @testset "V-type: no kf_T/kr_T for r_only cat steps" begin
        v_type = first(filter(
            r -> !isempty(r.r_only_cat_steps), allo_specs))
        m = AllostericEnzymeMechanism(v_type)
        params = parameters(m)
        @test length(params) == v_type.param_count
        t_k_params = filter(
            p -> contains(string(p), "f_T") ||
                 contains(string(p), "r_T"), params)
        @test isempty(t_k_params)
    end
end

@testset "Metabolite overlap: substrate as dead-end inhibitor" begin
    rxn_overlap = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
    end
    specs = EnzymeRates.init_mechanisms(rxn_overlap)
    @test !isempty(specs)

    spec = first(specs)
    de_specs = EnzymeRates._expand_add_dead_end_regulator(
        spec, rxn_overlap)
    @test !isempty(de_specs)

    # S-as-regulator uses __reg suffix
    for s in de_specs
        reg_syms = [sym for st in s.steps
            for sym in Iterators.flatten(
                (st.reactants, st.products))
            if contains(string(sym), "S__reg")]
        @test !isempty(reg_syms)
    end

    # All compile correctly
    for s in de_specs
        m = EnzymeMechanism(s)
        @test m isa EnzymeMechanism
    end
end

@testset "Metabolite overlap: substrate as allosteric regulator" begin
    rxn_allo_overlap = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: S
        oligomeric_state: 2
    end
    specs = EnzymeRates.init_mechanisms(rxn_allo_overlap)
    @test !isempty(specs)
    spec = first(specs)

    # Allosteric conversion works
    allo_specs = EnzymeRates._expand_to_allosteric(
        spec, rxn_allo_overlap)
    @test !isempty(allo_specs)

    # Add S as allosteric regulator
    allo = first(allo_specs)
    reg_specs = EnzymeRates._expand_add_allosteric_regulator(
        allo, rxn_allo_overlap)
    @test !isempty(reg_specs)

    # S appears in allosteric_reg_sites
    for r in reg_specs
        has_s = any(:S in site for site in r.allosteric_reg_sites)
        @test has_s
    end

    # All compile correctly
    for r in reg_specs
        m = AllostericEnzymeMechanism(r)
        @test m isa AllostericEnzymeMechanism
    end

    # TR equiv removal produces separate results for
    # S-as-substrate and S-as-regulator
    tr_spec = first(filter(
        r -> :S in r.tr_equiv_metabolites, reg_specs))
    result = EnzymeRates._expand_remove_tr_equiv(
        tr_spec, rxn_allo_overlap)
    # S as catalytic met and S as regulator are both
    # in tr_equiv_metabolites. Removing each should
    # produce separate variants.
    @test length(result) >= 2
end

@testset "Base-level moves on allosteric specs" begin
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
    spec = mechanism_spec_from_mechanism(m_uu, uni_uni_allo)
    @test EnzymeMechanism(spec) === m_uu
    allo_specs = EnzymeRates._expand_to_allosteric(
        spec, uni_uni_allo)
    allo = first(allo_specs)

    @testset "RE→SS on allosteric" begin
        result = EnzymeRates._expand_re_to_ss(allo)
        @test !isempty(result)
        for r in result
            @test r isa EnzymeRates.AllostericMechanismSpec
            @test r.catalytic_n == allo.catalytic_n
            @test r.allosteric_reg_sites ==
                allo.allosteric_reg_sites
            @test r.r_only_metabolites ==
                allo.r_only_metabolites
            @test r.param_count > allo.param_count
        end
    end

    @testset "Remove constraint on allosteric" begin
        m_bb = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        bi_bi_allo_rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        bb_spec = mechanism_spec_from_mechanism(
            m_bb, bi_bi_allo_rxn)
        @test EnzymeMechanism(bb_spec) === m_bb
        bb_spec_c = EnzymeRates.MechanismSpec(
            bb_spec.reaction, bb_spec.steps,
            EnzymeRates._max_equivalence_constraints(
                bb_spec),
            bb_spec.param_count)
        bb_allo = first(
            EnzymeRates._expand_to_allosteric(
                bb_spec_c, bi_bi_allo_rxn))
        @test !isempty(bb_allo.base.param_constraints)
        result = EnzymeRates._expand_remove_constraint(
            bb_allo)
        @test !isempty(result)
        for r in result
            @test r isa
                EnzymeRates.AllostericMechanismSpec
            @test r.catalytic_n == bb_allo.catalytic_n
            @test r.r_only_metabolites ==
                bb_allo.r_only_metabolites
        end
    end

    @testset "Add dead-end reg on allosteric" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
            oligomeric_state: 2
        end
        spec_i = mechanism_spec_from_mechanism(
            m_uu, rxn)
        allo_i = first(
            EnzymeRates._expand_to_allosteric(
                spec_i, rxn))
        result =
            EnzymeRates._expand_add_dead_end_regulator(
                allo_i, rxn)
        @test !isempty(result)
        for r in result
            @test r isa
                EnzymeRates.AllostericMechanismSpec
            @test r.catalytic_n == allo_i.catalytic_n
            @test r.r_only_metabolites ==
                allo_i.r_only_metabolites
        end
    end
end

end # top-level testset
