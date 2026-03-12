# ABOUTME: Tests for the staged mechanism enumeration pipeline
# ABOUTME: Organized by stage with hand-calculated expected values

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
            ) >= 1
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
            ) >= 1
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
            ) >= 1
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
        @test length(topos) == 19

        # Topo 16: classic ping-pong with Estar
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
        @test compile_mechanism(topos[16]) === m_pp

        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps
            ) >= 1
        end
    end

    @testset "Ter-Ter" begin
        topos = EnzymeRates._catalytic_topologies(
            ter_ter)
        @test length(topos) >= 1
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps
            ) >= 1
        end
    end

end

@testset "Stage 2: RE/SS expansion" begin
    @testset "Uni-Uni: 2^3 - 1 = 7 variants" begin
        topo = EnzymeRates._catalytic_topologies(uni_uni)[1]
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
        topo = EnzymeRates._catalytic_topologies(uni_bi)[1]
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
        topo = EnzymeRates._catalytic_topologies(bi_bi)[1]
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
        topo = EnzymeRates._catalytic_topologies(bi_bi)[end]
        result_default = EnzymeRates._expand_ress_variants(
            [topo], bi_bi)
        # With strict max_re_groups=2, fewer variants survive
        result_strict = EnzymeRates._expand_ress_variants(
            [topo], bi_bi; max_re_groups=2)
        @test length(result_strict) <=
            length(result_default)
    end
end

@testset "Stage 2.5: Substrate/product dead-end expansion" begin

    @testset "Uni-Uni: passthrough (no off-cycle forms)" begin
        # Uni-Uni: 3 forms (E, E_S, E_P), all on-cycle.
        # E_S has all subs, E_P has all prods, and E
        # can only bind S→E_S or P→E_P (both catalytic).
        # → 0 dead-end forms, passthrough.
        topo = EnzymeRates._catalytic_topologies(
            uni_uni)[1]
        ress = EnzymeRates._expand_ress_variants(
            [topo], uni_uni)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                ress, uni_uni)
        @test length(result) == length(ress)
    end

    @testset "Bi-Bi random: 4 dead-end forms" begin
        # Fully-random Bi-Bi (topo 9, 7 forms):
        # E, E_A, E_B, E_A_B, E_P, E_Q, E_P_Q
        # Eligible forms (not all-subs or all-prods):
        #   E: all 4 mets → catalytic forms (no dead-ends)
        #   E_A: +P→E_A_P, +Q→E_A_Q (valid)
        #   E_B: +P→E_B_P, +Q→E_B_Q (valid)
        #   E_P: +A→E_A_P, +B→E_B_P (valid, not
        #         all subs since 2 subs needed)
        #   E_Q: +A→E_A_Q, +B→E_B_Q (valid)
        # 4 unique dead-end forms → 2^4 = 16 variants
        topo = EnzymeRates._catalytic_topologies(
            bi_bi)[end]
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [topo], bi_bi)
        @test length(result) == 16
    end

    @testset "Uni-Bi ordered: 1 dead-end form" begin
        # Uni-Bi topo 1 (ordered Q-first release):
        # E, E_S, E_P_Q, E_Q (4 forms)
        # Eligible: E (not all subs/prods)
        #   E+S→E_S (catalytic), E+Q→E_Q (catalytic)
        #   E+P→E_P: not catalytic, 1 prod bound,
        #   not all prods → valid dead-end
        # E_Q: +S→all subs+prod → invalid
        #       +P→E_P_Q catalytic
        # 2^1 = 2 variants per input spec
        topo = EnzymeRates._catalytic_topologies(
            uni_bi)[1]
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [topo], uni_bi)
        @test length(result) == 2
        # The original (no dead-ends) must be included
        original_forms = EnzymeRates.all_form_names(topo)
        has_original = any(result) do spec
            EnzymeRates.all_form_names(spec) ==
                original_forms
        end
        @test has_original
    end
end

@testset "Stage 3: Regulator dead-end expansion" begin
    @testset "No regulators: passthrough" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni)[1]
        result = EnzymeRates._expand_dead_end(
            [topo], uni_uni; dead_end_regs=Symbol[])
        @test length(result) == 1
    end

    @testset "Uni-Uni + I: eligible forms" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I)[1]
        result = EnzymeRates._expand_dead_end(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        # Only E eligible (ES fully occupied with all
        # substrates, EP fully occupied with all
        # products) → (2^1)^1 = 2 variants
        @test length(result) == 2
    end

    @testset "Uni-Bi + I" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_bi_dead_end_I)[1]
        result = EnzymeRates._expand_dead_end(
            [topo], uni_bi_dead_end_I;
            dead_end_regs=[:I])
        # Topo[1] Q-first: E, ES, EQ, EPQ (4 forms).
        # ES has all subs, EPQ has all prods
        # → eligible = {E, EQ} = 2 forms
        # Sub/prod dead-ends: E_P (from E)
        # Regulator dead-ends: E_I__reg1, E_I__reg1_Q
        # Combined: 3 dead-end forms → 2^3 = 8
        @test length(result) == 8
    end

    @testset "Bi-Bi + I" begin
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_dead_end_I)[1]
        result = EnzymeRates._expand_dead_end(
            [topo], bi_bi_dead_end_I;
            dead_end_regs=[:I])
        # Topo[1] sequential (5 forms: E, EA, EAB,
        # EQ, EPQ). EAB has all subs, EPQ all prods
        # → eligible = {E, EA, EQ} = 3 forms
        # Sub/prod dead-ends: E_B, E_P, E_A_P,
        #   E_A_Q, E_B_Q (5 forms)
        # Regulator dead-ends: E_I__reg1,
        #   E_A_I__reg1, E_I__reg1_Q (3 forms)
        # Combined: 8 dead-end forms → 2^8 = 256
        @test length(result) == 256
    end

    @testset "2 inhibitors: Uni-Uni + I, J" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I_J)[1]
        result = EnzymeRates._expand_dead_end(
            [topo], uni_uni_dead_end_I_J;
            dead_end_regs=[:I, :J])
        # Only E eligible → (2^2)^1 = 4
        @test length(result) == 4
    end

    @testset "Allosteric-only: passthrough" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_allosteric_R)[1]
        result = EnzymeRates._expand_dead_end(
            [topo], uni_uni_allosteric_R;
            dead_end_regs=Symbol[])
        @test length(result) == 1
    end

    @testset "Properties" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I)[1]
        result = EnzymeRates._expand_dead_end(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        for s in result
            @test length(s.steps) >= length(topo.steps)
            @test s.param_count >= topo.param_count
        end
    end
end

@testset "Stage 4: Equivalence constraints" begin
    @testset "No equiv groups: passthrough" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni)[1]
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
        topo = EnzymeRates._catalytic_topologies(
            bi_bi)[1]
        de = EnzymeRates._expand_dead_end(
            [topo], bi_bi; dead_end_regs=Symbol[])
        # Find a spec with exactly 1 equiv group
        function _n_equiv_groups(spec)
            gs = Dict{Tuple{Symbol,Bool},
                Vector{Int}}()
            for (j, st) in enumerate(spec.steps)
                met = EnzymeRates.step_metabolite(st)
                met === nothing && continue
                key = (met, st.is_equilibrium)
                push!(get!(gs, key, Int[]), j)
            end
            count(v -> length(v) >= 2,
                values(gs))
        end
        s = first(
            x for x in de
            if _n_equiv_groups(x) == 1)
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
        topo = EnzymeRates._catalytic_topologies(
            bi_bi)[1]
        de = EnzymeRates._expand_dead_end(
            [topo], bi_bi; dead_end_regs=Symbol[])
        # Find spec with ≥ 2 equiv groups
        function _count_equiv_groups(spec)
            gs = Dict{Tuple{Symbol,Bool},
                Vector{Int}}()
            for (j, st) in enumerate(spec.steps)
                met = EnzymeRates.step_metabolite(st)
                met === nothing && continue
                key = (met, st.is_equilibrium)
                push!(get!(gs, key, Int[]), j)
            end
            count(v -> length(v) >= 2,
                values(gs))
        end
        multi = filter(
            x -> _count_equiv_groups(x) >= 2, de)
        @test !isempty(multi)
        s = first(multi)
        n_groups = _count_equiv_groups(s)
        eq = EnzymeRates._expand_equivalence_constraints(
            [s], bi_bi)
        @test length(eq) == 2^n_groups
    end

    @testset "Substrate/regulator same metabolite" begin
        # Dummy names (:I__reg1) prevent grouping
        # regulator binding steps with catalytic steps
        # that bind a different real metabolite. With
        # Uni-Uni + I, the dead-end spec has :S and
        # :I__reg1 as distinct metabolites → no equiv
        # group → only 1 variant (passthrough).
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I)[1]
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
        topo = EnzymeRates._catalytic_topologies(
            bi_bi)[1]
        de = EnzymeRates._expand_dead_end(
            [topo], bi_bi; dead_end_regs=Symbol[])
        # Find specs with substrate/product dead-ends
        has_de = filter(
            x -> length(x.steps) > length(topo.steps),
            de)
        if !isempty(has_de)
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
    end

    @testset "Properties" begin
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_dead_end_I)[1]
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
        @test minimum(
            s.param_count for s in constrained) <
            maximum(
                s.param_count
                for s in unconstrained)
    end
end

@testset "Stage 5: Deduplication" begin
    @testset "Uni-Uni dedup" begin
        topo = EnzymeRates._catalytic_topologies(uni_uni)[1]
        ress = EnzymeRates._expand_ress_variants(
            [topo], uni_uni)
        deduped = EnzymeRates._deduplicate(ress, uni_uni)
        # All RE/SS variants of single Uni-Uni topology
        # produce same concentration fingerprint → dedup to 1
        @test length(deduped) == 1
        @test deduped[1].param_count <=
            minimum(s.param_count for s in ress)
    end

    @testset "Duplicate removal" begin
        topo = EnzymeRates._catalytic_topologies(uni_uni)[1]
        result = EnzymeRates._deduplicate(
            [topo, deepcopy(topo)], uni_uni)
        @test length(result) == 1
    end

    @testset "Bi-Bi ordered: single-SS fingerprints" begin
        # For an ordered Bi-Bi topology (5 edges), create
        # single-SS variants by toggling each RE edge to SS
        # one at a time. The topo from _catalytic_topologies has
        # 4 RE edges and 1 SS (isomerization). Only 4 variants
        # are built here because _expand_ress_variants has a bug
        # where non-binding (isomerization) edges cannot be
        # toggled between RE and SS — all 5 edges should be
        # individually toggleable.
        topo = EnzymeRates._catalytic_topologies(bi_bi)[1]
        n = length(topo.edges)
        n_forms = length(Set(Iterators.flatten(topo.edges)))
        n_thermo = n - n_forms + 1
        single_ss = EnzymeRates.MechanismSpec[]
        for i in 1:n
            if topo.equilibrium_steps[i]
                eq = copy(topo.equilibrium_steps)
                eq[i] = false
                n_re = count(eq)
                n_ss = n - n_re
                pc = n_re + 2 * n_ss - n_thermo + 2
                push!(single_ss, EnzymeRates.MechanismSpec(
                    topo.reaction, topo.edges,
                    topo.n_catalytic_edges, eq,
                    topo.param_constraints,
                    pc))
            end
        end
        # Bug: _expand_ress_variants doesn't toggle non-binding
        # edges (documented in Stage 2 tests), so only 4 RE
        # edges are toggleable instead of all 5 edges.
        @test length(single_ss) == 4
        deduped = EnzymeRates._deduplicate(
            single_ss, bi_bi)
        # 4 single-SS variants dedup to 3: one pair
        # shares a concentration fingerprint.
        # Correct: all 5 single-SS variants (including
        # isom) have distinct fingerprints → dedup to 5.
        @test_broken length(deduped) == 5
        @test length(deduped) == 3
    end

    @testset "Keeps lower param_count" begin
        topo = EnzymeRates._catalytic_topologies(uni_uni)[1]
        topo2 = EnzymeRates.MechanismSpec(
            topo.reaction, topo.edges,
            topo.n_catalytic_edges,
            topo.equilibrium_steps,
            topo.param_constraints,
            topo.param_count + 5)
        result = EnzymeRates._deduplicate(
            [topo, topo2], uni_uni)
        @test length(result) == 1
        @test result[1].param_count == topo.param_count
    end
end

@testset "Stage 6: Allosteric expansion" begin
    @testset "No allosteric regs: passthrough" begin
        topo = EnzymeRates._catalytic_topologies(uni_uni)[1]
        dd = EnzymeRates._deduplicate([topo], uni_uni)
        result = EnzymeRates._expand_allosteric(
            dd, uni_uni; catalytic_n=1,
            allosteric_regs=Symbol[])
        @test length(result) == 1
    end

    @testset "1 reg, catalytic_n=1" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_allosteric_R)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        result = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R; catalytic_n=1,
            allosteric_regs=[:R])
        @test length(result) == 1  # m=1 only
    end

    @testset "1 reg, catalytic_n=2" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_bi_allosteric_R)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_bi_allosteric_R)
        result = EnzymeRates._expand_allosteric(
            dd, uni_bi_allosteric_R; catalytic_n=2,
            allosteric_regs=[:R])
        @test length(result) == 2  # m=1, m=2
    end

    @testset "1 reg, catalytic_n=3" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_allosteric_R)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        result = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R; catalytic_n=3,
            allosteric_regs=[:R])
        @test length(result) == 3  # m=1,2,3
    end

    @testset "2 regs, catalytic_n=1" begin
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_allosteric_R1_R2)[1]
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
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_allosteric_R1_R2)[1]
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
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I)[1]
        de = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        dd = EnzymeRates._deduplicate(
            de, uni_uni_dead_end_I)
        result = EnzymeRates._expand_allosteric(
            dd, uni_uni_dead_end_I; catalytic_n=1,
            allosteric_regs=Symbol[])
        @test length(result) == length(dd)
    end

    @testset "Dead-end I + allosteric R" begin
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_dead_end_I_allosteric_R)[1]
        # Expand dead-end inhibitors first so edges exist
        de = EnzymeRates._expand_dead_end_inhibitors(
            [topo], bi_bi_dead_end_I_allosteric_R;
            dead_end_regs=[:I])
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
        topo = EnzymeRates._catalytic_topologies(
            uni_bi_allosteric_R)[1]
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
    @testset "Uni-Uni + R: 2^3 = 8 variants" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_allosteric_R)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        allo = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R; catalytic_n=1,
            allosteric_regs=[:R])
        tr = EnzymeRates._expand_tr_equivalence(
            allo, uni_uni_allosteric_R)
        # 3 metabolites with T-state params: S, P, R
        @test length(tr) == 8
    end

    @testset "Uni-Bi + R: 2^4 = 16 variants" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_bi_allosteric_R)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_bi_allosteric_R)
        allo = EnzymeRates._expand_allosteric(
            dd, uni_bi_allosteric_R; catalytic_n=1,
            allosteric_regs=[:R])
        tr = EnzymeRates._expand_tr_equivalence(
            allo, uni_bi_allosteric_R)
        # 4 metabolites: S, P, Q, R
        @test length(tr) == 16
    end

    @testset "Bi-Bi + R1, R2: 2^6 = 64 per spec" begin
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_allosteric_R1_R2)[1]
        dd = EnzymeRates._deduplicate(
            [topo], bi_bi_allosteric_R1_R2)
        allo = EnzymeRates._expand_allosteric(
            dd, bi_bi_allosteric_R1_R2;
            catalytic_n=1, allosteric_regs=[:R1, :R2])
        # Use first allosteric spec only
        tr = EnzymeRates._expand_tr_equivalence(
            [allo[1]], bi_bi_allosteric_R1_R2)
        # 6 metabolites: A, B, P, Q, R1, R2
        @test length(tr) == 64
    end

    @testset "Properties" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_allosteric_R)[1]
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
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_allosteric_R1_R2)[1]
        dd = EnzymeRates._deduplicate(
            [topo], bi_bi_allosteric_R1_R2)
        allo = EnzymeRates._expand_allosteric(
            dd, bi_bi_allosteric_R1_R2;
            catalytic_n=1, allosteric_regs=[:R1, :R2])
        # Use first allosteric spec
        tr = EnzymeRates._expand_tr_equivalence(
            [allo[1]], bi_bi_allosteric_R1_R2)
        @test length(tr) == 64  # 2^6
        deduped = EnzymeRates._deduplicate_allosteric(
            tr, bi_bi_allosteric_R1_R2)
        # Bug: _allosteric_canonical_key uses sort(tr_equiv_metabolites)
        # directly, so complementary TR-equiv sets (T↔R mirrors) get
        # different keys and both survive dedup. Mirrors should be removed.
        @test_broken length(deduped) < length(tr)
        @test length(deduped) == length(tr)
    end

    @testset "Uni-Uni + R: no mirrors (odd metabolites)" begin
        # 3 metabolites (S, P, R): complementary subsets always
        # differ in size → different param counts → no true mirrors.
        # (Additionally, _allosteric_canonical_key never maps
        # complements to the same key, but even with correct dedup,
        # odd metabolite counts produce no equal-size mirrors.)
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_allosteric_R)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        allo = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R;
            catalytic_n=1, allosteric_regs=[:R])
        tr = EnzymeRates._expand_tr_equivalence(
            allo, uni_uni_allosteric_R)
        deduped = EnzymeRates._deduplicate_allosteric(
            tr, uni_uni_allosteric_R)
        @test length(deduped) == length(tr)
    end

    @testset "Keeps lower param_count on mirror" begin
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_allosteric_R1_R2)[1]
        dd = EnzymeRates._deduplicate(
            [topo], bi_bi_allosteric_R1_R2)
        allo = EnzymeRates._expand_allosteric(
            dd, bi_bi_allosteric_R1_R2;
            catalytic_n=1, allosteric_regs=[:R1, :R2])
        tr = EnzymeRates._expand_tr_equivalence(
            [allo[1]], bi_bi_allosteric_R1_R2)
        deduped = EnzymeRates._deduplicate_allosteric(
            tr, bi_bi_allosteric_R1_R2)
        # Bug: mirror dedup is broken — _allosteric_canonical_key
        # never maps complementary TR-equiv sets to the same key,
        # so no mirrors are removed and param_count selection
        # is untestable. Already covered by @test_broken in
        # "T/R mirrors dedup" above.
        for s in deduped
            @test s isa EnzymeRates.AllostericMechanismSpec
        end
    end

    @testset "Different base mechanisms survive" begin
        topos = EnzymeRates._catalytic_topologies(
            uni_bi_allosteric_R)
        @test length(topos) >= 2
        dd = EnzymeRates._deduplicate(
            topos[1:2], uni_bi_allosteric_R)
        @test length(dd) >= 2
        if length(dd) >= 2
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
end

@testset "Cross-stage properties" begin
    @testset "RE partition bounds" begin
        catalytic = EnzymeRates._catalytic_topologies(
            uni_bi)
        spec = catalytic[1]
        ress = EnzymeRates._expand_ress_variants(
            [spec], uni_bi)
        for s in ress
            partition =
                EnzymeRates._compute_re_partition(
                    s.edges, s.equilibrium_steps)
            @test 2 <= length(partition) <= 7
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
            EnzymeRates._expand_dead_end_inhibitors(
                ress_r, uni_bi_dead_end_I;
                dead_end_regs=[:I])
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
        for rxn in [uni_uni, bi_bi]
            all_specs = collect(
                EnzymeRates.enumerate_mechanisms(rxn))
            cat_specs = filter(
                s -> s isa EnzymeRates.MechanismSpec, all_specs)
            for s in first(cat_specs, 3)
                m = compile_mechanism(s)
                @test m isa EnzymeMechanism
                s2 = mechanism_spec_from_mechanism(
                    m, rxn; n_catalytic_edges=s.n_catalytic_edges)
                @test s2.edges == s.edges
                @test s2.equilibrium_steps == s.equilibrium_steps
                @test s2.param_count == s.param_count
            end
        end

        @testset "allosteric compilation" begin
            all_specs = collect(
                EnzymeRates.enumerate_mechanisms(
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
        # 1 topology, all RE/SS variants dedup to 1
        result = collect(
            EnzymeRates.enumerate_mechanisms(uni_uni))
        @test length(result) == 1
    end

    @testset "Uni-Bi, no regs" begin
        # 3 topologies × RE/SS variants, dedup to 9
        result = collect(
            EnzymeRates.enumerate_mechanisms(uni_bi))
        @test length(result) == 9
    end

    @testset "Bi-Bi, no regs" begin
        # 9 topologies → 81 unconstrained + 126
        # constrained = 207 after dedup + equiv expansion
        result = collect(
            EnzymeRates.enumerate_mechanisms(bi_bi))
        @test length(result) == 207
    end

    @testset "Bi-Bi Ping-Pong, no regs" begin
        # 10 topologies → 84 unconstrained + 126
        # constrained = 210 after dedup + equiv expansion
        result = collect(
            EnzymeRates.enumerate_mechanisms(
                bi_bi_ping_pong))
        @test length(result) == 210
    end

    @testset "Uni-Uni + 1 unknown reg" begin
        # 2 partitions (dead-end vs allosteric):
        # dead-end → 17 deduped, allosteric → 4 TR-deduped
        # = 21 total
        result = collect(
            EnzymeRates.enumerate_mechanisms(
                uni_uni_reg_unknown))
        @test length(result) == 21
    end

    @testset "Uni-Bi + 1 unknown reg" begin
        # 2 partitions (dead-end vs allosteric):
        # dead-end → 580 deduped, allosteric → 88
        # TR-deduped = 668 total
        result = collect(
            EnzymeRates.enumerate_mechanisms(
                uni_bi_reg_unknown))
        @test length(result) == 668
    end
end

@testset "param_count accuracy" begin
    @testset "All Uni-Bi specs" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(uni_bi))
        @test length(all_specs) > 0
        for s in all_specs
            m = compile_mechanism(s)
            @test s.param_count ==
                length(parameters(m))
        end
    end

    @testset "Sampled Bi-Bi specs (unconstrained)" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(bi_bi))
        base = filter(
            s -> s isa EnzymeRates.MechanismSpec &&
                isempty(s.param_constraints),
            all_specs)
        rng = Random.MersenneTwister(42)
        n = min(20, length(base))
        sample =
            base[randperm(rng, length(base))[1:n]]
        for s in sample
            m = compile_mechanism(s)
            @test s.param_count ==
                length(parameters(m))
        end
    end

    # AllostericMechanismSpec.base.param_count only covers
    # catalytic parameters, not total allosteric params.
    # Accuracy check is limited to verifying compilation
    # succeeds and produces parameters.
    @testset "Allosteric specs (compilation)" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(
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

    # BUG: param_count formula doesn't account for
    # parameter equivalence constraints making
    # Wegscheider constraints redundant.
    # 36 of 126 constrained Bi-Bi MechanismSpecs
    # have param_count off by 1.
    @testset "Bi-Bi constrained param_count" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(bi_bi))
        constrained = filter(
            s -> s isa EnzymeRates.MechanismSpec &&
                !isempty(s.param_constraints),
            all_specs)
        @test length(constrained) == 126
        n_match = count(constrained) do s
            m = compile_mechanism(s)
            s.param_count == length(parameters(m))
        end
        # 90 of 126 match; 36 are off by 1
        @test n_match == 90
        @test_broken n_match == 126
    end
end

end # outer testset
