@testset "Mechanism Enumeration" begin
    @testset "Uni-Uni default max_sites" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        forms = enumerate_enzyme_forms(r)
        names = Set(f.name for f in forms)
        # Sites: S₁, P₁ → 2² = 4 minus 1 excluded (E_S_P) = 3
        @test length(forms) == 3
        @test :E_0_0 ∈ names
        @test :E_S_0 ∈ names
        @test :E_0_P ∈ names
    end

    @testset "Uni-Uni with max_sites(S)=2" begin
        r = @enzyme_reaction begin
            substrates: S[C, 2]
            products:   P[C]
        end
        forms = enumerate_enzyme_forms(r)
        names = Set(f.name for f in forms)
        # Sites: S₁, P₁, S₂ → 2³ = 8 minus 3 excluded = 5
        @test length(forms) == 5
        @test :E_0_0_0 ∈ names   # free enzyme
        @test :E_S_0_0 ∈ names   # S in core site
        @test :E_0_0_S ∈ names   # S in extra site
        @test :E_S_0_S ∈ names   # S in both sites
        @test :E_0_P_0 ∈ names   # P in core site
    end

    @testset "Bi-Bi default max_sites" begin
        r = @enzyme_reaction begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end
        forms = enumerate_enzyme_forms(r)
        names = Set(f.name for f in forms)
        # Sites: A₁, B₁, P₁, Q₁ → 2⁴ = 16 minus 5 excluded = 11
        @test length(forms) == 11
        @test :E_0_0_0_0 ∈ names  # free enzyme
        @test :E_A_0_0_0 ∈ names  # only A bound
        @test :E_0_0_P_Q ∈ names  # both products (but not all substrates)
    end

    @testset "Ping-Pong Bi-Bi" begin
        r = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products:   P[C], Q[NX]
        end
        forms = enumerate_enzyme_forms(r)
        names = Set(f.name for f in forms)
        # Sites: A₁, B₁, P₁, Q₁ → 16 standard + 8 ping-pong - 7 excluded = 17
        @test :E_0_0_0_0 ∈ names     # free enzyme
        @test :E_A_0_0_0 ∈ names     # A bound
        @test :E_X_0_0_0 ∈ names     # ping-pong intermediate (residual X in A₁)
        @test :E_X_B_0_0 ∈ names     # intermediate with B bound
        @test :E_X_0_P_0 ∈ names     # intermediate with P bound
        @test :E_X_0_0_Q ∈ names     # intermediate with Q bound
        @test length(forms) == 17
    end

    @testset "Regulators add sites" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            regulators: I[N]
        end
        forms = enumerate_enzyme_forms(r)
        names = Set(f.name for f in forms)
        # Sites: S₁, P₁, I₁ → 2³ = 8 minus 2 excluded (E_S_P_0, E_S_P_I) = 6
        @test length(forms) == 6
        @test :E_0_0_0 ∈ names
        @test :E_S_0_I ∈ names   # S + inhibitor
        @test :E_0_P_I ∈ names   # P + inhibitor
    end

    @testset "No-atom species" begin
        r = @enzyme_reaction begin
            substrates: S
            products:   P
        end
        forms = enumerate_enzyme_forms(r)
        names = Set(f.name for f in forms)
        # Sites: S₁, P₁ → 4 minus 1 excluded (E_S_P) = 3
        @test length(forms) == 3
        @test :E_0_0 ∈ names
        @test :E_S_0 ∈ names
        @test :E_0_P ∈ names
    end

    @testset "Backward compat: 2-element tuples" begin
        # Direct constructor with 2-element tuples should auto-normalize
        r = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
        )
        forms = enumerate_enzyme_forms(r)
        @test length(forms) == 3
    end

    @testset "DSL produces 3-element tuples" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        subs = EnzymeRates.substrates(r)
        @test length(subs[1]) == 3
        @test subs[1][3] == 1  # default max_sites

        r2 = @enzyme_reaction begin
            substrates: S[C, 3]
            products:   P[C]
        end
        subs2 = EnzymeRates.substrates(r2)
        @test subs2[1][3] == 3
    end

    @testset "Site ordering" begin
        r = @enzyme_reaction begin
            substrates: A[C, 2], B[N]
            products:   P[C], Q[N]
        end
        forms = enumerate_enzyme_forms(r)
        # Site order: A₁(core), B₁(core), P₁(core), Q₁(core), A₂(extra)
        f = first(f for f in forms if f.name == :E_0_0_0_0_0)
        @test length(f.sites) == 5
        @test f.sites[1].metabolite == :A && f.sites[1].index == 1
        @test f.sites[2].metabolite == :B && f.sites[2].index == 1
        @test f.sites[3].metabolite == :P && f.sites[3].index == 1
        @test f.sites[4].metabolite == :Q && f.sites[4].index == 1
        @test f.sites[5].metabolite == :A && f.sites[5].index == 2
    end

    @testset "Multi-product ping-pong residual" begin
        # A[C2X] with products P1[C] and P2[C]: removing both gives residual {X}
        r = @enzyme_reaction begin
            substrates: A[C2X], B[N2]
            products:   P1[C], P2[C], Q[N2X]
        end
        forms = enumerate_enzyme_forms(r)
        names = Set(f.name for f in forms)
        # A[C2X] - P1[C] = {C,X} → valid residual
        # A[C2X] - P2[C] = {C,X} → same residual (dedup)
        # A[C2X] - P1[C] - P2[C] = {X} → another valid residual
        @test :E_CX_0_0_0_0 ∈ names  # residual {C,X} after one P release
        @test :E_X_0_0_0_0 ∈ names   # residual {X} after both P releases
    end

    @testset "Compilation time: 7-site reaction" begin
        # Guards against specialization explosion (e.g., Iterators.product splat
        # creating 2^n tuple types). Runs in a fresh subprocess to measure cold
        # compilation, then warm run in the same process.
        script = """
        using EnzymeRates
        rxn = @enzyme_reaction begin
            substrates: Glu[C6H12O6], ATP[C10H16N5O13P3]
            products: G6P[C6H13O9P], ADP[C10H15N5O10P2]
            regulators: Phosphate[PO4], G6P[C6H13O9P], G6P[C6H13O9P]
        end
        t1 = @elapsed enumerate_enzyme_forms(rxn)
        t2 = @elapsed enumerate_enzyme_forms(rxn)
        print(t1, " ", t2)
        """
        proj = dirname(@__DIR__)
        output = read(`$(Base.julia_cmd()) --project=$proj -e $script`, String)
        t_cold, t_warm = parse.(Float64, split(output))
        @test t_cold < 1.0
        @test t_warm < 0.1
    end

    # ─── enumerate_mechanisms tests ──────────────────────────────────────────

    @testset "Uni-Uni topologies" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        iter = enumerate_mechanisms(r)
        mechs = collect(iter)

        # Only one topology: E→ES→EP→E (3 steps: bind S, isomerize, release P)
        # RE/SS: 2³ - 1 = 7 valid assignments (at least one SS)
        @test length(mechs) == 7

        # Verify all convert to EnzymeMechanism
        for spec in mechs
            m = EnzymeMechanism(spec)
            @test m isa EnzymeMechanism
        end
    end

    @testset "Uni-Uni length matches iteration" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        iter = enumerate_mechanisms(r)
        @test length(iter) == length(collect(iter))
    end

    @testset "Bi-Bi topologies" begin
        r = @enzyme_reaction begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end
        iter = enumerate_mechanisms(r; max_forms=5)
        mechs = collect(iter)

        # Should include ordered and random-order Bi-Bi
        @test length(mechs) > 0

        # All mechanisms should convert without error
        for spec in mechs
            m = EnzymeMechanism(spec)
            @test m isa EnzymeMechanism
        end
    end

    @testset "Bi-Bi includes ordered and random-order" begin
        r = @enzyme_reaction begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end
        iter = enumerate_mechanisms(r; max_forms=6)
        mechs = collect(iter)

        # Check for ordered Bi-Bi: 5 forms (E, EA, EAB, EPQ, EQ)
        has_ordered = any(mechs) do spec
            length(spec.forms) == 5 && !any(spec.equilibrium_steps) # all SS
        end
        @test has_ordered

        # Check for random-order: 6+ forms (branched binding)
        has_random = any(mechs) do spec
            length(spec.forms) >= 6 && !any(spec.equilibrium_steps)
        end
        @test has_random
    end

    @testset "Regulator dead-end inhibitor" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            regulators: I[N]
        end
        iter = enumerate_mechanisms(r; max_forms=6)
        mechs = collect(iter)

        # Should include mechanisms with dead-end inhibitor complexes
        has_dead_end = any(mechs) do spec
            any(f -> occursin("I", string(f)), spec.forms)
        end
        @test has_dead_end

        # Should include EI (inhibitor bound to free enzyme)
        has_ei = any(mechs) do spec
            any(f -> f == :E_0_0_I, spec.forms)
        end
        @test has_ei

        # All convert without error
        for spec in mechs
            m = EnzymeMechanism(spec)
            @test m isa EnzymeMechanism
        end
    end

    @testset "Activator mechanism" begin
        # Regulators can participate in the catalytic cycle (essential activator)
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            regulators: A[N]
        end
        iter = enumerate_mechanisms(r; max_forms=8)
        mechs = collect(iter)

        # Should include activator cycle: E→EA→EAS→EAP→EA→... or similar
        # where A is bound during catalysis
        has_activator_cycle = any(mechs) do spec
            # Forms include both EA (activator bound) and EAS (activator + substrate)
            form_set = Set(spec.forms)
            :E_0_0_A ∈ form_set && :E_S_0_A ∈ form_set
        end
        @test has_activator_cycle

        # Also includes dead-end inhibitor variants (same regulator as inhibitor)
        has_inhibitor = any(mechs) do spec
            # Has inhibitor complex but NOT as part of the activator cycle
            form_set = Set(spec.forms)
            :E_0_0_A ∈ form_set && :E_S_0_A ∉ form_set
        end
        @test has_inhibitor
    end

    @testset "Product inhibition / abortive complexes" begin
        r = @enzyme_reaction begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end
        iter = enumerate_mechanisms(r; max_forms=8)
        mechs = collect(iter)

        # Should include dead-end product inhibition (P binding to E when
        # the cycle normally only has P released from EPQ)
        has_product_inhibition = any(mechs) do spec
            form_set = Set(spec.forms)
            :E_0_0_P_0 ∈ form_set || :E_0_0_0_Q ∈ form_set
        end
        @test has_product_inhibition

        # Should include abortive complexes (substrate + product simultaneously bound)
        has_abortive = any(mechs) do spec
            form_set = Set(spec.forms)
            :E_A_0_0_Q ∈ form_set || :E_0_B_P_0 ∈ form_set
        end
        @test has_abortive
    end

    @testset "All Uni-Uni mechanisms compile rate equations" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        for spec in enumerate_mechanisms(r)
            m = EnzymeMechanism(spec)
            s = rate_equation_string(m)
            @test s isa String
            @test !isempty(s)
        end
    end

    @testset "Stoichiometry: random-order is valid" begin
        r = @enzyme_reaction begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end
        # All generated mechanisms should pass constructor validation,
        # including random-order Bi-Bi which has futile cycles
        iter = enumerate_mechanisms(r; max_forms=6)
        for spec in iter
            @test EnzymeMechanism(spec) isa EnzymeMechanism
        end
    end

    @testset "max_forms limit" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            regulators: I[N]
        end
        small = collect(enumerate_mechanisms(r; max_forms=4))
        large = collect(enumerate_mechanisms(r; max_forms=8))
        @test length(small) <= length(large)

        # All mechanisms should respect the form limit
        for spec in small
            @test length(spec.forms) <= 4
        end
        for spec in large
            @test length(spec.forms) <= 8
        end
    end

    @testset "Equivalent step constraints" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            regulators: I[N]
        end
        iter = enumerate_mechanisms(r; max_forms=6)
        mechs = collect(iter)

        # Some mechanisms should have parameter constraints (equivalent step detection)
        has_constrained = any(spec -> !isempty(spec.param_constraints), mechs)
        has_unconstrained = any(spec -> isempty(spec.param_constraints), mechs)

        # Should have both constrained and unconstrained variants
        @test has_constrained
        @test has_unconstrained

        # Constrained mechanisms should still convert
        for spec in mechs
            if !isempty(spec.param_constraints)
                m = EnzymeMechanism(spec)
                @test m isa EnzymeMechanism
                break
            end
        end
    end

    @testset "Dead-end multi-child correctness" begin
        # Verify that dead-end nodes with multiple children generate all
        # valid downward-closed subsets (regression test for sibling bug)
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            regulators: I[N]
        end
        forms = enumerate_enzyme_forms(r)
        edges = EnzymeRates._build_reaction_graph(forms, r)
        adj = EnzymeRates._build_adjacency(forms, edges)
        cycles = EnzymeRates._find_valid_cycles(forms, edges, adj, r)
        topos = EnzymeRates._combine_cycles(cycles, forms, edges, 18, r)

        # Find a topology where a dead-end node has multiple children
        found_multi = false
        for topo in topos
            trees = EnzymeRates._build_dead_end_trees(topo, forms, adj, edges)
            for form_trees in trees
                for tree in form_trees
                    if length(tree.children) > 1
                        # This tree root has multiple children
                        subsets = EnzymeRates._dc_subsets(tree)
                        # Must include: empty, root only, root+child1, root+child2, root+both
                        @test length(subsets) >= 4  # at least: {}, {root}, {root,c1}, {root,c2}
                        # Check that root+both-children is present
                        root_idx = tree.form_idx
                        c1_idx = tree.children[1].form_idx
                        c2_idx = tree.children[2].form_idx
                        has_both = any(s -> root_idx ∈ s && c1_idx ∈ s && c2_idx ∈ s, subsets)
                        @test has_both
                        found_multi = true
                    end
                end
            end
        end
        @test found_multi  # confirm we actually tested a multi-child case
    end

    @testset "n_sites helper" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        @test n_sites(r) == 2

        r2 = @enzyme_reaction begin
            substrates: A[C, 2], B[N]
            products:   P[C], Q[N]
            regulators: I[N2, 3]
        end
        @test n_sites(r2) == 2 + 1 + 1 + 1 + 3  # A(2) + B(1) + P(1) + Q(1) + I(3)
    end
end
