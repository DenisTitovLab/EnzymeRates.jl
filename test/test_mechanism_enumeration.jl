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
end
