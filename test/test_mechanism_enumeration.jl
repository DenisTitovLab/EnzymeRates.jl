@testset "Mechanism Enumeration" begin
    @testset "Uni-Uni default max_sites" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        forms = enumerate_enzyme_forms(r)
        names = Set(f.name for f in forms)
        # Sites: S₁, P₁ → 2² = 4 standard forms, no ping-pong (same atoms)
        @test length(forms) == 4
        @test :E_0_0 ∈ names
        @test :E_S_0 ∈ names
        @test :E_0_P ∈ names
        @test :E_S_P ∈ names
    end

    @testset "Uni-Uni with max_sites(S)=2" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            max_sites:  S => 2
        end
        forms = enumerate_enzyme_forms(r)
        names = Set(f.name for f in forms)
        # Sites: S₁, P₁, S₂ → 2³ = 8 standard forms
        @test length(forms) == 8
        @test :E_0_0_0 ∈ names   # free enzyme
        @test :E_S_0_0 ∈ names   # S in core site
        @test :E_0_P_0 ∈ names   # P in core site
        @test :E_0_0_S ∈ names   # S in extra site
        @test :E_S_P_0 ∈ names   # S + P in core sites
        @test :E_S_0_S ∈ names   # S in both sites
        @test :E_0_P_S ∈ names   # S in both sites
        @test :E_S_P_S ∈ names   # all occupied
    end

    @testset "Bi-Bi default max_sites" begin
        r = @enzyme_reaction begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end
        forms = enumerate_enzyme_forms(r)
        names = Set(f.name for f in forms)
        # Sites: A₁, B₁, P₁, Q₁ → 2⁴ = 16 standard forms
        # No ping-pong: A[C]→P[C] leaves empty residual, A[C]→Q[N] atoms don't match
        @test length(forms) == 16
        @test :E_0_0_0_0 ∈ names  # free enzyme
        @test :E_A_B_P_Q ∈ names  # all occupied
        @test :E_A_0_0_0 ∈ names  # only A bound
    end

    @testset "Ping-Pong Bi-Bi" begin
        r = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products:   P[C], Q[NX]
        end
        forms = enumerate_enzyme_forms(r)
        names = Set(f.name for f in forms)
        # Sites: A₁, B₁, P₁, Q₁ → 16 standard + ping-pong intermediates
        # A[CX] - P[C] = {X} → valid residual for A₁
        @test :E_0_0_0_0 ∈ names     # free enzyme
        @test :E_A_0_0_0 ∈ names     # A bound
        @test :E_X_0_0_0 ∈ names     # ping-pong intermediate (residual X in A₁)
        @test :E_X_B_0_0 ∈ names     # intermediate with B bound
        @test :E_X_0_P_0 ∈ names     # intermediate with P bound
        @test :E_X_0_0_Q ∈ names     # intermediate with Q bound

        # Count: 16 standard + 8 ping-pong (A₁={X}, others 2³)
        @test length(forms) == 24
    end

    @testset "max_ping_pong_intermediates=0 disables ping-pong" begin
        r = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products:   P[C], Q[NX]
        end
        forms_no_pp = enumerate_enzyme_forms(r; max_ping_pong_intermediates=0)
        names_no_pp = Set(f.name for f in forms_no_pp)
        # Without ping-pong: 4 sites → 2⁴ = 16 standard forms only
        @test length(forms_no_pp) == 16
        @test :E_X_0_0_0 ∉ names_no_pp  # no ping-pong intermediates

        # With ping-pong (default): should have more forms
        forms_with_pp = enumerate_enzyme_forms(r)
        @test length(forms_with_pp) > 16
        @test :E_X_0_0_0 ∈ Set(f.name for f in forms_with_pp)
    end

    @testset "max_total_bound limits forms" begin
        r = @enzyme_reaction begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end
        # max_total_bound=1 means at most 1 site occupied
        forms = enumerate_enzyme_forms(r; max_total_bound=1)
        names = Set(f.name for f in forms)
        # Free enzyme + 4 single-bound forms = 5
        @test length(forms) == 5
        @test :E_0_0_0_0 ∈ names
        @test :E_A_0_0_0 ∈ names
        @test :E_0_B_0_0 ∈ names
        @test :E_0_0_P_0 ∈ names
        @test :E_0_0_0_Q ∈ names
    end

    @testset "max_total_bound=0 gives free enzyme only" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        forms = enumerate_enzyme_forms(r; max_total_bound=0)
        @test length(forms) == 1
        @test forms[1].name == :E_0_0
    end

    @testset "Regulators add sites" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            regulators: I[N]
        end
        forms = enumerate_enzyme_forms(r)
        names = Set(f.name for f in forms)
        # Sites: S₁, P₁, I₁ → 2³ = 8
        @test length(forms) == 8
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
        # Sites: S₁, P₁ → 4 standard, no ping-pong (no atoms)
        @test length(forms) == 4
        @test :E_0_0 ∈ names
        @test :E_S_0 ∈ names
        @test :E_0_P ∈ names
        @test :E_S_P ∈ names
    end

    @testset "SiteState and total_atoms" begin
        r = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products:   P[C], Q[NX]
        end
        forms = enumerate_enzyme_forms(r)
        # Find the ping-pong intermediate E_X_0_0_0
        f = first(f for f in forms if f.name == :E_X_0_0_0)
        @test total_atoms(f) == [:X => 1]

        # Find full A bound
        f_a = first(f for f in forms if f.name == :E_A_0_0_0)
        @test total_atoms(f_a) == [:C => 1, :X => 1]

        # Free enzyme has no atoms
        f_free = first(f for f in forms if f.name == :E_0_0_0_0)
        @test total_atoms(f_free) == Pair{Symbol,Int}[]
    end

    @testset "max_binding_sites accessor" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            max_sites:  S => 2
        end
        @test max_binding_sites(r, :S) == 2
        @test max_binding_sites(r, :P) == 1
        @test_throws ErrorException max_binding_sites(r, :X)
    end

    @testset "Backward compat: 2-element tuples" begin
        # Direct constructor with 2-element tuples should auto-normalize
        r = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
        )
        @test max_binding_sites(r, :S) == 1
        @test max_binding_sites(r, :P) == 1
        forms = enumerate_enzyme_forms(r)
        @test length(forms) == 4
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
            substrates: S[C]
            products:   P[C]
            max_sites:  S => 3
        end
        subs2 = EnzymeRates.substrates(r2)
        @test subs2[1][3] == 3
    end

    @testset "Site ordering" begin
        r = @enzyme_reaction begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
            max_sites:  A => 2
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

    @testset "show method" begin
        r = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        @test sprint(show, r) == "EnzymeReaction: S ⇌ P"  # no max_sites shown (all default)

        r2 = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            max_sites:  S => 2
        end
        s = sprint(show, r2)
        @test contains(s, "max_sites: S=2")
    end
end
