@testset "Types" begin
    @testset "EnzymeMechanism struct + accessors (new design)" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end

        @test EnzymeRates.substrates(m) == (:S,)
        @test EnzymeRates.products(m) == (:P,)
        @test EnzymeRates.regulators(m) == ()
        @test EnzymeRates.metabolites(m) == (:S, :P)
        @test EnzymeRates.reactions(m) == (
            ((:E, :S), (:ES,), true,  1),
            ((:ES,),   (:EP,), false, 2),
            ((:E, :P), (:EP,), true,  3),
        )
        @test EnzymeRates.equilibrium_steps(m) == (true, false, true)
        @test EnzymeRates.n_steps(m) == 3
        @test EnzymeRates.kinetic_group(m, 1) == 1
        @test EnzymeRates.kinetic_group(m, 2) == 2
        @test EnzymeRates.kinetic_group(m, 3) == 3
        @test EnzymeRates.kinetic_groups(m) == (1, 2, 3)
        @test EnzymeRates.steps_in_group(m, 1) == (1,)
        @test EnzymeRates.enzyme_forms(m) == (:E, :ES, :EP)
        @test EnzymeRates.n_states(m) == 3

        # Shared kinetic-group: two steps in group 1 (regulator R binds
        # both E and E(S) sharing one K).
        m2 = @enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: R
            steps: begin
                (E + R ⇌ E(R), E(S) + R ⇌ E(S, R))
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        @test EnzymeRates.kinetic_group(m2, 1) == 1
        @test EnzymeRates.kinetic_group(m2, 2) == 1
        @test EnzymeRates.steps_in_group(m2, 1) == (1, 2)
    end

    @testset "metabolites() lift covers all three loops + dedup" begin
        # Substrate, product, AND regulator so the lift exercises every loop
        # (substrates, then products, then regulators) and the dedup `seen`
        # set — coverage the plain S/P accessor tests do not reach.
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: I
            steps: begin
                E + S ⇌ E(S)
                E(S) + I ⇌ E(S, I)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        @test EnzymeRates.metabolites(m) == (:S, :P, :I)
    end

    @testset "EnzymeMechanism Sig repack" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end

        @test m isa EnzymeRates.EnzymeMechanism
        Sig = typeof(m).parameters[1]
        @test Sig isa Tuple
        @test length(Sig) == 2

        @test EnzymeRates.substrates(m) == (:S,)
        @test EnzymeRates.products(m)   == (:P,)
        @test EnzymeRates.n_steps(m)    == 3
    end

    @testset "_to_sig metabolite encoding round-trips" begin
        for M in (EnzymeRates.Substrate, EnzymeRates.Product,
                  EnzymeRates.AllostericRegulator, EnzymeRates.CompetitiveInhibitor)
            sig = EnzymeRates._to_sig(M(:X))
            @test sig == (nameof(M), :X)
            @test EnzymeRates._metabolite_from_sig(sig) == M(:X)
        end
    end

    @testset "EnzymeMechanism constructor" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        @test EnzymeRates.n_steps(m) == 3
        @test EnzymeRates.substrates(m) == (:S,)

        # Same-kinetics group test: regulator R binds E and E(S) sharing one K
        m_g = @enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: R
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                (E + R ⇌ E(R), E(S) + R ⇌ E(S, R))
            end
        end
        @test EnzymeRates.kinetic_group(m_g, 4) == EnzymeRates.kinetic_group(m_g, 5)
    end

    @testset "AllostericEnzymeMechanism (new design)" begin
        cm = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        # Dense format: (multiplicity, cat_allo_states) and (ligands, mult, reg_allo_states)
        cat_sites = (2, (:NonequalAI, :OnlyA, :NonequalAI))
        reg_sites = ((((:I,), 2, (:OnlyI,)),),)
        m = EnzymeRates.AllostericEnzymeMechanism{typeof(cm), cat_sites, reg_sites[1]}()

        @test EnzymeRates.catalytic_mechanism(m) === cm
        @test EnzymeRates.catalytic_multiplicity(m) == 2
        @test EnzymeRates.cat_allo_state(m, 1) == :NonequalAI
        @test EnzymeRates.cat_allo_state(m, 2) == :OnlyA
        @test EnzymeRates.regulatory_sites(m) == reg_sites[1]
    end

    @testset "AllostericEnzymeMechanism constructor + DSL" begin
        cm = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end

        # Single-ligand :EqualAI reg site → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalAI, :NonequalAI, :NonequalAI)),
            (((:I,), 2, (:EqualAI,)),),
        )

        # Catalytic group :OnlyI → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalAI, :OnlyI, :NonequalAI)), (),
        )

        # Build via DSL
        m = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: I::OnlyI

            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)     :: EqualAI
                E(S) <--> E(P)   :: OnlyA
                E(P) ⇌ E + P     :: EqualAI
            end
        end
        @test EnzymeRates.catalytic_multiplicity(m) == 2
        @test EnzymeRates.cat_allo_state(m, 1) == :EqualAI
        @test EnzymeRates.cat_allo_state(m, 2) == :OnlyA
        @test EnzymeRates.allosteric_regulators(m) == ((:I, :OnlyI),)
    end

    @testset "Pretty printing" begin
        # Linear mechanism: compact chain.
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
            end
        end
        @test sprint(show, m) ==
            "EnzymeMechanism: E + S <--> ES <--> E + P"

        # Branched mechanism: multi-line with header summary.
        m_b = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> E(A)
                E + B <--> E(B)
                E(A) + B <--> E(A, B)
                E(B) + A <--> E(A, B)
                E(A, B) <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end
        s = sprint(show, m_b)
        @test startswith(s, "EnzymeMechanism (7 steps, 6 enzyme forms):")
        @test contains(s, "E + A <--> EA")
        @test contains(s, "E + Q <--> EQ")

        # Linear chain with canonical RE binding: chain-walk renders the
        # release step in reverse so the chain stays linear.
        m_re = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        @test sprint(show, m_re) ==
            "EnzymeMechanism: E + S ⇌ ES <--> EP ⇌ E + P"

        # Mechanism with regulators: appended at end.
        m_reg = @enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: I
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
                E + I <--> E(I)
            end
        end
        @test contains(sprint(show, m_reg), "| regulators: I")

        # EnzymeReaction with oligomeric_state > 1.
        rxn_oligo = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            oligomeric_state: 4
        end
        @test sprint(show, rxn_oligo) ==
            "EnzymeReaction: S ⇌ P | oligomeric_state: 4"

        # AllostericEnzymeMechanism (smoke).
        m_allo = @allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            allosteric_regulators: I::OnlyI
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + F6P ⇌ E(F6P)         :: EqualAI
                E(F6P) <--> E(F16BP)     :: EqualAI
                E(F16BP) ⇌ E + F16BP     :: EqualAI
            end
        end
        s_allo = sprint(show, m_allo)
        @test contains(s_allo, "AllostericEnzymeMechanism (cat_n=2")
        @test contains(s_allo, "reg sites")
        @test contains(s_allo, "I::OnlyI")
    end

    @testset "EnzymeMechanism different orderings produce valid mechanisms" begin
        m1 = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
            end
        end
        m2 = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E(S) <--> E + P
                E + S <--> E(S)
            end
        end
        @test m1 isa EnzymeMechanism
        @test m2 isa EnzymeMechanism
    end

    @testset "EnzymeMechanism error cases" begin
        # Empty steps → error (re-pointed to _assert_mechanism_invariants:
        # the decomposed Mechanism with no steps errors on `isempty(flat)`).
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        m_empty = EnzymeRates.Mechanism(rxn, Vector{Vector{EnzymeRates.Step}}())
        @test_throws ErrorException EnzymeRates._assert_mechanism_invariants(m_empty)

        # Net stoichiometry mismatch: a declared substrate that no step binds.
        # Re-pointed to _assert_mechanism_invariants (substrate coverage).
        rxn_unused = @enzyme_reaction begin
            substrates: S[C], T[N]
            products:   P[CN]
        end
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)
        m_unused = EnzymeRates.Mechanism(rxn_unused, [
            [EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)],
            [EnzymeRates.Step(e_s, e_p, nothing, false)],
            [EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)],
        ])
        @test_throws ErrorException EnzymeRates._assert_mechanism_invariants(m_unused)

        # NOTE: "Duplicate reactions" and "Unreachable enzyme form" tests
        # were dropped — the new design accepts both. Two reactions with
        # the same (lhs, rhs) but distinct kinetic_groups are valid (they
        # represent dead-end mirrors with different parameters in the OLD
        # design's terms; in the new design, the kinetic_group integer
        # disambiguates). And the new constructor doesn't enforce a
        # connectivity invariant — enzyme forms are inferred from steps,
        # so an "unreachable" form simply has its own steps in isolation,
        # which is structurally valid (graph connectivity is a downstream
        # concern caught by Wegscheider analysis if it matters).
    end

    @testset "EnzymeMechanism valid with reachable enzyme forms" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
            end
        end
        @test m isa EnzymeMechanism
    end

    @testset "Kinetic-group validator error paths" begin
        # Group binding different metabolites → error (re-pointed to
        # _assert_mechanism_invariants over a hand-built decomposed Mechanism).
        rxn_two = @enzyme_reaction begin
            substrates: S[C], A[N]
            products:   P[CN]
        end
        e    = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s  = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_a  = EnzymeRates.Species([EnzymeRates.Substrate(:A)], :E)
        e_p2 = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)
        g1_s = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        g1_a = EnzymeRates.Step(e, e_a, EnzymeRates.Substrate(:A), true)
        g2_iso = EnzymeRates.Step(e_s, e_p2, nothing, false)
        m_diffmet = EnzymeRates.Mechanism(rxn_two, [[g1_s, g1_a], [g2_iso]])
        @test_throws ErrorException EnzymeRates._assert_mechanism_invariants(m_diffmet)

        # Group mixing RE and SS → error (re-pointed). Same metabolite, one
        # RE binding step and one SS binding step share a kinetic group.
        rxn_uni = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        s_re = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        s_ss = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), false)
        s_rel = EnzymeRates.Step(e, e_p2, EnzymeRates.Product(:P), true)
        m_mix = EnzymeRates.Mechanism(rxn_uni, [[s_re, s_ss], [s_rel]])
        @test_throws ErrorException EnzymeRates._assert_mechanism_invariants(m_mix)
    end

    @testset "AllostericEnzymeMechanism constructor validators" begin
        cm = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        # Wrong-length cat_allo_states (4 entries for 3 kinetic groups) → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalAI, :NonequalAI, :NonequalAI, :OnlyA)), ())

        # Invalid allo state value → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NotAState, :NonequalAI, :NonequalAI)), ())

        # Reg site with no ligands → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalAI, :NonequalAI, :NonequalAI)), (((), 2, ()),))

        # Reg site with all-:EqualAI ligands → error (cancels identically)
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalAI, :NonequalAI, :NonequalAI)),
            (((:I, :J), 2, (:EqualAI, :EqualAI)),))

        # Invalid reg-site ligand allo state → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalAI, :NonequalAI, :NonequalAI)),
            (((:I,), 2, (:NotAState,)),))
    end

    @testset "Base.show displays all dense states" begin
        m = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: I::NonequalAI, J::OnlyI
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)     :: NonequalAI
                E(S) <--> E(P)   :: OnlyA
                E(P) ⇌ E + P     :: EqualAI
            end
        end
        s = sprint(show, m)
        # Every catalytic state appears in cat_allo_states line
        cm_inner = EnzymeRates.catalytic_mechanism(m)
        n_groups = length(unique(EnzymeRates.kinetic_group(cm_inner, i)
                                 for i in 1:EnzymeRates.n_steps(cm_inner)))
        for g in 1:n_groups
            @test occursin(string(EnzymeRates.cat_allo_state(m, g)), s)
        end
        # No :NonequalAI ligand silently hidden from reg-site display
        for (i, _) in enumerate(EnzymeRates.regulatory_sites(m))
            for lig in EnzymeRates.regulatory_sites(m)[i][1]
                state = EnzymeRates.reg_allo_state(m, i, lig)
                @test occursin("$lig::$state", s)
            end
        end
    end

    @testset "EnzymeMechanism: regulator binding" begin
        # All regulators bound -> ok
        @test (@enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: A
            steps: begin
                E + S ⇌ E(S)
                E(S) + A ⇌ E(S, A)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end) isa EnzymeRates.EnzymeMechanism
        # No regulators -> ok
        @test (@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end) isa EnzymeRates.EnzymeMechanism
    end

    @testset "AllostericEnzymeMechanism display format" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        base = first(EnzymeRates.init_mechanisms(rxn))
        cat_allo_states = Symbol[]
        for g in EnzymeRates.kinetic_groups(base)
            rep = EnzymeRates.rep_step(base, g)
            met = EnzymeRates.bound_metabolite(rep)
            tag = (met isa EnzymeRates.Reactant &&
                   EnzymeRates.name(met) in (:S, :P)) ?
                  :EqualAI : :NonequalAI
            push!(cat_allo_states, tag)
        end
        site = EnzymeRates.RegulatorySite(
            [EnzymeRates.AllostericRegulator(:R)], 2, [:OnlyI])
        am = EnzymeRates.AllostericMechanism(
            EnzymeRates.reaction(base), copy(EnzymeRates.steps(base)),
            cat_allo_states, 2, [site])
        m = EnzymeRates.AllostericEnzymeMechanism(am)
        s = repr(m)

        # Old summary line gone:
        @test !occursin("cat_allo_states:", s)
        # Inline ::Tag annotations on each step or step group:
        @test occursin(":: EqualAI", s)
        # Multi-line catalytic display (no chain shortcut):
        n_steps_re = count(c -> c == '\n', s)
        @test n_steps_re >= 3   # header + ≥3 step lines
    end

    @testset "AllostericEnzymeMechanism display: shared kinetic group" begin
        cm = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                (E + S ⇌ E(S), E(P) + S ⇌ E(P, S))
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        am = EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:EqualAI, :EqualAI, :EqualAI)), ())
        s = repr(am)

        # The parenthesized-group branch must execute: look for the
        # exact "(...) :: EqualAI" shape with a comma-separated body.
        paren_group_match = match(
            r"\([^()]*,[^()]*\) :: EqualAI", s)
        @test paren_group_match !== nothing
        # Format: lhs/rhs joined with " + " (no brackets).
        @test occursin("ES <--> EP :: EqualAI", s)
        @test occursin(":: EqualAI", s)
        @test !occursin("cat_allo_states:", s)
    end

    # ─── Concrete type hierarchy ──────────────────────────────────────

    @testset "Metabolite hierarchy: Substrate / Product / Regulators" begin
        s = EnzymeRates.Substrate(:ATP)
        p = EnzymeRates.Product(:ADP)
        a = EnzymeRates.AllostericRegulator(:cAMP)
        i = EnzymeRates.CompetitiveInhibitor(:I)

        @test s isa EnzymeRates.Substrate
        @test s isa EnzymeRates.Reactant
        @test s isa EnzymeRates.Metabolite
        @test p isa EnzymeRates.Product
        @test p isa EnzymeRates.Reactant
        @test a isa EnzymeRates.AllostericRegulator
        @test a isa EnzymeRates.Regulator
        @test a isa EnzymeRates.Metabolite
        @test i isa EnzymeRates.CompetitiveInhibitor
        @test i isa EnzymeRates.Regulator

        @test EnzymeRates.name(s) === :ATP
        @test EnzymeRates.name(p) === :ADP
        @test EnzymeRates.name(a) === :cAMP
        @test EnzymeRates.name(i) === :I

        @test EnzymeRates.Substrate(:X) == EnzymeRates.Substrate(:X)
        @test EnzymeRates.Substrate(:X) != EnzymeRates.Substrate(:Y)
        # Distinct subtypes with same name are NOT equal (struct identity matters).
        @test EnzymeRates.Substrate(:X) != EnzymeRates.Product(:X)
        @test hash(EnzymeRates.Substrate(:X)) == hash(EnzymeRates.Substrate(:X))
    end

    @testset "Residual: empty default + canonical ordering" begin
        empty_r = EnzymeRates.Residual()
        @test isempty(empty_r)
        @test EnzymeRates.added(empty_r) == EnzymeRates.Substrate[]
        @test EnzymeRates.subtracted(empty_r) == EnzymeRates.Product[]

        r1 = EnzymeRates.Residual(
            [EnzymeRates.Substrate(:B), EnzymeRates.Substrate(:A)],
            [EnzymeRates.Product(:Q), EnzymeRates.Product(:P)],
        )
        r2 = EnzymeRates.Residual(
            [EnzymeRates.Substrate(:A), EnzymeRates.Substrate(:B)],
            [EnzymeRates.Product(:P), EnzymeRates.Product(:Q)],
        )
        @test r1 == r2
        @test hash(r1) == hash(r2)
        @test !isempty(r1)
        @test EnzymeRates.added(r1) ==
              [EnzymeRates.Substrate(:A), EnzymeRates.Substrate(:B)]
        @test EnzymeRates.subtracted(r1) ==
              [EnzymeRates.Product(:P), EnzymeRates.Product(:Q)]
    end

    @testset "Species: canonical bound ordering + accessors + name" begin
        s1 = EnzymeRates.Species(
            EnzymeRates.Metabolite[
                EnzymeRates.Substrate(:B), EnzymeRates.Substrate(:A)],
            :E,
        )
        s2 = EnzymeRates.Species(
            EnzymeRates.Metabolite[
                EnzymeRates.Substrate(:A), EnzymeRates.Substrate(:B)],
            :E,
        )
        @test s1 == s2
        @test hash(s1) == hash(s2)
        @test EnzymeRates.conformation(s1) === :E
        @test EnzymeRates.residual(s1) == EnzymeRates.Residual()
        @test !EnzymeRates.has_residual(s1)
        @test EnzymeRates.bound(s1) ==
              EnzymeRates.Metabolite[
                  EnzymeRates.Substrate(:A), EnzymeRates.Substrate(:B)]

        # Empty bound, :E conformation → :E
        s_e = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        @test EnzymeRates.name(s_e) === :E

        # Bound metabolites are appended in canonical order
        s_es = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:A),
                                   EnzymeRates.Substrate(:B)],
            :E,
        )
        @test EnzymeRates.name(s_es) === :EAB

        # Estar conformation
        s_estar = EnzymeRates.Species(EnzymeRates.Metabolite[], :Estar)
        @test EnzymeRates.name(s_estar) === :Estar

        # Three-arg constructor exposes residual
        res = EnzymeRates.Residual(
            [EnzymeRates.Substrate(:A)],
            [EnzymeRates.Product(:P)],
        )
        s_res = EnzymeRates.Species(EnzymeRates.Metabolite[], :Estar, res)
        @test EnzymeRates.has_residual(s_res)
        @test EnzymeRates.residual(s_res) == res

        # Same metabolite name bound in two roles (product + competitive
        # inhibitor) canonicalizes regardless of construction order: same
        # Species, same hash, same form name (non-inhibitor segment first).
        d1 = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Product(:G6P),
                                   EnzymeRates.CompetitiveInhibitor(:G6P)], :E)
        d2 = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.CompetitiveInhibitor(:G6P),
                                   EnzymeRates.Product(:G6P)], :E)
        @test d1 == d2
        @test hash(d1) == hash(d2)
        @test EnzymeRates.name(d1) === EnzymeRates.name(d2)
        @test EnzymeRates.name(d1) === :EG6PG6Pinh
    end

    @testset "RegulatorySite: validation + accessors" begin
        lig_a = EnzymeRates.AllostericRegulator(:A)
        lig_b = EnzymeRates.AllostericRegulator(:B)
        site = EnzymeRates.RegulatorySite(
            [lig_a, lig_b], 4, [:OnlyA, :NonequalAI],
        )
        @test EnzymeRates.ligands(site) == [lig_a, lig_b]
        @test EnzymeRates.multiplicity(site) == 4
        @test EnzymeRates.allo_states(site) == [:OnlyA, :NonequalAI]

        # Mismatched ligand / allo_state length → error
        @test_throws ErrorException EnzymeRates.RegulatorySite(
            [lig_a, lig_b], 4, [:OnlyA])

        # Multiplicity < 1 → error
        @test_throws ErrorException EnzymeRates.RegulatorySite(
            [lig_a], 0, [:OnlyA])

        # Invalid allo state → error
        @test_throws ErrorException EnzymeRates.RegulatorySite(
            [lig_a], 1, [:NotAState])

        # All four allowed states accepted
        for st in (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)
            @test EnzymeRates.RegulatorySite([lig_a], 1, [st]) isa
                  EnzymeRates.RegulatorySite
        end

        # Equality / hash
        site2 = EnzymeRates.RegulatorySite(
            [lig_a, lig_b], 4, [:OnlyA, :NonequalAI])
        @test site == site2
        @test hash(site) == hash(site2)
        # Order-sensitive: ligand ordering is parallel to allo_states,
        # so [A,B] / [OnlyA,NonequalAI] != [B,A] / [OnlyA,NonequalAI].
        site_reordered = EnzymeRates.RegulatorySite(
            [lig_b, lig_a], 4, [:OnlyA, :NonequalAI])
        @test site != site_reordered
    end

    @testset "RegulatorySite: A/I allo-state symbols accepted; R/T rejected" begin
        for st in (:EqualAI, :OnlyA, :OnlyI, :NonequalAI)
            @test EnzymeRates.RegulatorySite(
                [EnzymeRates.AllostericRegulator(:G6P)], 1, [st]) isa EnzymeRates.RegulatorySite
        end
        for st in (:OnlyR, :OnlyT, :EqualRT, :NonequalRT)
            @test_throws ErrorException EnzymeRates.RegulatorySite(
                [EnzymeRates.AllostericRegulator(:G6P)], 1, [st])
        end
    end

    @testset "Step fields + accessors" begin
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Product(:P)], :E)

        s1 = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        s2 = EnzymeRates.Step(e_s, e_p, nothing, false)
        s3 = EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)

        @test fieldnames(EnzymeRates.Step) ==
              (:from_species, :to_species, :bound_metabolite,
               :is_equilibrium)

        @test EnzymeRates.from_species(s1) === e
        @test EnzymeRates.to_species(s1) === e_s
        @test EnzymeRates.bound_metabolite(s1) ==
              EnzymeRates.Substrate(:S)
        @test EnzymeRates.is_equilibrium(s1)
        @test EnzymeRates.is_binding(s1)
        @test !EnzymeRates.is_iso(s1)
        @test EnzymeRates.direction(s1) === :binding

        @test EnzymeRates.bound_metabolite(s2) === nothing
        @test EnzymeRates.is_iso(s2)
        @test !EnzymeRates.is_binding(s2)
        @test EnzymeRates.direction(s2) === :iso

        @test EnzymeRates.is_binding(s3)
    end

    @testset "Step `==` / hash are structural" begin
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)], :E)
        s  = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        s2 = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        @test s == s2
        @test hash(s) == hash(s2)
    end

    @testset "Step canonicalizes binding direction" begin
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)], :E)

        # User authored release direction (E_S → E + S, metabolite on RHS).
        # Constructor swaps to binding direction (E + S → E_S). Both RE and
        # SS binding canonicalize this way.
        re_released = EnzymeRates.Step(e_s, e, EnzymeRates.Substrate(:S), true)
        re_bound    = EnzymeRates.Step(e,   e_s, EnzymeRates.Substrate(:S), true)
        @test re_released == re_bound
        @test hash(re_released) == hash(re_bound)
        @test EnzymeRates.from_species(re_released) === e
        @test EnzymeRates.to_species(re_released) === e_s

        ss_released = EnzymeRates.Step(e_s, e, EnzymeRates.Substrate(:S), false)
        ss_bound    = EnzymeRates.Step(e,   e_s, EnzymeRates.Substrate(:S), false)
        @test ss_released == ss_bound
        @test hash(ss_released) == hash(ss_bound)
        @test EnzymeRates.from_species(ss_released) === e
        @test EnzymeRates.to_species(ss_released) === e_s
    end

    @testset "Step preserves iso direction (canonicalized in Mechanism ctor)" begin
        e_s = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Product(:P)], :E)

        # The Step constructor does NOT canonicalize iso steps (RE or SS) —
        # iso direction depends on the reaction's substrate/product sets and
        # is decided by `_canonical_iso_direction` in the Mechanism / Allosteric
        # Mechanism constructor. At the bare-Step level, direction is preserved.
        re_fwd = EnzymeRates.Step(e_s, e_p, nothing, true)
        re_rev = EnzymeRates.Step(e_p, e_s, nothing, true)
        @test re_fwd != re_rev
        @test EnzymeRates.from_species(re_fwd) === e_s
        @test EnzymeRates.from_species(re_rev) === e_p

        ss_fwd = EnzymeRates.Step(e_s, e_p, nothing, false)
        ss_rev = EnzymeRates.Step(e_p, e_s, nothing, false)
        @test ss_fwd != ss_rev
        @test EnzymeRates.from_species(ss_fwd) === e_s
        @test EnzymeRates.from_species(ss_rev) === e_p
    end

    @testset "Parameter family: step-bound, Kreg, mechanism-level" begin
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)], :E)
        step = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)

        kd_none = EnzymeRates.Kd(step, :None)
        kd_i    = EnzymeRates.Kd(step, :I)
        @test kd_none isa EnzymeRates.Kd
        @test kd_none isa EnzymeRates.Parameter
        @test EnzymeRates.governing_step(kd_none) === step
        @test !EnzymeRates.is_i_state(kd_none)
        @test EnzymeRates.is_i_state(kd_i)
        @test kd_none == EnzymeRates.Kd(step, :None)
        @test kd_none != kd_i

        for T in (EnzymeRates.Kiso, EnzymeRates.Kon, EnzymeRates.Koff,
                  EnzymeRates.Kfor, EnzymeRates.Krev)
            p = T(step, :None)
            @test p isa EnzymeRates.Parameter
            @test EnzymeRates.governing_step(p) === step
            @test !EnzymeRates.is_i_state(p)
            @test EnzymeRates.is_i_state(T(step, :I))
        end

        lig_a = EnzymeRates.AllostericRegulator(:A)
        site = EnzymeRates.RegulatorySite([lig_a], 2, [:OnlyA])
        kr = EnzymeRates.Kreg(site, lig_a, :A)
        @test kr isa EnzymeRates.Parameter
        @test EnzymeRates.is_i_state(EnzymeRates.Kreg(site, lig_a, :I))
        @test !EnzymeRates.is_i_state(kr)
        @test kr == EnzymeRates.Kreg(site, lig_a, :A)

        # Mechanism-level scalars: singletons
        @test EnzymeRates.Keq() == EnzymeRates.Keq()
        @test EnzymeRates.Etot() == EnzymeRates.Etot()
        @test EnzymeRates.Lallo() == EnzymeRates.Lallo()
        @test EnzymeRates.Keq() isa EnzymeRates.Parameter
        @test EnzymeRates.Etot() isa EnzymeRates.Parameter
        @test EnzymeRates.Lallo() isa EnzymeRates.Parameter
    end

    @testset "ReactantAtoms canonicalizes atom ordering" begin
        ra1 = EnzymeRates.ReactantAtoms(
            EnzymeRates.Substrate(:ATP),
            [:C => 10, :H => 16, :N => 5],
        )
        ra2 = EnzymeRates.ReactantAtoms(
            EnzymeRates.Substrate(:ATP),
            [:N => 5, :H => 16, :C => 10],
        )
        @test ra1 == ra2
        @test hash(ra1) == hash(ra2)
        @test EnzymeRates.metabolite(ra1) == EnzymeRates.Substrate(:ATP)
        @test EnzymeRates.atoms(ra1) == [:C => 10, :H => 16, :N => 5]

        # Distinct metabolite kinds: ATP-Substrate != ATP-Product
        ra_p = EnzymeRates.ReactantAtoms(
            EnzymeRates.Product(:ATP), [:C => 10, :H => 16, :N => 5])
        @test ra1 != ra_p
    end

    @testset "ReactantAtoms validation" begin
        # Mandatory atoms: empty atom list is rejected.
        @test_throws ErrorException EnzymeRates.ReactantAtoms(
            EnzymeRates.Substrate(:S), Pair{Symbol,Int}[])
        # Positive counts: zero / negative rejected.
        @test_throws ErrorException EnzymeRates.ReactantAtoms(
            EnzymeRates.Substrate(:S), [:C => 0])
        @test_throws ErrorException EnzymeRates.ReactantAtoms(
            EnzymeRates.Substrate(:S), [:C => -1])
        # Bool is not a valid count. Use untyped vector so Bool survives to
        # the constructor (Pair{Symbol,Int}[:C => true] would convert true→1).
        @test_throws ErrorException EnzymeRates.ReactantAtoms(
            EnzymeRates.Substrate(:S), [:C => true])
        # Valid construction still works.
        ra = EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 6, :H => 12])
        @test EnzymeRates.atoms(ra) == [:C => 6, :H => 12]
    end

    @testset "EnzymeReaction validation" begin
        S = EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1])
        P = EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])
        P2 = EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])
        Punbal = EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 2])
        noregs = EnzymeRates.RegulatorMults[]
        regS = EnzymeRates.RegulatorMults(EnzymeRates.CompetitiveInhibitor(:S), [1])
        regA = EnzymeRates.RegulatorMults(EnzymeRates.CompetitiveInhibitor(:A), [1])
        regA2 = EnzymeRates.RegulatorMults(EnzymeRates.CompetitiveInhibitor(:A), [1])

        # Empty substrate set rejected (only a product in reactants).
        @test_throws ErrorException EnzymeRates.EnzymeReaction([P], noregs, Int[1])
        # Empty product set rejected (only a substrate in reactants).
        @test_throws ErrorException EnzymeRates.EnzymeReaction([S], noregs, Int[1])
        # Duplicate product names rejected.
        @test_throws ErrorException EnzymeRates.EnzymeReaction([S, P, P2], noregs, Int[1])
        # Atom imbalance rejected (S has 1 C, Punbal has 2 C).
        @test_throws ErrorException EnzymeRates.EnzymeReaction([S, Punbal], noregs, Int[1])
        # Duplicate regulator names rejected.
        @test_throws ErrorException EnzymeRates.EnzymeReaction(
            [S, P], [regA, regA2], Int[1])
        # A regulator MAY share a substrate/product name (`::Inh` role tag:
        # one metabolite binds as a CompetitiveInhibitor under its real name).
        rxn_inh = EnzymeRates.EnzymeReaction([S, P], [regS], Int[1])
        @test :S in EnzymeRates.name.(EnzymeRates.substrates(rxn_inh))
        # Valid balanced reaction with a distinct-named regulator still constructs.
        rxn = EnzymeRates.EnzymeReaction([S, P], [regA], Int[1])
        @test EnzymeRates.name.(EnzymeRates.substrates(rxn)) == [:S]
        @test EnzymeRates.name.(EnzymeRates.products(rxn)) == [:P]
        # Balance dispatches by metabolite TYPE, not name: a substrate and a
        # product sharing a name (X→X, 2 C each) balances and is not misrouted.
        Xs = EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:X), [:C => 2])
        Xp = EnzymeRates.ReactantAtoms(EnzymeRates.Product(:X), [:C => 2])
        @test EnzymeRates.EnzymeReaction([Xs, Xp], noregs, Int[1]) isa
              EnzymeRates.EnzymeReaction
    end

    @testset "RegulatorMults canonicalizes ordering + validates" begin
        rm1 = EnzymeRates.RegulatorMults(
            EnzymeRates.AllostericRegulator(:A), [4, 1, 2])
        rm2 = EnzymeRates.RegulatorMults(
            EnzymeRates.AllostericRegulator(:A), [1, 2, 4])
        @test rm1 == rm2
        @test hash(rm1) == hash(rm2)
        @test EnzymeRates.regulator(rm1) ==
              EnzymeRates.AllostericRegulator(:A)
        @test EnzymeRates.allowed_multiplicities(rm1) == [1, 2, 4]

        # Multiplicity < 1 → error
        @test_throws Exception EnzymeRates.RegulatorMults(
            EnzymeRates.AllostericRegulator(:A), [0, 1])
        @test_throws Exception EnzymeRates.RegulatorMults(
            EnzymeRates.AllostericRegulator(:A), [-1])

        # CompetitiveInhibitor also accepted
        rm_ci = EnzymeRates.RegulatorMults(
            EnzymeRates.CompetitiveInhibitor(:I), [1])
        @test EnzymeRates.regulator(rm_ci) ==
              EnzymeRates.CompetitiveInhibitor(:I)
    end

    @testset "EnzymeReaction (new concrete)" begin
        r = EnzymeReaction(
            [EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:ATP), [:C => 10]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Product(:ADP), [:C => 10])],
            [EnzymeRates.RegulatorMults(
                 EnzymeRates.AllostericRegulator(:cAMP), [2])],
            [1, 2],
        )

        @test length(EnzymeRates.reactants(r)) == 2
        @test EnzymeRates.allowed_catalytic_multiplicities(r) == [1, 2]
        @test length(EnzymeRates.regulators(r)) == 1
        @test EnzymeRates.substrates(r) == [EnzymeRates.Substrate(:ATP)]
        @test EnzymeRates.products(r) == [EnzymeRates.Product(:ADP)]
    end

    @testset "EnzymeReaction canonicalizes reactant + regulator ordering" begin
        r1 = EnzymeReaction(
            [EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:B), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:A), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Product(:P), [:C => 2])],
            [EnzymeRates.RegulatorMults(
                 EnzymeRates.AllostericRegulator(:Y), [2]),
             EnzymeRates.RegulatorMults(
                 EnzymeRates.AllostericRegulator(:X), [2])],
            [3, 1, 2],
        )
        r2 = EnzymeReaction(
            [EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:A), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:B), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Product(:P), [:C => 2])],
            [EnzymeRates.RegulatorMults(
                 EnzymeRates.AllostericRegulator(:X), [2]),
             EnzymeRates.RegulatorMults(
                 EnzymeRates.AllostericRegulator(:Y), [2])],
            [1, 2, 3],
        )
        @test r1 == r2
        @test hash(r1) == hash(r2)
    end

    @testset "EnzymeReaction rejects multiplicity < 1" begin
        @test_throws ErrorException EnzymeReaction(
            EnzymeRates.ReactantAtoms[
                EnzymeRates.ReactantAtoms(
                    EnzymeRates.Substrate(:S), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [0, 1],
        )
    end

    @testset "Mechanism (non-parametric)" begin
        r = EnzymeReaction(
            [EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Product(:P), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [1],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

        s_bind = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        s_iso  = EnzymeRates.Step(e_s, e_p, nothing, false)
        s_rel  = EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)

        m = EnzymeRates.Mechanism(r, [[s_bind], [s_iso], [s_rel]])
        @test EnzymeRates.reaction(m) == r
        @test EnzymeRates.steps(m) == [[s_bind], [s_iso], [s_rel]]
        @test EnzymeRates.kinetic_groups(m) == 1:3
        @test EnzymeRates.n_steps(m) == 3
        @test EnzymeRates.rep_step(m, 2) == s_iso

        m2 = EnzymeRates.Mechanism(r, [[s_bind], [s_iso], [s_rel]])
        @test m == m2
        @test hash(m) == hash(m2)
    end

    @testset "iso canonicalization (RE + SS, all tiers)" begin
        # Tier 1 (SS iso, score differs): forward = substrate-bound -> product-bound.
        m_fwd = @enzyme_mechanism begin
            substrates: S; products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
            end
        end
        m_rev = @enzyme_mechanism begin
            substrates: S; products: P
            steps: begin
                E + S <--> E(S)
                E(P) <--> E(S)        # SS iso written backwards
                E(P) <--> E + P
            end
        end
        @test EnzymeRates.Mechanism(m_fwd) == EnzymeRates.Mechanism(m_rev)

        # Tier 2 (the Segel Iso Uni Uni case): pure conformational F <--> E,
        # tied on Tier 1, decided by entry_kind. Two opposite-source-direction
        # mechanisms canonicalize to the same form.
        s1 = @enzyme_mechanism begin
            substrates: A; products: P
            steps: begin
                E + A <--> E(A); E(A) <--> E(P); E(P) <--> F + P; F <--> E
            end
        end
        s2 = @enzyme_mechanism begin
            substrates: A; products: P
            steps: begin
                E + A <--> E(A); E(A) <--> E(P); E(P) <--> F + P; E <--> F   # last step flipped
            end
        end
        @test EnzymeRates.Mechanism(s1) == EnzymeRates.Mechanism(s2)
        iso = only(s for grp in EnzymeRates.steps(EnzymeRates.Mechanism(s2))
                       for s in grp
                       if !EnzymeRates.is_binding(s) &&
                          EnzymeRates.bound_metabolite(s) === nothing &&
                          EnzymeRates.name(EnzymeRates.from_species(s)) in (:E, :F))
        @test EnzymeRates.name(EnzymeRates.from_species(iso)) == :F  # product-exit
        @test EnzymeRates.name(EnzymeRates.to_species(iso))   == :E  # substrate-entry
    end

    @testset "Mechanism stores steps in flat order" begin
        r = EnzymeReaction(
            [EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Product(:P), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [1],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

        s1 = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        s2 = EnzymeRates.Step(e_s, e_p, nothing, false)
        s3 = EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)

        m = EnzymeRates.Mechanism(r, [[s1], [s2], [s3]])
        flat = EnzymeRates._flat_steps(m)
        @test [s for (s, _) in flat] == [s1, s2, s3]
        @test [g for (_, g) in flat] == [1, 2, 3]
    end

    @testset "AllostericMechanism (non-parametric)" begin
        r = EnzymeReaction(
            [EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Product(:P), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [2],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)
        s_bind = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        s_iso  = EnzymeRates.Step(e_s, e_p, nothing, false)
        s_rel  = EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)

        site = EnzymeRates.RegulatorySite(
            [EnzymeRates.AllostericRegulator(:I)], 1, [:OnlyI])
        m = EnzymeRates.AllostericMechanism(
            r, [[s_bind], [s_iso], [s_rel]],
            [:EqualAI, :NonequalAI, :OnlyA], 2,
            [site])

        @test EnzymeRates.reaction(m) == r
        @test EnzymeRates.steps(m) == [[s_bind], [s_iso], [s_rel]]
        @test EnzymeRates.cat_allo_state(m, 1) == :EqualAI
        @test EnzymeRates.cat_allo_state(m, 3) == :OnlyA
        @test EnzymeRates.catalytic_multiplicity(m) == 2
        @test EnzymeRates.regulatory_sites(m) == [site]
        @test EnzymeRates.kinetic_groups(m) == 1:3
        @test EnzymeRates.n_steps(m) == 3
        @test EnzymeRates.rep_step(m, 2) == s_iso
        @test EnzymeRates.allosteric_regulators(m) ==
              [EnzymeRates.AllostericRegulator(:I)]

        m2 = EnzymeRates.AllostericMechanism(
            r, [[s_bind], [s_iso], [s_rel]],
            [:EqualAI, :NonequalAI, :OnlyA], 2,
            [site])
        @test m == m2
        @test hash(m) == hash(m2)
    end

    @testset "AllostericMechanism validation errors" begin
        r = EnzymeReaction(
            [EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Product(:P), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [1],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        s_bind = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        cat_steps = [[s_bind]]

        # :OnlyI for catalytic group is rejected (R-state-active convention)
        @test_throws ErrorException EnzymeRates.AllostericMechanism(
            r, cat_steps, [:OnlyI], 1,
            EnzymeRates.RegulatorySite[])

        # Length mismatch
        @test_throws ErrorException EnzymeRates.AllostericMechanism(
            r, cat_steps, [:EqualAI, :NonequalAI], 1,
            EnzymeRates.RegulatorySite[])

        # catalytic_multiplicity < 1
        @test_throws ErrorException EnzymeRates.AllostericMechanism(
            r, cat_steps, [:EqualAI], 0,
            EnzymeRates.RegulatorySite[])

        # Unknown allo-state symbol
        @test_throws ErrorException EnzymeRates.AllostericMechanism(
            r, cat_steps, [:Bogus], 1,
            EnzymeRates.RegulatorySite[])
    end

    @testset "AllostericMechanism(::AllostericEnzymeMechanism) converter" begin
        cm = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        aem = EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:EqualAI, :NonequalAI, :OnlyA)),
            (((:A, :B), 1, (:OnlyA, :NonequalAI)),),
        )

        am = EnzymeRates.AllostericMechanism(aem)
        @test am isa EnzymeRates.AllostericMechanism

        # Catalytic side lifted via Mechanism(CM())
        @test EnzymeRates.steps(am) == EnzymeRates.Mechanism(cm).steps
        @test EnzymeRates.reaction(am) == EnzymeRates.Mechanism(cm).reaction

        # cat_allo_states and multiplicity extracted from CS
        @test EnzymeRates.cat_allo_state(am, 1) === :EqualAI
        @test EnzymeRates.cat_allo_state(am, 2) === :NonequalAI
        @test EnzymeRates.cat_allo_state(am, 3) === :OnlyA
        @test EnzymeRates.catalytic_multiplicity(am) == 2

        # Regulatory sites: ligand Symbols wrapped as AllostericRegulator
        sites = EnzymeRates.regulatory_sites(am)
        @test length(sites) == 1
        @test sites[1].ligands ==
              [EnzymeRates.AllostericRegulator(:A),
               EnzymeRates.AllostericRegulator(:B)]
        @test sites[1].multiplicity == 1
        @test sites[1].allo_states == [:OnlyA, :NonequalAI]

        # Idempotent: two calls give equal results
        @test EnzymeRates.AllostericMechanism(aem) == am

        # _to_mechanism bridge dispatches to AllostericMechanism
        @test EnzymeRates._to_mechanism(aem) == am
        @test EnzymeRates._to_mechanism(aem) isa EnzymeRates.AllostericMechanism
        # And matches the non-allosteric side too
        @test EnzymeRates._to_mechanism(cm) == EnzymeRates.Mechanism(cm)
    end

    @testset "name(p, ::AllostericEnzymeMechanism) chokepoint overloads" begin
        cm = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        aem = EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalAI, :EqualAI, :NonequalAI)),
            (((:R,), 1, (:NonequalAI,)),),
        )
        am = EnzymeRates.AllostericMechanism(aem)

        # Step-bound parameters: AEM dispatch matches AM dispatch
        rep_bind = first(EnzymeRates.steps(am)[1])
        @test EnzymeRates.name(EnzymeRates.Kd(rep_bind, :None), aem) ==
              EnzymeRates.name(EnzymeRates.Kd(rep_bind, :None), am)
        @test EnzymeRates.name(EnzymeRates.Kd(rep_bind, :I), aem) ==
              EnzymeRates.name(EnzymeRates.Kd(rep_bind, :I), am)
        @test EnzymeRates.name(EnzymeRates.Kd(rep_bind, :I), aem) === :K_I_S_E

        rep_iso  = first(EnzymeRates.steps(am)[2])
        @test EnzymeRates.name(EnzymeRates.Kfor(rep_iso, :None), aem) ==
              EnzymeRates.name(EnzymeRates.Kfor(rep_iso, :None), am)
        @test EnzymeRates.name(EnzymeRates.Kfor(rep_iso, :None), aem) === :k_ES_to_EP

        # Kreg: AEM dispatch matches AM dispatch
        site = EnzymeRates.regulatory_sites(am)[1]
        lig  = first(site.ligands)
        @test EnzymeRates.name(EnzymeRates.Kreg(site, lig, :A), aem) ==
              EnzymeRates.name(EnzymeRates.Kreg(site, lig, :A), am)
        @test EnzymeRates.name(EnzymeRates.Kreg(site, lig, :I), aem) ===
              :K_I_Rreg
    end

    @testset "_sig_of / _mechanism_from_sig roundtrip" begin
        # Multi-element atom data exercises type-parameter validity:
        # Pair{Symbol,Int} is NOT a valid type-parameter value; encoding
        # must use Tuple{Symbol,Int} leaves. Substrate and product carry the
        # same atom totals so the declared reaction balances (this test checks
        # Sig round-trip + type-param leaves, not reaction chemistry).
        atoms = [:C => 10, :H => 16, :N => 5, :O => 13, :P => 3]
        r = EnzymeRates.EnzymeReaction(
            [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:ATP), atoms),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:ADP), copy(atoms))],
            EnzymeRates.RegulatorMults[],
            [1],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:ATP)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:ADP)], :E)

        m = EnzymeRates.Mechanism(r, [
            [EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:ATP), true)],
            [EnzymeRates.Step(e_s, e_p, nothing, false)],
            [EnzymeRates.Step(e, e_p, EnzymeRates.Product(:ADP), true)],
        ])

        sig = EnzymeRates._sig_of(m)
        @test sig isa Tuple

        m_recon = EnzymeRates._mechanism_from_sig(sig)
        @test m_recon == m   # roundtrip

        # CRITICAL: sig MUST be usable as a type parameter. Throws TypeError
        # if any leaf is invalid (Pair, Vector, DataType inside value-tuple).
        em_type = EnzymeRates.EnzymeMechanism{sig}
        @test em_type <: EnzymeRates.EnzymeMechanism
        em_inst = em_type()
        @test em_inst isa EnzymeRates.EnzymeMechanism

        # Roundtrip through the type-parameter form preserves everything.
        @test EnzymeRates.Mechanism(em_inst) == m
    end

    @testset "Mechanism <-> EnzymeMechanism converters" begin
        r = EnzymeRates.EnzymeReaction(
            [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [1],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

        m = EnzymeRates.Mechanism(r, [
            [EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)],
            [EnzymeRates.Step(e_s, e_p, nothing, false)],
            [EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)],
        ])

        em = EnzymeMechanism(m)
        @test em isa EnzymeMechanism

        m_back = EnzymeRates.Mechanism(em)
        @test m_back == m
    end

    @testset "name(p::Parameter, m) chokepoint" begin
        r = EnzymeRates.EnzymeReaction(
            [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [1],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

        step1 = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        step2 = EnzymeRates.Step(e_s, e_p, nothing, false)
        step3 = EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)

        m = EnzymeRates.Mechanism(r, [[step1], [step2], [step3]])

        # Structural naming: binding params encode metabolite + pre-binding form;
        # iso params encode directed species pair.
        @test EnzymeRates.name(EnzymeRates.Kd(step1, :None), m) === :K_S_E
        @test EnzymeRates.name(EnzymeRates.Kd(step1, :I),    m) === :K_I_S_E
        @test EnzymeRates.name(EnzymeRates.Kon(step2, :None), m) === :k_ES_to_EP
        @test EnzymeRates.name(EnzymeRates.Koff(step2, :None), m) === :k_EP_to_ES
        @test EnzymeRates.name(EnzymeRates.Kfor(step2, :None), m) === :k_ES_to_EP
        @test EnzymeRates.name(EnzymeRates.Krev(step2, :None), m) === :k_EP_to_ES
        @test EnzymeRates.name(EnzymeRates.Kd(step3, :None), m) === :K_P_E

        # I-state token on SS step
        @test EnzymeRates.name(EnzymeRates.Kon(step2, :I),  m) === :k_I_ES_to_EP
        @test EnzymeRates.name(EnzymeRates.Koff(step2, :I), m) === :k_I_EP_to_ES

        # Kiso uses Kiso_ prefix (RE iso)
        @test EnzymeRates.name(EnzymeRates.Kiso(step2, :None), m) === :Kiso_ES_to_EP
        @test EnzymeRates.name(EnzymeRates.Kiso(step2, :I),    m) === :Kiso_I_ES_to_EP

        # Mechanism-level scalars
        @test EnzymeRates.name(EnzymeRates.Keq(),   m) === :Keq
        @test EnzymeRates.name(EnzymeRates.Etot(),  m) === :E_total
        @test EnzymeRates.name(EnzymeRates.Lallo(), m) === :L

        # Same names resolve via EnzymeMechanism(m) (the parametric form).
        em = EnzymeMechanism(m)
        @test EnzymeRates.name(EnzymeRates.Kd(step1, :None), em) === :K_S_E
        @test EnzymeRates.name(EnzymeRates.Kon(step2, :None), em) === :k_ES_to_EP
        @test EnzymeRates.name(EnzymeRates.Keq(),   em) === :Keq
        @test EnzymeRates.name(EnzymeRates.Etot(),  em) === :E_total
        @test EnzymeRates.name(EnzymeRates.Lallo(), em) === :L
    end

    @testset "name(p::Parameter, m) rep_idx for shared kinetic group" begin
        # Group with 2 steps: rep is the first step's position in the
        # flattened step list. If group 1 contains steps at positions 1
        # and 2, rep_idx is 1; if group 2 starts at position 3 with two
        # steps, rep_idx for group 2 is 3.
        r = EnzymeRates.EnzymeReaction(
            [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [1],
        )
        e    = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s  = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p  = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)
        e_sp = EnzymeRates.Species(
            EnzymeRates.Metabolite[
                EnzymeRates.Substrate(:S), EnzymeRates.Product(:P)], :E)

        step_a = EnzymeRates.Step(e,   e_s,  EnzymeRates.Substrate(:S), true)
        step_b = EnzymeRates.Step(e_p, e_sp, EnzymeRates.Substrate(:S), true)
        step_c = EnzymeRates.Step(e_s, e_p,  nothing, false)
        step_d = EnzymeRates.Step(e,   e_p,  EnzymeRates.Product(:P), true)

        m = EnzymeRates.Mechanism(r, [[step_a, step_b], [step_c], [step_d]])

        # Group 1: both steps bind S; rep = step_a. Both yield the same name.
        @test EnzymeRates.name(EnzymeRates.Kd(step_a, :None), m) === :K_S_E
        @test EnzymeRates.name(EnzymeRates.Kd(step_b, :None), m) === :K_S_E
        # Group 2: SS iso ES → EP.
        @test EnzymeRates.name(EnzymeRates.Kon(step_c, :None), m) === :k_ES_to_EP
        # Group 3: RE binding P from E.
        @test EnzymeRates.name(EnzymeRates.Kd(step_d, :None), m) === :K_P_E
    end

    @testset "name(p::Kreg, m) chokepoint" begin
        r = EnzymeRates.EnzymeReaction(
            [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
            [EnzymeRates.RegulatorMults(
                 EnzymeRates.AllostericRegulator(:A), [2])],
            [2],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

        cat_steps = [
            [EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)],
            [EnzymeRates.Step(e_s, e_p, nothing, false)],
            [EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)],
        ]
        site_a = EnzymeRates.RegulatorySite(
            [EnzymeRates.AllostericRegulator(:A)], 2, [:NonequalAI])
        am = EnzymeRates.AllostericMechanism(
            r, cat_steps,
            [:EqualAI, :EqualAI, :EqualAI], 2, [site_a])

        @test EnzymeRates.name(
            EnzymeRates.Kreg(site_a, EnzymeRates.AllostericRegulator(:A), :A),
            am) === :K_A_Areg
        @test EnzymeRates.name(
            EnzymeRates.Kreg(site_a, EnzymeRates.AllostericRegulator(:A), :I),
            am) === :K_I_Areg

        # Step-bound parameters also resolve via AllostericMechanism.
        rep = first(cat_steps[1])
        @test EnzymeRates.name(EnzymeRates.Kd(rep, :None), am) === :K_S_E
        @test EnzymeRates.name(EnzymeRates.Kd(rep, :I),    am) === :K_I_S_E

        # Iso step in second kinetic group
        iso_step = first(cat_steps[2])
        @test EnzymeRates.name(EnzymeRates.Kiso(iso_step, :None), am) === :Kiso_ES_to_EP
        @test EnzymeRates.name(EnzymeRates.Kon(iso_step, :None),  am) === :k_ES_to_EP

        # Scalars also dispatch on AllostericMechanism
        @test EnzymeRates.name(EnzymeRates.Keq(),   am) === :Keq
        @test EnzymeRates.name(EnzymeRates.Etot(),  am) === :E_total
        @test EnzymeRates.name(EnzymeRates.Lallo(), am) === :L
    end

    @testset "structural parameter names" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
            end
        end
        ps = EnzymeRates.parameters(m)
        # SS iso forward step E(S) → E(P) must appear as :k_ES_to_EP
        @test :k_ES_to_EP in ps
        # No positional index names such as k2f or K1.
        @test !any(p -> occursin(r"^k[0-9]", String(p)), ps)
        @test !any(p -> occursin(r"^K[0-9]", String(p)), ps)
    end

    @testset "synth-dep I-state names consistent with chokepoint (NonequalAI)" begin
        # PK-like mechanism: NonequalAI PEP binding, EqualAI catalysis.
        # k5r is a Haldane dep whose RHS references K_PEP_E (NonequalAI),
        # so a synthesized I-state dep is produced. The synth-dep name must
        # be what name(_flip_to_inactive(_param_for_symbol(am, active)), am)
        # returns, not string(active) * "_T".
        m = @allosteric_mechanism begin
            substrates: PEP, ADP
            products:   Pyruvate, ATP
            allosteric_regulators: ATP::OnlyI, F16BP::OnlyA
            catalytic_multiplicity: 4
            catalytic_steps: begin
                (E + PEP ⇌ E(PEP),
                 E(ADP) + PEP ⇌ E(PEP, ADP))                          :: NonequalAI
                (E + ADP ⇌ E(ADP),
                 E(PEP) + ADP ⇌ E(PEP, ADP))                          :: EqualAI
                E(PEP, ADP) <--> E(Pyruvate, ATP)                     :: EqualAI
                (E(Pyruvate, ATP) ⇌ E(ATP) + Pyruvate,
                 E(Pyruvate) ⇌ E + Pyruvate)                          :: EqualAI
                (E(Pyruvate, ATP) ⇌ E(Pyruvate) + ATP,
                 E(ATP) ⇌ E + ATP)                                    :: EqualAI
            end
            regulatory_site(multiplicity = 2): begin
                ligands: ATP
            end
            regulatory_site(multiplicity = 4): begin
                ligands: F16BP
            end
        end
        # parameters(m, Reduced) must advertise the I-state variant of PEP
        # binding using the structural mid-name I_ token, not a _T suffix.
        ps = EnzymeRates.parameters(m)
        ps_str = String.(collect(ps))
        # The I-state PEP binding param must use mid-name I_ token
        @test any(s -> startswith(s, "K_I_"), ps_str)
        # No param name should use _T suffix
        @test !any(s -> endswith(s, "_T"), ps_str)
        # Calling rate_equation must not error (proves indep names and rate
        # body names are consistent — the intermediate state would KeyError here).
        rng = Random.MersenneTwister(1234)
        met_names = [:PEP, :ADP, :Pyruvate, :ATP, :F16BP]
        concs_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in met_names)
        concs = NamedTuple{Tuple(met_names)}(concs_vals)
        indep = EnzymeRates.fitted_params(m)
        param_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in indep)
        params = NamedTuple{indep}(param_vals)
        params = merge(params, (Keq=1.0, E_total=1.0))
        v = EnzymeRates.rate_equation(m, concs, params)
        @test isfinite(v)
    end

    @testset "_force_inactive forces :I regardless of tag" begin
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Product(:P)], :E)
        s = EnzymeRates.Step(e_s, e_p, nothing, false)

        # An :EqualAI parameter has no :I variant under _flip_to_inactive
        # (returns itself); _force_inactive must return the explicit :I variant.
        p_eq = EnzymeRates.Krev(s, :EqualAI)
        @test EnzymeRates._flip_to_inactive(p_eq) === p_eq          # unchanged
        @test EnzymeRates._force_inactive(p_eq) == EnzymeRates.Krev(s, :I)
        p_a = EnzymeRates.Krev(s, :A)
        @test EnzymeRates._force_inactive(p_a) == EnzymeRates.Krev(s, :I)
    end
end

# Chokepoint guard: no `Symbol("[KkVL]...")` literal is constructed outside
# parameter-name rendering bodies (the `name(::Parameter, m)` chokepoint).

const _CHOKEPOINT_PREFIX = r"^[KkVL][_a-zA-Z0-9]"

# Extract the function-name symbol from a signature expression.
# Handles `name(...)`, `name(...) where T`, `name(...)::Ret`, etc.
function _sig_fn_name(sig)
    while sig isa Expr && sig.head === :where
        sig = sig.args[1]
    end
    if sig isa Expr && sig.head === :(::)
        sig = sig.args[1]
    end
    return sig isa Expr && sig.head === :call ? sig.args[1] : nothing
end

# Extract the first positional arg type-annotation as a String.
function _sig_first_arg_str(sig)
    while sig isa Expr && sig.head === :where
        sig = sig.args[1]
    end
    if sig isa Expr && sig.head === :(::)
        sig = sig.args[1]
    end
    sig isa Expr && sig.head === :call || return ""
    args = sig.args[2:end]
    pos_args = filter(a -> !(a isa Expr && a.head === :parameters), args)
    isempty(pos_args) && return ""
    return string(pos_args[1])
end

# A method definition is a chokepoint body iff it is a `name` method
# dispatching on a Parameter subtype value.
function _is_chokepoint_def(expr)
    expr isa Expr || return false
    sig = if expr.head === :function && length(expr.args) >= 1
        expr.args[1]
    elseif expr.head === :(=) && expr.args[1] isa Expr &&
           expr.args[1].head in (:call, :where)
        expr.args[1]
    else
        return false
    end
    fn_name = _sig_fn_name(sig)
    fn_name === :name || return false
    arg_str = _sig_first_arg_str(sig)
    return occursin(
        r"Parameter|::(Kd|Kiso|Kon|Koff|Kfor|Krev|Kreg|Keq|Etot|Lallo)\b",
        arg_str)
end

# Reconstruct the string content of a `Symbol("...")` call. Supports
# both literal Strings and `:string` interpolation expressions like
# `Symbol("K\$idx")` → `Expr(:string, "K", :idx)`.
function _symbol_call_pattern(expr)
    expr isa Expr && expr.head === :call &&
        length(expr.args) >= 2 && expr.args[1] === :Symbol || return nothing
    arg2 = expr.args[2]
    if arg2 isa String
        return arg2
    elseif arg2 isa Expr && arg2.head === :string
        # Concatenate; non-String parts become a placeholder so the
        # prefix regex can still match (e.g., "K\$idx" → "K_"). The
        # placeholder must be a character class matched by
        # `_CHOKEPOINT_PREFIX` (`[_a-zA-Z0-9]`); underscore qualifies.
        return join(p isa String ? p : "_" for p in arg2.args)
    end
    return nothing
end

function _walk_violations!(expr, in_chokepoint::Bool, out::Vector{String})
    expr isa Expr || return
    if _is_chokepoint_def(expr)
        for child in expr.args
            _walk_violations!(child, true, out)
        end
    else
        pat = _symbol_call_pattern(expr)
        if pat !== nothing && occursin(_CHOKEPOINT_PREFIX, pat) && !in_chokepoint
            push!(out, "Symbol(\"$pat\")")
        end
        for child in expr.args
            _walk_violations!(child, in_chokepoint, out)
        end
    end
end

@testset "chokepoint: no Symbol(\"[KkVL]...\") outside parameter-name renderers" begin
    src_dir = joinpath(dirname(@__DIR__), "src")
    for f in readdir(src_dir; join=true)
        endswith(f, ".jl") || continue
        src = read(f, String)
        expr = Meta.parseall(src; filename=f)
        violations = String[]
        _walk_violations!(expr, false, violations)
        if !isempty(violations)
            @info "chokepoint violations" file=basename(f) violations
        end
        @test isempty(violations)
    end
end
