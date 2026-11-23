# ABOUTME: Tests of src/Derivation/thermodynamic_constraints.jl: the joint solve, the
# ABOUTME: parameter names and units, parameters and dependent_parameter_values.

# ── Parameter names and units ─────────────────────────────────────────────────────────

@testset "_constant_columns: the joint solve's columns and their names" begin
    M = ER.Mechanisms
    D = ER.Derivation
    col = D._StepParameter
    # Uni-uni: E + P → EP and E + S → ES at rapid equilibrium, ES → EP steady
    # state, each its own group. The isomerization's columns come first, by its shape,
    # its reverse rate, which runs against the reaction, before its forward one; then the
    # bindings, the one whose forms come later by name first (ES after EP).
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    c = D._constant_columns(m)
    @test c.columns == [col(3, 2, false), col(3, 1, false), col(2, 1, false),
                        col(1, 1, false)]
    @test c.names == [:k_EP_to_ES, :k_ES_to_EP, :K_ES_to_E_S, :K_EP_to_E_P]
    @test c.in_law == trues(4)
    @test c.flat == [1, 2, 3]

    # Allosteric uni-uni: E + P → EP :EqualAI, E + S → ES :NonequalAI and
    # ES → EP :OnlyA, one group each. The :NonequalAI and :OnlyA groups get inactive
    # columns; the inactive copies of the :OnlyA step lie outside the inactive law
    # (block 2) and come first, then the inactive binding (block 1), then the active
    # columns (block 0). An :OnlyA step's inactive copy takes its own name.
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: OnlyA
            EP ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    a = D._constant_columns(am)
    @test a.columns == [col(3, 2, true), col(3, 1, true), col(2, 1, true),
                        col(3, 2, false), col(3, 1, false), col(2, 1, false),
                        col(1, 1, false)]
    @test a.names == [:k_I_EP_to_ES, :k_I_ES_to_EP, :K_I_ES_to_E_S, :k_A_EP_to_ES,
                      :k_A_ES_to_EP, :K_A_ES_to_E_S, :K_EP_to_E_P]
    @test a.in_law == BitVector([0, 0, 1, 1, 1, 1, 1])
    @test a.flat == [1, 2, 3]
end

@testset "_regulator_constant: a regulator ligand's constant per conformation" begin
    M = ER.Mechanisms
    D = ER.Derivation
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: NonequalAI
            EP ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::OnlyI
        regulatory_site(2): X::NonequalAI
        regulatory_site(2): Y::EqualAI
        regulatory_site(2): Z::OnlyA
    end
    # Each ligand with its own tag, in the active and then the inactive conformation: an
    # :OnlyI ligand is absent from the active one and an :OnlyA ligand from the inactive
    # one; an :EqualAI ligand has one bare name in both.
    @test [D._regulator_constant(am, lig, tag, inactive)
           for (lig, tag) in ((:R, :OnlyI), (:X, :NonequalAI), (:Y, :EqualAI), (:Z, :OnlyA))
           for inactive in (false, true)] ==
          [nothing, :K_I_Rreg, :K_A_Xreg, :K_I_Xreg, :K_Yreg, :K_Yreg, :K_A_Zreg, nothing]
end

@testset "_Parameter family: step, regulator and allosteric parameters" begin
    # A step parameter serves every kind of edge: a slot of an edge of a mechanism's
    # catalytic graph in one conformation.
    P = ER.Derivation._StepParameter
    p = P(1, 1, false)
    @test p isa ER.Derivation._Parameter
    @test (p.step, p.slot, p.inactive) == (1, 1, false)
    @test p == P(1, 1, false)
    @test p != P(2, 1, false) && p != P(1, 2, false) && p != P(1, 1, true)

    # A regulator constant is its ligand's in one conformation; it needs no site, since a
    # ligand sits at one site.
    R = ER.Derivation._RegulatorParameter
    kr = R(:A, false)
    @test kr isa ER.Derivation._Parameter
    @test (kr.ligand, kr.inactive) == (:A, false)
    @test kr == R(:A, false)
    @test kr != R(:A, true) && kr != R(:B, false)

    @test ER.Derivation._AllostericParameter() isa ER.Derivation._Parameter
    @test ER.Derivation._AllostericParameter() == ER.Derivation._AllostericParameter()
end

@testset "parameter_units gives each constant's units" begin
    D = ER.Derivation
    # Each mechanism's fitted constants, then their units: row i is parameters(m)[i], its
    # power of the rate and its power of each concentration, in metabolites(m) order.

    # Uni-uni, rapid-equilibrium bindings and a steady-state isomerization: a binding's
    # dissociation constant has the units of what binds, the rate constant of a step that
    # takes nothing up those of the rate.
    uni_uni = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    @test parameters(uni_uni) == (:K_EP_to_E_P, :K_ES_to_E_S, :k_ES_to_EP)
    @test D.parameter_units(uni_uni) ==
          (rate = [0, 0, 1], concentrations = [0 1; 1 0; 0 0])

    # Steady-state bindings: a rate constant is the rate over what it takes up in its own
    # direction, so the binding of S has [S]⁻¹ and each release none.
    ss_bindings = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            ES <--> E + P
        end
    end
    @test parameters(ss_bindings) == (:k_ES_to_E_P, :k_ES_to_E_S, :k_E_S_to_ES)
    @test D.parameter_units(ss_bindings) ==
          (rate = [1, 1, 1], concentrations = [0 0; 0 0; -1 0])
    # A release's reverse rate constant binds the product: [P]⁻¹.
    ss_uni_uni = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP <--> E + P
        end
    end
    @test parameters(ss_uni_uni) ==
          (:k_EP_to_E_P, :k_ES_to_EP, :k_ES_to_E_S, :k_E_P_to_EP, :k_E_S_to_ES)
    @test D.parameter_units(ss_uni_uni) ==
          (rate = [1, 1, 1, 1, 1], concentrations = [0 0; 0 0; 0 0; 0 -1; -1 0])

    # A rapid-equilibrium Theorell–Chance step's constant is what it gives off over what
    # it takes up, [P]/[B]; the binding of Q to F is steady state.
    theorell_chance = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EQ + P
            EQ <--> FQ
            FQ <--> F + Q
            F <--> E
        end
    end
    @test parameters(theorell_chance) == (:K_EA_B_to_EQ_P, :K_EA_to_E_A, :k_EQ_to_FQ,
                                          :k_E_to_F, :k_FQ_to_F_Q, :k_F_Q_to_FQ, :k_F_to_E)
    @test D.parameter_units(theorell_chance) ==
          (rate = [0, 0, 1, 1, 1, 1, 1],
           concentrations = [0 -1 1 0; 1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 -1;
                             0 0 0 0])

    # At steady state the Theorell–Chance step's reverse rate constant takes up P, what
    # the step gives off as stored: the rate over [P]. The isomerization EQ <--> FQ
    # carries the Haldane relation.
    steady_theorell_chance = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EQ + P
            EQ <--> FQ
            FQ <--> F + Q
            F <--> E
        end
    end
    @test parameters(steady_theorell_chance) ==
          (:K_EA_to_E_A, :k_EA_B_to_EQ_P, :k_EQ_P_to_EA_B, :k_EQ_to_FQ, :k_E_to_F,
           :k_FQ_to_F_Q, :k_F_Q_to_FQ, :k_F_to_E)
    @test D.parameter_units(steady_theorell_chance) ==
          (rate = [0, 1, 1, 1, 1, 1, 1, 1],
           concentrations = [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 0; 0 0 0 0; 0 0 0 0;
                             0 0 0 -1; 0 0 0 0])

    # Two metabolites bound in one rapid-equilibrium step: the constant is the complex
    # over what it takes up, [A]⁻¹[B]⁻¹.
    two_in_rapid = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A + B ⇌ EAB
            EAB <--> EPQ
            EPQ <--> EQ + P
            EQ <--> E + Q
        end
    end
    @test first(parameters(two_in_rapid)) === :K_E_A_B_to_EAB
    units = D.parameter_units(two_in_rapid)
    @test (units.rate[1], units.concentrations[1, :]) == (0, [-1, -1, 0, 0])

    # A dead-end rapid-equilibrium isomerization's constant has no units.
    dead_end = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
            ES ⇌ FS
        end
    end
    @test parameters(dead_end) == (:K_EP_to_E_P, :K_ES_to_E_S, :K_ES_to_FS, :k_ES_to_EP)
    @test D.parameter_units(dead_end) ==
          (rate = [0, 0, 0, 1], concentrations = [0 1; 1 0; 0 0; 0 0])

    # A reactant's modifier copy takes the reactant's column.
    modifier_copy = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
            EP + S ⇌ EPSmod
        end
    end
    @test parameters(modifier_copy) ==
          (:K_EPSmod_to_EP_S, :K_EP_to_E_P, :K_ES_to_E_S, :k_ES_to_EP)
    @test D.parameter_units(modifier_copy) ==
          (rate = [0, 0, 0, 1], concentrations = [1 0; 0 1; 1 0; 0 0])

    # An allosteric mechanism ends with its regulator constants, which have the units of
    # their ligand, the last metabolite, and L, which has none.
    allosteric = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: EqualAI
            EP ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::NonequalAI
    end
    @test metabolites(allosteric) == [:S, :P, :R]
    @test parameters(allosteric) ==
          (:K_A_ES_to_E_S, :K_EP_to_E_P, :k_ES_to_EP, :K_A_Rreg, :K_I_Rreg, :L)
    @test D.parameter_units(allosteric) ==
          (rate = [0, 0, 1, 0, 0, 0],
           concentrations = [1 0 0; 0 1 0; 0 0 0; 0 0 1; 0 0 1; 0 0 0])
end

@testset "parameter_units gives each constant the units its name reads" begin
    # Each constant's power of the rate and of each concentration, exactly as the rule
    # reads them from the constant's name (`_testhelper_units_from_name`). Each call returns
    # arrays of its own.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        units = ER.Derivation.parameter_units(m)
        n, k = length(parameters(m)), length(metabolites(m))
        @test units.rate isa Vector{Int} && size(units.rate) == (n,)
        @test units.concentrations isa Matrix{Int} && size(units.concentrations) == (n, k)
        @test [(units.rate[i], units.concentrations[i, :]) for i in 1:n] ==
              [_testhelper_units_from_name(x, metabolites(m)) for x in parameters(m)]
        fill!(ER.Derivation.parameter_units(m).rate, 7)
        @test ER.Derivation.parameter_units(m) == units
    end
end

@testset "the fitted names of every spec" begin
    # parameters(m) of every MECHANISM_TEST_SPECS entry, written out: the names a fit
    # estimates, in parameters order.
    pinned = Dict(
        "Uni-Uni" =>
            (:k_ES_to_E_P, :k_ES_to_E_S, :k_E_S_to_ES),
        "Segel Uni Uni" =>
            (:k_EA_to_EP, :k_EA_to_E_A, :k_EP_to_E_P, :k_E_A_to_EA, :k_E_P_to_EP),
        "Segel Iso Uni Uni" =>
            (:k_EA_to_EP, :k_EA_to_E_A, :k_EP_to_F_P, :k_E_A_to_EA, :k_E_to_F,
             :k_F_P_to_EP, :k_F_to_E),
        "Segel Ordered Uni Bi" =>
            (:k_EA_to_EQ_P, :k_EA_to_E_A, :k_EQ_to_E_Q, :k_E_A_to_EA, :k_E_Q_to_EQ),
        "Segel Ordered Bi Bi" =>
            (:k_EAB_to_EA_B, :k_EAB_to_EQ_P, :k_EA_B_to_EAB, :k_EA_to_E_A, :k_EQ_to_E_Q,
             :k_E_A_to_EA, :k_E_Q_to_EQ),
        "Segel Theorell-Chance Bi Bi" =>
            (:k_EA_B_to_EQ_P, :k_EA_to_E_A, :k_EQ_to_E_Q, :k_E_A_to_EA, :k_E_Q_to_EQ),
        "Segel Ping Pong Bi Bi" =>
            (:k_EA_to_E_A, :k_EA_to_F_P, :k_E_A_to_EA, :k_E_Q_to_FB, :k_FB_to_E_Q,
             :k_F_B_to_FB, :k_F_P_to_EA),
        "Segel Ordered Ter Bi" =>
            (:k_EABC_to_EAB_C, :k_EABC_to_EQ_P, :k_EAB_C_to_EABC, :k_EAB_to_EA_B,
             :k_EA_B_to_EAB, :k_EA_to_E_A, :k_EQ_to_E_Q, :k_E_A_to_EA, :k_E_Q_to_EQ),
        "Segel Ordered Ter Ter" =>
            (:k_EABC_to_EAB_C, :k_EABC_to_EQR_P, :k_EAB_C_to_EABC, :k_EAB_to_EA_B,
             :k_EA_B_to_EAB, :k_EA_to_E_A, :k_EQR_P_to_EABC, :k_EQR_to_ER_Q, :k_ER_to_E_R,
             :k_E_A_to_EA, :k_E_R_to_ER),
        "Segel Bi Uni Uni Uni PP Ter Bi" =>
            (:k_EAB_to_EA_B, :k_EAB_to_F_P, :k_EA_B_to_EAB, :k_EA_to_E_A, :k_E_A_to_EA,
             :k_E_Q_to_FC, :k_FC_to_E_Q, :k_F_C_to_FC, :k_F_P_to_EAB),
        "Random-order Bi-Bi" =>
            (:k_EAB_to_EA_B, :k_EAB_to_EPQ, :k_EA_B_to_EAB, :k_EA_to_E_A, :k_EB_A_to_EAB,
             :k_EB_to_E_B, :k_EPQ_to_EQ_P, :k_EQ_P_to_EPQ, :k_EQ_to_E_Q, :k_E_A_to_EA,
             :k_E_B_to_EB, :k_E_Q_to_EQ),
        "Numerator: random chem-SS (RE/SS)" =>
            (:K_EP_to_E_P, :K_ES1_to_E_S1, :K_ES2_to_E_S2, :k_ES1S2_to_EP,
             :k_ES1S2_to_ES1_S2, :k_ES1_S2_to_ES1S2, :k_ES2_S1_to_ES1S2),
        "Numerator: RE-chemistry" =>
            (:K_EP1_to_E_P1, :K_EP2_to_E_P2, :K_ES1_to_E_S1, :K_ES2_to_E_S2,
             :k_EP1P2_to_EP1_P2, :k_EP1P2_to_EP2_P1, :k_EP1_P2_to_EP1P2,
             :k_ES1S2_to_ES1_S2, :k_ES1_S2_to_ES1S2, :k_ES2_S1_to_ES1S2),
        "Numerator: redundant SS-bind" =>
            (:K_EP_to_E_P, :K_ES1S2_to_ES1_S2, :K_ES1_to_E_S1, :K_ES2_to_E_S2,
             :k_ES1S2_to_EP, :k_ES2_S1_to_ES1S2),
        "Segel Bi Uni Uni Bi PP Ter Ter" =>
            (:k_EAB_to_EA_B, :k_EAB_to_F_P, :k_EA_B_to_EAB, :k_EA_to_E_A, :k_ER_Q_to_FC,
             :k_ER_to_E_R, :k_E_A_to_EA, :k_E_R_to_ER, :k_FC_to_ER_Q, :k_F_C_to_FC,
             :k_F_P_to_EAB),
        "Segel Bi Bi Uni Uni PP Ter Ter" =>
            (:k_EAB_to_EA_B, :k_EAB_to_FQ_P, :k_EA_B_to_EAB, :k_EA_to_E_A, :k_E_A_to_EA,
             :k_E_R_to_FC, :k_FC_to_E_R, :k_FC_to_F_C, :k_FQ_P_to_EAB, :k_FQ_to_F_Q,
             :k_F_C_to_FC),
        "Segel Hexa Uni Ping Pong" =>
            (:k_EA_to_E_A, :k_EA_to_F_P, :k_E_A_to_EA, :k_E_R_to_GC, :k_FB_to_F_B,
             :k_FB_to_G_Q, :k_F_B_to_FB, :k_F_P_to_EA, :k_GC_to_E_R, :k_G_C_to_GC,
             :k_G_Q_to_FB),
        "RE Uni-Uni" =>
            (:K_EA_to_E_A, :k_EA_to_E_P),
        "RE Ordered Bi-Bi" =>
            (:K_EAB_to_EA_B, :K_EA_to_E_A, :k_EAB_to_EQ_P, :k_EQ_to_E_Q, :k_E_Q_to_EQ),
        "RE Random Bi-Bi" =>
            (:K_EA_to_E_A, :K_EB_to_E_B, :k_EAB_to_EA_B, :k_EAB_to_EPQ, :k_EA_B_to_EAB,
             :k_EB_A_to_EAB, :k_EPQ_to_EQ_P, :k_EQ_P_to_EPQ, :k_EQ_to_E_Q, :k_E_Q_to_EQ),
        "Competitive Inhibitor" =>
            (:K_EP_to_E_P, :K_ER_to_E_R, :K_ES_to_E_S, :k_ES_to_EP),
        "Non-competitive Inhibitor" =>
            (:K_EP_to_E_P, :K_ER_to_E_R, :K_ES_to_E_S, :k_ES_to_EP),
        "Uncompetitive Inhibitor" =>
            (:K_EP_to_E_P, :K_ERS_to_ES_R, :K_ES_to_E_S, :k_ES_to_EP),
        "Essential Activator" =>
            (:K_EPR_to_ER_P, :K_ERS_to_ER_S, :K_ER_to_E_R, :k_ERS_to_EPR),
        "Non-essential Activator" =>
            (:K_EP_to_E_P, :K_ER_to_E_R, :K_ES_to_E_S, :k_ERS_to_EPR, :k_ES_to_EP),
        "Activator + Competitive Inhibitor" =>
            (:K_EA_to_E_A, :K_EI_to_E_I, :K_EP_to_E_P, :K_ES_to_E_S, :k_EAS_to_EAP,
             :k_ES_to_EP),
        "MWC Dimer" =>
            (:K_A_EP_to_E_P, :K_A_ES_to_E_S, :K_I_EP_to_E_P, :K_I_ES_to_E_S,
             :k_A_ES_to_EP, :k_I_ES_to_EP, :L),
        "MWC Homodimer + Non-competitive Inhibitor" =>
            (:K_A_EP_to_E_P, :K_A_ES_to_E_S, :K_I_EP_to_E_P, :K_I_ES_to_E_S,
             :k_A_ES_to_EP, :k_I_ES_to_EP, :K_A_Ireg, :K_I_Ireg, :L),
        "Two Competitive Inhibitors" =>
            (:K_EI1_to_E_I1, :K_EI2_to_E_I2, :K_EP_to_E_P, :K_ES_to_E_S, :k_ES_to_EP),
        "Two Non-competitive Inhibitors" =>
            (:K_EI1_to_E_I1, :K_EI2_to_E_I2, :K_EP_to_E_P, :K_ES_to_E_S, :k_ES_to_EP),
        "Non-competitive + Competitive Inhibitor" =>
            (:K_EI1_to_E_I1, :K_EI2_to_E_I2, :K_EP_to_E_P, :K_ES_to_E_S, :k_ES_to_EP),
        "Uncompetitive + Competitive Inhibitor" =>
            (:K_EI1S_to_ES_I1, :K_EI2_to_E_I2, :K_EP_to_E_P, :K_ES_to_E_S, :k_ES_to_EP),
        "Two Same-site Inhibitors" =>
            (:K_EI1_to_E_I1, :K_EI2_to_E_I2, :K_EP_to_E_P, :K_ES_to_E_S, :k_ES_to_EP),
        "MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites" =>
            (:K_A_EP1_to_E_P1, :K_A_EP2_to_E_P2, :K_A_ES1_to_E_S1, :K_A_ES2_to_E_S2,
             :K_I_EP1_to_E_P1, :K_I_EP2_to_E_P2, :K_I_ES1_to_E_S1, :K_I_ES2_to_E_S2,
             :k_A_ES1S2_to_EP1P2, :k_I_ES1S2_to_EP1P2, :K_A_R1reg, :K_A_R2reg, :K_A_R3reg,
             :K_I_R1reg, :K_I_R2reg, :K_I_R3reg, :L),
        "PFK-1" =>
            (:K_A_EF6P_to_E_F6P, :K_EADP_to_E_ADP, :K_EATP_to_E_ATP, :K_EF16BP_to_E_F16BP,
             :k_A_EATPF6P_to_EADPF16BP, :K_A_ADPreg, :K_Pireg, :K_A_F26BPreg, :K_I_ATPreg,
             :K_I_Citratereg, :K_I_F26BPreg, :L),
        "HK" =>
            (:K_A_EATP_to_E_ATP, :K_EADP_to_E_ADP, :K_EG6P_to_E_G6P, :K_EG6Pmod_to_E_G6P,
             :K_EGlucose_to_E_Glucose, :k_A_EATPGlucose_to_EADPG6P, :K_Pireg, :K_I_G6Preg,
             :L),
        "PK" =>
            (:K_A_EPEP_to_E_PEP, :K_EADP_to_E_ADP, :K_EATP_to_E_ATP,
             :K_EPyruvate_to_E_Pyruvate, :k_A_EADPPEP_to_EATPPyruvate, :K_A_F16BPreg,
             :K_I_ATPreg, :L),
        "m_all" =>
            (:K_A_EP1_to_E_P1, :K_A_ES1_to_E_S1, :K_EP2_to_E_P2, :K_ES2_to_E_S2,
             :K_I_EP1_to_E_P1, :K_I_ES1_to_E_S1, :k_A_ES1S2_to_EP1P2, :k_I_ES1S2_to_EP1P2,
             :K_A_R1reg, :K_R2reg, :K_I_R1reg, :L),
        "m_OnlyA_prod" =>
            (:K_A_EP_to_E_P, :K_ES_to_E_S, :k_A_ES_to_EP, :L),
        "LDH i-state NonequalAI 6-group, steady-state NAD group" =>
            (:K_EPyruvate_to_E_Pyruvate, :k_A_ELactate_to_E_Lactate,
             :k_A_ENADHPyruvate_to_ELactateNAD, :k_A_ENADH_to_E_NADH,
             :k_A_E_Lactate_to_ELactate, :k_A_E_NADH_to_ENADH, :k_ENAD_to_E_NAD,
             :k_E_NAD_to_ENAD, :k_I_ELactate_to_E_Lactate, :L),
        "LDH i-state NonequalAI 5-group" =>
            (:K_A_ELactate_to_E_Lactate, :K_ENADHPyruvate_to_ENADH_Pyruvate,
             :K_ENAD_to_E_NAD, :K_I_ELactate_to_E_Lactate,
             :k_A_ENADHPyruvate_to_ELactateNAD, :k_A_ENADH_to_E_NADH,
             :k_A_E_NADH_to_ENADH, :L),
        "LDH i-state EqualAI-NonequalAI 6-group" =>
            (:K_ELactateNAD_to_ENAD_Lactate, :K_ENADH_to_E_NADH, :K_ENAD_to_E_NAD,
             :K_EPyruvate_to_E_Pyruvate, :k_A_ENADPyruvate_to_EPyruvate_NAD,
             :k_ENADHPyruvate_to_ELactateNAD, :k_I_ENADPyruvate_to_EPyruvate_NAD, :L)
    )
    @test Set(keys(pinned)) == Set(spec.name for spec in MECHANISM_TEST_SPECS)
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test parameters(spec.mechanism) == pinned[spec.name]
    end
end

@testset "_name(p, m) chokepoint on an allosteric mechanism" begin
    am = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: NonequalAI
            EP ⇌ E + P      :: NonequalAI
        end
        regulatory_site(1): R::NonequalAI
    end

    # Edge parameters, edges of the catalytic graph picked by content.
    P = ER.Derivation._StepParameter
    rep_bind = _testhelper_edge(am, "E + S ⇌ ES")
    @test ER.Derivation._name(P(rep_bind, 1, false), am) === :K_ES_to_E_S
    @test ER.Derivation._name(P(rep_bind, 1, true), am) === :K_I_ES_to_E_S

    # An :EqualAI isomerization, steady-state and at rapid equilibrium, has no state tag.
    ss_iso = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: EqualAI
            EP ⇌ E + P      :: NonequalAI
        end
        regulatory_site(1): R::NonequalAI
    end
    rep_iso = _testhelper_edge(ss_iso, "ES <--> EP")
    @test ER.Derivation._name(P(rep_iso, 1, false), ss_iso) === :k_ES_to_EP
    @test ER.Derivation._name(P(rep_iso, 2, false), ss_iso) === :k_EP_to_ES
    re_iso = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + S <--> ES   :: NonequalAI
            ES ⇌ EP       :: EqualAI
            EP <--> E + P   :: NonequalAI
        end
        regulatory_site(1): R::NonequalAI
    end
    @test ER.Derivation._name(P(_testhelper_edge(re_iso, "ES ⇌ EP"), 1, false), re_iso) ===
          :K_ES_to_EP

    lig = first(am.allosteric.regulatory_sites[1].ligands)
    R = ER.Derivation._RegulatorParameter
    @test ER.Derivation._name(R(lig, false), am) === :K_A_Rreg
    @test ER.Derivation._name(R(lig, true), am) === :K_I_Rreg
    @test ER.Derivation._name(ER.Derivation._AllostericParameter(), am) === :L

end

@testset "_name(p::_Parameter, m) chokepoint" begin
    m = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    # The binding of S, the isomerization and the binding of P, as edges of `m`.
    bind, iso, rel = (_testhelper_edge(m, text)
                      for text in ("E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"))

    P = ER.Derivation._StepParameter

    # Structural naming: every step constant encodes its reaction's two sides;
    # a binding K reads in the release direction, iso params in the stored one.
    @test ER.Derivation._name(P(bind, 1, false), m) === :K_ES_to_E_S
    @test ER.Derivation._name(P(iso, 1, false), m) === :k_ES_to_EP
    @test ER.Derivation._name(P(iso, 2, false), m) === :k_EP_to_ES
    @test ER.Derivation._name(P(rel, 1, false), m) === :K_EP_to_E_P

    # A binding's rate pair: the forward rate binds, the reverse rate releases
    ss_bind = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    k = _testhelper_edge(ss_bind, "E + S <--> ES")
    @test ER.Derivation._name(P(k, 1, false), ss_bind) === :k_E_S_to_ES
    @test ER.Derivation._name(P(k, 2, false), ss_bind) === :k_ES_to_E_S

    # I-state token on RE and SS steps
    am_ss = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: NonequalAI
            EP ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::NonequalAI
    end
    k = _testhelper_edge(am_ss, "E + S ⇌ ES")
    @test ER.Derivation._name(P(k, 1, true), am_ss) === :K_I_ES_to_E_S
    k = _testhelper_edge(am_ss, "ES <--> EP")
    @test ER.Derivation._name(P(k, 1, true), am_ss) === :k_I_ES_to_EP
    @test ER.Derivation._name(P(k, 2, true), am_ss) === :k_I_EP_to_ES

    # The RE iso's equilibrium constant, named in the stored direction
    re_iso = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> ES
            ES ⇌ EP
            EP <--> E + P
        end
    end
    @test ER.Derivation._name(P(_testhelper_edge(re_iso, "ES ⇌ EP"), 1, false), re_iso) ===
          :K_ES_to_EP
    am_re = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + S <--> ES   :: EqualAI
            ES ⇌ EP       :: NonequalAI
            EP <--> E + P   :: EqualAI
        end
        regulatory_site(2): R::NonequalAI
    end
    @test ER.Derivation._name(P(_testhelper_edge(am_re, "ES ⇌ EP"), 1, true), am_re) ===
          :K_I_ES_to_EP
end

@testset "_name(p::_Parameter, m) for the steps of a shared kinetic group" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E + S ⇌ ES, EP + S ⇌ EPS)
            ES <--> EP
            EP ⇌ E + P
        end
    end

    # Both steps bind S; rep = bind. Both yield the same name.
    @test ER.Derivation._name(
        ER.Derivation._StepParameter(_testhelper_edge(m, "E + S ⇌ ES"), 1, false), m) ===
          :K_ES_to_E_S
    @test ER.Derivation._name(
        ER.Derivation._StepParameter(_testhelper_edge(m, "EP + S ⇌ EPS"), 1, false), m) ===
          :K_ES_to_E_S
end

@testset "fused steps are named by their sides" begin
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EAB
            EAB <--> EQ + P
            EQ <--> E + Q
        end
    end
    names = _testhelper_step_constant_names(m)
    @test :k_EAB_to_EQ_P in names && :k_EQ_P_to_EAB in names
    # The fused release EAB → EQ + P is stored as the binding of P it
    # reverses, EQ + P → EAB, so its forward rate constant binds P.
    k = _testhelper_edge(m, "EQ + P <--> EAB")
    kon, koff = _testhelper_step_parameters(m, k)
    @test (kon.slot, koff.slot) == (1, 2)
    @test ER.Derivation._name(kon, m) == :k_EQ_P_to_EAB &&
          ER.Derivation._name(koff, m) == :k_EAB_to_EQ_P
    @test :k_E_A_to_EA in names && :k_EQ_to_E_Q in names
end

@testset "step constants are named by their reaction" begin
    # Michaelis–Menten, every step steady-state: each binding has one rate constant
    # per direction, named after the reaction it drives; the isomerization has k both ways.
    mm_ss = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP <--> E + P
        end
    end
    @test issubset([:k_E_S_to_ES, :k_ES_to_E_S, :k_E_P_to_EP, :k_EP_to_E_P],
                   _testhelper_step_constant_names(mm_ss))
    # The Haldane reduction leaves the SS isomerization constant fitted.
    @test :k_ES_to_EP in ER.parameters(mm_ss)
    # An RE isomerization has one equilibrium constant, in its canonical direction.
    re_iso = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES ⇌ EP
            EP <--> E + P
        end
    end
    @test :K_ES_to_EP in _testhelper_step_constant_names(re_iso)
    # An RE fused release is stored as the binding of P it reverses, so its
    # constant is a dissociation constant, [EQ]·[P]/[EAB], named in the release
    # direction.
    re_fused = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB ⇌ EQ + P
            EQ <--> E + Q
        end
    end
    @test :K_EAB_to_EQ_P in _testhelper_step_constant_names(re_fused)
    # An RE fused binding takes a dissociation constant too, [E]·[S]/[EP], named in
    # the release direction like a plain binding's.
    re_fused_binding = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ EP
            EP <--> E + P
        end
    end
    @test :K_EP_to_E_S in _testhelper_step_constant_names(re_fused_binding)
    # A competitive-inhibitor copy of A and A itself bind E without colliding.
    inh = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A ⇌ EA
            EA <--> EP
            EP ⇌ E + P
            E + A ⇌ EAmod
        end
    end
    @test issubset([:K_EAmod_to_E_A, :K_EA_to_E_A],
                   _testhelper_step_constant_names(inh))
    _testhelper_check_law_against_oracle(inh)
    # The allosteric state tag follows the prefix: a :NonequalAI binding has an
    # A-state and an I-state constant, an :EqualAI binding one shared constant.
    allo = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES     :: NonequalAI
            ES <--> EP   :: EqualAI
            EP ⇌ E + P     :: EqualAI
        end
    end
    # Every constant the rate equation names, fitted or dependent.
    dep, indep = _testhelper_dependent_param_exprs(allo)
    allo_names = union(keys(dep), indep)
    @test issubset([:K_A_ES_to_E_S, :K_I_ES_to_E_S], allo_names)
    @test !(:K_ES_to_E_S in allo_names)
    @test :K_EP_to_E_P in ER.parameters(allo) && !(:K_A_EP_to_E_P in allo_names)
    @test :K_I_EP_to_E_P ∉ ER.parameters(allo)
    # A ping-pong residual form: B binds E_res_plusA_minusP, whose name is
    # E_res_plusA_minusP, so the release-direction K reads
    # EB_res_plusA_minusP → E_res_plusA_minusP + B.
    # E and the covalent E(; residual) share one RE segment whose weights vanish only
    # at B = Q = 0 (mixed), so the mechanism is accepted.
    pingpong = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            E + Q ⇌ EQ
        end
    end
    pp_names = _testhelper_step_constant_names(pingpong)
    @test :K_EB_res_plusA_minusP_to_E_res_plusA_minusP_B in pp_names
    _testhelper_check_law_against_oracle(pingpong)
    # One name per class: a rapid-equilibrium group's K, a steady-state group's kf, kr.
    @test length(pp_names) ==
          sum(pingpong.catalytic.steps[first(ks)].rapid_equilibrium ? 1 : 2
              for ks in _testhelper_group_members(pingpong.catalytic))
end

@testset "a modifier-only name is bare, a reactant's copy is Xmod in its form's name" begin
    # I binds only as a modifier, so nothing else carries its name.
    modifier_only = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
            E + I ⇌ EI
        end
    end
    @test _testhelper_step_constant_names(modifier_only) ==
          Set([:K_ES_to_E_S, :k_ES_to_EP, :k_EP_to_ES, :K_EP_to_E_P, :K_EI_to_E_I])
    _testhelper_check_law_against_oracle(modifier_only)
    # P is a product that also binds as a modifier: the form its modifier copy binds
    # into reads Pmod, and the constant names the metabolite it takes up, P.
    reactant_copy = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EPQ ⇌ EQ + P
            EQ ⇌ E + Q
            EQ + P ⇌ EPmodQ
        end
    end
    @test _testhelper_step_constant_names(reactant_copy) ==
          Set([:K_EA_to_E_A, :K_EAB_to_EA_B, :k_EAB_to_EPQ, :k_EPQ_to_EAB,
               :K_EPQ_to_EQ_P, :K_EQ_to_E_Q, :K_EPmodQ_to_EQ_P])
    _testhelper_check_law_against_oracle(reactant_copy)
end

@testset "a residual reads _res, then _plus per added and _minus per subtracted name" begin
    # The residual adds A and B and subtracts P; each run lists its names in order.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EAB
            EAB <--> EP_res_plusA_plusB_minusP
            EP_res_plusA_plusB_minusP <--> E_res_plusA_plusB_minusP + P
            E_res_plusA_plusB_minusP <--> EQ
            EQ <--> E + Q
        end
    end
    @test _testhelper_step_constant_names(m) == Set([
        :k_E_A_to_EA, :k_EA_to_E_A, :k_EA_B_to_EAB, :k_EAB_to_EA_B,
        :k_EAB_to_EP_res_plusA_plusB_minusP, :k_EP_res_plusA_plusB_minusP_to_EAB,
        :k_E_res_plusA_plusB_minusP_P_to_EP_res_plusA_plusB_minusP,
        :k_EP_res_plusA_plusB_minusP_to_E_res_plusA_plusB_minusP_P,
        :k_E_res_plusA_plusB_minusP_to_EQ, :k_EQ_to_E_res_plusA_plusB_minusP,
        :k_E_Q_to_EQ, :k_EQ_to_E_Q])
    # Every name is a Julia identifier, so a parameter tuple needs no var"…".
    @test all(Base.isidentifier, ER.parameters(m))
    _testhelper_check_law_against_oracle(m)
end

@testset "an isomerization runs the way the reaction turns over through it" begin
    # P leaves at F (the fused release ES → F + P) and S enters at E, so the cycle
    # E → ES → F → E turns the reaction over and the isomerization between F and E runs
    # F → E however either step is written.
    forward = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            ES <--> F + P
            F <--> E
        end
    end
    backward = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            F + P <--> ES
            E <--> F
        end
    end
    @test backward == forward
    k = only(k for (k, e) in enumerate(forward.catalytic.steps) if ER.Mechanisms._is_iso(e))
    c = forward.catalytic
    @test c.enzyme_forms[c.steps[k].from] === :F && c.enzyme_forms[c.steps[k].to] === :E
    kf, kr = _testhelper_step_parameters(forward, k)
    @test ER.Derivation._name(kf, forward) == :k_F_to_E &&
          ER.Derivation._name(kr, forward) == :k_E_to_F
end

# Chokepoint guard: no `Symbol("[KkVL]...")` literal is constructed outside
# parameter-name rendering bodies (the `_name(::_Parameter, m)` chokepoint).

const _CHOKEPOINT_PREFIX = r"^[KkVL][_a-zA-Z0-9]"

# Extract the function-name symbol from a signature expression.
# Handles `name(...)`, `name(...) where T`, `name(...)::Ret`, etc.
function _testhelper_sig_fn_name(sig)
    while sig isa Expr && sig.head === :where
        sig = sig.args[1]
    end
    if sig isa Expr && sig.head === :(::)
        sig = sig.args[1]
    end
    return sig isa Expr && sig.head === :call ? sig.args[1] : nothing
end

# Extract the first positional arg type-annotation as a String.
function _testhelper_sig_first_arg_str(sig)
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

# A method definition is a chokepoint body iff it is a `_name` method
# dispatching on a _Parameter subtype value.
function _testhelper_is_chokepoint_def(expr)
    expr isa Expr || return false
    sig = if expr.head === :function && length(expr.args) >= 1
        expr.args[1]
    elseif expr.head === :(=) && expr.args[1] isa Expr &&
           expr.args[1].head in (:call, :where)
        expr.args[1]
    else
        return false
    end
    fn_name = _testhelper_sig_fn_name(sig)
    fn_name === :_name || return false
    arg_str = _testhelper_sig_first_arg_str(sig)
    return occursin("Parameter", arg_str)
end

# Reconstruct the string content of a `Symbol("...")` call. Supports
# both literal Strings and `:string` interpolation expressions like
# `Symbol("K\$idx")` → `Expr(:string, "K", :idx)`.
function _testhelper_symbol_call_pattern(expr)
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

function _testhelper_walk_violations!(expr, in_chokepoint::Bool, out::Vector{String})
    expr isa Expr || return
    if _testhelper_is_chokepoint_def(expr)
        for child in expr.args
            _testhelper_walk_violations!(child, true, out)
        end
    else
        pat = _testhelper_symbol_call_pattern(expr)
        if pat !== nothing && occursin(_CHOKEPOINT_PREFIX, pat) && !in_chokepoint
            push!(out, "Symbol(\"$pat\")")
        end
        for child in expr.args
            _testhelper_walk_violations!(child, in_chokepoint, out)
        end
    end
end

"""The number of `_name(::_Parameter, m)` method definitions in `expr`."""
_testhelper_count_chokepoint_defs(expr) =
    expr isa Expr ?
    Int(_testhelper_is_chokepoint_def(expr)) +
    sum(_testhelper_count_chokepoint_defs, expr.args; init = 0) : 0

@testset "chokepoint: no Symbol(\"[KkVL]...\") outside parameter-name renderers" begin
    src_dir = normpath(joinpath(@__DIR__, "..", "..", "src"))
    names_file = joinpath(src_dir, "Derivation", "thermodynamic_constraints.jl")
    for f in sort!([joinpath(d, f) for (d, _, fs) in walkdir(src_dir) for f in fs])
        endswith(f, ".jl") || continue
        src = read(f, String)
        expr = Meta.parseall(src; filename=f)
        violations = String[]
        _testhelper_walk_violations!(expr, false, violations)
        if !isempty(violations)
            @info "chokepoint violations" file=basename(f) violations
        end
        @test isempty(violations)
        # The renderers are in src/Derivation/thermodynamic_constraints.jl and nowhere else.
        @test (_testhelper_count_chokepoint_defs(expr) > 0) == (f == names_file)
    end
end

# ── The joint solve ───────────────────────────────────────────────────────────────────

const _testhelper_JS = EnzymeRates

"""
A parameter point of `m` built from its joint solve: every column that is not a pivot
takes a random log value, `log Keq` one more, and each pivot column the value its row of
the reduced system gives. Returns the steps of each conformation as
`_testhelper_steady_state_rate` reads them (the steps of `_testhelper_graph`, whose order
is the solve's flat order, `_flat_edges`), with `L` and the regulator constants drawn at
random.
"""
function _testhelper_solve_point(m, rng)
    g = _testhelper_graph(m)
    sys = _testhelper_JS.Derivation._joint_system(m)
    law = _testhelper_JS.Mechanisms.steps_running_in_allosteric_conformations(m).inactive
    tags = _testhelper_JS.Derivation._solve_tags(m)
    edges = _testhelper_JS.Derivation._flat_edges(m.catalytic)
    M, pivots = _testhelper_JS.Derivation._reduce_joint_system(sys.A, sys.rhs)
    n = length(sys.columns)
    x = zeros(BigFloat, n)
    free = setdiff(1:n, pivots)
    foreach(f -> x[f] = 4 * big(rand(rng)) - 2, free)
    log_keq = 4 * big(rand(rng)) - 2
    for (r, p) in enumerate(pivots)
        x[p] = M[r, n + 1] * log_keq -
               sum((M[r, f] * x[f] for f in free); init = big(0))
    end
    value = Dict(c => exp(x[i]) for (i, c) in enumerate(sys.columns))
    allo = g.oligomeric_state !== nothing
    function _testhelper_conformation(inactive)
        out = []
        for (j, s) in enumerate(g.steps)
            k = edges[j]
            inactive && !law[k] && continue
            own = inactive && _testhelper_JS.Derivation._has_inactive_columns(tags[s.group])
            col(slot) = value[_testhelper_JS.Derivation._StepParameter(k, slot, own)]
            if s.rapid_equilibrium
                K = col(1)
                push!(out, (s, _testhelper_is_binding(s) ? 1 / K : K, nothing))
            else
                push!(out, (s, col(1) / col(2), col(1)))
            end
        end
        out
    end
    vals = Dict{Symbol, BigFloat}()
    L = big(0)
    if allo
        L = exp(4 * big(rand(rng)) - 2)
        for K in [_testhelper_regulator_constants(m, :A);
                  _testhelper_regulator_constants(m, :I)]
            vals[K] = exp(4 * big(rand(rng)) - 2)
        end
    end
    (; vals, Keq = exp(log_keq), active = _testhelper_conformation(false),
     inactive = allo ? _testhelper_conformation(true) : [], L)
end

# Reversible Michaelis–Menten with every step at steady state.
const _testhelper_mm_ss = @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S <--> ES
        ES <--> EP
        EP <--> E + P
    end
end

@testset "Joint thermodynamic solve" begin
    @testset "the priority key orders the columns" begin
        sys = EnzymeRates.Derivation._joint_system(_testhelper_mm_ss)
        describe(c) = (_testhelper_edge_text(_testhelper_mm_ss, c.step), c.slot)
        # The isomerization goes first, its reverse rate (against the reaction) before
        # its forward one; then the two bindings' constants that run against the
        # reaction, the release of S and the uptake of P (a product binding's uptake runs
        # against it), the step whose forms come later by name first (the binding of S,
        # as ES comes after EP); then the others.
        @test describe.(sys.columns) == [("ES <--> EP", 2), ("ES <--> EP", 1),
                                         ("E + S <--> ES", 2), ("E + P <--> EP", 1),
                                         ("E + S <--> ES", 1), ("E + P <--> EP", 2)]
        # One Haldane row: the first column, the isomerization's reverse rate, is the
        # dependent one.
        M, pivots = EnzymeRates.Derivation._reduce_joint_system(sys.A, sys.rhs)
        @test pivots == [1]
        @test M[1, end] == -1
    end

    @testset "shape, free E, flag, direction, then forms and metabolites by name" begin
        # Every step at steady state: a binding at free E, a binding off free E, an
        # isomerization, a Theorell–Chance exchange, the releases of P and Q, each stored as
        # its uptake, and an isomerization of free E.
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> EA
                EA + B <--> EAB
                EAB <--> EPQ
                EA + B <--> EQ + P
                EPQ <--> EQ + P
                EQ <--> E + Q
                E <--> F
            end
        end
        # The binding off free E at rapid equilibrium.
        m_re = @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A <--> EA
                EA + B ⇌ EAB
                EAB <--> EP
                EP <--> E + P
            end
        end
        _testhelper_key(x, text, slot, block) =
            EnzymeRates.Derivation._priority_key(
                x.catalytic, _testhelper_edge(x, text), slot, block)
        key(text, slot, block) = _testhelper_key(m, text, slot, block)
        at_free_E, off_free_E = "E + A <--> EA", "EA + B <--> EAB"
        # block first: an inactive copy outside the law before any other column
        @test key(at_free_E, 1, 2) > key(off_free_E, 2, 1) > key(off_free_E, 2, 0)
        # shape next: an isomerization, then a step with more than one metabolite, then a
        # binding. The three tie on the next three fields, so shape alone orders them.
        iso, exchange = "EAB <--> EPQ", "EA + B <--> EQ + P"
        @test key(iso, 2, 0) > key(exchange, 2, 0) > key(off_free_E, 2, 0)
        @test key(iso, 2, 0)[3:5] == key(exchange, 2, 0)[3:5] == key(off_free_E, 2, 0)[3:5]
        # shape before the constant against the reaction and before free E
        @test key(iso, 1, 0) > key(exchange, 2, 0)
        @test key("E <--> F", 1, 0) > key(off_free_E, 1, 0)
        # a step neither of whose forms is free E before one that touches it
        @test key(off_free_E, 1, 0) > key(at_free_E, 2, 0)
        # steady state before rapid equilibrium, then the constant against the reaction: a
        # binding's release, or the uptake of a product
        @test key(off_free_E, 1, 0) > _testhelper_key(m_re, "EA + B ⇌ EAB", 1, 0)
        @test key(off_free_E, 2, 0) > key(off_free_E, 1, 0)
        @test key("E + Q <--> EQ", 1, 0) > key("E + Q <--> EQ", 2, 0)
        # A release of two products is stored as its uptake, E + P + Q → EAB, which
        # runs against the reaction, so its forward rate is the constant against it.
        two = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A + B <--> EAB
                EAB <--> E + P + Q
            end
        end
        @test _testhelper_key(two, "E + P + Q <--> EAB", 1, 0)[5] &&
              !_testhelper_key(two, "E + P + Q <--> EAB", 2, 0)[5]
        # then the indices of the step's two forms, which are stored E first and then by
        # name (EA, EAB, EPQ, EQ, F), the later forms first
        @test key("E + Q <--> EQ", 1, 0)[1:5] == key(at_free_E, 2, 0)[1:5]
        @test key("E + Q <--> EQ", 1, 0)[6:7] == (1, 5)
        @test key(at_free_E, 2, 0)[6:7] == (1, 2)
        @test key("E + Q <--> EQ", 1, 0) > key(at_free_E, 2, 0)
        # then, between two steps of the same forms, the indices of their metabolites,
        # stored as in metabolites(m), the later first: in the steady-state uni-uni the
        # uptake of P into ES, which runs against the reaction, goes before the release of
        # S
        uni = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                ES <--> E + P
            end
        end
        p_uptake = _testhelper_key(uni, "E + P <--> ES", 1, 0)
        s_release = _testhelper_key(uni, "E + S <--> ES", 2, 0)
        @test p_uptake[1:7] == s_release[1:7]
        @test (p_uptake[8], s_release[8]) == ([2], [1]) && p_uptake > s_release
        # the stored place comes last
        @test key(at_free_E, 1, 0)[end] == _testhelper_edge(m, at_free_E)
    end

    @testset "a spec written in another order of forms, steps and groups is itself" begin
        # Each spec written with its steps reversed, the forms after free E reversed and its
        # kinetic groups numbered the other way round, each group's allosteric state
        # moving with it, is the same mechanism: the constructors store every writing in
        # one canonical order. The solve reads that stored order (`_priority_key` reads
        # the stored indices), so the fitted constants do not follow the writing.
        for spec in MECHANISM_TEST_SPECS
            m = spec.mechanism
            cm = m.catalytic
            nf = length(cm.enzyme_forms)
            G = maximum(e -> e.kinetic_group, cm.steps)
            form(f) = f == 1 ? 1 : nf + 2 - f
            steps = [EnzymeRates.Mechanisms.Step(form(e.from), form(e.to), e.taken_up,
                                                 e.given_off, e.rapid_equilibrium,
                                                 G + 1 - e.kinetic_group)
                     for e in reverse(cm.steps)]
            reordered = EnzymeRates.Mechanisms.CatalyticScheme(
                cm.metabolites, [cm.enzyme_forms[1]; reverse(cm.enzyme_forms[2:end])],
                cm.substrates, cm.products, steps)
            x = _testhelper_is_allosteric(m) ?
                EnzymeRates.Mechanism(
                    reordered,
                    EnzymeRates.Mechanisms.Allostery(
                        reverse(m.allosteric.allosteric_states),
                        m.allosteric.oligomeric_state, m.allosteric.regulatory_sites)) :
                _testhelper_plain_of(reordered)
            @test x == m
        end
    end

    @testset "a step on no cycle runs away from free E (C9)" begin
        # D <--> E lies on no cycle, so it runs away from free E, however it is written:
        # its forward rate constant is k_E_to_D.
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                ES <--> EP
                EP <--> E + P
                D <--> E
            end
        end
        iso = only(k for k in eachindex(m.catalytic.steps)
                   if occursin("D", _testhelper_edge_text(m, k)))
        @test _testhelper_edge_text(m, iso) == "E <--> D"
        kf, kr = _testhelper_step_parameters(m, iso)
        @test EnzymeRates.Derivation._name(kf, m) === :k_E_to_D
        @test EnzymeRates.Derivation._name(kr, m) === :k_D_to_E
        # Written E <--> D, the step is stored the same way.
        @test (@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                ES <--> EP
                EP <--> E + P
                E <--> D
            end
        end) == m
        _testhelper_check_stored_orientation(m)
        # At rapid equilibrium the step's one constant is the equilibrium constant of that
        # direction, [D]/[E].
        re = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                ES <--> EP
                EP <--> E + P
                D ⇌ E
            end
        end
        @test parameters(re) == (:K_E_to_D, :k_EP_to_E_P, :k_ES_to_EP,
                                                :k_ES_to_E_S, :k_E_P_to_EP, :k_E_S_to_ES)
        _testhelper_check_law_against_oracle(re)
        _testhelper_check_stored_orientation(re)
    end

    @testset "a kinetic group turns with its member nearest free E (C9)" begin
        # E <--> D and ES <--> DS share their constants. The reaction turns over
        # through them in opposite directions, E → D (D takes up S) and DS → ES, so on
        # its own ES <--> DS would run from DS. The group turns as one with its lead,
        # E → D, the member nearest free E, and ES <--> DS keeps the direction it is
        # written in relative to it. So both steps run from the E side to the D side.
        tied_mm = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                ES <--> EP
                EP <--> E + P
                D + S <--> DS
                (D <--> E, DS <--> ES)
            end
        end
        # Groups and steps in canonical order: the steps at free E first, then by form.
        @test [_testhelper_edge_text(tied_mm, k)
               for k in EnzymeRates.Derivation._flat_edges(tied_mm.catalytic)] ==
              ["E <--> D", "ES <--> DS", "E + P <--> EP", "E + S <--> ES", "D + S <--> DS",
               "ES <--> EP"]
        # The conformation a form's name gives: its first letter.
        label(f) = first(String(tied_mm.catalytic.enzyme_forms[f]))
        tied =
            only(ks for ks in _testhelper_group_members(tied_mm.catalytic)
                 if length(ks) == 2)
        @test all(k -> label(tied_mm.catalytic.steps[k].from) === 'E' &&
                       label(tied_mm.catalytic.steps[k].to) === 'D', tied)
        # Written E <--> D and ES <--> DS, the group is stored the same way.
        @test (@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                ES <--> EP
                EP <--> E + P
                D + S <--> DS
                (E <--> D, ES <--> DS)
            end
        end) == tied_mm
        _testhelper_check_stored_orientation(tied_mm)
        _testhelper_check_law_against_oracle(tied_mm)
    end

    @testset "fundamental cycles carry their multiple of the net reaction" begin
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A ⇌ EA
                E + B ⇌ EB
                EA + B ⇌ EAB
                EB + A ⇌ EAB
                EAB <--> EP
                EP ⇌ E + P
            end
        end
        Z, multiples = EnzymeRates.Derivation._fundamental_cycles(
            m.catalytic, EnzymeRates.Derivation._flat_edges(m.catalytic))
        @test size(Z) == (2, 6)
        @test sort(abs.(multiples)) == [0, 1]
        # A cycle that turns A into P without B and Q is no multiple of the net reaction.
        # The guard rejects such a mechanism, so the macro would throw before the solve:
        # the check is an internal assertion, reached here through the unguarded build.
        slip = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> EA
                EA + B <--> EAB
                EAB <--> EPQ
                EPQ <--> EQ + P
                EQ <--> E + Q
                EA <--> EP
                EP <--> E + P
            end
        end
        err = @test_throws ErrorException EnzymeRates.Derivation._fundamental_cycles(
            slip.catalytic, EnzymeRates.Derivation._flat_edges(slip.catalytic))
        @test occursin("produces metabolite change not proportional to net reaction",
                       err.value.msg)
        @test !EnzymeRates.Derivation.check_valid_rate_equation(Bool, slip)
        # A cycle that takes up a modifier and gives it off nowhere, E + I ⇌ EI then
        # EI <--> E, is no multiple of the net reaction either, which runs no I: the
        # guard rejects it (its message is pinned in
        # test_rate_equation_correctness_guard.jl), and the assertion reads every
        # metabolite.
        consumes_i = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: I
            steps: begin
                E + S ⇌ ES
                ES <--> E + P
                E + I ⇌ EI
                EI <--> E
            end
        end
        @test !EnzymeRates.Derivation.check_valid_rate_equation(Bool, consumes_i)
        err = @test_throws ErrorException EnzymeRates.Derivation._fundamental_cycles(
            consumes_i.catalytic, EnzymeRates.Derivation._flat_edges(consumes_i.catalytic))
        @test occursin("produces metabolite change not proportional to net reaction",
                       err.value.msg)
    end

    @testset "the reduction runs in Int64 and retries in BigInt on overflow" begin
        M, pivots = EnzymeRates.Derivation._reduce_joint_system([1 -1; 1 1], [1, 0])
        @test eltype(M) == Rational{Int} && pivots == [1, 2]
        A = [2^40 1; 3 2^40]
        M, pivots = EnzymeRates.Derivation._reduce_joint_system(A, [1, 0])
        @test eltype(M) == Rational{BigInt} && pivots == [1, 2]
        @test Rational{BigInt}.(A) * M[:, 3] == [1, 0]
        err = @test_throws ErrorException EnzymeRates.Derivation._reduce_joint_system(
            reshape([1, 1], 2, 1), [0, 1])
        @test occursin("Thermodynamically contradictory mechanism", err.value.msg)
    end

    @testset "a class keeps its last column, a member of the inactive law first" begin
        # I binds EA and EP under one :NonequalAI group. In the inactive
        # conformation EA and EAI lie behind the :OnlyA binding of A, so only the
        # EP + I member is in that conformation's law.
        am = @enzyme_mechanism begin
            substrates: A
            products: P
            modifiers: I
            oligomeric_state: 1
            steps: begin
                E + A <--> EA                                   :: OnlyA
                EA <--> EP                                    :: OnlyA
                E + P ⇌ EP                                      :: EqualAI
                (EA + I ⇌ EAI, EP + I ⇌ EIP)          :: NonequalAI
            end
        end
        kept = EnzymeRates.Derivation._kept_members(am)
        group = [e.kinetic_group for e in am.catalytic.steps]
        g = only(g for (g, members) in enumerate(
            _testhelper_group_members(am.catalytic)) if length(members) == 2)
        # Active: EA + I sorts before EP + I (A before P), so it is kept; inactive:
        # the member in the law is kept.
        @test _testhelper_edge_text(am, kept[(g, 1, false)]) == "EA + I ⇌ EAI"
        @test _testhelper_edge_text(am, kept[(g, 1, true)]) == "EP + I ⇌ EIP"
        sys = EnzymeRates.Derivation._joint_system(am)
        for (key, edge) in kept
            last_col = findlast(c -> c.slot == key[2] && c.inactive == key[3] &&
                                     group[c.step] == key[1], sys.columns)
            @test sys.columns[last_col].step == edge
        end
        _testhelper_check_law_against_oracle(am)
    end

    @testset "the inactive cycle row runs over the whole graph, islands included" begin
        # EA and EAI lie behind the :OnlyA binding of A, so the inactive law holds
        # only E, EP and EIP. The cycle through EA still gets an inactive row,
        # which reads the inactive copies of the two :OnlyA steps, outside the law.
        am = @enzyme_mechanism begin
            substrates: A
            products: P
            modifiers: I
            oligomeric_state: 1
            steps: begin
                E + A <--> EA                                   :: OnlyA
                EA <--> EP                                    :: OnlyA
                E + P ⇌ EP                                      :: EqualAI
                (EA + I ⇌ EAI, EP + I ⇌ EIP)          :: NonequalAI
            end
        end
        sys = EnzymeRates.Derivation._joint_system(am)
        # The cycle once per conformation, then the second member of the :NonequalAI
        # group's active class and of its inactive class.
        @test size(sys.A) == (4, 13)
        # Row 1 is the active cycle, row 2 the inactive one. The :EqualAI binding of P has
        # no inactive constant, so the inactive row reads its active one.
        reads = Set((_testhelper_edge_text(am, c.step), c.slot, c.inactive)
                    for c in sys.columns[findall(!iszero, sys.A[2, :])])
        @test reads == Set([("E + A <--> EA", 1, true), ("E + A <--> EA", 2, true),
                            ("EA <--> EP", 1, true), ("EA <--> EP", 2, true),
                            ("E + P ⇌ EP", 1, false)])
    end

    @testset "every spec's solve counts its declared dimension" begin
        for spec in MECHANISM_TEST_SPECS
            @testset "$(spec.name)" begin
                mech = spec.mechanism
                sys = EnzymeRates.Derivation._joint_system(mech)
                law = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(
                    mech).inactive
                _, pivots = EnzymeRates.Derivation._reduce_joint_system(sys.A, sys.rhs)
                excluded = (EnzymeRates.Derivation._excluded_steps(mech, 1),
                            EnzymeRates.Derivation._excluded_steps(mech, 2))
                class(c) = EnzymeRates.Derivation._class(mech, excluded, c)
                last_of = Dict(class(c) => i for (i, c) in enumerate(sys.columns))
                n_law = count(i -> (c = sys.columns[i]; (!c.inactive || law[c.step]) &&
                                                        last_of[class(c)] == i &&
                                                        i ∉ pivots), eachindex(sys.columns))
                n_reg = _testhelper_is_allosteric(mech) ?
                    sum((t === :NonequalAI ? 2 : 1
                         for site in mech.allosteric.regulatory_sites
                         for t in site.allosteric_states);
                        init = 0) + 1 : 0
                @test n_law + n_reg == _testhelper_declared_dimension(spec.mechanism)
            end
        end
    end

    @testset "every spec keeps the last column of each class (E3)" begin
        for spec in MECHANISM_TEST_SPECS
            @testset "$(spec.name)" begin
                mech = spec.mechanism
                sys = EnzymeRates.Derivation._joint_system(mech)
                tags = EnzymeRates.Derivation._solve_tags(mech)
                # In elimination order the last column of a class overwrites the earlier
                # ones. The inactive copy of an :OnlyA step is a class of its own, which
                # `_kept_members` leaves out.
                group = [e.kinetic_group for e in mech.catalytic.steps]
                last_of_class = Dict{Tuple{Int, Int, Bool}, Int}()
                for c in sys.columns
                    g = group[c.step]
                    c.inactive && tags[g] === :OnlyA && continue
                    last_of_class[(g, c.slot, c.inactive)] = c.step
                end
                @test EnzymeRates.Derivation._kept_members(mech) == last_of_class
            end
        end
    end

    @testset "every spec stores its steps in the direction the reaction graph gives" begin
        # The solve reads each stored step as the forward direction of its constants, and
        # the priority key reads the step as stored. The constructor orients every step
        # from the graph and turns each kinetic group with its lead, so orienting the
        # stored steps again leaves each one as stored, and writing every step the other
        # way round, in the reverse order, builds the same mechanism.
        # The C9 testsets above and the exchange testsets of test_mechanism_type.jl check
        # the same on steps that free E's distance orients and on kinetic groups.
        for spec in MECHANISM_TEST_SPECS
            @testset "$(spec.name)" begin
                _testhelper_check_stored_orientation(spec.mechanism)
            end
        end
    end

    @testset "every spec's solve satisfies detailed balance in the oracle" begin
        for (i, spec) in enumerate(MECHANISM_TEST_SPECS)
            @testset "$(spec.name)" begin
                rng = Random.MersenneTwister(i)
                mech = spec.mechanism
                point = setprecision(BigFloat, 256) do
                    _testhelper_solve_point(mech, rng)
                end
                subs, prods = _testhelper_substrates(mech), _testhelper_products(mech)
                concs = Dict{Symbol, BigFloat}(
                    x => exp(2 * big(rand(rng)) - 1)
                    for x in EnzymeRates.metabolites(spec.mechanism))
                p1 = first(prods)
                concs[p1] = point.Keq * prod(concs[x] for x in subs) /
                            prod((concs[x] for x in prods if x != p1); init = big(1))
                off = copy(concs)
                off[p1] /= 10
                v_eq = _testhelper_steady_state_rate(mech, point, concs)
                v_off = _testhelper_steady_state_rate(mech, point, off)
                @test abs(v_eq) <= big(10)^-30 * abs(v_off)
            end
        end
    end

end

@testset "_solve's half-integer dependent has keq -1//2 and exact exponents" begin
    # The conformational square ties two isomerizations of one group, so the Haldane
    # relation fixes k_F_to_E through square roots: Keq = (k_E_to_F / k_F_to_E)² ·
    # (k_F_S_to_ES / k_ES_to_F_S) · (k_FS_to_E_P / k_E_P_to_FS).
    half = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (ES <--> FS, E <--> F)
            FS <--> E + P
            F + S <--> ES
        end
    end
    s = EnzymeRates.Derivation._solve(half)
    @test issorted([d.col for d in s.dependents])
    @test sort([d.label for d in s.dependents]) ==
          [:haldane, :kinetic_group, :kinetic_group]
    d = only(d for d in s.dependents if d.label === :haldane)
    @test s.names[d.col] === :k_F_to_E && d.keq == -1//2
    @test Dict(s.names[f] => e for (f, e) in d.factors) ==
          Dict(:k_E_to_F => 1//1, :k_E_P_to_FS => -1//2, :k_FS_to_E_P => 1//2,
               :k_ES_to_F_S => -1//2, :k_F_S_to_ES => 1//2)
    # The rendered exponents are exact `Rational{Int}` literals (`_testhelper_power_expr`).
    @test _testhelper_dependent_param_exprs(s)[1] ==
          Dict(:k_F_to_E => :((Keq ^ $(-1//2)) * (k_ES_to_F_S ^ $(-1//2)) *
                              (k_E_P_to_FS ^ $(-1//2)) * k_E_to_F *
                              (k_FS_to_E_P ^ $(1//2)) * (k_F_S_to_ES ^ $(1//2))))
end

@testset "_solve labels kinetic-group, Haldane and Wegscheider dependents" begin
    # A binds E and EB under one group (a kinetic group: its occupied member's
    # constant equals the kept one); B's two bindings close a binding square with A's,
    # which ties them (Wegscheider); the cycle through the catalytic step fixes its
    # reverse rate (Haldane).
    folded = @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB)
            E + B ⇌ EB
            EA + B ⇌ EAB
            EAB <--> EP
            EP ⇌ E + P
        end
    end
    s = EnzymeRates.Derivation._solve(folded)
    described = [(s.names[d.col], d.label, d.keq, [s.names[f] => e for (f, e) in d.factors])
                 for d in s.dependents]
    @test described == [
        (:k_EP_to_EAB, :haldane, -1//1, [:k_EAB_to_EP => 1//1, :K_EP_to_E_P => 1//1,
                                         :K_EB_to_E_B => -1//1, :K_EA_to_E_A => -1//1]),
        (:K_EA_to_E_A, :kinetic_group, 0//1, [:K_EA_to_E_A => 1//1]),
        (:K_EAB_to_EA_B, :wegscheider, 0//1, [:K_EB_to_E_B => 1//1])]
    @test s.fitted == [:K_EA_to_E_A, :K_EB_to_E_B, :K_EP_to_E_P, :k_EAB_to_EP]
    @test s.names[s.independents] == s.fitted
    @test _testhelper_dependent_param_exprs(s) ==
          (Dict(:K_EAB_to_EA_B => :K_EB_to_E_B,
                :k_EP_to_EAB => :((1 / Keq) * (1 / K_EA_to_E_A) * (1 / K_EB_to_E_B) *
                                  K_EP_to_E_P * k_EAB_to_EP)),
           (:K_EA_to_E_A, :K_EB_to_E_B, :K_EP_to_E_P, :k_EAB_to_EP))
end

@testset "_solve maps an :EqualAI edge's inactive constant to its active column" begin
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> E + P     :: NonequalAI
        end
    end
    s = EnzymeRates.Derivation._solve(am)
    C = EnzymeRates.Derivation._StepParameter
    binding = _testhelper_edge(am, "E + S ⇌ ES")
    release = _testhelper_edge(am, "E + P <--> ES")
    @test length(s.columns) == 5
    @test s.column[C(binding, 1, true)] == s.column[C(binding, 1, false)]
    @test s.columns[s.column[C(binding, 1, false)]] == C(binding, 1, false)
    for slot in 1:2, inactive in (false, true)
        @test s.columns[s.column[C(release, slot, inactive)]] == C(release, slot, inactive)
    end
    @test s.fitted == [:K_ES_to_E_S, :k_A_ES_to_E_P, :k_I_ES_to_E_P, :L]
end

@testset "_solve lists the fitted names in parameter order" begin
    for spec in MECHANISM_TEST_SPECS
        @test Tuple(
            EnzymeRates.Derivation._solve(spec.mechanism).fitted) ==
              parameters(spec.mechanism)
    end
end

# ── Pivot priority ────────────────────────────────────────────────────────────────────

module PivotPriorityRegressionTests
using Test, EnzymeRates

@testset "pivot priority: non-allo Wegscheider tie enforced (no over-count)" begin
    # A non-allosteric LDH mechanism reached in enumeration. A closed all-RE cycle
    # — E→ELactate→ELactateNAD and E→ENAD→ELactateNAD with shared K_ENAD_to_E_NAD —
    # forces the Wegscheider tie K_ELactate_to_E_Lactate =
    # K_ELactateNADH_to_ENADH_Lactate. The solve must make one of the two dependent,
    # or parameters(m) over-counts by 1 (7 instead of 6). A split that recovers the
    # tie would then be a -1 parameter-count edge.
    m = @enzyme_mechanism begin
        substrates: NADH, Pyruvate
        products: Lactate, NAD
        steps: begin
            ELactate ⇌ E + Lactate
            (ENAD ⇌ E + NAD, ELactateNAD ⇌ ELactate + NAD)
            (E + NADH <--> ENADH, ELactate + NADH <--> ELactateNADH)
            (ELactateNAD ⇌ ENAD + Lactate, ELactateNADH ⇌ ENADH + Lactate)
            (ENAD + Pyruvate ⇌ ENADPyruvate,
             ENADH + Pyruvate ⇌ ENADHPyruvate)
            ENADHPyruvate <--> ELactateNAD
        end
    end
    fp = parameters(m)
    @test length(fp) == 6
    # the tie collapses one
    @test !(:K_ELactate_to_E_Lactate in fp && :K_ELactateNADH_to_ENADH_Lactate in fp)
end
end # module

# ── Collapse of :EqualAI constants ────────────────────────────────────────────────────

module AllostericCollapseTests
using Test, EnzymeRates, Random
const ER = EnzymeRates
# A step written from the form named `f` to the form named `t`, taking up the metabolites
# named `up` and giving off those named `off`, at rapid equilibrium when `re` holds.
_testhelper_step(f, t, up, off, re) = (; f, t, up, off, re)
# The allosteric mechanism over substrates `subs` and products `prods`, catalytic
# multiplicity 2, whose steps are `steps`, step `j` alone in group `j`, tagged `states[j]`.
function _testhelper_one_step_groups(subs, prods, steps, states)
    mets = [subs; prods]
    forms = unique([n for s in steps for n in (s.f, s.t)])
    form(f) = findfirst(==(f), forms)
    indices(xs) = Int[findfirst(==(x), mets) for x in xs]
    ER.Mechanism(
        ER.Mechanisms.CatalyticScheme(
            mets, forms, collect(eachindex(subs)),
            collect(length(subs) .+ eachindex(prods)),
            [ER.Mechanisms.Step(form(s.f), form(s.t), indices(s.up), indices(s.off), s.re,
                                j)
             for (j, s) in enumerate(steps)]),
        ER.Mechanisms.Allostery(collect(Symbol, states), 2,
                                ER.Mechanisms.RegulatorySite[]))
end
# `re` flags which of the S-binding, chemistry and P-release steps are rapid
# equilibrium; a steady-state binding carries kon/koff.
function _testhelper_uni(states; re=(true,false,true))
    none = Symbol[]
    steps=[_testhelper_step(:E,:ES,[:S],none,re[1]),
           _testhelper_step(:ES,:EP,none,none,re[2]),
           _testhelper_step(:EP,:E,none,[:P],re[3])]
    _testhelper_one_step_groups([:S], [:P], steps, states)
end
# Random-order bi-bi with two Wegscheider boxes; `states` tags its nine one-step groups.
function _testhelper_ro_bibi(states)
    none = Symbol[]
    sd=[_testhelper_step(:E,:EA,[:A],none,true),_testhelper_step(:E,:EB,[:B],none,true),
        _testhelper_step(:EB,:EAB,[:A],none,true),_testhelper_step(:EA,:EAB,[:B],none,true),
        _testhelper_step(:EAB,:EPQ,none,none,false),
        _testhelper_step(:EP,:EPQ,[:Q],none,true),
        _testhelper_step(:EQ,:EPQ,[:P],none,true),_testhelper_step(:E,:EP,[:P],none,true),
        _testhelper_step(:E,:EQ,[:Q],none,true)]
    _testhelper_one_step_groups([:A, :B], [:P, :Q], sd, states)
end
function _testhelper_evalrate(am; seed=1, split=nothing)
    cem=am; rng=MersenneTwister(seed)
    rate=ER.compile_rate_equation(cem)
    fp=parameters(am)
    mets=collect(ER.metabolites(cem))
    base=[(k===:L ? 0.6 : 0.4+2rand(rng)) for k in fp]
    split!==nothing && (base=[(fp[i]===split[1] ? split[2] : base[i]) for i in 1:length(fp)])
    prm=NamedTuple{fp}(Tuple(base))
    c=NamedTuple{Tuple(mets)}(ntuple(i->0.4+2rand(rng),length(mets)))
    v=real(rate(c,prm,3.0))
    ec=NamedTuple{Tuple(mets)}(ntuple(i->(mets[i]===:P ? 3.0 : 1.0),length(mets)))
    veq=real(rate(ec,prm,3.0))
    (fp, v, veq)
end

@testset "strict :EqualAI collapse" begin
    # baseline
    fp0,_,_ = _testhelper_evalrate(_testhelper_uni([:EqualAI,:EqualAI,:EqualAI]))

    @testset "single NonequalAI binding + EqualAI catalysis -> full collapse" begin
        fp,v,veq = _testhelper_evalrate(_testhelper_uni([:NonequalAI,:EqualAI,:EqualAI]))
        @test isfinite(v); @test abs(veq) < 1e-8
        @test !(:K_I_ES_to_E_S in fp)           # I-twin dropped (collapsed to a mirror)
        s = ER.rate_equation_string(_testhelper_uni([:NonequalAI,:EqualAI,:EqualAI]))
        @test occursin("K_I_ES_to_E_S=K_A_ES_to_E_S", replace(s," "=>""))  # explicit mirror
        @test !occursin("k_I_", s)              # catalysis not silently un-shared
    end

    @testset "two NonequalAI bindings + EqualAI catalysis -> 1 honorable DOF" begin
        am = _testhelper_uni([:NonequalAI,:EqualAI,:NonequalAI])
        fp,v,veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        # exactly one free I-split survives; the other is a derived mirror.
        nI = count(p->startswith(String(p),"K_I_"), fp)
        @test nI == 1
        s = replace(ER.rate_equation_string(am), " "=>"")
        @test occursin("k_I_", s) == false                    # catalysis stays shared
        # the surviving split moves the rate (identifiable)
        freeI = fp[findfirst(p->startswith(String(p),"K_I_"), fp)]
        v1 = _testhelper_evalrate(am; split=(freeI,1.3))[2]
        v2 = _testhelper_evalrate(am; split=(freeI,5.0))[2]
        @test !isapprox(v1, v2)
    end

    @testset "catalysis NonequalAI -> native, no collapse, no mirror" begin
        am = _testhelper_uni([:NonequalAI,:NonequalAI,:EqualAI])
        fp,v,veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        @test :K_I_ES_to_E_S in fp              # binding split free
        @test any(p->startswith(String(p),"k_I_"), fp)  # catalysis split free (native)
    end

    @testset "binding-Wegscheider single inner edge -> full collapse" begin
        # random-order bi-bi (two Wegscheider boxes); tag ONE inner box-independent
        # edge (EB+A->EAB) :NonequalAI, rest :EqualAI -> its split is forbidden ->
        # collapses.
        st=fill(:EqualAI,9); st[3]=:NonequalAI
        am=_testhelper_ro_bibi(st)
        fp,v,veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8             # thermo-consistent
        @test !(:K_I_EAB_to_EB_A in fp)                     # forbidden split collapsed
        dep,_=Main._testhelper_dependent_param_exprs(am)
        @test dep[:K_I_EAB_to_EB_A] == dep[:K_A_EAB_to_EB_A]  # one value in both states
    end

    # ── Steady-state binding: affinity/speed decomposition (Option 3) ──
    # A steady-state :NonequalAI binding carries two rate constants: its affinity
    # (kon/koff) can be forbidden and collapse (deriving the reverse), while its
    # speed (the forward kon) stays free. Only the affinity collapses — never the
    # whole binding. (The steady-state Wegscheider-box case is covered by `m_ro`
    # in test_rate_equation_derivation.jl.)
    @testset "SS binding + EqualAI catalysis -> affinity collapses, speed free" begin
        am = _testhelper_uni([:NonequalAI,:EqualAI,:EqualAI]; re=(false,false,false))
        fp,v,veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        @test !(:k_I_ES_to_E_S in fp)           # affinity collapsed: reverse derived
        @test :k_I_E_S_to_ES in fp              # speed (forward) stays free
        s = replace(ER.rate_equation_string(am), " "=>"")
        @test occursin("k_I_ES_to_E_S=", s)     # explicit reverse-rate mirror
        # the surviving speed split moves the rate (identifiable)
        v1 = _testhelper_evalrate(am; split=(:k_I_E_S_to_ES,1.3))[2]
        v2 = _testhelper_evalrate(am; split=(:k_I_E_S_to_ES,5.0))[2]
        @test !isapprox(v1, v2)
    end

    @testset "SS binding + NonequalAI catalysis -> not forbidden, stays free" begin
        am = _testhelper_uni([:NonequalAI,:NonequalAI,:EqualAI]; re=(false,false,false))
        fp,v,veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        # both free (affinity honorable)
        @test (:k_I_E_S_to_ES in fp) && (:k_I_ES_to_E_S in fp)
    end

    # ── Mixed-type + Wegscheider-coupled :NonequalAI bindings. Two mechanisms that
    #    were silently mis-derived before the combined solve: the constraint matrix
    #    must use one uniform binding-K sign, and an all-RE Wegscheider box with
    #    :NonequalAI edges must derive its tied binding-K with no undefined symbol.
    @testset "mixed RE/SS coupled :NonequalAI bindings -> consistent (uniform sign)" begin
        # S-binding RE, P-release SS, both :NonequalAI, catalysis :EqualAI: the two
        # coupled affinities are *different* step types, so the constraint matrix must
        # use one uniform sign — a per-type flip inverts the coupling → nonzero flux.
        am=_testhelper_uni([:NonequalAI,:EqualAI,:NonequalAI]; re=(true,false,false))
        fp,v,veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
    end

    @testset ":NonequalAI RE Wegscheider box binding -> consistent, no undefined symbol" begin
        # random-order bi-bi, all-RE bindings, SS :EqualAI catalysis; tag two box
        # edges :NonequalAI so one binding-K is Wegscheider-derived from the other.
        # The combined solve must express every derived symbol purely in free columns
        # (no circular reference ⇒ no UndefVarError) and hold detailed balance.
        st=fill(:EqualAI,9); st[2]=:NonequalAI; st[3]=:NonequalAI
        # must not throw UndefVarError
        fp,v,veq = _testhelper_evalrate(_testhelper_ro_bibi(st))
        @test isfinite(v); @test abs(veq) < 1e-8
    end

    @testset "dead-I NonequalAI binding -> K_I identifiable, NOT collapsed" begin
        # I state cannot turn over (OnlyA catalysis) but binds S with its own
        # affinity: K_A_ES_to_E_S and K_I_ES_to_E_S are BOTH identifiable (a dead-end
        # E_I·S is in no cycle, so nothing pins K_I to K_A). HEAD over-collapses this.
        am = _testhelper_uni([:NonequalAI, :OnlyA, :EqualAI])
        fp, v, veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        @test :K_I_ES_to_E_S in fp                    # NOT collapsed
        v1 = _testhelper_evalrate(am; split=(:K_I_ES_to_E_S, 1.3))[2]
        v2 = _testhelper_evalrate(am; split=(:K_I_ES_to_E_S, 5.0))[2]
        @test !isapprox(v1, v2)                       # identifiable (moves the rate)
        s = replace(ER.rate_equation_string(am), " " => "")
        @test !occursin("K_I_ES_to_E_S=K_A_ES_to_E_S", s)  # no collapse mirror
    end
end
end # module

# ── parameters and dependent_parameter_values ─────────────────────────────────────────

@testset "parameters API symmetry" begin
    # A bi-bi seed whose catalytic cycle's Haldane relation makes a constant dependent;
    # B and P also bind as dead ends, at EQ and EA.
    m_mono = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            E + Q ⇌ EQ
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            (EA + P ⇌ EAP, EQ + P ⇌ EPQ)
            EAB <--> EPQ
        end
    end
    params_mono = parameters(m_mono)
    @test isdisjoint(params_mono, (:Keq, :E_total))
    # dependent_parameter_values(m, values, Keq) is exactly the dependent constants to the
    # parameters.
    values_mono = NamedTuple{params_mono}(ntuple(_ -> 1.0, length(params_mono)))
    dep_mono, _ = _testhelper_dependent_param_exprs(m_mono)
    @test !isempty(dep_mono)
    @test Set(keys(dependent_parameter_values(m_mono, values_mono, 1.0))) ==
          Set(keys(dep_mono))
    @test isdisjoint(keys(dependent_parameter_values(m_mono, values_mono, 1.0)),
                     params_mono)

    base = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    allosteric_states = fill(:NonequalAI, length(_testhelper_group_members(base.catalytic)))
    site = EnzymeRates.Mechanisms.RegulatorySite([:R], 2, [:NonequalAI])
    am = EnzymeRates.Mechanism(
        base.catalytic, EnzymeRates.Mechanisms.Allostery(allosteric_states, 2, [site]))
    @test am.allosteric.allosteric_states == allosteric_states
    params_allo = parameters(am)
    @test any(occursin("_I_", string(p)) for p in params_allo)
    @test any(occursin("reg", string(p)) for p in params_allo)
    @test :L in params_allo
    @test isdisjoint(params_allo, (:Keq, :E_total))
    values_allo = NamedTuple{params_allo}(ntuple(_ -> 1.0, length(params_allo)))
    dep_allo, _ = _testhelper_dependent_param_exprs(am)
    @test !isempty(dep_allo)
    @test Set(keys(dependent_parameter_values(am, values_allo, 1.0))) ==
          Set(keys(dep_allo))
    @test isdisjoint(keys(dependent_parameter_values(am, values_allo, 1.0)),
                     params_allo)
end

@testset "fitted parameters of a mechanism equal its compiled rate law's" begin
    for spec in MECHANISM_TEST_SPECS
        @test parameters(spec.mechanism) ===
              parameters(compile_rate_equation(spec.mechanism))
    end
end

@testset "dependent_parameter_values returns the dependent constants" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    p = (K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0)
    deps = dependent_parameter_values(m, p, 4.0)
    # The Haldane relation: k_EP_to_ES = K_EP_to_E_P · k_ES_to_EP / (K_ES_to_E_S · Keq).
    @test deps == (k_EP_to_ES = 1.0 * 3.0 / (2.0 * 4.0),)
    @test keys(deps) == (:k_EP_to_ES,)
    # A dependent constant already in params is ignored and recomputed.
    @test dependent_parameter_values(m, merge(p, (k_EP_to_ES = 99.0,)), 4.0) == deps
    @test dependent_parameter_values(m, merge((k_EP_to_ES = 99.0,), p), 4.0) == deps
    # The concrete mechanism gives the same values, the argument Keq wins over a Keq in
    # params, and BigFloat values stay BigFloat.
    @test dependent_parameter_values(m, p, 4.0) == deps
    @test dependent_parameter_values(m, merge(p, (Keq = 99.0, E_total = 99.0)), 4.0) ==
          deps
    @test dependent_parameter_values(m, map(big, p), big(4.0)).k_EP_to_ES == big(3) / 8
    @test dependent_parameter_values(m, map(big, p), big(4.0)).k_EP_to_ES isa BigFloat
    @test_throws ErrorException(
        "dependent_parameter_values: params must hold every name of " *
        "parameters(m); it lacks K_ES_to_E_S",
    ) dependent_parameter_values(m, (K_EP_to_E_P = 1.0, k_ES_to_EP = 3.0), 4.0)
end

@testset "dependent_parameter_values reproduces the dependents of a consistent point" begin
    # The point's constants come from free energies, independently of the solve; the
    # fitted ones and Keq determine the rest. `half` has a dependent with
    # half-integer exponents. `cube` has one with third-integer exponents (k_F_to_E
    # through Keq^-1/3), which a Float64 exponent would miss by about 1e-16, far outside
    # the 1e-40 tolerance.
    half = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (ES <--> FS, E <--> F)
            FS <--> E + P
            F + S <--> ES
        end
    end
    cube = @enzyme_mechanism begin
        substrates: S, B
        products: P
        steps: begin
            (E <--> F, ES <--> FS, EP <--> FP)
            F + S <--> ES
            FS + B <--> EP
            FP <--> E + P
        end
    end
    rng = Random.MersenneTwister(12)
    for m in [[s.mechanism for s in MECHANISM_TEST_SPECS]; half; cube]
        point = setprecision(() -> _testhelper_consistent_point(m, rng), BigFloat, 256)
        fitted = parameters(m)
        p = NamedTuple{fitted}(getindex.(Ref(point.vals), fitted))
        deps = setprecision(() -> dependent_parameter_values(m, p, point.Keq), BigFloat,
                            256)
        all_p = merge(p, deps)
        @test all(isapprox(all_p[k], point.vals[k]; rtol = big(10)^-40)
                  for k in keys(all_p))
        @test issorted(collect(keys(all_p))[length(p) + 1:end])
        # Every step constant the law reads, in each conformation, is among the keys; the
        # regulator constants and L are fitted, so p holds them.
        states = _testhelper_is_allosteric(m) ? (:A, :I) : (:None,)
        @test issubset((n for st in states
                        for ns in _testhelper_edge_names(m, st) for n in ns),
                       keys(all_p))
    end
end

@testset "dependent_parameter_values computes rational powers of Keq from the solve" begin
    # The dependent of this conformational square is k_F_to_E = Keq^-1/2 ·
    # k_ES_to_F_S^-1/2 · k_E_P_to_FS^-1/2 · k_E_to_F · k_FS_to_E_P^1/2 · k_F_S_to_ES^1/2
    # (test_thermodynamic_constraints.jl): 1/4 at Keq = k_ES_to_F_S = 4 and every other
    # constant 1, and 1 once k_F_S_to_ES is 16. BigFloat inputs give the exact value.
    half = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (ES <--> FS, E <--> F)
            FS <--> E + P
            F + S <--> ES
        end
    end
    p = (k_E_to_F = 1.0, k_E_P_to_FS = 1.0, k_FS_to_E_P = 1.0, k_ES_to_F_S = 4.0,
         k_F_S_to_ES = 1.0)
    @test dependent_parameter_values(half, p, 4.0).k_F_to_E == 0.25
    @test dependent_parameter_values(half, merge(p, (k_F_S_to_ES = 16.0,)), 4.0).k_F_to_E ==
          1.0
    @test dependent_parameter_values(half, map(big, p), big(4.0)).k_F_to_E == big(0.25)
end

# ── Constant names of inhibitor complexes ─────────────────────────────────────────────

@testset "EI and EIinh, EPmod and EPinh name their constants apart" begin
    # EI and EIinh render EI and EIinh: the mechanism is accepted, and each constant
    # takes a name of its own.
    iinh = @enzyme_mechanism begin
        substrates: Iinh
        products: P
        modifiers: I
        steps: begin
            E + Iinh ⇌ EIinh
            EIinh <--> EP
            EP ⇌ E + P
            E + I ⇌ EI
        end
    end
    @test _testhelper_step_constant_names(iinh) ==
          Set([:K_EIinh_to_E_Iinh, :k_EIinh_to_EP, :k_EP_to_EIinh, :K_EP_to_E_P,
               :K_EI_to_E_I])
    _testhelper_check_law_against_oracle(iinh)
    # So with a reactant's modifier copy: EPmod and EPinh render EPmod and EPinh.
    pinh = @enzyme_mechanism begin
        substrates: Pinh
        products: P
        steps: begin
            E + Pinh ⇌ EPinh
            EPinh <--> EP
            EP ⇌ E + P
            E + P ⇌ EPmod
        end
    end
    @test _testhelper_step_constant_names(pinh) ==
          Set([:K_EPinh_to_E_Pinh, :k_EPinh_to_EP, :k_EP_to_EPinh, :K_EP_to_E_P,
               :K_EPmod_to_E_P])
    _testhelper_check_law_against_oracle(pinh)
end

# ── Fitted and dependent constants of single mechanisms ───────────────────────────────

@testset "a fused release keeps its catalytic rate constant fitted" begin
    # ES <--> E + P gives off P while turning S into P. It is stored as the binding it
    # reverses, E + P → ES, whose forward constant k_E_P_to_ES runs against the
    # reaction. The Haldane relation makes that constant dependent, as k₋₂ of the
    # textbook mechanism, and the catalytic constant k_ES_to_E_P stays fitted.
    ss = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S <--> ES
            ES <--> E + P
        end
    end
    @test Set(parameters(ss)) ==
          Set([:k_E_S_to_ES, :k_ES_to_E_S, :k_ES_to_E_P])
    re = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
        end
    end
    @test Set(parameters(re)) == Set([:K_ES_to_E_S, :k_ES_to_E_P])
end

function _testhelper_regulator_names(am)
    Symbol[EnzymeRates.Derivation._name(p, am)
           for p in EnzymeRates.Derivation._regulator_parameters(am)]
end

# ── Degenerate constraint handling ────────────────────────────────────────────

@testset "Degenerate constraint handling" begin
    # ── Unit tests: a rendered dependent's types ─────────────────

    @testset "a dependent renders as a Symbol, an Expr or 1.0 with exact exponents" begin
        render = _testhelper_power_expr
        R = Rational{BigInt}
        # Single symbol factor (exp=1): returns bare Symbol
        @test render(R(0), [(:k1f, R(1))]) === :k1f

        # Single inverse factor: returns Expr
        @test render(R(0), [(:k1f, R(-1))]) isa Expr

        # Keq-only: returns Symbol
        @test render(R(1), Tuple{Symbol, R}[]) === :Keq

        # Multiple factors: returns Expr
        @test render(R(0), [(:k1f, R(1)), (:k2f, R(1))]) isa Expr

        # No factors and zero Keq: returns the Float64 literal 1.0
        @test render(R(0), Tuple{Symbol, R}[]) === 1.0

        # Keq times a factor is a valid AST node
        @test render(R(1), [(:k1f, R(1))]) isa Union{Float64, Symbol, Expr}

        # A fractional exponent is an exact Rational{Int} literal, so the expression
        # evaluates exactly in BigFloat.
        @test render(R(-1//2), [(:k1f, R(1//2))]) ==
              Expr(:call, :*, Expr(:call, :^, :Keq, -1//2), Expr(:call, :^, :k1f, 1//2))
        @test render(R(-1//2), [(:k1f, R(1//2))]).args[2].args[3] isa Rational{Int}
        half = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                (ES <--> FS, E <--> F)
                FS <--> E + P
                F + S <--> ES
            end
        end
        dep, fitted = _testhelper_dependent_param_exprs(half)
        setprecision(BigFloat, 256) do
            p = NamedTuple{Tuple(fitted)}(
                Tuple(big(1) + big(i) / 7 for i in 1:length(fitted)))
            Keq = big(1) + big(length(fitted) + 1) / 7
            env = Dict{Symbol, BigFloat}(k => v for (k, v) in pairs(merge(p, (; Keq))))
            k_F_to_E = dependent_parameter_values(half, p, Keq).k_F_to_E
            @test abs(_testhelper_printed_value(dep[:k_F_to_E], env) - k_F_to_E) <=
                  big(10)^-60 * k_F_to_E
        end
    end

    @testset "a constant thermodynamics pins to 1 prints as 1" begin
        # I binds E and F under one dissociation constant, so the cycle E → EI → F → E
        # forces K_E_to_F = 1 and leaves no other constant in its row.
        m = @enzyme_mechanism begin
            substrates: A
            products: P
            modifiers: I
            steps: begin
                E + A <--> EA
                EA <--> EP
                EP <--> E + P
                E ⇌ F
                (E + I ⇌ EI, F + I ⇌ EI)
            end
        end
        @test occursin("# Wegscheider constraints:\nK_E_to_F = 1\n",
                       rate_equation_string(m))
        @test :K_E_to_F ∉ parameters(m)
        _testhelper_check_law_against_oracle(m)
        # The generated rate equation folds the pinned constant into its coefficients
        # and stays allocation-free.
        syms = parameters(m)
        params = NamedTuple{syms}(ntuple(i -> 0.5 + i / 10, length(syms)))
        @test dependent_parameter_values(m, params, 1.5).K_E_to_F == 1
        allocs, _ = _testhelper_test_rate_equation_performance(
            compile_rate_equation(m), params, 1.5, (A = 1.0, P = 0.5, I = 0.2))
        @test allocs == 0
    end

end

@testset "Parameter-struct allosteric helpers" begin
    @testset "_testhelper_edge_names and the regulator constants" begin
        flat_names(am, state) = reduce(vcat, _testhelper_edge_names(am, state))
        # :NonequalAI cat groups + :NonequalAI reg ligand → both conformations named.
        am = @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: NonequalAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: NonequalAI
            end
            regulatory_site(1): R::NonequalAI
        end
        # Catalytic side: every step of the inactive law takes its inactive name, the
        # :EqualAI iso its one shared name.
        @test Set(flat_names(am, :I)) ==
              Set([:K_I_ES_to_E_S, :k_ES_to_EP, :k_EP_to_ES, :K_I_EP_to_E_P])
        @test Set(flat_names(am, :A)) ==
              Set([:K_A_ES_to_E_S, :k_ES_to_EP, :k_EP_to_ES, :K_A_EP_to_E_P])
        @test _testhelper_regulator_constants(am, :A) == [:K_A_Rreg]
        @test _testhelper_regulator_constants(am, :I) == [:K_I_Rreg]
        @test _testhelper_regulator_names(am) == [:K_A_Rreg, :K_I_Rreg]

        # Balanced :OnlyA bindings (both S and P) keep the :NonequalAI iso
        # Haldane-valid. The inactive conformation keeps no form but E, so its law has
        # no step and the :OnlyA reg ligand no inactive constant.
        am_skip = @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: OnlyA
                ES <--> EP    :: NonequalAI
                EP ⇌ E + P      :: OnlyA
            end
            regulatory_site(1): R::OnlyA
        end
        @test all(isempty, _testhelper_edge_names(am_skip, :I))
        @test isempty(_testhelper_regulator_constants(am_skip, :I))
        @test _testhelper_regulator_names(am_skip) == [:K_A_Rreg]
        @test _testhelper_regulator_constants(am_skip, :A) ==
              [:K_A_Rreg]

        # An :OnlyI reg ligand has an inactive constant and no active one; an :EqualAI
        # ligand one bare constant, listed with the active ones.
        am_onlyi = @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: NonequalAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: NonequalAI
            end
            regulatory_site(1): R::OnlyI
        end
        @test isempty(_testhelper_regulator_constants(am_onlyi, :A))
        @test _testhelper_regulator_names(am_onlyi) == [:K_I_Rreg]
        @test _testhelper_regulator_constants(am_onlyi, :I) ==
              [:K_I_Rreg]
        am_equal = @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: NonequalAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: NonequalAI
            end
            regulatory_site(1): R::EqualAI
        end
        @test _testhelper_regulator_constants(am_equal, :A) ==
              [:K_Rreg]
        @test isempty(_testhelper_regulator_constants(am_equal, :I))
        @test _testhelper_regulator_names(am_equal) == [:K_Rreg]
    end
end

@testset "structural names + synth-dep routing: allosteric NonequalAI" begin
    # NonequalAI substrate binding (PEP), EqualAI catalysis.
    # k_cat_rev (Haldane dep) references the NonequalAI PEP binding K,
    # so an I-state dep name is produced. It must use the structural
    # mid-name I_ token, not string(active) * "_T"; a rate polynomial with
    # structural I-names beside _T-suffixed independent names would make
    # the compiled rate law error with a KeyError.
    m = @enzyme_mechanism begin
        substrates: PEP, ADP
        products:   Pyruvate, ATP
        oligomeric_state: 4
        steps: begin
            (E + PEP ⇌ EPEP,
             EADP + PEP ⇌ EADPPEP)                          :: NonequalAI
            (E + ADP ⇌ EADP,
             EPEP + ADP ⇌ EADPPEP)                          :: EqualAI
            EADPPEP <--> EATPPyruvate                     :: EqualAI
            (EATPPyruvate ⇌ EATP + Pyruvate,
             EPyruvate ⇌ E + Pyruvate)                          :: EqualAI
            (EATPPyruvate ⇌ EPyruvate + ATP,
             EATP ⇌ E + ATP)                                    :: EqualAI
        end
        regulatory_site(2): ATP::OnlyI
        regulatory_site(4): F16BP::OnlyA
    end
    ps_str = String.(collect(EnzymeRates.parameters(m)))
    # I-state param names use mid-name I_ token, not _T suffix.
    @test any(s -> startswith(s, "K_I_"), ps_str)
    @test !any(s -> endswith(s, "_T"), ps_str)
    # The compiled rate law must not error (intermediate state KeyErrors here).
    rng = Random.MersenneTwister(9999)
    met_names = [:PEP, :ADP, :Pyruvate, :ATP, :F16BP]
    concs_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in met_names)
    concs = NamedTuple{Tuple(met_names)}(concs_vals)
    indep = parameters(m)
    param_vals = Tuple(0.1 + 9.9 * rand(rng) for _ in indep)
    params = NamedTuple{indep}(param_vals)
    v = compile_rate_equation(m)(concs, params, 1.0)
    @test isfinite(v)
end

@testset "an allosteric Wegscheider tie leaves the grouped mechanism's model" begin
    # A binds E and EB under one shared K while B's two bindings keep their own, so
    # the RE binding square ties B's two K's (K_EAB_to_EA_B = K_EB_to_E_B) in each
    # conformation. The tied constants are dependents, so the mechanism fits the same
    # constants as the one that groups B's two bindings outright, and its law takes the
    # same values.
    tie = @enzyme_mechanism begin
        substrates: A, B
        products: P
        oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB) :: NonequalAI
            E + B ⇌ EB                       :: NonequalAI
            EA + B ⇌ EAB                 :: NonequalAI
            EAB <--> EP                  :: NonequalAI
            EP ⇌ E + P                       :: NonequalAI
        end
    end
    grouped = @enzyme_mechanism begin
        substrates: A, B
        products: P
        oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB) :: NonequalAI
            (E + B ⇌ EB, EA + B ⇌ EAB) :: NonequalAI
            EAB <--> EP                  :: NonequalAI
            EP ⇌ E + P                       :: NonequalAI
        end
    end
    @test parameters(tie) == parameters(grouped) ==
        (:K_A_EA_to_E_A, :K_A_EB_to_E_B, :K_A_EP_to_E_P, :K_I_EA_to_E_A,
         :K_I_EB_to_E_B, :K_I_EP_to_E_P, :k_A_EAB_to_EP, :k_I_EAB_to_EP, :L)
    rng = Random.MersenneTwister(4)
    tie_rate, grouped_rate = compile_rate_equation(tie), compile_rate_equation(grouped)
    for _ in 1:5
        params, Keq = _testhelper_random_independent_params(tie; rng)
        concs = (A = 0.5 + rand(rng), B = 0.5 + rand(rng), P = 0.5 + rand(rng))
        @test tie_rate(concs, params, Keq) ≈ grouped_rate(concs, params, Keq) rtol = 1e-12
    end
end

# ── The constants of tagged groups ────────────────────────────────────────────────────

@testset "Tagged groups exclude inactive-state params" begin
# A K-type uni-uni: the P binding and the chemistry :OnlyA, the S binding :EqualAI.
am = @enzyme_mechanism begin
    substrates: S
    products: P
    oligomeric_state: 2
    steps: begin
        E + S ⇌ ES    :: EqualAI
        ES <--> EP  :: OnlyA
        EP ⇌ E + P    :: OnlyA
    end
end
# The kinetic group of `am` that holds the edge written `text`.
_testhelper_group_of(text) = am.catalytic.steps[_testhelper_edge(am, text)].kinetic_group

# The :OnlyA group g of am has no inactive-state constant (`K_I_…`/`k_I_…`, rendered
# through the naming chokepoint) among the parameters. With every :OnlyA group of am
# tagged :NonequalAI those names do appear, so the first assertion can fail.
function _testhelper_test_no_inactive_constants(am, g)
    cm = am.catalytic
    i_names = EnzymeRates.Derivation._name.(
        _testhelper_step_parameters(am, first(_testhelper_group_members(cm)[g]), true),
        Ref(am))
    @test isempty(intersect(i_names, parameters(am)))
    nonequal = EnzymeRates.Mechanism(
        cm, EnzymeRates.Mechanisms.Allostery(
                replace(am.allosteric.allosteric_states, :OnlyA => :NonequalAI),
                am.allosteric.oligomeric_state, am.allosteric.regulatory_sites))
    @test !isempty(
        intersect(i_names, parameters(nonequal)))
end

@testset ":OnlyA binding group: no K_I_ param" begin
    # Must NOT be an iso-only group (iso `:OnlyA` is just a relabel — the test wants a
    # binding group whose K param disappears in the inactive state).
    _testhelper_test_no_inactive_constants(am, _testhelper_group_of("E + P ⇌ EP"))
end

@testset ":OnlyA iso group: no k_I_ params" begin
    _testhelper_test_no_inactive_constants(am, _testhelper_group_of("ES <--> EP"))
end

@testset "zero inactive-state numerator keeps the :NonequalAI binding constants" begin
    # K-type allosteric uni-uni: catalytic step is :OnlyA (so the I-state
    # numerator is zero), but binding steps are :NonequalAI.
    # When the I-state numerator is zero, the binding partition function
    # for :NonequalAI groups must still emit their I-state constants in the
    # I-state denominator so they appear in the rate-equation body and in
    # parameters(m).
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            E + P ⇌ EP      :: NonequalAI
            ES <--> EP    :: OnlyA
        end
    end
    @test iszero(EnzymeRates.Derivation._memo_law(m).inactive.Q)
    params = parameters(m)
    # K_I_ES_to_E_S and K_I_EP_to_E_P are referenced in the I-state
    # denominator of the body (the binding partition function for :NonequalAI
    # groups is built whether or not the I-state numerator is zero, since the
    # I-state denominator always appears in the rate's denominator).
    @test :K_I_ES_to_E_S in params
    @test :K_I_EP_to_E_P in params
    _testhelper_check_law_against_oracle(m)
end
end


@testset "fitted_parameter_count matches the fitted parameters on the test specs" begin
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        @test EnzymeRates.Derivation.fitted_parameter_count(m) ==
              _testhelper_fitted(m)
    end
end
