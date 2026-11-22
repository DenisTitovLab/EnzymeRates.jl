# ABOUTME: Tests of the parameter names and units: the _name(p, m) chokepoint, the constant
# ABOUTME: columns, the regulator constants and parameter_units.

const ER = EnzymeRates

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
    # A fit starts each constant at r^rate × Π s_j^power, r the typical rate and s_j the
    # typical concentration of metabolite j; in logs, at the rate 3 and concentrations
    # 10^-j of `_testhelper_typical_table`, the centre the rule reads from the constant's
    # name (`_testhelper_typical_log_centre`). Each call returns arrays of its own.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        units = ER.Derivation.parameter_units(m)
        n, k = length(parameters(m)), length(metabolites(m))
        @test units.rate isa Vector{Int} && size(units.rate) == (n,)
        @test units.concentrations isa Matrix{Int} && size(units.concentrations) == (n, k)
        expected = [_testhelper_typical_log_centre(x, metabolites(m))
                    for x in parameters(m)]
        log_s = [log(10.0^-j) for j in 1:k]
        @test units.rate .* log(3.0) .+ units.concentrations * log_s ≈ expected atol = 1e-12
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
    names_file = joinpath(src_dir, "Derivation", "names.jl")
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
        # The renderers are in src/Derivation/names.jl and nowhere else.
        @test (_testhelper_count_chokepoint_defs(expr) > 0) == (f == names_file)
    end
end
