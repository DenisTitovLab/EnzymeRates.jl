# ABOUTME: Three Haldane-valid :OnlyA-catalysis LDH mechanisms guarding sound
# ABOUTME: derivation of the dependent-parameter partition (EqualAI-shared release
# ABOUTME: rate, NonequalAI K_I, SS binding rate k_I — the shapes that once UndefVar'd).

const ALLOSTERIC_UNDEFVAR_REPRODUCERS = [
    m for m in (
        # k_ENADPyruvate_to_ENAD_Pyruvate family (EqualAI-shared dependent-param partition)
        @enzyme_mechanism(begin
            substrates: NADH, Pyruvate
            products: Lactate, NAD
            oligomeric_state: 4
            steps: begin
                (E + NAD ⇌ ENAD,
                 EPyruvate + NAD ⇌ ENADPyruvate) :: EqualAI
                (E + NADH ⇌ ENADH,
                 EPyruvate + NADH ⇌ ENADHPyruvate) :: OnlyA
                (E + Pyruvate <--> EPyruvate,
                 ENADH + Pyruvate <--> ENADHPyruvate) :: EqualAI
                (ENAD + Lactate ⇌ ELactateNAD,
                 ENADH + Lactate ⇌ ELactateNADH) :: EqualAI
                ENAD + Pyruvate <--> ENADPyruvate :: EqualAI
                ENADHPyruvate <--> ELactateNAD :: OnlyA
            end
        end),
        # K_ELactateNAD_to_ELactate_NAD family (NonequalAI dependent-param partition)
        @enzyme_mechanism(begin
            substrates: NADH, Pyruvate
            products: Lactate, NAD
            oligomeric_state: 4
            steps: begin
                (E + Lactate ⇌ ELactate,
                 ENADH + Lactate ⇌ ELactateNADH) :: NonequalAI
                E + NAD ⇌ ENAD :: NonequalAI
                (E + NADH ⇌ ENADH,
                 ELactate + NADH ⇌ ELactateNADH,
                 EPyruvate + NADH ⇌ ENADHPyruvate) :: OnlyA
                (E + Pyruvate <--> EPyruvate,
                 ENAD + Pyruvate <--> ENADPyruvate,
                 ENADH + Pyruvate <--> ENADHPyruvate) :: EqualAI
                ELactate + NAD ⇌ ELactateNAD :: EqualAI
                ENAD + Lactate ⇌ ELactateNAD :: NonequalAI
                ENADHPyruvate <--> ELactateNAD :: OnlyA
                EPyruvate + NAD ⇌ ENADPyruvate :: NonequalAI
            end
        end),
        # k_I_E_NAD_to_ENAD family (SS-speed dependent-param partition), with the Lactate
        # binding at ENAD steady-state: at rapid equilibrium beside the steady-state NAD
        # binding it traps ENAD at zero products, which the guard rejects
        # (test_rate_equation_correctness_guard.jl).
        @enzyme_mechanism(begin
            substrates: NADH, Pyruvate
            products: Lactate, NAD
            oligomeric_state: 4
            steps: begin
                E + Lactate ⇌ ELactate :: EqualAI
                E + NAD <--> ENAD :: NonequalAI
                E + NADH ⇌ ENADH :: OnlyA
                (E + Pyruvate <--> EPyruvate,
                 ENADH + Pyruvate <--> ENADHPyruvate) :: EqualAI
                ELactate + NAD ⇌ ELactateNAD :: EqualAI
                (ELactate + NADH ⇌ ELactateNADH,
                 EPyruvate + NADH ⇌ ENADHPyruvate) :: OnlyA
                ENAD + Lactate <--> ELactateNAD :: EqualAI
                ENAD + Pyruvate <--> ENADPyruvate :: EqualAI
                ENADH + Lactate ⇌ ELactateNADH :: EqualAI
                ENADHPyruvate <--> ELactateNAD :: OnlyA
                EPyruvate + NAD ⇌ ENADPyruvate :: EqualAI
            end
        end),
    )
]
