# ABOUTME: Three Haldane-valid :OnlyA-catalysis LDH mechanisms guarding sound
# ABOUTME: derivation of the dependent-parameter partition (koff EqualAI-shared,
# ABOUTME: K_I NonequalAI, kon_I SS-speed — the shapes that once UndefVar'd).

const ALLOSTERIC_UNDEFVAR_REPRODUCERS = [
    typeof(m) for m in (
        # koff_Pyruvate_ENAD family (EqualAI-shared dependent-param partition)
        @allosteric_mechanism(begin
            substrates: NADH, Pyruvate
            products: Lactate, NAD
            catalytic_multiplicity: 4
            catalytic_steps: begin
                (E + NAD ⇌ E(NAD),
                 E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate)) :: EqualAI
                (E + NADH ⇌ E(NADH),
                 E(Pyruvate) + NADH ⇌ E(NADH, Pyruvate)) :: OnlyA
                (E + Pyruvate <--> E(Pyruvate),
                 E(NADH) + Pyruvate <--> E(NADH, Pyruvate)) :: EqualAI
                (E(NAD) + Lactate ⇌ E(Lactate, NAD),
                 E(NADH) + Lactate ⇌ E(Lactate, NADH)) :: EqualAI
                E(NAD) + Pyruvate <--> E(NAD, Pyruvate) :: EqualAI
                E(NADH, Pyruvate) <--> E(Lactate, NAD) :: OnlyA
            end
        end),
        # K_NAD_ELactate family (NonequalAI dependent-param partition)
        @allosteric_mechanism(begin
            substrates: NADH, Pyruvate
            products: Lactate, NAD
            catalytic_multiplicity: 4
            catalytic_steps: begin
                (E + Lactate ⇌ E(Lactate),
                 E(NADH) + Lactate ⇌ E(Lactate, NADH)) :: NonequalAI
                E + NAD ⇌ E(NAD) :: NonequalAI
                (E + NADH ⇌ E(NADH),
                 E(Lactate) + NADH ⇌ E(Lactate, NADH),
                 E(Pyruvate) + NADH ⇌ E(NADH, Pyruvate)) :: OnlyA
                (E + Pyruvate <--> E(Pyruvate),
                 E(NAD) + Pyruvate <--> E(NAD, Pyruvate),
                 E(NADH) + Pyruvate <--> E(NADH, Pyruvate)) :: EqualAI
                E(Lactate) + NAD ⇌ E(Lactate, NAD) :: EqualAI
                E(NAD) + Lactate ⇌ E(Lactate, NAD) :: NonequalAI
                E(NADH, Pyruvate) <--> E(Lactate, NAD) :: OnlyA
                E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate) :: NonequalAI
            end
        end),
        # kon_I_NAD_E family (SS-speed dependent-param partition)
        @allosteric_mechanism(begin
            substrates: NADH, Pyruvate
            products: Lactate, NAD
            catalytic_multiplicity: 4
            catalytic_steps: begin
                E + Lactate ⇌ E(Lactate) :: EqualAI
                E + NAD <--> E(NAD) :: NonequalAI
                E + NADH ⇌ E(NADH) :: OnlyA
                (E + Pyruvate <--> E(Pyruvate),
                 E(NADH) + Pyruvate <--> E(NADH, Pyruvate)) :: EqualAI
                E(Lactate) + NAD ⇌ E(Lactate, NAD) :: EqualAI
                (E(Lactate) + NADH ⇌ E(Lactate, NADH),
                 E(Pyruvate) + NADH ⇌ E(NADH, Pyruvate)) :: OnlyA
                E(NAD) + Lactate ⇌ E(Lactate, NAD) :: EqualAI
                E(NAD) + Pyruvate <--> E(NAD, Pyruvate) :: EqualAI
                E(NADH) + Lactate ⇌ E(Lactate, NADH) :: EqualAI
                E(NADH, Pyruvate) <--> E(Lactate, NAD) :: OnlyA
                E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate) :: EqualAI
            end
        end),
    )
]
