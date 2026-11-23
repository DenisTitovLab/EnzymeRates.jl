# ABOUTME: Tests of src/Derivation/allosteric_rate_equation_correctness_guard.jl: the
# ABOUTME: readings, admission's interior and face conditions, its search and messages.

# ── Readings, admission's interior and face conditions, and their messages ────────────

using Random

@testset "a bridge's message names the tagged steps it needs" begin
    # E + A ⇌ EP_res_plusA_minusP binds A into a form the inactive conformation cuts
    # off, so it loses its affinity and stays fast: the cheapest spanning in-tree of that
    # conformation's fed forms needs only the vanishing :OnlyA binding of A at EQ, and
    # the message names it alone (the trap test's mechanism names a bridge of two). At zero
    # P and Q the ε-path converges: the rejection is conservative, as the message allows.
    bridge2 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                          :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP   :: OnlyA
            E_res_plusA_minusP + B <--> EQ                     :: NonequalAI
            EP_res_plusA_minusP + Q ⇌ EAQ                    :: EqualAI
            EQ + A ⇌ EAQ                                      :: OnlyA
        end
    end
    v = EnzymeRates.Derivation._guard(bridge2)
    @test (v.fault, v.face, v.conformation, v.edge) ==
          (EnzymeRates.Derivation._BRIDGE, 1, 2, _testhelper_edge(bridge2, "EQ + A ⇌ EAQ"))
    @test EnzymeRates.Derivation._guard_message(bridge2) ==
          "Mechanism: the :OnlyA step EQ + A ⇌ EAQ bridges forms of the " *
          "inactive conformation that its other steps keep apart, so admission cannot " *
          "show that the limit does not depend on how fast that step vanishes. Make the " *
          "kinetic group of EQ + A ⇌ EAQ steady-state."
    @test _testhelper_epath_verdict(bridge2; rng = Random.MersenneTwister(2),
                                    face = Set([:P, :Q])) === :converges
    # With the binding of A into EP_res_plusA_minusP steady-state, the guard's first
    # built reading still finds a bridge at zero P and Q; the search admits the mechanism
    # under another reading, and the oracles agree on every face.
    leading_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EP_res_plusA_minusP                       :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP   :: OnlyA
            E_res_plusA_minusP + B <--> EQ                     :: NonequalAI
            EP_res_plusA_minusP + Q ⇌ EAQ                    :: EqualAI
            EQ + A ⇌ EAQ                                      :: OnlyA
        end
    end
    w = EnzymeRates.Derivation._verdict_under(
        leading_ss, EnzymeRates.Derivation._built_reading(leading_ss))
    @test (w.fault, w.face, w.conformation) == (EnzymeRates.Derivation._BRIDGE, 2, 2)
    @test EnzymeRates.Derivation._guard(leading_ss) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(leading_ss; rng = Random.MersenneTwister(3),
                                              corners = true))
    _testhelper_check_law_against_limit(leading_ss)
end

# ── Tag sets of the admission tests ─────────────────────────────────────────

"""
Tag sets by name, each built without the guard (`@_testhelper_unguarded`): those that pass
every face under the guard's first built reading, those a face rejects under it
(`_testhelper_face_rejected`), and those whose inactive conformation excludes a
rapid-equilibrium step that reading cannot switch off, so admission's interior rejects it
(`_testhelper_pinned_tag_sets`); the search admits some of the last two kinds. The
ping-pong tag sets write their
covalent form `E_res_plusA_minusP`, which gives the graph a free form `F` would;
AD2a_iso_chemNonequal is R1, AD2e is R2 and AD2c is R3 written again, and AD6 holds one
reaction in two groups, which the constructor rejects.
"""
const _testhelper_tag_sets = Dict{String, EnzymeRates.Mechanism}()

_testhelper_tag_sets["R1"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: NonequalAI
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: NonequalAI
    end
end
_testhelper_tag_sets["A5"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: NonequalAI
        E + P ⇌ EP        :: EqualAI
        F ⇌ E               :: NonequalAI
    end
end
_testhelper_tag_sets["AD2a"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: OnlyA
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: EqualAI
    end
end
_testhelper_tag_sets["AD2d"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: EqualAI
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["AD3a"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: NonequalAI
        E + P ⇌ EP        :: EqualAI
        E + I ⇌ EI        :: OnlyA
    end
end
_testhelper_tag_sets["AD3b"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES                                       :: EqualAI
        ES <--> EP                                     :: EqualAI
        E + P ⇌ EP                                       :: EqualAI
        (E + S ⇌ ESmod, ES + S ⇌ ESSmod) :: OnlyA
        ESmod + S ⇌ ESSmod                       :: EqualAI
    end
end
_testhelper_tag_sets["AD3d"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES                    :: EqualAI
        ES <--> EP                  :: EqualAI
        E + P ⇌ EP                    :: EqualAI
        E + S <--> ESmod       :: OnlyA
        ES + S <--> ESSmod :: OnlyA
        ESmod + S ⇌ ESSmod    :: EqualAI
    end
end
_testhelper_tag_sets["AD4"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P, Q
    steps: begin
        E + S ⇌ ES            :: EqualAI
        ES <--> EPQ       :: EqualAI
        EPQ <--> E + P + Q  :: OnlyA
    end
end
_testhelper_tag_sets["AD4b"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P, Q
    steps: begin
        E + S ⇌ ES            :: EqualAI
        ES <--> EPQ       :: OnlyA
        EPQ ⇌ E + P + Q     :: EqualAI
    end
end
_testhelper_tag_sets["AD5"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E <--> F                        :: EqualAI
        (E + S ⇌ ES, F + S ⇌ FS)    :: OnlyA
        ES <--> FS                  :: EqualAI
        FS <--> FP                  :: OnlyA
        F + P ⇌ FP                    :: EqualAI
    end
end
_testhelper_tag_sets["AD7a"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A
    products: P
    modifiers: I
    steps: begin
        E + A <--> EA                 :: EqualAI
        (E + I ⇌ EI, F + I ⇌ EI)    :: OnlyA
        E + P <--> EP                 :: EqualAI
        E ⇌ F                           :: EqualAI
        EA <--> EP                  :: EqualAI
    end
end
_testhelper_tag_sets["AD7b"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A
    products: P
    modifiers: I
    steps: begin
        E + A <--> EA                 :: EqualAI
        (E + I ⇌ EI, F + I ⇌ EI)    :: EqualAI
        E + P <--> EP                 :: EqualAI
        E ⇌ F                           :: OnlyA
        EA <--> EP                  :: EqualAI
    end
end
_testhelper_tag_sets["AD7c"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A
    products: P
    modifiers: I
    steps: begin
        E + A <--> EA                 :: EqualAI
        (E + I ⇌ EI, F + I ⇌ EI)    :: EqualAI
        E + P <--> EP                 :: EqualAI
        E ⇌ F                           :: EqualAI
        EA <--> EP                  :: OnlyA
    end
end
_testhelper_tag_sets["AD8"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ EA            :: EqualAI
        EA + B <--> EQ + P  :: OnlyA
        E + Q ⇌ EQ            :: EqualAI
        EQ + B ⇌ EBQ      :: EqualAI
    end
end
_testhelper_tag_sets["AD9"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ EA            :: EqualAI
        EA + B ⇌ EAB      :: EqualAI
        EAB ⇌ EPQ       :: OnlyA
        EQ + P <--> EPQ   :: EqualAI
        E + Q <--> EQ         :: EqualAI
        EA + Q <--> EAQ   :: OnlyA
    end
end
_testhelper_tag_sets["AD10"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> EA         :: OnlyA
        E + B ⇌ EB            :: EqualAI
        EA + B ⇌ EAB      :: EqualAI
        EB + A <--> EAB   :: OnlyA
        EAB <--> EPQ    :: OnlyA
        EQ + P ⇌ EPQ      :: EqualAI
        E + Q ⇌ EQ            :: EqualAI
    end
end
# AD13: a second free form joined to E at rapid equilibrium, each label sorting before or
# after E, with the chemistry :NonequalAI or :OnlyA.
for label in (:D, :Ca, :Ea, :F, :Z), chem in (:NonequalAI, :OnlyA)
    block = quote
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: $chem
            E + P ⇌ EP        :: EqualAI
            $label ⇌ E          :: NonequalAI
        end
    end
    _testhelper_tag_sets["AD13_$(chem)_$label"] =
        EnzymeRates.DSL._construct(EnzymeRates.DSL.raw_mechanism(block))
end
_testhelper_tag_sets["AD14"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES                            :: EqualAI
        E + P ⇌ EP                            :: EqualAI
        (ES <--> EP, EIS <--> EIP)  :: OnlyA
        EI + S ⇌ EIS                      :: OnlyA
        EI + P ⇌ EIP                      :: EqualAI
        E + I ⇌ EI                            :: EqualAI
    end
end
_testhelper_tag_sets["AD14b"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES                            :: EqualAI
        E + P ⇌ EP                            :: EqualAI
        (ES <--> EP, EIS <--> EIP)  :: OnlyA
        EI + S <--> EIS                   :: OnlyA
        EI + P ⇌ EIP                      :: EqualAI
        E + I ⇌ EI                            :: EqualAI
    end
end
_testhelper_tag_sets["AD14c"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES                            :: EqualAI
        E + P ⇌ EP                            :: EqualAI
        (ES <--> EP, EIS <--> EIP)  :: OnlyA
        EI + S ⇌ EIS                      :: EqualAI
        EI + P ⇌ EIP                      :: EqualAI
        E + I ⇌ EI                            :: EqualAI
    end
end
_testhelper_tag_sets["AD15"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        (E + A ⇌ EA, EB + A ⇌ EAB)  :: OnlyA
        EA + B ⇌ EAB                  :: OnlyA
        E + B ⇌ EB                        :: EqualAI
        EAB <--> EPQ                :: OnlyA
        EQ + P ⇌ EPQ                  :: EqualAI
        E + Q ⇌ EQ                        :: EqualAI
    end
end
_testhelper_tag_sets["R2P_sealed"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> E + P     :: EqualAI
        E + P ⇌ EP        :: OnlyA
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["B_uphill"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: X
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: EqualAI
        E + P ⇌ EP        :: EqualAI
        E + X ⇌ EX        :: OnlyA
        F + X <--> EX     :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["TC1"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ EA            :: EqualAI
        EA + B ⇌ EQ + P     :: OnlyA
        EQ <--> E + Q         :: EqualAI
        EA + B <--> EAB   :: EqualAI
        EAB <--> EPQ    :: EqualAI
        EPQ ⇌ EQ + P      :: EqualAI
    end
end
_testhelper_tag_sets["pp_V_half1"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> EA                                       :: EqualAI
        EA <--> E_res_plusA_minusP + P                   :: OnlyA
        E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
        EB_res_plusA_minusP <--> E + Q                     :: EqualAI
    end
end
_testhelper_tag_sets["pp_V_half2"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> EA                                       :: EqualAI
        EA <--> E_res_plusA_minusP + P                   :: EqualAI
        E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
        EB_res_plusA_minusP <--> E + Q                     :: OnlyA
    end
end
_testhelper_tag_sets["pp_K_Abind_half1"] =
    @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> EA                                       :: OnlyA
        EA <--> E_res_plusA_minusP + P                   :: OnlyA
        E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
        EB_res_plusA_minusP <--> E + Q                     :: EqualAI
    end
end
_testhelper_tag_sets["PP1"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ EA                                          :: EqualAI
        EA <--> E_res_plusA_minusP + P                   :: OnlyA
        E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP    :: EqualAI
        EB_res_plusA_minusP <--> E + Q                     :: EqualAI
    end
end
_testhelper_tag_sets["PP2"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    modifiers: I
    steps: begin
        E + A <--> EA                                       :: EqualAI
        EA <--> E_res_plusA_minusP + P                   :: EqualAI
        E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP    :: EqualAI
        EB_res_plusA_minusP <--> E + Q                     :: OnlyA
        E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP    :: OnlyA
    end
end
# Census ter-bi d75d2415361c4dae: its faces with P zero pass by T0+ alone.
_testhelper_tag_sets["terbi_d75d2415361c4dae"] =
    @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B, D
    products: P, Q
    steps: begin
        (E + A ⇌ EA, EBP + A ⇌ EABP, EP + A ⇌ EAP)       :: OnlyA
        (E + P ⇌ EP, EA + P ⇌ EAP, EAB + P ⇌ EABP,
         EQ + P ⇌ EPQ)                                               :: EqualAI
        (E + Q <--> EQ, EP + Q <--> EPQ)                           :: OnlyA
        (EA + B <--> EAB, EAP + B <--> EABP,
         EP + B <--> EBP)                                            :: EqualAI
        (EAB + D <--> EPQ, EQ + D <--> EDQ)                  :: OnlyA
    end
end
# An :OnlyA group mixing B's fused binding into EPQ with its plain binding at EQ,
# which no untagged step joins to EBQ: its limit exists.
_testhelper_tag_sets["mixed_group_limit"] =
    @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ EA                                    :: EqualAI
        (EA + B <--> EPQ, EQ + B <--> EBQ)  :: OnlyA
        EQ + P ⇌ EPQ                              :: EqualAI
        E + Q ⇌ EQ                                    :: EqualAI
    end
end
_testhelper_tag_sets["R2"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: OnlyA
        ES <--> EP      :: EqualAI
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["BS1"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> EP      :: EqualAI
        EP ⇌ E + P        :: EqualAI
        E <--> F            :: OnlyA
        E + I ⇌ FI        :: OnlyA
        F + I ⇌ FI        :: EqualAI
        F + S ⇌ FS        :: EqualAI
        FS <--> FP      :: EqualAI
        FP ⇌ F + P        :: EqualAI
    end
end
_testhelper_tag_sets["BS2"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: OnlyA
        ES <--> EP      :: EqualAI
        F + P ⇌ EP        :: EqualAI
        F ⇌ E               :: OnlyA
    end
end
_testhelper_tag_sets["BS3"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S <--> ES     :: OnlyA
        ES <--> EP      :: EqualAI
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["BS4"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES            :: EqualAI
        ES <--> E + P         :: EqualAI
        E + I ⇌ EI            :: OnlyA
        EI + S <--> EIS   :: EqualAI
        EIS <--> FI + P   :: EqualAI
        FI ⇌ F + I            :: EqualAI
        F <--> E                :: OnlyA
    end
end
_testhelper_tag_sets["ADbridge"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S <--> ES                     :: EqualAI
        ES <--> E + P                     :: EqualAI
        (E + I ⇌ EI, ES + I ⇌ EIS)  :: OnlyA
        EI + S ⇌ EIS                  :: EqualAI
    end
end
_testhelper_tag_sets["ADbridge_noturnover"] =
    @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S <--> ES                     :: EqualAI
        ES <--> E + P                     :: OnlyA
        (E + I ⇌ EI, ES + I ⇌ EIS)  :: OnlyA
        EI + S ⇌ EIS                  :: EqualAI
    end
end
_testhelper_tag_sets["bridge_drain"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ ES        :: EqualAI
        ES <--> F + P     :: EqualAI
        F <--> E            :: OnlyA
        E + I ⇌ EI        :: OnlyA
        F + I ⇌ FI        :: OnlyA
        FI ⇌ EI         :: EqualAI
    end
end
_testhelper_tag_sets["bridge_deadend"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I, J
    steps: begin
        E + S <--> ES         :: EqualAI
        ES <--> E + P         :: EqualAI
        E + J <--> EJ         :: EqualAI
        E + I ⇌ EI            :: OnlyA
        EJ + I ⇌ EIJ      :: OnlyA
        EI + J ⇌ EIJ      :: EqualAI
    end
end
_testhelper_tag_sets["R2funnel"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ ES        :: OnlyA
        ES <--> EP      :: EqualAI
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: OnlyA
        F + S ⇌ FS        :: EqualAI
        FS <--> FP      :: EqualAI
        F + P ⇌ FP        :: EqualAI
        FP ⇌ EP         :: EqualAI
    end
end
_testhelper_tag_sets["R2funnel_SS"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S <--> ES     :: OnlyA
        ES <--> EP      :: EqualAI
        F + P ⇌ EP        :: EqualAI
        F <--> E            :: OnlyA
        F + S ⇌ FS        :: EqualAI
        FS <--> FP      :: EqualAI
        F + P ⇌ FP        :: EqualAI
        FP ⇌ EP         :: EqualAI
    end
end
_testhelper_tag_sets["S1"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: A, B
    products: P
    steps: begin
        E + A <--> EA                 :: OnlyA
        EA + B <--> EAB           :: EqualAI
        EAB <--> EP               :: OnlyA
        EP <--> E + P                 :: EqualAI
        E + A ⇌ EAmod          :: OnlyA
        EAmod + B ⇌ EAmodB    :: EqualAI
        EB + A ⇌ EAmodB    :: EqualAI
        F + B <--> EB                 :: EqualAI
    end
end
# Fuzz mechanism fz101_10232 (seed 101).
_testhelper_tag_sets["fz101_10232"] = @_testhelper_unguarded @enzyme_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S <--> ES                             :: OnlyA
        ES <--> EP                              :: OnlyA
        EP <--> F + P                             :: EqualAI
        F <--> E                                    :: OnlyA
        (ES + I <--> EIS, F + I <--> FI)    :: OnlyA
        EP + S ⇌ ES + P                         :: EqualAI
    end
end

@testset "how an excluded step is switched off" begin
    D = EnzymeRates.Derivation
    rng = Random.MersenneTwister(2611)
    # The S binding takes up one metabolite into ES, which the inactive conformation's
    # P binding no longer joins to E: it loses its affinity. The chemistry takes up
    # nothing, so it has no affinity to lose; at steady state it is stopped.
    k_type = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: OnlyA
            E + P ⇌ EP        :: EqualAI
        end
    end
    off = D._built_reading(k_type)[:, 2]
    @test off[_testhelper_edge(k_type, "E + S ⇌ ES")] === D._LOSES_FORWARD
    @test off[_testhelper_edge(k_type, "ES <--> EP")] === D._STOPPED
    @test off[_testhelper_edge(k_type, "E + P ⇌ EP")] === D._RUNS
    @test all(==(D._RUNS), D._built_reading(k_type)[:, 1])
    @test D._guard(k_type).fault === D._NO_FAULT
    @test all(_testhelper_admission_agreement(k_type; rng))
    # With the chemistry untagged, ES stays, made from P: thermodynamics pins the S
    # binding's equilibrium constant, and a rapid-equilibrium step cannot turn slow.
    pinned = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
    end
    v = D._guard(pinned)
    @test (v.fault, v.conformation, v.edge) ==
          (D._FAST_STEP_PINNED, 2, _testhelper_edge(pinned, "E + S ⇌ ES"))
    @test _testhelper_epath_verdict(pinned; rng) === :no_path
    # The same binding at steady state is stopped: S neither binds nor leaves in the
    # inactive conformation, and ES is made from P there.
    shut_door = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES     :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
    end
    @test D._STOPPED in D._built_reading(shut_door)[:, 2]
    @test D._guard(shut_door).fault === D._NO_FAULT
    @test all(_testhelper_admission_agreement(shut_door; rng))
    # A rapid-equilibrium step that binds B and reacts at once, into EPQ, which the
    # P binding still makes: it cannot be switched off.
    fused_re = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA         :: EqualAI
            E + Q ⇌ EQ            :: EqualAI
            EA + B ⇌ EPQ      :: OnlyA
            EQ + P <--> EPQ   :: EqualAI
        end
    end
    @test D._guard(fused_re).fault === D._FAST_STEP_PINNED
    # At steady state the same step is stopped.
    fused_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA         :: EqualAI
            E + Q ⇌ EQ            :: EqualAI
            EA + B <--> EPQ   :: OnlyA
            EQ + P <--> EPQ   :: EqualAI
        end
    end
    @test D._guard(fused_ss).fault === D._NO_FAULT
    @test all(_testhelper_admission_agreement(fused_ss; rng))
    # The package's reading agrees with the test's own statement of the rule.
    for x in (k_type, pinned, shut_door, fused_re, fused_ss)
        @test _testhelper_reading(D._built_reading(x)) ==
              _testhelper_reading(x, :all_lose_affinity)
    end
end

@testset "a fast isomerization excluded between two groups loses a direction" begin
    # The isomerization takes up nothing, so no built reading switches it off. But the
    # steps the inactive conformation keeps do not join its two sides: EP lies in a group
    # of its own, and the isomerization loses its rate into EP, staying fast. The search
    # admits the mechanism, and the oracles agree on every face.
    D = EnzymeRates.Derivation
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES ⇌ EP         :: OnlyA
            E + P <--> EP     :: OnlyA
        end
    end
    iso = _testhelper_edge(am, "ES ⇌ EP")
    @test D._admission_interior(am, D._built_reading(am)) ==
          D._Verdict(D._FAST_STEP_PINNED, D._LIMIT_UNTESTED, 0, 2, iso)
    @test D._guard(am) == D._WELL_POSED
    @test D._guard_with_reading(am).reading[iso, 2] === D._LOSES_FORWARD
    @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(2615),
                                              corners = true))
    _testhelper_check_law_against_limit(am)
end

@testset "a group that loses one affinity names its step that cannot be switched off" begin
    # EA is a dead end, cut off in the inactive conformation, so the A binding at E
    # loses its affinity; its groupmate binds A into EPQ, which the P binding still
    # makes, at rapid equilibrium: the verdict names the groupmate.
    D = EnzymeRates.Derivation
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EPQ)  :: OnlyA
            E + B ⇌ EB                        :: EqualAI
            EQ + P <--> EPQ               :: EqualAI
            E + Q ⇌ EQ                        :: EqualAI
        end
    end
    lost, pinned = _testhelper_edge(am, "E + A ⇌ EA"), _testhelper_edge(am, "EB + A ⇌ EPQ")
    off = D._built_reading(am)[:, 2]
    @test (off[lost], off[pinned]) == (D._LOSES_FORWARD, D._CANNOT_SWITCH_OFF)
    v = D._guard(am)
    @test (v.fault, v.conformation, v.edge) == (D._FAST_STEP_PINNED, 2, pinned)
end

@testset "an uptake from one cut-off form into another loses its affinity" begin
    # Every step :OnlyA: the inactive conformation keeps only E. The B binding takes up
    # one metabolite from EA, cut off, into EAB, cut off too, so it loses its
    # affinity, and free E reaches both along the bindings.
    D = EnzymeRates.Derivation
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA            :: OnlyA
            EA + B ⇌ EAB      :: OnlyA
            EAB <--> EPQ    :: OnlyA
            EPQ ⇌ EQ + P      :: OnlyA
            EQ ⇌ E + Q            :: OnlyA
        end
    end
    off = D._built_reading(am)[:, 2]
    @test off[_testhelper_edge(am, "E + A ⇌ EA")] === D._LOSES_FORWARD
    @test off[_testhelper_edge(am, "EA + B ⇌ EAB")] === D._LOSES_FORWARD
    @test D._guard(am).fault === D._NO_FAULT
    @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(2612)))
end

@testset "how an excluded step is switched off depends on neither subunits nor sites" begin
    # The pinned S binding and its steady-state variant of "how an excluded step is
    # switched off", with two catalytic subunits and a regulatory site.
    D = EnzymeRates.Derivation
    pinned = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    shut_door = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S <--> ES     :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    v = D._guard(pinned)
    @test (v.fault, v.conformation, v.edge) ==
          (D._FAST_STEP_PINNED, 2, _testhelper_edge(pinned, "E + S ⇌ ES"))
    @test D._guard(shut_door).fault === D._NO_FAULT
end

@testset "a steady-state uptake into a cut-off form may be stopped" begin
    # The guard first reads every excluded steady-state step that takes up one metabolite
    # into a cut-off form as losing its affinity; when that leaves an admission fault, it
    # checks again with such a step losing its affinity only when it leads away from free
    # E, then with every such step stopped. Four ter-ter mechanisms of the gate's
    # populations, built unguarded since the first reading rejects each: one that the
    # second reading admits, one that only the third admits, and two whose bridges are
    # rapid-equilibrium steps, which no reading stops.
    D = EnzymeRates.Derivation
    rng = Random.MersenneTwister(2613)
    second_try_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + D ⇌ ED, EBQ + D ⇌ EBDQ, EPQ + D ⇌ EDPQ,
             EQ + D ⇌ EDQ) :: EqualAI
            (E + Q ⇌ EQ, EBD + Q ⇌ EBDQ, ED + Q ⇌ EDQ,
             EDP + Q ⇌ EDPQ) :: OnlyA
            EBD + A <--> EPQR :: OnlyA
            (ED + B ⇌ EBD, EDQ + B ⇌ EBDQ, EQ + B ⇌ EBQ) :: OnlyA
            (ED + P <--> EDP, EDQ + P <--> EDPQ,
             EQ + P <--> EPQ, EQR + P <--> EPQR) :: EqualAI
            (EPQ + R ⇌ EPQR, EQ + R ⇌ EQR) :: EqualAI
        end
    end
    third_try = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + D ⇌ ED, EAQ + D ⇌ EADQ, EPQ + D ⇌ EDPQ,
             EQ + D ⇌ EDQ) :: OnlyA
            (E + Q ⇌ EQ, EAD + Q ⇌ EADQ, ED + Q ⇌ EDQ,
             EDP + Q ⇌ EDPQ) :: EqualAI
            (EAD + B ⇌ EPQR, ED + B ⇌ EBD) :: OnlyA
            (EBD + A <--> EPQR, ED + A <--> EAD,
             EDQ + A <--> EADQ, EQ + A <--> EAQ) :: EqualAI
            (ED + P ⇌ EDP, EDQ + P ⇌ EDPQ, EQ + P ⇌ EPQ) :: EqualAI
            EPQ + R <--> EPQR :: OnlyA
        end
    end
    # second_try_a, inactive conformation: the chemistry step joins two cut-off forms.
    k = _testhelper_edge(second_try_a, "EBD + A <--> EPQR")
    losing, stopped = D._built_reading(second_try_a)[:, 2],
                      D._built_reading(second_try_a, D._ALL_STOPPED)[:, 2]
    @test losing[k] === D._LOSES_FORWARD
    @test stopped[k] === D._STOPPED
    @test [losing[j] for j in eachindex(losing) if j != k] ==
          [stopped[j] for j in eachindex(stopped) if j != k]
    @test _testhelper_reading(second_try_a, :all_lose_affinity)[k, 2] === :loses_forward
    @test _testhelper_reading(second_try_a, :all_stopped)[k, 2] === :stopped
    # The two conformations are read alike.
    for r in (D._ALL_LOSE_AFFINITY, D._ALL_STOPPED)
        @test D._built_reading(_testhelper_mirror(second_try_a), r)[:, 1] ==
              D._built_reading(second_try_a, r)[:, 2]
    end
    # The first reading alone finds a bridge at zero substrates in the inactive
    # conformation.
    v = D._verdict_under(second_try_a, D._built_reading(second_try_a))
    @test (_testhelper_faces(second_try_a)[v.face], v.conformation, v.edge) ==
          (Set([:A, :B, :D]), 2, _testhelper_edge(second_try_a, "EDP + Q ⇌ EDPQ"))
    for am in (second_try_a, third_try)
        # The first try alone.
        @test D._verdict_under(am, D._built_reading(am)).fault === D._BRIDGE
        @test D._verdict_under(am, D._built_reading(am, D._ALL_STOPPED)) == D._WELL_POSED
        @test D._guard(am) == D._WELL_POSED
        @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(1)))
        # The compiled law, derived through the guard, against the oracle's own deletion
        # limit.
        _testhelper_check_law_against_limit(am)
    end
    # The outward reading admits second_try_a but not third_try, so each later try is
    # needed.
    @test D._verdict_under(second_try_a,
                           D._built_reading(second_try_a, D._OUTWARD_LOSE_AFFINITY)) ==
          D._WELL_POSED
    @test D._verdict_under(
        third_try, D._built_reading(third_try, D._OUTWARD_LOSE_AFFINITY)).fault ===
          D._BRIDGE
    # bridge_fast_a: its bridges, the rapid-equilibrium :OnlyA bindings such as the P
    # binding EDQ + P ⇌ EDPQ, lose their affinity under every reading, and stopping its
    # steady-state uptakes leaves forms that free E does not reach. The guard names the
    # first bridge in its search order, the B binding E_res_plusA_plusD_minusR + B ⇌ EPQ.
    # At zero P the ε-path oracle finds that the limit depends on the direction: the
    # rejection is real.
    bridge_fast_a =@_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + D ⇌ ED, EAP + D ⇌ EADP, EP + D ⇌ EDP,
             EPQ + D ⇌ EDPQ) :: EqualAI
            (E + P ⇌ EP, ED + P ⇌ EDP, EDQ + P ⇌ EDPQ,
             ER_res_plusA_plusD_minusR + P ⇌ EADP) :: OnlyA
            E_res_plusA_plusD_minusR + B ⇌ EPQ :: OnlyA
            E_res_plusA_plusD_minusR + R ⇌ ER_res_plusA_plusD_minusR :: EqualAI
            (EAP + Q <--> EAPQ, ED + Q <--> EDQ,
             EDP + Q <--> EDPQ, EP + Q <--> EPQ,
             ER_res_plusA_plusD_minusR + Q <--> EADQ) :: OnlyA
            (ED + A <--> ER_res_plusA_plusD_minusR, EDP + A <--> EADP,
             EDQ + A <--> EADQ, EP + A <--> EAP,
             EPQ + A <--> EAPQ) :: EqualAI
        end
    end
    v = D._guard(bridge_fast_a)
    @test v == D._verdict_under(bridge_fast_a, D._built_reading(bridge_fast_a))
    @test (v.fault, v.conformation, v.edge) ==
          (D._BRIDGE, 2,
           _testhelper_edge(bridge_fast_a, "E_res_plusA_plusD_minusR + B ⇌ EPQ"))
    @test D._verdict_under(
        bridge_fast_a, D._built_reading(bridge_fast_a, D._ALL_STOPPED)).fault === D._ISLAND
    @test D._verdict_under(
        bridge_fast_a, D._built_reading(bridge_fast_a, D._OUTWARD_LOSE_AFFINITY)).fault ===
          D._BRIDGE
    @test _testhelper_epath_verdict(bridge_fast_a; rng, face = Set([:P])) ===
          :direction_dependent
    # bridge_fast_b: every excluded uptake into a cut-off form is a rapid-equilibrium step,
    # so there is no steady-state one to read again.
    bridge_fast_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A ⇌ EA, EBQ + A ⇌ EABQ,
             EP_res_plusA_plusB_minusP + A ⇌ EAP_res_plusA_plusB_minusP,
             EQ + A ⇌ EAQ) :: OnlyA
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ,
             EP_res_plusA_plusB_minusP + Q ⇌ EABQ) :: OnlyA
            (E_res_plusA_plusB_minusP + P ⇌ EP_res_plusA_plusB_minusP,
             EA + P ⇌ EAP,
             EP_res_plusA_plusB_minusP + P ⇌ EABP) :: EqualAI
            E_res_plusA_plusB_minusP + D ⇌ EQR :: EqualAI
            (EA + B <--> EP_res_plusA_plusB_minusP, EAP + B <--> EABP,
             EAQ + B <--> EABQ,
             EP_res_plusA_plusB_minusP + B <--> EBP_res_plusA_plusB_minusP,
             EQ + B <--> EBQ) :: EqualAI
            (EQ + R <--> EQR,
             EQR + R <--> EDR_res_plusA_plusB_minusP) :: EqualAI
        end
    end
    @test D._guard(bridge_fast_b).fault === D._BRIDGE
    for c in 1:2
        @test D._built_reading(bridge_fast_b)[:, c] ==
              D._built_reading(bridge_fast_b, D._ALL_STOPPED)[:, c]
    end
    # The first try decides when it passes: with its S binding stopped, the K-type
    # uni-uni mechanism would leave ES out of free E's reach.
    k_type = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES     :: OnlyA
            ES <--> EP      :: OnlyA
            EP ⇌ E + P        :: EqualAI
        end
    end
    @test D._guard(k_type) == D._WELL_POSED
    @test D._admission_interior(k_type, D._built_reading(k_type, D._ALL_STOPPED)).fault ===
          D._ISLAND
end

@testset "a mechanism may need one uptake to lose its affinity and another to stop" begin
    # The guard's second try (`_OUTWARD_LOSE_AFFINITY`): an excluded steady-state step that
    # takes up one metabolite into a cut-off form loses its affinity only when it leads
    # away from free E, its bound side's part farther from free E's than its unbound
    # side's; otherwise it is stopped. A part's distance counts the excluded uptakes into
    # cut-off forms on the way to it from free E's part. The fixtures come from the gate's
    # populations.
    D = EnzymeRates.Derivation
    # both_readings, inactive conformation: the steps it keeps join EA to EAP and
    # EP to EPQ, each part one uptake from free E's. The B binding, from one part to
    # the other, must stop, or the steps that lose their affinity close a cycle with the
    # A binding at EP; the P binding at E must then lose its affinity, or no such step
    # leads to EP and EPQ.
    both_readings = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP)  :: OnlyA
            E + P <--> EP                     :: OnlyA
            EA + P <--> EAP               :: EqualAI
            EA + B <--> EPQ               :: OnlyA
            EP + Q ⇌ EPQ                  :: EqualAI
        end
    end
    p = _testhelper_edge(both_readings, "E + P <--> EP")      # leads away: distance 0 to 1
    b = _testhelper_edge(both_readings, "EA + B <--> EPQ")    # distance 1 to 1
    a = _testhelper_edge(both_readings, "EP + A ⇌ EAP")       # rapid equilibrium, 1 to 1
    L, S = D._LOSES_FORWARD, D._STOPPED
    for (reading, name, at_p, at_b) in (
            (D._ALL_LOSE_AFFINITY, :all_lose_affinity, L, L),
            (D._OUTWARD_LOSE_AFFINITY, :outward_lose_affinity, L, S),
            (D._ALL_STOPPED, :all_stopped, S, S))
        off = D._built_reading(both_readings, reading)[:, 2]
        @test off[p] === at_p && off[b] === at_b
        @test off[a] === D._LOSES_FORWARD                     # a fast step is read one way
        @test _testhelper_reading(D._built_reading(both_readings, reading))[:, 2] ==
              _testhelper_reading(both_readings, name)[:, 2]
    end
    # The same need with the substrates' roles swapped and the P binding at EB :OnlyI.
    both_readings_only_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + B ⇌ EB, EP + B ⇌ EBP)  :: OnlyA
            E + P <--> EP                     :: OnlyA
            EB + P <--> EBP               :: OnlyI
            EB + A <--> EPQ               :: OnlyA
            EP + Q ⇌ EPQ                  :: EqualAI
        end
    end
    for am in (both_readings, both_readings_only_i)
        @test D._verdict_under(am, D._built_reading(am)).fault === D._ARC_CYCLE
        @test D._verdict_under(am, D._built_reading(am, D._ALL_STOPPED)).fault === D._ISLAND
        @test D._verdict_under(am, D._built_reading(am, D._OUTWARD_LOSE_AFFINITY)) ==
              D._WELL_POSED
        @test D._guard(am) == D._WELL_POSED
        @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(1)))
        # The compiled law, derived through the guard, against the oracle's own deletion
        # limit.
        _testhelper_check_law_against_limit(am)
    end
    # outward_pumped: the first try admits it, and the outward reading alone would lose
    # it.
    outward_pumped = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + A ⇌ EA                                    :: OnlyA
            E + R <--> ER                                 :: OnlyA
            (EA + B <--> EAB, ER + B <--> EBR)  :: OnlyA
            (EA + P ⇌ EAP, ER + P ⇌ EPR)        :: EqualAI
            EAB + D <--> EPR + Q                    :: EqualAI
        end
    end
    @test D._verdict_under(outward_pumped, D._built_reading(outward_pumped)) ==
          D._WELL_POSED
    outward = D._built_reading(outward_pumped, D._OUTWARD_LOSE_AFFINITY)
    @test D._verdict_under(outward_pumped, outward).fault === D._PUMPED
    @test D._guard(outward_pumped) == D._WELL_POSED
    @test all(_testhelper_admission_agreement(outward_pumped;
                                              rng = Random.MersenneTwister(1)))
    # In the inactive conformation no uptake leads from free E to E_res_plusA_minusP,
    # which only the chemistry and the P release reach: its B binding leaves a part no
    # uptake reaches, so the outward reading stops it.
    unreached = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P      :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP   :: OnlyA
            EB_res_plusA_minusP <--> EQ                        :: OnlyA
            EQ ⇌ E + Q                                            :: OnlyA
        end
    end
    k = _testhelper_edge(unreached, "E_res_plusA_minusP + B <--> EB_res_plusA_minusP")
    @test D._built_reading(unreached)[k, 2] === D._LOSES_FORWARD
    @test D._built_reading(unreached, D._OUTWARD_LOSE_AFFINITY)[k, 2] === D._STOPPED
    # The package's reading agrees with the test's own statement of the rule, and the two
    # conformations are read alike.
    for x in (both_readings, both_readings_only_i, outward_pumped, unreached),
        (name, r) in pairs(_testhelper_guard_readings)
        @test _testhelper_reading(D._built_reading(x, r)) == _testhelper_reading(x, name)
        @test D._built_reading(_testhelper_mirror(x), r)[:, 1] ==
              D._built_reading(x, r)[:, 2]
    end
end

@testset "a reading may lose a step in either direction" begin
    # A ping-pong with both chemistry steps :OnlyA. In the
    # inactive conformation the covalent forms F, FB and FP form one cut-off group
    # (F is E_res_plusA_minusP). EA <--> FP loses its stored direction, into that
    # group, and FB ⇌ EQ its reverse, also into that group: both point uphill from
    # free E's group, and the cut-off group empties through both halves. No built reading
    # loses a step in reverse.
    D = EnzymeRates.Derivation
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
        end
    end
    chem_a = _testhelper_edge(am, "EA <--> EP_res_plusA_minusP")
    chem_b = _testhelper_edge(am, "EB_res_plusA_minusP ⇌ EQ")
    r = fill(D._RUNS, length(am.catalytic.steps), 2)
    r[chem_a, 2], r[chem_b, 2] = D._LOSES_FORWARD, D._LOSES_REVERSE
    @test D._admission_interior(am, r) == D._WELL_POSED
    eps = D._eps_network(am, D._re_flags(am.catalytic), r)
    fwd, rev = (only(t for t in eachindex(eps.from)
                     if eps.edge[t] == chem_b && eps.conformation[t] == 2 &&
                        eps.forward[t] == forward)
                for forward in (true, false))
    # FB ⇌ EQ in copy 2: the reverse (EQ → FB) vanishes, the forward does not,
    # both fast.
    @test [eps.vanishing[t] for t in (fwd, rev)] == [false, true]
    @test eps.fast[fwd] && eps.fast[rev]
    # Exchanged, both steps point into free E's group, below which nothing lies: no
    # energies exist.
    flipped = copy(r)
    flipped[chem_a, 2], flipped[chem_b, 2] = D._LOSES_REVERSE, D._LOSES_FORWARD
    # The verdict names the first of them in edge order.
    @test D._admission_interior(am, flipped) ==
          D._Verdict(D._ARC_CYCLE, D._LIMIT_UNTESTED, 0, 2, chem_a)
    # A step keeps its kind, so a reading that stops the rapid-equilibrium FB ⇌ EQ is
    # refused as one that cannot switch it off.
    stopped_fast = copy(r)
    stopped_fast[chem_b, 2] = D._STOPPED
    @test D._admission_interior(am, stopped_fast) ==
          D._Verdict(D._FAST_STEP_PINNED, D._LIMIT_UNTESTED, 0, 2, chem_b)
    stopped_names = _testhelper_reading(stopped_fast)
    @test !_testhelper_reading_is_consistent(am, stopped_names)
    @test _testhelper_exclusion(am, 2; reading = stopped_names).failure === :fast
    # A step the conformation keeps joins its two sides into one group, so a direction it
    # loses points uphill inside that group and closes a cycle alone.
    b_binding = _testhelper_edge(am, "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP")
    inside = copy(r)
    inside[b_binding, 2] = D._LOSES_FORWARD
    @test D._admission_interior(am, inside) ==
          D._Verdict(D._ARC_CYCLE, D._LIMIT_UNTESTED, 0, 2, b_binding)
    # The tests' own names for the reading's entries.
    expected = fill(:runs, length(am.catalytic.steps), 2)
    expected[chem_a, 2], expected[chem_b, 2] = :loses_forward, :loses_reverse
    @test _testhelper_reading(r) == expected
    # The ε-path oracle, fed the same reading, finds the limit on every face; under the
    # guard's first built reading FB ⇌ EQ takes up nothing, so it cannot be switched
    # off, and there is no path.
    rng = Random.MersenneTwister(2614)
    for face in _testhelper_faces(am)
        @test _testhelper_epath_verdict(am; rng, face, corners = true,
                                        reading = _testhelper_reading(r)) === :converges
    end
    @test _testhelper_epath_verdict(
        am; rng, reading = _testhelper_reading(am, :all_lose_affinity)) === :no_path
end

@testset "no form fills faster than it empties" begin
    D = EnzymeRates.Derivation
    NO, FI, U = D._NO_FAULT, D._FILLED, D._LIMIT_UNTESTED
    rng = Random.MersenneTwister(2650)
    # `filled`: a ping-pong whose B binding at the covalent form F (E_res_plusA_minusP)
    # and first chemistry step are :OnlyA. In the inactive conformation F and FP form a
    # cut-off group: EA <--> FP loses its stored direction, into it, and F + B ⇌ FB
    # its reverse, the release of B, also into it. With B zero the uptake of B dies, and
    # the release, fast though it vanishes, runs one way from the law's form FB into F,
    # whose group empties only through the slow chemistry FP → EA: F fills faster
    # than it empties.
    filled = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: EqualAI
        end
    end
    chem_a = _testhelper_edge(filled, "EA <--> EP_res_plusA_minusP")
    b_binding = _testhelper_edge(filled, "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP")
    r = fill(D._RUNS, length(filled.catalytic.steps), 2)
    r[chem_a, 2], r[b_binding, 2] = D._LOSES_FORWARD, D._LOSES_REVERSE
    faces = _testhelper_faces(filled)
    # The first face with B zero zeroes A as well; B alone is refused the same way.
    v = D._verdict_under(filled, r)
    @test v == D._Verdict(FI, U, 3, 2, b_binding)
    @test faces[3] == Set([:A, :B])
    @test faces[5] == Set([:B])
    @test _testhelper_admission_faults(filled; reading = _testhelper_reading(r)) ==
          [NO, NO, FI, NO, FI, NO, NO]
    # The message names the face, the step, the form it fills and the group's slow way
    # out, and offers the fix the guard admits: F + B <--> FB steady-state.
    @test D._steady_state_passes(filled, b_binding)
    @test D._admission_message(filled, v, r) ==
          "Mechanism: at zero A and B, in the inactive conformation the :OnlyA " *
          "rapid-equilibrium step E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP " *
          "runs one way into E_res_plusA_minusP, at a rate that vanishes yet grows " *
          "with the rapid-equilibrium speed, and the forms EP_res_plusA_minusP, " *
          "E_res_plusA_minusP leave only through the slower step " *
          "EA <--> EP_res_plusA_minusP, so E_res_plusA_minusP fills faster than " *
          "it empties and its weight depends on how fast the tagged step vanishes. Make " *
          "the kinetic group of E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP " *
          "steady-state."
    filled_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP   :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: EqualAI
        end
    end
    @test D._guard(filled_ss) == D._WELL_POSED
    @test all(_testhelper_admission_agreement(filled_ss; rng = Random.MersenneTwister(2653),
                                              corners = true))
    _testhelper_check_law_against_limit(filled_ss)
    # The harsh ε-path oracle under the same reading: the limit depends on the direction
    # on both faces with B zero, and exists on the interior.
    for i in (3, 5)
        @test _testhelper_epath_verdict(filled; rng, face = faces[i], corners = true,
                                        reading = _testhelper_reading(r)) ===
              :direction_dependent
    end
    @test _testhelper_epath_verdict(filled; rng, corners = true,
                                    reading = _testhelper_reading(r)) === :converges

    # The six-step ping-pong with both chemistry steps :OnlyA, under the reading that
    # admits it: neither excluded step takes up a metabolite, so no face kills the
    # direction it keeps, and nothing fills one way.
    pingpong = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
        end
    end
    p = fill(D._RUNS, length(pingpong.catalytic.steps), 2)
    p[_testhelper_edge(pingpong, "EA <--> EP_res_plusA_minusP"), 2] = D._LOSES_FORWARD
    p[_testhelper_edge(pingpong, "EB_res_plusA_minusP ⇌ EQ"), 2] = D._LOSES_REVERSE
    @test D._verdict_under(pingpong, p) == D._WELL_POSED

    # `slow_source`: a ping-pong whose B binding at F gives EQ at once, every step
    # :OnlyA but F + P ⇌ FP. In the inactive conformation free E stands alone and the
    # steps that leave it are steady-state uptakes, which lose their uptake. With B zero
    # F + B ⇌ EQ runs one way, fast though it vanishes, from EQ into F; but enzyme
    # reaches EQ only through steady-state steps, so the flow into F vanishes, and the
    # faces with B zero pass.
    slow_source = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EP_res_plusA_minusP, EQ + A <--> EAQ)    :: OnlyA
            (E + Q <--> EQ, EP_res_plusA_minusP + Q <--> EAQ)    :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP            :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                              :: OnlyA
        end
    end
    b_step = _testhelper_edge(slow_source, "E_res_plusA_minusP + B ⇌ EQ")
    s = fill(D._RUNS, length(slow_source.catalytic.steps), 2)
    for text in ("E + A <--> EP_res_plusA_minusP", "E + Q <--> EQ",
                 "EP_res_plusA_minusP + Q <--> EAQ", "EQ + A <--> EAQ")
        s[_testhelper_edge(slow_source, text), 2] = D._LOSES_FORWARD
    end
    s[b_step, 2] = D._LOSES_REVERSE
    @test D._verdict_under(slow_source, s) == D._WELL_POSED
    eps = D._eps_network(slow_source, D._re_flags(slow_source.catalytic), s)
    law_nodes = D._law_network(slow_source).present
    sets = D._face_sets(eps, law_nodes, D._faces(slow_source.catalytic)[3])
    fast_from_law = D._reached(D._successors(eps, sets.alive .& eps.fast),
                               findall(law_nodes))
    uptake, release = (only(t for t in eachindex(eps.from)
                            if eps.edge[t] == b_step && eps.conformation[t] == 2 &&
                               eps.forward[t] == forward)
                       for forward in (true, false))
    # At zero A and B, the release EQ → F + B is fast, vanishing and alive while the
    # uptake is dead; its source is fed, but not along fast transitions from the law's
    # forms.
    @test eps.fast[release] && eps.vanishing[release] && !eps.vanishing[uptake]
    @test sets.alive[release] && !sets.alive[uptake]
    @test sets.fed[eps.from[release]] && !fast_from_law[eps.from[release]]
    for face in _testhelper_faces(slow_source)
        @test _testhelper_epath_verdict(slow_source; rng, face, corners = true,
                                        reading = _testhelper_reading(s)) === :converges
    end

    # No built reading raises the fault. Under each of them `filled`'s excluded steps are
    # stopped or cannot be switched off, so none keeps a direction; `slow_source`'s
    # F + B ⇌ EQ loses its affinity, and the direction it keeps is the release of B,
    # which takes up nothing, so no face kills it.
    for am in (filled, slow_source)
        law, built_faces = D._law_network(am), D._faces(am.catalytic)
        for b in instances(D._SlowUptakeReading)
            eps = D._eps_network(am, D._re_flags(am.catalytic), D._built_reading(am, b))
            @test all(i -> D._admission_face(am, eps, law, built_faces, i).fault !== FI,
                      eachindex(built_faces))
        end
    end
    @test D._built_reading(slow_source)[b_step, 2] === D._LOSES_FORWARD
end

@testset "cut-off forms that cannot matter are set aside" begin
    D = EnzymeRates.Derivation
    NO, BR, U = D._NO_FAULT, D._BRIDGE, D._LIMIT_UNTESTED
    rng = Random.MersenneTwister(2651)
    # Whether face `i` of `am` under reading `r` is refused by the face condition alone and
    # passed by admission's face check, which sets groups aside.
    function _testhelper_passes_only_set_aside(am, r, i)
        eps = D._eps_network(am, D._re_flags(am.catalytic), r)
        law, faces = D._law_network(am), D._faces(am.catalytic)
        v = D._face_condition(eps, law.present, faces[i], first(am.catalytic.substrates))
        v.fault === BR && D._admission_face(am, eps, law, faces, i) == D._WELL_POSED
    end
    # `set_aside`: a ping-pong whose chemistry steps and the bindings of A and of P are
    # :OnlyA. In the inactive conformation the law's forms are free E and EQ, and the
    # cut-off forms EA, FP, F and FB (F is E_res_plusA_minusP), which the excluded
    # steps join, are entered from free E by the A binding and from EQ by
    # FB ⇌ EQ, both losing their direction into them. Their vanishing steps open a
    # route from free E to EQ on which the reaction turns over, so the face condition
    # refuses the interior as a bridge. But the law's own rapid-equilibrium step
    # E + Q ⇌ EQ holds the two forms they touch in equilibrium, their fast classes drain,
    # and none fills faster than it empties: they are set aside.
    set_aside = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: OnlyA
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: OnlyA
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
        end
    end
    r = fill(D._RUNS, length(set_aside.catalytic.steps), 2)
    r[_testhelper_edge(set_aside, "E + A ⇌ EA"), 2] = D._LOSES_FORWARD
    r[_testhelper_edge(set_aside, "EA <--> EP_res_plusA_minusP"), 2] = D._LOSES_FORWARD
    r[_testhelper_edge(set_aside, "E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP"), 2] =
        D._LOSES_REVERSE
    r[_testhelper_edge(set_aside, "EB_res_plusA_minusP ⇌ EQ"), 2] = D._LOSES_REVERSE
    @test D._verdict_under(set_aside, r) == D._WELL_POSED
    @test _testhelper_admission_faults(set_aside; reading = _testhelper_reading(r)) ==
          fill(NO, 7)
    # The face condition alone refuses the interior and the face with A zero. There the
    # uptake E + A → EA is dead, so that pair runs one way, through its release, which
    # does not vanish and only empties EA: it is safe, and the cut-off forms go again.
    @test _testhelper_faces(set_aside)[4] == Set([:A])
    let eps = D._eps_network(set_aside, D._re_flags(set_aside.catalytic), r),
        law = D._law_network(set_aside), faces = D._faces(set_aside.catalytic)
        x = first(set_aside.catalytic.substrates)
        @test [D._face_condition(eps, law.present, f, x).fault for f in faces] ==
              [BR, NO, NO, BR, NO, NO, NO]
    end
    @test _testhelper_passes_only_set_aside(set_aside, r, 1)
    @test _testhelper_passes_only_set_aside(set_aside, r, 4)
    for face in _testhelper_faces(set_aside)
        @test _testhelper_epath_verdict(set_aside; rng, face, corners = true,
                                        reading = _testhelper_reading(r)) === :converges
    end
    # `kept`: the same ping-pong with E + Q <--> EQ steady-state, so free E and EQ are
    # no longer held in equilibrium, and EA ⇌ FP at rapid equilibrium, so the route
    # through the cut-off forms runs fast, its rates vanishing only as the reading says.
    # In the interior every transition at the cut-off forms runs both ways, and each of
    # them empties into the law's forms along fast steps that do not vanish, so none lies
    # in a fast class that could fail to drain; but the route joins two law forms whose
    # ratio the law sets otherwise: nothing is set aside, the interior stays a bridge, and
    # the limit there depends on how the two limits are taken.
    kept = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: OnlyA
            EA ⇌ EP_res_plusA_minusP                           :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
        end
    end
    k = fill(D._RUNS, length(kept.catalytic.steps), 2)
    k[_testhelper_edge(kept, "E + A ⇌ EA"), 2] = D._LOSES_FORWARD
    k[_testhelper_edge(kept, "EA ⇌ EP_res_plusA_minusP"), 2] = D._LOSES_FORWARD
    k[_testhelper_edge(kept, "E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP"), 2] =
        D._LOSES_REVERSE
    k[_testhelper_edge(kept, "EB_res_plusA_minusP ⇌ EQ"), 2] = D._LOSES_REVERSE
    eps = D._eps_network(kept, D._re_flags(kept.catalytic), k)
    law = D._law_network(kept)
    sets = D._face_sets(eps, law.present, UInt64(0))
    @test all(D._fast_law_reach(eps, law.present, sets).to_law[.!law.present])
    @test D._without_set_aside_groups(eps, law.present, UInt64(0)) === eps
    @test D._verdict_under(kept, k) ==
          D._Verdict(BR, U, 1, 2, _testhelper_edge(kept, "EB_res_plusA_minusP ⇌ EQ"))
    @test _testhelper_epath_verdict(kept; rng, corners = true,
                                    reading = _testhelper_reading(k)) ===
          :direction_dependent

    # One touching law form counts: in this ping-pong's inactive conformation free E stands
    # alone, and the cut-off forms EA, FP, F, FB and EQ touch nothing else.
    one_touch = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: OnlyA
            E + Q ⇌ EQ                                            :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: OnlyA
            EA <--> EP_res_plusA_minusP                        :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                           :: EqualAI
        end
    end
    o = fill(D._RUNS, length(one_touch.catalytic.steps), 2)
    o[_testhelper_edge(one_touch, "E + A ⇌ EA"), 2] = D._LOSES_FORWARD
    o[_testhelper_edge(one_touch, "E + Q ⇌ EQ"), 2] = D._LOSES_FORWARD
    o[_testhelper_edge(one_touch, "E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP"), 2] =
        D._LOSES_REVERSE
    @test D._verdict_under(one_touch, o) == D._WELL_POSED
    @test _testhelper_admission_faults(one_touch; reading = _testhelper_reading(o)) ==
          fill(NO, 7)
    # The face condition alone refuses five faces. With P zero the uptake F + P → FP is
    # dead, and F + P ⇌ FP runs one way between two cut-off forms, through its release,
    # which vanishes and is fast; it is safe because no fast transition from the law's
    # forms reaches its source FP, so the flow through it vanishes.
    @test _testhelper_faces(one_touch)[6] == Set([:P])
    let eps = D._eps_network(one_touch, D._re_flags(one_touch.catalytic), o),
        law = D._law_network(one_touch), faces = D._faces(one_touch.catalytic)
        x = first(one_touch.catalytic.substrates)
        @test [D._face_condition(eps, law.present, f, x).fault for f in faces] ==
              [BR, BR, NO, BR, NO, BR, BR]
    end
    @test _testhelper_passes_only_set_aside(one_touch, o, 1)
    @test _testhelper_passes_only_set_aside(one_touch, o, 6)
    for face in _testhelper_faces(one_touch)
        @test _testhelper_epath_verdict(one_touch; rng, face, corners = true,
                                        reading = _testhelper_reading(o)) === :converges
    end

    # `first_try`, a ter-ter mechanism of the gate's populations: under the guard's first
    # reading the face condition alone refuses five faces as bridges of the inactive
    # conformation, and admission passes every one by setting cut-off forms aside, so the
    # first try admits it. The oracles agree on every face under that reading, and the
    # compiled law equals the oracle's deletion limit.
    first_try = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A ⇌ EA, EDP + A ⇌ EADP, EP + A ⇌ EAP,
             EPR + A ⇌ EAPR) :: OnlyA
            (E + P ⇌ EP, EA + P ⇌ EAP, EAD + P ⇌ EADP,
             EAR + P ⇌ EAPR) :: OnlyA
            (EA + D ⇌ EAD, EAP + D ⇌ EADP, EAR + D ⇌ EADR,
             EP + D ⇌ EDP, EPR + D ⇌ EDPR) :: EqualAI
            (EA + R <--> EAR, EAD + R <--> EADR,
             EAP + R <--> EAPR, EDP + R <--> EDPR,
             EP + R <--> EPR) :: EqualAI
            EAD + B <--> EPQR :: OnlyA
            EPR + Q ⇌ EPQR :: OnlyA
        end
    end
    let reading = D._built_reading(first_try)
        eps = D._eps_network(first_try, D._re_flags(first_try.catalytic), reading)
        law, faces = D._law_network(first_try), D._faces(first_try.catalytic)
        x = first(first_try.catalytic.substrates)
        @test [D._face_condition(eps, law.present, f, x).fault for f in faces] ==
              [BR, NO, BR, BR, BR, BR, NO, NO, NO]
        @test D._verdict_under(first_try, reading) == D._WELL_POSED
    end
    @test D._guard(first_try) == D._WELL_POSED
    @test all(_testhelper_admission_agreement(first_try; rng = Random.MersenneTwister(1)))
    _testhelper_check_law_against_limit(first_try)

    # R2funnel and R2funnel_SS under the guard's first built reading, at zero P. Their
    # cut-off forms touch free E alone, and each reaches it along transitions that do not
    # vanish; but F and FS form a fast class whose one slow way out, to FP, leads back
    # into it along fast steps, while the slow EP <--> ES leaves EP, which a fast
    # step empties into F. The class does not drain, nothing is set aside, the face stays
    # a bridge, and the limit there depends on the direction.
    for name in ("R2funnel", "R2funnel_SS")
        am = _testhelper_tag_sets[name]
        @test _testhelper_faces(am)[2] == Set([:P])
        eps = D._eps_network(am, D._re_flags(am.catalytic), D._built_reading(am))
        law, faces = D._law_network(am), D._faces(am.catalytic)
        @test !any(D._face_sets(eps, law.present, faces[2]).stranded)
        @test D._without_set_aside_groups(eps, law.present, faces[2]) === eps
        @test D._admission_face(am, eps, law, faces, 2).fault === BR
        @test _testhelper_epath_verdict(am; rng, face = _testhelper_faces(am)[2],
                                        corners = true) === :direction_dependent
    end

    # The safety clause deciding. `stopped_exit` (a gate entry; the mechanism of
    # `slow_source` in "no form fills faster than it empties") under a reading that stops
    # E + A <--> FP, F being E_res_plusA_minusP, rather than letting it lose its
    # uptake. In the inactive conformation free E stands alone, so the cut-off forms
    # touch one law form, and their fast classes drain through FP + Q <--> EAQ,
    # EAQ → EQ + A and EQ → E + Q. With A and B zero the stopped step runs one
    # way, from FP into free E, along its vanishing direction: a pair at a law form with
    # only its vanishing direction alive, which is not safe. It alone keeps them: nothing
    # is set aside, the face stays a bridge, and the limit there depends on the direction.
    stopped_exit = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EP_res_plusA_minusP, EQ + A <--> EAQ)    :: OnlyA
            (E + Q <--> EQ, EP_res_plusA_minusP + Q <--> EAQ)    :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP            :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                              :: OnlyA
        end
    end
    se = fill(D._RUNS, length(stopped_exit.catalytic.steps), 2)
    se[_testhelper_edge(stopped_exit, "E + A <--> EP_res_plusA_minusP"), 2] = D._STOPPED
    for text in ("E + Q <--> EQ", "EQ + A <--> EAQ")
        se[_testhelper_edge(stopped_exit, text), 2] = D._LOSES_FORWARD
    end
    for text in ("E_res_plusA_minusP + B ⇌ EQ", "EP_res_plusA_minusP + Q <--> EAQ")
        se[_testhelper_edge(stopped_exit, text), 2] = D._LOSES_REVERSE
    end
    @test _testhelper_faces(stopped_exit)[3] == Set([:A, :B])
    let eps = D._eps_network(stopped_exit, D._re_flags(stopped_exit.catalytic), se),
        law = D._law_network(stopped_exit), faces = D._faces(stopped_exit.catalytic)
        @test count(law.present[(length(stopped_exit.catalytic.enzyme_forms) + 1):end]) == 1
        @test !any(D._face_sets(eps, law.present, faces[3]).stranded)
        @test D._without_set_aside_groups(eps, law.present, faces[3]) === eps
    end
    @test D._verdict_under(stopped_exit, se) ==
          D._Verdict(BR, U, 3, 2, _testhelper_edge(stopped_exit, "E + Q <--> EQ"))
    @test _testhelper_epath_verdict(stopped_exit; rng,
                                    face = _testhelper_faces(stopped_exit)[3],
                                    corners = true,
                                    reading = _testhelper_reading(se)) ===
          :direction_dependent

    # The conservative rejections the suite pins stay rejected under the guard's first
    # built reading. R-B1's cut-off forms EI and EIJ touch free E and the dead end
    # EJ, which only the steady-state E + J <--> EJ joins; R-B2's basins reach the
    # law's forms only through vanishing steps, so they are stranded. (The search admits
    # R-B2 under other readings: "conservative rejections the oracle admits stay rejected
    # under their reading".)
    for (name, fault) in (("bridge_deadend", BR), ("BS1", D._PUMPED),
                          ("fz101_10232", D._PUMPED))
        am = _testhelper_tag_sets[name]
        @test (name, D._verdict_under(am, D._built_reading(am)).fault) == (name, fault)
    end
    @test D._guard(_testhelper_tag_sets["bridge_deadend"]).fault === BR
end

@testset "the guard searches for a reading" begin
    D = EnzymeRates.Derivation
    rng = Random.MersenneTwister(2652)
    # A six-step ping-pong with both chemistry steps :OnlyA (F is
    # E_res_plusA_minusP). No built reading switches off FB ⇌ EQ, a
    # rapid-equilibrium step that takes up nothing, so the guard searches. The first
    # reading it finds loses EA <--> FP's stored direction and FB ⇌ EQ's reverse,
    # both into the covalent group: the reading of "a reading may lose a step in either
    # direction".
    pingpong = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
        end
    end
    pingpong_dimer = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
        end
    end
    chem_a = _testhelper_edge(pingpong, "EA <--> EP_res_plusA_minusP")
    chem_b = _testhelper_edge(pingpong, "EB_res_plusA_minusP ⇌ EQ")
    r = fill(D._RUNS, length(pingpong.catalytic.steps), 2)
    r[chem_a, 2], r[chem_b, 2] = D._LOSES_FORWARD, D._LOSES_REVERSE
    @test D._built_reading(pingpong)[chem_b, 2] === D._CANNOT_SWITCH_OFF
    for am in (pingpong, pingpong_dimer)
        @test D._guard_with_reading(am) == (; verdict = D._WELL_POSED, reading = r)
        _testhelper_check_law_against_limit(am)
        @test all(_testhelper_admission_agreement(am; rng, corners = true))
    end
    # Its mirror image tags both chemistry steps :OnlyI, so the active conformation keeps
    # no chemistry: no reading helps.
    @test D._guard(_testhelper_mirror(pingpong)).fault === D._NO_ACTIVE_TURNOVER
    # The verdict does not depend on how the steps are written: the macro stores them in
    # canonical order. Nor on the order in which they are stored, which only orders the
    # search: stored backward, the same reading admits the mechanism.
    pingpong_rewritten = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E + A ⇌ EA                                            :: EqualAI
        end
    end
    @test pingpong_rewritten == pingpong
    @test D._guard(pingpong_rewritten) == D._WELL_POSED
    let cm = pingpong.catalytic, back = reverse(eachindex(cm.steps))
        backward = EnzymeRates.Mechanism(
            EnzymeRates.Mechanisms._unchecked_scheme(cm.metabolites, cm.enzyme_forms,
                                                     cm.substrates, cm.products,
                                                     cm.steps[back]),
            EnzymeRates.Mechanisms.Allostery(pingpong.allosteric.allosteric_states, 1,
                                             pingpong.allosteric.regulatory_sites))
        @test D._guard_with_reading(backward) ==
              (; verdict = D._WELL_POSED, reading = r[back, :])
    end

    # The three-step ping-pong: the inactive conformation keeps only E + Q <--> EQ. The
    # built readings stop both chemistry steps, which leaves F out of free E's reach; the
    # first reading the search finds loses the first half's stored direction and the
    # second half's reverse, both into F.
    three_step = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> F + P        :: OnlyA
            E + Q <--> EQ         :: EqualAI
            F + B <--> EQ         :: OnlyA
        end
    end
    half_a = _testhelper_edge(three_step, "E + A <--> F + P")
    half_b = _testhelper_edge(three_step, "F + B <--> EQ")
    t = fill(D._RUNS, length(three_step.catalytic.steps), 2)
    t[half_a, 2], t[half_b, 2] = D._LOSES_FORWARD, D._LOSES_REVERSE
    @test D._admission_interior(three_step, D._built_reading(three_step)).fault ===
          D._ISLAND
    @test D._guard_with_reading(three_step) == (; verdict = D._WELL_POSED, reading = t)
    _testhelper_check_law_against_limit(three_step)
    @test all(_testhelper_admission_agreement(three_step; rng, corners = true))

    # A ping-pong with no limit under any reading: the covalent forms empty only by
    # binding B. The guard finds no reading, and the ε-path oracle converges on every face
    # under none of the consistent readings.
    no_limit = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EP_res_plusA_minusP                       :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                        :: OnlyA
        end
    end
    @test D._guard(no_limit).fault !== D._NO_FAULT
    consistent = _testhelper_consistent_readings(no_limit)
    @test length(consistent) == 2
    for reading in consistent
        @test any(face -> _testhelper_epath_verdict(no_limit; rng, face, corners = true,
                                                    reading) !== :converges,
                  _testhelper_faces(no_limit))
    end

    # `filled` of "no form fills faster than it empties": its other readings fail too.
    filled = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: EqualAI
        end
    end
    @test D._guard(filled).fault !== D._NO_FAULT

    # A mechanism a built reading admits is admitted with it: `both_readings` of "a
    # mechanism may need one uptake to lose its affinity and another to stop".
    both_readings = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP)  :: OnlyA
            E + P <--> EP                     :: OnlyA
            EA + P <--> EAP               :: EqualAI
            EA + B <--> EPQ               :: OnlyA
            EP + Q ⇌ EPQ                  :: EqualAI
        end
    end
    @test D._guard_with_reading(both_readings) ==
          (; verdict = D._WELL_POSED,
           reading = D._built_reading(both_readings, D._OUTWARD_LOSE_AFFINITY))

    # A rapid-equilibrium step inside one group is the fault no reading changes. The
    # message offers the two fixes, each of whose variants passes.
    pinned = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
    end
    v = D._guard(pinned)
    @test v == D._Verdict(D._FAST_STEP_PINNED, D._LIMIT_UNTESTED, 0, 2,
                          _testhelper_edge(pinned, "E + S ⇌ ES"))
    @test D._guard_message(pinned) ==
          "Mechanism: in the inactive conformation the rapid-equilibrium step " *
          "E + S ⇌ ES is tagged :OnlyA, but the steps that conformation keeps still " *
          "join its two sides. Thermodynamics then fixes its equilibrium constant, and a " *
          "step keeps its kind in both conformations, so it cannot be switched off. Make " *
          "its kinetic group steady-state, or tag :OnlyA the steps that join its sides: " *
          "E + P ⇌ EP; ES <--> EP."
    @test _testhelper_consistent_readings(pinned) == Matrix{Symbol}[]

    # A cycle of rapid-equilibrium steps that turns the reaction over, each conformation
    # excluding two of them between its two groups. The built readings cannot switch off
    # the isomerizations, so under them the cycle is slow in each conformation; every
    # consistent reading lets each excluded step lose one direction, so the cycle runs
    # fast in both, and class 2 rejects it under the first, naming how each excluded step
    # of the cycle is switched off. The ε-path oracle finds the rate diverging.
    fast_cycle = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyI
            ES ⇌ EstarS     :: OnlyA
            EstarS ⇌ EP     :: OnlyI
            E + P ⇌ EP        :: OnlyA
        end
    end
    @test D._verdict_under(fast_cycle, D._built_reading(fast_cycle)).fault ===
          D._FAST_STEP_PINNED
    v = D._guard(fast_cycle)
    @test v.fault === D._RE_TURNOVER
    reported = _testhelper_reading(D._guard_with_reading(fast_cycle).reading)
    @test _testhelper_reading_is_consistent(fast_cycle, reported)
    msg = D._guard_message(fast_cycle)
    @test occursin("(whose affinity vanishes in the ", msg)
    @test occursin("(whose rate one way vanishes in the ", msg)
    @test _testhelper_epath_verdict(fast_cycle; rng, reading = reported) === :diverges

    # A mechanism and its mirror image, every :OnlyA and :OnlyI swapped, get one verdict
    # when both conformations catalyze: the six-step ping-pong with a ternary-complex path
    # beside it, which both conformations keep. The search admits both, the mirror with
    # the same reading in the other conformation, and the oracles agree on each.
    pp_ternary = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: EqualAI
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
            EA + B ⇌ EAB                                      :: EqualAI
            EAB <--> EPQ                                    :: EqualAI
            EQ + P ⇌ EPQ                                      :: EqualAI
        end
    end
    mirror = _testhelper_mirror(pp_ternary)
    @test D._admission_interior(pp_ternary, D._built_reading(pp_ternary)).fault ===
          D._FAST_STEP_PINNED
    p = fill(D._RUNS, length(pp_ternary.catalytic.steps), 2)
    p[_testhelper_edge(pp_ternary, "EA <--> EP_res_plusA_minusP"), 2] = D._LOSES_FORWARD
    p[_testhelper_edge(pp_ternary, "EB_res_plusA_minusP ⇌ EQ"), 2] = D._LOSES_REVERSE
    @test D._guard_with_reading(pp_ternary) == (; verdict = D._WELL_POSED, reading = p)
    @test D._guard_with_reading(mirror) ==
          (; verdict = D._WELL_POSED, reading = p[:, [2, 1]])
    for am in (pp_ternary, mirror)
        @test D._guard(am) == D._WELL_POSED
        _testhelper_check_law_against_limit(am)
        @test all(_testhelper_admission_agreement(am; rng, corners = true))
    end

    # The search is skipped whole for a mechanism too large to search: three ter-ter
    # taggings of the gate's populations whose inactive conformation has more ways to
    # switch off its excluded steps than the guard searches, counted by the tests' own
    # rule. The built readings admit `large_admitted`. They reject `large_island` and
    # `large_pumped`, and so does the guard, though a reading it never tries passes
    # every check: under no built reading are `large_island`'s forms all reached, so the
    # guard reports that there are too many readings; `large_pumped`'s first built reading
    # is consistent, so the guard reports its verdict.
    # The ways conformation `c` of `am` may switch off its excluded steps between two
    # groups, by the tests' own rule: each loses toward either side but free E's group,
    # or, at steady state, stops.
    function _testhelper_choice_count(am, c = 2)
        group = _testhelper_kept_components(am, c)
        prod((group[e.to] != group[1]) + (group[e.from] != group[1]) + !e.rapid_equilibrium
             for (k, e) in enumerate(am.catalytic.steps)
             if _testhelper_excludes(am, k, c) && group[e.from] != group[e.to]; init = 1)
    end
    large_admitted = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + D ⇌ ED :: EqualAI
            E + R ⇌ ER :: OnlyA
            (EAR + P <--> EAPR, EQR + P <--> EPQR,
             ER + P <--> EPR) :: OnlyA
            (EBD + A <--> EPQR, EPR + A <--> EAPR,
             ER + A <--> EAR) :: OnlyA
            (EBR + Q ⇌ EBQR, EPR + Q ⇌ EPQR,
             ER + Q ⇌ EQR) :: OnlyA
            (ED + B ⇌ EBD, EQR + B ⇌ EBQR, ER + B ⇌ EBR) :: OnlyA
        end
    end
    large_island = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ) :: OnlyA
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ,
             EP + D <--> EDP,
             EQ_res_plusB_plusD_minusQ + D <--> EDQ_res_plusB_plusD_minusQ) :: OnlyA
            (E + P <--> EP, ED + P <--> EDP, ER + P <--> EPR) :: OnlyA
            (E + R <--> ER, EP + R <--> EPR) :: EqualAI
            E_res_plusB_plusD_minusQ + A ⇌ EPR :: OnlyA
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             ED + Q ⇌ EDQ) :: OnlyA
        end
    end
    large_pumped = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A <--> EA, EB + A <--> EP_res_plusA_plusB_minusP,
             ER + A <--> EAR) :: OnlyA
            (E + B <--> EB, EA + B <--> EP_res_plusA_plusB_minusP) :: OnlyA
            (E + R <--> ER, EA + R <--> EAR,
             EQR + R <--> EDR_res_plusA_plusB_minusP) :: OnlyA
            (E_res_plusA_plusB_minusP + P ⇌ EP_res_plusA_plusB_minusP,
             EQR + P ⇌ EDP_res_plusA_plusB_minusP) :: EqualAI
            (E_res_plusA_plusB_minusP + D <--> EQR,
             EP_res_plusA_plusB_minusP + D <--> EDP_res_plusA_plusB_minusP,
             ER + D <--> EDR) :: OnlyA
            ER + Q ⇌ EQR :: OnlyA
        end
    end
    for am in (large_admitted, large_island, large_pumped)
        @test _testhelper_choice_count(am) > D._MAX_CHOICES_PER_CONFORMATION
    end
    @test D._guard_with_reading(large_admitted) ==
          (; verdict = D._WELL_POSED, reading = D._built_reading(large_admitted))
    @test D._guard(large_island) ==
          D._Verdict(D._TOO_MANY_READINGS, D._LIMIT_UNTESTED, 0, 2, 0)
    @test D._guard_message(large_island) ==
          "Mechanism: the inactive conformation can switch off the tagged " *
          "steps in more than 20000 ways, too many for the guard to search, and none of " *
          "the three ways it tries first fits the rule for switching them off, so it " *
          "cannot show that the rate law is the limit."
    @test D._guard(large_pumped).fault === D._PUMPED
    @test D._guard(large_pumped) ==
          D._verdict_under(large_pumped, D._built_reading(large_pumped))
    names = Dict(name => o for (o, name) in _testhelper_switch_off_names)
    LF, LR = :loses_forward, :loses_reverse
    for (am, column) in ((large_island, [LF, LF, LF, :runs, LF, LF, LF, LF, LF, :runs, LF,
                                         LF, LR, LR]),
                         (large_pumped, [LF, LF, LF, LF, LF, LF, LR, :runs, LF, LF, LF, LF,
                                         :runs, LR]))
        reading = [fill(:runs, length(column)) column]
        @test _testhelper_reading_is_consistent(am, reading)
        @test D._verdict_under(am, map(name -> names[name], reading)) == D._WELL_POSED
        @test D._guard(am) != D._WELL_POSED
    end

    # The two other ways a search can be too large, each first in a ter-ter tagging that
    # mixes :OnlyA and :OnlyI, found by sampling the gate's ter-ter init mechanisms. In
    # `too_many_pairs` each conformation stays within `_MAX_CHOICES_PER_CONFORMATION`, but
    # its consistent columns, counted by the tests' own rule and equal to the package's,
    # pair up more than `_MAX_READING_PAIRS` times; in `too_many_active` it is the active
    # conformation that has too many choices. No built reading of either is consistent.
    # Neither active conformation turns the reaction over, which no reading changes, so
    # the guard reports that and not the cap: the cap does not hide it.
    function _testhelper_consistent_column_count(am, c)
        n = length(am.catalytic.steps)
        group = _testhelper_kept_components(am, c)
        between = [k for (k, e) in enumerate(am.catalytic.steps)
                   if _testhelper_excludes(am, k, c) && group[e.from] != group[e.to]]
        r = fill(:runs, n, 2)
        for (k, e) in enumerate(am.catalytic.steps)
            _testhelper_excludes(am, k, c) && group[e.from] == group[e.to] &&
                (r[k, c] = e.rapid_equilibrium ? :cannot_switch_off : :stopped)
        end
        ways(k) = [w for w in (:loses_forward, :loses_reverse, :stopped)
                   if !(w === :loses_forward &&
                        group[am.catalytic.steps[k].to] == group[1]) &&
                      !(w === :loses_reverse &&
                        group[am.catalytic.steps[k].from] == group[1]) &&
                      !(w === :stopped && am.catalytic.steps[k].rapid_equilibrium)]
        count(Iterators.product((ways(k) for k in between)...)) do choice
            r[between, c] .= collect(choice)
            ex = _testhelper_exclusion(am, c; reading = r)
            ex.failure === :none && isempty(ex.lost)
        end
    end
    too_many_pairs = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A <--> EA, ED + A <--> EP_res_plusA_plusD_minusP,
             EQ + A <--> EAQ) :: OnlyI
            (E + D <--> ED, EA + D <--> EP_res_plusA_plusD_minusP,
             ER + D <--> EDR) :: OnlyI
            (E + Q <--> EQ, EA + Q <--> EAQ,
             EQR + Q <--> EBQ_res_plusA_plusD_minusP, ER + Q <--> EQR) :: OnlyI
            (E + R <--> ER, ED + R <--> EDR, EQ + R <--> EQR) :: OnlyA
            (E_res_plusA_plusD_minusP + P ⇌ EP_res_plusA_plusD_minusP,
             EQR + P ⇌ EBP_res_plusA_plusD_minusP) :: OnlyA
            (E_res_plusA_plusD_minusP + B ⇌ EQR,
             EP_res_plusA_plusD_minusP + B ⇌ EBP_res_plusA_plusD_minusP,
             EQ + B ⇌ EBQ) :: EqualAI
        end
    end
    too_many_active = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A ⇌ EA, EBR + A ⇌ EABR, ER + A ⇌ EAR) :: OnlyI
            (E + R <--> ER, EA + R <--> EAR,
             EQ_res_plusA_plusB_minusQ + R <--> EABR) :: OnlyI
            E_res_plusA_plusB_minusQ + D <--> EPR :: OnlyI
            E_res_plusA_plusB_minusQ + Q ⇌ EQ_res_plusA_plusB_minusQ :: OnlyI
            (EA + B <--> EQ_res_plusA_plusB_minusQ, EAR + B <--> EABR,
             ER + B <--> EBR) :: OnlyI
            ER + P ⇌ EPR :: OnlyI
        end
    end
    pair_columns = [_testhelper_consistent_column_count(too_many_pairs, c) for c in 1:2]
    @test pair_columns == [length(D._consistent_columns(too_many_pairs, c).columns)
                           for c in 1:2]
    @test all(c -> _testhelper_choice_count(too_many_pairs, c) <=
                   D._MAX_CHOICES_PER_CONFORMATION, 1:2)
    @test prod(pair_columns) > D._MAX_READING_PAIRS
    @test _testhelper_choice_count(too_many_active, 1) > D._MAX_CHOICES_PER_CONFORMATION
    for am in (too_many_pairs, too_many_active)
        @test !any(r -> _testhelper_reading_is_consistent(
                       am, _testhelper_reading(am, r)), keys(_testhelper_guard_readings))
    end
    for am in (too_many_pairs, too_many_active)
        @test !EnzymeRates.Mechanisms.has_turnover_cycle(
            am.catalytic,
            EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(am).active)
        @test D._guard(am) == D._Verdict(D._NO_ACTIVE_TURNOVER, D._LIMIT_UNTESTED, 0, 1, 1)
    end
    @test D._guard_message(too_many_pairs) ==
          "Mechanism: without its :OnlyI steps E + A <--> EA; " *
          "E + D <--> ED; E + Q <--> EQ; EA + Q <--> EAQ; " *
          "EA + D <--> EP_res_plusA_plusD_minusP; " *
          "ED + A <--> EP_res_plusA_plusD_minusP; EQ + A <--> EAQ; " *
          "EQR + Q <--> EBQ_res_plusA_plusD_minusP; ER + D <--> EDR; " *
          "ER + Q <--> EQR the active conformation turns no substrate into " *
          "product, but it must catalyze the reaction."
    @test D._guard_message(too_many_active) ==
          "Mechanism: without its :OnlyI steps E + A ⇌ EA; E + R <--> ER; " *
          "EA + R <--> EAR; EA + B <--> EQ_res_plusA_plusB_minusQ; " *
          "EAR + B <--> EABR; EBR + A ⇌ EABR; " *
          "EQ_res_plusA_plusB_minusQ + R <--> EABR; ER + A ⇌ EAR; " *
          "ER + B <--> EBR; ER + P ⇌ EPR; " *
          "E_res_plusA_plusB_minusQ + D <--> EPR; " *
          "E_res_plusA_plusB_minusQ + Q ⇌ EQ_res_plusA_plusB_minusQ the active " *
          "conformation turns no " *
          "substrate into product, but it must catalyze the reaction. Tag their kinetic " *
          "groups :EqualAI or :NonequalAI."

    # The same two ways with an active conformation that turns the reaction over, built
    # on random-order ter-ter: one full path of bindings, chemistry and releases kept by
    # both conformations; the other nine bindings, a second chemistry through EstarABD
    # and a second release order of Q and R tagged. In `cap_active` all three are :OnlyI,
    # which leaves the active conformation more than `_MAX_CHOICES_PER_CONFORMATION` ways;
    # in `cap_pairs` the bindings are :OnlyI and the rest :OnlyA, and the two
    # conformations' consistent columns pair up more than `_MAX_READING_PAIRS` times. Each
    # built reading leaves a group that no switched-off step reaches (EstarABD under
    # the two that lose affinity, EB or EQ when all stop), so none is consistent, and
    # the guard reports the cap.
    cap_active = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + A <--> EA                 :: EqualAI
            EA + B <--> EAB           :: EqualAI
            EAB + D <--> EABD     :: EqualAI
            E + B <--> EB                 :: OnlyI
            E + D <--> ED                 :: OnlyI
            EA + D <--> EAD           :: OnlyI
            EB + A <--> EAB           :: OnlyI
            EB + D <--> EBD           :: OnlyI
            ED + A <--> EAD           :: OnlyI
            ED + B <--> EBD           :: OnlyI
            EAD + B <--> EABD     :: OnlyI
            EBD + A <--> EABD     :: OnlyI
            EABD <--> EPQR      :: EqualAI
            EABD <--> EstarABD  :: OnlyI
            EstarABD <--> EPQR  :: OnlyI
            EPQR <--> EQR + P     :: EqualAI
            EQR <--> ER + Q           :: EqualAI
            ER <--> E + R                 :: EqualAI
            EQR <--> EQ + R           :: OnlyI
            EQ <--> E + Q                 :: OnlyI
        end
    end
    cap_pairs = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + A <--> EA                 :: EqualAI
            EA + B <--> EAB           :: EqualAI
            EAB + D <--> EABD     :: EqualAI
            E + B <--> EB                 :: OnlyI
            E + D <--> ED                 :: OnlyI
            EA + D <--> EAD           :: OnlyI
            EB + A <--> EAB           :: OnlyI
            EB + D <--> EBD           :: OnlyI
            ED + A <--> EAD           :: OnlyI
            ED + B <--> EBD           :: OnlyI
            EAD + B <--> EABD     :: OnlyI
            EBD + A <--> EABD     :: OnlyI
            EABD <--> EPQR      :: EqualAI
            EABD <--> EstarABD  :: OnlyA
            EstarABD <--> EPQR  :: OnlyA
            EPQR <--> EQR + P     :: EqualAI
            EQR <--> ER + Q           :: EqualAI
            ER <--> E + R                 :: EqualAI
            EQR <--> EQ + R           :: OnlyA
            EQ <--> E + Q                 :: OnlyA
        end
    end
    @test _testhelper_choice_count(cap_active, 1) > D._MAX_CHOICES_PER_CONFORMATION
    cap_columns = [_testhelper_consistent_column_count(cap_pairs, c) for c in 1:2]
    @test cap_columns ==
          [length(D._consistent_columns(cap_pairs, c).columns) for c in 1:2]
    @test all(c -> _testhelper_choice_count(cap_pairs, c) <=
                   D._MAX_CHOICES_PER_CONFORMATION, 1:2)
    @test prod(cap_columns) > D._MAX_READING_PAIRS
    for am in (cap_active, cap_pairs)
        @test EnzymeRates.Mechanisms.has_turnover_cycle(
            am.catalytic,
            EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(am).active)
        @test !any(r -> _testhelper_reading_is_consistent(
                       am, _testhelper_reading(am, r)), keys(_testhelper_guard_readings))
    end
    @test D._guard(cap_active) ==
          D._Verdict(D._TOO_MANY_READINGS, D._LIMIT_UNTESTED, 0, 1, 0)
    @test D._guard_message(cap_active) ==
          "Mechanism: the active conformation can switch off the tagged steps " *
          "in more than 20000 ways, too many for the guard to search, and none of the " *
          "three ways it tries first fits the rule for switching them off, so it cannot " *
          "show that the rate law is the limit."
    @test D._guard(cap_pairs) ==
          D._Verdict(D._TOO_MANY_READINGS, D._LIMIT_UNTESTED, 0, 0, 0)
    @test D._guard_message(cap_pairs) ==
          "Mechanism: the two conformations together can switch off the tagged " *
          "steps in more than 5000 ways, too many for the guard to search, and none of " *
          "the three ways it tries first fits the rule for switching them off, so it " *
          "cannot show that the rate law is the limit."
end

@testset "the graph test equals the exact LP on the uni-uni grid" begin
    # A limit exists exactly when no excluded rapid-equilibrium step is pinned (the reading
    # cannot switch it off, or stops it, though a step keeps its kind) and the
    # steps that lose a direction, each pointing uphill, close no cycle and never point
    # into free E's component; the test oracle's own reading of the rule and its exact
    # simplex over the cycle rows decide the same question, and its exclusion structure
    # names the same failure. So under the guard's first built reading, and under every
    # reading of the inactive conformation's excluded steps, the search's among them.
    grid = _testhelper_uniuni_grid()
    @test length(grid) == 108
    for am in grid
        @test _testhelper_existence_agrees(am)
        n = length(am.catalytic.steps)
        excluded = [k for k in 1:n if _testhelper_excludes(am, k, 2)]
        ways = (:loses_forward, :loses_reverse, :stopped, :cannot_switch_off)
        readings = map(Iterators.product(fill(ways, length(excluded))...)) do o
            r = fill(:runs, n, 2)
            r[excluded, 2] .= collect(o)
            r
        end
        @test all(r -> _testhelper_existence_agrees(am; reading = r), readings)
    end
end

@testset "a rapid-equilibrium step whose sides stay joined cannot be switched off" begin
    rng = Random.MersenneTwister(31)
    FP, U = EnzymeRates.Derivation._FAST_STEP_PINNED, EnzymeRates.Derivation._LIMIT_UNTESTED
    # The uni-uni S binding :OnlyA: the chemistry and the P binding join E to ES in the
    # inactive conformation, so the S binding cannot lose its affinity there, and at rapid
    # equilibrium it cannot be stopped either.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
    end
    m = am.catalytic
    bind, release, chem = (_testhelper_edge(am, t) for t in
                           ("E + S ⇌ ES", "E + P ⇌ EP", "ES <--> EP"))
    # The chemistry is stored substrate side first, so it runs toward products as stored.
    @test m.enzyme_forms[[m.steps[chem].from, m.steps[chem].to]] == [:ES, :EP]
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(FP, U, 0, 2, bind)
    @test EnzymeRates.Derivation._admission_interior(
        am, EnzymeRates.Derivation._built_reading(am)) == v
    @test _testhelper_epath_verdict(am; rng) === :no_path
    line(k) = EnzymeRates.Mechanisms.step_text(m, k)
    @test EnzymeRates.Derivation._guard_message(am) ==
          _testhelper_pinned_cause(am, bind, 2) * " Make its kinetic group steady-state, " *
          "or tag :OnlyA the steps that join its sides: $(line(release)); $(line(chem))."
    # The fixes pass the guard: the steps that join its sides tagged, and the binding
    # steady-state ("how an excluded step is switched off"); so does the binding untagged,
    # which the message does not offer.
    untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
    end
    joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: OnlyA
            E + P ⇌ EP        :: OnlyA
        end
    end
    @test EnzymeRates.Derivation._guard(untagged) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(joins_tagged) == EnzymeRates.Derivation._WELL_POSED
    # AD1 and AD11: a mixed :OnlyA group that binds B into EPQ and into EBQ, forms
    # the untagged steps still make. At steady state the group is stopped, and the law is
    # the limit on every face; at rapid equilibrium (AD1 RE) it cannot be switched off.
    # AD3c: a ::Mod copy's rapid-equilibrium binding whose sides untagged steps join.
    ad1_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                    :: EqualAI
            (EA + B <--> EPQ, EQ + B <--> EBQ)  :: OnlyA
            EQ + P ⇌ EPQ                              :: EqualAI
            E + Q ⇌ EQ                                    :: EqualAI
            E + B ⇌ EB                                    :: EqualAI
            EB + Q ⇌ EBQ                              :: EqualAI
        end
    end
    ad1_re = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)    :: OnlyA
            EQ + P <--> EPQ                       :: EqualAI
            E + Q ⇌ EQ                                :: EqualAI
            E + B ⇌ EB                                :: EqualAI
            EB + Q ⇌ EBQ                          :: EqualAI
        end
    end
    ad11 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            (EA + B <--> EPQ, EQ + B <--> EBQ,
             EP + B <--> EBP)                                 :: OnlyA
            EQ + P ⇌ EPQ                                      :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E + P ⇌ EP                                            :: EqualAI
            E + B ⇌ EB                                            :: EqualAI
            EB + Q ⇌ EBQ                                      :: EqualAI
        end
    end
    ad3c = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES                    :: EqualAI
            ES <--> EP                  :: OnlyA
            E + P ⇌ EP                    :: EqualAI
            E + S ⇌ ESmod          :: OnlyA
            ES + S ⇌ ESSmod    :: EqualAI
            ESmod + S ⇌ ESSmod    :: EqualAI
        end
    end
    for x in (ad1_ss, ad11)
        @test EnzymeRates.Derivation._guard(x) == EnzymeRates.Derivation._WELL_POSED
        @test all(_testhelper_admission_agreement(x; rng))
    end
    for (x, text) in ((ad1_re, "EA + B ⇌ EPQ"), (ad3c, "E + S ⇌ ESmod"))
        @test EnzymeRates.Derivation._guard(x) ==
              EnzymeRates.Derivation._Verdict(FP, U, 0, 2, _testhelper_edge(x, text))
        @test _testhelper_exclusion(x, 2).failure === :fast
        @test _testhelper_epath_verdict(x; rng) === :no_path
    end
    # On AD3c the message's steps that join the sides run from the binding's unbound side
    # to its bound side, and that fix passes there. On AD1 RE both fixes pass.
    message(x) = EnzymeRates.Derivation._guard_message(x)
    joins = join((EnzymeRates.Mechanisms.step_text(ad3c.catalytic,
                                                   _testhelper_edge(ad3c, t))
                  for t in ("E + S ⇌ ES", "ES + S ⇌ ESSmod", "ESmod + S ⇌ ESSmod")),
                 "; ")
    @test endswith(message(ad3c), "Tag :OnlyA the steps that join its sides: $joins.")
    ad1_re_joins = join((EnzymeRates.Mechanisms.step_text(ad1_re.catalytic,
                                                         _testhelper_edge(ad1_re, t))
                         for t in ("E + A ⇌ EA", "E + Q ⇌ EQ", "EQ + P <--> EPQ")), "; ")
    @test endswith(message(ad1_re), "Make its kinetic group steady-state, or tag :OnlyA " *
                                    "the steps that join its sides: $ad1_re_joins.")
    ad1_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                    :: EqualAI
            (EA + B <--> EPQ, EQ + B <--> EBQ)  :: EqualAI
            EQ + P ⇌ EPQ                              :: EqualAI
            E + Q ⇌ EQ                                    :: EqualAI
            E + B ⇌ EB                                    :: EqualAI
            EB + Q ⇌ EBQ                              :: EqualAI
        end
    end
    ad1_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                    :: EqualAI
            (EA + B <--> EPQ, EQ + B <--> EBQ)  :: OnlyA
            EQ + P ⇌ EPQ                              :: EqualAI
            E + Q ⇌ EQ                                    :: OnlyA
            E + B ⇌ EB                                    :: OnlyA
            EB + Q ⇌ EBQ                              :: OnlyA
        end
    end
    ad1_re_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)    :: EqualAI
            EQ + P <--> EPQ                       :: EqualAI
            E + Q ⇌ EQ                                :: EqualAI
            E + B ⇌ EB                                :: EqualAI
            EB + Q ⇌ EBQ                          :: EqualAI
        end
    end
    ad1_re_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)    :: OnlyA
            EQ + P <--> EPQ                       :: EqualAI
            E + Q ⇌ EQ                                :: OnlyA
            E + B ⇌ EB                                :: OnlyA
            EB + Q ⇌ EBQ                          :: OnlyA
        end
    end
    ad3c_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES                    :: EqualAI
            ES <--> EP                  :: OnlyA
            E + P ⇌ EP                    :: EqualAI
            E + S ⇌ ESmod          :: EqualAI
            ES + S ⇌ ESSmod    :: EqualAI
            ESmod + S ⇌ ESSmod    :: EqualAI
        end
    end
    ad3c_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES                    :: OnlyA
            ES <--> EP                  :: OnlyA
            E + P ⇌ EP                    :: EqualAI
            E + S ⇌ ESmod          :: OnlyA
            ES + S ⇌ ESSmod    :: OnlyA
            ESmod + S ⇌ ESSmod    :: OnlyA
        end
    end
    # Each variant of AD1, AD1 RE, AD3c and AD11 with the group untagged, or with the steps
    # that join the sides of its B or S binding tagged, passes the guard.
    ad11_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            (EA + B <--> EPQ, EQ + B <--> EBQ,
             EP + B <--> EBP)                                 :: EqualAI
            EQ + P ⇌ EPQ                                      :: EqualAI
            E + Q ⇌ EQ                                            :: EqualAI
            E + P ⇌ EP                                            :: EqualAI
            E + B ⇌ EB                                            :: EqualAI
            EB + Q ⇌ EBQ                                      :: EqualAI
        end
    end
    ad11_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            (EA + B <--> EPQ, EQ + B <--> EBQ,
             EP + B <--> EBP)                                 :: OnlyA
            EQ + P ⇌ EPQ                                      :: EqualAI
            E + Q ⇌ EQ                                            :: OnlyA
            E + P ⇌ EP                                            :: EqualAI
            E + B ⇌ EB                                            :: OnlyA
            EB + Q ⇌ EBQ                                      :: OnlyA
        end
    end
    for y in (ad1_untagged, ad1_joins_tagged, ad1_re_untagged, ad1_re_joins_tagged,
              ad3c_untagged, ad3c_joins_tagged, ad11_untagged, ad11_joins_tagged)
        @test EnzymeRates.Derivation._guard(y) == EnzymeRates.Derivation._WELL_POSED
    end
    # The message prints a fix only when its variant passes. A dead-end I joins the S
    # binding's sides a second way, through EI and EIS, so tagging the steps of the
    # path the message would name leaves that binding inside one component, and the guard
    # rejects its steady-state variant too: the message offers no fix. (Untagging the
    # binding passes, but is not one of the two fixes the message offers.)
    deadend = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
            E + I ⇌ EI        :: EqualAI
            ES + I ⇌ EIS  :: EqualAI
            EI + S ⇌ EIS  :: EqualAI
        end
    end
    deadend_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
            E + I ⇌ EI        :: EqualAI
            ES + I ⇌ EIS  :: EqualAI
            EI + S ⇌ EIS  :: EqualAI
        end
    end
    deadend_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: OnlyA
            E + P ⇌ EP        :: OnlyA
            E + I ⇌ EI        :: EqualAI
            ES + I ⇌ EIS  :: EqualAI
            EI + S ⇌ EIS  :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(deadend)
    @test v == EnzymeRates.Derivation._Verdict(FP, U, 0, 2,
                                               _testhelper_edge(deadend, "E + S ⇌ ES"))
    @test _testhelper_epath_verdict(deadend; rng) === :no_path
    @test EnzymeRates.Derivation._guard_message(deadend) ==
          _testhelper_pinned_cause(deadend, v.edge, 2)
    @test EnzymeRates.Derivation._guard(deadend_untagged) ==
          EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(deadend_joins_tagged).fault === FP
    # I's binding at E :OnlyA too: untagging it leaves the S binding inside one component,
    # while tagging the steps that join its sides separates every binding, so the message
    # offers only the joins fix.
    two_bindings = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
            E + I ⇌ EI        :: OnlyA
            ES + I ⇌ EIS  :: EqualAI
            EI + S ⇌ EIS  :: EqualAI
        end
    end
    two_bindings_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
            E + I ⇌ EI        :: EqualAI
            ES + I ⇌ EIS  :: EqualAI
            EI + S ⇌ EIS  :: EqualAI
        end
    end
    two_bindings_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: OnlyA
            E + P ⇌ EP        :: OnlyA
            E + I ⇌ EI        :: OnlyA
            ES + I ⇌ EIS  :: OnlyA
            EI + S ⇌ EIS  :: OnlyA
        end
    end
    v = EnzymeRates.Derivation._guard(two_bindings)
    @test v == EnzymeRates.Derivation._Verdict(FP, U, 0, 2,
                                               _testhelper_edge(two_bindings, "E + I ⇌ EI"))
    @test _testhelper_epath_verdict(two_bindings; rng) === :no_path
    joins = join((EnzymeRates.Mechanisms.step_text(two_bindings.catalytic,
                                                   _testhelper_edge(two_bindings, t))
                  for t in ("E + P ⇌ EP", "ES <--> EP", "ES + I ⇌ EIS", "EI + S ⇌ EIS")),
                 "; ")
    @test EnzymeRates.Derivation._guard_message(two_bindings) ==
          _testhelper_pinned_cause(two_bindings, v.edge, 2) *
          " Tag :OnlyA the steps that join its sides: $joins."
    @test EnzymeRates.Derivation._guard(two_bindings_untagged).fault === FP
    @test EnzymeRates.Derivation._guard(two_bindings_joins_tagged) ==
          EnzymeRates.Derivation._WELL_POSED
end

@testset "a tagged binding pinned inside one group, or closing a cycle of a reading" begin
    rng = Random.MersenneTwister(37)
    FP, U = EnzymeRates.Derivation._FAST_STEP_PINNED, EnzymeRates.Derivation._LIMIT_UNTESTED
    # The ter cube of seven :OnlyA edges: the untagged EA + C, EC + B and EB + A
    # bindings join EA to EAB, so its rapid-equilibrium :OnlyA binding of B stays
    # inside one component and cannot be switched off.
    cube = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, C
        products: P
        steps: begin
            E + A ⇌ EA             :: OnlyA
            E + B ⇌ EB             :: OnlyA
            E + C ⇌ EC             :: OnlyA
            EA + B ⇌ EAB       :: OnlyA
            EA + C ⇌ EAC       :: EqualAI
            EB + A ⇌ EAB       :: EqualAI
            EB + C ⇌ EBC       :: EqualAI
            EC + A ⇌ EAC       :: EqualAI
            EC + B ⇌ EBC       :: EqualAI
            EAB + C ⇌ EABC :: OnlyA
            EAC + B ⇌ EABC :: OnlyA
            EBC + A ⇌ EABC :: OnlyA
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
        end
    end
    @test _testhelper_exclusion(cube, 2).failure === :fast
    @test EnzymeRates.Derivation._guard(cube) ==
          EnzymeRates.Derivation._Verdict(FP, U, 0, 2,
                                          _testhelper_edge(cube, "EA + B ⇌ EAB"))
    @test _testhelper_epath_verdict(cube; rng) === :no_path
    # The message offers three fixes, each passing the guard: the B binding at EA
    # untagged, the five untagged bindings that join its sides (EA + C, EC + A,
    # EC + B, EB + C, EB + A) tagged, and the B binding steady-state.
    cube_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, C
        products: P
        steps: begin
            E + A ⇌ EA             :: OnlyA
            E + B ⇌ EB             :: OnlyA
            E + C ⇌ EC             :: OnlyA
            EA + B ⇌ EAB       :: EqualAI
            EA + C ⇌ EAC       :: EqualAI
            EB + A ⇌ EAB       :: EqualAI
            EB + C ⇌ EBC       :: EqualAI
            EC + A ⇌ EAC       :: EqualAI
            EC + B ⇌ EBC       :: EqualAI
            EAB + C ⇌ EABC :: OnlyA
            EAC + B ⇌ EABC :: OnlyA
            EBC + A ⇌ EABC :: OnlyA
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
        end
    end
    cube_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, C
        products: P
        steps: begin
            E + A ⇌ EA             :: OnlyA
            E + B ⇌ EB             :: OnlyA
            E + C ⇌ EC             :: OnlyA
            EA + B ⇌ EAB       :: OnlyA
            EA + C ⇌ EAC       :: OnlyA
            EB + A ⇌ EAB       :: OnlyA
            EB + C ⇌ EBC       :: OnlyA
            EC + A ⇌ EAC       :: OnlyA
            EC + B ⇌ EBC       :: OnlyA
            EAB + C ⇌ EABC :: OnlyA
            EAC + B ⇌ EABC :: OnlyA
            EBC + A ⇌ EABC :: OnlyA
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(cube_untagged) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(cube_joins_tagged) ==
          EnzymeRates.Derivation._WELL_POSED
    joins = join((EnzymeRates.Mechanisms.step_text(cube.catalytic,
                                                   _testhelper_edge(cube, t))
                  for t in ("EA + C ⇌ EAC", "EC + A ⇌ EAC", "EC + B ⇌ EBC",
                            "EB + C ⇌ EBC", "EB + A ⇌ EAB")), "; ")
    @test endswith(
        EnzymeRates.Derivation._guard_message(cube),
        "Make its kinetic group steady-state, or tag :OnlyA the steps that join its " *
        "sides: $joins.")
    # Both bindings of the ordered bi-uni :OnlyA: the chemistry and the P binding join
    # EAB back to E, so EA alone is cut off. No built reading switches off the B
    # binding, which takes up B into a form free E's group holds; the one consistent
    # reading loses the A binding's uptake and the B binding's release, both into EA.
    # Their vanishing directions then bridge E and EAB, and the limit depends on the
    # direction.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ EA            :: OnlyA
            EA + B ⇌ EAB      :: OnlyA
            EAB <--> EP       :: EqualAI
            E + P ⇌ EP            :: EqualAI
        end
    end
    b = _testhelper_edge(am, "EA + B ⇌ EAB")
    @test _testhelper_exclusion(am, 2).failure === :fast
    excluded(c) = EnzymeRates.Derivation._excluded_steps(am, c)
    @test EnzymeRates.Derivation.exclusion_components(am.catalytic, excluded(2)) ==
          [1, 2, 1, 1]  # E, EA, EAB, EP
    @test EnzymeRates.Derivation.exclusion_components(am.catalytic, excluded(1)) ==
          [1, 1, 1, 1]
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._BRIDGE, U, 1, 2, b)
    @test _testhelper_epath_verdict(am; rng) === :no_path
    reported = _testhelper_reading(EnzymeRates.Derivation._guard_with_reading(am).reading)
    @test _testhelper_reading_is_consistent(am, reported)
    @test _testhelper_epath_verdict(am; rng, corners = true, reading = reported) ===
          :direction_dependent
    @test EnzymeRates.Derivation._guard_message(am) ==
          "Mechanism: the :OnlyA step EA + B ⇌ EAB bridges forms of the " *
          "inactive conformation that its other steps keep apart, so admission cannot " *
          "show that the limit does not depend on how fast that step vanishes. Make the " *
          "kinetic group of EA + B ⇌ EAB steady-state."
    # Both bindings untagged, or the chemistry and the P binding tagged too, the mechanism
    # passes the guard; untagging the B binding alone leaves the A binding's sides joined.
    untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ EA            :: EqualAI
            EA + B ⇌ EAB      :: EqualAI
            EAB <--> EP       :: EqualAI
            E + P ⇌ EP            :: EqualAI
        end
    end
    joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ EA            :: OnlyA
            EA + B ⇌ EAB      :: OnlyA
            EAB <--> EP       :: OnlyA
            E + P ⇌ EP            :: OnlyA
        end
    end
    one_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ EA            :: OnlyA
            EA + B ⇌ EAB      :: EqualAI
            EAB <--> EP       :: EqualAI
            E + P ⇌ EP            :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(untagged) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(joins_tagged) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(one_untagged).fault === FP
    # With I forms, EA is still cut off alone, and the one consistent reading loses both
    # bindings into it, as above. Admission passes the faces with I present, and class 3
    # rejects zero I, where the law splits a rapid-equilibrium segment, as the RE-speed
    # oracle finds. Tagging the chemistry and the P binding too leaves the law a second
    # closed class on a face.
    with_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        modifiers: I
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B ⇌ EAB          :: OnlyA
            EAB <--> EP           :: EqualAI
            E + P ⇌ EP                :: EqualAI
            E + I ⇌ EI                :: EqualAI
            EAB + I ⇌ EABI    :: EqualAI
            EI + A ⇌ EAI          :: EqualAI
            EAI + B ⇌ EABI    :: EqualAI
        end
    end
    with_i_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        modifiers: I
        steps: begin
            E + A ⇌ EA                :: EqualAI
            EA + B ⇌ EAB          :: EqualAI
            EAB <--> EP           :: EqualAI
            E + P ⇌ EP                :: EqualAI
            E + I ⇌ EI                :: EqualAI
            EAB + I ⇌ EABI    :: EqualAI
            EI + A ⇌ EAI          :: EqualAI
            EAI + B ⇌ EABI    :: EqualAI
        end
    end
    with_i_joins_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        modifiers: I
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B ⇌ EAB          :: OnlyA
            EAB <--> EP           :: OnlyA
            E + P ⇌ EP                :: OnlyA
            E + I ⇌ EI                :: EqualAI
            EAB + I ⇌ EABI    :: EqualAI
            EI + A ⇌ EAI          :: EqualAI
            EAI + B ⇌ EABI    :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(with_i)
    @test (v.fault, _testhelper_faces(with_i)[v.face]) ==
          (EnzymeRates.Derivation._SPLIT, Set([:I]))
    @test _testhelper_epath_verdict(with_i; rng) === :no_path
    @test _testhelper_agrees(v.fault, _testhelper_respeed_verdicts(
        with_i; rng = Random.MersenneTwister(38), law = :limit)[v.face])
    @test EnzymeRates.Derivation._guard_message(with_i) ==
          "Mechanism: at zero I, the forms {E (active), EA (active), " *
          "EAB (active), EP (active), E (inactive), EP (inactive)} and " *
          "{EAB (inactive)} of one rapid-equilibrium segment are joined only by " *
          "steps that take up a zeroed metabolite (E + I ⇌ EI and " *
          "EAB + I ⇌ EABI), so the rate law cannot weight one against the other."
    @test EnzymeRates.Derivation._guard(with_i_untagged) ==
          EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(with_i_joins_tagged).fault ===
          EnzymeRates.Derivation._NONUNIQUE
    # Under the guard's built readings the arcs close a cycle between two cut-off parts:
    # EA with EAQ, and EQ with EAB and EPQ. The B binding at EA enters
    # the second, the A binding at EQ the first, so the two affinities cannot vanish
    # together. A reading in which one of them loses its release instead has no cycle,
    # and the search admits the mechanism under it, the oracles agreeing.
    cycle = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA              :: OnlyA
            EA + Q ⇌ EAQ        :: EqualAI
            EA + B ⇌ EAB        :: OnlyA
            EAB <--> EPQ      :: EqualAI
            EQ + P ⇌ EPQ        :: EqualAI
            EQ + A ⇌ EAQ        :: OnlyA
            E + Q ⇌ EQ              :: OnlyA
        end
    end
    line_c(k) = EnzymeRates.Mechanisms.step_text(cycle.catalytic, k)
    b_c, a_c = (_testhelper_edge(cycle, t) for t in ("EA + B ⇌ EAB", "EQ + A ⇌ EAQ"))
    @test _testhelper_exclusion(cycle, 2).failure === :cycle
    @test EnzymeRates.Derivation._admission_interior(
        cycle, EnzymeRates.Derivation._built_reading(cycle)) ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._ARC_CYCLE, U, 0, 2, a_c)
    @test _testhelper_epath_verdict(cycle; rng) === :no_path
    @test EnzymeRates.Derivation._guard(cycle) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(cycle; rng, corners = true))
    _testhelper_check_law_against_limit(cycle)
end

@testset "forms the built readings reach only through tagged chemistry" begin
    rng = Random.MersenneTwister(41)
    D = EnzymeRates.Derivation
    IS, U = EnzymeRates.Derivation._ISLAND, EnzymeRates.Derivation._LIMIT_UNTESTED
    FP, PU = EnzymeRates.Derivation._FAST_STEP_PINNED, EnzymeRates.Derivation._PUMPED
    # Each mechanism here the guard's built readings reject, with an island, a pumped
    # basin or a step they cannot switch off; most the search admits, the oracles agreeing
    # on every face under the reading it admits with.
    function _testhelper_admitted_by_search(x)
        @test D._guard(x) == D._WELL_POSED
        @test !(D._guard_with_reading(x).reading in
                [D._built_reading(x, r) for r in instances(D._SlowUptakeReading)])
        @test all(_testhelper_admission_agreement(x; rng, corners = true))
        _testhelper_check_law_against_limit(x)
    end
    _testhelper_edge_between(x, a, b) = only(k for (k, e) in enumerate(x.catalytic.steps)
        if Set(String.(x.catalytic.enzyme_forms[[e.from, e.to]])) == Set((a, b)))
    _testhelper_island_verdict(x) =
        _testhelper_epath_verdict(x; rng) in (:direction_dependent, :wrong_limit)
    # Both half-reactions of the ping-pong :OnlyA. Under the first built reading half 2
    # takes up Q into EB_res_plusA_minusP, which the inactive conformation cuts off, so
    # it loses its affinity and free E reaches the covalent forms along it: no island. At
    # zero A and B, B's binding feeds E_res_plusA_minusP, whose only way out, half 1,
    # is stopped: a pumped basin, as the ε-path oracle finds under that reading. The
    # search admits the ping-pong under a reading in which half 1 loses its rate into the
    # covalent forms instead.
    both = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: OnlyA
        end
    end
    half1 = _testhelper_edge(both, "E_res_plusA_minusP + P <--> EA")
    feed = _testhelper_edge(both, "E_res_plusA_minusP + B <--> EB_res_plusA_minusP")
    @test both.catalytic.enzyme_forms[[both.catalytic.steps[half1].from,
                                       both.catalytic.steps[half1].to]] ==
          [:E_res_plusA_minusP, :EA]
    built = D._built_reading(both)
    v = D._verdict_under(both, built)
    @test v == EnzymeRates.Derivation._Verdict(PU, U, 3, 2, feed)
    @test isempty(_testhelper_exclusion(both, 2).lost)
    @test _testhelper_admission_agrees(
        PU, _testhelper_epath_verdict(both; rng, face = _testhelper_faces(both)[3]))
    @test D._admission_message(both, v, built) ==
          "Mechanism: at zero A and B, " *
          "$(EnzymeRates.Mechanisms.step_text(both.catalytic, feed)) feeds the forms " *
          "E_res_plusA_minusP of the inactive conformation, which leave only through " *
          "the vanishing step $(EnzymeRates.Mechanisms.step_text(both.catalytic, half1)) " *
          "while nothing runs back the way it came, so admission cannot show that their " *
          "weight does not depend on how fast that step vanishes."
    _testhelper_admitted_by_search(both)
    @test EnzymeRates.Derivation._guard(_testhelper_tag_sets["pp_V_half2"]) ==
          EnzymeRates.Derivation._WELL_POSED
    # R3 and R4: a second free form behind a tagged isomerization; untagging it passes.
    # Under the built readings R3's isomerization is steady-state, so it is stopped and
    # F, FI are an island; R4's is at rapid equilibrium, which they cannot switch off.
    # Each isomerization joins F's group to free E's, and the search lets it lose its rate
    # into F.
    r3 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
            E <--> F            :: OnlyA
            F + I ⇌ FI        :: EqualAI
        end
    end
    r3_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
            E <--> F            :: EqualAI
            F + I ⇌ FI        :: EqualAI
        end
    end
    r4 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: NonequalAI
            E + P ⇌ EP        :: EqualAI
            F ⇌ E               :: OnlyA
        end
    end
    r4_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: NonequalAI
            E + P ⇌ EP        :: EqualAI
            F ⇌ E               :: EqualAI
        end
    end
    iso = _testhelper_edge_between(r3, "E", "F")
    # An isomerization of free forms ties on both sides; it is stored from E.
    @test r3.catalytic.steps[iso].from == 1
    @test D._admission_interior(r3, D._built_reading(r3)) ==
          EnzymeRates.Derivation._Verdict(IS, U, 0, 2, iso)
    @test _testhelper_island_verdict(r3)
    _testhelper_admitted_by_search(r3)
    iso = _testhelper_edge_between(r4, "E", "F")
    @test r4.catalytic.steps[iso].from == 1
    @test D._admission_interior(r4, D._built_reading(r4)) ==
          EnzymeRates.Derivation._Verdict(FP, U, 0, 2, iso)
    @test _testhelper_epath_verdict(r4; rng) === :no_path
    _testhelper_admitted_by_search(r4)
    for fixed in (r3_untagged, r4_untagged)
        @test EnzymeRates.Derivation._guard(fixed) == EnzymeRates.Derivation._WELL_POSED
    end
    # AD2b's P binding at F is at rapid equilibrium and joins F to EP, in free E's group,
    # so no built reading switches it off. Its consistent readings let it lose its release
    # into F; with P zero its uptake then dies, and it runs one way into F, which empties
    # only through the slow F <--> E: F fills faster than it empties, as the ε-path oracle
    # finds, and the message offers the steady-state fix, which the search admits. In the
    # ping-pong with B's binding and half 1 :OnlyA, the covalent form is an island behind
    # half 1 under the built readings; untagged, half 1 joins it to E, and B's steady-state
    # binding, its sides joined, is stopped, a fix that passes every face, the oracles
    # agreeing.
    ad2b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: EqualAI
            F + P ⇌ EP        :: OnlyA
            F <--> E            :: OnlyA
        end
    end
    bbind = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: OnlyA
            EB_res_plusA_minusP <--> E + Q                     :: EqualAI
        end
    end
    bbind_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: EqualAI
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: OnlyA
            EB_res_plusA_minusP <--> E + Q                     :: EqualAI
        end
    end
    ad2b_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> EP      :: EqualAI
            F + P <--> EP     :: OnlyA
            F <--> E            :: OnlyA
        end
    end
    pbind = _testhelper_edge(ad2b, "F + P ⇌ EP")
    @test D._admission_interior(ad2b, D._built_reading(ad2b)) ==
          EnzymeRates.Derivation._Verdict(FP, U, 0, 2, pbind)
    @test _testhelper_epath_verdict(ad2b; rng) === :no_path
    v = EnzymeRates.Derivation._guard(ad2b)
    @test v == EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._FILLED, U, 2, 2,
                                               pbind)
    reported = _testhelper_reading(D._guard_with_reading(ad2b).reading)
    @test _testhelper_reading_is_consistent(ad2b, reported)
    @test _testhelper_admission_agrees(
        v.fault, _testhelper_epath_verdict(ad2b; rng, face = _testhelper_faces(ad2b)[2],
                                           corners = true, reading = reported))
    @test EnzymeRates.Derivation._guard_message(ad2b) ==
          "Mechanism: at zero P, in the inactive conformation the :OnlyA " *
          "rapid-equilibrium step F + P ⇌ EP runs one way into F, at a rate that " *
          "vanishes yet grows with the rapid-equilibrium speed, and F leaves only " *
          "through the slower step F <--> E, so F fills faster than it empties and its " *
          "weight depends on how fast the tagged step vanishes. Make the kinetic group " *
          "of F + P ⇌ EP steady-state."
    _testhelper_admitted_by_search(ad2b_ss)
    half1 = _testhelper_edge(bbind, "E_res_plusA_minusP + P <--> EA")
    @test D._admission_interior(bbind, D._built_reading(bbind)) ==
          EnzymeRates.Derivation._Verdict(IS, U, 0, 2, half1)
    @test _testhelper_island_verdict(bbind)
    _testhelper_admitted_by_search(bbind)
    @test EnzymeRates.Derivation._guard(bbind_untagged) ==
          EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(bbind_untagged; rng))
    # Half 1 and its copy at EAI share one :OnlyA group, and I's binding at E is :OnlyA
    # too. As in the ping-pong above, under the first built reading, at zero A and B, the
    # covalent forms are fed, by B's binding and by the copy of half 1 at EAI running
    # forward, and their only ways out, the two copies of half 1 taking P back up, vanish.
    # The guard names the feed whose reverse is dead, B's binding: nothing runs back
    # through it. The search admits it as it admits that ping-pong. Untagging that group
    # joins EI to E, so I's rapid-equilibrium binding cannot be switched off.
    shared_half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A <--> EA                                       :: EqualAI
            (EA <--> E_res_plusA_minusP + P,
             EAI <--> EI_res_plusA_minusP + P)             :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: OnlyA
            E + I ⇌ EI                                          :: OnlyA
            EI + A <--> EAI                                 :: EqualAI
            E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP    :: EqualAI
        end
    end
    shared_half1_untagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A <--> EA                                       :: EqualAI
            (EA <--> E_res_plusA_minusP + P,
             EAI <--> EI_res_plusA_minusP + P)             :: EqualAI
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: OnlyA
            E + I ⇌ EI                                          :: OnlyA
            EI + A <--> EAI                                 :: EqualAI
            E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP    :: EqualAI
        end
    end
    shared_line(t) = EnzymeRates.Mechanisms.step_text(shared_half1.catalytic,
                                                      _testhelper_edge(shared_half1, t))
    built = D._built_reading(shared_half1)
    v = D._verdict_under(shared_half1, built)
    @test v == EnzymeRates.Derivation._Verdict(
        PU, U, 3, 2,
        _testhelper_edge(shared_half1, "E_res_plusA_minusP + B <--> EB_res_plusA_minusP"))
    @test _testhelper_admission_agrees(
        PU, _testhelper_epath_verdict(shared_half1; rng,
                                      face = _testhelper_faces(shared_half1)[3]))
    @test D._admission_message(shared_half1, v, built) ==
          "Mechanism: at zero A and B, " *
          "$(shared_line("E_res_plusA_minusP + B <--> EB_res_plusA_minusP")) feeds the " *
          "forms EI_res_plusA_minusP, E_res_plusA_minusP of the inactive " *
          "conformation, which leave only through the vanishing steps " *
          "$(shared_line("EI_res_plusA_minusP + P <--> EAI")); " *
          "$(shared_line("E_res_plusA_minusP + P <--> EA")) while nothing runs back " *
          "the way it came, so admission cannot show that their weight does not depend " *
          "on how fast those steps vanish."
    _testhelper_admitted_by_search(shared_half1)
    @test EnzymeRates.Derivation._guard(shared_half1_untagged).fault === FP
end

"""The tag sets of `_testhelper_tag_sets` whose inactive conformation excludes a
rapid-equilibrium step that the guard's first built reading cannot switch off, so
admission's interior rejects that reading (`_FAST_STEP_PINNED`). In AD7b, AD9 and TC1 the
step lies inside one group, so no reading switches it off; BS2's isomerization joins two
groups, and the search admits BS2."""
const _testhelper_pinned_tag_sets = Set(["AD7b", "AD9", "BS2", "TC1"])

@testset "admission's interior rejects only the tag sets with a pinned fast step" begin
    rng = Random.MersenneTwister(2613)
    @testset "$name" for name in sort!(collect(keys(_testhelper_tag_sets)))
        am = _testhelper_tag_sets[name]
        v = EnzymeRates.Derivation._admission_interior(
            am, EnzymeRates.Derivation._built_reading(am))
        if name in _testhelper_pinned_tag_sets
            @test (v.fault, v.conformation) == (EnzymeRates.Derivation._FAST_STEP_PINNED, 2)
            @test _testhelper_epath_verdict(am; rng) === :no_path
        else
            @test v == EnzymeRates.Derivation._WELL_POSED
        end
    end
end

@testset "a mixed :OnlyA group whose limit exists is admitted" begin
    # B's fused binding into EPQ and its plain binding at EQ share one :OnlyA group;
    # in the inactive conformation no untagged step joins EQ to EBQ, so the plain
    # binding's affinity can vanish, and the limit is deletion plus pruning.
    am = _testhelper_tag_sets["mixed_group_limit"]
    @test _testhelper_exclusion(am, 2).failure === :none
    @test _testhelper_lp_existence(am, 2)
    @test EnzymeRates.Derivation._admission_interior(
        am, EnzymeRates.Derivation._built_reading(am)) ==
          EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(am) == EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(43)) === :converges
end

"""The tag sets of `_testhelper_tag_sets` that a face rejects under the guard's first built
reading; every other one passes every face under it. The search admits R2, BS1, BS2, BS3,
BS4, R2funnel, R2funnel_SS and fz101_10232 under other readings."""
const _testhelper_face_rejected = Set(["R2", "BS1", "BS2", "BS3", "BS4", "ADbridge",
    "ADbridge_noturnover", "bridge_drain", "bridge_deadend", "R2funnel", "R2funnel_SS",
    "S1", "AD7b", "fz101_10232"])

@testset "R2's limit at zero product depends on how fast its isomerization vanishes" begin
    rng = Random.MersenneTwister(43)
    NO, PU = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._PUMPED
    # Classes 1-3 and admission's interior accept R2, so under the guard's first built
    # reading the face condition alone rejects it: at P zero, F's exits vanish (F <--> E,
    # stopped) or take up P (F + P ⇌ EP), and EP feeds it by a release whose reverse is
    # dead. The search admits R2 under a reading in which F <--> E loses its rate into F
    # instead, so F drains into free E, and the oracles agree on every face under it.
    am = _testhelper_tag_sets["R2"]
    m = am.catalytic
    @test _testhelper_face_faults(am) == fill(NO, 3)
    @test EnzymeRates.Derivation._admission_interior(
        am, EnzymeRates.Derivation._built_reading(am)) ==
          EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_admission_faults(am) == [NO, PU, NO]
    @test _testhelper_epath_verdict(am; rng, face = _testhelper_faces(am)[2]) ===
          :direction_dependent
    k, iso = _testhelper_edge(am, "F + P ⇌ EP"), _testhelper_edge(am, "F <--> E")
    built = EnzymeRates.Derivation._built_reading(am)
    v = EnzymeRates.Derivation._verdict_under(am, built)
    @test v == EnzymeRates.Derivation._Verdict(
        PU, EnzymeRates.Derivation._LIMIT_UNTESTED, 2, 2, k)
    line(j) = EnzymeRates.Mechanisms.step_text(m, j)
    @test EnzymeRates.Derivation._admission_message(am, v, built) ==
          "Mechanism: at zero P, $(line(k)) feeds the forms F of the inactive " *
          "conformation, which leave only through the vanishing step $(line(iso)) while " *
          "nothing runs back the way it came, so admission cannot show that their weight " *
          "does not depend on how fast that step vanishes."
    @test EnzymeRates.Derivation._guard(am) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard_with_reading(am).reading[iso, 2] !==
          EnzymeRates.Derivation._STOPPED
    @test all(_testhelper_admission_agreement(am; rng, corners = true))
    _testhelper_check_law_against_limit(am)
end

@testset "a stranded basin with two vanishing exits names both steps in its message" begin
    # Under the guard's first built reading, at P zero, ES feeds F by a release whose
    # reverse takes up P; F and FI leave the inactive conformation's basin only through
    # the :OnlyA steps EI ⇌ FI and F <--> E, which vanish there, so admission's face
    # condition rejects that face. But that reading cannot switch off EI ⇌ FI, at
    # rapid equilibrium, which consistent readings let lose its rate into FI; under the
    # first of them class 3 rejects zero I, where the law splits a rapid-equilibrium
    # segment, as the RE-speed oracle finds.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + I ⇌ EI        :: NonequalAI
            E + S ⇌ ES        :: OnlyA
            EI ⇌ FI         :: OnlyA
            F <--> E            :: OnlyA
            F + P <--> ES     :: NonequalAI
            F + I ⇌ FI        :: OnlyA
        end
    end
    rng = Random.MersenneTwister(43)
    k = _testhelper_edge(am, "FI ⇌ EI")
    @test EnzymeRates.Derivation._admission_interior(
        am, EnzymeRates.Derivation._built_reading(am)) ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._FAST_STEP_PINNED,
                                          EnzymeRates.Derivation._LIMIT_UNTESTED, 0, 2, k)
    @test _testhelper_epath_verdict(am; rng, face = _testhelper_faces(am)[2]) === :no_path
    v = EnzymeRates.Derivation._guard(am)
    @test (v.fault, _testhelper_faces(am)[v.face], v.conformation) ==
          (EnzymeRates.Derivation._SPLIT, Set([:I]), 0)
    @test _testhelper_agrees(v.fault, _testhelper_respeed_verdicts(
        am; rng = Random.MersenneTwister(44), law = :limit)[v.face])
    @test EnzymeRates.Derivation._guard_message(am) ==
          "Mechanism: at zero I, the forms {E (active), ES (active), " *
          "E (inactive)} and {F (active)} of one rapid-equilibrium segment are joined " *
          "only by steps that take up a zeroed metabolite (E + I ⇌ EI and " *
          "F + I ⇌ FI), so the rate law cannot weight one against the other."
    faces = EnzymeRates.Derivation._faces(am.catalytic)
    built = EnzymeRates.Derivation._built_reading(am)
    eps = EnzymeRates.Derivation._eps_network(
        am, EnzymeRates.Derivation._re_flags(am.catalytic), built)
    w = EnzymeRates.Derivation._admission_face(am, eps,
                                               EnzymeRates.Derivation._law_network(am),
                                               faces, 2)
    @test (w.fault, w.face, w.conformation) == (EnzymeRates.Derivation._PUMPED, 2, 2)
    @test EnzymeRates.Derivation._admission_message(am, w, built) ==
          "Mechanism: at zero P, F + P <--> ES feeds the forms F, FI of " *
          "the inactive conformation, which leave only through the vanishing steps " *
          "F <--> E; FI ⇌ EI while nothing runs back the way it came, so admission " *
          "cannot show that their weight does not depend on how fast those steps vanish."
    # With EI <--> FI steady-state the first built reading stops it, the interior
    # passes, and the face condition rejects the pumped basin at zero P, the ε-path
    # oracle agreeing. The search admits it under a reading in which the :OnlyA steps
    # lose their rates into F's group instead, both oracles agreeing on every face.
    stranded_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + I ⇌ EI        :: NonequalAI
            E + S ⇌ ES        :: OnlyA
            EI <--> FI      :: OnlyA
            F <--> E            :: OnlyA
            F + P <--> ES     :: NonequalAI
            F + I ⇌ FI        :: OnlyA
        end
    end
    built_ss = EnzymeRates.Derivation._built_reading(stranded_ss)
    @test EnzymeRates.Derivation._admission_interior(stranded_ss, built_ss) ==
          EnzymeRates.Derivation._WELL_POSED
    v = EnzymeRates.Derivation._verdict_under(stranded_ss, built_ss)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._PUMPED, 2, 2)
    @test _testhelper_admission_agrees(
        v.fault, _testhelper_epath_verdict(stranded_ss; rng,
                                           face = _testhelper_faces(stranded_ss)[2]))
    @test EnzymeRates.Derivation._admission_message(stranded_ss, v, built_ss) ==
          "Mechanism: at zero P, F + P <--> ES feeds the forms F, FI of " *
          "the inactive conformation, which leave only through the vanishing steps " *
          "F <--> E; FI <--> EI while nothing runs back the way it came, so " *
          "admission cannot show that their weight does not depend on how fast those " *
          "steps vanish."
    @test EnzymeRates.Derivation._guard(stranded_ss) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(stranded_ss; rng, corners = true))
    _testhelper_check_law_against_limit(stranded_ss)
end

@testset "admitted tag sets converge on every face" begin
    rng = Random.MersenneTwister(47)
    names = sort!([n for n in keys(_testhelper_tag_sets)
                   if n ∉ _testhelper_face_rejected && n ∉ _testhelper_pinned_tag_sets])
    @test length(names) == 37
    @testset "$name" for name in names
        am = _testhelper_tag_sets[name]
        @test EnzymeRates.Derivation._guard(am) == EnzymeRates.Derivation._WELL_POSED
        @test all(==(EnzymeRates.Derivation._NO_FAULT), _testhelper_admission_faults(am))
        @test all(_testhelper_admission_agreement(am; rng))
    end
    # The search admits these of the tag sets the first built reading rejects.
    @test Set(n for n in union(_testhelper_face_rejected, _testhelper_pinned_tag_sets)
              if EnzymeRates.Derivation._guard(_testhelper_tag_sets[n]) ==
                 EnzymeRates.Derivation._WELL_POSED) ==
          Set(["R2", "BS1", "BS2", "BS3", "BS4", "R2funnel", "R2funnel_SS", "fz101_10232"])
    # AD9's chemistry and TC1's exchange of B for P, both at rapid equilibrium and :OnlyA,
    # join forms the inactive conformation keeps joined, so it cannot switch them off.
    @testset "$name" for name in ("AD9", "TC1")
        am = _testhelper_tag_sets[name]
        @test EnzymeRates.Derivation._guard(am).fault ===
              EnzymeRates.Derivation._FAST_STEP_PINNED
        @test _testhelper_epath_verdict(am; rng) === :no_path
    end
end

@testset "a sealed basin drains, and a pumped one does not" begin
    rng = Random.MersenneTwister(53)
    NO, PU = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._PUMPED
    function _testhelper_stranded(am, i)
        face = EnzymeRates.Derivation._faces(am.catalytic)[i]
        s = EnzymeRates.Derivation._face_sets(
            EnzymeRates.Derivation._eps_network(
                am, EnzymeRates.Derivation._re_flags(am.catalytic),
                EnzymeRates.Derivation._built_reading(am)),
            EnzymeRates.Derivation._law_network(am).present, face)
        any(s.stranded)
    end
    # A stranded basin entered only reversibly from the law network is sealed: its weight
    # vanishes with the tagged constants, and the face passes. R2P_sealed at P zero,
    # B_uphill and BS4 with X or I zero.
    for (name, sealed) in (("R2P_sealed", [2]), ("B_uphill", [4, 5, 6]), ("BS4", [4, 5, 6]))
        am = _testhelper_tag_sets[name]
        @test all(i -> _testhelper_stranded(am, i), sealed)
        @test all(==(NO), _testhelper_admission_faults(am)[sealed])
    end
    # Under the guard's first built reading, which stops F <--> E, a basin fed by a step
    # whose reverse is dead is pumped: BS2, BS3 and BS4 at P zero. The search admits BS3
    # and BS4 under a reading in which F <--> E loses its rate into F, so F drains into
    # free E; the oracles agree on every face under it.
    @test _testhelper_admission_faults(_testhelper_tag_sets["BS2"]) == [NO, PU, NO]
    @test _testhelper_admission_faults(_testhelper_tag_sets["BS3"]) == [NO, PU, NO]
    @test _testhelper_admission_faults(_testhelper_tag_sets["BS4"]) ==
          [NO, PU, NO, NO, NO, NO]
    for name in ("R2P_sealed", "B_uphill")
        @test all(_testhelper_admission_agreement(_testhelper_tag_sets[name]; rng))
    end
    for name in ("BS3", "BS4")
        am = _testhelper_tag_sets[name]
        @test (name, EnzymeRates.Derivation._guard(am)) ==
              (name, EnzymeRates.Derivation._WELL_POSED)
        @test all(_testhelper_admission_agreement(am; rng, corners = true))
        _testhelper_check_law_against_limit(am)
    end
    # BS2's isomerization F ⇌ E is at rapid equilibrium, and the built readings cannot
    # switch it off. It joins F's group to free E's, and the search lets it lose its rate
    # into F.
    bs2 = _testhelper_tag_sets["BS2"]
    @test EnzymeRates.Derivation._admission_interior(
        bs2, EnzymeRates.Derivation._built_reading(bs2)).fault ===
          EnzymeRates.Derivation._FAST_STEP_PINNED
    @test _testhelper_epath_verdict(bs2; rng) === :no_path
    @test EnzymeRates.Derivation._guard(bs2) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(bs2; rng, corners = true))
    _testhelper_check_law_against_limit(bs2)
    # With F <--> E steady-state (R2's graph) the first built reading stops it, the
    # interior passes, and the face condition rejects the pumped basin at zero P; the
    # search admits it as it admits R2, both oracles agreeing on every face.
    bs2_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES <--> EP      :: EqualAI
            F + P ⇌ EP        :: EqualAI
            F <--> E            :: OnlyA
        end
    end
    v = EnzymeRates.Derivation._verdict_under(bs2_ss,
                                              EnzymeRates.Derivation._built_reading(bs2_ss))
    @test (v.fault, v.face, v.conformation) == (PU, 2, 2)
    @test _testhelper_admission_faults(bs2_ss) == [NO, PU, NO]
    @test EnzymeRates.Derivation._guard(bs2_ss) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(bs2_ss; rng, corners = true))
    _testhelper_check_law_against_limit(bs2_ss)
end

@testset "a tagged binding that bridges two basins under a reading" begin
    rng = Random.MersenneTwister(59)
    NO, BR = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._BRIDGE
    SP, TR = EnzymeRates.Derivation._SPLIT, EnzymeRates.Derivation._TRAP
    # ADbridge: the :OnlyA bindings of I join E and ES of the inactive conformation, which
    # turns over through its steady-state S binding. With them the fed forms' cheapest
    # in-tree needs one slow step fewer (R1 fails), and the turnover is no potential (T0+
    # fails). With I zero the bridge is dead, and class 3 finds the inhibitor square split.
    am = _testhelper_tag_sets["ADbridge"]
    m = am.catalytic
    @test _testhelper_admission_faults(am) == [BR, BR, BR, SP, SP, TR]
    @test all(_testhelper_admission_agreement(am; rng))
    k = _testhelper_edge(am, "ES + I ⇌ EIS")
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(
        BR, EnzymeRates.Derivation._LIMIT_UNTESTED, 1, 2, k)
    line = EnzymeRates.Mechanisms.step_text(m, k)
    @test EnzymeRates.Derivation._guard_message(am) ==
          "Mechanism: the :OnlyA step $line bridges forms of the inactive " *
          "conformation that its other steps keep apart, so admission cannot show that " *
          "the limit does not depend on how fast that step vanishes. Make the kinetic " *
          "group of $line steady-state."
    fixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> ES                         :: EqualAI
            ES <--> E + P                         :: EqualAI
            (E + I <--> EI, ES + I <--> EIS) :: OnlyA
            EI + S ⇌ EIS                      :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(fixed) == EnzymeRates.Derivation._WELL_POSED
    # The same bridge with the chemistry :OnlyA turns nothing over in the inactive
    # conformation: R1 still fails, but the first substrate's uptake is a potential on
    # every strongly connected part, so T0+ admits the first three faces.
    noturn = _testhelper_tag_sets["ADbridge_noturnover"]
    @test _testhelper_admission_faults(noturn) == [NO, NO, NO, SP, SP, TR]
    @test all(_testhelper_admission_agreement(noturn; rng))
    for (x, potential) in ((am, false), (noturn, true))
        eps = EnzymeRates.Derivation._eps_network(
            x, EnzymeRates.Derivation._re_flags(x.catalytic),
            EnzymeRates.Derivation._built_reading(x))
        law = EnzymeRates.Derivation._law_network(x)
        row = first(x.catalytic.substrates)
        for face in EnzymeRates.Derivation._faces(x.catalytic)[1:3]
            s = EnzymeRates.Derivation._face_sets(eps, law.present, face)
            nodes = s.fed .& .!s.stranded
            (; van, held) = EnzymeRates.Derivation._copy_transitions(eps, s.alive, nodes, 2)
            @test EnzymeRates.Derivation._shortens_in_tree(eps, nodes, s.n0, van)
            @test EnzymeRates.Derivation._potential_on_sccs(eps, held, row) == potential
        end
    end
    # bridge_drain: :OnlyA bindings of I let the inactive F return to E.
    @test _testhelper_admission_faults(_testhelper_tag_sets["bridge_drain"]) ==
          [BR, BR, BR, SP, NO, SP]
    @test all(_testhelper_admission_agreement(_testhelper_tag_sets["bridge_drain"]; rng))
    drain = _testhelper_tag_sets["bridge_drain"]
    k = _testhelper_edge(drain, "F + I ⇌ FI")
    v = EnzymeRates.Derivation._guard(drain)
    @test v == EnzymeRates.Derivation._Verdict(
        BR, EnzymeRates.Derivation._LIMIT_UNTESTED, 1, 2, k)
    @test EnzymeRates.Derivation._guard_message(drain) ==
          "Mechanism: the :OnlyA step F + I ⇌ FI bridges forms of the " *
          "inactive conformation that its other steps keep apart, so admission cannot " *
          "show that the limit does not depend on how fast that step vanishes. Make the " *
          "kinetic group of F + I ⇌ FI steady-state."
    drain_fixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES        :: EqualAI
            ES <--> F + P     :: EqualAI
            F <--> E            :: OnlyA
            E + I ⇌ EI        :: OnlyA
            F + I <--> FI     :: OnlyA
            FI ⇌ EI         :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(drain_fixed) == EnzymeRates.Derivation._WELL_POSED
    # R2funnel: under the guard's first built reading, which stops F <--> E, at P zero the
    # trap {F, FS} drains only through a slow step out of a fast transient, and the
    # bridge stays with the S binding steady-state (R2funnel_SS); the limit there depends
    # on the direction. The search admits both under a reading in which F <--> E loses
    # its rate into F instead, the oracles agreeing on every face.
    for name in ("R2funnel", "R2funnel_SS")
        x = _testhelper_tag_sets[name]
        @test _testhelper_admission_faults(x) == [NO, BR, NO]
        w = EnzymeRates.Derivation._verdict_under(
            x, EnzymeRates.Derivation._built_reading(x))
        @test (w.fault, w.face, w.conformation) === (BR, 2, 2)
        @test _testhelper_epath_verdict(x; rng, face = _testhelper_faces(x)[2],
                                        corners = true) === :direction_dependent
        @test EnzymeRates.Derivation._guard(x) == EnzymeRates.Derivation._WELL_POSED
        @test all(_testhelper_admission_agreement(x; rng, corners = true))
        _testhelper_check_law_against_limit(x)
    end
end

@testset "the census ter-bi mechanism passes by a potential on every component" begin
    rng = Random.MersenneTwister(61)
    am = _testhelper_tag_sets["terbi_d75d2415361c4dae"]
    m = am.catalytic
    @test EnzymeRates.Derivation._guard(am) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(am; rng))
    # With P zero its :OnlyA bindings shorten the inactive conformation's in-tree (R1
    # fails), yet the first substrate's uptake is a potential on each strongly connected
    # part of its kept and fast transitions (T0+): the census's case for T0+.
    eps = EnzymeRates.Derivation._eps_network(am, EnzymeRates.Derivation._re_flags(m),
                                              EnzymeRates.Derivation._built_reading(am))
    law = EnzymeRates.Derivation._law_network(am)
    row = first(m.substrates)
    faces = EnzymeRates.Derivation._faces(m)
    by_potential = map(faces) do face
        s = EnzymeRates.Derivation._face_sets(eps, law.present, face)
        nodes = s.fed .& .!s.stranded
        (; van, held) = EnzymeRates.Derivation._copy_transitions(eps, s.alive, nodes, 2)
        EnzymeRates.Derivation._shortens_in_tree(eps, nodes, s.n0, van) &&
            EnzymeRates.Derivation._potential_on_sccs(eps, held, row)
    end
    @test [Set(EnzymeRates.Derivation._face_names(m, f)) for f in faces[by_potential]] ==
          [Set([:P, :Q]), Set([:P])]
end

@testset "conservative rejections the oracle admits stay rejected under their reading" begin
    rng = Random.MersenneTwister(67)
    NO, BR, PU = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._BRIDGE,
        EnzymeRates.Derivation._PUMPED
    _testhelper_corner_verdict(am, i) =
        _testhelper_epath_verdict(am; rng, face = _testhelper_faces(am)[i], corners = true)
    # R-B1, an equilibrium bridge: the bridge joins E and the dead end EJ, which deletion
    # already holds at its equilibrium ratio, so the limit is deletion's, but no potential
    # proves it. The steady-state fix passes and the message offers it.
    am = _testhelper_tag_sets["bridge_deadend"]
    @test _testhelper_admission_faults(am) == [BR, BR, BR, NO, NO, NO]
    @test all(i -> _testhelper_corner_verdict(am, i) === :converges, 1:3)
    @test all(_testhelper_admission_agreement(am; rng)[4:6])
    k = _testhelper_edge(am, "EJ + I ⇌ EIJ")
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(
        BR, EnzymeRates.Derivation._LIMIT_UNTESTED, 1, 2, k)
    @test endswith(EnzymeRates.Derivation._guard_message(am),
                   "Make the kinetic group of " *
                   "$(EnzymeRates.Mechanisms.step_text(am.catalytic, k)) steady-state.")
    fixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S <--> ES         :: EqualAI
            ES <--> E + P         :: EqualAI
            E + J <--> EJ         :: EqualAI
            E + I ⇌ EI            :: OnlyA
            EJ + I <--> EIJ   :: OnlyA
            EI + J ⇌ EIJ      :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(fixed) == EnzymeRates.Derivation._WELL_POSED
    # R-B2, under the guard's first built reading, which stops E <--> F: a stranded basin
    # fed only through vanishing sources, BS1 with I zero and P or S zero, fz101_10232 on
    # its two faces with S zero. Every corner direction converges. The search admits both
    # under a reading in which E <--> F loses its rate into F, the oracles agreeing on
    # every face.
    bs1 = _testhelper_tag_sets["BS1"]
    @test _testhelper_admission_faults(bs1) == [NO, NO, NO, NO, PU, PU]
    @test all(i -> _testhelper_corner_verdict(bs1, i) === :converges, 5:6)
    # The step that rejects BS1 under that reading has both its forms in the basin: it
    # runs one way among them rather than feeding them.
    built = EnzymeRates.Derivation._built_reading(bs1)
    v = EnzymeRates.Derivation._verdict_under(bs1, built)
    @test (v.fault, v.face, v.conformation) == (PU, 5, 2)
    @test EnzymeRates.Derivation._admission_message(bs1, v, built) ==
          "Mechanism: at zero I and P, F + P ⇌ FP runs one way among the " *
          "forms F, FP, FS of the inactive conformation, which leave only through " *
          "the vanishing step E <--> F while nothing runs back the way it came, so " *
          "admission cannot show that their weight does not depend on how fast that " *
          "step vanishes."
    fz = _testhelper_tag_sets["fz101_10232"]
    @test _testhelper_admission_faults(fz) == [NO, PU, PU, NO, PU, PU]
    @test _testhelper_corner_verdict(fz, 3) === :converges
    # On both faces with P zero, I zero or not (faces 2 and 5), that reading's rejection
    # is right: at face 2 F's only alive exit, the blocked F <--> E, can vanish faster
    # than its inflow, which one corner direction shows.
    @test _testhelper_corner_verdict(fz, 2) !== :converges
    # The faces with I zero as well come last, so the earlier draws of `rng` stay put.
    @test _testhelper_corner_verdict(fz, 5) !== :converges
    @test _testhelper_corner_verdict(fz, 6) === :converges
    for (x, corners) in ((bs1, false), (fz, true))
        @test EnzymeRates.Derivation._guard(x) == EnzymeRates.Derivation._WELL_POSED
        @test all(_testhelper_admission_agreement(x; rng, corners))
        _testhelper_check_law_against_limit(x)
    end
    # R-B3, the set-aside rule's equilibrium condition. In this ping-pong's inactive
    # conformation the cut-off forms EA, FP, F and FB (F is E_res_plusA_minusP)
    # touch free E and EQ, which only the steady-state E + Q <--> EQ joins, so they are
    # not set aside and the interior stays a bridge. Yet the route through them crosses the
    # slow chemistry EA <--> FP, and every corner direction converges there: a
    # conservative rejection.
    D = EnzymeRates.Derivation
    rb3 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA                                            :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP      :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP      :: OnlyA
            EA <--> EP_res_plusA_minusP                        :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                           :: OnlyA
        end
    end
    r = fill(D._RUNS, length(rb3.catalytic.steps), 2)
    r[_testhelper_edge(rb3, "E + A ⇌ EA"), 2] = D._LOSES_FORWARD
    r[_testhelper_edge(rb3, "EA <--> EP_res_plusA_minusP"), 2] = D._LOSES_FORWARD
    r[_testhelper_edge(rb3, "E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP"), 2] =
        D._LOSES_REVERSE
    r[_testhelper_edge(rb3, "EB_res_plusA_minusP ⇌ EQ"), 2] = D._LOSES_REVERSE
    # The bridge is led by the last of the two vanishing steps on a cheapest tree, in
    # step order.
    @test D._verdict_under(rb3, r) ==
          D._Verdict(BR, D._LIMIT_UNTESTED, 1, 2,
                     _testhelper_edge(rb3, "E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP"))
    @test _testhelper_epath_verdict(rb3; rng, corners = true,
                                    reading = _testhelper_reading(r)) === :converges
    # R-B4, the set-aside rule's safety clause. In this random-order mechanism's inactive
    # conformation, with P and Q zero, the stopped E + P <--> EP runs one way from EP
    # into free E along its vanishing direction, a pair at a law form that is not safe;
    # the cut-off forms touch free E alone and their fast classes drain, so that pair
    # alone keeps them, and the face stays a bridge. Every corner direction converges
    # there: a conservative rejection.
    rb4 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP)        :: OnlyA
            (E + P <--> EP, EA + P <--> EAP)  :: OnlyA
            EA + B <--> EPQ                      :: EqualAI
            EP + Q ⇌ EPQ                         :: OnlyA
        end
    end
    r = fill(D._RUNS, length(rb4.catalytic.steps), 2)
    for text in ("E + A ⇌ EA", "EA + P <--> EAP")
        r[_testhelper_edge(rb4, text), 2] = D._LOSES_FORWARD
    end
    r[_testhelper_edge(rb4, "E + P <--> EP"), 2] = D._STOPPED
    for text in ("EP + A ⇌ EAP", "EP + Q ⇌ EPQ")
        r[_testhelper_edge(rb4, text), 2] = D._LOSES_REVERSE
    end
    @test _testhelper_faces(rb4)[2] == Set([:P, :Q])
    let eps = D._eps_network(rb4, D._re_flags(rb4.catalytic), r),
        law = D._law_network(rb4), faces = D._faces(rb4.catalytic)
        @test D._without_set_aside_groups(eps, law.present, faces[2]) === eps
    end
    @test D._verdict_under(rb4, r) ==
          D._Verdict(BR, D._LIMIT_UNTESTED, 2, 2, _testhelper_edge(rb4, "EP + Q ⇌ EPQ"))
    @test _testhelper_epath_verdict(rb4; rng, face = _testhelper_faces(rb4)[2],
                                    corners = true, reading = _testhelper_reading(r)) ===
          :converges
end

@testset "S1 passes admission's interior and fails class 3, as AD7b's faces do" begin
    rng = Random.MersenneTwister(71)
    NO, NU = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._NONUNIQUE
    # S1's stray free form F binds B at steady state: with A or B zero the law network has
    # a second closed class, which class 3 reports before admission reads those faces.
    s1 = _testhelper_tag_sets["S1"]
    @test EnzymeRates.Derivation._admission_interior(
        s1, EnzymeRates.Derivation._built_reading(s1)) ==
          EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_admission_faults(s1) == [NO, NO, NU, NU, NU]
    @test EnzymeRates.Derivation._guard(s1).fault === NU
    @test all(_testhelper_admission_agreement(s1; rng))
    # AD7b: with I zero no spanning in-tree reaches the fed forms, so the face condition
    # passes and class 3 decides. Its isomerization E ⇌ F is at rapid equilibrium, and
    # the inactive conformation cannot switch it off, so admission rejects the interior.
    ad7b = _testhelper_tag_sets["AD7b"]
    @test EnzymeRates.Derivation._verdict_under(
        ad7b, EnzymeRates.Derivation._built_reading(ad7b)).fault ===
          EnzymeRates.Derivation._FAST_STEP_PINNED
    @test _testhelper_admission_faults(ad7b) == [NO, NO, NO, NU, NU, NU]
    @test _testhelper_epath_verdict(ad7b; rng) === :no_path
    # Class 3's faces agree with the RE-speed oracle.
    @test all(_testhelper_admission_agreement(ad7b; rng)[4:6])
    # With E <--> F steady-state it is stopped: admission's interior passes, and the guard
    # reaches class 3's second closed class with I zero; both oracles agree on every face.
    ad7b_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I
        steps: begin
            E + A <--> EA                 :: EqualAI
            (E + I ⇌ EI, F + I ⇌ EI)    :: EqualAI
            E + P <--> EP                 :: EqualAI
            E <--> F                        :: OnlyA
            EA <--> EP                  :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._admission_interior(
        ad7b_ss, EnzymeRates.Derivation._built_reading(ad7b_ss)) ==
          EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_admission_faults(ad7b_ss) == [NO, NO, NO, NU, NU, NU]
    v = EnzymeRates.Derivation._guard(ad7b_ss)
    @test (v.fault, v.face, v.conformation) == (NU, 4, 2)
    @test all(_testhelper_admission_agreement(ad7b_ss; rng))
end

@testset "joint R1 catches a tree that only both conformations' bridges shorten" begin
    # The in-tree costs where each conformation passes R1 and both together fail it,
    # written in the guard's own type: the smallest network that shows the pattern, for
    # the in-tree arithmetic; "joint R1 rejects a mechanism whose bridges cut the tree
    # only together" tests the check on a mechanism.
    # Copy A holds E and X (nodes 1, 2), copy I holds E, Y and Z (nodes 4, 5, 6), and the
    # flip joins the two copies of E. Alone, copy A's vanishing E → X and copy I's
    # vanishing Y → E each leave the cheapest spanning in-tree at one slow transition;
    # together they make it free, rooted at X. Transitions in order: E→X slow, X→E slow,
    # E→X fast and vanishing (copy A); E→Y fast, Y→Z slow, Z→Y fast, Z→E slow, Y→E fast and
    # vanishing (copy I); the flip. Only the transitions these costs need are present, so
    # most lack the reverse partner a mechanism's would have; the face condition reads the
    # partners only on a stranded basin, which this network has none of.
    net = EnzymeRates.Derivation._Network(3, 2, BitVector([1, 1, 0, 1, 1, 1]),
        [1, 2, 1, 4, 5, 6, 6, 5, 1, 4],                  # from
        [2, 1, 2, 5, 6, 5, 4, 4, 4, 1],                  # to
        [1, 1, 2, 3, 4, 5, 6, 7, 0, 0],                  # edge
        [1, 1, 1, 2, 2, 2, 2, 2, 0, 0],                  # conformation
        BitVector([1, 0, 1, 1, 1, 1, 1, 1, 1, 0]),       # forward
        BitVector([0, 0, 1, 1, 0, 1, 0, 1, 1, 1]),       # fast
        BitVector([0, 0, 1, 0, 0, 0, 0, 1, 0, 0]),       # vanishing
        zeros(Int, 1, 10), zeros(UInt64, 10))
    s = EnzymeRates.Derivation._face_sets(net, net.present, UInt64(0))
    nodes = s.fed .& .!s.stranded
    @test nodes == net.present
    @test EnzymeRates.Derivation._fewest_slow(net, nodes, s.n0) == 1
    for c in 1:2
        (; van) = EnzymeRates.Derivation._copy_transitions(net, s.alive, nodes, c)
        @test count(van) == 1
        @test !EnzymeRates.Derivation._shortens_in_tree(net, nodes, s.n0, van)
    end
    @test EnzymeRates.Derivation._fewest_slow(net, nodes, s.alive) == 0
    @test EnzymeRates.Derivation._face_condition(net, net.present, UInt64(0), 1) ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._BRIDGE,
                                          EnzymeRates.Derivation._LIMIT_UNTESTED, 0, 0, 0)
    # The in-tree test reports -1 when no in-tree spans the nodes: nodes 2 and 6 share no
    # transition.
    @test EnzymeRates.Derivation._fewest_slow(
        net, BitVector([0, 1, 0, 0, 0, 1]), trues(10)) == -1
    # T0+'s potential: the copy-I cycle E → Y → Z → E taking up one first substrate is no
    # potential; without its first transition the cycle is broken and one exists.
    change = zeros(Int, 1, 10)
    change[1, 4] = 1
    turning = EnzymeRates.Derivation._Network(
        net.nforms, net.ncopies, net.present, net.from, net.to, net.edge, net.conformation,
        net.forward, net.fast, net.vanishing, change, net.takes)
    @test EnzymeRates.Derivation._potential_on_sccs(net, trues(10), 1)
    @test !EnzymeRates.Derivation._potential_on_sccs(turning, trues(10), 1)
    broken = BitVector([1, 1, 1, 0, 1, 1, 1, 1, 1, 1])
    @test EnzymeRates.Derivation._potential_on_sccs(turning, broken, 1)
end

@testset "a second free form keeps its law exact once the conformations flip at E" begin
    # Every admitted tag set of `_testhelper_tag_sets` that holds a free form besides E:
    # R1, AD2a, AD2d, AD5, R2P_sealed and B_uphill join it to E by a steady-state or
    # tagged step, A5, AD7a, AD7c and the AD13 family in rapid equilibrium (its label
    # sorting before or after E). The conformations flip at E alone, so every one
    # derives through its derivation entry, and each law equals the exact limit of the
    # oracle's own network, which prunes from free E by itself (R2P_sealed and B_uphill
    # lose their second free form in the inactive conformation).
    tag_sets = [["R1", "AD2a", "AD2d", "AD5", "R2P_sealed", "B_uphill", "A5", "AD7a",
                 "AD7c"];
                ["AD13_$(chem)_$label" for chem in (:NonequalAI, :OnlyA)
                 for label in (:D, :Ca, :Ea, :F, :Z)]]
    for (i, name) in enumerate(tag_sets)
        @testset "$name" begin
            _testhelper_check_law_against_limit(
                _testhelper_tag_sets[name]; seed = i)
        end
    end
end

# ── :OnlyI on modifier-binding groups ─────────────────────────────────────────

@testset "an :OnlyI modifier binding mirrors :OnlyA in the active conformation" begin
    fx = _testhelper_onlyi_fixtures()
    # The fixtures but the inhibitor bridge and the essential activator pass the guard,
    # and the oracles agree on every face of four of them: steady-state and
    # rapid-equilibrium steps, a substrate's modifier copy and an activator route.
    admitted = [k for k in keys(fx) if !(k in (:copy_group_bridge, :essential_activator))]
    @test length(admitted) == 11
    for k in admitted
        @test (k, EnzymeRates.Derivation._guard(fx[k]).fault) ==
              (k, EnzymeRates.Derivation._NO_FAULT)
    end
    for (i, k) in enumerate((:competitive_SS, :uncompetitive_RE, :mod_copy_substrate,
                             :nonessential_activator))
        agree = _testhelper_admission_agreement(
            fx[k]; rng = Random.MersenneTwister(500 + i))
        @test (k, all(agree)) == (k, true)
    end
    # The competitive inhibitor tagged :OnlyI is the :OnlyA one with the conformations
    # swapped: each conformation's law keeps what the other's drops.
    onlya = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + I ⇌ EI            :: OnlyA
        end
    end
    onlyi = fx.competitive_RE
    runs_i = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(onlyi)
    runs_a = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(onlya)
    @test runs_i.active == runs_a.inactive
    @test runs_i.inactive == runs_a.active
    @test !all(runs_i.active)
end

@testset "a fast cycle that only the inactive conformation keeps is named there" begin
    # The two :OnlyI steps each take up I and J at once, so no built reading switches them
    # off; but they join the active conformation's group of EIJS and EIJP to
    # free E's, and consistent readings let them lose their rates into it. With the
    # product binding :OnlyA, the inactive conformation keeps the I, J steps, which still
    # join EP to free E: the product binding lies inside one group, the fault no reading
    # changes. No fast cycle is left to class 2.
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES                    :: EqualAI
            ES + I + J ⇌ EIJS       :: OnlyI
            EIJS ⇌ EIJP         :: EqualAI
            EP + I + J ⇌ EIJP       :: OnlyI
            E + P ⇌ EP                    :: OnlyA
        end
    end
    k = _testhelper_edge(m, "EP + I + J ⇌ EIJP")
    @test EnzymeRates.Derivation._admission_interior(
        m, EnzymeRates.Derivation._built_reading(m)) ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._FAST_STEP_PINNED,
                                          EnzymeRates.Derivation._LIMIT_UNTESTED, 0, 1, k)
    p = _testhelper_edge(m, "E + P ⇌ EP")
    v = EnzymeRates.Derivation._guard(m)
    @test (v.fault, v.face, v.conformation, v.edge) ==
          (EnzymeRates.Derivation._FAST_STEP_PINNED, 0, 2, p)
    @test EnzymeRates.Derivation._guard_message(m) == _testhelper_pinned_cause(m, p, 2)
    # The same with the chemistry steady-state.
    iso_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES                    :: EqualAI
            ES + I + J ⇌ EIJS       :: OnlyI
            EIJS <--> EIJP      :: EqualAI
            EP + I + J ⇌ EIJP       :: OnlyI
            E + P ⇌ EP                    :: OnlyA
        end
    end
    @test EnzymeRates.Derivation._guard(iso_ss).fault ===
          EnzymeRates.Derivation._FAST_STEP_PINNED
    # With the product binding untagged, the inactive conformation keeps every step: the
    # only cycle of rapid-equilibrium steps, which turns the reaction over, runs there, and
    # class 2, which runs before admission, names that conformation.
    kept = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES                    :: EqualAI
            ES + I + J ⇌ EIJS       :: OnlyI
            EIJS ⇌ EIJP         :: EqualAI
            EP + I + J ⇌ EIJP       :: OnlyI
            E + P ⇌ EP                    :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(kept)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._RE_TURNOVER, 0, 2)
    @test EnzymeRates.Derivation._guard_message(kept) ==
          "Mechanism: the rapid-equilibrium steps EP + I + J ⇌ EIJP, " *
          "E + P ⇌ EP, E + S ⇌ ES, ES + I + J ⇌ EIJS and " *
          "EIJS ⇌ EIJP close a cycle in the inactive conformation that turns " *
          "the reaction over, so the mechanism has no finite rate. Make one step of the " *
          "cycle steady-state."
end

@testset "the ε-network's second closed class rejects the face" begin
    # fz101_1964: at zero products the :OnlyA steps vanish in the inactive conformation
    # and the :OnlyI binding in the active one, and the alive transitions of the
    # ε-network leave two closed classes, so the steady state depends on where the enzyme
    # starts; the same at zero products and modifiers. Class 3 passes every face: the law
    # networks keep one class each.
    fz = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: OnlyA
            ES <--> E + P         :: OnlyA
            ES + I ⇌ EIS      :: OnlyI
            EIS ⇌ EI + P      :: NonequalAI
        end
    end
    v = EnzymeRates.Derivation._guard(fz)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._EPS_NONUNIQUE, 2, 1)
    # The class the message names, EI of the active conformation, which only the
    # vanishing I binding at ES feeds, holds no form of the law networks.
    @test EnzymeRates.Derivation._guard_message(fz) ==
          "Mechanism: at zero P, the forms EI (active) hold the enzyme they " *
          "receive and are reached only through tagged steps, so the steady state " *
          "depends on where the enzyme starts."
    eps = EnzymeRates.Derivation._eps_network(
        fz, EnzymeRates.Derivation._re_flags(fz.catalytic),
        EnzymeRates.Derivation._built_reading(fz))
    law = EnzymeRates.Derivation._law_network(fz)
    s = EnzymeRates.Derivation._face_sets(eps, law.present,
                                          EnzymeRates.Derivation._faces(fz.catalytic)[2])
    named = EnzymeRates.Derivation._class_outside_law(
        EnzymeRates.Derivation._closed_classes(eps, s.alive), law.present)
    @test !any(view(law.present, named))
    NO, NU = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._EPS_NONUNIQUE
    @test _testhelper_admission_faults(fz) == [NO, NU, NO, NO, NU, NO]
    @test all(==(NO), _testhelper_face_faults(fz))
    @test _testhelper_epath_verdict(fz; rng = Random.MersenneTwister(65),
                                    face = Set([:P])) === :nonunique
    @test _testhelper_epath_verdict(fz; rng = Random.MersenneTwister(66),
                                    face = Set([:P, :I])) === :nonunique
end

@testset "an :OnlyI bridge is rejected" begin
    # copy_group_bridge and B2b: the :OnlyI bindings of I at E and ES bridge, in the
    # active conformation, forms its other steps keep apart, so the limit depends on how
    # fast the vanishing affinities fall; the ε-path settles on another value. With the
    # bridging group at steady state the guard admits the mechanism, the fix the message
    # prints. At zero modifiers class 3 also rejects copy_group_bridge.
    bridge = _testhelper_onlyi_fixtures().copy_group_bridge
    b2b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S <--> ES                       :: EqualAI
            ES <--> E + P                       :: NonequalAI
            (E + I ⇌ EI, ES + I ⇌ EIS)    :: OnlyI
            EI + S ⇌ EIS                    :: EqualAI
            E + J ⇌ EJ                          :: OnlyA
        end
    end
    for m in (bridge, b2b)
        v = EnzymeRates.Derivation._guard(m)
        @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._BRIDGE, 1, 1)
        @test _testhelper_admission_faults(m)[1:3] ==
              fill(EnzymeRates.Derivation._BRIDGE, 3)
        @test all(!=(EnzymeRates.Derivation._NO_FAULT), _testhelper_face_faults(m)[4:6])
    end
    msg =
        EnzymeRates.Derivation._guard_message(bridge)
    @test startswith(msg, "Mechanism: the :OnlyI step ") &&
          occursin(" bridges forms of the active conformation that its other steps " *
                   "keep apart, so admission cannot show that the limit does not depend " *
                   "on how fast that step vanishes. Make the kinetic group of ", msg) &&
          endswith(msg, " steady-state.")
    # Along the oracle's five orders the ε-path settles on another value or depends on
    # the order.
    @test _testhelper_epath_verdict(bridge; rng = Random.MersenneTwister(71)) in
          (:wrong_limit, :direction_dependent)
    fixed = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> ES                           :: EqualAI
            ES <--> E + P                           :: EqualAI
            (E + I <--> EI, ES + I <--> EIS)  :: OnlyI
            EI + S ⇌ EIS                        :: EqualAI
        end
    end
    @test all(_testhelper_admission_agreement(fixed; rng = Random.MersenneTwister(72)))
    # B2a bridges in both conformations, yet admission passes every face; class 3
    # rejects the three faces with zero modifiers, where the limit differs from the law.
    b2a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S <--> ES                       :: EqualAI
            ES <--> E + P                       :: EqualAI
            (E + I ⇌ EI, ES + I ⇌ EIS)    :: OnlyA
            EI + S ⇌ EIS                    :: EqualAI
            (E + J ⇌ EJ, ES + J ⇌ EJS)    :: OnlyI
            EJ + S ⇌ EJS                    :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(b2a)
    @test v.face == 4 && v.limit === EnzymeRates.Derivation._LIMIT_DIFFER
    @test _testhelper_admission_faults(b2a)[1:3] ==
          fill(EnzymeRates.Derivation._NO_FAULT, 3)
    @test _testhelper_respeed_verdicts(b2a; rng = Random.MersenneTwister(73),
                                       law = :limit) ==
          [:agree, :agree, :agree, :mismatch, :mismatch, :mismatch]
end

@testset "an essential activator tagged :OnlyI stops catalysis" begin
    # X binds only the inactive conformation, and every turnover passes through EX:
    # the active conformation, which must catalyze, turns nothing over. Tagging the X
    # binding :EqualAI or :NonequalAI, the message's fix, passes.
    ea = _testhelper_onlyi_fixtures().essential_activator
    v = EnzymeRates.Derivation._guard(ea)
    @test (v.fault, v.face, v.conformation) ==
          (EnzymeRates.Derivation._NO_ACTIVE_TURNOVER, 0, 1)
    @test EnzymeRates.Derivation._guard_message(ea) ==
          "Mechanism: without its :OnlyI steps E + X ⇌ EX the active " *
          "conformation turns no substrate into product, but it must catalyze the " *
          "reaction. Tag their kinetic groups :EqualAI or :NonequalAI."
    for (i, tag) in enumerate((:EqualAI, :NonequalAI))
        fixed = eval(:(@enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: X
            steps: begin
                E + X ⇌ EX            :: $tag
                EX + S ⇌ ESX      :: EqualAI
                ESX <--> EX + P   :: EqualAI
            end
        end))
        @test all(_testhelper_admission_agreement(fixed;
                                                  rng = Random.MersenneTwister(80 + i)))
    end
end

@testset "an :OnlyI activator whose retags fail is offered no fix" begin
    # X binds E and ES in the inactive conformation only, and P leaves only from
    # ESX, so the active conformation turns nothing over. Tagging the X bindings
    # :EqualAI or :NonequalAI joins E, EX, ESX and ES in one rapid-equilibrium
    # segment that only the X bindings hold together, so each variant splits it at zero
    # X and the message offers no fix.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> ES         :: EqualAI
            ES + X ⇌ ESX      :: OnlyI
            ESX ⇌ EX + P      :: EqualAI
            E + X ⇌ EX            :: OnlyI
        end
    end
    v = EnzymeRates.Derivation._guard(am)
    @test (v.fault, v.face, v.conformation) ==
          (EnzymeRates.Derivation._NO_ACTIVE_TURNOVER, 0, 1)
    ks = findall(EnzymeRates.Derivation._excluded_steps(am, 1))
    @test !EnzymeRates.Derivation._retag_passes(am, ks, :EqualAI)
    @test !EnzymeRates.Derivation._retag_passes(am, ks, :NonequalAI)
    @test EnzymeRates.Derivation._guard_message(am) ==
          "Mechanism: without its :OnlyI steps E + X ⇌ EX; ES + X ⇌ " *
          "ESX the active conformation turns no substrate into product, but it " *
          "must catalyze the reaction."
end

@testset "the active conformation's interior faults name the :OnlyI step" begin
    # The mirror of the :OnlyA interior faults, in conformation 1. Each message offers the
    # fixes whose variants pass. The oracles agree with each verdict.
    rng = Random.MersenneTwister(91)
    FP, U = EnzymeRates.Derivation._FAST_STEP_PINNED, EnzymeRates.Derivation._LIMIT_UNTESTED
    # I binds E and ES; without E + I ⇌ EI, the untagged I binding at ES and S
    # binding at EI still join E to EI, so the I binding cannot lose its affinity,
    # and at rapid equilibrium the active conformation cannot switch it off. Neither fix
    # the message offers passes; untagging it would.
    no_limit = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + I ⇌ EI            :: OnlyI
            EI + S ⇌ EIS      :: EqualAI
            ES + I ⇌ EIS      :: EqualAI
        end
    end
    k = _testhelper_edge(no_limit, "E + I ⇌ EI")
    v = EnzymeRates.Derivation._guard(no_limit)
    @test v == EnzymeRates.Derivation._Verdict(FP, U, 0, 1, k)
    @test EnzymeRates.Derivation._guard_message(no_limit) ==
          _testhelper_pinned_cause(no_limit, k, 1)
    @test all(t -> EnzymeRates.Derivation._retag_passes(no_limit, [k], t),
              (:EqualAI, :NonequalAI))
    @test _testhelper_epath_verdict(no_limit; rng) === :no_path
    # I binds E and J binds EI, both :OnlyI; the untagged J binding at E and I binding
    # at EJ join EIJ back to E, so the J binding at EI cannot lose its affinity,
    # and no built reading switches it off. But EI is cut off alone, and the search lets
    # both :OnlyI bindings lose their rates into it, the oracles agreeing.
    square = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + I ⇌ EI            :: OnlyI
            EI + J ⇌ EIJ      :: OnlyI
            E + J ⇌ EJ            :: EqualAI
            EJ + I ⇌ EIJ      :: EqualAI
        end
    end
    k = _testhelper_edge(square, "EI + J ⇌ EIJ")
    @test EnzymeRates.Derivation._admission_interior(
        square, EnzymeRates.Derivation._built_reading(square)) ==
          EnzymeRates.Derivation._Verdict(FP, U, 0, 1, k)
    @test _testhelper_epath_verdict(square; rng) === :no_path
    @test EnzymeRates.Derivation._guard(square) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(square; rng, corners = true))
    _testhelper_check_law_against_limit(square)
    # EI is reached only as EIJ gives off J, and the active conformation's other
    # steps make EIJ: no built reading switches that J binding off, though the search
    # lets it lose its rate into EI. But the inactive conformation keeps it, and at zero
    # I and J EI there holds the enzyme it receives: class 3 finds a second closed
    # class, as the RE-speed oracle does.
    j_release = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + J ⇌ EJ            :: EqualAI
            EJ + I ⇌ EIJ      :: EqualAI
            EI + J ⇌ EIJ      :: OnlyI
        end
    end
    k = _testhelper_edge(j_release, "EI + J ⇌ EIJ")
    @test EnzymeRates.Derivation._admission_interior(
        j_release, EnzymeRates.Derivation._built_reading(j_release)) ==
          EnzymeRates.Derivation._Verdict(FP, U, 0, 1, k)
    @test _testhelper_epath_verdict(j_release; rng) === :no_path
    v = EnzymeRates.Derivation._guard(j_release)
    @test (v.fault, _testhelper_faces(j_release)[v.face], v.conformation) ==
          (EnzymeRates.Derivation._NONUNIQUE, Set([:I, :J]), 2)
    @test _testhelper_agrees(v.fault, _testhelper_respeed_verdicts(
        j_release; rng = Random.MersenneTwister(92), law = :limit)[v.face])
    @test EnzymeRates.Derivation._guard_message(j_release) ==
          "Mechanism: at zero I and J in the inactive conformation, the forms " *
          "{EI} keep enzyme away from free E: every step that leaves them takes up a " *
          "zeroed metabolite (EI + J ⇌ EIJ), so the steady state is not unique."
    # At steady state the built readings stop that J binding, so EI is reached only
    # through a step that vanishes both ways: an island on the :OnlyI binding's unbound
    # side, which collects the active enzyme. The inactive conformation's second closed
    # class at zero I and J rejects it under every reading, as above.
    island = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + J ⇌ EJ            :: EqualAI
            EJ + I ⇌ EIJ      :: EqualAI
            EI + J <--> EIJ   :: OnlyI
        end
    end
    k = _testhelper_edge(island, "EI + J <--> EIJ")
    @test EnzymeRates.Derivation._admission_interior(
        island, EnzymeRates.Derivation._built_reading(island)) ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._ISLAND, U, 0, 1, k)
    @test _testhelper_epath_verdict(island; rng) in (:direction_dependent, :wrong_limit)
    v = EnzymeRates.Derivation._guard(island)
    @test (v.fault, _testhelper_faces(island)[v.face], v.conformation) ==
          (EnzymeRates.Derivation._NONUNIQUE, Set([:I, :J]), 2)
    @test _testhelper_agrees(v.fault, _testhelper_respeed_verdicts(
        island; rng = Random.MersenneTwister(93), law = :limit)[v.face])
    @test EnzymeRates.Derivation._guard_message(island) ==
          "Mechanism: at zero I and J in the inactive conformation, the forms " *
          "{EI} keep enzyme away from free E: every step that leaves them takes up a " *
          "zeroed metabolite (EI + J <--> EIJ), so the steady state is not unique."
end

@testset "a pinned step and its mirror get one message, but for conformation and tag" begin
    # Either conformation is offered the same fixes, each with that conformation's tag,
    # when both conformations catalyze in every fix's variant, as here, where the S path
    # stays untagged. I binds E :OnlyI, and the J bindings join EI back to E, so the
    # active conformation cannot switch the I binding off. Tagging :OnlyI the steps that
    # join its sides passes (and so would untagging it, which the message does not
    # offer). The mirror image, the I binding :OnlyA, gets the same message in the
    # inactive conformation, with :OnlyA.
    D = EnzymeRates.Derivation
    swapped(text) = replace(text, "the active conformation" => "the inactive conformation",
                            ":OnlyI" => ":OnlyA")
    square = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + I ⇌ EI            :: OnlyI
            EI + J ⇌ EIJ      :: EqualAI
            E + J ⇌ EJ            :: EqualAI
            EJ + I ⇌ EIJ      :: EqualAI
        end
    end
    k = _testhelper_edge(square, "E + I ⇌ EI")
    joins = [_testhelper_edge(square, t) for t in ("E + J ⇌ EJ", "EJ + I ⇌ EIJ",
                                                   "EI + J ⇌ EIJ")]
    v = D._guard(square)
    @test (v.fault, v.conformation, v.edge) == (D._FAST_STEP_PINNED, 1, k)
    msg = D._guard_message(square)
    @test msg == _testhelper_pinned_cause(square, k, 1) * " Tag :OnlyI the steps that " *
                 "join its sides: E + J ⇌ EJ; EJ + I ⇌ EIJ; EI + J ⇌ EIJ."
    @test all(t -> D._retag_passes(square, [k], t), (:EqualAI, :NonequalAI))
    @test D._retag_passes(square, joins, :OnlyI)
    mirror = _testhelper_mirror(square)
    w = D._guard(mirror)
    @test (w.fault, w.conformation, w.edge) == (D._FAST_STEP_PINNED, 2, k)
    @test D._guard_message(mirror) == swapped(msg)
    # The one rule that tells the states apart: the active conformation must catalyze.
    # The steps that join the sides of this I binding include the uptake of S, which
    # joins E to EP in one step, so tagging them :OnlyI leaves the active conformation no
    # turnover: only the mirror image is offered that fix.
    no_limit = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> EP         :: EqualAI
            EP ⇌ E + P            :: EqualAI
            E + I ⇌ EI            :: OnlyI
            EI + P ⇌ EIP          :: EqualAI
            EP + I ⇌ EIP          :: EqualAI
        end
    end
    k = _testhelper_edge(no_limit, "E + I ⇌ EI")
    e = no_limit.catalytic.steps[k]
    joins = D._kept_path(no_limit, 1, e.from, e.to)
    retagged = copy(no_limit.allosteric.allosteric_states)
    retagged[[no_limit.catalytic.steps[j].kinetic_group for j in joins]] .= :OnlyI
    @test D._guard(EnzymeRates.Mechanism(
              no_limit.catalytic, EnzymeRates.Mechanisms.Allostery(retagged, 1,
              EnzymeRates.Mechanisms.RegulatorySite[]))).fault === D._NO_ACTIVE_TURNOVER
    msg = D._guard_message(no_limit)
    mirror = _testhelper_mirror(no_limit)
    @test D._guard_message(mirror) ==
          swapped(msg) * " Tag :OnlyA the steps that join its sides: " *
          "E + S <--> EP; EP + I ⇌ EIP; EI + P ⇌ EIP."
    # The same graph with the substrate and product roles the other way round. The path
    # the message offers is the first breadth-first in step order (`_kept_path`), here
    # through the S binding E + S ⇌ ES: retagged :OnlyI it is pinned in its turn, and
    # retagged :OnlyA it pins the mirror's, so neither image is offered a fix, though the
    # path through the P release would pass for the mirror (ROADMAP, known issue).
    p_release_path = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> E + P         :: EqualAI
            E + I ⇌ EI            :: OnlyI
            EI + S ⇌ EIS          :: EqualAI
            ES + I ⇌ EIS          :: EqualAI
        end
    end
    k = _testhelper_edge(p_release_path, "E + I ⇌ EI")
    e = p_release_path.catalytic.steps[k]
    @test [EnzymeRates.Mechanisms.step_text(p_release_path.catalytic, j)
           for j in D._kept_path(p_release_path, 1, e.from, e.to)] ==
          ["E + S ⇌ ES", "ES + I ⇌ EIS", "EI + S ⇌ EIS"]
    msg = D._guard_message(p_release_path)
    @test msg == _testhelper_pinned_cause(p_release_path, k, 1)
    @test D._guard_message(_testhelper_mirror(p_release_path)) == swapped(msg)
end

@testset ":OnlyI on a step that does more than bind a modifier" begin
    # Any group may be :OnlyI while the active conformation still catalyzes the reaction.
    # Each mechanism the macro admits has a law equal to the limit on every face, by the
    # ε-path oracle and by the RE-speed oracle against the package's law.
    D = EnzymeRates.Derivation
    rng = Random.MersenneTwister(101)
    # Two isomerization paths from ES to EP; the direct one runs in the inactive
    # conformation alone. Its mirror image, the direct path :OnlyA, is the same model
    # with L inverted, and the guard admits it too.
    two_paths = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES              :: EqualAI
            ES <--> EP            :: OnlyI
            ES ⇌ EstarS           :: EqualAI
            EstarS <--> EstarP    :: EqualAI
            EstarP ⇌ EP           :: EqualAI
            E + P ⇌ EP              :: EqualAI
        end
    end
    @test all(_testhelper_admission_agreement(two_paths; rng))
    @test all(==(:agree), _testhelper_respeed_verdicts(two_paths; rng))
    mirror = _testhelper_mirror(two_paths)
    @test D._guard(mirror) == D._WELL_POSED
    @test all(_testhelper_admission_agreement(mirror; rng))
    @test all(==(:agree), _testhelper_respeed_verdicts(mirror; rng))
    # An exchange step beside a ternary-complex path, the exchange :OnlyI.
    exchange = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA              :: EqualAI
            EA + B <--> EQ + P    :: OnlyI
            EA + B ⇌ EAB        :: EqualAI
            EAB <--> EPQ      :: EqualAI
            EQ + P ⇌ EPQ        :: EqualAI
            E + Q ⇌ EQ              :: EqualAI
        end
    end
    @test all(_testhelper_admission_agreement(exchange; rng))
    @test all(==(:agree), _testhelper_respeed_verdicts(exchange; rng))
    # The only catalytic step :OnlyI: the active conformation cannot catalyze.
    no_turnover = "the active conformation turns no substrate into product"
    @test_throws no_turnover @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: OnlyI
            E + P ⇌ EP      :: EqualAI
        end
    end
    # A rapid-equilibrium S binding :OnlyI whose bound form the chemistry and the P
    # binding keep joined to free E: the active conformation cannot switch it off. The
    # message offers no fix: each would leave the active conformation no turnover
    # (untagging the binding passes, but is not one of the fixes it offers).
    pinned_active = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: OnlyI
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: EqualAI
        end
    end
    k = _testhelper_edge(pinned_active, "E + S ⇌ ES")
    v = D._guard(pinned_active)
    @test (v.fault, v.conformation, v.edge) == (D._FAST_STEP_PINNED, 1, k)
    msg = D._guard_message(pinned_active)
    @test occursin("in the active conformation", msg)
    @test occursin(":OnlyI", msg)
    @test msg == _testhelper_pinned_cause(pinned_active, k, 1)
    @test all(t -> D._retag_passes(pinned_active, [k], t), (:EqualAI, :NonequalAI))
    @test _testhelper_epath_verdict(pinned_active; rng) === :no_path
end

@testset "sampled taggings: the rule agrees with the mass-action oracles" begin
    # An aggregate check over a population: 60 taggings of the plain init mechanisms of
    # the bi-bi and ping-pong reactions, sorted by their text so that the sample does not
    # depend on the enumeration order, each built without the guard. 40 tag one random
    # nonempty set of groups :OnlyA; 20 draw a state for each group, at least one :OnlyI.
    # An admitted law agrees with both oracles on every face. A mechanism rejected for a
    # pinned step has a step that the tests' own statement of the rule cannot switch off,
    # and the ε-path oracle finds no path to the limit.
    D = EnzymeRates.Derivation
    bibi = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    pingpong = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    plain = sort!(unique([m for rxn in (bibi, pingpong)
                          for m in _testhelper_mechanisms(
                              EnzymeRates.Enumeration.init_mechanisms(rxn))
                          if !_testhelper_is_allosteric(m)]); by = m -> sprint(show, m))
    rng = Random.MersenneTwister(20261119)
    function _testhelper_random_tagging(holds_only_i::Bool)
        m = rand(rng, plain)
        G = length(_testhelper_group_members(m.catalytic))
        states = if holds_only_i
            s = Symbol[]
            while !(:OnlyI in s)
                s = rand(rng, [:EqualAI, :OnlyA, :OnlyI], G)
            end
            s
        else
            bits = rand(rng, 1:(2^G - 1))
            [isodd(bits >> (g - 1)) ? :OnlyA : :EqualAI for g in 1:G]
        end
        EnzymeRates.Mechanism(m.catalytic, EnzymeRates.Mechanisms.Allostery(states, 1,
                                        EnzymeRates.Mechanisms.RegulatorySite[]))
    end
    sample = [[_testhelper_random_tagging(false) for _ in 1:40];
              [_testhelper_random_tagging(true) for _ in 1:20]]
    faults = [D._guard(am).fault for am in sample]
    for (am, fault) in zip(sample, faults)
        if fault === D._NO_FAULT
            @test all(_testhelper_admission_agreement(am; rng))
        elseif fault === D._FAST_STEP_PINNED
            @test :cannot_switch_off in _testhelper_reading(am, :all_lose_affinity)
            @test _testhelper_epath_verdict(am; rng) === :no_path
        end
    end
    @test count(==(D._NO_FAULT), faults) >= 10
    @test count(==(D._FAST_STEP_PINNED), faults) >= 3
end

@testset "sampled taggings: the energy-level rule agrees with the oracles" begin
    # An aggregate check over a population: the nonempty :OnlyA sets of the plain init
    # mechanisms of 8 steps or fewer of the bi-bi and ping-pong reactions, sorted by their
    # text so that the sample does not depend on the enumeration order, each built without
    # the guard. 40 random ones that the guard admits only by a reading it searched, none
    # of its three built readings passing, agree with both oracles on every face, the
    # ε-path oracle along the corner directions under the reading the guard admitted with.
    # At least 5 of them are admitted with a step that loses its reverse direction, which
    # no built reading does.
    D = EnzymeRates.Derivation
    bibi = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    pingpong = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    plain = sort!(unique([m for rxn in (bibi, pingpong)
                          for m in _testhelper_mechanisms(
                              EnzymeRates.Enumeration.init_mechanisms(rxn))
                          if !_testhelper_is_allosteric(m) &&
                             length(m.catalytic.steps) <= 8]);
                  by = m -> sprint(show, m))
    built(am) = [D._built_reading(am, r) for r in instances(D._SlowUptakeReading)]
    rng = Random.MersenneTwister(20261120)
    sample = EnzymeRates.Mechanism[]
    # The 40 take 1,448 draws; the bound only stops a guard that admits fewer from
    # drawing forever.
    for _ in 1:10_000
        length(sample) == 40 && break
        m = rand(rng, plain)
        G = length(_testhelper_group_members(m.catalytic))
        bits = rand(rng, 1:(2^G - 1))
        am = EnzymeRates.Mechanism(
            m.catalytic,
            EnzymeRates.Mechanisms.Allostery(
                [isodd(bits >> (g - 1)) ? :OnlyA : :EqualAI for g in 1:G], 1,
                EnzymeRates.Mechanisms.RegulatorySite[]))
        g = D._guard_with_reading(am)
        searched = g.verdict === D._WELL_POSED && !(g.reading in built(am))
        searched && !(am in sample) && push!(sample, am)
    end
    @test length(sample) == 40
    for am in sample
        @test all(_testhelper_admission_agreement(am; rng, corners = true))
    end
    admitting(am) = _testhelper_reading(D._guard_with_reading(am).reading)
    @test count(am -> :loses_reverse in admitting(am), sample) >= 5
end

@testset "joint R1 rejects a mechanism whose bridges cut the tree only together" begin
    # The hand-built joint-R1 network's pattern on a mechanism. Half 2 :OnlyA strands the
    # inactive conformation's covalent forms at zero P, so no in-tree of that
    # conformation reaches its free E without half 2's vanishing steps. M binds EA and
    # the covalent form :OnlyI, and M-bound half 1 runs at rapid equilibrium in the active
    # conformation, where the vanishing M binding at EA lets free E's fast class leave,
    # one way, for the M-bound covalent form. Each conformation's vanishing transitions
    # alone leave the cheapest spanning in-tree at 4 slow transitions; together they cut
    # it to 3. The guard rejects the interior first: the M bindings, which the inactive
    # conformation keeps, join the two sides of M-bound half 1, an :OnlyA step at rapid
    # equilibrium, so that conformation cannot switch it off.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: M
        steps: begin
            E + A ⇌ EA                                             :: EqualAI
            EA <--> E_res_plusA_minusP + P                      :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: EqualAI
            EB_res_plusA_minusP <--> E + Q                        :: OnlyA
            EA + M ⇌ EAM                                       :: OnlyI
            E_res_plusA_minusP + M <--> EM_res_plusA_minusP    :: OnlyI
            EAM ⇌ EM_res_plusA_minusP + P                     :: OnlyA
        end
    end
    NO, BR, U = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._BRIDGE,
        EnzymeRates.Derivation._LIMIT_UNTESTED
    k = _testhelper_edge(am, "EM_res_plusA_minusP + P ⇌ EAM")
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._FAST_STEP_PINNED, U,
                                               0, 2, k)
    @test EnzymeRates.Derivation._guard_message(am) ==
          _testhelper_pinned_cause(am, k, 2) * " Make its kinetic group steady-state."
    @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(95)) === :no_path
    @test _testhelper_admission_faults(am) == [[NO, BR, NO, NO, NO, BR, NO]; fill(NO, 7)]
    eps = EnzymeRates.Derivation._eps_network(
        am, EnzymeRates.Derivation._re_flags(am.catalytic),
        EnzymeRates.Derivation._built_reading(am))
    law = EnzymeRates.Derivation._law_network(am)
    faces = EnzymeRates.Derivation._faces(am.catalytic)
    for i in (2, 6)   # zero P and Q, zero P
        s = EnzymeRates.Derivation._face_sets(eps, law.present, faces[i])
        nodes = s.fed .& .!s.stranded
        @test EnzymeRates.Derivation._fewest_slow(eps, nodes, s.n0) == 4
        for c in 1:2
            (; van) = EnzymeRates.Derivation._copy_transitions(eps, s.alive, nodes, c)
            @test EnzymeRates.Derivation._fewest_slow(eps, nodes, s.n0 .| van) == 4
        end
        @test EnzymeRates.Derivation._fewest_slow(eps, nodes, s.alive) == 3
        @test EnzymeRates.Derivation._face_condition(
            eps, law.present, faces[i],
            first(am.catalytic.substrates)) ==
              EnzymeRates.Derivation._Verdict(BR, U, 0, 0, 0)
    end
    # The face condition's message names the vanishing steps a cheapest in-tree needs: the
    # M binding at EA of the active conformation and half 2 of the inactive one.
    w = EnzymeRates.Derivation._admission_face(am, eps, law, faces, 2)
    @test w == EnzymeRates.Derivation._Verdict(BR, U, 2, 0, 0)
    @test EnzymeRates.Derivation._admission_message(
        am, w, EnzymeRates.Derivation._built_reading(am)) ==
          "Mechanism: at zero P and Q, the tagged steps EA + M ⇌ EAM; " *
          "E + Q <--> EB_res_plusA_minusP of both conformations together bridge forms " *
          "that each conformation's other steps keep apart, so admission cannot show " *
          "that the limit does not depend on how fast those steps vanish."
end

# ── Admission on written mechanisms ───────────────────────────────────────────────────

@testset "a second free form is admitted where a reading drains it" begin
    # The conformations flip at free E alone. A second free form joined to E by a
    # steady-state step is an ordinary form, whatever its label (R1, and its twin
    # labelled D, which sorts before E). Where tags strand it under the guard's built
    # readings, which stop a tagged isomerization and cannot switch off a
    # rapid-equilibrium one, admission's verdicts under them stand: every step at E
    # :OnlyA leaves F fed at zero product with a dead way back (`_PUMPED`); a tagged
    # steady-state isomerization alone joining F to E makes F an island (`_ISLAND`);
    # a tagged rapid-equilibrium one cannot be switched off (`_FAST_STEP_PINNED`). The
    # search admits each under a reading in which the isomerization loses its rate
    # into F, so F drains into free E, the oracles agreeing on every face. Each also
    # has an admitted variant under the built readings: the S binding untagged (AD2d),
    # or the catalysis tagged in place of the isomerization.
    r1 = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: NonequalAI
            F + P ⇌ EP      :: EqualAI
            F <--> E          :: NonequalAI
        end
    end
    r1_d = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: NonequalAI
            D + P ⇌ EP      :: EqualAI
            D <--> E          :: NonequalAI
        end
    end
    # R1's oracle is the RE-speed one (no tag excludes a step); the law also equals the
    # exact limit at one subunit, where L means [E_I]/[E_A] for D too.
    for (i, m) in enumerate((r1, r1_d))
        @test all(==(:agree), _testhelper_respeed_verdicts(
            m; rng = Random.MersenneTwister(20 + i)))
        _testhelper_check_law_against_limit(m; seed = 20 + i)
    end
    r2 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            F + P ⇌ EP      :: EqualAI
            F <--> E          :: OnlyA
        end
    end
    v = ER.Derivation._verdict_under(r2, ER.Derivation._built_reading(r2))
    @test (v.fault, v.face) == (ER.Derivation._PUMPED, 2)
    @test _testhelper_epath_verdict(r2; rng = Random.MersenneTwister(23),
                                    face = Set([:P])) === :direction_dependent
    @test @enzyme_mechanism(begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            F + P ⇌ EP      :: EqualAI
            F <--> E          :: OnlyA
        end
    end) == r2
    ad2d = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
            F + P ⇌ EP      :: EqualAI
            F <--> E          :: OnlyA
        end
    end
    @test _testhelper_epath_verdict(ad2d;
                                    rng = Random.MersenneTwister(24)) === :converges
    r3 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: EqualAI
            E <--> F          :: OnlyA
            F + I ⇌ FI      :: EqualAI
        end
    end
    r4 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: NonequalAI
            E + P ⇌ EP      :: EqualAI
            F ⇌ E             :: OnlyA
        end
    end
    @test ER.Derivation._admission_interior(
        r3, ER.Derivation._built_reading(r3)).fault === ER.Derivation._ISLAND
    @test _testhelper_epath_verdict(r3; rng = Random.MersenneTwister(26)) in
          (:direction_dependent, :wrong_limit)
    v = ER.Derivation._admission_interior(r4, ER.Derivation._built_reading(r4))
    @test (v.fault, v.conformation) == (ER.Derivation._FAST_STEP_PINNED, 2)
    @test ER.Mechanisms.step_text(r4.catalytic, v.edge) == "E ⇌ F"
    @test _testhelper_epath_verdict(r4; rng = Random.MersenneTwister(27)) === :no_path
    for (i, x) in enumerate((r2, r3, r4))
        @test ER.Derivation._guard(x) == ER.Derivation._WELL_POSED
        @test all(_testhelper_admission_agreement(
            x; rng = Random.MersenneTwister(30 + i), corners = true))
        _testhelper_check_law_against_limit(x)
    end
    @test @enzyme_mechanism(begin
        substrates: S
        products: P
        modifiers: I
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: EqualAI
            E <--> F          :: OnlyA
            F + I ⇌ FI      :: EqualAI
        end
    end) == r3
    @test @enzyme_mechanism(begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: NonequalAI
            E + P ⇌ EP      :: EqualAI
            F ⇌ E             :: OnlyA
        end
    end) == r4
    r3_base = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: OnlyA
            E + P ⇌ EP      :: EqualAI
            E <--> F          :: EqualAI
            F + I ⇌ FI      :: EqualAI
        end
    end
    r4_base = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: OnlyA
            E + P ⇌ EP      :: EqualAI
            F ⇌ E             :: EqualAI
        end
    end
    for (i, m) in enumerate((r3_base, r4_base))
        @test _testhelper_epath_verdict(
            m; rng = Random.MersenneTwister(27 + i)) ===
            :converges
    end
    # F in E's rapid-equilibrium segment, untagged, is an ordinary form too (A5).
    @test _testhelper_is_allosteric(@enzyme_mechanism(begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: NonequalAI
            E + P ⇌ EP      :: EqualAI
            F ⇌ E             :: NonequalAI
        end
    end))
end

@testset "admission rejects a one-sided :OnlyA binding and admits balanced ones" begin
    # Admission's fault for `am`'s graph tagged `tags`, and the metabolite of the
    # binding it names (`nothing` when it names none).
    function _testhelper_admission_reads(am, tags)
        cm = am.catalytic
        tagged = ER.Mechanism(cm, ER.Mechanisms.Allostery(tags, 1,
                                                          ER.Mechanisms.RegulatorySite[]))
        v = ER.Derivation._verdict_under(tagged, ER.Derivation._built_reading(tagged))
        e = v.edge == 0 ? nothing : cm.steps[v.edge]
        (v.fault, e === nothing || !ER.Enumeration._is_binding(e) ? nothing :
                  cm.metabolites[only(e.taken_up)])
    end
    NO, FP = ER.Derivation._NO_FAULT, ER.Derivation._FAST_STEP_PINNED

    # Uni-uni S -> P. Tags: (S binding, chemical step, P binding). Admission's fault,
    # and the metabolite of the binding it names. Each binding is at rapid equilibrium,
    # so the inactive conformation can switch off an :OnlyA binding only by its affinity,
    # which needs its two sides apart.
    function _testhelper_uni(s_tag, cat_tag, p_tag)
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: EqualAI
            end
        end
        firsts = _testhelper_first_steps(_testhelper_graph(m))
        # tag each group by the metabolite its step binds
        tags = [!_testhelper_is_binding(s) ? cat_tag :
                only(s.consumed) === :S ? s_tag : p_tag for s in firsts]
        _testhelper_admission_reads(m, tags)
    end

    # no :OnlyA anywhere -> admitted
    @test _testhelper_uni(:EqualAI, :EqualAI, :EqualAI) == (NO, nothing)
    # :OnlyA on the substrate only, catalysis :EqualAI: catalysis and the P binding
    # join E to ES in the inactive conformation -> the S binding cannot be switched off
    @test _testhelper_uni(:OnlyA, :EqualAI, :EqualAI) == (FP, :S)
    # :OnlyA on the product only, catalysis :EqualAI -> the P binding cannot be switched
    # off
    @test _testhelper_uni(:EqualAI, :EqualAI, :OnlyA) == (FP, :P)
    # :OnlyA on the substrate and the catalysis: ES lies behind the binding
    @test _testhelper_uni(:OnlyA, :OnlyA, :EqualAI) == (NO, nothing)
    # balanced: :OnlyA on both bindings, catalysis :EqualAI -> admitted
    @test _testhelper_uni(:OnlyA, :EqualAI, :OnlyA) == (NO, nothing)
    # V-system: :OnlyA chemical step only -> admitted
    @test _testhelper_uni(:EqualAI, :OnlyA, :EqualAI) == (NO, nothing)
    # :NonequalAI catalysis joins E to ES too -> same verdict
    @test _testhelper_uni(:OnlyA, :NonequalAI, :EqualAI) == (FP, :S)

    # Whether a binding is rapid-equilibrium or steady-state does not change
    # admission's verdict: with both bindings :OnlyA, ES and EP lie behind them.
    # Both orientations of the mix are admitted.
    mixed_uni_re_ss = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
            EP <--> E + P   :: EqualAI
        end
    end
    mixed_uni_ss_re = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + S <--> ES   :: EqualAI
            ES <--> EP    :: EqualAI
            EP ⇌ E + P      :: EqualAI
        end
    end
    # every binding :OnlyA, catalysis :EqualAI -> balanced -> valid
    function _testhelper_both_bindings_onlya(m)
        tags = [_testhelper_is_binding(s) ? :OnlyA : :EqualAI
                for s in _testhelper_first_steps(_testhelper_graph(m))]
        _testhelper_admission_reads(m, tags)
    end
    @test _testhelper_both_bindings_onlya(mixed_uni_re_ss) == (NO, nothing)
    @test _testhelper_both_bindings_onlya(mixed_uni_ss_re) == (NO, nothing)

    # A random-order binding square: E binds A and B, and each of EA and EB binds
    # the other. A lone :OnlyA edge has both sides joined by the other three, so the
    # inactive conformation cannot switch it off; tagging P's binding too separates them.
    @testset "random-order binding square" begin
        biuni = @enzyme_mechanism begin
            substrates: A, B
            products:   P
            oligomeric_state: 2
            steps: begin
                E + A ⇌ EA          :: EqualAI
                E + B ⇌ EB          :: EqualAI
                EA + B ⇌ EAB    :: EqualAI
                EB + A ⇌ EAB    :: EqualAI
                EAB <--> EP     :: EqualAI
                EP ⇌ E + P          :: EqualAI
            end
        end
        # Tag :OnlyA the groups named by (free form, bound metabolite); every
        # other group is :EqualAI. The chemical step's key is (:EAB, nothing).
        bu = _testhelper_graph(biuni)
        # Each group's key: the form its first step starts from, and the metabolite it
        # binds (`nothing` for the chemical step).
        bu_keys = [(Symbol(bu.forms[s.from]),
                    _testhelper_is_binding(s) ? only(s.consumed) : nothing)
                   for s in _testhelper_first_steps(bu)]
        _testhelper_tags(onlya_keys...) =
            [key in onlya_keys ? :OnlyA : :EqualAI for key in bu_keys]
        _testhelper_verdict(t) = _testhelper_admission_reads(biuni, t)

        # the square alone trips nothing: no :OnlyA anywhere -> admitted
        @test _testhelper_verdict(_testhelper_tags()) == (NO, nothing)
        # one square edge :OnlyA: the square's other three bindings join its sides
        @test _testhelper_verdict(_testhelper_tags((:E, :A))) == (FP, :A)
        # ...and the macro rejects that tagging: untagged steps join the two sides of the
        # rapid-equilibrium :OnlyA binding of A, so the inactive conformation cannot
        # switch it off (`_FAST_STEP_PINNED`).
        lone = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: A, B
            products:   P
            oligomeric_state: 2
            steps: begin
                E + A ⇌ EA          :: OnlyA
                E + B ⇌ EB          :: EqualAI
                EA + B ⇌ EAB    :: EqualAI
                EB + A ⇌ EAB    :: EqualAI
                EAB <--> EP     :: EqualAI
                EP ⇌ E + P          :: EqualAI
            end
        end
        v = ER.Derivation._guard(lone)
        @test v.fault === ER.Derivation._FAST_STEP_PINNED
        @test ER.Mechanisms.step_text(lone.catalytic, v.edge) == "E + A ⇌ EA"
        @test_throws ErrorException(
            ER.Derivation._guard_message(lone)) @enzyme_mechanism begin
            substrates: A, B
            products:   P
            oligomeric_state: 2
            steps: begin
                E + A ⇌ EA          :: OnlyA
                E + B ⇌ EB          :: EqualAI
                EA + B ⇌ EAB    :: EqualAI
                EB + A ⇌ EAB    :: EqualAI
                EAB <--> EP     :: EqualAI
                EP ⇌ E + P          :: EqualAI
            end
        end
        # both A-side bindings :OnlyA: the B bindings, the catalysis and the P binding
        # still join E to EA -> they cannot be switched off
        @test _testhelper_verdict(_testhelper_tags((:E, :A), (:EB, :A)))[1] === FP
        # ...and :OnlyA on P cuts the last join -> admitted
        @test _testhelper_verdict(
            _testhelper_tags((:E, :A), (:EB, :A), (:E, :P))) == (NO, nothing)
        # :OnlyA on the chemical step leaves the square, whose lone :OnlyA edge still
        # has its sides joined -> that edge cannot be switched off
        @test _testhelper_verdict(_testhelper_tags((:EAB, nothing), (:E, :A))) ==
              (FP, :A)
    end
end

@testset "admission rejects a ter-substrate cube joining a tagged binding's two sides" begin
    # A full random-order ter binding cube with seven :OnlyA edges and the cross bindings
    # :EqualAI. In the inactive conformation the untagged bindings join forms across the
    # cube, so a tagged binding's two sides share a component: no path takes its affinity to
    # zero. Admission finds untagged steps that join the two sides of an :OnlyA binding at
    # rapid equilibrium, which the inactive conformation cannot switch off
    # (`_FAST_STEP_PINNED`).
    cube = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, C
        products:   P
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA             :: OnlyA
            E + B ⇌ EB             :: OnlyA
            E + C ⇌ EC             :: OnlyA
            EA + B ⇌ EAB       :: OnlyA
            EA + C ⇌ EAC       :: EqualAI
            EB + A ⇌ EAB       :: EqualAI
            EB + C ⇌ EBC       :: EqualAI
            EC + A ⇌ EAC       :: EqualAI
            EC + B ⇌ EBC       :: EqualAI
            EAB + C ⇌ EABC :: OnlyA
            EAC + B ⇌ EABC :: OnlyA
            EBC + A ⇌ EABC :: OnlyA
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
        end
    end
    v = ER.Derivation._guard(cube)
    @test v.fault === ER.Derivation._FAST_STEP_PINNED
    @test_throws ErrorException(
        ER.Derivation._guard_message(cube)) @enzyme_mechanism begin
        substrates: A, B, C
        products:   P
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA             :: OnlyA
            E + B ⇌ EB             :: OnlyA
            E + C ⇌ EC             :: OnlyA
            EA + B ⇌ EAB       :: OnlyA
            EA + C ⇌ EAC       :: EqualAI
            EB + A ⇌ EAB       :: EqualAI
            EB + C ⇌ EBC       :: EqualAI
            EC + A ⇌ EAC       :: EqualAI
            EC + B ⇌ EBC       :: EqualAI
            EAB + C ⇌ EABC :: OnlyA
            EAC + B ⇌ EABC :: OnlyA
            EBC + A ⇌ EABC :: OnlyA
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
        end
    end
end

# The rejection testset above pins admission's reject direction on a ter-substrate
# binding cube; this one pins its accept direction there. An over-rejection does not
# fail loudly: it silently shrinks the enumeration's search space by dropping valid
# mechanisms before they are fitted.
@testset "admission admits feasible ter-substrate mechanisms" begin
    # A full random-order A/B/C binding lattice, every substrate binding and the
    # chemical step :OnlyA: in the inactive conformation every form stands alone but E
    # and EP, the bindings point away from free E and close no cycle, so admission
    # accepts the cube.
    cube = @enzyme_mechanism begin
        substrates: A, B, C
        products:   P
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA             :: OnlyA
            E + B ⇌ EB             :: OnlyA
            E + C ⇌ EC             :: OnlyA
            EA + B ⇌ EAB       :: OnlyA
            EA + C ⇌ EAC       :: OnlyA
            EB + A ⇌ EAB       :: OnlyA
            EB + C ⇌ EBC       :: OnlyA
            EC + A ⇌ EAC       :: OnlyA
            EC + B ⇌ EBC       :: OnlyA
            EAB + C ⇌ EABC :: OnlyA
            EAC + B ⇌ EABC :: OnlyA
            EBC + A ⇌ EABC :: OnlyA
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
        end
    end
    @test _testhelper_is_allosteric(cube)

    # Ordered ter with an :OnlyA B binding and an :OnlyA chemical step: EAB and
    # EABC lie behind the B binding, on its bound side, so this is admitted.
    ordered = @enzyme_mechanism begin
        substrates: A, B, C
        products:   P
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA             :: EqualAI
            EA + B ⇌ EAB       :: OnlyA
            EAB + C ⇌ EABC :: EqualAI
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
        end
    end
    @test _testhelper_is_allosteric(ordered)
end

@testset "a ping-pong with both fused releases :OnlyA, pumped as first read, passes" begin
    # EA → E_res_plusA_minusP + P runs the first half-reaction and releases P in
    # one step, EB_res_plusA_minusP → E + Q the second and releases Q. In the inactive
    # conformation the covalent forms are cut off. Under the guard's first built reading
    # the second half, read as Q taken up into EB_res_plusA_minusP, loses its affinity,
    # and the first half, whose far side EA stays, is stopped. At zero A and B the B
    # binding feeds E_res_plusA_minusP, which leaves only through the vanishing first
    # half, so admission cannot show that its weight does not depend on how fast that
    # step vanishes (`_PUMPED`), and the ε-path oracle finds a limit that depends on the
    # direction there. The search admits the ping-pong under a reading in which the first
    # half loses its rate into the covalent forms instead, both oracles agreeing on every
    # face. Either half alone :OnlyA is admitted too.
    both = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: OnlyA
        end
    end
    v = ER.Derivation._verdict_under(both, ER.Derivation._built_reading(both))
    @test (v.fault, v.face, v.conformation) == (ER.Derivation._PUMPED, 3, 2)
    @test _testhelper_epath_verdict(both; rng = Random.MersenneTwister(6),
                                    face = Set([:A, :B])) in
          (:direction_dependent, :wrong_limit)
    @test ER.Derivation._guard(both) == ER.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(both; rng = Random.MersenneTwister(9),
                                              corners = true))
    _testhelper_check_law_against_limit(both)
    @test @enzyme_mechanism(begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: OnlyA
        end
    end) == both
    half1 = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: EqualAI
        end
    end
    half2 = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: EqualAI
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: OnlyA
        end
    end
    # Each face: the RE-speed oracle (class 3) and the ε-path (admission) agree with the
    # guard, which accepts it.
    for (i, half) in enumerate((half1, half2))
        @test all(_testhelper_admission_agreement(half;
                                                  rng = Random.MersenneTwister(6 + i)))
    end
end

@testset "admission admits an abortive segment an :OnlyA group strands" begin
    # The A group mixes the fused binding into EPQ with A's abortive bindings at
    # EP and EQ. Tagged :OnlyA, it strands {EAP, EAQ, EAPQ} in the
    # inactive conformation: a rapid-equilibrium segment whose every form holds a
    # product, which a mechanism may not have. The product releases still make EPQ
    # there, so no built reading switches off the rapid-equilibrium fused binding. It
    # joins EB, cut off, to free E's group, and the first consistent reading lets it
    # lose its rate into EB; at zero A and B its uptake dies, and it runs one way into
    # EB, which empties only through the slow E + B <--> EB: EB fills faster than it
    # empties (`_FILLED`), as the ε-path oracle finds under that reading.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + B <--> EB                                              :: OnlyA
            (EB + A ⇌ EPQ, EP + A ⇌ EAP, EQ + A ⇌ EAQ) :: OnlyA
            EPQ <--> EP + Q                                        :: EqualAI
            EPQ <--> EQ + P                                        :: EqualAI
            (E + P ⇌ EP, EAQ + P ⇌ EAPQ)                     :: EqualAI
            (E + Q ⇌ EQ, EAP + Q ⇌ EAPQ)                     :: EqualAI
        end
    end
    k = _testhelper_edge(am, "EB + A ⇌ EPQ")
    @test ER.Derivation._admission_interior(am, ER.Derivation._built_reading(am)) ==
          ER.Derivation._Verdict(ER.Derivation._FAST_STEP_PINNED,
                                 ER.Derivation._LIMIT_UNTESTED, 0, 2, k)
    @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(3)) === :no_path
    v = ER.Derivation._guard(am)
    @test (v.fault, _testhelper_faces(am)[v.face], v.conformation, v.edge) ==
          (ER.Derivation._FILLED, Set([:A, :B]), 2, k)
    reported = _testhelper_reading(ER.Derivation._guard_with_reading(am).reading)
    @test _testhelper_admission_agrees(
        v.fault, _testhelper_epath_verdict(am; rng = Random.MersenneTwister(4),
                                           face = Set([:A, :B]), reading = reported))
    @test ER.Derivation._guard_message(am) ==
          "Mechanism: at zero A and B, in the inactive conformation the :OnlyA " *
          "rapid-equilibrium step EB + A ⇌ EPQ runs one way into EB, at a rate " *
          "that vanishes yet grows with the rapid-equilibrium speed, and EB leaves " *
          "only through the slower step E + B <--> EB, so EB fills faster than it " *
          "empties and its weight depends on how fast the tagged step vanishes."
    # With the fused binding steady-state in a group of its own, the inactive
    # conformation stops it, and admission reads each edge: the abortive bindings reach
    # {EAP, EAQ, EAPQ} from the product forms on their unbound side only,
    # so they lose their affinity, the inactive conformation holds no enzyme there, and
    # every face passes. Both oracles agree on every face.
    fused_ss = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + B <--> EB                                              :: OnlyA
            EB + A <--> EPQ                                        :: OnlyA
            (EP + A ⇌ EAP, EQ + A ⇌ EAQ)                     :: OnlyA
            EPQ <--> EP + Q                                        :: EqualAI
            EPQ <--> EQ + P                                        :: EqualAI
            (E + P ⇌ EP, EAQ + P ⇌ EAPQ)                     :: EqualAI
            (E + Q ⇌ EQ, EAP + Q ⇌ EAPQ)                     :: EqualAI
        end
    end
    @test all(_testhelper_admission_agreement(fused_ss; rng = Random.MersenneTwister(4)))
end

# ── The steps a conformation excludes ─────────────────────────────────────────────────

@testset "excluded steps of a plain and an allosteric mechanism" begin
    D = ER.Derivation
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES          :: EqualAI
            ES <--> EP        :: OnlyA
            EP ⇌ E + P          :: EqualAI
            E + I ⇌ EI          :: OnlyI
            EI + P ⇌ EIP    :: EqualAI
        end
    end
    m = _testhelper_plain_of(am)
    # A plain mechanism excludes no step. The steps are stored as in the
    # steps_running_in_allosteric_conformations testset (test_mechanism_type.jl): the active
    # conformation excludes the :OnlyI binding of I, the inactive one the :OnlyA chemistry.
    @test D._excluded_steps(m, 1) == falses(5)
    @test D._excluded_steps(m, 2) == falses(5)
    @test D._excluded_steps(am, 1) == BitVector([1, 0, 0, 0, 0])
    @test D._excluded_steps(am, 2) == BitVector([0, 0, 0, 0, 1])
end

@testset "the excluded steps and exclusion_components: an :OnlyA dead-end binding" begin
    M = ER.Mechanisms
    D = ER.Derivation
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                    :: EqualAI
            ES + S ⇌ ESSmod    :: OnlyA
            ES <--> EP                  :: EqualAI
            EP ⇌ E + P                    :: EqualAI
        end
    end
    m = am.catalytic
    # Canonical form numbers the forms E, EP, ES, ESSmod and stores the steps
    # E + P → EP, E + S → ES, ES → EP, ES + S → ESSmod. The inactive
    # conformation excludes the :OnlyA binding, which leaves ESSmod alone.
    @test M.step_text(m, 4) == "ES + S ⇌ ESSmod"
    @test D._excluded_steps(am, 1) == [false, false, false, false]
    @test D._excluded_steps(am, 2) == [false, false, false, true]
    @test D.exclusion_components(m, D._excluded_steps(am, 1)) == [1, 1, 1, 1]
    @test D.exclusion_components(m, D._excluded_steps(am, 2)) == [1, 1, 1, 4]
    # The rule by tag, one tag per group: an :OnlyA step leaves the inactive conformation,
    # an :OnlyI step the active one; an :EqualAI or :NonequalAI step stays in both.
    tags = [:OnlyA, :OnlyI, :EqualAI, :NonequalAI]
    every_tag = ER.Mechanism(m, ER.Mechanisms.Allostery(tags, 2, M.RegulatorySite[]))
    step_tags = [tags[e.kinetic_group] for e in m.steps]
    @test step_tags[D._excluded_steps(every_tag, 1)] == [:OnlyI]
    @test step_tags[D._excluded_steps(every_tag, 2)] == [:OnlyA]
end

# ── The gate: a dead-inactive ping-pong is rejected; its half-2 variant derives ──────
# With the F6P binding and both chemical steps :OnlyA the inactive conformation runs
# neither half-reaction. The half-2 variant, its last step steady-state, derives: a finite
# rate as a product goes to zero, a rate of 0 with both products at zero (its inactive
# conformation holds the enzyme in the covalent forms there), a forward rate far from
# equilibrium, and the non-allosteric twin's rate at L = 0.
@testset "a dead-inactive ping-pong is rejected; its half-2 variant is its L = 0 twin" begin
    # With the F6P binding and both half-reactions :OnlyA, the inactive conformation runs
    # neither half: the covalent forms between the halves join free E only through tagged
    # chemistry. Under the reading the guard reports, at zero products F16BP's release
    # feeds E_res_plusATP_minusF16BP, which leaves only through the vanishing F6P
    # binding (`_PUMPED`), and the ε-path oracle finds the limit direction-dependent there.
    # The guard finds no consistent reading that passes; under one of them it refuses a
    # filling at zero substrates on which the ε-path oracle converges, a conservative
    # rejection. With the second half alone :OnlyA, EF6P_res_plusATP_minusF16BP ⇌
    # EADP is a rapid-equilibrium step whose sides the inactive conformation's other
    # steps join, which no reading switches off; with that step steady-state the mechanism
    # is admitted (half 1 alone fails class 3, as it does for the ping-pong above).
    dead = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        oligomeric_state: 1
        steps: begin
            E + ATP ⇌ EATP                                                      :: EqualAI
            EATP <--> EF16BP_res_plusATP_minusF16BP                          :: OnlyA
            E_res_plusATP_minusF16BP + F16BP ⇌ EF16BP_res_plusATP_minusF16BP:: EqualAI
            E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP    :: OnlyA
            EF6P_res_plusATP_minusF16BP ⇌ EADP                               :: OnlyA
            E + ADP ⇌ EADP                                                      :: EqualAI
        end
    end
    reported = ER.Derivation._guard_with_reading(dead)
    v = reported.verdict
    @test (v.fault, v.face, v.conformation) == (ER.Derivation._PUMPED, 2, 2)
    msg = ER.Derivation._guard_message(dead)
    @test occursin("feeds the forms E_res_plusATP_minusF16BP of the inactive " *
                   "conformation, which leave only through the vanishing step " *
                   "E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP",
                   msg)
    @test _testhelper_epath_verdict(dead; rng = MersenneTwister(997),
                                    face = _testhelper_faces(dead)[v.face], corners = true,
                                    reading = _testhelper_reading(reported.reading)) ===
          :direction_dependent
    # The conservative rejection: of the three consistent readings, one lets the ε-path
    # oracle converge on every face, and under it the face checks refuse the filling of
    # E_res_plusATP_minusF16BP at zero substrates.
    readings = _testhelper_consistent_readings(dead)
    @test length(readings) == 3
    converging = filter(readings) do reading
        all(face -> _testhelper_epath_verdict(dead; rng = MersenneTwister(998), face,
                                              corners = true, reading) === :converges,
            _testhelper_faces(dead))
    end
    @test length(converging) == 1
    @test _testhelper_admission_faults(dead; reading = only(converging))[3] ===
          ER.Derivation._FILLED
    @test_throws ErrorException(msg) @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        oligomeric_state: 1
        steps: begin
            E + ATP ⇌ EATP                                                      :: EqualAI
            EATP <--> EF16BP_res_plusATP_minusF16BP                          :: OnlyA
            E_res_plusATP_minusF16BP + F16BP ⇌ EF16BP_res_plusATP_minusF16BP:: EqualAI
            E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP    :: OnlyA
            EF6P_res_plusATP_minusF16BP ⇌ EADP                               :: OnlyA
            E + ADP ⇌ EADP                                                      :: EqualAI
        end
    end
    pinned = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        oligomeric_state: 1
        steps: begin
            E + ATP ⇌ EATP                                                      :: EqualAI
            EATP <--> EF16BP_res_plusATP_minusF16BP                          :: EqualAI
            E_res_plusATP_minusF16BP + F16BP ⇌ EF16BP_res_plusATP_minusF16BP:: EqualAI
            E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP    :: EqualAI
            EF6P_res_plusATP_minusF16BP ⇌ EADP                               :: OnlyA
            E + ADP ⇌ EADP                                                      :: EqualAI
        end
    end
    @test ER.Derivation._guard(pinned) ==
          ER.Derivation._Verdict(ER.Derivation._FAST_STEP_PINNED,
                                 ER.Derivation._LIMIT_UNTESTED, 0, 2,
                                 _testhelper_edge(pinned,
                                                  "EF6P_res_plusATP_minusF16BP ⇌ EADP"))
    half2 = @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        oligomeric_state: 1
        steps: begin
            E + ATP ⇌ EATP                                                      :: EqualAI
            EATP <--> EF16BP_res_plusATP_minusF16BP                          :: EqualAI
            E_res_plusATP_minusF16BP + F16BP ⇌ EF16BP_res_plusATP_minusF16BP:: EqualAI
            E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP    :: EqualAI
            EF6P_res_plusATP_minusF16BP <--> EADP                            :: OnlyA
            E + ADP ⇌ EADP                                                      :: EqualAI
        end
    end
    _testhelper_check_law_against_oracle(half2)
    afp = parameters(half2)
    # Non-allosteric twin: the same six steps with no tags. At L = 0 the inactive
    # conformation is unpopulated, so the allosteric rate must reduce to the twin's, an
    # independent re-derivation through the plain King–Altman path. An :OnlyA constant
    # carries the A_ tag its twin lacks.
    nonallo = @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        steps: begin
            E + ATP ⇌ EATP
            EATP <--> EF16BP_res_plusATP_minusF16BP
            E_res_plusATP_minusF16BP + F16BP ⇌ EF16BP_res_plusATP_minusF16BP
            E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP
            EF6P_res_plusATP_minusF16BP <--> EADP
            E + ADP ⇌ EADP
        end
    end
    nfp = parameters(nonallo)
    half2_rate = ER.compile_rate_equation(half2)
    nonallo_rate = ER.compile_rate_equation(nonallo)
    twin_of(s) = Symbol(replace(String(s), r"^([Kk])_A_" => s"\1_"))
    @test Set(twin_of.(filter(!=(:L), collect(afp)))) == Set(nfp)
    rng = MersenneTwister(20260717)
    for _ in 1:6
        vals = Dict(s => 0.5 + 2rand(rng) for s in nfp)
        Keq = 2.0 + 2rand(rng)
        ATP = 0.5+2rand(rng); F6P = 0.5+2rand(rng)
        ADP = 0.5+2rand(rng); F16BP = 0.5+2rand(rng)
        concs = (ATP=ATP, F6P=F6P, ADP=ADP, F16BP=F16BP)
        prm = NamedTuple{afp}(Tuple(s === :L ? 0.7 : vals[twin_of(s)] for s in afp))
        @test isfinite(real(half2_rate(concs, prm, Keq)))
        @test isfinite(real(half2_rate(
            (ATP=ATP, F6P=F6P, ADP=ADP, F16BP=1e-8), prm, Keq)))
        # With both products at zero the inactive conformation holds its enzyme in the
        # covalent forms (half 2 cannot run, half 1 runs back only with F16BP): rate 0.
        @test real(half2_rate(merge(concs, (ADP=0.0, F16BP=0.0)), prm, Keq)) == 0
        @test real(half2_rate(
            (ATP=10.0, F6P=10.0, ADP=1e-3, F16BP=1e-3), prm, Keq)) > 0
        # v = 0 at the equilibrium metabolite ratio ADP·F16BP = Keq·ATP·F6P.
        @test isapprox(real(half2_rate(
            (ATP=ATP, F6P=F6P, ADP=Keq*ATP*F6P/F16BP, F16BP=F16BP), prm, Keq)),
            0.0; atol=1e-8)
        # L = 0 cross-check: the allosteric rate reduces to the non-allosteric twin.
        pn = NamedTuple{nfp}(Tuple(vals[s] for s in nfp))
        p0 = NamedTuple{afp}(Tuple(s === :L ? 0.0 : vals[twin_of(s)] for s in afp))
        @test isapprox(real(half2_rate(concs, p0, Keq)),
                       real(nonallo_rate(concs, pn, Keq)); rtol=1e-4)
    end
end
