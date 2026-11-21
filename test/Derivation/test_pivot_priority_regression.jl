# ABOUTME: Regression: a Wegscheider tie between a free-enzyme binding K and an occupied
# ABOUTME: one makes one of them dependent, so parameters(m) counts the tie once.
module PivotPriorityRegressionTests
using Test, EnzymeRates
const ER = EnzymeRates

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
