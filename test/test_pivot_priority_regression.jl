# ABOUTME: Regression for the shared -1 pivot-priority sentinel — a Wegscheider row
# ABOUTME: whose only pivot is a free-enzyme binding K must still be pivoted, not dropped.
module PivotPriorityRegressionTests
using Test, EnzymeRates
const ER = EnzymeRates

@testset "pivot priority: non-allo Wegscheider tie enforced (no over-count)" begin
    # A non-allosteric LDH mechanism reached in enumeration. A closed all-RE cycle
    # — E→ELactate→ELactateNAD and E→ENAD→ELactateNAD with shared K_NAD_E — forces
    # the Wegscheider tie K_Lactate_E = K_Lactate_ENADH. The tie's only pivot column
    # is a free-enzyme binding K (priority -1); the sentinel never selected it, so
    # the row was dropped as "redundant" and fitted_params over-counted by 1
    # (7 on the buggy kernel, 6 correct). This fixture is the O(1) distillation of
    # the split-monotonicity property (a split that recovers the dropped tie is a
    # -1 param-count edge, impossible on a correct kernel).
    m = Core.eval(EnzymeRates, Meta.parse(raw"EnzymeMechanism{(((((:Product, :Lactate), ((:C, 3), (:H, 6), (:O, 3))), ((:Product, :NAD), ((:C, 21), (:H, 27), (:N, 7), (:O, 14), (:P, 2))), ((:Substrate, :NADH), ((:C, 21), (:H, 29), (:N, 7), (:O, 14), (:P, 2))), ((:Substrate, :Pyruvate), ((:C, 3), (:H, 4), (:O, 3)))), (), (4,)), (((((), :E, ((), ())), (((:Product, :Lactate),), :E, ((), ())), (:Product, :Lactate), true),), ((((), :E, ((), ())), (((:Product, :NAD),), :E, ((), ())), (:Product, :NAD), true), ((((:Product, :Lactate),), :E, ((), ())), (((:Product, :Lactate), (:Product, :NAD)), :E, ((), ())), (:Product, :NAD), true)), ((((), :E, ((), ())), (((:Substrate, :NADH),), :E, ((), ())), (:Substrate, :NADH), false), ((((:Product, :Lactate),), :E, ((), ())), (((:Product, :Lactate), (:Substrate, :NADH)), :E, ((), ())), (:Substrate, :NADH), false)), (((((:Product, :NAD),), :E, ((), ())), (((:Product, :Lactate), (:Product, :NAD)), :E, ((), ())), (:Product, :Lactate), true), ((((:Substrate, :NADH),), :E, ((), ())), (((:Product, :Lactate), (:Substrate, :NADH)), :E, ((), ())), (:Product, :Lactate), true)), (((((:Product, :NAD),), :E, ((), ())), (((:Product, :NAD), (:Substrate, :Pyruvate)), :E, ((), ())), (:Substrate, :Pyruvate), true), ((((:Substrate, :NADH),), :E, ((), ())), (((:Substrate, :NADH), (:Substrate, :Pyruvate)), :E, ((), ())), (:Substrate, :Pyruvate), true)), (((((:Substrate, :NADH), (:Substrate, :Pyruvate)), :E, ((), ())), (((:Product, :Lactate), (:Product, :NAD)), :E, ((), ())), nothing, false),)))}"))()
    fp = ER.fitted_params(m)
    @test length(fp) == 6                                   # was 7 (over-counted)
    @test !(:K_Lactate_E in fp && :K_Lactate_ENADH in fp)   # the tie collapses one
end
end # module
