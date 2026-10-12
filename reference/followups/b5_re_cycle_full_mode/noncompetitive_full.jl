# ABOUTME: Prints the Full-mode law of the Non-competitive Inhibitor spec, whose open RE cycle
# ABOUTME: the single-symbol Wegscheider fold closes inside the Full polynomial.
using EnzymeRates
m = @enzyme_mechanism begin
    substrates: S
    products: P
    regulators: R
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E + P ⇌ E(P)
        (E + R ⇌ E(R), E(S) + R ⇌ E(S, R))
        E(R) + S ⇌ E(S, R)
    end
end
println("parameters(m, Full) = ", parameters(m, Full))
println(rate_equation_string(m, Full))
