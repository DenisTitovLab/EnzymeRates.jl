# ABOUTME: Probe the shipped allosteric ping-pong :NonequalAI mechanism: names + rate eq string.
using EnzymeRates
const ER = EnzymeRates

m = @allosteric_mechanism begin
    substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
    catalytic_steps: begin
        E + A <--> E(A)                                       :: EqualAI
        E(A) <--> E(; residual = A - P) + P                   :: NonequalAI
        E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
        E(B; residual = A - P) <--> E + Q                     :: NonequalAI
    end
end

println("=== fitted_params ===")
println(ER.fitted_params(m))
println("\n=== rate_equation_string (Reduced) ===")
println(ER.rate_equation_string(m))
println("\n=== rate_equation_string (Full) ===")
println(ER.rate_equation_string(m, ER.Full))
