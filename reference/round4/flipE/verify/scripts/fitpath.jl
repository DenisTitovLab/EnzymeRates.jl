# ABOUTME: What the fitting path does with a newly admitted law that vanishes at zero products:
# ABOUTME: its loss on product-free data, and the post-fit kcat rescale (rescale_parameter_values).
using EnzymeRates, Random, Printf
const ER = EnzymeRates

m = @allosteric_mechanism begin
    substrates: S
    products: P
    catalytic_multiplicity: 1
    catalytic_steps: begin
        E + S <--> E(S) :: EqualAI
        E(S) <--> Eprime(P) :: EqualAI
        Eprime + P <--> Eprime(P) :: EqualAI
        Eprime <--> E :: OnlyA
    end
end
println("package: ", pathof(EnzymeRates))
println("constructed; fitted_params = ", ER.fitted_params(m))
# Michaelis–Menten-like initial-rate data at zero product plus a product series.
S = [0.1, 0.3, 1.0, 3.0, 10.0, 1.0, 1.0, 1.0]
P = [0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 1.0, 3.0]
Rate = [1.0 * s / (1.0 + s + p / 2) for (s, p) in zip(S, P)]
table = (S = S, P = P, Rate = Rate, group = fill(1, length(S)))
fp = FittingProblem(m, table; Keq = 10.0)
x = zeros(length(ER.fitted_params(m)))
@printf("loss at log-params 0: %.4g (8 points, 5 at P = 0)\n", ER.loss!(x, fp))
prm = merge(NamedTuple{ER.fitted_params(m)}(Tuple(ones(length(x)))), (Keq = 10.0, E_total = 1.0))
for (s, p) in zip(S, P)
    @printf("  S=%-5g P=%-4g rate %.4g\n", s, p, rate_equation(m, (S = s, P = p), prm))
end
try
    rescale_parameter_values(m, prm)
    println("rescale_parameter_values: returned")
catch err
    println("rescale_parameter_values threw: ", first(sprint(showerror, err), 200))
end
