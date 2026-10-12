# ABOUTME: Skeptic of probe P2: Haldane-multiple gcd on the adversarial set, and a group whose two
# ABOUTME: SS constant classes take their representatives from different members.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/adv/AdvLib.jl")
using Serialization
const K = joinpath(PROTO, "skeptic_P2-solve")
for m in [deserialize(joinpath(K, "data", "adv_pkg.jls")); deserialize(joinpath(K, "data", "adv_pkg2.jls"))]
    sol = Thermo.solve(m)
    sol.haldane_gcd == 1 || println(rpad(m.source, 28), " gcd of Haldane cycle multiples = ", sol.haldane_gcd)
end

M(n, r) = Pop.PMet(n, r)
S_ = M(:S, :S); P_ = M(:P, :P); I_ = M(:I, :Inh); J_ = M(:J, :Inh); NONE = Pop.PMet[]
sp(conf, b...) = Pop.PSpecies(collect(b), conf, Symbol[], Symbol[],
                              Symbol(conf, join(string.(getfield.(collect(b), :name)))))
species = [sp(:E), sp(:E, S_), sp(:E, P_), sp(:F), sp(:E, J_), sp(:F, J_, S_), sp(:F, I_), sp(:E, I_, S_)]
steps = [Pop.PStep(1, 2, [S_], NONE, true), Pop.PStep(2, 3, NONE, NONE, false),
         Pop.PStep(1, 3, [P_], NONE, true), Pop.PStep(1, 4, NONE, NONE, true),
         Pop.PStep(1, 5, [J_], NONE, true), Pop.PStep(4, 7, [I_], NONE, true),
         Pop.PStep(5, 6, [S_], NONE, false),            # E(J) + S <--> F(J, S)   (group)
         Pop.PStep(7, 8, [S_], NONE, false)]            # F(I) + S <--> E(I, S)   (group)
m = Pop.PMech(; id = "hand:crossed_reps", source = "hand:crossed_reps", substrates = [:S],
              products = [:P], inhibitors = [:I, :J], metabolites = [:S, :P, :I, :J],
              species, steps, groups = [[1], [2], [3], [4], [5], [6], [7, 8]])
sol = Thermo.solve(m)
println("\n== crossed representatives: group {E(J) + S <--> F(J,S), F(I) + S <--> E(I,S)}")
println("   independents: ", Thermo.independent_names(sol))
for d in sol.deps
    println("   dep ", Thermo.colname(sol, d.col), " = ", join([string(Thermo.colname(sol, f), "^", e) for (f, e) in d.terms], " * "),
            d.keq == 0 ? "" : " * Keq^$(d.keq)", "  [", d.label, d.fold ? ", fold" : "", "]")
end
