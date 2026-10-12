# ABOUTME: Debugs the hand-built exchange-in-Wegscheider-cycle case: solve and per-parameter derivatives.
# ABOUTME: Throwaway P2 probe script; stdlib only.
include(joinpath(@__DIR__, "..", "common", "Oracle.jl"))
include(joinpath(@__DIR__, "..", "common", "Thermo.jl"))
using Random, Printf
const T = Thermo
setprecision(BigFloat, 512)
sp(conf, bound...) = Pop.PSpecies(collect(bound), conf, Symbol[], Symbol[],
                                  Symbol(conf, join(string.(getfield.(collect(bound), :name)))))
S = Pop.PMet(:S, :S); P = Pop.PMet(:P, :P); I = Pop.PMet(:I, :Inh); J = Pop.PMet(:J, :Inh)
species = [sp(:E), sp(:E, S), sp(:E, P), sp(:E, I), sp(:E, J)]
steps = [Pop.PStep(1, 2, [S], Pop.PMet[], true), Pop.PStep(2, 3, Pop.PMet[], Pop.PMet[], false),
         Pop.PStep(1, 3, [P], Pop.PMet[], true), Pop.PStep(1, 4, [I], Pop.PMet[], true),
         Pop.PStep(4, 5, [J], [I], false), Pop.PStep(1, 5, [J], Pop.PMet[], true)]
m = Pop.PMech(; id = "w", source = "w", substrates = [:S], products = [:P], inhibitors = [:I, :J],
              metabolites = [:S, :P, :I, :J], species = species, steps = steps,
              groups = [[j] for j in eachindex(steps)])
sol = T.solve(m)
for (i, c) in enumerate(sol.cols); println(i, " ", c.prio, " ", c.name, i in sol.indep ? " [indep]" : ""); end
for d in sol.deps
    println("dep ", T.colname(sol, d.col), " = Keq^", d.keq, " ", join([string(T.colname(sol, f), "^", e) for (f, e) in d.terms], " "), " ", d.label)
end
println("zero flux: ", T.zero_flux_ss(m))
rng = Xoshiro(5)
keq = Oracle.logu(rng); θ = Dict(i => Oracle.logu(rng) for i in sol.indep)
c = Dict{Symbol, BigFloat}(x => Oracle.logu(rng) for x in m.metabolites)
h = big(10)^-30
for i in sol.indep
    vs = map((1, -1)) do sg
        θp = copy(θ); θp[i] *= exp(sg * h)
        Oracle.steady_state_v(m, T.stored_constants(sol, T.all_values(sol, θp, keq)), c)
    end
    @printf("d v / d log %-32s = %.3e\n", T.colname(sol, i), Float64((vs[1] - vs[2]) / (2h)))
end
k = T.stored_constants(sol, T.all_values(sol, θ, keq))
println("stored constants: ", [Float64.(x) for x in k])
x = Oracle.steady_state_populations(m, k, c)
println("populations: ", Float64.(x))
# flux through the exchange step
f5 = k[5][1] * c[:J] * x[4] - k[5][2] * c[:I] * x[5]
println("net flux exchange: ", Float64(f5), "  v: ", Float64(Oracle.steady_state_v(m, k, c)))
