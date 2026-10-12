# ABOUTME: Smoke test of Thermo.solve: every plain mechanism of every population, counts vs the package.
# ABOUTME: Throwaway P2 probe script; stdlib only.
include(joinpath(@__DIR__, "..", "common", "Thermo.jl"))
using Printf
# checked Rational{Int64}?
try
    x = (typemax(Int64) - 1) // 1 + 3 // 1
    println("UNCHECKED rational add: ", x)
catch e
    println("rational add overflow -> ", typeof(e))
end
tot = 0; diff = 0; errs = 0
for pop in Pop.populations()
    ms = Pop.load_population(pop)
    for m in ms
        m.allosteric && continue
        global tot += 1
        try
            sol = Thermo.solve(m)
            n = length(sol.indep)
            if n != length(m.fitted_params) || n != m.independent_param_count
                global diff += 1
                println("DIFF ", m.id, " ", m.source, ": mine=", n, " fitted=", length(m.fitted_params), " ipc=", m.independent_param_count)
            end
        catch e
            global errs += 1
            println("ERR ", m.id, " ", m.source, ": ", sprint(showerror, e)[1:min(end, 300)])
        end
    end
end
println("plain total=$tot diff=$diff errs=$errs")
