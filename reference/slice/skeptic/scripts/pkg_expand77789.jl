# ABOUTME: Reproduces the package's expand_mechanisms failure on init_terter/77789 and prints the
# ABOUTME: stack frames inside the package, to see which move built the rejected mechanism.
using EnzymeRates
const ER = EnzymeRates
const TESTDIR = joinpath(dirname(dirname(pathof(EnzymeRates))), "test")
ast = Meta.parseall(read(joinpath(TESTDIR, "test_mechanism_enumeration.jl"), String))
for ex in ast.args
    ex isa Expr && ex.head === :const && ex.args[1] isa Expr &&
        ex.args[1].args[1] === :ter_ter_rxn && Core.eval(Main, ex)
end
init = unique(ER.init_mechanisms(ter_ter_rxn))
m = init[77789]
println(m)
try
    ER.expand_mechanisms([m], ter_ter_rxn)
    println("no error")
catch e
    println("ERROR: ", first(sprint(showerror, e), 400))
    for fr in stacktrace(catch_backtrace())
        s = string(fr)
        occursin("EnzymeRates", s) && println("  ", first(s, 200))
    end
end
