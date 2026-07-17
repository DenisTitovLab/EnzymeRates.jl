# ABOUTME: Inspect per-conformation free-enzyme weights + I-state polys for both mechs.
using EnzymeRates
const ER = EnzymeRates
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/fb96a560-d2b8-4694-a474-ff2f47968c52/scratchpad"
load(p) = Core.eval(ER, Meta.parse(strip(read(p,String))))()

function inspect(tag, m)
    am = ER.AllostericMechanism(m)
    numA,denA,dfA = ER._state_rate_polys(am, :A)
    numI,denI,dfI = ER._state_rate_polys(am, :I)
    println("\n===== $tag =====")
    println("  d_free_A = ", ER._poly_to_expr(dfA, Set{Symbol}(), Set{Symbol}()))
    println("  d_free_I = ", ER._poly_to_expr(dfI, Set{Symbol}(), Set{Symbol}()))
    println("  d_free_A == d_free_I ? ", dfA == dfI)
    println("  num_I is zero (dead I-cycle)? ", numI == ER.poly_zero())
    println("  num_A terms=", length(numA), "  den_A terms=", length(denA),
            "  num_I terms=", length(numI), "  den_I terms=", length(denI))
end

inspect("ERR1 (g3 :OnlyA)", load("$SP/pfkp_err1.txt"))
inspect("ERR2 (g4 :OnlyA)", load("$SP/pfkp_err2.txt"))
