# ABOUTME: Prints a mechanism's steps, its IndependentParameters law, dependents and Hill quotient.
# ABOUTME: Usage: julia peek_law.jl <population/index>...; throwaway P4 helper.
include(joinpath(@__DIR__, "P4.jl"))
using .Poly: Mono, mono, mmul, mpow, mexp, P, nterms
for id in ARGS
    pop, i = split(id, "/")
    m = Pop.load_population(pop)[parse(Int, i)]
    println("== ", m.id, " ", m.source, " ", m.family)
    for (j, s) in enumerate(m.steps)
        println("  $j ", Pop.step_text(m, s), "  ", m.step_names[j])
    end
    println("  groups ", m.groups)
    L = Engine.bare_rate_law(m; check = false)
    dpre, shift = P4.den_unshifted(L)
    println("  shift ", shift, " criteria ", P4.criteria(m))
    NL = P4.law(m)
    for d in NL.sol.deps
        println("  dep ", NL.sol.cols[d.col].name, " = Keq^", d.keq, " * ",
                join([string(NL.sol.cols[f].name, "^", e) for (f, e) in d.terms], " * "), "  ", d.label,
                " class ", P4.depclass(NL.sol, d))
    end
    W = NL.vars
    D = P4.cmap(NL, NL.den); N = P4.cmap(NL, NL.num)
    println("  D:")
    for (mo, c) in sort(collect(D); by = x -> (P4.degree(x[1]), P4.fmt_mono(x[1], W)))
        println("    [", P4.fmt_mono(mo, W), "] ", P4.fmt_poly(c, W))
    end
    println("  N:")
    for (mo, c) in sort(collect(N); by = x -> (P4.degree(x[1]), P4.fmt_mono(x[1], W)))
        println("    [", P4.fmt_mono(mo, W), "] ", P4.fmt_poly(c, W))
    end
    h = P4.hill(NL)
    println("  Q = ", h.Q === nothing ? "none" : P4.fmt_poly(h.Q, W), "  keq_free=", h.keq_free)
end
