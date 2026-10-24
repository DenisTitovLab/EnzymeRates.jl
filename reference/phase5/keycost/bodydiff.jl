# ABOUTME: Prints the generated rate_equation body of two specs whose body hash differs between
# ABOUTME: Julia 1.12 and 1.13, to locate the order that depends on the Julia version.
include(joinpath(@__DIR__, "klib.jl"))
for r in tsv(joinpath(@__DIR__, "data", "x5_small.tsv"))
    r[3] in ("62569106e04b1827", "446a2c9db9b9c356") || continue
    m = parse_neutral(r[end], names_of(r[4]), names_of(r[5]))
    sol = S.solve(m)
    p = S.program(m, sol, S.rate_law(m, sol))
    println("# ", r[2], " (Julia ", VERSION, ")")
    for (v, q) in p.defs
        println("def ", p.vars.names[v], " = ", sort!([string(mo, "=>", c) for (mo, c) in q.terms]))
    end
    println(S.body(p))
end
