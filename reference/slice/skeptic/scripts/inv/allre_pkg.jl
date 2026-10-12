# ABOUTME: What today's package does with the all-rapid-equilibrium uni-uni mechanism and with
# ABOUTME: metabolites named N, E_total and A_B (the printed-name collisions found in the slice).
using EnzymeRates
const ER = EnzymeRates
function try_mech(label, src)
    print(rpad(label, 22))
    m = try
        Core.eval(Main, Meta.parse(src))
    catch e
        println("rejected: ", first(sprint(showerror, e), 220)); return
    end
    try
        em = m isa ER.EnzymeMechanism ? m : ER.compile_mechanism(m)
        mets = ER.metabolites(em)
        ps = (ER.fitted_params(em)..., :Keq, :E_total)
        v = ER.rate_equation(em, NamedTuple{Tuple(mets)}(Tuple(1.5 for _ in mets)),
                             NamedTuple{ps}(Tuple(1.5 for _ in ps)))
        println("derived; rate at all 1.5 = ", v)
    catch e
        println("derive error: ", first(sprint(showerror, e), 220))
    end
end
try_mech("all_re", """@enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)
        E(S) ⇌ E(P)
        E(P) ⇌ E + P
    end
end""")
for name in ("N", "E_total", "Keq")
    try_mech("inhibitor $name", """@enzyme_mechanism begin
    substrates: S
    products: P
    regulators: $name
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E + P
        E + $name ⇌ E($name)
    end
end""")
end
