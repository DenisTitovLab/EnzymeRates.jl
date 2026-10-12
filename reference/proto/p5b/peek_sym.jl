# ABOUTME: Prints the A and I polynomials of mechanisms whose tag-stripped polynomials match but whose
# ABOUTME: printed text differs (probe P5b symmetric-printing check).
include(joinpath(@__DIR__, "Allo.jl"))
pop = filter(m -> m.allosteric, Pop.load_population("handwritten"))
for m in pop
    md = Allo.derive(m)
    isconc(s) = s in m.metabolites
    for k in 1:2
        strip(p) = (S = Allo.Vars(); (Allo.Poly.symbolic(Allo.subst(p, v -> Allo.mono(Allo.var!(S, Symbol(Allo.rolekey(md.T.names[v])[1])) => 1)), S)))
        strip(md.A[k]) == strip(md.I[k]) || continue
        tA, tI = (replace(Allo.poly_text(p, md.T, isconc), r"\b([kK])_[AI]_" => s"\1_") for p in (md.A[k], md.I[k]))
        tA == tI && continue
        println("== $(m.id) poly $k\n A: ", Allo.poly_text(md.A[k], md.T, isconc), "\n I: ", Allo.poly_text(md.I[k], md.T, isconc))
    end
end
