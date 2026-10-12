# ABOUTME: Shows that Full mode drops the K closing an open RE cycle and that the dropped K
# ABOUTME: depends on metabolite names (canonical step order), while Reduced mode does not.
using EnzymeRates
const ER = EnzymeRates

# expand_bibi/265: random-order RE binding square, SS chemistry, ordered RE product release.
mA = @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)
        E + B ⇌ E(B)
        E(A) + B ⇌ E(A, B)
        E(B) + A ⇌ E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) ⇌ E(P) + Q
        E(P) ⇌ E + P
    end
end
# Same mechanism, steps written in another order.
mA_shuffled = @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E(P) ⇌ E + P
        E(B) + A ⇌ E(A, B)
        E(A, B) <--> E(P, Q)
        E + B ⇌ E(B)
        E(P, Q) ⇌ E(P) + Q
        E(A) + B ⇌ E(A, B)
        E + A ⇌ E(A)
    end
end
# Same chemistry, substrate A renamed Z.
mZ = @enzyme_mechanism begin
    substrates: Z, B
    products: P, Q
    steps: begin
        E + Z ⇌ E(Z)
        E + B ⇌ E(B)
        E(Z) + B ⇌ E(B, Z)
        E(B) + Z ⇌ E(B, Z)
        E(B, Z) <--> E(P, Q)
        E(P, Q) ⇌ E(P) + Q
        E(P) ⇌ E + P
    end
end

function vline_missing(m)
    s = rate_equation_string(m, Full)
    vline = last(split(s, "\n"))
    toks = Set(Symbol.(eachmatch(r"[A-Za-z_][A-Za-z0-9_+\-]*", vline) .|> x -> x.match))
    [p for p in parameters(m, Full) if p ∉ toks]
end

println("mA === mA_shuffled: ", mA === mA_shuffled,
        "   Mechanism equal: ", ER.Mechanism(mA) == ER.Mechanism(mA_shuffled))
for (lbl, m) in (("A-named", mA), ("Z-named", mZ))
    println("\n================ ", lbl)
    println("stored flat step order:")
    for (j, (s, g)) in enumerate(ER._flat_steps(ER.Mechanism(m)))
        println("  [", j, "] g", g, "  ", s)
    end
    println("parameters(m, Full) = ", parameters(m, Full))
    println("--- rate_equation_string(m, Full):")
    println(rate_equation_string(m, Full))
    println("listed in Full but absent from the v line: ", vline_missing(m))
    println("--- parameters(m) (Reduced) = ", parameters(m))
    println("--- rate_equation_string(m) (Reduced):")
    println(rate_equation_string(m))
end
