# ABOUTME: Evaluates the Full- and Reduced-mode rate of one random-order RE square under two
# ABOUTME: metabolite namings at the same physical constants, with and without the RE cycle relation.
using EnzymeRates

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
rrb = @enzyme_mechanism begin     # spec "RE Random Bi-Bi": inner bindings are SS
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)
        E + B ⇌ E(B)
        E(A) + B <--> E(A, B)
        E(B) + A <--> E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) <--> E(Q) + P
        E(Q) <--> E + Q
    end
end
function absent_from_v(m)
    vline = last(split(rate_equation_string(m, Full), "\n"))
    toks = Set(Symbol(x.match) for x in eachmatch(r"[A-Za-z_][A-Za-z0-9_+\-]*", vline))
    [p for p in parameters(m, Full) if p ∉ toks]
end
println("RE Random Bi-Bi spec: Full params absent from v = ", absent_from_v(rrb),
        "  (", length(parameters(rrb, Full)) - 1, " listed constants)")

# Physical constants: KA_E = K(A|E), KB_E = K(B|E), KB_EA = K(B|EA), KA_EB = K(A|EB).
fullA(KA_EB) = (K_EA_to_E_A = 1.0, K_EB_to_E_B = 2.0, K_EP_to_E_P = 0.5, K_EAB_to_EA_B = 3.0,
                k_EAB_to_EPQ = 10.0, k_EPQ_to_EAB = 1.0, K_EAB_to_EB_A = KA_EB,
                K_EPQ_to_EP_Q = 0.7, E_total = 1.0)
fullZ(KA_EB) = (K_EB_to_E_B = 2.0, K_EP_to_E_P = 0.5, K_EZ_to_E_Z = 1.0, K_EBZ_to_EB_Z = KA_EB,
                k_EBZ_to_EPQ = 10.0, k_EPQ_to_EBZ = 1.0, K_EPQ_to_EP_Q = 0.7,
                K_EBZ_to_EZ_B = 3.0, E_total = 1.0)
cA = (A = 2.0, B = 1.0, P = 0.3, Q = 0.2)
cZ = (B = 1.0, Z = 2.0, P = 0.3, Q = 0.2)
@assert keys(fullA(1.0)) == parameters(mA, Full) && keys(fullZ(1.0)) == parameters(mZ, Full)

println("\nFull mode, K(A|EB) inconsistent with K(A|E)K(B|EA)/K(B|E) = 1.5:")
for KA_EB in (4.0, 1e6)
    println("  K(A|EB) = ", KA_EB, ":  v_Full(A-named) = ", rate_equation(mA, cA, fullA(KA_EB), Full),
            "   v_Full(Z-named) = ", rate_equation(mZ, cZ, fullZ(KA_EB), Full))
end
println("Full mode, K(B|EA) varied with K(A|EB) = 4 (the Z-named law ignores K(B|EA)):")
for KB_EA in (3.0, 1e6)
    pa = merge(fullA(4.0), (K_EAB_to_EA_B = KB_EA,)); pz = merge(fullZ(4.0), (K_EBZ_to_EZ_B = KB_EA,))
    println("  K(B|EA) = ", KB_EA, ":  v_Full(A-named) = ", rate_equation(mA, cA, pa, Full),
            "   v_Full(Z-named) = ", rate_equation(mZ, cZ, pz, Full))
end

# Consistent constants: K(A|EB) = 1*3/2 = 1.5, Keq from the Haldane relation with kr = 1.
Keq = 10.0 * 0.7 * 0.5 / (1.0 * 1.0 * 3.0)
redA = (K_EAB_to_EB_A = 1.5, K_EA_to_E_A = 1.0, K_EB_to_E_B = 2.0, K_EPQ_to_EP_Q = 0.7,
        K_EP_to_E_P = 0.5, k_EAB_to_EPQ = 10.0, Keq = Keq, E_total = 1.0)
redZ = (K_EBZ_to_EZ_B = 3.0, K_EB_to_E_B = 2.0, K_EPQ_to_EP_Q = 0.7, K_EP_to_E_P = 0.5,
        K_EZ_to_E_Z = 1.0, k_EBZ_to_EPQ = 10.0, Keq = Keq, E_total = 1.0)
@assert keys(redA) == parameters(mA) && keys(redZ) == parameters(mZ)
println("\nConsistent constants (K(A|EB) = 1.5, Keq = ", Keq, "):")
println("  v_Full(A-named)    = ", rate_equation(mA, cA, fullA(1.5), Full))
println("  v_Full(Z-named)    = ", rate_equation(mZ, cZ, fullZ(1.5), Full))
println("  v_Reduced(A-named) = ", rate_equation(mA, cA, redA))
println("  v_Reduced(Z-named) = ", rate_equation(mZ, cZ, redZ))
