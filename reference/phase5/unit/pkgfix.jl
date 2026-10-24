# ABOUTME: Checks whether the package (census env) accepts the unit-rule fixtures written in its
# ABOUTME: DSL: the ratio-form and fractional-λ two-cycle laws and the RE Theorell–Chance law.
using EnzymeRates

fixtures = [
    "F1 two-cycle uni-uni" => :(@enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
        end
    end),
    "F2 F1 + F + A + A ⇌ F(A, A)" => :(@enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A + A ⇌ F(A, A)
        end
    end),
    "F3 RE Theorell–Chance" => :(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(Q) + P
            E(Q) <--> E + Q
        end
    end),
    "F4 F1 + F(A) + F(A, A)" => :(@enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A ⇌ F(A)
            F(A) + A ⇌ F(A, A)
        end
    end),
]
for (name, ex) in fixtures
    println("=== ", name)
    try
        m = eval(ex)
        println(rate_equation_string(m))
    catch err
        println("REJECTED: ", first(sprint(showerror, err), 600))
    end
end
