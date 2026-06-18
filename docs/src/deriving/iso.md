# Iso mechanisms

In an iso mechanism the enzyme passes through two distinct conformations during
a single turnover. The substrate binds one conformation; catalysis converts the
enzyme–substrate complex into the *other* conformation of the enzyme–product
complex; the product leaves from that conformation; and the free enzyme must
then isomerize back to the first conformation before it can bind substrate
again. That free-enzyme isomerization — the "iso" step — is the defining
feature.

## An iso uni-uni example

A second enzyme conformation is written with a different **conformation label** —
`Eprime` in place of `E`. Bound forms carry their conformation (`Eprime(P)`),
and the free enzyme moves between conformations through an isomerization step
such as `Eprime <--> E`. Here `S` binds `E`, catalysis turns `E(S)` into
`Eprime(P)`, `P` leaves `Eprime`, and the free enzyme isomerizes from `Eprime`
back to `E`:

```@example iso
using EnzymeRates
m = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> Eprime(P)
        Eprime(P) ⇌ Eprime + P
        Eprime <--> E
    end
end
print(rate_equation_string(m))
```

The derived equation carries the signature of an iso mechanism: an `S·P` cross
term in the denominator (the `… P * S / (K_P_Eprime * K_S_E)` terms above). An
ordinary single-conformation Michaelis–Menten enzyme does not have an `S·P`
term. Saturating an iso enzyme with substrate and product together can
therefore slow turnover in a way a regular Michaelis–Menten enzyme cannot. Each
conformation appears in the parameter
names: `K_S_E` is `S` binding to `E`, `K_P_Eprime` is `P` binding to `Eprime`,
and `k_E_to_Eprime` / `k_Eprime_to_E` are the isomerization rate constants.
