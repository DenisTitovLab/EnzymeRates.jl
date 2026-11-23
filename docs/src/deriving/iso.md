# Iso mechanisms

In an iso mechanism the enzyme passes through two distinct conformations during
a single turnover. The substrate binds one conformation; catalysis converts the
enzyme–substrate complex into the *other* conformation of the enzyme–product
complex; the product leaves from that conformation; and the free enzyme must
then isomerize back to the first conformation before it can bind substrate
again. That free-enzyme isomerization — the "iso" step — is the defining
feature.

## An iso uni-uni example

A second enzyme conformation is a second free form, written by its name: `Eprime`
beside `E`. Its bound forms take names of their own (`EprimeP`), and the free enzyme
moves between conformations through an isomerization step such as `Eprime <--> E`. Here
`S` binds `E`, catalysis turns `ES` into `EprimeP`, `P` leaves `Eprime`, and the free
enzyme isomerizes from `Eprime` back to `E`:

```@example iso
using EnzymeRates
m = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S ⇌ ES
        ES <--> EprimeP
        EprimeP ⇌ Eprime + P
        Eprime <--> E
    end
end
print(rate_equation_string(m))
```

The derived equation carries the signature of an iso mechanism: an `S·P` cross
term in the denominator, the term of `D` above that holds both `S` and `P`.
An ordinary single-conformation Michaelis–Menten enzyme does not
have an `S·P` term. Saturating an iso enzyme with substrate and product
together can therefore slow turnover in a way a regular Michaelis–Menten
enzyme cannot. The forms' names make the parameter names: `K_ES_to_E_S`
is `S` binding to `E`, `K_EprimeP_to_Eprime_P` is `P` binding to `Eprime`, and
`k_E_to_Eprime` / `k_Eprime_to_E` are the isomerization rate constants.
