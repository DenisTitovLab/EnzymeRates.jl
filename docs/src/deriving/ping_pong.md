# Ping-pong mechanisms

In a ping-pong (double-displacement) mechanism the enzyme carries a covalent
fragment between two half-reactions. The first substrate binds and transfers part
of its atoms to the enzyme; the modified enzyme releases the first product, then
binds the second substrate and completes the transfer to form the second product.
The enzyme returns to its original form only at the end of the full cycle.

## A ping-pong bi-bi example

A covalently modified enzyme form is written by giving the enzyme a `residual = …`
keyword — an arithmetic expression over the declared metabolites that lists the
atoms it has **gained** (with `+`) and is **committed to release** (with `-`).
Here substrate `A` binds, product `P` leaves a modified enzyme
`E(; residual = A - P)` carrying `A`'s atoms, substrate `B` binds the modified
form, and product `Q` leaves to regenerate free `E`:

```@example pingpong
using EnzymeRates
m = @enzyme_mechanism begin
    substrates: A, B
    products:   P, Q
    steps: begin
        E + A <--> E(A)
        E(A) <--> E(; residual = A - P) + P
        E(; residual = A - P) + B <--> E(B; residual = A - P)
        E(B; residual = A - P) <--> E + Q
    end
end
print(rate_equation_string(m))
```

The derived equation has the signature of a ping-pong mechanism: its denominator has **no
constant term**. A sequential mechanism always carries a `1` in its denominator for the
free enzyme; a ping-pong enzyme never sits idle as free `E` during turnover (it is always
either loaded with substrate or carrying the covalent residual), so that term is absent.

Without a constant term, this law divides by Cleland's normalizer ([Normalization](@ref)).
Its denominator holds a term in `B` and a term in `A` (the product of the substrates over
each substrate) and a term in `A * B` (their product), so the normalizer is the product of
the coefficients of `B` and `A` divided by the coefficient of `A * B`. The substrate terms
then read `A/K_a`, `B/K_b` and `A·B/(K_a·K_b)`, `K_a` and `K_b` the Michaelis constants. The
normalizer's part that is not a single term prints once as `N = …`, and every denominator
term is dimensionless. A denominator that lacks one of these terms takes the next rule that
applies ([Normalization](@ref)): a virtual constant term on the first affinely independent
terms of a fixed order, which changes how `N` reads but never the law's value, or, when no
weights give the denominator's unit, the ratio form.

The coefficient form names each coefficient by the concentrations it multiplies:

```@example pingpong
print(rate_equation_string(m; coefficients = true))
```

The residual-bearing forms take descriptive names such as `E_res_plusA_minusP` (the
modified free enzyme) and `EB_res_plusA_minusP` (with `B` additionally bound).
