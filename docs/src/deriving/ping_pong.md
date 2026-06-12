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

The derived equation has the signature of a ping-pong mechanism: the numerator is
the difference `A·B − P·Q`, and the denominator has **no constant term**. A
sequential mechanism always carries a `1` in its denominator for the free enzyme;
a ping-pong enzyme never sits idle as free `E` during turnover — it is always
either loaded with substrate or carrying the covalent residual — so that constant
term is absent. The residual-bearing forms take descriptive names such as
`E_res_+A_-P` (the modified free enzyme) and `EB_res_+A_-P` (with `B`
additionally bound).
