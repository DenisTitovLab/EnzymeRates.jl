# Ping-pong mechanisms

In a ping-pong (double-displacement) mechanism the enzyme carries a covalent
fragment between two half-reactions. The first substrate binds and transfers part
of its atoms to the enzyme; the modified enzyme releases the first product, then
binds the second substrate and completes the transfer to form the second product.
The enzyme returns to its original form only at the end of the full cycle.

## Writing the covalent intermediate

A covalently modified enzyme form is written with a conformation label — any name
other than `E` — carrying a `residual = …` keyword. The residual lists the atoms
the enzyme has **gained** (with `+`) and the atoms it is **committed to release**
(with `-`), as an arithmetic expression over the declared metabolites:

```julia
Estar(; residual = A - P)    # free modified enzyme: carries A's atoms, owes P
Estar(B; residual = A - P)   # the same intermediate, with substrate B bound
```

Any label works (`Estar`, `F`, …); it is the `residual` payload, not the label,
that marks the form as a covalent intermediate.

## A ping-pong bi-bi example

Substrate `A` binds, product `P` leaves a covalently modified enzyme, substrate
`B` binds the modified form, and product `Q` leaves to regenerate free `E`:

```@example pingpong
using EnzymeRates
m = @enzyme_mechanism begin
    substrates: A, B
    products:   P, Q
    steps: begin
        E + A <--> E(A)
        E(A) <--> Estar(; residual = A - P) + P
        Estar(; residual = A - P) + B <--> Estar(B; residual = A - P)
        Estar(B; residual = A - P) <--> E + Q
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
`Estar_res_+A_-P` (the modified free enzyme) and `EstarB_res_+A_-P` (with `B`
additionally bound).
