# The Cha / King–Altman algorithm

EnzymeRates derives rate equations at compile time using the King–Altman
schematic method [KingAltman1956](@cite) extended by Cha's rapid-equilibrium
treatment [Cha1968](@cite).
Each unique mechanism type triggers a full symbolic derivation inside
`@generated` functions, so the derived equation is embedded directly
into the compiled code and the runtime call is allocation-free.

## Overview of the algorithm

The mechanism graph is partitioned into **rapid-equilibrium (RE) segments**
— sets of enzyme forms linked by RE steps — connected to each other by SS
steps.

**Within each RE segment**, the package picks the form with the fewest bound
metabolites as the **segment root** (the free-enzyme reference form), with
ties broken toward no covalent residual and then by name.
Every other form in the segment gets an alpha factor — its relative abundance
expressed as a ratio to the root.
Referencing to the free-enzyme form is what produces the readable
`1 + S / K_S_E + P / K_P_E` denominator in the final equation.

<!-- verified: _segment_root at src/rate_eq_derivation.jl:263-267 -->
<!-- verified: _compute_alpha at src/rate_eq_derivation.jl:277-322 -->

**Across segments**, the SS steps form an inter-segment rate matrix.
The Cha numerator and denominator come from cofactor determinants of that
matrix, with each diagonal entry weighted by the sigma sum of alpha factors
for its segment.

<!-- verified: _raw_symbolic_rate_polys at src/rate_eq_derivation.jl:339-399 -->
<!-- verified: sym_det at src/sym_poly_for_rate_eq_derivation.jl:88-114 -->

## Symbolic algebra: Laurent polynomials

The entire derivation works in **Laurent polynomial** arithmetic:
a monomial is a sorted list of `(Symbol, Int)` exponent pairs (exponents may
be negative), and a polynomial is a `Dict` mapping monomials to
`Rational{Int}` coefficients.

<!-- verified: MONO, POLY at src/sym_poly_for_rate_eq_derivation.jl:9-10 -->

Alpha factors and the numerator/denominator are maintained as fractional
Laurent POLYs throughout — the package never forms a common denominator or
expands the full rational expression.
This keeps intermediate expression size small for mechanisms with many
enzyme forms.

## Keeping the equation finite at zero concentration

A naive segment-referenced derivation can produce intermediate terms of the
form `conc / (conc * K)`, which simplify to `1/K` — a concentration-free
constant in the polynomial but numerically safe.
Worse, some coupling patterns leave a genuine `1/conc` factor in the raw
numerator or denominator, which would cause division by zero when a substrate
concentration is zero.

The package resolves this with `_reduce_conc_lowest_terms`: after building the
raw polynomials, it scans every monomial in num ∪ den, finds the minimum
exponent for each concentration symbol across all monomials (counting zero for
monomials where the symbol is absent), and shifts every exponent by that
minimum.
If the minimum is negative, the shift raises it to zero, clearing the `1/conc`
coupling.
If the minimum is positive (every monomial contains that concentration),
the shift factors it out, giving a cleaner expression.
The shift acts **only on concentration symbols**, never on parameter symbols,
so it cannot drop a fitted parameter.

<!-- verified: _reduce_conc_lowest_terms at src/sym_poly_for_rate_eq_derivation.jl:56-81 -->
<!-- verified: _concentration_symbols at src/rate_eq_derivation.jl:249-255 -->

## Seeing it in action

The textbook reversible Michaelis–Menten mechanism illustrates the result:

```@example cha
using EnzymeRates
m = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
    end
end
print(rate_equation_string(m))
```

Look at the denominator: `1 + P / K_P_E + S / K_S_E`.
The leading `1` comes from the free-enzyme term `E` (alpha = 1 at the segment
root) plus the concentration-GCD step that clears any residual `1/conc`
coupling.
The result is division-free and finite at every concentration, including zero.

The following call evaluates the rate at zero substrate and zero product:

```@example cha
rate_equation(m, (; S = 0.0, P = 0.0),
              (; K_S_E = 1.0, K_P_E = 1.0, k_ES_to_EP = 10.0, Keq = 2.0, E_total = 1.0))
```

The result is `0.0` — zero net rate at zero concentrations, and no
division-by-zero error.
That is the entire point of the concentration-GCD step.

## Compile-time derivation and the term limit

Derivation runs at compile time inside [`rate_equation`](@ref) and
[`rate_equation_string`](@ref), so each unique mechanism type
pays the symbolic cost once at first call and is free at runtime.
Very large mechanisms can be slow to compile, exhaust memory, or overflow the
stack; this is inherent to the type-parameter-based architecture.

To prevent unbounded compilation, the package aborts if the intermediate
determinant expansion exceeds `MAX_RATE_EQUATION_TERMS = 5000` polynomial
terms.

<!-- verified: MAX_RATE_EQUATION_TERMS at src/sym_poly_for_rate_eq_derivation.jl:7 -->
<!-- verified: sym_det term-count check at src/sym_poly_for_rate_eq_derivation.jl:102-111 -->

## Zero-allocation runtime

The emitted expression is a balanced binary `+`/`*` tree: every arithmetic
node has exactly two operands.
Julia inlines binary `+(::Float64, ::Float64)` into fused scalar arithmetic,
but falls back to a varargs path that boxes the operand tuple once an operator
chain exceeds roughly 30 terms.
The balanced tree avoids that path entirely, keeping every [`rate_equation`](@ref)
call allocation-free and sub-100 ns regardless of mechanism size.

<!-- verified: _nest_binary at src/sym_poly_for_rate_eq_derivation.jl:162-170 -->
