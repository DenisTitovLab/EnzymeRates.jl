# Derivation tutorial

This tutorial shows how to define the reversible Michaelis–Menten mechanism
and derive its rate equation as a `String`. After reading it, you will
recognize the package's parameter names as the textbook constants you already
know.

## Defining the mechanism

Use [`@enzyme_mechanism`](@ref) to declare the mechanism.
The `⇌` arrow marks a step as rapid equilibrium (RE);
`<-->` marks it as steady state (SS).

```jldoctest textbook
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products:   P
           steps: begin
               E + S ⇌ E(S)
               E(S) <--> E(P)
               E(P) ⇌ E + P
           end
       end;

julia> m isa EnzymeMechanism
true
```

The two binding steps are rapid equilibrium; the central isomerization
`E(S) <--> E(P)` is steady state.
That single SS step is the rate-limiting catalytic step.

## Deriving the rate equation

[`rate_equation_string`](@ref) **returns** the symbolic rate equation as a
multi-line `String` (it does not print).
Call `print` to display it without escaped newlines:

```jldoctest textbook
julia> print(rate_equation_string(m))
(; K_EP_to_E_P, K_ES_to_E_S, k_ES_to_EP, Keq, E_total) = params
(; S, P) = concs
# Haldane constraints:
k_EP_to_ES = K_EP_to_E_P * k_ES_to_EP / (K_ES_to_E_S * Keq)
v = E_total * (k_ES_to_EP / K_ES_to_E_S) * (S - P / Keq) / (1 + S / K_ES_to_E_S + P / K_EP_to_E_P)
```

The string reads, in order:

1. a `params` destructuring line: the independent fitted constants, then `Keq` and
   `E_total`;
2. a `concs` destructuring line: the concentrations;
3. the reference lines under `# Haldane constraints:` (after `# Wegscheider constraints:`
   when a cycle of bindings ties constants together): each dependent constant in `Keq` and
   the independent constants, for reference, since the law reads only the independent
   ones;
4. when the normalizer has a factor that is not a single term, or an integer to a
   fractional power, the line `N = …` ([Normalization](@ref));
5. the law, `v = E_total * (Q) * (H) / (D)`. `H = S - P / Keq` vanishes at equilibrium,
   and `Q` and `D` have positive terms.

## Normalization

Every term of `D` is dimensionless unless the law has no normalizer (rule 4): the package
divides numerator and denominator by a normalizer with `D`'s unit, chosen by the first rule
that applies.

1. **Constant term.** When `D` has a term free of concentrations, divide by that term's
   coefficient, the part of `D` free of concentrations; every term is then dimensionless.
   In a rapid-equilibrium law whose only form that binds nothing is free E, as in
   `1 + S / K_ES_to_E_S + P / K_EP_to_E_P` above, that part is free E's weight, and `D`
   sums the enzyme forms' abundances relative to free E.
2. **Cleland's normalizer.** When `D` holds the term `ΠS/X` for every substrate `X` and
   the term `ΠS`, `ΠS` the product of the substrates, divide by
   `N = Π_X c(ΠS/X) / c(ΠS)^(n−1)`, `c(m)` the coefficient of the term `m` and `n` the
   number of substrates. Each substrate term then reads as its concentrations over their
   Michaelis constants, `A·B·C/(K_a·K_b·K_c)`. The ping-pong bi-bi of
   [Ping-pong mechanisms](@ref) takes this rule.
3. **Virtual constant term.** Otherwise take `D`'s terms in a fixed order (lowest total
   degree first, then fewest metabolites other than substrates, then metabolite order) and
   weigh the first affinely independent ones so that the product of their coefficients,
   each to its weight, has `D`'s unit. Fractional weights print as roots.
4. **No normalizer.** When no such weights exist, no product of the coefficients has
   `D`'s unit: the law reads some concentrations only through their ratio, as a
   rapid-equilibrium Theorell–Chance step `E(A) + B ⇌ E(Q) + P` makes every term carry `B`
   or `P`. The law then prints as derived, and its terms keep their units.

A factor of the normalizer that is not a single term, or an integer to a fractional power,
prints once, as the line `N = …`, and the terms it divides read `/ N`.

## The coefficient form

`coefficients = true` prints the same law with one name per coefficient. It too runs from
top to bottom: the `params` and `concs` lines, a comment line, `N` when it prints, one line
per coefficient in the fitted constants, for an allosteric law its polynomials in those
coefficients and its site factors, then `v`. The reference lines print only in the full
law.

```jldoctest textbook
julia> print(rate_equation_string(m; coefficients = true))
(; K_EP_to_E_P, K_ES_to_E_S, k_ES_to_EP, Keq, E_total) = params
(; S, P) = concs
# The coefficients are functions of the fitted constants, not parameters:
n_0 = k_ES_to_EP / K_ES_to_E_S
d_S = 1 / K_ES_to_E_S
d_P = 1 / K_EP_to_E_P
v = E_total * n_0 * (S - P / Keq) / (1 + d_S * S + d_P * P)
```

A coefficient is named by the concentrations it multiplies: `d_A_B` for `A * B` and
`d_A²_B` for `A^2 * B`; `n_0` is the numerator's constant term. The fit still runs over the
fitted constants.

## Meaning of parameter names

| Package symbol | Role | Units |
|:---|:---|:---|
| `K_ES_to_E_S` | Dissociation constant of substrate `S` from the `E`–`S` complex. For a one-substrate, one-product Michaelis–Menten enzyme it is analogous to the Michaelis constant `Km` under the rapid-equilibrium approximation, but it is not the same as `Km` for more complex mechanisms. | M |
| `K_EP_to_E_P` | Dissociation constant of product `P` from the `E`–`P` complex. | M |
| `k_ES_to_EP` | Rate of conversion of enzyme species `ES` to enzyme species `EP`. | 1/s |
| `k_EP_to_ES` | Rate of the reverse conversion, `EP` to `ES`. | 1/s |
| `k_E_S_to_ES` | Rate constant for binding of `S` to `E`. | 1/(s·M) |
| `k_ES_to_E_S` | Rate of dissociation of `S` from the `E`–`S` complex. | 1/s |
| `Keq` | Equilibrium constant of the overall reaction; always user-supplied. | unitless |
| `E_total` | Active-site (protomer) concentration. For a monomer this equals the total enzyme concentration; for an allosteric enzyme with `N` catalytic subunits it is `N` × the oligomer concentration, so allosteric `kcat` is reported per active site. | M |

## Inspecting parameters and metabolites

The `@example` block below shows the fitted-parameter tuple and the
concentration symbols returned by [`parameters`](@ref) and
[`metabolites`](@ref):

```@example textbook_ex
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
(parameters(m), metabolites(m))
```

`parameters(m)` returns the symbols the fitter expects: the three independent
constants plus `Keq` and `E_total`.
`metabolites(m)` returns the concentration symbols the rate equation reads
at runtime.

## Every constant of the rate law

[`all_parameters`](@ref) adds the dependent constants to a set of independent values,
each computed from its Haldane or Wegscheider relation:

```@example textbook_ex
all_parameters(m, (K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0, Keq = 4.0))
```

Here the Haldane relation fixes the reverse catalytic rate `k_EP_to_ES` from `Keq`, the
two binding constants and the forward rate. See
[Rapid equilibrium vs steady state](@ref) for the contrast between RE and SS
parameters, and how adding SS steps changes the parameter list.
