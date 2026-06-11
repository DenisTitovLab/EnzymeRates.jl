# Rate equations from textbooks

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
(; K_P_E, K_S_E, k_ES_to_EP, Keq, E_total) = params
(; S, P) = concs
# Haldane constraints:
k_EP_to_ES = (1 / Keq) * K_P_E * (1 / K_S_E) * k_ES_to_EP
v = E_total * (k_ES_to_EP * S / K_S_E - k_EP_to_ES * P / K_P_E) / (1 + P / K_P_E + S / K_S_E)
```

The default mode is `Reduced`.
In `Reduced` mode the string has four sections:

1. A `params` destructuring line listing the independent fitted parameters
   plus `Keq` and `E_total`.
2. A `concs` destructuring line listing the concentration symbols.
3. An optional `# Haldane constraints:` (or `# Wegscheider constraints:`)
   section where each dependent rate constant is expressed in terms of `Keq`
   and the independent parameters.
4. The final `v = E_total * (num) / (den)` line.

## Meaning of parameter names

| Package symbol | Role | Units |
|:---|:---|:---|
| `K_S_E` | Dissociation constant of substrate `S` from the `E`–`S` complex. For a one-substrate, one-product Michaelis–Menten enzyme it is analogous to the Michaelis constant `Km` under the rapid-equilibrium approximation, but it is not the same as `Km` for more complex mechanisms. | M |
| `K_P_E` | Dissociation constant of product `P` from the `E`–`P` complex. | M |
| `k_ES_to_EP` | Rate of conversion of enzyme species `ES` to enzyme species `EP`. | 1/s |
| `k_EP_to_ES` | Rate of the reverse conversion, `EP` to `ES`. | 1/s |
| `kon_S_E` | Rate constant for binding of `S` to `E`. | 1/(s·M) |
| `koff_S_E` | Rate of dissociation of `S` from the `E`–`S` complex. | 1/s |
| `Keq` | Equilibrium constant of the overall reaction; always user-supplied. | unitless |
| `E_total` | Total enzyme concentration. | M |

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

## Full mode

Passing `Full` lists all raw rate constants, before any thermodynamic
reduction:

```@example textbook_ex
parameters(m, Full)
```

`Full` mode includes both `k_ES_to_EP` and `k_EP_to_ES` as independent symbols,
so there is no constraint section.
See [Rapid equilibrium vs steady state](@ref) for the contrast between RE and
SS parameters, and how adding SS steps changes this list.
