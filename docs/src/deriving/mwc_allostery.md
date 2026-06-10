# MWC allostery

Many enzymes exist in two conformations whose relative population regulators can
shift. The Monod–Wyman–Changeux (MWC) model [Monod1965](@cite) describes this
with two states: an **active** state (A) and an **inactive** state (I). In the
classic MWC notation R ≡ A and T ≡ I; the package uses A/I throughout, so
`K_A_…` names correspond to the literature's R-state constants.

Both conformations run the same catalytic cycle. A coupling constant `L`
weights the inactive-state population relative to the active state for the bare
enzyme. Regulatory ligands bind in rapid equilibrium and shift the A/I balance
by binding preferentially to one conformation, embodying the classical
partition-function logic.

## Setting up an MWC mechanism

Use [`@allosteric_mechanism`](@ref) and tag each catalytic step with its
allosteric-state behavior:

```@example mwc
using EnzymeRates
allo = @allosteric_mechanism begin
    substrates: S
    products:   P
    catalytic_multiplicity: 2
    allosteric_regulators: A::OnlyA, I::OnlyI
    catalytic_steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)       :: NonequalAI
        E(P) ⇌ E + P        :: EqualAI
    end
    regulatory_site(multiplicity = 2): begin
        ligands: A
    end
    regulatory_site(multiplicity = 2): begin
        ligands: I
    end
end
parameters(allo)
```

```@example mwc
print(rate_equation_string(allo))
```

## The four allosteric-state tags

Each catalytic kinetic group and each regulatory ligand carries one of four
tags that controls how the package renders A- and I-state symbols.

| Tag | Meaning | Example symbol |
|-----|---------|----------------|
| `:EqualAI` | One shared symbol for both states | `K_S_E` |
| `:NonequalAI` | Independent A and I symbols | `k_A_ES_to_EP`, `k_I_ES_to_EP` |
| `:OnlyA` | Present in the active state only; zeroed in the inactive polynomial | `K_A_Areg` |
| `:OnlyI` | Present in the inactive state only; zeroed in the active polynomial | `K_I_Ireg` |

In the example above:

- The two binding steps are tagged `:EqualAI`, so they share `K_S_E` and
  `K_P_E` — no A/I token appears.
- The isomerization step is tagged `:NonequalAI`, producing two independent
  rate constants `k_A_ES_to_EP` and `k_I_ES_to_EP`.
- The activator `A` is `:OnlyA`: it contributes `K_A_Areg` and appears only
  in the active-state polynomial.
- The inhibitor `I` is `:OnlyI`: it contributes `K_I_Ireg` and appears only
  in the inactive-state polynomial.
- `L` is the MWC coupling constant for the bare enzyme equilibrium.

## Catalytic groups cannot be `:OnlyI`

The active state must always be capable of catalysis. Tagging a catalytic
step `:OnlyI` is therefore a hard constructor error. Only regulatory ligands
may carry `:OnlyI`.

```@example mwc
try
    @allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_steps: begin
            E + S ⇌ E(S)   :: OnlyI
            E(S) <--> E(P)  :: EqualAI
            E(P) ⇌ E + P   :: EqualAI
        end
    end
catch err
    showerror(stdout, err)
end
```

## Partition-function structure

For a mechanism with `catalytic_multiplicity` equal to `cat_n`, the package
assembles the rate equation as

```
v = E_total * num / den
```

where

```
num = cat_n * (N_A * Q_A^(cat_n-1) * W_A + L * N_I * Q_I^(cat_n-1) * W_I)
den =          Q_A^cat_n * W_A           + L * Q_I^cat_n           * W_I
```

- `N_A`, `Q_A` are the active-state catalytic numerator and partition function
  (the King–Altman/Cha polynomials for the catalytic cycle).
- `N_I`, `Q_I` are the corresponding inactive-state quantities, obtained by
  zeroing `:OnlyA` symbols and renaming `:NonequalAI` symbols to their
  I-state counterparts.
- `W_A` and `W_I` are products of regulatory-site factors, each raised to
  that site's multiplicity. For a site with two ligands at multiplicity `m`,
  the factor is `(1 + lig1/K_lig1 + lig2/K_lig2)^m`.
- `L` weights the inactive-state branch throughout, encoding both the
  inactive-enzyme population and the allosteric coupling.

In the rendered equation above, `cat_n = 2`. The active branch carries
`(1 + A / K_A_Areg) ^ 2` (activator favors A-state), the inactive branch
carries `(1 + I / K_I_Ireg) ^ 2` (inhibitor favors I-state), and each
partition function appears squared in the denominator.

When any catalytic group is `:OnlyA`, the inactive catalytic cycle cannot
close. The package forces `N_I = 0` in that case to maintain Haldane
thermodynamic consistency — the inactive branch still contributes to the
denominator as enzyme mass, but carries no forward flux.

## Known rendering quirk in Haldane constraints

In the Haldane-constraint lines at the top of the rendered equation, the
right-hand side may reference a bare active-state base name rather than the
A-suffixed symbol listed in `params`. In the example above:

```
k_EP_to_ES = (1 / Keq) * K_P_E * (1 / K_S_E) * k_ES_to_EP
```

Here `k_ES_to_EP` is the un-prefixed base name, even though the independent
parameter list uses `k_A_ES_to_EP`. The numerical result is identical because
the assignment block sets `k_ES_to_EP` from `k_A_ES_to_EP` before evaluating
this line. This is confirmed behavior; the rendered string is correct.

## See also

- [`@allosteric_mechanism`](@ref) — DSL reference
- [`parameters`](@ref) — list independent parameter names
- [`rate_equation_string`](@ref) — inspect the derived equation
