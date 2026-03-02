# EnzymeMechanism Redesign: Per-Site Description with Conformational States

## Motivation

The current `EnzymeMechanism{Species, Reactions, EquilibriumSteps, ParamConstraints}`
requires every enzyme form and every reaction to be written explicitly. For oligomeric
enzymes this is infeasible:

- PFK tetramer (2 substrates, 2 products, 3 allosteric effectors per subunit):
  `(3×3×2×2×2×3)^4 × 2 ≈ 2.2 billion forms`
- Even a per-subunit expansion (54 forms) is intractable for King-Altman (O(n!))

The redesign separates **what each independent site does** from **how many copies
of each site exist** and **how many conformational states the enzyme has**. The existing
`EnzymeMechanism` is unchanged; `OligomericEnzymeMechanism` is a new wrapper around
it that composes multiple site mechanisms.

---

## Core Types

```julia
# Existing EnzymeMechanism type: UNCHANGED
# struct EnzymeMechanism{Species, Reactions, EquilibriumSteps, ParamConstraints}

# New wrapper for multi-site / multi-conformation mechanisms
# CatalyticMech : EnzymeMechanism type describing one catalytic subunit/domain
# CatalyticN    : Int — number of catalytic sites per enzyme molecule
# RegSites      : Tuple of ((ligand_symbols...,), multiplicity) pairs, one per
#                 regulatory site type; empty tuple for non-allosteric enzymes
# NConf         : Int — number of conformational states (1 = non-cooperative, 2 = two-state)
struct OligomericEnzymeMechanism{Metabolites, CatalyticMech, CatalyticN, RegSites, NConf}
end
```

`Metabolites` is a tuple of `(name, atoms)` pairs declared in the macro's `metabolites:`
block, shared across all site blocks. `CatalyticMech` is a full `EnzymeMechanism` encoding
enzyme states, reactions, RE/SS flags, and parameter constraints for the catalytic site.
`RegSites` stores regulatory sites as simple ligand lists — each entry is
`((:LigandA, :LigandB, ...), multiplicity)`. No King-Altman is needed for regulatory
sites: `Q_reg = 1 + Σ [Li]/K_Li` is assembled directly from the ligand names.

---

## Unified Rate Equation

The same assembly formula covers all values of `NConf`. There is no separate code
path — NConf=1 simply has no L parameter and no `_T` parameters.

```
Q_cat_R, N_cat_R  ← King-Altman on CatalyticMech with R-state parameters
                     N_cat_R = net-flux numerator (weighted spanning-tree sum for
                     all SS steps — the kcat contribution to the rate)
Q_cat_T, N_cat_T  ← King-Altman on CatalyticMech with T-state (_T suffix) parameters
                     T-state is always structurally present; set k_cat_T = 0 to
                     make it kinetically inactive (gives N_cat_T = 0 automatically)
                     NConf=1: no T-state terms; see simplified form below

For each regulatory site type i with multiplicity n_reg_i:
  Q_reg_i_R = 1 + Σ_ligands [Li]/K_Li      (direct sum from ligand list)
  Q_reg_i_T = 1 + Σ_ligands [Li]/K_Li_T    (NConf=2 only)

n_cat = CatalyticN

Q_R = Q_cat_R^n_cat × ∏_i Q_reg_i_R^n_reg_i    (per-enzyme partition function, R)
Q_T = Q_cat_T^n_cat × ∏_i Q_reg_i_T^n_reg_i    (per-enzyme partition function, T)
Z   = Q_R + L × Q_T                              (NConf=1: L absent → Z = Q_R)

Rate = E_total × n_cat × (N_cat_R × Q_cat_R^(n_cat-1) × ∏_i Q_reg_i_R^n_reg_i
                        + L × N_cat_T × Q_cat_T^(n_cat-1) × ∏_i Q_reg_i_T^n_reg_i) / Z
```

**NConf=1 has no L parameter and no `_T` parameters.** With no T-state terms and no
regulatory sites, the formula reduces algebraically:

```
Z = Q_cat_R^n_cat
Rate = E_total × n_cat × N_cat_R × Q_cat_R^(n_cat-1) / Q_cat_R^n_cat
     = E_total × n_cat × N_cat_R / Q_cat_R    (standard King-Altman for n_cat=1)
```

Implementation always produces the `N_cat × Q^(n_cat-1) / Q^n_cat` form; the Q
cancellation is a mathematical consequence, not a code branch.

### Generalized multiplicities (n_cat ≠ n_reg)

Catalytic and regulatory site multiplicities are independent. For ATCase
(6 catalytic c-chains, 6 regulatory r-chains): `CatalyticN=6`, `RegSites=(((:CTP,:ATP),6),)` —
same multiplicity, so `Q_R = (Q_cat × Q_reg)^6`. For a hypothetical enzyme with
2 catalytic and 1 regulatory subunit: `CatalyticN=2`, `RegSites=(((:I,),1),)` —
`Q_R = Q_cat^2 × Q_reg^1`, which cannot be written as `Q^n` for a single `n`.

### Why regulatory sites matter only for NConf=2

With NConf=1 (no L): `Z = Q_cat^n_cat × ∏ Q_reg_i^n_reg_i` and the `Q_reg` factors
cancel identically from numerator and denominator — regulatory sites have no effect.
With L≠0 (NConf=2): `Q_reg` appears asymmetrically in `Q_R` vs `Q_T` and does not
cancel. An activator binding preferentially to the R state increases `Q_R` more than
`Q_T`, shifting `Z` toward R and increasing the rate.

**Warning**: for NConf=1, all `RegSites` entries have no effect — Q_reg factors cancel
exactly from numerator and denominator. Dead-end inhibitors that bind the catalytic
site must be placed in `CatalyticMech` as dead-end complexes; they affect the rate
through Q_cat for any NConf.

### Site independence and binding within a site

**Same site, multiple ligands**: a regulatory site that binds AMP or ADP has
states `{E_r, E_AMP, E_ADP}` with **star topology** — all ligands bind the same
free enzyme form with no direct transitions between bound forms. King-Altman on
a star topology always gives:

```
Q_reg = 1 + [AMP]/K_AMP + [ADP]/K_ADP
```

In star topology, k_on and k_off cancel from every spanning tree ratio so parameters
are always Kd binding constants — never individual k_on/k_off pairs. The RE/SS
designation of regulatory binding steps is irrelevant (same formula either way).

Regulatory sites are declared with a `ligands:` list, not `steps:`. Q_reg is assembled
as a direct sum — no King-Altman machinery needed. Non-star topology (transitions
between bound forms in a regulatory site) is not supported and is a DSL error.

The outer ∏ in the assembly formula is over **different site types** (each a separate
entry in `RegSites`); the inner sum within one site is over mutually exclusive
occupancy states at that site.

**Different site types**: a catalytic site and a regulatory site can be simultaneously
occupied (they are on independent subunits/domains). This independence is what makes
the `∏_i Q_reg_i` product formula exact.

### Assumptions: when the factoring holds

The entire approach rests on two assumptions:

**Assumption 1 — Site independence within each conformation**: the rate constants at
every site (catalytic and regulatory, on every subunit) are unaffected by the occupancy
state of every other site, within the same conformational macrostate. This is one
statement that covers all forms of within-conformation kinetic coupling:

- Induced fit / heterocooperativity: binding at the regulatory site changes K_S at the
  catalytic site within the same conformation → **factoring fails**
- Sequential (KNF) conformational model: one subunit's conformational rate constants
  depend on a neighboring subunit's occupancy → **Q^n_cat factoring fails**
- Direct regulatory site–site coupling: binding at one regulatory site changes affinity
  at another regulatory site within the same conformation → **∏ Q_reg factoring fails**
- Half-of-the-sites reactivity: one subunit's catalytic constants depend on an adjacent
  subunit's occupancy → **factoring fails**

When Assumption 1 holds, the product distribution theorem (Cha 1968) gives exactly:

```
Q_enzyme_conformation = Q_cat^n_cat × ∏_i Q_reg_i^n_reg_i
```

The SS/RE classification of steps, mechanism topology (ordered, random, ping-pong),
and number of consecutive SS steps do not affect this factoring. For catalytic sites,
King-Altman handles any combination of SS and RE steps exactly. For regulatory sites,
the universal star-topology binding means Q_reg = 1 + Σ [Li]/K_Li — a direct sum,
no King-Altman required.

**Assumption 2 — R↔T is rapid equilibrium**: the conformational switch rate must be
fast relative to all catalytic and regulatory binding events so that L is an
equilibrium constant. If the switch is slow (rate comparable to kcat), L becomes a
steady-state ratio rather than an equilibrium constant and the two-state partition
function `Z = Q_R + L × Q_T` is no longer valid.

These are the only two conditions. Heterogeneous subunit types (different catalytic
kinetics on different subunit types) are handled by defining separate site types with
their own multiplicities — not a breakdown of the model.

### Parameters

Parameters use metabolite names in the **display layer only** (`rate_equation_string`,
`parameters`). Internal computation uses step-indexed names (K1, k1f, k1r, ...) as
today. The name conversion happens at the end, not in the core derivation.

RE steps (rapid equilibrium) → binding constants (Kd):

| Step | Parameter |
|------|-----------|
| `[E, F6P] ⇌ [EF6P]` | `K_F6P` |
| `[EF6P, ATP] ⇌ [ESATP]` | `K_ATP` |
| `[EPADP] ⇌ [EADP, FBP]` | `K_FBP` |

SS steps (steady state) → kinetic rate constants:

| Step | Parameters |
|------|-----------|
| `[ESATP] <--> [EPADP]` | `k_on_for_F6P_ATP_to_FBP_ADP`, `k_off_for_F6P_ATP_to_FBP_ADP` |

For regulatory sites, add `_reg{n}` suffix to distinguish from catalytic site params:
`K_AMP_reg1`, `K_ADP_reg1`.

For T-state (NConf=2), append `_T` to all parameters above:
`K_F6P_T`, `K_AMP_reg1_T`, etc.

**Name collision rule**: if the same metabolite appears in multiple RE steps of the
same site (e.g., two ATP binding events), append the enzyme state name of the complex
formed: `K_ATP_ESATP`, `K_ATP_ES2ATP`. If the same metabolite appears in different
site types, the `_reg{n}` suffix already disambiguates.

**Shared parameters**: `Keq` (same for all conformations), `E_total`, `L` (NConf=2
only).

### Reference enzyme form for L

`L = [E_T] / [E_R]` where `E` is the **bare enzyme with nothing bound** — the state
with all sites unoccupied in all site mechanisms simultaneously. This is the reference
form for King-Altman normalization in each conformation. For ping-pong, `E` is the
unmodified apo form (not the F modified form).

### When to use NConf=2

Two-state allosteric enzymes: PFK, hemoglobin, ATCase, pyruvate dehydrogenase,
aspartate transcarbamylase. Valid for any multiplicity including n_cat=1 (monomeric
enzyme with a conformational change that modulates activity). Competitive inhibitors
that bind the active site belong in the catalytic `EnzymeMechanism` as dead-end
complexes, not in `RegSites` (they affect the rate through `Q_cat`, which enters `Z`
regardless of conformation).

---

## Ping-Pong Mechanisms with NConf=2

Ping-pong mechanisms (e.g., pyruvate dehydrogenase) are compatible with NConf=2.
They have two ligand-free enzyme forms: **E** (unmodified) and **F** (modified).
Both exist in R and T conformations.

### L is well-defined

`L = [E_T]/[E_R]` — the conformational equilibrium of the reference free-enzyme form.
The conformational ratio for every other free-enzyme form is not independent: each
is determined by L and the parameters of the half-reactions connecting it to E.

For a single modified form F (standard ping-pong):

```
E_R ⇌ F_R    (K_ping_R: ping half-reaction equilibrium in R-state)
 ↕ L            ↕ L_F
E_T ⇌ F_T    (K_ping_T: ping half-reaction equilibrium in T-state)

L_F = L × K_ping_T / K_ping_R
```

(Wegscheider clockwise: `K_ping_R × L_F × (1/K_ping_T) × (1/L) = 1`.)

For mechanisms with n modified free-enzyme forms F₁, F₂, ... Fₙ (multi-step ping-pong),
each has a derived conformational ratio:

```
L_Fk = L × ∏ᵢ₌₁ᵏ (K_ping_i_T / K_ping_i_R)
```

All conformational ratios L_Fk are derived from L and mechanism parameters — L
remains the only conformational user parameter.

**Reference form**: for mechanisms with multiple free enzyme forms (all with zero
metabolites bound), the reference form E must be explicitly designated. By convention,
it is the first free-enzyme state listed in the catalytic site `states:` block. King-Altman
normalizes to this reference in both R and T states.

### Q_cat is a rational function for ping-pong

Normalizing relative to `E_R` (the reference form):

```
Q_cat_R = 1 + S1/K_S1 + S1/(K_S1·K_P1·P1) + S1·S2/(K_S1·K_P1·K_S2·P1)
               ↑E_R       ↑ES1_R               ↑F_R              ↑FS2_R
```

P1 appears in the denominator of the F_R and FS2_R terms. This is a rational function
of concentrations. King-Altman produces this correctly.

### All ping-pong conformational ratios are automatically handled

The constraints `L_F = L × K_ping_T / K_ping_R`, `L_G = L × K_ping1_T × K_ping2_T /
(K_ping1_R × K_ping2_R)`, etc. do not require any cross-conformation analysis —
they are automatically satisfied by `Z = Q_R + L × Q_T`.

King-Altman on the T-state mechanism normalizes to E_T = 1. For any free enzyme form
X reached via n ping half-reactions, the T-state K-A result is:

```
[X_T]/[E_T] = ∏ᵢ K_ping_i_T × f(concentrations)
```

Multiplying by L:

```
L × [X_T]/[E_T] = L × ∏ᵢ K_ping_i_T × f(concentrations)
```

From the thermodynamic cycle:

```
[X_T]/[E_R] = L_X × [X_R]/[E_R]
            = L × ∏ᵢ (K_ping_i_T/K_ping_i_R) × ∏ᵢ K_ping_i_R × f(concentrations)
            = L × ∏ᵢ K_ping_i_T × f(concentrations)   ✓
```

The T-state King-Altman values implicitly encode every L_X/L through the T-state
rate constants. Multiplying by L recovers the correct weight of every free-enzyme
form relative to E_R, for any number of ping-pong intermediates. No combined R+T
graph analysis is needed. Per-conformation Haldane with shared Keq is complete.

---

## Haldane and Wegscheider Analysis

### Per-conformation analysis is fully independent

For NConf=1, the analysis is identical to the existing `EnzymeMechanism` workflow:
apply Haldane/Wegscheider to the catalytic site to eliminate one dependent parameter
using Keq. Regulatory sites (star topology, no closed cycles) have no Haldane
conditions.

For NConf=2, apply the existing Haldane/Wegscheider machinery **twice** on the same
catalytic `EnzymeMechanism` — once with R-state parameters and once with T-state
parameters — both using the **same Keq**:

```
Haldane(CatalyticMech, R-state params, Keq)  →  eliminate one R-state rate constant
Haldane(CatalyticMech, T-state params, Keq)  →  eliminate one T-state rate constant
```

Sharing Keq across conformations is the only coupling between the two analyses. Keq
is thermodynamically fixed — the equilibrium constant for S→P is the same regardless
of which conformational state the enzyme is in when catalysis occurs.

### Cross-conformation Wegscheider is automatic for non-ping-pong mechanisms

Cross-conformation cycles appear to add constraints beyond per-conformation Haldane.
For example, the cycle connecting R and T through substrate binding:

```
E_R + S ⇌ ES_R  (K_S_R)
ES_R    ⇌ ES_T  (L_ES — conformational equilibrium for ES form)
ES_T    ⇌ E_T + S  (K_S_T)
E_T     ⇌ E_R   (1/L)
```

Wegscheider gives `L_ES = L × K_S_T / K_S_R`. But L_ES never appears as a user
parameter — it is encoded implicitly in `Z = Q_R + L × Q_T`. The formula already
places state-specific binding constants in Q_R and Q_T; L multiplies Q_T globally.
There is no independent equation for L_ES to violate.

The same reasoning applies to the SS catalytic step. The cross-conformation
Wegscheider condition for that cycle is:

```
K_eq_cat_R / K_eq_cat_T = (K_S_T × K_P_R) / (K_S_R × K_P_T)
```

But applying Haldane to both conformations with shared Keq already enforces:

```
Keq = K_S_R × K_eq_cat_R / K_P_R    (R-state Haldane)
Keq = K_S_T × K_eq_cat_T / K_P_T    (T-state Haldane)
```

Dividing these two equations yields exactly the cross-conformation condition above.
**Per-conformation Haldane with shared Keq is therefore sufficient for all
non-ping-pong mechanisms.** No additional cross-conformation analysis is needed.

### L is thermodynamically unconstrained

`L = [E_T]/[E_R]` represents the intrinsic free energy difference between
conformations in the absence of ligands. No Haldane condition constrains it — it is
an independent equilibrium constant, as free as any binding constant.

### Ping-pong is also automatic

As shown in the Ping-Pong section, every derived conformational ratio L_X is implicitly
encoded in [X_T]/[E_T] by T-state King-Altman. Per-conformation Haldane with shared
Keq is complete for any number of ping-pong intermediates; no additional analysis
is needed.

### Implementation

| Analysis | Applied to | When |
|---|---|---|
| Per-conformation Haldane (R-state) | `CatalyticMech` with R params | NConf ≥ 1 |
| Per-conformation Haldane (T-state) | Same `CatalyticMech` with T params | NConf = 2 |
| Regulatory sites | None (no cycles in star topology) | Always |
| Cross-conformation (non-ping-pong) | Automatic: shared Keq is sufficient | NConf = 2 |
| Cross-conformation (ping-pong) | Automatic: L_F implicit in T-state K-A | NConf = 2 |
| L | None (free parameter) | N/A |

---

## Validity Summary

All rows use the same assembly formula. `NConf` determines whether L and `_T`
parameters exist, not which formula to use. The formula is **exact** whenever the two
assumptions hold (site independence within each conformation; R↔T is RE).

### When the formula holds

| Mechanism | NConf | Notes |
|---|---|---|
| Any mechanism, n_cat=1 | 1 | Standard King-Altman; Q^n cancels to Q^1 |
| Any mechanism, n_cat>1 | 1 | n_cat× scaling, no cooperativity |
| Any mechanism, n_cat=1 | 2 | Monomeric allosteric; Z = Q_R + L×Q_T |
| Any mechanism, n_cat>1 | 2 | Cooperative; any SS/RE mix is exact |
| Ping-pong, any n_cat | 2 | Q_cat rational in concentrations; L = [E_T]/[E_R] |

The specific SS/RE assignment (which step is substrate binding, which is chemical
conversion, whether product release is SS) only determines what concentrations appear
in Q_cat and N_cat — it does not affect whether the formula is exact.

### When the formula breaks

| Violation | Effect |
|---|---|
| Rate constants at one site depend on occupancy at another site (induced fit, KNF, direct site–site cooperativity, half-of-the-sites) | Factoring invalid; need joint mechanism with all coupled states explicit |
| R↔T switch rate comparable to kcat (not RE) | L is not an equilibrium constant; two-state partition function is invalid |
| Ping-pong Q_cat contains product concentrations in denominator | Formula and numerics are correct; `rate_equation_string` display and symbolic factoring need extension for rational functions |

---

## DSL

Single macro `@enzyme_mechanism`. Old syntax (no `site:` or `conformations:`) is
unchanged and produces `EnzymeMechanism` directly.

```julia
# --- Old syntax: unchanged ---
m = @enzyme_mechanism begin
    species: begin
        substrates: S[C], ATP[CN5P3]
        products:   P[C], ADP[CN5P2]
        enzymes:    E, ES, EP
    end
    steps: begin
        [E, S]  ⇌ [ES]
        [ES]   <--> [EP]
        [EP]    ⇌ [E, P]
    end
end

# --- New syntax: homotetramer, two-substrate catalysis + allosteric regulatory site ---
m = @enzyme_mechanism begin
    metabolites: F6P[C6P], ATP[C10N5P3], FBP[C6P2], ADP[C10N5P2], AMP[C10N5P]

    conformations: 2    # NConf=2: two-state (R active, T)

    site :catalytic multiplicity=4 begin
        states: E_c, EF6P_c, EF6PATP_c, EFBPADP_c, EADP_c
        steps: begin
            [E_c, F6P]      ⇌ [EF6P_c]
            [EF6P_c, ATP]   ⇌ [EF6PATP_c]
            [EF6PATP_c]    <--> [EFBPADP_c]
            [EFBPADP_c]     ⇌ [EADP_c, FBP]
            [EADP_c]        ⇌ [E_c, ADP]
        end
    end

    site :regulatory multiplicity=4 begin
        ligands: AMP, ADP
    end
end
```

Names listed in the `metabolites:` block are recognized as global metabolites across
all site blocks. All other names in `states:` and `steps:` expressions are local enzyme
states. `site :catalytic` produces a full `EnzymeMechanism` (King-Altman on the
`steps:` block). `site :regulatory` uses a `ligands:` list — the macro generates a
`RegSites` entry directly; no King-Altman is run. Exactly one `:catalytic` site is
required; `:regulatory` sites are optional and may repeat for different site types.

`conformations: 2` triggers two-conformation assembly. Omitting `conformations:` (or
writing `conformations: 1`) gives NConf=1. `conformations: N` for N > 2 is not
currently supported.

---

## rate_equation Assembly

One code path for all `OligomericEnzymeMechanism` types:

```
rate_equation(::OligomericEnzymeMechanism{Mets, CatMech, CatN, RegSites, NConf}, conc, params)
  └─ King-Altman on CatMech with R-state params → Q_cat_R, N_cat_R
     NConf=2: King-Altman on CatMech with _T params → Q_cat_T, N_cat_T; read L
     NConf=1: Z = Q_cat_R^CatN; return E_total × CatN × N_cat_R / Q_cat_R
     For each ((Li₁, Li₂, ...), n_reg) in RegSites:
       Q_reg_R = 1 + [Li₁]/K_Li₁ + [Li₂]/K_Li₂ + ...
       Q_reg_T = 1 + [Li₁]/K_Li₁_T + [Li₂]/K_Li₂_T + ...    (NConf=2 only)
     Q_R = Q_cat_R^CatN × ∏ Q_reg_R^n_reg
     Q_T = Q_cat_T^CatN × ∏ Q_reg_T^n_reg
     Z   = Q_R + L × Q_T
     return E_total × CatN × (N_cat_R × Q_cat_R^(CatN-1) × ∏ Q_reg_R^n_reg
                            + L × N_cat_T × Q_cat_T^(CatN-1) × ∏ Q_reg_T^n_reg) / Z
```

`NConf` controls whether L, `_T` parameters, and RegSites T-state entries exist.
The assembly formula is the same for all `NConf` values; the NConf=1 early-return
above is an optimization (avoids the n_cat-1 exponent of a rational Q_cat).

---

## Open Design Questions

1. **Ping-pong in `rate_equation_string`**: Q_cat terms contain product concentrations
   in denominators (rational functions, not polynomials). The polynomial string builder
   and factoring pipeline assume non-negative exponents and will need extension.
