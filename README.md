# EnzymeRates.jl

Identify the best enzyme rate equation from kinetic data.

Given a reaction definition and experimental rate measurements at varying
substrate, product, and regulator concentrations, EnzymeRates enumerates all
biochemically valid mechanisms, fits each to the data, and selects the
mechanism with the fewest parameters that adequately describes the data based
on cross-validation.

```julia
# Define the reaction
rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products:   P[C], Q[N]
end

# Construct the selection problem (enumerates all mechanisms)
prob = IdentifyRateEquationProblem(rxn, data; Keq=5.0)

# Identify the best rate equation
results = identify_rate_equation(prob)

# Inspect the winner
rate_equation_string(results.best.mechanism)
rate_equation(results.best.mechanism, results.best.params, concentrations)
```

In addition to automated model selection, the package can be used to:

- Define enzyme mechanisms using a concise DSL and derive their QSSA rate
  equations (compiled into zero-allocation numeric functions)
- Fit a known mechanism's rate equation to kinetic data
- Analyze structural identifiability of a mechanism
- Enumerate all biochemically valid mechanisms for a given reaction

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/DenisTitovLab/EnzymeRates.jl")
```

## Quick Start

```julia
using EnzymeRates

# Define a reversible Uni-Uni mechanism
m = @enzyme_mechanism begin
    species: begin
        substrates: S[C]
        products:   P[C]
        enzymes:    E, ES[C]
    end
    steps: begin
        [E, S] <--> [ES]
        [ES] <--> [E, P]
    end
end

# Independent rate constants + Keq + E_total
params = (k1f=3.2, k2f=2.5, k2r=1.1, Keq=5.0, E_total=1.0)
concs = (S=0.7, P=0.3)

# Compute rate: zero allocations, compiled at first call
v = rate_equation(m, params, concs)

# Human-readable equation
rate_equation_string(m)
```

## Defining Mechanisms

### `@enzyme_mechanism` macro

`@enzyme_mechanism` requires explicit `species:` and `steps:` blocks.
Species include substrates, products, regulators, and enzyme forms. Atoms
use bracket syntax (`S[C6H12O6]`); multiple atoms are concatenated
(`A[C2H3]` means 2 C and 3 H).

Steps use `<-->` for steady-state or `⇌` for rapid-equilibrium:

```julia
m = @enzyme_mechanism begin
    species: begin
        substrates: S[C6H12O6]
        products:   P[C6H12O6]
        enzymes:    E, ES[C6H12O6]
    end
    steps: begin
        [E, S] ⇌ [ES]          # rapid-equilibrium binding
        [ES] <--> [E, P]        # steady-state catalysis
    end
end
```

An optional `constraints:` block can constrain rate parameters:

```julia
m = @enzyme_mechanism begin
    species: begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
        enzymes:    E, EA[C], EB[N], EAB[CN], EPQ[CN], EQ[N]
    end
    steps: begin
        [E, A] <--> [EA]
        [E, B] <--> [EB]
        [EA, B] <--> [EAB]
        [EB, A] <--> [EAB]
        [EAB] <--> [EPQ]
        [EPQ] <--> [EQ, P]
        [EQ] <--> [E, Q]
    end
    constraints: begin
        k4f = k3f           # same forward rate for A binding
        k4r = k3r           # same reverse rate for A binding
    end
end
```

### `@enzyme_reaction` macro

Define an overall reaction (substrates, products, optional regulators) without
specifying the mechanism steps. Used as input for mechanism enumeration:

```julia
rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products:   P[C], Q[N]
    regulators: I[P2]
end
```

### Programmatic construction

```julia
species = (
    ((:S, ((:C, 1),)),),       # substrates
    ((:P, ((:C, 1),)),),       # products
    (),                        # regulators
    ((:E, ()), (:ES, ((:C, 1),))),  # enzyme forms
)
reactions = (
    ((:E, :S), (:ES,)),
    ((:ES,), (:E, :P)),
)

m = EnzymeMechanism(species, reactions)
```

## Rate Equations

`rate_equation(m, params, concs)` computes the steady-state rate (net
consumption of the first substrate, normalized by its stoichiometric
coefficient):

- `params`: `NamedTuple` of independent rate constants (`k1f`, `k1r`, ...),
  `Keq` (equilibrium constant), and `E_total` (total enzyme concentration).
  Dependent rate constants (determined by Haldane/Wegscheider constraints)
  are computed internally and should not be included.
- `concs`: `NamedTuple` of metabolite concentrations (`S`, `P`, ...)

Use `parameters(m)` to get the list of required parameter names and
`fitted_params(m)` to get just the independent rate constants (excluding
`Keq` and `E_total`).

## Mechanism Enumeration

Given an `EnzymeReaction`, `enumerate_mechanisms` generates all biochemically
valid mechanisms as a lazy iterator of `MechanismSpec` instances. Each spec
can be converted to an `EnzymeMechanism` via `EnzymeMechanism(spec)`.

```julia
rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products:   P[C], Q[N]
end

# Lazy iterator over all valid mechanisms
iter = enumerate_mechanisms(rxn)
length(iter)  # total count (O(1), precomputed)

for spec in iter
    m = EnzymeMechanism(spec)
    # use m...
end
```

### Pipeline stages

The enumeration runs in five stages, exposed via `enumerate_mechanism_stages`:

| Stage | Description |
|-------|-------------|
| **1. Enzyme forms** | Enumerate all distinguishable enzyme states: free enzyme, substrate-bound, product-bound, and ping-pong residual intermediates. |
| **2. Catalytic topologies** | Build minimal catalytic cycles (sequential and ping-pong), then combine into multi-cycle unions. A purity filter removes hybrid topologies (see below). |
| **3. Activator configurations** | For each regulator, decide if it is an activator (essential or non-essential) or an inhibitor. Activators add new catalytic forms; inhibitors are handled in the next stage. |
| **4. Dead-end complexes** | Each inhibitor can independently form a dead-end complex at each topology form, giving `(2^r_inh)^n_topo` configurations per activator setup. |
| **5. RE/SS + constraints** | For each dead-end topology, enumerate all rapid-equilibrium (RE) vs. steady-state (SS) assignments for each elementary step, plus optional parameter constraints for equivalent binding steps. This stage is lazy — variants are generated on iteration. |

### Enzyme forms

Each binding site on the enzyme is distinguishable. A site can be:

- **Empty**: the metabolite is not bound (`atoms = nothing`)
- **Fully occupied**: the metabolite is fully bound (`atoms = full_atoms`)
- **Residual** (ping-pong only): partial atom content remains after a product
  has been released from the substrate site

For example, in a Bi-Bi ping-pong reaction `A[CX], B[N] → P[C], Q[NX]`:

| Form | A-site | B-site | P-site | Q-site | Description |
|------|--------|--------|--------|--------|-------------|
| `E` | — | — | — | — | Free enzyme |
| `E_A` | `[C,X]` | — | — | — | A bound |
| `E_A_B` | `[C,X]` | `[N]` | — | — | Both substrates bound (ternary complex) |
| `E_X` | `[X]` | — | — | — | Free intermediate (residual from A after P released) |
| `E_X_B` | `[X]` | `[N]` | — | — | Intermediate with B bound |
| `E_X_P` | `[X]` | — | `[C]` | — | Intermediate before P release |

### Allowed elementary steps

An elementary step connects two enzyme forms that differ by exactly one event.
The `edge_class` function classifies each valid transition:

| Step type | Rule | Example |
|-----------|------|---------|
| **Binding** | Exactly one site changes from empty to fully occupied. All other sites unchanged. | `E → E_A` (A binds) |
| **Release** | Exactly one site changes from occupied (full or residual) to empty. All other sites unchanged. For residual sites, the released metabolite is identified by matching the residual atoms to a product's atom signature. | `E_A → E + A` (full release), `E_X → E` (residual release) |
| **Isomerization** | Two or more core sites change simultaneously, with total atom balance preserved. Two sub-types: | |
| — *Standard* | All substrate sites switch to products (or vice versa). Every core site must differ between the two forms. No substrate site may remain occupied. | `E_A_B → E_P_Q` |
| — *Ping-pong* | A substrate site undergoes partial transformation (producing a residual). Non-differing substrate sites may remain occupied. Only total atom balance across all core sites is required. | `E_A → E_X_P` (A partially transforms, P appears) |

**Invalid transitions** (no elementary step exists):

- Two or more sites change but atom balance is violated
- A substrate site goes from one partial occupancy to another without a matching
  product change
- A regulatory or extra site changes simultaneously with a core site
- A site changes from full to residual without a corresponding product appearance
  (residual can only appear via isomerization)

### Pure topology filter

After combining individual cycles into multi-cycle topologies, a purity filter
removes biochemically implausible hybrids. A topology must be either:

- **Pure sequential**: no residual forms at all. The enzyme follows a standard
  ternary-complex pathway where all substrates bind before products are released.
- **Pure ping-pong**: has a free enzyme intermediate (a form carrying only
  residual atoms with all other core sites empty) AND does **not** contain the
  all-substrates-fully-bound form (ternary complex).

Topologies that mix both patterns — e.g., having both a ternary complex and
a residual intermediate, or having residuals without a free intermediate — are
rejected. These hybrids are biochemically implausible because ping-pong enzymes
by definition release the first product before the second substrate binds.

### Ping-pong catalytic cycle example

The classic Bi-Bi ping-pong mechanism (`A[CX], B[N] → P[C], Q[NX]`):

```
E + A ⇌ E_A ↔ E_X_P ⇌ E_X + P
                         ↓
                    E_X + B ⇌ E_X_B ↔ E_Q ⇌ E + Q
```

Key features:
- `E_A ↔ E_X_P` is a ping-pong isomerization: A partially transforms, P appears
  on the enzyme, residual X remains on the A-site
- `E_X_P → E_X + P` releases the first product
- `E_X` is the free intermediate (only residual X, everything else empty)
- `E_X_B ↔ E_Q` is another isomerization: X and B combine to form Q
- The ternary complex `E_A_B` never appears in this cycle

## Querying a Mechanism

Mechanisms are validated at construction (elementary-step structure, atomic
conservation, regulator balance). There is no separate `validate` API.

```julia
substrates(m)              # substrates with stoichiometric multiplicity
products(m)                # products with stoichiometric multiplicity
regulators(m)              # regulators
enzyme_forms(m)            # distinct enzyme states
metabolites(m)             # distinct metabolite names
metabolite_names(m)        # same as metabolites, as a tuple of Symbols
reactions(m)               # reaction steps as tuples of (lhs, rhs)
n_states(m)                # number of enzyme states
n_steps(m)                 # number of mechanism steps
graph(m)                   # (SimpleDiGraph, enzyme_forms)
stoich_matrix(m)           # metabolites x steps matrix
parameters(m)              # independent k's + Keq + E_total
parameters(m, RAW)         # all 2N k's + E_total
fitted_params(m)           # independent k's only (for fitting)
rate_equation_string(m)    # human-readable rate equation
```

## Identifiability

```julia
structural_identifiability_deficit(m)  # deficit (<=0 means identifiable)
is_identifiable(m)                     # true if deficit <= 0
```

## Parameter Fitting

```julia
fp = FittingProblem(m, data_table; Keq=5.0)
result = fit_rate_equation(fp, optimizer; n_restarts=10, maxtime=60.0)
```

The data table must have columns `Article`, `Fig`, `Rate`, and one column per
metabolite in `metabolites(m)`. Fitting operates in log-space on the
independent rate constants returned by `fitted_params(m)`.

## API Reference

### Types

| Type | Description |
|------|-------------|
| `EnzymeReaction{S,P,R}` | Overall reaction specification (substrates, products, regulators encoded in type parameters). |
| `EnzymeMechanism{Species,Reactions,EquilibriumSteps,ParamConstraints}` | Full mechanism with species, elementary steps, RE/SS flags, and parameter constraints encoded in type parameters. |
| `MechanismSpec` | Lightweight runtime description of a mechanism. Convert to `EnzymeMechanism` via `EnzymeMechanism(spec)`. |
| `FittingProblem` | Wraps a mechanism + data table for parameter fitting. |
| `SiteState` | State of a single binding site (metabolite, atoms, role). |
| `EnzymeFormSpec` | Specification of an enzyme form with named binding sites. |

### Macros

| Macro | Description |
|-------|-------------|
| `@enzyme_reaction` | Create an `EnzymeReaction` from a DSL block. |
| `@enzyme_mechanism` | Create an `EnzymeMechanism` from species + steps DSL blocks. |

### Functions

| Function | Description |
|----------|-------------|
| `rate_equation(m, params, concs)` | Compiled QSSA rate equation. Zero allocations. |
| `rate_equation_string(m)` | Human-readable rate equation string. |
| `parameters(m)` | Parameter names for the default (Haldane/Wegscheider) mode. |
| `parameters(m, RAW)` | All raw rate constant names + `E_total`. |
| `fitted_params(m)` | Independent rate constant names (excludes `Keq`, `E_total`). |
| `enumerate_mechanisms(rxn; max_forms)` | Lazy iterator over all valid mechanisms for a reaction. |
| `enumerate_mechanism_stages(rxn; max_forms)` | Run enumeration pipeline, returning intermediate results at each stage. |
| `enumerate_enzyme_forms(rxn)` | Enumerate all possible enzyme forms for a reaction. |
| `n_sites(rxn)` | Total number of binding sites for a reaction. |
| `is_identifiable(m)` | Whether the mechanism is structurally identifiable. |
| `structural_identifiability_deficit(m)` | Identifiability deficit (non-positive = identifiable). |
| `FittingProblem(m, table; Keq)` | Construct a fitting problem from mechanism + data. |
| `fit_rate_equation(fp, optimizer; ...)` | Fit rate constants via multi-start optimization. |
| `substrates(m)` | Substrates (with stoichiometric multiplicity). |
| `products(m)` | Products (with stoichiometric multiplicity). |
| `regulators(m)` | Regulators. |
| `enzyme_forms(m)` | Distinct enzyme states. |
| `metabolites(m)` | Distinct metabolites. |
| `metabolite_names(m)` | Metabolite names as a tuple of Symbols. |
| `reactions(m)` | Reaction steps as `(lhs, rhs)` tuples. |
| `n_states(m)` | Number of enzyme states. |
| `n_steps(m)` | Number of mechanism steps. |
| `graph(m)` | Enzyme-form connectivity graph. |
| `stoich_matrix(m)` | Stoichiometry matrix (metabolites x steps). |

## Running Tests

```julia
] test EnzymeRates
```
