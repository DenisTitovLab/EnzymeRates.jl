# Ping-pong mechanisms

In a ping-pong (double-displacement) mechanism, the enzyme carries a covalent
fragment between two half-reactions. Substrate A binds and transfers part of its
atoms to the enzyme; the modified enzyme then binds substrate B and completes the
transfer to form product Q. The enzyme returns to its original form only at the
end of the full cycle.

## The `Residual` type

The package represents the covalent fragment with a `Residual`:

```julia
struct Residual
    added::Vector{Substrate}     # atoms gained from consumed substrates
    subtracted::Vector{Product}  # atoms committed to released products
end
```

A non-empty `Residual` means atoms remain on the enzyme between the producing
and consuming steps. `has_residual(s::Species)` returns `true` when the
residual is non-empty. The rendered name of a residual-bearing form encodes both
the conformation and the residue: for example, `:E_res_+A_-Q` names an enzyme
on conformation `:E` that has gained the atoms of substrate `A` and is committed
to releasing product `Q`.

## Enumerated ping-pong intermediates

When the package enumerates mechanisms via `init_mechanisms`, every covalent
intermediate lives on conformation `:E` carrying a `Residual` — never on a
separate conformation label. The enumerator builds every enzyme form on conformation `:E` and computes each
form's residual from the consumed-substrate and released-product history. The
backtracking engine rejects a degenerate "ping-pong" step that would return
the enzyme to apo `E` with an empty residual mid-cycle: such a step would
split the reaction into two disconnected half-cycles, which is not a valid
mechanism.

This `:E`-only invariant is the **enumerator's** convention. The DSL is more
permissive: hand-written mechanisms via [`@enzyme_mechanism`](@ref) may use a
separate conformation label (for example, `:Estar`) for a covalent intermediate
provided the `residual:` field is supplied. The enumerator's output never uses
separate conformations for residual-bearing forms.

## Atom inventories

The `Residual` is computed by atom bookkeeping, so substrates and products must
declare real atom counts in the [`@enzyme_reaction`](@ref) block using the
`[atoms]` bracket syntax. For example, `A[C2N1]` declares substrate `A` with
two carbons and one nitrogen. The residual at an intermediate form is the
multiset difference between the atoms consumed from solution and the atoms
currently bound or already released.

## A concrete example

The following reaction transfers a C2N1 group from substrate `A` to substrate
`B` via the enzyme:

```@example pingpong
using EnzymeRates
rxn = @enzyme_reaction begin
    substrates: A[C2N1], B[C1]
    products:   P[C2], Q[C1N1]
end
mechs = EnzymeRates.init_mechanisms(rxn)
pp = nothing
for m in mechs, grp in EnzymeRates.steps(m), s in grp
    for sp in (EnzymeRates.from_species(s), EnzymeRates.to_species(s))
        EnzymeRates.has_residual(sp) && (global pp = m)
    end
end
residual_forms = Symbol[]
for grp in EnzymeRates.steps(pp), s in grp
    for sp in (EnzymeRates.from_species(s), EnzymeRates.to_species(s))
        EnzymeRates.has_residual(sp) &&
            (EnzymeRates.name(sp) in residual_forms ||
             push!(residual_forms, EnzymeRates.name(sp)))
    end
end
residual_forms
```

The three residual-bearing forms all sit on conformation `:E`:
`E_res_+A_-Q` (the free enzyme carrying the residual), `EB_res_+A_-Q`
(substrate `B` additionally bound — the moment just before the second
half-reaction), and `EQ_res_+A_-Q` (product `Q` bound). The `+A_-Q` suffix
records that the enzyme carries the atoms of substrate `A` and is committed to
releasing product `Q`.

Note that `init_mechanisms`, `steps`, `from_species`, `has_residual`, and `name`
are reached as `EnzymeRates.<name>` because they are internal-but-usable
functions, not part of the exported public API.

For how the enumerator generates these forms and which catalytic topologies
contain a ping-pong step, see the mechanism enumeration page.
