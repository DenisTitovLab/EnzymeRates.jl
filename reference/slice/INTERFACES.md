# Slice interfaces

The contract between the stages of the vertical slice (`module Slice`, sources in
`slice/src/`). Each stage documents the interfaces it adds or changes under its own heading,
and only edits its own files. Below, `S` is
`/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad`.

## Loading

```julia
include("S/proto/common/Population.jl")   # Main.Pop, only to read the probes' Pop records
include("S/slice/src/Slice.jl")           # Main.Slice: core files, then later files present
using Serialization
conv = deserialize("S/slice/core/data/init_bibi.jls")::Vector{Slice.Converted}
```

- `Slice.jl` includes `types.jl`, `names.jl`, `poly.jl`, `convert.jl`, then each of
  `thermo.jl`, `api.jl`, `engine.jl`, `normalize.jl`, `print.jl`, `kcat.jl`,
  `allosteric.jl`, `admission.jl`, `print_allo.jl`, `codegen.jl`, `enum.jl`, `key.jl`,
  `pipeline.jl` that exists, in that order. A later file adds its own `using` lines (stdlib
  only); it may call anything an earlier file defines.
- Serialized data names its types `Main.Slice.*`: include `Slice.jl` into `Main` before
  deserializing.

## Core (types.jl, names.jl, poly.jl, convert.jl)

### The mechanism value type (decision 1)

```julia
struct Met;  name::Symbol; modifier::Bool; end          # modifier: binds without reacting
struct Form; bound::UInt128; res::UInt128; conf::Int32; end          # isbits, 48 bytes
struct Edge; a::Int32; b::Int32; cons::UInt128; rel::UInt128; re::Bool; end   # isbits
struct Site; ligands::Vector{Symbol}; states::Vector{Symbol}; multiplicity::Int; end
struct Mechanism
    substrates::Vector{Symbol}; products::Vector{Symbol}   # sorted
    mets::Vector{Met}          # sorted table; a reactant's modifier copy is its own entry
    confs::Vector{Symbol}      # sorted conformation table
    forms::Vector{Form}        # sorted by formkey: bound count, residual count, content
    edges::Vector{Edge}        # sorted by (a, b, cons, rel, re)
    group::Vector{Int32}       # kinetic group of each edge, numbered by first edge
    tags::Vector{Symbol}       # one per group (:OnlyA, :EqualAI, :NonequalAI); empty: plain
    multiplicity::Int          # catalytic multiplicity (1 for a plain mechanism)
    sites::Vector{Site}        # sorted; ligands sorted within a site
end
```

- **Packs.** `bound`, `res`, `cons`, `rel` count table entries, four bits per entry
  (`count_at(p, i)`, `nibsum(p)`, `mets_of(m, p)`). Caps: 32 table entries and 15 copies of
  one entry per form or step side, with errors that say so.
- **Residual.** `res` counts the substrates added and the products subtracted; the role
  gives the sign.
- **Form order.** `formkey(f) = (nibsum(bound), nibsum(res), bound, res, conf)`, so free
  forms (no bound metabolite, no residual) come first, and each rapid-equilibrium segment's
  lowest form holds the fewest bound metabolites. `findfirst(is_free, m.forms)` is the first
  free form by content.
- **Edge storage.** `a < b`, except a binding-like edge (one metabolite in all, on one
  side), which is stored taking its metabolite up: `is_binding(e) = e.rel == 0 &&
  nibsum(e.cons) == 1` holds for every binding.
- **Identity.** `==` and `hash` compare every field and are defined explicitly for every
  struct, so hashes are stable across processes (checked on 4,917 mechanisms).
- **Derived, not stored:** `inhibitors(m)` (names that bind as modifiers), `regulators(m)`
  (site ligands), `metabolites(m)` (substrates, products, then every other name sorted: the
  concentration order).
- **Other accessors:** `ngroups`, `members(m, g)`, `is_allosteric`, `is_reactant(m, x)`,
  `is_dual(m, met)` (a reactant's modifier copy), `is_free(form)`, `is_plain_binding(m, e)`
  (the bound form holds the free form's metabolites plus the bound one, same residual; the
  conformation may change).

**Construction** (structure only):

```julia
m, edge_of, flipped = canonicalize(substrates, products, steps::Vector{RawStep},
                                   group::Vector{Int}; tags, multiplicity, sites)
m = Mechanism(substrates, products, steps, group; tags, multiplicity, sites)
RawStep(from::RawForm, to::RawForm, consumed::Vector{Met}, released::Vector{Met}, re)
RawForm(conf, bound::Vector{Met}, added::Vector{Symbol}, subtracted::Vector{Symbol})
(substrates, products, steps, group), kw = raw(m)      # back to raw input (moves)
```

`edge_of[j]` is the edge that step `j` became; `flipped[j]` says the edge runs against the
step as written. The constructor rejects: missing or overlapping substrates and products;
a residual that adds a product or subtracts a substrate; a non-reactant name not flagged
`modifier`; the caps; a self-loop; one metabolite on both sides of a step; one reaction in
two steps; group labels not `1:G` or a group without a step; tags not one per group; sites
or multiplicity on a plain mechanism; a malformed site or a ligand at two sites; group
members that, in their physical-forward direction, take up or give off different
metabolites or differ in RE/SS; group members that make different conformational
changes; steps that do not form one connected graph. The default field constructor
`Mechanism(subs, prods, mets, confs, forms, edges, group, tags, multiplicity, sites)` is
trusted (used by `decode`); call it only with canonical fields.

**Boundary checks** (separate functions, for the DSL and enumeration):
- `check_boundary(m)` errors on a bottomless rapid-equilibrium segment (`bottomless(m,
  mask)` returns its forms), a binding that changes a metabolite's role (`role_change(m)`,
  `E + P::Mod → E(P)`), and colliding rendered names (`check_names(m)`).
- `check_atoms(m, atoms)`: every step conserves atoms; `atoms[x]` lists `element => n`.
- `:OnlyA` admission is the admission stage's.

### Orientation and kinetic groups (decisions 1 and 2)

- `orient(m) -> fwd::BitVector`; `orientation(m) -> (fwd, tier)`. Today's tiers: a binding
  takes its metabolite up (0); more substrates, then fewer products, bound or free, first
  (1); from a product-only form to a substrate-only form (2); lower form index first (3).
  Tied members of a kinetic group turn with the group's lead.
- `oriented(m, k, fwd)` is edge `k` in its physical-forward direction: `o.a` plus `o.cons`
  becomes `o.b` plus `o.rel`.
- Group alignment is forward with forward: the constructor has checked that every member
  of a group, in its physical-forward direction, takes up and gives off the same
  metabolites with one RE/SS flag. Ties are therefore slot by slot (below); the solve never
  re-derives alignment.
- On all 38,263 package-built steps the physical-forward direction equals the package's
  stored direction.

### Bare per-step constants

Edge `k` carries the constants `constant_names(m)[k]` names, in that order. With
`o = oriented(m, k, orient(m))`, and concentrations read by name (a modifier copy reads its
reactant's concentration):

| Edge | Slot 1 | Slot 2 |
|---|---|---|
| steady state | `kf`: rate `kf·[o.a]·Π c(o.cons)` of `o.a → o.b` | `kr`: rate `kr·[o.b]·Π c(o.rel)` of `o.b → o.a` |
| RE binding (`is_binding(o)`) | `K`, the dissociation constant: `[o.a]·c(M) = K·[o.b]` | – |
| any other RE edge | `K` of the forward direction: `[o.b]·Π c(o.rel) = K·[o.a]·Π c(o.cons)` | – |

Thermodynamically, the forward equilibrium ratio `exp(G[o.a] + Σμ(o.cons) − G[o.b] −
Σμ(o.rel))` equals `kf/kr`, `1/K` (binding) or `K` (other RE); checked against the oracle's
consistent constants on all 38,449 edges (worst 3.1e-154). Slot 2 (`kr`) is the constant
that drives against the physical-forward direction (the solve priority's `against`).

**Columns and classes for the joint solve (decision 4).** A bare constant is
`(k, slot, inactive::Bool)`. Active state (and plain mechanisms): the members of a group tie
slot by slot. Inactive copies exist only for allosteric mechanisms: none for an `:EqualAI`
group (its constant is the active one), one class per `:NonequalAI` group (members tied),
and one class per member of an `:OnlyA` group (untied).

**Priority inputs** the type provides: kind (isomerization `o.cons == o.rel == 0`; plain
binding `is_plain_binding(m, o)`; else other chemistry), occupied (no form of the step has
`bound == 0`), ss (`!o.re`), against (slot 2). Content names for ties: `step_text(m, o)`
(DSL text from content), which is invariant under storage order, direction and grouping,
not under renaming a metabolite.

### Names (decision 11)

- `constant_names(m; inactive = false) -> names[k]`: per edge, its constants' names in slot
  order. Every edge renders its group's naming representative (`naming_reps(m, fwd)`,
  today's rule: free-enzyme metabolite steps first, isomerizations last, ties by rendered
  form names); in the inactive state an `:OnlyA` edge renders itself. State tags `A_`/`I_`
  except for `:EqualAI` groups and plain mechanisms.
- A rate constant is named after the direction it drives (`k_EA_to_EP`), a binding's `K`
  after its release (`K_EA_to_E_A`), any other `K` after the forward direction.
- `regulator_names(m) -> [ligand => names]`: `K_A_Rreg`, `K_I_Rreg`, both, or the bare
  `K_Rreg` for an `:EqualAI` site. Fixed names: `:L`, `:Keq`, `:E_total`.
- Renderers: `form_name(m, i)`, `side_name(m, i, pack)`, `label(m, i)`; DSL text for
  messages: `form_text`, `step_text`, `met_text`.
- **Changed from today** (decided): residuals print `_res_plusA_minusP` (today
  `_res_+A_-P`); a reactant's modifier copy prints `Amod` (today `Ainh`, DSL `A::Mod`); a
  modifier-only name prints bare (`K_EI_to_E_I`, today `K_EIinh_to_E_Iinh`); `:EqualAI`
  regulators are bare (`K_Rreg`); inactive `:OnlyA` members name themselves. Nothing else
  differs from the package on 77,289 compared names.

### Rapid-equilibrium forest

`re_forest(m, mask = all edges) -> (segment, parent, order)`: segments of the RE edges in
`mask`, numbered by their root (the segment's lowest form index); `parent[i]` is `+k`
(form `i` is `edges[k].b`) or `-k` (`edges[k].a`), 0 at a root; a breadth-first tree that
scans edges in index order (content order, never storage or written order). `order`
lists forms root first, each after its parent. Going `a → b` along edge `k`, form `b`'s
weight is form `a`'s times its constant factor (table above) and `Π c(cons) / Π c(rel)`.

### The Sig (decision 1)

```julia
struct EnzymeMechanism{S} end        # S::NTuple{N, Int}, the only type parameter
S = sig(m)                           # one homogeneous integer tuple
m = decode(S...)                     # inside a generator: decode by splatting
em = compile(m)                      # EnzymeMechanism{sig(m)}()
m = mechanism(em)                    # reads S at run time; compiles once
```

`decode(sig(m)...) == m` for all 4,917 mechanisms. A `@generated` function over
`EnzymeMechanism{S}` that calls `decode(S...)` compiles exactly one method instance per
mechanism (its own); `decode` and the generator compile once for every Sig length. Equal
mechanisms give the identical type.

### poly.jl

```julia
t = Vars(); v = var!(t, :A; conc = true); k = var!(t, :k_EA_to_E_A)
const Exp = Rational{Int}                 # integer for concentrations, rational for constants
const Mono = Vector{Pair{Int, Exp}}       # sorted by variable, no zero exponent
struct Poly{T <: Integer}; terms::Dict{Mono, T}; end   # T = Int (checked) or BigInt
```

Operations: `mono(v => e, ...)`, `mmul`, `mpow(m, e)`, `mexp(m, v)`; `term(T, m, c)`,
`one`, `zero`, `+`, `-`, `*` (by a polynomial or a monomial), `scale(p, c)`, `addto!`,
`nterms`, `iszero`, `==`, `hash`; `convert(Poly{BigInt}, p)` (Int arithmetic throws
`OverflowError`: retry in BigInt); `subst(p, v -> mono or nothing)` (power products,
rational exponents); `evaluate(p, vals)` in Float64 or BigFloat (`pow` takes half powers
through `sqrt`); `divexact(p, d)` for a two-term `d` (the Haldane factor; `ArgumentError`
when it does not divide); `conc_shift((num, den), t)` (the concentration monomial that
brings both to lowest terms: multiply both by it); `split_conc(p, t)` (`concentration
monomial => coefficient polynomial`).

### convert.jl and the converted data

- `Converted(id, source, index, mech, edge_of, flipped)`: `index` is the record's position
  in its Pop file; `edge_of`/`flipped` map each Pop step to its edge.
- `Converted(p, index)` converts a `Pop.PMech` (read through its fields); `pop_raw(p)`
  gives the raw input; `canonicalize(raw::NamedTuple)`.
- `pop_constants(c, x)`: per-edge constants (`constant_names` order) to per Pop step
  constants in the oracle's convention (`Oracle.steady_state_v(p, k, concs)`; per state for
  `MWC.jl`); `edge_constants(c, k)` the other way. Package-built records store each step in
  its physical-forward direction; hand-built ones may not, and the map handles both.
- `S/slice/core/data/<set>.jls` holds `Vector{Slice.Converted}`; `MANIFEST.tsv` lists each
  set, its counts and its Pop source file (`deserialize(source)[c.index]` is the record);
  `<set>.failed.txt` lists what does not convert. `pop_p1_direct.jls` and `pop_p2_hand.jls`
  are the Pop records of the hand-built sets.
- Sets: specs 42, handwritten 34, init_uniuni 1, init_unibi 7, init_bibi 239, init_bibi_pp
  264, expand_bibi 700, allo_uniuni 35, allo_bibi 1477; allo_uniuni_inh 123, allo_unibi 133,
  allo_bibi_pp 400, allo_terbi 192; p5a_adv_cases 22, p5a_adv_cases2 12, p5a_adv_cases3 3,
  p5a_adv_cases4 1, p5a_cases 116; p1_adv_pkg 25, p1_adv_pkg2 2, p1_enum_inh 1041, p1_direct
  5 (D1–D5); p2_adv_pkg 22, p2_adv_pkg2 3, p2_fuzz_dump 5, p2_fuzz_dump_flip 10, p2_hand 3
  (2 more are rejected by the constructor, as the package rejects them). 4,917 in all.
- Ids repeat across sets (`p1_adv_pkg2` reuses `adv_pkg/1`): key records by
  `(set, index)`.
- Records without package outputs: p1_direct, p2_hand, p2_fuzz_dump, p2_fuzz_dump_flip
  (and p1_enum_inh has names but no polynomials or samples).
- Rejected cases (`ok = false`) in the P5a sets were exported from a `:NonequalAI` twin:
  their Pop names and outputs are the twin's. `p5a_cases` Pop records of rejected grid
  cases carry the twin's tags; the converted mechanisms carry the restored ones
  (`S/slice/core/retag.jl`).

## Expected outputs of the later stages (core's proposal; each stage documents its own)

The shapes below are proposals so that stages can work in parallel; a stage that changes
one documents the change under its own heading.

### thermo.jl: the joint solve (decision 4)

```julia
struct Col; edge::Int; slot::Int; inactive::Bool; end       # a bare constant
struct Dep
    col::Int                                  # the dependent column
    terms::Vector{Pair{Int, Rational{Int}}}   # independent column => exponent
    keq::Rational{Int}                        # exponent of Keq
    label::Symbol                             # :kinetic_group, :haldane or :wegscheider
end
struct Solution
    cols::Vector{Col}        # every column (active, then inactive copies)
    indep::Vector{Int}       # independent columns, in parameter order
    deps::Vector{Dep}        # every other column, as a power product of independents and Keq
    fold::Dict{Int, Int}     # dependent => independent where terms == [i => 1], keq == 0;
                             # one table for both conformations
end
solve(m) -> Solution
```

Rows: fundamental cycles of a spanning forest (Haldane rows with integer multiples of
`log Keq`, Wegscheider rows), kinetic-group equalities, `:EqualAI` ties, RE-only cycles.
Priority: `(kind, occupied, ss, against)`, then the content name (`step_text` of the
oriented edge with its slot). Exact RREF over checked Int64, retried in BigInt on
overflow. A column's name is `constant_names(m; inactive)[edge][slot]`.

### api.jl (decision 5)
`parameters(m)` = the independents' names, then `:Keq`, `:E_total`; `fitted_params(m)` =
the independents (regulator constants and `:L` included for allosteric mechanisms);
`all_parameters(m, params)` = every class's constant by name, dependents computed.

### engine.jl: the law (decision 3)

```julia
struct Law
    vars::Vars          # concentrations (metabolites(m) order, conc = true), independents
                        # by name, :Keq
    H::Poly{Int}        # Π c(S) − Π c(P)·Keq^-1
    Q::Poly             # positive quotient: v = E_total·H·Q/D
    D::Poly             # Q and D multiplied by conc_shift((Q, D), vars)
end
```

The bare law (before the solve) has one variable per column (internal names, such as
`Symbol("c", col)`): `num`, `den` by Cha contraction on `re_forest`, whole-graph spanning
trees oriented per root, the edge-sum numerator, and a cap on the tree count. Substituting
the dependents (`subst`, rational exponents) and dividing by `H` (`divexact`) gives the
`Law`. For codegen the engine also exposes its node structure (each RE segment's weights
and sum, `H`).

### normalize.jl, print.jl (decision 8)
Printed law `(ΠS − ΠP/Keq)·Q / D`, numerator and denominator divided by the constant-term
coefficient `d0`, or by a virtual constant term `N` (Cleland's, then the general
construction, then plain coefficient division); the coefficient form with one legend line
per coefficient. Coefficients from `split_conc`; names from the chokepoint.

### kcat.jl (decision 7)
The saturating limit: products at 0, substrates to infinity, the ratio of the leading
coefficients of `Q` and `D` (`split_conc`), on the mechanism without its dead-end
branches; the maximum over regulator corners; an error naming the problem when nothing
saturates.

### allosteric.jl, admission.jl, print_allo.jl (decision 6)

```julia
struct AlloLaw
    vars::Vars                       # as Law, plus :L and the regulator constants
    H::Poly{Int}
    QA::Poly; DA::Poly               # active state, normalized by the flip form's weight
    QI::Poly; DI::Poly               # inactive state (edge mask pruned from the flip form)
    sites::Vector{Tuple{Poly, Poly, Int}}   # per site: W_A base, W_I base, multiplicity
    n::Int                           # catalytic multiplicity
end
# v = E_t·H·(QA·DA^(n-1)·WA + L·QI·DI^(n-1)·WI) / (DA^n·WA + L·DI^n·WI),
# WA = Π_sites (W_A base)^multiplicity, likewise WI.
admit(m) -> nothing or the rejection message    # G1 cone test and island rule
```

The flip form is chosen by content (the first free form, `findfirst(is_free, m.forms)`,
is the natural candidate); another free form is allowed only in its RE segment. Names and
orientation always come from the full mechanism.

### codegen.jl (decision 9)
`rate_equation(em::EnzymeMechanism{S}, concs, params)`: the generator decodes
`decode(S...)`, derives the `Law` or `AlloLaw` (no global memo), loads positionally in
`parameters(m)` and `metabolites(m)` order, computes each dependent and each reciprocal
once, keeps node locals (`H`, segment sums and weights, MWC state, site and cross-weight
terms) and emits greedy Horner with binary `+` and `*`; `sqrt` for half powers.

### key.jl (decision 10)
Input: the law's components in a fixed role order (`Q`, `D`, `H`; per conformation for an
`AlloLaw`, then the site bases with multiplicities), the fixed variables (concentrations,
`Keq`, `L`) and the parameters. Output: the SHA-256 hex of the canonical text, parameters
labelled by colour refinement, computed before code generation.

### enum.jl, pipeline.jl
Moves rebuild children with `raw`/`canonicalize` and filter them with `check_boundary`,
`check_atoms` and `admit`. The pipeline chains type, solve, engine, normalization, codegen
and key per mechanism.

## Solve (thermo.jl, api.jl)

These replace the core's proposal for `thermo.jl` and `api.jl` above where they differ.

### The joint solve (decision 4)

```julia
struct Col; edge::Int; slot::Int; inactive::Bool; end        # a bare constant
struct Dep
    col::Int                             # the dependent column
    terms::Vector{Pair{Int, Exp}}        # independent column => exponent (Rational{Int})
    keq::Exp                             # exponent of Keq
    label::Symbol                        # :kinetic_group, :haldane or :wegscheider
end
struct Solution
    cols::Vector{Col}        # every column, in elimination order
    index::Dict{Col, Int}    # every bare constant => its column; in an allosteric
                             # mechanism an :EqualAI edge's inactive one => its active one
    indep::Vector{Int}       # independents of the law, in parameter order
    deps::Vector{Dep}        # dependents of the law, in pivot order
    fold::Dict{Int, Int}     # dependent => independent when its right side is that one
end
sol = solve(m)
mask = inactive_mask(m)                # edges of the inactive law (BitVector)
cols, A, b = constraint_system(m)      # A·log(constants) = b·log(Keq), elimination order
piv = rref!(M, n)                      # in place on Matrix{Rational{T}}, first n columns
zs, multiples = cycles(m, fwd)         # fundamental cycles (physical-forward), Keq multiples
class(m, c)                            # kinetic class of a column (ties and labels)
```

- **Columns.** Active: one per `(edge, slot)`. Inactive (allosteric only): one per
  `(edge, slot)` of every edge outside `:EqualAI` groups; each `:OnlyA` member is its own
  class there, `:NonequalAI` members are tied by rows. An `:EqualAI` edge has no inactive
  column: its inactive constant is its active column (the `:EqualAI` tie, implicit; `index`
  maps it).
- **Rows.** Each fundamental cycle of a spanning tree of the whole graph, in the active state
  and, for an allosteric mechanism, again over the inactive constants (whole graph, `:OnlyA`
  edges included); then each class member equal to its class's first column. A cycle row sums
  `z·log(forward equilibrium ratio)` (`kf/kr`, `1/K` for a binding, `K` otherwise) and its right
  side is the cycle's Keq multiple. **RE-only cycle rows are not separate rows**: every RE-only
  cycle lies in the span of the fundamental cycles (adding them changed the rank in 0 records).
- **Elimination order.** Inactive columns outside the law first, then the other inactive
  columns, then the active ones; within a block the key `(kind, occupied, ss, against)` larger
  first, ties by the oriented edge's content `(o.a, o.b, o.cons, o.rel)` larger first (form
  indices are content-ordered, so this is a content key with the invariances of a content name).
  **`against` is not "slot 2"**: it marks the rate constant that drives against the net
  reaction, so for a binding of a product (physical-forward = taking it up) it is slot 1, the
  binding. With "slot 2" the Haldane dependent becomes the catalytic product release
  (Uni-Uni: `k_ES_to_E_P`), the reverse of P2's design.
- **Results.** Only columns of the law enter `indep`, `deps` and `fold`: columns outside it
  (inactive `:OnlyA` columns, pruned islands) are eliminated first, so no kept right side holds
  one, and a free one is not fitted. Labels: `:kinetic_group` when the dependent is not its
  class's last (kept) column, else `:haldane` when `keq != 0`, else `:wegscheider`.
- **Arithmetic.** Checked `Rational{Int}`, retried in `Rational{BigInt}` on `OverflowError`
  (exponents converted back to `Rational{Int}`). Errors: a cycle whose net change is no
  multiple of the net reaction (the package rejects the same, `p1_adv_pkg` A9 slip cycle);
  contradictory relations (`0 = log Keq`).
- **For the cone test (admission).** `constraint_system(m)` holds every column, the law's and
  the others; the inactive `:OnlyA` columns are `Col(k, s, true)` with
  `m.tags[m.group[k]] === :OnlyA`. An RE `:OnlyA` chemistry edge has one `K` column (its
  equilibrium ratio, sign-free under the kinetic-block reading); P5a's steady-state copy enters
  the rows only through its ratio, so feasibility is the same.
- **`inactive_mask(m)`**: the edges outside `:OnlyA` groups that the flip form (the first free
  form by content) reaches through such edges; all false for a plain mechanism; an error when
  an allosteric mechanism has no free form. It does not check decision 6's rule that another
  free form lies in the flip form's RE segment (admission's job).

### The parameter API (decision 5)

```julia
column_names(m, sol) -> Vector{Symbol}      # each column's constant_names in its state
fitted_params(m) -> Tuple    # the law's independents, then regulator constants, then :L
parameters(m) -> Tuple       # (fitted_params(m)..., :Keq, :E_total)
all_parameters(m, params::NamedTuple) -> NamedTuple   # params plus each dependent name once
```

- **Parameter order** (content order, not today's alphabetical one): active independents by
  `(edge, slot)`, inactive independents by `(edge, slot)`, regulator constants in
  `regulator_names` order, `:L`.
- `all_parameters` reads every name of `parameters(m)` from `params` (`:Keq` included) and adds
  each dependent class's name, `Keq^keq·Π x^e` evaluated with `pow` (half powers by `sqrt`) in
  the values' type (Float64 or BigFloat). A per-step constant reads its value by its
  `constant_names` name (an `:EqualAI` edge's inactive constant by its active name).
- **For the engine:** edge `k`'s slot `s` in a state is column `sol.index[Col(k, s, inactive)]`;
  substitute each `Dep` (`terms`, `keq`) with `subst`; an independent `i` is named
  `column_names(m, sol)[i]`.

### Notes for integration

- `naming_reps` (today's rule) and the solve's kept class member can differ. The solve is
  identical under same-model regroupings, but in 57 of 12,608 same-model splits a fitted
  *name* changes (init_bibi_pp/252: `K_EQ_to_E_res_plusA_minusP_B` becomes
  `K_EBP_res_plusA_minusP_to_EP_res_plusA_minusP_B`). Naming a class by its kept member would
  remove this but changes today's names (decision 11).
- `enum.jl` (parallel stage) already calls `solve` and `inactive_mask`; both are stable.
- Checks (a)–(f): scripts in `S/slice/solve/` (`run_all.sh` reruns them under the lock,
  about 17 minutes), logs in `S/slice/solve/logs/`.

## Enum (enum.jl)

The three ported moves. Each takes a parent and returns its children as `Mechanism`s, in no
particular order; a parent must satisfy the package's emission rules (every steady-state
group carries flux, no redundant copy group), as `expand_mechanisms` asserts today.

```julia
flip_children(m) -> Vector{Mechanism}                 # _expand_re_to_ss
split_children(m) -> Vector{Mechanism}                # _expand_split_kinetic_group
dead_end_children(m, declared::Vector{Symbol}) -> Vector{Mechanism}
                    # _expand_add_dead_end_regulator; declared = the reaction's dead-end
                    # inhibitor names (a substrate or product name binds as its ::Mod copy)
```

Shared helpers a later file may call:

```julia
regroup(m, labels, tags, re)   # m with new group labels (renumbered), tags and RE flags,
with_re(m, re)                 # through the trusted constructor: a flag or label never
                               # changes form or edge order; each group keeps one flag
minimal_sets(gains, partners)  # Apriori search of inclusion-minimal gaining unit sets
all_subsets(v)                 # every subset, the empty one first
segment_offsets(m, mask) -> (segment, names, ex)  # re_forest segments; ex[name, form] =
                               # copies held above the segment's fewest (a ::Mod copy by name)
edge_flux(m), group_flux(m)    # net steady-state flux per edge / per group (Tarjan blocks)
unbalanced_blocks(nv, ends, weights) -> (block, unbalanced)
chain_flanks(m)                # groups the flip skips (chain lemma; none if allosteric)
hyperbolic(m), reaches_all(n, arcs, root, fixed)  # degree ≤ 1 test of the flip (n > 1)
modifier_mask(m), taken(m, o), copy_entry(m, o)
twin_finder(m, mask), copy_gauge(m, fwd, mask, g, twin, label), redundant_copies(m)
bipartitions(m, fwd, g)        # the split's context divisions of group g
```

- **Dependencies.** The split's gain test is `length(solve(c).indep) > length(solve(m).indep)`
  (thermo.jl); the copy rule's inactive graph is `inactive_mask(m)` (thermo.jl). Children of
  the flip and the split come from `regroup`/`with_re`; dead-end children add forms and a table
  entry, so they go through `raw` and the constructor.
- **Same children as the package** on every plain parent (1,104: init_bibi_pp 264, dead1 300,
  expand_bibi 240, p1_enum_inh 300) for all three moves, and on every allosteric parent
  (760) for the flip and the dead-end move. The split differs on 143 of 760 allosteric
  parents; every difference is a gain verdict where the package's allosteric count exceeds
  the joint solve's (P5b's count agrees with the solve on all of them): ROADMAP known
  issues 1 and 2. Logs: `slice/enum/logs/compare_final.log`, `diag_split_b.log`.
- **Deliberate differences from the package code:** the copy rule runs on every copy group of
  a dead-end child (the package tests only the new group of a plain child; no difference on
  the data); a complex with several twins and none of its composition takes the first by
  content, not by rendered name (no difference on the data).

### Notes for integration (enum)

1. **thermo.jl `cycles` is slow by a Julia pitfall.** `isassigned(path, e.a)` with
   `e.a::Int32` takes Base's exception-based `isassigned(::AbstractArray, ::Integer...)`:
   40 µs per unassigned check against 6 ns with an `Int` index (`enum/isassigned_check.jl`).
   `cycles` then takes 250–340 µs on 7–10-edge mechanisms, about 95% of `solve`.
   `isassigned(path, Int(e.a))` (and `e.b`) fixes it. The split calls `solve` per candidate:
   with the fix it runs at 0.8–1.0× the package's speed on plain parents and 6–12× faster on
   allosteric ones (same children on all 1,864 parents); without it, 0.1–0.6×.
2. `bottomless` (types.jl) repeats `segment_offsets`'s exponent walk; it can be written on
   `segment_offsets` (about 10 code lines less).
3. **Admission.** The moves do not call `admit` (admission.jl did not exist). A flip or a
   split changes neither the `:OnlyA` members nor the inactive mask, so by the definitions it
   should not change the G1/island verdict (inference, untested); the dead-end move adds
   edges and forms and can, so the driver (or `dead_end_children`) must filter dead-end
   children of an allosteric parent with `admit`.
4. **Atoms and names.** The type has no reaction, so the driver runs `check_atoms(c, atoms)`
   and `check_names` (part of `check_boundary`) on children; every child in the test sets
   passed both, and `canonicalize(raw(c)...)` reproduced every child.
5. **Re-running the checks.** `slice/enum/pkg_moves.jl <set>` (er_oracle project) saves the
   package's children to `slice/enum/data/pkg_<set>.jls`; `slice/enum/compare.jl` (stdlib,
   loads only types/names/poly/convert/thermo/enum through `enum/load.jl`) diffs them;
   `coverage.jl` switches each filter off in turn; `diag_split_{a,b}.jl` with
   `pkg_counts.jl`/`pkg_issues.jl` attribute the allosteric split differences.
6. Top-level names enum.jl defines (a later file must not reuse them): `regroup`, `with_re`,
   `all_subsets`, `minimal_sets`, `modifier_mask`, `segment_offsets`, `unbalanced_blocks`,
   `edge_flux`, `group_flux`, `chain_flanks`, `reaches_all`, `hyperbolic`, `taken`,
   `copy_entry`, `twin_finder`, `copy_gauge`, `redundant_copies`, `flip_children`,
   `bipartitions`, `split_children`, `dead_end_children`.

## Engine (engine.jl, normalize.jl, print.jl, kcat.jl)

These replace the core's proposal for `engine.jl`, `normalize.jl`, `print.jl` and `kcat.jl`
above where they differ.

### The law (decision 3)

```julia
struct Law{T <: Integer}        # T = Int, or BigInt when Int overflows (retried in `law`)
    vars::Vars                  # concentrations (metabolites(m) order, conc = true), then
                                # every independent of the solve (sol.indep order, both
                                # conformations, named by column_names), then :Keq
    H::Poly{T}                  # ΠS − ΠP·Keq^-1
    Q::Poly{T}                  # positive coefficients; v = E_total·H·Q/D
    D::Poly{T}
    sigma::Vector{Poly{T}}      # per RE segment of the law's graph: the sum of its weights
    trees::Vector{Poly{T}}      # per segment: its spanning-tree sum
    shift::Mono                 # D == shift·Σ_g sigma[g]·trees[g]; Q carries the same shift
end
L = law(m, sol = solve(m), mask = all edges; inactive = false)
b = bare_law(m, mask = all edges, T = Int)  # (; vars, cols, num, den, sigma, trees)
complexity(m, mask = all edges) -> BigInt   # V×τ, today's _eq_complexity, for any mask
spanning_trees(n, pairs)       # errors when n × (tree count) exceeds MAX_TERMS = 5000
segment_nodes(m, mask) -> (node, G, parent, order)
conc_vars(m), conc_mono(vars, names)
```

- **Bare law.** Variables: the concentrations, then `Symbol("c", k, "_", s)` for slot `s` of
  each edge `k` in `mask` (`b.cols[j] = (k, s)` is variable `nc + j`). RE segments and their
  trees come from `re_forest(m, mask)` (lowest form, edges in content order). The numerator
  is the edge sum: `Σ_SS ω_e·(trees[a]·f_e − trees[b]·r_e)`, `ω_e` the exponent of
  `m.substrates[1]` in `f_e/r_e`.
- **Substitution.** Bare variable `(k, s)` becomes column `sol.index[Col(k, s, inactive)]`:
  an independent's variable, or a dependent's power product (rational exponents). The law is
  entirely in independents and `Keq`; dependents are not kept as names. Then `Q = num ÷ H`
  (`divexact`) and `(Q, D)` are shifted to lowest terms (`conc_shift`).
- **Masks.** `mask` must hold at least one edge. The allosteric stage calls
  `law(m, sol, inactive_mask(m); inactive = true)`; both conformations then share `vars`.
  An inactive state with no edge (the flip form alone) is the allosteric stage's case
  (`Q_I = 0`, `D_I = 1`). `kcat` calls `law(m, sol, catalytic_core(m))`.
- **Cap.** The engine expands `G` × (spanning trees of the simple segment graph) products and
  errors beyond `MAX_TERMS`; `complexity` counts parallel steady-state edges too (it is
  today's number, for the identify filter; sum it over both masks for an allosteric law).
- **For codegen.** `L.sigma`, `L.trees` and `L.shift` give the node structure; `L.H` is the
  Haldane node. The fold table and dependents are in `sol` (the law has them substituted).

### Normalization and printing (decision 8)

```julia
rule, lambda = normalizer(m, L)    # N = Π_{(mono, e) in lambda} c(mono)^e over L.D's
                                   # coefficients; rule: :constant, :cleland, :affine or
                                   # :coefficient (the fallback; terms not dimensionless)
P = printed(m, L)                  # (; rule, N, q, d): q/d are [(conc mono, (c, mo, rest, over))]
rate_equation_string(m, L = law(m), sol = solve(m)) -> String   # the full law
coefficient_form(m, L = law(m)) -> String
content(p) -> (c, mu, rest)        # integer and monomial content, rest or nothing
```

- **Full law.** `(; params...) = params`, `(; concs...) = concs`, `# Wegscheider constraints:`
  and `# Haldane constraints:` (each dependent of `sol` that is not a kinetic-group tie, as a
  power product of independents and `Keq`), `N = ...` when the normalizer has a non-monomial
  factor, then `v = E_total * (ΠS - ΠP / Keq) * (Q / N) / (D / N)`. The normalizer's integer
  and monomial contents divide every printed coefficient; a coefficient whose polynomial
  factor is `N` itself cancels to a monomial (ordered bi-bi prints `A * k_E_A_to_EA /
  k_EA_to_E_A`). The text is valid Julia: evaluated with `params` and `concs` it equals the
  law (checked on all 1,921 plain laws).
- **Coefficient form.** `v = E_total * (H) * (n_0 + n_A * A) / (1 + d_A * A + d_A²_B * A^2 * B
  + ...)`: coefficients named by the concentrations they multiply, a power as a superscript
  (`d_A²_B`; valid identifiers), then a legend line per coefficient (after the `N` line).
  Without a constant term there is no leading `1`.
- Monomials print as `k1 * k2^2 * sqrt(Keq) / (k3 * Keq)`: exponent 1 bare, integer powers
  `^n`, `1/2` as `sqrt`, other rationals `^(p/q)`; no `1 ^ 2` artifacts.

### kcat (decision 7, plain mechanisms)

```julia
catalytic_core(m) -> BitVector                 # edges outside the dead-end branches
vars, terms = kcat_terms(m, sol = solve(m), mask = catalytic_core(m))  # [(num, den)]
kcat(m, params::NamedTuple) -> the largest candidate at params
```

- Dead-end branches: the biconnected blocks without a turnover cycle, removed while one meets
  the rest at one form at most. The candidates: products at 0; for each set of modifiers
  (non-reactant names) one term holds, the other modifiers at 0; substrates grow, then that
  set's modifiers; the ratio on the denominator's leading face. Errors: the rate grows
  without bound ("never saturates"), or the face's ratio is not one number (depends on how
  the concentrations grow), or `D` vanishes; each suggests absolute mode.
- **Dependency on a later file.** `catalytic_core` calls enum.jl's `unbalanced_blocks`
  (resolved at call time). Integration should move `unbalanced_blocks` before kcat.jl
  (types.jl or a graph-helper file).

### Names engine-stage files define
`MAX_TERMS`, `Law`, `conc_vars`, `conc_mono`, `segment_nodes`, `spanning_trees`, `toward`,
`bare_law`, `law`, `complexity`, `normalizer`, `power_text`, `ratio_text`, `terms_of`,
`term_text`, `sum_text`, `content`, `conc_key`, `printed`, `group_sum`, `v_text`,
`rate_equation_string`, `SUPERSCRIPTS`, `conc_label`, `coefficient_form`,
`catalytic_core`, `kcat_terms`, `kcat`.

Checks and logs: `S/slice/engine/` (`run_all.sh` reruns them under the lock).

## Allosteric (allosteric.jl, admission.jl, print_allo.jl)

These replace the core's proposal for the allosteric files above where they differ.

### The MWC law (decision 6, formulation 1)

```julia
struct AlloLaw{T <: Integer}
    vars::Vars                    # active.vars (concentrations, independents of both states,
                                  # :Keq), then the regulator constants (regulator_names
                                  # order), then :L
    H::Poly{T}                    # ΠS − ΠP/Keq
    active::Law{T}                # law(m, sol, mask)
    inactive::Union{Nothing, Law{T}}  # law(m, sol, inactive_mask(m) .& mask; inactive = true);
                                  # nothing when the inactive state is the flip form alone
    fA::Poly{T}                   # the flip form's weight in each law over their common
    fI::Poly{T}                   # factor; both 1 when the weights are equal
    sites::Vector{Tuple{Poly{T}, Poly{T}, Int}}  # per site: 1 + Σ[l]/K in A, in I; multiplicity
    n::Int                        # catalytic multiplicity
end
A = allo_law(m, sol = solve(m), mask = all edges)   # error if no edge of mask holds forms[1]
```

    v = E_total·H·(fI^n·QA·DA^(n−1)·WA + L·fA^n·QI·DI^(n−1)·WI) / (fI^n·DA^n·WA + L·fA^n·DI^n·WI)

with `QA, DA = A.active.Q, A.active.D`, `QI, DI = A.inactive.Q, A.inactive.D` (`0` and `1` when
`A.inactive === nothing`), `WA = Π_s sites[s][1]^sites[s][3]`, `WI` likewise with `sites[s][2]`,
`L` the variable `A.vars.index[:L]`. Checked against P5a's two-conformation oracle on 2,995
allosteric records (`slice/allo/check_oracle.jl`).

- **One variable table.** `A.active.vars` and `A.inactive.vars` hold the same names in the same
  order and are a prefix of `A.vars`, so every polynomial of `A` is over `A.vars`.
- **The flip form** is `m.forms[1]`: free forms sort first, so it is the first free form by
  content, it is its rapid-equilibrium segment's root (weight 1 in the segment), and it is the
  `a` end of every edge it lies on. Its weight in a law is `L.trees[1] * L.shift`, on `L.D`'s
  scale.
- **Cross weights.** `cross_weights(fa, fi)` divides both flip weights by their common integer
  and monomial content and, when their polynomial parts are equal, by that part too. Both are
  1 in 2,003 of the 2,865 admitted laws. Concentration exponents of `fA`, `fI` are ≥ 0;
  constants may carry negative or rational exponents, as in any law.
- **Lowest terms.** Each law is in lowest terms in the concentrations, `WA`, `WI` have the
  constant term 1, and `cross_weights` leaves no common concentration factor in `fA`, `fI`, so
  the denominator has no common concentration monomial.
- **Masks.** `inactive_mask(m)` (thermo.jl) is deletion plus pruning from the flip form; the
  inactive law is the `:OnlyA` limit exactly when `admit(m) === nothing`.

### Admission (decision 6)

```julia
admit(m) -> nothing or the rejection message      # nothing for a plain mechanism
```

On the components of the forms that the edges outside `:OnlyA` groups join (the inactive
state before pruning), each `:OnlyA` plain binding (`is_plain_binding` of its physical-forward
edge) is an arc from its unbound to its bound side's component. Rejected, in this order:
1. no free form, or a free form outside the first one's rapid-equilibrium segment of the edges
   outside `:OnlyA` groups (the segment both states share, where another flip only rescales L);
2. G1: an arc inside one component, or arcs closing a directed cycle (the message names the
   bindings of the cycle);
3. islands: a component the flip form's component does not reach along the arcs.

G1 is exact: with every constant outside the inactive `:OnlyA` columns fixed, only the inactive
cycle rows move, so a direction's equilibrium-ratio exponents are differences of component
potentials; plain bindings need a positive difference, every other `:OnlyA` step (SS: both
rates → 0; RE: K free) none. An exact rational LP on `constraint_system(m)` agrees on all 2,995
records (`slice/allo/check_g1.jl`, `g1_lp.jl`). Admission takes about 1 µs.

**Enumeration.** The three ported moves keep admission (24,503 children of 1,345 admitted
parents, none rejected: `slice/allo/enum_admit.jl`); the package's tagging move
(`_expand_to_allosteric`, not ported) is what produces the 81 rejected `:OnlyA` ping-pong
mechanisms, so the driver must filter every allosteric child with `admit`.

### kcat (decision 7, allosteric)

```julia
kcat_mask(m) -> BitVector     # catalytic_core(m) plus the RE steps between free forms
vars, terms = kcat_terms(m, allo_law(m, solve(m), kcat_mask(m)))   # [(num, den)]
```

A method of `kcat_terms` (second argument an `AlloLaw`). Corners: each set of the names that are
no reactant (catalytic modifiers, regulators) at saturation with at most one ligand per
regulatory site (a saturated site holds one ligand, and the rate between two such corners is
monotone, so the maximum is the same), the others at 0; products at 0; substrates grow, then
the corner's names. Each factor of the MWC law keeps its leading face; the leading branches'
faces give the ratio; errors as in kcat.jl. **Integration:** `kcat(m, params)` must use this
method for an allosteric `m` (today it calls `kcat_terms(m)`, the plain one). `kcat_mask`
keeps the flip form when it is a free isomer on a dead-end branch.

### Printing (decision 8, per conformation)

```julia
rate_equation_string(m, A::AlloLaw, sol = solve(m)) -> String   # new methods (2nd arg AlloLaw)
coefficient_form(m, A::AlloLaw) -> String
```

- Each conformation `s` prints over its own normalizer (`printed`): `N_s` (when it prints
  apart), `Q_s` (unless it does not turn over), `D_s`; the flip weights over those normalizers,
  `F_A` and `F_I`, only when they differ and are not 1; `W_A`, `W_I` (regulatory sites, bare
  `K_Rreg` for `:EqualAI`); `v = E_total * (ΠS - ΠP / Keq) * (F_I^n * Q_A * D_A^(n-1) * W_A +
  L * F_A^n * Q_I * D_I^(n-1) * W_I) / (F_I^n * D_A^n * W_A + L * F_A^n * D_I^n * W_I)`
  without factors of 1. The full form starts with the active law's header (parameters,
  concentrations, Wegscheider and Haldane lines) and reads in order; the coefficient form puts
  `v` first, then `Q_A = nA_0 + …`, `D_A = 1 + dA_S * S + …`, then the legend.
- Both forms evaluate to the law on all 2,865 admitted laws (worst 9.8e-154); every printed `D`
  term and flip-weight term is dimensionless; no `^1`, `1 ^`, `* 1`, `L * 1` (`check_print.jl`).

### Notes for integration (allosteric)

1. **Shared printing helpers.** `mwc_lines` repeats coefficient_form's label/unit logic (about 6
   lines) and takes the header from `rate_equation_string(m, A.active, sol)` by splitting its
   lines; print.jl could expose `header(m, L, sol)` and `coefficient_sum(prefix, groups,
   names)` instead. The normalizer's printed-apart factor is passed to `term_text` as an extra
   variable `N_A` (its `over` flag prints a literal `N`).
2. **Printed locals** `N_s`, `Q_s`, `D_s`, `F_s`, `W_s` and coefficient names (`nA_0`, `dI_S`) can
   collide with a metabolite of that name; nothing checks it (as for the plain `N`).
3. **kcat semantics to decide:** PFK (specs/35) errors because ATP is a substrate and an
   `:OnlyI` regulator, so the saturating rate depends on ATP:F6P (the package treats the two
   roles separately). AD3b/AD3d give kcat 0: their substrate-inhibition copy binds E and E(S),
   closing a square with the catalytic forms, so `catalytic_core` keeps it.
4. A fitted constant can be absent from every polynomial (3 in specs/40, specs/42): the speed of
   a steady-state step that carries no flux in its conformation (decision 12 keeps such steps).
5. **Shared kcat helper.** `kcat_terms(m, ::AlloLaw)` repeats the plain `kcat_terms`' leading-face
   test and its three error messages (about 12 lines); one `face_ratio(num, den, deg)` could
   serve both.
6. **Dispatch.** `rate_equation_string(m)` and `coefficient_form(m)` (print.jl defaults) and
   `kcat(m, params)` must route an allosteric `m` to the `AlloLaw` methods.

### Names allosteric-stage files define
`AlloLaw`, `allo_law`, `cross_weights`, `kcat_mask`, `admit`, `conformation`, `mwc_lines`, and
methods `kcat_terms(m, ::AlloLaw)`, `rate_equation_string(m, ::AlloLaw[, sol])`,
`coefficient_form(m, ::AlloLaw)`.

## Codegen (codegen.jl)

This replaces the core's proposal for `codegen.jl` above where it differs.

### The generated rate equation (decisions 1 and 9)

```julia
rate_equation(em::EnzymeMechanism{S}, concs::NamedTuple, params::NamedTuple)  # @generated
rate_equation(m::Mechanism, concs, params)   # compiles the singleton first (allocates)
fitted_params(em), parameters(em), metabolites(em)   # @generated: constant tuples
```

- **Generator.** `body(program(decode(S...)))`: the Sig is the only type parameter; the
  generator decodes it, solves, derives the `Law` (`law`) or `AlloLaw` (`allo_law`) and emits
  the body; no memo. It does not call `admit` or `check_boundary` (boundary checks). One method
  instance per mechanism and argument types; `decode` and the generator's helpers compile once.
- **Inputs.** Read by name with `getfield` (positional loads once compiled): a `NamedTuple` in
  any field order works; fields the law does not read are ignored. Generic in the number type
  (integer literals, `inv`, `sqrt`), so BigFloat inputs give a BigFloat rate.
- **The generated accessors** `fitted_params(em)`, `parameters(em)`, `metabolites(em)` return
  the tuples of `fitted_params(m)` etc. as constants, so a fitter's `NamedTuple{names}` is
  type-stable (the package has the same `@generated` accessors). One instance each per
  mechanism.

### The program (codegen's IR)

```julia
struct Program{T}
    vars::Vars                         # the law's inputs, then the locals #1, #2, …
    defs::Vector{Pair{Int, Poly{T}}}   # local variable => its polynomial (earlier variables)
    num::Poly{T}; den::Poly{T}         # v = E_total·num/den
end
p = program(m)                         # plain or allosteric
p = program(m, sol, L::Law)            # num = H·Q, den = shift·Σ_g σ_g·T_g
p = program(m, sol, A::AlloLaw)        # MWC, below
Q, D = state!(vars, defs, m, sol, L, mask, inactive)  # one conformation's terms
local!(vars, defs, q) -> variable      # the local holding q (deduplicated by polynomial)
factor!(vars, defs, p) -> Poly         # p as one term: concentration content × local
ex = body(p)                           # the Expr the generator returns
```

- **Dependents are locals.** Each distinct power product of independents and `Keq` that a
  dependent equals is one local, computed once (class members with one right side share it); a
  dependent equal to one symbol (a fold) or to 1 (a pin) reads it directly. These locals appear
  only in the segment sums `σ_g` and tree sums `T_g`, which `state!` takes from
  `bare_law(m, mask, T)` (bare constants mapped to independents or dependent locals). `Q` and `H`
  come from the `Law` (in independents: `Q` exists only after the substitution). Substituting
  the dependents inside `σ_g`/`T_g` instead costs 17–48% more time per call (median 32%) on the
  40 largest plain laws (`slice/codegen/logs/variants2.log`).
- **Lowest terms per local.** `factor!` divides each multi-term polynomial by its concentration
  content and keeps that monomial in the term that uses it, so no local divides by a
  concentration and every term of `num`/`den` has nonnegative concentration exponents: the code
  is finite at zero concentrations wherever the law is.
- **MWC.** Locals: each state's `Q`, `σ_g`, `T_g` (`state!`, inactive with
  `inactive_mask(m)`), the cross weights `fA`, `fI`, the state terms `xA = fI·DA`,
  `xI = fA·DI`, each site's two bases. `num = H·(fI·QA·xA^(n−1)·WA + L·fA·QI·xI^(n−1)·WI)`,
  `den = xA^n·WA + L·xI^n·WI` (`AlloLaw`'s formula with `QI = 0`, `DI = 1` when
  `A.inactive === nothing`).
- **Body.** In variable order, each variable read: its load or its local's Horner form, then its
  reciprocal (`x5m1 = inv(x5)`) and roots (`x5p2 = sqrt(x5)`, `x5m2 = inv(sqrt(x5))`) if read.
  `x^e` reads `x` to the integer part of `e` and `x^(±1/q)` to the rest. Greedy Horner (the local
  most terms hold first, ties by name), binary `+`/`*`, a polynomial's negative terms
  subtracted. `return E_total * num / den` is the only division by a sum.
- **Names.** Locals are positional (`x<i>`); variables of locals are named `#k` in `vars`, which no
  metabolite or parameter name can be: no clash is possible.

### Notes for integration (codegen)

1. `state!` recomputes `bare_law` (the engine already computed it inside `law`); a `Law` that
   kept its bare law or its column substitution would save that (µs to 3 ms per generator run).
   `Law.sigma`/`Law.trees` (in independents) are not read by codegen; `allosteric.jl` reads
   `L.trees[1]` (the flip weight).
2. The fitter needs the three generated accessors; with a loss like the package's `loss!`, a
   fresh mechanism compiles about 21 method instances (3 of them codegen's) and 24–38 KB
   (`slice/codegen/logs/budget.log`).
3. Names codegen.jl defines: `Program`, `local!`, `factor!`, `state!`, `program`, `loc`,
   `power`, `literal`, `horner`, `emit`, `body`, `rate_equation`, and methods of
   `fitted_params`, `parameters`, `metabolites` for `EnzymeMechanism`.
4. Checks: `slice/codegen/run_all.sh` (values `check_values.jl`, `check_allo.jl`; benchmarks
   `bench.jl`, the package's `bench_pkg.jl`; `budget.jl`; `summary.jl` prints the tables).

## Integration (key.jl, pipeline.jl, and the reconciling edits)

### The dedup key (decision 10)

```julia
struct Canonical; order::Vector{Symbol}; text::String; key::String; end
c = canonical(L)                      # L a Law or an AlloLaw, before any code generation
comps, units, header = key_parts(L)   # the law's polynomials in role order, their units
```

- **Components and units.** Plain: `[H, Q, D]` with units `[H]`, `[Q, D]`. MWC: `[H, QA, DA, fA,
  QI, DI, fI, (WA_s, WI_s)…]` with units `[H]`, `[QA, DA, fA]`, `[QI, DI, fI]` and one per
  site: each unit is a set of components that rescale together without changing the rate.
  The header names the shape (`plain`, or `mwc|n=…|sites=…`), then `|P=` (the fitted count)
  and `|absent=` (fitted constants no term holds).
- **Normalization.** Each unit is multiplied by `lowest_shift` over every variable (rational
  exponents) and divided by its integer content. A common factor that is no monomial stays.
- **Labels.** Fixed variables: concentrations **by position** in `metabolites(m)` (`c1`, `c2`,
  …), `Keq`, `L`. Every other variable is a fitted constant. Terms are coloured by
  `(component, coefficient, fixed factor)`; colour refinement on the term–constant incidence
  (exponents as edge labels); while colours tie, individualization keeps the smallest text (a
  safety net: P8 never needed it; the toy 6-cycle needs it).
- **Invariances.** Renaming or reordering the fitted constants; storage order, directions and
  group labels (the mechanism is identical); same-model regroupings (the solve is identical);
  renaming metabolites without changing their sorted order. A renaming that reorders the
  metabolites changes positions, and can change the solve's tie-breaks (decision 4), so it can
  change the key.
- **Fits by canonical position.** `c.order` lists the fitted constants in canonical order (any
  absent last, by name): member `b` of a class reads representative `a`'s value of
  `a.order[i]` as `b.order[i]`.

### The pipeline

```julia
d = derive(m)   # (; mechanism, solution, law, normalizers, printed, coefficients, key, code,
                #    compiled, kcat)
```

`derive` runs `check_boundary`, `admit` (an error with its message when rejected), `solve`,
`rate_law`, the normalizer of each conformation, both printed forms, `canonical`, `body(program)`
and `compile`, then `kcat_candidates` (the `kcat: …` error message in place of the candidates
when there is no kcat). Atom conservation stays with the caller (the type has no reaction).

### Edits that reconcile the stages

- `types.jl` now holds `segment_offsets` and `unbalanced_blocks` (moved from `enum.jl`, so
  `kcat.jl` no longer calls a later file), and `bottomless` is written on `segment_offsets`
  (verdicts identical on all 4,917 records).
- `poly.jl`: `lowest_shift(ps, counted = every variable)`; `conc_shift(ps, t)` is its
  concentration case.
- `kcat.jl`: `kcat_error(what)` and `face_ratio(n, d, growing)` serve both kcat methods;
  `kcat(m, params)` moved to `allosteric.jl`.
- `allosteric.jl`: `inactive_terms(A)` (`(QI, DI)`, `0` and `1` for the flip form alone),
  `rate_law(m, sol = solve(m))` (`Law` or `AlloLaw`), `kcat_candidates(m, sol = solve(m))` (plain
  or MWC), and `kcat(m, params)` routes through it.
- `print.jl`: `coefficient_sum(prefix, gs, names)` serves both coefficient forms;
  `rate_equation_string(m, L::Law, sol)` and `coefficient_form(m, L::Law)` take the law
  explicitly; `print_allo.jl` adds the one-argument `rate_equation_string(m)` and
  `coefficient_form(m)`, plain or MWC.
- `codegen.jl`: `program(m)` uses `rate_law`.
- `engine.jl`: the cap is `complexity(m, mask) > MAX_TERMS`, checked in `bare_law` before the
  trees are enumerated: segments times the spanning trees of the segment *multigraph*
  (parallel steady-state steps counted apart), today's `_eq_complexity` guard. The simple-tree
  count it replaces let ter-ter laws with V×τ up to 61,600 through (10,960 D terms, about
  20 s to derive, 42,672 operations); `spanning_trees` no longer checks a cap.
- `enum.jl`: `expand(m, declared, atoms)`, one round of the three moves: children once each,
  `check_boundary` and `check_atoms` asserted, `admit` filtering.
- Every stage output is unchanged on all 4,917 converted records (`integrate/regress.jl`,
  `compare_regress.jl`): solve, names, laws, both printed forms, kcat candidates, bodies,
  boundary and admission verdicts, flux, and the moves on the 4,914 records with at most six
  modifiers. The dead-end move enumerates every subset of the bound modifiers (as the
  package's does), so `p1_adv_pkg`'s 20-inhibitor record exhausts memory there.

### Names integration files define
`Canonical`, `key_parts`, `ranks`, `refine`, `encode`, `canonical`, `derive`; and in other files
`lowest_shift`, `kcat_error`, `face_ratio`, `inactive_terms`, `rate_law`, `kcat_candidates`,
`coefficient_sum`, `expand`.

## Skeptic (no source changes)

The skeptic stage changed no file in `slice/src`. Its scripts in `slice/skeptic/scripts/`
are reusable checks:

- `indep.jl`: an independent route. `my_system(p, tags)` builds the constraint system over
  Pop steps itself (own cycle basis, ties aligned by metabolites, then bound-content change,
  then conformational change); `my_solution`/`my_sample` give consistent constants by an
  exact RREF in random column order; `my_network`/`my_mwc` are steady-state and MWC oracles;
  `my_core` recomputes the dead-end filter; `oracle_corners` gives kcat by saturation.
- `spotlib.jl` `check(r, d, rng)`: one derived law against the route (law, compiled code in
  Float64 and BigFloat, both printed forms, `all_parameters`, kcat).
- `advlib.jl` `build(...)` turns a step spec into a Pop record; `eps_limit` is decision 6's
  limit as one two-conformation network with the `:OnlyA` steps vanishing like `ε`.

Notes for integration (each with a reproducer in `slice/skeptic/scripts`):
1. `canonical` labels concentrations by position, so models that differ only in which
   inhibitor or regulator binds get one key (`keytest.jl`); label them by name.
2. No reserved-name check: a metabolite `E_total` or `N` misprints the law, `Keq` and `L` fail
   inside `var!` (`adv_cases.jl`); the package rejects `Keq` and `E_total` today.
3. An all-RE mechanism derives `Q = 0` and prints `* () /` (the package rejects it); a
   declared substrate no step uses fails in `divexact` (`inv/allre_slice.jl`).
4. Coefficient names collide for metabolites such as `A`, `B` and `A_B` (`adv_cases.jl`).
5. `power` emits `x^(1/q)` with a Float64 exponent: BigFloat results lose precision for
   `q = 3` (`cube_root`).
6. kcat is silently 0 for about 7% of ter-ter laws (`kcat_zero.jl`).

## Critic (no source changes)

The critic stage changed no file in `slice/src`. Its checks are in `slice/critic/`
(`timing_check.jl`, `kept_names.jl`, `density.py`; logs in `critic/logs/`). Entries above
that no longer describe the code:
- Core, "Bare per-step constants": "Slot 2 (`kr`) is ... the solve priority's `against`" and
  "Content names for ties: `step_text`". The solve's `against` is not slot 2 for a product
  binding, and its ties use the packed content `(o.a, o.b, o.cons, o.rel)`, not `step_text`.
- Engine, "The law": `spanning_trees` no longer errors at `MAX_TERMS`; `bare_law` checks
  `complexity(m, mask)` per mask (an allosteric law's two conformations are capped apart,
  not summed).
- Enum note 1: `cycles` no longer calls `isassigned`; the slowdown it describes is gone.
- Integration, "The dedup key": concentrations labelled by position give false merges
  (Skeptic note 1); P8 labelled them by name.
