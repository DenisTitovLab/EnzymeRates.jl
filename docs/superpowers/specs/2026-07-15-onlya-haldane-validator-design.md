# Reject thermodynamically impossible `:OnlyA` mechanisms

## The problem

Two HPC runs (`docs/hpc_results/ldh_hpc_results/2026_07_14_results_2`,
`docs/hpc_results/pfkp_hpc_results/2026_07_14_results`) produced 1354 identical
errors:

```
_kcat_forward: AllostericEnzymeMechanism produced no kcat components —
saturating-substrate pattern not found in numerator
```

The errors are the visible edge of a larger defect. Every allosteric mechanism
carrying an `:OnlyA` binding in either run is thermodynamically impossible:
200 of 200 in LDH iteration 1, 276 of 359 in the PFKP seed tier. Of the LDH
allosteric rows, 6024 of 6401 (94.1%) carry the defect; 1236 crash and 4788 fit
silently and compete in the beam. Their rate equations are wrong, not merely
unscaled.

## Root cause

An `:OnlyA` binding asserts `K_I = ∞`. Thermodynamics propagates that assertion
to the catalytic rate constants. Both conformations share `Keq`, so each cycle
obeys its own Haldane relation:

```
A-cycle:  Keq = (k_A_f/k_A_r) · ∏K_A(prod) / ∏K_A(subs)
I-cycle:  Keq = (k_I_f/k_I_r) · ∏K_I(prod) / ∏K_I(subs)
```

Dividing one by the other, and writing `K_I = K_A/ε` with `ε → 0⁺`:

```
k_I_f/k_I_r  =  (k_A_f/k_A_r) · ∏ε_p / ∏ε_s
```

An `:EqualAI` catalytic tag asserts `k_I_f/k_I_r = k_A_f/k_A_r`, which requires
`∏ε_p/∏ε_s = 1`. When `:OnlyA` appears on exactly one side of the reaction, that
ratio is `0` or `∞`, and the assertion is false. The inactive conformation cannot
run catalysis at the active rate toward a complex it cannot form.

`docs/src/deriving/mwc_allostery.md:207-209` already states the rule:

> its catalytic step should be `:OnlyA` too, since a state that cannot bind the
> substrate cannot catalyze.

Nothing enforces it.

### Why the existing machinery misses it

The constraint solver derives a constraint only from a catalytic thermodynamic
cycle. `rate_eq_derivation.jl:1485` states the principle: "Regulator-site
affinities complete no catalytic thermodynamic cycle, so they are independent."

`_state_allo_mechanism` implements `:OnlyA` as a **graph edit** — it deletes the
binding group from the I state. That edit severs the I-cycle, and a severed cycle
generates no constraint. The tag destroys the evidence of its own inconsistency.

The same edit fragments the I-graph into two rapid-equilibrium segments, so
`d_free_I ≠ 1`. The free-enzyme segment weight then carries the bridging form's
population weight — a metabolite. `_kcat_forward` cross-weights the A-numerator
by that factor, every saturating key acquires the metabolite, and the
products-are-zero filter rejects all of them. The crash and the missing
constraint are one cut seen from two sides.

The `:NonequalAI` case keeps its I-cycle, so the constraint fires and the docs
report the collapse `K_I_S_E = K_A_S_E`. Only `:OnlyA` deletes the evidence.

## The rule

**Validate the Haldane relation directly. Do not count tags.**

For each catalytic thermodynamic cycle, keep every `:OnlyA` binding in the graph
with a formal `K_I = K_A/ε`, build the cycle's Haldane relation, and ask whether
the tagged rate constants can satisfy it as `ε → 0⁺`.

Track the `ε` exponents per ligand symbolically. The monomial `∏ε_i^{a_i}` can be
made `O(1)` exactly when the exponents carry both signs: each `ε_i` is
independent, so writing `ε_i = t^{c_i}` with `c_i > 0` and `t → 0⁺` leaves
`Σ a_i c_i = 0` solvable. This yields:

| `:OnlyA` bindings in the cycle | `∏ε_p/∏ε_s` | `:EqualAI` catalysis |
|---|---|---|
| none on either side | `1` | valid |
| substrates only | `∞` | invalid — forces `k_I_r = 0` |
| products only | `0` | invalid — forces `k_I_f = 0` |
| at least one on each side | tunable to `1` | valid |

The criterion is exclusive-or, not equality. A cycle with two `:OnlyA` substrates
and one `:OnlyA` product is satisfiable: choose `ε_p = ε_s1·ε_s2`. Both sides
carrying `:OnlyA` is the thermodynamically consistent K-system — the affinities
diverge together and their ratio stays free.

**The validator must never substitute a numeric `∞` or compare two infinities.**
Doing so discards the ratio's freedom, which is the whole K-system degree of
freedom, and would reject valid mechanisms that #68's multi-`:OnlyA` move
reaches.

Bindings that complete no catalytic cycle — competitive inhibitors, dead-end
complexes, regulator sites — contribute nothing to the Haldane relation and take
no part in the check. This is the exclusion `rate_eq_derivation.jl:1485` already
applies. Keying on the metabolite name instead of the cycle would misclassify
LDH's `Lactate` and `NAD` competitive inhibitors, which share names with the
products.

## Where the rule lives

### Validator — the `AllostericMechanism` constructor

Throw when a cycle's Haldane relation is unsatisfiable under the tagged rate
constants. The message names the ligand, the tag thermodynamics forces, and the
balanced alternative.

The constructor is the single chokepoint: the enumerator, `compile_mechanism`,
and the DSL all route through it. `dsl.jl:1207-1211` says so directly — "Route
through `AllostericMechanism` so catalytic steps and their allosteric-state tags
canonicalize together."

This departs from the `:NonequalAI` precedent, which collapses the forbidden
degree of freedom instead of rejecting it. The departure is deliberate, and the
precedent is itself slated to change: a silent correction leaves the user holding
a mechanism they did not write and cannot see. The follow-up below brings
`:NonequalAI` onto the same footing.

### Move closure — `_expand_to_allosteric` and `_expand_promote_catalytic_to_onlya`

Both moves write `:OnlyA`, and both must close over the promotions the Haldane
relation forces. Each emits every minimal valid completion:

- promote the chemical step the Haldane forces, or
- promote a balancing binding on the opposite side.

Both completions are enumerated. In the balanced region each move emits both the
`:OnlyA`- and `:EqualAI`-catalysis children. The two tags render the same equation
today, but a `split` can repopulate the inactive forms and distinguish them, so
the search keeps both branches.

**`_expand_to_allosteric` is the primary generator and matters more.** It builds
`base_tags = [:EqualAI …]`, sets one group to `:OnlyA`, and — when that group is
a binding — emits the mechanism with the chemical step still `:EqualAI`:

```julia
for g in 1:n_g
    new_tags = copy(base_tags)
    new_tags[g] = :OnlyA
    if is_iso(rep_step(m, g))
        ...                                    # chemical step: pairs with a regulator
    else
        push!(results, AllostericMechanism(    # binding: catalysis stays :EqualAI
            reaction(m), copy(steps(m)), new_tags, cn, RegulatorySite[]))
    end
end
```

Every binding group yields one invalid mechanism at the root of the search. The
class is not drifted into, it is born: LDH is invalid 200 of 200 at generation 1,
and PFKP's seed tier fails 118 of 359 before any expansion runs.

The moves and the validator must land together. A validator alone would throw
uncaught during seeding.

### The derivation does not change

Once the constructor rejects invalid mechanisms, the graph-edit semantics are
correct for every mechanism that survives, and `d_free_I` returns to `1`. An
earlier candidate — pruning I-state forms that bear an `:OnlyA` ligand — is
unnecessary and is dropped.

## Reachability

`_expand_to_allosteric` and `_expand_promote_catalytic_to_onlya` are the only
moves that write `:OnlyA`. `change_allo_state` writes only `:NonequalAI`;
`add_dead_end` appends `:EqualAI`; `split` copies a tag. Today exactly one of
2170 generation-1 LDH mechanisms reaches a docs-compliant sibling, and that one
mechanism is itself invalid — so rejecting without fixing the moves would make
the valid mechanism unreachable.

Closing the move at the point of promotion produces the valid mechanism directly:

```
all-:EqualAI --promote S--> family A --promote P--> family B
```

Family A (`S::OnlyA`, catalysis `::OnlyA`, `P::EqualAI`) and family B
(`S::OnlyA`, `P::OnlyA`) are distinct hypotheses, not duplicates. Their L-terms
differ — `L·(1 + P/K_P_E)²` against `L·1²` — because family A's inactive
conformation still binds product and family B's binds nothing. The second step is
Δ0 (both mechanisms fit four parameters), so the beam explores it within a tier.

## Scope

**Out of scope: the fourteen chemical-step errors.** Fourteen PFKP errors carry
`:OnlyA` on a chemical step with no `:OnlyA` binding anywhere — the V-system that
`_expand_to_allosteric`'s `is_iso` branch builds. A V-system is thermodynamically
sound, so the validator calls them valid and they still fail. They are a separate
defect. This change addresses 1340 of 1354 errors (99.0%). Whether the
cycle-based validator catches them anyway is worth checking during
implementation; do not assume it will.

**Planned follow-up: reject a collapsed `:NonequalAI`.** The same principle
applies one rung in: a `:NonequalAI` tag whose degree of freedom thermodynamics
annihilates is also a lie, and it is also silently corrected today. Measured on
60 real LDH types, 33 (55%) collapse fully; the rest survive because a
steady-state binding or a chemical step keeps the speed (`kon·koff`) while only
the affinity ratio collapses. About 11% of allosteric rows carry `:NonequalAI`,
so roughly 6% would newly error.

It ships separately because its surface is disjoint. `change_allo_state` writes a
single `:NonequalAI`, and a K-system "needs two coupled `:NonequalAI` bindings,
not one" — so that move's arity changes rather than merely closing over a forced
promotion. And `mwc_allostery.md:157-172` is a live `@example` that exists to
demonstrate the collapse; under the new rule it throws and the doc build fails,
so the teaching example must be rewritten.

**Known and unaddressed.** In the balanced region every catalytic tag yields the
same equation — `:EqualAI`, `:OnlyA`, and `:NonequalAI` agree modulo the
`k_ES_to_EP`/`k_A_ES_to_EP` rename, and the solver already collapses
`:NonequalAI`. `eq_hash` hashes the rendered string and will not merge them, so
the search fits duplicates. We enumerate them anyway: a `split` move can
repopulate the inactive forms and make the tag identifiable again, and that
route is untested. Merging renaming duplicates belongs to the existing `eq_hash`
work.

## Cost

**Ground-truth gates.** All three `:OnlyA` gates in
`test/allosteric_ground_truth.jl` encode the invalid combination. Each carries a
phantom I-state species (`:ES_I`, `:EAB_I`) reached only by inactive catalysis,
and each omits the corresponding flip. They will not construct under the
validator. Rewriting each one deletes the phantom species, its edges, and its
`cat_edges` entry, and retags the chemical step `:OnlyA`. The rewritten gates
match the derivation to 4.5e-7 or better.

The rewrite retires most of #69's own coverage. Gates 1 and 2 exist to exercise
the fragmenting-I normalization; once `d_free_I` returns to `1` they take the raw
branch and no longer reach it. Gate 3 keeps the cross-weight branch covered
(`d_free_A = koff + k·B/K_B ≠ 1`), and validates under the corrected physics.

**Mechanism specs.** Two entries in
`test/mechanism_definitions_for_test_enzyme_derivation.jl` carry the invalid
combination: "LDH i-state NonequalAI 6-group" and "LDH i-state NonequalAI
5-group". Their `expected_n_wegscheider_constraints` and
`expected_n_independent_params` move. Both run without an ODE test or an
analytical rate function. PFK-1, HK, PK, and `m_OnlyA_prod` already comply.

**Goldens.** `test/reference/allosteric_golden_reference.txt` regenerates for
those two specs. The rest stays byte-identical.

**Enumeration tests.** The largest surface —
`test_mechanism_enumeration.jl` mentions `OnlyA` 125 times. `biuni_seed()`
(`:4058-4070`) is itself an instance of the invalid combination and seeds four
promote-move testsets; it will not construct. Both moves now emit more children
per parent, so every child-count assertion for `to_allosteric` and
`promote_catalytic_to_onlya` changes.

**Docstrings.** `_expand_promote_catalytic_to_onlya` claims Δ0 and claims every
promotion is distinguishable. Both claims are false today: the real LDH type
promotes 6 parameters to 6, 7, 7, and 8, and 37 HPC rows show +1. The Δ0 test
passes only because `biuni_seed()` happens to satisfy it. Rewrite the docstring
to state that the parameter count varies, and drop the Δ0 test or rescope it to
the seed it actually describes.

**Coverage.** `_kcat_forward` is tested only through `analytical_kcat_fn`, set on
the four already-compliant specs. The invalid class has no kcat coverage at all.
That is how 1354 codegen errors reached the cluster. The rewritten gates assert a
finite kcat.

**Performance.** Unaffected. The change touches construction and derivation only,
and `_kcat_forward` runs post-optimization from `rescale_parameter_values`. The
constructor gains a bounded cycle walk on the enumeration path.

**Data.** Every prior allosteric LDH and PFKP result is invalid, not only the
1354 errored fits.

## Accepted risks

**`n > 1` has no ground truth, and we ship without one.**
`test/allosteric_ground_truth.jl` solves a single protomer. Every number here and
in #69 comes from `n = 1`; the real LDH mechanisms use `n = 4`. The `^n` structure
has no mass-action reference. The Haldane argument is `n`-independent and the fix
was confirmed structurally on the real `n = 4` type, but no `n = 4` value check
exists and building one is not tractable. Accepted.

## Open questions

1. **Refitted loss and model selection are unmeasured.** Every distortion figure
   comes from random parameter values. Equations change for 200 of 200 LDH types
   and 27 of 71 gain parameters, so the selected mechanism will move. By how much
   is unknown, and only an HPC rerun answers it.
2. **The fourteen chemical-step errors** need their own diagnosis.
3. **Ping-pong fails silently.** Its `d_free_I` carries a substrate, and the
   products-are-zero filter fires only on products. The mechanism derives, cross-
   weights, and returns a wrong number.

## Evidence

The Haldane algebra above carries the argument. A mass-action ODE sweep confirms
it numerically and independently of the derivation code: with `K_I_F16BP` taken
to infinity along a Haldane-consistent path (`k_I_f = k_f·K_F16BP/K_I_F16BP`),
the steady-state flux at zero product converges to the promoted mechanism's
answer.

```
K_I_F16BP     k_I_f     flux at ADP=0     F16BP-bearing I island
1             10        8.910             5.3e-3
1e2           0.1       4.500             5.0e-6
1e6           1e-5      4.45506           4.5e-10
1e9           1e-8      4.455053          4.5e-13

mechanism as tagged (:OnlyA binding, :EqualAI catalysis)  =  0.0
promoted mechanism (:OnlyA catalysis)                     =  4.45505
```

The mechanism as tagged gives flux exactly proportional to ADP — 0 at zero
product, with the entire enzyme pool in the F16BP-bearing inactive island. That
is a faithful simulation of what the tags currently mean, not a second proof:
it reproduces the defect rather than the physics.

A third line of argument — that a zero-population inactive complex cannot absorb
the flux of reverse inactive catalysis — is **unsound and was dropped**. It
assumes the complex is unpopulated, which is a consequence of the Haldane limit
rather than an independent premise. Do not reintroduce it.
