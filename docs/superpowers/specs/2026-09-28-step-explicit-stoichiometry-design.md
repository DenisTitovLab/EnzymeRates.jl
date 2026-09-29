# Steps with explicit stoichiometry (sub-project A)

Date: 2026-09-28. Status: design agreed with Denis section by section on 2026-09-28; Denis
asked for the implementation plan and its execution to follow without a further review stop.
Part of the A + B + C plan in `2026-09-28-identifiable-enumeration-findings-and-plan.md`, whose
point 11 holds the evidence cited here.

## Why

A `Step` today stores at most one free metabolite (`bound_metabolite`) and not the side it is
on. `_step_sides` guesses the side from the shapes of the two species in five branches, and the
numerator's cut search trusts the direction a step was written in. On 43 fused and merged
constructs this gives 3 silent −v results, a silently wrong law when a cut mixes a reversed and
a forward step, 4 BoundsErrors, 4 thermodynamic errors and 1 missing cut. A Theorell–Chance step
(EA + B ⇌ EQ + P) cannot be written at all: the DSL silently drops the second metabolite. A
fused binding and a pure binding of the same metabolite from the same form get the same
parameter names, which silently ties their constants. Sub-project C needs Theorell–Chance and
merged central complexes in the enumerator, so these must work first.

## Goals

1. Every step states its own stoichiometry; the derivation never infers a side and gives the
   same rate law however a step is written.
2. Fused steps, Theorell–Chance steps and steps with several metabolites on one side derive
   correctly and can be written in the DSL.
3. Every mechanism the enumerator emits today keeps its fitted parameter names, its rate
   equations and its place in the enumeration.

## Non-goals

- The enumerator emitting or expanding fused steps (C), the zero-flux and duplicate-inhibitor
  filters (B), MERGE and ELIM moves, and rank on the beam's axis.
- `_kcat_forward` emitting an empty `max()` for a mechanism that never saturates, and
  residual-form parameter names left unquoted in the printed `v =` line. Both are recorded in
  the findings document for separate fixes.

## Design

### 1. Data model

```julia
struct Step
    from_species::Species
    to_species::Species
    consumed::Vector{Metabolite}   # taken up from solution going from → to
    released::Vector{Metabolite}   # given off to solution going from → to
    is_equilibrium::Bool
end
```

- This is the only constructor. The four-argument `Step(from, to, metabolite_or_nothing,
  is_equilibrium)` is removed, with no convenience form (Denis, 2026-09-28).
- Both lists are sorted like `Species.bound`: by metabolite name, a competitive-inhibitor copy
  after the metabolite of the same name.
- The constructor rejects a step whose two forms are equal and a step whose two lists share a
  metabolite.
- Accessors: `from_species`, `to_species`, `consumed`, `released`, `is_equilibrium`, and three
  kinds derived from the fields only:
  - `ligand(s)` returns M when the step is a **pure binding**: it consumes exactly M, releases
    nothing, and `bound(to)` is `bound(from)` plus M with the same residual. The conformation
    may change, as today. Otherwise it returns `nothing`. `is_binding(s) = ligand(s) !== nothing`.
  - `is_iso(s)` is true when both lists are empty. It keeps today's meaning (chemistry or
    conformational isomerization), so every move behaves exactly as today; B and C switch
    chemistry detection to a correct predicate.
  - A **transformation** is every other step: fused chemistry plus release, fused binding plus
    chemistry, Theorell–Chance, several metabolites on one side.
- `bound_metabolite` is removed; its callers use `ligand` or the lists.

### 2. Canonical orientation

- The `Step` constructor stores a pure binding with its metabolite consumed: a release written
  EA → E + A is stored as E + A → EA, as today.
- The `Mechanism` and `AllostericMechanism` constructors orient every other step with the three
  tiers of today's `_canonical_iso_direction`. Tier 1 scores each side by (number of substrates,
  −number of products) over the side's bound metabolites plus its free metabolites, so
  EA + B | EQ + P scores (2, 0) against (0, −2) and `from` becomes the substrate side. Tiers 2
  and 3 are unchanged. Isomerizations have no free metabolites, so their orientation does not
  change.
- Reversing a step swaps its two forms and its two lists.
- The function is renamed to say what it now does (it orients every non-binding step).

### 3. Identity and order

- `==` and `hash` include both lists.
- The sort key becomes (from name, to name, consumed names, released names, flag). For every
  existing step the third field equals today's (the bound metabolite's name, or "") and the
  fourth is empty, so step order, group order, Haldane pivots and fitted names do not change.
  CLAUDE.md marks this order as load-bearing.
- The `_to_sig` / `_step_from_sig` encoding of a step becomes (from, to, consumed tuple,
  released tuple, flag). Consequence: the `mechanism_type` strings in result CSVs written before
  this change no longer decode. There is no compatibility layer.

### 4. Derivation

The derivation reads each step's lists and is built so that reversing a step changes nothing but
which constant is called forward.

- **RE weights.** One rule replaces the two branches of `_compute_alpha`:
  - pure binding: [to] = [from]·[M]/K, with K the dissociation constant (`Kd`, named `K_M_F`);
  - every other RE step: [to]·Π[released] = K·[from]·Π[consumed] (`Kiso`).
  This fixes the BoundsError on RE steps that release a metabolite.
- **King–Altman edges.** Forward rate kf·Π[consumed]·w(from), reverse kr·Π[released]·w(to).
- **Thermodynamic constraints.** Each step's stoichiometry column comes from its lists. The sign
  of an RE constant in the log rows is −1 for a pure binding's `Kd` and +1 for every `Kiso`,
  matching the weight rule. The existing check that every cycle's net change is a whole number
  of turnovers, with no net change of inhibitors or regulators, is unchanged.
- **Rate.** The cut search in `_compute_numerator` is replaced by

  v · den = Σ over SS steps e of ω_e · (kf_e·Π[consumed_e]·w(from_e)·D[seg(from_e)] −
  kr_e·Π[released_e]·w(to_e)·D[seg(to_e)]),

  where ω_e = (copies of the first substrate consumed by e − copies released by e) + u(from_e) −
  u(to_e), and u(f) is the first substrate's exponent in f's RE weight. Flux conservation at
  every form makes this equal to the net consumption of the first substrate for every parameter
  value, not only thermodynamically consistent ones; any complete cut gives the same function.
  So the numerator polynomial is identical to today's wherever today's cut is correct, and a
  dead inactive allosteric state yields the zero polynomial exactly. Only steps with ω_e ≠ 0
  contribute, so the number of terms stays close to today's cuts.
- **Infinite rate.** An RE cycle with nonzero net stoichiometry (an all-RE catalytic cycle) is
  detected while the RE weights are built and raises an error saying the mechanism has no finite
  rate.
- **Deleted**: `_step_sides`, the cut search with its preference order, the "ambiguous
  central-complex cut" error, and the `allow_dead` keyword. Every former caller of `_step_sides`
  reads the lists: `_compute_alpha`, `_raw_symbolic_rate_polys`, `_segment_graph_terms`,
  `_re_segment_extras`, `_bottomless_re_segment`, `_hyperbolic_catalysis`, `_context_form`,
  `_thermodynamic_constraints`, `reactions` and `show`.
- **Unchanged**: the Cha and King–Altman structure, the polynomial machinery, code generation,
  the allosteric per-state derivation (it inherits the lists) and the zero-allocation, sub-120-ns
  `rate_equation` contract.

### 5. Parameters and names

- Pure bindings keep `Kd` / `Kon` / `Koff`, rendered as today: `K_M_F`, `kon_M_F`, `koff_M_F`,
  with the `inh` marker for a competitive-inhibitor copy.
- Every other step uses `Kiso` / `Kfor` / `Krev`, named by form pair: `Kiso_F1_to_F2`,
  `k_F1_to_F2`, `k_F2_to_F1`, with F1 the canonical `from` (Denis, 2026-09-28). A transformation's
  `Kiso` can carry concentration units. `rescale_parameter_values` classifies constants by
  parameter type and scales every SS rate constant alike, so it is unaffected.
- The `Mechanism` and `AllostericMechanism` constructors reject two kinetic groups whose
  representatives render the same name, with an error naming both steps. Today such a collision
  silently ties two constants.
- The structural heuristics that pick group representatives, Haldane pivots and isomerization
  directions (`_free_enz_set`, `_step_priority`, `_group_rep`, `_entry_kind`) read the lists: a
  step with free metabolites plays the part today's metabolite steps play, and a step that
  consumes a metabolite marks its `to` form as bound. For every pure binding and isomerization
  this reproduces the current behaviour, so group representatives and Haldane pivots do not
  move for any mechanism the enumerator emits.

### 6. DSL

- `_split_side` returns a side's enzyme form and all its metabolite terms; the left-hand terms
  become `consumed` and the right-hand terms `released`. Nothing is inferred or dropped.
- `E(A) + B <--> E(Q) + P` and `E + A + B <--> E(A, B)` parse. A side with no enzyme form, or
  with two, is still an error.

### 7. Display

`show` and `reactions` print each stored step with all its metabolites, for example
`EA + B <--> EQ + P`.

### 8. Enumeration code

Mechanical changes only; the enumerator's output must not change.

- Every `Step(` call switches to the list form: `_catalytic_topologies`, `_release_products!`,
  the init dead ends, the dead-end move and its mirror steps, `_flip_group_to_ss`, splits, and
  the allosteric paths.
- Readers of `bound_metabolite` switch to `ligand` or the lists: `_forms_with_binding_step_native`,
  `_bound_at_forms`, `_drop_unbound_regulators`, `_assert_mechanism_invariants`, the regulator
  checks.
- `_assert_step_atom_conserving` checks atoms(from) + atoms(consumed) = atoms(to) +
  atoms(released).
- `_apply_equivalence_grouping` keys on (consumed, released, RE/SS) instead of the bound
  metabolite.
- `_assert_chemistry_is_iso` keeps its role: a parent passed to `expand_mechanisms` may contain
  only pure bindings and isomerizations. C removes it when the moves learn fused steps.
- `_flux_carrying_groups` and the allosteric moves keep detecting chemistry through `is_iso`
  (B and C replace this).

## Test changes

- Tests that call `bound_metabolite` switch to `ligand` or the lists; test helpers that build
  steps switch to the five-argument constructor.
- "Numerator: ambiguous central cut (regulator sibling) raises"
  (`test/test_rate_eq_derivation.jl`) expects the "ambiguous central-complex cut" error. Its
  mechanism has an all-RE catalytic cycle (E(R) + S ⇌ E(S, R) ⇌ E(P, R) ⇌ E(R) + P), so its rate
  is infinite; the test keeps its mechanism and expects the infinite-rate error instead.
- Fixtures with fused steps (for example Segel ordered bi-bi, `E(A, B) <--> E(Q) + P`) change
  parameter names, e.g. `kon_P_EAB` → `k_EAB_to_EQ`. Their textbook-formula tests keep passing
  with the renamed parameters.

## Verification

Every change is made test-first.

1. **Nothing changes for today's mechanisms.** Before any code change, snapshot
   `fitted_params` and `rate_equation_string` (Reduced and Full) for every mechanism in the
   export (`init_mechanisms` plus two levels of the three catalytic moves on reactions R1–R6 of
   the findings document, including R6's level-2 sample) and for `MECHANISM_TEST_SPECS`, keyed by
   the steps written with explicit lists. After the change, regenerate the export with the new
   code and compare: the set of mechanisms, and every name and string, must be identical, except
   the listed fixtures with fused steps. Any other difference is checked against the mass-action
   oracle below; a difference where the new law matches the oracle and the old one does not is a
   fixed bug and is recorded in the PR. Record derivation time per mechanism in both runs and
   investigate if the median grows by more than 20%. This check is a script whose result is
   reported in the PR, not a committed golden file.
2. **Orientation property test.** For every mechanism in `MECHANISM_TEST_SPECS` and in the fused
   suite below, reversing any subset of the written steps gives an `==` mechanism and an
   identical `rate_equation_string`.
3. **Fused and Theorell–Chance correctness against an independent oracle.** A test helper
   (`_testhelper_mass_action_rate`) computes v by brute-force mass action over all forms in
   BigFloat, with no Cha segments and RE steps run at rate scale 1e60. Its parameter points
   come from free energies assigned to forms and metabolites, so every cycle is consistent by
   construction and independently of the package's constraint solver; the package's
   `rate_equation` receives the fitted values and computes its dependents itself. The suite
   covers: uni-uni and ordered fused steps written in both orientations and under both flags;
   merged ordered mechanisms with mixed RE/SS; ping-pong variants; merged random bi-bi, including
   the former missing-cut cases; Theorell–Chance with SS and RE steps; a two-metabolite binding;
   and dead ends on merged forms. Every case must match.
4. **Unit tests**: the `Step` constructor (binding orientation, reversal swapping the lists,
   rejection of a shared metabolite and of equal end forms); `ligand`, `is_binding`, `is_iso`;
   the DSL (Theorell–Chance, two metabolites on one side, the remaining errors); form-pair names
   for transformations; the name-collision error.
5. **Unchanged guards**: the `rate_equation` performance test (zero allocations, under 120 ns,
   every spec) and the parameter-naming chokepoint test pass without edits.
6. The full suite runs before every commit.

## Risks

- **Numerator size.** The telescoped sum could produce more intermediate terms than today's
  minimal cut on some mechanisms. Only steps with ω_e ≠ 0 enter, and the regression run
  measures derivation time.
- **A wrong numerator today.** If the regression finds a mechanism whose old law disagrees with
  the oracle, the change fixes it; the PR lists every such mechanism.
- **Saved CSVs.** `mechanism_type` strings from earlier runs stop decoding (section 3).
