# Design: enumeration combinatorics doc, shared-catalytic-site constraint, raw-loss LOOCV

Date: 2026-07-21

Three independent changes ship on one branch as three commits, one PR:

1. A documentation subsection that quantifies how fast the enumeration engine's
   candidate count grows, so users see why filters are mandatory.
2. A `shared_catalytic_site` constraint on `EnzymeReaction` that forbids a
   substrate/product pair from co-occupying the catalytic site, pruning
   mechanisms known to be wrong (ATP/ADP, NAD/NADH).
3. Removing the redundant outer `log` from the LOOCV `cv_score`, mean ± SE, and
   permutation test, so model selection compares raw per-fold losses.

The three touch disjoint code, so their commits stay separate and reviewable.

---

## Change 1 — Enumeration combinatorics subsection

### Goal

The enumeration page (`docs/src/identify/enumeration_engine.md`) explains each
move but never states how many candidates the moves produce together. A reader
cannot tell whether the search visits dozens or millions of mechanisms. The new
subsection makes the growth concrete and ties it to the filters documented on
the model-selection page.

### Design

Append a new `##` section after "The seven expansion moves" (currently ends at
line 181), titled **"How many mechanisms? The combinatorics of enumeration."**

The section derives the candidate count as a product of independent
combinatorial factors, each a function of the reaction's counts —
`n_substrates`, `n_products`, `n_competitive_inhibitors`,
`n_allosteric_regulators`:

- **Catalytic binding topologies** — ordered binding sequences plus the
  random-order topology, growing with `n_substrates!` and `n_products!`.
- **Substrate/product dead-end competition patterns** — the number of bipartite
  graphs on `n_substrates × n_products` in which every reactant has degree at
  least one (the source of catalytic-site co-occupancy, and the exact quantity
  Change 2 prunes).
- **Per-move multipliers** — each of the seven moves contributes a factor:
  rapid-equilibrium/steady-state flips (up to `2^{groups}`), kinetic-group
  splits, competitive-inhibitor sites (an inhibitor competition-pattern count),
  allosteric promotion, allosteric ligands and their states, `:EqualAI`
  relaxations, and regulatory-site merges.

The section then plugs in a **random-order bi-bi with two regulators**
(`n_substrates = 2`, `n_products = 2`, one competitive inhibitor, one allosteric
regulator) and multiplies the factors to a concrete total in the tens of
thousands. It closes by naming the filters that tame this total — the beam,
`eq_complexity_filter`, required-regulator seeding, and the new
`shared_catalytic_site` — with a cross-reference to the model-selection page.

Format matches the page: prose plus one small factor→count table, no runnable
block. The page ships zero runnable blocks today, and a live `@example` would
couple the doc to internal signatures and pay a doc-build cost for counts that
by design explode.

### Honesty and verification

The section presents the product as a growth estimate, not the engine's exact
output: structural deduplication and validity pruning trim the raw product, but
the growth rate is the point. Every factor formula is cross-checked against
actual `EnzymeRates.init_mechanisms` counts for uni-uni, uni-bi, and bi-bi
reactions while writing, so the numbers in the prose are trustworthy. The
verification runs offline and ships no code in the doc.

### Files

- `docs/src/identify/enumeration_engine.md` — new `##` subsection at end.

No `make.jl` change: the subsection lives inside an existing page.

---

## Change 2 — `shared_catalytic_site` constraint

### Goal

Let a user declare that a substrate/product pair binds the same catalytic site,
so no enumerated mechanism ever has both bound to the catalytic site at once.
The constraint targets physiological pairs (ATP/ADP, NAD/NADH). It must not
touch a metabolite acting as a competitive inhibitor at a non-physiological
site — that binding is a separate enumeration path and stays enabled.

### API

Declared on the reaction, since which reactants share a site is a biochemical
fact about the enzyme, not a per-search preference.

Constructor keyword:

```julia
EnzymeReaction(reactants, regulators, mults;
               shared_catalytic_site = [(:ATP, :ADP), (:NAD, :NADH)])
```

DSL label, each pair parenthesized and unordered (the constructor infers roles):

```julia
@enzyme_reaction begin
    substrates: ATP[...], NAD[...]
    products:   ADP[...], NADH[...]
    shared_catalytic_site: (ATP, ADP), (NAD, NADH)
end
```

The constructor keyword and the DSL both take `(:sub, :prod)`-shaped tuples, so
the two surfaces read identically.

### Storage and validation

- **Field** (fourth on `EnzymeReaction`, `src/types.jl:329`):
  `shared_catalytic_site::Vector{Tuple{Symbol,Symbol}}`, normalized to
  `(substrate, product)` order and sorted. Canonical storage matches how
  `reactants`, `regulators`, and `allowed_catalytic_multiplicities` are already
  sorted in the constructor; the `==`/`hash` dedup at `src/types.jl:401` depends
  on it.
- **Keyword** on the inner constructor (`src/types.jl:334`), default
  `Tuple{Symbol,Symbol}[]`. Existing three-positional calls — the two reaction
  macros, `IdentifyRateEquationProblem`, and the tests — keep working unchanged.
- **Validation**, added right after `sub_names`/`prod_names` are built
  (`src/types.jl:351-352`): each pair holds exactly two names, one a declared
  substrate and one a declared product. Reject, each with a specific message: a
  name that is neither substrate nor product, a name that is only a regulator,
  two substrates, two products, and duplicate pairs.
- **`==`/`hash`** (`src/types.jl:401-404`) extended to include the new field.
- **Accessor** `shared_catalytic_site(r::EnzymeReaction)` added beside the other
  field accessors (`src/types.jl:392-395`).
- **Struct docstring** updated to document the keyword.

### DSL parsing

Add `:shared_catalytic_site` to `_VALID_REACTION_LABELS` (`src/dsl.jl:41`) and a
branch in `_parse_reaction_block` (`src/dsl.jl:73-95`). Each value arrives from
`_parse_labeled_line` as `Expr(:tuple, name_a, name_b)`; a new helper reads the
two symbols per tuple and emits `(:name_a, :name_b)`. Any value that is not a
two-symbol tuple errors with a message naming the offending entry. The macro
passes the collected pairs to the constructor as the `shared_catalytic_site`
keyword.

### Enforcement

The only place a substrate and a product can end up co-bound to the catalytic
site is `_expand_substrate_product_dead_ends(topos, r)`
(`src/mechanism_enumeration.jl:964`); the main catalytic cycle never co-binds a
substrate and a product. That function gates dead-end forms with
`_competition_patterns(sub_names, prod_names)`
(`src/mechanism_enumeration.jl:869`), which enumerates every bipartite
"cannot-co-occupy" graph.

Enforcement filters that pattern list to patterns whose forbidden-edge set is a
superset of the declared pairs. Because the complete-bipartite pattern contains
every edge and always survives the filter, the retained list is never empty. No
retained pattern permits a dead-end form binding both members of a declared
pair, so no such form is ever built.

The reaction `r` is already in scope at the call site, so the filter reads
`shared_catalytic_site(r)` with no new plumbing. Competitive-inhibitor binding is
added later on a separate path (`_expand_add_dead_end_regulator_native`,
`src/mechanism_enumeration.jl:1607`) using `CompetitiveInhibitor` structs, which
this substrate/product-keyed filter never touches.

### Tests (TDD)

New tests in `test/test_types.jl` (constructor), `test/test_dsl.jl` (DSL
parsing), and `test/test_mechanism_enumeration.jl` (enforcement):

- Constructor rejects: unknown name, regulator-only name, two-substrate pair,
  two-product pair, duplicate pair — each asserting the specific error text.
- Constructor normalizes `(:ADP, :ATP)` to `(:ATP, :ADP)` and sorts pairs.
- DSL `shared_catalytic_site: (ATP, ADP), (NAD, NADH)` produces the same field as
  the constructor keyword; a malformed entry (bare symbol, three-tuple) errors.
- No mechanism from `init_mechanisms` for a bi-bi with a declared pair has both
  members bound to the catalytic site.
- Declaring a pair strictly reduces the `init_mechanisms` count versus the same
  reaction without the constraint.
- A metabolite declared both as a substrate and as a competitive inhibitor
  still yields its inhibitor dead-end form under the constraint.

### Files

- `src/types.jl` — field, keyword, validation, `==`/`hash`, accessor, docstring.
- `src/dsl.jl` — label, parse branch, helper, macro wiring.
- `src/mechanism_enumeration.jl` — pattern filter in
  `_expand_substrate_product_dead_ends`.
- `test/test_types.jl`, `test/test_dsl.jl`,
  `test/test_mechanism_enumeration.jl` — tests.

---

## Change 3 — Raw-loss LOOCV

### Goal

The per-fold loss is already a mean squared log-ratio, so the LOOCV aggregation
takes `log` of an already-logarithmic quantity — a double log. Remove the outer
`log` so `cv_score`, the mean ± SE, and the permutation test all work from the
raw per-fold loss.

### Design

Two source lines carry the outer log:

- `src/identify_rate_equation.jl:1081` — `log.(row.cv_fold_scores)` becomes
  `row.cv_fold_scores`. This one change propagates raw loss into the `n_min`
  selection, the mean paired difference, the SE, and the permutation-test diffs.
- `src/identify_rate_equation.jl:1189` — `mean(log.(v))` becomes `mean(v)`.

Downstream cleanup:

- Rename the diagnostic column `mean_log_loss_diff` to `mean_loss_diff`
  everywhere it appears (result construction, docstrings, tests). After the
  change, "log" in the name is false.
- Remove the `eps(Float64)` floor at `src/identify_rate_equation.jl:1267` and its
  comment. Its sole stated purpose is keeping `log(score)` finite; raw loss of
  zero flows through `mean`, `std`, and the permutation test without trouble.
- Update the `cv_score` definition and the "work in log space" rationale in
  `docs/src/identify/model_selection.md` (lines 82-90, 99-111, 142-145) and the
  docstrings at `src/identify_rate_equation.jl:152-160`, `185-203`, `1036-1070`.

### Behavior change

Removing the outer log shifts the scale of the SE and the permutation test from
log-of-loss to raw loss. Within a single reaction's model comparison the fold
losses sit at a similar scale, so the 1-SE and permutation decisions stay
well-behaved; the change is deliberate and requested.

### Tests

`test/test_identify_rate_equation.jl` recalibration:

- The `_select_best_n_params` fixtures build `cv_fold_scores = exp.([...])` so
  the internal `log` recovers round numbers. With the log removed, the fixtures
  become the round numbers directly — simpler, not more complex. Update every
  `exp.(...)` fixture and its hardcoded diff/SE expectation (lines 671-880).
- Rename `mean_log_loss_diff` to `mean_loss_diff` in the assertions (lines
  299-302, 329-331, and the value assertions).
- The `_onesided_permutation_p` tests (608-669) are scale-agnostic and need no
  change.
- The `_cv_fold_loss` finiteness tests (428-509) survive; drop any assertion
  that leans on the removed eps-floor rationale.

### Files

- `src/identify_rate_equation.jl` — two log removals, rename, eps-floor removal,
  docstrings.
- `docs/src/identify/model_selection.md` — cv_score definition and rationale.
- `test/test_identify_rate_equation.jl` — fixture recalibration, rename.

---

## Out of scope

- No change to the inner `loss!` log-ratio (a different, load-bearing log).
- No generalization of `shared_catalytic_site` beyond substrate/product pairs
  (no substrate/substrate grouping, no site partitions).
- No change to `optional_allosteric_regulators` / `optional_competitive_inhibitors`.
