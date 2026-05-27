# Structural Parameter Names — Design

Date: 2026-05-27
Branch base: `refactor-to-concrete-types-instead-of-symbols`

## Goal

Replace index-based parameter symbols (`:K1`, `:k6f`, `:K_G6P_reg1`) with
**structural, direction-independent** names derived from the species and
metabolites a step touches (`:K_ATP_E`, `:kon_ATP_E`, `:k_ES_to_EP`,
`:K_G6Preg`, `:K_G6Pinh`). Names become a pure function of the mechanism's
chemistry, independent of step write-direction and of any positional index.

This lets us:

- Delete `Step.source_idx` and `_rep_idx_for_step` — the only reason they
  exist is to render positional names.
- Delete the type/index-context chokepoint companion `name(::Type{P}, idx)` —
  the `@generated` derivation can render from the (reconstructed) Step value.
- Make cross-mechanism parameter identity visible at the **name** level
  (`:K_ATP_E` is the same symbol in every mechanism that has that step), so the
  rate-equation canonical-hash machinery that bridges rendered-name ↔
  structural token can collapse toward "rendered name *is* the canonical key."
- Remove the "SS steps are not canonicalized" special case in the `Step`
  constructor and dedup (see Decision 2).

## Locked decisions

These were resolved in brainstorming with Denis:

1. **Names-only on the elimination side.** The Haldane/Wegscheider pivot-priority
   scheme (`_dependent_param_exprs_kernel`, prefers keeping free-enzyme binding
   and forward-over-reverse) is **unchanged**. Which parameter becomes dependent
   does not change. Only its *name* changes.
2. **Canonicalize SS steps too.** The "SS steps preserve source direction"
   invariant (CLAUDE.md "Canonical Step Form") existed because the old `:kNf`/
   `:kNr` labels were tied to write-direction. Structural species-pair names
   remove that reason, so SS step storage direction is canonicalized like RE
   steps, and the special case is deleted from the `Step` constructor and dedup.
3. **Hybrid test migration.** Golden output assertions (`expected_factored_num`/
   `expected_factored_denom`, `rate_equation_string`, `parameters()` literal
   tuples) are **regenerated** to structural names — they are outputs we want to
   assert. Hand-derived numerical physics oracles (`analytical_rate_fn`,
   `rate_uni_uni`, `analytical_kcat_fn`) are kept **byte-for-byte unchanged**;
   a permanent positional-remap helper maps `parameters(m)`-ordered structural
   names back to the positional `k1f, k1r, …` the oracles destructure. The
   numerical comparison is the safety net (see Risk R1).

## Naming scheme

`name(s::Species)` already renders structural form names (`:E`, `:E_ATP`,
`:Estar`, `:E_c`). Parameter names compose those with the bound metabolite.

| Parameter | Step kind | Rendered name | Example |
|---|---|---|---|
| `Kd`   | RE binding | `K_<met>_<fromform>` | `:K_ATP_E`, `:K_Glc_E_ATP` |
| `Kiso` | RE isomerization (no metabolite) | `Kiso_<from>_to_<to>` | `:Kiso_EA_to_EstarA` |
| `Kon`  | SS binding, forward | `kon_<met>_<fromform>` | `:kon_ATP_E` |
| `Koff` | SS binding, reverse | `koff_<met>_<fromform>` | `:koff_ATP_E` |
| `Kfor` | SS isomerization, forward | `k_<from>_to_<to>` | `:k_ES_to_EP` |
| `Krev` | SS isomerization, reverse | `k_<to>_to_<from>` | `:k_EP_to_ES` |
| `Kreg` | allosteric regulator | `K_<lig>reg` | `:K_G6Preg` |
| `Kd` (CompetitiveInhibitor role) | dead-end RE binding | `K_<met>inh[_<fromform>]` | `:K_G6Pinh` |
| `Keq`  | — | `:Keq` | |
| `Etot` | — | `:E_total` | |
| `Lallo`| — | `:L` | |

Conventions:

- **`<fromform>`** is the pre-binding enzyme form (the `from_species` of the
  canonicalized RE/SS binding step). It is **always included** so a metabolite
  that binds more than one form (random mechanisms) stays unambiguous. For the
  free enzyme this is just `E`, giving `:K_ATP_E`.
- **`<from>`/`<to>`** for isomerizations are the canonicalized step's species
  pair. Each rate constant is named by its **actual directed transition**, so
  the name is unambiguous regardless of storage order. After SS canonicalization
  (Decision 2) the stored direction is deterministic (lex on species name).
- **T-state suffix** `_T` is appended unchanged: `:K_ATP_E_T`, `:k_ES_to_EP_T`,
  `:K_G6Preg_T`.
- **Uniqueness within a mechanism** is guaranteed because a step is identified
  by `(from_species, to_species, bound_metabolite, is_equilibrium)`, and the
  name encodes enough of that tuple to be injective. See Risk R2 for the one
  edge case (same ligand at two distinct regulatory sites).

## Architecture changes

### `src/types.jl`

- **`Step`**: drop the `source_idx` field. Drop the SS-direction special case
  in the constructor so SS steps canonicalize like RE steps. `==`/`hash` are
  unaffected (already ignore `source_idx`).
- **`Mechanism` / `AllostericMechanism` constructors**: stop assigning
  `source_idx`. Remove the density-of-`source_idx` validation.
- **Chokepoint** (`_param_symbol` family): rewrite the `_param_symbol` bodies to
  emit structural names from a Step value (and ligand/site for `Kreg`). Delete
  `name(::Type{P}, idx)` and `name(::Type{P}, idx, state)` (the index-context
  companion) and delete `_rep_idx_for_step`. The single remaining entry point is
  value-context `name(p::Parameter, m)`.

### Kinetic-group representative

The per-group representative survives (two physically-distinct steps in one
group do **not** share a structural name). Rep selection moves from "lowest
`source_idx`" to "first step in the **canonicalized** group storage order"
(dedup already sorts groups by `_step_canonical_key`, so this is deterministic
and structural). The group's parameter name = the rep step's structural name;
the Pass-1 rename map maps the other members' names to it, exactly as today —
only keyed on structural names instead of `K9 => K4`.

### `src/rate_eq_derivation.jl` and `src/thermodynamic_constr_for_rate_eq_derivation.jl`

- `@generated` callers that currently build symbols via `name(Kd, idx)` switch
  to constructing the rep step's `Parameter` and calling value-context
  `name(p, m)`. The mechanism is reconstructible from the `Sig` at compile time
  (`_mechanism_from_sig`), so the Step values are available.
- Display strings that interpolate indices (e.g. `"K$idx = K$rep"` in the
  Wegscheider/group annotation block) switch to rendering both sides through
  `name(p, m)`.

### Dedup / canonical hash (`src/mechanism_enumeration.jl`, `src/identify_rate_equation.jl`)

- Mechanism struct `==`/`hash` and `_canonicalize_mechanism!` are already
  position-independent and need no change for *mechanism* dedup.
- The rate-equation canonical-hash machinery (`_parameter_canonical_key`,
  the rendered-name → token substitution) can be simplified: with globally
  structural names, the rendered name is already a canonical, position-independent
  key. The cleanup pass collapses the token-substitution layer where it now
  only re-derives what the rendered name already encodes. Exact extent
  determined during the cleanup phase, guarded by existing hash-partition tests.

### Test support

- New permanent helper (test-only): `positional_params(m, structural_nt)` →
  NamedTuple re-keyed to `k1f, k1r, k2f, …` (and `K1, K2, …` for RE) in
  `parameters(m)` order, so unchanged numerical oracles keep working.
- Regenerate all golden-string and `parameters()`-tuple assertions.
- `test/test_chokepoint.jl`: the AST walker still forbids raw `K…`/`k…` symbol
  literals outside the chokepoint; update any allowances tied to the deleted
  index-context entry point.

## Risks and interactions

- **R1 — SS canonicalization × oracle shim.** Canonicalizing SS direction can
  change which physical direction is "forward" relative to what a hand-written
  oracle assumed, and can change group ordering. The positional-remap shim is
  the single adjustment point; if a mapping is wrong the **numerical** oracle
  test fails loudly (it can never silently pass with wrong physics). Resolution
  policy: when a numerical test goes red, fix the *shim mapping*, never the
  oracle. Oracles stay the trusted source of truth.
- **R2 — Name collision: same ligand at two distinct regulatory sites.** `Kreg`
  drops the site index (`:K_G6Preg`). If a mechanism has the same ligand at two
  structurally-distinct sites, the bare name collides. This is rare/possibly
  nonexistent in the enumerated space. Policy: detect at name time and error
  with a clear message (rather than silently merge); add a site discriminator
  only if a real mechanism needs it (YAGNI until proven).
- **R3 — `rate_equation` performance.** Names are compile-time constants; the
  refactor must not introduce allocations or per-call slowdown. The existing
  `test_rate_equation_performance` (`allocs == 0`, `< 100 ns`) is the gate and
  is non-negotiable.

## Sequencing (high level; detailed plan via writing-plans)

1. Rewrite the chokepoint to structural names (keep `source_idx` temporarily so
   the build stays green); regenerate golden strings + add the oracle shim.
2. Move kinetic-group rep to structural-first; switch `@generated` callers to
   value-context `name(p, m)`; delete `name(::Type{P}, idx)` and
   `_rep_idx_for_step`.
3. Remove `Step.source_idx` and its constructor/validation references.
4. Canonicalize SS steps; delete the "SS not canonicalized" special case;
   re-green oracles via shim mapping.
5. Cleanup sweep: simplify canonical-hash token layer, dedup, derivation, and
   enumeration code paths that referenced indices; remove now-dead comments.

## Success criteria

- Full test suite green (incl. Aqua, JET, compile-budget, chokepoint, and the
  `rate_equation` perf gate).
- `Step.source_idx`, `_rep_idx_for_step`, and `name(::Type{P}, idx)` are gone.
- Rendered parameter names are structural for every mechanism in
  `MECHANISM_TEST_SPECS`.
- Net source LOC does not increase (this is a simplification; "less code ==
  simpler code").
