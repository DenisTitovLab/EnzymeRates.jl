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
2. **Canonicalize all iso steps (RE and SS) in the physical-forward direction.**
   Iso steps are stored with `from` = the side further along the catalytic
   cycle's substrate→product progression. The direction question is identical
   for RE iso (`K`) and SS iso (`kf`/`kr`); only the parameter count differs,
   not the storage direction, so one canonicalization rule covers both.
   Iso canonicalization moves out of the `Step` constructor (no reaction
   context) into the `Mechanism` / `AllostericMechanism` constructor (which
   carries the `EnzymeReaction`'s substrate/product sets). RE *binding*
   canonicalization stays in the `Step` constructor (metabolite on `from`
   side — needs no reaction context). The derivation's existing assumption
   "LHS = kf-side" (`rate_eq_derivation.jl:356-360`) becomes physical-forward-
   aligned and *more* correct than today's source-order behavior.

   Algorithm (Mechanism constructor, per iso step):
   - **Tier 1:** `(n_subs_bound, -n_prods_bound)` — higher = more "from".
     Captures atom-balance progression. Decides the 95% case.
   - **Tier 2 (1-hop graph context):** For each species classify by which
     RE binding steps touch it as the free side: `:substrate_only` /
     `:product_only` / `:both` / `:neither`. If from=`:product_only` and
     to=`:substrate_only`, "product-exit → substrate-entry" is forward
     (the cycle-closing/regeneration step). Handles the Segel Iso Uni Uni
     `F ⇌ E` case fully source-direction-independently.
   - **Tier 3 (lex on conformation name):** deterministic fallback for the
     truly-symmetric tail (`:both` ↔ `:both` random-mechanism rearrangement,
     or `:neither` ↔ `:neither` theoretical conformational change with no
     metabolite traffic on either side).

   Residual atoms are NOT a separate tier: by atom conservation, an iso step
   cannot change residual content without also changing bound metabolites,
   so a Tier-1 tie implies a residual tie automatically — residual would
   never be the discriminator.

3. **Active/Inactive (A/I) everywhere.** Rename the allosteric-state notion
   from R/T to A (active) / I (inactive) throughout: taxonomy symbols
   (`:OnlyR→:OnlyA`, `:OnlyT→:OnlyI`, `:EqualRT→:EqualAI`, `:NonequalRT→:NonequalAI`),
   the user-facing DSL annotations, CLAUDE.md docs, and rendered names. The
   state token sits **right after the type prefix**, not at the end, and
   `EqualAI` carries **no token**: `K_ATP_E` (equal), `K_A_ATP_E` (OnlyA),
   `K_I_ATP_E` (OnlyI), both for NonequalAI; `kon_A_ATP_E`, `k_A_ES_to_EP`.

4. **Single priority function for rep + pivot.** Extract a pure
   `_step_priority(step, m)::Int`. The kinetic-group **name representative** is
   `argmin` within the group (least-eliminable = most primary, lexical
   tiebreak); the Haldane elimination **pivot** keeps its existing `argmax`
   behavior among reps. Kinetic groups are kind-homogeneous (all binding or all
   iso), so within a group this picks the free-enzyme-binding member — the
   intuitive name. Naming and elimination thus share one notion of primacy.
5. **Hybrid test migration.** Golden output assertions (`expected_factored_num`/
   `expected_factored_denom`, `rate_equation_string`, `parameters()` literal
   tuples) are **regenerated** to structural names — they are outputs we want to
   assert. Hand-derived numerical physics oracles (`analytical_rate_fn`,
   `rate_uni_uni`, `analytical_kcat_fn`) are kept **byte-for-byte unchanged**;
   a permanent positional-remap helper maps `parameters(m)`-ordered structural
   names back to the positional `k1f, k1r, …` the oracles destructure. The
   numerical comparison is the safety net (see Risk R1).

## Naming scheme

`name(s::Species)` renders structural form names. It is **changed** to drop the
underscore between conformation and bound metabolites (full concat), so form
names read `:E`, `:ES`, `:EP`, `:EATP`, `:EATPGlc`, `:EstarA`. Underscore then
serves purely as the parameter-field separator, removing the
`K_ATP_E_Glc`-is-it-`ATP`/`E_Glc`-or-`ATP_E`/`Glc` ambiguity. Residual markers
(`res`, `+X`, `-Y`) keep their existing form (e.g. `:Estarres+P`).

Parameter names compose the (concatenated) form name with the bound metabolite
and the optional A/I state token. State token placement: right after the type
prefix; `EqualAI` emits no token.

| Parameter | Step kind | Rendered name | Example |
|---|---|---|---|
| `Kd`   | RE binding | `K[_<state>]_<met>_<fromform>` | `:K_ATP_E`, `:K_A_Glc_EATP` |
| `Kiso` | RE isomerization (no metabolite) | `Kiso[_<state>]_<from>_to_<to>` | `:Kiso_EA_to_EstarA` |
| `Kon`  | SS binding, forward | `kon[_<state>]_<met>_<fromform>` | `:kon_ATP_E`, `:kon_A_ATP_E` |
| `Koff` | SS binding, reverse | `koff[_<state>]_<met>_<fromform>` | `:koff_ATP_E` |
| `Kfor` | SS isomerization, forward | `k[_<state>]_<from>_to_<to>` | `:k_ES_to_EP`, `:k_A_ES_to_EP` |
| `Krev` | SS isomerization, reverse | `k[_<state>]_<to>_to_<from>` | `:k_EP_to_ES` |
| `Kreg` | allosteric regulator | `K[_<state>]_<lig>reg` | `:K_G6Preg`, `:K_A_G6Preg` |
| `Kd` (CompetitiveInhibitor role) | dead-end RE binding | `K_<met>inh[_<fromform>]` | `:K_G6Pinh` |
| `Keq`  | — | `:Keq` | |
| `Etot` | — | `:E_total` (unchanged; baked into ~15 sites including the generated rate body and fitter — pure churn to rename) | |
| `Lallo`| — | `:L` | |

Conventions:

- **`<state>`** ∈ {`A`, `I`} appears only for `OnlyA` (`A`), `OnlyI` (`I`), and
  `NonequalAI` (both `A` and `I` variants emitted). `EqualAI` and non-allosteric
  parameters emit no state token.
- **`<fromform>`** is the pre-binding enzyme form (the `from_species` of the
  canonicalized RE/SS binding step), rendered via the concatenated
  `name(Species)`. **Always included** so a metabolite that binds more than one
  form stays unambiguous. For the free enzyme it is `E`, giving `:K_ATP_E`.
- **`<from>`/`<to>`** for isomerizations are the canonicalized step's species
  pair. Each rate constant is named by its **actual directed transition**, so
  the name is unambiguous regardless of storage order. After SS canonicalization
  (Decision 2) the stored direction is deterministic (lex on species name).
- **Uniqueness within a mechanism** is guaranteed because a step is identified
  by `(from_species, to_species, bound_metabolite, is_equilibrium)`, and the
  name encodes enough of that tuple to be injective. See Risk R2 for the one
  edge case (same ligand at two distinct regulatory sites).

## Architecture changes

### `src/types.jl`

- **`name(s::Species)`**: drop the `_` separator between conformation and bound
  metabolites (full concat). Residual markers (`res`/`+`/`-`) unchanged. This
  affects `enzyme_forms` rendering and any golden string showing form names.
- **`Step`**: drop the `source_idx` field. Drop the SS-direction special case
  in the constructor so SS steps canonicalize like RE steps. `==`/`hash` are
  unaffected (already ignore `source_idx`).
- **`Mechanism` / `AllostericMechanism` constructors**: stop assigning
  `source_idx`. Remove the density-of-`source_idx` validation.
- **Allosteric-state taxonomy**: rename `:OnlyR/:OnlyT/:EqualRT/:NonequalRT` →
  `:OnlyA/:OnlyI/:EqualAI/:NonequalAI` in the `RegulatorySite` validator, the
  `AllostericMechanism` validation, the symbolic R/T-rename machinery
  (`_onlyR_syms`, `_T_rename`, `_rename_symbols`, `_zero_symbols_in_poly`
  callers), and the DSL annotation parser. CLAUDE.md updated to match.
- **Chokepoint** (`_param_symbol` family): rewrite the `_param_symbol` bodies to
  emit structural names from a Step value (and ligand/site for `Kreg`), with the
  A/I state token placed after the type prefix and omitted for `EqualAI`. Delete
  `name(::Type{P}, idx)` and `name(::Type{P}, idx, state)` (the index-context
  companion) and replace `_rep_idx_for_step` with `_rep_step` returning the rep
  `Step`. The single remaining entry point is value-context `name(p::Parameter, m)`.
- **Inactive-state name helper `_inactive_name(sym)`**: returns the structural
  inactive (`I_`) name for a given active-state name by inserting `I_` after the
  first underscore. Every synthesized-dep `Symbol(string(k) * "_T")` site (~11
  across `rate_eq_derivation.jl` and `mechanism_enumeration.jl`) routes through
  this helper so the chokepoint's mid-name token placement and the
  synth-dep machinery use one structural convention.

### Kinetic-group representative (`_step_priority`)

The per-group representative survives (two physically-distinct steps in one
group do **not** share a structural name). Extract a pure
`_step_priority(step, m)::Int` from the existing Haldane pivot logic. The group
name representative is `argmin(_step_priority)` within the group (least-eliminable
= most primary), lexical tiebreak on the structural name for determinism. The
Haldane elimination pivot keeps its `argmax` selection among reps. The group's
parameter name = the rep step's structural name; the Pass-1 rename map maps the
other members' names to it, exactly as today — only keyed on structural names
instead of `K9 => K4`.

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
