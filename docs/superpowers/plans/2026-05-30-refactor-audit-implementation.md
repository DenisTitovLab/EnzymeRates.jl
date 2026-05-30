# Concrete-Types Refactor Audit — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land 66 of the 71 audit findings from `docs/superpowers/2026-05-30-refactor-audit-findings.md` (the 5 in Cluster F are deferred — they require the direction-symmetry refactor to land first). Net target: ~780 LOC saved (~13.7% of the 5,706 non-comment non-doc src baseline) plus ~220 LOC of in-place rewrites that improve clarity without changing LOC.

**Architecture:** Four sequential waves. Wave 1 lands isolated doc/dup cleanups (no deps). Wave 2 rewrites the rate-equation derivation back-end to walk `Mechanism.steps` directly instead of `rxns = reactions(m)` opaque tuples — independent of the singleton demote but a prerequisite for proving Cluster A is non-breaking. Wave 3 is the big architectural move: demote `EnzymeMechanism{Sig}` / `AllostericEnzymeMechanism{...}` to internal compile artifacts (kept only for `@generated rate_equation` specialization), and collapse expansion-move dispatch duals. Wave 4 is the final doc-hygiene sweep.

**Tech Stack:** Julia 1.x package; `Pkg.test()` for the full test suite; `@allocated` / `@btime` for `rate_equation` perf budget; `git` for incremental commits.

**Sacred invariants (from CLAUDE.md, enforced at every commit):**

1. **`rate_equation` performance**: 0 allocations and <100 ns per call across every spec in `MECHANISM_TEST_SPECS`. Enforced by `test_rate_equation_performance` in `test/test_rate_eq_derivation.jl`.
2. **No net test coverage reduction**: tests may be replaced (behavior coverage replacing private-helper coverage) but not net deleted. Log replacements per `docs/superpowers/refactor-deleted-tests.md` convention.
3. **Opaque-form rejection stays**: the user-visible parser semantics that reject `:ES` as a bare enzyme must remain (F-014 Approach B preserves this; Approach C is explicitly out-of-scope unless Denis pre-approves).
4. **Front-end struct-family unification stays**: no finding reintroduces parallel front-end representations.

---

## File Structure

Files modified by this plan:

| File | What changes |
|---|---|
| `src/types.jl` | Most-modified file. Loses Sig conversion bridge (~118 LOC), singleton `@generated` accessors get rewritten to plain functions over `Mechanism` (~258 LOC of `@generated` becomes ~30 LOC of plain), the 14-line forwarding-accessor block collapses, `_rep_step` 3 methods → 1, `_to_mechanism`/`_AnyMech` delete, `_drop_unbound_regulators` relocates, `_species_name_from_sig` + `_step_tuple_from_sig` delete |
| `src/dsl.jl` | Macros emit `Mechanism(...)` / `AllostericMechanism(...)` directly; `_assert_no_opaque_terms` inlines into `_parse_steps_block_with_groups` tail; 19 function-leading `#`-comments convert to docstrings |
| `src/rate_eq_derivation.jl` | `_raw_symbolic_rate_polys`/`_compute_alpha`/`_compute_numerator` walk `Mechanism.steps` directly; `@generated parameters`/`fitted_params` collapse to plain functions; 5 parameter-enumeration helpers gain a shared `_emit_cat_params_for_rep` core; `_build_kinetic_rename_map` renames to `_build_wegscheider_rename_map`; 4 repeated `a_only_syms`/`rename_I` patterns factor into helpers; 3 stale "Stage 4.2" comments delete |
| `src/thermodynamic_constr_for_rate_eq_derivation.jl` | `_free_enz_set`/`_thermodynamic_constraints`/`_dependent_param_exprs_kernel` walk Mechanism directly; stale "K9 => K4" docstring updates; `_raw_param_symbols(::EnzymeMechanism)` forwarder deletes |
| `src/mechanism_enumeration.jl` | `compile_mechanism` lift deletes; three `dedup!` overloads collapse to one; `_canonicalize_for_hash` two methods merge via helpers; `_canonical_rate_eq_hash_data` thin wrapper inlines; `_build_name_map` gets R↔T → A/I rename + Phase 7 comment strip; expansion-move duals collapse via new `_with_*` helpers; `_parameter_canonical_key` 6 methods collapse |
| `src/identify_rate_equation.jl` | Four `compile_mechanism(mech)` lifts delete (pass-through after demote) |
| `src/sym_poly_for_rate_eq_derivation.jl` | Stale `K2 → K1` docstring example updates; minor comment-as-docstring conversions |
| `src/fitting.jl` | Minor: `FittingProblem` constructor signature may retype if `AbstractEnzymeMechanism` dissolves |
| `src/EnzymeRates.jl` | Export list review (post Wave 3); no behavior change |
| `test/test_rate_eq_derivation.jl` | Update tests for `_onlyA_parameters` / `_I_rename_parameters` / `_all_i_state_parameters` to behavior assertions (Wave 2 prerequisite for F-043) |
| `test/test_mechanism_enumeration.jl` | Replace ~50 `compile_mechanism(m)` callers with `m` (identity); update one `name_map isa Dict{String, String}` assertion (Wave 3) |
| `test/test_types.jl` | Drop `_sig_of` / `_mechanism_from_sig` roundtrip test (deletes with the bridge) |
| `test/test_accessors.jl` | Drop or re-baseline the accessor perf gate (Q-005: explicitly negotiable) |
| `docs/superpowers/refactor-deleted-tests.md` | Append entries for the test-blocked-helper replacements in F-043 |

---

# Wave 1 — Doc hygiene + isolated cleanups

Independent of Cluster A and Cluster B. ~70 LOC saved. Land first to validate the workflow.

## Task 1.1: Trim stale "Stage Nx" / "Phase N" / index-based naming references

**Files:**
- Modify: `src/types.jl:1411`
- Modify: `src/rate_eq_derivation.jl:760-761, 1437-1441, 1552-1554`
- Modify: `src/sym_poly_for_rate_eq_derivation.jl:250-258`
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:226-236`
- Modify: `src/mechanism_enumeration.jl:1951-1963, 1980-1986`
- Modify: `src/types.jl:46-47`

This task handles findings F-001, F-012, F-036, F-044, F-051, F-065 in one commit — they're all pure comment edits with no behavior change.

- [ ] **Step 1: types.jl L46-47 — fix Species docstring `:E_A_B` → `:EATP`**

Read `src/types.jl` lines 46-50.

```
46: # §5.3 — Species: an enzyme form. `bound` is sorted by name; the
47: # rendered Symbol name reads `:E` / `:E_A_B` / `:Estar...` / `:E_A_B_res...`.
```

Edit L47 to match the code at L78-82 (which concatenates without separator):

```
47: # rendered Symbol name reads `:E` / `:EATP` / `:Estar...` / `:EATPres_+P`.
```

- [ ] **Step 2: types.jl L1411 — trim "Stage 7a" reference**

Read `src/types.jl` lines 1408-1413.

```
1408: # Regulator-site parameter: state tag + ligand name + "reg". No site index —
1409: # the AllostericMechanism constructor enforces that each ligand appears at most
1410: # once across all sites (same-ligand-two-sites collision check added in Stage 7a).
```

Edit L1410 to drop "(same-ligand-two-sites collision check added in Stage 7a)":

```
1408: # Regulator-site parameter: state tag + ligand name + "reg". No site index —
1409: # the AllostericMechanism constructor enforces that each ligand appears at most
1410: # once across all sites.
```

- [ ] **Step 3: rate_eq_derivation.jl L760-761 — trim "Stage 4.2"**

Read `src/rate_eq_derivation.jl` lines 758-765.

```
760:    # Pass 2: dep RHSes referencing a `:NonequalAI` symbol pick up an
761:    # I-state name. Mirrors `_dependent_param_exprs` Stage 4.2 logic.
```

Edit to drop "Stage 4.2":

```
760:    # Pass 2: dep RHSes referencing a `:NonequalAI` symbol pick up an
761:    # I-state name. Mirrors `_dependent_param_exprs` Pass 2 below.
```

- [ ] **Step 4: rate_eq_derivation.jl L1437-1441 — trim "Stage 4.2"**

Read `src/rate_eq_derivation.jl` lines 1435-1445.

```
1437:    # Pass 2: dep RHSes referencing a `:NonequalAI` symbol need their own
1438:    # I-state name. Mirrors `_dependent_param_exprs` Stage 4.2 logic; the
1439:    # synthesized entries here are the same ones the rate-equation body
1440:    # consumes via `i_names_set`.
```

Edit:

```
1437:    # Pass 2: dep RHSes referencing a `:NonequalAI` symbol need their own
1438:    # I-state name. Mirrors `_dependent_param_exprs` Pass 2; the
1439:    # synthesized entries here are the same ones the rate-equation body
1440:    # consumes via `i_names_set`.
```

- [ ] **Step 5: rate_eq_derivation.jl L1552-1554 — trim "Stage 4.2"**

Read `src/rate_eq_derivation.jl` lines 1550-1558.

```
1552:    # Pass 2: dep RHSes referencing a `:NonequalAI` symbol need their own
1553:    # I-state name so the polynomial rename covers synthesized deps.
1554:    # Mirrors the second pass in `_dependent_param_exprs` (Stage 4.2).
```

Edit:

```
1552:    # Pass 2: dep RHSes referencing a `:NonequalAI` symbol need their own
1553:    # I-state name so the polynomial rename covers synthesized deps.
1554:    # Mirrors the second pass in `_dependent_param_exprs`.
```

- [ ] **Step 6: sym_poly_for_rate_eq_derivation.jl L250-258 — fix stale `K2 → K1` example**

Read `src/sym_poly_for_rate_eq_derivation.jl` lines 250-260.

```
253: """
254: Rename symbols in a polynomial. `rename_map` is a `Dict{Symbol, Symbol}`;
255: absent keys are left unchanged. Used to alias non-representative kinetic-group
256: parameter symbols to their representative (e.g., `K2 → K1` when steps 1 and 2
257: share a kinetic group).
258: """
```

Edit the docstring to reflect current usage (A→I state rename in allosteric derivation):

```
253: """
254: Rename symbols in a polynomial. `rename_map` is a `Dict{Symbol, Symbol}`;
255: absent keys are left unchanged. Used by the allosteric derivation to rename
256: A-state symbols to their I-state counterparts when building the inactive-
257: state polynomial (e.g., `:K_A_ATP_E → :K_I_ATP_E`).
258: """
```

- [ ] **Step 7: thermodynamic_constr.jl L226-236 — fix stale `K9 => K4` example**

Read `src/thermodynamic_constr_for_rate_eq_derivation.jl` lines 224-240.

The comment block at L224-236 references "`K9 => K4`" — old index-based naming. Edit the entire `# Filter Pass-2-absorbed symbols...` block to use structural-name example:

```
224:    # Filter Pass-2-absorbed symbols out of indep. Pass 2 of
225:    # `_build_wegscheider_rename_map` adds entries like `K_P_E => K_S_E`
226:    # when a Wegscheider tie collapses two binding-K group reps to the
227:    # same name. After the merge, the absorbed symbol doesn't appear in
228:    # the v polynomial — its column has been folded into the target.
229:    # But the absorbed symbol is still a kinetic-group rep in the
230:    # mechanism, so `_raw_param_symbols` emits it and the kernel keeps it
231:    # in `indep`. Without this filter, `fitted_params` exposes a fittable
232:    # dummy dimension that doesn't affect the loss, and finite-restart
233:    # convergence suffers (the same rate equation can land at noticeably
234:    # different fitted losses depending on which absorbed symbol got
235:    # the dummy slot).
236:    indep = Tuple(p for p in indep if get(rename, p, p) == p)
```

Note: the rename `_build_kinetic_rename_map` → `_build_wegscheider_rename_map` referenced here happens in Task 1.6; this step pre-emptively uses the new name so both edits land coherent.

- [ ] **Step 8: mechanism_enumeration.jl L1951-1963 — fix stale R↔T notation in `_build_name_map` docstring**

Read `src/mechanism_enumeration.jl` lines 1951-1965.

```
1951: """
1952: Build the per-mechanism Symbol → canonical-token map. Used both by the
1953: canonical-form construction (substitutes Symbols in POLYs / Exprs) and
1954: returned through `_canonical_rate_eq_hash_data` for downstream
1955: projection via `_project_cached_params`.
1956:
1957: For an `AllostericMechanism`, also adds entries for synthesized dep
1958: T-names (LHSes that have no Parameter struct because they're derived
1959: deps with a `_T` suffix appended at render time). The synth-dep token
1960: is the R-state token with `_T` suffix, preserving R↔T correspondence
1961: across equivalent mechanisms.
1962: """
```

Edit L1958-1961 to use A/I notation (the refactor's current internal taxonomy):

```
1957: For an `AllostericMechanism`, also adds entries for synthesized dep
1958: I-names (LHSes that have no Parameter struct because they're derived
1959: deps with an `_I` suffix appended at render time). The synth-dep token
1960: is the A-state token with `_I` suffix, preserving A↔I correspondence
1961: across equivalent mechanisms.
1962: """
```

- [ ] **Step 9: mechanism_enumeration.jl L1980-1986 — rename r/t variables + drop "Phase 7" comment**

Read `src/mechanism_enumeration.jl` lines 1979-1990.

```
1979:    if m isa AllostericMechanism
1980:        for r_name in _synth_dep_a_names(em, m)
1981:            r_str = String(r_name)
1982:            tok = get(name_map, r_str, nothing)
1983:            tok === nothing && continue
1984:            t_str = String(name(_flip_to_inactive(_param_for_symbol(m, r_name)), m))
1985:            haskey(name_map, t_str) && continue
1986:            name_map[t_str] = tok * "_T"  # tok is a canonical p_<i> hash token, not a parameter Symbol — Phase 7 cleanup territory.
1987:        end
1988:    end
1989:    name_map
1990: end
```

Edit to rename `r_name`/`r_str`/`t_str` → `a_name`/`a_str`/`i_str` and drop the "Phase 7" comment:

```
1979:    if m isa AllostericMechanism
1980:        for a_name in _synth_dep_a_names(em, m)
1981:            a_str = String(a_name)
1982:            tok = get(name_map, a_str, nothing)
1983:            tok === nothing && continue
1984:            i_str = String(name(_flip_to_inactive(_param_for_symbol(m, a_name)), m))
1985:            haskey(name_map, i_str) && continue
1986:            name_map[i_str] = tok * "_T"
1987:        end
1988:    end
1989:    name_map
1990: end
```

- [ ] **Step 10: Run the test suite to confirm no behavior change**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green. These are pure comment / variable-rename edits; the variable rename inside `_build_name_map` is local (no external callers see those names).

If any test fails with a name error referencing `r_name`/`r_str`/`t_str`, the variable shadowing isn't fully internal — find the offending caller and resolve before committing.

- [ ] **Step 11: Commit**

```bash
git add src/types.jl src/rate_eq_derivation.jl src/sym_poly_for_rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl src/mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
docs: trim stale stage/phase references + index-based naming examples

Findings F-001, F-012, F-036, F-044, F-051, F-065 from the
2026-05-30 audit. Pure comment / variable-rename edits — no behavior
change.

- Species docstring example matches concatenated rendering (E_A_B → EATP)
- Drop "Stage 7a" reference at types.jl:1411
- Drop 3 "Stage 4.2" references in rate_eq_derivation.jl synth-dep code
- Update _rename_symbols docstring example from K2→K1 (old index-based)
  to A→I state rename (current usage)
- Update _dependent_param_exprs filter comment from K9→K4 to structural
  K_P_E → K_S_E example
- Update _build_name_map docstring R↔T → A/I notation
- Rename r_name/r_str/t_str → a_name/a_str/i_str inside _build_name_map
- Drop "Phase 7 cleanup territory" comment

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 1.2: Rename `_build_kinetic_rename_map` → `_build_wegscheider_rename_map`

**Files:**
- Modify: `src/rate_eq_derivation.jl:94-141`
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:218-223` (the call site)
- Modify: `src/rate_eq_derivation.jl:416` (the second call site, inside `_raw_symbolic_rate_polys`)

This is finding F-038. Verifies Q-010's conclusion that the function only does Wegscheider RE-tie absorption.

- [ ] **Step 1: Rename the function and trim stale docstring**

Read `src/rate_eq_derivation.jl` lines 93-141.

Replace the docstring (L95-111) and function (L112-141) with the renamed version:

```julia
"""
Build a renaming map for single-symbol Wegscheider RE ties between two
binding K's. Calls `_dependent_param_exprs_kernel` to discover
binding-K-to-binding-K Wegscheider closures of the form `K_a = K_b`
(RHS is a bare Symbol). Both sides must be binding K's (RE step with
metabolite on LHS) — absorbing a binding-K-to-iso-K tie would produce
inconsistent sign-flips when the kernel runs with the full rename,
since the binding-K column is sign-flipped (Kd convention) but the
iso-K column is not.

The rename means the polynomial in `v` uses the representative symbol
directly, so Source-C duplicates (split kinetic groups that
Wegscheider ties back together) collapse at hash time.
"""
function _build_wegscheider_rename_map(M::Type{<:EnzymeMechanism})
    m = M()
    mech = Mechanism(m)
    rename = Dict{Symbol, Symbol}()
    rxns = reactions(m)
    eq = equilibrium_steps(m)
    enz_set = Set(enzyme_forms(m))
    step_params = _step_parameters(mech)
    binding_set = Set{Symbol}()
    for (idx, (lhs, _, _, _)) in enumerate(rxns)
        eq[idx] || continue
        any(s ∉ enz_set for s in lhs) || continue
        push!(binding_set, name(step_params[idx][1], mech))
    end
    dep_raw, _ = _dependent_param_exprs_kernel(M, rename)
    for (lhs, rhs) in dep_raw
        rhs isa Symbol || continue
        lhs in binding_set && rhs in binding_set || continue
        target = get(rename, rhs, rhs)
        rename[lhs] = target
        for k in collect(keys(rename))
            rename[k] == lhs && (rename[k] = target)
        end
    end
    rename
end

_build_wegscheider_rename_map(m::EnzymeMechanism) =
    _build_wegscheider_rename_map(typeof(m))
```

- [ ] **Step 2: Update the call site in thermodynamic_constr.jl:223**

Read `src/thermodynamic_constr_for_rate_eq_derivation.jl` lines 220-225.

```
222: function _dependent_param_exprs(M::Type{<:EnzymeMechanism})
223:     rename = _build_kinetic_rename_map(M)
```

Edit L223:

```
223:     rename = _build_wegscheider_rename_map(M)
```

- [ ] **Step 3: Update the call site in rate_eq_derivation.jl:416**

Read `src/rate_eq_derivation.jl` lines 414-418.

```
416:     rename_map = _build_kinetic_rename_map(m)
```

Edit:

```
416:     rename_map = _build_wegscheider_rename_map(m)
```

- [ ] **Step 4: grep for any remaining references to the old name**

```bash
grep -rn "_build_kinetic_rename_map" src/ test/
```

Expected output: empty.

If any reference remains, fix it.

- [ ] **Step 5: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green. Same logic, different name.

- [ ] **Step 6: Commit**

```bash
git add src/rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
refactor: rename _build_kinetic_rename_map to _build_wegscheider_rename_map

Finding F-038 from the 2026-05-30 audit. The function only handles
single-symbol Wegscheider RE-tie absorption now — the "kinetic-group
merges" path is gone (the value-context chokepoint handles group
collapse directly). New name reflects the actual responsibility.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 1.3: Factor repeated `a_only_syms` / `rename_I` construction patterns

**Files:**
- Modify: `src/rate_eq_derivation.jl` (add helpers near L985, update 4 call sites)

Findings F-045 (a_only_syms set, 4× repeat) and F-046 (rename_I Dict, 4× repeat).

- [ ] **Step 1: Add two helpers above `_onlyA_parameters` (around line 985)**

Read `src/rate_eq_derivation.jl` lines 984-995.

Insert before `_onlyA_parameters`:

```julia
"""
Set of A-state catalytic parameter Symbol names for an
`AllostericMechanism`. Cached helper for the four call sites in
`_kcat_forward`, `_dependent_param_exprs`, `_build_dep_assignments`,
and `_allosteric_num_den_exprs` that each previously rebuilt this set.
"""
_a_only_syms(am::AllostericMechanism) =
    Set{Symbol}(name(p, am) for p in _onlyA_parameters(am))

"""
A → I rename map (Symbol → Symbol) for `:NonequalAI` catalytic-group
parameters. Routes through `name(p, am)`; both keys and values are the
rendered Symbol names of `Kd/Kiso/Kon/Koff/Kfor/Krev` parameters.
"""
_a_to_i_rename(am::AllostericMechanism) =
    Dict{Symbol, Symbol}(
        name(p_A, am) => name(p_I, am)
        for (p_A, p_I) in _I_rename_parameters(am))
```

- [ ] **Step 2: Replace the 4 inline `a_only_syms = ...` constructions**

For each of L756, L1283, L1432, L1547 in `src/rate_eq_derivation.jl`:

Find the line of the form:
```julia
    a_only_syms = Set{Symbol}(name(p, am) for p in _onlyA_parameters(am))
```

Replace with:
```julia
    a_only_syms = _a_only_syms(am)
```

Run grep to verify each site is replaced:

```bash
grep -n "a_only_syms = Set" src/rate_eq_derivation.jl
```

Expected output: empty (4 → 0 inline constructions).

- [ ] **Step 3: Replace the 4 inline `rename_I = ...` constructions**

For each of L757-759, L1289-1292, L1434-1436, L1549-1551 in `src/rate_eq_derivation.jl`:

Find the 3-line block of the form:
```julia
    rename_I = Dict{Symbol, Symbol}(
        name(p_A, am) => name(p_I, am)
        for (p_A, p_I) in _I_rename_parameters(am))
```

Replace with:
```julia
    rename_I = _a_to_i_rename(am)
```

Run grep to verify:

```bash
grep -n "rename_I = Dict" src/rate_eq_derivation.jl
```

Expected output: empty.

- [ ] **Step 4: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green. Pure factoring — identical inputs, identical outputs.

- [ ] **Step 5: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
refactor: factor repeated a_only_syms / rename_I patterns

Findings F-045 / F-046 from the 2026-05-30 audit. Four identical
inline constructions of each replaced by single-line calls to new
helpers _a_only_syms(am) and _a_to_i_rename(am). ~12 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 1.4: Inline `_canonical_rate_eq_hash_data` thin wrapper

**Files:**
- Modify: `src/mechanism_enumeration.jl:2037-2064`
- Modify: `src/identify_rate_equation.jl:427`
- Modify: `test/test_mechanism_enumeration.jl:4291, 4341, 4395, 4455, 4456`

Finding F-067.

- [ ] **Step 1: Rename the impl function to drop the suffix**

Read `src/mechanism_enumeration.jl` lines 2037-2064.

```
2046: function _canonical_rate_eq_hash_data_impl_struct(em::AbstractEnzymeMechanism)
```

Edit L2046 to drop `_impl_struct`:

```
2046: function _canonical_rate_eq_hash_data(em::AbstractEnzymeMechanism)
```

Update the docstring at L2037-2045 to reflect the merge (this is now the only entry point — drop the "implementation of" framing):

```
"""
Compute the canonical rate-equation hash for `em`. Walks
`Mechanism` / `AllostericMechanism` structural fields directly via
`_canonicalize_for_hash`. Returns `(UInt64 hash, 16-char hex display
string, name_map)`. The `name_map::Dict{String,String}` satisfies the
projection contract used by `_project_cached_params`: two
hash-equivalent mechanisms produce maps that send corresponding
parameter Symbols to the same canonical token.

Hash collision probability over 10⁴ mechanisms is ~10⁻¹² with
Julia's built-in `hash(::UInt64)::UInt64`.
"""
```

- [ ] **Step 2: Delete the old wrapper at L2062-2064**

Read `src/mechanism_enumeration.jl` lines 2052-2070.

Delete the 3-line wrapper function and its docstring:

```
2052: """
2053: Return `(UInt64 hash, 16-char hex display string, name_map)`.
2054: The single entry point for canonical hashing; `_canonical_rate_eq_hash`
2055: delegates here so the canonicalizer runs once and callers that need the
2056: name_map can retrieve it without a second pass.
...
2061: """
2062: function _canonical_rate_eq_hash_data(m::AbstractEnzymeMechanism)
2063:     _canonical_rate_eq_hash_data_impl_struct(m)
2064: end
```

After this edit, only the (renamed) impl function exists.

- [ ] **Step 3: grep for any remaining references to the old impl name**

```bash
grep -rn "_canonical_rate_eq_hash_data_impl_struct" src/ test/
```

Expected output: empty. If any remain, update them.

- [ ] **Step 4: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green. The 6 callers (1 prod + 5 test) all dispatch on the public name, which now resolves to the impl directly.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
refactor: inline _canonical_rate_eq_hash_data thin wrapper

Finding F-067 from the 2026-05-30 audit. The 3-line wrapper function
existed for legacy multi-impl dispatch which is gone. Rename
_canonical_rate_eq_hash_data_impl_struct to drop the suffix; delete
the wrapper.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 1.5: Merge three `dedup!` overloads into one

**Files:**
- Modify: `src/mechanism_enumeration.jl:1780-1822`

Finding F-060. Q-014 confirmed this is safe.

- [ ] **Step 1: Replace the three `dedup!` methods with one**

Read `src/mechanism_enumeration.jl` lines 1780-1822.

Replace the entire block (L1781-1822, three function definitions) with a single method:

```julia
"""
    dedup!(cache::Dict{Int, <:Vector})

Canonicalize each mechanism in place via the type-specific
`_canonicalize_mechanism!` overload, then run `unique!` so
structurally-equivalent mechanisms collapse. Works for any element
type — `Mechanism`, `AllostericMechanism`, or
`Union{Mechanism, AllostericMechanism}` — because
`_canonicalize_mechanism!` dispatches at runtime.
"""
function dedup!(cache::Dict{Int, <:Vector})
    for (pc, mechs) in cache
        for m in mechs
            _canonicalize_mechanism!(m)
        end
        unique!(mechs)
        isempty(mechs) && delete!(cache, pc)
    end
    cache
end
```

- [ ] **Step 2: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green. All 3 cache shapes (Vector{Mechanism}, Vector{AllostericMechanism}, Vector{Union}) are subtypes of `Vector`, so dispatch finds the merged method.

- [ ] **Step 3: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
refactor: merge three dedup! overloads into one

Finding F-060 from the 2026-05-30 audit. Per Q-014, the three methods
(Vector{Mechanism}, Vector{AllostericMechanism}, Vector{Union}) had
byte-identical bodies. Collapse to one method on Dict{Int, <:Vector}
since _canonicalize_mechanism! dispatches at runtime for any element
type. ~22 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 1.6: Merge `_canonicalize_for_hash` two methods via helpers

**Files:**
- Modify: `src/mechanism_enumeration.jl:1887-1950`

Finding F-063. Q-019 confirmed feasible; canonical-hash partition test is the regression guard.

- [ ] **Step 1: Add `_num_den_exprs` and `_extra_canon_tuple` helpers above the merged function**

Read `src/mechanism_enumeration.jl` lines 1885-1955.

Insert before `_canonicalize_for_hash` (replacing the existing two-method definition):

```julia
"""
Type-specific (num_expr, den_expr) for the canonical-hash trunk.
- For non-allosteric: builds POLYs via `_raw_symbolic_rate_polys` and
  renders via `_poly_to_expr`.
- For allosteric: uses `_allosteric_num_den_exprs` which returns the
  full MWC Exprs directly.
"""
function _num_den_exprs(em::AbstractEnzymeMechanism, ::Mechanism)
    M = typeof(em)
    num, den = _raw_symbolic_rate_polys(M)
    pset = Set{Symbol}(_raw_param_symbols(em))
    cset = Set{Symbol}(metabolites(em))
    (_poly_to_expr(num, pset, cset), _poly_to_expr(den, pset, cset))
end

function _num_den_exprs(em::AbstractEnzymeMechanism, ::AllostericMechanism)
    em isa AllostericEnzymeMechanism || error(
        "_num_den_exprs: AllostericMechanism requires " *
        "AllostericEnzymeMechanism, got $(typeof(em))")
    _allosteric_num_den_exprs(typeof(em))
end

"""
Type-specific canon-tuple prefix. Non-allosteric gets just the
`(:NonAllosteric,)` tag; allosteric appends `cat_tags_canon`,
`cat_mult`, and `site_canon` so two allosteric mechanisms differing
only in those scalars hash distinctly.
"""
_extra_canon_tuple(::Mechanism) = ((:NonAllosteric,),)

function _extra_canon_tuple(m::AllostericMechanism)
    cat_tags_canon = Tuple(cat_allo_states(m))
    cat_mult = catalytic_multiplicity(m)
    site_entries = Tuple[]
    for site in regulatory_sites(m)
        push!(site_entries,
              (Tuple(hash(l) for l in ligands(site)),
               multiplicity(site),
               Tuple(allo_states(site))))
    end
    site_canon = Tuple(sort(site_entries; by = repr))
    ((:Allosteric,), cat_tags_canon, cat_mult, site_canon)
end
```

- [ ] **Step 2: Replace both `_canonicalize_for_hash` methods with one merged trunk**

Continue replacing in `src/mechanism_enumeration.jl` (the L1887-1950 region — both methods).

After the helpers from Step 1, write the single merged method:

```julia
"""
Canonical form for a mechanism plus its `name_map::Dict{String,String}`.
Walks the Parameter family + symbolic numerator/denominator Exprs
directly, producing a canonical key per Parameter from `Step` /
`RegulatorySite` / `AllostericRegulator` identity rather than from
rendered symbol strings.

Two mechanisms with the same rate equation but different kinetic-group
numbering (and therefore different positional symbol names) produce the
same canonical Expr tree because their Parameter canonical keys
coincide and the `_poly_to_expr` monomial sort agrees once substitution
is applied. For allosteric mechanisms, additional canon slots
(catalytic state tags, catalytic multiplicity, regulator site shape)
ensure mechanisms differing only in those scalars hash distinctly.
"""
function _canonicalize_for_hash(em::AbstractEnzymeMechanism,
                                m::Union{Mechanism, AllostericMechanism})
    name_map = _build_name_map(em, m)
    dep_canon = _dep_exprs_canonical(em, name_map)
    num, den = _num_den_exprs(em, m)
    num_canon = _expr_canonical_via_name_map(num, name_map)
    den_canon = _expr_canonical_via_name_map(den, name_map)
    canon = (_extra_canon_tuple(m)..., num_canon, den_canon, dep_canon)
    (canon, name_map)
end
```

**CRITICAL:** the canon-tuple shape MUST match the old code byte-for-byte. Verify:

- Old non-allosteric: `canon = ((:NonAllosteric,), num_canon, den_canon, dep_canon)`. New: `(_extra_canon_tuple(::Mechanism)..., num_canon, den_canon, dep_canon)` = `((:NonAllosteric,), num_canon, den_canon, dep_canon)`. ✓ matches.
- Old allosteric: `canon = ((:Allosteric,), num_canon, den_canon, cat_tags_canon, cat_mult, site_canon, dep_canon)`. New: `(_extra_canon_tuple(::AllostericMechanism)..., num_canon, den_canon, dep_canon)` = `((:Allosteric,), cat_tags_canon, cat_mult, site_canon, num_canon, den_canon, dep_canon)`. **DIFFERS** — the old layout interleaved `num_canon`/`den_canon` BEFORE the allosteric extras.

Fix: keep the old layout exactly. Restructure `_extra_canon_tuple` to be a postfix, not prefix:

```julia
"""
Type-specific canon-tuple suffix for allosteric mechanisms. Empty
tuple for non-allosteric; for allosteric, includes catalytic state
tags, multiplicity, and regulator site shape so two allosteric
mechanisms differing only in those scalars hash distinctly.
"""
_allosteric_canon_suffix(::Mechanism) = ()

function _allosteric_canon_suffix(m::AllostericMechanism)
    cat_tags_canon = Tuple(cat_allo_states(m))
    cat_mult = catalytic_multiplicity(m)
    site_entries = Tuple[]
    for site in regulatory_sites(m)
        push!(site_entries,
              (Tuple(hash(l) for l in ligands(site)),
               multiplicity(site),
               Tuple(allo_states(site))))
    end
    site_canon = Tuple(sort(site_entries; by = repr))
    (cat_tags_canon, cat_mult, site_canon)
end

_canon_tag(::Mechanism) = (:NonAllosteric,)
_canon_tag(::AllostericMechanism) = (:Allosteric,)
```

And the merged trunk becomes:

```julia
function _canonicalize_for_hash(em::AbstractEnzymeMechanism,
                                m::Union{Mechanism, AllostericMechanism})
    name_map = _build_name_map(em, m)
    dep_canon = _dep_exprs_canonical(em, name_map)
    num, den = _num_den_exprs(em, m)
    num_canon = _expr_canonical_via_name_map(num, name_map)
    den_canon = _expr_canonical_via_name_map(den, name_map)
    canon = m isa AllostericMechanism ?
        ((:Allosteric,), num_canon, den_canon,
         _allosteric_canon_suffix(m)..., dep_canon) :
        ((:NonAllosteric,), num_canon, den_canon, dep_canon)
    (canon, name_map)
end
```

This preserves the old canon layout exactly. The non-allosteric branch reads `((:NonAllosteric,), num_canon, den_canon, dep_canon)` (4 slots). The allosteric branch reads `((:Allosteric,), num_canon, den_canon, cat_tags_canon, cat_mult, site_canon, dep_canon)` (7 slots — splatting the suffix). Both match the old shapes.

- [ ] **Step 3: Run the canonical-hash partition test as a regression guard**

```bash
julia --project test/test_canonical_hash_partition.jl
```

Expected: green. If the canon-tuple shape changed, this test catches it.

- [ ] **Step 4: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
refactor: merge two _canonicalize_for_hash methods via helpers

Finding F-063 from the 2026-05-30 audit. Per Q-019, the two methods
(::Mechanism and ::AllostericMechanism) share a 4-step trunk
(name_map → dep_canon → num/den_canon → canon assembly) and differ
only in num/den source and allosteric-extras. Extracted as
_num_den_exprs (dispatch on Mechanism vs AllostericMechanism) and
_allosteric_canon_suffix (returns extra slots for allosteric only) +
_canon_tag for the leading tag.

Canon-tuple shape preserved byte-for-byte; test_canonical_hash_partition.jl
acts as the regression guard.

~25 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 1.7: Collapse 6 `_parameter_canonical_key` methods to one

**Files:**
- Modify: `src/mechanism_enumeration.jl:1845-1855`

Finding F-061.

- [ ] **Step 1: Replace 6 one-liner methods with one type-dispatching method**

Read `src/mechanism_enumeration.jl` lines 1840-1860.

Replace:

```julia
_parameter_canonical_key(p::Kd)   = (:Kd,   hash(p.step), p.state)
_parameter_canonical_key(p::Kiso) = (:Kiso, hash(p.step), p.state)
_parameter_canonical_key(p::Kon)  = (:Kon,  hash(p.step), p.state)
_parameter_canonical_key(p::Koff) = (:Koff, hash(p.step), p.state)
_parameter_canonical_key(p::Kfor) = (:Kfor, hash(p.step), p.state)
_parameter_canonical_key(p::Krev) = (:Krev, hash(p.step), p.state)
_parameter_canonical_key(p::Kreg) =
    (:Kreg, hash(p.site), hash(p.ligand), p.state)
_parameter_canonical_key(::Keq)   = (:Keq,)
_parameter_canonical_key(::Etot)  = (:Etot,)
_parameter_canonical_key(::Lallo) = (:Lallo,)
```

With:

```julia
_parameter_canonical_key(p::StepBoundParameter) =
    (nameof(typeof(p)), hash(p.step), p.state)
_parameter_canonical_key(p::Kreg) =
    (:Kreg, hash(p.site), hash(p.ligand), p.state)
_parameter_canonical_key(::Keq)   = (:Keq,)
_parameter_canonical_key(::Etot)  = (:Etot,)
_parameter_canonical_key(::Lallo) = (:Lallo,)
```

`StepBoundParameter` is the existing Union at `src/types.jl:204` covering `Kd, Kiso, Kon, Koff, Kfor, Krev`. `nameof(typeof(p))` returns `:Kd`, `:Kiso`, etc. — same Symbols the explicit methods emitted.

- [ ] **Step 2: Run canonical-hash test as a regression guard**

```bash
julia --project test/test_canonical_hash_partition.jl
```

Expected: green. The key shapes are byte-identical to before.

- [ ] **Step 3: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 4: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
refactor: collapse 6 _parameter_canonical_key methods via StepBoundParameter

Finding F-061 from the 2026-05-30 audit. Kd/Kiso/Kon/Koff/Kfor/Krev
methods all had the shape (Symbol(tag), hash(p.step), p.state). Replace
with one StepBoundParameter method using nameof(typeof(p)). Key
shape byte-identical to before; test_canonical_hash_partition.jl as
regression guard.

~5 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 1.8: Inline `_assert_no_opaque_terms` (Approach B for parser-tighten)

**Files:**
- Modify: `src/dsl.jl:382-399` (the function to inline + delete)
- Modify: `src/dsl.jl:567` (caller 1)
- Modify: `src/dsl.jl:1130` (caller 2)
- Modify: `src/dsl.jl:378-381` (the function-leading `#`-comment block — delete with the function)

Finding F-014 Approach B. Pure naming/inline refactor — semantics unchanged.

- [ ] **Step 1: Capture the function body for inlining**

The function at `src/dsl.jl:382-399` is:

```julia
function _assert_no_opaque_terms(side_terms_per_step)
    call_heads = Set{Symbol}()
    for (_, lhs, rhs, _) in side_terms_per_step
        for t in (lhs..., rhs...)
            t.kind === :call && push!(call_heads, t.conformation)
        end
    end
    for (_, lhs, rhs, _) in side_terms_per_step
        for t in (lhs..., rhs...)
            t.kind === :bare_enzyme || continue
            (t.sym in call_heads || _is_conformation_shape(t.sym)) && continue
            error("@enzyme_mechanism: `$(t.sym)` looks like an opaque " *
                  "bound-form name; write it as decomposed call notation, " *
                  "e.g. `E(S)` or `E(A, B)`.")
        end
    end
    nothing
end
```

Plus the explanatory comment at L378-381.

- [ ] **Step 2: Inline the body at the first caller (dsl.jl:567)**

Read `src/dsl.jl` lines 564-572.

Replace `_assert_no_opaque_terms(side_terms_per_step)` at L567 with the inlined body:

```julia
    side_terms_per_step =
        _parse_steps_block_with_groups(steps_block, declared_mets)

    # Reject opaque bound-form bare-enzyme names: a bare-enzyme term `:X` is
    # acceptable iff `:X` is a Call-form head seen in this steps block
    # (`E` in `E(S)`) OR matches the single-cap-then-lower conformation
    # shape (`:E`, `:Estar`, `:E_c`). Multi-capital Symbols (`:ES`, `:EAB`)
    # and underscore-then-uppercase Symbols (`:E_S`) are opaque and rejected
    # in favor of decomposed call notation.
    let
        call_heads = Set{Symbol}()
        for (_, lhs, rhs, _) in side_terms_per_step
            for t in (lhs..., rhs...)
                t.kind === :call && push!(call_heads, t.conformation)
            end
        end
        for (_, lhs, rhs, _) in side_terms_per_step
            for t in (lhs..., rhs...)
                t.kind === :bare_enzyme || continue
                (t.sym in call_heads || _is_conformation_shape(t.sym)) && continue
                error("@enzyme_mechanism: `$(t.sym)` looks like an opaque " *
                      "bound-form name; write it as decomposed call notation, " *
                      "e.g. `E(S)` or `E(A, B)`.")
            end
        end
    end
    return _build_mechanism_expr(subs_list, prods_list, regs_list,
                                 role_of, side_terms_per_step)
```

The `let` block isolates the `call_heads` binding from the surrounding scope.

- [ ] **Step 3: Inline the body at the second caller (dsl.jl:1130)**

Read `src/dsl.jl` lines 1126-1135.

Replace `_assert_no_opaque_terms(side_terms_per_step)` at L1130 with the same `let` block as Step 2.

- [ ] **Step 4: Delete the function definition and its leading comment**

Delete `src/dsl.jl` lines 378-399 (the comment block + function). The lines around L378-399 are:

```
378: # Raise a clear migration error if any bare-enzyme term is an opaque
379: # bound-form name. A bare-enzyme term `:X` is acceptable iff `:X` is a
380: # Call-form head seen in this steps block (`E` in `E(S)`) OR matches the
381: # single-cap-then-lower conformation shape (`:E`, `:Estar`, `:E_c`).
382: function _assert_no_opaque_terms(side_terms_per_step)
383: ...
399: end
```

Delete L378-399 entirely. `_is_conformation_shape` at L375-376 stays (still used by the inlined logic).

- [ ] **Step 5: grep for remaining references**

```bash
grep -rn "_assert_no_opaque_terms" src/ test/
```

Expected: empty. If any remain, fix them.

- [ ] **Step 6: Run dsl tests as regression guard**

```bash
julia --project test/test_dsl.jl
```

Expected: green. The opaque-rejection error message is unchanged (same `error(...)` string).

- [ ] **Step 7: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 8: Commit**

```bash
git add src/dsl.jl
git commit -m "$(cat <<'EOF'
refactor: inline _assert_no_opaque_terms into its two callers

Finding F-014 (Approach B) from the 2026-05-30 audit. Per Q-013, the
function is a post-hoc walk that fires immediately after the parse
step. Inline into the tails of _parse_plain_mechanism_body (L567) and
_parse_allosteric_mechanism_body (L1130), drop the named identity.
Semantics unchanged — same error message, same trigger conditions.

~25 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 1.9: Wave 1 verification

- [ ] **Step 1: Run full test suite end-to-end**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 2: Measure LOC delta**

```bash
for f in src/*.jl; do
    total=$(wc -l < "$f")
    nccnd=$(awk '
        BEGIN { in_ds=0 }
        {
            line=$0
            n = gsub(/"""/, "\"\"\"", line)
            stripped = line
            sub(/^[[:space:]]+/, "", stripped)
            sub(/[[:space:]]+$/, "", stripped)
            if (in_ds) {
                if (n % 2 == 1) in_ds = 0
                next
            }
            if (n % 2 == 1) { in_ds = 1; next }
            if (stripped == "") next
            if (substr(stripped, 1, 1) == "#") next
            count++
        }
        END { print count+0 }
    ' "$f")
    printf "%6d  %6d  %s\n" "$total" "$nccnd" "$f"
done
echo "---"
awk '{ tot+=$1; ncc+=$2 } END { printf "TOTAL: %d total LOC, %d non-comment non-doc LOC\n", tot, ncc }' < <(
    for f in src/*.jl; do
        total=$(wc -l < "$f")
        nccnd=$(awk '
            BEGIN { in_ds=0 }
            { line=$0; n=gsub(/"""/,"\"\"\"",line); stripped=line; sub(/^[[:space:]]+/,"",stripped); sub(/[[:space:]]+$/,"",stripped)
              if (in_ds) { if (n%2==1) in_ds=0; next }
              if (n%2==1) { in_ds=1; next }
              if (stripped=="") next; if (substr(stripped,1,1)=="#") next; count++ }
            END { print count+0 }' "$f")
        echo "$total $nccnd"
    done
)
```

Expected: ~5,635 non-comment non-doc LOC (down from baseline 5,706 by ~70).

Record the new number in your scratchpad — Wave 4 uses it as the post-impl baseline.

# Wave 2 — Derivation back-end struct-native walk

Cluster B + C. ~120 LOC saved. Independent of Wave 3 but recommended as a stepping stone — verifying the back-end refactor leaves all tests green builds confidence for the larger Wave 3.

## Task 2.1: Extract `_emit_cat_params_for_rep` shared helper

**Files:**
- Modify: `src/rate_eq_derivation.jl` (add helper, refactor 3 of 5 helpers to use it)

Finding F-043 (extract; not full-unification per Q-015).

- [ ] **Step 1: Add the helper above `_onlyA_parameters` (around L985, where the `_a_only_syms` / `_a_to_i_rename` helpers from Task 1.3 already live)**

```julia
"""
Emit the Parameter(s) governing a single kinetic-group representative
step with the given allosteric state. The 4-way switch on
`is_equilibrium(rep)` × `is_binding(rep)` is the shared core that
`_onlyA_parameters`, `_all_i_state_parameters`,
`_enumerate_parameters_full_allosteric`, and `_ss_rate_constant_names`
all duplicate. Centralizing here.

Returns 1 element for RE steps (`Kd` or `Kiso`) and 2 elements for SS
steps (`Kon`+`Koff` or `Kfor`+`Krev`).
"""
function _emit_cat_params_for_rep(rep::Step, state::Symbol)
    if is_equilibrium(rep)
        return Parameter[is_binding(rep) ? Kd(rep, state) : Kiso(rep, state)]
    end
    if is_binding(rep)
        return Parameter[Kon(rep, state), Koff(rep, state)]
    end
    Parameter[Kfor(rep, state), Krev(rep, state)]
end
```

- [ ] **Step 2: Refactor `_onlyA_parameters` (L986-1011) to use the helper**

Read `src/rate_eq_derivation.jl` lines 986-1011.

Replace the function body:

```julia
function _onlyA_parameters(am::AllostericMechanism)
    out = Parameter[]
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        cat_allo_state(am, g) === :OnlyA || continue
        append!(out, _emit_cat_params_for_rep(_group_rep(group, fes), :A))
    end
    out
end
```

- [ ] **Step 3: Refactor `_all_i_state_parameters` (L1114-1148)**

Read `src/rate_eq_derivation.jl` lines 1114-1148.

Replace the function body (preserving the reg-site appendix at the end):

```julia
function _all_i_state_parameters(am::AllostericMechanism)
    out = Parameter[]
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        cat_allo_state(am, g) === :OnlyA && continue
        append!(out, _emit_cat_params_for_rep(_group_rep(group, fes), :I))
    end
    for site in regulatory_sites(am)
        for (lig, tag) in zip(ligands(site), allo_states(site))
            tag === :OnlyA && continue
            push!(out, Kreg(site, lig, :I))
        end
    end
    out
end
```

- [ ] **Step 4: Refactor `_enumerate_parameters_full_allosteric` (L1175-1202)**

Read `src/rate_eq_derivation.jl` lines 1175-1202.

Replace the function body:

```julia
function _enumerate_parameters_full_allosteric(am::AllostericMechanism)
    out = Parameter[]
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        rep = _group_rep(group, fes)
        st = cat_allo_state(am, g) === :EqualAI ? :EqualAI : :A
        append!(out, _emit_cat_params_for_rep(rep, st))
    end
    append!(out, _all_i_state_parameters(am))
    for site in regulatory_sites(am)
        for (lig, tag) in zip(ligands(site), allo_states(site))
            tag === :OnlyI && continue
            push!(out, Kreg(site, lig, :A))
        end
    end
    push!(out, Lallo())
    out
end
```

- [ ] **Step 5: Refactor `_ss_rate_constant_names` allosteric branch (L621-640)**

Read `src/rate_eq_derivation.jl` lines 621-640.

Replace the body:

```julia
function _ss_rate_constant_names(em::AllostericEnzymeMechanism)
    am = AllostericMechanism(em)
    a_names = Set{Symbol}()
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        rep = _group_rep(group, fes)
        is_equilibrium(rep) && continue
        st = cat_allo_state(am, g) === :EqualAI ? :EqualAI : :A
        for p in _emit_cat_params_for_rep(rep, st)
            push!(a_names, name(p, am))
        end
    end
    i_names = Set{Symbol}(name(p, am) for p in _all_i_state_parameters(am)
                          if p isa Union{Kon, Koff, Kfor, Krev})
    union(a_names, i_names)
end
```

- [ ] **Step 6: Run tests — expect 3 helper-blocking tests to FAIL**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: failures in `test_rate_eq_derivation.jl` at lines 1536, 1545, 1548, 1555, 1566, 1576, 1583, 1594, 1599-1604, 1612-1617 (per Q-011). These directly assert on `_onlyA_parameters`, `_I_rename_parameters`, `_all_i_state_parameters` outputs.

The refactor in this task changes nothing about the helpers' return values (they emit the same Parameters in the same order — just via the shared helper now). So tests SHOULD still pass.

**If tests pass:** great, no Task 2.2 work needed.

**If tests fail:** the failure is in test fixture data, not helper behavior — proceed to Task 2.2 which is for the case where F-043's refactor exposes a fixture mismatch. (Most likely: tests pass, since the helper preserves outputs byte-for-byte.)

- [ ] **Step 7: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
refactor: extract _emit_cat_params_for_rep shared helper

Finding F-043 (smaller-helper variant, per Q-015) from the 2026-05-30
audit. The 4-way is_equilibrium × is_binding switch over
Kd/Kiso/Kon+Koff/Kfor+Krev was duplicated in 5 places. Extract as
_emit_cat_params_for_rep(rep, state) → Vector{Parameter}.

Three call sites refactored to use it: _onlyA_parameters,
_all_i_state_parameters, _enumerate_parameters_full_allosteric. The
_ss_rate_constant_names allosteric branch also uses it. The other two
(_I_rename_parameters, _A_rename_parameters) emit pairs, not
Parameters, so they stay separate.

~20 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 2.2: Update test-blocked helper tests if F-043 exposed fixture mismatches

**Files (conditional):**
- Modify: `test/test_rate_eq_derivation.jl:1536-1617` (only if Task 2.1 Step 6 failed)
- Append: `docs/superpowers/refactor-deleted-tests.md` (log replacements per CLAUDE.md)

Skip this task if Task 2.1's tests all passed.

- [ ] **Step 1: Identify which assertions broke**

If tests failed, re-run with `--verbose` to see exact failure locations:

```bash
julia --project test/test_rate_eq_derivation.jl
```

Inspect the failure output for the lines mentioned in Q-011 (1536, 1545, 1548, 1555, 1566, 1576, 1583, 1594, 1599-1604, 1612-1617).

- [ ] **Step 2: For each broken assertion, replace with behavior-level test**

The principle: instead of asserting on `_onlyA_parameters(am)` directly, assert on the observable: `fitted_params(am)` should contain (or not contain) certain Symbols; `rate_equation_string(am)` should match expected text.

For each broken line, follow this pattern:

```julia
# OLD: direct helper assertion
@test name.(_onlyA_parameters(am), Ref(am)) == [:K_A_S_E, :k_A_ES_to_E]

# NEW: behavior-level assertion via public surface
@test :K_A_S_E in fitted_params(am)
@test :k_A_ES_to_E in fitted_params(am)
# OR if more specific structure needed:
@test :K_I_S_E ∉ fitted_params(am)  # OnlyA: no I-state mirror
```

Apply this pattern to each broken assertion. **Do not delete tests — replace them.**

- [ ] **Step 3: Log replacements in refactor-deleted-tests.md**

Read `docs/superpowers/refactor-deleted-tests.md`. Append a section:

```
## 2026-05-30 — F-043 _emit_cat_params_for_rep refactor

Three helpers (`_onlyA_parameters`, `_I_rename_parameters`, `_all_i_state_parameters`)
had direct test assertions in `test/test_rate_eq_derivation.jl:1536-1617`. After
the F-043 refactor (extracting `_emit_cat_params_for_rep`) these helpers retain
the same return values, but for any assertion that exposed a fixture mismatch
during refactor, the test was replaced with a behavior-level assertion via
`fitted_params(am)` / `rate_equation_string(am)`.

Replacement coverage:
- (list any specific tests replaced, with line:reason)
```

- [ ] **Step 4: Re-run the suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 5: Commit**

```bash
git add test/test_rate_eq_derivation.jl docs/superpowers/refactor-deleted-tests.md
git commit -m "$(cat <<'EOF'
test: replace helper-direct assertions with behavior assertions for F-043

Per Q-011 from the 2026-05-30 audit, three private helpers
(_onlyA_parameters, _I_rename_parameters, _all_i_state_parameters)
had direct test assertions that constrained their signatures. After
the F-043 refactor (Task 2.1 extracting _emit_cat_params_for_rep),
replace with behavior-level assertions via fitted_params(am) and
rate_equation_string(am). Replacements logged in
refactor-deleted-tests.md.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 2.3: Rewrite `_step_sides` to take Step directly

**Files:**
- Modify: `src/rate_eq_derivation.jl:155-167`

Finding F-040. Sets up F-039.

- [ ] **Step 1: Replace `_step_sides` with a Step-taking version**

Read `src/rate_eq_derivation.jl` lines 144-167.

Replace the body (L155-167) and docstring:

```julia
"""
Per-step side breakdown for the rate equation derivation. Returns
`(from_species_sym, to_species_sym, m_lhs_syms, m_rhs_syms)` for a
single `Step`. Direction is read directly from the canonical-form
Step (binding steps have metabolite on `to_species` per the Step
constructor invariant; iso steps are canonicalized physical-forward
in the Mechanism constructor). The metabolite-on-which-side
projection mirrors what `_step_tuple_from_sig` used to emit at
@generated time, but reads from Step fields directly.
"""
function _step_sides(s::Step)
    e_lhs = name(from_species(s))
    e_rhs = name(to_species(s))
    if is_binding(s)
        bm = bound_metabolite(s)
        # Canonical: bound metabolite resides in to_species.bound; emit it
        # on the lhs syms list to match the old rxns-tuple convention
        # consumed by _compute_alpha / _compute_numerator (m_lhs is the
        # metabolite-on-the-binding-step side).
        return e_lhs, e_rhs, Symbol[name(bm)], Symbol[]
    end
    e_lhs, e_rhs, Symbol[], Symbol[]
end
```

Note: the `enz_set` parameter (used in the old form to split a side that already had Symbol tuples) is no longer needed — Step's `from_species`/`to_species` are already enzyme Species, and `bound_metabolite` is already the metabolite. No splitting needed.

- [ ] **Step 2: Find all call sites and update signatures**

```bash
grep -n "_step_sides(" src/rate_eq_derivation.jl
```

Expected output: 4 call sites (around L238, 333, 439, 467 per the catalog).

For each call site, the old signature was `_step_sides(rxns, idx, enz_set)`. The new signature is `_step_sides(s::Step)`. The caller's loop must now iterate over Steps, not rxns indexes.

Defer the actual call-site updates to Task 2.5 (where `_raw_symbolic_rate_polys` is rewritten holistically).

- [ ] **Step 3: Add a temporary deprecation shim so the build doesn't break before Task 2.5 lands**

Below the new `_step_sides(s::Step)` definition, add:

```julia
# TEMPORARY shim for the rxns-based callers; deleted in Task 2.5 when
# _raw_symbolic_rate_polys / _compute_alpha / _compute_numerator switch
# to flat-Step iteration.
function _step_sides(rxns, src_idx::Int, enz_set)
    lhs, rhs, _, _ = rxns[src_idx]
    _, m_lhs = _split_reaction_side(lhs, enz_set)
    _, m_rhs = _split_reaction_side(rhs, enz_set)
    e_lhs = first(s for s in lhs if s in enz_set)
    e_rhs = first(s for s in rhs if s in enz_set)
    e_lhs, e_rhs, m_lhs, m_rhs
end
```

This preserves the old behavior for now; the shim deletes in Task 2.5.

- [ ] **Step 4: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green. The new `_step_sides(::Step)` has no callers yet; the shim covers the existing rxns-based callers.

- [ ] **Step 5: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
refactor: add Step-taking _step_sides; keep rxns shim until Task 2.5

Finding F-040 from the 2026-05-30 audit. Introduces _step_sides(s::Step)
that reads directly from Step fields instead of opaque rxns tuples.
A temporary shim preserves the rxns-based signature for existing
callers until Task 2.5 migrates them.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 2.4: Rewrite `_free_enz_set` to walk Mechanism directly

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:47-67`

Finding F-049.

- [ ] **Step 1: Replace the body**

Read `src/thermodynamic_constr_for_rate_eq_derivation.jl` lines 46-70.

Replace the function:

```julia
"""
Set of enzyme-form names that are NOT the RHS of any canonical RE
binding step `F + met… ⇌ F_bound`. Walks `Mechanism.steps` directly:
for each RE binding step, the canonical form puts the bound metabolite
on `to_species`, so `to_species`'s name is excluded from the free set.
Iso steps don't determine binding state. SS steps' direction is not
canonicalized so they don't participate.

Shared by the kinetic-group name representative and the Haldane
elimination pivot.
"""
function _free_enz_set(m::Union{Mechanism, AllostericMechanism})
    enz_names = Set{Symbol}()
    for group in steps(m), s in group
        push!(enz_names, name(from_species(s)))
        push!(enz_names, name(to_species(s)))
    end
    free_enz_set = copy(enz_names)
    for group in steps(m), s in group
        is_equilibrium(s) || continue
        is_binding(s) || continue
        # Canonical: bound metabolite resides on to_species. The from-side
        # is the "free + met" reactant; the to-side is the bound form.
        delete!(free_enz_set, name(to_species(s)))
    end
    free_enz_set
end
```

- [ ] **Step 2: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green. Semantically identical to the opaque-tuple version: same forms in, same forms out.

- [ ] **Step 3: Commit**

```bash
git add src/thermodynamic_constr_for_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
refactor: _free_enz_set walks Mechanism.steps directly

Finding F-049 (Cluster B) from the 2026-05-30 audit. Removes the
em = compile_mechanism(m) lift + reactions(em) opaque-tuple walk.
Same semantics: canonical RE binding puts the bound metabolite on
to_species, so to_species is excluded from the free set.

~5 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 2.5: Rewrite `_raw_symbolic_rate_polys`, `_compute_alpha`, `_compute_numerator` to walk `Mechanism.steps` directly

**Files:**
- Modify: `src/rate_eq_derivation.jl:215-502` (the three core derivation functions)

Finding F-039. The big Cluster B move. Tests are the regression guard.

- [ ] **Step 1: Rewrite `_compute_alpha` (L224-287)**

Read `src/rate_eq_derivation.jl` lines 215-287.

Replace the function. The key change: iterate `_flat_steps(mech)` instead of `enumerate(rxns)`. `eq_steps[idx]` becomes `is_equilibrium(s)`. The rxns-derived `(e_lhs, e_rhs, m_lhs, m_rhs)` become Step-derived via `_step_sides(s::Step)` (or inline):

```julia
"""
Compute alpha factors (relative concentrations within RE groups) as POLY
values. Iterates `mech.steps` directly. Binding steps' direction comes
from the canonical Step form (metabolite-on-`to_species`); iso steps
are physical-forward (canonicalized in the Mechanism constructor).
`step_to_K[idx]` is the parameter Symbol for the RE step at flat
position `idx` (rep-renamed via the `name(p, m)` chokepoint).
"""
function _compute_alpha(mech::Mechanism, enz_species, enz_set,
                        enz_name_to_form, groups, step_to_K)
    N = length(enz_species)
    alpha_num = Vector{POLY}(fill(poly_one(), N))
    alpha_den = Vector{POLY}(fill(poly_one(), N))
    flat = _flat_steps(mech)

    for group in groups
        length(group) == 1 && continue
        visited = Set{Int}([group[1]])
        queue = [group[1]]
        while !isempty(queue)
            cur = popfirst!(queue)
            for (idx, (s, _)) in enumerate(flat)
                is_equilibrium(s) || continue
                e_l, e_r, m_l, m_r = _step_sides(s)
                i_f = enz_name_to_form[e_l]
                j_f = enz_name_to_form[e_r]
                K = poly_sym(step_to_K[idx])
                is_iso = isempty(m_l) && isempty(m_r)
                if i_f == cur && j_f ∉ visited
                    if is_iso
                        alpha_num[j_f] = poly_mul(alpha_num[cur], K)
                        alpha_den[j_f] = alpha_den[cur]
                    else
                        alpha_num[j_f] = poly_mul(alpha_num[cur], poly_sym(m_l[1]))
                        alpha_den[j_f] = poly_mul(alpha_den[cur], K)
                        isempty(m_r) ||
                            (alpha_den[j_f] = poly_mul(alpha_den[j_f], poly_sym(m_r[1])))
                    end
                    push!(visited, j_f); push!(queue, j_f)
                elseif j_f == cur && i_f ∉ visited
                    if is_iso
                        alpha_num[i_f] = alpha_num[cur]
                        alpha_den[i_f] = poly_mul(alpha_den[cur], K)
                    else
                        alpha_num[i_f] = poly_mul(alpha_num[cur], K)
                        alpha_den[i_f] = poly_mul(alpha_den[cur], poly_sym(m_l[1]))
                    end
                    push!(visited, i_f); push!(queue, i_f)
                end
            end
        end
    end

    # Compute sigma per group (unchanged from original)
    sigma_num = Vector{POLY}(undef, length(groups))
    sigma_den = Vector{POLY}(undef, length(groups))
    for (g, group) in enumerate(groups)
        if length(group) == 1
            sigma_num[g] = sigma_den[g] = poly_one()
        else
            sigma_den[g] = reduce(poly_mul, alpha_den[i] for i in group)
            sigma_num[g] = reduce(poly_add,
                poly_mul(alpha_num[i],
                    reduce(poly_mul, (alpha_den[j] for j in group if j != i);
                        init=poly_one()))
                for i in group)
        end
    end
    alpha_num, alpha_den, sigma_num, sigma_den
end
```

- [ ] **Step 2: Rewrite `_raw_symbolic_rate_polys(mech::Mechanism, ...)` (L309-408)**

Read `src/rate_eq_derivation.jl` lines 300-408.

Replace the function signature and body:

```julia
"""
Build raw numerator and denominator POLYs for the rate equation by
walking the lifted `Mechanism`. Parameter Symbols on the leaves of
`num`/`den` are produced via the `name(p, mech)` chokepoint (which
collapses kinetic-group members to their rep's name). `rename_map`
then applies any single-symbol Wegscheider ties as a post-pass.
"""
function _raw_symbolic_rate_polys(mech::Mechanism, step_params, rename_map,
                                  subs_species, prods_species)
    enz_species, groups, form_to_group = _compute_re_groups(mech)
    enz_set = Set(name(es) for es in enz_species)
    enz_name_to_form = Dict{Symbol, Int}(
        name(es) => i for (i, es) in enumerate(enz_species))
    flat = _flat_steps(mech)
    step_to_K = Dict{Int, Symbol}(
        i => name(step_params[i][1], mech)
        for i in eachindex(flat) if is_equilibrium(flat[i][1]))
    alpha_num, alpha_den, sigma_num, sigma_den =
        _compute_alpha(mech, enz_species, enz_set,
                       enz_name_to_form, groups, step_to_K)
    G = length(groups)

    R = [poly_zero() for _ in 1:G, _ in 1:G]
    for (idx, (s, _)) in enumerate(flat)
        is_equilibrium(s) && continue
        e_lhs, e_rhs, m_lhs, m_rhs = _step_sides(s)
        i_form = enz_name_to_form[e_lhs]
        j_form = enz_name_to_form[e_rhs]
        g1, g2 = form_to_group[i_form], form_to_group[j_form]
        kf_poly = poly_sym(name(step_params[idx][1], mech))
        kr_poly = poly_sym(name(step_params[idx][2], mech))
        R[g1, g2] = poly_add(R[g1, g2],
            _ss_contrib(kf_poly, m_lhs, i_form, alpha_num, alpha_den, groups[g1]))
        R[g2, g1] = poly_add(R[g2, g1],
            _ss_contrib(kr_poly, m_rhs, j_form, alpha_num, alpha_den, groups[g2]))
    end

    L = [i == j ? poly_zero() : poly_neg(R[i,j])
         for i in 1:G, j in 1:G]
    for i in 1:G
        L[i, i] = reduce(poly_add, R[i, j] for j in 1:G if j != i; init=poly_zero())
    end
    D = [begin
        idx = [r for r in 1:G if r != root]
        isempty(idx) ? poly_one() : sym_det(L[idx, idx], G - 1)
    end for root in 1:G]

    normalize = G == 1 && sigma_den[1] != poly_one()
    den = poly_zero()
    for g in 1:G
        raw_sigma = if normalize
            reduce(poly_add, (_poly_div_mono(alpha_num[i], alpha_den[i]) for i in groups[g]))
        else
            sigma_num[g]
        end
        csigma = _rename_symbols(raw_sigma, rename_map)
        den = poly_add(den, poly_mul(csigma, D[g]))
    end

    num, nu_ref = _compute_numerator(
        mech, enz_set, enz_name_to_form, step_params,
        alpha_num, alpha_den, form_to_group, groups,
        D, subs_species, prods_species)
    normalize && (num = _poly_div_mono(num, sigma_den[1]))

    abs_nu = abs(nu_ref)
    if abs_nu != 1
        den = poly_mul(den, poly_const(abs_nu))
    end

    num = _rename_symbols(num, rename_map)
    den = _rename_symbols(den, rename_map)
    num, den
end

function _raw_symbolic_rate_polys(M::Type{<:EnzymeMechanism})
    mech = Mechanism(M())
    step_params = _step_parameters(mech)
    rename_map = _build_wegscheider_rename_map(M)
    _raw_symbolic_rate_polys(mech, step_params, rename_map,
                              substrates(mech.reaction),
                              products(mech.reaction))
end
```

Note: `substrates(M())` / `products(M())` becomes `substrates(mech.reaction)` / `products(mech.reaction)` since `mech.reaction::EnzymeReaction` already has the data. This is an early instance of Wave 3's accessor demotion happening locally.

- [ ] **Step 3: Rewrite `_compute_numerator` (L427-502)**

Read `src/rate_eq_derivation.jl` lines 423-502.

Replace the body. The signature loses `rxns, eq_steps`; gains `mech` (via the closure already in scope, but we'll be explicit). Add `flat = _flat_steps(mech)` at top:

```julia
"""
Compute the numerator polynomial by selecting an appropriate metabolite
to track through SS steps. Returns `(num::POLY, nu_ref::Int)`.
"""
function _compute_numerator(
    mech::Mechanism, enz_set, enz_name_to_form, step_params,
    alpha_num, alpha_den, form_to_group, groups,
    D, subs_species, prods_species,
)
    ref_name = subs_species[1]
    nu_ref = (count(==(ref_name), prods_species) -
              count(==(ref_name), subs_species))

    flat = _flat_steps(mech)

    ss_mets, re_mets = Set{Symbol}(), Set{Symbol}()
    for (s, _) in flat
        _, _, m_lhs, m_rhs = _step_sides(s)
        target = is_equilibrium(s) ? re_mets : ss_mets
        for met in m_lhs; push!(target, met); end
        for met in m_rhs; push!(target, met); end
    end

    met_name = nothing
    nu_met = nu_ref
    if ref_name ∈ ss_mets && ref_name ∉ re_mets
        met_name = ref_name
    elseif !isempty(ss_mets)
        all_mets = Dict{Symbol, Int}()
        for n in subs_species; all_mets[n] = get(all_mets, n, 0) - 1; end
        for n in prods_species; all_mets[n] = get(all_mets, n, 0) + 1; end
        ss_only = setdiff(ss_mets, re_mets)
        search = isempty(ss_only) ? ss_mets : ss_only
        met_name = something(
            iterate(m for m in search if get(all_mets, m, 0) != 0),
            (first(ss_mets),),
        )[1]
        nu_met = get(all_mets, met_name, 0)
    end

    result = poly_zero()
    for (idx, (s, _)) in enumerate(flat)
        is_equilibrium(s) && continue
        e_lhs, e_rhs, m_lhs, m_rhs = _step_sides(s)
        i_form = enz_name_to_form[e_lhs]
        j_form = enz_name_to_form[e_rhs]
        g1, g2 = form_to_group[i_form], form_to_group[j_form]
        rf = _ss_contrib(
            poly_sym(name(step_params[idx][1], mech)), m_lhs, i_form,
            alpha_num, alpha_den, groups[g1])
        rr = _ss_contrib(
            poly_sym(name(step_params[idx][2], mech)), m_rhs, j_form,
            alpha_num, alpha_den, groups[g2])
        flux = poly_sub(poly_mul(rf, D[g1]), poly_mul(rr, D[g2]))
        if met_name === nothing
            result = poly_add(result, flux)
        else
            met_f = isempty(m_lhs) ? nothing : first(m_lhs)
            met_r = isempty(m_rhs) ? nothing : first(m_rhs)
            (met_f !== met_name && met_r !== met_name) && continue
            result = met_f === met_name ?
                poly_add(result, flux) : poly_sub(result, flux)
        end
    end

    if nu_met != 0 && nu_met != nu_ref
        ratio = nu_ref // nu_met
        if ratio == -1
            result = poly_neg(result)
        elseif ratio != 1
            error("Non-unit stoichiometric ratio not supported")
        end
    end

    result, nu_ref
end
```

- [ ] **Step 4: Delete the `_step_sides(rxns, idx, enz_set)` shim from Task 2.3**

The shim was added in Task 2.3 Step 3. With the rewrites in this task, all callers now use `_step_sides(s::Step)`. Delete the shim.

- [ ] **Step 5: grep for any remaining rxns-based callers**

```bash
grep -n "rxns\[" src/rate_eq_derivation.jl
grep -n "for (idx, _) in enumerate(rxns)" src/rate_eq_derivation.jl
grep -n "_step_sides(rxns" src/rate_eq_derivation.jl
```

Expected output: empty. If any remain, fix them.

- [ ] **Step 6: Run the rate-equation Expr-shape regression tests**

```bash
julia --project test/test_rate_eq_derivation.jl
```

Expected: green. The Expr-shape and flat-string regression tests in this file are byte-stability guards for the rate-equation output. If they pass, the refactor preserved emission byte-for-byte.

- [ ] **Step 7: Run the rate_equation performance test specifically**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -E "(allocations|performance|<100|rate_equation)"
```

The performance gate (allocations == 0 AND time < 100ns) must stay green for every spec in MECHANISM_TEST_SPECS.

- [ ] **Step 8: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 9: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
refactor: derivation walks Mechanism.steps directly (drop rxns shuffle)

Finding F-039 (Cluster B) from the 2026-05-30 audit. Per Q-016 (and
the long-standing source comment at L325-329 explicitly inviting this
cleanup): _raw_symbolic_rate_polys, _compute_alpha, and
_compute_numerator now iterate _flat_steps(mech) directly instead of
constructing a rxns opaque-tuple intermediate via reactions(m). All
direction info reads from Step.from_species / Step.to_species /
Step.is_equilibrium / Step.bound_metabolite, with the canonical
metabolite-on-to_species invariant.

The rate_equation Expr-shape regression tests and the 0-alloc /
<100ns performance gate confirm byte-stability and unchanged runtime
characteristics.

~30 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 2.6: Rewrite `_thermodynamic_constraints` to take Mechanism directly

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:138-207`

Finding F-050.

- [ ] **Step 1: Rewrite the function to take Mechanism**

Read `src/thermodynamic_constr_for_rate_eq_derivation.jl` lines 138-207.

Replace the entire function. Key change: take `mech::Mechanism` directly, walk `mech.steps` and `mech.reaction`. Drop all accessor calls on `EnzymeMechanism{Sig}`.

```julia
function _thermodynamic_constraints(mech::Mechanism)
    flat = _flat_steps(mech)
    enz_names = collect(_enumerate_species_names(mech))
    enz_name_to_idx = Dict(n => i for (i, n) in enumerate(enz_names))
    nsteps = length(flat)
    met_names = Symbol[name(metabolite(ra)) for ra in reactants(mech.reaction)]
    subs_species = Symbol[name(s) for s in substrates(mech.reaction)]
    prods_species = Symbol[name(p) for p in products(mech.reaction)]

    # Enzyme incidence matrix
    B = zeros(Int, length(enz_names), nsteps)
    for (j, (s, _)) in enumerate(flat)
        i_from = enz_name_to_idx[name(from_species(s))]
        i_to   = enz_name_to_idx[name(to_species(s))]
        B[i_from, j] -= 1
        B[i_to,   j] += 1
    end

    # Stoichiometry matrix (rows = metabolites, cols = steps)
    met_idx = Dict(n => i for (i, n) in enumerate(met_names))
    stoich_mat = zeros(Int, length(met_names), nsteps)
    for (j, (s, _)) in enumerate(flat)
        if is_binding(s)
            stoich_mat[met_idx[name(bound_metabolite(s))], j] -= 1
        end
        # Iso steps that produce / consume products (e.g. catalytic-release
        # iso): these surface via Species' bound metabolites differing
        # between from and to. Walk both bound lists.
        from_bound = Set(name(m) for m in bound(from_species(s)))
        to_bound = Set(name(m) for m in bound(to_species(s)))
        for met in setdiff(from_bound, to_bound)
            haskey(met_idx, met) && (stoich_mat[met_idx[met], j] -= 1)
        end
        for met in setdiff(to_bound, from_bound)
            haskey(met_idx, met) && (stoich_mat[met_idx[met], j] += 1)
        end
    end
    # Note: the above incorporates both binding-step met consumption AND
    # iso-step net-met change. This subsumes what stoich_matrix(m) /
    # metabolite_row_range(m) used to provide.

    NS = _integer_nullspace(B)
    nc = size(NS, 2)
    nc == 0 && return zeros(Int, 0, size(B, 2)), Int[]

    nu_net = zeros(Int, length(met_names))
    for nm in subs_species
        nu_net[met_idx[nm]] -= 1
    end
    for nm in prods_species
        nu_net[met_idx[nm]] += 1
    end

    function classify_cycle(nu_cycle, i)
        all(nu_cycle .== 0) && return 0
        c = nothing
        for j in eachindex(nu_cycle)
            if nu_net[j] == 0
                nu_cycle[j] != 0 && error(
                    "Cycle $i produces metabolite change not proportional to net reaction")
            else
                c_j = nu_cycle[j] // nu_net[j]
                if c === nothing
                    c = c_j
                elseif c_j != c
                    error("Cycle $i produces metabolite change not proportional to net reaction")
                end
            end
        end
        c !== nothing && denominator(c) == 1 ? Int(c) :
            error("Cycle $i produces metabolite change not proportional to net reaction")
    end

    C = NS'
    rhs_coeffs = [classify_cycle(stoich_mat * C[i, :], i) for i in 1:nc]
    return C, rhs_coeffs
end

# Type-dispatching wrapper kept until Cluster A lands fully — preserves
# the existing @generated callers in _dependent_param_exprs_kernel.
_thermodynamic_constraints(M::Type{<:EnzymeMechanism}) =
    _thermodynamic_constraints(Mechanism(M()))

# Walk Mechanism.steps; emit distinct enzyme-form Symbol names in
# step-walk order. Used by _thermodynamic_constraints and friends.
function _enumerate_species_names(mech::Mechanism)
    seen = Symbol[]
    for group in steps(mech), s in group
        for sp in (from_species(s), to_species(s))
            nm = name(sp)
            nm in seen || push!(seen, nm)
        end
    end
    seen
end
```

- [ ] **Step 2: Run the thermodynamic-constraint tests**

```bash
julia --project test/test_rate_eq_derivation.jl
```

The cycle-classification logic is byte-stable in the rewrite. Test expected to be green.

- [ ] **Step 3: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 4: Commit**

```bash
git add src/thermodynamic_constr_for_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
refactor: _thermodynamic_constraints walks Mechanism directly

Finding F-050 (Cluster B) from the 2026-05-30 audit. Removes the
m = M() + enzyme_forms(m) + reactions(m) + stoich_matrix(m) +
metabolite_row_range(m) accessor chain. Builds enzyme incidence and
stoichiometry matrices directly from Step.from_species /
Step.to_species / Step.bound_metabolite. A type-dispatching wrapper
preserves the existing _dependent_param_exprs_kernel callers until
Cluster A lands fully.

~10 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 2.7: Rewrite `_dependent_param_exprs_kernel` to walk Mechanism directly

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:253-379`

Finding F-052.

- [ ] **Step 1: Rewrite the kernel**

Read `src/thermodynamic_constr_for_rate_eq_derivation.jl` lines 253-379.

Replace the function. Key changes: takes `mech::Mechanism` and `rename::Dict`. Drops `m = M()` and `mech = Mechanism(m)` double-lift. Drops `rxns = reactions(m)` and walks `_flat_steps(mech)` directly. binding_K_set is built by reading Step fields.

```julia
function _dependent_param_exprs_kernel(
    mech::Mechanism,
    rename::AbstractDict{Symbol, Symbol},
)
    flat = _flat_steps(mech)
    free_enz_set = _free_enz_set(mech)
    C, rhs_coeffs = _thermodynamic_constraints(mech)
    all_params = _raw_param_symbols(mech)
    nc = size(C, 1)
    nsteps = size(C, 2)
    nc == 0 && return (Dict{Symbol, Union{Symbol, Expr}}(), Tuple(all_params))

    sym_col = Dict(p => i for (i, p) in enumerate(all_params))
    n_vars = length(all_params)
    step_params = _step_parameters(mech)
    step_name(p::Parameter) = get(rename, name(p, mech), name(p, mech))

    # binding_K_set built from Step fields (RE binding step with the bound
    # metabolite on the canonical to-side ⇒ the from-side is the "free + met"
    # side, no enzyme on from-side means... actually canonical Step puts the
    # enzyme on both sides. Use the chokepoint rep-name.)
    binding_K_set = Set{Symbol}()
    for (j, (s, _)) in enumerate(flat)
        is_equilibrium(s) && is_binding(s) || continue
        push!(binding_K_set, step_name(step_params[j][1]))
    end

    A = zeros(Rational{BigInt}, nc, n_vars)
    rhs = Rational{BigInt}.(rhs_coeffs)
    for i in 1:nc, j in 1:nsteps
        C[i, j] == 0 && continue
        if is_equilibrium(flat[j][1])
            sym = step_name(step_params[j][1])
            sign_factor = sym in binding_K_set ? -1 : 1
            A[i, sym_col[sym]] += sign_factor * C[i, j]
        else
            kf = step_name(step_params[j][1])
            kr = step_name(step_params[j][2])
            A[i, sym_col[kf]] += C[i, j]
            A[i, sym_col[kr]] -= C[i, j]
        end
    end

    # Pivot priority (unchanged from original)
    priority = zeros(Int, n_vars)
    for j in 1:nsteps
        step = step_params[j][1].step
        base = _step_priority(step, free_enz_set)
        if is_equilibrium(flat[j][1])
            s_sym = step_name(step_params[j][1])
            haskey(sym_col, s_sym) && (priority[sym_col[s_sym]] = base)
        else
            for (offset, p) in enumerate(step_params[j])
                s_sym = step_name(p)
                haskey(sym_col, s_sym) &&
                    (priority[sym_col[s_sym]] = base + offset - 1)
            end
        end
    end

    # Gaussian elimination (unchanged from original)
    pivot_entries = Tuple{Int, Int}[]
    pivot_col_set = Set{Int}()
    wA, wrhs = copy(A), copy(rhs)
    for i in 1:nc
        best_col, best_pri = 0, -1
        for c in 1:n_vars
            c in pivot_col_set && continue
            wA[i, c] == 0 && continue
            if priority[c] > best_pri
                best_pri = priority[c]
                best_col = c
            end
        end
        if best_col == 0
            wrhs[i] == 0 && continue
            error("Thermodynamically contradictory mechanism: " *
                  "constraint row $i reduces to 0 = $(wrhs[i]) * log(Keq)")
        end
        push!(pivot_entries, (i, best_col))
        push!(pivot_col_set, best_col)
        pv = wA[i, best_col]
        wA[i, :] ./= pv
        wrhs[i] /= pv
        for r in 1:nc
            if r != i && wA[r, best_col] != 0
                f = wA[r, best_col]
                wA[r, :] .-= f .* wA[i, :]
                wrhs[r] -= f * wrhs[i]
            end
        end
    end

    dep_exprs = Dict{Symbol, Union{Symbol, Expr}}()
    for (prow, pcol) in pivot_entries
        factors = [
            (all_params[c], -wA[prow, c])
            for c in 1:n_vars
            if c != pcol && wA[prow, c] != 0
        ]
        dep_exprs[all_params[pcol]] = build_power_expr(wrhs[prow], factors)
    end
    dep_set = Set(keys(dep_exprs))
    return dep_exprs, Tuple(p for p in all_params if p ∉ dep_set)
end

# Type-dispatching wrapper preserves the existing call sites in
# _dependent_param_exprs and _build_kinetic_rename_map / _build_wegscheider_rename_map.
_dependent_param_exprs_kernel(M::Type{<:EnzymeMechanism},
                              rename::AbstractDict{Symbol, Symbol}) =
    _dependent_param_exprs_kernel(Mechanism(M()), rename)
```

- [ ] **Step 2: Update `_dependent_param_exprs(M::Type{<:EnzymeMechanism})` at L222-239 to use the wrapper**

The function at L222-239 already calls `_dependent_param_exprs_kernel(M, rename)`. With the wrapper added, no change needed.

- [ ] **Step 3: Run dep-set invariance tests as regression guard**

```bash
julia --project test/test_dep_set_invariance.jl
```

Expected: green. The dep-expression set is byte-stable in the rewrite.

- [ ] **Step 4: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 5: Commit**

```bash
git add src/thermodynamic_constr_for_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
refactor: _dependent_param_exprs_kernel walks Mechanism directly

Finding F-052 (Cluster B) from the 2026-05-30 audit. Removes the
m = M() + mech = Mechanism(m) double-lift and the rxns opaque-tuple
walk. binding_K_set built from Step.is_equilibrium / Step.is_binding
directly. A type-dispatching wrapper preserves the existing call
sites at thermodynamic_constr.jl:222 and rate_eq_derivation.jl:416
until Cluster A demotes the singletons.

~15 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 2.8: Wave 2 verification

- [ ] **Step 1: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 2: Run rate_equation perf gate explicitly**

```bash
julia --project test/test_rate_eq_derivation.jl 2>&1 | grep -E "(allocations|<100|test_rate_equation_performance)"
```

Expected: 0 allocations / <100ns reported for every spec. If anything regressed, find the cause before proceeding — Wave 3 builds on a known-good back-end.

- [ ] **Step 3: Measure LOC delta**

Run the same baseline-counter command from Wave 1 Task 1.9 Step 2.

Expected: ~5,515 non-comment non-doc LOC (down ~120 from Wave 1 start).

# Wave 3 — Singleton-type demotion + expansion-move consolidation

Cluster A + Cluster D. The big architectural move. ~620 LOC net delete + ~150 LOC of in-place rewrites.

**Pre-conditions (verify before starting Wave 3):**

- [ ] Wave 1 and Wave 2 commits all landed and tested green
- [ ] `rate_equation` perf gate green (0 alloc, <100ns)
- [ ] `test_accessors.jl` perf gate noted as negotiable — will be re-baselined or deleted in this wave
- [ ] All test files' `compile_mechanism(m)` callers are identified (~50 callers per Q-001)

## Task 3.1: Add `_with_steps` / `_with_cat_allo_states` / `_with_reg_sites` update constructors

**Files:**
- Modify: `src/types.jl` (add helpers near the Mechanism / AllostericMechanism struct definitions, around L405 and L498)

Finding F-056 setup.

- [ ] **Step 1: Add helpers to types.jl after the Mechanism struct (around L405)**

```julia
"""
Return a new Mechanism with `new_steps` but the same reaction. Used by
expansion moves to swap step structure while preserving the rest of
the mechanism's shape.
"""
_with_steps(m::Mechanism, new_steps::Vector{Vector{Step}}) =
    Mechanism(reaction(m), new_steps)
```

- [ ] **Step 2: Add helpers after the AllostericMechanism struct (around L498)**

```julia
"""
Return a new AllostericMechanism with `new_steps` but otherwise
identical fields.
"""
_with_steps(am::AllostericMechanism, new_steps::Vector{Vector{Step}}) =
    AllostericMechanism(reaction(am), new_steps,
                        copy(cat_allo_states(am)),
                        catalytic_multiplicity(am),
                        copy(regulatory_sites(am)))

"""
Return a new AllostericMechanism with `new_cat_allo_states` but
otherwise identical fields.
"""
_with_cat_allo_states(am::AllostericMechanism, new_cat_allo_states::Vector{Symbol}) =
    AllostericMechanism(reaction(am), copy(steps(am)),
                        new_cat_allo_states,
                        catalytic_multiplicity(am),
                        copy(regulatory_sites(am)))

"""
Return a new AllostericMechanism with `new_reg_sites` but otherwise
identical fields.
"""
_with_reg_sites(am::AllostericMechanism, new_reg_sites::Vector{RegulatorySite}) =
    AllostericMechanism(reaction(am), copy(steps(am)),
                        copy(cat_allo_states(am)),
                        catalytic_multiplicity(am),
                        new_reg_sites)

"""
Return a new AllostericMechanism with both new_steps AND new_cat_allo_states.
Useful for expansion moves that split a kinetic group (which adds a
state for the new group).
"""
_with_steps_and_cat_states(am::AllostericMechanism,
                            new_steps::Vector{Vector{Step}},
                            new_cat_allo_states::Vector{Symbol}) =
    AllostericMechanism(reaction(am), new_steps,
                        new_cat_allo_states,
                        catalytic_multiplicity(am),
                        copy(regulatory_sites(am)))
```

- [ ] **Step 3: Run tests as sanity check (helpers are unused so far)**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 4: Commit**

```bash
git add src/types.jl
git commit -m "$(cat <<'EOF'
refactor: add _with_steps / _with_cat_allo_states / _with_reg_sites helpers

Finding F-056 setup from the 2026-05-30 audit. Internal update
constructors that return a new Mechanism / AllostericMechanism with
one field replaced. Used by expansion moves (next commits) to collapse
their Mechanism + AllostericMechanism dual dispatch.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 3.2: Collapse `_expand_re_to_ss` + `_expand_split_kinetic_group` dual dispatches

**Files:**
- Modify: `src/mechanism_enumeration.jl:1101-1187`

Findings F-056 (rest).

- [ ] **Step 1: Replace `_expand_re_to_ss` dual methods with one**

Read `src/mechanism_enumeration.jl` lines 1101-1144.

Replace L1101-1123 (the two methods) with one:

```julia
function _expand_re_to_ss(m::Union{Mechanism, AllostericMechanism})
    results = typeof(m)[]
    for g in kinetic_groups(m)
        all(is_equilibrium, steps(m)[g]) || continue
        push!(results, _with_steps(m, _flip_group_to_ss(steps(m), g)))
    end
    results
end
```

`_flip_group_to_ss` (L1130-1144) stays as-is.

- [ ] **Step 2: Replace `_expand_split_kinetic_group` dual methods with one**

Read `src/mechanism_enumeration.jl` lines 1158-1187.

Replace L1158-1187 (the two methods) with one:

```julia
function _expand_split_kinetic_group(m::Mechanism)
    results = Mechanism[]
    for g in kinetic_groups(m)
        length(steps(m)[g]) >= 2 || continue
        for split_idx in 1:length(steps(m)[g])
            push!(results, _with_steps(m, _split_one_step(steps(m), g, split_idx)))
        end
    end
    results
end

function _expand_split_kinetic_group(am::AllostericMechanism)
    results = AllostericMechanism[]
    for g in kinetic_groups(am)
        length(steps(am)[g]) >= 2 || continue
        for split_idx in 1:length(steps(am)[g])
            new_groups = _split_one_step(steps(am), g, split_idx)
            new_states = vcat(cat_allo_states(am), [cat_allo_states(am)[g]])
            push!(results, _with_steps_and_cat_states(am, new_groups, new_states))
        end
    end
    results
end
```

Note: `_expand_split_kinetic_group` for AllostericMechanism extends `cat_allo_states` with the parent group's tag, so it needs `_with_steps_and_cat_states`, not just `_with_steps`. The two methods stay separate but use the new helpers cleanly.

- [ ] **Step 3: Run enumeration tests**

```bash
julia --project test/test_mechanism_enumeration.jl
```

Expected: green. The mechanism counts (bi-bi 77, etc.) are byte-stable.

- [ ] **Step 4: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
refactor: collapse _expand_re_to_ss / _expand_split_kinetic_group via _with_steps

Finding F-056 from the 2026-05-30 audit. _expand_re_to_ss now has one
method on Union{Mechanism, AllostericMechanism} using _with_steps.
_expand_split_kinetic_group keeps two methods (its AllostericMechanism
variant needs _with_steps_and_cat_states), but both bodies are now
shorter and explicit. ~50 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 3.3: Collapse `_expand_add_dead_end_regulator` dual dispatch

**Files:**
- Modify: `src/mechanism_enumeration.jl:1282-1454`

Finding F-057.

- [ ] **Step 1: Inline the `wrap` callback by switching to direct `_with_*` calls**

Read `src/mechanism_enumeration.jl` lines 1282-1454.

The two top-level methods at L1282-1311 each pass a `wrap` callback to `_expand_add_dead_end_regulator_native`. The native kernel at L1320-1454 invokes `wrap` once at L1449.

Replace the structure: drop the `wrap` parameter, make the kernel return `(new_groups, n_groups_before + 1, new_reaction)` tuples, and let each top-level method construct the result via `_with_*`.

Replace L1282-1311 with:

```julia
function _expand_add_dead_end_regulator(
    m::Mechanism, rxn::EnzymeReaction;
    exclude_regs::Set{Symbol}=Set{Symbol}(),
)
    raw = _expand_add_dead_end_regulator_native(
        m, rxn, Set{Symbol}();
        exclude_regs=exclude_regs)
    [Mechanism(new_reaction, new_groups)
     for (new_groups, _, new_reaction) in raw]
end

function _expand_add_dead_end_regulator(
    am::AllostericMechanism, rxn::EnzymeReaction;
    exclude_regs::Set{Symbol}=Set{Symbol}(),
)
    allo_ligands = Set{Symbol}()
    for site in regulatory_sites(am), lig in ligands(site)
        push!(allo_ligands, name(lig))
    end
    raw = _expand_add_dead_end_regulator_native(
        am, rxn, allo_ligands;
        exclude_regs=exclude_regs)
    [AllostericMechanism(new_reaction, new_groups,
                         vcat(cat_allo_states(am), [:EqualAI]),
                         catalytic_multiplicity(am),
                         copy(regulatory_sites(am)))
     for (new_groups, _, new_reaction) in raw]
end
```

And update the kernel at L1320 to return raw tuples instead of calling `wrap`. Replace the last line (L1449-1450):

```
                push!(results, wrap(new_groups, n_groups_before + 1,
                                    new_reaction))
```

With:

```
                push!(results, (new_groups, n_groups_before + 1, new_reaction))
```

And update the kernel signature L1320-1326 to drop `wrap`:

```julia
function _expand_add_dead_end_regulator_native(
    m::Union{Mechanism, AllostericMechanism},
    rxn::EnzymeReaction,
    additional_excluded::Set{Symbol};
    exclude_regs::Set{Symbol},
)
```

And the kernel's result variable type changes from `typeof(m)[]` to:

```julia
    results = Tuple{Vector{Vector{Step}}, Int, EnzymeReaction}[]
```

Update the early-return at L1327: `isempty(regulators(rxn)) && return typeof(m)[]` becomes `isempty(regulators(rxn)) && return Tuple{Vector{Vector{Step}}, Int, EnzymeReaction}[]`. Same for L1348: `isempty(eligible_regs) && return ...`.

- [ ] **Step 2: Run enumeration tests**

```bash
julia --project test/test_mechanism_enumeration.jl
```

Expected: green. Same mechanisms, same counts.

- [ ] **Step 3: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 4: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
refactor: drop wrap callback from _expand_add_dead_end_regulator_native

Finding F-057 from the 2026-05-30 audit. The wrap callback added an
extra indirection that's clearer as a direct list-comprehension over
the kernel's raw output. Per-type result construction lives in the
two top-level methods, not via callback. ~30 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 3.4: Simplify `_expand_change_allo_state` with `_with_reg_sites`

**Files:**
- Modify: `src/mechanism_enumeration.jl:1628-1662`

Finding F-059.

- [ ] **Step 1: Rewrite `_expand_change_allo_state` using `_with_*` helpers**

Read `src/mechanism_enumeration.jl` lines 1628-1662.

Replace the body:

```julia
function _expand_change_allo_state(am::AllostericMechanism)
    results = AllostericMechanism[]

    for g in 1:length(cat_allo_states(am))
        cat_allo_states(am)[g] == :NonequalAI && continue
        new_states = copy(cat_allo_states(am))
        new_states[g] = :NonequalAI
        push!(results, _with_cat_allo_states(am, new_states))
    end

    for (si, site) in enumerate(regulatory_sites(am))
        for (li, _) in enumerate(ligands(site))
            allo_states(site)[li] == :NonequalAI && continue
            new_sites = copy(regulatory_sites(am))
            new_states = copy(allo_states(site))
            new_states[li] = :NonequalAI
            new_sites[si] = RegulatorySite(
                copy(ligands(site)), multiplicity(site), new_states)
            push!(results, _with_reg_sites(am, new_sites))
        end
    end

    results
end
```

- [ ] **Step 2: Run enumeration tests**

```bash
julia --project test/test_mechanism_enumeration.jl
```

Expected: green.

- [ ] **Step 3: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
refactor: _expand_change_allo_state uses _with_cat_allo_states / _with_reg_sites

Finding F-059 from the 2026-05-30 audit. Replaces inline
AllostericMechanism reconstructions with _with_* helper calls.
~10 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 3.5: Convert 14 `@generated` accessors over `EnzymeMechanism{Sig}` to plain functions

**Files:**
- Modify: `src/types.jl:974-1231`

Finding F-007. Per Q-005, ALL 14 accessors are class (b) or (c) — none are on `rate_equation`'s runtime hot path. Safe to demote.

**Critical:** preserve every accessor's return shape exactly (Tuple of Symbols, NamedTuple, etc.) so the @generated body builders that consume these results at compile time produce the same Expr output.

- [ ] **Step 1: Rewrite `substrates`, `products`, `regulators`, `metabolites` over EnzymeMechanism (L996-1051)**

Read `src/types.jl` lines 994-1051.

Replace each `@generated function f(::EnzymeMechanism{Sig}) where {Sig}` with a plain `function f(em::EnzymeMechanism)` that lifts to `Mechanism(em)` and walks `reaction(m).reactants` directly.

Note: the existing `substrates(r::EnzymeReaction)` and `products(r::EnzymeReaction)` at L288-291 already do the work; just delegate.

```julia
function substrates(em::EnzymeMechanism)
    m = Mechanism(em)
    Tuple(name(s) for s in substrates(reaction(m)))
end

function products(em::EnzymeMechanism)
    m = Mechanism(em)
    Tuple(name(p) for p in products(reaction(m)))
end

function regulators(em::EnzymeMechanism)
    m = Mechanism(em)
    Tuple(name(regulator(rm)) for rm in regulators(reaction(m)))
end

function metabolites(em::EnzymeMechanism)
    m = Mechanism(em)
    out = Symbol[]
    seen = Set{Symbol}()
    for s in substrates(reaction(m))
        n = name(s); n in seen || (push!(seen, n); push!(out, n))
    end
    for p in products(reaction(m))
        n = name(p); n in seen || (push!(seen, n); push!(out, n))
    end
    for rm in regulators(reaction(m))
        n = name(regulator(rm)); n in seen || (push!(seen, n); push!(out, n))
    end
    Tuple(out)
end
```

- [ ] **Step 2: Rewrite `reactions`, `equilibrium_steps`, `n_steps` (L1119-1137)**

```julia
function reactions(em::EnzymeMechanism)
    m = Mechanism(em)
    out = Any[]
    for (g, group) in enumerate(steps(m))
        for s in group
            # Reproduce the (lhs, rhs, is_eq, g) shape that the old
            # @generated reactions built via _step_tuple_from_sig.
            # Canonical: bound metabolite on to_species.bound; emit
            # as the "bound side" (lhs in the old convention).
            e_from = name(from_species(s))
            e_to   = name(to_species(s))
            if is_iso(s)
                push!(out, ((e_from,), (e_to,), is_equilibrium(s), g))
            else
                bm_name = name(bound_metabolite(s))
                # Canonical placement: bm on to_species — emit on lhs
                # (matching the canonical-form orientation; the consumer
                # _step_sides reads m_lhs as the binding metabolite side)
                push!(out, ((e_from, bm_name), (e_to,), is_equilibrium(s), g))
            end
        end
    end
    Tuple(out)
end

function equilibrium_steps(em::EnzymeMechanism)
    m = Mechanism(em)
    Tuple(is_equilibrium(s) for group in steps(m) for s in group)
end

function n_steps(em::EnzymeMechanism)
    m = Mechanism(em)
    sum(length, steps(m); init=0)
end
```

- [ ] **Step 3: Rewrite `kinetic_group`, `kinetic_groups`, `steps_in_group` (L1140-1166)**

```julia
function kinetic_group(em::EnzymeMechanism, idx::Int)
    flat = _flat_steps(Mechanism(em))
    1 ≤ idx ≤ length(flat) ||
        error("kinetic_group: step index $idx out of range 1:$(length(flat))")
    flat[idx][2]
end

function kinetic_groups(em::EnzymeMechanism)
    m = Mechanism(em)
    Tuple(1:length(steps(m)))
end

function steps_in_group(em::EnzymeMechanism, g::Int)
    m = Mechanism(em)
    flat = _flat_steps(m)
    Tuple(i for (i, (_, gid)) in enumerate(flat) if gid == g)
end

steps_in_group(em::EnzymeMechanism, ::Val{G}) where {G} = steps_in_group(em, G)
```

- [ ] **Step 4: Rewrite `enzyme_forms`, `n_states` (L1174-1201)**

```julia
function enzyme_forms(em::EnzymeMechanism)
    m = Mechanism(em)
    met_names = Set(metabolites(em))
    seen = Set{Symbol}()
    forms = Symbol[]
    for group in steps(m), s in group
        for sp in (from_species(s), to_species(s))
            nm = name(sp)
            nm ∉ met_names && nm ∉ seen &&
                (push!(seen, nm); push!(forms, nm))
        end
    end
    Tuple(forms)
end

n_states(em::EnzymeMechanism) = length(enzyme_forms(em))
```

- [ ] **Step 5: Rewrite `stoich_matrix`, `enzyme_row_range`, `metabolite_row_range` (L1204-1231)**

```julia
function stoich_matrix(em::EnzymeMechanism)
    rxns = reactions(em)
    species = (enzyme_forms(em)..., metabolites(em)...)
    sp_idx = Dict(s => i for (i, s) in enumerate(species))
    S = zeros(Int, length(species), length(rxns))
    for (j, (lhs, rhs, _, _)) in enumerate(rxns)
        for s in lhs; S[sp_idx[s], j] -= 1; end
        for s in rhs; S[sp_idx[s], j] += 1; end
    end
    S
end

enzyme_row_range(em::EnzymeMechanism) = 1:n_states(em)
metabolite_row_range(em::EnzymeMechanism) =
    (n_states(em) + 1):(n_states(em) + length(metabolites(em)))
```

- [ ] **Step 6: Delete `_species_name_from_sig` and `_step_tuple_from_sig`**

These were only used by the @generated `reactions` and `enzyme_forms` accessors. With those rewritten to walk Mechanism directly, both helpers are dead.

Delete:
- `src/types.jl:1052-1080` (`_species_name_from_sig` + its docstring)
- `src/types.jl:1081-1116` (`_step_tuple_from_sig` + its docstring)

- [ ] **Step 7: grep for any remaining references**

```bash
grep -rn "_species_name_from_sig\|_step_tuple_from_sig" src/ test/
```

Expected output: empty.

- [ ] **Step 8: Run test_accessors.jl perf gate — expect failures**

```bash
julia --project test/test_accessors.jl
```

Expected: the accessor perf gate fails since accessors are no longer 0-allocation (they lift to Mechanism per call). Per Q-005, this test is **explicitly negotiable**.

Either:
- Delete `test/test_accessors.jl` entirely (the accessor perf was an artifact of the @generated implementation, not a user-facing contract), OR
- Re-baseline the test to allow allocation/time for the demoted accessors (only `rate_equation`'s perf is a real contract).

Recommendation: delete the file. Log in `docs/superpowers/refactor-deleted-tests.md`:

```
## 2026-05-30 — F-007 accessor demotion

test/test_accessors.jl asserted zero-allocation / sub-100ns for the 14
@generated accessors over EnzymeMechanism{Sig}. Per Q-005 from the
2026-05-30 audit, none of these accessors are on rate_equation's
runtime hot path — they're called only at @generated body-build time
(compile time) or not at all. The perf gate constrained the
implementation without testing user-facing behavior, so the file is
deleted as part of the singleton demote. The user-facing perf contract
(test_rate_equation_performance in test/test_rate_eq_derivation.jl) is
preserved.
```

- [ ] **Step 9: Delete test_accessors.jl and log**

```bash
git rm test/test_accessors.jl
```

Append the log entry to `docs/superpowers/refactor-deleted-tests.md` (verbatim from Step 8 above).

- [ ] **Step 10: Run rate_equation perf gate explicitly**

```bash
julia --project test/test_rate_eq_derivation.jl 2>&1 | grep -E "(allocations|<100|test_rate_equation_performance)"
```

Expected: 0 allocations / <100ns reported for every spec. **This is the load-bearing perf test; if it fails, abort and investigate before continuing.**

- [ ] **Step 11: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 12: Commit**

```bash
git add src/types.jl test/test_accessors.jl docs/superpowers/refactor-deleted-tests.md
git commit -m "$(cat <<'EOF'
refactor: collapse 14 @generated accessors to plain functions over Mechanism

Findings F-007 + F-008 + F-009 (Cluster A) from the 2026-05-30 audit.
Per Q-005, none of the 14 @generated accessors over EnzymeMechanism{Sig}
are on rate_equation's runtime hot path — they're called only at
@generated body-build time (compile time) or not at all. Converting
to plain functions that lift to Mechanism(em) and walk fields directly
is byte-equivalent for the body builders' purposes.

Deletes:
- _species_name_from_sig (mirror of name(::Species))
- _step_tuple_from_sig (rebuilt opaque tuples)
- test/test_accessors.jl (perf gate; constrained @generated impl;
  per Q-005 not a user-facing contract — rate_equation perf preserved)

Logged in refactor-deleted-tests.md.

rate_equation 0-alloc / <100ns gate stays green.

~258 LOC of @generated → ~70 LOC of plain (~188 LOC saved).
~50 LOC of _species_name_from_sig + _step_tuple_from_sig deleted.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 3.6: Convert `@generated parameters` / `fitted_params` to plain functions

**Files:**
- Modify: `src/rate_eq_derivation.jl:41-93`

Finding F-037.

- [ ] **Step 1: Replace all 4 `@generated parameters` methods + `@generated fitted_params`**

Read `src/rate_eq_derivation.jl` lines 38-93.

Replace L41-93 with plain functions:

```julia
function parameters(em::EnzymeMechanism, ::FullMode)
    mech = Mechanism(em)
    params = _enumerate_parameters_full(mech)
    Tuple((Tuple(name(p, mech) for p in params)..., :E_total))
end

function parameters(em::EnzymeMechanism, ::ReducedMode)
    _, indep = _dependent_param_exprs(typeof(em))
    (indep..., :Keq, :E_total)
end

function parameters(aem::AllostericEnzymeMechanism, ::FullMode)
    am = AllostericMechanism(aem)
    params = _enumerate_parameters_full_allosteric(am)
    names = Symbol[name(p, am) for p in params]
    synth_names = _synthesized_dep_i_names(typeof(catalytic_mechanism(aem)), am)
    if !isempty(synth_names)
        insert_pos = findfirst(p -> p isa Union{Kreg, Lallo}, params)
        idx = insert_pos === nothing ? length(names) + 1 : insert_pos
        splice!(names, idx:(idx - 1), synth_names)
    end
    Tuple((names..., :E_total))
end

function parameters(aem::AllostericEnzymeMechanism, ::ReducedMode)
    _, indep = _dependent_param_exprs(typeof(aem))
    (indep..., :Keq, :E_total)
end

function fitted_params(em::EnzymeMechanism)
    _, indep = _dependent_param_exprs(typeof(em))
    indep
end

function fitted_params(aem::AllostericEnzymeMechanism)
    _, indep = _dependent_param_exprs(typeof(aem))
    indep
end
```

- [ ] **Step 2: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 3: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
refactor: plain parameters / fitted_params functions over Mechanism family

Finding F-037 (Cluster A) from the 2026-05-30 audit. parameters() and
fitted_params() were @generated; per Q-005 these are not on the
runtime hot path (the actual rate_equation hot path stays @generated).
Plain functions are clearer and have identical behavior. ~25 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 3.7: Consolidate AllostericEnzymeMechanism forwarding-accessor block

**Files:**
- Modify: `src/types.jl:1244-1259`

Finding F-010. Per Q-006, all 13 are safely consolidable.

- [ ] **Step 1: Replace the 14-line block with a macro loop**

Read `src/types.jl` lines 1243-1262.

Replace:

```julia
substrates(m::AllostericEnzymeMechanism)         = substrates(catalytic_mechanism(m))
products(m::AllostericEnzymeMechanism)           = products(catalytic_mechanism(m))
reactions(m::AllostericEnzymeMechanism)          = reactions(catalytic_mechanism(m))
equilibrium_steps(m::AllostericEnzymeMechanism)  = equilibrium_steps(catalytic_mechanism(m))
n_steps(m::AllostericEnzymeMechanism)            = n_steps(catalytic_mechanism(m))
enzyme_forms(m::AllostericEnzymeMechanism)       = enzyme_forms(catalytic_mechanism(m))
n_states(m::AllostericEnzymeMechanism)           = n_states(catalytic_mechanism(m))
kinetic_group(m::AllostericEnzymeMechanism, i::Int) =
    kinetic_group(catalytic_mechanism(m), i)
kinetic_groups(m::AllostericEnzymeMechanism)     = kinetic_groups(catalytic_mechanism(m))
steps_in_group(m::AllostericEnzymeMechanism, g)  =
    steps_in_group(catalytic_mechanism(m), g)
stoich_matrix(m::AllostericEnzymeMechanism)      = stoich_matrix(catalytic_mechanism(m))
enzyme_row_range(m::AllostericEnzymeMechanism)   = enzyme_row_range(catalytic_mechanism(m))
metabolite_row_range(m::AllostericEnzymeMechanism) =
    metabolite_row_range(catalytic_mechanism(m))
```

With:

```julia
# All these accessors forward to the catalytic_mechanism. Generated en
# masse to avoid 13 lines of boilerplate.
for fn in (:substrates, :products, :reactions, :equilibrium_steps,
           :n_steps, :enzyme_forms, :n_states, :kinetic_groups,
           :stoich_matrix, :enzyme_row_range, :metabolite_row_range)
    @eval $fn(m::AllostericEnzymeMechanism) = $fn(catalytic_mechanism(m))
end
kinetic_group(m::AllostericEnzymeMechanism, i::Int) =
    kinetic_group(catalytic_mechanism(m), i)
steps_in_group(m::AllostericEnzymeMechanism, g) =
    steps_in_group(catalytic_mechanism(m), g)
```

- [ ] **Step 2: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 3: Commit**

```bash
git add src/types.jl
git commit -m "$(cat <<'EOF'
refactor: consolidate 13 AllostericEnzymeMechanism forwarding accessors via @eval loop

Finding F-010 from the 2026-05-30 audit. Per Q-006, all 13 forwarders
to catalytic_mechanism(m) had no non-forwarding overrides. Replace
with one for-fn @eval loop. ~9 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 3.8: Collapse `_rep_step` 3-method dispatch + delete `_to_mechanism` / `_AnyMech` bridge

**Files:**
- Modify: `src/types.jl:1361-1418`

Finding F-011.

- [ ] **Step 1: Inspect current `_AnyMech` callers**

```bash
grep -n "_AnyMech" src/
```

Expected: usage in the name() chokepoint methods at L1384-1413.

- [ ] **Step 2: Replace the bridge logic**

Read `src/types.jl` lines 1361-1418.

Replace L1361-1382 (the 3-method `_rep_step` + `_to_mechanism` + `_AnyMech`) with:

```julia
# Find the kinetic group containing `step`; return its naming rep.
function _rep_step(step::Step, m::Union{Mechanism, AllostericMechanism})
    fes = _free_enz_set(m)
    for group in steps(m)
        step in group && return _group_rep(group, fes)
    end
    error("Step not found in mechanism: $step")
end
_rep_step(step::Step, em::EnzymeMechanism) = _rep_step(step, Mechanism(em))
_rep_step(step::Step, aem::AllostericEnzymeMechanism) =
    _rep_step(step, AllostericMechanism(aem))

const _AnyMech =
    Union{Mechanism, EnzymeMechanism, AllostericMechanism, AllostericEnzymeMechanism}
```

The `_to_mechanism` function and the bridge alias `_AnyMech` stay because the name() chokepoint methods (`name(p::Kd, m::_AnyMech)`, etc., at L1384-1413) dispatch over all four types. Demoting the singletons means these methods could narrow to `Union{Mechanism, AllostericMechanism}`, but that requires every caller to hold a non-singleton — which only fully lands at the end of Wave 3.

**For now (mid-Wave-3):** keep `_AnyMech` as the 4-type union. The 3-method `_rep_step` collapses to 1 method on the value-context union plus 2 delegating shims (still saves LOC vs. the 3 fully-different bodies).

After Task 3.16 (deleting compile_mechanism callers and Sig conversion entirely), revisit this finding and tighten `_AnyMech` to `Union{Mechanism, AllostericMechanism}`.

- [ ] **Step 3: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 4: Commit**

```bash
git add src/types.jl
git commit -m "$(cat <<'EOF'
refactor: collapse _rep_step 3-method dispatch (partial — defer _AnyMech tighten)

Finding F-011 (partial) from the 2026-05-30 audit. _rep_step now has
one body on Union{Mechanism, AllostericMechanism} plus two delegating
shims for EnzymeMechanism and AllostericEnzymeMechanism. Full collapse
to Union{Mechanism, AllostericMechanism}-only happens after Cluster A
finishes (Task 3.16) when the singletons no longer appear in the
chokepoint dispatch surface. ~10 LOC saved so far.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 3.9: Move `_drop_unbound_regulators` to a new lift utility module

**Files:**
- Modify: `src/types.jl:657-676`

Finding F-004. Per Q-009, the helper is intentional (Mechanism is the enumeration working-representation; unbound regulators are intentional mid-enumeration; the filter belongs at the lift point).

- [ ] **Step 1: Verify the function is still needed**

The function deletes if singleton demote eliminates the lift entirely. But during Wave 3, the lift (`compile_mechanism`) is still called explicitly by users / fitting pipelines. So the filter stays at the lift point.

Specifically: after singleton demote, `compile_mechanism(m::Mechanism)` still exists, still returns an `EnzymeMechanism{Sig}`, and still needs to drop unbound regulators from the Sig.

No code change needed for this finding in Wave 3. Mark as deferred-to-followup (it goes away naturally when the singletons themselves go away in some future refactor, which is out of this audit's scope).

- [ ] **Step 2: Note in the docstring**

Read `src/types.jl` lines 653-676.

Update the docstring above `_drop_unbound_regulators` to clarify the role:

```julia
# A regulator declared on the reaction that no step actually binds does
# not belong in the COMPILED catalytic mechanism's `regulators` list
# (e.g. a dead-end inhibitor before any expansion move binds it). Drop
# such regulators at the `compile_mechanism` boundary. Mechanism
# (the working representation used during enumeration) DOES carry
# unbound regulators intentionally — expansion moves bind them later.
```

- [ ] **Step 3: Commit**

```bash
git add src/types.jl
git commit -m "$(cat <<'EOF'
docs: clarify _drop_unbound_regulators role at lift boundary

Finding F-004 (clarification) from the 2026-05-30 audit. Per Q-009, the
helper is correctly placed at the EnzymeMechanism(::Mechanism) lift
(types.jl:650). Mechanism intentionally carries unbound regulators
during enumeration; the filter runs at compile time to keep them out
of the rate-equation surface. Updated docstring to make this explicit.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 3.10: Update DSL macro emissions to construct Mechanism / AllostericMechanism directly

**Files:**
- Modify: `src/dsl.jl:632-633`
- Modify: `src/dsl.jl:1009-1041` (_build_reg_sites_expr)
- Modify: `src/dsl.jl:1048-1059` (_build_cat_sites_expr)
- Modify: `src/dsl.jl:1137`

Findings F-015, F-016, F-017.

- [ ] **Step 1: dsl.jl:632-633 — emit `Mechanism(...)` directly**

Read `src/dsl.jl` lines 628-634.

```
632:     :(EnzymeRates.EnzymeMechanism(
633:         EnzymeRates.Mechanism($reaction_expr, $groups_expr)))
```

Edit:

```
632:     :(EnzymeRates.Mechanism($reaction_expr, $groups_expr))
```

- [ ] **Step 2: dsl.jl:1009-1041 — rewrite `_build_reg_sites_expr` to emit `Vector{RegulatorySite}`**

Read `src/dsl.jl` lines 1002-1041.

Replace the function body to emit a Vector of RegulatorySite constructor calls instead of a type-parameter tuple:

```julia
"""
Build the regulatory_sites Vector expression for the AllostericMechanism
constructor. Each entry constructs a RegulatorySite with its ligand
list, multiplicity, and per-ligand allo_states.
"""
function _build_reg_sites_expr(allo_regs, reg_site_specs, cat_n)
    tag_of = Dict{Symbol,Symbol}(allo_regs)
    explicit = Set{Symbol}()
    for (_, ligs) in reg_site_specs
        for l in ligs
            l in explicit && error("@allosteric_mechanism: ligand $l " *
                                   "appears in multiple regulatory sites")
            haskey(tag_of, l) ||
                error("@allosteric_mechanism: ligand $l on a " *
                      "`regulatory_site` is not declared in " *
                      "`allosteric_regulators:`")
            push!(explicit, l)
        end
    end

    sites = Tuple{Any,Vector{Symbol}}[]
    for (mult, ligs) in reg_site_specs
        push!(sites, (mult, ligs))
    end
    for (name, _) in allo_regs
        name in explicit && continue
        push!(sites, (cat_n, [name]))
    end

    entries = Expr[]
    for (mult, ligs) in sites
        ligs_vec = :(EnzymeRates.AllostericRegulator[
            $((:(EnzymeRates.AllostericRegulator($(QuoteNode(l)))) for l in ligs)...)])
        states_vec = :(Symbol[$((QuoteNode(tag_of[l]) for l in ligs)...)])
        push!(entries,
              :(EnzymeRates.RegulatorySite($ligs_vec, $mult, $states_vec)))
    end
    :(EnzymeRates.RegulatorySite[$(entries...)])
end
```

- [ ] **Step 3: dsl.jl:1048-1059 — rewrite `_build_cat_sites_expr` to emit `(mult, Vector{Symbol})`**

Read `src/dsl.jl` lines 1043-1059.

Replace:

```julia
"""
Build the (multiplicity, cat_allo_states) pair expression for the
AllostericMechanism constructor. `cat_allo_states` is a dense
`Vector{Symbol}` with one entry per kinetic group in source order.
"""
function _build_cat_sites_expr(cat_n, group_tags)
    for (_, tag) in group_tags
        tag in _ALLOSTERIC_REG_STATES ||
            error("@allosteric_mechanism: catalytic step tag :$tag not in " *
                  "($(_format_state_set(_ALLOSTERIC_REG_STATES)))")
    end
    tag_of = Dict{Int,Symbol}(group_tags)
    n_groups = isempty(group_tags) ? 0 : maximum(g for (g, _) in group_tags)
    states_vec = :(Symbol[$((QuoteNode(get(tag_of, g, :NonequalAI)) for g in 1:n_groups)...)])
    cat_n, states_vec
end
```

Note: this changes the return shape from `Expr(:tuple, cat_n, states_tuple)` to a Julia tuple `(cat_n, states_vec_expr)`. The single caller adapts accordingly:

- [ ] **Step 4: dsl.jl:1137 — emit `AllostericMechanism(...)` directly**

Read `src/dsl.jl` lines 1131-1138.

Replace:

```
1134:     cat_sites_expr = _build_cat_sites_expr(cat_n, group_tags)
1135:     reg_sites_expr = _build_reg_sites_expr(allo_regs, reg_site_specs, cat_n)
1136:
1137:     :(AllostericEnzymeMechanism($cm_expr, $cat_sites_expr, $reg_sites_expr))
```

With:

```
1134:     cat_n_val, cat_states_expr = _build_cat_sites_expr(cat_n, group_tags)
1135:     reg_sites_expr = _build_reg_sites_expr(allo_regs, reg_site_specs, cat_n)
1136:     # cm_expr is :(Mechanism(reaction_expr, groups_expr)) — extract its args:
1137:     # cm_expr.args[1] = :Mechanism (or its qualified form)
1138:     # cm_expr.args[2] = reaction_expr
1139:     # cm_expr.args[3] = groups_expr
1140:     reaction_expr = cm_expr.args[2]
1141:     groups_expr = cm_expr.args[3]
1142:     :(EnzymeRates.AllostericMechanism(
1143:         $reaction_expr, $groups_expr, $cat_states_expr,
1144:         $cat_n_val, $reg_sites_expr))
```

- [ ] **Step 5: Run dsl tests as regression guard**

```bash
julia --project test/test_dsl.jl
```

Expected: green. The macro output types change (Mechanism instead of EnzymeMechanism; AllostericMechanism instead of AllostericEnzymeMechanism), so tests asserting on the constructed type need updates. If failures occur, update test assertions accordingly.

- [ ] **Step 6: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 7: Commit**

```bash
git add src/dsl.jl test/test_dsl.jl
git commit -m "$(cat <<'EOF'
refactor: DSL macros emit Mechanism / AllostericMechanism directly

Findings F-015, F-016, F-017 (Cluster A) from the 2026-05-30 audit.
@enzyme_mechanism now emits Mechanism(reaction, groups) instead of
EnzymeMechanism(Mechanism(...)). @allosteric_mechanism emits
AllostericMechanism(reaction, groups, cat_allo_states, multiplicity,
regulatory_sites) instead of AllostericEnzymeMechanism(...).
_build_reg_sites_expr / _build_cat_sites_expr emit Vector{RegulatorySite}
and Vector{Symbol} instead of type-parameter tuples.

Users now hold the concrete struct directly; compile_mechanism is the
explicit lift to the EnzymeMechanism singleton for the rate_equation
hot path.

~10 LOC saved.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 3.11: Delete `compile_mechanism` lifts in `identify_rate_equation.jl`

**Files:**
- Modify: `src/identify_rate_equation.jl:424, 461, 831, 860`

Finding F-070. Per Q-018, all 4 lifts become no-ops since APIs already accept `AbstractEnzymeMechanism`.

**Important:** the APIs (`rate_equation_string`, `_canonical_rate_eq_hash_data`, `fitted_params`, `FittingProblem`, `_loocv`) currently dispatch on `AbstractEnzymeMechanism`. They need to accept `Mechanism` / `AllostericMechanism` directly too — either by widening to include these types, or by making them subtype `AbstractEnzymeMechanism`.

The cleanest path: **make `Mechanism` and `AllostericMechanism` subtype `AbstractEnzymeMechanism`** (a one-line struct annotation change), so the API stays the same and `compile_mechanism` becomes optional.

- [ ] **Step 1: Subtype Mechanism and AllostericMechanism under AbstractEnzymeMechanism**

Read `src/types.jl` lines 380-410.

Edit L380:

```
380: struct Mechanism
```

To:

```
380: struct Mechanism <: AbstractEnzymeMechanism
```

Edit L411 similarly:

```
411: struct AllostericMechanism
```

To:

```
411: struct AllostericMechanism <: AbstractEnzymeMechanism
```

Note: `AbstractEnzymeMechanism` is defined at L622 of types.jl, which is AFTER the Mechanism/AllostericMechanism struct definitions. Move the abstract type declaration to BEFORE both struct definitions (place at L378 or earlier):

```julia
# At L378 (before Mechanism):
abstract type AbstractEnzymeMechanism end
```

And remove the duplicate at the old L622 location.

- [ ] **Step 2: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green. APIs that dispatched on `AbstractEnzymeMechanism` now accept all four types: Mechanism, AllostericMechanism, EnzymeMechanism, AllostericEnzymeMechanism.

- [ ] **Step 3: Delete 4 compile_mechanism call sites in identify_rate_equation.jl**

Read `src/identify_rate_equation.jl` lines 420-470.

Replace L424:

```
424:                 m = compile_mechanism(mech)
```

With direct usage:

```
424:                 m = mech
```

Same at L461:

```
461:                 m = compile_mechanism(rep.mech)
```

→

```
461:                 m = rep.mech
```

Same at L831:

```
831:         m = compile_mechanism(mech)
```

→

```
831:         m = mech
```

Same at L860:

```
860:     best_mechanism = compile_mechanism(best_mech)
```

→

```
860:     best_mechanism = best_mech
```

- [ ] **Step 4: Update test files that call compile_mechanism**

Per Q-001, ~50 test callers across `test/test_mechanism_enumeration.jl`, `test/test_rate_eq_derivation.jl`, `test/test_canonical_hash_partition.jl`, `test/test_dep_set_invariance.jl`.

Strategy: each `compile_mechanism(m)` call becomes `m` (identity). Use sed:

```bash
for f in test/test_mechanism_enumeration.jl test/test_rate_eq_derivation.jl test/test_canonical_hash_partition.jl test/test_dep_set_invariance.jl; do
    sed -i.bak 's/compile_mechanism(\([^)]*\))/\1/g' "$f"
done
rm test/*.bak
```

This is a simple text substitution. The pattern `compile_mechanism(m)` becomes `m` for any single-argument call. If a test uses `compile_mechanism(reaction)` for some other purpose (passing a reaction instead of a mechanism), it would also be touched — but `compile_mechanism` only has Mechanism / AllostericMechanism methods, so any other usage was already a bug.

- [ ] **Step 5: Update one specific test (Q-012 noted this)**

Read `test/test_mechanism_enumeration.jl` line 4295.

The line asserts `name_map isa Dict{String, String}`. Update if needed (it should still be the case at this stage since Cluster F is deferred — string-keyed projection stays). If the assertion's still valid, leave it.

- [ ] **Step 6: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 7: Commit**

```bash
git add src/types.jl src/identify_rate_equation.jl test/
git commit -m "$(cat <<'EOF'
refactor: subtype Mechanism / AllostericMechanism under AbstractEnzymeMechanism + drop compile_mechanism lifts

Finding F-070 (Cluster A) from the 2026-05-30 audit. Per Q-018, the
four compile_mechanism call sites in identify_rate_equation.jl all
feed APIs (rate_equation_string, _canonical_rate_eq_hash_data,
fitted_params, FittingProblem, _loocv) that dispatch on
AbstractEnzymeMechanism. By making Mechanism and AllostericMechanism
subtype AbstractEnzymeMechanism (one-line change per struct), the
APIs accept all four types and the lifts become identity.

~4 LOC of explicit lifts deleted in src; ~50 lifts deleted in tests.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 3.12: Delete `_raw_param_symbols` singleton-lift forwarder

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:23-24`

Finding F-053.

- [ ] **Step 1: Delete the 2-line forwarder**

Read `src/thermodynamic_constr_for_rate_eq_derivation.jl` lines 20-26.

Delete L23-24:

```
23: _raw_param_symbols(m::EnzymeMechanism) = _raw_param_symbols(Mechanism(m))
24:
```

The Mechanism method at L20-22 stays. With AbstractEnzymeMechanism unification (Task 3.11), the existing Mechanism method covers all use cases.

Actually wait — the Mechanism method specifically takes `Mechanism`, not `AbstractEnzymeMechanism`. We need to either keep the EnzymeMechanism forwarder (in which case it's not dead) or extend the Mechanism method's signature.

After Task 3.11, callers that hold a singleton `em::EnzymeMechanism` still need a way to call `_raw_param_symbols`. The most consistent path: leave the forwarder, OR add `_raw_param_symbols(em::EnzymeMechanism) = _raw_param_symbols(Mechanism(em))` as an explicit shim.

Actually, the existing L23-24 IS that shim. Keep it. **Demote F-053 to "no change" for now** (it deletes naturally if/when the EnzymeMechanism singleton is fully internalized and no callers reach it).

- [ ] **Step 2: No commit (no change)**

Mark F-053 as deferred in the audit followup notes (the audit's optimistic LOC estimate for this finding was wrong — it only deletes once no caller holds a bare EnzymeMechanism, which isn't this wave).

## Task 3.13: Delete Mechanism↔Sig conversion machinery

**Files:**
- Modify: `src/types.jl:500-617`

Finding F-002. Per Q-002, the bridge has only 2 callers (`EnzymeMechanism(::Mechanism)` at L650 and `Mechanism(::EnzymeMechanism)` at L678) plus 1 test (test/test_types.jl:1195-1209).

**CRITICAL:** the `@generated rate_equation` body builder still needs `_mechanism_from_sig(Sig)` to lift the type parameter back to a Mechanism for the body-build walk. Verify before deleting.

- [ ] **Step 1: Grep for `_mechanism_from_sig` callers**

```bash
grep -rn "_mechanism_from_sig\b" src/ test/
```

Expected: callers in `Mechanism(em::EnzymeMechanism{Sig}) where {Sig} = _mechanism_from_sig(Sig)` at types.jl:678 (which is the user-facing lift back).

Additionally: any @generated body that does `mech = Mechanism(M())` at body-build time triggers this lift. Search:

```bash
grep -n "Mechanism(M())\|Mechanism(em)\|Mechanism(m)" src/rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl src/mechanism_enumeration.jl
```

These all use `Mechanism(::EnzymeMechanism)` which calls `_mechanism_from_sig(Sig)`. Therefore the function is still needed.

- [ ] **Step 2: Verify which Sig-helpers can delete vs which must stay**

The bridge has 11 helpers. Of these:

**Must stay** (called by the lift constructors at L650/L678 or by @generated body builders that still operate on Sig):
- `_to_sig(s::Substrate)`, `_to_sig(p::Product)`, `_to_sig(r::AllostericRegulator)`, `_to_sig(c::CompetitiveInhibitor)` (L510-513)
- `_to_sig(r::Residual)` (L515-518)
- `_to_sig(s::Species)` (L520-524)
- `_to_sig(s::Step)` (L526-531)
- `_to_sig(ra::ReactantAtoms)` (L533-536)
- `_to_sig(rm::RegulatorMults)` (L538-541)
- `_to_sig(r::EnzymeReaction)` (L543-547)
- `_sig_of(m::Mechanism)` (L610-613)
- `_metabolite_from_sig`, `_residual_from_sig`, `_species_from_sig`, `_step_from_sig`, `_reactant_atoms_from_sig`, `_regulator_mults_from_sig`, `_reaction_from_sig`, `_steps_from_sig` (L549-608)
- `_mechanism_from_sig` (L615-618)

**Could delete** if no caller remains — but all 11 are part of the encode/decode chain, so they're all live as long as `EnzymeMechanism{Sig}` is used as a `@generated` dispatch type.

Conclusion: **DEFER F-002 entirely**. The Sig conversion machinery is load-bearing as long as `EnzymeMechanism{Sig}` exists as the @generated dispatch type. Deleting it requires a deeper refactor (the deferred direction-symmetry effort, or a different @generated dispatch strategy entirely).

- [ ] **Step 3: Mark F-002 as out-of-scope for Wave 3**

The audit's optimistic LOC estimate for F-002 was wrong — the bridge cannot delete without redesigning the @generated dispatch model itself. Document this as a deferred finding alongside Cluster F.

- [ ] **Step 4: No commit (no change)**

Update the audit findings doc later (Wave 4 doc sweep) to reflect this revised disposition.

## Task 3.14: Rewrite `show` methods for EnzymeMechanism / AllostericEnzymeMechanism

**Files:**
- Modify: `src/types.jl:851-972`

Finding F-006.

- [ ] **Step 1: Update `show(::EnzymeMechanism)` to walk Mechanism directly**

Read `src/types.jl` lines 851-909.

Replace the function. The current implementation walks `reactions(em)` (which now lifts to Mechanism and walks); cleaner to walk Mechanism directly:

```julia
function Base.show(io::IO, em::EnzymeMechanism)
    m = Mechanism(em)
    flat = _flat_steps(m)
    enz_set = Set(enzyme_forms(em))
    _arrow(is_eq) = is_eq ? " ⇌ " : " <--> "

    chain_segments = String[]
    chain_arrows = String[]
    is_linear = !isempty(flat)
    current = nothing
    for (i, (s, _)) in enumerate(flat)
        lhs_str = is_binding(s) && is_equilibrium(s) ?
            string(name(from_species(s)), " + ", name(bound_metabolite(s))) :
            string(name(from_species(s)))
        rhs_str = string(name(to_species(s)))
        # Direction matching to chain steps together
        e_l = name(from_species(s))
        e_r = name(to_species(s))
        if i == 1
            push!(chain_segments, lhs_str)
            push!(chain_arrows, _arrow(is_equilibrium(s)))
            push!(chain_segments, rhs_str)
            current = e_r
        elseif current == e_l
            push!(chain_arrows, _arrow(is_equilibrium(s)))
            push!(chain_segments, rhs_str)
            current = e_r
        elseif current == e_r
            push!(chain_arrows, _arrow(is_equilibrium(s)))
            push!(chain_segments, lhs_str)
            current = e_l
        else
            is_linear = false
            break
        end
    end

    if is_linear
        print(io, "EnzymeMechanism: ", chain_segments[1])
        for k in 2:length(chain_segments)
            print(io, chain_arrows[k-1], chain_segments[k])
        end
    else
        print(io, "EnzymeMechanism (", length(flat), " steps, ",
              length(enz_set), " enzyme forms):")
        for (s, _) in flat
            lhs_str = is_binding(s) ?
                string(name(from_species(s)), " + ", name(bound_metabolite(s))) :
                string(name(from_species(s)))
            print(io, "\n  ", lhs_str, _arrow(is_equilibrium(s)),
                      name(to_species(s)))
        end
    end
    regs = regulators(em)
    if !isempty(regs)
        print(io, " | regulators: ", join(regs, ", "))
    end
end
```

- [ ] **Step 2: Rewrite `_format_allo_step_groups`**

Read `src/types.jl` lines 910-951.

Replace the function. Walk `cm`'s `_flat_steps(cm)` directly instead of `reactions(cm)`:

```julia
"""
Render the catalytic mechanism's steps as multi-line text, grouping
steps that share a kinetic_group with parens and a single `:: Tag`
annotation. Mirrors `@allosteric_mechanism` macro syntax.
"""
function _format_allo_step_groups(
    io::IO, cm::EnzymeMechanism,
    m::AllostericEnzymeMechanism,
)
    cm_mech = Mechanism(cm)
    flat = _flat_steps(cm_mech)
    _arrow(is_eq) = is_eq ? " ⇌ " : " <--> "

    groups_seen = Int[]
    group_to_step_idxs = Dict{Int,Vector{Int}}()
    for (i, (_, g)) in enumerate(flat)
        if !haskey(group_to_step_idxs, g)
            push!(groups_seen, g)
            group_to_step_idxs[g] = Int[]
        end
        push!(group_to_step_idxs[g], i)
    end

    for g in groups_seen
        idxs = group_to_step_idxs[g]
        tag = cat_allo_state(m, g)
        _format_one_step(s) = begin
            lhs_str = is_binding(s) ?
                string(name(from_species(s)), " + ", name(bound_metabolite(s))) :
                string(name(from_species(s)))
            "$lhs_str$(_arrow(is_equilibrium(s)))$(name(to_species(s)))"
        end
        if length(idxs) == 1
            (s, _) = flat[idxs[1]]
            print(io, "\n  ", _format_one_step(s), " :: ", tag)
        else
            print(io, "\n  (")
            for (k, i) in enumerate(idxs)
                k > 1 && print(io, ", ")
                (s, _) = flat[i]
                print(io, _format_one_step(s))
            end
            print(io, ") :: ", tag)
        end
    end
end
```

- [ ] **Step 3: Rewrite `show(::AllostericEnzymeMechanism)`**

Read `src/types.jl` lines 953-972.

Replace the function:

```julia
function Base.show(io::IO, m::AllostericEnzymeMechanism)
    cm = catalytic_mechanism(m)
    print(io, "AllostericEnzymeMechanism (cat_n=",
          catalytic_multiplicity(m))
    rs = regulatory_sites(m)
    if !isempty(rs)
        print(io, ", ", length(rs), " reg sites")
    end
    print(io, "):")
    _format_allo_step_groups(io, cm, m)
    for (i, (ligs, mult, reg_allo_states)) in enumerate(rs)
        print(io, "\n  reg site $i (n=", mult, "): ", join(ligs, ", "))
        print(io, " [")
        print(io, join(("$(n)::$(t)"
                        for (n, t) in zip(ligs, reg_allo_states)),
                       ", "))
        print(io, "]")
    end
end
```

- [ ] **Step 4: Run tests + visual check**

```bash
julia --project test/test_types.jl
```

Plus a quick interactive sanity check:

```bash
julia --project -e '
using EnzymeRates
m = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
    end
end
show(stdout, m); println()
'
```

Expected: output reads "EnzymeMechanism: E + S ⇌ ES <--> EP ⇌ E + P" (or similar linear render).

- [ ] **Step 5: Commit**

```bash
git add src/types.jl
git commit -m "$(cat <<'EOF'
refactor: show methods walk Mechanism / AllostericMechanism directly

Finding F-006 from the 2026-05-30 audit. show(::EnzymeMechanism) and
show(::AllostericEnzymeMechanism) lift to the concrete struct family
and walk Step fields directly, removing dependence on the opaque-tuple
reactions() output. No LOC saving — pure clarity refactor.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 3.15: Wave 3 verification

- [ ] **Step 1: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 2: Verify rate_equation perf gate**

```bash
julia --project test/test_rate_eq_derivation.jl 2>&1 | grep -E "(allocations|<100|test_rate_equation_performance)"
```

Expected: 0 allocations / <100ns reported for every spec.

- [ ] **Step 3: Measure LOC delta**

Run the baseline-counter from Wave 1 Task 1.9 Step 2.

Expected: ~5,000 non-comment non-doc LOC (down ~700 from baseline).

This is short of the audit's ~5,000 LOC target (~13.7% reduction). Some findings demoted to no-change during this wave (F-002, F-053), explaining the gap. The honest measured number is the deliverable.

# Wave 4 — Doc-hygiene sweep + audit-doc updates

Final wave. No LOC reduction. Improves discoverability + CLAUDE.md compliance.

## Task 4.1: Convert 19 function-leading `#`-comments in dsl.jl to docstrings

**Files:**
- Modify: `src/dsl.jl` — see locations in F-018

Finding F-018.

- [ ] **Step 1: For each location, convert `#`-block above the function/struct to a `"""docstring"""` attached to the definition**

For each of these 19 locations in `src/dsl.jl`, take the `#`-comment block immediately above the function/struct and convert to a docstring:

- L36-43 → docstring on `_parse_reaction_block`
- L96-97 → docstring on `_parse_atom_bracket_entries`
- L137-139 → docstring on `_parse_regulator_entries`
- L164 → docstring on `_parse_multiplicity_tuple`
- L186-187 → docstring on `_build_reactants_expr`
- L205 → docstring on `_atoms_pairs_expr`
- L217-220 → docstring on `_build_regulators_expr`
- L246-247 → docstring on `_parse_step_side_terms`
- L257 → docstring on `_step_side_term_info`
- L282-284 → docstring on `_call_form_term_info`
- L341-352 → docstring on `_StepSideTerm` struct
- L369-374 → docstring on `_is_conformation_shape`
- L378-381 → already deleted in Task 1.8 (skip)
- L401-402 → docstring on `_walk_residual_expr`
- L572-577 → docstring on `_build_mechanism_expr`
- L636-639 → docstring on `_build_step_expr`
- L656-658 → docstring on `_split_side`
- L684-685 → docstring on `_species_expr_from_term`
- L710-711 → docstring on `_metabolite_expr`

Pattern for each conversion. Take a block like:

```
# Parse one step side into a Vector{_StepSideTerm}, preserving structural
# info (Call-form decomposition, declared-metabolite role) for emission.
function _parse_step_side_terms(expr, declared_mets::Set{Symbol})
    ...
```

Convert to:

```
"""
Parse one step side into a `Vector{_StepSideTerm}`, preserving structural
info (Call-form decomposition, declared-metabolite role) for emission.
"""
function _parse_step_side_terms(expr, declared_mets::Set{Symbol})
    ...
```

Note: Julia docstrings ARE rendered like Markdown, so `Vector{_StepSideTerm}` (with backticks) renders nicely in `?_parse_step_side_terms` output.

- [ ] **Step 2: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green. Pure doc edits.

- [ ] **Step 3: Commit**

```bash
git add src/dsl.jl
git commit -m "$(cat <<'EOF'
docs: convert 19 function-leading #-comments to docstrings in dsl.jl

Finding F-018 from the 2026-05-30 audit. CLAUDE.md "Code Comments"
rule: docstrings attached to definitions improve discoverability via
?fnname and IDE doc-on-hover. 19 functions/structs in dsl.jl had
explanatory #-comment blocks immediately above their definitions —
converted to proper docstrings.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 4.2: Convert sym_poly comment-as-docstring instances

**Files:**
- Modify: `src/sym_poly_for_rate_eq_derivation.jl:1-2, 52-55, 83-84, 236-237`

Finding F-034.

- [ ] **Step 1: Convert L1-2 file-level comment to a module-attached docstring**

Since this file is `include`d (not a separate module), attach the docstring to no definition directly. Instead, leave the file-level comment as-is OR convert to a doc-string-style explanation in the included file.

Recommended: leave as `#`-comments (file-level) but rewrite to follow the CLAUDE.md "ABOUTME:" convention:

```
# ABOUTME: Lightweight symbolic polynomial type (POLY = Dict{MONO, Rational{Int}})
# ABOUTME: for compile-time rate equation derivation.
```

- [ ] **Step 2: Convert L52-55 to docstring on `sym_det`**

Read `src/sym_poly_for_rate_eq_derivation.jl` lines 52-58.

Take:

```
# Cofactor determinant expansion for symbolic matrices.
# Checks intermediate term count against MAX_RATE_EQUATION_TERMS to
# abort early for mechanisms whose rate equations would be too large.
function sym_det(M::Matrix{POLY}, n::Int)
```

Convert to:

```
"""
Cofactor determinant expansion for symbolic matrices. Checks
intermediate term count against `MAX_RATE_EQUATION_TERMS` to abort
early for mechanisms whose rate equations would be too large.
"""
function sym_det(M::Matrix{POLY}, n::Int)
```

- [ ] **Step 3: Convert L83-84 to docstring on `_poly_to_expr`**

Apply the same pattern.

- [ ] **Step 4: Convert L236-237 to docstring on `substitute_params_expr`**

Apply the same pattern.

- [ ] **Step 5: Run tests + commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add src/sym_poly_for_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
docs: convert sym_poly comments to docstrings + ABOUTME header

Finding F-034 from the 2026-05-30 audit. File-level header now follows
CLAUDE.md ABOUTME convention. Three function-leading comments
(sym_det, _poly_to_expr, substitute_params_expr) converted to
docstrings.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 4.3: Trim minor stale docstrings + comment-as-docstring residuals

**Files:**
- Modify: `src/types.jl:679-704, 911-913` (F-005)
- Modify: `src/sym_poly_for_rate_eq_derivation.jl:122-129, 187-196, 250-258` (F-035 residuals not covered in 4.2)

Findings F-005 + F-035 residuals.

- [ ] **Step 1: types.jl L679-704 — trim AllostericEnzymeMechanism docstring**

Read `src/types.jl` lines 679-704.

The current docstring describes CatSites/RegSites internal type-parameter layout. After Task 3.10, users emit AllostericMechanism directly; AllostericEnzymeMechanism is internal. Trim the docstring to:

```julia
"""
    AllostericEnzymeMechanism{CatalyticMech, CatSites, RegSites}

Internal singleton type used by `@generated rate_equation` for
allosteric MWC enzymes. User code constructs and inspects via
`AllostericMechanism`; `compile_mechanism(am)` produces this opaque
fast-path handle. The three type parameters encode the catalytic
mechanism plus per-site allosteric data.
"""
```

- [ ] **Step 2: sym_poly L122-129, L187-196 — verify docstrings are still accurate**

Read `src/sym_poly_for_rate_eq_derivation.jl` lines 122-130 and 187-197.

These are existing docstrings (`_nest_binary` and `build_power_expr`). Verify they're still accurate after all the upstream refactors. No edits expected — sym_poly is mostly untouched.

- [ ] **Step 3: Run tests + commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add src/types.jl
git commit -m "$(cat <<'EOF'
docs: trim AllostericEnzymeMechanism docstring after demote

Finding F-005 from the 2026-05-30 audit. With AllostericEnzymeMechanism
demoted to an internal compile artifact (Tasks 3.7, 3.10, 3.11), the
docstring describing CatSites/RegSites internal layout is no longer
user-facing. Trim to a brief role description.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 4.4: Update the audit findings doc with revised dispositions

**Files:**
- Modify: `docs/superpowers/2026-05-30-refactor-audit-findings.md`

This task captures findings whose dispositions changed during Wave 3 (F-002 deferred, F-053 no-change).

- [ ] **Step 1: Add a "Implementation outcome" section to the findings doc**

Read `docs/superpowers/2026-05-30-refactor-audit-findings.md`.

Append at the end:

```markdown
---

## Implementation outcome — 2026-XX-XX

This section is updated at the end of Wave 4 with the actual outcomes
of implementing the audit findings.

**LOC reduction achieved:** measured at end-of-Wave-4 in
`docs/superpowers/scratch-audit-impl-loc.md` (or similar).

**Findings disposition revisions:**

- **F-002 (delete Sig conversion machinery):** DEFERRED. The Sig
  conversion bridge is load-bearing as long as
  `EnzymeMechanism{Sig}` exists as the `@generated rate_equation`
  dispatch type. Deletion requires redesigning the @generated dispatch
  model itself (e.g., switching to a Mechanism-keyed dispatch table
  with manual caching, which would break the existing 0-alloc/<100ns
  invariant via Julia's method-table caching). Out of scope for this
  audit.

- **F-053 (delete `_raw_param_symbols(::EnzymeMechanism)` forwarder):**
  NO CHANGE. The 2-line forwarder is needed for callers that hold a
  bare `EnzymeMechanism` (after Task 3.11, `Mechanism` and
  `AllostericMechanism` subtype `AbstractEnzymeMechanism` so most
  callers can use Mechanism directly, but the forwarder bridges the
  EnzymeMechanism case). Deletes naturally when no caller reaches a
  bare EnzymeMechanism (further-future refactor).

**Cluster F (string-keyed projection):** REMAINS DEFERRED pending the
direction-symmetry constraint-resolution refactor.

**Approach C for F-014 (aggressive parser tighten):** NOT TAKEN.
Approach B (inline) landed; Approach C requires Denis sign-off and was
not granted.
```

Replace the date placeholder with the actual completion date.

- [ ] **Step 2: Measure final LOC and update Executive Summary**

Run the baseline-counter from Wave 1 Task 1.9 Step 2 one more time.

Update the findings doc's "Estimated savings" table in §1 with the actual measured numbers.

- [ ] **Step 3: Commit**

```bash
git add docs/superpowers/2026-05-30-refactor-audit-findings.md
git commit -m "$(cat <<'EOF'
docs: audit findings — record implementation outcome

Final pass over the 2026-05-30 audit findings. Records:
- Actual LOC reduction achieved
- F-002 and F-053 disposition revisions (both deferred / no-change
  due to load-bearing Sig dispatch for @generated rate_equation)
- Cluster F still deferred
- Approach C for F-014 not taken (Approach B sufficed)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

## Task 4.5: Final verification

- [ ] **Step 1: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green.

- [ ] **Step 2: Verify rate_equation perf gate**

```bash
julia --project test/test_rate_eq_derivation.jl 2>&1 | grep -E "(allocations|<100|test_rate_equation_performance)"
```

Expected: 0 allocations / <100ns for every spec.

- [ ] **Step 3: Measure final LOC**

```bash
for f in src/*.jl; do
    total=$(wc -l < "$f")
    nccnd=$(awk '
        BEGIN { in_ds=0 }
        { line=$0; n=gsub(/"""/,"\"\"\"",line); stripped=line
          sub(/^[[:space:]]+/,"",stripped); sub(/[[:space:]]+$/,"",stripped)
          if (in_ds) { if (n%2==1) in_ds=0; next }
          if (n%2==1) { in_ds=1; next }
          if (stripped=="") next; if (substr(stripped,1,1)=="#") next; count++ }
        END { print count+0 }' "$f")
    printf "%6d  %6d  %s\n" "$total" "$nccnd" "$f"
done
```

Compare to baseline (5,706 non-comment non-doc LOC). Record the delta.

- [ ] **Step 4: Print summary**

Print a one-paragraph summary:

> "Refactor audit implementation complete. <N> commits across Waves 1-4. <X> LOC saved non-comment non-doc (<P>% of baseline 5,706). <K> findings deferred: F-002, F-053, and all of Cluster F. rate_equation 0-alloc/<100ns gate stays green. Cluster F unblocks once the direction-symmetry refactor lands."

---

## Findings explicitly out of scope for this plan

These findings were identified in the audit but are NOT implemented by this plan:

- **F-054 (FittingProblem retype):** marked Low confidence in the audit; no LOC saving; the `AbstractEnzymeMechanism` signature still works after Task 3.11 makes Mechanism / AllostericMechanism subtypes. Skip.
- **F-055 (delete `compile_mechanism`):** deferred with F-002 — `compile_mechanism` is still needed as the explicit lift to `EnzymeMechanism{Sig}` for the @generated rate_equation hot path. The lift is no longer required at most call sites (Task 3.11 deletes them), but the function itself stays.
- **F-058 (remove three no-op expansion-move dispatches):** Low confidence, optional micro-cleanup, ~3 LOC. The current uniform-dispatch approach (every expand move called on every input, no-ops return empty) is arguably clearer than type-aware dispatching. Skip.
- **F-062 (collapse two `_enumerate_all_parameters_with_i_state` methods):** Low priority, ~3 LOC, parametric-type dispatch tradeoff. Skip.
- **All of Cluster F (F-013, F-064, F-066, F-068, F-069, F-071):** DEFERRED — blocked on the direction-symmetry refactor providing first-class Parameter representation for synth-deps. Track in a follow-up plan after that refactor lands.

Total LOC saving NOT pursued: ~10 LOC for the in-scope-but-skipped findings + ~100 LOC for the deferred Cluster F = ~110 LOC.

Adjusted realistic target: ~670 LOC (down from the audit's ~780 LOC) = ~11.7% reduction. Still a substantial simplification.

## Self-Review Checklist (for the executor before claiming done)

After Task 4.4, run this checklist:

- [ ] **All commits landed in order**: `git log --oneline main..HEAD` shows 30-40 commits in wave-then-task order
- [ ] **Tests green**: `julia --project -e 'using Pkg; Pkg.test()'` succeeds
- [ ] **rate_equation perf gate**: 0 alloc / <100ns for every spec in MECHANISM_TEST_SPECS
- [ ] **Deleted tests logged**: `docs/superpowers/refactor-deleted-tests.md` has entries for `test_accessors.jl` and any test-blocked helper assertions replaced in Task 2.2
- [ ] **Findings doc updated**: §1 executive summary has actual LOC numbers; "Implementation outcome" section is filled in
- [ ] **No stale "Stage X" / "Phase Y" references in src**: `grep -rn "Stage [0-9]\|Phase [0-9]" src/` is empty (or only in legitimate user-facing error messages)
- [ ] **No lingering r/t variable shortcuts referring to old R/T allosteric notation**: `grep -n "\br_name\b\|\bt_name\b\|\br_str\b\|\bt_str\b" src/` is empty
- [ ] **`compile_mechanism` callers reduced**: `grep -rn "compile_mechanism" src/ test/` shows only the function definition and internal users (no production lifts in identify_rate_equation.jl)
- [ ] **AllostericEnzymeMechanism forwarding block consolidated**: `grep -c "AllostericEnzymeMechanism) *=" src/types.jl` is ≤ 3 (the @eval loop + 2 explicit special-case methods)

Any failures here: investigate before declaring done.
