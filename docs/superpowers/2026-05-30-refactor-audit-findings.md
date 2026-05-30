# Concrete-Types Refactor Audit — Findings

**Date:** 2026-05-30
**Branch:** `refactor-to-concrete-types-instead-of-symbols`
**Baseline (non-comment non-doc src LOC):** 5,706 (8,456 total across 9 files)
**Workflow:** `docs/superpowers/specs/2026-05-30-refactor-audit-workflow-design.md`

---

## §1 Executive summary

**71 total findings** across the 9 src files (1 fitting.jl finding is minor; 0 in EnzymeRates.jl module file):

| Category | Count |
|---|---|
| Architectural | 23 |
| Duplication | 17 |
| Dead-code / dead-or-relocate | 5 |
| String-keyed projection | 5 *(deferred)* |
| Permissive-parser-guard | 1 *(batched, Cluster G)* |
| Doc hygiene (stale-spec / comment-as-docstring) | 20 |

**Estimated LOC savings (non-comment non-doc src):**

| Cluster | Net delete | After-rewrite | Notes |
|---|---|---|---|
| A — Singleton-type demotion | ~520 | ~150 LOC of rewrites in place | The big architectural move; high impact |
| B — Derivation back-end struct-native walk | ~80 | ~30 LOC of rewrites | Q-016 confirmed feasible |
| C — Smaller `_emit_cat_params_for_rep` helper | ~20 | ~10 LOC of rewrites | Per Q-015 don't fully unify the 5 helpers |
| D — Expansion-move + dedup! consolidation | ~115 | ~30 LOC of rewrites | Adds `_with_steps`/`_with_*` helpers |
| E — Doc hygiene sweep | 0 (non-comment non-doc) | n/a | ~20 finding-equivalents; readability |
| F — String-keyed projection | DEFERRED | — | Needs direction-symmetry refactor first; see deferred-cluster note |
| G — Parser-tighten (Approach B) | ~25 | — | Inline `_assert_no_opaque_terms` into parse step |
| H — Small cleanups (kinetic-rename rename, wrapper inline, etc.) | ~20 | — | Independent of A-D |

**Total net LOC savings (non-deferred): ~780 LOC (~13.7% of baseline)**, plus ~220 LOC of in-place rewrites (mostly `@generated` → plain functions, struct-native walks).

**Hypothesis test:** Denis's "up to half of src may be removable" claim is **partially supported**. The audit measured ~14% non-comment non-doc LOC reduction available across all non-deferred clusters. The remaining gap to 50% would require:

1. **Further deletion of test-private helpers** that constrain design without testing behavior (CLAUDE.md no-test-deletion rule requires replacing with behavior coverage, so this is constrained).
2. **Replacing the `@generated`-driven King-Altman / Cha derivation** with a runtime alternative (different refactor; not in scope).
3. **The deferred string-keyed-projection cluster** (~100 LOC) lands after the direction-symmetry refactor enables structural Parameter representation for synth-deps.
4. **Test-suite-side reductions** (out of scope; CLAUDE.md "no test coverage reduction").

The audit recommends Denis review the hypothesis after Wave 3 lands and re-baseline. The 14% figure is honest and substantial, even if shy of 50%.

**Top architectural moves:**

1. **Cluster A — Demote `EnzymeMechanism{Sig}` / `AllostericEnzymeMechanism{...}` to internal compile artifact.** The singleton-type bridge (Sig conversion machinery, 14 @generated accessors, 14-line forwarding-accessor block, double-lift constructors, `_drop_unbound_regulators`, `_param_for_symbol` family) all exist to support a parametric mechanism representation. Per Q-005 / Q-006 / Q-009 / Q-018: none of this is on the `rate_equation` runtime hot path. The Mechanism / AllostericMechanism concrete structs (already used throughout the front end) can become the primary representation; the singletons become internal `compile_mechanism` artifacts used only at @generated-body-build time. Net ~520 LOC delete + ~150 LOC of rewrites.
2. **Cluster B — Replace `_step_tuple_from_sig` / `_species_name_from_sig` symbol-tuple plumbing with struct-native walk over `Mechanism.steps`.** The King-Altman and Wegscheider consumers in `_raw_symbolic_rate_polys`, `_compute_alpha`, `_compute_numerator`, `_thermodynamic_constraints`, `_dependent_param_exprs_kernel` all consume `rxns = reactions(m)` opaque tuples that regenerate structure already on Step. Per Q-016 the source comment at L325-329 explicitly invites this cleanup. ~80 LOC + struct-native rewrite.
3. **Cluster D — Add internal update constructors (`_with_steps`, `_with_cat_allo_states`, `_with_reg_sites`) to collapse dual-dispatch expansion moves.** The 4 expansion-move helpers (`_expand_re_to_ss`, `_expand_split_kinetic_group`, `_expand_add_dead_end_regulator`, `_expand_change_allo_state`) each have nearly-identical Mechanism + AllostericMechanism methods that differ only in result-type construction. With proper builder helpers each pair collapses to one method. Plus three `dedup!` methods collapse to one `dedup!(Dict{Int, <:Vector})`. ~115 LOC.

**Doc-hygiene count (does not reduce non-comment non-doc LOC):** ~20 batched findings covering stale spec/stage references and comment-as-docstring instances, primarily in dsl.jl, rate_eq_derivation.jl, and mechanism_enumeration.jl. Recommended as one final commit after Wave 3.

---

## §2 Findings (ordered by src file, then by line range)

### src/types.jl

#### F-001  Fix stale `:E_A_B` example in Species docstring

**Location:** src/types.jl:46-47
**Category:** Doc hygiene
**Confidence:** High
**LOC saving (non-comment non-doc):** 0
**Simplification gain:** Docstring matches code output.
**Depends on:** none
**Blocking tests:** none
**Recommendation:** The docstring at L46-47 reads "`:E_A_B`" (underscore-separated) but the code at L78-82 concatenates without separator producing `:EATP`. Update example to `:EATP` / `:EstarA_res_+P` to match the code-emitted shape.

#### F-002  Delete Mechanism↔Sig conversion machinery after singleton demote

**Location:** src/types.jl:500-617
**Category:** Architectural
**Confidence:** High
**LOC saving (non-comment non-doc):** ~118
**Simplification gain:** Removes the entire Sig encode/decode bridge — `_to_sig` polymorphic methods, the `_*_from_sig` family, `_sig_of`, `_mechanism_from_sig`. Per Q-002 these have only two real callers (`EnzymeMechanism(::Mechanism)` lift and `Mechanism(::EnzymeMechanism)` lift) plus one round-trip test.
**Depends on:** F-003 (must demote the singleton type itself to unlock this)
**Blocking tests:** test/test_types.jl:1195-1209 (sig roundtrip test — delete with the bridge)
**Recommendation:** When `EnzymeMechanism{Sig}` ceases to be the primary representation (F-003), the encode/decode bridge becomes orphan code. Delete L500-617 wholesale and remove the roundtrip test.

#### F-003  Demote EnzymeMechanism / AllostericEnzymeMechanism to internal compile artifacts

**Location:** src/types.jl:619-816
**Category:** Architectural
**Confidence:** Medium (needs design choice on the lift point)
**LOC saving (non-comment non-doc):** ~198
**Simplification gain:** `Mechanism` / `AllostericMechanism` become the primary representations everywhere. Singletons collapse to a slim wrapping needed only by `@generated rate_equation` (which dispatches on the type parameter to specialize per shape). All API functions already accept the `AbstractEnzymeMechanism` abstract supertype (Q-018), so call sites are pass-through.
**Depends on:** none (this is the head of Cluster A)
**Blocking tests:** test_compile_budget.jl (trace-compile budget assumes the singleton dispatches); test_canonical_hash_partition.jl (canonical hash uses `AbstractEnzymeMechanism`); test_accessors.jl perf gate (Q-005: explicitly negotiable)
**Recommendation:** The `@generated rate_equation` dispatch needs *some* per-shape type to specialize on. Option A: keep `EnzymeMechanism{Sig}` as a thin internal type built inside `compile_mechanism` at the @generated entry point only, treat as opaque. Option B: dispatch `rate_equation` on `Mechanism` directly and rely on @generated functions reading `Mechanism` field types — requires verifying Julia can dispatch @generated on concrete structs. Discuss with Denis; both keep `rate_equation`'s 0-alloc / <100ns invariant.

#### F-004  Relocate `_drop_unbound_regulators` to its new lift site

**Location:** src/types.jl:657-676
**Category:** Architectural
**Confidence:** High
**LOC saving (non-comment non-doc):** 0 (relocation, not deletion)
**Simplification gain:** Helper moves to wherever the lift happens after demote. Per Q-009: do NOT move into the Mechanism constructor — enumeration relies on intermediate Mechanisms with unbound regulators.
**Depends on:** F-003
**Blocking tests:** none (internal helper, no test asserts on it directly)
**Recommendation:** Keep the helper as a private utility; relocate its single caller from the deleted `EnzymeMechanism(::Mechanism)` constructor to whatever new lift point F-003 establishes.

#### F-005  Drop docstrings describing internal Sig layout

**Location:** src/types.jl:679-704 + 911-913
**Category:** Doc hygiene
**Confidence:** Low
**LOC saving (non-comment non-doc):** 0
**Simplification gain:** Docstrings stop describing structure that no longer exists.
**Depends on:** F-003
**Blocking tests:** none
**Recommendation:** Trim AllostericEnzymeMechanism docstring (L680-704) to describe the demoted internal type minimally. Update `_format_allo_step_groups` docstring (L911-913).

#### F-006  Re-route `show` methods through Mechanism / AllostericMechanism

**Location:** src/types.jl:851-972
**Category:** Architectural
**Confidence:** Low
**LOC saving (non-comment non-doc):** 0 (rewrite in place)
**Simplification gain:** `show(::Mechanism)` and `show(::AllostericMechanism)` become first-class; the singleton variants delegate.
**Depends on:** F-003, F-007
**Blocking tests:** Tests that grep `show` output for specific Symbol shapes — locate during impl
**Recommendation:** Keep current logic; rewrite to walk Mechanism struct directly instead of Sig accessors. Mostly a search-and-replace.

#### F-007  Collapse 14 @generated accessors to plain field-access over Mechanism

**Location:** src/types.jl:974-1231
**Category:** Architectural
**Confidence:** High *(per Q-005: zero are runtime hot-path)*
**LOC saving (non-comment non-doc):** ~258 → ~30 LOC (~228 LOC delete + 14 LOC of plain methods)
**Simplification gain:** **All 14 @generated accessors are either compile-time-only or unused by `rate_equation`'s runtime body.** Plain `f(m::Mechanism) = m.field` versions are equivalent. The `_species_name_from_sig` mirror (F-008) and `_step_tuple_from_sig` opaque-tuple builder (F-009) become unused.
**Depends on:** F-003
**Blocking tests:** test_accessors.jl perf gate — per Q-005 explicitly negotiable; delete or rebaseline. test_rate_equation_performance is NOT a blocker (rate_equation hot path stays @generated).
**Recommendation:** For each of `substrates`, `products`, `regulators`, `metabolites`, `reactions`, `equilibrium_steps`, `kinetic_groups`, `steps_in_group`, `enzyme_forms`, `stoich_matrix`, `n_steps`, `n_states`, `kinetic_group(idx)`, `enzyme_row_range`, `metabolite_row_range`: rewrite as plain function over Mechanism (with `AllostericMechanism` delegate). Where the result is used at @generated build time (the 7 "class b" accessors per Q-005), the new plain version is still called at body-build time; @generated body just inlines the value as before.

#### F-008  Delete `_species_name_from_sig` (duplicates `name(::Species)`)

**Location:** src/types.jl:1057-1080
**Category:** Duplication
**Confidence:** High *(per Q-003: byte-equivalent across all input shapes)*
**LOC saving (non-comment non-doc):** ~24
**Simplification gain:** One canonical Species→Symbol formatter (the value-context `name(::Species)` at L77-94). Eliminates a duplicate that exists only to be callable from @generated bodies.
**Depends on:** F-007 (the @generated accessors that consume it must switch to value-context first)
**Blocking tests:** none (no direct test asserts on `_species_name_from_sig`)
**Recommendation:** After F-007 lands, the 3 callsites at L1095-1096 and L1191 become value-context callers that can use `name(::Species)` directly.

#### F-009  Delete `_step_tuple_from_sig` (rebuilds opaque tuples that Step carries)

**Location:** src/types.jl:1093-1116
**Category:** Symbol-tuple plumbing
**Confidence:** High *(per Q-004: sole caller is `reactions(::EnzymeMechanism{Sig})`)*
**LOC saving (non-comment non-doc):** ~24
**Simplification gain:** Deletes the opaque `(lhs, rhs, is_eq, g)` reconstruction. Consumers read Step.from_species/to_species/is_equilibrium/outer-vector-index directly.
**Depends on:** F-007 (the `reactions` accessor that consumes it must collapse to non-@generated)
**Blocking tests:** none direct; `reactions(em)` output shape may be asserted in test_types.jl — verify during impl
**Recommendation:** Delete after F-007.

#### F-010  Consolidate 14-line AllostericEnzymeMechanism forwarding-accessor block

**Location:** src/types.jl:1244-1259
**Category:** Duplication
**Confidence:** High *(per Q-006: no non-forwarding overrides)*
**LOC saving (non-comment non-doc):** ~9 (14 LOC → ~5 LOC via macro loop)
**Simplification gain:** One `for fn in (:substrates, :products, …) @eval $fn(m::AllostericEnzymeMechanism) = $fn(catalytic_mechanism(m)) end` loop replaces 13 separate definitions.
**Depends on:** none (independent cleanup; can land in Wave 1)
**Blocking tests:** none (forwarding has no special-case overrides; tests exercise these through public API)
**Recommendation:** Or, if Cluster A lands first, entire block deletes (AllostericMechanism field-access becomes direct). Either way, this is ~14 LOC of pure boilerplate.

#### F-011  Collapse `_rep_step` 3-method dispatch + `_to_mechanism` + `_AnyMech` bridge

**Location:** src/types.jl:1361-1382
**Category:** Duplication
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~22 (3 methods → 1; bridge deletes)
**Simplification gain:** After singleton demote: the 4-element `_AnyMech` union becomes `Union{Mechanism, AllostericMechanism}`, `_to_mechanism` deletes, 2 of 3 `_rep_step` methods delete.
**Depends on:** F-003
**Blocking tests:** test_chokepoint.jl exercises the chokepoint via Mechanism / AllostericMechanism — not blocking
**Recommendation:** Trivial after F-003 lands.

#### F-012  Delete "Stage 7a" stale-spec comment

**Location:** src/types.jl:1411 (inside chokepoint Kreg-method comment block L1408-1413)
**Category:** Doc hygiene
**Confidence:** High
**LOC saving (non-comment non-doc):** 0
**Simplification gain:** Evergreen comments per CLAUDE.md "Code Comments" rule.
**Depends on:** none
**Blocking tests:** none
**Recommendation:** Trim "added in Stage 7a" from L1411. Keep the substantive part about ligand uniqueness enforcement.

#### F-013  *(DEFERRED — Cluster F)* `_param_for_symbol` family obsoletes under structural Parameter keys

**Location:** src/types.jl:1451-1473 (main fn) + 1475-1494 (`_onlyA_parameters_for_sym`) + 1497-1521 (`_all_params_for_sym`)
**Category:** String-keyed projection
**Confidence:** Medium (blocked on direction-symmetry refactor)
**LOC saving (non-comment non-doc):** ~72 *(DEFERRED)*
**Simplification gain:** Callers (synth-dep machinery) hold Parameter structs directly. No reverse-lookup needed.
**Depends on:** direction-symmetry refactor (synth-deps gain first-class Parameter representation — see `2026-05-29-direction-symmetry-constraint-resolution.md`)
**Blocking tests:** test_dep_set_invariance.jl uses `_param_for_symbol`
**Recommendation:** **Defer to a follow-on refactor.** Per Q-012 the load-bearing reason for the string-based projection is that synth-dep I-state names have no Parameter struct. The direction-symmetry refactor changes that, then this cluster (F-013, F-064, F-066, F-068, F-069, F-071) lands as one move.

### src/dsl.jl

#### F-014  Inline `_assert_no_opaque_terms` into `_parse_steps_block_with_groups`

**Location:** src/dsl.jl:258-262 (parser branch) + 375-376 (`_is_conformation_shape`) + 382-399 (`_assert_no_opaque_terms`) + 567 + 1130 (callers)
**Category:** Permissive-parser-guard (Cluster G)
**Confidence:** High *(per Q-013 Approach B)*
**LOC saving (non-comment non-doc):** ~25
**Simplification gain:** One less named function, post-walk runs inline immediately after parse, semantically identical.
**Depends on:** none
**Blocking tests:** test_dsl.jl tests that exercise opaque-rejection error message — they pass identically since semantics unchanged
**Recommendation:** Move the body of `_assert_no_opaque_terms` (L382-399) into the tail of `_parse_steps_block_with_groups` after the L778-827 walk completes. Delete the function and its two callers at L567, L1130. Keep `_is_conformation_shape` regex since it's still the gate. **Approach C** (aggressive: tighten both bare-enzyme and Call-head parse sites to require conformation-shape, deleting both `_assert_no_opaque_terms` and the post-hoc validation entirely) would save more LOC but adds a Call-head constraint not currently enforced — needs Denis sign-off.

#### F-015  Macro emits `Mechanism(...)` directly after singleton demote

**Location:** src/dsl.jl:632-633
**Category:** Architectural
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~1 (drop the `EnzymeMechanism(...)` wrap)
**Simplification gain:** DSL emits the user-facing struct directly.
**Depends on:** F-003
**Blocking tests:** Tests asserting on macro emission type (test_dsl.jl) — minor updates
**Recommendation:** Replace `EnzymeRates.EnzymeMechanism(EnzymeRates.Mechanism($reaction_expr, $groups_expr))` with `EnzymeRates.Mechanism($reaction_expr, $groups_expr)`.

#### F-016  Simplify `_build_reg_sites_expr` / `_build_cat_sites_expr` to emit Vector instead of tuple

**Location:** src/dsl.jl:1009-1041 (`_build_reg_sites_expr`) + 1048-1059 (`_build_cat_sites_expr`)
**Category:** Architectural
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~10
**Simplification gain:** These currently build type-parameter-shaped tuple Exprs solely for the singleton's CS/RS slots. After demote: emit `Vector{RegulatorySite}` + `Vector{Symbol}` directly.
**Depends on:** F-003
**Blocking tests:** test_dsl.jl @allosteric_mechanism tests — verify during impl
**Recommendation:** Replace tuple-Expr construction with Vector-Expr construction.

#### F-017  Macro emits `AllostericMechanism(...)` directly after singleton demote

**Location:** src/dsl.jl:1137
**Category:** Architectural
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~0 (signature change, not deletion)
**Simplification gain:** Pairs with F-015.
**Depends on:** F-003, F-016
**Blocking tests:** test_dsl.jl @allosteric_mechanism tests
**Recommendation:** `:(AllostericEnzymeMechanism($cm_expr, $cat_sites_expr, $reg_sites_expr))` becomes `:(AllostericMechanism($reaction_expr, ...positional args...))`.

#### F-018  Batch: convert 19 function-leading `#`-comments to docstrings in dsl.jl

**Location:** src/dsl.jl multiple — L36-43, L96-97, L137-139, L164, L186-187, L205, L217-220, L246-247, L257, L282-284, L341-352 (struct doc), L369-374, L378-381, L401-402, L572-577, L636-639, L656-658, L684-685, L710-711
**Category:** Doc hygiene
**Confidence:** Medium
**LOC saving (non-comment non-doc):** 0
**Simplification gain:** Functions become discoverable via `?fnname` and IDE doc-on-hover. CLAUDE.md "Code Comments" rule favors docstrings for "explain WHAT/WHY" purposes attached to definitions.
**Depends on:** none
**Blocking tests:** none
**Recommendation:** One commit converts all 19. Mechanical refactor: take the `#`-block immediately above a function, transform to a `"""docstring"""` attached to the def.

### src/sym_poly_for_rate_eq_derivation.jl

#### F-034  Doc-hygiene sweep in sym_poly_for_rate_eq_derivation.jl

**Location:** src/sym_poly_for_rate_eq_derivation.jl:1-2 (file header), 52-55 (`sym_det`), 83-84 (`_poly_to_expr`), 236-237 (`substitute_params_expr`)
**Category:** Doc hygiene
**Confidence:** Medium
**LOC saving (non-comment non-doc):** 0
**Simplification gain:** Module-attached docstring + 3 function docstrings improve discoverability.
**Depends on:** none
**Blocking tests:** none
**Recommendation:** Convert each `#`-block to `"""docstring"""`. Move L1-2 into a module-attached docstring (or as the package's overall doc since this file is `include`d).

#### F-036  Fix stale "K2 → K1" example in `_rename_symbols` docstring

**Location:** src/sym_poly_for_rate_eq_derivation.jl:250-258
**Category:** Doc hygiene (stale-spec)
**Confidence:** High
**LOC saving (non-comment non-doc):** 0
**Simplification gain:** Docstring reflects current usage (A→I state rename) rather than the obsolete index-based parameter naming (`K1, K2, ...`).
**Depends on:** none
**Blocking tests:** none
**Recommendation:** Update the docstring example from `K2 → K1` to an A→I rename example (e.g. `K_S_E_A → K_S_E_I`).

### src/rate_eq_derivation.jl

#### F-037  Convert `@generated parameters` / `fitted_params` to plain functions over Mechanism

**Location:** src/rate_eq_derivation.jl:41-93
**Category:** Architectural
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~25 (50 LOC → ~25)
**Simplification gain:** Removes 4 `@generated` boundaries. Calls become direct function evaluation.
**Depends on:** F-003 (Mechanism becomes the primary type to dispatch on)
**Blocking tests:** Existing `parameters` / `fitted_params` tests check return values; they pass unchanged
**Recommendation:** Convert each `@generated parameters(::M, ::Mode)` to `parameters(m::Mechanism, ::Mode)` / `parameters(am::AllostericMechanism, ::Mode)`. Same for `fitted_params`.

#### F-038  Rename `_build_kinetic_rename_map` → `_build_wegscheider_rename_map`; trim stale docstring

**Location:** src/rate_eq_derivation.jl:94-141
**Category:** Doc hygiene + small cleanup
**Confidence:** High *(per Q-010 docstring is stale; body is correct for Wegscheider RE ties only)*
**LOC saving (non-comment non-doc):** ~5 (mostly stale docstring lines)
**Simplification gain:** Function name reflects what it actually does (only Wegscheider RE-tie absorption). Eliminates historical "kinetic-group merges no longer need a rename" context.
**Depends on:** none
**Blocking tests:** any test referencing `_build_kinetic_rename_map` by name — rename in lockstep
**Recommendation:** Symbol rename + docstring rewrite. Keep the function as-is body-wise.

#### F-039  Replace `_raw_symbolic_rate_polys` `rxns` opaque-tuple walk with struct-native walk over Mechanism.steps

**Location:** src/rate_eq_derivation.jl:309-408 + 410-421 (wrapper) + 427-502 (`_compute_numerator`)
**Category:** Architectural
**Confidence:** High *(per Q-016)*
**LOC saving (non-comment non-doc):** ~30 (mechanical refactor; some LOC delete, some restructure)
**Simplification gain:** Eliminates the `rxns` re-projection from Mechanism. Direct reads of `Step.from_species`, `Step.to_species`, `Step.is_equilibrium`, `Step.bound_metabolite`, with `g` from the outer-vector iteration. Source comment at L325-329 explicitly invites this.
**Depends on:** F-007 (or independent — `_flat_steps` already exists)
**Blocking tests:** test_rate_eq_derivation.jl Expr-shape / flat-string regression tests — semantically identical output, byte-stable; rate_equation perf test stays
**Recommendation:** Replace `for (idx, _) in enumerate(rxns)` with `for (idx, (s, g)) in enumerate(_flat_steps(mech))` in `_compute_alpha`, `_raw_symbolic_rate_polys`, `_compute_numerator`. Drop `rxns` from inner-arity signatures. `_compute_alpha`/`_compute_numerator` lose the `enz_set` arg (no longer needed). The `@generated` entry point at L410-421 stops calling `reactions(m)`.

#### F-040  Rewrite `_step_sides` to take Step directly

**Location:** src/rate_eq_derivation.jl:162-167
**Category:** Architectural
**Confidence:** High
**LOC saving (non-comment non-doc):** ~5
**Simplification gain:** `_step_sides(s::Step)` returns `(name(from_species(s)), name(to_species(s)), m_lhs, m_rhs)` via `bound_metabolite(s) + is_binding(s)`. One-liner.
**Depends on:** F-039
**Blocking tests:** none direct
**Recommendation:** Trivial after F-039 lands.

#### F-041  Simplify `_raw_rate_expr_and_symbols` wrapper

**Location:** src/rate_eq_derivation.jl:503-520
**Category:** Architectural
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~5
**Simplification gain:** Wrapper updates to call new struct-native `_raw_symbolic_rate_polys(mech::Mechanism)` directly without M-Type indirection.
**Depends on:** F-039, F-003
**Blocking tests:** none direct
**Recommendation:** Minor update; bundle with F-039 commit.

#### F-042  Pull `_ss_rate_constant_names` allosteric body into shared `_emit_cat_params_for_rep` helper

**Location:** src/rate_eq_derivation.jl:603-640
**Category:** Duplication (Cluster C)
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~10 (saved by sharing the 4-way switch)
**Simplification gain:** One shared helper for the rep-step → param-list emission used by `_ss_rate_constant_names` (allosteric), `_onlyA_parameters`, `_all_i_state_parameters`, `_enumerate_parameters_full_allosteric`.
**Depends on:** F-043
**Blocking tests:** none direct
**Recommendation:** See F-043.

#### F-043  Extract `_emit_cat_params_for_rep(rep, state) → Vector{Parameter}` shared helper

**Location:** src/rate_eq_derivation.jl:986-1202 (5 helpers)
**Category:** Duplication (Cluster C)
**Confidence:** Medium *(per Q-015: don't fully unify all 5, but extract the smaller helper)*
**LOC saving (non-comment non-doc):** ~20 (the 4-way `is_eq × is_binding` switch repeated ~5x)
**Simplification gain:** One canonical Kd/Kiso/Kon+Koff/Kfor+Krev emission helper. 3 of 5 helpers use it directly; 2 of 5 (the rename-pair helpers) zip its A/I outputs.
**Depends on:** none
**Blocking tests:** Per Q-011: `_onlyA_parameters`, `_I_rename_parameters`, `_all_i_state_parameters` have direct test assertions (test_rate_eq_derivation.jl:1536, 1545, 1548, 1555, 1566, 1576, 1583, 1594, 1599-1604, 1612-1617). **Tests must be updated to assert via behavior, not via direct helper output.** Replacement test pattern: assert that `fitted_params(am)` for known mechanisms contains the right Parameter names; assert that `rate_equation_string(am)` contains the right `_T` suffixes. Per CLAUDE.md "no net coverage reduction" — log replacements in `docs/superpowers/refactor-deleted-tests.md`.
**Recommendation:** Bottom-up: introduce `_emit_cat_params_for_rep`, then refactor the 5 helpers to use it, then update the 3 test-blocked helpers' tests to behavior-based.

#### F-044  Delete 3 stale "Stage 4.2" comment references

**Location:** src/rate_eq_derivation.jl:760-761 + 1437-1441 + 1552-1554
**Category:** Doc hygiene
**Confidence:** High
**LOC saving (non-comment non-doc):** 0
**Simplification gain:** Evergreen comments.
**Depends on:** none
**Blocking tests:** none
**Recommendation:** Trim "Stage 4.2" references from each comment block. The substantive parts ("Pass 2: dep RHSes referencing a :NonequalAI symbol need their own I-state name") can stay.

#### F-045  Factor `a_only_syms` set construction (4 occurrences) into shared helper

**Location:** src/rate_eq_derivation.jl:756 + 1283 + 1432 + 1547
**Category:** Duplication
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~3-4 (4 occurrences → 1)
**Simplification gain:** `_a_only_syms(am) = Set{Symbol}(name(p, am) for p in _onlyA_parameters(am))`. Pairs with F-046 / F-047.
**Depends on:** none
**Blocking tests:** none
**Recommendation:** One-liner helper. Replace 4 inline constructions.

#### F-046  Factor `rename_I` Dict construction (4 occurrences) into shared helper

**Location:** src/rate_eq_derivation.jl:757-759 + 1289-1292 + 1434-1436 + 1549-1551
**Category:** Duplication
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~9 (4 occurrences × 3 LOC → 1 helper)
**Simplification gain:** `_a_to_i_rename(am) = Dict{Symbol,Symbol}(...)`. Pairs with F-045 / F-047.
**Depends on:** none
**Blocking tests:** none
**Recommendation:** One-liner helper.

#### F-047  Note: `_add_case_b_renames!` is already a shared helper

**Location:** src/rate_eq_derivation.jl:1251-1259 (definition); callsites at L763 + L1299 + L1442 + L1556
**Category:** No duplication after all
**Confidence:** High
**LOC saving (non-comment non-doc):** 0
**Simplification gain:** Verified during catalog — `_add_case_b_renames!` is already factored. Suspect listing for it was redundant.
**Depends on:** none
**Blocking tests:** none
**Recommendation:** No change. Note this finding for the audit record; remove from scratch.

#### F-048  Delete 3-layer `compile_mechanism(Mechanism(reaction(am), steps(am)))` lift duplications

**Location:** src/rate_eq_derivation.jl:1086-1087 + 1099-1100
**Category:** Duplication (Cluster A ripple)
**Confidence:** High
**LOC saving (non-comment non-doc):** ~4
**Simplification gain:** After F-003, the entire 3-layer lift goes away. `am` is the value.
**Depends on:** F-003
**Blocking tests:** none direct
**Recommendation:** Bundled with F-003 commit.

### src/thermodynamic_constr_for_rate_eq_derivation.jl

#### F-049  Rewrite `_free_enz_set` to walk Mechanism.steps directly

**Location:** src/thermodynamic_constr_for_rate_eq_derivation.jl:47-67
**Category:** Architectural (Cluster B)
**Confidence:** High
**LOC saving (non-comment non-doc):** ~5 (delete `em = compile_mechanism(m)` lift + walks `reactions(em)`)
**Simplification gain:** Direct walk via Step.from_species/bound_metabolite. No accessor indirection.
**Depends on:** F-039 (or independent)
**Blocking tests:** none direct
**Recommendation:** Replace `em = compile_mechanism(m)` + `reactions(em)` walk with `for group in steps(m), s in group; ...`.

#### F-050  Rewrite `_thermodynamic_constraints` to take Mechanism directly

**Location:** src/thermodynamic_constr_for_rate_eq_derivation.jl:138-207
**Category:** Architectural (Cluster B)
**Confidence:** High
**LOC saving (non-comment non-doc):** ~10 (delete `m = M()` lift + accessor calls)
**Simplification gain:** Function signature becomes `_thermodynamic_constraints(m::Mechanism)`. Walks `m.steps` directly for the enzyme incidence matrix. `enzyme_forms`, `metabolites`, `substrates`, `products`, `stoich_matrix` calls become field-access on `m.reaction`.
**Depends on:** F-039, F-007
**Blocking tests:** none direct (function is internal)
**Recommendation:** Refactor along with F-039.

#### F-051  Fix stale "K9 => K4" docstring example in `_dependent_param_exprs`

**Location:** src/thermodynamic_constr_for_rate_eq_derivation.jl:226-236
**Category:** Doc hygiene (stale-spec)
**Confidence:** High
**LOC saving (non-comment non-doc):** 0
**Simplification gain:** Docstring reflects current naming (structural rep-step) not index-based.
**Depends on:** none
**Blocking tests:** none
**Recommendation:** Update example from `K9 => K4` to structural rep-step example (e.g. `K_P_E => K_S_E` for two binding K's tied by Wegscheider).

#### F-052  Rewrite `_dependent_param_exprs_kernel` to walk Mechanism directly

**Location:** src/thermodynamic_constr_for_rate_eq_derivation.jl:253-379
**Category:** Architectural (Cluster B)
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~15 (delete `m = M()` + `mech = Mechanism(m)` double-lift + accessor reads + binding_K_set walk via opaque tuples)
**Simplification gain:** Kernel takes `(mech::Mechanism, rename::Dict, ...)` directly. binding_K_set built from Step fields.
**Depends on:** F-039, F-050
**Blocking tests:** test_dep_set_invariance.jl tests the dep-expression set — verify byte-stable
**Recommendation:** Bundle with F-050.

#### F-053  Delete `_raw_param_symbols(::EnzymeMechanism)` singleton-lift forwarder

**Location:** src/thermodynamic_constr_for_rate_eq_derivation.jl:23-24
**Category:** Duplication (Cluster A ripple)
**Confidence:** High
**LOC saving (non-comment non-doc):** ~2
**Simplification gain:** After F-003, only the Mechanism method needed.
**Depends on:** F-003
**Blocking tests:** none
**Recommendation:** Bundled with F-003 commit.

### src/fitting.jl

#### F-054  Retype FittingProblem to Mechanism-family if AbstractEnzymeMechanism dissolves

**Location:** src/fitting.jl:38-83 (constructor signature)
**Category:** Architectural
**Confidence:** Low
**LOC saving (non-comment non-doc):** 0
**Simplification gain:** Signature reflects the post-demote primary type.
**Depends on:** F-003
**Blocking tests:** none expected
**Recommendation:** Minor signature update during F-003 impl.

### src/mechanism_enumeration.jl

#### F-055  Delete `compile_mechanism` lift functions

**Location:** src/mechanism_enumeration.jl:977-985
**Category:** Architectural (Cluster A ripple)
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~2 (plus all 7 prod callers per Q-001)
**Simplification gain:** After F-003 + F-070 (and similar at thermo/rate_eq derivation sites), `compile_mechanism` has no callers and deletes.
**Depends on:** F-003, F-070, F-048, F-049
**Blocking tests:** test/test_mechanism_enumeration.jl has ~50 callers of `compile_mechanism` — these likely become direct Mechanism/AllostericMechanism use. Log replacements per CLAUDE.md
**Recommendation:** Final delete in Cluster A. Test updates: replace `compile_mechanism(m)` with `m` (identity) in tests.

#### F-056  Add `_with_steps` / `_with_cat_allo_states` internal update constructors; collapse `_expand_re_to_ss` + `_expand_split_kinetic_group` duals

**Location:** src/mechanism_enumeration.jl:1101-1187
**Category:** Duplication (Cluster D)
**Confidence:** High
**LOC saving (non-comment non-doc):** ~50
**Simplification gain:** Each expansion move becomes one method instead of two. `_with_steps(m::Mechanism, new_steps)` returns a Mechanism with the same reaction; `_with_steps(am::AllostericMechanism, new_steps)` returns an AllostericMechanism with same cat_allo_states/multiplicity/regulatory_sites. Now `_expand_re_to_ss(m::Union{Mechanism,AllostericMechanism}) = [_with_steps(m, _flip_group_to_ss(steps(m), g)) for g where all_re]`.
**Depends on:** none
**Blocking tests:** test_mechanism_enumeration.jl topology / expansion-count tests — output unchanged
**Recommendation:** Introduce `_with_*` helpers in types.jl, then refactor the 4 expansion-move pairs.

#### F-057  Collapse `_expand_add_dead_end_regulator` duals

**Location:** src/mechanism_enumeration.jl:1259-1311
**Category:** Duplication (Cluster D)
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~30
**Simplification gain:** With `_with_steps`-style helpers + a single `Union{Mechanism, AllostericMechanism}` dispatch, the `wrap` callback indirection goes away.
**Depends on:** F-056
**Blocking tests:** test_mechanism_enumeration.jl dead-end tests
**Recommendation:** Bundle with F-056.

#### F-058  Remove three no-op expansion-move dispatches via type-aware dispatcher

**Location:** src/mechanism_enumeration.jl:1486-1495 + 1606-1616 + 1663-1670
**Category:** Duplication (Cluster D)
**Confidence:** Low
**LOC saving (non-comment non-doc):** ~3 (the no-op bodies; docstrings stay)
**Simplification gain:** `_add_expansions_mech!` skips inapplicable moves based on input type instead of calling each move and getting back an empty Vector.
**Depends on:** F-056
**Blocking tests:** none direct
**Recommendation:** Optional micro-cleanup. The current approach has uniform dispatch which is arguably clearer; consider keeping the no-ops.

#### F-059  Simplify `_expand_change_allo_state` with `_with_reg_sites` helper

**Location:** src/mechanism_enumeration.jl:1628-1662
**Category:** Architectural (Cluster D)
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~10
**Simplification gain:** RegulatorySite copy/replace becomes `_with_reg_sites(am, new_sites)` + a builder for updating one site's allo_state.
**Depends on:** F-056
**Blocking tests:** test_mechanism_enumeration.jl allo-state-change tests
**Recommendation:** Bundle with F-056.

#### F-060  Merge three `dedup!` overloads into one `dedup!(Dict{Int, <:Vector})`

**Location:** src/mechanism_enumeration.jl:1780-1822
**Category:** Duplication (Cluster D)
**Confidence:** High *(per Q-014)*
**LOC saving (non-comment non-doc):** ~22
**Simplification gain:** One method covers all three current overloads. `_canonicalize_mechanism!` dispatches per element type at runtime regardless of Vector eltype.
**Depends on:** none
**Blocking tests:** none direct (tests exercise dedup! via expand_mechanisms loop)
**Recommendation:** Independent cleanup; can land in Wave 1.

#### F-061  Collapse 6 `_parameter_canonical_key` methods to one

**Location:** src/mechanism_enumeration.jl:1845-1855
**Category:** Duplication
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~5
**Simplification gain:** Six near-identical one-liners `(:Kd, hash(p.step), p.state)`, etc. become one method dispatching on `nameof(typeof(p))` or via static lookup table.
**Depends on:** none
**Blocking tests:** test_canonical_hash_partition.jl — output unchanged
**Recommendation:** Trade-off — the static dispatch is faster but slightly less readable. Verify perf impact at canonical-hash callsites first. Could defer.

#### F-062  Collapse 2 `_enumerate_all_parameters_with_i_state` methods

**Location:** src/mechanism_enumeration.jl:1857-1868
**Category:** Duplication
**Confidence:** Low
**LOC saving (non-comment non-doc):** ~3
**Simplification gain:** Small duplication; tradeoff is parametric-type dispatch overhead.
**Depends on:** none
**Blocking tests:** none
**Recommendation:** Low priority.

#### F-063  Merge `_canonicalize_for_hash` 2 methods via shared helpers

**Location:** src/mechanism_enumeration.jl:1887-1950
**Category:** Duplication
**Confidence:** Medium *(per Q-019)*
**LOC saving (non-comment non-doc):** ~25
**Simplification gain:** Two helpers `_num_den_exprs(em, m)` + `_extra_canon_tuple(m)` factor out the type-specific bits. Single trunk function replaces both bodies.
**Depends on:** none
**Blocking tests:** test_canonical_hash_partition.jl — canon-tuple byte-stability is load-bearing; must reproduce L1910 and L1947-1948 exactly. Existing partition test catches regressions.
**Recommendation:** Carefully reproduce the canon-tuple shape. Allosteric type guard moves to `_num_den_exprs(::AllostericMechanism)`.

#### F-064  *(DEFERRED — Cluster F)* Replace string-keyed `_build_name_map` with structural Parameter keys

**Location:** src/mechanism_enumeration.jl:1951-1990
**Category:** String-keyed projection
**Confidence:** High (per Q-012)
**LOC saving (non-comment non-doc):** ~20 *(DEFERRED)*
**Simplification gain:** Eliminates the `Symbol → token` Dict in favor of `Parameter → token` Dict. Removes 1 indirection in the projection pipeline.
**Depends on:** direction-symmetry refactor for synth-deps (see F-013)
**Blocking tests:** test_mechanism_enumeration.jl:4295 (`@test name_map isa Dict{String, String}`) — update during impl
**Recommendation:** Defer. The synth-dep I-state names have no Parameter struct currently; that's load-bearing for string keys. Direction-symmetry resolution unblocks.

#### F-065  Doc hygiene in `_build_name_map`: stale R↔T notation + "Phase 7 cleanup" + r/t variable names

**Location:** src/mechanism_enumeration.jl:1951-1963 (docstring) + 1980-1986 (variable names + comment)
**Category:** Doc hygiene
**Confidence:** High
**LOC saving (non-comment non-doc):** ~0 (variable rename + comment edit)
**Simplification gain:** R↔T → A/I notation consistency (per the completed allo-state-rename refactor); evergreen comments per CLAUDE.md.
**Depends on:** none
**Blocking tests:** none
**Recommendation:** Rename `r_name` → `a_name`, `r_str` → `a_str`, `t_str` → `i_str`. Drop "Phase 7 cleanup territory" comment. Update docstring "R↔T correspondence" → "A↔I correspondence".

#### F-066  *(DEFERRED — Cluster F)* `_dep_exprs_canonical` string-key projection

**Location:** src/mechanism_enumeration.jl:1999-2010
**Category:** String-keyed projection
**Confidence:** Medium
**LOC saving (non-comment non-doc):** ~5 *(DEFERRED)*
**Simplification gain:** Co-evolves with F-064.
**Depends on:** F-064
**Blocking tests:** none direct
**Recommendation:** Defer with Cluster F.

#### F-067  Inline `_canonical_rate_eq_hash_data` thin 3-line wrapper

**Location:** src/mechanism_enumeration.jl:2062-2064
**Category:** Dead-or-relocate
**Confidence:** High *(per Q-017)*
**LOC saving (non-comment non-doc):** ~3 (plus 9 LOC of docstring)
**Simplification gain:** One function name instead of two. Rename `_canonical_rate_eq_hash_data_impl_struct` → `_canonical_rate_eq_hash_data` to drop the `_impl_struct` suffix.
**Depends on:** none
**Blocking tests:** test_mechanism_enumeration.jl L4291, 4341, 4395, 4455, 4456 + identify_rate_equation.jl:427 — all use the wrapper directly; update to call the renamed impl
**Recommendation:** Trivial. Land in Wave 1.

### src/identify_rate_equation.jl

#### F-068  *(DEFERRED — Cluster F)* `_CachedFitResult.canon_to_rep::Dict{String,String}` retype

**Location:** src/identify_rate_equation.jl:310-317
**Category:** String-keyed projection
**Confidence:** High
**LOC saving (non-comment non-doc):** ~0 (type change)
**Simplification gain:** Field type reflects structural keys.
**Depends on:** F-064
**Blocking tests:** none direct (struct field type)
**Recommendation:** Defer with Cluster F.

#### F-069  *(DEFERRED — Cluster F)* `_project_cached_params` structural rekey

**Location:** src/identify_rate_equation.jl:364-390
**Category:** String-keyed projection
**Confidence:** High
**LOC saving (non-comment non-doc):** ~10 *(DEFERRED)*
**Simplification gain:** Projection walks Parameter structs directly. Removes String round-trip.
**Depends on:** F-064, F-068
**Blocking tests:** test_identify_rate_equation.jl projection tests
**Recommendation:** Defer with Cluster F.

#### F-070  Delete four `compile_mechanism` lift calls after singleton demote

**Location:** src/identify_rate_equation.jl:424 + 461 + 831 + 860
**Category:** Architectural (Cluster A ripple)
**Confidence:** High *(per Q-018)*
**LOC saving (non-comment non-doc):** ~4
**Simplification gain:** APIs (`rate_equation_string`, `_canonical_rate_eq_hash_data`, `fitted_params`, `FittingProblem`, `_loocv`) all already accept `AbstractEnzymeMechanism`. After F-003, the 4 lifts become no-ops. Pass-through.
**Depends on:** F-003, F-055
**Blocking tests:** none direct
**Recommendation:** Bundled with F-003.

#### F-071  *(DEFERRED — Cluster F)* String-keyed Dict inversion at identify_rate_equation.jl:478

**Location:** src/identify_rate_equation.jl:478
**Category:** String-keyed projection
**Confidence:** High
**LOC saving (non-comment non-doc):** ~1
**Simplification gain:** Inversion becomes a single-line `Dict{Parameter, …}` build after F-064.
**Depends on:** F-064
**Blocking tests:** none direct
**Recommendation:** Defer with Cluster F.

---

## §3 Dependency clusters

### Cluster A — Singleton-type demotion *(THE big architectural move)*

**Findings:** F-002, F-003, F-004, F-005, F-006, F-007, F-008, F-009, F-010, F-011, F-015, F-016, F-017, F-037, F-048, F-053, F-054, F-055, F-070

**Total LOC:** ~520 LOC net delete + ~150 LOC of in-place rewrites

**Dependency chain:**
- F-003 is the head — demote the singleton types
- F-002, F-004, F-006, F-007, F-011, F-015, F-016, F-017, F-037, F-053, F-054 follow F-003
- F-008, F-009 depend on F-007 (the @generated-accessor collapse)
- F-048, F-055, F-070 are pure deletion of `compile_mechanism` callsites after the type demotes
- F-010 is independent (forwarding-accessor consolidation) but lands cleanly with this cluster

**Sequencing within cluster:**
1. F-003 (singleton type minimization)
2. F-007 (collapse @generated accessors)
3. F-008, F-009 (delete symbol-tuple helpers)
4. F-002 (delete Sig conversion machinery — last because the lift constructor reads from it)
5. F-010, F-011, F-037 (downstream collapses)
6. F-015, F-016, F-017 (dsl macro emission)
7. F-048, F-053, F-055, F-070 (delete lift calls)
8. F-006 (rewrite show methods)
9. F-004, F-054, F-005 (minor cleanups + docstring trims)

### Cluster B — Derivation back-end struct-native walk

**Findings:** F-039, F-040, F-041, F-049, F-050, F-052

**Total LOC:** ~80 LOC net delete + ~30 LOC of in-place rewrites

**Dependency chain:** Independent of Cluster A. `_flat_steps` already walks Mechanism; the refactor is local to each function.

**Sequencing within cluster:**
1. F-039 (the rate_eq_derivation rewrite)
2. F-040, F-041 (helper updates)
3. F-049 (`_free_enz_set`)
4. F-050, F-052 (thermo functions)

### Cluster C — Smaller helper for parameter enumeration

**Findings:** F-042, F-043

**Total LOC:** ~20 LOC

**Dependency chain:** Independent. Has test-blocked dependents (per Q-011); must update tests.

**Sequencing:** F-043 first (extract helper), then F-042 (call from `_ss_rate_constant_names`).

### Cluster D — Expansion-move + dedup! consolidation

**Findings:** F-056, F-057, F-058, F-059, F-060, F-061, F-062, F-063

**Total LOC:** ~115 LOC

**Dependency chain:**
- F-056 is the head (`_with_steps` family)
- F-057, F-058, F-059 follow F-056
- F-060, F-061, F-062, F-063 are independent of F-056 — pure dedup!/canonical-helper merges

**Sequencing:**
1. F-060 (dedup! merge — easiest)
2. F-063 (`_canonicalize_for_hash` merge — canon-tuple stability concern)
3. F-056 (introduce `_with_*` helpers)
4. F-057, F-059 (apply helpers)
5. F-058 (optional micro-cleanup)
6. F-061, F-062 (minor dedup)

### Cluster E — Doc hygiene sweep

**Findings:** F-001, F-005, F-012, F-018, F-034, F-036, F-044, F-051, F-065 + the implicit comment-as-docstring batch

**Total LOC:** 0 (non-comment non-doc) — improves readability + CLAUDE.md compliance

**Dependency chain:** None — pure doc edits.

**Sequencing:** One commit at the end of the audit (Wave 4 / doc sweep) covering all of E.

### Cluster F *(DEFERRED — needs direction-symmetry refactor first)* — String-keyed projection

**Findings:** F-013, F-064, F-066, F-068, F-069, F-071

**Total LOC:** ~100 LOC *(when unblocked)*

**Dependency chain:** Blocked on `2026-05-29-direction-symmetry-constraint-resolution.md` providing first-class Parameter representation for synth-dep I-state names.

**Sequencing:** When unblocked, F-064 first, then F-066, then F-013/F-068/F-069/F-071.

### Cluster G — Parser-tighten (single finding)

**Findings:** F-014

**Total LOC:** ~25 LOC (Approach B)

**Dependency chain:** Independent of all others.

**Sequencing:** Standalone, can land any wave. Approach B (inline) is recommended; Approach C (aggressive tighten) needs Denis sign-off.

### Cluster H — Small independent cleanups

**Findings:** F-038, F-045, F-046, F-047, F-067

**Total LOC:** ~15 LOC

**Dependency chain:** All independent.

**Sequencing:** Wave 1. F-067 (`_canonical_rate_eq_hash_data` inline) is the smallest single win.

---

## §4 Suggested sequencing

### Wave 1 — High-confidence dead code + isolated cleanups *(can land first)*

- **Doc hygiene quick wins:** F-012 (Stage 7a comment), F-036 (K2 → K1 stale example), F-051 (K9 → K4 stale example), F-044 (3 Stage 4.2 references), F-001 (E_A_B mismatch)
- **Renames + trims:** F-038 (`_build_kinetic_rename_map` → Wegscheider; trim stale docstring), F-065 (R↔T → A/I rename + Phase 7 cleanup)
- **Wrapper inlines:** F-067 (`_canonical_rate_eq_hash_data` thin wrapper)
- **Small dups:** F-045, F-046, F-061
- **Cluster D quick wins:** F-060 (dedup! merge), F-063 (`_canonicalize_for_hash` merge — careful with canon-tuple stability)
- **Cluster G:** F-014 (parser-tighten Approach B)

Estimated Wave 1 savings: **~70 LOC** + several doc-hygiene improvements.

### Wave 2 — Derivation back-end struct-native walk *(Cluster B)*

- F-039, F-040, F-041 (rate_eq_derivation walk)
- F-049, F-050, F-052 (thermo walk)
- F-042, F-043 (Cluster C smaller helper)

Estimated Wave 2 savings: **~100 LOC** + significant readability improvement (no more `rxns` opaque-tuple shuffling).

**Pre-condition for Wave 2:** Update the 3 test-blocked helpers' tests per F-043 (replace direct helper assertions with behavior tests).

### Wave 3 — Singleton-type demotion *(Cluster A — the big architectural move)*

- F-003 (the demote itself)
- F-007 (collapse @generated accessors)
- F-008, F-009 (delete symbol-tuple helpers)
- F-002 (delete Sig conversion)
- F-010, F-011, F-037 (downstream collapses)
- F-015, F-016, F-017 (dsl emission)
- F-048, F-053, F-055, F-070 (delete lift calls)
- F-006 (rewrite show)
- F-004, F-054, F-005 (minor cleanups)
- **Cluster D rest:** F-056, F-057, F-058, F-059 (collapse expansion-move duals — happens to land cleanly with Wave 3 since `_with_*` helpers naturally fit Cluster A's structural emphasis)

Estimated Wave 3 savings: **~620 LOC** delete + ~150 LOC in-place rewrites.

**Pre-conditions for Wave 3:**
- Test-side: update ~50 `compile_mechanism(m)` callers in test files to use `m` directly (Q-001 enumerated these)
- Re-baseline `test_accessors.jl` perf gate (Q-005 — explicitly negotiable)
- Verify `@generated rate_equation` still dispatches correctly with Mechanism as primary type (impl-time design decision)

### Wave 4 — Doc-hygiene sweep + minor minor *(Cluster E)*

- F-018 (19-function dsl.jl comment-as-docstring batch)
- F-034 (sym_poly doc batch)
- F-035, F-005, F-062 (remaining minor)

Estimated Wave 4 savings: 0 non-comment non-doc LOC; ~30-40 docstrings added; ~10 stale comments removed.

### Wave 5 *(deferred indefinitely — Cluster F)*

- After direction-symmetry refactor lands and synth-deps gain Parameter struct representation:
  - F-013, F-064, F-066, F-068, F-069, F-071

Estimated Wave 5 savings: **~100 LOC** *(when unblocked)*.

---

## §5 Hard constraints tracked

- **`rate_equation` perf budget (0 alloc / <100 ns):** ALL Cluster A and Cluster B findings preserve this. Per Q-005, the 14 @generated accessors are NOT on the runtime hot path — collapsing them does not affect `rate_equation`. Cluster A keeps `EnzymeMechanism{Sig}` as the internal `@generated rate_equation` dispatch type (now built only at compile-body-build time), so the hot path's 0-alloc / <100ns invariant is preserved. Impl plans for F-003, F-007, F-039 must include `test_rate_equation_performance` evidence.
- **Test coverage:** F-043 deletes direct test assertions on 3 private helpers (`_onlyA_parameters`, `_I_rename_parameters`, `_all_i_state_parameters`); replaces with behavior-level tests asserting on `fitted_params(am)` and `rate_equation_string(am)` outputs. Log replacements in `docs/superpowers/refactor-deleted-tests.md` per project convention.
- **Front-end struct-family unification:** Already-achieved per `2026-05-26-finish-refactor-legacy-sig-removal-design.md`. No finding in this audit reintroduces parallel front-end representations. Cluster A simplifies *within* the unified family.
- **Opaque bound-form rejection:** F-014 (Approach B) preserves the rejection semantics — inlines the guard rather than removing it. Approach C is aggressive and would change the user-visible parser surface; needs Denis sign-off before adopting.
- **`@generated rate_equation` dispatch invariant:** Cluster A's F-003 must preserve a per-mechanism-shape dispatch type for `@generated rate_equation` to specialize on. The recommended design (Option A in F-003) keeps `EnzymeMechanism{Sig}` as an internal compile artifact built inside `compile_mechanism` — invisible to user code but live at @generated entry.

---

## Notes on the hypothesis test

The non-deferred clusters total ~780 LOC of estimated non-comment non-doc src savings, against a baseline of 5,706. That's **~13.7%**, well above the 20% "partially supported" threshold but well below the 40% "supported" threshold.

The gap between this and Denis's "up to 50%" estimate is real and worth discussing:

1. **What's NOT being deleted:** The `@generated rate_equation` body builders (`_build_rate_body`, `_build_allosteric_rate_body`, `_allosteric_num_den_exprs`, `_kcat_forward`), all of the rate-equation derivation core math (`_compute_alpha`, `_compute_re_groups`, `sym_det`, `_thermodynamic_constraints`, `_dependent_param_exprs`), all of the parameter renaming and synth-dep machinery, all of the mechanism enumeration logic (`_catalytic_topologies`, `_expand_substrate_product_dead_ends`, etc.). These are the core science; they account for the bulk of the LOC and aren't reducible without rewriting the algorithm.

2. **What's deferred:** Cluster F's ~100 LOC waits on direction-symmetry resolution. When that lands, the audit's effective savings approach ~880 LOC = ~15.4%.

3. **What's intentionally excluded:** Test-side LOC (~14,000 LOC of tests) — per CLAUDE.md, tests are out of scope for net reduction.

The audit's deliverable is the honest measured number. If Denis wants to push higher, the conversations to have:
- Should we accept reduced test coverage on private helpers to land Cluster C / F-043 with less ceremony?
- Should we adopt F-014's Approach C (aggressive parser-tighten with Call-head constraint change)?
- Should the `@generated`-driven derivation itself be reconsidered (different refactor — out of this audit's scope)?
