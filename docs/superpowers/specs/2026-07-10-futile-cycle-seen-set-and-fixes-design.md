# Beam termination via a structural seen-set; idempotent canonicalization and RE→SS enumeration fix

## Context

The LDH search with four competitive inhibitors
(`docs/ldh_hpc_results/2026_07_09_results_2`) finishes the useful search by ~iteration 7, then
never terminates: it re-emits a fixed ~5000-equation set indefinitely. The evidence is direct —
the iteration CSV file sizes form a **period-6 limit cycle** from iteration 20 on
(iters 20–25 sizes `{22597665, 25872302, 24941089, 23496220, 24973675, 25790549}` repeat
byte-for-byte at 26–31, 32–37, …; iter 20 and 26 both have 5059 rows). The equation-complexity
filter (#63) stopped the earlier segfault; this is the remaining non-termination.

**Root cause (verified in code).** The beam's `frontier` (`identify_rate_equation.jl:613–632`)
holds *every structurally-distinct mechanism with no global seen-set*. The `memo`
(`_process_batch`) deduplicates only *fitting* (a repeated `eq_hash` returns "inherited"); an
inherited mechanism still joins the `frontier` and is expanded again. `expanded_by_count` is a
*count*, not a set. So any expansion move that produces a same-parameter-count,
structurally-distinct child feeds that child back into the frontier, and it is expanded forever.

Three move classes produce such Δ≤0 children (all reproduced firsthand under `c201d50`):

- **A — non-idempotent `_canonical_mechanism`.** A split no-op canonicalizes to 9 groups on one
  pass but 8 (the parent) on a second pass, so `_expand_split_kinetic_group`'s self-loop guard
  `child == _canonical_mechanism(parent)` misses it and keeps a non-canonical mechanism. This is
  the only move that puts a *non-canonical* mechanism into the frontier, breaking the invariant the
  whole "monotone drain → terminates" argument rests on.
- **B2 — over-parameterized split duplicate.** A split of an over-parameterized parent yields a
  child with a different `eq_hash` but the **identical** rate function
  (`max rel |v_parent − v_child| = 6.9e-16` over 3000 points); the textual `eq_hash` cannot see the
  equivalence. Root cause is an enumeration/canonicalization pair, detailed under the
  `_expand_re_to_ss` fix and Defect 2.
- **C — change_allo Δ0.** `_expand_change_allo_state` relaxes a tag toward `:NonequalAI` and can
  free no identifiable parameter, producing a same-count, genuinely-different child. This move is
  monotone (tags only move toward `:NonequalAI`), so it cannot sustain a cycle alone; it adds
  same-count volume that the split classes recirculate.

The fix that guarantees termination is a **structural seen-set**: a mechanism, once produced, is
expanded at most once. Termination then holds regardless of whether we have found every Δ≤0 edge.
The two confirmed correctness bugs behind A and B2 are fixed separately, on their own merits.

## Goals

1. **Guarantee termination** of the beam search at a bounded parameter count.
2. **Fix the two confirmed correctness bugs**: non-idempotent canonicalization (A) and
   `_expand_re_to_ss` emitting type-inconsistent, never-identifiable mechanisms (Defect 1 under B2).
3. **Keep Δ≤0 expansion edges visible** as an ongoing bug detector.

## Non-goals

- The near-degenerate loss-landscape / beam-width flood (the 200→5000 fan-out). That is a separate
  problem calling for a structural-identifiability prune, out of scope here.
- **The change_allo Δ0 prune** — deferred. Measured volume is small (≈2–7% of children at the
  parameter cap) and it drops genuinely-different, often fully-identifiable allosteric variants. The
  seen-set removes its termination relevance. Revisit during a dedicated Δ≤0-edge audit.
- **Defect 2 (multi-symbol Wegscheider merge)** — deferred (see the `_expand_re_to_ss` fix). Fixing
  Defect 1 resolves the B2 instance; independent reachability is unverified.

## The structural seen-set — primary termination mechanism

Give the beam a persistent set of the structural identities it has already produced, and never
process the same structure twice. Filter it **in `_process_batch`, alongside the existing
`max_params` and `complexity` filters**, so all three "skipped, and why" reasons are counted and
reported in one place.

**Design**

- `seen::Set{UInt64}` of `hash(mech)`. `Mechanism`/`AllostericMechanism` `==`/`hash` are structural
  (`types.jl:571`, `:664`) and construction is canonical, so the hash is a sound structural key.
  `UInt64` matches the `memo`'s existing eq_hash keying, stays light, and does not pin mechanism
  objects in memory.
- Thread `seen` through `_beam_search` like `memo` (one instance for the whole search, seeded empty;
  the base tier populates it on its first `_process_batch` call).
- In `_process_batch`, **master-side, before the parallel PASS-1 compile**: partition the input into
  fresh (`hash(m) ∉ seen`) and already-seen; add every fresh structure to `seen`; run PASS-1/PASS-2
  on the fresh set only; return `n_seen_skip` alongside `n_param_skip` and `n_complexity_skip`. The
  check is master-side (not inside the `pmap`) because `seen` is shared search history; broadcasting
  it to workers each batch would be wasteful, and the check is a cheap hash lookup.
- A structure is added to `seen` on **first production regardless of outcome** (fit, cap-skip,
  complexity-skip, or compile error). Repeats are then cheap pre-compile seen-skips — no wasted
  derivation, which the profiling showed is the dominant cost.
- `_batch_summary` gains a sixth bucket, `skipped (already seen)`, so the buckets still partition the
  child count: `new + inherited + already-seen + >params + >complexity + errored = total`.

**Why this terminates.** At a bounded `max_param_count` the set of structurally-distinct mechanisms
is finite. Each enters `_process_batch` once, is ingested into the frontier at most once, and is
expanded at most once. The frontier therefore drains. No dependence on having enumerated every Δ≤0
edge.

**Why the selected model is unchanged.** Expansion is deterministic, so re-expanding a structure
reproduces the same children. The reachable-mechanism set is identical with or without the seen-set;
the seen-set only removes repeated work. Every distinct rate function is still fit once (the first
time its structure appears), so cross-validation selection is unaffected. The one behavioral change
is that a borderline structure now gets a single (its earliest, best-budgeted) beam-selection
opportunity instead of a fresh one every sweep — which is the correct, well-defined behavior for a
terminating search.

**Δ≤0 detection is preserved.** First-production child rows already carry `parent_n_params`
(post-#64), so a Δ≤0 edge is any CSV row with `n_params ≤ parent_n_params`. No extra instrumentation;
this is the ongoing bug-detector for Defect 2, the change_allo prune, and any future Δ≤0 class.

**Tests**

- A constructed pair of parents whose expansions overlap: the second parent's repeat children are
  reported as `already seen` and do not re-enter the frontier.
- A bounded local LDH four-inhibitor run terminates by draining the frontier (today it does not).
- The `skipped (already seen)` count appears in `progress.log` and the buckets sum to the child count.
- An existing end-to-end `identify_rate_equation` selection test is unchanged (same selected model).

## Idempotent `_canonical_mechanism` — correctness (class A)

`_canonical_mechanism` must reach a fixed point. `_merge_tied_kinetic_groups` is single-pass and not
idempotent: merging one tie can expose a second tie (the reproducer's split no-op merges 9→9 on pass
one, 9→8 on pass two). Iterate the merge inside `_canonical_mechanism` (both the `Mechanism` and
`AllostericMechanism` methods) until the partition stops changing, with a small max-iteration safety
bound (convergence is two passes in every case observed). Keep `_merge_tied_kinetic_groups`
single-pass and testable; iterate at the caller so every caller benefits.

This restores the guard `child == _canonical_mechanism(parent)` for the split no-ops it currently
misses (class A), and it collapses renaming variants more broadly. With the seen-set already
guaranteeing termination, this is a correctness/quality fix, not a termination requirement.

**Re-baseline.** Some mechanisms' canonical form changes, so allosteric golden fixtures and eq_hashes
need a reviewed re-baseline.

**Tests**

- `_canonical_mechanism(_canonical_mechanism(m)) == _canonical_mechanism(m)` for every
  `MECHANISM_TEST_SPECS` entry.
- The reproducer A split no-op is dropped by `_expand_split_kinetic_group`.

## `_expand_re_to_ss` RE-only intent — correctness (Defect 1 under B2)

`_expand_add_dead_end_regulator` creates competitive-inhibitor bindings as rapid-equilibrium
(`mechanism_enumeration.jl:1646`, `is_equilibrium=true`) and builds each catalytic mirror step the
same type. `_expand_re_to_ss` (`:1207`) then flips *any* all-RE group to SS with no guard, so it
emits mechanisms the enumeration never intended. Two invariants restore the intent:

1. **Competitive-inhibitor bindings are RE-only.** A binding onto a dead-end inhibitor complex has
   only its dissociation constant identifiable; its speed is never identifiable, and it is always a
   dead end. `_expand_re_to_ss` must skip any group whose binding step binds a `CompetitiveInhibitor`.
2. **A catalytic step and its inhibitor-bound mirror share RE/SS type.** The same-group case is
   already correct (the whole group flips at once). When a split has separated the mirror into its
   own group, `_expand_re_to_ss` must flip the mirror group(s) together with the base — never one
   without the other. Find the mirror by reusing `_expand_add_dead_end_regulator`'s
   `de_species_map` machinery (`:1636–1660`): map a species to its inhibitor-added / -removed
   counterpart and find the step between them.

**Design decision (Denis).** The type-lock is on **step type**, not identifiability. Biochemically
the binding mechanism (RE vs SS) is independent of whether an inhibitor is bound elsewhere, so a base
step legitimately flipped to SS carries its dead-end mirror to SS too. The resulting "both-SS"
mechanism may not be fully identifiable, but it is parameter-dominated and is not a Δ0 cycle edge, so
it is accepted. Do **not** force the mirror to stay RE — that would fuse inhibitor-bound and
inhibitor-free forms into one RE group, which is biochemically wrong.

**Why this also fixes B2's Δ0 duplicate (Defect 2 is downstream here).** The B2 parent (`eq_hash
7f09e180b56580aa`, 10 fitted params) is over-parameterized with identifiable **rank 7, deficiency
3** — verified by complex-step sensitivity SVD (clean cliff; see appendix). The three null directions
are the two dead-end-SS binding *speeds* (`kon/koff_A_Pyruvateinh_E`,
`kon/koff_A_Pyruvate_EPyruvateinh`) and `L` (a literal zero column, absent from `v`, entangled with
those speeds). Those dead-end SS bindings exist only because `_expand_re_to_ss` flipped RE inhibitor
bindings to SS. Their `koff/kon` ratio is also what makes the A-state Wegscheider box tie
**multi-symbol** (`K_A_Pyruvate_E` and `K_A_Pyruvate_ENAD` both resolve to that ratio), which the
single-symbol-only kinetic-group merge cannot fold (Defect 2) — leaving the two Pyruvate groups
separate in the parent while a perturbed split fuses them in the child, the Δ0 duplicate. Converting
the dead-end SS bindings back to RE (this fix) removes the `koff/kon` ratio: the tie collapses to the
single-symbol `K_A_Pyruvate_E = K_A_Pyruvate_ENAD` the merge already folds, and the parent becomes
fully identifiable (**np 10→8, rank 8**, verified). So the RE-only fix resolves B2 end to end.

**Re-baseline.** The move emits different mechanisms; allosteric golden fixtures need a reviewed
re-baseline.

**Tests**

- `_expand_re_to_ss` never yields a steady-state inhibitor binding.
- After any `re_to_ss` flip, a catalytic step and its inhibitor-bound mirror share type (no
  split-separated divergence).
- The B2 reproducer parent, with its dead-end SS bindings RE, is fully identifiable (`np 10→8`,
  rank 8) and its split no longer produces a distinct-text duplicate.
- The equilibrium-flux oracle (`v = 0` at `Q = Keq`) still holds.

## Deferred, with rationale

- **Defect 2 — multi-symbol Wegscheider merge.** `_merge_tied_kinetic_groups` /
  `_state_wegscheider_rename_map` fold only *single-symbol* binding-K ties (`rhs isa Symbol`,
  `rate_eq_derivation.jl:1302–1304`); documented tech debt (`2026-07-09-wegscheider-rename-phantom-filter-brittleness.md`).
  In B2 the unfolded tie is in the **A-state** (verified: the A-state rename map is empty, the I-state
  map folds `K_I_Pyruvate_E → K_I_Pyruvate_ENAD` correctly), and the multi-symbol form is caused by the
  dead-end SS binding, so the RE-only fix resolves it. Whether a multi-symbol tie can arise **without** a dead-end
  SS binding is unverified; if it can, the merge needs a multi-symbol fix (the brittleness spec's
  "apply the rename when building the column set" is the clean lever) or the A-state pivot should prefer
  single-symbol binding-K ties. Deferred to a Δ≤0 audit; the seen-set backstops any residual duplicate.
- **The change_allo Δ0 prune.** Deferred; see Non-goals.

## Verification (whole change)

- All three reproductions still hold before the fixes and are resolved after: class A dropped by the
  idempotent guard; B2 parent fully identifiable and no split duplicate.
- A bounded local LDH four-inhibitor run terminates by draining the frontier rather than hitting
  wall-clock.
- Full `Pkg.test()` green, including the allosteric golden re-baseline, the `rate_equation`
  0-allocation / sub-120 ns performance gate, and the parameter-naming chokepoint guard.

## Sequencing

Three independent commits. **Seen-set first** (unblocks the search; behavior-preserving for the
selected model, no golden re-baseline). Then **idempotent canonicalization** and **the
`_expand_re_to_ss` invariants**, each as a correctness commit with its own reviewed golden re-baseline.

## Appendix — reproduced facts (firsthand, `c201d50`)

- **Cycle:** period-6 CSV file-size cycle from iteration 20; iter 20 and 26 both 5059 rows.
- **A:** `_canonical_mechanism(raw)` = 9 groups ≠ parent's 8 (guard keeps it);
  `_canonical_mechanism(_canonical_mechanism(raw))` = 8 = parent.
- **B2 function identity:** split child `eq_hash` ≠ parent yet `max rel |v_parent − v_child| = 6.9e-16`
  over 3000 wide-range points.
- **B2 identifiability:** complex-step ∂v/∂logθ SVD (no column-normalization — column-normalizing
  rescues near-zero columns and inflates rank), parent np=10 rank 7 deficiency 3; null directions =
  `{kon/koff_A_Pyruvateinh_E}`, `{kon/koff_A_Pyruvate_EPyruvateinh}`, `L` (zero column). Dead-end
  SS→RE: np=8, rank 8, deficiency 0.
- **B2 merge:** `_state_wegscheider_rename_map(:A)` = `Dict()` (multi-symbol tie skipped);
  `(:I)` = `Dict(K_I_Pyruvate_E → K_I_Pyruvate_ENAD)`; `_merge_tied_kinetic_groups` 9 groups → 9 (no
  merge), blocked on the A-state key component.
- **C:** `A_parent` (np=11, deficiency 0) yields a Δparams=0 change_allo child
  (`:OnlyA → :NonequalAI`) with a different `eq_hash` and `max rel |v_p − v_c| = 1.0` (genuinely
  different function).
- **Volume (120 cycling parents, pre-dedup):** re_to_ss 24%, split 27%, dead_end 23%, change_allo
  26%. Δ0 fraction (bounded compile sample, np=13 cap): change_allo ≈7%, split ≈0%.

Reproducers: `docs/superpowers/specs/2026-07-10-futile-cycle-reproducers.jl` (A and B2).
