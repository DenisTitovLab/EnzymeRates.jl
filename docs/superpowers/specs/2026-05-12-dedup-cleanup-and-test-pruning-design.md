# Dedup cleanup and test pruning — design

Date: 2026-05-12
Branch: `improve-deduplication-of-ident-equations`
Status: design (awaiting plan)

## Problem

Deduplication of equivalent rate equations is currently working correctly
on this branch — the recent commits fixed the canonicalizer, Pass 1 + Pass 2
absorption chain, Wegscheider section emission, and the LDH 4-mechanism
cluster collapse. What remains is code and test cleanup:

1. **Tests are scattered and slow.** Dedup assertions live in four files:
   - `test/test_mechanism_enumeration.jl` (spec-level `_canonicalize!`, `_dedup_key`, `dedup!`)
   - `test/test_identify_rate_equation.jl` (canonical-hash basic + Pattern-A LDH)
   - `test/test_eq_hash_dedup.jl` (4-mech LDH cluster + section labels + fitted_params shape — depends on `dedup_investigation/cv_results.csv`)
   - `test/test_dedup_csv_replay.jl` (replays ~22k rows from `params_estimate_{5,6,7,8}.csv` — takes >30 min because each unique mechanism type triggers a fresh `@generated rate_equation_string` compile)

2. **Source code is scattered.** The rate-equation canonical hash
   (`_canonicalize_rate_eq_with_map`, `_canonical_rate_eq_hash`,
   `_canonical_rate_eq_hash_data`, `_factor_sort_key`, `_sort_run_factors`)
   lives in `src/identify_rate_equation.jl` even though it is pure dedup
   logic with no beam-search dependency.

## Goal

Consolidate. No behavior changes.

- Move the rate-equation canonicalizer into `src/mechanism_enumeration.jl`
  alongside the existing spec-level `_canonicalize!` / `_dedup_key` /
  `dedup!` block.
- Replace the slow CSV-replay test and the cv_results.csv–dependent cluster
  test with targeted, minimal hand-synthesized exemplars — one test per
  specific dedup behavior.
- Consolidate all dedup tests under one top-level testset in
  `test/test_mechanism_enumeration.jl`.
- Drop `dedup_investigation/` (CSVs + status docs + diagnostic script).

Test-suite runtime for the dedup block goes from multi-minute (CSV replay)
to <30 s (handful of minimal mechanism compilations).

## Non-goals

- No change to the canonicalizer algorithm, sort keys, regex set, or
  section-stripping logic. Pure file moves.
- No change to `_build_kinetic_rename_map` Pass 1 / Pass 2 / transitive
  closure. That code stays in `src/rate_eq_derivation.jl` — it is
  intertwined with rate-equation rendering (it folds equivalent params
  into v and emits the constraint sections at the same time).
- No change to `ANNOTATION_SUBSTITUTED`. Stays in `rate_eq_derivation.jl`
  as the display constant it is; cross-module reference from the
  canonicalizer is fine.
- No change to `_project_cached_params`. Stays in
  `src/identify_rate_equation.jl` — its only call site is beam-search
  stage 1, and pulling it into mechanism_enumeration would split it from
  its single caller for organizational tidiness only.
- No package version bump. No API change. Exports list unchanged.

## Code moves

### From `src/identify_rate_equation.jl` → `src/mechanism_enumeration.jl`

Move these definitions verbatim:

- `_canonicalize_rate_eq_with_map(m::AbstractEnzymeMechanism)` — rate-equation-string → (canonical text, name_map). Strips destructure header lines, section header lines, and single-symbol `(substituted into v)` equality lines. Renames parameters by first-appearance order. Sorts factors within each multiplicative run.
- `_canonical_rate_eq_hash_data(m::AbstractEnzymeMechanism)` — returns (UInt64 hash, 16-char hex, name_map).
- `_canonical_rate_eq_hash(m::AbstractEnzymeMechanism)` — UInt64 hash.
- `_sort_run_factors(run::AbstractString)` — private helper.
- `_factor_sort_key(f::AbstractString)` — private helper.

Insert location: at the bottom of `src/mechanism_enumeration.jl`, after
the existing `dedup!` definition. The file's "Dedup" comment header
becomes the natural anchor for both spec-level and rate-equation
dedup logic.

The canonicalizer references `ANNOTATION_SUBSTITUTED` (defined in
`rate_eq_derivation.jl`). Since both files are inside the
`EnzymeRates` module, the cross-file reference works without any
import change. Confirm by running the test suite after the move.

### Stays in `src/identify_rate_equation.jl`

- `_project_cached_params(...)` and all beam-search machinery.

### Stays in `src/rate_eq_derivation.jl`

- `_build_kinetic_rename_map` (Pass 1 user-defined kinetic-group merges +
  Pass 2 single-symbol Wegscheider absorption + transitive closure).
- User-defined / Wegscheider / Haldane section emission inside both
  `rate_equation_string` methods (EnzymeMechanism + AllostericEnzymeMechanism).
- `ANNOTATION_SUBSTITUTED` constant.

## Test consolidation

### New top-level testset in `test/test_mechanism_enumeration.jl`

Insert as a new top-level `@testset "Rate-equation canonical hash dedup"`
block. Names describe the equivalence being asserted; no
"Source A / B / C" jargon.

```
@testset "Rate-equation canonical hash dedup"
    @testset "_factor_sort_key sorts p_i atoms numerically"
    @testset "_sort_run_factors reorders mixed run by p_i"
    @testset "Hash is deterministic across repeated calls"
    @testset "Hash hex string is 16 lowercase hex chars"
    @testset "Distinct mechanisms produce distinct hashes"
    @testset "Steps sharing a kinetic_group → only representative K_i in hash"
    @testset "Wegscheider single-symbol tie K_i = K_j → K_i absent from hash"
    @testset "Chained kinetic-group + Wegscheider renames close transitively"
    @testset "Allosteric T-state K_i_T renamed away in canonical hash"
    @testset "rate_equation_string emits User-defined / Wegscheider / Haldane section labels"
    @testset "Hash-equivalent mechanisms share fitted_params shape"
end
```

### Per-test plan

**`_factor_sort_key`** — pure unit. Feed `"p_1"`, `"p_10"`, `"p_2 ^ 3"`,
`"E_total"`. Assert returned tuples sort in the expected order
(numerical-by-p_i for `p_*`, lex for everything else). No mechanism
compile needed.

**`_sort_run_factors`** — pure unit. `"p_3 * p_1 * p_2"` →
`"p_1 * p_2 * p_3"`. Single string transform.

**Hash is deterministic** — one minimal `@enzyme_mechanism` uni-uni,
call `_canonical_rate_eq_hash` twice, assert equal.

**Hash hex string** — same minimal uni-uni; call
`_canonical_rate_eq_hash_data`; assert `length(hex) == 16` and
`all(c -> c in "0123456789abcdef", hex)`.

**Distinct mechanisms produce distinct hashes** — two `@enzyme_mechanism`
variants with genuinely different v polynomials (e.g., different
catalytic topologies). Assert unequal hashes.

**Steps sharing a kinetic_group → only representative K_i in hash** —
uni-uni with a dead-end inhibitor mirror step. The mirror step shares
the catalytic binding step's `kinetic_group`, so K_mirror gets absorbed
into K_rep at Pass 1. Construct two specs: one where the partition is
expressed via mirror-step kinetic_group sharing, one where the same v
is produced via a different step ordering with the same partition.
Assert equal hashes AND assert the canonical text contains the rep
symbol once (not the absorbed one). This is also the regression
covering what was previously called "Pattern-A LDH duplicates".

**Wegscheider single-symbol tie K_i = K_j → K_i absent from hash** —
pick the smallest topology where the Haldane/Wegscheider closure
produces a single-symbol equality (likely 2-substrate ordered-binding
with iso). Two specs that produce the same v after the absorption (one
with the tie collapsed at spec time, one expressing it via the
Wegscheider rename map) → equal hashes.

**Chained kinetic-group + Wegscheider renames close transitively** —
one mechanism designed so Pass 1 creates `K_i → K_j` AND Pass 2
creates `K_j → K_k`. Assert the canonical text contains only the
chain-end symbol — no residual K_j. The failure mode is "K_j leaks
into v", so a single mechanism + canonicalize is sufficient (no pair
needed). Finding the smallest topology that exhibits the chain is
done during implementation; falling back to 3-substrate if uni/bi-uni
doesn't produce one.

**Allosteric T-state K_i_T renamed away** — move the canonicalizer-
invariant assertion from the existing `t_state_dead with :NonequalRT`
testset (currently around `test/test_mechanism_enumeration.jl:4470`)
into this block:

```julia
canon, _ = EnzymeRates._canonicalize_rate_eq_with_map(m)
@test !occursin(r"\bK\d+_T\b", canon)
@test !occursin(r"\bk\d+[fr]_T\b", canon)
```

Leave the `parameters(m, Full)` assertions (`:K1_T in params_full`,
etc.) in their original block — those are parameters-API regressions,
not dedup regressions.

**Section labels render correctly** — two `@enzyme_mechanism` cases:
one with a user-defined kinetic-group merge that emits
`# User defined constraints:`, one whose Wegscheider closure emits
`# Wegscheider constraints:`. Assert both header strings appear via
`occursin`. Haldane closure is exercised implicitly by the
Wegscheider mechanism (binding-K mechanisms always have a Haldane
constraint involving `:Keq`).

**Hash-equivalent mechanisms share fitted_params shape** — reuse one
of the equal-hash pairs from above. Assert
`length(fitted_params(m1)) == length(fitted_params(m2))` AND
`sort(_fp_kind.(fitted_params(m1))) == sort(_fp_kind.(fitted_params(m2)))`.
The `_fp_kind` helper currently lives in
`test/test_dedup_csv_replay.jl:15-31` — move into this test block
as a local helper.

### Deletions

- `test/test_eq_hash_dedup.jl` (whole file).
- `test/test_dedup_csv_replay.jl` (whole file).
- The two canonical-hash testsets in `test/test_identify_rate_equation.jl`
  at lines 492-559:
  - `@testset "canonical rate-equation hash: basic"` — covered by
    Hash deterministic / hex / Distinct tests.
  - `@testset "canonical hash collapses Pattern-A LDH duplicates"` —
    covered by the "Steps sharing a kinetic_group" minimal-exemplar test.
- `dedup_investigation/` directory and all its contents
  (LDH_data.csv, cv_results.csv, params_estimate_{5,6,7,8}.csv,
  refit_diagnostic.jl, status_2026-05-11.md,
  investigation_eq_hash_duplication.md) — ~22 MB freed.

### `test/runtests.jl` edits

Remove these two lines:
```
include("test_eq_hash_dedup.jl")
include("test_dedup_csv_replay.jl")
```

## Verification

Implementation is done when:

1. `julia --project -e 'using Pkg; Pkg.test()'` passes (including Aqua + JET).
2. The new `@testset "Rate-equation canonical hash dedup"` block runs
   in under ~30 s in isolation (vs the multi-minute CSV replay it
   replaces).
3. `git status` shows the four target files removed (`test_eq_hash_dedup.jl`,
   `test_dedup_csv_replay.jl`, the deleted canonical-hash testsets in
   `test_identify_rate_equation.jl`, and the entire `dedup_investigation/`
   directory).
4. `grep -r "_canonical_rate_eq_hash\|_canonicalize_rate_eq_with_map" src/`
   shows the canonicalizer definitions in `mechanism_enumeration.jl`
   and references (if any) in `identify_rate_equation.jl`.
5. No new exports. No API change. Public `parameters`, `fitted_params`,
   `rate_equation`, `rate_equation_string`, `EnzymeMechanism`,
   `AllostericEnzymeMechanism`, beam-search interface all behave
   identically to pre-cleanup.

## Risks

- **Cross-file reference for `ANNOTATION_SUBSTITUTED`:** the canonicalizer
  in its new home references a constant defined in another file. Both
  files share the `EnzymeRates` module, so the reference resolves at
  load time regardless of include order. If `mechanism_enumeration.jl`
  is included *before* `rate_eq_derivation.jl` in
  `src/EnzymeRates.jl`, the `const ANNOTATION_SUBSTITUTED = ...`
  binding won't exist yet when the canonicalizer source is parsed —
  but Julia binds constants at module-finalization time, not at
  parse time, so this is fine for function bodies. Confirmed during
  implementation by reading `src/EnzymeRates.jl`'s include order; if
  needed, swap include order (cheap and safe).

- **Finding the minimal chained-rename exemplar:** designing the
  smallest mechanism that produces both a Pass 1 kinetic-group share
  AND a Pass 2 Wegscheider single-symbol tie chaining through that
  share takes more thought than the per-source minimal cases. If
  uni/bi-uni topologies don't naturally produce the chain, fall back
  to a 3-substrate case — still minimal, still hand-written, still
  fast.

- **Test pruning losing real-world coverage:** the CSV replay was a
  noisy but real signal that canonicalization worked on 22k mechanisms
  pulled from actual fits. Replacing with targeted minimal exemplars
  trades broad coverage for focused coverage. Acceptable because each
  targeted test directly asserts the behavior the CSV replay was
  indirectly checking, and a regression that escapes the targeted
  tests would necessarily be a bug in a code path none of them
  exercise — which means the targeted set is incomplete and should
  be extended, not "we needed CSV replay".
