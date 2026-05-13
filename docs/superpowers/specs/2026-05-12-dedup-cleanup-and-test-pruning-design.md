# Dedup cleanup and test pruning — design

Date: 2026-05-12 (revised after REPL investigation + reviewer audit)
Branch: `improve-deduplication-of-ident-equations`
Status: design (awaiting plan)

## Problem

Deduplication of equivalent rate equations works correctly on this branch.
What remains is code and test cleanup:

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

## What the canonicalizer actually does (investigation finding)

REPL investigation across nine candidate topologies established what the
canonicalizer's job actually is and what's testable on minimal exemplars:

- **Pass 1 (user-defined kinetic-group merges)** — fires whenever a
  mechanism has multiple steps sharing a `kinetic_group` integer.
  Renders a `# User defined constraints:` section with
  `K_dup = K_rep (substituted into v)` lines. Verified on a 6-step
  random bi-uni with substrate-side mirror sharing.
- **Polynomial-equivalence detection** — the canonicalizer's main job:
  two mechanisms with different enzyme-form graphs may produce
  equivalent v polynomials after Pass-1 absorption. The canonicalizer
  detects equivalence via section-label stripping + first-appearance
  parameter rename + multiplicative-factor sorting. Verified only on
  the LDH Pattern-A 11-step literal pair (smaller hand-tuned topologies
  did not produce graph-distinct yet v-equivalent pairs in nine attempts).
- **Pass 2 (single-symbol Wegscheider absorption + transitive closure)**
  — present in code at `src/rate_eq_derivation.jl:152-167` as a defensive
  pipeline. **Does not fire on real or minimal mechanisms.** Even the
  LDH Pattern-A literal — the canonicalizer's original motivating case
  — emits no `# Wegscheider constraints:` section; all its single-symbol
  equalities come from Pass 1. Gaussian elimination of thermodynamic
  constraints on hand-synthesizable topologies always produces
  multi-symbol Wegscheider RHSes (or no Wegscheider closure at all when
  Pass-1 sharing has collapsed the cycle constraints).
- **Allosteric T-state K_i_T renaming** — separate codepath for allosteric
  mechanisms with `:NonequalRT` site tags. Canonicalizer must rename
  K_i_T tokens away alongside K_i tokens.

This finding overrides the original "Sources A / B / C" framing: there
is **one** practical dedup source (Pass 1) plus the canonicalizer's
polynomial-equivalence machinery. Tests should target what's testable.

## Goal

Consolidate. No behavior changes.

- Move the rate-equation canonicalizer into `src/mechanism_enumeration.jl`
  alongside the existing spec-level `_canonicalize!` / `_dedup_key` /
  `dedup!` block.
- Replace the slow CSV-replay test and the cv_results.csv–dependent
  cluster test with targeted exemplars — one minimal Pass-1 unit test
  + the LDH Pattern-A literal pair for polynomial-equivalence regression.
- Consolidate all dedup tests under one new top-level `@testset`
  (peer to the file's other top-level testsets, NOT nested inside the
  existing `@testset "Mechanism Enumeration"` block at line 487).
- Drop `dedup_investigation/` (CSVs + status docs + diagnostic script).

Test-suite runtime for the dedup block goes from multi-minute (CSV replay)
to <30 s (handful of mechanism compilations).

## Non-goals

- No change to the canonicalizer algorithm, sort keys, regex set, or
  section-stripping logic. Pure file moves.
- No change to `_build_kinetic_rename_map`. The Pass-2 transitive-closure
  code stays as defensive infrastructure even though it doesn't fire on
  realistic mechanisms.
- No change to `ANNOTATION_SUBSTITUTED`. Stays in `rate_eq_derivation.jl`
  as the display constant it is; cross-module reference from the
  canonicalizer is fine — `mechanism_enumeration.jl` is included after
  `rate_eq_derivation.jl` in `src/EnzymeRates.jl`, and the function-body
  reference resolves at call time (not parse time).
- No change to `_project_cached_params`. Stays in
  `src/identify_rate_equation.jl` — its only call site is beam-search
  stage 1.
- No package version bump. No API change. Exports list unchanged.
- This is a **pure refactor**, not new behavior. The "Write failing
  test first" TDD ordering does not apply; tests are added against
  existing-and-correct code and remain passing throughout.

## Code moves

### From `src/identify_rate_equation.jl` → `src/mechanism_enumeration.jl`

Move these definitions verbatim:

- `_canonicalize_rate_eq_with_map(m::AbstractEnzymeMechanism)` — rate-equation-string → (canonical text, name_map). Strips destructure header lines, section header lines, and single-symbol `(substituted into v)` equality lines. Renames parameters by first-appearance order. Sorts factors within each multiplicative run.
- `_canonical_rate_eq_hash_data(m::AbstractEnzymeMechanism)` — returns (UInt64 hash, 16-char hex, name_map).
- `_canonical_rate_eq_hash(m::AbstractEnzymeMechanism)` — UInt64 hash.
- `_sort_run_factors(run::AbstractString)` — private helper.
- `_factor_sort_key(f::AbstractString)` — private helper.

Insert location: at the bottom of `src/mechanism_enumeration.jl`, after
the existing `dedup!` definition. Move uses function-boundary text
matching (start at the `_canonicalize_rate_eq_with_map` docstring, end at
the closing `end` of `_canonical_rate_eq_hash`), not absolute line
numbers.

### Stays in `src/identify_rate_equation.jl`

- `_project_cached_params(...)` and all beam-search machinery.

### Stays in `src/rate_eq_derivation.jl`

- `_build_kinetic_rename_map` (Pass 1 user-defined kinetic-group merges +
  Pass 2 defensive transitive closure).
- User-defined / Wegscheider / Haldane section emission inside both
  `rate_equation_string` methods (EnzymeMechanism + AllostericEnzymeMechanism).
- `ANNOTATION_SUBSTITUTED` constant.

## Test consolidation

### New top-level testset (peer, not nested)

Insert as a NEW top-level `@testset` block AFTER the closing
`end # top-level testset` of the file's existing
`@testset "Mechanism Enumeration"` (currently line 4537). The new
block is a peer to the existing top-level testsets at lines 139, 182,
250, 487.

```julia
@testset "Rate-equation canonical hash dedup"
    @testset "_factor_sort_key sort order"
    @testset "_sort_run_factors sort order"
    @testset "Hash is deterministic across repeated calls"
    @testset "Hash hex string is 16 lowercase hex chars"
    @testset "Distinct mechanisms produce distinct hashes"
    @testset "Pass-1 kinetic-group merge: User-defined section + canonical text invariant"
    @testset "LDH Pattern-A: graph-distinct mechanisms with equivalent v hash equally"
    @testset "Allosteric T-state K_i_T renamed away in canonical hash"
    @testset "rate_equation_string emits section labels"
    @testset "Hash-equivalent mechanisms share fitted_params shape"
end
```

10 testsets. Pure-unit testsets run in microseconds. The mechanism-driven
testsets compile a handful of distinct `EnzymeMechanism` types — total
dedup-block runtime target: **<30 s**.

### Per-test plan

**`_factor_sort_key sort order`** — pure unit. Assert observable ordering
only (not tuple-shape internals):
```julia
@test _factor_sort_key("p_1") < _factor_sort_key("p_2")
@test _factor_sort_key("p_2") < _factor_sort_key("p_10")
@test _factor_sort_key("p_99") < _factor_sort_key("E_total")
```
No assertions on individual tuple slot values — those lock to a specific
encoding and would break on an internal refactor without a real
regression.

**`_sort_run_factors sort order`** — pure unit. Three string-transform
assertions (numerical p_i ordering, exponent preservation, non-p atoms
sort to end).

**Hash is deterministic** — one minimal `@enzyme_mechanism` uni-uni,
call `_canonical_rate_eq_hash` twice, assert equal.

**Hash hex string** — same uni-uni; call `_canonical_rate_eq_hash_data`;
assert length=16 + all lowercase hex chars.

**Distinct mechanisms** — uni-uni vs bi-bi `@enzyme_mechanism`, hashes
differ.

**Pass-1 kinetic-group merge** — verified random bi-uni exemplar:

```julia
m = EnzymeMechanism(
    ((:A, :B), (:P,), ()),
    (((:E, :A), (:E_A,), true, 1),
     ((:E_B, :A), (:E_A_B,), true, 1),
     ((:E, :B), (:E_B,), true, 2),
     ((:E_A, :B), (:E_A_B,), true, 2),
     ((:E_A_B,), (:E_P,), false, 3),
     ((:E, :P), (:E_P,), true, 4)))
```

Assert (a) `rate_equation_string(m)` contains `# User defined constraints:`
and `(substituted into v)`, (b) `_canonicalize_rate_eq_with_map(m)`
canonical text contains no raw `K\d+\b` or `k\d+[fr]\b` tokens.
**This pair was verified to produce the expected rendering in REPL.**

**LDH Pattern-A polynomial-equivalence** — copy the `m_a` / `m_b` pair
verbatim from `test/test_identify_rate_equation.jl:528-559`. Assert
`_canonical_rate_eq_hash(m_a) == _canonical_rate_eq_hash(m_b)`. **Both
hashes confirmed equal in REPL (`4711f0996b051276`).** This is the
real-world regression the canonicalizer was built for; smaller
hand-tuned topologies tried in investigation did not exhibit the
graph-distinct-but-v-equivalent property. Acknowledged in the test
comment that 11 steps is the minimal known case.

**Allosteric T-state K_i_T renamed away** — move the canonicalizer-
invariant assertion from `test/test_mechanism_enumeration.jl:4496-4502`
into this block:
```julia
canon, _ = EnzymeRates._canonicalize_rate_eq_with_map(m)
@test !occursin(r"\bK\d+_T\b", canon)
@test !occursin(r"\bk\d+[fr]_T\b", canon)
```
Leave the `parameters(m, Full)` assertions in their original allosteric
block (parameters-API regression, not dedup regression).

**`rate_equation_string` emits section labels** — two cases:
- The Pass-1 random bi-uni exemplar emits `# User defined constraints:`.
- The same exemplar (or any RE binding uni-uni) emits `# Haldane constraints:`.

`# Wegscheider constraints:` is **NOT** asserted here — investigation
confirmed that no minimal hand-synthesizable mechanism produces a
`# Wegscheider constraints:` section with single-symbol absorptions.
The LDH Pattern-A test is the indirect regression for the
section-stripping pipeline as a whole; if the rendering side ever
emits malformed section headers, the LDH hash equality will break.

**Hash-equivalent mechanisms share fitted_params shape** — use the LDH
Pattern-A `m_a`/`m_b` pair (the only verified hash-equivalent pair).
Assert `length(fp1) == length(fp2)` and `sort(_fp_kind.(fp1)) ==
sort(_fp_kind.(fp2))`. The `_fp_kind` helper lifts from
`test/test_dedup_csv_replay.jl:15-31`. Wrap inside a `let` block to
keep it scope-local to the testset (Julia's `function f(...) end`
inside a `@testset begin ... end` block leaks the binding to the
file's module scope).

### Deletions

- `test/test_eq_hash_dedup.jl` (whole file).
- `test/test_dedup_csv_replay.jl` (whole file).
- The two canonical-hash testsets in `test/test_identify_rate_equation.jl`
  at lines 492-559:
  - `@testset "canonical rate-equation hash: basic"` — covered by
    Hash deterministic / hex / Distinct tests.
  - `@testset "canonical hash collapses Pattern-A LDH duplicates"` —
    the m_a/m_b literal moves into the new dedup block under the
    "LDH Pattern-A" test.
- `dedup_investigation/` directory and all its contents (including
  the currently-untracked `refit_results.csv`). Confirmed OK to delete
  with user; ~22 MB freed.

### `test/runtests.jl` edits

Remove the two `include` lines for the deleted test files.

## Verification

Implementation is done when:

1. `julia --project -e 'using Pkg; Pkg.test()'` passes (including Aqua + JET).
2. The new `@testset "Rate-equation canonical hash dedup"` block runs
   in under ~30 s in isolation.
3. `git status` shows the four target files removed.
4. `grep -n "^function _canonical\|^function _factor_sort_key\|^function _sort_run_factors\|^function _canonicalize_rate_eq" src/` shows the definitions in `mechanism_enumeration.jl` (and not in `identify_rate_equation.jl`).
5. No new exports. No API change.

## Risks

- **Cross-file reference for `ANNOTATION_SUBSTITUTED`:** the canonicalizer
  in its new home references a constant defined in another file. Both
  files share the `EnzymeRates` module; function-body references resolve
  at call time, and `rate_eq_derivation.jl` is included before
  `mechanism_enumeration.jl` in `src/EnzymeRates.jl`, so even
  parse-time resolution would work. No action needed.

- **LDH literal as "minimal" exemplar:** 11 steps is not minimal in
  absolute terms, but it is the smallest verified Pattern-A
  graph-distinct-but-v-equivalent case after nine investigation attempts
  on smaller topologies. The test comment documents this empirical
  finding explicitly.

- **Test pruning losing CSV-replay coverage:** the CSV replay was a
  noisy but real signal that canonicalization worked on 22k mechanisms
  pulled from actual fits. The targeted tests directly assert the
  pipeline (Pass 1, polynomial-equivalence, allosteric T-state,
  section labels, fitted_params shape). The class of regression the
  CSV replay would catch but targeted tests miss is "different
  mechanism with same loss yields different hash" — i.e., a polynomial-
  equivalence detection gap on some untested topology. Accepted: a
  regression there would also necessarily mean an untested code path,
  which means the targeted test set should be extended, not "we needed
  CSV replay."
