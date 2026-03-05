# Implementation Plan: Enumeration Pipeline Restructure

## Overview

Restructure `enumerate_mechanisms` to:
1. Remove activator expansion (stage 2) — replaced by MWC conformational switching
2. Change RE/SS default from all-SS to all-RE-except-first-catalytic-SS
3. Cap King-Altman rate matrix dimension at G ≤ 7 (random-order bi-bi ceiling)
4. Add regulator partitioning (dead-end vs allosteric)
5. Produce both `EnzymeMechanism` and `OligomericEnzymeMechanism` candidates
6. Add `catalytic_n` parameter, remove `max_forms`

## Key Design Decisions

- **G = number of RE groups** (rate matrix dimension in `_raw_symbolic_rate_polys`). Cap at 7.
- **Baseline**: all steps RE except first isomerization in canonical edge order → SS. This guarantees at least 1 SS step. Stage 4 adds more SS steps.
- **G ≥ 2 required**: mechanisms where all forms merge into 1 RE group are invalid (SS self-loops don't contribute to spanning trees).
- **Activator stage removed**: essential/non-essential activators are now handled by `OligomericEnzymeMechanism` with `NConf=2` (R/T conformational switching).
- **Regulator partitioning**: each regulator independently assigned to dead-end (in CatalyticMech) or allosteric (RegSite on OligomericEnzymeMechanism). Enumerate all 2^n_reg partitions.
- **OligomericEnzymeMechanism**: `NConf=2` always (NConf=1 is redundant with plain EnzymeMechanism). `CatalyticN` user-provided. Allosteric regulator multiplicity ∈ 1:CatalyticN. Produced even for no-regulator reactions (substrates/products as allosteric ligands).
- **catalytic_n=1** still valid for OligomericEnzymeMechanism — produces essential/non-essential activator behavior via R/T states.

## New Pipeline

```
Stage 1: Catalytic topologies (unchanged algorithm)
Stage 2: Dead-end inhibitors for a SUBSET of regulators (parameterized by which regs are dead-end)
Stage 3: RE/SS expansion with G≤7 cap + equivalence constraints
         Baseline: all RE except first canonical isomerization → SS
         Expand: try all additional RE→SS conversions where 2 ≤ G ≤ 7
Stage 4: Regulator partition + Oligomeric expansion
         For each 2^n_reg partition of regulators into {dead-end, allosteric}:
           - Run stages 1-3 with dead-end regulators only
           - Produce EnzymeMechanism candidates
           - Produce OligomericEnzymeMechanism candidates (NConf=2, CatalyticN=user)
             with allosteric regulators as RegSites (multiplicity ∈ 1:CatalyticN)
```

## Agent Breakdown (5 sequential agents)

---

### Agent 1: Remove Activator Stage

**Goal**: Remove `WithActivator` stage and activator expansion code. Simplify pipeline to Catalytic → DeadEnd → FullEnumeration. All existing tests must still pass after adjusting expected counts.

**Files modified**:
- `src/mechanism_enumeration.jl`
- `test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl`
- `test/test_mechanism_enum_of_enz_reaction.jl`

#### RED: Update tests first

1. **Update `EnumerationTestSpec`**: Remove `expected_n_cat_with_act` field.

2. **Update `test_mechanism_enum_of_enz_reaction.jl`**:
   - Remove the `WithActivator()` stage test block (lines 13-18).
   - Remove `with_act` from the dead-end independent verification (replace with `catalytic` in `_compute_expected_dead_end_count` call).
   - Remove the `with_act` subset check in "Stage subset" testset.
   - Keep catalytic, dead-end, and full enumeration tests.

3. **Recompute `expected_n_cat_act_de`** (now just `expected_n_cat_de`):
   Without activator expansion, dead-end is computed directly from catalytic topologies:
   - Formula per catalytic topology: `(2^n_reg)^n_topo_forms`
   - Uni-Uni: `(2^0)^3 = 1` (no regulators) → 1
   - Uni-Uni 1 Reg: `(2^1)^3 = 8` (1 topology, 3 forms, 1 reg) → 8
   - Uni-Uni 2 Reg: `(2^2)^3 = 64` (1 topology, 3 forms, 2 regs) → 64
   - Uni-Bi 1 Reg: sum over 3 catalytic topologies:
     - 2 sequential (4 forms each): `2 × (2^1)^4 = 32`
     - 1 random (5 forms): `(2^1)^5 = 32`
     - Total: 64
   - Bi-Bi 1 Reg: sum over 9 catalytic topologies (need to compute form counts per topology)
   - Bi-Bi PP: 10 (no regulators, same as before)
   - Bi-Bi Budget: 4 (no regulators, same as before)

4. **Recompute `expected_n_total`**: The RE/SS counts change because the dead-end base specs change (no activator-expanded topologies). The formula `_compute_expected_n_total` remains valid but is applied to fewer base specs. New totals must be computed.

5. **Update `_compute_expected_dead_end_count`**: Input is now `catalytic_specs` instead of `activator_specs` — same formula applies, just different input.

**IMPORTANT**: Run tests → they should FAIL (RED) because the implementation hasn't changed yet.

#### GREEN: Implement changes

1. **Remove from `src/mechanism_enumeration.jl`**:
   - Delete `WithActivator` struct (line 78).
   - Delete `_build_activator_options` function (lines 377-408).
   - Delete `_expand_activators` function (lines 419-447).

2. **Update `enumerate_mechanisms`**:
   - Remove `with_activators` step (lines 590-591).
   - Pass `catalytic` directly to `_expand_inhibitors` (instead of `with_activators`).
   - Remove `stage isa WithActivator && return with_activators` check.

3. Run tests → should PASS (GREEN).

#### Verification
- All 8 enumeration test specs pass with new expected counts.
- `test_enzyme_derivation.jl` still passes (mechanism definitions unchanged).
- `test_aqua_jet.jl` passes (no stale exports).

---

### Agent 2: New RE/SS Expansion with G ≤ 7 Cap

**Goal**: Change stage 4 from "enumerate all 2^n-1 RE/SS masks" to "baseline = all RE except first canonical isomerization SS, expand with G ≤ 7 cap." Only valid masks with 2 ≤ G ≤ 7 are kept.

**Files modified**:
- `src/mechanism_enumeration.jl`
- `test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl`
- `test/test_mechanism_enum_of_enz_reaction.jl`

#### RED: Update tests

1. **Add `_compute_re_group_count` test helper** in the test reaction definitions file:
   ```julia
   """Compute G (number of RE groups) for a given RE/SS mask on edges."""
   function _compute_re_group_count(edges, eq_steps, forms)
       # Union-find: merge form indices connected by RE edges
       ...
       return number_of_groups
   end
   ```

2. **Add targeted unit tests** for the G computation:
   - Uni-Uni 3 edges, all SS → G=3
   - Uni-Uni 3 edges, edge 2 SS only → G=1 (invalid)
   - Uni-Uni 3 edges, edges 1,2 SS → G=2 (valid)

3. **Update `_compute_expected_n_total`**: Replace with a new version that:
   - Identifies the first isomerization edge (baseline SS)
   - Iterates over subsets of remaining edges to make SS
   - For each subset, computes G via union-find
   - Keeps only masks with 2 ≤ G ≤ 7
   - For each valid mask, counts equivalence constraint combos
   - This is now a brute-force counting function (slow but correct for tests)

4. **Recompute all `expected_n_total` values** using the new counting function.
   Note: For reactions with no isomerization edge (shouldn't happen for valid catalytic cycles — every cycle has at least one isomerization), fall back to first edge as SS.

5. **Tests should FAIL** because old `_ress_variants` produces different counts.

#### GREEN: Implement

1. **Add `_compute_re_group_count` to `src/mechanism_enumeration.jl`**:
   ```julia
   """Compute G (number of RE groups) for given edges and eq_steps."""
   function _compute_re_group_count(edges, eq_steps, forms)
       form_indices = Set(Iterators.flatten(edges))
       parent = Dict(i => i for i in form_indices)
       function find(x)
           while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end
           x
       end
       for (idx, (a, b)) in enumerate(edges)
           eq_steps[idx] || continue  # only RE steps merge
           ra, rb = find(a), find(b)
           ra != rb && (parent[ra] = rb)
       end
       length(Set(find(i) for i in form_indices))
   end
   ```

2. **Add `_find_first_isomerization` helper**:
   ```julia
   """Find index of first isomerization edge in canonical order."""
   function _find_first_isomerization(edges, adj)
       for (i, (a, b)) in enumerate(edges)
           info = adj[minmax(a, b)]
           info.type == :isomerization && return i
       end
       return 1  # fallback: first edge
   end
   ```

3. **Rewrite `_ress_variants`**:
   ```julia
   function _ress_variants(spec, adj, forms; max_re_groups::Int=7)
       edges = spec.edges
       n = length(edges)
       iso_idx = _find_first_isomerization(edges, adj)
       equiv_groups = _find_equivalent_groups(edges, adj, forms)

       # Iterate over subsets of OTHER edges to make SS (iso_idx always SS)
       other_indices = [i for i in 1:n if i != iso_idx]
       Iterators.flatmap(0:(1 << length(other_indices)) - 1) do ss_mask
           eq_steps = fill(true, n)  # start all RE
           eq_steps[iso_idx] = false  # first isomerization always SS
           for (bit, idx) in enumerate(other_indices)
               if (ss_mask >> (bit - 1)) & 1 == 1
                   eq_steps[idx] = false  # make SS
               end
           end

           # Check G constraint
           G = _compute_re_group_count(edges, eq_steps, forms)
           (G < 2 || G > max_re_groups) && return Iterators.map(identity, ())

           # Equivalence group constraints (same logic as before)
           valid_groups = [g for g in equiv_groups
               if all(eq_steps[s] == eq_steps[g[1]] for s in g)]
           Iterators.map(0:(1 << length(valid_groups)) - 1) do constraint_mask
               constraints = ParamConstraint[...]  # same as current
               MechanismSpec(spec.reaction, edges, eq_steps, constraints)
           end
       end
   end
   ```

4. **Update `enumerate_mechanisms`**: pass `max_re_groups=7` to `_ress_variants`. Update the O(1) count formula (may need brute-force counting or a new closed-form formula for the G-capped case).

5. **Update stages 1-3 MechanismSpec construction**: Change from `MechanismSpec(reaction, edges)` (all SS default) to setting the baseline RE/SS mask with first isomerization as SS.

6. **O(1) count for MechanismIterator**: The old closed-form formula `2^(n-Σgᵢ) × ∏(2^gᵢ+2) - 2^k` no longer applies because of the G≤7 filter. Options:
   - Compute exact count by iterating over masks (may be slow for large n)
   - Use lazy iterator without precomputed count (change `MechanismIterator` to `SizeUnknown`)
   - Precompute count via brute force during `enumerate_mechanisms`

   **Decision**: For now, precompute count via brute force during `enumerate_mechanisms`. The mask iteration is 2^(n-1) which for n≤~20 edges is feasible. For very large mechanisms, the G≤7 cap means n is bounded anyway. If the number of edges is large, we can add a shortcut: if all-SS gives G ≤ 7, use the old formula.

#### Verification
- All enumeration test specs pass with new expected counts.
- New unit tests for `_compute_re_group_count` pass.
- Rate equation tests still pass (mechanism definitions unchanged).

---

### Agent 3: Regulator Partitioning + Dead-End Subset

**Goal**: Modify dead-end expansion to work with regulator subsets. Enumerate all 2^n_reg partitions of regulators into {dead-end, allosteric}. Track partition info in `MechanismSpec`.

**Files modified**:
- `src/mechanism_enumeration.jl`
- `test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl`
- `test/test_mechanism_enum_of_enz_reaction.jl`

#### RED: Update tests

1. **Extend `MechanismSpec`** (or add a wrapper): Add field to track which regulators are dead-end vs allosteric:
   ```julia
   struct MechanismSpec
       reaction::Any
       edges::Vector{Tuple{Int,Int}}
       equilibrium_steps::Vector{Bool}
       param_constraints::Vector{ParamConstraint}
       dead_end_regulators::Vector{Symbol}      # NEW
       allosteric_regulators::Vector{Symbol}     # NEW
   end
   ```

2. **Add new `EnumerationTestSpec` fields**:
   ```julia
   expected_n_partitions::Int  # 2^n_reg
   expected_n_total_all_partitions::Int  # sum over all partitions
   ```

3. **Add tests**: For Uni-Uni 1 Reg:
   - 2 partitions: {R dead-end} and {R allosteric}
   - {R dead-end}: dead-end expansion as Agent 1 computed (8 topologies), then RE/SS
   - {R allosteric}: no dead-end forms (just catalytic topologies), then RE/SS
   - Total = sum of both partitions

4. **Tests should FAIL** because implementation doesn't partition yet.

#### GREEN: Implement

1. **Update `MechanismSpec`** with new fields and backward-compatible constructor.

2. **Modify `_expand_inhibitors`**: Accept a set of `dead_end_regulators` (subset of all regulators). Only create dead-end forms for regulators in this subset.

3. **Modify `enumerate_mechanisms`**: Loop over all 2^n_reg partitions:
   ```julia
   for reg_mask in 0:(1 << n_reg) - 1
       dead_end_regs = [regs[i] for i in 1:n_reg if (reg_mask >> (i-1)) & 1 == 0]
       allosteric_regs = [regs[i] for i in 1:n_reg if (reg_mask >> (i-1)) & 1 == 1]
       # Run stages 1-3 with dead_end_regs only
       # Stage 4: RE/SS expansion
       # Tag MechanismSpecs with partition info
   end
   ```

4. **Optimization**: Catalytic topologies (stage 1) are the same for all partitions — compute once, reuse.

5. **`enumerate_enzyme_forms` change**: Currently generates forms for ALL regulators. For dead-end subset enumeration, we need forms with only the dead-end regulators. Options:
   - Filter existing forms to exclude allosteric regulator positions
   - Pass regulator subset to `enumerate_enzyme_forms`

   **Decision**: Pass dead-end regulator subset to `_expand_inhibitors`. The form enumeration stays the same (all forms exist), but `_expand_inhibitors` only uses dead-end regulator positions.

#### Verification
- Partition counts correct for all test specs.
- MechanismSpec carries correct partition info.
- Total mechanism count = sum over all partitions.

---

### Agent 4: OligomericEnzymeMechanism Enumeration

**Goal**: Produce `OligomericEnzymeMechanism` candidates from enumerated `MechanismSpec`s. Add `catalytic_n` parameter. Enumerate allosteric regulator multiplicities.

**Files modified**:
- `src/mechanism_enumeration.jl`
- `src/types.jl` (if needed for OligomericEnzymeMechanism constructor)
- `test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl`
- `test/test_mechanism_enum_of_enz_reaction.jl`

#### RED: Update tests

1. **Add `OligomericEnzymeMechanism` enumeration tests**:
   - New test specs or extension of existing ones with `catalytic_n=2`.
   - Uni-Uni no reg, catalytic_n=2: OligomericEnzymeMechanism candidates with NConf=2, CatalyticN=2, no RegSites. Count = same as EnzymeMechanism count (each EnzymeMechanism → one OligomericEnzymeMechanism).
   - Uni-Uni 1 Reg, catalytic_n=2: For allosteric partition, each catalytic MechanismSpec → multiplicity ∈ {1, 2} → 2 OligomericEnzymeMechanism per MechanismSpec. For dead-end partition, each → 1 OligomericEnzymeMechanism (no RegSites).
   - Verify constructed `OligomericEnzymeMechanism` has correct NConf, CatalyticN, RegSites.

2. **Add `OligomericEnzymeMechanism(spec::MechanismSpec; catalytic_n, n_conf)` constructor test**:
   - Construct from a known MechanismSpec, verify type parameters match.

3. **Tests should FAIL** because constructor doesn't exist yet.

#### GREEN: Implement

1. **Add `OligomericEnzymeMechanism(spec::MechanismSpec; ...)` constructor** in `src/mechanism_enumeration.jl`:
   ```julia
   function OligomericEnzymeMechanism(
       spec::MechanismSpec;
       catalytic_n::Int,
       allosteric_regulators::Vector{Tuple{Symbol, Int}}=Tuple{Symbol,Int}[],
   )
       cm = EnzymeMechanism(spec)
       mets = metabolites(cm)
       reg_sites = Tuple(
           (Tuple(reg), mult)
           for (reg, mult) in allosteric_regulators
       )
       OligomericEnzymeMechanism{
           mets, typeof(cm), catalytic_n, reg_sites, 2
       }()
   end
   ```

2. **Extend `enumerate_mechanisms`**:
   - Add `catalytic_n::Int=1` parameter.
   - In the regulator partition loop, for each partition:
     - Produce `EnzymeMechanism` specs (Path A, as before).
     - Produce `OligomericEnzymeMechanism` specs (Path B):
       - For each allosteric regulator, enumerate multiplicity 1:catalytic_n.
       - Cartesian product of multiplicities across allosteric regulators.
       - Each combo → one OligomericEnzymeMechanism candidate per base MechanismSpec.

3. **Iterator output type**: The iterator now yields a union or a wrapper that indicates whether to construct EnzymeMechanism or OligomericEnzymeMechanism. Options:
   - **Option A**: Separate iterators (one for each type).
   - **Option B**: Yield `MechanismSpec` with oligomeric fields — caller decides constructor.
   - **Option C**: New wrapper `MechanismCandidate` that holds either an EnzymeMechanism or OligomericEnzymeMechanism.

   **Decision**: Option B — `MechanismSpec` already has `allosteric_regulators` from Agent 3. Add `catalytic_n` and `n_conf` fields. The `identify_rate_equation` (future) checks these fields to call the right constructor. Add a convenience `compile_mechanism(spec)` function that dispatches.

4. **Count update**: Total count now includes multiplicity expansion. For each base spec with k allosteric regulators and catalytic_n=N:
   - EnzymeMechanism: 1
   - OligomericEnzymeMechanism: N^k (Cartesian product of multiplicities 1:N)
   - But if k=0 and catalytic_n > 0: still 1 OligomericEnzymeMechanism (no RegSites, just cooperative)
   - If catalytic_n == 0 or not provided: no OligomericEnzymeMechanism

   Wait, we said catalytic_n=1 is valid for OligomericEnzymeMechanism. So:
   - catalytic_n=1: OligomericEnzymeMechanism with CatalyticN=1, NConf=2
   - allosteric reg multiplicity: only 1 (since CatalyticN=1)
   - So N^k = 1^k = 1 per partition with allosteric regs

#### Verification
- OligomericEnzymeMechanism construction from MechanismSpec works.
- Correct NConf, CatalyticN, RegSites in constructed mechanisms.
- `rate_equation` works on constructed OligomericEnzymeMechanism.
- Total counts correct.

---

### Agent 5: API Cleanup + Integration

**Goal**: Clean up API, remove `max_forms`, add `catalytic_n` to public interface, update SPEC.md, ensure full test suite passes.

**Files modified**:
- `src/mechanism_enumeration.jl`
- `src/EnzymeRates.jl`
- `SPEC.md`
- `.claude/CLAUDE.md`
- All test files as needed

#### Tasks

1. **Remove `max_forms` parameter** from `enumerate_mechanisms`. The G≤7 cap is the sole complexity control. Add `max_re_groups::Int=7` as the replacement parameter.

2. **Export new symbols** if needed: `compile_mechanism` or similar.

3. **Update `EnumerationStage` types**:
   - `Catalytic()` — unchanged
   - `WithDeadEnd()` — works with regulator subset
   - `FullEnumeration()` — includes RE/SS + oligomeric expansion
   - Remove `WithActivator` if not done in Agent 1.

4. **Update SPEC.md**: Reflect new pipeline, new parameters, OligomericEnzymeMechanism enumeration.

5. **Update `.claude/CLAUDE.md`**: Update architecture decisions, source layout, enumeration description.

6. **Integration tests**:
   - Round-trip: enumerate → compile → rate_equation for both mechanism types.
   - Verify no mechanism has G > 7.
   - Verify every mechanism has at least 1 SS step.
   - Verify OligomericEnzymeMechanism candidates have NConf=2.

7. **Run full test suite**: `julia --project -e 'using Pkg; Pkg.test()'`

8. **Performance check**: Enumeration time for Bi-Bi + 1 Reg should be reasonable (< 30s).

#### Verification
- Full test suite passes.
- Aqua + JET pass (no stale exports, no type instabilities).
- No regressions in existing mechanism derivation tests.

---

### Agent 6: Simplify, Commit, PR, Review

**Goal**: Run code simplification, commit all changes, push, create PR, and run code review.

#### Tasks

1. **Run `/simplify`**: Review all changed code for reuse, quality, efficiency. Fix any issues found. Ensure 92-char line length, 4-space indentation, no dead code.

2. **Commit**: Stage all changes and commit with a descriptive message summarizing the pipeline restructure:
   ```
   Restructure mechanism enumeration: remove activators, add G≤7 cap, oligomeric expansion

   - Remove activator stage (essential/non-essential) — replaced by MWC NConf=2
   - Change RE/SS baseline to all-RE-except-first-catalytic-SS
   - Cap King-Altman rate matrix dimension at G≤7 (random-order bi-bi ceiling)
   - Add regulator partitioning (dead-end vs allosteric)
   - Enumerate OligomericEnzymeMechanism candidates (NConf=2, variable CatalyticN)
   - Replace max_forms with max_re_groups parameter
   - Add catalytic_n parameter to enumerate_mechanisms
   ```

3. **Push** to the current branch (`add-identify_rate_equation-functionality`).

4. **Create PR** against `main` with:
   - Title matching the commit summary
   - Body describing the motivation (too many complex mechanisms for finite data), the design decisions, and the new pipeline structure

5. **Run `/code-review`**: Review the PR for correctness, style, and adherence to project conventions.

---

## Expected Count Changes Summary

The exact expected counts for each test spec will need to be computed during implementation. Key factors:

| Spec | Old total | Direction | Reason |
|------|-----------|-----------|--------|
| Uni-Uni | 7 | ↓ | G≤7 filters some masks |
| Uni-Uni 1 Reg | 1,779 | ↓↓ | No activator expansion, G cap |
| Uni-Uni 2 Reg | 24,646,535 | ↓↓↓ | No activator expansion, G cap |
| Uni-Bi 1 Reg | 435,521 | ↓↓ | No activator expansion, G cap |
| Bi-Bi 1 Reg | 114,684,452 | ↓↓↓ | No activator expansion, G cap |
| Bi-Bi PP | 2,157 | ↓ | G cap |
| Bi-Bi Budget | 124 | removed | max_forms removed |

With OligomericEnzymeMechanism expansion (catalytic_n > 1), total counts increase again but from a much smaller base.

## Test File Changes Summary

| File | Changes |
|------|---------|
| `reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl` | Rewrite `EnumerationTestSpec`, recompute all counts, add oligomeric specs |
| `test_mechanism_enum_of_enz_reaction.jl` | Remove activator stage test, add G-cap tests, add oligomeric tests |
| `mechanism_definitions_for_test_enzyme_derivation.jl` | No changes (mechanism definitions are manual, not from enumeration) |
| `test_enzyme_derivation.jl` | No changes expected |
| Other test files | No changes expected |

## Risk Assessment

- **Highest risk**: Computing correct expected counts for the new pipeline. Mitigated by brute-force counting in test helpers.
- **Medium risk**: `OligomericEnzymeMechanism` constructor from `MechanismSpec` — need to ensure type parameters match DSL-constructed ones. Mitigated by comparing against existing test mechanisms.
- **Low risk**: Removing activator expansion — straightforward deletion.

---

## Progress Log

### Agent 1: Remove Activator Stage — COMPLETED

**Changes made:**
- `src/mechanism_enumeration.jl`:
  - Removed `WithActivator` struct from `EnumerationStage` hierarchy
  - Removed `_build_activator_options` function (was lines 377-408)
  - Removed `_expand_activators` function (was lines 419-447)
  - Updated `enumerate_mechanisms` to pass `catalytic` directly to `_expand_inhibitors` (skipping activator stage)
  - Updated pipeline docstring to reflect 3-stage pipeline (Catalytic → WithDeadEnd → FullEnumeration)
  - Simplified `_expand_inhibitors`: removed `activator_positions` logic (now dead since catalytic topologies never have regulators bound), renamed `inhibitor_positions` to `reg_positions`
  - Updated `_expand_inhibitors` docstring

- `test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl`:
  - Removed `expected_n_cat_with_act` field from `EnumerationTestSpec`
  - Renamed `expected_n_cat_act_de` to `expected_n_cat_de`
  - Recomputed all expected counts (dead-end and total RE/SS)

- `test/test_mechanism_enum_of_enz_reaction.jl`:
  - Removed `WithActivator()` stage test block
  - Updated dead-end verification to use `catalytic` instead of `with_act`
  - Removed `with_act` subset check in "Stage subset" testset

**New expected counts (activator stage removed):**

| Spec | cat | de (was cat_act_de) | total |
|------|-----|---------------------|-------|
| Uni-Uni | 1 | 1 (was 1) | 7 (unchanged) |
| Uni-Uni 1 Reg | 1 | 8 (was 10) | 808 (was 1779) |
| Uni-Uni 2 Regs | 1 | 64 (was 228) | 3089511 (was 24646535) |
| Uni-Bi 1 Reg | 3 | 64 (was 70) | 212624 (was 435521) |
| Bi-Bi + 1 Reg | 9 | 512 (was 530) | 63632894 (was 114684452) |
| Bi-Bi 1 Reg | 9 | 512 (was 530) | 63632894 (was 114684452) |
| Bi-Bi PP | 10 | 10 (unchanged) | 2157 (unchanged) |
| Bi-Bi Budget | 4 | 4 (unchanged) | 124 (unchanged) |

**Tests:** All 2064 tests pass. 1 pre-existing error in "Large equation compilation (<20s)" (StackOverflow in large mechanism — not related to this change).

### Agent 2: New RE/SS Expansion with G ≤ 7 Cap — COMPLETED

**Changes made:**
- `src/mechanism_enumeration.jl`:
  - Added `_compute_re_group_count(edges, eq_steps)`: union-find to compute G (number of RE groups / connected components when only RE edges merge forms)
  - Added `_find_first_isomerization(edges, adj)`: finds first isomerization edge in canonical edge order (fallback to index 1)
  - Rewrote `_ress_variants(spec, adj, forms; max_re_groups=7)`: baseline = all RE except first isomerization (always SS), then enumerate subsets of remaining edges to make SS, keeping only masks with 2 ≤ G ≤ max_re_groups
  - Added `_count_ress_variants(spec, adj, forms; max_re_groups=7)`: same logic as `_ress_variants` but returns count only (no materialization)
  - Updated `enumerate_mechanisms` FullEnumeration block: replaced closed-form O(1) count formula with brute-force `_count_ress_variants` summed over dead-end specs

- `test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl`:
  - Added `_compute_re_group_count_test(edges, eq_steps)`: test-local union-find reimplementation
  - Added `_find_first_isomerization_test(edges, forms)`: test-local isomerization detection using multi-site diff heuristic (no adjacency dict needed)
  - Rewrote `_compute_expected_n_total(spec, forms; max_re_groups=7)`: brute-force enumeration matching new `_ress_variants` logic for independent verification
  - Updated all `expected_n_total` values

**New expected counts (G ≤ 7 cap applied):**

| Spec | cat | de | total (old → new) |
|------|-----|----|--------------------|
| Uni-Uni | 1 | 1 | 7 → 3 |
| Uni-Uni 1 Reg | 1 | 8 | 808 → 338 |
| Uni-Uni 2 Regs | 1 | 64 | 3,089,511 → 1,245,541 |
| Uni-Bi 1 Reg | 3 | 64 | 212,624 → 92,136 |
| Bi-Bi + 1 Reg | 9 | 512 | 63,632,894 → 28,004,728 |
| Bi-Bi 1 Reg | 9 | 512 | 63,632,894 → 28,004,728 |
| Bi-Bi PP | 10 | 10 | 2,157 → 989 |
| Bi-Bi Budget | 4 | 4 | 124 → 60 |

**Key design details:**
- Brute-force counting for `MechanismIterator.total` iterates 2^(n-1) masks per dead-end spec, where n = edge count. For the largest case (Bi-Bi 1 Reg, 512 specs), this completes in ~10s — acceptable.
- The `max_re_groups` parameter defaults to 7 but is configurable via keyword argument. Agent 5 will expose it in the public `enumerate_mechanisms` API.
- For Uni-Uni (3 edges, 3 forms), the G cap eliminates 4 of the old 7 masks: the all-RE-except-iso baseline has G=1 (all forms in one RE group), so additional SS edges are needed. Only 3 masks give valid G ∈ [2,7].
- Test helpers use independent implementations (no `EnzymeRates._*` calls for G computation or isomerization detection).

**Tests:** All 2064 tests pass. 681 mechanism enumeration tests pass. Pre-existing StackOverflow error in large mechanism compilation unchanged.

### Agent 3: Regulator Partitioning + Dead-End Subset — COMPLETED

**Changes made:**
- `src/mechanism_enumeration.jl`:
  - Extended `MechanismSpec` with `dead_end_regulators::Vector{Symbol}` and `allosteric_regulators::Vector{Symbol}` fields
  - Added backward-compatible constructors: 2-arg `(reaction, edges)` and 4-arg `(reaction, edges, eq_steps, constraints)` both default to empty regulator vectors
  - Modified `_expand_inhibitors`: added `dead_end_regs` and `allosteric_regs` keyword arguments; only uses regulator positions matching `dead_end_regs` for dead-end expansion; tags resulting `MechanismSpec`s with partition info
  - Modified `enumerate_mechanisms`: loops over all 2^n_reg partitions of regulators into {dead-end, allosteric}; catalytic topologies computed once and reused; dead-end and RE/SS results concatenated across all partitions
  - `_ress_variants` propagates `dead_end_regulators` and `allosteric_regulators` from input spec to output specs

- `test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl`:
  - Updated `_compute_expected_dead_end_count` to sum over all 2^n_reg partitions
  - Updated all `expected_n_cat_de` and `expected_n_total` values

- `test/test_mechanism_enum_of_enz_reaction.jl`:
  - Added "Regulator partitioning" testset: verifies each dead-end spec carries valid partition info (dead_end + allosteric = all regulators), and that 2^n_reg distinct partitions are present

**New expected counts (with regulator partitioning):**

| Spec | cat | de (old → new) | total (old → new) |
|------|-----|----------------|-------------------|
| Uni-Uni | 1 | 1 → 1 | 3 → 3 |
| Uni-Uni 1 Reg | 1 | 8 → 9 | 338 → 341 |
| Uni-Uni 2 Regs | 1 | 64 → 81 | 1,245,541 → 1,246,220 |
| Uni-Bi 1 Reg | 3 | 64 → 67 | 92,136 → 92,177 |
| Bi-Bi + 1 Reg | 9 | 512 → 521 | 28,004,728 → 28,005,686 |
| Bi-Bi 1 Reg | 9 | 512 → 521 | 28,004,728 → 28,005,686 |
| Bi-Bi PP | 10 | 10 → 10 | 989 → 989 |
| Bi-Bi Budget | 4 | 4 → 4 | 60 → 60 |

**Key design details:**
- Reactions without regulators have exactly 1 partition (the empty partition); counts unchanged.
- For 1 regulator: 2 partitions. The "allosteric" partition has no dead-end expansion (just the bare catalytic topology), contributing a small additional count.
- For 2 regulators: 4 partitions. Symmetry: both single-regulator-dead-end partitions contribute equally.
- The `dead_end_regs` filter in `_expand_inhibitors` uses `s.metabolite in dead_end_regs` on regulator site positions — only sites matching dead-end regulators participate in dead-end form generation.

**Tests:** All 3297 tests pass. 1912 mechanism enumeration tests pass. Pre-existing "Rate equation too large error" failure (1 test) and segfault in large mechanism compilation unchanged.
