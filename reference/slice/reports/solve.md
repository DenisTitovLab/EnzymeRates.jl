I built the solve stage (`thermo.jl`, `api.jl`), and checks (a) to (f) all pass on every population. I did not write `slice/solve/REPORT.md`: the harness forbids report files, as it did for the core stage, so this message is the report. The detailed numbers are in `slice/solve/logs/`.

**Files** (code lines from `proto/p5b/codelines.py`; I hand-counted `api.jl` and got the same 23):

| File | Code lines | Contents |
|---|---|---|
| `slice/src/thermo.jl` | 182 | Solve types, inactive mask, priority, cycles, the row builder `constraint_system`, `rref!`, `solve` |
| `slice/src/api.jl` | 23 | `column_names`, `fitted_params`, `parameters`, `all_parameters` |
| **Total** | **205** | The critic projected 150–200, plus this stage's roughly 14-line inactive mask |

`INTERFACES.md` has a new "Solve" section. The check scripts are in `slice/solve/`, and `run_all.sh` reruns them all under the lock in about 17 minutes.

**Results**
- **(a) Relation space:** on all 806 plain mechanisms my relations span exactly the package's (dependents, monomials, name ties and folds, compared by exact rank). The same holds on 51 more package-built adversarial mechanisms (857/857). That includes `iso_group_opposite_chem` and `iso_group_mixed_tier`, where P2's alignment failed (F1).
  - My dependent set differs from the package's in 564 of 806.
  - Package samples rebuilt through `all_parameters`: 3,428 samples, worst error 4.0e-154.
- **(b) Detailed balance:** 3,842 draws on all 1,921 solvable plain records across 16 sets. Worst |v(Keq)|/|v(Keq/10)| is 4.5e-105, and none is above 1e-30. One record is rejected, as the package rejects it: `p1_adv_pkg` A9 slip cycle, whose cycle is no multiple of the net reaction.
- **(c) Determinism:** the solve differs in 0 cases across:
  - 24,580 rewritten inputs (step order, direction, group labels);
  - 14,748 row changes that keep the row space;
  - 12,608 same-model splits (4,313 plain, 8,295 allosteric);
  - 45 merges along folds.

  Labels change as expected. **But the fitted names change in 57 of those splits:** today's naming rule names a group after a different member than the one the solve keeps (e.g. init_bibi_pp/252).
- **(d) Rational exponents:**
  - `half_keq` gives `Keq^(-1/2)`, the same as the package.
  - `half_exponent` gives ±1/2 exponents.
  - A10g gets an integer-exponent dependent from my priority, where the package needs half powers.
  - Across all records, 5 have non-integer exponents (3 with Keq^(p/q)); the largest denominator is 2.
- **(e) Allosteric, the 1,993:**
  - Phantom names: 0 (package 142 in 115 mechanisms). Fold hazards: 0 (package: 9 opposite-direction folds, plus 30 mechanisms with one-state folds).
  - Fitted length against the package's `fitted_params`: 1,844 equal, 123 one fewer, 26 two fewer. These are P5b's numbers exactly. Every difference is a package phantom name (110), a package fold hazard (38) or both (1).
  - Same model where counts are equal: the package's samples satisfy my relations to 3.0e-154. Where counts differ, the residuals are 0.91 to 7.2e6.
- **(e) Allosteric, the 848:** 791 equal, 47 one fewer, 10 two fewer. All 57 differences are package phantom names (67 names).
- **(e) Both sets:**
  - Two-state detailed balance in P5a's oracle: worst 1.6e-148 and 7.8e-148.
  - My inactive mask equals the package's pruning in all 2,841.
  - Giving each `:OnlyA` member its own inactive class versus tying them changed the fitted length in 0 of 2,841, so the per-step classes matter only for the cone test. Island cycles changed no count either.
  - Identifiable rank, by finite differences at Λ = 1e60 on 322 mechanisms (report only):
    - Where my count is below the package's, count equals rank in 195 of 206. The 10 n = 1 deficits are the known one-subunit L phantom.
    - Across specs, handwritten and allo_uniuni, 14 of 56 are rank-deficient (MWC dimer, collapse fixtures) and have equal counts. By (e)'s same-model check the package's model is the same, so I infer the package has those deficits too.
- **(f) Arithmetic and cost:**
  - Checked Int64 throughout: 0 overflows. The largest numerator or denominator written is 2, and the result equals the BigInt RREF in 4,916 of 4,916.
  - Systems are at most 34 rows × 84 columns.
  - A solve takes a median 13.5 µs (p90 36, max 185). `fitted_params` takes 37 µs and `all_parameters` 38 µs; rendering names is most of that.
  - Separate RE-only cycle rows are redundant: adding 7,150 of them changed no rank.

**Deviations from the decisions**
1. **"Against" follows the reaction direction, not "slot 2".** For a binding of a product, whose physical-forward direction is taking the product up, the binding rate is "against". Under the literal slot-2 rule the Haldane dependent is a catalytic rate constant in 228 of 806 mechanisms (Uni-Uni: `k_ES_to_E_P`); under mine, in 0. P2's design never eliminates a forward constant. Please confirm this choice.
2. **No separate RE-only cycle rows or `:EqualAI` tie rows.** The fundamental cycles already span the RE-only cycles. An `:EqualAI` edge has no inactive column; it reads its active one.
3. **Ties are broken by the oriented edge's content, not by `step_text`.** It has the same invariances and needs no rendering.
4. **Parameter order is content order:** active constants by edge and slot, then inactive ones, then regulator constants, then L. The package sorts names alphabetically.
5. **Inactive cycle rows cover the whole graph.** Columns outside the law are eliminated first, so islands' cycles still bind, which is what the limit semantics need. This changed no count here.
6. **`inactive_mask` lives in `thermo.jl`.** The flip form is the first free form by content.

**My mistakes, all fixed**
- `isassigned` with an `Int32` index throws and catches an exception internally, which made a solve take 1.4 ms instead of 12 µs.
- My first "against" rule was literal slot 2.
- My first allosteric sample check skipped the `:EqualAI` inactive tie, so the unequal-count mechanisms passed trivially.
- At Λ = 1e40, an O(1/Λ) artifact made expand_bibi/94 look full-rank (kept pivot 5.4e-35, just above the 1e-35 cutoff); at Λ = 1e60 it has the deficit of 1.
- A `@printf` string-concatenation bug in a check script.

No Julia process of mine is left running.

**What later stages must know**
- **Engine:** edge `k`, slot `s` in a state is column `sol.index[Col(k, s, inactive)]`. Substitute `sol.deps` (terms, Keq exponent) and name independents with `column_names`.
- **Admission:**
  - `constraint_system(m)` holds the full system, including each `:OnlyA` member's inactive columns. An RE `:OnlyA` chemistry edge has one K column with no sign constraint.
  - `inactive_mask` does not reject a second free form outside the flip form's RE segment; admission must.
- **Names:** the 57 name changes under same-model splits need a decision: name a group after the member the solve keeps, or keep today's names and accept the changes.
- `enum.jl` (a parallel stage) already calls `solve` and `inactive_mask`.

Files are in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/`:
- `src/thermo.jl`
- `src/api.jl`
- `INTERFACES.md`
- `solve/run_all.sh`
- `solve/logs/`