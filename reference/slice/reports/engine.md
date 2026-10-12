**Engine stage: the four files are built and checks (a)–(g) all pass on every plain population. kcat matches the package on 818 of 825 mechanisms, and I can explain each of the 7 differences.** I did not write `slice/engine/REPORT.md` because subagents may not write report files, so this message is the report. The printed example laws are in `slice/engine/logs/show_laws.log`.

## What I built
Code lines are counted with `proto/p5b/codelines.py`; I hand-counted `normalize.jl` and got the same 24.

| File | Code lines | Contents |
|---|---|---|
| `src/engine.jl` | 150 | `Law{T}` (`H`, `Q`, `D`, plus `sigma`, `trees` and `shift` for codegen), `bare_law`, `law` (substitutes the solve, divides exactly by H, shifts to lowest terms, retries in BigInt on overflow), `spanning_trees` with the cap, `complexity(m, mask)` |
| `src/normalize.jl` | 24 | `normalizer`: constant term, then Cleland, then affine basis, then plain coefficient |
| `src/print.jl` | 92 | `rate_equation_string`, `coefficient_form`, `printed`, `content` |
| `src/kcat.jl` | 56 | `catalytic_core` (removes dead-end branches), `kcat_terms`, `kcat` |
| **Total** | **322** | |

I added an "Engine" section to `INTERFACES.md`. The checks are in `slice/engine/` (`run_all.sh` reruns them all in about 5 minutes).

## Results
- **(a) Oracle:** 1,921 plain laws at 3 random consistent points each (all solvable records; A9 is rejected by the solve). Worst relative error 5.1e-33, none above 1e-30. Q is positive in all 1,921, no law needed BigInt, and the largest coefficient is 16.
- **(b) Package:** on its exact polynomials (857 laws × 2 points) the worst error is 1.7e-150; on its 3,428 samples, 7.9e-14. The other plain records have no package polynomials, so only (a) covers them.
- **(c) Determinism:** 5,766 rewritten inputs (step order, direction, group labels) give identical bare and final polynomials, printed strings, kcat terms and complexity. Two processes also agree on all 1,922. The rewrites canonicalize to the same mechanism, so this mainly shows the derivation is deterministic.
- **(d) Rational exponents:** 5 laws have half exponents (adv/9, adv/10, half_exponent, half_keq and the hand-built half_keq). They pass (a), and (b) where package outputs exist, and print with `sqrt`. A10g gets integer exponents from the solve.
- **(e) Units:** every printed denominator term is dimensionless and every numerator term is 1/time in all 1,919 laws that use a normalizer:

  | Normalizer | Laws (of the 806) |
  |---|---|
  | Constant term | 1,862 (757) |
  | Cleland | 42 (40) |
  | Affine basis | 15 (8) |
  | Plain coefficient (fallback) | 2: handwritten/21 and tc_pair_group |

  Both printed forms are valid Julia; evaluated, they equal the law to at most 1.9e-153 on all 1,921.
- **(f) kcat:** 818 of 825 equal the package to 1e-10.
  - 101 of those match only because dead-end branches are removed; with them kept, the limit is 0.
  - The 7 differences are expand_bibi/388, 442, 600, 609, 638, 663 and A7. In each, `_kcat_forward`'s maximum over patterns picks a sub-saturating pattern (A·B where the leading face is A²B), so it overstates kcat. My value, the leading face's ratio, is 0.02–0.91 of the package's.
  - I rebuilt the package's rule on my law and it reproduces the package in 825 of 825, which validates the comparison.
  - All 7 package errors give the clear error: the 4 hand-written "no kcat components" cases and 3 adversarial ones. Example: "kcat: the rate grows without bound with A, B, so the mechanism never saturates; fit in absolute mode (scale_k_to_kcat = nothing)".
- **(g) Cost (806 mechanisms, package re-timed on this machine):**
  - solve + law takes a median 70.7 µs against the package's 141 µs: ratio 0.58 at the median, 0.97 at p90, 1.30 at worst. The largest laws (expand_bibi/626, 642, 658) take 2.9 ms.
  - Terms: D has 30,199 against the package's 30,501 (never more). Q has 4,017 against a package numerator of 8,052.
- **Extras:**
  - `complexity` equals the package's `_eq_complexity` on 3,853 of 3,853 records.
  - The cap triggers on a complete 6-node segment graph (7,776 > 5,000) and not on 5 nodes (625).
  - Inactive-state laws derive on all 2,440 allosteric records with a non-empty mask. Whether the inactive state turns over agrees with the package on 2,960 of 2,995; the 35 differences are all P5a hand-built cases.

## Printed examples
RE random bi-bi (init_bibi/55):
```
v = E_total * (A * B - P * Q / Keq) * k_EAB_to_EPQ / (K_EA_to_E_A * K_EB_to_E_B) / (1 + A / K_EA_to_E_A + B / K_EB_to_E_B + P / K_EP_to_E_P + Q / K_EQ_to_E_Q + A * B / (K_EA_to_E_A * K_EB_to_E_B) + P * Q / (K_EP_to_E_P * K_EQ_to_E_Q))
```
- **Steady-state uni-uni:** prints `N = k_ES_to_E_P + k_ES_to_E_S` on its own line, then `… / (1 + S * k_E_S_to_ES / N + …)`.
- **Ordered bi-bi:** coefficients come out in Cleland's form, e.g. `d_A = k_E_A_to_EA / k_EA_to_E_A`.
- **Classic ping-pong:** Cleland's N, and no leading 1.
- **half_keq:** `k_F_to_E = … / (… * sqrt(Keq))`.
- **Squared concentrations** are named with superscripts: `d_A²_B * A^2 * B`.

## Deviations and choices to confirm
1. **Printed laws are in the independents.** Dependents are substituted, and the Haldane and Wegscheider lines are documentation only. Today's strings use the dependents' names inside v.
2. **How N prints.** The normalizer's integer and monomial content is divided into each term; any remaining polynomial factor prints as its own `N = …` line. A coefficient whose polynomial part equals N cancels to a monomial.
3. **kcat semantics:**
   - Plain-mechanism modifiers are treated as corners, which matches the package on the activator specs.
   - Substrates grow first, then the corner's modifiers.
   - Corners are the modifier sets that the law's terms hold. This gives the same maximum as all 2^k corners and keeps A12 to 211 corners.
   - An error is raised if the leading face is not one ratio; this never happened on the data.
4. **The cap counts simple trees × segments** (at most 5,000). `complexity` keeps today's count, which includes parallel edges.

## My mistakes, all fixed
- A vcat bug in `normalize.jl`.
- A print check that ignored E_total.
- A variable clash in a check script.
- Enumerating all 2^20 kcat corners on A12, which ran for 500 s.
- A Julia process that ignored SIGTERM held the lock for about 20 minutes; I then used `timeout -s KILL` and killed by PID.
- I once ran `pgrep -f` with a pattern that also matched my own shell, so it killed my shell along with my Julia process; nothing else was touched.

No Julia process of mine is left running.

## For later stages
- **Allosteric:** call `law(m, sol, inactive_mask(m); inactive = true)`; both conformations share `vars`. 555 allosteric records have an empty inactive mask, and you must handle them (Q_I = 0, D_I = 1).
- **Codegen:** the law is entirely in independents; the dependents are in `sol`. Keq appears in Q in 577 laws.
- **Integration:**
  - `catalytic_core` calls enum.jl's `unbalanced_blocks`; move it to an earlier file.
  - The printed local `N` could collide with a metabolite named N; nothing checks for that.

Everything is in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/`:
- `src/engine.jl`
- `src/normalize.jl`
- `src/print.jl`
- `src/kcat.jl`
- `INTERFACES.md`
- `engine/run_all.sh`
- `engine/logs/`