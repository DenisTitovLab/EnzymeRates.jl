Denis, I could not refute the first-passage argument. With no product present, the textbook mechanism never runs faster than v_sat, whatever the rate constants, and I found a second proof that agrees with it. Your intuition is right once a product is present. Q (or P) can then push the rate above v_sat. The maximum sits at saturating A and B, zero P and an intermediate Q. The substrates never peak; that needs a substrate binding a product-bound form, which the textbook mechanism doesn't have. The earlier probe held P = Q = 0, where no parameter set can do it. Its 91/300 were peaks along fixed-B slices, which stay below v_sat.

## 1. Solver and cross-check
- **Solver:** an independent steady state over the enzyme forms, using an elimination that never subtracts and so stays accurate for very stiff rates. Float64 for the search, 256-bit BigFloat for checking. Each step gets a forward constant and each form a free energy, with Keq fixed, so every Haldane and Wegscheider condition holds automatically. Every forward and reverse constant is kept in 1e-6..1e6.
- **Against `rate_equation`:** the package accepted the textbook and all 8 variants. Over 40 parameter sets × 4 product settings per mechanism (half spread over ±6 decades), the worst relative difference was 3.6e-15 in Reduced mode and 2.1e-15 in Full mode. The free-parameter counts match (15 for the textbook). At the optima in §3 and §5, the package and the solver agree to 10 digits.
- **A numerical trap you should know about:** with both products present, the flux through a single step can cancel badly in Float64 (relative error up to 85 on random points). My first run reported R = 1.0003 for variant (g), which BigFloat puts at 1 − 1e-11. The solver now picks the flux expression with the least cancellation (worst error 2.6e-5), and every R below is recomputed in BigFloat.

## 2. Search
- **Optimizer:** CMA-ES with restarts, over rate constants and log-concentrations together. Bounds: A, B in 1e-8..1e8; P, Q either 0 or 1e-8..1e8. Keq = 1e-3, 1 and 1e3 (also 1e±6 for d1, d2).
- **Checks on each best point:** a polish that keeps the constraints, a grid search over concentrations at the best constants, and a BigFloat recheck.
- **Extra pressure:** 20 independent seeds for every case that stayed at or below 1, and boxes capped at A, B ≤ 1e2 and ≤ 1e4 so the saturated corner is out of reach.
- **Also run:** 7 targeted textbook sets (slow E(A)+B with fast E(B)+A and fast E(A,B) dissociation, tight A binding, near-equilibrium catalysis, and others) and a census of 5000 random feasible textbook sets.
- **v_sat:** computed at A, B = 1e48 along three directions, and at 1e80. For textbook, a, b, c and g it is the same in every direction.

## 3. Results
R = sup v / v_sat, with v_sat at A, B → ∞ and P = Q = 0, checked in BigFloat.

| mechanism | P = Q = 0 | products allowed |
|---|---|---|
| textbook | 1 − 4e-18, reached only as A, B → 1e8 | **1.7155e11** at A = B = 1e8, P = 0, Q = 1.1e4 |
| (a) E(A,B) ⇌ E + P + Q | ≤ 1 | ≤ 1 |
| (b) rapid-equilibrium release | ≤ 1 | ≤ 1 |
| (c) E(A,B) ⇌ Estar(A,B) ⇌ E(P,Q) | ≤ 1 | **1.339e11** at A = B = 1e8, Q = 2.3e3 |
| (d1) X1 from E(A)+B, X2 from E(B)+A, both ⇌ E(P,Q) | ≤ 1 | **1.7155e11** |
| (d2) X1, X2 each ⇌ E + P + Q | ≤ 1, against the larger directional v_sat | ≤ 1 |
| (e) dead end E(Q)+A ⇌ E(A,Q) | v_sat = 0 for every set, so a peak always | same |
| (f) E(Q)+A ⇌ E(A,Q) ⇌ E(A)+Q | **3.33e11** at A ≈ 1e4, B = 1e8 | same |
| (g) random binding, ordered release (package fixture 11) | ≤ 1 | ≤ 1 |

- "≤ 1" means the best R lay between 1 − 6e-16 and 1 − 2e-19, always at the saturated corner. No seed exceeded 1 + 1e-9 in any configuration.
- For (d2), v_sat depends on the direction of saturation because the first binding fixes the route (spread up to 7.5%).
- The values around 1e11 are set by the box on rate constants; they are not intrinsic.

## 4. Why no substrate peak at zero product (second proof)
At P = Q = 0, release from E(P,Q) is irreversible. All forms except E(P) and E(Q) form a reversible network whose only exit is E(P,Q), at rate k_rel = k6 + k7. Split each turnover into three parts:
- τ_in: from E to the first arrival at E(P,Q). It is at least 1/k5, because E(P,Q) is entered only from E(A,B).
- t_S: from E(P,Q) to release, back-excursions included. A renewal argument gives t_S = Σ π_j / (π_E(P,Q) · k_rel), where π is the equilibrium distribution of that network at the current A, B. The terms for E, E(A) and E(B) fall as A and B rise and vanish at saturation; π_E(A,B)/π_E(P,Q) = 1/K5 is fixed.
- τ_rel: the dwell in E(P) or E(Q), which does not depend on A or B.

Saturation reaches both lower bounds, so v(A, B) < v_sat at every finite A, B and for every rate-constant set.
- **(d1):** the same argument with max(kc1, kc2) in place of k5.
- **(d2):** each turnover has exactly one completion, whose rate never exceeds max kc.
- **(a), (b) with products:** the steady state solves in closed form, J = (k5 α f_E − k−5 β f_EPQ) / [(1 + c1 f_E)(β + k5) + (1 + c2)(α f_E + k−5 f_EPQ)]. Here α, β, c1 and c2 depend only on A and B. Product lowers f_E and raises f_EPQ, so it can only slow the enzyme.

## 5. What a product changes: rescue of a slow release branch
Steps: 5 is E(A,B)→E(P,Q), 6 is E(P,Q)→E(Q)+P, 7 is E(P,Q)→E(P)+Q, 8 is E(Q)→E+Q, 9 is E(P)→E+P. At P = 0 and saturating A, B, the turnover time as a function of Q is:

T(Q) = 1/k5 + (1 + k−5/k5 + k6/k8 + k7·q/k9) / (k6 + k7·q), where q = k9 / (k9 + k−7·Q)

- **Condition:** Q speeds the enzyme exactly when k9 < k6 / (1 + k−5/k5 + k6/k8). In words: P leaves E(P) more slowly than the enzyme can finish through the E(Q) branch.
- **Mechanism:** Q rebinds E(P), returning it to E(P,Q), which then exits as E(Q) + P → E + Q. The gain is T(0)/T(∞).
- **Why Q is intermediate:** at P = 0, E + Q ⇌ E(Q) is a dead end that competes with A and B for free enzyme, so the rate falls again at high Q.
- **P works the same way** when step 8 is the slow one.
- **Substrates never peak:** at P = 0 and any fixed Q, v rises monotonically to its A, B → ∞ limit (same gateway argument as §4).

**Example** (Keq = 1, all constants 1 except E(P) ⇌ E + P at 1e-3 both ways): v_sat = 1.990e-3. At A = B = 1e8 and P = 0:

| Q | v / v_sat |
|---|---|
| 1e-2 | 5.8 |
| 1 | 100.6 |
| 1e4 | 125.62 (maximum) |
| 1e8 | 111.7 |

The formula gives 502.5/4 = 125.625; the package's `rate_equation` gives 125.6202969.

**Box optimum** (Keq = 1): E(P) → E + P at the 1e-6 floor (reverse 1.1e-6), E(P,Q) → E(P) + Q at 1e6 with Q rebinding at 1e6, E(P,Q) → E(Q) + P at 7.1e5, E(Q) → E + Q at 1e6, catalysis 1e6 (reverse 90). v_sat is 1.71e-6 because 59% of turnovers wait in E(P). At Q = 1.1e4 the rate is 2.93e5. The formula predicts 1.71566e11 at these constants; the search reached 1.71547e11, the gap coming from A and B stopping at 1e8.

## 6. Which structures allow a peak
- **No peak above v_sat at zero product, for any constants:** textbook, (a), (b), (c), (d1), (d2), (g).
- **Peak in a substrate at zero product:** needs a substrate binding a form that exists only after product release.
  - (e): v → 0 as A → ∞ for every set. With all constants 1, v peaks at A ≈ 1 and then falls as 1/A.
  - (f): with all constants 1 except E(A,Q) ⇌ E(A) + Q at 1e-3, v peaks at 96.06 × v_sat at A = 2.2, B = 1e8 (the package agrees).
- **Peak in a product:** needs steady-state random release from E(P,Q) (textbook, c, d1). None with rapid-equilibrium release (b), lumped release (a, d2) or ordered release (g).
- **The earlier 91/300:** in my census, 3665 of 15000 fixed-B slices (24%) peak in A, up to 2.2e6 × v(A = 1e8) at that B. All of them stay below v_sat. The census's largest sup v/v_sat is 1 to 15 digits; the targeted guesses reach at most 0.99999999998.

## 7. Not checked
- Topologies beyond these nine: ordered binding, ping-pong, other abortive complexes (E(P)+B, E(A)+P), oligomeric or allosteric forms. Rapid-equilibrium substrate binding was not run (the §4 bound covers it as a limit).
- Ranges beyond rate constants 1e-6..1e6, concentrations 1e-8..1e8, and Keq 1e-3..1e3 (1e±6 for d1, d2).
- (g) and (d2) with both products: search only, no proof. The claim that A, B never peak when both P and Q are present is also search only.
- Whether the package's fitter or identify pipeline ever proposes, or can fit, product-activated parameter sets.

Total Julia compute was about 9 minutes; no package files were touched.

Files are in /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/peaksearch/:
- solver.jl
- mechs.jl
- pkgdefs.jl
- search_lib.jl
- crosscheck.jl, crosscheck.log
- crosscheck_best.jl, crosscheck_best.log
- search.jl, full.log, d1focus.log, results_full.jls
- search2.jl, multiseed.log, results_multiseed.jls
- guesses.jl, guesses.log
- profiles.jl, profiles.log
- selfcheck_re.jl, selfcheck_re.log