Denis, the search agent's main conclusions hold. Three of its statements are wrong or stated too broadly.

I rebuilt every claimed peak with my own exact rational-arithmetic solver, and each R matches to 10 digits. Every parameter set is thermodynamically consistent and inside the 1e-6..1e6 box. My own 22.5-minute search of the textbook mechanism at zero product found no rate above v_sat. The first-passage identities hold exactly.

The three corrections:
- **"Rises monotonically" with Q present** is false along fixed-B slices.
- **"No peak" for d1 and d2** holds only against the best saturation direction. At a fixed A/B ratio both peak far above that ratio's own saturated rate.
- **The (f) v_sat** is 6e-5 too high.

## Verdicts

**1. Textbook, P = Q = 0: no rate constants give v > v_sat — CONFIRMED**
- **Proof check:** I checked the first-passage argument exactly on 400 random consistent sets (kf, K, Keq in 1e-6..1e6; A, B in 1e-8..1e8). Each identity held in all 400 sets:
  - mean time in E(A,B) before the first catalysis = 1/k5;
  - t_S from the fundamental matrix = Σπ_j/(π_EPQ·k_rel);
  - 1/v = τ_in + t_S + τ_rel;
  - v < v_sat (largest ratio 0.99999999998963).
- **Search:** 4.2e8 joint evaluations plus 753,000 sampled rate-constant sets, each with its own sweep over A and B.
  - No point reached R > 1 + 1e-12. The best Float64 R was 1 + 1.3e-15, which is rounding.
  - The 25 best points, rechecked exactly, all have 1 − R between 6.1e-20 and 6.2e-17.
  - With k in 1e-3..1e3 and A, B ≤ 1, the best was 0.999998999979 at A = B = 1.
- **Rays:** along any fixed ratio A/B the rate rises monotonically.
  - Argument: τ_in = H_E + (k−3/k5)·H_EA + (k−4/k5)·H_EB + 1/k5, where H_E = N(t)/D(t) with N falling and D rising in t; t_S also falls along the ray.
  - The ray search found no decrease (largest drop 1 + 1.6e-15).
- **Slices:** substrate peaks occur only at fixed B. The fall from the peak reaches 1.7e11×, but the peaks stay below v_sat (largest slice maximum = v_sat × (1 − 4.8e-13)).

**2. Textbook with a product: R = 1.7155e11 at A = B = 1e8, P = 0, Q = 1.1e4 — CONFIRMED**
- **Exact value:** R = 1.7154678312e11 (v = 2.9285121312e5, v_sat = 1.7071215664e-6). All 9 T rows (Q, P, PQ × Keq 1e-3, 1, 1e3) reproduce.
- **Thermodynamics:** the three cycle products deviate by less than 1e-150 in log10. The checker does detect errors: changing any one reverse constant by 1% gives a deviation of 0.00432.
- **Package:** `rate_equation` matches my exact solver within 2.2e-16 at 10 points, in both Full mode (18 constants) and Reduced mode (15 + Keq).
- **Set by the box:** my DE in the same box reached 1.7156575e11; the analytic supremum is 0.171573 × kmax/kmin = 1.71573e11. With k in 1e-9..1e9 it reached 1.7157280e17.
- **Wording (PARTIAL):** "maximum at intermediate Q" depends on A and B being finite.
  - At A = B = 1e8, R peaks near Q = 1.1e4 and falls to 1.2689e11 at Q = 1e8.
  - At A = B = 1e30, R is still rising at Q = 1e8 (1.7156644e11), matching the closed-form v_sat(Q).
- **§5 example** (all k = 1, E(P) ⇌ E + P at 1e-3 both ways):
  - v_sat = 2/1005 exactly.
  - At A = B = 1e8, v/v_sat = 5.7703 (Q = 1e-2), 100.6004 (Q = 1), 125.6203 (Q = 1e4), 111.6667 (Q = 1e8).
  - The true maximum is 125.62056 at Q ≈ 1.4e4; the limit is 1005/8 = 125.625.
  - The T(Q) closed form matches the exact limit to 5.8e-134 on 100 random sets. The condition k9 < k6/(1 + k−5/k5 + k6/k8) checks out algebraically.

**3. (c) 1.339e11 — CONFIRMED.** Exact 1.3386629083e11; all 9 rows reproduce.

**4. (d1) 1.7155e11 — CONFIRMED.** Exact 1.7154678052e11; all 9 rows reproduce.

**5. (e) v_sat = 0 for every set — CONFIRMED.**
- All-ones set: v·A = 2.000000 from A = 1e20 to 1e200. The rate is 0.1667 at A = 1 and 0.1212 at A = 10 (B = 1e8).
- 20 random sets: v·A at 1e100 and at 1e200 agree to 3.8e-84.

**6. (f) — CONFIRMED, with a correction.**
- Exact R = 3.3327073865e11; the agent reported 3.3325055390e11.
- The agent's v_sat is the maximum over three directions at 1e48. One of them, B ≫ A, uses A = 1e24, which does not yet saturate this set, so its v_sat is 6.1e-5 too high. At 1e100–1e300 every direction gives 1.000000e-6.
- All 12 f rows carry the same 5e-5 to 2e-4 offset.
- Example 96.06484288 at A = 2.24, B = 1e8: CONFIRMED.

**7. (a), (d2), (g) ≤ 1 — CONFIRMED at the agent's best points.**
- Exact 1 − R ranges from 4e-18 to 4e-14.
- The g "numerical trap" is CONFIRMED: with both products, the Float64 R of 1.0003 is exactly 1 − 1.3e-11, 1 − 3.45e-8 and 1 − 1.2e-13.
- My DE for g with one product (9e7 evaluations each for Q-only and P-only) found a maximum R of 1 + 1.1e-15, which is rounding. g shows no product activation.

**8. §6 "no peak at zero product" for d1 and d2 — PARTIAL.**
- True against the best saturation direction:
  - d2: largest v/max(kc1, kc2) = 1 − 9e-13;
  - d1: all 200 exact sets satisfy time in X1 + X2 ≥ 1/max(kc), with the largest ratio 1 − 4e-11.
- False against the saturated rate at the same A/B ratio (exact):
  - **d2:** v/v_sat = 8.22e10 with kc1 = 1e-6, kc2 = 9.0e5, k1 = 698, k2 = 1e6, at A = 5.7e3, B = 2.5.
  - **d1:** v/v_sat = 6.0e10 with kc1 = 1e6, kc2 = 1e-6, k1 = 1e6, k2 = 165. Along the ray the ratio rises from 3.1 (t = 1e-6) to 6.0e10 (t = 1), then falls to 1.28 (t = 1e12).
- So diluting A and B at a fixed ratio gives a peak whenever the binding order picks ternary complexes with different kcat. No product-bound form is involved, which contradicts the agent's general claim that a substrate peak needs one.

**9. "At P = 0 and any fixed Q, v rises monotonically to its A, B → ∞ limit" — PARTIAL.**
- **The bound holds:**
  - exact v < v_sat(Q) in 400 of 400 random sets;
  - a DE search in two boxes (5.5e7 evaluations) found best points with exact 1 − R of 1.6e-17 and 6.7e-18.
- **"Monotonically" is refuted:** exact example at Q = 1, P = 0, B = 1, all constants in 1e-4..1e4, Keq = 100.
  - kf = [5.82e-4, 1e4, 1e-4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4]
  - kr = [1e-4, 2170, 1e-4, 7911, 2487, 1e4, 2341, 2341, 1e4]
  - v = 907.9 at A = 1, 285.3 at A = 10, and 1.34e-4 at A = 1e8 — a 6.8e6-fold fall, all below v_sat(Q).

**10. v_sat is the true limit — CONFIRMED, except (f) (item 6).**
- T, a, c, g: the 1e48 value and the exact 1e200 value differ by at most 5e-55, and all directions agree to 1e-107.
- d1: differs by at most 1.2e-31; the spread across directions is up to 2.2e-8 at its optima.
- d2: the agent's 7.5% direction spread reproduces.
- The textbook closed form matches the exact limit to 4.2e-86 (100 sets × 3 directions).

## How I checked
- **Solver:** exact rational steady state. Net P, net Q, A uptake and B uptake must agree exactly, and did.
- **Rate evaluator:** King–Altman spanning trees (48 for T), with no subtraction when one product is absent. It agrees with exact to 2.1e-15.
- **Thermodynamics:** a cycle-product checker that works from the rate constants alone, independent of the agent's free-energy parameterization.
- **Main search (1352 s), using a different method from the agent's:**
  - differential evolution with about 8e4 restarts on stagnation, plus Nelder–Mead polish;
  - my own parameterization: forward constants, spanning-tree reverse constants, chords from Wegscheider/Haldane, Keq free;
  - boxes: k in 1e-12..1e12, 1e-9..1e9 or 1e-3..1e3; A, B ≤ 1e12, ≤ 1, or 1e-6..1; Keq free or pinned at 1e±6;
  - the agent's 7 guesses as seeds;
  - bilevel sampling over the 753,000 rate-constant sets.
- **Secondary searches (695 s):** Q > 0, rays, slices, g with one product, d1/d2 at fixed ratio.
- **Total Julia:** about 35 minutes, one process at a time under the lock. No package files were touched.

## Not checked
- (b) was not rebuilt, and the (a)/(b) closed-form J with products was not verified.
- No independent search with both P and Q > 0; I only rechecked the agent's P,Q optima exactly.
- Rays with Q > 0, and rays for c and g: argued, not searched.
- The earlier probe's 91/300 (no files).
- Topologies beyond the nine, including dead-end substrate complexes; the package's fitting and enumeration.
- Rate constants outside 1e-12..1e12.

Everything is in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/peaksearch/verify/`:
- `vlib.jl`, `vopt.jl`
- `check_points.jl` / `check_points.log`
- `fp.jl` / `fp.log`
- `vsearch.jl` / `vsearch.log` / `vsearch_results.jls`
- `vextra.jl` / `vextra.log`
- `d1ray.jl`, `pkgcheck.jl`
- `followup.jl` / `followup.log`
- `runall.jl`