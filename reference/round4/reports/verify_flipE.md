**Verdict:** three of the four claims I tested hold. The fourth, the justification for reading B, fails for one class of binding. The probe also missed one thing that changes "adopt flip-at-E now": some of the newly admitted laws are exact but have zero rate at zero products, and the kcat rescale then throws.

All runs used my own scripts in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/flipE/verify/`. "Prototype" is the probe's worktree `wt`; "PR #86" is `wt_base`.

## Claim 1 (Task 2): flip-at-E laws are exact — CONFIRMED

I wrote an oracle that does not use the package's internals:
- It maps parameter names to steps from the naming convention.
- It solves the Haldane relation itself and asserts every other constant is known.
- The flip is a bridge edge at finite speed, so only its ratio L matters.
- Rapid-equilibrium steps use Λ = 1e40 with extrapolation from Λ and 2Λ.

Coverage: 12 shapes, n = 1, 2 or 3, 20 (shape, n) pairs and 208 law-vs-oracle points (`logs/exact_v_proto.log`). Three shapes are new: three free forms on one cycle, a steady-state dead-end `D <--> E`, and one rapid-equilibrium segment rooted at `Ca`.
- Worst error is 3.2e-69. For all-steady-state shapes it is ≤ 3.4e-136. The oracle's rate at equilibrium is ≤ 1.4e-111 of a typical rate.
- In all 20 pairs, flipping at another free form gives a different rate: ≥ 6.3e-4, and ≥ 0.089 in 19 of 20. So the test can tell the flip form apart.
- PR #86 rejects 15 of the 20 pairs at construction (`logs/exact_v_base.log`).

Examples:
- `E+S<-->E(S); E(S)<-->Eprime(P); Eprime+P<-->Eprime(P)` (all EqualAI) plus `Eprime<-->E :: NonequalAI`, n = 3: error 3.4e-136. A flip at Eprime would differ by 0.21.
- `E+S⇌E(S); E(S)<-->F(P); F+P⇌F(P)` (EqualAI), `F<-->G :: NonequalAI`, `G<-->E :: EqualAI`, n = 2: error 2.5e-80. Flips at G and F would differ by 6.3e-4 and 2.3e-3.
- `E+S⇌E(S) :: EqualAI; E(S)<-->E(P) :: NonequalAI; E+P⇌E(P) :: EqualAI; Ca⇌E :: NonequalAI`, n = 2: error 5.4e-80.

## Claim 2 (Task 4): nothing changes for enumerated mechanisms — CONFIRMED, on a wider set than the probe's

- **Census.** I checked 9,308 distinct enumerated mechanisms from reactions the probe did not run:
  - ping-pong bi-bi: 3,427
  - ter-bi: 5,636 (all 4,885 init; gen1 from 30 seeds)
  - uni-bi: 167
  - two regulator/inhibitor uni-uni reactions: 36 and 42

  Every one has exactly one free form, and every form has conformation `:E`. 1,160 of them carry residual forms (`logs/census_proto.log`, `census_proto_terbi.log`, `census_base.log`).
- **Law comparison.** 223 sampled allosteric laws (60 ping-pong, 60 ter-bi, 53 uni-bi, 25 + 25 regulator) give identical rows in both checkouts: fitted params, rate-string hash, three rates, kcat (`data/census_{proto,base}.tsv`). One ter-bi row fails identically in both at the 5,000-term guard.
- **Probe's own data.** Its two TSVs are byte-identical.
- **Code reading.** `_segment_root` ranks by bound count, then residual, then name, so a lone free E roots its segment and its alpha is 1. `_reduce_conc_lowest_terms` keeps the free-enzyme weight out of the minimum, so num/den cannot change.
- **Examples** (`logs/examples.log`, identical in PR #86):
  - Ping-pong `E+A⇌EA :: OnlyA; EP_res+A⇌EAP_res :: OnlyA; E+Q⇌EQ; EB_res+Q⇌EBQ_res; EA+P⇌EAP; E_res+P⇌EP_res; EA<-->EP_res :: OnlyA; EB_res⇌EQ :: OnlyA; EQ+B⇌EBQ; E_res+B⇌EB_res` (unmarked steps are EqualAI). The inactive state keeps E, EQ, EBQ.
  - Ping-pong `E+A⇌EP_res :: OnlyA; E+Q⇌EQ :: OnlyA; E_res+P<-->EP_res :: EqualAI; E_res+B<-->EQ :: OnlyA`. The inactive state is fully pruned.

## Claim 3 (Task 1): PR #86 flips at D, and the alpha factor moves the flip to E — CONFIRMED, but PR #86's law is not wrong

- PR #86's law at L equals the prototype's law at L·(K_I/K_A)^n to all 16 printed digits, in 16 of 16 comparisons: 2 shapes × n = 1, 2 × 2 draws × 2 points (`logs/reparam_{base,proto}.log`).
  - Example: `E+S<-->E(S); E(S)<-->E(P) :: NonequalAI; E+P<-->E(P); D⇌E :: NonequalAI`, n = 2, S = 0.7, P = 1.3. PR #86 gives 8.443163169904860e-03 at L; the prototype gives the same value at L·(K_I/K_A)², and 1.2303e-02 at L.
- So for free forms joined to E in rapid equilibrium, the two laws are the same model family. Only the meaning of L differs, and in PR #86 that meaning depends on whether the other form's name sorts before "E".
- Flip-at-E fixes laws only for steady-state-joined free forms, which PR #86 rejects today.

## Claim 4 (Task 3 and recommendation 2B): `:OnlyA` islands — PARTIAL

**Confirmed for R4** (`logs/block.log`). R4 is `E+S⇌E(S); E(S)<-->E(P) :: NonequalAI; E+P⇌E(P); F⇌E :: OnlyA`.
- I kept `F⇌E` in the inactive state at rates ε·(ratio, 1).
- The result equals the EqualAI twin to ≤ 2.5e-80 at ε = 1, 1e-6 and 1e-14. It does not depend on ε at all, so under reading A, `:OnlyA` on a dead-end step means the same as `:EqualAI`.
- The deletion law differs by 1.4e-2 to 5.2e-1 over 8 rows (example: n = 1, S = 3, P = 0 gives 0.1723 for deletion against 0.0820 for the twin).

**Confirmed for R2.** R2 is `E+S⇌E(S) :: OnlyA; E(S)<-->E(P); F+P⇌E(P); F<-->E :: OnlyA`.
- At P = 1 the kinetic-block network approaches deletion as ε falls (a = 1: 1.6e-5 at ε = 1e-6, 1.6e-11 at 1e-12).
- At P = 0 the gap stays constant as ε shrinks: 0.76 for a = 0.5 and 1.21 for a = 1. For a = 2 it grows like 1/ε.

**Refuted: "B is the reading `:OnlyA` bindings already have (K_I → ∞)."** Counterexample: `E+S⇌E(S); E(S)<-->E(P); E+P⇌E(P)` (EqualAI) plus `F+P<-->E(P) :: OnlyA`, in steady-state and rapid-equilibrium variants (`logs/away_counter.log`).
- The prototype accepts it and prunes F from the inactive state.
- Deletion equals the "away" limit: error ∝ ε, from 3.7e-9–2.5e-7 at ε = 1e-6 down to 3.7e-15–2.5e-13 at 1e-12.
- The K_I → ∞ limit (P uptake vanishes, release kept) drives the rate to 0, with all inactive enzyme in F (relative error 1.0).
- Kinetic block is off by 3.7e-3 to 2.0e-1.
- So for a binding whose bound side lies in free E's component, reading B means K_I → 0, the opposite of the binding convention.
- At exactly P = 0 every reading's network is singular (F's only exit binds P, in both conformations). The package's value there is the P → 0 limit.

## What the probe missed

1. **Degenerate laws are newly admitted, and kcat throws on them.** Shapes such as `E+S<-->E(S); E(S)<-->Eprime(P); Eprime+P<-->Eprime(P)` (EqualAI) with `Eprime<-->E :: OnlyA` (also with RE bindings, a catalysis-NonequalAI variant, the bi-bi analogue, and three free forms with `G<-->E :: OnlyA`) are rejected by PR #86 and accepted by the prototype.
   - They give a rate of exactly 0 at zero products, at any substrate level including 1e9. At products 1e-6 the rate scales like Pⁿ (about 1e-7 for n = 1, about 1e-27 for bi-bi n = 2).
   - `_kcat_forward` throws "no kcat components" in 24 of 24 runs (`logs/degenerate_proto.log`).
   - The oracle confirms the zero: at P = 0 all enzyme sits in the inactive sink, inactive mass 1.000000 (`logs/debug_sink.log`).
   - The loss does flag them: 104.8 on 8 points, 5 of them product-free. But `rescale_parameter_values` throws (`logs/fitpath.log`), so `fit_rate_equation` would throw after optimizing. From code reading, identify would record a FitFailure after a full fit.
   - The probe's kcat check left these shapes out. "Adopt now" needs a fourth addition: reject these, or define kcat for them.
   - Enumeration does not produce them. All 1,954 gen1 allosteric ping-pong laws have a nonzero rate at zero products and a finite kcat (`logs/pp_scan.log`).
2. **The flip form is chosen by name.** The DSL documents `Estar` and `E_c` as valid enzyme names. Under the prototype these are rejected, and swapping the names E and Eprime changes the law: at equal parameters it differs by 21–100% (the flip-at-Eprime column above). The proposed fallback only covers "no E and one free form".
3. **The probe's oracle was not independent.** It takes constant names and values from the package (`_state_parts`, `_dependent_param_exprs`), so it could not catch a wrong inactive-state constant. Mine found no such error.
4. **Unrelated, pre-existing, same in both checkouts.** `expand_mechanisms` throws on 1 of 25 sampled ter-bi allosteric parents ("more than 5000 polynomial terms"). The throw comes from `_split_gain_test` → `_dependent_param_exprs`, so one oversized split candidate aborts that parent's whole expansion. It is not in ROADMAP.md.

## Not checked

- Regulatory sites combined with a second free form.
- n > 3; n = 3 was tested on one shape only.
- Ter-ter enumeration.
- Generation 2 beyond 25 parents per reaction.
- A real `fit_rate_equation` or `identify` run.
- R3 numerically (it has the same dead-end structure as R4).
- Whether swapping E and Eprime gives a different law family, not just different values at equal parameters.
- The slice's `admit`.
- I did not rerun the test suite; the probe's log does show 32,824 pass, 5 fail, 1 error (E_c), and Aqua and JET passing.

Every Julia run went through the lock with a kill timeout, about 23 minutes in total, and none is left running. Nothing in the main checkout, `slice/src` or the probe's worktrees was modified.