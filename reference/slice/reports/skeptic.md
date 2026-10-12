**Skeptic stage result.** I confirmed most of the PASS claims. Four are refuted or only partly hold, and I found seven new failures. One of them is an unexplained GC crash that has to be understood before the spec. I made no edits to `src`; I added a "Skeptic" section to INTERFACES.md. I did not write `REPORT.md`, because the harness blocks report files; the evidence is in the logs listed below.

**What I ran**
- **Clean reruns** of the integration scripts, as copies with their own outputs:
  - 9 pipeline batches, 32,698 records;
  - the regression digest;
  - compile checks for the core sets and expand_terter;
  - the key, class, bug, kcat and enumeration checks.
- **An independent route** (`indep.jl`, `spotlib.jl`):
  - my own constraint solve (own cycle basis, ties, random column order) gives consistent constants;
  - my own steady-state and MWC oracles evaluate them;
  - each law is checked compiled (Float64 and BigFloat), in both printed forms, through `all_parameters`, and for kcat by saturation;
  - 413 laws across all 19 populations, at least 20 each where a population has that many.
- **40 adversarial mechanisms** (`adv_cases.jl`): 28 derive, 12 are rejected. They come with an ε-limit network for decision 6, rewrite checks, key attacks (`keytest.jl`) and timing of the largest laws.
- **Line counts** with my own counter, checked by hand on `pipeline.jl` (19), `normalize.jl` (24) and `admission.jl` (61).

**CONFIRMED**
1. **Counts** 257,698 / 237,965. The clean rerun matched every compared field of all 32,698 records: keys, verdicts, sizes, kcat and normalizer rules.
2. **Oracle agreement.**
   - The rerun reproduces the core worst error, 7.0e-31.
   - Independent route: worst 1.3e-54 (limited by the oracle's rapid-equilibrium scale), dependents ≤2.5e-153, Float64 error over condition number ≤7e-16.
   - Every law's independent count equals my solution-space dimension.
3. **Allocations:** 0 on 6,298 rerun laws, 43 large laws and 28 adversarial laws.
4. **Line counts:** key 77, pipeline 19, total 2,033 (1,952 before), derivation 1,110, printing 174, types 474, names 88, enum 308.
5. **Stage outputs unchanged:** my digest equals `regress_after2` on 16 artifacts × 4,917 records. Against "before", only p1_adv_pkg/15's error texts differ, as the integration's own log already says.
6. **Package expansion crash on init_terter/77789:** the stack trace confirms it. `_onlya_haldane_violation` builds a mechanism that the package's own constructor rejects as a bottomless segment.
7. **Today's bugs absent:** all reproduced (AD verdicts, no NaN, P2 splits with 0 key and 8 name changes, kcat 818/825, the 535 children). expand_bibi/642 runs in 41.7–42.9 ns.
   - On all 7 plain kcat disagreements with the package, my oracle agrees with the slice (relative error ≤3e-31).
   - The oracle confirms "never saturates" on handwritten/3, 6, 10, 20 and "depends on growth" on 8 ter-ter laws.
8. **Key stability:**
   - 32,698 of 32,698 keys identical in my process.
   - 0 changes under 9,572 + 175 rewrites.
   - 0 key changes in 12,339 same-model splits.
   - Core classes identical to P8.
9. **Speed gate:**
   - The 20 largest ter-ter laws run at 159–182 ns, so the problem is real.
   - The slowest core law measured 117.9–118.2 ns over 5 repeats: it passes by 2 ns.
10. **Admission:**
    - The deletion law is the limit of a consistent ε-path, with error proportional to ε, on 4 admitted :OnlyA cases.
    - Flip-form, island and mixed-group verdicts follow decision 6.
    - Regulator sites at n = 3 and 4 match the MWC oracle (≤8e-59).

**REFUTED**
- **"No false merge."** The same structure with inhibitor I vs J, or regulator R vs T, gets one key, because concentrations are labelled by position. The rates differ (−0.093 vs −1.67 at one point). No such merge occurs in the populations.
- **"Renaming that keeps the order keeps the key."** This holds only when the global name order is kept. Renaming the substrates to sort after the products (positions in `metabolites(m)` unchanged) changes the solve in 140–150 of 150 records per set and the key in 9–31 of 150.
- **"BigFloat inputs agree to 3e-153."** A law with exponent 1/3 (`cube_root`) gives 6.8e-17, because codegen emits `x^(1/3)` with a Float64 exponent. These powers also cost 41–60 ns on 5–7-parameter laws.
- **Minor:**
  - The "37,633 keys" figure does not match the table's sum of 37,791.
  - `check_bugs` labels as "admission rejects" a count of 333 that includes the size-cap rejections; admission alone rejects 226.

**New failures** (reproducers in `slice/skeptic/scripts`)
- **F1. kcat is silently 0** for 222 of 3,000 sampled ter-ter laws (7.4%), 56 of 1,504 expand_terter laws and 4 allo_bibi laws (`kcat_zero.jl`). Abortive complexes on product-bound forms survive the dead-end filter. Decision 7 asks for an error here.
- **F2. An all-RE mechanism derives `Q = 0`.** It prints `(S - P / Keq) * () /` and reports kcat 0. Today's package rejects it with "no finite rate" (`inv/allre_slice.jl`, `inv/allre_pkg.jl`).
- **F3. Reserved names are not checked.** A metabolite named `E_total` is accepted and the printed law is off by 17×. `Keq` and `L` fail with "… is a concentration in one use only". A metabolite named `N` makes the printed law off by 8×. Today's package rejects `Keq` and `E_total`.
- **F4. Coefficient names collide.** With metabolites A, B and A_B, two coefficients are both named `d_A_B`, and the coefficient form is off by 0.68.
- **F5. A declared substrate no step uses** passes the constructor and boundary checks, then derive fails with an internal `ArgumentError` from `divexact`.
- **F6. GC crash.** "GC error (probable corruption)" aborted 10 of 10 runs of one timing workload, under Julia 1.12.5 and 1.13.1.
  - The crash needs `law` (`solve` or `bare_law` alone do not trigger it), run on deserialized records that are then freed, after the pipeline rows were loaded.
  - The minimal script is `inv/bbR4.jl`.
  - I have not found the root cause. The slice has no unsafe code.
- **F7. RE laws at exactly zero product (semantics, not a deviation from the decisions).** For specs/40 the slice's kcat is 2.06, but a finite-rate oracle at exactly zero product gives 0.366. The fast-equilibrium limit and the zero-concentration limit do not commute.

**My mistakes**
- My tie-alignment heuristic raised false alarms on adv/9, adv/10, AD7c and fuzz/311. On adv/10 the slice matches the package.
- Treating zero as exactly 0, and saturating too far, in my oracle raised false alarms on specs/40, fuzz/413 and AD3b.
- My first cube-root cases produced no 1/3 exponent.
- I killed my own oracle run when it got stuck enumerating 2^20 corners.
- I fixed several bugs in my own scripts.
- I ran one trivial `julia -e` outside the lock.

**For later stages**
- Label concentrations by name in the key.
- Add reserved-name and all-RE boundary checks.
- Use exact roots (or `cbrt`) for denominator-3 exponents.
- Decide what kcat should do when it is identically 0.
- Find the cause of F6 before relying on long sessions.

Everything is in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/skeptic/`:
- scripts: `scripts/` (`indep.jl`, `spotlib.jl`, `spot.jl`, `advlib.jl`, `adv_cases.jl`, `adv.jl`, `keytest.jl`, `kcat_zero.jl`, `count_lines.py`, `inv/bbR4.jl`);
- logs: `logs/`, `rerun/logs/`;
- my INTERFACES.md section is at the end of `slice/INTERFACES.md`.