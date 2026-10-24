# Claim audit: law.md, unit.md, phantom.md, keycost.md

Paths are relative to `scratchpad/phase45/research5/` unless they start with `slice/`, `phase3/`
or `phase45/`. "Recount" means I recomputed the number from the run's own data or logs.

## 0. Answer

- 0 blocking, 3 major, 15 minor findings.
- Most numbers reproduce. Recounted and confirmed:
  - hoistscan counts for every variant and version (law: 710/670/612/841; 861/921/850/935; U's 2 hoists in both);
  - the timing ratios, R5.5 percentages, the unit census tallies, the 42-spec printed lengths, and the phantom totals.
- The major findings:
  - M1: the law type that passed the gates is not the one L1 proposes.
  - M2: phantom's ping-pong round 1 duplicates bi-bi round 1, so the totals double-count.
  - M3: no report maps phase 4's active mask (catalytic `:OnlyI`).

## 1. Blocking

None. No ruling rests on a wrong number that would change it.

## 2. Major

**M1. L1's law type was never built. The variant that passed the scan and timing keeps the parametric types and the Union field.**
- `law/variants/C` differs from the slice only by the deleted retry (`diff`: `C/engine.jl` vs `slice/src/engine.jl:175-180`; `C/allosteric.jl` vs `slice/src/allosteric.jl:40-42`).
- So C keeps `Poly{T}`, `Law{T}`, `AlloLaw{T}` and `inactive::Union{Nothing, Law{T}}` (`slice/src/allosteric.jl:15-24`).
- C's `rate_law` returns `Union{Law{Int}, AlloLaw{Int}}`, a Union of parametric structs. R18 forbids that by its letter.
- `scan.jl:49` reads the `Union{Nothing, Law}` field inside the loop.
- L1 (`law.md:128-131, 186-187`) proposes non-parametric `Law`/`AlloLaw` with the empty inactive conformation stored as (0, 1). It cites section 1.2's gates, which ran on the other layout.
- The evidence still points the same way: the stricter layout removes Unions rather than adding them.
- Proposed ruling:
  - Ruling: build the proposed layout (C′) and rerun `run_scans2.sh` on it (about 2 min per version) before the plan; otherwise L1 says "tested as Int-only parametric C".
  - Why: the claim and its evidence must describe the same code.
  - Cost if wrong: an untested layout reaches the PR, and only L7's package scan catches it.

**M2. Phantom's ping-pong round 1 duplicates bi-bi round 1, and every total double-counts it.**
- Recount from `phantom/data/facts_{bibi,pingpong}.tsv` and `allokeys_*.tsv`: ping-pong round 1 flags 188 variants, all `to_allosteric`, on 76 keys. These are the same 76 keys bi-bi round 1 flags (overlap 76 of 76).
- Both populations admit 1,833 variants. The 331 ping-pong-specific variants are all inadmissible (`phantom/logs/summary.log:3,5`).
- The certificate counts also repeat: 716 phantoms and 528 certificate-only phantoms in both rounds (`phantom/logs/extra.log:2,4`), and 512 removed fits in both (`phantom/logs/fits.log:2,5`).
- Distinct counts:
  - flagged: 738, not 926;
  - `to_allosteric`: 472, not 660;
  - certificate phantoms: 3,102, not 3,818;
  - "76 of 2,535" fits removed in ping-pong (`phantom.md:14`) is bi-bi's 76 again.
- `phantom.md:14, 57` presents ping-pong as an independent third population.
- No ping-pong-specific one-subunit variant was tested. Under phase 4's tagging (P4-1) some may become admissible. `phantom.md` §5 notes that P4-1 is not covered, but not this.
- Proposed ruling:
  - Ruling: report distinct counts, and state that ping-pong-specific phantom coverage is zero.
  - Ruling: rerun R5.2 on ping-pong once the package's P4-1 tagging exists (phase-5 Task 13 gate).
  - Why: the R5.2 pass claim covers bi-bi only.
  - Cost if wrong: a ping-pong phantom class that bi-bi lacks goes unseen.

**M3. No report maps phase 4's active mask.**
- Phase 4 makes catalytic `:OnlyI` legal on modifier-binding groups and adds `_active_mask(am)` (`phase45/research/phase4_research.md:234-253`; P4-16, P4-17).
- The phase-5 plan, step 5, requires "both masks" (`phase45/understand/SUMMARY.md`). Engine critic E3 says the same (`phase45/understand/engine.md`).
- Where the reports miss it:
  - `keycost.md:246` sets the cap and filter to `complexity(catalytic(am), trues) + …`, which counts `:OnlyI` edges the active conformation lacks.
  - `keycost.md:252`'s codegen mapping names only the inactive mask.
  - `law.md:121-131`'s `law(am)` mapping names only `_inactive_mask`.
  - `phantom.md:158` takes the L→0 law to be `law(catalytic(am))`. That is false when a group is `:OnlyI`. It is harmless in practice: the pre-filter excludes `:OnlyI` variants (they bind in I), and the certificate fails on them (X is not a sub-sum of Y). The report should say so.
- Proposed ruling:
  - Ruling: replace `trues` with `_active_mask(am)` in the law, the cap, the filter, codegen and the key's MWC parts. Phase 4's `:OnlyI` fixtures join Task 5's gate.
  - Why: the phase-4 rulings P4-16 and P4-17.
  - Cost if wrong: an `:OnlyI` mechanism gets a wrong active law and a wrong cap count.

## 3. Minor

1. **The overflow proof covers less than "Int64 cannot overflow".** `law.md` §1.4 bounds one conformation's law only.
   - Integer polynomial products also occur in codegen (`fI*DA`, `fA*DI`; `slice/src/codegen.jl:81-85`) and in the phantom certificate (`H*fI*QA*WA`; `phantom/scripts/plib.jl`, `certificate`).
   - These stay below about τ_I·F·τ_A ≈ 4e8, so they are safe, but no proof states it.
   - The bound also needs the cap check to stay inside the law entry (`slice/src/engine.jl:147-150`). L1 and L2 do not say so.
2. **L3 leaves `complexity`'s return type open after the BigInt fallback.** A tree count above 2^63 does not fit an `Int`.
   - `keycost.md:241` ports the BigInt determinant "line for line", which contradicts L3.
3. **Q11's byte identity covers integer exponents only.** The run had 0 non-integer exponents (`law/logs/q11.log:3`).
   - `build_power_expr`'s Float-exponent path and L4's `Rational{Int}` narrowing of fractional exponents are untested.
4. **The R5.1 speed gate was replaced.** The gate was ≤ 1.2 × 70.7 µs; the report judges a same-run ratio to U instead (`law.md` §1.3, which says so).
   - U itself runs 55-120 µs on these sets (`law/logs/timing_*.log`, "sl all" lines), so the absolute gate was not evaluated. Denis should see that the gate changed.
5. **The identify workload is thin.**
   - `NFIT = 10` is reached in pass 1, so no child of the beam is ever fitted (`law/scripts/scan.jl:108`; log tuple `(10, 56, 10, …)`).
   - The sample's children are the smallest-id ones, not random (`law/scripts/extract_sample.py`, PLAN). The report does not say so.
6. **`phantom.md:10` overstates its containment numbers.**
   - It says every flagged variant lies inside the L→0 family "to ≤2.7e-33". Numeric containment ran on 43 of 926.
   - The 40 certified variants sampled reach 1.5e-31 (`phantom/logs/contain_summary.log:2,4`).
   - Harmless against the 1e-12 pass mark.
7. **"Strictly larger family" for 3 rank-equal variants rests on a local optimizer failing to fit** (residuals 1.8e-1, 1.4e-1 and 9.8e-3; `contain_summary.log:14,15,18`).
   - Optimizer non-convergence does occur in this run: the L-free fit is worse than the fixed-L fit, for example 1.5e-16 against 8.0e-6 (`contain_summary.log:52,57`).
   - Label the result as numeric.
8. **Phantom ruling 1's case for the certificate is weaker than stated.**
   - 25 µs against 1.3 ms (`phantom/logs/timing.log:2-3`): both are under 1% of one fit, which takes 0.21-1.33 s (`keycost/logs/key_vs_fit.log:1-5`).
   - A seeded modular rank is also deterministic.
   - Inside the pre-filter, rank equality and family equality never disagreed on the populations; the 3 counterexamples lie outside it.
   - The real trade-off: a proof against 3 of 926 flags lost.
9. **`keycost.md:18-19` compares with the wrong stratum.** "Cost what dropped plain laws cost" cites admitted plain laws (0.091 ns per V×τ).
   - Dropped plain laws cost 0.073 (`keycost/logs/r55_cost_ratio.log:4`).
   - Per candidate the claim holds: 33.5 against 34.8 ns (`keycost/logs/r55_cost.log:4-5`).
   - Per V×τ unit the allosteric laws cost 26% more, so the sum understates their cost. The conclusion stands.
10. **F9's "package helper" is an emulation.** The timing was `klib.jl:193-212`'s copy of the helper's loop on slice-generated bodies, not the helper on package code.
    - The 4 hoisted specs (`keycost/logs/c_calib.log:11`) predict phase-5 behaviour; they say nothing measured about today's package.
11. **F3 still lacks a time bound.** At 1,000 leaves, a large law can still take about 1,000 × 0.35 s ≈ 6 min.
    - Slowest single-leaf key: 351 ms (`keycost/logs/kt_terter_summary.log:11`).
    - Propose a time budget alongside the leaf cap.
12. **R5.5's ping-pong populations use phase-2 tagging, not P4-1's.** P4-1's inactive conformation runs a half-reaction, which raises V×τ(inactive).
    - The conclusion is likely robust: the largest admitted sum is 50-72 against 337 (`keycost/logs/r55_table.log:20,26`). The report does not state the caveat.
13. **U2 cites SPEC:38-39 ("Cleland-style coefficient form") in support of Cleland's N.** That line is about the coefficient form, which the graded basis also produces.
    - The real evidence is ter-ter: 257 of 292 eligible conformations differ (`unit/logs/summary.log:133`).
14. **`unit.md:18` mislabels a cause, and U3 leaves a direction unstated.**
    - "A square does in 27" mislabels 23 ter-ter cases that `unit/logs/why.log:5` classes as "lower-degree substrate term + a power".
    - U3 does not state its exponent-vector direction: larger exponent of the earlier metabolite first (`unit/ulib.jl:167-180`, `Tuple(-ex)`). That direction puts A² before AB.
15. **Cross-report inconsistencies and a shortfall.**
    - R5.4 asked for 2,000 ter-ter laws; 1,865 were derived because 135 hit the cap (`unit/logs/summary.log:125-128`).
    - `keycost.md:216` uses the old 2,630-character median ter-ter law; unit's layout gives 4,534 (`unit/logs/summary.log:168`).
    - `law.md:191`'s "a non-escaping read would hoist as U's does" for B is an inference that was not run.

## 4. Checked and confirmed

- Positive control: 2 dangerous hoists and a GC abort (exit 134) on 1.12 and 1.13 (`law/logs/remarks_repro_11{2,3}.log`).
- Variants C, B and A: 0 hoists. Rounds identical, 0 GC errors, both versions.
- The 995-s and 744-s run times include waits for the lock, since `phase3/jl.sh:6-10` starts the clock before `flock`. Every Julia run used a timeout of at most 600 s.
- R5.6: 1,390 ter-ter laws, 1 leaf each; key time 1.21% of fit time at the 5,000 cap and 0.34% at 337 (`keycost/logs/key_vs_fit.log:6-7`).
- R5.7: 0 allocations in 578 laws (`keycost/logs/c_calib.log:1`). c = 2 × 0.108 × 4.58, rounded up to 1.0. The CI readings match `graph_mechanism_wt/test/test_rate_eq_derivation.jl:1309-1310`.
- Unit:
  - 9,406 laws, 0 unit failures;
  - N differs from Cleland's in 264 of 382 eligible conformations;
  - 29 of 42 specs print shorter, 82,422 → 70,011 characters (`unit/logs/examples.log:48-89`);
  - U7's names cannot clash, since R10 reserves names containing `_` (`phase3/RULINGS.md:62-65`).
- Phantom:
  - 18,915 variants;
  - modular rank equals BigFloat rank on 37,830 ranks;
  - the site clause changes no verdict (recount: the 3 rank-equal regulator variants bind in the catalytic graph);
  - the certificate math holds, with `Keq` and `L` weighted 0 (`phantom/scripts/plib.jl`).
- No report touches the version (R19).
