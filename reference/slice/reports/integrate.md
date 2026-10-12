I've integrated and run the whole slice: 257,698 mechanisms, 237,965 derived laws. Every law agrees with its oracle, the compiled code equals the law with 0 allocations, and the key had no false merge. Three results need your call, Denis: 343 large ter-ter laws run over 120 ns, 2.8% of ter-ter laws have no unique kcat, and fitted names change in 9.4% of same-model regroupings. `slice/REPORT.md` is not written because the harness blocks report files; the full tables are in `logs/final_summary.log`.

**What I built** (code lines from `codelines.py`; `linecounts.py` gives the same)
- **`key.jl`, 77 lines** (critic projected 75–87). It brings the law to lowest terms, labels the fitted constants by colour refinement, falls back to searching tied colours, and hashes with SHA-256. Rational exponents work. Concentrations are labelled by their position in the metabolite list, not by name.
- **`pipeline.jl`, 19 lines.** `derive(m)` runs the checks, solve, law, normalizers, both printed forms, key, generated code and kcat. It errors when a mechanism is rejected and returns the kcat error message when there is no kcat.
- **Edits to reconcile the stages:**
  - Moved two graph helpers into `types.jl` and rewrote `bottomless` on them.
  - Shared the kcat and coefficient-printing helpers between plain and allosteric code.
  - Routed allosteric mechanisms correctly through `kcat(m, params)` and one-argument printing.
  - Added an `expand` driver that drops rejected tag sets with `admit`.

  Every stage output is unchanged on all 4,917 core records (16 artifacts per record compared before and after).
- Slice total goes from 1,952 to 2,033 lines. The derivation parts sum to 1,110 lines, against 926 today and the critic's ~970. Printing (174 lines against 72 today) causes most of the overrun: the coefficient form, the Cleland/affine normalizer and the allosteric printing are new. `types.jl` is 474 (746 today), names 88, enum 308.

**Populations**
- Core: 4,917 records.
- New: `init_mechanisms(ter_ter_rxn)` gives 250,855 mechanisms (exported in 30 s). One expansion round of 60 sampled parents gives 1,926 children (535 plain, 1,391 allosteric).
- The package's own expansion crashed on one parent: its `:OnlyA` guard builds an invalid mechanism.

**Results**

| Check | Core | Ter-ter |
|---|---|---|
| Oracle, worst relative error | 7.0e-31 | 5.9e-26 (see below) |
| Laws compiled | 4,786 (all) | 26,679 (largest tenth + 2,000 random + all expansion children) |
| Allocations | 0 | 0 |
| Slowest call | 113 ns | 184 ns |
| Over 120 ns | 0 | 343, all with more than 1,000 denominator terms |

- **Oracle:** 189 ter-ter laws are above 1e-30. On the 12 worst, the error shrinks 10^20-fold when the oracle's rapid-equilibrium scale goes from 1e40 to 1e60, so it is the oracle's error, not the law's.
- **Generated code vs law:** Float64 error at most 3.4e-15 after allowing for cancellation near equilibrium; BigFloat inputs agree to 3e-153; every zero concentration gives a finite, correct value.
- **Key invariance:** 0 changes in all of these:
  - a second process (37,633 keys);
  - 21,908 rewritten inputs (order, directions, group labels);
  - 37,849 same-model regroupings;
  - renaming metabolites without changing their order.

  Renaming that reorders them does change the key, because positions move.
- **Key classes:** identical to P8's redesign classes in all six shared pools. All 997 merges were checked through the oracles (worst 3.9e-80), and deliberately swapped name maps were caught in 751 of 752 cases.

**Today's bugs are absent**
- Mixed groups: AD1 ×3 and AD11 are rejected; AD14, AD14b and AD15 are admitted.
- Two free forms: AD2a–e and AD5 are rejected with a message; no NaN.
- The four "no kcat components" cases now get "never saturates; fit in absolute mode".
- expand_bibi/642 runs in 42 ns.
- P2's 663 same-model splits: 0 key changes (package 327) and 8 fitted-name changes (package 312).
- Admission and kcat reproduce the earlier stages exactly (kcat matches the package on 818/825 plain and 303/305 allosteric mechanisms).
- `expand` reproduces the package's 535 plain ter-ter children exactly.

**Ter-ter scaling** (largest size bin, 2,001–5,000 denominator terms)
- Derivation 0.86 s, compile 1.0 s, 29 KB of native code, 94 ns per call, 277 of 1,023 laws over 120 ns.
- Growth with law size is mildly super-linear: law 1.23, key 1.34, compile 1.24 (log-log slopes).
- Largest law: 3,236 terms, 2.9 s to derive.
- Today's package on the 25 slowest laws takes 2.55 µs per call against 177 ns here. On random laws it is 233 ns against 26 ns.

**Failures, all fixed**
- **Engine size cap:** it counted spanning trees without parallel steps, so 77 ter-ter laws got through at up to 12× today's limit (10,960 terms, 20 s each). I changed the cap to today's count (`complexity`). Core outputs are unchanged and 19,376 ter-ter laws are now refused.
- My kcat refactor called `zero` on a polynomial value, which has no method; a stage-by-stage probe caught it.
- Check-script mistakes: per-mechanism NamedTuples made the oracle checks recompile, a 0/0 zero check gave 2 false failures, there were scope and name clashes, and NaN broke a maximum. Separately, my own `pgrep` pattern killed my shell once; nothing else was touched.

**Deviations**
- The key labels concentrations by position, as noted above.
- Ter-ter key invariances are sampled 1 record in 50.
- Ter-ter compile checks skip the zero-concentration and BigFloat-input runs.

**Questions for you**
1. The 343 laws over 120 ns: lower the size cap, speed up code generation, or accept them for ter-ter?
2. 6,476 ter-ter laws have no unique kcat. A substrate binds a product-release form (D on E(R)), and Q rebinding closes a cycle, so the dead-end filter keeps it. Should the filter look for dead ends with the products at zero?
3. Fitted names change in 3,561 of 37,849 same-model regroupings (12% on ter-ter). Should each group be named after the member the solve keeps?
4. Should printed local names (`N`, `Q_A`, `d_A`, …) be reserved, like `Keq`?
5. The tagging move must call `admit`: today it produces 226 rejected ter-ter tag sets out of 1,391.

No Julia process of mine is left running.

All files are in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/`:
- `src/key.jl`
- `src/pipeline.jl`
- `INTERFACES.md` (new Integration section)
- `integrate/logs/final_summary.log`
- `integrate/logs/check_bugs.log`
- `integrate/logs/compare_regress.log`
- `integrate/logs/compare_regress2.log`
- `integrate/logs/pkg_bench_terter.log`
- `integrate/logs/linecounts.log`