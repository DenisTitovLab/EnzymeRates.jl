I changed no files in `slice/src`. I did not write `REPORT.md` because the harness blocks report files, so this message is the report. I added a "Critic" section to `INTERFACES.md` that lists the entries the code no longer matches. My scripts and logs (`timing_check.jl`, `kept_names.jl`, `density.py`, `logs/`) are in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/critic/`.

**What I checked myself**
- **The timing protocol measures real work.** Fixed inputs and 64 changing inputs agree within 2%:

  | Law | Fixed | Changing |
  |---|---|---|
  | expand_bibi/642 | 39.5 ns | 39.6 ns |
  | bi_bi_with_inh/1013 | 112.2 ns | 111.2 ns |
  | init_terter/62921 | 163.1 ns | 166.5 ns |

- **Large laws are limited by how fast the CPU does floating-point (FP) arithmetic.** LLVM keeps 629 of 642's 1,436 counted operations as FP instructions, and 2,593 of 7,670 for 62921. All three laws run at 15–16 FP instructions per ns, so 120 ns buys about 1,900; 62921 needs 2,593.
- **Line counts reproduce.** The slice's code is not denser than the package's (31 against 35 non-blank characters per code line), so its counts are comparable.
- **Naming each fitted constant after its own edge** (the member the solve keeps) changes 554 of 27,649 fitted names (2.0%). Names then no longer depend on grouping.

**What the spec can rely on**
- **Values:** all 237,965 laws agree with the oracles. For the 12 worst of 189 ter-ter outliers, the error is the oracle's. The relation space equals the package's on 857 of 857 plain laws, and no compiled law allocates.
- **Determinism:** keys are stable across processes and about 31,700 input rewrites. The solve is identical across 37,849 same-model regroupings.
- **Admission is exact:** it matches an exact LP on 2,995 records and the skeptic's ε-limit checks, and rejects the 81 bad ping-pong tag sets.
- **Speed and compile cost:** expand_bibi/642 runs in 40 ns against the package's 121 ns. Ter-ter laws run 8–15× faster than today. A fresh mechanism costs at most 126 method instances and 76 KB.
- **Enumeration:** the three ported moves give the package's children; every difference comes from known issues 1 or 2.
- **Size cap:** it uses today's measure, V×τ ≤ 5,000.

**What the spec cannot rely on (fix in the spec)**
1. **Key:** it labels concentrations by position, so a model with inhibitor I and one with inhibitor J share a key. Label them by name. The key also changes when metabolites are renamed.
2. **Missing boundary checks:**
   - reserved names (E_total, Keq, L, N, the printed MWC locals);
   - mechanisms where every step is rapid-equilibrium (they give Q = 0);
   - a declared substrate that no step uses;
   - coefficient names that collide (A, B, A_B).
3. **kcat:** it is silently 0 for 7.4% of sampled ter-ter laws and for AD3b/d, because substrate inhibition inside a cycle gets past the dead-end filter. Another 2.8% of ter-ter laws raise "depends on how A, B, D grow".
4. **Exponents p/q with q ≥ 3** compile to Float64 powers. Use `cbrt` or a rational power.
5. **Two choices nobody decided:** parameters come out in content order (the package sorts by name), and ties are broken by packed content rather than decision 4's content name.
6. **Allosteric size cap:** it applies to each conformation separately, not to their sum.
7. **F6** (the GC corruption) is still unexplained.

**Contradictions between stages**
- Integration labels concentrations by position; P8 labelled them by name, and the skeptic found false merges.
- Integration says renames that keep the metabolite order keep the key; the skeptic saw 9–31 of 150 keys change.
- The ROADMAP says the rewrite removes name churn, yet 9.4% of same-model splits still change fitted names.
- Stale INTERFACES entries: "against = slot 2", the cap inside `spanning_trees`, and the `isassigned` slowdown.
- Counts: 37,633 vs 37,791 keys. The 333 "admission rejects" include 188 size-cap rejects; admission alone rejects 226.
- "types.jl 746 → 474" ignores that EnzymeReaction (143 lines) stays.

**Code lines, today (HEAD 1907f9e) against the slice**

| Component | Today | Slice |
|---|---|---|
| Engine, polynomials, cap | 193 | 302 |
| Joint solve, API | 202 | 205 |
| Allosteric law | 96 | 48 |
| :OnlyA admission, kcat | 155 | 166 |
| Codegen | 172 | 119 |
| Printing, normalization | 72 | 174 |
| Key | 6 | 77 |
| API, pipeline | 30 | 19 |
| **Derivation** | **926** | **1,110** |
| Type layer and names | ≈708* | 562 |
| Three enumeration moves | ≈411* | 308 |
| **Total** | **≈2,045** | **1,980** |

\*These two come from my own region counts, so they are estimates.

The derivation grows 20%, and printing alone adds 102 lines. Overall the redesign comes out about the same size (−3%), before roughly 30–40 more lines for the fixes above.

**Risks and probes**
- **No new vertical probe is needed.**
- **F6 needs one probe of at most half a day before PR (2).** It matters because identify's workers deserialize mechanisms. The probe should:
  - rebuild the records in the process instead of deserializing them;
  - run on Julia 1.10 and on x86;
  - report the bug upstream if Julia is at fault.
- **Time the speed gate on CI hardware within PR 0.** bi_bi_with_inh/1013 (112–118 ns here) has no margin.

**PR sequence (adjusted)**
- **Fix now:** `_onlya_haldane_violation` builds a mechanism that its own constructor rejects, so `expand_mechanisms` crashes on init_terter/77789.
- **PR 0:** Horner code generation with locals for the Haldane factor, the RE segments and the MWC terms (decision 9).
- **PR 0b:** tests only — rewrite test assertions to check counts and values rather than literal names and strings.
- **PR (1):** joint solve, single-graph allosteric and `all_parameters` on today's types. Its tie-break must not depend on which type represents the mechanism (for example, the step's DSL text); otherwise PR (2) changes the dependents again.
- **PR (2):** graph Mechanism, DSL, show/CSV, names, boundary checks, Sig.
- **PR (4) moves before (3):** :OnlyA admission is 61 lines, does not need the engine, and fixes verdicts that are wrong today. The flip-form rule is already in commit 1907f9e.
- **PR (3):** engine, normalization, printing, key, kcat.

**Test pins that will move**
- Dependents differ from the package's in 564 of 806 plain mechanisms.
- Allosteric fitted counts drop in 132 of 1,477 (allo_bibi), 29 of 400 (allo_bibi_pp) and 28 of 192 (allo_terbi) mechanisms.
- Test lines affected:
  - 55 lines with `_res_` and 98 with `inh`;
  - at least 15 string assertions, 2 golden files and 11 factored-form pairs;
  - 30 lines on Full/Reduced and 61 on eq_hash;
  - 58 lines on kcat (7 of 825 plain values change, AD3b/d change, PFK now errors);
  - 20 lines on the guards (AD1, AD11 and AD15 flip verdict).
- These may move: the Haldane, Wegscheider and mirror counts in MECHANISM_TEST_SPECS.
- These will be rewritten: the Expr-shape tests and the compile-budget Sig numbers.

**Questions for Denis**
1. **Speed gate:** a hard 120 ns gate for MECHANISM_TEST_SPECS and expand_bibi/626, 642 and 658, plus a bound that scales with law size for the 343 ter-ter laws over 120 ns? I recommend this.
2. **kcat at 0, or kcat that depends on growth:** raise an error and fall back to absolute mode (I recommend this), or leave out abortive bindings?
3. **Fitted names:** name each class after its kept member? About 2% of names change, and the renaming on regrouping stops. I recommend this.
4. **Parameter order:** keep alphabetical? I recommend this.
5. **Metabolite names:** reject `_` in them, and reserve N and the printed locals? I recommend this.
6. **F6:** do you approve the probe?