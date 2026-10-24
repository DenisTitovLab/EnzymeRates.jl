# Law type and engine output (R5.1; engine.md Q1, Q9, Q10, Q11, Q14)

Paths are relative to `scratchpad/phase45/research5/`. `law/` holds this study's scripts,
variants and logs (`law/logs/`); `slice/src/` is `scratchpad/slice/src/`; WT3 is
`scratchpad/phase3/proto3_wt`; GM is `scratchpad/graph_mechanism_wt` (7205ba4), which built
the inputs (`law/corpus/`, 81,498 mechanisms, as base6; `law/logs/corpus.log:8`).

## 0. Answer

- Law type: **C, one concrete law type over checked `Int64` with no `BigInt` retry**; an
  overflow is an internal error. Candidates A (widen to `BigInt`) and B (`Dict` alias) also
  pass R5.1's gates; C passes them with less code and no extra live bytes.
- GC gate, 1,535 census mechanisms (642 allosteric), derivation and identify workloads,
  `--pkgimages=no`, LICM remarks, `hoistscan.py`:
  - Julia 1.12.5: C 0 hoists of 710 escape analyses, B 0/670, A 0/612; the slice (U) 2
    dangerous of 841 (`law/logs/remarks_scanall_CBA_112.log`, `remarks_scan_U_112.log:10226,10262`).
  - Julia 1.13.1: C 0/861, B 0/921, A 0/850; U 2 dangerous of 935
    (`law/logs/remarks_scanall_CBAU_113.log:56242,56279`).
  - Every variant, both versions: 3 derivation and 3 identify rounds identical, 0 GC errors.
  - Positive control `repro_upstream.jl`: 2 dangerous hoists and a GC abort (exit 134) on
    both versions (`remarks_repro_112.log:36,72,137`; `remarks_repro_113.log:37,74,75`).
- Speed, 2,030 laws, variants interleaved per mechanism: law time over U's, median (p90):
  C 0.991-1.003 (1.019-1.033), B 0.997-1.005, A 1.030-1.038 (1.059-1.076)
  (`law/logs/timing_*.log`, "law all" lines). Same law, printed string, key and code as U
  in all 2,030 for A, B and C ("differences against U: none").
- `Int64` cannot overflow under the 5,000 cap (proof, section 1.4). Measured: largest
  coefficient 36 in 2,642 laws; `BigInt` needed 0 times (`law/logs/checks.log:13-14`).
- A triples live bytes after a full collection: law 2.22 vs 0.71 KiB per mechanism,
  law+print+key 12.83 vs 4.19 KiB (bi-bi; `timing_bibi.log:63-70`).
- Q10: checked-`Int` Bareiss determinant, `BigInt` only after an overflow: equal on
  2,647/2,647 pairs; 0.00 vs 0.39 KiB live bytes per call; "together" = the sum.
- Q11: a structured solve output rebuilds `_dependent_param_exprs` byte for byte on
  1,974/1,974 mechanisms (811 allosteric, 42 specs).
- Q9: the guard's class 2 owns RE turnover; the engine keeps a cheap assertion, since
  the slice engine returns a finite, wrong law on both RE-turnover hand cases.
- Q14: class 5 on the new law costs 0.4 µs median and rejects exactly the two face cases.

## 1. R5.1: law type

### 1.1 Variants (`law/variants/mk.py` patches the slice's files; each patch matches once)

| Variant | Law type | Change |
|---|---|---|
| U | `Union{Law{Int}, Law{BigInt}}` from `law`; `AlloLaw{T}` with `active::Law{T}`, `inactive::Union{Nothing, Law{T}}` | none (`slice/src/engine.jl:173-181`, `allosteric.jl:15-24, 40-42`) |
| A | `Law{BigInt}` only: the `Int` result widened at the boundary | `law/variants/A/engine.jl` (`widen_law`) |
| B | `Poly{T} = Dict{Mono, T}`; `law` still returns the Union | `law/variants/B/poly.jl` (from `research5_interrupted`), `.terms` reads in 3 files |
| C | `Law{Int}` only: `law` = `law(Int, …)`, no retry | `law/variants/C/engine.jl`, retry deleted in `C/allosteric.jl` |

### 1.2 GC scan (`law/scripts/scan.jl`, driver `scan_all.jl`, splitter `split_scan.py`)

- Derivation workload: mechanisms deserialized each round; the loop reads `L.D`,
  `L.active.D`, `L.sigma`, `L.sites`, `L.inactive.D` as the probe's loop did, then prints,
  coefficient form, key, codegen body, class 5 (`scan.jl:36-61`).
- Identify workload: 60 seeds (bi-bi init, ping-pong init, allosteric ping-pong children);
  complexity filter 337, law, key, printed law, fitted names; 10 fits through the compiled
  `rate_equation`; beam of 3 keys expanded by flip and split; CV list deduplicated by key.
- Results (escape analyses / hoisted / hoisted with a bits store into a reference field):

  | | U | A | B | C |
  |---|---|---|---|---|
  | 1.12.5 | 841 / 2 / 2 | 612 / 0 / 0 | 670 / 0 / 0 | 710 / 0 / 0 |
  | 1.13.1 | 935 / 2 / 2 | 850 / 0 / 0 | 921 / 0 / 0 | 861 / 0 / 0 |

- U's two hoists lie outside the library functions: staged calls of each one on concrete
  inputs hoist nothing (`law/logs/remarks_locate_U_112.log`), and the hoist remarks sit
  among the derivation loop's (`scan.jl:52, 55`), the caller-loop pattern of the probe.
  Each hoisted box has one 8-byte reference field written by an integer store, the
  `Poly{T}` layout (`slice/src/poly.jl:63-65`; probe `phase3/gcprobe/report.txt:93-96`).
- U ran 3+3 clean rounds without crashing; the crash needs GC timing, so the scan is the gate.
- B stays clean only because its loop's boxes escape: U and B allocate a box of the
  Union law at `scan.jl:52,55` (escaped: 1; `remarks_scanall_CBA_112.log:22000-22112`,
  `remarks_scanall_CBAU_113.log:27403-27579`); A and C allocate none there.
- The run with the original try/catch loop also showed U's 2 hoists
  (`remarks_scan_U_112_trydriver.log:8707,8743`).

### 1.3 Speed and live bytes (`law/scripts/timing.jl`; best of 5, U run in rotation)

| Set (laws) | U law median µs | A | B | C | U solve+law µs |
|---|---|---|---|---|---|
| specs, two-site, uni-uni (135) | 50.3 | 1.031 | 1.002 | 0.991 | 55.2 |
| bi-bi init + L1 (600) | 78.8 | 1.038 | 0.997 | 1.003 | 100.8 |
| ping-pong init, L1, L2 (800) | 97.8 | 1.030 | 1.005 | 1.000 | 110.5 |
| ter-bi, prefix (495) | 97.1 | 1.031 | 0.997 | 0.995 | 119.7 |

- Ratios are medians of the per-law ratio to U (`timing_spec_twosites_uniuni.log:121-124`,
  `timing_bibi.log:51-54`, `timing_pingpong.log:52-55`, `timing_terbi_prefix.log:87-90`).
- The 70.7 µs reference was solve+law on 806 other plain mechanisms (lost with the outage).
  On the 42 specs C takes 64.1 µs (`timing_spec_twosites_uniuni.log:8`). I judge the 1.2×
  gate by the same-run ratio: C's worst median ratio is 1.005 (`timing_bibi.log:18`).
- Downstream (print + key + codegen body), median ratio: A 1.062-1.087, B 0.996-1.021,
  C 1.000-1.005 ("down all" lines).
- Live bytes after two full collections, KiB per mechanism, law / law+print+key:
  U 0.37-1.14 / 1.50-5.85; C and B equal U; A 0.98-3.15 / 6.01-17.28
  (`timing_*.log` "live bytes" lines). U's law-only bytes come from `complexity`'s
  `BigInt` determinant (0.39 KiB per call, `checks.log:12`; see Q10).

### 1.4 Why `Int64` cannot overflow under the cap

With τ the segment multigraph's spanning trees (parallel steady-state steps counted, as
`complexity` counts them), F the forms the masked steps touch, w_k the copies of the
first substrate steady-state step k moves (`slice/src/engine.jl:141`):

1. Each `W[(a, b)]` sums one coefficient-1 flux per parallel step, so a tree's product has
   |·|₁ = Π multiplicities and `trees[g]` has |·|₁ = τ for every root g (matrix-tree theorem).
2. `den = Σ σ_g·trees[g]` has positive coefficients and |den|₁ = F·τ.
3. `num = Σ w_k(trees[a]·f_k − trees[b]·r_k)`: every partial sum is at most 2τ·Σ|w_k|.
4. Substituting the solve maps monomials to monomials and the shift multiplies by one, so
   no coefficient exceeds the |·|₁ before them.
5. Dividing by H (coefficients ±1) is the chain recurrence b_k = ±(a_k + b_{k−1}), so
   |b_k| ≤ |num|₁.

So every coefficient is at most max(F·τ, 2τ·Σ|w_k|) with τ ≤ 5,000. Overflow would need
about 10¹³ steps. Measured on 2,642 laws (both conformations): |bare den|₁ = F·τ in all; every
coefficient within the bound, at most half of it; largest coefficient 36; at most 15 touched
forms (`checks.log:13-14`). The one exponentiation, codegen raising a single-term state
factor to the power n, saw only coefficient 1 in 76 cases (`checks.log:15`). The slice's
237,965 laws never took the `BigInt` branch (spec line 931; `phase45/understand/engine.md:50-51`).

### 1.5 On phase 3's types (S3_design.md sections 1, 7)

- `law(m::Mechanism, sol, mask::BitVector, inactive::Bool)::Law` reads `m.forms`
  (`bound`/`added`/`subtracted`), `m.edges` (`from`, `to`, `consumed`, `released`,
  `rapid_equilibrium`), `m.group`, `m.table`, where the slice reads `a`, `b`, `cons`, `rel`, `re`
  and `mets` (S3_design.md:101-139).
- `law(am::AllostericMechanism)::AlloLaw` reads `catalytic(am)`, `cat_allo_states`,
  `catalytic_multiplicity`, `regulatory_sites`; the inactive mask is `_inactive_mask(am)` and
  the forest `_segments(m, mask, re)` (S3_design.md:977, 982). The flip form is `forms[1]`.
- Under C, `Poly` holds `Dict{Mono, Int}` and `Law`, `AlloLaw` take no type parameter. A
  population loop over `Union{Mechanism, AllostericMechanism}` gets `Union{Law, AlloLaw}`,
  whose shared fields (`vars`, `H`) have one type each, so no field read boxes.
- A flip-form-only inactive conformation is stored as Q = 0, D = 1 (as
  `slice/src/allosteric.jl:75-76` computes), not as a `Union{Nothing, Law}` field.

## 2. Q10: cap arithmetic (`law/scripts/checks.jl`)

- Checked-`Int` Bareiss equals the `BigInt` determinant on 2,647 of 2,647
  (mechanism, conformation) pairs (`checks.log:2`).
- Median 0.4 vs 0.5 µs per call; p90 0.5 vs 1.9; max 1.0 vs 16.8 (`checks.log:3`).
- Live bytes after a full collection: 0.00 vs 0.39 KiB per call over 20,000 calls
  (`checks.log:11-12`).
- Stress, 390 random multigraph Laplacians (2-14 nodes, multiplicity 1-3): 351 exact,
  0 wrong, 39 overflow; the smallest overflowing tree count is 1.17e10, far above the cap
  (`law/logs/q11.log:4`). On overflow the `BigInt` determinant takes over.
- "Both conformations together" read as V×τ(active) + V×τ(inactive): 0 of 749 sampled
  allosteric mechanisms move across 337 or 5,000 (`checks.log:4-5`). R5.5 counts the
  filter at scale.

## 3. Q11: structured solve output (`law/scripts/q11.jl`, run on GM)

- Per dependent: name, `Keq` exponent and factor exponents as `Rational{Int}`, from the same
  reduction as `_dependent_param_exprs` (GM `src/thermo.jl:367-393`; WT3 `src/thermo.jl:341-369`).
- Rendered with `build_power_expr`, it equals `_dependent_param_exprs` (`==`) on 1,974 of
  1,974 mechanisms: 42 specs, bi-bi and ping-pong init with the children of 25 parents
  each, every 20th ter-bi init (811 allosteric; `q11.log:2-3`).
- Values: worst relative difference 6.2e-16 against `_evaluate_dependent`; 2,028
  dependents, 0 with a fractional exponent; 0 reductions needed `Rational{BigInt}` (`q11.log:3`).
- Consumers to switch: the engine (substitute per column, as `slice/src/engine.jl:185-192`),
  `all_parameters` (drop the `Expr` interpreter `_evaluate_dependent`, WT3
  `rate_eq_derivation.jl:68-90`), codegen's `_dep_assignments` (WT3
  `rate_eq_derivation.jl:379`), the printer's reference lines.

## 4. Q9: where the RE-turnover check lives

- 0 of 2,035 census mechanisms and 2 of the 62 hand cases that build (`re_turnover_beside_ss`,
  `re_turnover_flip2_453`) close an RE cycle with net uptake (`checks.log:20-21, 79`).
- The slice engine, which has no check, returns a finite law for both (Q has 1 term;
  `checks.log:20-21`); the true rate is unbounded (spec lines 429-430).
- The check costs 4.6 µs median standalone, against 81.9 µs for the law (`checks.log:80`).
- The derivation entry runs the guard: `_derivation_input` → `_well_posed`
  (WT3 `src/admission.jl:82-88`), and phase 4 moves the check there as class 2
  (`phase45/research/phase4_research.md:122-136, 303-305`).

## 5. Q14: class 5 on the new law (`law/scripts/lib.jl` `class5`)

- Rule (phase4_research.md:222-230): plain, D in lowest terms keeps a term free of the
  face's names; MWC, fI·DA or fA·DI keeps one (all coefficients positive, sites hold 1).
- 2,092 laws: fails only `ss_dead_end_at_face` and `bottomless_entered`, both at products = 0
  (`checks.log:22, 25, 81`), as phase 4 expects (phase4_research.md:676).
- All 42 specs pass, the two trap specs' hand forms included (`checks.log:30-31`); class 3
  rejects those.
- 0.4 µs median, 1.4 µs max (`checks.log:82`), against phase 4's 1.4-2.3 µs on POLY
  (phase4_research.md:228-231).

## 6. Proposed rulings (ruling — why — cost if wrong)

- **L1 (Q1, R5.1).** Law type C: `Poly` over `Dict{Mono, Int}` with checked arithmetic,
  non-parametric `Law` and `AlloLaw`, the empty inactive conformation stored as (0, 1), no
  `BigInt` retry; an `OverflowError` surfaces as an internal error. — Gates pass on 1.12 and
  1.13 (section 1.2); speed equals the slice (1.3); overflow is impossible under the cap
  (1.4); it deletes the retry and the only Union. A passes but triples live bytes; B passes
  but keeps the Union return R18 forbids, and its loops box the Union, which a non-escaping
  read would hoist as U's does. — A law above 2⁶³ would error; that needs about 10¹³ steps.
- **L2 (spec line 477).** Amend to "checked `Int64`; under the cap no coefficient leaves
  `Int64` (proof in the PR), so an overflow is an internal error". The solve keeps its
  `Rational{BigInt}` retry (spec line 280; matrices are heap objects, GC-safe). — The spec
  is the authority; Denis approves. — If he declines, take A and accept 3× live bytes.
- **L3 (Q10).** `complexity` by checked-`Int` Bareiss, the `BigInt` determinant only after
  an `OverflowError`; the allosteric cap is V×τ(active) + V×τ(inactive) ≤ 5,000. — Equal on
  2,647 pairs; faster; removes the determinant's live bytes (ROADMAP "BigInt arithmetic
  inflates…", `_eq_complexity` part); overflow only at about 1e10 trees. — About 25 lines.
- **L4 (Q11).** The solve returns a concrete struct: column names, independents, and per
  dependent (column, `keq::Rational{Int}`, `factors::Vector{Pair{Int, Rational{Int}}}`, label),
  narrowed to `Rational{Int}` at the boundary; `_dependent_param_exprs` renders it;
  `all_parameters` and the engine read it. — Byte-identical on 1,974 mechanisms; drops the
  `Expr` interpreter. — One struct and a renderer.
- **L5 (Q9).** Class 2 at the derivation entry stays the user-facing check; the engine adds
  an internal assertion while it walks the RE forest (each non-tree RE edge closes with zero
  net uptake per name). — Without it the engine returns a finite wrong law on RE turnover
  (section 4); about 5% of law time standalone, less fused. — About 10 redundant lines.
- **L6 (Q14).** Class 5 ships as phase 4 specified, on the new law's D, fA, fI. — Exact on
  the hand cases, 0.4 µs. — None found.
- **L7 (R18).** The phase-5 PR reruns this scan on the package's own workloads,
  `--pkgimages=no`, 1.12 and 1.13, with `repro_upstream.jl` as the positive control.
  — The scan covers only compiled code; the control shows it still detects (section 0).
  — About 2 minutes per version.

## 7. Limits

- Prototype on slice types; the package port must rerun the scan (L7).
- The identify workload fits with a loss stand-in, not an optimizer; it exercises the key
  cache, beam per key and CV deduplication, not CMA-ES or Distributed workers.
- Rounds differ between 1.12 and 1.13 (code text 1,889,258 vs 1,889,261 characters; 10 vs
  14 children), so each version's three rounds are compared only with each other.
- The 70.7 µs population was lost; section 1.3 compares against U on the same laws instead.
- Reproduce: `law/variants/mk.py`; `law/scripts/run_scans.sh`, `run_scans2.sh`,
  `run_rest.sh`, `locate_U.jl` through `phase45/jl_research.sh`.
