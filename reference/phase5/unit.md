# R5.4 and engine.md Q2-Q7, Q12, Q13: unit rule and printing

Paths are relative to `scratchpad/phase45/research5/unit/` (scripts, `data/`, `logs/`). SPEC =
`phase3/proto3_wt/docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`. SUM =
`logs/summary.log`. Prototype on the slice types (`slice/src/`); section 4 maps it onto phase 3.

## 0. Answer

- R5.4 passes: under the spec's unit rule every printed denominator term is dimensionless in
  all 9,406 laws, in the full law and in the coefficient form (SUM:22, 60, 102, 154).
- Printed laws equal the derived laws to at most 3.1e-152 (512-bit) in both layouts; reference lines
  match to at most 4.7e-154 (SUM:24-31, 62-69, 104-111, 156-163).
- The spec layout prints 0 artifacts. The slice layout prints `* 1 / (…)` in 8 ping-pong laws
  (SUM:114-115; `logs/artprobe.log`).
- The spec rule's N differs from Cleland's in 264 of 382 Cleland-eligible conformations: 3 of 5
  specs, 4 of 85 ping-pong and 257 of 292 ter-ter (SUM:10, 89, 133).
- Cause: a lower-degree substrate term (`C`, `D`) enters the graded basis before Cleland's
  monomials in 237 (all 3 specs, 234 ter-ter); a square (`A²`) does in 27 (`logs/why.log:1-6`).
- No population law needs fractional λ or the ratio form (SUM:20-21, 58-59, 100-101, 152-153).
  Package-accepted hand fixtures exercise both (sections 2-3).
- Main proposal: restore Cleland-first (U2). It reverses Denis's ledger item 11 and needs his
  call. Rulings U1-U9 are in section 3.

## 1. Inputs and method

- Environment: `env/` develops `graph_mechanism_wt` at 7205ba4, whose src is #92 (census mkenv
  recipe; `logs/mkenv.log`).
- Populations were regenerated with the package as census neutral text (`export_pops.jl`):
  - 42 specs (`logs/export_small.log:2-3`).
  - bi-bi (census reaction, competitive A and I): init 239 plus every level-1 child gives 5,671
    distinct mechanisms (`logs/export_small.log:5-6`).
  - ping-pong (census reaction: R, n = 2, shared site): init 191 plus level 1 gives 1,828
    (`logs/export_small.log:8-9`).
  - ter-ter: a seeded sample of 2,000 from the 250,855 init mechanisms (34 s, 2.2 GB;
    `logs/export_terter.log:2-3`).
- Laws: the slice derives each law (`Slice.solve`, `rate_law`). 135 ter-ter mechanisms hit the
  5,000-product cap, leaving 9,406 laws, 2,778 of them allosteric (SUM:4, 45, 83, 125-127).
- Library: `ulib.jl` comes from `research5_interrupted/ulib.jl`. I read and checked it and fixed
  one bug: in `allo_texts` the nested `build` closure rebound the outer `defs` and `legend`, so
  the full allosteric law printed coefficient names (now `local`; ulib.jl:476-477, 511-515).
- Rules per conformation (ulib.jl:194-222):
  - the spec rule: constant term, else the graded basis with fractional λ, else the ratio form;
  - five order variants (ulib.jl:176-189);
  - the slice rule (`slice/src/normalize.jl:17-40`: constant, Cleland, integer affine,
    coefficient);
  - Cleland's N.
- N equality is compared on log N at two random points to 1e-100 (ulib.jl:228-234).
- Checks per law (ulib.jl:707-748):
  - value: the printed text is parsed and evaluated in BigFloat against `E_total·H·Q/D`;
  - units: every top-level denominator term must be (0, 0), and v must be E_total/time;
  - reference lines are checked against `all_parameters`;
  - artifact regexes (ulib.jl:649-653).
- Cost: peak 2.26 GB (ter-ter export). Census time per law: median 0.011-0.153 s, max 13.6 s
  (SUM:41, 79, 122, 173).

## 2. R5.4 census

| Population | laws | conformations | constant | graded λ (integer) | fractional λ | ratio | unit failures | max value error | artifacts (spec layout) |
|---|---|---|---|---|---|---|---|---|---|
| specs | 42 | 53 | 48 | 5 | 0 | 0 | 0 | 4.4e-154 | 0 |
| bi-bi | 5,671 | 7,063 | 7,063 | 0 | 0 | 0 | 0 | 2.9e-152 | 0 |
| ping-pong | 1,828 | 2,593 | 2,423 | 170 | 0 | 0 | 0 | 2.4e-153 | 0 |
| ter-ter | 1,865 | 1,865 | 1,141 | 724 | 0 | 0 | 0 | 3.1e-152 | 0 |

Sources: SUM:2-41, 43-79, 81-122, 124-173.

Where N differs from Cleland's (no constant term, all of Cleland's monomials present):

| Population | eligible | spec rule = Cleland | differs | cause (`logs/why.log`) | slice rule (Cleland-first) differs from spec rule |
|---|---|---|---|---|---|
| specs | 5 | 2 | 3 | lower-degree term C: λ = C:2, AB:1, AC:−1, BC:−1 (SUM:11-13) | 3 (SUM:14) |
| ping-pong | 85 | 81 | 4 | power: λ = A:2, A²:−1 (SUM:90-93) | 4 (SUM:94) |
| ter-ter | 292 | 35 | 257 | 234 lower-degree term (D:2, AB, AD, BD), 23 power (why.log:5-6) | 257 (SUM:146) |

- The 3 specs are Segel Bi Uni Uni Uni PP Ter Bi, Bi Uni Uni Bi PP Ter Ter and Bi Bi Uni Uni PP
  Ter Ter (SUM:11-13). These are the same 3 of 40 that C4 found (`followups/results/C4.md:77`).
- Why: E carries a degree-1 term `d_C·C` (`logs/examples.log:20`). SPEC:502 takes the lowest
  degree first, so the basis {C, AB, AC, BC} solves the units equation before ABC enters.
- Spec rule's N: `(k_EAB_to_F_P + k_EAB_to_EA_B) * (k_FC_to_E_Q + k_FC_to_F_C)`
  (`logs/examples.log:21`).
- Cleland's N: the same divided by `(k_FC_to_E_Q + k_EAB_to_F_P)^2`, which gives Segel's
  Ka·Kb·Kc form (`logs/examples.log:25`).
- Ping-pong power case (6e46f55b0b68e5d9): A² precedes AB in metabolite order. The spec rule's N
  is Cleland's polynomial squared (`logs/examples.log:36, 40`).
- Ordered bi-bi forms (Segel Ping Pong Bi Bi, Hexa Uni Ping Pong) agree with Cleland (SUM:10;
  `logs/examples.log:4-10`).
- Order variants: none reaches Cleland in the 3 specs (SUM:15-19). Squarefree-first and C4's
  text order reach it in 4 of 4 ping-pong cases (SUM:95, 98) but only 37 of 292 ter-ter
  (SUM:147, 150).

Fixtures (slice `fixtures.jl`; the package accepts all four, `logs/pkgfix.log:4-27`):

| Fixture | Rule | Printed | Value, units |
|---|---|---|---|
| F1: `E + A <--> F + P`, `F + A <--> E + P` | ratio, divide by A | `(n_0 + n_P_A⁻¹ * P / A) … / (1 + d_P_A⁻¹ * P / A)`, `N = k_E_A_to_F_P + k_F_A_to_E_P` (fixtures.log:9-16) | 4.8e-154, pass (fixtures.log:4) |
| F2: F1 + `F + A + A ⇌ F(A, A)` | fractional, λ = A:3/2, A³:−1/2 | `N = (k_E_A_to_F_P + k_F_A_to_E_P)^(3/2)`, `sqrt(…)` in coefficients (fixtures.log:20-37) | 3.4e-154, pass (fixtures.log:21) |
| F3: RE Theorell–Chance `E(A) + B ⇌ E(Q) + P` | ratio, divide by B | `v = E_total * (n_B⁻¹ / B) * (A * B - P * Q / Keq) / (1 + d_P_B⁻¹ * P / B + …)`, `n_B⁻¹ = k_E_A_to_EA` (fixtures.log:48-51) | 2.5e-154, pass (fixtures.log:41) |
| F4: F1 + `F(A)`, `F(A, A)` | integer, λ = A:2, A²:−1 | `N = (…)^2` (fixtures.log:64) | 2.6e-154, pass (fixtures.log:59) |

- The slice's coefficient fallback leaves `A` and `B` with units on F1-F3 (fixtures.log:4, 21,
  41). This is SEL's "2 laws keep dimensional terms", and the ratio form fixes it.
- `d_P_B⁻¹`, `n_B⁻¹`, `d_A⁻²_B` and `d_A²_B` are valid identifiers and parse as names
  (`logs/smoke.log:1-10`).

## 3. Proposed rulings (ruling — why — cost if wrong)

**U1 (R5.4). The unit rule ships as in SPEC:490-514, with U2's amendment.**
- Why: 0 unit failures in 9,406 laws plus 4 fixtures; the units follow from homogeneity, since
  Σλ = 1 and Σλ·a = 0 give N the unit of D.
- Cost if wrong: none found.

**U2 (Q4). Restore Cleland-first, as the slice does. Denis decides, because this reverses ledger
item 11's choice of the graded basis.**
- Rule: when ΠS/X for every substrate X and ΠS are all denominator terms, N is Cleland's.
  Otherwise use the graded basis.
- Why:
  - Dropping Cleland-first changes N in 264 of 382 eligible conformations: 3/5 specs, 4/85
    ping-pong and 257/292 ter-ter, which is 13.8% of all ter-ter laws (SUM:133).
  - Ledger 11 rested on C4's 37 of 40 bi-bi-scale laws (C4.md:77). Ter-ter overturns that scale.
  - SPEC:38-39 promises a "Cleland-style coefficient form", and `deriving/ping_pong.md` teaches
    Cleland's K's.
  - The slice does it in 4 lines (`slice/src/normalize.jl:20-23`).
  - Values never change.
- Cost if wrong: one extra rule and its test. Dropping it later re-pins N in the affected prints.

**U3 (Q3). Keep SPEC:502-503's order.**
- Order: degree; then distinct non-substrate count; then metabolite order, read as the
  lexicographic exponent vector over `metabolites(m)` (substrates, products, modifiers, each
  sorted).
- Why:
  - Under U2, alternative tie-breaks change N only in 8-17 of 432 ter-ter graded-basis
    conformations: text 8, squarefree 8, multiplicity count 10, no non-substrate key 17
    (why.log:9). They change none in specs or ping-pong (why.log:7-8).
  - None has a ground truth.
  - Metabolite order matches the `concs` line and phase 3's `metabolites`. C4's text order
    depends on spelling.
  - If Denis declines U2, add "fewest copies of one metabolite" after the non-substrate count. It
    fixes the 4 ping-pong power cases (SUM:98) and nothing else of note (37 of 292 ter-ter,
    SUM:150).
- Cost if wrong: one sort key, and N moves in at most 17 ter-ter laws.

**U4 (Q5). Keep fractional λ.**
- Why:
  - No population law needs it, but F2 is a package-accepted mechanism that does
    (pkgfix.log:10-15).
  - The solve is exact rational anyway.
  - Roots print through the existing `power_text` (`sqrt`, `^(3/2)`).
  - Integer content raised to a fractional power goes into N (ulib.jl:327-358).
  - F2 is the test fixture.
- Cost if wrong: about 5 lines kept for a case no enumerated mechanism reaches.

**U5 (Q6). Ratio form as in SPEC:506-509.**
- Rules:
  - Divide numerator and denominator by the first term in U3's order, its concentration
    monomial included.
  - Name coefficients by the shifted monomial: positive exponents first, then negative, each in
    metabolite order (`n_B⁻¹`, `d_P_B⁻¹`, `d_A_P_B⁻¹`).
  - Print no `n_0` when the numerator has no term at the origin.
  - Print `N` only when the first term's coefficient has a polynomial factor (F1 has one, F3 does
    not).
- Fixtures: F1 and F3.
- Why: both pass values and units, while the slice fallback fails units (section 2).
- Cost if wrong: names change on a class that no enumerated law reaches.

**U6 (Q2). Print the coefficient form with `rate_equation_string(m; coefficients = true)`.**
- The default stays the full law, and identify's CSV keeps the full law only
  (`identify_rate_equation.jl:254`).
- Why:
  - The coefficient form is about as long as the full law, so carrying both doubles the CSV
    cell. Median lengths: specs 881 vs 898, bi-bi 1,353 vs 1,341, ping-pong 814 vs 841, ter-ter
    4,993 vs 4,534. Maxima: 406,617 vs 394,061 (SUM:36-37, 74-75, 117-118, 168-169).
  - A keyword adds no export, so no new docstring page is needed (SPEC:740-741).
- Cost if wrong: moving a keyword into a function later.

**U7 (Q7). Allosteric names follow the slice.**
- Per conformation s ∈ {A, I}:
  - `Q_s`: numerator;
  - `D_s`: partition function;
  - `F_s`: flip-form weight, printed only when F_A ≠ F_I and neither is 1;
  - `W_s`: site factors;
  - `N_s`: the normalizer's polynomial factor.
- Coefficients are `nA_0`, `dA_S`, `fA_…`. In v the factors of 1 are omitted
  (`logs/smoke.log:121-144`).
- Why:
  - `N` already means the normalizer in plain laws (SPEC:511), and `n_`/`d_` are the plain
    coefficient prefixes. The docs' `N_A` means numerator, which would give N two meanings.
  - Every name holds `_`, which E5 bans from metabolite names, so none can clash.
  - 0 value errors above 1e-100, 0 unit failures and 0 artifacts in 2,778 allosteric laws.
- Cost if wrong: renames in the printer, the golden file and `deriving/mwc_allostery.md`, which
  phase 5 rewrites anyway.

**U8 (Q12). `v = E_total * (Q) * (H) / (D)`.**
- Q comes first, as in SPEC:521. It is parenthesized unless it is one name (`n_0`) and dropped
  when it is 1. H and D are always parenthesized.
- Why:
  - The slice's `E_total * (H) * Q / (D)` prints `* 1 / (K_EA_to_E_A * …)` in 8 ping-pong laws
    (artprobe.log), and this layout prints none.
  - Both layouts are value-equal (SUM:24-27).
  - The TR:683-700 extractor reads Q, H and D as three groups.
- Cost if wrong: re-pinning the 11 factored pairs again.

**U9 (Q13). Reference lines print as a ratio of products.**
- They use the term printer that v and the legend use: `k_E_P_to_ES = k_ES_to_E_P *
  k_E_S_to_ES / (k_ES_to_E_S * Keq)` (`logs/smoke.log:59`).
- Exponents print exact (`sqrt(x)`, `x^(p/q)`), a tie to 1 prints `= 1`, and each section is
  sorted by name.
- Why:
  - Both styles are value-equal: 4.7e-154 against 3.9e-154 (SUM:30-31, 68-69).
  - Today's style prints Float exponents (`x ^ 0.5`).
  - Every population dependent has an integer exponent (SUM:32, 70, 112, 164), so only
    hand-written laws show roots.
- Cost if wrong: re-pin TR:1461 (`0.5` becomes `sqrt`) and TR:1563 (`1.0` becomes `1`).

## 4. Mapping onto phase 3's types

- The unit rule reads only the law, the substrates and the concentration order, so it ports as
  `normalizer(D, vars, substrate indices)`:
  - D comes from the law split by concentration monomial (`Slice.split_conc`).
  - The order comes from phase 3's `metabolites(m)` (S3_design.md:932-937). Site ligands append
    last but never enter a conformation's D.
  - The slice's `metabolites` differs only by sorting regulators together with modifiers
    (`slice/src/types.jl:299-300`).
- λ: the RREF has at most (number of metabolites + 1) rows. `Rational{Int}` fits, and its
  overflow throws.
  - R18: `Norm` is one concrete struct, and the normalizer returns it, never a Union.
- Printer inputs:
  - the column names, through the naming chokepoint (`_edge_names`, S3_design.md table §7);
  - the structured solve output (engine.md Q11: `Dep(col, terms, keq, label)`) for reference
    lines;
  - `parameters(m)`.
- Allosteric: the conformation tags A and I. The flip-form weight is `forms[1]` = free E (phase
  3 stores free E first, S3_design.md:17, 126; phase 4 flips at E).
- Unit oracle for tests:
  - a rate constant's unit is conc^(−copies taken up) per time, from `total(e.consumed)` or
    `total(e.released)`;
  - an RE binding K is a concentration, and another RE K is conc^(released − consumed);
  - this follows ulib.jl:275-295.

## 5. Notes

- Length against today's printer on the 42 specs: today 82,422 characters, new full law 70,011,
  coefficient form 66,921. The new law is shorter for 29 specs. Extremes: LDH 6-group 15,300 →
  6,164; Random-order Bi-Bi 9,387 → 12,016 (`logs/examples.log:48-89`).
- The ter-ter full law reaches 394,061 characters (SUM:168). The CSV-size issue predates phase
  5, and the key map should note it.
- Files:
  - code: `ulib.jl`, `pops.jl`, `export_pops.jl`, `census.jl`, `run_census.sh`,
    `summarize.jl`, `fixtures.jl`, `pkgfix.jl`, `examples.jl`, `why.jl`, `artprobe.jl`,
    `smoke.jl`;
  - data: `data/*.tsv` (neutral text) and `data/*.jls` (census rows).
