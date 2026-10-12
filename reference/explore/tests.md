# The test suite as the specification and safety net for a derivation rewrite

All facts come from reading the files at 277a785. I ran no Julia. Counts are grep-based static counts and are approximate where marked "≈". Statements marked "(inference)" are my reading of the code, not something I measured.

## 1. Test inventory: what each file pins

`runtests.jl:8` includes `mechanism_definitions_for_test_enzyme_derivation.jl` before the top-level testset. That fixtures file is therefore a single point of failure. Its source-capture macros call `EnzymeRates._parse_mechanism_body` directly (`fixtures:64-73`), so a change to that parser stops the whole suite.

| File (testsets / static `@test` lines) | What it pins |
|---|---|
| fixtures file (0/0) | The 42 specs (section 2), `_testhelper_allo_from_source` (83-96) and `_testhelper_random_reduced_params` (144-152). |
| `test_rate_eq_derivation.jl` (47/139) | The per-spec harness `_testhelper_run_all_tests` (1103-1119) over all specs; the spec count is pinned at 1125. See the list below. |
| `test_allosteric_golden.jl` (1/6) | Byte-identical Reduced strings for the 11 allosteric specs against `test/reference/allosteric_golden_reference.txt` (22 lines). Also `parameters(m, Reduced) == (fitted..., :Keq, :E_total)` and the fitted order: sorted catalytic names, then regulator constants, then `L` (41-49). The golden file contains the goal-14 artifacts `1 ^ 4` (PFK-1, line 8) and `1 ^ 2`/`1 ^ 4` (PK, line 12). |
| `allosteric_ground_truth.jl` (14/58) | Mass-action networks (oracle O6 in section 3). Each gate also pins an exact ordered `fitted_params` tuple (92, 149, 235, 553-558, 673-676, 727, 941-942, 1006-1010) and internal free-enzyme weights `d_free` (120-124, 776-777, 1015-1016). |
| `test_allosteric_collapse.jl` (10/27) | Which inactive-state splits collapse or stay free (literal names), mirror lines in the string, `veq≈0`, and identifiability by perturbation. Builds mechanisms from raw `Step`/`Species`/`EnzymeReaction` constructors (9-31). |
| `test_types.jl` (72/411) | Constructors, canonical form, `==`/`hash`, 49 `@test_throws` validations (bottomless RE segments, kinetic-group rules, the `:OnlyA` Stiemke guard at 2190-2453). Exact names through the chokepoint (1091-1217, 1536-1695), `_sig_of` round trip (1122-1174), the AST chokepoint guard (2175-2188), `_rational_nullspace`/`_solve_dependent_set` (2455-2517), and naming-cache/serialization identity (2519-2578). |
| `test_dsl.jl` (35/136) | Derivation-facing parts only: `::Inh` keeps the concentration name `P` in `(; S, P) = concs` (882-903), and dual-role binding types (905-952). |
| `test_compile_budget.jl` (5/8) | Two fresh subprocesses. `init_mechanisms` trace-compile ≤200 (line 15). First `rate_equation` call ≤6 s (16). Warm/cold reuse ratio <0.01. Uni-uni and ter-ter dispatch to the same method. |
| `test_accessors.jl` (1/12) | `Full` vs `Reduced` parameter lists; `parameters(allo, Full)` raises `MethodError` (34). |
| `test_pivot_priority_regression.jl` (1/2) | A raw `EnzymeMechanism{…}` Sig literal; fitted count = 6; one of two tied names is absent. |
| `test_fitting.jl` (17/47) | A literal uni-uni `true_params` (22-23) used by 7 testsets. `fitted_params == parameters` minus `Keq`/`E_total` (66-72). `loss!` allocates nothing (203-245) and takes <70 µs for 500 points (301). |
| `test_identify_rate_equation.jl` (41/295) | A literal allosteric `true_params` (59-63). `eq_hash` = hash of `rate_equation_string` with comment lines and `ANNOTATION_SUBSTITUTED` lines removed (`src/identify_rate_equation.jl:262-268`). 239 classes over bi-bi `init_mechanisms` (692-734). Twins share key and fitted tuple (1230-1239). Split/merge parameter-count difference (1365-1415). The `mechanism_type` column round-trips through `Core.eval(Meta.parse(…))` (250-260, 350, 418-432). |
| `test_mechanism_enumeration.jl` (237/962) | `_testhelper_identifiable_rank` (72-99) used 95 times. `_testhelper_fitted` used 93 times: 27 literal count pins, 28 `fitted == rank`, 10 Δ-param pins. `_independent_param_count == fitted` on all 42 specs (6558-6566). `_hyperbolic_catalysis` checked against the derived denominator degree over >100 mechanisms (8480-8554). Division-freeness (11165-11183) and a structural init golden (11137-11163). Mechanism equality through 68 `Set(…) == Set(…)` asserts and 95 `_testhelper_on_reaction` rebuilds (152-153). |

The standalone testsets in `test_rate_eq_derivation.jl` pin:
- The BigFloat mass-action oracle on 25 further mechanisms (1264-1328, 1693).
- Byte-identical strings after reversing steps or group order (1330-1376).
- `allocs == 0` and `t < 120e-9` per call (874-888).
- Expression shape: only 2-argument `+`/`*` calls (1503-1529).
- At most 6 parentheses in the Full-mode string (1597-1613).
- Exact lines for folded ties, `= 1.0`, and residual names (1531-1595, 1675-1700).
- One byte-identical allosteric string (2012-2032).
- `_eq_complexity` values 25, 336, 49 and 9 (1789-1852).
- The "no finite rate" and "polynomial terms" errors (1744-1787, 1854-1882).
- Compile time under 20 s (1706-1742).
- kcat contracts (2205-2360).
- Literal fitted tuples and names (1240-1262, 1615-1646, 2034-2094, 2143-2176).

**Two tests pass without testing anything (inference).** `:OnlyA binding group: no K_T param` and `:OnlyA iso group: no kf_T/kr_T param` (enumeration 6487 and 6505) look for a `_T` suffix. Current names put the state in the middle (`K_I_…`), and `test_rate_eq_derivation.jl:2127-2129` asserts that no name ends in `_T`.

## 2. MECHANISM_TEST_SPECS anatomy

**Fields** (`fixtures:13-57`):
- `name` and `mechanism`.
- `metabolite_names`: pins the order of `metabolites(m)`.
- `expected_n_states` and `expected_n_steps`.
- `expected_n_haldane_constraints` (the dependent expression's right side mentions `Keq`), `expected_n_mirror_constraints` (the right side is a single Symbol; the file itself calls this vestigial at 26-27) and `expected_n_wegscheider_constraints` (any other right side).
- `expected_n_independent_params`: equals `length(fitted_params)`. The comment at line 29, "2*n_steps - n_constraints", is wrong for RE steps; RE Uni-Uni pins 2.
- `run_ode_test`, `reference_rtol` (1e-8) and `ode_rtol` (1e-6).
- `analytical_rate_fn` and `analytical_kcat_fn`.
- `expected_factored_num` and `expected_factored_denom`: exact strings.
- `source_steps` and `source_reg_sites`: the as-written order. Positional oracles number constants `k1f, k2f, …` in written order, and the constructor reorders steps, so the bridge needs it.

**The 42 specs** (count pinned at `test_rate_eq_derivation.jl:1125`):

| Family | Specs |
|---|---|
| Uni-uni / iso | Uni-Uni; Segel Uni Uni (IX-8); Segel Iso Uni Uni (IX-45); RE Uni-Uni |
| Ordered | Segel Ordered Uni Bi (IX-60), Bi Bi (IX-87), Ter Bi (IX-195), Ter Ter (IX-261); RE Ordered Bi-Bi |
| Theorell–Chance | Segel TC Bi Bi (IX-122) |
| Ping-pong (conformation labels F/G, no residuals) | Segel IX-140, IX-228, IX-278, IX-288, IX-308 |
| Random order / Wegscheider / RE-SS mixes | Random-order Bi-Bi (all SS); RE Random Bi-Bi; three "Numerator:" specs |
| Inhibitors / activators | 11 specs, each with exact factored strings |
| Allosteric MWC | MWC Dimer; Homodimer + NC inhibitor; MWC Tetramer; PFK-1; HK; PK; m_all; m_OnlyA_prod; LDH i-state ×3 |

Totals:
- 31 non-allosteric and 11 allosteric.
- 14 all-SS; 14 with multi-step kinetic groups.
- 34 with an analytical oracle, 11 with a kcat oracle, 11 with factored strings.
- 11 skip the ODE check, and all 11 are the allosteric specs.
- 36 are built with the source-capture macros.

What the specs do not cover:
- Covalent-residual forms (`E(; residual = A - P)`).
- Fused bindings into a merged complex (`E(A) + B <--> E(P, Q)`).
- Two-metabolite step sides (`E + A + B`).

These appear only in the mass-action cases and in `allosteric_ground_truth.jl`.

**How each spec is checked** (1103-1119):
1. Structure (687-700): numbers of forms and steps, metabolite order, all parameter names unique.
2. Constraints (765-793): the dependent graph is sound (730-763), the dependent expressions are classified by the shape of their right side, and the independent count matches.
3. Reference QSSA for the 14 all-SS specs: 20 draws, rtol 1e-8.
4. Analytical rate for 34 specs: 20 draws, rtol 1e-10.
5. Haldane equilibrium for all 42: `|v| < 1e-10`.
6. Performance for all 42: 0 allocations and under 120 ns.
7. Rate-equation string (916-952): 8 format `occursin` checks; the string is evaluated as Julia and must match `rate_equation` (10 draws, rtol 1e-10); exact factored numerator and denominator for 11 specs.
8. Zero-metabolite finiteness (1081-1097): `v` stays finite, and stays nonzero when a substrate or product is zeroed.
9. ODE steady state for the 31 non-allosteric specs.
10. Analytical kcat for 11 specs.
11. kcat rescaling for all 42 (998-1068).

Draws are uniform on 0.1–10.

## 3. Oracle inventory

| Oracle | How it is invoked | Coverage | Independent of names? Of the dependent-parameter choice? |
|---|---|---|---|
| O1 Textbook / hand closed forms | `_testhelper_test_analytical_rate` (818-845), through the bridge `_testhelper_analytical_oracle_params` (279-451). The bridge maps chokepoint names to positional `K1/k1f/k1f_T/K_<lig>_reg<n>`, renumbers to written order, and swaps f/r for releases that are stored as the reverse binding. | 34 specs: 12 Segel, 11 inhibitor, 2 Cha-style RE, Uni-Uni, 8 hand-derived MWC | Names: yes, through `name(p,m)`. Dependent choice: **no**. The dependent constants come from the package's own `_dependent_param_exprs` (485-489), so a wrong constraint expression passes. Exceptions: PFK, HK, PK, m_all and m_OnlyA_prod recompute k_r from `Keq` (fixtures 1919, 2039, 2154, 2238, 2304). |
| O2 Reference QSSA | Laplacian cofactors (18-97) | 14 all-SS specs | Same caveat as O1; one metabolite per step side |
| O3 ODE steady state | RadauIIA9 at 1e-13 (628-672); RE steps become 1e6 rates (530-577) | 31 non-allosteric specs | Same caveat; RE steps only to rtol 1e-3 |
| O4 BigFloat mass-action at thermodynamically consistent points | Builds every constant from random free energies, with kinetic groups projected (108-154); 512-bit linear solve with RE at Λ=1e60 (158-185); feeds the package only the fitted subset (187-201) | 25 mechanisms: fused steps, TC, `E+A+B`, residual ping-pong, dead ends, dual role, iso pinned to 1 (1138-1217, 1275-1312, 1693) | **Yes and yes.** The only oracle that checks the constraint solution. Non-allosteric only, and not run on the specs. |
| O5 Haldane equilibrium | 847-872 | 42 specs | Yes and yes, but it only checks the overall equilibrium point |
| O6 Two-conformation mass-action networks | `_testhelper_mwc_ground_truth_flux` (ground truth 16-27) on hand-built networks, formulation 1, FAST 1e7–1e9; the harness is self-validated (315, 470, 643, 844) by L=0 limits, identical conformations, FAST invariance and Segel IX-140 | OnlyA uni/bi-uni, metabolite-bearing D, NonequalAI at 1–3 protomers, NonequalAI ping-pong (rtol 1e-10), an RE segment holding E and F, an inert inactive state, dead-inactive ping-pong with an L=0 twin | Values yes; the gates pin literal fitted tuples |
| O7 Detailed-balance sweeps | Reproducers (1464-1499), collapse `veq`, `m_ro` (1965-2008), ground truth 1066-1068 | LDH-like and RO bi-bi tag mixes | Yes and yes |
| O8 Finite-difference identifiable rank | Enumeration 72-99 | 95 uses | Yes and yes; pins counts |
| O9 String evaluates to the same value as the code | 679-683 | 42 specs | Self-consistency only; impossible for residual names |
| O10 Invariances | Step reversal (1330-1376), tie equals group (2143-2176), eq_hash partition (identify 692-734) | Specs plus 25 cases | Self-comparisons |
| O11 kcat | Closed forms, numeric vmax sweeps (2205-2342), rescaling | — | Goes through internal `_kcat_forward` |
| O12 Structural cross-checks | `_independent_param_count` (enumeration 6558), `_hyperbolic_catalysis` (8480), `_eq_complexity`, undefined-symbol scan (1415-1456) | — | Coupled to internals |

## 4. Breakage matrix (counts are testsets)

| Change | What breaks | ≈ count |
|---|---|---|
| **(a1)** Residual names without `+`/`-` | Ground-truth gates 542, 660 and 991 (51 literal lines); "prints residual names bare" (1577; intended); types 1661; enumeration 8738 only if species names change too. No spec has residuals. `eq_hash` values move (ROADMAP). | 5–6 |
| **(a2)** Any general renaming scheme | Golden file (11 asserts); 11 factored-string specs; `test_rate_eq_derivation.jl` 10 (1240, 1314, 1531, 1577, 1615, 1675, 1886, 2012, 2034, 2143); `test_types.jl` 12; ground truth 9; collapse 6; fitting 7; identify 1; pivot 1; accessors 1; enumeration 2; 3 doctests. O1–O5, O8 and O10 survive because they key through the chokepoint. | ≈ 61 |
| **(b)** A different, equivalent dependent-parameter choice | The Constraints counts of the 15 specs with ≥2 dependents (see note). Also: 11 factored strings (they print the dependent `k_EP_to_ES`), golden, byte-identical, folded ties, `=1.0`, fused-release policy (1240), single-symbol tie, 9 ground-truth tuples, collapse 6, types "tied steps" (`:K_F_to_E in fitted_params`), fitting 7, identify 2, doctests. Invariant: `expected_n_independent_params` and all counts. | up to ≈ 55 |
| **(c)** Normalized or dimensionless form, names kept | No numeric test breaks, provided no `1/[metabolite]` factor appears. The zero-metabolite checks (42 specs) and division-freeness (11165) enforce that, so a ping-pong form must divide by parameter monomials, not concentration monomials. Internal `d_free` pins break (1219-1238, ground truth 120-124, 776-777, 1015-1016, 2372-2373; goal 7 removes them anyway), plus the parenthesis budget and `_kcat_forward` (18 call sites). The factored strings are RE and already `1 + …`, so they likely survive (inference). | 5–8 |
| **(d)** Different string layout | All 42 string testsets: 8 checks each, string evaluation needs valid Julia with the `(; …) = params`/`concs` lines, and the numerator/denominator extraction needs `E_total * (` on the last line. Also factored strings, golden, byte-identical, folded-tie lines, the `# Wegscheider constraints:` header, residual regex, §5a parser, edge cases, collapse ≈5, types ≈4, DSL 1, doctests 3. `eq_hash` depends on the comment and annotation-line convention. | 42 + ≈ 25 |
| **(e)** Remove `Full` / rename `Reduced` | `Full`: accessors, 1326-1327, 1611, types 1547 and 1560-1662. Those name pins use `parameters(m, Full)` as the "every constant" list, so they need a replacement accessor. **No test evaluates Full-mode `rate_equation` numerically.** `Reduced`: accessors, golden, 2202, 2231, 2338. All mechanical. | 9 |
| **(f)** New Mechanism representation or constructor signature | Call sites: `Step` 97, `Species` ≈78, 2-arg `Mechanism` ≥32 plus 95 `_testhelper_on_reaction`, 5-arg `AllostericMechanism` ≥38, `EnzymeReaction` 30. Testsets: types ≥33 of 72; enumeration ≥50 directly, plus every mechanism set-equality through `==`/`hash`; identify 5 plus 3 Sig round trips; all collapse tests; the reversal test; the 42-spec bridge (`Step ==`, `_canonicalize_step_directions`, `_group_rep`, `_free_enz_set`, the stored-direction swap); pivot; `_sig_of`; serialization. Dropping the `reaction` field (goal 8) hits fixtures, O4 and `_testhelper_on_reaction`. | ≥ 100 |
| **(g)** One DSL syntax | 224 `@allosteric_mechanism` call sites in 9 files. Aliases: `dead_end_inhibitors` 36 / `competitive_inhibitors` 17; `oligomeric_state` 67 / `allowed_catalytic_multiplicities` 11; `regulators` ≈150 / `catalytic_inhibitors` 36. `::Inh` has 283 uses (goal 16). Changing `_parse_mechanism_body` breaks the 36 source-captured specs, and with them the whole suite. `test_dsl.jl` has 35 testsets with 59 exact-message `@test_throws`. | DSL file + fixtures |

**Note on (b): why the constraint counts are not invariant.** The Haldane/Wegscheider/mirror split depends on the chosen basis and on which single-symbol ties are folded into renames. Non-essential Activator pins 2 Haldane and 0 Wegscheider, because K4=K1 and K6=K3 are folded. A cycle-type classification would instead give 1 reaction cycle and 3 net-zero cycles (inference).

**Internal-API unit tests will be rewritten or deleted along with their functions:**
- `build_power_expr` (1653-1673) and the POLY helpers (2362-2374).
- `_num_den_exprs` (1511).
- `_raw_symbolic_rate_polys` and `_state_rate_polys` (rate_eq 1236; ground truth ×4; enumeration ×4).
- `_dependent_param_exprs` (5 uses).
- `_rational_nullspace`, `_has_strict_positive_combination` and `_solve_dependent_set` (types 2455-2517).

## 5. Recommended verification strategy

**What stays fixed and what gets regenerated:**
1. **Keep unchanged (semantic contract):** numeric agreement with O1–O7, Haldane equilibrium, finiteness and the no-`1/[metabolite]` rule, kcat values, the allocation-free sub-120 ns contract and `loss!` speed, the expression-shape guard, the compile budgets, and the O10 invariances.
2. **Keep (counts that do not depend on names or basis):** `expected_n_independent_params`; `fitted == identifiable rank` wherever it holds today; `_independent_param_count == fitted`; numbers of forms and steps; 239 eq_hash classes; twins with equal fitted tuples (goal 5).
3. **Regenerate after one hand review:**
   - the golden file, which changes anyway under goal 14;
   - the byte-identical allosteric fixture and the factored strings;
   - the literal fitted tuples and `true_params`;
   - the chokepoint name pins and the 3 doctests;
   - the Haldane/Wegscheider/mirror counts, redefined as invariants: Wegscheider = rank of the net-zero cycles after kinetic-group identification, Haldane = the rest.

**Strengthen the harness before writing the new derivation:**
- Run O4 on the 31 non-allosteric specs. It is cheap, and unlike O1–O3 it does not depend on the package's dependent-parameter choice.
- Generalize O6 so it builds the formulation-1 network automatically from an `AllostericMechanism` and its tags. Then run it on the 11 allosteric specs; the three LDH i-state specs currently have no value oracle at all.
- Re-key the oracles and ground-truth gates by role (for example, "the constant of the step binding S to E") instead of literal names, so a rename touches only one naming-pin testset.
- Add `@test_broken` acceptance tests for the three ROADMAP known issues; the redesign should flip them. Retarget or remove the two empty `_T` tests.
- Add property tests for goals 5 and 9. Random step reversal, group reordering, conformation relabeling, and thermodynamically equivalent groupings must all give identical strings.
- Add a codegen and memory budget for goals 11 and 18: compile N distinct mechanisms in a subprocess and bound the growth in resident memory or method instances. No current test guards this.
- Replace the direction-swap logic in the bridge (406-451) with a "chemically forward" step accessor.
- Write a migration script outside `test/`. Over the specs plus `init_mechanisms` and two expansion rounds (bi-bi, ping-pong, uni-bi, allosteric), compare old and new derivations:
  - fitted counts: equal, except where a known-issue fix applies, where the new count should equal the identifiable rank;
  - values at consistent points, with names mapped by role;
  - eq_hash partition sizes.

**Infrastructure constraints:**
- The full `Pkg.test` takes 5 min 39 s and runs 32,779 assertions (277a785 commit message).
- Use `--heap-size-hint=2500M` and run one Julia process at a time. The machine has 7.7 GB and no swap, and the full suite needs about 4 GB free (CLAUDE.md).
- Focused TestEnv runs skip Aqua and JET. The enumeration file alone takes about 2.5 min (memory notes).
- The compile-budget test starts two fresh Julia processes.
- The 120 ns gate is sensitive to machine load; it takes the minimum of 5×10k calls.
- O3 needs OrdinaryDiffEqFIRK. O4's BigFloat solve grows cubically with the number of enzyme forms, which is fine for the specs but slow beyond about 30 forms (inference).