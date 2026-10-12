# Slice report: documented semantics and promises (docs and public API)

**Sources.** I read `docs/src/**`, `README.md`, `ROADMAP.md`, the exported docstrings in `src/*.jl`, and the tests that pin outputs. Rendered example outputs come from the local `docs/build/`, which git ignores (`.gitignore:7`). It was built 2026-10-08 23:59, about 25 minutes before squash commit 277a785, so it is close to current but not guaranteed identical. Items marked **(inference)** are my reading of the code; I did not run them.

## 1. Public API and its documented semantics

**Exports** (`src/EnzymeRates.jl:5-21`, 18 names):
- `EnzymeReaction`, `EnzymeMechanism`, `AllostericEnzymeMechanism`
- `FittingProblem`, `IdentifyRateEquationProblem`, `IdentifyRateEquationResults`
- the three DSL macros
- `Full`, `Reduced`
- `rate_equation`, `rate_equation_string`, `parameters`, `metabolites`, `rescale_parameter_values`
- `fit_rate_equation`, `identify_rate_equation`

`docs/make.jl` sets `doctest = true` and `checkdocs = :exports`. Every export therefore needs a docstring, and every `jldoctest` output is pinned. An `@example` block that errors, for example because a parameter name changed, breaks the docs build.

### `rate_equation(m, concs, params, [mode])`
Docstring at `src/rate_eq_derivation.jl:390-433`.
- **What it returns:** the net quasi-steady-state turnover of the whole mechanism. Rapid-equilibrium (RE) steps collapse to binding constants; steady-state (SS) steps go through the King–Altman determinant.
- **Arguments:** `concs` and `params` are NamedTuples, destructured by name (`_rate_body`, 379-386). Extra keys are ignored and missing keys error. The default mode is `Reduced`.
- **Performance:** the body is generated once as "a single arithmetic expression with no allocations, loops, or matrix operations" (406-408).
- **Concrete-type methods:** methods on `Mechanism` and `AllostericMechanism` lift through `compile_mechanism`. They allocate and are documented as not the hot path (52-59, 415-419).

### Modes (`src/types.jl:1213-1243`)
- **`Full`:** "all 2N microscopic rate constants… No thermodynamic constraints applied. Parameters: all k's + E_total". It takes no `Keq` and exists only for plain `EnzymeMechanism`. The docs repeat that an allosteric mechanism has neither a Full rate equation nor a Full parameter list (`docs/src/deriving/textbooks.md:110-111`).
- **`Reduced`:** dependent constants are substituted from the independents and `Keq`. Parameters are the independents plus `Keq` plus `E_total`.

### `parameters(m, [mode])` (`rate_eq_derivation.jl:18-44`)
- **Reduced** returns `(indep..., :Keq, :E_total)`, with the independents sorted by name (`thermodynamic_constr…jl:548-553`).
- **Allosteric Reduced** orders them as: sorted catalytic independents, then active-state regulator Ks site by site, then inactive-state ones, then `L` (`rate_eq_derivation.jl:943-950, 1003-1004`).
- **Full** returns the raw constants in canonical group order plus `:E_total`. Rendered: `(:K_EP_to_E_P, :K_ES_to_E_S, :k_ES_to_EP, :k_EP_to_ES, :E_total)`.
- **Inconsistent prose:** `getting_started.md:71` calls the Reduced output "the independent parameter symbols", although it also holds `Keq` and `E_total`.

### `fitted_params(m)` is not exported, but users see it
- Docstring at `rate_eq_derivation.jl:46-50`: "Independent rate constant names for fitting (excludes Keq, E_total)".
- The user docs call it as `EnzymeRates.fitted_params(generator)` (`identify/tutorial.md:52`).
- `fit_rate_equation`'s documented default bounds use it (`fitting.jl:184-185`), and `loss!` reads it as a compile-time constant (`fitting.jl:121-129`).
- A test pins `parameters(m, Reduced) == (fitted_params(m)..., :Keq, :E_total)` (`test/test_allosteric_golden.jl:26-31`).

### `rate_equation_string(m, [mode])` (`rate_eq_derivation.jl:437-537`)
It returns a multi-line String that "spells out how `rate_equation` evaluates". The user docs (`textbooks.md:50-59`) describe four sections: the `(; …) = params` line, the `(; …) = concs` line, the optional `# Wegscheider constraints:` / `# Haldane constraints:` sections, and `v = E_total * (num) / (den)`.

Details the user docs do not mention:
1. **Section labels** are assigned by whether the right-hand side mentions `Keq` (Haldane) or not (Wegscheider) (`_equation_text`, 498-511). They do not follow the cycle class the docs define.
2. **Sorting:** each section is sorted by name, which is "load-bearing" for `eq_hash` (500-503).
3. **Substituted ties:** the Reduced string for an `EnzymeMechanism` solves without the Wegscheider rename. It marks single-symbol ties with `  (substituted into v)` (`ANNOTATION_SUBSTITUTED`, 10-16, 524-532).
4. **Pinned constants:** a constant pinned by a loop prints as `= 1.0` (`thermodynamic_constraints.md:39-43`).
5. **Invalid identifiers:** names that are not Julia identifiers "print bare, so such a string reads as text, not code" (456-458).
6. **v-line names:** the v-line uses raw names, dependents included, except that a single-symbol binding-K tie is renamed into it.

### `metabolites(m)` (`types.jl:1298-1339`)
- Returns substrates, then products, then regulators, each sorted by name and deduplicated. Allosteric mechanisms append the regulatory ligands.
- It stays `@generated` because `loss!` needs the tuple as a compile-time constant (1322-1326).
- Regulators that no step binds are dropped from the compiled type (`_sig_of`, 1042-1066; `developer.md:45-47`).

### `rescale_parameter_values(m, params; scale_k_to_kcat=1.0)` (`rate_eq_derivation.jl:730-746`)
- It multiplies every SS rate constant (`_ss_rate_constant_names`, 541-555) by `scale/kcat`. Ks, `Keq`, `E_total`, `L` and regulatory Ks are untouched.
- kcat is the maximum forward rate at saturating substrates, zero products and `E_total = 1`, taking the maximum over paths (577-592). For allosteric mechanisms it is the maximum over saturating patterns and regulator corners, per active site (627-638).
- The docs justify the rescaling by time units (`fitting/normalized_vs_absolute.md:59-81`).
- It accepts `AbstractEnzymeMechanism` only. `ROADMAP.md:90-97` proposes replacing `_kcat_forward`.

### Types and `show`
- `EnzymeReaction` docstring: `types.jl:332-374`, with a doctest that prints `EnzymeReaction: S ⇌ P`.
- `EnzymeMechanism{Sig}`: "exact layout is internal — users construct via `EnzymeMechanism(::Mechanism)`" (1102-1118).
- **Contradiction:** `AllostericEnzymeMechanism` is exported, but its docstring calls it an "Internal singleton…opaque fast-path handle". It sends users to `AllostericMechanism`, which is not exported (1142-1155).
- `show(::EnzymeMechanism)` prints steps such as `E + P <--> ES`. It uses the stored binding direction and concatenated form names, and shows no groups (1247-1265). The allosteric `show` prints groups in parentheses with their tags (1267-1285).
- `getting_started.md:39-41` explains that the docs check `isa` because "the Sig … is an unreadable string".

### Identify outputs (`identify/model_selection.md:124-154`)
- **`mechanism_type`:** the docs call it the "Julia type name of the compiled mechanism". The code calls it "the round-trippable parametric `EnzymeMechanism{Sig}` string" (`identify_rate_equation.jl:329-336, 464-467`).
- **`eq_hash`:** the docs promise "Two mechanisms with the same eq_hash compute the same rate function" (model_selection.md:144).
- **Parameter columns:** one column per fitted parameter name.
- **How `eq_hash` is computed:** it hashes the string minus `#` lines and minus `(substituted into v)` lines (`identify_rate_equation.jl:252-268`). The params line and the multi-symbol dependent lines are therefore part of the key.

## 2. Kinetic groups, thermodynamics and allosteric tags as documented

### Kinetic groups (`deriving/kinetic_groups.md`)
- **Meaning:** "a set of steps that share one kinetic parameter", a modeling choice for chemically identical steps (3-9). It is "distinct from the thermodynamic reduction … a decision you make: it ties together steps the thermodynamics leaves independent" (33-37).
  - **Imprecise:** an SS group shares two constants, forward and reverse (`types.jl:1437-1439`).
- **Syntax:** parentheses inside `steps:`; by default each step is its own group (11-20).
- **Rules** (22-31; enforced at `types.jl:703-777`):
  - all steps have the same `(consumed, released)`;
  - all steps share one RE/SS flag;
  - fused and plain bindings of one metabolite may share a group;
  - a Theorell–Chance step may not share a group with a binding;
  - conformational changes in one group run between the same two conformations (`_orient_tied_steps`, 538-574);
  - each reaction appears once.
- **Interaction with thermodynamics:** the docs say grouping turns Wegscheider relations into identities that "drop out" (111-116). In the code, grouping and the constraint solve happen in the same pass:
  - each group's cycle columns are merged into one before elimination, "equivalent to a kinetic-group equality constraint" (`thermo…jl:226-230, 352-355`);
  - a second solve, the Pass-2 rename (`rate_eq_derivation.jl:66-103`), collapses "split kinetic groups that Wegscheider ties back together" so they hash the same (85-87).
- **Naming representative:** a group's constants are named after one step, the argmin of (`_step_priority`, canonical key). The order is: free-enzyme RE binding (−1), free-enzyme SS metabolite step (0), other metabolite step (10), isomerization (20) (`thermo…jl:83-110`).
  - The user docs show only the result (`K_EA_to_E_A`, kinetic_groups.md:101-103).
  - Known issue (ROADMAP.md:10-24): pruning the inactive state can remove the representative and give a shared constant a second, fitted name.

### Thermodynamic constraints (`deriving/thermodynamic_constraints.md`)
- **Definitions** (15-37): a Haldane constraint comes from a cycle whose net change matches the reaction, and its dependent carries `Keq`. A Wegscheider constraint comes from a cycle with zero net change.
- **How they are found:** as the exact integer null space of the form incidence matrix, classified by net metabolite change. They are linearly independent by construction and computed in exact rationals (45-54; code at `thermo…jl:171-221`).
- **Infeasible mechanisms:** a contradiction `0 = c·log Keq` raises an error (56-59; `thermo…jl:539-541`).
- **Documented choice of dependent constants** (61-82):
  - free-enzyme binding Ks stay independent;
  - bindings on occupied forms are eliminated next;
  - isomerization rates are eliminated first;
  - within an SS step, the reverse rate is eliminated first;
  - for a release that carries chemistry and is stored as a binding (`E(S) <--> E + P`), the rebinding rate is eliminated first.
- **Determinism promise:** "two mechanisms with the same structure written in a different step order produce the identical reduced equation" (79-82).

Three behaviors the code adds beyond the docs:
- **(a) "Reverse" means stored-reverse.** For a plain SS product release, which is stored as a binding, the release constant is the one preferred as dependent.
- **(b) Group priority comes from the last step.** A grouped column takes the priority of the *last* step in canonical order, "not the highest such score" (`thermo…jl:89-93, 376-388`).
- **(c) Ties break by name.** Equal priorities fall back to ascending column order, which is canonical step order, which is sorted by rendered names (`thermo…jl:513-515, 537`; `types.jl:619-630, 663-666`).

**Consequence (inference):** renaming a conformation can change which constant is dependent. In the iso example (`iso.md:20-33`), two columns tie at priority 21 and `k_EprimeP_to_ES` wins because it comes first. Naming the second conformation `D…` instead would sort the D→E isomerization first and make `k_E_to_D` dependent. The rate function is the same, but the string and `eq_hash` differ. The docs promise invariance only under step order and direction, not under relabeling or regrouping. This is the class of problem goal 5 targets.

**What the user sees:** the labels follow §1's rule. In the MWC "collapse" example, the cross-state tie `K_I_ES_to_E_S = K_A_ES_to_E_S` is printed under `# Wegscheider constraints:`. The prose credits it to the Haldane relation (`mwc_allostery.md:180-184`); the code builds it as "I-row − A-row (the log Keq cancels)" (`rate_eq_derivation.jl:970`).

### Parameter naming (`re_vs_ss.md:18-44`, `developer.md:92-98`, `types.jl:1341-1411`)
- **Name formats:**
  - RE binding: `K_<bound>_to_<free…>`, a dissociation constant named in the release direction.
  - RE isomerization or other RE step: `K_<from>_to_<to>`, products over reactants in the stored direction, so not a Kd.
  - SS step: `k_<from>_to_<to>` plus its reverse.
- **Sides and forms:** a side is the form name plus its free metabolites, joined by `_` (`_side_label`, 679-685). A form name is the conformation plus the bound labels with no separator (`_species_name`, 102-116). Colliding names are rejected (`getting_started.md:43-50`), and adding a separator was decided against (ROADMAP.md:156-158).
- **Orientation:** the direction a step is written in does not matter. Tier 1 is substrate/product progression, Tier 2 is entry/exit context, and Tier 3 is name order (`re_vs_ss.md:35-44`; `types.jl:498-532`). Tier 3 makes some orientations depend on names.

### Allosteric MWC (`deriving/mwc_allostery.md`)
- **States and `L`:** A corresponds to R and I to T. `L` is the inactive/active ratio of free enzyme (9-18; `identify/tutorial.md:59`). `L` multiplies the inactive branch (`rate_eq_derivation.jl:1042-1045`).
- **Assumptions** (84-107): regulators bind at rapid equilibrium, and the enzyme changes conformation only while free, keeping it for a whole cycle. The alternative model differs by "a few percent" under `:NonequalAI` catalysis.
- **Formula** (109-116):
  - `num = (N_A/D_A)(Q_A/D_A)^(n−1)·W_A + L·(N_I/D_I)(Q_I/D_I)^(n−1)·W_I`
  - `den = (Q_A/D_A)^n·W_A + L·(Q_I/D_I)^n·W_I`
  - `W = Π(1 + lig/K)^m` over regulatory sites.
  - The rate is per catalytic site; catalytic multiplicity enters only through exponents (79-82). `E_total` is the protomer concentration (`textbooks.md:72`).
- **Normalization** (129-142): `D = 1` for a single RE segment. Otherwise the code divides by `D` when it is a single rate constant and cross-weights when it carries a metabolite. The code has three cases (`_mwc_state_polys`, `rate_eq_derivation.jl:1079-1115`).
- **Tags** (54-70):
  - Catalytic groups: `:EqualAI` gives one shared bare symbol; `:NonequalAI` gives `k_A_…` and `k_I_…`; `:OnlyA` is pruned from the inactive state and named with `A_`. `:OnlyI` is rejected for catalytic groups (`types.jl:918, 941-946`).
  - Regulators: `K_A_<lig>reg` and `K_I_<lig>reg`.
  - **(inference)** An `:EqualAI` regulator renders as `K_A_<lig>reg` in both states plus a printed mirror line `K_I_<lig>reg = K_A_<lig>reg` (999-1002, 1023-1026). That differs from the bare catalytic `:EqualAI` symbol the table describes.
- **Inactive state** (122-142):
  - The docs say `:OnlyA` groups are pruned, and that "the package sets N_I = 0".
  - The code also prunes every form disconnected from free E (`_inactive_groups`, 767-845). It derives `N_I` on the pruned graph and keeps it when it is nonzero ("live redundant-path", 1157-1166).
- **Thermodynamics across states** (149-249):
  - A lone `:NonequalAI` binding collapses to `:EqualAI`.
  - Two coupled bindings form a K-system with a fixed ratio.
  - An `:OnlyA` binding forces `:OnlyA` catalysis unless it is balanced. The constructor checks this (`thermo…jl:427-509`, called from `types.jl:963-965`).
  - A steady-state `:NonequalAI` binding constrains the affinity but leaves the speed free.
- **`1 ^ 2`** is documented as intended (144-147). It comes from `_reg_site_expr` emitting a literal `1`, which `_power_expr` then raises to the multiplicity (1028-1040, 1145-1146).
- **Existing design work:** ROADMAP.md:45-63 records a design that derives the catalytic equation once and takes limits as leading-order terms. It is on branch `mwc-solve-then-limit` (spec b067aa0/651a821, prototype 6f0011e), with `test/allosteric_ground_truth.jl` as the oracle.

### Regulators and competitive inhibitors
- **Docs' view** (`regulators.md:1-17`): a regulator's kind "is not a setting you choose but an emergent consequence of where it binds". Identification adds only dead-end inhibitors and allosteric regulators.
- **Four labels for one concept**, all producing `CompetitiveInhibitor` (`dsl.jl:228-231`):
  - `competitive_inhibitors:` and `dead_end_inhibitors:` in `@enzyme_reaction` (`dsl.jl:28-32, 79-81`);
  - `regulators:` in `@enzyme_mechanism` (249-252);
  - `catalytic_inhibitors:` in `@allosteric_mechanism` (568-569).
- **`::` means three different things:**
  - `X::Inh` is the inhibitor copy of a name that is also a substrate or product (253-255, 487-517).
  - `A::Activator` / `A::Inhibitor` in `@enzyme_reaction` pins the enumeration type (33-40; `identify/tutorial.md:174-180`).
  - `A::OnlyA` and step `:: EqualAI` in `@allosteric_mechanism` set the allosteric state (570-577).
- **The `inh` label:** every competitive inhibitor carries an `inh` label (`types.jl:51-57`). `regulators.md:41-44` explains `K_EIinh_to_E_Iinh` as "because it occupies a regulator binding site". In the same page, the essential activator renders as `K_ERinh_to_E_Rinh` and `k_ERinhS_to_EPRinh` (build output), while the prose says "Nothing in the DSL marks R as an activator" (70-71).

### Ping-pong, iso, RE vs SS, Cha/King–Altman
- **Ping-pong** (`ping_pong.md`):
  - The residual syntax is `E(; residual = A - P)`, listing atoms gained (+) and committed to release (−) (11-17).
  - Forms are named `E_res_+A_-P` (38-40). Constants such as `k_EA_to_E_res_+A_-P_P` are not identifiers and need `var"…"` to pass in a NamedTuple.
  - The docstring example `k_EA_to_EP_res_+A_-P` (`rate_eq_derivation.jl:456-458`) is not a name this page produces.
  - ROADMAP.md:98-104 plans `_res_pA_mP`; that change moves `eq_hash` and requires banning `_` in metabolite names.
  - Documented signature: numerator `A·B − P·Q` and a denominator with "no constant term" (33-37).
- **Iso** (`iso.md`): signature is an S·P cross term (35-42). Iso mechanisms are not enumerated (`roadmap.md:18-24`).
- **RE vs SS** (`re_vs_ss.md`): `⇌` costs 1 constant and `<-->` costs 2 (3-44). A chain with RE flanks gives every law that SS flanks give, so the flank K behaves as an apparent Km (126-140).
- **Cha/King–Altman** (`cha_king_altman.md`):
  - About 100 ns per call (3-13).
  - Each RE segment is rooted at its least-bound form, which yields the `1 + S/K…` denominator (17-27).
  - The output is "division-free in the concentrations" and finite at zero (36-48; `sym_poly…jl:45-74`).
  - No dimensionless form is documented. The iso and all-SS examples print multi-term constant parts (e.g. `k_E_to_Eprime + k_Eprime_to_E + …`), and ping-pong prints none.

## 3. Textbook bridge today

**User-facing, prose only:**
- `textbooks.md` (titled "Derivation tutorial") promises readers will "recognize the package's parameter names as the textbook constants you already know" (3-6).
- Its table (61-72) says `K_ES_to_E_S` is "analogous to Km under the rapid-equilibrium approximation, but… not the same as Km for more complex mechanisms", and lists units for the rate constants, `Keq` and `E_total`.
- `re_vs_ss.md:126-140` covers apparent constants.
- `enumeration_engine.md:43-47` cites the Cleland forms (`V·A·B/(Kia·Kb + …)`) only to motivate the merged and Theorell–Chance seeds.

**No code maps a textbook equation onto a derived one.** Nothing translates Km, Ki, Kia or Vmax, and nothing rewrites a derived equation into apparent-constant form. The only textbook quantity computed is kcat (`_kcat_forward`, internal).

**Test side (the working bridge):**
- `test/mechanism_definitions_for_test_enzyme_derivation.jl:185-650` transcribes Segel's equations (IX-8, -45, -60, -87, -122, -140, -195, -228, -261) in positional microscopic constants `k1f`, `k1r`, …, not Km or Ki.
- The mapping runs in four steps:
  1. Evaluate the package's dependent expressions numerically (`_testhelper_compute_all_params`, `test_rate_eq_derivation.jl:479-489`).
  2. Re-key structural names to positional `K{i}`, `k{i}f`, `k{i}r` (plus `_T` variants and `K_<lig>_reg{site}`) by matching each canonical step to its as-written step (205-404).
  3. Swap f/r for product releases stored as bindings (406-451).
  4. Compare at random points with `rtol 1e-10` (818-845).
- The as-written order comes from `@_testhelper_enzyme_mechanism_src` (fixtures 59-73). That macro relies on `_parse_mechanism_body` returning source-order groups (`dsl.jl:653-663`), a DSL feature that exists only for tests.
- **Other oracles:**
  - QSSA/ODE (`test_rate_eq_derivation.jl:795-914`);
  - mass action (158-201);
  - the allosteric ground truth;
  - `eval(Meta.parse(string))` equivalence (679-683, 916-951);
  - exact factored num/den strings (e.g. mechanism_definitions:1094-1098);
  - the byte-identical `test/reference/allosteric_golden_reference.txt`.

The bridge proves that two **functions** agree numerically. Nothing relates their **forms** or named apparent constants.

## 4. Doc examples affected by likely changes

Key: N = names, D = choice of dependent, Z = normalization, L = string layout and labels, P = `1 ^ 2`, T = type/`show` output. A star (★) marks a pinned doctest, or an `@example` with hard-coded names that breaks the build if names change.

| Location | Content | Affected by |
|---|---|---|
| textbooks.md:41-48 ★; rate_eq_derivation.jl:476-481 ★ | uni-uni string doctest | N D Z L |
| types.jl:1318-1319 ★, 371-372 ★ | `(:S, :P)`; `EnzymeReaction: S ⇌ P` | order/show |
| textbooks.md:50-59, 61-72, 91-96, 104-111 | section list, name table, params, Full | N D L |
| getting_started.md:59-69, 155 | params and equation prints | N D Z L |
| getting_started.md:84-95 ★; README.md:42-45 (not executed, goes stale silently) | hard-coded `K_ES_to_E_S`… | N D |
| thermodynamic_constraints.md:23-25, 35-37, 39-43 | static Haldane, Wegscheider and `= 1.0` lines | N D L |
| thermodynamic_constraints.md:90-118 | random-order string and params | N D Z L |
| kinetic_groups.md:51-116 | both strings and params; prose lists 4 Wegscheider-fixed names (75-80) and denominator terms (109-113) | N D Z L |
| re_vs_ss.md:18-44, 56-112 | naming rules; params ("Six"/"Nine" include `Keq`/`E_total`, while kinetic_groups.md:75's "nine" excludes them); strings | N D Z L |
| cha_king_altman.md:21-27, 36-48 | `1 + S/K` claim; division-free promise | Z |
| ping_pong.md:18-40 | string; "no constant term"; `E_res_+A_-P` | N Z |
| iso.md:20-42 | string; cross-term snippet; names | N D Z |
| regulators.md:25-71 | two strings; `K_EIinh_to_E_Iinh` | N Z |
| mwc_allostery.md:26-70 | string (`1 ^ 2`, I-branch terms ordered differently: `S * k_I_ES_to_EP…`, from 1129-1131), params, tag table | N D Z L P |
| mwc_allostery.md:109-147 | formula, D normalization, `1 ^ 2` prose | Z P |
| mwc_allostery.md:162-249 | collapse/coupled strings; stated relations at 180, 207, 229-231; ss params | N D L |
| fitting/tutorial.md:58-60 ★, 71, 107 | hard-coded names; `FittingProblem` shows the full Sig; result params | N D T |
| fitting/normalized_vs_absolute.md:95-106 ★ | hard-coded names; rescaled output | N D |
| identify/tutorial.md:52-60, 79-80 ★, 190, 197, 200-205, 212 | fitted params; true params; `best` show; string with `1 ^ 2`; dependent `k_A_EP_to_ES`; `cv_results` (Sig-string `mechanism_type`, parameter columns, `eq_hash`) | N D Z L P T |
| identify/model_selection.md:139-146 | column definitions | T L |
| dsl.jl:268-275 | `:ES`, `:EPS`, `:Estar_res_+A_-P` | N |

**Broken docstring example (inference):** the `@enzyme_mechanism` example (`dsl.jl:239-244`) is all-RE with a closed catalytic cycle. It would fail in `_compute_alpha` (`rate_eq_derivation.jl:183-195`).

**Test pins (not docs) that would also change:**
- the golden file;
- the factored strings;
- `test_pivot_priority_regression.jl:16`, which `Meta.parse`s an `EnzymeMechanism{Sig}` literal.

## 5. Documented rationale and promises

### Rationale on record
- **Two representations on a performance boundary** (`developer.md:7-41`):
  - Singleton types move derivation to compile time, so `rate_equation` is a flat expression, allocation-free and under 120 ns. A test enforces this (`test_rate_eq_derivation.jl:874-888`; CLAUDE.md).
  - The cost: every new type pays a full derivation, and large mechanisms can be slow, exhaust memory, or overflow the stack (36-41).
  - Guards: `MAX_RATE_EQUATION_TERMS = 5000` (`sym_poly…jl:4-9`; `rate_eq_derivation.jl:303-347`) and the default `eq_complexity_filter = 337`.
- **Concrete value types for enumeration:** building and deduplicating candidates compiles nothing, and only fitted candidates are lifted (`developer.md:122-130`).
- **Keeping per-type codegen small:**
  - Only the `@generated` methods do per-type work. Helpers are `@nospecialize`, the lifts read `Sig` at run time, and `FittingProblem`/`loss!` specialize on purpose (`developer.md:49-63`; `types.jl:1023-1027`).
  - The lift is deliberately not the inverse of construction: it drops `shared_catalytic_site`, `reg_type` and unbound regulators "so that equivalent mechanisms share one compile" (1134-1138).
  - `AllostericEnzymeMechanism` has three type parameters because a DataType cannot sit in a value tuple (1151-1154).
  - `metabolites` stays `@generated` for type stability (1322-1326).
  - `_nest_binary` builds balanced trees to avoid varargs boxing (`sym_poly…jl:136-151`).
  - Compile budgets are tested (`test/test_compile_budget.jl:13-14`).
- **Canonical Step Form is load-bearing** (CLAUDE.md; `developer.md:96-121`):
  - Dependents are chosen by step order, and deduplication by `==`/`hash` works only because construction is canonical.
  - Canonical keys use rendered names, not `hash`, because the default struct hash mixes in `objectid` and is not stable across precompile sessions (`types.jl:619-626`).
  - The segment-root rule and the ω_e numerator weights make the derivation independent of step order and of the direction a step is written in (`rate_eq_derivation.jl:116-121, 276-282`).
- **String determinism for deduplication:**
  - independents sorted by name (`thermo…jl:548-552`);
  - dependent lines sorted (`rate_eq_derivation.jl:500-503`);
  - factors sorted (`sym_poly…jl:216-218`);
  - `eq_hash` is exact-render identity only (`identify_rate_equation.jl:259-260`).
- **Naming chokepoint:** all parameter names go through `name(p, m)`, enforced by an AST-walker test (CLAUDE.md; `types.jl:1341-1346`).

### Promises to keep, or retire explicitly
**Keep:**
- 0 allocations and under 120 ns per `rate_equation` call; `loss!` allocation-free.
- `Keq` always supplied by the user.
- `E_total` is the active-site concentration, and the rate is reported per site.
- Step order and direction never change the type, names, dependents or equation (`developer.md:18-21`; `thermodynamic_constraints.md:79-82`; `re_vs_ss.md:38-43`).
- One dependent per independent cycle, expressed only through `Keq` and independents (`rate_eq_derivation.jl:364-369`); infeasible mechanisms raise an error.
- Every fitted equation satisfies detailed balance (`thermodynamic_constraints.md:116-118`).
- Equations are division-free in concentrations and finite at zero.
- Equal `eq_hash` implies the same rate function.
- Free-enzyme binding Ks are preferred as independents.
- `:OnlyI` catalysis is rejected, and inconsistent `:OnlyA` combinations are rejected.
- The docs build passes doctests and `@example` blocks.

**Currently documented, candidates to retire or revise:**
- `Full`/`Reduced` modes (goal 12).
- Section labels decided by whether `Keq` appears, and the `(substituted into v)` annotation (goal 12; they also contradict the MWC prose).
- The `1 ^ 2` prose (goal 14).
- The "Sig is unreadable, use `isa`" note, and `mechanism_type` as a round-trippable Sig string (goals 13 and 18).
- `regulators:` meaning competitive inhibitor, the four labels, the three meanings of `::`, and the `inh` label in names (goals 15 and 16).
- `_res_+A_-P` names and the "print bare" docstring (goal 17).
- The "N_I = 0" and pruning description (goal 7).
- Today's normalization, which is only "RE root weight 1 plus a concentration shift" (goal 6).
- `Mechanism` holding a `reaction` field (`developer.md:69-71`; goal 8).
- `fitted_params` being unexported but used in the docs.

**Missing promise:** the docs promise invariance under step order and direction only. Goal 5's stronger invariance, under relabeling and regrouping, does not hold today, because ties break by name and a group takes its last step's priority (§2).