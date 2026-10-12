# DSL slice: inventory and audit of `src/dsl.jl`

## Headline findings

- **Size.** dsl.jl is 752 lines: 421 code, 287 docstring, 13 comment and 32 blank (my count). About 99 of the code lines are error-message text spread over 54 `error(` calls.
- **Denis's ">1200 LOC" is the file before v0.9.0.** It was 1256 lines at e59ebd6 and stayed between 1177 and 1259 lines from 80367d0 to e59ebd6. v0.9.0 (277a785, merged today) cut it to 752: it merged the two mechanism parsers (a1fe2b3) and emitted Exprs while parsing (9c3f270). Over the same change test/test_dsl.jl grew from 772 to 1163 lines. DSL source plus tests is now 1915 lines.
- **No legacy syntax is left.** Opaque bound-form names are rejected (dsl.jl:476-485), and the legacy emitter was deleted in a203856. What remains:
  - about 8 duplicate surface forms;
  - one inert feature (per-regulator multiplicities);
  - three label names for one concept;
  - an allosteric-regulator declaration split across two places;
  - heavy prose.
- **Removing every duplicate saves about 110–140 code lines (25–33%).** The core (steps, species, metabolites, residuals, groups and tags) is about 150 code lines and stays. Reaching half needs three more things:
  - Denis's goal 8 (no `EnzymeReaction` inside `Mechanism`);
  - one grammar documented once;
  - trimmed error prose.

## 1. Syntax inventory

Counts are lines containing the form (grep). "docs" means docs/src + README + docs/hpc_results. Outside dsl.jl, src uses the DSL only in three docstring examples (types.jl:363, types.jl:1309, rate_eq_derivation.jl:466), all basic forms.

**`@enzyme_reaction`** lowers to `EnzymeReaction(ReactantAtoms[...], RegulatorMults[...], mults; shared_catalytic_site)` (dsl.jl:203-217).

| # | Form | Example | Lowers to | Uses: test / docs | Verdict |
|---|---|---|---|---|---|
| R1 | reactant + atoms | `S[C6H12O6]` | `ReactantAtoms(Substrate(:S), [:C=>6,…])` (dsl.jl:204-207) | 164 / 4+1+3 (README.md:58, pgk.jl:50) | core |
| R2 | multi-arg bracket | `B[N, P]` | same as `B[NP]` (dsl.jl:119-120) | 2, only test_dsl.jl:469-470 / 0 | duplicate, test-only |
| R3 | `competitive_inhibitors:` | `I` | `RegulatorMults(CompetitiveInhibitor(:I), cat_mults, :unspecified)` | 17 / 2 (pgk.jl:53, ldh.jl:57) | core (enumeration input) |
| R4 | `dead_end_inhibitors:` | `I` | same branch as R3 (dsl.jl:79-81) | 36 (31 in test_mechanism_enumeration.jl, e.g. :219) / 0 | duplicate, test-only |
| R5 | `allosteric_regulators:` ± `::Activator`/`::Inhibitor` | `ATP::Inhibitor` | `RegulatorMults(AllostericRegulator, …, :inhibitor)`, read by `_respects_reg_type` (mechanism_enumeration.jl:1866) | ~43 / 3 (identify/tutorial.md:34,179; pfkp.jl:56) | core |
| R6 | per-regulator multiplicities | `A(1, 2)`, `X(1,2)::Inhibitor` | `RegulatorMults.allowed_multiplicities`; its only reader is `_sig_of` (types.jl:1062) | 8 (5 in test_dsl) / 0, though advertised in dsl.jl:9 | inert; ROADMAP.md:81-89 plans to make it work |
| R7 | `allowed_catalytic_multiplicities:` | `(1, 2, 4)` | the field (read at mechanism_enumeration.jl:1581) | 11 (8 in test_dsl) / 0 | duplicate of R8, test-only |
| R8 | `oligomeric_state:` | `4`, or `(4,)` | same field | 67 / 4 | core; the `(4,)` variant is test-only (test_dsl.jl:778) |
| R9 | `shared_catalytic_site:` | `(ATP, ADP)` | kwarg normalized at types.jl:410-429, read at mechanism_enumeration.jl:323 | 10 / 1 (pgk.jl:54) | core |

**Mechanism macros.** Both share `_parse_mechanism_body(block, allosteric::Bool)` (dsl.jl:664).

| # | Form | Example | Lowers to | Uses: test / docs | Verdict |
|---|---|---|---|---|---|
| M1 | bare reactants | `substrates: A, B` | names plus a placeholder-atom `EnzymeReaction` (dsl.jl:228-231) | 652 / 22 | core |
| M2 | `regulators:` (plain macro) | `I` | role `CompetitiveInhibitor`; the reaction gets `RegulatorMults(…, [1])` | 161 / 3 (deriving/regulators.md:30,56) | core concept |
| M3 | `catalytic_inhibitors:` (allosteric macro) | `I` | identical to M2 (dsl.jl:668-669) | 36 (e.g. test_types.jl:428) / 0 | duplicate label, test-only |
| M4 | `steps: begin…end` | | `EnzymeMechanism(Mechanism(rxn, groups))` (dsl.jl:748-749) | 483 / 16 | core |
| M5 | `catalytic_steps:` | | `AllostericEnzymeMechanism(AllostericMechanism(…))` (dsl.jl:742-745) | 223 / 5 (mwc_allostery.md:33) | M4 under another name |
| M6 | `catalytic_multiplicity:` | `2` | field; defaults to 1 (dsl.jl:737) | 194 / 5 | core |
| M7 | `allosteric_regulators:` (allosteric macro, tag required) | `A::OnlyA` | per-ligand state for a `RegulatorySite`; never put on the reaction (mechanism_enumeration.jl:1641-1643) | 67 / 5 (mwc_allostery.md:32) | core content, declared apart from its site |
| M8 | `regulatory_site(multiplicity = N): begin ligands: … end` | | `RegulatorySite(AllostericRegulator[…], N, states)`; a ligand with no site gets a one-ligand site at the catalytic multiplicity (dsl.jl:643) | 38 blocks / 5 | see below |

M8 in detail:
- Needed for multi-ligand sites (15 test blocks; the merge move emits them) and for site multiplicities that differ from the catalytic one (14 test blocks).
- All 5 docs blocks (mwc_allostery.md:38,41,173,200,244) and 7 test blocks restate the default, so they duplicate it.
- 38 allosteric fixtures rely on the implicit default instead.

**Steps**

| # | Form | Example | Lowers to | Uses: test / docs | Verdict |
|---|---|---|---|---|---|
| S1 | arrows | `E + S ⇌ E(S)`, `E(S) <--> E(P)` | `Step(from, to, consumed, released, op == :⇌)` (dsl.jl:433-434); no other operator is accepted (dsl.jl:409-410) | ⇌ 2427, `<-->` 1297 / 82, 52 | core |
| S2 | enzyme form | `E`, `E(A, B)` | `Species(Metabolite[…], :E, Residual())` | everywhere | core |
| S3 | other conformations | `Eprime(P)`, `Estar`, `F` | `Species(…, :Eprime, …)`; shape regex at dsl.jl:484-485 | ~128 lines with Estar; Segel ping-pong written with F and G (mechanism_definitions…:884-899) / iso.md:27-29 | core for hand-written mechanisms; the enumerator emits only `:E` (mechanism_enumeration.jl:126,180) |
| S4 | residual | `E(B; residual = A - P)` | `Residual(Substrate[A], Product[P])` (dsl.jl:456-463, 524-543) | 239, all a single `X - Y` / ping_pong.md:15-27 | core; multi-term residuals are buildable by ter ping-pong seeds (mechanism_enumeration.jl:122-126) but used nowhere |
| S5 | any number of free terms per side | `E(A) + B <--> E(Q) + P`; `E + A + B ⇌ E(A, B)` | `consumed` / `released` | Theorell–Chance common (test_types.jl:210); multi-metabolite sides 3 (test_rate_eq_derivation.jl:1207,1307,1309) | core; one code path |
| S6 | untagged group | `(E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))` | one inner `Vector{Step}` (dsl.jl:342-347) | 456 / 5 (kinetic_groups.md:17) | core |
| S7 | tagged group | `(…, …) :: NonequalAI` | group plus a `cat_allo_states` entry (dsl.jl:332-340) | 64 (test_rate_eq_derivation.jl:2154) / 0 | core |
| S8 | tagged step | `E(S) <--> E(P) :: EqualAI` | one-step group; the tag is peeled off the right side's last term (dsl.jl:384-392) | 826 / 15 | core |
| S9 | parenthesized tagged step | `(E(P) ⇌ E + P) :: EqualAI` | same as S8 (dsl.jl:332-337); test_dsl.jl:421 pins them equal | 2 (test_dsl.jl:392) / dsl.jl:555 docstring | duplicate |
| S10 | role copy | `E + P::Inh ⇌ E(P::Inh)` | `CompetitiveInhibitor(:P)`; requires P in M2/M3 (dsl.jl:508-512) | 695 tokens, 593 of them on reactant names (A, Q, S, P, B) in enumeration fixtures / developer.md:85 (prose) | core |
| S11 | `I` vs `I::Inh` for an inhibitor-only name | `E(I)` ≡ `E(I::Inh)` | both give `CompetitiveInhibitor(:I)` (dsl.jl:512) | `I::Inh` 26× in enumeration tests | duplicate |

**Other duplicates and quirks**
- `::` carries four meanings:
  - step or group state (`:: EqualAI`);
  - regulator state (`A::OnlyA`);
  - role copy (`X::Inh`);
  - reaction regulator type (`A::Activator`).
- `::Inh` and `::Inhibitor` are near-homonyms with unrelated meanings.
- List labels can be repeated and their values append (dsl.jl:680-684), but `steps:` and `catalytic_multiplicity:` cannot be repeated (dsl.jl:686-702).

## 2. Duplicated or dead code paths

1. **Label aliases (R4, R7/R8).** dsl.jl:79-81, 84-91 and 189-192, about 10 code lines. test_dsl.jl:752-797 exists only because two labels set one value.
2. **Inert per-regulator multiplicities (R6).** dsl.jl:106-107 and 176-183 (about 10 code lines), plus the `RegulatorMults` field (types.jl:311-329) and its Sig encoding and decoding (types.jl:1062, 1085).
   - Inference: a fixture's inhibitor carries `[1]` (dsl.jl:231), while an enumerated one carries the user's catalytic multiplicities (dsl.jl:107). The same mechanism therefore encodes two different Sigs.
3. **Mode switch between the two mechanism grammars (M2/M3, M4/M5).** dsl.jl:668-669, 705-709 and the parametrized message at 508-511, about 8 code lines. The larger cost is two 46-line macro docstrings that restate one grammar (dsl.jl:233-278, 545-590), plus the per-macro error testsets.
4. **Split allosteric declaration (M7/M8).** Tag list, site blocks and default sites are reconciled across four places, about 56 code lines:
   - dsl.jl:600-624 (24 lines);
   - dsl.jl:633-651 (19);
   - dsl.jl:689-695 (7);
   - dsl.jl:729-735 (6).

   One line per site would need about 15.
5. **One label, two parsers.** `allosteric_regulators:` goes through `_parse_regulator_entries` in the reaction (dsl.jl:163-185) and through dsl.jl:689-695 in the allosteric mechanism.
6. **S9.** A `:call` head is accepted inside `::` (dsl.jl:332-337), about 2 lines.
7. **`_peel_step_tag!`** (dsl.jl:373-392; 9 code, 11 doc).
   - It exists because `::` binds tighter than `+` and `⇌`, so `s :: Tag` attaches to the right side's last term.
   - It must exempt `::Inh` (dsl.jl:389), which needs a 76-line testset (test_dsl.jl:1087-1162).
   - Requiring parentheses on every tagged step would delete it, but would rewrite 826 fixture lines.
8. **Placeholder reaction** (dsl.jl:219-231). It exists only because `Mechanism` stores an `EnzymeReaction` (types.jl:902-910). Knock-on effects:
   - 95 `_testhelper_on_reaction*` calls rebuild fixtures on the real reaction (test_mechanism_enumeration.jl:152-161);
   - fixtures cannot run `_assert_atom_conserving` (mechanism_enumeration.jl:76-84);
   - fixture reactions omit allosteric regulators (mechanism_enumeration.jl:1641-1643).
9. **Return value kept for tests.** `_parse_mechanism_body` returns `(mech, groups, sites)` exprs (dsl.jl:653-662, 751) only for two test macros (mechanism_definitions_for_test_enzyme_derivation.jl:59-73). Those macros have 42 call sites feeding source-order oracles.
10. **No DSL printer.** The only code that renders a form back to DSL text is a local `text(sp)` in an error path (types.jl:712-719).
    - `show` prints rendered names (`EPinh`, `E_res_+A_-P`; types.jl:1250-1285) that the macros reject.
    - The identify CSVs store `string(typeof(em))` (identify_rate_equation.jl:466).
    - A printer plus a round-trip test would serve goal 13 and could replace many single-form pins.

## 3. Inhibitor and regulator semantics today, and the `::Inh`-as-suffix proposal

**How a type is chosen**
- **`@enzyme_reaction`: by label only.** `competitive_inhibitors:` and `dead_end_inhibitors:` give `CompetitiveInhibitor`; `allosteric_regulators:` gives `AllostericRegulator` with a reg_type (dsl.jl:79-83). One name may hold both kinds and also be a reactant (types.jl:431-435; test_dsl.jl:562-585).
- **Mechanism macros: through `role_of`**, built in the order AllostericRegulator, CompetitiveInhibitor, Substrate, Product, with later entries winning (dsl.jl:717-724):
  - a bare name binds by its reactant role first, then as an inhibitor;
  - an allosteric-only name is rejected in steps (dsl.jl:512-515);
  - `X::Inh` forces `CompetitiveInhibitor` and requires X in `regulators:`/`catalytic_inhibitors:` (dsl.jl:508-511);
  - inhibitors enter the synthesized reaction at multiplicity `[1]` (dsl.jl:228-231);
  - allosteric ligands go only into `RegulatorySite` (dsl.jl:633-651).
- **Naming.** `_met_label` appends `inh` to every `CompetitiveInhibitor` (types.jl:56-57).
  - Even a pure activator gets it: the essential activator renders as `K_ERinh_to_E_Rinh` and `k_ERinhS_to_EPRinh` (mechanism_definitions…:1231-1233), and regulators.md:40-43 documents this.
  - Site ligands are named `K_<state>_<name>reg` (types.jl:1410-1411).

**Why `::Inh` exists.** The dead-end move binds `CompetitiveInhibitor(name)` for every declared competitive inhibitor (mechanism_enumeration.jl:1463-1525, line 1507). The PGK and LDH HPC runs declare only reactants as competitive inhibitors (pgk.jl:53, ldh.jl:57). E(P), the product at its catalytic site, and E(P::Inh), P at a dead-end site, are different complexes. The role is kept apart in seven places:
- `Species` identity (types.jl:83-87);
- the sort key (types.jl:49);
- names (types.jl:56-57);
- the Step role check, which rejects `E + P::Inh → E(P)` (types.jl:226-236);
- Canonical Step Form, which counts a copy on neither side (types.jl:514-518, 606-612);
- the form-name collision guard (types.jl:709-728);
- the Sig type tags (types.jl:1035, 1069-1071).

**What the derivation needs.** It never mentions "competitive": grep finds no `CompetitiveInhibitor` in rate_eq_derivation.jl, thermodynamic_constr_for_rate_eq_derivation.jl or sym_poly_for_rate_eq_derivation.jl. It needs three things:
- **The metabolite name**, for the concentration. A copy and its reactant share it:
  - stoichiometry keyed by name (thermodynamic_constr:195-205);
  - numerator by first-substrate name (rate_eq_derivation.jl:283-290);
  - RE weights by name (rate_eq_derivation.jl:137-143);
  - test_dsl.jl:899-902 pins `(; S, P) = concs`.
- **Distinct identity of copy-bound forms.** Today this comes from rendered names: forms are keyed by `name(sp)` (thermodynamic_constr:176-178) and in `_reachable_from_free` (rate_eq_derivation.jl:779-795). So the `inh` suffix is load-bearing for correctness, not cosmetic, while forms are keyed by name.
- **One role-sensitive tie-break**, `bound_metabolite(s) isa Product` in the pivot priority (thermodynamic_constr:385).

The kind of inhibition is emergent (regulators.md:3-9).

**What enumeration needs.** It reads the role on each bound entry:
- productive forms (mechanism_enumeration.jl:1234);
- sites by role (:1395-1404, :1494);
- the inhibiting-names map (:1481);
- copy redundancy (:1358-1362);
- RE→SS eligibility (:450-451);
- the catalytic-site filter (:926-927);
- seed bookkeeping (:2029-2035; identify_rate_equation.jl:661-670).

**Assessment of the proposal.** At the reaction level, "competitive inhibitor" is already enumeration-only: `competitive_inhibitors:` is read only by the dead-end move (mechanism_enumeration.jl:1470-1472), the seed requirement and `_sig_of`. At the mechanism level, the four-subtype hierarchy can become `(name, label)` on bound entries:
- Substrate and Product come from the mechanism's lists (goal 8);
- position separates catalytic-site regulators (in steps) from site ligands.

This works if the label is data in identity, not only a suffix. What breaks or needs care:
1. **Dual-role names.** The label must take part in `==`, `hash`, the Sig and the sort key, or E(P) and E(P::Inh) merge. Declaring `regulators: P::Inh` once is not enough, because each occurrence in steps must still say which role it means. The per-occurrence `P::Inh` stays.
2. **Canonical Step Form is load-bearing.** If a copy were stored as a labelled Product, Tier 2 would mark the forms where it binds as product-entry forms. Inference: that could flip isomerization directions and so move dependent parameters, `fitted_params` names and `eq_hash`. The label rule must reproduce the current "neither side" rule exactly.
3. **Renaming.** A user-chosen suffix (`R::Act`, or none) renames every regulator constant:
   - golden strings (mechanism_definitions…:1099-1347);
   - regulators.md;
   - HPC CSV columns.

   Dropping the suffix is safe only for names no reactant shares, and `_assert_distinct_form_names` would catch collisions.
4. **Tag peeling.** dsl.jl:389 exempts only `Inh`. With open labels, `E(R::Act) ⇌ E + R::Act` reads like a state tag unless the parser checks the closed set {OnlyA, OnlyI, EqualAI, NonequalAI} or the declared labels. A misspelled state would then surface as an unknown label.
5. **Vocabulary clash.** `::Activator`/`::Inhibitor` already pin allosteric states in `@enzyme_reaction` (dsl.jl:166-174; mechanism_enumeration.jl:1850-1880). `::Act` and `::Inh` would be confusing near-homonyms.
6. **Allosteric regulator that also binds the catalytic site** (test_dsl.jl:1004-1021): position plus the `reg` suffix already separate the two roles. A label is needed only if the name is also a reactant.
7. **Code changes.** About 12 enumeration sites and the Sig encoder change from `isa CompetitiveInhibitor`/`isa Regulator` to label checks. In the derivation only thermodynamic_constr:385 changes.

## 4. What test_dsl.jl pins

The file has 34 testsets, 59 `@test_throws` (about 46 pin message text, 10 accept any `Exception`, 2 match an exact `ErrorException`) and 77 positive `@test`s.

- **Accepted forms:**
  - decomposed species, groups and residuals (4-136);
  - S9 equals S8 (380-434);
  - reaction labels and kinds (436-631, 728-835);
  - dual roles and `::Inh` (882-1055);
  - writing a binding as a release (1087-1150).
- **DSL errors:**
  - atoms in a mechanism (113);
  - allosteric-only labels in the plain macro (124, 293-314);
  - repeated or missing labels (263-355);
  - `regulatory_site` arguments (357-378);
  - formulas and multiplicities (487-536);
  - `regulators:` in a reaction (632-641);
  - opaque names (643-678);
  - undeclared names reported before the opaque-name check (680-703);
  - macro self-naming (705-726);
  - multiplicities set once (752-797);
  - malformed shared-site pairs (815-834);
  - enzyme-form count per side (852-864, 1057-1067);
  - non-step calls (866-880);
  - regulator misuse in steps (954-1055);
  - residual roles (1069-1085);
  - tags that are not Symbols (1153).
- **Constructor errors reached through the macros:** 152, 219, 232, 244, 451-465.
- **Weak pins** that pass on any error: 140, 165, 178, 580, 614-629, 816-834.

Under the grammar in section 5:
- **Deleted** (about 150–185 lines):
  - 380-435 (S9);
  - 752-797 (R7/R8);
  - 538-561 (R4; its assertions fold into 436);
  - 293-315 and 705-727, if the two mechanism macros merge (otherwise they shrink).
- **Changed (about 12):** 4-137 (labels), 138-262 and 357-379 (site syntax), 316-356, 467-486 (R2), 487-537 (R6), 562-631 (R6/R4), 728-751, 954-1055 (label in messages).
- **Unchanged (about 15).** The rest, if bare step tags stay.

Mechanical fixture renames outside test_dsl.jl: about 700 lines in total.
- `catalytic_steps:` 223 lines;
- `catalytic_inhibitors:` 36 and `dead_end_inhibitors:` 36;
- `allowed_catalytic_multiplicities:` 11;
- about 74 allosteric fixtures with regulators (38 of them rely on default sites).

## 5. Minimal grammar sketch, LOC estimate, risks

```
@enzyme_reaction begin                       # enumeration input
    substrates: A[C3H4O3], B[...]            # Name[Formula], one bracket arg
    products:   P[...], Q[...]
    competitive_inhibitors: B, P, I          # dead-end move candidates
    allosteric_regulators: X, Y::Activator, Z::Inhibitor
    oligomeric_state: 4 | (1, 2, 4)
    shared_catalytic_site: (A, P), (B, Q)
end
@enzyme_mechanism begin                      # allosteric iff catalytic_multiplicity/sites present
    substrates: A, B
    products:   P, Q
    regulators: I, P                         # catalytic-site regulators
    catalytic_multiplicity: 2
    steps: begin
        E + A ⇌ E(A)                                   :: EqualAI
        (E(A) + B ⇌ E(A, B), E(B) + A ⇌ E(A, B))       :: NonequalAI
        E(A, B) <--> E(Q; residual = A - P) + P        :: EqualAI
        E + P::Inh ⇌ E(P::Inh)                         :: EqualAI
    end
    regulatory_site(2): X::OnlyA, Y::OnlyI           # one line per site, no default
end
```

Grammar rules:
- `side := term {+ term}` with exactly one enzyme form per side.
- `form := Conf | Conf(Met, …; residual = ±Name …)`, where `Conf` matches dsl.jl:484-485.
- `Met := Name | Name::Inh`.
- A group is a step or a parenthesized tuple of steps, each with an optional `:: State`.

Decisions for Denis:
- one macro, or two macros sharing these labels;
- the name of the mechanism regulator label;
- keep R6 for the ROADMAP item (about 10 lines) or drop it;
- explicit sites only, or default sites only;
- whether goal 8 lands first.

**LOC (code lines)**

| Area | Today | Minimal |
|---|---|---|
| reaction parsing and emission | 120 | ~80 |
| label loop, roles, emission, macros | 97 | ~60 |
| steps block and tag peel | 51 | ~35 |
| step, species, metabolite, residual, shape | 106 | ~92 |
| allosteric sites | 43 | ~14 |
| placeholder reaction | 4 | 0 with goal 8 |
| **total code** | **421** | **~280** |

With docstrings at about 120 (grammar documented once) and about 28 blank lines, the total is about 430, 43% below today. Trimming the error prose (about 99 lines down to about 55) brings it to about 380, roughly half. Removing syntax alone does not get there.

A further option to prototype: build the mechanism value at macro-expansion time instead of emitting constructor Exprs. That would remove the Expr boilerplate. Inference: constructor errors would then fire at expansion and break the `@test_throws … @enzyme_mechanism` sites that are not wrapped in `eval` (test_types.jl:330, 366, 1773, 1859, 2315, 2351; test_dsl.jl:457, 461).

**Risks**
1. **Canonical Step Form and parameter names**, but only if the role model changes (section 3). Pure syntax cuts do not touch types.
2. **About 700 lines of mechanical fixture churn**, mostly in the 11,713-line enumeration tests, where CLAUDE.md requires inline fixtures.
3. **46 error-message pins** in test_dsl.jl.
4. **42 source-order test-macro uses** depend on the 3-tuple return (dsl.jl:751).
5. **Docs.** 32 macro examples across docs/src, README and HPC scripts need updating, and doctests run (docs/make.jl:11).
6. **Parse shape of the site line.** `regulatory_site(2): X::OnlyA, …` should parse like the existing `allosteric_regulators: A::OnlyA` and call-label lines (dsl.jl:689-698). Inference, not verified, since I was told not to run Julia.
7. **Open labels versus tag peeling** (section 3, item 4).
8. **Two ways to write ping-pong at the modeling level.** Hand-written fixtures use conformations (F, G); the enumerator uses residuals (mechanism_enumeration.jl:122-126). They are different `Mechanism` values. Inference: they also have different parameter names. "One way" here is a separate decision from the parser.