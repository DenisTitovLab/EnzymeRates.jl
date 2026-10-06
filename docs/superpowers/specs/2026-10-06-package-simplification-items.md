# Package simplification — approved item catalog
Generated from the review workflows of 2026-10-05/06 and Denis's sign-off page. Each item keeps the reviewer's change text, evidence and test impact. The design spec (2026-10-06-package-simplification-design.md) decides order and gates; this file supplies the detail. Line numbers refer to main at e59ebd6.

## src/types.jl
1. About 260 lines of src are reached only by tests. They are the dead Parameter flip and lookup family (87 lines), the Parameter predicates and Keq/Etot (11), and a tuple-accessor layer on the singleton types (~165) whose only src consumer is show. Each of those accessors re-lifts the whole Mechanism on every call.
2. Boilerplate that a loop or a generic form says once: 10 hand-written ==/hash pairs (69 lines, every hash stays bit-identical), a 16-function Sig codec (it becomes 4 encoders plus 1 decoder with a byte-identical layout and no throw-away Mechanism), and five _with_* rebuild helpers.
3. Validation runs two or three times. The 3-arg AllostericEnzymeMechanism constructor re-checks what AllostericMechanism already checked, at a cost of 2n+2 Mechanism re-lifts per allosteric compile. Both mechanism constructors also spell out the same Canonical Step Form sequence.
4. show walks name tuples, lifts the mechanism 5 times and regroups steps through a Dict. It can instead walk the Steps of one lifted mechanism and print the same strings.
5. Several fixes cut test time. metabolites(::AllostericEnzymeMechanism) has a runtime body on the loss! hot path. _rep_step rebuilds the free-enzyme set on every name(p, m) call. The canonical sort renders its key on every comparison. The Tier-2 direction rule does O(n^2) Species scans.
6. Behaviour-changing items are kept separate for Denis to decide. They are:
   - the 3-arg constructor removal;
   - deleting the test-oracle accessors versus compacting them;
   - RegulatorMults multiplicities (DSL syntax and Sig shape);
   - the bottomless-segment side-map bug fix;
   - canonical ligand order in RegulatorySite;
   - the Parameter type collapse from 6 to 3;
   - dropping the chain show.
7. Net effect: about 530 lines without any decision, about 720 if every decision is yes (1945 → ~1225 lines).

### types/T1 — Delete the dead Parameter flip and lookup helpers
- Net lines removed (this file group): 87
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** Delete types.jl:1819-1905 together with their docstrings:
- _flip_to_inactive (2 methods) and _force_inactive (2 methods);
- _param_for_symbol (2 methods);
- _onlyA_parameters_for_sym and _all_params_for_sym.
No src code calls any of them; the only src hits are their own definitions and internal calls. Their docstrings name a 'synth-dep machinery' that no longer exists. In the _emit_cat_params_for_rep docstring (1927-1930), remove the two deleted walkers from the caller list. This also removes _all_params_for_sym's near-duplicate of _all_i_state_parameters (rate_eq_derivation.jl:1247-1260).

**Other files and tests.** src: none.

Tests:
- Delete testset test_types.jl:1770-1786 ('_force_inactive forces :I'); it tests only the dead helper.
- Rewrite _dep_struct_key (test_rate_eq_derivation.jl:2421-2432) as `only(q for q in EnzymeRates._enumerate_parameters_full(mech) if EnzymeRates.name(q, mech) == sym)`. It is only ever called with a plain Mechanism, so its catch branch (_I_ to _A_) can never run.
- Reword the comments at test_types.jl:1723 and test_rate_eq_derivation.jl:2338.

Optional test-time saving, which needs Denis's sign-off under the test-deletion rule: test_types.jl:1719-1768 duplicates test_rate_eq_derivation.jl:2333-2378. Both use the same PK mechanism and the same assertions and differ only in RNG seed, so the test_types.jl copy can go.

**Evidence.** I read 1819-1905; it is 87 lines. Four verifiers ran rg over src, test and docs and found only the definitions, the error strings inside them and one docstring mention (1929). Tests reach them only through the dead-helper testset and _dep_struct_key, whose only caller is gated on `spec.mechanism isa EnzymeMechanism`.

### types/T2 — Delete the test-only Parameter API: governing_step, is_i_state, the two Union aliases, and the Keq/Etot structs
- Net lines removed (this file group): 11
- Risk: low
- Behaviour change: none (Keq(), Etot(), governing_step and is_i_state are unexported and only tests call them)
- Depends on: —

**Change.** Delete:
- types.jl:266-267 (the Keq and Etot structs). Keep Lallo at 268, which rate_eq_derivation.jl:1302 constructs, and reword the comment at 265 so it describes Lallo alone.
- types.jl:270-277: StepBoundParameter, StatefulParameter, governing_step, is_i_state and their comment.
- types.jl:1815-1816: name(::Keq, _) and name(::Etot, _).
src writes :Keq and :E_total as literal Symbols everywhere and never constructs Keq() or Etot(). The is_i_state hits in thermodynamic_constr_for_rate_eq_derivation.jl:397-454 and rate_eq_derivation.jl:1327 are an unrelated keyword argument.

**Other files and tests.** src: none.

Tests: delete the asserts at test_types.jl:925-945 (governing_step/is_i_state), 949-953, 1601-1602, 1609-1610 and 1696-1697 (Keq()/Etot()). Keq and Etot can stay in the AST-walker regex (test_types.jl around 2417); leaving them is harmless.

**Evidence.** grep for governing_step|is_i_state(|StepBoundParameter|StatefulParameter|Keq()|Etot() over src finds only the definitions. In test/ they appear only at test_types.jl:925-954, 1601-1611 and 1696-1697.

### types/T3 — Make the name(p, m) chokepoint concrete-only and compute the free-enzyme set in _rep_step only when needed
- Net lines removed (this file group): 5
- Risk: low
- Behaviour change: name(p, ::EnzymeMechanism) and name(p, ::AllostericEnzymeMechanism) now raise a MethodError. `name` is unexported and only tests call these methods. Rendered names are unchanged. Derivation is faster, because name() sits in the inner loops of _assemble_constraints (generator time only; rate_equation runtime is untouched).
- Depends on: —

**Change.** (a) Delete the two singleton _rep_step methods (types.jl:1784-1786). Narrow `_AnyMech` (1789-1790) to a one-line `Union{Mechanism, AllostericMechanism}`, and narrow the Kreg name method's union (1811) to `::AllostericMechanism`. Every src call of name(p, x) passes a concrete Mechanism or AllostericMechanism: rate_eq_derivation.jl:48, 70, 637, 697, 710, 713, 1196 and thermodynamic_constr:24, 410.
(b) In _rep_step (1777-1783), find the step's group first and compute the free-enzyme set only for multi-step groups: `step in group || continue; return length(group) == 1 ? only(group) : _group_rep(group, _free_enz_set(m))`. _group_rep on a one-element group returns that element, so the returned step is identical. Today every name(p, m) call rebuilds _free_enz_set: three passes over all steps, rendering species names as Strings.

**Other files and tests.** src: none; every caller already passes a concrete mechanism.

Tests: rewrite test_types.jl:1448-1490 ('name(p, ::AllostericEnzymeMechanism) chokepoint overloads') and the em lines at 1605-1611 to pass AllostericMechanism(aem) or Mechanism(em). The rendered names are the same, so coverage of the rendering is kept.

**Evidence.** _group_rep (thermodynamic_constr:120) is an argmin over the group, so a single-element group returns its only step whatever fes is. I read the chokepoint at 1776-1811. The AST walker checks Symbol literals only, and the signature change does not touch rendering.

### types/T4 — Compact the parameter-emit helpers in types.jl
- Net lines removed (this file group): 12
- Risk: low
- Behaviour change: none (same Parameter types and order, which feeds the column order of the constraint solve)
- Depends on: T1

**Change.** Three rewrites:
- _emit_cat_params_for_rep (1937-1945) becomes one ternary expression: `is_equilibrium(rep) ? Parameter[(is_binding(rep) ? Kd : Kiso)(rep, state)] : is_binding(rep) ? Parameter[Kon(rep, state), Koff(rep, state)] : Parameter[Kfor(rep, state), Krev(rep, state)]`. The element order stays identical.
- _enumerate_parameters_full (1915-1922) becomes `fes = _free_enz_set(m); Parameter[p for g in steps(m) for p in _emit_cat_params_for_rep(_group_rep(g, fes), :None)]`.
- _state_tag (1760-1767) becomes a ternary chain that keeps the error for unknown states.
Keep all docstrings; T1 already updates the caller list. Optional: once _emit_cat_params_for_rep serves any step, not only group representatives, rename it to `_step_params`. That rename crosses files.

**Other files and tests.** This lets the derivation and thermo groups write _step_parameters (thermodynamic_constr:36-46, about -9) and _state_step_params (rate_eq_derivation.jl:1139-1152, about -4) as comprehensions over this helper; they contain the same four-way switch token for token. Use typed `Vector{Parameter}[...]` comprehensions there. No test changes.

**Evidence.** thermodynamic_constr:39-42 and rate_eq_derivation.jl:1145-1148 repeat the same RE/SS × binding switch. The verifier counted about 13 lines in types.jl, and I get 12 after allowing for line wraps at 92 characters.

### types/T5 — Generate ==/hash for the value structs with one @eval loop, keeping the explicit Species pair
- Net lines removed (this file group): 54
- Risk: low
- Behaviour change: none (== semantics are identical and every hash value is unchanged)
- Depends on: —

**Change.** Replace these hand-written pairs: Residual (46-49), RegulatorySite (146-151), Step (234-240), the Kd..Krev loop and Kreg (278-287), ReactantAtoms (311-314), RegulatorMults (336-343), EnzymeReaction (483-491), Mechanism (850-853) and AllostericMechanism (948-959).

Put one loop after the AllostericMechanism struct (after line 928):
`for T in (Residual, RegulatorySite, Step, Kd, Kiso, Kon, Koff, Kfor, Krev, Kreg, ReactantAtoms, RegulatorMults, EnzymeReaction, Mechanism, AllostericMechanism); fs = fieldnames(T); @eval Base.:(==)(a::$T, b::$T) = $(foldr((f, acc) -> :(a.$f == b.$f && $acc), fs; init = true)); @eval Base.hash(x::$T, h::UInt) = $(foldl((acc, f) -> :(hash(x.$f, $acc)), fs; init = :(hash($(QuoteNode(nameof(T))), h)))); end`
Add a 2-line comment explaining why: these types hold Vector fields, so the default == would compare identity.

Keep the explicit Species pair (86-90). Its hash nests residual innermost and its == compares the cheap fields first. Keeping it makes every hash value bit-identical to today's. Use this unrolled @eval form rather than the ntuple/Union `_fields` variant, so the generated code on the enumeration hot path is the same as today's hand-written code.

**Other files and tests.** src: none; no ==, hash or isequal is defined outside types.jl. Tests: none. 'Step == / hash are structural' (test_types.jl:848) and the dedup tests are the gate.

**Evidence.** I read all ten hand-written pairs. Each compares and hashes every field. Each hash folds fields in declaration order on a hash(:T, h) seed, except Species. No struct has a cache field to exclude, and the package has no precompile workload that would run before the loop is defined.

### types/T6 — Compact the Sig codec and fold _drop_unbound_regulators into _sig_of (Sig layout byte-identical)
- Net lines removed (this file group): 85
- Risk: medium
- Behaviour change: none (every compiled type is identical). compile_mechanism no longer builds a throw-away EnzymeReaction plus a second Mechanism (full canonicalization and asserts) whenever a declared regulator is unbound, which is common for init seeds that declare dead-end inhibitors.
- Depends on: —

**Change.** Replace the 7 _to_sig methods (1014-1054), the 9 _*_from_sig decoders (1056-1116), _sig_of/_mechanism_from_sig (1118-1126) and _drop_unbound_regulators (1163-1193) with the following.

Four encoder methods:
- Metabolite leaf: `(nameof(typeof(m)), name(m))`;
- Vector{<:Metabolite} → Tuple;
- Species: `(bound, conformation, (added, subtracted))`;
- Step: the existing 5-slot tuple.

`_sig_of(m)` builds the reaction tuple inline and takes over _drop_unbound_regulators' docstring:
- reactants as `(leaf, Tuple(Tuple.(atoms)))`;
- only the regulators bound by some step, as `(leaf, Tuple(allowed_multiplicities))`;
- the catalytic multiplicities.

One `_mechanism_from_sig((reaction_sig, steps_sig))` decoder with local met/mets/species closures. Its NamedTuple tag lookup keeps today's 4-type whitelist.

`EnzymeMechanism(m::Mechanism) = EnzymeMechanism{_sig_of(m)}()`.

Keep the leaf-rules header comment (1004-1012) but drop the temporal 'replaces' wording (1015-1018). Use typed comprehensions (Metabolite[...], Vector{Step}[...]) so empty tuples never produce Vector{Union{}}.

**Other files and tests.** src: update the comment at mechanism_enumeration.jl:3361 that names _drop_unbound_regulators.

Tests:
- test_types.jl:94-100 (_to_sig and _metabolite_from_sig) and 1887 (_step_from_sig) call deleted helpers. Fold them into the round-trip testset at 1492-1534.
- Extend that fixture with a Theorell–Chance step, a competitive inhibitor, a bound allosteric regulator and an unbound regulator, so every leaf kind stays covered.
- The 7 pinned Sig strings (test_identify_rate_equation.jl:1261-1480, test_pivot_priority_regression.jl:17) must still eval.

**Evidence.** The verifier checked byte-identity leaf by leaf. `Tuple(p::Pair)` gives (first, second) with Int64 leaves. Re-canonicalizing an already canonical mechanism is idempotent, because direction tiers read only substrates/products and the stable re-sort of sorted keys is the identity. No src code indexes into a Sig.

### types/T7 — Delete the 3-arg AllostericEnzymeMechanism constructor; build the singleton type directly in the AllostericMechanism lift
- Net lines removed (this file group): 52
- Risk: low
- Behaviour change: The 3-arg method of the exported type AllostericEnzymeMechanism goes away. It is undocumented, and the type's docstring says to construct through AllostericMechanism. An empty-ligand RegulatorySite now errors when it is constructed instead of at compile_mechanism. Compiled types are identical, and every allosteric compile skips 2n+2 Mechanism reconstructions, which should help allosteric test time.
- Depends on: —
- Approved by Denis on the sign-off page.

**Change.** Delete AllostericEnzymeMechanism(cm, cat_sites, reg_sites) (1215-1263). Its only src caller is the 1-arg lift at 1312; dsl.jl:1248 and mechanism_enumeration.jl:1183 already go through the 1-arg lift.

Why it is redundant:
- AllostericMechanism already enforces the multiplicity, cat_allo_states length and tags, including :OnlyI (881-894).
- RegulatorySite already enforces lengths, multiplicity and states (129-138).
- The 'group numbers are 1..n consecutive' check cannot fail.
- It costs 2n+2 full Mechanism re-lifts per call, because each kinetic_group(cm, i) and n_steps(cm) call decodes the Sig and reruns the constructor.

Rewrite AllostericEnzymeMechanism(am) (1302-1313) as `cm = EnzymeMechanism(Mechanism(reaction(am), steps(am)))`, returning `AllostericEnzymeMechanism{typeof(cm), (catalytic_multiplicity(am), Tuple(cat_allo_states(am))), Tuple((Tuple(name.(ligands(s))), multiplicity(s), Tuple(allo_states(s))) for s in regulatory_sites(am))}()`. The type parameters are identical to today's.

Move the one check found nowhere else into the RegulatorySite constructor (127-140): `isempty(ligands) && error("RegulatorySite: needs at least one ligand")`.

Fallback if Denis keeps the method: replace 1222-1229 with `n_groups = length(steps(Mechanism(cm)))`. That saves about 7 lines and removes the quadratic re-lifts and the check that cannot fail.

**Other files and tests.** src: none.

Tests: rebuild the 11 3-arg call sites (test_types.jl:169, 175, 514, 518, 522, 526, 532, 635, 1417, 1458; test_rate_eq_derivation.jl:2255) through AllostericMechanism or @allosteric_mechanism; the allo_from_source helper already does this.
- Error-path cases retarget the existing AllostericMechanism and RegulatorySite tests; they still raise ErrorException.
- The empty-ligand cases (test_types.jl:522, test_rate_eq_derivation.jl:2255) move to the RegulatorySite testset (around 764).
- 635, the parenthesised-group show test, is rebuilt via @allosteric_mechanism so it still exercises that branch.

Type literals `AllostericEnzymeMechanism{CM,CS,RS}()` (CSV mechanism_type strings) are unaffected.

**Evidence.** I read 1215-1263 (49 lines) and the lift at 1302-1313. The verifier confirmed that the only check found nowhere else is the empty-ligand list (1248-1250). My count: -50 for the constructor and its blank line, -5 for the lift, +1 for the RegulatorySite check.

### types/T8 — Rewrite show for EnzymeMechanism and AllostericEnzymeMechanism over the Steps of one lifted mechanism
- Net lines removed (this file group): 80
- Risk: medium
- Behaviour change: None for every pinned case. Output differs only when an enzyme-form name equals a metabolite name; there today's enz_set lookup can throw or pick the wrong item. Printing lifts the mechanism once instead of five times.
- Depends on: —

**Change.** Replace 1348-1510 with the following.

`_step_text(s::Step)` prints the entry side, the arrow and the exit side. It replaces the twice-defined _arrow and the three join triples.

`_step_chain(flat::Vector{Step}, subs::Set{Symbol})` returns the chain string or nothing, keeping today's comments:
- empty flat → nothing;
- any form with degree > 2 → nothing;
- start at the first degree-1 form in from1, to1, from2, ... order, else the form named :E, else forms[1];
- the first step prefers one whose entry metabolites include a substrate;
- a later step with a non-empty entry side that is _is_chemistry → nothing;
- walk the steps with popat!.

show(em) lifts once with `m = Mechanism(em)`. It prints the chain or the '(n steps, k enzyme forms):' listing via _step_text, then '| regulators:' from regulators(reaction(m)).

show(aem) lifts once with `am = AllostericMechanism(aem)`. It prints each group as `body :: tag` or `(body) :: tag` straight from steps(am) and cat_allo_state(am, g), which deletes _format_allo_step_groups' Dict regrouping. Sites print from the RegulatorySite fields.

The output strings are identical for every pinned case. (After T7 every aem comes from a canonical AllostericMechanism, so the lift re-sorts nothing.)

**Other files and tests.** src: none. Tests: none needed; the pins at test_types.jl:205-330 and 537-650 are the gate. If no pin exists for the 0-step mechanism (constructible per test_types.jl:359), add one: 'EnzymeMechanism (0 steps, 0 enzyme forms):'.

**Evidence.** show(em) calls reactions, regulators, enzyme_forms, substrates and steps(Mechanism(m)), and each starts by re-lifting the Sig. I traced each tie-break (degree counting, start form, substrate preference, chemistry break) against the Step walk. The rewrite comes to about 72 lines, against 163 today.

### types/T9 — Delete the accessors no src code calls and rewrite their test uses as one-liners
- Net lines removed (this file group): 78
- Risk: low
- Behaviour change: none for package behaviour (unexported internals, not mentioned in docs)
- Depends on: T7, T8

**Change.** Delete:
- n_steps and rep_step on Mechanism and AllostericMechanism (847-848, 937-938). rep_step also misleads: it returns first(group), not the naming representative _group_rep.
- On EnzymeMechanism: equilibrium_steps, n_steps, kinetic_group, kinetic_groups and both steps_in_group methods (1613-1645). The Val{G} method has no caller anywhere.
- n_states (1670-1672) and cat_allo_state(::AllostericEnzymeMechanism) (1678-1683).
- The kinetic_group and steps_in_group forwards, plus the :equilibrium_steps/:n_steps/:n_states/:kinetic_groups entries of the forwarder loop (1686-1693).
- allosteric_regulators, regulatory_sites and reg_allo_state on AllostericEnzymeMechanism (1724-1745).

Compact, and keep, allosteric_regulators(::AllostericMechanism) (940-946) as `unique(AllostericRegulator[l for s in regulatory_sites(m) for l in ligands(s)])`. It has no caller today, but the derivation and enumeration groups propose reusing it (_kcat_forward's all_ligs loop, the regulator moves).

Keep kinetic_groups on the concrete types (src calls at mechanism_enumeration.jl:1273 and 2003), catalytic_mechanism/catalytic_multiplicity(aem) and both metabolites methods.

**Other files and tests.** src: none.

Tests: about 65 mechanical rewrites.
- n_steps (25 uses) → sum(length, steps(...)).
- rep_step (15) → first(steps(m)[g]).
- kinetic_group (6), mostly `length(unique(kinetic_group(m,i) for i in 1:n_steps(m)))` → length(steps(Mechanism(m))).
- kinetic_groups on em (test_types.jl:28).
- equilibrium_steps (3), n_states (4), steps_in_group (1), reg_allo_state (1), allosteric_regulators(aem) (3).
- cat_allo_state and regulatory_sites on aem (test_types.jl:151-153, 555-559).

test_structure (test_rate_eq_derivation.jl:800-801) reads n_steps and n_states for every MECHANISM_TEST_SPEC, allosteric ones included. Add a 3-line _testhelper_ that lifts an em or aem to its catalytic Mechanism. No coverage of package behaviour is lost.

**Evidence.** I grepped src outside types.jl. None of these is called, except cat_allo_state and regulatory_sites, which are always called on a concrete `am` (rate_eq_derivation.jl:708, 1103, 1143, 1251, 1292, 1401, 1428, 1499, 1550). kinetic_group(em) and n_steps(em) feed only the 3-arg constructor (T7), and cat_allo_state/regulatory_sites(aem) feed only show (T8). I left out the concrete kinetic_groups deletion: 2 src lines are not worth about 17 test edits. Counted line by line, blank lines included.

### types/T10 — Delete the test-oracle tuple accessors: reactions, enzyme_forms, and substrates/products/regulators on the singletons
- Net lines removed (this file group): 86
- Risk: medium
- Behaviour change: none for package behaviour (unexported, and no docs page mentions them)
- Depends on: T8
- Approved by Denis on the sign-off page.

**Change.** After T8, the only src caller left is products(M()) at rate_eq_derivation.jl:784. Replace it with `name.(products(reaction(Mechanism(M()))))`.

Then delete:
- substrates, products and regulators on EnzymeMechanism (1532-1549);
- reactions (1595-1612) and enzyme_forms (1646-1669);
- the remaining AllostericEnzymeMechanism forwarder loop (1684-1689);
- regulators(::AllostericEnzymeMechanism) (1695-1711), whose docstring falsely says rate-equation code reads it.
Shrink the false accessor header comment (1512-1517: 'Matrix{Int} stoichiometry', 'contract consumed by the @generated body builders').

Fallback if Denis keeps them: rewrite the bodies as one-liners with unique/setdiff (same tuples, same order) and correct the regulators(aem) docstring. That saves about 40 lines.

**Other files and tests.** src: a one-line swap at rate_eq_derivation.jl:784 (0 net).

Tests (about 45 sites):
- The QSSA/ODE oracles reference_qssa and build_ode_rhs (test_rate_eq_derivation.jl:23-30, 610-640, 684-730) and test_kcat_rescaling (1097-1099) walk _flat_steps(Mechanism(m)), or use about 20 lines of _testhelper_reactions/_testhelper_enzyme_forms.
- The accessor asserts in test_types.jl:4-55 and test_dsl.jl:88, 143, 159-160, 540, 645, 718 are rewritten against Mechanism(m).

**Evidence.** grep -rnw over src after T8 finds only products(M()) at rate_eq_derivation.jl:784. The names are unexported, and docs/README contain no references. Line spans counted from the file: 18+18+24+6+17, plus 4 from the header.

### types/T11 — One metabolite-name list per reaction, and make metabolites(::AllostericEnzymeMechanism) @generated
- Net lines removed (this file group): 16
- Risk: low
- Behaviour change: Same tuple values in the same order. Allosteric loss! should become type-stable and allocation-free. This is unmeasured, so the TDD allocation test is the gate; it may cut allosteric fit and identify test time.
- Depends on: —

**Change.** Add `_metabolite_names(rxn::EnzymeReaction) = unique!(Symbol[name.(substrates(rxn)); name.(products(rxn)); [name(regulator(rm)) for rm in regulators(rxn)]])`. It is an internal name, so the exported metabolites gets no new method.

The @generated metabolites(::EnzymeMechanism) body (1579-1592) becomes `Tuple(_metabolite_names(reaction(Mechanism(EnzymeMechanism{Sig}()))))`. Keep its comment, minus the 'unlike the other demoted accessors' wording.

Make metabolites(::AllostericEnzymeMechanism) (1712-1722) @generated as well: `Tuple(unique!(Symbol[metabolites(CM())..., (l for (ligs, _, _) in RS for l in ligs)...]))`. Today it returns `(cat_mets..., extra...)` with extra::Vector built through Set/push!. So loss! (fitting.jl:126-138) builds NamedTuple{MetNames} from a runtime-length tuple at every data point of every allosteric fit.

**Other files and tests.** src:
- rate_eq_derivation.jl:202-207: _concentration_symbols becomes `Set(_metabolite_names(reaction(mech)))` (about -5).
- identify_rate_equation.jl:47-51: mnames becomes `_metabolite_names(reaction)` (about -3). This adds dedup to a check that only tests whether columns are present, which is harmless.

Tests (TDD): first add a failing `@allocated loss!(x, fp) == 0` for an allosteric FittingProblem in test_fitting.jl. Today only uni_uni is checked (270-290).

**Evidence.** Read 1574-1593 and 1712-1722. The EnzymeMechanism method's own comment says a runtime body would make loss!'s NamedTuple type-unstable, and the AllostericEnzymeMechanism method is exactly that kind of runtime body. Line count: helper +4, em body -13, aem body -8.

### types/T12 — Replace the five _with_* rebuild helpers with one keyword _with
- Net lines removed (this file group): 40
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** Replace _with_steps(::Mechanism) (855-861) and _with_steps/_with_cat_allo_states/_with_reg_sites/_with_steps_and_cat_states (961-1002) with one docstring and two methods:
- `_with(m::Mechanism; groups = steps(m)) = Mechanism(reaction(m), groups)`;
- `_with(am::AllostericMechanism; groups = steps(am), states = cat_allo_states(am), sites = regulatory_sites(am)) = AllostericMechanism(reaction(am), groups, states, catalytic_multiplicity(am), sites)`.
The keyword names must not shadow accessor names: a default like `steps = steps(am)` hits an undefined local. The defensive copy() calls go. The constructor never aliases its inputs: _canonicalize_step_directions builds fresh vectors, states go through permute!(copy(...)), and sites through sort(...).

**Other files and tests.** src, mechanism_enumeration.jl:
- Rename the call sites at 1293, 2143, 2148, 2905, 2920, 2932 and 3020.
- The manual rebuild at 2816-2818 becomes `_with(am; sites = new_sites)` (about -2).
- Keep _dead_end_child's explicit constructor at 2425, which passes a different reaction.
- The _make_am_with_added_reg rewrite belongs to the enumeration group.

Tests: rename test_mechanism_enumeration.jl:9434 and 9459.

**Evidence.** The helpers take 49 lines (7+42), and the replacement takes 7. Verifiers confirmed the constructor never aliases its inputs (types.jl:895-907). They also refuted replacing 2425, because it uses a different reaction, and caught the keyword-shadowing bug in one sketch.

### types/T13 — One canonical-groups routine for both mechanism constructors, with precomputed sort keys
- Net lines removed (this file group): 12
- Risk: high
- Behaviour change: None. The keys are the same, sortperm is stable on Julia ≥1.9 (the package's compat), and the operations run in the same order. An empty group now raises an explicit error instead of a BoundsError; no test pins that.
- Depends on: —

**Change.** Fold _canonical_group_order! (629-641), the canonicalize/permute!/assert sequence of both constructors (826-830, 895-904) and _assert_re_segments_have_bottom (835-842) into `_canonical_groups(reaction, groups) -> (groups, perm)`. Its body, in today's order:
- an explicit 'empty kinetic group' error;
- `gs = [g[sortperm(_step_canonical_key.(g))] for g in _canonicalize_step_directions(reaction, groups)]`;
- `perm = sortperm([_step_canonical_key(first(g)) for g in gs]); gs = gs[perm]`;
- _assert_uniform_groups, _assert_each_reaction_once, then the bottomless-segment error with the same message.

The Mechanism constructor becomes `new(reaction, first(_canonical_groups(reaction, steps)))`. The AllostericMechanism constructor validates tags first, as today, then runs `cat_steps, perm = _canonical_groups(...)` and `cat_allo_states = cat_allo_states[perm]`.

Keep _canonicalize_step_directions as a separate function, because test_rate_eq_derivation.jl:237 calls it. Merge the docstrings.

**Other files and tests.** src: none. Tests: none. The gates are the canonicalization tests (test_types.jl:1234-1306, 2264), the fitted_params pins and the enumeration count pins.

**Evidence.** types.jl:826-830 and 895-904 run the same five calls. sort!(by=) re-renders the key, including String allocations from name(::Species), on every comparison. Precomputing renders each key once, and the constructor runs for every enumeration child and every Mechanism(em) lift. Rated high only because it touches invariant 1; it is pure code motion plus key caching.

### types/T14 — Precompute the Tier-2 entry kinds once and delete _entry_kind
- Net lines removed (this file group): 10
- Risk: medium
- Behaviour change: none (each form gets the same booleans; a Species-keyed Dict uses the same ==/hash as `form == sp`)
- Depends on: —

**Change.** In _canonicalize_step_directions (597-604), build a `kind::Dict{Species, Tuple{Bool,Bool}}` once over all steps. For each step and each pair (form, free) in ((from, consumed), (to, released)), OR in (any substrate in free, any product in free). Pass `kind` to _canonical_step_direction in place of all_steps. Tier 2 (581-584) keeps s when `get(kind, f, (false,false)) == (false, true)` and `get(kind, t, (false,false)) == (true, false)`, and flips in the mirrored case. Delete _entry_kind (536-552) and its four-symbol encoding. Move its docstring (516-535), including the 'Why ALL steps' paragraph, onto the table construction.

**Other files and tests.** src: none. Tests: none; the same gates as T13 apply.

**Evidence.** _entry_kind (536-552) scans every step with Species == on each call, and Tier 2 calls it twice per Tier-1 tie. That is O(n^2) in every Mechanism/AllostericMechanism constructor call. The only callers are at types.jl:581-582. This subsumes part (c) of the Species/Step proposal (_entry_kind returning a tuple).

### types/T15 — Small rewrites in the Species/Step layer
- Net lines removed (this file group): 16
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** (a) One `name(m::Metabolite) = m.name` replaces the four per-subtype methods (24-27). All four subtypes have a name::Symbol field. Confirm that JET (test_aqua_jet.jl) stays clean.
(b) name(::Species) (100-116) becomes `parts = [String(conformation(s)) * join(_met_label.(bound(s)))]`. When there is a residual, push "res", append `"+" .* String.(name.(added(r)))` and `"-" .* String.(name.(subtracted(r)))`, and join with "_". The string is the same.
(c) `_flat_steps(m::Mechanism) = [(s, g) for (g, group) in enumerate(steps(m)) for s in group]` replaces 1519-1529. It returns the same Vector{Tuple{Step,Int}}, even when empty. Fix its docstring, which says 'Mechanism(em)' although the function takes a Mechanism, and drop the doubled blank line at 1530-1531.

**Other files and tests.** none (_flat_steps has 18 src callers, and its return type is unchanged)

**Evidence.** Read each function. The verifier confirmed that (a) holds for all four subtypes, that (b) produces the same parts joined by '_', and that (c) infers the same element type. Part (c) of that proposal moves to T14 and part (d) to T17.

### types/T16 — Simplify constructor bodies and validation loops, keeping every error message
- Net lines removed (this file group): 28
- Risk: low
- Behaviour change: none (validation conditions and message texts are kept; the reviewer's signed-balance variant would have changed the atom-imbalance text)
- Depends on: —

**Change.** (a) EnzymeReaction constructor (404-467):
- build sub_names/prod_names directly as typed comprehensions and drop subs/prods;
- use allunique for the three duplicate checks;
- compute per-side atom totals with `mergewith(+, Dict{Symbol,Int}(), (Dict(atoms(ra)) for ra in sorted_reactants if metabolite(ra) isa T)...)` and keep today's atom-imbalance message word for word.
(b) show(::EnzymeReaction) line 505: drop the unreachable `|| (length(mults) == 1 && mults[1] != 1)` disjunct (0 lines).
(c) Kreg-collision check (911-921): `ligs = [name(l) for s in regulatory_sites for l in ligands(s)]; dup = findfirst(i -> ligs[i] in view(ligs, 1:i-1), eachindex(ligs))`. This reports the same ligand as today.
(d) AllostericMechanism(aem) (1277-1292): build sites with a typed `RegulatorySite[...]` comprehension. An untyped one gives Vector{Union{}} for regulator-free mechanisms and then a MethodError.
(e) Define reaction once, and kinetic_groups once, over `Union{Mechanism, AllostericMechanism}` (844-846, 930, 936).
(f) Turn the RegulatorySite state loop (134-138), the ReactantAtoms count loop (300-304) and the AllostericMechanism tag loop (888-894) into findfirst plus one error that keeps the offending value or group index.

**Other files and tests.** none

**Evidence.** Read 404-467, 503-508, 888-921 and 1277-1292. mults is sorted, unique and ≥1, so the second disjunct is dead. allunique is equivalent to length(x) == length(Set(x)). No test pins any of these messages, but I keep them identical anyway.

### types/T17 — Add a _forms(groups) helper and return a segment index from _re_segment_extras (types.jl side of cross-file merges)
- Net lines removed (this file group): -3
- Risk: low
- Behaviour change: none (unique keeps first occurrences, and Species has ==/hash)
- Depends on: —

**Change.** Add `_forms(groups)`: the distinct Species of all steps in first-seen, from-then-to order, written `unique(sp for g in groups for s in g for sp in (from_species(s), to_species(s)))`. Place it next to _re_segment_extras. Use it for _re_segment_extras' species list (726-730), which removes an O(n^2) `in` scan per constructor, and for T8's enzyme-form count.

Extend _re_segment_extras to return each form's segment number and `idx` as well. Fill a `segment_of` vector during the BFS; it can replace the `haskey(expo, v)` visited test. This lets enumeration delete _indexed_re_segments. _bottomless_re_segment keeps destructuring three names. Update the docstring, and keep BFS order.

**Other files and tests.** The gain lands in other groups:
- mechanism_enumeration.jl:1393-1404 deletes _indexed_re_segments, and its callers at 1430, 1835 and 2205 call _re_segment_extras (about -12).
- _forms replaces the distinct-form loops in thermodynamic_constr (_free_enz_set 72-76, _enumerate_species_names 260-273) and rate_eq_derivation.jl (1056-1060, 1106-1109, and _enumerate_species 184-193 unless the RE-segment unification deletes it), about -10.
Tests: none.

**Evidence.** _re_segment_extras already builds idx (731), and its BFS knows each form's segment. In types.jl this adds about 3 lines (+4 for the helper, -4 for the loop, +3 for segment_of and the docstring) so that about 22 lines can go elsewhere.

### types/T18 — Fix the remaining false comments and document the lossy Mechanism(em) lift
- Net lines removed (this file group): 0
- Risk: low
- Behaviour change: none (comments only)
- Depends on: T6, T9

**Change.** - The Kreg docstring (255-258) uses the retired 'R or T branch' wording; change it to the A/I terminology.
- Give `Mechanism(em::EnzymeMechanism)` (1195, currently undocumented) a short docstring. It should say the lift is not the inverse of EnzymeMechanism(m): shared_catalytic_site, RegulatorMults.reg_type and unbound regulators are not encoded, on purpose, so that equivalent mechanisms share one compile.

Other items remove the rest of the falsehoods: 1351 'both Sig shapes' (T8), 1015-1018 temporal 'replaces' (T6), 1519 'Mechanism(em)' (T15), 1927-1930 walker list (T1), and 1512-1517 plus 1695-1700 (T10). If T10 is declined, fix those last two here.

**Other files and tests.** none

**Evidence.** Each cited comment is provably false or temporal. For example, rg 'stoichiometry' src finds only the 1515 comment. _to_sig(::EnzymeReaction) (1050-1054) omits shared_catalytic_site, and _to_sig(::RegulatorMults) (1045-1048) omits reg_type.

### types/T20 — Fix _bottomless_re_segment's side map so a reactant's inhibitor copy cannot hide a bottomless segment
- Net lines removed (this file group): 5
- Risk: medium
- Behaviour change: The constructor and the enumeration filters now reject mechanisms whose RE segment is bottomless only through a reactant's inhibitor copy, and the verdict no longer depends on order. Enumeration counts may drop in substrate-inhibition families. The 128,493-mechanism rxn6 depth-2 population passes today, so this is likely rare.
- Depends on: —
- Approved by Denis on the sign-off page.

**Change.** Replace the side-map loop (791-797) with `side = Dict(name(m) => typeof(m) for g in steps for s in g for m in Iterators.flatten((consumed(s), released(s))) if m isa Reactant)`. This classifies each name by its reactant role over all steps, RE and SS.

Today `side[name(m)] = typeof(m)` runs over RE steps with last-write-wins. A reactant's competitive-inhibitor copy has the same name, so if it binds at RE in a group that sorts later, it retypes the name as CompetitiveInhibitor. _re_segment_extras is keyed by name, so it then stops counting that reactant's extras. A segment that is bottomless at zero substrate is accepted, giving a 0/0 rate law. The verdict also depends on group order: enumeration calls it on non-canonical orders (mechanism_enumeration.jl:1288, 1784, 2012), while the constructor calls it on canonical ones.

The 'm isa Reactant over RE steps only' variant is incomplete. It misses a reactant whose own binding is steady state.

**Other files and tests.** src: none; the enumeration filters at mechanism_enumeration.jl:1288, 1784 and 2012 call this function.

Tests:
- TDD: first add a failing constructor fixture. The verifier's example: substrates A and B; E+A and E+B steady state; RE steps EA+B⇌EAB, EB+A⇌EAB and EB+A::Inh⇌EAinhB.
- Rerun the enumeration count pins.
- Rerun the allosteric derivation tests, because _state_allo_mechanism builds AllostericMechanisms for pruned I-state graphs.

**Evidence.** Read 790-807. extras is keyed by name (737-738), and name(::CompetitiveInhibitor) returns the bare name (27). The verifier built a concrete false negative in which the canonical key ('EB','EAB') < ('EB','EAinhB') lets the copy write last. Every declared reactant appears as a Reactant in some step, so the all-steps map is complete.

### types/T21 — Store RegulatorySite ligands in canonical (name) order
- Net lines removed (this file group): -1
- Risk: medium
- Behaviour change: Fewer duplicate allosteric children, so counts can drop. DSL sites written out of name order get reordered, which changes the RS type parameter and possibly the Kreg order in parameters(m). Kreg names do not change.
- Depends on: —
- Approved by Denis on the sign-off page.

**Change.** In RegulatorySite's inner constructor, use `p = sortperm(ligands; by = name); new(ligands[p], multiplicity, allo_states[p])`, and rewrite its docstring ('ordering is meaningful (canonicalize at the call site if needed)').

Today _make_am_with_added_reg appends the new ligand (mechanism_enumeration.jl:2805-2806), while the merge move sorts (3008-3012). Different add or merge routes therefore produce AllostericMechanisms that are physically identical but differ under ==. Both pass the parent dedup and the hash skip (identify_rate_equation.jl:521-527, 725) and are compiled and fit.

**Other files and tests.** src: the sortperm block in the merge (mechanism_enumeration.jl:3008-3012) is deleted, about -3.

Tests:
- TDD: add a test that the two routes give == mechanisms.
- Re-pin allosteric golden and parameter-order tests whose DSL sites list ligands out of name order.

An alternative with no DSL effect: sort only in _make_am_with_added_reg (enumeration group).

**Evidence.** types.jl:127-139 does no ordering, and RegulatorySite ==/hash compare ligands by position (146-151). Memory notes list 'within-site canon pre-existing' as an open flag. How many duplicates actually arise in a beam run is unmeasured.

### types/T22 — Collapse the six step-bound Parameter types into three (RE constant, forward rate, reverse rate)
- Net lines removed (this file group): 8
- Risk: medium
- Behaviour change: none (rendered names are identical)
- Depends on: T2, T4
- Approved by Denis on the sign-off page.

**Change.** Merge Kd and Kiso into one RE-constant type, Kon and Kfor into a forward rate, and Koff and Krev into a reverse rate (245-253). Within a pair the types differ only by is_binding(step), which _assert_uniform_groups keeps uniform within a group. The k-constant name methods already share methods (1797-1800). The RE constant renders reverse(sides) when is_binding(rep) and forward sides otherwise. The emit switch becomes two-way. Every rendered name stays identical, and the name(p, m) chokepoint remains the only renderer.

**Other files and tests.** src: the SS unions at rate_eq_derivation.jl:698 and 714, and thermodynamic_constr:39-42 unless T4's reuse lands. Tests: about 100 references (Kd 39, Kiso 13, Kon 14, Koff 10, Kfor 12, Krev 13). The developer docs mention the Kd/Kon/Koff naming.

**Evidence.** The type-level consumers are only the name methods (1797-1804) and the SS filter (rate_eq_derivation.jl:698, 714). The AST walker does not need the new type names. The gain is small once T2 and T4 have landed.

### types/T23 — Optional: drop the chain rendering of linear mechanisms in show
- Net lines removed (this file group): 33
- Risk: low
- Behaviour change: REPL and show output of linear mechanisms changes from the one-line chain 'EnzymeMechanism: E + S ⇌ ES <--> EP ⇌ E + P' to the multi-line listing.
- Depends on: T8
- Approved by Denis on the sign-off page.

**Change.** Delete T8's _step_chain and always print the '(n steps, k enzyme forms):' multi-line listing for EnzymeMechanism.

**Other files and tests.** Tests: update the chain-format pins at test_types.jl:216 and 250. The multi-line pins (233, 265) are unaffected.

**Evidence.** After T8 the chain walk is a self-contained helper of about 35 lines including comments. The multi-line fallback already renders every mechanism correctly.

## src/dsl.jl
src/dsl.jl (1256 lines) runs only at macro-expansion time: nothing in src calls it and it is off the hot path. test_dsl.jl takes 3.2 s of the 1152 s suite, so this group shortens code but does nothing for the 2x test-time goal.
Where the excess lines are:
(1) The _StepSideTerm round trip, about 450 lines (308-506, 661-833, 857-982). Each side term is parsed into an 8-field record, the records are walked again for the opaque-name check, regrouped through a Dict that does nothing, and only then turned into Exprs. Emitting the Exprs while parsing removes about 300 lines and constructs identical values.
(2) Reaction-Expr emission is written three times, plus an atoms Expr that is built only to be taken apart again (about 92 lines).
(3) The DSL repeats checks the constructors already make: the state-tag sets, ligand collisions and a duplicated _ALLOSTERIC_REG_STATES. The DSL copy has already drifted: it lets :OnlyI through on catalytic steps (about 44 lines).
(4) The label shape is parsed by three separate routines (_line_label_expr/_reject_allosteric_syntax!, _parse_labeled_line, _match_regulatory_site_line), and the positive-Int check is written out six times (about 95 lines).
Parse bugs found along the way:
- Repeated or conflicting declarations silently keep the last value: multiplicity labels, steps:, catalytic_multiplicity:, a name listed twice in allosteric_regulators:, and a repeated regulatory_site kwarg.
- The docstring's own `(step) :: Tag` example does not parse.
- A right-hand `X::Inh` release is rejected.
- An allosteric-only regulator written into a catalytic step is silently accepted.
- A wrong-role residual fails at run time with a MethodError.
Canonical Step Form is not affected: the DSL still hands source-order groups to the constructors, and the test helpers' source-order return values keep their shape.
Size: about -600 if every item is accepted (file to about 657 lines). About -195 with no decision needed. About -500 if Denis also approves the step-parser rewrite.
The one real API change is removing the per-regulator multiplicity syntax `A(1, 2)`, which nothing ever reads. It spans this file and types.jl, and it stops old regulator-bearing CSV mechanism_type strings from evaluating.

### dsl/D1 — One reaction-Expr emitter for all three macros; parse atoms straight to element => count pairs
- Net lines removed (this file group): 92
- Risk: low
- Behaviour change: None. The constructors receive equal values:
- EnzymeReaction sorts reactants and regulators and uniques the mults.
- ReactantAtoms sorts atoms.
- The mechanism macros now pass `:unspecified` and an empty shared_catalytic_site explicitly; those are exactly the constructor defaults (types.jl:322, 396).
- Depends on: —

**Change.** Replace the three reaction-Expr emission paths with two functions:
- `_reaction_expr(reactants, regulators, mults, shared)`, about 17 lines;
- `_mechanism_reaction_expr(subs, prods, regs)`, about 13 lines, which carries the placeholder-carbon docstring moved from dsl.jl:669-674.

Delete:
- `_build_shared_site_expr` (150-154);
- `_build_reactants_expr` (244-263);
- `_atoms_pairs_expr` (265-275);
- `_build_regulators_expr` (277-301), including its unreachable 'internal: unknown regulator kind' branch;
- `_build_catalytic_mults_expr` (303-306);
- the hand-rolled reactants/regulators/EnzymeReaction block inside `_build_mechanism_expr` (679-706), which becomes one call to `_mechanism_reaction_expr`.

Rewrite the parsers:
- `_parse_chemical_formula` (164-176) checks the whole string with one anchored element-token regex and returns a `Pair{Symbol,Int}` comprehension over `eachmatch`. It accepts the same strings as today's matched-length test.
- `_parse_atom_bracket_entries` (116-131) becomes a `map` that returns `(name, Pair{Symbol,Int}[...])`. The vector is typed, so `S[]` still reaches ReactantAtoms' no-atoms error.
- `_parse_reaction_block` collects `reactants` as `(type, name, atoms)` tuples.
- It passes the type names `:CompetitiveInhibitor`/`:AllostericRegulator` to `_parse_regulator_entries` in place of `:competitive`/`:allosteric`.
- It resolves the default multiplicities once, at its end: `something(mults, [1])` for the reaction and `something(ms, mults)` per regulator.
- The macro body (33-43) becomes two lines.

Emission splices the type name: `EnzymeRates.$t(...)`. Sketch: scratchpad/draft_dsl_all.jl:33-215.

**Other files and tests.** None. Every deleted helper is referenced only inside dsl.jl. test_dsl.jl:270-363 already covers atoms, multi-atom formulas, regulator kinds and shared sites.

**Evidence.** I read dsl.jl:1-307 and 661-733. The verifier confirmed P2 at 100 (range 90-105). I mapped current lines onto the draft:
- macro: -7
- atom brackets: -8
- shared-site emitter: -6
- formula: -9
- emitter block 243-307 (65 lines down to 33): -32
- mechanism reaction block: -33
- parse block: +3
Net: -92. This covers P27's _build_regulators_expr dead branch and P29(e,f)/P30(h), which become moot.

### dsl/D2 — One `_positive_int` check for the six positive-integer literals; compact the multiplicity and regulator parsers
- Net lines removed (this file group): 30
- Risk: low
- Behaviour change: Valid input parses the same. Only errors change:
- `catalytic_multiplicity: 1, 2` fails with an ArgumentError from `only` instead of the custom 'takes a single Int' message.
- The positive-Int messages become uniform.
- Depends on: D1

**Change.** Add `_positive_int(v, what) = v isa Integer && v >= 1 ? Int(v) : error(...)` with a one-line docstring. Use it at all six copies (dsl.jl:98, 208, 228, 234, 1081, 1195):
- `_parse_multiplicity_tuple` (223-242) shrinks to its length check plus an `Int[_positive_int(a, ...) for a in (tuple ? v.args : (v,))]` comprehension.
- `_parse_regulator_entries` (185-220) becomes a `map(values) do v ... end`; its call-form branch becomes one comprehension.
- `catalytic_multiplicity:` (1190-1198) becomes `cat_n = _positive_int(only(values), ...)`.
- The regulatory_site multiplicity check (1080-1084) uses the helper too.
Keep the `isa Integer` guard: it stops a literal like 2.0 from converting silently.

**Other files and tests.** None. test_dsl.jl:365-408 and 466-488 cover these paths, and no test pins their messages.

**Evidence.** `grep -n '>= 1' src/dsl.jl` finds 98, 208, 228, 234, 1081 and 1195. The verifier confirmed it at 28. Mapped onto the draft:
- helper: +4
- regulator entries: -13
- multiplicity tuple: -13
- catalytic_multiplicity: -6
- regulatory_site multiplicity: -3
Net: about -30.

### dsl/D3 — _parse_reaction_block: drop the label Set and merge the two multiplicity labels (fixes order-dependent last-wins)
- Net lines removed (this file group): 11
- Risk: low
- Behaviour change: - `oligomeric_state:` written after `allowed_catalytic_multiplicities:` already errors.
- The reverse order, and a repeated `allowed_catalytic_multiplicities:` line, silently keep the last value today; both would now error.
- `oligomeric_state: (2,)` becomes accepted as [2].
- The list of valid labels becomes a hardcoded string, a minor drift risk.
- Depends on: D1
- Approved by Denis on the sign-off page.

**Change.** - Delete `_VALID_REACTION_LABELS` (45-51) and its membership check (77-79).
- Add a final `else error(...)` that names the unknown label and lists the valid ones.
- Merge the `allowed_catalytic_multiplicities:` branch (89-90) and the `oligomeric_state:` branch (91-101) into one: `mults === nothing || error(...)`, then `mults = _parse_multiplicity_tuple(values, label)`, then `label === :oligomeric_state && length(mults) != 1 && error(...)`.

**Other files and tests.** Add TDD tests to test_dsl.jl:
- `oligomeric_state: 2` followed by `allowed_catalytic_multiplicities: (1, 2)` errors;
- two `allowed_catalytic_multiplicities:` lines error.
test_dsl.jl:410-418 (unknown label) still passes.

**Evidence.** Read 89-101: only the oligomeric_state branch checks `mults === nothing`. The verifier confirmed it at 12, and no fixture combines the two labels. P30(i), letting oligomeric_state reuse _parse_multiplicity_tuple, is this same change.

### dsl/D5 — _parse_shared_site_pairs as a map
- Net lines removed (this file group): 2
- Risk: low
- Behaviour change: None: a label always carries at least one value.
- Depends on: —

**Change.** Replace the push loop in `_parse_shared_site_pairs` (138-148) with `map(values) do v <same check>; (v.args[1]::Symbol, v.args[2]::Symbol) end`.

**Other files and tests.** None. test_dsl.jl:490-524 covers valid and malformed pairs.

**Evidence.** Single caller at line 103. The function goes from 11 lines to 9.

### dsl/D6 — Delete `_reject_allosteric_syntax!` and `_line_label_expr`; one label-shape parser; `only(values)` for step blocks
- Net lines removed (this file group): 47
- Risk: low
- Behaviour change: The same inputs are rejected; only errors change:
- Misplaced allosteric labels get different wording.
- Without the whole-block pre-scan, an earlier malformed line is reported before a later allosteric label.
- `catalytic_steps: a, b` raises an ArgumentError from `only`.
- Depends on: —

**Change.** - Delete `_reject_allosteric_syntax!` (557-572), `_line_label_expr` (574-590) and the call at 553.
- The plain body's final `else` (630) becomes one message: the unknown label, plus a hint that `allosteric_regulators:`, `catalytic_steps:`, `regulatory_site(...)` and the like belong in @allosteric_mechanism. Each of those already reaches this branch.
- Route `steps:` through `_parse_labeled_line` (`elseif label == :steps; steps_block = only(values)`) instead of the hand-matched special case at 618-621.
- Shorten `_parse_labeled_line` (597-611) to one shape test: split a leading tuple into head and rest, require `head` to be a `:(:)` call, and return `head.args[2], Any[head.args[3], rest...]`.
- In the allosteric body, `catalytic_steps:` (1203-1206) becomes `cat_steps_block = only(values)`; the 'multiple blocks' check stays.

**Other files and tests.** Delete test/mechanism_definitions_for_test_enzyme_derivation.jl:66, the `_reject_allosteric_syntax!` call in @enzyme_mechanism_src. test_dsl.jl:179 and 190 use `@test_throws Exception` and still pass.

**Evidence.** `_parse_plain_mechanism_body` accepts only substrates, products, regulators and steps. All four allosteric labels, and the regulatory_site(...) call label, therefore fall to the `else` at 630.

References: `_line_label_expr` has one caller. `_reject_allosteric_syntax!` is called at 553 and by the test helper.

The verifier confirmed P3 at 42. P26's tailored-message variant saves only 28.

Mapped onto the draft:
- 556-611 (56 lines down to 14): -42
- macro: -1
- catalytic_steps: -4
Net: -47.

### dsl/D7 — Step parser emits Step/Species/Metabolite Exprs while parsing; delete the _StepSideTerm round trip
- Net lines removed (this file group): 300
- Risk: medium
- Behaviour change: The constructed Mechanism and AllostericMechanism values are identical:
- canonicalization stays in the constructors;
- Species, Residual and Step sort their own contents;
- `Residual()` is `Residual(Substrate[], Product[])`;
- the current regroup is the identity.

Only errors change:
- Errors in @allosteric_mechanism steps now name @allosteric_mechanism; today they hardcode @enzyme_mechanism.
- A `::Tag` on a plain step gets a role-tag or no-enzyme-form error instead of 'tag annotation not allowed'.
- A comma-form `E(S, residual = ...)` reports 'not a declared metabolite' instead of 'invalid entry in species'.
- Errors fire in source order, not after the whole block has been parsed.

As drafted, the rewrite also fixes D9 and D10.
- Depends on: D1
- Approved by Denis on the sign-off page.

**Change.** Delete:
- 308-506: `_parse_step_side_terms`, `_step_side_term_info`, `_call_form_term_info`, `_StepSideTerm`, `_term_metabolite`/`_term_bare_enzyme`, `_reject_opaque_bound_forms`, `_walk_residual_expr`;
- the rest of `_build_mechanism_expr`: its Dict regroup and wrapper (661-668, 675-678, 707-733);
- 735-833: `_build_step_expr`, `_split_side`, `_species_expr_from_term`, `_metabolite_expr`;
- 857-953: `_parse_steps_block_with_groups`, `_step_struct_info`.

Replace them with:
- `_parse_steps_block(steps_block, role_of, macro_name; allow_tag)`. It returns the Vector{Vector{Step}} Expr and `tags::Vector{Symbol}` in source order. The return shape is fixed, with no Ref counter and no regroup.
- `_step_expr`. It checks for ⇌ / <--> and builds `Step(from, to, Metabolite[consumed], Metabolite[released], is_eq)`.
- `_side_exprs`. It splits the `+` terms and requires exactly one enzyme form. It keeps the pinned 'no enzyme-form term' and 'more than one enzyme-form term' messages.
- A single `_metabolite_expr(t, role_of, macro_name)` for both `X` and `X::Inh`. It splices the role in as the type name (`EnzymeRates.$type(:X)`). That removes the duplicated ::Inh validation and the unreachable 'not declared' and 'unknown role' errors.
- `_species_expr` for bare and call-form terms. It always emits `Residual(Substrate[...], Product[...])` and does no parser-side sorting, because Species, Residual and Step sort in their constructors.
- `_walk_residual_expr`, which pushes `_metabolite_expr` results.

In both mechanism bodies:
- Build `role_of` in one expression: `Dict{Symbol,Symbol}([allo .=> :AllostericRegulator; inhibitors .=> :CompetitiveInhibitor; subs .=> :Substrate; prods .=> :Product])`. This keeps the same last-wins precedence; keep the precedence comment.
- Drop `declared_mets`.
- Keep allosteric regulators in role_of. That is behaviour-preserving; D12 decides whether to reject them.

Keep today's opaque-name rule exactly: a bare name that is not conformation-shaped is accepted only if the same name is also a call head somewhere in the block. That takes about 6 lines collecting call heads and bare names and checking them after the loop. D8 is the simpler alternative.

Also:
- `_peel_step_tag!` stays, but is called only in allow_tag mode.
- Every error message carries `macro_name`.
- `_is_conformation_shape` and every docstring of surviving code stay.
- The bodies' extra return values (source-order groups and reg sites) keep their shape for the @..._src test helpers.

Sketch: scratchpad/draft_dsl_all.jl:216-395. It has never been run, and it differs from this item in three ways:
- it moves the @allosteric_mechanism docstring to the end of the file, detached from the macro; keep it attached, because the docs build uses checkdocs=:exports;
- it drops allosteric regulators from role_of (that is D12's decision);
- it shape-checks bare names only (that is D8's decision).

**Other files and tests.** No src changes outside dsl.jl, and no test edits are required:
- the test helper macros keep their return shapes;
- all pinned messages still hold: 'opaque bound-form name', 'no enzyme-form term', 'more than one enzyme-form term', and the @allosteric_mechanism-names-itself test at test_dsl.jl:446-464.
Run the full suite, since every DSL fixture is built through this code.

**Evidence.** I read dsl.jl:308-506, 613-982 and 1164-1256. What the code shows:
- Every deleted helper is referenced only inside dsl.jl.
- In both bodies, role_of's keys equal declared_mets, and its values are exactly the four type names.
- `_StepSideTerm.conformation` equals `.sym` at all three construction sites (403, 431, 433).
- Group numbers are contiguous, so the regroup at 708-730 is the identity.

The verifier confirmed P1 at 260 with the opaque rule kept exactly (266 in the draft), and P12 at 12.

My own count:
- 308-506, 735-833 and 856-982 (425 lines) become 180;
- the remainder of _build_mechanism_expr: -41;
- body call sites, declared_mets and role_of: -22;
- keeping the opaque rule exactly: +6.
Net: about -300.

### dsl/D8 — Apply the opaque-name shape rule to every conformation label (bare and call head)
- Net lines removed (this file group): 6
- Risk: low
- Behaviour change: Newly rejected:
- call heads that are not conformation-shaped, such as `ES(P)`;
- a bare `ER` written alongside `ER(S)`, which the exemption accepts today.
No fixture is affected.
- Depends on: D7
- Approved by Denis on the sign-off page.

**Change.** In D7's `_species_expr`, check `_is_conformation_shape(conformation)` for every enzyme-form term, whether a bare name or a call head. Drop the exemption that accepts a bare name when it is a call head elsewhere, and the ~6 lines that collect call heads and bare names for it.

**Other files and tests.** None. Every call head on a step line in test/, docs/src, src and README is conformation-shaped (E 6424, Estar 80, F 32, E_c 4, G 2, Eprime 2, E_alt 2).

**Evidence.** Today `_reject_opaque_bound_forms` (456-472) exempts bare names that are call heads, and call heads themselves are never shape-checked. The verifier noted that the P1 text and its draft disagree on this rule; the draft's bare-only check is a third variant that is inconsistent (it accepts `ER(S)` but rejects a bare `ER`).

### dsl/D9 — Accept `(step) :: Tag`, the parenthesized single-step form the @allosteric_mechanism docstring shows
- Net lines removed (this file group): 0
- Risk: low
- Behaviour change: Newly accepted: a parenthesized single step followed by `:: Tag`. Today it reaches 'Expected step or step-group' (927). Unparenthesized single steps are unaffected.
- Depends on: D7
- Approved by Denis on the sign-off page.

**Change.** In allow_tag mode, accept any line whose head is `:(::)`, and take `steps = inner is a tuple ? inner.args : [inner]`. D7's `_parse_steps_block` already has this shape. In today's code the fix is the same one-condition change at 883-884.

**Other files and tests.** Add a TDD test to test_dsl.jl checking that `(E(F16BP) ⇌ E + F16BP) :: EqualAI` parses as a one-step group.

**Evidence.** A parenthesized single expression is not a tuple, so the line parses as `Expr(:(::), call, :EqualAI)` and matches none of the branches at 883, 897 or 909. grep finds this form only in the docstring, which is not a jldoctest, so CI never runs it.

### dsl/D10 — Accept a right-hand inhibitor-copy release in @enzyme_mechanism: `E(X::Inh) ⇌ E + X::Inh`
- Net lines removed (this file group): 0
- Risk: low
- Behaviour change: Newly accepted. Today plain mode peels `::Inh` off the last right-hand term and then errors 'tag annotation ... is not allowed'. The Step constructor already stores a release as the binding it reverses, so the result equals the left-to-right form.
- Depends on: D7
- Approved by Denis on the sign-off page.

**Change.** Call `_peel_step_tag!` only in allow_tag mode. D7's `_parse_steps_block` already does this; in today's code it means moving the call at 913 under `allow_tag`.

**Other files and tests.** Add a TDD test to test_dsl.jl.

**Evidence.** `_peel_step_tag!` runs on every single step whatever `allow_tag` is (913). For `E + X::Inh`, the last `+` term is a `::` Expr (976-981). No fixture writes the right-hand form.

### dsl/D11 — `_peel_step_tag!`: one code path that never peels `::Inh`
- Net lines removed (this file group): 8
- Risk: low
- Behaviour change: Errors only. An untagged allosteric step ending in `X::Inh` now reports a missing annotation. Today `:Inh` is taken as the group state and fails as an invalid state; after D14 that failure would come from the AllostericMechanism constructor.
- Depends on: D7

**Change.** Rewrite `_peel_step_tag!` (965-982):
- Take `holder, i` as `(rhs, lastindex)` when rhs is a `+` call, otherwise `(step_expr, 3)`.
- Let `t = holder.args[i]`.
- Return `nothing` unless `t` is a `::` Expr whose tag is not `:Inh`.
- Otherwise put `t.args[1]` back and return the tag.
The `tag isa Symbol` check then runs once, in the caller (D7's missing-annotation error). Update the docstring to say `::Inh` is never taken as a step tag.

**Other files and tests.** Add two tests:
- `E + X::Inh :: EqualAI` with a dual-role X still peels EqualAI (`::` is left-associative, so this parses as `(X::Inh)::EqualAI`);
- an untagged allosteric step ending in `X::Inh` reports a missing annotation.

**Evidence.** The two branches differ only in which args vector holds the trailing term. There is one caller (913). The verifier confirmed it at 9.

### dsl/D12 — Reject an allosteric-only regulator written into a catalytic step
- Net lines removed (this file group): -2
- Risk: low
- Behaviour change: Newly rejected: R inside `catalytic_steps:` when R is declared only in `allosteric_regulators:`, whether written `E + R ⇌ E(R)`, as a bound `E(R)`, or in a residual. Today this silently builds a catalytic-site AllostericRegulator, contrary to the docstring (1022-1025); no constructor checks step metabolites against the reaction.
- Depends on: D7
- Approved by Denis on the sign-off page.

**Change.** Keep allosteric regulators in role_of, so the conformation-collision check still sees them. Add one line to D7's `_metabolite_expr`: if `type === :AllostericRegulator`, error with '$macro_name: `$name` is an allosteric regulator; it binds only at its regulatory site'. Dual-role names are unaffected, because their substrate, product or inhibitor role wins in role_of. The path that emits AllostericRegulator into catalytic steps disappears.

**Other files and tests.** Add a rejection test to test_dsl.jl. A scan of all 162 @allosteric_mechanism(_src) blocks in test/ and docs/src found none that writes an allosteric-only regulator into a catalytic step.

**Evidence.** `role_of[r] = :AllostericRegulator` (1229), declared_mets includes the allosteric regulators (1221), and `_metabolite_expr` emits AllostericRegulator (828-829). The verifier found the bug real. It recommended this explicit error over P13's fix of removing allosteric regulators from role_of, which would let the name serve as a conformation label and give a misleading 'more than one enzyme-form term' error.

### dsl/D13 — Residual role check: a product added (or substrate subtracted) errors at macro expansion
- Net lines removed (this file group): -2
- Risk: low
- Behaviour change: Errors only. `residual = P - A` now gets a clear error at macro expansion instead of a run-time `MethodError: Cannot convert Product to Substrate`. No valid input changes.
- Depends on: D7

**Change.** In the Symbol branch of D7's `_walk_residual_expr`, add `get(role_of, e, nothing) === (sign_positive ? :Substrate : :Product) || error(...)`, with the message '$macro_name: a residual adds substrates and subtracts products; got `$e`'.

**Other files and tests.** Add a rejection test to test_dsl.jl.

**Evidence.** The constructor signature is `Residual(added::Vector{Substrate}, subtracted::Vector{Product})` (types.jl:37). Today the walker (479-506) checks only that names are declared, and lines 798-803 emit role-typed entries.

### dsl/D14 — Drop DSL re-checks of allosteric state tags and ligand collisions; reject duplicate allosteric_regulators names
- Net lines removed (this file group): 44
- Risk: low
- Behaviour change: - Invalid state tags are rejected by RegulatorySite/AllostericMechanism when the expanded code runs, with the constructor's message, instead of during macro expansion.
- A name listed twice in `allosteric_regulators:` is always rejected. Today it errors without an explicit site (two default sites collide) but silently keeps the last tag when it has one.
- Depends on: —
- Approved by Denis on the sign-off page.

**Change.** Delete these DSL copies of constructor checks:
- `_ALLOSTERIC_REG_STATES` (1031), a copy of `_VALID_REG_ALLO_STATES` (types.jl:867);
- `_format_state_set` (1059-1060);
- `_tagged_symbols_from_values` (1033-1057). Inline its shape check in the `allosteric_regulators:` branch: `v isa Expr && v.head == :(::) && all(a -> a isa Symbol, v.args) || error(...)`, then push `v.args[1] => v.args[2]`.
- `_build_cat_allo_states_expr` (1148-1162). Its tag check is weaker than AllostericMechanism's, which also rejects :OnlyI, and its :NonequalAI default never fires. Emit `Symbol[...]` of the group tags inline; this works with either tag format.
- the 'ligand appears in multiple regulatory sites' check (1119-1120), which AllostericMechanism already makes (types.jl:911-921).

Add `allunique(first.(allo_regs)) || error(...)`. It is needed so that no input becomes newly accepted: `A::Bogus, A::OnlyA` with an explicit site would otherwise drop :Bogus silently.

**Other files and tests.** None in src; the constructors already make these checks. Add a TDD rejection test for an allosteric_regulators name listed twice with an explicit regulatory_site. The existing tag-rejection tests (test_dsl.jl:233, 245; test_rate_eq_derivation.jl:2070-2104) use `@test_throws Exception` and still pass.

**Evidence.** The constructors already check everything the DSL re-checks:
- RegulatorySite checks the states (types.jl:134-138);
- AllostericMechanism checks catalytic tags, including :OnlyI (888-894), and ligand collisions (911-921).

The verifier only partially confirmed P4 (43-47 lines) because of the :Bogus edge case, which P17's allunique check closes.

Mapped onto the draft:
- 1030-1060: -31
- 1147-1162: -16
- inline shape check: +3
- collision check: -2
- allunique: +2
Net: -44.

### dsl/D15 — Parse `regulatory_site(...)` through `_parse_labeled_line`
- Net lines removed (this file group): 17
- Risk: low
- Behaviour change: `regulatory_site(multiplicity = 2, multiplicity = 4)` now errors instead of silently using 4. The error for a missing kwarg is reworded.
- Depends on: D2
- Approved by Denis on the sign-off page.

**Change.** Delete `_match_regulatory_site_line` (1062-1105).

In the allosteric label loop, call `_parse_labeled_line` first. Dispatch a label that is a `regulatory_site` call to `_parse_regulatory_site(label, only(values))`, which:
- requires exactly one `:kw` argument named `multiplicity` (the `:(=)` head is dead, because call arguments always parse as `:kw`);
- keeps the existing `ligands:` body loop.

Type the specs as `Tuple{Int,Vector{Symbol}}` (1129, 1169).

**Other files and tests.** None. Many fixtures use regulatory_site, so a regression would fail loudly.

**Evidence.** `_match_regulatory_site_line` re-tests the labeled-line shape (1068-1069), and its kwarg loop keeps the last value (1075-1085). It has one caller (1173). The verifier confirmed it at 18. On the draft the function goes from 45 lines to 28 and the dispatch saves 3; D2 claims 3 of that.

### dsl/D16 — _build_reg_sites_expr: build the default sites with a comprehension
- Net lines removed (this file group): 15
- Risk: low
- Behaviour change: None.
- Depends on: D14

**Change.** Rewrite `_build_reg_sites_expr` (1114-1146):
- `explicit = Set(l for (_, ligs) in reg_site_specs for l in ligs)`;
- `sites = [reg_site_specs; [(cat_n, [n]) for (n, _) in allo_regs if n ∉ explicit]]`;
- one `map(sites) do (mult, ligs)` that keeps the `haskey(tag_of, l)` declared-ligand check and emits the RegulatorySite Expr.

Site order stays as today: explicit sites first, then default sites in allo_regs order. Ligand order within a site also stays; it matters because RegulatorySite does not sort ligands and the test helpers return source-order reg sites.

**Other files and tests.** None.

**Evidence.** The function has copy loops at 1129-1136 and a separate validation loop; it has a single caller (1243). The verifier confirmed it at 13. On the draft the function goes from 33 lines to 15. The 2-line collision check belongs to D14; if D14 is declined, keep it here, which makes this item -13.

### dsl/D17 — _bare_symbols_from_values: one check whose message names no macro
- Net lines removed (this file group): 11
- Risk: low
- Behaviour change: Errors only. Today the messages hardcode '@enzyme_mechanism' and say tags are 'only valid in @allosteric_mechanism', even when the function is called from @allosteric_mechanism (1100, 1180, 1182, 1185).
- Depends on: —

**Change.** Replace the three error branches (839-855) with `all(v -> v isa Symbol, values) || error(...)`, then return `Symbol[values...]`. The message reads: '`$label:` expects bare names; got ...; atom brackets belong in @enzyme_reaction'.

**Other files and tests.** None. test_dsl.jl:179 uses `@test_throws Exception`.

**Evidence.** 8 references, all in dsl.jl. The verifier confirmed it at 11.

### dsl/D18 — Reject repeated `steps:` and `catalytic_multiplicity:` lines
- Net lines removed (this file group): -5
- Risk: low
- Behaviour change: A repeated `steps:` line (plain) or `catalytic_multiplicity:` line (allosteric) now errors instead of silently keeping the last one. No fixture repeats either label.
- Depends on: D6, D2
- Approved by Denis on the sign-off page.

**Change.** - Plain body: add `steps_block === nothing || error(...)` before `steps_block = only(values)`.
- Allosteric body: start with `cat_n = nothing`, guard with `cat_n === nothing || error(...)`, and default to `cat_n = something(cat_n, 1)` after the loop.
`catalytic_steps:` already has this guard (1200-1202). If D19 merges the bodies, one shared check covers both step labels.

**Other files and tests.** Add two TDD rejection tests to test_dsl.jl.

**Evidence.** Read 618-621 (steps_block is overwritten) and 1190-1198 (cat_n is overwritten). Only catalytic_steps has a `=== nothing ||` guard.

### dsl/D19 — Optional: merge the plain and allosteric mechanism-body parsers
- Net lines removed (this file group): 15
- Risk: low
- Behaviour change: Errors only:
- the plain macro's 'not specified' messages gain the macro prefix;
- a non-block macro argument gets a clear error instead of a FieldError on `block.args`.
- Depends on: D6, D7, D14, D15, D18
- Approved by Denis on the sign-off page.

**Change.** After D6, D7, D14, D15 and D18, replace `_parse_plain_mechanism_body` and `_parse_allosteric_mechanism_body` with one `_parse_mechanism_body(block, allosteric::Bool)`:
- label names switch on the mode: `regulators:`/`steps:` versus `catalytic_inhibitors:`/`catalytic_steps:`;
- the three allosteric-only labels are guarded by `allosteric &&`;
- the 'not specified' checks, role_of, the steps parse and the reaction emission are shared;
- it always returns three values, with `reg_sites_expr = nothing` in plain mode, to avoid an arity switch;
- it adds a `block isa Expr && block.head === :block` check, which neither macro has today.

**Other files and tests.** Edit the two test-helper call sites (test/mechanism_definitions_for_test_enzyme_derivation.jl:67, 73).

**Evidence.** After the other items, the plain body is about 37 lines and the allosteric body about 68. They share the label-loop skeleton, the checks, role_of, the steps parse and the reaction emission. The verifier partially confirmed it at 15-20 and noted the cost of an arity switch.

### dsl/D20 — Fix docstrings that are false
- Net lines removed (this file group): 0
- Risk: low
- Behaviour change: None.
- Depends on: D4, D7, D9

**Change.** Five corrections:
(a) @enzyme_mechanism docstring, 538-540: a bare `ES` is rejected as an opaque name; it does not become a conformation-only species.
(b) Line 533: `site(...)` should read `regulatory_site(...)`.
(c) @allosteric_mechanism docstring, 1013-1014: catalytic step tags are OnlyA, EqualAI and NonequalAI; :OnlyI is rejected on catalytic groups, so they are not 'the same set'.
(d) @enzyme_reaction docstring, lines 9 and 23-25: a bare `A` is accepted and takes allowed_catalytic_multiplicities. Fix this only if D4 is declined, since D4 rewrites these lines.
(e) The line 994 example works once D9 lands; if D9 is declined, rewrite it unparenthesized.
The `_StepSideTerm` 'reclassified' docstring (413-416) disappears with D7; if D7 is declined, fix it in place.

**Other files and tests.** None. docs/src api.md renders these exported docstrings.

**Evidence.** - test_dsl.jl:421-430 asserts that `ES` is rejected.
- types.jl:888-894 rejects :OnlyI on catalytic groups.
- `_parse_regulator_entries` accepts bare allosteric entries (201-202, 294), and docs/src/identify/tutorial.md:34 uses one.

### dsl/D21 — test_dsl.jl: fold four testsets that rebuild the same uni-uni mechanism
- Net lines removed (this file group): 0
- Risk: low
- Behaviour change: None.
- Depends on: —
- Approved by Denis on the sign-off page.

**Change.** Fold four testsets into one, keeping each distinct assertion once:
- test_dsl.jl:76-89 ('+ step-side syntax');
- 133-144 (m_call, which repeats 9-17);
- 528-541 (the first block of '@enzyme_mechanism');
- 599-611 ('No-atom species').
Reword the temporal comment at line 76 ('New form: + separator').

**Other files and tests.** Test-only: removes about 30 test lines with no measurable change in test time (test_dsl.jl takes 3.2 s of 1152 s). No coverage is lost, because the removed assertions are exact duplicates.

**Evidence.** The verifier read all four testsets. Each builds the same 2- or 3-step uni-uni mechanism and makes near-identical isa/n_steps/enzyme_forms assertions. Their one rate_equation call compiles a law that MECHANISM_TEST_SPECS already compiles.

## src/rate_eq_derivation.jl + src/sym_poly_for_rate_eq_derivation.jl
1) The allosteric path re-implements the plain path as _state_* twins: a second rename-map builder, a second step-param builder, column-set helpers and an A/I split of dependent assignments. It rebuilds the inactive AllostericMechanism about 6 times per constraint solve and re-runs Pass-1 solves 3-4 times.
2) Rate-body and equation-string assembly is written five times (Full/Reduced bodies in thermo, the allosteric body, three string methods).
3) Allosteric _kcat_forward (183 lines) carries hand-built k/metabolite name sets, dead branches (isempty(RS), i_state_dead, length==1 ternaries), a full I-state re-derivation, and a copy of the MWC normalization and of the regulator-K tag rule.
4) The King–Altman engine duplicates work: RE segmentation (a union-find next to _re_segment_extras), a hand-rolled Bareiss det, an R→L conversion, a numerator helper that recomputes every flux, and a post-hoc rename pass.
5) sym_poly has single-use helpers and hand-rolled monomial canonicalization.
The fix: one engine entry with defaults, one per-state view (_state_parts), one rename builder, one dependent-assignment emitter, one body builder and one text builder for both kinds, one MWC normalization, and one regulator-K rule shared by the rate law and kcat.
How the counts were measured: every item was applied in this order to a never-run draft that keeps the docstrings of surviving code (scratchpad/draft_deriv/re.jl and sp.jl). The behaviour-preserving items V1-V22 remove 605 lines (about 30%). The decision items would add -11 (V23), -50 (V27), +5/+1/+1 (V24/V25/V26). Expect ±10%, because the draft was never run.
Compile-time work also drops: fewer I-state rebuilds, fewer constraint solves, no I-state KA re-derivation in kcat, and one combined solve less per allosteric generator. This shortens allosteric test compile and _independent_param_count(am) in the split move, but the gain is unmeasured.

### deriv/V1 — Delete the _AnyMechanism alias and the duplicate concrete-mechanism default-mode forwarders
- Net lines removed (this file group): 4
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** Delete `const _AnyMechanism = AbstractEnzymeMechanism` (rate_eq:6). Write `AbstractEnzymeMechanism` in fitted_params (rate_eq:82) and rescale_parameter_values (1001). Widen the three no-mode generics, parameters(m) (40), rate_equation(m, concs, params) (546) and rate_equation_string(m) (616), to `Union{AbstractEnzymeMechanism, Mechanism, AllostericMechanism}`, and delete their concrete copies (98, 553-554, 622-623). Keep the moded concrete lifts (95-97, 550-552, 619-621) and fitted_params(::Union{Mechanism, AllostericMechanism}).

**Other files and tests.** None. The alias is used only in this file; no tests or docs reference it.

**Evidence.** grep _AnyMechanism across src/test/docs finds 6 hits, all in rate_eq. Mechanism and AllostericMechanism are not subtypes of AbstractEnzymeMechanism, so the Union methods cannot be ambiguous with the 4-argument generated methods. types.jl defines a different `_AnyMech` union, which makes the alias easy to confuse. Measured -4 on the draft.

### deriv/V2 — _poly_to_expr: drop the unused conc_syms argument, the unreachable tail and _nest_binary's redundant n==2 line
- Net lines removed (this file group): 4
- Risk: low
- Behaviour change: none (identical Expr trees)
- Depends on: —

**Change.** sym_poly:113: change the signature to `_poly_to_expr(p::POLY, param_syms::Set{Symbol} = Set{Symbol}())`; conc_syms is never read. Collapse the return tail (143-148) to `isempty(neg) && return _nest_binary(:+, pos); ne = _nest_binary(:+, neg); isempty(pos) ? :(-$ne) : :($(_nest_binary(:+, pos)) - $ne)`. The final `return 0` is unreachable after the isempty(p) early return. Delete `_nest_binary`'s `n == 2` line (161): mid = 1 builds the same 2-operand Expr. Callers drop the third argument, and the kcat calls become `_poly_to_expr(x)`. Most call sites are rewritten by later items anyway.

**Other files and tests.** Test: allosteric_ground_truth.jl:1204 drops its third argument.

**Evidence.** conc_syms appears only in the signature (sym_poly:113). There are 17 src callers and 1 test caller. The verifier confirmed it. Measured -4 in sym_poly on the draft.

### deriv/V3 — _mentions(ex, s) replaces _expr_references_any; drop the redundant first line of _expr_to_string
- Net lines removed (this file group): 8
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** Replace sym_poly:258-266 with the one-liner `_mentions(ex, s::Symbol) = ex === s || (ex isa Expr && any(a -> _mentions(a, s), ex.args))`; every caller passes Set([:Keq]). Delete sym_poly:186 (`x isa Union{Number, Symbol} && return string(x)`), which line 187 already covers.

**Other files and tests.** Test: test_rate_eq_derivation.jl:816 becomes `EnzymeRates._mentions(expr, :Keq)`. The two keq_set call sites in rate_eq (652-654, 1676-1682) are rewritten by V18.

**Evidence.** All three callers (rate_eq:654, 1682, test:816) pass Set([:Keq]). A Number or Symbol is not an Expr, so line 187 already returns string(x). Measured -8 in sym_poly on the draft.

### deriv/V4 — build_power_expr: ternary _pa and one-line product assembly (emitted Expr unchanged)
- Net lines removed (this file group): 18
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** Replace the nested if-chain of `_pa` (sym_poly:226-240) with the ternary `exp == 1 ? sym : exp == -1 ? :(1 / $sym) : !isinteger(exp) ? :($sym ^ $(Float64(exp))) : exp > 0 ? :($sym ^ $(Int(exp))) : :(1 / $sym ^ $(Int(-exp)))`. Replace 249-255 with `isempty(terms) ? 1 : length(terms) == 1 ? only(terms) : Expr(:call, :*, terms...)`. Keep the varargs `*`; the switch to _nest_binary is in dropped.

**Other files and tests.** None. The build_power_expr return-type tests (test_rate_eq_derivation.jl:1751-1779) still hold.

**Evidence.** isinteger(::Rational) holds exactly when denominator == 1. `:(1)` is the Int 1. Callers never pass exponent 0 (lines 242 and 247 filter it). The single src caller is thermo:697. Measured -18 in sym_poly on the draft.

### deriv/V5 — Monomial helpers: fold _mono_op into _mono_mul, shift by poly_mul in _reduce_conc_lowest_terms, drop _re_weight_ratio's inverse keyword
- Net lines removed (this file group): 14
- Risk: low
- Behaviour change: none (identical canonical POLYs)
- Depends on: —

**Change.** sym_poly:39-46: `_mono_mul(a, b)` keeps the Dict accumulate/filter/sort body but loses the always-1 sign argument and the wrapper line. Keep two plain loops, with no varargs splat inside poly_mul's inner loop. sym_poly:74-86: replace the hand-written `shift` closure with `shift = POLY(_mono((s => -mn for (s, mn) in mins)...) => 1); poly_mul(num, shift), poly_mul(den, shift), poly_mul(weight, shift)`, with a comment that multiplying by a monomial is injective. rate_eq:222-238: drop the `inverse` keyword of `_re_weight_ratio` (and its sgn variable); its single inverse caller in `_compute_alpha` (270-271) uses `_invert_monomial(_re_weight_ratio(s, K))`.

**Other files and tests.** None. Tests: none (the `_mono` varargs signature that tests call is unchanged).

**Evidence.** `_mono_op` has one caller, _mono_mul, always with sign 1. shift() multiplies every monomial by a fixed monomial and keeps the coefficients, which equals poly_mul by that monomial. The inverse=true caller is unique (rate_eq:270-271). Following the verifier of #19, the `_mono(a..., b...)` splat form is avoided in the sym_det hot path. Measured sym_poly -11 and rate_eq -3.

### deriv/V6 — One RE-segment routine: the derivation reads _indexed_re_segments; delete _compute_re_groups and _enumerate_species
- Net lines removed (this file group): 54
- Risk: medium
- Behaviour change: none. Species order and segment order are identical, since both routines number segments by their lowest-index form. Only member order inside a segment changes (ascending to BFS), and every consumer ignores it: _segment_root is an argmin over a total key, sigma is a POLY sum, D and alpha are indexed per segment or form, and _poly_to_expr sorts by a total monomial key.
- Depends on: —

**Change.** Delete `_compute_re_groups` (rate_eq:149-181, a second union-find) and `_enumerate_species` (183-193). In `_raw_symbolic_rate_polys`, write `enz_species, groups, _, idx, form_to_group = _indexed_re_segments(steps(mech))` and index forms by Species (`idx[from_species(s)]`). This replaces the name-keyed `enz_name_to_form` Dict (321-322) and the paired i_f/j_f lookups. `_compute_alpha` and the numerator take `idx`. `_eq_complexity` (V7) uses the same call.

**Other files and tests.** mechanism_enumeration.jl:1898 (`_re_segment_count`) and 2068 (`_group_re_segments`) call `_compute_re_groups`. Either switch them to `_indexed_re_segments(steps(cm))` (enum group, about 0 net lines), or keep a 3-line `_compute_re_groups` shim here (which costs 3 of these lines). Consider moving `_indexed_re_segments` (mechanism_enumeration.jl:1396) next to `_re_segment_extras` in types.jl. thermo's `_enumerate_species_names` can index the same way (thermo group). Tests: none; both deleted functions have 0 test refs.

**Evidence.** I read types.jl:725-768. Species are walked in the same from-then-to, first-seen order. Each BFS root is the smallest unvisited index, and the union-find creates a group at each first-seen root, so segment order matches. The guards are the golden file, exact-string tests and Expr-shape tests. The enumeration-side `_re_segment_count_after_flip` swap is left to the enum group (it needs a timing). Measured -54 on the draft.

### deriv/V7 — _eq_complexity: build the segment Laplacian directly and use LinearAlgebra.det; delete _segment_graph_terms, _spanning_tree_count and _bareiss_det
- Net lines removed (this file group): 41
- Risk: low
- Behaviour change: none (exact BigInt determinant; same V×τ value and type; a disconnected graph still gives 0)
- Depends on: V6
- Approved by Denis on the sign-off page.

**Change.** Replace rate_eq:399-453 with one `_eq_complexity(mech::Mechanism)` of about 13 lines: `_indexed_re_segments`, a `G <= 1 && return big(G)` early return, a BigInt Laplacian over the SS steps that skips self-loops, and `G * det(lap[2:end, 2:end])`. Keep `_eq_complexity(m::AllostericMechanism) = _eq_complexity(_state_mechanism(m, :A))`. Fold the Matrix–Tree sentence into the kept docstring (385-398). Add `using LinearAlgebra: det` to EnzymeRates.jl.

**Other files and tests.** EnzymeRates.jl +1 line. LinearAlgebra is already in [deps] but no src file loads it. Test: delete the `_bareiss_det exact integer determinant` testset (test_rate_eq_derivation.jl:1946-1955, 10 lines). The `_eq_complexity` testsets keep the V×τ coverage: test_rate_eq_derivation.jl:1881-1944, test_identify_rate_equation.jl:1040 and 1065, and test_mechanism_enumeration.jl:9382.

**Evidence.** The verifier checked the installed 1.12 stdlib: generic.jl:1834 has `det(A::AbstractMatrix{BigInt}) = det_bareiss(A)`. The effective Julia floor is 1.11, set by the Dates/Statistics compat entries. A direct Laplacian that skips g1 == g2 equals today's adjacency→Laplacian path. `_segment_graph_terms` and `_spanning_tree_count` each have a single caller. Measured -41 on the draft.

### deriv/V8 — King–Altman assembly: build the Laplacian directly, keep each SS step's fluxes, inline _compute_numerator
- Net lines removed (this file group): 19
- Risk: low
- Behaviour change: none (same Laplacian entries; self-loops still enter the numerator only)
- Depends on: V6

**Change.** In `_raw_symbolic_rate_polys` (rate_eq:331-353), replace the R matrix, the R→L conversion and the D special case with one pass over the SS steps. Each step computes fwd and rev once, pushes `(s, a, b, fwd, rev)` to `ss`, and accumulates L directly: `L[g1,g1]+=fwd; L[g1,g2]-=fwd; L[g2,g2]+=rev; L[g2,g1]-=rev`, skipping g1 == g2. Then `D = [sym_det(L[o, o]) for g in 1:G for o in (setdiff(1:G, g),)]`. Inline `_compute_numerator` (470-502) as a loop over `ss`, with its docstring kept as a comment. `_ss_contrib` takes metabolites and applies `name.`. sym_poly: delete `poly_neg` (28) and `poly_const` (18), which have one use each, and let `sym_det(M)` read `n = size(M, 1)`. Keep the `i_free === nothing ? poly_one()` fallback (see dropped).

**Other files and tests.** None. Tests: none (test_rate_eq_derivation.jl:1304-1305 and 2576-2579 cover d_free).

**Evidence.** R[i,i] is filled but never read, because L[i,i] sums only j≠i (345-349). _compute_numerator recomputes the same kf/kr poly_sym, _ss_contrib values and form indices (486-498 vs 331-343). sym_det returns poly_one() at n == 0, so a 0×0 minor needs no special case. Measured rate_eq -18 and sym_poly -1.

### deriv/V9 — Apply the Pass-2 Wegscheider rename at the RE leaves; delete _rename_symbols
- Net lines removed (this file group): 31
- Risk: medium
- Behaviour change: none. Rename keys and targets are RE binding Kd names, and these enter the polys only through step_to_K and alpha. Substitution is a ring homomorphism. The L entries are single-signed, so no entry can vanish. _compute_alpha's cycle check compares concentration parts only.
- Depends on: —

**Change.** Rename while building step_to_K (rate_eq:324-326): `K = name(step_params[i][1], mech); get(rename_map, K, K)`. Delete the four post-passes: `_rename_symbols` on D[g_free] (360), on csigma (365) and on num/den (371-372). Delete `_rename_symbols` with its section header (sym_poly:268-294); its docstring is stale, since no A→I renaming happens. Update the engine docstring to say rename_map is applied at the leaves. Pass 2 (invariant 4) is kept: the same map is built and applied, only earlier.

**Other files and tests.** None. Tests: the exact-string tests over the non-competitive-inhibitor and non-essential-activator specs (where the rename is non-empty) and test_mechanism_enumeration.jl:9255 and 9336 must pass unchanged.

**Evidence.** _rename_symbols has callers only at rate_eq:360, 365, 371 and 372 (0 in tests and docs). The rename keys come only from binding_set (rate_eq:132, 1228). The csigma rename is already redundant, because the whole den is renamed again at 372. Measured rate_eq -3 and sym_poly -28.

### deriv/V10 — _concentration_symbols as one expression
- Net lines removed (this file group): 4
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** Replace rate_eq:202-208 with `_concentration_symbols(mech) = (rxn = reaction(mech); Set{Symbol}([name.(substrates(rxn)); name.(products(rxn)); [name(regulator(rm)) for rm in regulators(rxn)]]))` and keep the docstring. If the types group adds a shared `_metabolite_names(rxn)` (used by metabolites(em) and IdentifyRateEquationProblem), use `Set(_metabolite_names(reaction(mech)))` instead.

**Other files and tests.** Optional types.jl helper `_metabolite_names(rxn)` (types group, from proposals #84/#103/#1). Tests: test_mechanism_enumeration.jl:96 and 131 call it unchanged.

**Evidence.** The same substrates/products/regulators walk appears in types.jl:1574-1593 and identify_rate_equation.jl:45-52. Measured -4 on the draft.

### deriv/V11 — One inactive-conformation graph (_inactive_groups) for _state_allo_mechanism and the copy rule; _reachable_from_free seeds directly
- Net lines removed (this file group): 2
- Risk: low
- Behaviour change: none. `_reachable_from_free` runs to a fixpoint where no step has exactly one endpoint in reach, so 'neither endpoint stranded' holds exactly when 'both endpoints in reach', which is the enum rule. Group order and tags are preserved, and the constructor still canonicalizes.
- Depends on: —

**Change.** Add `_inactive_groups(am)`, a list aligned with steps(am): `:OnlyA` groups are emptied, and every other group keeps only the steps whose two forms are in `_reachable_from_free` of the non-OnlyA groups. `_state_allo_mechanism(:I)` (1101-1121) becomes `groups = _inactive_groups(am); keep = findall(!isempty, groups); AllostericMechanism(reaction(am), groups[keep], cat_allo_states(am)[keep], catalytic_multiplicity(am), regulatory_sites(am))`, which drops the all_forms/stranded detour. `_reachable_from_free` (1055-1062) builds its root set in one comprehension instead of an O(N²) `forms` vector; its fixpoint loop is unchanged.

**Other files and tests.** mechanism_enumeration.jl:2346-2351 (`_redundant_copy_groups(am)`) becomes `inactive = _inactive_groups(am)` (-5 enum lines). Without the shared helper, this file alone would save about 12, but the inactive graph would keep two definitions that must agree.

**Evidence.** I read rate_eq:1055-1121 and mechanism_enumeration.jl:2342-2351. Both drop the :OnlyA groups and then prune by reachability from free E. This path is codegen-time only (the I-state polys). Measured rate_eq -2 (-12 of simplification, +10 for the documented shared helper).

### deriv/V12 — One parameters(Full) and one parameters(Reduced) method via a _concrete lift; _ss_rate_constant_names filters the full enumeration
- Net lines removed (this file group): 38
- Risk: low
- Behaviour change: none (same symbols in the same order; same SS set)
- Depends on: —

**Change.** Add `_concrete(m::EnzymeMechanism) = Mechanism(m)` and `_concrete(m::AllostericEnzymeMechanism) = AllostericMechanism(m)`. Rename `_enumerate_parameters_full_allosteric` to the method `_enumerate_parameters_full(am::AllostericMechanism)`. Replace the four @generated parameters methods (42-79) with one Full method over AbstractEnzymeMechanism, `(Tuple(name(p, m) for p in _enumerate_parameters_full(m))..., :E_total)`, keeping the collapse-mirror comment, and one Reduced method (the two Reduced bodies are identical). Replace both `_ss_rate_constant_names` methods (695-716) with one filter `p isa Union{Kon, Koff, Kfor, Krev}` over `_enumerate_parameters_full(_concrete(em))`. `_dependent_param_exprs(::Type{<:AllostericEnzymeMechanism})` (1415-1416) lifts via `_concrete`.

**Other files and tests.** thermo:313 `_dependent_param_exprs(M::Type{<:EnzymeMechanism})` can merge into one `Type{<:AbstractEnzymeMechanism}` method via `_concrete` (thermo -1). If the types group adds the same lift (#88), define it once there; `_concrete` is unused today. Tests: none. The golden PARAMS_FULL/PARAMS_REDUCED lines and the `_ss_rate_constant_names` testset (test_rate_eq_derivation.jl:1713-1744) pin the outputs.

**Evidence.** The full allosteric enumeration is the A part, then _all_i_state_parameters, then Kreg A, then Lallo. Filtering it by SS type drops the RE reps (Kd/Kiso), Kreg and Lallo, which leaves exactly today's allosteric _ss_rate_constant_names set. `_enumerate_parameters_full_allosteric` has a single caller (69). Measured -38 on the draft.

### deriv/V13 — Two allosteric parameter walkers: _cat_params(am, state) and _kreg_params(am, state)
- Net lines removed (this file group): 13
- Risk: low
- Behaviour change: none (same Parameters, same order)
- Depends on: V12

**Change.** Add `_cat_params(am, state)`, which emits group reps: in :A, an :EqualAI group is tagged :EqualAI and every other group :A; in :I, :OnlyA groups are skipped and the rest tagged :I. Add `_kreg_params(am, state)`, which skips :OnlyI ligands in :A and :OnlyA ligands in :I. Then `_all_i_state_parameters(am) = [_cat_params(am, :I); _kreg_params(am, :I)]` (1247-1261) and `_enumerate_parameters_full(am) = [_cat_params(am, :A); _all_i_state_parameters(am); _kreg_params(am, :A); Lallo()]` (1287-1304; keep its order docstring). In `_dependent_param_exprs(am)` (1399-1410), `reg_params_a = [name(p, am) for p in _kreg_params(am, :A)]`. The I-side loop stays, because :EqualAI ligands become dep mirrors.

**Other files and tests.** Tests: the `_all_i_state_parameters` testset (test_rate_eq_derivation.jl:2296-2330) is unchanged. The types.jl lookup walkers (_onlyA_parameters_for_sym/_all_params_for_sym) belong to the types group's dead-code deletion.

**Evidence.** The tag rules match today's code: the I walker tags every non-OnlyA group :I (EqualAI included); the A walker maps EqualAI to :EqualAI; Kreg A skips OnlyI and Kreg I skips OnlyA. The site/ligand order is the same. Measured -13 on the draft.

### deriv/V14 — One Wegscheider rename-map builder; delete _state_wegscheider_rename_map, the Type/instance wrappers and the dead chain-compression loop
- Net lines removed (this file group): 47
- Risk: medium
- Behaviour change: none
- Depends on: —

**Change.** Rewrite `_build_wegscheider_rename_map(mech::Mechanism; step_params = _step_parameters(mech), all_params = _raw_param_symbols(mech))` with three parts: binding_set as a comprehension; one kernel call with an empty rename that forwards step_params and all_params; and `Dict{Symbol,Symbol}(lhs => rhs for (lhs, rhs) in dep_raw if rhs isa Symbol && lhs in binding_set && rhs in binding_set)`. Delete `_state_wegscheider_rename_map` (1202-1236); the allosteric callers pass `step_params = sp, all_params = _state_all_params(cm, sp)`. Delete the ::Type and ::EnzymeMechanism wrappers (142-145); the src caller at 381 already holds mech. Merge the two docstrings, adding why no tie can chain onto another (full Gauss–Jordan). Keep all_params explicit at the allosteric call sites, so the non-allosteric column list is untouched.

**Other files and tests.** test_rate_eq_derivation.jl:619 passes `mech` (already in hand at 611) instead of `m`. test_types.jl:1307-1316 ('wegscheider rename map on Mechanism') compares the deleted EnzymeMechanism wrapper with the Mechanism method. Rewrite it as `_build_wegscheider_rename_map(Mechanism(compile_mechanism(m))) == _build_wegscheider_rename_map(m)`, which keeps the compile round-trip coverage with no deletion. thermo: the stale `_build_kinetic_rename_map` docstring sits on the kernel Type wrapper that V18 deletes. If the thermo group drops the all_params kwarg (#46), the allosteric call sites shrink further.

**Evidence.** rate_eq:116-140 and 1212-1236 are identical line for line apart from step_params and all_params. _solve_dependent_set (thermo:679-686) eliminates each pivot column from every other row, so a dependent's RHS uses only independent columns. Therefore `get(rename, rhs, rhs) == rhs`, and the chain loops (133-137, 1229-1233) never fire. Pass-2 absorption (invariant 4) is kept, in one copy. Measured -47 on the draft.

### deriv/V15 — _state_parts(am, state) builds each conformation's (Mechanism, step_params) from one _state_allo_mechanism; the engine takes defaults and the V×τ guard
- Net lines removed (this file group): 20
- Risk: low
- Behaviour change: none for derivable mechanisms. Direct 3-argument engine calls now also enforce MAX_RATE_EQUATION_TERMS, and only tests make such calls.
- Depends on: V14

**Change.** Replace `_state_step_params` (1123-1152) with `_state_parts(am, state) -> (cm, sp)`. It makes one `_state_allo_mechanism` call, sets `cm = Mechanism(reaction(sam), steps(sam))`, and builds `sp = Vector{Parameter}[_emit_cat_params_for_rep(s, tag(g)) for (g, group) in enumerate(steps(sam)) for s in group]`; the 4-way Kd/Kiso/Kon/Koff/Kfor/Krev switch is `_emit_cat_params_for_rep`. Keep `_state_mechanism` (used by enum, `_eq_complexity` and tests) and the `@assert` alignment check. Write `_state_all_params(cm, sp) = unique(name(p, cm) for group in sp for p in group)`. Give the engine defaults, `_raw_symbolic_rate_polys(mech, step_params = _step_parameters(mech), rename_map = _build_wegscheider_rename_map(mech; step_params))`, with `_assert_derivable(mech)` as its first line. The Type method becomes one line, and `_state_rate_polys` becomes `cm, sp = _state_parts(am, state)` plus the engine call.

**Other files and tests.** thermo `_step_parameters` (thermo:36-46) can also use `_emit_cat_params_for_rep` (thermo group, #48/#90). Tests: none. The direct 3-argument engine calls (test_mechanism_enumeration.jl:9255, 9334) are unchanged and now also pass through the V×τ guard.

**Evidence.** Every `_state_step_params` call is paired with `_state_mechanism` on the same state (864, 1178-1180, 1213-1214, 1323-1324, 1496-1497, 1555-1556). Each pair rebuilds the inactive AllostericMechanism twice, and every rebuild re-runs canonicalization plus `_onlya_haldane_violation`. Lines 1145-1148 match `_emit_cat_params_for_rep` (types.jl:1937-1945) token for token. Measured -20 on the draft.

### deriv/V16 — Allosteric _dependent_param_exprs: inline _combined_state_dependent_exprs, build each state's parts, rename and system once, stack by block assignment
- Net lines removed (this file group): 26
- Risk: medium
- Behaviour change: none (same columns, priorities, rows and solve)
- Depends on: V14, V15

**Change.** Fold `_combined_state_dependent_exprs` (1308-1356) into `_dependent_param_exprs(am)`: `systems = map((:A, :I)) do state; cm, sp = _state_parts(am, state); cols = _state_all_params(cm, sp); rename = _build_wegscheider_rename_map(cm; step_params = sp, all_params = cols); rename, _assemble_constraints(cm, rename; step_params = sp, all_params = cols, is_i_state = state === :I) end`. Then set `columns = unique([cols_A; cols_I])` and `priority = merge(Dict(zip(cols_I, pri_I)), Dict(zip(cols_A, pri_A)))` (A wins on shared :EqualAI columns), fill a zero matrix with two block assignments, and make one `_solve_dependent_set` call. Reuse `merge(rename_A, rename_I)` for the indep filter (1384-1385) instead of re-deriving both maps. Merge the two docstrings, keeping the I-above-A pivot rationale.

**Other files and tests.** The comment at test_rate_eq_derivation.jl:1636 names the deleted helper and needs rewording. Recommended, to cover the touched branch: add one `@allosteric_mechanism` fixture with a grouped RE binding square, so the `!isempty(rename)` reference guard (1386-1394) runs; the #54 verifier found no test that reaches it.

**Evidence.** Today one call builds the inactive AllostericMechanism about 6 times and runs 4 Pass-1 solves plus the combined solve; afterwards it is 1 rebuild and 2 + 1 solves. The call runs once per split candidate via `_independent_param_count(am)` (mechanism_enumeration.jl:2058-2062). unique() reproduces the A-then-I-only column order, and block assignment equals the nonzero-only copy loops. Measured -26 on the draft.

### deriv/V17 — Delete the A/I split of dependent assignments: one _dep_assignments(dep)
- Net lines removed (this file group): 51
- Risk: low
- Behaviour change: Inside the generated allosteric rate_equation/_kcat_forward bodies, the assignment statements come out in one global name order instead of A-block then I-block. Values, allocations and rate_equation_string (which re-sorts its lines) are unchanged.
- Depends on: —

**Change.** Delete `_i_state_symbol_set` (1488-1506) and `_build_dep_assignments` (1508-1537). Add `_dep_assignments(dep) = [Expr(:(=), sym, rhs) for (sym, rhs) in sort!(collect(dep); by = first)]` with a docstring stating the full Gauss–Jordan argument. Allosteric `_kcat_forward` uses `dep, indep = _dependent_param_exprs(M_type)` and `_dep_assignments(dep)...`, which deletes the alias lines 886-891. The rate bodies and strings switch over in V18, and the non-allosteric kcat in V23.

**Other files and tests.** The same sort-and-assign comprehension at thermo:769-770 disappears with `_build_rate_body` in V18. Tests: none (0 refs to either helper).

**Evidence.** _solve_dependent_set is a full Gauss–Jordan elimination, so no RHS reads another dependent. The :EqualAI regulator mirror reads K_A_reg, which is in indep (1403, 1411). Every consumer splices a..., i... back to back. Each `_build_dep_assignments` call also re-ran the combined solve plus 4 state rebuilds. Measured -51 on the draft.

### deriv/V18 — One rate-body builder and one equation-text builder for both mechanism kinds and modes
- Net lines removed (this file group): 91
- Risk: medium
- Behaviour change: none. The arithmetic Exprs come from the same `_poly_to_expr`/`_nest_binary` calls, and the non-allosteric body gains an explicit `return`. Non-allosteric constraint RHS are now rendered by `_expr_to_string`, which prints every shape build_power_expr emits exactly as Base `string` does: Symbol, 1/x, x^n, 1/x^n, x^f, and a flat `*` with parenthesized `/` factors. The jldoctest line is unchanged. Sorting by LHS equals the allosteric whole-line sort, because a space sorts below every identifier character.
- Depends on: V1, V2, V3, V12, V17

**Change.** Add `_num_den_exprs(M::Type{<:EnzymeMechanism})`, which builds the num/den Exprs from `_raw_symbolic_rate_polys` and `_poly_to_expr`, and rename `_allosteric_num_den_exprs` to its `Type{<:AllostericEnzymeMechanism}` method. Add `_rate_body(M, param_syms, dep)`, a block of the params destructure, the `metabolites(M())` destructure, `_dep_assignments(dep)...` and `return E_total * (num) / (den)`. Two @generated rate_equation methods: FullMode for EnzymeMechanism (`(_raw_param_symbols(M())..., :E_total)`, empty dep) and ReducedMode for `M <: AbstractEnzymeMechanism` (`dep, indep = _dependent_param_exprs(M)`). They replace 556-567 and the allosteric 1654-1661. Add `_equation_text(M, param_syms, dep; annotate = false)`: the destructure lines; the dep lines sorted by name, with RHS rendered by `_expr_to_string` and split by `_mentions(rhs, :Keq)` into Wegscheider/Haldane sections; and the v line from `_num_den_exprs`. Keep the load-bearing-sort note in its docstring. Three one-line string methods: Full; non-allosteric Reduced (`first(_dependent_param_exprs_kernel(Mechanism(m), Dict()))`, annotate = true, with a comment that the Pass-1 solve keeps absorbed ties visible); allosteric Reduced (`first(_dependent_param_exprs(typeof(m)))`). Delete `_raw_rate_expr_and_symbols` (504-520), `_rate_v_line`, `_partition_constraint_lines!`, `_append_constraint_sections!`, the old string methods (625-682), `_build_allosteric_rate_body` (1630-1652) and the allosteric string (1663-1702). Reword the ANNOTATION_SUBSTITUTED docstring, since one display site remains.

**Other files and tests.** thermo: delete `_sorted_raw_param_symbols` and both `_build_rate_body` methods (thermo:747-776), plus the `_dependent_param_exprs_kernel(::Type)` wrapper with its stale `_build_kinetic_rename_map` docstring (thermo:731-738). That is about -38 thermo lines, which the thermo group must not count again (#50, #57, #36, #52). Test: test_rate_eq_derivation.jl:1675 walks both Exprs of `EnzymeRates._num_den_exprs(typeof(spec.mechanism))` instead of `_raw_rate_expr_and_symbols`.

**Evidence.** Today the E_total*(num)/(den) body is built twice, the v line twice, the Haldane/Wegscheider partition twice and the Reduced tuple 4 times. rate_equation_string runs per candidate at identify time (identify_rate_equation.jl:542), and after this item does 2 Pass-1 solves instead of 3 plus a Pass-2 solve. This touches invariant 2: run test_rate_equation_performance, the Expr-shape and flat-string tests, the jldoctest and test_allosteric_golden.jl. Measured -91 in this group.

### deriv/V19 — One MWC free-enzyme normalization (_mwc_state_polys) for the rate law and kcat
- Net lines removed (this file group): 16
- Risk: low
- Behaviour change: none (both consumers receive identical polys and Exprs)
- Depends on: V2, V15

**Change.** Add `_mwc_state_polys(am, mets) -> (num_A, den_A, num_I, den_I, d_A, d_I)`. It derives both states (`_state_rate_polys`) and then normalizes: with equal weights, both d's become poly_one(); with two metabolite-free monomials, each state's polys are divided by that state's inverse monomial and both d's become poly_one(); otherwise the raw d's are returned for cross-weighting. `_num_den_exprs(::Type{AEM})` renders `D_X_expr = _poly_to_expr(d_X, cat_params)`; poly_one() renders as 1, so `_mwc_cross_weight` stays a no-op. Allosteric `_kcat_forward` groups `poly_mul(num_A, d_I)`, `poly_mul(den_A, d_I)` and the I-state counterparts. Move the two explanatory comment blocks (832-841, 1559-1571) into the helper's docstring.

**Other files and tests.** None. Tests: none (allosteric_ground_truth d_free checks, the golden file and the analytical-kcat specs guard it).

**Evidence.** kcat (842-855) and the num/den Exprs (1574-1585) use the same three-way branch with the same `_is_metabolite_free_monomial` test, and kcat's comment says it must match rate_equation. This area has had real bugs (the MWC ping-pong d_free fix 43f4d23). Measured -16 on the draft.

### deriv/V20 — kcat grouping splits monomials by metabolite set; delete the hand-built a/i_param_names
- Net lines removed (this file group): 33
- Risk: low
- Behaviour change: none
- Depends on: V19

**Change.** Change `_kcat_groups_from_polys(num, den, mets::Set{Symbol})` (718-753) to one local `group(p, forward_only)` that keys by `filter(q -> q.first in mets, mono)` and accumulates the `∉ mets` part. The filters keep MONO sorted, and Keq/E_total never occur in a POLY. Allosteric kcat passes `cat_mets` and deletes `a_param_names`/`i_param_names` (857-869); non-allosteric kcat passes `Set(metabolites(M()))`.

**Other files and tests.** None. Tests: none (test_analytical_kcat over MECHANISM_TEST_SPECS and test_kcat_rescaling guard it).

**Evidence.** a_param_names = cols ∪ (poly symbols \ cat_mets), so classifying a symbol as k exactly when it is ∉ cat_mets reproduces today's split. Non-allosteric polys hold only raw parameter symbols (the Wegscheider targets are raw binding K's) and metabolites. Measured -33 on the draft.

### deriv/V21 — Allosteric _kcat_forward: delete dead branches and the redundant I-state re-derivation; delete _i_state_num_zero
- Net lines removed (this file group): 46
- Risk: low
- Behaviour change: kcat values unchanged. One Expr-shape change: when the I state is live but lacks the pattern in its numerator, A_I is the literal 0 instead of `0 * y^(n-1)`, which has the same value for finite y.
- Depends on: V19, V20

**Change.** (a) `_mwc_power_pair(x, y, n) = (n == 1 || x == 0 ? x : …)`. A missing I-numerator pattern or a dead I cycle then gives A_I = 0 without `i_state_dead = _i_state_num_zero(am)` (894), which re-ran the whole inactive King–Altman derivation. Collapse 908-925 to two `_mwc_power_pair` calls, the I one guarded by `haskey(den_I_groups, k) ? … : (0, 0)`. (b) Delete the `if isempty(RS)` branch (927-930): the mask loop over 0:0 emits the same Expr. Keep the `isempty(W_A_factors) && isempty(W_I_factors)` special case. (c) Drop the `length(x) == 1 ? x[1] : _nest_binary(...)` ternaries (955-956, 960-961, 969-974); `_nest_binary` returns a lone term itself. (d) `all_ligs = allosteric_regulators(am)` replaces 897-902. (e) Delete `_i_state_num_zero` (1031-1041) and keep its dead-cycle explanation in the num/den comment at 1619-1622.

**Other files and tests.** Tests: test_rate_eq_derivation.jl:1154 and test_mechanism_enumeration.jl:7374 become `isempty(first(EnzymeRates._state_rate_polys(am, :I)))`. Reword the comment mentions at mechanism_definitions_for_test_enzyme_derivation.jl:2411 and test_mechanism_enumeration.jl:7360-7361. types.jl `allosteric_regulators(::AllostericMechanism)` gains its first src caller, so the types/dead-code group must not delete it; otherwise inline `unique(lig for site in regulatory_sites(am) for lig in ligands(site))`.

**Evidence.** Multiplying by a nonzero monomial or POLY preserves emptiness, so isempty(num_I) == _i_state_num_zero(am). RS empty means n_ligs = 0, which gives one mask with no site factors. _nest_binary with n == 1 returns terms[1] (sym_poly:160). The constructor rejects repeated ligand names, so the all_ligs dedup never fires. Measured -46 on the draft.

### deriv/V22 — One regulator-K rule (_reg_K) for the rate law and kcat corners; one MWC term builder
- Net lines removed (this file group): 25
- Risk: low
- Behaviour change: none (identical Exprs with the same factor order and _nest_binary trees)
- Depends on: V21

**Change.** Add `_reg_K(am, site, lig, tag, inactive)`, which returns `nothing` when `tag === (inactive ? :OnlyA : :OnlyI)` and otherwise `name(Kreg(site, lig, inactive && tag !== :EqualAI ? :I : :A), am)`; its docstring carries the comment that :EqualAI ligands share the A name. `_reg_site_expr(am, site, inactive)` (1427-1443) takes the site and pushes `lig / K` for each non-nothing K. The kcat corners (938-952) loop over `((sat_terms_A, false), (sat_terms_I, true))` with `_reg_K`. In `_num_den_exprs(::Type{AEM})`, precompute `reg_A`/`reg_I = Any[_power_expr(_reg_site_expr(am, s, inactive), multiplicity(s)) for s in RS]`, and replace make_num_term/make_den_term (1595-1612) with two one-liners that keep the factor order `[N, Q^(CatN-1)?, reg...]` and `[Q^CatN, reg...]`. Delete `_power_expr`'s unreachable `n == 0` line and reword its docstring.

**Other files and tests.** None. Tests: none (the golden file byte-checks the rate strings).

**Evidence.** Case check: in the active state NonequalAI/EqualAI/OnlyA map to :A and OnlyI is skipped; in the inactive state EqualAI maps to :A, NonequalAI/OnlyI to :I and OnlyA is skipped. _reg_site_expr (1431-1438) and the corners (940-951) agree on every case. Every _power_expr caller has n ≥ 1. Measured -25 on the draft.

### deriv/V23 — Non-allosteric _kcat_forward shares pattern selection (_kcat_keys) and fails loudly when no candidate exists
- Net lines removed (this file group): 11
- Risk: low
- Behaviour change: A non-allosteric mechanism with no product-free saturating pattern now errors at generation (the first _kcat_forward call) with 'no kcat components…' instead of `MethodError: max()`. The allosteric message loses its 'AllostericEnzymeMechanism produced' prefix unless a label is passed. No values change otherwise.
- Depends on: V17, V20
- Approved by Denis on the sign-off page.

**Change.** Add `_kcat_keys(num_groups, den_groups, prods)`, which returns the sorted keys present in both maps with no product, errors 'no kcat components…' when none remain, and carries the products = 0 comment. Add `_max_expr(cands)`. The non-allosteric generator (771-806) becomes: polys, then `_kcat_groups_from_polys(…, Set(metabolites(M())))`, then candidates over `_kcat_keys`, then a block of the params destructure `(indep..., :Keq)`, `_dep_assignments(dep)...` and `return _max_expr(candidates)`. Allosteric kcat uses `_kcat_keys` for a_keys (878-883) and `_max_expr`.

**Other files and tests.** Optional new test: a mechanism with no product-free saturating pattern errors with the descriptive message.

**Evidence.** rate_eq:799-801 builds Expr(:call, :max) with zero candidates and no check, while 881-883 checks for the allosteric method. Both use the same substrate-only filter (789 vs 880), and their sort orders agree (MONO lexicographic). No test checks the message text. Memory records the 'no kcat components' failure family from HPC runs. Measured -11 on the draft.

### deriv/V24 — PLAUSIBLE BUG: allosteric kcat regulator corners treat :OnlyA/:OnlyI ligands as unit concentration
- Net lines removed (this file group): -5
- Risk: medium
- Behaviour change: kcat changes for allosteric mechanisms that have :OnlyA/:OnlyI regulators and a live inactive saturating pattern, and with it the scale rescale_parameter_values applies to fitted SS constants (fitting.jl:267-270). Losses and fits are unaffected.
- Depends on: V21, V22
- Approved by Denis on the sign-off page.

**Change.** TDD: first add a failing numerical test, a :NonequalAI-catalysis mechanism with an :OnlyA regulator, comparing rate_equation at large S and large X against `_kcat_forward`. Then track per-state regulator degrees in the mask loop: d_A = Σ n_site·[site has a saturating ligand that binds A], and d_I likewise. When d_A ≠ d_I, keep only the dominant state's terms (A_A/B_A or A_I/B_I), and skip the corner when that state's pattern is dead (0/0). Equal degrees keep today's weighted form.

**Other files and tests.** One new numerical test in test_rate_eq_derivation.jl (kcat section). Existing specs are unaffected: PFK/HK/PK use :OnlyA catalysis, so B_I = 0.

**Evidence.** I re-derived the limit from the code. With X → ∞ and an :OnlyA ligand, the A-site factor grows like X^n while the I factor stays 1, so v → A_A/B_A. The code instead returns (A_A·W_A + L·A_I)/(B_A·W_A + L·B_I) with W_I = 1. The verifier agreed; this is not confirmed numerically (no Julia run).

### deriv/V25 — Annotate '(substituted into v)' only on ties actually folded by Pass 2
- Net lines removed (this file group): -1
- Risk: low
- Behaviour change: Lines for single-symbol ties that are not folded (iso-K = iso-K, binding-K = iso-K, or a bare `k = Keq` line) lose the suffix, so `_rate_eq_dedup_key` (identify_rate_equation.jl:319-325) stops stripping them. rate_equation_string and eq_hash change for those mechanisms, which are rare.
- Depends on: V18
- Approved by Denis on the sign-off page.

**Change.** In the non-allosteric Reduced string, pass the Pass-2 rename keys (`keys(_build_wegscheider_rename_map(Mechanism(m)))`) to `_equation_text`. Append ANNOTATION_SUBSTITUTED only when `sym in` that set, instead of whenever the RHS is a bare Symbol.

**Other files and tests.** New test: a mechanism with a non-binding single-symbol tie emits no suffix, and a bare `k = Keq` Haldane line likewise.

**Evidence.** rate_eq:655 annotates every Symbol RHS, but only binding-K/binding-K ties are folded (132). The verifier found a second case: build_power_expr returns a bare :Keq. In both cases a runtime-assigned constraint drops out of the dedup key.

### deriv/V26 — Sort the allosteric catalytic independent parameters by name, as the non-allosteric path does
- Net lines removed (this file group): -1
- Risk: medium
- Behaviour change: The order of allosteric fitted_params changes, along with the params destructure line of rate_equation_string, eq_hash and the golden PARAMS_REDUCED lines. Rate-equivalent allosteric mechanisms with different column orders would then dedup; whether such pairs occur is unmeasured.
- Depends on: V16
- Approved by Denis on the sign-off page.

**Change.** In `_dependent_param_exprs(am)`, sort the catalytic part of `indep` by string before appending the regulator and L parameters (+1 line), matching thermo:305-309.

**Other files and tests.** Regenerate test/reference/allosteric_golden_reference.txt and any pins on the order of allosteric fitted_params. If accepted, #31 (classify concentrations by metabolite set in _poly_to_expr, -4 lines) can ride on the same golden regeneration.

**Evidence.** thermo:305-309 sorts the non-allosteric indep 'so rate-equivalent mechanisms produce the identical string'. rate_eq:1411-1412 does not sort. The allosteric params line feeds `_rate_eq_dedup_key`.

### deriv/V27 — Drop allosteric Full mode: parameters(::AllostericEnzymeMechanism, Full) and the allosteric full enumeration
- Net lines removed (this file group): 50
- Risk: medium
- Behaviour change: `parameters(aem, Full)` becomes a MethodError. It is an exported function, but no rate_equation or rate_equation_string Full method exists for allosteric mechanisms, and the docs use Full only for plain mechanisms (textbooks.md:99-108).
- Depends on: V12, V13
- Approved by Denis on the sign-off page.

**Change.** Restrict the Full parameters method to EnzymeMechanism and drop the allosteric over-emission comment. Delete `_enumerate_parameters_full(::AllostericMechanism)` with its 24-line docstring, and `_all_i_state_parameters`, which then has no src caller. Give `_ss_rate_constant_names` an allosteric method that filters `[_cat_params(am, :A); _cat_params(am, :I)]`. Trim the parameters docstring (28-36).

**Other files and tests.** Tests: delete the allosteric Full asserts (test_accessors.jl:32-39) and the 'parameters(Full) injective for collapsed allosteric shapes' testset (test_rate_eq_derivation.jl:1465-1468). Drop the allosteric PARAMS_FULL golden lines; test_allosteric_golden.jl:17 must skip allosteric specs. Rewrite the `_all_i_state_parameters` testset (test_rate_eq_derivation.jl:2296-2330) onto `_cat_params`/`_kreg_params`.

**Evidence.** grep finds no src call of parameters(…, Full), and the docs use Full only on plain mechanisms. Measured -50 on a copy of the draft taken after V12/V13.

## src/thermodynamic_constr_for_rate_eq_derivation.jl
Excess code in this file comes from five sources.
(1) Three exact Gauss-Jordan routines: _rref_partition, _rational_nullspace, and the priority-pivot loop in _solve_dependent_set. One RREF can serve all three (T1, T2).
(2) Hand-rolled bookkeeping in _thermodynamic_constraints: a 27-line classifier that writes one error string three times, a separate species walker, and separate substrate/product vectors (T3).
(3) The constraint sign convention is written out in _assemble_constraints and again in the split counter, then undone in _onlya_haldane_violation. On top of that, a redundant all_params kwarg must agree with step_params, and names are re-rendered for every nonzero cycle entry (T6, T7).
(4) Forwarding layers with stale docstrings (the kernel and its Type wrapper, which names a function that no longer exists), plus rate-body builders that belong next to @generated rate_equation, not here (T8, T9).
(5) Small re-implementations: the 4-way parameter switch, and a copy-then-delete _free_enz_set (T4, T5).
Every behaviour-preserving item keeps C, the column order, the priorities and the generated Expr identical. The Pass-2 rename is not touched here.
Two items change outputs and need Denis's decision: whether a group's pivot priority comes from its last step or its naming rep (T10), and whether allosteric independent parameters are name-sorted (T11). T2 also changes an error message's wording.
About 263 of 776 lines go; in package terms about 30 of those move to rate_eq_derivation.jl (see T8/T9).

### thermo/T1 — One exact RREF: build _rational_nullspace and _integer_nullspace on _rref_partition
- Net lines removed (this file group): 28
- Risk: low
- Behaviour change: none. RREF is unique, so pivots, free columns and both nullspace matrices are bit-identical, and so is the cycle basis C that feeds pivot choice.
- Depends on: —

**Change.** Widen `_rref_partition(A::Matrix{Int})` (thermo:131) to `AbstractMatrix`. Its first line `Matrix{Rational{BigInt}}(A)` already copies Int or Rational input. Rewrite `_rational_nullspace` (thermo:466-494, a line-for-line second Gauss-Jordan) as the basis build on `_rref_partition`: free column k gets 1 at its free coordinate, and pivot row r gets `-R[r, fc]`. Rewrite `_integer_nullspace` (151-171) as a column-wise scaling of `_rational_nullspace(A)`: multiply by the lcm of the denominators, divide by the gcd, make the first nonzero positive. Delete its `isempty(free_cols)` early return and the dead `g > 0 &&` / `fnz !== nothing &&` guards; every basis column holds a 1 at its free coordinate, so neither guard can fire. Optionally move `_rational_nullspace` next to `_rref_partition`. Update the 'Used by' line of the `_rref_partition` docstring. Keep the `_rational_nullspace` name and signature, since tests call it.

**Other files and tests.** No other src file. Tests: none edited. The test_types.jl:2728-2760 `_rational_nullspace` testset (size and M*N==0 checks) keeps its API. `_partition_independent_count` keeps calling `_rref_partition` on an Int matrix.

**Evidence.** Read thermo:125-171 and 461-494. Both routines do the same pivot search (first nonzero row), swap, normalize and eliminate-all-other-rows. `r > m && break` in `_rational_nullspace` is equivalent to an empty `row:m` findfirst. The basis build at 485-492 equals NS at 155-159. Callers: `_integer_nullspace` only at thermo:211; `_rational_nullspace` at thermo:630 plus tests; `_rref_partition` at 153 and 368. Verifiers agree on about 28 lines. The Bareiss/det half of the 'exact linear algebra' proposal is in rate_eq_derivation.jl and belongs to the derivation group.

### thermo/T2 — _solve_dependent_set as column-order RREF on priority-sorted columns
- Net lines removed (this file group): 32
- Risk: high
- Behaviour change: Same dependent set and the same dependent expressions (matroid-greedy argument, independently re-verified). `build_power_expr` sorts factors by name, so the Expr is identical. Every dep consumer sorts by key or is order-independent. Only change: the 'Thermodynamically contradictory mechanism' error can no longer name the original row index or the unnormalized coefficient.
- Depends on: T1
- Approved by Denis on the sign-off page.

**Change.** Replace the hand-rolled row-greedy priority pivoting in `_solve_dependent_set` (thermo:647-701) with: `order = sortperm(priority; rev = true)`; `pivots, _, R = _rref_partition([A[:, order] rhs])`; error when `n + 1 in pivots`, which is today's `0 = c·log Keq` contradiction (rank[A rhs] > rank A). Then for each `(r, p)` in `enumerate(pivots)`, set `dep_exprs[columns[order[p]]] = build_power_expr(R[r, n+1], [(columns[order[c]], -R[r, c]) for c in 1:n if c != p && R[r, c] != 0])`. indep is the columns not in dep, in column order, as today. Rewrite the docstring to state why this matches the old rule: the per-row highest-priority unused column is the lexicographic-greedy basis, which is the column-order RREF pivot set on columns sorted by (priority desc, index asc). Stable sortperm breaks ties by ascending index, as the strict `>` scan does. Bool-first tuples sort with true first.

**Other files and tests.** None in src. The allosteric combined solve (rate_eq:1355) calls it with the same signature. Tests: no edits; nothing calls it directly or greps its error text. Oracles: the full suite, especially the fitted_params pins over MECHANISM_TEST_SPECS, test_pivot_priority_regression.jl, test_allosteric_golden.jl and the jldoctest.

**Evidence.** Read thermo:647-701 and sym_poly:225-256 (factors are sorted by name). Every pivot eliminates its column from all other rows, so the final row space is in RREF and each dependent row is unique. Zero rows keep the rhs they had when processed, so the error fires iff the rhs column is a pivot. `grep contradictory|Thermodynamically` in test/ and docs/: 0 hits. Every iteration over a dep map in src sorts by key, except the rename fold (rate_eq:130, 1226), which is order-independent because its lhs are pivots and its rhs are free columns. The verifier confirmed 33 lines; my draft comes to about 32, with a slightly longer docstring.

### thermo/T3 — Rewrite _thermodynamic_constraints compactly; delete _enumerate_species_names
- Net lines removed (this file group): 48
- Risk: low
- Behaviour change: none. B has the same rows and columns, so C is unchanged, and the error text is byte-identical. Accept and reject are the same: an all-zero nu gives c = 0 and passes; a nonzero entry where nu_net = 0 fails `nu == c .* nu_net`; unequal ratios fail; a non-integer ratio fails `isinteger`.
- Depends on: —

**Change.** In `_thermodynamic_constraints` (thermo:173-258):
(a) Build the form index inline, in flat-step order, with `get!(form, name(sp), length(form) + 1)`, and delete `_enumerate_species_names` (260-273); its single caller is 175.
(b) Merge the consumed and released loops (202-209) into one loop over `((consumed(s), -1), (released(s), 1))`. Keep the `haskey(met_idx, …)` guard, because competitive-inhibitor copies are consumed but are not reactants.
(c) Build `nu_net` in one loop over `reactants(mech.reaction)`, adding ±1 by `metabolite(ra) isa Substrate`, and drop `subs_species`/`prods_species` (179-180, 215-221). This is the same name-indexed accumulation. The `_reactant_signs` variant was rejected because it lives in a later file and indexes differently.
(d) Drop the `nc == 0` early return (211-213). `C = _integer_nullspace(B)'` is already 0×nsteps, and `Int[classify(i) for i in axes(C, 1)]` gives `Int[]`.
(e) Replace the 27-line `classify_cycle` closure (227-253) with `j0 = findfirst(!iszero, nu_net)` and a closure: `nu = stoich_mat * C[i, :]; c = nu[j0] // nu_net[j0]; nu == c .* nu_net && isinteger(c) || error("Cycle $i produces metabolite change not proportional to net reaction"); Int(c)`.
Keep the 9-line stoichiometry comment (191-199) and the 4-line classification comment (223-226).

**Other files and tests.** None in src. Tests: no edits. test_identify_rate_equation.jl:1532-1562 (m_bad) still raises the same error, which surfaces as a FitFailure.

**Evidence.** Read thermo:173-273. `grep -w _enumerate_species_names` over src/test/docs finds only the call at 175 and its own definition. `get!` evaluates the default before inserting, so the first-seen index order is the same. `reactants` holds only Substrate/Product (types.jl:478-481). Duplicate substrate or product names are rejected (types.jl:415-418), so the name-indexed nu_net equals today's. The verifier counted 45 lines; my draft is about 52 lines against today's 101 plus a blank (about 48). The other proposals' _forms, RE-segment and unique versions touch this file only through _enumerate_species_names, which goes here; the rest of those proposals belongs to the types, derivation and enumeration groups.

### thermo/T4 — _step_parameters as one comprehension over _emit_cat_params_for_rep
- Net lines removed (this file group): 9
- Risk: low
- Behaviour change: none. The same Vector{Vector{Parameter}} is produced, with the same element order (Kon,Koff / Kfor,Krev), each anchored on its own step.
- Depends on: —

**Change.** Replace the body of `_step_parameters` (thermo:36-46), which hand-writes the same RE/SS × binding switch as `_emit_cat_params_for_rep` (types.jl:1937-1945), with `_step_parameters(m::Mechanism) = Vector{Parameter}[_emit_cat_params_for_rep(s, :None) for (s, _) in _flat_steps(m)]`. Keep the docstring, which says each Parameter is anchored on its own step. Use a typed comprehension so the element type does not depend on inference.

**Other files and tests.** rate_eq_derivation.jl:1139-1152: `_state_step_params` can use the same call with its state tag (about -4 lines; derivation group). Optional in types.jl: the name `_emit_cat_params_for_rep` and its docstring wording ('representative') become misleading once it is applied to every step; a rename such as `_step_constants` touches about 6 call sites (types group). Tests: none. They index or forward the result (test_rate_eq_derivation.jl:140, test_types.jl:2095/2292, test_mechanism_enumeration.jl:9255/9335).

**Evidence.** thermo:39-42 matches types.jl:1938-1944 token for token apart from the anchor. The helper reads only is_equilibrium and is_binding of its argument, never whether the step is a rep. Several verifiers confirmed. The rename-map part of the plumbing-unification proposals lives in rate_eq_derivation.jl (derivation group); only this switch is in this file.

### thermo/T5 — _free_enz_set as set comprehensions; reword its false SS-direction clause
- Net lines removed (this file group): 15
- Risk: medium
- Behaviour change: none. This relies on enzyme-form names being injective, which the whole derivation already assumes since it keys forms by name.
- Depends on: —

**Change.** Rewrite `_free_enz_set` (thermo:71-95), which today copies a set and then deletes from it over three passes, as: `flat = [s for g in steps(m) for s in g]`; `into(re) = Set{Symbol}(name(to_species(s)) for s in flat if !isempty(consumed(s)) && (!re || is_equilibrium(s)))`; `bound_in, re_bound_in = into(false), into(true)`; return `Set{Symbol}(name(sp) for s in flat for sp in (from_species(s), to_species(s)) if name(sp) ∉ re_bound_in && (isempty(bound(sp)) || name(sp) in bound_in))`. The result must be `Set{Symbol}`, because `_step_priority`'s signature requires it. Keep the docstring. Its clause 'SS steps' direction is not canonicalized so they don't participate' (55-56) is provably false (types.jl:556-560 orients RE and SS steps alike), so reword it to state the rule without the false reason.

**Other files and tests.** None. Tests: none edited. Parameter-name pins across the suite are the oracle, because this set picks every kinetic group's naming rep (`_group_rep`, `_rep_step` at types.jl:1778).

**Evidence.** Today a name is free iff (a) it is not the to-side of an RE consuming step, and (b) every occurrence is unbound or the to-side of some consuming step. Under name injectivity, 'every occurrence' equals 'any occurrence'. Callers: types.jl:1778, 1882, 1893, 1917; rate_eq:704, 1249, 1289; thermo:400; test_rate_eq_derivation.jl:297. The verifier confirmed about 15 lines. The _forms proposals' use of this function is subsumed, since the comprehension needs no enz_names list.

### thermo/T6 — _onlya_haldane_violation reads the ε exponents straight off the cycle basis
- Net lines removed (this file group): 24
- Risk: medium
- Behaviour change: none. M has the same rows, the same columns in the same order (sorted column index into cm's group order) and the same values: RE Kd is −ΣC times −1, SS Kon is +ΣC times +1. Offender strings are identical, so the FM blow-up cap behaves the same.
- Depends on: —

**Change.** In `_onlya_haldane_violation` (thermo:577-635), stop building `_step_parameters(cm)`, running `_assemble_constraints(cm, Dict())`, rendering names to find columns and undoing the Kd sign flip with a per-column multiplier (589-605). Instead:
1. `onlya(g) = cat_allo_states[g] === :OnlyA`.
2. Build `keep` as today, then `any(onlya, keep) || return nothing`.
3. `cm = Mechanism(rxn, kept)`, `flat = _flat_steps(cm)`.
4. `groups = unique(g for (s, g) in flat if s in onlyA_steps)` gives the :OnlyA binding groups in cm's group order. Add `isempty(groups) && return nothing`, which keeps today's accept when no column survives.
5. `C, _ = _thermodynamic_constraints(cm)`.
6. `M = Rational{BigInt}[sum((C[i, j] for (j, (_, g)) in enumerate(flat) if g == k); init = 0) for i in axes(C, 1), k in groups]`.
7. Offender label of column k: `string(name(first(_emit_cat_params_for_rep(<first step of group k>, :None)), cm))`.
8. Per-row test: the nonzero columns of a row are all one sign → message.
9. Then `_has_strict_positive_combination(_rational_nullspace(M)) === false || return nothing` and the message over all columns.
The docstring paragraph at 568-572 (Kd/Kon normalization) becomes false: reword it to say an :OnlyA group's ε exponent is the sum of its steps' cycle entries for RE and SS alike. Drop the stale 'line 366' mirror comment (594-598). Keep the two stage comments (607-609, 622-627).

**Other files and tests.** None in src. Callers are unchanged: types.jl:922 and enumeration 2679/2902/2918, all of which test `=== nothing`. Tests: no edits. test_types.jl:2490-2724 pins both stages, including the mixed RE/SS balanced accept. Removes one `_assemble_constraints`, one `_step_parameters` and O(nc·nsteps) `name` renders from every AllostericMechanism construction and enumeration-move check. This lowers generator-time and test time, but the gain is unmeasured.

**Evidence.** Read thermo:577-635 together with _assemble_constraints' sign rule (424-430). Binding steps are never re-oriented (types.jl:566), so `s in onlyA_steps` survives `Mechanism(rxn, kept)`, as today's code already assumes. Message text is never grepped in tests; they check `isa String` or `=== nothing`. allequal needs Julia ≥ 1.8 (compat effective 1.11). The verifier counted 25; my draft is about 37 lines against today's 59, plus 2 docstring lines saved (about 24).

### thermo/T7 — _assemble_constraints: columns from step_params, one render per step, one sign-convention fold shared with the split counter
- Net lines removed (this file group): 31
- Risk: medium
- Behaviour change: none. Columns, column order, A entries and priorities are identical, so the pivots and fitted_params are unchanged.
- Depends on: —

**Change.** In `_assemble_constraints` (thermo:392-459):
(1) Drop the `all_params` kwarg. Set `columns = _param_columns(mech, step_params)`, defined in this file as `unique(Symbol[name(p, m) for ps in sp for p in ps])`. This is `_state_all_params` (rate_eq:1192-1200) moved and renamed. It equals `_raw_param_symbols(mech)` element for element and in order, because groups are uniform in kind and RE/SS (`_assert_uniform_groups`) and no reaction appears twice (`_assert_each_reaction_once`).
(2) Delete `binding_K_set` (412-416). The sign key is `_count_kind(s)`, which is `:binding_K` exactly for an RE binding, since the rename maps binding K to binding K only.
(3) One loop over the flat steps that renders each step's names once (`cols = [sym_col[get(rename, n, n)] for n in (name(p, mech) for p in step_params[j])]`), folds the cycle column through a new helper `_add_step_column!(A, C, j, kind, cols)` (binding K negated, iso K as is, SS +kf/−kr) and sets the priority. This keeps the last-write-wins order. `against` now also requires `!is_equilibrium(s)`, so RE steps keep rank 0. Delete the dead `haskey(sym_col, s)` guards (444, 454). Keep the two priority comments (434-437, 446-449).
(4) The `_partition_independent_count` counter's fill loop (357-367) calls the same `_add_step_column!`. Keep its two length guards (348-351) for robustness. Make the counter docstring's sign-convention sentence (339-341) point at the helper. If T6 has not landed, also repoint the mirror comment in `_onlya_haldane_violation`.

**Other files and tests.** rate_eq_derivation.jl: delete `_state_all_params` (−13; moved here as the one-line `_param_columns`). Rename its callers at 864, 1496-1497 and 1555, and drop the `all_params = …` kwarg at 1225 and 1326 (−2). That is about −15 there. It also simplifies the derivation group's rename-map merge, which no longer threads all_params. Tests: no edits. test_mechanism_enumeration.jl:9744-9806 (counter == _independent_param_count) guards the shared fold.

**Evidence.** Read thermo:343-459 and rate_eq:1192-1236 and 1321-1356. Today `step_name` renders `name(p, mech)` twice per nonzero C[i,j], and each render rebuilds `_free_enz_set` inside `_rep_step` (types.jl:1777-1782). The merged loop renders once per parameter. The flat step s and `step_params[j][1].step` are structurally equal for both builders; a direction mismatch would already break `_rep_step`. `_step_priority` and `against` are symmetric under reversal for non-bindings. No tests call `_assemble_constraints`, the kernel, or the counter guard messages. Verifier estimate: 35. Mine is about 31 because I keep the counter's length guards and add the 4-line `_param_columns`.

### thermo/T8 — Delete _dependent_param_exprs_kernel and its stale Type wrapper
- Net lines removed (this file group): 37
- Risk: low
- Behaviour change: none
- Depends on: T7

**Change.** After T7 the kernel (thermo:703-730) only forwards: `_solve_dependent_set(_assemble_constraints(mech, rename; step_params)...)`. Its docstring describes a 'Pass-1-only / user-defined kinetic-group rename' that no longer exists; every such caller passes an empty Dict. The Type wrapper (731-738) has one caller, and its docstring names `_build_kinetic_rename_map`, which exists nowhere. Delete both. Write the solve inline at the four call sites: thermo:291, rate_eq:129, rate_eq:673 (pass `Mechanism(m)`) and rate_eq:1223. Move the one still-true docstring fact to the Reduced `rate_equation_string` (rate_eq:669): the display path solves with the empty rename so absorbed single-symbol ties stay visible under `# Wegscheider constraints:`. Reword the `_dependent_param_exprs` docstring (284-287), which mentions forwarding to the kernel. An alternative if Denis wants a named entry point is to keep a 1-line kernel with a 3-line docstring (net about 32).

**Other files and tests.** rate_eq_derivation.jl: the three call sites are rewritten in place (0 net), and about 3 docstring lines move onto rate_equation_string (+3). Package net is about −34. Tests: none; nothing in test/ or docs/ calls the kernel.

**Evidence.** `grep _dependent_param_exprs_kernel` finds thermo:291 and rate_eq:129, 673 and 1223 only; the Type method is reached only from 673. `grep _build_kinetic_rename_map` over src/test/docs finds only thermo:734. The rename map at rate_eq:117 is empty when the kernel runs, so 'Pass-1-only' is false today. The Pass-2 explanation already lives in the `_build_wegscheider_rename_map` docstring, so nothing load-bearing is lost.

### thermo/T9 — Move rate-body assembly out of this file; delete _sorted_raw_param_symbols
- Net lines removed (this file group): 38
- Risk: medium
- Behaviour change: none (byte-identical generated Expr and strings)
- Depends on: —

**Change.** Delete the 'Preamble Building Helpers' section (thermo:740-776):
- `_sorted_raw_param_symbols` (751-753). Its single caller is rate_eq:518. It is misnamed: it keeps group order and equals `(_raw_param_symbols(m)..., :E_total)`, which is also `parameters(m, Full)`. Inline that at rate_eq:518.
- The two `_build_rate_body` methods (755-776). Each has one caller, the @generated `rate_equation` methods at rate_eq:556-567. Inline them there, or into the derivation group's unified body builder if that item lands. Share one `_dep_assignments(dep)` for the sorted assignment comprehension duplicated at thermo:769-770 and rate_eq:797-798.
- `_destructuring_expr` (742-745) moves beside its other callers in rate_eq_derivation.jl (803, 986, 1647-1648).
Update ABOUTME line 2 and the module docstring (4-13), which say this file builds rate-equation preambles. The generated Expr must stay byte-identical: same destructure order, the same sorted assignments, and the same `:(E_total * (num) / (den))`, with no added `return`.

**Other files and tests.** rate_eq_derivation.jl gains about +13: `_destructuring_expr` (+4), the two bodies inlined into the @generated methods (about +12), and one shared `_dep_assignments` (−3 at `_kcat_forward`). Package net is about −25; less is added if the derivation group's body unification (one front end for plain and allosteric) lands instead. The docstring mention of `_build_rate_body` at rate_eq:648 needs rewording. `_raw_param_symbols` (thermo:17-27) then has no caller in this file; leave it or move it beside `_enumerate_parameters_full` (0 net). Tests: none, unless the derivation group changes the return shape of `_raw_rate_expr_and_symbols` (test_rate_eq_derivation.jl:1675). Gate: test_rate_equation_performance (0 allocs, <120 ns) and the Expr-shape and flat-string tests.

**Evidence.** `grep -rn` over src/test/docs: `_build_rate_body` has its definitions at thermo:756/765, calls at rate_eq:559/566 and a docstring mention at 648. `_sorted_raw_param_symbols` has one caller (rate_eq:518) and no tests. `_destructuring_expr` has callers at thermo:759-773 and rate_eq:803, 986, 1647-1648. This file's lines are counted here; the num/den sharing, the parameters(Reduced) merge and the v-line unification in the overlapping proposals are rate_eq_derivation.jl edits for the derivation group.

### thermo/T11 — Decision: name-sort the allosteric independent parameters (move the sort into _solve_dependent_set)
- Net lines removed (this file group): 1
- Risk: medium
- Behaviour change: If approved, it reorders allosteric fitted_params, the NamedTuple field order of the params line, the rate_equation_string params line and the eq_hash of allosteric candidates. Non-allosteric output is unchanged.
- Depends on: —
- Approved by Denis on the sign-off page.

**Change.** `_dependent_param_exprs(::Mechanism)` sorts `indep` by string (thermo:305-309) so that rate-equivalent mechanisms render the same params line and eq_hash. The allosteric `_dependent_param_exprs` (rate_eq:1374-1413) returns the combined solve's column order (A-state columns, then I-only columns), so the params line is not canonical, for example `K_A_EP…, K_A_ES…, k_A_ES…, K_I_EP…`. If Denis approves, move the sort and its comment into `_solve_dependent_set`'s return. Both paths then get it, and the plain path's separate sort line goes. Sorting before the Pass-2 filter keeps the order. Regulator Ks and :L stay appended after, as today. Land after T2, since both edit the same return line.

**Other files and tests.** No code change in rate_eq_derivation.jl, but the allosteric output order changes. Regenerate all 12 REDUCED_STRING entries in test/reference/allosteric_golden_reference.txt and any allosteric fitted_params order pins. Separate coverage gap from the same proposal, for the tests group: the allosteric Pass-2 reference guard (rate_eq:1384-1394) runs only for a non-empty state rename, which no current fixture produces. Add a DSL fixture with a grouped RE binding square inside @allosteric_mechanism.

**Evidence.** thermo:309 sorts; rate_eq:1411-1412 does not. The golden file's first params line, `K_A_EP_to_E_P, K_A_ES_to_E_S, k_A_ES_to_EP, K_I_EP_to_E_P, …`, is not string-sorted, because 'K' < 'k' would group all the K_ names first. The verifier marked it needs_measurement: whether rate-equivalent allosteric mechanisms with permuted columns actually collide in a run is unmeasured. The params line feeds `_rate_eq_dedup_key` (identify_rate_equation.jl:319-325).

## src/mechanism_enumeration.jl
1. Seed generation (lines 1-1231) holds most of the excess code. The backtracker has four branches and repeats the bind block four times and the isomerization blocks twice. Eleven helpers thread state that could be derived from the arguments. The dead-end builder makes three passes, and its last pass (the connectivity fill) already adds every step the first two add. A group-id array runs alongside the steps but carries no information, and five bitmask subset loops are written out by hand.
2. Graph primitives are written several times. RE connectivity has two union-finds, the `_re_segment_extras` BFS, the `_indexed_re_segments` wrapper and a node BFS. Cycle balance has Tarjan, a second BFS for potentials and the screen's union-find. The minimal-flip search appears twice, once in `_expand_re_to_ss` and once in `_seed_variants`.
3. The regulator and allosteric moves copy vectors by hand and recompute values that do not change inside the loop. There are four no-op wrong-kind methods, three functions implementing one reg-type rule, several single-use helpers, and about 70 lines of src that only tests use (`_assert_mechanism_invariants`).
4. The proposals overlap heavily. 98 proposals collapse into 26 items.
- Items that preserve behaviour and need no permission save about 697 lines (668 if S1 is declined, because S2 then saves 21, not 50).
- The two rewrites need Denis's go-ahead under CLAUDE.md: S1 (walker, -555) and S4 (dead-end fill, -175).
- Behaviour-changing or perf-gated items add 133 lines net.
- Total: about 1,560 of 3,415 lines (46%).
Every load-bearing invariant stays untouched: the Mechanism constructor still canonicalizes every child, and the rate_equation codegen and naming are unchanged.

### enum/S0 — One `_subsets` helper replaces the bitmask subset loops in the competition-pattern generators
- Net lines removed (this file group): 56
- Risk: low
- Behaviour change: none (same subsets in the same order)
- Depends on: —

**Change.** Add `_subsets(v)` near the topology helpers: a docstring and a 2-line definition giving every subset in binary-counting order, empty first, with bit i-1 selecting v[i].
- Rewrite `_competition_patterns` (886-919) as `[Set(pat) for pat in _subsets(edges)[2:end] if <every substrate covered> && <every product covered>]`.
- Rewrite `_inhibitor_competition_patterns` (921-972) as `[(Set(s), Set(p), Set(i)) for s in _subsets(subs)[2:end] for p in _subsets(prods)[2:end] for i in _subsets(inhs)]`. The substrate loop stays outermost and the inhibitor loop innermost, and the inhibitor subsets still start from the empty set.
- Keep both signatures and docstrings. Do not inline `_inhibitor_competition_patterns`; that would drop its count testset.

The other mask loops move onto `_subsets` elsewhere: `_wo_recurse!` in S1 and the K-type loop in A6.

**Other files and tests.** No other src lines. The derivation group may reuse `_subsets` for the kcat regulator corners (rate_eq_derivation.jl:931-939). Tests: none. The order and count pins at test_mechanism_enumeration.jl:713-831 still hold.

**Evidence.** Four of five verifiers confirmed that mask order and comprehension nesting match the loops at 903-909 and 945-965, with an estimate of 56. My count is 26+42 code lines down to about 9, plus the 3-line helper.

### enum/S1 — Rewrite `_catalytic_topologies` as one walker whose state is derived from its arguments
- Net lines removed (this file group): 555
- Risk: medium
- Behaviour change: None. The topologies are identical. Only the internal order in which paths are generated changes, and the Set unions and sorted outputs erase it. A latent tie exists today and stays: two distinct forms that render the same name.
- Depends on: S0
- Approved by Denis on the sign-off page.

**Change.** Delete `_make_species` and `_residual_for` (7-45), `_can_pingpong`, `_subtract_atoms` and `_add_atoms` (64-99), `_combinations`, `_release_products!`, `_binding_order`, `_release_order` and `_linearizes` (198-300), and the backtracker body (309-821). Replace them with three pieces (draft: scratchpad/enum_seeds_sketch.jl:1-133):
(a) `_net_atoms(rxn, plus, minus)`: the signed atom multiset over reactant names, with zero counts dropped.
(b) A top-level `_weak_orderings(items)` built on `_subsets`, in the same order as today.
(c) A recursive `walk!(cur, consumed, released, on_subs, on_prods, path)`. The atoms on the enzyme are atoms(consumed) - atoms(released); `post_final` is `!isempty(on_prods)`; the ping-pong flag is implied.

Per node, `walk!`:
1. Releases each bound product.
2. Otherwise binds each remaining substrate.
3. Stops if no substrate is bound or more than 3 are (C7, C6).
4. Otherwise isomerizes to each nonempty subset of the remaining products. It skips a subset whose `residue = _net_atoms(rxn, consumed, [released; subset])` is negative. It keeps a subset only if `isempty(residue) == isempty(remaining_subs)` (the admissible-residual rule) and `length(subset) + !isempty(residue) <= 3`.

Post-processing:
- Drop the `unique_paths` pass; it is a no-op.
- Compute the weak orderings once per reaction, with a rank Dict per ordering, and each path's binding and release orders once.
- Keep exactly: the iso-pattern sort, the weak-ordering iteration order, the topology step sort, and the rule that the first isomerization is steady state.

Keep the union-of-paths and determinism comments and the C6, C7, C8 and admissible-residual comments. Drop the vacuous C5 guard and `max_bound`. Keep `_atoms_dict` until S2.

Before deleting the old code, run a one-off script asserting old(r) == new(r), as `Vector{Vector{Step}}` equality, for uni-uni, uni-bi, bi-bi, bi-bi ping-pong, ter-bi, ter-ter, the admissible-residual ter-ter, pyruvate carboxylase, PDH and quad-quad.

**Other files and tests.** No other src. Tests: no change is needed, because tests call only `_catalytic_topologies` (15 references) with the same signature. The pins are the counts 1/3/9/10/223/45/312/334, the init counts 239/264/250,855 and the bi-bi structural golden at test 12115-12149. The tests group may delete the vacuous C5 assertion loop (test_mechanism_enumeration.jl:547-555).

**Evidence.** All four rewrite proposals were confirmed, with estimates of 555, 270, 200 and 78 depending on scope. I traced every branch against the sketch. The admissible-residual rule differs from today's only on routes that dead-end anyway. Atom balance (types.jl:449-467) and reactant sorting make the outputs equal. I count 692 lines deleted and about 131 added.

### enum/S2 — Rebuild the atom-conservation check on `_net_atoms`
- Net lines removed (this file group): 50
- Risk: low
- Behaviour change: none
- Depends on: S1

**Change.** Rewrite `_assert_step_atom_conserving` (133-151) as one `_net_atoms` call:
- from-side names: bound and residual-added names of from_species, residual-subtracted names of to_species, and consumed;
- minus to-side names: the mirror set.

Delete `_atoms_dict` (47-63), `_accumulate_atoms!` (103-109), `_nonzero_atoms` (111-112) and `_species_atoms` (114-131), and fold `_species_atoms`'s description into the remaining docstring.

Keep `_assert_atom_conserving` and all four call sites: init, `_seed_variants` bases, the per-child check in expand_mechanisms, and `_base_tier`. Keep the error message byte-identical, including its 'atom-non-conserving step X → Y' prefix. Its label states from - to while the printed diff is to - from; fixing that is a separate one-line message fix.

Full deletion of the guard (from #11) is not adopted:
- It is the only error-injection route for the `_expand_parent` and `_base_tier` failure-isolation tests (identify 1580-1615 and 1777-1805).
- enumeration_engine.md:145 documents it.
- Denis asked for more robustness.

**Other files and tests.** None. Tests at test_identify_rate_equation.jl:1611 and 1803 check only the message prefix. Without S1 this item saves about 21 lines: `_net_atoms` is added here and `_atoms_dict` stays for the old walker.

**Evidence.** The verifiers confirmed the guard cannot fire in production, but they also found two tests that use it as their only error trigger. #11's own fallback is to keep the guard and rewrite it on `_net_atoms`. `_atoms_dict` matches by name, so `_net_atoms` gives the same result except when a substrate and a product share a name, which is already rejected downstream.

### enum/S5 — Fail loudly when a reaction has no admissible catalytic cycle
- Net lines removed (this file group): -2
- Risk: low
- Behaviour change: For reactions with no route, such as S → P+Q+R+T or a quad-uni with no ping-pong route, init_mechanisms, seed_mechanisms and identify_rate_equation now throw a clear error. Today identify silently returns an empty result, and seed_mechanisms raises a misleading 'no mechanism binds every required regulator'.
- Depends on: —
- Approved by Denis on the sign-off page.

**Change.** At the end of `_catalytic_topologies`, add `isempty(result) && error("no catalytic cycle for $(reaction): every route needs an isomerization converting more than three substrates or products at once")`, wrapped at 92 characters. Write the failing `@test_throws` on a uni-quad `@enzyme_reaction` first (TDD). The check fits either the old backtracker or the S1 walker.

**Other files and tests.** Tests: +4 lines (one `@test_throws`). No test expects an empty result; quad-quad still has ping-pong topologies.

**Evidence.** Confirmed. Traced uni-quad and quad-uni to an empty result. Reactions with at most three substrates and at most three products always have the sequential route, so the error fires only in unsupported cases.

### enum/S3 — Group seed steps by kind: `_seed_groups` replaces `_apply_equivalence_grouping` + `_to_group_list`
- Net lines removed (this file group): 40
- Risk: low
- Behaviour change: none (identical Mechanisms in identical order)
- Depends on: —

**Change.** Replace `_to_group_list` (1187-1203) and `_apply_equivalence_grouping` (1205-1230) with `_seed_groups(steps)`. It groups steps in a Dict keyed by `(_step_kind(s), is_equilibrium(s))`, with every isomerization under a key of its own. The Mechanism constructor canonicalizes group and step order, so Dict order does not leak.

In init_mechanisms (3184-3201), build `mechs = [Mechanism(r, _seed_groups(steps)) for (steps, _) in _expand_substrate_product_dead_ends(_catalytic_topologies(r), r)]` and then `foreach(_assert_atom_conserving, mechs)`. Keep the seen/out variant loop as it is; #47's `unique!` would change the count if duplicate seeds exist, which has not been checked for ter-ter.

Reword the init docstring, which names `_apply_equivalence_grouping` at 3178. Replace the orphan `# --- Dedup ---` header (3166-3167) with a `# ─── Entry points ───` section header in the file's style.

**Other files and tests.** Tests: the `_topo_mech` helper (test_mechanism_enumeration.jl:34-36) becomes `EnzymeRates.Mechanism(rxn, [[s] for s in t])`, and the comment at 1404 names `_seed_groups`. The bi-bi structural golden (12115-12149) and 'Same-metabolite RE bindings share kinetic_group' (1403-1432) gate the change.

**Evidence.** Confirmed: the parallel id array carries no information. A mirror always shares its parent's key, a singleton keeps a unique id, and an isomerization is never mirrored. The final partition is therefore key-based. Canonical group order is total.

### enum/S4 — Rewrite `_expand_substrate_product_dead_ends` as dead-end species plus one connectivity fill
- Net lines removed (this file group): 175
- Risk: medium
- Behaviour change: None at the Mechanism level. Only the flat pre-constructor step order changes, and the constructor canonicalizes it.
- Depends on: S3
- Approved by Denis on the sign-off page.

**Change.** Follow the #10 sketch (scratchpad/enum_seeds_sketch.jl:167-220). Per topology:
- Compute the catalytic Species once, as unique from/to species.
- Compute the dead-end Species directly: a form plus one reactant it does not hold, kept when the result holds both substrates and products but not all of either.
- Per competition pattern, filter the active dead ends and dedup by the Species vector.
- Run one fill over [catalytic; active]. For each ordered pair of forms one metabolite apart, with the same conformation and residual, that is not already a topology edge, push an RE plain-binding Step.
- Return `Vector{Vector{Step}}`.

This deletes `_substrate_product_dead_end_opportunities` (825-884), the explicit dead-end and mirror passes, the name-keyed de_forms/de_bound/present/have_edge bookkeeping, the `_add` closure and the groups/next_g threading.

State both preconditions in the docstring: every input binding is RE, and every isomerization joins a substrate-only form to a product-only form. Keep the shared-site and connectivity comments.

Gate: old_init(r) == new_init(r), order-sensitive, for uni-uni, uni-bi, bi-bi, bi-bi ping-pong, ter-bi and the shared_catalytic_site reaction, plus the ter-ter count 250,855. Keying by Species also removes a latent name collision (E(NAD,H) and E(NADH) render the same name).

**Other files and tests.** Tests:
- The `_expand_substrate_product_dead_ends` tests change `for (steps, _groups) in result` to `for steps in result` (about 10 sites, test_mechanism_enumeration.jl:1020-1438).
- The `_substrate_product_dead_end_opportunities` testset (929-1004) is rewritten through `_expand_substrate_product_dead_ends`: a union of 27 forms, and exactly 12 forms on the diagonal pattern.
- The test helpers `_boundmap` and `_form_species` become unused.

**Evidence.** Confirmed: the fill subsumes the dead-end bindings and mirrors. Every topology binding is RE, there are no isomerization mirrors, and the Step constructor stores single-metabolite releases as bindings, smaller form first. Verifier estimate 175 (258 lines deleted, about 79 added).

### enum/T1 — Move `_assert_mechanism_invariants` out of src into a `_testhelper_`
- Net lines removed (this file group): 70
- Risk: low
- Behaviour change: none for the package (no src caller)
- Depends on: —

**Change.** Delete both methods (3346-3415). Add `_testhelper_assert_mechanism_invariants(m)` (about 20 lines). It checks that the mechanism has steps and that every declared substrate and product appears in some step, for both mechanism types through `steps(m)` and `reaction(m)`.

Put it in test/mechanism_definitions_for_test_enzyme_derivation.jl, because test_types.jl and test_rate_eq_derivation.jl also call it and run before the enumeration file. The allosteric-only checks are not carried over; the constructor or the field type already guarantees each one: tag vector length and validity, multiplicity ≥ 1, and `isa Vector{RegulatorySite}`.

**Other files and tests.** Tests:
- Rename 83 call sites (77 in the enumeration tests, 5 in test_types, 1 in test_rate_eq_derivation).
- Allosteric call sites become stricter, since the coverage check now applies to them too; they are expected to pass.
- test_types.jl:352-379 then tests a test helper. The tests group may reword or retarget it.

**Evidence.** Five proposals were confirmed: 0 src callers and 83 test callers. The check 'every substrate/product appears in some step' is the only non-tautological one, so it moves to the helper. #73's other parts (types.jl and identify) belong to other groups.

### enum/B1 — One DFS for flux blocks: merge `_edge_blocks`, `_unbalanced_blocks` and `_block_balanced`
- Net lines removed (this file group): 35
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** Give the Tarjan DFS in `_edge_blocks` a potential array: on a tree edge e from u to w, set `phi[w] = phi[u] + (edges[e][1] == u ? weights[e] : -weights[e])`. Return `block, Set(block[e] for (e, (u, v)) in enumerate(edges) if weights[e] != phi[v] - phi[u])` as the single `_unbalanced_blocks(nv, edges, weights)`.

Delete `_edge_blocks` (1345-1384), the `members` grouping and `_block_balanced` (1478-1507). Put the telescoping-defect proof in the docstring and reword the `_flux_carrying_steps` docstring, which names `_edge_blocks`.

**Other files and tests.** None (0 test or docs references). The flux testsets (test 7416-7760) and the depth-2 aggregate pins cover it.

**Evidence.** Confirmed. A non-tree edge's fundamental cycle lies in its own block, and the tree edges span each block, so a nonzero defect exists exactly when the block is unbalanced. Self-loops never reach the function. About 81 lines become about 45.

### enum/C1 — Small compactions in the flank, flux, vmax, hyperbolic and gauge helpers
- Net lines removed (this file group): 15
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** (a) `_chain_flank_groups`: replace `group_of` with `lone = Dict(only(group) => g for (g, group) in enumerate(groups) if length(group) == 1)`. `flank` returns the group index, and the loop runs over `lone`.
(b) Delete the method `_flux_carrying_groups(m)`. Its only caller, `_assert_emission_rules`, passes `steps(m), reaction(m)`.
(c) `_re_turnover_cycle(groups, rxn, fast = is_equilibrium)` gains a step filter. `_has_vmax` calls it with `s -> is_equilibrium(s) || touches(s)` instead of rebuilding Steps through `_with_equilibrium`.
(d) `_hyperbolic_catalysis`: write `name.(consumed(s))` and `name.(released(s))`, build `mets` with one comprehension, and drop `filter(!isempty, …)`.
(e) `_gauge_rescaling`: replace the `shared` tally with `count(==(t), values(twin_of)) == 1`.

**Other files and tests.** Tests: the three `_flux_carrying_groups(m)` calls at test_mechanism_enumeration.jl:7538, 7582 and 7613 pass `steps(m), reaction(m)`.

**Evidence.** Confirmed by the verifiers, about 14 + 4 lines. #24's items (e) and (f), on the split move, are done in SP. Item (c) also makes `_has_vmax` cheaper, which matters for F.

### enum/R1 — `_re_segment_extras` returns the form index and segment map; delete `_indexed_re_segments` and `_group_re_segments`
- Net lines removed (this file group): 28
- Risk: low
- Behaviour change: none (segment ids are used only for intersection tests)
- Depends on: —

**Change.** `_re_segment_extras` (types.jl:725-768) also returns `index`, the idx Dict it already builds, and `segment_of`, filled during the BFS.
- Delete `_indexed_re_segments` (1393-1404) and call `_re_segment_extras` at 1430, 1835 and 2205.
- Inline `_group_re_segments` (2064-2074) into `_expand_split_kinetic_group`: `_, _, _, idx, seg = _re_segment_extras(groups)`, then per group `Set{Int}(seg[idx[from_species(s)]] for s in group)` for an RE group and an empty Set for an SS group. Keep its reason as a comment. This drops the `_state_mechanism(m, :A)` construction and the per-step `findfirst`.
- In `_chemistry_equilibrates_both_sides`, replace the `node` BFS closure (1619-1629) with the segments of `_re_segment_extras([re_iso])`. A form outside every RE isomerization is a node of its own.

**Other files and tests.** types.jl: +3 lines and a docstring update. `_bottomless_re_segment` destructures 3 of the 5 values, which Julia allows. The derivation group's RE-segment unification should call the 5-tuple `_re_segment_extras` instead of `_indexed_re_segments`. After R1 and R2 this file no longer calls `_compute_re_groups` or `_state_mechanism`. Tests: none (0 references).

**Evidence.** All three were confirmed: 9, 6 and 9 lines. `_state_mechanism(am, :A)` has the steps of `am`, and uniform groups make the first step representative of its group.

### enum/R2 — One minimal-flip search over groups (`_minimal_flips`) and a groups-level `_re_segment_count`
- Net lines removed (this file group): 66
- Risk: medium
- Behaviour change: none (same children in the same order)
- Depends on: —

**Change.** Add `_minimal_flips(admissible, base, rxn, eligible = _ -> true)`, about 20 lines with its docstring. It precomputes `steady = _all_steady_state(base)` and the flux units (RE, flux-carrying on `steady`, eligible). It builds each candidate by reference from `steady`. It accepts a candidate when `admissible(gs) && _bottomless_re_segment(gs) === nothing && all(_flux_carrying_groups(gs, rxn)[units[sel]])`, and returns the flipped groups of every minimal set.

Other changes:
- Replace both `_re_segment_count` methods and `_re_segment_count_after_flip` (1895-1929) with `_re_segment_count(groups) = length(_re_segment_extras(groups)[2])`.
- `_expand_re_to_ss` (1270-1295) shrinks to about 9 lines. The admissible test is a segment-count gain on the groups, the eligible test is the existing flank and regulator test, and the hyperbolic filter is unchanged. It also loses the Core.Box created by reassigning the captured `flux`.
- `_seed_variants`: build `bases = [[merged]; filter(!isnothing, [_eliminate_form(merged, x) for x in complexes])]`, use `enumerate` in place of the is_merged flag, and run each base through `_minimal_flips(screen, base, rxn)`. `_seed_candidate_screen` takes the candidate's groups and reads its RE flags from them, which deletes the `group` index array and the mask plumbing.
- Delete `_flip_group_to_ss` (1332-1338).

The order of tests inside each predicate and the unit and set order are unchanged.

Perf gate: time the `_expand_re_to_ss` testsets (12.1 s and 14.1 s) and the ter-ter expand budget test (16.8 s). If calling `_re_segment_extras` per candidate is measurably slower, make `_re_segment_count(groups)` a Dict-indexed union-find of about 10 lines instead.

**Other files and tests.** Tests:
- The `_flip_group_to_ss` uses (test 2754, 9417, 9435, and `_testhelper_flip_groups` at 9454-9460, still used at 10255) become an inline comprehension.
- About 11 `_re_segment_count(m)` calls become `_re_segment_count(steps(m))`.
- The screen agreement test (8270) passes groups instead of a mask.
- Delete the testset '_re_segment_count_after_flip agrees with the built child' (9453-9491). It tests only the deleted function and takes 0.1 s; Denis signs off.

**Evidence.** The verifiers estimated 34 and 32. The flip move needs the groups-level count. Precomputing `steady` avoids rebuilding Steps per candidate, the cost `_re_segment_count_after_flip` existed to avoid. The verifiers flagged the per-candidate cost of `_re_segment_extras` as needing measurement. The derivation-side parts of #63/#74/#90 belong to the derivation group.

### enum/F — Delete `_seed_candidate_screen` and give seed candidates the group predicates (perf-gated)
- Net lines removed (this file group): 78
- Risk: high
- Behaviour change: No output change (the agreement test shows identical verdicts on 1,456 candidates). Only speed changes: a rejected candidate now pays for the Dict/DFS predicates instead of array passes.
- Depends on: R2
- Approved by Denis on the sign-off page.

**Change.** Split `_degenerate` into `_degenerate(groups, rxn)` and `_degenerate(m) = _degenerate(steps(m), reaction(m))`. `_seed_variants` then calls `_minimal_flips(gs -> !_degenerate(gs, rxn), base, rxn)`.
- Delete the screen (1652-1724).
- Inline `_fused_substrate_binding` (1591-1594) and `_crosses_sides` (1596-1600) into `_chemistry_equilibrates_both_sides`, now their only caller.
- Reword the `_seed_variants` docstring (1735-1736) and developer.md:151-153.

If F is rejected, the fallback is #84: one weighted union-find shared by the screen and `_re_turnover_cycle`, saving about 8 lines.

**Other files and tests.** docs/src/developer.md:151-153. Tests: delete '_seed_candidate_screen agrees with the group predicates on every candidate' (test_mechanism_enumeration.jl:8207-8316, about 110 lines, 1.1 s). It exists only to keep the two copies in agreement.

**Evidence.** Both verifiers marked it needs_measurement. They confirmed the equivalence line by line: a base holds no isomerization and every base step is RE. Nothing in the repo records how much time the screen saves, and ter-ter init_mechanisms runs in the default test suite.

### enum/SP — Split-move compactions
- Net lines removed (this file group): 34
- Risk: low
- Behaviour change: none (same bipartitions, order and origin)
- Depends on: —

**Change.** - `_context_form` becomes one ternary.
- `_context_bipartitions`: build ligands with a `unique!` comprehension filtered by `b != own`. Replace the seen/out loop with a comprehension plus `unique!(first, out)`.
- `_bipartitioned_groups`: two comprehensions over `get(Dict(selection), g, (group,))`.
- One `_apply_bipartitions(m::Union{Mechanism, AllostericMechanism}, selection)` method.
- `_revert_zero_flux_parts`: `map` over `bp` with a Step→flag Dict. The caller tests value equality, `parts != bp`, and the RE-group early return goes.
- `selection(sel)` becomes `units[sel]`.
- `_split_gain_test(::Mechanism)`: replace the `next_id` bookkeeping with `enumerate(sel)`.
- At 2003, `kinetic_groups(m)` becomes `eachindex(groups)`, which decouples this file from the types group deleting that accessor.

**Other files and tests.** Tests: test_mechanism_enumeration.jl:9741 changes `=== bp` to `== bp`. The `_context_bipartitions`/`_context_form` testsets pin the output.

**Evidence.** Confirmed by the verifiers at 14, 12 and 4 lines, plus #24's (e) and (f) at 4. Julia 1.9 compat covers `unique!(f, A)`. `_with_equilibrium(s, true)` returns an RE step unchanged, so value equality holds.

### enum/RC — One `_redundant_copy_groups` method on a shared inactive-state graph
- Net lines removed (this file group): 9
- Risk: medium
- Behaviour change: none (both graphs are 'both ends reachable from free E', proven equal)
- Depends on: —

**Change.** Merge the two methods. They share the competitive-inhibitor guard and `_productive_twin(active)`; after that prefix, `m isa Mechanism && return [g for g in eachindex(active) if _gauge_rescaling(active, g, twin, identity) !== nothing]`.

Replace the hand-built inactive graph (2346-2351) with `_inactive_groups(am)`, a documented helper in rate_eq_derivation.jl. It returns groups aligned with `steps(am)`: Step[] for :OnlyA groups, and otherwise the steps with both ends in `_reachable_from_free`. `_state_allo_mechanism` uses the same helper, then builds from `keep = findall(!isempty, groups)` and the matching tags.

**Other files and tests.** rate_eq_derivation.jl: +8 lines for `_inactive_groups` and about -10 for `_state_allo_mechanism`'s pruning. This is the codegen-time derivation path, so rerun the allosteric golden, collapse and ground-truth files. Optionally, in the derivation group, rewrite `_reachable_from_free`'s fixpoint as one BFS. Tests: none.

**Evidence.** Both were confirmed: 'not stranded' equals 'in reach', because reach is a subset of all forms. The change removes a definition that could drift between the copy rule and the derivation.

### enum/D1 — Refactor `_expand_add_dead_end_regulator`: hoist per-regulator invariants, fuse the mirror loops, delete `_add_competitive_inhibitor` and `exclude_regs`
- Net lines removed (this file group): 95
- Risk: low
- Behaviour change: none (same children in the same order: regulator-major, then pattern order)
- Depends on: —

**Change.** Restructure the move (2465-2596, 132 lines, down to about 52):
- `bound_regs`: the consumed Regulator names, the same predicate as today's `existing_regs`. This keeps behaviour identical even for a DSL fixture whose catalytic step consumes an AllostericRegulator.
- `eligible_regs`: a comprehension over the name-sorted `regulators(rxn)` that excludes `bound_regs ∪ _bound_allo_regs(m)`, with no `sort!`.
- `form_sp`: a comprehension.
- Placements are computed once per parent: each competition pattern's sorted active forms, deduped in pattern order, with `eligible_forms` folded into the filter.
- `for reg_name in eligible_regs, active in placements`: build the complex map; build the mirror groups as typed `Step[group; …]`; put the copy group last; call `_dead_end_child(m, groups, rxn, twin)`.

Also:
- Delete `_add_competitive_inhibitor` (2389-2402); it returns `rxn` unchanged on every src path.
- Delete the `exclude_regs` keyword, which no caller passes.
- Correct the false docstring paragraph at 2459-2463: the child's reaction is `rxn` itself.
- Keep the in-body comments.

**Other files and tests.** Tests:
- Delete the `_add_competitive_inhibitor` testset (test_mechanism_enumeration.jl:1280-1287) and the `exclude_regs` testset (4836-4864). Both test deleted code; Denis's rule is to remove unused features entirely.
- Replace `EnzymeRates._add_competitive_inhibitor(rxn, :I)` with `rxn` at 3388 and 11921.
- The exact-children dead-end tests (3375-3425, 11905-11945) are the gate.

**Evidence.** The rewrite (#38) was confirmed at 75. The `_add_competitive_inhibitor` identity was confirmed at 15, and `exclude_regs` has 0 src callers. `existing_inhibitors` equals the bound set because reg_name is never bound, and `name(bm) == reg_name` can never fire. The init-side part of #85 is subsumed by S4.

### enum/D2 — Dead-end children: build, then apply `_redundant_copy_groups`; delete `_dead_end_child` and `_all_twin`
- Net lines removed (this file group): 30
- Risk: medium
- Behaviour change: The move judges copy redundancy on the built canonical child, over all copy groups. That is the same test `_assert_emission_rules` applies to every parent. A child whose older copy group became redundant is now dropped; today it is emitted and fails when expanded. The verifier found no change over the 128,493 bi-bi depth-2 mechanisms.
- Depends on: D1
- Approved by Denis on the sign-off page.

**Change.** In D1's loop, build the child: `m isa Mechanism ? Mechanism(rxn, groups) : AllostericMechanism(rxn, groups, [cat_allo_states(m); :EqualAI], catalytic_multiplicity(m), regulatory_sites(m))`. Keep it if and only if `isempty(_redundant_copy_groups(child))`.
- Delete both `_dead_end_child` methods and their docstring (2404-2429), and the `twin` precompute.
- Delete `_all_twin` (2224-2233) and inline its test as the first line of `_gauge_rescaling`.
- Keep `rxn` for the child, not reaction(m).
- Reword the docstrings at 2241, 2310 and 2453-2456.

**Other files and tests.** docs/src/developer.md:117. Tests: 8651 and 8694 call `_all_twin`; switch them to the inline predicate. Gate: time the dead-end testsets (38.8 s) and the depth-2 aggregate (70.5 s).

**Evidence.** Confirmed at 33. The verifier found that the deleted docstring's claim ('the placement cannot make an older one redundant') is unproven. Today the Mechanism and allosteric paths test different things.

### enum/A1 — Delete the four no-op wrong-kind move methods and `_add_expansions_mech!`
- Net lines removed (this file group): 40
- Risk: low
- Behaviour change: None for the package. Calling a move on the wrong kind now raises a MethodError instead of returning []; only tests do that.
- Depends on: —

**Change.** Delete these methods and their docstrings:
- `_expand_to_allosteric(::AllostericMechanism)` (2701-2709)
- `_expand_add_allosteric_regulator(::Mechanism)` (2821-2830)
- `_expand_change_allo_state(::Mechanism)` (2939-2945)
- `_expand_merge_regulatory_sites(::Mechanism)` (3026-3033)
- `_add_expansions_mech!` (3153-3164)

expand_mechanisms appends re_to_ss, split and dead_end, then `m isa Mechanism ? to_allosteric : (add_regulator, change_allo_state, merge_sites)` with multi-argument `append!`. Per-parent child order is unchanged. `_seed_children` branches the same way. Reword 'Mechanism-native overload' in the docstrings.

**Other files and tests.** Tests:
- Delete the no-op assertions at test_mechanism_enumeration.jl:5152-5162, 5907-5915, 6331-6333 and the Mechanism half of 6531-6533. Keep 5324, 5383 and 5524, which call the real moves.
- Rewrite 7083-7088 (`_add_expansions_mech!`) and raw6 at 7140-7148 to call the applicable moves.

**Evidence.** Confirmed by three verifiers: 46 lines replaced by about 6. Julia 1.9 supports multi-argument `append!`. Per #68's verifier, 3 of its 6 listed test lines test real moves, so they stay.

### enum/A2 — Inline the seed-BFS helpers `_seed_children`, `_has_nonequalai` and `_binds_all_required`
- Net lines removed (this file group): 26
- Risk: low
- Behaviour change: none
- Depends on: A1

**Change.** - Inline `_seed_children` (3278-3297) into the `pmap` closure in seed_mechanisms, with A1's kind branch.
- Fold `_has_nonequalai` (3299-3304) and the single-ligand check into one `m isa AllostericMechanism && (...) && return false` in `_is_seed_node`. `only(allo_states(site))` is safe once a site holds one ligand.
- Inline `_binds_all_required` (3341-3344) into `consider!`.
- Keep `_is_seed_node`, `_bound_allo_regs` and `_bound_comp_inhibitors`. The seed_mechanisms docstring already documents the child moves.

**Other files and tests.** Tests: rewrite the 'seed_mechanisms wave-parallel equivalence' serial reference (test_mechanism_enumeration.jl:12450-12485) with the inline expressions. At 12296, replace `_has_nonequalai` with a direct tag check.

**Evidence.** Confirmed at 25: each helper has one src caller, and the folded predicate is equivalent.

### enum/A3 — Compact `_expand_add_allosteric_regulator` and `_make_am_with_added_reg`
- Net lines removed (this file group): 40
- Risk: low
- Behaviour change: none (children in the same order: regulators, then tags, then sites)
- Depends on: —

**Change.** In `_expand_add_allosteric_regulator`:
- Take the bound names from `_bound_allo_regs(am)`.
- Iterate the name-sorted `regulators(rxn)` directly, with `reg isa AllostericRegulator && !(name(reg) in …) || continue`. This drops `new_regs`, `sort!` and the early return.
- Write the nested loops as `for tag in (:OnlyA, :OnlyI, :NonequalAI), site_idx in 0:length(sites)`.
- Write the :EqualAI pass as one loop.

`_make_am_with_added_reg` (2790-2819, 30 lines) becomes about 10: copy the sites, then either push a new single-ligand site or replace site i with `[ligands(s); lig]` and `[allo_states(s); tag]`, then `_with_reg_sites(am, sites)` (or the types group's keyword `_with`).

Correct the provably false claims at 2724-2726 and 2770-2772. RegulatorySite has no all-:EqualAI rule; the move's own guard prevents it. Also correct the rxn paragraph at 2732-2734. Do not name a local `bound`, which would shadow `bound(::Species)`.

**Other files and tests.** types.jl: the `_with_*` to keyword-`_with` consolidation from #65/#91 is counted in the types group; this file only uses whichever helper exists. Tests: none. The add-regulator exact-children tests (5427-6034) are the gate.

**Evidence.** Confirmed at 38: regulators are name-sorted (types.jl:396), and the constructor does not alias its inputs. #79's other parts are distributed: (c) to D1, (e) to A6, (f) to A4. Its `_regulator_names` helper is unnecessary once each site is a one-line comprehension. (b) is skipped because the types group may delete `allosteric_regulators(m)`.

### enum/A4 — One relaxation loop in `_expand_change_allo_state`
- Net lines removed (this file group): 19
- Risk: low
- Behaviour change: none (same candidates in the same order)
- Depends on: —

**Change.** Build the relaxation sets first:
- `[g]` for each non-chemistry group that is not :NonequalAI;
- then `findall(chem)` once, if any chemistry group is not :NonequalAI.

One loop copies the tags, sets the set to :NonequalAI with `new_states[gs] .= :NonequalAI`, applies `_onlya_haldane_violation` and `_partial_onlya_catalysis`, and pushes `_with_cat_allo_states`. The regulatory-ligand loop becomes `for (si, site) in enumerate(regulatory_sites(am)), li in eachindex(ligands(site))`. Keep the chemistry-relaxes-together comment.

**Other files and tests.** None. The exact-children tests at 6055-6333 and 12554-12716 are unchanged.

**Evidence.** Confirmed at 18. The two blocks have identical bodies. Keep the 5-line in-body comment.

### enum/A5 — Inline `_partial_onlya_catalysis` as a condition on the chemistry relaxation
- Net lines removed (this file group): 16
- Risk: low
- Behaviour change: A hand-built partial parent (some :OnlyA tag plus a live chemistry group) passed straight to expand_mechanisms would now also get its partial binding-relaxation children. Enumeration never builds such a parent.
- Depends on: A4
- Approved by Denis on the sign-off page.

**Change.** Delete `_partial_onlya_catalysis` (2832-2850). Enqueue the all-chemistry relaxation only when no binding group is :OnlyA; on that branch this is exactly the partial test. Drop the check on binding singletons; for every enumerated parent it is vacuous there, because every enumerated parent is non-partial (by induction over all seven moves). Move the rationale into the move docstring, which already states it at 2867-2879, and update the mention at types.jl:227.

**Other files and tests.** types.jl:227 (docstring mention). Tests: comment wording at 11291 and 12589.

**Evidence.** No verifier checked this proposal. I checked it myself. After the chemistry relaxation every chemistry group is :NonequalAI, so 'partial' holds exactly when an :OnlyA tag remains. For non-partial parents a binding relaxation cannot create a partial mechanism.

### enum/A6 — Compact `_expand_to_allosteric`
- Net lines removed (this file group): 18
- Risk: low
- Behaviour change: None: same children in the same order. The Haldane check is no longer repeated for each multiplicity.
- Depends on: S0

**Change.** - Make `chem` a Bool vector, with `vtags = [c ? :OnlyA : :EqualAI for c in chem]`.
- Build `regs` with one comprehension over the name-sorted regulators, without `sort!`.
- Compute the K-type tag vectors once, before the multiplicity loop: `for sel in _subsets(findall(!, chem))[2:end]`, with the Haldane prefilter on `reaction(m)`. That matches the constructor's own check and equals rxn on every src path.
- For each multiplicity cn, append the K-types, then build each V-type directly: `AllostericMechanism(reaction(m), steps(m), vtags, cn, [RegulatorySite([AllostericRegulator(reg)], cn, [tag])])`. This removes the throwaway `am_cat`.
- Drop the no-op `unique!` and the redundant `copy(steps(m))`.
- Update 'duplicate variants are removed' (2642).

**Other files and tests.** None.

**Evidence.** Confirmed at 15. The prefilter does not depend on cn, the constructor copies the tags, and the multiplicities are distinct. #58: the V-type site multiplicity is cn either way.

### enum/A7 — Inline `_merged_site_state_assignments` into the merge move
- Net lines removed (this file group): 16
- Risk: low
- Behaviour change: none (same order: all-keep first, then retags in ligand order)
- Depends on: —

**Change.** Inline the helper (3035-3057) into `_expand_merge_regulatory_sites`:
- push the all-keep assignment unless `redundant`;
- then push each assignment that retags exactly one non-:EqualAI ligand to :EqualAI, unless the result is all-:EqualAI.

The merge docstring already documents both kinds, and the `drop_all_keep` keyword disappears with the helper.

**Other files and tests.** Tests: delete the helper's unit testset (test_mechanism_enumeration.jl:6352-6365). The merge-move exact-children tests (7106-7136, `count(is_merged, children) == 2`) subsume it; Denis signs off.

**Evidence.** Confirmed at 16: one src caller, and the keyword exists only for that call.

### enum/A8 — One reg-type predicate: fold `_declared_reg_type`, `_state_respects_reg_type` and `_filter_by_reg_type` into `_respects_reg_type`
- Net lines removed (this file group): 40
- Risk: medium
- Behaviour change: none
- Depends on: —

**Change.** In `_respects_reg_type(am, rxn)`:
- Build `declared = Dict(name(regulator(rm)) => reg_type(rm) for rm in regulators(rxn) if regulator(rm) isa AllostericRegulator)` once.
- Per site, `types = [get(declared, name(l), :unspecified) for l in ligands(site)]`.
- A ligand passes if its type is :unspecified, or if none of these holds: activator and :OnlyI; inhibitor and :OnlyA; :EqualAI with `count(==(t), types) > 1`.

Merge the three docstrings into one, keeping the activator, inhibitor and antagonist rules. expand_mechanisms uses `filter!(c -> _respects_reg_type(c, rxn), result)`. Reword the mentions at 3001 and 3226. This replaces about 67 lines (3059-3125) with about 25.

**Other files and tests.** Tests:
- Rewrite the `_declared_reg_type`/`_state_respects_reg_type` unit block (test_mechanism_enumeration.jl:6576-~6730, about 28 asserts) as a `_respects_reg_type` table over one- and two-ligand sites, with no coverage loss.
- 11 `_filter_by_reg_type` calls become `filter(c -> EnzymeRates._respects_reg_type(c, rxn), xs)`.

**Evidence.** #46 was confirmed at 12; #56, which subsumes it, was not checked by a verifier. Checked: the early return guarantees t is not :unspecified at the :EqualAI test, so `count(==(t), types) > 1` means 'another ligand at the site has the same type'. Risk is medium only because of the test rewrite.

### enum/A9 — Canonical ligand order inside a regulatory site
- Net lines removed (this file group): 3
- Risk: medium
- Behaviour change: Allosteric enumeration emits fewer duplicate children, so allosteric child counts can drop. DSL sites written out of name order get reordered, which can reorder entries in parameters(m) and fitted_params (the names do not change) and changes their Sig.
- Depends on: —
- Approved by Denis on the sign-off page.

**Change.** The RegulatorySite inner constructor (types.jl:127-140) sorts the ligands by name and applies the same permutation to allo_states. Delete the explicit sortperm in `_expand_merge_regulatory_sites` (3010-3012).

Today `_make_am_with_added_reg` appends the new ligand. So 'add Y at X's site' gives [X, Y], and 'add X at Y's site' or a merge gives [Y, X] or [X, Y]. These are the same mechanism with different ==/hash; both pass the beam's structural dedup and are compiled, and fit too unless eq_hash catches them.

TDD: a failing test that the add-then-merge routes give == mechanisms.

**Other files and tests.** types.jl +2 lines; the types group input lists this as #5, so count types.jl lines there. Sites in DSL fixtures written out of name order get reordered. Expect possible re-pins of allosteric golden or parameter-order tests, and a changed Sig for such inputs.

**Evidence.** Confirmed. RegulatorySite ==/hash compare ligands by position, and only the merge move sorts them. Denis's memory already lists 'within-site canon pre-existing' as an open flag. How many duplicates a beam run produces has not been measured.

## src/identify_rate_equation.jl (1147) + src/fitting.jl (277)
Excess code in this group comes from six sources.
(1) _beam_search writes the compile → prefit log → fit → CSV → ingest → postfit log pipeline twice (base tier and iteration, about 90 lines). The compile/fit split hands a 6-tuple between two single-use functions, plus a test-only _process_batch, a duplicated `orig` field and mixed skip sentinels.
(2) Single-use helpers and wrappers: _select_beam and _parsimony_cutoff around _select_count!, _expand_parent, _subset_data, _evaluate_loss, _scatter_fold_scores, and two CSV one-liners.
(3) Validation is written twice in the two problem constructors, with hand-written typed-constructor and `sk` boilerplate.
(4) Some invariants are enforced twice: the LOOCV eq_hash re-dedup repeats what _offer_cv! already guarantees, and the pre_best/improved recomputation repeats what _ingest! already knows.
(5) Redundant data shapes: parallel fitted_param_names/values tuples and a hand-typed parent NamedTuple.
(6) Unused public knobs: lb/ub, abstol/reltol/callback, and the Mechanism overload of FittingProblem.
About 183 lines go with no observable change and no decision. About 39 more need only Denis's sign-off on deleting tests. The remaining ~109 change log text, the CSV schema or public keywords and are flagged.
Robustness items: reject non-finite Rate and Keq ≤ 0, keep Inf-loss fits out of LOOCV, and validate n_cv_candidates ≥ 1.
Test time: merge duplicate end-to-end runs and kcat fits, and derive the 240 bi-bi equations once instead of twice. Also add an allosteric zero-allocation loss! test that gates the types group's @generated metabolites(::AllostericEnzymeMechanism) fix, likely the largest speed win for the allosteric identify tests.
loss! stays allocation-free throughout. No load-bearing invariant (canonical form, rate_equation perf, name chokepoint, Pass-2 rename) lives in this group.

### identify/V1 — One shared rate-table validator for FittingProblem and IdentifyRateEquationProblem; get! grouping; drop the sk/Float64 boilerplate
- Net lines removed (this file group): 48
- Risk: low
- Behaviour change: none. Same checks in the same order, same messages, same group order.
- Depends on: —

**Change.** Add `_rate_table(table, mnames, scale_k_to_kcat)` to fitting.jl: a 5-line docstring and about 14 lines of code. It runs the same checks in the same order with byte-identical messages:
- scale_k_to_kcat is positive or nothing;
- `Tables.columntable`;
- the :group/:Rate columns exist ('Missing required column: $req');
- every mnames column exists ('Missing metabolite column: $m');
- `i = findfirst(iszero, data.Rate)` finds no zero ('Zero rate at row $i: log(0) is undefined').
It returns `data`.

FittingProblem (fitting.jl:46-94, 49 lines → ~12):
- `data = _rate_table(table, metabolites(mechanism), scale_k_to_kcat)`;
- the group map becomes `for (i, g) in enumerate(data.group); push!(get!(() -> Int[], group_map, g), i); end`. Keep the Dict so group_point_indexes order is unchanged;
- keep the typed call `FittingProblem{typeof(mechanism), typeof(data)}(mechanism, data, collect(values(group_map)), Keq, scale_k_to_kcat, log.(abs.(data.Rate)), Vector{Float64}(undef, length(data.Rate)))`. The inner constructor converts Keq and scale (Base's convert for T>:Nothing), so the `sk` line 89 and the Float64(Keq) call go, along with the log_abs_rates and buffer temporaries.

IdentifyRateEquationProblem constructor (identify:32-77, 46 lines → ~14): a 2-line mnames tuple, then `_rate_table`, the existing ≥2-groups check, and the typed call without `sk`.

Do NOT switch to the default outer constructor: it requires the exact field types and would raise a MethodError on `Keq=2` (Int) or a Float32 Rate column.

Net: fitting.jl about -17, identify about -32.

**Other files and tests.** No src lines elsewhere. Tests unchanged: test_fitting.jl:203-209 and 505-521 and test_identify_rate_equation.jl:86-144 assert ErrorException only.

**Evidence.** Read fitting.jl:46-94 and identify_rate_equation.jl:32-77: the same five checks with identical messages. The verifier of proposal 11 showed the default outer constructor breaks Int Keq; the inner typed constructor converts. I kept Dict grouping instead of the findall/unique variant so loss summation order stays bit-identical.

### identify/V2 — Drop the useless R type parameter of IdentifyRateEquationProblem
- Net lines removed (this file group): 2
- Risk: low
- Behaviour change: The exported type's parameter list changes, so the printed type of a problem object changes from IdentifyRateEquationProblem{EnzymeReaction, NamedTuple{…}} to IdentifyRateEquationProblem{NamedTuple{…}}.
- Depends on: V1
- Approved by Denis on the sign-off page.

**Change.** Change identify_rate_equation.jl:23-25 to `struct IdentifyRateEquationProblem{D<:NamedTuple}` (one line). Change the docstring header at line 10 from `{R, D}` to `{D}`, and the constructor call to `IdentifyRateEquationProblem{typeof(data)}(reaction, data, Keq, scale_k_to_kcat)`. EnzymeReaction is a concrete, non-parametric type (types.jl:387), so R carries no information.

**Other files and tests.** None. `IdentifyRateEquationProblem{` appears nowhere in test/ or docs/src.

**Evidence.** `grep IdentifyRateEquationProblem{` finds only the definition and the constructor call at identify:10/23/74. EnzymeReaction is non-parametric (types.jl:387).

### identify/V3 — Reject non-finite/missing Rate and Keq ≤ 0 at problem construction
- Net lines removed (this file group): -3
- Risk: low
- Behaviour change: Tables with NaN, ±Inf or missing Rate values and Keq ≤ 0 now error at construction. Today:
- a NaN rate is silently scored as a sign mismatch, and its sentinel shifts the group mean;
- an Inf rate makes every centered loss NaN, so every fit ends :NoFiniteLoss;
- a missing rate throws a cryptic 'non-boolean (Missing)' error;
- Keq ≤ 0 yields NaN or Inf Haldane constants.
- Depends on: V1
- Approved by Denis on the sign-off page.

**Change.** Add two checks to `_rate_table`:
- `i = findfirst(r -> !(r isa Real && isfinite(r)), data.Rate); i === nothing || error("Rate at row $i must be a finite number; got $(data.Rate[i])")`, placed before the zero check so the zero message stays unchanged;
- a `Keq` argument with `Keq > 0 || error("Keq must be positive; got $Keq")`.
Both constructors pass Keq. Optionally apply the same finite check to the metabolite columns.

**Other files and tests.** Add @test_throws cases for NaN, Inf and missing Rate and for Keq=0 to test_fitting.jl 'Validation errors' (505-521) and test_identify 'construction' (86-144).

**Evidence.** fitting.jl:66-67 checks only `== 0`. In loss!, sign(NaN) never equals pred's sign, so the 10.0 sentinel is written (141-143), and an Inf in buf gives NaN after centering (165-168). Keq is never validated anywhere.

### identify/V4 — IdentifyRateEquationProblem uses the shared reaction metabolite-name helper; add an allosteric zero-alloc loss! test
- Net lines removed (this file group): 2
- Risk: low
- Behaviour change: none in values. The identify column check now dedups names, which only affects which duplicate name is reported first. It may make allosteric loss! allocation-free.
- Depends on: V1

**Change.** Once the types group adds `_metabolite_names(rxn)` (substrates, products, then regulators, made unique) and makes `metabolites(::AllostericEnzymeMechanism)` @generated, replace the 2-line mnames tuple in the IdentifyRateEquationProblem constructor with `_rate_table(table, _metabolite_names(reaction), …)`. In test_fitting.jl, first add a failing `@allocated loss!(x, fp) == 0` case for an allosteric FittingProblem (TDD gate). fitting.jl:126 keeps calling `metabolites(fp.mechanism)` unchanged.

**Other files and tests.** The src change lives in types.jl (helper plus @generated AEM metabolites; owned by the types group), rate_eq_derivation.jl:195-208 (_concentration_symbols). test_fitting.jl gains about 8 lines in 'Zero allocations'. If the type instability is confirmed, allosteric fits in test_identify (main pipeline, 'mechanism recovery' refit of 3×10 s) should speed up.

**Evidence.** loss! builds `NamedTuple{metabolites(fp.mechanism)}` per point (fitting.jl:126-139). The AEM method (types.jl:1712-1722) splats a runtime Vector. The only zero-alloc loss! test uses uni_uni (test_fitting.jl:271-290). Verifiers rated it very likely real but unmeasured.

### identify/L1 — loss!: count sign mismatches in pass 1; drop the redundant `pred == 0.0` test and the sentinel-rescan loop
- Net lines removed (this file group): 2
- Risk: low
- Behaviour change: Only a correctly signed point whose log-ratio is exactly 10.0 is no longer miscounted as a mismatch (measure-zero bug fix). All other loss values are bit-identical, and loss! stays allocation-free.
- Depends on: —

**Change.** In loss! (fitting.jl:123-184):
- add `n_mismatch = 0` before pass 1;
- write the branch as `if sign(pred) != meas_sign` with `buf[i] = 10.0; n_mismatch += 1`. `|| pred == 0.0` is redundant: meas_sign is ±1 for validated nonzero rates, and sign(±0.0) and sign(NaN) never equal it;
- delete the third loop (178-181) that recounts by `buf[i] == 10.0`;
- keep the comment 173-177 above the return and adjust the docstring wording to 'counted in the per-point loop and added after'.
Keep the explicit uncentered loop (154-157). `sum(abs2, buf)` would save 2 more lines but changes absolute-mode losses in the last ULP (pairwise summation). Optional, 0 lines: drop the unused `where {M,D}`.

**Other files and tests.** None. 'Zero allocations' (test_fitting.jl:271-290), 'Speed' (293-347) and the exact-100.0 penalty tests (221-268) must stay green.

**Evidence.** fitting.jl:142-146 writes the sentinel and 178-181 rescans it with float equality. Zero rates are rejected at construction. An Int counter does not allocate.

### identify/FR1 — fit_rate_equation: compute pnames once, pass loss! directly, index the rescaled NamedTuple
- Net lines removed (this file group): 5
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** In fitting.jl:222-277:
- put `pnames = fitted_params(fp.mechanism); np = length(pnames)` at the top and delete line 265;
- write `Optimization.OptimizationFunction(loss!)` instead of the `(x, p) -> loss!(x, p)` lambda at 233;
- replace the rescale block 267-275 with `result_params = rescale_parameter_values(fp.mechanism, merge(result_params, (Keq = fp.Keq, E_total = 1.0)); scale_k_to_kcat = fp.scale_k_to_kcat)[pnames]`, inside the existing `if`.
rescale_parameter_values returns `NamedTuple{keys(params)}` (rate_eq_derivation.jl:1006), and NamedTuple indexing by a Symbol tuple exists from Julia 1.7 (compat is 1.9).

**Other files and tests.** None. 'fit_rate_equation kcat rescaling' asserts `keys(res.params) == fitted_params(uni_uni)`.

**Evidence.** Read fitting.jl:222-277. Verifier confirmed iip defaults and key order.

### identify/C1 — Delete the two CSV writer one-liners and compact the save_dir guard
- Net lines removed (this file group): 17
- Risk: low
- Behaviour change: none (same files, same guard message)
- Depends on: —

**Change.** Delete `_save_initial_csv` and `_save_iteration_csv` (identify:255-266). Call `_write_rows_csv(save_dir, "initial_mechanisms.csv", rows)` at 844 and 851, and `_write_rows_csv(save_dir, "equation_search_iteration_$(iteration).csv", rows)` at 923. Move the 'iteration is a 1-based sequential counter, not a parameter count' note into the existing gap-free comment at 901-902.

Compact the save_dir guard (222-230, 9 lines) to `isdir(save_dir) && any(f -> endswith(f, ".csv"), readdir(save_dir)) && error(<same message>)`, about 3 lines.

Optional, 0 lines: `isdir(save_dir) || mkpath(save_dir)` becomes `mkpath(save_dir)` at 251, 452 and 1086, since mkpath is idempotent.

**Other files and tests.** Rewrite test 'csv writers' (test_identify 500-521) to call _write_rows_csv; its dir-creation sub-case stays valid because _write_rows_csv still creates the directory.

**Evidence.** _save_initial_csv is called at 844/851 and _save_iteration_csv at 923; the only test uses are at 509, 511 and 518. The guard is 9 lines of one-word-per-line formatting.

### identify/C2 — save_dir guard also refuses a stale progress.log
- Net lines removed (this file group): 0
- Risk: low
- Behaviour change: A save_dir containing only progress.log, left by a run that crashed before its first CSV, is now rejected. Today a new run silently appends its log to the old one.
- Depends on: C1
- Approved by Denis on the sign-off page.

**Change.** Extend the compacted guard predicate to `f -> endswith(f, ".csv") || f == "progress.log"` and reword the message to 'save_dir already contains results (CSV files or progress.log). Use an empty directory to avoid mixing results.'

**Other files and tests.** Add one test: a directory holding only progress.log is rejected. 'save_dir non-empty check' (test_identify 404-420) still passes.

**Evidence.** The guard (identify:222-230) checks only .csv files, while _progress appends to progress.log (452-455).

### identify/R1 — Rows carry `params::NamedTuple` instead of parallel fitted_param_names/fitted_param_values tuples
- Net lines removed (this file group): 15
- Risk: low
- Behaviour change: none. CSV and cv_results columns and values are unchanged. Mechanisms sharing an eq_hash share the fitted-param tuple and order, because the hashed text includes the `(; indep…, Keq, E_total) = params` line (rate_eq_derivation.jl:677), and test 1301 asserts it.
- Depends on: —

**Change.** - _fit_batch row (identify:625-626): `params = fit.params`.
- _failure_row (436-437): `params = (;)`.
- Drop `fitted_param_names = fkeys` from the PASS-1 record (546). `n_params = length(fkeys)` stays.
- _rows_to_dataframe (273-307): replace the 7-line Set loop with `sorted_pnames = sort!(unique(k for r in rows for k in keys(r.params)))`, and the O(P²·rows) findfirst loop with `df[!, pn] = [get(r.params, pn, missing) for r in rows]`.
Keep the explicit fixed-column DataFrame(...) construction; see dropped items.

**Other files and tests.** In test_identify_rate_equation.jl, about 9 hand-built row literals change to `params = (a = 1.0, b = 2.0)` or `params = (;)`: lines 148-159, 179-189, 501-507, 888-890, 1082-1085, 1821-1824, 1840-1844, 1869-1873. Assertions at 935 and 1329-1334 read `e.row.params`.

**Evidence.** Read identify:273-307, 424-440 and 615-629. fitted_param_values is rebuilt from a NamedTuple already in hand, then turned back into columns by findfirst.

### identify/B1 — _ingest! returns the set of counts whose best loss improved
- Net lines removed (this file group): 5
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** In _ingest! (identify:657-669), keep `improved = Set{Int}()`, push `e.n_params` where the best is set (661-663), and return `improved`. Both call sites become `improved = _ingest!(...)`: the base tier at 854-858 and the iteration at 927-932, which also drops the `!isempty(child_entries) &&` guard because ingesting an empty vector is a no-op. Extend the docstring by one clause.

**Other files and tests.** '_ingest! and cv pool' (test 1077-1102) can additionally assert the returned set.

**Evidence.** The same pre_best/improved comprehension is written at 854-858 and 927-932. _ingest! makes the exact comparison at 661-662. Fit losses are never NaN.

### identify/B2 — Merge _select_beam, _select_count! and _parsimony_cutoff into one _select_count!
- Net lines removed (this file group): 45
- Risk: low
- Behaviour change: none (the only src call always passed best_override = best_loss_by_count[c])
- Depends on: —

**Change.** Replace _select_beam (identify:327-371), _select_count! (373-393) and _parsimony_cutoff (961-971) with one function:
`_select_count!(expanded::Dict{Int,Int}, best_loss_by_count::Dict{Int,Float64}, c::Int, losses; loss_rel_threshold, loss_abs_threshold, loss_parsimony_threshold, min_beam_width)`.
Its body:
- `cutoff = loss_rel_threshold * best_loss_by_count[c] + loss_abs_threshold`;
- `smaller = [l for (k, l) in best_loss_by_count if k < c]`, and if it is non-empty, `cutoff = min(cutoff, loss_parsimony_threshold * minimum(smaller))`;
- `budget = max(0, min_beam_width - get(expanded, c, 0))`;
- `perm` = finite indices stably sorted by loss;
- `sel = sort!([i for (rank, i) in enumerate(perm) if losses[i] <= cutoff || rank <= budget])`;
- update `expanded[c]` and return `sel`.
Merge the docstrings and keep all their content: the cutoff rule, min over ALL smaller counts, the cumulative budget spent once, non-finite losses excluded, input-order return. The call at 876-881 shrinks to about 3 lines.

**Other files and tests.** Rewrite test_identify testsets 523-569, 766-777, 779-808, 810-835 and 837-843 (18 _select_beam, 4 _select_count! and 4 _parsimony_cutoff calls) to call the merged function with an explicit `best_loss_by_count` and a fresh Dict. Every asserted case stays expressible: the floor, the abs term, Inf/NaN exclusion, input order, the best override, parsimony tightening and min over all smaller counts.

**Evidence.** grep: _select_beam is called only from _select_count!, and both _select_count! and _parsimony_cutoff only from _beam_search:876-881. Verifier checked the sketch line by line: 83 old lines become about 40.

### identify/B3 — Sweep frontier buckets directly instead of pooling into `swept` and regrouping
- Net lines removed (this file group): 7
- Risk: low
- Behaviour change: none. Sorting the counts would give a deterministic ascending order but would change parent_of attribution and CSV child order, so it is not proposed.
- Depends on: B2

**Change.** Replace identify:864-883 and 950-951 with:
- `target = 0` before the loop;
- `target = max(target + 1, minimum(keys(frontier)))` at the top of the loop;
- then `for c in [k for k in keys(frontier) if k <= target]; entries_at_count = pop!(frontier, c); sel = _select_count!(…); append!(to_expand, entries_at_count[sel]); end`.
Iterate the keys in Dict order and do NOT sort. That reproduces today's order exactly, because `unique(e.n_params for e in swept)` followed the Dict order of popped buckets and each bucket keeps its insertion order. n_params ≥ 1, so the first target equals the minimum key.

**Other files and tests.** none

**Evidence.** Buckets hold exactly one count (`get!(frontier, e.n_params, …)` at 660). I read 864-883 and confirmed that the unsorted key iteration preserves the current order.

### identify/B4 — _expand_parents: inline _expand_parent into the pmap; parent_of maps child → parent BatchEntry
- Net lines removed (this file group): 10
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** Delete _expand_parent (identify:692-706). Write the pmap at 716-717 as a do-block that runs the same try/catch and returns `(children, nothing)` or `(Union{…}[], FitFailure(pe.mech, _exc_string(e)))`, under a 2-line comment that keeps the docstring's rationale. Make `parent_of = Dict{Union{Mechanism, AllostericMechanism}, BatchEntry}()` with `parent_of[child] = pe` (718-719, 726-727). The row builder reads `parent.n_params` and `parent.row.mechanism_type`. parent_of stays on the main process, so nothing extra is serialized.

**Other files and tests.** Test changes:
- '_expand_parent records an expansion error' (test 1581-1617): call `_expand_parents([BatchEntry(m_bad, 0, 0.0, :Success, UInt64(0), (mechanism_type="M",))], rxn_bad)` and assert on expand_failures.
- '_expand_parents parallel equivalence' (1619-1657) compares against a hand-copied serial loop built on _expand_parent. Rewrite it to assert directly against `expand_mechanisms` output: children in first-parent order, and parent_of[child] === the first parent producing it. This keeps coverage without a copy of the code under test.

**Evidence.** _expand_parent has one src caller (716). parent_of is read only at 614 (`get(parent_of, c.orig, nothing)`).

### identify/B5 — Compress _required_regulators
- Net lines removed (this file group): 5
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** Replace the two near-identical setdiff blocks (identify:776-784) with a local `required(T, optional) = setdiff(Set{Symbol}(name(regulator(rm)) for rm in regulators(rxn) if regulator(rm) isa T), optional)` and return `(required(AllostericRegulator, optional_allosteric_regulators), required(CompetitiveInhibitor, optional_competitive_inhibitors))`. If the enumeration group adds `_regulator_names(rxn, T)`, use `Set(_regulator_names(rxn, T))` instead. Keep the init/seed branch at 822-825 (see dropped items).

**Other files and tests.** None. '_required_regulators' (test 1174-1210) keeps its exact set assertions.

**Evidence.** The two blocks at 776-783 differ only in the type and the optional list. The single caller is at 819.

### identify/B6 — One batch runner: merge _compile_batch + _fit_batch (+ test-only _process_batch) into _process_batch, and run base tier and iterations through one closure
- Net lines removed (this file group): 75
- Risk: medium
- Behaviour change: Search results, CSV contents and file numbering are unchanged. Only progress-log text changes:
(1) the base post-fit line loses its 'Base tier: ' prefix;
(2) the iteration header loses '(child n_params lo-hi)' and prints the tentative iteration number;
(3) an all-skipped expansion prints the standard header plus the prefit line instead of 'Expanded P parents → C children | all skipped (…)';
(4) a non-empty to_expand that yields no children now prints a header and an all-zero prefit line;
(5) in the all-base-failed path, a postfit line prints before the error.
- Depends on: C1, R1, B1, B4
- Approved by Denis on the sign-off page.

**Change.** (a) Merge _compile_batch, _fit_batch and _process_batch (identify:501-648) into one `_process_batch(mechs, prob; optimizer, max_param_count, eq_complexity_filter, memo, seen, parent_of, log = msg -> nothing, kwargs...)`. It returns the current 5-tuple `(entries, failures, n_param_skip, n_cx_skip, n_seen_skip)`. Changes inside it:
- the 13-line seen loop becomes 4 lines and keeps the in-batch duplicate drop;
- the PASS-1 sentinels become `:param_skip` and `:complexity_skip`, and the row loop skips with `c isa Symbol && continue`;
- the duplicate `orig = m` field is dropped; use `c.mech` at 608 and 614;
- rep_idx is built with `get!`;
- `log(_prefit_summary(...))` is called before PASS 2;
- KEEP the `reps = [(mech, eq_hash) …]` vector so each pmap task does not serialize the whole `compiled` batch;
- rename `fitted` to `seen` here and in _beam_search, since it holds every produced structure, not only fitted ones;
- merge the two docstrings, keeping all their content.
(b) In _beam_search, add `progress(msg)` and a closure `run_batch!(mechs, extra_failures, header, csv; parent_of = Dict())`. It calls _process_batch with `log = msg -> progress(header * "\n  " * msg)`, appends extra_failures, returns `nothing` when there are no rows, and otherwise writes the CSV, calls `improved = _ingest!(…)` and prints the postfit and best-loss lines. The closure captures only variables that are never reassigned, so there is no Core.Box.
(c) The base tier becomes `run_batch!(base, base_expand_failures, "Fitting N initial mechanisms…", "initial_mechanisms.csv")`. If it returns `nothing`, return empty (as today). Keep the 'Every base-tier fit failed …' error verbatim.
(d) Each iteration becomes `run_batch!(children, expand_failures, "Iteration $(iteration + 1): P parents → C children", "equation_search_iteration_$(iteration + 1).csv"; parent_of) === nothing || (iteration += 1)`. The special 'all skipped (…)' branch (936-947) is deleted, because the prefit summary already prints the same three skip counts.
If Denis says no: keep the two functions, extract only the shared post-fit tail into a closure with byte-identical text (about -20), and move _process_batch to test/ as `_testhelper_process_batch` (src -13). Total about -33.

**Other files and tests.** test_identify:
- the M2 regex at 1145-1147 changes to match the prefit summary, e.g. r"0 new fits \+ 0 inherited \+ \d+ skipped \(already fit\) \+ \d+ skipped \(>3 params\)";
- the one test passing `fitted` by name (971-981) passes `seen`;
- the 18 _process_batch call sites keep the same return shape;
- the 'CSV output (new schema)' log assertions (342-344) still hold.

**Evidence.** I read 501-648 and 816-952. The two blocks share every call in sequence. _process_batch has no src caller (18 test calls). `orig` is always `mech` (543/608/614). I drafted the merged function at about 102 lines against 146 today, and the beam at -45 after the closure. Verifiers gave -48 (unify) and -35 (merge) separately. Keeping `reps` follows the proposal 9 verifier's serialization warning.

### identify/B8 — Keep non-finite fits out of the frontier, best-loss map and LOOCV pool; validate n_cv_candidates ≥ 1
- Net lines removed (this file group): -2
- Risk: low
- Behaviour change: (a) An Inf-loss fit (:NoFiniteLoss) no longer enters LOOCV. Today, at a sparse count, it fills a pool slot, and its fold can abort model selection after the whole search. It also no longer shows as 'c:Inf' in the best-loss line. CSV rows are unchanged.
(b) n_cv_candidates=0 now errors immediately instead of failing after the search.
- Depends on: B2
- Approved by Denis on the sign-off page.

**Change.** - Add `isfinite(e.loss) || continue` at the top of the _ingest! loop.
- Add `n_cv_candidates >= 1 || error("n_cv_candidates must be ≥ 1; got $n_cv_candidates")` at the top of identify_rate_equation.
- Delete `n == 0 && return pool` in _offer_cv! (677).
- In the merged _select_count!, `perm = sortperm(losses)` replaces the finite filter (0 lines).
- Update the fitting.jl:247-248 comment ('dropped by the beam's non-finite filter' becomes '_ingest! drops it').

**Other files and tests.** Test changes:
- move the Inf/NaN cases of the selection tests (test 547-559) to an _ingest! test asserting that non-finite entries are skipped;
- replace the 'n=0 must not panic' assert (1099-1101) with an identify-boundary @test_throws, which is cheap because it raises before any fit.

**Evidence.** _offer_cv! pushes unconditionally while `length(pool) < n` (683-684). fit_rate_equation returns loss=Inf when no restart is finite (fitting.jl:245-263). _cv_fold_loss errors on a non-finite fold (1117-1119). Verifier corrected the saving to roughly 0 lines.

### identify/CV1 — Inline _subset_data and _evaluate_loss into _cv_fold_loss
- Net lines removed (this file group): 22
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** Delete _subset_data (identify:973-981) and _evaluate_loss (983-994). _cv_fold_loss (1099-1121) becomes:
- `held = prob.data.group .== held_out`;
- `fold(mask) = (idx = findall(mask); FittingProblem(mechanism, map(col -> view(col, idx), prob.data); Keq = prob.Keq, scale_k_to_kcat = prob.scale_k_to_kcat))`, with a 2-line comment keeping the 'views, so no fold copies the data (O(G·N·ncols))' rationale;
- `fit = fit_rate_equation(fold(.!held), optimizer; kwargs...)`;
- `test_loss = loss!([log(v) for v in fit.params], fold(held))`.
The isfinite guard and its comment are kept. fit.params is `NamedTuple{fitted_params(mechanism)}`, so iteration order matches what _evaluate_loss rebuilt.

**Other files and tests.** None. Tests call only _cv_fold_loss (test 422-496, 1884).

**Evidence.** grep: each helper has one caller and no test or docs use. Keys are preserved through the rescale step (fitting.jl:265-274).

### identify/CV2 — Delete the LOOCV eq_hash re-dedup in _cv_model_selection (the cv_pool already guarantees it)
- Net lines removed (this file group): 17
- Risk: low
- Behaviour change: none on every identify_rate_equation path: the candidate set and order are identical. A direct call with duplicate eq_hash rows, which only that test makes, would now cross-validate all of them.
- Depends on: —
- Approved by Denis on the sign-off page.

**Change.** Replace identify:1027-1044, the groupby/sort/seen_hashes loop, with `candidate_indices = sortperm(collect(zip(df.n_params, df.loss)))`. This keeps today's n_params-then-loss row order, which matches DataFrames' current grouping of Int columns and the stable per-group sort. Drop `n_cv_candidates` from _cv_model_selection's signature (1018) and from its call (242). Keep the comment explaining that candidates are distinct because _offer_cv! keeps the top n distinct eq_hash per count.

**Other files and tests.** Test changes in test_identify:
- delete the first half of 'LOOCV eq_hash-uniqueness guard (§4)' (1807-1834), which feeds two same-eq_hash rows straight into _cv_model_selection; its _offer_cv! half (1836-1851) stays;
- tests at 1826 and 1876 must drop `n_cv_candidates=5`, or it would flow into fit_rate_equation's kwargs and throw;
- the comment at test 287-292 should name _offer_cv!.

**Evidence.** cv_pool is written only via _offer_cv! (665), which keeps the top n distinct eq_hash per count. _beam_search returns exactly the pool (954-958), and _cv_model_selection has one src caller (240). The proposal 38 verifier required the sortperm form to keep row order.

### identify/CV3 — Replace _scatter_fold_scores with a reshape of the order-preserving pmap result
- Net lines removed (this file group): 22
- Risk: low
- Behaviour change: none (same columns, same values, same order)
- Depends on: —
- Approved by Denis on the sign-off page.

**Change.** In _cv_model_selection (identify:1050-1081):
- `scores = pmap([(m, g) for m in candidate_mechs for g in groups]) do (m, g); _cv_fold_loss(compile_mechanism(m), prob, g; optimizer, kwargs...); end`, with a comment that pmap keeps grid order so column ci holds candidate ci;
- `S = reshape(Float64.(scores), length(groups), :)`;
- `cv_df.cv_score = [mean(c) for c in eachcol(S)]` and `cv_df.cv_score_se = [std(c) / sqrt(length(groups)) for c in eachcol(S)]`. Use eachcol, not `mean(S; dims=1)`, so values stay bit-identical to today's per-vector mean/std;
- `cv_df[!, Symbol("cv_fold_$g")] = S[i, :]`.
Delete _scatter_fold_scores (1123-1135), the temporary :cv_fold_scores column (1063) and `select!(…, Not(:cv_fold_scores))` (1081). The final column order is unchanged.

**Other files and tests.** Test changes:
- delete '_scatter_fold_scores places each fold' (test 1890-1898);
- extend '_cv_model_selection flatten reproduces serial LOOCV' (1854-1888), which today has ONE candidate and so cannot detect a transposed reshape, to two candidates with distinct equations and fold scores. Give it the exotic group labels ("a=b", "c,d", "x y") and a CSV round-trip;
- replace the 'exotic group labels' testset (667-691), which exercises a hand-copied loop and calls no package code (CLAUDE.md violation).

**Evidence.** _scatter_fold_scores has one call (1059). pmap order is already relied on (_expand_parents docstring 709-713). The tasks vary g fastest. The flatten test uses one candidate, so I added the 2-candidate extension to keep placement coverage.

### identify/M1 — Remove the unused `using Random` from the module
- Net lines removed (this file group): 0
- Risk: low
- Behaviour change: none
- Depends on: —

**Change.** Delete `using Random` at src/EnzymeRates.jl:29 and move Random from [deps] to [extras] and [targets].test in Project.toml, keeping `Random = "1"` in [compat]. fitting.jl:255's `randn` resolves to Base.randn, whose methods Random defines in the sysimage.

**Other files and tests.** src/EnzymeRates.jl -1 line and Project.toml edits. Tests that `using Random` (runtests and 5 files) need it in the test target, or Aqua's stale-deps check and the test imports fail.

**Evidence.** grep src for Random./shuffle/MersenneTwister/seed!/rand( finds only fitting.jl:255 randn (Base).

### identify/T1 — Cheaper same-coverage test runs for fitting/identify (all assertions kept; time each before and after)
- Net lines removed (this file group): 0
- Risk: low
- Behaviour change: none (test-only)
- Depends on: —

**Change.** (a) test_fitting.jl: merge 'scale_k_to_kcat normalization' (350-391) and 'fit_rate_equation kcat rescaling' (394-424) into one testset with n_restarts=1 and maxtime≈1.0. Keep every assertion: kcat≈1 (default), ≈42, ≈7, keys, finite loss, retcode, and the `nothing` raw path. The anchoring is exact whatever the fit quality, because rescale_parameter_values divides by the current kcat. The old budget was up to about 62 s of maxtime.
(b) test_identify: fold 'identify runs on a solver that rejects popsize' (1104-1125), 'all-cap-skipped expansion batch is reported (M2)' (1127-1150) and 'loss_parsimony_threshold threads through' (1152-1172) into ONE run: max_param_count=3, loss_parsimony_threshold=2.0, show_progress=true. This saves two full searches plus LOOCV.
(c) Compute the ping-pong bi-bi `init_mechanisms` + `_base_tier` once at file level for 1659-1680 and 1746-1775. This is an aggregate seed-set pin, an allowed exception.
(d) Measure first: the 'mechanism recovery' refit (255-268) uses n_restarts=3 and maxtime=10; cut it to the smallest budget that passes reliably, especially after V4.
(e) Measure first: the main pipeline's max_param_count=8 against the 5-parameter generating mechanism; lower it only if the recovery assert still passes.

**Other files and tests.** Test-only: about 40-60 test lines removed and no assertion lost. If B6 lands, the merged run's M2 regex follows B6's text.

**Evidence.** The three uni-uni runs share the rxn, data and optimizer settings. The two kcat testsets assert the same properties with 5 BBO fits. Times are unmeasured (no Julia): wrap each top-level testset in @time during one focused run before acting.

### identify/T2 — Delete the subsumed '_cv_fold_loss: one fold, finite' testset
- Net lines removed (this file group): 0
- Risk: low
- Behaviour change: none (test-only)
- Depends on: —
- Approved by Denis on the sign-off page.

**Change.** Delete test_identify_rate_equation.jl:477-496. Its fixture, settings and assertions (Float64, ≥0, finite on group 2) are a strict subset of '_cv_fold_loss over all groups: per-fold scores, finite' (422-452), which scores groups 1, 2 and 3 with identical settings.

**Other files and tests.** Test-only: -20 test lines, one fewer pair of 2-restart × 2 s fits.

**Evidence.** Read both testsets: the same DataFrame, Keq=10, CMA-ES n_restarts=2, maxtime=2, maxiters=500, and a superset of groups in 422-452.

### identify/T3 — 'rate-eq dedup-key partition stability': derive each equation once
- Net lines removed (this file group): 0
- Risk: low
- Behaviour change: none (test-only)
- Depends on: —
- Approved by Denis on the sign-off page.

**Change.** In test_identify_rate_equation.jl:744-763, compute `h` once per mechanism. Today line 759 calls `rate_equation_string(em)` a second time for every one of the 240 mechanisms (1 uni-uni and 239 bi-bi), and rate_equation_string re-derives at runtime (rate_eq_derivation.jl:669-683). Keep the determinism re-check only for uni-uni and the first few bi-bi mechanisms. The partition-count pin (239 distinct classes) is kept in full.

**Other files and tests.** Test-only. Roughly halves this testset's time, likely one of the heaviest in the file (unmeasured).

**Evidence.** Lines 756 and 759 both call rate_equation_string(em). test_mechanism_enumeration.jl:1466-1471 derives only 5 bi-bi equations because 'full derivation is slow'.

## Tests

### Test time plan across the whole suite (1152 s of test time, target 576 s or less)

Denis, the plan below does not reach 576 s. All time figures are estimates: I did not run Julia. I took the 1152 s baseline and per-testset rows from the timing table. Where proposals conflicted I read the code: the hyperbolic sweep, division-freeness, depth-2, the compile-budget subprocesses, the rank oracle, the BBO fits, the derivation helpers, the Sig decode and to_allosteric.

Where the time stands, in three tiers:
- **Tier 1, confirmed or exact-by-construction: about 313 s.** Main items: TT01 drops the bi-bi branch of the hyperbolic sweep (90 s). TT02 checks division-freeness on the polynomials instead of compiling 264 types (60 s). TT03 removes the duplicate ter-ter run in compile_budget (50 s). Then BBO maxtime (25 s), mass-action replacing QSSA/ODE (18 s), depth-2 run once (13 s), rank-oracle seed (12 s), the ter-ter split run once (10 s) and the identify popsize run (9.5 s). This takes the suite to about 839 s.
- **Tier 2, coverage-neutral but each needs a stated measurement: +82 s.** The two big ones are the hyperbolic memo (30 s) and a compile-free rank oracle for plain mechanisms (25 s).
- **Tier 3, lowered thresholds or deletions with a stated coverage loss: +38 s.** The biggest are the BFS generation cap and the identify cap (needs measurement and is risky).
- **Test-side total: about 433 s, giving about 719 s (1.6x).**
- **Two src items, both unmeasured: +55 s, giving about 664 s.** One stops per-type JIT by not specializing the Sig decode and derivation entry points; rate_equation's generated body is untouched. The other is profile-guided enumeration speedups.

Overlaps resolved:
- **ter-ter runs.** Only one of the two can go. Keep the in-process 250,855 pin (CLAUDE.md puts init tests in the enumeration file) and drop the compile_budget subprocess run. Proposal #0, which does the reverse, is dropped.
- **Ordering.** Item times are incremental under a fixed order:
  - TT01 before TT04 and TT14 (the two hyperbolic memos).
  - TT07, then TT16, then TT22 (derivation spec tests).
  - TT09, then TT06, then TT32/TT33/TT34 (rank oracle).
  - TT02, TT09 and TT20 before TT35.

  The rank-oracle deletions (TT32/TT33/TT34) are worth about 0.3 s each once TT06 lands. If Denis rejects TT06, they are worth about 9 s together.
- **Compile cost that only moves is not counted.** For example, #168's renders reappear in the identify partition testset.

Every deletion and every lowered threshold is flagged for Denis's sign-off. So is every change that drops incidental codegen coverage: TT02 (about 257 init types) and TT06 (plain rank fixtures).

#### time/TT01 — Hyperbolic exactness sweep: drop the bi-bi iteration (the ping-pong-capable bi-bi population contains it)
- Lever: delete_redundant; time saved ≈ 90 s; test lines removed ≈ 15; risk low
- Sign-off: approved by Denis
- Measure first: Before deleting, in one focused Julia process (no other Julia running), enumerate both populations exactly as the testset does and assert issubset(Set(key.(bibi_mechs)), Set(key.(pp_mechs))), with key(m) = m isa Mechanism ? steps(m) : (steps(m), cat_allo_states(m), catalytic_multiplicity(m)). Record @elapsed of the bi-bi branch (enumeration plus derivation loop) to price the saving; the verifier's range is 80-120 s.

**Testsets.** test/test_mechanism_enumeration.jl:9319-9401 "_hyperbolic_catalysis matches the derived denominator" (271.3 s)

**Change.** In test_mechanism_enumeration.jl:9375 loop over ((unibi, 2), (pingpong, 2)). Delete the bibi reaction (9363-9367), _testhelper_allosteric_children (9349-9356) and the append line (9379-9380). Rewrite the header comment (9320-9327): the sequential seeds of the ping-pong-capable reaction are the bi-bi seeds, so its depth-2 population holds every bi-bi mechanism and their allosteric children.

**Coverage.** - Seeds: the 55 sequential topologies are identical under both reactions. Ping-pong paths carry residual isomerizations, so they never enter the sequential group. merged-3 and TC tallies match exactly (pins 1486-1505).
- Moves, the predicate and the derivation read metabolite names, never atoms. So ping-pong level 1 contains bi-bi level 1, and ping-pong level 2 = expand(level 1) contains the allosteric children of the level-1 allosteric mechanisms.
- Each bi-bi assertion is therefore repeated verbatim by a ping-pong assertion.
- n_checked > 100 and 0 < n_nonhyperbolic < n_checked stay meaningful.
- No compile moves; nothing here compiles a @generated function.

_Merged from:_ tests_all#47 (enum_hyper)

#### time/TT02 — init division-freeness: check the reduced denominator symbolically instead of compiling 264 rate equations
- Lever: cheaper_same_coverage; time saved ≈ 60 s; test lines removed ≈ 0; risk medium
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:12166-12184 "init division-freeness (bi_bi_pp)" (74.4 s)

**Change.** Rewrite test_mechanism_enumeration.jl:12166-12184. For each bi_bi_pp init mechanism:
- keep `compile_mechanism(m) isa EnzymeMechanism`;
- derive `_, den, _ = _raw_symbolic_rate_polys(m, _step_parameters(m), _build_wegscheider_rename_map(m))`;
- for each of A, B, P, Q, assert that some den monomial is free of it (1056 tests kept);
- add a closure check: every free symbol of num, den and each dependent RHS of `_dependent_param_exprs(m)` lies in indep ∪ {keys assigned earlier in sorted order} ∪ metabolites ∪ {:Keq}.

Keep the current compiled rate_equation check for the 7 ping-pong seeds only. Rewrite the 1/conc comment, which becomes false. This also removes the file's only use of random_reduced_params, the focused-run UndefVarError artifact.

**Coverage.** - The compiled body is E_total·num/den plus concentration-free dependent assignments built from the same _raw_symbolic_rate_polys.
- den has positive coefficients (matrix-tree cofactors), and exponents are non-negative after _reduce_conc_lowest_terms. So 'finite at x=0' is exactly 'den has an x-free monomial'.
- The closure check keeps the §5a UndefVarError class that the numeric call caught incidentally.
- Lost: LLVM codegen and execution on about 257 init types (7 kept). This is the same kind of cap Denis accepted at 1446-1471.
- No later test compiles these types, so nothing moves.
- The remaining cost is about 264 symbolic derivations (2-15 s), so 60 s is conservative.

_Merged from:_ tests_all#54 (enum_C); tests_all#157 (suite_level); src_enum#16

#### time/TT03 — compile_budget: drop the cold ter-ter warm-up subprocess and time warm uni-uni inside the bi-bi trace-compile subprocess
- Lever: cheaper_same_coverage; time saved ≈ 50 s; test lines removed ≈ 25; risk medium
- Sign-off: approved by Denis
- Measure first: One Julia process at a time: run the new traced script (bi-bi init, GC.gc(), timed uni-uni) and the unchanged cold uni-uni script 3 times each. Accept if max(UNI_WARM/UNI_COLD) <= 2e-3 and the trace count stays <= 100. If the ratio fails, keep a separate warm subprocess but warm it with bi_bi_pp (A[CX],B[N] -> P[C],Q[NX]) instead of ter-ter; that variant saves ~45 s.

**Testsets.** test/test_compile_budget.jl:92-116 'trace-compile: init_mechanisms (bi-bi)', 164-208 'compile reuse: ter-ter warms all of uni-uni' (file 78.2 s)

**Change.** - In the trace-compile script (test_compile_budget.jl:93-104), after init_mechanisms(bi-bi), build uni-uni with the direct EnzymeReaction constructor (no DSL, so no parser trace lines), call GC.gc(), then print UNI_WARM from @elapsed init_mechanisms(r_uni).
- Make the trace helper also return stdout, merging _count_relevant_precompiles with _measure_labeled_subprocess.
- Delete the warm ter-ter script and its subprocess (174-191), plus `@test isfinite(t_ter)` and `@test t_ter < 150` (196-200).
- Keep the cold uni-uni subprocess and the `warm/cold < 0.01` gate.
- Update the ABOUTME line, the 155-163 comment and the testset name.

Three subprocesses instead of four, and no ter-ter.

**Coverage.** - The ter-ter time ceiling stays enforced by test_mechanism_enumeration.jl:1568-1575: same structure, same 150 s budget, plus the 250,855 count pin. Only the cold-JIT share of the ter-ter call is lost, a few seconds against 150 s, and the trace-count gate already bounds it.
- The reuse gate keeps its assertion and threshold; only the warmer changes. EnzymeReaction is non-parametric, and dispatch identity stays asserted at 213-225.
- Warm uni-uni adds no trace lines if reuse holds. If reuse breaks, the trace gate only gets stricter.
- Rule-compliant: the init_mechanisms pin stays in the enumeration file. Do not also apply #0.

_Merged from:_ tests_all#158 (suite_level); tests_all#73 (derivation_tests); tests_all#166 (suite_level); tests_all#24 (enum_A)

#### time/TT04 — Hyperbolic sweep: derive each allosteric state once per distinct (state, state Mechanism)
- Lever: cheaper_same_coverage; time saved ≈ 30 s; test lines removed ≈ -4; risk low
- Sign-off: not required (coverage unchanged)
- Measure first: Proposal #50, extended by its verifier. One focused process, nothing else running: per reaction, @elapsed the enumeration (_testhelper_levels) and the derivation loop separately. Print length(plain), length(allo), the distinct _state_mechanism(m,:A) and (m,:I) counts, and how many A-keys are in Set(plain). Also time _expand_split_kinetic_group and _expand_to_allosteric on the level-1 parents. If enumeration dominates, this item is worth less than 30 s.

**Testsets.** test/test_mechanism_enumeration.jl:9381-9395 (allosteric branch of the hyperbolic sweep)

**Change.** In the allosteric branch (9386-9391) memoize:
```
derived = Dict{Tuple{Symbol, EnzymeRates.Mechanism}, Bool}()
state_hyp(m, s) = get!(derived, (s, EnzymeRates._state_mechanism(m, s))) do
    poly_hyperbolic(EnzymeRates._state_rate_polys(m, s)[2], mets)
end
```
Keep `structural = _hyperbolic_catalysis(m)` and both @test lines for every mechanism.

**Coverage.** - Two mechanisms with equal _state_mechanism(m,s) give _state_rate_polys the same graph, step and group order, RE/SS flags and grouping. Tags change parameter names only.
- The pivot choice is name-blind (is_i_state is never passed; _step_priority is structural).
- den's concentration support is a positive matrix-tree sum, so it is independent of symbol names. The numerator shift is decided by den.
- The predicate still runs and is asserted on every mechanism. Each distinct derivation input still goes through the production _state_rate_polys, including _assert_derivable. No coverage is lost.

_Merged from:_ tests_all#48 (enum_hyper); tests_all#50 (enum_hyper)

#### time/TT05 — BBO kcat-anchoring fits: merge the two testsets and cut every fit to maxtime 0.5 s
- Lever: lower_threshold; time saved ≈ 25 s; test lines removed ≈ 37; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_fitting.jl:349-424 (inside 'Fitting', 43.4 s)

**Change.** Replace test_fitting.jl 'scale_k_to_kcat normalization' (350-391) and 'fit_rate_equation kcat rescaling' (393-424) with one testset. Build the data once and run three BBO fits, all at maxtime=0.5:
- default target 1.0 with n_restarts=2;
- scale_k_to_kcat=7.0 with n_restarts=1;
- scale_k_to_kcat=nothing with n_restarts=1.

Keep every assertion: keys(params) == fitted_params, isfinite(loss), retcode isa Symbol, _kcat_forward ≈ target (rtol 0.01), and nothing-mode keys. Keep fit_capturing_convergence.

**Coverage.** - rescale_parameter_values (src/rate_eq_derivation.jl:1000-1010) scales every SS constant by target/kcat_current. uni_uni's kcat is degree-1 homogeneous in those constants, so kcat ≈ target holds for any finite fitted point. No assertion depends on fit quality.
- Targets 42.0 and 7.0 run the same path, so one custom target suffices.
- Multi-restart selection stays exercised by n_restarts=2.
- MaxTime maps to ReturnCode.MaxTime without a warning, so the output stays pristine.
- BBO/Optimization first-use compile stays here, so nothing moves.

_Merged from:_ tests_all#138 (identify_fit_tests); tests_all#159 (suite_level); tests_all#167 (suite_level); src_enum#60(a); src_identify_fit#31(a)

#### time/TT06 — Rank and degeneracy oracles: evaluate the derived polynomials for plain Mechanisms instead of compiling rate_equation
- Lever: cheaper_same_coverage; time saved ≈ 25 s; test lines removed ≈ -15; risk medium
- Sign-off: approved by Denis
- Measure first: Warm focused process: for about 10 representative plain fixtures (uni-uni, ordered bi-bi with dead-end copies, the ter-ter split cube, ping-pong), time the current oracle against the polynomial version and assert both give identical ranks and verdicts. Count plain vs allosteric oracle calls in the file. Then run the whole enumeration file once to confirm every (fitted, rank) pin.

**Testsets.** test/test_mechanism_enumeration.jl helpers 92-160. 117 call sites, concentrated in _expand_add_dead_end_regulator (38.8 s), _expand_split_kinetic_group (context bipartitions) (20.4 s), _expand_re_to_ss (group-set flips) (14.1 s), _redundant_copy_groups (12.8 s), _expand_re_to_ss (12.1 s) and _seed_variants (about 8.5 s)

**Change.** In _testhelper_identifiable_rank (92-122) and _testhelper_degenerate (124-160), handle `m isa Mechanism` without compiling:
1. Take `num, den, _ = _raw_symbolic_rate_polys(m, _step_parameters(m), _build_wegscheider_rename_map(m))` and `dep, indep = _dependent_param_exprs(m)`.
2. Evaluate the dependent constants with a small recursive Expr evaluator (*, /, ^, +, - over Symbols and numbers).
3. Evaluate num/den from monomials pre-indexed into Vector{Float64} slots, using indep as the θ coordinates.

Keep the compiled path for AllostericMechanism, so the MWC combination is not re-implemented in test code. The helper grows by about 15 lines.

**Coverage.** - The asserted quantity, the numerical rank of ∂v/∂logθ (or the degeneracy verdict), is a property of the reduced rational function. That function is exactly E_total·num/den with the dependent constants substituted, the same objects _build_rate_body compiles (thermodynamic_constr:765-776).
- Float64 roundoff differences (~1e-16) are far below the 1e-7 SVD cut.
- Lost: incidental LLVM codegen of each plain oracle fixture. Codegen stays guarded by the 43-spec perf, Expr-shape and flat-string tests, the 25 fused/TC mass-action cases and TT02's 7 compiled seeds.
- The saving is incremental to TT09. It shrinks TT32, TT33 and TT34 to about 0.3 s each.

_Merged from:_ tests_all#178 (suite_level); tests_all#71 (enum_C); tests_all#26 (enum_A, superseded for plain fixtures)

#### time/TT08 — Depth-2 catalytic aggregate: enumerate each reaction once and split the levels by holds_iso
- Lever: cheaper_same_coverage; time saved ≈ 13 s; test lines removed ≈ 0; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_mechanism_enumeration.jl:11169-11240 'catalytic moves on bi-bi to depth 2' (70.5 s)

**Change.** Replace the 5 levels() calls (test_mechanism_enumeration.jl:11210-11212, 11235-11236) with `bibi = levels(rxn)` and `copies = levels(rxn6)`.
- Inside levels, accumulate one Bool for `holds_iso(c) == holds_iso(m)` over every (parent, child) pair and @test it once. Do not emit 135k tests.
- Assert `[count(holds_iso, l) for l in bibi] == [62, 349, 1134]`, the !holds_iso counts == [202, 669, 1237], and `[count(holds_iso, l) for l in copies] == [62, 1629, 31382]`.
- Add the rxn6 non-iso pin [202, 5109, 97111] for free.

#55's halves variant with an isdisjoint assertion is equivalent.

**Coverage.** - `keep` filters only level 0, and moves never add or remove an isomerization; the new accumulated assertion now checks that invariant explicitly.
- So an iso mechanism's BFS level in the full run equals its level in the iso-only run, and every pinned number is reproduced exactly.
- The obeys_rules sweep still covers every mechanism. It stays an aggregate pin, already labelled as Denis's rule allows.
- Removed work: one full rxn growth and rxn6's iso half, about 24% of rxn6's expansions.

_Merged from:_ tests_all#163 (suite_level); tests_all#55 (enum_C); src_enum#34

#### time/TT09 — Rank oracle: seed the RNG without a full rate_equation_string derivation
- Lever: cheaper_same_coverage; time saved ≈ 12 s; test lines removed ≈ 0; risk low
- Sign-off: not required (coverage unchanged)
- Measure first: In one focused run, accumulate @elapsed rate_equation_string(em) inside the helper and print the sum (that is the saving). Switch the seed, rerun the whole file once and confirm no rank pin moves at the 1e-7 cut.

**Testsets.** test/test_mechanism_enumeration.jl:92-122 _testhelper_identifiable_rank (102 call sites)

**Change.** At test_mechanism_enumeration.jl:99 replace `MersenneTwister(hash(rate_equation_string(em)) % 2^31)` with `MersenneTwister(1)`, or with `MersenneTwister(hash(string(typeof(em))) % 2^31)`. Reword the comment: a fixed or Sig-derived seed reproduces a failure in a later session.

**Coverage.** - The seed only picks draws. A phantom is an exact rank deficiency and full rank is generic over 3×60 points, so any deterministic seed tests the same property.
- With a fixed seed, paired fixtures with the same fp and mets still get identical draws.
- Lost: incidental rate_equation_string rendering on about 200 oracle shapes. rate_equation_string stays tested on the 43 specs, the allosteric golden file and the 239 bi-bi init mechanisms (identify:751-757).
- Each call today pays a per-type rate_equation_string specialization plus an uncached derivation.

_Merged from:_ tests_all#1 (enum_A); tests_all#30 (enum_B); tests_all#70 (enum_C); tests_all#164 (suite_level)

#### time/TT10 — bi-uni+R reachability BFS: cap at 3 generations instead of 14
- Lever: lower_threshold; time saved ≈ 11 s; test lines removed ≈ 0; risk low
- Sign-off: approved by Denis
- Measure first: Fresh process: print length(frontier) and @elapsed per generation of the current loop, then choose the largest cap that fits. The 11 s assumes the frontier grows combinatorially past generation 3.

**Testsets.** test/test_mechanism_enumeration.jl:7154-7192 'expand_mechanisms reaches ≥2 distinct-metabolite catalytic :OnlyA' (~15 of the 16.6 s 'expand_mechanisms' row)

**Change.** In test_mechanism_enumeration.jl:7176 change `gen < 14` to `gen < 3`; leave the rest unchanged. The cap value follows the measurement.

**Coverage.** - The explicit assertion is met at generation 1, because _expand_to_allosteric emits every binding subset.
- Three generations keep the only fuzz of the regulator moves on a bi-substrate reaction. That covers add_allosteric_regulator on K- and V-types and change_allo_state, RE→SS and split on R-bearing allosteric bi-uni mechanisms, under expand_mechanisms' emission and atom asserts and the constructor validators.
- Generations 4-13 only re-apply the same moves to larger mechanisms. That fuzz is what is lost.
- Chosen over #28 (stop at generation 1) and #44 (delete the testset), which drop that fuzz entirely.

_Merged from:_ tests_all#43 (enum_B); tests_all#28 (enum_B); tests_all#44 (enum_B)

#### time/TT11 — Ter-ter 55-step worst-case seed: run the split once and keep both budgets
- Lever: cheaper_same_coverage; time saved ≈ 10 s; test lines removed ≈ 33; risk medium
- Sign-off: approved by Denis
- Measure first: Warm focused run: @elapsed _expand_split_kinetic_group(worst) gives the saving. Confirm length(rest) == 69 once.

**Testsets.** test/test_mechanism_enumeration.jl:11007-11045 (nested in the 20.4 s split testset) and 11128-11167 (16.8 s)

**Change.** Merge the two testsets into one with one copy of the seed literal:
```
t_split = @elapsed split = _expand_split_kinetic_group(worst)
t_rest = @elapsed rest = vcat(_expand_re_to_ss(worst), _expand_add_dead_end_regulator(worst, terter), _expand_to_allosteric(worst, terter), _expand_add_allosteric_regulator(worst, terter), _expand_change_allo_state(worst), _expand_merge_regulatory_sites(worst))
@test length(split) == 12
@test t_split < 60
@test length(rest) == 69
@test t_split + t_rest < 120
```

**Coverage.** - Both counts stay pinned: 12 splits, plus 69 = 6 flips + 63 K-types (81 in all), and both budgets stay.
- expand_mechanisms re-runs the split internally with no cache, so today the 55-step split runs twice.
- expand_mechanisms' composition (filtered union of the moves) stays pinned at 7151. Its parent emission assert and child atom asserts run on every expand_mechanisms call elsewhere.
- Lost: a wall-clock bound on expand_mechanisms' own overhead on this seed, which is negligible.

_Merged from:_ tests_all#56 (enum_C)

#### time/TT12 — Main allosteric identify run: max_param_count 8 → 6
- Lever: lower_threshold; time saved ≈ 10 s; test lines removed ≈ 0; risk medium
- Sign-off: approved by Denis
- Measure first: Focused run of the file with the change, then read progress.log. Confirm all of:
- iterations 1 and 2 are written;
- an 'all skipped (…>6 params…)' line appears;
- results.best and at least one LOOCV candidate are AllostericEnzymeMechanism;
- recovery loss < 0.01.
Record the testset time.

**Testsets.** test/test_identify_rate_equation.jl:237-253 main run (testset 29.8 s), assertions 92/360/415

**Change.** test_identify_rate_equation.jl:249 and 415 max_param_count=8→6; 360 `<=(8)`→`<=(6)`. Add an assertion that at least one fitted CSV row is allosteric. Do not change n_groups.

**Coverage.** - Every assertion keeps its meaning: gap-free iteration CSVs 1..2, the 'skipped (>' line, parent columns, cv_results invariants and 1-SE selection, and the CSV roundtrip.
- Lost: fits of seven distinct 7-8-parameter types, two LOOCV candidates, and the mixed fitted/cap-skipped batch.
- Risk: results.best and the 'mechanism recovery' refit may stop being allosteric. The file header promises an allosteric path, so adopt only if the measurement shows it survives.

_Merged from:_ tests_all#141 (identify_fit_tests); src_enum#60(h)

#### time/TT13 — Identify: delete the uni-uni popsize run and fold loss_parsimony_threshold into the all-cap-skipped run
- Lever: delete_redundant; time saved ≈ 9.5 s; test lines removed ≈ 42; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_identify_rate_equation.jl:1104-1172 (9.8 s + 0.1 s)

**Change.** Delete 'identify runs on a solver that rejects popsize' (test_identify_rate_equation.jl:1104-1125) and 'loss_parsimony_threshold threads through identify_rate_equation' (1152-1172). Add `loss_parsimony_threshold=2.0` to the M2 run at 1139-1143, with a one-line comment.

**Coverage.** - Popsize contract (runs on CMA with default solver_kwargs): subsumed by the main run (245-253) and its results-type assertion (271-272). An injected option would make every base fit throw.
- Parsimony: its only claim is that the keyword is accepted. The M2 run proves that, and its assertions are unchanged (cutoff nothing at the base count, all children cap-skipped).
- Both must go together: the two runs are identical, so deleting one alone moves its 9.8 s of per-type compile onto the other. No later testset uses the uni-uni 4-6 parameter types.

_Merged from:_ tests_all#140 (identify_fit_tests); tests_all#169 (suite_level); src_identify_fit#33; src_enum#60(c)

#### time/TT14 — Hyperbolic sweep: reuse the denominator verdict by flattened step graph (plain and A-state)
- Lever: cheaper_same_coverage; time saved ≈ 7 s; test lines removed ≈ -6; risk medium
- Sign-off: approved by Denis
- Measure first: Same run as TT04's measurement: count distinct flattened-step graphs against n_checked and against the exact-Mechanism key counts.

**Testsets.** test/test_mechanism_enumeration.jl:9381-9395

**Change.** Applied after TT01 and TT04.
- Plain branch: keep `_build_wegscheider_rename_map(m)` per mechanism (the kernel smoke check), but memoize the hyperbolic verdict from `_raw_symbolic_rate_polys` by `Set(Iterators.flatten(steps(m)))`.
- Allosteric branch: look up the A-state verdict by `graph(_state_mechanism(m, :A))` in the same dict.
- `structural` stays per mechanism.

**Coverage.** - den = Σσ_g·D_g with positive matrix-tree terms, so its concentration-exponent support depends only on forms, RE/SS flags and consumed/released metabolites. Grouping and renaming only merge symbols, and nothing cancels.
- Numerator monomials dominate den monomials, so the lowest-terms shift and the verdict are functions of the graph.
- A grouping-dependent predicate bug is still caught, because `structural` is compared on every mechanism against a provably grouping-invariant verdict.
- Lost: repeated smoke executions of _raw_symbolic_rate_polys and _state_rate_polys(:A) on duplicate graphs.

_Merged from:_ tests_all#53 (enum_hyper); tests_all#49 (enum_hyper)

#### time/TT15 — 'reversing written steps changes nothing': assert compiled-type identity instead of re-deriving strings
- Lever: cheaper_same_coverage; time saved ≈ 6 s; test lines removed ≈ 0; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_rate_eq_derivation.jl:1395-1439 (7.5 s)

**Change.** Keep `@test m2 == m` and `am2 == am`. Replace both per-variant rate_equation_string comparisons with `@test ER.compile_mechanism(m2) === em` (and `=== spec.mechanism`). Keep one rate_equation_string(em) render per non-spec mechanism (the 20 fused cases, tc_* and two_in), the only place that renders them. Replace the tc_backwards string check with `tc_backwards === _testhelper_tc_cases.tc_ss`.

**Coverage.** - rate_equation_string(::M, mode) is a deterministic function of the singleton type M, so `===` implies equal strings. It is also stronger: it would catch a Sig difference that does not change the printed equation.
- Keeping about 25 single renders preserves the only renderer coverage of the fused and two-metabolite shapes.
- Spec renders stay covered by test_rate_equation_string.

_Merged from:_ tests_all#76 (derivation_tests); tests_all#179 (suite_level)

#### time/TT16 — Derivation spec helpers: one derivation per spec instead of three per random draw
- Lever: cheaper_same_coverage; time saved ≈ 5 s; test lines removed ≈ 10; risk low
- Sign-off: not required (coverage unchanged)
- Measure first: #88: in one focused process, include the definitions and the helper section of test_rate_eq_derivation.jl. Per spec, time warm _dependent_param_exprs(typeof(spec.mechanism)), and time each check function over two passes (compile+run vs run). On the first run, confirm the allosteric dep exprs contain only :call nodes.

**Testsets.** test/test_rate_eq_derivation.jl helpers 516-599 (used by every per-spec oracle in 'Enzyme Derivation Tests')

**Change.** - `_get_independent_params(m) = EnzymeRates.fitted_params(m)` (@generated, cached per type).
- Compute `dep, _ = _dependent_param_exprs(typeof(m))` once per spec in run_all_tests and pass it down.
- Evaluate dependent constants with a small Expr walker over :*, :/, :^ instead of eval(Meta.parse("let …")).
- Hoist analytical_oracle_params' Mechanism/AllostericMechanism decompile per spec.
- Delete _get_dependent_params and _eval_dep_expr.

**Coverage.** - Test plumbing only. The same expressions are evaluated at the same parameter values, so the oracles get identical inputs and no assertion changes.
- The package derivation stays exercised uncached by test_structure, test_constraint_counting, rate_equation_string and every @generated body.
- The time is incremental after TT07 and TT22 (draws fall from ~1806 to ~277). Standalone it would be 10-20 s.

_Merged from:_ tests_all#162 (suite_level); tests_all#75 (derivation_tests)

#### time/TT17 — Identify small cuts: cap-stop the random bi-bi, fold the original-mechanism check, assert the beam's recorded loss
- Lever: cheaper_same_coverage; time saved ≈ 5 s; test lines removed ≈ 20; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_identify_rate_equation.jl 1004-1048 (3.3 s), 1290-1352 + 1532-1579 (2.0 s), 255-268

**Change.** Three changes:
1. #143: in 'eq_complexity_filter' (1004-1048), call _process_batch with max_param_count=0. Assert `_eq_complexity(Mechanism(random_bibi)) == 336`, n_cx_skip == 1, n_param_skip == 1, and that entries and failures are empty. Keep cx_tight == 2.
2. #144: add `@test [f.mech for f in bad_failures] == [m1, m2]` after 1351. Delete the LDH-split half (1564-1578).
3. #145: replace the 'mechanism recovery' refit (255-268) with `best_row.loss < 0.01` on the _select_best_row row.

**Coverage.** - The complexity check precedes the parameter cap, so n_param_skip == 1 proves the 336-term bi-bi passed the default filter. The old assertion never asserted a successful fit.
- Fit-stage failures are built from c.orig whether the mechanism is a representative or a duplicate, and pmap preserves order, so the new check is a superset of the deleted LDH one.
- cv_results.loss comes from the same FittingProblem and fit_rate_equation path on the same data; baseline losses are ~1e-13. FittingProblem on an allosteric singleton stays exercised by _cv_fold_loss.

_Merged from:_ tests_all#143; tests_all#144; tests_all#145; src_enum#60(b)

#### time/TT18 — Count fitted parameters with _independent_param_count where only the count is needed (Integration buckets and the bi-bi init tally)
- Lever: cheaper_same_coverage; time saved ≈ 5 s; test lines removed ≈ 45; risk medium
- Sign-off: approved by Denis
- Measure first: Fresh process per run: time `[length(fitted_params(compile_mechanism(m))) for m in ms]` against `[_independent_param_count(m) for m in ms]` over fresh types (the uni-uni+R closure at cap 8). Record a Set of the types compiled only for counting in one focused run.

**Testsets.** test/test_mechanism_enumeration.jl Integration 7197-7305 (19.5 s) and the count sites listed in #161

**Change.** - In enumerate_all_mechanism (test_mechanism_enumeration.jl:252-253) use `actual(m) = EnzymeRates._independent_param_count(m)`, memoized per call because add! and the sweep both call it.
- Do the same in the bi-bi init tally (7241-7243) and in `deltas` comprehensions that have no rank-oracle neighbour.
- Keep the per-member `length(fitted_params(compile_mechanism(m))) == pc` loops at 7205-7208 and 7273-7276. They become real cross-checks instead of tautologies.

**Coverage.** - fitted_params' generator is _dependent_param_exprs(Mechanism(M())), the same computation as _independent_param_count(m) after the Sig round trip. That round trip is pinned on all 43 specs (7395-7403) and in test_types.
- The kept loops turn tautologies into real cross-checks.
- Several sites call the compiled count 'ground truth' (2323, 3221, 3281, 3341, 4444, 4693). Denis has a history of 'three parameter counts disagree', so this needs his explicit OK. If declined, the sites stay as they are.

_Merged from:_ tests_all#161 (suite_level); tests_all#31 (enum_B)

#### time/TT20 — Dedup-key partition: render each mechanism once; re-render a fixed sample for determinism
- Lever: lower_threshold; time saved ≈ 5 s; test lines removed ≈ 8; risk low
- Sign-off: approved by Denis
- Measure first: Warm session: for each bi-bi init mechanism, time the first and second render. The saving is about Σ of the second renders minus the 10 kept.

**Testsets.** test/test_identify_rate_equation.jl:707-764 'rate-eq dedup-key partition stability' (55.9 s)

**Change.** In test_identify_rate_equation.jl:752-759 render once to bucket. Re-render and assert `=== h` only for `i % 25 == 1` (10 of 240). Drop the vacuous uni_uni entry (1 mechanism, 1 class), so the test is a single bi-bi block. Keep the 239-class pin.

**Coverage.** - The 239-class pin over all bi-bi init mechanisms stays. That is the regression value.
- src holds no global cache, RNG or objectid hashing, so a second in-process render can differ only through global state, which would not depend on the mechanism.
- Remaining confidence: a defect hitting a fraction f of mechanisms is caught with probability 1-(1-f)^10, against 1-(1-f)^239 today.
- The second render is warm (~25 ms), so only that part is saved. The first-render specialization is TT35's target.

_Merged from:_ tests_all#139 (identify_fit_tests); tests_all#165 (suite_level)

#### time/TT21 — compile_budget: run the cold uni-uni init and the first-call rate_equation timing in one subprocess (only if they do not warm each other)
- Lever: cheaper_same_coverage; time saved ≈ 5 s; test lines removed ≈ 15; risk medium
- Sign-off: not required (coverage unchanged)
- Measure first: Three runs of each order. Accept only if UNI_COLD stays in its separate-process band (1-2 s) and the first rate_equation call stays near its separate-process baseline (~1 s).

**Testsets.** test/test_compile_budget.jl:118-153 and the cold_script of 'compile reuse'

**Change.** After TT03, merge the wall-clock script (118-153) and the cold uni-uni script into one subprocess that prints UNI_COLD and ELAPSED, in whichever order the measurement accepts.

**Coverage.** The same two quantities and thresholds are measured. Contested: the #166 verifier argues the two cold timings must stay separate, because shared Mechanism/Step JIT would weaken the 6 s first-call gate or distort the reuse denominator. Adopt only if the measurement shows neither run warms the other.

_Merged from:_ tests_all#89 (derivation_tests)

#### time/TT22 — Derivation spec oracles: fewer random draws (analytical 20→3, string 10→2)
- Lever: lower_threshold; time saved ≈ 4 s; test lines removed ≈ 0; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_rate_eq_derivation.jl:868-893, 960-986

**Change.** In test_rate_eq_derivation.jl, set test_analytical_rate n_trials=3 (868) and test_rate_equation_string `all(1:2)` (976). Leave the 40-seed Fix B and the reproducer loops alone.

**Coverage.** - Every analytical oracle is a branch-free rational function, apart from one single-draw max().
- Two different rational functions disagree almost everywhere, so 3 or 2 generic points detect any wrong term, sign, name or dependent value at rtol 1e-10.
- Remaining exposure: a draw landing on a near-root, which is negligible.
- The time is incremental after TT07 and TT16; the string-eval share does not shrink with caching.

_Merged from:_ tests_all#78 (derivation_tests)

#### time/TT23 — Enumeration file: small subsumed checks with measurable time
- Lever: delete_redundant; time saved ≈ 3.3 s; test lines removed ≈ 144; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl 440-709, 1368-1472, Dedup (3.4 s)

**Change.** - #168/#4: delete 'Init compiles for all small reactions' (1446-1458), 'bi-bi exit gate' (1460-1472) and the repeated connectivity asserts (1390-1391, 1409-1410, 1441-1442).
- #5: compute ter-ter topologies once inside the _catalytic_topologies testset, and delete 'weak-ordering combining' (522-538).
- #6: delete the duplicate connectivity testset in init_mechanisms.
- #25: delete the duplicate full-rank check on parent 3952.
- #27: one init_mechanisms call per reaction in that testset.
- #33: delete the five Dedup testsets subsumed by the allunique, equal-hash and canonical tests.

**Coverage.** - Each deleted assertion is a strict subset of a full-set assertion on an identical reaction: connectivity at 1368-1379, the 239 pin at 1504, compile on every bi_bi_pp init (12172), and rate_equation_string on all 239 bi-bi inits at identify:751-757.
- #168: commit message must name identify:751-757 as a subsumer, because trimming that test later would drop this coverage silently.
- The 5 bi-bi renders' JIT moves to the identify partition, so only about 0.5 s of #168 is saved.

_Merged from:_ tests_all#168; tests_all#4; tests_all#5; tests_all#6; tests_all#25; tests_all#27; tests_all#33; tests_all#23

#### time/TT24 — Build MECHANISM_TEST_SPECS from top-level statements instead of one 2,460-line function
- Lever: cheaper_same_coverage; time saved ≈ 3 s; test lines removed ≈ 0; risk low
- Sign-off: not required (coverage unchanged)
- Measure first: Fresh process: after `using EnzymeRates` and one @enzyme_mechanism plus one @allosteric_mechanism warm-up, include a scratch copy without the final const line. Compare `@time precompile(build_mechanism_test_specs, ())` with `@time build_mechanism_test_specs()`. Adopt only if the precompile share is several seconds.

**Testsets.** test/mechanism_definitions_for_test_enzyme_derivation.jl:103-2568 (DEFS_LOADED 18.6 s)

**Change.** `const MECHANISM_TEST_SPECS = MechanismTestSpec[]`, then each existing let block at top level with push!. Hoist any helper closure that several lets share.

**Coverage.** Same specs, values and order; only how the fixture table is evaluated changes. Most of the 18.6 s is first JIT of the DSL and constructors, which would move to the next caller, so only the giant-function inference and codegen is saved.

_Merged from:_ tests_all#77 (derivation_tests)

#### time/TT25 — Derivation file: delete duplicate specs and same-branch error checks
- Lever: delete_redundant; time saved ≈ 2.6 s; test lines removed ≈ 227; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_rate_eq_derivation.jl 1957-2018, 1830-1879, group-rep and §5a testsets; definitions file; reference/allosteric_golden_reference.txt

**Change.** - #81/#96: delete spec 'MWC Dimer + Independent Inhibitor', an exact duplicate of 'Homodimer + Non-competitive Inhibitor', and its 4 golden lines.
- #90: delete m_cyclic from 'Rate equation too large error'.
- #91: merge the two 'all-RE catalytic cycle raises' testsets.
- #82: fold the group-rep invariance tests into the reversing test.
- #85: consolidate the dependent-parameter soundness checks.

**Coverage.** - The spec is character-identical to another one.
- Both 'too large' mechanisms are all-SS and reach the same _assert_derivable → _segment_graph_terms branch, which has no shape-dependent path. The surviving @test_throws keeps the message check.
- The all-RE testsets raise on the same line.
- The group-rep and soundness checks are restated by the reversing test and test_constraint_counting.
- #90 and #91 are unverified (spotted by verifiers), so read them once before deleting.

_Merged from:_ tests_all#81; tests_all#96; tests_all#90; tests_all#91; tests_all#82; tests_all#85

#### time/TT26 — Spec metabolite order: canonical HK/PFK-1/PK name lists, so rate_equation compiles once per spec
- Lever: cheaper_same_coverage; time saved ≈ 2.5 s; test lines removed ≈ 43; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/mechanism_definitions_for_test_enzyme_derivation.jl; test/test_rate_eq_derivation.jl test_structure and the six test_* users

**Change.** Reorder the metabolite_names lists of HK (defs :2204), PFK-1 (:2071) and PK (:2305) to metabolites(m) order. In test_structure assert `metabolites(m) == Tuple(spec.metabolite_names)`. Delete the now-implied expected_n_metabolites field and its 43 lines.

**Coverage.** Today two concs NamedTuple orders reach the @generated allosteric rate_equation for these three specs, so each compiles twice. The new assertion pins names and order, which is stronger than the count it replaces. By-name concs access stays covered by 'structural names + synth-dep routing', which passes a PEP-first order.

_Merged from:_ tests_all#79 (derivation_tests, verifier variant)

#### time/TT27 — Delete 'a base-tier row of a degenerate seed's child carries no parent'
- Lever: delete_redundant; time saved ≈ 2.5 s; test lines removed ≈ 30; risk medium
- Sign-off: approved by Denis

**Testsets.** test/test_identify_rate_equation.jl:1746-1775 (2.5 s)

**Change.** Delete test_identify_rate_equation.jl:1746-1775.

**Coverage.** - `child in base && !(child in seeds)` follows from 1677-1679 (|base| = 278 = 257 + 21, after unique!).
- The test bypasses _beam_search. Parentless base rows are asserted on the real base tier at 365-367.
- mechanism_type is covered at 310 and n_params at 935.
- Only an incidental ping-pong child compile is lost, and no compile moves.

_Merged from:_ tests_all#142 (identify_fit_tests); src_enum#60(f)

#### time/TT28 — Derivation-guard testset: a smaller over-limit mechanism
- Lever: cheaper_same_coverage; time saved ≈ 2 s; test lines removed ≈ 0; risk low
- Sign-off: not required (coverage unchanged)
- Measure first: Warm @elapsed of the testset with each fixture.

**Testsets.** test/test_identify_rate_equation.jl:1050-1075 (3.4 s)

**Change.** Replace the 25-step random ter-ter (test_identify_rate_equation.jl:986-1002) with the 16-step bi-bi+R1 m_manual from test_rate_eq_derivation.jl:1962-1982, written inline with an R1 data column. Keep all four assertions, including `_eq_complexity(m) > MAX_RATE_EQUATION_TERMS` and the 'polynomial terms' message.

**Coverage.** The same render-stage throw inside _compile_batch's try, with the same guard message; the guard is a pure V×τ threshold.

_Merged from:_ tests_all#153 (identify_fit_tests)

#### time/TT29 — Depth-2 aggregate: read the no-inhibitor reaction's counts off rxn6's inhibitor-free mechanisms
- Lever: cheaper_same_coverage; time saved ≈ 2 s; test lines removed ≈ 8; risk medium
- Sign-off: approved by Denis
- Measure first: Measure the rxn share together with TT08's timing.

**Testsets.** test/test_mechanism_enumeration.jl:11186-11214

**Change.** After TT08:
```
free(m) = isempty(_bound_comp_inhibitors(m))
@test [count(free, l) for l in iso6] == [62, 349, 1134]
@test [count(free, l) for l in rest6] == [202, 669, 1237]
```
Delete the rxn reaction and its growth.

**Coverage.** - Moves never remove steps, so an inhibitor-free mechanism is reached only through inhibitor-free ones, at the same BFS distance.
- _redundant_copy_groups returns early when no inhibitor is bound.
- Lost: a direct check of the moves on a reaction that declares no inhibitors. That path stays covered by the bi_bi_rxn aggregate pins (9835, 10914) and the hyperbolic sweep.
- Unverified (spotted by a verifier).

_Merged from:_ tests_all#72 (enum_C)

#### time/TT30 — Derivation file small cheaper checks: pin the large-equation mechanism inline; derive the spec string once; evaluate with fitted parameters only
- Lever: cheaper_same_coverage; time saved ≈ 1.8 s; test lines removed ≈ 9; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_rate_eq_derivation.jl:1785-1828 and the per-spec string and factored tests

**Change.** - #80: write the mechanism 'Large equation compilation' selects inline, keep the t_compile < 20 gate, and fix the false comment at 1786-1788.
- #84: derive rate_equation_string once per spec and fold test_factored_form into test_rate_equation_string.
- #95: evaluate the whole string, constraint lines plus v, with only the fitted parameters.

**Coverage.** The timed compile of the same large equation is kept, and dead-end enumeration stays covered in the enumeration file. The string checks see the same string, now derived once.

_Merged from:_ tests_all#80; tests_all#84; tests_all#95

#### time/TT31 — Integration: enumerate the uni-uni closure once and drop the tautologies
- Lever: test_code_simplification; time saved ≈ 1.5 s; test lines removed ≈ 15; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:7197-7305

**Change.** Compute `enumerate_all_mechanism(uni_uni_rxn; max_params=8)` once in the outer Integration testset and reuse it at 7200, 7291 and 7301. Merge 'Multiple levels populated' into it. Drop `issorted(sort(...))`. Drop the recount loops only if TT18 is declined; under TT18 they become real cross-checks.

**Coverage.** Identical results are reused. issorted(sort(x)) cannot fail, and the recount loops recompute the bucket key with the same expression.

_Merged from:_ tests_all#175 (suite_level); tests_all#34 (enum_B)

#### time/TT32 — Dead-end section: delete rank assertions that _redundant_copy_groups makes on the same steps
- Lever: delete_redundant; time saved ≈ 0.4 s; test lines removed ≈ 17; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:3629-3795, 4600-4617

**Change.** Delete test_mechanism_enumeration.jl:3705, 3794, 3795, 4615 and 4617, and the now-dead relaxed_s fixture (4600-4611).

**Coverage.** - Step-identical twins carry the same (fitted, rank) pins: b1 at 8632, b14 at 8673, case3 at 8653, and relaxed_s at 8956-8957 (stronger, `== [copy_group]`).
- The local package decisions (`!(x in kids)`, exact Set(kids)) stay.
- Saves 2.5 s if TT06 is rejected (4 withrxn compiles), and about 0.4 s under TT06.

_Merged from:_ tests_all#2 (enum_A)

#### time/TT33 — Rank oracle on absent children: one representative per symmetry orbit
- Lever: delete_redundant; time saved ≈ 0.4 s; test lines removed ≈ 45; risk medium
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl 9844-9883, 10277-10282, 10285-10360, 10584-10641, 10921-10953, 11373-11377

**Change.** Keep one fixture per orbit:
- ping-pong flipA (drop flipB, flipP, flipQ);
- one uni-uni flip;
- flip(A1,B2) only for the ter-ter rank line;
- the A group's two bipartitions in the rejected-split loop;
- plus the zero-flux-pair and ordered-TC mirrors the verifier found.

**Coverage.** - Each dropped child is the image of a kept one under a relabeling or reversal symmetry of the parent, and identifiable rank is invariant under those.
- Lost: derivation in the dropped orientations. Relabeling changes canonical step order and so the Haldane pivot, so a pivot bug in one orientation would go unseen here.
- Saves about 3.3 s if TT06 is rejected, about 0.4 s under TT06. Recommend skipping it if TT06 lands.

_Merged from:_ tests_all#58 (enum_C)

#### time/TT34 — _seed_variants tests: numeric degeneracy probe on one candidate per failure kind per base
- Lever: lower_threshold; time saved ≈ 0.3 s; test lines removed ≈ 0; risk medium
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:7809-7817, 8021-8029, 8119-8126, 8201-8204

**Change.** Keep `ER._degenerate` on every rejected candidate. Run _testhelper_degenerate only on a representative set:
- shared-B: 6 of 13;
- ping-pong: 10 of 20, keeping [binds(x), binds(y)] per TC base as the verifier corrected;
- ordered/random: 4 of 8;
- C on ping-pong: :A and :B.

The full 14-candidate probe table at 7716-7725 stays.

**Coverage.** - The structural predicate the enumerator uses is still asserted on every candidate.
- The numeric oracle still confirms every failure kind (V substrates, V products, C) on every base, but no longer on each individual subset.
- Saves about 3 s if TT06 is rejected, about 0.3 s under TT06. Recommend skipping it if TT06 lands.

_Merged from:_ tests_all#32 (enum_B)

#### time/TT35 — src: stop per-Sig JIT in the Sig decode and the Type{M} derivation entry points (rate_equation body untouched)
- Lever: src_change; time saved ≈ 30 s; test lines removed ≈ 0; risk medium
- Sign-off: approved by Denis
- Measure first: Warm session, 30 fresh bi-bi init types: record @timed compile_time and time of compile_mechanism plus fitted_params, of the first and second rate_equation_string, and of the first rate_equation call. Apply the @nospecialize patch and repeat. Then run test_rate_eq_derivation.jl's performance gate and the partition testset.

**Testsets.** Suite-wide: the identify partition (55.9 s, about 200 ms of each mechanism's ~233 ms is first-render specialization), the Integration closures, the rank oracle's allosteric calls, the identify compile batches and the derivation spec loop

**Change.** Mark the Sig decoders `@nospecialize`: _mechanism_from_sig, _reaction_from_sig, _steps_from_sig, _step_from_sig, _species_from_sig and _metabolite_from_sig in src/types.jl:1060-1127. Also mark the non-generated entry points that take `M::Type{<:EnzymeMechanism}`: _dependent_param_exprs(M), _dependent_param_exprs_kernel(M, …), _raw_symbolic_rate_polys(M) and _rate_v_line(M). Route rate_equation_string(::M, mode) through an unspecialized implementation on Mechanism(em).

Today every new mechanism shape compiles fresh specializations of these functions, because the nested Sig tuple type differs per shape. That cost lands on compile_mechanism, every generator call and every runtime rate_equation_string. Output is unchanged.

**Coverage.** - No test changes and identical outputs. The byte-identical golden file, the factored-form checks, the 239-class partition pin and the fitted_params pins guard equality.
- rate_equation's generated body and its 0-alloc/<120 ns contract are untouched, because the generators run these functions only at generation time. Re-run test_rate_equation_performance to confirm.
- The time is incremental after TT02, TT09 and TT20, which already remove some per-type work. Range 10-60 s. It also speeds production identify, which renders the string per candidate.
- Needs Denis's approval as a src change.

_Merged from:_ tests_all#180 (suite_level); src_derivation#107; src_derivation#33; src_types#4

#### time/TT36 — src: profile-guided enumeration speedups (expand_mechanisms / init_mechanisms), output-identical
- Lever: src_change; time saved ≈ 25 s; test lines removed ≈ 0; risk medium
- Sign-off: approved by Denis
- Measure first: Profile.@profile (one Julia process, nothing else running) on the three workloads above, and act on the top frames. Re-time 'init_mechanisms on ter-ter within 150 s' and the hyperbolic sweep before and after each change.

**Testsets.** Enumeration-dominated rows after the plan: the hyperbolic sweep (~145 s left), depth-2 (~57 s), in-process ter-ter init (45 s), Integration, the BFS and the split/flip testsets (~30 s)

**Change.** Profile first: Profile.@profile on init_mechanisms(ter_ter_rxn), on _testhelper_levels(pingpong, 2) and on the rxn6 depth-2 BFS. Candidates already found in code reading:
- _expand_to_allosteric runs _onlya_haldane_violation explicitly (mechanism_enumeration.jl:2679), then again inside the AllostericMechanism constructor (types.jl:922), for every mask.
- _canonical_group_order! re-renders String sort keys on every comparison (types.jl:615-640). Precompute them; src_types#9 shows the order is identical.
- _split_gain_test runs a full _independent_param_count(am) kernel solve per candidate set. The Pass-1 kernel is also repeated (src_derivation#33).
- _rep_step recomputes _free_enz_set on every name(p, m) call (src_derivation#8).
- _state_view memoization (src_derivation#81).

**Coverage.** - Output-identical refactors gated by every enumeration count pin and canonicalization test.
- src_types#9 touches Canonical Step Form, and is labelled high risk only for that reason; it is behaviour-identical (stable sort, same keys).
- 25 s is a placeholder for an unprofiled 10-20% speedup of roughly 280 s of enumeration. A 1.5-2× speedup would be worth 90-140 s.

_Merged from:_ src_enum#33; src_types#9; src_derivation#8; src_derivation#81; src_derivation#32; src_enum#57; src_enum#97

### TEST CODE SIZE of test/test_mechanism_enumeration.jl (12,734 lines). Merged items from scopes enum_A, enum_B, enum_C, enum_hyper and suite_level, plus src items that delete tests in this file. Every item respects Denis's enumeration-test rules (fixtures stay inline in their testset, exact child sets stay).

I read the code for every item that overlaps another or that the verifiers disputed, and I did not run Julia. The ~70 proposals that touch this file merge into 26 items.

**Line counts**
- Test-only items T1-T24 remove about 1,584 lines, 12.4% of the file.
- Two conditional items remove about 353 more. They apply only if Denis accepts the matching src deletions (T25: test-only src helpers; T26: the dead-end opportunities inline and the _seed_candidate_screen deletion).
- The largest blocks:
  - Property-only move tests superseded by later exact-children tests on the same seed: T6, T8, T12, T18, T24.
  - Copy-pasted testset bodies turned into tables, with fixtures still inline: T7, T9, T13.
  - Repeated test-local closures turned into _testhelper_ helpers: T15.
  - Tautological asserts that the constructor or element type already guarantees: T14.

**Conflicts I resolved by reading the code**
- **relaxed_s duplicate.** Tests#2 and tests#35 each delete the other's copy of the duplicate relaxed_s block. I kept #35 (T18): it saves 55 lines against 17, and T18 moves the only exact (fitted, rank) == (6, 5) pin to line 4616 so that pin survives.
- **Init connectivity check.** Tests#6 and tests#168 each delete a different copy of the same check. T5 keeps one loop.
- **Ter-ter budget tests.** The time version of tests#56 depends on line 7151, which tests#40 (in T13) deletes. So T21 is the merge with no coverage loss: 33 lines, 0 s. The time-saving variant is left to the time concern.
- **exclude_regs.** Tests#13 asserts the I-only count through exclude_regs, which src_enum#49/#66 deletes. T10 and T25 say how to re-home that assert.
- **Duplicate subsumers.** Tests#23 is a subset of tests#33, and tests#16 duplicates tests#67.

**Additions of mine, not checked by any verifier**
- The `_inhibitor_competition_patterns` table in T2 (about 15 lines).
- Two extra deletions in T24 that 284-325 subsumes (about 15 lines).

**Rule note for Denis**
Many move tests here still assert only properties, which breaks CLAUDE.md rule 2 (exact children). T7, T9 and T13 keep those gaps exactly as they are. Writing out the expected children would add roughly 60-200 lines per area, so Denis needs to decide that separately.

**Time**
The items total about 105 s. About 96 s of that is the same saving the time concern counts (T19 = tests#47's 90 s, T11, T20). Pure time items with no line change are listed under dropped.

#### loc/T1 — _catalytic_topologies: one count table, is_iso instead of the adapters, one iso predicate, delete the duplicate weak-ordering testset
- Lever: test_code_simplification; time saved ≈ 0.7 s; test lines removed ≈ 110; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:8-30, 416-713

**Change.** Replace the six count testsets (418-488) with one loop over [(uni_uni,1),(uni_bi,3),(bi_bi,9),(bi_bi_pp,10),(ter_ter,223),(ter_bi,45)]. Each entry asserts the count, all-connected and one SS step. Keep the uni-bi EnzymeMechanism line and the per-reaction count comments.

Also:
- Delete 'weak-ordering combining' (522-538).
- Compute ter-ter topologies once inside the _catalytic_topologies testset; today they are computed at 466, at 502 (an identical ter_ter re-declared at 498-501), at 534 and at 541. Drop the repeated connectivity asserts at 503 and 543.
- Replace _t_reactants/_t_products (22-30) with EnzymeRates.is_iso at 561-562 and 704-705. Rewrite the adapter comment at 8-13; _iso_orient stays.
- Turn each symmetric 8-line `any(spec) do` block in the pyruvate tests into one line through a local iso_between closure.
- Fold the per-element @test loops (C5, C8, conformation) into @test all(...).

**Coverage.** Every predicate is unchanged. Topology steps are only plain bindings and isomerizations, so the adapter test is exactly is_iso. iso_between is the same orientation-agnostic check. The deleted 522-538 re-pins 9 and 223 on reactions identical to those at 441-452 and 465-474. _catalytic_topologies is deterministic (it sorts the iso groups), so a shared vector gives the same assertions. Folding the loops cuts the @test record count from about 15.5k to a few hundred: reporting granularity only.

_Merged from:_ tests#19; tests#5; src_enum#15

#### loc/T2 — _competition_patterns and _inhibitor_competition_patterns as tables with one cover predicate
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 45; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_mechanism_enumeration.jl:714-838

**Change.** Define `covers(pat,S,P)` once. Loop over [(2x2,7),(3x3,265),(2x3,25),(3x2,25)], asserting the count and all(covers). This replaces the four copy-pasted cover loops (737-747, 756-767, 774-781, 788-795). Keep the 1x1/1x2/2x1 exact-pattern asserts and the invalid-2x2 check.

Also drive 800-837 (_inhibitor_competition_patterns counts 9/49/18/36/72) from one (S, P, existing, n) table. Keep the uni-uni exact pattern and the 2^k comments.

**Coverage.** The same counts and the same cover property run on every pattern, and the exact-pattern and invalid-pattern asserts stay. The test count drops from about 1,890 to about 20 records; no predicate is lost. The inhibitor-pattern table (about 15 lines) is my own addition, not reviewed by a verifier.

_Merged from:_ tests#22

#### loc/T3 — Dead-end opportunity and shape tests: exact sets instead of a filter written in the test and repeated calls
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 65; risk low
- Sign-off: approved by Denis
- Measure first: One focused run to confirm the expected form names and composition sets.

**Testsets.** test/test_mechanism_enumeration.jl:929-1018, 1046-1086, 1169-1225

**Change.** In 'Bi-Bi random: 4 dead-end forms' (1046), assert the exact set of new-form sets. That set is the complements of the 7 edge covers of K(2,2): {AQ,BP}, {AP,BQ}, {BQ}, {BP}, {AQ}, {AP}, {}. Keep `length(result) == 7` (1085). Delete the 'complete competition -> bare topology' and 'diagonal has exactly 2 dead-end forms' subtests (1195-1225), which rerun the same call; keep the ter-ter per-topology subtest.

In _substrate_product_dead_end_opportunities (929-1018):
- Drop the diagonal filter written in the test (983-997) and the vacuous '∉ allowed' checks (1015-1017).
- Assert the exact composition of all 27 forms (1S+1P, 2S+1P, 1S+2P over [:A,:B,:D] × [:P,:Q,:R]).
- Keep a names assertion: either Set(de_form_names) == the 27 names, or the 12 lines rewritten as `in de_form_names`.

**Coverage.** Stronger than today. The composition check fixes all 27 forms, where today only the 12 that pass the test's own filter are fixed. Equality of the exact set plus the count of 7 implies one bare variant and two 2-form variants, which is what the two deleted subtests assert on the same seed and reaction. The names assertion keeps the rendered-name coverage. If src_enum#10 lands, 929-1018 is replaced anyway (T26) and the 25 lines of this half are moot.

_Merged from:_ tests#20; tests#21

#### loc/T4 — Delete the fixture sanity testsets for atom balance and uni-uni invariants
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 23; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:327-334, 1296-1310

**Change.** Delete 'test reaction atom balance' (1296-1310) and '_assert_mechanism_invariants on uni-uni init' (327-334).

**Coverage.** The EnzymeReaction constructor errors on any per-element atom imbalance (src/types.jl:449-466). The two pyruvate consts are built at file load, and test_types.jl:1011 tests the rejection. For 327: 345-347 runs _assert_mechanism_invariants on the same uni-uni init steps, and the unbound inhibitor declared there is skipped by the check.

_Merged from:_ tests#15

#### loc/T5 — init_mechanisms testset: one connectivity loop, call init once per reaction, delete the subset compile and derive checks
- Lever: delete_redundant; time saved ≈ 1.5 s; test lines removed ≈ 48; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:1366-1507

**Change.** Inside 'init_mechanisms':
- Compute bb = init_mechanisms(bi_bi_rxn) and pp = init_mechanisms(bi_bi_pp_rxn) once.
- Keep the connectivity assert only in the 1386 loop over the four reactions. Delete 'init mechanisms are connected' (1368-1379) and the repeats at 1409-1410 and 1441-1442.
- Fold 'Uni-uni: exactly 1 init mechanism' (1436-1444) into that loop as `rxn === uni_uni_rxn && @test length(specs) == 1`, keeping its comment.
- Delete 'Init compiles for all small reactions' (1446-1458) and 'bi-bi exit gate: init mechanisms derive (subset)' (1460-1473).

**Coverage.** Each deleted assert is a strict subset (first 5 or smallest 5) of a full-set check on the identical reaction:
- uni-uni compile: 1347-1348.
- All 239 bi_bi init mechanisms through compile_mechanism and rate_equation_string: test_identify_rate_equation.jl:751-757 (239 distinct keys at 759). This subsumer is in another file; name it in the commit.
- Every bi_bi_pp init through compile_mechanism: 12172-12175. Any rewrite of init division-freeness (tests#54/#157) must keep `compile_mechanism(m) isa EnzymeMechanism` per mechanism.
- Count 239 and allunique: 1504.
- Connectivity: the kept loop over the same four reactions.

_Merged from:_ tests#4; tests#168; tests#6; tests#27

#### loc/T6 — Delete the old _expand_re_to_ss property tests that exact-children tests on the same seed subsume
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 124; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:2037-2151, 2732-2759; subsumers 1741-1813, 9844-9883, 9885-9947

**Change.** - Delete 'bi-bi sequential: 4 RE binding groups -> 2 variants' (2037-2104).
- Delete 'bi-bi multi-step kinetic group: atomic conversion' (2106-2151).
- Delete 'uni-uni: both RE binding groups are chain flanks' (2732-2744), after moving its `_chain_flank_groups(m) == Set(non-iso groups)` assert into 9844.
- Write 2746-2759's seed inline instead of pulling it from init_mechanisms (rule 1). This changes no line count.

**Coverage.** - 2037: its seed is the 1747-1756 seed with steps in a different order. Canonical construction makes them equal (tested at 284-325). 1778-1779 asserts length 2 and Set == {a_flip, q_flip}, which implies every property 2037 checks.
- 2106: its seed is identical to 9893-9897. 9939-9940 asserts exactly {flipA, flipB, flipP, flipQ}, each flipping both steps of one shared group.
- 2732: 9877 asserts isempty(kids) on the same three steps, and the flank-set assert moves over. Neither _chain_flank_groups nor the move reads atoms.

_Merged from:_ tests#8; tests#9; tests#16; tests#67

#### loc/T7 — Four uni-uni allosteric RE→SS testsets as one table, each seed still written inline
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 120; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_mechanism_enumeration.jl:2299-2475, 2761-2808

**Change.** Merge 2299-2353, 2355-2417, 2419-2475 and 2761-2808 into one testset that loops over (@allosteric_mechanism seed written inline, deltas). The deltas are [1,1], [1,1], [1,2] and [2,2]. Keep each seed's explanatory comment, including 2378-2385.

The shared body asserts:
- length 2;
- the sorted fitted deltas;
- invariants on every child;
- exactly one newly all-SS group per child;
- cat_allo_states, multiplicity, regulatory_sites and reaction preserved.

**Coverage.** I read all four bodies. They differ only in the seed tags and the deltas. Every assertion of every body runs on every seed; 2761 gains the invariant check, which is strictly more. fitted() goes through compile_mechanism, so the compile check is implied. Rule 1 holds because each macro stays inline in the testset. The tests remain property-only (a rule 2 gap that predates this item). Writing out the 8 expected children would add about 60 lines, so Denis needs to decide that.

_Merged from:_ tests#17

#### loc/T8 — Delete the two old split tests that duplicate a parent tested elsewhere
- Lever: delete_redundant; time saved ≈ 0.25 s; test lines removed ≈ 49; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:2925-2951, 10984-11005; subsumers 10706-10754, 2858-2923

**Change.** - Delete 'Mechanism — bi-bi random: single splits tied, a pair frees a constant' (2925-2951).
- Delete '_expand_split_kinetic_group: allosteric parent' (10984-11005).

**Coverage.** - 2925: its seed is identical to 10710-10719. 10747-10753 asserts exactly {splitAB, splitPQ}, base count 5, each child 6 with 7 groups. That implies non-empty, every child gaining and ≥ +2 groups.
- 10984: its parent is line for line the em_seed at 2862-2873. Its four asserts match 2881 (non-empty), 2882-2883 (gain), 2897 (isa) and 2908 (tag length, which the constructor enforces anyway).
- 2858 itself stays property-only, a rule 2 gap that predates this item.

_Merged from:_ tests#10; tests#64; tests#11; tests#61

#### loc/T9 — Three bi-bi + I dead-end testsets as one table; fold the I-binding-steps-share-one-group testset into it
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 133; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:3193-3369, 4752-4795

**Change.** Replace 3193-3251, 3253-3310 and 3312-3369 with one testset. It loops over (seed written inline, reaction written inline, count): sequential 4, random 9, ping-pong 3.

The shared body asserts:
- the count;
- Δfitted == 1 against the compiled count;
- the children compile;
- I steps are present and all in one kinetic group.

Delete 'I-binding steps share one kinetic group' (4752-4795). Add its one existence claim to the random entry: some child has ≥2 I-binding steps.

**Coverage.** The three bodies differ only in seed, reaction and count, and each fixture stays inline in the table. 4752 uses the same seed and reaction as 3256-3273. The random entry asserts the one-group property on every variant, which is a superset of 4752's multi-binding subset, and the added line keeps its non-emptiness claim. The tests remain property-only, a rule 2 gap that predates this item; exact children would add about 190 lines.

_Merged from:_ tests#18; tests#12

#### loc/T10 — Delete 'uni-uni + I: 1 variant' and 'dead-end on allosteric base'
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 57; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:4797-4834, 4836-4864, 4875-4894

**Change.** - Delete 4797-4834. It also breaks rule 1 by pulling its fixture from init_mechanisms.
- In 4836, add `@test length(EnzymeRates._expand_add_dead_end_regulator(m, rxn_ij; exclude_regs=Set([:J]))) == 1`.
- Delete 4875-4894, which is weak, property-only and also pulls its base from init_mechanisms.

**Coverage.** - 4797: the I-only count becomes a direct assert through the added line. The compile of uni-uni+I children is covered at 1592-1601, which also asserts :I in some child's regulators.
- 4875: the closest subsumers are 4445-4498 and 4500-4541, exact children on :OnlyA bases at multiplicity 2 with multiplicity carried. The foreign-inhibitor placement on a K-type base is covered only through the shared placement path (the copy filter is skipped), so Denis must accept that.
- If src_enum#49/#66 deletes exclude_regs (T25), 4836 goes with it. The I-only count must then stay as a 6-line inline uni-uni+I count check, and this item saves about 32 lines instead.

_Merged from:_ tests#13; tests#14

#### loc/T11 — Delete the rank assertions in the dead-end section that are repeated on identical step sets
- Lever: delete_redundant; time saved ≈ 3 s; test lines removed ≈ 5; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:3705, 3794-3795, 4097

**Change.** Delete these lines:
- 3705: `rank(at_E) == rank(m)`. Keep its comment beside `!(at_E in kids)`.
- 3794-3795: the at_E_EQ and at_E rank pins. Keep r0 == fitted(m).
- 4097: the rank of m against its fitted count.

**Coverage.** Each deleted claim is asserted on the same steps elsewhere:
- at_E (3676) has b1's steps: rank 5 at 8632. m (3629) is the 1747 seed: rank 5 at 1783.
- at_E_EQ (3775) has b14's steps: fitted 6, rank 5 at 8673.
- at_E (3763) has case3's steps: fitted 6, rank 6 at 8653.
- 4097's m has at_E_EQ's six steps, whose fitted == rank is asserted at 3702.
Rank does not depend on atoms or on unbound declared regulators. The package decisions `!(x in kids)` stay local.

The 4615/4617 half of tests#2 is not taken; see T18.

_Merged from:_ tests#2; tests#25

#### loc/T12 — _expand_to_allosteric: delete four subsumed testsets; replace the two uni-uni tests with one exact-children test
- Lever: delete_redundant; time saved ≈ 0.1 s; test lines removed ≈ 138; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:4904-5211, 5045-5101, 12486-12539

**Change.** Delete these four testsets:
- 4904-4931 (oligomeric_state from reaction);
- 5140-5149 (single-valued guard);
- 5164-5181 (no all-:EqualAI baseline);
- 5183-5211 (catalysis-:OnlyA is V-type only).
Add the mandatory assert to the V-type testset (5213): every V-type has its catalytic group(s) :OnlyA.

Replace 5045-5101 and 12486-12539 with one testset:
- The parent is written inline and reattached: `Mechanism(rxn, steps(Mechanism(em)))`.
- The expected children {S}, {P} and {S,P} are written with @allosteric_mechanism and lifted onto rxn with AllostericMechanism(rxn, steps, tags, 2, RegulatorySite[]).
- Assert `length(kids) == 3 && Set(kids) == Set(expected)`, and keep the invariant and compile smoke checks.

**Coverage.** Multiplicity from the reaction is pinned at 5113-5114 ({2,4}) and by the exact children (2); the oligomeric_state:4 DSL path is tested at test_dsl.jl:466-475. No all-:EqualAI child: 5126-5131 for cn=4 and 5295 for the regulator reaction. 'No bare catalysis-:OnlyA' is subsumed by 5299-5300. The one claim nothing else subsumes, a V-type with :OnlyA catalysis, becomes the added assert.

The exact set equality implies everything the two uni-uni tests check: multiplicity, empty sites, reaction, tags, :OnlyA count, uniqueness and shapes. It fixes the rule 1 violation (5055) and the rule 2 violations in both tests.

_Merged from:_ tests#36; tests#62

#### loc/T13 — _expand_add_allosteric_regulator: merge the identical-seed testsets, a table for R/S/P, and delete two hand-rebuilt duplicates
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 230; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:5430-5905, 7138-7151

**Change.** (1) 5527-5586, 5588-5632 and 5752-5791 have the same seed. The reaction is uni_uni_allo_2reg in all three (5752-5791 re-declares it locally). Merge them into one testset with one call and the union of assertions: count 6, deltas [1,1,1,1,2,2], R2 in every child, an R2 :EqualAI child at site 1, 3 new-site and 3 existing-site children, preservation.

(2) Table-drive 5430-5502, 5690-5750 and 5793-5843 over inline reactions declaring R, S or P. Assert the exact triple set {(x,:OnlyA),(x,:OnlyI),(x,:NonequalAI)}, deltas [1,1,2], invariants and preservation. Move the dual-role comments into the table.

(3) Delete 'uni-uni + R: enumerate variants' (5870-5905).

(4) Delete the hand-rebuilt second half of 'Merge move wired' (7138-7151).

**Coverage.** - (1) and (2): each call is deterministic, so one call carries every assert the copies made separately. The exact triple set is stronger than the 'ligand present' checks.
- (3): the per-child site and multiplicity checks are covered by 5430-5502 (the move never branches on catalytic tags), and K-type non-emptiness by 6031-6041.
- (4): merge being a no-op on m is pinned at 6531-6534; expand == filter(_add_expansions_mech!) at 7097-7103; which moves are wired, as a Set, at 11271-11273. The one thing lost is the pin that RE→SS children come before to_allosteric children in the output; Denis must accept that.

The tests stay property-only (rule 2). Writing out the 6 + 9 expected children would add lines.

_Merged from:_ tests#37; tests#39; tests#40

#### loc/T14 — Delete tautological type, Haldane and tag-length asserts on move results
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 34; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:4982-6295 (isa lines), 9695, 12234-12236, 12550-12634

**Change.** Delete the asserts that are always true:
- `@test r isa EnzymeRates.AllostericMechanism` at 4982, 5038, 5671, 5952, 6100, 6148, 6205, 6263 and 6295. The moves return Vector{AllostericMechanism}.
- `compile_mechanism(r) isa AllostericEnzymeMechanism` at 4984, 5040, 5673, 6150 and 6207, where a delta loop already compiles every r (4975, 5032, 5665, 6139, 6196).
- The haldane_ok helper (12550-12552) and its uses (12572, 12576, 12604, 12626).
- The tag-length and Haldane asserts at 12234-12236, 12632-12634 and 9695.
Fold any loop left holding only `_assert_mechanism_invariants` into its neighbour. Rename '…emits only Haldane-valid mechanisms' to describe what remains, and keep the intent comments.

The other lines tests#38/#65 list sit inside T12, T13 and T8 and are counted there.

**Coverage.** - The isa lines are true by the element type of the moves' return vectors (src/mechanism_enumeration.jl:2661/2706/2754/2892).
- The compile isa lines are implied by fitted_params(compile_mechanism(r)) running on the same r.
- AllostericMechanism has one inner constructor that errors on a tag-length mismatch (types.jl:880-883) and on a Haldane violation (922-924). The struct is immutable and is never mutated, so these asserts cannot fail on any existing object. The guard itself is tested at test_types.jl:2598ff.

_Merged from:_ tests#38; tests#65

#### loc/T15 — File-level _testhelper_ helpers for the repeated test-local closures
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 85; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_mechanism_enumeration.jl file-wide (1782-11661 fitted closures; 3241-4859 I/J loops; 3626-4351 withrxn; 7632-8213 seed closures; 12544)

**Change.** Add file-level helpers that replace the repeated closures:
- `_testhelper_fitted(x)` replaces 23 identical fitted closures.
- `_testhelper_param_deltas(parent, kids)` replaces the remaining delta sites (4973, 5030, 5663, 6137, 6194).
- Use the `fitted_params(::Mechanism/AllostericMechanism)` convenience (src/rate_eq_derivation.jl:99-100) instead of fitted_params(compile_mechanism(x)).
- `_testhelper_holds_iso` (5 copies), `_testhelper_iso_groups` (3 isochem copies), `_testhelper_on_reaction(rxn, em)` (7 withrxn copies with their 2-line comment), a lift helper (3), and `_testhelper_groups_binding(m, x)` for the copy-pasted I/J detection loops.
- In the seed tests: mech/groups (9 copies), binds (5), is_tc (3) and flip(gs, preds) (2). 7791's flip closes over its seed and stays local.
Drop the inner testset that repeats its parent's name at 6340; at 6370 keep the shared setup in scope. Rename the file-level tags_by_bound_metabolite to the _testhelper_ prefix (rule). Fixtures stay inline at their call sites. The body of _testhelper_fitted stays the compiled count.

**Coverage.** This is a pure refactor: the same expressions are called through one definition, and the convenience method forwards to fitted_params(compile_mechanism(m)). The count of 85 nets out the sites already removed by T7, T9, T12 and T13 and adds about 10 lines for the helper definitions. Swapping the helper body to _independent_param_count (tests#161) is a separate time decision and changes the stated 'compiled count is ground truth' intent.

_Merged from:_ tests#63; tests#41; tests#42; tests#18

#### loc/T16 — Integration: enumerate uni-uni once, drop the tautological recounts and the sorted-is-sorted check
- Lever: delete_redundant; time saved ≈ 1.5 s; test lines removed ≈ 22; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:7195-7305

**Change.** - Compute `uni = enumerate_all_mechanism(uni_uni_rxn; max_params=8)` once in the Integration testset and reuse it at 7200, 7291 and 7301.
- Move the ≥2-buckets and gap≤4 asserts from 'Multiple levels populated' (7297-7305) into 7197, then delete 7297-7305.
- Delete `@test issorted(pcs)` (7203), which tests the output of sort.
- Delete the per-member recount loops at 7204-7209 and 7272-7277.

**Coverage.** The same deterministic call backs every assert. The recount loops evaluate the exact expression enumerate_all_mechanism used for the bucket key (252-258, 266) on the same type, so they cannot fail. If the time concern adopts tests#31 (bucket by _independent_param_count), keep the loops, because they then become real cross-checks; the saving drops to about 10 lines.

_Merged from:_ tests#34; tests#175; tests#45

#### loc/T17 — Delete the duplicate bi-bi reaction constant
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 4; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_mechanism_enumeration.jl:7390-7393

**Change.** Delete `_testhelper_bibi_rxn` (7390-7393) and use bi_bi_rxn (167-170) at 9477, 9745, 9835 and 10914. Keep the section comment at 7386-7389.

**Coverage.** The two reactions are byte-identical (A[C], B[N] → P[C], Q[N]), and EnzymeReaction == is structural.

_Merged from:_ tests#177

#### loc/T18 — _productive_twin/_redundant_copy_groups: delete the duplicated allosteric copy fixtures and twin_only; keep the one exact pin
- Lever: delete_redundant; time saved ≈ 0.2 s; test lines removed ≈ 55; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:8464-8478, 8918-8957, 4616

**Change.** - Delete the allosteric kept/relaxed_s block 8918-8957 (rxn, lift, fixtures and the asserts at 8955-8957).
- Delete twin_only and its asserts (8464-8478); keep 8462-8463.
- Rewrite 4616 as `@test (fitted(relaxed_s), _testhelper_identifiable_rank(relaxed_s)) == (6, 5)`, so the only exact pin survives.

This replaces the opposite deletion in tests#2, which removed 4600-4617 instead; the two cannot both be applied.

**Coverage.** - kept (8931) has the same steps and tags as onlya_with_equalai_copy (8516). 8528-8530 asserts it is not redundant and fitted == rank, and 4538 has the identical lifted kept. _redundant_copy_groups reads only steps and tags.
- relaxed_s (8943) is the same lifted type as 4599. 4614 asserts !isempty, which is equivalent to == [copy_group] with one copy group. The (6,5) pin moves to 4616.
- twin_only equals b1 (8619), where 8631 asserts == [copy_group(b1)].

_Merged from:_ tests#35; tests#2

#### loc/T19 — Hyperbolic exactness sweep: drop the bi-bi iteration, which the ping-pong population contains, and fold the local helpers
- Lever: delete_redundant; time saved ≈ 90 s; test lines removed ≈ 28; risk low
- Sign-off: approved by Denis
- Measure first: A one-off check before deleting (not a permanent test): enumerate both populations and assert issubset of the bi-bi keys (steps, tags, multiplicity, ignoring reaction) in the ping-pong keys.

**Testsets.** test/test_mechanism_enumeration.jl:9319-9401

**Change.** - Loop over ((unibi,2), (pingpong,2)) only.
- Delete the bibi reaction (9363-9367), _testhelper_allosteric_children (9349-9356) and the append/unique lines (9379-9380).
- Fold _testhelper_mets, _testhelper_den_hyperbolic and _testhelper_levels into local closures and a 5-line population loop.
- Rewrite the header comment (9320-9327): the ping-pong-capable reaction's sequential seeds are the bi-bi seeds.

**Coverage.** The verifier traced every divergence point:
- The 55 sequential topologies are identical under both reactions (residual-free isos).
- Dead-end decoration, _seed_variants, the moves, _hyperbolic_catalysis and the derivation read names only; atoms appear only in the topology builder and the conservation asserts.
- Both reactions are oligomeric_state 2.
So ping-pong level 1 ⊇ bi-bi level 1, and the allosteric children of level-1 allosteric mechanisms ⊆ ping-pong level 2. The n_checked and n_nonhyperbolic guards stay meaningful.

The 90 s is the same saving the time concern counts under tests#47; do not double count.

_Merged from:_ tests#47; tests#51

#### loc/T20 — Absent-child rank oracles: keep one representative per symmetry orbit
- Lever: delete_redundant; time saved ≈ 3.3 s; test lines removed ≈ 55; risk medium
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:9844-9883, 10213-10360, 10494-10641, 10921-10953

**Change.** - Ping-pong (10285): keep flipA and delete the flipP, flipB and flipQ fixtures (10318-10353).
- Uni-uni (9844): keep one of the two mirrored absent flips (9867-9875).
- Split ter-ter (10278): run the rank line on flip(A1,B2) only, and keep the `!(absent in kids)` and segment-count checks for both.
- Rejected split (10943): loop over the A group's two bipartitions only.
- Zero-flux pair (10584-10641): drop one mirrored absent fixture.

**Coverage.** The verifier checked each orbit: each parent is invariant under its relabeling or reversal, and identifiable rank is invariant under relabeling. The cost is that the dropped orientations are derived nowhere else. Relabeling changes the canonical step order and with it the Haldane/Wegscheider pivot, so a pivot bug confined to one orientation would go uncaught here. 'reversing written steps' does not cover relabeling. Denis must accept that.

_Merged from:_ tests#58

#### loc/T21 — Merge the two ter-ter worst-case budget testsets into one seed and one testset
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 33; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_mechanism_enumeration.jl:11007-11045, 11128-11167

**Change.** Write the 55-step seed once, in one testset:
```
t_split = @elapsed split = _expand_split_kinetic_group(worst)
t = @elapsed kids = expand_mechanisms([worst], terter)
```
Keep `length(split) == 12`, `t_split < 60`, `length(kids) == 81` and `t < 120`. Delete the second 27-line seed copy and its testset.

**Coverage.** Every assertion and both budgets stay, on the identical seed. This is the line-count half of tests#56 with no coverage loss.

The time variant drops the duplicate split run by timing the split and the six other moves separately. It relies on 7151 to pin expand_mechanisms' composition, which T13 deletes. If the time concern takes that variant, T13 part (4) must stay unapplied.

_Merged from:_ tests#56

#### loc/T22 — Depth-2 test: read the no-inhibitor reaction's counts off rxn6's inhibitor-free mechanisms
- Lever: cheaper_same_coverage; time saved ≈ 2 s; test lines removed ≈ 8; risk medium
- Sign-off: approved by Denis
- Measure first: One run to confirm the inhibitor-free counts per level equal [62,349,1134] and [202,669,1237].

**Testsets.** test/test_mechanism_enumeration.jl:11169-11240

**Change.** After the halves restructure (tests#55/#163, time concern), assert the five rxn pins as counts of inhibitor-free mechanisms per level of rxn6's iso and non-iso families:
```
free(m) = isempty(_bound_comp_inhibitors(m))
```
Then delete the rxn reaction and its growth (11186-11214).

**Coverage.** Moves never remove steps, so an inhibitor-free mechanism is reached only through inhibitor-free mechanisms, at the same BFS distance. Flips and splits read the reaction only through reactant signs, and _redundant_copy_groups returns early when no inhibitor is bound. The rule sweep over rxn6 covers the same structures. What is lost is a direct run on a reaction that declares no inhibitors at all; the bi_bi_rxn aggregate pins at 9835 and 10914 and the hyperbolic sweep still exercise that path.

_Merged from:_ tests#72

#### loc/T23 — Late-file cleanups: hoist the seed_mechanisms reaction, drop the false 'shared exemplars' let, trim golden asserts the equality implies
- Lever: test_code_simplification; time saved ≈ 0.4 s; test lines removed ≈ 25; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:11959-12113, 12115-12150, 12256-12401

**Change.** - seed_mechanisms (12256): define rxn and seeds once in the outer testset, replacing the 3 copies at 12280, 12309 and 12352.
- Rate-equation dedup key (11959): remove the `let` and the false comment at 11961, move each fixture into the subtest that uses it, and bind ldh key_a once.
- Golden (12145-12149): delete `length(init) == length(golden)` and the holds_iso subset check, which `mech_keys == golden` (12144) implies.

**Coverage.** seed_mechanisms is deterministic (asserted at 12482), and no subtest mutates the seeds. The 11961 comment is provably false: each exemplar is used once, which allows deleting it under CLAUDE.md. mech_keys has one key per init mechanism, so equality to golden implies the length and the subset checks. The only loss is extra diagnostics on failure.

_Merged from:_ tests#59; tests#60; tests#68

#### loc/T24 — Dedup and canonical-construction testsets: keep the strongest of each
- Lever: delete_redundant; time saved ≈ 0.5 s; test lines removed ≈ 88; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:6741-6908, 12152-12164

**Change.** Delete these Dedup testsets:
- 6812-6824 (reversed groups);
- 6826-6834 (bi-bi distinct types);
- 6836-6842 (idempotent);
- 6844-6852 (no-op on init);
- 6854-6872 (site permutation).
Delete 'mechanism dedup via unique!' (12152-12164). Delete Test 1 of 'Mechanism — canonical by construction' (6742-6749), keeping the site-permutation Test 2. Delete the Mechanism half of '_dedup_key' (6780-6789), keeping the multiplicity-inequality half.

Keep 'expansion-path overlap' (6874) and 'permuted groups collapse' (6889).

**Coverage.** - Reversal, equality, equal hash and unique! collapse without mutation are asserted at 284-325, on a construction that is both reversed and step-swapped (stronger). The non-trivial-permutation collapse is asserted at 6889-6908.
- 239 distinct and allunique: 1504 and test_identify_rate_equation.jl:741/752-758. They imply idempotence, the no-op and distinct types.
- Site swap: 6751-6769, a stronger fixture. Equal hashes on different construction routes: 6495-6496.
- Deepcopy structural equality: every exact Set(kids) == Set(expected) comparison, plus the depth-2 `c in seen`.
The 6742-6749 and 6780-6789 deletions (about 15 lines) are my additions; they are subsumed by 284-325.

_Merged from:_ tests#33; tests#23; tests#66

#### loc/T25 — CONDITIONAL on src deletions: tests of test-only or dead src helpers go with them
- Lever: src_change; time saved ≈ 0.3 s; test lines removed ≈ 183; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_mechanism_enumeration.jl:336-366, ~46 invariant calls, 1280-1288, 4836-4864, 5152-5162, 5907-5916, 6331-6334, 6353-6366, 6531-6534, 9462-9491

**Change.** These apply only if the src items are accepted.
- _assert_mechanism_invariants deleted (src_enum#37/#3/#89): about 46 bare call lines remain after T4-T24, plus the 'ported coverage' testset (336-366). About 75 lines.
- exclude_regs kwarg removed (src_enum#49/#66): 4836-4864, 29 lines, minus about 6 to re-home T10's I-only count.
- _add_competitive_inhibitor deleted (src_enum#23/#43): 1280-1288, and pass rxn at 3388 and 11921. About 10 lines.
- _merged_site_state_assignments inlined (src_enum#45): 6353-6366, 14 lines.
- Wrong-type no-op move methods removed (src_enum#39/#94): 5152-5162, 5907-5916, 6331-6334 and the Mechanism half of 6531-6534. About 31 lines. Keep 5324, 5383 and 5524, which call the real moves.
- _re_segment_count_after_flip deleted (src_enum#21): 9462-9491, 30 lines. _testhelper_flip_groups stays, because 10255 uses it.

**Coverage.** Each deleted test exercises only the src helper that goes away. The src verifiers confirmed that each helper is test-only, dead, or always the identity on src paths.

The invariant checks are constructor-guaranteed, except 'every substrate/product appears in a step'. That check is vacuous on move children, because moves never remove steps.

The exclude_regs, no-op and _merged_site_state assertions test behaviour that no longer exists. If _assert_mechanism_invariants is moved into test/ instead of deleted (src_enum#62/#73), this item's line count drops to about 0.

_Merged from:_ src_enum#37; src_enum#3; src_enum#62; src_enum#89; src_enum#49; src_enum#66; src_enum#23; src_enum#43; src_enum#45; src_enum#39; src_enum#94; src_enum#68; src_enum#21

#### loc/T26 — CONDITIONAL on larger src rewrites: the dead-end opportunities testset and the seed-screen agreement sweep
- Lever: src_change; time saved ≈ 1.1 s; test lines removed ≈ 170; risk medium
- Sign-off: approved by Denis
- Measure first: src_enum#17: time the seed construction (ter-ter init) with and without the screen before deleting it.

**Testsets.** test/test_mechanism_enumeration.jl:44-53, 929-1018, 8207-8316

**Change.** - If src_enum#10 lands (inline _substrate_product_dead_end_opportunities into the dead-end fill): replace 929-1018 with about 15 lines through `_expand_substrate_product_dead_ends([random_topo], ter_ter_rxn)`. Assert 27 added form names in the union, and the 12 forms of the diagonal pattern. Delete _boundmap and _form_species (44-53), which become unused. This is about 60 lines beyond T3.
- If src_enum#17 lands (delete _seed_candidate_screen; give candidates the group predicates): delete the agreement testset 8207-8316 (110 lines). The group-predicate tests at 7631-7760 stay.

**Coverage.** - src_enum#10 was confirmed: the connectivity fill subsumes the explicit dead-end and mirror steps. The replacement tests the public composition through the remaining function.
- src_enum#17's equivalence was checked line by line, but its verdict is needs_measurement, because the screen exists for speed. The agreement sweep only pins that the screen equals the predicates it would be replaced by.

_Merged from:_ src_enum#10; src_enum#17

### TEST CODE SIZE of the non-enumeration test files (test_types, test_dsl, test_rate_eq_derivation, mechanism_definitions_for_test_enzyme_derivation, allosteric_ground_truth, test_allosteric_collapse/golden, reference/*, test_identify_rate_equation, test_fitting, test_compile_budget, test_accessors, test_pivot_priority_regression)

Denis, these 13,101 test lines can shrink by about 2,580 lines (20%) without losing coverage. Every deletion below names the test that already covers it. The work is in 61 deduplicated items, built from the reviewer and verifier proposals. Where those disagreed, I read the code and decided myself (I did not run Julia).

Biggest single item, DR1 (about 330 lines): replace the QSSA and stiff-ODE oracles with the BigFloat mass-action oracle that already exists. That oracle is stronger: independent Haldane values, rtol 1e-8, and it handles two-metabolite steps. This also drops OrdinaryDiffEqFIRK from the test dependencies.

Biggest clusters:
- test_types.jl: about 700 lines of duplicate construction, accessor, Step and naming tests.
- test_rate_eq_derivation.jl and the defs file: about 900 lines, mostly dead helpers, default-valued spec keywords, duplicate dependent-graph checks and duplicate reject cases.
- allosteric_ground_truth.jl: about 250 lines (unused EqualAI networks and duplicated self-validation loops).
- identify and fitting: about 380 lines.

About 123 s of time saving comes with these items. The main ones are CB2 (45 s), FI1 (25 s), DR1 (18 s) and ID1 (9.5 s). They overlap the time concern's list.

Code checks that changed the proposals:
- Single-ligand :EqualAI acceptance: I keep test_types.jl:167-172 and delete the DSL copy instead (DR13). The DSL copy's comment says it derives, which is false.
- 'Unused declared product' (test_rate_eq_derivation.jl:2120-2137): this is a third copy of the unused-reactant check, so I added it to TY12.
- DR3 removes the only caller of `_param_for_symbol`, which matches the src proposal to delete that helper. SRC1 lists the tests that go with that src change.
- DR1 plus TY2 leave the instance method `_build_wegscheider_rename_map(m::EnzymeMechanism)` (src/rate_eq_derivation.jl:144-145) with no caller. That is a src dead-code item for the src sweep.

Cross-file conflict: of the two ter-ter enumeration runs, exactly one may go. I recommend the compile_budget one (CB2), because it follows CLAUDE.md. Proposal T#0 must then not be adopted. needs_denis_signoff is true for every deletion, lowered threshold and style call.

#### loc/CB1 — test_compile_budget.jl: delete the never-called _measure_elapsed_subprocess; the wall-clock testset reuses _measure_labeled_subprocess
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 29; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_compile_budget.jl:48-64, 118-153

**Change.** Delete _measure_elapsed_subprocess (test_compile_budget.jl:48-64). Nothing in the repo calls it. In 'wall-clock: rate_equation body-build', replace the inline run/IOBuffer/regex block (136-148) with `t_first = _measure_labeled_subprocess(script, ["ELAPSED"])[1]`. Keep `@test isfinite(t_first)` and `@test t_first < RATE_EQUATION_WALLCLOCK_BUDGET_S`.

**Coverage.** Same subprocess, same script, same two assertions. A failed subprocess now returns NaN, which fails `isfinite(t_first)`; `NaN < 6.0` is false as well. Today the same failure hits `@test false`. The deleted helper has no caller anywhere (grep).

_Merged from:_ T#171; T#86

#### loc/CB2 — test_compile_budget.jl: warm the compile-reuse gate on the ping-pong bi-bi instead of a second cold ter-ter enumeration; drop the duplicate ter-ter ceiling
- Lever: cheaper_same_coverage; time saved ≈ 45 s; test lines removed ≈ 5; risk low
- Sign-off: approved by Denis
- Measure first: Run one Julia process at a time. Run the new warm script and the cold script 3 times each. Accept if max(warm/cold) ≤ 2e-3, the band the current comment records. A higher ratio means a uni-uni-only specialization was missed: report it, do not loosen the gate.

**Testsets.** test/test_compile_budget.jl:155-208 'compile reuse: ter-ter warms all of uni-uni'; subsumer test/test_mechanism_enumeration.jl:1568-1575 'init_mechanisms on ter-ter within 150 s'

**Change.** In warm_script (174-188):
- Replace the ter-ter reaction with `substrates: A[CX], B[N]; products: P[C], Q[NX]`.
- Add `GC.gc()` before the warm `@elapsed init_mechanisms(r_uni)`.
- Drop the TER_COLD print and parse.

Also:
- Delete `@test isfinite(t_ter)`, `@test t_ter < 150.0` and the ceiling comment (196-200).
- Rewrite the testset title, the 155-163 comment, the @info line and ABOUTME line 2 for the new warm-up reaction.
- Keep the warm/cold < 1e-2 gate unchanged.
- Keep test_mechanism_enumeration.jl:1568-1575 (in-process ter-ter, 250,855 pin, < 150 s). Enumeration proposal T#0, which deletes that run instead, must not also be adopted.

**Coverage.** - Ter-ter ceiling: still enforced in-process by test_mechanism_enumeration.jl:1568-1575. That test uses a structurally identical reaction (A[C],B[N],D[X] against A[C],B[N],C[O]), the same 150 s budget, and also pins the count. The only extra in the subprocess run is the cold JIT of init_mechanisms (about 1-2 s), which the bi-bi trace-compile count (≤ 100) already gates.
- Reuse gate: keeps its assertion and threshold. The ping-pong bi-bi reaches the sequential, merged, Theorell–Chance and residual paths, a superset of what ter-ter reaches (ter-ter admits no ping-pong).
- Dispatch identity (213-225) is untouched.
- Rule: this is the CLAUDE.md-compliant direction. T#0 would move the 250,855 init pin out of the enumeration file.

_Merged from:_ T#158; T#24; T#73

#### loc/CB3 — test_compile_budget.jl: time the warm uni-uni inside the trace-compile subprocess (one subprocess fewer)
- Lever: cheaper_same_coverage; time saved ≈ 4 s; test lines removed ≈ 8; risk medium
- Sign-off: not required (coverage unchanged)
- Measure first: Run the merged traced script 3 times. The trace count must stay ≤ 100 (baseline 47) and warm/cold < 1e-2. If plain bi-bi misses a uni-uni path, warm on the ping-pong bi-bi instead, as in CB2.

**Testsets.** test/test_compile_budget.jl:26-46, 92-108, 164-208

**Change.** After CB2:
- In the trace-compile script (93-104), after `init_mechanisms(r)` on bi-bi, build a uni-uni EnzymeReaction with the direct constructor. Run `GC.gc()`, then print `UNI_WARM:` from `@elapsed init_mechanisms(r_uni)`.
- Let _count_relevant_precompiles capture stdout (`run(pipeline(cmd; stdout=buf))`) and return `(n, out)`.
- Delete the separate warm subprocess.
- Keep the cold uni-uni subprocess and the rate_equation wall-clock subprocess separate.

**Coverage.** Both gates keep their assertions and thresholds. If reuse holds, the warm uni-uni compiles nothing and adds no trace lines. If reuse breaks, the ratio gate fails as before and the trace count only rises, so the trace gate gets stricter, never looser. The cost is coupled diagnosis: one regression can trip both gates. The cold and wall-clock runs stay in their own processes so that shared JIT cannot distort their timings (as the verifier noted on T#166).

_Merged from:_ T#166; T#73

#### loc/AC1 — test_accessors.jl: assert cat_allo_states on the identical mechanism in 'parameters API symmetry'; delete the rebuilt copy
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 14; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_accessors.jl:3-42, 44-58

**Change.** Add `@test EnzymeRates.cat_allo_states(am) == cat_allo_states` after test_accessors.jl:30. Delete 'added field accessors: cat_allo_states' (44-58).

**Coverage.** 44-58 rebuilds exactly what 18-30 builds: the same reaction, the same `first(init_mechanisms)` base, all-:NonequalAI tags, the same R site at multiplicity 2 and the same AllostericMechanism. It does so only to assert cat_allo_states. That assertion moves onto the first testset's `am`.

_Merged from:_ T#83; T#176

#### loc/TY1 — test_types.jl: delete five duplicate EnzymeMechanism construction/accessor testsets
- Lever: delete_redundant; time saved ≈ 0.05 s; test lines removed ≈ 106; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:73-92, 103-132, 331-350, 392-402, 566-589

**Change.** Delete these testsets in test_types.jl:
- 'EnzymeMechanism Sig repack' (73-92)
- 'EnzymeMechanism constructor' (103-132)
- 'EnzymeMechanism different orderings produce valid mechanisms' (331-350)
- 'EnzymeMechanism valid with reachable enzyme forms' (392-402)
- 'EnzymeMechanism: regulator binding' (566-589)

Move `@test length(sig) == 2` (the (reaction_sig, steps_sig) layout) into '_sig_of / _mechanism_from_sig roundtrip' near 1519.

**Coverage.** - 73-92: the same RE Michaelis–Menten as 4-31, which asserts substrates, products and n_steps. `Sig isa Tuple` is checked at 1519, and the moved `length == 2` keeps the layout check.
- 103-132: `m` repeats lines 15/27. `m_g` canonicalizes to the same type as m2 at 35-45, with the same 4-group assertion (48-49). Tuple groups written after other lines are still parsed elsewhere (1917, 2020).
- 331-350 and 392-402: the same type is built at 207-214, in the Uni-Uni spec and at test_dsl.jl:528-541. Order invariance is pinned at 1302-1304 and in test_rate_eq_derivation.jl:1427-1435.
- 566-589: mechanism 1 is 55-71 with the regulator renamed (asserted at 70); mechanism 2 is the RE MM of 4-31.

_Merged from:_ T#119

#### loc/TY2 — test_types.jl: merge the raw-Step Mechanism and converter testsets into the canonical-order and Sig-roundtrip tests; drop the rename-map round trip
- Lever: delete_redundant; time saved ≈ 0.02 s; test lines removed ≈ 64; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:1199-1232, 1278-1305, 1307-1316, 1526-1560

**Change.** - Fold 'Mechanism (non-parametric)' (1199-1232) into 'Mechanism canonicalizes step order' (1278-1305). Carry over `kinetic_groups(m) == 1:3` and `n_steps(m) == 3` (1224-1225), and add `ER.reaction(m) == r` and `m == m_perm && hash(m) == hash(m_perm)`.
- Delete 'wegscheider rename map on Mechanism' (1307-1316).
- Delete 'Mechanism <-> EnzymeMechanism converters' (1535-1560), and add `@test EnzymeMechanism(m) === em_inst` to the roundtrip testset (1526-1532).

**Coverage.** - Step set and singleton reps: implied by the `_flat_steps` assertions at 1300-1301. Equality and hash become stronger, because they now run on permuted input. kinetic_groups and n_steps on a Mechanism are kept explicitly (separate methods, src/types.jl:846-847).
- Rename map: `_build_wegscheider_rename_map(::EnzymeMechanism)` is literally the Mechanism method applied to `Mechanism(M())` (src/rate_eq_derivation.jl:142-145). Its equality with the Mechanism result therefore follows from the round trip `Mechanism(EnzymeMechanism(m)) == m` asserted at 1532.
- Converters: `EnzymeMechanism(m) = EnzymeMechanism{_sig_of(_drop_unbound_regulators(m))}()`, and the fixture has no regulators. The added `===` plus 1532 cover the converters.
- Interaction: after DR1 and this item, the instance method at src/rate_eq_derivation.jl:144-145 has no caller (src dead code).

_Merged from:_ T#120

#### loc/TY3 — test_types.jl 'step constants are named by their reaction': drop the cases other tests already pin
- Lever: delete_redundant; time saved ≈ 0.02 s; test lines removed ≈ 58; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:2101-2262

**Change.** Delete the mm_re block (2105-2115), moving its explanatory comment (2103-2104) above mm_ss. Delete the tc (2128-2138), mm_ss_backward and tc_backward (2221-2241) and re_backward (2242-2256) blocks. Keep mm_ss, re_iso, re_fused, re_fused_binding, inh, allo and pingpong.

**Coverage.** - mm_re: all four names are pinned on the public path by the Competitive Inhibitor factored strings (defs:1102-1105, test_factored_form), and on the renderer by test_types.jl:1584-1590.
- tc: the same two names are asserted on the identical tc_ss at test_rate_eq_derivation.jl:1391-1392.
- tc_backward: covered by the raw reversal loop (test_rate_eq_derivation.jl:1427-1435), which includes tc_ss and tc_backwards, plus DSL parsing of release-written steps (1216-1219, test_types.jl:2281-2283).
- mm_ss_backward and re_backward: the reversal loop asserts m2 == m for every step of Segel Uni Uni and 'Numerator: RE-chemistry'. re_iso (kept) pins K_ES_to_EP, and the parameter count is pinned by expected_n_independent_params.

_Merged from:_ T#121

#### loc/TY4 — test_types.jl: delete 'AllostericEnzymeMechanism display format' (subsumed by the shared-group display test)
- Lever: delete_redundant; time saved ≈ 0.04 s; test lines removed ≈ 32; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:591-648

**Change.** Delete test_types.jl:591-623. Add `@test count(==('\n'), s) >= 3` to 'AllostericEnzymeMechanism display: shared kinetic group' (625-648).

**Coverage.** 591-623 asserts three things:
- No 'cat_allo_states:' line: also asserted at 647.
- Inline ':: EqualAI' tags: asserted at 645-646.
- A multi-line display: moved. _format_allo_step_groups prints one newline per group (src/types.jl:1470-1477), so the 3-group fixture at 625 satisfies it.

Reg-site and NonequalAI/OnlyA tag display stay asserted at 313-328 and 537-564.

_Merged from:_ T#122

#### loc/TY5 — test_types.jl: merge the two allosteric name-chokepoint testsets, drop 'structural parameter names' and delete a provably false comment
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 62; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:1448-1490, 1614-1618, 1652-1717

**Change.** - In 1448-1490, add literal assertions on `am`: Kreg :A → :K_A_Rreg; Kd → :K_ES_to_E_S; Kiso → :K_ES_to_EP; Kon → :k_ES_to_EP. Then delete 'name(p::Kreg, m) chokepoint' (1652-1699).
- Delete 'structural parameter names' (1701-1717), but move `@test :k_ES_to_EP in ER.parameters(mm_ss)` into the mm_ss block near 2126.
- Delete the comment at 1615-1618, which says the representative is the first step's flat position, and retitle the testset that names rep_idx.

**Coverage.** name(::Kd/Kiso/Kon, ::_AnyMech) is one method for all four mechanism kinds, and Kreg has one method for AM and AEM (src/types.jl:1797-1811). After the merge, every literal from 1678-1698 is still asserted on an AllostericMechanism, along with the AEM == AM dispatch checks. The scalar names stay asserted at 1601-1611.

The Reduced-mode `:k_ES_to_EP` assertion on the all-SS MM is kept, because it pins which constant the Haldane reduction leaves fitted. The positional `k1`/`K1` negatives are tautological: names come only from the renderers, and the AST guard (2457-2470) forbids such literals in src.

The deleted comment is false: the representative is argmin `_step_priority` (src/thermodynamic_constr_for_rate_eq_derivation.jl:120-121).

_Merged from:_ T#123

#### loc/TY6 — test_types.jl: shared fixtures (ER alias at top, a raw-Step uni-uni builder, one RE Michaelis–Menten constant)
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 80; risk medium
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl 'Types' block (3-1787): 156-165, 503-512, 887-955, 1278-1405, 1407-1457, 1492-1533, 1562-1650, 1770-1786

**Change.** - Move `const ER = EnzymeRates` and `_testhelper_sp` (1789-1791) to the top of the file, and write `ER.` throughout the 'Types' block.
- Add `_testhelper_uniuni(; mults=[1], s=:S, p=:P)` returning `(; rxn, E, ES, EP, bind, iso, rel)` for the repeated reaction + Species + Step blocks.
- Add `const _testhelper_re_mm = @enzyme_mechanism ...` for the RE MM copies at 156-165, 503-512, 1407-1416 and 1448-1457.
- Keep the 4-13 copy inline, because its assertions enumerate that fixture.

**Coverage.** No assertion changes, and every testset builds the same objects. The RE MM copies are already the same singleton type, because the Sig is canonical. The count already excludes lines that TY1, TY2, TY5 and TY9 delete. Denis's inline-fixture rule is written for test_mechanism_enumeration.jl only, but this is still a style call.

_Merged from:_ T#124

#### loc/TY7 — test_types.jl: fold the Step unit testsets into 'Step: explicit consumed/released lists'
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 79; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_types.jl:815-885, 1803-1906

**Change.** Delete 'Step fields + accessors' (815-846), 'Step == / hash are structural' (848-858), 'Step canonicalizes binding direction' (860-885) and 'a step that takes up one metabolite and gives off none binds it' (1891-1906). Move their unique assertions into 1808-1839:
- fieldnames(Step) and is_equilibrium at 1808.
- `===` species identity at 1810 and 1821.
- An RE release-vs-binding twin with equality and hash.
- Hash equality at 1817.
- `_is_chemistry` for the iso (1827) and TC (1835) steps.
- `to_species == EAB` at 1832.
- The forward fused binding `Step(EA, EPQ, [B], [], true)`: bound_metabolite == B and _is_chemistry.

**Coverage.** Every assertion is kept, or strengthened by comparing a release-written step with its binding twin instead of two identically built steps. These are pure runtime unit tests with no generated calls.

_Merged from:_ T#125

#### loc/TY8 — test_types.jl: delete the duplicate PK :NonequalAI synth-dep testset
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 50; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:1719-1768; subsumer test/test_rate_eq_derivation.jl:2333-2382

**Change.** Delete 'synth-dep I-state names consistent with chokepoint (NonequalAI)' (test_types.jl:1719-1768).

**Coverage.** test_rate_eq_derivation.jl:2333-2382 builds the byte-identical @allosteric_mechanism and makes the same three assertions: some name starts with K_I_, no name ends in _T, and rate_equation is finite. Only the RNG seed differs. It is the same @generated specialization, so the suite's first allosteric rate_equation JIT just moves to that file.

_Merged from:_ T#126; src_types#1 (test impact)

#### loc/TY9 — test_types.jl: small collapses (RegulatorySite pair, constructor tautologies, hand-built AEM accessors, Kd loop, two in-file duplicates, two re-lifts)
- Lever: delete_redundant; time saved ≈ 0.15 s; test lines removed ≈ 36; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:134-154, 652-679, 764-813, 914-955, 1417-1420, 1644-1649, 2653-2699

**Change.** (a) Merge 804-813 into 'RegulatorySite: validation + accessors' (764-802), keeping one 'all four states accepted' loop.
(b) Delete the four isa tautologies at 658, 661, 663 and 666; keep the supertype checks.
(c) Replace 'AllostericEnzymeMechanism struct + accessors' (134-154), a hand-built type that bypasses validation, with four accessor assertions on the validated aem at 1417-1420: catalytic_mechanism, multiplicity, cat_allo_state for g = 1:3, and regulatory_sites. Its tags (:EqualAI/:NonequalAI/:OnlyA) are distinct, so an indexing error is still caught.
(d) Fold Kd into the `for T in (Kiso, Kon, Koff, Kfor, Krev)` loop (931-938), keeping its equality and inequality lines.
(e) Delete the group-2/3 name assertions at 1646-1649.
(f) Delete the re-lift lines `ER.AllostericMechanism(cube/ordered) isa …` at 2681 and 2698.

**Coverage.** - (a), (d): every assertion survives in merged form.
- (b): each value was just built by that constructor, so the check cannot fail.
- (c): the same accessors are asserted on a validated instance whose distinct tags still catch a group-index error.
- (e): step_c and step_d are the same steps as step2 and step3, whose names are pinned at 1586 and 1590.
- (f): @allosteric_mechanism already ran the identical AllostericMechanism constructor, Stiemke guard included (src/dsl.jl:1248-1251), and 2680/2697 assert that it accepted. Re-lift idempotence is covered for every allosteric spec at test_rate_eq_derivation.jl:1436-1446.
- The single-ligand :EqualAI acceptance at 167-172 is deliberately kept (see DR13).

_Merged from:_ T#127; T#137; T#117

#### loc/TY10 — test_types.jl iso canonicalization: drop the cases the step-reversal test already proves
- Lever: delete_redundant; time saved ≈ 0.01 s; test lines removed ≈ 25; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:1234-1276

**Change.** Delete the Tier 1 block (1235-1252) and the s1 mechanism with its equality check (1257-1262, 1269). Keep s2 and the F→E direction assertions (1270-1275), and trim the 1254-1256 comment to the direction claim.

**Coverage.** m_fwd is the Segel Uni Uni spec renamed, and s1 is exactly the Segel Iso Uni Uni spec. test_rate_eq_derivation.jl:1427-1435 reverses every step, and random subsets, of both specs and asserts m2 == m, which proves both tier equalities. The canonical F→E direction, which that loop does not pin, stays at 1270-1275.

_Merged from:_ T#128

#### loc/TY11 — test_types.jl bottomless-RE-segment testset: build the rejected mechanism once; drop the duplicate ping-pong accept case
- Lever: delete_redundant; time saved ≈ 0.01 s; test lines removed ≈ 24; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:404-501

**Change.** Delete the @test_throws at 410-421. Turn the try/catch at 422-438 into `catch e; e end` and assert `err isa ErrorException` plus the two message checks on `sprint(showerror, err)`. Delete m_pingpong (488-500), moving its two-line rationale (486-487) next to the identical pingpong at 2206.

**Coverage.** The rejection is still asserted for both error type and message, on the same mechanism. The ping-pong accept case is the identical mechanism built at test_types.jl:2206-2217 (names asserted) and at test_rate_eq_derivation.jl:1293-1303 (d_free asserted). A wrong rejection fails both.

_Merged from:_ T#129

#### loc/TY12 — Delete two copies of the 'declared reactant appears in no step' check (test_types, test_rate_eq_derivation) and a false comment
- Lever: delete_redundant; time saved ≈ 0.02 s; test lines removed ≈ 40; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:352-390; test/test_rate_eq_derivation.jl:2120-2137; subsumer test/test_mechanism_enumeration.jl:349-364

**Change.** Delete the unused-substrate case in test_types.jl 'EnzymeMechanism error cases' (362-379). Delete the false first half of the NOTE (381-384), which says duplicate reactions in distinct kinetic groups are valid; keep its still-true unreachable-forms half. Delete the 'Stoichiometric infeasibility' block in test_rate_eq_derivation.jl 'Allosteric edge cases' (2120-2137, product Q never in a step).

**Coverage.** test_mechanism_enumeration.jl:349-364 asserts ErrorException from _assert_mechanism_invariants on a mechanism whose declared substrate T is never bound. The check is one loop over `(substrates..., products...)` (src/mechanism_enumeration.jl:3374-3377), so the substrate and product copies take the identical path. The comment is false: the constructor calls _assert_each_reaction_once (src/types.jl:829), and test_types.jl:2316-2324 asserts that the duplicate throws. Interaction: the enumeration testset at 349-364 must stay. Enumeration proposal T#15 deletes 327-334, a different testset.

_Merged from:_ T#130

#### loc/TY13 — Move the fused-plus-plain-binding fitted count from test_types.jl next to its mass-action case
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 12; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_types.jl:1908-1923; test/test_rate_eq_derivation.jl:1279-1284, 1332-1336

**Change.** Delete 'a fused and a plain binding of one metabolite share a kinetic group' (test_types.jl:1908-1923). In test_rate_eq_derivation.jl 'fused steps derive the mass-action rate', add the `length(fitted_params(..)) == 5` count with its one-line rationale. Look the case up by content or a named binding, not `_testhelper_fused_cases[end]`.

**Coverage.** It is the identical mechanism, the last fused case, and test_rate_eq_derivation.jl:1332-1336 already proves it derives the mass-action rate. The count assertion moves rather than disappears, and it is the same fitted_params instance.

_Merged from:_ T#131

#### loc/TY14 — test_types.jl: delete 'AllostericMechanism rejects an unsatisfiable Haldane'
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 13; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:2598-2609

**Change.** Delete test_types.jl:2598-2609.

**Coverage.** Two existing tests cover it:
- The verdict for this exact tagging (S :OnlyA, catalysis :EqualAI, P :EqualAI, on the same three steps) is asserted at test_types.jl:2499.
- The constructor's only verdict-to-error path (src/types.jl:922-924) is exercised by the @allosteric_mechanism @test_throws at 2619-2639.

_Merged from:_ T#132

#### loc/TY15 — test_types.jl: delete 'EnzymeReaction struct + accessors'
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 17; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:1140-1156

**Change.** Delete test_types.jl:1140-1156.

**Coverage.** Each assertion is made elsewhere:
- substrates/products on a constructor-built reaction: test_types.jl:1021-1023.
- reactants: test_dsl.jl:278-281.
- allowed_catalytic_multiplicities with several values: test_dsl.jl:482-487; default and oligomeric_state: 466-479.
- regulators count: test_types.jl:1036-1037.
- RegulatorMults accessors: 1109-1111.

_Merged from:_ T#135

#### loc/TY16 — test_types.jl: drop the 44-spec kon/koff/Kiso prefix sweep
- Lever: delete_redundant; time saved ≈ 0.25 s; test lines removed ≈ 5; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:2257-2261

**Change.** Delete the `for spec in MECHANISM_TEST_SPECS … @test !any(n -> occursin(r"^(kon|koff|Kiso)_", n), names)` loop and its comment (test_types.jl:2257-2261).

**Coverage.** Every Full-mode name goes through name(p, m). The renderers (src/types.jl:1797-1817) can emit only k_/K_ + tag, the Kreg form or the three scalars, and each renderer's literal output is pinned:
- Kd/Kon/Koff/Kfor/Krev/Kiso with :None and :I: 1584-1598.
- Kreg: 1488-1489 and 1680-1683.
- The full PARAMS_FULL bytes of every allosteric spec: test_allosteric_golden.jl:17.

The AST guard forbids other literals in src. The implicit 'Full enumeration does not throw' smoke check survives, because the derivation calls _enumerate_parameters_full for every spec.

_Merged from:_ T#116

#### loc/DS1 — test_dsl.jl: remove ten duplicates of DSL behaviour asserted elsewhere
- Lever: delete_redundant; time saved ≈ 0.1 s; test lines removed ≈ 170; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_dsl.jl:76-107, 130-160, 205-218, 284-293, 432-464, 543-632

**Change.** Delete or merge in test_dsl.jl:
(1) 76-89 '+ step-side syntax'.
(2) 91-107 'multi-product balances placeholders'.
(3) The m_call and m_multi blocks (130-160) of the decomposed-Species grammar test.
(4) The positive half of '@allosteric_mechanism (parsing & validation)' (205-218); keep the four rejections.
(5) spec2 in '@enzyme_reaction' (284-293).
(6) Fold 432-442 into 444-464 with `@test occursin("opaque bound-form name", msg)`.
(7) The numeric spot check and the 6-state m2 (543-562).
(8) Everything after the no-enzyme-form rejection in 'Elementary steps' (578-596).
(9) 'No-atom species' (599-611).
(10) 'Constraint DSL parsing' (613-632).

Also reword the temporal 'New form' comment at 76 wherever it survives.

**Coverage.** Each deleted piece and its subsumer:
(1) and (9) build the same mechanism as 528-541, which makes the same and more assertions.
(2) The Segel Ordered Uni Bi spec declares the same `substrates: A; products: P, Q` through the same parser; bi-uni is covered at 44-56.
(3) test_types.jl:4-31 (call-form species), test_dsl.jl:57-65 (multi-bound parse) and test_types.jl:280 (EPQ).
(4) test_types.jl:313-328 builds the identical F6P/I::OnlyI mechanism, and 202 asserts the regulators exactly.
(5) Multi-digit formulas: test_mechanism_enumeration.jl:228-237 and 373-374; competitive_inhibitors: test_dsl.jl:316-338.
(6) The merged test asserts a superset.
(7) The Uni-Uni spec's analytical oracle plus test_rate_eq_derivation.jl:1316-1317; the F(P) call form is covered by Segel Ping Pong.
(8) test_types.jl:283-293 m_reg.
(10) test_dsl.jl:115-128, test_types.jl:33-52 and 2009-2025.

_Merged from:_ T#118; src_dsl#19

#### loc/DS2 — test_dsl.jl: delete 'fused catalytic release: metabolite in neither bound list' and the positive halves of 'several metabolites on a step side'
- Lever: delete_redundant; time saved ≈ 0.03 s; test lines removed ≈ 54; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_dsl.jl:705-761

**Change.** Delete the testset at test_dsl.jl:705-724. In 726-771, delete only_transformation, the tc block and the two_in block (727-760). Keep the `@test_throws "more than one enzyme-form term"` rejection.

**Coverage.** - 705-724: the mechanism is the Segel Ordered Uni Bi spec (defs:220-228), built by the same parser and checked against Segel IX-60. Storing a fused release as the binding it reverses is pinned on the DSL path at test_types.jl:2078-2099 and at Step level at 1829-1832. The reactions() layout is pinned at test_types.jl:21-25.
- tc: identical to tc_ss, whose from_species and consumed are asserted at test_rate_eq_derivation.jl:1386-1390. The k_EA_B_to_EQ_P / k_EQ_P_to_EA_B names at 1391-1392 pin `released == [P]` through _forward_sides.
- two_in: its first step `E + A + B <--> E(A, B)` is the first step of _testhelper_tc_cases.two_in, whose mass-action check (1383-1385) fails on any mis-parsed consumed list. Sorted consumed lists are pinned at Step level (test_types.jl:1836-1838).

_Merged from:_ T#133; T#134

#### loc/DR2 — test_rate_eq_derivation.jl helpers: derive the dependent parameters once per spec and evaluate them without eval(Meta.parse)
- Lever: cheaper_same_coverage; time saved ≈ 4 s; test lines removed ≈ 15; risk low
- Sign-off: not required (coverage unchanged)
- Measure first: Warm `@elapsed EnzymeRates._dependent_param_exprs(typeof(spec.mechanism))` per spec, and time 'Enzyme Derivation Tests' before and after. Confirm on the first run that the allosteric combined-solve dependent Exprs contain only :*, :/ and :^ calls. The saving drops to about 2 s if the time concern also cuts the draw counts (T#78).

**Testsets.** test/test_rate_eq_derivation.jl:513-599 (helpers), call sites 868-893, 960-986, 895-934, 1030-1048

**Change.** - Replace `_get_independent_params(m)` with `EnzymeRates.fitted_params(m)`; the @generated method returns exactly `indep`.
- Compute `dep_exprs, indep = _dependent_param_exprs(typeof(m))` once per spec in run_all_tests, or in a type-keyed IdDict cache, and pass it down.
- Delete _get_dependent_params and _eval_dep_expr. Evaluate the dependent Exprs with a 4-line evaluator: Symbol → p[x], Real → x, `:call` → `getfield(Base, op)(args...)`.
- Hoist the per-draw Mechanism/AllostericMechanism decompile in analytical_oracle_params/positional_params to once per spec.

**Coverage.** Only test-helper plumbing changes. The oracles get identical values: the same expressions, evaluated at the same draws, minus the string round trip. The package derivation is still exercised uncached by test_structure, test_constraint_counting, rate_equation_string and every @generated body. No assertion is removed.

_Merged from:_ T#75; T#162

#### loc/DR3 — test_rate_eq_derivation.jl: fold the two group-rep invariance testsets into 'reversing written steps' as a within-group-reversal variant
- Lever: delete_redundant; time saved ≈ 0.4 s; test lines removed ≈ 84; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_rate_eq_derivation.jl:1395-1439, 2384-2469

**Change.** Extend `variants(groups)` (1400-1401) with a third variant, `[reverse(g) for g in groups]`, so both the EnzymeMechanism and allosteric loops assert `m2 == m` for it. Delete 'kinetic-group name rep is structurally primary' (2384-2401), the invariance comment block, the helpers _dep_struct_key and _dep_struct_key_set (2403-2452) and 'dependent-param choice invariant to group-rep' (2454-2469).

**Coverage.** _canonical_group_order! sorts the steps of each group (src/types.jl:636-641), so `Mechanism(rxn, [reverse(g) …]) == mech` exactly. Both deleted testsets therefore compare a mechanism with itself. The new variant fails if within-group canonicalization is lost, and it also runs on the allosteric specs.

The name-rep assertion (K_ES_to_E_S, never K_EI1inhS_to_EI1inh_S) is pinned by spec 29's factored strings (defs:1733-1736).

Interaction: _dep_struct_key is the only caller of `_param_for_symbol` anywhere. That matches the src proposals deleting it (src_types#1/#87); see SRC1.

_Merged from:_ T#82

#### loc/DR4 — test_rate_eq_derivation.jl: run the dependent-graph soundness check once per spec in test_constraint_counting; delete the subset loops and the §5a perf duplicate
- Lever: delete_redundant; time saved ≈ 0.2 s; test lines removed ≈ 27; risk low
- Sign-off: approved by Denis
- Measure first: One focused run confirming that clauses (b) and (c) hold on the 31 non-allosteric specs. A failure is a real finding.

**Testsets.** test/test_rate_eq_derivation.jl:823-845, 1510-1553, 1555-1612

**Change.** Change `_dep_graph_is_sound` to take `(dep, indep)`, move it and rhs_syms above test_constraint_counting, and add `@test _dep_graph_is_sound(dep_exprs, indep)` there, reusing the call at 826. Move the LDH-leak rationale (1536-1540) there as a comment.

Then delete:
- 'indep ∩ keys(dep) == ∅ (allosteric i-state mechanisms)' (1535-1546).
- '… (all MECHANISM_TEST_SPECS)' (1548-1553).
- The MECHANISM_TEST_SPECS half of 'allosteric dependent-param graph is sound' (1608-1611); keep the reproducer loop.
- The §5a finite-at-1.5 check (1518-1519) and the §5a perf block (1526-1531).

Keep the §5a products=0 finiteness and _undefined_rhs_symbols checks.

**Coverage.** - _LDH_ISTATE_MECHS (1481-1482) is a filter of MECHANISM_TEST_SPECS, so 1535-1546 is a strict subset of 1548-1553.
- Clause (a) of _dep_graph_is_sound is that same intersection test, so running it on every spec subsumes both loops and adds the acyclic and closed checks for the 31 non-allosteric specs.
- The §5a perf block calls the same helper with the same 0-alloc / 120 ns bounds on the same three specs that test_performance gates in run_all_tests.
- The finite-at-1.5 evaluation compiles the same body as the kept products=0 evaluation, so an UndefVarError still surfaces.

_Merged from:_ T#85; T#103; T#173; T#114

#### loc/DR5 — Delete dead test helpers and explicit default-valued spec keywords; collapse a self-duplicating assertion loop; fold the pivot-regression file in
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 59; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_rate_eq_derivation.jl:534-548, 1769-1778; test/mechanism_definitions_for_test_enzyme_derivation.jl:639-2581 (listed lines); test/test_pivot_priority_regression.jl; test/runtests.jl

**Change.** Deletions:
- make_independent_params (test_rate_eq_derivation.jl:534-548, no caller).
- _spec_by_name, pfk_mechanism, pfk_rate_analytical, hk_mechanism and hk_rate_analytical (defs:2570-2581, no caller).
- The five `analytical_rate_fn=nothing` lines (defs:639, 669, 702, 732, 1055).
- The twelve `expected_factored_num/denom=nothing` lines (defs:2001-2002, 2086-2087, 2219-2220, 2319-2320, 2405-2406, 2471-2472). Keep `analytical_kcat_fn=nothing` at 2404, whose comment explains it.

Replacement: swap the build_power_expr loop (1769-1778) for `@test bpe(R(1), [(:k1f, R(1))]) isa Union{Int, Symbol, Expr}`.

Optional: move the single testset of test_pivot_priority_regression.jl into test_rate_eq_derivation.jl, dropping its module wrapper, `using`, `const ER`, ABOUTME and the runtests include.

**Coverage.** - Dead helpers and default-valued keywords carry no assertions; grep confirms each helper appears only at its definition.
- In the build_power_expr loop, four of five entries are already asserted more strictly at 1755-1767 (=== :k1f, isa Expr, === :Keq, isa Expr). Only (1,[k1f]) was unique, and it is kept.
- The pivot regression moves unchanged.

_Merged from:_ T#86; T#174; T#107

#### loc/DR6 — MechanismTestSpec: pin metabolite names exactly and drop expected_n_metabolites; fix the HK, PFK-1 and PK lists to canonical order
- Lever: cheaper_same_coverage; time saved ≈ 2.5 s; test lines removed ≈ 43; risk low
- Sign-off: approved by Denis
- Measure first: One focused run: `metabolites(m) == Tuple(spec.metabolite_names)` for all 43 specs. Any other mismatch means fixing that list's order.

**Testsets.** test/test_rate_eq_derivation.jl:797-809; test/mechanism_definitions_for_test_enzyme_derivation.jl struct + 43 `expected_n_metabolites=` lines

**Change.** Reorder the hand-written metabolite_names of HK (defs:2204), PFK-1 (2071) and PK (2305) to metabolites(m) order (sorted substrates, sorted products, then regulators). In test_structure replace `@test length(metabolites(m)) == spec.expected_n_metabolites` with `@test metabolites(m) == Tuple(spec.metabolite_names)`. Delete the expected_n_metabolites field and its 43 spec lines.

**Coverage.** Exact equality of the metabolite tuple implies the count, so the deleted field is redundant. It also pins the names and order, which the count never did (the verifier's correction of T#79, which would have lost the only name pin).

Today the three out-of-order lists give the concs NamedTuple a second type. test_kcat_rescaling and test_zero_metabolite_finite use metabolites(m), so each of those allosteric specs compiles @generated rate_equation twice; one order removes the duplicate compile. By-name concs access stays covered by 'structural names + synth-dep routing' (2367-2373), which passes a PEP-first non-canonical order.

_Merged from:_ T#79

#### loc/DR7 — MechanismTestSpec: default the mirror and Wegscheider constraint counts to 0 (style call)
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 75; risk low
- Sign-off: approved by Denis

**Testsets.** test/mechanism_definitions_for_test_enzyme_derivation.jl:22-28 and 75 `=0,` lines

**Change.** Give expected_n_mirror_constraints and expected_n_wegscheider_constraints a @kwdef default of 0 and delete the 40 + 36 zero-valued keyword lines.
- The commented zero at defs:1553 goes with DR8.
- The comment at defs:1997 ('site-independence constraints are in param_constraints of CM') is provably false, because no `param_constraints` exists in src, so its line can go.

**Coverage.** test_constraint_counting (test_rate_eq_derivation.jl:840-843) still asserts both counts on every spec, and a default of 0 asserts exactly what the deleted lines assert. The struct comment already calls the mirror count 'largely vestigial'. Explicit zeros do document intent, which is why this is Denis's style call.

_Merged from:_ T#94

#### loc/DR8 — Delete spec 'MWC Dimer + Independent Inhibitor', an exact duplicate of 'Homodimer + Non-competitive Inhibitor', and its 4 golden lines
- Lever: delete_redundant; time saved ≈ 0.4 s; test lines removed ≈ 62; risk low
- Sign-off: approved by Denis

**Testsets.** test/mechanism_definitions_for_test_enzyme_derivation.jl:1504-1564; test/reference/allosteric_golden_reference.txt:9-12

**Change.** Remove the 26B `let` block (defs:1504-1564) and its golden reference lines 9-12 (test/reference/allosteric_golden_reference.txt) in the same commit. Fold the 26B header sentence (Wegscheider closure is automatic) into the surviving spec's comment so both readings stay documented. Denis picks the surviving name.

**Coverage.** The two @allosteric_mechanism_src bodies (defs:1450-1463 and 1511-1524) are identical, so both build the same singleton type. The oracle bodies and every expected_* field are identical too.

Every per-spec check, the golden test, the reversal loop, the dependent-graph check, test_types.jl:2258 and test_mechanism_enumeration.jl:7396 therefore lose only a repeat iteration on the same type. No test counts MECHANISM_TEST_SPECS.

_Merged from:_ T#81; T#96

#### loc/DR9 — Drop the hand-written factored-form strings of the allosteric dimer specs; the golden file pins them byte for byte
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 12; risk medium
- Sign-off: approved by Denis

**Testsets.** test/mechanism_definitions_for_test_enzyme_derivation.jl:1430-1435, 1495-1500; subsumer test/test_allosteric_golden.jl:9-31

**Change.** After DR8, remove expected_factored_num/denom from 'MWC Dimer' (defs:1430-1435) and 'Homodimer + Non-competitive Inhibitor' (1495-1500). Non-allosteric specs keep theirs.

**Coverage.** test_factored_form slices the v-line of `rate_equation_string(m)` (Reduced is the default). test_allosteric_golden.jl compares that same string byte for byte for every allosteric spec, and a script check confirmed that each golden v-line contains exactly these num/denom strings.

What is lost: an independently hand-written expectation, which matters only if someone regenerates the golden file without review.

_Merged from:_ T#107; T#93; T#172

#### loc/DR10 — test_rate_eq_derivation.jl 'Rate equation too large error': keep one of the two @test_throws that exercise the identical guard
- Lever: delete_redundant; time saved ≈ 1 s; test lines removed ≈ 32; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_rate_eq_derivation.jl:1957-2018

**Change.** Delete m_cyclic and its comment (1986-2017). Keep m_manual and its `@test_throws "polynomial terms"`.

**Coverage.** Both mechanisms are all-SS and reach the same function by the same path: rate_equation_string → _raw_symbolic_rate_polys → _assert_derivable (src/rate_eq_derivation.jl:461-467) → _segment_graph_terms with G > 1. That function has no shape-dependent branch. The _bareiss_det zero-pivot and singular paths are tested directly (1946-1955), _eq_complexity values are pinned on 4 shapes (1881-1944), and the guard is also covered end to end by test_identify_rate_equation.jl 'derivation guard'. The m_cyclic 'reachable two ways' history predates the upfront V×τ guard.

_Merged from:_ T#90

#### loc/DR11 — test_rate_eq_derivation.jl: keep one of the two 'all-RE catalytic cycle raises' testsets
- Lever: delete_redundant; time saved ≈ 0.3 s; test lines removed ≈ 25; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_rate_eq_derivation.jl:1830-1879

**Change.** Delete 'Numerator: all-RE catalytic cycle raises' (1830-1853). Keep the activator-route testset (1855-1879), whose cycle does not pass through free E.

**Coverage.** Both assert an ErrorException containing 'no finite rate'. Its only source is the RE-ratio consistency loop in _compute_alpha (src/rate_eq_derivation.jl:282-291), which checks every RE step with no branch on cycle topology, and both reach it through the same rate_equation_string path. The 'Numerator:' prefix is historical; there is no separate numerator-side check any more. The enumeration-side _re_turnover_cycle predicate has its own tests.

_Merged from:_ T#91

#### loc/DR12 — test_rate_eq_derivation.jl: fold 'allosteric reproducers: rate_equation is callable' into the detailed-balance testset
- Lever: delete_redundant; time saved ≈ 0.03 s; test lines removed ≈ 13; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_rate_eq_derivation.jl:1614-1663; fixtures test/reference/allosteric_undefvar_reproducers.jl

**Change.** Delete 1614-1631 and move its rationale (1615-1620: a structurally sound dependent graph can still leave a forward reference) into the detailed-balance comment at 1634-1639.

**Coverage.** The detailed-balance testset calls rate_equation 100 times on the same three reproducer types, with the same params/concs key sets and Float64 values. An UndefVarError in the straight-line @generated body errors that testset. @test_nowarn adds nothing: src has no @warn/@info/@debug in the derivation path (grep). The 3-type compile moves, so only lines are saved.

_Merged from:_ T#99; T#92

#### loc/DR13 — 'Allosteric edge cases': delete the three :OnlyI rejections, the single-ligand :EqualAI acceptance and the empty-ligand rejection that test_types and test_dsl already make
- Lever: delete_redundant; time saved ≈ 0.1 s; test lines removed ≈ 64; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_rate_eq_derivation.jl:2068-2118, 2245-2259

**Change.** In test_rate_eq_derivation.jl 'Allosteric edge cases', delete:
- The three `@test_throws Exception eval(:(@allosteric_mechanism … :: OnlyI …))` blocks (2068-2103).
- The single-ligand `I::EqualAI` acceptance (2106-2118). Its comment says the derivation runs, but the test never derives.
- The cm_simple / empty-ligand reg-site block (2245-2259).

**Coverage.** - :OnlyI: rejection is one loop over all catalytic tags in the AllostericMechanism inner constructor (:OnlyI ∉ _VALID_CAT_ALLO_STATES, src/types.jl:866), with no step-kind or position branch. It is asserted via the DSL on the catalytic step (test_dsl.jl:232-242), via the 3-arg constructor on the iso group (test_types.jl:174-177), and via the raw constructor on an S-binding group (test_types.jl:1386-1389).
- Single-ligand :EqualAI: the identical site is accepted at test_types.jl:167-172 (kept), and RegulatorySite has no ligand-count or :EqualAI branch (src/types.jl:127-140). DSL parsing of `::EqualAI` regulators is derived for real in the PFK-1, HK and m_all specs (defs:2015, 2113, 2336).
- Empty ligand: test_types.jl:521-523 is the identical 3-arg call.

_Merged from:_ T#111; T#112

#### loc/DR14 — 'Allosteric edge cases': drop the m_mixed equilibrium check (same catalytic core as the collapse test)
- Lever: delete_redundant; time saved ≈ 0.15 s; test lines removed ≈ 27; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_rate_eq_derivation.jl:2172-2198; subsumers test/test_allosteric_collapse.jl:39-46 and m_ro at test_rate_eq_derivation.jl:2200-2244

**Change.** Delete test_rate_eq_derivation.jl:2172-2198.

**Coverage.** test_allosteric_collapse.jl:39-46 builds the identical catalytic core: uni dimer, S binding :NonequalAI RE, SS chemistry :EqualAI, P binding :EqualAI RE. It asserts |v| < 1e-8 at equilibrium, K_I_ES_to_E_S absent from fitted_params, and the explicit mirror in the string.

Collapse together with a regulator is still covered by m_ro in the same testset: the same `(((:I,), 2, (:NonequalAI,)),)` site, a forbidden-split collapse, and v = 0 at equilibrium with I = 0.5. The regulator's equilibrium property is also asserted per spec by 'Haldane Equilibrium' on the regulated dimer specs.

_Merged from:_ T#98

#### loc/DR15 — test_rate_eq_derivation.jl: fold test_factored_form into test_rate_equation_string (render once)
- Lever: test_code_simplification; time saved ≈ 0.3 s; test lines removed ≈ 5; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_rate_eq_derivation.jl:960-1028, 1185-1186

**Change.** Move the body of test_factored_form (1012-1028) into test_rate_equation_string so it reuses `s` (line 964). Delete the separate function and its call at 1186.

**Coverage.** The same string and the same assertions, computed once instead of twice for the specs that carry factored strings.

_Merged from:_ T#84

#### loc/DR17 — _eval_rate_string: evaluate the whole rendered string (constraint lines + v) from the fitted parameters only
- Lever: test_code_simplification; time saved ≈ 0.5 s; test lines removed ≈ 4; risk medium
- Sign-off: not required (coverage unchanged)
- Measure first: One full run. A failure is a real string defect to report to Denis, not a reason to revert.

**Testsets.** test/test_rate_eq_derivation.jl:773-793, 960-986

**Change.** In _eval_rate_string (773-793):
- Strip `EnzymeRates.ANNOTATION_SUBSTITUTED`, as _undefined_rhs_symbols does (1492).
- Keep the constraint lines with the v line, and evaluate them in a `let` that binds only `new_params` and `concs`.
- In test_rate_equation_string, pass `new_params` instead of `all_params` (981), so the string oracle no longer needs compute_all_params.

**Coverage.** Stronger than today. Today `eq_line = last(split(code_body, "v = "))` discards the rendered Haldane and Wegscheider lines and feeds the v line the package's own dependent values, so a mis-rendered constraint line is never executed. With the change, a wrong constraint line, or a v line that references an undefined dependent, fails the check. Reduced mode renders constraint lines through a different call than the v line (src/rate_eq_derivation.jl:669-689), so this is a real gap.

_Merged from:_ T#95

#### loc/DR18 — test_zero_metabolite_finite: delete the trap_mets exemption if no spec triggers it, and drop the false trap sentence from its docstring
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 13; risk low
- Sign-off: approved by Denis
- Measure first: In one focused session, for each AllostericEnzymeMechanism spec run the code at 1151-1161 and print the spec name and trap_mets whenever it is non-empty.

**Testsets.** test/test_rate_eq_derivation.jl:1122-1171

**Change.** Measure whether trap_mets (1151-1161) is ever non-empty.
- If it is empty for every spec, delete the exemption block and the `!(zeroed in trap_mets)` clause (1168), and rewrite the docstring (1122-1137) without the trap case.
- If it is non-empty, keep the code and fix only the false docstring sentences. AG3 fixes the gate citation.

**Coverage.** If trap_mets is empty for every spec, the exemption never changes which assertions run, so deleting it changes nothing. A future trap would then surface as a failing `v != 0` instead of being silently exempted. allosteric_ground_truth.jl:571 records that the trap→0 was fake.

_Merged from:_ T#113

#### loc/AG1 — allosteric_ground_truth.jl: delete the three all-:EqualAI networks that no gate uses, and fold the remaining L=0 harness checks into their gates' draw loops
- Lever: delete_redundant; time saved ≈ 0.05 s; test lines removed ≈ 108; risk low
- Sign-off: approved by Denis

**Testsets.** test/allosteric_ground_truth.jl:52-160, 167-226, 257-350

**Change.** - Delete uni_equalAI_flux (52-74), multi_equalAI_flux (100-125) and metab_dfree_equalAI_flux (257-282), with their headers and their (b) assertions (137-141, 155-159, 308-312).
- Move each remaining (a) assertion (the :OnlyA oracle at L=0 equals its closed form or base network) into the matching gate's draw loop, as the inert gate already does at 972-975. Then delete the three emptied self-validation testsets (127-142, 145-160, 299-313).
- Keep the gate names; defs:2561 cites one.
- Edit the comment at 285, which becomes false.

**Coverage.** The three networks are called only by their own (b) checks (grep) and are never a gate's oracle.
- metab_dfree_equalAI_flux is the same 8-species network as biuni_nonequalAI_flux with k_A = k_I, whose 5-draw check at 432-438 is a superset.
- The shared solver stays anchored at 130, 148, 719-720, 875-881 and 1082-1092.
- Each moved (a) assertion still compares the oracle with an independent closed form, never with derivation code. It now runs at 5 random draws instead of one fixed point.

The only loss is diagnostic: a broken oracle and a broken derivation now fail in the same testset.

_Merged from:_ T#100; T#104

#### loc/AG2 — allosteric_ground_truth.jl: per-form-flip reference as the oligomer oracle at nprot=1; one self-validation loop; the n=1 :NonequalAI gate joins the ^n gate
- Lever: delete_redundant; time saved ≈ 0.02 s; test lines removed ≈ 95; risk medium
- Sign-off: approved by Denis

**Testsets.** test/allosteric_ground_truth.jl:352-546, 1016-1167

**Change.** (1) Delete biuni_nonequalAI_flux (352-387) and use `biuni_mwc_oligomer_flux(1, …; freeflip=false)` in its place.
(2) Delete the discriminator at 494-503, which is argument for argument identical to 1117-1118.
(3) Move the per-form (a)/(b) checks into the oligomer self-validation via `for nprot in (1,2,3), freeflip in (true,false)`.
(4) Merge the two bi-uni self-validation loops (397-440 and 464-492), which draw from the same seed in the same order.
(5) Run the :NonequalAI gate as `for nprot in (1, 2, 3)` against the oligomer oracle, removing the n=1 copy (514-546). Build the three fixtures with a literal ternary, not @eval.

Keep biuni_nonequalAI_freeflip_flux and the 1115-1116 formulation-1 pin.

**Coverage.** - Network identity was checked line by line. At nprot=1 the oligomer species order and per-protomer edges, the Haldane kr formula and the flipping of every form under freeflip=false match 366-386, so the deleted per-form assertions reappear for nprot 1-3 as a superset.
- The merged loops reproduce the same parameter points exactly.
- Gate coverage: the same mechanism types, the same fitted_params pin (525-526 = 1149-1150), the same rtol 1e-4 and 6 draws per n.
- The independent formulation-1 transcription is deliberately kept (T#102 part 1 is rejected), so no oracle is left without an independent cross-check, in line with the policy at 679-685.

_Merged from:_ T#101; T#102

#### loc/AG3 — Fold the 'metabolite in D[g_free] at zero concentration' testset and two uni-OnlyA testsets into the gates that declare the same mechanisms; fix the false docstring
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 40; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/allosteric_ground_truth.jl:167-192, 321-350, 514-591; test/test_rate_eq_derivation.jl:1131-1136, 2560-2580, 2596-2607

**Change.** (1) Move the metabD B=0 block (allosteric_ground_truth.jl:565-572) to the end of the metabD gate, and the :NonequalAI B=0 block (581-590) to the end of the :NonequalAI gate (the nprot == 1 iteration if AG2 lands). Carry the 548-555 rationale and delete the testset (548-591).
(2) Move 'D[g_free] surfaced per allosteric state' (test_rate_eq_derivation.jl:2560-2580) and 'kcat consistent … (uni-OnlyA)' (2596-2607) into the OnlyA gate (allosteric_ground_truth.jl:167-192). Keep the 2561-2566 comment.
(3) Rewrite test_rate_eq_derivation.jl:1131-1136, which wrongly says the gate checks a trap → 0 (allosteric_ground_truth.jl:571 records that the trap was fake).

**Coverage.** The same assertions run on the same mechanism types with the same fixed parameter values:
- At B=0: oracle match, v0 < 0, finite kcat and |vN| > 1e-3.
- For uni-OnlyA: D_A == D_I == 1 and the rescale to kcat 5.0.

The redeclared mechanisms (557-563 ≡ 322-330, 574-580 ≡ 515-523, and three copies of uni-OnlyA) disappear. Compile is already shared with the gates.

_Merged from:_ T#106; T#108

#### loc/AG4 — allosteric_ground_truth.jl: remove guard assertions the preceding constructor already makes; fold one duplicated B point
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 5; risk low
- Sign-off: approved by Denis

**Testsets.** test/allosteric_ground_truth.jl:922-925, 1000-1004, 1198-1202

**Change.** Delete `@test ER._onlya_haldane_violation(…) === nothing` at 1003-1004 and 1201-1202; keep a comment that construction proves acceptance. Put B=1.0 into the loop tuple at 923 and delete the separate assertion at 922.

**Coverage.** `ER.AllostericMechanism(err1/dead)` on the line before lifts into the inner constructor. That constructor calls _onlya_haldane_violation on the same canonical fields and errors on any violation (src/types.jl:922-924), so the deleted assertion cannot fail on its own. The B=1.0 assertion is identical to a loop entry.

_Merged from:_ T#110

#### loc/CO1 — test_allosteric_collapse.jl: one uni builder with an RE mask, one random-order bi-bi builder, evalrate for the bi-bi equilibrium
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 30; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_allosteric_collapse.jl:9-34, 71-96, 130-169

**Change.** (1) Merge uni and uni_ss (9-22) into `uni(states; re=(true,false,true))`, and build 'mixed RE/SS coupled' (134-139) with it.
(2) Hoist the bi-bi construction, duplicated at 74-85 and 149-160, into `ro_bibi(states)`.
(3) In both bi-bi testsets, replace the hand-rolled params, concs and equilibrium code with `fp, v, veq = evalrate(am)` and `@test abs(veq) < 1e-8`. Keep the UndefVarError comment.

**Coverage.** The mechanisms are identical: same Steps, flags and tags. Every assertion stays: the collapsed forbidden split, the explicit mirror string, a finite rate and zero at equilibrium. evalrate's point (P=3, others 1, Keq=3) also satisfies P·Q/(A·B) = Keq for the bi-bi, and it adds an off-equilibrium finiteness check. Only the arbitrary draw changes.

_Merged from:_ T#105

#### loc/GD1 — Golden reference: drop the PARAMS_REDUCED lines and assert the parameters wrapper exactly in the golden loop instead
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 13; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_allosteric_golden.jl:9-31; test/reference/allosteric_golden_reference.txt (12 lines)

**Change.** Remove the PARAMS_REDUCED push (test_allosteric_golden.jl:18-19) and its 12 reference lines. Inside the golden testset add `@test parameters(m, EnzymeRates.Reduced) == (fitted_params(m)..., :Keq, :E_total)` per allosteric spec. Reword defs:2479-2480 ('Their golden PARAMS_REDUCED confirm…'), which DR1 also touches.

**Coverage.** A script confirmed that all 12 PARAMS_REDUCED tuples equal the `(; …) = params` header of the same spec's REDUCED_STRING, which is pinned byte for byte, so indep is pinned exactly. The added line keeps exactness for the exported parameters(m, Reduced) wrapper itself (the verifier's condition). test_accessors checks membership only.

_Merged from:_ T#109

#### loc/FI1 — test_fitting.jl: merge the two BBO kcat-rescaling testsets and cut their fits to maxtime 0.5 s
- Lever: lower_threshold; time saved ≈ 25 s; test lines removed ≈ 37; risk low
- Sign-off: approved by Denis
- Measure first: One run of test_fitting.jl: confirm no convergence warning leaks and kcat ≈ target on all three fits.

**Testsets.** test/test_fitting.jl:349-425

**Change.** Replace 'scale_k_to_kcat normalization' (350-391) and 'fit_rate_equation kcat rescaling' (394-424) with one testset. It builds the 8-point data once and runs 3 BBO fits at maxtime=0.5: default target 1.0 with n_restarts=2 (keeps the best-of-restarts loop), scale_k_to_kcat=7.0 with n_restarts=1, and `nothing` with n_restarts=1. Assert keys == fitted_params, isfinite(loss), retcode isa Symbol and `_kcat_forward ≈ target rtol=0.01`, plus the `nothing` keys check.

**Coverage.** No assertion depends on fit quality. rescale_parameter_values (src/rate_eq_derivation.jl:1000-1010) scales every SS constant by target/kcat_current, and uni_uni's kcat is degree-1 homogeneous in those constants, so kcat ≈ target holds at any finite fitted point.

The 42.0 and 7.0 targets take the same path (fp.scale_k_to_kcat). BBO always takes at least one step before checking a stop condition, and MaxTime maps to ReturnCode.MaxTime with no warning, so fit_capturing_convergence stays pristine. Fit quality stays covered by identify's recovery and the CMA tests.

_Merged from:_ T#138; T#159; T#167; src_enum#60; src_identify_fit#31

#### loc/FI2 — test_fitting.jl: delete 'Loss at true params is zero' and 'Centering invariance'
- Lever: delete_redundant; time saved ≈ 0 s; test lines removed ≈ 52; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_fitting.jl:84-136; subsumers 139-170, 173-200

**Change.** Delete test_fitting.jl:84-104 and 106-136.

**Coverage.** - 103 makes the same assertion on the same true_params, 5-point concs, default FittingProblem and x_true as 'Absolute mode uncentered loss' (185, 189).
- Single-group centering invariance is the one-iteration case of the per-group loop in loss! (src/fitting.jl:159-170), which 'Multi-group centering invariance' (139-170) tests with unequal groups and scales at 10 random x.
- The single-group offset at x_true stays at 195-198.

_Merged from:_ T#147

#### loc/FI3 — test_fitting.jl: hoist the shared uni-uni true_params and 5-point concentrations
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 12; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_fitting.jl:63-200, 271-290, 427-502

**Change.** Define `Keq_val = 2.0`, `true_params` and `concs5` once after uni_uni (line 18), and use them in the testsets that repeat them.

**Coverage.** Code moves only; every assertion keeps its values. The count already excludes the copies that FI1 and FI2 delete.

_Merged from:_ T#148

#### loc/ID1 — test_identify_rate_equation.jl: delete the uni-uni popsize and loss_parsimony identify runs; pass the parsimony kwarg to the all-cap-skipped run
- Lever: delete_redundant; time saved ≈ 9.5 s; test lines removed ≈ 42; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_identify_rate_equation.jl:1104-1172

**Change.** Delete 'identify runs on a solver that rejects popsize' (1104-1125) and 'loss_parsimony_threshold threads through identify_rate_equation' (1152-1172). Add `loss_parsimony_threshold=2.0` to the identify call in 'all-cap-skipped expansion batch is reported (M2)' (1139-1143), with a one-line comment that an explicit non-default value proves the keyword is accepted.

**Coverage.** - Popsize contract (no injected solver option): every CMA identify run with default solver_kwargs covers it. That includes the main run (245-253, asserting `isa IdentifyRateEquationResults`) and the M2 run. An injected option would make every base fit throw and the run raise.
- Parsimony kwarg: the M2 run proves it is accepted. Its assertions are unchanged, because the base count has no smaller count, so the cutoff is `nothing`.
- Both deletions are needed: with only one gone, its 9.8 s of per-type compile moves to the other. No later testset uses the 4-6-parameter uni-uni children.

_Merged from:_ T#140; T#169; src_identify_fit#33; src_enum#60

#### loc/ID2 — test_identify_rate_equation.jl: delete 'a base-tier row of a degenerate seed's child carries no parent'
- Lever: delete_redundant; time saved ≈ 2.5 s; test lines removed ≈ 30; risk medium
- Sign-off: approved by Denis

**Testsets.** test/test_identify_rate_equation.jl:1746-1775

**Change.** Delete test_identify_rate_equation.jl:1746-1775.

**Coverage.** - `child in base && !(child in seeds)`: follows from 1677-1679. `unique!(base)` plus |base| = 278 = 257 + 21 means the cures are disjoint from the seeds.
- Parentless base rows: the test calls _process_batch directly with the default `parent_of = Dict()`, so it never tested base-tier wiring. Real base-tier rows are asserted parentless in the main run (365-367).
- Row fields: mechanism_type at 310, n_params at 935, and the loss copy uses the same row code.

The deleted cost (re-running init_mechanisms + _base_tier and compiling one new ping-pong child) is unique to this test.

_Merged from:_ T#142

#### loc/ID3 — test_identify_rate_equation.jl: assert the original-mechanism property inside fit-dedup; delete the LDH-split half
- Lever: cheaper_same_coverage; time saved ≈ 1.7 s; test lines removed ≈ 14; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_identify_rate_equation.jl:1290-1352, 1532-1579

**Change.** Add `@test [f.mech for f in bad_failures] == [m1, m2]` after line 1351 in 'fit-dedup by eq_hash in _process_batch'. Delete the LDH-split half (1564-1578) of '_process_batch failures report the ORIGINAL mechanism', keeping the m_bad compile-stage half.

**Coverage.** A fit-stage FitFailure is built from c.orig (src/identify_rate_equation.jl:599-603), whether the mechanism is the representative or an inherited duplicate. pmap preserves order, so the added check covers both representative and duplicate, a superset of the deleted single check. The LDH structure is irrelevant to this code, and its FittingProblem/Optimization specializations are unique to the deleted half.

_Merged from:_ T#144

#### loc/ID4 — test_identify_rate_equation.jl 'mechanism recovery': assert the beam's recorded loss for the selected row instead of refitting
- Lever: cheaper_same_coverage; time saved ≈ 1 s; test lines removed ≈ 6; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_identify_rate_equation.jl:255-272

**Change.** Replace the refit at 255-268 (n_restarts=3, maxtime=10) with `best_row = results.cv_results[EnzymeRates._select_best_row(results.cv_results), :]; @test best_row.loss < 0.01`, sharing best_row with 'results structure'.

**Coverage.** cv_results.loss is the beam's training loss for that exact row, on the same full data, through the same FittingProblem/fit_rate_equation path (src/identify_rate_equation.jl:586-588). Baseline losses are about 1e-13, 11 orders of magnitude under the threshold.

A public FittingProblem on an AllostericEnzymeMechanism singleton is still exercised by every LOOCV fold. The assertion now rests on a single-restart fit with maxtime 1, which is slightly weaker. The change also removes a 30 s worst-case tail.

_Merged from:_ T#145; src_enum#60

#### loc/ID5 — test_identify_rate_equation.jl: merge the three _cv_fold_loss testsets; delete the 'one fold' duplicate
- Lever: delete_redundant; time saved ≈ 0.5 s; test lines removed ≈ 35; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_identify_rate_equation.jl:422-496

**Change.** Make one '_cv_fold_loss' testset that builds rxn, m, data and prob once. Keep the all-groups success block (422-452) and the loud-failure @test_throws (454-475). Delete 'one fold, finite' (477-496).

**Coverage.** 'over all groups' makes the identical call for g = 2 with identical kwargs. Its Vector{Float64}, ≥ 0 and isfinite checks on scores[2] imply the deleted Float64, ≥ 0 and finite checks. One duplicate CMA fold fit goes.

_Merged from:_ T#146; T#170; src_enum#60

#### loc/ID6 — test_identify_rate_equation.jl LOOCV tests: exotic labels through the real flatten, one twin fixture for dedup and flatten, _offer_cv! checks in the _ingest! testset
- Lever: delete_redundant; time saved ≈ 0.2 s; test lines removed ≈ 45; risk low
- Sign-off: approved by Denis
- Measure first: One focused run of the merged testset.

**Testsets.** test/test_identify_rate_equation.jl:667-691, 1077-1102, 1807-1888

**Change.** - Delete 'cv_results: exotic group labels survive CSV roundtrip' (667-691), which calls no package code and re-implements the src flatten loop.
- Merge 'LOOCV eq_hash-uniqueness guard' (1807-1852) and '_cv_model_selection flatten reproduces serial LOOCV' (1854-1888) into one testset: the twins [m1, m2], 3-group data with exotic labels ("a=b", "c,d", "x y"), two rows (0.5, 0.2) and one _cv_model_selection call. Assert nrow == 1, eq_hash == h, loss == 0.2, flat == serial (serial on the kept twin), and `all("cv_fold_$g" in names(CSV.read(joinpath(dir, "loocv_results.csv"), DataFrame)) for g in groups)`.
- Move the _offer_cv! repeat-hash block (1838-1851) into '_ingest! and cv pool', reusing its `mk` builder.

**Coverage.** Every assertion is kept, and exotic-label coverage goes up: the real flatten (src/identify_rate_equation.jl:1072-1075) and the real loocv_results.csv write now run on exotic labels instead of a copy of the logic.

The group count (2 against 3) does not affect dedup or flatten. The twins render the same equation, and the stub sets every parameter to the same value, so the fold losses do not depend on which twin is kept. The _offer_cv! and _ingest! checks are complementary and both stay.

_Merged from:_ T#149; T#154

#### loc/ID7 — test_identify_rate_equation.jl: replace the six serialized-Sig string constants with DSL-built fixtures
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 55; risk low
- Sign-off: not required (coverage unchanged)
- Measure first: While converting, check once (not committed) that each DSL mechanism == its Sig reconstruction.

**Testsets.** test/test_identify_rate_equation.jl:1254-1530, 1566, 1811, 1857

**Change.** Rewrite _DEDUP_SIG1/2 (1254-1288), _CANON_SIG_MERGED/SPLIT (1355-1419) and _ALLO_SIG_SPLIT/MERGED (1443-1510) as @enzyme_mechanism / @allosteric_mechanism fixtures:
- Use parenthesized step tuples for shared groups and `:: EqualAI` tags with oligomeric_state 4 for the allosteric pair.
- Use the attach(em) pattern (1694) where an atom-bearing reaction is needed.
- Delete the recon/recon_am closures (1291, 1422, 1513, 1566, 1811, 1857).

**Coverage.** These are the same mechanisms. Each testset already re-checks its precondition at run time (1299-1301 collision, 1429-1431 distinct keys, 1519-1528), so a mismatch fails loudly. The DSL fixtures also survive future canonical-order changes that would invalidate serialized Sigs.

_Merged from:_ T#150

#### loc/ID8 — test_identify_rate_equation.jl: one shared uni-uni identify fixture (reaction, 4-row data, problem) (style call)
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 22; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_identify_rate_equation.jl:579-587, 886-887, 915-925, 1079-1080, 1131-1137, 1217-1225, 1838-1839

**Change.** Add file-level `_testhelper_uni_rxn` (S[C] → P[C]) and `_testhelper_uni_prob(data_form)`, keeping both the NamedTuple/String-group and the DataFrame/Int-group forms. Use them at the sites that survive ID1 and ID9 (579-587, 915-925, 1131-1137, 1217-1225), and use `first(init_mechanisms(_testhelper_uni_rxn))` in the three BatchEntry builders.

**Coverage.** Same reaction, data and assertions; code only moves. Both data forms stay exercised. Denis's inline-fixture rule is written for enumeration tests only, but this is still his style call.

_Merged from:_ T#151

#### loc/ID9 — test_identify_rate_equation.jl: fold the fitted-set testset into _process_batch, reusing its fixture
- Lever: test_code_simplification; time saved ≈ 0.3 s; test lines removed ≈ 12; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_identify_rate_equation.jl:914-982

**Change.** Move the two fitted-set `_process_batch([m], …; fitted)` calls and their assertions (957-982) into the '_process_batch' testset (914-955) with `m = first(ms)`, and delete the duplicated rxn/data/prob block (958-969). Optionally assert solve counts with _CountingStubOpt (its definition must then move above 914).

**Coverage.** The fitted-set skip (src/identify_rate_equation.jl:518-529) runs before any optimizer and is asserted unchanged (ss1 == 0 with one entry; ss2 == 1 with nothing emitted). Real CMA fitting of the same uni-uni type stays asserted in _process_batch.

_Merged from:_ T#155

#### loc/ID10 — test_identify_rate_equation.jl: merge the _rows_to_dataframe pair and the _select_beam trio into single testsets
- Lever: test_code_simplification; time saved ≈ 0 s; test lines removed ≈ 8; risk low
- Sign-off: not required (coverage unchanged)

**Testsets.** test/test_identify_rate_equation.jl:147-200, 523-569, 766-808

**Change.** Merge '_rows_to_dataframe' (147-176) with '_rows_to_dataframe with failure row' (178-200), and 'beam selection' (523-569) with '_select_beam best_override' (766-777) and '_select_beam parsimony_cutoff' (779-808). Keep the lean search kwargs in 'save_dir non-empty check': dropping them would turn a regression of the up-front check into an hours-long run.

**Coverage.** Only testset scaffolding goes; every assertion is kept.

_Merged from:_ T#152

#### loc/ID11 — test_identify_rate_equation.jl dedup-key partition: drop the vacuous uni-uni entry and check render determinism on a sample
- Lever: lower_threshold; time saved ≈ 4.5 s; test lines removed ≈ 8; risk low
- Sign-off: approved by Denis
- Measure first: In a warm session, record the first and second render times per bi-bi mechanism. The saving ≈ Σ(second render) × 0.96. It is estimated at 3.5-6 s, because the first render is mostly per-type specialization.

**Testsets.** test/test_identify_rate_equation.jl:707-764

**Change.** Turn the loop over test_reactions into one bi-bi block, since uni-uni has one init mechanism, so one class. Render each of the 239 bi-bi mechanisms once for bucketing. Re-render for the `=== h` determinism check only on a fixed sample (`i % 25 == 1`, 10 mechanisms). Keep the 239-class pin.

**Coverage.** The 239-distinct-classes pin and derivability of every bi-bi equation stay fully covered, because every mechanism is still rendered once. That single render is also what test_mechanism_enumeration.jl 'bi-bi exit gate' (T#168) relies on.

The uni-uni entry asserts 1 class from 1 mechanism and can fail only by crashing; uni-uni rendering and the eq_hash format stay covered at 936.

Only the in-process determinism re-render is sampled. A defect affecting a fraction f of mechanisms is caught with probability 1-(1-f)^10 instead of 1-(1-f)^239; src has no global caches, RNG or objectid hashing (grep).

_Merged from:_ T#139; T#165

#### loc/SRC1 — With the src deletion of the test-only Parameter helpers (src_types#1/#87), delete the tests of that dead code
- Lever: src_change; time saved ≈ 0 s; test lines removed ≈ 30; risk low
- Sign-off: approved by Denis

**Testsets.** test/test_types.jl:914-955, 1584-1612, 1770-1786

**Change.** If Denis approves src_types#1/#87 (delete _flip_to_inactive, _force_inactive, _param_for_symbol, its walkers, governing_step/is_i_state and the Keq/Etot singletons):
- Delete test_types.jl '_force_inactive forces :I regardless of tag' (1770-1786).
- Delete the governing_step/is_i_state asserts (925-927, 935-937, 944-945) and the Keq()/Etot() singleton asserts (949-950, 952-953).
- Delete the Keq()/Etot() name asserts at 1601-1602 and 1609-1610; 1696-1697 go with TY5.
- Keep the Lallo lines, since src constructs Lallo.
- DR3 already removes _dep_struct_key, the only other caller of _param_for_symbol.

**Coverage.** These tests exercise only code that the src proposal deletes. A repo-wide grep finds no src caller of these functions or of Keq()/Etot() construction, so no package behaviour loses coverage. The :Keq and :E_total names users see stay pinned through parameters() in the specs and the golden file.

_Merged from:_ src_types#1; src_types#87; src_derivation#0; src_derivation#93

## Cross-file critic

### Ordering

These are the cross-group constraints only; the order inside each group follows its own depends_on.
1. Pure deletions and compactions, in any order, each followed by a full suite run: types T1, T2, T15, T16, T18; thermo T1, T3, T5; deriv V1-V5, V10; dsl D5, D17; enum T1, A1→A2, A7, A8, B1, C1, SP; identify V1, L1, FR1, C1, R1, B1, B2, B4, B5, CV1, CV3, M1.
2. types foundations, one commit each: T5 → T13 → T14 (Canonical Step Form gate) → T6 → T12 → T17 → T4 (fix the helper's name here) → T3 → T7 → T8 → T9 → T10 → T11.
3. Enumeration items that consume the types foundations:
   - after T17: R1 → R2 → F, A10;
   - after T12: A3, A4 → A5, A6 (S0 first);
   - S0 → S1 → S2; S3 → S4.
4. thermo T4 and T6 (after types T4); thermo T1 → T2; thermo T7 before deriv V14-V16.
5. Derivation:
   - V6/V7/V8 only after types T17 and enum R1+R2;
   - V11 together with enum RC, then enum D1 → D2;
   - V12 (lift via `_concrete`) → V13;
   - V14 → V15 → V16, written without all_params and without the kernel; thermo T8 lands here;
   - V17 → V18, together with thermo T9 (V18 owns the body);
   - V19 → V20 → V21 → V22 → V23 (absorbs types T10's line-784 swap).
6. dsl D1 → D2 → D6 → D7 → D8-D13 → D14 → D15 → D16 → D18 → D19.
7. identify:
   - V1 → V2, V3, and V4 (V4 after types T11);
   - C1, R1, B1, B4 → B6 → B7;
   - B2 → B3, B8.
8. Decisions last, each in its own commit with re-pins:
   - types T19 + dsl D4 (after T6 and D1);
   - types T20;
   - types T21 + enum A9, as one change;
   - thermo T11 or deriv V26 (after thermo T2);
   - thermo T10, deriv V24, V25 (after V18), V27 (before the allosteric half of V12/V13 if adopted);
   - enum D2, S5, A5;
   - types T22, T23;
   - identify FR2, FR3, FP1, V2, V3, B7, B8, C2.

### Conflicts and resolutions

- **types/T19 + dsl/D4 (+ enum/D1)** — Both items make the same change: remove RegulatorMults.allowed_multiplicities and the `A(1, 2)` / `I(2)` syntax. Their types.jl estimates disagree: T19 says 7 lines, counted after T6; D4 says about 12, counted before T6. D4's depends_on omits types/T6, even though T6 rewrites `_to_sig(::RegulatorMults)` and `_regulator_mults_from_sig`, which D4 edits. Both items also edit mechanism_enumeration.jl:2398 (`Int[1]`), but that line is inside `_add_competitive_inhibitor`, which enum/D1 deletes. **Resolution:** Treat this as one decision for Denis and land it as one commit, after types/T6 and dsl/D1. Count 7 types.jl lines (T19) and 10 dsl.jl lines (D4). If enum/D1 lands first, the edit at enumeration:2398 disappears.

- **types/T21 + enum/A9 (+ dsl/D16)** — Both items make the same change: the RegulatorySite constructor sorts its ligands. T21's 'other' field also claims the enumeration -3 (the sortperm at 3008-3012) that A9 counts. D16 justifies keeping source ligand order with 'RegulatorySite does not sort ligands' and with test helpers returning source-order sites. That reason becomes false. The oracle itself (`_positional_site_idx`, test_rate_eq_derivation.jl:262) matches sites by ligand-name Set, so it keeps working. **Resolution:** Make it one item needing Denis's decision. Count it once: types about +1, enumeration -3. Reword D16's comment.

- **deriv/V26 + thermo/T11** — Both items make the same decision, sorting the allosteric independent parameters by name, but implement it in two places. V26 adds a sort in `_dependent_param_exprs(am)`; T11 moves the sort into `_solve_dependent_set`'s return. The net counts are -1 and +1. **Resolution:** Adopt one. T11 is cleaner: it lands after thermo/T2 because both edit the same return line, and it deletes the plain path's separate sort. Regenerate the allosteric golden file once.

- **deriv/V11 + enum/RC** — Both items introduce the same `_inactive_groups(am)` helper. V11 counts the rate_eq_derivation.jl side (+8/-10) and RC counts the mechanism_enumeration.jl side. **Resolution:** V11 defines the helper and RC consumes it, in one commit. Their nets do not overlap. enum/D2 must land after RC, because D2 calls `_redundant_copy_groups`.

- **types/T17 + enum/R1/R2 + deriv/V6/V7/V8** — V6 and V7 read `_indexed_re_segments`, which R1 deletes. V6 deletes `_compute_re_groups`, but enumeration:1898 (`_re_segment_count`) and :2068 (`_group_re_segments`) still call it until R2 and R1 land. R1's 'other' field claims the same types.jl +3 that T17 counts. **Resolution:** Land T17 first, with `_re_segment_extras` returning a 5-tuple, then R1 and R2. After that, V6, V7 and V8 destructure `species, segments, _, idx, seg = _re_segment_extras(steps(mech))`, and no `_compute_re_groups` shim is needed. Count the types.jl +3 once, in T17.

- **thermo/T7 vs deriv/V14, V15, V16** — T7 drops `_assemble_constraints`' `all_params` keyword and moves `_state_all_params` into thermo as `_param_columns`, claiming -15 in rate_eq_derivation.jl. V14 says to keep all_params explicit at the allosteric call sites, V15 rewrites `_state_all_params`, and V16 passes `all_params = cols`. All three assume the keyword survives. T7's -15 also falls on lines (864, 1496-1497, 1555, 1225, 1326) that V15, V16, V17, V19 and V20 delete or rewrite. **Resolution:** Land T7 before V14-V16, and write those three items without all_params. Give the column helper one name in one file. Do not count T7's claimed -15 in rate_eq_derivation.jl.

- **thermo/T8 vs deriv/V14, V18** — T8 deletes `_dependent_param_exprs_kernel` and its Type wrapper. V14's rename builder still makes 'one kernel call', and V18's non-allosteric Reduced string uses `first(_dependent_param_exprs_kernel(Mechanism(m), Dict()))`. T8's call-site edits (rate_eq:129, 673, 1223) land in code that V14 and V18 rewrite. **Resolution:** Write V14 and V18 against the inlined `_solve_dependent_set(_assemble_constraints(mech, rename; step_params)...)`, and land T8 together with them or after them.

- **thermo/T9 vs deriv/V17/V18** — T9 inlines both `_build_rate_body` methods into rate_eq_derivation.jl (+~12 lines) and requires the generated Expr to stay byte-identical with no added `return`. V18 replaces those bodies with one `_rate_body` that does add `return`. V17 and T9 each propose the same `_dep_assignments` helper. **Resolution:** V18 owns the body builder. T9 shrinks to deleting `_sorted_raw_param_symbols` and `_build_rate_body` and moving `_destructuring_expr`. The explicit `return` is harmless: the allosteric body already has one, and the Expr-shape test (test_rate_eq_derivation.jl:1666) walks only the num/den Expr. Define `_dep_assignments` once, in V17. T9's +4 line move and T8's +3 docstring move are not netted anywhere, so subtract 7.

- **types/T3 vs deriv/V12** — V12's merged `parameters(::AbstractEnzymeMechanism, Full)` sketch calls `name(p, m)`. After T3, `name(p, ::EnzymeMechanism)` and `name(p, ::AllostericEnzymeMechanism)` raise a MethodError. Today every src caller passes a concrete mechanism (verified). **Resolution:** In the generator, V12 must lift first: `cm = _concrete(M())`, then `name(p, cm)`. Define `_concrete` once.

- **types/T4 (optional rename) vs thermo/T4, thermo/T6, deriv/V13, deriv/V15** — Five items call `_emit_cat_params_for_rep` on arbitrary steps, not only on group representatives. T4 optionally renames the function to `_step_params` (and thermo/T4 suggests `_step_constants`). **Resolution:** Decide the name once, in the types/T4 commit, before the thermo and derivation items use it.

- **types/T12 vs enum/R2, SP, A3, A4, merge move** — T12 renames call sites at enumeration lines 1293, 2143, 2148, 2905, 2920, 2932, 3020 and 2816-2818. R2, SP, A4, A3 and the merge move rewrite those same lines. T12's -2 at 2816-2818 lies inside A3's rewrite of `_make_am_with_added_reg`. **Resolution:** Land T12 first. The enumeration items then call `_with(am; groups=…/states=…/sites=…)`. The -2 belongs to A3; neither group currently counts it twice.

- **types/T10 vs deriv/V23** — Both items edit rate_eq_derivation.jl:784 (`products(M())`): T10 swaps the line in place, and V23 rewrites the whole non-allosteric kcat block around it. **Resolution:** Sequence the two. Whichever lands second uses `name.(products(reaction(Mechanism(M()))))`.

- **identify/B9 vs enum/S2** — B9 applies only if the enumeration group drops the per-child `_assert_atom_conserving` checks. S2 explicitly keeps them: they are the only error-injection route for the failure-isolation tests, and the docs rely on them. **Resolution:** Drop B9, which removes 1 line from the identify total.

- **types/T22 vs deriv/V12, V27, thermo/T6** — T22 collapses Kd/Kiso, Kon/Kfor and Koff/Krev into three types. V12 and V27 filter with `p isa Union{Kon, Koff, Kfor, Krev}`. The chokepoint recognizer regex at test_types.jl:2417 lists the six type names. T22 needs about 100 test edits to save 8 src lines. **Resolution:** Recommend dropping T22. If Denis keeps it, land it after V12 and V27 and update the AST-walker regex.

- **deriv/V27 vs V12/V13** — V27 deletes the allosteric Full enumeration that V12 and V13 have just restructured, and it re-adds an allosteric `_ss_rate_constant_names` method that V12 had merged away. **Resolution:** Get Denis's decision on V27 before doing the allosteric half of V12 and V13.

- **types/T9 + types/T10** — Both items claim the AllostericEnzymeMechanism forwarder loop and forwards at types.jl:1684-1693. **Resolution:** About 2 lines are double counted. Subtract them.

- **dsl/D6, D19 vs test/mechanism_definitions_for_test_enzyme_derivation.jl:65-75** — The `@enzyme_mechanism_src` and `@allosteric_mechanism_src` test macros call `_reject_allosteric_syntax!` (deleted by D6) and `_parse_plain_mechanism_body` / `_parse_allosteric_mechanism_body` (merged by D19). **Resolution:** Edit the test macros in the same commit (test side only).

- **ID collisions across groups** — T1-T23 (types), T1-T11 (thermo), T1 (enum) and T1-T3 (identify) share IDs, as do D1/D2 (dsl) and D1/D2 (enum). **Resolution:** Use group-prefixed IDs (types/T1, thermo/T1, …) in the plan.

### Understated behaviour changes

- **types/T8 combined with types/T21 or enum/A9** — `show(aem)` now lifts through `AllostericMechanism(aem)`. Once the RegulatorySite constructor sorts ligands, a site written out of name order prints its ligands, and possibly its sites, in a different order. T8's 'identical for every pinned case' holds only if T21/A9 are rejected.

- **types/T21 + enum/A9** — The effect goes beyond Kreg order inside a site. `_regulatory_site_canonical_key` sorts sites by their ligand-name tuple, so sorting ligands can also reorder the sites. That changes the order of reg Ks in fitted_params, the params destructure line, rate_equation_string, eq_hash, and the RS type parameter of compiled allosteric types. Old CSV mechanism_type strings then no longer equal freshly compiled ones. This needs Denis's decision.

- **enum/S4 (also deriv/V6)** — Form identity switches from name-keyed to Species-keyed. `name(::Species)` concatenates labels with no separator, so E with {NAD, H} bound and E with {NADH} bound both render :ENADH. Today such forms are silently merged. Afterwards both forms are kept, and the constructor (`_assert_each_reaction_once`) will likely error. This is an edge case, but it is not 'none'.

- **enum/D2** — Today the allosteric path tests redundancy only on the new group, and only when `_all_twin` holds. D2 tests every copy group on the built canonical child. The verifier checked only the non-allosteric bi-bi depth-2 population, so allosteric and ter-ter child counts can drop. It is also slower, because every candidate child is now constructed before the test. This should be a needs_denis_decision item with count re-pins and a timing gate.

- **deriv/V27** — Beyond the MethodError on the exported `parameters(aem, Full)`, which the parameters docstring documents today, V27 deletes the allosteric PARAMS_FULL golden lines and two testsets. That is a coverage reduction, which needs explicit sign-off.

- **test deletions inside 'behaviour none' items (deriv/V7, enum/F, enum/D1, enum/A1, identify/CV2, identify/CV3)** — Each item deletes a testset or assertions:
- V7: the `_bareiss_det` testset;
- F: the 110-line screen-agreement testset;
- enum D1: the `_add_competitive_inhibitor` and `exclude_regs` testsets;
- A1: the no-op assertions;
- CV2 and CV3: halves of testsets.
Under Denis's rule every one of these needs his sign-off. A7, R2 and types T1 already say so; these items do not.

- **identify/B8** — Inf-loss fits leave the LOOCV pool. At sparse parameter counts this can change which candidates are cross-validated, and therefore which model is selected. The item describes this as a robustness fix, but it can change the outcome.

- **enum/A6** — The Haldane prefilter moves from the `rxn` argument to `reaction(m)`. The two are equal on identify paths. A direct caller (a test) that passes a different rxn, such as a reaction with extra declared regulators, gets different filtering.

- **deriv/V18** — The non-allosteric generated body gains an explicit `return`, so it is no longer byte-identical, as thermo/T9 demands. Values and performance are unchanged. Only the claim of byte identity is understated.

### Invariant risks

- **thermo/T2 (Haldane pivot rule)** — This replaces the solver that picks dependent constants. I verified the equivalence argument. In today's Gauss-Jordan, once row r takes its pivot p_r, every nonzero at a non-pivot column of row r stays ≤ p_r in the (priority, then lowest index) order: later eliminations add only entries ≤ p_k, and p_k ≤ p_r. So today's basis is the unique max-weight basis, which is what column-order RREF finds on stably sorted columns. The remaining risk is in the details: tuple priorities with rev=true, and how the rhs pivot is detected. Gate on the fitted_params pins, test_pivot_priority_regression.jl and the allosteric golden file.

- **types/T13, T14, T16 (Canonical Step Form)** — These rewrite the canonicalizing constructors. The equivalence is plausible: the keys are distinct within a group (each reaction appears once), sorts are stable on Julia ≥1.9, and T14's table is built from the same original steps that `_entry_kind` scans today. They still need their own commits, each gated by the canonicalization, fitted_params and count pins.

- **types/T21 + enum/A9** — They add a new canonicalization inside a constructor, which changes ==, hash and the Sig of allosteric mechanisms. Land them only on Denis's decision.

- **types/T5** — I verified that every generated hash nests fields in today's order for all 15 types, with Species kept explicit, so the hash values are identical. The @eval loop must be defined before any top-level code builds a Dict or Set over these types.

- **types/T6 (+T19)** — T6 produces an identical Sig (verified): the second Mechanism that `_drop_unbound_regulators` builds re-canonicalizes steps in which no unbound regulator appears. T19 changes the Sig of every mechanism that binds a regulator.

- **types/T3, T22 (name chokepoint)** — T3 is safe, because every src caller of name(p, ·) passes a concrete mechanism (verified). V12 is the one new caller that must lift first. T22 renames the Parameter subtypes that the chokepoint recognizer regex lists (test_types.jl:2417).

- **types/T12** — T12 drops the defensive copies. The constructors copy steps and states (Step, Species and Residual sort copies, and `permute!(copy)` for states), but RegulatorySite objects, which hold Vectors, stay shared between parent and child. That is safe only while no move mutates a site in place; today none does.

- **deriv/V9 (Pass-2 absorption)** — Applying the rename at the RE leaves gives the same polynomials: substitution is a ring homomorphism, and the only branch on intermediate polynomials, `_compute_alpha`'s cycle check, reads concentrations only. But this is the load-bearing absorption, so gate on the exact strings for the non-competitive-inhibitor and non-essential-activator specs and on test_mechanism_enumeration.jl:9255 and 9336.

- **deriv/V14, V16 + thermo/T7, T8 (Pass-2 builder)** — The chain-compression loop is provably dead: a bare-Symbol dependent's right-hand side is a non-pivot column, so it is never a rename key. Column equality under T7 rests on `_assert_uniform_groups` and `_assert_each_reaction_once`. Four items edit the same call graph, so the order in which they land matters.

- **deriv/V25, thermo/T10, thermo/T11 or V26** — These deliberately change which constant becomes dependent, which ties are annotated, or the parameter order. Their outputs and eq_hash change, which is why each is a decision.

- **deriv/V17, V18, V22, V2, V4 (rate_equation 0-alloc/<120 ns)** — These rebuild the generated body. They must keep the `_nest_binary` 2-operand trees, `build_power_expr`'s varargs `*`, and every factor order. test_rate_equation_performance is the gate. V17 reorders the allosteric assignments, which is safe because full Gauss-Jordan leaves no dependent referencing another.

- **enum/S1, S3, S4 (seed identity)** — The equality scripts need Julia. S3's 'every isomerization under its own key' depends on S4's precondition that no isomerization has a mirror step. Today a mirror inherits its catalytic step's group id (mechanism_enumeration.jl:1124-1126), and an isomerization mirror would break the equivalence.

- **dsl/D7** — D7 removes the parser-side sorting and relies on the constructors to sort. I verified that the Species, Residual and Step constructors sort into copies. Every DSL fixture passes through this code, so run the full suite.
