# Derivation bug history and its causes

Sources: the main history (120 commits), the squash messages of #14–#83, sub-commits on the `package-simplification` and `mwc-solve-then-limit` branches, `.superpowers/harness/record/*`, `.superpowers/sdd/2026-10-08-merge-decisions/items.json` and `ROADMAP.md`. Line numbers refer to main at 277a785. I did not run Julia. Anything marked *(inference)* comes from reading code, not from a commit or a measurement.

## Headline

- **Allosteric axis: the record supports the premise.** The largest allosteric bug families come from three things. The inactive state was derived on a separate, pruned, re-canonicalized graph. `:OnlyA` was applied by deleting steps before the constraints were solved. And the two states were joined by rendered name. The incidents: §5a's 54,674 errors (#56), #70's 1,354 errors, and the open 766 + 768 cases (ROADMAP:10-33).
- **Kinetic-group axis: the record contradicts a literal reading of the premise.** The first implementation solved thermodynamics without the user equalities and got inconsistent laws. dec0707 (#16, Feb 2026): "the Haldane/Wegscheider thermodynamic constraint matrix was ignoring user-defined parameter constraints (e.g., K3=4K1 …), causing constrained rate equations to differ". The fix merged them into one matrix. The bugs that came later come from how they are joined, not from solving them together:
  - by rendered name;
  - with a column priority that depends on which group a step is in;
  - with a second rename before the solve (Pass 2) that leaves zero columns.
- **So the workable form is:** one linear system over bare step constants (cycle rows + group-equality rows + A/I tie rows), solved once with a fixed priority. Then fold the pure equalities, then normalize, then take the limits. "Solve, then substitute" only works if every substitution triggers a re-solve, which makes it a joint solve anyway.

## 1. Bug catalog

| Class | Examples (commits) | Root cause (file:line) | Pinned by | Open? | Verdict under the new design |
|---|---|---|---|---|---|
| **A. Kinetic groups join the solve by rendered name; the pivot depends on the group** | dec0707 (#16: matrix ignored `K3=4K1`); a06a051 (constraint K symbols were position-based, wrong after a sort); 885ace1 (#14: groups create `0=0` / `0=c·logKeq` rows); ee462a1 (group pins K=1 → MethodError); a39cc93 (re-canonicalized order changed the Haldane dependents → dedup keys diverged) | thermo:376-389: each step's constants go into `sym_col[get(rename, name(p,mech))]`, and a column takes the priority of the *last* step written into it (documented at thermo:84-94). types.jl:1386-1403: name = the group representative's sides | test_rate_eq_derivation.jl:1531-1575, 1675-1700, 765-792 (per-spec constraint counts) | Crashes closed. Group-dependent pivot still live *(inference)* | **Eliminated** only with one joint solve and a membership-independent priority. Pins and contradictions created by groups are real thermodynamics and remain. A literal "solve, then substitute" repeats dec0707 |
| **B. Pass-2 fold before the solve leaves zero columns** | 7638558 (Source C); 4fb462e (#36: chained renames were not idempotent); c8a3302, 8b4a524 (phantom `K_T`); 3f8289b (#62: allosteric phantom filter; the naive version threw UndefVarError); 58dd6cb (V25: unfolded single-symbol lines stripped from eq_hash → false merges); Known issue 2 (256 opposite folds + 512 one-state folds) | rate_eq:89-103; thermo:237-253 (plain filter); rate_eq:977-993 (allosteric filter re-derives both states); rate_eq:984 (`merge(rename_A, rename_I)`) | test_rate_eq_derivation.jl:1531-1547, 2143-2176 | **Open**: ROADMAP:25-33, 112-121 | **Eliminated** if folding becomes a post-solve rename of dependents whose right-hand side is one symbol. Pass 2's role inside the solve was a workaround for the pre-#62 sentinel (bcf11d8 note) |
| **C. Free-enzyme columns never pivoted** | 226e2cb (free-set heuristic marked every form free → rows silently dropped, `K6=K1` left free); 92e26ad/#62 (−1 sentinel: LDH over-count by 1, broken allosteric detailed balance) | thermo:95-101, 530-554 (fixed) | test_pivot_priority_regression.jl:7-22; test_allosteric_collapse.jl:135-144 | Closed | Irrelevant. Keep the sentinel-free RREF |
| **D. Conformations derived on separate graphs, joined by name** | Pre-#58, derive in A names then rename to I: §5a, 54,674 UndefVarErrors in one LDH run (#56); #54 Fix A; c921dc3, e2bfbb6 (Case-B names). #61 bug 2 (A-state Haldane dependent leaked through the dead I state). #60 (review caught an RE/SS sign inversion and a circular RE pivot). #64 (pruned-graph free set flipped a representative name: +2/+3 instead of +1 in 44 of 276 LDH children). Known issue 1 (766 of 15,219 mechanisms, 1–3 extra parameters) | rate_eq:839-845 (I graph rebuilt through the constructor, so re-oriented); rate_eq:866-874 (named against the pruned `cm`); rate_eq:965 (columns unioned by name); thermo:48-56, 68-69 (free-set patch that exists only for the pruned graph) | test_allosteric_collapse.jl:45-160; test_rate_eq_derivation.jl:1438-1499, 2096-2141, 2178-2203; test_mechanism_enumeration.jl:11185-11248; test/reference/allosteric_golden_reference.txt | **Open**: ROADMAP:10-33 | **Eliminated**, if the I state is the same step list with relabeled columns and the relabeling happens before the solve |
| **E. `:OnlyA` as deletion before the solve** | #70 (1,354 errors; 94% of LDH allosteric mechanisms thermodynamically impossible, the rest silently wrong); #71 (Stiemke completion of the guard; 112 of 691 PFK-P mechanisms failing); #68 + 13a0020 (L-term leak, fix reverted); 5b8bc6d (pruned graph makes a segment bottomless) | rate_eq:805-845; guard thermo:474-509, run in the constructor at types.jl:963-965 | test_types.jl:2190-2431, 418-440; allosteric_ground_truth.jl:82-222 | Guard is complete; deletion remains (ROADMAP:45-63) | **Root cause eliminated** (limits are taken after the solve). New limit risks, see §3 |
| **F. Free-enzyme normalization (`d_free`, cross-weighting, L)** | #68 leak; #69 (formulation 1, three renderings); #71 err1 (ping-pong covalent island → `d_free_I = 0` → NaN); 43f4d23 (#81: weight not shifted by the concentration GCD, 707 laws wrong); 0b42228 (inert I state rendered `L·0ⁿ`) | rate_eq:265-268 (`findfirst` free form); rate_eq:779-796, 1054-1059, 1102-1115, 1152-1154 | allosteric_ground_truth.jl:223, 369-541, 660-757, 991; test_rate_eq_derivation.jl:1219-1238 | **Open**: which form carries L when there are two free forms (ROADMAP:57-59) | **Partial.** No pruning means no fragmentation, but normalization becomes a mandatory step before the limit |
| **G. Orientation, Canonical Step Form, naming representative** | #48 (a hash containing `objectid` made `fitted_params` flip between precompiles); 5ceff74 (SS iso direction rule changed); a93fa56 (tied steps oriented one by one: silent `K_iso=1`, crash at RE); 4435eac, 7da2cf7; d4b8408 (an inhibitor copy read as its reactant changed which constant is dependent) | types.jl:510-617, 628-630; thermo:381-388 (reverse constant eliminated first) | test_types.jl:963-1090, 1665-1900; test_rate_eq_derivation.jl:1330-1376 | Closed for full graphs | **Not addressed.** Orientation still sets names, the reverse-first pivot, and what a group equality equates |
| **H. Dedup / eq_hash** | 7638558; #59 (66% of fits were re-fits; constructor-time merge reverted because enumeration got 15× slower); #61 (self-loop, beam never terminated); #66 (limit cycle); 4dc4877 (canonicalization pruned reachable models); #78 (deleted canonicalization); 58dd6cb; aa79ee0 (39 of 1,075 hashes moved) | identify:262-268 (string key); rate_eq:505-532 | test_rate_eq_derivation.jl:1531-1575, 2143-2176 | **Open**: ROADMAP:72-80 | **Partial.** A deterministic joint solve plus general folding makes regroupings of one graph render identically. Different graphs with the same rate law still need an algebraic key |
| **I. kcat** | #54 B/B′ (single-pattern assert; product monomials up to 22× off); 9478c49; 3779154; #56 §5b → #79 (61 errors); #70/#71 ("no kcat components"); d868889 (V24, up to 9× off); 1dd5516 | rate_eq:593-726 | test_rate_eq_derivation.jl:2205-2361 | **Open**: ROADMAP:90-97 | Mostly irrelevant. The `d_free`-induced share goes with F |
| **J. RE-segment / Cha pathologies** | #52 (numerator 0.60×, 0.48×, v≡0); #48 (65 of 69 non-finite); #79 (bottomless segments, 0/0); #15, #63 (size limits) | rate_eq:154-197, 276-294; types.jl:865-876 | test_types.jl:301, 418; test_rate_eq_derivation.jl:1264-1328, 1744-1788 | Closed | Irrelevant, but limits can create new 0/0 and bottomless states |
| **K. Structural phantoms** | One-subunit L (ROADMAP:64-71); dead-end SS step (ROADMAP:34-41); iso dwell gauge and copy twins (ROADMAP:122-143); `:NonequalAI` degeneracy (#58 §4a, settled by #60's collapse) | Properties of the rate law itself | ROADMAP only | Open | Irrelevant. Limits add one new candidate (R6) |
| **L. Other** | `var"…"` residual names (aa79ee0; ROADMAP:98-104); `1 ^ 2` / `1 ^ 4` in printed MWC laws (deferred since #58; rate_eq:1040, 1145-1150; 2 golden lines); DSL role typing (#66, #80, ec8b025, 0cc13dc, 03f1cab) | — | — | Names and `1^n` open | Irrelevant |

## 2. Deep dives

### 2.1 Kinetic groups × thermodynamics (classes A, B, C)

**What the code does.**
- `_assemble_constraints` makes one column per rendered name. Each step's constants enter the column `get(rename, name(p, mech))` (thermo:378), and `name` renders the group representative's sides (types.jl:1398-1403).
- So a kinetic group is an equality imposed by names coinciding, and an `:EqualAI` constant shared across states works the same way (rate_eq:965).
- A column's pivot priority is that of the last step written into it, "not the highest such score" (thermo:84-94, 386-388). *(Inference)* Splitting or merging a group, or a Pass-2 fold, can therefore change which constant a Wegscheider row eliminates. No commit measures this.
- Pass 2 runs in this order:
  1. Solve with no rename.
  2. Keep the single-symbol ties between two binding K's (rate_eq:99-102).
  3. Rename and solve again.

  The folded symbol is still a group representative, so it stays as a zero column. Both `_dependent_param_exprs` methods strip it afterwards (thermo:252; rate_eq:984-993).

**What went wrong.**
1. **Thermodynamics solved without the groups gave inconsistent laws** (dec0707).
2. **Groups have thermodynamic consequences that only a joint solve finds:**
   - redundant and contradictory rows (885ace1);
   - a constant pinned to 1 (ee462a1);
   - a Theorell–Chance constant pinned to Keq (test_rate_eq_derivation.jl:1565-1574);
   - a silent `K_iso = 1` from tied steps oriented in opposite directions (a93fa56).
3. **Pivot bugs silently dropped rows** (226e2cb, 92e26ad).
4. **The Pass-2 zero column produced phantoms** in plain mechanisms (c8a3302, 8b4a524) and allosteric ones (3f8289b). Fixing them needed a reference guard. Known issue 2 is what remains.
5. *(Inference, unverified)* **The equation string comes from two solves.** `rate_equation_string` prints dependents from the unrenamed solve and parameters from the renamed solve (rate_eq:526-532). A fold also moves the folded column's priority (thermo:89-93), so the two solves could choose different pivots.

**Verdict.**
- The redesign removes 4, 5 and the name-coincidence part of 1–3, if it is one solve over bare constants with group equalities entered as rows, followed by a post-solve fold of every dependent whose right-hand side is a single symbol.
- It cannot remove 2. Those are real consequences of the user's grouping.
- The enumeration already uses the right kernel shape: `_partition_independent_count` computes the cycle basis of the step graph once and regroups the columns on each call (thermo:303-326).
- Pass 2 rarely fires for plain enumerated mechanisms: 58dd6cb found no single-symbol constraint line in any of the harness's 4,850 distinct plain mechanisms. It does fire for allosteric ones (768 cases).

### 2.2 Two conformations, two graphs (classes D, E)

**Three architectures so far:**

- **Before #58: derive in A names, rename to I.** Every I-state relation had to be rebuilt afterwards. The #58 spec calls this "the source of the §5a-class desync bugs (four functions that must stay in lockstep)". Symptoms:
  - 54,674 UndefVarErrors in one LDH run (#56);
  - undefined I-state parameters (#54);
  - a phantom `K_T` (8b4a524);
  - Case-B names for `:EqualAI` dependents whose Haldane relation involves a `:NonequalAI` constant (c921dc3, e2bfbb6).
- **#58–#64: native derivation per state.** The I state uses a pruned graph rebuilt through the constructor (rate_eq:839-845). The two systems are stacked over columns unioned by name (#62, rate_eq:953-975). Bugs:
  - Haldane-dependent leak (#61);
  - name mismatch from the pruned graph's free set (#64);
  - deletion severing cycles (#70);
  - both open Known issues: pruning changes the group representative or the orientation (766 cases), and the per-state renames disagree (768 cases).
- **July: solve-then-limit.** Designed and prototyped on `mwc-solve-then-limit` (b067aa0, 651a821, 9c95f12, 6f0011e). It matched the mass-action oracle on all six uni-uni tag combinations. #71 declined it in favour of the dead-inactive enumeration rule.

**What the redesign changes:**
- One step list, one orientation, one naming. An I-state column is either a relabeled copy (`:NonequalAI`, `:OnlyA`) or the same column (`:EqualAI`).
- Known issues 1 and 2, #64, #61 and #70's root cause can no longer arise.
- The October prototype showed that naming against the whole mechanism alone removes all 766 cases. 14,524 of 15,422 mechanisms stayed byte-identical, and the aggregate pin went from 68 failures to 0 (equalai-prototype-report.md).

**What it must keep:**
- The stacked rows with I-above-A priority (thermo:372-376).
- The collapse cases:
  - With `:EqualAI` catalysis, a `:NonequalAI` binding has no freedom left and collapses to `K_I = K_A` (test_allosteric_collapse.jl:49-56).
  - A steady-state binding collapses its affinity but keeps its speed (test_allosteric_collapse.jl:100-112).
- These are new constraints created by imposing the tags. They come out of the joint solve, not out of relabeling. Relabeling a solved expression instead of the columns rebuilds Case-B.

### 2.3 Free-enzyme normalization (class F)

**What the code does.**
- The MWC combine needs both conformations on one basis (formulation 1).
- The code picks one reference form per state: the first unbound, residual-free form in species order (rate_eq:267).
- It divides by that form's spanning-tree weight only when the weight is a metabolite-free monomial. Otherwise it cross-weights each state by the other state's weight (rate_eq:1102-1115, 1152-1154).

**Failures, all since July:**
- **Pruning fragmented the I graph** (err1): a covalent island gave `d_free_I = 0`, so the rate was NaN. A product-bearing weight emptied the kcat keys (6f0011e findings).
- **The weight missed the concentration-GCD shift**: 707 laws wrong (43f4d23).
- **An inert I state rendered `L·0ⁿ`** (0b42228).
- **A cross-weighting fix regressed LDH** and was reverted (13a0020).

**What the redesign changes:**
- Without pruning there is no fragmentation, so the err1 class goes away.
- But the July spec found that normalization is required before the limit. Its words: "The `:OnlyA` limit is clean only on a dimensionless denominator … skipping normalization makes `v` concentration-dependently wrong by 12–43%". Per-state normalization equals today's cross-weighting algebraically; it is re-expressed, not deleted.
- The choice of reference form when two free forms sit in different steady-state segments is still open (ROADMAP:57-59; the fixture at test_types.jl:418-440 has both E and F).
- *(Inference)* Denis's "1 in the denominator" is the same operation. For a steady-state mechanism the free-enzyme weight is a sum of spanning-tree products, so a literal 1 needs apparent constants (Km-like ratios), not a monomial division.

## 3. Risks the new design introduces

- **R1. Substitution creates constraints.** An imposed group or `:EqualAI` equality can turn a dependent's expression into a new relation among independents: a collapse, K=1, K=Keq, or 0 = log Keq. Evidence: dec0707, 885ace1, ee462a1, test_rate_eq_derivation.jl:1565-1574, test_allosteric_collapse.jl. A two-phase pipeline has to re-pivot each time; one joint solve avoids that.
- **R2. Relabeling after the solve brings back Case-B and the §5a desync.** This is the pre-#58 architecture.
- **R3. 0/0 at bridges.** Sending both constants of a steady-state step to zero zeroes every spanning tree through it (ROADMAP:55-56). A bridge lies in every spanning tree, so the numerator and denominator vanish together. Getting the leading order needs an explicit scaling parameter, applied after normalization.
- **R4. Partner limits and their direction.** An RE binding's `K→∞` forces partner limits through the Haldane relations (thermo:431-441). Infeasible tag sets must be rejected (Stiemke check, thermo:501-508). *(Inference, untested)* From ter-substrate up, where the per-row test and the Stiemke check part (thermo:462-463), the feasible exponent cone can have more than one direction, and the leading-order law may depend on which one is taken.
- **R5. Semantics diverge from deletion.**
  - An SS `:OnlyA` step (both k→0, K finite) and an RE `:OnlyA` binding (K→∞) are different limits.
  - Today both mean deletion plus pruning of what free E cannot reach, and that keeps forms repopulated by reverse catalysis (rate_eq:828-829).
  - The July attack found a pendant form whose mass the limit drops. The guard rejects those mechanisms.
  - *(Inference)* test_types.jl:418-440 is a case where pruning raises a bottomless-segment error that a K→∞ limit may instead derive.
  - Expect fitted lists and laws to move for some mechanisms, each needing an oracle check.
- **R6. A new phantom.** *(Inference)* The constants of an `:OnlyA` steady-state step drop out of the leading-order law. They still have to be solve columns, because their ratio is what keeps the inactive Haldane relation satisfiable (thermo:436-441). If the solve leaves one of them independent, it is a Pass-2-style dummy.
- **R7. Determinism.**
  - Priority must not depend on group membership, storage order or folds, and keys must not use `hash` (#48).
  - *(Inference)* Splitting dependents into Wegscheider and Haldane (goal 12) depends on the basis after RREF. Today's rule, "the right-hand side mentions Keq" (rate_eq:510), is the only basis-free one.
- **R8. Cost.**
  - Thermodynamic work in a constructor once made ter-ter enumeration 15× slower and was reverted (#59).
  - The `AllostericMechanism` constructor already runs a rational nullspace for every `:OnlyA` binding group (types.jl:963 → thermo:482-483).
  - Columns per step instead of per group make the RREF larger.
- **R9. Hashes and names move.** Every eq_hash and many names change. CSVs from older versions will not compare, and loaders for them were declined (ROADMAP:162-163).

## 4. Requirements the new design must satisfy

1. **One joint solve (A, B, D, R1).** One log-space linear system over bare step constants:
   - Haldane/Wegscheider rows from the full graph's integer nullspace (thermo:171-221);
   - one equality row per tied constant;
   - `:EqualAI` tie rows and regulator mirrors.

   Solve it once with exact RREF and no sentinel. A row reducing to 0 = log Keq is an error. A pinned constant becomes a dependent equal to 1.0 or `Keq^c`.
2. **Priority.** A total order over bare constants, keyed on content and independent of group membership, storage order and folds. A group column takes a fixed aggregate score. I-above-A stays.
3. **Folding (B, H).** After the solve, fold every dependent whose right-hand side is a single symbol with exponent 1 into v, using the class representative chosen by that order. This covers binding-K, iso-K and SS ties.
   - Folded symbols are dependents, never columns.
   - Annotate only lines that were actually folded (V25).
   - `rate_equation` and `rate_equation_string` read the same solve.
4. **One step list for both states (D, R2).** The I state is the same step list. Relabel columns before the solve. Never re-canonicalize or re-name a subgraph.
5. **Limits (E, F, R3–R6).**
   - Normalize each state first.
   - Use an explicit scaling parameter and take the leading order of numerator and denominator together.
   - Completeness check: the solved constraints stay finite under the limit (651a821).
   - Write down RE versus SS `:OnlyA` semantics and the behaviour at bridges, pendant forms and two free forms.
   - Keep or move the #70 guard as a cheap enumeration pre-filter.
6. **Orientation (G).** Orient once, on the full graph, by metabolite role (d4b8408). Tied steps of a group share one orientation (a93fa56). Keys stay string-based.
7. **Regression oracles to carry over:**
   - test_pivot_priority_regression.jl;
   - the collapse suite;
   - the §5a closed-transcript check, and detailed balance at Q = Keq for random parameters;
   - the mass-action oracle for 1–3 protomers and ping-pong;
   - "reversing written steps changes nothing";
   - dependent ∩ independent = ∅.

   Add pins for:
   - fitted names ⊆ declared constants (the prototype's T5: 68 → 0);
   - opposite and one-state folds: 256 → 0 and 512 → 0;
   - split-versus-merged regroupings rendering identical strings.
8. **Performance.** No solves added to constructors. Cache the cycle basis per graph for the enumeration moves. Keep `rate_equation` allocation-free and under 120 ns.