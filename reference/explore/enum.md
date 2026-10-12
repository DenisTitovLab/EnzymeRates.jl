# How enumeration and the mechanism types couple to the derivation

I read all of `src/mechanism_enumeration.jl` (enum). I also read the parts of `src/types.jl` (types), `src/thermodynamic_constr_for_rate_eq_derivation.jl` (thermo), `src/rate_eq_derivation.jl` (rate_eq) and `src/identify_rate_equation.jl` (ident) that it reaches, plus the docs, the HPC progress logs and the tests. I did not run Julia. All costs come from reading the code, test comments and logs. "Inference" marks my own reasoning rather than something I read.

## 1. Dependency table

Hot-path codes:
- **R**: once per run
- **P**: once per parent or move
- **C**: once per candidate inside a minimal-set search, before any child is built
- **K**: once per constructed child
- **F**: once per fresh child in identify

| Function | Defined | Call sites | Why enumeration needs it | Hot | Fate with a graph-centric Mechanism |
|---|---|---|---|---|---|
| `Mechanism(rxn, groups)` → `_canonical_groups` | types:906, 658 | enum:885, 1423, 1933; through `_with` (types:1008) from enum:455, 1187; thermo:482; rate_eq:868; types:1095, 1202 (lifts) | Canonical storage that ==/hash dedup, naming and group validation depend on | K, plus every lift | A cheap builder: species sorted by content, integer edges. Orientation is computed for naming only. Children built by moves skip re-validation. |
| `AllostericMechanism(...)` | types:929 | enum:1427, 1604, 1611; through `_with` (types:1009) from enum:455, 1188, 1693, 1755, 1763, 1844; rate_eq:843; types:1181 | Does the above, then permutes the tags with the groups (951), sorts the sites (952), and runs `_onlya_haldane_violation` (963) | K | Tags become per-group labels. The Haldane check runs only when the catalytic tags change. |
| `_with` | types:1008–1011 | 15 sites (listed above) | Rebuilds a child with one field replaced | K | A plain field replacement. Moves that change only tags or sites (enum:1693, 1763, 1844) stop re-canonicalizing a catalytic graph they did not touch. |
| `_bottomless_re_segment` | types:865 | enum:523 (flips and seed variants), enum:1081 (split); types:670 (every constructor) | The constructor rejects a bottomless segment, so moves test candidates first | C, then K again on a candidate that already passed | A predicate on (graph, RE mask) over cached segments |
| `_re_segment_extras` | types:796 | enum:552, 775, 929, 1085, 1233; types:867; rate_eq:234, 320 | RE segments plus the exponent offset of each form | P; C through the flux and bottomless checks | The core cached analysis. Today it allocates a `Dict{Symbol,Int}` per RE step and per form (types:803, 824). |
| `_forms`, `_step_kind`, `_flat_steps` | types:781, 701, 1296 | enum:330; enum:404; enum:1109 | Species list; seed grouping key; flat positions for the parameter counter | R, P | Fields or edge ids |
| `is_iso`, `is_binding`, `bound_metabolite`, `_is_chemistry` | types:250–267 | About 20 sites, e.g. enum:404, 483, 679, 745, 1583, 1740; thermo:385, 479 | Step classes; tells chemistry from binding when assigning allosteric tags | Everywhere | Per-edge flags, computed once |
| `_partition_independent_count`, `_count_kind` | thermo:303, 263 | enum:1108, 1113, 1122 | Gain test of the split move for a plain Mechanism: cycle basis once per parent, then a rank per candidate | P (basis), C (rank, `Rational{BigInt}` RREF) | Already shaped as "graph plus labels" (`group_of_step`, `kind_of_step`). It should become the single parameter counter, allosteric included. |
| `_independent_param_count` | thermo:257, which calls rate_eq:953 | enum:1130–1131 | Gain test of the split move for an allosteric mechanism | C, and the most expensive path (§2) | Replace with the counter run on one catalytic system holding both A and I columns, without building the child |
| `_onlya_haldane_violation` | thermo:474 | enum:1597, 1754; types:963 | Rejects `:OnlyA` patterns that make a binding-only cycle unsatisfiable | C, then K again | Compute from the cached cycle basis restricted to the kept edges, cached per (graph, groups, tags) |
| `_thermodynamic_constraints` | thermo:171 | Indirect only: thermo:304, 367, 483 | Cycle basis by exact rational nullspace (thermo:121–169), plus the Haldane/Wegscheider class of each cycle | P, C, K | Fundamental cycles of a spanning forest (entries ±1/0), cached on the graph |
| `_inactive_groups` (through `_reachable_from_free`) | rate_eq:805, 779 | enum:1367 | The inactive-state graph that the copy rule reads | P, plus K in some moves | An edge mask over the parent's edges, by BFS on integer adjacency. Not a new mechanism. |
| `_redundant_copy_groups` (defined in enum, built on the rows above) | enum:1358 | enum:66, 1098, 1429, 1765 | Emission rule 2 | P, K | Same logic on cheaper inputs |
| `_free_enz_set`, `_group_rep`, `name(p, m)`, `_step_constants` | thermo:76, 109; types:1398–1403, 1437 | Indirect only: the Haldane labels (thermo:490, rendered eagerly even when nothing is violated); `_independent_param_count` through thermo:366 and rate_eq:95 | Naming and pivot priority inside the counts | C | Names only at render time. A count never needs a name. |
| `_eq_complexity` | rate_eq:319, 332 | ident:460 | V×τ gate before compiling | F | Read from the cached segment graph. Today the allosteric case builds a Mechanism through `_state_mechanism` (rate_eq:868). |
| `compile_mechanism`, which calls `EnzymeMechanism` / `AllostericEnzymeMechanism` | enum:388–389; types:1130, 1201 | ident:342, 461, 491, 846, 855 | Singleton type for code generation | F; also once per LOOCV fold | Outside this slice. Note that each lift rebuilds the mechanism (§2). |
| `fitted_params`, `rate_equation_string` | rate_eq:47, 486–537 | ident:462, 464 | Parameter count and the `eq_hash` dedup key | F | Concrete functions on the graph, cached |
| `==` / `hash` | types:974–982, 83–87 | ident:450, 614, 646, 731–734; enum:343, 352, 1935–1938, 1985 | Dedup | K, F | See §5 |
| Reaction accessors (`substrates`, `products`, `regulators`, `reg_type`, `atoms`, `allowed_catalytic_multiplicities`, `shared_catalytic_site`) | types:459–467, 309, 330 | `reaction(m)` is read at enum:59, 87, 453, 789, 865, 940, 1070, 1081, 1597, 1604, 1611, 1754 | Atoms, which side a metabolite is on, declared regulators, multiplicities | P | Pass `rxn` explicitly. Most moves already take it (enum:1463, 1579, 1646, 1891). |

## 2. Mechanism life cycle, construction counts and cost

**Stages of a search:**
1. **`init_mechanisms` (enum:1931).**
   - Topologies (enum:114) and dead-end seeds (enum:309) are built as flat `Step` lists. Each seed becomes one `Mechanism` (1933), passes the atom assertion (1934), and is deduplicated with `Set`/`unique` (1935–1938).
   - The merged and Theorell–Chance variants are screened on index arrays (`_seed_candidate_screen`, enum:808). A `Mechanism` is built only for a variant that passes (885).
2. **`seed_mechanisms` (enum:1976).** This runs only when regulators are required. It is a breadth-first search through the to-allosteric, add-regulator and dead-end moves, run on workers with `pmap` (1996). It deduplicates by `hash(m)` alone, into a `Set{UInt64}` (1978, 1984–1991).
3. **`_base_tier` (ident:632).** `_degenerate` runs on every seed. A degenerate seed is replaced by its flip children, then `unique!` (646).
4. **`_process_batch` (ident:439).**
   - It skips anything already seen, again by hash alone (449–452).
   - On workers it runs `_eq_complexity`, `compile_mechanism`, `fitted_params` and `rate_equation_string` to get the `eq_hash` (458–471).
   - It fits one mechanism per new `eq_hash`, compiling again (489–498).
5. **`_expand_parents` (ident:596).** `expand_mechanisms` runs per parent under `pmap`. Children are deduplicated in the `parent_of` Dict, where the first parent wins (608–617).
6. **`expand_mechanisms` (enum:1891).** It checks `_assert_emission_rules` once per parent (1896), applies the seven moves, filters with `_respects_reg_type` (1906), and runs `_assert_atom_conserving` on each child (1907–1909).
7. **LOOCV.** One compile per candidate per fold (ident:846).

**Counts on record:**
- **Starting sets:** 239 for bi-bi, 264 when the atoms allow ping-pong, 7 for uni-bi (enum tests :1168, :1175; combinatorics.md). Ter-ter gives 250,855 in about 26 s locally (enum tests :1239–1246), so roughly 0.1 ms per emitted mechanism, screening included.
- **Expansion depth:** bi-bi with every reactant also a dead-end inhibitor, three catalytic moves, two levels: 264 → 6,738 → 128,493 distinct mechanisms (enum tests :10242). The copy rule was checked on "about 218,000 enumerated parents" (enum:1414).
- **HPC logs (my sums over `docs/hpc_results/*/*/progress.log`):** between 4,643 and 102,898 children per run.
  - LDH 2026-07-17: 2,729 parents → 102,898 children.
  - PGK 2026-07-12: 2,411 parents → 55,252 children.
  - PFKP, iteration 1: 276 parents → 9,467 children, from 1,104 base-tier mechanisms.
  - All of these runs predate 0.8.0 (ROADMAP:72–73).

**Constructions per child (inference).** The logged children are only the survivors. On top of them, the code builds:
- children that are built and then dropped (enum:456, 1098, 1430, 1765, 1906; ident:614);
- internal sub-mechanisms (thermo:482 for each Haldane check; rate_eq:843 and :868 for each state part);
- lifts during identify. Each lift reruns the full constructor (types:1095):
  - `compile_mechanism(am)` builds a Mechanism (types:1202);
  - `rate_equation_string` lifts the Sig at least four times (rate_eq:527 through :1007, :528, :301 through :354, and thermo:26) and computes about five cycle bases;
  - the generated `rate_equation`, `metabolites` and `_kcat_forward` lift again for every fitted representative (rate_eq:431, :598; types:1326).

That is roughly 5–6 constructions per fresh plain child before fitting, about 8 more for each fitted representative, and more for allosteric mechanisms. Over a run this comes to roughly 10⁵–10⁶ constructions.

**Constructor work paid on every `Mechanism`:**
- `_canonicalize_step_directions`: a `Dict{Species,…}` that hashes whole bound vectors whose element type is abstract (types:83–87). `_reversed` rebuilds `Step` values (types:535).
- `_orient_tied_steps`.
- String sort keys, four strings per step, for both the step sort and the group sort (types:628–630, 663–666).
- `_assert_distinct_form_names`, `_assert_uniform_groups`, and `_assert_each_reaction_once`, which joins strings for every step (types:709–777).
- `_bottomless_re_segment`.

**Extra work on `AllostericMechanism`:** site sorting with string keys, then `_onlya_haldane_violation`. Whenever a kept group is `:OnlyA`, which holds for every K-type lineage, that check:
- builds a second `Mechanism`;
- computes an exact rational nullspace (thermo:121–169);
- may run Fourier–Motzkin elimination (thermo:402–418);
- renders labels eagerly (thermo:490).

**Duplicated work:**
- A candidate's checks run again in the constructor: the bottomless-segment test at enum:523 and again at types:670; the Haldane check at enum:1597/1754 and again at types:963.
- `_expand_to_allosteric` runs the Haldane check once per binding subset (2^b − 1 subsets), then again in the constructor for each multiplicity (enum:1594–1605).

**The most expensive per-candidate path is the allosteric split gain test (enum:1129–1132).** For each candidate set it builds the child and runs `_dependent_param_exprs(am)`, which does all of this:
- two state parts, where the I state builds both an `AllostericMechanism` and a `Mechanism` (rate_eq:843, 868);
- two Wegscheider renames and two constraint assemblies, each computing a cycle basis;
- one combined solve;
- when either rename is non-empty, a full symbolic derivation of both states (rate_eq:985–993).

**Perspective.** Fitting dominates the run time: one fit takes seconds to minutes (combinatorics.md:91). Memory growth comes from compiled code, about 0.47 MB per distinct fitted mechanism (`docs/superpowers/specs/2026-07-03-ldh-memory-findings.md` at 241f410). Enumeration values can be garbage-collected. Inference: constructor speed matters at ter-ter scale, in the allosteric split and K-type paths, and in test time. It is not what limits HPC runs. The larger win is reusing analyses that a move does not change (§3), not shaving the constructor.

## 3. Graph reasoning duplicated in enumeration

By my classification, about 440 of enum's roughly 1,040 code lines are graph analysis.

1. **Species indexing and incidence are rebuilt on every call.** Sites:
   - `_re_segment_extras` (types:797–809), `_re_turnover_cycle` (enum:722–725), `_seed_candidate_screen` (814–820), `_re_segment_count` (990–997);
   - `_chain_flank_groups` (471–474), `_eliminate_form` (676), `lumping_twin` (873);
   - `_forms_where_free` (1395–1403) and the dead-end move's name maps (1475–1482);
   - `_thermodynamic_constraints`, which keys forms by name (thermo:176–184);
   - `_reachable_from_free`, a fixed-point loop over names (rate_eq:779–796);
   - the position map in `_split_gain_test` (enum:1110–1111).

   About 80 lines.
2. **Three cycle implementations decide the same question.** Each asks whether a cycle runs the net reaction:
   - DFS potentials in `_unbalanced_blocks` (enum:603–639);
   - a weighted union-find in `_has_unbalanced_cycle` (691–710), used by `_re_turnover_cycle` and the seed screen;
   - an exact rational nullspace plus a stoichiometric class (thermo:121–221).

   A cycle has nonzero uptake weight (enum:712–717) exactly when its Haldane multiple `c ≠ 0` (thermo:213–219). Inference: one spanning-forest utility could feed the flux test, the turnover test, conditions V and C, the Haldane check, the parameter counter and the derivation's constraint rows, replacing about 60–80 lines. RE segments are likewise computed three ways: BFS in `_re_segment_extras`, union-find in `_re_segment_count`, and isomerization-only segments for condition C (enum:775).
3. **Candidates are evaluated over flags and labels.** `_seed_candidate_screen` (40 lines) is a hand-indexed fast copy of `_re_turnover_cycle`, `_has_vmax` and `_chemistry_equilibrates_both_sides`. `_minimal_flips` builds candidate `Vector{Vector{Step}}` values (enum:521). The counter already takes label vectors. Inference: with edges, a `BitVector` RE mask and group labels as the representation, the fast path becomes the normal path.
4. **Four definitions of "free enzyme" or "bottom form":**
   - `_compute_free_enz_set` (thermo:61–70);
   - apo-form roots, at rate_eq:780–782 and as `i_free` at rate_eq:267;
   - `_segment_root` (rate_eq:122–126);
   - the segment bottom read from the offsets (types:829–837, 865–876).
5. **Subgraphs are built as full, validated mechanisms** (thermo:482; rate_eq:843, 868). Re-canonicalizing a subgraph re-orients its steps. This is the cause of ROADMAP known issue 1 (`F ⇌ E` names `K_F_to_E` in the active state and `K_E_to_F` in the inactive one; ROADMAP:10–24). Inference: treating subgraphs as edge masks over the parent removes this class of bug.
6. **Names act as identity.** Forms are keyed by rendered name at thermo:65–69 and 176–179, rate_eq:780–796, and enum:1395–1403 and 1475–1511. `_assert_distinct_form_names` (types:709) exists because of this.
7. **Groups are tracked by position.** Tags run parallel to groups and are permuted with them (types:950–951; enum:1188). `_dead_end_child` relies on the new group being last before construction (enum:1422).
8. **Dead-end lattice.** An O(F²) `issubset` scan joins forms that differ by one bound metabolite (enum:366–372).

Estimate (inference): about 150–250 lines across enum, types, thermo and rate_eq are plumbing that a graph core would provide once. The domain rules stay: gauge, twins, conditions V/C, hyperbolic scoring, minimal sets and competition patterns. The biggest reuse comes from invariances:
- the cycle basis does not change under flips, splits or tag changes;
- RE segments do not change under splits or tag changes;
- nothing graph-level changes under site moves.

## 4. Features the derivation must support because enumeration emits them

- **Multi-step groups that share constants:**
  - bindings of one metabolite at several forms, RE or SS (enum:401–407);
  - a plain and a fused binding of one metabolite in one group (enumeration_engine.md, "Every step but a plain binding is chemistry");
  - mirror steps that inherit their counterpart's group and RE/SS flag (enum:1512–1520). Mirrors include SS steps and isomerizations (a ping-pong half-reaction), so a group's steps must share one orientation (types:555–574).
- **Chemistry in all its enumerated forms:**
  - isomerizations, where a seed's first one is steady state and the rest are at rapid equilibrium (enum:228–233), two of them in ping-pong;
  - fused bindings into a merged complex (enum:659–666);
  - Theorell–Chance steps in a group of their own (enum:675–685), RE or SS.
- **Ping-pong residual forms** (enum:122–127).
- **Substrate/product dead-end complexes**, joined by every binding that differs by one metabolite and folded into that metabolite's group (enum:355–372).
- **Competitive inhibitors**, each in one fresh RE group. This includes a reactant's copy, which shares the reactant's name and concentration but has a different role (types:49–57; enum:1386–1394).
- **Context splits** by ligand, conformation or residual (enum:1161–1176). A zero-flux part reverts to RE (enum:1197–1200), so one SS group can become an SS part plus an RE part.
- **Allosteric:**
  - catalytic tags `:OnlyA`, `:EqualAI`, `:NonequalAI` (types:918), with K-type and V-type variants (enum:1593–1614);
  - chemistry groups relax to `:NonequalAI` together, binding groups one at a time (enum:1747–1749);
  - tags on SS groups after flips; the copy group tagged `:EqualAI` (enum:1427); split parts inherit their tag (enum:1188);
  - regulatory sites with any of the four states, merged multi-ligand sites (enum:1824–1848), each site at the catalytic multiplicity (enum:1687);
  - one name both an allosteric regulator and a dead-end inhibitor (enum:1623–1626);
  - catalytic multiplicity above 1 only with a hyperbolic scheme (enum:1580–1581, 456).
- **The inactive state as forms reachable from free E** once the `:OnlyA` groups are removed (rate_eq:805–813). The copy rule reads this independently of how the derivation gets the I-state law. If the derivation moves to limits, enumeration still needs an explicit "live in I" mask.
- **Hand-written only, never emitted:**
  - conformations other than `:E` (enumeration builds only `:E`: enum:126, 180; complexes inherit theirs at enum:339, 1509–1510);
  - a merged complex on the substrate side;
  - a step with more than one metabolite on a side, other than Theorell–Chance;
  - SS dead-end-only steps (ROADMAP known issue 3).
- **Structural predicates pinned against the derivation, which make good oracles for a new one:**
  - `_hyperbolic_catalysis` (developer.md:213);
  - counter equals `fitted_params` (enum tests :6558, :8873–8880);
  - V×τ (rate_eq:304–318).

## 5. Equality, hash and canonical form

**What == means today.** `==`/`hash` compare fields in order (types:974–982); `Species` compares (bound, conformation, residual) (types:83–87); the naming cache takes no part. This works only because the constructor canonicalizes:
- a binding is stored with its metabolite consumed (types:223–225);
- every other step is oriented by three tiers, and Tier 2 reads the whole graph (types:510–617);
- steps and groups are sorted by string keys (types:663–666);
- sites and ligands are sorted (types:142–143, 952–953);
- metabolite lists and residuals are sorted (types:35, 72, 213–214).

The property is pinned at `test_types.jl`:980–993.

**Who relies on it:**
- `unique!`, `Set` and `Dict` (enum:1935–1938; ident:614, 646, 731–734), and the tests' `Set(children) == Set(expected)`;
- two hash-only `Set{UInt64}` dedups (enum:1978–1991; ident:449–452). Hash quality is therefore a correctness property, not just a performance one;
- `Step` and `Species` hashes inside the moves (enum:183, 215, 343–353, 1174).

**Why step order is load-bearing today:**
- `_solve_dependent_set` breaks priority ties by column order, which is flat step order (stable `sortperm`, thermo:537);
- a column takes the priority of the last step in canonical order (thermo:377–389, documented at thermo:89–93);
- the cycle-basis rows depend on order (thermo:176–211).

**What a graph-based equality must preserve:**
1. **No isomorphism test is needed.** Every form has an intrinsic label (role-aware bound multiset, conformation, residual), and enumeration builds forms from that content (inference). The keys must stay role-aware: an inhibitor copy is a different vertex from the reactant-bound form (types:49).
2. **Edges are identical up to reversal.** Orientation exists for naming (k forward vs reverse, the direction of K; types:1398–1403). It must be one function of the whole graph, computed once and carried into subgraphs.
3. **Groups compare as a set of edge sets**, with flags and tags attached to the group, not to a position.
4. **Regulatory sites compare as a multiset.**
5. **The reaction should leave ==.** This deletes about 95 `_testhelper_on_reaction` lines in the enumeration tests (definition at 152–153). Within a run every mechanism carries `prob.reaction` anyway (inference from ident:732–734 and enum:1463–1522).

**Two options:**
- **(A) Canonical integer order.** Sort species by content key (Symbols compare lexicographically), edges by integer tuples, and groups by their smallest edge. Then == is array equality and the hash can be computed once and cached. Determinism of the derivation comes for free.
- **(B) Insertion order with an order-insensitive ==/hash.** Construction is cheaper, but then the derivation must not depend on order: name-keyed tie-breaks and content-keyed group priority, not last-write-wins. Otherwise two equal mechanisms derive different-looking equations and get different `eq_hash` values, which is exactly the problem in goal 5. A commutative hash must also stay strong, because two of the dedup sets are hash-only.

**Either way:**
- no ordering may use `objectid`-based hashes (types:619–626);
- emission order must stay deterministic, because it decides CSV row order, `_select_count!` ties (ident:299–304), first-parent-wins (ident:614) and which mechanism represents each `eq_hash` (ident:474–477).

## 6. Risks and constraints

1. **Parameter-count drift.** The beam buckets mechanisms by fitted count, the split move's gain uses the same count, and tests pin both. Fixing ROADMAP known issues 1–2 changes allosteric counts, so the tests need re-pinning and search trajectories will change.
2. **Hash-only dedup** (enum:1978; ident:449) breaks silently if the hash gets weaker.
3. **Orientation depends on the whole graph** (Tier 2, types:605–617). Recomputing it on a subgraph causes known issue 1.
4. **Mirror groups need one orientation shared by the whole group.** Goal 1 ("substitute the groups afterwards") needs the same: a group shares its (k forward, k reverse) only in one physical direction.
5. **Goal 16 (`::Inh` as a naming suffix).** The inhibitor role must stay part of vertex identity. It is read by:
   - the twin, gauge and copy rules;
   - the hyperbolic filter (enum:926–927) and flip eligibility (enum:449–451);
   - Tier 1 scoring (types:516–518).
6. **Subgraph validation hazard.** The Haldane check builds a fully validated Mechanism from the kept groups (thermo:482). A pruned graph could violate a constructor invariant such as the bottomless-segment rule. I did not verify that this happens (inference). Masks avoid the question.
7. **Removing the reaction from Mechanism.**
   - Enumeration still needs atoms (enum:11–44, 124), declared regulators, multiplicities and the shared site (enum:323, 1470, 1581, 1652).
   - Inference: the derivation needs only roles and the net reaction (thermo:195–205; rate_eq:114, 283), which can be read from the graph.
   - The Sig embeds atoms and allowed multiplicities (types:1061, 1064), so the same mechanism declared on different reactions compiles twice (inference).
8. **Mechanisms cross `pmap`** (enum:1996; ident:458, 598). Caches must either be excluded from serialization or be cheap.
9. **Test coupling.**
   - The 442 inline fixtures in the enumeration tests must build values that are == to the enumerator's, through the single DSL.
   - Tests round-trip `mechanism_type` Sig strings (`test_identify_rate_equation.jl`:257, :434; `test_pivot_priority_regression.jl`:17). Goals 13 and 18 need a replacement round-trip format.
10. **Where the Haldane check lives.** Under the limit-based derivation (goal 7), `_onlya_haldane_violation` is the feasibility condition of the limits. It should become one predicate on (graph, groups, tags), shared by the derivation and enumeration.