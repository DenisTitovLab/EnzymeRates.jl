# Zero-flux steady-state groups (Denis's point 1)

Denis, this is the full report for track T1. It works in Python only, outside the repository; I did not change any package file. Evidence labels: PROVED (the proof is below), EXACT (exact computation over a named set), SAMPLED, CONJECTURE.

## 0. Answer in brief

- A steady-state (SS) **kinetic group** in which every step carries zero net flux is pure waste:
  - the rate depends on its kf and kr only through kf/kr;
  - it adds exactly one unidentifiable parameter ("phantom");
  - converting it to rapid equilibrium (RE) gives the same family with one fewer parameter.
- A zero-flux **step** inside a group that also holds a flux-carrying step costs nothing. The rule must therefore be group-level.
- The exact test is structural and costs O(steps + forms):
  1. Contract each RE segment to one vertex.
  2. Make each SS step an edge weighted by its net uptake of the first substrate, including the RE offsets of its two end forms.
  3. A step carries flux iff its biconnected block contains a cycle of nonzero weight.
- On the export this test matches exact flux evaluation on all 25,850 SS steps, with 0 mismatches.
- The package's `_flux_carrying_groups` is sound but runs in the wrong place: on the parent's full, uncontracted graph. As a result:
  - It passed 148 distinct R4 level-2 flip children that each carry 2 zero-flux groups. Each one is exactly its parent plus 2 phantoms.
  - The split move, which has no screen, created 40 more.
- Rejecting all 188 loses no family; each one's RE-converted twin is already enumerated.

## 1. Definitions

- **Steps and forms.** A mechanism has forms F and steps S. Step e goes from a_e to b_e and has uptake vector σ_e = consumed − released. Its type is RE or SS; group g(e) shares one K (RE) or one pair kf, kr (SS), applied in the written direction.
- **Validity** (as asserted by the export validator):
  - every cycle of the form graph has Σ±σ = n·ρ and Σ±log K = n·log Keq, where ρ = (+1 substrates, −1 products, 0 others);
  - every all-RE cycle has n = 0 (otherwise the rate is infinite).
- **Segments and offsets.** RE segments are the connected components of the RE steps. u(f) is the uptake along an RE path from the segment root to f; it is well defined because RE cycles are balanced.
- **Segment multigraph G.** Vertices are the segments. Each SS step e is an edge seg(a_e)→seg(b_e) with weight ν_e = σ_e + u(a_e) − u(b_e).
  - If both ends lie in one segment, the edge is a self-loop.
  - The weights around a cycle of G sum to the net stoichiometry of its lift to the form graph, which is n·ρ.
  - A cycle is **balanced** iff n = 0. The first-substrate component ν̂_e suffices, because its cycle sum is n.
- **Steady state (Cha, RE limit).**
  - Forms: [f] = y_{seg f}·w(f).
  - One-way rates: α_e = kf·Πx^cons·w(a_e) and β_e = kr·Πx^rel·w(b_e).
  - Flux: J_e = α_e y_{s(a)} − β_e y_{s(b)}.
  - Rate: v = Σ_SS ν̂_e J_e.
  - Around any cycle C of G, Π(α/β)^{±1} = (Keq/Γ)^{n_C}, where Γ is the mass-action ratio.
- **Flux-carrying SS step:** J_e is not identically zero. **Zero-flux SS group:** every step has J ≡ 0. **Mixed group:** it has both kinds of step.
- **Package predicate** (`_flux_carrying_groups`): a group is flagged if one of its steps shares a biconnected block of the **full** form graph (RE and SS steps alike) with an isomerization step (nothing consumed or released).

## 2. Theory (task A)

**Lemma 0.** At steady state the SS non-loop fluxes form a circulation on G: RE fluxes and self-loops cancel inside each segment. The cycle space splits into block cycle spaces, so the restriction of this circulation to each block is again a circulation.

**Theorem 1 (balanced block ⇒ zero flux, any parameters, grouping allowed). PROVED.**

Let B be a block of G all of whose cycles are balanced.

- *Self-loop e.* J_e = y_s(α_e − β_e), and α_e/β_e = (Keq/Γ)^0 = 1, so J_e = 0.
- *Non-loop block.* Every cycle C in B has Π_C(α/β)^{±} = 1. So there are potentials φ on B with log(α_e/β_e) = φ_a − φ_b.
  - Set ψ_s = e^{φ_s} y_s. Then J_e = c_e(ψ_a − ψ_b) with c_e = β_e e^{−φ_b} > 0.
  - By Lemma 0, Σ_{e∈B} J_e(ψ_a − ψ_b) = Σ_s ψ_s·(net outflow at s) = 0.
  - The same sum equals Σ c_e(ψ_a − ψ_b)², so every J_e = 0. ∎

**Lemma 2 (theta). PROVED.** In a 2-connected block containing an unbalanced cycle C, every edge e lies on an unbalanced cycle.

- By Menger, e = (x, y) has two disjoint paths to distinct vertices p, q of C (trivial paths if x or y lies on C).
- Together with the two arcs C1, C2 of C these give cycles Θ1, Θ2 through e.
- Their weights differ by ±weight(C) ≠ 0, so one of them is unbalanced.
- Parallel edges form 2-cycles, and the edge-stack block algorithm handles them.

**Theorem 3 (unbalanced block ⇒ nonzero flux). PROVED under either condition below.**

Take e in a block with an unbalanced cycle. By Lemma 2, e lies on an unbalanced simple cycle. Hill's cycle-flux formula on G gives

J_e·Σ_t W_t A_t = Σ_{simple cycles c ∋ e} (Π⁺_c − Π⁻_c)·Φ_c,

with Φ_c > 0 (spanning forests rooted on c) and Π⁺_c/Π⁻_c = (Keq/Γ)^{n_c}.

- **(i) e is one-way** (all unbalanced simple cycles through e have n of one sign, oriented along e). Every term then has the same sign, so J_e ≠ 0 at every non-equilibrium point, whatever the grouping. For a self-loop, J_e = y_s·β_e((Keq/Γ)^n − 1) ≠ 0.
- **(ii) SS constants are step-independent.** Write kr_e = b_e and kf_e = b_e·K_e, with the b_e independent.
  - Pick an unbalanced cycle C through e and a spanning forest F rooted on C.
  - The monomial b_C·b_F has edge set C ∪ F, a unicyclic spanning graph whose only cycle is C, so it occurs only in C's term.
  - Its coefficient contains ((Keq/Γ)^{n_C} − 1) ≠ 0, so J_e ≢ 0. ∎
- **Coverage.** Condition (i) holds for all 24,124 flux-carrying SS steps of the export (7,014 self-loops and 17,110 edges; 0 two-way steps). EXACT.
- **Open case.** For grouped mechanisms with a two-way step, a Wheatstone-bridge balance forced by group ties could in principle give J_e ≡ 0 inside an unbalanced block. This is a CONJECTURE (no counterexample known); none occurs in the export. Theorem 1 does not depend on it, so a rule built on the predicate is always sound.

**Lemma 4 (lift). PROVED.** Let M* be M with some RE groups flipped to SS, so its segments refine M's.

- Take a simple cycle of G(M). Inside each segment it visits (each at most once), join its entry and exit points by a simple path of flipped edges. Different segments are disjoint, so the result is a simple cycle of G(M*) with the same weight.
- **Corollary (flips are monotone).** A flux-carrying SS step of M stays flux-carrying in every flip-descendant.
- **Corollary (package screen is sound).** Take M* = all-SS; then G(M*) is the full form graph.
  - A flux-carrying step lies on an unbalanced full-graph simple cycle.
  - That cycle must contain an isomerization, because binding-only cycles are balanced (EXACT: 0 exceptions in 12,556 mechanisms).
  - So the package flag is True. Package flag False ⇒ zero flux in M and in all its flip-descendants.
  - EXACT: the package flag equals the all-SS structural status on all 67,235 groups of the 10,071 distinct exported mechanisms.

**Theorem 5 = (a). PROVED, grouping allowed.** Suppose every step of SS group g has J ≡ 0.

1. **Scaling invariance.** Let x be the steady state. Scaling kf_g and kr_g by λ keeps every g-flux at 0 and leaves all other equations unchanged. The Cha system has a unique solution (connected segment graph, positive rates), so x is the steady state for every λ and v is invariant. Hence v depends on the pair only through kf_g/kr_g.
2. **The RE version M′ is valid.** Suppose M′ had an unbalanced all-RE cycle. In G it is a closed walk of g-edges; decompose it and take one unbalanced simple cycle of g-edges. Zero g-fluxes force α_e y_a = β_e y_b along it, so Π(α/β) = 1 = (Keq/Γ)^n. That forces n = 0, a contradiction.
3. **Same steady state.** x solves M′ (g at equilibrium with K_g = kf_g/kr_g, with its free RE flux set to 0). By uniqueness, v_M = v_M′ ∘ π with π(kf, kr) = kf/kr.
4. **Counts.** The constraints involve log kf_g − log kr_g only, so fitted(M) = fitted(M′) + 1. π is a surjective submersion between the constraint manifolds, so the generic ranks are equal, there is exactly one extra phantom, and the images (families) are equal.
5. **Several groups.** Convert them one after another; the steady state x is shared throughout. ∎

**Theorem 6 = (b). PROVED under condition (i) or (ii) of Theorem 3.**

- **Near-equilibrium expansion.** Let X = log(Keq/Γ). Linearizing J_e = β_e(e^{A_e} − 1) gives v = X·G(c) + O(X²), where G(c) = min_φ Σ_e c_e(ν̂_e − Δφ_e)² and c_e = one-way equilibrium flux (Thomson/Kirchhoff; circulations annihilate Δφ).
- **Effect of scaling.** Scaling g multiplies c_h by λ for h ∈ g, and the equilibrium occupancies do not change. By the envelope theorem, dG/dlog λ = Σ_{h∈g} i_h²/c_h, where i_h is the linear-response current.
- **Consequence.** If the scaling were a kernel direction, v would be λ-invariant, hence so would G, and every i_h would have to be 0.
- For a one-way step h, i_h = Σ_c n_c Π_c Φ_c / Σ W A > 0. With step-independent constants, making the conductances of one unbalanced cycle dominant gives i_h ≠ 0.
- **Coverage.** EXACT on the export: all 19,862 flux-carrying SS groups have a nonzero exact scaling derivative at a GF(2^61−1) point, and all 739 zero-flux groups have a zero one.

## 3. Computation over the whole export (task B)

**Method (files in `t1_flux/`).**

- `fluxlib.flux_values` computes each SS step's flux J from engine A's Cha solution over GF(2^61−1), at 3 random points of free-parameter generators and concentrations. A nonzero value certifies J ≢ 0.
- `structural_flux` implements the Section 1 predicate independently: segments, offsets, `_edge_blocks` port, with self-loops handled separately.
- `package_flux_groups` is a line-by-line port of `_flux_carrying_groups` and `_edge_blocks`.
- `group_scaling_derivative` computes dv/dlog λ for (kf_g, kr_g) → λ(kf_g, kr_g), using engine A's reverse-mode gradient. The code checks that the free-coordinate direction maps exactly onto the joint scaling.

**Result:** 0 step mismatches out of 25,850. The zero-flux steps are zero by Theorem 1 and the rest carry a nonzero certificate, so the classification is exact with no probabilistic step.

**Package screen versus truth, over 20,601 SS groups (pkg flag, truly flux-carrying):**

| package flag | truly flux-carrying | groups |
|---|---|---|
| True | True | 19,862 |
| True | False | 654 |
| False | False | 85 |
| False | True | 0 |

**How the package screen differs from the true predicate:**

1. It works on the full form graph with RE edges uncontracted, so it cannot see RE shunts. On the export it equals the all-SS status.
2. It is evaluated on the parent, before the flip.
3. It is only a unit pre-filter in `_expand_re_to_ss`; split, dead-end and init never call it.
4. It shares the needed group-level "any step" semantics.

**Where zero-flux SS groups occur.** "zero steps" and "zero groups" are zero-flux SS steps and groups; "pkg-flagged" counts zero groups the package flags as flux-carrying; "mechs w/ zero group" counts mechanisms with at least one zero group.

| reaction | level | move | mechs | SS steps | zero steps | SS groups | zero groups | pkg-flagged | mechs w/ zero group |
|---|---|---|---|---|---|---|---|---|---|
| R4 | 2 | re_to_ss | 842 | 3226 | 488 | 2486 | 296 | 296 | 148 |
| R4 | 2 | split | 546 | 1466 | 88 | 1202 | 40 | 0 | 40 |
| R5 | 2 | re_to_ss | 842 | 3226 | 488 | 2486 | 296 | 296 | 148 |
| R5 | 2 | split | 939 | 1860 | 88 | 1596 | 40 | 0 | 40 |
| R6 | 2 | re_to_ss (sample) | 150 | 594 | 99 | 441 | 62 | 62 | 31 |
| R6 | 2 | split (sample) | 301 | 429 | 15 | 394 | 5 | 0 | 5 |
| all other strata | | | 9,936 | 15,049 | 448 | 11,996 | 0 | 0 | 0 |
| **all** | | | **12,556** | **25,850** | **1,726** | **20,601** | **739** | **654** | **412** |

- No zero-flux groups occur in R1–R3, at any level 0 or 1, or in any dead_end child. The zero-flux steps in the "all other strata" row all sit inside flux-carrying (mixed) groups.
- The 412 mechanisms are 188 distinct step structures, all of them R4 structures: R5 repeats them and the R6 sample holds 36. The full list is in `zero_flux_mechanisms.tsv`.

**Shape of every case** (`patterns.py`, all 188 distinct structures). The zero-flux steps attach a piece of the segment graph to exactly one other segment, and that attachment is balanced.

- **Flip cases (148).** Two context-split groups are flipped together. This isolates a form or a small piece ({E}, {EA}, {EAQ}, {EP,EAP}, {E,EA,EP}, …), whose SS edges all run to one neighbouring segment with equal weight.
- **Split cases (40).** A single abortive-complex step (EBQ, EAP, EBP, EAQ, or a residual version) is a pendant bridge.

**Example E1 (flip).** Parent R4_00078 (a level-1 split of seed R4_00003):

- E+A⇌EA (g1), E+Q⇌EQ (g2), EA+B⇌EAB (g3), EA+Q⇌EAQ (g4), EQ+A⇌EAQ (g6), EQ+P⇌EPQ (g7), all RE.
- EAB→EPQ (g5, SS).
- fitted 6, rank 6.

Child R4_00470 flips {g1, g2}: fitted 8, rank 6.

- The segments become {E} and X = {EA, EAB, EAQ, EQ, EPQ}.
- E+A and E+Q are two parallel SS edges from {E} to X, both of weight 1. The 2-cycle (bind A, bind Q, release A, release Q) is balanced.
- The chemistry step is an unbalanced self-loop inside X, so all turnover runs through the abortive route EQ+A⇌EAQ⇌EA+Q.
- Both flipped groups therefore carry zero flux. Converting them back to RE returns the parent.
- The package screen passed both groups, because their steps share the full-graph block with EAB→EPQ.

**Example E2 (split).** R4_00065 has the mixed SS group g3 = {EA+B→EAB (catalytic), EQ+B→EBQ (pendant abortive complex)}: fitted 6, rank 6.

- Splitting by context (Q bound) gives R4_00442. EQ+B→EBQ is now its own SS group g6, which carries zero flux: fitted 8, rank 7.
- Its RE twin is the exported level-2 flip mechanism R4_00446 (fitted 7, rank 7).

**Mixed groups are harmless.** 819 mechanisms have zero-flux steps only inside flux-carrying groups, and 743 of them are fully identifiable. A step-level rule would be wrong.

## 4. Converting zero-flux groups to RE (task C)

- **Counts.** 412/412 mechanisms lose exactly #zero-groups fitted parameters, keep the same rank and lose #zero-groups phantoms; 739/739 single-group conversions give fitted − 1 with the same rank. The fitted and rank values agree between engine A, engine B and the package's fitted count, both for the originals and for the equivalents.
- **Exact parameter maps (`exactmap.py`).**
  - Z→C uses K_g = kf_g/kr_g; C→Z uses kf_g = K_g, kr_g = 1.
  - Parameters are rational points generated from each model's free generators.
  - All thermodynamic constraints were checked exactly, and v is equal in exact rationals at 3 concentration points in each direction: **412/412**, so the families are EQUAL.
  - The evaluator agrees exactly with engine A on 150 random exported mechanisms.
- **family_fit (SAMPLED).** 30 pairs, a seeded sample of 19 flip-origin and 11 split-origin R4 cases, 5 laws each way.
  - With 8 starts: 280 of 300 law fits are below 1e-20; the worst is 0.0805.
  - Re-fitting with 32 starts: 7 of the 13 failing directions reached 6.6e-41 to 7.2e-37. For example, R4_01421 Z→C went from 0.0717 to 2.3e-41.
  - R4_01270 Z→C overflowed mpmath memory under all 3 seeds. The other failing directions were stopped to protect shared memory.
  - Final: 26/30 pairs are below 1e-20 in both directions. The 4 unresolved split-origin pairs are optimizer failures, since the exact maps prove equality.
- **Twins in the enumeration** (canonical key invariant to group relabeling):
  - all 327 flip cases (148 + 148 + 31) equal their own parent;
  - the 81 R4/R5 split cases, and 1 R6 split case, equal exported level-2 re_to_ss mechanisms;
  - the 4 remaining R6 split cases equal R4/R5 structures (R4_00655, R4_00761, R4_01004, R4_01179). Their R4 parents' copies are in R6 level 1 (R6_00369, R6_00510, R6_00850, R6_01061). The flip move depends only on the steps, so R6's full level 2 contains them; the R6 level-2 export is a sample.
  - No converted mechanism has a bottomless RE segment (Python port of `_bottomless_re_segment`).
  - **Rejecting these children loses nothing.**

## 5. How flux status evolves under the moves (task D)

| move | effect on SS flux status | proof | export check |
|---|---|---|---|
| RE→SS flip | never flux→zero for existing SS steps; only the flipped groups can be zero in the child | Lemma 4 | 4,172 parent SS steps, 0 violations (0 zero→flux either) |
| split | step statuses unchanged; zero-flux group only when a mixed SS group is bipartitioned | graph and weights unchanged | 4,097 steps, 0 changes; 85 created groups |
| dead-end | all preserved; SS mirrors inherit | a mirror is a parallel edge with equal weight (u(XI) = u(X) + I); new forms join their base segment | 11,846 steps, 0 changes; 557 SS mirrors all match |
| MERGE, RE isomerization | nothing changes | u(X1) = u(X2), so G is unchanged | — |
| MERGE, SS isomerization | edge contraction of G: blocks only split, so steps can only lose flux, and only when a second X1–X2 route with the same weight exists | lift argument | single-step: 16/11,332 flux→zero, all chemistry groups with a dead-end mirror (e.g. R5_04414, where the unmerged mirror EAIinh→EIinhP_res becomes a balanced self-loop); group-wide merge: 0/39 (27 others infinite rate) |
| ELIM, fused type SS if either input is SS | all preserved | series reduction, or re-anchoring inside a segment | 40,058 applications, 0 changes |
| ELIM, RE fused type | edge contraction, like an SS MERGE | — | — |

**Where the check must run:**

- the flip move (on the child, flipped groups only);
- the split move (unit filter);
- MERGE of SS isomerization groups (merge whole groups, then check all SS groups);
- seeds (a test only; 0 failures today).

The check is not needed after dead-end additions or ELIM with the SS-if-either rule.

## 6. Recommendation (task E)

**Rule.** Emit a mechanism only if every SS kinetic group has at least one flux-carrying step.

**Helper sketch** (Julia; reuses `_re_segment_extras`, `_step_sides`, `_edge_blocks`):

```julia
function _flux_carrying_steps(groups::Vector{Vector{Step}}, x::Symbol)
    species, segments, extras = _re_segment_extras(groups)
    idx = Dict(sp => i for (i, sp) in enumerate(species))
    seg = zeros(Int, length(species))
    for (k, members) in enumerate(segments), i in members; seg[i] = k; end
    flags = [falses(length(g)) for g in groups]
    edges = Tuple{Int,Int}[]; weight = Int[]; owner = Tuple{Int,Int}[]
    for (g, group) in enumerate(groups), (j, s) in enumerate(group)
        is_equilibrium(s) && continue
        _, _, lhs, rhs = _step_sides(s)
        a, b = idx[from_species(s)], idx[to_species(s)]
        w = count(==(x), lhs) - count(==(x), rhs) +
            get(extras[a], x, 0) - get(extras[b], x, 0)
        seg[a] == seg[b] ? (flags[g][j] = w != 0) :
            (push!(edges, (seg[a], seg[b])); push!(weight, w); push!(owner, (g, j)))
    end
    block = _edge_blocks(length(segments), edges)
    # per block: BFS potentials; the block is unbalanced iff some edge's weight is not a
    # potential difference; flag every edge of an unbalanced block
    ...
end
```

**Where to use it:**

- `_expand_re_to_ss`: add `all(u -> any(flags[units[u]]), sel)`, with the flags computed on `flipped_groups(sel)`, to `gains(sel)`. A zero-flux set then counts as failed, like the bottomless rule. Keep `_flux_carrying_groups` as the (proven-sound) unit pre-filter.
- `_expand_split_kinetic_group`: compute the flags once per parent and drop every SS-group bipartition in which one part has no flux-carrying step.

**Cost.** O(steps + forms) per call, the same order as the segment-count and bottomless probes already inside `gains`.

**Validated effect.** A Python port of `_expand_re_to_ss` reproduces the export on all 431 R4 level-0/1 parents (1,777 children, 0 differences). With the rule it emits 1,629: 148 removed, 0 added, whether zero-flux sets are extended or dropped. Removing a split unit removes only the sets that contain it.

| scope | mechanisms before → after | phantoms before → after | mechanisms with phantoms before → after |
|---|---|---|---|
| all exported | 12,556 → 12,144 | 2,388 → 1,619 | 1,909 → 1,497 |
| R4 | 1,819 → 1,631 | 616 → 266 | 418 → 230 |
| R4 level 2 | 1,388 → 1,200 | 579 → 229 | 383 → 195 |

**Next generation (R4 port).**

- 64 distinct flux-valid level-3 flip mechanisms are currently reached only through a zero-flux level-2 "bridge", which has its parent's loss but 2 extra nominal parameters.
- Under the rule, all 64 remain reachable, via 3 flips from the surviving level-1 ancestor through identifiable intermediates, i.e. one move later. None is lost.

## 7. Files (`t1_flux/`)

- `fluxlib.py`: flux, predicate, package port and scaling derivative.
- `run_export.py` → `flux_export.jsonl`.
- `convert.py` → `converted.jsonl`, `zero_flux_mechanisms.tsv`.
- `exactmap.py`, `exactmap_sanity.py` → `exactmap.jsonl`.
- `ff_pairs*.py` → `ff_pairs.jsonl`; `ff_retry_all.py` → `ff_retry_all.jsonl`.
- `moves_check.py`, `merge_elim.py` → `merge_elim.jsonl`; `merge_group.py`.
- `oneway.py`.
- `sim_flip.py`, `sim_level3.py`, `sim_lost.py`, `sim_reach.py`.
- `patterns.py`.
- `tables.py` → `tables.json`.