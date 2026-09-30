# Identifiable mechanism enumeration: findings and plan

Date: 2026-09-28. Status: findings computed by seven research tracks with two independent exact
engines; 77 of the tracks' 147 claims rechecked by independent skeptics; every sentence below
fact-checked against the track reports. Scope decided (A + B + C); A designed, its spec follows;
B and C not yet designed. This document supersedes
`2026-09-24-structural-identifiability-of-enumerated-mechanisms.md`: Denis asked that nothing be
taken from that survey, so every finding below was rederived with new engines and new code, and
no track used the survey as evidence.

## Goal

The enumerator emits mechanisms with more fitted parameters than steady-state data can
identify. The beam, the parsimony filter and the 1-SE rule charge such a mechanism for
parameters its rate law cannot see. Denis's concern is that such a mechanism then loses to
cheaper rivals, so the search misses mechanisms that would win; this study did not measure how
often that changes a selected model. Denis's goal: the enumerator emits mostly identifiable
mechanisms, or at least stops emitting the clearly non-identifiable classes, and it reaches the
textbook families (a single central complex, Theorell–Chance, ping-pong) at their true parameter
counts. For ping-pong that is possible only as two identifiable 6/6 families plus a 5/5 boundary,
or as the non-identifiable 7/6 textbook form (point 6).

## Terms

- **Fitted**: the constants left free after the Haldane and Wegscheider constraints, with Keq
  and E_total known.
- **Rank**: the generic rank of ∂v/∂(free log-parameters) over all positive concentrations.
  **Phantoms** = fitted − rank. **Identifiable**: no phantoms.
- **f/r**: fitted/rank.
- **Family**: the set of rate laws a mechanism reaches over positive parameters. Two families
  are **EQUAL** when each mechanism reaches every law of the other at finite parameters.
- **Type strings** list each step's type in order: R for rapid equilibrium (RE), S for steady
  state (SS). Ordered bi-bi RRSRR is E+A⇌EA (R), EA+B⇌EAB (R), EAB⇌EPQ (S), EPQ⇌EQ+P (R),
  EQ⇌E+Q (R), the package's RE seed.
- **MERGE(X1, X2)** contracts a ligand-free step X1⇌X2 into one form X, where the step is alone in
  its kinetic group and is the only step joining X1 and X2. Every other step keeps its
  metabolites, type and group. (Point 1's single-step merges also allow group-mates, which stay;
  point 10's R5 counts merge such a group whole.)
- **ELIM(X)** replaces a form with exactly two steps, Y⇌X⇌Z, each in its own group, by one fused
  step Y⇌Z in a new group carrying both steps' metabolites. The fused step is SS if either step
  was SS and RE otherwise, unless a variant sets its type explicitly. It is a Theorell–Chance
  step when one metabolite is consumed and another released.
- **Clean merge** (point 10 onward): every kinetic group with a step at X1 or X2 holds only steps
  at that same end, all oriented the same way. In point 9, **clean** means non-degenerate.
- **Evidence labels**: PROVED (a proof exists in a track report or a skeptic's verdict), EXACT
  (exact computation over a named finite set), SAMPLED (random laws fitted numerically),
  CONJECTURE (not yet checked).

## How the findings were verified

- **Two independent engines**, written in Python with no shared code:
  - A: Cha segments plus King–Altman, with an exact rank over GF(2^61−1) at random points and an
    exact symbolic rank over Q.
  - B: brute-force mass action over all enzyme forms, with the exact rapid-equilibrium limit for
    its laws, its own cycle enumeration, and a 150-digit numeric rank at Λ = 10^60. Its
    singular-value gap is at least 10^52.9 on the export and never ambiguous.

  They agree on 77 curated mechanisms (69 valid; both reject the other 8) and 500 random ones
  (384 valid; both reject the other 116). On the valid ones they agree on fitted counts, ranks,
  term supports and exact rate values at 5 points. They also agree on the fitted count and rank
  of all 12,556 exported mechanisms, and both reproduce Segel's textbook laws (IX-8, IX-45,
  IX-87, IX-140).
- **The package agrees with both**: fitted counts on all 12,556 exported mechanisms, ranks on the
  250 the package itself ranked.
- **The export**: `init_mechanisms` plus two levels of the three catalytic moves
  (`_expand_re_to_ss`, `_expand_split_kinetic_group`, `_expand_add_dead_end_regulator`) on six
  reactions:
  - R1: uni-uni S ⇌ P.
  - R2: R1 plus a competitive inhibitor I.
  - R3: R1 with S and P also declared competitive inhibitors.
  - R4: bi-bi A[CX] + B[N] ⇌ P[C] + Q[NX], which admits ping-pong.
  - R5: R4 plus a competitive inhibitor I.
  - R6: R4 with A, B, P and Q declared competitive inhibitors (level 2 sampled: 4,169 of
    28,304).

  "Bi-bi" below means R4: 1,819 mechanisms (62, 369 and 1,388 at levels 0, 1 and 2), 418 of
  them with phantoms, 616 phantoms in all.
- **Seven research tracks.** For each track an independent skeptic wrote its own code and
  rechecked claims: all 17 for point 1 and 10 selected claims for each other track (77 of 147).
  It confirmed 63, weakened 12, refuted 1 and did not rerun 1. Where a skeptic weakened a claim,
  this document states the weaker form.

## Findings, by Denis's point

### 1. Zero-flux SS groups

- An SS kinetic group none of whose steps carries net steady-state flux adds exactly one
  phantom: the rate law depends on its kf and kr only through kf/kr. Making the group RE gives
  the same family with one fewer parameter. PROVED; EXACT on all 412 affected exported
  mechanisms.
- **The rule is group-level.** A zero-flux step inside a group that also holds a flux-carrying
  step adds no phantom: scaling the group's kf and kr together changes v (PROVED when the
  flux-carrying step is one-way, as all exported ones are; EXACT on all 19,862 flux-carrying SS
  groups of the export). A step-level rule would reject the 819 exported mechanisms whose
  zero-flux steps all sit in such mixed groups, and 743 of them are identifiable.
- **Exact structural test**, O(steps + forms): contract each RE segment to a node; give each SS
  step the net uptake of the first substrate, including the RE offsets of its two end forms; a
  step carries flux iff its biconnected block contains a cycle with nonzero net uptake.
  - The zero-flux direction is PROVED for any parameters and grouping (a potential/Thomson
    argument).
  - The flux-carrying direction is PROVED (Hill's cycle fluxes and a theta argument) when the
    step is one-way or the SS constants are step-independent, and is a CONJECTURE for a two-way
    step in a shared group. All 24,124 flux-carrying SS steps of the export are one-way.
  - A rule built on the test can therefore miss a rejection but never rejects wrongly. The test
    matches exact flux on all 25,850 SS steps of the export (EXACT).
- **Today's enumeration emits 188 zero-flux bi-bi mechanisms.** `_flux_carrying_groups` is sound
  for today's step kinds, but it runs on the parent's uncontracted graph before the flip, and the
  split move never calls it. It lets through 148 flip children, each exactly its parent plus 2
  phantoms; the split move adds 40 children that isolate a dead-end abortive binding in its own
  SS group. Rejecting all 188 takes the bi-bi phantoms from 616 to 266 and loses no family: each
  one's RE-converted form is its parent or another exported mechanism (EXACT).
- **The screen breaks on fused steps.** It flags only groups that share a form-graph block with an
  `is_iso` step, so a mechanism with no isomerization step left, such as a fully merged or
  Theorell–Chance mechanism, gets no flip children (two exact counterexamples: 0 flip children,
  against 2 under the correct test).
- **Reachability.** In the bi-bi flip-and-split move graph, 64 identifiable level-3 mechanisms
  (all 9/9) are reached today only through a zero-flux level-2 child. Counting a zero-flux flip set
  as failed inside the existing Apriori search gives exactly the graph of rejecting zero-flux
  children: the 64 appear only at level 4, and 272 identifiable level-4 mechanisms are delayed the
  same way. Searching the minimal flux-carrying supersets of each failing set delays none (0
  identifiable mechanisms lost at levels 3 and 4, with or without a cap on extra groups). At depth
  4, 72 dropped split children have RE-reverted forms that the moves never reach; emitting the
  reverted form instead of dropping the child keeps them.
- **Flux under the moves.** A flip never turns a flux-carrying step zero-flux (PROVED). Splits
  and dead-end additions leave each step's status unchanged. A single-step MERGE of an SS
  isomerization can lose flux: 16 of the 3,818 valid ones on the export, each in a chemistry group
  with an unmerged dead-end mirror; merging whole groups lost no flux in 39 valid cases. The check
  must therefore run on every child of every move. The four allosteric moves were not analyzed,
  and whether the check runs on the A-state projection or on the combined A/I graph is still
  open.

### 2. Substrates and products declared as competitive inhibitors

- A dead-end copy L\* of a substrate or product binds with L's own concentration. At depth ≤ 2 it
  adds a phantom only when it creates no new complex, that is, when every complex X·L\* it makes
  has the composition of an existing form: the productive complex, an init abortive complex, or
  another copy. Copies that create at least one new complex never added a phantom there (0 of
  26,292 R6 events, EXACT). Copies of foreign inhibitors never add one (PROVED).
- **Denis's example is the textbook case** (non-productive binding): K_S and K_S,inh enter the
  denominator only as K_S + K_S,inh, and K_S enters the numerator only through k·K_S, which a
  rescaled kcat absorbs (K_S' = K_S + K_S,inh, k' = k·K_S/(K_S + K_S,inh)). The child's family
  equals the parent's (PROVED; 4/3 against the parent's 3/3).
- **Not every duplicate is a phantom.** Among the children that create no new complex, 29% are
  phantoms (3,082 of 10,552 R6 events). The other 7,470 are identifiable. For the 5,900 inside the
  proof's scope (all sites RE-bind L through one group, no SS mirror), a shared kinetic group pins
  the productive complex ("tie-only"), and their families lie inside a kinetic-group split of the
  parent (PROVED). The remaining 1,570 have no proof.
- **Rule "a dead-end copy must create a new complex"**: at depth ≤ 2 it removes 10,552 of the
  36,844 R6 dead-end events (28.6%), every phantom-creating one among them. 7,470 of the removed
  events are identifiable. Some removed phantom children reach laws their parent cannot (PROVED for
  R6_03189 against R6_00179 and for a hand example), and up to 278 phantom events lie outside
  every containment proof. At level 3 the rule also discards 12 duplicate-only children that are
  more identifiable than their parents.
- **Uni-uni.** Every S or P copy the dead-end move can build today (only at E; it cannot build
  second-site complexes such as ES·S) is a phantom, so under the rule a required S or P inhibitor
  leaves no seed. Merged uni-uni seeds from C admit an identifiable placement (E+S→X SS, X⇌E+P RE,
  E+S⇌ESinh RE is 3/3).
- **Two gaps the skeptic found:**
  - The split move can split a kept copy's group so that one part's sites all duplicate existing
    complexes, which recreates the phantom (42 of the 184 exported level-2 split children of
    dead-end parents the rule keeps, engine B; about 285 across all of level 2 by extrapolation,
    SAMPLED). Either the check must hold for every child, or the split move must not split copy
    groups.
  - Composition misses a complex that matches an existing form only through an RE isomerization
    (ping-pong EB_res ⇌ EQ): 3 new-complex phantoms among the 70,687 level-3 dead-end children of
    the 4,169 sampled level-2 mechanisms, one of which (R6_28303 + Q\*, 8/8 → 9/8) has a lineage
    the rule keeps. Also counting a site as a duplicate when it matches an existing form's (RE
    segment, net uptake) catches every phantom-creating dead-end event at depth ≤ 2 (3,082,
    exhaustive) and in the level-3 sample (5,737), at the cost of 14 identifiable events beyond
    those the composition key already removes. Both keys are needed: the uptake key alone misses
    the SS-binding duplicates.

### 3–5 and 8. Merging the central complex; Theorell–Chance

- **MERGE and ELIM keep thermodynamic consistency** (PROVED). With private groups, MERGE removes
  1 parameter from an RE isomerization and 2 from an SS one. With the default fused type, ELIM
  removes 2 when both flanking steps are SS and 1 otherwise (PROVED; EXACT on the catalog); other
  choices of fused type remove 3 or 0 parameters and change which forms share an RE segment.
- **Chain lemma** (PROVED; EXACT on 227 pairs). Let X1's only other step bind a metabolite into X1
  and X2's only other step release one out of X2, with all three steps in private groups. The
  rest of the mechanism then sees the chain only through four numbers: two flux coefficients and
  two occupancy coefficients. (The skeptic added the orientation condition: a free-enzyme
  isomerization, as in iso uni-uni, does not qualify.) Consequences:
  - Of the chain's eight RE/SS types, only RE binding / SS isomerization / RE release (the
    sequential seeds' pattern, RSR) has no phantom. SSS carries 2 phantoms; SSR, RSS, SRS, RRS,
    SRR and RRR carry 1.
  - Each type equals exactly one merged form, its **class representative**: SSS, RSR, SSR, RSS
    and SRS merge with both flanks SS; RRS with RE entry and SS exit; SRR with SS entry and RE
    exit; RRR with both RE. A merge that keeps the original flank types is exact only for SSS
    chains and RE isomerizations.
- **Denis's examples** (EXACT, family relations PROVED):
  - Uni-uni all SS, 5/3, equals the merged E+S⇌X⇌E+P, 3/3.
  - Ordered bi-bi all SS with two central complexes, 9/7, equals the merged form, 7/7, with the
    same 11 denominator terms.
  - Theorell–Chance all SS, 5/5, lacks exactly the ABP and BPQ terms. It is a limit of the
    merged family, never reached at finite parameters.
- **ELIM of the merged central complex**, with the default fused type, removes exactly the
  denominator terms that only the central complex's occupancy produced (PROVED; the skeptic
  reproduced the listed cases): {ABP, BPQ} for ordered all SS, {AP} or {BQ} for a ping-pong
  half-reaction.
- **RE sandwich** (PROVED for chains): Y⇌X1 RE, X1⇌X2 SS, X2⇌Z RE equals one merged X with SS entry
  and SS exit, at the same count. It does not hold in general when X1 or X2 has more than two
  steps: random SSRRSRRSS is 11/11 while its SS-flank merge is 13/13, and the random RE seed
  RRRRSRRRR (7/7) equals its SS-flank merge only at 9 fitted.
- **Today's enumeration duplicates families.** The ordered RE seed RRSRR (R4_00004, 5/5), its flip
  children RRSSR and RSSRR (R4_00081, R4_00082, 6/5) and its grandchild RSSSR (R4_00481, 7/5) are
  one family. 108 of the 1,819 bi-bi mechanisms equal a cheaper exported mechanism (EXACT).
- **On the export**, replacing every qualifying chain by its class representative removes 159 of
  the 616 bi-bi phantoms and makes 127 mechanisms identifiable, with no family lost (EXACT).

### 6. Ping-pong

- The six-step all-SS form (E+A⇌EA⇌FP⇌F+P, F+B⇌FB⇌EQ⇌E+Q), 11/6, equals the textbook four-step
  form, 7/6 (PROVED). None of the 63 valid six-step RE/SS assignments is identifiable.
- **The textbook form is not identifiable.** Let kf3 be the forward constant of (EA=FP)⇌F+P and
  kf6 that of (FB=EQ)⇌E+Q. The law fixes 1/kf3 + 1/kf6 (that is, 1/kcat in the forward
  direction) but not its split: which half-reaction limits turnover is invisible (PROVED; v is
  bit-identical along the null direction).
- **No mechanism in the MERGE/ELIM space is both identifiable and EQUAL to it.** Its family is
  exactly the union of two identifiable 6/6 families, RSSS and SSSR (A binding RE, or Q release
  RE, the rest SS), and a 5/5 boundary, RSSR. The sign of D_AB·D_P − D_AP·D_B decides which piece
  a law belongs to (PROVED; the skeptic built exact maps for 19 random textbook laws, 12 into RSSS
  and 7 into SSSR, and 10 converse maps from RSSS).
- **The textbook family is already reached today** as R4_00429 (level 1, 7/6), the lumping twin
  of the four-step all-SS form.
- **All 7 ping-pong seeds are degenerate**: at zero products their rate does not depend on B
  (PROVED). All seven have an RE second chemistry step EB_res⇌EQ. In R4_00061 and R4_00062 (6/5)
  every step at EB_res and EQ owns its group, and merging that RE step removes the phantom exactly
  (5/5, EQUAL). In R4_00056–60 shared binding groups make the seed 6/6, and the merge loses rank.

### 7. Random order

- Fully SS random bi-bi, 15/15, merges to 13/13 with identical terms (18 in the numerator, 48 in
  the denominator) and consistent thermodynamics (EXACT). The family is not EQUAL: merging an SS
  isomerization is the same as making it RE and then merging, so the result is the parent's
  fast-chemistry limit. It lies in the closure of the parent's family and has rank 13 against the
  parent's 15, a strictly smaller, lower-dimensional family (PROVED).
- Of the 256 valid random assignments with an SS isomerization (EXACT):
  - 49 merges are invalid (an all-RE catalytic cycle);
  - of the 207 valid merges, the rank drops by 2 in 176, by 1 in 30 and by 0 in 1;
  - the terms stay identical in 56.
- **Merging an RE isomerization is exact** (−1, EQUAL) whenever every step touching X1 or X2 owns
  its group, and every one of the 207 valid random assignments with an RE isomerization carries
  a phantom (EXACT). A shared group can break it: in package ping-pong seeds R4_00056–60 a binding
  group joins a step at EB_res or EQ to a step elsewhere, the seeds are 6/6, and their merges
  (5/5) lose rank.

### 9. New seeds

MERGE of a seed's chemistry leaves four binding groups, all RE, and an infinite rate. With k of
them SS the count is 3 + k; sequential seeds have 5 and ping-pong seeds 6.

- **k = 1** (4 parameters): all 248 variants are degenerate; the flipped ligand never saturates
  (EXACT on all 248; the proof covers only seeds whose flipped group holds no abortive step).
- **k = 2** (5 parameters, the sequential seed count): 130 of 372 variants are clean, all
  identifiable.
  - 8 are the seed itself: when each group next to X holds only its step at X and both are SS,
    the merged form equals RE entry, SS chemistry and RE exit (lumping, PROVED; seeds 00003,
    00004, 00006, 00007, 00018, 00022, 00025, 00026). In the other 12 ordered seeds an adjacent
    group also holds an abortive step, and the same SS pair gives a new family (EXACT).
  - The other 122 give 114 new families, none EQUAL to any mechanism at levels 0–2 (EXACT:
    signatures, plus exact containment tests on the 20 variants whose signature matches their own
    seed). Some lie inside larger mechanisms: one merged ping-pong variant (PP-BP) inside the 7/6
    level-1 R4_00429 at finite parameters, and the ordered representatives inside 6-parameter
    mechanisms only as limits (both SAMPLED).
  - Example: ordered A-first/P-off-first (R4_00004) with SS A binding and SS P release has, at
    zero products, the full steady-state ordered forward law V·A·B/(Kia·Kb + Kb·A + Ka·B + A·B)
    at 5 parameters (PROVED).
- **Random/random** has no clean 5-parameter variant; the minimum is 6 (three SS groups). Each of
  the 28 ordered/random seeds also has one clean three-group set at 6 parameters.
- **Theorell–Chance** is clean only with all three groups SS (EXACT). That gives 5 parameters for
  the 8 ordered seeds whose central-complex groups hold one step each, and 6 or 7 for the other
  12. Among the undecorated seeds, the Theorell–Chance form of A-first/P-off-first equals that of
  B-first/Q-off-first (PROVED, explicit map), and A-first/Q-off-first equals B-first/P-off-first
  (EXACT on 3 random laws each way), so the 8 five-parameter variants give 6 families.
- **Merged ping-pong** with A binding and Q release SS, or B binding and P release SS, is clean at
  5/5, one parameter fewer than today's degenerate seeds (14 variants, 2 per seed). For the
  undecorated seed R4_00062 both have the classic ping-pong zero-product form
  V·A·B/(Kb·A + Ka·B + A·B) and the textbook denominator terms, but neither is the textbook family
  (7/6): PP-BP lies inside it, PP-AQ only partly (SAMPLED). Four half-Theorell–Chance variants
  (R4_00057 X1, R4_00061 X2, R4_00062 X1 and X2) are also clean at 5/5.
- **Degeneracy check.** Three structural conditions (V: a Vmax exists in both directions; C: no
  chemistry node is left by RE releases of both a substrate and a product; I: every ligand's
  binding is visible) match exact degeneracy flags on 3,056 mechanisms with 0 disagreements. With
  the zero-flux rule of point 1 in force, V and C gave the same verdict as V, C and I on every
  flux-valid mechanism tested (about 5,600; EXACT, not proved). V wrongly rejects clean mechanisms
  that carry dead-end copies of substrates or products (100 merged R3 and R6 variants), so it
  cannot be applied to those as it stands.
- **Proposed bi-bi seed set: 179** (today 62), all 5/5, clean and pairwise distinct (the skeptic
  checked distinctness with exact containment tests): the 55 sequential seeds, 114 merged
  two-group variants, 6 ordered Theorell–Chance and 4 ping-pong half-Theorell–Chance variants. An
  undecorated ordered topology with n substrates and m products gains n + m − 2 merged seeds
  (EXACT from uni-uni to ter-ter); a bi-bi ordered seed whose groups next to X also hold an
  abortive step gains n + m − 1 = 3.
- **Dropping the degenerate ping-pong seeds loses reachability.** All 132 clean ping-pong
  mechanisms at levels 1–2 descend only from them. 83 are identifiable, and none of these equals a
  merged variant without splits (no shared signature, EXACT); among the non-identifiable ones,
  R4_00429 (7/6) equals R4_00062 merged with all four groups SS.

### 10. Merging identifiable mechanisms

- MERGE of an SS isomerization always removes 2 fitted parameters (EXACT on 4,412 valid R4 and R5
  merges). On an identifiable parent the rank always drops, by at least 2, so the merge is never
  EQUAL to its parent. A clean merge lies in the closure of the parent's family (PROVED), often
  only as a limit. For obstructed merges (1,072 of the 1,189 R4 merges of identifiable parents),
  whether the merge lies even in that closure is open.
- The terms stay identical in 124 of the 1,464 valid bi-bi merges: all 24 whose X1 and X2 carry
  only SS steps, and 100 of the 1,440 others.
- A clean merge is exactly the parent with its central isomerization at rapid equilibrium
  (PROVED). Example: R4_00079 (6/6) merges to the rate-limiting A-binding law at 4/4. That law is
  degenerate: at zero products v = kf·A, with no Vmax and no dependence on B. So are all 220
  four-parameter merges of level-1 parents (point 9, k = 1).
- As an extra child that keeps its parent, MERGE adds 1,464 bi-bi mechanisms at depth ≤ 2 (+80%),
  1,400 of them with a (terms, rank) signature found nowhere else (1,283 with term supports found
  nowhere else). How many of these a merged seed set plus today's moves would reach anyway was not
  measured.
- The beam already handles a child with fewer parameters than its parent: it is swept as a
  straggler in the next iteration.

### 11. Derivation support

- **43 fused and merged constructs** were checked against exact laws: 29 derive exactly (SAMPLED:
  3 random laws each, fitted to engine A's exact law below 1e-20); 3 silently return −v; 4 raise a
  BoundsError (an RE step that releases its metabolite); 4 raise a thermodynamic error (the side
  was guessed wrong); 1 finds no reaction cut; and 2 (Theorell–Chance) cannot be written, because
  the DSL silently drops the second metabolite (`src/dsl.jl:744-745`).
- **Causes.** `Step` stores one metabolite and no side, and `_step_sides` guesses the side from
  the species' shapes in five branches. `_compute_alpha` assumes an RE step consumes its
  metabolite (the BoundsError). The numerator's cut search trusts the written direction of its
  steps, so a cut that mixes a reversed and a forward step gives a law that is neither v nor −v.
  The search also builds its central forms only from isomerization endpoints, so it can find no
  cut.
- **A merged complex written product-side** (E(A)+B⇌E(P,Q), then an ordinary release) derives
  exactly today for every RE/SS combination of uni-uni (3), ordered (15) and ping-pong (15)
  mechanisms, and for 161 of the 176 merged random patterns that have a finite rate and pass the
  bottomless-segment check. The other 15 hit the no-cut bug. Another 31 of the 255 random patterns
  have finite rates but are rejected at construction by the bottomless-segment check. Merging
  alone therefore does not require a new `Step`; Theorell–Chance does.
- **Silent name collision.** When both are SS, a fused binding and a pure binding of the same
  metabolite from the same form get the same parameter names (kon_B_EA, koff_B_EA), and
  `fitted_params` lists each name twice without an error.
- **Enumeration blockers for fused steps**: the atom check (for substrate-side fused releases);
  `_flux_carrying_groups` and the allosteric moves, which detect chemistry through `is_iso`;
  `_assert_chemistry_is_iso`, for covalent ping-pong fused steps (bindings such as
  E + A → E(P; residual = A − P) and releases such as EA → F + P); dead-end placement next to fused
  releases; and deduplication, because fused steps have no canonical orientation.
- **Other bugs found**: `_kcat_forward` emits an empty `max()` when no product-free numerator term
  matches a denominator term (in laws that never saturate, and after the −v bug), and residual-form
  parameter names are not quoted in the printed `v =` line.

### 12. Which unmerged isomerizations to keep; the cost of rank

- **Series rule.** Replace an SS isomerization X1⇌X2 by its merge when X1's only other step binds a
  metabolite into X1 and X2's only other step releases one out of X2, both SS and each in a private
  group. The merge has the same family and 2 fewer parameters (PROVED). At depth ≤ 2 this applies
  to 10 bi-bi mechanisms and removes 20 phantoms net. Another 12 mechanisms (22 in all, with 33
  phantoms) have series ends whose groups are shared with abortive steps. Their merge loses a rank,
  so they can only be dropped, and that is lossless at depth ≤ 2 only because both single-flank
  parents and the seed are emitted (PROVED for R4_00441, SAMPLED for 4 more). The orientation
  condition is needed: the free-enzyme isomerization of all-SS iso uni-uni meets the other
  conditions, yet its merge (5/4 → 3/3) loses the P·S term. The fully SS unmerged ordered bi-bi is
  the textbook case of the rule.
- **The broader rule** "every other step at X1 and X2 is SS" loses 2 ping-pong families
  (R4_01650, 8/8; R4_01788, 8/7).
- **Keep every other unmerged SS isomerization**: the merge lowers the rank in 1,450 of the 1,464
  valid bi-bi merges.
- **Rank is cheap.** An exact modular rank prototyped in Julia on the package's own polynomials
  takes 3.5 ms per bi-bi mechanism and matches engine A on 12,534 mechanisms; compile plus
  `fitted_params` takes 0.2–0.6 s.
- **Rank filters lose models.** Rejecting a flip or split child whose rank does not exceed its
  parent's removes 328 bi-bi mechanisms (85% of the phantoms), but a non-identifiable mechanism
  can reach rate laws its identifiable parent cannot (PROVED for R4_00413 against seed R4_00059 and
  R4_00419 against seed R4_00060). In a sample, 8 of 25 removed mechanisms had random laws that
  the union of their kept equal-rank ancestors did not fit (residuals 5e-3 to 0.27; SAMPLED and
  weak, because family fits miss proven containments at such residuals). "Non-identifiable, hence
  useless" is false in general.

### Phantoms that remain

The 616 bi-bi phantom directions fall into six classes (EXACT): zero-flux groups 336, absorbed
transit forms 147, RE isomerizations between forms with the same monomial 60, parallel channels
42, series blocks 28 and a two-segment ping-pong gauge 3. After B and C, about 117 remain at depth
≤ 2. This is computed from the class table and the chain predicate, not measured on implemented
rules: rejecting the 188 zero-flux mechanisms removes 350 phantoms, and class-representative chain
merges in the rest remove 149. By class: absorbed transits whose groups are shared with abortive
steps (about 47), parallel channels (42), RE isomerizations descending from ping-pong seeds 56–60
(about 21), series blocks with shared groups (4) and the two-segment ping-pong gauge (3), which is
the textbook ping-pong's phantom. Allosteric mechanisms were not analyzed. Mechanisms with
competitive inhibitors (R2, R3, R5, R6) were checked alongside R4 for zero flux, duplicate
inhibitors, merges and rank policies (R5), phantom classes (R2, R5) and the modular rank (R5, R6);
the class counts and estimates in this section are for R4 only.

## Decisions (Denis, 2026-09-28)

- **Scope: A + B + C**, as three specs in the order A, B, C.
- **Deferred:**
  - D, MERGE and ELIM as parent-keeping moves anywhere in the search: until an HPC run on C and a
    measurement of what D reaches beyond C.
  - E, charging rank instead of fitted count on the beam's complexity axis: a possible small
    follow-up.
- **What deferring costs**: the unmeasured part of point 10, and the remaining phantoms above,
  which stay charged at their fitted count.

## Plan

### A. Steps with explicit stoichiometry (design agreed; its own spec follows)

- **Data model**: `Step(from_species, to_species, consumed, released, is_equilibrium)`, the only
  constructor; the old four-argument constructor is removed. Kinds: pure binding
  (`bound_metabolite(s)` returns its metabolite), isomerization (`is_iso`), and transformation
  (everything else).
- **Canonical orientation**: the `Step` constructor stores a pure binding with its metabolite
  consumed; the `Mechanism` constructor orients every other step by the existing three tiers, with
  Tier 1 counting free metabolites. Reversing a step swaps its forms and its lists. The sort key
  keeps today's step order, so Haldane pivots and fitted names should not move (CONJECTURE, checked
  by A's regression run).
- **Derivation**: reads the lists and ignores orientation. `_step_sides` is deleted; one RE
  weight rule replaces two; v is the net consumption of the first substrate, summed over the SS
  steps, which replaces the cut search.
- **Names**: pure bindings keep `K_`, `kon_` and `koff_`; every other step is named by form pair
  (`k_EA_to_EQ`, `Kiso_EA_to_EQ`). The `Mechanism` constructor rejects two kinetic groups with the
  same rendered name.
- **DSL**: any number of metabolites on either side of a step.
- **Enumeration code**: mechanical changes only; fused parents stay out of `expand_mechanisms`
  until C.
- **Verification**: identical fitted names and rate equations for all 12,556 exported mechanisms
  and `MECHANISM_TEST_SPECS`, except intended renames of fixtures with fused steps; an orientation
  property test; a fused and Theorell–Chance suite checked against a brute-force mass-action
  oracle; unchanged performance and naming-chokepoint tests.

### B. Exact filters

- **Zero-flux rule**: every SS kinetic group needs a flux-carrying step, checked on every child of
  every move with the structural test of point 1. It replaces `_flux_carrying_groups`. The flip
  move extends a zero-flux set to its minimal flux-carrying supersets instead of dropping it. The
  test is proven for non-allosteric mechanisms; its use on allosteric moves and A/I state graphs
  has not been analyzed.
- **Duplicate-inhibitor rule**: a dead-end copy of a substrate or product must create a new
  complex, judged by composition and by (RE segment, net uptake), on every child (or the split move
  must not split copy groups).
- **To decide**: whether the split move drops a bipartition that leaves an SS part without a
  flux-carrying step or emits its RE-reverted form; required S or P inhibitors in uni-uni, which
  have no identifiable placement in today's unmerged topology; whether to keep tie-only duplicate
  children.

### C. Merged and Theorell–Chance seeds; canonical merges

- **New seeds**: MERGE each seed's chemistry and emit every inclusion-minimal set of groups to
  make SS that gives a valid, flux-carrying, non-degenerate mechanism; ELIM a central complex with
  exactly two steps into a Theorell–Chance step, under the same filter; drop lumping twins and
  duplicates. Variants with more parameters than the seed enter at level 1: the random/random ones
  at 6, and one 6-parameter three-group set for each of the 28 ordered/random seeds.
- **Canonical merges**: replace every qualifying chain (the orientation and private-group
  conditions of points 3–5) by its class-representative merge, whose flank types follow the
  chain's type. The series rule is the SSS case.
- **Moves learn fused steps**: chemistry detection in the flip and allosteric moves, dead-end
  placement, and removal of `_assert_chemistry_is_iso`.
- **`:OnlyA` fused chemistry**: `_onlya_haldane_violation` rejects a hand-written allosteric
  mechanism whose `:OnlyA` chemistry is a fused step, such as a ping-pong with
  `E(A) <--> E(; residual = A - P) + P :: OnlyA`, with a message asking to tag the chemical step
  `:OnlyA`; the same mechanism with the chemistry as an isomerization is accepted, because the
  check drops only `:OnlyA` isomerizations from its cycle graph. The behaviour predates
  sub-project A; C must decide how the check treats fused chemistry once the moves emit it.
- **To decide**: which species represents a merged complex (substrate side or product side); what
  to do with the degenerate ping-pong seeds (drop, keep as non-fitted parents, or replace); how to
  apply the degeneracy check when dead-end copies of substrates or products are present.

## Reproducing

The full track reports, with their proofs and tables, and the skeptics' verdicts are in
`2026-09-28-identifiability-evidence/`. The engines, the Julia export script and the per-track
scripts are session work outside the repository. Everything above can be reproduced from
the descriptions here: build the export with the three catalytic moves, compute ranks exactly at
random points over a large prime field, and test family relations by exact parameter maps or exact
coefficient matching rather than by numerical fits, which often miss proven containments.
