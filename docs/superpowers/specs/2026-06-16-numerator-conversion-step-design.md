# Numerator = net flux through the product-creation cut

ABOUTME: Design for fixing `_compute_numerator` so the rate-equation numerator is
ABOUTME: the net reaction flux, not a partial substrate-binding flux.

## Problem

`_compute_numerator` (src/rate_eq_derivation.jl:418) computes the rate-equation
numerator by tracking ONE substrate metabolite through **SS steps only**, summing
the flux of the steps where that metabolite is bound/released as a *free* species.

This is wrong whenever the tracked metabolite reaches the catalytic core by more
than one route, or via a rapid-equilibrium step:

- **Random-order bi-bi (`m_RE`, G=2):** the substrate binds via two parallel SS
  steps; tracking captures only one branch → numerator is a constant fraction of
  the truth (`0.5977×` in the worked example). Right shape, right equilibrium null,
  wrong magnitude.
- **Single-RE-segment + redundant SS binding (`example_1/2`, G=1):** the tracked
  branch is a Wegscheider-closed cycle → fraction is 0 → **v ≡ 0**.
- The genuine rate-carrying step is the **substrate→product conversion** (iso
  chemistry) step, which carries no free metabolite and is therefore invisible to
  the current tracker.

Confirmed numerically: the package rate disagrees with the exact full-network
steady state (RE limit), while "net flux through the conversion step" matches it to
7 digits for G=1, G=2, and the ordered control. The **denominator is correct**
(expands to the printed King-Altman denominator); only the numerator is wrong.

Scope: the allosteric numerator path reuses `_compute_numerator` (inherits the
bug), and `_kcat_forward` consumes the derived numerator POLY (also corrupted).
Fixing `_compute_numerator` fixes all three.

## The rule

The steady-state velocity is the **King-Altman net cycle flux**: net forward-reaction
flux across a directed cut of SS steps that the catalytic cycle crosses exactly once
per turnover. A correct numerator must

- **sum parallel routes** — random-order binding: `{bind S1, bind S2}` together = v,
  either alone = a fraction of v; and
- **count series routes once** — ordered all-SS bi-bi: `bind A` and `bind B` are each
  = v, so summing them double-counts.

The current code (track one metabolite through SS steps) gets series right, parallel
wrong. A "net flux through the conversion step" rule gets random-order right only when
chemistry is SS; it fails when the catalytic isomerization is at rapid equilibrium
(a physically real case: binding/release SS, `E(S1,S2) ⇌ E(P1,P2)` RE). Neither is
general. Only a true directed-cut / net-cycle-flux numerator is.

The building blocks already exist and are correct: per-SS-step `flux(s) = rf·D[g1] −
rr·D[g2]`, the segment cofactors `D[g]`, the α-weights, and the denominator. Because
segments are RE-connected components, **every inter-segment edge is SS**, so any
segment-boundary cut is computable. What is missing is choosing the cut and the
forward-reaction orientation of each crossing edge so that parallel edges sum and
series edges count once.

### Orientation of an SS step (forward reaction)

- substrate-binding step → forward = canonical from→to (substrate consumed)
- product-release step (stored canonically as product-binding `E+P→EP`) → forward =
  canonical to→from (product released)
- chemistry/iso step → forward = canonical from→to (physical-forward, already canon)

### The cut algorithm: reaction-progress (ρ) boundaries

A catalytic turnover is a sequence of elementary forward events
(bind substrate, isomerize, release product) returning free E → free E. Assign each
enzyme form a **reaction-progress level** ρ(f) = number of forward events from free E
to f. Then the forward-oriented form graph is *graded*: every forward edge goes from
level k to level k+1, and one turnover advances ρ from 0 to N (= total events) — so
**every ρ-boundary k|k+1 is crossed exactly once per turnover.**

Therefore, for ANY boundary k|k+1, the net forward flux across it = v. We need a
boundary whose crossing steps are **all SS** (so every flux is computable; segments
are RE components, so this is automatically satisfiable when such a boundary exists):

Implementation note: this is computed WITHOUT a global ρ assignment. Equivalent and
simpler (avoids the gradedness assumption) is the **step-type / parallel-RE** form:

```
A reaction CUT is a set of SS steps that the turnover crosses exactly once. Steps are
PARALLEL (belong to one cut, summed) iff they share the same to-form (converging binds,
ES1→ES1S2 & ES2→ES1S2) or the same from-form (diverging releases/chem, EP1P2→EP1 &
EP1P2→EP2). SERIES steps (to-form of one = from-form of next, E→EA→EAB) are NEVER summed.
An SS step "has a parallel RE step" iff one of its parallel routes is RE → that cut is
unusable (summing only SS routes would undercount).

Candidate cuts (each is a complete cut crossed once per turnover ⇒ flux sum = v):
  (a) SS iso/chemistry steps that convert ONE chosen substrate→product, no parallel RE
  (b) SS substrate-binding steps for one progress level, no parallel RE
  (c) SS product-release steps for one progress level, no parallel RE
NUM = Σ flux(s) over the qualifying cut with the FEWEST steps (simplest numerator —
they all give the identical equation); tie-break toward the chemistry cut (a)
(standard Vmax·(S−P/Keq) form). No stoichiometric normalisation — the cut is crossed
once per turnover, so the sum IS v·DEN/E_total.
If NO cut qualifies (every progress level has a parallel RE) → raise: equivalent to "a
complete all-RE catalytic cycle exists ⇒ Vmax→∞ ⇒ no finite RE rate" (proven below).
```

The ρ-levels above are the conceptual model (each cut = one ρ-boundary); the step-type
form computes the same cuts via local shared-form checks, with no graded-graph
requirement — robust to interleaved-release / Theorell-Chance mechanisms.

### Construction guard (separate, small): forbid a reaction being both RE and SS

`@enzyme_mechanism` currently accepts the same physical reaction twice with different
arrows, e.g. `E + S <--> E(S)` AND `E + S ⇌ E(S)` — representable only because `Step`
carries `is_equilibrium` in its identity, so the duplicates aren't `==` and both survive
`unique!`. Such a mechanism is nonsensical (a step can't be slow and fast). The numerator
algorithm already raises on it (all-RE cycle ⇒ rule (4)), but it should be rejected at
construction: the `Mechanism`/`AllostericMechanism` ctor must reject any two steps with the
same `(from_species, to_species, bound_metabolite)` differing only in `is_equilibrium`.
This is a small standalone validation, distinct from the numerator fix; broader
"mirror via different bystanders" handling stays with #1.

Why this is exactly right: parallel routes (random order) live on the SAME boundary →
summed. Series routes (ordered) live on DIFFERENT boundaries → never summed (we use
one boundary). Redundant internal-SS binding steps (the `example_1/2` degeneracy) sit
on a *mixed* boundary (one route SS, its parallel route RE) which step 4 skips,
falling through to the all-SS chemistry boundary — yielding the correct non-zero rate.

### Per-class walkthrough (ρ-levels; **bold** = chosen all-SS boundary)

| Mechanism | ρ-levels (forms) | boundaries & kinds | NUM | vs exact |
|---|---|---|---|---|
| **Uni-Uni RE** (bind RE, chem SS, rel RE) | E0, ES1, EP2 | 0\|1 bind RE; **1\|2 chem SS**; 2\|0 rel RE | flux(chem) | ✓ (= today) |
| **Ordered Bi-Bi all-SS** | E0, EA1, EAB2, EPQ3… | **0\|1 bindA SS**; 1\|2 bindB SS; 2\|3 chem SS… | flux(bindA) | ✓ (= today's substrate-A track) |
| **RE Ordered Bi-Bi** (bind RE, chem+rel SS) | E0, EA1, EAB2, EPQ3… | 0\|1,1\|2 bind RE; **2\|3 chem SS**; rel SS | flux(chem) | ✓ value (raw string may move from rel→chem) |
| **Ping-pong Bi-Bi** (all SS, single cycle) | E0, EA1, F2, FB3 | **0\|1 bindA SS**; 1\|2 chem1 SS; 3\|0 chem2 SS | flux(bindA) | ✓ (= today) |
| **`m_RE` random, chem SS** (G=2) | E0, ES1/ES2 =1, ES1S2=2, EP=3 | 0\|1 bind RE; **1\|2 {bindS2,bindS1} SS**; 2\|3 chem SS | flux(3)+flux(4) | ✓ (today: only one → 0.598×) |
| **Denise RE-chem** (G=2, iso RE) | E0, ES·=1, ES1S2=2, EP1P2=3, EP·=4 | 0\|1 bind RE; **1\|2 {3,4} SS**; 2\|3 chem RE; 3\|4 {6,7} SS; 4\|0 rel RE | flux(3)+flux(4) | ✓ (today: 0.478×) |
| **`example_1/2` degenerate** (G=1) | E0, ES1/ES2=1, ES1S2=2, EP=3 | 0\|1 **mixed** (S1 SS, S2 RE) → skip; 1\|2 bind RE; **2\|3 chem SS** | flux(chem) | ✓ non-zero (today: 0) |
| **mixed-binding, all-SS release** (G=1) | E0, ES·=1, ES1S2=2, EP1P2=3, EP·=4 | 1\|2 **mixed** (S2 SS, S1 RE) → skip; **3\|4 {relP2,relP1} SS** | flux(6)+flux(7) | ✓ (today: one release → 0.48×; denom faithful, "no P1·P2" correct) |
| **no all-SS boundary** (all-RE catalytic cycle) | — | every boundary has an RE crossing | **raise** | rate diverges (Vmax→∞); today silently emits v≡0 |

### "No all-SS boundary → raise" is correct, not too aggressive (resolved)

**No all-SS ρ-boundary ⟺ a complete all-RE catalytic cycle exists ⟺ Vmax → ∞.**
Proof: if every boundary has an RE crossing, choosing the RE edge at each level 0→N
gives an all-RE path around one full turnover — an all-RE catalytic cycle. In the
RE limit that cycle carries flux ∝ M (the RE rate), so turnover diverges and the
mechanism has no finite rapid-equilibrium rate. Verified numerically: a mechanism
with binding boundary mixed AND release boundary mixed AND RE chemistry gives
`v/M = 0.071055` constant from M=1e4..1e8 (v = 0.071·M → ∞). Today's code derives it
and prints `v ≡ 0` (numerator cancels after Wegscheider). So `raise` is strictly
correct and strictly better than the current silent-wrong-answer. The ρ-cut raises
intrinsically (it has no boundary to sum); this is part of the numerator fix, distinct
from #1's enumeration pruning of these all-RE-cycle / redundant-SS mechanisms.

### Why existing oracles are still expected to be preserved

All oracle specs are ordered / single-cycle, where the net flux through any series
cut equals J. A correct directed-cut numerator yields J for them — same value, and
(for single-cycle) the same numerator polynomial — so they should stay green. The
specs that expose the bug (random-order, RE-chemistry) have no analytical oracle,
which is why it went uncaught; the exact-steady-state solve is their oracle.

### Faithful-but-degenerate mechanisms are out of scope (resolved)

A uni-uni with SS substrate binding and RE chemistry+release
(`E+S <--> E(S)`, `E(S) ⇌ E(P)`, `E(P) ⇌ E+P`) derives to
`v = (kon·S − koff·P/(K_P·Kiso)) / (1 + P/(K_P·Kiso) + P/K_P)` — **no S in the
denominator, no substrate saturation**. This was verified to equal the exact
full-network steady state in the RE limit to 6 digits: it is **faithful, not a bug**.
The mechanism is degenerate (the lone SS step is a redundant binding internal to one
RE segment; chemistry+release drain `E(S)` infinitely fast so it never accumulates →
Vmax → ∞). The ρ-cut algorithm reproduces this faithful rate (boundary 0|1 is the only
all-SS boundary → flux of the binding step). **The derivation must NOT special-case or
raise on these** — they are removed by the enumeration guard (#1, separate work).

This sharpens the taxonomy — two different problems, two different fixes:
- **Genuine numerator bug** (`m_RE` random chem-SS; RE-chemistry bi-bi; `example_1/2`):
  current code ≠ exact. Cause: parallel SS routes with only one tracked, or the
  numerator locking onto a dead binding cycle. → **ρ-cut numerator fixes these.**
- **Faithful-but-degenerate** (the uni-uni above): current code = exact, but the
  mechanism has no saturation. → **#1 prunes; derivation already correct.**

### Oracle safety of the ρ-cut (gradedness + side-branches)

- **Gradedness** (well-defined ρ) holds for every spec that has an analytical oracle:
  they are ordered/single-route (trivially graded) or RE-binding-random with a single
  SS chemistry step (graded: both binding orders reach the central complex in equal
  events). The only gradedness risk — random-order with product release *interleaved*
  before all substrates bind — has **no analytical oracle**; there the exact-steady-
  state solve is the oracle, and a genuinely non-graded mechanism raises.
- **Side-branches** (dead-end inhibitors, bystanders): inhibitor/regulator binding is
  NOT a forward event, so those forms stay intra-level (ρ unchanged) and their RE
  dead-end steps never define a boundary. The inhibitor specs all reduce to a single
  SS chemistry boundary → correct.

### Edge cases / open risks (validate in TDD)

1. **Raw-string regression.** Numeric rate unchanged for ordered oracles, but the raw
   numerator POLY may be written via different rate constants. If the flat-string /
   Expr-shape snapshots in test/test_rate_eq_derivation.jl move, re-baseline with
   justification (rate is still correct).
2. **Laurent / `_reduce_conc_lowest_terms`.** Cut sums add/subtract whole `flux(s)`
   POLYs; result stays a POLY, no new division. Confirm on multi-product specs.
3. **Performance.** Change is in the @generated derivation (compile time), not the
   emitted body — runtime unaffected — but re-run test_rate_equation_performance.

## TDD plan

1. Add the three exact-steady-state cases as oracles (ground truth = exact
   King-Altman solve of the full network in the RE limit, already prototyped in
   scratch_g2_compare.jl / scratch_exact_steadystate.jl): `m_RE` (G=2, factor 0.5977),
   `example_1`/`example_2` (G=1, v≡0). Confirm RED against current code.
2. Implement the δ rule in `_compute_numerator`; add `residual`-aware `π`.
3. GREEN on the new cases; run the FULL suite — every Segel/RE/ping-pong/allosteric
   analytical oracle and every kcat oracle must stay green; JET/Aqua clean;
   chokepoint + perf tests green.
4. If raw-string snapshots move, review each and re-baseline with justification.
5. Drop the enumeration redundancy guard (fix #1) into a separate follow-up — it is
   independent and only prunes non-identifiable candidates; it is not needed for
   correctness once the numerator is fixed.
