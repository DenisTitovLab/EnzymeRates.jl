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

### The cut algorithm: SS steps grouped by mirror, RE-mirror check

A **reaction cut** is a set of SS steps that one turnover crosses exactly once. The
numerator is the oriented-flux sum over one complete cut. Two reaction steps are
**mirrors** (belong to the same cut) iff they have the same `(level, type)`, where
`type ∈ {:bind, :chem, :release}` and `level` is a *local* property of the step's
reaction-from form:

```
level(bind step)    = #substrates committed in its from-form
level(chem step)    = #substrates committed in its from-form
level(release step) = #products  committed in its from-form
  #substrates committed(form) = count(Substrate in bound(form)) + length(residual.added)
  #products  committed(form)  = count(Product  in bound(form)) + length(residual.subtracted)
```

(Reaction-from form: bind/chem use canonical `from`; release uses canonical `to` —
the fuller complex — because release is stored as product *binding*.)

`level` is just "how far into the reaction this step's starting form is" without any
graph traversal. Mirrors are the parallel alternatives at one stage: random-order
first binding `E+S1` / `E+S2` are both `(0, :bind)`; second binding `ES1+S2` / `ES2+S1`
are both `(1, :bind)`; the two product releases from `E(P1,P2)` are both
`(2, :release)`. The `residual` terms keep ping-pong's second substrate (its from-form
carries the first substrate covalently in `residual.added`) from colliding with the
first.

```
For each reaction step compute (level, type); record RE keys and SS keys separately.
A cut (level, type) is USABLE iff it has SS steps and NO RE step shares its key
(no parallel rapid-equilibrium "mirror" route — else the SS-only sum undercounts).
NUM = Σ oriented flux over the USABLE cut with the FEWEST steps (simplest numerator;
all usable cuts give the identical equation); tie-break toward chemistry (standard
Vmax·(S−P/Keq) form). No stoichiometric normalisation — the cut is crossed once per
turnover, so the sum IS v·DEN/E_total.
If NO cut is usable (every reaction step has an RE mirror) → raise: proven below to be
equivalent to "a complete all-RE catalytic cycle exists ⇒ Vmax→∞ ⇒ no finite RE rate".
```

This needs no reaction-progress BFS and makes no gradedness assumption — `level` is a
direct metabolite count on each form, so it is robust to interleaved-release /
Theorell-Chance / ping-pong mechanisms. `example_1`'s redundant `E+S1<-->E(S1)` (SS)
is killed because the RE `E+S2⇌E(S2)` shares its `(0, :bind)` key, so the algorithm
falls through to the chemistry cut and yields the correct non-zero rate.

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

Why this is exactly right: parallel routes (random order) share a mirror key → summed.
Series routes (ordered) have different keys (different committed-metabolite counts) →
never summed. A redundant internal-SS binding (the `example_1/2` degeneracy) shares its
key with a parallel RE step (same stage, other metabolite) → that cut is skipped,
falling through to the chemistry cut → correct non-zero rate.

### Per-class walkthrough (`(level, type)` mirror keys; **bold** = chosen cut)

| Mechanism | reaction-step keys (SS unless noted RE) | chosen cut | NUM | vs exact |
|---|---|---|---|---|
| **Uni-Uni RE** (bind RE, chem SS, rel RE) | (0,bind)RE; **(1,chem)**; (1,rel)RE | (1,chem) | flux(chem) | ✓ (= today) |
| **Ordered Bi-Bi all-SS** | (0,bind); (1,bind); **(2,chem)**; (2,rel) | (2,chem) — tie-break | flux(chem) | ✓ value (string moves bindA→chem) |
| **RE Ordered Bi-Bi** (bind RE, chem+rel SS) | (0,bind)RE,(1,bind)RE; **(2,chem)**; (2,rel) | (2,chem) | flux(chem) | ✓ value (string may move rel→chem) |
| **Ping-pong Bi-Bi** (all SS, concerted release) | **(0,bind)**; (0,rel); (1,bind); (1,rel) | (0,bind) — no chem, bind<rel | flux(bindA) | ✓ (= today) |
| **`m_RE` random, chem SS** (G=2) | (0,bind)RE×2; (1,bind)×2; **(2,chem)** | (2,chem) — 1 step < 2 | flux(chem) | ✓ (today: one branch → 0.60×) |
| **Denise RE-chem** (G=2, iso RE) | (1,bind): S2 SS+S1 RE → skip; (2,chem)RE; **(2,rel)×2** | (2,rel) {6,7} | flux(6)+flux(7) | ✓ (today: 0.48×) |
| **`example_1/2` degenerate** (G=1) | (0,bind): S1 SS+S2 RE → skip; (1,bind)RE; **(2,chem)** | (2,chem) | flux(chem) | ✓ non-zero (today: 0) |
| **mixed-binding, all-SS release** (G=1) | (1,bind): S2 SS+S1 RE → skip; (2,chem)RE; **(2,rel)×2** | (2,rel) {6,7} | flux(6)+flux(7) | ✓ (today 0.48×; denom faithful) |
| **no usable cut** (all-RE catalytic cycle) | every key has an RE mirror | — | **raise** | diverges (Vmax→∞); today v≡0 |

### "No usable cut → raise" is correct, not too aggressive (resolved)

**No usable cut ⟺ a complete all-RE catalytic cycle exists ⟺ Vmax → ∞.**
Proof: if every reaction step has an RE mirror (same `(level, type)`), then at every
reaction stage some route is RE; choosing the RE route at each stage gives an all-RE
path around one full turnover — an all-RE catalytic cycle. In the RE limit that cycle
carries flux ∝ M (the RE rate), so turnover diverges and the mechanism has no finite
rapid-equilibrium rate. Verified numerically: a mechanism with binding stage mixed AND
release stage mixed AND RE chemistry gives `v/M = 0.071055` constant from M=1e4..1e8
(v = 0.071·M → ∞). Today's code derives it and prints `v ≡ 0` (numerator cancels after
Wegscheider). So `raise` is strictly correct and strictly better than the current
silent-wrong-answer; the algorithm raises intrinsically (no usable cut to sum). This is
part of the numerator fix, distinct from #1's enumeration pruning of these
all-RE-cycle / redundant-SS mechanisms.

### Why existing oracles are still expected to be preserved

All oracle specs are ordered / single-cycle, where the net flux through any series cut
equals J. A correct reaction-cut numerator yields J for them — same value, and (for
single-cycle) the same numerator polynomial for the chosen cut — so they should stay
green. The specs that expose the bug (mixed RE/SS binding, RE-chemistry) are not in the
suite; the exact-steady-state solve is their oracle.

### Faithful-but-degenerate mechanisms are out of scope (resolved)

A uni-uni with SS substrate binding and RE chemistry+release
(`E+S <--> E(S)`, `E(S) ⇌ E(P)`, `E(P) ⇌ E+P`) derives to
`v = (kon·S − koff·P/(K_P·Kiso)) / (1 + P/(K_P·Kiso) + P/K_P)` — **no S in the
denominator, no substrate saturation**. This was verified to equal the exact
full-network steady state in the RE limit to 6 digits: it is **faithful, not a bug**.
The mechanism is degenerate (the lone SS step is a redundant binding internal to one
RE segment; chemistry+release drain `E(S)` infinitely fast so it never accumulates →
Vmax → ∞). The reaction-cut algorithm reproduces this faithful rate: its only usable cut
is the `(0,bind)` binding step (no parallel RE mirror), so NUM = flux of that step.
**The derivation must NOT special-case or raise on these** — they are removed by the
enumeration guard (#1, separate work).

This sharpens the taxonomy — two different problems, two different fixes:
- **Genuine numerator bug** (`m_RE` random chem-SS; RE-chemistry bi-bi; `example_1/2`):
  current code ≠ exact. Cause: parallel SS routes with only one tracked, or the
  numerator locking onto a dead binding cycle. → **reaction-cut numerator fixes these.**
- **Faithful-but-degenerate** (the uni-uni above): current code = exact, but the
  mechanism has no saturation. → **#1 prunes; derivation already correct.**

### Oracle safety (no graph traversal; side-branches)

- **No gradedness assumption.** `level` is a direct metabolite count on each form, so
  the rule needs no reaction-progress BFS and is robust to interleaved-release /
  Theorell-Chance / ping-pong mechanisms. The committed-count (with the `residual`
  term) separates series steps for every single-cycle mechanism.
- **Side-branches** (dead-end inhibitors, bystanders): inhibitor/regulator binding is
  not a reaction step (`_reaction_step` returns `nothing`), so those forms never form a
  cut. The inhibitor specs reduce to a single SS chemistry cut → correct.

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
