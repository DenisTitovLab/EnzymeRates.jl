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
- product-release step → forward depends on storage form: a product-*binding* storage
  (`E+P→EP`, product on `m_lhs`) stores the reverse reaction → forward = canonical
  to→from; an SS-*dissociation* storage (`EA→E+P`, product on `m_rhs`) already stores
  the release → forward = canonical from→to. (The endpoints `(ff,ft)` used for
  central-species cuts AND the flux sign must both respect this.)
- chemistry/iso step → forward = canonical from→to. An iso step counts as chemistry
  ONLY if it changes bound substrate/product content or the covalent `Residual`; a pure
  conformational iso (identical bound metabolites and residual on both endpoints) is not
  a reaction-progress step and generates no cut (its endpoints are not central species).

### The cut algorithm: metabolite cuts + central-species cuts

A **reaction cut** is a set of SS steps that one turnover crosses exactly once; the
numerator is the oriented-flux sum over one complete cut. Each conserved per-turnover
"event" is a candidate cut (its steps all carry that event once, so the net flux sums
to v). A candidate is **usable** iff *all* its forward reaction steps are SS (any RE
step in the event would make the SS-only sum an undercount).

**Exclude dead-end steps** when collecting a candidate's steps (from both the all-SS
check and the flux sum): a step touching a *substrate–product mixed complex* — a form
with **both** a bound `Substrate` and a bound `Product` — is off the catalytic path
(reached only by product rebinding to a substrate complex, or substrate binding to a
product complex) and carries zero net flux, so a dead-end RE step must not disqualify an
otherwise-all-SS cut. **Regulator-bound forms are NOT mixed** (a regulator is neither a
substrate nor a product), so steps binding S / releasing P on a regulator-bound enzyme
are kept — they are genuine parallel catalytic routes, and their RE-ness *should* count.
The central complexes stay valid (`E(S1,S2)` is pure-substrate, `E(P1,P2)` pure-product);
ping-pong intermediates carry a `Residual`, not bound substrate+product, so are not flagged.

**Why two kinds of candidate are needed.** A single metabolite's consumption can span
several stages (`bind S1 to E` and `bind S1 to E(S2)`) — a *metabolite* cut. A single
stage can span several metabolites (two different products leaving `E(P1,P2)`) — a
*species* cut. Neither subsumes the other (verified by counterexample), so we generate
both. Restricting species cuts to **iso-step-endpoint** forms keeps them valid (those
forms lie on every turnover's path) and avoids the *branch trap* (e.g. "produce
`E(S1)` = {one SS step}" is all-SS but wrong, because `E(S1)` is a bypassed side branch
— `E(S1)` is not an iso endpoint, so it is never offered as a cut).

Candidate cuts (forward-oriented; usable iff every step is SS):

```
Metabolite cuts
  per substrate S:  {steps that BIND S}            or {iso steps that CONVERT S → residual/product}
  per product   P:  {steps that RELEASE P}         or {iso steps that PRODUCE P from substrate/residual}
Central-species cuts   (X = from- or to-form of any ISO step; the chemistry endpoints)
  per such X:       {steps that PRODUCE X}          or {steps that CONSUME X}
```

(Forward orientation: bind/chem = canonical from→to; release = canonical to→from, since
release is stored as product binding.) Central species are the substrate/residual-side
and product/residual-side forms of iso steps — this covers ping-pong covalent
intermediates, where there is no "all-substrate-bound" form but the iso step's
endpoints still sit on every turnover.

```
NUM = Σ oriented flux over the chosen usable candidate. Selection order:
  (1) prefer a METABOLITE cut over a central-species cut, then
  (2) fewest steps, then
  (3) a chemistry/iso cut (standard Vmax·(S−P/Keq) form), then
  (4) sorted step indices (deterministic across precompile sessions).
A metabolite cut is ALWAYS complete (its metabolite is consumed/produced once per
turnover, summed over all parallel routes). A central-species cut can be INCOMPLETE
when parallel routes pass through different iso endpoints (e.g. a free `E(P)` and a
regulator-bound `E(P,R)`) — "produce `E(P)`" would then capture only the free route —
so central cuts are the fallback, used only when chemistry is RE and no metabolite cut
is all-SS. No stoichiometric normalisation — the cut is crossed once per turnover, so
the sum IS v·DEN/E_total. Orientation sign: bind/chem keep the canonical from→to flux;
product release flips ONLY when stored as product-binding (`E+P→EP`, product on
`m_lhs`); a release stored as SS-dissociation (`EA→E+P`, product on `m_rhs`) is already
forward, no flip.
If NO candidate is usable → raise: equivalent to "a complete all-RE catalytic cycle
exists ⇒ Vmax→∞ ⇒ no finite RE rate" (proven below). No free-E reference anywhere, so
the rule is robust to multiple free-enzyme conformations.
```

Worked cases: `example_1` → convert-S1 / consume-`E(S1,S2)` = {chem} (metabolite-S1
*binding* is mixed RE/SS → skipped) → correct non-zero; divergence (S1 binds SS to both
`E` and `E(S2)`) → bind-S1 = {step1, step4}; the RE-chemistry mechanism (no metabolite
all-SS) → consume-`E(P1,P2)` = {step6, step7}; ping-pong → bind-A or its iso-endpoint
central species. Completeness caveat: a mechanism whose only SS bottleneck is neither a
metabolite's full bind/release/convert set nor an iso-endpoint's full produce/consume
set is not provably covered; it would `raise` rather than mis-derive (none constructed).

### Construction guard (separate, small): forbid a reaction being both RE and SS

`@enzyme_mechanism` currently accepts the same physical reaction twice with different
arrows, e.g. `E + S <--> E(S)` AND `E + S ⇌ E(S)` — representable only because `Step`
carries `is_equilibrium` in its identity, so the duplicates aren't `==` and both survive
`unique!`. Such a mechanism is nonsensical (a step can't be slow and fast). The numerator
algorithm already raises on it (no usable cut ⇒ all-RE cycle), but it should be rejected at
construction: the `Mechanism`/`AllostericMechanism` ctor must reject any two steps with the
same `(from_species, to_species, bound_metabolite)` differing only in `is_equilibrium`.
This is a small standalone validation, distinct from the numerator fix; broader
"mirror via different bystanders" handling stays with #1.

Why this is exactly right: each candidate is a per-turnover-conserved event, so its
SS steps' fluxes sum to v (parallel routes within the event are summed; the event is
crossed once, so series steps live in *different* candidates and are never summed). A
redundant internal-SS binding (the `example_1/2` degeneracy) has an RE step in its
metabolite-binding event (`bind S1` includes an RE route) → that candidate is skipped,
falling through to the convert-S1 / consume-`E(S1,S2)` candidate → correct non-zero.

### Per-class walkthrough (usable candidate cut in **bold**)

| Mechanism | usable candidate cuts | chosen | NUM | vs exact |
|---|---|---|---|---|
| **Uni-Uni RE** (bind RE, chem SS, rel RE) | **convert-S = {chem}** | {chem} | flux(chem) | ✓ (= today) |
| **Ordered Bi-Bi all-SS** | bind-A, bind-B, **convert = {chem}**, release | {chem} — tie-break | flux(chem) | ✓ value (string moves bindA→chem) |
| **RE Ordered Bi-Bi** (bind RE, chem+rel SS) | **convert = {chem}**, release-P | {chem} | flux(chem) | ✓ value (string may move) |
| **Ping-pong Bi-Bi** (concerted release) | **bind-A**, bind-B, convert(iso) | bind-A (fewest, no separate chem) | flux(bindA) | ✓ (= today) |
| **`m_RE` random, chem SS** (G=2) | bind-S1/S2 (2 steps), **convert = {chem}** | {chem} — 1 step | flux(chem) | ✓ (today one branch → 0.60×) |
| **Denise RE-chem** (G=2, iso RE) | bind/release all mixed; convert(iso) RE; **consume-`E(P1,P2)` = {6,7}** | {6,7} | flux6+flux7 | ✓ (today 0.48×) |
| **`example_1/2` degenerate** (G=1) | bind-S1 mixed → skip; **convert-S1 = {chem}** | {chem} | flux(chem) | ✓ non-zero (today 0) |
| **mixed-binding, all-SS release** (G=1) | bind mixed; convert RE; **consume-`E(P1,P2)` = {6,7}** | {6,7} | flux6+flux7 | ✓ (today 0.48×; denom faithful) |
| **divergence** (S1 binds SS to E & ES2) | **bind-S1 = {step1,step4}** (convert/consume all RE) | {step1,step4} | flux1+flux4 | ✓ |
| **no usable candidate** (all-RE catalytic cycle) | — | — | **raise** | diverges (Vmax→∞); today v≡0 |

### "No usable candidate → raise" is correct, not too aggressive (resolved)

**No usable candidate ⟺ a complete all-RE catalytic cycle exists ⟺ Vmax → ∞.**
Proof: if every candidate cut contains an RE step, then every conserved event has an
RE route; stringing those RE routes together gives an all-RE
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
Vmax → ∞). The reaction-cut algorithm reproduces this faithful rate: its only usable
candidate is bind-S = {the lone SS binding step} (no RE step in that event), so NUM =
flux of that step. **The derivation must NOT special-case or raise on these** — they
are removed by the enumeration guard (#1, separate work).

This sharpens the taxonomy — two different problems, two different fixes:
- **Genuine numerator bug** (`m_RE` random chem-SS; RE-chemistry bi-bi; `example_1/2`):
  current code ≠ exact. Cause: parallel SS routes with only one tracked, or the
  numerator locking onto a dead binding cycle. → **reaction-cut numerator fixes these.**
- **Faithful-but-degenerate** (the uni-uni above): current code = exact, but the
  mechanism has no saturation. → **#1 prunes; derivation already correct.**

### Oracle safety (no free-E split; side-branches)

- **No free-enzyme reference.** Candidates are anchored at metabolites and iso-step
  endpoints, never at "free E about to bind vs having released", so the rule is robust
  to mechanisms with multiple free-enzyme conformations (`E`, `E*`).
- **Branch trap avoided.** Species cuts are restricted to iso-step-endpoint forms,
  which lie on every turnover's path; a bypassed side branch like `E(S1)` is never an
  iso endpoint, so "produce `E(S1)` = {one SS step}" is never offered as a candidate.
- **Side-branches** (dead-end inhibitors, bystanders): inhibitor/regulator binding is
  not a reaction step (`_reaction_step` returns `nothing`), so those forms never form a
  cut. The inhibitor specs reduce to a single SS chemistry/convert cut → correct.

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
