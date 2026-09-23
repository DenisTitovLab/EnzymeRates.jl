# Iso mechanisms and single-subunit kinetic cooperativity: findings

Date: 2026-09-23. Status: research findings; no design decided yet.

This note records what a research pass established before any design work on two goals:

1. **Iso mechanisms in enumeration.** Chemistry turns the substrate complex of conformation
   `E` into the product complex of a second conformation `Eprime`; products leave `Eprime`;
   free `Eprime` returns to `E` through a slow steady-state (SS) step.
2. **Kinetic cooperativity in single-subunit enzymes.** Two slowly interconverting
   conformations both catalyze, so competing routes stay out of equilibrium and the rate law
   gains power terms, ideally with fewer parameters than fully SS random bi-bi.

The evidence comes from three independent sources: an exact sympy engine (Cha segments plus
the Markov-chain tree theorem, checked against brute-force mass action), the package's own
derivation checked against a BigFloat mass-action ground truth, and two independent
re-derivations of every theoretical claim. Hand-written PGK mechanisms used Keq = 760.

## Summary

- A slow conformational step that lies on no catalytic cycle carries zero net flux at steady
  state. Its two conformations then sit at their equilibrium ratio, and the mechanism is
  equivalent to one in which that step is rapid equilibrium (RE).
- With RE binding, conformations stay out of equilibrium only when **chemistry itself crosses
  conformations**, for example `E(BPG, ADP) <--> Eprime(ThreePG, ATP)`.
- The "every step duplicated in `E` and `Eprime`" model therefore equals RE random bi-bi: it
  has 16 fitted parameters and 7 identifiable ones, the same 7 as RE random bi-bi.
- Cleland iso mechanisms stay hyperbolic in every metabolite; they add substrate × product
  cross terms only. They serve goal 1, not goal 2.
- The symmetric **flip-flop** (chemistry in `E` lands in `Eprime` and chemistry in `Eprime`
  lands in `E`) is non-hyperbolic in all four metabolites of a bi-bi reaction.
- In the PGK data, most of the RE→SS random bi-bi loss gap comes from the ADP
  product-inhibition pattern, not from curvature. The genuine curvature is apparent negative
  cooperativity (Hill h ≈ 0.5–0.7), mainly in reverse-direction ATP titrations.

## 1. Which conformational mechanisms collapse

**Theorem (proved; exhaustively checked).** Contract each RE segment (a connected component of
RE steps) to one node. Keep each SS step that joins two different segments as an edge labelled
with its reaction extent: 1 for chemistry, 0 for binding or isomerization. An SS step whose
biconnected block contains only zero-extent cycles (a bridge qualifies) carries zero net flux
at steady state for every positive concentration and every admissible parameter set. Making
that step RE leaves the rate law unchanged.

Proof sketch: in the gauge `x_X = ψ·w_X`, the flux through each step of such a block is
`g·Δψ` with `g > 0`. Steady-state flows restricted to the block form a circulation, so
`Σ g·Δψ² = 0` and every `Δψ` vanishes.

The rate law was identical before and after this collapse for all 22,560 two-conformation
uni-uni mechanisms built from an 11-step pool and for 9,000 random two-conformation bi-bi
mechanisms.

Consequences:

- **Duplicated model.** When binding is RE inside each conformation and chemistry stays
  inside a conformation, every set of pure isomerizations `E(X) <--> Eprime(X)` carries zero
  flux. That holds whether the link sits at free enzyme only or at every bound form. The
  duplicated uni-uni variants reduce to Michaelis–Menten (3 identifiable parameters of 8–10).
  The duplicated PGK bi-bi has identifiable rank 7 of 16 (links at free enzyme only) or 7 of
  22 (links everywhere), equal to RE random bi-bi's 7 of 7.
- **Power terms that cancel.** Adding SS links at bound forms, such as
  `Eprime(S) <--> E(S)`, puts S² terms (ADP²·BPG² for PGK) into the raw equation. The
  Wegscheider constraints on the link cycles force a common factor into numerator and
  denominator: the summed rate of the conformation-changing steps. The package never cancels
  polynomial factors, so the printed equation looks cooperative while `x/v` stays linear in
  `x` to about 1e-16. The structural `_hyperbolic_catalysis` predicate calls these models
  non-hyperbolic.
- **RE links.** With RE binding, one RE binding or isomerization step between the
  conformations merges them into one RE segment, and the rate law becomes a Cha RE law. An RE
  *chemistry* step between conformations is the one exception: it can stay non-hyperbolic when
  product is present.
- **SS binding changes the picture.** With SS binding, isomerizations at bound forms start to
  carry flux (for example, the duplicated uni-uni with links at `ES`: rank 3 with RE binding,
  8 with SS binding). A bridge still carries none.

## 2. What produces genuine curvature

With RE binding in both conformations, a metabolite `m` gets a degree-2 denominator only when
a conformation that holds an `m`-bound form is entered by a cross step whose rate depends on
`m`. A 2:2 (cooperative-type) law also needs numerator degree 2: either internal chemistry in a
conformation entered by an `m`-dependent step, or two cross chemistry steps of opposite
direction. One catalytic route alone gives at most 1:2, the substrate-inhibition shape.

The mnemonic model as reviews describe it (Ricard, Meunier & Buc 1974; the primary text was
not checked) ties the conformational change to binding, not chemistry. With RE binding it
collapses to Michaelis–Menten; it needs SS binding to curve. Denis's variant ties the change
to chemistry and curves with RE binding.

Uni-uni, RE binding (reduced degrees after all constraints and cancellation):

| Mechanism | Degree in S (num:den) | Degree in P | Identifiable / fitted |
|---|---|---|---|
| Iso: `E(S) <--> Eprime(P)`, `Eprime <--> E` | 1:1 | 1:1 | 4 / 5 |
| Duplicated (any links, chemistry inside each conformation) | 1:1 | 1:1 | 3 / 8–10 |
| Mnemonic, change tied to chemistry | 2:2 | 1:1 | 6 / 7 |
| Mnemonic, change tied to binding (classic) | 1:1 | 1:1 | 3 / 7 |
| Flip-flop: `E(S) <--> Eprime(P)`, `Eprime(S) <--> E(P)`, `Eprime <--> E` | 2:2 | 2:2 | 7 / 8 |
| Flip-flop without `Eprime <--> E` | 1:1 at P = 0; 2:2 with P | 2:2 | 5 / 7 |

PGK bi-bi (BPG + ADP ⇌ 3PG + ATP), RE binding:

| Mechanism | Non-hyperbolic in | Identifiable / fitted |
|---|---|---|
| RE random bi-bi | none | 7 / 7 |
| Iso random bi-bi | none | 8 / 9 |
| Duplicated random bi-bi | none | 7 / 16 |
| Mnemonic random bi-bi | BPG, ADP | 12 / 13 |
| Flip-flop random bi-bi | all four | 15 / 16 |
| Flip-flop, ordered binding | all four | 11 / 12 |
| Flip-flop, binding constants shared across `E` and `Eprime` | all four | 10 / 10 |
| Mnemonic, binding constants shared | BPG, ADP | 10 / 10 |
| Flip-flop, bindings and chemistry both shared | none (Wegscheider forces Kiso = 1) | 7 / 8 |
| SS random bi-bi | all four | 15 / 15 |

The random-binding flip-flop matches SS random bi-bi in degree but not in economy. Savings
come only from sharing binding constants across conformations or from ordered binding. The
shared-binding flip-flop has an exact closed form: the RE random bi-bi law times
`(1 + a(AB + PQ/Keq)) / (1 + c·AB + d·PQ)`. At initial rates its activation stays below
2-fold, while stronger inhibition is allowed.

Size of the effect: sampled apparent Hill coefficients span about 0.05–1.9 across these
models, so each allows both positive and negative apparent cooperativity. Sampling gives lower
bounds on the extremes, not maxima.

## 3. Cleland iso mechanisms (goal 1)

Iso uni-uni adds one denominator term, S·P. Iso random bi-bi adds ADP·ATP·BPG,
ADP·ATP·3PG, ADP·BPG·3PG, ATP·BPG·3PG and ADP·ATP·BPG·3PG. Every iso variant stays degree 1
in each metabolite, so iso mechanisms change product-inhibition patterns and never produce
curvature on their own. Each carries one structurally unidentifiable fitted parameter (4 of 5
uni-uni, 8 of 9 random bi-bi).

## 4. PGK data and search results

**Data.** 386 points in 15 groups. Six groups run forward (BPG + ADP) and nine run reverse
(3PG + ATP); no group mixes directions. At most one product is present at a time: 3PG in
forward groups (Lavoinne 1983 Fig4/5) and ADP in reverse groups (Lavoinne 1979 Fig7/8, Lee
1975 Fig3A/B). The data mix rat-liver PGK (Lavoinne) and human-erythrocyte PGK (Lee, Ali).

**Where the curvature is.** Per-series Hill fits give apparent negative cooperativity in
Lavoinne 1979 Fig4 (ATP over 160-fold, h = 0.49–0.60), Fig7 (h = 0.58–0.73), and Ali 1976
Fig1A/3A (BPG and ADP, h = 0.43–0.69). No series is sigmoidal (h > 1.3). Letting every series
have its own V and K, Michaelis–Menten leaves loss 0.0055 and a Hill equation leaves 0.0009,
so within-series curvature accounts for at most about 0.0046 of loss. Ali 1976 conflicts with
Lavoinne 1983 over matched forward ranges (h ≈ 0.5 against h ≈ 1.0), so no single mechanism
fits both.

**Search baselines** (2026-09-17 runs, Keq = 760, best training loss by fitted-parameter
count):

| Params | No inhibitor modes | With inhibitor modes |
|---|---|---|
| 5 | 0.0615 | 0.0615 |
| 6 | 0.0407 | 0.0287 |
| 7 | 0.0303 | 0.0276 |
| 8 | 0.0261 | 0.0272 |
| 9 | 0.0243 | 0.0250 |
| 10 | 0.0232 | 0.0194 |
| 11 | 0.0228 | 0.0187 |
| 12 | 0.0215 | 0.0174 |
| 13 | 0.0196 | 0.0165 |

- RE random bi-bi scores 0.0615 (5 parameters). The best SS random bi-bi variants in the
  search score 0.0303 (7) and 0.0273 (9). The fully SS random bi-bi with independent
  constants (15 parameters) exceeds `max_param_count = 13` and was never fitted.
- Almost all of the 0.06 → 0.03 gap sits in the ADP-containing reverse groups (Lavoinne 1979
  Fig7/8, Lee 1975 Fig3B): the ADP product-inhibition pattern. An RE model with one extra ADP
  binding mode closes it (0.0287 at 6 parameters).
- Models that reproduce the ATP curvature are n = 1 `AllostericMechanism`s with SS binding
  (0.0194 at 10, 0.0165 at 13) and one 10-parameter single-conformation SS model with
  context-dependent binding constants (0.0198). In LOOCV (with inhibitor modes) the scores are
  0.0383 at 7 parameters (RE plus inhibitor modes), 0.0291 at 10 and 0.0271 at 12.
- The SS random bi-bi variants in the search share binding constants across contexts. That
  sharing makes them exactly hyperbolic in every substrate when products are absent; their
  fitted Hill coefficient on Lavoinne 1979 Fig4 is about 1.0 against the data's 0.5.

**Keq.** `identify_pgk.jl` uses Keq = 760 and `fit_pgk.jl` uses 1000. RE-model losses do not
depend on Keq, because no group contains reactants of both directions and the loss is centred
per group. SS and iso losses do: the 7-parameter SS model scores 0.0303 at 760, 0.0323 at
1000 and 0.0635 at 3000. Comparisons need one Keq.

**Literature.** PGK's hinge-bending domain closure is fast, and product release limits the
rate (Banks et al. 1979; Varga et al. 2006; Schmidt, Travers & Barman 1995). The source papers
attribute the curvature to anion or second-nucleotide sites (Scopes 1978; Lavoinne 1979;
Szilágyi & Vas 1998), and 20–40 mM sulfate linearizes the plots (Scopes 1978). No paper
proposes an iso, mnemonic or hysteretic mechanism for PGK. PGK is therefore a numerical
benchmark for goal 2, not mechanistic validation.

## 5. Package state

The derivation already handles a second conformation label. Thirteen hand-written variants
(iso, mnemonic, flip-flop, duplicated, linked, and PGK versions of each) matched a BigFloat
mass-action ground truth to relative error ≤ 3.4e-14, with `v` at equilibrium ≤ 6e-14.
`rate_equation` stays allocation-free at about 10–26 ns for the iso family, against 97 ns for
SS random bi-bi; V×τ is 4–6 against 336.

Gaps:

- **Enumeration builds one conformation.** `_catalytic_topologies` hardcodes `:E` in
  `_make_species` (`src/mechanism_enumeration.jl:13-21`), the start form (`:682`) and the
  cycle-closure test (`:353`). The rest of the pipeline keys species on conformation already.
- **Conformational change counts as chemistry.** `is_iso` means "no bound metabolite"
  (`src/types.jl:179`), so `_flux_carrying_groups` and `_expand_to_allosteric` treat a pure
  conformational step as chemistry. Only `_compute_numerator` separates them.
- **Numerator cut.** `_compute_numerator` skips pure conformational steps
  (`src/rate_eq_derivation.jl:538-542`). A mechanism whose only SS steps are conformational
  (slow induced fit plus slow relaxation) fails with a false "all-RE catalytic cycle" error.
- **kcat.** `_kcat_forward` (`src/rate_eq_derivation.jl:848-930`) overstates kcat by up to
  43× the true peak for mnemonic and flip-flop mechanisms. It affects post-fit rescaling of
  returned parameters only, not the loss.
- **Phantom parameter.** Every iso, mnemonic and flip-flop variant with independent bindings
  carries one structurally unidentifiable fitted parameter, which the beam would charge as a
  real one. The shared-binding variants carry none.
- **Dependent-parameter choice.** The Haldane-dependent parameter of an iso cycle (chemistry
  reverse rate or `k_Eprime_to_E`) and the canonical direction of `Eprime <--> E` fall to
  name ordering when both conformations bind substrate.
- **Allosteric plus iso.** Wrapping an iso scheme in an `AllostericMechanism` works under the
  implicit rule that only free `E` flips between A and I. Nobody chose that rule, and no test
  covers it.

## 6. Caveats

- **RE limit versus initial rate.** The RE (Cha) law evaluated as products → 0 and the
  mass-action limit at exactly zero product disagree for about 4% of the uni-uni pool. It
  happens when an RE product-binding step is the only fast link between two populated forms:
  at zero product the enzyme is trapped. Example:
  `E + S ⇌ E(S); E(S) <--> Eprime(P); Eprime(P) ⇌ Eprime + P; Eprime(P) ⇌ E(P); E(P) ⇌ E + P`
  predicts a finite initial rate, while mass action gives zero. An enumerator must not count an
  RE route through product complexes as a return path from `Eprime` to `E`.
- **Unproven converse.** "Dead blocks collapse" is proved. The converse, that every live
  multi-segment mechanism differs from all single-segment RE laws, is supported by sampling
  but unproven. The graph-based degree predictor held in every exhaustive scan but is also
  unproven.
- **Practical identifiability.** Structural rank overstates what data pin down. Sampled
  singular values of `d ln v / d ln θ` give 9–10 effective parameters for the random
  flip-flop (of 15) and 9–13 for SS random bi-bi (of 15).

## 7. Next step: PGK feasibility spike

Question: can a conformational kinetic-cooperativity mechanism beat the current
curvature-capable models on PGK with fewer parameters?

Probe: fit hand-written mechanisms to the PGK data with the HPC settings (Keq = 760,
`scale_k_to_kcat = 1.0`, CMA-ES). Baselines: RE random bi-bi, RE plus the ADP inhibitor mode,
and SS random bi-bi. Candidates: iso random bi-bi; mnemonic in both orientations; the random,
ordered and shared-binding flip-flops; and the most promising candidates again with the ADP
inhibitor mode. Success means beating about 0.0194 at 9 parameters or fewer with inhibitor
modes (or 0.026 at 7 or fewer without), while reproducing the ATP curvature of Lavoinne 1979
Fig4 (h ≈ 0.5–0.6).

## References

- Ali M, Brownstone YS (1976) Biochim Biophys Acta 445:89. doi:10.1016/0005-2744(76)90162-5
- Banks RD et al. (1979) Nature 279:773. doi:10.1038/279773a0
- Cleland WW (1963) Biochim Biophys Acta 67:104. doi:10.1016/0006-3002(63)91800-6
- Cornish-Bowden A, Cárdenas ML (1987) J Theor Biol 124:1
- Ferdinand W (1966) Biochem J 98:278. doi:10.1042/bj0980278
- Lavoinne A et al. (1979) Biochimie 61:1043. PMID 534662
- Lavoinne A et al. (1983) Biochimie 65:211. doi:10.1016/s0300-9084(83)80086-8
- Lee CS, O'Sullivan WJ (1975) J Biol Chem 250:1275. PMID 1112804
- Ricard J, Meunier JC, Buc J (1974) Eur J Biochem 49:195.
  doi:10.1111/j.1432-1033.1974.tb03825.x
- Schmidt PM, Travers F, Barman T (1995) Biochemistry 34:824
- Scopes RK (1978) Eur J Biochem 85:503. doi:10.1111/j.1432-1033.1978.tb12266.x
- Storer AC, Cornish-Bowden A (1977) Biochem J 165:61. doi:10.1042/bj1650061
- Szilágyi AN, Vas M (1998) Biochemistry 37:8551. doi:10.1021/bi973072k
- Varga A et al. (2006) FEBS Lett 580:2698. doi:10.1016/j.febslet.2006.04.024
