<!-- ABOUTME: Known issues, planned work, ideas and declined proposals for the package. -->
<!-- ABOUTME: The developer backlog; docs/src/roadmap.md holds user-facing directions. -->

# Roadmap

Update an entry in the same commit as the change that fixes, adds or retires it.

## Known issues

- **A steady-state step that only reaches dead-end forms fits both rate constants.**
  When a steady-state step's only role is to reach forms that lead nowhere else, as
  `E(I) + J <--> E(J) + I` does when `E(J)` binds nothing more, no net flux crosses it, so
  only the ratio of its two rate constants enters the rate law. `fitted_params` lists
  both, and the law carries their common factor uncancelled in numerator and
  denominator. One constant, the step's equilibrium constant, describes it. Only
  hand-written mechanisms meet this: enumerated dead-end bindings stay at rapid
  equilibrium.
- **A fit is not thread-safe.**
  `loss!` (`src/fitting.jl`) writes its problem's shared `log_ratios_buffer`, and
  `_fit_loss` copies each point into the one `x_buf` of its `FitObjective`, so an
  optimizer that scores points on several threads at once, such as CMA-ES with
  `multi_threading` passed through `solver_kwargs`, would race on both buffers. Fits run
  in parallel only across `pmap` workers, each of which builds its own problem.
- **`eq_complexity_filter` counts only the active conformation of an allosteric
  mechanism.** `_eq_complexity` (`src/rate_eq_derivation.jl`) measures an allosteric
  mechanism's V×τ on its catalytic graph, the active conformation, alone, while the
  generated MWC rate equation evaluates the polynomials of both conformations, so the filter
  undercounts what an allosteric equation costs and can admit one that evaluates more
  terms than its limit.
- **`BigInt` arithmetic inflates the GC's live-byte count.**
  Three functions that compute in `BigInt` raise `Base.gc_live_bytes` by bytes that no
  live object holds: `_onlya_offenders` (`src/thermo.jl`; a `Rational{BigInt}` null
  space), the tagging and change-tag moves' `:OnlyA` test, by about 15 KiB per call; the
  counter that `_split_gain_test` uses for a plain mechanism
  (`_partition_independent_count`; a `Rational{BigInt}` row reduction) by 18 KiB; and
  `_eq_complexity` (a `BigInt` determinant, skipped for a single rapid-equilibrium
  segment) by 15 KiB. A 6×6 determinant outside the package does the same in `BigInt`
  and `Rational{BigInt}`, not in `Rational{Int}`. The moves run the `:OnlyA` test on
  every tag set they try: 100,000 calls add 1.46 GiB and 12 full sweeps. The joint solve
  adds nothing, on a plain mechanism or on a concrete allosteric one; it runs in
  `Rational{Int}` on 10,226 enumerated bi-bi and ter-bi mechanisms without overflow. The
  count stays up after a full collection: 100,000 calls that start 1.47 GiB high make 52
  full sweeps in 10.0 s, against 12 in 6.5 s from a clean start. The count holds no
  memory and raises no peak: a process whose count reaches 2.8 GiB peaks at 1.2 GB RSS,
  and the enumeration tests' count after a full collection stays near 1.8 GiB whether
  their peak RSS is 2.2 or 3.5 GB. The higher peak came from mechanisms that kept vectors
  at the capacity their builds grew to; `_canonicalize` stores each at its length, and
  `test/test_reaction_graph.jl` pins it. Filtering the moves through admission replaces
  the `:OnlyA` test, not the other two. A loop that calls one of these and prints
  `Base.gc_live_bytes()` and `Base.gc_num().full_sweep` shows it. The figures come from a
  random-order bi-bi of nine steady-state steps, plain and as a two-subunit allosteric
  mechanism with `:OnlyA` on both bindings of A and on the chemical step and `:EqualAI`
  elsewhere (Julia 1.12.5, aarch64, `--heap-size-hint=2500M`).

- **The tagging and change-tag moves misjudge `:OnlyA` groups that mix a binding and a
  chemistry step.** Their filter, `_onlya_offenders` (`src/thermo.jl`), skips every
  `:OnlyA` group that holds a chemistry step, including one that also holds a plain
  binding of the same metabolite (enumeration builds these by merging and tags every
  allosteric one `:OnlyA`), and gives the bindings of one group a single rate of
  divergence. The moves therefore emit tag sets that admission rejects, such as a mixed
  `:OnlyA` group whose plain member sits on a square of `:EqualAI` bindings (the exit
  drops them), and miss tag sets that admission admits, such as a split `:OnlyA` binding
  group whose limit exists. Filtering the moves through admission replaces the test.
- **Splitting a tied kinetic group changes `eq_hash` without changing the model.**
  When a square of rapid-equilibrium bindings already ties two binding constants,
  splitting their kinetic group leaves the rate function and the fitted names unchanged,
  but the split part's constant becomes a dependent with its own line under
  `# Wegscheider constraints:`, and the law names it, so the printed equation and the
  dedup key (`_rate_eq_dedup_key`, `src/identify_rate_equation.jl`) differ, and identify
  fits both. The rename-proof key under Planned removes this.
- **A rapid-equilibrium isomerization of a covalent intermediate can be unidentifiable.**
  In 5 of 751 enumerated plain mechanisms in the derivation probes, an isomerization such
  as `E(B; residual = A - P) ⇌ E(Q)` joins two forms whose rapid-equilibrium weights
  carry the same concentrations, so the data determine only a combination of its
  constant and its neighbour's, and `fitted_params` lists one constant more than the
  identifiable rank. The exact rank under Ideas would count it.
- **A law the guard admits can be 0/0 where one substrate and every product are zero.**
  The guard checks the faces of the rewrite spec's section 8 (every product zero, every
  substrate zero, each single reactant zero, each also with every modifier zero), never a
  corner that zeroes one substrate with every product. Of 65,423 census mechanisms, 357
  pass the guard yet have a law that is 0/0 at such a corner: 293 plain ter-bi, 38 plain
  ping-pong and 26 allosteric ter-bi mechanisms; there 264 of these laws have no limit and
  93 a removable 0/0, census mechanism `d7a09031d5dd2fc6` at B = P = Q = 0 among them.
  Another 12 allosteric ping-pongs that pass the guard have a non-unique steady state at
  such a corner. `rate_equation` returns NaN there, where a well-posed law returns 0. Fits
  never meet it: `FittingProblem` refuses a zero rate, and the loss scores a NaN prediction
  as a sign mismatch, the penalty a 0 prediction gets (`test/test_admission.jl`). Adding
  these corners to the faces would drop at least these 369 mechanisms, 93 of whose laws
  have a limit there.
- **On Julia 1.10 the full suite can peak above 4 GB.**
  With the full-suite command in `.claude/CLAUDE.md`, Julia 1.10.12 peaked at 3.30-4.36 GB
  in seven runs, and an eighth was killed at 4.49 GB, most likely by the OOM killer of a
  7.7 GB machine; Julia 1.12 peaks at 2.6-3.0 GB and Julia 1.13 at 2.6-2.7 GB. The ter-ter
  `init_mechanisms` testset (`test/test_mechanism_enumeration.jl`) leaves about 1 GB of
  garbage. On 1.10 the full suite still counts it after the enumeration file, 1.83-1.84 GiB
  live there against 0.55-0.65 GiB with one `GC.gc()` after the testset; with the file run
  alone, the count drops at the ter-ter seed testset. No test value retains the garbage: a
  full collection after the testset frees it (live bytes fall from 0.94 to 0.055 GiB), and
  that `GC.gc()` brings the suite to 2.93-3.05 GB. CI runs Julia 1 alone, so nothing checks
  the compat floor against the 4 GB that `.claude/CLAUDE.md` sets aside for the suite.

## Planned

- **Rewrite the derivation around the reaction graph.**
  `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md` records the design,
  decided between 2026-10-09 and 2026-10-12. Two PRs remain: admission, the
  well-posedness guard and the flip at free E; then the engine, printing, a rename-proof
  dedup key, code-generation locals and the CSV columns that go with them (`key`,
  `fit_source_id`, `L_phantom`). The table in its section 16 lists the entries that each
  PR retires. Priority: high.
- **Take `:OnlyA` steps to their limit in the inactive conformation.**
  Both conformations come from one joint solve on one graph, but the inactive law still
  deletes every `:OnlyA` group and prunes the forms that deletion strands. An `:OnlyA`
  binding should lose its affinity for the inactive conformation instead (`K → ∞`, or
  `kon → 0` with `koff` fixed for a steady-state binding), and an `:OnlyA` chemistry step
  be blocked (both rate constants → 0, their ratio set by thermodynamics); each step of an
  `:OnlyA` group already keeps its own class. A tag set is admitted when an exact
  feasibility test passes and every form the free enzyme cannot reach in the inactive
  conformation lies behind a binding, on its bound side; the limit then equals deletion
  and pruning. A ping-pong with both half-reactions `:OnlyA` fails this, so the tagging
  move will stop one half per V-type child instead (spec section 7, and the entry
  below). The conformations will flip at free E, which every mechanism holds; until
  then the boundary check admits a second free form of an allosteric mechanism only in
  rapid equilibrium with the first (`_assert_one_flip_segment`, `src/admission.jl`),
  where the choice only rescales `L`. The rewrite's admission PR brings the limits,
  admission and the flip at free E.
  Branch `mwc-solve-then-limit` holds a July design of this approach (spec b067aa0 and
  651a821, plan 9c95f12, prototype 6f0011e), whose prototype matched the mass-action
  ground truth on all six uni-uni tag combinations and would delete the `d_free`
  cross-weighting; `test/allosteric_ground_truth.jl` is the oracle. Priority: high.
- **Let an allosteric ping-pong stop both half-reactions in the inactive conformation.**
  PDH and OGDH are ping-pong and allosteric, so allosteric ping-pong must stay in the
  search. With the conformations flipping only at free E, a ping-pong whose inactive
  conformation runs neither half-reaction strands the covalent forms between the halves:
  their weight, and with it the rate, depends on how fast each half stops (by up to 100%
  at one catalytic subunit), so admission rejects it. The admission PR keeps allosteric
  ping-pong reachable by stopping one half per V-type child, which covers PDH, whose
  phosphorylation stops E1 and leaves E2 and E3 running. Stopping both halves needs the
  conformations to flip at the covalent form as well, with an equilibrium constant of
  its own, and the multi-subunit rate formula rederived for two flip forms. Take it up
  with iso mechanisms, which add free forms of their own. Priority: medium.
- **Flag one-subunit K-type variants whose L is a phantom instead of fitting them.**
  At multiplicity 1, `_expand_to_allosteric` (`src/mechanism_enumeration.jl`) emits
  K-type variants whose L the data cannot determine, so the fitter optimizes a flat
  direction and the beam counts the model one parameter too high. This hits every
  reaction that allows one catalytic subunit, which is the default. The rewrite keeps the
  variants in the frontier, since their children can carry laws that the parent's
  children lack, but stops fitting them: a variant whose inactive conformation binds
  nothing and does not catalyze, and whose L an exact rank test shows adds no
  identifiable rank, takes its parent's fit with L = 0 and is flagged.
  `docs/src/identify/enumeration_engine.md` documents the phantom. Priority: medium.
- **Profile the first v0.9.0 HPC runs.**
  Every search-cost number on hand comes from runs made before 0.8.0: 66% of equations
  duplicated under renaming, `to_allosteric` producing 47% of LDH children, 94% of fits
  non-parsimonious. On the next LDH, PGK and PFKP runs, measure fits whose parent fails
  the final parsimony cutoff, children and survivors per move, same-function equations
  that `eq_hash` misses (`_rate_eq_dedup_key`, `src/identify_rate_equation.jl`),
  zero-gain `change_allo_state` children, and how often PGK meets the one-subunit L
  phantom. These numbers rank the open search-speed and deduplication ideas.
  Priority: medium.
- **Make per-regulator multiplicities such as `A(1, 2)` work.**
  The DSL parses and stores a regulator's multiplicity list (`RegulatorMults`,
  `src/types.jl`), but no enumeration move reads it: `_expand_to_allosteric` and
  `_expand_add_allosteric_regulator` give every regulatory site the catalytic
  multiplicity, so `A(1, 2)` does nothing. The derivation already powers each site by its
  own multiplicity, so only enumeration changes: those two moves, and
  `_expand_merge_regulatory_sites`, which must merge only sites of equal multiplicity.
  Keep a bare `A` at the catalytic multiplicity, widen only on an explicit list, and
  reject `I(m)` for a competitive inhibitor. Priority: medium-low.
- **Drop `change_allo_state` relaxations that free no parameter, if measurement finds any.**
  `_expand_change_allo_state` (`src/mechanism_enumeration.jl`) has no gain test, unlike
  the split and RE→SS moves, and in a past LDH run it produced 20 of 33 zero-gain edges,
  each costing a fit. First count, over the bi-bi allosteric depth-2 population, the
  children whose `_independent_param_count` does not exceed their parent's, and confirm
  them with `_testhelper_identifiable_rank`. A naive filter can drop a different model
  with the same count, so check reachability before filtering. Priority: medium-low.
- **Enumerate Cleland iso mechanisms.**
  The derivation handles iso mechanisms, but no enumeration move generates them, and each
  carries a dwell-time-gauge phantom unless the exact rank (see Ideas) or a
  "switch-product = 1" constraint row comes first. A `Form` carries its conformation
  label, and `_requires_hyperbolic_catalysis` (`src/mechanism_enumeration.jl`) is the
  extension hook. Two-conformation kinetic cooperativity stays parked: on PGK data no
  conformational model beat MWC with steady-state steps at equal parameter count, so
  reopen it only for a dataset whose curvature MWC cannot fit. Revisit allosteric
  ping-pong with this work: a flip at a second free form is the same question (see
  "Let an allosteric ping-pong stop both half-reactions"). Priority: low.

## Ideas

- **Count identifiable parameters with an exact modular rank.**
  Several phantom classes make the beam count a model above its true dimension: the
  one-subunit L, chain flanks kept after a later split, copy-copy twin pairs,
  `:NonequalAI` copies and the iso dwell-time gauge. The rank of ∂v/∂log θ modulo a prime
  at random points of the rate polynomials (`src/sym_poly_for_rate_eq_derivation.jl`) is
  exact up to a negligible failure probability, took 3.5 ms per mechanism in a prototype,
  and has a finite-difference oracle in `_testhelper_identifiable_rank`
  (`test/test_mechanism_enumeration.jl`). Run it over several saved runs to find the major
  phantom classes; performance may keep it out of `src/`, where it would otherwise replace
  the beam's parameter count. It must never filter children, since a phantom child can
  carry families its parent lacks.
- **Seed each identify fit from its mechanism's content.**
  `fit_rate_equation` (`src/fitting.jl`) draws each restart's starting point with `randn`
  from the running task's random stream, and `_process_batch`
  (`src/identify_rate_equation.jl`) fits a batch's representatives in the iteration order
  of a `Dict` keyed by `eq_hash`, filled in the order the moves emit the mechanisms. So a
  change of emission order or of a law's printed text moves the starting points of every
  later fit, and a fit-dependent test can flip with no change to any model. A seed per fit
  from a process-stable content key, mixed with one master seed drawn per
  `identify_rate_equation` call, would make each fit independent of the processing order
  while `Random.seed!` still reproduces a run; every saved fit would move once.
- **Test the compat bounds in a scheduled CI job.**
  `.github/workflows/CI.yml` tests only the latest Julia release with the newest
  dependencies, while `Project.toml` declares `julia = "1.10"` and a lower bound for every
  dependency, so a wrong Julia or dependency floor surfaces only in a user's environment,
  and a dependency release that breaks the package only at the next CI run. A scheduled
  workflow that runs the suite at the lower and upper bounds of Julia and the dependencies
  would catch both without adding its run time to every push or pull request.
- **Give the thermodynamics one row reduction and one cycle basis.**
  `src/thermo.jl` holds two exact row reductions, `_rref!` (the joint solve, in
  `Rational{Int}` retried in `Rational{BigInt}`) and `_rref_partition` (in
  `Rational{BigInt}`), and two cycle-basis builders, `_fundamental_cycles` and
  `_thermodynamic_constraints`, which raise the same error. The `:OnlyA` guard's use of
  `_rref_partition` and `_thermodynamic_constraints` goes when admission replaces the
  guard, but the split move's counter (`_partition_independent_count`) keeps them;
  computing its rank with `_rref!` over the joint solve's cycles would leave one of each.
- **Admit the equilibrium bridges and the vanishing-fed basins that admission rejects.**
  Admission's face condition (`_face_condition`, `src/admission.jl`) is sound but not
  exact: it rejects two rare classes whose limit exists. In an equilibrium bridge (R-B1) a
  tagged binding joins forms its conformation's other steps keep apart, yet those forms
  already sit at equilibrium, so the limit is deletion's (`bridge_deadend` in
  `test/test_admission.jl`). In a stranded basin fed only through vanishing sources (R-B2)
  the basin's weight vanishes along every path, but no structural test proves it (`BS1`).
  The phase-4 research found R-B1 on 332 faces and R-B2 on 21, all in fuzzed or hand-built
  mechanisms, none among the census mechanisms or their `:OnlyI` variants. No exact
  narrower structural rule is known: the one candidate tried, T0++, admits faces the
  ε-path oracle rejects. The exact test compares Freidlin–Wentzell exponents over the
  admissible directions: the basin vanishes when, for every rise of the tagged constants,
  every per-member rate and every order of the two limits, its cheapest in-tree costs more
  than the law network's; each comparison is one linear-programming feasibility test, not
  yet prototyped.

## Decided against

- **Move dead-end decorations out of `init_mechanisms` into an expansion move.** It is not
  clear this saves time, since the move would have to run again on every child.
- **Separate the metabolite names inside enzyme-form names (`Ac + CoA ⇌ AcCoA`).**
  `@enzyme_reaction` rejects a reaction whose names spell one name two ways and names both
  spellings, and renaming one metabolite removes the collision, so no separator is needed.
- **Map BlackBoxOptim's convergence stop to `Success` upstream in OptimizationBase.** Only
  BlackBoxOptim users meet the warning, `test/test_fitting.jl` captures it, and the docs
  and HPC scripts use CMA-ES.
- **Load identify CSVs written by EnzymeRates 0.7 or earlier, or that name a mechanism by
  its type string (`mechanism_type`) instead of its text and id.** The loader would be a
  backward-compatibility shim; rerunning is cleaner, and the CSVs read as plain tables.
- **Replay the 0.8 LOOCV random stream so seeded runs reproduce its fold scores.** Nothing
  needs this backward-compatibility code, and unseeded workers never made HPC fold scores
  bit-reproducible.
- **Re-rescale the rate constants of saved allosteric fits.** Only CSVs written before
  0.9.0 carry the mis-scaled constants, current output is correct, and losses and model
  selection never depended on the scale.
