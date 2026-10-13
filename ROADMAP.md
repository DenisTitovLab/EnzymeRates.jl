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
  mechanism's V×τ on its active-state graph (`steps(am)`) alone, while the generated MWC
  rate equation evaluates the polynomials of both conformations, so the filter
  undercounts what an allosteric equation costs and can admit one that evaluates more
  terms than its limit.
- **`BigInt` arithmetic inflates the GC's live-byte count.**
  Three functions that compute in `BigInt` raise `Base.gc_live_bytes` by bytes that no
  live object holds: `_onlya_haldane_violation` (`src/thermo.jl`; a `Rational{BigInt}`
  null space) by about 16 KiB per call, the counter that `_split_gain_test` uses for a
  plain mechanism (`_partition_independent_count`; a `Rational{BigInt}` row reduction) by
  14 KiB and `_eq_complexity` (a `BigInt` determinant, skipped for a single
  rapid-equilibrium segment) by 8 KiB. A 6×6 determinant outside the package does the same
  in `BigInt` and `Rational{BigInt}`, not in `Rational{Int}`. An `AllostericMechanism`
  with an `:OnlyA` group of plain bindings runs the guard when built, so every lift of
  such a compiled mechanism pays it, as does each such candidate of the allosteric
  `_split_gain_test`. 100,000 such lifts add 1.5 GiB and 16 full sweeps. The joint solve
  adds nothing, on a plain mechanism or on a concrete allosteric one; it runs in
  `Rational{Int}` on 10,226 enumerated bi-bi and ter-bi mechanisms without overflow. The
  count stays up after a full collection: a lift loop that starts 1.4 GiB high makes 115
  full sweeps in 29 s, against 16 in 20 s from a clean start. Admission under Planned
  replaces the guard, not the other two. A loop that calls one of these and prints
  `Base.gc_live_bytes()` and `Base.gc_num().full_sweep` shows it. The figures come from a
  random-order bi-bi of nine steady-state steps, plain and as a two-subunit allosteric
  mechanism with `:OnlyA` on both bindings of A and on the chemical step and `:EqualAI`
  elsewhere (Julia 1.12.5, aarch64, `--heap-size-hint=2500M`).

- **The `:OnlyA` guard misjudges groups that mix a binding and a chemistry step.**
  `_onlya_haldane_violation` (`src/thermo.jl`)
  skips every `:OnlyA` group that holds a chemistry step, including one that also holds a
  plain binding of the same metabolite (enumeration builds these by merging and tags
  every allosteric one `:OnlyA`), and gives the bindings of one group a single rate of
  divergence. Hand-built mechanisms show both errors: it admits a mixed `:OnlyA` group
  whose plain member sits on a square of `:EqualAI` bindings, whose inactive conformation
  is not a limit of the model, and rejects a split `:OnlyA` binding group whose limit
  exists. The `:OnlyA` semantics under Planned replace the guard.
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

## Planned

- **Rewrite the derivation around the reaction graph.**
  `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md` records the design,
  decided between 2026-10-09 and 2026-10-12. Three PRs remain: the graph mechanism type
  with two macros, `show` and a new CSV format; admission, the well-posedness guard and
  the flip at free E; then the engine, printing, a rename-proof dedup key and
  code-generation locals. The table in its section 16 lists the entries that each PR
  retires. Priority: high.
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
  below). The conformations will flip at free E, which every mechanism must then
  contain; until then the constructor admits a second free form only in rapid
  equilibrium with the first (`_assert_one_flip_segment`), where the choice only rescales
  `L`. The rewrite's admission PR brings the limits, admission and the flip at free E.
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
- **Give ping-pong residual forms names that are Julia identifiers.**
  `_species_name` (`src/types.jl`) writes residuals as `_res_+A_-P`, so a name such as
  `k_EA_to_EP_res_+A_-P` needs `var"…"` to pass as a parameter, and the destructuring
  line of `rate_equation_string` is not valid Julia. An encoding such as `_res_pA_mP`
  fixes this but renames every ping-pong parameter and moves its `eq_hash`, and it stays
  unambiguous only if `EnzymeReaction` also rejects `_` in metabolite names.
  Priority: medium-low.
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
  "switch-product = 1" constraint row comes first. `Species` carries a conformation
  field, and `_requires_hyperbolic_catalysis` (`src/mechanism_enumeration.jl`) is the
  extension hook. Two-conformation kinetic cooperativity stays parked: on PGK data no
  conformational model beat MWC with steady-state steps at equal parameter count, so
  reopen it only for a dataset whose curvature MWC cannot fit. Revisit allosteric
  ping-pong with this work: a flip at a second free form is the same question (see
  "Let an allosteric ping-pong stop both half-reactions"). Priority: low.

## Ideas

- **Encode a mechanism's Sig without compiling per mix of step shapes.**
  `_sig_of` (`src/types.jl`) builds the Sig tuples from generators, so Base's collect
  widening compiles once per new mix of step shapes. Over a 6,475-mechanism bi-bi run in one
  process (the seeds, then two expansion rounds; 5,045 allosteric) that totals about 1,500
  method instances, 0.5 MB of native code and 1 s, 85% of it within the first 250
  mechanisms, so a search pays it early and once per worker. Building each tuple from a
  `Vector{Any}` should remove it. `test/test_compile_budget.jl` gates the total for the
  239 bi-bi seeds.

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

## Decided against

- **Move dead-end decorations out of `init_mechanisms` into an expansion move.** It is not
  clear this saves time, since the move would have to run again on every child.
- **Separate the metabolite names inside enzyme-form names (`Ac + CoA ⇌ AcCoA`).** Renaming
  one metabolite removes any collision, and the error names both forms; the same answer
  covers a separator only for colliding reactions and a collision check at reaction build.
- **Map BlackBoxOptim's convergence stop to `Success` upstream in OptimizationBase.** Only
  BlackBoxOptim users meet the warning, `test/test_fitting.jl` captures it, and the docs
  and HPC scripts use CMA-ES.
- **Load identify CSVs written by EnzymeRates 0.7 or earlier.** The loader would be a
  backward-compatibility shim; rerunning is cleaner, and the CSVs read as plain tables.
- **Replay the 0.8 LOOCV random stream so seeded runs reproduce its fold scores.** Nothing
  needs this backward-compatibility code, and unseeded workers never made HPC fold scores
  bit-reproducible.
- **Re-rescale the rate constants of saved allosteric fits.** Only CSVs written before
  0.9.0 carry the mis-scaled constants, current output is correct, and losses and model
  selection never depended on the scale.
- **Have the `Mechanism` constructors check step metabolites against the reaction.** The
  DSL rejects undeclared names and the enumerator builds only declared ones, so the check
  would catch nothing and cost every enumerated mechanism.
