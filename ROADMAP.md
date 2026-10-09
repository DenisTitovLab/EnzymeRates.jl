<!-- ABOUTME: Known issues, planned work, ideas and declined proposals for the package. -->
<!-- ABOUTME: The developer backlog; docs/src/roadmap.md holds user-facing directions. -->

# Roadmap

Update an entry in the same commit as the change that fixes, adds or retires it.

## Known issues

- **The inactive state can give an `:EqualAI` constant a second, fitted name.**
  `_state_parts` (`src/rate_eq_derivation.jl`) names the inactive conformation's
  constants after its pruned graph (`_state_allo_mechanism`), so a shared constant takes
  a second name when pruning removes its group's naming representative (`_group_rep`)
  but not the whole group, or reverses a rapid-equilibrium step that metabolite
  progression leaves unoriented (`F ⇌ E` gives `K_F_to_E` in the active state and
  `K_E_to_F` in the inactive one). The second name enters `fitted_params` as a free
  constant, so the rate law describes a larger model than the one declared and the beam
  overcounts its parameters. Enumerated mechanisms meet the first case: when
  `E(A) + P ⇌ E(A, P)` and `E(Q) + P ⇌ E(P, Q)` share `K_EAP_to_EA_P` and the inactive
  state cannot reach `E(A)`, that state names the constant `K_EPQ_to_EQ_P`, and 766 of
  15,219 bi-bi and ping-pong allosteric mechanisms sampled within two expansion rounds
  gain one to three such parameters. The single-equation derivation under Planned
  removes the cause. A narrower fix, prototyped and measured, names and orients every
  inactive-state constant by the full mechanism's own step and removes all 766.
- **The two conformations fold Wegscheider ties separately.**
  `_dependent_param_exprs(am)` (`src/rate_eq_derivation.jl`) solves each conformation's
  single-symbol rename on its own and merges the two maps. When the conformations fold a
  tied pair in opposite directions, or one folds a constant the other still uses,
  `fitted_params` keeps a constant that thermodynamics ties to another: the model is
  larger than declared, though identifiable. Within two expansion rounds of the bi-bi and
  ping-pong reactions, 256 enumerated allosteric mechanisms fold a pair in opposite
  directions and 512 more keep a folded constant in the other conformation. The
  single-equation derivation under Planned removes the cause.
- **A kinetic group cannot tie one metabolite exchange at two conformations.**
  A group holding `E(I) + J ⇌ E(J) + I` and the same exchange at conformation `F` raises
  "change different pairs of conformations", although neither step changes conformation;
  the same exchange at two forms of one conformation passes. `_orient_tied_steps`
  (`src/types.jl`) compares the conformation pair of every step that metabolite
  progression leaves tied with the lead step's pair, and a step that keeps its
  conformation has the pair `(E, E)` or `(F, F)`. The enumerator never builds such a
  group, so only hand-written mechanisms meet it. A fix would compare pairs only for
  tied steps that change conformation and orient the others by their consumed and
  released metabolites.

## Planned

- **Derive both conformations from one catalytic rate equation.**
  Each conformation is derived on its own graph: `_state_parts` rebuilds the inactive one
  as a pruned mechanism (`_state_allo_mechanism`), and each solves its own Wegscheider
  rename, which causes the first two known issues. Derive the catalytic mechanism once
  instead, with one symbol per constant and the constraints solved on that one system,
  and obtain each conformation from it: relabel a `:NonequalAI` constant for the inactive
  conformation, and remove an `:OnlyA` or `:OnlyI` step by a limit (`K → ∞` for a
  rapid-equilibrium binding, both `k → 0` for a steady-state step). Haldane and
  Wegscheider relations make such a limit come with partner limits, so both conformations
  stay thermodynamically consistent by construction. Take each limit as the leading-order
  term in a scaling parameter, not by substitution: zeroing both rate constants of a step
  that holds the graph together zeroes every spanning tree and leaves 0/0. The design
  should also settle which free form carries `L` when a hand-written mechanism has two
  residual-free free forms in different steady-state segments; the first in species
  order does today. Branch `mwc-solve-then-limit` holds a July design of this approach
  (spec b067aa0 and 651a821, plan 9c95f12, prototype 6f0011e), whose prototype matched
  the mass-action ground truth on all six uni-uni tag combinations and would delete the
  `d_free` cross-weighting; `test/allosteric_ground_truth.jl` is the oracle.
  Priority: high.
- **Stop emitting K-type variants whose L is a phantom at one catalytic subunit.**
  At multiplicity 1, `_expand_to_allosteric` (`src/mechanism_enumeration.jl`) emits
  K-type variants whose L the data cannot determine, so the fitter optimizes a flat
  direction and the beam counts the model one parameter too high. This hits every
  reaction that allows one catalytic subunit, which is the default. No exact criterion
  exists yet: derive a structural test on the dead-inactive law `v = N_A/(Q_A + L·Q_I)`,
  or apply the exact rank (see Ideas) to this move's children;
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
- **Anchor a measured kcat at its assay concentrations.**
  `rescale_parameter_values` (`src/rate_eq_derivation.jl`) divides by `_kcat_forward`, an
  analytic saturating limit behind a run of NaN, "no kcat components" and mis-scaling
  bugs. Evaluating `rate_equation` at the assay concentrations is unambiguous and mirrors
  the experiment: add a `kcat_concs` keyword and thread it like `scale_k_to_kcat` through
  `FittingProblem` and `identify_rate_equation`. A wider version anchors the relative
  default at a reference point and deletes `_kcat_forward` (~150 lines), but changes
  every reported rescaled constant. Priority: medium-low.
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
- **Build the Wegscheider rename into the constraint columns.**
  Pass 2 of `_build_wegscheider_rename_map` folds one binding constant into another, but
  the absorbed symbol stays a constraint column, so both `_dependent_param_exprs` methods
  (`src/thermodynamic_constr_for_rate_eq_derivation.jl`, `src/rate_eq_derivation.jl`)
  strip the resulting dummy parameter afterwards; the allosteric one re-derives both
  conformations' polynomials to keep any symbol they still reference. Dropping renamed
  symbols in `_assemble_constraints` would delete both filters; the hard case is an
  `:EqualAI` symbol folded in one conformation and used in the other. Pass 2 itself is
  load-bearing and stays, and the derivation output must stay byte-identical. The
  single-equation derivation above may absorb this. Priority: low.
- **Enumerate Cleland iso mechanisms.**
  The derivation handles iso mechanisms, but no enumeration move generates them, and each
  carries a dwell-time-gauge phantom unless the exact rank (see Ideas) or a
  "switch-product = 1" constraint row comes first. `Species` carries a conformation
  field, and `_requires_hyperbolic_catalysis` (`src/mechanism_enumeration.jl`) is the
  extension hook. Two-conformation kinetic cooperativity stays parked: on PGK data no
  conformational model beat MWC with steady-state steps at equal parameter count, so
  reopen it only for a dataset whose curvature MWC cannot fit. Priority: low.

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
- **Test the compat bounds in a scheduled CI job.**
  `.github/workflows/CI.yml` tests only the latest Julia release with the newest
  dependencies, while `Project.toml` declares `julia = "1.10"` and a lower bound for every
  dependency, so a wrong Julia or dependency floor surfaces only in a user's environment,
  and a dependency release that breaks the package only at the next CI run. A scheduled
  workflow that runs the suite at the lower and upper bounds of Julia and the dependencies
  would catch both without adding its run time to every push or pull request.

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
