<!-- ABOUTME: Known issues, planned work, ideas and declined proposals for the package. -->
<!-- ABOUTME: The developer backlog; docs/src/roadmap.md holds user-facing directions. -->

# Roadmap

Update an entry in the same commit as the change that fixes, adds or retires it.

## Known issues

- **The stored facts make the ter-ter `init_mechanisms` slower and larger.**
  Each `Mechanism` stores its modifiers, groups, step kinds, directions toward products and
  step sides. Across the 251,107 mechanisms the ter-ter `init_mechanisms` testset retains,
  the live heap grew 48% (1,272 to 1,883 MiB) and the peak resident set from 2.3 to 3.1 GB.
  The extra garbage collection makes the call about 20% slower (94 s to 117 s on Julia 1.13);
  its time limit of 400 s excludes garbage collection, about 23% of the call on CI, since a
  slow macOS runner reached 430 s with it and a fast one took 77 s. The step sides cost most:
  two strings per step, 668 of the 5,408 retained bytes per mechanism. Denis raised the
  limit instead of taking one of three measured options: (A) keep already-sorted input
  vectors instead of copying them, which lets a mechanism share the caller's vectors (72 s,
  24% under the time before the facts were stored); (B) compute `step_sides` on demand
  (102 s, 13% under the time with the facts, 8% over the time before them); (D) reuse the
  variants the flip search in `_seed_variants` already builds (with A and B, 51 s, 45%
  under the time before the facts). No benchmark or pin yet covers the live bytes of a held
  result set.
  - **On Julia 1.11 the full suite can peak above 4 GB.**
    With the full-suite command in `.claude/CLAUDE.md`, Julia 1.11.9 peaked at 4.12-4.15 GB (3.84-3.86
    GiB) in two runs; Julia 1.12 and 1.13 peaked at 4.0 and 3.8 GB in one run each, after the stored facts
    grew the retained mechanisms (2.6-3.0 and 2.6-2.7 GB before). The ter-ter
    `init_mechanisms` testset (`test/Enumeration/test_mechanism_enumeration.jl`) leaves about
    1 GB of garbage. No test value retains the garbage: a full collection after the testset
    frees it (live bytes fall from 0.94 to 0.055 GiB). CI runs Julia 1.11 only in the
    lower-bound job, on a runner with more memory, so nothing checks the compat floor against
    the 4 GB that `.claude/CLAUDE.md` sets aside for the suite.
  - **macOS CI takes 23-39 minutes for the tests, up from about 8.**
    Ubuntu and Windows stay at 15-17 minutes. The macOS runner has 7 GiB against their 16:
    in one instrumented run its heap reached 5.5 GiB with 3.0 GiB resident, and the
    enumeration test file took 1,165 s, 600 s of it in garbage collection, against 273 s on
    Ubuntu. A `--heap-size-hint=2500M` on macOS gave 25.6 and 32.4 minutes in two runs,
    inside the range without it, and slows Ubuntu and Windows to about 20 minutes. macOS
    runner speed varies about 30% between runs, so one run proves nothing. The next step is
    to measure which stored fields hold the memory; a change to the `Mechanism`
    representation needs the maintainer's approval.
- **A steady-state step that only reaches dead-end forms fits both rate constants.**
  When a steady-state step's only role is to reach forms that lead nowhere else, as
  `E(I) + J <--> E(J) + I` does when `E(J)` binds nothing more, no net flux crosses it, so
  only the ratio of its two rate constants enters the rate law.
  `parameters(m)` lists both, and the law carries their common factor
  uncancelled in numerator and denominator. One constant, the step's equilibrium
  constant, describes it. Only hand-written mechanisms meet this: enumerated dead-end
  bindings stay at rapid equilibrium.
- **A fit is not thread-safe.**
  `loss!` (`src/Fitting/fit_rate_equation.jl`) writes its problem's shared
  `log_ratios_buffer`, and `_fit_loss` copies each point into the one `x_buf` of its
  `FitObjective`, so an optimizer that scores points on several threads at once, such as
  CMA-ES with `multi_threading` passed through `solver_kwargs`, would race on both buffers.
  Fits run in parallel only across `pmap` workers, each of which builds its own problem.
- **A rapid-equilibrium isomerization of a covalent intermediate can be unidentifiable.**
  In 5 of 751 enumerated plain mechanisms in the derivation probes, an isomerization such
  as `E(B; residual = A - P) ⇌ E(Q)` joins two forms whose rapid-equilibrium weights
  carry the same concentrations, so the data determine only a combination of its
  constant and its neighbour's, and `parameters(m)` lists one constant
  more than the identifiable rank. The exact rank under Ideas would count it.
- **A law the guard admits can be 0/0 where one substrate and every product are zero.**
  The guard checks the faces of the rewrite spec's section 8 (every product zero, every
  substrate zero, each single reactant zero, each also with every modifier zero), never a
  corner that zeroes one substrate with every product. Measured with the guard before
  0.11.0, 357 of 65,423 census mechanisms pass it yet have a law that is 0/0 at such a
  corner: 293 plain ter-bi, 38 plain ping-pong and 26 allosteric ter-bi mechanisms; there
  264 of these laws have no limit and 93 a removable 0/0, census mechanism
  `d7a09031d5dd2fc6` at B = P = Q = 0 among them. Another 12 allosteric ping-pongs that pass
  the guard have a non-unique steady state at such a corner. `rate_equation` returns NaN
  there, where a well-posed law returns 0. A fit scores the NaN as it scores the 0 a
  well-defined law predicts there: `FittingProblem` refuses a zero rate, and the loss scores
  a NaN prediction as a sign mismatch, the penalty a 0 prediction gets
  (`test/Derivation/test_admission.jl`). Adding these corners to the faces would drop at
  least these 369 mechanisms, 93 of whose laws have a limit there.
- **The guard rejects some allosteric mechanisms whose law is the limit.**
  A gate run on the guard as built (see the next entry) found three kinds:
  - **Too large to search.** When a conformation has more than 20,000 ways to switch off its
    excluded steps, or the two conformations' consistent ways form more than 5,000 pairs,
    the guard does not search (`_MAX_CHOICES_PER_CONFORMATION`, `_MAX_READING_PAIRS`,
    `src/Derivation/guard.jl`) and tries only its three built readings. 441 gate mechanisms
    are too large to search. The 33 among them with no consistent built reading are rejected
    with `_TOO_MANY_READINGS`, whatever their law: 27 of the 2,000 ter-ter taggings, 2 of
    3,000 ter-bi, 1 of 3,000 bi-ter and 3 of the admission tests' mechanisms. Without a
    regulator, 20 of the 4,885 census ter-bi starting mechanisms keep no allosteric child,
    since their only K-type is too large to search. Raising the caps costs time: the slowest
    guard call, on a ter-bi mechanism of 25 steps with 2,157 pairs, takes 0.3 s.
  - **Refused as a bridge.** The face condition refuses vanishing transitions that open a
    faster route between two of the law's forms when it cannot prove the route harmless.
    Gate entry 53246, a ter-ter tagging that holds an excluded rapid-equilibrium step that
    binds and reacts at once, is refused so, yet it has one consistent reading, under which
    the ε-path oracle converges on every face along standard directions. A wider set-aside
    rule would admit it: one that takes a cut-off form the enzyme does not reach on the face
    along with a group set aside when any of its alive transitions leads into it, where the
    rule as built wants all of them to (`_without_set_aside_groups`,
    `src/Derivation/admission.jl`). It would also set aside unfed forms that empty into a
    set-aside group only through one another, which the rule as built keeps. One argument
    makes it sound: the face has one closed class, which free E reaches, so a form free E
    does not reach holds no enzyme at any ε and any speed. The wider rule was not built; it
    changes verdicts, so it needs the gate. Entry 53145, of the same kind and refused as a
    bridge too, is not shown valid: under the one reading tested, the ε-path oracle finds a
    limit that depends on the path at zero products and at zero P.
    The suite pins two more, `bridge_deadend` and R-B3 (`test/Derivation/test_admission.jl`,
    "conservative rejections the oracle admits stay rejected under their reading"). A rule
    that also set aside groups touching forms in equilibrium for another reason, such as a
    dead-end branch, admitted nothing in any population.
  - **Refused by a face check though the rate equals the law's.** Of the mechanisms of
    8 steps or fewer that a face check refuses, sampled 40 per population, these have a
    limit along corner directions under some consistent reading: 9 of 40 `:OnlyA` taggings
    of the uni-uni, bi-bi and ping-pong starting mechanisms (an estimated 127 of 566), 16 of
    40 of the admission tests' mechanisms (47 of 118), 7 of 40 ter-bi (16 of 89), 4 of the 8
    ping-pong taggings with both halves' chemistry `:OnlyA` (below), 1 of 40 bi-ter (2 of
    86), and none of the `:OnlyI`, ter-ter, bi-bi split and inhibitor ones. The guard
    reports each with its first consistent reading's `_PUMPED` or `_FILLED`. Nearly all are
    ping-pongs. Under the reading the oracles accept, the guard refuses 8 of the 9 `:OnlyA`
    taggings by the filling check and 1 as a bridge, 8 times at B's binding to the covalent
    form and once at P's. On such a face the cut-off forms hold most of the enzyme, yet the
    rate equals the law's because the inactive conformation's kept chemistry has the active
    one's constants. No simple rule found separates these from the wrong ones.
- **The guard's face checks under a searched reading rest on numerical evidence.**
  The filling check and the set-aside rule (Developer page, "Well-posedness guard") have
  arguments, not proofs, that the face condition with them suffices under every reading. A
  wrong law the guard admitted would not show in a fit: it is a smooth function that can fit
  well, with a wrong mechanism label and unreliable extrapolation toward the face where it
  fails. The protection is a gate run once on the guard as built. It ran the guard on 78,554
  mechanisms in ten populations: `:OnlyA` taggings of the uni-uni, bi-bi and ping-pong
  starting mechanisms (9,616), taggings with an `:OnlyI` group (43,311), ter-ter taggings
  (2,000), K-type taggings of the bi-bi split children (9,878), ping-pong taggings (50),
  bi-bi with a competitive inhibitor, ter-bi and bi-ter taggings (3,000 each), one expansion
  round's allosteric children (1,990) and the admission tests' mechanisms (2,709). None of
  the 29,517 the built readings admit is rejected. It ran the ε-path and RE-speed oracles on
  the 1,575 that only the search or the set-aside rule admits, all of them of 8 steps or
  fewer and samples of the rest: all 591 small ones agree along corner directions, 598
  larger ones agree along standard directions, and 100 ter-ter ones agree along corner
  directions; at two subunits the steady state equals the law for 40 of 40 sampled. With the
  filling check off, the guard admits 96 small `:OnlyA` taggings more (138 in all), and 18
  of them fail; standard directions miss these, so corner directions are the test.
  `test/Derivation/test_admission.jl` keeps a sampled population test with corner
  directions. No hand-written test shows a one-way fast vanishing transition that passes
  the filling check because its target drains fast, along rapid-equilibrium steps that do
  not vanish, into the law's forms (the check's third clause); only the sampled test covers
  it. The gate's scripts, logs and report (`report.md`) are on the maintainer's
  machine, outside the repository, in `claude-scratchpad/lean-mechanism-spike/gate/`.
- **An allosteric ping-pong whose covalent forms empty only by binding B cannot stop both
  half-reactions.**
  PDH and OGDH are ping-pong and allosteric, so allosteric ping-pong must stay in the
  search. The guard admits a ping-pong whose inactive conformation runs neither
  half-reaction when the covalent forms can lie higher in energy there: each half's step
  into them loses its rate and its step out keeps it, so they empty whatever the ratio of
  the two vanishing rates. Of the 25 bi-bi ping-pong `init_mechanisms` with a covalent form,
  each with both halves' chemistry `:OnlyA` and with every step `:OnlyA`, it admits 40 of
  the 50. The 3 rejected with only their chemistry `:OnlyA` (`_FILLED`) have no limit under
  any reading: their covalent forms empty only by binding B. Stopping both halves there
  needs the conformations to flip at the covalent form as well, with an equilibrium constant
  of its own and the multi-subunit rate formula rederived for two flip forms; that belongs
  with iso mechanisms. The 7 rejected with every step `:OnlyA` (`_PUMPED`) hold the 4 valid
  ping-pongs counted under "The guard rejects some allosteric mechanisms whose law is the
  limit". What the search does today:
  - **The tagging move stops both half-reactions where the guard admits it.**
    Its first candidate base tags every chemistry group `:OnlyA` (`_tagging_bases`,
    `src/Enumeration/moves.jl`). On the six-step and three-step bi-bi ping-pongs, and on the
    six-step one whose A and P groups each hold an abortive binding, the guard admits it, so
    it is the one base. On the six-step ping-pong every allosteric child kept at the exits,
    and every allosteric child of those on the same graph, turns over at zero products: 44
    with a regulator declared, at one subunit and again at two, and 18 without one. A
    one-half base arises only where the guard rejects the every-chemistry set, as on the
    four-step bi-bi ping-pong whose A binding is fused at steady state and whose B binding
    is fused at rapid equilibrium; several bases occur in 48 of the 4,885 census ter-bi
    parents.
  - **Change-tag cannot give the six-step bi-bi ping-pong a live inactive conformation.**
    On the V-type that stops both halves, change-tag's relaxation of both half-reactions to
    `:NonequalAI` excludes nothing, but class 3 rejects it: at zero substrates, and at zero
    B, the covalent forms of both conformations lie on one rapid-equilibrium segment, joined
    only through free E, so the rate's limit depends on how fast each conformation's steps
    run (`_SPLIT`, `_LIMIT_SPEED`). `expand_mechanisms` drops that child, and none of the 44
    children above has an inactive conformation that turns the reaction over. On the
    three-step ping-pong, whose steps are all steady-state, the same relaxation passes.
  - **A bare tagging base is fitted only beside a regulator.**
    The tagging move emits a base only beside a declared regulator, and no move reaches the
    bare base: change-tag relaxes tags only to `:NonequalAI`, and no move removes a site.
    Yet a bare base can be well-posed with its `L` identifiable. The three-step ping-pong's
    every-chemistry base, its only base, fits 6 parameters at rank 6 over the parent's 5 at
    two subunits (at rank 5 at one, where `L` is a phantom). The four-step ping-pong's
    one-half base, which stops the first half-reaction, fits 6 at rank 6 at one and at two
    subunits. Both turn over at zero products. Whether the enumerator should emit a bare
    base, at least above one subunit, is the maintainer's decision.
- **A batch's progress-log header can print twice.**
  When an L phantom's source's fit throws in the phantom's own batch, `run_batch!`'s log
  closure (`src/identify_rate_equation.jl`) prints the batch header again above the demotion
  lines, and the pre-fit summary has already counted those phantoms as inherited. Rows and
  CSVs are right. A phantom's source is usually fitted a batch earlier, so no identify run
  in the suite reaches this path; the fix (print the header before `_process_batch` and let
  `log` only indent) needs an identify fixture that does.
- **A hand-built step that does not balance fails only inside the joint solve.**
  A step whose sides take up and give off different net metabolites, such as `E(P) → E`
  taking up `P`, passes the `Mechanism` constructor and the guard. It fails later, in the
  joint solve, with "Cycle 1 produces metabolite change not proportional to net
  reaction". The constructor or the guard should reject it with a message that names the
  step. Only hand-built mechanisms meet this: the macros and the enumerator write balanced
  steps.
- **A regulator measured only at saturating levels gets an arbitrary K.**
  With an activator at 0 and 100× K, restarts that fit perfectly reported 0.015× to 2× the
  true K. The data do not constrain K there, so the fit returns any value in the box and
  lists it in `at_bound` only when it lands at an edge.
- **`loss!` runs 2-5× slower in about one process in nine.**
  `loss!` stores each row's concentrations and the parameters in one block of its stack
  frame and passes them by pointer to `rate_equation`, which is not inlined (for every gated
  law but Uni-Uni). When a 4 KB page boundary of the stack falls inside that block, one of
  `rate_equation`'s 16- or 32-byte loads straddles the page on every row and costs 20-30 ns
  more per row, measured in a Lima VM on Apple Silicon with 4 KB pages. The stack's offset
  within its page is random per process and fixed for its life, so the slow mode looks
  random: 29 of 256 offsets slow at least one of the 42 gated laws. With address-space
  randomization off, one 2 KB environment variable reproduces it on unchanged code, and the
  same laws fail at the same offset at an older commit (47e581ce). The `loss!` gate in
  `test/Derivation/test_rate_eq_derivation.jl` therefore fails falsely in about one full
  suite in nine, on several laws at once, each at 26-38 ns per row: its 20 batches use fresh
  problems, which vary the heap, but share one stack depth. Rerun such a suite. Two fixes
  are open, for the maintainer to choose: run each batch a few frames deeper (test only: 0
  false failures in 10,752 (law, offset) pairs, against 147 now); or inline `rate_equation`
  at `loss!`'s call site, which makes `loss!` 15-35% faster per row on small and medium laws
  and cuts the failing pairs from 147 to 10, at a compile cost not measured (see the Idea
  "Inline `rate_equation` in `loss!` and emit `muladd` in `_horner`"). Not measured: 16 KB
  pages (macOS on Apple Silicon), x86, and the cost to a whole fit; that a fit in such a
  process evaluates `loss!` slower for its whole life is inferred. The earlier observations
  may be this effect: three large laws 15-25% over the gate in one local process, and one
  `FittingProblem` up to about 40% slower than a fresh one. Windows CI running up to 5% over
  on six laws is more likely coverage counters, unproven. The investigation, with its
  scripts and scans, is `claude-scratchpad/loss-timing-investigation/report.md` on the
  maintainer's machine.

## Planned

- **Shrink the module interfaces.**
  The package aims at deep modules: few public names, simple arguments, the complexity
  inside (`.claude/CLAUDE.md`). Four places fall short. Priority: high.
  - **`Mechanisms` has 25 public names.**
    Some serve one caller and expose how it works: `segment_walk` and `present_species`
    are used only by `Derivation`, `check_conformation_label` and `check_decodable_names`
    only by `DSL`, and `has_unbalanced_cycle` takes six vectors. Each either moves to its
    one caller or folds into a function that hides it.
  - **`identify_rate_equation` has at least 18 keyword arguments,** such as
    `eq_complexity_filter = 337`. Each becomes an internal constant unless users have a
    real reason to change it, as they do for optimizer settings.
  - **`DSL.parse_mechanism` is public and no other module in `src/` uses it.** Make it
    private if nothing needs it.
  - **The guard's internals carry names and arguments nothing needs.**
    `_ARC_CYCLE` and `_ISLAND` are never reported, so `_interior_verdict` could return one
    "inconsistent" outcome from a Boolean acyclicity test, and `_Fault`'s docstring could
    become a list. `exclusion_components(am, c::Symbol)` takes a Symbol that its three
    callers make from an Int, and lacks the leading underscore. `_admitting_reading` and
    `_admission` are called only by tests. Of `_eps_network`'s callers, only
    `_every_step_cures` passes a `re` other than `_re_flags(m)`, and it could build its
    variant with `_with_steps` instead. Three words, component, group and part, name the
    forms a conformation's kept steps join.
- **Decide what the guard keeps that no measured mechanism needs.**
  Whether each part below goes is the maintainer's decision: he chose the three built
  readings, and guard code does not go on a search alone. Priority: medium.
  - **The second and third built readings.** Once the guard searches the readings,
    `_ALL_STOPPED` changes no verdict in any of the gate's ten populations, and
    `_OUTWARD_LOSE_AFFINITY` was not measured. Measure both against the mechanisms the
    guard before 0.11.0 admits, those too large to search in particular, and remove what
    changes nothing: about 45 lines of `src/Derivation/conformations.jl` and their mirror
    in `test/mechanism_definitions_for_test_enzyme_derivation.jl`.
  - **Joint R1.** The branch of `_face_condition` (`src/Derivation/admission.jl`) that
    refuses a bridge the two conformations' vanishing transitions open only together is
    reached through `_guard` by no mechanism searched: 20,096 variants of its fixture and
    the 54,927 taggings of the gate's `:OnlyA`, `:OnlyI` and ter-ter populations.
- **Log the time of each search phase.**
  `progress.log` carries no timestamps, so a cluster run cannot say where its wall time goes.
  In the 0.10.0 runs the average worker was idle for 30% of PFKP's time and 75% of PGK's
  (Developer page, What a fit costs). Workers connect within 15 s and the master's CPU covers
  under half of that idle time, so workers waiting for the slowest fit of a batch are the
  leading suspect. A timestamp on each line, and each batch's slowest fit, would rank all the
  search-speed ideas below. Priority: high.
- **Profile the search on the HPC runs.**
  In the 0.10.0 runs, 64% of PFKP's and 75% of PGK's fits reach a loss within 10⁻⁸ of another
  fit at the same parameter count. These may be laws the key misses, equal only after a
  reparametrization (`_canonical`, `src/Derivation/key.jl`), or laws that differ only outside
  the data. Separating the two would size a saving of up to 3–4×. Still unmeasured on these
  runs: fits whose parent fails the final parsimony cutoff, children and survivors per move,
  zero-gain `change_allo_state` children, and how many PGK rows identify flags `L_phantom`.
  Priority: medium.
- **Check the energy-level rule on ter-ter mechanisms beyond the gate's sample.**
  The gate (see the Known issue "The guard's face checks under a searched reading rest on
  numerical evidence") drew 2,000 of the ter-ter `:OnlyA` taggings of `init_mechanisms`. The
  guard admits 253 of them only through the search or the set-aside rule: the 19 of 8 steps
  or fewer agree with both oracles along corner directions, 150 of the 234 larger ones along
  standard directions, and 100 of the 253, evenly spaced, along corner directions (a median
  6.9 s each). 223 of the 2,000 are too large to search. The other ter-ter taggings, and
  ter-ter mechanisms after an expansion round, are unchecked. Priority: medium.
- **Make per-regulator multiplicities such as `A(1, 2)` work.**
  `@enzyme_reaction` rejects a regulator's multiplicity list such as `A(1, 2)`; the
  syntax returns with this feature. `_expand_to_allosteric` and
  `_expand_add_allosteric_regulator` give every regulatory site the catalytic
  multiplicity. The derivation already powers each site by its own multiplicity, so only
  the DSL (`_parse_regulator_entries`, `src/DSL/read.jl`; `Regulator`,
  `src/Mechanisms/types.jl`) and enumeration change: those two moves, and
  `_expand_merge_regulatory_sites`, which must merge only sites of equal multiplicity.
  Keep a bare `A` at the catalytic multiplicity, widen only on an explicit list, and
  reject `I(m)` for a competitive inhibitor. Priority: medium-low.
- **Drop `change_allo_state` relaxations that free no parameter, if measurement finds any.**
  `_expand_change_allo_state` (`src/Enumeration/moves.jl`) has no gain test, unlike
  the split and RE→SS moves, and in a past LDH run it produced 20 of 33 zero-gain edges.
  A child whose law equals one already fitted, up to renaming its constants, shares that
  fit through the key (`_canonical`, `src/Derivation/key.jl`) and costs no fit; a zero-gain
  child with a law of its own still costs one. First count, over the bi-bi allosteric
  depth-2 population, the children whose `fitted_parameter_count` does not exceed their
  parent's and whose key is new, and confirm them with `_testhelper_identifiable_rank`. A
  naive filter can drop a different model with the same count, so check reachability
  before filtering. Priority: medium-low.
- **Decide whether to keep the hyperbolic drops above one subunit.**
  Above one subunit the tagging and change-tag moves drop a child whose inactive law is not
  hyperbolic (`_hyperbolic_where_required`, `src/Enumeration/moves.jl`). Neither drop
  removes a candidate the guard admits in any population searched: 81,687 tagging and
  378,543 change-tag candidates at two subunits on bi-bi and census ping-pong parents, and
  1,046,178 tag vectors on hyperbolic parents of at most 7 groups. On the six-step bi-bi
  ping-pong whose A and P groups each hold an abortive binding, the tagging move's one base
  stops both half-reactions; at two subunits the drop removes none of its 15 K-types, and
  the exits keep 5. Why the drops find nothing, as an argument rather than a proof: a
  tagging or change-tag child has no `:OnlyI` catalytic group, so its active law is the
  parent's, which the move requires to be hyperbolic. Its inactive law keeps only forms of
  free E's group. An excluded rapid-equilibrium step between two of them would be pinned,
  which the guard rejects, so each rapid-equilibrium segment of the inactive law is a
  segment of the parent cut down to the kept forms, and the law's degree in a metabolite
  cannot exceed the parent's. The exception is a segment the cut splits in two: the cut-off
  forms between the parts then drain into both, which the face condition refuses as a
  bridge, unless the set-aside rule sets them aside; that case the argument does not cover.
  Tests pin both drops on candidates the guard rejects anyway. Priority: low.
- **Enumerate Cleland iso mechanisms.**
  The derivation handles iso mechanisms, but no enumeration move generates them, and each
  carries a dwell-time-gauge phantom unless the exact rank (see Ideas) or a
  "switch-product = 1" constraint row comes first. An `EnzymeForm` carries its conformation
  label, and `_requires_hyperbolic_catalysis` (`src/Enumeration/filters.jl`) is the
  extension hook. Two-conformation kinetic cooperativity stays parked: on PGK data no
  conformational model beat MWC with steady-state steps at equal parameter count, so
  reopen it only for a dataset whose curvature MWC cannot fit. Revisit allosteric
  ping-pong with this work: a flip at a second free form is the same question (see the
  Known issue "An allosteric ping-pong whose covalent forms empty only by binding B cannot
  stop both half-reactions"). Priority: low.

## Ideas

- **Extend `plot_rate_equation_on_data`.**
  Considered and left out: overlaying the fits of two mechanisms, a fit summary (measured
  against predicted, error per group), a tick rule for log axes, and a keyword to override
  the inferred x metabolite.
- **Spend less time compiling each law.**
  Compiling a fresh law takes 0.10-0.19 s, 9-12% of the PFKP and PGK searches' fit time and
  30-50% on PGK's laws of 6-8 parameters. `loss!`'s second pass does not depend on the mechanism,
  so moving it behind a function that compiles once may save part of this; the saving is not
  measured.
- **Reduce what CMA-ES itself costs per evaluation.**
  At 0.10.0, `loss!` no longer takes nearly all of a restart: an evaluation took 22 µs against
  10 µs in `loss!` on a 20-parameter PFKP law, and 8.1 µs against 6.0 µs at the median of 60
  PGK laws, so any faster `rate_equation` gains only on the `loss!` part. The solver's
  eigendecomposition takes 5-14% of a warm fit, and the solver with the rest of the fit up
  to 20% of PFKP's. It belongs to CMAEvolutionStrategy, so this means a solver option or a
  change upstream; neither is explored.
- **Count identifiable parameters with an exact modular rank.**
  Several phantom classes make the beam count a model above its true dimension: chain
  flanks kept after a later split, copy-copy twin pairs, `:NonequalAI` copies and the
  iso dwell-time gauge. The one-subunit L has an exact test of its own: a rescaling
  certificate (`_phantom_weights`, `src/Derivation/key.jl`) finds an L phantom whose
  inactive conformation turns nothing over, and identify lets it read the fit of its law at
  L = 0 when that law holds its other fitted constants and a mechanism with it was fitted in
  an earlier batch or is fitted in the phantom's own.
  The rank of ∂v/∂log θ modulo a prime at random points of the derived law (`_law`,
  `src/Derivation/engine.jl`) is exact up to a negligible failure probability, took 3.5 ms
  per mechanism in a prototype, and has a finite-difference oracle in
  `_testhelper_identifiable_rank`
  (`test/mechanism_definitions_for_test_enzyme_derivation.jl`). Run it over several
  saved runs to find the other phantom classes; performance may keep it out of `src/`,
  where it would otherwise replace the beam's parameter count. It must never filter
  children, since a phantom child can carry families its parent lacks.
- **Apply the chain-flank exclusion at one catalytic subunit.**
  The flip move keeps a chain's flanks at rapid equilibrium only in plain mechanisms
  (`_chain_flank_groups`, `src/Enumeration/moves.jl`): an allosteric mechanism has
  none, since above one subunit a flank's flip can show through the two
  conformations. At one subunit an L phantom's flank flip stays invisible, so the move
  emits a phantom whose law at L = 0 no plain mechanism of the search holds, and
  identify fits it as an ordinary candidate (census #3268; the phase-5 research found
  62 such bi-bi variants on 16 keys in a sampled round 2). Applying the exclusion to
  one-subunit variants would drop these fits; first check that no other law becomes
  unreachable.
- **Choose the number of restarts per fit.**
  Over 120 sampled fits per enzyme, the best of 10 restarts took half the time of 20 and came
  within 1% of the best of 20 in 98% of PFKP's fits (worst 10.7% above it) and in all of PGK's.
  Stopping once the best loss recurs 4 times gives 1.7× and 2.4×, and 2% of PFKP's fits end
  more than 5% high. A loss above the beam's cutoff drops a mechanism for good, so test a rule
  by comparing whole search results, not single fits.
  More restarts have a case too: on five PFKP laws of 16 or more parameters whose fits had
  missed, 20, 40 and 80 restarts gave 13, 3 and 0 misses in 60 fits, and one such miss
  dropped the line to PFKP's 18-parameter law from the beam. No whole search has been run
  at 10 or at 40 restarts, and which way to move the default is undecided.
- **Start a child's restarts from its parent's fit.**
  A child usually holds its parent as a limit, so one restart could begin at the parent's
  parameters. Not measured.
- **Pin each law's operation counts.**
  The `rate_equation` time bound tolerates the spread between machines, so undoing a 3× gain
  on a small law passes it: the MWC tetramer's former 20.3 ns sits under its 34 ns bound. A
  count of multiplications, additions and power calls per law does not depend on the machine
  and would catch this, at the cost of re-pinning when a law changes on purpose.
- **Inline `rate_equation` in `loss!` and emit `muladd` in `_horner`.**
  A call-site `@inline` lets LLVM hoist the arithmetic on parameters alone out of the row loop:
  `loss!` runs 9–35% faster over seven sampled laws, with up to 75% more first-call compile
  time. `muladd` gives 2–18%, and the generated code then breaks the binary `+` and `*` shape
  that the tests assert. Both leave the loss unchanged to 12 digits. They were measured
  separately on aarch64, and neither has been measured on the cluster's x86 cores.
- **Seed each identify fit from its mechanism's content.**
  `fit_rate_equation` (`src/Fitting/fit_rate_equation.jl`) draws each restart's starting
  point with `randn` from the running task's random stream, and `_process_batch`
  (`src/identify_rate_equation.jl`) fits a batch's representatives in the iteration order of
  a `Dict` keyed by the rate-law key (`_canonical`), filled in the order the moves emit the
  mechanisms. So a change of emission order or of a law's key moves the starting points of
  every later fit, and a fit-dependent test can flip with no change to any model. A seed per
  fit from a process-stable content key, mixed with one master seed drawn per
  `identify_rate_equation` call, would make each fit independent of the processing order
  while `Random.seed!` still reproduces a run; every saved fit would move once.
- **Give the thermodynamics one row reduction.**
  `src/Derivation/thermo.jl` holds two exact row reductions, `_rref!` (the joint solve, in
  `Rational{Int}` retried in `Rational{BigInt}`) and `_rref_partition` (in
  `Rational{BigInt}`). `_rref_partition` serves only class 3's model point (`_model_point`,
  `src/Derivation/limit_test.jl`), whose integer null space (`_integer_nullspace`) goes
  through `_rational_nullspace`. Computing that null space with `_rref!` would leave one.
- **Decide admission's face condition exactly.**
  Admission's face condition (`_face_condition`, `src/Derivation/admission.jl`), with the
  filling check and the set-aside rule, is sound in every test but not exact: it rejects
  mechanisms whose limit exists (Known issues, "The guard rejects some allosteric
  mechanisms whose law is the limit"). In an equilibrium bridge (R-B1) a tagged binding
  joins forms its conformation's other steps keep apart, yet those forms already sit at
  equilibrium, so the limit is deletion's (`bridge_deadend` in
  `test/Derivation/test_admission.jl`); the set-aside rule admits such a bridge only when
  the law's own rapid-equilibrium steps hold the forms it touches in equilibrium. No exact
  narrower structural rule is known: the one candidate tried, T0++, admits faces the
  ε-path oracle rejects. The exact test compares Freidlin–Wentzell exponents over the
  admissible directions: the cut-off forms vanish when, for every rise of the tagged
  constants, every per-member rate and every order of the two limits, their cheapest
  in-tree costs more than the law network's; each comparison is one linear-programming
  feasibility test, not yet prototyped.
## Decided against

- **Scheduled CI and Downgrade runs.**
  A weekly or monthly run would catch a dependency release that breaks the package or stops
  the lower bounds resolving, with no push here. CI runs on pull requests and pushes to `main`,
  and `.github/workflows/Downgrade.yml` on pull requests that change `Project.toml`; the next
  such run surfaces either break.
- **A global ±ln(1e9) box, physical ranges or a plausibility warning.** Either needs a
  declared unit. `lb` and `ub` in `fit_rate_equation` stay, and the box follows the data.
- **`lb` and `ub` in `identify_rate_equation`.** The search sets each fit's box from its
  data.
- **Move dead-end decorations out of `init_mechanisms` into an expansion move.** It is not
  clear this saves time, since the move would have to run again on every child.
- **Separate the metabolite names inside enzyme-form names (`Ac + CoA ⇌ AcCoA`).**
  `@enzyme_reaction` rejects a reaction whose names spell one name two ways and names both
  spellings, and renaming one metabolite removes the collision, so no separator is needed.
- **Map BlackBoxOptim's convergence stop to `Success` upstream in OptimizationBase.** Only
  BlackBoxOptim users meet the warning, `test/Fitting/test_fitting.jl` captures it, and the
  docs and HPC scripts use CMA-ES.
