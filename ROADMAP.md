<!-- ABOUTME: Known issues, planned work, ideas and declined proposals for the package. -->
<!-- ABOUTME: The developer backlog; docs/src/roadmap.md holds user-facing directions. -->

# Roadmap

Update an entry in the same commit as the change that fixes, adds or retires it.

## Known issues

- **The ter-ter starting set is large to hold.**
  The ter-ter `init_mechanisms` testset (`test/Enumeration/test_starting_mechanisms.jl`)
  keeps 251,107 mechanisms, which hold 723 MiB of live heap (3,019 bytes each). Its call
  peaks at 1.83 GiB resident (1.82-1.83) and takes 61.9 s (60.4-62.0), 27.6 s of it in
  garbage collection, against a time limit of 400 s that excludes garbage collection.
  `candidate_mechanisms` on the ter-ter reaction returns 271,621 candidates in 104.0 s
  (103.2-104.5), 56.9 s of it in garbage collection, peaks at 2.09 GiB resident
  (2.09-2.10), and its result holds 930 MiB. These are medians of three runs, ranges in
  brackets, by `claude-scratchpad/lean-mechanism-spike/gate/measure_terter.jl` on the
  maintainer's machine (a fresh process per run under `--startup-file=no` and
  `--heap-size-hint=2500M`, Julia 1.13.4, a warm-up on ter-bi). The testset "the
  enumeration holds and allocates a bounded number of bytes" in the same file bounds the
  bytes each candidate of the bi-bi starting set holds (`Base.summarysize`), and the bytes
  that set and one expansion allocate; the ter-ter testset bounds the bytes its call
  allocates. Two options may cut the time, by amounts not
  measured on the mechanisms as stored: keep already-sorted input vectors instead of
  copying them, so that a mechanism shares the caller's vectors; and reuse the variants the
  flip search in `_seed_variants` already builds.
  - **On Julia 1.11 the full suite can peak above 4 GB.**
    With the full-suite command in `.claude/CLAUDE.md`, Julia 1.11.9 peaked at 4.12-4.15 GB
    (3.84-3.86 GiB) in two runs; Julia 1.12 and 1.13 peaked at 4.0 and 3.8 GB in one run
    each, when the stored facts grew the retained mechanisms (2.6-3.0 and 2.6-2.7 GB before
    them). The ter-ter `init_mechanisms` testset
    (`test/Enumeration/test_starting_mechanisms.jl`) leaves about 1 GB of garbage. No test
    value retains the garbage: a full collection after the testset frees it (live bytes fall
    from 0.94 to 0.055 GiB). CI runs Julia 1.11 only in the lower-bound job, on a runner
    with more memory, so nothing checks the compat floor against the 4 GB that
    `.claude/CLAUDE.md` sets aside for the suite.
  - **macOS CI takes 23-39 minutes for the tests, up from about 8.**
    Ubuntu and Windows stay at 15-17 minutes. The macOS runner has 7 GiB against their 16.
    Before pull request #109, in one instrumented run its heap reached 5.5 GiB with 3.0 GiB
    resident, and the enumeration test file took 1,165 s, 600 s of it in garbage
    collection, against 273 s on Ubuntu; a `--heap-size-hint=2500M` on macOS gave 25.6 and
    32.4 minutes in two runs, inside the range without it, and slowed Ubuntu and Windows to
    about 20 minutes. Four runs on the macOS runner timed the ter-ter `init_mechanisms`
    call and the suite:

    | Run | Commit | Ter-ter call | Of it in garbage collection | Suite |
    |:--|:--|--:|--:|--:|
    | `main` | 261c7eb4 | 318.9 s | 171.0 s | 37 min |
    | 44986745521 | 61f3bf5b | 187.1 s | 76.6 s | about 22 min |
    | 44990337442 | 6b37a084 | 78.2 s | 13.9 s | 10.5 min |
    | 45001777789 | dc936ab2 | 286.0 s | 116.1 s | 35 min |

    The last three ran pull request #109 on nearly the same code and took 10.5 to 35
    minutes: no gain on macOS is shown. `test/Enumeration/test_starting_mechanisms.jl`
    bounds the bytes each candidate holds and the bytes the enumeration allocates. In the
    last run each candidate held 1,594 bytes on Ubuntu, Windows and macOS, and the bytes
    allocated agreed to within 2.5%, so these bytes do not explain a slow macOS run. In that
    run the ter-ter call spent 170 s outside garbage collection and 116 s in it on macOS,
    against 83 s and 21 s on Ubuntu, and macOS peaked at 2.8 GiB resident against
    3.0-3.6 GiB elsewhere. The cause of the slow macOS runs is not known.
- **A steady-state step that only reaches dead-end forms fits both rate constants.**
  When a steady-state step's only role is to reach forms that lead nowhere else, as
  `EI + J <--> EJ + I` does when `EJ` binds nothing more, no net flux crosses it, so
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
  as `EB_res_plusA_minusP ⇌ EQ` joins two forms whose rapid-equilibrium weights
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
  the guard have a non-unique steady state at such a corner. The compiled rate law returns
  NaN there, where a well-posed law returns 0. A fit scores the NaN as it scores the 0 a
  well-defined law predicts there: `FittingProblem` refuses a zero rate, and the loss scores
  a NaN prediction as a sign mismatch, the penalty a 0 prediction gets
  (`test/Fitting/test_fit_rate_equation.jl`). Adding these corners to the faces would drop
  at least these 369 mechanisms, 93 of whose laws have a limit there.
- **The guard rejects some allosteric mechanisms whose law is the limit.**
  A gate run on the guard as built (see the next entry) found three kinds:
  - **Too large to search.** When a conformation has more than 20,000 ways to switch off its
    excluded steps, or the two conformations' consistent ways form more than 5,000 pairs,
    the guard does not search (`_MAX_CHOICES_PER_CONFORMATION`, `_MAX_READING_PAIRS`,
    `src/Derivation/rate_equation_correctness_guard.jl`) and tries only its three built
    readings. 441 gate mechanisms are too large to search. The 33 among them with no
    consistent built reading are rejected with `_TOO_MANY_READINGS`, whatever their law: 27
    of the 2,000 ter-ter taggings, 2 of 3,000 ter-bi, 1 of 3,000 bi-ter and 3 of the
    admission tests' mechanisms. No starting mechanism is left without an allosteric child:
    with no regulator declared and one subunit, each of the 5,158 starting mechanisms of the
    census ter-bi reaction `A[C] + B[N] + D[X] → P[CN] + Q[X]` (`candidate_mechanisms(rxn)`,
    which replaces each degenerate one of its 4,885 `init_mechanisms` by its flip children)
    gets at least one from `candidate_mechanisms(rxn, parent)`, a median of 26 (measured
    2026-11-25). Raising the caps costs time: the slowest guard call, on a ter-bi mechanism
    of 25 steps with 2,157 pairs, takes 0.3 s.
  - **Refused as a bridge.** The face condition refuses vanishing transitions that open a
    faster route between two of the law's forms when it cannot prove the route harmless.
    Gate entry 53246, a ter-ter tagging that holds an excluded rapid-equilibrium step that
    binds and reacts at once, is refused so, yet it has one consistent reading, under which
    the ε-path oracle converges on every face along standard directions. A wider set-aside
    rule would admit it: one that takes a cut-off form the enzyme does not reach on the face
    along with a group set aside when any of its alive transitions leads into it, where the
    rule as built wants all of them to (`_without_set_aside_groups`,
    `src/Derivation/allosteric_rate_equation_correctness_guard.jl`). It would also set aside
    unfed forms that empty into a set-aside group only through one another, which the rule
    as built keeps. One argument makes it sound: the face has one closed class, which free E
    reaches, so a form free E does not reach holds no enzyme at any ε and any speed. The
    wider rule was not built; it changes verdicts, so it needs the gate. Entry 53145, of the
    same kind and refused as a bridge too, is not shown valid: under the one reading tested,
    the ε-path oracle finds a limit that depends on the path at zero products and at zero P.
    The suite pins two more, `bridge_deadend` and R-B3
    (`test/Derivation/test_allosteric_rate_equation_correctness_guard.jl`, "conservative
    rejections the oracle admits stay rejected under their reading"). A rule that also set
    aside groups touching forms in equilibrium for another reason, such as a dead-end
    branch, admitted nothing in any population. A wider rule could also give the test for
    hyperbolic laws something to drop. Above one subunit the enumerator returns no mechanism
    whose conformations' laws are not hyperbolic (`_hyperbolic_where_required`,
    `src/Enumeration/candidate_mechanism_type.jl`), and today that drops nothing that
    `_enumerator_keeps` keeps: every one of the 13,678 state assignments it keeps on the
    bi-bi starting mechanisms with hyperbolic schemes passes `_hyperbolic_catalysis` at two
    subunits, and every one of 7,846 on mechanisms with an inhibitor, but only because the
    guard refuses the assignments that would fail it.
  - **Refused by a face check though the rate equals the law's.** Of the mechanisms of 8
    steps or fewer that a face check refuses, sampled 40 per population, these have a limit
    along corner directions under some consistent reading: 9 of 40 `:OnlyA` taggings of the
    uni-uni, bi-bi and ping-pong starting mechanisms (an estimated 127 of 566), 16 of 40 of
    the admission tests' mechanisms (47 of 118), 7 of 40 ter-bi (16 of 89), 4 of the 8
    refused ping-pong taggings (the gate's 25 bi-bi ping-pongs with a covalent form, each
    with both halves' chemistry `:OnlyA` and with every step `:OnlyA`; all 4 tag every
    step), 1 of 40 bi-ter (2 of 86), and none of the `:OnlyI`, ter-ter, bi-bi split and
    inhibitor ones. The guard reports each with its first consistent reading's `_PUMPED` or
    `_FILLED`. Nearly all are ping-pongs. Under the reading the oracles accept, the guard
    refuses 8 of the 9 `:OnlyA` taggings by the filling check and 1 as a bridge, 8 times at
    B's binding to the covalent form and once at P's. On such a face the cut-off forms hold
    most of the enzyme, yet the rate equals the law's because the inactive conformation's
    kept chemistry has the active one's constants. No simple rule found separates these from
    the wrong ones.
- **The guard's face checks under a searched reading rest on numerical evidence.**
  The filling check and the set-aside rule (Developer page, "Well-posedness guard") have
  arguments, not proofs, that the face condition with them suffices under every reading. A
  wrong law the guard admitted would not show in a fit: it is a smooth function that can fit
  well, with a wrong mechanism label and unreliable extrapolation toward the face where it
  fails. The protection is a gate run once on the guard as built. It ran the guard on 78,554
  mechanisms in ten populations: `:OnlyA` taggings of the uni-uni, bi-bi and ping-pong
  starting mechanisms (9,616), taggings with an `:OnlyI` group (43,311), ter-ter taggings
  (2,000), taggings of the bi-bi split children with the chemistry and some bindings
  `:OnlyA` (9,878), ping-pong taggings (50), bi-bi with a competitive inhibitor, ter-bi and
  bi-ter taggings (3,000 each), one expansion round's allosteric children (1,990) and the
  admission tests' mechanisms (2,709). None of
  the 29,517 the built readings admit is rejected. It ran the ε-path and RE-speed oracles on
  the 1,575 that only the search or the set-aside rule admits, all of them of 8 steps or
  fewer and samples of the rest: all 591 small ones agree along corner directions, 598
  larger ones agree along standard directions, and 100 ter-ter ones agree along corner
  directions; at two subunits the steady state equals the law for 40 of 40 sampled. With the
  filling check off, the guard admits 96 small `:OnlyA` taggings more (138 in all), and 18
  of them fail; standard directions miss these, so corner directions are the test.
  `test/Derivation/test_allosteric_rate_equation_correctness_guard.jl` keeps a sampled
  population test with corner directions. No hand-written test shows a one-way fast
  vanishing transition that passes the filling check because its target drains fast, along
  rapid-equilibrium steps that do not vanish, into the law's forms (the check's third
  clause); only the sampled test covers it. The gate's scripts, logs and report
  (`report.md`) are on the maintainer's machine, outside the repository, in
  `claude-scratchpad/lean-mechanism-spike/gate/`.
- **An allosteric ping-pong whose covalent forms empty only by binding B cannot stop both
  half-reactions.**
  PDH and OGDH are ping-pong and allosteric, so allosteric ping-pong must stay in the
  search. On 10 of the 39 ping-pong starting mechanisms of the bi-bi reaction
  `A[CX] + B[N] → P[C] + Q[NX]`, the guard rejects the version that tags both
  half-reactions `:OnlyA` and nothing else (3 as `_FILLED`, 7 as `_BRIDGE`), and the test
  suite's oracle finds no mass-action limit under any reading, at zero B and at zero A and
  B. In all ten the first half-reaction is at steady state and the second half's chemistry
  is a rapid-equilibrium step, and the enumerator keeps the version with the first half
  alone `:OnlyA`. Stopping both halves there needs the conformations to interconvert at the
  covalent form as well, with an equilibrium constant of its own and the multi-subunit rate
  formula rederived for two such forms; that belongs with iso mechanisms.
- **Some `:NonequalAI` mechanisms fit a constant the data cannot determine at one subunit.**
  The search tags a set of kinetic groups `:NonequalAI` when that raises the fitted count
  (`Derivation.fitted_parameter_count`), which counts constants the thermodynamic
  equations leave free, not constants the data determine. On the six-step bi-bi ping-pong
  with both half-reactions at steady state, the bindings of A and Q tagged together, or of
  B and P, fit two constants more than the plain mechanism; the data determine none more at
  one subunit and one more at two (`_testhelper_identifiable_rank`, `test/testhelpers.jl`).
  The L-phantom certificate (`_phantom_weights`,
  `src/Derivation/rate_equation_dedup_key.jl`) covers only mechanisms whose inactive
  conformation turns nothing over, so identify fits these. How many allosteric versions of
  the starting mechanisms are like them was not counted.
- **Some `:OnlyI` inhibitors fit a constant the data cannot determine at one subunit.**
  `_expand_to_allosteric` (the tagging move, `src/Enumeration/expansion_moves.jl`) can tag
  a competitive inhibitor's group `:OnlyI`, so that only the inactive conformation binds
  it. Of 122 sampled such children with no `:NonequalAI` group, of the bi-bi starting
  mechanisms of `A[CX] + B[N] → P[C] + Q[NX]` with a competitive inhibitor `I`, 47 fit more
  constants than the data determine at one subunit (2 of them two more), and 3 at two
  subunits. Where no group is `:OnlyA`, in 4 of 5 sampled, the rate equation is the plain
  mechanism's with the inhibitor's constant `K` replaced by `K·(1 + L)/L`. The L-phantom
  certificate covers 4 of the 47, whose inhibitor's terms lie in the inactive
  conformation's denominator alone, so identify fits 43 with a constant the data cannot
  determine. With the inhibitor's group `:OnlyA` instead, 28 of the same 122 fit one
  constant too many, and 19 of those pass the certificate.
- **Some `:NonequalAI` versions are reached only through an L phantom.**
  Where a set of kinetic groups gains a fitted constant tagged `:NonequalAI` and the
  enumerator also returns it tagged `:OnlyA`, the tagging move emits only the `:OnlyA`
  child, and `_expand_change_allo_state` (the relaxing move,
  `src/Enumeration/expansion_moves.jl`) reaches the `:NonequalAI` version from it. At one
  subunit that child is an L phantom for 227 of the 539 such sets on the 239 sequential
  bi-bi starting mechanisms of `A[CX] + B[N] → P[C] + Q[NX]`, and for 11 of the 124 on the
  39 ping-pong ones. An L phantom reads its parent's fit, so it enters the beam only
  through the beam's minimum width, and the `:NonequalAI` version is fitted only if it does.
  Above one subunit there are no L phantoms. Emitting these sets tagged `:NonequalAI` from
  the tagging move as well would end this, for a median of one more child per sequential
  mechanism and five per ping-pong one.
- **The flip move skips allosteric flank flips that the data can see.**
  `_expand_re_to_ss` never flips a flank of a steady-state chain (`_chain_flank_groups`,
  `src/Enumeration/expansion_moves.jl`), in allosteric mechanisms as in plain ones. For a
  plain mechanism a proof shows the flip adds no rate equation; for an allosteric one the
  rule is a modeling choice. Of the 278 bi-bi starting mechanisms of
  `A[CX] + B[N] → P[C] + Q[NX]`, 15 hold a chain. On four allosteric versions of each, a
  flip of flanks alone adds a constant the data determine (`_testhelper_identifiable_rank`)
  in 20 of 96 cases at two subunits and in 2 of 128 at one; the two at one subunit come
  from a version with an `:OnlyA` and a `:NonequalAI` group. So the search misses a slow
  binding beside a slow isomerization where the two conformations treat the chain
  differently.
- **A steady-state step beside an `:OnlyI` inhibitor can fit two constants for one ratio.**
  `_expand_split_kinetic_group` can give a step between two of a competitive inhibitor's
  complexes a kinetic group of its own, and `_expand_re_to_ss` can then flip that group to
  steady state. With the inhibitor's group `:OnlyI` the step runs in the inactive
  conformation alone, and where that conformation turns nothing over, its two rate
  constants enter the rate equation only as their ratio. The check that every steady-state
  group carries net flux (`_assert_emission_rules`, `src/Enumeration/expansion_moves.jl`)
  reads the whole catalytic scheme and does not see it. How often the search builds such a
  mechanism was not counted.
- **Requiring a competitive inhibitor and an allosteric regulator together makes 236,444
  starting mechanisms.**
  On the bi-bi reaction `A[CX] + B[N] → P[C] + Q[NX]` at oligomeric states 1 and 2, with a
  competitive inhibitor `I` and an allosteric regulator `R` both required,
  `candidate_mechanisms(rxn)` returns 236,444 starting mechanisms in 173 s, peaking at
  2.29 GiB resident, against 69,178 in 36 s and 1.19 GiB at 0.11.0; identify fits all of
  them in its first batch. With `R` alone required it returns 14,078 (6,428 at 0.11.0), and
  with `I` alone 1,548 (unchanged). The tagging move emits every set of `:OnlyA` groups and
  gives an inhibitor's group three states, and the dead-end move gives a new inhibitor
  three states on an allosteric parent. Marking either regulator optional
  (`optional_competitive_inhibitors`, `optional_allosteric_regulators`) avoids the cost.
  Measured once each in a fresh process (compilation included) on the maintainer's
  machine, Julia 1.13.4, `--heap-size-hint=2500M`.
- **One parent's allosteric children double with each kinetic group.**
  `_expand_to_allosteric` (the tagging move, `src/Enumeration/expansion_moves.jl`) emits
  every set of `:OnlyA` groups at each oligomeric state, `2^g − 1` for `g` kinetic groups,
  and the design measured how many of them the enumerator returns only up to seven groups
  (90 at seven, on ter-bi). Take the ordered ter-ter mechanism of
  `A[X] + B[Y] + C[Z] → P[X] + Q[Y] + R[Z]` with the abortive complexes EAR, EAQR and
  EABR, every step at rapid equilibrium but the chemistry: ten groups. At oligomeric
  states 1 and 2 its children (`candidate_mechanisms(rxn, parent)`) number 1,838 and take
  0.4 s, against 1,026 at 0.11.0. With a competitive inhibitor bound (eleven groups) they
  number 5,488 and take 3.0 s, and the call allocates 21.4 GB in passing, against 2,052,
  0.9 s and 6.3 GB at 0.11.0. The design forbids a bound on these children: any bound
  would hide a feasible version, and the argument that the tagging and relaxing moves
  reach every irreducible version starts from every set of `:OnlyA` groups. The relaxing
  move's search for minimal gaining sets (`_expand_change_allo_state`) grows the same way
  where `:OnlyA` groups cut others off from free E: tagging such a group `:NonequalAI`
  gains nothing, so the search extends every set that holds it. It tries 48 of 63 sets at
  six groups and 768 of 1,023 at ten (40 ms). At that rate it is harmless to about 15
  groups, an estimate and not a measurement, and exponential beyond. Measured on the
  maintainer's machine, Julia 1.13.4, `--heap-size-hint=2500M`; the times of the children
  are of a second call in one process.
- **No test shows the inhibitor checks deciding a `:NonequalAI` search.**
  Inside its search for `:NonequalAI` sets, `_expand_to_allosteric` counts a set whose child
  binds a competitive inhibitor in neither conformation or redundantly
  (`_inhibitors_bind`, `_redundant_copy_groups`) as not gaining and extends it, and asks the
  same of the set tagged `:OnlyA`. No known mechanism of six or fewer kinetic groups makes
  these checks change which children the search emits, so no testset pins that path.
  `_expand_change_allo_state` likewise counts a set whose child binds an inhibitor
  redundantly as not gaining and extends it. The testset "Allosteric mechanism — no
  relaxation makes a kept copy redundant" (`test/Enumeration/test_expansion_moves.jl`)
  pins such a refusal, but no testset shows a larger set emitted after one. The
  reachability testset ("to_allosteric and change_allo_state reach every irreducible
  version", same file) would show a version either search hides on its samples.
- **An allosteric mechanism is reached only through a plain mechanism the search expands.**
  `identify_rate_equation` fits the plain starting mechanisms first and expands the best of
  them, and a mechanism becomes allosteric only as a child of an expanded plain one
  (`_expand_to_allosteric`), unless the reaction requires an allosteric regulator. Data
  from an allosteric enzyme fit no plain mechanism well, so the plain mechanisms with the
  fewest constants win, and the true mechanism's plain parent can stay unexpanded.
  Measured once, on the bi-bi reaction `A[CX] + B[N] → P[C] + Q[NX]` at two subunits: 84
  noiseless rates generated from the six-step ping-pong with both half-reactions at steady
  state and the bindings of A and B `:NonequalAI` (9 fitted constants, of which the data
  determine 8), searched with `min_beam_width = 10` and `max_param_count = 9`. Of the 278
  starting mechanisms the search expanded 10, all with 5 fitted constants (losses 0.030 to
  0.045). None of the 64 with 6 constants or the 25 with 7 was expanded: their best loss
  is 0.033. The generator's plain parent has 7 constants and a loss of 0.042, ninth of
  those 25, so the generator was never built. The search still selected a mechanism that
  fits the data to a loss of 3.7e-13: an 8-constant ping-pong in which each of A's and
  B's bindings is one step with its chemistry, those two steps `:NonequalAI`. So this run
  recovered the rate equation and which steps differ between the conformations, and not
  the mechanism. Whether a rate equation can be missed this way is not known. A second
  run, with all chemistry `:OnlyA` and a required regulator, selected its generator, which
  is a starting mechanism there.
- **A batch's progress-log header can print twice.**
  When an L phantom's source's fit throws in the phantom's own batch, `run_batch!`'s log
  closure (`src/identify_rate_equation.jl`) prints the batch header again above the demotion
  lines, and the pre-fit summary has already counted those phantoms as inherited. Rows and
  CSVs are right. A phantom's source is usually fitted a batch earlier, so no identify run
  in the suite reaches this path; the fix (print the header before `_process_batch` and let
  `log` only indent) needs an identify fixture that does.
- **A regulator measured only at saturating levels gets an arbitrary K.**
  With an activator at 0 and 100× K, restarts that fit perfectly reported 0.015× to 2× the
  true K. The data do not constrain K there, so the fit returns any value in the box and
  lists it in `at_bound` only when it lands at an edge.
- **Where their data fall in memory slows the compiled law and `loss!`.**
  Measured on aarch64 (a Lima VM on Apple Silicon, 4 KB pages); each placement is fixed for
  the life of a process or of a `FittingProblem`.
  - The compiled law takes a call's concentrations and parameters by pointer. Called alone,
    it runs up to 7 ns slower per call (18 ns under coverage counters) where that block
    straddles a 4 KB page on the heap, and up to 17 ns slower (29 ns under coverage) at
    some stack offsets whose cause is not known. That touches `_reference_scale` (once per
    fit), `plot_rate_equation_on_data` and user code. The law's speed gate runs each of its
    five batches at its own depth with its own copy of the inputs, which removed every
    placement effect measured on aarch64: no reading came 5 ns over the law's median in
    any scan, under coverage included.
  - `loss!` inlines the law and puts no row on the stack. In two scans, each one process
    stepping through 48 stack depths for 12 laws, 3 of 1,152 readings ran over 1.5 times
    the law's median (HK at three adjacent depths, in one scan only), against 23 with the
    law called once per row. Where a problem's arrays fall on the heap still matters: of
    300 fresh problems of each of four small laws, 7-26% ran `loss!` over 1.15 times the
    law's median and 2-10% over 1.5 times, the worst 1.8-2.2 times (with the law called
    once per row, 6-27%, 5-25% and 2.6-4.6 times). The cause is not known.
  - Over 96 placements of whole processes (environment pads under `setarch -R`), 97 of
    1,728 warm fits of 18 laws ran over 1.15 times their law's median and 5 over 1.5 times,
    the worst 1.71 times, and the worst placement's 18 fits took 1.053 times the median
    placement's. With the law called once per row the same scan gave 75, 46, 4.65 and 1.81.
  On the slowest class of CI runner the law's gate failed in two jobs of five, each on one
  law, 22 and 31 ns over that law's usual time; those failures resemble the first effect,
  and their cause on x86 is not measured. Those jobs timed the law in the test process,
  under coverage counters and `--check-bounds=yes`. The gate reads its times in a child
  process without either (Developer page, What a fit costs). Without coverage the effects
  of the first bullet measured 7 and 17 ns rather than 18 and 29, in scans in one process
  on aarch64, not in the child; whether reading the times in the child ends the failures
  on CI is not known. Not measured: x86 and 16 KB pages (macOS on Apple Silicon). The
  investigations, with their scripts and scans, are on the maintainer's machine:
  `claude-scratchpad/loss-timing-investigation/`, and `lawgate/`, `inline-memory/`,
  `compiled-loss/`, `loss-two-loops/`, `loss-guard/` and `nocov-gate/` in
  `claude-scratchpad/lean-mechanism-spike/`.
- **The fit-cost pin failed once in a full suite.**
  "A fit costs about its loss! calls" (`test/Fitting/test_fit_rate_equation.jl`) failed once
  in a full suite on the allosteric bi-bi mechanism, at 0.387 s against 0.265 s allowed; its
  focused run gave 0.185 s against 0.253 s and the next full suite passed. With the law
  called once per row (5e1da4e5), the 96-placement scan of the previous entry reproduced
  that signature: the test failed at 2 of 96 placements, where the fit took 3.4-4.0 times
  evaluations × the standalone `loss!`, which read its usual 5.1 or 7.0 µs. With the
  present `loss!` it fails at none of 96. A fit there takes 1.23-1.46 times evaluations ×
  `loss!` (medians 1.36 and 1.32 on the two laws) against the bound of 1.5 times plus 3 ms
  per restart: the solver's 1.0-1.4 µs per evaluation weighs more beside a faster `loss!`.
  In one CI run, under coverage counters, the two fits took 0.67-0.81 of what the test
  allows in seven readings of eight and 0.96 in the eighth: the ordered bi-bi on an Ubuntu
  runner where `loss!` took 10 µs per call and the fit 1.53 times evaluations × `loss!`.
  A faster runner would fail the test with nothing wrong.
- **A pumped-basin message can say "nothing runs back" of a feed that runs back.**
  The guard's `_PUMPED` verdict names a feed into a stranded basin. It prefers a feed whose
  reverse is dead; when every feed has a live reverse and enters from outside the basin and
  the law network, it names one of those, and the message still says "while nothing runs
  back the way it came" (`_admission_message`), which is then false.
- **The guard's offered fix follows step order.**
  For a rapid-equilibrium step that its conformation cannot switch off, the message offers
  to tag the steps that join its sides, the path `_kept_path` finds first breadth-first in
  step order, and offers the fix only when that path's retag passes. It tries no other
  path. In `E + S ⇌ ES`, `ES <--> E + P`, `E + I ⇌ EI` (:OnlyI), `EI + S ⇌ EIS`,
  `ES + I ⇌ EIS` neither the mechanism nor its mirror image is offered a fix, though the
  path through the P release passes for the mirror
  (`test/Derivation/test_allosteric_rate_equation_correctness_guard.jl`, "the one rule that
  tells the states apart").
- **A rejection reports the fault of the first consistent reading in search order.**
  When no built reading is consistent and the consistent readings fail in different ways,
  the guard reports the verdict of the first one its search tries (`_guard_with_reading`),
  and the search walks the excluded steps in stored step order (`_consistent_columns`), so
  the fault and face a rejection reports depend on that order. The step a verdict names
  follows the stored step order too and, for `_PUMPED`, the rule that prefers a feed whose
  reverse is dead. Which mechanisms the guard admits does not depend on the order. The
  admitting reading could, since the guard admits under the first passing reading in
  search order, but on the gate's 78,554 mechanisms it does not: the guard at commit
  56bea2e6, which stores the steps in another order and differs in code too (no
  `_SHARED_NAME` check, another `_PUMPED` feed rule), gives them the same admissions and
  admitting readings as the package's (`report-lean.md`, beside `report.md`). Nine ter-bi
  rejections report `_FILLED` under the package's guard where that one reports `_PUMPED`,
  and 31 rejections, the nine among them, report another face, because the first
  consistent reading differs; for the nine and for 30 of the 31, every consistent reading
  gets the same verdict from both guards. 800 rejections name another step, 279 of them by
  the feed rule.
- **The search solves each candidate twice.**
  Under identify's `max_param_count`, `candidate_mechanisms` solves every candidate at its
  exit to count its fitted constants (`Enumeration._within_size_limits`), and identify's
  first pass solves it again to derive its law; for the starting set those counts run
  serially on the master. A count that throws on one child fails its parent's whole
  expansion, which identify records as one failure (`_expand_parents`). The reaction's
  name-decodability error (`_check_decodable_names`) fires after identify, with
  `show_progress` on (the default), has created `save_dir` and its `progress.log`.
- **The fitting tutorial shows a `FittingProblem`'s uninitialised buffer.**
  The block of `docs/src/fitting/tutorial.md` that builds `fp` ends with the assignment,
  so the page prints the problem with Julia's default `show`, which includes
  `log_ratios_buffer`, allocated `undef`, at whatever values its memory held.
- **The two guard files take names from each other.**
  `allosteric_rate_equation_correctness_guard.jl` takes 31 names from
  `rate_equation_correctness_guard.jl`: the verdict, its faults and the limit verdict
  (`_Verdict`, `_WELL_POSED`, `_NO_FAULT`, `_ARC_CYCLE`, `_BRIDGE`, `_EPS_NONUNIQUE`,
  `_FAST_STEP_PINNED`, `_FILLED`, `_ISLAND`, `_NO_ACTIVE_TURNOVER`, `_PUMPED`,
  `_TOO_MANY_READINGS`, `_LIMIT_UNTESTED`), the search caps
  (`_MAX_CHOICES_PER_CONFORMATION`, `_MAX_READING_PAIRS`), the networks and their graph
  walks (`_Network`, `_eps_network`, `_law_network`, `_node`, `_form_of`, `_copy_of`,
  `_out_transitions`, `_successors`, `_reached`, `_strong_components`, `_closed_classes`,
  `_aggregated_chain`), the faces (`_faces`, `_face_phrase`), `_guard` and
  `_steady_state_passes`. The guard takes eight back: `_SlowUptakeReading`,
  `_built_reading`, `_consistent_columns`, `_admission_interior`, `_admission_face`,
  `_active_turnover_verdict`, `_admission_message` and `_conformation_suffix`. Those eight
  are the use that would go: the allosteric guard builds on the guard's verdicts and
  networks, while the guard reaches admission only from its reading search
  (`_guard_with_reading`, `_verdict_under`, `_faces_verdict`) and its messages
  (`_guard_message`, `_face_message`). Merging the two files (2,551 lines), or cutting them
  so that the guard's entry, its reading search and the messages of both files, which
  rerun the guard to check an offered fix, sit above admission, would end the loop.
  `test/test_boundaries.jl` sets the guard's use of the allosteric guard aside until then.
- **The guard and the thermodynamic solve take names from each other.**
  `rate_equation_correctness_guard.jl` takes `_shared_constant_names` from
  `thermodynamic_constraints.jl`, to reject two constants with one name, and
  `thermodynamic_constraints.jl` takes `check_valid_rate_equation` from the guard:
  `parameters`, `parameter_units` and `dependent_parameter_values` call it to refuse a law
  that is not defined. That use is the one that would go: moving those three functions to a
  file that uses both would end the loop; moving `_shared_constant_names` alone would not,
  since it reads `_constant_columns` and `_regulator_names`, which the solve uses too.
  `test/test_boundaries.jl` sets the solve's use of the guard aside until then.

## Planned

- **Shrink the module interfaces.**
  The package aims at deep modules: few public names, simple arguments, the complexity
  inside (`.claude/CLAUDE.md`). The lean mechanism took `Mechanisms` from 25 public names
  to 15, and one function for both conformations' steps to 14; it took `Derivation` from
  16 to 9 and `Enumeration` to two, `CandidateMechanism` and `candidate_mechanisms`. Seven
  places still fall short. Priority: high.
  - **`has_nonzero_weight_cycle` takes six vectors.**
    `Mechanisms.has_nonzero_weight_cycle(from, to, weight, on, parent, offset)` is a
    union-find whose two callers outside `Mechanisms` pass the working arrays it fills. The
    enumerator's seed screen (`_seed_candidate_screen`,
    `src/Enumeration/starting_mechanisms.jl`) asks the turnover question for many step
    subsets and reuses its buffers. The guard calls it at one place, on its own network of
    transitions with a metabolite's uptake as the weight (`_potential_on_sccs`,
    `src/Derivation/allosteric_rate_equation_correctness_guard.jl`). A signature on a
    catalytic scheme and a mask, as `has_turnover_cycle` has, would fit the seed screen but
    not the guard, which passes transitions, not a scheme's steps. The guard also has its
    own walk over that network that finds a cycle whose offsets disagree (`_fast_segments`
    with `_tree_cycle`, `src/Derivation/rate_equation_correctness_guard.jl`), which its
    check that every cycle balances each metabolite uses (`_unbalanced_cycle`). Unverified:
    whether `_potential_on_sccs` could ask its question through that walk. If it can, the
    seed screen is the only caller outside `Mechanisms`, and with buffers passed to
    `has_turnover_cycle`, or an allocation accepted there, the routine becomes private: one
    public name and its six-argument signature fewer. Considered and not recommended:
    merging the two functions under one name (two call shapes remain); dropping
    `has_turnover_cycle` (the rule that turnover is uptake of the first substrate would then
    be written at four call sites); dropping the general routine outright (the guard would
    need a copy).
  - **Three modules count which way a step runs.**
    The count of substrates taken up and products given off that orients a step is made
    in `Mechanisms` (`_orientation`, the second row of the table of stored directions),
    again in `Derivation` (`_priority_key`, "stored against the reaction") and a third time
    in `Enumeration` (`_uptake_weight`). One owner needs a fifteenth public name in
    `Mechanisms`, or `_priority_key` reading the direction another way, which moves fitted
    names.
  - **Small facts of a mechanism are read in several modules.**
    Whether a mechanism is allosteric: `_is_allosteric` in `Derivation` and in
    `Enumeration`, inline in identify's `_fit_order` and in `Mechanisms`. A scheme's
    kinetic groups: `Mechanisms._groups` and `_kinetic_groups` in `Derivation` and in
    `Enumeration`. Its rapid-equilibrium flags: `Derivation._re_flags`,
    `Enumeration._rapid_equilibrium_flags` and the default argument of
    `rapid_equilibrium_segments`. An isomerization: `_is_iso` in `Mechanisms` and in
    `Enumeration`; a binding's shape is written four times in `Derivation` and as
    `Enumeration._is_binding`. The steps a conformation excludes, those of its groups'
    excluded allosteric state: `Mechanisms._steps_in_conformation`
    (`src/Mechanisms/mechanism_type.jl`) and `Derivation._excluded_steps`
    (`src/Derivation/Derivation.jl`). A mechanism's id, the first 16 hex digits of the
    SHA-256 of its text: identify's `_mechanism_id` and the capped key of `_canonical`
    (`src/Derivation/rate_equation_dedup_key.jl`); keying a capped law by the whole text
    instead would change every capped key once. Each copy is a line or two; one owner of
    each would cost a public name.
  - **Store allosteric states as what runs.** Deferred by the maintainer on 2026-11-22.
    `steps_running_in_allosteric_conformations` derives which steps each conformation runs
    from the written states and from which forms the running steps join to free E: a step
    between two cut-off forms does not run, whatever its state. The constructor could
    store states that are already true of what runs (a cut-off `:EqualAI` or `:NonequalAI`
    step stored as `:OnlyA` or `:OnlyI`), so that the masks are a read of the states, with
    no function shared across modules and no stored mask. Checked on one pair only, at
    oligomeric state 2: `E + S ⇌ ES :: OnlyA`, `ES <--> EX` and `EX <--> E + P :: OnlyA`,
    the steady-state step `ES <--> EX`, between two forms cut off in the inactive
    conformation, written `:: EqualAI` in one and `:: OnlyA` in the other. The guard
    admits both, the laws are equal where evaluated, the dedup keys (`key`, `fit_key`) are
    equal, and one constant's name differs (`k_ES_to_EX` against `k_A_ES_to_EX`). The design
    must answer three questions. States are stored per kinetic group, and a step is cut off
    on its own: in a group such as `(E + I ⇌ EI, ES + I ⇌ EIS) :: EqualAI` with `ES` cut
    off, one member runs and one does not, and both share constants. The written state
    carries intent the result does not: whether a step is cut off depends on other groups'
    states, so a move that reconnects the forms needs the written state; the enumerator
    would keep it in its record, which reopens the decision to keep proposing same-law
    taggings because their children differ; and a user who writes `:: EqualAI` would see
    `:: OnlyA` printed. The guard treats "excluded by its state" (`_excluded_steps`) and
    "does not run" in separate code, so its verdicts would need checking again on the gate's
    mechanisms. Storing whether each conformation turns the reaction over would serve one of
    the three callers of `has_turnover_cycle`, admission's `_active_turnover_verdict`; the
    guard's class 1 asks about every step, and the enumerator's maximal-rate test
    (`_has_vmax`) about a hypothetical subset of steps, so the function stays shared.
  - **`rate_equation_panels` takes a `Keq` it only passes on.**
    `Plotting.rate_equation_panels` hands `Keq` to `Data.RateData` and reads it nowhere
    else, and the method of `plot_rate_equation_on_data` for a rate function passes a
    placeholder `1.0`.
  - **`identify_rate_equation` has at least 18 keyword arguments,** such as
    `eq_complexity_filter = 337`. Each becomes an internal constant unless users have a
    real reason to change it, as they do for optimizer settings.
  - **The guard's internals carry names and arguments nothing needs.**
    `_ARC_CYCLE` and `_ISLAND` are never reported, so `_interior_verdict` could return one
    "inconsistent" outcome from a Boolean acyclicity test, and `_Fault`'s docstring could
    become a list. `exclusion_components` lacks the leading underscore. Of
    `_eps_network`'s callers, only `_every_step_cures` passes a `re` other than
    `_re_flags(m)`, and it could build its variant with `_with_steps` instead. The docs call
    the forms a conformation's kept steps join an island; the code calls them a component, a
    group and a part.
- **Decide what the guard keeps that no measured mechanism needs.**
  Whether each part below goes is the maintainer's decision: he chose the three built
  readings, and guard code does not go on a search alone. Priority: medium.
  - **The second and third built readings.** Once the guard searches the readings,
    `_ALL_STOPPED` changes no verdict in any of the gate's ten populations, and
    `_OUTWARD_LOSE_AFFINITY` was not measured. Measure both against the mechanisms the
    guard before 0.11.0 admits, those too large to search in particular, and remove what
    changes nothing: about 45 lines of
    `src/Derivation/allosteric_rate_equation_correctness_guard.jl` and their mirror in
    `test/testhelpers.jl`.
  - **Joint R1.** The branch of `_face_condition`
    (`src/Derivation/allosteric_rate_equation_correctness_guard.jl`) that refuses a bridge
    the two conformations' vanishing transitions open only together is reached through
    `_guard` by no mechanism searched: 20,096 variants of its fixture and the 54,927
    taggings of the gate's `:OnlyA`, `:OnlyI` and ter-ter populations.
- **Log the time of each search phase.**
  `progress.log` carries no timestamps, so a cluster run cannot say where its wall time
  goes. In the searches of 2026-11-02 the average worker was idle for 30% of PFKP's time and
  75% of PGK's (Developer page, What a fit costs). Workers connect within 15 s and the
  master's CPU covers under half of that idle time, so workers waiting for the slowest fit
  of a batch are the leading suspect. A timestamp on each line, and each batch's slowest
  fit, would rank all the search-speed ideas below. Priority: high.
- **Profile the search on the HPC runs.**
  In the searches of 2026-11-02, 64% of PFKP's and 75% of PGK's fits reach a loss within
  10⁻⁸ of another fit at the same parameter count. These may be laws the key misses, equal
  only after a reparametrization (`_canonical`,
  `src/Derivation/rate_equation_dedup_key.jl`), or laws that differ only outside the data.
  Separating the two would size a saving of up to 3–4×. Still unmeasured on these runs: fits
  whose parent fails the final parsimony cutoff, children and survivors per move, and how
  many PGK rows identify flags `L_phantom`. Priority: medium.
- **Tighten the speed gate's two bounds once CI's readings without coverage are in.**
  The bounds on the compiled law (50 ns + 0.5 ns per denominator term) and on `loss!` (1.25
  times the law plus eight `log` calls per row) were set from timings taken under coverage
  counters, and the gate reads its times in a child process without them (Developer page,
  What a fit costs). In ten such processes on aarch64 no law took more than 0.22 of its
  bound and `loss!` no more than 0.30 of what it is allowed, so both bounds pass a slowdown
  of three times or more. One CI run has been read. The worst law took 0.24 of its bound
  on Ubuntu, 0.35 on Windows, 0.43 on Ubuntu with Julia 1.11 and 0.50 on macOS; `loss!`
  took at most 0.12, 0.17, 0.16 and 0.29 of what it is allowed. macOS decides how tight
  either bound can go: there a 22-term law read 30 ns and three laws of 6 to 9 terms 12 to
  15 ns, against a median of 4 ns for laws of 16 terms or fewer, so half the law's bound
  would have left one law at 0.99 of it; and on 20 of the 48 laws the batch the gate
  reports for `loss!` was a slow one, with an allowance three to four times the usual.
  Read the gate's `@info` lines from a few more runs, then propose bounds that leave room
  for the slowest runner, or a statistic for `loss!` that such batches do not move
  (`_testhelper_loss_over_law` takes the best batch of each), for the maintainer's
  approval. The `loss!` guard read 0.16, 0.19, 0.24 and 0.42 against its bound of 0.7.
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
  `_expand_add_allosteric_regulator` give every regulatory site the mechanism's
  oligomeric state as its multiplicity. The derivation already powers each site by its
  own multiplicity, so only the DSL (`_parse_regulator_entries`,
  `src/DSL/enzyme_reaction_macro.jl`; `Regulator`, `src/Mechanisms/enzyme_reaction_type.jl`)
  and enumeration change: those two moves, and `_expand_merge_regulatory_sites`, which must
  merge only sites of equal multiplicity. Keep a bare `A` at the oligomeric state, widen
  only on an explicit list, and reject `I(m)` for a competitive inhibitor. Priority:
  medium-low.
- **Enumerate Cleland iso mechanisms.**
  The derivation handles iso mechanisms, but no enumeration move generates them, and each
  carries a dwell-time-gauge phantom unless the exact rank (see Ideas) or a
  "switch-product = 1" constraint row comes first. The enumerator's record of a form
  (`_FormRecord`, `src/Enumeration/candidate_mechanism_type.jl`) carries its conformation
  label, and `_requires_hyperbolic_catalysis`
  (`src/Enumeration/candidate_mechanism_type.jl`) is the extension hook. The enumerator
  builds the label `E` only, so until iso mechanisms are enumerated only tests reach the
  code of its moves for other labels, such as the split of a kinetic group by conformation.
  Two-conformation kinetic cooperativity stays parked: on PGK data no conformational model
  beat MWC with steady-state steps at equal parameter count, so reopen it only for a dataset
  whose curvature MWC cannot fit. Revisit allosteric ping-pong with this work: a flip at a
  second free form is the same question (see the Known issue "An allosteric ping-pong whose
  covalent forms empty only by binding B cannot stop both half-reactions"). Priority: low.

## Ideas

- **Remove the DSL's hints for the notation it no longer reads.**
  `@enzyme_mechanism` answers notation the lean mechanism dropped with a hint of what to
  write instead (`src/DSL/enzyme_mechanism_macro.jl`): a form in call notation such as
  `E(A, B)` (`_called_form_name` and the call branch of `_raw_term`), a reactant's modifier
  copy `X::Mod` and a count such as `2S` (the `::Mod` and count branches of `_raw_term`),
  and a tagged name in a name list such as `P::Mod` (`_bare_name`, `_untagged` and the hint
  of `_bare_symbols_from_values`). The maintainer kept them for this release. Removing them
  in a later one takes about 30 lines of source and 60 of
  `test/DSL/test_enzyme_mechanism_macro.jl`, and old notation then gets the generic error.
- **Extend `plot_rate_equation_on_data`.**
  Considered and left out: overlaying the fits of two mechanisms, a fit summary (measured
  against predicted, error per group), a tick rule for log axes, and a keyword to override
  the inferred x metabolite.
- **Spend less time compiling each law.**
  Compiling a fresh law takes 0.10-0.19 s, 9-12% of the PFKP and PGK searches' fit time and
  30-50% on PGK's laws of 6-8 parameters. Inlining the law into `loss!` adds 4-10% to a
  fresh small law's first `loss!`, law call and short fit, and 77-117% to a law of 396-858
  terms (aarch64, medians of 2 processes). The loops of `loss!` after its first do not
  depend on the mechanism, so moving them behind a function that compiles once may save
  part of this; the saving is not measured.
- **Reduce what CMA-ES itself costs per evaluation.**
  At 0.10.0, `loss!` no longer takes nearly all of a restart: an evaluation took 22 µs
  against 10 µs in `loss!` on a 20-parameter PFKP law, and 8.1 µs against 6.0 µs at the
  median of 60 PGK laws, so any faster rate law gains only on the `loss!` part. With the
  present `loss!`, the solver and the rest of the fit take 1.0-1.4 µs of 3.9-5.1 µs per
  evaluation on the fit-cost test's two laws (aarch64). The solver's eigendecomposition
  takes 5-14% of a warm fit, and the solver with the rest of the fit up to 20% of PFKP's.
  It belongs to CMAEvolutionStrategy, so this means a solver option or a change upstream;
  neither is explored.
- **Count identifiable parameters with an exact modular rank.**
  Several phantom classes make the beam count a model above its true dimension: chain
  flanks kept after a later split, copy-copy twin pairs, `:NonequalAI` copies and the
  iso dwell-time gauge. The one-subunit L has an exact test of its own: a rescaling
  certificate (`_phantom_weights`, `src/Derivation/rate_equation_dedup_key.jl`) finds an L
  phantom whose inactive conformation turns nothing over, and identify lets it read the fit
  of its law at L = 0 when that law holds its other fitted constants and a mechanism with it
  was fitted in an earlier batch or is fitted in the phantom's own.
  The rank of ∂v/∂log θ modulo a prime at random points of the derived law (`_law`,
  `src/Derivation/rate_equation_derivation.jl`) is exact up to a negligible failure
  probability, took 3.5 ms per mechanism in a prototype, and has a finite-difference oracle
  in `_testhelper_identifiable_rank`
  (`test/testhelpers.jl`). Run it over several
  saved runs to find the other phantom classes; performance may keep it out of `src/`,
  where it would otherwise replace the beam's parameter count. It must never filter
  children, since a phantom child can carry families its parent lacks.
- **Replace allosteric mechanisms at one subunit by a regulator that binds forms of the
  catalytic scheme.**
  At one catalytic subunit the equilibrium between the two conformations only reweights the
  enzyme forms, so the rate equation is a plain one written with `L` and two sets of
  constants. That is the common cause of the one-subunit Known issues: L phantoms, and
  `:NonequalAI` and `:OnlyI` versions with a constant the data cannot determine or reached
  only through an L phantom. The maintainer keeps allosteric mechanisms at one subunit: they
  are the only way a monomer binds an allosteric regulator, and the only allosteric versions
  of a catalytic scheme that is not hyperbolic. Their replacement would be a regulator that
  binds forms of the catalytic scheme and leaves them active with constants of their own, a
  family of mechanisms that needs a design of its own.
- **Choose the number of restarts per fit.**
  Over 120 sampled fits per enzyme, the best of 10 restarts took half the time of 20 and
  came within 1% of the best of 20 in 98% of PFKP's fits (worst 10.7% above it) and in all
  of PGK's. Stopping once the best loss recurs 4 times gives 1.7× and 2.4×, and 2% of PFKP's
  fits end more than 5% high. A loss above the beam's cutoff drops a mechanism for good, so
  test a rule by comparing whole search results, not single fits. More restarts have a case
  too: on five PFKP laws of 16 or more parameters whose fits had missed, 20, 40 and 80
  restarts gave 13, 3 and 0 misses in 60 fits, and one such miss dropped the line to PFKP's
  18-parameter law from the beam. No whole search has been run at 10 or at 40 restarts, and
  which way to move the default is undecided.
- **Start a child's restarts from its parent's fit.**
  A child usually holds its parent as a limit, so one restart could begin at the parent's
  parameters. Not measured.
- **Pin each law's operation counts.**
  The rate law's time bound tolerates the spread between machines, so undoing a 3× gain
  on a small law passes it: the MWC tetramer's former 20.3 ns sits under its 34 ns bound. A
  count of multiplications, additions and power calls per law does not depend on the machine
  and would catch this, at the cost of re-pinning when a law changes on purpose.
- **Emit `muladd` in `_horner`.**
  `muladd` made `loss!` 2–18% faster over seven sampled laws, and the generated code then
  breaks the binary `+` and `*` shape that the tests assert. It leaves the loss unchanged to
  12 digits, not to the bit. It was measured on aarch64 with a `loss!` that called the law
  once per row, and has not been measured with the present `loss!` or on the cluster's x86
  cores.
- **Compute the law's parameter-only coefficients once per `loss!` call.**
  A prototype splits the law into coefficients that depend only on the constants and
  `Keq`, computed once per call, and a polynomial in the concentrations for each row, in a
  vectorized row loop. Warm fits took 0.39 of the one-loop `loss!`'s time on a 96-term law
  and 0.20 on laws of 396 to 858 terms, where the present `loss!` gives 0.55 at 96 terms,
  0.48 at 396 and 0.39 at 858; small allosteric laws gain nothing. It costs a second
  generator, one new `Derivation` name (a method filling one prediction per row of a
  table), and a loss that changes in the last digits (at most 3.4e-15 relative). The
  prototype is `claude-scratchpad/lean-mechanism-spike/compiled-loss/` on the
  maintainer's machine.
  Separately, the scalar `log` is the largest per-row cost left for small laws (about 3 ns
  of 4 to 6); a vectorized `log` is untried.
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
  Derivation holds two exact row reductions, `_rref!` (the joint solve, in
  `Rational{Int}` retried in `Rational{BigInt}`) and `_rref_partition` (in
  `Rational{BigInt}`). `_rref_partition` serves only class 3's model point (`_model_point`,
  `src/Derivation/rate_equation_correctness_guard.jl`), whose integer null space
  (`_integer_nullspace`) goes through `_rational_nullspace`. Computing that null space with
  `_rref!` would leave one.
- **Decide admission's face condition exactly.**
  Admission's face condition (`_face_condition`,
  `src/Derivation/allosteric_rate_equation_correctness_guard.jl`), with the filling check
  and the set-aside rule, is sound in every test but not exact: it rejects mechanisms whose
  limit exists (Known issues, "The guard rejects some allosteric mechanisms whose law is the
  limit"). In an equilibrium bridge (R-B1) a tagged binding joins forms its conformation's
  other steps keep apart, yet those forms already sit at equilibrium, so the limit is
  deletion's (`bridge_deadend` in
  `test/Derivation/test_allosteric_rate_equation_correctness_guard.jl`); the set-aside rule
  admits such a bridge only when the law's own rapid-equilibrium steps hold the forms it
  touches in equilibrium. No exact narrower structural rule is known: the one candidate
  tried, T0++, admits faces the ε-path oracle rejects. The exact test compares
  Freidlin–Wentzell exponents over the admissible directions: the cut-off forms vanish when,
  for every rise of the tagged constants, every per-member rate and every order of the two
  limits, their cheapest in-tree costs more than the law network's; each comparison is one
  linear-programming feasibility test, not yet prototyped.
## Decided against

- **Scheduled CI and Downgrade runs.**
  A weekly or monthly run would catch a dependency release that breaks the package or stops
  the lower bounds resolving, with no push here. CI runs on pull requests and pushes to
  `main`, and `.github/workflows/Downgrade.yml` on pull requests that change `Project.toml`;
  the next such run surfaces either break.
- **A global ±ln(1e9) box, physical ranges or a plausibility warning.** Either needs a
  declared unit. `lb` and `ub` in `fit_rate_equation` stay, and the box follows the data.
- **`lb` and `ub` in `identify_rate_equation`.** The search sets each fit's box from its
  data.
- **Move dead-end decorations out of `init_mechanisms` into an expansion move.** It is not
  clear this saves time, since the move would have to run again on every child.
- **Separate the metabolite names inside enzyme-form names (`Ac + CoA ⇌ AcCoA`).**
  The enumerator rejects a reaction whose names spell one name two ways and names both
  spellings, and renaming one metabolite removes the collision, so no separator is needed.
- **Map BlackBoxOptim's convergence stop to `Success` upstream in OptimizationBase.** Only
  BlackBoxOptim users meet the warning, `test/Fitting/test_fit_rate_equation.jl` captures
  it, and the docs and HPC scripts use CMA-ES.
