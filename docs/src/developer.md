# Developer / Architecture

This page is for maintainers and contributors. It documents the internal
architecture of EnzymeRates.jl; the public API lives on the
[API Reference](@ref) page.

The package keeps two parallel mechanism representations, split along a
deliberate performance boundary. **Singleton types** drive the compile-time
rate-equation derivation; **concrete types** drive the mechanism enumeration.
The two halves and the optimizer layer are described below.

## Derivation architecture

A mechanism's rate equation is derived once, at compile time, by the
King–Altman/Cha method in `@generated` functions (`src/codegen.jl`).
To make that possible, a mechanism is encoded as a **singleton type** — the
exported `EnzymeMechanism{Sig}` or `AllostericEnzymeMechanism{Sig}`
(`src/types.jl`) — whose one type parameter, `Sig`, encodes the mechanism's canonical
reaction graph as a tuple of integers (`_sig_of`), so the order or direction in which
the steps were written does not change the resulting type. The compiler
specializes the derivation per type, so the symbolic King–Altman work runs
during precompilation rather than at every call.

This is the most important architectural decision in the package, and it is
deliberate. Moving the derivation to compile time leaves
`rate_equation(m, conc, params)` as straight-line arithmetic, written in Horner
form (`_horner`), that must be
**allocation-free and fast**: under 120 ns per call for every fixture mechanism, and
under `max(120 ns, 0.8 ns × denominator terms)` for the three largest enumerated bi-bi
laws (396 denominator terms each) and three ter-bi laws (589, 799 and 858 terms).
`_testhelper_test_rate_equation_performance` (`test/test_rate_eq_derivation.jl`)
measures each law's allocations and time per call, and the tests assert
`allocs == 0` and the law's bound. `rate_equation_string` (`src/print.jl`)
prints the law rather than the code: each denominator term over the unit rule's normalizer,
dimensionless whenever the law has one. `rate_equation`'s speed is the binding
constraint on the whole package: the fitter is a multi-start, global,
gradient-free optimizer that
evaluates `rate_equation` millions of times per fit, and a single rate equation
can take minutes to fit, so any per-call allocation or microsecond-scale
overhead would make fitting — and therefore `identify_rate_equation`, which fits
thousands of candidates — impractical.
`loss!` is held to the same standard: `FittingProblem` pre-allocates its
`log_ratios_buffer` once and `loss!` reuses it, so the inner optimization loop
allocates nothing. The cost is the flip side of the benefit: each unique
singleton type triggers a full symbolic derivation, and a mechanism with many
enzyme forms or steps can be slow to compile, exhaust memory, or `StackOverflow`.

`compile_mechanism` (internal) is the boundary between the two representations:
`compile_mechanism(m::Mechanism) = EnzymeMechanism(m)` and
`compile_mechanism(am::AllostericMechanism) = AllostericEnzymeMechanism(am)`. The
lifts back, `Mechanism(em)` and `AllostericMechanism(aem)`, decode the Sig
(`_mechanism_from_sig`, `_allosteric_mechanism_from_sig`) into the mechanism it encodes,
so a mechanism survives the round trip unchanged. Decoding trusts the Sig, which only
`_sig_of` writes: a plain Sig's decoding checks nothing, and an allosteric one's runs only
the composition constructor's tag and site checks. A mechanism of either representation
prints as the macro block that builds it (`show`,
`src/dsl.jl`), from the concrete mechanism, so the compiled types' `show` takes its
argument `@nospecialize` and `@nospecializeinfer` and compiles once. `_parse_mechanism`
reads that text back through the macros' walker and checks, never evaluating it, and
identify's CSVs store it with a short id (`_mechanism_id`).

Only the `@generated` methods (`rate_equation`, `parameters`, `fitted_params`,
`metabolites`) do work per singleton type: `rate_equation` runs the derivation,
`fitted_params` and `parameters` run its constraint solve, and `metabolites` reads names off
the lifted mechanism. The lifts read the Sig at run time, and the derivation helpers that
take a singleton or its type take it `@nospecialize`, so they compile once for all
mechanisms. Each mechanism type then pays for its generated bodies and the derivation
they run, and for the methods that specialize on it: `all_parameters`, Base's
`text/plain` display, `LogConcentrationTurnover`, and `FittingProblem` and `loss!` below.
`FittingProblem` and `loss!` specialize on the type on purpose: `loss!` then reads
`fitted_params` and `metabolites` as constants and calls the generated `rate_equation`
without dispatch. The optimizer never sees that type. `fit_rate_equation` takes its
`FittingProblem` `@nospecialize` and hands `Optimization.solve` a `FitObjective`, which
holds `x -> loss!(x, fp)` behind a `FunctionWrapper`, so every mechanism gives the solver
one problem type: `solve` and the solver compile once per optimizer, not once per mechanism.
The peak search that scales a relative fit (`_peak_turnover`) wraps the mechanism's turnover
the same way, so its search loop compiles once.

### The joint thermodynamic solve

`_solve` (`src/thermo.jl`) builds one integer system per mechanism
(`_joint_system`): a column per constant of every edge (`_Column`), with an inactive
copy for each edge of a `:NonequalAI`, `:OnlyA` or `:OnlyI` group; a row per fundamental
cycle and conformation, and a row tying each further member of a kinetic class to its first.
It orders the columns by `_priority_key` (`_ordered_columns`) and reduces the system
exactly (`_reduce_joint_system`: `Rational{Int}`, redone in `Rational{BigInt}` on
overflow). Pivot columns are dependents. `_solve`
returns them with the independents and the fitted names as a `_Solve`, which the engine,
`parameters`, `fitted_params`, `all_parameters`, code generation and the printer read.
Columns outside a conformation's law (an `:OnlyA` edge's inactive copy, an `:OnlyI` edge's
active copy, an edge on a form that conformation cannot reach) come first, so no kept
dependent reads them. Each class is named after its member
that comes last in that order (`_kept_members`), through the `name(p, m)` chokepoint. The
engine (`_law`, `src/engine.jl`) derives each conformation's
law on a mask of the edges it keeps (`_active_mask(am)`, `_inactive_mask(am)`; every edge
of a plain mechanism). Every edge keeps its stored direction, so
both conformations share the mechanism's orientation and names. The canonical form
decides the order and the direction it reads: free E is form 1 and the other forms follow
in content order (conformation label, bound metabolites, residual); an edge that neither
chemistry nor graph context orients runs from the earlier form; and a tied group's lead
turns the group's other tied edges to match it (`_align_tied_groups!`). The order's last
tie-break reads each edge's content as stored: its forms' places in content order, then
the metabolites it takes up and gives off (`_edge_base_rank`). The canonical form thus
decides the dependents and so the fitted names: a change to it changes fitted names.

### The law, its code, its text and its key

`src/poly.jl` holds sparse polynomials over integer-indexed variables (`_Poly`), with
rational exponents for constants and checked `Int64` coefficients: under the term cap no
coefficient leaves `Int64`, so an overflow is an internal error. `_law(m, s, mask)`
(`src/engine.jl`) derives one conformation's law, `_Law`, as `v = E_total·H·Q/D`. It
contracts each rapid-equilibrium segment to a node rooted at its first form (free E roots
its own), sums the spanning trees of the segment graph, builds the numerator from each
steady-state edge's net flux, substitutes the solve, divides out `H = Π[S] − Π[P]/Keq`
exactly and brings `Q` and `D` to lowest terms in the concentrations. An MWC law
(`_AlloLaw`) holds both conformations' laws, their free-E weights over a common factor and
the regulatory sites' polynomials. `_eq_complexity` counts V×τ, both conformations of an
MWC law together, by a checked Bareiss determinant; the cap (`_assert_derivable`) and
identify's `eq_complexity_filter` read it before any polynomial is built. `_Poly`, `_Law`
and `_AlloLaw` are mutable structs with `const` fields, never an immutable wrapper of one
reference: Julia 1.12 and 1.13 can hoist such a wrapper's box out of a loop that merges it
with another type, without a write barrier.

`_rate_body` (`src/codegen.jl`) emits the law as straight-line locals: the dependents in
solve-column order, each segment's weight and tree sums, `Q`, `H` and the MWC terms, all in
greedy Horner form with binary `+` and `*`. Each input is read by name into a positional
local, so a `NamedTuple` in any field order works, and reciprocals and roots (`sqrt`, or
`x^(one(x)/q)`) are computed once.

`rate_equation_string` divides the law by the unit rule's normalizer (`_normalizer`,
`src/normalize.jl`: the constant term, else Cleland's normalizer when the denominator holds
`ΠS/X` for every substrate `X` and `ΠS`, else a virtual constant term on the first
affinely independent terms), so every denominator term is dimensionless; a law with no such
normalizer prints as derived. With `coefficients = true` it prints the law in named
coefficients.

`_canonical` (`src/key.jl`) labels a law canonically: each part to lowest terms, the fitted
constants coloured by refinement on the term–constant incidence and individualized while
colours tie, the concentrations by name. The SHA-256 of its text is the law's key, by which
identify fits one mechanism per law and reads that fit through each mechanism's own names
(the canonical order). At 100 labellings the search stops, and the law takes its
mechanism's own key. The L-phantom certificate (`_phantom_weights`) solves exactly for a
rescaling of the other constants that absorbs `L` in a one-subunit variant whose inactive
conformation turns nothing over, and it alone decides an L phantom (`_is_phantom`). Such a
variant reads the fit of its active law when that law holds every fitted constant but `L`
and a mechanism with that law was fitted in an earlier batch or is fitted in the variant's
own (`_process_batch`).

## Well-posedness guard

`src/admission.jl` decides whether a mechanism's rate law is the unique fast-equilibrium
limit of mass-action kinetics wherever data sit. `_guard(m)` reads the reaction graph
through masks and flag overrides and never constructs a mechanism, so a move can judge a
candidate on its flags (`_faces_pass(m, re)`). It returns an isbits `_Verdict` (fault,
limit-test verdict, face, conformation, edge), never a `Union` that a loop reads: Julia
1.12 and 1.13 can hoist the box of a `Union` value out of a loop without a write barrier,
which corrupts the heap. `_check_well_posed(m)` throws `_guard_message(m, v)`, which
names each step by its DSL line (`_step_text`) and the face; `_is_well_posed(m)` is the
enumerator's drop test.

The guard builds two networks of one catalytic subunit whose conformations flip at free E
(`_Network`): the **law network**, each conformation's kept steps on the forms free E
reaches through them (`_law_mask`), and the **ε-network**, every step of both
conformations, each excluded step marked vanishing and an excluded rapid-equilibrium plain
binding kept fast (`_eps_network`). The faces (`_faces`) are bit masks over the metabolite
names. In order it checks class 1, a cycle that turns the reaction over
(`_has_turnover_cycle`); class 2, no fast cycle of the ε-network that takes up a name net;
admission's interior (`_admission_interior`: existence by components and arcs, islands,
the active conformation's turnover); then on each face class 3 (`_check_face`: closed
classes, bottomless segments, the aggregated chains, and a limit test modulo 2^61 − 1 where
structure decides nothing) and admission's face condition (`_admission_face`: the
ε-network's closed classes, stranded basins, and the fewest-slow-transitions in-tree test
per excluding conformation and for both together). Class 5, a denominator that keeps no
monomial on a face, runs on the derived law (`_assert_support`) wherever a law is
derived for use: the generated `rate_equation`, identify's PASS1 and
`rate_equation_string`. The macros, the parser
and the derivation entry throw; `init_mechanisms`, `seed_mechanisms`, `expand_mechanisms`
and identify's `_base_tier` drop; the flip move runs classes 2 and 3 on each candidate's
flags; the tagging and change-tag moves filter by admission alone (`_admits`) and leave
class 3 to the exit.

## Enumeration engine architecture

Mechanism enumeration uses the **concrete types** `Mechanism` and
`AllostericMechanism` (`src/types.jl`) — internal, not exported. A `Mechanism` is a
reaction graph in canonical form: its sorted `substrates` and `products`; `table`, each
metabolite entry its forms and steps hold, once, sorted by name and then by whether it
binds as a modifier (`Met`), a reactant and its modifier copy being two entries;
`conformations`, the sorted conformation labels; `forms` (`Form`: a label's index, and the
entries it binds, adds and subtracts, each a packed count of table entries, `Counts`);
`edges` (`Edge`: two form indices, the entries the step takes up and gives off, and its
rapid-equilibrium flag); and `group[k]`, edge `k`'s kinetic group. It holds no reaction,
and `==` and `hash` read its fields alone, alike in every process. An
`AllostericMechanism` holds a catalytic `Mechanism`, one allosteric tag per kinetic
group, the catalytic multiplicity and the regulatory sites.

The macros, the parser and the moves that change flags, groups or steps write a mechanism
as raw input (`RawMechanism` of `RawStep`s between `RawForm`s, and
`RawAllostericMechanism`): its steps in written order and direction. `Mechanism(raw)` and
`AllostericMechanism(raw)` build its canonical graph (`_canonicalize`) and check it. The
moves that tag a mechanism or change its regulatory sites compose the parent's catalytic
graph with other tags or sites through the composition constructor,
`AllostericMechanism(m, tags, multiplicity, sites)` (directly or by `_with`), which runs
only its tag and site checks. Going `from → to`, an edge takes
up the entries in `consumed` and gives off those in `released`. Steps are classified by
these counts. A binding takes up exactly one metabolite and gives off none
(`_is_binding`). It is plain when its target form is its source form with that metabolite
added, the residual unchanged and the conformation free to change (`_is_plain_binding`),
and fused otherwise, as when the last substrate binds straight into the product-bound
form (`E(A) + B → E(P, Q)`). The constructor rejects a binding whose forms match a plain
binding by name but not by role, such as `E + P::Mod ⇌ E(P)`. An isomerization takes up
and gives off nothing (`_is_iso`); a Theorell–Chance step takes up one metabolite and
gives off another. `_is_chemistry` is true for every edge but a plain binding; the
allosteric moves (`_chemistry_groups`), admission, the joint
solve's priority (`_priority_key`) and `_merged_complexes` read it to tell catalysis from
binding. The derivation reads every binding alike, so a fused binding's rapid-equilibrium
constant is a dissociation constant and its steady-state pair a binding and a release
rate constant. Every step's constants are named by the two sides of its kinetic class's
kept member — each side's enzyme form followed by its free metabolites (`K_ES_to_E_S` for
`E + S ⇌ E(S)`, `K_EPQ_to_EA_B` for `E(A) + B ⇌ E(P, Q)`, `k_EA_B_to_EQ_P` for
`E(A) + B <--> E(Q) + P`).

Like the singleton types, the graph is canonical, so the order or direction in which
steps are written does not change the resulting mechanism. `_canonicalize` stores free E
first and the other forms in content order, and stores every edge in its
physical-forward direction (`_orientation`): a step with one metabolite takes it up, so a
release is stored as the binding it reverses; otherwise the side holding more substrates,
then fewer products, comes first, then the form where only products enter or leave, and
failing those the earlier form. It sorts the edges and numbers the kinetic groups by
their first edge. The edges of one group that metabolite progression leaves tied, such as
a conformational change and its mirror at inhibitor-bound forms, turn as one
(`_align_tied_groups!`): each takes up and gives off the same metabolites, and each that
changes conformation runs between the same two conformations in the same order, so the
group's shared constants describe one physical direction. A tied edge that keeps its
conformation, such as one metabolite exchange at two conformations
(`E(I) + J <--> E(J) + I` and `F(I) + J <--> F(J) + I`), turns by its metabolites alone.
A group whose tied edges share no such direction is rejected. The `RegulatorySite`
constructor sorts a site's ligands by name, and the `AllostericMechanism` constructor
sorts the sites.

The constructors check representation and the `:OnlyI` rule only (`_check_declarations`,
`_check_representation`, and for an allosteric mechanism `_check_tags`,
`_check_onlyi_groups` and `_check_ligands`): names a metabolite can take and conformation
labels shaped as labels,
none a metabolite's name; every name the steps hold declared, and every declared name
taken up, given off or bound by some step; residuals that add substrates and subtract
products; the free form E (nothing bound, no residual), and edges that join every form to
it. They also enforce the kinetic-group rules. A group's edges take up and give off the
same entries and carry one RE/SS flag (`_check_uniform_groups`), so a fused and a plain
binding of one metabolite may share a group and a Theorell–Chance step cannot share one
with a binding. A reaction — the pair of an edge's two sides — appears in one edge
(`_check_each_reaction_once`). No two distinct forms render one name
(`_check_distinct_form_names`): a form's name joins its conformation and bound
metabolites without a separator, so E with a metabolite `star` bound and the free
conformation `Estar` would both be `Estar` and share their constants' names. Names that
spell one string two ways, such as `NAD`, `P` and `NADP`, never reach the constructor:
`@enzyme_reaction` and the mechanism macros reject them (`_check_uniquely_decodable`).
Whether a mechanism's rate law is defined is the well-posedness guard's question
(`src/admission.jl`; see [Well-posedness guard](@ref)). The macros, the parser and every
derivation entry throw on a mechanism it rejects (`_check_well_posed`), and identify's base
tier and the enumerator's exits drop one (`_is_well_posed`), except that init's exit runs
the guard on its seeds and checks class 1 on their variants, whose classes 2 and 3 the flip
search ran (`_faces_pass`). The constructors and the Sig decoder do not run it.

These are ordinary value types to avoid excessive precompilation costs. The enumeration builds,
expands, and deduplicates many thousands of candidate mechanisms (see
[The enumeration engine](@ref)). If each candidate were a singleton parametric
type, merely constructing it would trigger compiler specialization — the same
precompilation cost the derivation pays — multiplied across thousands of
mechanisms, which would be prohibitive. Carrying mechanism structure as runtime
values instead means enumeration and deduplication cost no compilation at all.
Only the candidates the search actually fits are lifted to singleton types, one
at a time, through `compile_mechanism`.

The moves never emit a child that is provably a reparameterization
of its parent. The split move divides a group by binding context and accepts a
set of splits only if the independent-parameter count rises; the count comes from
the thermodynamic constraint solve, and for a `Mechanism` it is evaluated without
building the child, from a cycle basis computed once per parent
(`_partition_independent_count`). For an `AllostericMechanism` it is the joint solve's
count on the built child, which derives no rate polynomial (`_split_gain_test`).
A split never changes admission's interior verdict, since both parts keep the group's tag,
and a child the guard rejects, inside or on a face, is dropped where `expand_mechanisms`
returns it. The search extends a failed set by the splits of the groups that share its
rapid-equilibrium segments. The RE→SS move flips whole groups and
accepts a set only if the rapid-equilibrium segment count rises and every flipped group
carries net flux in the child. In a `Mechanism` it never offers as a unit a flank
of a qualifying chain whose isomerization is steady state (`_chain_flank_groups`):
by the chain lemma that flip adds parameters and no rate law. An
`AllostericMechanism` keeps those flips, since over more than one catalytic
subunit a flank's flip can raise the rank. Flux is decided on the graph of
rapid-equilibrium segments (`_flux_carrying_groups`): each steady-state step is an
edge weighted by
its net uptake of substrates minus products, segment offsets included, and a step
carries flux exactly when its biconnected block holds a cycle of nonzero weight. A
steady-state group with no such step exposes only the ratio of its constants, so a
flip set that leaves one counts as failed and is extended, and a split part with
none is emitted at rapid equilibrium. A dead-end copy of a substrate or product
is redundant, and never emitted, when in every conformation where it binds each of
its complexes has a productive twin (a form that holds no modifier, in
the complex's rapid-equilibrium segment with its offsets, so the two weights are
proportional: `_productive_twin`, tested over a group by `_all_twin`; every
conformation holds its free enzyme, `_redundant_copy_groups`), and the dwell gauge,
which rescales each twin by
its own factor and all of the copy's complexes by one factor, is consistent with
every kinetic group's shared constants, the copy's own included, an `:EqualAI`
group taking one rescaling in both conformations (`_gauge_rescaling`, bookkeeping
over the edges whose classes, `_GaugeClass`, are a kinetic group's or a twin form's own,
per conformation); a copy whose complexes all have twins but whose gauge fails is
kept, since only the gauge proves the copy's constant invisible. The
dead-end move skips a redundant placement, the split drops every child that holds
a redundant copy group (a split can complete a gauge, and splitting a group that
forms a copy's twins can break one, so only the whole child decides), and
`_expand_change_allo_state` and the parent check reject a redundant group
(`_redundant_copy_groups`). A complex
that duplicates only a copy's complex never rejects, since the sites that pin
either copy may keep both constants separable.
`expand_mechanisms` asserts both rules on every parent (`_assert_emission_rules`).
Both refinement moves share one minimal-set search (`_minimal_gaining_sets`).
Duplicate laws that survive these proofs share one fit, an L phantom's key aside:
identify fits one mechanism per key (`_canonical`, `src/key.jl`). A numerical
identifiability rank exists only in the test suite, as an oracle for the proofs and the
L-phantom certificate; nothing in `src/` estimates identifiability numerically.

`init_mechanisms` appends to the seeds their merged and Theorell–Chance variants
(`_seed_variants`, plain `Mechanism`s only). `_merge_isomerization` removes a
chemistry isomerization and moves every other step at its substrate side onto its
product side, so the last substrate's binding becomes fused; `_eliminate_form`
replaces a form whose only two steps are bindings into it by one Theorell–Chance
step in a group of its own. Every step of a base starts at rapid equilibrium. The
variants are the inclusion-minimal sets of the base's flux-carrying groups
(`_minimal_flips`) whose flip to steady state leaves no rapid-equilibrium
turnover cycle (`_re_turnover_cycle`: each RE step weighted by its uptake of
substrates minus products, a turnover cycle being a cycle of nonzero weight),
keeps a maximal rate both ways (`_has_vmax`, condition V), keeps chemistry out of
equilibrium with both sides (`_chemistry_equilibrates_both_sides`, condition C),
passes the guard's face check and keeps every steady-state group flux-carrying.
`_seed_candidate_screen` runs the first three tests on index arrays with a
weighted union-find, so the other two run only on a candidate that passes them. A
merged variant whose every merged complex has two steady-state steps, each alone
in its group, has the family of the unmerged chain with rapid-equilibrium flanks
at the same count and is skipped. `_degenerate`, the failure of V or C on the
active-state graph, marks the starting mechanisms that `_base_tier`
(`src/identify_rate_equation.jl`) replaces by their flip children that are not
degenerate.

Conformational mechanism types declare `_requires_hyperbolic_catalysis` (true for
an `AllostericMechanism` whose catalytic multiplicity is above 1), and the moves
that can give such a mechanism a non-hyperbolic catalytic scheme consult
`_hyperbolic_catalysis`: `_expand_to_allosteric` promotes a non-hyperbolic
parent at multiplicity 1 only, and above one subunit the tagging and change-tag
moves (`_expand_to_allosteric`, `_expand_change_allo_state`), like the flip move
(`_expand_re_to_ss`), drop a child that fails the predicate. On an
`AllostericMechanism` the predicate reads each conformation's law (`_active_mask`,
`_inactive_mask`), since a conformation that runs part of the catalytic graph can
carry a power the whole graph lacks. The predicate scores
the rapid-equilibrium segment graph of the catalytic scheme, leaving out every
edge that touches a form holding a modifier, a reactant's modifier copy included: for
each substrate and product, a directed steady-state edge scores one if it binds the
metabolite in that direction plus the count the edge's source form carries beyond its
segment's bottom form (the offsets of `_segments`, which the split move's
`_bottomless_segment` also reads), and a segment scores the most any of its forms
carries. The equation's degree in the
metabolite is the best score over root segments and spanning arborescences toward
the root; the predicate decides "at most 1" with reachability checks
(`_all_reach`) rather than a max-arborescence solve, since a weight-2 pattern is
one segment or edge scoring 2, a scoring root plus a compatible scoring edge, or
two compatible scoring edges. The exactness of the predicate against the derived
denominator is pinned by a test over the uni-bi enumeration and the
ping-pong-capable bi-bi enumeration, which contains every bi-bi mechanism. A future
conformational type (KNF, mnemonic, slow isomerization) adds one
`_requires_hyperbolic_catalysis` method and calls the predicate in its promotion
move.

## Optimization algorithm architecture

Fitting reaches solvers only through Optimization.jl; the package ships no solver of
its own. `fit_rate_equation` (`src/fitting.jl`) builds an `OptimizationProblem`
whose objective, `_fit_loss`, scores each point with `loss!` through a
`FitObjective` (see [Derivation architecture](@ref)), and calls
`Optimization.solve` with whatever optimizer the caller passes. This gives the
package access to the global, gradient-free optimizers that non-convex
rate-equation fitting needs — CMA-ES (`OptimizationCMAEvolutionStrategy`) and
BBO differential evolution (`OptimizationBBO`) are the tested choices — without
taking on a solver dependency or locking users into one algorithm. Optimization.jl
common options (`maxtime`, `maxiters`, …) are named keyword arguments, and
solver-specific options pass through a `solver_kwargs` named tuple; see
[Loss & optimizers](@ref) for the user-facing view.
