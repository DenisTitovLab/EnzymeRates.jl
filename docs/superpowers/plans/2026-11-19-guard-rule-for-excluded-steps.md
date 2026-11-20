# Guard Rule for Excluded Steps Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the well-posedness guard read how an excluded step is switched off from the
reaction graph alone (stage 1 of the lean-mechanism design), on today's data model.

**Architecture:** One internal function in `Derivation` classifies each excluded step of a
conformation as losing its affinity, stopped, or impossible to switch off; the guard's
existing checks read that classification where they read `step_kinds` today. The
constructor's `:OnlyI` check goes, and `Enumeration` gains one filter. No public name
changes.

**Tech Stack:** Julia 1.11+, the package's own test oracles
(`test/mechanism_definitions_for_test_enzyme_derivation.jl`), Documenter doctests.

**Spec:** `docs/superpowers/specs/2026-11-18-lean-mechanism-design.md`, decision 2,
"What Enumeration does not build" in decision 5, "Tests" and "Documentation". Read those
before starting. Stages 2 and 3 of the spec are not in this plan.

## Global Constraints

- **The energy-level rule (Tasks 10 to 13), copied from the spec's decision 2.** The steps a
  conformation keeps join the forms into groups; every group but free E's is cut off and
  lies higher in energy.

  | Excluded step | How it is switched off |
  |:--|:--|
  | Between two groups | Loses its uphill rate and keeps its downhill rate. A rapid-equilibrium step stays fast. |
  | Between two groups, steady state | May be stopped instead: both rates vanish |
  | Inside one group, steady state | Stopped |
  | Inside one group, rapid equilibrium | Fault |

  A reading chooses, for each excluded step between two groups, which side is higher, or,
  for a steady-state step, that it stops. It is consistent when the steps that lose, each
  pointing uphill, close no cycle, never point into free E's group, and lead from free E's
  group to every group. The guard admits a mechanism when some consistent reading passes
  class 2 and the face checks. It tries the three built readings first, then searches the
  rest, unless the mechanism is too large to search. Read the spec's subsections "The
  search for a reading" and "The face checks" before Tasks 11 and 12: they state the
  filling check and the set-aside rule in full.
- **The three built readings (Tasks 2, 8 and 9).** In a conformation, for each step it
  excludes:

  | Excluded step | How it is switched off |
  |:--|:--|
  | Rapid equilibrium, takes up one metabolite into a cut-off form | Loses its affinity: stays fast, uptake vanishes |
  | Rapid equilibrium, anything else | Fault |
  | Steady state, takes up one metabolite into a cut-off form | Loses its affinity (uptake rate vanishes, release keeps its rate), or is stopped: whichever gives a limit |
  | Steady state, anything else | Stopped: both rates vanish |

  The guard tries three readings of the steady-state uptakes into cut-off forms (Tasks 8
  and 9), in this order: every one loses its affinity; one loses its affinity only when
  it leads away from free E, the others stopped; every one is stopped. It admits at the
  first try that finds no fault and otherwise reports the first try's verdict. No other
  mixture of the two readings.

  "Leads away from free E": the steps the conformation keeps join the forms into parts. A
  part's distance is the fewest excluded steps that take up one metabolite into a cut-off
  form, rapid-equilibrium or steady-state, each followed from its unbound side to its
  bound side, that lead to it from free E's part. A step leads away when its bound side's
  part is farther than its unbound side's; a part no such path reaches is infinitely far.

  "Takes up one metabolite into a cut-off form" means: `length(e.taken_up) == 1`,
  `isempty(e.given_off)`, the steps the conformation keeps do not join `e.from` and
  `e.to`, and they do not join `e.to` to free E (form 1).
- `Derivation` reads no `step_kinds` and no form contents for this rule. Other reads of
  `step_kinds` in `src/Derivation/` (names, thermo, engine) stay until stage 2.
- No `public` statement and no export changes. The module-map test must still pass.
- A test may change only in the four ways the spec lists under "Verdicts that change".
  The tests of a law that is not the limit pass unchanged.
- **Stop and report to the maintainer** when: a mechanism admitted before this work is
  rejected for any reason but the built guard's fault for a rapid-equilibrium step it
  cannot switch off (the maintainer accepted two exceptions, gate
  entries 53145 and 53246, rejected with `_BRIDGE`: each holds an excluded
  rapid-equilibrium step that binds and reacts at once, which the old guard read as
  stopped); a newly admitted law fails an oracle; a
  hand-written `MECHANISM_TEST_SPECS` entry or an identify test's true mechanism is now
  rejected; or a test would have to change in a fifth way.
- One Julia process at a time. Focused run of one test file (from `.claude/CLAUDE.md`):
  `julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<Module>/<file>.jl")'`.
  Full suite, in the background:
  `timeout 1200 julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`.
- Work on branch `guard-rule-excluded-steps`. Stage files by path; never `git add -A`
  (`simulations/` is untracked and not part of this work). Never skip a hook.
- Version 0.11.0, bumped once, in Task 7.
- Spike scripts and logs for reference only:
  `/home/denis.linux/claude-scratchpad/lean-mechanism-spike/`. They patch a running
  session; do not copy their patching into `src/` or `test/`.

## Review Focus

Inputs the spec implies and that are most likely to go wrong; each has a test in the task
named.

1. An excluded rapid-equilibrium isomerization or exchange step whose far side is cut off:
   the rule rejects it for now, and the message must say why and offer the steady-state
   fix. (Task 2)
2. A kinetic group that mixes a step losing its affinity with a rapid-equilibrium step
   that cannot be switched off: the fault names the second step, not the first. (Task 2)
3. A one-metabolite uptake from one cut-off form into another: it loses its affinity, and
   free E must still reach both along such steps. (Task 2)
4. A rapid-equilibrium step excluded from the active conformation by `:OnlyI`: the fault
   and its message name the active conformation and `:OnlyI`. (Task 3)
5. The verdict does not depend on `catalytic_multiplicity` or on regulatory sites: the
   same catalytic graph and states give the same fault at 1 and 2 subunits with a site.
   (Task 2)

---

### Task 1: Record today's verdicts on the gate populations

**Files:**
- Create (outside the repository): `/home/denis.linux/claude-scratchpad/lean-mechanism-spike/gate/populations.jl`
- Create (outside the repository): `/home/denis.linux/claude-scratchpad/lean-mechanism-spike/gate/baseline.jls`

**Interfaces:**
- Produces: `gate_populations() -> Vector{NamedTuple{(:name, :catalytic, :states)}}`, a
  deterministic list (fixed seeds), each entry a plain `Mechanism` and one allosteric
  state per kinetic group; and `baseline.jls`, a serialized
  `Dict{String, Tuple{Bool, Symbol}}` from `"$name|$(sprint(show, catalytic))|$states"` to
  `(constructor accepts, Symbol(_guard(am).fault))`.

- [ ] **Step 1: Create the branch**

Run: `git switch -c guard-rule-excluded-steps`

- [ ] **Step 2: Write `gate_populations()`**

Four populations, seeds fixed with `MersenneTwister`:
`:onlyA`, every nonempty set of groups `:OnlyA`, the rest `:EqualAI`, for every plain
`init_mechanisms` seed of uni-uni (`S[C]` to `P[C]`), bi-bi (`A[C], B[N]` to `P[C], Q[N]`)
and ping-pong (`A[CX], B[N]` to `P[C], Q[NX]`): 9,616 entries;
`:onlyI`, for 300 random plain mechanisms from init plus one `expand_mechanisms` round of
each of bi-bi, ping-pong and bi-bi with `competitive_inhibitors: A, I`, 100 random state
vectors (each group `:OnlyI` with probability 0.25, else `:OnlyA` with 0.2, else
`:EqualAI`) that hold an `:OnlyI`, duplicates dropped;
`:terter`, 2,000 random plain init seeds of a ter-ter reaction, one random nonempty
`:OnlyA` set each. `spike_rule18.jl` and `spike_rule20.jl` show the loops.

- [ ] **Step 3: Write the baseline**

For each entry build `AllostericMechanism(catalytic, states, 1, RegulatorySite[])`, catching
the constructor's `:OnlyI` error as `(false, :none)`. Serialize the Dict.

- [ ] **Step 4: Verify**

Run the script. Expected: it prints the number of entries per population (`:onlyA` 9,616)
and, for `:onlyA`, 6,748 entries with `(true, :_NO_FAULT)`.

Nothing to commit: both files live outside the repository.

### Task 2: The rule in the guard

**Files:**
- Modify: `src/Derivation/conformations.jl` (new `_SwitchOff`, `_switch_offs`)
- Modify: `src/Derivation/guard.jl` (`_Fault` enum and its docstring, `_ADMISSION_FAULTS`,
  `_eps_network`, the `_guard` docstring)
- Modify: `src/Derivation/admission.jl` (`_arcs`, `_island_edge`, `_admission_interior`,
  `_admission_message`)
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (the five oracle
  helpers)
- Test: `test/Derivation/test_admission.jl`

**Interfaces:**
- Produces, all internal to `Derivation`:
  - `@enum _SwitchOff::UInt8 _RUNS _LOSES_AFFINITY _STOPPED _CANNOT_SWITCH_OFF`
  - `_switch_offs(am::Mechanisms.AllostericMechanism, c::Int) -> Vector{_SwitchOff}`, one
    per catalytic step for conformation `c` (1 active, 2 inactive), by the table in Global
    Constraints; `_RUNS` for a step `c` does not exclude.
  - fault `_FAST_STEP_PINNED`, which replaces `_NO_LIMIT` at the same place in the enum
    and in `_ADMISSION_FAULTS`.
- Produces, in the test helpers: `_testhelper_switch_off(am, c::Int, k::Int) -> Symbol`
  (`:runs`, `:loses_affinity`, `:stopped`, `:cannot_switch_off`), written from the table
  with the helper file's own union-find, not by calling `Derivation`.

- [ ] **Step 1: Give the oracles their own statement of the rule**

Add `_testhelper_switch_off`. In `_testhelper_exclusion`, `_testhelper_network`,
`_testhelper_epsilon_span`, `_testhelper_corner_directions` and
`_testhelper_lp_existence`, replace each `EnzymeRates.Mechanisms._is_plain_binding(...)`
read by `_testhelper_switch_off(...) === :loses_affinity`. `_testhelper_exclusion` reports
`failure = :fast` when any excluded step is `:cannot_switch_off` (this replaces
`:inside`), so `_testhelper_epath_verdict` returns `:no_path` for it.
`_testhelper_admission_agrees` already maps a fault to `:no_path`.

- [ ] **Step 2: Write the failing tests**

A testset `"how an excluded step is switched off"` in `test_admission.jl`, each mechanism
written inline with `@_testhelper_unguarded @allosteric_mechanism`, unlisted steps
`:EqualAI`, one catalytic subunit:

```julia
D = EnzymeRates.Derivation
# 1. RE binding into a cut-off form: admitted.
#    E + S ⇌ E(S) :: OnlyA;  E(S) <--> E(P) :: OnlyA;  E + P ⇌ E(P)
off = D._switch_offs(k_type, 2)
@test off[_testhelper_edge(k_type, "E + S ⇌ ES")] === D._LOSES_AFFINITY
@test off[_testhelper_edge(k_type, "ES <--> EP")] === D._STOPPED
@test off[_testhelper_edge(k_type, "E + P ⇌ EP")] === D._RUNS
@test all(==(D._RUNS), D._switch_offs(k_type, 1))
@test D._guard(k_type).fault === D._NO_FAULT
@test all(_testhelper_admission_agreement(k_type; rng))
# 2. RE binding whose bound form stays: fault.
#    E + S ⇌ E(S) :: OnlyA;  E(S) <--> E(P);  E + P ⇌ E(P)
v = D._guard(pinned)
@test (v.fault, v.conformation, v.edge) ==
      (D._FAST_STEP_PINNED, 2, _testhelper_edge(pinned, "E + S ⇌ ES"))
@test _testhelper_epath_verdict(pinned; rng) === :no_path
# 3. The same with a steady-state binding: admitted, stopped.
#    E + S <--> E(S) :: OnlyA;  E(S) <--> E(P);  E + P ⇌ E(P)
@test D._STOPPED in D._switch_offs(shut_door, 2)
@test D._guard(shut_door).fault === D._NO_FAULT
@test all(_testhelper_admission_agreement(shut_door; rng))
# 4. RE step that binds and reacts at once, its bound form still reached: fault.
#    E + A <--> E(A);  E + Q ⇌ E(Q);  E(A) + B ⇌ E(P, Q) :: OnlyA;  E(Q) + P <--> E(P, Q)
@test D._guard(fused_re).fault === D._FAST_STEP_PINNED
# 5. The same with that step steady-state (`E(A) + B <--> E(P, Q) :: OnlyA`): admitted.
@test D._guard(fused_ss).fault === D._NO_FAULT
@test all(_testhelper_admission_agreement(fused_ss; rng))
```

Add one testset per Review Focus line 1, 2, 3 and 5:
(1) `E + S ⇌ E(S)`, `E(S) ⇌ E(P) :: OnlyA`, `E + P <--> E(P) :: OnlyA` gives
`_FAST_STEP_PINNED` on the isomerization, and its message ends with the steady-state fix
when `D._steady_state_passes` holds;
(2) the group `(E + A ⇌ E(A), E(B) + A ⇌ E(P, Q)) :: OnlyA` in a random bi-bi whose
`E(P, Q)` stays: the verdict's `edge` is the second step;
(3) a chain `E + A ⇌ E(A) :: OnlyA`, `E(A) + B ⇌ E(A, B) :: OnlyA` with the chemistry and
product steps also `:OnlyA`: both bindings are `_LOSES_AFFINITY`, fault `_NO_FAULT`,
agreement on every face;
(5) fixture 2 and fixture 3 rebuilt with `catalytic_multiplicity: 2` and
`regulatory_site(2): R::OnlyI` give the same faults.

- [ ] **Step 3: Run the tests to verify they fail**

Run the focused command on `test/Derivation/test_admission.jl`.
Expected: the new testsets fail with `UndefVarError: _switch_offs` (and
`_FAST_STEP_PINNED`).

- [ ] **Step 4: Implement `_switch_offs` and wire it in**

`_switch_offs` computes `exclusion_components` once per call. Then:
- `_arcs(am, c)`: the steps with `_LOSES_AFFINITY`.
- `_eps_network`: `plain` is `_switch_offs(am, c)[k] === _LOSES_AFFINITY`; compute the
  classification once per conformation, not per step.
- `_admission_interior`: where it returns `_NO_LIMIT` today, return
  `_Verdict(_FAST_STEP_PINNED, _LIMIT_UNTESTED, 0, c, k)` for the first step with
  `_CANNOT_SWITCH_OFF`. The arc-inside-a-component test goes: `_LOSES_AFFINITY` steps
  never lie inside one component.
- `_island_edge`: rank 3 is a step with `length(taken_up) == 1 && isempty(given_off)`
  instead of `step_kinds[k] === :binding`.
- `_admission_message`, the `_FAST_STEP_PINNED` branch. Cause, exact text, with `$conf`,
  `$tag` and `$(line(k))` as the other branches use them:
  `"AllostericMechanism: in the $conf conformation the rapid-equilibrium step $(line(k)) is tagged $tag but cannot be switched off. A step keeps its kind in both conformations, so a rapid-equilibrium step can only lose its affinity, and it can do that only when it takes up one metabolite into a form the other steps no longer reach."`
  Then the fixes, each only when its variant passes the guard, as the branch does today:
  the steady-state fix (`_steady_state_passes`), and for `c == 2` the existing
  `fixes([k], "its group", _kept_path(...), "its sides")`.
- In the `_ISLAND` branch, the ping-pong hint reads `step_kinds` and form contents. Keep it
  only if a fixture in `test_admission.jl` still reaches it; otherwise delete the hint and
  its test assertion, and say so in the commit message.
- Update the docstrings of `_Fault`, `_eps_network`, `_arcs`, `_admission_interior` and
  `_guard` to the rule; no "spec section" reference names the old rule.

- [ ] **Step 5: Bring the file's existing testsets to the rule**

Run the focused command. Nine existing testsets assert the old verdicts or messages
(found in the spike): "a trap across both conformations labels each form", "a bridge of
several tagged steps names each one", "a tagged binding that untagged steps join on both
sides has no limit", "tagged bindings that close a directed cycle have no limit", "forms
reached only through tagged chemistry are an island", "the macros reject what the guard
rejects when the expanded code runs", "a fast cycle that only the inactive conformation
keeps is named there", "the active conformation's interior faults name the :OnlyI step",
and any that read `.failure === :inside`.
For each failing assertion apply exactly one of the spec's four allowed changes. A fixture
that becomes admitted gets `@test all(_testhelper_admission_agreement(x; rng))` in place
of its fault assertion. Rename a testset whose title states the old rule. Keep a list of
every changed assertion (file line, old expectation, new expectation, which of the four)
in `/home/denis.linux/claude-scratchpad/lean-mechanism-spike/gate/changed-tests.md` for
the pull request.

- [ ] **Step 6: Run the tests to verify they pass**

Run the focused command on `test/Derivation/test_admission.jl`. Expected: all pass, no
warnings.

- [ ] **Step 7: Commit**

```bash
git add src/Derivation/conformations.jl src/Derivation/guard.jl src/Derivation/admission.jl \
        test/mechanism_definitions_for_test_enzyme_derivation.jl test/Derivation/test_admission.jl
git commit -m "Read how an excluded step is switched off from the graph"
```

### Task 3: `:OnlyI` on any group

**Files:**
- Modify: `src/Mechanisms/types.jl` (delete `_check_onlyi_groups` and its call at line 516)
- Modify: `src/DSL/macros.jl` (the `@allosteric_mechanism` docstring, the two paragraphs
  on tags at lines 199 to 213)
- Test: `test/Derivation/test_admission.jl`, `test/Mechanisms/test_types.jl`,
  `test/Mechanisms/test_reaction_graph.jl`, `test/DSL/test_dsl.jl`

**Interfaces:**
- Consumes: `_FAST_STEP_PINNED`, `_switch_offs`, `_testhelper_admission_agreement` from
  Task 2.
- Produces: `AllostericMechanism` constructors accept `:OnlyI` on every group.

- [ ] **Step 1: Write the failing tests**

In `test_admission.jl`, a testset `":OnlyI on a step that does more than bind a modifier"`
with these mechanisms inline through `@allosteric_mechanism` (guarded), one subunit:

```julia
# Two isomerization paths; the second belongs to the inactive conformation alone.
#   E + S ⇌ E(S) :: EqualAI;  E(S) <--> E(P) :: OnlyI;  E(S) ⇌ Estar(S) :: EqualAI
#   Estar(S) <--> Estar(P) :: EqualAI;  Estar(P) ⇌ E(P) :: EqualAI;  E + P ⇌ E(P) :: EqualAI
@test all(_testhelper_admission_agreement(two_paths; rng))
@test all(_testhelper_admission_agreement(_testhelper_mirror(two_paths); rng))
# An exchange step beside a ternary-complex path, the exchange :OnlyI.
#   E + A ⇌ E(A);  E(A) + B <--> E(Q) + P :: OnlyI;  E(A) + B ⇌ E(A, B)
#   E(A, B) <--> E(P, Q);  E(Q) + P ⇌ E(P, Q);  E + Q ⇌ E(Q)
@test all(_testhelper_admission_agreement(exchange; rng))
# The only catalytic path :OnlyI: the active conformation cannot catalyze.
#   E + S ⇌ E(S) :: EqualAI;  E(S) <--> E(P) :: OnlyI;  E + P ⇌ E(P) :: EqualAI
@test_throws "the active conformation turns no substrate into product" @allosteric_mechanism ...
# Review Focus 4: a rapid-equilibrium reactant binding :OnlyI whose bound form stays.
#   E + S ⇌ E(S) :: OnlyI;  E(S) <--> E(P) :: EqualAI;  E + P ⇌ E(P) :: EqualAI
v = D._guard(pinned_active)   # built unguarded
@test (v.fault, v.conformation) == (D._FAST_STEP_PINNED, 1)
@test occursin("in the active conformation", D._guard_message(pinned_active, v))
@test occursin(":OnlyI", D._guard_message(pinned_active, v))
```

`_testhelper_mirror(am)` goes in the shared helper file: the same catalytic graph,
multiplicity and sites with every `:OnlyA` and `:OnlyI` group state swapped.

- [ ] **Step 2: Run to verify they fail**

Expected: the first three error with "is tagged :OnlyI but does more than bind a modifier".

- [ ] **Step 3: Delete `_check_onlyi_groups` and its call; rewrite the docstring paragraphs**

The docstring says: a catalytic statement takes any of the four tags; the active
conformation must still catalyze; how an excluded step is switched off, in three
sentences from the spec's table; and a link to the docs section of Task 7,
`[How an excluded step is switched off](@ref)`.

- [ ] **Step 4: Bring the constructor's old tests to the rule**

Run the focused command on each of `test/Mechanisms/test_types.jl`,
`test/Mechanisms/test_reaction_graph.jl` and `test/DSL/test_dsl.jl`. Each assertion of the
constructor's `:OnlyI` error becomes an assertion of the guard's verdict on that mechanism
(the spec's fourth allowed change): the guard's fault and message if it rejects, oracle
agreement if it admits. Append each to `changed-tests.md`.

- [ ] **Step 5: Run to verify they pass**

Focused runs of the four test files. Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add src/Mechanisms/types.jl src/DSL/macros.jl test/Derivation/test_admission.jl \
        test/Mechanisms/test_types.jl test/Mechanisms/test_reaction_graph.jl \
        test/DSL/test_dsl.jl test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "Let any group be :OnlyI while the active conformation catalyzes"
```

### Task 8: The second try for a steady-state uptake (runs after Task 3, before Task 4)

The gate (Task 4) found two ter-ter mechanisms, admitted before this work, that the guard
now rejects with `_BRIDGE` though their law is the limit. In each, a steady-state step
that binds and reacts at once leads into a cut-off form; the old guard read it as stopped,
the rule reads it as losing its affinity. Both readings keep the step steady-state, so the
maintainer decided the guard tries both.

**Files:**
- Modify: `src/Derivation/conformations.jl` (`_switch_offs`), `src/Derivation/guard.jl`
  (`_guard`, `_eps_network`), `src/Derivation/admission.jl` (`_arcs`,
  `_admission_interior`, `_admission`), and any docstring in `src/` that states the rule
  for a steady-state step (`_SwitchOff`, `_Fault`, `@allosteric_mechanism`)
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/Derivation/test_admission.jl`

**Interfaces:**
- Consumes: `_switch_offs(am, c)`, `_SwitchOff`, `_testhelper_switch_off(am, c, k)`,
  `_testhelper_admission_agreement(am; rng)`, `_testhelper_mirror(am)` from Tasks 2 and 3.
- Produces:
  - `_switch_offs(am, c::Int, slow_uptake::_SwitchOff = _LOSES_AFFINITY)`: `slow_uptake`
    is how a steady-state step that takes up one metabolite into a cut-off form is switched
    off, `_LOSES_AFFINITY` or `_STOPPED`. Nothing else in the result depends on it.
  - The same trailing argument, with the same default, on `_arcs`, `_admission_interior`,
    `_admission` and `_eps_network`, so every existing call keeps the first reading.
  - `_guard(am)`: the first try as now. If its fault is one of `_ADMISSION_FAULTS`, the
    interior and face checks run again with `_STOPPED`; `_WELL_POSED` if that try finds no
    fault, else the first try's verdict unchanged. Classes 1 to 3 do not depend on the
    reading. A second try that would read no step differently is not run.
  - `_testhelper_switch_off(am, c, k; slow_uptake = :loses_affinity)` (`:stopped` for the
    second reading), threaded through the oracle helpers that read it, and
    `_testhelper_admission_agreement(am; rng)` judging under the reading the guard admitted
    with: `:stopped` when the guard admits and the first reading alone does not, else the
    first. The helpers keep their own rule; they do not call `_switch_offs`.

The fixtures, each written out in full in its testset with `@allosteric_mechanism` (they
are entries of the gate's recorded populations; their text is in
`/home/denis.linux/claude-scratchpad/lean-mechanism-spike/gate/report.md`, "Violations").
Build them unguarded where the guard would refuse them.

- `second_try_a` (entry 54443) and `second_try_b` (entry 54593): admitted by the second
  try alone.
- `bridge_fast_a` (entry 53145) and `bridge_fast_b` (entry 53246): stay rejected.

- [ ] **Step 1: Write the failing tests**

A testset `"a steady-state uptake into a cut-off form may be stopped"` in
`test/Derivation/test_admission.jl`:

```julia
# second_try_a, inactive conformation: the chemistry step joins two cut-off forms.
k = _testhelper_edge(second_try_a, "E(B, D) + A <--> E(P, Q, R)")
first, second = D._switch_offs(second_try_a, 2), D._switch_offs(second_try_a, 2, D._STOPPED)
@test first[k] === D._LOSES_AFFINITY
@test second[k] === D._STOPPED
@test [first[j] for j in eachindex(first) if j != k] ==
      [second[j] for j in eachindex(second) if j != k]
@test _testhelper_switch_off(second_try_a, 2, k) === :loses_affinity
@test _testhelper_switch_off(second_try_a, 2, k; slow_uptake = :stopped) === :stopped

for am in (second_try_a, second_try_b)
    @test D._admission(am).fault === D._BRIDGE            # the first try alone
    @test D._admission(am, D._STOPPED) == D._WELL_POSED
    @test D._guard(am) == D._WELL_POSED
    @test all(_testhelper_admission_agreement(am; rng = MersenneTwister(1)))
end
```

and, in the same testset:
- the first try's `_BRIDGE` on `second_try_a` is at the face with A, B and D zero, in
  conformation 2, on the edge `E(D, P) + Q ⇌ E(D, P, Q)`;
- the law of `second_try_a` and of `second_try_b` equals the oracle's own deletion limit
  (`_testhelper_check_law_against_limit`), so the admitted law is compared with a
  separately derived one;
- `bridge_fast_a`: `D._guard(am) == D._admission(am)`, fault `_BRIDGE`, conformation 2,
  edge `E(D, Q) + P ⇌ E(D, P, Q)`; `D._admission(am, D._STOPPED).fault === D._ISLAND`; the
  ε-path oracle under the first reading is `:direction_dependent` at the face with P zero;
- `bridge_fast_b`: fault `_BRIDGE`, and `D._switch_offs(am, c) ==
  D._switch_offs(am, c, D._STOPPED)` for both `c` (no steady-state uptake to read again);
- the first try decides when it passes: for the K-type uni-uni mechanism
  `E + S <--> E(S) :: OnlyA`, `E(S) <--> E(P) :: OnlyA`, `E(P) ⇌ E + P`, written in the
  testset, `D._guard(am) == D._WELL_POSED` and
  `D._admission_interior(am, D._STOPPED).fault === D._ISLAND`;
- the two conformations are read alike: `D._switch_offs(_testhelper_mirror(am), 1, r) ==
  D._switch_offs(am, 2, r)` for `am = second_try_a` and both readings `r`.

- [ ] **Step 2: Run them and watch them fail**

Focused run of `test/Derivation/test_admission.jl`. Expected: a `MethodError` on the
three-argument `_switch_offs`.

- [ ] **Step 3: Implement**

As the Interfaces block says. Update each docstring that states how a steady-state step is
switched off; the guard's error messages stay the first try's.

- [ ] **Step 4: Run the tests and the guards**

- Focused runs: `test/Derivation/test_admission.jl`, `test/Mechanisms/test_types.jl`,
  `test/Mechanisms/test_reaction_graph.jl`, `test/DSL/test_dsl.jl`,
  `test/Enumeration/test_mechanism_enumeration.jl`. No existing assertion may change. **If
  one would have to, a population pin included, stop and report it.**
- `test/test_type_stability.jl` and `test/test_boundaries.jl` pass (the default `julia`
  is 1.13.4; `julia +1.13` is not installed).
- Report the wall time of the enumeration file's focused run before and after the change:
  enumeration calls the guard for every tagging it proposes.

- [ ] **Step 5: Commit**

```bash
git add src/Derivation test/Derivation/test_admission.jl test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "Let a steady-state uptake into a cut-off form be stopped on a second try"
```

### Task 9: The third try, for a mechanism that needs both readings (runs before Task 5)

Task 5 found 64 mechanisms, admitted before this work and outside the gate's populations
(K-type taggings of the bi-bi split children), that the two tries of Task 8 reject: one
excluded steady-state binding must lose its affinity and another must stop. The maintainer
decided the guard also tries the reading in which such a step loses its affinity only
when it leads away from free E (spec, decision 2). No single reading replaces the three.

**Files:**
- Modify: `src/Derivation/conformations.jl` (`_switch_offs`), `src/Derivation/guard.jl`
  (`_guard`, `_eps_network`), `src/Derivation/admission.jl` (`_arcs`,
  `_admission_interior`, `_admission`), and each docstring in `src/` that states the two
  tries
- Modify (outside the repository):
  `/home/denis.linux/claude-scratchpad/lean-mechanism-spike/gate/compare.jl`, `report.md`
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/Derivation/test_admission.jl`

**Interfaces:**
- Consumes: Task 8's `slow_uptake::_SwitchOff` argument and its test helpers, which this
  task replaces.
- Produces:
  - `@enum _SlowUptakeReading::UInt8 _ALL_LOSE_AFFINITY _OUTWARD_LOSE_AFFINITY _ALL_STOPPED`
    in `conformations.jl`: which excluded steady-state steps that take up one metabolite
    into a cut-off form lose their affinity, the rest being stopped. Its order is the
    order of the guard's tries.
  - `_switch_offs(am, c::Int, reading::_SlowUptakeReading = _ALL_LOSE_AFFINITY)`, and the
    same trailing argument, with the same default, in place of `slow_uptake` on `_arcs`,
    `_admission_interior`, `_admission` and `_eps_network`. No `slow_uptake` argument
    remains. A rapid-equilibrium step is read as before under every reading.
  - `_guard(am)`: the first try as now. When its fault is one a reading can change, the
    interior and face checks run under `_OUTWARD_LOSE_AFFINITY`, then under `_ALL_STOPPED`;
    `_WELL_POSED` at the first try that finds no fault, else the first try's verdict
    unchanged. A try that reads no step differently from an earlier try is not run. Task 8
    retries on `_ARC_CYCLE`, `_ISLAND`, `_PUMPED` and `_BRIDGE` only; check against the
    code that the faults it skips cannot change under the outward reading either, and say
    in the report how you checked.
  - `_testhelper_switch_off(am, c, k; reading = :all_lose_affinity)`
    (`:outward_lose_affinity`, `:all_stopped`), in place of the `slow_uptake` keyword, with
    its own search for the distances; it does not call `_switch_offs`. An unknown symbol
    is an error.
  - `_testhelper_admission_agreement(am; rng)` judges under the reading the guard admitted
    with: the first of the three, in the guard's order, under which `D._admission(am, r)`
    finds no fault; the first reading when none does.

The fixtures, each written out in full in its testset with `@allosteric_mechanism`, built
unguarded where the guard would refuse them. Their text is in
`/home/denis.linux/claude-scratchpad/lean-mechanism-spike/gate/fixtures3.out` under the
entry number given; `second_try_a`, `second_try_b`, `bridge_fast_a` and `bridge_fast_b`
are Task 8's.

- `both_readings` (split entry 5644; admitted before this work). The inactive conformation
  keeps `E(A) + P <--> E(A, P)` and `E(P) + Q ⇌ E(P, Q)`, so its cut-off parts are
  `{E(A), E(A, P)}` and `{E(P), E(P, Q)}`, each at distance 1.
- `both_readings_only_i` (`:onlyI` entry 13305; not admitted before).
- `outward_pumped` (ter-ter entry 53328; admitted before).

- [ ] **Step 1: Write the failing tests**

A testset `"a mechanism may need one uptake to lose its affinity and another to stop"` in
`test/Derivation/test_admission.jl`:

```julia
am = both_readings
p = _testhelper_edge(am, "E + P <--> E(P)")              # leads away: distance 0 to 1
b = _testhelper_edge(am, "E(A) + B <--> E(P, Q)")        # distance 1 to 1
a = _testhelper_edge(am, "E(P) + A ⇌ E(A, P)")           # rapid equilibrium, 1 to 1
for (reading, name, at_p, at_b) in (
        (D._ALL_LOSE_AFFINITY, :all_lose_affinity, D._LOSES_AFFINITY, D._LOSES_AFFINITY),
        (D._OUTWARD_LOSE_AFFINITY, :outward_lose_affinity, D._LOSES_AFFINITY, D._STOPPED),
        (D._ALL_STOPPED, :all_stopped, D._STOPPED, D._STOPPED))
    off = D._switch_offs(am, 2, reading)
    @test off[p] === at_p && off[b] === at_b
    @test off[a] === D._LOSES_AFFINITY                    # a fast step is read one way
    @test all(k -> Symbol(lowercase(string(off[k])[2:end])) ===
                   _testhelper_switch_off(am, 2, k; reading = name), eachindex(off))
end
@test D._admission(am).fault === D._ARC_CYCLE
@test D._admission(am, D._ALL_STOPPED).fault === D._ISLAND
@test D._admission(am, D._OUTWARD_LOSE_AFFINITY) == D._WELL_POSED
@test D._guard(am) == D._WELL_POSED
@test all(_testhelper_admission_agreement(am; rng = MersenneTwister(1)))
```

and, in the same testset:
- the law of `both_readings` equals the oracle's own deletion limit
  (`_testhelper_check_law_against_limit`);
- `both_readings_only_i`: the same three verdicts per reading as `both_readings`,
  `D._guard(am) == D._WELL_POSED`, the agreement, and the law against the limit;
- `outward_pumped`: `D._admission(am) == D._WELL_POSED`,
  `D._admission(am, D._OUTWARD_LOSE_AFFINITY).fault === D._PUMPED`, `D._guard(am) ==
  D._WELL_POSED`, and the agreement: the first try decides, and the outward reading alone
  would lose this mechanism;
- `second_try_a`: `D._admission(am, D._OUTWARD_LOSE_AFFINITY) == D._WELL_POSED`;
  `second_try_b`: `D._admission(am, D._OUTWARD_LOSE_AFFINITY).fault === D._BRIDGE` and
  `D._admission(am, D._ALL_STOPPED) == D._WELL_POSED`, so each later try is needed;
- `bridge_fast_a`: rejected by all three, and `D._guard(am) == D._admission(am)` (the
  first try's verdict);
- the two conformations are read alike: `D._switch_offs(_testhelper_mirror(am), 1, r) ==
  D._switch_offs(am, 2, r)` for `am = both_readings` and every reading `r`;
- a step out of a part no path reaches is stopped under the outward reading: assert it on
  a mechanism written in the testset, with the step named, if `both_readings` and Task 8's
  fixtures hold none.

Task 8's testset keeps every assertion, with `D._STOPPED` as a reading written
`D._ALL_STOPPED` and `slow_uptake = :stopped` written `reading = :all_stopped`. No
assertion is weakened or removed.

- [ ] **Step 2: Run them and watch them fail**

Focused run of `test/Derivation/test_admission.jl`. Expected: `UndefVarError` on
`_ALL_LOSE_AFFINITY`.

- [ ] **Step 3: Implement**

As the Interfaces block says. Update each docstring and comment that states two tries;
the guard's error messages stay the first try's. `docs/src/` waits for Task 7.

- [ ] **Step 4: Run the tests and the guards**

- Focused runs: `test/Derivation/test_admission.jl`, `test/Mechanisms/test_types.jl`,
  `test/Mechanisms/test_reaction_graph.jl`, `test/DSL/test_dsl.jl`. No existing assertion
  may change but Task 8's renames. **If one would have to, a population pin included, stop
  and report it.**
- `test/Enumeration/test_mechanism_enumeration.jl` fails at this commit (32,144 pass, 41
  fail, 16 error at 14f8513a; Tasks 5 and 6 bring it to the rule). Report the failing
  testsets before and after this change and the wall time of both runs; a difference in
  which fail is expected only where the guard now admits a mechanism.
- `test/test_type_stability.jl` and `test/test_boundaries.jl` pass (the default `julia`
  is 1.13.4; `julia +1.13` is not installed).

- [ ] **Step 5: Commit**

```bash
git add src/Derivation test/Derivation/test_admission.jl test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "Let the guard read an outward steady-state uptake apart from the others"
```

- [ ] **Step 6: Rerun the gate with the split children added**

In `compare.jl`, add a fourth population, `:split`: the 9,878 mechanisms whose text is in
`gate/split_population.jls` (field `allparents`; rebuild each as `level_rule.jl`'s
`unshow` does), with the guard's faults at the branch base in `gate/split_base_faults.jls`
(a vector of symbols in the same order; admitted before means `:_NO_FAULT`). The checks of
Task 4 step 1 apply to it as to `:onlyA`. Keep the last run's files as
`*-before-third-try.*`. Rerun, and write `report.md`.

Expected, from the spike (`gate/level_rule.out`, `fixtures3.out`): admitted now 7,520
(`:onlyA`), 2,200 (`:onlyI`), 1,350 (`:terter`), 7,474 (`:split`); the only lost entries
without the new fault are 53145 and 53246; the outward try admits 5 `:onlyI`, 2 `:terter`
and 64 `:split` entries the other two do not, and every one must give
`all(_testhelper_admission_agreement(am; rng))`. **If a count differs, an entry admitted
before is lost for another reason, or an oracle disagrees, stop and report it.**

### Task 10: A reading names how every excluded step is switched off (runs after Task 9)

A refactor with no change of verdict: the guard's checks take a whole reading, in which a
step may lose in either direction, so that Tasks 11 and 12 can feed them any reading.

**Files:**
- Modify: `src/Derivation/conformations.jl`, `src/Derivation/admission.jl`,
  `src/Derivation/guard.jl`, `src/Derivation/messages.jl`
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/Derivation/test_admission.jl`, and any test that names `_LOSES_AFFINITY` or
  `_switch_offs`

**Interfaces:**
- Produces, internal to `Derivation`:
  - `@enum _SwitchOff::UInt8 _RUNS _LOSES_FORWARD _LOSES_REVERSE _STOPPED _CANNOT_SWITCH_OFF`.
    `_LOSES_FORWARD`: the step's stored direction is uphill and vanishes, its reverse keeps
    its rate (what `_LOSES_AFFINITY` meant; that name goes). `_LOSES_REVERSE`: the reverse.
  - `const _Reading = Matrix{_SwitchOff}`: `reading[k, c]` for catalytic step `k` in
    conformation `c` (1 active, 2 inactive).
  - `_built_reading(am, r::_SlowUptakeReading = _ALL_LOSE_AFFINITY) -> _Reading`, in place
    of `_switch_offs(am, c, r)`.
  - `_admission_interior(am, reading::_Reading)`, `_admission(am, reading::_Reading)` and
    `_eps_network(am, re::BitVector, reading::_Reading)`: no `_SlowUptakeReading` argument
    remains on them, and no default reading. The interior follows each losing step from
    its lower side's group to its higher side's, and also refuses a losing step that points
    into free E's group (`_ARC_CYCLE`, since no energies exist). `_eps_network` marks the
    forward transition of a `_LOSES_FORWARD` step and the reverse of a `_LOSES_REVERSE`
    step as vanishing; either stays fast when the step is rapid-equilibrium.
  - `_guard(am)` tries `_built_reading(am, r)` for the three `r` exactly as now.
- Produces, in the test helpers (their own statement of the rule; they do not call the
  package's):
  - `_testhelper_reading(am, name::Symbol) -> Matrix{Symbol}` for the three built names,
    with entries `:runs`, `:loses_forward`, `:loses_reverse`, `:stopped`,
    `:cannot_switch_off`; and `_testhelper_reading(reading::D._Reading)`, which only
    renames the entries.
  - The keyword `reading` of `_testhelper_exclusion`, `_testhelper_network`,
    `_testhelper_epath_verdict`, `_testhelper_corner_directions`,
    `_testhelper_standard_directions`, `_testhelper_epsilon_span` and
    `_testhelper_admission_faults` takes a `Matrix{Symbol}` (a built name still works).
    Arcs point from the lower group to the higher. A steady-state step that loses keeps
    the rate out of its higher side whichever way it is stored: with `d` the potential's
    rise along the stored direction, `(f ε^d, r)` when `d >= 0` and `(f, r ε^-d)` when
    `d < 0`. `_testhelper_switch_off(am, c, k; reading)` goes or becomes an index into the
    matrix.

- [ ] **Step 1: Write the failing tests**

A testset `"a reading may lose a step in either direction"` in `test_admission.jl`, on the
six-step ping-pong of the spec's fourth example, written in the testset with both
chemistry steps `:OnlyA` and built unguarded. With `r` a `_Reading` that is `_RUNS`
everywhere except `_LOSES_FORWARD` at `EA <--> FP` and `_LOSES_REVERSE` at `FB ⇌ EQ` in
conformation 2:

```julia
@test D._admission_interior(am, r) == D._WELL_POSED
eps = D._eps_network(am, D._re_flags(am.catalytic), r)
# FB ⇌ EQ in copy 2: the reverse (EQ → FB) vanishes, the forward does not, both fast.
@test [eps.vanishing[t] for t in (fwd, rev)] == [false, true]
@test eps.fast[fwd] && eps.fast[rev]
# Both steps pointing into the covalent group's lower neighbour: no energies exist.
@test D._admission_interior(am, flipped).fault === D._ARC_CYCLE
# The tests' own network agrees on which rates vanish.
@test _testhelper_reading(r) == expected_symbols
```

with `fwd`, `rev` the two transitions of that step in copy 2 found through `eps.edge` and
`eps.conformation`, and `flipped` the reading with the two directions exchanged. Add: the
ε-path oracle under `_testhelper_reading(r)` converges on every face with
`corners = true`, and under `_testhelper_reading(am, :all_lose_affinity)` it is `:no_path`.

- [ ] **Step 2: Run them and watch them fail** (`UndefVarError: _LOSES_REVERSE`).

- [ ] **Step 3: Implement**, migrating every caller and removing `_switch_offs` and
`_LOSES_AFFINITY` in the same change.

- [ ] **Step 4: Run the tests and the guards**

Focused runs of `test/Derivation/test_admission.jl`, `test/Mechanisms/test_types.jl`,
`test/Mechanisms/test_reaction_graph.jl`, `test/DSL/test_dsl.jl`; `test/test_boundaries.jl`
and `test/test_type_stability.jl`. No assertion changes but the renames. The enumeration
file fails at this commit exactly as before it (32,144 pass, 41 fail, 16 error): report
its failing testsets and wall time before and after. Rerun
`gate/compare.jl` and confirm every verdict is unchanged. **If a verdict changes, stop.**

- [ ] **Step 5: Commit**

`git commit -m "Let a reading switch a step off in either direction"`

### Task 11: The two additions to the face checks

**Files:**
- Modify: `src/Derivation/admission.jl`, `src/Derivation/guard.jl`,
  `src/Derivation/messages.jl`
- Test: `test/Derivation/test_admission.jl`

**Interfaces:**
- Consumes: Task 10's `_Reading`, `_eps_network(am, re, reading)`, `_admission(am, reading)`.
- Produces, internal to `Derivation`:
  - A fault `_FILLED` in `_Fault` and in `_ADMISSION_FAULTS`: on a face, a rapid-equilibrium
    transition that vanishes is alive while its reverse, which does not vanish, is dead,
    its source is reached by the enzyme and from the law's forms along alive
    rapid-equilibrium transitions, and its target reaches no law form along alive
    rapid-equilibrium transitions that do not vanish (the spec's "No form fills faster
    than it empties"). `_face_condition` reports it before `_PUMPED`, with the
    conformation and the step. Its message names the face, the step, the form it fills
    and the form's slow way out, and offers the steady-state fix where that variant passes.
  - The set-aside rule exactly as the spec's "Groups that cannot matter are set aside"
    states it: when `_face_condition` refuses a face with `_PUMPED` or `_BRIDGE`, it runs
    once more on the network without the groups that qualify, and the face passes if that
    run does. `_FILLED` and `_EPS_NONUNIQUE` are not retried.
- Reference only, not to copy: `pump_ok` and `without_appendages` in
  `/home/denis.linux/claude-scratchpad/lean-mechanism-spike/gate/final.jl` are the
  versions the spec's numbers come from (use its `:h1` branch; the `:h2` branch is not
  built). They patch nothing but are written as a spike.

- [ ] **Step 1: Write the failing tests**

A testset `"no form fills faster than it empties"`, each mechanism written out and built
unguarded, each reading written out as a `_Reading`:

- `filled` (gate entry 8524 of the `:onlyA` population; text and reading in
  `gate/pump_diag3.out`): `D._admission(am, r).fault === D._FILLED` on the face with B
  zero, in conformation 2, at `F + B ⇌ FB`; the ε-path oracle under
  `_testhelper_reading(r)` with `corners = true` is not `:converges` on that face, and is
  `:converges` on the interior.
- The six-step ping-pong of Task 10 under its reading: `D._admission(am, r) ==
  D._WELL_POSED` (a one-way transition never arises: neither excluded step takes up a
  metabolite).
- `slow_source` (ping-pong case 5 of `gate/pump_diag2.out`, its first reading): a fast
  one-way filling whose source is fed through a steady-state step passes, and the oracle
  with `corners = true` converges on every face.
- No built reading can raise the fault: for the three `_built_reading`s of `filled`, no
  face gives `_FILLED`.

A testset `"cut-off forms that cannot matter are set aside"`:

- `set_aside` (a mechanism the rule admits only through it: take the smallest of the 35
  in the `:onlyA` population, found by running `gate/final.jl`'s search with and without
  `SHARP[]`; write its text and reading in the testset): `_admission(am, r) ==
  _WELL_POSED`; with the group's touching forms not in one rapid-equilibrium set (the same
  mechanism with the joining step steady-state) the face stays `_BRIDGE`; the oracle with
  `corners = true` agrees with both verdicts.
- The two conservative rejections the suite already pins
  (`"conservative rejections the oracle admits stay rejected"`) stay rejected.

- [ ] **Step 2: Run them and watch them fail** (`UndefVarError: _FILLED`).

- [ ] **Step 3: Implement.**

- [ ] **Step 4: Run the tests and the guards**

As Task 10's step 4. `gate/compare.jl`: no mechanism it admits may be lost; list each
mechanism newly admitted (the set-aside rule can pass a built reading) and check every
one with `_testhelper_admission_agreement`, `corners = true` for 8 steps or fewer. **If
one is lost or one disagrees, stop.**

- [ ] **Step 5: Commit**

`git commit -m "Refuse a form that fills faster than it empties; set harmless groups aside"`

**Split in execution.** The filling check landed alone (`852d69ab`, "Refuse a form that
fills faster than it empties"): the set-aside rule as first specified admitted the
R2funnel family with laws the harsh test rejects. The set-aside rule lands as Task 11b,
with the spec's drain condition as now stated (the group's fast classes drain in the
aggregated chain) and these additions to the testset above:

- `R2funnel` and `R2funnel_SS` under the first built reading: the group at zero P is not
  set aside (its only way out that does not vanish leaves a form with a fast way out), the
  face stays `_BRIDGE`, and the oracle with `corners = true` does not converge there.
- A counter-case for the equilibrium condition in which the oracle does not converge (the
  one first planned, entry 8561 with `E + Q <--> E(Q)`, stays `_BRIDGE` though the oracle
  converges: a conservative rejection, pinned as such).
- Fixtures the rule moves to another verdict are brought to it in the spec's allowed
  ways: the trap fixture is rejected as `_TRAP`, not `_BRIDGE`; `second_try_a` and
  `second_try_b`, now admitted on the first try, are replaced by mechanisms from the
  gate's later-try entries so that the second and third tries stay exercised.

Reference only: `without_appendages` with `REPAIR[] = :agg` in `gate/repair.jl`, and its
output `gate/repair.out`. Commit: `"Set aside cut-off groups that cannot matter"`.

### Task 12: The search for a reading

**Files:**
- Modify: `src/Derivation/guard.jl`, `src/Derivation/admission.jl`,
  `src/Derivation/conformations.jl`, `src/Derivation/messages.jl`, the
  `@allosteric_mechanism` docstring in `src/DSL/macros.jl`
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/Derivation/test_admission.jl`, `test/Mechanisms/test_types.jl`,
  `test/DSL/test_dsl.jl`

**Interfaces:**
- Consumes: Tasks 10 and 11.
- Produces, internal to `Derivation`:
  - `_admitting_reading(am) -> Union{_Reading, Nothing}`: the first reading under which
    every check passes, `nothing` when none does. Order: the three built readings, then
    the other consistent readings in a fixed order that does not depend on how the steps
    are numbered beyond their stored order. `_guard(am)` admits exactly when it returns a
    reading.
  - The search as the spec's "The search for a reading": per conformation, each excluded
    step between two groups takes `_LOSES_FORWARD` (unless its `to` side is in free E's
    group), `_LOSES_REVERSE` (unless its `from` side is), or, when steady-state,
    `_STOPPED`; an excluded step inside one group is `_STOPPED` or, rapid-equilibrium,
    `_CANNOT_SWITCH_OFF`. The choices of a conformation that pass the interior are paired
    with the other's and each pair is checked by class 2 and the faces.
  - `const _MAX_CHOICES_PER_CONFORMATION = 20_000` and `const _MAX_READING_PAIRS = 5_000`:
    the search is skipped whole when either count is exceeded. Set them from the timings
    if the slowest guard call in the gate's populations is far from 0.1 s, and report it.
  - `_FAST_STEP_PINNED` now means a rapid-equilibrium step excluded inside one group; it
    is reported before any other admission fault, since no reading changes it. Its message
    says the step's two sides stay joined by the steps the conformation keeps, and offers
    the two fixes (make the step steady-state; exclude the steps that keep its sides
    joined), each only where its variant passes. Any other rejection reports the first
    built reading's verdict when that reading is consistent, else the first consistent
    reading's. No message says a step cannot be switched off when a reading switches it
    off.
- Produces, in the test helpers:
  - `_testhelper_reading_is_consistent(am, reading::Matrix{Symbol}) -> Bool`, the tests'
    own statement of the four-row table and the three conditions, from
    `_testhelper_kept_components`.
  - `_testhelper_admission_agreement(am; rng, corners = false)`: judges under
    `_testhelper_reading(D._admitting_reading(am))` when the guard admits, after asserting
    that reading consistent by the tests' own rule; under the first built reading
    otherwise.

- [ ] **Step 1: Write the failing tests**

A testset `"the guard searches for a reading"` in `test_admission.jl`, mechanisms written
out:

- the six-step ping-pong with both chemistry steps `:OnlyA`: `D._guard(am) ==
  D._WELL_POSED`; `D._admitting_reading(am)` is Task 10's reading; the law equals the
  oracle's deletion limit (`_testhelper_check_law_against_limit`);
  `all(_testhelper_admission_agreement(am; rng, corners = true))`; the same at
  `catalytic_multiplicity: 2`, and the same verdict for its mirror image;
- the three-step ping-pong `E + A <--> F + P :: OnlyA`, `E + Q <--> EQ`,
  `F + B <--> EQ :: OnlyA`: admitted, with the same four checks;
- a ping-pong with no limit under any reading (case `init 13` of `gate/levels.out`:
  `E + A <--> FP :: OnlyA`, `E + Q <--> EQ`, `F + P ⇌ FP`, `F + B ⇌ EQ :: OnlyA`):
  `D._admitting_reading(am) === nothing`, and the oracle converges under none of its
  consistent readings;
- `filled` of Task 11: rejected, since its other readings fail too;
- a mechanism admitted by a built reading returns that built reading (the search does not
  run): `both_readings` of Task 9 returns `_built_reading(am, _OUTWARD_LOSE_AFFINITY)`;
- a rapid-equilibrium step inside one group is still the fault, with its message: the
  spec's second example;
- a mechanism too large to search that the built readings reject is rejected, and one
  they admit is admitted: assert on `_MAX_CHOICES_PER_CONFORMATION` with a ter-ter
  mechanism from `gate/levels.out` ("too large") written in the testset;
- step order: the six-step ping-pong written with its steps in another order gets the
  same verdict.

- [ ] **Step 2: Run them and watch them fail** (`UndefVarError: _admitting_reading`).

- [ ] **Step 3: Implement.** Reference only: `choices`, `interior_ok`,
`interior_columns` and `search` in `gate/levels.jl`.

- [ ] **Step 4: Bring the tests to the rule**

Focused runs of the four test files of Task 10. A fixture rejected before and admitted
now is the spec's first allowed change: the test then asserts
`all(_testhelper_admission_agreement(am; rng, corners = true))` (`corners = false` above
8 steps) and its law against the limit. A fixture rejected before and rejected now for
another fault, or with another message, is the third. Append each to
`gate/changed-tests.md`. **Stop and report** if a fixture admitted before is rejected, a
newly admitted fixture fails an oracle, or a `MECHANISM_TEST_SPECS` entry changes verdict.
`test/test_boundaries.jl`, `test/test_type_stability.jl` and `test/test_compile_budget.jl`
pass. Report the enumeration file's failing testsets and wall time before and after.

- [ ] **Step 5: Commit**

`git commit -m "Search for a reading of the excluded steps that gives a limit"`

### Task 13: The gate for the energy-level rule

**Files:**
- Modify (outside the repository): `gate/populations.jl`, `gate/compare.jl`,
  `gate/report.md`; keep the last run's files as `*-before-energy-levels.*`
- Test: `test/Derivation/test_admission.jl`

**Interfaces:**
- Consumes: `D._admitting_reading`, `_testhelper_admission_agreement(am; rng, corners)`.

- [ ] **Step 1: Widen the populations.** To the four of `compare.jl` add, each recorded
with its fingerprint: the 50 ping-pong taggings of `gate/final.jl`; every nonempty
`:OnlyA` set on the plain init mechanisms of a bi-bi reaction with a competitive
inhibitor, of a ter-bi and of a bi-ter reaction (a spaced 3,000 of each when there are
more); one expansion round of the census reactions (uni-uni, bi-bi, ping-pong), the
allosteric children as emitted; and the distinct allosteric mechanisms
`test/Derivation/test_admission.jl` puts through the guard (the spike's way of collecting
them is in `lean-mechanism-spike/`). "Before" for the new populations is the guard at
9850682d, computed in a worktree at that commit.

- [ ] **Step 2: Write and run the checks**, each a line of `report.md` with its counts:

1. No mechanism admitted at 9850682d is rejected.
2. Every mechanism of 8 steps or fewer admitted now and not at 9850682d gives
   `all(_testhelper_admission_agreement(am; rng, corners = true))`.
3. Of the larger newly admitted ones, a spaced 150 per population with `corners = false`.
4. Every mechanism whose admitting reading has, on some face, a vanishing transition
   alive while its reverse is dead is listed with that count per population; each of 8
   steps or fewer is in check 2.
5. With the filling check switched off (a patched session, not a flag in `src/`), check 2
   on the `:onlyA` population must fail for at least 14 mechanisms, entries 6580, 6595 and
   8524 among them.
6. With the set-aside rule switched off, every verdict outside the mechanisms it admits
   is unchanged.
7. At two catalytic subunits: for a spaced 40 of the newly admitted mechanisms of 8 steps
   or fewer, the verdict equals the one-subunit verdict, and where
   `_testhelper_protomer_network` applies its steady state agrees with the law.
8. Reported, not required: how many of the mechanisms lost between the branch base and
   9850682d are admitted again; the slowest guard call and the guard's total time per
   population against 9850682d; how many mechanisms are too large to search.

- [ ] **Step 3: The long run.** In the background, alone: check 2's assertion with
`corners = true` on a spaced 100 of the newly admitted ter-ter mechanisms, under
`timeout` of 8 hours, writing each result as it finishes. Report what completed.

- [ ] **Step 4: The permanent sampled test.** A testset
`"sampled taggings: the energy-level rule agrees with the oracles"` beside Task 4's, with
a comment that it is an aggregate check over a population: from the ping-pong and bi-bi
init mechanisms of 8 steps or fewer, `MersenneTwister(20261120)`, 40 random nonempty
`:OnlyA` sets drawn among those the guard admits only by a searched reading; each gives
`all(_testhelper_admission_agreement(am; rng, corners = true))`. It asserts that at least
5 of the 40 have an admitting reading with a `:loses_reverse` entry. Under two minutes.

- [ ] **Step 5: Commit** the test. **If any of checks 1, 2, 3, 5, 6 or 7 fails, stop and
report the mechanism, its reading and both verdicts.**

`git commit -m "Pin the energy-level rule against the oracles on sampled taggings"`

### Task 4: The gate

**Files:**
- Create (outside the repository): `/home/denis.linux/claude-scratchpad/lean-mechanism-spike/gate/compare.jl`,
  `.../gate/report.md`
- Test: `test/Derivation/test_admission.jl`

**Interfaces:**
- Consumes: `gate_populations()` and `baseline.jls` from Task 1; the guard of Tasks 2 and 3.

- [ ] **Step 1: Write `compare.jl`**

For every entry of `gate_populations()`, build the mechanism and take the guard's fault
now. Classify against the baseline: admitted before and now; admitted before, rejected
now; rejected or refused before, admitted now; rejected both. Then:
- every "admitted before, rejected now" entry must have fault `_FAST_STEP_PINNED`, except
  entries 53145 and 53246 of the recorded populations, which must have `_BRIDGE` and hold
  an excluded rapid-equilibrium step with `step_kinds[k] === :fused_binding` that
  `_switch_offs` reads as `_LOSES_AFFINITY`; the report lists both;
- every "admitted now, not before" entry of `:onlyA` and `:terter`, and 400 spaced ones of
  `:onlyI`, must give `all(_testhelper_admission_agreement(am; rng))`;
- 200 spaced "admitted before and now" entries must give the same;
- for every `:onlyI` entry whose two conformations both turn the reaction over
  (`has_turnover_cycle` on `active_steps` and on `inactive_steps`), the entry and its
  mirror (`:OnlyA` and `:OnlyI` swapped) must get the same verdict.
Write the counts per population as a table to `report.md`, with the reference counts from
the spec's evidence table beside them (`:onlyA`: 6,748 admitted before, 258 lost, 1,030
gained).

- [ ] **Step 2: Run it**

Expected: no violation. Counts for `:onlyA` equal the spec's. **If a violation is found,
stop: do not go on to Task 5. Report the mechanism and both verdicts to the maintainer.**

- [ ] **Step 3: Add the permanent sampled test**

A testset `"sampled taggings: the rule agrees with the mass-action oracles"` in
`test_admission.jl`, with a comment that it is an aggregate check over a population. From
the plain init mechanisms of bi-bi and ping-pong, `MersenneTwister(20261119)`: 40 random
mechanisms with one random nonempty `:OnlyA` set each and 20 with a random state vector
holding an `:OnlyI`, built unguarded. For each: when the fault is `_NO_FAULT`,
`all(_testhelper_admission_agreement(am; rng))`; when it is `_FAST_STEP_PINNED`,
`any(c -> any(k -> _testhelper_switch_off(am, c, k) === :cannot_switch_off, steps), 1:2)`
and `_testhelper_epath_verdict(am; rng) === :no_path`. Assert that at least 10 of the 60
are admitted and at least 3 hold the new fault, so the sample exercises both.

- [ ] **Step 4: Run and commit**

Focused run of `test/Derivation/test_admission.jl`; the new testset takes under a minute.

```bash
git add test/Derivation/test_admission.jl
git commit -m "Pin the guard's rule against the oracles on sampled taggings"
```

### Task 5: What Enumeration does not build (runs after Task 13)

**Files:**
- Modify: `src/Enumeration/Enumeration.jl` (`_guard_accepts`)
- Modify: `src/Enumeration/filters.jl` (new `_kept_forms_can_be_made`)
- Test: `test/Enumeration/test_mechanism_enumeration.jl`

**Interfaces:**
- Produces, internal to `Enumeration`:
  `_kept_forms_can_be_made(m::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}) -> Bool`:
  true for a plain mechanism; for an allosteric one, false when in some conformation a
  form other than free E is touched by a step of that conformation's law
  (`active_steps`, `inactive_steps`), has at least one step into it that takes up one
  metabolite and gives off none, and every such step is excluded there
  (`excluded_active`, `excluded_inactive`).
  `_guard_accepts(m) = Derivation.check_valid_rate_equation(Bool, m) && _kept_forms_can_be_made(m)`.

- [ ] **Step 1: Write the failing tests**

A testset `"the enumerator does not build a form nothing can still bind to make"`, each
mechanism inline:

```julia
En = EnzymeRates.Enumeration
# S binding steady-state :OnlyA, the rest :EqualAI: the guard admits it, E(S) stays.
#   E + S <--> E(S) :: OnlyA;  E(S) <--> E(P) :: EqualAI;  E + P ⇌ E(P) :: EqualAI
@test EnzymeRates.Derivation.check_valid_rate_equation(Bool, shut_door)
@test !En._kept_forms_can_be_made(shut_door)
@test !En._guard_accepts(shut_door)
# The chemistry :OnlyA too: E(S) is cut off.
@test En._kept_forms_can_be_made(k_type) && En._guard_accepts(k_type)
# A merged complex still made by a product binding stays.
#   E + A ⇌ E(A);  E(A) + B <--> E(P, Q) :: OnlyA;  E(Q) + P ⇌ E(P, Q);  E + Q ⇌ E(Q)
@test En._kept_forms_can_be_made(merged_v_type)
@test En._kept_forms_can_be_made(EnzymeRates.Mechanism(plain_uni_uni))
```

- [ ] **Step 2: Run to verify they fail**

Expected: `UndefVarError: _kept_forms_can_be_made`.

- [ ] **Step 3: Implement both functions**

- [ ] **Step 4: Bring the enumeration tests to the rule**

Run the focused command on `test/Enumeration/test_mechanism_enumeration.jl`. For each
failing exact-children test, compute the children lost and gained. A lost child must hold
`_FAST_STEP_PINNED`; a gained child must give
`all(_testhelper_admission_agreement(child; rng, corners = true))` (`corners = false`
above 8 steps) and pass `_kept_forms_can_be_made`.
Then update the expected set and count, keeping every expected child written out. For an
aggregate count pin, record old and new counts. Append all to `changed-tests.md`. Any
other difference is a stop condition.
A lost child may also be one the guard still admits and the enumerator no longer proposes
(its tagging base is rejected or fails the filter, so the search builds on another base):
record it with the base that replaced it. A move's raw output may lose a child the guard
rejected before and rejects now, because its base is gone: record it as "rejected before
and now; base gone".
Reference from the spike: one expansion round of the census ping-pong reaction went from
1,593 to 1,559 children, bi-bi from 2,060 to 2,048, uni-uni unchanged.

- [ ] **Step 5: Run to verify they pass, and commit**

```bash
git add src/Enumeration/Enumeration.jl src/Enumeration/filters.jl \
        test/Enumeration/test_mechanism_enumeration.jl
git commit -m "Keep the enumerator from building a form nothing can still bind to make"
```

### Task 6: The rest of the suite

**Files:**
- Test: every other file under `test/` that fails

- [ ] **Step 1: Run the full suite in the background and read all of its output**

- [ ] **Step 2: Bring each failing test to the rule**

By the four allowed changes only, appended to `changed-tests.md`. Stop conditions from
Global Constraints apply, in particular: a `MECHANISM_TEST_SPECS` entry the macro now
refuses, and an identify test whose true or expected mechanism is now rejected. For an
identify test whose selected mechanism changed without either, record both selections in
`changed-tests.md` and report it at the end.

- [ ] **Step 3: Run the guards**

Full suite again; then `julia +1.13 --project -e '...'` focused run of
`test/test_type_stability.jl`, and the focused run of `test/test_compile_budget.jl` and
`test/test_boundaries.jl`. Expected: all pass, output pristine.

- [ ] **Step 4: Commit**

Stage the changed test files by path.
`git commit -m "Bring the remaining tests to the guard's rule for excluded steps"`

### Task 7: Documentation, roadmap and version

**Files:**
- Modify: `docs/src/deriving/mwc_allostery.md`
- Modify: `docs/src/developer.md` (the guard and admission sections)
- Modify: `ROADMAP.md`, `docs/src/roadmap.md` if it lists a changed entry
- Modify: `Project.toml` (`version = "0.11.0"`)

- [ ] **Step 1: Rewrite the user section**

In `mwc_allostery.md`, replace the bullets under "Tags as limits, and admission" with a
subsection titled exactly `### How an excluded step is switched off`. It holds, in the
user's terms and in this order: the assumption (an excluded step has a very high
dissociation constant or very low rate constants, and every thermodynamic constraint holds
in each conformation); that a step keeps its kind in both conformations, with the
maintainer's reason (rapid equilibrium and steady state differ more than 100-fold); the
energy picture and the four-row table of the spec's decision 2, and under them, in two
or three sentences, that the package looks for a choice of energies under which the rate
equation is the limit and accepts the mechanism when one exists; the spec's four worked
examples as `@allosteric_mechanism` blocks, the ping-pong among them; that the last check
is validated numerically, with what the gate tested; what
the error for a rapid-equilibrium step inside one group means and its two fixes (make the
step steady-state, or exclude the steps that keep its two sides joined).
Then correct what the page says elsewhere: the paragraph at lines 66 to 70 (any group may
be `:OnlyI` while the active conformation catalyzes), the numbered "Existence" and
"Uniqueness" items and the example after them, and the K-system paragraphs near line 280
(the catalytic step must be `:OnlyA` too when the substrate binding is at rapid
equilibrium; a steady-state binding is stopped instead).

- [ ] **Step 2: Rewrite the developer page's guard section** to the energy-level rule:
readings and their three conditions, the search with the three built readings first and
its two size limits, the filling check and the set-aside rule with their arguments, and
the faults, with no mention of the rule before this work.

- [ ] **Step 3: Update `ROADMAP.md`**

Rewrite "Let an allosteric ping-pong stop both half-reactions in the inactive
conformation" for the rule as built: build the examples under "What the search does
today" and record what the guard and the enumerator now do with each; turn "The
energy-level rule" bullets into what remains open (the ping-pongs with no limit under any
reading, the valid mechanisms still rejected and why, the wider set-aside variant that
was not built, mechanisms too large to search, and that the last check rests on the
gate's numerical evidence); move the entry out of Planned if nothing planned is left. Add under Planned the two things the spec leaves unchecked or rejected for
now: ter-ter beyond the gate's sample. Add a Known issue for the valid mechanisms the
guard still rejects (too large to search; refused as a bridge; refused by the filling
check though the rate equals the law's), with the gate's counts. Update any Known issue that names `_NO_LIMIT` or the plain-binding rule.

- [ ] **Step 4: Bump the version and build the docs**

Run: `julia --project=docs docs/make.jl`
Expected: the build and every doctest pass.

- [ ] **Step 5: Commit**

```bash
git add docs/src/deriving/mwc_allostery.md docs/src/developer.md ROADMAP.md Project.toml
git commit -m "Document how an excluded step is switched off; bump to 0.11.0"
```

(Add `docs/src/roadmap.md` if it changed.)

- [ ] **Step 6: Report**

Give the maintainer `gate/report.md` and `gate/changed-tests.md` for the pull request
description, with any identify selections that changed. Do not open or merge a pull
request unless asked.
