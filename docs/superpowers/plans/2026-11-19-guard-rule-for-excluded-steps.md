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

- The rule, copied from the spec. In a conformation, for each step it excludes:

  | Excluded step | How it is switched off |
  |:--|:--|
  | Rapid equilibrium, takes up one metabolite into a cut-off form | Loses its affinity: stays fast, uptake vanishes |
  | Rapid equilibrium, anything else | Fault |
  | Steady state, takes up one metabolite into a cut-off form | Loses its affinity: uptake rate vanishes, release keeps its rate |
  | Steady state, anything else | Stopped: both rates vanish |

  "Takes up one metabolite into a cut-off form" means: `length(e.taken_up) == 1`,
  `isempty(e.given_off)`, the steps the conformation keeps do not join `e.from` and
  `e.to`, and they do not join `e.to` to free E (form 1).
- `Derivation` reads no `step_kinds` and no form contents for this rule. Other reads of
  `step_kinds` in `src/Derivation/` (names, thermo, engine) stay until stage 2.
- No `public` statement and no export changes. The module-map test must still pass.
- A test may change only in the four ways the spec lists under "Verdicts that change".
  The tests of a law that is not the limit pass unchanged.
- **Stop and report to the maintainer** when: a mechanism admitted before this work is
  rejected for any reason but the new fault; a newly admitted law fails an oracle; a
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
- every "admitted before, rejected now" entry must have fault `_FAST_STEP_PINNED`;
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

### Task 5: What Enumeration does not build

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
`all(_testhelper_admission_agreement(child; rng))` and pass `_kept_forms_can_be_made`.
Then update the expected set and count, keeping every expected child written out. For an
aggregate count pin, record old and new counts. Append all to `changed-tests.md`. Any
other difference is a stop condition.
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
four-row table; the spec's three worked examples as `@allosteric_mechanism` blocks; what
the error for a rapid-equilibrium step means and its two fixes (make the step
steady-state, or exclude the steps that keep its bound form reachable).
Then correct what the page says elsewhere: the paragraph at lines 66 to 70 (any group may
be `:OnlyI` while the active conformation catalyzes), the numbered "Existence" and
"Uniqueness" items and the example after them, and the K-system paragraphs near line 280
(the catalytic step must be `:OnlyA` too when the substrate binding is at rapid
equilibrium; a steady-state binding is stopped instead).

- [ ] **Step 2: Rewrite the developer page's guard section** to the rule and the new
fault, with no mention of the old one.

- [ ] **Step 3: Update `ROADMAP.md`**

Re-examine "Let an allosteric ping-pong stop both half-reactions in the inactive
conformation": build its examples, record what the guard now says of each, and rewrite the
entry to match. Add under Planned the two things the spec leaves unchecked or rejected for
now: an excluded rapid-equilibrium isomerization that could stay fast, and ter-ter beyond
the gate's sample. Update any Known issue that names `_NO_LIMIT` or the plain-binding rule.

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
