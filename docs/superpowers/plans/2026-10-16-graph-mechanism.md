# Graph Mechanism (Phase 3) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land phase 3 of the reaction-graph rewrite on `graph-mechanism`: `Mechanism` and
`AllostericMechanism` become the reaction graph, both macros speak the modifier grammar,
`show` prints the macro text, an eval-free parser reads it back, and identify's CSVs store
that text and its id.

**Architecture:** The plan replays the reviewed prototype (branch `proto3`, 233 commits,
`4f9a34a..7c2a657`) onto `graph-mechanism` at `7205ba4` (R60). Each task applies its
prototype commits from the shared object store, test paths first, so the red run comes
before the source and the green run after it; code is written out only where the replay
departs from the prototype. Each stage ends at a checkpoint that diffs the branch against
the prototype's stage tip.

**Tech Stack:** Julia 1.12 (local runs), 1.13 (CI), 1.10 (compat floor); Test, TestEnv,
Aqua, JET 0.9/0.11/0.12; the SHA and Unicode stdlibs; git; the phase-3 census tools.

**Spec:** `/home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/graph_mechanism_wt/docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`
(phase-3 sections 3, 4, 5, 13-16). The rulings in `$P/RULINGS.md`, R1 to the last, settle
what it leaves open; the prototype's reports (`$P/proto/S1_report.md` ... `S5_report.md`)
say what each commit does.

**Paths** (every command below uses them):

- `R=/home/denis.linux/claude-scratchpad/EnzymeRates-12f42995`
- `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad` (a symlink to `$R`)
- `P=$S/phase3`; `W=$R/graph_mechanism_wt` (the target worktree, branch `graph-mechanism`;
  every command runs in `$W` unless it says otherwise); `L=$P/replay/logs`
- The prototype worktree `$P/proto3_wt` is read only: never edit, commit, reset or check
  out there. `git show <sha>` reads its commits from `$W`, since both share one object store.
- Abbreviations: FX = `test/mechanism_definitions_for_test_enzyme_derivation.jl`, TT =
  `test/test_types.jl`, TD = `test/test_dsl.jl`, TE = `test/test_mechanism_enumeration.jl`,
  TRD = `test/test_rate_eq_derivation.jl`, TRG = `test/test_reaction_graph.jl`, TJS =
  `test/test_joint_solve.jl`, TO = `test/test_oracles.jl`, TID =
  `test/test_identify_rate_equation.jl`, TCB = `test/test_compile_budget.jl`.

## Task index

Units: each entry is one commit, the main prototype commit and the number of commits folded
into it (`+n`; a piece of a commit the map splits counts as one). ROADMAP entries written
out as departures (Tasks 34, 36, 42, 48, 54, 60) are not counted.

| Task | Stage | Title | Units |
|---|---|---|---|
| 0 | setup | Commit the plan and prepare the replay | none |
| 1 | S5 test environment | Julia 1.13 and 1.10: JET 0.12 and m_OnlyA_prod's peak turnover | 5c08ede, 4a96d9a +4 |
| 2 | S1 | Check the laws behind the string pins phase 3 renames | a3de3c4 +3 |
| 3 | S1 | C9: orient tier-3 steps and tied leads by form index | 2bd4b5d +2, 43d7af8 |
| 4 | S1 | Render phase 3's constant names | 8b7fb88 +10, c86f61c |
| 5 | S2a | Raw input types | 8d5e465 +1 |
| 6 | S2a | One walker reads every new-grammar block into raw input | b3eabc4 +3 |
| 7 | S2a | oligomeric_state tuples and a metabolite named M | 7cf6b0f +1, 6075bba |
| 8 | S2a | Reserved metabolite names and names that hold _ | 51fad19 +2 |
| 9 | S2a | Reject a reaction whose names spell one name two ways (R38) | 410f7c0 +4 |
| 10 | S2b | Write the allosteric test fixtures inline | c80e857 +1 |
| 11 | S2b | Migrate the tests and docs to the modifier grammar | c2cd9b3 |
| 12 | S2b | Hand-migrated tests | bb9b2e3, b0b3a4e, 96846b8, 9e6e64c, 9227cbe, bc804d6 |
| 13 | S2b | Describe the modifier grammar in titles, comments, docstrings and docs | 3c0c040, 588e99f, fcfd5e6, 9e41af7 +1, 756abee |
| 14 | S2c | Reject the inhibitor grammar's mechanism forms; read every block with the walker | 3b28fca, 68a4961, fa09735, 9399c58 +1, 822d590, 9736179, b361684 |
| 15 | S2c | @enzyme_reaction: one form per label | 961ec20, 37c6a8b +1, f8d0199, 6c85d75 |
| 16 | S2c | Messages in the modifier grammar | 59136b3, c457930, 8fc17da +1, 257e81b, be61164 |
| 17 | S2c | Print a reaction as its block; summary modifiers; the documented grammar | 2493151 +1, 67dc293 +2, 3c44dd3, 7adea3e, b3431bf, 60bd37a |
| 18 | S3a | A mechanism's identity and names without its reaction | b588424, db5b3ff, 10527e7 |
| 19 | S3a | Labels and metabolite names from the mechanism | 31ff003, 6ba4fed |
| 20 | S3a | The boundary check | 2b79bd3 +2 |
| 21 | S3a | Naming without a cache; Symbol ligands | f7bd003, 4d13802 |
| 22 | S3b | Position-picked fixtures written out; populations sampled by content | ad4aca2, cd380a2 |
| 23 | S3b | Steps, groups and forms in content order (rows 1-2) | 1362b5b, 9916684 |
| 24 | S3b | Segment roots, cycle rows and form identity by content (rows 3-5) | 77b7b61 +1, eb508a5, 9c6d4cb |
| 25 | S3b | Enumeration ties by content (rows 6-9) | 1a33d0f, 7f1aa95, af01ff6, 0266794 |
| 26 | S3b | Seed-base order, the R32 Idea, forms that share an internal name (rows 10-11) | 74a8fc6, 61a1256, a2bbf86 |
| 27 | S3c | Tests read mechanisms through one reader | 8107bc0, 990cb6a, 279855d, f57370e, 0e051b4 |
| 28 | S3c | R48's test items | 56f7871, bdd16e2, 1cd8cd5, d52040f, 380a8c4, e080b63 |
| 29 | S3d1 | The reaction graph's value types, content order and identity | e12daba |
| 30 | S3d1 | Canonical form and checked constructors | 2db9d74, 584d7bf +1, 43c50fc |
| 31 | S3d1 | Sig, compiled types, lifts and _raw | 6177781 +1, f8eb0df |
| 32 | S3d1 | The graph against the Step types; the constructor's speed gate | bf41ad0, 6ca71bd, 1bf2fee, 03fad30, 94eaa5d, 5298b84 |
| 33 | S3d2 | Before the swap: tier-3 texts, names.jl, the pivot block | f7c1587, 6214bb1, b0a4b29 |
| 34 | S3d2 | The swap | f1aa056 +10 |
| 35 | S3d2 | After the swap: the bottomless check typed, views built once | 92ba4f1, 86d21f6, 1894d39 +1 |
| 36 | S3e | Row 25: free E | 85e0790 |
| 37 | S3e | Rows 26-28: one connected graph, every declared name used | b9c535e, 450a713 +2 |
| 38 | S3e | Rows 29-31: names and labels in the constructors, R38 at the macros | a76ddf2 +3, 2aa3cbd +2 |
| 39 | S3e | Multiplicities, typed inhibitors, elements, Bool | ad74058 +1, e7698be, 56dbf00, 5b72e27 |
| 40 | S3e | Test cleanups, tighter Sig gates, the review's pins | c6d5a47, 9cbb9af, 44661cc, dfc5f79 +1, ff5e59c |
| 41 | S3f | The boundary check on the graph | f3d5e2c |
| 42 | S3f | Constants, the joint solve and the derivation over edges | 19e3fe8 +1 |
| 43 | S3f | Decode once per generator; reject labels A and I | 43519c4, a75ef55 +1 |
| 44 | S3g | Segment sides and cycle rows; A/I pins | d427228, ee46733 |
| 45 | S3g | Flip on the graph | 8752036 |
| 46 | S3g | Copy rule and dead-end move on the graph | 0f795c4 |
| 47 | S3g | Split on the graph | 5d1197d +1, 7d69b90 |
| 48 | S3g | Tag, relax and check on the graph | 2c40783 +4 |
| 49 | S3g | Hotspot fixtures inline; the census population freed | 612a678 +1, 698fd6a |
| 50 | S3h | Seeds from raw input | 9dee115 |
| 51 | S3h | Seed variants and degeneracy on the graph | 2395401 +1 |
| 52 | S3i | The compiled types' summary from the graph | 6ada4b2 |
| 53 | S3i | Every test reads the graph | 36384ef +4 |
| 54 | S3i | Delete the Step views, Step, Species and Residual | 2692bb7 +3 |
| 55 | S3i | Minors: raw iso test, seed bases, merged complexes | 14ef525 +4 |
| 56 | S4 | Count a repeated element once | b784fdc |
| 57 | S4 | @enzyme_reaction resolves its names; the tag hints name a bare metabolite | a520f8e, 3f6c43c +5 |
| 58 | S4 | Reject names Julia reads as another name | 13f202c +1 |
| 59 | S4 | Print a mechanism as its macro text | 3464d63 +4 |
| 60 | S4 | Read a mechanism back without eval | 0ec7d51 +5 |
| 61 | S4 | identify's CSVs store the text and its id | e59fb5b +5 |
| 62 | S5 | CLAUDE.md facts for Denis | ea6bf86 |
| 63 | S5 | Final: the 1.10 Known issue, suites, census, docs, PR | 7395164 +2 |

## Global Constraints

- Version stays `0.10.0`; no task touches `version` in Project.toml.
- Julia compat stays `julia = "1.10"`; CI runs Julia 1 (1.13). The final task runs the full
  suite on 1.12 and 1.13, and on 1.10 when at least 5 GB is available (R64).
- `rate_equation`: 0 allocations and under 120 ns per call for every `MECHANISM_TEST_SPECS`
  mechanism (`_testhelper_test_rate_equation_performance`, TRD). Any task that slows it:
  stop and ask Denis.
- R18: no function that callers loop over returns a Union of immutable, reference-holding
  parametric structs; never read, in a loop, a field whose type differs across a Union's
  members. The final task reruns the remark scan (`--pkgimages=no`) on derivation and
  identify workloads on 1.12 and 1.13.
- R52 (gate at Task 55, the end of S3i; reported from Task 34 on): full-suite maxrss at most
  3.1 GB, ter-ter `init_mechanisms` at most 31 s, `test/test_mechanism_enumeration.jl` alone
  at most 145 s.
- R17: the median `Mechanism(raw)` time over the specs stays under 10 µs (TRG's gate from
  Task 32 on; the prototype measures 3.2-3.8 µs).
- Runner: Julia runs only through `$P/jl.sh <log> <timeout_s> <julia args...>`, one process
  at a time (7.7 GB, no swap); logs go to `$L/T<NN>_<what>.log`. Bash calls time out at
  600 s: start a long run with `run_in_background: true` and poll with
  `timeout 580 bash -c 'until grep -q "^exit " LOG; do sleep 20; done'; tail -40 LOG`.
- Focused run: `focused <log> <test file>...` from `env.sh` (Task 0 Step 2), or the Gt and
  Rt runners of Tasks 44-63, which run the same Julia code on file names relative to
  `test/`. Either runs CLAUDE.md's focused command with each named file in a testset of its
  own, so a failing file does not stop the files after it and the summary prints one count
  per file. No task includes the test files one after another.
- Counts: an expected count is a report's figure at the same tree (with the offsets the
  Replay procedure names) or derived from the commits' added and removed assertions; each
  task says which. One rule holds for every task and both kinds: a count that differs from
  its expectation, whether or not anything fails, stops the task. Check that every path of
  the unit's commits applied once and in order (`git show --stat --format= <sha>` for each
  commit against `git diff --stat` of the unit). A path missing, doubled or out of order:
  redo the unit from its first apply. Every path applied once: the expectation was wrong;
  record the measured count, carry the difference into the task's later expectations,
  report both to the controller (DONE_WITH_CONCERNS) and go on. Never change a test to meet
  a count.
- Shell helpers: every shell block first sources `$P/replay/env.sh`, which Task 0 Step 2
  writes. It sets the Paths, leaves the shell in `$W` and defines `focused` (the focused
  run, each test file in a testset of its own, so a failing file does not stop the files
  after it), `full`, `wait_log`, `summary`, `pristine`, `census_run`, `census_check` and
  `chg`. The Runners section defines the shorthands the steps use.
- Full suite: `$P/jl.sh $L/T<NN>_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`;
  output pristine (no warning or error beyond the base's 12 `[ Info:` lines).
- Census (from `$P`; corpus about 2.5 min, compare 5 s):
  `./jl.sh census/logs/replay_T<NN>.log 1800 --project=census/env_replay --heap-size-hint=2500M census/corpus.jl $W census/replay_T<NN> census/rename_phase3.toml`,
  then `./jl.sh census/logs/cmp_replay_T<NN>.log 600 --project=census/env --heap-size-hint=2500M census/compare.jl census/<ref> census/replay_T<NN> [census/rename_phase3.toml]`.
  The reference and the map per task are in each task's gate. Tasks 34-53 prefix the corpus
  command with `CENSUS_ADAPTER=neutral_graph.jl` (the Step views exist beside the graph);
  from Task 54 `corpus.jl` picks it alone.
- Commits: one per commit unit (below). Message: the unit's main commit's subject and body
  (`git show -s --format=%B <sha>`), edited only where CLAUDE.md's evergreen rule or a fold
  requires (a fold that changes what the commit does adds its sentence to the body). Every
  message ends with the two lines
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE`
  (45 prototype commits, the last of them 8107bc0, lack the second; add it).
- No whitespace edits by hand, except R55 m2 (Task 49, listed there).
- Push only in the final task; open a draft PR against `main` and never merge it. The PR body
  is `$P/proto/drafts/pr_body_phase3.md` with the edits Task 63 Step 11 writes out, which
  make it true of this branch (the prototype's hashes, when the memory rose and fell, the
  rulings' range, run times from Task 63's own logs); R64 approved them, on the draft
  that carries R62 r1's and R66 m2's wording. No session link.
- Every CLAUDE.md rule applies to anything written out: ABOUTME headers, `_testhelper_`
  prefixes, 92 columns, evergreen comments, no mocks.

## Runners

Every task's steps run through the helpers and shorthands below; an implementer gets this
section, the Global Constraints and the Replay procedure with its task.

**The shell.** The Bash tool keeps no shell state between calls, so every shell block from
Task 0 Step 3 on starts by sourcing `env.sh`, which Task 0 Step 2 writes (Task 0's
Interfaces list its helpers). It sets the Paths (`R`, `S`, `P`, `W`, `L`), leaves the
shell in `$W`, and defines `focused <log> <test files>`, `full <log>`, `wait_log <log>`,
`summary <log>`, `pristine <log>`, `masked <log> <reference log>`, `census_run <NN> <ref>
[map]`, `census_check <NN>` and `chg`; a `<log>` without a `/` names `$L/<log>.log`. Prose
writes the source line `. .../replay/env.sh`; in a command it is
`. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh`.
"Background: `. .../replay/env.sh; X`" means one call with `run_in_background: true` that
runs X after that line, written out in full; a runner written in prose (F, C, SC, PR, CD,
RS) likewise runs after it in the same call. Commands run in `$W` unless they start with
`cd`.

**Background runs.** Every Julia run goes through `$P/jl.sh`, one at a time (7.7 GB, no
swap). Start every focused, named, full-suite, census, seed-count, TE-alone, JET, docs,
version-suite, measurement and archive run with `run_in_background: true`, then call
`wait_log <log>` until the log ends with an `exit` line; for any other log, poll
`timeout 580 bash -c 'until grep -q "^exit " LOG; do sleep 20; done'; tail -40 LOG`.

**Shorthands.** Tasks 14-28 and Tasks 44-63 run their focused runs differently, so their
green and red runners carry different names: G and Rd, Gt and Rt.

- **G `<log>` `<files>`** (green run, Tasks 14-28): background `. .../replay/env.sh;
  focused <log> test/<file> ...`, one `test/` path per file the step names, then
  `wait_log <log>` and `summary <log>`. `focused` runs each file in a testset of its own
  inside one outer testset, so a failing file does not stop the files after it. Passed:
  under `focused`, one line per file whose Pass equals its Total and the count the step
  gives (`test/test_dsl.jl | 233  233`), no Fail or Error column, and the log's last line
  `exit 0 elapsed <s>`.
- **Rd `<log>` `<files>`** (red run, Tasks 14-28): the same call. The log ends `exit 1`.
  The step names every testset expected to fail and how. Read the `Test Failed` and
  `Error During Test` blocks: a failure or error in any other testset means a step was
  skipped or reordered; stop and compare with the task's commit list.
- **Gt `<log>` `<files>`** (green run, Tasks 44-63): each named file inside one top-level
  testset of its own, so each file prints one Test Summary line with its count:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh <log> 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); for f in ARGS; @testset "$f" begin include("test/$f") end; end' <files>
```

  Passed: one line per file whose Pass equals its Total and the count the step gives
  (`test_dsl.jl | 262 262`), and the log's last line `exit 0 elapsed <s>`.
- **Rt `<log>` `<files>`** (red run, Tasks 44-63): Gt with `try ... catch; end` around
  each file's testset, so every testset of every file runs and the log ends `exit 0`:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh <log> 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); for f in ARGS; try; @testset "$f" begin include("test/$f") end; catch; end; end' <files>
```

  The step names every testset expected to fail and how. Count them with
  `grep -E '(Test Failed|Error During Test) at' <log> | sed -E 's/: (Test Failed|Error During Test) at .*\/(test_[a-z_]+\.jl):[0-9]+/ | \1 | \2/' | sort | uniq -c`.
  A failure or error in any other testset means a step was skipped or reordered: stop and
  compare with the task's commit list.
- **N `<log>` `<file>` `"<title>"` ...** (named run): evaluates each file's top-level
  definitions, then only the named testsets, and prints `PASS` or `FAIL: <message>` per
  testset under `==== run_named summary ====`; `<log>` is a path (`$L/T20_red.log`) and
  `<file>` the path from `$W` (`test/test_dsl.jl`). Used where a unit's test edits also
  call the API the unit introduces, so a whole-file red would error far past the tests
  the unit adds; the prototype took these reds this way (`$P/proto/scripts/S3a_named.sh`).
  Background, then `wait_log <log>`; Tasks 29-43 write the command out:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh <log> 1800 --project=. --heap-size-hint=2500M $P/s3c/scripts/run_named.jl <file> "<title>" [<file> "<title>" ...]
```

- **F `<NN>`** (full suite; Tasks 0-13 and 29-43 write `full T<NN>_full` out): background
  `. .../replay/env.sh; full T<NN>_full` (the Global Constraints' command,
  `$P/jl.sh $L/T<NN>_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`),
  then `wait_log T<NN>_full` until it ends, `summary T<NN>_full; pristine T<NN>_full`.
  Passed: `EnzymeRates.jl | <count>  <count>`, `Testing EnzymeRates tests passed`,
  `exit 0`; `pristine` prints `12` and nothing else (the base's 12 `[ Info:` lines:
  `rate_equation` per-call figures and the compile budgets). Record the log's `maxrss_kb`
  line: R52's full-suite figure is that number / 10^6, in GB.
- **C `<NN>` [`<ref>`]** (census, always with the rename map; `<ref>` defaults to
  `S3b_end`): background `. .../replay/env.sh; census_run <NN> <ref>
  census/rename_phase3.toml` (the corpus of `$W` into `census/replay_T<NN>`, then
  `compare.jl census/<ref> census/replay_T<NN> census/rename_phase3.toml`), then
  `wait_log $P/census/logs/cmp_replay_T<NN>.log` until it ends, and `census_check <NN>`.
  From Task 34 through Task 53 the census goes through the graph adapter: prefix the call
  with `CENSUS_ADAPTER=neutral_graph.jl` (Tasks 34-43 write it out), since from the swap a
  mechanism is a graph while `EnzymeRates.Species` lives on as the Step views' type, so
  `corpus.jl` would pick `neutral_steps.jl`, which reads the reaction a graph no longer
  holds (census/README.md); from Task 54 `corpus.jl` picks `neutral_graph.jl` by itself
  and the variable is left out. "Zero differences": `census_check` prints `RESULT: PASS`
  and nothing else (init, seed and child sets equal and in order in all six reactions,
  every law string exact, values max rel err 0.0, no flipped stored step). Against
  `S3b_end` the corpus holds 81,498 mechanisms and 542 laws (S3g-S5 §5).
- **K `<sha>` `<paths>`** (commit a unit whose message is its main commit's, trailers
  added):

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short
git add <paths>
git show -s --format=%B <sha> | git interpret-trailers --if-exists addIfDifferent --trailer "Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE" > $P/replay/msg
git commit -F $P/replay/msg
git status --short; git log --oneline -1
```

  Expected: the first status lists exactly `<paths>`; after the commit no status line, and
  the log line ends with the main commit's subject. The command copies the prototype's
  message exactly and adds the Claude-Session line where it is missing (the main commits
  from b588424 to 8107bc0; every later prototype commit carries both trailers). A unit
  whose message changes writes it out in full in its step and commits with
  `git commit -F $P/replay/msg` after `git add`. Tasks 0-13 write each commit out with
  its expected output. Tasks 29-43 write every commit out: stage the unit's paths by
  name, print `git status --short`, commit with the message written out, then print
  `git status --short; git log --oneline -1`. Expected: the first status lists only
  staged entries (`M `, `A `, `D `), each one of the unit's paths; after the commit no
  status line, and the log line ends with the subject. Every message is the prototype's,
  verbatim unless the unit says it is edited.
- **SC `<NN>`** (seed counts; Tasks 34-43 write "Seed counts (`<NN>` = n)"):

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh $L/T<NN>_seed_counts.log 1800 --project=. --heap-size-hint=2500M $P/proto/scripts/S3a_seed_counts.jl
grep -E '^(pp|bb|ub|tt):|^exit' $L/T<NN>_seed_counts.log
```

  Passed (as `$P/proto/logs/S3f_seed_counts_head.log`): `pp: 264 mechanisms in <t> s`,
  `bb: 239 mechanisms in <t> s`, `ub: 7 mechanisms in <t> s`, `tt: 250855 mechanisms in
  <t> s`, `exit 0`; the `tt:` time is R52's ter-ter init figure.
- **TEA `<NN>`** (the enumeration file alone, R52's third figure; Tasks 34-43 write "TE
  alone (`<NN>` = n)"):

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh $L/T<NN>_TE.log 1800 --project=. --heap-size-hint=2500M -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); t = @elapsed include("test/test_mechanism_enumeration.jl"); println("TE_ELAPSED ", t)'
```

  Passed: `TE_ELAPSED <s>`, `exit 0`. TE has no outer testset, so the log holds one
  summary per top-level testset and no line with the file's count (as
  `$P/proto/logs/S3f_TE_head.log`); a failing testset throws and ends the run with `exit
  1`. The file's count comes from the task's focused run and full suite.
- **R52 record `<NN>`** (Tasks 34-43), after the task's full suite, SC and TEA:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
grep -H maxrss_kb $L/T<NN>_full.log $L/T<NN>_seed_counts.log $L/T<NN>_TE.log
grep -h '^tt:' $L/T<NN>_seed_counts.log; grep -h TE_ELAPSED $L/T<NN>_TE.log
```

  R52's three figures: the full suite's `maxrss_kb` / 10^6 (R52's GB), the ter-ter init
  seconds, the TE seconds. R52 gates only at Task 55 (3.1 GB, 31 s, 145 s); earlier tasks
  report them beside the prototype's figures the task names. Those lack 050e16f's copies,
  which Task 34 folds into the swap, so lower memory is expected. Do not stop on a figure
  above the gate before Task 55.
- **JET on ReactionGraph `<NN>`** (Tasks 29-32): the suite's JET names `EnzymeRates` as its
  only target, so it skips the submodule until Task 34 moves the module into `types.jl`.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh $L/T<NN>_jet.log 1800 --project=. $P/s3d/scripts/jet.jl
grep '^targets' $L/T<NN>_jet.log
```

  Expected (as `$P/proto/logs/S3d_c1_jet.log`): `targets (EnzymeRates,
  EnzymeRates.ReactionGraph): 0 reports` and `targets (EnzymeRates,): 0 reports`.
- **PR `<log>` `<script>`** (a probe of the prototype's, run on `$W`):
  `$P/jl.sh <log> 900 --project=. --heap-size-hint=2500M <script>`.
- **CD `<NN>`** (copy-rule drops over the census parents):
  `cd $P && ./jl.sh $L/T<NN>_copydrops.log 1800 --project=census/env_replay --heap-size-hint=2500M s3g/scan/copydrops.jl $P/census census/S3b_end`.
  Passed (`$P/proto/logs/S3g_copydrops_head.log`): the `bibi` line and the last line read
  `["change_state" => 2, "dead_end.allo" => 44, "dead_end.plain" => 81, "other" => 86, "split.allo" => 13, "split.plain" => 27]`,
  every other reaction `Pair{String, Int64}[]`, `exit 0`.
- **RS `<NN>` `<env>`** (resolve an environment after a new stdlib dependency;
  `<env>` is `.`, `census/env_replay` from `$P`, or `$P/replay/docenv`):
  `$P/jl.sh $L/T<NN>_resolve_<name>.log 900 --project=<env> -e 'using Pkg; Pkg.resolve()'`.
  Passed: `exit 0`. The Manifests already list the stdlib (other packages use it); in
  an environment that develops `$W` the resolve records it as a dependency of
  EnzymeRates, in the `deps` line under `[[deps.EnzymeRates]]`.
- **ENT `<title>`** (a ROADMAP entry against its final text): prints nothing when the
  entry of `ROADMAP.md` that starts `- **<title>` equals the one in
  `$P/replay/roadmap_final.md`:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
entry() { awk -v t="- **$1" 'index($0, t) == 1 {p = 1; print; next} p && (/^- \*\*/ || /^## /) {p = 0} p' "$2"; }
diff <(entry '<title>' ROADMAP.md) <(entry '<title>' $P/replay/roadmap_final.md)
```

- **CK `<tip>` `<ops>`** (stage checkpoint, Tasks 44-63): `$P/replay/ckpt_p4.sh <tip>
  <ops>` builds the stage's expected tree in a scratch index (the tip, D1-D4 and the ops:
  `c:<sha>` a whole commit, `p:<sha>:<path>` one path of one, `e:<title>=<file|delete>` a
  ROADMAP entry, `i` R55 m2's indentation) and prints `git diff --cached --stat HEAD`.
  Passed: it prints nothing. Task 44 Step 1 writes it.

**Counts.** Julia 1.12; under `focused`, `summary` prints one line per test file. Every
count, measured or derived, follows the Global Constraints'
Counts rule: a count that differs from its expectation, whether or not anything fails,
stops the task until the unit's paths are checked against its commits; then redo the
unit, or record the count, report it and go on. Never change a test to meet a count.

- Tasks 0-13: the full-suite counts of Tasks 4, 9 and 13 are the reports' counts at the
  same tree plus Task 1's 2 (S1 r1 38,573, S2a r1 38,842, S2b r1 38,849); Task 0's is
  derived (S1 §1: 38,163 outside Aqua and JET, plus their 11); Tasks 1 and 12 add what
  their commits add (+2; +6 and +1); and Tasks 2, 3 and 5-8 carry the counts the fix
  round measured on the dry replay's trees of those tasks
  (`$P/plan/fix_counts/logs/T0<N>_full.log`, Julia 1.12, the prototype's Manifest; the
  same run gave Task 4's tree 38,575, the reports' figure). The per-file counts in the
  focused runs of Tasks 4 and 8-13, and the full-suite counts of Tasks 10 and 11, are
  derived from S2b §4's per-file table (taken at 62ab877, 7b86c41 and c34c525), on the
  assumption that each commit between those points keeps a file's count unless the task
  says otherwise.
- Tasks 14-28: full-suite counts are each task's gate; per-file counts are the
  prototype's (S2_report.md section 4, S3_report.md sections 1 and 4,
  `$P/proto/logs/S3c_counts_base.log`) moved by each commit's added assertions, with
  test_rate_eq_derivation.jl 2 higher than the prototype's throughout (Task 1's 68f46f0
  testset). Start of Task 14: test_dsl.jl 235, test_types.jl 618,
  test_mechanism_enumeration.jl 32,534, test_rate_eq_derivation.jl 3,184.
- Tasks 29-43: at Task 28's end: TT 662, TD 251, TRD 3,186, TO 662, TJS 453, TE 32,543,
  TA 155 (`test/test_accessors.jl`), TID 446, TCB 16, TFS 104
  (`test/test_fitting_scale.jl`), `test/test_pivot_priority_regression.jl` 2. TRD runs 2
  above the prototype throughout (Task 1's testset); TRG runs 45 above it from Task 34 on
  (050e16f's testset, which the prototype adds at S3i). A red step gives pass, fail and
  error per file where the prototype's log fixes them, adjusted for those offsets and for
  the folds the unit names. A failure or error outside the testsets a step names means a
  step was skipped or reordered: stop and compare with the unit's commit list.
- Tasks 44-63: at Task 44's start (Task 43's tree): test_mechanism_enumeration.jl 32,545,
  test_types.jl 694, test_reaction_graph.jl 956, test_rate_eq_derivation.jl 3,193,
  test_joint_solve.jl 455, test_oracles.jl 667, test_dsl.jl 257,
  test_identify_rate_equation.jl 446, test_compile_budget.jl 16, test_accessors.jl 155,
  test_allosteric_collapse.jl 36, test_fitting.jl 88. They are the prototype's counts at
  7fe5c6d (TRD, TRG, TD and TT 5,053 in `$P/proto/logs/S3f_c5_green.log`; TE from
  `S3g_c3a.log` and `S3g_red_c6.log`; the rest from `S3i_counts_6ada4b2.log`) plus Task
  1's +2 in TRD and 050e16f's +45 in TRG, which the prototype adds only at S3i. TRD stays
  2 above the prototype's to the end; TRG equals the prototype's from Task 52 on. Each
  step's count follows from the commits' added and removed assertions (S3g-S4 §4, S5 §4,
  S5 minors round 2 for ad3ca91, round 3 for eef3af7 and round 4 for 7c2a657).

## Review Focus

Inputs the spec implies that no task's tests exercise, most likely to bite first. The plan
replays the prototype exactly (R60), so it adds no test for them; the final task's
whole-branch review checks each, and Denis sees them in the PR body.

1. **identify on real worker processes.** `identify_rate_equation` ships `Mechanism`,
   `AllostericMechanism` and compiled values through `pmap` (docs/src/identify/parallel.md,
   SlurmClusterManager); no test starts a worker (`addprocs`), so R19's "structs serialize"
   is unchecked. Expected: the same CSV rows as one process.
2. **A `mechanism_id` read as a number.** `CSV.read(path, DataFrame)` without
   `types = String` turns an id of decimal digits, or digits around one `e`, into an
   `Int64` or `Float64`, so `parent_id` joins break; every test passes `types = String`
   (D-S4-11; PR question 5). Expected: an id reads back as written.
3. **REPL display of a compiled mechanism.** `display(m)` goes through `text/plain`, which
   no test calls; it prints the macro text but compiles Base's method once per Sig (D-S4-6,
   m9; PR question 6). Expected: the same text as `show(io, m)`, with no compile per
   mechanism in a loop.
4. **Mechanism text a person edits.** `_parse_mechanism` is tested on `show` output and on
   surrounding whitespace (aa158dd); not on CRLF line endings, comments, reordered lines or
   a block for the other macro. Expected: the same mechanism, or an error that names the
   line.
5. **The suite on Julia 1.10.** CI runs Julia 1 only; on 1.10 the full suite peaks at
   3.30-4.36 GB, above CLAUDE.md's 4 GB budget in four of seven runs (ROADMAP Known issue,
   Task 63; PR question 9). Expected: the suite fits beside a language server on every
   supported version.

## File Structure

Net effect of the plan (`git diff 7205ba4 HEAD --stat` at the end: 43 files besides the plan):

- `src/types.jl` (modify): `Met`, the raw input (`RawForm`, `RawStep`, `RawMechanism`,
  `RawSite`, `RawAllostericMechanism`), the graph (`Counts`, `Form`, `Edge`, `Mechanism`,
  `AllostericMechanism`), `_canonicalize`, the representation checks, the Sig and its
  decoders, the lifts and compiled types, `EnzymeReaction` and its name checks.
- `src/names.jl` (create, Task 33): the naming chokepoint `name(p, m)`, the label renderers,
  and from Task 59 the DSL text of metabolites, forms and steps.
- `src/admission.jl` (create, Task 20): `_check_well_posed`, the boundary check the macros,
  the parser, the enumerator's exits and every derivation entry call.
- `src/reaction_graph.jl` (create, Task 29; delete, Task 34): `module ReactionGraph`, the
  graph built beside the Step types until the swap moves it into `types.jl` and `names.jl`.
- `src/dsl.jl` (modify): the three macros, the walker `_raw_mechanism`, `show` for every
  mechanism type, `_checked_mechanism`, `_parse_mechanism`.
- `src/thermo.jl`, `src/rate_eq_derivation.jl` (modify): the joint solve and the derivation
  over edges and masks.
- `src/mechanism_enumeration.jl` (modify): moves and seeds on the graph and raw input.
- `src/identify_rate_equation.jl` (modify): CSV columns `mechanism`, `mechanism_id`,
  `parent_id`.
- `src/EnzymeRates.jl`, `Project.toml` (modify): includes; `SHA` and `Unicode` deps; JET
  compat `0.9, 0.11, 0.12`.
- `test/test_reaction_graph.jl`, `test/fixtures/phase3_canonical_golden.txt` (create); FX
  gains each spec's raw value (`@_testhelper_with_raw`) and the reader `_testhelper_graph`;
  every other test file is migrated (S2b) and ported (S3c-S3i).
- `docs/src/**`, `ROADMAP.md` (modify, in the task that makes each text false);
  `.claude/CLAUDE.md` (Task 62, for Denis); `docs/superpowers/plans/2026-10-16-graph-mechanism.md`
  (Task 0).

## Replay procedure

Every task follows these steps; its own section lists only its units, files, red, green,
gate and departures.

A **commit unit** is one commit on `graph-mechanism`: a main prototype commit plus the
commits folded into it, listed in prototype order (`[main + fold + ...]`). For each unit:

1. Apply the test paths of each commit in the unit, in order:
   `git show --format= --binary <sha> -- <test/ paths of that commit> | git apply`
   (skip a commit with no test path). `git show --stat --format= <sha>` lists its paths.
2. Red: the focused run of the unit's test files. Expected: the task's red line. A unit
   that changes only tests says "green before and after" and runs once.
3. Apply the remaining paths (src, docs, ROADMAP.md, Project.toml, .claude) of each commit,
   in order: `git show --format= --binary <sha> -- <paths> | git apply`. A piece of 216be7f
   names its path; a ROADMAP entry operation is written out in the task.
4. Departures, if the task lists any, exactly as written there.
5. Green: the focused run again; expected counts as stated.
6. Commit: `git add` exactly the unit's paths;
   `git show -s --format=%B <main> | git interpret-trailers --if-exists addIfDifferent --trailer "Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE" > $P/replay/msg`,
   edit `$P/replay/msg` where the task says, then `git commit -F $P/replay/msg`.

At the task's end: the full suite (expected count in the gate), the census when the task
changed `src/`, and the task's other gates. A `git apply` that fails means a step was
skipped or reordered: every apply in this plan was checked, in this order, on a scratch
index built from 7205ba4 (`$S/p3skel/sim.py`, `spec.py`); stop and compare with the task's
commit list instead of editing a patch. The pieces the fix round moved (the `split` rows of
the commit map, and 4d2af5a and 8c867d3 into Tasks 38-39) were checked the same way, in
order, on the dry replay's commits (`$P/plan/fix_counts/`), with every checkpoint they
touch. The R62-R64 update's moves (a494c9c into Task 10, ad3ca91 into Task 57, da2d8c7
into Task 63) were checked the same way, with the S4 and S5 checkpoints
(`$P/plan/evidence/r62_update/`), and so were the R66 update's (eef3af7 into Task 57,
2af6a00 into Task 63; `$P/plan/evidence/r66_update/`) and the R67 update's (7c2a657 into
Task 57; `$P/plan/evidence/r67_update/`). Expected full-suite counts are the prototype's
count at the same point (its report's table) plus what the replay moved ahead: +2 from
Task 1 (68f46f0's testset) everywhere, +45 from Task 34 to Task 51 (050e16f's testset,
which the prototype adds at S3i), and the folded S5 pins in S4 (+16 from Task 57, +1
from Task 60).
Task 0 records the base count (38,174 expected, derived: the S1 report's 38,163 outside Aqua
and JET plus their 11; 1.10 ran 38,167 at the base and runs 7 fewer than 1.12 at the tip).

## Commit map

Every commit of `git log --reverse --format='%h %s' 4f9a34a..7c2a657`, in order: the stage
whose report lists it (`r1`/`r2` = that stage's review rounds; `S1 r2` = S1's minors, made at
S2a's start; `S2 r1` = the whole-S2 review; `S5m` = the S5 minors; `S5m2` = their second
round, R62; `S5m3` = their third, R66; `S5m4` = their fourth, R67), and the task that
replays it. "Task N" alone means the commit is a unit's main commit, or its own commit;
"folded into <sha>" means it is applied in that unit and squashed into its commit.

Checked by script: the table's 233 rows list every commit in the range exactly once and in
order. A row's Task column names one task, or starts `split` for a commit applied in pieces
and names each piece's task: 216be7f by file and ROADMAP entry across eleven tasks (S5 §9
item 3; its windows are the first item of Open for Denis), and 2b29fd2, 4d2af5a, 3962124,
b42cc90, 4952d79, 919bc4c and cd7a4ae by path, each piece in the task its report's section 9
or fold note names (Tasks 53 and 55 say where a piece cannot follow its note). The check
sets its own paths, so it runs in a fresh shell, before Task 0 too:

```bash
P=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3
W=/home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/graph_mechanism_wt
F=$P/plan/2026-10-16-graph-mechanism.md
rows() { awk -F'|' '/^\| [0-9]+ \| [0-9a-f]{7,} \|/{gsub(/ /, "", $3); print $3}' "$F"; }
test -r "$F" && [ "$(git -C $W rev-list --count 4f9a34a..7c2a657)" = 233 ] &&
    [ "$(rows | wc -l)" = 233 ] &&
    diff <(git -C $W log --reverse --format=%h 4f9a34a..7c2a657) <(rows) &&
    echo "map: all 233 commits, each once, in order" || echo "MAP CHECK FAILED"
test -r "$F" && awk -F'|' '/^\| [0-9]+ \| [0-9a-f]{7,} \|/ { r++
    if ($6 ~ /^ split/) {
        if ($6 !~ /Task [0-9]+.*Task [0-9]+/) { print "split into one task:", $3; bad++ }
    } else if (gsub(/Task [0-9]+/, "&", $6) != 1) { print "not one task:", $3; bad++ } }
    END { if (r == 233 && !bad) print "tasks: one per commit, two or more per split commit"
          else print "TASK CHECK FAILED" }' "$F" || echo "TASK CHECK FAILED"
```

Expected: `map: all 233 commits, each once, in order`, then `tasks: one per commit, two or
more per split commit`. Any other line, or a missing one, means the map is wrong.

| # | Commit | Subject | Stage | Task |
|---|---|---|---|---|
| 1 | a3de3c4 | Check the laws behind the string pins that phase 3 renames | S1 | Task 2 |
| 2 | 2bd4b5d | Orient tier-3 steps and tied-group leads by form index (C9) | S1 | Task 3 |
| 3 | 43d7af8 | Drop _solve_forward: the joint solve reads the stored direction | S1 | Task 3 |
| 4 | 8b7fb88 | Render phase 3's constant names: Xmod, bare modifiers, _res_plusA_minusP | S1 | Task 4 |
| 5 | a048534 | Check that a printed ping-pong law runs as written | S1 | Task 4, folded into 8b7fb88 |
| 6 | c7698fd | Pin the C9 fixtures' stored steps and fitted names exactly | S1 | Task 3, folded into 2bd4b5d |
| 7 | 7133720 | Check the kept-members fixture's law against the oracle | S1 | Task 2, folded into a3de3c4 |
| 8 | c86f61c | Check stored orientation on the Tier-3 and tied fixtures | S1 | Task 4 |
| 9 | 5e04242 | Name an internal-name clash as one, not as a rendered-name clash | S1 | Task 4, folded into 8b7fb88 |
| 10 | 1910833 | Pass _compute_alpha the concentration symbols its caller builds | S1 | Task 4, folded into 8b7fb88 |
| 11 | e74c8af | Reflow the tied-group docstring; use valid residuals in name examples | S1 | Task 4, folded into 8b7fb88 |
| 12 | e51c2dd | Point the golden test at the per-spec oracle instead of repeating it | S1 | Task 2, folded into a3de3c4 |
| 13 | 55e310d | Keep the invariants check in the dual-role testset | S1 | Task 4, folded into 8b7fb88 |
| 14 | adeb136 | Retire the residual-names ROADMAP entry; correct two texts it outdated | S1 | Task 4, folded into 8b7fb88 |
| 15 | cf1eb53 | Say what the golden test's oracle comment covers | S1 r2 | Task 2, folded into a3de3c4 |
| 16 | c378047 | Define the reactant-copy predicate once | S1 r2 | Task 4, folded into 8b7fb88 |
| 17 | df09b1d | Show a dead-end form's internal name as the code renders it | S1 r2 | Task 4, folded into 8b7fb88 |
| 18 | a8884ce | Keep test_dsl.jl within 92 columns | S1 r2 | Task 4, folded into 8b7fb88 |
| 19 | 8d5e465 | Add the raw input types the DSL walker returns | S2a | Task 5 |
| 20 | b3eabc4 | Read mechanism blocks into raw input with one walker | S2a | Task 6 |
| 21 | 7cf6b0f | Let oligomeric_state: take a tuple of multiplicities | S2a | Task 7 |
| 22 | 6075bba | Let a metabolite take the name M | S2a | Task 7 |
| 23 | 51fad19 | Reject reserved metabolite names and names that hold `_` | S2a | Task 8 |
| 24 | 410f7c0 | Reject a reaction whose metabolite names spell one name two ways | S2a | Task 9 |
| 25 | 9336851 | Name the raw input types in types.jl's header | S2a | Task 5, folded into 8d5e465 |
| 26 | 157ff8c | Test the constant-name clash a conformation label A or I causes | S2a r1 | Task 8, folded into 51fad19 |
| 27 | 6f9511f | Reject an empty list of catalytic multiplicities | S2a r1 | Task 7, folded into 7cf6b0f |
| 28 | 62280d3 | Pin the full walker messages for a missing and a forbidden step tag | S2a r1 | Task 6, folded into b3eabc4 |
| 29 | 10a974b | Say where a single step's tag goes without claiming groups need two steps | S2a r1 | Task 6, folded into b3eabc4 |
| 30 | fa2cc67 | Read a tag after a parenthesized right-hand side as the step's tag | S2a r1 | Task 6, folded into b3eabc4 |
| 31 | 03e2143 | Say that the metabolite-name guards cannot fire for constructor-built input | S2a r1 | Task 8, folded into 51fad19 |
| 32 | 62ab877 | Name the spelled-twice check in EnzymeReaction's docstring; balance the example | S2a r1 | Task 9, folded into 410f7c0 |
| 33 | c80e857 | Write the allosteric test fixtures inline, with tags on their steps | S2b | Task 10 |
| 34 | c2cd9b3 | Write the tests and docs in the modifier grammar | S2b | Task 11 |
| 35 | bb9b2e3 | Migrate the DSL error tests whose message or point moves with the grammar | S2b | Task 12 |
| 36 | b0b3a4e | Check the zero-inactive-numerator fixture's law against the oracle | S2b | Task 12 |
| 37 | 96846b8 | Rename the zero-inactive-numerator fixture's free form E_c to E | S2b | Task 12 |
| 38 | 9e6e64c | Migrate the eval-built DSL loops that interpolate role-tagged lines | S2b | Task 12 |
| 39 | 3c0c040 | Write the tests' titles and comments in the modifier grammar | S2b | Task 13 |
| 40 | 588e99f | Describe the modifier grammar in the mechanism macros' docstrings | S2b | Task 13 |
| 41 | fcfd5e6 | Quote the modifier grammar in the docs' prose | S2b | Task 13 |
| 42 | 9e41af7 | Name the raw-value test macro's argument instead of sketching a block | S2b | Task 13 |
| 43 | a494c9c | Restore the space after `p_mix =` in the mixed-tag regression test | S2b r1 | Task 10, folded into c80e857 |
| 44 | 9227cbe | Test that AllostericMechanism rejects a ligand on two regulatory sites | S2b r1 | Task 12 |
| 45 | bc804d6 | End the plain release in "::Mod is never a step tag" with a role tag | S2b r1 | Task 12 |
| 46 | 756abee | Restore the tests of three @enzyme_reaction forms the migration rewrote | S2b r1 | Task 13 |
| 47 | c34c525 | Name the raw-value test macro's argument as the macro does | S2b r1 | Task 13, folded into 9e41af7 |
| 48 | 3b28fca | Reject `catalytic_steps:`, `catalytic_inhibitors:` and plain `regulators:` | S2c | Task 14 |
| 49 | 68a4961 | Reject the role tag `::Inh` | S2c | Task 14 |
| 50 | fa09735 | Reject a tag on an `allosteric_regulators:` entry | S2c | Task 14 |
| 51 | 9399c58 | Reject a keyword multiplicity on a `regulatory_site` line | S2c | Task 14 |
| 52 | 822d590 | Reject a tag after a parenthesized single step | S2c | Task 14 |
| 53 | 9736179 | Read every mechanism block with the walker | S2c | Task 14 |
| 54 | b361684 | Delete the inhibitor-grammar parser | S2c | Task 14 |
| 55 | 961ec20 | Reject `dead_end_inhibitors:` in @enzyme_reaction | S2c | Task 15 |
| 56 | 1748df6 | Keep the multiplicity test within 92 columns | S2c | Task 14, folded into 9399c58 |
| 57 | 37c6a8b | Reject `allowed_catalytic_multiplicities:` in @enzyme_reaction | S2c | Task 15 |
| 58 | f8d0199 | Take one chemical formula per atom bracket in @enzyme_reaction | S2c | Task 15 |
| 59 | 59136b3 | Advise `X::Mod` in the Step role-change error | S2c | Task 16 |
| 60 | c457930 | Write a reactant's modifier copy as `X::Mod` in the form-name error | S2c | Task 16 |
| 61 | 8fc17da | Say why `_peel_tag!` reads `::Inh` as a role tag | S2c | Task 16 |
| 62 | 2493151 | List a mechanism's modifiers as `modifiers:` in its summary | S2c | Task 17 |
| 63 | 67dc293 | Print a reaction as the @enzyme_reaction block that builds it | S2c | Task 17 |
| 64 | 3c44dd3 | Document the macros' grammar, each with a checked example | S2c | Task 17 |
| 65 | 9838d79 | Sort a reaction's regulators by name, then competitive inhibitor first | S2 r1 | Task 17, folded into 67dc293 |
| 66 | 5cc7f8e | Name the reaction in one line in the no-catalytic-cycle error | S2 r1 | Task 17, folded into 67dc293 |
| 67 | a421d95 | List in a mechanism's summary the names its steps bind as modifiers | S2 r1 | Task 17, folded into 2493151 |
| 68 | 6c85d75 | Reject an empty oligomeric_state: tuple while the macro expands | S2 r1 | Task 15 |
| 69 | 257e81b | Name the step tags as the grammar writes them in the missing-tag error | S2 r1 | Task 16 |
| 70 | be61164 | Hint where a bracket or a role tag goes only for that value | S2 r1 | Task 16 |
| 71 | 7adea3e | Test that each macro's list labels accumulate over repeated lines | S2 r1 | Task 17 |
| 72 | b3431bf | Name the MWC dimer fixture's forms as its block writes them | S2 r1 | Task 17 |
| 73 | 60bd37a | Say "modifier copy" in test titles, comments and locals | S2 r1 | Task 17 |
| 74 | 4a22318 | Say that `_peel_tag!` reads `::Inh` the way it reads `::Mod` | S2 r1 | Task 16, folded into 8fc17da |
| 75 | 0a13357 | Say in Getting Started that @enzyme_reaction rejects names spelled two ways | S2 r1 | Task 9, folded into 410f7c0 |
| 76 | b588424 | Read the init golden's groups through steps(m) | S3a | Task 18 |
| 77 | db5b3ff | Compare mechanisms by their steps, not their reactions | S3a | Task 18 |
| 78 | 10527e7 | Name a mechanism's substrates, products and modifiers | S3a | Task 18 |
| 79 | 31ff003 | Render constant and form labels from the mechanism | S3a | Task 19 |
| 80 | 6ba4fed | Read metabolite names from the mechanism in the derivation and moves | S3a | Task 19 |
| 81 | 2b79bd3 | Check a mechanism's rate law at the DSL, the enumerator's exits and derivation | S3a | Task 20 |
| 82 | f7bd003 | Name constants without a cache inside the mechanism | S3a | Task 21 |
| 83 | 4d13802 | Hold a regulatory site's ligands as names | S3a | Task 21 |
| 84 | e4b5927 | Name the macro in the bottomless-segment testset's title | S3a | Task 20, folded into 2b79bd3 |
| 85 | ad4aca2 | Write out the test mechanisms picked as an enumeration's first | S3b | Task 22 |
| 86 | cd380a2 | Sample the oracle population and the fit census by content | S3b | Task 22 |
| 87 | 1362b5b | Order steps and kinetic groups by content | S3b | Task 23 |
| 88 | 9916684 | List a mechanism's forms in the canonical form order | S3b | Task 23 |
| 89 | 77b7b61 | Break a segment root's tie by the canonical form order | S3b | Task 24 |
| 90 | eb508a5 | Key the cycle basis's form rows by form, in the canonical order | S3b | Task 24 |
| 91 | 9c6d4cb | Identify forms by content in the reach, flip and dead-end sets | S3b | Task 24 |
| 92 | 1a33d0f | Pick a tied productive twin by the canonical form order | S3b | Task 25 |
| 93 | 7f1aa95 | Order catalytic topologies and their steps by form content | S3b | Task 25 |
| 94 | af01ff6 | Order a split's residual contexts by content | S3b | Task 25 |
| 95 | 0266794 | Sort the dead-end move's placements in the canonical form order | S3b | Task 25 |
| 96 | 74a8fc6 | Order each seed base's groups canonically before flipping them | S3b | Task 26 |
| 97 | 61a1256 | Record that identify's fit starting points follow processing order | S3b | Task 26 |
| 98 | a2bbf86 | Accept forms that share an internal name but render apart | S3b | Task 26 |
| 99 | 8107bc0 | Read test mechanisms through one adapter | S3c | Task 27 |
| 100 | 990cb6a | Read the test helpers' mechanisms through the adapter | S3c | Task 27 |
| 101 | 279855d | Take each test's reaction or names, not a mechanism's reaction | S3c | Task 27 |
| 102 | f57370e | Build enumeration fixtures with the macros alone | S3c | Task 27 |
| 103 | 0e051b4 | Read a mechanism's kinetic groups through steps(m) | S3c | Task 27 |
| 104 | 56f7871 | Pin the one-flip message for free forms D and E | S3c | Task 28 |
| 105 | bdd16e2 | Pin the children of a seed with substrate Iinh and inhibitor I | S3c | Task 28 |
| 106 | 1cd8cd5 | Rename the segment-root testset's copy local | S3c | Task 28 |
| 107 | d52040f | Read form labels where enumeration tests read internal names | S3c | Task 28 |
| 108 | 380a8c4 | Drop the ter-ter init timing figure from its test comment | S3c | Task 28 |
| 109 | e080b63 | Read isomerizations off the steps in the population predicate | S3c | Task 28 |
| 110 | e12daba | Add the reaction graph's value types, content order and identity | S3d1 | Task 29 |
| 111 | 2db9d74 | Canonicalize raw input into the reaction graph | S3d1 | Task 30 |
| 112 | 584d7bf | Build and check the reaction graph from raw input | S3d1 | Task 30 |
| 113 | 43c50fc | Reject a ligand written twice on one regulatory site | S3d1 | Task 30 |
| 114 | 6177781 | Encode the reaction graph as one integer tuple, its compiled type's parameter | S3d1 | Task 31 |
| 115 | f8eb0df | Write a reaction graph back as raw input | S3d1 | Task 31 |
| 116 | bf41ad0 | Check the reaction graph against today's canonical form | S3d1 | Task 32 |
| 117 | 6ca71bd | Gate the graph constructor's speed at 10 µs per call | S3d1 | Task 32 |
| 118 | 1bd88e8 | Lift Sigs of lengths no other test decodes in the compile-once test | S3d1 | Task 31, folded into 6177781 |
| 119 | 1bf2fee | Compare tier-3 and tied-group orientation with today's constructor | S3d1 | Task 32 |
| 120 | 03fad30 | Name the graph tests' comparison types instead of a point in time | S3d1 r1 | Task 32 |
| 121 | d1fb2df | Test that one reaction stored once each way in two groups is rejected | S3d1 r1 | Task 30, folded into 584d7bf |
| 122 | 94eaa5d | Test the allosteric graph's per-group tag read | S3d1 r1 | Task 32 |
| 123 | 5298b84 | Read a binding's entry and forms once for both binding checks | S3d1 r1 | Task 32 |
| 124 | f7c1587 | Pin the canonical text of the tier-3 and tied-group cases | S3d2 | Task 33 |
| 125 | 6214bb1 | Move the naming chokepoint and label renderers to names.jl | S3d2 | Task 33 |
| 126 | b0a4b29 | Write the pivot regression's mechanism as an @enzyme_mechanism block | S3d2 | Task 33 |
| 127 | f1aa056 | Store mechanisms as reaction graphs, with Step views for unported code | S3d2 | Task 34 |
| 128 | 92ba4f1 | Type the bottomless-segment check's mechanism and keep it with the boundary check | S3d2 | Task 35 |
| 129 | 86d21f6 | Check a step's shared metabolites without building sets | S3d2 | Task 35 |
| 130 | 1894d39 | Build a mechanism's Step views once per move and once per naming | S3d2 | Task 35 |
| 131 | 924353b | Pin the graph types' fields beside the Step type's | S3d2 | Task 34, folded into f1aa056 |
| 132 | 69964b0 | Wrap the canonical-order comment at 92 columns | S3d2 r1 | Task 34, folded into f1aa056 |
| 133 | 8567c69 | Say what the Sig round-trip fixture holds | S3d2 r1 | Task 34, folded into f1aa056 |
| 134 | 7c9dd07 | Describe the renderers' context and the modifiers label as they are | S3d2 r1 | Task 34, folded into f1aa056 |
| 135 | 1ecb534 | Refuse a written mechanism other than m in the textbook bridge | S3d2 r1 | Task 34, folded into f1aa056 |
| 136 | ddf94d9 | Drop the second blank line before the Sig lifts | S3d2 r1 | Task 34, folded into f1aa056 |
| 137 | e5ef9c5 | Read Step input over a reaction through one _raw_of method | S3d2 r1 | Task 34, folded into f1aa056 |
| 138 | 38c9c75 | Build the Step views once in show and in the allosteric boundary check | S3d2 r1 | Task 35, folded into 1894d39 |
| 139 | 85e0790 | Reject a mechanism without free E | S3e | Task 36 |
| 140 | b9c535e | Reject steps that do not form one connected graph | S3e | Task 37 |
| 141 | 450a713 | Reject a declared name that no step exchanges | S3e | Task 37 |
| 142 | a76ddf2 | Check metabolite names in the constructors, decodability at the macros | S3e | Task 38 |
| 143 | 2aa3cbd | Check conformation labels in the constructors | S3e | Task 38 |
| 144 | ad74058 | Reject an empty or repeated multiplicity list | S3e | Task 39 |
| 145 | e7698be | Reject a regulator type on a competitive inhibitor | S3e | Task 39 |
| 146 | 56dbf00 | Reject an atom element that no formula writes | S3e | Task 39 |
| 147 | 5b72e27 | Reject Bool multiplicities | S3e | Task 39 |
| 148 | c6d5a47 | Build test fixtures without reaction wrappers and drop an unread helper method | S3e | Task 40 |
| 149 | 9cbb9af | Name only kept steps' inactive constants in the test reader, and pin a Theorell–Chance step | S3e | Task 40 |
| 150 | 44661cc | Pin a step that shares a metabolite in unequal copies, and require each spec's raw value | S3e | Task 40 |
| 151 | 4d2af5a | Wrap the constructor checks' lines at 92 columns | S3e | split by path (S3e §9: the wraps fold into tasks 4-6): Task 38 src/dsl.jl and TT, folded into a76ddf2; Task 38 TRG, folded into 2aa3cbd; Task 39 TD, folded into ad74058 |
| 152 | dfc5f79 | Tighten the Sig-encoding budgets to a few times their baseline | S3e | Task 40 |
| 153 | 86d0c31 | Test the invariant helper's and the oracle's failure paths again | S3e r1 | Task 37, folded into 450a713 |
| 154 | ff5e59c | Pin that free E holds no residual | S3e r1 | Task 40 |
| 155 | 8c867d3 | Call a ligand a ligand when a conformation label takes its name | S3e r1 | Task 38, folded into 2aa3cbd |
| 156 | 9221bfc | Pin the allosteric macro's decodability check | S3e r1 | Task 38, folded into a76ddf2 |
| 157 | 3bb8d94 | Pin contextual keywords as metabolite names | S3e r1 | Task 38, folded into a76ddf2 |
| 158 | 88f896e | Say which Sig budget catches which per-seed compile | S3e r1 | Task 40, folded into dfc5f79 |
| 159 | f3d5e2c | Check a mechanism's rate law on its reaction graph | S3f | Task 41 |
| 160 | 19e3fe8 | Name constants and solve the joint system over edges | S3f | Task 42 |
| 161 | f6be872 | Derive the rate law from the reaction graph through edge masks | S3f | Task 42, folded into 19e3fe8 |
| 162 | 43519c4 | Decode a compiled mechanism once per generator and derive from the value | S3f | Task 43 |
| 163 | a75ef55 | Reject the conformation labels A and I | S3f | Task 43 |
| 164 | 7fe5c6d | Name the tied group's lead by the function that picks it | S3f | Task 43, folded into a75ef55 |
| 165 | d427228 | Read segment sides through the name index and cycle rows through the touched forms | S3g | Task 44 |
| 166 | ee46733 | Pin labels A and I at the plain macro and the allosteric constructor | S3g | Task 44 |
| 167 | 8752036 | Flip rapid-equilibrium groups to steady state on the reaction graph | S3g | Task 45 |
| 168 | 0f795c4 | Decide the copy rule and the dead-end move on the reaction graph | S3g | Task 46 |
| 169 | 5d1197d | Split kinetic groups on the reaction graph | S3g | Task 47 |
| 170 | 2c40783 | Tag, relax and check mechanisms on the reaction graph | S3g | Task 48 |
| 171 | 612a678 | Write the enumeration hotspot fixtures inline with their exact children | S3g | Task 49 |
| 172 | 7d69b90 | Pin the turn of a tied member a split frees, and every bi-bi child's raw round trip | S3g | Task 47 |
| 173 | 698fd6a | Release the graph tests' census population after its two testsets | S3g | Task 49 |
| 174 | 78d8376 | Wrap the regrouping oracle's loop header at 92 columns | S3g r1 | Task 47, folded into 5d1197d |
| 175 | 01e2dbd | Name the seed variants as the readers of the Step flip helpers and content keys | S3g r1 | Task 48, folded into 2c40783 |
| 176 | 704bac2 | Read the names a count holds through one helper | S3g r1 | Task 48, folded into 2c40783 |
| 177 | 9cece42 | Ask whether a group holds an edge through one helper | S3g r1 | Task 48, folded into 2c40783 |
| 178 | 801cd5c | Drop a dead-end copy, an unread argument name and a single-use helper | S3g r1 | Task 48, folded into 2c40783 |
| 179 | 5fc1f44 | Write out every expected child of the first-regulator and relaxation testsets | S3g r1 | Task 49, folded into 612a678 |
| 180 | 9dee115 | Build the seeds from raw input | S3h | Task 50 |
| 181 | 2395401 | Build the seed variants and judge degeneracy on the reaction graph | S3h | Task 51 |
| 182 | 050e16f | Store a mechanism's table, labels and forms at their length | S3i | Task 34, folded into f1aa056 |
| 183 | 6ada4b2 | Print the compiled types' summary from the reaction graph | S3i | Task 52 |
| 184 | 36384ef | Read the reaction graph in every test, and build test input as raw input | S3i | Task 53 |
| 185 | 2692bb7 | Delete the Step views, Step, Species and Residual | S3i | Task 54 |
| 186 | 14ef525 | Share the raw isomerization test, the seed bases and the bound-inhibitor read | S3i | Task 55 |
| 187 | e69cd27 | Pin the seed screen's message for a graph that holds an isomerization | S3i | Task 55, folded into 14ef525 |
| 188 | 3962124 | Describe three tests by what they check | S3i r1 | split by path (S3i round-1 note): Task 34 TRG (m6), folded into f1aa056; Task 53 FX (m7), folded into 36384ef; Task 54 TT (B1), folded into 2692bb7 |
| 189 | b42cc90 | Read isomerizations, Theorell–Chance steps and group flags through the helpers | S3i r1 | split by path (S3i round-1 note): Task 53 TD (m4) and TE (B2, m2), folded into 36384ef; Task 55 TT (m1, and m3 in the same file), folded into 14ef525 |
| 190 | 1182f84 | Find a graph's merged complexes in one function | S3i r1 | Task 55, folded into 14ef525 |
| 191 | b784fdc | Count an element listed twice in a formula once, with the summed count | S4 | Task 56 |
| 192 | a520f8e | Let @enzyme_reaction resolve the package's names itself | S4 | Task 57 |
| 193 | 3f6c43c | Show a site tag, not the user's, in the tagged-regulator error's example | S4 | Task 57 |
| 194 | 13f202c | Reject a metabolite name that Julia reads as another name | S4 | Task 58 |
| 195 | 3464d63 | Print a mechanism as the text of the macro that builds it | S4 | Task 59 |
| 196 | 0ec7d51 | Read a mechanism back from its text without evaluating it | S4 | Task 60 |
| 197 | b8a604f | Compile a compiled mechanism's show once, whatever the Sig | S4 | Task 59, folded into 3464d63 |
| 198 | e59fb5b | Store each mechanism in identify's CSVs as its text and a stable id | S4 | Task 61 |
| 199 | aa158dd | Read a mechanism's text past the whitespace around it | S5 | Task 60, folded into 0ec7d51 |
| 200 | 9578d2a | Name the bare metabolite in a list-label hint under nested tags | S5 | Task 57, folded into 3f6c43c |
| 201 | 8aed057 | Claim the text read-back only for mechanisms the macros accept | S5 | Task 60, folded into 0ec7d51 |
| 202 | 02214d5 | Write the Unicode names the rejection tests pin as escapes | S5 | Task 58, folded into 13f202c |
| 203 | cd7a4ae | Share one reactant-copy test and one edge text | S5 | split by path (S5 §9 item 5): Task 59 src/names.jl, FX, TE and TO, folded into 3464d63; Task 60 TD and TRG, which edit 0ec7d51's tests, folded into 0ec7d51 |
| 204 | 4952d79 | Keep two printing tests' reasons and name the reach oracle's argument | S5 | split by path (S5 §9 items 5-6): Task 53 TRG (R58 m1), folded into 36384ef; Task 59 TT (R59 m5), folded into 3464d63 |
| 205 | 919bc4c | Name two multi-line test closures and the CSV testset by what they do | S5 | split by path (S5 §9 item 5): Task 59 TT (m7), folded into 3464d63; Task 61 TID (m7, m8), folded into e59fb5b |
| 206 | a8f6f51 | Build the allosteric uni-uni cases with one helper and pin each message | S5 | Task 53, folded into 36384ef |
| 207 | 2b29fd2 | Write the docs' constant names and form order as the package renders them | S5 | split by path: Task 3 re_vs_ss.md (the tier-3 tie-break C9 changes), folded into 2bd4b5d; Task 4 ping_pong.md and regulators.md (the names), folded into 8b7fb88 |
| 208 | be65542 | Say in the docs where reserved and doubly spelled names are rejected | S5 | Task 9, folded into 410f7c0 |
| 209 | 216be7f | Describe the reaction-graph mechanism in the docs and ROADMAP | S5 | split by file and ROADMAP entry (S5 §9): Task 20 mwc_allostery.md; Task 24 cha_king_altman.md; Task 34 getting_started.md and the ROADMAP Idea "Encode a mechanism's Sig ..." and Decided-against "Have the `Mechanism` constructors check ..." deleted; Task 36 ROADMAP "Take `:OnlyA` steps ..." (merged with P4-1 and dc4035b); Task 37 src/dsl.jl; Task 42 ROADMAP "`eq_complexity_filter` ..."; Task 48 ROADMAP "The `:OnlyA` guard misjudges ..."; Task 51 enumeration_engine.md; Task 54 re_vs_ss.md and ROADMAP "Enumerate Cleland iso ..."; Task 55 developer.md; Task 60 ROADMAP "`BigInt` arithmetic ..." |
| 210 | 81f9f1e | Describe the printed mechanism text and the CSV's text and id columns | S5 | Task 61, folded into e59fb5b |
| 211 | ea6bf86 | Correct three statements of fact in CLAUDE.md for Denis to take or drop | S5 | Task 62 |
| 212 | fbdb9aa | Call the test reader by its name in two test comments | S5 | Task 54, folded into 2692bb7 |
| 213 | a84d194 | Name the reaction label the allosteric promotion reads its multiplicities from | S5 | Task 15, folded into 37c6a8b |
| 214 | 5c08ede | Allow JET 0.12, the only JET release that runs on Julia 1.13 | S5 | Task 1 |
| 215 | 4a96d9a | Expect m_OnlyA_prod's true peak turnover, which a product can raise | S5 r1 | Task 1 |
| 216 | 26f7323 | Name the parser among the boundary check's callers | S5 r1 | Task 60, folded into 0ec7d51 |
| 217 | 5cfd1e6 | Say which slice of the cv table the identify tutorial shows | S5 r1 | Task 61, folded into e59fb5b |
| 218 | 32ea3c9 | Say the one-flip and :OnlyA conditions bind allosteric mechanisms only | S5 r2 | Task 60, folded into 0ec7d51 |
| 219 | 1b228ef | Name both id shapes CSV.jl reads as a number | S5 r2 | Task 61, folded into e59fb5b |
| 220 | 0b5e637 | Name m_OnlyA_prod's peak-turnover formula as a test helper | S5 r2 | Task 1, folded into 4a96d9a |
| 221 | 68f46f0 | Check m_OnlyA_prod's peak at a nonzero product on every Julia version | S5 r2 | Task 1, folded into 4a96d9a |
| 222 | 33b94d8 | Say what the analytical kcat function receives | S5m | Task 1, folded into 4a96d9a |
| 223 | dc4035b | Say the ROADMAP's one-flip condition binds allosteric mechanisms only | S5m | Task 36, folded into the ROADMAP entry text that task writes |
| 224 | 09aec13 | Name both id shapes CSV.jl reads as a number in the CSV test | S5m | Task 61, folded into e59fb5b |
| 225 | 0b2021d | Name the reserved metabolite names in the identify tutorial | S5m | Task 9, folded into 410f7c0 |
| 226 | 761337e | Name only a bare metabolite in the tag hints, past a multiplicity call | S5m | Task 57, folded into 3f6c43c |
| 227 | cfe9211 | Draw the kcat check's parameters only for a spec with a kcat | S5m | Task 1, folded into 4a96d9a |
| 228 | 7395164 | Record the full suite's peak above 4 GB on Julia 1.10 | S5m | Task 63 |
| 229 | ad3ca91 | Unwrap only a named call in the tag hints, never an operator call | S5m2 | Task 57, folded into 3f6c43c |
| 230 | da2d8c7 | Keep the 1.10 Known issue to what the traces measured | S5m2 | Task 63, folded into 7395164 |
| 231 | eef3af7 | Unwrap only a multiplicity call in the tag hints | S5m3 | Task 57, folded into 3f6c43c |
| 232 | 2af6a00 | Scope the 1.10 Known issue's garbage to what the suite traces measured | S5m3 | Task 63, folded into 7395164 |
| 233 | 7c2a657 | Unwrap a multiplicity call of several Integers in the tag hints | S5m4 | Task 57, folded into 3f6c43c |

## Stage checkpoints

After each stage's last task, from `$W`:

```bash
git diff --stat <tip> HEAD -- src test docs ROADMAP.md Project.toml .claude
```

must list exactly the files below, and `git diff <tip> HEAD -- <file>` of each must hold only
the changes named for it. The file lists come from a replay of every task's git steps on a
scratch index (`$S/p3skel/sim.py`, log `$S/p3skel/sim_final.log`), which applied every patch
in this order without a 3-way merge.

Standing departures, in every checkpoint up to S4:

- **D1** `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`: 7205ba4's edits
  (`git diff 4f9a34a 7205ba4 -- docs/superpowers/specs`).
- **D2** `docs/superpowers/plans/2026-10-16-graph-mechanism.md`: this plan (Task 0).
- **D3** `ROADMAP.md`: 7205ba4's P4-1 edits (`git diff 4f9a34a 7205ba4 -- ROADMAP.md`), plus
  the entries the stage list says were written ahead.
- **D4** Task 1's test environment: `Project.toml` (5c08ede's JET compat), FX (4a96d9a,
  0b5e637, 33b94d8) and TRD (4a96d9a, 68f46f0, cfe9211).
- **Docs ahead:** each S5 commit or piece, and 3962124's TRG piece, applied before its
  prototype position, cumulatively, as listed.

| Stage | After | Tip | Files besides D1-D4 | Changes ahead of the tip |
|---|---|---|---|---|
| S1 | Task 4 | a8884ce | docs/src/deriving/{ping_pong,re_vs_ss,regulators}.md | 2b29fd2 |
| S2a | Task 9 | 62ab877 | + docs/src/developer.md, docs/src/getting_started.md, docs/src/identify/tutorial.md | + 0a13357, be65542, 0b2021d |
| S2b | Task 13 | c34c525 | as S2a | as S2a |
| S2c | Task 17 | 0a13357 | getting_started.md leaves (0a13357 is the tip); + docs/src/identify/enumeration_engine.md | 2b29fd2, be65542, 0b2021d, a84d194 |
| S3a | Task 21 | e4b5927 | + docs/src/deriving/mwc_allostery.md | + 216be7f (mwc_allostery.md) |
| S3b | Task 26 | a2bbf86 | + docs/src/deriving/cha_king_altman.md | + 216be7f (cha_king_altman.md) |
| S3c | Task 28 | e080b63 | as S3b | as S3b |
| S3d1 | Task 32 | 5298b84 | as S3b | as S3b |
| S3d2 | Task 35 | 38c9c75 | + src/types.jl, TRG, docs/src/getting_started.md | + 050e16f (the copies and TRG's testset); 3962124 (TRG: that testset's helper names its argument `m`); 216be7f (getting_started.md); ROADMAP: the Sig Idea and the reaction-check Decided-against entry deleted |
| S3e | Task 40 | 88f896e | + src/dsl.jl | + 216be7f (the `@enzyme_mechanism` docstring); ROADMAP "Take `:OnlyA` steps ..." (Task 36's text, with dc4035b) |
| S3f | Task 43 | 7fe5c6d | as S3e | + ROADMAP "`eq_complexity_filter` ..." |
| S3g | Task 49 | 5fc1f44 | + TE | + R55 m2 (TE lines 6697-6723 re-indented); ROADMAP "The `:OnlyA` guard misjudges ..." |
| S3h | Task 51 | 2395401 | as S3g | + 216be7f (enumeration_engine.md) |
| S3i | Task 55 | 1182f84 | src/types.jl leaves (050e16f and 3962124 are in the tip); TRG stays, 4 lines (4952d79); + TT; re_vs_ss.md and developer.md also carry 216be7f's pieces | 050e16f and 3962124 leave the list; + a8f6f51 (TT), fbdb9aa (TE, TRD), 4952d79 (TRG); 216be7f (re_vs_ss.md, developer.md); ROADMAP "Enumerate Cleland iso ..." |
| S4 | Task 61 | e59fb5b | + src/admission.jl, src/names.jl, TD, TID, TO, TRG (back: 02214d5, cd7a4ae, 4952d79), docs/src/identify/model_selection.md | + aa158dd, 9578d2a, 8aed057, 02214d5, cd7a4ae, 4952d79, 919bc4c, 26f7323, 32ea3c9, 5cfd1e6, 1b228ef, 09aec13, 761337e, ad3ca91, eef3af7, 7c2a657, 81f9f1e; ROADMAP "`BigInt` arithmetic ..." |
| S5 | Task 63 | 7c2a657 | only ROADMAP.md, TE, D1, D2 | none; ROADMAP.md differs by exactly P4-1's three hunks (Task 36's merged paragraph, the new entry "Let an allosteric ping-pong stop both half-reactions ...", the iso entry's "Revisit allosteric ping-pong ..." sentence); TE by R55 m2's 27 re-indented lines |

S4's full file list (24 files): Project.toml, ROADMAP.md, docs/src/deriving/{cha_king_altman,
mwc_allostery, ping_pong, re_vs_ss, regulators}.md, docs/src/{developer,getting_started}.md,
docs/src/identify/{enumeration_engine,model_selection,tutorial}.md, the spec, src/admission.jl,
src/dsl.jl, src/names.jl, FX, TD, TID, TE, TO, TRD, TRG, TT (plus the plan).

## Tasks

Each task lists its commit units (`[main + folds]`, one commit each, applied by the replay
procedure), its files, interfaces, red, green and gate, the S5 docs edits it takes and its
departures. Counts are Julia 1.12 full-suite passes unless marked focused. "Census = base6"
means `compare.jl census/base6 census/replay_T<NN>` (with the map from Task 4 on, without it
before); "census = S3b_end" means `compare.jl census/S3b_end census/replay_T<NN>
census/rename_phase3.toml`; both expect `RESULT PASS` with zero differences unless the task
says otherwise.

The tasks fall into four groups. The Runners section defines every shorthand and count
convention the steps use; give an implementer the Global Constraints, the Runners and
the Replay procedure with its task.

Section 9's analysis probes are not rerun: the S1 flip scan, S3a's call-site mutation
check, S3b's K10 scans, S3h's exit-check search (D-S3h-9), S2c's grep for removed
tokens and S4's 115,209-mechanism round-trip probe. Each stage checkpoint shows the
branch's tree equals the prototype's stage tip, on which they ran, but for the listed
departures.

## Tasks 0-13: setup, test environment, S1 names and C9, S2a walker and names, S2b migration

### Task 0: Commit the plan and prepare the replay [setup]

**Replays:** nothing.

**Files:**
- Create: `docs/superpowers/plans/2026-10-16-graph-mechanism.md` (a copy of
  `$P/plan/2026-10-16-graph-mechanism.md`)
- Create (scratch, not committed): `$R/phase3/replay/env.sh`, `$P/replay/mutate.py`, `$L/`,
  `$P/replay/trees/`, `$P/census/env_replay/`, `$P/census/replay_T00/`; `$W/Manifest.toml`
  (git-ignored, a copy of the prototype's)

**Interfaces:**
- Produces: `env.sh`, which every later shell block of this plan sources first. It
  sets the Paths' `R`, `S`, `P`, `W`, `L`, leaves the shell in `$W`, and defines:
  `focused <log> <test file>...` (the focused run, Julia 1.12), `full <log>` (the full
  suite, Julia 1.12), `wait_log <log>` (waits up to 580 s, prints the tail),
  `summary <log>`, `pristine <log>`, `masked <log> <reference log>` (two suite logs
  compared with every number masked), `census_run <NN> <ref> [map]`, `census_check <NN>`
  and `chg` (the changed lines of a diff). A `<log>` without a `/` means `$L/<log>.log`.
- Produces: `$P/replay/mutate.py <name> [--revert]`, which applies one named source mutant
  for a red that only a mutant shows, and restores the file exactly (Step 8). Tasks 32, 37,
  40, 46 and 47 call it.
- Produces: the base count (38,174) and `census/replay_T00`.
- Consumes: `$P/jl.sh`, `census/{mkenv,corpus,compare}.jl`, `census/base6`, `census/env`.

**Gate:** R62-R67 recorded (Step 1); full suite 38,174 pass, pristine; census of
the base against `census/base6` without the map: `RESULT: PASS`, zero differences;
`mutate.py`'s sum (Step 8); JET 0.12 at every stage tip (Step 9).

- [ ] **Step 1: Check the worktree and Denis's rulings**

Rulings R62-R67 in `$P/RULINGS.md` answer what this plan asked Denis (R63 the questions
before Task 0, R64 those before Task 63), set how SDD runs it (R65) and rule on the
prototype's last five commits (R62, R66, R67); this plan carries each answer. Check that
they are recorded before anything runs.

```bash
cd /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/graph_mechanism_wt
git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD; git status --short
grep -c '^R6[2-7]\.' /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/RULINGS.md
```

Expected: `graph-mechanism`, then `7205ba4`, then nothing, then `6` (R62 to R67).
Anything else: stop and ask.

- [ ] **Step 2: Write `env.sh` and make the log directories**

```bash
mkdir -p /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/logs \
         /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/trees
cat > /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh <<'EOF'
# ABOUTME: Paths and run helpers for replaying the phase-3 prototype onto graph-mechanism.
# ABOUTME: Source it at the start of every shell; it leaves the shell in the worktree $W.
R=/home/denis.linux/claude-scratchpad/EnzymeRates-12f42995
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates
S=$S/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
P=$S/phase3
W=$R/graph_mechanism_wt
L=$P/replay/logs

# The focused run's Julia code: the fixtures file, then each test file named on the
# command line in a testset of its own, as runtests.jl includes them, so a failing file
# does not stop the files after it; the outer testset prints each file's count.
FOCUSED_JL='using TestEnv; TestEnv.activate()
using Test, EnzymeRates, LinearAlgebra, Random
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
@testset verbose = true "focused" begin
    for f in ARGS
        @testset "$f" begin include(f) end
    end
end'

# A log argument without a `/` names $L/<log>.log.
_log() { [[ $1 == */* ]] && echo $1 || echo $L/$1.log; }

# Usage: focused <log> <test file>...   (Julia 1.12, in $W)
focused() {
    local log=$(_log $1)
    shift
    (cd $W && $P/jl.sh $log 1800 --project=. -e "$FOCUSED_JL" "$@")
}

# Usage: full <log>   (the full suite on Julia 1.12, in $W)
full() {
    (cd $W && $P/jl.sh $(_log $1) 3600 --project=. \
        -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])')
}

# Usage: wait_log <log>   (waits up to 580 s for the run to end, then prints its tail)
wait_log() {
    local f=$(_log $1)
    timeout 580 bash -c "until grep -qs '^exit ' $f; do sleep 20; done"
    tail -40 $f
}

# Usage: summary <log>   (the test summary table and the run's exit line)
summary() { sed -n '/^Test Summary/,/^exit /p' $(_log $1); }

# Usage: pristine <log>   (a full-suite log: the number of `[ Info:` lines between
# "Testing Running tests" and the summary, then every other line there but precompilation)
pristine() {
    local body
    body=$(sed -n '/Testing Running tests/,/^Test Summary/p' $(_log $1))
    grep -c '^\[ Info:' <<< "$body"
    grep -v '^\[ Info:\|Testing Running tests\|^Test Summary\|^Precompiling' <<< "$body" |
        grep -v '✓ EnzymeRates\|dependency successfully precompiled'
}

# Usage: masked <log> <reference log>   (the two suite logs from "Testing Running tests" on,
# every number masked as N by s3h/scripts/mask.sh, precompilation lines dropped and the
# warm/cold ratio's number format ignored; prints nothing when they match line for line)
_masked() {
    bash $P/s3h/scripts/mask.sh $1 | grep -v 'Precompiling\|✓ EnzymeRates\|precompiled in' |
        sed 's/Ne-N/N/g'
}
masked() { diff <(_masked $(_log $1)) <(_masked $2); }

# Usage: census_run <NN> <reference> [census/rename_phase3.toml]
# Builds census/replay_T<NN> from $W, then compares census/<reference> with it, renaming
# the reference's names by the map when a third argument gives it.
census_run() {
    rm -rf $P/census/replay_T$1 $P/census/logs/replay_T$1.log \
           $P/census/logs/cmp_replay_T$1.log
    (cd $P && ./jl.sh census/logs/replay_T$1.log 1800 --project=census/env_replay \
         --heap-size-hint=2500M census/corpus.jl $W census/replay_T$1 \
         census/rename_phase3.toml &&
     ./jl.sh census/logs/cmp_replay_T$1.log 600 --project=census/env \
         --heap-size-hint=2500M census/compare.jl census/$2 census/replay_T$1 $3)
}

# Usage: census_check <NN>   (prints "RESULT: PASS", then every difference: none expected)
census_check() {
    local c=$P/census/logs/cmp_replay_T$1.log
    grep -x "RESULT: PASS" $c
    grep -E "order changed|MULTIPLICITY|-only [1-9]|differs [1-9]" $c
    grep -E "[1-9][0-9]* only in [AB]|inexact [1-9]|NOTE" $c
    grep "rel err" $c | grep -v "rel err 0.0 (tol"
    grep "directions:" $c | grep -v "directions: 0 of"
}

# Usage: <a diff> | chg   (its changed lines, without the file headers)
chg() { grep '^[-+]' | grep -v '^+++ \|^--- '; }

cd $W
EOF
```

Expected: no output.

- [ ] **Step 3: Check the runner, TestEnv and the three Julia versions**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh $L/T00_v112.log 300 -e 'using TestEnv; println(VERSION)'
$P/jl.sh $L/T00_v113.log 300 +1.13 -e 'println(VERSION)'
$P/jl.sh $L/T00_v110.log 300 +1.10 -e 'println(VERSION)'
grep -h "^1\.\|^exit" $L/T00_v112.log $L/T00_v113.log $L/T00_v110.log
```

Expected: `1.12.5`, `exit 0 ...`, `1.13.1`, `exit 0 ...`, `1.10.12`, `exit 0 ...`.

- [ ] **Step 4: Commit the plan**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
test -f $P/plan/2026-10-16-graph-mechanism.md
cp $P/plan/2026-10-16-graph-mechanism.md docs/superpowers/plans/
git add docs/superpowers/plans/2026-10-16-graph-mechanism.md
git commit -F - <<'EOF'
Add the phase-3 implementation plan

The plan replays the reviewed phase-3 prototype onto graph-mechanism task by
task: each task applies its prototype commits, test paths first, so the red
run comes before the source and the green run after it, and each stage ends at
a checkpoint against the prototype's stage tip.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: `git status --short` prints nothing; the log line ends `Add the phase-3
implementation plan`.

- [ ] **Step 5: Instantiate the worktree and build the census environment**

The worktree takes the prototype's 1.12 Manifest (a read-only copy from `$P/proto3_wt`,
resolved 2026-10-16 17:54, which the prototype's last 1.12 suites ran on). The General
registry has moved since (20:26), so a fresh resolve could pick other dependency versions
than the ones behind the reports' counts and pristine output.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
cp $P/proto3_wt/Manifest.toml $W/
$P/jl.sh $L/T00_instantiate.log 1800 --project=. -e 'using Pkg; Pkg.instantiate()'
(cd $P && ./jl.sh census/logs/mkenv_replay.log 1800 --heap-size-hint=2500M \
    census/mkenv.jl $W census/env_replay)
tail -2 $L/T00_instantiate.log; tail -3 $P/census/logs/mkenv_replay.log
git status --short
```

Expected: both logs end `exit 0` (the instantiate log may warn that the project changed
since the manifest was resolved: the manifest is the prototype tip's); the census log's
tail holds `ok 0.10.0 from .../graph_mechanism_wt/src/EnzymeRates.jl` (then `maxrss_kb`
and `exit 0`); `git status --short` prints nothing (Manifest.toml is ignored).

- [ ] **Step 6: Full suite at the base**

Background: `. .../replay/env.sh; full T00_full`. Then, until the log ends:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
wait_log T00_full
```

Then:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
summary T00_full; pristine T00_full; grep -o "JET v[0-9.]*" $L/T00_full.log | head -1
```

Expected: `EnzymeRates.jl | 38174  38174`, `Testing EnzymeRates tests passed`, `exit 0`;
`pristine` prints `12` and nothing else; `JET v0.11.x`. Every later count in this plan is
derived from 38,174: if the count differs, stop and report it.

- [ ] **Step 7: Census of the base**

Background: `. .../replay/env.sh; census_run 00 base6`. Then
`wait_log $P/census/logs/cmp_replay_T00.log` until it ends, and:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
census_check 00
```

Expected: `RESULT: PASS` and nothing else.

- [ ] **Step 8: Write the mutant helper**

Some tests pass on arrival: the prototype showed them failing only on a mutant of the code
they pin (a report's section 9 Start line names it). Tasks 32, 37, 40, 46 and 47 rerun those
reds with this helper. Each mutant is the prototype's (`$P/s3d/scripts/mutants.py`,
`$P/s3d1r1/scripts/mutants.py`, and the diffs of the snapshots `$P/s3e/snap/{mut_occupied,
mut_sig2, r1_red_mutants}` and `$P/s3g/snap/{mut_k3, mut_class, mut_tied}` against their
commits); each text occurs once in its file at the step that uses it, checked on the dry
replay's trees. It edits the file under the current directory, the worktree `env.sh`
leaves the shell in, and keeps the original beside itself in `$P/replay`; `--revert`
writes it back only if the file holds the mutant and nothing else changed.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
cat > $P/replay/mutate.py <<'EOF'
# ABOUTME: Applies or reverts one named mutant of the worktree's code, for the red runs of
# ABOUTME: tests that pass on arrival. Usage: mutate.py <name> [--revert].
import os
import sys

W = os.getcwd()  # the worktree: every step runs this from the shell env.sh leaves there
R = os.path.dirname(os.path.abspath(__file__))  # $P/replay, which keeps the saved file
FX = "test/mechanism_definitions_for_test_enzyme_derivation.jl"
# name: (path, text, mutant text); each text occurs once in its file when its task runs it
MUTANTS = {
    "tier3": ("src/reaction_graph.jl",
              "fwd[k], tier[k] = e.from < e.to, 3",
              "fwd[k], tier[k] = e.from > e.to, 3"),
    "no_alignment": ("src/reaction_graph.jl",
                     "    for l in unique(labels)\n        tied =",
                     "    for l in Int[]\n        tied ="),
    "first_tag": ("src/reaction_graph.jl",
                  "g::Int) = am.cat_allo_states[g]",
                  "g::Int) = am.cat_allo_states[1]"),
    "invariants_blind": (FX,
                         "function _testhelper_assert_mechanism_invariants(m)\n"
                         "    g = _testhelper_graph(m)",
                         "function _testhelper_assert_mechanism_invariants(m)\n"
                         "    return nothing\n    g = _testhelper_graph(m)"),
    "oracle_blind": (FX,
                     "        length(reached) == n ||\n"
                     "            error(\"_testhelper_steady_state_rate: the network",
                     "        true ||\n"
                     "            error(\"_testhelper_steady_state_rate: the network"),
    "occupied": ("src/types.jl",
                 "shared = _occupied(e.consumed) & _occupied(e.released)",
                 "shared = e.consumed & e.released"),
    "sig_per_mechanism": ("src/types.jl",
                          "_sig_of(m::Mechanism) = Tuple(_sig_mechanism!(Int[], m))\n",
                          "function _sig_of(m::Mechanism)\n"
                          "    sig = Tuple(_sig_mechanism!(Int[], m))\n"
                          "    metabolites(EnzymeMechanism{sig}())\n"
                          "    sig\n"
                          "end\n"),
    "residual_blind": ("src/types.jl",
                       "&& m.forms[1] == Form(E, 0, 0, 0) ||",
                       "&& m.forms[1].conformation == E &&\n"
                       "        m.forms[1].bound == Counts(0) ||"),
    "k3": ("src/mechanism_enumeration.jl",
           "    twin = _productive_twin(cm, trues(length(cm.edges)))",
           "    twin = _productive_twin(m, trues(length(m.edges)))"),
    "one_index": ("src/mechanism_enumeration.jl",
                  "_one_state_class(h::Int, t::Int) = _GaugeClass(h, t, 1)",
                  "_one_state_class(h::Int, t::Int) = _GaugeClass(h + t, 0, 1)"),
    "stale_direction": ("src/mechanism_enumeration.jl",
                        "    m isa Mechanism ? Mechanism(_raw(m, re, group)) :",
                        "    m isa Mechanism ? Mechanism(m.substrates, m.products, "
                        "m.table, m.conformations, m.forms,\n"
                        "        [Edge(e.from, e.to, e.consumed, e.released, re[k]) "
                        "for (k, e) in enumerate(m.edges)],\n"
                        "        _number_groups(group, maximum(group))[1]) :"),
}

name, revert = sys.argv[1], sys.argv[2:] == ["--revert"]
path, old, new = MUTANTS[name]
f, saved = f"{W}/{path}", f"{R}/mutant_{name}.orig"
text = open(f, encoding="utf-8").read()
if revert:
    orig = open(saved, encoding="utf-8").read()
    if text.count(new) != 1 or text.replace(new, old) != orig:
        sys.exit(f"{name}: {path} changed while mutated; its text is in {saved}")
    open(f, "w", encoding="utf-8").write(orig)
    os.remove(saved)
    print(f"restored {path}")
else:
    if text.count(old) != 1:
        sys.exit(f"{name}: {text.count(old)} matches in {path}, expected 1")
    open(saved, "w", encoding="utf-8").write(text)
    open(f, "w", encoding="utf-8").write(text.replace(old, new))
    print(f"mutated {path} ({name})")
EOF
md5sum $P/replay/mutate.py
```

Expected: `22217bbb69e412c83b11b2a9451011f4`. A different sum means the heredoc was
mistyped: write the file again.

- [ ] **Step 9: JET 0.12 at the prototype's stage tips**

From Task 1 on, the 1.12 environment resolves JET 0.12, which the prototype ran only at its
tip, so a JET report at an intermediate task would stop it with no prototype code to fall
back on. This pre-flight runs `test/test_aqua_jet.jl` under JET 0.12 on an archive of each
stage tip the checkpoints compare against (5c08ede's compat applied, the Manifest of Step
5), and JET with `ReactionGraph` targeted at S3d1's tip. The fix round ran it on
2026-10-17 with this Manifest (`$P/plan/fix_counts/jet/`): JET 0.12.3, Aqua and JET 11
of 11 at all 15 tips, and 0 reports for both target sets at 5298b84. The result depends on
the code and the Manifest alone, so read those logs unless the Manifest differs:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
cmp $W/Manifest.toml $P/proto3_wt/Manifest.toml && echo same-manifest
J=$P/plan/fix_counts/jet
for t in a8884ce 62ab877 c34c525 0a13357 e4b5927 a2bbf86 e080b63 5298b84 38c9c75 \
         88f896e 7fe5c6d 5fc1f44 2395401 1182f84 e59fb5b; do
    echo "$t $(grep -h 'JET version' $J/jet_$t.log)" \
         "$(grep -h 'test_aqua_jet.jl |' $J/jet_$t.log | tr -s ' ')"
done
grep -h '^targets' $J/jet_rg_5298b84.log
```

Expected: `same-manifest`; 15 lines `<tip> JET version 0.12.3 test/test_aqua_jet.jl | 11 11
<time>`; `targets (EnzymeRates, EnzymeRates.ReactionGraph): 0 reports` and `targets
(EnzymeRates,): 0 reports`. If the Manifest differs (Step 5 copies it, so it should not),
rerun the pre-flight with this Manifest, in the background (about 10 minutes once the
dependencies are compiled), and expect the same lines from `$L/T00_jet_<tip>.log`:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
JET_JL='using TestEnv; TestEnv.activate()
using Test, EnzymeRates, LinearAlgebra, Random
import JET; println("JET version ", pkgversion(JET))
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
@testset verbose = true "focused" begin
    @testset "test/test_aqua_jet.jl" begin include("test/test_aqua_jet.jl") end
end'
git show --format= 5c08ede -- Project.toml > $P/replay/jet_compat.patch
for t in a8884ce 62ab877 c34c525 0a13357 e4b5927 a2bbf86 e080b63 5298b84 38c9c75 \
         88f896e 7fe5c6d 5fc1f44 2395401 1182f84 e59fb5b; do
    D=$P/replay/trees/jet_$t
    rm -rf $D && mkdir -p $D && git archive $t | tar -x -C $D
    (cd $D && patch -s -p1 < $P/replay/jet_compat.patch && cp $W/Manifest.toml . &&
     $P/jl.sh $L/T00_jet_$t.log 1800 --project=. -e "$JET_JL")
    [ $t = 5298b84 ] && (cd $D && $P/jl.sh $L/T00_jet_rg_$t.log 1800 --project=. \
        $P/s3d/scripts/jet.jl)
    rm -rf $D
done
```

A JET report at any tip: stop and ask Denis (pin JET 0.11 on Julia 1.12 until Task 63, or
look at the report at that stage first).

---

### Task 1: Julia 1.13 and 1.10: JET 0.12 and m_OnlyA_prod's peak turnover [S5 test environment]

**Replays:** `[5c08ede]`, `[4a96d9a + 0b5e637 + 68f46f0 + 33b94d8 + cfe9211]` (S5 §9
item 7: first in the plan, since no phase-3 task makes it true or false).

**Files:**
- Modify: `Project.toml` (JET compat `0.9, 0.11, 0.12`)
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (`m_OnlyA_prod`'s peak
  formula and the `analytical_kcat_fn` field comment)
- Modify: `test/test_rate_eq_derivation.jl` (`_testhelper_test_analytical_kcat`, the
  fixed-constant testset)

**Interfaces:**
- Produces: `_testhelper_m_OnlyA_prod_peak_turnover(p)` (FX, inside `m_OnlyA_prod`'s `let`);
  `_testhelper_test_analytical_kcat(spec; params = nothing)` (TRD; draws seed 42's
  constants only for a spec with a kcat); TRD testset "m_OnlyA_prod's peak turnover at a
  nonzero product".
- Consumes: `_peak_turnover` (src/fitting.jl), `_testhelper_written_params`,
  `_testhelper_random_independent_params`.

**Gate:** 1.12 full 38,176 pass, pristine, JET 0.12.x; 1.13 full 38,176 pass; 1.10 full
38,168 pass, 1 broken (the fit-compile census's `@test_skip`), 0 fail. No census (no src).

The version suites run on `git archive` copies under `$P/replay/trees/`, each with the
Manifest S5 resolved for its version (`$P/s5/trees/HEAD_113/Manifest.toml`,
`$P/s5/trees/p2_110/Manifest.toml`), so `$W` keeps its 1.12 Manifest (D-S5-8).

- [ ] **Step 1: Red on Julia 1.13 and 1.10, at the base**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
for v in 113 110; do
    D=$P/replay/trees/T01_red_$v
    rm -rf $D && mkdir -p $D && git archive HEAD | tar -x -C $D
done
cp $P/s5/trees/HEAD_113/Manifest.toml $P/replay/trees/T01_red_113/
cp $P/s5/trees/p2_110/Manifest.toml $P/replay/trees/T01_red_110/
```

Background (both runs, one after the other):

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
(cd $P/replay/trees/T01_red_113 && $P/jl.sh $L/T01_red_113.log 900 +1.13 --project=. \
    -e 'using Pkg; Pkg.instantiate(); Pkg.test(julia_args=["--heap-size-hint=2500M"])')
(cd $P/replay/trees/T01_red_110 && $P/jl.sh $L/T01_red_110.log 3600 +1.10 --project=. \
    -e 'using Pkg; Pkg.instantiate(); Pkg.test(julia_args=["--heap-size-hint=2500M"])')
```

`wait_log T01_red_113`, then `wait_log T01_red_110` until it ends, then:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
grep -m1 "Unsatisfiable requirements detected for package JET" $L/T01_red_113.log
grep "Some tests did not pass" $L/T01_red_110.log
grep -A2 "Test Failed" $L/T01_red_110.log
```

Expected (`$P/proto/logs/S5_full_113_a84d194.log`, `S5_full_110_p2.log`): the 1.13 log
holds `ERROR: Unsatisfiable requirements detected for package JET [c3a54625]:` and ends
`exit 1`; the 1.10 log holds `Some tests did not pass: 38165 passed, 1 failed, 0 errored,
1 broken.`, and its one failure is `Peak turnover approaches the analytical kcat: Test
Failed at .../test/test_rate_eq_derivation.jl:720` with `Evaluated: 0.09908900906283286 ≈
0.09782929093792445 (rtol=1.0e-9)`. This 1.10 run is also unit 2's red: unit 2 changes
only test paths.

- [ ] **Step 2: Unit 1, apply 5c08ede**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 5c08ede -- Project.toml | git apply
git diff --stat
```

Expected: `Project.toml | 2 +-`.

- [ ] **Step 3: Unit 1 green on 1.12: Aqua and JET resolve and pass**

Background: `. .../replay/env.sh; focused T01_jet test/test_aqua_jet.jl`. Then
`wait_log T01_jet`, `summary T01_jet`.

Expected: `focused | 11  11` (`test/test_aqua_jet.jl | 11  11`), `exit 0`.

- [ ] **Step 4: Commit unit 1**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add Project.toml
git commit -F - <<'EOF'
Allow JET 0.12, the only JET release that runs on Julia 1.13

JET 0.11 requires Julia 1.12 and JET 0.12 requires 1.12 or later, so with
JET = "0.9, 0.11" the test environment cannot resolve on Julia 1.13, CI's
"Julia 1": Pkg.test stops with "Unsatisfiable requirements detected for
package JET". With 0.12 allowed the suite passes on 1.13 with JET 0.12.3.
On 1.12 the resolver picks 0.12 too; 1.10 keeps 0.9.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Allow JET 0.12, the only JET release that
runs on Julia 1.13`.

- [ ] **Step 5: Unit 2, apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 4a96d9a -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_rate_eq_derivation.jl | git apply
git show --format= --binary 0b5e637 -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl | git apply
git show --format= --binary 68f46f0 -- test/test_rate_eq_derivation.jl | git apply
git show --format= --binary 33b94d8 -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl | git apply
git show --format= --binary cfe9211 -- test/test_rate_eq_derivation.jl | git apply
git diff --stat
```

Expected: no `git apply` output; the stat lists FX and TRD only.

- [ ] **Step 6: Unit 2 green on 1.12**

Background: `. .../replay/env.sh; focused T01_trd test/test_rate_eq_derivation.jl`. Then
`wait_log T01_trd`, `summary T01_trd`.

Expected: `focused` passes with 0 fail, 0 error (record TRD's count). On 1.12 seed 42's
draw puts the peak at zero products, so the seeded kcat checks pass as before. The summary
lists files, not the passing testsets inside them, so the new testset "m_OnlyA_prod's peak
turnover at a nonzero product" shows only in the counts: its 2 assertions are Step 8's +2
(38,176 against Task 0's 38,174).

- [ ] **Step 7: Commit unit 2**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/mechanism_definitions_for_test_enzyme_derivation.jl \
        test/test_rate_eq_derivation.jl
git commit -F - <<'EOF'
Expect m_OnlyA_prod's true peak turnover, which a product can raise

P binds only the active conformation of m_OnlyA_prod, so at saturating S
adding P draws the enzyme out of the inactive state and raises the
turnover above the zero-product kcat whenever k2f(L - 1) > k2r(1 + L).
The peak-turnover test expected kcat, which held only for some parameter
draws: 643 of 1,000 draws exceed it, and MersenneTwister(42) draws one of
them on Julia 1.10, the compat floor, where the test failed (0.099089
against 0.097829). The search was right; the expectation was wrong.

The spec now expects the supremum along saturating S, g(c) at the
positive root of g'(c) or at c = 0, which matches the search to 1.3e-11
over the same 1,000 draws; rtol stays 1e-9. The other ten specs' kcats
hold on all 1,000 draws. The comments that claimed every peak sits at
zero products, and two false comments in the spec (kcat 2·k2f/(1+L),
group 2 EqualAI), now state the facts.

The formula is a named test helper, _testhelper_m_OnlyA_prod_peak_turnover,
beside the spec's analytical rate function. Seed 42's draw puts the peak at
zero products on Julia 1.12 and 1.13, so a testset at fixed constants
(k2f·(L − 1) > k2r·(1 + L)) checks the peak at a nonzero product on every
version; the kcat check takes its constants as a keyword and draws seed 42's
only for a spec with a kcat. The spec field's comment says what the function
receives: _testhelper_written_params' constants and Et.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Expect m_OnlyA_prod's true peak turnover,
which a product can raise`.

- [ ] **Step 8: Full suite on 1.12**

Background: `. .../replay/env.sh; full T01_full`. Then `wait_log T01_full` until it ends,
and:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
summary T01_full; pristine T01_full; grep -o "JET v[0-9.]*" $L/T01_full.log | head -1
```

Expected: `EnzymeRates.jl | 38176  38176`, `exit 0`; `pristine` prints `12` only;
`JET v0.12.x` (S5 resolved 0.12.3).

- [ ] **Step 9: Full suites on 1.13 and 1.10**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
for v in 113 110; do
    D=$P/replay/trees/T01_$v
    rm -rf $D && mkdir -p $D && git archive HEAD | tar -x -C $D
done
cp $P/s5/trees/HEAD_113/Manifest.toml $P/replay/trees/T01_113/
cp $P/s5/trees/p2_110/Manifest.toml $P/replay/trees/T01_110/
```

Background:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
(cd $P/replay/trees/T01_113 && $P/jl.sh $L/T01_full_113.log 3600 +1.13 --project=. \
    -e 'using Pkg; Pkg.instantiate(); Pkg.test(julia_args=["--heap-size-hint=2500M"])')
(cd $P/replay/trees/T01_110 && $P/jl.sh $L/T01_full_110.log 3600 +1.10 --project=. \
    -e 'using Pkg; Pkg.instantiate(); Pkg.test(julia_args=["--heap-size-hint=2500M"])')
```

`wait_log T01_full_113` and `wait_log T01_full_110` until both end; then:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
summary T01_full_113; summary T01_full_110
pristine T01_full_113; masked T01_full_113 $L/T01_full.log
pristine T01_full_110; masked T01_full_110 $P/proto/logs/S5m_full_110_cfe9211_run1.log
```

Expected: 1.13 `EnzymeRates.jl | 38176  38176`, `exit 0`; 1.10
`EnzymeRates.jl | 38168  1  38169` (pass, broken, total), `Testing EnzymeRates tests
passed`, `exit 0`. Output pristine on both, as S5 checked it: `pristine` prints `12` for
1.13 and `9` for 1.10 (1.10 skips three compile-budget `@info` lines), and nothing else;
each `masked` prints nothing (1.13 against this tree's 1.12 run; 1.10 against the
prototype tip's pristine 1.10 run, whose `[ Info:` lines phase 2 prints too). 1.13 has
never run the base code before (its test environment could not resolve): any other
failure there is a finding; stop and debug it (superpowers:systematic-debugging) before
Task 2.

---

### Task 2: Check the laws behind the string pins phase 3 renames [S1]

**Replays:** `[a3de3c4 + 7133720 + e51c2dd + cf1eb53]` (S1 §9 task 1; the commit
map's set, applied in prototype order).

**Files:**
- Modify: `test/test_allosteric_golden.jl` (a comment citing the per-spec oracle and string
  check), `test/test_dsl.jl`, `test/test_mechanism_enumeration.jl`, `test/test_types.jl`,
  `test/test_joint_solve.jl` (oracle checks)

**Interfaces:**
- Produces: `_testhelper_check_law_against_oracle` calls behind the pins Task 4 renames: TT
  (`inh`, `pingpong`, the plain, allosteric and mixed exchange fixtures), TD (the `::Inh`
  role-tag testset), TE (the dual-role testset), TJS (the kept-members testset).
- Consumes: `_testhelper_check_law_against_oracle` (FX); TRD's
  `_testhelper_test_mass_action_oracle` and `_testhelper_test_rate_equation_string` (named
  in the golden test's comment).

**Gate:** full suite 38,236 pass, pristine. No census (no src).

Test-only: green before and after. a3de3c4's oracles pass on the base code (S1 §9 task 1:
"Starts green"); so does 7133720's kept-members oracle (S1 §4: 24 of 24 in
`S1r1_base_oracles.log` with the two C9 oracles Task 3 carries; the log was lost on
2026-10-14). e51c2dd replaces a3de3c4's golden loop with a comment; cf1eb53 says what that
comment covers.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary a3de3c4 -- test/test_allosteric_golden.jl test/test_dsl.jl \
    test/test_mechanism_enumeration.jl test/test_types.jl | git apply
git show --format= --binary 7133720 -- test/test_joint_solve.jl | git apply
git show --format= --binary e51c2dd -- test/test_allosteric_golden.jl | git apply
git show --format= --binary cf1eb53 -- test/test_allosteric_golden.jl | git apply
git diff --stat
```

Expected: no `git apply` output; the stat lists the five test files.

- [ ] **Step 2: Green (test-only)**

Background: `. .../replay/env.sh; focused T02 test/test_types.jl test/test_dsl.jl
test/test_joint_solve.jl test/test_allosteric_golden.jl test/test_mechanism_enumeration.jl`.
Then `wait_log T02` until it ends, `summary T02`.

Expected: every file passes, 0 fail, 0 error; `test/test_allosteric_golden.jl | 57  57`.

- [ ] **Step 3: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_allosteric_golden.jl test/test_dsl.jl \
        test/test_mechanism_enumeration.jl test/test_types.jl test/test_joint_solve.jl
git commit -F - <<'EOF'
Check the laws behind the string pins that phase 3 renames

Phase 3 renames modifier copies and residuals in every constant name, so these
tests change pins that check only a name or a printed law. Each now also checks
its law against the exact steady state of the mass-action network, so a re-pin
cannot hide a wrong law.

The allosteric golden test does not repeat the per-spec checks of
test_rate_eq_derivation.jl; its comment says that the mass-action oracle checks
each spec's law with one catalytic subunit and that the string check ties the
printed string to the compiled law. The kept-members testset of
test_joint_solve.jl pins its kept members by step text, which the names change,
so it checks its law against the oracle too.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Check the laws behind the string pins that
phase 3 renames`.

- [ ] **Step 4: Full suite**

Background: `. .../replay/env.sh; full T02_full`. Then `wait_log T02_full` until it ends,
`summary T02_full; pristine T02_full`.

Expected: `EnzymeRates.jl | 38236  38236`, `exit 0`: Task 1's 38,176 + 60, ten oracle
calls of six assertions each (a3de3c4's golden loop and e51c2dd's removal of it cancel),
as the fix round measured on the dry replay's tree (`$P/plan/fix_counts/logs/T02_full.log`);
`pristine` prints `12` only.

---

### Task 3: C9: orient tier-3 steps and tied leads by form index [S1]

**Replays:** `[2bd4b5d + c7698fd + 2b29fd2: docs/src/deriving/re_vs_ss.md]` with departure
P1-1, `[43d7af8]` (S1 §9 tasks 2 and 3; 2b29fd2's re_vs_ss.md piece is S5's correction of
the tier-3 sentence C9 makes false, so it lands with C9). c86f61c, the start of S1 §9 task
3, lands in Task 4: its stored-orientation helper calls the three-argument
`_orient_tied_steps(reaction, oriented, tiers)` that 8b7fb88 writes, and its TJS hunk
removes the check as 8b7fb88's version of that file holds it, so it cannot run or apply
before Task 4.

**Files:**
- Modify: `src/types.jl` (`_canonical_step_direction`, `_orient_tied_steps`, docstrings),
  `src/thermo.jl` (`_priority_key`, `_joint_system`, `_kept_members`; `_solve_forward`
  deleted)
- Modify: `docs/src/deriving/re_vs_ss.md` (tier 3's last tie-break: the form order)
- Modify: `test/test_joint_solve.jl`

**Interfaces:**
- Consumes: `_form_index_key`, `_step_index_key` (thermo.jl, phase 2).
- Produces: `_canonical_step_direction` orients tier 3 by `_form_index_key` (free E first);
  `_orient_tied_steps(oriented, tiers)` ranks the lead by `_step_index_key` (no
  `lead_key`); `_priority_key(s, slot, block)` reads the stored direction;
  `_solve_forward` is gone.

**Gate:** full suite 38,496 pass, pristine; census against `base6` without
the map: `RESULT: PASS`, zero differences (C9 moves no census mechanism, S1 §5).

**Departure P1-1.** c7698fd's pin calls `_step_text(s, rxn)`, which Task 4 (8b7fb88)
introduces; at this task `_step_text` takes the step alone. The pin is applied with
`_step_text` called on the step alone; Task 4 restores c7698fd's text (P1-2). The expected
strings are the same either way: the fixture's forms bind no modifier and hold no
residual.

**Known, left to Task 4:** 43d7af8 leaves one short line in `_orient_tied_steps`'s
docstring (`its two forms). Every other tied step`); e74c8af reflows it in Task 4. Do not
reflow it here.

- [ ] **Step 1: Unit 1, apply the test paths and departure P1-1**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2bd4b5d -- test/test_joint_solve.jl | git apply
git show --format= --binary c7698fd -- test/test_joint_solve.jl | git apply
git apply <<'EOF'
--- a/test/test_joint_solve.jl
+++ b/test/test_joint_solve.jl
@@ -182,8 +182,7 @@ end
             end
         end
         m = EnzymeRates.Mechanism(tied_mm)
-        @test EnzymeRates._step_text.(first.(EnzymeRates._flat_steps(m)),
-                                      Ref(EnzymeRates.reaction(m))) ==
+        @test EnzymeRates._step_text.(first.(EnzymeRates._flat_steps(m))) ==
               ["D + S <--> DS", "E <--> D", "ES <--> DS", "E + P <--> EP", "E + S <--> ES",
                "ES <--> EP"]
         tied = only(g for g in EnzymeRates.steps(m) if length(g) == 2)
EOF
git hash-object test/test_joint_solve.jl
```

Expected: no `git apply` output; the hash is `6123cc5fffbcf25ea70821f77e0ed1ec3269d59c`.

- [ ] **Step 2: Unit 1 red**

Background: `. .../replay/env.sh; focused T03_red1 test/test_joint_solve.jl`. Then
`wait_log T03_red1`, `summary T03_red1`.

Expected (S1 §9 task 2): 6 failures, 0 errors, all in the two C9 testsets: the step
`D <--> E` is stored D → E (`from_species` is D, `k_D_to_E`/`k_E_to_D` swapped), `re`'s
`fitted_params` hold `K_D_to_E`, and the tied fixture's flat-step list reads `"D <--> E"`
and `"DS <--> ES"` and its tied steps run D → E. The two oracle checks pass.

- [ ] **Step 3: Unit 1, apply the src paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2bd4b5d -- src/thermo.jl src/types.jl | git apply
git show --format= --binary 2b29fd2 -- docs/src/deriving/re_vs_ss.md | git apply
git diff --stat
```

Expected: no `git apply` output; the stat lists docs/src/deriving/re_vs_ss.md (6 lines),
src/thermo.jl, src/types.jl and TJS.

- [ ] **Step 4: Unit 1 green**

Background: `. .../replay/env.sh; focused T03_green1 test/test_joint_solve.jl`. Then
`wait_log T03_green1`, `summary T03_green1`.

Expected: TJS passes, 0 fail, 0 error.

- [ ] **Step 5: Commit unit 1**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add docs/src/deriving/re_vs_ss.md src/thermo.jl src/types.jl test/test_joint_solve.jl
git commit -F - <<'EOF'
Orient tier-3 steps and tied-group leads by form index (C9)

The constructor stored a step that neither chemistry nor graph context orients by
the rendered names of its forms, and led a tied kinetic group by the rendered-name
key, while the joint solve orients such steps from the form earlier in the
canonical form order. Phase 3 stores what the solve reads, so the constructor now
uses the solve's rule: free E first, then the forms' base-name keys, and the tied
step ranked first by _step_index_key leads its group.

Only hand-written mechanisms change: a test-suite and docs scan of 650,092
constructions finds two, both in test_joint_solve.jl. A rapid-equilibrium step
that turns has its K renamed and inverted (K_D_to_E becomes K_E_to_D).

The two fixtures pin their whole stored flat-step list and fitted_params tuple,
not only the turned steps and the renamed K.

Docs: re_vs_ss.md's last tier-3 tie-break is the form that comes first in the
form order (free E, then conformation label, bound metabolites and residual),
not the form whose name sorts first.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Orient tier-3 steps and tied-group leads by
form index (C9)`.

- [ ] **Step 6: Unit 2, apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 43d7af8 -- test/test_joint_solve.jl | git apply
```

Expected: no output.

- [ ] **Step 7: Unit 2 red**

Background: `. .../replay/env.sh; focused T03_red2 test/test_joint_solve.jl`. Then
`wait_log T03_red2`, `summary T03_red2`.

Expected: "kind, occupancy, flag and direction rank in that order" errors on its 7 `@test`
lines, which hold the 16 `_priority_key` calls (MethodError: `_priority_key` still takes
`fwd`, so the three-argument calls match no method); the new testset "every spec stores its
steps in the direction the orientation tiers give" passes for every spec (it pins what
`_solve_forward` guaranteed, on the code before the deletion). Nothing else fails.

- [ ] **Step 8: Unit 2, apply the src paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 43d7af8 -- src/thermo.jl src/types.jl | git apply
grep -c "_solve_forward" src/thermo.jl src/types.jl
```

Expected: no `git apply` output; `src/thermo.jl:0`, `src/types.jl:0`.

- [ ] **Step 9: Unit 2 green**

Background: `. .../replay/env.sh; focused T03_green2 test/test_joint_solve.jl`. Then
`wait_log T03_green2`, `summary T03_green2`.

Expected: TJS passes, 0 fail, 0 error.

- [ ] **Step 10: Commit unit 2**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/thermo.jl src/types.jl test/test_joint_solve.jl
git commit -F - <<'EOF'
Drop _solve_forward: the joint solve reads the stored direction

With C9 the constructor stores every step in the direction the solve's tiers
give, so _solve_forward returned true for every step of every constructed
mechanism. The priority key now reads each step as stored, and the tied-group
lead key is no longer a parameter. A test pins, on every spec, that the stored
steps are a fixed point of the orientation tiers and that writing every step
the other way round, in reverse order, builds the same mechanism.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Drop _solve_forward: the joint solve reads
the stored direction`.

- [ ] **Step 11: Full suite**

Background: `. .../replay/env.sh; full T03_full`. Then `wait_log T03_full` until it ends,
`summary T03_full; pristine T03_full`.

Expected: `EnzymeRates.jl | 38496  38496`, `exit 0`: Task 2's 38,236 + 260 (the C9
testsets and the stored-orientation check over every spec), as the fix round measured on
the dry replay's tree (`$P/plan/fix_counts/logs/T03_full.log`); `pristine` prints `12`
only.

- [ ] **Step 12: Census without the map**

Background: `. .../replay/env.sh; census_run 03 base6`. Then
`wait_log $P/census/logs/cmp_replay_T03.log` until it ends, and `census_check 03`.

Expected: `RESULT: PASS` and nothing else (0 flipped stored steps).

---

### Task 4: Render phase 3's constant names [S1]

**Replays:** `[8b7fb88 + a048534 + 5e04242 + 1910833 + e74c8af + 55e310d + adeb136 + c378047
+ df09b1d + a8884ce + 2b29fd2: ping_pong.md and regulators.md]` with departure P1-2,
`[c86f61c]` (S1 §9 task 4, S1 review rounds 1 and 2, S5's 2b29fd2 but for its re_vs_ss.md,
which Task 3 takes). c86f61c is S1 §9 task 3's start; it lands here, after 8b7fb88, because
its helper calls `_orient_tied_steps(reaction, oriented, tiers)` and its TJS hunk edits
8b7fb88's text (Task 3 says so too).

**Files:**
- Modify: `src/types.jl`, `src/rate_eq_derivation.jl`, `src/mechanism_enumeration.jl`
- Modify: `ROADMAP.md` (the Planned entry for residual names that are identifiers retires)
- Modify: `docs/src/deriving/ping_pong.md`, `docs/src/deriving/regulators.md`
- Modify: `test/allosteric_ground_truth.jl`,
  `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/reference/allosteric_golden_reference.txt`, `test/test_dsl.jl`,
  `test/test_fitting_scale.jl`, `test/test_joint_solve.jl`,
  `test/test_mechanism_enumeration.jl`, `test/test_oracles.jl`,
  `test/test_rate_eq_derivation.jl`, `test/test_types.jl`

**Interfaces:**
- Produces (signatures at S1's tip a8884ce): `_is_reactant_copy(x::Metabolite, rxn)`,
  `_met_label(x, rxn)` (`Xmod` for a reactant's copy, bare otherwise),
  `_form_label(sp, rxn)` (residual `_res_plusA_minusP`), `_side_label(form, mets, rxn)`,
  `_forward_sides(s, rxn)`, `_assert_distinct_form_names(reaction, steps)` with two
  messages ("both render the name X", "share the internal name X, which orders steps and
  identifies forms"), `_step_text(s, rxn)`, `_orient_tied_steps(reaction, oriented, tiers)`,
  `_compute_alpha(rxn, flat, conc_set, species, idx, segments, step_to_K)`; constant names
  `K_EI_to_E_I`, `K_EPmodQ_to_EQ_Pmod`, `_res_plusA_minusP`. FX:
  `_testhelper_check_stored_orientation(m)` (c86f61c moves Task 3's per-spec check there
  and also runs it on the C9 and the three exchange fixtures). S3a (Task 19) later moves
  these renderers to read the mechanism instead of `rxn`.
- Consumes: Task 3's `_orient_tied_steps`, `_canonical_step_direction`.

**Gate:** full 38,575 pass (S1 r1 38,573 + 2), pristine; census against `base6` with the
map: `RESULT: PASS`, zero differences (S1 §5). Then the S1 checkpoint (Step 11).

**Departure P1-2.** Undoes P1-1: after the unit's test paths, TJS's flat-step pin calls
`_step_text(s, rxn)` again, as c7698fd wrote it. TJS then equals the prototype's blob at
c86f61c^ (8710538).

**Known:** e74c8af reflows the `_orient_tied_steps` docstring line Task 3 left short.

- [ ] **Step 1: Unit 1, apply the test paths and departure P1-2**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 8b7fb88 -- test/allosteric_ground_truth.jl \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/reference/allosteric_golden_reference.txt test/test_dsl.jl \
    test/test_fitting_scale.jl test/test_joint_solve.jl \
    test/test_mechanism_enumeration.jl test/test_oracles.jl \
    test/test_rate_eq_derivation.jl test/test_types.jl | git apply
git show --format= --binary a048534 -- test/test_rate_eq_derivation.jl | git apply
git show --format= --binary 5e04242 -- test/test_types.jl | git apply
git show --format= --binary 55e310d -- test/test_mechanism_enumeration.jl | git apply
git show --format= --binary adeb136 -- test/test_dsl.jl | git apply
git show --format= --binary df09b1d -- test/test_mechanism_enumeration.jl | git apply
git show --format= --binary a8884ce -- test/test_dsl.jl | git apply
git apply <<'EOF'
--- a/test/test_joint_solve.jl
+++ b/test/test_joint_solve.jl
@@ -183,7 +183,8 @@
             end
         end
         m = EnzymeRates.Mechanism(tied_mm)
-        @test EnzymeRates._step_text.(first.(EnzymeRates._flat_steps(m))) ==
+        @test EnzymeRates._step_text.(first.(EnzymeRates._flat_steps(m)),
+                                      Ref(EnzymeRates.reaction(m))) ==
               ["D + S <--> DS", "E <--> D", "ES <--> DS", "E + P <--> EP", "E + S <--> ES",
                "ES <--> EP"]
         tied = only(g for g in EnzymeRates.steps(m) if length(g) == 2)
EOF
git hash-object test/test_joint_solve.jl; git rev-parse c86f61c^:test/test_joint_solve.jl
```

Expected: no `git apply` output; both hashes are
`8710538ac4e35d8ea8af3ab414f1bc307bf6766a`.

- [ ] **Step 2: Unit 1 red**

Background: `. .../replay/env.sh; focused T04_red test/test_types.jl test/test_dsl.jl
test/test_rate_eq_derivation.jl test/test_oracles.jl test/test_joint_solve.jl
test/test_allosteric_golden.jl test/allosteric_ground_truth.jl test/test_fitting_scale.jl
test/test_mechanism_enumeration.jl`. Then `wait_log T04_red` until it ends,
`summary T04_red`.

Expected (S1 §9 task 4): TD passes (its changes are comments, a line split and a
`!occursin` pin that holds on the old names); every other file fails or errors, each at a
renamed name, a re-pinned string, or a call of the new two-argument `_step_text` or
three-argument `_orient_tied_steps` (MethodError). Among them: TRD "rate_equation_string
prints residual names that run as written" (a048534's title) fails its regex pin on the
`_res_plusA_minusP` names and errors where it runs the printed ping-pong law ("invalid
assignment location"); TT "a modifier-only name is bare in constant names, a reactant's
copy is Xmod" and "a residual reads _res, then _plus per added and _minus per subtracted
name" fail their exact name sets; TT "two forms whose constants would share a name are
rejected" fails (no exception: the old internal names `EPinh` and `EPmod` differ); TT "two
forms with one internal name are rejected" (substrate `Iinh` beside modifier `I`) fails on
its message (the old text says "both render the name EIinh"); the allosteric golden fails
at the line holding `G6Pmod`.

- [ ] **Step 3: Unit 1, apply the src, ROADMAP and docs paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 8b7fb88 -- src/mechanism_enumeration.jl \
    src/rate_eq_derivation.jl src/types.jl | git apply
git show --format= --binary 5e04242 -- src/types.jl | git apply
git show --format= --binary 1910833 -- src/rate_eq_derivation.jl | git apply
git show --format= --binary e74c8af -- src/types.jl | git apply
git show --format= --binary adeb136 -- ROADMAP.md src/rate_eq_derivation.jl | git apply
git show --format= --binary c378047 -- src/types.jl | git apply
git show --format= --binary 2b29fd2 -- docs/src/deriving/ping_pong.md \
    docs/src/deriving/regulators.md | git apply
git diff --stat | tail -1
```

Expected: no `git apply` output; `16 files changed` (3 src, ROADMAP.md, 2 docs, 10 test
paths).

- [ ] **Step 4: Unit 1 green**

Background: `. .../replay/env.sh; focused T04_green` with the same nine files as Step 2.
Then `wait_log T04_green` until it ends, `summary T04_green`.

Expected: every file passes, 0 fail, 0 error; TO 634, `test/test_allosteric_golden.jl`
57, `test/allosteric_ground_truth.jl` 282, `test/test_fitting_scale.jl` 104 (S2b §4's
table at 62ab877; Tasks 5-9 touch none of these four files); record the other counts (TJS
reaches 453 with unit 2).

- [ ] **Step 5: Commit unit 1**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/mechanism_enumeration.jl src/rate_eq_derivation.jl src/types.jl ROADMAP.md \
    docs/src/deriving/ping_pong.md docs/src/deriving/regulators.md \
    test/allosteric_ground_truth.jl \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/reference/allosteric_golden_reference.txt test/test_dsl.jl \
    test/test_fitting_scale.jl test/test_joint_solve.jl \
    test/test_mechanism_enumeration.jl test/test_oracles.jl \
    test/test_rate_eq_derivation.jl test/test_types.jl
git commit -F - <<'EOF'
Render phase 3's constant names: Xmod, bare modifiers, _res_plusA_minusP

Phase 3's constant names (spec section 5). A competitive-inhibitor copy of a
substrate or product renders with suffix mod (K_EPmodQ_to_EQ_Pmod); a name that
binds only as a modifier has no other role, so it renders bare (K_EI_to_E_I); a
residual renders _res_plusA_minusP, so every residual name is a Julia identifier.

The label needs the mechanism's reactants, so the chokepoint renders forms with
_form_label from the reaction, and so do the user-visible texts that name forms:
show, and the constructor, derivation and enumeration errors. The internal name
cached in Species keeps its spelling: it orders steps and identifies forms, so
this change moves no step or group order. The distinct-form check compares both,
so no two forms share a constant name or an internal name.

Every re-pinned name, factored law and the allosteric golden equal the mechanical
rename of the old text byte for byte; the census passes with the rename map.

The printed law of a ping-pong now runs as written, its destructuring line
included, and a test runs it against the compiled law. Each half of the
distinct-form check names its own clash: two forms that render one name, or two
that share the internal name (substrate Iinh beside the modifier-only I). The
docs write the new names (regulators.md, ping_pong.md); the ROADMAP's Planned
entry for residual names that are identifiers is retired, and the
rate_equation_string docstring no longer says residual names print as text.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Render phase 3's constant names: Xmod, bare
modifiers, _res_plusA_minusP`.

- [ ] **Step 6: Unit 2, apply c86f61c (test-only)**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary c86f61c -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_joint_solve.jl \
    test/test_types.jl | git apply
```

Expected: no output.

- [ ] **Step 7: Unit 2 green (test-only: green before and after)**

Background: `. .../replay/env.sh; focused T04u2 test/test_types.jl
test/test_joint_solve.jl`. Then `wait_log T04u2`, `summary T04u2`.

Expected: both pass, 0 fail, 0 error; `test/test_joint_solve.jl | 453  453`.

- [ ] **Step 8: Commit unit 2**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/mechanism_definitions_for_test_enzyme_derivation.jl \
        test/test_joint_solve.jl test/test_types.jl
git commit -F - <<'EOF'
Check stored orientation on the Tier-3 and tied fixtures

The per-spec check that stored steps are a fixed point of the orientation
tiers, and that the reversed writing builds the same mechanism, never reached
a step Tier 3 orients or a tied group: no spec holds one. The check moves to
a shared helper, which the C9 fixtures (in both writings) and the exchange
fixtures of test_types.jl now call. With Tier 3 mutated to keep the written
direction, the helper still passes on all 42 specs and fails on each of
these fixtures.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Check stored orientation on the Tier-3 and
tied fixtures`.

- [ ] **Step 9: Full suite**

Background: `. .../replay/env.sh; full T04_full`. Then `wait_log T04_full` until it ends,
`summary T04_full; pristine T04_full`.

Expected: `EnzymeRates.jl | 38575  38575`, `exit 0`; `pristine` prints `12` only.

- [ ] **Step 10: Census with the map**

Background: `. .../replay/env.sh; census_run 04 base6 census/rename_phase3.toml`. Then
`wait_log $P/census/logs/cmp_replay_T04.log` until it ends, and `census_check 04`.

Expected: `RESULT: PASS` and nothing else (every law exact after the rename, `rename
inexact 0`).

- [ ] **Step 11: S1 checkpoint against a8884ce**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
FX=test/mechanism_definitions_for_test_enzyme_derivation.jl
TRD=test/test_rate_eq_derivation.jl
diff <(git diff --name-only a8884ce HEAD -- src test docs ROADMAP.md Project.toml .claude) - <<'EOF'
Project.toml
ROADMAP.md
docs/src/deriving/ping_pong.md
docs/src/deriving/re_vs_ss.md
docs/src/deriving/regulators.md
docs/superpowers/plans/2026-10-16-graph-mechanism.md
docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md
test/mechanism_definitions_for_test_enzyme_derivation.jl
test/test_rate_eq_derivation.jl
EOF
diff <(git diff a8884ce HEAD -- Project.toml | chg) \
     <(git show --format= 5c08ede -- Project.toml | chg)
diff <(git diff a8884ce HEAD -- ROADMAP.md | chg) \
     <(git diff 4f9a34a 7205ba4 -- ROADMAP.md | chg)
diff <(git diff a8884ce HEAD -- docs/superpowers/specs | chg) \
     <(git diff 4f9a34a 7205ba4 -- docs/superpowers/specs | chg)
diff <(git diff a8884ce HEAD -- $FX $TRD | chg) \
     <(git diff 4a96d9a^ cfe9211 -- $FX $TRD | chg)
diff <(git diff a8884ce HEAD -- docs/src | chg) \
     <(git show --format= 2b29fd2 -- docs/src | chg)
git diff a8884ce HEAD -- docs/superpowers/plans | grep '^-[^-]'
```

Expected: no output. The branch differs from S1's tip only by D1 (the spec), D2 (the
plan), D3 (7205ba4's ROADMAP edits), D4 (Task 1: Project.toml, FX, TRD) and 2b29fd2 ahead
(the three deriving pages). Departures P1-1 and P1-2 cancel: TJS is not in the list.

---

### Task 5: Raw input types [S2a]

**Replays:** `[8d5e465 + 9336851]` (S2a §9 task 2).

**Files:**
- Modify: `src/types.jl` (the raw input types; the header names them)
- Modify: `test/test_types.jl`

**Interfaces:**
- Produces: `Met(name::Symbol, modifier::Bool)` (ordered by name, a modifier after the
  reactant of its name); `RawForm(conformation, bound[, added, subtracted])` and
  `RawStep(from, to, consumed, released, rapid_equilibrium)`, whose inner constructors sort
  their lists; `RawMechanism(substrates, products, modifiers, steps, group)`,
  `RawSite(ligands, multiplicity, allo_states)`,
  `RawAllostericMechanism(catalytic, cat_allo_states, catalytic_multiplicity, sites)`;
  content `==`, and a `hash` that folds the type's name, then its fields. TT:
  `_testhelper_raw_examples()`.

**Gate:** full suite 38,602 pass, pristine; census against `base6` with the
map: `RESULT: PASS`, zero differences.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 8d5e465 -- test/test_types.jl | git apply
```

Expected: no output.

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T05_red test/test_types.jl`. Then
`wait_log T05_red`, `summary T05_red`.

Expected (S2a §9 task 2): 2 errors, `UndefVarError` for `Met` in `EnzymeRates`, one in "raw
input: sorted lists, content equality and hash" and one in "a raw value's hash reads its
type's name and its fields only"; nothing else fails.

- [ ] **Step 3: Apply the src paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 8d5e465 -- src/types.jl | git apply
git show --format= --binary 9336851 -- src/types.jl | git apply
```

Expected: no output.

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T05_green test/test_types.jl`. Then
`wait_log T05_green`, `summary T05_green`.

Expected: TT passes, 0 fail, 0 error.

- [ ] **Step 5: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/types.jl test/test_types.jl
git commit -F - <<'EOF'
Add the raw input types the DSL walker returns

Met, RawForm, RawStep, RawMechanism, RawSite and RawAllostericMechanism hold
a mechanism as written: steps in written order and direction, one group label
per written statement, the declared names and the regulatory sites. The
mechanism macros will read a block into these plain values and the
constructors will build from them, so the macros and an eval-free parser can
share one walker, and tests can read the written steps.

Each constructor sorts its lists, so two writings of one form or step are
equal. == and hash compare content; Met hashes its fields, since a struct's
default hash mixes in its type's identity, which changes between precompile
sessions.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Add the raw input types the DSL walker
returns`.

- [ ] **Step 6: Full suite**

Background: `. .../replay/env.sh; full T05_full`. Then `wait_log T05_full` until it ends,
`summary T05_full; pristine T05_full`.

Expected: `EnzymeRates.jl | 38602  38602`, `exit 0`: Task 4's 38,575 + 27 (the raw input
types' tests), as the fix round measured on the dry replay's tree
(`$P/plan/fix_counts/logs/T05_full.log`); `pristine` prints `12` only.

- [ ] **Step 7: Census with the map**

Background: `. .../replay/env.sh; census_run 05 base6 census/rename_phase3.toml`. Then
`wait_log $P/census/logs/cmp_replay_T05.log` until it ends, and `census_check 05`.

Expected: `RESULT: PASS` and nothing else.

---

### Task 6: One walker reads every new-grammar block into raw input [S2a]

**Replays:** `[b3eabc4 + 62280d3 + 10a974b + fa2cc67]` (S2a §9 task 3; S2a review round 1
M6, M3, M4).

**Files:**
- Modify: `src/dsl.jl` (the walker, the routing in both mechanism macros)
- Modify: `src/types.jl` (`Mechanism(raw)`, `AllostericMechanism(raw)`,
  `_placeholder_reaction`, `_step_groups`)
- Modify: `test/test_dsl.jl`

**Interfaces:**
- Consumes: Task 5's raw types.
- Produces: `_raw_mechanism(block, allosteric::Bool)` and its helpers `_raw_site`,
  `_check_site_ligands`, `_raw_steps`, `_step_parts`, `_side_terms`, `_peel_tag!`,
  `_raw_step`, `_is_declared`, `_declared_names`, `_raw_met`, `_raw_form`,
  `_read_residual!`, `_check_conformation_label`; `Mechanism(raw::RawMechanism)` and
  `AllostericMechanism(raw::RawAllostericMechanism)`, which build the Step types through
  `_placeholder_reaction(raw)` and `_step_groups(raw)`; `_compiled_mechanism(raw)`; the
  routing `_uses_inhibitor_grammar(block, allosteric)` and `_holds_inh_tag` (Task 14
  deletes both).

**Gate:** full suite 38,680 pass, pristine; census against `base6` with the
map: `RESULT: PASS`, zero differences (the specs' blocks without inhibitor-grammar tokens
now go through the walker; S2a §5).

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary b3eabc4 -- test/test_dsl.jl | git apply
git show --format= --binary 62280d3 -- test/test_dsl.jl | git apply
git show --format= --binary 10a974b -- test/test_dsl.jl | git apply
git show --format= --binary fa2cc67 -- test/test_dsl.jl | git apply
```

Expected: no output.

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T06_red test/test_dsl.jl`. Then
`wait_log T06_red`, `summary T06_red`.

Expected (S2a §9 task 3): TD records 1 failure and 1 error. The failure: the re-pinned
missing-`steps:` test still gets "`catalytic_steps:` not specified". The error: the new
testset "DSL: modifiers, ::Mod and regulatory_site(n)" fails to expand, `LoadError` with
"@enzyme_mechanism: unknown label `modifiers:`", so none of its 10 inner testsets (which
hold the folds' pins: the full missing- and forbidden-tag messages, the single step's tag,
the tag after a parenthesized right-hand side) runs.

- [ ] **Step 3: Apply the src paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary b3eabc4 -- src/dsl.jl src/types.jl | git apply
git show --format= --binary 10a974b -- src/dsl.jl | git apply
git show --format= --binary fa2cc67 -- src/dsl.jl | git apply
```

Expected: no output.

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T06_green test/test_dsl.jl`. Then
`wait_log T06_green`, `summary T06_green`.

Expected: TD passes, 0 fail, 0 error (record the count).

- [ ] **Step 5: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/dsl.jl src/types.jl test/test_dsl.jl
git commit -F - <<'EOF'
Read mechanism blocks into raw input with one walker

_raw_mechanism walks a mechanism macro's block into RawMechanism or
RawAllostericMechanism in the grammar the rewrite settles on: modifiers: for
names that bind without reacting, X::Mod for a reactant's modifier copy, steps:
in both macros, one tag per allosteric statement, and regulatory_site(n):
lines that carry each ligand's tag. It builds no code and evaluates nothing,
so the macros and an eval-free parser can share it; its syntax errors throw at
expansion and describe this grammar.

Mechanism(raw) and AllostericMechanism(raw) build today's types with the
placeholder reaction the macros build today, and the macros embed the raw
value and call _compiled_mechanism at run time, so constructor errors still
throw when the expanded code runs.

A block with a token only the inhibitor grammar has (catalytic_steps:,
catalytic_inhibitors:, the plain regulators:, ::Inh, a tagged
allosteric_regulators: entry, a keyword regulatory_site, a parenthesized
single step with a tag) keeps today's expansion; every other block goes
through the walker. One test re-pins: an allosteric block that writes no
steps label now reports `steps:` missing.

A one-step tuple group, (step,) :: Tag, is a group, and the error for
(step) :: Tag says only where a single step's tag goes. A tag after a
parenthesized right-hand side, E(P) ⇌ (E + P) :: EqualAI, tags the step,
whose side terms then replace the sum. The tests pin the walker's full
messages for a missing and a forbidden step tag.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Read mechanism blocks into raw input with one
walker`.

- [ ] **Step 6: Full suite**

Background: `. .../replay/env.sh; full T06_full`. Then `wait_log T06_full` until it ends,
`summary T06_full; pristine T06_full`.

Expected: `EnzymeRates.jl | 38680  38680`, `exit 0`: Task 5's 38,602 + 78 (the walker's
testsets), as the fix round measured on the dry replay's tree
(`$P/plan/fix_counts/logs/T06_full.log`); `pristine` prints `12` only.

- [ ] **Step 7: Census with the map**

Background: `. .../replay/env.sh; census_run 06 base6 census/rename_phase3.toml`. Then
`wait_log $P/census/logs/cmp_replay_T06.log` until it ends, and `census_check 06`.

Expected: `RESULT: PASS` and nothing else.

---

### Task 7: oligomeric_state tuples and a metabolite named M [S2a]

**Replays:** `[7cf6b0f + 6f9511f]`, `[6075bba]` (S2a §9 tasks 4 and 5; S2a review round 1
M2).

**Files:**
- Modify: `src/dsl.jl` (`_parse_multiplicity_tuple`, the `@enzyme_reaction` docstring)
- Modify: `src/types.jl` (the EnzymeReaction constructor's multiplicity check)
- Modify: `src/rate_eq_derivation.jl` (`rate_equation`'s static parameter)
- Modify: `test/test_dsl.jl`, `test/test_types.jl`, `test/test_rate_eq_derivation.jl`

**Interfaces:**
- Produces: `_parse_multiplicity_tuple(values, label)` (`oligomeric_state:` takes one Int
  or one tuple); the EnzymeReaction constructor requires a nonempty list of Ints ≥ 1
  ("allowed_catalytic_multiplicities must be a nonempty list of Ints ≥ 1, got ...");
  `rate_equation(m::M_type, concs, params) where {M_type <: AbstractEnzymeMechanism}` (R33).

**Gate:** full suite 38,694 pass, pristine; census against `base6` with the
map: `RESULT: PASS`, zero differences.

- [ ] **Step 1: Unit 1, apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 7cf6b0f -- test/test_dsl.jl | git apply
git show --format= --binary 6f9511f -- test/test_dsl.jl test/test_types.jl | git apply
```

Expected: no output.

- [ ] **Step 2: Unit 1 red**

Background: `. .../replay/env.sh; focused T07_red1 test/test_types.jl test/test_dsl.jl`.
Then `wait_log T07_red1`, `summary T07_red1`.

Expected (S2a §9 task 4; review round 1 M2): TD errors while loading, `LoadError` with
"@enzyme_reaction: `oligomeric_state:` takes a single Int." (the testset holding
`oligomeric_state: (4, 1, 2)` and `oligomeric_state: ()` fails to expand); TT's
"EnzymeReaction takes one or more catalytic multiplicities, each ≥ 1" fails twice (the old
message for `[0, 1]`; no exception for `Int[]`).

- [ ] **Step 3: Unit 1, apply the src paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 7cf6b0f -- src/dsl.jl | git apply
git show --format= --binary 6f9511f -- src/types.jl | git apply
```

Expected: no output.

- [ ] **Step 4: Unit 1 green**

Background: `. .../replay/env.sh; focused T07_green1 test/test_types.jl
test/test_dsl.jl`. Then `wait_log T07_green1`, `summary T07_green1`.

Expected: both pass, 0 fail, 0 error.

- [ ] **Step 5: Commit unit 1**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/dsl.jl src/types.jl test/test_dsl.jl test/test_types.jl
git commit -F - <<'EOF'
Let oligomeric_state: take a tuple of multiplicities

The reaction macro reads oligomeric_state: as one Int or a tuple such as
(1, 2, 4), so a search can still span several catalytic multiplicities once
allowed_catalytic_multiplicities: goes. Two values after the label, not one
tuple, stay an error.

The EnzymeReaction constructor also requires at least one multiplicity, so
oligomeric_state: () and an empty list built directly are rejected with one
message; the old test of the multiplicity check, which passed a reaction
without products, now builds a valid reaction and pins both messages.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Let oligomeric_state: take a tuple of
multiplicities`.

- [ ] **Step 6: Unit 2, apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 6075bba -- test/test_rate_eq_derivation.jl | git apply
```

Expected: no output.

- [ ] **Step 7: Unit 2 red**

Background: `. .../replay/env.sh; focused T07_red2 test/test_rate_eq_derivation.jl`. Then
`wait_log T07_red2`, `summary T07_red2`.

Expected (S2a §9 task 5): "a metabolite named M or m reads as a concentration" errors in
the `M` mechanism's oracle check, "local variable name M conflicts with a static
parameter"; the `m` mechanism passes; nothing else fails.

- [ ] **Step 8: Unit 2, apply the src paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 6075bba -- src/rate_eq_derivation.jl | git apply
```

Expected: no output.

- [ ] **Step 9: Unit 2 green**

Background: `. .../replay/env.sh; focused T07_green2 test/test_rate_eq_derivation.jl`.
Then `wait_log T07_green2`, `summary T07_green2`.

Expected: TRD passes, 0 fail, 0 error.

- [ ] **Step 10: Commit unit 2**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl
git commit -F - <<'EOF'
Let a metabolite take the name M

rate_equation's generated body binds every metabolite name as a local, and a
local may not reuse the generator's static parameter, so a metabolite named M
broke the compiled law ("local variable name M conflicts with a static
parameter"). The static parameter is now M_type, a name with `_`, which the
reaction's name check bars from metabolites as it lands, so no name joins the
reserved list. A metabolite named after the first argument, m, already worked
and is tested beside it.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Let a metabolite take the name M`.

- [ ] **Step 11: Full suite**

Background: `. .../replay/env.sh; full T07_full`. Then `wait_log T07_full` until it ends,
`summary T07_full; pristine T07_full`.

Expected: `EnzymeRates.jl | 38694  38694`, `exit 0`: Task 6's 38,680 + 14, as the fix
round measured on the dry replay's tree (`$P/plan/fix_counts/logs/T07_full.log`);
`pristine` prints `12` only.

- [ ] **Step 12: Census with the map**

Background: `. .../replay/env.sh; census_run 07 base6 census/rename_phase3.toml`. Then
`wait_log $P/census/logs/cmp_replay_T07.log` until it ends, and `census_check 07`.

Expected: `RESULT: PASS` and nothing else.

---

### Task 8: Reserved metabolite names and names that hold _ [S2a]

**Replays:** `[51fad19 + 157ff8c + 03e2143]` (S2a §9 task 6; S2a review round 1 B1, M1).

**Files:**
- Modify: `src/types.jl` (`_RESERVED_NAMES`, `_check_metabolite_names`, the EnzymeReaction
  and AllostericMechanism constructors)
- Modify: `src/dsl.jl` (the `@enzyme_reaction` docstring)
- Modify: `src/thermo.jl` (`_dependent_param_exprs`' clash message and docstring)
- Modify: `src/rate_eq_derivation.jl`, `src/identify_rate_equation.jl` (docstrings and a
  comment: the name guards cannot fire for constructor-built input)
- Modify: `test/test_types.jl`, `test/test_rate_eq_derivation.jl`,
  `test/test_identify_rate_equation.jl`, `test/test_mechanism_enumeration.jl`

**Interfaces:**
- Produces: `_RESERVED_NAMES = (:E_total, :Keq, :L, :N, :params, :concs, :v)`;
  `_check_metabolite_names(who::String, names)` (R10, R37, R24), called by the
  EnzymeReaction constructor on every name and by the AllostericMechanism constructor on
  its site ligands; the clash error of `_dependent_param_exprs` names the conformation
  label as a fix.

**Gate:** full suite 38,839 pass, pristine; census against `base6` with the
map: `RESULT: PASS`, zero differences.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 51fad19 -- test/test_identify_rate_equation.jl \
    test/test_mechanism_enumeration.jl test/test_rate_eq_derivation.jl \
    test/test_types.jl | git apply
git show --format= --binary 157ff8c -- test/test_rate_eq_derivation.jl | git apply
```

Expected: no output.

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T08_red test/test_types.jl
test/test_rate_eq_derivation.jl test/test_mechanism_enumeration.jl
test/test_identify_rate_equation.jl`. Then `wait_log T08_red` until it ends,
`summary T08_red`.

Expected (S2a §9 task 6; review round 1 B1): TT's reserved-name and `_` assertions in its
EnzymeReaction validation testset fail (no exception, or the old "Keq and E_total are
rate-equation parameters" text), and "a regulatory-site ligand takes no reserved name and
no `_`" fails twice; TRD's "a metabolite named after a constant is rejected" and "a ligand
name holding `_` is rejected" fail (no exception), "constants that share a name raise"
fails 8 assertions on its message (S2a r1 B1: `logs/S2a_r1_red_clash.log`, lost), and "the
reserved names cover every non-metabolite name of a law" passes (R37 pins existing state);
TID's `L` assertion fails (no exception); TE passes (`LAla`, `DAla` build either way).

- [ ] **Step 3: Apply the src paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 51fad19 -- src/dsl.jl src/types.jl | git apply
git show --format= --binary 157ff8c -- src/thermo.jl | git apply
git show --format= --binary 03e2143 -- src/identify_rate_equation.jl \
    src/rate_eq_derivation.jl | git apply
```

Expected: no output.

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T08_green` with the same four files as Step 2.
Then `wait_log T08_green` until it ends, `summary T08_green`.

Expected: every file passes, 0 fail, 0 error; `test/test_rate_eq_derivation.jl | 3184
3184` and `test/test_identify_rate_equation.jl | 442  442` (S2b §4's table at 62ab877,
TRD + 2 for Task 1; Task 9 touches neither file).

- [ ] **Step 5: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/dsl.jl src/types.jl src/thermo.jl src/identify_rate_equation.jl \
    src/rate_eq_derivation.jl test/test_identify_rate_equation.jl \
    test/test_mechanism_enumeration.jl test/test_rate_eq_derivation.jl test/test_types.jl
git commit -F - <<'EOF'
Reject reserved metabolite names and names that hold `_`

A rate equation binds its constants, the locals of its code and the
concentrations in one scope, and constant names join form and metabolite
names with `_`. The EnzymeReaction constructor now rejects the names a rate
equation keeps (E_total, Keq, L, N, params, concs, v) and any name holding
`_`, so @enzyme_reaction and, through their placeholder reaction, both
mechanism macros reject them before enumeration or derivation. A mechanism
macro's reaction lists no allosteric regulator, so the AllostericMechanism
constructor checks its ligands' names the same way.

A test pins the reserved list: besides the metabolites, every spec's
generated rate_equation body and printed law bind or read only reserved names
or names holding `_`. The racemase fixture's L_Ala and D_Ala become LAla and
DAla; the testsets whose subject was a clash at derivation (L, k_EP_to_ES,
A_R, I_X) and the identify problem's L check now assert the constructor's
error.

One way remains for two constants of a law to share a name: a conformation
labelled A or I renders the state letter at the start of a step side, so an
:EqualAI step's bare constant can spell a :NonequalAI step's state-tagged one.
The guard in _dependent_param_exprs catches it; a testset builds both clashes
with @allosteric_mechanism and asserts the exact error from every entry point,
and the message names the conformation label as a fix. The docstrings of
_assert_no_name_clash, _rate_body and _equation_text, and the comment on the
identify problem's L check, say that these guards cannot fire for input the
constructors build.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends ``Reject reserved metabolite names and names
that hold `_` ``.

- [ ] **Step 6: Full suite**

Background: `. .../replay/env.sh; full T08_full`. Then `wait_log T08_full` until it ends,
`summary T08_full; pristine T08_full`.

Expected: `EnzymeRates.jl | 38839  38839`, `exit 0`: Task 7's 38,694 + 145, as the fix
round measured on the dry replay's tree (`$P/plan/fix_counts/logs/T08_full.log`);
`pristine` prints `12` only.

- [ ] **Step 7: Census with the map**

Background: `. .../replay/env.sh; census_run 08 base6 census/rename_phase3.toml`. Then
`wait_log $P/census/logs/cmp_replay_T08.log` until it ends, and `census_check 08`.

Expected: `RESULT: PASS` and nothing else.

---

### Task 9: Reject a reaction whose names spell one name two ways (R38) [S2a]

**Replays:** `[410f7c0 + 62ab877 + 0a13357 + be65542 + 0b2021d]` (S2a §9 task 7; S2a review
round 1 M5; S5's 0a13357, be65542, 0b2021d).

**Files:**
- Modify: `src/types.jl` (`_check_uniquely_decodable`, `_spelled_twice_error`, the
  EnzymeReaction constructor and docstring)
- Modify: `src/dsl.jl` (the `@enzyme_reaction` docstring and its balanced example)
- Modify: `ROADMAP.md` (the Decided-against AcCoA entry: the reaction rejects such names)
- Modify: `docs/src/getting_started.md`, `docs/src/developer.md`,
  `docs/src/identify/tutorial.md`
- Modify: `test/test_types.jl`, `test/test_mechanism_enumeration.jl`

**Interfaces:**
- Produces: `_check_uniquely_decodable(spellings::Vector{Pair{String, String}})` and
  `_spelled_twice_error(a, b)` (Sardinas-Patterson with witness); the EnzymeReaction
  constructor runs it on substrate, product and competitive-inhibitor names plus `Xmod` for
  a reactant also declared a competitive inhibitor; error text "EnzymeReaction: Ac + CoA and
  AcCoA both spell AcCoA: ...".

**Gate:** full 38,844 pass (S2a r1 38,842 + 2), pristine; census against `base6` with the
map: `RESULT: PASS`, zero differences. Then the S2a checkpoint (Step 8).

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 410f7c0 -- test/test_mechanism_enumeration.jl \
    test/test_types.jl | git apply
```

Expected: no output.

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T09_red test/test_types.jl
test/test_mechanism_enumeration.jl`. Then `wait_log T09_red` until it ends,
`summary T09_red`.

Expected (S2a §9 task 7): 8 failures. TT "a reaction's metabolite names spell each name one
way" fails its 6 rejections (no exception); TT "two forms whose constants would share a
name are rejected" gets the constructor's "both render the name EPmod" instead of "P::Mod
and Pmod both spell Pmod"; TE "a reaction whose names spell one name two ways never reaches
enumeration" gets no exception. TT "two distinct forms that render one name are rejected"
(`E(star)` beside `Estar`) passes: the constructor's check catches it today.

- [ ] **Step 3: Apply the src, ROADMAP and docs paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 410f7c0 -- ROADMAP.md src/dsl.jl src/types.jl | git apply
git show --format= --binary 62ab877 -- src/dsl.jl src/types.jl | git apply
git show --format= --binary 0a13357 -- docs/src/getting_started.md | git apply
git show --format= --binary be65542 -- docs/src/developer.md \
    docs/src/identify/tutorial.md | git apply
git show --format= --binary 0b2021d -- docs/src/identify/tutorial.md | git apply
```

Expected: no output.

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T09_green test/test_types.jl
test/test_mechanism_enumeration.jl`. Then `wait_log T09_green` until it ends,
`summary T09_green`.

Expected: `test/test_types.jl | 617  617`, `test/test_mechanism_enumeration.jl | 32528
32528` (S2b §4's table at 62ab877), 0 fail, 0 error.

- [ ] **Step 5: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add ROADMAP.md src/dsl.jl src/types.jl docs/src/getting_started.md \
    docs/src/developer.md docs/src/identify/tutorial.md \
    test/test_mechanism_enumeration.jl test/test_types.jl
git commit -F - <<'EOF'
Reject a reaction whose metabolite names spell one name two ways

Enzyme forms and their constants join the names of the metabolites they bind
without a separator, so with Ac, CoA and AcCoA the forms E(Ac, CoA) and
E(AcCoA) would share a name. Bare modifier names add such clashes (A, I and
AI), and the mechanism constructor's throw aborted a parent's whole expansion
during enumeration. The EnzymeReaction constructor now runs the
Sardinas-Patterson test on the labels forms are built from: substrate,
product and competitive-inhibitor names, plus Xmod for a reactant also
declared a competitive inhibitor; allosteric regulators bind no form. The
error spells both splits ("Ac + CoA and AcCoA both spell AcCoA").

The AcCoA enumeration test becomes a reaction-level rejection test. The
constructor keeps its form-name check for clashes that come from conformation
labels, so its test now builds E(star) beside the conformation Estar; the
reactant-copy clash test (P::Mod beside Pmod) asserts the reaction's error.

EnzymeReaction's docstring names the check, and the @enzyme_reaction example
balances its atoms (ADP joins the products). The docs say where names are
rejected: Getting Started says @enzyme_reaction rejects names spelled two ways
and names both spellings; the identify tutorial credits the L rule to the
reaction and lists the other reserved names; the developer page shows the
form-name check on a conformation label (E(star) beside Estar).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Reject a reaction whose metabolite names
spell one name two ways`.

- [ ] **Step 6: Full suite**

Background: `. .../replay/env.sh; full T09_full`. Then `wait_log T09_full` until it ends,
`summary T09_full; pristine T09_full`.

Expected: `EnzymeRates.jl | 38844  38844`, `exit 0`; `pristine` prints `12` only.

- [ ] **Step 7: Census with the map**

Background: `. .../replay/env.sh; census_run 09 base6 census/rename_phase3.toml`. Then
`wait_log $P/census/logs/cmp_replay_T09.log` until it ends, and `census_check 09`.

Expected: `RESULT: PASS` and nothing else (the six census reactions spell every name one
way; census/README.md "Reactions").

- [ ] **Step 8: S2a checkpoint against 62ab877**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
FX=test/mechanism_definitions_for_test_enzyme_derivation.jl
TRD=test/test_rate_eq_derivation.jl
diff <(git diff --name-only 62ab877 HEAD -- src test docs ROADMAP.md Project.toml .claude) - <<'EOF'
Project.toml
ROADMAP.md
docs/src/deriving/ping_pong.md
docs/src/deriving/re_vs_ss.md
docs/src/deriving/regulators.md
docs/src/developer.md
docs/src/getting_started.md
docs/src/identify/tutorial.md
docs/superpowers/plans/2026-10-16-graph-mechanism.md
docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md
test/mechanism_definitions_for_test_enzyme_derivation.jl
test/test_rate_eq_derivation.jl
EOF
diff <(git diff 62ab877 HEAD -- Project.toml | chg) \
     <(git show --format= 5c08ede -- Project.toml | chg)
diff <(git diff 62ab877 HEAD -- ROADMAP.md | chg) \
     <(git diff 4f9a34a 7205ba4 -- ROADMAP.md | chg)
diff <(git diff 62ab877 HEAD -- docs/superpowers/specs | chg) \
     <(git diff 4f9a34a 7205ba4 -- docs/superpowers/specs | chg)
diff <(git diff 62ab877 HEAD -- $FX $TRD | chg) \
     <(git diff 4a96d9a^ cfe9211 -- $FX $TRD | chg)
diff <(git diff 62ab877 HEAD -- docs/src/deriving | chg) \
     <(git show --format= 2b29fd2 -- docs/src/deriving | chg)
diff <(git diff 62ab877 HEAD -- docs/src/getting_started.md | chg) \
     <(git show --format= 0a13357 -- docs/src/getting_started.md | chg)
diff <(git diff 62ab877 HEAD -- docs/src/developer.md | chg) \
     <(git show --format= be65542 -- docs/src/developer.md | chg)
diff <(git diff 62ab877 HEAD -- docs/src/identify/tutorial.md | chg) - <<'EOF'
-- No substrate, product, or regulator may be named `L`: the search fits allosteric
-  mechanisms, whose rate equations take the conformational constant `L`.
+The reaction checks the names before a problem exists: [`@enzyme_reaction`](@ref)
+rejects a metabolite named `L`, the conformational constant of the allosteric
+mechanisms the search fits, the other names a rate equation keeps (`E_total`, `Keq`,
+`N`, `params`, `concs`, `v`) and any name holding `_`.
+
EOF
git diff 62ab877 HEAD -- docs/superpowers/plans | grep '^-[^-]'
```

Expected: no output. The branch differs from S2a's tip only by D1-D4, 2b29fd2 (Task 4), and
0a13357, be65542 and 0b2021d ahead (Task 9; the tutorial lines above are be65542 and
0b2021d together).

---

### Task 10: Write the allosteric test fixtures inline [S2b]

**Replays:** `[c80e857 + a494c9c]` (S2b §9 task 1; a494c9c is review round 1's B4, folded
in by R63 item 2).

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
  (`_testhelper_allo_from_source` deleted), `test/test_types.jl` (nine fixtures inline;
  `_testhelper_aem` deleted), `test/test_rate_eq_derivation.jl` (six fixtures inline)

**Interfaces:**
- Produces: nothing new; `_testhelper_allo_from_source` and `_testhelper_aem` are gone.
- Consumes: Task 6's walker (the inline blocks write tags on their steps).

**Gate:** full 38,844 pass (c80e857 is count-neutral: S2b §4's table), pristine. No census
(test-only).

Test-only: green before and after. Each inline block builds the helper's mechanism (11
equal with one compiled type) or both throw (4) (S2b §4; its log, `s2b/probe/inline_allo.log`, no longer exists).

**Fold:** c80e857 drops the space in `p_mix = (` (TRD, the mixed-tag regression test);
a494c9c restores it in the same commit, so no commit carries the slip (R63 item 2). A
whitespace fix: counts and c80e857's message stand.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary c80e857 -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_rate_eq_derivation.jl test/test_types.jl | git apply
git show --format= --binary a494c9c -- test/test_rate_eq_derivation.jl | git apply
grep -c "_testhelper_allo_from_source\|_testhelper_aem" \
    test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_types.jl \
    test/test_rate_eq_derivation.jl
grep -n "p_mix =" test/test_rate_eq_derivation.jl
```

Expected: no `git apply` output; each count is `0`; every `p_mix =` line reads
`p_mix = (`.

- [ ] **Step 2: Green (test-only)**

Background: `. .../replay/env.sh; focused T10 test/test_types.jl
test/test_rate_eq_derivation.jl`. Then `wait_log T10` until it ends, `summary T10`.

Expected: `test/test_types.jl | 617  617`, `test/test_rate_eq_derivation.jl | 3184  3184`,
0 fail, 0 error.

- [ ] **Step 3: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/mechanism_definitions_for_test_enzyme_derivation.jl \
        test/test_rate_eq_derivation.jl test/test_types.jl
git commit -F - <<'EOF'
Write the allosteric test fixtures inline, with tags on their steps

Fifteen tests built an allosteric mechanism by handing a plain mechanism's groups and a
list of tags to `_testhelper_allo_from_source`, directly or through `_testhelper_aem`,
the tags in written or in canonical group order. Each is now an `@allosteric_mechanism`
block with the tag on its step, so it reads without counting groups, and the three
plain blocks that fed the helper no longer need the source-order macro. Both helpers
are gone.

Each block builds the helper's mechanism (`==`, one compiled type); each error case
throws as before. A tag list longer than the groups cannot be written in the DSL, so
that case calls the constructor, as the length-mismatch test beside it does.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Write the allosteric test fixtures inline,
with tags on their steps`.

- [ ] **Step 4: Full suite**

Background: `. .../replay/env.sh; full T10_full`. Then `wait_log T10_full` until it ends,
`summary T10_full; pristine T10_full`.

Expected: `EnzymeRates.jl | 38844  38844`, `exit 0`; `pristine` prints `12` only.

---

### Task 11: Migrate the tests and docs to the modifier grammar [S2b]

**Replays:** `[c2cd9b3]` (S2b §9 task 2; a494c9c folds into Task 10, R63 item 2). c2cd9b3
is the output of the prototype's migration script; the replay applies the commit and
runs no script.

**Files:**
- Modify: `src/dsl.jl` (docstring examples only)
- Modify: `docs/src/deriving/mwc_allostery.md`, `docs/src/deriving/regulators.md`,
  `docs/src/identify/tutorial.md`
- Modify: `test/allosteric_ground_truth.jl`,
  `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/reference/allosteric_undefvar_reproducers.jl`, `test/test_dsl.jl`,
  `test/test_fitting.jl`, `test/test_fitting_scale.jl`,
  `test/test_identify_rate_equation.jl`, `test/test_joint_solve.jl`,
  `test/test_mechanism_enumeration.jl`, `test/test_oracles.jl`,
  `test/test_rate_eq_derivation.jl`, `test/test_types.jl`

**Interfaces:**
- Produces: FX's `@_testhelper_with_raw call` (returns `(mechanism, raw)`) and the spec
  field `raw::Union{RawMechanism, RawAllostericMechanism, Nothing}`, which replaces
  `source_steps` and `source_reg_sites`; TRD's `_testhelper_written_params(m, all_params;
  raw)`, the textbook bridge, reads written steps through `_step_groups(raw)` and written
  sites from `raw.sites` (R14).
- Consumes: Task 6's `_raw_mechanism` and `_step_groups`.

**Gate:** full 38,844 pass (count-neutral: S2b §4's table), pristine; census against
`base6` with the map: `RESULT: PASS`, zero differences (FX changes, and the census
evaluates the spec file's own `@_testhelper_with_raw`, which calls the tree's walker;
census/README.md).

Green before and after: the walker already reads every new-grammar block, and the src and
docs changes are docstrings and pages. The bridge compares the same 346 constants over 34
specs before and after (S2b §1; `s2b/bridge_before.tsv` equals `bridge_after.tsv`).

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary c2cd9b3 -- test/allosteric_ground_truth.jl \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/reference/allosteric_undefvar_reproducers.jl test/test_dsl.jl \
    test/test_fitting.jl test/test_fitting_scale.jl \
    test/test_identify_rate_equation.jl test/test_joint_solve.jl \
    test/test_mechanism_enumeration.jl test/test_oracles.jl \
    test/test_rate_eq_derivation.jl test/test_types.jl | git apply
```

Expected: no `git apply` output.

- [ ] **Step 2: Green before the src and docs paths**

Background: `. .../replay/env.sh; focused T11 test/test_types.jl test/test_dsl.jl
test/test_rate_eq_derivation.jl test/test_oracles.jl test/test_joint_solve.jl
test/allosteric_ground_truth.jl test/test_fitting.jl test/test_fitting_scale.jl
test/test_mechanism_enumeration.jl test/test_identify_rate_equation.jl`. Then
`wait_log T11` until it ends (several polls), `summary T11`.

Expected: 0 fail, 0 error; TT 617, TD 235, TRD 3,184, TO 634, TJS 453,
`allosteric_ground_truth.jl` 282, `test_fitting.jl` 88, `test_fitting_scale.jl` 104, TE
32,528, TID 442 (S2b §4's table; TRD + 2 for Task 1).

- [ ] **Step 3: Apply the src and docs paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary c2cd9b3 -- docs/src/deriving/mwc_allostery.md \
    docs/src/deriving/regulators.md docs/src/identify/tutorial.md src/dsl.jl | git apply
git diff --stat src/dsl.jl
```

Expected: no `git apply` output; `src/dsl.jl | 16 ++++++----------` (6 insertions, 10
deletions), all inside the two mechanism macros' docstrings.

- [ ] **Step 4: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/dsl.jl docs/src/deriving/mwc_allostery.md docs/src/deriving/regulators.md \
    docs/src/identify/tutorial.md test/allosteric_ground_truth.jl \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/reference/allosteric_undefvar_reproducers.jl test/test_dsl.jl \
    test/test_fitting.jl test/test_fitting_scale.jl \
    test/test_identify_rate_equation.jl test/test_joint_solve.jl \
    test/test_mechanism_enumeration.jl test/test_oracles.jl \
    test/test_rate_eq_derivation.jl test/test_types.jl
git commit -F - <<'EOF'
Write the tests and docs in the modifier grammar

scripts/migrate_dsl.py (phase-3 scratchpad) rewrote every DSL block in test/,
test/reference/, docs/src/ and the src docstrings, except those its skip list names:
`catalytic_steps:` becomes `steps:`; `regulators:` and `catalytic_inhibitors:` become
`modifiers:`; `X::Inh` becomes `X::Mod` on a substrate or product and bare `X` on a
modifier-only name; a regulator's tag moves from `allosteric_regulators:` onto a
`regulatory_site(n):` line, a default site becoming a line of its own;
`dead_end_inhibitors:` becomes `competitive_inhibitors:`,
`allowed_catalytic_multiplicities:` `oligomeric_state:`, and `B[N, P]` becomes `B[NP]`.
Every converted block builds what its old text builds, or both throw
(check_migration.jl: 0 mismatches over 1,026 blocks).

The skip list keeps four testsets whose subject is today's grammar, for the commit that
removes it, and seven blocks that the next commits migrate by hand.

The spec fixtures' source-order macros read only today's grammar, so they give way to
`@_testhelper_with_raw`, which returns the compiled mechanism and the DSL walker's raw
value (R14). The specs store the raw value, and the textbook bridge reads the written
steps through `_step_groups(raw)` and the written sites from it.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Write the tests and docs in the modifier
grammar`.

- [ ] **Step 5: Full suite**

Background: `. .../replay/env.sh; full T11_full`. Then `wait_log T11_full` until it ends,
`summary T11_full; pristine T11_full`.

Expected: `EnzymeRates.jl | 38844  38844`, `exit 0`; `pristine` prints `12` only.

- [ ] **Step 6: Census with the map**

Background: `. .../replay/env.sh; census_run 11 base6 census/rename_phase3.toml`. Then
`wait_log $P/census/logs/cmp_replay_T11.log` until it ends, and:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
census_check 11
sed -n 3p $P/census/logs/cmp_replay_T11.log
```

Expected: `RESULT: PASS`, then the compare header's B line, which ends `specs loaded 42,
skipped 0)`; nothing else.

---

### Task 12: Hand-migrated tests [S2b]

**Replays:** `[bb9b2e3]`, `[b0b3a4e]`, `[96846b8]`, `[9e6e64c]`, `[9227cbe]`, `[bc804d6]`,
each its own commit (S2b §9 tasks 3-7; review round 1 B2, B3).

**Files:**
- Modify: `test/test_dsl.jl` (bb9b2e3, 9e6e64c, bc804d6),
  `test/test_mechanism_enumeration.jl` (b0b3a4e, 96846b8), `test/test_types.jl` (9e6e64c,
  9227cbe)

**Interfaces:**
- Produces: exact pins of existing behavior: the walker's messages for five error tests,
  the zero-inactive-numerator fixture's oracle and its free form `E` (`K_I_ES_to_E_S`,
  `K_I_EP_to_E_P`), the `P::Mod`/`X::Mod` eval-built loops, the AllostericMechanism
  constructor's two-site message, the plain release ending in a role tag.

**Gate:** full 38,851 pass (38,844 + 6 from b0b3a4e + 1 from 9227cbe), pristine. No census
(test-only).

Each unit is test-only and pins existing behavior, so each is green on arrival. 9227cbe and
bc804d6 are each red only under a mutation of src (`$P/proto/logs/S2b_r1_red.log`); the
replay repeats both mutations and reverts them.

- [ ] **Step 1: bb9b2e3**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary bb9b2e3 -- test/test_dsl.jl | git apply
```

Background: `. .../replay/env.sh; focused T12_bb9b2e3 test/test_dsl.jl`; `wait_log
T12_bb9b2e3`; `summary T12_bb9b2e3`. Expected: `test/test_dsl.jl | 235  235`.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_dsl.jl
git commit -F - <<'EOF'
Migrate the DSL error tests whose message or point moves with the grammar

Five error tests in test_dsl.jl lose their meaning or their pinned text once their
blocks are written in the modifier grammar, so they are migrated by hand:

- an untagged allosteric regulator is, in the new grammar, one on no site line;
  the test now pins the walker's "sits at no regulatory site" error;
- a ligand on two sites is now a walker error, so the test pins its text, and the
  testset's title says that the macro checks names and sites;
- the repeated label is `steps:`;
- `R::Inh` on an allosteric-only name becomes `R::Mod`, which the walker rejects as
  neither substrate nor product; stripping the tag would only repeat the testset's
  first assertion;
- the step that ends in a role tag and has no state now ends in a reactant's `S::Mod`,
  the one role tag the grammar keeps, and the test pins the step's text, so it still
  shows that the role tag is not read as the state.

Each old block and its new form both throw. The old blocks raise the inhibitor
grammar's errors: an `allosteric_regulators:` entry needs a `name::Tag` annotation;
ligand A appears in two distinct regulatory sites; `catalytic_steps:` given more than
once; `R::Inh` names `R`, which `catalytic_inhibitors:` does not declare; step
`E(I::Inh) ⇌ E + I::Inh` is missing its state annotation. The new blocks raise the
walker messages the tests now pin.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Migrate the DSL error tests whose message or
point moves with the grammar`.

- [ ] **Step 2: b0b3a4e**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary b0b3a4e -- test/test_mechanism_enumeration.jl | git apply
```

Background: `. .../replay/env.sh; focused T12_b0b3a4e test/test_mechanism_enumeration.jl`;
`wait_log T12_b0b3a4e` until it ends; `summary T12_b0b3a4e`. Expected:
`test/test_mechanism_enumeration.jl | 32534  32534` (the oracle adds 6: 3 points, 2 tests
each).

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_mechanism_enumeration.jl
git commit -F - <<'EOF'
Check the zero-inactive-numerator fixture's law against the oracle

The testset pins two constant names and an empty inactive numerator, and the next
commit renames the fixture's free form, which renames both constants. Before a test
that checks only names changes them, it gains a value check (spec section 14): the
fixture's law now agrees with the exact steady state of its mass-action network.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Check the zero-inactive-numerator fixture's
law against the oracle`.

- [ ] **Step 3: 96846b8**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 96846b8 -- test/test_mechanism_enumeration.jl | git apply
grep -c "E_c" test/test_mechanism_enumeration.jl
```

Expected: no `git apply` output; the count is `0`. Background: `. .../replay/env.sh;
focused T12_96846b8 test/test_mechanism_enumeration.jl`; `wait_log T12_96846b8` until it
ends; `summary T12_96846b8`. Expected: `test/test_mechanism_enumeration.jl | 32534
32534`.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_mechanism_enumeration.jl
git commit -F - <<'EOF'
Rename the zero-inactive-numerator fixture's free form E_c to E

A conformation label holds no `_` (R40): constant names join labels with `_`, so a
label `E_c` could make two steps' constants spell one name. The walker rejects it, so
the fixture, the one mechanism in the tests whose free form is not E, writes its free
form as E (R28), in the modifier grammar. Its two pinned constants become
K_I_ES_to_E_S and K_I_EP_to_E_P. Renaming E_c to E in the old mechanism's parameters
and fitted parameters gives the new mechanism's, and the two laws agree at five
points to a relative difference of 0. The law check added before the rename still
passes.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Rename the zero-inactive-numerator fixture's
free form E_c to E`.

- [ ] **Step 4: 9e6e64c**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9e6e64c -- test/test_dsl.jl test/test_types.jl | git apply
```

Background: `. .../replay/env.sh; focused T12_9e6e64c test/test_types.jl
test/test_dsl.jl`; `wait_log T12_9e6e64c`; `summary T12_9e6e64c`. Expected:
`test/test_types.jl | 617  617`, `test/test_dsl.jl | 235  235`.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_dsl.jl test/test_types.jl
git commit -F - <<'EOF'
Migrate the eval-built DSL loops that interpolate role-tagged lines

Two loops build their blocks with `eval` and splice in a line that writes `X::Inh`,
which the migration script cannot see. The role-change loop of test_types.jl now
declares P under `modifiers:` and splices `P::Mod` lines; the step constructor raises
the same role-change error for each. The X::Inh loop of test_dsl.jl becomes the X::Mod
loop and pins the walker's message for a substrate or product written `X::Mod` that
`modifiers:` does not list. Each old line and its new form both throw: each of the
role-change loop's three steps raises the same message in both forms, and each X::Mod
line raises the walker's message for an undeclared modifier where its X::Inh line
raised the inhibitor grammar's message for an undeclared inhibitor.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Migrate the eval-built DSL loops that
interpolate role-tagged lines`.

- [ ] **Step 5: 9227cbe, with its red under a mutation**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9227cbe -- test/test_types.jl | git apply
sed -i 's/^        dup === nothing ||$/        true ||/' src/types.jl
git diff --stat src/types.jl
```

Expected: `src/types.jl | 2 +-`. Background: `. .../replay/env.sh; focused
T12_9227cbe_red test/test_types.jl`; `wait_log T12_9227cbe_red`; `summary
T12_9227cbe_red`. Expected: 1 failure, in "AllostericEnzymeMechanism lift validators":
"No exception thrown" for the two-site construction (`S2b_r1_red.log`). Then revert the
mutation and run green:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git checkout -- src/types.jl
git status --short
```

Expected: ` M test/test_types.jl` only. Background: `. .../replay/env.sh; focused
T12_9227cbe test/test_types.jl`; `wait_log T12_9227cbe`; `summary T12_9227cbe`.
Expected: `test/test_types.jl | 618  618`.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_types.jl
git commit -F - <<'EOF'
Test that AllostericMechanism rejects a ligand on two regulatory sites

The constructor rejects a ligand that sits at two sites, whose Kreg names would
collide. Its one test wrote the two sites in the DSL; in the modifier grammar the
walker rejects such a block while the macro expands, so the test now pins the
walker's message and nothing reached the constructor's check. A direct construction
with ligand A at sites of multiplicity 2 and 4 now pins the constructor's message.
With the check disabled, the assertion fails: no exception is thrown.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Test that AllostericMechanism rejects a
ligand on two regulatory sites`.

- [ ] **Step 6: bc804d6, with its red under a mutation**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary bc804d6 -- test/test_dsl.jl | git apply
sed -i 's/t.args\[2\] ∉ (:Mod, :Inh) || return nothing/t.args[2] ∉ (:Inh,) || return nothing/' \
    src/dsl.jl
git diff --stat src/dsl.jl
```

Expected: `src/dsl.jl | 2 +-`. Background: `. .../replay/env.sh; focused T12_bc804d6_red
test/test_dsl.jl`; `wait_log T12_bc804d6_red`; `summary T12_bc804d6_red`. Expected: TD
errors with "@enzyme_mechanism: tag annotation on `E(S::Mod) ⇌ E + S::Mod` is not allowed"
(`S2b_r1_red.log`: with `_peel_tag!` reading `::Mod` as a step tag, the release no longer
expands). Then revert the mutation and run green:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git checkout -- src/dsl.jl
git status --short
```

Expected: ` M test/test_dsl.jl` only. Background: `. .../replay/env.sh; focused T12_bc804d6
test/test_dsl.jl`; `wait_log T12_bc804d6`; `summary T12_bc804d6`. Expected:
`test/test_dsl.jl | 235  235`.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_dsl.jl
git commit -F - <<'EOF'
End the plain release in "::Mod is never a step tag" with a role tag

The testset's plain pair checked that a release that ends in a role tag equals the
binding it reverses, so the walker never reads the tag as a step tag. The migration
wrote its modifier-only I bare, which leaves no role tag to misread. The pair now
binds and releases the reactant copy S::Mod, as the allosteric pair above it does.
If `_peel_tag!` read `::Mod` as a tag, the release would now fail to expand ("tag
annotation on `E(S::Mod) ⇌ E + S::Mod` is not allowed"); the bare pair still built
and compared equal.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `End the plain release in "::Mod is never a
step tag" with a role tag`.

- [ ] **Step 7: Full suite**

Background: `. .../replay/env.sh; full T12_full`. Then `wait_log T12_full` until it ends,
`summary T12_full; pristine T12_full`.

Expected: `EnzymeRates.jl | 38851  38851`, `exit 0`; `pristine` prints `12` only.

---

### Task 13: Describe the modifier grammar in titles, comments, docstrings and docs [S2b]

**Replays:** `[3c0c040]`, `[588e99f]`, `[fcfd5e6]`, `[9e41af7 + c34c525]`, `[756abee]`
(S2b §9 task 8; review round 1 m1, m3).

**Files:**
- Modify: `src/dsl.jl` (the two mechanism macros' docstrings)
- Modify: `docs/src/deriving/mwc_allostery.md`, `docs/src/deriving/regulators.md`,
  `docs/src/developer.md`
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`, `test/test_dsl.jl`,
  `test/test_mechanism_enumeration.jl`, `test/test_types.jl`

**Interfaces:**
- Produces: docstrings and docs in the modifier grammar; `@_testhelper_with_raw`'s docstring
  names its argument `call`; TD's tests of the `@enzyme_reaction` forms Task 15 removes
  (`dead_end_inhibitors:`, `B[N,P]`, `allowed_catalytic_multiplicities:`) are back as
  written.

**Gate:** full 38,851 pass (S2b r1 38,849 + 2), pristine. No census (prose and tests). Then
the S2b checkpoint (Step 7).

Prose units are green before and after; 756abee restores three `@enzyme_reaction` forms'
tests the migration rewrote (count-neutral) and comes after the prose commits whose TD lines
it reads. Message fix: 9e41af7 named the argument `macro_call`; with c34c525 folded in, the
docstring says `call`, and so does the message.

- [ ] **Step 1: 3c0c040**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 3c0c040 -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_dsl.jl \
    test/test_mechanism_enumeration.jl test/test_types.jl | git apply
```

Background: `. .../replay/env.sh; focused T13_3c0c040 test/test_types.jl test/test_dsl.jl
test/test_mechanism_enumeration.jl`; `wait_log T13_3c0c040` until it ends; `summary
T13_3c0c040`. Expected: TT 618, TD 235, TE 32,534, 0 fail, 0 error.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_dsl.jl \
        test/test_mechanism_enumeration.jl test/test_types.jl
git commit -F - <<'EOF'
Write the tests' titles and comments in the modifier grammar

Test titles and comments that quote DSL named the labels and the role tag of the
blocks beside them before the migration: `X::Inh`, `catalytic_inhibitors:`,
`dead_end_inhibitors:`, `allowed_catalytic_multiplicities:`,
`regulatory_site(multiplicity = N)`. They now quote the modifier grammar those blocks
are written in. Comments that name the competitive-inhibitor role, which the types
still carry, are unchanged.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Write the tests' titles and comments in the
modifier grammar`.

- [ ] **Step 2: 588e99f**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 588e99f -- test/test_dsl.jl | git apply
git show --format= --binary 588e99f -- src/dsl.jl | git apply
```

Background: `. .../replay/env.sh; focused T13_588e99f test/test_dsl.jl`; `wait_log
T13_588e99f`; `summary T13_588e99f`. Expected: `test/test_dsl.jl | 235  235`.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/dsl.jl test/test_dsl.jl
git commit -F - <<'EOF'
Describe the modifier grammar in the mechanism macros' docstrings

The migration rewrote the examples of @enzyme_mechanism and @allosteric_mechanism in
the modifier grammar, while the bullets under them still described `regulators:`,
`catalytic_inhibitors:`, `catalytic_steps:`, `X::Inh`, tagged `allosteric_regulators:`
entries, keyword sites with default sites, and the label `E_c`. They now describe the
grammar the examples use: `modifiers:`, `X::Mod` on a substrate or product only,
`steps:`, bare regulator names with each tag on its `regulatory_site(n):` line, one
site per regulator, and conformation labels without `_`.

A test of the parenthesized one-step group cited the docstring's example, which no
longer holds that form; its comment now says what it checks.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Describe the modifier grammar in the
mechanism macros' docstrings`.

- [ ] **Step 3: fcfd5e6 (docs only, no run)**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary fcfd5e6 -- docs/src/deriving/mwc_allostery.md \
    docs/src/deriving/regulators.md docs/src/developer.md | git apply
git add docs/src/deriving/mwc_allostery.md docs/src/deriving/regulators.md \
        docs/src/developer.md
git commit -F - <<'EOF'
Quote the modifier grammar in the docs' prose

Three sentences quote DSL beside examples that the migration rewrote: how to declare a
regulator (`regulators: <name>`), what a site's multiplicity means
(`regulatory_site(multiplicity = m)`), and the binding the step constructor rejects
(`E + P::Inh ⇌ E(P)`). They now write `modifiers: <name>`, `regulatory_site(m)` and
`E + P::Mod ⇌ E(P)`.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no `git apply` or status output; the log line ends `Quote the modifier grammar in
the docs' prose`.

- [ ] **Step 4: 9e41af7 + c34c525**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9e41af7 -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl | git apply
git show --format= --binary c34c525 -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl | git apply
grep -c "macro_call" test/mechanism_definitions_for_test_enzyme_derivation.jl
```

Expected: no `git apply` output; the count is `0`. Background: `. .../replay/env.sh;
focused T13_9e41af7 test/test_dsl.jl` (FX loads); `wait_log T13_9e41af7`; `summary
T13_9e41af7`. Expected: `test/test_dsl.jl | 235  235`.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -F - <<'EOF'
Name the raw-value test macro's argument instead of sketching a block

The usage lines of `@_testhelper_with_raw` wrote `@enzyme_mechanism begin ... end`,
which tools that scan the test files for DSL blocks (the migration script's block
finder among them) count as two more blocks. The usage now names the argument as the
macro's parameter does, `call`, and the text says which calls it takes.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Name the raw-value test macro's argument
instead of sketching a block`.

- [ ] **Step 5: 756abee**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 756abee -- test/test_dsl.jl | git apply
git hash-object test/test_dsl.jl; git rev-parse 756abee:test/test_dsl.jl
```

Expected: no `git apply` output; both hashes are
`25a8bbd4dde162d798ebf4194e21afc0454b07cc` (TD equals the prototype's). Background:
`. .../replay/env.sh; focused T13_756abee test/test_dsl.jl`; `wait_log T13_756abee`;
`summary T13_756abee`.
Expected: `test/test_dsl.jl | 235  235`.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_dsl.jl
git commit -F - <<'EOF'
Restore the tests of three @enzyme_reaction forms the migration rewrote

Four @enzyme_reaction tests have the reaction grammar as their subject: that
`dead_end_inhibitors:` emits competitive inhibitors, that it rejects a type tag, that
an atom bracket may hold several formulas (`B[N,P]`), and that
`allowed_catalytic_multiplicities:` takes a tuple. The migration rewrote them in the
labels and brackets that replace these forms, so the forms lost their tests while the
parser still reads them, and the regulator-kinds test came to pin that two
`competitive_inhibitors:` lines accumulate. They are back as written, beside the
mechanism-grammar testsets kept for the commit that removes the forms, which turns
them into rejection tests.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: no status output; the log line ends `Restore the tests of three @enzyme_reaction
forms the migration rewrote`.

- [ ] **Step 6: Full suite**

Background: `. .../replay/env.sh; full T13_full`. Then `wait_log T13_full` until it ends,
`summary T13_full; pristine T13_full`.

Expected: `EnzymeRates.jl | 38851  38851`, `exit 0`; `pristine` prints `12` only.

- [ ] **Step 7: S2b checkpoint against c34c525**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
FX=test/mechanism_definitions_for_test_enzyme_derivation.jl
TRD=test/test_rate_eq_derivation.jl
diff <(git diff --name-only c34c525 HEAD -- src test docs ROADMAP.md Project.toml .claude) - <<'EOF'
Project.toml
ROADMAP.md
docs/src/deriving/ping_pong.md
docs/src/deriving/re_vs_ss.md
docs/src/deriving/regulators.md
docs/src/developer.md
docs/src/getting_started.md
docs/src/identify/tutorial.md
docs/superpowers/plans/2026-10-16-graph-mechanism.md
docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md
test/mechanism_definitions_for_test_enzyme_derivation.jl
test/test_rate_eq_derivation.jl
EOF
diff <(git diff c34c525 HEAD -- Project.toml | chg) \
     <(git show --format= 5c08ede -- Project.toml | chg)
diff <(git diff c34c525 HEAD -- ROADMAP.md | chg) \
     <(git diff 4f9a34a 7205ba4 -- ROADMAP.md | chg)
diff <(git diff c34c525 HEAD -- docs/superpowers/specs | chg) \
     <(git diff 4f9a34a 7205ba4 -- docs/superpowers/specs | chg)
diff <(git diff c34c525 HEAD -- $FX $TRD | chg) \
     <(git diff 4a96d9a^ cfe9211 -- $FX $TRD | chg)
diff <(git diff c34c525 HEAD -- docs/src/deriving | chg) \
     <(git show --format= 2b29fd2 -- docs/src/deriving | chg)
diff <(git diff c34c525 HEAD -- docs/src/getting_started.md | chg) \
     <(git show --format= 0a13357 -- docs/src/getting_started.md | chg)
diff <(git diff c34c525 HEAD -- docs/src/developer.md | chg) \
     <(git show --format= be65542 -- docs/src/developer.md | chg)
diff <(git diff c34c525 HEAD -- docs/src/identify/tutorial.md | chg) - <<'EOF'
-- No substrate, product, or regulator may be named `L`: the search fits allosteric
-  mechanisms, whose rate equations take the conformational constant `L`.
+The reaction checks the names before a problem exists: [`@enzyme_reaction`](@ref)
+rejects a metabolite named `L`, the conformational constant of the allosteric
+mechanisms the search fits, the other names a rate equation keeps (`E_total`, `Keq`,
+`N`, `params`, `concs`, `v`) and any name holding `_`.
+
EOF
git diff c34c525 HEAD -- docs/superpowers/plans | grep '^-[^-]'
```

Expected: no output: as at S2a, D1-D4, 2b29fd2, 0a13357, be65542 and 0b2021d. The S2b
units add no departure (a494c9c's fold leaves TRD equal to the tip's apart from D4).

## Tasks 14-28: S2c grammar removal, S3a boundary check, S3b content order, S3c test reader

Every `git apply` and checkpoint block of Tasks 14-28 ran as written, in order, on a
scratch index built from 7205ba4 through Task 13 (`$P/plan/check_p2/sim_p2.py`): every
apply was clean, each unit applied every path of its commits once, and the four
checkpoint checks printed nothing.

### Task 14: Reject the inhibitor grammar's mechanism forms; read every block with the walker [S2c]

**Replays:** `[3b28fca]`, `[68a4961]`, `[fa09735]`, `[9399c58 + 1748df6]`, `[822d590]`,
`[9736179]`, `[b361684]`: seven commits, in this order (S2c D-S2c-1: a form's test goes red
only while the router still holds that form; the parser goes last, dead).

**Files:**
- Modify: `src/dsl.jl` (the router `_uses_inhibitor_grammar`, `_holds_inh_tag`, the routing
  line in each mechanism macro, and the inhibitor-grammar parser go)
- Test: `test/test_dsl.jl`

**Interfaces:**
- Consumes: `_raw_mechanism(block, allosteric::Bool)`, `_compiled_mechanism(raw)` (Task 6).
- Produces: `@enzyme_mechanism` and `@allosteric_mechanism` expand every block to
  `:($_compiled_mechanism($(_raw_mechanism(block, allosteric))))`. Gone:
  `_uses_inhibitor_grammar`, `_holds_inh_tag`, `_parse_mechanism_body`,
  `_parse_steps_block`, `_peel_step_tag!`, `_step_expr`, `_species_expr`,
  `_is_conformation_shape`, `_metabolite_expr`, `_walk_residual_expr`,
  `_parse_regulatory_site`, `_build_reg_sites_expr`, `_mechanism_reaction_expr`. Kept:
  `_reaction_expr`, `_parse_labeled_line`, `_bare_symbols_from_values`, `_positive_int`.

**Gate:** F 14 = 38,851; C 14 base6: zero differences.

The S2c red and green logs were lost in the 2026-10-14 outage; the reds below are S2_report.md
section 2's.

**Unit 1: 3b28fca** (`catalytic_steps:`, `catalytic_inhibitors:`, plain `regulators:`)

- [ ] **Step 1: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 3b28fca -- test/test_dsl.jl | git apply
```

- [ ] **Step 2: Red.** Rd `$L/T14_1_red.log` `test_dsl.jl`. Expected: 5 failures, all in
  "labels: the valid ones, each single-valued label once with one value" (the plain
  `regulators:` block builds a mechanism; the other four blocks raise the inhibitor
  grammar's errors).
- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 3b28fca -- src/dsl.jl | git apply
```

- [ ] **Step 4: Green.** G `$L/T14_1_green.log` `test_dsl.jl`: 233.
- [ ] **Step 5: Commit.** K `3b28fca` `src/dsl.jl test/test_dsl.jl` ("Reject
  `catalytic_steps:`, `catalytic_inhibitors:` and plain `regulators:`").

**Unit 2: 68a4961** (`::Inh`)

- [ ] **Step 6: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 68a4961 -- test/test_dsl.jl | git apply
```

- [ ] **Step 7: Red.** Rd `$L/T14_2_red.log` `test_dsl.jl`. Expected: 2 failures in
  "::Mod marks a substrate or product under modifiers: only" (both blocks reach the
  inhibitor-grammar parser, which rejects `modifiers:`).
- [ ] **Step 8: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 68a4961 -- src/dsl.jl | git apply
```

- [ ] **Step 9: Green.** G `$L/T14_2_green.log` `test_dsl.jl`: 234.
- [ ] **Step 10: Commit.** K `68a4961` `src/dsl.jl test/test_dsl.jl` ("Reject the role tag
  `::Inh`").

**Unit 3: fa09735** (a tagged `allosteric_regulators:` entry)

- [ ] **Step 11: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary fa09735 -- test/test_dsl.jl | git apply
```

- [ ] **Step 12: Red.** Rd `$L/T14_3_red.log` `test_dsl.jl`. Expected: 2 failures in
  "every allosteric regulator sits at exactly one regulatory site" (the parser rejects
  their `steps:` label).
- [ ] **Step 13: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary fa09735 -- src/dsl.jl | git apply
```

- [ ] **Step 14: Green.** G `$L/T14_3_green.log` `test_dsl.jl`: 235.
- [ ] **Step 15: Commit.** K `fa09735` `src/dsl.jl test/test_dsl.jl` ("Reject a tag on an
  `allosteric_regulators:` entry").

**Unit 4: 9399c58 + 1748df6** (a keyword multiplicity; 1748df6 rewraps the 93-column
message line 9399c58 wrote)

- [ ] **Step 16: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9399c58 -- test/test_dsl.jl | git apply
git show --format= --binary 1748df6 -- test/test_dsl.jl | git apply
```

- [ ] **Step 17: Red.** Rd `$L/T14_4_red.log` `test_dsl.jl`. Expected: 3 failures in
  "regulatory_site takes exactly one multiplicity" (its three keyword shapes reach the
  parser, which rejects their untagged `allosteric_regulators:` entry).
- [ ] **Step 18: Apply the source** (1748df6 has none)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9399c58 -- src/dsl.jl | git apply
```

- [ ] **Step 19: Green.** G `$L/T14_4_green.log` `test_dsl.jl`: 235.
- [ ] **Step 20: Commit.** K `9399c58` `src/dsl.jl test/test_dsl.jl` ("Reject a keyword
  multiplicity on a `regulatory_site` line"; the fold changes nothing the message says).

**Unit 5: 822d590** (a tag after a parenthesized single step)

- [ ] **Step 21: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 822d590 -- test/test_dsl.jl | git apply
```

- [ ] **Step 22: Red.** Rd `$L/T14_5_red.log` `test_dsl.jl`. Expected: 1 failure in
  "allosteric steps carry one tag per statement, plain steps none" (the block now expands
  the macro and reaches the parser, which rejects `steps:`).
- [ ] **Step 23: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 822d590 -- src/dsl.jl | git apply
```

- [ ] **Step 24: Green.** G `$L/T14_5_green.log` `test_dsl.jl`: 235.
- [ ] **Step 25: Commit.** K `822d590` `src/dsl.jl test/test_dsl.jl` ("Reject a tag after a
  parenthesized single step").

**Unit 6: 9736179** (the `begin ... end` site body, the router)

- [ ] **Step 26: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9736179 -- test/test_dsl.jl | git apply
```

- [ ] **Step 27: Red.** Rd `$L/T14_6_red.log` `test_dsl.jl`. Expected: 1 failure in
  "every allosteric regulator sits at exactly one regulatory site" (`_testhelper_sites`
  through the macro; the body block reaches the parser, which rejects its untagged
  `allosteric_regulators:` entries).
- [ ] **Step 28: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9736179 -- src/dsl.jl | git apply
```

- [ ] **Step 29: Green.** G `$L/T14_6_green.log` `test_dsl.jl test_types.jl`: 235, 618.
- [ ] **Step 30: Commit.** K `9736179` `src/dsl.jl test/test_dsl.jl` ("Read every mechanism
  block with the walker").

**Unit 7: b361684** (the dead parser; green before and after)

- [ ] **Step 31: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary b361684 -- src/dsl.jl | git apply
git diff --stat
```

  Expected: `src/dsl.jl | 408 ---...` and `1 file changed, 408 deletions(-)`.
- [ ] **Step 32: Green.** G `$L/T14_7_green.log` `test_dsl.jl test_types.jl`: 235, 618.
- [ ] **Step 33: Commit.** K `b361684` `src/dsl.jl` ("Delete the inhibitor-grammar parser").

**Task end**

- [ ] **Step 34: Full suite.** F 14: 38,851.
- [ ] **Step 35: Census.** C 14 `base6`: zero differences (the walker builds the same
  mechanisms the parser built).

### Task 15: @enzyme_reaction: one form per label [S2c]

**Replays:** `[961ec20]`, `[37c6a8b + a84d194]`, `[f8d0199]`, `[6c85d75]`.

**Files:**
- Modify: `src/dsl.jl` (`_parse_reaction_block`, `_parse_atom_bracket_entries`,
  `_parse_multiplicity_tuple`, the `@enzyme_reaction` docstring),
  `docs/src/identify/enumeration_engine.md` (a84d194)
- Test: `test/test_dsl.jl`

**Interfaces:**
- Consumes: `_parse_multiplicity_tuple` (Task 7).
- Produces: `@enzyme_reaction` takes exactly the labels `substrates`, `products`,
  `competitive_inhibitors`, `allosteric_regulators`, `oligomeric_state`,
  `shared_catalytic_site`; any other is reported as "@enzyme_reaction: unknown label
  `<label>:`. Valid labels: ..."; a second `oligomeric_state:` line reports
  "@enzyme_reaction: `oligomeric_state:` given more than once."; a bracket holds exactly one
  formula ("`substrates:` lists each reactant with one chemical formula in brackets, such
  as `S[C6H12O6]`; got `...`."); `oligomeric_state: ()` is rejected while the macro expands
  ("`oligomeric_state:` takes one Int or a tuple of one or more Ints, such as `(1, 2, 4)`;
  got ().").

**Gate:** F 15 = 38,854 (Task 14's 38,851: 961ec20 +1, 37c6a8b −1, f8d0199 +3); C 15 base6:
zero differences.

**Unit 1: 961ec20** (`dead_end_inhibitors:`)

- [ ] **Step 1: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 961ec20 -- test/test_dsl.jl | git apply
```

- [ ] **Step 2: Red.** Rd `$L/T15_1_red.log` `test_dsl.jl`. Expected: 2 failures in
  "@enzyme_reaction names its labels for any other" (`dead_end_inhibitors:` builds a
  reaction; the `regulators:` message still lists it as valid).
- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 961ec20 -- src/dsl.jl | git apply
```

- [ ] **Step 4: Green.** G `$L/T15_1_green.log` `test_dsl.jl`: 236.
- [ ] **Step 5: Commit.** K `961ec20` `src/dsl.jl test/test_dsl.jl` ("Reject
  `dead_end_inhibitors:` in @enzyme_reaction").

**Unit 2: 37c6a8b + a84d194** (`allowed_catalytic_multiplicities:`; the engine page names
the label)

- [ ] **Step 6: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 37c6a8b -- test/test_dsl.jl | git apply
```

- [ ] **Step 7: Red.** Rd `$L/T15_2_red.log` `test_dsl.jl`. Expected: 4 failures (the
  report's 3 + 1): 3 in "@enzyme_reaction names its labels for any other" (both messages
  still list `allowed_catalytic_multiplicities`; that label builds a reaction) and 1 in
  "@enzyme_reaction sets the catalytic multiplicities once" (the old "sets the catalytic"
  message).
- [ ] **Step 8: Apply the source and the docs**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 37c6a8b -- src/dsl.jl | git apply
git show --format= --binary a84d194 -- docs/src/identify/enumeration_engine.md | git apply
```

- [ ] **Step 9: Green.** G `$L/T15_2_green.log` `test_dsl.jl`: 235.
- [ ] **Step 10: Commit** (message edited: one paragraph for a84d194)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M docs/src/identify/enumeration_engine.md, M src/dsl.jl, M test/test_dsl.jl
git add docs/src/identify/enumeration_engine.md src/dsl.jl test/test_dsl.jl
cat > $P/replay/msg <<'EOF'
Reject `allowed_catalytic_multiplicities:` in @enzyme_reaction

`oligomeric_state:` takes one Int or a tuple, so it sets every catalytic multiplicity
the search may assign; `allowed_catalytic_multiplicities:` is now an unknown label.
With one label left, a second `oligomeric_state:` line reports "`oligomeric_state:`
given more than once", as the mechanism macros do for their single-valued labels.

Tests: the label test adds `allowed_catalytic_multiplicities:`; before this change it
was read as a multiplicity label. The three tests of a second multiplicity line, which
mixed the two labels in either order, become one test of a second `oligomeric_state:`
line. The tuple test writes `oligomeric_state: (1, 2, 4)`.

Docs: the enumeration engine page names `oligomeric_state:` as the label whose
multiplicities the allosteric promotion runs over; it named only the field,
`allowed_catalytic_multiplicities`, which reads like the label the macro now rejects.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```


**Unit 3: f8d0199** (one formula per bracket)

- [ ] **Step 11: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary f8d0199 -- test/test_dsl.jl | git apply
```

- [ ] **Step 12: Red.** Rd `$L/T15_3_red.log` `test_dsl.jl`. Expected: 3 failures in
  "multi-atom metabolites" (`B[N, P]` builds; `B[]` fails only when the reaction is built;
  a bare `B` gets the terser message).
- [ ] **Step 13: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary f8d0199 -- src/dsl.jl | git apply
```

- [ ] **Step 14: Green.** G `$L/T15_3_green.log` `test_dsl.jl`: 238.
- [ ] **Step 15: Commit.** K `f8d0199` `src/dsl.jl test/test_dsl.jl` ("Take one chemical
  formula per atom bracket in @enzyme_reaction").

**Unit 4: 6c85d75** (an empty `oligomeric_state:` tuple; S2 review m3, its own commit)

- [ ] **Step 16: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 6c85d75 -- test/test_dsl.jl | git apply
```

- [ ] **Step 17: Red.** Rd `$L/T15_4_red.log` `test_dsl.jl`. Expected: 2 failures in
  "@enzyme_reaction sets the catalytic multiplicities once" (the empty tuple fails at run
  time with the constructor's message; the two-values message is the old wording).
- [ ] **Step 18: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 6c85d75 -- src/dsl.jl | git apply
```

- [ ] **Step 19: Green.** G `$L/T15_4_green.log` `test_dsl.jl`: 238.
- [ ] **Step 20: Commit.** K `6c85d75` `src/dsl.jl test/test_dsl.jl` ("Reject an empty
  oligomeric_state: tuple while the macro expands").

**Task end**

- [ ] **Step 21: Full suite.** F 15: 38,854.
- [ ] **Step 22: Census.** C 15 `base6`: zero differences.

### Task 16: Messages in the modifier grammar [S2c]

**Replays:** `[59136b3]`, `[c457930]`, `[8fc17da + 4a22318]`, `[257e81b]`, `[be61164]`.

**Files:**
- Modify: `src/types.jl` (the `Step` role-change error and docstring,
  `_assert_distinct_form_names`' form text), `src/dsl.jl` (`_peel_tag!`'s docstring, the
  missing-tag error, `_bare_symbols_from_values`)
- Test: `test/test_types.jl`, `test/test_dsl.jl`

**Interfaces:**
- Consumes: `_peel_tag!`, `_raw_met` (Task 6); `_assert_distinct_form_names` (Task 4).
- Produces: the `Step` role-change error advises "write a modifier copy as X::Mod on both
  sides"; a reactant's modifier copy prints `P::Mod` in the form-name error; the
  missing-tag error names "`:: EqualAI`, `:: NonequalAI` or `:: OnlyA`"; a list label's
  value other than a bare name gets "Atom brackets belong in @enzyme_reaction." only for
  `X[...]` and "List `X` bare; a role tag goes on the metabolite in a step." only for
  `X::Tag`.

**Gate:** F 16 = 38,858 (c457930 +1, be61164 +3); C 16 base6: zero differences.

**Unit 1: 59136b3** (the role-change error)

- [ ] **Step 1: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 59136b3 -- test/test_types.jl | git apply
```

- [ ] **Step 2: Red.** Rd `$L/T16_1_red.log` `test_types.jl`. Expected: 7 failures in
  "a binding keeps each metabolite's role" (the message still advises `X::Inh`).
- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 59136b3 -- src/types.jl | git apply
```

- [ ] **Step 4: Green.** G `$L/T16_1_green.log` `test_types.jl`: 618.
- [ ] **Step 5: Commit.** K `59136b3` `src/types.jl test/test_types.jl` ("Advise `X::Mod`
  in the Step role-change error").

**Unit 2: c457930** (the form-name error)

- [ ] **Step 6: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary c457930 -- test/test_types.jl | git apply
```

- [ ] **Step 7: Red.** Rd `$L/T16_2_red.log` `test_types.jl`. Expected: 1 failure in
  "two forms with one internal name are rejected" (the error writes `E(P::Inh)`).
- [ ] **Step 8: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary c457930 -- src/types.jl | git apply
```

- [ ] **Step 9: Green.** G `$L/T16_2_green.log` `test_types.jl`: 619.
- [ ] **Step 10: Commit.** K `c457930` `src/types.jl test/test_types.jl` ("Write a
  reactant's modifier copy as `X::Mod` in the form-name error").

**Unit 3: 8fc17da + 4a22318** (`_peel_tag!`'s docstring; docstring only)

- [ ] **Step 11: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 8fc17da -- src/dsl.jl | git apply
git show --format= --binary 4a22318 -- src/dsl.jl | git apply
```

- [ ] **Step 12: Green.** G `$L/T16_3_green.log` `test_dsl.jl`: 238.
- [ ] **Step 13: Commit.** K `8fc17da` `src/dsl.jl` ("Say why `_peel_tag!` reads `::Inh` as
  a role tag"; 4a22318 rewords the same docstring, which the message's "the docstring now
  says so" still covers).

**Unit 4: 257e81b** (the missing-tag error; S2 review m5)

- [ ] **Step 14: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 257e81b -- test/test_dsl.jl | git apply
```

- [ ] **Step 15: Red.** Rd `$L/T16_4_red.log` `test_dsl.jl`. Expected: 3 failures, one
  each in "@allosteric_mechanism (parsing & validation)" (the group case),
  "::Mod is never a step tag" and "allosteric steps carry one tag per statement, plain
  steps none" (the message names the tags in Symbol syntax).
- [ ] **Step 16: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 257e81b -- src/dsl.jl | git apply
```

- [ ] **Step 17: Green.** G `$L/T16_4_green.log` `test_dsl.jl`: 238.
- [ ] **Step 18: Commit.** K `257e81b` `src/dsl.jl test/test_dsl.jl` ("Name the step tags as
  the grammar writes them in the missing-tag error").

**Unit 5: be61164** (the list-label hints; S2 review m6)

- [ ] **Step 19: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary be61164 -- test/test_dsl.jl | git apply
```

- [ ] **Step 20: Red.** Rd `$L/T16_5_red.log` `test_dsl.jl`. Expected: 3 failures in "a
  mechanism macro's list labels take bare names" (every value gets the bracket hint).
- [ ] **Step 21: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary be61164 -- src/dsl.jl | git apply
```

- [ ] **Step 22: Green.** G `$L/T16_5_green.log` `test_dsl.jl test_types.jl`: 241, 619.
- [ ] **Step 23: Commit.** K `be61164` `src/dsl.jl test/test_dsl.jl` ("Hint where a bracket
  or a role tag goes only for that value").

**Task end**

- [ ] **Step 24: Full suite.** F 16: 38,858.
- [ ] **Step 25: Census.** C 16 `base6`: zero differences.

### Task 17: Print a reaction as its block; summary modifiers; the documented grammar [S2c]

**Replays:** `[2493151 + a421d95]`, `[67dc293 + 9838d79 + 5cc7f8e]`, `[3c44dd3]`,
`[7adea3e]`, `[b3431bf]`, `[60bd37a]`.

**Files:**
- Modify: `src/types.jl` (`show(::EnzymeMechanism)`, `show(::EnzymeReaction)`, the
  `EnzymeReaction` constructor's regulator sort, its docstring and doctest), `src/dsl.jl`
  (the three macro docstrings with jldoctests, the ABOUTME),
  `src/mechanism_enumeration.jl` (the no-catalytic-cycle error)
- Test: `test/test_types.jl`, `test/test_dsl.jl`, `test/test_mechanism_enumeration.jl`,
  `test/test_rate_eq_derivation.jl`, `test/mechanism_definitions_for_test_enzyme_derivation.jl`

**Interfaces:**
- Consumes: `EnzymeReaction`, `RegulatorMults`, `_forms`, `bound` (types.jl).
- Produces: `show(io::IO, rxn::EnzymeReaction)` prints the `@enzyme_reaction begin ... end`
  block (one line per non-empty label, 4-space indent; a regulator prints its
  multiplicities only when they differ from the catalytic ones); regulators sort by
  `(name, kind)` with a competitive inhibitor first; `show(::EnzymeMechanism)` ends
  `" | modifiers: <names>"`, the names bound as modifiers in the mechanism's forms, each
  once, sorted; the no-cycle error reads "no catalytic cycle for S ⇌ P + Q + R + T: ...";
  `$P/replay/docenv`, a docs environment that develops `$W` (Tasks 61 and 63 may reuse it
  after `Pkg.resolve()` in it, since S4 adds the SHA and Unicode deps).

**Gate:** F 17 = 38,872; C 17 base6: zero differences; the docs build passes; checkpoint S2c.

**Unit 1: 2493151 + a421d95** (the summary's `modifiers:`; S2 review m2 folded)

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2493151 -- test/test_types.jl | git apply
git show --format= --binary a421d95 -- test/test_types.jl | git apply
```

- [ ] **Step 2: Red.** Rd `$L/T17_1_red.log` `test_types.jl`. Expected: 2 failures in
  "Pretty printing": the `"| modifiers: I"` pin and the dual-role `endswith(...,
  " | modifiers: I")` (the summary still ends "| regulators: I, I"). The prototype's red
  was "1 + 3", since a421d95's red ran on 2493151's label: here its two
  `!contains(..., "modifiers")` assertions pass on the old "regulators" label.
- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2493151 -- src/types.jl | git apply
git show --format= --binary a421d95 -- src/types.jl | git apply
```

- [ ] **Step 4: Green.** G `$L/T17_1_green.log` `test_types.jl`: 622.
- [ ] **Step 5: Commit** (message edited: the fold's list)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M src/types.jl, M test/test_types.jl
git add src/types.jl test/test_types.jl
cat > $P/replay/msg <<'EOF'
List a mechanism's modifiers as `modifiers:` in its summary

`show` on an `EnzymeMechanism` ended with "| regulators: I", the label of the plain
macro that `modifiers:` replaced. It now ends with "| modifiers: I" and lists the names
the steps bind as modifiers in the mechanism's forms, each once, in name order. The
regulators of the reaction the Sig keeps would print "| modifiers: S" for a plain
mechanism whose substrate S the reaction also declares as an allosteric regulator or a
competitive inhibitor, though S binds only as itself, and "I, I" for a name in both
regulator roles. The pin and the dual-role assertion failed on the old text first.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

**Unit 2: 67dc293 + 9838d79 + 5cc7f8e** (the printed reaction; S2 review B1 and m1 folded)

- [ ] **Step 6: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 67dc293 -- test/test_dsl.jl test/test_types.jl | git apply
git show --format= --binary 9838d79 -- test/test_dsl.jl | git apply
git show --format= --binary 5cc7f8e -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 7: Red** (background). Rd `$L/T17_2_red.log`
  `test_dsl.jl test_types.jl test_mechanism_enumeration.jl`. Expected:
  - test_dsl.jl: "a reaction prints as the @enzyme_reaction text that builds it", 2 fail
    (the two pinned texts) and 4 errors (each of the four reactions' summary text fails to
    evaluate; the prototype's red had three reactions, 9838d79 adds the fourth);
    "@enzyme_reaction dual-role regulator (allosteric + competitive)", 2 fail (the
    allosteric-first reaction is not `==` and hashes apart).
  - test_types.jl: "Pretty printing", 2 fail (the two reaction pins).
  - test_mechanism_enumeration.jl: "no admissible catalytic cycle throws", 1 fail (the
    message holds "EnzymeReaction: S ⇌ P + Q + R + T").
- [ ] **Step 8: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 67dc293 -- src/types.jl | git apply
git show --format= --binary 9838d79 -- src/types.jl | git apply
git show --format= --binary 5cc7f8e -- src/mechanism_enumeration.jl | git apply
```

- [ ] **Step 9: Green** (background). G `$L/T17_2_green.log`
  `test_dsl.jl test_types.jl test_mechanism_enumeration.jl`: 249, 622, 32,534.
- [ ] **Step 10: Commit** (message edited: one paragraph each for 9838d79 and 5cc7f8e)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M src/mechanism_enumeration.jl, M src/types.jl, M test/test_dsl.jl,
                     # M test/test_mechanism_enumeration.jl, M test/test_types.jl
git add src/mechanism_enumeration.jl src/types.jl test/test_dsl.jl \
    test/test_mechanism_enumeration.jl test/test_types.jl
cat > $P/replay/msg <<'EOF'
Print a reaction as the @enzyme_reaction block that builds it

`show(::EnzymeReaction)` printed a summary, "EnzymeReaction: S ⇌ P | regulators: I",
that no macro reads and that named `allowed_catalytic_multiplicities:`, a label
@enzyme_reaction rejects. It now prints the block: each reactant with its formula,
the competitive inhibitors and allosteric regulators with their types, the catalytic
multiplicities unless they are the default 1, and the shared catalytic sites. A
regulator prints its multiplicities only when they differ from the catalytic ones,
which a bare name takes.

The constructor sorts regulators by name and then kind, a competitive inhibitor first,
the order `show` prints. A stable sort by name alone kept a name declared as an
allosteric regulator before its competitive-inhibitor line in that order, so its text
built a reaction that was not == to it, and the two declaration orders built unequal
reactions with different hashes and Sigs. Only a name in both roles moves.

The no-catalytic-cycle error writes the reaction as its substrates and products,
`S ⇌ P + Q + R + T`, instead of interpolating the multi-line block into its sentence.

A new test evaluates the printed text of four reactions, which cover the shapes the
tests use, the allosteric-first dual-role reaction among them, and checks that each
builds an equal reaction; it pins two texts. The dual-role testset checks that both
label orders build one reaction, and the no-cycle test pins the whole message. The two
summary pins and EnzymeReaction's doctest now show the block.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

**Unit 3: 3c44dd3** (the macros' docstrings; docs only)

- [ ] **Step 11: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 3c44dd3 -- src/dsl.jl | git apply
```

- [ ] **Step 12: Green.** G `$L/T17_3_green.log` `test_dsl.jl`: 249.
- [ ] **Step 13: Docs environment** (once; skip if `$P/replay/docenv/Manifest.toml` exists).
  Background, then `wait_log T17_docenv`:

The environment starts from S5's docs Manifest (`$P/s5/docenv`, resolved for the
prototype's docs builds; `$W/docs/Project.toml` equals its Project.toml), and Pkg runs
offline and keeps every other version, so the build runs on the Documenter and
DocumenterCitations versions behind S5's warning counts, as Task 0 Step 5 does for the
worktree.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
mkdir -p $P/replay/docenv
cp $W/docs/Project.toml $P/s5/docenv/Manifest.toml $P/replay/docenv/
$P/jl.sh $L/T17_docenv.log 1800 --project=$P/replay/docenv -e "using Pkg; Pkg.offline(true); Pkg.develop(path = \"$W\"; preserve = Pkg.PRESERVE_ALL); Pkg.instantiate()"
grep -A5 '^\[\[deps.Documenter' $P/replay/docenv/Manifest.toml | grep -E '^\[\[|^version'
```

  Expected: `exit 0`; `$P/replay/docenv/Manifest.toml` lists `EnzymeRates` with
  `path = "$W"`; the `grep` prints `Documenter` 1.19.0 and `DocumenterCitations` 1.5.0,
  S5's; afterwards `git status --short` lists only ` M src/dsl.jl` (Pkg writes nothing
  into `$W`).
- [ ] **Step 14: Docstring doctests.** Background, then `wait_log T17_3_doctest`:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh $L/T17_3_doctest.log 1800 --project=$P/replay/docenv -e 'using Documenter, EnzymeRates; doctest(EnzymeRates; manual = false)'
```

  Expected: `Doctests: EnzymeRates` passes, `exit 0` (S2c checked the same run fails when
  an expected output changes).
- [ ] **Step 15: Commit.** K `3c44dd3` `src/dsl.jl` ("Document the macros' grammar, each
  with a checked example").

**Unit 4: 7adea3e** (repeated list labels; test-only, green on arrival; S2 review m4)

- [ ] **Step 16: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 7adea3e -- test/test_dsl.jl | git apply
```

- [ ] **Step 17: Green.** G `$L/T17_4_green.log` `test_dsl.jl`: 252.
- [ ] **Step 18: Mutation check** (the test can fail: every list label made to replace the
  earlier line's entries). Mutate `src/dsl.jl`:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/s2c_r1/mutate_m4.py $W
git status --short
```

  Expected: `mutated 7 sites`; ` M src/dsl.jl` and ` M test/test_dsl.jl`. Then N
  `$L/T17_4_mut.log` `test/test_dsl.jl` `"list labels may repeat; their names accumulate"`
  (background, `wait_log T17_4_mut`). Expected under `==== run_named summary ====`:
  `test/test_dsl.jl | list labels may repeat; their names accumulate | FAIL: ...` (the
  prototype's run: two expansion errors, A and X undeclared, and one false comparison).
  Then revert the mutation, whatever the run printed:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git checkout -- src/dsl.jl
git diff --stat src; git status --short
```

  Expected: no diff output; ` M test/test_dsl.jl` only.
- [ ] **Step 19: Commit.** K `7adea3e` `test/test_dsl.jl` ("Test that each macro's list
  labels accumulate over repeated lines").

**Unit 5: b3431bf** (a fixture comment; test-only)

- [ ] **Step 20: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary b3431bf -- test/mechanism_definitions_for_test_enzyme_derivation.jl | git apply
```

- [ ] **Step 21: Green.** G `$L/T17_5_green.log` `test_types.jl`: 622 (the run loads the
  fixtures file).
- [ ] **Step 22: Commit.** K `b3431bf`
  `test/mechanism_definitions_for_test_enzyme_derivation.jl` ("Name the MWC dimer fixture's
  forms as its block writes them").

**Unit 6: 60bd37a** ("modifier copy" in titles, comments and locals; test-only)

- [ ] **Step 23: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 60bd37a -- test/test_dsl.jl test/test_mechanism_enumeration.jl test/test_rate_eq_derivation.jl test/test_types.jl | git apply
```

- [ ] **Step 24: Green** (background). G `$L/T17_6_green.log`
  `test_dsl.jl test_types.jl test_rate_eq_derivation.jl test_mechanism_enumeration.jl`:
  252, 622, 3,184, 32,534.
- [ ] **Step 25: Commit.** K `60bd37a` `test/test_dsl.jl test/test_mechanism_enumeration.jl
  test/test_rate_eq_derivation.jl test/test_types.jl` ("Say "modifier copy" in test titles,
  comments and locals").

**Task end**

- [ ] **Step 26: Full suite.** F 17: 38,872.
- [ ] **Step 27: Census.** C 17 `base6`: zero differences.
- [ ] **Step 28: Docs build.** Background, then `wait_log T17_docs`:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh $L/T17_docs.log 1800 --project=$P/replay/docenv --heap-size-hint=2500M docs/make.jl; rm -rf $W/docs/build
```

  Expected: `exit 0`; doctests and `checkdocs = :exports` pass. At most two warnings, both
  in S2c's build: an `@example` block on identify/tutorial.md renders a large `text/html`
  output, and deploydocs cannot auto-detect the building environment and skips deployment.
- [ ] **Step 29: Checkpoint S2c** (tip 0a13357). Run

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git diff --stat 0a13357 HEAD -- src test docs ROADMAP.md Project.toml .claude
```

  Expected: exactly these 12 files, with the line counts of the dry replay
  (`$S/p3skel/sim_final.log`): Project.toml 2, ROADMAP.md 20,
  docs/src/deriving/ping_pong.md 2, docs/src/deriving/re_vs_ss.md 6,
  docs/src/deriving/regulators.md 6, docs/src/developer.md 6,
  docs/src/identify/enumeration_engine.md 5, docs/src/identify/tutorial.md 7,
  docs/superpowers/plans/2026-10-16-graph-mechanism.md (the plan, D2),
  docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md 48,
  test/mechanism_definitions_for_test_enzyme_derivation.jl 35,
  test/test_rate_eq_derivation.jl 31. No `src/` file; getting_started.md is absent
  (the tip holds 0a13357's paragraph). Then check that each holds only its named changes:
  build the expected tree in a scratch index and compare it with HEAD.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
export GIT_INDEX_FILE=$P/replay/ckpt.idx
git read-tree 0a13357
git diff 4f9a34a 7205ba4 -- docs/superpowers/specs ROADMAP.md | git apply --cached   # D1, D3
for c in 5c08ede 4a96d9a 0b5e637 68f46f0 33b94d8 cfe9211; do git show --format= --binary $c | git apply --cached; done   # D4
for c in 2b29fd2 be65542 0b2021d a84d194; do git show --format= --binary $c | git apply --cached; done   # docs ahead
git update-index --add --cacheinfo 100644,$(git rev-parse HEAD:docs/superpowers/plans/2026-10-16-graph-mechanism.md),docs/superpowers/plans/2026-10-16-graph-mechanism.md   # D2
git diff --cached --stat HEAD -- src test docs ROADMAP.md Project.toml .claude
rm -f $GIT_INDEX_FILE; unset GIT_INDEX_FILE
```

  Expected: every `git apply --cached` silent, and the last `git diff` prints nothing.
  Departures in S2c: none in code; the messages of Task 15 unit 2 and Task 17 units 1-2 are
  edited as written there.

### Task 18: A mechanism's identity and names without its reaction [S3a]

**Replays:** `[b588424]`, `[db5b3ff]`, `[10527e7]`.

**Files:**
- Modify: `src/types.jl` (the `==`/`hash` loop, the name accessors, the compiled
  `metabolites` generator, `show(::EnzymeMechanism)`)
- Test: `test/test_mechanism_enumeration.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: `reaction(m)`, `steps(m)`, `_forms`, `bound` (types.jl).
- Produces: `==` and `hash` of `Mechanism`/`AllostericMechanism` skip the `reaction` field;
  `substrates(m::_AnyMech)`, `products(m::_AnyMech)` (sorted `Vector{Symbol}`),
  `modifiers(m::_AnyMech)` (every name bound as a `CompetitiveInhibitor` in any form,
  sorted), `metabolites(m::Mechanism)` (substrates, products, modifiers, each once),
  `metabolites(am::AllostericMechanism)` (then site ligands, site by site);
  `metabolites(::EnzymeMechanism)` reads `metabolites(Mechanism(...))`.

**Gate:** F 18 = 38,886; C 18 base6: zero differences.

The S3a red and green logs were lost in the outage; the reds are S3_report.md S3a
section 2's.

**Unit 1: b588424** (the init golden's key reads `steps(m)`; test-only)

- [ ] **Step 1: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary b588424 -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 2: Green** (background). G `$L/T18_1_green.log`
  `test_mechanism_enumeration.jl`: 32,534.
- [ ] **Step 3: Commit.** K `b588424` `test/test_mechanism_enumeration.jl` ("Read the init
  golden's groups through steps(m)").

**Unit 2: db5b3ff** (identity without the reaction)

- [ ] **Step 4: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary db5b3ff -- test/test_types.jl | git apply
```

- [ ] **Step 5: Red.** Rd `$L/T18_2_red.log` `test_types.jl`. Expected: 2 failures (of 4)
  in "a mechanism's identity leaves out its reaction" (the two equalities; the two
  reaction inequalities pass).
- [ ] **Step 6: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary db5b3ff -- src/types.jl | git apply
```

- [ ] **Step 7: Green.** G `$L/T18_2_green.log` `test_types.jl`: 626.
- [ ] **Step 8: Commit.** K `db5b3ff` `src/types.jl test/test_types.jl` ("Compare mechanisms
  by their steps, not their reactions").

**Unit 3: 10527e7** (name accessors)

- [ ] **Step 9: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 10527e7 -- test/test_types.jl | git apply
```

- [ ] **Step 10: Red.** Rd `$L/T18_3_red.log` `test_types.jl`. Expected: 9 errors in "a
  mechanism names its substrates, products and modifiers" (no method
  `substrates(::Mechanism)`, `metabolites(::AllostericMechanism)`, ...).
- [ ] **Step 11: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 10527e7 -- src/types.jl | git apply
```

- [ ] **Step 12: Green.** G `$L/T18_3_green.log` `test_types.jl`: 636.
- [ ] **Step 13: Commit.** K `10527e7` `src/types.jl test/test_types.jl` ("Name a
  mechanism's substrates, products and modifiers").

**Task end**

- [ ] **Step 14: Full suite.** F 18: 38,886.
- [ ] **Step 15: Census.** C 18 `base6`: zero differences.

### Task 19: Labels and metabolite names from the mechanism [S3a]

**Replays:** `[31ff003]`, `[6ba4fed]`. Both refactors: green before (Task 18's full suite)
and after, 38,886.

**Files:**
- Modify: `src/types.jl`, `src/thermo.jl`, `src/rate_eq_derivation.jl`,
  `src/mechanism_enumeration.jl`, `src/identify_rate_equation.jl`
- Test: `test/test_joint_solve.jl`, `test/test_mechanism_enumeration.jl` (calls follow the
  moved signatures)

**Interfaces:**
- Consumes: Task 18's `substrates(m)`, `products(m)`, `modifiers(m)`, `metabolites(m)`.
- Produces: `_is_reactant(m::_AnyMech, n::Symbol)`; the label renderers take `ctx`, a
  mechanism or the reaction of one under construction: `_is_reactant_copy(x, ctx)`,
  `_met_label(x, ctx)`, `_form_label(sp, ctx)`, `_side_label`, `_forward_sides(s, ctx)`,
  `_step_text(s, ctx)`; name readers on `m`: `_fundamental_cycles(m::_AnyMech, groups)`,
  `_joint_system(m::_AnyMech)`, `_raw_symbolic_rate_polys(m::_AnyMech, groups, names)`,
  `_compute_alpha(m, flat, conc_set, species, idx, segments, step_to_K)`,
  `_flux_carrying_steps(groups, m)`, `_flux_carrying_groups(groups, m)`,
  `_reactant_signs(m)`, `_re_turnover_cycle(groups, m, ...)`, `_has_vmax(groups, m, side)`,
  `_chemistry_equilibrates_both_sides(groups, m)`, `_seed_candidate_screen(groups,
  m::Mechanism)`, `_minimal_flips(admissible, base, m, ...)`; atoms from the caller:
  `_assert_atom_conserving(m, rxn::EnzymeReaction)`, `_seed_variants(m::Mechanism,
  rxn::EnzymeReaction)`; `_concentration_symbols` deleted.

**Gate:** F 19 = 38,886; C 19 base6: zero differences.

**Unit 1: 31ff003** (renderers take the mechanism; source only)

- [ ] **Step 1: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 31ff003 -- src/mechanism_enumeration.jl src/types.jl | git apply
```

- [ ] **Step 2: Green** (background). G `$L/T19_1_green.log`
  `test_types.jl test_joint_solve.jl test_mechanism_enumeration.jl`: 636, 453, 32,534.
- [ ] **Step 3: Commit.** K `31ff003` `src/mechanism_enumeration.jl src/types.jl` ("Render
  constant and form labels from the mechanism").

**Unit 2: 6ba4fed** (name readers take `m`). No red run: the test edits only follow the
moved helpers' signatures, so before the source they error on every call.

- [ ] **Step 4: Apply the test paths, then the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 6ba4fed -- test/test_joint_solve.jl test/test_mechanism_enumeration.jl | git apply
git show --format= --binary 6ba4fed -- src/identify_rate_equation.jl src/mechanism_enumeration.jl src/rate_eq_derivation.jl src/thermo.jl | git apply
```

- [ ] **Step 5: Green** (background). G `$L/T19_2_green.log`
  `test_joint_solve.jl test_mechanism_enumeration.jl`: 453, 32,534.
- [ ] **Step 6: Commit.** K `6ba4fed` `src/identify_rate_equation.jl
  src/mechanism_enumeration.jl src/rate_eq_derivation.jl src/thermo.jl
  test/test_joint_solve.jl test/test_mechanism_enumeration.jl` ("Read metabolite names from
  the mechanism in the derivation and moves").

**Task end**

- [ ] **Step 7: Full suite.** F 19: 38,886.
- [ ] **Step 8: Census.** C 19 `base6`: zero differences.

### Task 20: The boundary check [S3a]

**Replays:** `[2b79bd3 + e4b5927 + 216be7f: docs/src/deriving/mwc_allostery.md]`, one
commit.

**Files:**
- Create: `src/admission.jl`
- Modify: `src/EnzymeRates.jl` (includes admission.jl after thermo.jl), `src/types.jl`,
  `src/dsl.jl` (`_compiled_mechanism` wraps both constructors in `_well_posed`),
  `src/thermo.jl`, `src/rate_eq_derivation.jl`, `src/mechanism_enumeration.jl`,
  `src/identify_rate_equation.jl`, `docs/src/deriving/mwc_allostery.md` (lines 93 and 216)
- Test: `test/allosteric_ground_truth.jl` (comments), `test/test_identify_rate_equation.jl`,
  `test/test_mechanism_enumeration.jl`, `test/test_oracles.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: Task 19's name readers; `_bottomless_re_segment`, `_inactive_groups`.
- Produces (src/admission.jl): `_check_well_posed(m::Mechanism)`,
  `_check_well_posed(am::AllostericMechanism)` (bottomless on `steps(am)`, one flip
  segment, the `:OnlyA` check, bottomless on `_inactive_groups(am)`),
  `_assert_bottomed_re_segments(m::_AnyMech, groups::Vector{Vector{Step}})`,
  `_assert_one_flip_segment(am::AllostericMechanism)`, `_well_posed(m) = (_check_well_posed(m); m)`,
  `_derivation_input(M::Type)`; called by both macros (`_compiled_mechanism`),
  `init_mechanisms` and `seed_mechanisms` (on what they return), `expand_mechanisms` (each
  child), identify's `_base_tier` (each flip child) and the derivation entries (through the
  `Type` method of `_dependent_param_exprs`, `all_parameters`, `_ss_rate_constant_names`).
  thermo.jl: `_onlya_haldane_violation(m::_AnyMech, cat_steps, cat_allo_states) ->
  Vector{Int}` (offending groups, empty when none),
  `_onlya_violation_message(am::AllostericMechanism, groups::Vector{Int})`. types.jl:
  `_bottomless_re_segment(m, steps)`; the constructors check representation only.

**Gate:** F 20 = 38,902; C 20 base6: zero differences.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2b79bd3 -- test/allosteric_ground_truth.jl test/test_identify_rate_equation.jl test/test_mechanism_enumeration.jl test/test_oracles.jl test/test_types.jl | git apply
git show --format= --binary e4b5927 -- test/test_types.jl | git apply
```

- [ ] **Step 2: Red** (background; named run, S3a section 9 task 6). The rewritten tests
  that call `_check_well_posed` or the new `_onlya_haldane_violation` signature error in a
  whole-file run before the source; the four new tests carry the red.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh $L/T20_red.log 1800 --project=. --heap-size-hint=2500M $P/s3c/scripts/run_named.jl \
  test/test_types.jl "the inactive conformation's segments need a bottom form too" \
  test/test_types.jl "every derivation entry checks that the rate law is defined" \
  test/test_mechanism_enumeration.jl "expand_mechanisms checks that every child's rate law is defined" \
  test/test_identify_rate_equation.jl "_base_tier records a flip child whose rate law is not defined"
```

  Then `wait_log T20_red` until it ends. Expected under `==== run_named summary ====`,
  `exit 0`: the inactive-segment test `FAIL` (1 failed: no exception, the macro
  does not yet check the inactive conformation); the derivation-entry test `FAIL` (an
  error: the constructor throws the bottomless error while `Mechanism(raw)` builds the
  fixture); the `expand_mechanisms` test `FAIL` (an error at fixture construction, for the
  same reason); the `_base_tier` test `PASS` (the constructor already threw inside the
  flip).
- [ ] **Step 3: Apply the source and the docs**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2b79bd3 -- src/EnzymeRates.jl src/admission.jl src/dsl.jl src/identify_rate_equation.jl src/mechanism_enumeration.jl src/rate_eq_derivation.jl src/thermo.jl src/types.jl | git apply
git show --format= --binary 216be7f -- docs/src/deriving/mwc_allostery.md | git apply
git status --short
```

  Expected: `?? src/admission.jl` beside the 13 modified paths the commit step adds.
- [ ] **Step 4: Green** (background). G `$L/T20_green.log`
  `test_types.jl test_oracles.jl test_mechanism_enumeration.jl test_identify_rate_equation.jl`:
  646, 634, 32,536, 446.
- [ ] **Step 5: Commit** (message edited: one sentence each for e4b5927 and the docs)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/EnzymeRates.jl src/admission.jl src/dsl.jl src/identify_rate_equation.jl \
    src/mechanism_enumeration.jl src/rate_eq_derivation.jl src/thermo.jl src/types.jl \
    docs/src/deriving/mwc_allostery.md test/allosteric_ground_truth.jl \
    test/test_identify_rate_equation.jl test/test_mechanism_enumeration.jl \
    test/test_oracles.jl test/test_types.jl
git status --short   # every line staged: `M ` or `A ` (src/admission.jl), no ` M` or `??`
cat > $P/replay/msg <<'EOF'
Check a mechanism's rate law at the DSL, the enumerator's exits and derivation

The constructors rejected a bottomless rapid-equilibrium segment, free forms in two
segments and an unsatisfiable :OnlyA cycle, so every internal construction ran checks
meant for the mechanisms users and the enumerator hand out. src/admission.jl now holds
the boundary check, _check_well_posed, which the mechanism macros, init_mechanisms,
seed_mechanisms, expand_mechanisms, identify's base tier and every derivation entry
run; the constructors check representation only. The check also covers the inactive
conformation's graph, which only the derivation read before. The :OnlyA test returns
the offending groups, so the moves read its verdict without naming constants, and the
check names them from the mechanism itself instead of building a plain one.

The bottomless-segment testset's title names the macro, whose boundary check rejects
the mechanism the constructor builds. The MWC page says that @allosteric_mechanism, not
the constructor, rejects a second free form outside the first's segment and the
inconsistent :OnlyA combination.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

- [ ] **Step 6: Full suite.** F 20: 38,902 (Aqua and JET included: the two moved helpers are
  typed in admission.jl).
- [ ] **Step 7: Census.** C 20 `base6`: zero differences.

### Task 21: Naming without a cache; Symbol ligands [S3a]

**Replays:** `[f7bd003]`, `[4d13802]`.

**Files:**
- Modify: `src/types.jl`, `src/thermo.jl`, `src/rate_eq_derivation.jl`,
  `src/mechanism_enumeration.jl`
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`, `test/test_accessors.jl`,
  `test/test_dsl.jl`, `test/test_identify_rate_equation.jl`,
  `test/test_mechanism_enumeration.jl`, `test/test_rate_eq_derivation.jl`,
  `test/test_types.jl`

**Interfaces:**
- Consumes: `name(p, m)` (the chokepoint), `RegulatorySite`, `Kreg`.
- Produces: no `naming` field and no `_NamingCache`; `_class_sides(m)`,
  `_constant_sides(step, slot, state, m, sides)`, `name(p::Union{Kfor,Krev,Krapid}, m,
  sides = _class_sides(m))`, `_step_names(m, state, sides = _class_sides(m))`;
  `RegulatorySite(ligands::Vector{Symbol}, multiplicity, allo_states)` (ligands sorted),
  `Kreg(ligand::Symbol, state::Symbol)`, `allosteric_regulators(am) -> Vector{Symbol}`,
  `_reg_K(am::AllostericMechanism, lig::Symbol, tag::Symbol, inactive::Bool)`.

**Gate:** F 21 = 38,902; C 21 base6: zero differences, and
`compare.jl census/S3a_6 census/replay_T21 census/rename_phase3.toml`: zero differences;
K15: test_rate_eq_derivation.jl and test_mechanism_enumeration.jl within +20% of their
time before f7bd003; checkpoint S3a.

**Unit 1: f7bd003** (naming without a cache; refactor)

- [ ] **Step 1: K15 baseline** (background; before any apply)
  G `$L/T21_1_k15_before.log` `test_rate_eq_derivation.jl test_mechanism_enumeration.jl`:
  3,184, 32,536. Note each file's Time column.
- [ ] **Step 2: Apply the test path, then the source** (the test edits are a title and two
  comments)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary f7bd003 -- test/test_types.jl | git apply
git show --format= --binary f7bd003 -- src/rate_eq_derivation.jl src/thermo.jl src/types.jl | git apply
```

- [ ] **Step 3: Green and K15** (background). G `$L/T21_1_green.log`
  `test_rate_eq_derivation.jl test_mechanism_enumeration.jl test_types.jl`: 3,184, 32,536,
  646. Each of the first two files' Time is at most 1.2 times Step 1's (the prototype:
  61.1 s to 61.1 s, 130.2 s to 133.0 s). Over 1.2: stop and ask Denis (K15: the helpers
  would then pass `sides`).
- [ ] **Step 4: Commit.** K `f7bd003` `src/rate_eq_derivation.jl src/thermo.jl src/types.jl
  test/test_types.jl` ("Name constants without a cache inside the mechanism").

**Unit 2: 4d13802** (Symbol ligands)

- [ ] **Step 5: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 4d13802 -- test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_accessors.jl test/test_dsl.jl test/test_identify_rate_equation.jl test/test_mechanism_enumeration.jl test/test_rate_eq_derivation.jl test/test_types.jl | git apply
```

- [ ] **Step 6: Red** (background; named run: every file's site readers now take names, so
  a whole-file run errors wherever a site is built)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh $L/T21_2_red.log 1800 --project=. --heap-size-hint=2500M $P/s3c/scripts/run_named.jl \
  test/test_types.jl "RegulatorySite: validation + accessors" \
  test/test_types.jl "Parameter family: step constants and Kreg"
```

  Then `wait_log T21_2_red` until it ends. Expected: both `FAIL` with an error (no method
  `RegulatorySite(::Vector{Symbol}, ...)`; no method `Kreg(::Symbol, ::Symbol)`).
- [ ] **Step 7: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 4d13802 -- src/mechanism_enumeration.jl src/rate_eq_derivation.jl src/types.jl | git apply
```

- [ ] **Step 8: Green** (background). G `$L/T21_2_green.log`
  `test_accessors.jl test_types.jl test_dsl.jl test_rate_eq_derivation.jl test_mechanism_enumeration.jl test_identify_rate_equation.jl`:
  155, 647, 251, 3,184, 32,536, 446.
- [ ] **Step 9: Commit.** K `4d13802` `src/mechanism_enumeration.jl
  src/rate_eq_derivation.jl src/types.jl test/mechanism_definitions_for_test_enzyme_derivation.jl
  test/test_accessors.jl test/test_dsl.jl test/test_identify_rate_equation.jl
  test/test_mechanism_enumeration.jl test/test_rate_eq_derivation.jl test/test_types.jl`
  ("Hold a regulatory site's ligands as names").

**Task end**

- [ ] **Step 10: Full suite.** F 21: 38,902.
- [ ] **Step 11: Census.** C 21 `base6`: zero differences. Then compare the same corpus
  with the prototype's S3a corpus:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
(cd $P && ./jl.sh census/logs/cmp_replay_T21_S3a_6.log 600 --project=census/env --heap-size-hint=2500M census/compare.jl census/S3a_6 census/replay_T21 census/rename_phase3.toml)
census_check 21_S3a_6
```

  Expected: `RESULT: PASS` and nothing else (the replay equals the prototype's S3a corpus).
- [ ] **Step 12: Checkpoint S3a** (tip e4b5927). Run

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git diff --stat e4b5927 HEAD -- src test docs ROADMAP.md Project.toml .claude
```

  Expected: exactly these 13 files: Project.toml 2, ROADMAP.md 20,
  docs/src/deriving/mwc_allostery.md 8, docs/src/deriving/ping_pong.md 2,
  docs/src/deriving/re_vs_ss.md 6, docs/src/deriving/regulators.md 6,
  docs/src/developer.md 6, docs/src/identify/enumeration_engine.md 5,
  docs/src/identify/tutorial.md 7, the plan (D2), the spec 48 (D1), FX 35, TRD 31. Then
  the scratch-index check of Task 17 Step 29 against e4b5927,
  with 216be7f's mwc_allostery.md piece among the docs ahead (one block; the index file
  must be set in the same shell as every `git apply --cached`):

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
export GIT_INDEX_FILE=$P/replay/ckpt.idx
git read-tree e4b5927
git diff 4f9a34a 7205ba4 -- docs/superpowers/specs ROADMAP.md | git apply --cached   # D1, D3
for c in 5c08ede 4a96d9a 0b5e637 68f46f0 33b94d8 cfe9211; do git show --format= --binary $c | git apply --cached; done   # D4
for c in 2b29fd2 be65542 0b2021d a84d194; do git show --format= --binary $c | git apply --cached; done   # docs ahead
git show --format= --binary 216be7f -- docs/src/deriving/mwc_allostery.md | git apply --cached   # docs ahead
git update-index --add --cacheinfo 100644,$(git rev-parse HEAD:docs/superpowers/plans/2026-10-16-graph-mechanism.md),docs/superpowers/plans/2026-10-16-graph-mechanism.md   # D2
git diff --cached --stat HEAD -- src test docs ROADMAP.md Project.toml .claude
rm -f $GIT_INDEX_FILE; unset GIT_INDEX_FILE
```

  Expected: every `git apply --cached` silent, and the last `git diff` prints nothing.
  Departures in S3a: none in code; Task 20's message is edited as written there, and K
  adds the Claude-Session trailer to the other main commits' messages, which lack it.

### Task 22: Position-picked fixtures written out; populations sampled by content [S3b]

**Replays:** `[ad4aca2]`, `[cd380a2]`. Test-only: green before (Task 21) and after.

**Files:**
- Test: `test/test_accessors.jl`, `test/test_identify_rate_equation.jl`,
  `test/test_oracles.jl`, `test/test_compile_budget.jl`

**Interfaces:**
- Produces: TID's `_testhelper_uni_seed` (the uni-uni seed written out, shared by the
  identify tests); TO's `_testhelper_content_text(m)`, `_testhelper_population_one_round`,
  `_testhelper_population_sample` (120 mechanisms spaced through the population sorted by
  content text); TCB writes the fit census's five bi-bi seeds out.

**Gate:** F 22 = 38,905. No census (no src).

- [ ] **Step 1: Apply ad4aca2's test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary ad4aca2 -- test/test_accessors.jl test/test_identify_rate_equation.jl | git apply
```

- [ ] **Step 2: Green** (background). G `$L/T22_1_green.log`
  `test_accessors.jl test_identify_rate_equation.jl`: 155, 446.
- [ ] **Step 3: Commit.** K `ad4aca2` `test/test_accessors.jl
  test/test_identify_rate_equation.jl` ("Write out the test mechanisms picked as an
  enumeration's first").
- [ ] **Step 4: Apply cd380a2's test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary cd380a2 -- test/test_compile_budget.jl test/test_oracles.jl | git apply
```

- [ ] **Step 5: Green** (background; the budget file first, as in the suite). G
  `$L/T22_2_green.log` `test_compile_budget.jl test_oracles.jl`: 16, 637 (+3: three more
  same-model regroupings in the new sample; the gates hold: 5 Sig shapes, 5 name tuples,
  at least 30 allosteric, at least 20 and 50 regroupings).
- [ ] **Step 6: Commit.** K `cd380a2` `test/test_compile_budget.jl test/test_oracles.jl`
  ("Sample the oracle population and the fit census by content").
- [ ] **Step 7: Full suite.** F 22: 38,905.

### Task 23: Steps, groups and forms in content order (rows 1-2) [S3b]

**Replays:** `[1362b5b]`, `[9916684]`.

**Files:**
- Modify: `src/types.jl` (`_step_canonical_key`, `_forms`; `_met_names` deleted)
- Test: `test/test_joint_solve.jl`, `test/test_mechanism_enumeration.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: `_step_index_key`, `_form_index_key` (thermo.jl, phase 2).
- Produces: `_step_canonical_key(s::Step) = (_step_index_key(s)..., is_equilibrium(s))`;
  `_forms(groups)` lists the distinct forms sorted by `_form_index_key` (free E first).

**Gate:** F 23 = 38,907; C 23 `S3b_end`: `census_check 23` prints `RESULT: PASS` and one
more line, `  init: A 4885, B 4885; set equal; order changed` (ter-bi's init order, which
row 7 brings in Task 25; `$P/census/logs/rec_cmp_S3a_6_vs_S3b_end.log` shows it); the four
path laws of S3b section 5 and the order-only child counts already match.

The S3b red logs were lost in the outage; the reds are S3_report.md S3b section 2's, plus the
re-pins each commit carries.

**Unit 1: 1362b5b** (row 1)

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 1362b5b -- test/test_joint_solve.jl test/test_types.jl | git apply
```

- [ ] **Step 2: Red.** Rd `$L/T23_1_red.log` `test_types.jl test_joint_solve.jl`. Expected:
  test_types.jl, 2 failures in "sort key orders bindings and isomerizations" (string keys)
  and 1 in "Mechanism canonicalizes step order" (old flat order [1, 6, 2, 5, 4, 3]);
  test_joint_solve.jl, 1 failure in "a tied group turns with the lead the form order ranks
  first (C9)" (the re-pinned step list).
- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 1362b5b -- src/types.jl | git apply
```

- [ ] **Step 4: Green.** G `$L/T23_1_green.log` `test_types.jl test_joint_solve.jl`: 648,
  453.
- [ ] **Step 5: Commit.** K `1362b5b` `src/types.jl test/test_joint_solve.jl
  test/test_types.jl` ("Order steps and kinetic groups by content").

**Unit 2: 9916684** (row 2)

- [ ] **Step 6: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9916684 -- test/test_mechanism_enumeration.jl test/test_types.jl | git apply
```

- [ ] **Step 7: Red** (background). Rd `$L/T23_2_red.log`
  `test_types.jl test_mechanism_enumeration.jl`. Expected: test_types.jl, 1 failure in
  "Mechanism canonicalizes step order" (forms in first-seen order) and 1 in "the inactive
  conformation's segments need a bottom form too" (the re-pinned segment message);
  test_mechanism_enumeration.jl, 1 failure in "expand_mechanisms checks that every child's
  rate law is defined" (the re-pinned segment message).
- [ ] **Step 8: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9916684 -- src/types.jl | git apply
```

- [ ] **Step 9: Green** (background). G `$L/T23_2_green.log`
  `test_types.jl test_mechanism_enumeration.jl`: 649, 32,536.
- [ ] **Step 10: Commit.** K `9916684` `src/types.jl test/test_mechanism_enumeration.jl
  test/test_types.jl` ("List a mechanism's forms in the canonical form order").

**Task end**

- [ ] **Step 11: Full suite.** F 23: 38,907.
- [ ] **Step 12: Census.** C 23 `S3b_end`, as the Gate says.

### Task 24: Segment roots, cycle rows and form identity by content (rows 3-5) [S3b]

**Replays:** `[77b7b61 + 216be7f: docs/src/deriving/cha_king_altman.md]`, `[eb508a5]`,
`[9c6d4cb]`.

**Files:**
- Modify: `src/rate_eq_derivation.jl` (`_segment_root`, `_reachable_from_free`,
  `_inactive_groups`), `src/thermo.jl` (`_thermodynamic_constraints`), `src/admission.jl`
  (`_assert_one_flip_segment`), `src/mechanism_enumeration.jl` (`_forms_where_free`, the
  dead-end tables), `src/types.jl` (`Species` docstring), `docs/src/deriving/cha_king_altman.md`
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/test_rate_eq_derivation.jl`, `test/test_mechanism_enumeration.jl`,
  `test/test_types.jl`

**Interfaces:**
- Consumes: Task 23's `_forms`, `_form_index_key`.
- Produces: `_segment_root(segment, species)` breaks its last tie by `_form_index_key`;
  `_thermodynamic_constraints` keys form rows by `Species` in `_forms` order;
  `_reachable_from_free(groups)` and `_inactive_groups(am)` use `Set{Species}`;
  `_forms_where_free(m, kind, name) -> Set{Species}`; `_assert_one_flip_segment` lists free
  forms by `_form_label` in `_forms` order; FX's `_testhelper_flip_form(mech)` takes
  `argmin(_form_index_key, ...)`; TE's `_testhelper_form_labels(m, forms)`.

**Gate:** F 24 = 38,909; C 24 `S3b_end`: `census_check 24` prints `RESULT: PASS` and one
more line, `  init: A 4885, B 4885; set equal; order changed` (ter-bi's init order, which
row 7 brings in Task 25).

**Unit 1: 77b7b61 + cha_king_altman.md** (row 3)

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 77b7b61 -- test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_rate_eq_derivation.jl | git apply
```

- [ ] **Step 2: Red.** Rd `$L/T24_1_red.log` `test_rate_eq_derivation.jl`. Expected: 2
  failures in "a segment's root ties go to the canonical form order" (EAb picked). The
  fixtures' flip form already follows the key order; no spec's free forms order apart by
  name and key, so no oracle check moves.
- [ ] **Step 3: Apply the source and the docs**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 77b7b61 -- src/rate_eq_derivation.jl | git apply
git show --format= --binary 216be7f -- docs/src/deriving/cha_king_altman.md | git apply
```

- [ ] **Step 4: Green.** G `$L/T24_1_green.log` `test_rate_eq_derivation.jl`: 3,186.
- [ ] **Step 5: Commit** (message edited: one sentence for the docs)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M docs/src/deriving/cha_king_altman.md, M src/rate_eq_derivation.jl,
                     # M test/mechanism_definitions_for_test_enzyme_derivation.jl,
                     # M test/test_rate_eq_derivation.jl
git add docs/src/deriving/cha_king_altman.md src/rate_eq_derivation.jl \
    test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_rate_eq_derivation.jl
cat > $P/replay/msg <<'EOF'
Break a segment root's tie by the canonical form order

Among a rapid-equilibrium segment's forms with the fewest bound metabolites and no
residual, `_segment_root` took the one whose internal name sorts first, which spells a
modifier copy after a longer name (EAinh after EAb). It now takes the first in the
canonical form order (`_form_index_key`). Another root scales a law's numerator and
denominator by one monomial, so its printed form may change with equal values.

The oracle's flip form, which the law's segment root fixes, follows the same order. The
Cha / King–Altman page says that the last tie goes to the form that comes first in the
mechanism's form order.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

**Unit 2: eb508a5** (row 4; refactor, source only)

- [ ] **Step 6: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary eb508a5 -- src/thermo.jl | git apply
```

- [ ] **Step 7: Green.** G `$L/T24_2_green.log` `test_joint_solve.jl test_rate_eq_derivation.jl`:
  453, 3,186 (a row permutation leaves the cycle basis unchanged).
- [ ] **Step 8: Commit.** K `eb508a5` `src/thermo.jl` ("Key the cycle basis's form rows by
  form, in the canonical order").

**Unit 3: 9c6d4cb** (row 5)

- [ ] **Step 9: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9c6d4cb -- test/test_mechanism_enumeration.jl test/test_types.jl | git apply
```

- [ ] **Step 10: Red** (background). Rd `$L/T24_3_red.log`
  `test_types.jl test_mechanism_enumeration.jl`. Expected: test_types.jl passes (its edit is
  a comment); test_mechanism_enumeration.jl, 13 errors under "_forms_where_free": 6 in the
  inner "_forms_where_free", 3 in "a copy binds at its own sites, apart from its
  reactant's", 4 in "a Theorell–Chance step frees B at E(A) and P at E(Q)" (the forms come
  back as Symbols, which `_form_label` cannot render).
- [ ] **Step 11: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9c6d4cb -- src/admission.jl src/mechanism_enumeration.jl src/rate_eq_derivation.jl src/types.jl | git apply
```

- [ ] **Step 12: Green** (background). G `$L/T24_3_green.log`
  `test_types.jl test_mechanism_enumeration.jl`: 649, 32,536.
- [ ] **Step 13: Commit.** K `9c6d4cb` `src/admission.jl src/mechanism_enumeration.jl
  src/rate_eq_derivation.jl src/types.jl test/test_mechanism_enumeration.jl
  test/test_types.jl` ("Identify forms by content in the reach, flip and dead-end sets").

**Task end**

- [ ] **Step 14: Full suite.** F 24: 38,909.
- [ ] **Step 15: Census.** C 24 `S3b_end`: as Task 23 (`census_check 24` prints
  `RESULT: PASS` and the ter-bi init line only).

### Task 25: Enumeration ties by content (rows 6-9) [S3b]

**Replays:** `[1a33d0f]`, `[7f1aa95]`, `[af01ff6]`, `[0266794]`.

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_productive_twin`, `_catalytic_topologies`,
  `_context_bipartitions`, the dead-end placements), `src/types.jl` (texts about internal
  names)
- Test: `test/test_mechanism_enumeration.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: `_form_index_key`, Task 24's `Species`-keyed sets.
- Produces: `_productive_twin(groups)` breaks a tie among twins by
  `argmin(_form_index_key, twins)`; `_catalytic_topologies(reaction)` sorts isomerization
  patterns and each topology's steps by `_form_index_key` pairs (the first isomerization
  stays the steady-state step); `_context_bipartitions(group)` orders residual contexts by
  (added names, subtracted names); dead-end placements sort by `_form_index_key`.

**Gate:** F 25 = 38,912; C 25 `S3b_end`: zero differences, order included.

**Unit 1: 1a33d0f** (row 6)

- [ ] **Step 1: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 1a33d0f -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 2: Red** (background). Rd `$L/T25_1_red.log` `test_mechanism_enumeration.jl`.
  Expected: 1 failure in "_productive_twin" (`F2(S)` picked).
- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 1a33d0f -- src/mechanism_enumeration.jl | git apply
```

- [ ] **Step 4: Green** (background). G `$L/T25_1_green.log`
  `test_mechanism_enumeration.jl`: 32,537.
- [ ] **Step 5: Commit.** K `1a33d0f` `src/mechanism_enumeration.jl
  test/test_mechanism_enumeration.jl` ("Pick a tied productive twin by the canonical form
  order").

**Unit 2: 7f1aa95** (row 7)

- [ ] **Step 6: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 7f1aa95 -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 7: Red** (background). Rd `$L/T25_2_red.log` `test_mechanism_enumeration.jl`.
  Expected: 1 failure in "isomerization patterns come in the canonical form order" (the
  ter-bi ping-pong patterns in the other order).
- [ ] **Step 8: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 7f1aa95 -- src/mechanism_enumeration.jl | git apply
```

- [ ] **Step 9: Green** (background). G `$L/T25_2_green.log`
  `test_mechanism_enumeration.jl`: 32,538.
- [ ] **Step 10: Commit.** K `7f1aa95` `src/mechanism_enumeration.jl
  test/test_mechanism_enumeration.jl` ("Order catalytic topologies and their steps by form
  content").

**Unit 3: af01ff6** (row 8)

- [ ] **Step 11: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary af01ff6 -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 12: Red** (background). Rd `$L/T25_3_red.log` `test_mechanism_enumeration.jl`.
  Expected: 1 failure in "_context_bipartitions" (old order A + B − P, A − P, none).
- [ ] **Step 13: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary af01ff6 -- src/mechanism_enumeration.jl | git apply
```

- [ ] **Step 14: Green** (background). G `$L/T25_3_green.log`
  `test_mechanism_enumeration.jl`: 32,539.
- [ ] **Step 15: Commit.** K `af01ff6` `src/mechanism_enumeration.jl
  test/test_mechanism_enumeration.jl` ("Order a split's residual contexts by content").

**Unit 4: 0266794** (row 9; refactor: placements dedup as sets, children are canonicalized
when built)

- [ ] **Step 16: Apply the test path, then the source** (the test edit is a comment)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 0266794 -- test/test_types.jl | git apply
git show --format= --binary 0266794 -- src/mechanism_enumeration.jl src/types.jl | git apply
```

- [ ] **Step 17: Green** (background). G `$L/T25_4_green.log`
  `test_types.jl test_mechanism_enumeration.jl`: 649, 32,539.
- [ ] **Step 18: Commit.** K `0266794` `src/mechanism_enumeration.jl src/types.jl
  test/test_types.jl` ("Sort the dead-end move's placements in the canonical form order").

**Task end**

- [ ] **Step 19: Full suite.** F 25: 38,912.
- [ ] **Step 20: Census.** C 25 `S3b_end`: zero differences, order included (row 7 brings
  the ter-bi init order).

### Task 26: Seed-base order, the R32 Idea, forms that share an internal name (rows 10-11) [S3b]

**Replays:** `[74a8fc6]`, `[61a1256]`, `[a2bbf86]` (prototype order: the ROADMAP commit
precedes row 11, the stage's last commit, K4).

**Files:**
- Modify: `src/types.jl` (`_canonical_order`, `_canonical_groups`,
  `_assert_distinct_form_names`), `src/mechanism_enumeration.jl` (`_seed_variants`),
  `ROADMAP.md` (Ideas: "Seed each identify fit from its mechanism's content")
- Test: `test/test_mechanism_enumeration.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: `_orient_tied_steps`, `_canonical_step_direction` (Task 3).
- Produces: `_canonical_order(reaction::EnzymeReaction, groups::Vector{Vector{Step}}) ->
  (groups, perm)` (orientation and both sorts, no representation check);
  `_canonical_groups(reaction, groups)` = empty-group check, `_canonical_order`, then the
  representation checks; `_seed_variants(m, rxn)` takes each base through
  `_canonical_order` before `_minimal_flips`; `_assert_distinct_form_names(reaction, steps)`
  compares rendered names only (two forms may share an internal name).

**Gate:** F 26 = 38,925; C 26 `S3b_end`: zero differences, order included; seed counts
264 / 239 / 7 / 250,855; checkpoint S3b.

**Unit 1: 74a8fc6** (row 10)

- [ ] **Step 1: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 74a8fc6 -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 2: Red** (background). Rd `$L/T26_1_red.log` `test_mechanism_enumeration.jl`.
  Expected: 1 failure in "_seed_variants: the variants follow each base's canonical group
  order" (variants 4 and 5 swapped).
- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 74a8fc6 -- src/mechanism_enumeration.jl src/types.jl | git apply
```

- [ ] **Step 4: Green** (background). G `$L/T26_1_green.log`
  `test_mechanism_enumeration.jl`: 32,540.
- [ ] **Step 5: Commit.** K `74a8fc6` `src/mechanism_enumeration.jl src/types.jl
  test/test_mechanism_enumeration.jl` ("Order each seed base's groups canonically before
  flipping them").

**Unit 2: 61a1256** (R32's Idea; ROADMAP only, nothing to run)

- [ ] **Step 6: Apply**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 61a1256 -- ROADMAP.md | git apply
git diff --stat
```

  Expected: `ROADMAP.md | 10 ++++++++++`; the entry sits under Ideas, after
  "...It must never filter children, since a phantom child can carry families its parent
  lacks."
- [ ] **Step 7: Commit.** K `61a1256` `ROADMAP.md` ("Record that identify's fit starting
  points follow processing order").

**Unit 3: a2bbf86** (row 11)

- [ ] **Step 8: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary a2bbf86 -- test/test_types.jl | git apply
```

- [ ] **Step 9: Red.** Rd `$L/T26_3_red.log` `test_types.jl`. Expected: 1 error in "two
  forms that render apart are accepted, whatever their internal names" (the constructor
  rejects E(I) beside E(Iinh): "share the internal name EIinh").
- [ ] **Step 10: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary a2bbf86 -- src/types.jl | git apply
```

- [ ] **Step 11: Green.** G `$L/T26_3_green.log` `test_types.jl`: 661.
- [ ] **Step 12: Commit.** K `a2bbf86` `src/types.jl test/test_types.jl` ("Accept forms that
  share an internal name but render apart").

**Task end**

- [ ] **Step 13: Full suite.** F 26: 38,925. test_identify_rate_equation.jl passes, so R32's
  seeding stays an Idea (K5).
- [ ] **Step 14: Census.** C 26 `S3b_end`: zero differences, order included.
- [ ] **Step 15: Seed counts.** Background, then `wait_log T26_seed_counts`:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh $L/T26_seed_counts.log 1800 --project=. --heap-size-hint=2500M $P/proto/scripts/S3a_seed_counts.jl
```

  Expected: `pp: 264 mechanisms`, `bb: 239 mechanisms`, `ub: 7 mechanisms`, `tt: 250855
  mechanisms in <t> s` with `<t>` about 31 (S3b: 30.9 s warm), `exit 0`.
- [ ] **Step 16: Checkpoint S3b** (tip a2bbf86). Run

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git diff --stat a2bbf86 HEAD -- src test docs ROADMAP.md Project.toml .claude
```

  Expected: S3a's 13 files with their line counts, plus
  docs/src/deriving/cha_king_altman.md (3 lines), 14 in all. Then the scratch-index check
  against a2bbf86, with both 216be7f pieces among the docs ahead:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
export GIT_INDEX_FILE=$P/replay/ckpt.idx
git read-tree a2bbf86
git diff 4f9a34a 7205ba4 -- docs/superpowers/specs ROADMAP.md | git apply --cached   # D1, D3
for c in 5c08ede 4a96d9a 0b5e637 68f46f0 33b94d8 cfe9211; do git show --format= --binary $c | git apply --cached; done   # D4
for c in 2b29fd2 be65542 0b2021d a84d194; do git show --format= --binary $c | git apply --cached; done   # docs ahead
git show --format= --binary 216be7f -- docs/src/deriving/mwc_allostery.md | git apply --cached   # docs ahead
git show --format= --binary 216be7f -- docs/src/deriving/cha_king_altman.md | git apply --cached   # docs ahead
git update-index --add --cacheinfo 100644,$(git rev-parse HEAD:docs/superpowers/plans/2026-10-16-graph-mechanism.md),docs/superpowers/plans/2026-10-16-graph-mechanism.md   # D2
git diff --cached --stat HEAD -- src test docs ROADMAP.md Project.toml .claude
rm -f $GIT_INDEX_FILE; unset GIT_INDEX_FILE
```

  Expected: every `git apply --cached` silent, and the last `git diff` prints nothing.
  Departures in S3b: none in code; Task 24's message is edited as written there, and K
  adds the Claude-Session trailer to the other main commits' messages, which lack it.

### Task 27: Tests read mechanisms through one reader [S3c]

**Replays:** `[8107bc0]`, `[990cb6a]`, `[279855d]`, `[f57370e]`, `[0e051b4]`. Test-only;
each unit green before (the previous unit) and after.

**Files:**
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`, `test/test_oracles.jl`,
  `test/test_dsl.jl`, `test/test_joint_solve.jl`, `test/test_mechanism_enumeration.jl`,
  `test/test_rate_eq_derivation.jl`, `test/test_types.jl`, `test/test_accessors.jl`,
  `test/test_identify_rate_equation.jl`

**Interfaces:**
- Consumes: Task 18's name accessors, Task 23's `_forms` order, `_flat_steps`,
  `_step_constants`, `_class_sides`, `_inactive_groups`, `_form_label`.
- Produces (FX): `_testhelper_graph(m) -> NamedTuple` with `forms` (labels, free E first),
  `free`, `steps` (each `(; from, to, consumed, released, consumed_modifier,
  released_modifier, rapid_equilibrium, group, names, inactive_names)`), `tags`,
  `inactive`, `catalytic_multiplicity`, `sites` (`(; multiplicity, ligands, allo_states)`),
  `substrates`, `products`, `modifiers`; `_testhelper_is_binding(s)`,
  `_testhelper_is_iso(s)`, `_testhelper_first_steps(g)`, `_testhelper_step_text(g, s)`,
  `_testhelper_placeholder_reaction(m)`. TE: `_testhelper_keeps_reaction(child, parent)`.
  Gone: FX's `_testhelper_flat_steps`, `_testhelper_enzyme_forms`; TE's
  `_testhelper_on_reaction`, `_testhelper_on_reaction_mult2`. Task 34 points
  `_testhelper_graph` at the graph; Task 53 makes it the tests' only reader.

**Gate:** F 27 = 38,950. No census (no src).

990cb6a's red, the reader's `free` pinned before the field existed (1 fail, 4 errors), is
recorded in `$P/proto/logs/S3c_red_free.log`; it ran on a working tree the prototype later
split into commits, so the replay does not reproduce it.

- [ ] **Step 1: Apply 8107bc0**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 8107bc0 -- test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_oracles.jl | git apply
```

- [ ] **Step 2: Green.** G `$L/T27_1_green.log` `test_oracles.jl`: 656 (as
  `$P/proto/logs/S3c_c1_oracles.log`).
- [ ] **Step 3: Commit.** K `8107bc0` `test/mechanism_definitions_for_test_enzyme_derivation.jl
  test/test_oracles.jl` ("Read test mechanisms through one adapter"; its message lacks the
  Claude-Session line, which K adds).
- [ ] **Step 4: Apply 990cb6a**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 990cb6a -- test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_dsl.jl test/test_joint_solve.jl test/test_mechanism_enumeration.jl test/test_oracles.jl test/test_rate_eq_derivation.jl test/test_types.jl | git apply
```

- [ ] **Step 5: Green** (background). G `$L/T27_2_green.log`
  `test_types.jl test_dsl.jl test_rate_eq_derivation.jl test_oracles.jl test_joint_solve.jl test_mechanism_enumeration.jl`:
  661, 251, 3,186, 662, 453, 32,540 (the init golden, the oracle sample and the 42 spec
  checks unchanged).
- [ ] **Step 6: Commit.** K `990cb6a` `test/mechanism_definitions_for_test_enzyme_derivation.jl
  test/test_dsl.jl test/test_joint_solve.jl test/test_mechanism_enumeration.jl
  test/test_oracles.jl test/test_rate_eq_derivation.jl test/test_types.jl` ("Read the test
  helpers' mechanisms through the adapter").
- [ ] **Step 7: Apply 279855d**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 279855d -- test/test_accessors.jl test/test_dsl.jl test/test_identify_rate_equation.jl test/test_joint_solve.jl test/test_mechanism_enumeration.jl test/test_rate_eq_derivation.jl test/test_types.jl | git apply
```

- [ ] **Step 8: Green** (background). G `$L/T27_3_green.log`
  `test_accessors.jl test_types.jl test_dsl.jl test_rate_eq_derivation.jl test_joint_solve.jl test_mechanism_enumeration.jl test_identify_rate_equation.jl`:
  155, 661, 251, 3,186, 453, 32,540, 446.
- [ ] **Step 9: Commit.** K `279855d` `test/test_accessors.jl test/test_dsl.jl
  test/test_identify_rate_equation.jl test/test_joint_solve.jl
  test/test_mechanism_enumeration.jl test/test_rate_eq_derivation.jl test/test_types.jl`
  ("Take each test's reaction or names, not a mechanism's reaction").
- [ ] **Step 10: Apply f57370e**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary f57370e -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 11: Green** (background). G `$L/T27_4_green.log`
  `test_mechanism_enumeration.jl`: 32,540; `grep -c _testhelper_on_reaction
  test/test_mechanism_enumeration.jl` prints 0.
- [ ] **Step 12: Commit.** K `f57370e` `test/test_mechanism_enumeration.jl` ("Build
  enumeration fixtures with the macros alone").
- [ ] **Step 13: Apply 0e051b4**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 0e051b4 -- test/test_dsl.jl test/test_mechanism_enumeration.jl test/test_types.jl | git apply
```

- [ ] **Step 14: Green** (background). G `$L/T27_5_green.log`
  `test_types.jl test_dsl.jl test_mechanism_enumeration.jl`: 661, 251, 32,540.
- [ ] **Step 15: Commit.** K `0e051b4` `test/test_dsl.jl test/test_mechanism_enumeration.jl
  test/test_types.jl` ("Read a mechanism's kinetic groups through steps(m)").
- [ ] **Step 16: Full suite.** F 27: 38,950.

### Task 28: R48's test items [S3c]

**Replays:** `[56f7871]`, `[bdd16e2]`, `[1cd8cd5]`, `[d52040f]`, `[380a8c4]`, `[e080b63]`.
Test-only, each green on arrival: 56f7871 and bdd16e2 fail only on e4b5927's source
(`$P/proto/logs/S3c_red_e4b.log`: "the free forms D, E ..." expected "E, D"; the Iinh
expansion aborts on the internal-name check), which Tasks 24 and 26 replaced.

**Files:**
- Test: `test/test_types.jl`, `test/test_mechanism_enumeration.jl`,
  `test/test_rate_eq_derivation.jl`, `test/test_joint_solve.jl`

**Interfaces:**
- Consumes: Task 27's reader; Task 24's one-flip message; Task 26's acceptance of E(I)
  beside E(Iinh).
- Produces: TE's `_testhelper_src_forms(part, ctx)` and `_testhelper_form_called` read form
  labels; `_testhelper_holds_iso(m)` reads `steps(m)` (D-S3c-9: about a million calls in
  the depth-2 population testset).

**Gate:** F 28 = 38,954; checkpoint S3c. No census (no src).

- [ ] **Step 1: Apply 56f7871** (m1: the one-flip message for free forms D and E)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 56f7871 -- test/test_types.jl | git apply
```

- [ ] **Step 2: Green.** G `$L/T28_1_green.log` `test_types.jl`: 662.
- [ ] **Step 3: Commit.** K `56f7871` `test/test_types.jl` ("Pin the one-flip message for
  free forms D and E").
- [ ] **Step 4: Apply bdd16e2** (m2: the seed with substrate Iinh and inhibitor I, its four
  children written out)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary bdd16e2 -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 5: Green** (background). G `$L/T28_2_green.log`
  `test_mechanism_enumeration.jl`: 32,543.
- [ ] **Step 6: Commit.** K `bdd16e2` `test/test_mechanism_enumeration.jl` ("Pin the
  children of a seed with substrate Iinh and inhibitor I").
- [ ] **Step 7: Apply 1cd8cd5** (m3: the `copy` local becomes `a_copy`)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 1cd8cd5 -- test/test_rate_eq_derivation.jl | git apply
```

- [ ] **Step 8: Green.** G `$L/T28_3_green.log` `test_rate_eq_derivation.jl`: 3,186.
- [ ] **Step 9: Commit.** K `1cd8cd5` `test/test_rate_eq_derivation.jl` ("Rename the
  segment-root testset's copy local").
- [ ] **Step 10: Apply d52040f** (m4: labels where tests read internal names)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary d52040f -- test/test_joint_solve.jl test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 11: Green** (background). G `$L/T28_4_green.log`
  `test_joint_solve.jl test_mechanism_enumeration.jl`: 453, 32,543.
- [ ] **Step 12: Commit.** K `d52040f` `test/test_joint_solve.jl
  test/test_mechanism_enumeration.jl` ("Read form labels where enumeration tests read
  internal names").
- [ ] **Step 13: Apply 380a8c4** (m5: the ter-ter timing figure leaves the comment)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 380a8c4 -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 14: Green** (background). G `$L/T28_5_green.log`
  `test_mechanism_enumeration.jl`: 32,543.
- [ ] **Step 15: Commit.** K `380a8c4` `test/test_mechanism_enumeration.jl` ("Drop the
  ter-ter init timing figure from its test comment").
- [ ] **Step 16: Apply e080b63** (`_testhelper_holds_iso` reads the steps)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary e080b63 -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 17: Green** (background). G `$L/T28_6_green.log`
  `test_mechanism_enumeration.jl`: 32,543. Its Time should fall back near Task 26
  Step 4's (D-S3c-9: 151 s against 142 s at a2bbf86, where the reader in this predicate
  had slowed the depth-2 population testset from 49 s to 81 s); informational, not a
  gate.
- [ ] **Step 18: Commit.** K `e080b63` `test/test_mechanism_enumeration.jl` ("Read
  isomerizations off the steps in the population predicate").
- [ ] **Step 19: Full suite.** F 28: 38,954.
- [ ] **Step 20: Checkpoint S3c** (tip e080b63). Run

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git diff --stat e080b63 HEAD -- src test docs ROADMAP.md Project.toml .claude
```

  Expected: exactly these 14 files: Project.toml 2, ROADMAP.md 20,
  docs/src/deriving/cha_king_altman.md 3, docs/src/deriving/mwc_allostery.md 8,
  docs/src/deriving/ping_pong.md 2, docs/src/deriving/re_vs_ss.md 6,
  docs/src/deriving/regulators.md 6, docs/src/developer.md 6,
  docs/src/identify/enumeration_engine.md 5, docs/src/identify/tutorial.md 7, the plan
  (D2), the spec 48 (D1), FX 35, TRD 31. Then the scratch-index check against e080b63,
  with the docs ahead of Task 26 Step 16:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
export GIT_INDEX_FILE=$P/replay/ckpt.idx
git read-tree e080b63
git diff 4f9a34a 7205ba4 -- docs/superpowers/specs ROADMAP.md | git apply --cached   # D1, D3
for c in 5c08ede 4a96d9a 0b5e637 68f46f0 33b94d8 cfe9211; do git show --format= --binary $c | git apply --cached; done   # D4
for c in 2b29fd2 be65542 0b2021d a84d194; do git show --format= --binary $c | git apply --cached; done   # docs ahead
git show --format= --binary 216be7f -- docs/src/deriving/mwc_allostery.md | git apply --cached   # docs ahead
git show --format= --binary 216be7f -- docs/src/deriving/cha_king_altman.md | git apply --cached   # docs ahead
git update-index --add --cacheinfo 100644,$(git rev-parse HEAD:docs/superpowers/plans/2026-10-16-graph-mechanism.md),docs/superpowers/plans/2026-10-16-graph-mechanism.md   # D2
git diff --cached --stat HEAD -- src test docs ROADMAP.md Project.toml .claude
rm -f $GIT_INDEX_FILE; unset GIT_INDEX_FILE
```

  Expected: every `git apply --cached` silent, and the last `git diff` prints nothing.
  Departures in S3c: none in code; K adds the Claude-Session trailer to 8107bc0's
  message, the last prototype commit without it.

## Tasks 29-43: S3d1 graph core, S3d2 swap, S3e constructor checks, S3f derivation on the graph

Every git block of Tasks 29-43 ran as written, in order, on a scratch clone at Task 28's
end (`$P/plan/evidence/p3part/validate.py`): every apply was clean, the four checkpoint
checks printed nothing, and the 39 commits end in the tree the prototype holds at 7fe5c6d
plus the listed departures.

### Task 29: The reaction graph's value types, content order and identity [S3d1]

**Replays:** `[e12daba]` (S3d1 §9 task 1).

**Files:**
- Create: `src/reaction_graph.jl` (`module ReactionGraph`), `test/test_reaction_graph.jl`
- Modify: `src/EnzymeRates.jl` (includes `reaction_graph.jl` after `types.jl`),
  `test/runtests.jl` (includes `test_reaction_graph.jl` after `test_types.jl`)

**Interfaces:**
- Consumes: `Met`, `RegulatorySite`, `_VALID_CAT_ALLO_STATES` (types.jl),
  `_check_metabolite_names` (Task 8).
- Produces, in `EnzymeRates.ReactionGraph` beside the Step types: `Counts = UInt128` with
  `copies(c, i)`, `total(c)`, `entries(c)`; `Form(conformation, bound, added, subtracted)`;
  `Edge(from, to, consumed, released, rapid_equilibrium)`; `Mechanism(substrates, products,
  table, conformations, forms, edges, group)`; `AllostericMechanism(catalytic,
  cat_allo_states, catalytic_multiplicity, regulatory_sites)`, whose inner constructor runs
  `_check_tags` and `_check_ligands(sites)`, sorts the sites by `_site_key` and copies the
  tags; `_seq_less`, `isless(::Form, ::Form)`, `isless(::Edge, ::Edge)`, `_reversed(::Edge)`,
  `ngroups(m)`; content `hash` and `==` for forms, edges and both mechanisms. TRG: `const RG
  = EnzymeRates.ReactionGraph` and its first five testsets.

**Gate:** full 39,000 (38,954 + TRG's 45 + one chokepoint-guard assertion, since the guard
reads each src file), pristine; JET on ReactionGraph, 0 reports. No census: nothing outside
the module reads it.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary e12daba -- test/runtests.jl test/test_reaction_graph.jl | git apply
git status --short
```

Expected: ` M test/runtests.jl` and `?? test/test_reaction_graph.jl`.

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T29_red test/test_reaction_graph.jl`. Then
`wait_log T29_red`, `summary T29_red`.

Expected (`$P/proto/logs/S3d_c1_red.log`): `test/test_reaction_graph.jl` records 1 error
and no pass: the file stops at `const RG = EnzymeRates.ReactionGraph` ("Got exception
outside of a @test", `LoadError: UndefVarError: ReactionGraph not defined in EnzymeRates`),
so none of its five testsets runs; `exit 1`.

- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary e12daba -- src/EnzymeRates.jl src/reaction_graph.jl | git apply
git status --short
```

Expected: ` M src/EnzymeRates.jl`, `?? src/reaction_graph.jl` and the two test paths.

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T29_green test/test_reaction_graph.jl
test/test_types.jl`. Then `wait_log T29_green`, `summary T29_green`.

Expected: `test/test_reaction_graph.jl | 45  45` (S3d_c1_green.log), `test/test_types.jl |
663  663`, `exit 0`.

- [ ] **Step 5: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/EnzymeRates.jl src/reaction_graph.jl test/runtests.jl test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Add the reaction graph's value types, content order and identity

The graph core is built beside today's types, in module ReactionGraph, so a
differential test can compare the two before the swap. This first step holds
the types (packed metabolite counts, forms, edges, the plain and allosteric
mechanisms), the content order that canonical form sorts by (entry lists in
lexicographic order, not by their packed numbers) and process-stable content
hashes.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

- [ ] **Step 6: Full suite**

Background: `. .../replay/env.sh; full T29_full`. Then `wait_log T29_full` until it ends,
`summary T29_full; pristine T29_full`.

Expected: `EnzymeRates.jl | 39000  39000`, `exit 0`; `pristine` prints `12` only.

- [ ] **Step 7: JET on ReactionGraph** (Runners, `<NN>` = 29). Expected: both
  `targets` lines end `0 reports`.

### Task 30: Canonical form and checked constructors [S3d1]

**Replays:** `[2db9d74]`, `[584d7bf + d1fb2df]`, `[43c50fc]` (S3d1 §9 tasks 2-4; d1fb2df is
review round 1's m1, folded into the constructor task whose check it pins).

**Files:**
- Modify: `src/reaction_graph.jl` (canonical form, the constructors and their checks, the
  graph renderers), `src/types.jl` (`RegulatorySite` rejects a ligand written twice)
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (raw values of the six
  specs that lacked one: 11b, 11c, 11d and the three LDH i-state specs, which move into
  `let` blocks), `test/test_reaction_graph.jl`

**Interfaces:**
- Consumes: Task 29's types; `RawForm`, `RawMechanism`, `RawAllostericMechanism` (Task 5).
- Produces (ReactionGraph): `_canonicalize(raw::RawMechanism) -> (; mechanism, edge_of,
  flipped, group_of, form_of, misaligned)` with `_metabolite_table`, `_counts`, `_form`,
  `_sorted_forms`, `_role_masks`, `_orientation`, `_align_tied_groups!`, `_number_groups`;
  `Mechanism(raw::RawMechanism)` and `AllostericMechanism(raw::RawAllostericMechanism)`, each
  tag following its written label to its canonical group; `_check_declarations(raw)`,
  `_check_representation(m, misaligned, modifiers)` (rows 13-23) with `_occupied`,
  `_check_distinct_form_names`, `_check_uniform_groups`, `_check_each_reaction_once` (keyed
  by the unordered pair of an edge's sides); the reads `substrates`, `products`,
  `modifiers`, `_is_reactant`, `_role_masks(m)`, `catalytic`, `_group_members`,
  `_has_residual`, `_is_binding`, `_bound_entry`, `_is_plain_binding`, `_binds_by_name`;
  the renderers `_entry_label`, `_form_name`, `_form_text`, `_side_name`, `_sides`.
  types.jl: `RegulatorySite` errors "RegulatorySite: ligand R appears twice in one
  regulatory site" (row 2).

**Gate:** full 39,314 (39,000 + TRG 185 + 127 + 1 + d1fb2df's 1), pristine; JET on
ReactionGraph, 0 reports; census against `S3b_end` with the map, zero differences
(`RegulatorySite` in types.jl changed; no enumerated site repeats a ligand).

**Unit 1: 2db9d74** (canonical form)

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2db9d74 -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_reaction_graph.jl | git apply
```

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T30_1_red test/test_reaction_graph.jl`. Then
`wait_log T30_1_red`, `summary T30_1_red`.

Expected (S3d_c2_red.log): `test/test_reaction_graph.jl` 45 pass, 3 fail, 48 errors (96),
every one an `UndefVarError: _canonicalize not defined in EnzymeRates.ReactionGraph`:
"every writing of a mechanism gives one canonical graph" 42 errors (one per spec),
"canonical form maps each written step to its edge, group and forms" 1, "canonical form
runs each step physical-forward, tied steps as one" 3, and "a graph holds at most 32 table
entries and 15 copies of one entry" 2 errors and 3 fails (its `@test_throws` see the
UndefVarError). The five testsets of Task 29 pass.

- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2db9d74 -- src/reaction_graph.jl | git apply
```

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T30_1_green test/test_reaction_graph.jl
test/test_rate_eq_derivation.jl`. Then `wait_log T30_1_green`, `summary T30_1_green`.

Expected: TRG 230 (S3d_c2_green2.log), TRD 3,186 (S3d1 §4: no other file's count moves),
`exit 0`.

- [ ] **Step 5: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/reaction_graph.jl test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Canonicalize raw input into the reaction graph

Every writing of a mechanism (step order, written directions, group labels,
declaration order) now gives one graph with equal fields and hashes: forms in
content order with free E first, edges physical-forward by today's tiers with
tied groups turned as one, groups numbered by their first edge. Canonical form
also maps each written step to its edge, group and forms, which the textbook
bridge and the moves will read. It errors only where no encoding exists: more
than 32 table entries or 15 copies of one.

The six specs that lacked a raw value gain one, so every spec's written form
feeds the graph tests.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 2: 584d7bf + d1fb2df** (checked constructors; the once-each-way case)

- [ ] **Step 6: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 584d7bf -- test/test_reaction_graph.jl | git apply
git show --format= --binary d1fb2df -- test/test_reaction_graph.jl | git apply
```

Expected: no output.

- [ ] **Step 7: Red**

Background: `. .../replay/env.sh; focused T30_2_red test/test_reaction_graph.jl`. Then
`wait_log T30_2_red`, `summary T30_2_red`.

Expected: TRG 230 pass, 39 fail, 45 errors (314): S3d_c3_red.log's 38 fails and 45
errors, all from `MethodError: no method matching EnzymeRates.ReactionGraph.Mechanism(::
EnzymeRates.RawMechanism)` (or `AllostericMechanism`): "the constructors build one value
from every writing of a mechanism" 43 errors, "the constructors reject what no graph
represents, with today's messages" 38 fails and 1 error, "a form's name and text read as
constants and the macros write them" 1 error; plus d1fb2df's case, 1 more fail in the
rows testset ("kinetic groups 3 and 6 both hold the reaction E ⇌ F", thrown a
MethodError).

- [ ] **Step 8: Apply the source** (d1fb2df has none)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 584d7bf -- src/reaction_graph.jl | git apply
```

- [ ] **Step 9: Green**

Background: `. .../replay/env.sh; focused T30_2_green test/test_reaction_graph.jl`. Then
`wait_log T30_2_green`, `summary T30_2_green`.

Expected: TRG 358 (S3d_c3_green2.log's 357 + d1fb2df's 1), `exit 0`.

- [ ] **Step 10: Commit** (message edited: one paragraph for d1fb2df)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/reaction_graph.jl test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Build and check the reaction graph from raw input

Mechanism(raw) and AllostericMechanism(raw) check the declarations, then the
representation of the canonical graph: every name declared in its role,
residuals by role, no self-loop or metabolite on both sides, no role change in
a binding, tied steps aligned, forms that render apart, uniform kinetic
groups and each reaction once. Each tag follows its written label to its
canonical group. The messages are today's, checked against today's
constructors on the same raw values; a form renders by label, so a role
change into P's modifier copy names EPmod where today names EPinh.

The renderers that name forms and step sides from the graph land here, since
the messages read them.

The once-per-reaction check keys a reaction by the unordered pair of its
edge's sides, so a reaction one group stores E -> F and another F -> E is
caught too. A test ties E -> F with E(I) -> F(I) and keeps a lone F -> E,
which a key by the ordered pair would accept.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 3: 43c50fc** (row 2 in `RegulatorySite`)

- [ ] **Step 11: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 43c50fc -- test/test_reaction_graph.jl | git apply
```

- [ ] **Step 12: Red**

Background: `. .../replay/env.sh; focused T30_3_red test/test_reaction_graph.jl`. Then
`wait_log T30_3_red`, `summary T30_3_red`.

Expected (S3d_c4_red.log): TRG 358 pass, 1 fail: the row-2 case (a ligand written twice
on one site) throws `AllostericMechanism: ligand R appears in two distinct regulatory
sites; rendered Kreg names would collide. ...` where it expects RegulatorySite's message.

- [ ] **Step 13: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 43c50fc -- src/types.jl | git apply
```

- [ ] **Step 14: Green**

Background: `. .../replay/env.sh; focused T30_3_green test/test_reaction_graph.jl
test/test_types.jl`. Then `wait_log T30_3_green`, `summary T30_3_green`.

Expected: TRG 359, TT 663, `exit 0`.

- [ ] **Step 15: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/types.jl test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Reject a ligand written twice on one regulatory site

RegulatorySite accepted a repeated ligand, and the allosteric mechanism then
reported it as a ligand "in two distinct regulatory sites", which the site
did not have. The site now says what is wrong. The graph's constructors build
their sites with it, so they report it too.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Task end**

- [ ] **Step 16: Full suite.** Background: `. .../replay/env.sh; full T30_full`; `wait_log
  T30_full`; `summary T30_full; pristine T30_full`. Expected: `EnzymeRates.jl | 39314
  39314`, `exit 0`; `pristine` prints `12` only.
- [ ] **Step 17: JET on ReactionGraph** (`<NN>` = 30). Expected: both lines `0 reports`.
- [ ] **Step 18: Census.** Background: `. .../replay/env.sh; census_run 30 S3b_end
  census/rename_phase3.toml`; `wait_log $P/census/logs/cmp_replay_T30.log`; `census_check
  30`. Expected: `RESULT: PASS` and nothing else.

### Task 31: Sig, compiled types, lifts and _raw [S3d1]

**Replays:** `[6177781 + 1bd88e8]`, `[f8eb0df]` (S3d1 §9 tasks 5-6; §9 folds 1bd88e8 into
the Sig task).

**Files:**
- Modify: `src/reaction_graph.jl`
- Test: `test/test_reaction_graph.jl`

**Interfaces:**
- Consumes: Task 30's constructors and canonical form.
- Produces (ReactionGraph): `AbstractEnzymeMechanism`, `EnzymeMechanism{Sig}`,
  `AllostericEnzymeMechanism{Sig}` (module-local until the swap); `_TAG_CODES`,
  `_sig_name!`, `_sig_counts!`, `_sig_mechanism!`, `_sig_of(m)` and `_sig_of(am)` (one
  homogeneous `Int` tuple: names as code points, counts as two halves); `_SigReader` and
  its readers; `_mechanism_from_sig(sig::Int...)`, `_allosteric_mechanism_from_sig(sig::
  Int...)` (varargs, so each compiles once whatever the length); `EnzymeMechanism(m)`,
  `AllostericEnzymeMechanism(am)`, `compile_mechanism`, the lifts `Mechanism(em)` and
  `AllostericMechanism(aem)`, `_concrete`, `catalytic_mechanism`,
  `catalytic_multiplicity(aem)`; `_raw(m)` and `_raw(am)` (every edge a step in its stored
  direction, labelled with its group, every vector fresh). TRG: `_testhelper_long_named(n)`.

**Gate:** full 39,703 (39,314 + 220 + 169), pristine; JET on ReactionGraph, 0 reports. No
census: only the module changes, and nothing outside it reads it.

**Unit 1: 6177781 + 1bd88e8** (the Sig; the compile-once test on lengths no other test
decodes)

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 6177781 -- test/test_reaction_graph.jl | git apply
git show --format= --binary 1bd88e8 -- test/test_reaction_graph.jl | git apply
```

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T31_1_red test/test_reaction_graph.jl`. Then
`wait_log T31_1_red`, `summary T31_1_red`.

Expected (S3d_c5_red.log): TRG 359 pass, 45 errors, each an `UndefVarError` in
ReactionGraph (`_sig_of`, `compile_mechanism`, a decoder or a compiled type): "a graph's
Sig decodes to the graph, and equal graphs share one type" 44, "the Sig decoders compile
once, whatever a Sig's length" 1. 1bd88e8's evidence is recorded, not rerun: a decoder
that takes an `NTuple{N, Int}` fails the compile-once test (S3d_mut_ntuple_decoder3.log:
1 fail, and 33 errors in the calls that splat a Sig).

- [ ] **Step 3: Apply the source** (1bd88e8 has none)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 6177781 -- src/reaction_graph.jl | git apply
```

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T31_1_green test/test_reaction_graph.jl`. Then
`wait_log T31_1_green`, `summary T31_1_green`.

Expected: TRG 579 (S3d_c5_green2.log's 578 + d1fb2df's 1), `exit 0`.

- [ ] **Step 5: Commit** (message edited: its last sentence describes 1bd88e8's test)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/reaction_graph.jl test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Encode the reaction graph as one integer tuple, its compiled type's parameter

A graph's Sig is a homogeneous Int tuple (names as code points, packed counts
as two halves), the only parameter of EnzymeMechanism and of
AllostericEnzymeMechanism, so equal graphs share one compiled type. Decoding
trusts the Sig and runs no check but the composition constructor's tags and
sites. The lifts splat the Sig with its length left open, so each varargs
decoder compiles once whatever the Sig's length; a test lifts uni-uni
mechanisms whose substrate's name holds hundreds of characters, Sig lengths no
other test decodes, and counts no new method instance, which a decoder that
takes its Sig as an NTuple{N, Int} would add.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 2: f8eb0df** (`_raw`)

- [ ] **Step 6: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary f8eb0df -- test/test_reaction_graph.jl | git apply
```

- [ ] **Step 7: Red**

Background: `. .../replay/env.sh; focused T31_2_red test/test_reaction_graph.jl`. Then
`wait_log T31_2_red`, `summary T31_2_red`.

Expected (S3d_c6_red.log): TRG 579 pass, 43 errors, each `UndefVarError: _raw not defined
in EnzymeRates.ReactionGraph` ("a graph written back as raw input builds the same graph").

- [ ] **Step 8: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary f8eb0df -- src/reaction_graph.jl | git apply
```

- [ ] **Step 9: Green**

Background: `. .../replay/env.sh; focused T31_2_green test/test_reaction_graph.jl`. Then
`wait_log T31_2_green`, `summary T31_2_green`.

Expected: TRG 748 (S3d_c6_green.log's 747 + 1), `exit 0`.

- [ ] **Step 10: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/reaction_graph.jl test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Write a reaction graph back as raw input

The enumeration moves will build children from raw input, so a graph must
write itself back: every edge a step in its stored direction, labelled with
its group, the declared modifiers its modifier entries. Building the raw
value again gives the same graph, unchanged in edge order, direction and
group numbering, and every vector is fresh, so a move that edits the raw
value never edits its parent.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Task end**

- [ ] **Step 11: Full suite.** Background: `. .../replay/env.sh; full T31_full`; `wait_log
  T31_full`; `summary T31_full; pristine T31_full`. Expected: `EnzymeRates.jl | 39703
  39703`, `exit 0`; `pristine` prints `12` only.
- [ ] **Step 12: JET on ReactionGraph** (`<NN>` = 31). Expected: both lines `0 reports`.

### Task 32: The graph against the Step types; the constructor's speed gate [S3d1]

**Replays:** `[bf41ad0]`, `[6ca71bd]`, `[1bf2fee]`, `[03fad30]`, `[94eaa5d]`, `[5298b84]`,
each its own commit (S3d1 §9 tasks 7-8; 1bf2fee lands after 6ca71bd, whose TRG lines it
reads, and the review round's commits keep theirs).

**Files:**
- Create: `test/fixtures/phase3_canonical_golden.txt` (472 lines, the 42 specs' canonical
  text)
- Modify: `src/reaction_graph.jl`
- Test: `test/test_reaction_graph.jl`

**Interfaces:**
- Consumes: Tasks 29-31; the test reader `_testhelper_graph` (Task 27);
  `$P/replay/mutate.py` (Task 0 Step 8): `tier3`, `no_alignment`, `first_tag`.
- Produces (ReactionGraph): `metabolites(m)`, `metabolites(am)`, `cat_allo_states(am)`,
  `cat_allo_state(am, g)`, `catalytic_multiplicity(am)`, `regulatory_sites(am)`;
  `_binding_ends(m, e) -> Tuple{Int, UInt128, UInt128}`, read by `_is_plain_binding` and
  `_binds_by_name`. TRG: the differential against the Step types
  (`_testhelper_step_view`, `_testhelper_step_built`, 03fad30's names), the canonical
  golden, the tier-3 and tied-group writings, R17's gate (median `Mechanism(raw)` over the
  specs under 10 µs).

**Gate:** full 39,761 (39,703 + 50 + 1 + 6 + 0 + 1 + 0), pristine; JET on ReactionGraph, 0
reports; census against `S3b_end`, zero differences, order included; the constructor's
median over the specs, about 3.2 µs (S3d_measure.log: 3.17); checkpoint S3d1.

**Unit 1: bf41ad0** (the differential and the canonical golden)

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary bf41ad0 -- test/fixtures/phase3_canonical_golden.txt \
    test/test_reaction_graph.jl | git apply
git status --short
```

Expected: ` M test/test_reaction_graph.jl`, `?? test/fixtures/phase3_canonical_golden.txt`.

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T32_1_red test/test_reaction_graph.jl`. Then
`wait_log T32_1_red`, `summary T32_1_red`.

Expected: TRG 749 pass, 43 errors, each `UndefVarError` of `metabolites` or
`cat_allo_states` in ReactionGraph, in the differential testset. S3d_c7_red.log shows one
more failure, the golden text against a fixture with a trailing blank line, which the
prototype fixed before committing the fixture: here the golden comparison passes, and a
failure there means the fixture applied wrong.

- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary bf41ad0 -- src/reaction_graph.jl | git apply
```

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T32_1_green test/test_reaction_graph.jl`. Then
`wait_log T32_1_green`, `summary T32_1_green`.

Expected: TRG 798 (S3d_c7_green3.log's 797 + 1), `exit 0`. The run read the golden, not
wrote it (it bootstraps only a missing file): `git show
bf41ad0:test/fixtures/phase3_canonical_golden.txt | cmp - test/fixtures/phase3_canonical_golden.txt`
prints nothing.

- [ ] **Step 5: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/reaction_graph.jl test/fixtures/phase3_canonical_golden.txt \
    test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Check the reaction graph against today's canonical form

The graph must reproduce today's canonical mechanisms exactly before it can
replace them. For the 42 specs and for the census's first populations (the
init sets of its six reactions and the 5,505 children of the bi-bi seeds, each
written against its stored form), the graph's forms, edges, groups, tags,
sites and names equal today's, read through the test adapter, in order; each
graph also decodes from its Sig and rebuilds from its raw form. A free form
labelled before E is checked too, since no population holds one.

The specs' canonical text goes to test/fixtures/phase3_canonical_golden.txt,
so the comparison outlives today's types. The graph gains the concentration
names and the allosteric accessors the comparison reads.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 2: 6ca71bd** (R17's gate on existing code: green on arrival)

- [ ] **Step 6: Apply and run**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 6ca71bd -- test/test_reaction_graph.jl | git apply
```

Background: `. .../replay/env.sh; focused T32_2_green test/test_reaction_graph.jl`. Then
`wait_log T32_2_green`, `summary T32_2_green`. Expected: TRG 799, `exit 0`.

- [ ] **Step 7: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Gate the graph constructor's speed at 10 µs per call

Every enumerated child goes through the constructor, so its speed bounds the
enumeration's. The median over the spec raws of the time per call (each the
best of five batches) fails the test above 10 µs (R17); locally it is about
3 µs, against about 9 µs for today's constructor on the same mechanisms.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 3: 1bf2fee** (tier 3 and tied groups against the Step types' constructor; green on
arrival, red only on its mutants: S3d_mut2_tier3.log 8 fail, S3d_mut2_no_alignment.log 4
fail and 2 errors, both of 804)

- [ ] **Step 8: Apply and run**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 1bf2fee -- test/test_reaction_graph.jl | git apply
```

Background: `. .../replay/env.sh; focused T32_3_green test/test_reaction_graph.jl`. Then
`wait_log T32_3_green`, `summary T32_3_green`. Expected: TRG 805, `exit 0`.

- [ ] **Step 9: Red on its two mutants** (S3d1 §9 task 7: tier 3 run from the later form,
  the tied-group alignment skipped; the testset titles still say "today's" until unit 4)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py tier3
```

Expected: `mutated src/reaction_graph.jl (tier3)`. Background: `. .../replay/env.sh;
focused T32_3_mut_tier3 test/test_reaction_graph.jl`; `wait_log T32_3_mut_tier3`;
`summary T32_3_mut_tier3`. Expected (`$P/proto/logs/S3d_mut2_tier3.log`, there 8 of 804):
TRG 797 pass, 8 fail, 0 error: "canonical form runs each step physical-forward, tied steps
as one" 3, "the constructors reject what no graph represents, with today's messages" 1,
"the graph orients as today where tier 3 and tied groups decide" 4. Then:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py tier3 --revert
python3 $P/replay/mutate.py no_alignment
```

Expected: `restored src/reaction_graph.jl`, `mutated src/reaction_graph.jl
(no_alignment)`. Background: `. .../replay/env.sh; focused T32_3_mut_align
test/test_reaction_graph.jl`; `wait_log T32_3_mut_align`; `summary T32_3_mut_align`.
Expected (`$P/proto/logs/S3d_mut2_no_alignment.log`): TRG 799 pass, 4 fail, 2 errors: the
first two testsets above 3 fail and 1 fail, and "the graph orients as today where tier 3
and tied groups decide" 2 errors. Then, whatever the run printed:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py no_alignment --revert
git status --short
```

Expected: `restored src/reaction_graph.jl`, then ` M test/test_reaction_graph.jl` only.

- [ ] **Step 10: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Compare tier-3 and tied-group orientation with today's constructor

No spec and no census population holds a step that only form order orients,
or a tied kinetic group whose alignment turns a step: mutants that reverse
tier 3 or skip the alignment passed the differential and failed only the
graph's own orientation tests, whose expected flips were derived by hand.
The graph and today's constructor now orient the same raw values in both
writings: a tier-3 isomerization between substrate-bound forms, a mirrored
isomerization tied to F → E, and a metabolite exchange at two conformations,
plain and allosteric. Both mutants fail it.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 4: 03fad30** (names the comparison types; test-only, green before and after)

- [ ] **Step 11: Apply and run**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 03fad30 -- test/test_reaction_graph.jl | git apply
grep -ciE "today|\byet\b" test/test_reaction_graph.jl || true
```

Expected: `0`. Background: `. .../replay/env.sh; focused T32_4_green
test/test_reaction_graph.jl`. Then `wait_log T32_4_green`, `summary T32_4_green`.
Expected: TRG 805, `exit 0`.

- [ ] **Step 12: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Name the graph tests' comparison types instead of a point in time

The test file called the Step-based mechanisms "today's types" in its
header, helper names, docstrings and testset titles. Names and comments
must stay true as the code changes, so they now say what they compare
against: the Step types, the mechanisms built from Step and Species.
The row-check testset title also claimed every message matches the
Step types', which not every row does (row 19 names the form EPmod
where the Step types name EPinh); it now says each row has its message.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 5: 94eaa5d** (the per-group tag read; green on arrival, red only on its mutant:
S3d1r1_mut_first_tag.log, `[:EqualAI, :EqualAI, :EqualAI] == [:EqualAI, :EqualAI,
:OnlyA]`)

- [ ] **Step 13: Apply and run**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 94eaa5d -- test/test_reaction_graph.jl | git apply
```

Background: `. .../replay/env.sh; focused T32_5_green test/test_reaction_graph.jl`. Then
`wait_log T32_5_green`, `summary T32_5_green`. Expected: TRG 806 (S3d1r1_m12.log),
`exit 0`.

- [ ] **Step 14: Red on its mutant** (every group read with the first group's tag)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py first_tag
```

Expected: `mutated src/reaction_graph.jl (first_tag)`. Background: `. .../replay/env.sh;
focused T32_5_mut test/test_reaction_graph.jl`; `wait_log T32_5_mut`; `summary T32_5_mut`.
Expected (`$P/proto/logs/S3d1r1_mut_first_tag.log`): TRG 805 pass, 1 fail, in "an allosteric
graph sorts its sites and checks its tags and ligands", evaluated `[:EqualAI, :EqualAI,
:EqualAI] == [:EqualAI, :EqualAI, :OnlyA]`. Then, whatever the run printed:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py first_tag --revert
git status --short
```

Expected: `restored src/reaction_graph.jl`, then ` M test/test_reaction_graph.jl` only.

- [ ] **Step 15: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Test the allosteric graph's per-group tag read

cat_allo_state(am, g) had no test and no caller yet in the module; a
mutant that returns the first group's tag for every group passed the
graph tests. The composition-constructor testset now reads each group's
tag through it.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 6: 5298b84** (`_binding_ends`; a refactor, green before and after)

- [ ] **Step 16: Apply and run**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 5298b84 -- src/reaction_graph.jl | git apply
```

Background: `. .../replay/env.sh; focused T32_6_green test/test_reaction_graph.jl`. Then
`wait_log T32_6_green`, `summary T32_6_green`. Expected: TRG 806, `exit 0`
(S3d1r1_m5.log ran TRG and TT at 5298b84: 1,469 = 806 + 663).

- [ ] **Step 17: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Read a binding's entry and forms once for both binding checks

_is_plain_binding and _binds_by_name each read the bound entry, the two
forms and compared their residuals. _binding_ends does that once and
returns the entry (0 unless the edge is a binding that keeps its
residual) with the two forms' bound entries, a concrete tuple.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Task end**

- [ ] **Step 18: Full suite.** Background: `. .../replay/env.sh; full T32_full`; `wait_log
  T32_full`; `summary T32_full; pristine T32_full`. Expected: `EnzymeRates.jl | 39761
  39761`, `exit 0`; `pristine` prints `12` only.
- [ ] **Step 19: JET on ReactionGraph** (`<NN>` = 32). Expected: both lines `0 reports`.
- [ ] **Step 20: Census.** Background: `. .../replay/env.sh; census_run 32 S3b_end
  census/rename_phase3.toml`; `wait_log $P/census/logs/cmp_replay_T32.log`; `census_check
  32`. Expected: `RESULT: PASS` and nothing else.
- [ ] **Step 21: Constructor speed**

Background:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh $L/T32_measure.log 1800 --project=. --heap-size-hint=2500M $P/s3d/scripts/measure.jl
```

Then `wait_log T32_measure` and `grep -E 'constructor|_canonicalize' $L/T32_measure.log`.
Expected (S3d_measure.log): `specs constructor (µs): (median = <m>, ...)` with `<m>` near
3.2 and under 10 (R17); `bi-bi init constructor` near 3.7; `today's constructor` near 9.2;
`exit 0`. Report the medians.

- [ ] **Step 22: Checkpoint S3d1** (tip 5298b84)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git diff --stat 5298b84 HEAD -- src test docs ROADMAP.md Project.toml .claude
```

Expected: the 14 files of the S3c checkpoint (Task 28 Step 20) with the same line counts:
Project.toml 2, ROADMAP.md 20, docs/src/deriving/cha_king_altman.md 3,
docs/src/deriving/mwc_allostery.md 8, docs/src/deriving/ping_pong.md 2,
docs/src/deriving/re_vs_ss.md 6, docs/src/deriving/regulators.md 6, docs/src/developer.md
6, docs/src/identify/enumeration_engine.md 5, docs/src/identify/tutorial.md 7, the plan
(D2), the spec 48, test/mechanism_definitions_for_test_enzyme_derivation.jl 35,
test/test_rate_eq_derivation.jl 31. No src file and no TRG: S3d1's code equals the
prototype's. Then build the expected tree in a scratch index and compare:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
export GIT_INDEX_FILE=$P/replay/ckpt.idx
git read-tree 5298b84
git diff 4f9a34a 7205ba4 -- docs/superpowers/specs ROADMAP.md | git apply --cached   # D1, D3
for c in 5c08ede 4a96d9a 0b5e637 68f46f0 33b94d8 cfe9211; do git show --format= --binary $c | git apply --cached; done   # D4
for c in 2b29fd2 be65542 0b2021d a84d194; do git show --format= --binary $c | git apply --cached; done   # docs ahead
git show --format= --binary 216be7f -- docs/src/deriving/mwc_allostery.md docs/src/deriving/cha_king_altman.md | git apply --cached
git update-index --add --cacheinfo 100644,$(git rev-parse HEAD:docs/superpowers/plans/2026-10-16-graph-mechanism.md),docs/superpowers/plans/2026-10-16-graph-mechanism.md   # D2
git diff --cached --stat HEAD -- src test docs ROADMAP.md Project.toml .claude
rm -f $GIT_INDEX_FILE; unset GIT_INDEX_FILE
```

Expected: every `git apply --cached` silent and the `git diff --cached` prints nothing.
Departures in S3d1: none in code; the messages of Task 30 unit 2 and Task 31 unit 1 carry
their folds, as written there.

### Task 33: Before the swap: tier-3 texts, names.jl, the pivot block [S3d2]

**Replays:** `[f7c1587]`, `[6214bb1]`, `[b0a4b29]`, each its own commit (S3d2 §9 tasks 1-3,
which must precede the swap).

**Files:**
- Create: `src/names.jl` (the naming chokepoint and the label renderers, moved unchanged)
- Modify: `src/types.jl` (the moved code leaves), `src/EnzymeRates.jl` (includes
  `names.jl` after `types.jl`)
- Test: `test/test_reaction_graph.jl`, `test/test_pivot_priority_regression.jl`

**Interfaces:**
- Consumes: Task 32's graph and TRG.
- Produces: `src/names.jl` holds `_is_reactant_copy`, `_met_label`, `_form_label`,
  `_side_label`, `_forward_sides`, `_step_text`, `_state_tag`, `_render_reaction`,
  `_class_sides`, `_constant_sides`, the `name` methods and `_step_constants`, byte for
  byte as types.jl held them. TRG pins the canonical text of the six tier-3 and tied-group
  writings, so they stay checked once the swap makes the comparison with the Step types
  circular. The pivot regression builds its mechanism from an `@enzyme_mechanism` block, not
  a literal of the old Sig.

**Gate:** full 39,768 (39,761 + 6 + 1 chokepoint-guard assertion for `names.jl`), pristine;
census against `S3b_end`, zero differences (6214bb1 moves src).

**Unit 1: f7c1587** (pins texts the graph prints; green on arrival:
S3d2_probe_tier3_text.log holds the texts at 5298b84)

- [ ] **Step 1: Apply and run**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary f7c1587 -- test/test_reaction_graph.jl | git apply
```

Background: `. .../replay/env.sh; focused T33_1_green test/test_reaction_graph.jl`. Then
`wait_log T33_1_green`, `summary T33_1_green`. Expected: TRG 812 (S3d2_c1_trg.log), `exit 0`.

- [ ] **Step 2: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Pin the canonical text of the tier-3 and tied-group cases

The graph test compares these cases with the Step types' constructor, which the
swap replaces by the graph's own; the pinned text keeps the cases checked after
that comparison stops being independent.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 2: 6214bb1** (a pure move; green before and after)

- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 6214bb1 -- src/EnzymeRates.jl src/names.jl src/types.jl | git apply
git diff --stat; git status --short
```

Expected: `src/EnzymeRates.jl | 1 +` (`include("names.jl")`), `src/types.jl | 149 +----`
(1 insertion, 148 deletions: the moved code leaves and the second ABOUTME line is
rewritten), and `?? src/names.jl` (151 lines).

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T33_2_green test/test_types.jl
test/test_rate_eq_derivation.jl`. Then `wait_log T33_2_green`, `summary T33_2_green`.

Expected: TT 664 (the chokepoint guard reads `names.jl` too), TRD 3,186, `exit 0`.

- [ ] **Step 5: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/EnzymeRates.jl src/names.jl src/types.jl
git status --short
git commit -F - <<'EOF'
Move the naming chokepoint and label renderers to names.jl

The spec's source layout gives the parameter-naming chokepoint its own file; the
graph's renderers join it at the swap. The code moves unchanged.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 3: b0a4b29** (the block builds `==` the old Sig literal, with equal fitted names
and law: S3d2_probe_pivot.log; green on arrival)

- [ ] **Step 6: Apply and run**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary b0a4b29 -- test/test_pivot_priority_regression.jl | git apply
```

Background: `. .../replay/env.sh; focused T33_3_green
test/test_pivot_priority_regression.jl`. Then `wait_log T33_3_green`, `summary
T33_3_green`. Expected: `test/test_pivot_priority_regression.jl | 2  2`, `exit 0`.

- [ ] **Step 7: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_pivot_priority_regression.jl
git status --short
git commit -F - <<'EOF'
Write the pivot regression's mechanism as an @enzyme_mechanism block

The test held its mechanism as a literal of today's Sig layout, which the swap
to the reaction graph replaces; the block builds the same mechanism (checked
equal, with equal fitted names and law, on the tree before the swap).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Task end**

- [ ] **Step 8: Full suite.** Background: `. .../replay/env.sh; full T33_full`; `wait_log
  T33_full`; `summary T33_full; pristine T33_full`. Expected: `EnzymeRates.jl | 39768
  39768`, `exit 0`; `pristine` prints `12` only.
- [ ] **Step 9: Census.** Background: `. .../replay/env.sh; census_run 33 S3b_end
  census/rename_phase3.toml`; `wait_log $P/census/logs/cmp_replay_T33.log`; `census_check
  33`. Expected: `RESULT: PASS` and nothing else.

### Task 34: The swap [S3d2]

**Replays:** `[f1aa056 + 924353b + 69964b0 + 8567c69 + 7c9dd07 + 1ecb534 + ddf94d9 + e5ef9c5
+ 050e16f + 3962124: test/test_reaction_graph.jl + 216be7f: docs/src/getting_started.md]`
plus two ROADMAP deletions, one commit (S3d2 §9 task 4 with task 7, review round 1's six
fixes of the swap's own code, and S3i's 050e16f, R58: the swap writes `_canonicalize`, so
its copies land with it; 3962124's TRG piece is S3i round-1 m6, which the S3i note folds
into 050e16f's task: the testset helper's docstring names its argument `m`).

**Files:**
- Delete: `src/reaction_graph.jl` (its contents move into `types.jl` and `names.jl`)
- Modify: `src/types.jl`, `src/names.jl`, `src/dsl.jl`, `src/mechanism_enumeration.jl`,
  `src/EnzymeRates.jl`, `ROADMAP.md`, `docs/src/getting_started.md`
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/test_compile_budget.jl`, `test/test_dsl.jl` (a comment),
  `test/test_mechanism_enumeration.jl` (comments), `test/test_rate_eq_derivation.jl`,
  `test/test_reaction_graph.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: Tasks 29-33.
- Produces: `Mechanism` and `AllostericMechanism` are the graph (table, conformations,
  forms, edges, group; tags, multiplicity, sites), defined in types.jl, the graph
  renderers in names.jl; `EnzymeMechanism{Sig}` and `AllostericEnzymeMechanism{Sig}` take
  one `Int` tuple and are exported with public docstrings; one `metabolites` generator over
  `AbstractEnzymeMechanism`; `_check_aligned`, shared by the constructor and the adapter
  `_canonical_order(rxn, groups)` (groups only); `_canonicalize` stores the table, labels
  and forms as copies at their length (050e16f). Transitional, all deleted by Task 54: the
  Step views `steps(m)`, `steps(am)` (`_view_metabolite`, `_view_species`), and today's
  signatures `Mechanism(rxn, groups)` and `AllostericMechanism(rxn, groups, tags, n, sites)`
  through `_raw_of(rxn::EnzymeReaction, groups)` and `_allosteric_of`; `_with` over Step
  groups. `_compiled_mechanism` builds `_placeholder_reaction(raw)` before the constructor,
  so the macros keep the reaction-level name checks and messages (D-S3d2-1).
  `_expand_to_allosteric` builds children with the composition constructor. Deleted: the
  Step types' canonical form and checks (`_canonical_step_direction`,
  `_orient_tied_steps`, `_step_canonical_key`, `_assert_*`), today's Sig and lifts,
  `reaction(m)`, `_step_groups`. Tests: TRG reads `const ER = EnzymeRates`; FX's
  `_testhelper_check_stored_orientation` checks that `_canonicalize(_raw(m))` turns no
  edge; TRD's textbook bridge maps written step `i` to `_canonicalize(raw).edge_of[i]` and
  errors unless `raw` canonicalizes to `m`'s catalytic graph (1ecb534); TCB's fifth fit
  pick is the ordered B-first bi-bi seed (Sig length 103); TRG pins the graph types'
  fields (924353b) and "a mechanism stores each vector at its length" (050e16f).

**Gate:** full 39,808 (39,768 − 9 + 924353b's 4 + 050e16f's 45: the TT Sig test's 13
assertions become 7, the two `_step_canonical_key` pins go and the guard reads one src file
fewer), pristine; census through the graph adapter against `S3b_end`, zero differences,
order included; both goldens unchanged and passing; seed counts 264 / 239 / 7 / 250,855;
R52's three figures recorded.

**Departures (written out):**
1. ROADMAP.md: delete the Idea "**Encode a mechanism's Sig without compiling per mix of step
   shapes.**" with the blank line after it, and the Decided-against entry "**Have the
   `Mechanism` constructors check step metabolites against the reaction.**" (216be7f's two
   hunks; the Int Sig settles the first, and a graph holds no reaction, which makes the
   second moot). Step 4 applies them as a patch.
2. The commit message: f1aa056's subject and body, then one paragraph for the folds that add
   behavior (924353b's pins, 1ecb534's guard, the Getting Started sentence), then 050e16f's
   reason in GiB (S3i round-1 m5), then one sentence for the two retired ROADMAP entries
   (Step 6). 69964b0, 8567c69, 7c9dd07, ddf94d9 and e5ef9c5 fix the swap's own wrapping,
   comments, docstrings, spacing and one helper, and 3962124's piece the argument name in
   050e16f's helper docstring, so they add no sentence. f1aa056's body
   says the graph is the one "that task 1 built", S3d's task 1; on this branch Task 1 is
   the JET task, so the body names module ReactionGraph instead (the paragraph rewraps).
3. 050e16f's TRG testset was written at S3i, on S3h's moves; here the moves still run on
   the Step views. The constructors it reaches (`_canonicalize` through `Mechanism(raw)`,
   the composition constructor, `RegulatorySite`) store every vector at its length at this
   point too, so it passes; if it fails, stop and report: do not edit it.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary f1aa056 -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_compile_budget.jl \
    test/test_dsl.jl test/test_mechanism_enumeration.jl test/test_rate_eq_derivation.jl \
    test/test_reaction_graph.jl test/test_types.jl | git apply
git show --format= --binary 924353b -- test/test_reaction_graph.jl | git apply
git show --format= --binary 69964b0 -- test/test_types.jl | git apply
git show --format= --binary 8567c69 -- test/test_types.jl | git apply
git show --format= --binary 1ecb534 -- test/test_rate_eq_derivation.jl | git apply
git show --format= --binary 050e16f -- test/test_reaction_graph.jl | git apply
git show --format= --binary 3962124 -- test/test_reaction_graph.jl | git apply
git status --short
```

Expected: no `git apply` output; the status lists the seven test files, each ` M`.

- [ ] **Step 2: Red** (the ported tests on Task 33's source)

Background: `. .../replay/env.sh; focused T34_red test/test_types.jl
test/test_reaction_graph.jl test/test_joint_solve.jl test/test_rate_eq_derivation.jl`.
Then `wait_log T34_red`, `summary T34_red`.

Expected, as the prototype's three red logs, with the folded tests added:
- TT: 1 fail and 3 errors (S3d2_red_swap.log): "_sig_of / _mechanism_from_sig roundtrip"
  (1 fail, 1 error), "a kinetic group ties one metabolite exchange at two conformations"
  and "a tied exchange turns to take up its lead's metabolites" (1 error each).
- TRG: the file stops at its first top-level read of the graph's types ("Got exception
  outside of a @test", `UndefVarError: Counts not defined in EnzymeRates`). Before that,
  the ported testsets fail or error as in S3d2_red_swap.log (TRG's 52 fails and 187
  errors there, the stop included); "a mechanism's fields are its graph and nothing else"
  fails for `Mechanism` and `AllostericMechanism` and errors for `Form` and `Edge`, which
  EnzymeRates does not define yet (S3d2_fields_red.log: the two fails); "a mechanism stores
  each vector at its length" errors for each of the 42 specs and for both populations,
  since it rebuilds a mechanism from the Step type's fields (`copy` of a reaction, a
  missing `catalytic` field), while its check that some child is allosteric passes.
- TJS: 44 errors; TRD: 45 errors (S3d2_red_swap2.log), each `UndefVarError` of `_raw` or
  `_canonicalize` in EnzymeRates, from FX's stored-orientation helper and TRD's textbook
  bridge; TRD's other testsets, Task 1's included, pass.

- [ ] **Step 3: Apply the source and the docs**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary f1aa056 -- src/EnzymeRates.jl src/dsl.jl \
    src/mechanism_enumeration.jl src/names.jl src/reaction_graph.jl src/types.jl | git apply
git show --format= --binary 7c9dd07 -- src/dsl.jl src/names.jl | git apply
git show --format= --binary ddf94d9 -- src/types.jl | git apply
git show --format= --binary e5ef9c5 -- src/types.jl | git apply
git show --format= --binary 050e16f -- src/types.jl | git apply
git show --format= --binary 216be7f -- docs/src/getting_started.md | git apply
grep -n "copy(table), copy(conformations), copy(forms)" src/types.jl
```

Expected: no `git apply` output; one `grep` line in `_canonicalize`.

- [ ] **Step 4: Apply the ROADMAP departure**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
cat > $P/replay/T34_roadmap.patch <<'PATCH'
--- a/ROADMAP.md
+++ b/ROADMAP.md
@@ -158,15 +158,6 @@
 
 ## Ideas
 
-- **Encode a mechanism's Sig without compiling per mix of step shapes.**
-  `_sig_of` (`src/types.jl`) builds the Sig tuples from generators, so Base's collect
-  widening compiles once per new mix of step shapes. Over a 6,475-mechanism bi-bi run in one
-  process (the seeds, then two expansion rounds; 5,045 allosteric) that totals about 1,500
-  method instances, 0.5 MB of native code and 1 s, 85% of it within the first 250
-  mechanisms, so a search pays it early and once per worker. Building each tuple from a
-  `Vector{Any}` should remove it. `test/test_compile_budget.jl` gates the total for the
-  239 bi-bi seeds.
-
 - **Count identifiable parameters with an exact modular rank.**
   Several phantom classes make the beam count a model above its true dimension: the
   one-subunit L, chain flanks kept after a later split, copy-copy twin pairs,
@@ -222,6 +213,3 @@
 - **Re-rescale the rate constants of saved allosteric fits.** Only CSVs written before
   0.9.0 carry the mis-scaled constants, current output is correct, and losses and model
   selection never depended on the scale.
-- **Have the `Mechanism` constructors check step metabolites against the reaction.** The
-  DSL rejects undeclared names and the enumerator builds only declared ones, so the check
-  would catch nothing and cost every enumerated mechanism.
PATCH
git apply $P/replay/T34_roadmap.patch
git diff --stat -- ROADMAP.md; git status --short
```

Expected: `ROADMAP.md | 12 ------------`; the status lists the seven test files, `M`
src/EnzymeRates.jl, src/dsl.jl, src/mechanism_enumeration.jl, src/names.jl, src/types.jl,
ROADMAP.md and docs/src/getting_started.md, and ` D src/reaction_graph.jl`. The patch is
kept: the S3d2 to S3f checkpoints apply it again.

- [ ] **Step 5: Green**

Background: `. .../replay/env.sh; focused T34_green test/test_types.jl
test/test_reaction_graph.jl test/test_joint_solve.jl test/test_rate_eq_derivation.jl
test/test_dsl.jl test/test_compile_budget.jl test/test_pivot_priority_regression.jl
test/test_mechanism_enumeration.jl`. Then `wait_log T34_green` until it ends (about 6
min), `summary T34_green`.

Expected: TT 655, TRG 861, TJS 453, TRD 3,186, TD 251, TCB 16,
`test/test_pivot_priority_regression.jl` 2, TE 32,543 (its init golden included),
`exit 0`. TRG's 861 = the prototype's 816 at 924353b (S3d2_fields_trg.log) + 050e16f's 45.

- [ ] **Step 6: Commit** (message edited: departure 2, the first paragraph included)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/EnzymeRates.jl src/dsl.jl src/mechanism_enumeration.jl src/names.jl \
    src/reaction_graph.jl src/types.jl ROADMAP.md docs/src/getting_started.md \
    test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_compile_budget.jl \
    test/test_dsl.jl test/test_mechanism_enumeration.jl test/test_rate_eq_derivation.jl \
    test/test_reaction_graph.jl test/test_types.jl
git status --short
git commit -F - <<'EOF'
Store mechanisms as reaction graphs, with Step views for unported code

Mechanism and AllostericMechanism are now the reaction graph that module
ReactionGraph built beside them: a metabolite table, forms, edges stored
physical-forward, a group per edge, and the allosteric tags, multiplicity and
sites. Their compiled types take one Int tuple as their Sig. The canonical-form
and check code of the Step types goes; the graph's canonical form and
representation checks replace it.

Code not yet ported reads a mechanism through steps(m), which builds today's
kinetic groups of Steps from the graph, and builds one through today's
constructor signatures, which turn Step groups into raw input. The macros check
the names as a reaction's before building the mechanism, as they did through
the placeholder reaction. Tests that read the old Sig, the deleted orientation
helpers or written Steps read the graph's Sig and _canonicalize instead; the
compile-budget picks take five Sig lengths.

A test pins the graph types' fields beside the Step type's, and the textbook
bridge refuses a raw value whose catalytic steps canonicalize to another
mechanism than the one it reads. Getting Started says the Sig encodes the
mechanism as a long tuple of integers.

The mechanism stores its metabolite table, conformation labels and forms at
their length. _canonicalize builds them in scratch vectors that grow past their
length (append!, unique!, push!), and a search keeps every mechanism it emits,
so spare capacity would be paid once per mechanism: the 250,855 ter-ter init
mechanisms take 1.49 GiB (6,359 B each) at grown capacity and 0.57 GiB (2,450 B
each) at their length. test/test_reaction_graph.jl pins it.

ROADMAP: the Int Sig retires the Idea "Encode a mechanism's Sig without
compiling per mix of step shapes", and a graph holds no reaction, which
retires "Have the Mechanism constructors check step metabolites against the
reaction".

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

Expected: the first status shows `D  src/reaction_graph.jl` among the staged entries: the
deletion is staged by name, in this commit (S3d2 §7: a staged `git rm` once rode into the
wrong commit).

- [ ] **Step 7: Full suite.** Background: `. .../replay/env.sh; full T34_full`; `wait_log
  T34_full`; `summary T34_full; pristine T34_full`. Expected: `EnzymeRates.jl | 39808
  39808`, `exit 0`; `pristine` prints `12` only. Aqua and JET pass: from here the suite's
  JET reads the graph in types.jl.
- [ ] **Step 8: Census through the graph adapter.** Background: `. .../replay/env.sh;
  CENSUS_ADAPTER=neutral_graph.jl census_run 34 S3b_end census/rename_phase3.toml`;
  `wait_log $P/census/logs/cmp_replay_T34.log`; `census_check 34`. Expected: `RESULT:
  PASS` and nothing else (S3d2 §5: the swap's corpus equals S3b_end line for line).
- [ ] **Step 9: Goldens.** `git log --oneline 7205ba4..HEAD -- test/fixtures` lists one
  commit, Task 32 unit 1's ("Check the reaction graph against today's canonical form"):
  neither golden changed, and both passed in Steps 5 and 7.
- [ ] **Step 10: Seed counts** (`<NN>` = 34). Expected: 264, 239, 7, 250,855.
- [ ] **Step 11: TE alone** (`<NN>` = 34). Expected: `TE_ELAPSED <s>`, `exit 0`.
- [ ] **Step 12: R52 record** (`<NN>` = 34). Report the three figures. The prototype's
  swap, without the copies, measured a 4.10 GB full-suite peak (S3d2_swapD_full) and TE
  alone 231.3 s (S3d2 §1, §7); ter-ter init took about 55 s after the swap (86d21f6's
  message) and 50.7 s before Task 35 unit 3's view hoisting (S3d2_seed_counts_wip.log).
  With the copies the peak should fall near S3b's 2.84-3.06 GB (the S3i bisect's 050e16f
  row: TE peak 2.30 GB against 3.51 GB at the swap), while the Step views still cost time
  until S3f-S3h port their readers.

### Task 35: After the swap: the bottomless check typed, views built once [S3d2]

**Replays:** `[92ba4f1]`, `[86d21f6]`, `[1894d39 + 38c9c75]` (S3d2 §9 tasks 5-6; 38c9c75 is
review round 1's m6 for the same views).

**Files:**
- Modify: `src/admission.jl`, `src/types.jl`, `src/mechanism_enumeration.jl`,
  `src/names.jl`, `src/rate_eq_derivation.jl`, `src/thermo.jl`

**Interfaces:**
- Consumes: Task 34's views and constructors.
- Produces: `_bottomless_re_segment(m::_AnyMech, steps)` lives in admission.jl beside
  `_assert_bottomed_re_segments`; the `Step` constructor finds a shared metabolite with a
  scan (`any(in(r), c)`), building the sets only for its message; `_class_sides(m)` returns
  `(; groups, kept)` and `_constant_sides` searches `sides.groups`; the moves, the emission
  rules and `_seed_variants` build the views once; `_split_gain_test` reads `ngroups(m)`;
  `_solve_tags(m)` reads `ngroups(m)`, `_inactive_law(m)` reads `length(m.edges)`;
  `show(::EnzymeMechanism)` builds the views once; `_check_well_posed(am)` builds them once
  and passes them to `_assert_one_flip_segment(am, groups)`; `_inactive_groups` builds its
  own once.

**Gate:** full 39,808, pristine; census through the graph adapter against `S3b_end`, zero
differences; ter-ter init and TE timed after unit 3 against Task 34's records (R48 m5);
checkpoint S3d2. Refactors and speed: each unit is green before and after.

**Unit 1: 92ba4f1** (move and type the bottomless check)

- [ ] **Step 1: Apply and run**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 92ba4f1 -- src/admission.jl src/types.jl | git apply
git diff --stat
```

Expected: `src/admission.jl | 36 +` and `src/types.jl | 36 -`. Background: `.
.../replay/env.sh; focused T35_1_green test/test_types.jl`; `wait_log T35_1_green`;
`summary T35_1_green`. Expected: TT 655, `exit 0`.

- [ ] **Step 2: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/admission.jl src/types.jl
git status --short
git commit -F - <<'EOF'
Type the bottomless-segment check's mechanism and keep it with the boundary check

_bottomless_re_segment reads the mechanism's substrates and products, so it
takes a mechanism, and it serves the boundary check and the moves' pre-filters,
so it lives in admission.jl beside _assert_bottomed_re_segments.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 2: 86d21f6** (the Step constructor's shared-metabolite scan)

- [ ] **Step 3: Apply and run**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 86d21f6 -- src/types.jl | git apply
```

Background: `. .../replay/env.sh; focused T35_2_green test/test_types.jl`; `wait_log
T35_2_green`; `summary T35_2_green`. Expected: TT 655 (its Step testsets pin the
shared-metabolite message), `exit 0`.

- [ ] **Step 4: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/types.jl
git status --short
git commit -F - <<'EOF'
Check a step's shared metabolites without building sets

The Step views build every step of a mechanism again on each read, and the
constructor's intersect of the consumed and released lists built two sets per
step: 6 of the 55 s the ter-ter init took after the swap. A scan finds a shared
metabolite; the sets are built only for the error message.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 3: 1894d39 + 38c9c75** (the views once per move, per naming, in `show` and in the
allosteric boundary check)

- [ ] **Step 5: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 1894d39 -- src/mechanism_enumeration.jl src/names.jl \
    src/rate_eq_derivation.jl src/thermo.jl | git apply
git show --format= --binary 38c9c75 -- src/admission.jl src/rate_eq_derivation.jl \
    src/types.jl | git apply
```

- [ ] **Step 6: Green**

Background: `. .../replay/env.sh; focused T35_3_green test/test_types.jl
test/test_rate_eq_derivation.jl test/test_reaction_graph.jl
test/test_mechanism_enumeration.jl`. Then `wait_log T35_3_green` until it ends, `summary
T35_3_green`. Expected: TT 655, TRD 3,186, TRG 861, TE 32,543, `exit 0`.

- [ ] **Step 7: Commit** (message edited: one paragraph for 38c9c75)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/admission.jl src/mechanism_enumeration.jl src/names.jl \
    src/rate_eq_derivation.jl src/thermo.jl src/types.jl
git status --short
git commit -F - <<'EOF'
Build a mechanism's Step views once per move and once per naming

Each steps(m) builds the mechanism's Steps from its graph, so a function that
read it several times, or once per candidate or per constant, now builds them
again each time. The moves, the emission rules and the seed variants read one
copy; the naming context carries the groups it searches; group and edge counts
come from the graph. The bi-bi depth-2 workload drops from 21.8 s to 18.4 s.

show(::EnzymeMechanism) builds the views once, and the allosteric boundary
check builds them twice where it built them five times: the one-flip check
takes the groups its caller built, and _inactive_groups builds its own once.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Task end**

- [ ] **Step 8: Full suite.** Background: `. .../replay/env.sh; full T35_full`; `wait_log
  T35_full`; `summary T35_full; pristine T35_full`. Expected: `EnzymeRates.jl | 39808
  39808`, `exit 0`; `pristine` prints `12` only.
- [ ] **Step 9: Census through the graph adapter.** Background: `. .../replay/env.sh;
  CENSUS_ADAPTER=neutral_graph.jl census_run 35 S3b_end census/rename_phase3.toml`;
  `wait_log $P/census/logs/cmp_replay_T35.log`; `census_check 35`. Expected: `RESULT:
  PASS` and nothing else.
- [ ] **Step 10: Timing after the views' hoisting** (R48 m5): Seed counts and TE alone
  with `<NN>` = 35, then compare `grep -h '^tt:\|TE_ELAPSED'
  $L/T34_seed_counts.log $L/T34_TE.log $L/T35_seed_counts.log $L/T35_TE.log`. Expected:
  counts 264 / 239 / 7 / 250,855 and TE `exit 0`. Report both pairs: the prototype
  measured ter-ter init about 55 s at the swap, 50.7 s before the hoisting and 51.1 s
  after, and TE 231.3 s at the swap and 189.6 s after this task (S3d2 §1); a time that
  rises here is a finding to report, not a gate.
- [ ] **Step 11: Checkpoint S3d2** (tip 38c9c75)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git diff --stat 38c9c75 HEAD -- src test docs ROADMAP.md Project.toml .claude
```

Expected: exactly these 17 files: Project.toml 2, ROADMAP.md 32,
docs/src/deriving/cha_king_altman.md 3, docs/src/deriving/mwc_allostery.md 8,
docs/src/deriving/ping_pong.md 2, docs/src/deriving/re_vs_ss.md 6,
docs/src/deriving/regulators.md 6, docs/src/developer.md 6, docs/src/getting_started.md 4
(216be7f), docs/src/identify/enumeration_engine.md 5, docs/src/identify/tutorial.md 7, the
plan (D2), the spec 48 (D1), src/types.jl 5 (050e16f's copies and comment), FX 35, TRD 31,
test/test_reaction_graph.jl 29 (050e16f's testset, its helper's docstring as 3962124 writes
it). Then:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
export GIT_INDEX_FILE=$P/replay/ckpt.idx
git read-tree 38c9c75
git diff 4f9a34a 7205ba4 -- docs/superpowers/specs ROADMAP.md | git apply --cached   # D1, D3
for c in 5c08ede 4a96d9a 0b5e637 68f46f0 33b94d8 cfe9211; do git show --format= --binary $c | git apply --cached; done   # D4
for c in 2b29fd2 be65542 0b2021d a84d194; do git show --format= --binary $c | git apply --cached; done   # docs ahead
git show --format= --binary 216be7f -- docs/src/deriving/mwc_allostery.md docs/src/deriving/cha_king_altman.md docs/src/getting_started.md | git apply --cached
git show --format= --binary 050e16f -- src/types.jl test/test_reaction_graph.jl | git apply --cached
git show --format= --binary 3962124 -- test/test_reaction_graph.jl | git apply --cached
git apply --cached $P/replay/T34_roadmap.patch
git update-index --add --cacheinfo 100644,$(git rev-parse HEAD:docs/superpowers/plans/2026-10-16-graph-mechanism.md),docs/superpowers/plans/2026-10-16-graph-mechanism.md   # D2
git diff --cached --stat HEAD -- src test docs ROADMAP.md Project.toml .claude
rm -f $GIT_INDEX_FILE; unset GIT_INDEX_FILE
```

Expected: every `git apply --cached` silent and the `git diff --cached` prints nothing.
Departures in S3d2: Task 34's ROADMAP deletions, 050e16f's copies and testset, 3962124's
docstring of that testset's helper and 216be7f's Getting Started sentence, all ahead of
their prototype position; the messages of Task 34
(its folds, 050e16f's reason in GiB, the retired entries, and module ReactionGraph named
in place of S3d's "task 1") and Task 35 unit 3, as written there.

### Task 36: Row 25: free E [S3e]

**Replays:** `[85e0790]` plus the ROADMAP entry below, one commit (S3e §9 task 1; the entry
merges 216be7f's hunk, dc4035b's scoping to allosteric mechanisms and 7205ba4's P4-1
sentences). ff5e59c, the residual-only pin, lands in Task 40: its TRG lines sit among later
S3e commits'.

**Files:**
- Modify: `src/types.jl`, `ROADMAP.md`
- Test: `test/test_mechanism_enumeration.jl`, `test/test_reaction_graph.jl`,
  `test/test_types.jl`

**Interfaces:**
- Consumes: Task 34's `_check_representation`.
- Produces: `_check_free_E(m::Mechanism)` (row 25), run after the representation rows:
  "Mechanism: no free enzyme form E; every mechanism holds the free form E, ...". Free E is
  `E`, nothing bound, no residual, form 1. TT's zero-step block is a rejection now; TE's
  `_assert_atom_conserving` fixture is a whole LDH mechanism whose one offending step is
  the transmuting isomerization.

**Gate:** full 39,812 (39,808 + TRG 3 + TE 1), pristine; census through the graph adapter
against `S3b_end`, zero differences.

**Departure (written out):** ROADMAP.md: the entry "**Take `:OnlyA` steps to their limit in
the inactive conformation.**" takes this text (7c2a657's entry with 7205ba4's P4-1
sentences, resolved; 92 columns):

```markdown
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
```

Against the branch's entry only the
sentences on free E and the boundary check change; Step 4 applies them as a patch. The
commit message adds: "ROADMAP: every mechanism holds free E, and the boundary check, not
the constructor, admits a second free form."

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 85e0790 -- test/test_mechanism_enumeration.jl \
    test/test_reaction_graph.jl test/test_types.jl | git apply
```

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T36_red test/test_reaction_graph.jl
test/test_types.jl`. Then `wait_log T36_red`, `summary T36_red`.

Expected (S3e_c1_red.log, 4 fail): TRG 861 pass, 3 fail, the no-free-E cases of "the
constructors reject what no graph represents, each with its message" (no step; a free form
F; E only with S or P bound), each "No exception thrown"; TT 654 pass, 1 fail, "Pretty
printing" (the zero-step block, "No exception thrown").

- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 85e0790 -- src/types.jl | git apply
```

- [ ] **Step 4: Apply the ROADMAP departure**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
cat > $P/replay/T36_roadmap.patch <<'PATCH'
--- a/ROADMAP.md
+++ b/ROADMAP.md
@@ -91,10 +91,11 @@
   conformation lies behind a binding, on its bound side; the limit then equals deletion
   and pruning. A ping-pong with both half-reactions `:OnlyA` fails this, so the tagging
   move will stop one half per V-type child instead (spec section 7, and the entry
-  below). The conformations will flip at free E, which every mechanism must then
-  contain; until then the constructor admits a second free form only in rapid
-  equilibrium with the first (`_assert_one_flip_segment`), where the choice only rescales
-  `L`. The rewrite's admission PR brings the limits, admission and the flip at free E.
+  below). The conformations will flip at free E, which every mechanism holds; until
+  then the boundary check admits a second free form of an allosteric mechanism only in
+  rapid equilibrium with the first (`_assert_one_flip_segment`, `src/admission.jl`),
+  where the choice only rescales `L`. The rewrite's admission PR brings the limits,
+  admission and the flip at free E.
   Branch `mwc-solve-then-limit` holds a July design of this approach (spec b067aa0 and
   651a821, plan 9c95f12, prototype 6f0011e), whose prototype matched the mass-action
   ground truth on all six uni-uni tag combinations and would delete the `d_free`
PATCH
git apply $P/replay/T36_roadmap.patch
sed -n '/^- \*\*Take `:OnlyA` steps/,/Priority: high\./p' ROADMAP.md
```

Expected: no `git apply` output; the `sed` prints the Departure's block line for line
(20 lines, from "- **Take `:OnlyA` steps to their limit in the inactive conformation.**" to
"... `test/allosteric_ground_truth.jl` is the oracle. Priority: high.").

- [ ] **Step 5: Green**

Background: `. .../replay/env.sh; focused T36_green test/test_reaction_graph.jl
test/test_types.jl test/test_mechanism_enumeration.jl`. Then `wait_log T36_green` until it
ends, `summary T36_green`. Expected: TRG 864, TT 655, TE 32,544 (S3e_c1_green.log: 34,018
over the three), `exit 0`.

- [ ] **Step 6: Commit** (message edited: the ROADMAP sentence)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/types.jl ROADMAP.md test/test_mechanism_enumeration.jl \
    test/test_reaction_graph.jl test/test_types.jl
git status --short
git commit -F - <<'EOF'
Reject a mechanism without free E

Every mechanism holds the free form E, with nothing bound and no residual, as form 1
(spec section 4, R9, R28): the enumerator labels free E so, and a mechanism without it
has no form the conformations flip at and no reference for the dead-end placements. The
constructor now checks it after the representation rows. A block with no steps, which
the show test printed as a bare header, is rejected; the atom-conservation fixture
becomes a whole LDH mechanism whose one offending step is the transmuting isomerization.

ROADMAP: every mechanism holds free E, and the boundary check, not the constructor,
admits a second free form.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

- [ ] **Step 7: Full suite.** Background: `. .../replay/env.sh; full T36_full`; `wait_log
  T36_full`; `summary T36_full; pristine T36_full`. Expected: `EnzymeRates.jl | 39812
  39812`, `exit 0`; `pristine` prints `12` only.
- [ ] **Step 8: Census through the graph adapter.** Background: `. .../replay/env.sh;
  CENSUS_ADAPTER=neutral_graph.jl census_run 36 S3b_end census/rename_phase3.toml`;
  `wait_log $P/census/logs/cmp_replay_T36.log`; `census_check 36`. Expected: `RESULT:
  PASS` and nothing else (every enumerated mechanism holds free E).

### Task 37: Rows 26-28: one connected graph, every declared name used [S3e]

**Replays:** `[b9c535e]`, `[450a713 + 86d0c31 + 216be7f: src/dsl.jl]` (S3e §9 tasks 2-3;
86d0c31 is review round 1's B1, which restores the two test-side failure pins; 216be7f's
`@enzyme_mechanism` docstring states the three rules once all three hold).

**Files:**
- Modify: `src/types.jl`, `src/dsl.jl` (the `@enzyme_mechanism` docstring)
- Test: `test/test_reaction_graph.jl`, `test/test_oracles.jl`,
  `test/test_mechanism_enumeration.jl`, `test/test_types.jl`,
  `test/test_identify_rate_equation.jl`

**Interfaces:**
- Consumes: Task 36's `_check_free_E`; `$P/replay/mutate.py` (Task 0 Step 8):
  `invariants_blind`, `oracle_blind`.
- Produces: `_check_connected(m)` (row 26: "Mechanism: the steps do not form one connected
  graph; no steps lead from free E to EstarP"), `_check_declared_used(m, modifiers)` (rows
  27-28: "Mechanism: <kind> X is declared, but no step uses it", "Mechanism: modifier X is
  declared, but no step binds it"), both after `_check_free_E`. TO pins the oracle's own
  disconnected-network error and TE the invariant helper's unused-substrate error, each on
  `_canonicalize(_raw_mechanism(block, false)).mechanism`, beside the macro's rejection.
  Fixtures that broke a rule use every name and join every form: TE's
  representation-independence blocks, TT's mirrored isomerization and bottomless-segment
  fixtures, TE's conformation context (binds B at Estar), TID's three atom-breaking
  mechanisms (bind T in a dead end at free E).

**Gate:** full 39,821 (39,812 + 1 + 5 + 86d0c31's 3), pristine; census through the graph
adapter against `S3b_end`, zero differences.

**Unit 1: b9c535e** (row 26)

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary b9c535e -- test/test_mechanism_enumeration.jl \
    test/test_oracles.jl test/test_reaction_graph.jl test/test_types.jl | git apply
```

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T37_1_red test/test_reaction_graph.jl
test/test_oracles.jl`. Then `wait_log T37_1_red`, `summary T37_1_red`.

Expected (S3e_c2_red.log, 2 fail): TRG 864 pass, 1 fail (the `EstarP` island, "No
exception thrown"); TO 661 pass, 1 fail ("a disconnected network is no mechanism": the
two-enzyme network builds).

- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary b9c535e -- src/types.jl | git apply
```

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T37_1_green test/test_reaction_graph.jl
test/test_oracles.jl test/test_types.jl test/test_mechanism_enumeration.jl`. Then
`wait_log T37_1_green` until it ends, `summary T37_1_green`. Expected: TRG 865, TO 662, TT
655, TE 32,544 (S3e_c2_green.log: 34,681 over the four), `exit 0`.

- [ ] **Step 5: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/types.jl test/test_mechanism_enumeration.jl test/test_oracles.jl \
    test/test_reaction_graph.jl test/test_types.jl
git status --short
git commit -F - <<'EOF'
Reject steps that do not form one connected graph

A form no steps join to free E sits in a part of the network with its own steady
state, so the enzyme's split between the parts, and the rate, is undefined (spec
section 4). The constructor now rejects such a network and names the first form free E
does not reach. The oracle's two-enzyme network becomes a rejection test; the mirrored
isomerization and canonical-construction fixtures gain the steps that join their
forms, written as mechanism blocks.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 2: 450a713 + 86d0c31 + 216be7f: src/dsl.jl** (rows 27-28; the failure pins; the
docstring)

- [ ] **Step 6: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 450a713 -- test/test_identify_rate_equation.jl \
    test/test_mechanism_enumeration.jl test/test_reaction_graph.jl \
    test/test_types.jl | git apply
git show --format= --binary 86d0c31 -- test/test_mechanism_enumeration.jl \
    test/test_oracles.jl | git apply
```

- [ ] **Step 7: Red**

Background: `. .../replay/env.sh; focused T37_2_red test/test_reaction_graph.jl
test/test_mechanism_enumeration.jl test/test_identify_rate_equation.jl test/test_types.jl
test/test_oracles.jl`. Then `wait_log T37_2_red` until it ends, `summary T37_2_red`.

Expected (S3e_c3_red.log, 6 fail): TRG 865 pass, 5 fail (the unused substrate, product
and modifier, a substrate exchanged only as its copy, a reactant modifier no copy binds,
each "No exception thrown"); TE 32,544 pass, 1 fail, in
"_testhelper_assert_mechanism_invariants: coverage of every reactant" (the macro builds
the block that declares T); TID 446, TT 655 and TO 664 pass. 86d0c31's three pins pass
here (TE's helper error, TO's oracle error concrete and compiled): they fail only on
their mutants, which Step 10 runs (S3e_r1_red_mutants.log: TO:307 twice
`SingularException(4)`, TE:408 "No exception thrown").

- [ ] **Step 8: Apply the source and the docstring**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 450a713 -- src/types.jl | git apply
git show --format= --binary 216be7f -- src/dsl.jl | git apply
grep -n "Every mechanism holds the free enzyme" src/dsl.jl
```

Expected: no `git apply` output; one `grep` line in the `@enzyme_mechanism` docstring.

- [ ] **Step 9: Green**

Background: `. .../replay/env.sh; focused T37_2_green test/test_reaction_graph.jl
test/test_mechanism_enumeration.jl test/test_identify_rate_equation.jl test/test_types.jl
test/test_oracles.jl`. Then `wait_log T37_2_green` until it ends, `summary T37_2_green`.
Expected: TRG 870, TE 32,545, TID 446, TT 655, TO 664, `exit 0`.

- [ ] **Step 10: Red on the two test-side mutants** (S3e §9 tasks 2-3: the invariant
  helper and the oracle each lose their check)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py invariants_blind
python3 $P/replay/mutate.py oracle_blind
```

Expected: `mutated test/mechanism_definitions_for_test_enzyme_derivation.jl
(invariants_blind)` and the same for `oracle_blind`. Background: `. .../replay/env.sh;
focused T37_2_mut test/test_oracles.jl test/test_mechanism_enumeration.jl`; `wait_log
T37_2_mut` until it ends; `summary T37_2_mut`. Expected
(`$P/proto/logs/S3e_r1_red_mutants.log`): TO 662 pass, 2 fail, both in "rejects a
disconnected network" (concrete and compiled: thrown `SingularException(4)`, not the
oracle's disconnected-network error); TE 32,544 pass, 1 fail, in
"_testhelper_assert_mechanism_invariants: coverage of every reactant" ("No exception
thrown"). Then, whatever the run printed:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py oracle_blind --revert
python3 $P/replay/mutate.py invariants_blind --revert
git diff --stat -- test/mechanism_definitions_for_test_enzyme_derivation.jl
```

Expected: `restored ...` twice and no `git diff` output (the unit leaves FX alone).

- [ ] **Step 11: Commit** (message edited: one paragraph for 86d0c31 and the docstring)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/dsl.jl src/types.jl test/test_identify_rate_equation.jl \
    test/test_mechanism_enumeration.jl test/test_oracles.jl test/test_reaction_graph.jl \
    test/test_types.jl
git status --short
git commit -F - <<'EOF'
Reject a declared name that no step exchanges

A substrate or product that no step takes up or gives off, or a declared modifier that
no step binds as a modifier, takes no part in the rate law (spec section 4, R9). The
constructor now rejects each and names it; a reactant exchanged only as its modifier
copy counts as unused. Fixtures that declared such a name now use it: the identify
fixtures whose chemistry breaks atom balance bind T in a dead end at free E, the
bottomless-segment Step fixture gains P's release and the chemistry, and the
conformation-context fixture binds B. The invariant helper's negative case becomes the
constructor's rejection.

A mechanism the constructors reject still reaches the invariant helper and the oracle
as the canonical graph of the walker's raw value, which no constructor check reads, so
each testset keeps the macro's rejection and pins the check's own error on that graph:
the helper's unused-substrate message, and the oracle's disconnected-network message,
concrete and compiled. The @enzyme_mechanism docstring states the rules: every
mechanism holds free E, its steps join every form to it, and each declared name takes
part in some step.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Task end**

- [ ] **Step 12: Full suite.** Background: `. .../replay/env.sh; full T37_full`; `wait_log
  T37_full`; `summary T37_full; pristine T37_full`. Expected: `EnzymeRates.jl | 39821
  39821`, `exit 0`; `pristine` prints `12` only.
- [ ] **Step 13: Census through the graph adapter.** Background: `. .../replay/env.sh;
  CENSUS_ADAPTER=neutral_graph.jl census_run 37 S3b_end census/rename_phase3.toml`;
  `wait_log $P/census/logs/cmp_replay_T37.log`; `census_check 37`. Expected: `RESULT:
  PASS` and nothing else.

### Task 38: Rows 29-31: names and labels in the constructors, R38 at the macros [S3e]

**Replays:** `[a76ddf2 + 9221bfc + 3bb8d94 + 4d2af5a: src/dsl.jl test/test_types.jl]`,
`[2aa3cbd + 4d2af5a: test/test_reaction_graph.jl + 8c867d3]` (S3e §9 tasks 4-5; 9221bfc and
3bb8d94 are review round 1's m1 and m4, pins of a76ddf2's checks; 8c867d3 is its m3, which
R53 m2 and §9 task 5 fold into 2aa3cbd, "shown failing on the one-message version" in Step
9; 4d2af5a wraps at 92 columns four lines that a76ddf2, 2aa3cbd and ad74058 write, and §9
folds the wraps into those tasks, so its pieces land with them, the TD line in Task 39).
8c867d3's TRG pin sits on the label-collision pin as 4d2af5a wraps it, so the wrap of that
pin comes first.

**Files:**
- Modify: `src/types.jl`, `src/dsl.jl`, `src/mechanism_enumeration.jl` (the
  `_assert_atom_conserving` docstring)
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/test_rate_eq_derivation.jl`, `test/test_reaction_graph.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: Tasks 34-37's constructors; `_check_metabolite_names` (Task 8),
  `_check_uniquely_decodable` (Task 9).
- Produces: `_KEYWORDS`; `_check_metabolite_names(who, names::Vector{Symbol})` (typed for
  JET) also rejects `E` and any name that is no identifier or is a keyword, and runs in
  `EnzymeReaction`, `_check_declarations` (row 29) and `_check_ligands`;
  `_check_uniquely_decodable(reactants, modifiers, who)` and `_spelled_twice_error(a, b,
  who)` take the error's prefix; `_compiled_mechanism` calls
  `_check_uniquely_decodable(..., "Mechanism")` before the constructor, so
  `_placeholder_reaction` leaves src and FX's `_testhelper_placeholder_reaction` builds its
  own reaction; `_check_conformation_label(c, who)` and `_label_collision_error(c, who,
  kind)` (`kind`: `"declared metabolite"` for the walker and the catalytic check, `"ligand"`
  for `_check_ligands`; 8c867d3) move to types.jl, and `_check_representation` (rows 30-31,
  first) and `_check_ligands(sites, catalytic)` call them. Pins move with their checks:
  TRD's `k_EP_to_ES` and TT's `P::Mod`/`Pmod` now read "Mechanism:".

**Gate:** full 39,846 (39,821 + 13 + 1 + 7 + 4), pristine (JET included: a76ddf2's first
cut failed JET on an untyped `names`); census through the graph adapter against `S3b_end`,
zero differences.

**Unit 1: a76ddf2 + 9221bfc + 3bb8d94 + 4d2af5a's src/dsl.jl and TT** (row 29, R38 at the
macros, their pins, and the wraps of two lines a76ddf2 writes)

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary a76ddf2 -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_rate_eq_derivation.jl test/test_reaction_graph.jl \
    test/test_types.jl | git apply
git show --format= --binary 9221bfc -- test/test_types.jl | git apply
git show --format= --binary 3bb8d94 -- test/test_types.jl | git apply
git show --format= --binary 4d2af5a -- test/test_types.jl | git apply
```

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T38_1_red test/test_reaction_graph.jl
test/test_types.jl test/test_rate_eq_derivation.jl`. Then `wait_log T38_1_red`, `summary
T38_1_red`.

Expected (S3e_c4_red.log's 15, plus 9221bfc's 1): TRG 870 pass, 7 fail (the five
declared-name cases, which report row 27 first, "Thrown: Mechanism: modifier 1I is
declared, but no step binds it", and the two ligand cases, "No exception thrown"); TT 661
pass, 8 fail ("EnzymeReaction validation" 6: `E` by the constructor and by the macro,
`1A`, `A B`, `end`, `in`; "two forms whose constants would share a name are rejected" 2:
the `P::Mod`/`Pmod` pin and 9221bfc's allosteric block, both thrown with the placeholder
reaction's `EnzymeReaction:` prefix); TRD 3,185 pass, 1 fail ("a metabolite named after a
constant is rejected": `k_EP_to_ES` gets `EnzymeReaction:`). 3bb8d94's seven
contextual-keyword round trips pass: they fail only on a keyword list that names them
(S3e_r1_red_mutants.log, TT:1023).

- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary a76ddf2 -- src/dsl.jl src/mechanism_enumeration.jl \
    src/types.jl | git apply
git show --format= --binary 4d2af5a -- src/dsl.jl | git apply
! grep -n "_placeholder_reaction" src/*.jl
```

Expected: no output (`_placeholder_reaction` left src with its one call).

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T38_1_green test/test_reaction_graph.jl
test/test_types.jl test/test_rate_eq_derivation.jl test/test_dsl.jl
test/test_accessors.jl`. Then `wait_log T38_1_green`, `summary T38_1_green`. Expected:
TRG 877, TT 669, TRD 3,186, TD 251, TA 155 (S3e_c4_green.log's five files), `exit 0`.

- [ ] **Step 5: Commit** (message edited: one sentence for the two pins; 4d2af5a's wraps add
  none)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/dsl.jl src/mechanism_enumeration.jl src/types.jl \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_rate_eq_derivation.jl test/test_reaction_graph.jl test/test_types.jl
git status --short
git commit -F - <<'EOF'
Check metabolite names in the constructors, decodability at the macros

The constructor now rejects a declared name that is no identifier, a keyword, E, a
reserved name or one holding `_` (spec section 4, R10, R26, row 29), and so does the
reaction: E labels the free enzyme form, so a metabolite named E would make every
enumerated form read as that metabolite, and a name no DSL text can write would print
a block that does not read back (R46 m1, m4). The ligand check shares the rule.

The mechanism macros no longer build a placeholder reaction: the constructor checks the
names, and the macros check only that the names spell every string one way (R38),
before the constructor, with the mechanism's prefix. The placeholder reaction moves to
the test helper that rebuilds mechanisms over Step groups, its only user.

Tests pin the allosteric macro's decodability check, and the contextual keywords
(abstract, mutable, outer, primitive, public, type, var) as substrate names whose
reaction text reads back equal.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 2: 2aa3cbd + 4d2af5a's TRG + 8c867d3** (rows 30-31; the wrap of 2aa3cbd's
label-collision pin; a ligand that a label names is called a ligand)

- [ ] **Step 6: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2aa3cbd -- test/test_reaction_graph.jl | git apply
git show --format= --binary 4d2af5a -- test/test_reaction_graph.jl | git apply
git show --format= --binary 8c867d3 -- test/test_reaction_graph.jl | git apply
```

- [ ] **Step 7: Red**

Background: `. .../replay/env.sh; focused T38_2_red test/test_reaction_graph.jl`. Then
`wait_log T38_2_red`, `summary T38_2_red`.

Expected (S3e_c5_red.log): TRG 877 pass, 4 fail: the label shape `Ex_1`, the opaque label
`ES` (row 21's message comes first), a label `S` naming a metabolite and a label `F`
naming a ligand, each "No exception thrown" but `ES`. 8c867d3's re-pin of the last case
(`with ligand \`F\``) fails the same way here, so the count is the log's.

- [ ] **Step 8: Apply 2aa3cbd's source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2aa3cbd -- src/dsl.jl src/types.jl | git apply
```

- [ ] **Step 9: Red for 8c867d3** (the one-message version, as S3e §9 task 5 asks)

Background: `. .../replay/env.sh; focused T38_2_red2 test/test_reaction_graph.jl`. Then
`wait_log T38_2_red2`, `summary T38_2_red2`.

Expected (`$P/proto/logs/S3e_r1_red_mutants.log`, TRG:752): TRG 880 pass, 1 fail, in "the
constructors reject what no graph represents, each with its message": the label `F` naming
a ligand throws "AllostericMechanism: conformation label `F` collides with declared
metabolite `F`; choose a different conformation label." where the pin reads "with ligand
`F`".

- [ ] **Step 10: Apply 8c867d3's source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 8c867d3 -- src/dsl.jl src/types.jl | git apply
```

- [ ] **Step 11: Green**

Background: `. .../replay/env.sh; focused T38_2_green test/test_reaction_graph.jl
test/test_types.jl test/test_dsl.jl`. Then `wait_log T38_2_green`, `summary T38_2_green`.
Expected: TRG 881, TT 669, TD 251 (the walker's label messages), `exit 0`.

- [ ] **Step 12: Commit** (message edited: one paragraph for 8c867d3; the wrap adds none)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/dsl.jl src/types.jl test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Check conformation labels in the constructors

A conformation label must have a label's shape and name no metabolite, so `show`
prints every mechanism the constructors accept as text the DSL reads back (spec section
5, R26, R40, R46 m4; rows 30 and 31). The walker checked both for the macros; the
constructor now checks each label of the graph and the composition constructor each
ligand, with the walker's messages: the label rule and the collision message move to
types.jl, and the walker calls them.

The collision error takes what the name is: a declared metabolite (the walker and
the catalytic check) or a ligand (the composition constructor), so a site ligand F
that a label also names is reported as a ligand.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Task end**

- [ ] **Step 13: Full suite.** Background: `. .../replay/env.sh; full T38_full`; `wait_log
  T38_full`; `summary T38_full; pristine T38_full`. Expected: `EnzymeRates.jl | 39846
  39846`, `exit 0`; `pristine` prints `12` only.
- [ ] **Step 14: Census through the graph adapter.** Background: `. .../replay/env.sh;
  CENSUS_ADAPTER=neutral_graph.jl census_run 38 S3b_end census/rename_phase3.toml`;
  `wait_log $P/census/logs/cmp_replay_T38.log`; `census_check 38`. Expected: `RESULT:
  PASS` and nothing else (8c867d3 changes one message, which no census mechanism meets).

### Task 39: Multiplicities, typed inhibitors, elements, Bool [S3e]

**Replays:** `[ad74058 + 4d2af5a: test/test_dsl.jl]`, `[e7698be]`, `[56dbf00]`, `[5b72e27]`,
each its own commit (S3e §9 tasks 6-9; 4d2af5a's TD piece wraps the multiplicity pin
ad74058 writes, as §9 folds the wraps into tasks 4-6, and adds no sentence).

**Files:**
- Modify: `src/types.jl`, `src/dsl.jl`
- Test: `test/test_dsl.jl`, `test/test_mechanism_enumeration.jl` (`rxn_ci`),
  `test/test_reaction_graph.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: `RegulatorMults`, `EnzymeReaction`, `ReactantAtoms`, `RawSite`,
  `RawAllostericMechanism`, the walker's `_parse_multiplicity_tuple`, `_positive_int`,
  `_parse_regulator_entries`.
- Produces: `_check_multiplicities(ms, what)` (nonempty, each ≥ 1, distinct), called by
  `RegulatorMults` and `EnzymeReaction`, and the walker rejects a repeated
  `oligomeric_state:` entry at expansion; `RegulatorMults` rejects a regulator type on a
  `CompetitiveInhibitor` (the walker's own check goes; the macro reports it at run time,
  D-S3e-6); `ReactantAtoms` rejects an element that is not a capital letter and any
  lowercase letters; `_positive_int` rejects `Bool`, and `RawSite` and
  `RawAllostericMechanism` gain typed inner constructors.

**Gate:** full 39,866 (39,846 + 5 + 2 + 3 + 10), pristine; census through the graph adapter
against `S3b_end`.

**Unit 1: ad74058** (an empty or repeated multiplicity list)

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary ad74058 -- test/test_dsl.jl test/test_types.jl | git apply
git show --format= --binary 4d2af5a -- test/test_dsl.jl | git apply
```

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T39_1_red test/test_types.jl test/test_dsl.jl`.
Then `wait_log T39_1_red`, `summary T39_1_red`.

Expected (S3e_c6_red.log, 7 fail): TT 667 pass, 5 fail ("RegulatorMults canonicalizes
ordering + validates": the empty and `[2, 2]` lists; "EnzymeReaction takes one or more
catalytic multiplicities, each ≥ 1": `[0, 1]` and `Int[]`, whose message now says
"distinct", and
`[2, 2]`, accepted); TD 251 pass, 2 fail ("@enzyme_reaction sets the
catalytic multiplicities once": `oligomeric_state: (2, 2)` at expansion, `A(2, 2)` at run
time).

- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary ad74058 -- src/dsl.jl src/types.jl | git apply
```

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T39_1_green test/test_types.jl
test/test_dsl.jl`. Then `wait_log T39_1_green`, `summary T39_1_green`. Expected: TT 672,
TD 253, `exit 0`.

- [ ] **Step 5: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/dsl.jl src/types.jl test/test_dsl.jl test/test_types.jl
git status --short
git commit -F - <<'EOF'
Reject an empty or repeated multiplicity list

`oligomeric_state: (2, 2)` stored [2] on the reaction but gave each bare regulator
[2, 2], so the reaction differed from `oligomeric_state: 2` and printed a regulator list
no DSL text means (R46 m2). The walker now rejects a repeated catalytic multiplicity at
expansion, and the reaction and RegulatorMults reject what no DSL text writes: an empty
or repeated multiplicity list, through one check (R46 m1, R13).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 2: e7698be** (a regulator type on a competitive inhibitor)

- [ ] **Step 6: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary e7698be -- test/test_dsl.jl \
    test/test_mechanism_enumeration.jl test/test_types.jl | git apply
```

- [ ] **Step 7: Red**

Background: `. .../replay/env.sh; focused T39_2_red test/test_types.jl test/test_dsl.jl`.
Then `wait_log T39_2_red`, `summary T39_2_red`.

Expected (S3e_c7_red.log, 4 fail): TT 672 pass, 2 fail ("RegulatorMults reg_type": a typed
competitive inhibitor builds); TD 251 pass, 2 fail ("@enzyme_reaction regulator types":
`X::Activator` and `X::Inhibitor` on a competitive inhibitor get the walker's expansion
error, not RegulatorMults' message).

- [ ] **Step 8: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary e7698be -- src/dsl.jl src/types.jl | git apply
```

- [ ] **Step 9: Green**

Background: `. .../replay/env.sh; focused T39_2_green test/test_types.jl test/test_dsl.jl
test/test_mechanism_enumeration.jl`. Then `wait_log T39_2_green` until it ends, `summary
T39_2_green`. Expected: TT 674, TD 253, TE 32,545 (`rxn_ci` declares its inhibitor
untyped through the macro), `exit 0`.

- [ ] **Step 10: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/dsl.jl src/types.jl test/test_dsl.jl test/test_mechanism_enumeration.jl \
    test/test_types.jl
git status --short
git commit -F - <<'EOF'
Reject a regulator type on a competitive inhibitor

Only an allosteric regulator activates or inhibits through the conformation it
favors, so the walker refused `competitive_inhibitors: X::Activator`, while a reaction
built directly could hold a typed competitive inhibitor, which `show` printed as text
the macro rejects (R46 m1, R13). RegulatorMults now rejects the type, and the walker's
own check goes, so one check covers both paths; the macro reports it when its code
runs. The enumeration test that built such a reaction directly declares the inhibitor
untyped through the macro instead.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 3: 56dbf00** (an element no formula writes)

- [ ] **Step 11: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 56dbf00 -- test/test_types.jl | git apply
```

- [ ] **Step 12: Red**

Background: `. .../replay/env.sh; focused T39_3_red test/test_types.jl`. Then `wait_log
T39_3_red`, `summary T39_3_red`.

Expected (S3e_c8_red.log): TT 674 pass, 3 fail ("ReactantAtoms validation": the elements
`c`, `CL` and `C1` are accepted).

- [ ] **Step 13: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 56dbf00 -- src/types.jl | git apply
```

- [ ] **Step 14: Green**

Background: `. .../replay/env.sh; focused T39_3_green test/test_types.jl
test/test_dsl.jl`. Then `wait_log T39_3_green`, `summary T39_3_green`. Expected: TT 677,
TD 253 (S3e_c8_green.log: 922 over the two), `exit 0`.

- [ ] **Step 15: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/types.jl test/test_types.jl
git status --short
git commit -F - <<'EOF'
Reject an atom element that no formula writes

A formula in `@enzyme_reaction` writes each element as a capital letter and any
lowercase letters, but a reaction built directly could hold an element such as `c` or
`C1`, which `show` printed as a formula the macro rejects or reads as other elements
(R46 m1, R13). ReactantAtoms now rejects such an element and names it.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 4: 5b72e27** (Bool multiplicities)

- [ ] **Step 16: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 5b72e27 -- test/test_dsl.jl test/test_reaction_graph.jl \
    test/test_types.jl | git apply
```

- [ ] **Step 17: Red**

Background: `. .../replay/env.sh; focused T39_4_red test/test_reaction_graph.jl
test/test_dsl.jl test/test_types.jl`. Then `wait_log T39_4_red`, `summary T39_4_red`.

Expected (S3e_c9_red.log, 6 fail): TRG 883 pass, 2 fail (`RawSite` and
`RawAllostericMechanism` convert `true` to 1, "No exception thrown"); TD 253 pass, 4 fail
("chemical formulas and multiplicities are validated": `A(true)`, `oligomeric_state:
true`, `catalytic_multiplicity: true`, `regulatory_site(true)`); TT 679 pass (its two Bool
pins hold already: typed signatures).

- [ ] **Step 18: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 5b72e27 -- src/dsl.jl src/types.jl | git apply
```

- [ ] **Step 19: Green**

Background: `. .../replay/env.sh; focused T39_4_green test/test_reaction_graph.jl
test/test_dsl.jl test/test_types.jl`. Then `wait_log T39_4_green`, `summary T39_4_green`.
Expected: TRG 885, TD 257, TT 679, `exit 0`.

- [ ] **Step 20: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/dsl.jl src/types.jl test/test_dsl.jl test/test_reaction_graph.jl \
    test/test_types.jl
git status --short
git commit -F - <<'EOF'
Reject Bool multiplicities

Julia counts `true` an Integer, so `catalytic_multiplicity: true`,
`regulatory_site(true)`, `oligomeric_state: true` and `A(true)` built multiplicity 1,
and a raw site or allosteric mechanism built directly converted `true` to 1 (R46 m7).
`_positive_int` now rejects a Bool, and the two raw types take an Int without
converting, as RegulatorySite, AllostericMechanism, RegulatorMults and EnzymeReaction
already do; tests pin all six.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Task end**

- [ ] **Step 21: Full suite.** Background: `. .../replay/env.sh; full T39_full`; `wait_log
  T39_full`; `summary T39_full; pristine T39_full`. Expected: `EnzymeRates.jl | 39866
  39866`, `exit 0`; `pristine` prints `12` only.
- [ ] **Step 22: Census through the graph adapter.** Background: `. .../replay/env.sh;
  CENSUS_ADAPTER=neutral_graph.jl census_run 39 S3b_end census/rename_phase3.toml`;
  `wait_log $P/census/logs/cmp_replay_T39.log`; `census_check 39`. Expected: `RESULT:
  PASS` and nothing else.

### Task 40: Test cleanups, tighter Sig gates, the review's pins [S3e]

**Replays:** `[c6d5a47]`, `[9cbb9af]`, `[44661cc]`, `[dfc5f79 + 88f896e]`, `[ff5e59c]` (S3e
§9 tasks 10-13 and review round 1's B2 and m2; 4d2af5a's wraps and 8c867d3, round 1's m3,
are in Tasks 38-39, where §9 puts them). ff5e59c is §9 task 1's residual-only case (R53 m2),
which lands here because its TRG lines sit among later S3e commits' (Task 36); Step 15 shows
it failing on the residual-blind check that task 1's Start line names.

**Files:**
- Modify: `src/types.jl` (the lifts' comment)
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/test_compile_budget.jl`, `test/test_identify_rate_equation.jl`,
  `test/test_mechanism_enumeration.jl`, `test/test_oracles.jl`,
  `test/test_rate_eq_derivation.jl`, `test/test_reaction_graph.jl`

**Interfaces:**
- Consumes: Tasks 36-39; `$P/replay/mutate.py` (Task 0 Step 8): `occupied`,
  `sig_per_mechanism`, `residual_blind`.
- Produces: fixtures built from the macros alone (no reaction wrappers: TE `withrxn`, TE
  `_testhelper_on_rxn`, TID `attach`); FX's `_testhelper_regulators` keeps its plain method
  only; FX's `_testhelper_graph` names inactive constants only for steps the inactive law
  keeps, and TRD's `_testhelper_written_params` skips empty names; FX's spec struct takes a
  raw value always (no `Nothing`, no default); TRG's row-18 case with unequal copies, the
  benchmark comment citing the spec's 2 µs goal, the residual-only free-E case; TCB's Sig
  budgets 250 instances and 100 KB.

**Gate:** full 39,871 (39,866 + 3 + 1 + 1), pristine; census through the graph adapter
against `S3b_end`, zero differences; both goldens unchanged; seed counts; the
constructor's median over the specs, about 3.7 µs (S3e_measure.log: 3.66); R52's figures;
checkpoint S3e.

**Unit 1: c6d5a47** (fixtures without reaction wrappers; test-only, green before and
after)

- [ ] **Step 1: Apply and run**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary c6d5a47 -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_identify_rate_equation.jl test/test_mechanism_enumeration.jl | git apply
```

Background: `. .../replay/env.sh; focused T40_1_green test/test_mechanism_enumeration.jl
test/test_identify_rate_equation.jl`. Then `wait_log T40_1_green` until it ends, `summary
T40_1_green`. Expected: TE 32,545, TID 446, `exit 0`.

- [ ] **Step 2: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_identify_rate_equation.jl test/test_mechanism_enumeration.jl
git status --short
git commit -F - <<'EOF'
Build test fixtures without reaction wrappers and drop an unread helper method

A mechanism's identity leaves out its reaction, so the three local wrappers that rebuilt
a macro's mechanism over the test's reaction (two enumeration testsets, one identify
testset) build the same mechanism the macro does; the fixtures now come from the macro
alone (R49 m4). No test reads the allosteric method of `_testhelper_regulators`, which
goes with its half of the comment (R49 m1). The invariant helper's comment said init
mechanisms declare unbound inhibitors, which the constructor now rejects; it says what
the helper checks.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 2: 9cbb9af** (the reader's inactive names; a Theorell–Chance pin). The reader is test
code, so the red applies TO's pins first and the reader after them, as S3e_c11_red.log
ran them against the old reader.

- [ ] **Step 3: Apply TO's pins**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9cbb9af -- test/test_oracles.jl | git apply
```

- [ ] **Step 4: Red**

Background: `. .../replay/env.sh; focused T40_2_red test/test_oracles.jl`. Then `wait_log
T40_2_red`, `summary T40_2_red`.

Expected (S3e_c11_red.log): TO 666 pass, 1 fail ("the test helpers read a mechanism's
forms, steps, tags and sites": the three re-pinned steps' `inactive_names` read
`[:K_I_EI_to_E_I]`, `[:K_I_EP_to_E_P]` and `[:k_I_ES_to_EP, :k_I_EP_to_ES]` where the pin
expects empty). The Theorell–Chance pin passes on the old reader.

- [ ] **Step 5: Apply the reader and the bridge**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9cbb9af -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_rate_eq_derivation.jl | git apply
```

- [ ] **Step 6: Green**

Background: `. .../replay/env.sh; focused T40_2_green test/test_oracles.jl
test/test_rate_eq_derivation.jl`. Then `wait_log T40_2_green`, `summary T40_2_green`.
Expected: TO 667, TRD 3,186 (the 42 bridge checks unchanged), `exit 0`.

- [ ] **Step 7: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_oracles.jl \
    test/test_rate_eq_derivation.jl
git status --short
git commit -F - <<'EOF'
Name only kept steps' inactive constants in the test reader, and pin a Theorell–Chance step

`_testhelper_graph` named every step's constants in the inactive conformation, even for
steps the inactive law does not keep, which have none; those steps now read no inactive
names, and the textbook bridge skips them where it skipped names the law lacks (R49 m6).
The reader's testset re-pins the three such steps of its allosteric mechanism and pins a
Theorell–Chance mechanism, the first whose step gives off a metabolite, so `released`
carries a name (R49 m5).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 3: 44661cc** (row 18 with unequal copies; the raw field; green on arrival, red only
on the `_occupied`-free mutant: S3e_c12_mut_occupied.log, TRG:503 "No exception thrown")

- [ ] **Step 8: Apply and run**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 44661cc -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_reaction_graph.jl | git apply
```

Background: `. .../replay/env.sh; focused T40_3_green test/test_reaction_graph.jl`. Then
`wait_log T40_3_green`, `summary T40_3_green`. Expected: TRG 886, `exit 0`.

- [ ] **Step 9: Red on the `_occupied`-free mutant** (S3e §9 task 12)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py occupied
```

Expected: `mutated src/types.jl (occupied)`. Background: `. .../replay/env.sh; focused
T40_3_mut test/test_reaction_graph.jl`; `wait_log T40_3_mut`; `summary T40_3_mut`.
Expected (`$P/proto/logs/S3e_c12_mut_occupied.log`, TRG:503): TRG 885 pass, 1 fail, in
"the constructors reject what no graph represents, each with its message": row 18's
unequal-copies step, "No exception thrown". Then, whatever the run printed:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py occupied --revert
git status --short
```

Expected: `restored src/types.jl`, then ` M` FX and TRG only.

- [ ] **Step 10: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Pin a step that shares a metabolite in unequal copies, and require each spec's raw value

Row 18 rejects a step whose sides share a metabolite, which a plain bitwise test of the
two sides' counts misses when the copies differ; a step that takes up one B and gives
off two now pins the check (R50 m1). The constructor benchmark's comment cites the
spec's 2 µs goal (R50 m2). Every spec fills its raw value, so the field takes no
`nothing` and no default (R50 m3).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 4: dfc5f79 + 88f896e** (Sig budgets at a few times their baseline; green on arrival,
red only on a per-Sig-compile mutant: S3e_c14_mut_sig2.log, `1043.0 <= 250` and `141823.0
<= 102400`)

- [ ] **Step 11: Apply and run**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary dfc5f79 -- test/test_compile_budget.jl \
    test/test_reaction_graph.jl | git apply
git show --format= --binary 88f896e -- test/test_compile_budget.jl | git apply
git show --format= --binary dfc5f79 -- src/types.jl | git apply
```

Background: `. .../replay/env.sh; focused T40_4_green test/test_compile_budget.jl
test/test_reaction_graph.jl`. Then `wait_log T40_4_green`, `summary T40_4_green`, `grep -i
'sig' $L/T40_4_green.log | head`. Expected: TCB 16, TRG 886, `exit 0`; the Sig-encoding
line reads about 52 method instances and 34 KB, under 250 and 100 KB (S3e_c13_green.log).

- [ ] **Step 12: Red on the per-Sig-compile mutant** (S3e §9 task 13: encoding compiles
  the generated `metabolites` once per Sig value)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py sig_per_mechanism
```

Expected: `mutated src/types.jl (sig_per_mechanism)`. Background: `. .../replay/env.sh;
focused T40_4_mut test/test_compile_budget.jl`; `wait_log T40_4_mut`; `summary T40_4_mut`.
Expected (`$P/proto/logs/S3e_c14_mut_sig2.log`): TCB 14 pass, 2 fail, both in "Sig
encoding: total compile for the bi-bi seeds": the instances over 250 (the log: 1043.0) and
the native code over 100 KB (141823.0 > 102400). Then, whatever the run printed:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py sig_per_mechanism --revert
git status --short
```

Expected: `restored src/types.jl`, then ` M` src/types.jl (dfc5f79's comment), TCB and TRG.

- [ ] **Step 13: Commit** (message edited: one paragraph for 88f896e)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/types.jl test/test_compile_budget.jl test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Tighten the Sig-encoding budgets to a few times their baseline

Encoding the 239 bi-bi seeds' Sigs compiles 52 method instances and 34 KB, so budgets of
4,000 instances and 1,000 KB could not catch the regression their comment names, a
compile per mechanism, which adds at least one instance per seed (R51 m2). They are now
250 instances and 100 KB. The comments that said encoding and decoding compile "once for
every mechanism" or "for every length" now say once, whatever the mechanism or the Sig
length (R51 m1).

The budgets' comment says which catches which per-seed compile: one small method
instance per seed adds about 239 instances and almost no code, so only the instance
budget fails it; the native-code budget fails a compile that emits more than about
0.3 KB per seed.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 5: ff5e59c** (free E holds no residual; green on arrival, red only on the
residual-blind mutant: S3e_r1_red_mutants.log, TRG:690 "No exception thrown")

- [ ] **Step 14: Apply and run**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary ff5e59c -- test/test_reaction_graph.jl | git apply
```

Background: `. .../replay/env.sh; focused T40_5_green test/test_reaction_graph.jl`. Then
`wait_log T40_5_green`, `summary T40_5_green`. Expected: TRG 887, `exit 0`.

- [ ] **Step 15: Red on the residual-blind mutant** (S3e §9 task 1's Start: a free-E check
  that reads label and `bound` only)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py residual_blind
```

Expected: `mutated src/types.jl (residual_blind)`. Background: `. .../replay/env.sh;
focused T40_5_mut test/test_reaction_graph.jl`; `wait_log T40_5_mut`; `summary T40_5_mut`.
Expected (`$P/proto/logs/S3e_r1_red_mutants.log`, TRG:690): TRG 886 pass, 1 fail, in "the
constructors reject what no graph represents, each with its message": the ping-pong whose
free form is F and whose E holds only A's residual, "No exception thrown". Then, whatever
the run printed:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py residual_blind --revert
git status --short
```

Expected: `restored src/types.jl`, then ` M test/test_reaction_graph.jl` only.

- [ ] **Step 16: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Pin that free E holds no residual

Row 25's tests covered a mechanism with no step, one whose free form is F and one
whose E forms all hold a metabolite, so a check blind to the residual passed them
all. A ping-pong enzyme whose free form is F, and whose E holds nothing only with
A's residual, uses every declared name, so only row 25's residual condition
rejects it.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Task end**

- [ ] **Step 17: Full suite.** Background: `. .../replay/env.sh; full T40_full`; `wait_log
  T40_full`; `summary T40_full; pristine T40_full`. Expected: `EnzymeRates.jl | 39871
  39871`, `exit 0`; `pristine` prints `12` only.
- [ ] **Step 18: Census through the graph adapter.** Background: `. .../replay/env.sh;
  CENSUS_ADAPTER=neutral_graph.jl census_run 40 S3b_end census/rename_phase3.toml`;
  `wait_log $P/census/logs/cmp_replay_T40.log`; `census_check 40`. Expected: `RESULT:
  PASS` and nothing else (every enumerated init, seed and child mechanism passes rows
  25-31).
- [ ] **Step 19: Goldens.** `git log --oneline 7205ba4..HEAD -- test/fixtures` lists only
  Task 32 unit 1's commit; both goldens passed in Step 17.
- [ ] **Step 20: Seed counts** (`<NN>` = 40). Expected: 264, 239, 7, 250,855.
- [ ] **Step 21: Constructor speed**

Background:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh $L/T40_measure.log 1800 --project=. --heap-size-hint=2500M $P/s3d2/scripts/measure.jl
```

Then `wait_log T40_measure` and `grep -E 'constructor|_canonicalize|views' $L/T40_measure.log`.
Expected (S3e_measure.log): `specs constructor (µs): (median = <m>, ...)` with `<m>` near
3.7 and under 10 (R17; D-S3d-1's revisit mark is 5 µs); `bi-bi init constructor` near 4.1;
`exit 0`. Report the medians.

- [ ] **Step 22: TE alone** (`<NN>` = 40). Expected: `TE_ELAPSED <s>`, `exit 0`.
- [ ] **Step 23: R52 record** (`<NN>` = 40). Report the three figures. The prototype at
  S3e, without the copies: 4.23 GB at dfc5f79 and 4.28 GB at 88f896e (S3e §1, review
  round 1), 50.8 s, 186.3 s (S3e §1).
- [ ] **Step 24: Checkpoint S3e** (tip 88f896e)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git diff --stat 88f896e HEAD -- src test docs ROADMAP.md Project.toml .claude
```

Expected: exactly these 18 files: Project.toml 2, ROADMAP.md 39,
docs/src/deriving/cha_king_altman.md 3, docs/src/deriving/mwc_allostery.md 8,
docs/src/deriving/ping_pong.md 2, docs/src/deriving/re_vs_ss.md 6,
docs/src/deriving/regulators.md 6, docs/src/developer.md 6, docs/src/getting_started.md 4,
docs/src/identify/enumeration_engine.md 5, docs/src/identify/tutorial.md 7, the plan (D2),
the spec 48 (D1), src/dsl.jl 4 (216be7f's docstring paragraph), src/types.jl 5, FX 35,
TRD 31, TRG 29. Then:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
export GIT_INDEX_FILE=$P/replay/ckpt.idx
git read-tree 88f896e
git diff 4f9a34a 7205ba4 -- docs/superpowers/specs ROADMAP.md | git apply --cached   # D1, D3
for c in 5c08ede 4a96d9a 0b5e637 68f46f0 33b94d8 cfe9211; do git show --format= --binary $c | git apply --cached; done   # D4
for c in 2b29fd2 be65542 0b2021d a84d194; do git show --format= --binary $c | git apply --cached; done   # docs ahead
git show --format= --binary 216be7f -- docs/src/deriving/mwc_allostery.md docs/src/deriving/cha_king_altman.md docs/src/getting_started.md src/dsl.jl | git apply --cached
git show --format= --binary 050e16f -- src/types.jl test/test_reaction_graph.jl | git apply --cached
git show --format= --binary 3962124 -- test/test_reaction_graph.jl | git apply --cached
git apply --cached $P/replay/T34_roadmap.patch
git apply --cached $P/replay/T36_roadmap.patch
git update-index --add --cacheinfo 100644,$(git rev-parse HEAD:docs/superpowers/plans/2026-10-16-graph-mechanism.md),docs/superpowers/plans/2026-10-16-graph-mechanism.md   # D2
git diff --cached --stat HEAD -- src test docs ROADMAP.md Project.toml .claude
rm -f $GIT_INDEX_FILE; unset GIT_INDEX_FILE
```

Expected: every `git apply --cached` silent and the `git diff --cached` prints nothing.
Departures in S3e: Task 36's ROADMAP entry (with dc4035b's scoping and P4-1's sentences)
and 216be7f's `@enzyme_mechanism` paragraph, ahead of their prototype position, beside
S3d2's; the messages of Tasks 36, 37 unit 2, 38 units 1-2 and 40 unit 4, as written there.

### Task 41: The boundary check on the graph [S3f]

**Replays:** `[f3d5e2c]` (S3f §9 task 1, R52 first: the boundary check ran on Step views at
every enumerator exit and derivation entry).

**Files:**
- Modify: `src/types.jl` (the graph's reads and sub-graphs), `src/admission.jl`,
  `src/thermo.jl`
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/test_reaction_graph.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: Task 34's graph and Task 20's boundary check.
- Produces: types.jl `_is_chemistry(m, e)`, `_flat_edges(m)`, `_is_free(f)`, `_flags(m)`,
  `_present_forms(m, mask)`, `_name_index(m)`, `_segments(m, mask, re) -> (; segment_of,
  segments, offsets)` (breadth-first, adjacency in `_flat_edges` order, so the bottomless
  message keeps its form order), `_outside_onlya(am)`, `_inactive_mask(am)`; admission.jl
  `_check_well_posed(m)` and `(am)` on the graph, with `_bottomless_segment(m, mask, re)`,
  `_assert_bottomed_segments(m, mask)`, `_assert_one_flip_segment(am)` (the Step versions
  go; `_bottomless_re_segment(m, steps)` stays for the moves' pre-filters); thermo.jl
  `_thermodynamic_constraints(m, flat)` and `_cycle_basis(B, stoich, ns)`,
  `_onlya_offenders(m, group, tags)`, `_onlya_haldane_violation(am)`, sharing
  `_onlya_cone_offenders` with the Step method the moves keep. FX: `_testhelper_edge_text(m,
  k)`, `_testhelper_edge(m, text)`.

**Gate:** full 39,947 (39,871 + TRG 61 + TT 15), pristine; census through the graph adapter
against `S3b_end`, zero differences; R52's figures recorded (the prototype measured ter-ter
init 35.1 s after this step, against 50.7 s before it).

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary f3d5e2c -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_reaction_graph.jl test/test_types.jl | git apply
```

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T41_red test/test_reaction_graph.jl
test/test_types.jl`. Then `wait_log T41_red`, `summary T41_red`.

Expected (S3f_c1_red.log: 1,522 pass and 63 errors there, 1,567 pass here with TRG's 45):
TRG 48 errors ("a sub-graph reads the forms, segments and bottom of the edges a mask
keeps" 5, "a sub-graph's segments and inactive edges are its Step views'" 43); TT 15 errors
("OnlyA Haldane validator" 14, ":OnlyA guard admits a ping-pong whose fused releases are
:OnlyA" 1); no failure; every error an `UndefVarError` of `_segments`, `_onlya_offenders`,
`_flat_edges`, `_flags`, `_is_free` or `_present_forms`.

- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary f3d5e2c -- src/admission.jl src/thermo.jl src/types.jl | git apply
```

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T41_green test/test_reaction_graph.jl
test/test_types.jl`. Then `wait_log T41_green`, `summary T41_green`. Expected: TRG 948, TT
694 (S3f_c1_trg_tt.log: 1,597 over the two, TRG 45 fewer), `exit 0`.

- [ ] **Step 5: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/admission.jl src/thermo.jl src/types.jl \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_reaction_graph.jl test/test_types.jl
git status --short
git commit -F - <<'EOF'
Check a mechanism's rate law on its reaction graph

The boundary check ran on Step views: every call built the views and the
Step-based segment and cycle helpers, once per enumerator output and per
derivation entry, which cost the swap its speed and memory (R52). It now reads
the graph through edge masks: rapid-equilibrium segments (_segments), the
inactive conformation's edges (_inactive_mask), the cycle basis of a list of
edges and the :OnlyA cone test (_onlya_offenders). The moves keep their
Step-based pre-filters, which share the cycle classification and the cone test.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

- [ ] **Step 6: Full suite.** Background: `. .../replay/env.sh; full T41_full`; `wait_log
  T41_full`; `summary T41_full; pristine T41_full`. Expected: `EnzymeRates.jl | 39947
  39947`, `exit 0`; `pristine` prints `12` only.
- [ ] **Step 7: Census through the graph adapter.** Background: `. .../replay/env.sh;
  CENSUS_ADAPTER=neutral_graph.jl census_run 41 S3b_end census/rename_phase3.toml`;
  `wait_log $P/census/logs/cmp_replay_T41.log`; `census_check 41`. Expected: `RESULT:
  PASS` and nothing else.
- [ ] **Step 8: Seed counts and TE alone** (`<NN>` = 41). Expected: 264 / 239 / 7 /
  250,855; `TE_ELAPSED <s>`, `exit 0`.
- [ ] **Step 9: R52 record** (`<NN>` = 41). Report the three figures beside Task 40's. The
  prototype, without the copies: 4.31 GB, 35.1 s, 154.0 s after this step, against 4.18
  GB, 50.7 s, 173.0 s before it (S3f §1).

### Task 42: Constants, the joint solve and the derivation over edges [S3f]

**Replays:** `[19e3fe8 + f6be872]` plus 216be7f's ROADMAP entry below, one commit (S3f §9
task 2 folds task 3, so the four-line `_flat_names` adapter 19e3fe8 adds and f6be872
deletes never lands).

**Files:**
- Modify: `src/types.jl` (`Krapid`, `Kfor`, `Krev` hold `edge::Int`; edge reads and
  ranks), `src/names.jl`, `src/thermo.jl`, `src/rate_eq_derivation.jl`, `ROADMAP.md`
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/test_accessors.jl`, `test/test_fitting_scale.jl`, `test/test_joint_solve.jl`,
  `test/test_mechanism_enumeration.jl`, `test/test_oracles.jl`,
  `test/test_rate_eq_derivation.jl`, `test/test_reaction_graph.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: Task 41's sub-graphs and `_flat_edges`.
- Produces: types.jl `_is_iso(e)`, `_binds_product(m, e)`, `_base_rank(m)`,
  `_edge_base_rank(m)`; thermo.jl `_Column` with `edge`, `_priority_key(m, k, slot, block,
  edge_rank)`, `_ordered_columns(m)` shared by `_joint_system(m, cols)` and
  `_kept_members(m, cols)` (R20), `_fundamental_cycles(m, flat)`, `_column_parameter(m, c)`,
  `_dependent_param_exprs(m)`, `_count_kind(e::Edge)`, the `:OnlyA` message over edges;
  names.jl `_constant_sides(m, k, slot, state, kept)`, the step-constant `name` methods by
  edge, `_edge_constants(m, k, state)` (`_class_sides` and `_step_constants` go);
  rate_eq_derivation.jl `_state_mask`, `_edge_names(m, state, kept)` (one entry per edge,
  empty outside the state's law), `_segment_root(m, segment)`, `_re_weight_ratio`,
  `_compute_alpha(m, flat, conc_set, segments, names)`, `_raw_symbolic_rate_polys(m, mask,
  names)`, `_eq_complexity(m, mask)`, `_eq_complexity(m)`, `_assert_derivable(m, mask)`,
  `_state_rate_polys(am, state, kept = _kept_members(am))` (`_step_names`, `_state_groups`
  and the Step derivation methods go).

**Gate:** full 39,958 (39,947 + 8 + 3), pristine; census through the graph adapter against
`S3b_end`, zero differences (542 laws exact); TRD's Expr-shape and flat-string pins and
`rate_equation`'s 0 allocations under 120 ns per spec hold (Step 7).

**Departure (written out):** ROADMAP.md: the entry "**`eq_complexity_filter` counts only the
active conformation of an allosteric mechanism.**" takes its text in `git show
7c2a657:ROADMAP.md` (lines 24-29): `_eq_complexity` measures the catalytic graph, the
active conformation, not `steps(am)`. Step 4 applies it as a patch. The commit message is
19e3fe8's subject and body, f6be872's body, and one sentence for the entry.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 19e3fe8 -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_accessors.jl \
    test/test_fitting_scale.jl test/test_joint_solve.jl \
    test/test_mechanism_enumeration.jl test/test_oracles.jl \
    test/test_rate_eq_derivation.jl test/test_reaction_graph.jl \
    test/test_types.jl | git apply
git show --format= --binary f6be872 -- test/test_mechanism_enumeration.jl \
    test/test_rate_eq_derivation.jl | git apply
```

- [ ] **Step 2: Red** (two runs: the prototype took 19e3fe8's red over four files and
  f6be872's on one testset)

Background, one after the other:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
focused T42_red test/test_joint_solve.jl test/test_types.jl test/test_reaction_graph.jl \
    test/test_accessors.jl
$P/jl.sh $L/T42_red_root.log 1800 --project=. --heap-size-hint=2500M \
    $P/s3c/scripts/run_named.jl test/test_rate_eq_derivation.jl \
    "a segment's root ties go to the canonical form order"
```

Then `wait_log T42_red`, `summary T42_red`, `wait_log T42_red_root`.

Expected: `T42_red`: 1,685 pass, 3 fail, 257 errors over the four files (S3f_c2_red.log:
1,640 pass there, TRG 45 more here), among them `UndefVarError` of `_edge_constants`
(131), `_edge_base_rank` (9), `_base_rank` (2), `_is_iso` and `_binds_product`, and FX's
`_testhelper_edge_text` given a kept Step. `T42_red_root`: `test/test_rate_eq_derivation.jl
| a segment's root ties go to the canonical form order | FAIL: Some tests did not pass: 1
passed, 0 failed, 4 errored, 0 broken.` (S3f_c3_red.log: no
`_segment_root(::EnzymeRates.Mechanism, ...)`).

- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 19e3fe8 -- src/names.jl src/rate_eq_derivation.jl \
    src/thermo.jl src/types.jl | git apply
git show --format= --binary f6be872 -- src/rate_eq_derivation.jl | git apply
! grep -n "_flat_names" src/*.jl
```

Expected: no output: the `_flat_names` adapter never lands.

- [ ] **Step 4: Apply the ROADMAP departure**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
cat > $P/replay/T42_roadmap.patch <<'PATCH'
--- a/ROADMAP.md
+++ b/ROADMAP.md
@@ -23,8 +23,8 @@
   in parallel only across `pmap` workers, each of which builds its own problem.
 - **`eq_complexity_filter` counts only the active conformation of an allosteric
   mechanism.** `_eq_complexity` (`src/rate_eq_derivation.jl`) measures an allosteric
-  mechanism's V×τ on its active-state graph (`steps(am)`) alone, while the generated MWC
-  rate equation evaluates the polynomials of both conformations, so the filter
+  mechanism's V×τ on its catalytic graph, the active conformation, alone, while the
+  generated MWC rate equation evaluates the polynomials of both conformations, so the filter
   undercounts what an allosteric equation costs and can admit one that evaluates more
   terms than its limit.
 - **`BigInt` arithmetic inflates the GC's live-byte count.**
PATCH
git apply $P/replay/T42_roadmap.patch
diff <(sed -n '/^- \*\*`eq_complexity_filter`/,/terms than its limit\./p' ROADMAP.md) \
     <(git show 7c2a657:ROADMAP.md | sed -n 24,29p)
```

Expected: no output.

- [ ] **Step 5: Green**

Background: `. .../replay/env.sh; focused T42_green test/test_joint_solve.jl
test/test_types.jl test/test_reaction_graph.jl test/test_accessors.jl
test/test_rate_eq_derivation.jl test/test_oracles.jl test/test_fitting_scale.jl
test/test_mechanism_enumeration.jl`. Then `wait_log T42_green` until it ends (about 7
min), `summary T42_green`.

Expected: TJS 455, TT 694, TRG 954, TA 155, TRD 3,189, TO 667, TFS 104, TE 32,545,
`exit 0`.

- [ ] **Step 6: Commit** (message edited: f6be872's body and the ROADMAP sentence)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/names.jl src/rate_eq_derivation.jl src/thermo.jl src/types.jl ROADMAP.md \
    test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_accessors.jl \
    test/test_fitting_scale.jl test/test_joint_solve.jl \
    test/test_mechanism_enumeration.jl test/test_oracles.jl \
    test/test_rate_eq_derivation.jl test/test_reaction_graph.jl test/test_types.jl
git status --short
git commit -F - <<'EOF'
Name constants and solve the joint system over edges

Parameters name an edge of the catalytic graph instead of a Step, and the
joint solve's columns, cycles, priority keys and kept members read edges, so
naming a mechanism's constants no longer builds its Step views. The kept
members come from the solve's ordered columns (_ordered_columns), shared by
the solve and the names (R20); the content tie of the priority key is the
edges' base rank, the order of today's content key. The :OnlyA message and the
split's counter read the graph too.

The King–Altman derivation read a conformation's Step views and their
segments; it now reads the edges a mask keeps (_state_mask): the segments and
alphas over the edges group by group (_flat_edges), each segment rooted at its
form with the fewest bound metabolites, then without a residual, then the
smaller index. The size guard reads the same segments. The order of every
operation is the Step derivation's, so the polynomials are the same.

ROADMAP: the eq_complexity_filter entry says _eq_complexity measures the
catalytic graph, the active conformation, not steps(am).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

- [ ] **Step 7: Full suite.** Background: `. .../replay/env.sh; full T42_full`; `wait_log
  T42_full`; `summary T42_full; pristine T42_full`. Expected: `EnzymeRates.jl | 39958
  39958`, `exit 0`; `pristine` prints `12` only. The gate's `rate_equation` checks are TRD
  testsets this count includes: every spec's "Performance" (0 allocations, under 120 ns),
  "rate_equation time per call grows at most in proportion to the law's size", "the
  generated rate_equation body is in Horner form" and "rate_equation_string prints flat
  +/* sums (no nested parens)"; Step 5's TRD run passed them too.
- [ ] **Step 8: Census through the graph adapter.** Background: `. .../replay/env.sh;
  CENSUS_ADAPTER=neutral_graph.jl census_run 42 S3b_end census/rename_phase3.toml`;
  `wait_log $P/census/logs/cmp_replay_T42.log`; `census_check 42`. Expected: `RESULT:
  PASS` and nothing else.

### Task 43: Decode once per generator; reject labels A and I [S3f]

**Replays:** `[43519c4]`, `[a75ef55 + 7fe5c6d]` (S3f §9 tasks 4-5; 7fe5c6d fixes a TJS
comment the swap made false, folded into the last S3f code commit).

**Files:**
- Modify: `src/rate_eq_derivation.jl`, `src/types.jl` (`_check_conformation_label`),
  `src/thermo.jl` (the `_dependent_param_exprs` docstring), `src/dsl.jl` (the
  `@enzyme_mechanism` docstring)
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/test_joint_solve.jl` (a comment), `test/test_rate_eq_derivation.jl`,
  `test/test_reaction_graph.jl`

**Interfaces:**
- Consumes: Task 42's derivation over edges; `_derivation_input` (Task 20).
- Produces: `_num_den_exprs(m::Mechanism)`, `_num_den_exprs(am::AllostericMechanism)`,
  `_assert_no_name_clash(m, ...)`, `_rate_body(m, ...)`, `_equation_text(m, ...)` take the
  decoded mechanism; the `rate_equation` generator and `rate_equation_string` decode once
  through `_derivation_input`; `_mwc_state_polys(am, mets, kept)` computes `kept` once for
  both conformations; FX's oracle law value reads `_num_den_exprs` of the decoded
  mechanism. `_check_conformation_label` rejects `A` and `I` first (R43): "`A` cannot be a
  conformation label: constant names tag a conformation with A_ or I_ (`k_A_ES_to_EP`), so
  the labels A and I could give two constants one name; choose another label." The
  name-clash guard in `_dependent_param_exprs` stays, tested on a graph built without the
  label check (`_testhelper_relabelled_raw`, `_testhelper_unchecked`, TRD).

**Gate:** full 39,964 (39,958 + 6), pristine; census through the graph adapter against
`S3b_end`, zero differences; seed counts; R52's figures; the R18 probe; checkpoint S3f.

**Unit 1: 43519c4** (decode once)

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 43519c4 -- \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_rate_eq_derivation.jl | git apply
```

- [ ] **Step 2: Red**

Background: `. .../replay/env.sh; focused T43_1_red test/test_rate_eq_derivation.jl`. Then
`wait_log T43_1_red`, `summary T43_1_red`.

Expected (S3f_c4_red.log, plus Task 1's 2 passes): TRD 2,841 pass, 223 errors, no
failure: `MethodError`s of `_num_den_exprs(::EnzymeRates.Mechanism)` (each spec's
"Mass-action oracle", through FX's law value), `_num_den_exprs(::AllostericMechanism)`
and `_rate_body(::EnzymeRates.Mechanism, ...)` (the generated-body readers).

- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 43519c4 -- src/rate_eq_derivation.jl | git apply
```

- [ ] **Step 4: Green**

Background: `. .../replay/env.sh; focused T43_1_green test/test_rate_eq_derivation.jl
test/test_oracles.jl`. Then `wait_log T43_1_green`, `summary T43_1_green`. Expected: TRD
3,189, TO 667 (FX's law value), `exit 0`.

- [ ] **Step 5: Commit**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/rate_eq_derivation.jl \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_rate_eq_derivation.jl
git status --short
git commit -F - <<'EOF'
Decode a compiled mechanism once per generator and derive from the value

The rate_equation generator and rate_equation_string decoded the Sig once for
the solve, again for the polynomials and again for the metabolite names; the
MWC assembly decoded it twice more. They now decode and check it once
(_derivation_input) and pass the mechanism to the solve, the generated body,
the printed law and the MWC assembly, which computes the kept members once for
both conformations.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Unit 2: a75ef55 + 7fe5c6d** (labels A and I; the tied lead's comment)

- [ ] **Step 6: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary a75ef55 -- test/test_rate_eq_derivation.jl \
    test/test_reaction_graph.jl | git apply
git show --format= --binary 7fe5c6d -- test/test_joint_solve.jl | git apply
```

- [ ] **Step 7: Red**

Background: `. .../replay/env.sh; focused T43_2_red test/test_rate_eq_derivation.jl
test/test_reaction_graph.jl`. Then `wait_log T43_2_red`, `summary T43_2_red`.

Expected (S3f_c5_red.log, 4 fail): TRD 3,191 pass, 2 fail ("constants that would share a
name are rejected": the macro accepts labels A and I, "No exception thrown"); TRG 954
pass, 2 fail (rows 30-31: the constructor accepts labels A and I).

- [ ] **Step 8: Apply the source** (7fe5c6d has none)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary a75ef55 -- src/dsl.jl src/thermo.jl src/types.jl | git apply
```

- [ ] **Step 9: Green**

Background: `. .../replay/env.sh; focused T43_2_green test/test_rate_eq_derivation.jl
test/test_reaction_graph.jl test/test_dsl.jl test/test_types.jl test/test_joint_solve.jl`.
Then `wait_log T43_2_green`, `summary T43_2_green`. Expected: TRD 3,193, TRG 956, TD 257,
TT 694, TJS 455 (S3f_c5_green.log: 5,053 over the first four, TRD 2 and TRG 45 fewer),
`exit 0`.

- [ ] **Step 10: Commit** (message edited: one sentence for 7fe5c6d)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git add src/dsl.jl src/thermo.jl src/types.jl test/test_joint_solve.jl \
    test/test_rate_eq_derivation.jl test/test_reaction_graph.jl
git status --short
git commit -F - <<'EOF'
Reject the conformation labels A and I

Constant names tag a conformation with A_ or I_, and on the reaction graph an
:EqualAI step at a form labelled A still renders the tagged name of another
step (A(Eb) <--> A + Eb and the active constant of E(b) <--> A(Eb) are both
k_A_Eb_to_AEb), so the derivation-time guard stayed reachable (R43). The label
check, which the macros and the constructors share, now rejects A and I. The
guard stays and is tested on a graph built without the label check.

A test comment names the tied group's lead by the function that picks it,
_align_tied_groups!, in place of _step_index_key, which the swap deleted.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git status --short; git log --oneline -1
```

**Task end**

- [ ] **Step 11: Full suite.** Background: `. .../replay/env.sh; full T43_full`; `wait_log
  T43_full`; `summary T43_full; pristine T43_full`. Expected: `EnzymeRates.jl | 39964
  39964`, `exit 0`; `pristine` prints `12` only.
- [ ] **Step 12: Census through the graph adapter.** Background: `. .../replay/env.sh;
  CENSUS_ADAPTER=neutral_graph.jl census_run 43 S3b_end census/rename_phase3.toml`;
  `wait_log $P/census/logs/cmp_replay_T43.log`; `census_check 43`. Expected: `RESULT:
  PASS` and nothing else (S3f §5: the S3f corpus equals S3b_end).
- [ ] **Step 13: Seed counts and TE alone** (`<NN>` = 43). Expected: 264 / 239 / 7 /
  250,855; `TE_ELAPSED <s>`, `exit 0`.
- [ ] **Step 14: R52 record** (`<NN>` = 43). Report the three figures beside Tasks 40 and
  41's. The prototype at S3f's end, without the copies: 4.19 GB, 35.7 s, 151.5 s (S3f §1).
- [ ] **Step 15: R18 probe**

Background:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
$P/jl.sh $L/T43_r18.log 600 --project=. $P/s3f/probe/r18.jl
```

Then `wait_log T43_r18` and `grep -E '^NOT|^NO METHOD' $L/T43_r18.log`. Expected (S3f_r18.log,
64 methods S3f adds or ports): exactly five lines, the shapes they had before S3f:
`_column_parameter(Mechanism, _Column) -> Parameter`, `_dependent_param_exprs` of
`Mechanism` and of `AllostericMechanism` -> `Tuple{Dict{Symbol, Union{Float64, Expr,
Symbol}}, Tuple{Vararg{Symbol}}}`, `_num_den_exprs(Mechanism) -> Tuple{Any, Any}`,
`_num_den_exprs(AllostericMechanism) -> Tuple{Any, Expr}`; no `NO METHOD` line; `exit 0`.

- [ ] **Step 16: Checkpoint S3f** (tip 7fe5c6d)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git diff --stat 7fe5c6d HEAD -- src test docs ROADMAP.md Project.toml .claude
```

Expected: exactly these 18 files: Project.toml 2, ROADMAP.md 43,
docs/src/deriving/cha_king_altman.md 3, docs/src/deriving/mwc_allostery.md 8,
docs/src/deriving/ping_pong.md 2, docs/src/deriving/re_vs_ss.md 6,
docs/src/deriving/regulators.md 6, docs/src/developer.md 6, docs/src/getting_started.md 4,
docs/src/identify/enumeration_engine.md 5, docs/src/identify/tutorial.md 7, the plan (D2),
the spec 48 (D1), src/dsl.jl 4, src/types.jl 5, FX 35, TRD 31, TRG 29. Then:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
export GIT_INDEX_FILE=$P/replay/ckpt.idx
git read-tree 7fe5c6d
git diff 4f9a34a 7205ba4 -- docs/superpowers/specs ROADMAP.md | git apply --cached   # D1, D3
for c in 5c08ede 4a96d9a 0b5e637 68f46f0 33b94d8 cfe9211; do git show --format= --binary $c | git apply --cached; done   # D4
for c in 2b29fd2 be65542 0b2021d a84d194; do git show --format= --binary $c | git apply --cached; done   # docs ahead
git show --format= --binary 216be7f -- docs/src/deriving/mwc_allostery.md docs/src/deriving/cha_king_altman.md docs/src/getting_started.md src/dsl.jl | git apply --cached
git show --format= --binary 050e16f -- src/types.jl test/test_reaction_graph.jl | git apply --cached
git show --format= --binary 3962124 -- test/test_reaction_graph.jl | git apply --cached
git apply --cached $P/replay/T34_roadmap.patch
git apply --cached $P/replay/T36_roadmap.patch
git apply --cached $P/replay/T42_roadmap.patch
git update-index --add --cacheinfo 100644,$(git rev-parse HEAD:docs/superpowers/plans/2026-10-16-graph-mechanism.md),docs/superpowers/plans/2026-10-16-graph-mechanism.md   # D2
git diff --cached --stat HEAD -- src test docs ROADMAP.md Project.toml .claude
rm -f $GIT_INDEX_FILE; unset GIT_INDEX_FILE
```

Expected: every `git apply --cached` silent and the `git diff --cached` prints nothing.
Departures in S3f: Task 42's ROADMAP entry, ahead of its prototype position, beside S3e's;
the messages of Task 42 and Task 43 unit 2, as written there.

## Tasks 44-63: S3g moves, S3h seeds, S3i views removed, S4 show, parser and CSV, CLAUDE.md, final

### Task 44: Segment sides and cycle rows; A/I pins [S3g]

**Replays:** `[d427228]`, `[ee46733]` (S3g §9 tasks 1-2: R54 m4, R53 m1 and R54 m6; then
R54 m3, R54 m2 and R53 m3).

**Files:**
- Create (scratch, not committed): `$P/replay/roadmap_entry.py`, `$P/replay/r55m2.py`,
  `$P/replay/ckpt_p4.sh`, `$P/replay/p41_roadmap.patch`, `$P/replay/roadmap_final.md`
- Modify: `src/admission.jl` (`_bottomless_segment`), `src/thermo.jl`
  (`_thermodynamic_constraints(m, flat)`), `src/types.jl` (`_mechanism_from_sig`'s
  docstring)
- Test: `test/test_joint_solve.jl` (the tied-lead comment),
  `test/test_rate_eq_derivation.jl`

**Interfaces:**
- Consumes: `_segments`, `_name_index`, `_present_forms` (Tasks 41-42); TRD's
  `_testhelper_relabelled_block`, `_testhelper_relabelled_raw`, `_testhelper_unchecked`,
  `_testhelper_written_params`.
- Produces: `_bottomless_segment(m, mask, re)` reads each name's product and substrate side
  through `_name_index`; `_thermodynamic_constraints(m, flat)` takes its rows from
  `_present_forms`; `_testhelper_unchecked(raw)` keeps `raw.sites`; the textbook bridge
  skips an empty name list only in the inactive state; labels `A` and `I` pinned at
  `@enzyme_mechanism` and `AllostericMechanism(raw)`. The five helpers of Step 1, which the
  ROADMAP departures (Tasks 48, 54, 60) and every checkpoint of Tasks 44-63 call.

**Gate:** F 44 = 39,968; C 44 (adapter): zero differences.

- [ ] **Step 1: Helpers for Tasks 44-63** (once; scratch files outside the repository)

`roadmap_entry.py` sets one ROADMAP entry to its text in another file, or deletes it;
`r55m2.py` makes R55 m2's indentation; `ckpt_p4.sh` is the CK runner;
`roadmap_final.md` is the ROADMAP the plan ends with: 7c2a657's with 7205ba4's three P4-1
hunks (`p41_roadmap.patch`). Every ROADMAP departure of Tasks 44-63 takes its entry's text
from it.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
cat > $P/replay/roadmap_entry.py <<'EOF'
# ABOUTME: Sets one ROADMAP.md entry to the same entry's text in another file, or deletes it.
# ABOUTME: Usage: roadmap_entry.py <title start> <source file | delete> <ROADMAP.md>.
import sys


def block(lines, title):
    """First line, end of text and end of block of the one entry whose first line starts
    with "- **" + title; the block runs to the next entry or section, blank lines included."""
    starts = [i for i, l in enumerate(lines) if l.startswith("- **" + title)]
    if len(starts) != 1:
        sys.exit(f"{len(starts)} entries start with {title!r}")
    s = starts[0]
    e = s + 1
    while e < len(lines) and not lines[e].startswith(("- **", "## ")):
        e += 1
    t = e
    while lines[t - 1] == "":
        t -= 1
    return s, t, e


title, source, path = sys.argv[1:4]
text = open(path, encoding="utf-8").read()
assert text.endswith("\n"), path
lines = text[:-1].split("\n")
s, t, e = block(lines, title)
if source == "delete":
    lines[s:e] = []
else:
    new = open(source, encoding="utf-8").read().rstrip("\n").split("\n")
    fs, ft, _ = block(new, title)
    lines[s:t] = new[fs:ft]
open(path, "w", encoding="utf-8").write("\n".join(lines) + "\n")
EOF
cat > $P/replay/r55m2.py <<'EOF'
# ABOUTME: Indents the three expected children of TE's "uni-uni all-:EqualAI: 3 cat
# ABOUTME: relaxations" testset by four spaces, not their closing bracket (R55 m2).
import sys

path = sys.argv[1]
lines = open(path, encoding="utf-8").read().split("\n")
(t,) = [i for i, l in enumerate(lines) if "uni-uni all-:EqualAI: 3 cat relaxations" in l]
s = next(i for i in range(t, len(lines)) if lines[i] == "    expected = [")
e = next(i for i in range(s + 1, len(lines)) if lines[i] == "    ]")
assert all(lines[i] for i in range(s + 1, e))
lines[s + 1:e] = ["    " + l for l in lines[s + 1:e]]
open(path, "w", encoding="utf-8").write("\n".join(lines))
print(f"indented lines {s + 2}-{e}")
EOF
cat > $P/replay/ckpt_p4.sh <<'EOF'
#!/bin/bash
# ABOUTME: Builds a stage's expected tree in a scratch index (the prototype's stage tip, the
# ABOUTME: standing departures D1-D4 and the given ones) and prints how HEAD differs from it.
# Usage: ckpt_p4.sh <tip> <op>...   op: c:<sha> | p:<sha>:<path> | e:<title>=<file|delete> | i
set -e
P=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3
cd ${CK_REPO:-/home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/graph_mechanism_wt}
D=${CK_DIR:-$P/replay}
H=${CK_HEAD:-HEAD}
export GIT_INDEX_FILE=$D/ckpt.idx
T=$D/ckpt.tmp
rm -f $GIT_INDEX_FILE
git read-tree "$1"
shift
git diff 4f9a34a 7205ba4 -- docs/superpowers/specs ROADMAP.md | git apply --cached   # D1, D3
for c in 5c08ede 4a96d9a 0b5e637 68f46f0 33b94d8 cfe9211; do                      # D4
    git show --format= --binary $c | git apply --cached
done
plan=docs/superpowers/plans/2026-10-16-graph-mechanism.md                             # D2
if git cat-file -e $H:$plan 2>/dev/null; then
    git update-index --add --cacheinfo 100644,$(git rev-parse $H:$plan),$plan
fi
edit() {  # <path> <command...>: runs the command on a copy of <path>, which it takes last
    local path=$1
    shift
    git show :$path > $T
    "$@" $T > /dev/null
    git update-index --cacheinfo 100644,$(git hash-object -w $T),$path
}
for op in "$@"; do
    case $op in
        c:*) git show --format= --binary ${op#c:} | git apply --cached ;;
        p:*) s=${op#p:}; git show --format= --binary ${s%%:*} -- ${s#*:} | git apply --cached ;;
        e:*) s=${op#e:}; edit ROADMAP.md python3 $D/roadmap_entry.py "${s%%=*}" "${s#*=}" ;;
        i) edit test/test_mechanism_enumeration.jl python3 $D/r55m2.py ;;
        *) echo "unknown op: $op"; exit 1 ;;
    esac
done
git diff --cached --stat $H -- src test docs ROADMAP.md Project.toml .claude
rm -f $GIT_INDEX_FILE $T
EOF
chmod +x $P/replay/ckpt_p4.sh
cat > $P/replay/p41_roadmap.patch <<'EOF'
--- ROADMAP.md
+++ ROADMAP.md
@@ -107,8 +107,9 @@
   `:OnlyA` group already keeps its own class. A tag set is admitted when an exact
   feasibility test passes and every form the free enzyme cannot reach in the inactive
   conformation lies behind a binding, on its bound side; the limit then equals deletion
-  and pruning, and enumeration stops producing the `:OnlyA` ping-pong mechanisms this
-  rejects. The conformations will flip at free E, which every mechanism holds; until
+  and pruning. A ping-pong with both half-reactions `:OnlyA` fails this, so the tagging
+  move will stop one half per V-type child instead (spec section 7, and the entry
+  below). The conformations will flip at free E, which every mechanism holds; until
   then the boundary check admits a second free form of an allosteric mechanism only in
   rapid equilibrium with the first (`_assert_one_flip_segment`, `src/admission.jl`),
   where the choice only rescales `L`. The rewrite's admission PR brings the limits,
@@ -117,6 +118,17 @@
   651a821, plan 9c95f12, prototype 6f0011e), whose prototype matched the mass-action
   ground truth on all six uni-uni tag combinations and would delete the `d_free`
   cross-weighting; `test/allosteric_ground_truth.jl` is the oracle. Priority: high.
+- **Let an allosteric ping-pong stop both half-reactions in the inactive conformation.**
+  PDH and OGDH are ping-pong and allosteric, so allosteric ping-pong must stay in the
+  search. With the conformations flipping only at free E, a ping-pong whose inactive
+  conformation runs neither half-reaction strands the covalent forms between the halves:
+  their weight, and with it the rate, depends on how fast each half stops (by up to 100%
+  at one catalytic subunit), so admission rejects it. The admission PR keeps allosteric
+  ping-pong reachable by stopping one half per V-type child, which covers PDH, whose
+  phosphorylation stops E1 and leaves E2 and E3 running. Stopping both halves needs the
+  conformations to flip at the covalent form as well, with an equilibrium constant of
+  its own, and the multi-subunit rate formula rederived for two flip forms. Take it up
+  with iso mechanisms, which add free forms of their own. Priority: medium.
 - **Flag one-subunit K-type variants whose L is a phantom instead of fitting them.**
   At multiplicity 1, `_expand_to_allosteric` (`src/mechanism_enumeration.jl`) emits
   K-type variants whose L the data cannot determine, so the fitter optimizes a flat
@@ -159,7 +171,9 @@
   label, and `_requires_hyperbolic_catalysis` (`src/mechanism_enumeration.jl`) is the
   extension hook. Two-conformation kinetic cooperativity stays parked: on PGK data no
   conformational model beat MWC with steady-state steps at equal parameter count, so
-  reopen it only for a dataset whose curvature MWC cannot fit. Priority: low.
+  reopen it only for a dataset whose curvature MWC cannot fit. Revisit allosteric
+  ping-pong with this work: a flip at a second free form is the same question (see
+  "Let an allosteric ping-pong stop both half-reactions"). Priority: low.
 
 ## Ideas
 
EOF
git show 7c2a657:ROADMAP.md > $P/replay/roadmap_final.md
patch -s --no-backup-if-mismatch $P/replay/roadmap_final.md < $P/replay/p41_roadmap.patch
(cd $P/replay && md5sum roadmap_entry.py r55m2.py ckpt_p4.sh p41_roadmap.patch roadmap_final.md)
```

Expected: `patch` prints nothing, and

```
1217072d153848f1a971b408a0f74bfc  roadmap_entry.py
628903646a92e95120086a0384605507  r55m2.py
a57779cbdc9fcc8df09c1241dc16cc77  ckpt_p4.sh
05a47bfe7956fe9a8e9519cef5d6e24b  p41_roadmap.patch
720318593664be4e079dfbd63825d011  roadmap_final.md
```

A different sum means a heredoc was mistyped: write the file again.

**Unit 1: d427228** (R54 m4, R53 m1, R54 m6: a refactor; green before, at Task 43, and
after)

- [ ] **Step 2: Apply the test path** (a comment)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary d427228 -- test/test_joint_solve.jl | git apply
```

- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary d427228 -- src/admission.jl src/thermo.jl src/types.jl | git apply
```

- [ ] **Step 4: Green** (background). Gt `$L/T44_1_green.log`
  `test_joint_solve.jl test_types.jl test_reaction_graph.jl test_rate_eq_derivation.jl`:
  455, 694, 956, 3,193.
- [ ] **Step 5: Commit.** K `d427228` `src/admission.jl src/thermo.jl src/types.jl
  test/test_joint_solve.jl` ("Read segment sides through the name index and cycle rows
  through the touched forms").

**Unit 2: ee46733** (R54 m3, R54 m2, R53 m3: test-only, pins of existing behavior, green on
arrival; the label check rejects `A` and `I` since Task 43)

- [ ] **Step 6: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary ee46733 -- test/test_rate_eq_derivation.jl | git apply
```

- [ ] **Step 7: Green** (background). Gt `$L/T44_2_green.log` `test_rate_eq_derivation.jl`:
  3,197 (+4: each label now pinned at three entries, where one was).
- [ ] **Step 8: Commit.** K `ee46733` `test/test_rate_eq_derivation.jl` ("Pin labels A and I
  at the plain macro and the allosteric constructor").

**Task end**

- [ ] **Step 9: Full suite** (background). F 44: 39,968.
- [ ] **Step 10: Census** (background). C 44: zero differences.

### Task 45: Flip on the graph [S3g]

**Replays:** `[8752036]` (S3g §9 task 3).

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_expand_re_to_ss`, `_minimal_flips(admissible, m,
  eligible)`, `_chain_flank_groups(m)`, `_flux_carrying_edges(m, re)`,
  `_flux_carrying_groups(m, re)`, `_hyperbolic_catalysis`, `_re_segment_count(m, re)`;
  the Step `_re_segment_count` goes), `src/types.jl` (`_build`, `_modifier_mask`,
  `_raw_form`, the `_raw` overrides, `_raw_sites`)
- Test: `test/test_mechanism_enumeration.jl`

**Interfaces:**
- Consumes: `_segments`, `_flags`, `_bottomless_segment` (Tasks 41, 44); `_raw(m)`
  (Task 31); the checked constructors `Mechanism(raw)`, `AllostericMechanism(raw)`.
- Produces: `_build(raw)` (`Mechanism(raw)` or `AllostericMechanism(raw)` by the raw
  value's type); `_raw(m, re = _flags(m), group = m.group)`, `_raw(am, re, group, tags)`,
  `_raw_form(m, f, extra...)`, `_raw_sites(am)`; `_modifier_mask(m)`; the flip reads the
  graph, a candidate a flag per edge, each child `_build(_raw(m, re))`;
  `_flux_carrying_edges(m, re)`, `_flux_carrying_groups(m, re)`, `_re_segment_count(m,
  re)`. The seed variants keep the Step `_minimal_flips`, `_flux_carrying_steps`,
  `_flux_carrying_groups` and `_all_steady_state` (D-S3g-1) until Task 51. TE:
  `_testhelper_segment_count(m)`.

**Gate:** F 45 = 39,968; C 45 (adapter): zero differences.

- [ ] **Step 1: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 8752036 -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 2: Red** (background). Rt `$L/T45_red.log` `test_mechanism_enumeration.jl`.
  Expected (`$P/proto/logs/S3g_red_c3.log`): 32,494 pass and 15 errors, every one a
  `MethodError` for the graph signatures (`_re_segment_count(m, re)` 12,
  `_flux_carrying_groups(m, re)` 2, `_raw(m, re)` 1): "_re_segment_count" 5,
  "_expand_re_to_ss: a zero-flux pair is absent, its supersets reachable" 3,
  "Mechanism — flip past a bottomless RE segment" 2, "_expand_re_to_ss: a split metabolite
  flips as a pair, never singly" 2, "_flux_carrying_groups" 2, "_expand_re_to_ss: split
  ter-ter pairs the A and B parts" 1.
- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 8752036 -- src/mechanism_enumeration.jl src/types.jl | git apply
```

- [ ] **Step 4: Green** (background). Gt `$L/T45_green.log` `test_mechanism_enumeration.jl`:
  32,545.
- [ ] **Step 5: Commit.** K `8752036` `src/mechanism_enumeration.jl src/types.jl
  test/test_mechanism_enumeration.jl` ("Flip rapid-equilibrium groups to steady state on
  the reaction graph").
- [ ] **Step 6: Full suite** (background). F 45: 39,968.
- [ ] **Step 7: Census** (background). C 45: zero differences.

### Task 46: Copy rule and dead-end move on the graph [S3g]

**Replays:** `[0f795c4]` (S3g §9 task 4).

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_GaugeClass`, `_NO_RESCALING`,
  `_one_state_class`, `_productive_twin(m, mask)`, `_all_twin`, `_gauge_rescaling`,
  `_redundant_copy_groups`, `_forms_where_free(m, entry)`, `_dead_end_raw`, both
  `_dead_end_child`, `_expand_add_dead_end_regulator`), `src/rate_eq_derivation.jl`
  (`_inactive_groups` and `_reachable_from_free` go), `src/types.jl` (`_entry(m, x::Met)`)
- Test: `test/test_mechanism_enumeration.jl`, `test/test_reaction_graph.jl`

**Interfaces:**
- Consumes: Task 45's `_raw` overrides, `_build`, `_modifier_mask`; `_inactive_mask`
  (Task 42); `_segments`; `$P/replay/mutate.py` (Task 0 Step 8): `k3`, `one_index`.
- Produces: `_GaugeClass` (a group, a form or a state: group 2 and form 2 are two classes);
  `_productive_twin(m, mask)` returns `(site, ligand::Symbol) -> form` (0 for none,
  D-S3g-2); `_redundant_copy_groups` reads the inactive state from `_inactive_mask`;
  `_forms_where_free(m, entry) -> BitVector`; `_dead_end_raw(m, reg, sites)`;
  `_dead_end_child(m::Mechanism, raw, reg)` (canonicalizes, finds twins on the candidate,
  builds checked only when kept, K3) and the allosteric method (D-S3g-4);
  `_entry(m, x::Met)` (0 when absent). TE: `_testhelper_free_at(m, name, modifier)` (in
  place of `_testhelper_form_labels`); TRG: `_testhelper_inactive_steps(am)`, the Step-view
  reach that Task 53 rewrites.

**Gate:** F 46 = 39,980; C 46 (adapter): zero differences; CD 46.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 0f795c4 -- test/test_mechanism_enumeration.jl test/test_reaction_graph.jl | git apply
```

- [ ] **Step 2: Red** (background). Rt `$L/T46_red.log`
  `test_mechanism_enumeration.jl test_reaction_graph.jl`. Expected
  (`$P/proto/logs/S3g_red_c4.log`): 20 errors, all in TE (`_entry` undefined 15,
  `_productive_twin(m, mask)` 5): "_forms_where_free" 6, "_redundant_copy_groups: every
  complex a twin and a consistent dwell gauge" 4, "a Theorell–Chance step frees B at E(A)
  and P at E(Q)" 4, "a copy binds at its own sites, apart from its reactant's" 4,
  "_gauge_rescaling: group 2 and form 2 are two classes" 1, "_productive_twin" 1; TRG
  passes, 956. The two new testsets also fail on the K3 and one-index mutants of the source
  (`$P/proto/logs/S3g_mut_k3.log`, `S3g_mut_class.log`; S3g §1), which Step 5 runs.
- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 0f795c4 -- src/mechanism_enumeration.jl src/rate_eq_derivation.jl src/types.jl | git apply
```

- [ ] **Step 4: Green** (background). Gt `$L/T46_green.log`
  `test_mechanism_enumeration.jl test_reaction_graph.jl`: 32,557, 956.
- [ ] **Step 5: Red on the K3 and one-index mutants** (S3g §9 task 4: the plain dead-end
  candidate's twins found on the parent; a gauge class keyed by one index). The
  prototype ran them at 7d69b90, after the split and tag ports; here the moves around them
  are older, so the expectation names the testsets that must fail and bounds the rest.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py k3
```

  Expected: `mutated src/mechanism_enumeration.jl (k3)`. Background: Rt
  `$L/T46_mut_k3.log` `test_mechanism_enumeration.jl`. Expected
  (`$P/proto/logs/S3g_mut_k3.log`: 1 fail, 6 errors): "Mechanism — a copy whose table
  entry sorts before the other reactants'" (0f795c4's new testset) errors, and every
  testset that fails or errors is one of the log's seven: that one, "Mechanism — ordered
  bi-bi with A as its own inhibitor: three placements", "Mechanism — shared A group with
  abortive E(A, Q): the gauge decides", "Mechanism — a second copy may share an older
  copy's composition", "Mechanism — a copy that only matches a form through a
  conformational isomer", "catalytic moves on bi-bi to depth 2: counts and both rules on
  every child" and "seed_mechanisms errors when no mechanism binds every required
  regulator". Then:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py k3 --revert
python3 $P/replay/mutate.py one_index
```

  Expected: `restored src/mechanism_enumeration.jl`, `mutated
  src/mechanism_enumeration.jl (one_index)`. Background: Rt `$L/T46_mut_class.log`
  `test_mechanism_enumeration.jl`. Expected (`$P/proto/logs/S3g_mut_class.log`: 4 fail, 1
  error): "_gauge_rescaling: group 2 and form 2 are two classes" (0f795c4's other new
  testset) fails twice, and nothing outside it and "catalytic moves on bi-bi to depth 2:
  counts and both rules on every child" fails or errors. Then, whatever the runs printed:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py one_index --revert
git status --short
```

  Expected: `restored src/mechanism_enumeration.jl`, then ` M` for the five paths of Step 6.
- [ ] **Step 6: Commit.** K `0f795c4` `src/mechanism_enumeration.jl
  src/rate_eq_derivation.jl src/types.jl test/test_mechanism_enumeration.jl
  test/test_reaction_graph.jl` ("Decide the copy rule and the dead-end move on the reaction
  graph").
- [ ] **Step 7: Full suite** (background). F 46: 39,980.
- [ ] **Step 8: Census** (background). C 46: zero differences.
- [ ] **Step 9: Copy-rule drops.** CD 46: the counts in the runner's line (split 27 + 13,
  dead-end 81 + 44, change-state 2).

### Task 47: Split on the graph [S3g]

**Replays:** `[5d1197d + 78d8376]`, `[7d69b90]` (S3g §9 task 5, which folds 7d69b90's
tied-group split testset into the split; 7d69b90 keeps its own commit, since it also adds
the bi-bi round-trip aggregate).

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_expand_split_kinetic_group`,
  `_context_bipartitions(m, g)`, `_bipartitioned_groups(m, selection)`, both
  `_split_gain_test`, `_apply_bipartitions(m, selection, re)`,
  `_revert_zero_flux_parts(m, bp, flux)`, `_context_form(e::Edge)`), `src/thermo.jl`
  (`_count_kind(e, re)`; the Step method goes)
- Test: `test/test_mechanism_enumeration.jl`, `test/test_oracles.jl`

**Interfaces:**
- Consumes: Tasks 45-46 (`_raw`, `_build`, `_redundant_copy_groups`); `_onlya_offenders(m,
  group, tags)` (Task 41); `_partition_independent_count`, `_independent_param_count`;
  `$P/replay/mutate.py` (Task 0 Step 8): `stale_direction`.
- Produces: a split unit is `(g, (part1, part2))` of edge lists; `_bipartitioned_groups`
  numbers each second part right after its group (D-S3g-5); the allosteric gain runs
  `_onlya_offenders` on the parent's graph under the candidate's groups and tags; every
  child goes through the checked constructor; `_count_kind(e, re = e.rapid_equilibrium)`.
  TE: "_expand_split_kinetic_group: a split member of a tied group turns alone" and "every
  bi-bi init mechanism and child is what its raw input builds" (239, 2,831, 1,833).

**Gate:** F 47 = 39,982 (5d1197d −7, 7d69b90 +9); C 47 (adapter): zero differences.

**Unit 1: 5d1197d + 78d8376** (78d8376 wraps TO's 94-column loop header that 5d1197d
writes, S3g round-1 B1)

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 5d1197d -- test/test_mechanism_enumeration.jl test/test_oracles.jl | git apply
git show --format= --binary 78d8376 -- test/test_oracles.jl | git apply
```

- [ ] **Step 2: Red** (background). Rt `$L/T47_1_red.log`
  `test_mechanism_enumeration.jl test_oracles.jl`. Expected
  (`$P/proto/logs/S3g_red_c5.log`): 9 errors, each a `MethodError` for
  `_context_bipartitions(m, g)` or the graph forms after it: TE 8, one each in
  "_expand_split_kinetic_group: a rejected split is the parent's model",
  "_expand_split_kinetic_group: a four-step group splits two and two",
  "_context_bipartitions", "_context_bipartitions separates an inhibitor-bound mirror",
  "_apply_bipartitions", "_revert_zero_flux_parts", "_partition_independent_count agrees
  with _independent_param_count" and "AllostericMechanism — the split's :OnlyA test agrees
  with _check_well_posed"; TO 1, in "a same-model regrouping keeps every fitted name".
- [ ] **Step 3: Apply the source** (78d8376 has none)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 5d1197d -- src/mechanism_enumeration.jl src/thermo.jl | git apply
```

- [ ] **Step 4: Green** (background). Gt `$L/T47_1_green.log`
  `test_mechanism_enumeration.jl test_oracles.jl`: 32,550 (−7: the reverted-parts testset
  checks the reverted edge list and the child's flags, 15 → 8; S3g §4), 667.
- [ ] **Step 5: Commit.** K `5d1197d` `src/mechanism_enumeration.jl src/thermo.jl
  test/test_mechanism_enumeration.jl test/test_oracles.jl` ("Split kinetic groups on the
  reaction graph"; the fold changes nothing the message says).

**Unit 2: 7d69b90** (test-only; green on arrival. Its tied-group testset fails on the
stale-direction mutant, 2 fail and 1 error, `$P/proto/logs/S3g_mut_tied.log`, which Step
8 reruns; the aggregate pins every bi-bi child's raw round trip)

- [ ] **Step 6: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 7d69b90 -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 7: Green** (background). Gt `$L/T47_2_green.log`
  `test_mechanism_enumeration.jl`: 32,559.
- [ ] **Step 8: Red on the stale-direction mutant** (S3g §9 task 5: the split's child
  built from the parent's stored edges, regrouped without the constructor)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py stale_direction
```

  Expected: `mutated src/mechanism_enumeration.jl (stale_direction)`. Background: Rt
  `$L/T47_2_mut.log` `test_mechanism_enumeration.jl`. Expected
  (`$P/proto/logs/S3g_mut_tied.log`): 2 fail and 1 error, all in
  "_expand_split_kinetic_group: a split member of a tied group turns alone", and no other
  testset fails or errors (S3g §1: only the new tied-group fixture catches it). Then,
  whatever the run printed:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/mutate.py stale_direction --revert
git status --short
```

  Expected: `restored src/mechanism_enumeration.jl`, then ` M
  test/test_mechanism_enumeration.jl` only.
- [ ] **Step 9: Commit.** K `7d69b90` `test/test_mechanism_enumeration.jl` ("Pin the turn of
  a tied member a split frees, and every bi-bi child's raw round trip").

**Task end**

- [ ] **Step 10: Full suite** (background). F 47: 39,982.
- [ ] **Step 11: Census** (background). C 47: zero differences.

### Task 48: Tag, relax and check on the graph [S3g]

**Replays:** `[2c40783 + 01e2dbd + 704bac2 + 9cece42 + 801cd5c]` plus 216be7f's ROADMAP
entry below (S3g §9 task 6, with round 1's m1, m4 and m5).

**Departure from S3g's round-1 note** (Open for Denis): the note writes each fix where its
code is first written (m1's comments in tasks 3 and 6, `_entry_names` and `_groups_with`
from task 3 on, m5 in tasks 4 and 6), so that no task writes code a later one rewrites.
The four commits land together here instead, so Tasks 45-47 write code this task rewrites.
Each commit edits sites across the code of S3g's tasks 3, 4 and 6 in the same files
(`src/mechanism_enumeration.jl`, `src/types.jl`, `src/thermo.jl`, `src/names.jl`,
`src/rate_eq_derivation.jl`, `src/admission.jl`), so following the note means splitting
them hunk by hunk into written-out patches and building intermediate trees the prototype
never built or tested. After 2c40783 they apply whole, in the prototype's order, and the
prototype ran the full suite on each (704bac2, 9cece42, 801cd5c: 39,942).

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_expand_to_allosteric`,
  `_expand_change_allo_state`, `_assert_emission_rules`, `_assert_atom_conserving`,
  `_check_step_atoms`, `_bound_comp_inhibitors`, `_dead_end_raw`, the plain
  `_dead_end_child`), `src/types.jl` (`_chemistry_groups`, `_with`, `_entry_names`,
  `_groups_with`; `_allosteric_of` inlined), `src/thermo.jl` (`_onlya_offenders`; the Step
  `_onlya_haldane_violation(m, groups, tags)` and `_thermodynamic_constraints(m, groups)`
  go), `src/names.jl` (`_forward_sides`, `_side_label` go), `src/admission.jl`
  (`_bottomless_re_segment`'s docstring), `src/rate_eq_derivation.jl` (`_entry_names` at the
  flux terms), `ROADMAP.md`
- Test: `test/test_mechanism_enumeration.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: Tasks 45-47; `_onlya_offenders(m, group, tags)` (Task 41).
- Produces: `_chemistry_groups(m, group, n)`; `_with(am; states, sites)`, the composition
  constructor and `_with`'s only form; `_entry_names(m, c) -> Vector{Symbol}` (the names a
  count holds, at the nine sites S3g round-1 m4 lists); `_groups_with(pred, group, n) ->
  BitVector` (behind `_chemistry_groups`, `_flux_carrying_groups` and the flip's
  modifier-exchange test); `_dead_end_raw` extends the vectors `_raw` made fresh. Gone:
  the Step `_onlya_haldane_violation(m, groups, tags)`, the Step
  `_thermodynamic_constraints(m, groups)`, `_forward_sides`, `_side_label`,
  `_allosteric_of`. TT's `:OnlyA` validator helpers read `_onlya_offenders` alone.

**Gate:** F 48 = 39,968 (TT −14: the agreement assertions against the deleted Step method,
S3g §4); C 48 (adapter): zero differences; CD 48.

- [ ] **Step 1: Apply the test paths** (the four fold commits have none)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2c40783 -- test/test_mechanism_enumeration.jl test/test_types.jl | git apply
```

- [ ] **Step 2: Green before the source** (background; a refactor, so the new tests pass
  on Task 47's source, as on 2c40783's parent: `$P/proto/logs/S3g_red_c6.log`, 33,897).
  Gt `$L/T48_before.log` `test_mechanism_enumeration.jl test_types.jl test_oracles.jl`:
  32,559, 680, 667.
- [ ] **Step 3: Apply the source, in commit order**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2c40783 -- src/mechanism_enumeration.jl src/names.jl src/thermo.jl src/types.jl | git apply
git show --format= --binary 01e2dbd -- src/admission.jl src/mechanism_enumeration.jl src/thermo.jl | git apply
git show --format= --binary 704bac2 -- src/mechanism_enumeration.jl src/names.jl src/rate_eq_derivation.jl src/types.jl | git apply
git show --format= --binary 9cece42 -- src/mechanism_enumeration.jl src/types.jl | git apply
git show --format= --binary 801cd5c -- src/mechanism_enumeration.jl src/types.jl | git apply
```

- [ ] **Step 4: Departure: the ROADMAP entry** (216be7f's piece; S3g section 8 lists it).
  Set "The `:OnlyA` guard misjudges groups that mix a binding and a chemistry step." to
  its final text, the guard's test named as the code now holds it:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/roadmap_entry.py 'The `:OnlyA` guard misjudges' $P/replay/roadmap_final.md ROADMAP.md
git diff --stat -- ROADMAP.md
```

  Expected: `ROADMAP.md | 5 +++--`, and ENT with the title `` The `:OnlyA` guard
  misjudges `` prints nothing. The entry then reads:

```markdown
- **The `:OnlyA` guard misjudges groups that mix a binding and a chemistry step.**
  The guard's test, `_onlya_offenders` (`src/thermo.jl`), which the boundary check runs
  through `_onlya_haldane_violation`, skips every `:OnlyA` group that holds a chemistry
  step (`_chemistry_groups`), including one that also holds a
  plain binding of the same metabolite (enumeration builds these by merging and tags
  every allosteric one `:OnlyA`), and gives the bindings of one group a single rate of
  divergence. Hand-built mechanisms show both errors: it admits a mixed `:OnlyA` group
  whose plain member sits on a square of `:EqualAI` bindings, whose inactive conformation
  is not a limit of the model, and rejects a split `:OnlyA` binding group whose limit
  exists. The `:OnlyA` semantics under Planned replace the guard.
```

- [ ] **Step 5: Green** (background). Gt `$L/T48_green.log`
  `test_mechanism_enumeration.jl test_types.jl test_oracles.jl test_rate_eq_derivation.jl`:
  32,559, 680, 667, 3,197.
- [ ] **Step 6: Commit** (message: 2c40783's, plus one paragraph for the folds and one for
  the entry)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M ROADMAP.md, src/admission.jl, src/mechanism_enumeration.jl, src/names.jl, src/rate_eq_derivation.jl, src/thermo.jl, src/types.jl, test/test_mechanism_enumeration.jl, test/test_types.jl
git add ROADMAP.md src/admission.jl src/mechanism_enumeration.jl src/names.jl src/rate_eq_derivation.jl src/thermo.jl src/types.jl test/test_mechanism_enumeration.jl test/test_types.jl
cat > $P/replay/msg <<'EOF'
Tag, relax and check mechanisms on the reaction graph

The allosteric moves read Step views for their chemistry groups and ran the Step
form of the :OnlyA test on them. To-allosteric and change-state now read
_chemistry_groups and _onlya_offenders on the catalytic graph, and every
allosteric move builds its children through the composition constructor (_with
keeps only its tags-and-sites form). The emission rules name groups by their
edges' sides, the atom check reads each edge's forms and lists, and the bound
competitive inhibitors are the modifier entries some edge takes up. The Step forms
of the :OnlyA test and the cycle basis, _forward_sides and _side_label lose their
last callers and go; the two atom checks, the seed bases' and a mechanism's, share
one message.

_entry_names(m, c) reads the names a count holds, at the nine sites that spelled
the comprehension, and _groups_with(pred, group, n) marks the groups that hold an
edge meeting a predicate, for _chemistry_groups, _flux_carrying_groups and the
flip's modifier-exchange test. _dead_end_raw extends the vectors _raw made fresh
instead of copying them, the plain _dead_end_child leaves its unread regulator
argument unnamed, and _allosteric_of, left with one caller, is inlined there. The
comments on _bottomless_re_segment and _all_steady_state name the seed variants'
Step _minimal_flips as their reader, and those on the Species content keys name
_catalytic_topologies and _forms, which order forms by them.

ROADMAP: the misjudged-groups entry names the guard's test, _onlya_offenders,
which the boundary check runs through _onlya_haldane_violation, and
_chemistry_groups.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

- [ ] **Step 7: Full suite** (background). F 48: 39,968.
- [ ] **Step 8: Census** (background). C 48: zero differences.
- [ ] **Step 9: Copy-rule drops.** CD 48: the counts in the CD runner's line (split 27 +
  13, dead-end 81 + 44, change-state 2), unchanged from Task 46.

### Task 49: Hotspot fixtures inline; the census population freed [S3g]

**Replays:** `[612a678 + 5fc1f44]` with R55 m2's indentation, `[698fd6a]` (S3g §9 tasks 7-8;
S3g round-1 m3 is 5fc1f44).

**Files:**
- Test: `test/test_mechanism_enumeration.jl`, `test/test_reaction_graph.jl`

**Interfaces:**
- Consumes: the moves of Tasks 45-48.
- Produces: TE's hotspot testsets build every parent and every expected child as an
  `@enzyme_mechanism` or `@allosteric_mechanism` block and assert the child count and set
  (CLAUDE.md enumeration rules 1-3; S3g §4 lists them); `_testhelper_mkm` goes. TRG builds
  its census population in a `let` around the two testsets that read it, with their
  helpers above it.

**Gate:** F 49 = 39,989 (612a678 +21); R52's three figures recorded; checkpoint S3g.

**Unit 1: 612a678 + 5fc1f44, R55 m2** (test-only: every expected child was rendered from
the move's output and shown `==` to it, `$P/proto/logs/S3g_probe_hotspots.log`; green on
arrival)

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 612a678 -- test/test_mechanism_enumeration.jl | git apply
git show --format= --binary 5fc1f44 -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 2: Departure: R55 m2's indentation** (D-S3i-10; the S3g review's m2: "Indent
  the three `EnzymeRates.AllostericMechanism(...)` elements and the closing `]` by four
  spaces, as the first-regulator testset does"). In TE's testset "AllostericMechanism —
  uni-uni all-:EqualAI: 3 cat relaxations", 5fc1f44 wrote the three expected children at
  the indentation of `expected = [`; the script prefixes four spaces to every line after
  `    expected = [` up to the closing `    ]`, which stays at four spaces like every
  other `expected = [` list in TE (R63 item 3). No later commit touches those lines.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
sed -n '6673p;6696p;6724p' test/test_mechanism_enumeration.jl
python3 $P/replay/r55m2.py test/test_mechanism_enumeration.jl
diff <(git show 5fc1f44:test/test_mechanism_enumeration.jl) test/test_mechanism_enumeration.jl | grep -c '^>'
diff -w <(git show 5fc1f44:test/test_mechanism_enumeration.jl) test/test_mechanism_enumeration.jl && echo whitespace-only
```

  Expected: the testset's title line, `    expected = [` and `    ]`; then
  `indented lines 6697-6723`; then `27`; then `whitespace-only`. The file now equals
  5fc1f44's but for those 27 lines.
- [ ] **Step 3: Green** (background). Gt `$L/T49_1_green.log`
  `test_mechanism_enumeration.jl`: 32,580.
- [ ] **Step 4: Commit** (message: 612a678's, plus one paragraph for 5fc1f44 and the
  indentation)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M test/test_mechanism_enumeration.jl
git add test/test_mechanism_enumeration.jl
cat > $P/replay/msg <<'EOF'
Write the enumeration hotspot fixtures inline with their exact children

The move tests that built their parents from an init seed and today's Step
constructors, or that asserted properties of children only, now write every parent
and every expected child as an @enzyme_mechanism or @allosteric_mechanism block and
assert the child count and set (CLAUDE.md enumeration-test rules): the paired flip
of an inhibitor-bound mirror, the allosteric split, to-allosteric with and without
a bound competitive inhibitor, a first regulator, the skipped OnlyI append, the
change-state relaxations and the merges. The _respects_reg_type, canonical-form
and dedup-key fixtures become blocks too, and the remaining Step reads of these
testsets read the graph.

The first-regulator testset pairs each reaction with its three children and the
relaxation testset lists its three, each written out as a block rather than built
from one eval'd template with interpolated tags and ligands; the relaxation
testset's children sit one level inside their list, as the other lists' do.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

**Unit 2: 698fd6a** (R54 m5: test-only, green on arrival; frees TRG's population after
its two testsets)

- [ ] **Step 5: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 698fd6a -- test/test_reaction_graph.jl | git apply
```

- [ ] **Step 6: Green.** Gt `$L/T49_2_green.log` `test_reaction_graph.jl`: 956.
- [ ] **Step 7: Commit.** K `698fd6a` `test/test_reaction_graph.jl` ("Release the graph
  tests' census population after its two testsets").

**Task end**

- [ ] **Step 8: Full suite** (background). F 49: 39,989.
- [ ] **Step 9: R52's figures** (record; the gate is Task 55's). SC 49, then TEA 49, both
  background. Record F 49's maxrss, the ter-ter time and `TE_ELAPSED`. The prototype's S3g
  ran without 050e16f's copies (4.17-4.20 GB, 35.6-37.5 s, 117.8-136.4 s); this branch
  carries them since Task 34, so expect lower figures, near S3i's for memory. Seed counts
  264 / 239 / 7 / 250,855.
- [ ] **Step 10: Checkpoint S3g** (tip 5fc1f44). Run
  `git diff --stat 5fc1f44 HEAD -- src test docs ROADMAP.md Project.toml .claude`.
  Expected: exactly these 19 files, with these line counts (the dry replay in a scratch
  index, which matched `$S/p3skel/sim_final.log` but for ROADMAP.md, whose final
  newline that simulation dropped): Project.toml 2, ROADMAP.md 48,
  docs/src/deriving/cha_king_altman.md 3, docs/src/deriving/mwc_allostery.md 8,
  docs/src/deriving/ping_pong.md 2, docs/src/deriving/re_vs_ss.md 6,
  docs/src/deriving/regulators.md 6, docs/src/developer.md 6,
  docs/src/getting_started.md 4, docs/src/identify/enumeration_engine.md 5,
  docs/src/identify/tutorial.md 7, the plan (D2), the spec 48 (D1), src/dsl.jl 4,
  src/types.jl 5, FX 35, TE 54, TRD 31, TRG 29. Then check that each holds only its named
  changes:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
F=$P/replay/roadmap_final.md
$P/replay/ckpt_p4.sh 5fc1f44 c:2b29fd2 c:be65542 c:0b2021d c:a84d194 \
    p:216be7f:docs/src/deriving/mwc_allostery.md p:216be7f:docs/src/deriving/cha_king_altman.md \
    p:216be7f:docs/src/getting_started.md p:216be7f:src/dsl.jl c:050e16f \
    p:3962124:test/test_reaction_graph.jl \
    "e:Encode a mechanism's Sig=delete" "e:Have the \`Mechanism\` constructors check=delete" \
    "e:Take \`:OnlyA\` steps to their limit=$F" "e:\`eq_complexity_filter\` counts only=$F" \
    "e:The \`:OnlyA\` guard misjudges=$F" i
```

  Expected: no output. Departures in S3g: R55 m2's 27 lines (Task 49); the ROADMAP entry
  "The `:OnlyA` guard misjudges ..." (Task 48); the messages of Task 48 and of Task 49's
  first unit, written out there. Ahead of the tip, as the Stage checkpoints table lists: 050e16f
  (src/types.jl, TRG), 3962124 (TRG), the S1-S2 docs commits, 216be7f's mwc_allostery.md,
  cha_king_altman.md, getting_started.md and src/dsl.jl pieces, and the ROADMAP entries of
  Tasks 34, 36 and 42.

### Task 50: Seeds from raw input [S3h]

**Replays:** `[9dee115]` (S3h §9 task 1).

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_catalytic_topologies`,
  `_expand_substrate_product_dead_ends`, `_seed_groups`, `init_mechanisms`),
  `src/thermo.jl` (the Species-key header), `src/types.jl` (`_step_kind` goes)
- Test: `test/test_mechanism_enumeration.jl`

**Interfaces:**
- Consumes: `RawForm`, `RawStep`, `RawMechanism`, `Met` (Task 5), with content `==` and
  `hash`; `Mechanism(raw)` (Task 34).
- Produces: `_catalytic_topologies(rxn) -> Vector{Vector{RawStep}}`, a release written as
  the binding it reverses (`RawStep(to, cur, [Met(p, false)], Met[], true)`), patterns and
  steps sorted by a raw-form key in the canonical form order (D-S3h-1);
  `_expand_substrate_product_dead_ends` on raw topologies (forms in first-seen order,
  D-S3h-2); `_seed_groups(steps) -> Vector{Int}` (labels by first use); each seed
  `Mechanism(RawMechanism(subs, prods, Symbol[], steps, _seed_groups(steps)))`. TE:
  `_testhelper_iso_orient`, `_testhelper_topo_mech(rxn, t)`, `_testhelper_form_names(rxn,
  t)`, `_testhelper_flat_topo`, `_testhelper_raw_iso`, `_testhelper_raw_residual`.

**Gate:** F 50 = 39,991; C 50 (adapter): zero differences, order included; SC 50; the
init golden (TE "bi-bi init_mechanisms structural golden") passes in the suite.

- [ ] **Step 1: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9dee115 -- test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 2: Red** (background). Rt `$L/T50_red.log` `test_mechanism_enumeration.jl`.
  Expected (`$P/proto/logs/S3h_redA2.log`): 32,535 pass, 1 fail and 35 errors, all in the
  topology, dead-end and seed-group testsets: the fail in "a topology writes each release
  as the binding it reverses"; errors "topology counts" 13, "isomerization constraints" 3,
  "pyruvate carboxylase mechanism" 3, "pyruvate dehydrogenase mechanism" 3, "quad-quad: C6
  forces ping-pong" 2, "admissible-residual ping-pong" 2, "Ter-ter per-topology (OOM on
  full init)" 2, and 1 each in "Bi-Bi Ping-Pong: 6 dead-end forms → 7 variants", "Bi-Bi
  random: 4 dead-end forms", "Ter-ter random: 27 dead-end forms, 12 on the diagonal
  pattern", "Uni-Bi ordered: no dead-end forms", "Uni-Uni: no dead-end forms",
  "_seed_groups labels one metabolite's bindings at one flag as one group" and
  "isomerization patterns come in the canonical form order".
- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 9dee115 -- src/mechanism_enumeration.jl src/thermo.jl src/types.jl | git apply
```

- [ ] **Step 4: Green** (background). Gt `$L/T50_green.log` `test_mechanism_enumeration.jl`:
  32,582 (`$P/proto/logs/S3h_greenA.log`).
- [ ] **Step 5: Commit.** K `9dee115` `src/mechanism_enumeration.jl src/thermo.jl
  src/types.jl test/test_mechanism_enumeration.jl` ("Build the seeds from raw input").
- [ ] **Step 6: Full suite** (background). F 50: 39,991.
- [ ] **Step 7: Census** (background). C 50: zero differences, init and seed order
  included (init 1, 239, 191, 1, 4,885, 191; seeds 16, 1,355, 1,728, 4, none, 1,058).
- [ ] **Step 8: Seed counts** (background). SC 50: 264 / 239 / 7 / 250,855.

### Task 51: Seed variants and degeneracy on the graph [S3h]

**Replays:** `[2395401 + 216be7f: docs/src/identify/enumeration_engine.md]` (S3h §9 task 2;
S5 §9 item 3's S3h and S3a pieces of the engine page).

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_merge_isomerization`, `_eliminate_form`,
  `_re_turnover_cycle`, `_has_vmax`, `_fused_substrate_binding`, `_crosses_sides`,
  `_chemistry_equilibrates_both_sides`, `_degenerate`, `_seed_candidate_screen`,
  `_seed_variants`, `_entry_signs`, `_uptake_weight(e, sign)`, `_assert_atom_conserving`;
  the Step helpers go), `src/types.jl` (`_entry_mask`, `_names_mask`; `_canonical_order`
  goes), `src/admission.jl` (`_bottomless_re_segment` goes),
  `docs/src/identify/enumeration_engine.md`
- Test: `test/test_mechanism_enumeration.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: Task 50's raw seeds; `_canonicalize` and `_check_aligned` (K12); the graph
  `_minimal_flips`, `_flux_carrying_edges` (Task 45); `_segments`, `_bottomless_segment`.
- Produces: `_merge_isomerization(raw, s0)` and `_eliminate_form(raw, x)` on raw input
  (labels renumbered by first use, D-S3h-3); the degeneracy predicates take a graph and a
  flag vector (D-S3h-4); `_seed_candidate_screen(m)` returns `re -> Bool`; `_seed_variants`
  builds each variant `Mechanism(_raw(base, re))`; `_entry_signs`, `_uptake_weight(e,
  sign)`; `_entry_mask(pred, m)`, `_names_mask(m, names)`. Gone: `_canonical_order`,
  `_bottomless_re_segment`, `_with_equilibrium`, `_all_steady_state`, the Step
  `_minimal_flips`, `_flux_carrying_steps`, the Step `_flux_carrying_groups`,
  `_reactant_signs`, the Step `_uptake_weight`, `_any_named`,
  `_assert_step_atom_conserving`, `_check_step_atoms`. TE: `_testhelper_binds(m, xs...)`,
  `_testhelper_is_tc`, `_testhelper_flip_flags`, `_testhelper_flip_matching`,
  `_testhelper_all_re`, `_testhelper_seed_bases`, `_testhelper_base` (Task 55 replaces
  `_testhelper_all_re` and `_testhelper_seed_bases` with `_seed_bases`).

**Gate:** F 51 = 39,997; C 51 (adapter): zero differences, order included; SC 51; R52's
figures recorded; the R18 probe; checkpoint S3h.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2395401 -- test/test_mechanism_enumeration.jl test/test_types.jl | git apply
```

- [ ] **Step 2: Red** (background). Rt `$L/T51_red.log`
  `test_mechanism_enumeration.jl test_types.jl`. Expected (`$P/proto/logs/S3h_redB.log`):
  TE 32,362 pass and 24 errors, every one in the rewritten seed testsets: "C on the
  ping-pong seed and its single flips" 10, "a reactant's modifier copy counts as its
  reactant in the seed tests" 3, "MERGE, ELIM and the seed tests on the ordered bi-bi seed"
  2, "_re_turnover_cycle: a balanced rapid-equilibrium cycle runs no turnover" 2,
  "_seed_variants: the uni-uni seed has none" 2, and 1 each in "ELIM refuses the uni-uni
  merged complex", "_seed_candidate_screen agrees with the group predicates on every
  candidate", "_seed_variants: an ordered seed whose B group is shared", "_seed_variants:
  the ordered/random seed" and "_seed_variants: the undecorated ping-pong seed"; TT passes,
  680 (the log ran TE alone; TT's rewritten bottomless test reads `_bottomless_segment`,
  which Task 41 wrote, and the Step-input `Mechanism(rxn, groups)`, which lives until
  Task 54).
- [ ] **Step 3: Apply the source and the engine page**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2395401 -- src/admission.jl src/mechanism_enumeration.jl src/types.jl | git apply
git show --format= --binary 216be7f -- docs/src/identify/enumeration_engine.md | git apply
```

  `git diff --stat -- docs` shows `docs/src/identify/enumeration_engine.md | 8 +++++---`.
- [ ] **Step 4: Green** (background). Gt `$L/T51_green.log`
  `test_mechanism_enumeration.jl test_types.jl test_identify_rate_equation.jl`: 32,588,
  680, 446 (`$P/proto/logs/S3h_greenB.log`).
- [ ] **Step 5: Commit** (message: 2395401's, plus one paragraph for the engine page)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M docs/src/identify/enumeration_engine.md, src/admission.jl, src/mechanism_enumeration.jl, src/types.jl, test/test_mechanism_enumeration.jl, test/test_types.jl
git add docs/src/identify/enumeration_engine.md src/admission.jl src/mechanism_enumeration.jl src/types.jl test/test_mechanism_enumeration.jl test/test_types.jl
cat > $P/replay/msg <<'EOF'
Build the seed variants and judge degeneracy on the reaction graph

Merge and elimination run on the seed's raw input, and each base is read as
the graph _canonicalize gives it, unchecked, as bases were never built. The
candidate screen and the minimal-flip search read the base with flag
overrides, and each variant is built through the checked constructor; the
lumping twin finds the merged complexes in the base through form_of.
_degenerate and its tests (turnover cycle, maximal rate, chemistry node) read
the catalytic graph with its flags, a modifier copy counting as its reactant
by name as before. The Step copies of the minimal-flip search, the flux test,
the bottomless test, the canonical-order adapter and the seed atom check go
with their last caller.

Docs: the enumeration engine page says the seed tests run on a candidate's
flags over the base's graph, names `_bottomless_segment` and the boundary check
`_check_well_posed` as the rule they apply, and says the ligand move adds one of
the reaction's allosteric regulators, which a regulatory site stores by name.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

- [ ] **Step 6: Full suite** (background). F 51: 39,997.
- [ ] **Step 7: Census** (background). C 51: zero differences, order included.
- [ ] **Step 8: R52's figures** (record). SC 51 (264 / 239 / 7 / 250,855), then TEA 51,
  both background. The prototype's S3h, without 050e16f's copies: 4.22 GB, 20.8 s, 105.3 s.
- [ ] **Step 9: R18 probe.** PR `$L/T51_r18.log` `$P/s3h/probe/r18.jl`. Expected
  (`$P/proto/logs/S3h_r18_B.log`): 25 lines `concrete ...` and one `NOT`,
  `_eliminate_form(EnzymeRates.RawMechanism, EnzymeRates.RawForm) ->
  Union{Nothing, EnzymeRates.RawMechanism}` (a non-parametric struct beside `Nothing`,
  which its one caller filters once; R18 allows it).
- [ ] **Step 10: Checkpoint S3h** (tip 2395401).
  `git diff --stat 2395401 HEAD -- src test docs ROADMAP.md Project.toml .claude` lists
  exactly these 19 files: Project.toml 2, ROADMAP.md 48,
  docs/src/deriving/cha_king_altman.md 3, docs/src/deriving/mwc_allostery.md 8,
  docs/src/deriving/ping_pong.md 2, docs/src/deriving/re_vs_ss.md 6,
  docs/src/deriving/regulators.md 6, docs/src/developer.md 6,
  docs/src/getting_started.md 4, docs/src/identify/enumeration_engine.md 13,
  docs/src/identify/tutorial.md 7, the plan (D2), the spec 48 (D1), src/dsl.jl 4,
  src/types.jl 5, FX 35, TE 54, TRD 31, TRG 29. Then:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
F=$P/replay/roadmap_final.md
$P/replay/ckpt_p4.sh 2395401 c:2b29fd2 c:be65542 c:0b2021d c:a84d194 \
    p:216be7f:docs/src/deriving/mwc_allostery.md p:216be7f:docs/src/deriving/cha_king_altman.md \
    p:216be7f:docs/src/getting_started.md p:216be7f:src/dsl.jl c:050e16f \
    p:3962124:test/test_reaction_graph.jl \
    "e:Encode a mechanism's Sig=delete" "e:Have the \`Mechanism\` constructors check=delete" \
    "e:Take \`:OnlyA\` steps to their limit=$F" "e:\`eq_complexity_filter\` counts only=$F" \
    "e:The \`:OnlyA\` guard misjudges=$F" i p:216be7f:docs/src/identify/enumeration_engine.md
```

  Expected: no output. Departures in S3h: Task 51's message; ahead of the tip, S3g's list
  (3962124's TRG piece included) plus 216be7f's enumeration_engine.md.

### Task 52: The compiled types' summary from the graph [S3i]

**Replays:** `[6ada4b2]` (S3i §9 task 2; 050e16f, S3i's task 1, is in Task 34).

**Files:**
- Modify: `src/names.jl` (`_side_text(m, f, c)`, `_edge_text(m, k)`), `src/types.jl`
  (`show(::EnzymeMechanism)`, `show(::AllostericEnzymeMechanism)`)

**Interfaces:**
- Consumes: `_flat_edges`, `_group_members`, `_form_name`, `_entry_label`.
- Produces: `_side_text(m, f, c)`, `_edge_text(m, k)` (an edge as the summary prints it);
  both `show` methods read the graph. Task 59 deletes all four when `show` prints DSL text.

**Gate:** F 52 = 39,997; C 52 (adapter): zero differences.

- [ ] **Step 1: Apply the source** (no test path; TT "Pretty printing" pins the printed
  text, green before, at Task 51, and after: `$P/proto/logs/S3i_showA.log`)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 6ada4b2 -- src/names.jl src/types.jl | git apply
```

- [ ] **Step 2: Green.** Gt `$L/T52_green.log` `test_types.jl test_oracles.jl`: 680, 667.
- [ ] **Step 3: Commit.** K `6ada4b2` `src/names.jl src/types.jl` ("Print the compiled
  types' summary from the reaction graph").
- [ ] **Step 4: Full suite** (background). F 52: 39,997.
- [ ] **Step 5: Census** (background). C 52: zero differences.

### Task 53: Every test reads the graph [S3i]

**Replays:** `[36384ef + a8f6f51 + 3962124: FX + 4952d79: TRG + b42cc90: test/test_dsl.jl
test/test_mechanism_enumeration.jl]` (S3i §9 task 3 with the round-1 items its note puts
here: m7 is 3962124's FX piece, the reader's step order; B2 and m2 are b42cc90's TE piece,
m4 its TD piece. m3 sits in b42cc90's TT piece beside m1, which needs Task 55's
`_is_iso(::RawStep)`, so it lands there. R58 m2 is a8f6f51 and R58 m1 4952d79's TRG piece,
both S5 §9 item 6.)

**Files:**
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/test_accessors.jl`,
  `test/test_allosteric_collapse.jl`, `test/test_compile_budget.jl`, `test/test_dsl.jl`,
  `test/test_joint_solve.jl`, `test/test_mechanism_enumeration.jl`, `test/test_oracles.jl`,
  `test/test_rate_eq_derivation.jl`, `test/test_reaction_graph.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: the graph's read API (`_flat_edges`, `entries`, `_form_name`, `_is_free`,
  `_group_members`, `_bound_entry`, `_is_binding`, `_is_iso`, `_is_chemistry(m, e)`), the
  composition constructor, `_raw_mechanism`, Task 52's `_edge_text`.
- Produces: FX's `_testhelper_graph` reads the graph's fields (same fields and values);
  `_testhelper_groups(g)`, `_testhelper_re_segments(g, kept)` (the segment oracle,
  D-S3i-6); `_testhelper_edge_text` calls `_edge_text`; `_testhelper_placeholder_reaction`
  goes, and its docstring gives its step order (group by group, each group in edge
  order). TD: `_testhelper_binds_entry(x, t)`, which `_testhelper_binding_of` filters
  with; TE's inline Theorell–Chance and steady-state group counts read `_testhelper_is_tc`
  and `_testhelper_group_flags`. TRG: `_testhelper_inactive_reach` (on the reader's steps;
  its docstring names its argument `g`) in place of
  `_testhelper_inactive_steps`; the Step-view comparisons and the test-side converter go
  (D-S3i-5). TT: `_testhelper_uniuni_steps`, `_testhelper_raw_uniuni`,
  `_testhelper_tagged_uniuni(steps, tags, n, sites)` (a8f6f51); `_testhelper_sp` and
  `_testhelper_uniuni` go. TE: `_testhelper_group_flags(m)`. TD: `_testhelper_form_roles`,
  `_testhelper_role`, `_testhelper_binding_of`. Every Step-input construction in the tests
  is raw input, the composition constructor or a macro block.

**Gate:** F 53 = 39,933 (−64: S3i §4 lists each removed representation test and the test
that covers its subject). No census (no src).

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 36384ef -- test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_accessors.jl test/test_allosteric_collapse.jl test/test_compile_budget.jl test/test_dsl.jl test/test_joint_solve.jl test/test_mechanism_enumeration.jl test/test_oracles.jl test/test_rate_eq_derivation.jl test/test_reaction_graph.jl test/test_types.jl | git apply
git show --format= --binary a8f6f51 -- test/test_types.jl | git apply
git show --format= --binary 3962124 -- test/mechanism_definitions_for_test_enzyme_derivation.jl | git apply
git show --format= --binary 4952d79 -- test/test_reaction_graph.jl | git apply
git show --format= --binary b42cc90 -- test/test_dsl.jl test/test_mechanism_enumeration.jl | git apply
```

- [ ] **Step 2: Green** (background; test-only, green on Task 52's source:
  `$P/proto/logs/S3i_portB.log`, 38,941 over these ten files with the prototype's TRD).
  Gt `$L/T53_green.log` `test_accessors.jl test_allosteric_collapse.jl test_compile_budget.jl test_dsl.jl test_joint_solve.jl test_mechanism_enumeration.jl test_oracles.jl test_rate_eq_derivation.jl test_reaction_graph.jl test_types.jl`:
  155, 36, 16, 254, 455, 32,582, 667, 3,197, 948, 633 (38,943).
- [ ] **Step 3: Commit** (message: 36384ef's, plus one paragraph for a8f6f51 and one for
  the 3962124, 4952d79 and b42cc90 pieces)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M the eleven test files of Step 1
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_accessors.jl test/test_allosteric_collapse.jl test/test_compile_budget.jl test/test_dsl.jl test/test_joint_solve.jl test/test_mechanism_enumeration.jl test/test_oracles.jl test/test_rate_eq_derivation.jl test/test_reaction_graph.jl test/test_types.jl
cat > $P/replay/msg <<'EOF'
Read the reaction graph in every test, and build test input as raw input

The tests read mechanisms through the Step views (steps, Step, Species)
and built fixtures through the Step-input constructors, which go with the
views. _testhelper_graph now reads the graph's fields; the view readers
read edges, forms and groups; the Step-input constructions are raw input,
the composition constructor or macro blocks. The Step, Species and
Residual representation tests go: the graph's own tests cover each
subject. The segment and inactive-mask oracles, which read the views, are
rewritten on the test reader's steps, so they stay independent of
_segments and _inactive_mask. The test-side converter (Step views to raw)
and the population differential it fed go; the canonical golden and the
raw and Sig round trips over the census population stay.

The two allosteric uni-uni testsets of test_types.jl build their mechanisms
with one header helper, _testhelper_tagged_uniuni, in place of a closure each,
and the four validation cases pin the message of the check each names, where
they asserted only an ErrorException.

The reader's docstring gives its step order (group by group, each group in
edge order), and the inactive-reach oracle's names its argument g. Inline
copies of the Theorell–Chance edge test and the steady-state group count read
the test helpers, and the dual-role DSL test counts the edges that bind R as a
modifier through one predicate, _testhelper_binds_entry.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

- [ ] **Step 4: Full suite** (background). F 53: 39,933.

### Task 54: Delete the Step views, Step, Species and Residual [S3i]

**Replays:** `[2692bb7 + 3962124: test/test_types.jl + fbdb9aa + 216be7f:
docs/src/deriving/re_vs_ss.md]` plus the ROADMAP entry below (S3i §9 task 4 with round 1's
B1, 3962124's TT piece; its m6 and m7 pieces are in Tasks 34 and 53, where the S3i note
puts them; S5's fbdb9aa; 216be7f's S3i pieces).

**Files:**
- Modify: `src/types.jl` (`Residual`, `Species`, `Step`, the views, the Step-input
  constructors, `_raw_of`, `_forms`, `_re_segment_extras`, `_flat_steps` and their helpers
  go), `src/names.jl` (`_is_reactant_copy`, `_met_label`, `_form_label`, `_step_text` go),
  `src/thermo.jl` (the Species content keys go), `docs/src/deriving/re_vs_ss.md`,
  `ROADMAP.md`
- Test: `test/test_types.jl` (one testset's title and comments),
  `test/test_mechanism_enumeration.jl` and `test/test_rate_eq_derivation.jl` (one comment
  each)

**Interfaces:**
- Consumes: Tasks 50-53 (no src code or test reads a view).
- Produces: a mechanism is its reaction graph alone. Gone from src: everything the
  Files line lists; `steps(m)`, `_view_*`, both `_raw_of`, `Mechanism(rxn, groups)`,
  `AllostericMechanism(rxn, groups, tags, n, sites)`, `_is_reactant(::EnzymeReaction, ...)`,
  `_is_reactant(::AllostericMechanism, ...)`. Task 59 writes a new `_step_text(m, k)` and a
  new `_is_reactant_copy(m, i)`, both on the graph.

**Gate:** F 54 = 39,933; C 54 without the adapter: zero differences; the grep of Step 4.

- [ ] **Step 1: Apply the test paths** (a title and comments)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 3962124 -- test/test_types.jl | git apply
git show --format= --binary fbdb9aa -- test/test_mechanism_enumeration.jl test/test_rate_eq_derivation.jl | git apply
```

- [ ] **Step 2: Apply the source and the page** (a deletion: the tests pass before, at
  Task 53, and after, `$P/proto/logs/S3i_delC.log`)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 2692bb7 -- src/names.jl src/thermo.jl src/types.jl | git apply
git show --format= --binary 216be7f -- docs/src/deriving/re_vs_ss.md | git apply
git diff --stat -- src
```

  Expected: `src/names.jl | 43`, `src/thermo.jl | 30`, `src/types.jl | 394`, deletions
  but for 10 inserted lines (2692bb7: +10 −457).
- [ ] **Step 3: Departure: the ROADMAP entry.** Set "Enumerate Cleland iso mechanisms." to
  its final text: 216be7f's "A `Form` carries its conformation label" with 7205ba4's
  sentence "Revisit allosteric ping-pong with this work ...", which 216be7f's version
  lacks.

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/roadmap_entry.py 'Enumerate Cleland iso' $P/replay/roadmap_final.md ROADMAP.md
git diff --stat -- ROADMAP.md
```

  Expected: `ROADMAP.md | 4 ++--`, and ENT `Enumerate Cleland iso` prints nothing. The
  entry then reads:

```markdown
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
```

- [ ] **Step 4: No view left in src**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git grep -nE '\bStep\b|\bSpecies\b|\bResidual\b|(^|[^_a-z])steps\(|_view_|_raw_of\b|_flat_steps|(^|[^_a-z])_forms\(|_re_segment_extras|_step_text|_form_label' -- src
```

  Expected: exactly the three step-level messages D18 keeps, `src/types.jl:1069`,
  `:1073` and `:1077` (`error("Step: both ends ...`, `error("Step $(_form_name(...)) → ...`
  twice); `RawStep` does not match the pattern.
- [ ] **Step 5: Green** (background). Gt `$L/T54_green.log`
  `test_accessors.jl test_allosteric_collapse.jl test_compile_budget.jl test_dsl.jl test_joint_solve.jl test_mechanism_enumeration.jl test_oracles.jl test_rate_eq_derivation.jl test_reaction_graph.jl test_types.jl test_identify_rate_equation.jl test_fitting.jl`:
  155, 36, 16, 254, 455, 32,582, 667, 3,197, 948, 633, 446, 88 (S3i_delC's twelve files).
- [ ] **Step 6: Commit** (message: 2692bb7's, plus one paragraph for the folded tests and
  one for the docs)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M ROADMAP.md, docs/src/deriving/re_vs_ss.md, src/names.jl, src/thermo.jl, src/types.jl, TE, TRD, TT
git add ROADMAP.md docs/src/deriving/re_vs_ss.md src/names.jl src/thermo.jl src/types.jl test/test_mechanism_enumeration.jl test/test_rate_eq_derivation.jl test/test_types.jl
cat > $P/replay/msg <<'EOF'
Delete the Step views, Step, Species and Residual

No src code or test reads a mechanism through its Step views any more, and
no caller builds one from Step groups, so the transitional types go with
everything that served them: steps(m) and its view builders, the
Step-input constructors and _raw_of, _forms, _re_segment_extras,
_flat_steps, the Step and Species readers and role checks, the Species
content keys in thermo.jl, the Step renderers in names.jl (_form_label,
_step_text, _met_label, _is_reactant_copy) and the two _is_reactant
methods only they called. A mechanism is its reaction graph alone.

With Species gone forms have no internal name, so the testset on E(I) beside
E(Iinh) and E(P::Mod) beside E(Pinh) says what it checks: each pair renders apart
and every constant keeps its own name. Two comments call _testhelper_graph the
test reader.

Docs: re_vs_ss.md says the mechanism stores the distinction as each step's
rapid_equilibrium flag. ROADMAP: the iso entry says a Form carries its
conformation label.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

- [ ] **Step 7: Full suite** (background). F 54: 39,933.
- [ ] **Step 8: Census** (background). C 54 without `CENSUS_ADAPTER`: zero differences,
  order included. `corpus.jl` chose `neutral_graph.jl` by itself: the compare log's `B:`
  line reads `adapter neutral_graph.jl`, as `census/replay_T54/meta.toml` records (the
  corpus log does not print it).

### Task 55: Minors: raw iso test, seed bases, merged complexes [S3i]

**Replays:** `[14ef525 + e69cd27 + b42cc90: test/test_types.jl + 1182f84 + 216be7f:
docs/src/developer.md]` (S3i §9 task 5, test first, D-S3i-9; round 1's m1 and m8, which
the S3i note puts here, are b42cc90's TT piece and 1182f84. That piece also holds m3, a
test-local rename the note puts in task 3: it shares test_types.jl with m1, whose
`step_of(ER._is_iso)` reads raw steps through the `_is_iso(::RawStep)` 14ef525 adds, so
the file lands here, a departure Open for Denis lists. S5 §9 item 3's developer.md.)

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_catalytic_topologies` and `_seed_groups` call
  `_is_iso`; `_seed_candidate_screen`'s message; `_seed_bases`, `_lumping_twin`,
  `_merged_complexes`; `_seed_variants`; `_bound_comp_inhibitors`;
  `_expand_add_dead_end_regulator`; `_chemistry_equilibrates_both_sides`;
  `_fused_substrate_binding` folds into `_merged_complexes`), `src/types.jl`
  (`_is_iso(::RawStep)`), `docs/src/developer.md`
- Test: `test/test_mechanism_enumeration.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: Task 51's seed variants and screen; `_flat_edges`.
- Produces: `_is_iso(s::RawStep)` beside the `Edge` method; `_seed_bases(m) ->
  Vector{RawMechanism}` (the merged base, then the Theorell–Chance bases);
  `_lumping_twin(base, cx, re)`; `_merged_complexes(m, subs)` (the forms a fused substrate
  binding enters, each once), read by `_seed_variants`, `_seed_candidate_screen` and
  `_chemistry_equilibrates_both_sides`; `_bound_comp_inhibitors` reads the edges group by
  group and the dead-end move's `bound_regs` calls it; the screen's message "the seed base
  holds an isomerization; a base has every isomerization merged". TE: `_testhelper_twin`,
  `_testhelper_group_of`; `_testhelper_raw_iso`, `_testhelper_all_re` and
  `_testhelper_seed_bases` go. TT: `iso = step_of(ER._is_iso)`, `copy_binding`.

**Gate:** F 55 = 39,942; C 55: zero differences, order included; SC 55; the R18 probe;
the constructor median; **R52's gate**; checkpoint S3i.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 14ef525 -- test/test_mechanism_enumeration.jl | git apply
git show --format= --binary e69cd27 -- test/test_mechanism_enumeration.jl | git apply
git show --format= --binary b42cc90 -- test/test_types.jl | git apply
```

- [ ] **Step 2: Red** (the screen's message, e69cd27's testset alone; the other new
  tests call `_seed_bases` and `_is_iso(::RawStep)`, which do not exist yet). N
  `$L/T55_red.log` `test/test_mechanism_enumeration.jl`
  `"_seed_candidate_screen refuses a graph that holds an isomerization"`. Expected
  (`$P/proto/logs/S3i_screen_2692bb7.log`): one `Test Failed` with `Expected:
  ErrorException("_seed_candidate_screen: the seed base holds an isomerization; a base has
  every isomerization merged")` and `Thrown: ErrorException("_seed_candidate_screen: a seed
  base holds no isomerization")`, and the summary line `... | FAIL: Some tests did not
  pass: 0 passed, 1 failed, 0 errored, 0 broken.` If the named run stops before that
  testset (a top-level definition of TE that needs the new source), take the red as the
  prototype did, with the same testset in a file of its own: PR `$L/T55_red.log`
  `$P/s3i/probe/screen_message.jl`, the same failure expected.
- [ ] **Step 3: Apply the source and the developer page**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 14ef525 -- src/mechanism_enumeration.jl src/types.jl | git apply
git show --format= --binary 1182f84 -- src/mechanism_enumeration.jl | git apply
git show --format= --binary 216be7f -- docs/src/developer.md | git apply
```

  `git diff --stat -- docs` shows `docs/src/developer.md | 203`.
- [ ] **Step 4: Green** (background). Gt `$L/T55_green.log`
  `test_mechanism_enumeration.jl test_dsl.jl test_types.jl test_identify_rate_equation.jl test_reaction_graph.jl`:
  32,591 (+8 at 14ef525, +1 at e69cd27), 254, 633, 446, 948. The red testset of Step 2
  now passes.
- [ ] **Step 5: Commit** (message: 14ef525's without its first sentence, which names the
  prototype's review rounds; the folds' sentences merged in, since 1182f84 rewrites how
  `_seed_variants` finds the merged complexes)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M docs/src/developer.md, src/mechanism_enumeration.jl, src/types.jl, test/test_mechanism_enumeration.jl, test/test_types.jl
git add docs/src/developer.md src/mechanism_enumeration.jl src/types.jl test/test_mechanism_enumeration.jl test/test_types.jl
cat > $P/replay/msg <<'EOF'
Share the raw isomerization test, the seed bases and the bound-inhibitor read

_is_iso gains a RawStep method, which the topology search, the seed groups and
the seed variants call in place of three inline tests, and the tests'
_testhelper_raw_iso goes. _seed_bases(m) builds the merged and Theorell-Chance
bases that _seed_variants and four seed testsets read, so the tests' copies
(_testhelper_seed_bases, _testhelper_all_re) go. _merged_complexes(m, subs)
gives the forms a fused substrate binding enters: _seed_variants reads them
once, for the merged base, and skips the lumping twin through _lumping_twin,
and the seed screen and the degeneracy check read them in place of their own
copies, so the one-use predicate behind them folds in. _seed_variants'
docstring says a misaligned base errors. The screen's assertion message states
the failure, and a testset pins it; on the old message ("a seed base holds no
isomerization") it fails. _bound_comp_inhibitors reads the edges group by
group, and the dead-end move calls it. In the tests, one _testhelper_twin
serves its three users, and the move testsets that drew fixtures from
init_mechanisms or an earlier move write them inline as macro blocks, with
every expected child written out. An inline copy of the raw isomerization test
reads the helper, and a local no longer shadows Base.copy.

Docs: the developer page describes the graph's fields, raw input,
`_canonicalize`, the representation checks (free E, a connected graph, every
declared name used), the boundary check and `_merged_complexes`, where it
described Step and Species values, the reaction field, the naming cache, the
Step constructor and the content keys of the joint solve.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

- [ ] **Step 6: Full suite** (background). F 55: 39,942. Record its maxrss.
- [ ] **Step 7: Census** (background). C 55: zero differences, order included.
- [ ] **Step 8: R52's gate.** SC 55, then TEA 55, both background. Passed, each figure
  within its gate (S3i measured all three inside it, D-S3i-2's copies being in Task 34):
  - full-suite maxrss from F 55: at most 3.1 GB, `maxrss_kb` at most 3,100,000 (S3i
    2.86-2.94 GB);
  - ter-ter init, SC 55's `tt:` time: at most 31 s (S3i 13.1-13.7 s), with counts
    264 / 239 / 7 / 250,855;
  - TE alone, TEA 55's `TE_ELAPSED`: at most 145 s (S3i 82.9-83.9 s).
  A figure over its gate: stop and report it with the three logs (R52 asks for the cause,
  R57's method); do not go on to S4.
- [ ] **Step 9: R18 probe.** PR `$L/T55_r18.log` `$P/s3i/probe/r18.jl`. Expected
  (`$P/proto/logs/S3i_r18.log`): 18 lines `concrete ...` and one `NOT`,
  `show(IOBuffer, EnzymeMechanism{...}) -> Union{Nothing, Bool}` (the summary `show`'s last
  line, D-S3i-11; Task 59 replaces it).
- [ ] **Step 10: Constructor median.** PR `$L/T55_measure.log` `$P/s3i/scripts/measure.jl`.
  Expected: `specs constructor (µs): (median = <m>, ...)` with `<m>` about 3.5-3.9 (S3i
  3.66 in `$P/proto/logs/S3i_measure_14ef525.log`, S4 3.50, S5 3.82: the spread is the
  machine's) and `bi-bi init constructor (µs)` about 4.2. The gate is TRG's 10 µs median
  (R17), which passed in the suite; record the figures.
- [ ] **Step 11: Checkpoint S3i** (tip 1182f84).
  `git diff --stat 1182f84 HEAD -- src test docs ROADMAP.md Project.toml .claude` lists
  exactly these 19 files: Project.toml 2, ROADMAP.md 52,
  docs/src/deriving/cha_king_altman.md 3, docs/src/deriving/mwc_allostery.md 8,
  docs/src/deriving/ping_pong.md 2, docs/src/deriving/re_vs_ss.md 8,
  docs/src/deriving/regulators.md 6, docs/src/developer.md 209,
  docs/src/getting_started.md 4, docs/src/identify/enumeration_engine.md 13,
  docs/src/identify/tutorial.md 7, the plan (D2), the spec 48 (D1), src/dsl.jl 4, FX 35,
  TE 56, TRD 33, TRG 4, TT 42. src/types.jl leaves the list (050e16f and 3962124 are in the
  tip); TRG's 4 lines are 4952d79's docstring, ahead since Task 53. Then:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
F=$P/replay/roadmap_final.md
$P/replay/ckpt_p4.sh 1182f84 c:2b29fd2 c:be65542 c:0b2021d c:a84d194 \
    p:216be7f:docs/src/deriving/mwc_allostery.md p:216be7f:docs/src/deriving/cha_king_altman.md \
    p:216be7f:docs/src/getting_started.md p:216be7f:src/dsl.jl \
    "e:Encode a mechanism's Sig=delete" "e:Have the \`Mechanism\` constructors check=delete" \
    "e:Take \`:OnlyA\` steps to their limit=$F" "e:\`eq_complexity_filter\` counts only=$F" \
    "e:The \`:OnlyA\` guard misjudges=$F" i p:216be7f:docs/src/identify/enumeration_engine.md \
    p:216be7f:docs/src/deriving/re_vs_ss.md p:216be7f:docs/src/developer.md c:a8f6f51 c:fbdb9aa \
    p:4952d79:test/test_reaction_graph.jl "e:Enumerate Cleland iso=$F"
```

  Expected: no output. Departures in S3i: the ROADMAP iso entry (Task 54); the messages of
  Tasks 53-55; ahead of the tip, S3h's list (but 050e16f and 3962124) plus a8f6f51 (TT),
  fbdb9aa (TE, TRD), 4952d79 (TRG) and 216be7f's re_vs_ss.md and developer.md.

### Task 56: Count a repeated element once [S4]

**Replays:** `[b784fdc]` (S4 §9 task 1; S3e §7's `S[CHC]`, D-S4-13).

**Files:**
- Modify: `src/types.jl` (`ReactantAtoms` and the comment above it)
- Test: `test/test_types.jl`, `test/test_dsl.jl`

**Interfaces:**
- Consumes: `ReactantAtoms`, `@enzyme_reaction`'s atom brackets (Task 15).
- Produces: `ReactantAtoms` sums an element listed twice, so `S[CHC] == S[C2H]`, the two
  reactions hash alike and both print `S[C2H]`.

**Gate:** F 56 = 39,946; C 56: zero differences.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary b784fdc -- test/test_dsl.jl test/test_types.jl | git apply
```

- [ ] **Step 2: Red.** Rt `$L/T56_red.log` `test_types.jl test_dsl.jl`. Expected
  (`$P/proto/logs/S4_c1_red.log`): 4 failures, 2 in "ReactantAtoms canonicalizes atom
  ordering" (TT) and 2 in "multi-atom metabolites" (TD).
- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary b784fdc -- src/types.jl | git apply
```

- [ ] **Step 4: Green.** Gt `$L/T56_green.log` `test_types.jl test_dsl.jl`: 635, 256.
- [ ] **Step 5: Commit.** K `b784fdc` `src/types.jl test/test_dsl.jl test/test_types.jl`
  ("Count an element listed twice in a formula once, with the summed count").
- [ ] **Step 6: Full suite** (background). F 56: 39,946.
- [ ] **Step 7: Census** (background). C 56: zero differences.

### Task 57: @enzyme_reaction resolves its names; the tag hints name a bare metabolite [S4]

**Replays:** `[a520f8e]`, `[3f6c43c + 9578d2a + 761337e + ad3ca91 + eef3af7 + 7c2a657]`
(S4 §9 tasks 2-3; R46 m5 and m3; S5's R59 m4, D-S5-1, the S5 minors' m6, R62's r2, R66's
m1 and R67's m1-m3 fold into the second, test first).

**Files:**
- Modify: `src/dsl.jl` (`@enzyme_reaction`, the `_reaction_expr` docstring; the
  tagged-regulator message in `_raw_mechanism`, `_bare_symbols_from_values`'s hint,
  `_bare_name`, `_untagged`)
- Test: `test/test_dsl.jl`

**Interfaces:**
- Consumes: `_reaction_expr`, `_raw_mechanism`, `_bare_symbols_from_values`.
- Produces: `@enzyme_reaction` returns its expansion unescaped, so it works under
  `import EnzymeRates as ER` and `using EnzymeRates: @enzyme_reaction` (D-S4-9);
  `_bare_name(x)`: the name `x` writes under its `::` tags and a multiplicity call,
  `nothing` when no name stands there (a number, or any other call, such as `A + B`,
  `A in B` or `(-(2))`: a call unwraps only as a multiplicity call, a Symbol head that
  `Base.isoperator` does not flag applied to one or more Integers, as `X(2)` and
  `X(1, 2)`, 7c2a657; that test passes the word operators, so `in(2)` gives `in`);
  `_untagged(x)`: `x` without its `::` tags (9578d2a's helper, which 761337e widened into
  `_bare_name` and ad3ca91 brings back beside it); an atom bracket under a tag
  (`S[C]::Mod`) gets the atom-bracket hint; the tagged-regulator example writes
  `regulatory_site(2): <name>::OnlyA`, and leaves the example out without a name.

**Gate:** F 57 = 39,969 (a520f8e +6, 3f6c43c +1, 9578d2a +2, 761337e +4, ad3ca91 +3,
eef3af7 +5, 7c2a657 +2);
C 57: zero differences.

**Unit 1: a520f8e** (R46 m5)

- [ ] **Step 1: Apply the test path**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary a520f8e -- test/test_dsl.jl | git apply
```

- [ ] **Step 2: Red.** Rt `$L/T57_1_red.log` `test_dsl.jl`. Expected
  (`$P/proto/logs/S4_c2_red.log`): 2 errors in "the macros work in a module that names the
  package otherwise", each `UndefVarError: EnzymeRates not defined` from the
  `@enzyme_reaction` calls in the renamed-import and macro-only modules.
- [ ] **Step 3: Apply the source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary a520f8e -- src/dsl.jl | git apply
```

- [ ] **Step 4: Green.** Gt `$L/T57_1_green.log` `test_dsl.jl`: 262.
- [ ] **Step 5: Commit.** K `a520f8e` `src/dsl.jl test/test_dsl.jl` ("Let @enzyme_reaction
  resolve the package's names itself").

**Unit 2: 3f6c43c + 9578d2a + 761337e + ad3ca91 + eef3af7 + 7c2a657** (R46 m3, R59 m4
with D-S5-1, S5 minors m6, R62 r2, R66 m1, R67 m1-m3)

- [ ] **Step 6: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 3f6c43c -- test/test_dsl.jl | git apply
git show --format= --binary 9578d2a -- test/test_dsl.jl | git apply
git show --format= --binary 761337e -- test/test_dsl.jl | git apply
git show --format= --binary ad3ca91 -- test/test_dsl.jl | git apply
git show --format= --binary eef3af7 -- test/test_dsl.jl | git apply
git show --format= --binary 7c2a657 -- test/test_dsl.jl | git apply
```

- [ ] **Step 7: Red.** Rt `$L/T57_2_red.log` `test_dsl.jl`. Expected: 262 pass and 17
  failures, the six commits' reds together (`$P/proto/logs/S4_c3_red.log`,
  `S5_red_m3m4.log`, `S5_red_m4b.log`, `S5m_red_m6.log`, `S5m_red_m6_site2.log`,
  `S5m2_red_r2.log`, `S5m3_red_m1.log`, `S5m3_mutant_shape_only.log`, `S5m4_red_m3.log`;
  the plan's update ran this red on the dry replay's Task 57 tree,
  `$P/plan/evidence/r67_update/T57_2_red.log`):
  - "every allosteric regulator sits at exactly one regulatory site" 8: the example
    echoes `X::Activator`, `X::Activator::OnlyA`, `X(2)::OnlyA`, `X(1, 2)::OnlyA`,
    `1::OnlyA`, `(A + B)::OnlyA`, `(A in B)::OnlyA` and `(A isa B)::OnlyA` into
    `regulatory_site(2): ...` (3f6c43c, 9578d2a, 761337e twice, ad3ca91, eef3af7 twice,
    7c2a657);
  - "a mechanism macro's list labels take bare names" 9: the hint advises listing
    `P::Mod`, `P(2)`, `X(1, 2)`, `1`, `A + B`, `S[C]`, `A in B`, `A isa B` and `-2` bare
    (9578d2a, 761337e twice, ad3ca91 twice, eef3af7 three times, 7c2a657).
- [ ] **Step 8: Apply the source, in commit order**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 3f6c43c -- src/dsl.jl | git apply
git show --format= --binary 9578d2a -- src/dsl.jl | git apply
git show --format= --binary 761337e -- src/dsl.jl | git apply
git show --format= --binary ad3ca91 -- src/dsl.jl | git apply
git show --format= --binary eef3af7 -- src/dsl.jl | git apply
git show --format= --binary 7c2a657 -- src/dsl.jl | git apply
grep -c '_untagged' src/dsl.jl; grep -c '_bare_name' src/dsl.jl
```

  Expected: `3` (ad3ca91 brings `_untagged` back: its definition and its calls in
  `_bare_symbols_from_values` and `_bare_name`), then `3`.
- [ ] **Step 9: Green.** Gt `$L/T57_2_green.log` `test_dsl.jl`: 279 (the plan's update
  ran it on the dry replay's Task 57 tree, `$P/plan/evidence/r67_update/T57_2_green.log`).
- [ ] **Step 10: Commit** (message: 3f6c43c's, plus one paragraph for the five folds)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M src/dsl.jl, M test/test_dsl.jl
git add src/dsl.jl test/test_dsl.jl
cat > $P/replay/msg <<'EOF'
Show a site tag, not the user's, in the tagged-regulator error's example

A tagged `allosteric_regulators:` entry is rejected with an example of the
site line its tag belongs on. The example repeated the tag the user wrote,
so `X::Activator`, a regulator type of @enzyme_reaction and no site tag,
came back as `regulatory_site(2): X::Activator`. The example now writes the
ligand with `::OnlyA`.

The example and the list labels' role-tag hint name the bare metabolite
(`_bare_name`), past nested tags and a multiplicity call, a Symbol head that
`Base.isoperator` does not flag applied to one or more Integers:
`allosteric_regulators: X::A::B`, `X(2)::OnlyA` and `X(1, 2)::OnlyA` give
`regulatory_site(2): X::OnlyA`, and `modifiers: P::Mod::Mod`, `P(2)::Mod` and
`P(1, 2)::Mod` advise listing `P` bare. The head test passes the word
operators, so `in(2)` counts as a call. A value with no name there gets
neither: a number (`1::OnlyA`, `1::Mod`) or any other call, such as an
operator call (`(A + B)::Mod`, `(-(2))::Mod`) or the infix `(A in B)::Mod` and
`(A isa B)::Mod`, whose arguments are no Integers. An atom bracket under a tag
(`S[C]::Mod`) gets the atom-bracket hint.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

**Task end**

- [ ] **Step 11: Full suite** (background). F 57: 39,969.
- [ ] **Step 12: Census** (background). C 57: zero differences.

### Task 58: Reject names Julia reads as another name [S4]

**Replays:** `[13f202c + 02214d5]` (S4 §9 task 4, D-S4-5; R59 m1 is 02214d5).

**Files:**
- Modify: `Project.toml` (`Unicode` in `[deps]`, compat `Unicode = "1"`),
  `src/EnzymeRates.jl` (`import Unicode`), `src/types.jl` (`_check_metabolite_names`, the
  `EnzymeReaction` docstring)
- Test: `test/test_types.jl`, `test/test_reaction_graph.jl`

**Interfaces:**
- Consumes: `_check_metabolite_names(who, names)` (Tasks 8, 38).
- Produces: `_check_metabolite_names` rejects a non-ASCII name not in Julia's identifier
  normal form (NFC plus `Unicode.julia_chartransform`): "Julia reads the name "µM" as
  "μM", so no DSL text can write it; use "μM"". The Unicode stdlib joins the
  dependencies.

**Gate:** F 58 = 39,979; C 58: zero differences.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 13f202c -- test/test_reaction_graph.jl test/test_types.jl | git apply
git show --format= --binary 02214d5 -- test/test_reaction_graph.jl test/test_types.jl | git apply
grep -c 'u00b5M' test/test_types.jl test/test_reaction_graph.jl
```

  Expected: each file holds the escaped micro sign at least once (02214d5).
- [ ] **Step 2: Red.** Rt `$L/T58_red.log` `test_types.jl test_reaction_graph.jl`.
  Expected: 5 failures (`$P/proto/logs/S4_c4_red.log` for TT, `S4_c4_red_trg.log` for
  TRG): 2 in "EnzymeReaction validation" (µM and a decomposed é are accepted) and 3 in
  "the constructors reject what no graph represents, each with its message" (µM, a
  decomposed é and the ligand ɛ are accepted).
- [ ] **Step 3: Apply the source and resolve**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 13f202c -- Project.toml src/EnzymeRates.jl src/types.jl | git apply
$P/jl.sh $L/T58_resolve.log 900 --project=. -e 'using Pkg; Pkg.resolve()'
cd $P && ./jl.sh $L/T58_resolve_census.log 900 --project=census/env_replay -e 'using Pkg; Pkg.resolve()'
```

  Expected: both `exit 0`; `grep -A1 '^\[\[deps.EnzymeRates\]\]'
  $P/census/env_replay/Manifest.toml` shows `"Unicode"` in its `deps` line.
- [ ] **Step 4: Green.** Gt `$L/T58_green.log` `test_types.jl test_reaction_graph.jl`: 642,
  951 (`$P/proto/logs/S4_c4_green.log`).
- [ ] **Step 5: Commit** (message: 13f202c's, plus one paragraph for 02214d5)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M Project.toml, src/EnzymeRates.jl, src/types.jl, test/test_reaction_graph.jl, test/test_types.jl
git add Project.toml src/EnzymeRates.jl src/types.jl test/test_reaction_graph.jl test/test_types.jl
cat > $P/replay/msg <<'EOF'
Reject a metabolite name that Julia reads as another name

Julia's parser normalizes an identifier: NFC, and `julia_chartransform`,
which reads a micro sign as μ and a Latin open e as ε. A name built
directly with `Symbol("µM")` or a decomposed é passed `isidentifier`, so a
reaction, mechanism or ligand could hold it, and the text such a value
prints would read back as another name. The name check now rejects a name
not in that normal form and says which name Julia reads. It normalizes
non-ASCII names only, with the Unicode stdlib, which documents the
equivalence with the parser; an ASCII name costs nothing more.

The rejection tests write those names as escapes: the micro sign, the
decomposed é and the open epsilon look like the names Julia reads them as, so
written out each pin would read as two equal strings, and an editor that
normalizes the file would turn a rejection case into an accepted name.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

- [ ] **Step 6: Full suite** (background). F 58: 39,979; Aqua's stale-dependency check
  passes with `Unicode` used.
- [ ] **Step 7: Census** (background). C 58: zero differences.

### Task 59: Print a mechanism as its macro text [S4]

**Replays:** `[3464d63 + b8a604f + 4952d79: test/test_types.jl + 919bc4c: test/test_types.jl
+ cd7a4ae: src/names.jl FX TE TO]` (S4 §9 task 5, which takes b8a604f's
`@nospecializeinfer` and TCB window test first, and S5 §9 item 5's R59 minors for this
task: m5 is 4952d79's TT piece (its TRG piece, R58 m1, is in Task 53); m7's TT closure is
919bc4c's TT piece (its TID piece is Task 61's); m6's predicate and test-side edge text are
cd7a4ae's src/names.jl, FX, TE and TO pieces. cd7a4ae's TD and TRG pieces edit tests that
0ec7d51 writes, so they land with it in Task 60).

**Files:**
- Modify: `src/names.jl` (`_met_text`, `_form_text`, `_step_text` and
  `_is_reactant_copy(m, i)` written; `_side_text` and `_edge_text` go), `src/dsl.jl`
  (section "Mechanisms as text": the four `show` methods,
  `_print_line`, `_print_names`, `_print_steps`), `src/types.jl` (the "Pretty printing"
  section goes)
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
  (`_testhelper_edge_text`), `test/test_compile_budget.jl`,
  `test/test_mechanism_enumeration.jl` (two lookups), `test/test_oracles.jl`,
  `test/test_reaction_graph.jl`, `test/test_types.jl`

**Interfaces:**
- Consumes: the graph's table, forms and edges; `_form_name`, `_entry_label`,
  `_flat_edges`, `_group_members`, `modifiers(m)`, `allosteric_regulators(am)`.
- Produces: `show(io, m::Mechanism)` and `show(io, am::AllostericMechanism)` print the
  macro block (S4 §0.3: sorted names, `modifiers:` only when a step binds one, one line per
  kinetic group in stored order and direction, `:: Tag` on an allosteric group,
  `catalytic_multiplicity:` always, the sites); `Base.@nospecializeinfer show(io,
  @nospecialize(em::EnzymeMechanism))` and the allosteric one print the concrete
  mechanism's text; `_met_text(m, i)`, `_form_text(m, f)` (a form with nothing bound prints
  its bare label, `E`), `_step_text(m, k)`, and `_is_reactant_copy(m, i)`, which
  `_entry_label` and `_met_text` read (cd7a4ae). FX's `_testhelper_edge_text` renders the
  step text the tests pin (D-S4-1), each metabolite by its label (`Amod` for a copy,
  cd7a4ae), so two TE lookups read `E + Amod ⇌ EAmod`; TT's printing tests read back
  through `_testhelper_builds_back` (919bc4c). TCB counts the show instances PASS1 and
  printing a compiled mechanism add (gate 0); its `@info` line reads "show instances per
  fresh mechanism".

**Gate:** F 59 = 39,979 (TT's pins rewritten, net 0); TCB "fit compile: PASS1 and a compiled
mechanism's text compile no show" 0 instances; C 59: zero differences (the census reads no
`show`, S4 §5).

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 3464d63 -- test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_oracles.jl test/test_reaction_graph.jl test/test_types.jl | git apply
git show --format= --binary b8a604f -- test/test_compile_budget.jl | git apply
git show --format= --binary 4952d79 -- test/test_types.jl | git apply
git show --format= --binary 919bc4c -- test/test_types.jl | git apply
git show --format= --binary cd7a4ae -- test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_mechanism_enumeration.jl test/test_oracles.jl | git apply
```

- [ ] **Step 2: Red** (background). Rt `$L/T59_red.log`
  `test_compile_budget.jl test_types.jl test_reaction_graph.jl test_oracles.jl test_joint_solve.jl`.
  Expected: 26 failures. `$P/proto/logs/S4_c5_red.log`'s 25: TT 23 ("a mechanism prints as
  the text of the macro that builds it" 20, "an allosteric kinetic group prints as one
  tagged tuple" 1, "an allosteric mechanism prints every tag and every site" 1, "two
  distinct forms that render one name are rejected" 1, its `Estar()` message), TRG 2 ("a
  form's name and text read as constants and the macros write them" 1, "the constructors
  reject what no graph represents, each with its message" 1), TO 667 and TJS 455 passing;
  and TCB 1 in "fit compile: PASS1 and a compiled mechanism's text compile no show" (the
  summary `show` compiles per Sig; its `@info` line reads a nonzero count).
- [ ] **Step 3: Apply 3464d63's source alone**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 3464d63 -- src/dsl.jl src/names.jl src/types.jl | git apply
git show --format= --binary cd7a4ae -- src/names.jl | git apply
```

- [ ] **Step 4: Red for b8a604f** (background). Rt `$L/T59_red_show.log`
  `test_compile_budget.jl`. Expected (`$P/proto/logs/S4_c7a_red.log`): 15 pass, 1 failure
  in "fit compile: PASS1 and a compiled mechanism's text compile no show", the log's
  `[ Info: show instances per fresh mechanism: [1.0, 1.0, 1.0, 1.0] (budget: 0)`:
  `@nospecialize` alone lets `sprint(show, em)` infer one show per fresh Sig (D-S4-6).
- [ ] **Step 5: Apply b8a604f's source**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary b8a604f -- src/dsl.jl | git apply
```

- [ ] **Step 6: Green** (background). Gt `$L/T59_green.log`
  `test_compile_budget.jl test_types.jl test_reaction_graph.jl test_oracles.jl test_joint_solve.jl test_mechanism_enumeration.jl`:
  16, 642, 951, 667, 455, 32,591 (TE's two lookups read the copy's label, `Amod`); the log
  holds
  `[ Info: show instances per fresh mechanism: [0.0, 0.0, 0.0, 0.0] (budget: 0)`.
- [ ] **Step 7: Commit** (message: 3464d63's, with b8a604f's reason, 4952d79's comments
  and the cd7a4ae and 919bc4c pieces)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M src/dsl.jl, src/names.jl, src/types.jl, FX, TCB, TE, TO, TRG, TT
git add src/dsl.jl src/names.jl src/types.jl test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_compile_budget.jl test/test_mechanism_enumeration.jl test/test_oracles.jl test/test_reaction_graph.jl test/test_types.jl
cat > $P/replay/msg <<'EOF'
Print a mechanism as the text of the macro that builds it

A compiled mechanism printed a summary in internal form names (`ES <--> EP`,
`reg site 1 (n=3): ...`) that no macro reads, and a concrete one printed
Julia's default struct text. Both now print their @enzyme_mechanism or
@allosteric_mechanism block: the sorted names, a modifiers: line for the
names the steps bind as modifiers, every kinetic group on one line in
stored order and direction, the tags, the catalytic multiplicity and the
sites. The text depends on the mechanism alone, so identify can hash it
into the mechanism's id. The compiled types print through the concrete
mechanism and take their argument @nospecialize and @nospecializeinfer, so a
fresh Sig compiles no show method; with @nospecialize alone, a caller compiled
for one Sig, such as `sprint`, inferred a show for that Sig. The compile-budget
census prints each fresh compiled mechanism after the fit census and counts the
show instances that adds: 0 per fresh mechanism, against 1 without
@nospecializeinfer.

A form with nothing bound prints as its bare label (`E`), as the DSL
writes free E; the form-name collision message follows (`Estar`, not
`Estar()`). `_edge_text` and `_side_text`, whose one caller was the
summary, become the test helper that renders the internal-name step text
the tests pin. Each new pin also builds the mechanism back through the
macro, and the printing tests say why two texts read as they do: the 3-cycle
prints in stored order, not chain order, and a fused binding prints B on its
entry side because E(P, Q) does not name it. `_entry_label` and `_met_text`
share one reactant-copy predicate, `_is_reactant_copy(m, i)`. The test-side
edge text writes labels, so two TE lookups name a modifier copy `Amod`, and
the two-line read-back closure of the printing tests takes the `_testhelper_`
prefix as `_testhelper_builds_back`, since the macro, not the parser, reads
the text back.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

- [ ] **Step 8: Full suite** (background). F 59: 39,979, pristine: the 12 `[ Info:` lines,
  the compile gate's now "show instances per fresh mechanism".
- [ ] **Step 9: Census** (background). C 59: zero differences.

### Task 60: Read a mechanism back without eval [S4]

**Replays:** `[0ec7d51 + aa158dd + 8aed057 + cd7a4ae: test/test_dsl.jl
test/test_reaction_graph.jl + 26f7323 + 32ea3c9]` plus 216be7f's ROADMAP entry below (S4
§9 task 6; R59 m3 and m2, m6's TD and TRG pieces (they edit 0ec7d51's tests; its other
pieces are Task 59's), S5 round 1's m2 and round 2's B2 fold in, aa158dd test first).

**Files:**
- Modify: `src/dsl.jl` (`_checked_mechanism`, `_compiled_mechanism`, `_parse_mechanism`; the
  `show` and parser docstrings), `src/admission.jl` (`_check_well_posed`'s docstring, the
  ABOUTME), `docs/src/developer.md` (the boundary-check sentence), `ROADMAP.md`
- Test: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
  (`_testhelper_reads_back`), `test/test_dsl.jl`, `test/test_reaction_graph.jl`

**Interfaces:**
- Consumes: `_raw_mechanism(block, allosteric)` (Task 6), `_check_uniquely_decodable`,
  the checked constructors, `_check_well_posed` (Task 41), Task 59's `show`.
- Produces: `_checked_mechanism(raw::RawMechanism) -> Mechanism` and its allosteric
  method (decodability, constructor, boundary check); the macros expand to
  `_compiled_mechanism(raw) = compile_mechanism(_checked_mechanism(raw))`;
  `_parse_mechanism(text::AbstractString) -> Mechanism or AllostericMechanism` (one
  unqualified macro call, surrounding whitespace stripped, never evaluated; the macros'
  messages). FX: `_testhelper_reads_back(x)`. TRG:
  `_testhelper_census_reactions`, `_testhelper_census_allosteric`, and "a graph's DSL text
  reads back as the graph" (42 specs, 11,013 plain and 944 allosteric mechanisms).

**Gate:** F 60 = 40,093 (0ec7d51 +113, aa158dd +1); C 60: zero differences.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 0ec7d51 -- test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_dsl.jl test/test_reaction_graph.jl | git apply
git show --format= --binary aa158dd -- test/test_dsl.jl | git apply
git show --format= --binary 8aed057 -- test/test_dsl.jl | git apply
git show --format= --binary cd7a4ae -- test/test_dsl.jl test/test_reaction_graph.jl | git apply
```

- [ ] **Step 2: Red.** Rt `$L/T60_red.log` `test_dsl.jl test_reaction_graph.jl`. Expected:
  every failure and error in the parser's testsets, each error an `UndefVarError` for
  `_parse_mechanism` (or `_checked_mechanism`): `$P/proto/logs/S4_c6_red.log` had 9
  failures and 88 errors on 0ec7d51's tests; the folds make it 9 failures and 102 errors:
  - TD "a text the macros reject is rejected with their message" 4 failures;
  - TD "the text is one mechanism block, read without evaluation" 5 failures and 2 errors
    (aa158dd's trailing-blank-lines case adds one, `$P/proto/logs/S5_red_m3m4.log`);
  - TD "every mechanism the macros accept reads back" (8aed057's title) 14 errors: cd7a4ae
    moves the read-back inside `@test`, so each of the seven mechanisms errors twice where
    0ec7d51's loop stopped at its first;
  - TRG 86 errors in "a graph's DSL text reads back as the graph": two in each of its 42
    per-spec testsets ("Uni-Uni", "HK", ..., "m_all") and two at its own level.
  Any other failure or error: stop.
- [ ] **Step 3: Apply the source and the developer page, in commit order**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 0ec7d51 -- src/dsl.jl | git apply
git show --format= --binary aa158dd -- src/dsl.jl | git apply
git show --format= --binary 8aed057 -- src/dsl.jl | git apply
git show --format= --binary 26f7323 -- src/admission.jl | git apply
git show --format= --binary 32ea3c9 -- docs/src/developer.md src/admission.jl | git apply
```

- [ ] **Step 4: Departure: the ROADMAP entry** (216be7f's piece; it names the parser, so it
  lands with it). Set "`BigInt` arithmetic inflates the GC's live-byte count." to its final
  text:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
python3 $P/replay/roadmap_entry.py '`BigInt` arithmetic inflates' $P/replay/roadmap_final.md ROADMAP.md
git diff --stat -- ROADMAP.md
```

  Expected: `ROADMAP.md | 40 +++++++++++++++++++++++-----------------` (23 lines in, 17
  out), and ENT `` `BigInt` arithmetic inflates `` prints nothing. The entry then
  reads (a blank line follows it, as before):

```markdown
- **`BigInt` arithmetic inflates the GC's live-byte count.**
  Three functions that compute in `BigInt` raise `Base.gc_live_bytes` by bytes that no
  live object holds: the `:OnlyA` guard's test `_onlya_offenders` (`src/thermo.jl`; a
  `Rational{BigInt}` null space) by about 15 KiB per call, the counter that
  `_split_gain_test` uses for a plain mechanism (`_partition_independent_count`; a
  `Rational{BigInt}` row reduction) by 18 KiB and `_eq_complexity` (a `BigInt`
  determinant, skipped for a single rapid-equilibrium segment) by 15 KiB. A 6×6
  determinant outside the package does the same in `BigInt` and `Rational{BigInt}`, not in
  `Rational{Int}`. The boundary check (`_check_well_posed`, `src/admission.jl`) runs the
  guard on every allosteric mechanism the macros or the parser build, at the enumerator's
  exits and at every derivation entry, and the allosteric `_split_gain_test` runs its test
  on each candidate. 100,000 runs of the guard add 1.46 GiB and 12 full sweeps. The joint
  solve adds nothing, on a plain mechanism or on a concrete allosteric one; it runs in
  `Rational{Int}` on 10,226 enumerated bi-bi and ter-bi mechanisms without overflow. The
  count stays up after a full collection: 100,000 guard runs that start 1.47 GiB high make
  52 full sweeps in 10.0 s, against 12 in 6.5 s from a clean start. The count holds no
  memory and raises no peak: a process whose count reaches 2.8 GiB peaks at 1.2 GB RSS,
  and the enumeration tests' count after a full collection stays near 1.8 GiB whether
  their peak RSS is 2.2 or 3.5 GB. The higher peak came from mechanisms that kept vectors
  at the capacity their builds grew to; `_canonicalize` stores each at its length, and
  `test/test_reaction_graph.jl` pins it. Admission under Planned replaces the guard, not
  the other two. A loop that calls one of these and prints `Base.gc_live_bytes()` and
  `Base.gc_num().full_sweep` shows it. The figures come from a random-order bi-bi of nine
  steady-state steps, plain and as a two-subunit allosteric mechanism with `:OnlyA` on
  both bindings of A and on the chemical step and `:EqualAI` elsewhere (Julia 1.12.5,
  aarch64, `--heap-size-hint=2500M`).
```

- [ ] **Step 5: Green** (background). Gt `$L/T60_green.log`
  `test_dsl.jl test_reaction_graph.jl test_mechanism_enumeration.jl test_oracles.jl test_types.jl test_joint_solve.jl test_fitting_scale.jl`:
  306, 1,038, 32,591, 667, 642, 455, 104. TRG's canonical golden passes unchanged
  (cd7a4ae writes its text with `_testhelper_edge_text`).
- [ ] **Step 6: Commit** (message: 0ec7d51's, with a paragraph for the folds and one for
  the entry)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M ROADMAP.md, docs/src/developer.md, src/admission.jl, src/dsl.jl, FX, TD, TRG
git add ROADMAP.md docs/src/developer.md src/admission.jl src/dsl.jl test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_dsl.jl test/test_reaction_graph.jl
cat > $P/replay/msg <<'EOF'
Read a mechanism back from its text without evaluating it

`_parse_mechanism(text)` takes one @enzyme_mechanism or @allosteric_mechanism
block, such as `show` prints, and returns the concrete mechanism. It parses
the text, walks the expression with the macros' walker and runs the macros'
checks through `_checked_mechanism`, which the macros now call too, so the
parser accepts exactly what the macros accept, rejects the rest with their
messages, and evaluates nothing (an interpolation stays a walker error). It
strips the whitespace around the text first: a block copied from a CSV cell or
a file often ends in blank lines, which Meta.parse reads as an extra token. The
show and _parse_mechanism docstrings claim the read-back for the mechanisms the
macros accept; the constructors also accept an ill-posed mechanism, or one
whose names spell a string two ways, and the parser rejects its text.

The text reads back as the mechanism that printed it for the 42 specs, the
census's first populations (11,013 plain mechanisms) and its first 944
allosteric ones, and for contextual keywords, Greek names, labels at the
edge of the shape rule, residuals that only subtract, repeated modifier
copies, dual roles, a mechanism without sites and two-digit multiplicities.

_check_well_posed's docstring and admission.jl's header name the parser among
the boundary check's callers. The docstring and the Developer page state the
one-flip-segment and :OnlyA conditions for allosteric mechanisms only: a plain
iso mechanism, whose free forms E and Eprime meet at a steady-state step, is
valid. In the tests, TRG's canonical text writes each edge with
_testhelper_edge_text, as TE's lookups do, and TD's read-back loop calls
_testhelper_reads_back.

ROADMAP: the BigInt entry names the guard's test, _onlya_offenders, and the
boundary check, which runs it on every allosteric mechanism the macros or the
parser build; it carries figures re-measured on a guard loop and records that
the count raises no memory peak, the higher peak having come from vectors kept
at their grown capacity.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

- [ ] **Step 7: Full suite** (background). F 60: 40,093.
- [ ] **Step 8: Census** (background). C 60: zero differences.

### Task 61: identify's CSVs store the text and its id [S4]

**Replays:** `[e59fb5b + 919bc4c: test/test_identify_rate_equation.jl + 81f9f1e + 5cfd1e6 +
1b228ef + 09aec13]` (S4 §9 task 7; R59 m7's and m8's TID edits are 919bc4c's TID piece, its
TT piece is Task 59's; 81f9f1e lands whole here, where every sentence is true; S5
round 1's m3, round 2's m3 and the S5 minors' m4).

**Files:**
- Modify: `Project.toml` (`SHA` in `[deps]`, compat `SHA = "0.7, 1"`),
  `src/identify_rate_equation.jl` (`_mechanism_id`, `_write_rows_csv`,
  `_rows_to_dataframe`, `_failure_row`, `_process_batch`, `_beam_search`;
  `_mechanism_type_string` and `_print_type_value` go), `src/dsl.jl` (the macro
  docstrings), `ROADMAP.md`, `docs/src/developer.md`, `docs/src/getting_started.md`,
  `docs/src/identify/model_selection.md`, `docs/src/identify/tutorial.md`
- Test: `test/test_identify_rate_equation.jl`, `test/test_compile_budget.jl`

**Interfaces:**
- Consumes: Task 59's `show` (the text), Task 60's `_parse_mechanism` (the read-back).
- Produces: `_mechanism_id(text) = bytes2hex(sha256(text))[1:16]`;
  `_write_rows_csv(save_dir, filename, rows, texts)` errors before writing when a row's id
  names a second text; CSV columns `mechanism`, `mechanism_id`, `parent_id` in place of
  `mechanism_type` and `parent_mechanism_type`, every other column kept (R16, D-S4-7).
  TID: `_testhelper_read_csv` (919bc4c).

**Gate:** F 61 = 40,109; C 61: zero differences; R52's figures recorded; the R18 probe;
the docs build; checkpoint S4.

- [ ] **Step 1: Apply the test paths**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary e59fb5b -- test/test_compile_budget.jl test/test_identify_rate_equation.jl | git apply
git show --format= --binary 919bc4c -- test/test_identify_rate_equation.jl | git apply
git show --format= --binary 09aec13 -- test/test_identify_rate_equation.jl | git apply
```

- [ ] **Step 2: Red** (background). Rt `$L/T61_red.log`
  `test_identify_rate_equation.jl test_compile_budget.jl test_types.jl`. Expected
  (`$P/proto/logs/S4_c7_red.log`, the CSV testset under 919bc4c's title): TID 10 errors,
  `_mechanism_id` undefined ("a failure row holds its mechanism's text and id" 4, and 1
  each in "_rows_to_dataframe", "results structure", "CSV output: every file's columns,
  texts and ids", "_mechanism_id cuts the SHA-256 of a mechanism's text to 16 hex digits",
  "csv writer" and "_cv_model_selection: flatten reproduces serial LOOCV per candidate");
  TCB 8 failures, its traced fit census erroring on `_mechanism_id` ("fit compile: PASS1 +
  fit + LOOCV fold per fresh mechanism" 5, "Sig encoding: total compile for the bi-bi
  seeds" 2, "fit compile: PASS1 and a compiled mechanism's text compile no show" 1, which
  the log, taken before b8a604f renamed it, calls "fit compile: PASS1 compiles no show for
  a fresh Sig shape"); TT passes, 642.
- [ ] **Step 3: Apply the source, the docs and the ROADMAP edits; resolve**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary e59fb5b -- Project.toml src/identify_rate_equation.jl | git apply
git show --format= --binary 81f9f1e -- ROADMAP.md docs/src/developer.md docs/src/getting_started.md docs/src/identify/model_selection.md docs/src/identify/tutorial.md src/dsl.jl | git apply
git show --format= --binary 5cfd1e6 -- docs/src/identify/tutorial.md | git apply
git show --format= --binary 1b228ef -- docs/src/identify/model_selection.md | git apply
$P/jl.sh $L/T61_resolve.log 900 --project=. -e 'using Pkg; Pkg.resolve()'
cd $P && ./jl.sh $L/T61_resolve_census.log 900 --project=census/env_replay -e 'using Pkg; Pkg.resolve()'
$P/jl.sh $L/T61_resolve_docenv.log 900 --project=$P/replay/docenv -e 'using Pkg; Pkg.resolve()'
```

  Expected: the three resolves `exit 0`. Every Manifest already lists `SHA` and `Unicode`
  (other packages use them); what changes is the `deps` line under `[[deps.EnzymeRates]]`
  in `$P/census/env_replay/Manifest.toml` and `$P/replay/docenv/Manifest.toml`, which
  `grep -A1 '^\[\[deps.EnzymeRates\]\]' <Manifest>` shows holding `"SHA"` (and, in the
  docs environment, which Task 58 did not resolve, `"Unicode"`). Then, in a new call,
  `git diff --stat -- ROADMAP.md` shows 13 lines (81f9f1e:
  the graph-rewrite entry lists the two PRs that remain; "Decided against: loading old
  CSVs" covers CSVs that name a mechanism by its type string).
- [ ] **Step 4: Green** (background). Gt `$L/T61_green.log`
  `test_identify_rate_equation.jl test_compile_budget.jl test_types.jl`: 462, 16, 642.
- [ ] **Step 5: Commit** (message: e59fb5b's, with a paragraph for the test names and one
  each for the docs and ROADMAP)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M Project.toml, ROADMAP.md, docs/src/developer.md, docs/src/getting_started.md, docs/src/identify/model_selection.md, docs/src/identify/tutorial.md, src/dsl.jl, src/identify_rate_equation.jl, TCB, TID
git add Project.toml ROADMAP.md docs/src/developer.md docs/src/getting_started.md docs/src/identify/model_selection.md docs/src/identify/tutorial.md src/dsl.jl src/identify_rate_equation.jl test/test_compile_budget.jl test/test_identify_rate_equation.jl
cat > $P/replay/msg <<'EOF'
Store each mechanism in identify's CSVs as its text and a stable id

The `mechanism_type` and `parent_mechanism_type` columns held the compiled
type's text, which only `Core.eval` read back. The CSVs now hold
`mechanism`, the mechanism's DSL text, which `_parse_mechanism` reads back
without evaluation (median 342 characters over the bi-bi seeds' children,
against 516 for the type's text); `mechanism_id`, the first 16 hex digits
of the SHA-256 of that text (the SHA stdlib joins the dependencies), the
same in every process and Julia version; and `parent_id`, the parent's id.
Every other column stays. PASS1 renders the text from the concrete
mechanism, so it compiles nothing per Sig, and a failure row no longer
compiles to name its mechanism. A run keeps every id it wrote with its
text, and an id that names a second text is an error before the file is
written. `_mechanism_type_string` and `_print_type_value` had no other
reader and are deleted.

The CSV tests' two-line closure takes the `_testhelper_` prefix
(`_testhelper_read_csv`), and the CSV testset is named by what it checks: every
file's columns, texts and ids.

Docs: a mechanism prints as the macro block that builds it, and identify's CSVs
name it by that text (`mechanism`), its SHA-256 id (`mechanism_id`) and its
parent's id (`parent_id`). The model-selection table lists the three columns,
names `EnzymeRates._parse_mechanism` as the text's reader and says to read the
ids with `types = String`, since CSV.jl reads an id that looks like a number
(decimal digits, or digits around one `e`) as a number, as the CSV test's
comment says too; the macro docstrings, Getting Started and the identify
tutorial say how a mechanism displays; the developer page says `show` compiles
once and the parser never evaluates. The tutorial shows the first five rows of
`cv_results` in five of its columns, so its table stays under Documenter's size
threshold.

ROADMAP: the graph-rewrite entry lists the two PRs that remain, and "Decided
against: loading old CSVs" also covers CSVs that name a mechanism by its type
string.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
```

- [ ] **Step 6: Full suite** (background). F 61: 40,109.
- [ ] **Step 7: Census** (background). C 61: zero differences.
- [ ] **Step 8: R52's figures** (record). SC 61, then TEA 61, both background. S4 measured
  2.84 GB, 13.3 s and 83.7 s (all within R52).
- [ ] **Step 9: R18 probe.** PR `$L/T61_r18.log` `$P/s4/probe/r18.jl`. Expected
  (`$P/proto/logs/S4_r18_e59fb5b.log`): 16 lines `concrete ...` and three `NOT`:
  `_compiled_mechanism(EnzymeRates.RawMechanism,) -> EnzymeMechanism`, the allosteric one
  `-> AllostericEnzymeMechanism` (the Sig is a run-time value) and
  `_parse_mechanism(String,) -> Union{EnzymeRates.AllostericMechanism,
  EnzymeRates.Mechanism}` (two non-parametric structs; no package loop calls it). Every
  `show` returns `Nothing`.
- [ ] **Step 10: Docs build** (background)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
cd $W && $P/jl.sh $L/T61_docs.log 1800 --project=$P/replay/docenv --heap-size-hint=2500M docs/make.jl; rm -rf $W/docs/build
```

  Expected (`$P/proto/logs/S5_docs_m_7395164.log`): `exit 0`; doctests and
  `checkdocs = :exports` pass; one warning, `Documenter could not auto-detect the building
  environment. Skipping deployment.` The tutorial's size warning is gone (D-S5-3).
- [ ] **Step 11: Checkpoint S4** (tip e59fb5b).
  `git diff --stat e59fb5b HEAD -- src test docs ROADMAP.md Project.toml .claude` lists
  exactly these 25 files: Project.toml 2, ROADMAP.md 105,
  docs/src/deriving/cha_king_altman.md 3, docs/src/deriving/mwc_allostery.md 8,
  docs/src/deriving/ping_pong.md 2, docs/src/deriving/re_vs_ss.md 8,
  docs/src/deriving/regulators.md 6, docs/src/developer.md 215,
  docs/src/getting_started.md 5, docs/src/identify/enumeration_engine.md 13,
  docs/src/identify/model_selection.md 5, docs/src/identify/tutorial.md 19, the plan
  (D2), the spec 48 (D1), src/admission.jl 22, src/dsl.jl 60, src/names.jl 11, FX 50,
  TD 143, TID 16, TE 60, TO 4, TRD 33, TRG 32, TT 78. Then, every S5 commit folded into S4
  applied to the tip in prototype order:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
F=$P/replay/roadmap_final.md
$P/replay/ckpt_p4.sh e59fb5b c:aa158dd c:9578d2a c:8aed057 c:02214d5 c:cd7a4ae c:4952d79 \
    c:919bc4c c:a8f6f51 c:2b29fd2 c:be65542 \
    p:216be7f:docs/src/deriving/mwc_allostery.md p:216be7f:docs/src/deriving/cha_king_altman.md \
    p:216be7f:docs/src/getting_started.md p:216be7f:src/dsl.jl \
    p:216be7f:docs/src/identify/enumeration_engine.md p:216be7f:docs/src/deriving/re_vs_ss.md \
    p:216be7f:docs/src/developer.md c:81f9f1e c:fbdb9aa c:a84d194 c:26f7323 c:5cfd1e6 \
    c:32ea3c9 c:1b228ef c:09aec13 c:0b2021d c:761337e c:ad3ca91 c:eef3af7 c:7c2a657 \
    "e:Encode a mechanism's Sig=delete" "e:Have the \`Mechanism\` constructors check=delete" \
    "e:Take \`:OnlyA\` steps to their limit=$F" "e:\`eq_complexity_filter\` counts only=$F" \
    "e:The \`:OnlyA\` guard misjudges=$F" "e:Enumerate Cleland iso=$F" \
    "e:\`BigInt\` arithmetic inflates=$F" i
```

  Expected: no output. Departures in S4: the ROADMAP `BigInt` entry (Task 60); the
  messages of Tasks 57 (unit 2) and 58-61; ahead of the tip, the Stage checkpoints table's S4 row:
  aa158dd, 9578d2a, 8aed057, 02214d5, cd7a4ae, 4952d79, 919bc4c, 26f7323, 32ea3c9,
  5cfd1e6, 1b228ef, 09aec13, 761337e, ad3ca91, eef3af7, 7c2a657 and 81f9f1e, with S3i's
  list.

### Task 62: CLAUDE.md facts for Denis [S5]

**Replays:** `[ea6bf86]`: its own commit, never folded; the PR asks Denis to take it or
drop it (question 4). It corrects three statements of fact the code made false and
changes no rule: line 176 (the enumeration types are reaction graphs built from raw input,
not from `Step` / `Species`), line 189 (a focused run of the enumeration file takes about
1.5 minutes) and line 220 (edge direction, edge order and group order are canonicalized in
the `Mechanism` / `AllostericMechanism` constructors, `_canonicalize`; no `Step`
constructor exists). The texts are S3i section 8's, word for word.

**Files:**
- Modify: `.claude/CLAUDE.md`

**Interfaces:** none (no code).

**Gate:** `git diff --stat HEAD~1 HEAD` lists only `.claude/CLAUDE.md | 6 +++---`.

**Applied by the controller** (R63 item 6): a subagent may refuse to edit an instructions
file, so the controller runs this task's steps itself.

- [ ] **Step 1: Check that no earlier task touched CLAUDE.md**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git diff --stat 7205ba4 HEAD -- .claude
```

  Expected: no output.
- [ ] **Step 2: Apply**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary ea6bf86 -- .claude/CLAUDE.md | git apply
git diff --stat
```

  Expected: `.claude/CLAUDE.md | 6 +++---`.
- [ ] **Step 3: Commit.** K `ea6bf86` `.claude/CLAUDE.md` ("Correct three statements of
  fact in CLAUDE.md for Denis to take or drop").
- [ ] **Step 4: Gate.** `git diff --stat HEAD~1 HEAD`: `.claude/CLAUDE.md | 6 +++---`, one
  file. No suite: nothing the suite reads changed.

### Task 63: Final: the 1.10 Known issue, suites, census, docs, PR [S5]

**Replays:** `[7395164 + da2d8c7 + 2af6a00]` (the ROADMAP Known issue "On Julia 1.10 the
full suite can peak above 4 GB"; its figures describe the finished tree, S5 §9 item 8;
da2d8c7 keeps its cause to what the traces measured, R62 r1, and 2af6a00 scopes the
garbage to what the full-suite traces measured, R66 m2), then the final checks, the push
and the draft PR.

**Files:**
- Modify: `ROADMAP.md`
- Create (scratch, not committed): `$P/replay/trees/T63_{113,110,scan_112,scan_113}`,
  `$P/replay/scanenv/T63_{112,113}`, `$P/replay/pr_body.py`, `$P/replay/pr_body.md`

**Interfaces:**
- Consumes: every earlier task; `$P/proto/drafts/pr_body_phase3.md`; the remark-scan tools
  `$P/s5/probe/remarks_{derivation,identify}.jl`, `$P/gcprobe/repro_upstream.jl`,
  `$P/gcprobe/hoistscan.py`; S5's phase-2 scan trees and environments
  (`$P/s5/trees/p2_{112,113}`, `$P/s5/scanenv/p2_{112,113}`) and HEAD's environments
  (`$P/s5/scanenv/HEAD_{112,113}`), which seed this task's; the prototype's suite logs at
  its tip (`$P/proto/logs/S5m4_full_112_4b72021.log`; 4b72021 holds 7c2a657's tree) and at
  7395164, whose src and test differ from the tip's only by ad3ca91, eef3af7 and 7c2a657
  (`S5m_full_113_7395164.log`, `S5m_full_110_cfe9211_run{1,2}.log`; cfe9211 holds
  7395164's src and test).
- Produces: branch `graph-mechanism` on `origin`; a draft PR against `main`.

**Gate:** 1.12 40,109; 1.13 40,109; 1.10 40,101 pass and 1 broken when at least 5 GB is
available, else the prototype's 1.10 runs cited (R64); the docs build with one warning;
C 63 zero differences; the remark scan 0 hoists in every workload and 2 in each control;
checkpoint S5; the review's findings settled; the PR body with the edits R64 approved;
the draft PR open. Version stays 0.10.0.

- [ ] **Step 1: Apply 7395164, da2d8c7 and 2af6a00**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git show --format= --binary 7395164 -- ROADMAP.md | git apply
git show --format= --binary da2d8c7 -- ROADMAP.md | git apply
git show --format= --binary 2af6a00 -- ROADMAP.md | git apply
git diff --stat
diff ROADMAP.md $P/replay/roadmap_final.md && echo final-roadmap
```

  Expected: `ROADMAP.md | 11 +++++++++++`, then `final-roadmap`: the ROADMAP is the one
  the plan set out to reach (7c2a657's with 7205ba4's three P4-1 hunks).
- [ ] **Step 2: Commit** (message: 7395164's subject; its body in R62 r1's and R66 m2's
  wording, which da2d8c7 and 2af6a00 give the entry: 7395164's body names a cause no trace
  measured)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short   # M ROADMAP.md
git add ROADMAP.md
cat > $P/replay/msg <<'EOF'
Record the full suite's peak above 4 GB on Julia 1.10

On Julia 1.10, the compat floor, the suite peaked at 3.30-4.36 GB in seven
runs, against 2.6-3.0 GB on 1.12. The ter-ter init testset leaves about 1 GB
of garbage. On 1.10 the full suite still counts it after the enumeration
file, 1.83-1.84 GiB live there against 0.55-0.65 GiB with one GC.gc() after
the testset; with the file run alone, the count drops at the ter-ter seed
testset. No test value retains the garbage: a full collection after the
testset frees it (live bytes 0.94 to 0.055 GiB), and that GC.gc() brings the
suite to 2.93-3.05 GB. CI runs Julia 1 alone, so the entry records the peak
for a decision.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
git commit -F $P/replay/msg
git status --short; git log --oneline -1
```

  Expected: no status line after the commit; the log line ends `Record the full suite's
  peak above 4 GB on Julia 1.10`.
- [ ] **Step 3: Full suite on Julia 1.12** (background). F 63: 40,109, pristine; record its
  maxrss (the prototype's tip: 2.81 GB). `grep '^version' Project.toml` prints
  `version = "0.10.0"`, and `masked T63_full $P/proto/logs/S5m4_full_112_4b72021.log`
  prints nothing: the output equals the prototype tip's (4b72021 holds 7c2a657's tree),
  numbers masked.
- [ ] **Step 4: Full suites on Julia 1.13 and 1.10**, each on an archive of HEAD with the
  Manifest S5 resolved for its version (D-S5-8), so `$W` keeps its 1.12 Manifest

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
for v in 113 110; do
    D=$P/replay/trees/T63_$v
    rm -rf $D && mkdir -p $D && git archive HEAD | tar -x -C $D
    cp $P/s5/trees/HEAD_$v/Manifest.toml $D/
done
free -g
```

  The 1.13 suite always runs. The 1.10 suite runs only when `free -g` shows at least 5 GB
  available after the 1.13 run (it peaks at 3.3-4.4 GB, the Known issue Step 1 recorded);
  with less, it does not run, no `$L/T63_full_110.log` exists, and Step 11's PR body cites
  the prototype's 1.10 runs at 7395164 instead and says so (R64). Nobody asks Denis to
  close his language server. Background, in one call:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
(cd $P/replay/trees/T63_113 && $P/jl.sh $L/T63_full_113.log 3600 +1.13 --project=. \
    -e 'using Pkg; Pkg.instantiate(); Pkg.test(julia_args=["--heap-size-hint=2500M"])')
if [ "$(free -g | awk '/^Mem:/ {print $7}')" -ge 5 ]; then
    (cd $P/replay/trees/T63_110 && $P/jl.sh $L/T63_full_110.log 3600 +1.10 --project=. \
        -e 'using Pkg; Pkg.instantiate(); Pkg.test(julia_args=["--heap-size-hint=2500M"])')
else
    echo "1.10 not run: $(free -g | awk '/^Mem:/ {print $7}') GB available"
fi
```

  Expected (`$P/proto/logs/S5m_full_113_7395164.log`, `S5m_full_110_cfe9211_run1.log`,
  each 10 below this branch, which holds ad3ca91's three pins, eef3af7's five and
  7c2a657's two): 1.13 `EnzymeRates.jl | 40109  40109`, `exit 0`, JET 0.12; 1.10, when it
  ran, `EnzymeRates.jl | 40101  1  40102` (pass, broken, total: the fit-compile census's
  `@test_skip`, which needs 1.12's global method table), `Testing EnzymeRates tests passed`,
  `exit 0`, JET 0.9. Record each maxrss figure. Then the output check S5 ran, pristine on
  every version that ran:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
pristine T63_full_113; masked T63_full_113 $P/proto/logs/S5m_full_113_7395164.log
test -f $L/T63_full_110.log &&
    { pristine T63_full_110; masked T63_full_110 $P/proto/logs/S5m_full_110_cfe9211_run1.log; }
```

  Expected: `12`, then nothing; when 1.10 ran, `9` (1.10 skips three compile-budget
  `@info` lines), then nothing. Each log equals the prototype's on its version, numbers
  masked. A 1.10 run that ends with "received signal: KILL" is the Known issue's
  out-of-memory case, which proves nothing (R64): report it with the log, move the log to
  `$L/T63_full_110_killed.log` and go on as if 1.10 had not run.
- [ ] **Step 5: Docs build** (background)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
cd $W && $P/jl.sh $L/T63_docs.log 1800 --project=$P/replay/docenv --heap-size-hint=2500M docs/make.jl; rm -rf $W/docs/build
```

  Expected: as Task 61 Step 10: `exit 0`, the one deploydocs warning.
- [ ] **Step 6: Final census** (background). C 63, without `CENSUS_ADAPTER`: zero
  differences, order included (S5 minors: `census/logs/cmp_S5m_vs_S3b_end.log`).
- [ ] **Step 7: R18 remark scan.** R18 asks that no function callers loop over returns a
  Union of immutable, reference-holding parametric structs; the scan looks for the
  allocation hoist such a Union lets LLVM make (E6's bug), on derivation and identify
  workloads, at HEAD and at phase 2, on Julia 1.12 and 1.13, as
  `$P/s5/scripts/remarks.sh` does. Phase 2's trees and environments are S5's. HEAD's trees
  are made here, and HEAD's environments start from S5's (`$P/s5/scanenv/HEAD_{112,113}`,
  resolved with phase 2's on one registry state); Pkg runs offline and keeps every
  version, so both sides scan on S5's package versions:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
for v in 112 113; do
    D=$P/replay/trees/T63_scan_$v
    rm -rf $D && mkdir -p $D && git archive HEAD | tar -x -C $D
done
```

  Background (setup, then the scan; S5's runs took 20-60 s each, the setup longer when it
  precompiles):

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
for row in "112 +release" "113 +1.13"; do
    set -- $row
    E=$P/replay/scanenv/T63_$1
    rm -rf $E && mkdir -p $E
    cp $P/s5/scanenv/HEAD_$1/Project.toml $P/s5/scanenv/HEAD_$1/Manifest.toml $E/
    $P/jl.sh $L/T63_scanenv_$1.log 2400 $2 --project=$E -e "using Pkg; Pkg.offline(true); Pkg.develop(path=\"$P/replay/trees/T63_scan_$1\"; preserve=Pkg.PRESERVE_ALL); Pkg.instantiate()"
    diff <(grep -v '^path = ' $P/s5/scanenv/HEAD_$1/Manifest.toml) <(grep -v '^path = ' $E/Manifest.toml) &&
        echo "T63_$1: S5's versions"
done
export JULIA_LLVM_ARGS="-pass-remarks=julia-licm -pass-remarks-analysis=julia-licm"
for row in "HEAD_112 +release $P/replay/scanenv/T63_112 $P/replay/trees/T63_scan_112" \
           "p2_112 +release $P/s5/scanenv/p2_112 $P/s5/trees/p2_112" \
           "HEAD_113 +1.13 $P/replay/scanenv/T63_113 $P/replay/trees/T63_scan_113" \
           "p2_113 +1.13 $P/s5/scanenv/p2_113 $P/s5/trees/p2_113"; do
    set -- $row
    for w in derivation identify; do
        $P/jl.sh $L/T63_remarks_${w}_$1.log 3000 $2 --startup-file=no --pkgimages=no \
            --project=$3 --heap-size-hint=2500M $P/s5/probe/remarks_$w.jl $4
    done
done
for row in "112 +release" "113 +1.13"; do
    set -- $row
    $P/jl.sh $L/T63_remarks_control_$1.log 600 $2 --startup-file=no $P/gcprobe/repro_upstream.jl
done
unset JULIA_LLVM_ARGS
python3 $P/gcprobe/hoistscan.py $L/T63_remarks_*.log > $L/T63_remarks_scan.txt 2>&1 &&
    echo REMARKS_DONE >> $L/T63_remarks_scan.txt ||
    echo REMARKS_FAILED >> $L/T63_remarks_scan.txt
```

  Poll until `REMARKS_DONE` or `REMARKS_FAILED`, then `cat $L/T63_remarks_scan.txt`.
  `REMARKS_FAILED` means hoistscan.py failed, its error above the marker (a dry replay,
  which runs no Julia, gets it): stop and report it. The setup lines print
  `T63_112: S5's versions` and `T63_113: S5's versions` (the Manifests differ from S5's
  only in the developed tree's path); a `diff` line instead means Pkg moved a version:
  stop and report it. Expected
  (`$P/proto/logs/S5_remarks_scan.txt`): every derivation and identify log reads
  `hoisted allocations: 0  hoisted with a bits store into a reference field: 0` (escape
  analyses near S5's: HEAD 358 and 1,173 on 1.12, 776 and 1,779 on 1.13; phase 2 560 and
  1,457, 691 and 1,745); each control reads `hoisted allocations: 2  hoisted with a bits
  store into a reference field: 2` with `DANGER 2x`, its run aborting (`exit 134`): the
  positive control shows the scan can see a hoist. A hoist in a package workload: stop and
  report it with its remark lines.
- [ ] **Step 8: Checkpoint S5** (tip 7c2a657).
  `git diff --stat 7c2a657 HEAD -- src test docs ROADMAP.md Project.toml .claude` lists
  exactly: ROADMAP.md 20, the spec 48 (D1), the plan (D2), TE 54. Then:

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
export GIT_INDEX_FILE=$P/replay/ckpt.idx
rm -f $GIT_INDEX_FILE
git read-tree 7c2a657
git diff 4f9a34a 7205ba4 -- docs/superpowers/specs | git apply --cached                       # D1
git update-index --cacheinfo 100644,$(git hash-object -w $P/replay/roadmap_final.md),ROADMAP.md  # P4-1's hunks
git show :test/test_mechanism_enumeration.jl > $P/replay/ckpt.tmp
python3 $P/replay/r55m2.py $P/replay/ckpt.tmp                                                 # R55 m2
git update-index --cacheinfo 100644,$(git hash-object -w $P/replay/ckpt.tmp),test/test_mechanism_enumeration.jl
plan=docs/superpowers/plans/2026-10-16-graph-mechanism.md
git update-index --add --cacheinfo 100644,$(git rev-parse HEAD:$plan),$plan                    # D2
git diff --cached --stat HEAD -- src test docs ROADMAP.md Project.toml .claude
rm -f $GIT_INDEX_FILE $P/replay/ckpt.tmp; unset GIT_INDEX_FILE
git diff 7c2a657 HEAD -- ROADMAP.md | grep -c '^@@'
git diff --stat 7205ba4 HEAD | tail -1
```

  Expected: `indented lines 6977-7003` (the testset's place at the tip), no `git diff`
  output, `3` (P4-1's three hunks: Task 36's merged paragraph, the entry "Let an allosteric
  ping-pong stop both half-reactions ...", the iso entry's "Revisit allosteric ping-pong
  ..." sentence), and `44 files changed, ...` (the 43 of the File Structure and the plan).
  Departures over the whole replay: D1-D4, the ROADMAP entries written ahead, R55 m2's 27
  lines, and the commit messages each task writes out.
- [ ] **Step 9: Whole-branch final review.** Use superpowers:requesting-code-review over
  `7205ba4..HEAD`, giving the reviewer the spec, this plan's Review Focus (its five items
  are the inputs no test exercises) and its Global Constraints. Expected: no
  Critical or Important finding, since every commit is the reviewed prototype's. Any
  Critical or Important finding: stop before the push and report it to Denis with the
  review; this plan adds no code past the prototype's without him. List the Minor findings
  for Denis in the session; the PR body changes only by Step 11's edits.
- [ ] **Step 10: Push**

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
git status --short
git rev-parse --abbrev-ref HEAD
git push -u origin graph-mechanism
```

  Expected: no status output, `graph-mechanism`, and the push sets the upstream.
- [ ] **Step 11: Write the PR body**

The prototype's draft (`$P/proto/drafts/pr_body_phase3.md`) describes the prototype's
history: its Memory table cites local proto3 commits, it says the type swap raised the
peak (here Task 34's swap carries the fix), it names rulings R1-R59, it quotes the
prototype's run times, and it does not mention 7205ba4. `pr_body.py` writes the draft with
six edits that make it true of this branch, the times and peaks read from this task's
logs. R64 approved the six edits, on the draft that carries R62 r1's and R66 m2's
wording. When Step 4 did not run Julia 1.10, the closing figures say so and the
Verification line cites the prototype's two 1.10 runs at 7395164 and how its src and
test differ from this tip's (R64).

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
cat > $P/replay/pr_body.py <<'EOF'
# ABOUTME: Writes the phase-3 PR body: the prototype's draft with the edits that make it
# ABOUTME: true of the replayed branch, times and peaks read from Task 63's suite logs.
import os
import re
import sys

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # $P, two levels up
L = P + "/replay/logs"


def suite(path):
    """The Test Summary time and the peak in GB (maxrss_kb / 10^6) of a suite log."""
    text = open(path, encoding="utf-8").read()
    time = re.search(r"^EnzymeRates\.jl \|.* (\S+)$", text, re.M).group(1)
    peak = int(re.search(r"^maxrss_kb (\d+)$", text, re.M).group(1)) / 1e6
    return time, f"{peak:.2f}"


(t112, m112), (t113, m113) = [suite(f"{L}/{x}.log") for x in ("T63_full", "T63_full_113")]
rulings = open(P + "/RULINGS.md", encoding="utf-8").read()
last = max(int(n) for n in re.findall(r"^R(\d+)\.", rulings, re.M))  # the last ruling
if os.path.exists(f"{L}/T63_full_110.log"):  # Step 4 ran the 1.10 suite
    t110, m110 = suite(f"{L}/T63_full_110.log")
    peaks = f"peaked at {m112} GB on 1.12, {m113} GB on 1.13 and {m110} GB on 1.10."
    suites = (f"  {t112}; Julia 1.13.1, 40,109 pass in {t113}; Julia 1.10.12, 40,101"
              f" pass in {t110},\n  the fit-compile census skipped by design. Output"
              " pristine on all three.")
else:  # less than 5 GB free: the prototype's 1.10 runs at 7395164's src and test (R64)
    (ta, _), (tb, _) = [suite(f"{P}/proto/logs/S5m_full_110_cfe9211_run{k}.log")
                        for k in (1, 2)]
    peaks = (f"peaked at {m112} GB on 1.12 and {m113} GB on 1.13; 1.10 did not run here"
             " (Verification).")
    suites = (f"  {t112}; Julia 1.13.1, 40,109 pass in {t113}. Output pristine on both."
              " Julia 1.10.12\n  did not run here: less than 5 GB of memory was free. On"
              " the prototype at 7395164\n"
              f"  it passed 40,091 in {ta} and {tb}, the fit-compile census skipped by"
              " design,\n  output pristine; 7395164's src and test differ from this tip's"
              " only by ad3ca91,\n  eef3af7 and 7c2a657 (tag hints unwrap only a"
              " multiplicity call; ten more test_dsl.jl\n  pins) and the indentation of one"
              " list in test_mechanism_enumeration.jl.")
EDITS = [
    ("The type swap raised the full suite's peak RSS from 2.84-3.06 GB to 4.1-4.2 GB,"
     " above\nCLAUDE.md's budget of about 4 GB beside a resident language server. The"
     " cause: the\n",
     "In the prototype this branch replays, the type swap raised the full suite's peak"
     " RSS\n"
     "from 2.84-3.06 GB to 4.1-4.2 GB, above CLAUDE.md's budget of about 4 GB beside a\n"
     "resident language server. The cause: the\n"),
    ("and the heap hint were each measured and rejected as causes.\n",
     "and the heap hint were each measured and rejected as causes. This branch's swap\n"
     "commit already stores the vectors at their length, so none of its commits carries\n"
     "the higher peak. The bisect below ran on the prototype's local commits, which are\n"
     "not on this branch:\n"),
    ("Full suite on Julia 1.12 since the fix: 2.64-2.96 GB (S5 2.71-2.92 GB); on 1.13:\n"
     "2.59-2.65 GB; ter-ter init 13.0 s (31 s gate); the enumeration file alone 83.8 s"
     " at 2.28\nGB (145 s gate).",
     "In the prototype since the fix, the full suite peaked at 2.64-2.96 GB on Julia 1.12\n"
     "(S5 2.71-2.92 GB) and 2.59-2.65 GB on 1.13; ter-ter init took 13.0 s (31 s gate)"
     " and\nthe enumeration file alone 83.8 s at 2.28 GB (145 s gate). This branch's"
     " final suites\n" + peaks),
    ("- The prototype's rulings (R1-R59) record every decision",
     f"- The prototype's rulings (R1-R{last}) record every decision"),
    ("## Verification\n",
     "- The branch's first commit after main, 7205ba4 (\"Keep allosteric ping-pong"
     " reachable\n  under phase 4's admission\"), edits the spec and ROADMAP for phase 4:"
     " the tagging\n  move will stop one half-reaction per V-type child.\n\n"
     "## Verification\n"),
    ("- Full suite with Aqua and JET: Julia 1.12.5, 40,099 pass in 6m09s; Julia 1.13.1,\n"
     "  40,099 pass in 4m59s; Julia 1.10.12, 40,091 pass, the fit-compile census skipped"
     " by\n  design. Output pristine on all three.",
     "- Full suite with Aqua and JET on this branch's tip: Julia 1.12.5, 40,109 pass in\n"
     + suites),
]
body = open(P + "/proto/drafts/pr_body_phase3.md", encoding="utf-8").read()
for old, new in EDITS:
    if body.count(old) != 1:
        sys.exit(f"{body.count(old)} matches for {old[:60]!r}")
    body = body.replace(old, new)
sys.stdout.write(body)
EOF
md5sum $P/replay/pr_body.py
if python3 $P/replay/pr_body.py > $P/replay/pr_body.md; then
    diff $P/proto/drafts/pr_body_phase3.md $P/replay/pr_body.md
    grep -c 'claude.ai/code/session' $P/replay/pr_body.md
else
    echo "PR BODY FAILED"
fi
```

  (`pr_body.py` reads Steps 3-4's suite logs, so in a dry replay, which runs no Julia, it
  fails and the block prints `PR BODY FAILED` after the `md5sum`.)

  Expected: `785773365f5c6528ff7ddd8e641cfc49`, then the diff of the six edits (the Memory
  section's first paragraph, its note that the table is the prototype's local bisect, its
  closing figures with this branch's peaks, the rulings' range R1-R<n>, where R<n> is the
  last ruling in `$P/RULINGS.md` when the script runs, the 7205ba4 bullet, and the
  Verification suite line with this task's times; without a 1.10 run, the closing
  figures say 1.10 did not run here and the suite line cites the prototype's
  two 1.10 runs at 7395164), then `0` (no session link; the body closes with the
  "Generated with Claude Code" line). `PR BODY FAILED`, the script's error above it:
  stop and report it. Stop and ask Denis before Step 12 if a count of
  Steps 3-7 differs from the body's: 40,109 on 1.12 and 1.13, 40,101 pass and 1 broken on
  1.10 when it ran, the census's 81,498 mechanisms and 542 laws with zero differences,
  and 0 hoists in every workload of the remark scan. Times and peaks differ on every run;
  the body takes this task's. R64 approved these edits; a body that needs any other
  change goes to Denis before Step 12.
- [ ] **Step 12: Open the draft PR** (never merge it; the title is the one R64 confirms)

```bash
. /home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase3/replay/env.sh
gh pr create --draft --base main --head graph-mechanism \
    --title "Store each mechanism as a reaction graph that prints and reads back as its macro text" \
    --body-file $P/replay/pr_body.md
```

  Expected: the new PR's URL. To change the body later, use `gh api -X PATCH
  repos/{owner}/{repo}/pulls/<number> -F body=@<file>`: `gh pr edit` fails on the
  Projects-classic deprecation. Report the PR's URL.

## Open for Denis

### Decide before Task 0

Settled by R63; the plan carries each answer, and Task 0 Step 1 checks that R62-R67 are
recorded.

**1. Docs false for a window (S5 §9 item 3, R60).** The plan applies 216be7f by file and
ROADMAP entry, one task per piece, and S4's 81f9f1e whole in Task 61. A text that two
sub-stages make false waits for the later one, and a piece that names later code waits
for that code; S5 §9 item 3 asks for one piece per sub-stage, "each piece true at its
task", and R60 and CLAUDE.md want each doc edit in the commit that makes the old text
false. Settled by R63 item 1: the windows below stand. The final tree is exact, and Denis
squash-merges, so main never sees an intermediate tree; splitting 216be7f by hunk would
put text no review has read into Tasks 20-60 (its developer.md hunk at :71-155 alone
mixes the swap, S3a, S3e, S3f, S3g, S3h and S3i). Each window runs from the task that
makes the text false to the task that corrects it (S5 §3's site table, this plan's
placement):
- developer.md:14-22 (`AllostericEnzymeMechanism{CatalyticMech, CatSites, RegSites}`),
  :47-52 (the lift drops the regulators a reaction declares) and :57-61 (the lifts read
  type parameters; `catalytic_mechanism` and `catalytic_multiplicity` specialize): Tasks
  34-55.
- developer.md:71-95, the joint solve: `(reaction, groups, names)` from Task 42,
  `_inactive_groups` from Task 46, the Species content keys (`_form_index_key`,
  `_step_index_key`, `_form_base_key`, `_metabolites_key`) from Task 54; to Task 55.
- developer.md:97-152, the enumeration types: the `naming` cache from Task 21 (f7bd003);
  the `reaction` and `steps` fields, `Step` and `Species`, the `Step` constructor,
  `_canonical_step_direction`, `_orient_tied_steps` and the `_assert_*` checks from Task
  34; `_is_chemistry` of a Step from Task 42; `_step_kind` from Task 50; to Task 55.
- developer.md's boundary-check sentences (:104-106, :171-173) and :170-174's "the
  constructor rejects", from Task 20; :170-174's `_onlya_haldane_violation` on its groups
  from Task 48; :191 and :193-199 (a twin "bound to no competitive inhibitor", the gauge
  "bookkeeping over the steps") from Task 46; :242-247 (`_re_segment_extras`,
  `_bottomless_re_segment`) from Task 51; all to Task 55.
- developer.md's boundary-check sentence and `_check_well_posed`'s docstring
  (src/admission.jl) state the one-flip and `:OnlyA` conditions for every mechanism, not
  only allosteric ones: from Task 20 to Task 60 (32ea3c9).
- enumeration_engine.md:74 ("the rule the `Mechanism` constructor enforces" from Task 20,
  `_bottomless_re_segment` from Task 51) and :330 ("Adds one `AllostericRegulator`", from
  Task 21, 4d13802): to Task 51.
- ROADMAP's MWC entry ("the constructor admits a second free form", for every
  mechanism): Tasks 20-36 (216be7f's hunk with dc4035b's scoping).
- ROADMAP's `BigInt` entry: "runs the guard when built" from Task 20, and
  `_onlya_haldane_violation` named as the holder of the null space from Task 41; to Task
  60, where the entry, which names the parser, lands with S5's re-measured figures.
- Text ahead of its code: developer.md names "the parser" from Task 55, five tasks before
  Task 60 adds it.
- S4's docs (81f9f1e, Task 61): developer.md's "`show` ... still specializes" from Task
  59 (b8a604f), and how a mechanism displays in the identify tutorial (:193-216) and
  Getting Started (:41-42) from Task 59 (3464d63).
- The identify tutorial's `L` rule (be65542, 0b2021d): from Task 8 (51fad19, R10) to
  Task 9.
- CLAUDE.md:176 and :220 from Task 54 (Step and Species go), :189 from S3h-S3i's
  speedups; Task 62 corrects them for you to take or drop (PR question 4).
The re_vs_ss.md tier-3 sentence, false at C9, now lands with C9 in Task 3 (2b29fd2's
re_vs_ss.md piece), so it has no window.

**2. Fold placements kept against a report's note.** Settled by R63 item 2: each stands,
with the reason its task states, except a494c9c, which folds into Task 10.
- Task 48: S3g round 1's 01e2dbd, 704bac2, 9cece42 and 801cd5c land together after
  2c40783, though the note writes each fix where its code is first written, so Tasks 45-47
  write code that Task 48 rewrites. Following the note means hunk-level written-out
  patches and intermediate trees the prototype never built or tested.
- Task 55: b42cc90's test_types.jl piece holds m3 (the note: task 3) beside m1 (task 5),
  whose `step_of(ER._is_iso)` needs 14ef525's `_is_iso(::RawStep)`; m3 lands in Task 55.
  b42cc90's TD and TE pieces follow the note into Task 53.
- Task 40: ff5e59c (R53 m2: S3e §9 task 1, Task 36) lands in Task 40, since its TRG lines
  sit among later S3e commits'; its red now runs there on the residual-blind mutant.
- Task 60: cd7a4ae's TD and TRG pieces (S5 §9 item 5 puts m6 in the `show` task, Task
  59) stay with the parser, since they edit tests 0ec7d51 writes (TRG's apply fails at
  `_testhelper_census_allosteric` on Task 59's tree); its src/names.jl, FX, TE and TO
  pieces moved to Task 59.
- Task 4: c86f61c (S1 §9 task 3's start) lands after 8b7fb88, since its helper calls the
  three-argument `_orient_tied_steps` and its TJS hunk edits 8b7fb88's text. e74c8af's
  docstring reflow (S1 §9: task 3) lands there too, since its hunk needs 8b7fb88's
  signature; Task 3 marks the short line as Known.
- Task 10: a494c9c (the `p_mix =(` whitespace fix) folds into c80e857, which introduces
  the slip, so no commit carries it (R63 item 2); Task 11 applies c2cd9b3 alone.

**3. R55 m2.** Settled by R63 item 3: the S3g review's fix text indents the closing `]`
too, which would leave it at 8 spaces where every other `expected = [` list in TE closes
at 4, so the plan indents the 27 element lines only: `r55m2.py`, Task 49 Step 2
(6697-6723, 27), Task 63 Step 8 (6977-7003), the S3g and S5 checkpoint rows and every
TE checkpoint count (2 lower) follow.

### Decide before Task 63

Settled by R64; Task 63 carries each answer.

- **The PR body's edits** (Task 63 Step 11): `pr_body.py`'s six edits stand, on the draft
  that carries R62 r1's and R66 m2's wording; the suite line takes this branch's counts
  (40,109 on 1.12 and 1.13, 40,101 on 1.10).
- **Memory for the 1.10 suite**: Task 63 runs the 1.10 suite only when at least 5 GB is
  available, and otherwise cites the prototype's 1.10 runs at 7395164 (the same src and
  test but for ad3ca91, eef3af7, 7c2a657 and the listed departures) and says so in the
  PR body. Nobody asks Denis to close his language server.
- **PR title** (Task 63 Step 12): "Store each mechanism as a reaction graph that prints
  and reads back as its macro text".

### The PR's questions

The PR body (`$P/proto/drafts/pr_body_phase3.md`, with Task 63 Step 11's edits) asks them
in this order.

1. Delete three guards that cannot fire: `_assert_no_name_clash`, the `L` check in
   `IdentifyRateEquationProblem` and the clash check in `_dependent_param_exprs`? Each
   stays, with a docstring saying it cannot fire. (R42, D-S3f-6)
2. Delete the test-only accessors `catalytic_mechanism` and `catalytic_multiplicity` on
   `AllostericEnzymeMechanism`? (R54 m1)
3. CLAUDE.md enumeration rule 1: two boundary tests write an ill-posed mechanism inline
   through `Mechanism(_raw_mechanism(:(begin ... end), false))`. Permit it, or amend the
   rule? (D-S3a-10)
4. Take or drop Task 62's commit, which corrects CLAUDE.md lines 176, 189 and 220.
5. Prefix `mechanism_id` with a letter, so CSV.jl never reads an id as a number? The docs
   tell readers to pass `types = String`. (D-S4-11; Review Focus 2)
6. Add two `@nospecializeinfer` `text/plain` methods, so REPL display compiles once
   whatever the Sig? (D-S4-6; Review Focus 3)
7. The constructor gate: the median `Mechanism(raw)` takes 3.5-3.8 µs, under the test's
   10 µs but over the design's 3 µs trigger for a byte hash of form names, which was not
   built. (D-S3d-1)
8. JET 0.12 and the `m_OnlyA_prod` test fix (Task 1): main needs both. Keep them here or
   move them to a PR of their own?
9. Add one `GC.gc()` after the ter-ter init testset, so the suite stays under 4 GB on
   Julia 1.10? (Review Focus 5)

### What the plan leaves open

**Commits and messages.**

- Departures P1-1 and P1-2: c7698fd (Task 3) calls the two-argument `_step_text(s, rxn)`,
  which arrives with 8b7fb88 in Task 4, so Task 3 calls the one-argument form and Task 4
  restores it. They cancel at the S1 checkpoint. Settled by R63 item 4: no Review Focus
  entry.
- Task 34's message adds one paragraph for 924353b's field pins, 1ecb534's bridge guard
  and the Getting Started sentence, by the rule that a fold which changes what a commit
  does adds its sentence (it already names module ReactionGraph where f1aa056 said "task
  1"). Settled by R63 item 4: the paragraph stays.
- Tasks 14-28 and 44-63 copy most messages with `git show -s --format=%B <sha> | git
  interpret-trailers`, writing out only those a fold changes; Tasks 0-13 and 29-43 write
  every message out. Settled by R63 item 4: messages copied so are exact and suffice.

**Reds and counts.**

- One count rule now holds for every task (Global Constraints, Counts): a count that
  differs stops the task until its paths are checked, then the unit is redone or the count
  recorded and reported. Many counts stay derived: most red logs of S1-S3b were lost on
  2026-10-14, so the reds of Tasks 3-9 and the per-file counts of Tasks 4 and 8-28 come
  from reading the test diffs; Task 34's red for 050e16f's testset and Task 60's red (9
  failures, 102 errors) are inferred. The full-suite counts of Tasks 2, 3 and 5-8 were
  measured by the fix round on the dry replay's trees.
- Tasks 20, 21 (unit 2), 42 and 55 take their red as named runs (`run_named.jl`), as the
  prototype did: a whole-file red would error far past the tests the unit adds. Settled
  by R63 item 4: a named run is the red run R60 asks for.
- Reds that only a mutant shows now run with `$P/replay/mutate.py` (Tasks 32, 37, 40, 46
  and 47), and Task 38 shows 8c867d3's red on the one-message version, as S3e §9 asks.
  Two stay recorded, not rerun: Task 31's `NTuple{N, Int}` decoder (S3d1 §9 task 5; no
  script holds that mutant), and four of S3d1 §9 task 7's six mutants (forms ordered by
  packed numbers, free E not moved first, written group labels, a per-length decoder),
  whose `$P/proto/logs/S3d_mut_*.log` figures come from an earlier S3d1 tree, so they fix
  no expectation for Task 32. Settled by R63 item 4: they stay recorded from the
  prototype's logs.
- Task 27 does not replay 990cb6a's red: `S3c_red_free.log` came from an unsplit tree, and
  990cb6a calls FX helpers it adds itself, so a replayed red would differ from the log
  (R63 item 4).
- Task 40 unit 2 applies TO's path before FX and TRD to reproduce `S3e_c11_red.log` (the
  reader under test is test code), departing from "all test paths, then the red"; the
  order stays (R63 item 4).
- JET 0.12: from Task 1 the 1.12 environment resolves JET 0.12, which the prototype ran
  only at its tip. Task 0 Step 9's pre-flight shows it passing at all 15 stage tips
  (Aqua and JET 11 of 11; ReactionGraph targeted at S3d1's tip, 0 reports), so a JET
  report at an intermediate task would come from a commit between tips. Settled by R63
  item 4: JET 0.12 from Task 1; such a report stops that task for debugging.
- 050e16f's TRG testset lands at the swap (Task 34), before S3g and S3h port the moves; the
  prototype first ran it at S3i. The task stops if it fails. Its R52 figures will be lower
  than the prototype's, which lacked the copies; they are recorded, not gated, until Task
  55.

**Execution.**

- Runners: settled by R63 item 5. One Runners section after the Global Constraints
  defines every shorthand; Tasks 44-63's green and red runners, which differ from Tasks
  14-28's G and Rd, are Gt and Rt, and their TE-alone runner is TEA, since TA names
  `test/test_accessors.jl` in Tasks 29-43's counts.
- Gates beyond the commit map in Tasks 29-43: a census at Task 30 (it changes
  `RegulatorySite`) and JET on ReactionGraph at Tasks 30 and 32. Settled by R63 item 4:
  they stay.
- Task 62 edits `.claude/CLAUDE.md`. Settled by R63 item 6: the controller applies it,
  since a subagent may refuse to edit an instructions file.
- The evidence that the git blocks run clean lives under `$P/plan`: Tasks 14-28 in
  `$P/plan/check_p2`, Tasks 29-43 in `$P/plan/evidence/p3part` (`validate.py`,
  `val.log`, a private clone; moved there by R63 item 7), the fix round's moved pieces
  in `$P/plan/fix_counts` (`replay_applies.sh` on the dry replay's commits), and the
  R62-R64 update's in `$P/plan/evidence/r62_update`, the R66 update's in
  `$P/plan/evidence/r66_update` and the R67 update's in `$P/plan/evidence/r67_update`.
- A few prose lines run 93-100 columns, mostly long command or function names; left as
  is.

Settled by the assembly: the parts' task sections replace the skeleton's short task
entries, so each task's own text governs its red counts and signatures (Task 3 unit 2's
red, Task 2's unit order, Task 0's Manifest copy, Task 4's `_is_reactant_copy(x::Metabolite,
rxn)`, Task 17's reds, Task 34's message). Task 36 writes its ROADMAP text out in full.
The ROADMAP checkpoint counts follow Tasks 29-43's patches, which keep the file's final
newline: 32, 39 and 43 lines at S3d2, S3e and S3f; 48, 48, 52 and 105 at S3g, S3h, S3i
and S4.
