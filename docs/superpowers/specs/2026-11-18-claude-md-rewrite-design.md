# CLAUDE.md rewrite: design

## Intent

Agents working on EnzymeRates under-invest in interface design: they add shallow helpers
and widen module interfaces without looking for a simpler one. The current `CLAUDE.md`
states the deep-module goal in one paragraph, buried under generic rules, and some of its
rules pull against each other. The rewrite makes lowering complexity the file's organizing
idea, ahead of further refactoring toward deep modules.

The design ideas come from John Ousterhout, *A Philosophy of Software Design* (2018). Four
ideas are borrowed from the pstack plugin; pstack itself is not installed.

## Goals

- State the goal first: a package where one module can be changed or fixed knowing only the
  other modules' interfaces.
- Make agents design an interface before they implement it, and bring alternatives to the
  maintainer when the interface is shared.
- Give a short playbook for each kind of task: feature, refactor, bug fix, optimization.
- Remove contradictions and the shouted, person-specific rules.
- Stay under 200 lines in one file.

## Non-goals

- No new test machinery for type stability or boundaries. A new public name already fails
  `test/test_boundaries.jl` (developer-docs table) and `test/test_type_stability.jl` (no
  row) until both are updated, and JET's analysis follows each row's calls into internal
  functions.
- No refactoring of `src/`. `Mechanisms` keeps its 25 public names for now.
- No rename to `AGENTS.md`, no split into `.claude/rules/`, no project skills.

## Anti-goals

- A file that restates what a test enforces. A test-enforced rule gets one or two lines
  naming the test.
- A rigid checklist that makes small fixes slow. The design step scales with the change: a
  private helper needs a thought, a shared name needs written alternatives.

## Approaches considered

| | A. One file | B. Core file plus path-scoped rules | C. Core file plus project skills |
|---|---|---|---|
| Layout | Everything in `CLAUDE.md` | Test-fixture rules move to `.claude/rules/` | Playbooks become skills |
| Cost | Every session loads the fixture rules | Claude Code only | A skill that does not trigger skips the design step |

Chosen: **A**. The maintainer's reason: under 200 lines needs no split.

## Decisions made with the maintainer

- Features: interface design first, then failing tests against that interface, then the
  implementation. Bug fixes stay strictly test-first.
- pstack is not installed. Borrowed: a playbook per task type; the selection rule of its
  `architect` skill; its refactoring sequence; its benchmark questions, rewritten for Julia.
- One file (approach A).

## Content of the new file, by section

Target: about 150 lines. Calm wording throughout: no capitals for emphasis, no "Rule #1".

1. **Goal: low complexity** (about 12 lines). Complexity is whatever makes the package hard
   to understand and change; it comes from dependencies and from obscurity, and it grows by
   small additions. The structure that answers it: six modules, each with a small interface
   and the complexity inside, so work in one module needs only the others' interfaces. The
   two anti-patterns, named: functions calling each other across the package, so that any
   change needs knowledge of all of it; and many small shallow functions and modules, so
   that there are many interfaces to learn.

2. **Module map** (about 14 lines). A generated table, one row per module and one for the
   top level: what the module does (first sentence of its docstring), its public names, and
   the modules it may call. Followed by the facts the table cannot hold: exports are the
   only promise to users; `compile_mechanism`, `test_rate_equation`, `Mechanism` and
   `AllostericMechanism` are usable but not exported; a mechanism's stored fields are part
   of `Mechanisms`' interface and are read with dots.

3. **Design rules** (about 14 lines).
   - A function or module is deep: a simple interface over substantial work. Interface
     complexity counts names, arguments, argument complexity, errors and preconditions.
   - A simple interface matters more than a simple implementation: take the hard case
     inside rather than passing it to callers through an argument, an option or an error.
   - Each design decision lives in one module. Two modules that know the same format or
     rule should be merged, or the knowledge moved behind one interface.
   - Build what is needed now, behind an interface general enough not to name its one
     current caller. (Replaces "YAGNI" beside "architect for extensibility".)
   - No pass-through functions, no one-caller wrappers, no arguments that only travel
     through. Do not split a function for length; split only when the parts can be
     understood separately.
   - Extract repeated code when the snippet is substantial and the helper's signature is
     simple. (Qualifies "work hard to reduce duplication".)
   - Fewer error cases: define the operation so the odd input is a normal case; validate
     user input once, at the entry point, with a message in the user's terms; trust types
     inside.
   - A value is an internal constant unless users have a real reason to change it.
   - Remove unused features; no backward compatibility without approval.

4. **Changing an interface** (about 8 lines). Before adding or changing a function used
   outside its file: write its docstring first; a docstring that must be long, or must
   describe the implementation, means the design is shallow. For a public name, an argument
   of one, or an exported function: sketch two or three different interfaces, compare them
   from the caller's side, pick the one where each fact has one owner, each task has one way
   to do it and callers need to know least, and bring the comparison to the maintainer
   before implementing. Every change looks for a public name or argument it can remove;
   removing one needs no approval and is reported.

5. **Playbooks** (about 24 lines).
   - *Feature:* design the interface (section 4); write failing tests against it; implement;
     re-read for dead code and shallow helpers.
   - *Refactor:* pin the behaviour with tests first; delete before adding; migrate every
     caller and remove the old name in the same change; keep the result only if a reader
     needs to know less than before, measured by public names, arguments and cross-module
     calls first and line count second; commit separately from fixes.
   - *Bug fix:* follow `superpowers:systematic-debugging`; reproduce in the simplest failing
     test; fix the root cause, one change at a time; leave the touched code simpler.
   - *Optimization:* measure before changing. Name what limits the time. Time with
     alternating runs and report median and range; check allocations (`@allocated`), run
     the JET rows, and run `test/test_compile_budget.jl` for compile time. Revert a change
     that does not measurably help unless it simplifies. Simplify the hot path before
     adding machinery.

6. **Julia style** (about 8 lines). Dispatch on types, not `isa` or tags (walking an `Expr`
   is the exception). Type-stable code and concretely typed fields. Type parameters keep
   fields concrete; they do not specialize per mechanism beyond `EnzymeMechanism{Sig}` and
   `AllostericEnzymeMechanism{Sig}`, because of compile time. A trait only where a second
   implementation exists. Closures over `@nospecialize` arguments are untyped on Julia
   1.13; use `Base.Fix2`. Where this file is silent the SciML style guide applies, and repo
   conventions win where they differ (single `_` prefix, plain `@testset`). 92 characters,
   4 spaces, match the surrounding code, no manual whitespace changes.

7. **Names, comments and docstrings** (about 16 lines). The current naming rules and their
   four examples, kept. Comments say what the code cannot: why, a chemical or numerical
   reason, a trap. No comment about history or change. A docstring gives what a caller needs
   and nothing about the implementation. `ABOUTME:` two-line file header kept.

8. **Guards** (about 22 lines). One short paragraph each, naming the test:
   module boundaries; type stability rows; the parameter-name chokepoint; Canonical Step
   Form; `rate_equation` and `loss!` performance (the numeric bounds stay, since they are
   the contract, with "stop and discuss before a change that would allocate or slow it").

9. **Tests** (about 20 lines). Test output is pristine; no mocked behaviour; never delete a
   failing test; `_testhelper_` prefix. The two Derivation-test rules and the three
   Enumeration-test rules, compressed but complete, since no test enforces them.

10. **Commands and workflow** (about 20 lines). The two commands and the four facts about
    them, as now. Read and update `ROADMAP.md`. Delete specs and plans in the implementing
    PR. Git: branch for the task, commit often, never skip a hook, no `git add -A` without
    `git status`, ask about uncommitted changes found at the start.

11. **Working with the maintainer** (about 6 lines). Say so when you do not know. Push back
    with technical reasons, or say it is a gut feeling. Stop and ask for a design decision,
    a public name or anything hard to undo; otherwise choose, say what you chose and what it
    costs if wrong. Architectural changes are discussed before implementation.

## Removed from the current file

"Rule #1" and the capitalized imperatives; addressing the maintainer by name; the
sycophancy and "Circle K" lines; "you'll be replaced"; the five-step TDD list (replaced by
the feature and bug-fix playbooks); the debugging rules that restate
`superpowers:systematic-debugging`; "inline single-use helpers" and "prefer minimal code"
(replaced by the design rules); the long forms of the guard sections.

## The one code change

`CLAUDE.md` gets the module table between two marker comments. `test/module_map.jl` gains a
function that builds the table from `src/` (it already reads each module's docstring
sentence, public names and layer), `test/update_module_map.jl` rewrites the block along
with the Mermaid diagram, and `test/test_boundaries.jl` fails when the block differs from
`src/`. Built test-first, with a small-tree fixture like the existing diagram test.

## Calls made by Claude, for the maintainer to check

- Refactors are judged by interface size first and line count second. This changes the
  earlier working rule that less code is the measure.
- The personal working rules stay in the repo file, worded for "the maintainer", so that
  cloud sessions load them.
- The module table lists every public name (about 50 in all), so that a new public name
  shows in the `CLAUDE.md` diff.
- "A trait only where a second implementation exists" narrows the request to prefer Holy
  traits.
- The SciML guide is a tie-breaker, not a second rulebook.
- `docs/src/developer.md` is unchanged except where it quotes a rule this rewrite changes.

## Left to the builder

The exact wording of each section, within the line targets; the table's column layout; the
name of the new helper in `test/module_map.jl`.
