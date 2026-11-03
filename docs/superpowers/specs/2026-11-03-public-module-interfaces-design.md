# Public module interfaces and revised project rules

## Intent

Denis wants the package to get simpler with every change, and wants `.claude/CLAUDE.md` to say
so. Today the file pulls the other way: it asks for the smallest change, forbids removing
comments, tells agents to add names to an interface, and requires tools that no longer exist.
Each module's interface is also declared with a homemade `const API` tuple that nothing in
`src/` reads, where Julia has `public`.

## Goals

- Declare each module's interface with Julia's `public`, keeping every boundary check.
- Make simplification the default in `CLAUDE.md`: deep modules, few public names, simple
  arguments, and refactoring as part of each change.
- Require Denis's approval for any change that grows a shared interface or alters an exported
  function.
- State conventions for dispatch, comments, docstrings and names.
- Remove rules that are stale or that conflict with these.

## Non-goals

- Cleaning up existing comments, docstrings or names in a sweep. They change when touched.
- Removing `identify_rate_equation` keyword arguments. That becomes a `ROADMAP.md` entry.
- Making non-public access an error. Julia does not offer this; the boundary test does the job.

## Anti-goals

- Losing a check: the layer rule, the "Used by" columns and the module map diagram all stay.
- Rules so strict that every small fix stops to ask Denis.
- Using parametric types as a default, which would add compile cost per mechanism.

## Approaches

Measured on this branch's code with the `API` lists declared public at runtime.
`improper_qualified_accesses` listed 81 own-module uses before that and 0 after, because it
reports only improper uses. Every use is available only through internals, which
[ExplicitImports issue #186](https://github.com/JuliaTesting/ExplicitImports.jl/issues/186)
says ties a tool to an exact minor version. The scanner feeds five checks and ExplicitImports
could replace one.

| Approach | Result |
|---|---|
| **A. `public`, keep the scanner (Denis's choice)** | `const API` goes, no new dependency, every check stays; test tooling stays about 634 lines |
| B. `public`, ExplicitImports for rule 1 only | A dependency to save a few lines |
| C. `public`, ExplicitImports, no scanner | An estimated 150–200 fewer lines, but no layer check, no "Used by" columns, and the diagram needs another source |
| `export` in submodules | `using ..Mod` would bring names in bare, hiding the `Mod.` that shows a call crosses modules |

## Part 1: `public` replaces `const API`

Decided by Denis: approach A, and the minimum Julia version rises to 1.11, because `public` is
a syntax error on 1.10.

- `Project.toml`: `julia = "1.11"`.
- Each `src/<Module>/<Module>.jl` replaces its `const API = (...)` with a `public` statement
  listing the same names. Macros are written `public @enzyme_mechanism`.
- `_testhelper_api` in `test/module_map.jl` reads the `Expr(:public, ...)` statement, whose
  arguments are plain symbols, so `_testhelper_api_name` is deleted. A module with no `public`
  statement is an error. It keeps parsing `src/` without loading the package, because the
  fixture tests run it on temporary trees.
- Fixtures in `test/test_boundaries.jl` move to `public`, and testset names and messages say
  "public names" rather than "API list". Rule 3's fixture that lists a duplicate name
  (`:_guard, :ok, :_guard`) keeps its intent.
- Rules 1 to 4, the docs tables, the diagram and "Mechanisms API is the data model" stay with
  the same meaning.
- `docs/src/developer.md`: prose that names `const API` or an `API` list changes to match.
- `test/test_compile_budget.jl`: the comment's Julia 1.10 baseline becomes the figure measured
  on 1.11. The budget does not change.
- `ROADMAP.md`: "On Julia 1.10 the full suite can peak above 4 GB" is rewritten with the 1.11
  measurement or retired if 1.11 stays under 4 GB, and the compat-bounds entry quotes the new
  `julia = "1.11"`.
- `src/Mechanisms/names.jl` gains its missing `ABOUTME` header, since `CLAUDE.md` requires one.

The change adds and removes no interface name: each module's public names equal its old `API`.

## Part 2: `CLAUDE.md`

Denis decided the content of each rule. The wording is my proposal, and his approval of this
document approves it, after which the builder copies it verbatim. These lists are closed: every
other rule stays, and a builder who finds another conflict reports it to Denis instead of editing.

### Added

**In "Writing code", replacing "YOU MUST make the SMALLEST reasonable changes":**

> Every change, bug fixes included, leaves the code it touches simpler, in its interface and
> inside the module, rather than taking the smallest patch that adds complexity. Stay within
> the code the change touches, and commit a refactor separately from a fix where the two can
> be separated.

**In "Writing code", replacing "YOU MUST NEVER throw away or rewrite implementations without
EXPLICIT permission":**

> Rewriting a whole implementation or module needs Denis's explicit permission. Simplifying
> the code a change touches does not.

**In "Code Style":**

> Dispatch on types rather than branching with `isa` or tags; walking an `Expr` is the
> exception. Use type parameters to keep fields concretely typed, not to specialize per
> mechanism: `EnzymeMechanism{Sig}` and `AllostericEnzymeMechanism{Sig}` exist because
> `rate_equation` is generated, and specializing further costs compile time that
> `test/test_compile_budget.jl` guards.

**In "Naming", first in the section:**

> Choose names with care. A name says what a function does, so that a reader, human or agent,
> can follow a call without opening its docstring. Follow Julia conventions: `lowercase`
> functions, `CamelCase` types, `!` on a function that mutates its arguments, a leading `_` on
> an internal name.

**"Code Comments", replacing the bullets on when to remove comments:**

> Avoid comments; let the code speak for itself. A comment adds information at a different
> level of abstraction from the code, higher (why this exists, the idea behind it) or lower (a
> subtle numerical or chemical detail). It never restates the code. Delete a comment that
> restates the code when you touch it.

The bans on history in comments and the `ABOUTME` rule stay.

**A new "Function docs" paragraph after "Code Comments":**

> A docstring tells a caller what they need in order to use the function: its arguments,
> result, units and errors. It does not describe the implementation, so that the
> implementation stays hidden. An internal function needs no docstring when its name and
> arguments make its job completely obvious.

**"Module boundaries (guard)", replacing the paragraph:**

> `src/` is six submodules layered `Mechanisms ← Derivation ← {Fitting, Enumeration, DSL}` and
> `Data ← Fitting` (`DSL` holds the macros and `parse_mechanism`; `Data` holds `RateData`),
> each in `src/<Module>/` with its interface declared by the `public` statement in
> `src/<Module>/<Module>.jl`; the top level (identify, the exports) uses all of them.
>
> We want deep modules: the fewest public names and the complexity inside the module, so you
> can work on a module knowing only the others' interfaces. Complexity counts public names and
> the number and complexity of their arguments, with the cleanest, simplest arguments as the
> goal. A value is an internal constant unless users have a real reason to change it, as they
> do for optimizer settings.
>
> Denis's approval comes first for adding a public name, adding an argument to one, and any
> change to an exported function. Removing or shrinking a non-exported shared name goes in
> with the change that finds it, and the report says so.
>
> `test/test_boundaries.jl` parses `src/` without loading the package and fails on four rules:
> (1) a name of another module (`Mod.name`, or `using ..Mod: name`) that is not public there;
> (2) a module naming a module its layer forbids; (3) an underscore name in a `public`
> statement or used across modules; (4) an `export` in a submodule. Call another module's
> function qualified (`Derivation.fitted_parameter_count(m)`); outside `Mechanisms`, read a
> mechanism's stored fields with dots (`m.steps[k].from`, `m.groups`). The Developer page lists
> each module's interface in tables and draws the module map as a Mermaid diagram generated
> from `src/` (each box's sentence is the first sentence of its module's docstring), and the
> same test fails when they disagree with `src/`: a change to a `public` statement, an export,
> a module docstring or a cross-module call updates `docs/src/developer.md` in the same commit,
> and `julia test/update_module_map.jl` rewrites the diagram.

### Changed

| Current rule | New text or action | Why |
|---|---|---|
| "YOU MUST ALWAYS STOP and ask for clarification rather than making assumptions." | "Stop and ask when the decision is Denis's: design, a public name or exported function, or anything hard to undo. Otherwise choose, say what you chose, and report what it costs if wrong." | Conflicts with "just do it"; Denis asked for the rewrite |
| Proactiveness: "The action would delete or significantly restructure existing code" | "The action would rewrite a whole implementation or module" | Otherwise the simplification rule stops at every refactor |
| Proactiveness: "Multiple valid approaches exist and the choice matters" | "Multiple valid approaches exist and the choice is Denis's to make" | Matches the rewritten rule on stopping |
| "Fix broken things immediately when you find them. Don't ask permission to fix bugs." | Add: "Anything else unrelated goes in `ROADMAP.md`." | Takes over from the journal rule on unrelated fixes |
| "never park on a background monitor — poll a log in a bounded shell loop" | "run the full suite in the background under `timeout 1200`, so a hang also ends in an exit notice; wait for that notice and read the result before reporting" | The harness notifies on exit, but a hung process never exits |
| `rate_equation` section: "Locally the 120-ns bound leaves about 4× headroom…" | Delete that clause, keep the bounds and the CI rerun advice | The timings go stale |
| "Data tables use a `group` column (not `Article`+`Fig`)" | Drop the parenthesis | Neither name is in the tree |

### Deleted

- The journal rules: "Our relationship" (two bullets), "Commit your journal entries", and
  "Learning and Memory Management". No journal tool exists; the memory system does this job.
- "Issue tracking", which requires TodoWrite. That tool does not exist in current sessions.
- "YOU MUST NEVER remove code comments unless you can PROVE they are actively false" and
  "Comments should explain WHAT the code does or WHY it exists", which the new comment rule
  replaces.

### Elsewhere

- `.gitignore` gains `.claude/settings.local.json`. "STOP and ask about untracked files" stays,
  and this file no longer triggers it every session.
- `ROADMAP.md` gains two Ideas: review `identify_rate_equation`'s keyword arguments (at least
  18, such as `eq_complexity_filter = 337`) under the new rule, and `DSL.parse_mechanism`, a
  public name that no other module in `src/` uses.

## Verification

- `git diff` of the `API` tuples against the `public` statements shows the same names.
- The full suite passes on Julia 1.13 under the stated command, with a pristine log.
- The full suite runs once on Julia 1.11 (`juliaup add 1.11`), since CI and the installed
  versions would otherwise never exercise the new minimum. Record its peak memory. A failure
  there is fixed if it is a real bug and otherwise reported to Denis; no budget or timing bound
  is loosened.
- A deliberately private cross-module call fails rule 1, checked by the existing fixtures.
- The docs build passes.
- Each symbol and path `CLAUDE.md` names exists in the tree.

## Delivery

One PR against `main`, because the new text describes the new checks. The local branch
`public-module-interfaces` stacks on `data-based-fit-search` (PR #99), since both change the `API`
lists and the same part of `CLAUDE.md`. Decided by Denis: this PR opens after #99 merges, and the
branch rebases onto `main` first. This spec is deleted in the PR that implements it.

## Left to the builder

The order of commits, and the wrapping of `public` statements within 92 characters.
