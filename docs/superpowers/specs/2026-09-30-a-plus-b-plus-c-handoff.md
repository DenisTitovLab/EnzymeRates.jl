# A+B+C handoff: state after A and B

Date: 2026-09-30. Read this, then the findings document
(`2026-09-28-identifiable-enumeration-findings-and-plan.md`), before starting C.

## Where things stand

- Branch `step-explicit-stoichiometry` holds all of A+B+C. Denis wants A, B and C in one pull
  request from this branch (2026-09-30), so B and C are built here. Nothing is pushed.
- Done on the branch, full suite 54,642/54,642:
  - A, explicit-stoichiometry `Step` (spec `2026-09-28-step-explicit-stoichiometry-design.md`,
    plan `../plans/2026-09-28-step-explicit-stoichiometry.md`, regression record
    `2026-09-28-step-explicit-stoichiometry-regression.md`).
  - Denis's eight review follow-ups (plan `../plans/2026-09-29-step-refactor-review-followups.md`):
    - `bound_metabolite(s)` returns a pure binding's metabolite, `nothing` otherwise.
    - A kinetic group holds one kind of step with one RE/SS flag (`_assert_uniform_groups`).
    - A reaction appears in one step only (`_assert_each_reaction_once`).
    - Step constants are named by their reaction: `k_X_to_Y` for SS steps and `K_X_to_Y` for RE
      steps, a pure binding in the release direction (`K_ES_to_E_S`, a dissociation constant).
      The regression record shows this changed names only.
    - Internal functions are documented with docstrings.
- B, exact filters, is implemented (spec `2026-09-30-exact-filters-design.md`, plan
  `../plans/2026-09-30-exact-filters.md`, regression record
  `2026-09-30-exact-filters-regression.md`). Next: C, merged and Theorell–Chance seeds and
  canonical merges (findings plan section C), starting with its spec.

## Open questions parked for Denis

- `Kon`/`Koff` now render exactly like `Kfor`/`Krev`, and nothing dispatches on the difference.
  Merging the pairs would delete three copies of a four-way switch but changes the internal
  `Parameter` type family.
- `rate_equation_string` prints ping-pong residual names such as `K_EB_res_+A_-P_…`, which are
  not valid Julia. Older names had the same characters; src never parses the strings back.
- The `Estar`-style ping-pong test seeds (five copies in `test/test_mechanism_enumeration.jl`)
  are not atom-balanced. Direct move tests do not check atoms, so they pass.
- Group-rule errors number kinetic groups in canonical order, not the order written in the DSL.

## Working on this machine

- 7.7 GB RAM, no swap: run one Julia process at a time; the VS Code language-server process is
  fine.
- The full suite takes about 14 minutes. Start it in the background and poll its log in bounded
  foreground loops. A subagent's background job dies when the subagent hands back, so an agent
  that starts the suite must poll it to the end itself. Check that a process is alive before
  reporting it as running, and match `bin/julia` in `pgrep` so the pattern does not match the
  polling shell itself.
- A subagent can exhaust its context on a large change (one did at about 270k tokens). Split big
  tasks, and have agents read files by line range, not whole.
- Regression method for B and C: enumerate a fixed population, key each mechanism by its steps
  (a name-independent key), and compare fitted sets and Reduced rate strings. Evaluate the
  strings with a small expression interpreter rather than compiling each one. For code that
  predates a change, run it from a detached worktree. The scripts used for A lived in session
  scratch space and are gone; the regression record describes the method.
