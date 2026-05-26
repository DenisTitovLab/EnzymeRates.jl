# Finish the Concrete-Types Refactor — Roadmap + Phase 1 Detail

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. **This work is exploratory, not mechanical — do NOT delegate whole tasks to unsupervised subagents** (a prior attempt reverted the architecture to escape a blocker, continuation spec §11.4). Drive directly or with tightly-scoped, monitored single-step subagents; escalate architectural blockers rather than working around them.

**Goal:** Finish the Symbol→concrete-struct refactor: decompose every remaining opaque-form fixture, delete the legacy Sig path, rewrite the enumeration internals off the opaque-Symbol working representation, and ship the PR — meeting the continuation-spec §3 success criteria (one struct family, no Symbol-string dispatch, perf gates green, ≤8,200 src LOC).

**Architecture:** The branch already landed Stages 6β–7d (`EnzymeReactionLegacy` gone, spec types now `_Raw*` scratch, DSL emits decomposed Species, README/CLAUDE.md updated). Two coupled cleanups remain plus the high-risk enumeration-internals rewrite. The new Step-based `Sig` is **lossy for opaque Species** (continuation spec §11.1), so every fixture must be *decomposed* (or deleted) — there is no opaque shortcut.

**Tech Stack:** Julia 1.10+, Test/Aqua/JET. Non-negotiable gates: `rate_equation` 0-alloc/<100ns (`test_rate_equation_performance`), 3 compile-budget gates, `test_chokepoint.jl` param-naming chokepoint, `scripts/check_test_integrity.sh main`.

**Base:** branch `refactor-to-concrete-types-instead-of-symbols`, tip `407dc1a` (red by construction — opaque rejection on, opaque fixtures not yet migrated; continuation spec §11.1).

---

## Roadmap — three buckets to "done"

Current status: **9,363 src LOC** (`wc -l src/*.jl`; main = 7,136; gate ≤ 8,200, so **−1,163 net still to remove**). Most of that reduction lives in Phase 2.

| Bucket | What | Risk | Success criteria it closes (continuation spec §3) | Plan detail |
|---|---|---|---|---|
| **Phase 1** | Decompose opaque fixtures (HK via `::Inh`, 6 lumped Segel by rename, delete Theorell-Chance), re-enable opaque rejection. **Does NOT delete the legacy Sig path** (moved to Phase 2 — see ③). | Low–medium | §3.1 #1 (partial — fixtures), #3, #4 | **Fully detailed** below + Tasks 0–9 of [`2026-05-25-finish-phase1-multisite-and-opaque-removal.md`](2026-05-25-finish-phase1-multisite-and-opaque-removal.md) |
| **Phase 2** | Rewrite the enumeration internals (`_expand_*` moves + topology backtracker + dead-end enumeration) off the opaque-Symbol working-rep; delete the 7 parse-back helpers AND the `_Raw*` working-rep pipeline. **Then** delete the legacy Sig path (2-arg `EnzymeMechanism(metabolites, reactions)`, `_mechanism_from_legacy_sig`, `_is_new_sig`, `_legacy_step_tuple` — 25 refs), now safe because nothing emits opaque Sig. **Exit gate: end-to-end bi-bi `identify_rate_equation` derives without the "not proportional" error.** | **High** | §3.1 #1 (completes "one struct family"), #2 ("form names are not parsed back into structure"), §3.2 (the bulk of −1,163 LOC) | **Sketch + spike gate** below; full task detail deferred to its own brainstorm once Phase 1 lands |
| **Final** | Dead-code sweep, LOC-gate check, §3 success-criteria verification, PR + push. | Low | §3.2 (LOC), §5 (PR) | Sketched below |

**Sequencing rule (continuation spec §10.3):** Phase 2 is high-risk and its exact shape depends on what Phase 1 leaves behind, so it is sequenced **last** and re-planned in detail only after the Phase 1 green gate. Do not start Phase 2 against the pre-Phase-1 codebase.

**Decisions locked in this brainstorm (2026-05-25):**
- **Multi-site binding syntax: `G6P::Inh`** role tag → `CompetitiveInhibitor(:G6P)`, real metabolite name preserved (so `concs.G6P` drives it). Chosen over an explicit `CompetitiveInhibitor(G6P)` constructor because `::` is already the DSL annotation idiom (`:: EqualRT`/`:: OnlyR` on steps, `G6P::OnlyT` on regulators) and the parse is unambiguous (`Expr(:(::), …)` never collides with a Species call). Untagged G6P takes its declared role (product); the tag's presence is the disambiguator.
- **Delete Theorell-Chance, with a `refactor-deleted-tests.md` §2.1 entry** at the deletion commit. Confirmed on the record: it is tested on `main`, so the log entry is required (reason: the middle step `EA + B ⇌ EQ + P` binds *and* releases in one transition; `Step` has a single `bound_metabolite` field, so it is un-representable; adjacent coverage via Ordered Bi-Bi + Ping-Pong).
- **Legacy-Sig deletion is sequenced AFTER the Phase 2 enumerator rewrite, not in Phase 1 (③).** Reason: the multi-substrate `identify_rate_equation` path routes enumerator output (opaque central complexes) through the new Step-based Sig, which is lossy for opaque Species (§11.1). Deleting the legacy Sig before the enumerator stops emitting opaque risks a window with neither a working legacy path nor a fixed enumerator. *Honest caveat:* it is **not yet verified** whether the enumerator currently depends on the legacy Sig path or already routes through the (lossy) new Sig — two reviewers contested this. The reorder is safe under **both** readings; **Phase 2's opening step must verify the dependency empirically** (build a bi-bi `IdentifyRateEquationProblem`, observe whether opaque central complexes hit "Cycle N not proportional") before any deletion.
- **Ping-pong fixtures migrate by pure rename, formulas preserved (②, spike-verified 2026-05-25).** A standalone spike of #10 derived 5 states / 5 steps / 9 params and matched Segel Eq. IX-228 to `rtol=1e-10` over 20 trials. No analytical formulas are dropped; no step expansion. See Phase 1 plan Task 6.

---

## Phase 1 — fixtures + rejection (DETAILED)

**Execute Tasks 0–9 of [`2026-05-25-finish-phase1-multisite-and-opaque-removal.md`](2026-05-25-finish-phase1-multisite-and-opaque-removal.md)** (revised 2026-05-25 to absorb this review round — it is the single source of truth, not reproduced here to avoid drift). The accepted review findings are already folded into that plan:
- **Task 3** now patches `_species_name_from_sig` (`types.jl:1429`) as the *primary* form-name producer (finding ①), with the ripple grep broadened to all `regulators:`/`competitive_inhibitors:` forms (finding ⑥).
- **Task 6** records the spike-verified ping-pong rename (②); formulas are preserved. Theorell-Chance deletion has on-the-record approval; write the §2.1 log entry in the same commit (it exists on `main`).
- **Task 8** is now *re-enable rejection only* — the legacy Sig path deletion moved to Phase 2 (③).
- **Perf/compile-budget/chokepoint gates run per-commit** during the red window (finding ⑦), not only at the Task 7 gate.

**Phase 1 green gate (its Task 7):** full suite green, `check_test_integrity.sh main` EXIT=0, perf/compile/chokepoint gates green, zero opaque enzyme-form step entries remain in any test file. **Phase 2 does not begin until this gate is green.**

**Intentional red window:** while opaque fixtures are mid-migration the *full* suite is red (Task 0 disables rejection to unblock the file; it stays red until the Task 7 gate). This is documented and deliberate on this unpushed branch — verify each step with the relevant *per-file* test and commit granularly; do not paper over it.

---

## Phase 2 — enumeration internals Symbol→struct (SKETCH + SPIKE GATE)

> **Full task-level detail is deferred to its own brainstorm/spec/plan, started only after the Phase 1 green gate.** This section scopes the work and defines the spike that gates it — it is intentionally not broken into TDD steps yet (continuation spec §10.3 sequences it last because its shape depends on Phase 1's output).

**Problem.** After Phase 1, the public surface is fully decomposed, but the enumeration *pipeline* (`src/mechanism_enumeration.jl`) still manipulates an **opaque-Symbol working representation** internally and parses form-name Symbols back into structure. That violates continuation-spec §3.1 #1 (one struct family, no parallel representations) and #2 (form names are *not* parsed back into structure). The 7 live helpers that encode this representation (caller counts as of `407dc1a`, all in `src/mechanism_enumeration.jl`):

| Helper | Refs | Role |
|---|---|---|
| `_dead_end_form_name` | 7 | Synthesize opaque dead-end form Symbol |
| `_subtract_atoms` | 5 | Atom bookkeeping on opaque forms |
| `_can_pingpong` | 4 | Ping-pong eligibility from opaque atoms |
| `_atoms_dict` | 3 | Atom inventory keyed by form Symbol |
| `_parse_bound` | 2 | Parse bound metabolites *back out* of a form Symbol |
| `_bound_mets_from_form_name` | 2 | Same, higher-level |
| `_is_estar_form` | 2 | Estar-conformation test on a Symbol |

(`_form_name` is already dead — 0 refs — quick delete, can fold into Phase 1's final sweep or open Phase 2.)

**Target shape.** The `_expand_*` moves, the topology backtracker (`backtrack!` / `_catalytic_topologies`), and dead-end enumeration operate on `Step` / decomposed `Species` (bound lists of `Metabolite` subtypes) throughout, reading atom inventories from `ReactantAtoms` on the reaction — never synthesizing or re-parsing a form-name Symbol. The 7 helpers AND the `_Raw*` working-rep pipeline are deleted as their callers disappear.

**Phase 2 opening step — establish the baseline (③):** before any rewrite, build a bi-bi `IdentifyRateEquationProblem` and run `init_mechanisms`/`expand_mechanisms`/`compile_mechanism` on it. Record whether multi-substrate enumerator output (opaque central complexes) currently derives or fails with "Cycle N not proportional". This (a) confirms whether the bi-bi path is broken today, (b) tells us whether the enumerator depends on the legacy Sig path or already routes through the new (lossy) one — the contested question behind ③. **Turn that bi-bi case into a test: it is Phase 2's EXIT GATE** (must derive cleanly before Phase 2 is "done").

**Spike gate — scope is contested, resolve at the Phase 2 brainstorm (⑤):** reviewer B argued the `_expand_*` moves already dispatch on `Mechanism` structs and the opaque dependency is confined behind a `Set{Symbol}`/`_bound_at_forms` boundary, so rewriting one move (`_expand_add_dead_end_regulator`, highest fan-out — drives `_dead_end_form_name`) behind a shim is feasible. Reviewer A argued the move's *input* is opaque (from `init_mechanisms` via the `_Raw*` pipeline), so a single move can't be isolated without also changing the producer. **Both read the same code and disagreed — so the first design question of the Phase 2 brainstorm is: is the real spike boundary one move, or the `_catalytic_topologies`/`_stepspec_from_step` *production* boundary?** Whichever spike is chosen, prove on `test_mechanism_enumeration.jl`: (1) enumeration counts unchanged (bi-bi=11, ter-ter=283 — CLAUDE.md), (2) `dedup!` canonicalizes identically (struct `==`/`hash`), (3) no `init_mechanisms` perf regression (compile-budget gate 750). **If the working-rep proves load-bearing in a way structs can't express (as opaque Species were for the Sig, §11.1), STOP and escalate** — do not work around it.

**Dead-end form naming (⑥):** Phase 1's `::Inh` rename gives competitive-inhibitor-bound DSL/derivation forms an `inh` suffix, but `_dead_end_form_name` synthesizes enumeration forms *without* it. Unify the convention here (the enumeration should produce the same names the derivation expects) when the dead-end enumeration is rewritten.

**Then — delete the legacy Sig path (was Phase 1 Task 8):** once the enumerator emits decomposed Species and the bi-bi gate is green, delete the 2-arg `EnzymeMechanism(metabolites, reactions)` constructor + dead validators, `_mechanism_from_legacy_sig`, the 12 `_is_new_sig` accessor branches (collapse to the new-shape body), `_is_new_sig`, dead `_legacy_step_tuple`/`_species_sym`; rename `_legacy_step_tuple_from_sig` → `_step_tuple_from_sig`. Handle the 4 `@test_throws` on the 2-arg constructor in `test/test_types.jl` per §2.1.

**Expected LOC (④, corrected):** the −1,163-to-target reduction does **not** come from the 7 helpers alone (~100 LOC). The bulk is the `_Raw*` working-rep pipeline — `_RawSpec`/`_RawStep`/`_RawAllostericSpec` + `_init_raw_specs`/`_raw_from_mechanism`/`_mechanism_from_raw`/`_stepspec_from_step`/`_apply_equivalence_grouping` — plus the legacy Sig deletion (~310). **Re-baseline the true Phase 2 deletion LOC at the Phase 2 brainstorm and confirm ≤8,200 is reachable; if not, renegotiate the gate with Denis (continuation spec §3.2 is the only gating LOC target) rather than letting it force a hack.**

---

## Final — sweep, LOC gate, PR (SKETCH)

After the Phase 2 green gate:

- [ ] **Dead-code sweep.** `grep -rnE "_is_new_sig|_mechanism_from_legacy_sig|_legacy_step_tuple" src/` → empty. `grep -rn` each of the 7 Phase-2 helpers + `_form_name` → empty. Any `@generated` dual-Sig branch collapsed.
- [ ] **Success-criteria §3 verification** (continuation spec):
  - §3.1 #1 one struct family — no `_Raw*`-vs-`Mechanism` or opaque-vs-decomposed parallel reps remain in `src/`.
  - §3.1 #2 no Symbol-string dispatch — `test_chokepoint.jl` green; no form-name parse-back helpers.
  - §3.1 #3 test integrity — `bash scripts/check_test_integrity.sh main; echo "EXIT=$?"` → 0 (no pipe).
  - §3.1 #4 perf — `test_rate_equation_performance` 0-alloc/<100ns green; 3 compile-budget gates green.
  - §3.2 LOC — `wc -l src/*.jl` ≤ 8,200.
- [ ] **Full suite green** — `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/out.log 2>&1; echo "EXIT=$?"` then **read the `Test Summary` line** (the trailing `echo` swallows `Pkg.test`'s exit — trust the summary, not the notification exit code).
- [ ] **Docs final pass** — README + CLAUDE.md already updated (7e.2/7e.3); re-read for any opaque-form / legacy-Sig references that the deletions invalidated. Update continuation spec §11 status to "done".
- [ ] **PR.** Branch is unpushed (origin at `d638636`). Push, open PR with honest LOC numbers (celebrate the architectural simplification; the residual +LOC vs main is the cost of structured types — link the deferred parameter-naming refactor as the next step, continuation spec §7). End the PR body with the Claude Code attribution footer.

---

## Conventions (every commit, both phases)

- **TDD**: failing test → implement → green. **No `--amend`**; stay on the branch.
- **Test integrity** (continuation spec §2/§4): no deletion/weakening without a `docs/superpowers/refactor-deleted-tests.md` §2.1 entry *in the same commit*. `bash scripts/check_test_integrity.sh main; echo "EXIT=$?"` must be 0 (no pipe — the pipe masks the exit).
- **Full suite** check: `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/out.log 2>&1; echo "EXIT=$?"` then read the `Test Summary` line. Per-file iteration uses the temp-env recipe (resolve deps once):
  ```bash
  julia --project=. -e '
    using Pkg; Pkg.activate(temp=true); Pkg.develop(path=".")
    Pkg.add(["Test","Aqua","JET","OptimizationBBO","OptimizationPyCMA","OrdinaryDiffEqFIRK","Tables","DataFrames","Statistics","Optimization","Random","CSV"])
    using Test, EnzymeRates
    include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")'
  ```
- **Perf/chokepoint gates** — run `test_rate_equation_performance`, the 3 compile-budget gates, and `test_chokepoint.jl` **per-commit** (they build their own small mechanisms and don't need the red fixtures to load), so a regression is caught at its origin commit rather than bisected from a later gate (⑦).
- **No temporal-context comments** in code (no "Stage N", "previously", "legacy", "will be") — evergreen comments only. Plan/spec docs may reference stages.
- Commit footer: `src delta: -X / +Y net Z, cumulative: ±W` (`wc -l src/*.jl`; cumulative vs main 7,136). End messages with `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>`.

---

## Self-review notes
- **Riskiest work is Phase 2** (enumeration rewrite) — gated by a spike on the single highest-fan-out move before full commitment, with an escalation clause if the Symbol working-rep proves load-bearing.
- **Riskiest Phase 1 task is its Task 3** (role-distinct naming ripples to every competitive-inhibitor form); de-risked by its Task 1 `::Inh` spike.
- **Don't trust `Pkg.test` notification exit codes** — read the summary line.
- **Don't hand exploratory migration to unsupervised subagents** (continuation spec §11.4).
