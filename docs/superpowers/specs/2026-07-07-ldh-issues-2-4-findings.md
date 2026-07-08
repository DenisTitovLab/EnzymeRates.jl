# ABOUTME: Findings for LDH v0.1.6-run Issues 2-4 (serial bottleneck, parameter
# ABOUTME: counts, futile enumeration). Investigation record only — no fix in this PR.

# LDH v0.1.6 run — Issues 2, 3, 4 findings

Investigation record from the 2026-07-07 LDH HPC run
(`docs/ldh_hpc_results/2026_07_07_results/`). **Issue 1 (allosteric `UndefVarError`)
is fixed on branch `allosteric-combined-constraint-solve`; Issues 2-4 are documented
here for later work and are NOT addressed in that PR.** Evidence: the run CSVs, the
`csv_summary.json` aggregate, and the move-behaviour experiments in the session
scratchpad (`task1.jl`, `task2*.jl`).

## Headline: the non-monotonicity is real, not an Issue-1 artifact

The central open question was whether the beam's frozen low parameter tier — the same
6-param equations re-derived every iteration — is a genuine property of the expansion
moves, or a downstream artifact of Issue-1 corrupting `fitted_params` on allosteric
mechanisms. **It is genuine.** Measured on a 900-mechanism all-non-allosteric pool
(correct kernel, immune to Issue 1):

- `_expand_split_kinetic_group` delta histogram `{-1: 22, 0: 22, +1: 2484, +2: 272}` —
  **22 real −1 edges**, every one a 7-param parent → 6-param child on a plain
  `Mechanism`. The vanished parameter is always a free-enzyme binding K
  (`K_NADH_E`, `K_NAD_E`, …): the split creates a mechanism where a single-symbol
  Wegscheider RE tie forces e.g. `K_NADH_E == K_NADH_EPyruvate`, so the parameter is
  correctly absorbed (the documented Pass-2 absorption,
  `thermodynamic_constr_for_rate_eq_derivation.jl:277-289`). These are distinct, correct
  6-param models.
- `_expand_re_to_ss`: `{0: 108, +1: 2828}` — monotone.
- Allosteric `_expand_change_allo_state` −1/−2 edges are also genuine, not corruption: a
  delta −2 child evaluates `rate_equation` to a finite value with no `UndefVarError`; all
  16 `change_allo_state` and 8 `re_to_ss` negative edges scanned were clean.

**Consequence:** fixing Issue 1 corrects the allosteric `UndefVarError` crashes and the
off-by-one allosteric counts, but will **not** make the moves monotone and will **not**
drain the frozen tier. Draining it requires dedup/memoisation of re-derived equations in
the beam, not a `fitted_params` fix. (This also refutes the earlier "split can't reduce
params" intuition — splitting *can*, by forcing a Wegscheider tie.)

## Issue 3 — parameter-count questions

- **14 six-param initial mechanisms** are legitimate ping-pong / iso topologies
  (covalent-intermediate residual + a `Kiso` isomerisation group = 6 groups vs the
  sequential bi-bi's 5). Not a bug.
- **Iteration-1 child range 6-8 (the "+2/+3")** comes solely from `_expand_to_allosteric`
  (`{+1: 232, +2: 40, +3: 4}` on the base). For LDH — which declares no regulators — this
  is a bare K-type `:OnlyA` on a binding group, not a V-type-with-regulator move
  (`_expand_add_allosteric_regulator` and `_expand_add_dead_end_regulator` emit zero
  children here).
- **The frozen 6-param tier** is the non-monotone `split` (7→6) above, re-swept every
  iteration by the advancing-target sweep
  (`identify_rate_equation.jl:746`, `target = max(target+1, min(keys(frontier)))`, which
  only increases yet pops every count ≤ target). Genuine, per the headline.

## Issue 2 — per-iteration single-core stall

Denis observed a multi-minute single-core phase after the (fully parallel) fitting stage.
Measurement **refutes** the CSV-write hypothesis: the whole post-fit serial write path is
<1 s even at the 440 MB / 15,481×122 shape. The real costs:

- `expand_mechanisms(parents, rxn)` runs on the master with no `pmap`
  (`identify_rate_equation.jl:702-703`): ~17-34 s at iteration 15, growing with the parent
  count — this is the "between the progress line and the next fit" gap.
- PASS-1 re-renders `rate_equation_string` for **every** child every iteration to compute
  `eq_hash` (~1.15 s each at real depth), even for the ~92 % that are inherited duplicates.
  Only the *fit* is memoised, never the render or the expansion.

The single-node benchmark tops out near 34 s; the remaining gap to the observed ~3.5 min is
most likely master-side `pmap` marshalling at ~1000-worker scale, which a single node
cannot reproduce — flagged as needing cluster confirmation.

**Directional fixes** (future PR): `pmap` `expand_mechanisms`; and stop re-processing /
re-writing inherited duplicates (Issue 4) — which shrinks the PASS-1 render load, the
dominant systemic cost.

## Issue 4 — futile enumeration cycle

`_beam_search` has **no global cross-iteration guard**. The only dedup is `unique!` within
one iteration (`:702`) and the `eq_hash` memo that gates the *fit* only (`:641, :532`).
`_ingest!` (`:591-594`) pushes every entry — including inherited ones — back onto the
frontier, so equation-duplicate-but-structurally-distinct mechanisms become parents and
re-expand. With `expand_mechanisms` a heavy fan-in DAG, the same child is regenerated
across many iterations. Quantified from the run CSVs: **140,826 full row-processings of
~22,156 distinct mechanisms (6.36×, up to 15×)**; iterations 16 and 17 are byte-identical
(13,143 distinct equations, 0 new) — a genuine fixed point that `while !isempty(frontier)`
never breaks. The run was killed by wall-clock, not convergence.

**Fix-direction analysis (for a future brainstorm):**

- **Deduping the frontier by `eq_hash` is UNSAFE.** Adversarial check: 6 real collision
  pairs (two structurally-distinct mechanisms sharing one `eq_hash`), each expanded, reach
  5-8 child equations the other misses (6/6). Keeping one per `eq_hash` silently drops real
  equations.
- **A global STRUCTURAL seen-set is safe and lossless** (`expand_mechanisms` is
  deterministic): it keeps every structurally-distinct mechanism, expands each once, and
  terminates.
- **Denis's stated preference** is to fix at the move level (like the #61 split self-loop
  drop) rather than add a guard, *unless proven infeasible*. The monotonicity result
  sharpens the tension: the −1 `split` edges are legitimate and produce distinct correct
  6-param models, so a move-level rule that simply drops them would lose coverage. The
  unresolved question that decides feasibility: **is each −1-split 6-param model also
  reachable via a monotone (+1-from-5-param) path?** If yes, the −1 edge is redundant and a
  move-level drop is lossless; if no, only a seen-set / memoisation is safe. This needs one
  more targeted experiment before choosing.

## Recommended next step

Issues 2, 3, and 4 are one coupled problem: the beam re-derives (re-enumerates, re-renders,
re-writes) equations it has already found. A future brainstorm should decide between the
structural seen-set and a move-level redundancy rule, gated on the monotone-reachability
experiment above, and pair it with `pmap`-ing `expand_mechanisms`. None of this blocks the
Issue-1 fix.
