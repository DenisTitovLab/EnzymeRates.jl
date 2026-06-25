# Design note: enumeration redundancy and a loss-dominance pruning filter

Status: proposal, backed by the empirical investigation summarized below. Not yet
implemented.

## Summary

The mechanism enumeration generates many **redundant or structurally
non-identifiable** mechanisms — candidates that carry one or more parameters the
data cannot constrain. They are fit anyway (wasted optimizer effort, occasional
ill-conditioning) and then lose model selection on parsimony. This note
characterizes where the redundancy comes from, argues that the current
`fitted_params` behavior is correct and should be kept, and proposes a cheap
**loss-dominance** beam filter to stop expanding mechanisms that do not improve
on the best simpler model. The filter's key safety property — that it never
discards the best mechanism at any parameter level — held in every complete
landscape tested.

## Where the redundancy comes from

A mechanism is **structurally non-identifiable** when `rank(∂v/∂θ) < #fitted_params`
— some parameter (or parameter combination) leaves the rate function `v`
unchanged. Measured over the non-allosteric `bi_bi` enumeration (~564 deduped
mechanisms from `init_mechanisms` + 3 expansion rounds), **~36%** are
non-identifiable. None at the minimal parameter count (the init mechanisms are
all clean); the redundancy lives entirely in the expanded, higher-parameter
mechanisms (51% at nfit=6, 42% at 7, etc.) — exactly the candidates that are most
expensive to fit. For allosteric mechanisms the fraction is far higher
(~79% for `bi_bi`, ~98% for a `uni_uni` + competitive-inhibitor enumeration): MWC
variants carry an allosteric constant `L` and inactive-state parameter mirrors
that substrate-only data largely cannot resolve.

The non-identifiability splits into two kinds, in roughly even proportion:

1. **Futile-cycle / "zero-column" (~56%).** A parameter is *literally absent*
   from `v` (`∂v/∂p ≡ 0`). Every measured case (72/72) traces to an **SS step
   whose two endpoints lie in the same rapid-equilibrium group** — a futile cycle
   inside an RE segment. Cha lumping correctly places that step's flux on the
   discarded diagonal of the King–Altman matrix, so the step's rate constants
   drop out of `v`. The rate equation is *right*; `fitted_params` (derived from
   the thermodynamic Wegscheider/Haldane null space) reports the constant as
   "independent" even though the rate cannot depend on it. These mechanisms are
   exactly the ones reachable by a rapid-equilibrium → steady-state flip that
   does **not** increase the number of RE-connected enzyme-species groups
   (`ΔG = 0`).

2. **Steady-state-intermediate lumping / "lumped" (~44%).** Two consecutive SS
   steps share a steady-state intermediate (e.g. `EA + B → EAB` SS feeding
   `EAB → EPQ` SS). The intermediate's rate constants enter `v` only in a fixed
   combination — visible as a constant-sum signature in the denominator, e.g.
   `(k_EAB_to_EPQ + k_off)`. The parameters *are* in `v`, but one degree of
   freedom among them is unidentifiable.

A key non-result: **`ΔG` (the change in RE-group count) is not the right
criterion** for identifiability. `_expand_re_to_ss` flips a whole metabolite
kinetic group, which always splits an RE group (`ΔG ≥ 1`) — yet it can still
produce a non-identifiable mechanism of kind (2). The correct test is
`rank(∂v/∂θ)`, which is also what catches both kinds uniformly.

## `fitted_params` over-counting is correct — keep it

For a futile-cycle mechanism, `fitted_params` reports a parameter that `v` cannot
depend on. This is **not** a bug to fix: a non-flux-carrying SS step inside an RE
segment is a nonsensical, over-specified mechanism, and the extra parameter in
the count is the signal that flags this to the user. The decision is to keep the
thermodynamic parameter count as-is.

## A representability caveat that dissolved an apparent gap

An intermediate version of this analysis classified a degenerate mechanism's
rate *law* as "uncovered" when no lower-parameter mechanism shared its exact set
of concentration monomials, and flagged ~72 such `bi_bi` laws as a potential
enumeration-completeness gap. **This was a false positive.** Fitting is the
ground truth for representability, and on noiseless data every "uncovered"
mechanism tested was fit to machine precision (~1e-13) by a *lower*-parameter
mechanism — one with a *superset* of monomials whose extra terms vanish at the
fitted parameters. So those functions are representable more parsimoniously after
all; the genuine completeness gap is negligible. The lesson:
**fitting/representability, not monomial-set matching, is the right notion of
"covered."** It also means the redundant mechanisms are genuinely dominated by a
simpler fit — which is what makes the loss-based filter below safe.

## Proposal: a loss-dominance beam filter

Rather than structurally pruning non-identifiable mechanisms (which would require
a per-candidate `rank(∂v/∂θ)` and risks edge cases), keep every mechanism but add
a third beam threshold, alongside the existing `loss_rel_threshold` and
`loss_abs_threshold`, that stops expanding a mechanism which fails to improve on
the best model with one fewer parameter.

Proposed keyword: `loss_m1param_threshold`. A mechanism at parameter count `n`
qualifies for the next-level beam only if, **in addition to** the existing
criteria, its loss satisfies

```
loss ≤ loss_m1param_threshold * minimum(loss of all (n-1)-parameter mechanisms)
```

This is an **AND** condition layered onto the current beam selection
(`_select_beam` in `identify_rate_equation.jl`). The required quantity is already
tracked: `best_loss_by_count[n-1]`. The base parameter count (no `n-1` level) is
exempt. Suggested default: a no-op (`Inf`) so existing behavior is unchanged and
the filter is opt-in; the evidence below supports an active value near `1.0`
("an added parameter must earn its keep by beating the best simpler model").

Rationale: a redundant or non-identifiable mechanism achieves the same best loss
as its simpler, identifiable equivalent, so it cannot beat the best `(n-1)`-param
model and is pruned — without any rank computation, reusing losses already
computed. It also subsumes the "parameter that doesn't pay off" case, which a
pure identifiability gate would miss.

### Why not a structural `rank(∂v/∂θ)` gate

A `rank < nfit` pre-fit gate was considered. It is principled but has two
drawbacks the loss filter avoids: (a) it needs a per-candidate sensitivity
Jacobian; (b) a *blanket* rank prune would also remove the rare mechanism that is
the sole realization of a legitimate rate law — and the "uncovered" analysis that
motivated that concern turned out to be a false positive (see above). The
loss-dominance filter handles redundancy through parsimony instead, and never
removes a function that actually fits.

## Empirical support: the filter never loses the per-level best

The filter is only safe if pruning a level-`n` mechanism with `loss > best(n-1)`
never makes the global-best mechanism at some level unreachable — i.e. if a
"parallel all-downhill-in-loss" path to every good mechanism always exists. The
failure mode would be **parameter epistasis**: two structural additions that help
only together, each one alone fitting worse than `best(n-1)`, so both single-step
intermediates get pruned and the target becomes unreachable.

This was tested on **complete** mechanism landscapes (so that `true-best(k)` and
the expansion DAG's parent edges are real, not sampled), comparing
`true-best(k)` = min loss over *all* level-`k` mechanisms against
`reached-best(k)` = the min loss the filter actually sees. Divergence at any `k`
would be a counterexample.

| reaction | landscape | result |
| --- | --- | --- |
| `bi_uni` | complete, 70 mechanisms | holds, 8/8 true-mechanism landscapes |
| `uni_bi` | complete, 70 mechanisms | holds, 5/5 |
| allosteric `uni_uni` (+R) | complete, 197 mechanisms | holds, all tested |

`reached-best(k) == true-best(k)` at every level in every landscape — **zero
divergence across ~22 landscapes**. This includes the genuinely hard cases:
non-allosteric mechanisms whose true function requires the maximum parameter
count, with loss improving only gradually (e.g. `6e-3 → 3e-3 → 8e-6 → ~0` across
nfit 4→7). Fits reached machine precision (~1e-13–1e-16), so the verdicts are
trustworthy. No parameter epistasis was observed.

### Caveats

- Every *complete* landscape that can be fit exhaustively and reliably is a
  **small** reaction (≤ ~284 mechanisms). `bi_bi` and `ter_bi` complete
  landscapes to nfit 7 are too large to fit exhaustively; a *sampled* landscape
  was deliberately not used because missing parent edges would manufacture false
  divergences. Parameter epistasis on richer systems is not ruled out, only
  unobserved.
- The allosteric landscapes, while passing, are weak tests: the `uni_uni`
  enumeration has no identifiable high-parameter MWC mechanism (only 4/197 are
  identifiable, all low-parameter), so every allosteric rate law is simple and
  the landscapes are trivially easy. The meaningful stress tests are the deep
  non-allosteric ones.

## Related follow-ups surfaced by this investigation

- **Allosteric derivation bug.** Deriving/fitting at least one allosteric
  mechanism raises `UndefVarError: k_I_EP_to_ES not defined` — a generated rate
  equation references an undefined inactive-state rate constant. Independent of
  this proposal; worth a separate fix.
- **Allosteric over-generation.** ~98% of MWC mechanisms in the tested
  enumeration are non-identifiable. The loss filter would prune most of them, but
  it may also be worth not generating so many unfittable allosteric variants in
  the first place.
