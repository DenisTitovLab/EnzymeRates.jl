# Loss-parsimony beam filter — implementation design

**Date:** 2026-06-24
**Status:** Approved, not yet implemented.
**Companion:** `docs/design_notes/enumeration_redundancy_and_loss_dominance_pruning.md`
— the durable design note that characterizes the redundancy this filter prunes
and reports the complete-landscape validation. This spec covers only *how* to
build the filter.

## What this adds

A third loss knob, `loss_parsimony_threshold`, on `identify_rate_equation`. A
mechanism at parameter count `n` keeps expanding to count `n+1` only if its loss
stays within this factor of the best model with one fewer parameter. An added
parameter must earn its keep: if it cannot beat the best simpler model, the
beam stops growing it. The filter reuses losses already computed and needs no
sensitivity Jacobian.

The note proposed the keyword `loss_m1param_threshold` and an `Inf` (opt-in)
default. Two decisions override it: the keyword is **`loss_parsimony_threshold`**
and the default is **`1.01`** (active). The note's "AND onto all existing beam
criteria" is also superseded — see "The floor stays" below.

## The rule

Write `best(k)` for `best_loss_by_count[k]`, the lowest loss seen so far at
parameter count `k`. For a mechanism at count `n`, ranked by ascending loss:

```
keep i  ⇔  ( loss[i] ≤ loss_rel_threshold·best(n) + loss_abs_threshold
             AND loss[i] ≤ loss_parsimony_threshold·best(n-1) )
          OR  rank(i) ≤ min_beam_width
```

The two loss conditions are exactly `loss[i] ≤ min(loss_rel_threshold·best(n) +
loss_abs_threshold, loss_parsimony_threshold·best(n-1))`. The parsimony
threshold therefore only **tightens the existing cutoff**; everything else in
the beam logic stays unchanged.

## The floor stays

`min_beam_width` remains a guaranteed floor. The top `min_beam_width` mechanisms
by loss always expand, even when *zero* mechanisms clear the loss cutoff. This
is the deliberate override of the design note, which would have ANDed the
parsimony condition onto the whole selection (floor included) and let it shrink
the beam below `min_beam_width`.

Three consequences follow from keeping the floor:

- **The per-level best always survives.** Rank 1 always satisfies `rank ≤
  min_beam_width`, so the lowest-loss mechanism at every count always expands.
  The note's empirical "zero divergence between reached-best and true-best" is
  now a structural guarantee, not just an observation.
- **The filter is an expansion gate only.** A pruned mechanism is still fit,
  still written to the iteration CSV, and still eligible for the `cv_pool` and
  model selection. The filter chooses what to *grow*, nothing else. `cv_pool`
  is untouched.
- **The filter bites only past the floor.** With `min_beam_width = 50` it
  changes nothing for a tier of 50 or fewer candidates; it prunes the redundant
  tail on large reactions (`bi_bi`, `ter_bi`, hundreds per tier). This is the
  conservative trade the floor buys: it prunes where the cost lives while
  guarding against parameter epistasis on the large reactions the note could not
  validate exhaustively.

## Default value

`loss_parsimony_threshold = 1.01`. An added parameter must bring the loss within
1% of the best one-fewer-parameter model to keep expanding. The 1% slack absorbs
optimizer and floating-point noise, so the filter keeps a genuine tie. `Inf`
disables the filter.

## Exemptions

The parsimony term is dropped — leaving today's `loss_rel`/`loss_abs` behavior —
whenever `best(n-1)` is absent:

- **Base count.** The lowest parameter count has no `n-1` level. The init
  mechanisms are all identifiable, so nothing there needs pruning anyway.
- **Count gap.** If counts 3 and 5 exist but 4 does not, level 5 has no `best(4)`
  reference and falls through to the exempt branch. A missing reference only
  forgoes extra tightening, so the exemption is safe.

## Code changes

Four edits in `src/identify_rate_equation.jl`, each mirroring the existing
`loss_rel_threshold` / `loss_abs_threshold` plumbing.

1. **`_select_beam`** gains `parsimony_cutoff::Union{Nothing,Float64}=nothing`
   (the pre-multiplied cutoff, passed in like `best_override` is). One line
   after the existing `cutoff = …`:

   ```julia
   parsimony_cutoff !== nothing && (cutoff = min(cutoff, parsimony_cutoff))
   ```

   The rank/floor logic and the input-order return are unchanged.

2. **Call site** in `_beam_search` (the per-count `_select_beam` call): add

   ```julia
   parsimony_cutoff = haskey(best_loss_by_count, c - 1) ?
       loss_parsimony_threshold * best_loss_by_count[c - 1] : nothing,
   ```

   `best_loss_by_count` is already in scope at this call.

3. **Threading.** Add `loss_parsimony_threshold::Float64 = 1.01` to
   `identify_rate_equation`'s keyword list and to `_beam_search`'s signature,
   forwarded the same way the other two loss thresholds already are.

4. **Docstrings.** Update `identify_rate_equation`'s keyword list and "Beam
   selection" section, and `_select_beam`'s docstring, to describe the
   `min()` tightening and to state that `min_beam_width` stays a hard floor.

## Safety

The cutoff reads `best(n-1)` from a running minimum that only ever falls. At
selection time it therefore over-estimates the final `best(n-1)`, so the cutoff
is looser than a hindsight-perfect one. The filter can only ever be too
permissive, never too strict — it cannot wrongly make a level's true best
unreachable. The `min_beam_width` floor makes this airtight: the per-level best
expands regardless of the cutoff.

## Tests

Test-driven; written before the edits. Unit tests on `_select_beam`, matching
the existing `_select_beam best_override` testset:

1. **Floor guarantee.** With `parsimony_cutoff` below every loss and
   `min_beam_width = k`, exactly the top `k` by loss survive. This pins the
   explicit requirement that the floor outranks the loss filters.
2. **Tightening.** A `parsimony_cutoff` stricter than the rel/abs cutoff drops
   the mechanisms between the two cutoffs (those beyond the floor).
3. **No-op.** `parsimony_cutoff = nothing` reproduces current selection; the
   existing beam-selection testset stays green.
4. **Interaction.** Combined with `best_override`, `min()` picks the smaller of
   the two cutoffs.

One threading check: an end-to-end `identify_rate_equation` run on the existing
allosteric uni-uni fixture, passing an explicit `loss_parsimony_threshold`,
completes and selects a model — confirming that the keyword forwards through
`_beam_search` to the call site.
