# Parallelize seed enumeration and expansion

Date: 2026-07-18
Status: Approved, ready for planning
Branch: `parallelize-seed-and-expansion`

## Problem

The beam search runs two steps single-threaded on the main node while every
worker sits idle:

1. **Seed enumeration** (`seed_mechanisms`, `mechanism_enumeration.jl:2289`).
   The July 2026 PFKP run spent roughly 30 minutes here before any fitting
   began. LDH showed no such delay.
2. **Per-iteration expansion** (`for pe in to_expand`, `identify_rate_equation.jl:832`).
   Each iteration pauses on one core to expand every selected parent into its
   children before the next round of fitting starts.

Profiling (single-threaded, matching one HPC core) confirmed both costs:

| Seed enumeration, by required allosteric regulators | 1 | 2 | 3 | 4 | 5 (PFKP) |
|---|---|---|---|---|---|
| wall time | 49 s | 104 s | 238 s | 492 s | ~17–30 min (projected) |
| seed count | 1104 | 1104 | 1104 | 1104 | 1104 |

The seed count stays flat while the time grows ~2.1× per required regulator:
the BFS closure explores an exponentially growing graph of intermediate nodes
and keeps only the 1104 that bind every required regulator. PFKP requires five
allosteric regulators (`#67`, v0.2.0), so it pays the full cost; LDH marks all
inhibitors optional and takes the cheap `init_mechanisms` path.

Expansion measured ~2.7 s per parent. A real iteration expands 700–1000
parents, so one core spends 30–45 minutes per iteration while the workers wait.
Writing the iteration CSV adds ~33 s (586 MB at ~18 MB/s).

Fitting already runs across all workers (`pmap`, two passes in `_process_batch`).
Only these two enumeration steps remain serial.

## Goal

Move both steps onto the workers. Keep the results identical — this spec changes
*where* work runs, never *what* the search enumerates or fits. Reducing the
enumeration itself (pruning the seed closure, restricting expansion moves) is a
separate effort tracked under the "fit fewer mechanisms" theme.

## Design

### Correctness principle

`pmap` returns results in input order: `result[i]` holds the output for
`input[i]`, whatever worker computed it. Both rewrites keep every dedup and
ordering decision on the main node, driven by that preserved order, so each
produces byte-identical output to the serial code it replaces. When no worker
processes exist — the test suite's case — `pmap` runs on the main process, so
the parallel path stays exercised and correct there too.

### Seed enumeration: wave-parallel BFS

The serial BFS pops one node at a time from a FIFO queue, generates its
children, filters them, and enqueues the new ones. Because a FIFO queue drains
level by level, the enqueue order equals the level order.

Process one level at a time instead:

```
frontier = deduplicated init_mechanisms(rxn)
while frontier is non-empty
    childsets = pmap(frontier) do node          # parallel, pure, independent
        filter(child -> _is_seed_node(child, rxn, required_allo, required_comp),
               _seed_children(node, rxn, required_allo))
    end
    frontier = []                                # serial barrier: dedup + collect
    for childset in childsets, child in childset # frontier order = enqueue order
        if consider!(child)                      # new? mark visited, collect if a seed
            push!(frontier, child)
        end
    end
end
```

`pmap` distributes the expensive per-node work — `_seed_children` runs three
expansion moves, each constructing and canonicalizing mechanisms. The `visited`
set and seed collection stay serial and cheap. Iterating `childsets` in frontier
order visits children in the same order the serial BFS enqueued them, so the
`visited` set, the seed vector, and its order all match the current output.

### Expansion: parallel flat-map

Expand every selected parent on a worker, then run the existing serial merge
over the order-preserved results:

```
results = pmap(pe -> _expand_parent(pe.mech, prob.reaction), to_expand)
for (pe, (kids, failure)) in zip(to_expand, results)
    failure === nothing || push!(expand_failures, failure)
    for child in kids
        haskey(parent_of, child) && continue     # first parent wins, unchanged
        parent_of[child] = (mechanism_type = pe.row.mechanism_type,
                            n_params = pe.n_params)
        push!(children, child)
    end
end
```

`_expand_parent` already catches a per-parent expansion error and returns it as
a `FitFailure` carrying the exception text, so failures serialize back and land
in the same bucket as before. Iterating `zip(to_expand, results)` preserves the
first-parent-wins dedup, the child order, and the failure order.

Send only what the worker needs. Mapping over `to_expand` ships each
`BatchEntry`, including its row's long `mechanism_type` string; mapping over the
mechanisms alone (`pmap(m -> _expand_parent(m, rxn), [pe.mech for pe in to_expand])`)
keeps `parent_of`'s bookkeeping on the main node, where it already lives. The
same applies to the seed frontier — ship mechanisms, not richer wrappers.

## Scope

Out of scope, by decision:

- The CSV write (~33 s). Secondary, and streaming it is a separate change.
- The seed closure size and the expansion move set. Shrinking either changes
  what the search enumerates — the "fit fewer mechanisms" theme, not this spec.
- `rate_equation`. Untouched, so its 0-allocation / sub-120 ns gate is
  unaffected.

## Testing

1. **Equivalence.** An inline serial reference BFS and the wave-parallel
   `seed_mechanisms` return the same seeds, in the same order, for a reaction
   whose closure runs in well under a second — a uni-uni or bi-uni with one or
   two required regulators, not the slow PFKP case. The logic is size-independent,
   so a small closure proves equivalence; the large case belongs in the speedup
   benchmark below. The parallel flat-map and the serial expansion loop return the
   same children and the same `parent_of` map for a set of parents.
2. **Invariants.** Every seed binds all required regulators; the seed vector has
   no duplicates; the order is deterministic across runs.
3. **Regression.** The full suite stays green.
4. **Speedup.** A local `addprocs` benchmark confirms the wall-time drop — the
   empirical check belongs in implementation, not the suite.

## Roadmap

This is the first of a sequenced set the profiling motivated:

1. **This spec** — parallelize seed enumeration and expansion (ideas 1, 6).
2. Log the fit plan before fitting, results after (idea 2).
3. Fit fewer mechanisms — strict-sequential parameter advance, limit the
   dominant low-yield moves (`to_allosteric` for LDH at 47% of children,
   `merge_reg_sites` for PFKP at 64%), and recognize near-inseparable equations
   (ideas 3, 4, 5). Idea 3 counts only children of parents that a strict advance
   would never have expanded, not every non-parsimonious fit.
