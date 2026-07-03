# LDH memory-growth diagnostic — findings (preliminary)

Diagnoses the RAM climb (~10 GB → >50 GB by iteration 14) on the 2026-07-02 LDH
HPC run. Spec change 4 of `2026-07-03-ldh-search-improvements-design.md`.
**No fix is implemented here** — the fix is deferred pending an at-scale
confirmation run.

## Bottom line

- **Not type instability.** The growth is retained memory, not churn.
- **`GC.gc()` will not reduce RSS.** Measured directly: a full GC freed ~240 MB
  of heap (`gc_live` 326→85 MB) yet the process's *current* RSS did not move
  (952.3→952.4 MB). Freed memory is not returned to the OS, and compiled code is
  never freed.
- **Per-distinct-mechanism compilation is the plausible primary driver.**
  Running the same stack a real worker runs per mechanism (compile + rate
  equation + a fit), each distinct mechanism costs ~0.5 MB of non-returnable
  RSS on top of a ~267 MB one-time JIT of the optimizer stack. Spread across the
  run's tens of thousands of distinct mechanisms and a large worker pool, that
  reaches tens of GB — quantitatively consistent with the observed 50 GB.
- **The cumulative-floor fix already shipped in this branch mitigates it** by
  terminating the search far earlier (LDH: ~iteration 9 vs 22) and expanding far
  fewer mechanisms — fewer distinct compiles is directly less growth.

This **revises** the design doc's framing: the code/compilation hypothesis was
right. An earlier probe here wrongly dismissed it — it used `Sys.maxrss()`
(a monotonic high-water mark that cannot show reclaim) and skipped the fit
stack, undercounting per-mechanism cost ~10–25×. Both were corrected; the
numbers below are from the corrected probe.

## Method

This box is 7.7 GB — too small to reproduce a 50 GB climb, and each `@generated`
derivation for oligomeric LDH is seconds-slow, so a full-scale local run is out
of reach. Two experiments:

1. **Reduced full-pipeline run** (`docs/ldh_hpc_results/profile_memory.jl`,
   `max_param_count=7`, `min_beam_width=3`, 2 workers): per-process
   `gc_live_bytes`/`maxrss` sampled every 5 s, written incrementally.
2. **Single-process per-mechanism probe**
   (`docs/ldh_hpc_results/profile_codecache_probe.jl`): for each distinct
   mechanism, `compile_mechanism` + `rate_equation` + a tiny `fit_rate_equation`
   (the full stack a worker specializes per mechanism), tracking **current RSS
   (`/proc/self/status` `VmRSS`) — not `Sys.maxrss()`** — then a forced GC to
   measure reclaim.

## Evidence

**Per-mechanism probe** (40 distinct LDH init mechanisms, current RSS):

```
baseline:                gc_live= 48.7MB  rss= 658.4MB
after 10 compiled+fit:   gc_live=305.3MB  rss= 925.5MB    (+267 MB — one-time)
after 20 compiled+fit:   gc_live=190.7MB  rss= 942.9MB
after 40 compiled+fit:   gc_live=325.7MB  rss= 952.3MB    (~0.47 MB/mech, 20→40)
after forced GC:         gc_live= 85.1MB  rss= 952.4MB    (RSS unmoved)
```

- **GC does not shrink RSS.** The GC reclaimed ~240 MB of heap (`gc_live`
  326→85 MB) but current RSS held at 952 MB. So the retained ~294 MB over
  baseline is not returned to the OS — part compiled code (never freeable),
  part allocator-retained freed pages. Either way, `GC.gc()` is not a fix.
- **The unbounded term is compiled code.** The allocator-retained heap plateaus
  at the working set, but each *new* distinct mechanism adds fresh specialized
  code (derivation + rate equation + the optimizer's per-mechanism closure).
  Steady-state ~0.47 MB per distinct mechanism (the first ~10 also pay the
  ~267 MB one-time optimizer-stack JIT).

**Reduced full-pipeline run** — master `maxrss` flat at ~1168 MB after the first
sample, `gc_live` churning ~190–589 MB; a worker's `maxrss` climbs steadily to
~1187 MB by run's end (still rising, not plateaued). Too small to reach the
50 GB regime, but its steady worker climb is consistent with the per-mechanism
growth above. (This project resolves EnzymeRates to a non-dev build — visible in
its old progress-log format — irrelevant to the memory mechanics, which the
branch did not change.)

## Diagnosis

Per worker, RSS = baseline (~0.66 GB loaded Julia + EnzymeRates) + one-time
optimizer-stack JIT (~0.27 GB) + ~0.5 MB of non-returnable memory per distinct
mechanism it fits. The last term is unbounded in the number of distinct
mechanisms and is what makes RSS *climb* across iterations.

The run's 50 GB is the aggregate across the worker pool:
~50k distinct mechanisms / ~32 workers ≈ 1.5k per worker × ~0.5 MB ≈ 0.75 GB of
per-mechanism growth per worker, plus ~0.9 GB baseline+JIT, ≈ 1.6 GB/worker
× 32 ≈ 50 GB. Master-side retention (`frontier`/`memo`/`cv_pool` holding every
distinct mechanism with its full rendered equation string, plus the wide sparse
per-iteration DataFrame — iteration 9's CSV was 962 MB) adds to the master's
share but is a secondary term next to the worker aggregate.

The 50 GB estimate is inference — this box cannot reach that scale — but it is
now quantitatively consistent, not just qualitative.

## Recommendations (for the deferred fix)

1. **Confirm the split at scale.** Run `profile_memory.jl` on a full machine at
   the real worker count and `max_param_count=13`. It already samples every
   worker via `remotecall`, so the per-worker vs. master breakdown is directly
   readable — decide the fix from that split.
2. **Lean on the cumulative-floor fix first.** It sharply cuts the distinct-
   mechanism count and iteration count; re-measure peak RAM on the branch before
   adding anything else — it may bring the run under the RAM ceiling on its own.
3. **If workers dominate (expected):** recycle workers periodically (restart
   them to return memory to the OS and drop accumulated code cache), or run
   fewer workers. `GC.gc()` on workers reclaims only the working heap and does
   not shrink RSS, so it is not a substitute.
4. **If the master contributes materially:** stop storing the full equation
   string in every retained `BatchEntry.row` (needed only to write the CSV, then
   dropped), and write the per-iteration CSV without materializing the full wide
   `missing`-filled DataFrame in memory.

## Artifacts

- `docs/ldh_hpc_results/profile_memory.jl` — full-pipeline sampler (run at scale).
- `docs/ldh_hpc_results/profile_codecache_probe.jl` — the per-mechanism probe
  (current RSS + full fit stack).
- Raw samples from the reduced run are left untracked
  (`docs/ldh_hpc_results/memprofile_samples.csv`).
