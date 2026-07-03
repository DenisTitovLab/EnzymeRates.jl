# LDH memory-growth diagnostic — findings (preliminary)

Diagnoses the RAM climb (~10 GB → >50 GB by iteration 14) seen on the
2026-07-02 LDH HPC run. Spec change 4 of
`2026-07-03-ldh-search-improvements-design.md`. **No fix is implemented here**
— the fix is deferred pending an at-scale confirmation run.

## Bottom line

- **Not type instability.** The growth is retained memory, not churn.
- **`GC.gc()` will not fix it.** The two candidate stores are both immune:
  compiled code lives outside the GC heap, and the search's retained data is
  *live by design* (not garbage).
- **The `@generated`-per-mechanism code cache is real but small per process**
  — a one-time ~110 MB of infrastructure, then only ~20–90 KB per additional
  distinct mechanism. That caps it at single-digit GB even at the run's scale,
  so it is **not** the 50 GB driver on its own. (This refutes the design doc's
  leading hypothesis.)
- **The 50 GB is an aggregate**, most plausibly dominated by **per-worker
  baselines × many workers** (a fresh worker is ~0.65 GB of loaded Julia +
  EnzymeRates before it fits anything; 32–40 of them is 20–26 GB at startup,
  matching "~10 GB at start" for a smaller pool) **plus master-side data
  accumulation** (`frontier`/`memo`/`cv_pool` retain every distinct mechanism,
  each carrying its full rendered equation string, across all iterations; plus
  the wide sparse per-iteration DataFrame).
- **The cumulative-floor fix already shipped in this branch mitigates it** by
  terminating the search far earlier (LDH: ~iteration 9 vs 22) and expanding
  far fewer mechanisms — fewer distinct compiles *and* less retained data.

## Method

This box is 7.7 GB — too small to reproduce a 50 GB climb, and each `@generated`
derivation for oligomeric LDH is seconds-slow, so a full-scale local run is out
of reach. Two experiments instead:

1. **Reduced full-pipeline run** (`docs/ldh_hpc_results/profile_memory.jl`,
   `max_param_count=7`, `min_beam_width=3`, 2 workers): per-process
   `gc_live_bytes`/`maxrss` sampled every 5 s to a CSV written incrementally.
2. **Single-process code-cache probe**
   (`docs/ldh_hpc_results/profile_codecache_probe.jl`): compile + evaluate each
   distinct mechanism (69 init, then 30 diverse expanded children), forcing the
   `@generated` specialization, sampling `maxrss`; then force GC and measure
   reclaim. Run against the branch source (`--project=<repo>`).

## Evidence

**Reduced full-pipeline run** — master `maxrss` flat at **1168 MB** for the whole
~9 min while `gc_live` churned 207–532 MB; workers ~965 MB after ~123 compiled
mechanisms. At `max_param_count=7` the search is too small to accumulate, so it
plateaus — it neither reproduces nor rules out the at-scale growth. (It did
confirm the shipped progress-log format is absent here — this project resolves
EnzymeRates to a non-dev build; irrelevant to the memory mechanics, which the
branch did not change.)

**Code-cache probe** (the decisive one):

```
baseline (0 compiled):   gc_live= 47.6MB  maxrss= 654.7MB
init 20:                 gc_live=195.2MB  maxrss= 760.9MB     (+106 MB)
init 69:                 gc_live=164.5MB  maxrss= 765.2MB     (+4 MB over 20→69)
+30 diverse children:    gc_live= 58.5MB  maxrss= 765.9MB     (+0.7 MB)
after forced GC:         gc_live= 54.3MB  maxrss= 765.9MB     (maxrss unmoved)
```

Two facts:
- **Non-reclaimable:** a full GC returns `gc_live` to its ~50 MB baseline but
  leaves `maxrss` at 766 MB. The retained ~110 MB is compiled code
  (MethodInstances / CodeInstances), which the GC never frees.
- **Small and front-loaded:** the bulk is the first ~20 compiles (the symbolic
  derivation machinery + LLVM/type infrastructure). Each *additional* distinct
  mechanism — including structurally diverse expanded children — adds only
  ~20–90 KB. Extrapolated to the run's tens of thousands of distinct equations
  that is single-digit GB per process, not 50 GB.

## Diagnosis

Per process, memory = baseline (~0.65 GB) + one-time code infrastructure
(~0.11 GB) + small per-mechanism code cache (~KB each, GC-immune) + GC-able
working heap. None of these makes one process reach 50 GB.

The run's 50 GB is therefore an **aggregate across the worker pool plus the
master**. The two terms that scale are:

1. **Worker count.** Every worker pays the full baseline + infrastructure, plus
   its own code cache for the mechanisms it fits. 32–40 workers turns ~1 GB/worker
   into 30–40 GB with no single-process leak in sight. This is the most likely
   dominant term and is fully consistent with "10 GB at start, 50 GB later."
2. **Master-side retention.** `frontier` (all unexpanded structurally-distinct
   mechanisms), `memo` (every distinct equation's fit, kept for the whole run),
   and `cv_pool` accumulate across iterations, and each retained `BatchEntry.row`
   holds the **full rendered equation string** (multi-KB for LDH). The wide
   sparse per-iteration DataFrame (`_rows_to_dataframe`, one column per distinct
   param name, `missing`-filled — iteration 9's CSV was 962 MB) is a transient
   multi-GB spike per iteration.

Neither term yields to `GC.gc()`: worker code cache is GC-immune, and the
master's retained mechanisms are live references, not garbage.

## Recommendations (for the deferred fix)

1. **Confirm the split at scale.** Run `profile_memory.jl` on a full machine
   with the real worker count and `max_param_count=13`. It already samples every
   worker via `remotecall`, so the per-worker vs. master breakdown is directly
   readable. Decide the fix from that split, not from this constrained box.
2. **Lean on the cumulative-floor fix.** It already cuts the distinct-mechanism
   count and iteration count sharply — re-measure peak RAM on the branch before
   adding anything else; it may suffice.
3. **If workers dominate:** cap or recycle the worker pool (restart workers
   periodically to flush accumulated code cache), or trade breadth for fewer
   workers. `GC.gc()` on workers reclaims only the working heap, not the cache.
4. **If the master dominates:** stop storing the full equation string in every
   retained `BatchEntry.row` (it is only needed to write the CSV, then dropped)
   and write the per-iteration CSV without materializing the full wide
   `missing`-filled DataFrame in memory.

## Artifacts

- `docs/ldh_hpc_results/profile_memory.jl` — full-pipeline sampler (run at scale).
- `docs/ldh_hpc_results/profile_codecache_probe.jl` — the decisive code-cache probe.
- Raw samples from the reduced run are left untracked
  (`docs/ldh_hpc_results/memprofile_samples.csv`).
