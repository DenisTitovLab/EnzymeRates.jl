# LDH memory-growth diagnostic — findings (preliminary)

Diagnoses the RAM climb (~10 GB → >50 GB by iteration 14) on the 2026-07-02 LDH
HPC run — the memory item from the LDH search-improvements work.
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
  RSS on top of a ~267 MB one-time JIT of the derivation machinery and optimizer
  stack. Spread across the
  run's tens of thousands of distinct mechanisms and a large worker pool, that
  reaches tens of GB — quantitatively consistent with the observed 50 GB.
- **The cumulative-floor fix already shipped in this branch mitigates it** by
  terminating the search far earlier (LDH: ~iteration 9 vs 22) and expanding far
  fewer mechanisms — fewer distinct compiles is directly less growth.

This **revises** the design doc's framing: the code/compilation hypothesis was
right. An earlier probe here wrongly dismissed it — it used `Sys.maxrss()`
(a monotonic high-water mark that cannot show reclaim) and skipped the fit
stack, undercounting per-mechanism cost ~5–25×. Both were corrected; the
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
derivation + optimizer-stack JIT (~0.27 GB) + ~0.5 MB of non-returnable memory
per distinct mechanism it fits. The last term is unbounded in the number of
distinct mechanisms and is what makes RSS *climb* across iterations.

The run's 50 GB is the aggregate across the worker pool. Its iteration CSVs
total tens of thousands of distinct fitted equations; taking ~50k across a
~32-core node: ≈1.5k per worker × ~0.5 MB ≈ 0.75 GB of per-mechanism growth per
worker, plus ~0.9 GB baseline+JIT, ≈ 1.6 GB/worker × 32 ≈ 50 GB (the 50k and 32
are the rough inputs the estimate is most sensitive to — confirm at scale).
Master-side retention (`frontier`/`memo`/`cv_pool` holding every
distinct mechanism with its full rendered equation string, plus the wide sparse
per-iteration DataFrame — iteration 9's CSV was 962 MB) adds to the master's
share but is a secondary term next to the worker aggregate.

The 50 GB estimate is inference — this box cannot reach that scale — but it is
now quantitatively consistent, not just qualitative.

## Recommendations (for the deferred fix)

1. **Confirm the split at scale.** Run `profile_memory.jl` on a full machine at
   the real worker count and `max_param_count=13`. It already samples every
   worker via `remotecall`, so the per-worker vs. master breakdown is directly
   readable — decide the fix from that split. Add a same-mechanism control (fit
   one mechanism many times) to confirm the per-mechanism RSS growth tracks the
   distinct-type count rather than an allocator working-set ratchet.
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

## Artifacts — the probe scripts

Both scripts live in the gitignored `docs/ldh_hpc_results/` working area (which
holds the multi-GB result CSVs), so they are not tracked as files — they are
reproduced here so the repo retains them. Raw samples from the reduced run are
also left there untracked (`memprofile_samples.csv`).

### `profile_memory.jl` — full-pipeline sampler (run this at scale)

```julia
# ABOUTME: Local memory profiler for identify_rate_equation on reduced LDH data.
# ABOUTME: Samples master + per-worker gc_live_bytes / maxrss to a CSV (written
# ABOUTME: incrementally) to isolate code-cache growth from data retention.
using Pkg
Pkg.activate(joinpath(@__DIR__))
Pkg.instantiate()

using Distributed, Dates, CSV, DataFrames

addprocs(2; exeflags = ["--project"])   # 2 workers: box is memory-constrained
@everywhere using EnzymeRates, OptimizationCMAEvolutionStrategy

raw = CSV.read(joinpath(@__DIR__, "Enzyme data", "LDH_data.csv"), DataFrame)
filter!(row -> row.Rate != 0.0, raw)
data = (group = String.(raw.Article .* "_" .* raw.Fig),
        Rate = Float64.(raw.Rate), NADH = Float64.(raw.NADH),
        Pyruvate = Float64.(raw.Pyruvate), Lactate = Float64.(raw.Lactate),
        NAD = Float64.(raw.NAD))
rxn = @enzyme_reaction begin
    substrates:NADH[C21H29N7O14P2], Pyruvate[C3H4O3]
    products:Lactate[C3H6O3], NAD[C21H27N7O14P2]
    oligomeric_state:4
end
prob = IdentifyRateEquationProblem(rxn, data; Keq=20000.0, scale_k_to_kcat=1.0)

# Sampler: master + each worker, every 5 s, appended to a CSV as it goes so a
# partial or OOM-killed run still yields the growth trend.
@everywhere _mem() = (gc_live = Base.gc_live_bytes(), maxrss = Sys.maxrss())
samples_path = joinpath(@__DIR__, "memprofile_samples.csv")
io = open(samples_path, "w")
println(io, "t,who,gc_live,maxrss"); flush(io)
t0 = time()
sampling = Ref(true)
sample_row(t, who, m) = (println(io, "$t,$who,$(m.gc_live),$(m.maxrss)"); flush(io))
sampler = @async while sampling[]
    t = round(time() - t0; digits=1)
    sample_row(t, "master", _mem())
    for w in workers()
        sample_row(t, "worker$w", remotecall_fetch(_mem, w))
    end
    sleep(5)
end

# Reduced + FAST run. Memory growth tracks the count of distinct @generated
# mechanisms + iterations, not fit quality. For an at-scale profile, raise
# max_param_count to 13 and use the real worker count.
results = identify_rate_equation(prob;
    optimizer=CMAEvolutionStrategyOpt(),
    max_param_count=7, min_beam_width=3,
    n_restarts=2, maxtime=2.0,
    loss_rel_threshold=1.3, loss_abs_threshold=0.001,
    loss_parsimony_threshold=1.01,
    save_dir=joinpath(@__DIR__, "memprofile_results"))

sampling[] = false; wait(sampler)

# Does a forced GC reclaim? NOTE maxrss is a high-water mark; for a valid reclaim
# test read current RSS (see the probe below). This sampler is for the growth
# TREND and the per-worker vs. master split.
GC.gc(true); GC.gc(true)
sample_row(round(time() - t0; digits=1), "master_postGC", _mem())
for w in workers()
    sample_row(round(time() - t0; digits=1), "worker$(w)_postGC",
               remotecall_fetch(() -> (GC.gc(true); _mem()), w))
end
close(io)

allsamples = CSV.read(samples_path, DataFrame)
mb(x) = round(x / 2^20; digits=1)
println("=== peak maxrss (MB) ===")
for who in unique(allsamples.who)
    ws = allsamples[allsamples.who .== who, :]
    println("  $who: first=", mb(first(ws.maxrss)), " peak=", mb(maximum(ws.maxrss)),
            " gc_live first=", mb(first(ws.gc_live)), " peak=", mb(maximum(ws.gc_live)))
end
println("Selected: n_params=", results.cv_results.n_params[1])
for p in workers(); rmprocs(p); end
```

### `profile_codecache_probe.jl` — per-mechanism probe (current RSS + full fit stack)

Run against the LDH project (it has the optimizer):
`julia --project=docs/ldh_hpc_results docs/ldh_hpc_results/profile_codecache_probe.jl`.

```julia
# ABOUTME: Single-process probe for the LDH RAM growth. Per distinct mechanism it
# ABOUTME: runs the SAME stack a real worker compiles (compile + rate_equation +
# ABOUTME: a tiny fit), tracks CURRENT RSS (VmRSS), then GCs and measures reclaim.
# Run against the LDH project (has the optimizer): --project=docs/ldh_hpc_results.
# NOTE: uses current RSS from /proc/self/status, NOT Sys.maxrss() — maxrss is a
# monotonic high-water mark and cannot show reclaim.
using Pkg
Pkg.activate(@__DIR__)
using EnzymeRates, OptimizationCMAEvolutionStrategy, Printf

function current_rss()
    for line in eachline("/proc/self/status")
        startswith(line, "VmRSS:") && return parse(Int, split(line)[2]) * 1024
    end
    return 0
end
mb(x) = round(x / 2^20; digits=1)
sample(tag) = @printf("%-26s gc_live=%7.1fMB  rss(current)=%7.1fMB\n",
                      tag, mb(Base.gc_live_bytes()), mb(current_rss()))

rxn = @enzyme_reaction begin
    substrates:NADH[C21H29N7O14P2], Pyruvate[C3H4O3]
    products:Lactate[C3H6O3], NAD[C21H27N7O14P2]
    oligomeric_state:4
end
inits = unique!(collect(EnzymeRates.init_mechanisms(rxn)))
# Tiny synthetic dataset so a real fit runs (fit quality is irrelevant — we only
# want the per-mechanism compilation the optimizer stack does).
data = (group = ["a","a","b","b"], Rate = [1.0, 2.0, 1.5, 2.5],
        NADH = [1.0, 2.0, 1.0, 2.0], Pyruvate = [1.0, 1.0, 2.0, 2.0],
        Lactate = [0.5, 0.5, 0.5, 0.5], NAD = [0.5, 0.5, 0.5, 0.5])
opt = CMAEvolutionStrategyOpt()
function compile_and_fit(m)
    try
        em = EnzymeRates.compile_mechanism(m)
        fp = FittingProblem(em, data; Keq = 20000.0, scale_k_to_kcat = 1.0)
        fit_rate_equation(fp, opt; n_restarts = 1, maxtime = 0.3, maxiters = 100)
    catch
    end
end

n = min(40, length(inits))
GC.gc(true); GC.gc(true)
println("distinct init mechanisms: ", length(inits), " (probing ", n, ")")
sample("baseline:")
for (i, m) in enumerate(inits[1:n])
    compile_and_fit(m)
    (i % 10 == 0 || i == n) && sample("after $i compiled+fit:")
end

# Valid reclaim test: with CURRENT RSS, if it FALLS toward baseline after a full
# GC the memory was reclaimable (GC.gc helps); if it stays high the retained
# memory is non-reclaimable (compiled code) or allocator-retained.
GC.gc(true); GC.gc(true); GC.gc(true)
sample("after forced GC:")
```
