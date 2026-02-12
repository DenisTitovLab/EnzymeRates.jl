# Rule: Profile and Optimize Compilation Latency (Julia)

Use this workflow whenever first-run latency is high, especially when warm runs are much faster.

## Goals
- Identify whether time is spent in compilation vs runtime.
- Find the concrete method families causing compile cost.
- Apply the lowest-risk fixes first, then architectural changes if needed.
- Verify improvements with cold and warm benchmarks.

## 1) Reproduce and Separate Cold vs Warm
1. Run the exact user workload in a fresh Julia process.
2. Run the same workload again in the same process.
3. Record:
   - wall time
   - allocation volume
   - `% compilation`

Example:
```julia
@time workload()
@time workload()
```

Interpretation:
- High first run + low second run + high `% compilation` => compilation latency problem.

## 2) Capture Compile Trace
Run with:
```bash
julia --project=. --trace-compile=/tmp/trace.jl -e '...workload...'
```
Then inspect:
- number of trace lines
- repeated package methods
- explosion in specializations for same function family

Useful shell queries:
```bash
wc -l /tmp/trace.jl
rg -o "YourPkg\.[A-Za-z0-9_!]+" /tmp/trace.jl | sort | uniq -c | sort -nr | head
rg -n "YourPkg\.(functionA|functionB|...)" /tmp/trace.jl | head -n 100
```

## 3) Find Common Root Causes
Prioritize these patterns:
- Type-encoded data structures where each instance creates a new concrete type.
- Hot functions specialized on rich type parameters (hashing, dedup, grouping, metadata access).
- Heavy `@generated` methods triggered across many candidate types.
- Eager calls to expensive generated paths for all candidates, even before pruning.

## 4) Optimization Strategy (Order Matters)

### A. Reduce specialization fanout (first choice)
- Add function barriers and `@nospecialize` in high-fanout utility paths.
- Rewrite utility methods to operate on runtime values via abstract interfaces.
- Avoid methods typed as `f(x::TypeWithHugeParams{...})` in broad loops.

### B. Delay expensive generated work
- Do not call generated methods for every candidate upfront.
- Add cheap structural estimators for grouping/filtering.
- Only materialize full generated outputs for finalists.

### C. Cache deterministic results
- Cache expensive derived metadata keyed by stable canonical keys.
- Cache failures too when failure is deterministic for a structure.

### D. Consider representation split (higher impact)
- Enumerate/search using runtime structs.
- Convert to type-heavy representation only at API boundary or on-demand.

## 5) Guardrails
- Preserve external API behavior unless explicitly changing API.
- Keep correctness tests first-class; performance changes are not enough.
- Do not rely on exhaustive precompile of combinatorial user inputs.

## 6) Validation Checklist
1. Targeted correctness tests pass.
2. Cold timing improves materially on user's exact command.
3. Warm timing does not regress significantly.
4. Candidate count/output equivalence remains unchanged.
5. Any known unrelated test failures are called out explicitly.

## 7) Reporting Format
Always report:
- baseline cold/warm metrics
- post-change cold/warm metrics
- key code paths changed (file + line)
- residual risks and next optimization options

## 8) Practical Tactics
- Prefer `rg` for fast code/trace search.
- Use direct script-style Julia invocations for reproducible measurements.
- Benchmark in fresh process for cold numbers.
- Re-run same process for warm numbers.

