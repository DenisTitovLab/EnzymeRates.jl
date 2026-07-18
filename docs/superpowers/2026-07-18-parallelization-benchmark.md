# Seed-enumeration parallelization benchmark

Date: 2026-07-18
Branch: `parallelize-seed-and-expansion`

Confirms the two claims behind parallelizing `seed_mechanisms`: the wave-parallel
BFS produces byte-identical seeds across real worker processes, and it runs
faster.

## Setup

- Machine: 4 physical cores.
- Reaction: PFKP (`docs/hpc_results/pfkp_hpc_results/identify_pfkp.jl`),
  restricted to three required allosteric regulators (`ATP`, `ADP`,
  `Phosphate`) so the serial baseline stays tractable. The full five-regulator
  case is ~20 minutes serial.
- Both runs are JIT-warm (a one-regulator warm-up precedes each timed run, and
  the workers run one warm-up wave before the parallel timing).

## Result

| Run | Seeds | Wall time |
|---|---|---|
| Serial (no workers) | 1104 | 229.3 s |
| Parallel (3 workers) | 1104 | 86.7 s |

- **Seeds identical across workers (content and order): true.** This is the
  load-bearing check — it proves `Mechanism` objects serialize round-trip to the
  workers and back, and that the wave-parallel BFS is byte-identical in
  production, not only under `pmap`-on-the-main-process (the test suite's case).
- **Speedup: 2.64× on 3 workers** — near-linear, bounded by the serial
  `visited`-dedup barrier between waves plus worker dispatch overhead.

Three workers give 2.64×; the closure explored is exponential in the number of
required regulators, so the same wave-parallel structure across an HPC node's
hundreds of cores turns the observed ~30-minute five-regulator stall into
seconds. `_expand_parents` parallelizes the per-iteration expansion with the
identical `pmap`-over-mechanisms pattern, so the same serialization guarantee
this benchmark confirms covers it.
