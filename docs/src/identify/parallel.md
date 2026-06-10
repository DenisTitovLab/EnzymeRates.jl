# Running in parallel

The equation search is embarrassingly parallel. Every candidate equation is
compiled and fit independently of the others, and every leave-one-group-out
cross-validation fold is independent too. [`identify_rate_equation`](@ref)
farms both out through its `pmap_function` argument, which defaults to
`Distributed.pmap`. The same call therefore scales from one core to a whole
high-performance compute (HPC) cluster — you change nothing in the call, you
only add worker processes.

Because a realistic search fits thousands of candidates with a multi-start
optimizer, this is the difference between a run that finishes overnight on a
cluster and one that does not finish on a laptop.

## One machine, many cores

Add local worker processes with `addprocs`, load the package and your optimizer
on every worker with `@everywhere`, then call `identify_rate_equation` as
usual — `pmap` distributes the work across the workers automatically.

```julia
using Distributed
addprocs(8)                                        # 8 worker processes

@everywhere using EnzymeRates, OptimizationPyCMA   # load on every worker

results = identify_rate_equation(prob;
    optimizer = PyCMAOpt(),
    # pmap_function = pmap is the default: candidate fits and CV folds
    # are now spread across the 8 workers.
)
```

The `@everywhere` line is required: each worker compiles and fits mechanisms on
its own, so it needs both `EnzymeRates` and the optimizer package loaded. The
main process collects the results and writes the CSV outputs to `save_dir`;
workers do not write CSVs themselves.

## A cluster: many machines

On a cluster, use a cluster manager to start workers across the allocated nodes,
then keep the call identical. For Slurm, `SlurmClusterManager.jl` reads the
allocation and starts one worker per task:

```julia
using Distributed, SlurmClusterManager
addprocs(SlurmManager())                           # one worker per Slurm task

@everywhere using EnzymeRates, OptimizationPyCMA

results = identify_rate_equation(prob; optimizer = PyCMAOpt())
```

Run this from inside a Slurm allocation (`salloc`, or a script submitted with
`sbatch`); the number of workers equals the number of allocated tasks. Other
schedulers are reached the same way through `ClusterManagers.jl` (PBS, SGE,
LSF) or `MPIClusterManagers.jl` (MPI). In every case the `identify_rate_equation`
call is unchanged — only the `addprocs` line differs.

## Forcing serial execution

Pass `pmap_function = map` to run everything on the main process with no
workers. This is useful for debugging, for small problems where worker startup
is not worth it, and for the fast example in the [Identify tutorial](@ref).

```julia
results = identify_rate_equation(prob;
    optimizer = PyCMAOpt(),
    pmap_function = map,    # single process, no Distributed workers
)
```
