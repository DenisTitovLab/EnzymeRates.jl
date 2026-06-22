# Running in parallel

A reaction as ordinary as two substrates and two products has more than 100,000
biochemically valid rate equations: the order substrates bind and products
leave, the rapid-equilibrium-versus-steady-state choice at each step, and the
dead-end complexes that can form all combine into a vast space.
[`identify_rate_equation`](@ref) makes this tractable two ways. The
rate-equation search keeps only the promising candidates at each parameter count
instead of fitting the whole space, and the candidates it does fit run in
parallel. This page
covers the parallelism.

The search is embarrassingly parallel: every candidate equation is compiled and
fit independently, and every leave-one-group-out cross-validation fold is
independent too. `identify_rate_equation` runs both through `Distributed.pmap`.
The same call scales from one core to a whole high-performance compute (HPC)
cluster — you add worker processes and change nothing else. With no workers,
`pmap` runs every fit on the main process, so the call already works on a laptop.

## Running on a local machine with several cores

Add local workers with `addprocs`, load the package and your optimizer on every
worker with `@everywhere`, then load the data and run the search on the main
process:

```julia
using Distributed, CSV
addprocs(8)                                  # 8 local worker processes

@everywhere using EnzymeRates, OptimizationCMAEvolutionStrategy

# Read the data and build the problem on the main process. The workers receive
# the problem automatically when the search distributes its fits to them.
table = CSV.File("rate_data.csv")
data = (group = table.group, Rate = table.Rate,
        A = table.A, B = table.B, P = table.P, Q = table.Q)   # one column per metabolite

rxn = @enzyme_reaction begin
    substrates: A[C], B[C]
    products:   P[C], Q[C]
end
prob = IdentifyRateEquationProblem(rxn, data; Keq = 5.0)

results = identify_rate_equation(prob; optimizer = CMAEvolutionStrategyOpt())

# Clean up the workers when the run finishes
for p in workers()
    rmprocs(p)
end
```

The `@everywhere using` line is required: each worker compiles and fits
mechanisms on its own, so it needs both `EnzymeRates` and the optimizer package.
The data file is read once on the main process and converted to a plain
`NamedTuple`; `pmap` ships that problem to the workers with each batch of fits,
so the workers need no data-loading packages. The main process collects the
results and writes the CSV outputs to `save_dir`; the workers write nothing.

## Running on a compute cluster with many machines

On a cluster, swap the `addprocs` line for a cluster manager that starts workers
across the allocated nodes — the rest of the script is identical. For Slurm,
`SlurmClusterManager.jl` reads the allocation and starts one worker per task.
Save the script as `identify.jl`:

```julia
using Distributed, SlurmClusterManager, CSV
addprocs(SlurmManager())                     # one worker per Slurm task

@everywhere using EnzymeRates, OptimizationCMAEvolutionStrategy

table = CSV.File("rate_data.csv")
data = (group = table.group, Rate = table.Rate,
        A = table.A, B = table.B, P = table.P, Q = table.Q)

rxn = @enzyme_reaction begin
    substrates: A[C], B[C]
    products:   P[C], Q[C]
end
prob = IdentifyRateEquationProblem(rxn, data; Keq = 5.0)

results = identify_rate_equation(prob; optimizer = CMAEvolutionStrategyOpt())

# Clean up the workers when the run finishes
for p in workers()
    rmprocs(p)
end
```

Submit the script with `sbatch`: the batch script requests the allocation and
launches Julia, and `SlurmManager()` then starts one worker per allocated task.
This template uses [`juliaup`](https://github.com/JuliaLang/juliaup) to select
the Julia version:

```bash
#!/bin/bash
# Job name:
#SBATCH --job-name=YOUR_JOB_NAME
#
# Account:
#SBATCH --account=YOUR_HPC_ACCOUNT
#
# Partition:
#SBATCH --partition=YOUR_HPC_PARTITION_NAME
#
# Number of requested nodes:
#SBATCH --nodes=24
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=24:00:00

# Load software
module purge

# Run Julia script
/PATH/TO/USER/.juliaup/bin/julia YOUR_IDENTIFY_RATE_EQUATION_CODE.jl
```

Save this batch script with a `.slurm` extension and submit it with `sbatch
identify.slurm`. Replace the placeholders with your account, partition, and
paths, and point the last line at your Julia script — run Julia in the project
where `EnzymeRates` and your optimizer are installed (add
`--project=/path/to/project` if it is not the default). The worker count follows
the allocation, so widening the search means requesting more nodes. Other
schedulers work the same way through `ClusterManagers.jl` (PBS, SGE, LSF) or
`MPIClusterManagers.jl` (MPI); only the `addprocs` line changes.
