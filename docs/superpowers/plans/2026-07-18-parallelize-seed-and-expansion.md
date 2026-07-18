# Parallelize seed enumeration and expansion — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the two serial main-node steps of the beam search — seed
enumeration and per-iteration expansion — onto the workers via `pmap`, keeping
results byte-identical.

**Architecture:** Both steps run independent per-item work. `seed_mechanisms`
becomes a wave-parallel BFS: each level's per-node child generation runs on
workers, the `visited` dedup and seed collection stay serial on the main node.
The expansion loop moves into a new `_expand_parents` helper that maps
`_expand_parent` over workers, then runs the existing serial merge. `pmap`
preserves input order, so both keep every dedup and ordering decision on the
main node and produce identical output to the serial code.

**Tech Stack:** Julia, `Distributed` (`pmap`, already package-wide via
`EnzymeRates.jl:28`), the existing test suite (`Pkg.test()`).

## Global Constraints

- Results stay byte-identical. This changes where work runs, never what the
  search enumerates or fits.
- `rate_equation` is untouched. Its 0-allocation / sub-120 ns gate must stay
  green.
- `pmap` runs on the main process when no workers exist (the test suite's case),
  so the parallel path is exercised and correct there.
- Both rewrites are behavior-preserving refactors. Each is guarded by a
  characterization test that compares the new code to an inline serial reference
  and must stay green across the change.
- 92-character line limit, 4-space indent. Match surrounding style.

---

### Task 1: Wave-parallel `seed_mechanisms`

**Files:**
- Modify: `src/mechanism_enumeration.jl:2289-2311` (the `seed_mechanisms` body)
- Test: `test/test_mechanism_enumeration.jl`

**Interfaces:**
- Consumes: `init_mechanisms(rxn)`, `_seed_children(m, rxn, required_allo)`,
  `_is_seed_node(c, rxn, required_allo, required_comp)`,
  `_binds_all_required(m, required_allo, required_comp)` — all unchanged.
- Produces: `seed_mechanisms(rxn::EnzymeReaction, required_allo::Set{Symbol},
  required_comp::Set{Symbol}) -> Vector{Union{Mechanism, AllostericMechanism}}`
  — same signature, same output vector as before.

- [ ] **Step 1: Write the characterization test**

Add to `test/test_mechanism_enumeration.jl` (the `uni_uni_allo_reg` reaction —
substrates S, products P, one allosteric regulator R, oligomeric 2 — is defined
near the top of the file):

```julia
@testset "seed_mechanisms wave-parallel equivalence" begin
    # Inline serial FIFO BFS = the reference the parallel version must match.
    function serial_seed_reference(rxn, req_allo, req_comp)
        visited = Set{UInt64}()
        queue = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
        seeds = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
        enq(m) = begin
            h = hash(m)
            h in visited && return
            push!(visited, h); push!(queue, m)
            EnzymeRates._binds_all_required(m, req_allo, req_comp) && push!(seeds, m)
        end
        for m in EnzymeRates.init_mechanisms(rxn); enq(m); end
        while !isempty(queue)
            m = popfirst!(queue)
            for c in EnzymeRates._seed_children(m, rxn, req_allo)
                EnzymeRates._is_seed_node(c, rxn, req_allo, req_comp) && enq(c)
            end
        end
        seeds
    end

    req = Set([:R])
    empty = Set{Symbol}()
    got = EnzymeRates.seed_mechanisms(uni_uni_allo_reg, req, empty)
    ref = serial_seed_reference(uni_uni_allo_reg, req, empty)

    @test got == ref                                   # same seeds, same order
    @test !isempty(got)                                # the case is non-trivial
    @test allunique(hash.(got))                        # no duplicate structures
    @test all(m -> EnzymeRates._binds_all_required(m, req, empty), got)
    @test EnzymeRates.seed_mechanisms(uni_uni_allo_reg, req, empty) == got  # deterministic
end
```

- [ ] **Step 2: Run the test against the current serial code**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: the whole suite passes, including the new
`seed_mechanisms wave-parallel equivalence` testset. The current
`seed_mechanisms` is the serial FIFO BFS, so it matches the reference. This pins
current behavior before the refactor. (For faster iteration, run just this file
after loading fixtures: `julia --project -e 'using EnzymeRates, Test, Random,
LinearAlgebra; include("test/mechanism_definitions_for_test_enzyme_derivation.jl");
include("test/test_mechanism_enumeration.jl")'`.)

- [ ] **Step 3: Refactor `seed_mechanisms` to a wave-parallel BFS**

Replace the body at `src/mechanism_enumeration.jl:2289-2311` with:

```julia
function seed_mechanisms(rxn::EnzymeReaction, required_allo::Set{Symbol},
                         required_comp::Set{Symbol})
    visited = Set{UInt64}()
    seeds = Union{Mechanism, AllostericMechanism}[]
    # Dedup + collect on the main node. Returns true when `m` is new, so the
    # caller advances only genuinely-new nodes to the next wave. Called in
    # frontier order, which equals the serial BFS enqueue order, so `visited`
    # and `seeds` end byte-identical to the FIFO version.
    consider!(m) = begin
        h = hash(m)
        h in visited && return false
        push!(visited, h)
        _binds_all_required(m, required_allo, required_comp) && push!(seeds, m)
        true
    end
    frontier = Union{Mechanism, AllostericMechanism}[
        m for m in init_mechanisms(rxn) if consider!(m)]
    while !isempty(frontier)
        # Per-node child generation is pure and independent — distribute it.
        childsets = pmap(frontier) do m
            filter(c -> _is_seed_node(c, rxn, required_allo, required_comp),
                   _seed_children(m, rxn, required_allo))
        end
        next = Union{Mechanism, AllostericMechanism}[]
        for cs in childsets, c in cs
            consider!(c) && push!(next, c)
        end
        frontier = next
    end
    seeds
end
```

Update the docstring's closing line ("Returned deduped (each node is visited
once).") to note the wave-parallel structure if it now reads as stale; keep the
rest of the docstring, which still describes the seed set correctly.

- [ ] **Step 4: Run the test to confirm the refactor preserved behavior**

Run: the same command as Step 2.
Expected: PASS. `got == ref` proves the wave-parallel output matches the serial
FIFO output exactly.

- [ ] **Step 5: Run the full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: all green (including the `rate_equation` performance gate, untouched).

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Parallelize seed_mechanisms as a wave-parallel BFS"
```

---

### Task 2: Parallel `_expand_parents`

**Files:**
- Modify: `src/identify_rate_equation.jl:820-846` (extract + parallelize the
  expansion loop in `_beam_search`)
- Test: `test/test_identify_rate_equation.jl`

**Interfaces:**
- Consumes: `_expand_parent(m, reaction) -> (children::Vector, failure)`
  (`identify_rate_equation.jl:694`, unchanged); `BatchEntry` fields `.mech`,
  `.n_params`, `.row.mechanism_type`.
- Produces: `_expand_parents(to_expand::Vector{BatchEntry},
  reaction::EnzymeReaction) -> (children::Vector{Union{Mechanism,
  AllostericMechanism}}, parent_of::Dict, expand_failures::Vector{FitFailure})`.

- [ ] **Step 1: Write the characterization test**

Add to `test/test_identify_rate_equation.jl`. The test defines its own small
reaction — a uni-uni with a dead-end inhibitor, so expansion produces children —
so it does not depend on consts from other test files:

```julia
@testset "_expand_parents parallel equivalence" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: I
    end
    mechs = collect(EnzymeRates.init_mechanisms(rxn))
    to_expand = EnzymeRates.BatchEntry[
        EnzymeRates.BatchEntry(
            m, 2, 0.0, :Success, hash(m),
            (mechanism_type = string(typeof(EnzymeRates.compile_mechanism(m))),))
        for m in mechs]

    # Inline serial reference = the loop being replaced.
    function serial_expand_reference(to_expand, reaction)
        parent_of = Dict{Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism},
                         @NamedTuple{mechanism_type::String, n_params::Int}}()
        children = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
        fails = EnzymeRates.FitFailure[]
        for pe in to_expand
            kids, failure = EnzymeRates._expand_parent(pe.mech, reaction)
            failure === nothing || push!(fails, failure)
            for child in kids
                haskey(parent_of, child) && continue
                parent_of[child] = (mechanism_type = pe.row.mechanism_type,
                                    n_params = pe.n_params)
                push!(children, child)
            end
        end
        (children, parent_of, fails)
    end

    gc, gp, gf = EnzymeRates._expand_parents(to_expand, rxn)
    rc, rp, rf = serial_expand_reference(to_expand, rxn)

    @test gc == rc            # same children, same order, same first-parent dedup
    @test gp == rp            # same parent_of map
    @test length(gf) == length(rf)
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: the `_expand_parents parallel equivalence` testset errors —
`_expand_parents` does not exist yet (`UndefVarError: _expand_parents`). (Faster
iteration: `julia --project -e 'using EnzymeRates, Test, Random, LinearAlgebra;
include("test/mechanism_definitions_for_test_enzyme_derivation.jl");
include("test/test_identify_rate_equation.jl")'`.)

- [ ] **Step 3: Add the `_expand_parents` function**

Add near `_expand_parent` in `src/identify_rate_equation.jl` (after its
definition around line 701):

```julia
"""
Expand every selected parent into its children across the workers, then merge
serially. `pmap` preserves input order, so iterating `zip(to_expand, results)`
reproduces the serial loop's first-parent-wins dedup, child order, and failure
order exactly. Returns `(children, parent_of, expand_failures)`.
"""
function _expand_parents(to_expand::Vector{BatchEntry},
                         reaction::EnzymeReaction)
    results = pmap(m -> _expand_parent(m, reaction),
                   [pe.mech for pe in to_expand])
    parent_of = Dict{Union{Mechanism, AllostericMechanism},
                     @NamedTuple{mechanism_type::String, n_params::Int}}()
    children = Union{Mechanism, AllostericMechanism}[]
    expand_failures = FitFailure[]
    for (pe, (kids, failure)) in zip(to_expand, results)
        failure === nothing || push!(expand_failures, failure)
        for child in kids
            haskey(parent_of, child) && continue
            parent_of[child] = (mechanism_type = pe.row.mechanism_type,
                                n_params = pe.n_params)
            push!(children, child)
        end
    end
    (children, parent_of, expand_failures)
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: the same command as Step 2.
Expected: PASS.

- [ ] **Step 5: Call `_expand_parents` from `_beam_search`**

Replace the serial block at `src/identify_rate_equation.jl:828-841` (the
`parent_of` / `children` / `expand_failures` construction, from the `parent_of =
Dict...` line through the closing `end` of the `for pe in to_expand` loop) with:

```julia
            children, parent_of, expand_failures =
                _expand_parents(to_expand, prob.reaction)
```

Leave the surrounding code unchanged — the comment above the block, the
following `_process_batch(children, prob; ..., parent_of, ...)` call, and the
`append!(child_failures, expand_failures)` line all still apply.

- [ ] **Step 6: Run the full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: all green.

- [ ] **Step 7: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Parallelize per-iteration expansion via _expand_parents"
```

---

### Task 3: Confirm the speedup with workers

This is a manual verification — the empirical check the spec calls for, not an
automated test. It confirms both that the parallel paths serialize correctly
across real worker processes and that wall time drops.

**Files:**
- Create: `docs/hpc_results/parallelization_benchmark.md` (record the numbers)

- [ ] **Step 1: Benchmark seed enumeration serially vs across workers**

Write a throwaway script (scratchpad) that builds the PFKP reaction (from
`docs/hpc_results/pfkp_hpc_results/identify_pfkp.jl`) and times
`seed_mechanisms(pfk_rxn, Set([:ATP,:ADP,:Phosphate,:F26BP,:Citrate]),
Set{Symbol}())` twice: once with no workers, once after `addprocs(3);
@everywhere using EnzymeRates`. Confirm the seed vectors are equal (`==`) and
record both wall times. Use a smaller required set (e.g. 3 regulators, ~238 s
serial) if the 5-regulator case is too slow to iterate on.

- [ ] **Step 2: Record the result**

Write `docs/hpc_results/parallelization_benchmark.md` with the serial and
`addprocs` wall times for seed enumeration and a note confirming the seed
vectors matched across the two runs. State the worker count used.

- [ ] **Step 3: Commit**

```bash
git add docs/hpc_results/parallelization_benchmark.md
git commit -m "Record seed-enumeration parallelization benchmark"
```

---

## Notes for the implementer

- `pmap` inside `seed_mechanisms` and `_expand_parents` sends only mechanisms to
  workers; all bookkeeping (`visited`, `seeds`, `parent_of`) stays on the main
  node. Do not ship `BatchEntry` wrappers to workers.
- If a wave in `seed_mechanisms` is small, `pmap`'s per-task overhead may show;
  `pmap(...; batch_size=n)` amortizes it. Treat this as tuning, only if the
  benchmark shows it matters — it does not change results.
- Every equivalence test compares with `==` on `Vector`s of mechanisms and on
  the `parent_of` `Dict`; `Mechanism`/`AllostericMechanism` define structural
  `==`/`hash`, so these comparisons are order- and content-exact.
