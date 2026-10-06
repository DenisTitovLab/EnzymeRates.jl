# Package Simplification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Shrink EnzymeRates' source by about 4,070 lines and its suite to at most ~550 s of
test time, with no loss of coverage or documentation and no change in behaviour beyond the
changes Denis approved.

**Architecture:** An equivalence harness records enumeration output, rate strings, parameter
lists and kcat values before any change; every behaviour-preserving commit must reproduce it.
A lazy per-mechanism cache removes the free-enzyme-set rebuild that dominates enumeration and
derivation time. Then test-speed work, test-size work, pure deletions, consolidations in
dependency order, and the approved behaviour changes, each its own commit.

**Tech Stack:** Julia 1.12, EnzymeRates (Test, TestEnv, Aqua, JET, Documenter).

**Spec:** `docs/superpowers/specs/2026-10-06-package-simplification-design.md`, with the item
catalog `docs/superpowers/specs/2026-10-06-package-simplification-items.md`.

**Procedures:** `docs/superpowers/plans/2026-10-06-package-simplification-procedures.md` holds
the commands, rules and the four procedures (R refactor, B behaviour change, S test change,
T timing). Every task brief names its procedure; implementers read that file first.

## Global Constraints

- One Julia process at a time (7.7 GB, no swap).
- Never push, switch branches or rewrite history; commit on `package-simplification` only.
- `rate_equation`: 0 allocations and < 120 ns per call for every `MECHANISM_TEST_SPECS`
  mechanism (`test_rate_equation_performance`).
- Canonical Step Form, the `name(p, m)` chokepoint and Pass-2 Wegscheider absorption stay intact.
- Behaviour-preserving commits reproduce the harness byte for byte (`HARNESS OK`).
- Comments and docstrings of surviving code stay; docs pages that name a changed internal are
  updated in the same commit.
- Every file keeps its two `ABOUTME:` lines; 92-character lines; match surrounding style.
- Test deletions only as approved catalog items, or tests whose only subject is deleted code.
- Enumeration-engine tests keep CLAUDE.md's rules: inline fixtures, exact child sets and counts,
  chemistry written the enumerator's way.
- Commit trailers: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_015DG9FSn9ZSbZhzeSDkqqnS`.
- Declined items stay untouched: thermo/T10 (only its docstring is corrected, Task 38),
  identify/FP1, FR2, FR3, B7, enum/A10, types/T19, dsl/D4, time/TT07, loc/DR1, loc/DR16,
  time/TT19.

## Review Focus

1. Mechanisms built different ways (DSL, enumerator, reversed steps, a cache filled or empty)
   must stay `==`, hash-equal and compile to the same type: Task 2 adds the test.
2. A mechanism whose cache was filled in one process and is serialized to a Distributed worker
   must behave like a fresh one: Task 2's test checks that a deserialized copy names its
   parameters identically.
3. Two distinct species that render the same name (E·{NAD, H} vs E·{NADH}) must be refused
   loudly, never merged: Task 30 adds the assertion and its test.
4. A DSL block that repeats a label (`steps:`, multiplicity lines, `multiplicity =` in
   `regulatory_site`) must raise an error naming the label: Task 35's tests.
5. Rate data with NaN, ±Inf or missing `Rate`, or `Keq ≤ 0`, must be refused at problem
   construction with a message naming the column: Task 37's tests.

---

### Task 1: Equivalence harness and timed runner

**Procedure:** none of R/B/S/T — this task builds the tools. No repo commit (the files are
git-ignored).

**Files:**
- Create: `.superpowers/harness/.gitignore` (one line: `*`)
- Create: `.superpowers/harness/harness.jl`
- Create: `.superpowers/harness/timed_runtests.jl`
- Create (by running): `.superpowers/harness/snapshot/*.txt`

**Interfaces:**
- Produces: `julia --project --heap-size-hint=2500M .superpowers/harness/harness.jl check|record`
  (exit 0 and `HARNESS OK` on match) and `.superpowers/harness/timed_runtests.jl`, used by every
  later task.

- [ ] **Step 1: Write `.superpowers/harness/.gitignore`** containing the single line `*`.

- [ ] **Step 2: Write `.superpowers/harness/harness.jl`** exactly:

```julia
# ABOUTME: Equivalence harness for the package-simplification branch: records enumeration
# ABOUTME: output, rate strings, parameter lists and kcat values, or checks them.
using EnzymeRates
const ER = EnzymeRates
const SNAP = joinpath(@__DIR__, "snapshot")
const M = Union{ER.Mechanism, ER.AllostericMechanism}

const UNIUNI = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R2
    competitive_inhibitors: R1
    oligomeric_state: 2
end
const UNIBI = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    oligomeric_state: 2
end
const BIBI = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    oligomeric_state: 2
end
const PINGPONG = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    oligomeric_state: 2
end
const TERTER = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
end

canon(m) = string(typeof(ER.compile_mechanism(m)))

# Raw output of init_mechanisms and of each expansion level (each level expands the unique
# mechanisms of the previous one), keyed "<reaction>_<label>".
function enumerations()
    lists = Pair{String, Vector{M}}[]
    function grow!(name, rxn, depth)
        level = M[m for m in ER.init_mechanisms(rxn)]
        push!(lists, "$(name)_init" => level)
        for d in 1:depth
            level = ER.expand_mechanisms(unique(level), rxn)
            push!(lists, "$(name)_depth$d" => level)
        end
        level
    end
    grow!("uniuni", UNIUNI, 2)
    push!(lists, "uniuni_seeds" =>
        M[m for m in ER.seed_mechanisms(UNIUNI, Set([:R2]), Set([:R1]))])
    grow!("unibi", UNIBI, 2)
    lvl1 = grow!("bibi", BIBI, 1)
    allo = M[m for m in unique(lvl1) if m isa ER.AllostericMechanism]
    push!(lists, "bibi_allochildren" => ER.expand_mechanisms(allo, BIBI))
    grow!("pingpong", PINGPONG, 2)
    lists
end

# At most n evenly spaced elements of v, deterministic.
sample(v, n) = isempty(v) ? v : v[unique(round.(Int, range(1, length(v); length = n)))]

function equation_record(m)
    cm = ER.compile_mechanism(m)
    io = IOBuffer()
    println(io, "## ", canon(m))
    println(io, rate_equation_string(cm))
    println(io, "fitted: ", join(ER.fitted_params(cm), ", "))
    println(io, "reduced: ", join(parameters(cm, Reduced), ", "))
    cm isa EnzymeMechanism && println(io, "full: ", join(parameters(cm, Full), ", "))
    fp = Tuple(ER.fitted_params(cm))
    vals = merge(NamedTuple{fp}(Tuple(1.0 + 0.137 * i for i in eachindex(fp))), (Keq = 2.5,))
    rescaled = try
        join(round.(collect(values(rescale_parameter_values(cm, vals))); sigdigits = 10), ", ")
    catch e
        "ERROR " * string(nameof(typeof(e)))
    end
    println(io, "rescaled: ", rescaled)
    String(take!(io))
end

function snapshot()
    files = Pair{String, Vector{String}}[]
    eqs = String[]
    for (key, ms) in enumerations()
        push!(files, "$key.txt" => [canon(m) for m in ms])
        for m in sample(unique(ms), 25)
            push!(eqs, "# $key")
            push!(eqs, equation_record(m))
        end
    end
    push!(files, "equations.txt" => eqs)
    ter = ER.init_mechanisms(TERTER)
    push!(files, "terter_init_sample.txt" =>
        ["count=$(length(ter))"; [canon(m) for m in ter[1:250:end]]])
    files
end

function main(mode)
    t0 = time()
    files = snapshot()
    if mode == "record"
        rm(SNAP; force = true, recursive = true)
        mkpath(SNAP)
        for (f, lines) in files
            write(joinpath(SNAP, f), join(lines, "\n"), "\n")
        end
        println("HARNESS RECORDED ", length(files), " files in ",
                round(time() - t0; digits = 1), " s")
        return 0
    end
    bad = 0
    for (f, lines) in files
        path = joinpath(SNAP, f)
        old = isfile(path) ? split(chomp(read(path, String)), "\n") : String[]
        new = split(chomp(join(lines, "\n") * "\n"), "\n")
        old == new && continue
        bad += 1
        println("MISMATCH ", f, ": recorded ", length(old), " lines, now ", length(new))
        shown = 0
        for i in 1:max(length(old), length(new))
            o = get(old, i, "<missing>")
            n = get(new, i, "<missing>")
            o == n && continue
            println("  line ", i, "\n    was: ", first(o, 300), "\n    now: ", first(n, 300))
            (shown += 1) == 5 && break
        end
    end
    println(bad == 0 ? "HARNESS OK" : "HARNESS FAILED ($bad files differ)",
            " in ", round(time() - t0; digits = 1), " s")
    bad == 0 ? 0 : 1
end

exit(main(get(ARGS, 1, "check")))
```

- [ ] **Step 3: Write `.superpowers/harness/timed_runtests.jl`** exactly:

```julia
# ABOUTME: Runs the EnzymeRates test files in runtests.jl order, printing every nested
# ABOUTME: testset's time and each file's wall, compile and GC seconds.
using TestEnv; TestEnv.activate()
using Test, EnzymeRates, LinearAlgebra, Random
Base.cumulative_compile_timing(true)
const T = joinpath(pkgdir(EnzymeRates), "test")
const FILES = [m.captures[1] for m in eachmatch(r"include\(\"([^\"]+\.jl)\"\)",
                   read(joinpath(T, "runtests.jl"), String))
               if m.captures[1] != "mechanism_definitions_for_test_enzyme_derivation.jl"]
t0 = time()
include(joinpath(T, "mechanism_definitions_for_test_enzyme_derivation.jl"))
println("DEFS ", round(time() - t0; digits = 1))
rows = Tuple{String, Float64, Float64, Float64}[]
try
    @testset verbose = true "EnzymeRates.jl" begin
        for f in FILES
            c0 = Base.cumulative_compile_time_ns()[1]
            g0 = Base.gc_num().total_time
            w0 = time_ns()
            @testset verbose = true "$f" begin
                include(joinpath(T, f))
            end
            push!(rows, (f, (time_ns() - w0) / 1e9,
                         (Base.cumulative_compile_time_ns()[1] - c0) / 1e9,
                         (Base.gc_num().total_time - g0) / 1e9))
        end
    end
finally
    for (f, w, c, g) in rows
        println("FILE ", rpad(f, 40), " wall ", round(w; digits = 1), " compile ",
                round(c; digits = 1), " gc ", round(g; digits = 1))
    end
    println("TOTAL ", round(time() - t0; digits = 1))
end
```

- [ ] **Step 4: Record the baseline.** Run (background, then poll for `HARNESS_EXIT`; see the
procedures file) `julia --project --heap-size-hint=2500M .superpowers/harness/harness.jl record`.
Expected: `HARNESS RECORDED 15 files`. Then run `check` the same way. Expected: `HARNESS OK`.
Report both run times and `ls -la .superpowers/harness/snapshot`.

- [ ] **Step 5: Confirm it is not tracked.** `git status --short` prints nothing.

---

### Task 2: Lazy per-mechanism cache (spec §5a)

**Procedure:** R for the cache itself (outputs identical), with a TDD test for the invariance
properties; the speed gate below replaces "watch it fail".

**Files:**
- Modify: `src/types.jl` (`Species`, `name(::Species)`, `Mechanism`, `AllostericMechanism`,
  `_rep_step`)
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl` (`_free_enz_set`,
  `_assemble_constraints`' `step_name`)
- Test: `test/test_types.jl`

**Interfaces:**
- Produces: `Species` has a stored field `name::Symbol`; `name(s::Species)` returns it.
  `Mechanism` and `AllostericMechanism` each have a last field `naming::_NamingCache`.
  `_free_enz_set(m)` returns the cached `Set{Symbol}` (callers must not mutate it).
  `_group_reps(m)::Vector{Step}` returns each kinetic group's naming representative, cached.
  `_rep_step(step, m)` reads `_group_reps(m)`.

- [ ] **Step 1: Measure the baseline speed.** Write `.superpowers/harness/speed.jl`:

```julia
# ABOUTME: Times the ping-pong depth-2 enumeration, the workload the free-enzyme cache targets.
# ABOUTME: Prints seconds, GC seconds and gigabytes allocated.
using EnzymeRates
const ER = EnzymeRates
const M = Union{ER.Mechanism, ER.AllostericMechanism}
pingpong = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    oligomeric_state: 2
end
function grow(rxn, depth)
    level = M[m for m in ER.init_mechanisms(rxn)]
    for _ in 1:depth
        level = unique!(ER.expand_mechanisms(level, rxn))
    end
    level
end
grow(pingpong, 1)
t = @timed grow(pingpong, 2)
println("SPEED pingpong depth2: ", round(t.time; digits = 1), " s, gc ",
        round(t.gctime; digits = 1), " s, ", round(t.bytes / 1e9; digits = 1), " GB, n = ",
        length(t.value))
```

Run `julia --project --heap-size-hint=2500M .superpowers/harness/speed.jl`. Expected around
85 s. Record the line in the report.

- [ ] **Step 2: Write the invariance test** at the end of `test/test_types.jl`'s top-level
testsets:

```julia
@testset "the naming cache never changes a mechanism's identity" begin
    m1 = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) ⇌ E(Q) + P
            E(Q) ⇌ E + Q
        end
    end)
    m2 = EnzymeRates.Mechanism(EnzymeRates.reaction(m1), deepcopy(EnzymeRates.steps(m1)))
    names1 = [EnzymeRates.name(p, m1) for p in EnzymeRates._enumerate_parameters_full(m1)]
    @test m1 == m2
    @test hash(m1) == hash(m2)
    @test EnzymeRates.compile_mechanism(m1) === EnzymeRates.compile_mechanism(m2)
    @test names1 == [EnzymeRates.name(p, m2) for p in EnzymeRates._enumerate_parameters_full(m2)]
    io = IOBuffer(); serialize(io, m1); seekstart(io)
    m3 = deserialize(io)
    @test m3 == m1 && hash(m3) == hash(m1)
    @test names1 == [EnzymeRates.name(p, m3) for p in EnzymeRates._enumerate_parameters_full(m3)]
    sp = EnzymeRates.to_species(first(first(EnzymeRates.steps(m1))))
    sp2 = EnzymeRates.Species(reverse(EnzymeRates.bound(sp)), EnzymeRates.conformation(sp),
                              EnzymeRates.residual(sp))
    @test sp2 == sp && hash(sp2) == hash(sp) && EnzymeRates.name(sp2) === EnzymeRates.name(sp)
end
```

Add `using Serialization` at the top of `test/test_types.jl` if it is not already loaded there
(check `test/runtests.jl` and the test target in `Project.toml`; `Serialization` is a stdlib —
if it is not available in the test environment, add it to `[extras]` and the `test` target in
`Project.toml`). If `_enumerate_parameters_full` has a different name in the current code,
use the function that `parameters(m, Full)` calls for a `Mechanism`. Adjust the fixture syntax
only if the DSL rejects it (keep a 5-step ordered bi-bi). Run the focused `test_types.jl`.
Expected: PASS before the change too (it pins invariance, not speed).

- [ ] **Step 3: Store the species name.** In `src/types.jl`:
  - Add a fourth field `name::Symbol` to `Species`. In the inner constructor, sort the bound
    list as today, then compute the name from the sorted fields with the existing rendering code
    (move the body of today's `name(s::Species)` into a helper
    `_species_name(bound, conformation, residual)::Symbol`) and call
    `new(sorted, conformation, residual, _species_name(sorted, conformation, residual))`.
  - Replace `name(s::Species)`'s body with `s.name`, keeping its docstring (it describes the
    rendered format, which is unchanged) on the helper or the accessor — keep it exactly once.
  - `==` and `hash` for `Species` stay as they are (they read `bound`, `conformation`,
    `residual` only).
  - `grep -n 'Species(' src/*.jl` and confirm no code builds a `Species` through `new` outside
    the inner constructor or reads `fieldnames(Species)`; the Sig encoder (`_to_sig`) must keep
    encoding exactly the three semantic fields.

- [ ] **Step 4: Add the mechanism cache.** In `src/types.jl`, before `struct Mechanism`:

```julia
"""
Naming data derived from a mechanism's steps, filled on first use: the free-enzyme
form names and each kinetic group's naming representative. It never takes part in
`==`, `hash`, the compiled type or display.
"""
mutable struct _NamingCache
    free_enz::Union{Nothing, Set{Symbol}}
    reps::Union{Nothing, Vector{Step}}
end
_NamingCache() = _NamingCache(nothing, nothing)
```

  Add a last field `naming::_NamingCache` to `Mechanism` and to `AllostericMechanism`, and
  pass `_NamingCache()` as the last argument of each `new(...)`. Their `==` and `hash` stay
  as they are. `grep -n 'Mechanism(' src/*.jl | grep -n 'new('` must show no other `new` call.

- [ ] **Step 5: Read through the cache.** In
`src/thermodynamic_constr_for_rate_eq_derivation.jl`, rename today's computing function to
`_compute_free_enz_set(m)` (same body, same docstring) and add:

```julia
function _free_enz_set(m::Union{Mechanism, AllostericMechanism})
    c = m.naming
    c.free_enz === nothing && (c.free_enz = _compute_free_enz_set(m))
    c.free_enz
end
```

  Then `grep -n '_free_enz_set(' src/*.jl` and read every caller: none may mutate the returned
  set (`push!`, `delete!`, `union!`, `setdiff!`, `empty!` on it). If one does, make that caller
  `copy` the set first.

  In `src/types.jl`, replace `_rep_step` with:

```julia
"""Each kinetic group's naming representative, in group order; cached per mechanism."""
function _group_reps(m::Union{Mechanism, AllostericMechanism})
    c = m.naming
    if c.reps === nothing
        fes = _free_enz_set(m)
        c.reps = Step[_group_rep(group, fes) for group in steps(m)]
    end
    c.reps
end

"""Find the kinetic group containing `step`; return its naming rep."""
function _rep_step(step::Step, m::Union{Mechanism, AllostericMechanism})
    for (g, group) in enumerate(steps(m))
        step in group && return _group_reps(m)[g]
    end
    error("Step not found in mechanism: $step")
end
```

  (Keep the singleton `_rep_step` methods; Task 19 removes them.) `_group_rep` lives in the
  thermo file, which is included after `types.jl`; a call inside a function body resolves at
  run time, so the order is fine.

- [ ] **Step 6: Render each name once in `_assemble_constraints`.** Replace
`step_name(p::Parameter) = get(rename, name(p, mech), name(p, mech))` with:

```julia
    step_name(p::Parameter) = (s = name(p, mech); get(rename, s, s))
```

- [ ] **Step 7: Verify.** Run the focused `test_types.jl` and `test_rate_eq_derivation.jl`
(PASS), the harness `check` (`HARNESS OK`), `speed.jl` (gate: at most 45 s, record the line),
and the full suite (PASS, including `test_rate_equation_performance`). If a test or golden file
printed a `Species` or `Mechanism` with Julia's default `show` (which would now include the new
fields), find it with `grep -rn 'repr(\|string(m)\|sprint(show' test/` and report it rather
than re-pinning it; such a pin means the default display leaked into an expectation.

- [ ] **Step 8: Commit** `src/types.jl`, the thermo file and `test/test_types.jl` (and
`Project.toml` if Step 2 changed it) with subject
`Cache each mechanism's free-enzyme set and species names`. Body: the profile finding (≈65 % of
enumeration samples rebuilt the free-enzyme set per parameter name), the measured before/after
`speed.jl` lines, and that `==`, `hash` and compiled types are unchanged (harness OK).

---

### Task 3: Per-type JIT and profile-guided enumeration fixes (time/TT35, time/TT36), checkpoint

**Procedure:** R, plus procedure T at the end.

**Files:**
- Modify: `src/rate_eq_derivation.jl` (entry points of `rate_equation_string`, `parameters`,
  `fitted_params`, `_kcat_forward`, `_ss_rate_constant_names` and the helpers they reach), and
  whatever the profile in Step 4 names.
- Create: `.superpowers/harness/jit.jl`

**Interfaces:**
- Consumes: Task 2's cache.
- Produces: no new names. `rate_equation` stays `@generated` per type.

- [ ] **Step 1: Read catalog items time/TT35 and time/TT36.**

- [ ] **Step 2: Write `.superpowers/harness/jit.jl`** exactly, run it, record the baseline line:

```julia
# ABOUTME: Times the first calls of the derivation entry points on 50 fresh bi-bi mechanisms.
# ABOUTME: Prints total and compile seconds; per-type JIT work targets the compile part.
using EnzymeRates
const ER = EnzymeRates
bibi = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end
function touch_all(ms)
    for m in ms
        cm = ER.compile_mechanism(m)
        rate_equation_string(cm)
        parameters(cm, Reduced)
        fp = Tuple(ER.fitted_params(cm))
        vals = merge(NamedTuple{fp}(Tuple(1.0 + 0.137 * i for i in eachindex(fp))),
                     (Keq = 2.5,))
        try
            rescale_parameter_values(cm, vals)
        catch
        end
    end
end
all_ms = ER.init_mechanisms(bibi)
touch_all(all_ms[end-4:end])                     # compile the shared helpers once
Base.cumulative_compile_timing(true)
c0 = Base.cumulative_compile_time_ns()[1]
t = @elapsed touch_all(all_ms[1:50])
c = (Base.cumulative_compile_time_ns()[1] - c0) / 1e9
println("JIT total ", round(t; digits = 1), " s, compile ", round(c; digits = 1), " s")
```

- [ ] **Step 3: Remove per-type re-specialization.** Find the non-`@generated` methods that take
an `EnzymeMechanism{Sig}` / `AllostericEnzymeMechanism{…}` instance or type and only lift it to
the concrete mechanism (`Mechanism(em)`, `AllostericMechanism(aem)`) before doing runtime work.
Mark their mechanism argument `@nospecialize` (or lift once at the exported entry point and pass
the concrete mechanism down), so the helper chain compiles once rather than once per Sig.
`@generated` functions (`rate_equation`, and any other generator whose emitted body must stay
per type) keep their signatures. Run `jit.jl` after each change; keep only changes that cut
compile time. Gate: compile time at least 30 % lower than Step 2's baseline, harness OK,
`test_rate_equation_performance` green.

- [ ] **Step 4: Profile once more.** Re-run the hyperbolic sweep profile:
`julia --project --heap-size-hint=2500M` with `Profile.@profile` around the ping-pong depth-2
enumeration of `speed.jl`, and print `Profile.print(format = :flat, sortedby = :count,
C = false, mincount = 30)` filtered to `EnzymeRates/src`. If one package function still owns
more than 20 % of samples and an output-identical fix is local (caching a value already
computed per mechanism, hoisting a loop invariant, avoiding a String render in a comparison),
make it; otherwise record the top 10 frames in the report and stop. Do not restructure
algorithms here — the consolidation tasks own that.

- [ ] **Step 5: Verify and commit.** Harness OK; full suite PASS. Commit with subject
`Compile the derivation helpers once rather than per mechanism type` (and a separate commit for
any Step 4 fix). Bodies carry the before/after `jit.jl` and `speed.jl` lines.

- [ ] **Step 6: Checkpoint (procedure T).** Run the timed suite; put the `FILE` and `TOTAL`
lines and the slowest testsets in the report. The controller compares them with the spec's
§8 projection.

---

### Task 4: Compile-budget tests (time/TT03, loc/CB2, time/TT21, loc/CB1, loc/CB3)

**Procedure:** S. **Files:** `test/test_compile_budget.jl`.

- [ ] **Step 1: Read the five catalog sections.** time/TT03 and loc/CB2 describe the same
change two ways (drop the second cold ter-ter enumeration; warm the compile-reuse gate with a
cheaper reaction); implement it once, keeping every budget assertion and the trace-compile gate.
time/TT21 runs only if its condition holds (the two measurements must not interfere); measure
and report.
- [ ] **Step 2: Edit, run the focused file before/after, record times.**
- [ ] **Step 3: Full suite, commit** (`Time the compile-reuse gate without a second ter-ter
enumeration`).

---

### Task 5: Hyperbolic sweep and division-freeness (time/TT01, loc/T19, time/TT04, time/TT14, time/TT02)

**Procedure:** S. **Files:** `test/test_mechanism_enumeration.jl` (testsets
"_hyperbolic_catalysis matches the derived denominator" and "init division-freeness
(bi_bi_pp)").

- [ ] **Step 1: Read the catalog sections.** time/TT01 = loc/T19 (one change). time/TT04 and
time/TT14 memoize verdicts; keep each only if its measurement holds. time/TT02 replaces the 264
compiled rate equations with a symbolic check of the same property (no metabolite
concentration in any denominator term's divisor, exactly what the compiled check asserted).
- [ ] **Step 2: Implement TT01, measure; TT04, measure; TT14, measure; TT02, measure.** Record
each focused-file time.
- [ ] **Step 3: Full suite; two commits** (sweep; division-freeness).

---

### Task 6: Rank oracles and dead-end rank assertions (time/TT06, time/TT09, time/TT32, loc/T11, time/TT33, loc/T20, time/TT34)

**Procedure:** S. **Files:** `test/test_mechanism_enumeration.jl` (`_testhelper_identifiable_rank`
and the testsets the catalog names).

- [ ] **Step 1: Read the sections.** Pairs that are one change: TT32 = loc/T11, TT33 = loc/T20.
TT06 evaluates the derived polynomials instead of compiling `rate_equation` for plain
mechanisms; it must return the same rank on every fixture (run old and new oracle side by side
once over all call sites and report agreement before deleting the old path). TT09 seeds the RNG
from the type string.
- [ ] **Step 2: Implement in the order TT09, TT06, TT32, TT33, TT34; focused runs after each.**
- [ ] **Step 3: Full suite; commit** (`Rank oracle without compiling rate equations`) and a
second commit for the deletions.

---

### Task 7: Depth-2, integration, ter-ter and BFS tests (time/TT08, time/TT29, loc/T22, time/TT18, time/TT11, loc/T21, time/TT10, time/TT31, loc/T16, time/TT23)

**Procedure:** S. **Files:** `test/test_mechanism_enumeration.jl`.

- [ ] **Step 1: Read the sections.** Same-change pairs: TT29 = loc/T22, TT11 = loc/T21,
TT31 = loc/T16. TT10 lands only if the BFS target is still reached within 3 generations
(measure and report the generation it is reached at).
- [ ] **Step 2: Implement one item at a time with focused runs; record times.**
- [ ] **Step 3: Full suite; commit per item group** (depth-2; ter-ter; BFS; integration;
small checks).

---

### Task 8: Derivation test speed and helpers (time/TT15, TT16, TT22, TT24, TT25, TT26, TT30, loc/DR2, loc/DR15, loc/DR17)

**Procedure:** S. **Files:** `test/test_rate_eq_derivation.jl`,
`test/mechanism_definitions_for_test_enzyme_derivation.jl`.

- [ ] **Step 1: Read the sections.** TT24 restructures `MECHANISM_TEST_SPECS` construction;
the set of specs and every field value must stay identical (compare a printed list of
`(name, n_params, metabolites)` before and after). TT26 canonicalizes metabolite order in the
HK/PFK-1/PK specs so each compiles once; their assertions stay.
- [ ] **Step 2: Implement one item at a time with focused runs.**
- [ ] **Step 3: Full suite (the perf gate must stay green); commit per item group.**

---

### Task 9: Identify and fitting test speed (time/TT05, loc/FI1, time/TT12, time/TT13, loc/ID1, time/TT17, loc/ID3, loc/ID4, time/TT20, loc/ID11, time/TT27, loc/ID2, time/TT28, identify/T1, identify/T3, loc/FI3), checkpoint

**Procedure:** S, then T. **Files:** `test/test_identify_rate_equation.jl`, `test/test_fitting.jl`.

- [ ] **Step 1: Read the sections.** Same-change pairs: TT05 = FI1, TT13 = ID1,
TT17 = ID3 + ID4, TT20 = ID11 = identify/T3, TT27 = ID2. TT12 lands only if the allosteric
identify run still recovers the generating mechanism at `max_param_count = 6` (run it three
times with different seeds and report). TT05: run the merged testset ten times at maxtime 0.5 s
and report that every run passes; if any fails, use the smallest maxtime that passes ten of ten.
- [ ] **Step 2: Implement with focused runs.**
- [ ] **Step 3: Full suite; commit per item group.**
- [ ] **Step 4: Checkpoint (procedure T).**

---

### Task 10: test_types.jl and test_accessors.jl size (loc/TY1–TY16, loc/AC1)

**Procedure:** S. **Files:** `test/test_types.jl`, `test/test_accessors.jl`.

Items: loc/TY1, TY2, TY3, TY4, TY5, TY6, TY7, TY8, TY9, TY10, TY11, TY12, TY13, TY14, TY15,
TY16, AC1. TY12 also touches `test/test_rate_eq_derivation.jl`; TY13 moves an assertion into
the file the catalog names. Keep Task 2's invariance testset.

- [ ] **Step 1: Read the sections; implement one item at a time; focused runs.**
- [ ] **Step 2: Full suite; commit** (one commit, or two if the diff passes 600 lines).

---

### Task 11: test_dsl.jl size (loc/DS1, loc/DS2, dsl/D21)

**Procedure:** S. **Files:** `test/test_dsl.jl`.

- [ ] **Step 1: Read; implement; focused run.**
- [ ] **Step 2: Full suite; commit.**

---

### Task 12: Derivation and allosteric test size (loc/DR3–DR14, loc/DR18, loc/GD1, loc/AG1–AG4, loc/CO1)

**Procedure:** S. **Files:** `test/test_rate_eq_derivation.jl`,
`test/mechanism_definitions_for_test_enzyme_derivation.jl`, `test/allosteric_ground_truth.jl`,
`test/test_allosteric_collapse.jl`, `test/test_allosteric_golden.jl`, `test/reference/*`.

Items: loc/DR3, DR4, DR5, DR6, DR7, DR8, DR9, DR10, DR11, DR12, DR13, DR14, DR18, GD1, AG1,
AG2, AG3, AG4, CO1. (loc/DR1 and loc/DR16 are declined.)

- [ ] **Step 1: Read; implement in groups (derivation file; spec definitions; allosteric
files); focused runs after each.**
- [ ] **Step 2: Full suite; one commit per group.**

---

### Task 13: Enumeration test size, part 1 (loc/T1, T2, T3, T4, T5, T6, T7, T8, T15, T17)

**Procedure:** S. **Files:** `test/test_mechanism_enumeration.jl`.

Fixtures stay inline per testset; exact child sets and counts stay. loc/T15 introduces file-level
`_testhelper_` helpers only for closures the catalog lists.

- [ ] **Step 1: Read; implement one item at a time; focused runs (the file takes minutes; run
it after every two or three items).**
- [ ] **Step 2: Full suite; commit.**

---

### Task 14: Enumeration test size, part 2 (loc/T9, T10, T12, T13, T14, T18, T23, T24)

**Procedure:** S. **Files:** `test/test_mechanism_enumeration.jl`.

- [ ] **Step 1: Read; implement; focused runs.**
- [ ] **Step 2: Full suite; commit.**

---

### Task 15: Identify and fitting test size (loc/FI2, loc/ID5–ID10, identify/T2)

**Procedure:** S. **Files:** `test/test_identify_rate_equation.jl`, `test/test_fitting.jl`.

Items: loc/FI2, ID5, ID6, ID7, ID8, ID9, ID10, identify/T2.

- [ ] **Step 1: Read; implement; focused runs.**
- [ ] **Step 2: Full suite; commit.**

---

### Task 16: Pure deletions (types/T1, types/T2, enum/T1, enum/A1, deriv/V1, deriv/V2, deriv/V3, identify/M1, loc/SRC1)

**Procedure:** R. **Files:** `src/types.jl`, `src/mechanism_enumeration.jl`,
`src/rate_eq_derivation.jl`, `src/sym_poly_for_rate_eq_derivation.jl`, `src/EnzymeRates.jl`,
and the tests the catalog names. enum/T1 moves `_assert_mechanism_invariants` into
`test/test_mechanism_enumeration.jl` as `_testhelper_assert_mechanism_invariants`.

- [ ] **Step 1: Read; delete one item at a time; focused runs.**
- [ ] **Step 2: Harness OK; full suite; commit** (`Delete dead and test-only source code`).

---

### Task 17: Value-type boilerplate (types/T5, T15, T16, T18)

**Procedure:** R. **Files:** `src/types.jl`.

types/T5 generates `==`/`hash` with one `@eval` loop. It must skip `Species` (explicit pair
kept) and must not include Task 2's `naming` field or `Species.name` in any generated method;
list the generated types and their fields in the report. Every hash value stays identical
(the harness's dedup-sensitive lists prove it).

- [ ] **Step 1: Read; implement T5, then T15, T16, T18; focused `test_types.jl` after each.**
- [ ] **Step 2: Harness OK; full suite; commit.**

---

### Task 18: Canonical Step Form internals (types/T13, then types/T14)

**Procedure:** R, two commits. **Files:** `src/types.jl`.

These rewrite the canonicalizing constructors. Each lands alone: implement T13, run the focused
`test_types.jl`, harness OK, full suite, commit; then the same for T14.

- [ ] **Step 1: T13 (one canonical-groups routine with precomputed keys); verify; commit.**
- [ ] **Step 2: T14 (precomputed Tier-2 entry kinds); verify; commit.**

---

### Task 19: Sig codec, rebuild helper and the name chokepoint (types/T6, T12, T4, T3, T7, T11)

**Procedure:** R; types/T7 is an approved removal of an internal method (procedure B step 2:
assert the 3-argument method no longer exists). **Files:** `src/types.jl`, call sites in
`src/mechanism_enumeration.jl` (T12's `_with`), `src/rate_eq_derivation.jl`,
`src/thermodynamic_constr_for_rate_eq_derivation.jl`.

Order: T6 → T12 → T4 → T3 → T7 → T11 (catalog "Cross-file critic" ordering). T6's Sig must be
byte-identical (harness). T4 decides the final name of the parameter-emit helper once; later
tasks use that name. T3 makes `name(p, m)` concrete-only and keeps Task 2's cached
representatives (do not re-add a free-enzyme computation inside `_rep_step`).

- [ ] **Step 1: Implement in order; focused runs after each.**
- [ ] **Step 2: Harness OK; full suite; commits: T6; T12; T4 + T3; T7; T11.**

---

### Task 20: Display and test-oracle accessors (types/T8 + T23, types/T9, types/T10)

**Procedure:** B for T8 + T23 (approved: `show` of an `EnzymeMechanism` always prints the
step listing; the one-line chain display goes), R for T9 and T10 (approved deletions; their
test uses become one-liners over `Mechanism(m)`). **Files:** `src/types.jl`,
`test/test_types.jl`, any test that calls a deleted accessor, `docs/src` pages that show a
mechanism's display.

- [ ] **Step 1: T8 + T23 as one change**: write `show` over the steps of one lifted mechanism
without the chain walker; re-pin the display tests to the step listing; commit alone.
- [ ] **Step 2: T9, then T10; harness OK; full suite; commit.**

---

### Task 21: Thermo compactions and polynomial helpers (thermo/T1, T3, T4, T5, T6; deriv/V4, V5, V10)

**Procedure:** R. **Files:** `src/thermodynamic_constr_for_rate_eq_derivation.jl`,
`src/sym_poly_for_rate_eq_derivation.jl`, `src/rate_eq_derivation.jl`.

thermo/T5 rewrites `_compute_free_enz_set` (Task 2's name for the computing function). V4 must
emit the same `Expr` shape (the Expr-shape and flat-string tests and the perf gate guard it).

- [ ] **Step 1: Implement one item at a time; focused `test_rate_eq_derivation.jl` and
`test_types.jl`.**
- [ ] **Step 2: Harness OK; full suite; commit (thermo; polynomial helpers).**

---

### Task 22: Dependent-parameter solve as column-order RREF (thermo/T2)

**Procedure:** R (approved rewrite of a load-bearing routine). **Files:**
`src/thermodynamic_constr_for_rate_eq_derivation.jl`.

- [ ] **Step 1: Read thermo/T2 and the critic's invariant note on it.**
- [ ] **Step 2: Before deleting the old loop, run both solvers side by side** in a scratch
script over every mechanism in the harness lists and every `MECHANISM_TEST_SPECS` mechanism;
report that the dependent sets and expressions agree on all of them.
- [ ] **Step 3: Replace; harness OK; full suite; commit.**

---

### Task 23: One RE-segment routine and the seed screen (types/T17, enum/R1, enum/R2, deriv/V6, deriv/V7, deriv/V8, then enum/F)

**Procedure:** R; enum/F is perf-gated. **Files:** `src/types.jl`,
`src/mechanism_enumeration.jl`, `src/rate_eq_derivation.jl`, the tests the items name
(deriv/V7 deletes the `_bareiss_det` testset; enum/F deletes the screen-agreement testset,
loc/T26's second half).

- [ ] **Step 1: T17 → R1 → R2 → V6 → V7 → V8 (catalog conflict resolution: T17 first, then
R1/R2, then V6–V8 destructure `_re_segment_extras`'s 5-tuple); focused runs; harness OK; full
suite; commit.**
- [ ] **Step 2: enum/F.** Time `init_mechanisms` on the ter-ter reaction and the ping-pong
bi-bi reaction (oligomeric_state 2), three runs each, with and without the screen. Land the
deletion only if both medians are at most 5 % slower; otherwise keep the screen and apply the
fallback the catalog names (#84 shared weighted union-find). Report the medians. Harness OK;
full suite; commit.

---

### Task 24: Inactive-state graph, redundant copies, dead-end move (deriv/V11, enum/RC, enum/D1, then enum/D2)

**Procedure:** R for V11 + RC + D1 (one commit); B for enum/D2 (approved: the dead-end move
judges redundancy on the built child). **Files:** `src/rate_eq_derivation.jl`,
`src/mechanism_enumeration.jl`, tests the catalog names.

- [ ] **Step 1: V11 + RC + D1; harness OK; full suite; commit.**
- [ ] **Step 2: D2 by procedure B; its test is a dead-end placement whose older copy group
becomes redundant (catalog evidence) that today emits a child failing on expansion; report the
harness diff (expected: none on the recorded reactions); commit alone.**

---

### Task 25: Parameter lists (deriv/V27, then deriv/V12, deriv/V13)

**Procedure:** B for V27 (approved: `parameters(::AllostericEnzymeMechanism, Full)` is
removed; its tests and the golden `PARAMS_FULL` lines go; the `parameters` docstring and
docs pages say Full is for plain mechanisms), then R for V12 and V13. **Files:**
`src/rate_eq_derivation.jl`, `test/test_rate_eq_derivation.jl`,
`test/reference/allosteric_golden_reference.txt`, `test/test_allosteric_golden.jl`,
`test/test_accessors.jl`, `docs/src`.

- [ ] **Step 1: V27 (procedure B); commit alone.**
- [ ] **Step 2: V12 then V13 (V12's generator lifts with `_concrete` before calling
`name(p, cm)`); harness OK; full suite; commit.**

---

### Task 26: One Wegscheider rename builder and dependent assignments (thermo/T7, deriv/V14, V15, V16, thermo/T8, deriv/V9)

**Procedure:** R. **Files:** `src/thermodynamic_constr_for_rate_eq_derivation.jl`,
`src/rate_eq_derivation.jl`.

Order: thermo/T7 (it keeps Task 2's single render) → V14 → V15 → V16 (written without
`all_params` and without the kernel) → thermo/T8 → V9. Gate for V9 and V14: the exact
`rate_equation_string` of every `MECHANISM_TEST_SPECS` mechanism in the non-competitive-inhibitor
and non-essential-activator families is unchanged (the harness plus the spec string tests).

- [ ] **Step 1: Implement in order; focused runs after each.**
- [ ] **Step 2: Harness OK; full suite; commits: T7; V14–V16 + T8; V9.**

---

### Task 27: One rate-body and one equation-text builder (deriv/V17, deriv/V18, thermo/T9)

**Procedure:** R. **Files:** `src/rate_eq_derivation.jl`,
`src/thermodynamic_constr_for_rate_eq_derivation.jl`.

V18 owns the body builder; thermo/T9 shrinks to deleting `_sorted_raw_param_symbols` and the
old `_build_rate_body` methods and moving `_destructuring_expr`. The generated body keeps
`_nest_binary` two-operand trees, `build_power_expr`'s varargs `*` and every factor order; an
explicit `return` is allowed.

- [ ] **Step 1: V17, then V18 + T9; focused runs; perf gate green.**
- [ ] **Step 2: Harness OK; full suite; commit.**

---

### Task 28: MWC normalization and kcat (deriv/V19, V20, V21, V22, then V23, then V24)

**Procedure:** R for V19–V22 (one commit; the harness's kcat lines guard them); B for V23
(approved: a descriptive error instead of `MethodError max()`) and B for V24 (approved bug fix,
TDD: write the failing numerical test the catalog describes first; if it passes on the current
code, do not change the code, report that, and keep the test).

- [ ] **Step 1: V19 → V22; harness OK; full suite; commit.**
- [ ] **Step 2: V23; commit alone.**
- [ ] **Step 3: V24; commit alone (or the test alone if the bug does not reproduce).**

---

### Task 29: Catalytic topologies (enum/S0, enum/S1)

**Procedure:** R (S1 is an approved rewrite). **Files:** `src/mechanism_enumeration.jl`.

- [ ] **Step 1: S0 (`_subsets`); focused run; commit after harness OK.**
- [ ] **Step 2: S1.** Before deleting the backtracker, keep it under a temporary name and assert
in a scratch script that the walker returns the identical vector (same elements, same order) for
every reaction in the harness and in `test/test_mechanism_enumeration.jl`'s reaction constants
(including ter-ter, ter-bi and quad-quad). Report the comparison, delete the old code, harness
OK, full suite, commit.

---

### Task 30: Seed grouping, atom checks and the dead-end builder (enum/S2, S3, S4, then S5)

**Procedure:** R for S2, S3, S4 (S4 is an approved rewrite); B for S5 (approved: throw when a
reaction has no admissible catalytic cycle). **Files:** `src/mechanism_enumeration.jl`,
`test/test_mechanism_enumeration.jl` (loc/T26's dead-end opportunities part goes with S4).

- [ ] **Step 1: S2, S3; harness OK; commit.**
- [ ] **Step 2: S4.** Add the assertion the spec requires — the mechanism constructors (or the
dead-end builder, if the constructor route is not possible) raise an error naming both species
when two distinct `Species` in one mechanism render the same name — with a test that builds E
with {NAD, H} and E with {NADH} in one mechanism and expects the error. Compare
`init_mechanisms` output (vector equality, order included) old versus new on the six reactions
the catalog names before deleting the old builder. Harness OK; full suite; commit.
- [ ] **Step 3: S5 by procedure B; commit alone.**

---

### Task 31: Flux blocks and small enumeration compactions (enum/B1, C1, SP)

**Procedure:** R. **Files:** `src/mechanism_enumeration.jl`.

- [ ] **Step 1: Implement; focused runs; harness OK; full suite; commit.**

---

### Task 32: Allosteric moves (enum/A2, A3, A4, A6, A7, A8, then A5)

**Procedure:** R for A2–A4, A6–A8; B for A5 (approved). **Files:**
`src/mechanism_enumeration.jl`, tests the items name.

A6 moves the Haldane prefilter from the `rxn` argument to `reaction(m)`; on identify paths they
are equal, and the catalog's critic note says a test passing a different `rxn` would filter
differently — check the enumeration tests for such a call and report.

- [ ] **Step 1: A2 → A3 → A4 → A6 → A7 → A8; harness OK; full suite; commit.**
- [ ] **Step 2: A5 by procedure B; commit alone.**

---

### Task 33: DSL emitters and small parsers (dsl/D1, D2, D5, D6, D11, D13, D16, D17, D20)

**Procedure:** R; D13 is listed as behaviour-preserving in the catalog but adds an error for a
product added to (or a substrate subtracted from) a residual — if the harness or a test shows
it rejects an accepted input, treat it as procedure B and report. **Files:** `src/dsl.jl`,
`test/mechanism_definitions_for_test_enzyme_derivation.jl` (the `_src` test macros, D6).

- [ ] **Step 1: Implement in the order D1, D2, D5, D6, D11, D13, D16, D17, D20; focused
`test_dsl.jl` after each.**
- [ ] **Step 2: Harness OK; full suite; commit.**

---

### Task 34: DSL step parser emits Exprs directly (dsl/D7)

**Procedure:** R (approved rewrite). **Files:** `src/dsl.jl`.

- [ ] **Step 1: Before deleting `_StepSideTerm`, compare the macro expansions old versus new**
(`macroexpand` of every `@enzyme_mechanism` / `@allosteric_mechanism` block in
`test/test_dsl.jl` and the mechanism definitions file, evaluated to mechanisms and compared with
`==`) and report agreement.
- [ ] **Step 2: Delete the old path; harness OK; full suite; commit.**

---

### Task 35: DSL decisions (dsl/D3, D8, D9, D10, D12, D14, D15, D18, D19)

**Procedure:** B for each of D3, D8, D9, D10, D12, D14, D15, D18 (one commit each, test first);
R for D19 (parser merge). **Files:** `src/dsl.jl`, `test/test_dsl.jl`.

Review Focus item 4 lives here: D3, D15 and D18 tests assert that the error message names the
repeated label.

- [ ] **Step 1: One item at a time, each with its failing test first; commit each.**
- [ ] **Step 2: D19; harness OK; full suite; commit.**

---

### Task 36: Identify and fitting consolidations (identify/V1, V4, L1, FR1, C1, R1, B1, B2, B3, B4, B5, CV1, CV2, CV3)

**Procedure:** R; CV2 and CV3 delete the test halves the catalog names (approved). **Files:**
`src/identify_rate_equation.jl`, `src/fitting.jl`, `test/test_identify_rate_equation.jl`,
`test/test_fitting.jl`.

`lb`/`ub`, `abstol`/`reltol`/`callback`, the `FittingProblem(::Mechanism, …)` overload and the
`fit_inherited` column stay (declined items). V4 depends on types/T11 (Task 19).

- [ ] **Step 1: Implement in the order V1, V4, L1, FR1, C1, R1, B1, B2, B3, B4, B5, CV1, CV2,
CV3; focused runs; `loss!` must stay allocation-free (V4 adds the allosteric zero-alloc test
the catalog names).**
- [ ] **Step 2: Harness OK; full suite; commit (validation and loss; beam selection; CV).**

---

### Task 37: Identify decisions (identify/B6, V2, V3, B8, C2)

**Procedure:** B, one commit each, test first. **Files:** `src/identify_rate_equation.jl`,
`src/fitting.jl`, tests.

Review Focus item 5 lives here: V3's tests feed NaN, Inf, -Inf and missing `Rate` and
`Keq = 0` and `Keq = -1` to both problem constructors and assert an error naming the column.

- [ ] **Step 1: B6 (one batch runner; progress-log text changes as the catalog states), V2,
V3, B8, C2; one commit each.**

---

### Task 38: Remaining approved behaviour changes (types/T20, types/T21 + enum/A9, thermo/T11, deriv/V25, types/T22) and the thermo/T10 docstring

**Procedure:** B, one commit each. **Files:** `src/types.jl`, `src/mechanism_enumeration.jl`,
`src/thermodynamic_constr_for_rate_eq_derivation.jl`, `src/rate_eq_derivation.jl`, tests,
`test/reference/allosteric_golden_reference.txt`.

- types/T20: bottomless-segment side map (TDD with the masking inhibitor copy).
- types/T21 + enum/A9: `RegulatorySite` sorts its ligands; report how many harness strings and
  which golden lines move.
- thermo/T11: one name sort in `_solve_dependent_set` (replaces deriv/V26); regenerate the
  allosteric golden file once.
- deriv/V25: annotate only Pass-2-folded ties; report the moved strings.
- types/T22: three Parameter types (RE constant, forward rate, reverse rate); names unchanged
  (harness OK expected); update the chokepoint test's type regex.
- thermo/T10 declined: correct only the `_step_priority` docstring so it matches today's rule
  (the column keeps the priority of the last step of its group in canonical order); procedure R.

- [ ] **Step 1: One item at a time, each its own commit.**

---

### Task 39: Finish

**Procedure:** R + T. **Files:** `Project.toml`, `docs/src/**`, `README.md`, `.claude/CLAUDE.md`
(only where it names a changed internal).

- [ ] **Step 1: `grep -rn` every function name deleted or renamed on this branch
(`git diff e59ebd6 --stat` and the commit log list them) across `docs/`, `README.md` and
`.claude/CLAUDE.md`; update stale mentions.**
- [ ] **Step 2: Build the docs and run the doctests** (`julia --project=docs -e 'using Pkg;
Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl")'`, or the
command `docs/make.jl` documents); fix failures.
- [ ] **Step 3: Version 0.9.0** in `Project.toml` (breaking: `parameters(::AllostericEnzymeMechanism,
Full)` removed).
- [ ] **Step 4: Full suite (with Aqua and JET); timed suite (procedure T); `wc -l src/*.jl
test/*.jl` before (at `e59ebd6`) and after.** Put both tables in the report.
- [ ] **Step 5: Commit** (`Bump the version to 0.9.0 and update the docs for the simplification`).
