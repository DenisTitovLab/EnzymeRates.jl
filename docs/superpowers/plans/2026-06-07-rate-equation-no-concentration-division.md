# Representation-Correct, Division-Free Rate-Equation Derivation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make a mechanism's derived rate equation depend only on *what the mechanism is* (not how its steps are written) and stay finite at zero metabolite concentration, so `identify_rate_equation` stops silently dropping mechanisms.

**Architecture:** Four coupled changes — (1) enumeration represents ping-pong covalent intermediates with metabolite **residuals** on conformation `:E` (not empty `:Estar` conformations), validated atom-conservingly; (2) the Cha/King-Altman derivation references each rapid-equilibrium segment to its **free enzyme** and computes the polynomials directly in **fractional/Laurent** form, then reduces to lowest terms over concentrations (concentration-GCD); (3) mechanisms become **canonical by construction** (canonicalization moves into the `Mechanism`/`AllostericMechanism` constructors, `_dedup_flat!` becomes non-mutating), with a textbook-oracle permutation bridge; (4) regression tests lock in atom-conservation, division-freeness at zero concentration, and representation-independence. The PR **net-removes** code (`is_estar` plumbing, `_canonicalize_mechanism!`, the entire common-denominator apparatus).

**Tech Stack:** Julia, `@generated` rate-equation derivation, the in-house `POLY = Dict{MONO, Rational{Int}}` Laurent polynomial type, Optimization.jl (untouched here). Tests: Test.jl + Aqua + JET.

**Settled design decision (this session):** the **atom-aware** option. Atom-conservation is validated where the reaction's atoms are available (enumerator / `Mechanism` constructor); the enumerator emits apo `:E` (empty `Residual`) when a covalent residue's atoms cancel, closing the catalytic cycle. The `Step` constructor keeps only a *cheap metabolite-level* check for the cases it can decide (same conformation **and** equal residual on both sides) — which still rejects the original `ENADH → EstarLactate` residual-free transmuting iso. The residual-*consuming* close (`E(B)·res(+A−P) → E(Q)`, atom-balanced but metabolite-imbalanced) is validated only by the atom-level check.

---

## Worked residual example (reference for Part 1)

`bi_bi_pp_rxn`: substrates `A[C,X]`, `B[N]`; products `P[C]`, `Q[N,X]` (so `A+B` atoms `= P+Q` atoms). Residual threading rule: at any form, `residual = Residual(added = consumed_subs ∖ on_enzyme_subs, subtracted = released_prods ∪ on_enzyme_prods)`, **reduced to `Residual()` whenever its net atoms cancel** (the enumerator already tracks residue atoms in `acc_atoms`). One ping-pong path:

| form | bound | residual | units = bound ∪ added − subtracted | step into it | metabolite check |
|------|-------|----------|------|------|------|
| `E` | — | ∅ | ∅ | (start, apo) | — |
| `E(A)` | A | ∅ | {A} | `E + A` binding | to−from = {A} ✓ |
| `E(P)·res(+A−P)` | P | +A −P | {A} | iso | ∅ ✓ |
| `E·res(+A−P)` ≡ F | — | +A −P | {A}−{P} (atoms {X}) | release P (binding) | to−from = {P} ✓ |
| `E(B)·res(+A−P)` | B | +A −P | {A,B}−{P} | `+ B` binding | to−from = {B} ✓ |
| `E(Q)` | Q | ∅ (atoms cancel → reduced) | {Q} | **iso, residual-consuming** | **metabolite ✗ / atom ✓** |
| `E` | — | ∅ | ∅ | release Q (binding) | to−from = {Q} ✓ |

Only the single residual-consuming iso (`E(B)·res(+A−P) → E(Q)`) fails a metabolite-only check; it is atom-balanced (`{X}+B = Q`) and validated by the atom check. There is exactly **one** residual-bearing form per covalent intermediate (`F ≡ E·res(+A−P)`), matching the Segel `F`-conformation ping-pong.

---

## File Structure

| File | Responsibility | Change |
|------|----------------|--------|
| `src/types.jl` | concrete types + `Step`/`Mechanism`/`AllostericMechanism` ctors + `name` chokepoint | add cheap metabolite check to `Step` ctor; move step/group canonicalization into the two mechanism ctors; delete nothing structural in `Residual`/`Species` (already support residuals) |
| `src/mechanism_enumeration.jl` | `init_mechanisms`/`expand_mechanisms`/`_dedup_flat!`, `_make_species`, `backtrack!`, `_release_products!` | residual threading; drop `is_estar`; add atom-conservation validator; delete `_canonicalize_mechanism!`; `_dedup_flat!` → `unique!` |
| `src/rate_eq_derivation.jl` | Cha/King-Altman derivation, kcat, allosteric assembly | free-enzyme reference seed; Laurent refactor (single `alpha`, delete common-denominator apparatus); apply concentration-GCD |
| `src/sym_poly_for_rate_eq_derivation.jl` | `POLY`/`MONO` Laurent algebra | add `_reduce_conc_lowest_terms`; possibly delete `_poly_div_mono` if no remaining caller |
| `test/test_rate_eq_derivation.jl` | `run_all_tests`, oracle bridge, snapshots, perf | add zero-metabolite finiteness; oracle permutation bridge; regenerate ~14 factored-form snapshots |
| `test/test_mechanism_enumeration.jl` | enumeration unit/integration | add division + atom-conservation + representation-independence tests; update count assertions |
| `.claude/CLAUDE.md` | project notes | update "verified topology counts"; remove C4 "empty-residual ping-pong" + `is_estar` notes |

---

## Task 0 — De-risk spike: Laurent derivation feasibility (Part 2)

**Goal:** Before refactoring the derivation, confirm on real mechanisms that computing the Cha polynomials directly in Laurent form is correct and meets the guardrails. **Scratch only — not committed.** Fall back to an explicit normalization pass only if a guardrail can't be met (then STOP and raise with Denis).

**Files:**
- Create: `spike_laurent.jl` (repo root, deleted in Task 13)

- [ ] **Step 1: Write the spike script**

The spike monkeypatches nothing; it reproduces the math by hand on two mechanisms to validate the building blocks the Part 2 refactor relies on. Targets: (a) a sequential `G>1` (steady-state) mechanism that currently triggers the `normalize`-only-for-`G==1` bloat — e.g. `Segel Ordered Bi Bi` from `MECHANISM_TEST_SPECS`; (b) the `Segel Ping Pong Bi Bi` fixture (a correct division-free `G=4` ping-pong). Construct LDH (`substrates NADH,Pyruvate; products Lactate,NAD; oligomeric_state 4`) and take the first few `init_mechanisms` for the division spot-check.

```julia
using EnzymeRates, Random
const ER = EnzymeRates

# 1. sym_det on a hand-built Laurent matrix: confirm negative exponents survive
let
    a = ER.poly_sym(:S); k = ER.poly_sym(:K)                # S, K
    alpha = ER.poly_mul(a, ER.POLY(ER._mono(:K => -1) => 1)) # S * K^-1  (Laurent)
    M = ER.POLY[ alpha ER.poly_one(); ER.poly_zero() ER.poly_one() ]
    d = ER.sym_det(M, 2)
    @assert haskey(d, ER._mono(:K => -1, :S => 1)) "sym_det dropped a Laurent term"
    println("sym_det Laurent OK: ", d)
end

# 2. kcat invariance under a common monomial factor (the GCD multiplies num&den
#    by the same monomial; _kcat_groups_from_polys must give identical ratios).
let spec = first(s for s in MECHANISM_TEST_SPECS if s.name == "Segel Ordered Bi Bi")
    m = spec.mechanism
    p = ER.random_reduced_params(m; rng = MersenneTwister(1))   # (test helper)
    k0 = ER._kcat_forward(m, p)
    println("kcat baseline = ", k0, "  (compare to refactor output)")
end

# 3. perf: today's normalized form already contains 1/K divisions; confirm a
#    Laurent num/den still compiles 0-alloc / <100ns for a small mechanism.
#    (Use test_rate_equation_performance after the refactor; here just eyeball
#    rate_equation_string for the absence of K^3-style big-K bloat in a G>1.)
println(ER.rate_equation_string(
    first(s for s in MECHANISM_TEST_SPECS if s.name == "Segel Ordered Bi Bi").mechanism))
```

- [ ] **Step 2: Run the spike**

Run: `julia --project -e 'include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); const MECHANISM_TEST_SPECS = build_mechanism_test_specs(); include("test/test_rate_eq_derivation.jl")' ` is heavy; instead load the package and the spec builder minimally. Simplest: `julia --project spike_laurent.jl` after `include`ing the two test support files the spike needs. Expected: assertion in Step 1 passes; kcat prints a finite positive value; the `G>1` string shows `K^3`-style big-K bloat **today** (the thing Part 2 removes).

- [ ] **Step 3: Record findings, do NOT commit**

Write 3-5 bullet findings into the plan's Task 6-8 notes (inline) if anything differs from the spec's assumptions. If `sym_det`/kcat/perf can't be met on Laurent polys, STOP and raise with Denis before Part 2.

---

## Task 1 — Regression tests that fail today (Part 4, written first)

These encode the executable spec. They fail today; Parts 1–3 make them pass.

**Files:**
- Modify: `test/test_rate_eq_derivation.jl` (inside `run_all_tests`)
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Zero-metabolite finiteness inside `run_all_tests`**

Add a `test_zero_metabolite_finite(spec)` called from `run_all_tests` (`test/test_rate_eq_derivation.jl:945`). It piggybacks on the existing `rate_equation` compile. For each metabolite, set just that one concentration to 0 (others random-positive, params random-positive) and assert the rate is **finite and nonzero** (`isfinite` and `!= 0.0`; nonzero guards the spurious `0.0`-from-`1/Inf`). Covers allosteric fixtures.

```julia
function test_zero_metabolite_finite(spec::MechanismTestSpec)
    m = spec.mechanism
    @testset "Zero-metabolite finiteness" begin
        rng = Random.MersenneTwister(777 + hash(spec.name) % 1000)
        mets = collect(metabolites(m))
        params = random_reduced_params(m; rng)
        for zeroed in mets
            cvals = Tuple(n == zeroed ? 0.0 : 0.5 + rand(rng) for n in mets)
            concs = NamedTuple{Tuple(mets)}(cvals)
            v = rate_equation(m, concs, params)
            @test isfinite(v)
            @test v != 0.0
        end
    end
end
```
Add `test_zero_metabolite_finite(spec)` to the `run_all_tests` body. NOTE: some current division-free fixtures may already pass; the headline failure is the enumeration test below.

- [ ] **Step 2: Run — expect the enumeration-derived case to fail later; fixtures may pass now**

Run: `julia --project -e 'using Pkg; Pkg.test()'` (or scope to the derivation testset). Expected: green for division-free fixtures; this step is the harness for the real failures in Step 3.

- [ ] **Step 3: Enumeration division-freeness on `bi_bi_pp_rxn` (FAILS today)**

In `test/test_mechanism_enumeration.jl`, add a testset that enumerates the `init` set of `bi_bi_pp_rxn`, compiles each (`compile_mechanism`), and for each metabolite zeroed asserts finite+nonzero. Today this errors/НaNs on the 14 ping-pong artifacts.

```julia
@testset "init division-freeness (bi_bi_pp)" begin
    mets = [:A, :B, :P, :Q]
    for m in EnzymeRates._dedup_flat!(collect(EnzymeRates.init_mechanisms(bi_bi_pp_rxn)))
        cm = EnzymeRates.compile_mechanism(m)
        params = random_reduced_params(cm; rng = Random.MersenneTwister(1))
        for zeroed in mets
            cvals = Tuple(n == zeroed ? 0.0 : 1.0 for n in mets)
            concs = NamedTuple{Tuple(mets)}(cvals)
            v = rate_equation(cm, concs, params)
            @test isfinite(v)
        end
    end
end
```
(`random_reduced_params` lives in `test/test_rate_eq_derivation.jl`; ensure include order makes it visible, or inline a local generator.)

- [ ] **Step 4: Atom-conservation test (FAILS today)**

Add `_assert_atom_conserving(m, rxn)` test helper (test-only) that, for every `Step` in every `init`/`expand` mechanism of `uni_uni_rxn`, `bi_bi_rxn`, `bi_bi_pp_rxn`, `pyruvate_carboxylase_rxn`, asserts the signed-metabolite-unit / atom rule (Part 1). Today fails because enumeration emits atom-non-conserving `:Estar` forms.

- [ ] **Step 5: Representation-independence + dedup-non-mutating (FAILS today)**

Build one mechanism two ways (two step orders) via `EnzymeRates.Mechanism(rxn, steps_orderA)` and `…orderB`; assert `==` and equal `hash`. Assert `_dedup_flat!` does not mutate a deep-copied input's `steps` ordering. Today fails because canonicalization lives in `_dedup_flat!`, not the ctor.

- [ ] **Step 6: Commit the failing tests**

```bash
git add test/test_rate_eq_derivation.jl test/test_mechanism_enumeration.jl
git commit -m "test: division-freeness, atom-conservation, representation-independence (failing)"
```

---

## Task 2 — Atom-conservation validation (Part 1 core)

**Files:**
- Modify: `src/types.jl` (`Step` ctor; new metabolite-unit helper)
- Modify: `src/mechanism_enumeration.jl` (atom-level validator with reaction atoms; called from the `Mechanism` ctor path or a `_assert_atom_conserving`)

- [ ] **Step 1: Failing unit test for the cheap metabolite check**

In `test/test_mechanism_enumeration.jl`, assert that constructing a residual-free same-conformation transmuting iso `Step(E(NADH), E(Lactate), nothing, true)` **errors** (it changes bound `{NADH}→{Lactate}` with no residual). Run → fails (no check yet).

- [ ] **Step 2: Add the metabolite-unit helper + cheap `Step`-ctor check**

In `src/types.jl`, add:
```julia
# Signed metabolite multiset of a species: bound (+1 each) ∪ residual.added (+1)
# − residual.subtracted (−1). Pure metabolite counting — no reaction atoms.
function _metabolite_units(s::Species)
    u = Dict{Symbol,Int}()
    for m in bound(s);                 u[name(m)] = get(u, name(m), 0) + 1; end
    for m in added(residual(s));       u[name(m)] = get(u, name(m), 0) + 1; end
    for m in subtracted(residual(s));  u[name(m)] = get(u, name(m), 0) - 1; end
    filter!(p -> p.second != 0, u); u
end
```
In the `Step` inner constructor, **after** the binding canonical-direction swap, add a cheap check **gated** on same conformation AND equal residual on both sides (the only case decidable without atoms):
```julia
if conformation(from_species) == conformation(to_species) &&
   residual(from_species) == residual(to_species)
    diff = _metabolite_units(to_species)
    for (k, v) in _metabolite_units(from_species); diff[k] = get(diff, k, 0) - v; end
    filter!(p -> p.second != 0, diff)
    expected = bound_metabolite === nothing ? Dict{Symbol,Int}() :
               Dict(name(bound_metabolite) => 1)
    diff == expected || error(
        "Step violates metabolite conservation: to−from = $diff, expected $expected")
end
```
This rejects the `ENADH→ELactate` iso, leaves residual-changing steps (validated atomically later) untouched.

- [ ] **Step 3: Run → cheap-check test passes**

Run the unit test from Step 1. Expected: PASS (errors as asserted).

- [ ] **Step 4: Add the atom-level validator (has the reaction)**

In `src/mechanism_enumeration.jl`, add `_assert_atom_conserving(m::Mechanism)` (and an `AllostericMechanism` method) that, using `reactants(reaction(m))` atom payloads, computes each species' atom multiset (`bound + residual.added − residual.subtracted`, expanded to atoms) and asserts: binding step → `atoms(to) − atoms(from) == atoms(bound_metabolite)`; iso step → `atoms(to) == atoms(from)`. Wire it as a call inside the `Mechanism`/`AllostericMechanism` constructor (cheap relative to enumeration; the reaction is in scope). Iso steps that the cheap `Step` check skipped (residual-consuming close) are caught here.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "feat: atom-conservation validation (cheap metabolite check + atom-level ctor check)"
```

---

## Task 3 — Residual-bearing enumeration; drop `is_estar` (Part 1)

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_make_species`, `backtrack!`, `_release_products!`)

- [ ] **Step 1: Rewrite `_make_species` to carry a `Residual`, always `:E`**

New signature (`src/mechanism_enumeration.jl:13`):
```julia
function _make_species(bound_subs::Vector{Symbol},
                       bound_prods::Vector{Symbol},
                       residual::Residual)
    mets = Metabolite[Substrate.(bound_subs)...; Product.(bound_prods)...]
    Species(mets, :E, residual)
end
```
Remove the `is_estar` doc/param. Update the file ABOUTME if it mentions `is_estar`.

- [ ] **Step 2: Thread the metabolite residual through `backtrack!`**

Replace the `is_estar`/`has_residual::Bool` plumbing with the actual residual. `backtrack!` already threads `consumed_subs`, `released_prods`, `on_enzyme_subs`, `on_enzyme_prods`, and `acc_atoms`. Build the form's residual at each `_make_species` call:
```julia
# residue metabolites = consumed substrates not currently discretely bound,
# minus released + currently-bound (iso-formed) products. Reduced to empty
# when acc_atoms (the running atom residue) is empty.
_residual_for(consumed, on_subs, released, on_prods, acc_atoms) =
    isempty(acc_atoms) ? Residual() :
    Residual(Substrate.(setdiff(consumed, on_subs)),
             Product.(vcat(released, on_prods)))
```
Keep `has_residual` as a branch predicate but derive it from `!isempty(acc_atoms)` (the covalent-intermediate state), not a separately-passed flag. The cycle-complete check (`conformation === :E && isempty(bound) && consumed==subs && released==prods`) still works because the closing form has empty `acc_atoms` → `Residual()` → apo `:E`.

NOTE (TDD-driven): the exact placement of `_residual_for(...)` at each of the ~5 `_make_species` call sites in `backtrack!`/`_release_products!` is validated by the Task 1 atom-conservation test (Task 2's atom check is the guardrail — a mis-threaded residual errors at construction). Implement one call site at a time, re-running the atom-conservation test.

- [ ] **Step 3: Update `_release_products!`**

Replace its `is_estar::Bool`/`has_residual_atoms::Bool` params with the residual threading; each release `_make_species(Symbol[], new_unreleased, residual)` uses `_residual_for(...)` computed from the release recursion's running atom dict (it already maintains `cur_atoms`/`_subtract_atoms`).

- [ ] **Step 4: Run atom-conservation + division tests**

Run: the Task 1 atom-conservation testset and the `bi_bi_pp` enumeration testset. Expected: atom-conservation PASSES; division-freeness still partially fails (needs Part 2 derivation) but compiles now (no atom-non-conserving forms → no derivation blow-ups).

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "feat: enumerate ping-pong with residuals on :E (drop is_estar)"
```

---

## Task 4 — Recompute & update mechanism counts (Part 1)

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (count assertions ~1095–1233, ~305, ~928)
- Modify: `.claude/CLAUDE.md` (verified topology counts; C4 note)

- [ ] **Step 1: Measure**

Run a one-off: print `length(init_mechanisms(rxn))` and `length(_dedup_flat!(collect(init_mechanisms(rxn))))` for `uni_uni_rxn`, `bi_bi_rxn`, `bi_bi_pp_rxn`, `ter_ter_rxn`, `ter_bi_rxn`, `pyruvate_carboxylase_rxn`, `pyruvate_dehydrogenase_rxn`. Record actual numbers. **Do not assume** — the bi-bi init may stay ~69 or change.

- [ ] **Step 2: Update assertions to measured values**

Update every count assertion in `test/test_mechanism_enumeration.jl` that changed, and the `_catalytic_topologies` count expectations (bi-bi=11, ter-ter=283, etc.) **only if** they changed (topology count should be unaffected by the residual relabel; verify). Update the `(bi_bi_pp_rxn, 2, 2)` tuple and the init-count testset.

- [ ] **Step 3: Update CLAUDE.md**

Replace the "Verified topology counts" line with measured values; remove/replace the C4 "Empty-residual ping-pong" constraint note and the "Bystander…" line only if affected; remove the `has_residual … 'enzyme is in Estar conformation'` wording. Keep edits minimal and factual (no temporal/"changed" language per project rules).

- [ ] **Step 4: Run full enumeration testset; commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'   # scope: enumeration testset if iterating
git add test/test_mechanism_enumeration.jl .claude/CLAUDE.md
git commit -m "test: update mechanism counts for residual ping-pong; refresh CLAUDE.md notes"
```

---

## Task 5 — Free-enzyme reference seed (Part 2a)

**Files:**
- Modify: `src/rate_eq_derivation.jl` (`_compute_alpha`, `:250`)

- [ ] **Step 1: Failing/observational test**

Add a temporary assertion (or use `rate_equation_string`) on a sequential bi-bi `init` mechanism that the denominator reads as `1 + [S]/K + …` (free-enzyme-referenced, constant term present). Today, after dedup moves the reference onto a bound complex, it doesn't.

- [ ] **Step 2: Seed the per-segment BFS at the free enzyme**

In `_compute_alpha`, replace `group[1]` seeds (`visited = Set{Int}([group[1]])`, `queue = [group[1]]`, and `sigma` `length==1` handling) with the **free enzyme of the segment**: among the segment's forms, the one with the fewest bound metabolites, tie-broken toward no residual, then a deterministic key. Reuse the existing notion via `_free_enz_set(mech)` / form names, e.g.:
```julia
# enz_species[i] → Species; pick segment root by (n_bound, has_residual, name)
function _segment_root(group, enz_species)
    argmin(i -> (length(bound(enz_species[i])),
                 has_residual(enz_species[i]) ? 1 : 0,
                 string(name(enz_species[i]))), group)
end
```
Use `_segment_root(group, enz_species)` as the BFS seed in place of `group[1]`. This makes the derivation order-independent on its own.

- [ ] **Step 3: Run derivation tests**

Run: the `MECHANISM_TEST_SPECS` derivation testset (`test_reference_qssa`, `test_haldane_equilibrium`, kcat). Expected: numerically unchanged (reference choice doesn't change the rate value), denominator now interpretable. The 55/69 single-valley LDH mechanisms become division-free. Snapshot tests will change — that's Task 11.

- [ ] **Step 4: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "feat: reference each RE segment to its free enzyme in _compute_alpha"
```

---

## Task 6 — Laurent derivation refactor (Part 2b) — delete the common-denominator apparatus

**Files:**
- Modify: `src/rate_eq_derivation.jl` (`_compute_alpha`, `_ss_contrib`, `_raw_symbolic_rate_polys`, `_compute_numerator`)

This is the riskiest task (de-risked by Task 0). Compute `α_i`, `D_g`, `num`, `den` directly as Laurent `POLY`s; delete the linearization.

- [ ] **Step 1: `_compute_alpha` returns a single Laurent `alpha::Vector{POLY}`**

Replace the `alpha_num`/`alpha_den` split with one Laurent `alpha[i]` per form (`α_i = [bound metabolites] / [path big-K's]`, e.g. `S * K^-1`). Build it by `poly_mul` with `poly_sym(K)` for iso edges and `poly_mul(poly_sym(m_l[1]), POLY(_mono(K=>-1)=>1))` for binding edges (negative exponent on K — `POLY` supports it). Delete `sigma_num`/`sigma_den` and the per-group sigma block (lines ~295–309); the per-group σ is now just `reduce(poly_add, alpha[i] for i in group)`. Return `alpha` only.

- [ ] **Step 2: Simplify `_ss_contrib`**

It no longer multiplies by other forms' `alpha_den` to clear denominators. New body:
```julia
function _ss_contrib(k_poly, mets, i_form, alpha)
    r = isempty(mets) ? k_poly : poly_mul(k_poly, reduce(poly_mul, poly_sym.(mets)))
    poly_mul(r, alpha[i_form])
end
```
(Drops the `group`/`alpha_den` args and the `_poly_div_mono` clearing.)

- [ ] **Step 3: `_raw_symbolic_rate_polys` — delete `normalize`, build den/num from Laurent α**

Remove the `normalize = G == 1 && …` flag and both its branches; the denominator is `Σ_g σ_g · D_g` with `σ_g = reduce(poly_add, alpha[i] for i in group)` (Laurent). Remove `normalize && (num = _poly_div_mono(num, sigma_den[1]))`. Update `_ss_contrib`/`_compute_numerator` call sites to the new signature (pass `alpha` not `alpha_num,alpha_den,groups[g]`).

- [ ] **Step 4: `_compute_numerator` — same Laurent α**

Update its `_ss_contrib` calls to the new signature; the metabolite-tracking logic is unchanged. The `nu_ref`/`abs_nu` stoichiometric handling stays (it's a scalar, not the bloat).

- [ ] **Step 5: Run kcat + perf + derivation tests**

Run: `test_kcat_rescaling`, `test_analytical_kcat`, `test_rate_equation_performance`, `test_reference_qssa` over `MECHANISM_TEST_SPECS`. Expected: kcat invariant (common monomial cancels — Task 0 confirmed); perf 0-alloc/<100 ns; QSSA numerically identical. The `G>1` big-K bloat is gone from `rate_equation_string`.

- [ ] **Step 6: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "refactor: compute Cha polynomials in Laurent form; delete common-denominator apparatus"
```

---

## Task 7 — Concentration-GCD: reduce to lowest terms over concentrations (Part 2c)

**Files:**
- Modify: `src/sym_poly_for_rate_eq_derivation.jl` (new `_reduce_conc_lowest_terms`)
- Modify: `src/rate_eq_derivation.jl` (apply as final step of `_raw_symbolic_rate_polys`)

- [ ] **Step 1: Add `_reduce_conc_lowest_terms`**

```julia
"""
Reduce num/den to lowest terms over concentrations only. For each
concentration symbol, shift its exponent in every monomial of num ∪ den so
the minimum exponent across all those monomials is 0. Never touches
parameter symbols (would drop a fitted parameter). Clears the G=1 ping-pong
1/conc coupling; identity for sequential mechanisms (they have a constant
term).
"""
function _reduce_conc_lowest_terms(num::POLY, den::POLY, conc_set::Set{Symbol})
    mins = Dict{Symbol,Int}()
    for p in (num, den), mono in keys(p), (s, e) in mono
        s in conc_set || continue
        mins[s] = haskey(mins, s) ? min(mins[s], e) : e
    end
    isempty(mins) && return num, den
    shift(p) = POLY(
        sort!(MONO([s => (s in conc_set ? e - get(mins, s, 0) : e) for (s, e) in mono]);
              by=first) => v
        for (mono, v) in p)
    shift(num), shift(den)
end
```
(Shift is by the same per-concentration `min` across num∪den, so no term goes negative and the ratio is unchanged.)

- [ ] **Step 2: Apply it in `_raw_symbolic_rate_polys`**

Just before `return num, den` (after the `_rename_symbols` post-pass), call:
```julia
conc_set = Set{Symbol}(name(s) for s in vcat(subs_species_metas, prods_species_metas, reg_metas))
num, den = _reduce_conc_lowest_terms(num, den, conc_set)
```
Build `conc_set` from `metabolites(mech)` (the existing concentration-symbol source used by `_poly_to_expr`). Reduce over **concentrations only**.

- [ ] **Step 3: Run the division tests (Part 4 headline)**

Run: `test_zero_metabolite_finite` over `MECHANISM_TEST_SPECS` and the `bi_bi_pp` init division-freeness testset. Expected: now PASS — the G=1 ping-pong `1/conc` coupling is cleared; sequential mechanisms unchanged.

- [ ] **Step 4: Remove `_poly_div_mono` if now unused**

Grep `_poly_div_mono`; if no remaining caller (it was only the `normalize` clearing), delete it from `src/sym_poly_for_rate_eq_derivation.jl`. Otherwise leave it.

- [ ] **Step 5: Commit**

```bash
git add src/sym_poly_for_rate_eq_derivation.jl src/rate_eq_derivation.jl
git commit -m "feat: concentration-GCD reduces rate equation to lowest terms over concentrations"
```

---

## Task 8 — Verify allosteric path stays division-free (Part 2)

**Files:**
- Inspect: `src/rate_eq_derivation.jl` (`_allosteric_num_den_exprs` `:1518`, `_build_allosteric_rate_body` `:1597`, `_raw_symbolic_rate_polys_allosteric`)

- [ ] **Step 1: Run allosteric fixtures' zero-metabolite test**

`_allosteric_num_den_exprs` assembles from the per-state catalytic polys (now clean) by multiply/add only — no new concentration denominator should appear. Run `test_zero_metabolite_finite` for the allosteric specs in `MECHANISM_TEST_SPECS` and the `uni_uni_allo*` enumerations. Expected: PASS. If a `Q_cat^CatN` factor reintroduces a concentration denominator, the GCD must run on the assembled allosteric num/den too — only add that if a test fails (don't pre-engineer).

- [ ] **Step 2: Commit (only if a change was needed)**

---

## Task 9 — Canonical by construction + non-mutating dedup (Part 3a)

**Files:**
- Modify: `src/types.jl` (`Mechanism`/`AllostericMechanism` ctors)
- Modify: `src/mechanism_enumeration.jl` (delete `_canonicalize_mechanism!`; `_dedup_flat!` → `unique!`)

- [ ] **Step 1: Move canonicalization into the `Mechanism` ctor**

In `src/types.jl`, after `_canonicalize_iso_groups`, sort steps within each group by `_step_canonical_key` and sort groups by their first step's key (the body of `_canonicalize_mechanism!`). The `_step_canonical_key` helper currently lives in `mechanism_enumeration.jl:1660` — move it to `types.jl` (or make it visible to the ctor). Do the same in the `AllostericMechanism` ctor (also sort `regulatory_sites` + permute `cat_allo_states` in lockstep, per the existing `_canonicalize_mechanism!(am)`).

- [ ] **Step 2: Delete `_canonicalize_mechanism!`; make `_dedup_flat!` non-mutating**

```julia
function _dedup_flat!(mechs::Vector)
    unique!(mechs)   # mechanisms are canonical at construction
    mechs
end
```
Delete `_canonicalize_mechanism!` (both methods) and `_regulatory_site_canonical_key` if now only used by it (check — it's used in the allosteric ctor sort, keep it).

- [ ] **Step 3: Run representation-independence + dedup tests**

Run: Task 1 Step 5 testset + the full enumeration testset. Expected: representation-independence PASSES; dedup still collapses duplicates; counts unchanged from Task 4 (canonicalization is idempotent with what dedup did before).

- [ ] **Step 4: Commit**

```bash
git add src/types.jl src/mechanism_enumeration.jl
git commit -m "refactor: canonicalize in Mechanism/AllostericMechanism ctors; dedup becomes unique!"
```

---

## Task 10 — Oracle permutation bridge (Part 3b)

**Files:**
- Modify: `test/test_rate_eq_derivation.jl` (`positional_params` `:113`)
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (capture as-written steps per oracle-bearing spec)

Now that the ctor reorders steps, the textbook oracles (which key `k1f`, `k2f`, … by **as-written** step position) must map through the as-written→canonical permutation. **No oracle formula changes.**

- [ ] **Step 1: Capture as-written step order**

The `@enzyme_mechanism` macro emits `EnzymeMechanism(Mechanism(...))`, which now canonicalizes — losing source order. Add an optional `as_written_steps::Union{Vector{Vector{Step}},Nothing} = nothing` field to `MechanismTestSpec`, populated for the oracle-bearing fixtures by constructing the same step list a second time as a raw `Vector{Vector{Step}}` (not run through the canonicalizing ctor) OR by reading them off the `@enzyme_mechanism` source via a helper that returns the pre-canonical list. Simplest: a test-only `@enzyme_mechanism_steps` companion or a literal `Vector{Vector{Step}}` mirroring the block. (Pick the lower-churn option during implementation; the oracle fixtures are a fixed set.)

- [ ] **Step 2: `positional_params` applies the permutation**

Compute, for each as-written step `i`, the canonical group `g` and within-group slot by matching `as_written_steps[i] == steps(mech)[g][j]` structurally (`Step` compares structurally). Assign positional `K{i}`/`k{i}f`/`k{i}r` to the **as-written** index `i` while reading values from the canonical rep. Where `as_written_steps === nothing` (non-oracle specs), keep today's behavior.

- [ ] **Step 3: Run oracle tests**

Run: `test_analytical_rate`, `test_analytical_kcat`, `test_ode_steadystate` across the oracle-bearing specs. Expected: PASS (formulas untouched, indices bridged).

- [ ] **Step 4: Commit**

```bash
git add test/test_rate_eq_derivation.jl test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "test: bridge textbook oracles through as-written→canonical step permutation"
```

---

## Task 11 — Regenerate factored-form snapshots (Part 4)

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (~14 `expected_factored_num`/`expected_factored_denom` pairs)

- [ ] **Step 1: Print the new factored forms**

For each spec with a non-`nothing` `expected_factored_*`, print `rate_equation_string(spec.mechanism)`'s num/denom (via `_extract_num_denom`). 

- [ ] **Step 2: Eyeball each against the spec's intent**

Confirm each new denominator is the readable `1 + [S]/K + …` form, no concentration in a denominator, no big-K bloat. **Eyeball — do not blindly paste.** The spec example to match: `den = 1 + NADH/K_NADH_E + NAD/K_NAD_E + Lactate·NADH/(K_Lactate_ENAD·K_NADH_E) + …`.

- [ ] **Step 3: Update the snapshot strings; run `test_factored_form`**

Paste the eyeballed strings. Run the derivation testset. Expected: `test_factored_form` PASS. Also re-check the Expr-shape / flat-string regression tests in `test/test_rate_eq_derivation.jl` (same file) and update if they assert specific shapes.

- [ ] **Step 4: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "test: regenerate factored-form snapshots to reduced division-free equations"
```

---

## Task 12 — Full suite + LDH confirmation + cleanup (Part 4 / TDD #5)

**Files:**
- Delete: `spike_laurent.jl`
- Inspect: `.claude/CLAUDE.md` (final pass)

- [ ] **Step 1: Run the full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`. Expected: ALL green, including Aqua/JET, `test_rate_equation_performance`, kcat/rescaling/scale-invariance.

- [ ] **Step 2: LDH end-to-end confirmation**

One-off script: build the LDH reaction (`substrates NADH,Pyruvate; products Lactate,NAD; oligomeric_state 4` — supply atom brackets so `EnzymeReaction` mass-balances), run `identify_rate_equation` on a small synthetic dataset that includes reverse initial velocity (`NADH=Pyruvate=0` rows), and confirm reverse-direction mechanisms now fit (no NaN/Inf drops). Compare candidate-survival count to the spec's 69/69-then-65/69 baseline qualitatively. (Confirmation, not a committed test, unless a compact version fits the suite.)

- [ ] **Step 3: Delete the spike; final CLAUDE.md pass**

```bash
rm -f spike_laurent.jl
```
Ensure CLAUDE.md "Known Issues"/architecture notes reflect: no `is_estar`, residual ping-pong, free-enzyme reference, Laurent derivation, canonical-by-construction. No "changed/old/new" language.

- [ ] **Step 4: Final commit**

```bash
git add -A   # after git status
git commit -m "chore: LDH confirmation, remove spike, finalize CLAUDE.md notes"
```

---

## Deletion checklist (the "net-remove" Denis values)

- `is_estar` parameter/flag plumbing in `_make_species`, `backtrack!`, `_release_products!` (Task 3).
- `_canonicalize_mechanism!` (both methods) (Task 9).
- `_compute_alpha` `alpha_num`/`alpha_den` split → single `alpha`; `sigma_num`/`sigma_den` block (Task 6).
- `normalize` flag + its `G==1` special-case in `_raw_symbolic_rate_polys` (Task 6).
- `_poly_div_mono` clearing in `_ss_contrib` and the `num`-clearing line; delete `_poly_div_mono` itself if unused (Tasks 6–7).

## Guardrails (must stay green)

- `test_rate_equation_performance` — 0 alloc, <100 ns (most important test). Verify after Tasks 6, 7, 9.
- kcat / rescaling / scale-invariance (`test_kcat_rescaling`, `test_analytical_kcat`).
- Aqua / JET.
- Full suite before every commit that touches `src/`.

## Open measurement points (do not assume)

- `init`/`expand`/topology counts after Part 1 (Task 4) — measure, update assertions + CLAUDE.md.
- Whether the allosteric assembly needs the GCD (Task 8) — measure via the zero-metabolite test, add only if it fails.
- The ~14 factored-form snapshot strings (Task 11) — regenerate + eyeball.

---

## Self-review

**Spec coverage:** Part 1 → Tasks 2–4 (residuals, atom-conservation, counts). Part 2 → Tasks 5–8 (free-enzyme reference, Laurent, GCD, allosteric). Part 3 → Tasks 9–10 (canonical ctor, dedup→`unique!`, oracle bridge). Part 4 → Task 1 (failing tests up front) + Tasks 11–12 (snapshots, full suite, LDH). TDD order matches the spec's (tests first, then Part 1→2→3, then confirm). De-risk spike is Task 0, first.

**Type/name consistency:** `_make_species(bound_subs, bound_prods, residual::Residual)` used in Tasks 2–3. `_compute_alpha → alpha::Vector{POLY}`, `_ss_contrib(k_poly, mets, i_form, alpha)` consistent across Tasks 6–7. `_reduce_conc_lowest_terms(num, den, conc_set::Set{Symbol})` consistent Tasks 7–8. `_segment_root(group, enz_species)` Task 5. `_metabolite_units(s::Species)` Task 2.

**Risk notes:** The two intricate bodies — residual threading in `backtrack!` (Task 3) and the Laurent `_compute_alpha` (Task 6) — are TDD-gated: Task 3 by the atom-conservation constructor check (a mis-thread errors immediately), Task 6 by Task 0's spike plus kcat/perf/QSSA. Neither is "blind." If a guardrail can't be met, STOP and raise with Denis (per the spec's fallback clause).
