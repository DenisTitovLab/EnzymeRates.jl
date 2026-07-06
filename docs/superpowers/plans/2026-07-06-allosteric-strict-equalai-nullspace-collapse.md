# Strict `:EqualAI` Allosteric Derivation (No Case-B, Nullspace Collapse) — Implementation Plan

> **As-built note ("Option 3").** Tasks 1–2's mechanic was refined during implementation: the collapse operates on per-step **affinities** (with steady-state **speeds** always free), keyed on the base **free/derived (`indep_A`) partition** rather than the per-`:NonequalAI`-group `δ_g` / binding-vs-catalytic scheme this plan sketches. This fixed a steady-state Wegscheider-box case (`m_ro`) the original mechanic left thermodynamically inconsistent, and dropped the D2 work. See the spec's "Design update (Option 3)" banner and §2 for the mechanic as built; `_split_resolution` / `_collapse_mirror_exprs` in `src/rate_eq_derivation.jl` are the implementation. The task **structure** below (removal set, PK retag, reg-guard removal, test/gate strategy) is as executed.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `:EqualAI` mean genuinely shared between conformations, always: remove the "Case-B" silent promotion of `:EqualAI` dependents, and instead collapse the thermodynamically-forbidden part of a `:NonequalAI` split to `K_A=K_I` (keeping the honorable, nullspace part free), producing a valid degenerate equation and never an error.

**Architecture:** The honorable split space is `nullspace(C[:, N])` where `C` is the catalytic cycle-incidence matrix (already computed by `_thermodynamic_constraints`) and `N` the `:NonequalAI` group columns. A new `_split_resolution` computes, per mechanism, which `:NonequalAI` groups keep a free `K_A/K_I` split and which are *derived* (`K_I_g = K_A_g·∏(K_I_f/K_A_f)^{a}`). The per-state derivation consumes this instead of Case-B: `:EqualAI` dependents stay one shared symbol; derived `:NonequalAI` I-symbols become dependent mirror lines dropped from `fitted_params`. The whole Case-B cluster is deleted.

**Tech Stack:** Julia; the internal derivation engine (`_thermodynamic_constraints`, `_integer_nullspace`, `build_power_expr`, the per-state `_state_*` functions, the `name(p,m)` chokepoint); `Test`, `Random`.

**Spec:** `docs/superpowers/specs/2026-07-06-allosteric-strict-equalai-nullspace-collapse-design.md` — read §2 (the mechanic) and §3 (removal) before Tasks 1–2.

## Global Constraints

- **Thermodynamic consistency (hard test):** every derived mechanism must give equilibrium flux ≈ 0 (net rate = 0 at concentrations satisfying the mass-action ratio = `Keq`), for arbitrary positive parameters. A collapse that breaks this is wrong.
- **`rate_equation` perf:** allocation-free, < 100 ns/call (`test_rate_equation_performance`). The collapse is compile-time only; mirrors are as cheap as the existing `:EqualAI` regulator mirror. No runtime work added.
- **`fitted_params`/`metabolites` stay `@generated`** — the collapse changes *which* symbols they list at derivation time, never their `@generated`-ness.
- **Naming chokepoint:** every `K…`/`k…`/`V…`/`L…` symbol via `name(p, m)` (AST-walker guard in `test/test_types.jl`).
- **Golden churn is intended** for every mechanism where Case-B fired (only PK); regenerate and review `test/reference/allosteric_golden_reference.txt`.
- **Never reject.** A collapsing config derives a valid degenerate equation; it must never throw.
- 92-char lines, 4-space indent.

**Escalation note for the implementer of Tasks 1–2:** the split-resolution linear algebra (group-column merge of `C`, nullspace pivot/free partition, signed derived exponents) is the hard core. If after honest effort the exact derived-expression signs or the group-merge cannot be gotten right against the tests, report **BLOCKED** with what you tried — do NOT ship a version that passes the equilibrium-flux test by coincidence but mis-derives a coupled case. The controller will take it over.

---

### Task 1: `_split_resolution` — the honorable-split partition

**Files:**
- Modify: `src/rate_eq_derivation.jl` (add `_split_resolution` near the current Case-B block ~1256)
- Test: `test/test_split_resolution.jl` (new)
- Modify: `test/runtests.jl` (include it)

**Interfaces:**
- Consumes: `_thermodynamic_constraints(m::Mechanism) → (C::Matrix{Int}, rhs)` where `C` rows are cycles, columns are per-flat-step (`thermodynamic_constr_for_rate_eq_derivation.jl:147`); `_integer_nullspace(A::Matrix{Int})::Matrix{Int}` (rational nullspace basis, `thermodynamic_constr…:111`); `steps(am)`, `cat_allo_state(am, g)`, `_flat_steps`, `_state_mechanism(am, :A)` (`rate_eq_derivation.jl:1145`), `_group_rep`, `_free_enz_set`.
- Produces: `_split_resolution(am::AllostericMechanism) → SplitResolution` with fields `free::Vector{Int}` (indices, into `steps(am)`, of `:NonequalAI` groups that keep a free split) and `derived::Vector{Pair{Int, Vector{Pair{Int,Int}}}}` (each entry `g => [f1=>a1, f2=>a2, …]`: the derived `:NonequalAI` group `g`'s split equals `Σ aᵢ·δ_{fᵢ}` over free groups `fᵢ`, i.e. `K_I_g = K_A_g·∏ (K_I_{fᵢ}/K_A_{fᵢ})^{aᵢ}`). Consumed by Task 2.

- [ ] **Step 1: Write the failing tests**

Create `test/test_split_resolution.jl`. It builds mechanisms with the internal API and asserts the partition. (Builders are verbose; this is the correctness spec for the whole feature.)

```julia
# ABOUTME: Unit tests for _split_resolution — the honorable-split nullspace partition.
module SplitResolutionTests
using Test, EnzymeRates
const ER = EnzymeRates
const Sub=ER.Substrate; const Prd=ER.Product; const Sp=ER.Species
const St=ER.Step; const RA=ER.ReactantAtoms; const Met=ER.Metabolite
const S=Sub(:S); const P=Prd(:P); const A=Sub(:A); const B=Sub(:B)

# uni-uni: E+S->ES (bind), ES<->EP (SS catalysis), EP->E+P (release)
function uni(states)
    E=Sp(Met[],:E); ES=Sp(Met[S],:E); EP=Sp(Met[P],:E)
    rxn=ER.EnzymeReaction(RA[RA(S,[:C=>1]),RA(P,[:C=>1])], ER.RegulatorMults[], Int[2])
    steps=Vector{St}[[St(E,ES,S,true)],[St(ES,EP,nothing,false)],[St(EP,E,P,true)]]
    ER.AllostericMechanism(rxn, steps, collect(Symbol,states), 2, ER.RegulatorySite[])
end
gi(am) = [g for g in 1:length(ER.steps(am)) if ER.cat_allo_state(am,g)===:NonequalAI]

@testset "_split_resolution" begin
    # (a) single NonequalAI binding + EqualAI catalysis -> fully forbidden: NO free split.
    am = uni([:NonequalAI, :EqualAI, :EqualAI])
    r = ER._split_resolution(am)
    @test isempty(r.free)
    @test length(r.derived) == 1
    @test first(r.derived).first in gi(am)         # the S-binding group is derived (K_I=K_A)
    @test isempty(first(r.derived).second)          # derived from NO free split => K_I_S=K_A_S

    # (b) two NonequalAI bindings + EqualAI catalysis -> 1 honorable DOF.
    am2 = uni([:NonequalAI, :EqualAI, :NonequalAI])  # S-binding + P-binding NonequalAI
    r2 = ER._split_resolution(am2)
    @test length(r2.free) == 1
    @test length(r2.derived) == 1
    # the derived group's split is +1 or -1 times the free group's split (delta_P = delta_S).
    d = first(r2.derived)
    @test length(d.second) == 1
    @test d.second[1].first == r2.free[1]
    @test abs(d.second[1].second) == 1

    # (c) catalysis NonequalAI -> its reverse differs natively; the binding split is free, no collapse.
    am3 = uni([:NonequalAI, :NonequalAI, :EqualAI])
    r3 = ER._split_resolution(am3)
    @test length(r3.free) == 2      # both S-binding and catalysis keep free splits
    @test isempty(r3.derived)

    # (d) all EqualAI -> no NonequalAI groups -> empty resolution.
    r4 = ER._split_resolution(uni([:EqualAI,:EqualAI,:EqualAI]))
    @test isempty(r4.free) && isempty(r4.derived)
end
end # module
```

Add `include("test_split_resolution.jl")` to `test/runtests.jl` near the other allosteric includes.

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_split_resolution.jl")'`
Expected: `UndefVarError: _split_resolution` (function not defined).

- [ ] **Step 3: Implement `_split_resolution`**

Add to `src/rate_eq_derivation.jl` (near the Case-B block that Task 2 deletes). Algorithm:
1. `cm = _state_mechanism(am, :A)` — the full catalytic graph (both states share topology; use A).
2. `C, _ = _thermodynamic_constraints(cm)` — cycle incidence, per **flat-step** column.
3. **Merge flat-step columns to per-group.** Build `group_of[flat_index] → group g` by matching each flat step to the kinetic group whose steps contain it (`steps(am)` are the groups; `_flat_steps(cm)` the flat list, in the same canonical order the constraints use). Sum each group's flat columns into one group column, giving `Cg` (cycles × groups).
4. `N = [g for g in 1:length(steps(am)) if cat_allo_state(am,g)===:NonequalAI]`. If empty, return the empty resolution.
5. `CN = Cg[:, N]`. Compute `NS = _integer_nullspace(CN)` (columns = nullspace basis). The **free** split variables are `_integer_nullspace`'s free columns; the **pivot** variables are derived. Recover the free/pivot split: rerun the same RREF `_integer_nullspace` uses, OR read it from `NS` — a column of `NS` has a `1` at its free variable and the pivot coefficients elsewhere. Free variables = the set of column-indices `k` where `NS[k, :]` is a unit basis row; equivalently, extend `_integer_nullspace` to also return `(free_cols, pivot_cols, RREF)` and use those directly (preferred — add a 3-value method or an internal variant; keep the existing 1-value method for its other callers).
6. Map local indices back to group indices in `N`. `free = N[free_local]`. For each pivot group `p = N[pivot_local]`, its split relation from the RREF: `δ_p = Σ_f coeff·δ_f` (the negated RREF entries, as `_integer_nullspace` already computes at line 132: `NS[pc,k] = -R[r,fc]`). Emit `p => [ N[free_local_f] => coeff … ]` with nonzero `coeff` only.
7. Return `SplitResolution(free, derived)`.

Define `struct SplitResolution; free::Vector{Int}; derived::Vector{Pair{Int,Vector{Pair{Int,Int}}}}; end`.

*Note on the group-merge (Step 3):* the constraint columns are per flat-step and merge into groups by summing (a group's symbol appears once per traversal of any member step); this is the same merge the elimination kernel does at `thermodynamic_constr…:328-355` — mirror its `step_name`/`sym_col` grouping.

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_split_resolution.jl")'`
Expected: all 4 sub-testsets pass.

- [ ] **Step 5: Commit**

```bash
git add src/rate_eq_derivation.jl test/test_split_resolution.jl test/runtests.jl
git commit -m "Add _split_resolution: honorable-split nullspace partition"
```

---

### Task 2: Collapse in the derivation; delete Case-B

**Files:**
- Modify: `src/rate_eq_derivation.jl` — replace Case-B in `_state_dependent_exprs`/`_state_rate_polys`; emit collapse mirrors in `_dependent_param_exprs` (~1423) and `_build_dep_assignments` (~1525); delete `_case_b_rename_map` (~1277), `_i_nonequalai_syms` (~1256), `_state_i_case_b_renames` (~1290), `_dep_inactive_name` (~1386)
- Modify: `src/sym_poly_for_rate_eq_derivation.jl` — no change (keep `_expr_references_any`, `build_power_expr`)
- Test: `test/test_allosteric_collapse.jl` (new — the collapse-consistency regression suite)
- Modify: `test/runtests.jl` (include it)

**Interfaces:**
- Consumes: `_split_resolution(am) → SplitResolution` (Task 1); `name(p, am)`, `_emit_cat_params_for_rep(rep, state)`, `_group_rep`, `_free_enz_set`, `build_power_expr`, the reg-mirror emission pattern in `_dependent_param_exprs` (~1467-1476) and `_build_dep_assignments` (~1559-1566).
- Produces: strict-`:EqualAI` derivation. `fitted_params` for a `:NonequalAI` group in the resolution's `derived` list drops its `K_I_`/`k_I_` symbol; `rate_equation_string` shows a mirror line for it.

- [ ] **Step 1: Write the collapse-consistency tests (failing-first)**

Create `test/test_allosteric_collapse.jl`. Reuse the `uni`/`SplitResolutionTests`-style builders (copy them in; keep this file self-contained as a module). Assert:

```julia
# ABOUTME: Strict :EqualAI collapse regression — no Case-B promotion; forbidden splits
# ABOUTME: collapse to K_A=K_I, honorable splits stay free, all equations thermo-consistent.
module AllostericCollapseTests
using Test, EnzymeRates, Random
const ER = EnzymeRates
const Sub=ER.Substrate; const Prd=ER.Product; const Sp=ER.Species
const St=ER.Step; const RA=ER.ReactantAtoms; const Met=ER.Metabolite
const S=Sub(:S); const P=Prd(:P)
function uni(states)
    E=Sp(Met[],:E); ES=Sp(Met[S],:E); EP=Sp(Met[P],:E)
    rxn=ER.EnzymeReaction(RA[RA(S,[:C=>1]),RA(P,[:C=>1])], ER.RegulatorMults[], Int[2])
    steps=Vector{St}[[St(E,ES,S,true)],[St(ES,EP,nothing,false)],[St(EP,E,P,true)]]
    ER.AllostericMechanism(rxn, steps, collect(Symbol,states), 2, ER.RegulatorySite[])
end
function evalrate(am; seed=1, split=nothing)
    cem=ER.compile_mechanism(am); fp=ER.fitted_params(am); rng=MersenneTwister(seed)
    mets=collect(ER.metabolites(cem))
    base=[(k===:L ? 0.6 : 0.4+2rand(rng)) for k in fp]
    split!==nothing && (base=[(fp[i]===split[1] ? split[2] : base[i]) for i in 1:length(fp)])
    prm=NamedTuple{(fp...,:Keq,:E_total)}((base...,3.0,1.0))
    c=NamedTuple{Tuple(mets)}(ntuple(i->0.4+2rand(rng),length(mets)))
    v=real(ER.rate_equation(cem,c,prm))
    ec=NamedTuple{Tuple(mets)}(ntuple(i->(mets[i]===:P ? 3.0 : 1.0),length(mets)))
    veq=real(ER.rate_equation(cem,ec,prm))
    (fp, v, veq)
end

@testset "strict :EqualAI collapse" begin
    fp0,_,_ = evalrate(uni([:EqualAI,:EqualAI,:EqualAI]))     # baseline

    @testset "single NonequalAI binding + EqualAI catalysis -> full collapse" begin
        fp,v,veq = evalrate(uni([:NonequalAI,:EqualAI,:EqualAI]))
        @test isfinite(v); @test abs(veq) < 1e-8
        @test !(:K_I_S_E in fp)                 # I-twin dropped (collapsed to a mirror)
        s = ER.rate_equation_string(uni([:NonequalAI,:EqualAI,:EqualAI]))
        @test occursin("K_I_S_E = K_A_S_E", replace(s," "=>""))  # explicit mirror
        @test !occursin("k_I_", s)              # NO Case-B catalytic promotion
    end

    @testset "two NonequalAI bindings + EqualAI catalysis -> 1 honorable DOF" begin
        am = uni([:NonequalAI,:EqualAI,:NonequalAI])
        fp,v,veq = evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        # exactly one free I-split survives; the other is a derived mirror.
        nI = count(p->startswith(String(p),"K_I_"), fp)
        @test nI == 1
        s = replace(ER.rate_equation_string(am), " "=>"")
        @test occursin("k_I_", s) == false                    # catalysis stays shared
        # the surviving split moves the rate (identifiable)
        freeI = fp[findfirst(p->startswith(String(p),"K_I_"), fp)]
        v1 = evalrate(am; split=(freeI,1.3))[2]; v2 = evalrate(am; split=(freeI,5.0))[2]
        @test !isapprox(v1, v2)
    end

    @testset "catalysis NonequalAI -> native, no collapse, no mirror" begin
        am = uni([:NonequalAI,:NonequalAI,:EqualAI])
        fp,v,veq = evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        @test :K_I_S_E in fp                    # binding split free
        @test any(p->startswith(String(p),"k_I_"), fp)  # catalysis split free (native)
    end

    @testset "binding-Wegscheider single inner edge -> full collapse" begin
        # random-order bi-bi (two Wegscheider boxes); tag ONE inner box-independent
        # edge (EB+A->EAB) :NonequalAI, rest :EqualAI -> its split is forbidden -> collapses.
        A2=Sub(:A); B2=Sub(:B); Q2=Prd(:Q)
        E=Sp(Met[],:E); EA=Sp(Met[A2],:E); EB=Sp(Met[B2],:E); EAB=Sp(Met[A2,B2],:E)
        EPQ=Sp(Met[P,Q2],:E); EP=Sp(Met[P],:E); EQ=Sp(Met[Q2],:E)
        rxn=ER.EnzymeReaction(RA[RA(A2,[:C=>1]),RA(B2,[:N=>1]),RA(P,[:C=>1]),RA(Q2,[:N=>1])],
                              ER.RegulatorMults[], Int[2])
        sd=[St(E,EA,A2,true),St(E,EB,B2,true),St(EB,EAB,A2,true),St(EA,EAB,B2,true),
            St(EAB,EPQ,nothing,false),St(EP,EPQ,Q2,true),St(EQ,EPQ,P,true),
            St(E,EP,P,true),St(E,EQ,Q2,true)]
        st=fill(:EqualAI,9); st[3]=:NonequalAI
        am=ER.AllostericMechanism(rxn, Vector{St}[[s] for s in sd], st, 2, ER.RegulatorySite[])
        cem=ER.compile_mechanism(am); fp=ER.fitted_params(am)
        @test !(:K_I_A_EB in fp)                            # forbidden split collapsed
        s=replace(ER.rate_equation_string(am)," "=>"")
        @test occursin("K_I_A_EB=K_A_A_EB", s)              # explicit mirror
        rng=MersenneTwister(2)
        pv=Tuple((k===:L ? 0.6 : 0.4+2rand(rng)) for k in fp)
        prm=NamedTuple{(fp...,:Keq,:E_total)}((pv...,3.0,1.0))
        mets=collect(ER.metabolites(cem))
        ec=NamedTuple{Tuple(mets)}(ntuple(i->(mets[i] in (:A,:B) ? 1.0 : sqrt(3.0)),length(mets)))
        @test abs(real(ER.rate_equation(cem,ec,prm))) < 1e-8   # thermo-consistent
    end
end
end # module
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_allosteric_collapse.jl")'`
Expected: FAILs — current Case-B emits `k_I_EP_to_ES` (the `!occursin("k_I_", s)` assertions fail) and keeps `K_I_S_E` in `fitted_params`.

- [ ] **Step 3: Replace Case-B with the collapse in `_state_dependent_exprs` / `_state_rate_polys`**

Remove the Case-B rename from both `_state_dependent_exprs(am,:I)` (~1301-1308) and `_state_rate_polys(am,:I)` (~1169-1171) so the I-run produces its native names with NO promotion (both functions become single-path: `return dep, indep` / `return num, den`).

- [ ] **Step 4: Emit the collapse in `_dependent_param_exprs` and `_build_dep_assignments`**

Using `res = _split_resolution(am)`:
- For each `derived` entry `g => combo`: compute the I-symbol `K_I_g = name(_emit_cat_params_for_rep(_group_rep(steps(am)[g],fes), :I)…, am)` and A-symbol `K_A_g`, and the derived RHS via `build_power_expr(0//1, factors)` where `factors` are `(K_I_f => aᵢ, K_A_f => -aᵢ)` for each `f=>aᵢ` in `combo`, plus `K_A_g => 1` — i.e. `K_I_g = K_A_g·∏(K_I_f/K_A_f)^{aᵢ}`. Put `dep[K_I_g] = <RHS>` (drops it from `indep`/`fitted_params`), and in `_build_dep_assignments` `push!` the same as an `Expr(:(=), K_I_g, <RHS>)` ordered before any dependent that reads it (mirror the reg-`:EqualAI` ordering).
- `:NonequalAI` groups in `res.free` keep their native distinct `K_A_g`/`K_I_g` (unchanged).
- `:EqualAI` dependents are now a single shared symbol (the merge at ~1447-1453 no longer clobbers because there is no distinct I-name; guard the merge so a shared dependent takes its A-value once).

- [ ] **Step 5: Delete the Case-B cluster**

Grep-confirm no remaining callers, then delete `_case_b_rename_map`, `_i_nonequalai_syms`, `_state_i_case_b_renames`, `_dep_inactive_name` from `src/rate_eq_derivation.jl`. Remove the two now-dead Case-B unit tests: `test/test_rate_eq_derivation.jl` `_state_i_case_b_renames` assertion (~1660-1663) and the `_dep_inactive_name` testset (~1727-1762).

Run: `cd /home/denis.linux/.julia/dev/EnzymeRates && grep -rn "_case_b_rename_map\|_i_nonequalai_syms\|_state_i_case_b_renames\|_dep_inactive_name" src/ test/`
Expected: no matches.

- [ ] **Step 6: Run the collapse tests + the derivation suite**

Run: `julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_allosteric_collapse.jl")'`
Expected: all sub-testsets pass. Then run `test/test_split_resolution.jl` again — still green.

- [ ] **Step 7: Commit**

```bash
git add src/rate_eq_derivation.jl test/test_allosteric_collapse.jl test/runtests.jl test/test_rate_eq_derivation.jl
git commit -m "Strict :EqualAI: collapse forbidden splits, delete Case-B"
```

---

### Task 3: Remove the all-`:EqualAI` regulator guard

**Files:**
- Modify: `src/types.jl:965-969` (delete the guard)
- Modify: `test/test_types.jl` (flip two `@test_throws` guard assertions), `test/test_rate_eq_derivation.jl` (flip one)

**Interfaces:** none new.

- [ ] **Step 1: Flip the three guard-assertion tests to assert construction succeeds**

In `test/test_types.jl`, the single-ligand `:EqualAI` reg assertion (`@test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(cm, (2,(:NonequalAI,:NonequalAI,:NonequalAI)), (((:I,),2,(:EqualAI,)),))`) becomes:
```julia
        # Single-ligand :EqualAI reg site is allowed (degenerate but valid)
        @test EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalAI, :NonequalAI, :NonequalAI)),
            (((:I,), 2, (:EqualAI,)),),
        ) isa EnzymeRates.AllostericEnzymeMechanism
```
and the two-ligand one (`(((:I,:J),2,(:EqualAI,:EqualAI)),)`) becomes the same `@test … isa EnzymeRates.AllostericEnzymeMechanism` form. In `test/test_rate_eq_derivation.jl`, the DSL `@test_throws Exception eval(:(@allosteric_mechanism … allosteric_regulators: I::EqualAI …))` becomes `@test eval(:(@allosteric_mechanism …)) isa EnzymeRates.AllostericEnzymeMechanism` (the DSL returns an `AllostericEnzymeMechanism`).

- [ ] **Step 2: Run to verify they fail**

Run: `julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates; include("test/test_types.jl")'`
Expected: the two flipped `test_types.jl` assertions fail — the guard still throws.

- [ ] **Step 3: Delete the guard**

In `src/types.jl`, delete lines 965-969 (the `# All-:EqualAI site cancels identically — error` comment plus the `all(st === :EqualAI …) && error(...)` block). Leave every other per-reg-site validator and the trailing `AllostericEnzymeMechanism{…}()` construction intact.

- [ ] **Step 4: Run to verify they pass**

Run: `julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates; include("test/test_types.jl")'`
Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl test/test_rate_eq_derivation.jl
git commit -m "Allow all-:EqualAI regulator: remove degeneracy guard"
```

---

### Task 4: Retag PK to PEP `:OnlyA`; rewrite analytical + golden

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (PK spec: source tag, `pk_rate_analytical` ~2222, `analytical_kcat_fn` ~2268, `expected_n_independent_params` ~2265)
- Modify: `test/reference/allosteric_golden_reference.txt` (PK block, ~lines 25-28)

**Interfaces:** none new.

- [ ] **Step 1: Retag the PEP binding step in the PK source to `:OnlyA`**

In the PK `@allosteric_mechanism_src` block, change the PEP-binding **source** step's tag from `:NonequalAI` to `:OnlyA` (the source step, NOT a stored index — canonicalization reorders). Set `expected_n_independent_params` from 9 to 8.

- [ ] **Step 2: Run PK derivation to capture the new form**

Run: `julia --project -e 'using TestEnv; TestEnv.activate(); include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); using EnzymeRates; spec=first(s for s in MECHANISM_TEST_SPECS if s.name=="PK"); am=spec.mechanism; println(EnzymeRates.fitted_params(am)); println(EnzymeRates.rate_equation_string(am))'`
Expected: 8 fitted params, **no** `K_I_PEP_E` and **no** `k_I_EATPPyruvate_to_EADPPEP`; the T-state PEP terms and T-catalytic flux are absent (`N_T=0`).

- [ ] **Step 3: Rewrite `pk_rate_analytical` for the `:OnlyA` PEP form**

Update `pk_rate_analytical` (`~2222`) so the T-state contributes no catalytic flux (`N_T=0`) and its `Q_cat_T` carries no PEP/`K1_T` terms — the T-state cannot bind PEP. Keep `analytical_kcat_fn = p -> p.k5f` (OnlyA-substrate saturates to the R-state, `kcat=k5f`, matching HK/PFK-1). Match the exact reduced form emitted in Step 2.

- [ ] **Step 4: Regenerate PK's golden block and run the golden test**

Regenerate `test/reference/allosteric_golden_reference.txt`'s PK block (`### PK` … PARAMS lines) from the Step-2 output, then run:
`julia --project -e 'using TestEnv; TestEnv.activate(); include("test/test_allosteric_golden.jl")'`
Expected: PK passes byte-exact; the other 8 blocks unchanged.

- [ ] **Step 5: Run the PK analytical + kcat oracles**

Run the analytical-rate and analytical-kcat testsets for PK (in `test/test_rate_eq_derivation.jl`); confirm `pk_rate_analytical ≈ rate_equation` and `p->p.k5f ≈ _kcat_forward` for PK.

- [ ] **Step 6: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/reference/allosteric_golden_reference.txt
git commit -m "Retag PK PEP binding :OnlyA (strict-:EqualAI K-system); regen analytical + golden"
```

---

### Task 5: Full-suite / perf / golden / chokepoint gate

**Files:** none modified — the merge gate.

- [ ] **Step 1: Run the full suite**

Run (memory-heavy, ~11-13 min; do not run two at once): `julia --project -e 'using Pkg; Pkg.test()'`
Expected: all pass, including `test_split_resolution.jl`, `test_allosteric_collapse.jl`, the golden suite, and every existing allosteric derivation test.

- [ ] **Step 2: Confirm the perf contract held**

In the output, confirm `test_rate_equation_performance` passed (`allocs == 0`, `t < 100e-9` for every `MECHANISM_TEST_SPECS` mechanism) and `rate_equation` first-call is under its 6 s budget. If it regressed, STOP and report.

- [ ] **Step 3: Confirm golden only changed for PK**

Run: `git diff --stat main..HEAD -- test/reference/allosteric_golden_reference.txt`
Confirm the only changed golden block is PK (the eight other blocks byte-unchanged from `main`). If any other block changed, STOP and investigate — no other mechanism uses Case-B, so nothing else should move.

- [ ] **Step 4: Confirm the naming chokepoint guard passed**

Confirm the `test/test_types.jl` AST-walker testset passed (the derived-mirror `K_I_…` symbols all flow through `name(p,m)`).
