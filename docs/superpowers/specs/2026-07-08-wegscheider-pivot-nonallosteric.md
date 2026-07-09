# The `-1` pivot-priority bug also fires on non-allosteric mechanisms

**TL;DR.** The “`-1` pivot priority never fires” defect you found in the Wegscheider
Gaussian-elimination kernel is **not allosteric-only**. It drops real Wegscheider
constraints on plain non-allosteric `Mechanism`s, over-counting `fitted_params` by 1.
It is the root cause of the LDH beam’s non-monotone `split −1` edges (the Issue 2-4
“futile enumeration” symptom): the split move re-partitions a group and its canonical
child *happens* to recover the constraint the parent’s grouping lost, so the child has
one fewer fitted param than the parent — a spurious `−1` edge.

## The bug (same one you have)

`thermodynamic_constr_for_rate_eq_derivation.jl`, pivot loop in
`_dependent_param_exprs_kernel`:

```julia
best_col, best_pri = 0, -1              # ← init
for c in 1:n_vars
    c in pivot_col_set && continue
    wA[i, c] == 0 && continue
    priority[c] > best_pri && (best_pri = priority[c]; best_col = c)
end
if best_col == 0
    wrhs[i] == 0 && continue            # ← silently drops the constraint (0 = 0)
    error(“Thermodynamically contradictory …“)
end
```

Free-enzyme RE binding K’s get `_step_priority == -1`. With `best_pri` initialized to
`-1`, the test `priority[c] > best_pri` is `-1 > -1 == false`, so a `-1` column is
**never** eligible as a pivot. When a Wegscheider row’s only remaining pivotable column
is a free-enzyme binding K, `best_col` stays `0`; because a Wegscheider constraint has
`rhs == 0`, the row is skipped as “redundant” and the parameter reduction it encodes is
lost.

## Non-allosteric example — LDH (`NADH + Pyruvate ⇌ Lactate + NAD`)

A plain `Mechanism` (6 kinetic groups, all catalytic, no allostery) reached in the LDH
enumeration:

```
grp1: E→ELactate            [Lactate, RE]
grp2: E→ENAD, ELactate→ELactateNAD          [NAD, RE]   (shared K_NAD_E)
grp3: E→ENADH, ELactate→ELactateNADH        [NADH, SS]  (shared kon/koff_NADH_E)
grp4: ENAD→ELactateNAD, ENADH→ELactateNADH  [Lactate, RE] (shared K_Lactate_ENADH)
grp5: ENAD→ENADPyruvate, ENADH→ENADHPyruvate [Pyruvate, RE] (shared K_Pyruvate_ENADH)
grp6: ENADHPyruvate→ELactateNAD             [iso, SS]
```

- **`fitted_params` = 7 on main (v0.1.6) — should be 6.** The missing constraint is
  `K_Lactate_E = K_Lactate_ENADH`.
- It is required by a **closed all-RE cycle** already present in this mechanism:

  ```
  E ──Lactate(RE, K_Lactate_E)──:arrow_forward: ELactate ──NAD(RE, K_NAD_E)──:arrow_forward: ELactateNAD
  E ──NAD(RE, K_NAD_E)────────:arrow_forward: ENAD ────Lactate(RE, K_Lactate_ENADH)──:arrow_forward: ELactateNAD
  ```

  Both paths hit the same `ELactateNAD`; all four steps are RE; NAD binding is shared
  (`K_NAD_E` on both legs). Detailed balance ⇒
  `K_Lactate_E · K_NAD_E = K_NAD_E · K_Lactate_ENADH` ⇒ `K_Lactate_E = K_Lactate_ENADH`.
  `K_Lactate_E` is free-enzyme binding (priority `-1`), so the pivot can’t fire on it and
  the constraint is dropped.

- Probe on main v0.1.6: `_thermodynamic_constraints` finds **3 cycles (1 Haldane, 2
  Wegscheider)**; the kernel imposes only the Haldane; the single-symbol rename patch
  (`_build_wegscheider_rename_map`) does not recover it for this grouping ⇒
  `fitted_params = 7`, thermodynamically inconsistent.

## Fix confirmed on this non-allosteric case

Initializing `best_pri = typemin(Int)` (so a `-1`-priority column becomes eligible as a
last-resort pivot; higher priorities are still preferred, so free-enzyme K’s are only
made dependent when nothing else can carry the constraint):

- this parent: `fitted_params` **7 → 6**, and its rate equation now shows the
  `K_Lactate_E = K_Lactate_ENADH` Wegscheider line;
- LDH `split` delta histogram (~1.4k-mechanism pool):
  - **before:** `{-1: 64, 0: 54, +1: 4680, +2: 496}`
  - **after:**  `{+1: 3184, +2: 496}` — strictly monotone; the `−1` *and* `delta-0`
    edges are gone, and the `+1` count drops because many spurious split variants now
    collapse to self-loops (dropped). The fix also shrinks the canonical pool
    (1513 → 1389), i.e. fewer thermodynamically-inconsistent duplicates.

(`typemin(Int)` is the minimal probe that confirms root cause — use whatever pivot
handling your plan prefers; the point is that a Wegscheider row whose only pivot is a
`-1`-priority free-enzyme K must still be pivoted, not dropped as redundant.)

## Minimal repro (run against your branch)

Directly detects the symptom: a `split` that *reduces* the fitted-param count. On a
correct kernel no such edge exists for these all-catalytic topologies; on the buggy
kernel there are dozens.

```julia
using EnzymeRates; const ER = EnzymeRates
rxn = @enzyme_reaction begin
    substrates:NADH[C21H29N7O14P2], Pyruvate[C3H4O3]
    products:Lactate[C3H6O3], NAD[C21H27N7O14P2]
    oligomeric_state:4
end
np(m) = length(ER.fitted_params(m))

# enumerate two generations (re_to_ss + split), structural seen-set
base = unique!(collect(ER.init_mechanisms(rxn)))
seen = Set(base); pool = collect(base); frontier = copy(base)
for _ in 1:2
    kids = ER.Mechanism[]
    for m in frontier
        append!(kids, ER._expand_re_to_ss(m)); append!(kids, ER._expand_split_kinetic_group(m))
    end
    newk = ER.Mechanism[]
    for c in kids; c in seen && continue; push!(seen,c); push!(pool,c); push!(newk,c); end
    frontier = newk
end

# find a split that reduces params (a −1 edge) — should be impossible on a correct kernel
hist = Dict{Int,Int}()
for P in pool, c in ER._expand_split_kinetic_group(P)
    d = np(c) - np(P); hist[d] = get(hist, d, 0) + 1
    if d < 0
        println(“BUG: split “, np(P), ” → “, np(c), ” params”)
        println(”  parent fitted: “, ER.fitted_params(P))
        println(”  child  fitted: “, ER.fitted_params(c))
    end
end
println(“split delta histogram: “, sort(collect(hist)))
# buggy main v0.1.6:  {-1: 64, 0: 54, +1: 4680, +2: 496}
```