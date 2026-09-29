# T5: what the package derives today, and a Step with metabolites on both sides

**Setup.**
- Package: `/home/denis.linux/.julia/dev/EnzymeRates`, branch `parameter-identifiability`, HEAD `96cdf23`, not modified.
- Julia 1.12.5, one process at a time; peak RSS 1.1 GB.
- Everything is in `…/scratchpad/ident/t5_derivation/` (TD below).
- Evidence labels (PROVED / EXACT / SAMPLED / CONJECTURE) are as defined in the task.

## 0. Verdicts

1. **The Step representation.** Today's `Step` holds one metabolite (`bound_metabolite`) and no side. `_step_sides` guesses the side from the shapes of the two species, using five branches. So whether a fused step derives correctly depends on how it is written and on its RE/SS flag.
2. **Empirical results on 43 constructs:**
   - 29 give the exact law.
   - 3 silently return exactly −v.
   - 4 crash with `BoundsError`.
   - 4 fail with a thermodynamic error because the side was guessed wrong.
   - 1 fails because no flux cut is found.
   - 2 are Theorell–Chance, which cannot be written at all (the DSL drops P).
3. **Enumeration machinery.** It rejects or ignores fused steps:
   - the atom assertion rejects them;
   - the flux-carrying predicate finds no chemistry, so RE→SS gives 0 children;
   - the dead-end move loses patterns;
   - the allosteric moves find chemistry through `is_iso`;
   - deduplication fails because fused steps have no canonical orientation.
4. **Redesign.** `Step(from, to, consumed, released, is_equilibrium)` removes all side inference and fixes every derivation failure above. It touches about 25 functions (~300 lines) and should shrink `src` slightly.
5. **Rank cost.** An exact structural rank costs 3.5–6 ms per mechanism, measured with a Julia prototype. It agrees with engine A on 12,534 of 12,534 mechanisms and is about 100× cheaper than compiling a mechanism. A rank check inside enumeration is affordable.

## 1. Definitions

Write a step X→Y as consuming the multiset C and releasing R of free metabolites. Let content(X) be the bound metabolites of X plus its covalent residual.

| kind | C | R | content | example |
|---|---|---|---|---|
| pure binding | {M} | ∅ | bound(Y) = bound(X)+M | E + S → ES |
| conformational iso | ∅ | ∅ | equal | E → E* |
| chemistry iso | ∅ | ∅ | changes | EAB → EPQ |
| fused chemistry + release | ∅ | ≠∅ | changes | EAB → EQ + P |
| fused binding + chemistry | ≠∅ | ∅ | changes | EA + B → EPQ |
| Theorell–Chance (TC) | ≠∅ | ≠∅ | changes | EA + B → EQ + P |

Package terms:
- `_step_sides(s)` returns (from, to, m_lhs, m_rhs). Here m_lhs is what is consumed going from→to and m_rhs what is released.
- A "segment" is an RE-connected component (Cha).
- The numerator "cut" is the set of SS steps whose net flux is chosen as v in `_compute_numerator`.

## 2. Part A: code audit

### (a) A step carries at most one free metabolite

| where | assumption | effect on fused steps |
|---|---|---|
| types.jl:148-176 `Step` | one field, `bound_metabolite::Union{Metabolite,Nothing}` | TC cannot be represented |
| dsl.jl:759-786 `_split_side` | errors on a 2nd metabolite on one side | `E + A + B <--> E(A,B)` is rejected with "each elementary step binds at most one metabolite" |
| dsl.jl:738-757 `_build_step_expr` (744-745) | takes the LHS metabolite, else the RHS one | **TC is silently corrupted**: `E(A) + B <--> E(Q) + P` is stored as EA + B → EQ with P dropped. The error appears only at `fitted_params` ("Cycle 1 produces metabolite change not proportional to net reaction") |
| rate_eq_derivation.jl:269-315 `_compute_alpha` | uses `m_l[1]` at 295 and 305; ignores `m_r` in the reverse branch | an RE step whose metabolite is released gives `BoundsError: attempt to access 0-element Vector{Symbol} at index [1]` |
| thermo:35-46 `_step_parameters`, types.jl:1786 `_emit_cat_params_for_rep` | a step with a metabolite is a binding | fused steps get Kd, or kon/koff |
| mechanism_enumeration.jl:1205 `_apply_equivalence_grouping` | groups steps by (bm name, RE/SS) | a fused and a pure release of the same product would be forced to share constants if a merged topology went through init's grouping |
| mechanism_enumeration.jl:139-151 `_assert_step_atom_conserving` | atoms(to) − atoms(from) = atoms(bm), i.e. uptake | every fused release fails: "atom-non-conserving step EAB → EQ (bound P): Δatoms Dict(:C => -1) ≠ Dict(:C => 1)". `expand_mechanisms` therefore errors on a merged parent |
| types.jl:538 / 568 / 896 / 945 (key, RE/SS duplicate check, Sig encoding and decoding) | a single bm | fine as long as there is at most one |

### (b) Chemistry is a ligand-free isomerization

| where | assumption | effect |
|---|---|---|
| mechanism_enumeration.jl:1353-1370 `_flux_carrying_groups` (line 1361: `edge_is_chemistry = is_iso(s)`) | chemistry = isomerization | a merged mechanism has no flux-carrying group, so `_expand_re_to_ss` gives 0 children. Measured: merged ordered RE seed 0 vs unmerged seed 4; merged uni-uni 0 vs 2. Pure conformational isos are also counted as chemistry |
| rate_eq_derivation.jl:535-541 and 572-575 `_compute_numerator` | `:chem` = non-conformational iso; central-species cuts only on iso endpoints | merged random bi-bi with RE bindings has no cut (RB2) |
| mechanism_enumeration.jl:163-172 `_assert_chemistry_is_iso` | a binding step must not change the residual | a merged ping-pong EA → F + P is rejected as a parent; merged ordered passes |
| 2004-2052 `_expand_to_allosteric` (2007), 2197-2201 `_partial_onlya_catalysis`, 2237-2240 `_expand_change_allo_state`, thermo:561-565 `_onlya_haldane_violation` | chemistry groups = iso groups | on a merged mechanism: K-type variants keep catalysis live while bindings are OnlyA (partial catalysis); `_partial_onlya_catalysis` is always false; no V-type variant; fused chemistry is treated as an ε-binding |
| thermo:103-108 `_step_priority` | iso = 20, metabolite step = 10/0/−1 | only changes which parameter the Haldane pivot eliminates (names) |
| 214-800 `_catalytic_topologies` / `_release_products!`, 703 `_iso_pattern` | builds only unfused steps | the enumerator never emits fused steps |

### (c) The metabolite's side is inferred from species shape

`_step_sides` branches (rate_eq_derivation.jl:166-185):

1. **bm ∈ bound(to) → consumed.** Correct.
2. **bm ∈ bound(from) → released.** Unreachable after the constructor's swap.
3. **SS and a Product not bound on either side → released.** Wrong if the step is written with the product on the left: `E + P <--> E(S)` becomes E → ES + P (cases U4, B4, UB3).
4. **|bound(from)| > |bound(to)| → released.** The stoichiometry is right, but:
   - an RE step then crashes `_compute_alpha` (U3, U8, B5, P4);
   - a Substrate step (e.g. `E(P) <--> E + S`) is treated as a forward binding by `_compute_numerator` (543-544; line 623 negates only product steps stored as bindings). This gives −v when that step is the cut (U6c, B6c, B6d).
5. **Otherwise → consumed.** Wrong for an RE fused release between equal-size forms: `E(A) ⇌ E(Q) + P` becomes E(A) + P → E(Q) (UB2).

**No working orientation for both flags (PROVED).** U1 and U3 store the identical step ES → E + P and differ only in the flag: SS works, RE crashes. U4b and U4 store the identical step E + P → ES: RE works, SS errors. `_flip_group_to_ss` (1285-1300) keeps the orientation. So no stored orientation of a fused product release derives under both flags, and the RE→SS move turns a working RE fused step into a failing SS one.

Other readers of `_step_sides`, correct only when the inference is:
- `_re_segment_extras` (types.jl:590), `_bottomless_re_segment` (656);
- `_hyperbolic_catalysis` (1406), `_context_form` (1633);
- `_thermodynamic_constraints` stoichiometry (thermo:175);
- `reactions` / `show` (types.jl:1442, 1205).

Code that reads from/to directly:
- `_forms_with_binding_step_native` (1740-1751) returns `from` for every step carrying the metabolite. For EAB → EQ + P it returns EAB, which is ineligible, instead of EQ. **Under R6 the dead-end move gives the merged ordered mechanism 12 children versus 16 for the unmerged one**: every pattern that places an inhibitor at EQ is lost.
- `_free_enz_set` (thermo:71-95) and `_entry_kind` (types.jl:474-489) treat any step with a metabolite as a binding.

### (d) Canonical orientation

- The Step constructor (types.jl:160-170) swaps only when bm ∈ bound(from) \ bound(to), i.e. only for pure bindings.
- `_canonical_iso_direction` (495-517) handles only isomerizations.
- Fused steps are stored as written. **U5** (`E + S <--> E(P)`; `E(P) <--> E + P`) and **U6** (the same mechanism written as `E(P) <--> E + S`; `E + P <--> E(P)`):
  - `Mechanism` equality is false and the types differ;
  - fitted names are `(koff_S_E, kon_P_E, kon_S_E)` versus `(koff_S_EP, kon_P_E, kon_S_EP)`;
  - both laws are exact. EXACT.

### (e) Parameter naming

- `name(p, m)` (types.jl:1640-1664) dispatches on `is_binding(rep)`. `_render_binding` (1611-1617) builds prefix + metabolite + "_" + from; `_render_iso` (1620) builds prefix + from + "_to_" + to.
- An SS fused release EAB → EQ + P gets **kon_P_EAB** (the forward, releasing rate) and **koff_P_EAB**. This is consistent but misleading: P never binds EAB.
- `E(Q) + P ⇌ E(A,B)` gets K_P_EQ = [EQ][P]/[EAB].
- The Kd sign flip in `_assemble_constraints` (thermo:417) applies to every RE step that carries a metabolite. That is right only when the metabolite is consumed.
- **Pre-existing:** residual forms put `+`/`-` into names (e.g. `kon_B_E_res_+A_-P`). `rate_equation_string` quotes them with `var"…"` in the Haldane lines but not in the `v =` line, which is therefore not valid Julia. This affects every enumerated ping-pong mechanism.

### Other

**`_kcat_forward`** (rate_eq_derivation.jl:895-930). When no product-free numerator pattern matches a denominator pattern, it emits `max()` with no arguments, which raises `MethodError: no method matching max()` at run time.
- U4b: an SS binding inside an RE segment makes v linear in S.
- U6c, B6c, B6d: the numerator sign is flipped.
- D3 (abortive EBQ): kcat is the ratio of the A·B coefficients, while v → 0 at saturating B. This is pre-existing semantics, not specific to fused steps.

**Allosteric derivation** (`_state_allo_mechanism`, `_reachable_from_free`). These are graph-generic. The empirical check in §3 is correct.

## 3. Part B: empirical suite

**Method (SAMPLED).**
1. Each mechanism was built with `@enzyme_mechanism` from a string (`TD/suite_b.jl`, `suite_b.jl.cases`).
2. `fitted_params` and `rate_equation_string` were recorded.
3. There were 3 random laws per mechanism: free parameters and Keq log-uniform in 10^±1, E_total = 1.
4. Each law was evaluated at 3(F+2)+4 random concentrations (10^±1.5) in BigFloat with 256 bits.
5. The intended chemistry was written in the shared JSON format (`TD/json/`), and engine A's free log-parameters were fitted to the package values by mpmath LM at 50 digits (`compare_b.py`). A correct law fits below 1e-20.
6. `_kcat_forward` was compared with v at substrates = 1e30 and products = 0.
7. Six laws first stalled at 3e-3..8e-2 and were refitted with 80 starts, reaching ≤ 6.9e-37 (`refit.py`).

The "modrank" column is the Julia prototype of §6 applied to the package's derivation.

| case | steps as written | status | pkg fitted | A fitted/rank | modrank | worst fit | kcat |
|---|---|---|---|---|---|---|---|
| U1 | `E + S <--> E(S); E(S) <--> E + P` | ok | 3 | 3/3 | 3 | 2.2e-41 | ok |
| U2 | `E + S ⇌ E(S); E(S) <--> E + P` | ok | 2 | 2/2 | 2 | 2.5e-36 | ok |
| U3 | `E + S <--> E(S); E(S) ⇌ E + P` | **BoundsError** | 2 | 2/2 | – | – | – |
| U4 | `E + S <--> E(S); E + P <--> E(S)` | **thermo error** | – | 3/3 | – | – | – |
| U4b | `E + S <--> E(S); E + P ⇌ E(S)` | ok | 2 | 2/2 | 2 | 5.5e-39 | **max() error** |
| U5 | `E + S <--> E(P); E(P) <--> E + P` | ok | 3 | 3/3 | 3 | 3.4e-39 | ok |
| U6 | `E(P) <--> E + S; E + P <--> E(P)` | ok | 3 | 3/3 | 3 | 1.8e-37 | ok |
| U7 | `E + S ⇌ E(P); E(P) <--> E + P` | ok | 2 | 2/2 | 2 | 1.9e-36 | ok |
| U8 | `E(P) ⇌ E + S; E(P) <--> E + P` | **BoundsError** | 2 | 2/2 | – | – | – |
| U9 | U1 + `E + I ⇌ E(I)` | ok | 4 | 4/4 | 4 | 1.4e-37 | ok |
| U10 | U1 + `E(S) + I ⇌ E(S, I)` | ok | 4 | 4/4 | 4 | 1.4e-36 | ok |
| UB1 | `E + A <--> E(A); E(A) <--> E(Q) + P; E(Q) <--> E + Q` | ok | 5 | 5/5 | 5 | 1.6e-37 | ok |
| UB2 | same, middle step `⇌` | **thermo error** | – | 4/4 | – | – | – |
| UB3 | middle step `E(Q) + P <--> E(A)` | **thermo error** | – | 5/5 | – | – | – |
| B1 | Segel: `E+A<-->E(A); E(A)+B<-->E(A,B); E(A,B)<-->E(Q)+P; E(Q)<-->E+Q` | ok | 7 | 7/7 | 7 | 6.9e-37 | ok |
| B2 | B1 with outer steps RE (RE, SS, SS, RE) | ok | 5 | 5/5 | 5 | 7.2e-37 | ok |
| B3 | SS, RE, SS, RE | ok | 5 | 5/5 | 5 | 6.7e-38 | ok |
| B3b | RE, RE, SS, RE | ok | 4 | 4/4 | 4 | 1.1e-38 | ok |
| B4 | B1 with `E(Q) + P <--> E(A, B)` | **thermo error** | – | 7/7 | – | – | – |
| B5 | B1 with `E(A, B) ⇌ E(Q) + P` | **BoundsError** | 6 | 6/6 | – | – | – |
| B5b | B1 with `E(Q) + P ⇌ E(A, B)` | ok | 6 | 6/6 | 6 | 3.3e-37 | ok |
| B6 | `…; E(A)+B <--> E(P,Q); E(P,Q) <--> E(Q)+P; …` | ok | 7 | 7/7 | 7 | 6.4e-37 | ok |
| B6b | B6 with `E(P,Q) <--> E(A)+B` | ok | 7 | 7/7 | 7 | 6.7e-38 | ok |
| B7 | B6 with outer steps RE | ok | 5 | 5/5 | 5 | 4.2e-38 | ok |
| P1 | textbook ping-pong: `E+A<-->E(A); E(A)<-->F+P; F+B<-->F(B); F(B)<-->E+Q`, F = `E(; residual = A - P)` | ok | 7 | 7/**6** | 6 | 1.6e-38 | ok |
| P2 | P1 with RE, SS, RE, SS | ok | 5 | 5/5 | 5 | 4.3e-37 | ok |
| P3 | P1 with SS, SS, RE, SS | ok | 6 | 6/6 | 6 | 9.5e-37 | ok |
| P4 | P1 with `E(A) ⇌ F + P` | **BoundsError** | 6 | 6/6 | – | – | – |
| P5 | P1 with `F + B <--> E(Q)` (fused binding + chemistry) | ok | 7 | 7/**6** | 6 | 9.8e-38 | ok |
| RB1 | random, X = E(A,B): 4 bindings, `X<-->E(Q)+P`, `X<-->E(P)+Q`, 2 releases, all SS | ok | 13 | 13/13 | 13 | 2.4e-37 | ok |
| RB2 | RB1 with the 4 bindings and 2 releases RE | **no cut** | 7 | 7/**6** | – | – | – |
| RB3 | RB1 with the 4 bindings RE | ok | 9 | 9/9 | 9 | 6.0e-40 | ok |
| RB4 | random RE binding, `X <--> E(Q)+P`, `E(Q) ⇌ E+Q` | ok | 5 | 5/5 | 5 | 1.9e-38 | ok |
| D1 | B1 + `E + I ⇌ E(I)` | ok | 8 | 8/8 | 8 | 5.2e-40 | ok |
| D2 | B1 + `E(Q) + I ⇌ E(I, Q)` | ok | 8 | 8/8 | 8 | 2.6e-37 | ok |
| D3 | B1 + abortive `E(Q) + B ⇌ E(B, Q)` | ok | 8 | 8/8 | 8 | 2.2e-37 | differs (see §2 Other) |
| D4 | B1 + abortive `E(A) + P ⇌ E(A, P)` | ok | 8 | 8/8 | 8 | 9.1e-37 | ok |
| D5 | B3b + `E + I ⇌ E(I)` + `E(Q) + I ⇌ E(I, Q)` | ok | 6 | 6/6 | 6 | 1.2e-37 | ok |
| TC1 | `E+A<-->E(A); E(A)+B<-->E(Q)+P; E(Q)<-->E+Q` | **P dropped, thermo error** | – | 5/5 | – | – | – |
| TC2 | TC1 with outer steps RE | **P dropped, thermo error** | – | 3/3 | – | – | – |
| U6c | `E(P) <--> E + S; E(P) ⇌ E + P` | **returns −v** | 2 | 2/2 | 2 | fit(−v) ≤ 1.7e-36 | max() error |
| B6c | `E+A⇌E(A); E(P,Q)<-->E(A)+B; E(P,Q)<-->E(Q)+P; E(Q)⇌E+Q` | **returns −v** | 5 | 5/5 | 5 | fit(−v) ≤ 4.9e-40 | max() error |
| B6d | B6c with `E(P,Q) ⇌ E(Q)+P` | **returns −v** | 4 | 4/4 | 4 | fit(−v) ≤ 2.6e-39 | max() error |

**Error texts.**
- Thermo error: "Cycle 1 produces metabolite change not proportional to net reaction".
- No cut: "rate_equation: no rapid-equilibrium-consistent reaction cut — a complete all-RE catalytic cycle exists, so the mechanism has no finite rate."

**U6c law (as printed by the package).** `v = E_total*(kon_S_EP*P/K_P_E − koff_S_EP*S)/(1 + P/K_P_E)`. This is the true law with the sign reversed. In B6c a correct SS P-release cut also existed; the tie-break on step index chose the reversed B step instead.

**What derives correctly today:**
- SS fused release written substrate side → product side;
- fused binding + chemistry written with the substrate consumed;
- RE fused steps written with the metabolite consumed;
- fused ping-pong half-reactions, SS or with mixed RE/SS;
- merged random bi-bi, except with RE bindings and RE releases together;
- dead-end and abortive complexes on merged mechanisms.

Where a mechanism builds, the package's fitted count equals engine A's in every case, including the ones that later crash.

**Allosteric check (PROVED).** Mechanism: `E + S ⇌ E(S) :: NonequalAI`, `E(S) <--> E + P :: NonequalAI`, n = 2. The package gives
`v = ((kon_A_P_ES·S/K_A_S_E − koff_A_P_ES·P)(1+S/K_A_S_E) + L(kon_I_P_ES·S/K_I_S_E − koff_I_P_ES·P)(1+S/K_I_S_E)) / ((1+S/K_A_S_E)^2 + L(1+S/K_I_S_E)^2)`,
with `koff_X = kon_X/(K_X·Keq)`.

Proof: in state X the forms E and ES form one RE segment with [ES] = [E]·S/K_X. So the per-subunit flux is [E](kf S/K_X − kr P), normalised by 1 + S/K_X, and equilibrium requires kr = kf/(K_X·Keq). The MWC combination for n subunits is (N_A r^{n−1} + L N_I t^{n−1})/(r^n + L t^n), with r = 1 + S/K_A and t = 1 + S/K_I. With n = 2 this is exactly the package's expression.

**Enumeration probes** (`TD/suite_b2.jl`, `suite_b2.log`; real-atom reactions `S[C]→P[C]` and `A[CX]+B[N]→P[C]+Q[NX]`):

| mechanism | `_assert_atom_conserving` | `_flux_carrying_groups` | RE→SS children | `expand_mechanisms` |
|---|---|---|---|---|
| merged uni-uni U2 | errors | all false | 0 | 0 children (unmerged: 5) |
| merged ordered B3b | errors | all false | 0 | **errors** (unmerged: 19 children) |

The split move gives 0 children for both mechanisms, the same as their unmerged counterparts.

## 4. Part C: Theorell–Chance

- **DSL.** `E(A) + B <--> E(Q) + P` expands without error but is stored as `EA→EQ [bound B, SS]`: P is gone (dsl.jl:744-745). `fitted_params` then errors "Cycle 1 produces metabolite change not proportional to net reaction" (TC1, TC2).
- **Two metabolites on one side.** `E + A + B <--> E(A, B)` fails at macro expansion: "step side has more than one metabolite term (A, B); each elementary step binds at most one metabolite."
- **Direct construction.** `fieldnames(Step) = (:from_species, :to_species, :bound_metabolite, :is_equilibrium)`.
  - `Step(E(A), E(Q), Substrate(:B), false)` constructs, with `_step_sides = (:EA, :EQ, [:B], [])`, so P is not in it.
  - `Step(E(A), E(Q), [B, P], false)` raises a `MethodError`.
- **Correct laws (engine A).** TC all SS: 5 fitted, rank 5. TC with RE outer steps: 3 fitted, rank 3.

## 5. Part D: Step redesign (proposal; CONJECTURE until built)

```julia
struct Step
    from::Species
    to::Species
    consumed::Vector{Metabolite}   # free metabolites taken up going from → to (sorted)
    released::Vector{Metabolite}   # free metabolites given off going from → to (sorted)
    is_equilibrium::Bool
end
```

**Derived predicates (one place).**
- `ligand(s)`: the metabolite M if the step is a pure binding (C = {M}, R = ∅, bound(to) = bound(from) + M, same residual), else `nothing`.
- `is_chemistry(s)`: bound(from) ⊎ C ≠ bound(to) ⊎ R as metabolite multisets, or the residual changes.
- A pure conformational step has C = R = ∅ and equal content.

**Canonical orientation.**
1. **Pure binding or release.** The ligand is consumed (it sits on `to`), exactly as today. Every mechanism the enumerator emits now therefore keeps its stored steps, names and Haldane choices.
2. **Chemistry-bearing steps** (chemistry iso, fused, TC). The substrate side is `from`, chosen by Tier 1 extended to free metabolites: score(side) = (#substrates, −#products) over bound ∪ free on that side. This is strict for every step seen here: EAB|EQ+P (2,0)>(0,−2); EA+B|EQ+P; EA|F+P (1,0)>(0,−1); E+S|EP.
3. **Pure conformational steps.** Tiers 2 and 3, unchanged.

The group sort key becomes (from, to, consumed, released, flag).

**Parameters.** No new Parameter types are needed.
- Pure bindings keep Kd / kon / koff, rendered as today.
- Every other step uses Kiso / kf / kr, rendered by form pair: `k_EAB_to_EQ`, `k_EQ_to_EAB`, `Kiso_EAB_to_EQ`.
- **RE weight rule:** w(to)/w(from) = K^σ · Π[consumed] / Π[released]. Here σ = −1 for a binding Kd and +1 otherwise, and the thermodynamic column signs match.
- Name uniqueness needs one assertion: at most one non-binding step per form pair.
- Everything still flows through `name(p, m)`, so the chokepoint guard test is unaffected.

**Merged central complex.** Keep the substrate-side Species (the canonical `from` of the removed chemistry step) and attach fused steps to it: MERGE(EAB, EPQ) keeps EAB, and EPQ ⇌ EQ + P becomes EAB → EQ + P. No new Species kind is needed. Why:
- this representation already derives exactly (B1–B3b, RB1/3/4, P1–P3, D1–D5);
- the product-side alternative also derives (B6), so the choice is only a dedup convention;
- a display alias such as "EAB=EPQ" is optional (YAGNI).

**DSL.** `_split_side` returns every metabolite term. `consumed` = LHS metabolites and `released` = RHS metabolites, with no inference and no silent drop.

**Impact map** (rough lines changed):

| function / area | change | size |
|---|---|---|
| `Step`, constructor, predicates (types.jl:136-190) | new fields, classification, orientation | +40 / −15 |
| `_step_sides` (rate_eq:159-186) | **deleted**; callers read fields | −28 |
| `_compute_alpha` (269-315) | two branches → one rule; fixes BoundsError | −15 |
| `_compute_numerator` (518-628) | orient each cut step by its uptake; central forms from every chemistry step; fixes −v and RB2 | ±30 |
| `_thermodynamic_constraints`, `_assemble_constraints`, `_step_parameters`, `_free_enz_set`, `_step_priority`, `_partition_independent_count` | stoichiometry from fields; binding predicate is `ligand` | ±15 |
| `name(p::Kon/Koff)`, `_render_binding`, `_emit_cat_params_for_rep` | dispatch on `ligand` | ±10 |
| `_to_sig` / `_step_from_sig`, `_step_canonical_key`, `_assert_no_re_ss_duplicate`, `_step_lex_key` | tuples of consumed/released | ±15 |
| `_canonical_iso_direction`, `_entry_kind` | Tier 1 applied to chemistry-bearing steps | ±15 |
| DSL `_build_step_expr`, `_split_side` | multi-term sides | ±15 |
| `_assert_step_atom_conserving` | atoms(from)+atoms(C) = atoms(to)+atoms(R) | ±5 |
| `_assert_chemistry_is_iso` + 3 call sites | **deleted** | −18 |
| `_flux_carrying_groups`, `_expand_to_allosteric`, `_partial_onlya_catalysis`, `_expand_change_allo_state`, `_onlya_haldane_violation` | `is_iso` → `is_chemistry` | ±8 |
| `_forms_with_binding_step_native`, `_context_form`, `_context_bipartitions`, `_apply_equivalence_grouping`, dead-end mirror steps, `_assert_mechanism_invariants`, `_drop_unbound_regulators` | read consumed/released | ±30 |
| `_catalytic_topologies`, `_release_products!`, dead-end init (≈17 `Step(` calls) | mechanical | ±20 |
| `_kcat_forward` | guard against an empty candidate list (independent fix) | +3 |
| `show`, `reactions` | direct field reads | −5 |

- **Total:** about 25 functions and ~300 lines, net about −50 in `src`.
- **Tests:** 182 references to `bound_metabolite`, `is_binding`, `is_iso` or `Step(` — test_types 69, test_mechanism_enumeration 86, test_rate_eq_derivation 15, test_identify 7, test_dsl 5.
- **Convenience constructor.** A 4-argument `Step(from, to, ligand_or_nothing, eq)` would keep most tests untouched, but it is a backward-compatibility shim, so it needs Denis's approval.
- **Not included:** the MERGE and ELIM moves themselves.

## 6. Part E: timing and a prototype exact rank

**Prototype** (`TD/modrank.jl`, ~130 lines, no code generation):
1. Take the package's `_raw_symbolic_rate_polys(mech, …)` → N and D over the raw group parameters.
2. Run `_assemble_constraints` with the package's priority pivoting. This writes each dependent parameter as Keq^r · Π θ_j^{c_j}, with rational exponents.
3. Let L be the least common multiple of the exponent denominators. Set θ_j = φ_j^L and Keq = κ^L, which makes every exponent an integer.
4. Pick random φ and κ in GF(2^61−1). At F+5 random concentration points, form the row D·Σ_t x_tj N_t − N·Σ_t x_tj D_t for each free parameter j. This equals N·D times ∂log v/∂log φ_j.
5. Take the rank of these rows over GF(p). Repeat with 2 seeds and keep the maximum.

**Agreement (EXACT on finite sets).** Rank and fitted count equal engine A's (`cv/results/export_a.jsonl`) on:
- R4 levels 0–2: 1819 of 1819;
- R5: 4715 of 4715;
- R6 exported: 6000 of 6000 (ids reproduced by re-enumeration, 30,135 distinct).

Phantom mechanisms: 418, 567 and 906 respectively.

**Timings (warm):**

| operation | median | p90 / max |
|---|---|---|
| `init_mechanisms(R4)` (62 seeds) | 7 ms total | – |
| `expand_mechanisms`, per R4 level-0 parent (21 children median) | 26–30 ms | max 0.14 s |
| `expand_mechanisms`, per R4 level-1 parent (100 sampled) | 23 ms | max 0.15 s |
| R4 levels 0–2 under the 3 catalytic moves (1819 mechanisms) | 0.30 s total | – |
| `_independent_param_count` | 2.8 ms | max 0.10 s |
| **modular exact rank, R4** | **3.5 ms** (of which symbolic N/D 1.8 ms) | 5.9 ms / 112 ms |
| modular exact rank, R5 | 4.5 ms | 9.1 ms / 127 ms |
| modular exact rank, R6 | 6.2 ms | 13.8 ms / 199 ms |
| `compile_mechanism` + `fitted_params` (first call, 40 mechanisms) | 0.19 s | max 0.70 s |
| first `rate_equation` call (code generation) | 0.29 s | max 1.03 s |
| BigFloat-Jacobian rank (earlier export, 250 mechanisms) | 0.57 s | max 2.7 s |

**Conclusion.** An exact rank costs about twice `_independent_param_count`, which enumeration already computes for the split gain test.
- For R4 it adds about 75 ms per expanded parent (21 × 3.5 ms).
- Over all of R4 levels 0–2 it totals 7.7 s, against ~870 s to compile the same 1819 mechanisms.

A rank filter or rank annotation in enumeration is affordable. Caveat: the rank inherits the package's derivation, so the bugs of §3 would appear in it; fix those first.

## 7. Open questions for Denis

1. Canonical orientation for fused steps (substrate side proposed). It changes Haldane pivot names relative to hand-written product-side fixtures.
2. Prefix for RE non-binding constants that carry concentration units: keep `Kiso_`, or add a new prefix?
3. Keep a 4-argument convenience `Step` constructor? It counts as a backward-compatibility shim.
4. Allow two consumed metabolites in one step (E + A + B → EAB), or only the kinds the enumerator will emit?
5. Pre-existing: residual-form names break the printed `v` line; and `_kcat_forward` semantics under substrate inhibition (D3).
6. Should the rank be a filter or only an annotation? This is a policy question for the identifiability tracks.