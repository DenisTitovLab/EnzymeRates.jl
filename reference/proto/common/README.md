# Shared probe data: mechanism populations and a numeric oracle

Everything here is plain data plus stdlib code. Nothing outside `export/` loads
EnzymeRates.

| File | What it is |
|---|---|
| `Population.jl` | `module Pop`: the data model (`PMet`, `PSpecies`, `PStep`, `PSite`, `PSample`, `PMech`) and `load_population(name)`. |
| `Oracle.jl` | `module Oracle`: BigFloat steady-state rate of a plain mechanism, consistent constants, helpers. |
| `data/<name>.jls` | One population: a serialized `Vector{Pop.PMech}`. |
| `data/<name>.summary.txt` | Count, families, export time, notes and every stage error of that population. |
| `export/` | The exporter (needs EnzymeRates): `export_pop.jl`, `export_lib.jl`, `merge_parts.jl`, `run_all.sh`. |
| `validate/validate_oracle.jl` | The oracle validation (stdlib only). `validate/census.jl` summarizes every population (`logs/census.log`); `validate/peek.jl` prints records. |
| `logs/` | Export and validation logs. |

## Loading

`Serialization` resolves the structs as `Main.Pop.*`, so include `Population.jl` into
`Main` before you read a file. `Oracle.jl` includes it for you when it is missing.

```julia
using Random
include("/…/proto/common/Oracle.jl")          # defines Main.Pop and Main.Oracle
m = only(filter(m -> m.source == "spec:Segel Ordered Bi Bi", Pop.load_population("specs")))
s = m.samples[1]; vals = Pop.all_values(s)    # one package evaluation: fitted, Keq, E_total, dependents
k = Oracle.constants_from_names(m, vals; rename = Dict(m.pass2_rename))
v = Oracle.steady_state_v(m, k, Pop.concs(s)) # BigFloat
@assert abs(v - s.v) <= 1e-12 * abs(s.v)      # package Reduced rate (E_total = 1)
kc, _, _ = Oracle.consistent_constants(m, Random.Xoshiro(1); Keq = 3.0)
c = Oracle.concs_at_ratio(m, Random.Xoshiro(2), 3.0)
@assert abs(Oracle.steady_state_v(m, kc, c)) < 1e-30   # detailed balance at Γ = Keq
```

Run scripts as `flock -w 14400 …/proto/julia.lock julia --startup-file=no
--heap-size-hint=1500M script.jl`. `Pop.populations()` lists the files present.

## Populations

All populations are deterministic: fixed seeds, the reactions of
`test/test_mechanism_enumeration.jl`, and the package at the oracle worktree
(EnzymeRates 0.9.1, commit ffb36ad).

| Name | Count | Contents |
|---|---|---|
| `specs` | 42 | All `MECHANISM_TEST_SPECS` (31 plain, 11 allosteric); `source = "spec:<name>"`. |
| `init_uniuni` | 1 | `init_mechanisms(uni_uni_rxn)` (S[C] → P[C]). |
| `init_unibi` | 7 | `init_mechanisms(uni_bi_rxn)` (S[AB] → P[A] + Q[B]). |
| `init_bibi` | 239 | `init_mechanisms(bi_bi_rxn)` (A[C] + B[N] → P[C] + Q[N]). |
| `init_bibi_pp` | 264 | `init_mechanisms(bi_bi_pp_rxn)` (A[CX] + B[N] → P[C] + Q[NX], ping-pong capable). |
| `expand_bibi` | 700 | One round of `expand_mechanisms(·, bi_bi_pp_rxn)` over 60 parents sampled from `init_bibi_pp`: all 700 distinct children (240 plain, 460 allosteric at multiplicity 1 with `:OnlyA`/`:EqualAI` tags). The 2000 cap was not reached. |
| `allo_uniuni` | 35 | `uni_uni_allo_reg` (S → P, allosteric regulator R, oligomeric state 2): its 1 seed gives 5 allosteric children (`child[i]`), and these give 30 new allosteric grandchildren (`grandchild[i]`). 22 of the 35 have a regulatory site. |
| `allo_bibi` | 1477 | Bi-bi at oligomeric state 2: all 977 `_expand_to_allosteric` children of its 239 seeds (`to_allosteric[i]`), plus 500 of their 4910 distinct `expand_mechanisms` children (`expand[i]`; 187 of the population carry `:NonequalAI`). |
| `handwritten` | 34 | The mechanisms of the BigFloat mass-action testsets of `test/test_rate_eq_derivation.jl` (24 plain, listed below) and the ten tag patterns of `test/test_allosteric_collapse.jl` (10 allosteric). |

`handwritten` holds:
- 18 of the 20 `_testhelper_fused_cases`, the 4 `_testhelper_tc_cases`, and two testset mechanisms: the free-enzyme-weight ping-pong and the "pins to 1" iso mechanism;
- 5 structural duplicates were skipped (`fused_cases[4]` = `[1]`, `[6]` = `[5]`, the fused-release `ss`/`re` = `fused_cases[1]`/`[2]`, and `tc_backwards` = `tc_ss`);
- the collapse mechanisms, built by the file's own `_testhelper_uni`/`_testhelper_ro_bibi` helpers with the tags its testsets pass.

Every mechanism exported with `ok = true` (2799 in total: 806 plain, 1993 allosteric).
The only stage errors are four `_kcat_forward` failures ("no kcat components") on
handwritten mechanisms with a rapid-equilibrium fused release or Theorell–Chance outer
steps (`handwritten/3, 6, 10, 20`).

Each `.summary.txt` lists how the population was built, its family counts and every
stage error.

## Data model

### Structure (`PMech`)

- `id` is `"<population>/<index>"`. `source` says where the mechanism came from.
- `family` holds feature tags:
  - the reaction type, such as `:bi_bi`;
  - `:pingpong` (a residual form), `:theorell_chance`, `:fused_binding`,
    `:multi_metabolite_side`, `:inhibitor`, `:dead_end_leaf` (a form with one step),
    `:conformations` (a form not in conformation `:E`);
  - `:all_ss`, `:all_re` or `:mixed_re_ss`, and `:multi_step_group`;
  - `:tag_OnlyA`, `:tag_EqualAI`, `:tag_NonequalAI` and `:reg_sites`.
- `substrates`, `products`: in the package's order, which is sorted by name. The rate is
  the net consumption of `substrates[1]`.
- `inhibitors`: competitive inhibitors the steps bind. `regulators`: allosteric site
  ligands.
- `metabolites`: the concentration keys, in the order of `metabolites(em)`.
- `species[i]` is a form in the order of `_forms`. `bound` holds `PMet(name, role)`
  sorted as the package sorts. `role` is `:S`, `:P`, `:Inh` (a standalone inhibitor or a
  copy `A::Inh`, which keeps the reactant's name and concentration) or `:Reg`.
  `name` is the package's rendered name.
- `steps[j]` is the j-th step **in flat stored order**: group by group, and within a
  group in the package's stored order.
  - The step runs in the package's **stored direction**: `from` takes up `consumed` and
    gives off `released` to become `to`.
  - Every binding is stored with its metabolite consumed. A "binding" here means one
    metabolite consumed and none released, plain or fused, as `Pop.is_binding` tests.
  - `re` marks a rapid-equilibrium step.
- `groups[g]` lists the flat indices of kinetic group g. `Pop.group_of_step(m)` gives
  the inverse.
- Allosteric mechanisms add:
  - `group_tags[g]`, one of `:OnlyA`, `:EqualAI`, `:NonequalAI`;
  - `cat_multiplicity`;
  - `sites`, each with ligands, multiplicity and per-ligand states.

### Constant conventions

Constants per step follow the order `[K]` for an RE step and `[kf, kr]` for an SS step.

- `kf` multiplies the consumed concentrations and drives `from → to`.
- `kr` multiplies the released concentrations and drives `to → from`.
- An RE **binding** has `K` = dissociation constant, so `[to]/[from] = [M]/K`.
- Any other RE step has `K = [to]·Π[released] / ([from]·Π[consumed])`.

The names come from the package (`name(p, m)` of `_step_constants(s, state)`), and they
are **group-representative** names: every step of a group renders the same names.

- **Plain mechanisms:** `step_names[j]`.
- **Allosteric mechanisms:**
  - `step_names_A[j]` and `step_names_I[j]` render every step against the whole
    mechanism with the package's tags. An `:EqualAI` group keeps the bare shared name in
    both states; other groups render `K_A_…`/`k_A_…` and `K_I_…`/`k_I_…`, `:OnlyA` groups
    included.
  - `inactive_kept[j]`: the step survives `_inactive_groups`, the package's pruning of the
    inactive state.
  - `step_names_I_pkg[j]`: the names the package's derivation actually uses for that step
    in the inactive state. It renders them against the pruned, re-canonicalized graph
    (ROADMAP known issue 1), and it is empty for a pruned step.
  - `inactive_reversed[j]`: the pruned graph stores the step in the opposite direction.

### Package outputs

- `fitted_params`; `full_params` (plain only, `parameters(em, Full)`, ends with
  `:E_total`).
- `dependents` (name ⇒ rendered right-hand side) and `dependent_monomials`.
  - Each right-hand side is a power product; the monomial form is `symbol => exponent`
    with `:Keq` included.
  - **Use the monomials, not the strings, to compute.** Residual names such as
    `k_EA_to_EP_res_+A_-P` do not parse as Julia.
  - The package classifies a dependent as Haldane when its right-hand side mentions
    `Keq`, Wegscheider otherwise.
- `rate_equation_string` (Reduced); `eq_hash` = `_rate_eq_dedup_key` of that string;
  `independent_param_count` = `_independent_param_count`.
- `pass2_rename`: the single-symbol Wegscheider fold `_build_wegscheider_rename_map`
  (absorbed ⇒ target). For allosteric mechanisms the per-state maps are
  `pass2_rename_A` and `pass2_rename_I`.
  - **An absorbed name is neither fitted nor dependent.** Its value is its target's.
- Polynomials (`Pop.PPoly` = `[monomial => Rational coefficient]`, monomial =
  `[symbol => exponent]`, exponents may be negative for parameters):
  - **Plain:** `raw_num/raw_den/raw_dfree` come from
    `_raw_symbolic_rate_polys(mech, _step_parameters(mech), Dict())`, with no rename.
    `pkg_num/pkg_den/pkg_dfree` use the package's default arguments, which apply the
    Pass-2 rename. This is what `rate_equation` evaluates in **both** Full and Reduced
    mode. `pass2_changes_polys` flags a difference.
  - **Allosteric:** `state_polys_A`/`state_polys_I` = `[num, den, dfree]` from
    `_state_rate_polys(am, :A/:I)`.
  - `v = E_total·pkg_num/pkg_den` holds for plain mechanisms. Concentrations appear
    with non-negative exponents.
- `samples`: four `PSample`s.
  - `concs`, and `params` = fitted (log-uniform in [1e-2, 1e2], `L` included), `Keq`
    (log-uniform in [1e-2, 1e2]) and `E_total` = 1.
  - `dep_values`: every dependent evaluated in Float64 from `params`.
  - `v` = `rate_equation(em, concs, params)` (Reduced).
  - Plain only: `full_params` gives random raw group constants (log-uniform, `E_total` =
    1), and `v_full` = `rate_equation(em, concs, full_params, Full)` at the same `concs`.
  - `kcat` = `_kcat_forward(em, params)` for all specs and handwritten mechanisms and
    100 of each enumerated population; otherwise NaN.
- Costs:
  - `ctor_ns`: the median of 20 `_with(m)` rebuilds.
  - `first_call_s`: the first Reduced `rate_equation` call, which includes generation
    and compilation.
  - `eval_ns`: the per-call time in the package's performance-test protocol (minimum
    over 5 × 10,000 calls). `eval_allocs`: the bytes one call allocates.
  - `eq_complexity`: V×τ (`_eq_complexity`).
- `ok`/`error`/`stage_errors`: a stage that threw keeps its message. `ok = false` when
  the derivation or the Reduced evaluation failed. The structure is always present.

## The oracle (`Oracle.jl`)

`steady_state_v(m, k, concs; Λ = 10^40)` solves the full mass-action network of a
**plain** mechanism in 512-bit BigFloat:

- every form is its own state;
- an SS step uses `kf`/`kr`;
- an RE step is a fast pair at scale `Λ`: `Λ` and `Λ·K` for a binding, `Λ` and `Λ/K`
  otherwise.

It replaces one row of the generator with the normalization Σx = 1 and solves. It
returns the net consumption of `substrates[1]` per unit enzyme, counted over every step
by name (an `A::Inh` copy carries no net flux at steady state).

`k[j]` must follow `m.steps[j]`: `[K]` for RE and `[kf, kr]` for SS.

Helpers:
- `consistent_constants(m, rng; Keq)` returns bare per-step constants (no kinetic
  groups) that satisfy every Haldane and Wegscheider relation. They come from random form
  free energies (uniform in [-3, 3]), metabolite potentials (products shifted so
  `Σμ_S − Σμ_P = log Keq`, inhibitor copies sharing their reactant's potential) and
  log-uniform `kf` prefactors. It returns `(k, G, mu)`.
- `concs_at_ratio(m, rng, Γ)` returns random concentrations with `Π[P]/Π[S] = Γ`,
  adjusted on `products[1]`.
- `group_constants(m, rng)` returns one random value per distinct name in `step_names`,
  shared across each group. These are the Full-mode semantics without a Pass-2 rename.
- `constants_from_names(m, values; rename, names = m.step_names)` builds per-step
  constants from named values. Pass `names = m.step_names_A` for an allosteric active
  state.
- `eval_poly(p, values)` evaluates a polynomial in BigFloat.
- `steady_state_populations` returns the form populations.

`Pop.re_cycle_closure(m)` (plain mechanisms) returns the number of independent cycles
among the RE steps and how many of them the kinetic grouping leaves **open**, meaning
not closed identically by shared names. With an open RE cycle and independent constants,
the RE subsystem has no equilibrium. The package's Full mode then ignores the K of the
cycle-closing edge, since segment weights follow a BFS spanning tree, and the oracle and
Full mode legitimately disagree.

The oracle covers plain mechanisms only. It has no two-conformation (MWC) network.

## Validation (`validate/validate_oracle.jl`, log in `logs/validate_oracle.log`)

The validation covers 255 plain mechanisms: the 31 plain specs, the 24 plain
handwritten mechanisms, and 200 enumerated ones (all of `init_uniuni` and `init_unibi`,
60 from `init_bibi`, 60 from `init_bibi_pp`, 72 plain from `expand_bibi`).

| Check | Result |
|---|---|
| (i) bare consistent constants, `\|v(Γ=Keq)\| / \|v(Γ=Keq/10)\|` | max 7.2e-108; 0 of 255 above 1e-30 |
| (ii) Full mode at per-group random constants, 197 mechanisms with no Pass-2 fold and no RE cycle | oracle vs the exact raw polynomials ≤ 3.2e-35; oracle vs the package's Float64 `v_full` ≤ 7.0e-15; 0 above 1e-12 |
| (ii) the 58 mechanisms with an RE cycle | the 41 whose RE cycles the grouping closes agree (≤ 2.8e-15); all 17 with an open RE cycle disagree (rel 0.08 to 6e3). That includes the 3 specs with a Pass-2 fold (Non-competitive Inhibitor, Non-essential Activator, Activator + Competitive Inhibitor), whose Full-mode polynomials the fold changes. |
| (iii) Reduced mode at the package's own dependents (BigFloat from `dependent_monomials`) | oracle vs exact package polynomials ≤ 1.6e-34; vs the package's Float64 `v` ≤ 2.0e-14; 0 of 255 above 1e-12 |
| (iv) the package's own Reduced constants at Γ = Keq | oracle ≤ 1.2e-107, package numerator ≤ 5.2e-152 relative; 0 above 1e-30 |

The whole validation takes about 2.5 s of computation.
