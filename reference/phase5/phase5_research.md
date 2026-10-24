# Phase 5: law engine, unit rule, printer, key, identify cache, L phantoms, codegen, gates

Input to the phase-5 implementation plan (2026-10-14). It merges law.md, unit.md, phantom.md,
keycost.md and the two verifications, and adds one run (section 0, "law type"). The spec is
the authority; the rulings in section 2 settle what it leaves open, for Denis's review.

Paths are relative to the session scratchpad. Abbreviations:
- R5 = `phase45/research5/`. Reports: LAW, UNIT, PH, KC, VI, VC = R5 `law.md`, `unit.md`,
  `phantom.md`, `keycost.md`, `verify_independent.md`, `verify_claim.md`. Logs: LL, UL, PL,
  KL, VL = R5 `law/logs/`, `unit/logs/`, `phantom/logs/`, `keycost/logs/`, `verify/logs/`;
  SL5 = R5 `synth/logs/` (this study: `synth/mkF.py`, `agree_f.jl`, `scan_f.jl`, `run_f.sh`).
- SPEC = `phase3/proto3_wt/docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`;
  S3 = `phase3/proto/S3_design.md` (types §1, read API §7); P4R =
  `phase45/research/phase4_research.md`; EM, KM, SUM = `phase45/understand/{engine,key,SUMMARY}.md`;
  R# = `phase3/RULINGS.md`; P4-# = `phase45/RULINGS.md`.
- WT3 = `phase3/proto3_wt` at 1bf2fee (read only). GM = `graph_mechanism_wt` at 7205ba4 (src =
  4f9a34a). Tests are named by testset, with GM's line, because phase 3 moves lines. TR, TME,
  TID, MD = GM `test/` test_rate_eq_derivation.jl, test_mechanism_enumeration.jl,
  test_identify_rate_equation.jl, mechanism_definitions_for_test_enzyme_derivation.jl.
- SLC = `slice/src/`. Distinct counts replace phantom.md's totals (VI m1, VC M2).

## 0. Answer

**What phase 5 builds** (SPEC:469-575, 650-660; files per SPEC:113-117)
- `poly.jl`: exact sparse polynomials, checked `Int64`.
- `engine.jl`: the law of one conformation on a mask of phase 3's graph; the MWC law on phase
  4's two masks; `complexity` (V×τ) and the 5,000 cap summed over both conformations.
- `thermo.jl`: the solve returns a structured output; `_dependent_param_exprs` renders it.
- `normalize.jl`, `print.jl`: the unit rule; `rate_equation_string` in independents with
  reference lines, `N =`, the coefficient form, per-conformation MWC lines.
- `key.jl`: the canonical key and the L-phantom certificate.
- identify: fit once per key, rows read through their own labelling, `_offer_cv!`, beam and
  CV per key, phantoms inherit, CSV `key`, `fit_source_id`, `L_phantom`; the complexity filter
  counts both conformations.
- `codegen.jl`: locals for dependents, segment sums, H and MWC terms; roots once; the
  size-scaled gate re-calibrated.
- Deleted: `POLY`, `sym_det`, `_raw_symbolic_rate_polys`, `_num_den_exprs`, `_mwc_*`,
  `_invert_monomial`, `eq_hash` and `_rate_eq_dedup_key`, Float exponents in `build_power_expr`.

**Law type (P5-1): checked `Int64` only; `Poly`, `Law`, `AlloLaw` are `mutable struct`s with
`const` fields; a flip-only inactive conformation is a law with Q = 0, D = 1.**
- Outputs equal the slice's on 2,030 of 2,030 laws (law parts, printed law, coefficient form,
  key, code body); law time 1.002× median, p90 1.026× (SL5/agree_f.log:6-8).
- Typed GC workload with six merging loops, `--pkgimages=no`, LICM remarks: this layout (F)
  hoists nothing in 184 escape analyses on 1.12.5 and 214 on 1.13.1; three rounds identical
  (SL5/hoistscan.log:13-16; scan_F_11{2,3}.log "round 3: same").
- Positive control in the same workload: the immutable `Poly` (law.md's C) hoists 3 with a bits
  store and aborts the GC on both versions, in the loop merging a fresh `Law` with a `Poly`
  (SL5/hoistscan.log:1-8; scan_C_112.log "GC error", `fresh_active_or_poly`, scan_f.jl:87-95).
- Mutable `Poly` alone (E) is also clean, 0 hoists in 189 and 222 escape analyses
  (SL5/hoistscan.log:9-12), so the mutable `Law` and `AlloLaw` are margin, not need.
- No coefficient leaves `Int64` under the cap: bound max(F·τ, 2τ·Σ|w_k|), τ ≤ 5,000
  (LAW §1.4); largest coefficient 36 in 2,642 laws; `BigInt` needed 0 times (LL/checks.log:13-14).

**Unit rule (P5-9 to P5-13): constant term; else [Cleland's N, if Denis restores it]; else the
graded basis with fractional λ; else the ratio form.**
- 0 unit failures in 9,406 laws and 4 fixtures; printed laws equal the law to ≤ 3.1e-152
  (UNIT §0; UL/summary.log:22-31).
- The graded basis differs from Cleland's N in 264 of 382 eligible conformations: 3 of 5 specs,
  4 of 85 ping-pong, 257 of 292 ter-ter (UNIT §2). A fresh ter-ter draw: 23 of 28 (VI §4).
- No enumerated law needs fractional λ or the ratio form; four package-accepted fixtures do
  (UNIT §2 F1-F4).

**Phantom rule (P5-30 to P5-36): an exact rescaling certificate in `key.jl`; the phantom
inherits the fit of whoever holds its L→0 key; if none does, it is fitted as an ordinary
candidate.**
- The spec's rule (pre-filter + rank) flags 738 distinct bi-bi variants, 188 in round 1 and 550
  in round 2 (PL/summary.log:3-4). Every one has key(active law) = key(L→0 law) and
  np = np₀ + 1 (PL/summary.log:23). The certificate proves 735 (PL/summary.log:3-4).
- Modular and BigFloat ranks agree on 37,830 of 37,830 ranks (PL/summary.log:14).
- Fits removed, spec rule: 76 of 2,398 (3.2%) and 179 of 7,248 (2.5%); certificate alone: 512
  (21.4%) and 1,386 (19.1%) (PL/fits.log:2-3). A fresh sample: 2.8% and 17.2% (VI §4).
- Sources of the 738: 472 at the enumeration parent, 248 at a relative fitted earlier, 18
  never fitted (chain-flank `re_to_ss`, 8 keys) (PL/summary.log:38-45).
- Rank equality is not family equality: of 60 sampled rank-equal variants without a
  certificate, 3 have a family their L→0 law cannot reproduce numerically (residuals 1.8e-1,
  1.4e-1, 9.8e-3; PL/contain_summary.log:14,15,18; numeric, VC m7).

**Key (P5-18 to P5-21): the slice's canonical labelling with concentrations by name, full 64-hex
SHA-256; at 100 leaves the candidate gets a key of its own.**
- 1,390 ter-ter laws up to 3,438 key terms: 1 leaf each, 0 errors (KL/kt_terter_summary.log:1-9).
- Key time over a round: 1.21% of serial fit time at the 5,000 cap, 0.34% at the default 337
  (KL/key_vs_fit.log:6-7). Slowest key 351 ms (KL/kt_terter_summary.log:11).
- A 16-hex value reads back as a number from 2.1e-3 of one-row CSVs; 64-hex from 0 of 200,000
  (KL/csvkey.log:2,5).

**Gates (P5-22 to P5-25): 0 allocations; 120 ns on the specs; `max(120 ns, c·D)` with D the
denominator after substitution; c by phase 1's rule on CI.**
- Slice codegen: 0 allocations in 578 laws; specs ≤ 13.1 ns (package 41.6); bi-bi gate laws
  42.7 ns (package 61-63); worst 0.108 ns per term on the gate laws (KL/c_calib.log:1,5-11).
- D after substitution equals today's count on 67 of 67 laws, so [396, 396, 396] and
  [589, 799, 858] stay (KL/dterms_vs_pkg.log:1).
- Provisional c = 1.0 ns: local worst 0.108 × the worst CI/local ratio 4.58, doubled (KC §3).

**Complexity filter (P5-26): V×τ(active) + V×τ(inactive) against the threshold, default 337.**
- Census regulated reactions (n = 2): 0 newly filtered in 6-7 rounds (KL/r55_table.log:2-27).
- One-subunit laws: 0.05-0.48% per round sampled uniformly; 10-55% walking largest-first
  (KL/r55_table.log:29-54).

**Still uncertain**
- Ping-pong phantoms are untested: ping-pong round 1 duplicated bi-bi round 1 (VC M2). The plan
  reruns R5.2 on P4-1's ping-pong children (task 10).
- Every number is on slice types. The package port reruns the typed GC workload, key cost and
  gates on its own code (P5-7, tasks 6 and 14).
- c needs CI on three runners; the fixed 120 ns gate already flakes (KM §0 risk 3).
- Five calls wait for Denis (section 2, marked ★): which fit a phantom inherits, where the
  coefficient form lives, the key column's length, the complexity-filter rule, Cleland-first.
- The spec's one minute of population tests is already spent by phases 2 and 4 (31 s + ≤ 30 s,
  P4-28); phase 5 asks for ≤ 20 s more (P5-41).
- No identify run with a real optimizer exercised the key cache, beam per key or phantom CV
  (LAW §7; VC m5). Task 9 adds one.
- R5.5's ping-pong rows used phase 2's tagging; P4-1's half-reaction may raise V×τ(inactive)
  (VC m12). Task 4 reruns the count.

## 1. Decided designs, in pseudo-code over phase 3's graph

### 1.1 Inputs

```
m     = catalytic(am)                 # S3 §1: forms (forms[1] = free E), edges, group, table
e     = m.edges[k]                    # from, to, consumed, released, rapid_equilibrium (S3:101-139)
amask = _active_mask(am)              # phase 4 (P4R §1.7): edges outside :OnlyI groups on forms
                                      #   reached from form 1; plain mechanism: trues
imask = _inactive_mask(am)            # S3:982, pruned from form 1 only (P4R §1.8)
n     = catalytic_multiplicity(am); sites = regulatory_sites(am)
names = metabolites(am)               # S3:932-937: substrates, products, modifiers, site ligands
```
- Sub-graphs are read through masks and `_segments(m, mask, _flags(m))`, never by building a
  mechanism (S3 §7).
- Every derivation entry runs `_check_well_posed` first (R11, P4-23); the engine assumes classes
  1-4 passed.

### 1.2 Structured solve output (`thermo.jl`; P5-4)

```
struct Dependent
    col::Int                                  # the dependent column
    keq::Rational{Int}                        # exponent of Keq
    factors::Vector{Pair{Int, Rational{Int}}} # independent column => exponent
    label::Symbol                             # :kinetic_group, :haldane or :wegscheider
end
struct Solve
    names::Vector{Symbol}                     # every column's name, through _edge_names
    column::Dict{Tuple{Int, Int, Bool}, Int}  # (edge, slot, inactive) => column
    independents::Vector{Int}                 # parameter order
    dependents::Vector{Dependent}
end
_dependent_param_exprs(s) = Dict(s.names[d.col] => power_expr(s, d) for d in s.dependents)
```
- Built by the reduction phase 2 already runs (GM `src/thermo.jl:367-393`), narrowed from
  `Rational{BigInt}` to `Rational{Int}` at the boundary (overflow throws).
- Renders byte-identically on 1,974 of 1,974 mechanisms (LL/q11.log:2-3).
- Consumers: the engine (substitution per column), `all_parameters` (drops the `Expr`
  interpreter `_evaluate_dependent`), codegen's dependent locals, the printer's reference lines.

### 1.3 Polynomials and law types (`poly.jl`, `engine.jl`; P5-1, P5-2)

```
const Exp  = Rational{Int}
const Mono = Vector{Pair{Int, Exp}}           # sorted variable => exponent
mutable struct Poly
    const terms::Dict{Mono, Int}              # Base.checked_add / checked_mul
end
mutable struct Law                            # one conformation: v = E_total·H·Q/D
    const vars::Vars; const H::Poly; const Q::Poly; const D::Poly
    const sigma::Vector{Poly}; const trees::Vector{Poly}; const shift::Mono
end
mutable struct AlloLaw
    const vars::Vars; const H::Poly; const active::Law; const inactive::Law
    const fA::Poly; const fI::Poly
    const sites::Vector{Tuple{Poly, Poly, Int}}   # (W_A base, W_I base, site multiplicity)
    const n::Int
end
flip_only_law(LA) = Law(LA.vars, LA.H, zero(Poly), one(Poly), Poly[], [one(Poly)], Mono())
flip_only(L) = isempty(L.sigma)
```
- No type parameter, no `Union` field, no `BigInt` retry. An `OverflowError` escapes as an
  internal error naming the mechanism.
- `==` and `hash` run over the fields, as `Poly`'s already do (SLC/poly.jl:71-72); a mutable
  struct's default `==` is identity.
- The prototype kept `{T}` with only `Int` instantiated (synth/variantF, mkF.py); the layout is
  the same.
- No phase-5 struct is an immutable wrapper of one reference (the crash layout, LAW §1.2).

### 1.4 Plain engine (`engine.jl`)

```
law(m, s::Solve, mask; inactive = false)::Law
    forest   = RE forest of the masked edges: root = lowest form index, tree edges in edge
               order (free E, form 1, roots its segment)
    alpha[f] = Π equilibrium ratios along f's tree path      # assert: a non-tree RE edge
                                                             # closes with zero net uptake
    one bare variable per (edge k in mask, slot); an RE edge has one slot
    W[(a,b)] += forward flux monomial of SS edge k joining segments a ≠ b; W[(b,a)] += reverse
    T[g]     = Σ over spanning trees of the segment multigraph of Π W, oriented toward g
    sigma[g] = Σ_{f in g} alpha[f]
    den      = Σ_g sigma[g]·T[g]
    num      = Σ_SS k w_k·(T[a]·f_k − T[b]·r_k)    # w_k: net copies of substrates[1] moved,
                                                  # RE-weight exponents included (SPEC:481-486)
    substitute each bare variable: an independent, or Keq^keq·Π independent^e (s.dependents)
    H = Π[S] − Π[P]/Keq;  Q = divexact(num, H)    # throws when H does not divide
    shift = concentration monomial bringing Q, D to lowest terms
    Law(vars, H, Q·shift, D·shift, sigma, T, shift)

complexity(m, mask)::Int = V × τ          # V segments, τ spanning trees with parallel SS
    τ = det(reduced Laplacian) by checked-Int Bareiss; on OverflowError the BigInt det,
        saturated at typemax(Int)
```
- Reads `e.from`, `e.to`, `copies(e.consumed, i)`, `copies(e.released, i)`,
  `e.rapid_equilibrium`, `m.table[i].name` (a `::Mod` copy reads as its reactant).
- The slice's `bare_law` and `law` port line for line (SLC/engine.jl:114-198; LAW §1.5).

### 1.5 Allosteric law (P5-3, P5-8)

```
complexity(am) = complexity(m, amask) + (any(imask) ? complexity(m, imask) : 0)
law(am, s)::AlloLaw
    complexity(am) > MAX_RATE_EQUATION_TERMS && throw(RateEquationTooLargeError(...))
    LA = law(m, s, amask)
    LI = any(imask) ? law(m, s, imask; inactive = true) : flip_only_law(LA)
    flip(L) = L.trees[segment of form 1] · L.shift
    fA, fI  = cross_weights(flip(LA), flip(LI))  # common integer and monomial content out
    per site: W_A = 1 + Σ_{ligand X not :OnlyI} X/K_A(X);  W_I = 1 + Σ_{X not :OnlyA} X/K_I(X)
    # v = E_total·H·(fI^n·QA·DA^(n−1)·WA + L·fA^n·QI·DI^(n−1)·WI)
    #               / (fI^n·DA^n·WA + L·fA^n·DI^n·WI),   WA = Π_sites W_A^multiplicity
class5(law, faces)   # P4R §1.6: plain D, or fI·DA and fA·DI, keeps a term free of the face
```
- The cap and the filter call the same `complexity` (P5-3, P5-26).
- `:OnlyI` groups are outside `amask` (VC M3); the 13 `:OnlyI` fixtures of phase 4 gate this.

### 1.6 Unit rule (`normalize.jl`; P5-9 to P5-13)

```
normalizer(D, vars, substrates)::Norm      # per conformation; Norm is one concrete struct
    Dc = split_conc(D)                     # concentration monomial => coefficient polynomial
    Mono() in keys(Dc)            → Norm(:constant, [Mono() => 1], Mono())
    [P5-10] parts = ΠS/X for each substrate X; all of parts and ΠS in keys(Dc)
                                  → Norm(:cleland, parts .=> 1, ΠS => 1 − |S|)
    order = sort(keys(Dc); by = (degree, distinct non-substrates,
                                 exponent vector over names, larger exponent of an earlier
                                 name first))           # A² before A·B
    basis = first affinely independent [1; a_m] in order
    λ     = solve Σλ = 1, Σλ·a = 0 on basis    (Rational{Int} RREF, fractional allowed)
    consistent                    → Norm(:graded, basis .=> λ, Mono())
    otherwise                     → Norm(:ratio, [order[1] => 1], order[1])
N = Π c(m)^λ_m; integer and monomial content go into the coefficients; a non-monomial
factor prints once as `N = …`
```
- An allosteric law normalizes QA, DA, fA by A's rule and QI, DI, fI by I's (SPEC:512-514).

### 1.7 Printer and coefficient form (`print.jl`; P5-14 to P5-17)

```
rate_equation_string(m; coefficients::Bool = false)::String
  # params and concs lines as today (R37 locals)
  # reference lines: dependent = ratio of products, exact exponents (sqrt, ^(p/q)), sorted
  N = …                                              # only when non-monomial
  v = E_total * (Q) * (H) / (D)    # Q bare when one name, dropped when 1; H is always the
                                   # expression (S - P / Keq), never a name
  coefficients = true:
  v = E_total * n_0 * (H) / (1 + d_A * A + d_A²_B * A^2 * B + …)
  d_A = …                                            # one legend line per coefficient
  ratio form: d_P_B⁻¹ (P/B), n_B⁻¹; `n_0` absent when the numerator has no origin term
  MWC: per conformation s ∈ {A, I}: Q_s, D_s, N_s, F_s (only when F_A ≠ F_I and neither
       is 1), W_s; then the v line with factors of 1 left out; coefficients nA_0, dA_S, fA_…
```
- identify's CSV keeps the full law only (UNIT U6).

### 1.8 Key (`key.jl`; P5-18, P5-19)

```
key_parts(L::Law)     = [H, Q, D];  units [[1], [2, 3]]; header "plain"
key_parts(A::AlloLaw) = [H, QA, DA, fA, QI, DI, fI, (W_A, W_I) per site]
                        units [[1], [2,3,4], [5,6,7], one per site]; header "mwc|n=…|sites=…"
canonical(L)::Canonical            # (order::Vector{Symbol}, text, key)
    each unit × lowest_shift, ÷ integer content
    fixed variables: concentrations labelled c_<name>, Keq, L
    colour refinement on term–constant incidence; individualize the first tied class in
    turn; smallest text; at 100 leaves stop
    text = header|P=…|absent=…\n sorted term lines;  key = bytes2hex(sha256(text))  # 64 hex
leaf cap reached → key = bytes2hex(sha256("mechanism|" * mechanism_id)); log names it
```
- Port of SLC/key.jl:27-118 with `canonical_named`'s one-line change (KC §6).

### 1.9 identify: key cache, beam and CV (P5-27 to P5-29)

```
PASS1 (workers, per fresh mechanism; no compile)
    complexity(m) > eq_complexity_filter            → :complexity_skip
    np = _independent_param_count(m) > max_param_count → :param_skip
    L = law(m, s); class5(L) fails                  → FitFailure
    c = canonical(L); text = rate_equation_string(m, L, s)
    phantom = am with n == 1 && certificate(L)      # 1.10
    c0 = phantom ? canonical(L.active) : c           # the L→0 law's key and order; one
                                                     # concrete type, no Union (R18)
REPS: first member of each key not in memo; a phantom whose c0.key is in memo or among
      this batch's reps inherits; any other phantom is a plain candidate (phantom = false)
PASS2 (pmap over reps): fit compile_mechanism(rep)
    memo[key] = (values in c.order, loss, retcode, peak_at_range_edge,
                 source = rep, source_id = mechanism_id(rep))
ROWS: member  params = c.order[i] => memo[key].values[i]
      phantom params = c0.order[i] => memo[c0.key].values[i], plus L => 0.0
      fit_inherited = row is not the fit's source; fit_source_id = memo[...].source_id
      (a representative's own id); peak_at_range_edge from the source; L_phantom
BatchEntry gains key and fit_key (= key, or c0.key for a phantom)
_offer_cv!: one slot per key (the same code, `key` for `eq_hash`)
beam: _select_count! over distinct keys at a count (budget counts keys);
      selected = Set(keys): a late member of a selected key expands without budget;
      rejected = Set(keys): a late member of a rejected key stays out
CV: fold_losses[fit_key] computed once, by fitting memo[fit_key].source per fold;
    a phantom reads its source's fold losses
```
- Today's code: WT3 `src/identify_rate_equation.jl:304-321, 326-333, 488-585, 600-633,
  725-845, 871-924`.

### 1.10 L phantoms and CSV columns (P5-30 to P5-36)

```
certificate(A::AlloLaw)::Bool                  # phantom.md's rule, without the pre-filter
    A.n == 1 && iszero(QI) || return false
    N = H·fI·QA·WA;  Y = fI·DA·WA;  X = fA·DI·WI
    X ⊆ Y term by term, equal coefficients || return false
    rows: e(m)·w − s = 0 for every m in N and in Y − X;  e(m)·w − s = 1 for m in X
          (e over the fitted constants; concentrations, Keq and L weigh 0)
    return the system is consistent (Rational{Int}; overflow → false, so it is fitted)
# then v(θ, L) = v₀(θ ⊙ (1+L)^w) for all positive θ, L: the families are equal
```
- CSV: `key` (64 hex), `fit_inherited`, `fit_source_id` (16-hex `mechanism_id`, R16),
  `L_phantom`, `peak_at_range_edge`; `eq_hash` goes. A failure row has `missing` in each.
- A test-visible `_phantom_weights(A)::Vector{Rational{Int}}` (empty when none) backs
  `certificate`, so the oracle test can check the map (section 3.1).

### 1.11 Complexity filter (P5-26)

```
_eq_complexity(m) = complexity(m)          # plain: V×τ; allosteric: active + inactive
PASS1: complexity(m) > eq_complexity_filter → :complexity_skip     (default 337)
```
- The docstring states that a one-subunit variant of the 336-term random-order bi-bi counts
  about twice that (KC F1; task 4 checks the value).

### 1.12 Codegen locals and gates (`codegen.jl`; P5-21 to P5-25)

```
program(m, s, L::Law): per conformation, dependent locals numbered by solve column (P5-21),
    single symbols read directly; sigma_g and T_g locals per segment; D = shift·Σ σ_g·T_g; Q
program(am, s, A): both conformations, plus fA, fI, xA = fI·DA, xI = fA·DI, site bases;
    num = H·(fI·QA·xA^(n−1)·WA + L·fA·QI·xI^(n−1)·WI), den = xA^n·WA + L·xI^n·WI
body: inputs by getfield(concs, :A) into positional locals x<i>; x<i>m1 reciprocals and
    roots once (sqrt for q = 2, x^(one(x)/q) for q > 2); greedy Horner with binary + and *;
    E_total * num / den
```
- Gates: `allocs == 0`; `t < 120e-9` on `MECHANISM_TEST_SPECS`; `t < max(120e-9, c·D)` on
  the six gate laws with D = `nterms(law.D)` (MWC: DA + DI, a flip-only DI counts 1).
- c stays 0.8 ns until the CI runs of task 14 set it by phase 1's rule.

## 2. Proposed rulings (ruling — why — cost if wrong)

★ marks the five Denis should see first. ◆ marks other changes to the spec or to CLAUDE.md
that need his consent.

**Law, engine, solve**

P5-1 ◆ Law type F (§1.3): checked `Int64`, no `BigInt` retry; `Poly`, `Law`, `AlloLaw` mutable
  with `const` fields, no type parameter; a flip-only inactive conformation stored as
  Q = 0, D = 1. Amends SPEC:477 ("retried in BigInt") to "checked Int64; overflow is an
  internal error". — Equal outputs on 2,030 laws at 1.002× (SL5/agree_f.log:6-8); 0 hoists on
  1.12 and 1.13 where the immutable `Poly` crashes (SL5/hoistscan.log); overflow needs about
  10¹³ steps (LAW §1.4). Settles EM Q1, VI B1, VC M1. — A law past 2⁶³ would error; if Denis
  keeps the retry, widen at the boundary (LAW variant A, 3× live bytes).
P5-2 The overflow note in the PR covers products as well as laws: codegen's fI·DA, fA·DI and
  the certificate's products stay below τ_I·F·τ_A ≈ 4e8; the cap is checked at the law entry
  before any product. — VC m1. — None; it is a proof note and an ordering.
P5-3 `complexity` returns `Int`: checked Bareiss, the `BigInt` det only after an overflow,
  saturated at `typemax(Int)`; allosteric = V×τ(amask) + V×τ(imask), a flip-only inactive
  conformation counting 0. — Equal to the `BigInt` det on 2,647 pairs; 0.4 vs 0.5 µs and no live
  bytes (LL/checks.log:2-3,11-12); a concrete return type (VC m2); the active mask (VC M3).
  Settles EM Q10. — About 25 lines.
P5-4 The solve returns `Solve` (§1.2); `_dependent_param_exprs` renders it. A fixture with a
  half-integer dependent ("rate_equation_string runs a half-integer dependent as written", GM
  TR:1447) checks the `Rational` path. — Byte-identical on 1,974 mechanisms, all with integer
  exponents (LL/q11.log:2-3; VC m3). Settles EM Q11. — One struct and a renderer.
P5-5 RE turnover stays with phase 4's class 2 at the derivation entry; the engine asserts, while
  it builds the forest, that each non-tree RE edge closes with zero net uptake. — The slice
  engine returns a finite wrong law on both RE-turnover hand cases (LL/checks.log:20-21); 4.6 µs
  standalone. Settles EM Q9. — About 10 lines no public path reaches; tested on the internal
  function, as P4-15's class 5.
P5-6 Class 5 ports onto the new law as P4R §1.6 states. — Rejects exactly the two face cases;
  passes all 42 specs; 0.4 µs (LL/checks.log:22-31,82). Settles EM Q14. — None found.
P5-7 The PR's GC gate is a script outside the suite (R18) in `synth/scan_f.jl`'s pattern, run
  on the package: typed loops over a `Vector{Union{Mechanism, AllostericMechanism}}`, the six
  merging loops with full GCs inside, the derivation and an identify workload,
  `--pkgimages=no`, 1.12 and 1.13. Positive controls: the same loops on an immutable one-field
  wrapper defined in the script, and `repro_upstream.jl`. No phase-5 struct is an immutable
  wrapper of one reference. — An untyped workload reported 0 for code that crashes (VI M1);
  this one separates C from E and F on both versions (SL5/hoistscan.log). Settles VI M1, L7. —
  About 2 minutes per version.
P5-8 Every MWC use reads `_active_mask(am)`: the law, `complexity`, the filter, codegen and the
  key's parts. Phase 4's 13 `:OnlyI` fixtures (P4R §1.7) gate tasks 3-4. — P4-16, P4-17; the
  reports used `trues` (VC M3; EM E3, E6). — An `:OnlyI` mechanism gets a wrong active law.

**Unit rule and printing**

P5-9 The unit rule ships as SPEC:490-514, per conformation for MWC laws. Tests add allosteric
  fixtures on the graded and ratio rules. — 0 unit failures in 9,406 laws (UNIT U1); the
  populations reach graded MWC conformations only in ping-pong n = 2 (VI m2). — None found.
P5-10 ★ Restore Cleland-first: when ΠS/X for every substrate X and ΠS are all denominator terms,
  N is Cleland's. — The graded basis moves N in 3 of 5 specs and 257 of 292 eligible ter-ter
  conformations (13.8% of ter-ter laws) (UNIT U2; VC m13: the ter-ter count, not SPEC:38-39, is
  the evidence); Cleland's N is the textbook Ka·Kb form; 4 lines (SLC/normalize.jl:20-23); values
  never change. Reverses ledger item 11 (EM Q4). — If Denis declines, N stays graded and U3's
  "fewest copies" tie-break fixes the 4 ping-pong power cases.
P5-11 Term order: degree, then distinct non-substrates, then the exponent vector over
  `metabolites(m)` with the earlier name's larger exponent first (A² before A·B). — Other
  tie-breaks move N in at most 17 of 432 ter-ter graded conformations and none elsewhere
  (UL/why.log:7-9); VC m14 asked for the direction. Settles EM Q3. — One sort key.
P5-12 Fractional λ stays; it prints through `power_text` (`sqrt`, `^(3/2)`); fixture F2. —
  A package-accepted mechanism needs it (UL/pkgfix.log:10-15). Settles EM Q5. — About 5 lines.
P5-13 Ratio form as UNIT U5 (divide by the first term, monomial included; `d_P_B⁻¹`, `n_B⁻¹`;
  no `n_0` without an origin term). Fixtures F1, F3. — Values and units pass; the slice's
  coefficient fallback fails units (UNIT §2). Settles EM Q6. — Names change on a class no
  enumerated law reaches.
P5-14 ★ The coefficient form lives behind `rate_equation_string(m; coefficients = true)`; the
  default and identify's CSV keep the full law. — It is about as long as the full law (ter-ter
  medians: coefficient form 4,993, full law 4,534 characters, UL/summary.log:168-169), so both
  in the CSV double the cell; a keyword adds no export or docstring page. Settles EM Q2. —
  Moving it later.
P5-15 MWC names as UNIT U7 (`Q_s`, `D_s`, `F_s`, `W_s`, `N_s`; `nA_0`, `dA_S`); H never prints
  as a name. — `N` already means the normalizer; every name holds `_`, which no metabolite may
  (R10); 0 artifacts in 2,778 MWC laws (UNIT U7). Settles EM Q7, E4. — Renames in the golden
  file and `mwc_allostery.md`, which phase 5 rewrites anyway.
P5-16 `v = E_total * (Q) * (H) / (D)`. — The slice's layout prints `* 1 / (…)` in 8 ping-pong
  laws, this one in none (UNIT U8). Settles EM Q12. — Re-pinning the 11 factored pairs again.
P5-17 Reference lines print as ratios of products with exact exponents, `= 1` for a tie to 1,
  sorted by name. — Value-equal (UNIT U9); today's Float exponents (`^ 0.5`) go. Settles EM Q13.
  — Two re-pins (section 3.4 rows 6, 8).

**Key, gates, codegen**

P5-18 The key ports SLC/key.jl with concentrations by name; PASS1 computes it from the concrete
  law after the complexity filter, without compiling. — 0.34-1.21% of a round's fit time
  (KL/key_vs_fit.log:6-7); I vs J must differ (SPEC:542-546). Settles KM §2d, KC F2. — Per
  candidate it reaches 8-22% above V×τ 2,000.
P5-19 Leaf cap 100 (amends the slice's 100,000); at the cap the candidate's key is the SHA-256 of
  `"mechanism|"` and its `mechanism_id`, and the log names it. — No real law needs 2 leaves
  (1,800 of 1,800, KL/compile_and_rows.log:3); the cap bounds leaves, not time, so 1,000 still
  allows about 6 minutes on a large law (VC m11); a leaf count is deterministic, a time budget
  is not. Settles KM item 12. — A law with over 100 symmetric labellings loses dedup: one
  extra fit per twin.
P5-20 ★ The `key` column holds all 64 hex digits. `mechanism_id` and `parent_id` stay 16 (R16);
  the identify docs read the CSVs with `types = String` for the three id columns. — 16-hex
  values read back as numbers from 2.1e-3 of one-row CSVs, 64-hex from none (KL/csvkey.log:2-5);
  `best_equation.csv` is one row. Settles KM item 6, KC F4. — 48 more characters per row;
  Denis may lengthen the ids as well.
P5-21 Second process: one subprocess (as `test_compile_budget.jl` spawns) computes keys and body
  digests of the 42 specs and 6 gate laws and compares them with the parent's; about 5 key
  digests are pinned; dependent locals are numbered by solve column. — Keys agree across
  processes and 1.10/1.12/1.13; bodies differ on 2 of 78 under 1.13 by local numbering
  (KL/xversion.log:1-7). Settles KM item 8. — About 10 s of suite time.
P5-22 Roots: `sqrt` for q = 2, `x^(one(x)/q)` for q > 2. — Same 9.1 ns as the Float exponent and
  exact in BigFloat (KC §4). Settles KM item 11. — None.
P5-23 The size gate counts D after substitution (MWC: DA + DI). — Equals today's count on 67 of
  67, so the pins stay (KL/dterms_vs_pkg.log:1). Settles EM Q15, KM item 7. — A law whose
  operations outgrow D gets a looser bound.
P5-24 ◆ c follows phase 1's rule on the PR's CI (two runs, three runners, twice the worst
  reading, one significant figure up); 0.8 stays until then; provisional 1.0. CLAUDE.md names
  the constant, so Denis approves the edit. — Local worst 0.108 × 4.58 × 2 (KC F8). — Too loose
  misses a regression under 2×; too tight flakes.
P5-25 ◆ `_testhelper_test_rate_equation_performance` calls `rate_equation` through a
  `@noinline` barrier. — Its loop lets the compiler hoist the call: 4 of 42 specs time under
  0.5 ns (KL/c_calib.log:11); reproduced (VL/gate_loop.log:1-2). This is an emulation on slice
  bodies (VC m10); task 6 confirms it on package code before the change. CLAUDE.md names the
  helper. — About 2 lines; readings rise 1-2 ns.
P5-26 ★ The complexity filter compares V×τ(active) + V×τ(inactive) with the threshold; default
  337. — 0 newly filtered in the census regulated reactions; 0.05-0.48% per round uniformly;
  "distinct" equals "sum"; "max" undercounts; retuning to 674 admits up to 4,920 plain laws per
  round (KL/r55_table.log). Per V×τ unit the dropped MWC laws cost 26% more than dropped plain
  laws (VC m9), so the sum understates cost. Settles EM Q8. — A beam walking toward large
  one-subunit laws loses 10-55% of deep rounds; the user raises the keyword.

**identify and L phantoms**

P5-27 A fit is stored in canonical order with its source; each row reads it through its own
  order; `fit_source_id` is the source's `mechanism_id`, a representative's own; a row's
  `peak_at_range_edge` is its source's. — SPEC:550-553; a renamed law has the same peak (KM
  K8). — None found.
P5-28 `_offer_cv!` keeps one slot per key; LOOCV fits each fit key's source once per fold. —
  SPEC:553-554; duplicates fill CV slots today (KM row 7). — None found.
P5-29 Beam per key: `_select_count!` ranks distinct keys and charges the budget per key; selected
  and rejected keys persist across sweeps. — SPEC:557-558; sources precede phantoms (PH §3);
  the monotone cutoff (WT3 src/identify_rate_equation.jl:744-756). Settles KM item 4. — A subtree
  lost to a late member; unmeasured.
P5-30 ◆ The L-phantom test in `src/` is the certificate (§1.10), in `key.jl`; ranks and numeric
  containment are test oracles. Amends SPEC:651-653 ("exact rank test"). — Inheritance needs
  equal families, which the certificate proves and rank equality does not (PL/contain_summary
  .log:14,15,18, numeric per VC m7); deterministic; 25 µs vs 1.3 ms per modular draw
  (PL/timing.log:2-3), both under 1% of a fit (VC m8); about 45 lines against 135. Settles KM
  item 1, K2. — 3 of 738 spec-rule flags are lost (equal only as a union over L) and get fitted.
P5-31 ◆ No pre-filter: every admitted one-subunit variant with a certificate is a phantom.
  Amends SPEC:659-660 ("a variant whose inactive conformation binds something is fitted"). —
  19-21% of bi-bi fits removed instead of 2.5-3.2%, each with a proof (PL/fits.log:2-3); 17.2% vs
  2.8% on a fresh sample (VI §4). If Denis keeps the pre-filter, it is
  `!any(_inactive_mask(am))` with every site ligand `:OnlyA`, and the certificate runs inside it
  (PH 6). Settles KM item 9. — None found; the certificate's map holds to 8.6e-77 on 3,971
  certificates (PL/summary.log:24).
P5-32 ★ A phantom inherits the fit of whoever holds its L→0 key (the key of its active law);
  when nobody does, it is fitted as an ordinary candidate, unflagged (b2). — The enumeration
  parent is right for only 472 of 738 (PL/summary.log:38-45); b2 needs no fit outside the
  candidate set; it costs at most 18 fits (8 keys) per sampled round 2, or 62 (16 keys) under
  P5-31 (PL/extra.log:15,19). Settles KM item 2. — b1 (fit `catalytic(am)` as the key's
  representative) saves at most 10 fits per sampled round.
P5-33 The phantom's names come from the canonical order of its active law; L = 0. — Keys match
  on all 3,102 certificate phantoms (PH 3; PL/extra.log:2-3). Settles KM item 3. — None found.
P5-34 A phantom takes its source's fold losses; with its source in the pool it loses the 1-SE
  tie; chosen otherwise, it reports L = 0 and `fit_source_id`. — Equal families give equal fold
  optima (SPEC:657-658). Settles KM item 5. — A slot at n + 1 holds a phantom.
P5-35 The oracle rank runs in absolute mode. — Relative and absolute flags agree on 18,915
  variants (PL/summary.log:19). Settles KM item 10. — None.
P5-36 The PR body reports distinct counts (738 spec-rule flags, 3,102 certificate phantoms) and
  says ping-pong phantom coverage is zero; task 10 reruns R5.2 on P4-1's ping-pong children. —
  Ping-pong round 1 duplicated bi-bi round 1 (VI m1, VC M2). — A ping-pong-only class unseen.
P5-37 ◆ ROADMAP Idea: apply the chain-flank exclusion at one catalytic subunit, which removes
  the 18 (62) never-fitted phantoms (PH 2). Denis decides whether to add it. — Today
  `_chain_flank_groups` skips allosteric mechanisms (WT3 src/mechanism_enumeration.jl:459-468).
  — None now.

**Process**

P5-38 Phase-5 tests run on Julia 1.10 once before the PR; the site list has an explicit element
  type. — The slice failed 5 of 78 laws on 1.10 (an untyped empty `map` in a parametric
  constructor, KC §3); non-parametric types remove that path (inferred). — One extra run.
P5-39 A test passes `concs` and `params` as NamedTuples in shuffled field order, and a mechanism
  with a metabolite named `x1`, through the generated body: same value, 0 allocations; R37's
  pin test checks that the body binds only positional locals. — SPEC:567-568 has no test today
  (KM K5). — None.
P5-40 ROADMAP: retire the `eq_hash` half of the regrouping entry, the Planned L-phantom entry and
  `eq_complexity_filter`'s known issue; rewrite the exact-rank Idea (the certificate is in
  `src/`; other uses stay ideas); update the live-byte entry (`complexity` no longer allocates
  `BigInt`s); reword "Drop `change_allo_state` relaxations that free no parameter" (the key
  shares their fit, KM K7); retire the Planned graph rewrite (last phase, KM K1). — SPEC:782,
  843-846. — None.
P5-41 ◆ Phase 5's population tests stay within 20 s, reusing phase 2's and phase 4's samples.
  Denis sees the total, about 80 s against the spec's one minute. — Phases 2 and 4 take 31 s and
  ≤ 30 s (P4-28). — Smaller samples catch less.
P5-42 Version stays 0.10.0; draft PR stacked on phase 4's (R19). — Denis's override. — None.
P5-43 ◆ Before the session ends, R5's logs and data and `phase3/gcprobe/` go to the reference
  branch (SUM, "Both"). — Most cited evidence lives only in the scratchpad (EM E5, KM K6). —
  A few MB in a branch never merged.
P5-44 The PR body labels four evidence limits: R5.1's speed gate was judged by a same-run ratio,
  since the 70.7 µs population was lost (VC m4); the identify GC workload fitted only pass 1 with
  a stand-in loss (VC m5); containment was sampled (43 of 926) and numeric (VC m6, m7); ter-ter
  gave 1,865 laws of 2,000 asked (VC m15). — Honest records. — None.

**Coverage map**
- EM: Q1 → P5-1; Q2 → P5-14; Q3 → P5-11; Q4 → P5-10; Q5 → P5-12; Q6 → P5-13; Q7 → P5-15;
  Q8 → P5-26; Q9 → P5-5; Q10 → P5-3; Q11 → P5-4; Q12 → P5-16; Q13 → P5-17; Q14 → P5-6;
  Q15 → P5-23. E1 → P5-42; E2 → §0 file list; E3, E6 → P5-8; E4 → P5-15; E5 → P5-43;
  E7 → research only.
- KM §5: 1 → P5-30; 2 → P5-32; 3 → P5-33; 4 → P5-29; 5 → P5-34; 6 → P5-20; 7 → P5-23 and
  P5-24; 8 → P5-21; 9 → P5-31; 10 → P5-35; 11 → P5-22; 12 → P5-19. K1 → P5-40, P5-42;
  K2 → P5-30; K3, K4 → PH §2-3 (settled); K5 → P5-39; K6 → P5-43; K7 → P5-40; K8 → P5-27.
- VI: B1 → P5-1; M1 → P5-7; m1 → P5-36; m2 → P5-9.
- VC: M1 → P5-1 (the layout is now built and scanned); M2 → P5-36; M3 → P5-8; minors 1 → P5-2;
  2 → P5-3; 3 → P5-4; 4-7, 15 → P5-44; 8 → P5-30; 9 → P5-26; 10 → P5-25; 11 → P5-19;
  12 → task 4; 13 → P5-10; 14 → P5-11.

## 3. Test plan

### 3.1 Oracle helpers (shared fixtures file MD; every name `_testhelper_`)

| Helper | Source | Role and tolerance |
|---|---|---|
| `_testhelper_steady_state_rate` (exists, MD:301-358) | phase 2 | exact mass-action rate, 256-bit; the law agrees to 1e-20 relative (SPEC:681-682) |
| `_testhelper_law_value` (rebase, MD:367-385) | new `law` | evaluates `Law`/`AlloLaw` in BigFloat by name; replaces the `_num_den_exprs` reader |
| `_testhelper_check_law_against_oracle` (rebase, MD:390-414) | as is | law 1e-20; compiled body with Float64 inputs 1e-8; BigFloat inputs through the body 1e-20 |
| `_testhelper_detailed_balance` (rebase, MD:430-465) | as is | net rate 0 at equilibrium concentrations, per conformation and combined, exact in BigFloat |
| `_testhelper_eval_rate_string` (tighten, GM TR:402-406) | slice check | parses the printed law, full and coefficient forms and reference lines, evaluates in 256-bit; equals the oracle to 1e-20 (today 1e-10) |
| `_testhelper_term_units(m, term)` (new) | UNIT §4, `unit/ulib.jl:275-295` | units from the graph: a rate constant is conc^(−copies taken up)/time, an RE binding K a conc, another RE K conc^(released − consumed); every printed denominator term (0, 0), v in E_total/time; exact |
| `_testhelper_print_artifacts(s)` (new) | `unit/ulib.jl:649-653` | no `^1`, `1 ^`, `* 1`, `L * 1`, `(1)`; exact |
| `_testhelper_same_law_under_map(L1, L2, order1, order2)` (new) | KM §4 row 3 | both laws at 3 random points with constants mapped by canonical position, 1e-20; a swapped map must fail |
| `_testhelper_identifiable_rank` (exists, TME:71-99; moves to MD when a second file uses it) | as is | rank of the variant and of its L→0 law on phantom fixtures |
| `_testhelper_phantom_map_holds(am)` (new) | PL/summary.log:24 | v_oracle(am, θ, L) = v_oracle(catalytic(am), θ ⊙ (1+L)^w) at 3 points per L ∈ {0.1, 10}, 1e-20, with the exact solver on both sides |
| `_testhelper_declared_dimension` (exists, phase 2/4) | as is | `length(fitted_params)` pins |
| phase 4's RE-speed, ε-path and protomer oracles (P4R §3.1) | as is | the new law must keep phase 4's verdicts and values |

### 3.2 Fixtures and expected values (inline DSL; CLAUDE.md's enumeration rules where they apply)

| Fixture | Expected |
|---|---|
| Michaelis–Menten (SS) | constant rule; `v = E_total * (…) * (S - P / Keq) / (1 + …)`; value = oracle |
| Segel Ping Pong Bi Bi spec | graded λ = A:1, B:1, A·B:−1 under both rules; N = Cleland's (UNIT §2) |
| Segel Bi Uni Uni Uni PP Ter Bi spec | Cleland rule N = graded N / (k_FC_to_E_Q + k_EAB_to_F_P)² (UL/examples.log:21,25), names as phase 3 renders them |
| Ping-pong power case 6e46f55b0b68e5d9 (written inline) | graded N = Cleland's squared without P5-10; equal with it (UL/examples.log:36,40) |
| F1 `E + A <--> F + P`, `F + A <--> E + P` | ratio rule, divide by A; `n_P_A⁻¹`, `d_P_A⁻¹`; N = k_E_A_to_F_P + k_F_A_to_E_P (UL/fixtures.log:9-16) |
| F2 F1 + `F + A + A ⇌ F(A, A)` | fractional λ A:3/2, A³:−1/2; N = (…)^(3/2); `sqrt` in coefficients (UL/fixtures.log:20-37) |
| F3 RE Theorell–Chance `E(A) + B ⇌ E(Q) + P` | ratio rule, divide by B; `n_B⁻¹ = k_E_A_to_EA`; no N line (UL/fixtures.log:48-51) |
| F4 F1 + `F(A)`, `F(A, A)` | integer λ A:2, A²:−1; N = (…)² (UL/fixtures.log:64) |
| An MWC fixture per rule (n = 2 ping-pong graded; F3 tagged) | each conformation by its own rule; MWC ratios unchanged (VI m2) |
| Half-integer dependent (GM TR:1447 testset's mechanism) | reference line prints `sqrt(…)`; `Solve` keq 1//2 renders the same `Expr` as today |
| Cube-root fixture (slice `cube_root`) | body emits `x^(one(x)/3)`; BigFloat inputs agree to 1e-20 |
| `re_turnover_beside_ss` built through `_raw_mechanism` | the engine's assertion throws; the derivation entry throws class 2 first |
| Key toys (`proto/p8/test_canon.jl:14-94`) | renaming keeps the key and gives the induced map; monomial and integer normalization; symmetric pairs (4 leaves); 6-cycle vs two triangles differ; an extra declared constant changes the key |
| Key attacks (`slice/skeptic/scripts/keytest.jl:36-58`) | I vs J differ; R vs T differ; competitive vs uncompetitive differ |
| Leaf cap: 9 symmetric constants (`keycost/leafcap.jl`) | 9! leaves → the per-mechanism key; logged |
| uni-bi pair #9/#13 (`census/reports/small_report.md:113`) | one key; one fit (stub optimizer counts 1); each row its own names; rates equal to 1e-12 |
| Phantom fixtures (TME:10828-11250 testsets) | verdicts re-derived under P5-30/P5-31 and checked against the rank oracle and `_testhelper_phantom_map_holds` before pinning; the random-order merged variant (rank 6) stays unflagged; `p_relaxed` may flip under P5-31 |
| A `re_to_ss` chain-flank phantom (#3268 lineage, written inline), in a batch without its source | its L→0 key has no fit, so it is fitted unflagged (b2) |
| Filter: one-subunit K-type of the 336-term random-order bi-bi | `complexity` = active + inactive (task 4 pins the value); skipped at 337, kept at 700 |
| Cap: a one-subunit variant of a ter-ter law at V×τ 4,896 (KC §2 "40 largest"; task 4 confirms the sum) | active ≤ 5,000 < sum → `RateEquationTooLargeError`, today's text |
| Shuffled NamedTuples; a metabolite named `x1` | same value, 0 allocations (P5-39) |
| 13 `:OnlyI` fixtures of phase 4 | law vs exact limit 1e-20; detailed balance; declared dimension (P4R §1.7) |

### 3.3 Population tests (≤ 20 s added, P5-41)

| Test | Population | Check | Estimate |
|---|---|---|---|
| Q1 | phase 2's sample (TO:5-19, 120 mechanisms) | law = oracle 1e-20; printed full and coefficient forms = law 1e-20; units; artifacts; Q > 0 | ~4 s |
| Q2 | the same | rewrites, regroupings and renames change no law, fitted name or key | ~2 s |
| Q3 | the same, plus P4's 60 I-binding children | merged keys pass `_testhelper_same_law_under_map`; a swapped map fails | ~2 s |
| Q4 | init bi-bi seeds | 239 key classes (TID "rate-eq dedup-key partition stability") | ~2 s |
| Q5 | 200 one-subunit bi-bi round-1 variants (seeded) | every flagged variant passes the map oracle; every rank-gaining variant is unflagged | ~5 s |
| Q6 | ping-pong-capable bi-bi, n = 1: P4-1 tagging children of 30 seeded parents | the same as Q5 | ~2 s |
| Q7 | subprocess | keys and body digests of 42 specs and 6 gate laws equal the parent's | ~10 s, outside the 20 s (P5-21) |
| JIT | key, printer, certificate | | ~5 s |

### 3.4 Pins that change, and the value check before each re-pin (SPEC:697-699)

No test is deleted; a subject whose code goes is listed for Denis (CLAUDE.md).

| # | Testset (GM line) | Change | Value check first |
|---|---|---|---|
| 1 | 11 `expected_factored_num`/`den` pairs (MD:1532-2173), read by "Rate Equation String" (TR:644, extractor 683-700) | layout `(Q) * (H) / (D)` | `_testhelper_eval_rate_string` = oracle at 1e-20 on each spec |
| 2 | `test/reference/allosteric_golden_reference.txt` (test_allosteric_golden.jl:9-28) | regenerate | each of 11 specs: printed = oracle 1e-20 |
| 3 | "rate_equation_string allosteric byte-identical fixture" (TR:2006) | rewrite | printed = oracle |
| 4 | docstring doctest (GM src/rate_eq_derivation.jl:470-489); `docs/src/deriving/textbooks.md:41-60` | rewrite | the doctest's law evaluated against the oracle in a test |
| 5 | "rate_equation_string prints every tie as a constraint line" (TR:1377) | ratio-of-products lines | `all_parameters` values equal the printed lines |
| 6 | "rate_equation_string runs a half-integer dependent as written" (TR:1447; `0.5` at 1461) | `sqrt` | printed reference line evaluated = `all_parameters` |
| 7 | "rate_equation_string prints flat +/* sums (no nested parens)" (TR:1472; `count('(') <= 3` at 1485) | the new layout's invariant: flat sums inside Q, H, D | printed = oracle |
| 8 | "a constant thermodynamics pins to 1 is the dependent constant 1.0" (TR:1548; 1563) | `= 1` | value |
| 9 | "Allosteric cat_allo_state stays bound to its step" (TR:1100; regex at 1130) | the RHS-symbol regex skips `sqrt` and legend names | printed = oracle |
| 10 | "Allosteric edge cases" (TR:1882; 1958) and test_allosteric_collapse.jl:53-55, 106-107, 157-158 | likely survive (`D_I`, Wegscheider lines) | run; re-pin only after printed = oracle |
| 11 | Internal pins: "the free-enzyme weight is free E's weight…" (TR:934), "rendering helpers" (TR:2362), allosteric_ground_truth.jl:121-124, 779, 1021-1022 | rewrite on `Law`/`AlloLaw` fields | the law's flip weight = the oracle's free-E weight; subjects whose code goes listed for Denis |
| 12 | "_eq_complexity (V×τ term-count estimate)" (TR:1662; `== 9` at 1724) | allosteric = active + inactive | the value recomputed by the test's own Laplacian |
| 13 | TME "AllostericMechanism — the split count of a law over the term cap" (2654) and "…splits by its count" (2765) | children may change with the summed cap | exact children re-derived; each child's V×τ checked |
| 14 | "Rate equation too large error" (TR:1727); TID "derivation guard … → FitFailure" (1150) | type and text stay | — |
| 15 | "rate_equation time per call grows at most in proportion…" (TR:1317) | counter reads `nterms(law.D)`; pins stay; c per P5-24 | D equals today's count before the switch (67 of 67) |
| 16 | "rate_equation polynomial body uses 2-arg +/* calls" (TR:1224), "…body is in Horner form" (TR:1363) | walk every local of the new body | written first (task 6) |
| 17 | "Degenerate constraint handling" (TR:1523, `build_power_expr`) | power-emission tests | BigFloat agreement |
| 18 | TID "_rate_eq_dedup_key" (303), "rate-eq dedup-key partition stability" (866), "results structure" (359), "fit-dedup by eq_hash in _process_batch" (1402; summary strings 1432-1433, 1468-1469), "renaming-dup pair…different eq_hash" (1487) → one key, "_rows_to_dataframe" (229) | `key` for `eq_hash`; new columns | key unit tests and the same-law oracle first |
| 19 | test_compile_budget.jl "compile-budget" (PASS1 workload at 186-196) | key replaces `_rate_eq_dedup_key` (task 8); no compile in PASS1 once the printer reads the concrete law (task 12) | budgets measured before re-pin |
| 20 | TME "parent and children are rate-equation-distinct" (6306), "rate_equation_string emits section labels" (11633) | likely survive | run |

## 4. Task order, gates, timings and risks

Each task leaves the full suite green (about 5 min 40 s with the heap hint). The old engine and
printer stay beside the new code until task 13 deletes them, so each step is small. Julia runs
go one at a time through `phase3/jl.sh`.

| # | Task | Gate | Runs, suite added |
|---|---|---|---|
| 0 | Baseline at the phase-4 tip: law values at fixed points, printed strings, `eq_hash` classes, fitted names for specs and regenerated populations (census tooling, R22) | files written outside the repo | ~5 min |
| 1 | Structured solve output (P5-4) | `_dependent_param_exprs`, `parameters`, `all_parameters` byte-identical on specs and populations; half-integer fixture | +0 s |
| 2 | `poly.jl` (P5-1): `Poly`, `divexact`, `subst`, shifts, checked overflow | unit tests (rational Keq exponents, n_κ = 2, a forced overflow throws `OverflowError`); P5-7's merging loops on `Poly` alone, scanned on 1.12 and 1.13 outside the suite | ~2 min scan |
| 3 | Plain engine (§1.4), not wired; `complexity` (P5-3); RE-forest assertion (P5-5) | = task 0's values at 1e-20; Q > 0; H divides; detailed balance; rewrites change nothing; `complexity` = `_eq_complexity` on plain laws | +2 s |
| 4 | MWC law on both masks, summed cap and filter (§1.5, §1.11; P5-8, P5-26) | allosteric ground truth (1-3 protomers); 13 `:OnlyI` fixtures; cap and filter fixtures; pins 12-13; R5.5 recount on P4-1 tagging (VC m12) for the PR body | +3 s; recount ~10 min |
| 5 | Rebase test helpers and class 5 on the new law (P5-6) | phase 4's verdicts unchanged; RE-speed and ε-path oracles agree with the new law | +0 s |
| 6 | Codegen on the new law (§1.12; P5-21, P5-22, P5-39); confirm P5-25's hoisting on package code, then change the helper (Denis) | Expr-shape tests rewritten first (pin 16); 0 allocations; 120 ns specs; size gate with c = 0.8; compile budget; subprocess bodies equal; package GC gate P5-7 on 1.12 and 1.13 | +10 s subprocess; scan ~4 min |
| 7 | Key (§1.8; P5-18, P5-19) | toys, attacks, leaf cap; Q2-Q4 | +4 s |
| 8 | identify key cache (§1.9; P5-27, P5-28, P5-20): PASS1 derives law and key from the concrete mechanism, memo by key, rows through own order, `fit_source_id`, CV folds per fit key, CSV `key`; the CSV text stays today's printer until task 12 | TID green; pins 18-19; same-key fixture fits once | +2 s |
| 9 | Beam per key (P5-29) | crowding and late-member tests; one seeded identify run (min_beam_width 10) on synthetic data from a known bi-bi mechanism selects a law with that mechanism's key; fits counted against task 0's run | +2 s; run ~10 min |
| 10 | L phantoms (§1.10; P5-30 to P5-34): certificate, inheritance, `L_phantom` | phantom fixtures; Q5, Q6; R5.2 rerun on P4-1's ping-pong children (P5-36) | +7 s; rerun ~20 min |
| 11 | Unit rule (§1.6; P5-9 to P5-13) | fixtures F1-F4, MM, Segel PP, the power case, MWC per rule; Q1's units | +2 s |
| 12 | Printer (§1.7; P5-14 to P5-17): value checks added first, then the switch | pins 1-10, 17, 20 each after its value check; doctests | +2 s |
| 13 | Delete the old engine and printer (POLY, `sym_det`, `_raw_symbolic_rate_polys`, `_num_den_exprs`, `_mwc_*`, `eq_hash`); port remaining internal pins (11) | full suite, Aqua, JET; deleted subjects listed for Denis | — |
| 14 | c on CI (P5-24): two runs on three runners | CI green twice; c set by phase 1's rule; CLAUDE.md's figure (Denis) | CI |
| 15 | Docs (`textbooks.md`, `ping_pong.md`, `mwc_allostery.md`, identify pages, `developer.md`, README, getting_started), ROADMAP (P5-40), PR body (rulings, P5-36, P5-44, subjects), final GC gate, Julia 1.10 run (P5-38) | docs build with doctests | ~10 min |

Changes to SUM's order:
- The key cache (8-10) moves before the unit rule and printer (11-12), so `eq_hash` is gone before
  the printed text changes; otherwise the 239-class pin would follow the new text.
- The complexity filter joins the cap (task 4): both call one `complexity`.
- Test helpers rebase in task 5, after the law exists and before codegen switches
  `rate_equation`.
- The old engine stays until task 13, so tasks 3-12 each switch one consumer.
- Task 9 adds a real identify run; task 10 adds the ping-pong phantom rerun.

**Risks, ranked**
1. `rate_equation` speed and allocation (task 6), the contract CLAUDE.md calls non-negotiable.
   Slice codegen ran 0 allocations on 578 laws and 42.7 ns on the bi-bi gate laws against the
   package's 61-63 (KL/c_calib.log:1,5-7). If the package port allocates or slows, stop and ask
   Denis.
2. Re-pin volume: 20 pin groups (section 3.4). Each value check comes first; a pin that fails its
   value check is a bug, not a re-pin.
3. identify behaviour (tasks 8-10) changes with no real-optimizer evidence; task 9's run is the
   only end-to-end check.
4. Ping-pong phantoms (task 10) may show a class bi-bi lacks.
5. Suite time: about 20 s of population tests plus 10 s of subprocess over phase 4's total.
6. The package GC gate may find loops the slice lacked; P5-7's adversarial loops and the
   one-field-wrapper rule bound it.
7. c on CI may flake (P5-24, P5-25).
8. Phases 3 and 4 still move: `_active_mask`, `_guard` and S3's read API names may change; the
   plan rebases cites on phase 4's tip.
9. Julia 1.10 (P5-38).

**For Denis** (★ first): P5-10 Cleland-first, P5-14 coefficient form, P5-20 key length,
P5-26 filter rule, P5-32 phantom fit source; then P5-1, P5-24, P5-25, P5-30, P5-31, P5-37,
P5-41, P5-43.
