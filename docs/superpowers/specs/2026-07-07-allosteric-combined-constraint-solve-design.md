# ABOUTME: Design spec — replace the fragmented allosteric dependent-parameter
# ABOUTME: derivation with one combined constraint solve over A- and I-state cycles.

# Allosteric dependent-parameter derivation via one combined constraint solve

## Summary

The allosteric rate-equation derivation solves the A-state and I-state thermodynamic
constraints in two independent passes, resolves the cross-state affinity splits in a
third pass, then reconciles the three results with hand-written guards. The guards are
wrong: they leave parameters in a circular or undefined state, so `fitted_params`
returns the wrong set and the generated `rate_equation` throws `UndefVarError` at call
time. This is the second class of bugs found in this reconciliation code.

Replace all of it with a single linear solve. Assemble the A-state and I-state
thermodynamic constraints into one matrix over state-tagged parameter columns and run
the existing priority-pivoting Gaussian elimination once. A single solve cannot produce
a circular or undefined dependent set, and the cross-state ties the current code builds
by hand fall out of the solve for free.

## Goal and success criterion

The goal is **simplicity and robustness**, and the measure of success is **net code
removed**. This reconciliation code has now produced two separate classes of bugs; each
past fix added another guard to the seam and made the next bug easier to write. The
correct move is not another guard — it is deleting the seams. This change must remove far
more patchwork than it adds: the whole `_split_resolution` nullspace pass, the collapse
mirrors, the `S_I` referenced-symbol gate, the per-state merge, and `#61`'s filter all
collapse into "assemble the constraints, solve once." If the diff is not strongly
net-negative in lines, the rewrite has drifted back toward patching and should be
reconsidered. Correctness comes from the single solve being unable to express the failure
(no circular or undefined dependent is representable), not from new checks layered on top.

## Background: the bug

The v0.1.6 LDH HPC run errored on ~6.6% of children (10,180 rows across 17 iterations),
every one an `UndefVarError` on an `AllostericEnzymeMechanism`. Twenty-three distinct
symbols, three failure shapes, all from one function —
`_dependent_param_exprs(::Type{AllostericEnzymeMechanism})`
(`src/rate_eq_derivation.jl:1555-1610`):

- **Circular, `:EqualAI`-shared** (`koff_Pyruvate_E`/`koff_Pyruvate_ENAD`, `K_*`): the
  I-state solve picks the opposite Wegscheider pivot from the A-state solve on a shared
  symbol, so the merge marks *both* members of one cycle dependent, each defined from the
  other, neither fitted.
- **Circular, `:NonequalAI`** (`K_I_*`, `koff_I_*`): the merge admits an I-state
  dependent whose left-hand side is a free split referenced by a collapse mirror, forming
  a mirror↔dependent cycle.
- **Undefined, `:NonequalAI` steady-state speed** (`kon_I_*`): the `S_I` gate drops a
  mirror-referenced free speed from the independent set, but the mirror is still emitted,
  referencing a symbol that is never defined.

Reproduction (a minimal case, active-state `koff` cycle):

```julia
# derived rate-equation body for the failing mechanism
(; ..., kon_Pyruvate_E, kon_Pyruvate_ENAD, ...) = params   # neither koff is fitted
koff_Pyruvate_E    = koff_Pyruvate_ENAD * kon_Pyruvate_E / kon_Pyruvate_ENAD
koff_Pyruvate_ENAD = koff_Pyruvate_E  / kon_Pyruvate_E * kon_Pyruvate_ENAD   # circular
```

The base non-allosteric mechanism is correct — it emits a single Wegscheider line with
`koff_Pyruvate_ENAD` fitted. The defect is entirely in the allosteric reconciliation.

`#61`'s uniform `p ∉ keys(dep)` filter (line 1610) did not cause this; it *exposed* it by
turning a pre-existing silent over-parameterization into a hard crash. The dependent-set
undercount and the earlier `#58` overcount are two faces of the same fragile merge.

## Root cause

Two coupled linear subsystems cannot be solved independently and glued. The A-state and
I-state share parameters (`:EqualAI` groups) and are tied by cross-state thermodynamic
cycles. Solving each alone lets the two passes make inconsistent pivot choices on shared
columns; no after-the-fact guard can repair that without effectively re-solving. The
current code has three separate solves —
`_state_dependent_exprs(am, :A)`, `_state_dependent_exprs(am, :I)`, and `_split_resolution`
(a nullspace pass over the affinity splits) — plus a merge. Every bug lives in the seams
between them.

## The fix: one combined solve

The non-allosteric kernel `_dependent_param_exprs_kernel`
(`src/thermodynamic_constr_for_rate_eq_derivation.jl:315-433`) is already the "one clean
function" this needs: it builds a rational constraint matrix over the parameter columns,
runs Gaussian elimination with priority pivoting, and emits `(dependent expressions,
independent set)`. The allosteric path simply refuses to use it globally.

Assemble one matrix and call it once:

- **Columns** = the union of state-tagged parameters. An `:EqualAI` group contributes a
  single shared column used by both states. A `:NonequalAI` group contributes distinct
  `K_A`/`K_I` (or `k_A`/`k_I`) columns. `:OnlyA`/`:OnlyI` contribute a single-state column.
  `L` and regulator `Kreg`s are independent columns that no catalytic cycle constrains.
- **Rows** = the A-state thermodynamic cycles over the A-columns *plus* the I-state cycles
  over the I-columns. Both use the same cycle incidence from `_thermodynamic_constraints`;
  the Haldane right-hand side (in `log Keq`) is identical across conformations.

The cross-state ties need no separate construction. For a cycle `c`, the A-row is
`C_c · logθ_A = r_c` and the I-row is `C_c · logθ_I = r_c`; their difference is
`C_c · (logθ_I − logθ_A) = 0` — exactly the affinity-split constraint `_split_resolution`
computes by hand. The combined row space already contains it, so solving `[A-rows; I-rows]`
yields the collapse mirrors as ordinary dependent expressions, in a valid topological
order, with no circularity possible.

## Design

### Priority (canonical choice)

Re-canonicalize cleanly. Keep the existing non-allosteric column priority (internal
isomerization > metabolite step > free-enzyme binding) applied to the state-tagged
columns, with one allosteric rule: **a non-`:EqualAI` state-specific parameter
(`K_A`/`K_I`, `k_A`/`k_I`) is preferred to stay independent**; the cross-state /
thermodynamically-tied quantities are the pivots (dependents). Within each state the
Wegscheider/Haldane rank still forces some native parameters dependent, exactly as
non-allosteric — the "prefer independent" rule only governs which member of a coupled
pair survives.

### Naming

Allosteric equations are derived here, not taken from a textbook, so positional parameter
naming carries no correctness constraint. Name allosteric parameters by ligand
(`K_ATP`-style) and let the clean re-canonicalization choose the natural form. Positional
naming remains load-bearing only for the textbook-derived **non-allosteric** equations,
which this change does not touch.

### Kernel factoring

Split `_dependent_param_exprs_kernel` into:

- `_assemble_constraints(mech, rename; step_params, all_params)` → `(A, rhs, columns,
  priority)` — the matrix build and priority scoring (current lines 324-386).
- `_solve_dependent_set(A, rhs, columns, priority)` → `(dep_exprs, indep)` — the priority
  pivoting and dependent-expression build (current lines 388-433).

The non-allosteric path assembles one mechanism and solves. The allosteric path assembles
the stacked A/I matrix (shared `:EqualAI` columns, state priorities) and calls the same
`_solve_dependent_set`. Whether to keep one function with an optional coupling argument or
two named helpers is an implementation choice to settle against the code while writing the
plan; the requirement is that both paths run the *same* solver.

### Body assignments (second consumer)

`_build_dep_assignments(M_type)` (`:1674-1720`) builds the constraint lines emitted into
the rate-equation body and today reads the same fragmented sources
(`_state_dependent_exprs`, `_collapse_mirror_exprs`, `_i_state_referenced_syms`). It must
consume the combined solve instead, so the param list and the body can never disagree —
the invariant this bug violates. The body builder's remaining job is placement:
topologically order the combined `dep_exprs`, split them into the A-state and I-state
assignment blocks by where each symbol is referenced, and emit only the dependents the
retained body polynomials actually use. That "actually uses" set is computed from the
body Exprs directly, replacing the separate `_i_state_referenced_syms` gate (which cannot
then drift from real usage).

### Deletions

The combined solve subsumes, and this change removes:

- the merge loop and `S_I` handling in `_dependent_param_exprs(::AllostericEnzymeMechanism)`;
- `_collapse_mirror_exprs` and the whole `_split_resolution` / `SplitResolution` pass;
- `_i_state_referenced_syms` (replaced by direct body-symbol usage);
- the per-state *merge* consumers of `_state_dependent_exprs` (the per-state matrix build
  via `_state_mechanism` / `_thermodynamic_constraints` / column naming stays — it feeds
  the assembly);
- `#61`'s line-1610 `p ∉ keys(dep)` filter (the solve produces disjoint sets by
  construction).

### Fractional split coefficients

`_split_resolution` currently `error`s on a non-integer split coefficient (`1//2`, a
multiply-traversed cycle — a separate ~25-row error class in the run). The combined solve
works over `Rational{BigInt}`, so such a coefficient is a normal pivot and renders as a
rational power rather than crashing. **Open question for Denis:** is a fractional-exponent
dependent (an affinity as a geometric mean of others) a valid mechanism we want to keep,
or should it be rejected earlier as degenerate? The rewrite makes it *not crash*; whether
to accept it is a modeling decision, out of scope for the mechanical fix.

## Testing

- **Oracle / golden.** The existing derivation oracle and golden equations in
  `test/test_rate_eq_derivation.jl` guard non-allosteric output byte-for-byte — those must
  stay green (this change does not touch the non-allosteric path beyond the harmless
  kernel factoring). Allosteric goldens are **re-baselined** to the corrected equations
  and ligand-based names.
- **Thermodynamic consistency.** For every allosteric `MECHANISM_TEST_SPECS` variant,
  assert `v = 0` at `Q = Keq` (equilibrium flux). This is the ground-truth correctness
  check independent of naming.
- **Structural invariant (new regression test).** For every allosteric variant, assert the
  emitted dependent-assignment graph is acyclic and every right-hand-side symbol is defined
  — either fitted (destructured) or a preceding dependent. A combined triangular solve
  satisfies this by construction; the test fails loudly on any regression. Assert also
  `indep ∩ keys(dep) == ∅`.
- **Promote the failures.** Add the three LDH failure families (the saved reproducers for
  the `koff`, `K_I`, and `kon_I` shapes) as permanent specs.
- **Performance contract unaffected.** `_dependent_param_exprs` runs at derivation
  (compile) time, not in `rate_equation`; the 0-allocation / sub-120 ns `rate_equation`
  contract is untouched. Confirm derivation-time cost stays negligible (small rational
  matrices).

## Risks and open questions

1. **Cross-state emergence must be verified, not assumed.** The claim that `[A-rows;
   I-rows]` reproduces every collapse mirror is the load-bearing correctness argument. The
   equilibrium-flux test and the allosteric golden re-baseline are the checks; if any
   mechanism's corrected equation fails `v = 0` at `Q = Keq`, the assembly is incomplete.
2. **Priority tuning.** The exact priority scores that keep non-`:EqualAI` state params
   independent while reproducing sensible mirrors need working out against real mechanisms,
   guarded by the acyclicity and equilibrium tests.
3. **Golden / eq_hash churn.** Allosteric param counts change to the *true* identifiable
   dimension (correcting both the Issue-1 undercount and the `#58` overcount). Every
   allosteric golden and `eq_hash` must be re-baselined; this is expected, not a
   regression, but it is real churn to review.
4. **Fractional coefficients** — see above; a modeling decision for Denis.

## Out of scope

The enumeration and beam-search issues from the same run — the serial `expand_mechanisms`
bottleneck, the parameter-count questions, and the futile re-enumeration cycle — are
separate work, tracked apart from this derivation-correctness fix. A pending experiment is
testing whether the apparent parameter-count non-monotonicity that drives the futile cycle
is itself a downstream artifact of this `fitted_params` corruption; if so, this fix also
shrinks that problem, but that is not a dependency of shipping this one.
