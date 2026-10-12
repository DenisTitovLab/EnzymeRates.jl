# Joint Solve (Phase 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the per-conformation Haldane/Wegscheider solve and its second rename pass
with one joint thermodynamic solve over the constants of both conformations, name every
kinetic class after the member the solve keeps, derive the inactive conformation from the
one graph, export `all_parameters`, remove the `Full`/`Reduced` modes, and test every law
against an exact mass-action oracle.

**Architecture:** A new `src/thermo.jl` builds one integer system over bare step
constants (`_Column`): the fundamental cycles of the graph in each conformation, and one
row per further member of each kinetic class. It orders the columns by the spec's fixed
priority (`_priority_key`, C9 orientation for steps no chemistry orients) and reduces the
system exactly (`Rational{Int}`, retried in `Rational{BigInt}` on overflow). Pivots are
dependents. The naming chokepoint names each class after its kept member
(`_kept_members`, E3); the polynomial engine reads `(reaction, groups, names)`, so the
inactive conformation's law is derived on `_inactive_groups(am)` with the full mechanism's
orientation and names. Tests gain a 256-bit BigFloat mass-action oracle that builds its
network from any plain or two-conformation mechanism.

**Tech Stack:** Julia 1.10+ (CI runs the latest release and 1.13), Test, TestEnv for
focused runs, Documenter for the docs build.

**Spec:** `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`, sections
6 (thermodynamics, names, parameter API), 7 (the phase-2 part of the allosteric model),
14 (testing), 15 (phases), 16 (docs and ROADMAP) and 18 items 1, 3, 11 and 12. Read
sections 6 and 7 before starting any task.

## Global Constraints

- Branch `joint-solve` in the worktree
  `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/joint_solve_wt`
  (below: `$WT`), stacked on `fitting-scale` (PR #91, itself stacked on #85–#89). Never
  rebase or rewrite the stack. Never touch the main checkout
  `/home/denis.linux/.julia/dev/EnzymeRates`. Every command below runs from `$WT`.
- Run Julia only through the shared lock, one process at a time (7.7 GB, no swap):
  `timeout -s KILL <seconds> flock -w 14400 /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock julia --project=. --heap-size-hint=2500M <args>`
  (below: `$JL <args>`, with `<seconds>` = 1800 unless a step says otherwise). Never park
  on a background monitor; poll a log in a bounded shell loop.
- Focused run of one test file (1–3 minutes):
  `$JL -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")'`
  It skips Aqua and JET and sees only the shared fixtures file and the named file, so a
  helper used by two files lives in `test/mechanism_definitions_for_test_enzyme_derivation.jl`.
- Full suite (about 7 minutes), before every commit that changes `src/`, logged and
  polled:
  ```bash
  (JULIA_NUM_PRECOMPILE_TASKS=2 timeout -s KILL 3600 flock -w 14400 /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock julia --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])' > /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase2/logs/suite.log 2>&1; echo "exit=$?" >> /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase2/logs/suite.log) > /dev/null 2>&1 &
  for i in $(seq 1 9); do grep -q '^exit=' /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase2/logs/suite.log 2>/dev/null && break; sleep 60; done
  ```
  Repeat the polling loop until the log ends in `exit=…`; then read its `Test Summary`
  and every `Test Failed` / `Error During Test` block. Baseline at `$WT` HEAD e5fea31:
  33,120 tests pass in 6m11s.
- Code style: 92-character lines, 4-space indentation, match the surrounding code. Every
  new file starts with two `# ABOUTME: ` lines. Names say what code does; no "new",
  "old", "legacy", "improved" or history in names or comments. Comments are evergreen.
- Tests: every helper under `test/` carries the `_testhelper_` prefix (only a one-line
  closure local to a testset or helper keeps a bare name). Test output must be pristine.
  Never delete a failing test: a test of removed behavior is rewritten to test the
  behavior that replaces it, keeping its fixture. Enumeration tests keep CLAUDE.md's
  rules: fixtures inline, exact children with counts.
- `rate_equation` stays allocation-free, under 120 ns per call for every mechanism in
  `MECHANISM_TEST_SPECS` and the three largest enumerated bi-bi laws, within the
  size-scaled gate, and `test/test_compile_budget.jl`'s budgets hold.
- Decided behavior (spec section 6; controller rulings): one joint solve; columns ordered
  by block, then kind (isomerization 3, other chemistry 2, plain binding 1), then both
  forms holding a metabolite, then steady state before rapid equilibrium, then the
  constant against the reaction (E7: the reverse of the solve's forward direction, except
  that a product binding's uptake runs against it), then the base-name content key, the
  smaller key kept (C8); steps that orientation tier 3 decides take C9's form order (free
  E first) in the key only; each class is named after its kept member (E3); `:EqualAI`
  regulator ligands are bare `K_<ligand>reg` with no mirror (B10); parameter order stays
  alphabetical (E4); no second rename pass; `Full`/`Reduced` removed (B5);
  `all_parameters(m, params)` exported (spec section 18 item 1).
- A re-recorded pin must be backed by a value check in the same task: the oracle
  (`_testhelper_check_law_against_oracle`), the string-evaluation test, or an analytic
  value. Only names and the choice of dependents change, never a law's value.
- Commits: never `git add -A` without `git status` first; never skip hooks; every
  commit message ends with exactly these two lines (an implementer may name its own
  model in the first):
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`
  `Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE`
- PRs: never merge. A PR description never contains a claude.ai/code/session link; it
  ends with `🤖 Generated with [Claude Code](https://claude.com/claude-code)`. When
  `gh pr edit` fails, use `gh api -X PATCH repos/DenisTitovLab/EnzymeRates.jl/pulls/<n>`.
- Write prose (docstrings, docs, comments, commit and PR text) clearly and concisely:
  active voice, definite statements, no needless words.

## Review Focus

1. A `:NonequalAI` group that lies wholly behind an `:OnlyA` binding: its inactive class
   is outside the law, so no `K_I_…` name may appear in `fitted_params`, and the law must
   still match the oracle (Task 7, island fixture in `test/test_oracles.jl`).
2. A dependent with a half-integer exponent (`k_F_to_E = Keq ^ -0.5 · …`): the compiled
   law, the printed law and `all_parameters` must agree with the oracle (Task 1 oracle
   fixture; Task 3 `all_parameters` test).
3. `all_parameters` called with a `params` that lacks an independent constant or `Keq`:
   a clear error naming what is missing, never a `KeyError` from inside (Task 3).
4. A hand-written step that neither chemistry nor graph context orients, between free E
   and a form whose rendered name sorts first (`E <--> D`): the solve's key reads C9's
   form order, the stored direction and its names stay (Task 6).
5. A plain mechanism or allosteric lift called 100,000 times: no live-byte growth from
   the solve itself (Task 11 measures; the allosteric lift's `:OnlyA` guard still runs in
   `Rational{BigInt}` until phase 4).

---

## File Structure

- `src/thermo.jl` (create in Task 6; Task 7 adds the dependents, Task 8 folds the
  cycle-space helpers of the guard into it): the joint solve — `_Column`, the content
  keys, `_solve_forward`, `_priority_key`, `_inactive_law`, `_fundamental_cycles`,
  `_joint_system`, `_rref!`, `_reduce_joint_system`, `_kept_members`,
  `_column_parameter`, `_dependent_param_exprs`, `_independent_param_count`.
- `src/thermodynamic_constr_for_rate_eq_derivation.jl` (Task 7 deletes the old solve;
  Task 8 moves what remains into `src/thermo.jl` and deletes the file).
- `src/types.jl` (modify): `_entry_sides`, `_orient_tied_steps(...; lead_key)` (Task 6);
  the E3 chokepoint (`_class_sides`, `_constant_sides`), `_assert_bottomed_re_segments`
  (Task 7); mode types removed (Task 5).
- `src/rate_eq_derivation.jl` (modify): `all_parameters` (Task 3); the engine over
  `(reaction, groups, names)` (Task 4); modes removed (Task 5); `_step_names`,
  `_state_groups`, the single-graph inactive law, B10, no rename pass (Task 7).
- `src/EnzymeRates.jl`, `src/mechanism_enumeration.jl`,
  `src/identify_rate_equation.jl` (modify).
- `test/test_oracles.jl` (create in Task 1), `test/test_joint_solve.jl` (create in
  Task 6), `test/mechanism_definitions_for_test_enzyme_derivation.jl` (shared oracle,
  balance, dimension and population helpers), `test/runtests.jl`, and the existing test
  files each task names.
- `docs/src/deriving/{textbooks,thermodynamic_constraints,kinetic_groups,mwc_allostery}.md`,
  `docs/src/developer.md`, `.claude/CLAUDE.md`, `ROADMAP.md`, `Project.toml` (modify).

Line numbers below are those of `$WT` HEAD e5fea31 unless a task says otherwise; earlier
tasks shift them, so find each anchor by the function or testset it names.

---

### Task 1: Exact steady-state oracle

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (two `_testhelper_lift`
  methods; the oracle block before `# ── Mechanism test specifications`)
- Modify: `test/test_rate_eq_derivation.jl` (delete the 512-bit oracle at lines 99–201;
  call the new check; add `_testhelper_test_mass_action_oracle` to every spec)
- Create: `test/test_oracles.jl`
- Modify: `test/runtests.jl` (include `test_oracles.jl` after `test_rate_eq_derivation.jl`)

**Interfaces:**
- Consumes: `EnzymeRates._forms`, `_inactive_groups`, `_kreg_params`, `_reg_K`, `name`,
  `_num_den_exprs`, `_dependent_param_exprs` (today's, unchanged).
- Produces (later tasks rely on these exact names, all in the fixtures file):
  - `_testhelper_lift(m::Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}) = m`
  - `_testhelper_eval_expr(x, env)` — BigFloat evaluation of a law or dependent `Expr`
  - `_testhelper_single_subunit(m)` — an allosteric compiled mechanism rebuilt with
    catalytic multiplicity 1
  - `_testhelper_flip_form(mech)`
  - `_testhelper_consistent_point(m, rng) -> (; vals::Dict{Symbol,BigFloat}, Keq,
    active, inactive, L)`; `active`/`inactive` hold `(step, Ka, kf)` per law step,
    `kf === nothing` at rapid equilibrium
  - `_testhelper_steady_state_rate(m, point, concs::AbstractDict; Λ = big(10)^30) -> BigFloat`
  - `_testhelper_law_value(m, params::AbstractDict, concs::AbstractDict) -> BigFloat`
    (`m` compiled; per enzyme, no `E_total` factor)
  - `_testhelper_check_law_against_oracle(m; n_points = 3, seed = 1, compiled = true,
    rtol = big(10)^-20)`

The oracle builds the mass-action network from the mechanism's own steps and takes every
constant from free energies, so it never reads the package's solve: a wrong dependent
shows as a mismatch. Rapid-equilibrium steps run as fast pairs at `Λ = 1e30`; the error
falls as 1/Λ (measured on today's code: 8e-13 at 1e10, 8e-23 at 1e20, 8e-33 at 1e30), so
1e-20 agreement has eight decades of margin. Measured on today's code: every one of the
42 specs agrees to 3e-28 or better (all-steady-state specs to 1e-74). An allosteric
mechanism is checked with one catalytic subunit, where the MWC law is the exact steady
state; regulators bind each conformation at rapid equilibrium independently of the
catalytic form, so they rescale the flip ratio by the conformations' site factors.

- [ ] **Step 1: Add the shared oracle helpers**

In `test/mechanism_definitions_for_test_enzyme_derivation.jl`, after the two existing
`_testhelper_lift` methods (lines 119–121) add:

```julia
_testhelper_lift(m::Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}) = m
```

Then insert this block immediately before the line
`# ── Mechanism test specifications ───────────────────────────────────────────`:

```julia
# ── Exact steady-state oracle ───────────────────────────────────────────────

"""
Evaluate an expression of the rate law or of a dependent constant in BigFloat: a Symbol
reads `env`, a number becomes a BigFloat, and a call applies the `Base` function it
names to its evaluated arguments.
"""
_testhelper_eval_expr(x::Symbol, env) = env[x]
_testhelper_eval_expr(x::Real, env) = BigFloat(x)
_testhelper_eval_expr(e::Expr, env) =
    getfield(Base, e.args[1])((_testhelper_eval_expr(a, env) for a in e.args[2:end])...)

"""
`m` with one catalytic subunit: an allosteric mechanism rebuilt with catalytic
multiplicity 1, whose MWC law is the exact steady state of one protomer; a plain mechanism
unchanged.
"""
_testhelper_single_subunit(m::EnzymeRates.EnzymeMechanism) = m
function _testhelper_single_subunit(m::EnzymeRates.AllostericEnzymeMechanism)
    am = EnzymeRates.AllostericMechanism(m)
    EnzymeRates.AllostericEnzymeMechanism(EnzymeRates.AllostericMechanism(
        EnzymeRates.reaction(am), EnzymeRates.steps(am), EnzymeRates.cat_allo_states(am),
        1, EnzymeRates.regulatory_sites(am)))
end

"""
The free form at which the two conformations of an allosteric mechanism interconvert in
its law: the free form (no bound metabolite, no residual) with the smallest name, which
roots the rapid-equilibrium segment that holds every free form.
"""
_testhelper_flip_form(mech) =
    argmin(f -> string(EnzymeRates.name(f)),
           [f for f in EnzymeRates._forms(EnzymeRates.steps(mech))
            if isempty(EnzymeRates.bound(f)) && isempty(EnzymeRates.residual(f))])

"""
    _testhelper_consistent_point(m, rng) -> (; vals, Keq, active, inactive, L)

A thermodynamically consistent parameter point of `m`, built from free energies and so
independent of the package's thermodynamic solve. Each enzyme form has a free energy in
each conformation and each metabolite one; a step's association constant is
`exp(g_from + Σ g_consumed − g_to − Σ g_released)`, so every cycle multiplies to `Keq^n`.
The free energies are projected so that the steps of a kinetic group share their
constants (their metabolites cancel) and, in an allosteric mechanism, an `:EqualAI` step
has one constant in both conformations. A steady-state class draws one forward rate.

`vals` maps the name of every step constant of the law (`EnzymeRates.name(p, mech)`; the
inactive conformation's steps are those of `_inactive_groups`), every regulator constant
and `L` to its value. `active[j]` (and `inactive[j]` for the inactive law's steps) holds
a step as `(step, Ka, kf)`, `kf === nothing` for a rapid-equilibrium step.
"""
function _testhelper_consistent_point(m, rng)
    mech = _testhelper_lift(m)
    allo = mech isa EnzymeRates.AllostericMechanism
    groups = EnzymeRates.steps(mech)
    forms = EnzymeRates._forms(groups)
    rxn = EnzymeRates.reaction(mech)
    mets = EnzymeRates._metabolite_names(rxn)
    nf, nconf = length(forms), allo ? 2 : 1
    fi = Dict(f => i for (i, f) in enumerate(forms))
    mi = Dict(x => nconf * nf + i for (i, x) in enumerate(mets))
    n = nconf * nf + length(mets)
    # The form part of step s's free-energy drop in conformation c (1 or 2).
    function drop(c, s)
        v = zeros(BigFloat, n)
        v[(c - 1) * nf + fi[EnzymeRates.from_species(s)]] += 1
        v[(c - 1) * nf + fi[EnzymeRates.to_species(s)]] -= 1
        v
    end
    tag(g) = allo ? EnzymeRates.cat_allo_state(mech, g) : :None
    ties = Vector{BigFloat}[]
    for (g, group) in enumerate(groups), s in group
        s1 = first(group)
        s == s1 || push!(ties, drop(1, s) - drop(1, s1))
        allo || continue
        tag(g) === :NonequalAI && s != s1 && push!(ties, drop(2, s) - drop(2, s1))
        tag(g) === :EqualAI && push!(ties, drop(2, s) - drop(1, s))
    end
    basis = Vector{BigFloat}[]
    for v in ties
        for b in basis; v = v - (b' * v) * b; end
        nv = sqrt(v' * v)
        nv > big(10)^-30 && push!(basis, v / nv)
    end
    u = BigFloat[4 * big(rand(rng)) - 2 for _ in 1:n]
    for b in basis; u = u - (b' * u) * b; end
    gm(x) = u[mi[EnzymeRates.name(x)]]
    Keq = exp(sum(gm, EnzymeRates.substrates(rxn)) - sum(gm, EnzymeRates.products(rxn)))
    vals = Dict{Symbol, BigFloat}()
    kf_of = Dict{Symbol, BigFloat}()
    function conformation_steps(c, law_groups, state)
        out = []
        for (g, group) in enumerate(law_groups), s in group
            Ka = exp(u[(c - 1) * nf + fi[EnzymeRates.from_species(s)]] +
                     sum((gm(x) for x in EnzymeRates.consumed(s)); init = big(0)) -
                     u[(c - 1) * nf + fi[EnzymeRates.to_species(s)]] -
                     sum((gm(x) for x in EnzymeRates.released(s)); init = big(0)))
            st = tag(g) === :EqualAI ? :EqualAI : state
            if EnzymeRates.is_equilibrium(s)
                vals[EnzymeRates.name(EnzymeRates.Krapid(s, st), mech)] =
                    EnzymeRates.is_binding(s) ? 1 / Ka : Ka
                push!(out, (s, Ka, nothing))
            else
                kf = get!(() -> exp(4 * big(rand(rng)) - 2), kf_of,
                          EnzymeRates.name(EnzymeRates.Kfor(s, st), mech))
                vals[EnzymeRates.name(EnzymeRates.Kfor(s, st), mech)] = kf
                vals[EnzymeRates.name(EnzymeRates.Krev(s, st), mech)] = kf / Ka
                push!(out, (s, Ka, kf))
            end
        end
        out
    end
    active = conformation_steps(1, groups, allo ? :A : :None)
    inactive = allo ? conformation_steps(2, EnzymeRates._inactive_groups(mech), :I) : []
    L = big(0)
    if allo
        L = exp(4 * big(rand(rng)) - 2)
        vals[:L] = L
        for p in [EnzymeRates._kreg_params(mech, :A); EnzymeRates._kreg_params(mech, :I)]
            vals[EnzymeRates.name(p, mech)] = exp(4 * big(rand(rng)) - 2)
        end
    end
    (; vals, Keq, active, inactive, L)
end

"""
    _testhelper_steady_state_rate(m, point, concs; Λ = big(10)^30) -> BigFloat

The exact steady-state rate per enzyme of the mass-action network of `m` at the constants
of `point` (`_testhelper_consistent_point`) and the concentrations `concs`, in 256-bit
BigFloat: the net consumption of the first substrate. Each rapid-equilibrium step runs
as a fast pair at rate scale `Λ`, so the steady state approaches the rapid-equilibrium
limit with a relative error of about 1/Λ times the spread of the constants. An allosteric
mechanism (one catalytic subunit) adds the inactive conformation's forms and steps and one
flip between the two copies of the flip form (`_testhelper_flip_form`) at equilibrium
ratio `L` times the ratio of the conformations' regulatory-site factors. Errors when the
network is disconnected, since its steady state is then not unique.
"""
function _testhelper_steady_state_rate(m, point, concs::AbstractDict; Λ = big(10)^30)
    setprecision(BigFloat, 256) do
        mech = _testhelper_lift(m)
        node = Dict{Tuple{Symbol, EnzymeRates.Species}, Int}()
        id(c, f) = get!(node, (c, f), length(node) + 1)
        x1 = EnzymeRates.name(first(EnzymeRates.substrates(EnzymeRates.reaction(mech))))
        edges = Tuple{Int, Int, BigFloat, BigFloat, Int}[]
        for (c, list) in ((:A, point.active), (:I, point.inactive)), (s, Ka, kf) in list
            f, r = kf === nothing ? (Λ, Λ / Ka) : (kf, kf / Ka)
            rf = f * prod((concs[EnzymeRates.name(x)] for x in EnzymeRates.consumed(s));
                          init = big(1))
            rr = r * prod((concs[EnzymeRates.name(x)] for x in EnzymeRates.released(s));
                          init = big(1))
            ω = count(x -> EnzymeRates.name(x) == x1, EnzymeRates.consumed(s)) -
                count(x -> EnzymeRates.name(x) == x1, EnzymeRates.released(s))
            push!(edges, (id(c, EnzymeRates.from_species(s)),
                          id(c, EnzymeRates.to_species(s)), rf, rr, ω))
        end
        if mech isa EnzymeRates.AllostericMechanism
            # Regulators bind each conformation at rapid equilibrium, independently of
            # the catalytic form, so summed over regulator occupancy the flip ratio is
            # L times the ratio of the conformations' regulatory-site factors.
            function site_factor(inactive)
                prod((1 + sum((concs[EnzymeRates.name(lig)] / point.vals[K]
                               for (lig, tag) in zip(EnzymeRates.ligands(site),
                                                     EnzymeRates.allo_states(site))
                               for K in (EnzymeRates._reg_K(mech, site, lig, tag,
                                                            inactive),)
                               if K !== nothing); init = big(0)))^
                     EnzymeRates.multiplicity(site)
                     for site in EnzymeRates.regulatory_sites(mech); init = big(1))
            end
            flip = _testhelper_flip_form(mech)
            push!(edges, (id(:A, flip), id(:I, flip),
                          point.L * site_factor(true) / site_factor(false), big(1), 0))
        end
        n = length(node)
        reached, frontier = Set([1]), [1]
        while !isempty(frontier)
            u = pop!(frontier)
            for (a, b, _, _, _) in edges, (p, q) in ((a, b), (b, a))
                p == u && q ∉ reached && (push!(reached, q); push!(frontier, q))
            end
        end
        length(reached) == n ||
            error("_testhelper_steady_state_rate: the network of $n enzyme forms is " *
                  "disconnected, so its steady state is not unique")
        M = zeros(BigFloat, n, n)
        for (a, b, rf, rr, _) in edges
            M[b, a] += rf; M[a, a] -= rf; M[a, b] += rr; M[b, b] -= rr
        end
        M[n, :] .= 1
        rhs = zeros(BigFloat, n); rhs[n] = 1
        x = M \ rhs
        sum((ω * (rf * x[a] - rr * x[b]) for (a, b, rf, rr, ω) in edges); init = big(0))
    end
end

"""
The rate law of the compiled mechanism `m` per enzyme, evaluated in 256-bit BigFloat from
the polynomials `rate_equation` compiles (`_num_den_exprs`) at the parameter values
`params` (by name: the fitted constants, `Keq` and, for an allosteric mechanism, `L` and
the regulator constants) and the concentrations `concs`. The dependent constants are
computed from `_dependent_param_exprs`, as the compiled body computes them.
"""
function _testhelper_law_value(m, params::AbstractDict, concs::AbstractDict)
    setprecision(BigFloat, 256) do
        num, den = EnzymeRates._num_den_exprs(typeof(m))
        dep, _ = EnzymeRates._dependent_param_exprs(typeof(m))
        env = Dict{Symbol, BigFloat}(k => BigFloat(v) for (k, v) in params)
        for (k, ex) in dep
            env[k] = _testhelper_eval_expr(ex, env)
        end
        for (k, v) in concs
            env[k] = BigFloat(v)
        end
        _testhelper_eval_expr(num, env) / _testhelper_eval_expr(den, env)
    end
end

"""
Check that the rate law of `m` equals the exact steady state of its mass-action network
(`_testhelper_steady_state_rate`) to `rtol` at `n_points` consistent points
(`_testhelper_consistent_point`). The law reads the point's fitted constants and computes
its dependents from them, so a wrong dependent shows as a mismatch. With `compiled`, the
compiled `rate_equation` is checked too, in Float64 to 1e-8. An allosteric mechanism is
checked with one catalytic subunit.
"""
function _testhelper_check_law_against_oracle(m; n_points = 3, seed = 1, compiled = true,
                                              rtol = big(10)^-20)
    rng = Random.MersenneTwister(seed)
    em = _testhelper_single_subunit(m)
    fitted = EnzymeRates._dependent_param_exprs(typeof(em))[2]
    mets = EnzymeRates.metabolites(em)
    for _ in 1:n_points
        point = setprecision(() -> _testhelper_consistent_point(em, rng), BigFloat, 256)
        concs = Dict{Symbol, BigFloat}(x => exp(3 * big(rand(rng)) - big(1.5))
                                       for x in mets)
        params = Dict{Symbol, BigFloat}(k => point.vals[k] for k in fitted)
        params[:Keq] = point.Keq
        oracle = _testhelper_steady_state_rate(em, point, concs)
        @test abs(_testhelper_law_value(em, params, concs) - oracle) <=
              rtol * abs(oracle)
        compiled || continue
        p64 = NamedTuple{(fitted..., :Keq, :E_total)}(
            (Float64.(getindex.(Ref(params), fitted))..., Float64(point.Keq), 1.0))
        c64 = NamedTuple{mets}(Tuple(Float64(concs[x]) for x in mets))
        @test isapprox(rate_equation(em, c64, p64), Float64(oracle); rtol = 1e-8)
    end
end

```

- [ ] **Step 2: Write the oracle's own tests**

Create `test/test_oracles.jl`:

```julia
# ABOUTME: Tests of the exact steady-state oracle, and of the laws against it over
# ABOUTME: regenerated populations: agreement, detailed balance and declared dimension.
using Random

@testset "Exact steady-state oracle" begin
    @testset "matches the steady-state Michaelis–Menten law exactly" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
            end
        end
        point = setprecision(() -> _testhelper_consistent_point(
            m, Random.MersenneTwister(1)), BigFloat, 256)
        # Both steps are stored as bindings: E + S → E(S) and E + P → E(S).
        _, Ka_S, k1 = only(x for x in point.active
                           if EnzymeRates.consumed(x[1]) == [EnzymeRates.Substrate(:S)])
        _, Ka_P, km2 = only(x for x in point.active
                            if EnzymeRates.consumed(x[1]) == [EnzymeRates.Product(:P)])
        km1, k2 = k1 / Ka_S, km2 / Ka_P
        S, P = big"0.7", big"1.3"
        v = setprecision(BigFloat, 256) do
            (k1 * k2 * S - km1 * km2 * P) / (km1 + k2 + k1 * S + km2 * P)
        end
        oracle = _testhelper_steady_state_rate(m, point, Dict(:S => S, :P => P))
        @test abs(oracle - v) <= big(10)^-60 * abs(v)
    end

    @testset "a rapid-equilibrium step's error falls as 1/Λ" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        rng = Random.MersenneTwister(3)
        point = setprecision(() -> _testhelper_consistent_point(m, rng), BigFloat, 256)
        concs = Dict(:S => big"0.7", :P => big"1.3")
        params = Dict{Symbol, BigFloat}(k => point.vals[k]
                                        for k in EnzymeRates.fitted_params(m))
        params[:Keq] = point.Keq
        law = _testhelper_law_value(m, params, concs)
        err(Λ) = abs(_testhelper_steady_state_rate(m, point, concs; Λ) - law) / abs(law)
        e10, e20, e30 = err(big(10)^10), err(big(10)^20), err(big(10)^30)
        @test e20 < 1e-9 * e10
        @test e30 < 1e-9 * e20
        @test e30 < big(10)^-20
    end

    @testset "rejects a disconnected network" begin
        two_enzymes = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
                F + S <--> F(S)
                F(S) <--> F + P
            end
        end
        point = _testhelper_consistent_point(two_enzymes, Random.MersenneTwister(1))
        @test_throws ErrorException("_testhelper_steady_state_rate: the network of 4 " *
                                    "enzyme forms is disconnected, so its steady state " *
                                    "is not unique") _testhelper_steady_state_rate(
            two_enzymes, point, Dict(:S => big(1), :P => big(1)))
    end

    @testset "the two conformations combine at the flip ratio L" begin
        # One catalytic subunit, each conformation a rapid-equilibrium Michaelis–Menten
        # enzyme with its own S binding and chemistry: v = Σ_c w_c·(kf_c·S/KS_c −
        # kr_c·P/KP) / Σ_c w_c·(1 + S/KS_c + P/KP), with w = (1, L).
        m = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: NonequalAI
                E(P) ⇌ E + P      :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(m)
        point = setprecision(() -> _testhelper_consistent_point(
            m, Random.MersenneTwister(5)), BigFloat, 256)
        flat = first.(EnzymeRates._flat_steps(EnzymeRates.steps(am)))
        bind_s = only(s for s in flat if EnzymeRates.consumed(s) ==
                                         [EnzymeRates.Substrate(:S)])
        bind_p = only(s for s in flat if EnzymeRates.consumed(s) ==
                                         [EnzymeRates.Product(:P)])
        chem = only(s for s in flat if EnzymeRates.is_iso(s))
        val(p) = point.vals[EnzymeRates.name(p, am)]
        KP = val(EnzymeRates.Krapid(bind_p, :EqualAI))
        S, P = big"0.7", big"1.3"
        terms = map(((:A, big(1)), (:I, point.L))) do (state, w)
            KS = val(EnzymeRates.Krapid(bind_s, state))
            kf, kr = val(EnzymeRates.Kfor(chem, state)), val(EnzymeRates.Krev(chem, state))
            (w * (kf * S / KS - kr * P / KP), w * (1 + S / KS + P / KP))
        end
        v = setprecision(() -> sum(first, terms) / sum(last, terms), BigFloat, 256)
        oracle = _testhelper_steady_state_rate(m, point, Dict(:S => S, :P => P))
        @test abs(oracle - v) <= big(10)^-25 * abs(v)
    end

    @testset "a dependent with a half-integer exponent matches the oracle" begin
        # The conformational square ties two isomerizations of one group, so the
        # Haldane relation fixes k_F_to_E through square roots (Keq ^ -0.5 · …).
        half = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                (E(S) <--> F(S), E <--> F)
                F(S) <--> E + P
                F + S <--> E(S)
            end
        end
        dep, _ = EnzymeRates._dependent_param_exprs(typeof(half))
        @test any(ex -> occursin("^ -0.5", string(ex)), values(dep))
        _testhelper_check_law_against_oracle(half)
    end
end
```

In `test/runtests.jl`, add `include("test_oracles.jl")` on the line after
`include("test_rate_eq_derivation.jl")`.

- [ ] **Step 3: Run the oracle tests (they pass against today's code)**

Run the focused command for `test_oracles`.
Expected: `Exact steady-state oracle | 13 13` (every test passes). These are oracle checks,
not a red step: the oracle must agree with today's laws before it judges the new ones. If
the half-exponent assertion fails because today's solve picks an integer-exponent
dependent for `half`, report it with the printed `rate_equation_string(half)` instead of
changing the fixture.

- [ ] **Step 4: Replace the 512-bit oracle in the derivation tests**

In `test/test_rate_eq_derivation.jl`:

1. Delete the block from the comment line
   `# ── Brute-force mass-action oracle (no Cha segments) ────────────────────────`
   through the end of `_testhelper_check_against_mass_action` (lines 99–201): the old
   `_testhelper_consistent_point(em, rng)`, `_testhelper_mass_action_rate` and
   `_testhelper_check_against_mass_action`. The shared helpers replace them.
2. Replace every call `_testhelper_check_against_mass_action(em)` (three sites: the
   fused cases, the Theorell–Chance cases and the `K_E_to_F = 1` testset) with
   `_testhelper_check_law_against_oracle(em)`.
3. Immediately before the docstring that starts `Run all tests for a mechanism
   specification.` add:

```julia
"""The spec's rate law against the exact steady state of its mass-action network
(`_testhelper_check_law_against_oracle`)."""
function _testhelper_test_mass_action_oracle(spec::_testhelper_MechanismTestSpec)
    @testset "Mass-action oracle" begin
        _testhelper_check_law_against_oracle(spec.mechanism)
    end
end

```

4. In `_testhelper_run_all_tests`, add the line
   `        _testhelper_test_mass_action_oracle(spec)` right after
   `        _testhelper_test_constraint_counting(spec)`.

- [ ] **Step 5: Run the derivation tests**

Run the focused command for `test_rate_eq_derivation`.
Expected: all pass; each of the 42 spec testsets gains a `Mass-action oracle` testset
with 6 passing tests (3 points × law and compiled check).

- [ ] **Step 6: Run the full suite, then commit**

Run the full suite. Expected: all pass (about 33,400 tests).

```bash
git status
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_oracles.jl \
        test/test_rate_eq_derivation.jl test/runtests.jl
git commit -m "Check every spec's law against an exact mass-action oracle" \
  -m "The oracle builds the steady-state network of any plain or one-subunit allosteric
mechanism in 256-bit BigFloat, runs rapid-equilibrium steps as fast pairs at 1e30 and
takes its constants from free energies, so it never reads the package's solve. It
replaces the plain 512-bit solver, rejects a disconnected network, and agrees with every
spec's law to 1e-20." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---

### Task 2: Detailed balance over populations and the declared dimension

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (balance, dimension
  and population block after the oracle block)
- Modify: `test/test_oracles.jl` (population sample; detailed-balance testset)
- Modify: `test/test_rate_eq_derivation.jl` (`_testhelper_test_constraint_counting`)

**Interfaces:**
- Consumes: Task 1's `_testhelper_law_value`, `_testhelper_eval_expr`, `_testhelper_lift`.
- Produces (fixtures file):
  - `_testhelper_poly_value(p, env)`
  - `_testhelper_detailed_balance(m, rng) -> Vector{BigFloat}` (each entry is zero in
    exact arithmetic)
  - `_testhelper_rank(A::AbstractMatrix{<:Integer}) -> Int`
  - `_testhelper_cycle_basis(groups) -> Vector{Vector{Int}}`
  - `_testhelper_declared_dimension(m) -> Int` (fitted catalytic constants + regulator
    constants + `L`)
  - `_testhelper_population(rxn, rounds) -> Vector{Union{Mechanism, AllostericMechanism}}`
- Produces (`test/test_oracles.jl`): `_testhelper_pingpong_allo_rxn`,
  `_testhelper_population_sample` (120 mechanisms, the same on every run).

`_testhelper_declared_dimension` is spec success criterion 4's independent count: it
builds the classes (a group's constant per slot and conformation; one shared class for
an `:EqualAI` group; one class per `:OnlyA` step's inactive copy) and the cycle rows
itself, from its own cycle basis and exact rank, and keeps
`n_law − (rank(rows) − rank(rows over the classes outside the law))` dimensions. It
equals today's fitted count on all 42 specs. On populations it exposes ROADMAP known
issues 1 and 2 (today 71 of the 2,344 mechanisms in `_testhelper_population(_testhelper_pingpong_allo_rxn, 1)`
miss it), so the population dimension test arrives with the switch (Task 7).

- [ ] **Step 1: Add the shared balance, dimension and population helpers**

In the fixtures file, insert after the oracle block (before
`# ── Mechanism test specifications`):

```julia
# ── Detailed balance, declared dimension and populations ────────────────────

"""The value of the polynomial `p` at the values `env` (by symbol), in BigFloat."""
_testhelper_poly_value(p, env) =
    sum((c * prod((env[s]^e for (s, e) in mono); init = big(1)) for (mono, c) in p);
        init = big(0))

"""
    _testhelper_detailed_balance(m, rng) -> Vector{BigFloat}

Detailed-balance residuals of `m`'s rate law at a random parameter point, in 256-bit
BigFloat: the law's value at equilibrium concentrations (the first product scaled so
that Π[P]/Π[S] = Keq) over its value with that product ten times lower, then the same
ratio for each conformation's numerator that is not identically zero (`_state_rate_polys`;
for a plain mechanism the law's own numerator). Each is zero in exact arithmetic when the
law satisfies every Haldane and Wegscheider relation.
"""
function _testhelper_detailed_balance(m, rng)
    setprecision(BigFloat, 256) do
        em = m isa EnzymeRates.AbstractEnzymeMechanism ? m :
             EnzymeRates.compile_mechanism(m)
        mech = _testhelper_lift(em)
        params = Dict{Symbol, BigFloat}(k => exp(4 * big(rand(rng)) - 2)
                                        for k in EnzymeRates._dependent_param_exprs(typeof(em))[2])
        params[:Keq] = exp(4 * big(rand(rng)) - 2)
        rxn = EnzymeRates.reaction(mech)
        subs = EnzymeRates.name.(EnzymeRates.substrates(rxn))
        prods = EnzymeRates.name.(EnzymeRates.products(rxn))
        concs = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1)
                                       for x in EnzymeRates.metabolites(em))
        p1 = first(prods)
        concs[p1] = params[:Keq] * prod(concs[x] for x in subs) /
                    prod((concs[x] for x in prods if x != p1); init = big(1))
        off = copy(concs)
        off[p1] /= 10
        ratios = BigFloat[abs(_testhelper_law_value(em, params, concs)) /
                          abs(_testhelper_law_value(em, params, off))]
        mech isa EnzymeRates.AllostericMechanism || return ratios
        dep, _ = EnzymeRates._dependent_param_exprs(typeof(em))
        env = copy(params)
        for (k, ex) in dep
            env[k] = _testhelper_eval_expr(ex, env)
        end
        for state in (:A, :I)
            num = first(EnzymeRates._state_rate_polys(mech, state))
            isempty(num) && continue
            push!(ratios, abs(_testhelper_poly_value(num, merge(env, concs))) /
                          abs(_testhelper_poly_value(num, merge(env, off))))
        end
        ratios
    end
end

"""The rank of an integer matrix, by exact elimination over the rationals."""
function _testhelper_rank(A::AbstractMatrix{<:Integer})
    R = Rational{BigInt}.(A)
    rank = 0
    for c in axes(R, 2)
        p = findfirst(i -> !iszero(R[i, c]), (rank + 1):size(R, 1))
        p === nothing && continue
        p += rank
        R[rank + 1, :], R[p, :] = R[p, :], R[rank + 1, :]
        for i in (rank + 2):size(R, 1)
            iszero(R[i, c]) || (R[i, :] .-= (R[i, c] / R[rank + 1, c]) .* R[rank + 1, :])
        end
        rank += 1
    end
    rank
end

"""A basis of the cycle space of the step graph of `groups`: one integer vector over the
flat steps (each counted in its stored direction) per independent cycle, from the null
space of the form-by-step incidence matrix."""
function _testhelper_cycle_basis(groups)
    flat = [s for group in groups for s in group]
    forms = unique(f for s in flat for f in (EnzymeRates.from_species(s),
                                             EnzymeRates.to_species(s)))
    R = zeros(Rational{BigInt}, length(forms), length(flat))
    for (j, s) in enumerate(flat)
        R[findfirst(==(EnzymeRates.from_species(s)), forms), j] -= 1
        R[findfirst(==(EnzymeRates.to_species(s)), forms), j] += 1
    end
    pivots = Int[]
    for c in axes(R, 2)
        r = length(pivots) + 1
        p = findfirst(i -> !iszero(R[i, c]), r:size(R, 1))
        p === nothing && continue
        p += r - 1
        R[r, :], R[p, :] = R[p, :], R[r, :]
        R[r, :] ./= R[r, c]
        for i in axes(R, 1)
            i != r && !iszero(R[i, c]) && (R[i, :] .-= R[i, c] .* R[r, :])
        end
        push!(pivots, c)
    end
    map(setdiff(axes(R, 2), pivots)) do f
        z = zeros(Rational{BigInt}, size(R, 2))
        z[f] = 1
        for (k, c) in enumerate(pivots); z[c] = -R[k, f]; end
        Int.(z .* lcm(denominator.(z)))
    end
end

"""
The dimension of the model `m` declares, counted without the package's thermodynamic
solve: the classes of constants in its law (a kinetic group's constant, slot by slot, in
each conformation the law holds it; an `:EqualAI` group shares one class between the
two), less the rank the cycle rows of both conformations impose on them, plus each
regulator constant and `L`. A cycle row also reads the classes outside the law (the
inactive copy of each `:OnlyA` step, which is a class of its own, and the inactive
classes on pruned forms), which absorb what they can: the law keeps
`n_law − (rank(rows) − rank(rows over the classes outside the law))` dimensions.
"""
function _testhelper_declared_dimension(m)
    mech = _testhelper_lift(m)
    allo = mech isa EnzymeRates.AllostericMechanism
    groups = EnzymeRates.steps(mech)
    flat = [(s, g) for (g, group) in enumerate(groups) for s in group]
    law = allo ? Set(s for group in EnzymeRates._inactive_groups(mech) for s in group) :
                 Set{EnzymeRates.Step}()
    tag(g) = allo ? EnzymeRates.cat_allo_state(mech, g) : :None
    index = Dict{Any, Int}()
    in_law = Bool[]
    function class(key, inlaw)
        i = get!(index, key) do
            push!(in_law, false)
            length(in_law)
        end
        in_law[i] |= inlaw
        i
    end
    function classes(j, inactive)
        s, g = flat[j]
        slots = EnzymeRates.is_equilibrium(s) ? (1,) : (1, 2)
        !inactive || tag(g) === :EqualAI ? [class((:A, g, k), true) for k in slots] :
        tag(g) === :OnlyA ? [class((:I, j, k), false) for k in slots] :
                            [class((:I, g, k), s in law) for k in slots]
    end
    foreach(j -> classes(j, false), eachindex(flat))
    allo && foreach(j -> classes(j, true), eachindex(flat))
    rows = Vector{Int}[]
    for z in _testhelper_cycle_basis(groups), inactive in (allo ? (false, true) : (false,))
        row = zeros(Int, length(in_law))
        for (j, (s, _)) in enumerate(flat)
            z[j] == 0 && continue
            cs = classes(j, inactive)
            if EnzymeRates.is_equilibrium(s)
                row[cs[1]] += EnzymeRates.is_binding(s) ? -z[j] : z[j]
            else
                row[cs[1]] += z[j]
                row[cs[2]] -= z[j]
            end
        end
        push!(rows, row)
    end
    A = isempty(rows) ? zeros(Int, 0, length(in_law)) : permutedims(reduce(hcat, rows))
    n_law = count(in_law) - (_testhelper_rank(A) - _testhelper_rank(A[:, .!in_law]))
    allo || return n_law
    n_reg = sum((t === :NonequalAI ? 2 : 1 for site in EnzymeRates.regulatory_sites(mech)
                 for t in EnzymeRates.allo_states(site)); init = 0)
    n_law + n_reg + 1
end

"""The mechanisms `init_mechanisms(rxn)` seeds and `rounds` rounds of `expand_mechanisms`
reach, without duplicates."""
function _testhelper_population(rxn, rounds)
    level = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[
        m for m in EnzymeRates.init_mechanisms(rxn)]
    mechs = copy(level)
    for _ in 1:rounds
        level = unique!(EnzymeRates.expand_mechanisms(level, rxn))
        append!(mechs, level)
    end
    unique!(mechs)
end

```

- [ ] **Step 2: Add the population sample and the detailed-balance test**

In `test/test_oracles.jl`, after the line `using Random` add:

```julia

# The ping-pong-capable bi-bi reaction with two catalytic subunits: its seeds hold every
# bi-bi mechanism, and one expansion round adds the allosteric children.
const _testhelper_pingpong_allo_rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    oligomeric_state: 2
end

# Mechanisms sampled from the seeds and one expansion round of the reaction above, the
# same sample on every run.
const _testhelper_population_sample = let
    population = _testhelper_population(_testhelper_pingpong_allo_rxn, 1)
    population[randperm(Random.MersenneTwister(2026), length(population))[1:120]]
end
```

At the end of the file add:

```julia

@testset "Detailed balance over a regenerated population" begin
    rng = Random.MersenneTwister(7)
    @test count(m -> m isa EnzymeRates.AllostericMechanism,
                _testhelper_population_sample) >= 30
    for m in _testhelper_population_sample
        @test all(<(big(10)^-30), _testhelper_detailed_balance(m, rng))
    end
end
```

- [ ] **Step 3: Pin every spec's fitted count to its declared dimension**

In `test/test_rate_eq_derivation.jl`, in `_testhelper_test_constraint_counting`, after
`        @test length(indep) == spec.expected_n_independent_params` add:

```julia
        @test length(indep) == _testhelper_declared_dimension(m)
```

- [ ] **Step 4: Run the tests (they pass against today's code)**

Run the focused commands for `test_oracles` and `test_rate_eq_derivation`.
Expected: all pass. `Detailed balance over a regenerated population | 121 121` (about
10 s); each spec's `Constraints` testset gains one passing test. Today's laws satisfy
detailed balance; only their parameter counts miss the declared dimension, and only on
populations (Task 7 adds that test).

- [ ] **Step 5: Run the full suite, then commit**

Run the full suite. Expected: all pass.

```bash
git status
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_oracles.jl \
        test/test_rate_eq_derivation.jl
git commit -m "Test detailed balance over a population and count each spec's dimension" \
  -m "A sample of 120 mechanisms regenerated from the ping-pong-capable bi-bi seeds and one
expansion round must sit at equilibrium to 1e-30, per conformation and combined. Each
spec's fitted count must equal the dimension the test counts itself: the classes in the
law less the rank of the cycle rows over them." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---
### Task 3: Export `all_parameters`

**Files:**
- Modify: `src/rate_eq_derivation.jl` (after the two concrete forwarders of
  `parameters`/`fitted_params`, lines 60–64)
- Modify: `src/EnzymeRates.jl` (export; the ABOUTME count)
- Modify: `test/test_accessors.jl` (new testsets)
- Modify: `test/test_rate_eq_derivation.jl` (replace `_testhelper_compute_all_params`)

**Interfaces:**
- Consumes: `_dependent_param_exprs` (today's; the joint solve's after Task 7 — the
  function's contract does not change), `_concrete`.
- Produces: exported `all_parameters(m::Union{AbstractEnzymeMechanism, Mechanism,
  AllostericMechanism}, params::NamedTuple) -> NamedTuple` — `params` merged with every
  dependent constant, in name order after `params`' own entries; the internal
  `_evaluate_dependent(expr, params)`.

Spec section 6 and section 18 item 1: `all_parameters` replaces the test helper
`_testhelper_compute_all_params` and is public API.

- [ ] **Step 1: Write the failing tests**

Append to `test/test_accessors.jl`:

```julia

@testset "all_parameters adds every dependent constant" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    p = (K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0, Keq = 4.0, E_total = 1.0)
    all_p = all_parameters(m, p)
    # The Haldane relation: k_EP_to_ES = K_EP_to_E_P · k_ES_to_EP / (K_ES_to_E_S · Keq).
    @test all_p.k_EP_to_ES == 1.0 * 3.0 / (2.0 * 4.0)
    @test keys(all_p) == (keys(p)..., :k_EP_to_ES)
    @test all_p[keys(p)] == p
    # The concrete mechanism gives the same values, E_total may be left out, and the
    # values may be of any number type.
    @test all_parameters(EnzymeRates.Mechanism(m), p) == all_p
    @test all_parameters(m, Base.structdiff(p, (E_total = 0,))).k_EP_to_ES ==
          all_p.k_EP_to_ES
    @test all_parameters(m, map(big, p)).k_EP_to_ES == big(3) / 8
    @test_throws ErrorException("all_parameters: params lacks K_ES_to_E_S, Keq, which " *
                                "the dependent constants need") all_parameters(
        m, (K_EP_to_E_P = 1.0, k_ES_to_EP = 3.0))
end

@testset "all_parameters reproduces the dependent constants of a consistent point" begin
    # The point's constants come from free energies, independently of the solve; the
    # fitted ones and Keq determine the rest. `half` has a dependent with
    # half-integer exponents.
    half = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E(S) <--> F(S), E <--> F)
            F(S) <--> E + P
            F + S <--> E(S)
        end
    end
    rng = Random.MersenneTwister(12)
    plain = [s.mechanism for s in MECHANISM_TEST_SPECS
             if s.mechanism isa EnzymeRates.EnzymeMechanism]
    for m in [plain; half]
        point = setprecision(() -> _testhelper_consistent_point(m, rng), BigFloat, 256)
        fitted = EnzymeRates.fitted_params(m)
        p = NamedTuple{(fitted..., :Keq)}((getindex.(Ref(point.vals), fitted)...,
                                           point.Keq))
        all_p = setprecision(() -> all_parameters(m, p), BigFloat, 256)
        @test all(isapprox(all_p[k], point.vals[k]; rtol = big(10)^-40)
                  for k in keys(all_p) if k !== :Keq)
    end
end
```

- [ ] **Step 2: Run the tests to verify they fail**

Run the focused command for `test_accessors`.
Expected: the two new testsets error with `UndefVarError: all_parameters not defined`.

- [ ] **Step 3: Implement `all_parameters`**

In `src/rate_eq_derivation.jl`, immediately after

```julia
fitted_params(m::Union{Mechanism, AllostericMechanism}) =
    fitted_params(compile_mechanism(m))
```

insert:

````julia

"""
    all_parameters(m, params) -> NamedTuple

`params` with every dependent constant of `m`'s rate law added, each computed from the
independent constants and `Keq` in `params` by its Haldane or Wegscheider relation (the
relations [`rate_equation_string`](@ref) prints). `params` holds the values of
[`parameters`](@ref)`(m)` by name; `E_total` may be left out, and every entry is kept as
passed. The values may be of any number type.

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S ⇌ E(S)
               E(S) <--> E(P)
               E(P) ⇌ E + P
           end
       end;

julia> all_parameters(m, (K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0, Keq = 4.0))
(K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0, Keq = 4.0, k_EP_to_ES = 0.375)
```
"""
function all_parameters(m::Union{AbstractEnzymeMechanism, Mechanism, AllostericMechanism},
                        params::NamedTuple)
    dep, indep = _dependent_param_exprs(m isa AbstractEnzymeMechanism ? _concrete(m) : m)
    absent = [k for k in (indep..., :Keq) if !haskey(params, k)]
    isempty(absent) || error("all_parameters: params lacks $(join(absent, ", ")), " *
                             "which the dependent constants need")
    names = sort!(collect(keys(dep)))
    merge(params, NamedTuple{Tuple(names)}(
        Tuple(_evaluate_dependent(dep[k], params) for k in names)))
end

"""The value at `params` of a dependent constant's expression (`build_power_expr`): a
power product of `Keq` and independent constants, or the number it is pinned to."""
_evaluate_dependent(x::Symbol, params) = params[x]
_evaluate_dependent(x::Real, params) = x
function _evaluate_dependent(e::Expr, params)
    op, args = e.args[1], [_evaluate_dependent(a, params) for a in e.args[2:end]]
    op === :* ? prod(args) : op === :/ ? args[1] / args[2] : args[1]^args[2]
end
````

In `src/EnzymeRates.jl`, change
`export rate_equation, rate_equation_string, parameters, metabolites` to
`export rate_equation, rate_equation_string, parameters, all_parameters, metabolites`, and
the ABOUTME line `# ABOUTME: Exports the 18 public names and wires together src/ source files.`
to `# ABOUTME: Exports the 19 public names and wires together src/ source files.`

- [ ] **Step 4: Run the tests to verify they pass**

Run the focused command for `test_accessors`. Expected: all pass.

- [ ] **Step 5: Replace the test helper**

In `test/test_rate_eq_derivation.jl`:

1. Delete `_testhelper_compute_all_params` with its docstring (lines 479–489) and
   `_testhelper_eval_dep_expr` with its docstring (lines 491–500).
2. In `_testhelper_random_independent_params_concs`, replace
   `    all_params = _testhelper_compute_all_params(m, new_params)` with
   `    all_params = all_parameters(m, new_params)`.
3. In `_testhelper_test_analytical_kcat`, replace
   `                      m, _testhelper_compute_all_params(m, params);` with
   `                      m, all_parameters(m, params);`.

Run the focused command for `test_rate_eq_derivation`. Expected: all pass.

- [ ] **Step 6: Run the full suite, then commit**

Run the full suite (it runs Aqua's export checks). Expected: all pass.

```bash
git status
git add src/rate_eq_derivation.jl src/EnzymeRates.jl test/test_accessors.jl \
        test/test_rate_eq_derivation.jl
git commit -m "Export all_parameters: the dependent constants from the independent ones" \
  -m "all_parameters(m, params) merges into params every dependent constant of the rate
law, computed by its Haldane or Wegscheider relation from the independent constants and
Keq. It replaces the test helper that did the same, and it reproduces the dependents of
a free-energy point on every plain spec and on a dependent with half-integer exponents." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---

### Task 4: Derive the rate polynomials from a reaction, groups and names

**Files:**
- Modify: `src/rate_eq_derivation.jl` (`_concentration_symbols`, `_compute_alpha`,
  `_raw_symbolic_rate_polys`, `_eq_complexity`, `_assert_derivable`, `_num_den_exprs` for
  `EnzymeMechanism`, `_state_rate_polys`; new `_leaf_names`)
- Modify: `test/test_mechanism_enumeration.jl` (lines 73–76, 106–111, 2751, 8553–8554,
  8675–8677)

**Interfaces:**
- Consumes: today's `_step_parameters`, `_build_wegscheider_rename_map`, `_state_parts`.
- Produces (Task 7 relies on these exact signatures):
  - `_raw_symbolic_rate_polys(rxn::EnzymeReaction, groups::Vector{Vector{Step}},
    names::Vector{Vector{Symbol}}) -> (num::POLY, den::POLY, d_free::POLY)`, with
    `names[j]` the Symbols of flat step `j`'s constants (`[K]` or `[kf, kr]`)
  - `_raw_symbolic_rate_polys(m::Mechanism)` and
    `_raw_symbolic_rate_polys(M::Type{<:EnzymeMechanism})`
  - `_compute_alpha(flat, conc_set, species, idx, segments, step_to_K)`
  - `_concentration_symbols(rxn::EnzymeReaction)`
  - `_eq_complexity(groups::Vector{Vector{Step}})`, `_eq_complexity(m::_AnyMech)`
  - `_assert_derivable(groups::Vector{Vector{Step}})`
  - `_leaf_names(mech::Mechanism, step_params, rename)` (Task 7 deletes it)

A pure refactor: the engine stops reading a `Mechanism`, so Task 7 can derive the
inactive conformation on `_inactive_groups(am)` with the full mechanism's own steps and
names. No output changes; every test passes unchanged apart from the internal calls
listed in Step 2.

- [ ] **Step 1: Refactor the engine**

In `src/rate_eq_derivation.jl`:

1. Replace `_concentration_symbols(mech::Mechanism) = Set(_metabolite_names(reaction(mech)))`
   with `_concentration_symbols(rxn::EnzymeReaction) = Set(_metabolite_names(rxn))`.
2. Replace the docstring and signature head of `_compute_alpha`:

```julia
"""
Compute alpha factors (relative concentrations within RE segments) as POLY
values. Iterates the flat steps `flat` directly; each RE step's weight ratio comes
from `_re_weight_ratio`, entering `to_species` from `from_species` or the
reverse. `idx[sp]` is the position of form `sp` in `species`; `step_to_K[j]` is
the parameter Symbol for the RE step at flat position `j`, and `conc_set` the
concentration symbols. Raises when the RE steps close a catalytic cycle, which gives
the mechanism no finite rate.
"""
function _compute_alpha(flat, conc_set, species, idx, segments, step_to_K)
    N = length(species)
    alpha = Vector{POLY}(fill(poly_one(), N))
```

   deleting its old line `    flat = _flat_steps(mech)`, and in its RE-turnover check
   delete the line `    conc_set = _concentration_symbols(mech)`.
3. Replace `_raw_symbolic_rate_polys` (its docstring and both methods, lines 209–301)
   with:

```julia
"""
    _raw_symbolic_rate_polys(rxn, groups, names) -> (num, den, d_free)

Build raw numerator and denominator POLYs for the rate equation of the kinetic groups
`groups` of a mechanism of reaction `rxn`. `names[j]` holds the parameter Symbols of
flat step `j`'s constants (`[K]` or `[kf, kr]`), which are the leaves of `num`/`den`.
Aborts, before the King–Altman expansion, a graph whose denominator would be too large
to derive (`_assert_derivable`).

Also returns `d_free`, the weight in the returned denominator of the free
resting enzyme (the form with empty `bound` and empty `residual`): the
spanning-tree weight `D[g_free]` of its segment times the concentration monomial
that brings `num`/`den` to lowest terms. Free E roots its segment, so its own
`alpha` is 1.
"""
function _raw_symbolic_rate_polys(rxn::EnzymeReaction, groups::Vector{Vector{Step}},
                                  names::Vector{Vector{Symbol}})
    _assert_derivable(groups)
    species, segments, _, idx, seg = _re_segment_extras(groups)
    # A fully-inert conformation (every binding pruned) has no enumerated form; it
    # exists only as free enzyme — no flux, partition 1, D[g_free] 1.
    isempty(species) && return poly_zero(), poly_one(), poly_one()
    flat = _flat_steps(groups)
    conc_set = _concentration_symbols(rxn)
    step_to_K = Dict{Int, Symbol}(
        i => names[i][1] for i in eachindex(flat) if is_equilibrium(flat[i][1]))
    alpha = _compute_alpha(flat, conc_set, species, idx, segments, step_to_K)
    G = length(segments)

    # The Laplacian of the segment graph: each SS step joining two segments adds its
    # forward flux out of its source segment and its reverse flux out of its target. A
    # step within one segment enters only the numerator. `ss` keeps every SS step with
    # its two forms and fluxes for the numerator.
    L = [poly_zero() for _ in 1:G, _ in 1:G]
    ss = Tuple{Step, Int, Int, POLY, POLY}[]
    for (j, (s, _)) in enumerate(flat)
        is_equilibrium(s) && continue
        a, b = idx[from_species(s)], idx[to_species(s)]
        fwd = _ss_contrib(poly_sym(names[j][1]), consumed(s), a, alpha)
        rev = _ss_contrib(poly_sym(names[j][2]), released(s), b, alpha)
        push!(ss, (s, a, b, fwd, rev))
        g1, g2 = seg[a], seg[b]
        g1 == g2 && continue
        L[g1, g1] = poly_add(L[g1, g1], fwd); L[g1, g2] = poly_sub(L[g1, g2], fwd)
        L[g2, g2] = poly_add(L[g2, g2], rev); L[g2, g1] = poly_sub(L[g2, g1], rev)
    end
    # D[g], the spanning-tree weight of segment g: L without g's row and column.
    D = [(o = setdiff(1:G, g); sym_det(L[o, o])) for g in 1:G]

    # A graph with no unbound, residual-free form (e.g. a hand-written mechanism whose
    # every form carries a residual) has no free enzyme to weight, so `d_free` is 1.
    i_free = findfirst(f -> isempty(bound(f)) && isempty(residual(f)), species)
    d_free = i_free === nothing ? poly_one() : D[seg[i_free]]

    den = poly_zero()
    for g in 1:G
        sigma = reduce(poly_add, (alpha[i] for i in segments[g]); init=poly_zero())
        den = poly_add(den, poly_mul(sigma, D[g]))
    end

    # Numerator of the rate: v·den summed over steady-state steps. With u(f) the
    # first substrate's exponent in form f's RE weight, each SS step e contributes
    # ω_e·(forward − reverse flux) with ω_e = (copies of the first substrate e
    # consumes − copies it releases) + u(from_e) − u(to_e). Flux conservation at every
    # form makes the sum equal the net consumption of the first substrate for every
    # parameter value, so it needs no choice of reaction cut and does not depend on
    # the direction a step is written in. Steps with ω_e = 0 contribute nothing.
    x = name(first(substrates(rxn)))
    expo(p) = (mono = only(keys(p));
               k = findfirst(q -> q.first == x, mono);
               k === nothing ? 0 : mono[k].second)
    num = poly_zero()
    for (s, a, b, fwd, rev) in ss
        ω = count(m -> name(m) == x, consumed(s)) - count(m -> name(m) == x, released(s)) +
            expo(alpha[a]) - expo(alpha[b])
        ω == 0 && continue
        term = poly_sub(poly_mul(fwd, D[seg[a]]), poly_mul(rev, D[seg[b]]))
        num = poly_add(num, POLY(key => ω * c for (key, c) in term))
    end

    _reduce_conc_lowest_terms(num, den, d_free, conc_set)
end

_raw_symbolic_rate_polys(m::Mechanism) = _raw_symbolic_rate_polys(
    reaction(m), steps(m),
    _leaf_names(m, _step_parameters(m), _build_wegscheider_rename_map(m)))
_raw_symbolic_rate_polys(@nospecialize(M::Type{<:EnzymeMechanism})) =
    _raw_symbolic_rate_polys(Mechanism(M()))

"""
The parameter Symbols of each flat step's constants, `step_params` rendered through the
naming chokepoint on `mech`, a rapid-equilibrium step's constant through `rename` (the
single-symbol Wegscheider ties, `_build_wegscheider_rename_map`).
"""
_leaf_names(mech::Mechanism, step_params, rename) =
    [Symbol[(n = name(p, mech); p isa Krapid ? get(rename, n, n) : n) for p in ps]
     for ps in step_params]
```

4. In `_eq_complexity`, change the signature to
   `function _eq_complexity(groups::Vector{Vector{Step}})`, its first line to
   `    _, segments, _, idx, seg = _re_segment_extras(groups)` and its loop header to
   `    for group in groups, s in group`; replace the line
   `_eq_complexity(m::AllostericMechanism) = _eq_complexity(_state_mechanism(m, :A))`
   with `_eq_complexity(m::_AnyMech) = _eq_complexity(steps(m))`.
5. Change `_assert_derivable(mech::Mechanism)` to
   `_assert_derivable(groups::Vector{Vector{Step}})` and its body's call to
   `_eq_complexity(groups)`.
6. Replace `_num_den_exprs(@nospecialize(M::Type{<:EnzymeMechanism}); horner = false)`
   (the plain method, lines 365–370) with:

```julia
function _num_den_exprs(@nospecialize(M::Type{<:EnzymeMechanism}); horner = false)
    m = Mechanism(M())
    names = _leaf_names(m, _step_parameters(m), _build_wegscheider_rename_map(m))
    num, den, _ = _raw_symbolic_rate_polys(reaction(m), steps(m), names)
    horner && return _poly_to_horner(num), _poly_to_horner(den)
    param_syms = Set{Symbol}(n for ns in names for n in ns)
    _poly_to_expr(num, param_syms), _poly_to_expr(den, param_syms)
end
```

7. Replace `_state_rate_polys` (the one-line definition, lines 736–737; keep its
   docstring) with:

```julia
function _state_rate_polys(am::AllostericMechanism, state::Symbol)
    cm, sp = _state_parts(am, state)
    rename = _build_wegscheider_rename_map(cm; step_params = sp)
    _raw_symbolic_rate_polys(reaction(cm), steps(cm), _leaf_names(cm, sp, rename))
end
```

- [ ] **Step 2: Update the internal calls in the enumeration tests**

In `test/test_mechanism_enumeration.jl`:

1. In `_testhelper_identifiable_rank`, replace
   ```julia
       cm = m isa EnzymeRates.Mechanism ? m : EnzymeRates._state_mechanism(m, :A)
       mets = sort!(collect(EnzymeRates._concentration_symbols(cm)))
   ```
   with
   ```julia
       mets = sort!(collect(EnzymeRates._concentration_symbols(EnzymeRates.reaction(m))))
   ```
2. In `_testhelper_degenerate`, replace
   ```julia
       cm = m isa EnzymeRates.Mechanism ? m : EnzymeRates._state_mechanism(m, :A)
       subs = [EnzymeRates.name(s) for s in EnzymeRates.substrates(EnzymeRates.reaction(cm))]
       prods = [EnzymeRates.name(p) for p in EnzymeRates.products(EnzymeRates.reaction(cm))]
       mets = sort!(collect(EnzymeRates._concentration_symbols(cm)))
   ```
   with
   ```julia
       rxn = EnzymeRates.reaction(m)
       subs = [EnzymeRates.name(s) for s in EnzymeRates.substrates(rxn)]
       prods = [EnzymeRates.name(p) for p in EnzymeRates.products(rxn)]
       mets = sort!(collect(EnzymeRates._concentration_symbols(rxn)))
   ```
3. In the testset "AllostericMechanism — a split whose count needs a law over the term
   cap", replace `    @test EnzymeRates._eq_complexity(EnzymeRates._state_mechanism(am, :A)) >`
   with `    @test EnzymeRates._eq_complexity(am) >`.
4. Replace the two three-argument calls
   ```julia
           EnzymeRates._raw_symbolic_rate_polys(k, EnzymeRates._step_parameters(k),
               EnzymeRates._build_wegscheider_rename_map(k))[2]))
   ```
   with `        EnzymeRates._raw_symbolic_rate_polys(k)[2]))`, and
   ```julia
                   _, den, _ = EnzymeRates._raw_symbolic_rate_polys(
                       m, EnzymeRates._step_parameters(m),
                       EnzymeRates._build_wegscheider_rename_map(m))
   ```
   with `                _, den, _ = EnzymeRates._raw_symbolic_rate_polys(m)`.

- [ ] **Step 3: Run the tests**

Run the focused commands for `test_rate_eq_derivation`, `test_allosteric_golden` and
`test_mechanism_enumeration`. Expected: all pass, the golden strings byte-identical.

- [ ] **Step 4: Run the full suite, then commit**

Run the full suite. Expected: all pass.

```bash
git status
git add src/rate_eq_derivation.jl test/test_mechanism_enumeration.jl
git commit -m "Derive rate polynomials from a reaction, kinetic groups and their names" \
  -m "The King–Altman engine reads the groups and each flat step's parameter names instead
of a Mechanism, so a conformation's law can be derived on a subset of a mechanism's own
steps without rebuilding a mechanism. Output is unchanged." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---

### Task 5: Remove the `Full` and `Reduced` modes

**Files:**
- Modify: `src/types.jl` (delete the mode types and constants, lines 1250–1283)
- Modify: `src/rate_eq_derivation.jl` (`parameters`, `rate_equation`,
  `rate_equation_string` and their docstrings; `_raw_param_symbols` users)
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl` (delete
  `_raw_param_symbols`)
- Modify: `src/EnzymeRates.jl` (exports)
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
  (`_testhelper_step_constant_names`), `test/test_accessors.jl`,
  `test/test_allosteric_golden.jl`, `test/test_rate_eq_derivation.jl`,
  `test/test_types.jl`
- Modify: `docs/src/deriving/textbooks.md`, `docs/src/developer.md`, `ROADMAP.md`

**Interfaces:**
- Consumes: Task 3's `all_parameters` (the docs page shows it where the Full section
  stood).
- Produces: `parameters(m)`, `rate_equation(m, concs, params)`,
  `rate_equation_string(m)` with no mode argument;
  `_testhelper_step_constant_names(m) -> Set{Symbol}` in the fixtures file.

Spec B5: no all-parameters law. `all_parameters` lists every constant the law reads; the
naming tests that read `parameters(m, Full)` read each step constant's name through the
chokepoint instead, which names the same constants.

- [ ] **Step 1: Rewrite the tests of the removed modes**

1. In the fixtures file, after the concrete `_testhelper_lift` method add:

```julia

"""The name of every step constant of the plain mechanism `m`, fitted or dependent: each
flat step's constants as the naming chokepoint renders them."""
_testhelper_step_constant_names(m) =
    Set(EnzymeRates.name(p, mech) for mech in (_testhelper_lift(m),)
        for (s, _) in EnzymeRates._flat_steps(mech)
        for p in EnzymeRates._step_constants(s, :None))
```

2. `test/test_accessors.jl`, testset "parameters API symmetry": replace lines 10–16
   (`full_mono = …` through `@test length(full_mono) >= length(reduced_mono)`) with

```julia
    reduced_mono = parameters(m_mono)
    @test :E_total in reduced_mono
    @test :Keq in reduced_mono
    values_mono = NamedTuple{reduced_mono}(ntuple(_ -> 1.0, length(reduced_mono)))
    @test issubset(reduced_mono, keys(all_parameters(m_mono, values_mono)))
```

   and replace lines 33–35 (the comment, `@test_throws MethodError parameters(m_allo,
   Full)` and `reduced_allo = parameters(m_allo, Reduced)`) with

```julia
    reduced_allo = parameters(m_allo)
    values_allo = NamedTuple{reduced_allo}(ntuple(_ -> 1.0, length(reduced_allo)))
    @test issubset(reduced_allo, keys(all_parameters(m_allo, values_allo)))
```

3. `test/test_allosteric_golden.jl`: replace
   `EnzymeRates.rate_equation_string(m, EnzymeRates.Reduced)` with
   `EnzymeRates.rate_equation_string(m)`, and `@test parameters(m, EnzymeRates.Reduced) ==`
   with `@test parameters(m) ==`; in the comment above that test replace
   ``# `parameters(m, Reduced)` wrapper`` with ``# `parameters(m)` wrapper``.
4. `test/test_types.jl`:
   - testset "fused steps are named by their sides": replace
     `    names = Set(ER.parameters(m, ER.Full))` with
     `    names = _testhelper_step_constant_names(m)`, and
     `    kon, koff = ER._step_parameters(mech)[i]` with
     `    kon, koff = ER._step_constants(first(ER._flat_steps(mech)[i]), :None)`.
   - testset "step constants are named by their reaction": replace
     `    full_names(m) = Set(ER.parameters(m, ER.Full))` with
     `    full_names(m) = _testhelper_step_constant_names(m)`; replace the comment
     `    # Reduced mode: the Haldane reduction leaves the SS isomerization constant fitted.`
     with `    # The Haldane reduction leaves the SS isomerization constant fitted.`; replace
     `    @test length(pp_names) == length(ER.parameters(pingpong, ER.Full))` with

```julia
    # One name per class: a rapid-equilibrium group's K, a steady-state group's kf, kr.
    @test length(pp_names) == sum(ER.is_equilibrium(first(g)) ? 1 : 2
                                  for g in ER.steps(ER.Mechanism(pingpong)))
```

   - testset "Tier 2 reads the free metabolites at both ends of a step": replace
     `    kf, kr = ER._step_parameters(m)[i]` with
     `    kf, kr = ER._step_constants(first(ER._flat_steps(m)[i]), :None)`.
   - the testset holding `names = [ER.name(p, m) for p in ER._enumerate_parameters_full(m)]`
     (line 2094): replace that line with

```julia
    names = [ER.name(p, m) for g in ER.steps(m) for p in ER._step_constants(first(g), :None)]
```

     (one name per class, so `allunique(names)` still checks that no two classes
     render one name).
5. `test/test_rate_eq_derivation.jl`:
   - in `_testhelper_test_structure` replace the comment lines
     `        # paths the fitter consumes (Reduced + fitted_params). A collision` with
     `        # paths the fitter consumes (parameters + fitted_params). A collision`, and
     `        @test allunique(EnzymeRates.parameters(m))           # Reduced (default)`
     with `        @test allunique(EnzymeRates.parameters(m))`.
   - testset "Theorell–Chance and two-metabolite steps": replace
     `    @test :k_EA_B_to_EQ_P in ER.parameters(tc_ss, ER.Full)` and the next line with
     ```julia
         @test :k_EA_B_to_EQ_P in _testhelper_step_constant_names(tc_ss)
         @test :k_EQ_P_to_EA_B in _testhelper_step_constant_names(tc_ss)
     ```
   - testset "rate_equation_string prints flat +/* sums (no nested parens)": replace
     its comment and last two lines with
     ```julia
         # Orthogonal guard (NOT a perf-fix backstop): _expr_to_string is
         # precedence-aware and flattens nested +/* nodes
         # transparently: balanced +(+(a,b), +(c,d)) prints as "a + b + c + d"
         # with no added parens. If that flattening regresses, the printed
         # output would gain nested parens but still evaluate correctly —
         # easy to miss without an explicit guard. The `v =` line of Random-order
         # Bi-Bi holds three: the numerator, the denominator and the negative terms
         # inside the numerator.
         spec = only(s for s in MECHANISM_TEST_SPECS
                     if s.name == "Random-order Bi-Bi")
         s = last(split(rate_equation_string(spec.mechanism), '\n'))
         @test count(==('('), s) <= 3
     ```
   - in the clash testset, replace the comment
     `    # An inhibitor named after the reverse chemistry constant: in Reduced mode the`
     `    # Haldane assignment would overwrite its concentration, in Full mode its`
     `    # concentration would overwrite the constant.`
     with
     `    # An inhibitor named after the reverse chemistry constant: the Haldane assignment`
     `    # would overwrite its concentration.`
     and delete the three lines
     `    full = NamedTuple{parameters(uni_k, Full)}((1.0, 2.0, 3.0, 4.0, 5.0, 1.0))`,
     `    @test_throws clash_k rate_equation(uni_k, concs, full, Full)` and
     `    @test_throws clash_k rate_equation_string(uni_k, Full)` (the default-mode
     assertions on the same fixture stay and test the clash that remains).
   - testset "Fix A: dead-inactive-state allosteric body defines all I-state symbols":
     replace `rate_equation(m, concs, params, Reduced)` with
     `rate_equation(m, concs, params)`.

- [ ] **Step 2: Run the rewritten tests**

Run the focused commands for `test_accessors` and `test_types`.
Expected: all pass. The rewritten tests read only what stays when the modes go, so they
pass before Step 3 and after it; running them now shows they are right before the code
they replace disappears.

- [ ] **Step 3: Remove the modes**

1. `src/types.jl`: delete the block from `# --- Rate equation mode types ---` through
   `const Reduced = ReducedMode()` and its blank line (lines 1250–1283).
2. `src/thermodynamic_constr_for_rate_eq_derivation.jl`: delete `_raw_param_symbols`
   (both methods and the docstring, lines 16–26).
3. `src/EnzymeRates.jl`: delete the lines `# Rate equation modes` and
   `export Full, Reduced` and the blank line after them; change the ABOUTME count from
   19 to 17.
4. `src/rate_eq_derivation.jl`:
   - Replace the `parameters` docstring and its first four methods (from the docstring
     `    parameters(m::EnzymeMechanism, [mode])` through
     `@generated function parameters(::M, ::ReducedMode) … end`) with

```julia
"""
    parameters(m)

The names of the parameters the rate law of `m` takes, as a tuple of Symbols: its
independent constants (`fitted_params`), then `Keq` and `E_total`. The rate law computes
every other constant from these (see [`all_parameters`](@ref)).
"""
@generated function parameters(::M) where {M <: AbstractEnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    (indep..., :Keq, :E_total)
end
```

   - Replace the concrete forwarder
     `parameters(m::Union{Mechanism, AllostericMechanism}, mode::AbstractRateEquationMode) =`
     `    parameters(compile_mechanism(m), mode)` with
     `parameters(m::Union{Mechanism, AllostericMechanism}) = parameters(compile_mechanism(m))`.
   - Replace the block from the line `# ─── Mode-dispatched rate_equation ───…` up to the
     line `# ─── String Representation ───…` (the docstring, `function rate_equation
     end`, the default-mode method, the concrete forwarder and the `FullMode` and
     `ReducedMode` methods) with:

```julia
# ─── rate_equation ─────────────────────────────────────────────

"""
    rate_equation(m, concs, params)

Return the net reaction rate of mechanism `m` at the metabolite concentrations
`concs` and parameters `params` (both `NamedTuple`s). The rate equation itself
is derived from the mechanism's elementary steps by the King–Altman/Cha method:
rapid-equilibrium steps collapse to binding constants, steady-state steps are
assembled into the King–Altman determinant, and the two combine into the
quasi-steady-state flux through the whole mechanism — so one call returns the
overall turnover, not the rate of a single step.

`params` holds the independent constants (`parameters(m)`: the fitted constants, `Keq`
and `E_total`); the rate law computes each dependent constant from them by its
Haldane or Wegscheider relation.

The derivation runs once, at compile time: the body is generated as a single
arithmetic expression with no allocations, loops, or matrix operations. Use
[`rate_equation_string`](@ref) to inspect the derived equation symbolically.
"""
@generated function rate_equation(m::M, concs::NamedTuple,
                                  params::NamedTuple) where {M <: AbstractEnzymeMechanism}
    dep, indep = _dependent_param_exprs(M)
    _rate_body(M, (indep..., :Keq, :E_total), dep)
end

# Concrete-mechanism convenience: lift to the singleton (see the note on
# the parameters/fitted_params methods above — allocates, not the hot path).
rate_equation(m::Union{Mechanism, AllostericMechanism}, concs, params) =
    rate_equation(compile_mechanism(m), concs, params)
```

   - Replace the block after `# ─── String Representation ───…` up to the docstring of
     `_equation_text` (the `rate_equation_string` docstring, `function
     rate_equation_string end`, the default-mode method and the concrete forwarder) with:

````julia
"""
    rate_equation_string(m) -> String

Return the symbolic rate equation for mechanism `m` as a multi-line
`String` (it returns the text — it does not print). Pass a concrete `Mechanism` /
`AllostericMechanism` or its compiled [`EnzymeMechanism`](@ref) singleton.

The string spells out how [`rate_equation`](@ref)
evaluates: a `(; …) = params` destructure line, a `(; …) = concs`
destructure line, then the `v = E_total * (num) / (den)` line. Dependent rate
constants are listed first under `# Wegscheider constraints:` and
`# Haldane constraints:` headers — the thermodynamic identities that eliminate
parameters — and only the independent set appears in the `params` destructure.

Names that are not Julia identifiers, such as the residual forms of
ping-pong mechanisms (`k_EA_to_EP_res_+A_-P`), print bare, so such a
string reads as text, not code.

Use `print` on the result to see the multi-line layout without escaped
newlines.

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S ⇌ E(S)
               E(S) <--> E(P)
               E(P) ⇌ E + P
           end
       end;

julia> print(rate_equation_string(m))
(; K_EP_to_E_P, K_ES_to_E_S, k_ES_to_EP, Keq, E_total) = params
(; S, P) = concs
# Haldane constraints:
k_EP_to_ES = (1 / Keq) * K_EP_to_E_P * (1 / K_ES_to_E_S) * k_ES_to_EP
v = E_total * (k_ES_to_EP * S / K_ES_to_E_S - k_EP_to_ES * P / K_EP_to_E_P) / (1 + P / K_EP_to_E_P + S / K_ES_to_E_S)
```
"""
function rate_equation_string end

# Concrete-mechanism convenience: lift to the singleton.
rate_equation_string(m::Union{Mechanism, AllostericMechanism}) =
    rate_equation_string(compile_mechanism(m))
````

   - Replace the three methods after `_equation_text` (the `FullMode` method, the
     commented `ReducedMode` method for `EnzymeMechanism` and the `ReducedMode` method
     for `AllostericEnzymeMechanism`) with:

```julia
# The display solves with no Wegscheider rename, so an absorbed single-symbol tie stays
# visible under `# Wegscheider constraints:`, marked as substituted into v.
function rate_equation_string(@nospecialize(m::EnzymeMechanism))
    _, indep = _dependent_param_exprs(typeof(m))
    mech = Mechanism(m)
    dep, _ = _solve_dependent_set(_assemble_constraints(mech, Dict{Symbol, Symbol}())...)
    _equation_text(typeof(m), (indep..., :Keq, :E_total), dep;
                   substituted = keys(_build_wegscheider_rename_map(mech)))
end

function rate_equation_string(@nospecialize(m::AllostericEnzymeMechanism))
    dep, indep = _dependent_param_exprs(typeof(m))
    _equation_text(typeof(m), (indep..., :Keq, :E_total), dep)
end
```

   (Task 7 merges the two `rate_equation_string` methods once the rename pass goes.) In
   the docstring of `ANNOTATION_SUBSTITUTED`, replace
   ``Suffix the `Reduced` `rate_equation_string` of an `EnzymeMechanism` appends`` with
   ``Suffix the `rate_equation_string` of an `EnzymeMechanism` appends``.

- [ ] **Step 4: Update the docs**

1. `docs/src/deriving/textbooks.md`: replace the two lines

```text
The default mode is `Reduced`.
In `Reduced` mode the string has four sections:
```

   with the single line `The string has four sections:`. Replace the whole
   `## Full mode` section (from `## Full mode` through the paragraph that ends
   `SS parameters, and how adding SS steps changes this list.`) with:

````markdown
## Every constant of the rate law

[`all_parameters`](@ref) adds the dependent constants to a set of independent values,
each computed from its Haldane or Wegscheider relation:

```@example textbook_ex
all_parameters(m, (K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0, Keq = 4.0))
```

Here the Haldane relation fixes the reverse catalytic rate `k_EP_to_ES` from `Keq`, the
two binding constants and the forward rate. See
[Rapid equilibrium vs steady state](@ref) for the contrast between RE and SS
parameters, and how adding SS steps changes the parameter list.
````

2. `docs/src/developer.md`: in the paragraph that starts
   ``Only the `@generated` methods``, replace the lines

```text
`fitted_params` and the reduced `parameters` run its constraint solve, and `metabolites` and
the full `parameters` read names off the lifted mechanism. The lifts `Mechanism(em)` and
```

   with

```text
`fitted_params` and `parameters` run its constraint solve, and `metabolites` reads names off
the lifted mechanism. The lifts `Mechanism(em)` and
```

   and the lines

```text
run, plus a few thin methods that still specialize on it: the forwarders that supply the
default mode (such as `rate_equation(m, concs, params)` and `rate_equation_string(m)`),
`show`, `catalytic_mechanism` and `catalytic_multiplicity`. `FittingProblem` and `loss!`
```

   with

```text
run, plus a few thin methods that still specialize on it: `show`, `catalytic_mechanism`
and `catalytic_multiplicity`. `FittingProblem` and `loss!`
```

   then re-wrap the paragraph to 92 characters.
3. `ROADMAP.md`: delete the known-issue entry **Full mode lists rapid-equilibrium
   constants that never reach the rate law.** (lines 83–92; spec section 16, phase-2
   row).

- [ ] **Step 5: Run the tests**

Run the focused commands for `test_accessors`, `test_types`, `test_rate_eq_derivation`
and `test_allosteric_golden`. Expected: all pass. Then confirm no mode remains:
`grep -rn "FullMode\|ReducedMode\|AbstractRateEquationMode\|\bFull\b\|\bReduced\b\|_raw_param_symbols" src test docs/src README.md`
prints nothing but prose (`Full documentation`, `Full symbolic`, `Full-rank`,
`Full multi-tier`, `Reduced row echelon form`).

- [ ] **Step 6: Run the full suite, then commit**

Run the full suite. Expected: all pass (test_compile_budget's method-instance count per
fresh mechanism may drop by one: the mode forwarders no longer specialize).

```bash
git status
git add src/types.jl src/rate_eq_derivation.jl \
        src/thermodynamic_constr_for_rate_eq_derivation.jl src/EnzymeRates.jl \
        test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_accessors.jl \
        test/test_allosteric_golden.jl test/test_rate_eq_derivation.jl test/test_types.jl \
        docs/src/deriving/textbooks.md docs/src/developer.md ROADMAP.md
git commit -m "Remove the Full and Reduced modes" \
  -m "parameters, rate_equation and rate_equation_string take no mode: the rate law is the
thermodynamically reduced one, and all_parameters adds the dependent constants. The
naming tests read each step constant's name through the chokepoint, and the textbooks
page shows all_parameters where the Full section stood." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---
### Task 6: The joint solve as a unit

**Files:**
- Create: `src/thermo.jl`
- Modify: `src/EnzymeRates.jl` (include `thermo.jl` after
  `thermodynamic_constr_for_rate_eq_derivation.jl`)
- Modify: `src/types.jl` (`_entry_sides` out of `_canonicalize_step_directions`; a
  `lead_key` keyword on `_orient_tied_steps`)
- Create: `test/test_joint_solve.jl`
- Modify: `test/runtests.jl` (include `test_joint_solve.jl` after `test_oracles.jl`)

**Interfaces:**
- Consumes: `_canonical_step_direction`, `_reversed`, `_flat_steps`, `_forms`,
  `_inactive_groups`, `_is_chemistry`, `is_iso`, `bound_metabolite`; Task 1's
  `_testhelper_steady_state_rate`, Task 2's `_testhelper_declared_dimension`.
- Produces (Task 7 relies on these exact names):
  - `struct _Column; step::Int; slot::Int; inactive::Bool; end`
  - `_metabolites_key(mets)`, `_form_base_key(sp::Species)`, `_is_free_E(sp)`,
    `_form_index_key(sp)`, `_step_index_key(s::Step)`
  - `_solve_forward(groups::Vector{Vector{Step}}) -> Vector{Bool}` (per flat step)
  - `_priority_key(s::Step, fwd::Bool, slot::Int, block::Int)` (a tuple; larger is
    eliminated first)
  - `_solve_tags(m)`, `_inactive_law(m) -> BitVector`, `_has_inactive_columns(tag)`,
    `_class(c::_Column, g::Int, tag::Symbol)`
  - `_fundamental_cycles(rxn, groups) -> (Z::Matrix{Int}, multiples::Vector{Int})`
  - `_joint_system(m::_AnyMech) -> (; columns, A::Matrix{Int}, rhs::Vector{Int}, law,
    tags, flat)`, columns in elimination order
  - `_rref!(M::Matrix{<:Rational}, n::Int) -> pivots`
  - `_reduce_joint_system(A::Matrix{Int}, rhs::Vector{Int}) -> (M, pivots)`
  - `_kept_members(m::_AnyMech) -> Dict{Tuple{Int,Int,Bool}, Step}` keyed by
    `(group, slot, inactive)`
  - `_entry_sides(groups)` and `_orient_tied_steps(oriented, tiers; lead_key =
    _step_canonical_key)` in `src/types.jl`

The solve is built and tested on its own here; nothing calls it until Task 7. It follows
spec section 6 exactly:

- **Columns.** One per slot of every flat step (slot 1: a rapid-equilibrium step's
  equilibrium constant or a steady-state step's rate constant in its stored direction;
  slot 2: the reverse rate), and for an allosteric mechanism one more per slot of every
  step of a `:NonequalAI` or `:OnlyA` group, for its inactive copy; an `:EqualAI` step's
  inactive constant reads its active column.
- **Rows.** The fundamental cycles of a spanning forest, once per conformation (the
  inactive rows over the whole graph, islands included), with right side the cycle's
  multiple of the net reaction; then one row per further member of each kinetic class,
  equal to the class's first member, slot by slot. A class is a group's constant in one
  slot and conformation, except that each `:OnlyA` step's inactive copy is a class of
  its own (B2). Integer arithmetic throughout; the guard keeps today's
  `_thermodynamic_constraints` (spec: today's guards stay unchanged).
- **Order.** `_priority_key(s, fwd, slot, block)`, larger first:
  `(block, kind, occupied, ss, against, content)` with block 2 for an inactive copy
  outside the law (an `:OnlyA` step, a step on a pruned form), 1 for another inactive
  copy, 0 for an active column; kind 3 for an isomerization, 2 for other chemistry, 1
  for a plain binding; occupied when both forms hold a metabolite; ss for steady state;
  against for the reverse of the solve's forward direction, except a product binding's
  uptake (E7); content the step's base-name key in the solve's direction (C8: conformation
  label, bound `(name, role)` pairs and residual name sets of each form, then the
  consumed and released pairs; never rendered names), the smaller key kept.
- **Orientation (C9).** `_solve_forward` re-runs today's orientation tiers; a step that
  tier 3 decides (by rendered name in the constructor) runs instead from the form earlier
  in the canonical form order (`_form_index_key`: free E first, then the base-name key),
  and tied group members align to a lead ranked by that order. Only the key reads it:
  columns stay attached to the stored `Kfor`/`Krev`, so names and the engine are
  unaffected. Measured: 0 of the steps of 20,867 enumerated mechanisms differ from their
  stored direction; the hand-written `E <--> D` of the test below does.
- **Reduction.** Exact Gauss–Jordan over `Rational{Int}`, which throws `OverflowError`
  instead of wrapping; on overflow the system is reduced again in `Rational{BigInt}`.
  Rows that reduce to `0 = log Keq` raise today's error text.
- **Kept members (E3).** The member of a class with the smallest key comes last in
  elimination order, so the solve keeps it whenever it keeps one of the class; every
  non-last member is always a pivot (its group row with the last member's row isolates
  it), so `_kept_members` computes the names without the reduction.

- [ ] **Step 1: Write the failing tests**

Create `test/test_joint_solve.jl`:

```julia
# ABOUTME: Tests of the joint thermodynamic solve: the priority key and orientation, the
# ABOUTME: cycles, the exact reduction, kept members, and its constants against the oracle.
using Random

const _testhelper_JS = EnzymeRates

"""
A parameter point of `m` built from its joint solve: every column that is not a pivot
takes a random log value, `log Keq` one more, and each pivot column the value its row of
the reduced system gives. Returns the steps of each conformation as
`_testhelper_steady_state_rate` reads them, with `L` and the regulator constants drawn at
random.
"""
function _testhelper_solve_point(m, rng)
    mech = _testhelper_lift(m)
    sys = _testhelper_JS._joint_system(mech)
    M, pivots = _testhelper_JS._reduce_joint_system(sys.A, sys.rhs)
    n = length(sys.columns)
    x = zeros(BigFloat, n)
    free = setdiff(1:n, pivots)
    foreach(f -> x[f] = 4 * big(rand(rng)) - 2, free)
    log_keq = 4 * big(rand(rng)) - 2
    for (r, p) in enumerate(pivots)
        x[p] = M[r, n + 1] * log_keq -
               sum((M[r, f] * x[f] for f in free); init = big(0))
    end
    value = Dict(c => exp(x[i]) for (i, c) in enumerate(sys.columns))
    allo = mech isa _testhelper_JS.AllostericMechanism
    function conformation(inactive)
        out = []
        for (j, (s, g)) in enumerate(sys.flat)
            inactive && !sys.law[j] && continue
            own = inactive && _testhelper_JS._has_inactive_columns(sys.tags[g])
            col(slot) = value[_testhelper_JS._Column(j, slot, own)]
            if _testhelper_JS.is_equilibrium(s)
                K = col(1)
                push!(out, (s, _testhelper_JS.is_binding(s) ? 1 / K : K, nothing))
            else
                push!(out, (s, col(1) / col(2), col(1)))
            end
        end
        out
    end
    vals = Dict{Symbol, BigFloat}()
    L = big(0)
    if allo
        L = exp(4 * big(rand(rng)) - 2)
        for p in [_testhelper_JS._kreg_params(mech, :A); _testhelper_JS._kreg_params(mech, :I)]
            vals[_testhelper_JS.name(p, mech)] = exp(4 * big(rand(rng)) - 2)
        end
    end
    (; vals, Keq = exp(log_keq), active = conformation(false),
     inactive = allo ? conformation(true) : [], L)
end

# Reversible Michaelis–Menten with every step at steady state.
const _testhelper_mm_ss = @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S <--> E(S)
        E(S) <--> E(P)
        E(P) <--> E + P
    end
end

@testset "Joint thermodynamic solve" begin
    @testset "the priority key orders the columns" begin
        mech = EnzymeRates.Mechanism(_testhelper_mm_ss)
        sys = EnzymeRates._joint_system(mech)
        describe(c) = (EnzymeRates._step_text(sys.flat[c.step][1]), c.slot)
        # The isomerization goes first, its reverse rate (against the reaction) before
        # its forward one; then the two bindings' constants that run against the
        # reaction, the release of S and the uptake of P (a product binding's uptake runs
        # against it), the larger content key first (S sorts after P); then the others.
        @test describe.(sys.columns) == [("ES <--> EP", 2), ("ES <--> EP", 1),
                                         ("E + S <--> ES", 2), ("E + P <--> EP", 1),
                                         ("E + S <--> ES", 1), ("E + P <--> EP", 2)]
        # One Haldane row: the first column, the isomerization's reverse rate, is the
        # dependent one.
        M, pivots = EnzymeRates._reduce_joint_system(sys.A, sys.rhs)
        @test pivots == [1]
        @test M[1, end] == -1
    end

    @testset "kind, occupancy, flag and direction rank in that order" begin
        E = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        A, B = EnzymeRates.Substrate(:A), EnzymeRates.Substrate(:B)
        EA = EnzymeRates.Species(EnzymeRates.Metabolite[A], :E)
        EAB = EnzymeRates.Species(EnzymeRates.Metabolite[A, B], :E)
        plain = EnzymeRates.Step(E, EA, [A], EnzymeRates.Metabolite[], false)
        occupied = EnzymeRates.Step(EA, EAB, [B], EnzymeRates.Metabolite[], false)
        occupied_re = EnzymeRates.Step(EA, EAB, [B], EnzymeRates.Metabolite[], true)
        key = EnzymeRates._priority_key
        # block first: an inactive copy outside the law before any other column
        @test key(plain, true, 1, 2) > key(occupied, true, 2, 1) > key(occupied, true, 2, 0)
        # a binding between two occupied forms before one that touches free E
        @test key(occupied, true, 1, 0) > key(plain, true, 2, 0)
        # steady state before rapid equilibrium, then the constant against the reaction
        @test key(occupied, true, 1, 0) > key(occupied_re, true, 1, 0)
        @test key(occupied, true, 2, 0) > key(occupied, true, 1, 0)
    end

    @testset "steps no chemistry orients run from the earlier form (C9)" begin
        # E ⇌ D is decided by neither atom progression nor graph context. The
        # constructor stores it by rendered name, D before E; the solve orients it from
        # free E, the first form of the canonical form order.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
                E <--> D
            end
        end)
        fwd = EnzymeRates._solve_forward(EnzymeRates.steps(m))
        flat = EnzymeRates._flat_steps(m)
        j = only(j for (j, (s, _)) in enumerate(flat)
                 if EnzymeRates.name(EnzymeRates.from_species(s)) === :D)
        @test fwd[j] === false
        @test all(fwd[k] for k in eachindex(fwd) if k != j)
        # The stored forward rate D → E runs against the solve's direction.
        @test EnzymeRates._priority_key(first(flat[j]), false, 1, 0)[5] === true
        @test EnzymeRates._form_index_key(EnzymeRates.to_species(first(flat[j]))) <
              EnzymeRates._form_index_key(EnzymeRates.from_species(first(flat[j])))
    end

    @testset "fundamental cycles carry their multiple of the net reaction" begin
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A ⇌ E(A)
                E + B ⇌ E(B)
                E(A) + B ⇌ E(A, B)
                E(B) + A ⇌ E(A, B)
                E(A, B) <--> E(P)
                E(P) ⇌ E + P
            end
        end)
        Z, multiples = EnzymeRates._fundamental_cycles(EnzymeRates.reaction(m),
                                                       EnzymeRates.steps(m))
        @test size(Z) == (2, 6)
        @test sort(abs.(multiples)) == [0, 1]
        # A cycle that turns A into P without B and Q is no multiple of the net reaction.
        slip = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) <--> E + Q
                E(A) <--> E(P)
                E(P) <--> E + P
            end
        end)
        err = @test_throws ErrorException EnzymeRates._fundamental_cycles(
            EnzymeRates.reaction(slip), EnzymeRates.steps(slip))
        @test occursin("produces metabolite change not proportional to net reaction",
                       err.value.msg)
    end

    @testset "the reduction runs in Int64 and retries in BigInt on overflow" begin
        M, pivots = EnzymeRates._reduce_joint_system([1 -1; 1 1], [1, 0])
        @test eltype(M) == Rational{Int} && pivots == [1, 2]
        A = [2^40 1; 3 2^40]
        M, pivots = EnzymeRates._reduce_joint_system(A, [1, 0])
        @test eltype(M) == Rational{BigInt} && pivots == [1, 2]
        @test Rational{BigInt}.(A) * M[:, 3] == [1, 0]
        err = @test_throws ErrorException EnzymeRates._reduce_joint_system(
            reshape([1, 1], 2, 1), [0, 1])
        @test occursin("Thermodynamically contradictory mechanism", err.value.msg)
    end

    @testset "a class keeps its last column, a member of the inactive law first" begin
        # I binds E(A) and E(P) under one :NonequalAI group. In the inactive
        # conformation E(A) and E(A, I) lie behind the :OnlyA binding of A, so only the
        # E(P) + I member is in that conformation's law.
        am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: A
            products: P
            catalytic_inhibitors: I
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + A <--> E(A)                                   :: OnlyA
                E(A) <--> E(P)                                    :: OnlyA
                E + P ⇌ E(P)                                      :: EqualAI
                (E(A) + I ⇌ E(A, I), E(P) + I ⇌ E(I, P))          :: NonequalAI
            end
        end)
        kept = EnzymeRates._kept_members(am)
        g = only(g for (g, group) in enumerate(EnzymeRates.steps(am)) if length(group) == 2)
        # Active: E(A) + I sorts before E(P) + I (A before P), so it is kept; inactive:
        # the member in the law is kept.
        @test EnzymeRates._step_text(kept[(g, 1, false)]) == "EA + I ⇌ EAIinh"
        @test EnzymeRates._step_text(kept[(g, 1, true)]) == "EP + I ⇌ EIinhP"
        sys = EnzymeRates._joint_system(am)
        for (key, step) in kept
            last_col = findlast(c -> c.slot == key[2] && c.inactive == key[3] &&
                                     sys.flat[c.step][2] == key[1], sys.columns)
            @test sys.flat[sys.columns[last_col].step][1] == step
        end
    end

    @testset "every spec's solve counts its declared dimension" begin
        for spec in MECHANISM_TEST_SPECS
            mech = _testhelper_lift(spec.mechanism)
            sys = EnzymeRates._joint_system(mech)
            _, pivots = EnzymeRates._reduce_joint_system(sys.A, sys.rhs)
            class(c) = EnzymeRates._class(c, sys.flat[c.step][2], sys.tags[sys.flat[c.step][2]])
            last_of = Dict(class(c) => i for (i, c) in enumerate(sys.columns))
            n_law = count(i -> (c = sys.columns[i]; (!c.inactive || sys.law[c.step]) &&
                                                    last_of[class(c)] == i &&
                                                    i ∉ pivots), eachindex(sys.columns))
            n_reg = mech isa EnzymeRates.AllostericMechanism ?
                    sum((t === :NonequalAI ? 2 : 1
                         for site in EnzymeRates.regulatory_sites(mech)
                         for t in EnzymeRates.allo_states(site)); init = 0) + 1 : 0
            @test n_law + n_reg == _testhelper_declared_dimension(spec.mechanism)
        end
    end

    @testset "every spec's solve satisfies detailed balance in the oracle" begin
        rng = Random.MersenneTwister(11)
        for spec in MECHANISM_TEST_SPECS
            mech = _testhelper_lift(spec.mechanism)
            point = setprecision(() -> _testhelper_solve_point(mech, rng), BigFloat, 256)
            rxn = EnzymeRates.reaction(mech)
            subs = EnzymeRates.name.(EnzymeRates.substrates(rxn))
            prods = EnzymeRates.name.(EnzymeRates.products(rxn))
            concs = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1)
                                           for x in EnzymeRates.metabolites(spec.mechanism))
            p1 = first(prods)
            concs[p1] = point.Keq * prod(concs[x] for x in subs) /
                        prod((concs[x] for x in prods if x != p1); init = big(1))
            off = copy(concs)
            off[p1] /= 10
            v_eq = _testhelper_steady_state_rate(mech, point, concs)
            v_off = _testhelper_steady_state_rate(mech, point, off)
            @test abs(v_eq) <= big(10)^-30 * abs(v_off)
        end
    end
end
```

In `test/runtests.jl`, add `include("test_joint_solve.jl")` on the line after
`include("test_oracles.jl")`.

- [ ] **Step 2: Run the tests to verify they fail**

Run the focused command for `test_joint_solve`.
Expected: errors with `UndefVarError: _joint_system not defined` (and the other new
names).

- [ ] **Step 3: Factor the orientation helpers in `src/types.jl`**

1. In `_orient_tied_steps`, change the signature to

```julia
function _orient_tied_steps(oriented::Vector{Step}, tiers::Vector{Int};
                            lead_key = _step_canonical_key)
```

   and in its `lead = …` line replace `_step_canonical_key(oriented[i])` with
   `lead_key(oriented[i])`. In its docstring replace
   ``first), then by deciding tier, then by `_step_canonical_key`. Every other tied step``
   with the two lines
   ``first), then by deciding tier, then by `lead_key` (by default `_step_canonical_key`;``
   ``the joint solve passes a content key, `_solve_forward`). Every other tied step``.
2. In `_canonicalize_step_directions`, replace the loop that fills `entry_sides` (from
   `    entry_sides = Dict{Species, Tuple{Bool, Bool}}()` through its `end`) with
   `    entry_sides = _entry_sides(groups)`, and in its docstring replace
   ``The table `entry_sides` built here classifies each species`` with
   ``The table `entry_sides` (`_entry_sides`) classifies each species``.
3. Immediately before the comment line ``# Canonical key for a `Step`. Gives the sorts``
   add:

```julia
"""
Each enzyme form of `groups` with the metabolites that enter or leave solution at it, as
`(any substrate, any product)`: the consumed metabolites of every step leaving it and the
released metabolites of every step arriving at it (`_canonicalize_step_directions`).
"""
function _entry_sides(groups::Vector{Vector{Step}})
    entry_sides = Dict{Species, Tuple{Bool, Bool}}()
    for group in groups, s in group, (form, free) in ((from_species(s), consumed(s)),
                                                      (to_species(s), released(s)))
        has_sub, has_prod = get(entry_sides, form, (false, false))
        entry_sides[form] = (has_sub  || any(m -> m isa Substrate, free),
                             has_prod || any(m -> m isa Product, free))
    end
    entry_sides
end

```

- [ ] **Step 4: Write `src/thermo.jl`**

Create `src/thermo.jl`:

```julia
# ABOUTME: The joint thermodynamic solve: one exact row reduction over the constants of
# ABOUTME: every step in both conformations, with a fixed priority choosing the dependents.

"""
A bare constant of the joint solve: slot `slot` of flat step `step` (`_flat_steps`
order), in the inactive conformation when `inactive`. Slot 1 is a rapid-equilibrium
step's equilibrium constant or a steady-state step's rate constant in its stored
direction; slot 2 is a steady-state step's rate constant in the reverse direction.
"""
struct _Column
    step::Int
    slot::Int
    inactive::Bool
end

# ─── Content keys and orientation ──────────────────────────────────────────

"""A metabolite list as the content keys read it: each metabolite's name and whether it
binds as a modifier copy (a competitive inhibitor), in stored order."""
_metabolites_key(mets) = Tuple((String(name(x)), x isa CompetitiveInhibitor) for x in mets)

"""
An enzyme form's content key, built from base names: its conformation label, its bound
metabolites (`_metabolites_key`) and the names its residual adds and subtracts. Rendered
names never enter it.
"""
_form_base_key(sp::Species) =
    (String(conformation(sp)), _metabolites_key(bound(sp)),
     Tuple(String(name(x)) for x in added(residual(sp))),
     Tuple(String(name(x)) for x in subtracted(residual(sp))))

"""Whether `sp` is free E: conformation label `E`, nothing bound, no residual."""
_is_free_E(sp::Species) =
    conformation(sp) === :E && isempty(bound(sp)) && isempty(residual(sp))

"""
The canonical order of enzyme forms that orients the steps no chemistry orients (C9):
free E first, then by `_form_base_key`.
"""
_form_index_key(sp::Species) = (!_is_free_E(sp), _form_base_key(sp))

"""A step in a given direction as the canonical form order reads it: its two forms
(`_form_index_key`), then the metabolites it takes up and gives off."""
_step_index_key(s::Step) =
    (_form_index_key(from_species(s)), _form_index_key(to_species(s)),
     _metabolites_key(consumed(s)), _metabolites_key(released(s)))

"""
Whether each flat step of `groups`, as stored, runs in the joint solve's forward
direction. Chemistry (orientation tier 1) and graph context (tier 2) orient a step as
the `Mechanism` constructor does; a step neither orients (tier 3) runs from the form
earlier in the canonical form order (`_form_index_key`) to the later one, and the tied
steps of a kinetic group then turn with a lead ranked by `_step_index_key`
(`_orient_tied_steps`). The constructor breaks those last two ties by rendered name, so
for such a step the stored direction can differ from the solve's.
"""
function _solve_forward(groups::Vector{Vector{Step}})
    sides = _entry_sides(groups)
    fwd = Bool[]
    for group in groups
        decided = [_canonical_step_direction(s, sides) for s in group]
        tiers = [last(d) for d in decided]
        by_index = Step[t == 3 ? (_form_index_key(from_species(s)) <=
                                  _form_index_key(to_species(s)) ? s : _reversed(s)) : d
                        for (s, (d, t)) in zip(group, decided)]
        aligned = _orient_tied_steps(by_index, tiers; lead_key = _step_index_key)
        append!(fwd, [a == s for (a, s) in zip(aligned, group)])
    end
    fwd
end

"""
The priority key of constant `slot` of step `s`, whose stored direction is the solve's
forward one when `fwd`, in conformation block `block` (2: an inactive copy outside the
law, 1: another inactive copy, 0: active). The solve eliminates the larger key first, so
a column with a larger key becomes dependent sooner. In turn: the block; the step's kind
(3 an isomerization, 2 other chemistry, 1 a plain binding); whether both of its forms hold
a metabolite; steady state before rapid equilibrium; the constant that runs against the
reaction (the reverse of the solve's forward direction, except that for a product binding
the uptake runs against it); then the step's content in the solve's direction, its forms'
`_form_base_key`s and its metabolites, the smaller key kept.
"""
function _priority_key(s::Step, fwd::Bool, slot::Int, block::Int)
    kind = is_iso(s) ? 3 : _is_chemistry(s) ? 2 : 1
    occupied = !isempty(bound(from_species(s))) && !isempty(bound(to_species(s)))
    ss = !is_equilibrium(s)
    against = ss && ((slot == 2) == fwd) != (bound_metabolite(s) isa Product)
    a, b, c, r = fwd ? (from_species(s), to_species(s), consumed(s), released(s)) :
                       (to_species(s), from_species(s), released(s), consumed(s))
    (block, kind, occupied, ss, against,
     (_form_base_key(a), _form_base_key(b), _metabolites_key(c), _metabolites_key(r)))
end

# ─── The joint system ──────────────────────────────────────────────────────

"""Each kinetic group's allosteric tag, `:None` for a plain mechanism."""
_solve_tags(m::Mechanism) = fill(:None, length(steps(m)))
_solve_tags(am::AllostericMechanism) = cat_allo_states(am)

"""Whether each flat step lies in the inactive conformation's law (`_inactive_groups`);
false for every step of a plain mechanism."""
_inactive_law(m::Mechanism) = falses(length(_flat_steps(m)))
function _inactive_law(am::AllostericMechanism)
    kept = _inactive_groups(am)
    BitVector([s in kept[g] for (s, g) in _flat_steps(steps(am))])
end

"""Whether the steps of a group with allosteric tag `tag` take their own constants in the
inactive conformation: those of a `:NonequalAI` or `:OnlyA` group do; an `:EqualAI` group's
read the active ones."""
_has_inactive_columns(tag::Symbol) = tag in (:NonequalAI, :OnlyA)

"""
The kinetic class of column `c`, whose step lies in group `g` with tag `tag`: the group's
constant in the column's slot and conformation, except that the inactive copy of each
`:OnlyA` step is a class of its own.
"""
_class(c::_Column, g::Int, tag::Symbol) =
    c.inactive && tag === :OnlyA ? (c.step, c.slot, true, true) :
                                   (g, c.slot, c.inactive, false)

"""
    _fundamental_cycles(rxn, groups) -> (Z, multiples)

The fundamental cycles of a spanning forest of the graph whose nodes are the enzyme forms
of `groups` and whose edges are its flat steps: row `i` of `Z` counts each flat step of
cycle `i` with its sign in the step's stored direction, and `multiples[i]` is the multiple
of the net reaction the cycle runs (0 for a Wegscheider cycle). A metabolite counts by
name, so a reactant's modifier copy counts as the reactant. Errors on a cycle whose net
change is no multiple of the net reaction. Integer arithmetic throughout.
"""
function _fundamental_cycles(rxn::EnzymeReaction, groups::Vector{Vector{Step}})
    flat = _flat_steps(groups)
    forms = _forms(groups)
    idx = Dict(f => i for (i, f) in enumerate(forms))
    ends = [(idx[from_species(s)], idx[to_species(s)]) for (s, _) in flat]
    path = [zeros(Int, length(flat)) for _ in forms]   # signed tree steps from the root
    seen, tree = falses(length(forms)), falses(length(flat))
    for root in eachindex(forms)
        seen[root] && continue
        seen[root] = grown = true
        while grown
            grown = false
            for (k, (a, b)) in enumerate(ends)
                seen[a] == seen[b] && continue
                from, to, sgn = seen[a] ? (a, b, 1) : (b, a, -1)
                path[to] .= path[from]
                path[to][k] += sgn
                seen[to] = tree[k] = grown = true
            end
        end
    end
    net = Dict{Symbol, Int}()
    foreach(x -> net[name(x)] = -1, substrates(rxn))
    foreach(x -> net[name(x)] = 1, products(rxn))
    x1 = name(first(substrates(rxn)))
    rows, multiples = Vector{Int}[], Int[]
    for (k, (a, b)) in enumerate(ends)
        tree[k] && continue
        z = path[a] - path[b]
        z[k] += 1
        change = Dict{Symbol, Int}(x => 0 for x in keys(net))
        for (j, (s, _)) in enumerate(flat)
            z[j] == 0 && continue
            for x in released(s); haskey(net, name(x)) && (change[name(x)] += z[j]); end
            for x in consumed(s); haskey(net, name(x)) && (change[name(x)] -= z[j]); end
        end
        c = change[x1] * net[x1]
        all(x -> change[x] == c * net[x], keys(net)) || error(
            "Cycle $(length(rows) + 1) produces metabolite change not proportional to " *
            "net reaction")
        push!(rows, z)
        push!(multiples, c)
    end
    Z = isempty(rows) ? zeros(Int, 0, length(flat)) : permutedims(reduce(hcat, rows))
    Z, multiples
end

"""
    _joint_system(m) -> (; columns, A, rhs, law, tags, flat)

The thermodynamic system of `m`, `A · log(constants) = rhs · log(Keq)`, over one column
per bare constant (`_Column`): one per slot of every flat step, and for an allosteric
mechanism one more per slot of every step of a `:NonequalAI` or `:OnlyA` group, for its
inactive copy. `columns` lists them in elimination order: by `_priority_key`, larger
first, the steps oriented by `_solve_forward`.

Rows: every fundamental cycle (`_fundamental_cycles`) once per conformation, the inactive
rows over the whole graph (islands included) with an `:EqualAI` step reading its active
column; a cycle row sums, over its steps, the log of each step's forward equilibrium ratio
(`kf/kr` at steady state, `1/K` for a binding's dissociation constant, `K` otherwise), and
its right side is the cycle's multiple of the net reaction. Then one row per further
member of each kinetic class (`_class`), setting it equal to the class's first member,
slot by slot. `law[j]` tells whether flat step `j` lies in the inactive conformation's
law, `tags` holds each group's tag and `flat` the flat steps.
"""
function _joint_system(m::_AnyMech)
    groups = steps(m)
    flat = _flat_steps(groups)
    tags = _solve_tags(m)
    law = _inactive_law(m)
    fwd = _solve_forward(groups)
    allosteric = m isa AllostericMechanism
    columns = _Column[]
    for (j, (s, g)) in enumerate(flat), inactive in (false, true)
        inactive && !(allosteric && _has_inactive_columns(tags[g])) && continue
        for slot in 1:(is_equilibrium(s) ? 1 : 2)
            push!(columns, _Column(j, slot, inactive))
        end
    end
    block(c) = !c.inactive ? 0 : law[c.step] ? 1 : 2
    keys = [_priority_key(flat[c.step][1], fwd[c.step], c.slot, block(c)) for c in columns]
    columns = columns[sortperm(keys; rev = true)]
    index = Dict(c => i for (i, c) in enumerate(columns))
    column(j, slot, inactive) =
        index[_Column(j, slot, inactive && _has_inactive_columns(tags[flat[j][2]]))]
    Z, multiples = _fundamental_cycles(reaction(m), groups)
    rows, rhs = Vector{Int}[], Int[]
    for inactive in (allosteric ? (false, true) : (false,)), i in eachindex(multiples)
        row = zeros(Int, length(columns))
        for (j, (s, _)) in enumerate(flat)
            z = Z[i, j]
            z == 0 && continue
            row[column(j, 1, inactive)] += is_equilibrium(s) && is_binding(s) ? -z : z
            is_equilibrium(s) || (row[column(j, 2, inactive)] -= z)
        end
        push!(rows, row)
        push!(rhs, multiples[i])
    end
    first_of = Dict{NTuple{4, Any}, Int}()
    for (i, c) in enumerate(columns)
        g = flat[c.step][2]
        f = get!(first_of, _class(c, g, tags[g]), i)
        f == i && continue
        row = zeros(Int, length(columns))
        row[i], row[f] = 1, -1
        push!(rows, row)
        push!(rhs, 0)
    end
    A = isempty(rows) ? zeros(Int, 0, length(columns)) : permutedims(reduce(hcat, rows))
    (; columns, A, rhs, law, tags, flat)
end

# ─── Exact reduction ───────────────────────────────────────────────────────

"""
Bring `M` to reduced row echelon form in place over its first `n` columns and return the
pivot columns, in row order. `Rational{Int}` arithmetic throws `OverflowError` rather than
wrap, which `_reduce_joint_system` catches.
"""
function _rref!(M::Matrix{<:Rational}, n::Int)
    pivots = Int[]
    for col in 1:n
        r = length(pivots) + 1
        r > size(M, 1) && break
        p = findfirst(i -> !iszero(M[i, col]), r:size(M, 1))
        p === nothing && continue
        p += r - 1
        for k in axes(M, 2)
            M[r, k], M[p, k] = M[p, k], M[r, k]
        end
        f = inv(M[r, col])
        for k in axes(M, 2)
            M[r, k] *= f
        end
        for i in axes(M, 1)
            (i == r || iszero(M[i, col])) && continue
            g = M[i, col]
            for k in axes(M, 2)
                M[i, k] -= g * M[r, k]
            end
        end
        push!(pivots, col)
    end
    pivots
end

"""
    _reduce_joint_system(A, rhs) -> (M, pivots)

The reduced row echelon form `M` of `[A rhs]` over `A`'s columns and its pivot columns,
in checked `Rational{Int}` arithmetic, redone in `Rational{BigInt}` on overflow. A pivot
column is a dependent constant: row `r` of `M` reads `x_p + Σ_f M[r, f]·x_f =
M[r, end]·log Keq` over the columns `f` that are not pivots. Errors when the rows combine
to `0 = log Keq`.
"""
function _reduce_joint_system(A::Matrix{Int}, rhs::Vector{Int})
    n = size(A, 2)
    M = Matrix{Rational{Int}}([A rhs])
    pivots = try
        _rref!(M, n)
    catch err
        err isa OverflowError || rethrow()
        M = Matrix{Rational{BigInt}}([A rhs])
        _rref!(M, n)
    end
    all(iszero, view(M, (length(pivots) + 1):size(M, 1), n + 1)) || error(
        "Thermodynamically contradictory mechanism: a combination of its " *
        "constraint rows reduces to 0 = log(Keq)")
    M, pivots
end

# ─── Class names ───────────────────────────────────────────────────────────

"""
The kept member of each kinetic class of `m`, keyed by `(group, slot, inactive)`: the
member with the smallest `_priority_key`, which comes last in the joint solve's
elimination order, so the solve keeps it whenever it keeps one of the class (E3). A
`:NonequalAI` group has an inactive class, kept by a member in the inactive law when it has
one; an `:EqualAI` group's steps share their active class, and each `:OnlyA` step's
inactive copy is a class of its own, so neither appears.
"""
function _kept_members(m::_AnyMech)
    flat = _flat_steps(steps(m))
    tags = _solve_tags(m)
    law = _inactive_law(m)
    fwd = _solve_forward(steps(m))
    best = Dict{Tuple{Int, Int, Bool}, Tuple{Any, Step}}()
    for (j, (s, g)) in enumerate(flat), inactive in (false, true)
        inactive && tags[g] !== :NonequalAI && continue
        for slot in 1:(is_equilibrium(s) ? 1 : 2)
            key = _priority_key(s, fwd[j], slot, !inactive ? 0 : law[j] ? 1 : 2)
            cur = get(best, (g, slot, inactive), nothing)
            (cur === nothing || key < first(cur)) && (best[(g, slot, inactive)] = (key, s))
        end
    end
    Dict{Tuple{Int, Int, Bool}, Step}(k => last(v) for (k, v) in best)
end
```

In `src/EnzymeRates.jl`, add `include("thermo.jl")` on the line after
`include("thermodynamic_constr_for_rate_eq_derivation.jl")`.

- [ ] **Step 5: Run the tests to verify they pass**

Run the focused command for `test_joint_solve`.
Expected: `Joint thermodynamic solve | 113 113`. Then run the focused command for
`test_types` (the orientation refactor): all pass.

- [ ] **Step 6: Run the full suite, then commit**

Run the full suite. Expected: all pass (JET checks the new file).

```bash
git status
git add src/thermo.jl src/EnzymeRates.jl src/types.jl test/test_joint_solve.jl \
        test/runtests.jl
git commit -m "Add the joint thermodynamic solve over both conformations" \
  -m "One integer system holds every step constant of both conformations: the fundamental
cycles once per conformation and one row per further member of each kinetic class. A
fixed priority (block, kind, occupancy, flag, the constant against the reaction, then a
base-name content key) orders the columns, steps no chemistry orients take the canonical
form order, and an exact Int64 reduction, retried in BigInt on overflow, picks the
dependents. Nothing calls it yet; its constants satisfy detailed balance in the
mass-action oracle on every spec and its counts equal each spec's declared dimension." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---
### Task 7: Switch the derivation to the joint solve

**Files:**
- Modify: `src/thermo.jl` (append the dependents)
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl` (delete the old solve)
- Modify: `src/types.jl` (the E3 chokepoint; `_assert_bottomed_re_segments`)
- Modify: `src/rate_eq_derivation.jl` (no rename pass; `_step_names`, `_state_groups`;
  the inactive law on `_inactive_groups`; B10)
- Modify: `src/mechanism_enumeration.jl`, `src/identify_rate_equation.jl`
- Modify: `test/test_oracles.jl` (population laws), and re-record:
  `test/mechanism_definitions_for_test_enzyme_derivation.jl`,
  `test/test_rate_eq_derivation.jl`, `test/test_types.jl`,
  `test/test_allosteric_collapse.jl`, `test/test_allosteric_golden.jl`,
  `test/reference/allosteric_golden_reference.txt`, `test/allosteric_ground_truth.jl`,
  `test/test_fitting_scale.jl`, `test/test_mechanism_enumeration.jl`,
  `test/test_identify_rate_equation.jl`, `test/test_pivot_priority_regression.jl`

**Interfaces:**
- Consumes: everything Task 6 produces; Task 4's engine; Tasks 1–2's oracle helpers.
- Produces:
  - `_dependent_param_exprs(m::_AnyMech) -> (dep::Dict{Symbol, Union{Symbol, Expr,
    Float64}}, indep::Tuple)` and its `Type{<:AbstractEnzymeMechanism}` method, in
    `src/thermo.jl`; `_independent_param_count(m::_AnyMech)`; `_column_parameter`
  - `name(p::Kfor|Krev|Krapid, m)` renders the kept member of `p`'s class
    (`_class_sides(m)`, `_constant_sides(step, slot, state, m)`)
  - `_state_groups(m, state)`, `_step_names(m::_AnyMech, state::Symbol) ->
    Vector{Vector{Symbol}}` (`state` is `:None` for a plain mechanism, `:A` or `:I`)
  - `_state_rate_polys(am, state)` on `_state_groups(am, state)`;
    `_assert_bottomed_re_segments(reaction, groups)`
  - `_kreg_params(am, :A)` lists an `:EqualAI` ligand once, as `Kreg(site, lig,
    :EqualAI)`; `_reg_K` names it bare in both conformations (B10)
  - removed: `ANNOTATION_SUBSTITUTED`, `_build_wegscheider_rename_map`, `_leaf_names`,
    `_step_parameters`, `_free_enz_set`, `_step_priority`, `_group_rep`,
    `_group_reps`, `_rep_sides`, `_enumerate_parameters_full`, `_param_columns`,
    `_assemble_constraints`, `_solve_dependent_set`, `_state_allo_mechanism`,
    `_state_parts`, `_state_mechanism`, `_cat_params`

This task switches every law to the joint solve in one step and re-records the pins it
moves. A prototype of exactly this task (the same code, the same test edits) passed the
full suite at 33,299 tests in 6m54s, Aqua and JET included; `rate_equation` ran 77.3,
116.6 and 103.2 ns on the three largest bi-bi laws (baseline 77.5, 118.6, 103.8).

What changes, and why every changed pin is the same law:

- **Plain mechanisms.** Counts never move (measured on 20,867 enumerated mechanisms),
  but the dependents follow the new priority: 20 of the 42 specs and most enumerated
  mechanisms keep a different set fitted. The textbook case: in Segel's ordered bi-bi
  the product uptake `k_EQ_P_to_EAB` (it runs against the reaction, E7) becomes the
  Haldane dependent and the free-enzyme release `k_EA_to_E_A` is fitted. The 42 specs'
  `Mass-action oracle`, analytical-rate, QSSA, ODE and string-evaluation tests check the
  same law values under the new names.
- **No rename pass.** A single-symbol tie (`K_EAB_to_EA_B = K_EB_to_E_B`) is an ordinary
  dependent: it prints under `# Wegscheider constraints:` and the law reads the
  dependent class's own name. This moves the `expected_n_mirror_constraints` count of
  three inhibitor specs and their factored strings; the printed law of a same-model
  regrouping now differs (its `eq_hash` too), while its fitted names stay (E3).
- **Allosteric mechanisms.** The inactive law is derived on `_inactive_groups(am)` with
  `am`'s own steps, orientation and names, and both conformations' cycle rows sit in one
  system. ROADMAP known issues 1 and 2 go: 71 of the 2,344 mechanisms of
  `_testhelper_population(_testhelper_pingpong_allo_rxn, 1)` (818 of 12,985 after two
  rounds) lose one or two phantom or folded parameters and now equal the declared
  dimension. `:EqualAI` regulator ligands print bare (`K_Pireg`) with no mirror
  dependent (B10).
- **Enumeration.** The allosteric split count is the joint solve's and derives no rate
  polynomial, so a law over `MAX_RATE_EQUATION_TERMS` no longer refuses a count; the
  count also drops the phantoms above, so the split move's children change for parents
  that carried them. Among the enumeration tests only the two term-cap testsets move: a
  parent over the cap now splits by its count (two new children, written out below).
- **identify.** `eq_hash` keeps its role, recomputed from the new printed law; the
  dedup key strips only `#` lines. A same-model regrouping no longer prints the same law,
  so the twin fixture of the batch dedup test becomes a mechanism that compiles to the
  same type (an unbound inhibitor), and the allosteric split/merge pair, whose extra
  parameter was a fold hazard (known issue 2), becomes one model with one fitted set.

- [ ] **Step 1: Write the failing population tests**

Append to `test/test_oracles.jl`:

```julia

@testset "Laws over a regenerated population" begin
    @testset "every sampled law matches the exact oracle" begin
        for (i, m) in enumerate(_testhelper_population_sample)
            _testhelper_check_law_against_oracle(EnzymeRates.compile_mechanism(m);
                                                 n_points = 2, seed = i, compiled = false)
        end
    end

    @testset "fitted counts equal the declared dimension" begin
        population = _testhelper_population(_testhelper_pingpong_allo_rxn, 1)
        @test all(m -> EnzymeRates._independent_param_count(m) ==
                       _testhelper_declared_dimension(m), population)
        for m in _testhelper_population_sample[1:10]
            @test length(EnzymeRates.fitted_params(m)) == _testhelper_declared_dimension(m)
        end
    end

    @testset "every fitted catalytic name names a constant of the law" begin
        # Known issues 1 and 2 of the ROADMAP gave an :EqualAI constant a second fitted
        # name, the name of no step of the law.
        for m in _testhelper_population_sample
            m isa EnzymeRates.AllostericMechanism || continue
            law_names = Set(n for state in (:A, :I)
                            for ns in EnzymeRates._step_names(m, state) for n in ns)
            reg = Set(EnzymeRates.name(p, m) for state in (:A, :I)
                      for p in EnzymeRates._kreg_params(m, state))
            fitted = EnzymeRates._dependent_param_exprs(m)[2]
            @test issubset(setdiff(fitted, reg, (:L,)), law_names)
        end
    end

    @testset "a same-model regrouping keeps every fitted name" begin
        n_splits = 0
        for m in _testhelper_population_sample
            m isa EnzymeRates.Mechanism || continue
            counter = EnzymeRates._partition_independent_count(m)
            flat = EnzymeRates._flat_steps(m)
            position = Dict(s => j for (j, (s, _)) in enumerate(flat))
            ids = [g for (_, g) in flat]
            names = Set(EnzymeRates._dependent_param_exprs(m)[2])
            for (g, group) in enumerate(EnzymeRates.steps(m)),
                bp in EnzymeRates._context_bipartitions(group)
                split_ids = copy(ids)
                for s in bp[2]; split_ids[position[s]] = length(EnzymeRates.steps(m)) + 1; end
                counter(split_ids) == counter(ids) || continue
                child = EnzymeRates._apply_bipartitions(m, [(g, bp)])
                n_splits += 1
                @test Set(EnzymeRates._dependent_param_exprs(child)[2]) == names
            end
        end
        @test n_splits >= 20
    end

    @testset "a :NonequalAI group behind an :OnlyA binding takes no inactive name" begin
        # E(S) and E(I, S) lie behind the :OnlyA binding of S, so the inactive
        # conformation never reaches the I binding: its inactive class is outside the
        # law, and no K_I_ constant is fitted.
        island = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_inhibitors: I
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + S <--> E(S)          :: OnlyA
                E(S) <--> E(P)           :: OnlyA
                E + P ⇌ E(P)             :: EqualAI
                E(S) + I ⇌ E(I, S)       :: NonequalAI
            end
        end
        @test EnzymeRates.fitted_params(island) ==
              (:K_A_EIinhS_to_ES_Iinh, :K_EP_to_E_P, :k_A_ES_to_EP, :k_A_ES_to_E_S,
               :k_A_E_S_to_ES, :L)
        @test length(EnzymeRates.fitted_params(island)) ==
              _testhelper_declared_dimension(island)
        _testhelper_check_law_against_oracle(island)
    end
end
```

- [ ] **Step 2: Run them to verify they fail**

Run the focused command for `test_oracles`.
Expected: `fitted counts equal the declared dimension` fails (71 mechanisms carry a
phantom or folded parameter today), `every fitted catalytic name names a constant of the
law` errors (`_step_names` does not exist yet), `a same-model regrouping keeps every
fitted name` fails on many splits (today's names follow the grouping), and the oracle
and island testsets pass.

- [ ] **Step 3: Add the dependents to `src/thermo.jl`**

Append to `src/thermo.jl`:

```julia

"""The parameter of column `c` of `m`, whose flat steps are `flat`, in the state its name
takes: `:None` in a plain mechanism; `:I` for an inactive copy; `:EqualAI` or `:A` for an
active constant of an allosteric mechanism."""
function _column_parameter(m::_AnyMech, flat, c::_Column)
    s, g = flat[c.step]
    state = m isa Mechanism ? :None : c.inactive ? :I :
            cat_allo_state(m, g) === :EqualAI ? :EqualAI : :A
    c.slot == 2 ? Krev(s, state) : is_equilibrium(s) ? Krapid(s, state) : Kfor(s, state)
end

# ─── Dependents ────────────────────────────────────────────────────────────

"""
    _dependent_param_exprs(m) -> (dep, indep)

The dependent constants of `m`'s law and its independent ones, from one joint solve
(`_joint_system`, `_reduce_joint_system`). Every class of constants of the law is
independent or dependent, under the name of its kept member (`_kept_members`):

- `indep` holds the classes whose kept column is not a pivot, sorted by name (`isless`
  on the strings, so `K_…` precedes `k_…`); an allosteric mechanism's list continues with
  its regulator constants (`_kreg_params`: active, then inactive, site by site) and `L`.
- `dep` maps each class whose kept column is a pivot to `Keq^q · Π independent^e` with
  rational exponents (`build_power_expr`); the empty product is `1.0`.

Columns outside the law (an `:OnlyA` step's inactive copy, an inactive copy on a pruned
form) are eliminated first, so no dependent of the law reads them; a pivot that is not
its class's kept member equals the kept one and adds no entry. The `Type` method lifts a
compiled mechanism with `_concrete`.
"""
function _dependent_param_exprs(m::_AnyMech)
    sys = _joint_system(m)
    M, pivots = _reduce_joint_system(sys.A, sys.rhs)
    n = length(sys.columns)
    names = [name(_column_parameter(m, sys.flat, c), m) for c in sys.columns]
    class(c) = _class(c, sys.flat[c.step][2], sys.tags[sys.flat[c.step][2]])
    last_of = Dict(class(c) => i for (i, c) in enumerate(sys.columns))
    kept_in_law(i) = (c = sys.columns[i];
                      (!c.inactive || sys.law[c.step]) && last_of[class(c)] == i)
    free = setdiff(1:n, pivots)
    dep = Dict{Symbol, Union{Symbol, Expr, Float64}}()
    for (r, p) in enumerate(pivots)
        kept_in_law(p) || continue
        factors = [(names[f], -M[r, f]) for f in free if !iszero(M[r, f])]
        dep[names[p]] = build_power_expr(M[r, n + 1], factors)
    end
    indep = sort!([names[f] for f in free if kept_in_law(f)]; by = string)
    m isa Mechanism && return dep, Tuple(indep)
    reg = Symbol[name(p, m) for p in [_kreg_params(m, :A); _kreg_params(m, :I)]]
    dep, Tuple([indep; reg; :L])
end

_dependent_param_exprs(@nospecialize(M::Type{<:AbstractEnzymeMechanism})) =
    _dependent_param_exprs(_concrete(M()))

"""Number of independent (fitted) constants of a concrete mechanism, from the joint solve,
without compiling the mechanism."""
_independent_param_count(m::_AnyMech) = length(_dependent_param_exprs(m)[2])
```

- [ ] **Step 4: Delete the old solve**

In `src/thermodynamic_constr_for_rate_eq_derivation.jl`:

1. Replace the two ABOUTME lines and the module docstring under them with:

```julia
# ABOUTME: Cycle-space helpers for the :OnlyA guard and the split move's counter:
# ABOUTME: the exact null space, the cycle basis and the parameter count of a regrouping.
```

2. Delete the section `# ─── Shared Helpers ───…` with `_step_parameters`, and the
   section `# ─── Structural primacy: free-enzyme set + step priority ───…` with
   `_compute_free_enz_set`, `_free_enz_set`, `_step_priority` and `_group_rep` (all
   docstrings included), up to the line `# ─── Thermodynamic Constraint Infrastructure ───…`.
3. Replace the docstring of `_rref_partition` with:

```julia
"""
Reduced row echelon form over `Rational{BigInt}`. Returns the pivot and
free column indices (pivot_cols in row-pivot order, so pivot_cols[i] is the
pivot at reduced-matrix row i) plus the reduced matrix R. Used by
`_rational_nullspace` (nullspace basis) and `_partition_independent_count` (rank).
"""
```
4. Delete `_dependent_param_exprs(mech::Mechanism)` and `_independent_param_count` with
   their docstrings (`src/thermo.jl` holds both now), `_param_columns` and
   `_assemble_constraints` with their docstrings, and `_solve_dependent_set` with its
   docstring.
5. In the docstring of `_partition_independent_count`, replace the sentences from
   ``Equals `_independent_param_count` of the`` to the end of the docstring with:

```text
Equals `_independent_param_count` of the
constructed child: the joint solve's independent set is its columns of the law less
its pivots, and a kinetic group's constants form one column here and one class there.
```

6. In the docstring of `_onlya_haldane_violation`, replace
   ``not call `_state_allo_mechanism`, which would construct an `AllostericMechanism` and``
   ``recurse.`` (two lines) with
   ``not construct an `AllostericMechanism`, whose constructor calls it.``

- [ ] **Step 5: Name each class after its kept member in `src/types.jl`**

1. Replace the `_NamingCache` docstring and struct with:

```julia
"""
Naming data derived from a mechanism's steps, filled on first use: the rendered sides
of each kinetic class's kept member (`_class_sides`). It never takes part in `==`,
`hash`, the compiled type or display.
"""
mutable struct _NamingCache
    sides::Union{Nothing, Dict{Tuple{Int, Int, Bool}, Tuple{String, String}}}
end
_NamingCache() = _NamingCache(nothing)
```

2. Replace `_group_reps`, `_group_sides`, `_rep_sides`, the comment block above the
   `name` methods and the three methods `name(p::Kfor, …)`, `name(p::Krev, …)`,
   `name(p::Krapid, …)` (lines 1408–1444) with:

```julia
"""The `_forward_sides` of each kinetic class's kept member (`_kept_members`), keyed by
`(group, slot, inactive)`; cached per mechanism."""
function _class_sides(m::_AnyMech)
    c = m.naming
    sides = c.sides
    sides === nothing || return sides
    c.sides = Dict{Tuple{Int, Int, Bool}, Tuple{String, String}}(
        k => _forward_sides(s) for (k, s) in _kept_members(m))
end

"""
The sides that name constant `slot` of `step` in `state`: those of the kept member of
its kinetic class (`_class_sides`). An `:EqualAI` step's constant is one class in both
conformations; the inactive copy of an `:OnlyA` step is a class of its own, named after
its own step. A slot the step's group has no class for (the reverse rate of a
rapid-equilibrium step) takes the sides of its slot-1 class, so every parameter of every
step renders.
"""
function _constant_sides(step::Step, slot::Int, state::Symbol, m::_AnyMech)
    g = findfirst(group -> step in group, steps(m))
    g === nothing && error("Step not found in mechanism: $step")
    tag = m isa AllostericMechanism ? cat_allo_state(m, g) : :None
    state === :I && tag === :OnlyA && return _forward_sides(step)
    sides = _class_sides(m)
    inactive = state === :I && tag === :NonequalAI
    get(() -> sides[(g, 1, inactive)], sides, (g, slot, inactive))
end

# Every step constant is named after the reaction of its class's kept member: a rate
# constant (`k_`) after the direction it drives, an equilibrium constant (`K_`) after
# the direction whose products-over-reactants it is. A binding's K is named in the
# release direction, so it is a dissociation constant (`K_ES_to_E_S`); an
# isomerization's or other RE step's K in the stored direction (`K_ES_to_EP`).
name(p::Kfor, m::_AnyMech) =
    _render_reaction("k_", _constant_sides(p.step, 1, p.state, m), p.state)
name(p::Krev, m::_AnyMech) =
    _render_reaction("k_", reverse(_constant_sides(p.step, 2, p.state, m)), p.state)
function name(p::Krapid, m::_AnyMech)
    sides = _constant_sides(p.step, 1, p.state, m)
    _render_reaction("K_", is_binding(p.step) ? reverse(sides) : sides, p.state)
end
```

3. Delete `_enumerate_parameters_full` and its docstring. In the docstring of
   `_step_constants`, replace the lines

```text
`is_equilibrium(s)`. The walkers over group representatives
(`_enumerate_parameters_full`, `_cat_params`) apply it to each group's
representative step; `_step_parameters` applies it to every step, `_state_parts` to
every step of one conformation, with that conformation's tag, and
`_onlya_haldane_violation` to a group's first step to label the group.
```

   with

```text
`is_equilibrium(s)`. `_step_names` applies it to every step of one conformation, with
that conformation's tag, `_column_parameter` picks a joint-solve column's, and
`_onlya_haldane_violation` applies it to a group's first step to label the group.
```

4. In `_canonical_groups`, replace the lines from
   `    segment = _bottomless_re_segment(reaction, gs)` through the `error(…)` call that
   follows it with `    _assert_bottomed_re_segments(reaction, gs)`, and in its
   docstring replace ``(`_bottomless_re_segment`). The`` with
   ``(`_assert_bottomed_re_segments`). The``. After `_canonical_groups` add:

```julia

"""
Error when a rapid-equilibrium segment of `groups` has no bottom form
(`_bottomless_re_segment`): its rate is 0/0 where one side's metabolites are absent. The
mechanism constructors check their groups, and the inactive conformation's derivation
checks its pruned graph (`_state_rate_polys`).
"""
function _assert_bottomed_re_segments(reaction::EnzymeReaction,
                                      groups::Vector{Vector{Step}})
    segment = _bottomless_re_segment(reaction, groups)
    segment === nothing ||
        error("Mechanism: rapid-equilibrium segment {" *
              join(name.(segment), ", ") * "} has no form free of one side's " *
              "metabolites, so its rate is undefined when that side is absent " *
              "(0/0). Make one of the segment's binding steps steady-state.")
    nothing
end
```

5. In the docstring of `_is_chemistry`, replace the lines

```text
`_expand_change_allo_state`), `_onlya_haldane_violation`, the steady-state pivot
tie-break in `_assemble_constraints` and `_fused_substrate_binding` read it to tell
catalysis from binding."""
```

   with

```text
`_expand_change_allo_state`), `_onlya_haldane_violation`, the joint solve's priority
(`_priority_key`) and `_fused_substrate_binding` read it to tell catalysis from
binding."""
```

   and in the docstring of `_bottomless_re_segment` replace
   ``(the inactive conformation's pruned graph, `_state_allo_mechanism`).`` with
   ``(the inactive conformation's pruned graph, `_inactive_groups`).``

- [ ] **Step 6: Derive both conformations from the one graph in `src/rate_eq_derivation.jl`**

1. Delete `ANNOTATION_SUBSTITUTED` and its docstring (the first block after the
   `@nospecialize` comment), and `_build_wegscheider_rename_map` with its docstring.
2. Replace `_raw_symbolic_rate_polys(m::Mechanism) = …` (with its `_leaf_names` call)
   and `_leaf_names` with its docstring by:

```julia
_raw_symbolic_rate_polys(m::Mechanism) =
    _raw_symbolic_rate_polys(reaction(m), steps(m), _step_names(m, :None))
```

   and, after the `Type` method of `_raw_symbolic_rate_polys`, add:

```julia

"""The kinetic groups of `m`'s law in conformation `state`: every group for `:None` and
`:A`; for `:I`, the inactive conformation's pruned graph (`_inactive_groups`)."""
_state_groups(m::Mechanism, ::Symbol) = steps(m)
_state_groups(am::AllostericMechanism, state::Symbol) =
    state === :I ? _inactive_groups(am) : steps(am)

"""
The parameter Symbols of each flat step's constants in conformation `state` of `m` —
`[K]` at rapid equilibrium, `[kf, kr]` at steady state — through the naming chokepoint
(`name(p, m)`), over the steps of `_state_groups(m, state)`: `:None` for a plain
mechanism; `:A` or `:I` for an allosteric one, an `:EqualAI` group's constants taking its
one shared name in both.
"""
function _step_names(m::_AnyMech, state::Symbol)
    tag(g) = m isa AllostericMechanism && cat_allo_state(m, g) === :EqualAI ?
             :EqualAI : state
    [Symbol[name(p, m) for p in _step_constants(s, tag(g))]
     for (s, g) in _flat_steps(_state_groups(m, state))]
end
```

3. In `_num_den_exprs` for `EnzymeMechanism`, replace
   `    names = _leaf_names(m, _step_parameters(m), _build_wegscheider_rename_map(m))`
   with `    names = _step_names(m, :None)`.
4. Replace the docstring of `_dep_assignments` with:

```julia
"""
The assignments `sym = rhs` of the dependent parameters `dep`, in name order.
The joint solve is a full Gauss–Jordan elimination (`_reduce_joint_system`), so every
right-hand side reads only `Keq` and independent parameters, and the assignments may run
in any order.
"""
```

5. Replace `_equation_text` and its docstring with:

```julia
"""
The `rate_equation_string` text of `M`: the `params` destructure of `param_syms`, the
`concs` destructure, one line `sym = rhs` per dependent parameter of `dep`, and the
`v = E_total * (num) / (den)` line from `_num_den_exprs`. A dependent whose right-hand
side mentions `Keq` goes under `# Haldane constraints:`, every other one under
`# Wegscheider constraints:`. Each section is sorted by name, which is load-bearing: the
eq_hash dedup of rate-equivalent mechanisms compares these strings, so their line order
must not depend on the order the solve emits its dependents in. A metabolite named after
a parameter raises (`_assert_no_name_clash`).
"""
function _equation_text(@nospecialize(M::Type{<:AbstractEnzymeMechanism}), param_syms,
                        dep)
    _assert_no_name_clash(M, param_syms, dep)
    weg, hal = String[], String[]
    for (sym, rhs) in sort!(collect(dep); by = first)
        push!(_mentions(rhs, :Keq) ? hal : weg, "$sym = $(_expr_to_string(rhs))")
    end
    num, den = _num_den_exprs(M)
    lines = ["(; $(join(param_syms, ", "))) = params",
             "(; $(join(metabolites(M()), ", "))) = concs"]
    isempty(weg) || push!(lines, "# Wegscheider constraints:", weg...)
    isempty(hal) || push!(lines, "# Haldane constraints:", hal...)
    push!(lines, "v = E_total * ($(_expr_to_string(num))) / ($(_expr_to_string(den)))")
    join(lines, "\n")
end
```
6. Replace the two methods after `_equation_text` (the commented one for
   `EnzymeMechanism` and the one for `AllostericEnzymeMechanism`) with:

```julia
function rate_equation_string(@nospecialize(m::AbstractEnzymeMechanism))
    dep, indep = _dependent_param_exprs(typeof(m))
    _equation_text(typeof(m), (indep..., :Keq, :E_total), dep)
end
```

7. Replace `_ss_rate_constant_names` and its docstring with:

```julia
"""
The names of the steady-state rate constants of `em`'s law (`_step_names`): for an
`AllostericEnzymeMechanism`, those of both conformations, the inactive one's `I_` names
included. `_scale_rate_constants` scales exactly these, leaving the rapid-equilibrium
constants, `Keq`, `E_total`, `L` and the regulator constants alone.
"""
function _ss_rate_constant_names(@nospecialize(em::AbstractEnzymeMechanism))
    m = _concrete(em)
    states = m isa Mechanism ? (:None,) : (:A, :I)
    Set{Symbol}(n for st in states
                for ((s, _), ns) in zip(_flat_steps(_state_groups(m, st)),
                                        _step_names(m, st))
                if !is_equilibrium(s) for n in ns)
end
```

8. Replace the docstring of `_inactive_groups` with:

```julia
"""
The inactive conformation's step graph of `am`, aligned with `steps(am)`: an `:OnlyA`
group is empty, and every other group keeps the steps whose two forms are reachable
from the free enzyme over the non-`:OnlyA` groups (`_reachable_from_free`). The steps
are `am`'s own, in their stored direction, so the inactive law reads `am`'s orientation
and names. The derivation (`_state_rate_polys`), the joint solve (`_inactive_law`) and
the copy rule (`_redundant_copy_groups`) read the inactive state from it.
"""
```

9. Delete `_state_allo_mechanism`, `_state_parts` and `_state_mechanism` with their
   docstrings. Replace `_state_rate_polys` (its docstring and Task 4's body) with:

```julia
"""
    _state_rate_polys(am, state) -> (num, den, d_free)

The rate polynomials of `am`'s catalytic graph in conformation `state` (`:A` or `:I`),
in that state's names (`_step_names`): the active conformation on the full graph, the
inactive one on its pruned graph (`_inactive_groups`), each step in `am`'s stored
direction. King–Altman on the pruned graph derives the inactive conformation's binding
partition, and for a dead cycle its steady-state fluxes cancel, so `num` is 0. The
inactive graph is checked for a rapid-equilibrium segment with no bottom form
(`_assert_bottomed_re_segments`), as the mechanism constructors check the full graph.
"""
function _state_rate_polys(am::AllostericMechanism, state::Symbol)
    groups = _state_groups(am, state)
    state === :I && _assert_bottomed_re_segments(reaction(am), groups)
    _raw_symbolic_rate_polys(reaction(am), groups, _step_names(am, state))
end
```

10. Delete `_cat_params` and its docstring. Replace `_kreg_params` and its docstring
    with:

```julia
"""
The regulator constants of `am` listed with conformation `state`, site by site: for `:A`
every ligand the active conformation binds, an `:EqualAI` ligand under its one bare name
(`K_<ligand>reg`, state `:EqualAI`) and an `:OnlyA` or `:NonequalAI` one under its
`K_A_` name; for `:I` the `K_I_` constant of every `:OnlyI` or `:NonequalAI` ligand.
"""
_kreg_params(am::AllostericMechanism, state::Symbol) =
    state === :A ?
    Kreg[Kreg(site, lig, tag === :EqualAI ? :EqualAI : :A)
         for site in regulatory_sites(am)
         for (lig, tag) in zip(ligands(site), allo_states(site)) if tag !== :OnlyI] :
    Kreg[Kreg(site, lig, :I) for site in regulatory_sites(am)
         for (lig, tag) in zip(ligands(site), allo_states(site))
         if tag in (:OnlyI, :NonequalAI)]
```

11. Delete the section `# ─── Dependent parameter expressions ───…` up to
    `# ─── Rate body building helpers ───…`: the allosteric `_dependent_param_exprs`, the
    `Type` forwarder and the comment after it (`src/thermo.jl` holds both methods).
12. Replace `_reg_K` and its docstring with:

```julia
"""
The dissociation-constant name of ligand `lig` (allosteric tag `tag`) at regulatory
`site` of `am` in the active (`inactive = false`) or inactive conformation, or
`nothing` when the ligand is absent from that conformation (`:OnlyI` in the active
state, `:OnlyA` in the inactive one). An `:EqualAI` ligand has one bare name in both
conformations (`K_<ligand>reg`); an `:OnlyA`, `:OnlyI` or `:NonequalAI` ligand has its
conformation's `K_A_` or `K_I_` name. Renders the name via the `name(::Kreg, am)`
chokepoint.
"""
function _reg_K(am::AllostericMechanism, site, lig, tag::Symbol, inactive::Bool)
    tag === (inactive ? :OnlyA : :OnlyI) && return nothing
    name(Kreg(site, lig, tag === :EqualAI ? :EqualAI : inactive ? :I : :A), am)
end
```
13. In the docstring of `_mwc_state_polys`, replace the sentences from
    ``The I-state polys are always re-derived`` through ``so `num_I` is 0.`` with:

```text
The I-state polys come from the inactive conformation's pruned
graph (`_inactive_groups` drops `:OnlyA` groups and every form they disconnect from free
E). King–Altman on that subgraph derives the inactive conformation's binding partition,
and for a dead cycle the pruned graph's steady-state fluxes cancel exactly, so `num_I`
is 0.
```

14. In `_num_den_exprs` for `AllostericEnzymeMechanism`, replace

```julia
    # A-state catalytic param symbols (the tagged column set) drive `_poly_to_expr`'s
    # param/metabolite ordering split; the I-poly's `:I` symbols sort as non-params.
    cat_params = Set(_param_columns(_state_parts(am, :A)...))
```

    with

```julia
    # The active conformation's catalytic names drive `_poly_to_expr`'s
    # param/metabolite ordering split; the I-poly's `I_` names sort as non-params.
    cat_params = Set{Symbol}(n for ns in _step_names(am, :A) for n in ns)
```

- [ ] **Step 7: Enumeration and identify**

1. `src/mechanism_enumeration.jl`, docstring of `_expand_split_kinetic_group`: replace

```text
and by the combined state solve on the built child for an `AllostericMechanism`,
where a count the derivation refuses as too large gains nothing (`_split_gain_test`).
```

   with

```text
and by the joint solve on the built child for an `AllostericMechanism`
(`_split_gain_test`).
```

2. Replace the allosteric `_split_gain_test` and its docstring with:

```julia
"""
Gain test of the split move for an allosteric mechanism: `sel -> Bool`, true when the child
built from the selected units has more independent parameters than `am`, both counted by
the joint solve (`_independent_param_count`), which derives no rate polynomial.
"""
function _split_gain_test(am::AllostericMechanism, units)
    base = _independent_param_count(am)
    sel -> _independent_param_count(_apply_bipartitions(am, units[sel])) > base
end
```

3. In the docstring of `_redundant_copy_groups`, replace the lines

```text
is `:OnlyA`. The inactive state's graph is that of `_state_mechanism(am, :I)`: the steps
of the non-`:OnlyA` groups between forms reachable from the free enzyme
(`_inactive_groups`). The copy binds there only at the sites that graph keeps, and a
```

   with

```text
is `:OnlyA`. The inactive state's graph is `_inactive_groups(am)`: the steps of the
non-`:OnlyA` groups between forms reachable from the free enzyme. The copy binds there
only at the sites that graph keeps, and a
```

   (re-wrap the paragraph to 92 characters).
4. `src/identify_rate_equation.jl`: replace the docstring of `_rate_eq_dedup_key` and its
   filter with:

```julia
"""
Equation-identity key: the rendered rate-equation string without its `# …` header
lines. Two mechanisms with the same key compute the identical rate function. Used as a
CSV tag and the LOOCV distinct-equation key.

Exact-render identity only: it hashes the (provenance-stripped) equation string
verbatim, so it detects textually identical equations, not algebraic equivalence.
"""
function _rate_eq_dedup_key(eq_text::AbstractString)
    kept = Iterators.filter(split(eq_text, '\n')) do ln
        l = strip(ln)
        !startswith(l, "#")
    end
    hash(join(kept, '\n'))
end
```

- [ ] **Step 8: Adapt the test helpers that read removed internals**

1. `test/test_rate_eq_derivation.jl`:
   - In `_testhelper_positional_params`, delete `    fes = EnzymeRates._free_enz_set(mech)`
     and replace `        rep = EnzymeRates._group_rep(group, fes)` with
     `        rep = first(group)` (the chokepoint names every member of a class alike).
     In its regulator branch replace the `:EqualAI` case (from
     `                if tag === :EqualAI` through the `if haskey(nt, struct_key)` line)
     with

```julia
                if tag === :EqualAI
                    # An :EqualAI ligand has one bare constant in both conformations.
                    struct_key = EnzymeRates.name(
                        EnzymeRates.Kreg(site, lig, :EqualAI), mech)
                    if haskey(nt, struct_key)
```

   - In `_testhelper_raw_to_ode_params`, delete the three comment lines starting
     `    # Kinetic-group rename map:` and the line
     `    rename = EnzymeRates._build_wegscheider_rename_map(mech)`, and replace the
     comment line `    # Resolve a structural param key through the rename map if not present`
     and the three-line `_testhelper_lookup` definition under it with
     `    _testhelper_lookup(k) = Float64(raw_params[k])` (no constant is absorbed any
     more, so every key exists).
   - In `_testhelper_eval_rate_string`, replace
     `           replace(s, EnzymeRates.ANNOTATION_SUBSTITUTED => "") * "\nend"` with
     `           s * "\nend"`; in `_testhelper_undefined_rhs_symbols`, replace
     `        line = replace(raw, EnzymeRates.ANNOTATION_SUBSTITUTED => "")` with
     `        line = raw`.
2. `test/test_mechanism_enumeration.jl`:
   - In "Tagged groups exclude inactive-state params", replace
     `        EnzymeRates._step_constants(EnzymeRates._group_reps(am)[g], :I), Ref(am))`
     with
     `        EnzymeRates._step_constants(first(EnzymeRates.steps(am)[g]), :I), Ref(am))`.
   - In "_hyperbolic_catalysis matches the derived denominator", replace
     `        derived = Dict{Tuple{Symbol, EnzymeRates.Mechanism}, Bool}()` with
     `        derived = Dict{Tuple{Symbol, Vector{Vector{EnzymeRates.Step}}}, Bool}()` and
     `            derived, (state, EnzymeRates._state_mechanism(m, state))) do` with
     `            derived, (state, EnzymeRates._state_groups(m, state))) do`.
   - In "to_allosteric shares :EqualAI downstream constants across conformations",
     replace the first comment paragraph (from
     ``    # An inactive-conformation graph (`_state_mechanism(am, :I)`) drops each`` to
     ``    # solve fails to merge it and over-counts.``) with

```julia
    # The inactive conformation drops each `:OnlyA` binding step, so that ligand's
    # downstream complex is reached only by conformational flip or reverse catalysis.
    # A lumped `:EqualAI` group is one class in both conformations, named after the
    # kept member of the full mechanism, so the inactive conformation cannot give it a
    # second name and the joint solve counts it once.
```

3. `test/test_types.jl`:
   - Delete the testset "_solve_dependent_set rejects a contradictory system"; its
     fixture (rows `x = 0` and `x = log Keq`) and assertion test
     `_reduce_joint_system` in `test/test_joint_solve.jl` since Task 6.
   - In "the naming cache never changes a mechanism's identity", replace the two
     `_testhelper_param_names` methods with

```julia
    _testhelper_param_names(m::EnzymeRates.Mechanism) = EnzymeRates._step_names(m, :None)
    _testhelper_param_names(m::EnzymeRates.AllostericMechanism) =
        [[EnzymeRates._step_names(m, state) for state in (:A, :I)];
         [[EnzymeRates.name(p, m)] for state in (:A, :I)
          for p in EnzymeRates._kreg_params(m, state)]]
```

4. `test/test_pivot_priority_regression.jl`: replace the two ABOUTME lines with

```julia
# ABOUTME: Regression: a Wegscheider tie between a free-enzyme binding K and an occupied
# ABOUTME: one makes one of them dependent, so fitted_params counts the tie once.
```

   and the comment block in the testset (from
   `    # A non-allosteric LDH mechanism reached in enumeration.` through
   `    # -1 param-count edge, impossible on a correct kernel).`) with

```julia
    # A non-allosteric LDH mechanism reached in enumeration. A closed all-RE cycle
    # — E→ELactate→ELactateNAD and E→ENAD→ELactateNAD with shared K_ENAD_to_E_NAD —
    # forces the Wegscheider tie K_ELactate_to_E_Lactate =
    # K_ELactateNADH_to_ENADH_Lactate. The solve must make one of the two dependent,
    # or fitted_params over-counts by 1 (7 instead of 6). A split that recovers the
    # tie would then be a -1 parameter-count edge.
```

- [ ] **Step 9: Re-record the pins the switch moves**

Each edit below is backed by a value check that already runs: the spec's
`Mass-action oracle`, `Analytical Rate` and `Rate Equation String` testsets for the
fixtures file; the hand-built ground truth for the AGT gate; the peak values for the
fitting-scale tests; exact children for the enumeration tests.

1. Fixtures file:
   - In the struct, replace
     `    expected_n_mirror_constraints::Int = 0    # RHS is a single Symbol (allosteric :EqualAI rename)`
     and the two `# NOTE:` comment lines under it with
     `    expected_n_mirror_constraints::Int = 0    # RHS is a single Symbol: a tie K_a = K_b`.
   - Delete the line `        expected_n_mirror_constraints=1,` from the PFK-1, HK and
     m_all specs (B10: no mirror).
   - In "Non-competitive Inhibitor", "Non-essential Activator" and "Activator +
     Competitive Inhibitor", add after their `expected_n_haldane_constraints=…,` line
     `        expected_n_mirror_constraints=1,`, `…=2,` and `…=2,` respectively (the
     single-symbol ties, now dependents: `K_ERinhS_to_ERinh_S = K_ES_to_E_S`;
     `K_EPRinh_to_ERinh_P = K_EP_to_E_P` and `K_ERinhS_to_ERinh_S = K_ES_to_E_S`;
     `K_EAinhP_to_EAinh_P = K_EP_to_E_P` and `K_EAinhS_to_EAinh_S = K_ES_to_E_S`).
   - Replace the factored strings (the law reads each tied class by its own name):
     - Non-competitive Inhibitor, denominator second line:
       `        " + R * S / (K_ERinhS_to_ERinh_S * K_ERinh_to_E_Rinh)",`
     - Non-essential Activator:
       numerator `"k_ES_to_EP * S / K_ES_to_E_S + k_ERinhS_to_EPRinh * R * S / (K_ERinhS_to_ERinh_S * K_ERinh_to_E_Rinh) - (k_EP_to_ES * P / K_EP_to_E_P + k_EPRinh_to_ERinhS * P * R / (K_EP_to_E_P * K_ERinh_to_E_Rinh))"`;
       denominator `"1 + P / K_EP_to_E_P + R / K_ERinh_to_E_Rinh + S / K_ES_to_E_S + P * R / (K_EP_to_E_P * K_ERinh_to_E_Rinh) + R * S / (K_ERinhS_to_ERinh_S * K_ERinh_to_E_Rinh)"`
     - Activator + Competitive Inhibitor:
       numerator `"k_ES_to_EP * S / K_ES_to_E_S + k_EAinhS_to_EAinhP * A * S / (K_EAinhS_to_EAinh_S * K_EAinh_to_E_Ainh) - (k_EP_to_ES * P / K_EP_to_E_P + k_EAinhP_to_EAinhS * A * P / (K_EAinhP_to_EAinh_P * K_EAinh_to_E_Ainh))"`;
       denominator `"1 + A / K_EAinh_to_E_Ainh + I / K_EIinh_to_E_Iinh + P / K_EP_to_E_P + S / K_ES_to_E_S + A * P / (K_EAinhP_to_EAinh_P * K_EAinh_to_E_Ainh) + A * S / (K_EAinhS_to_EAinh_S * K_EAinh_to_E_Ainh)"`
2. `test/test_rate_eq_derivation.jl`:
   - Replace the testset "rate_equation_string marks only the ties folded into v as
     substituted" with:

```julia
@testset "rate_equation_string prints every tie as a constraint line" begin
    # A binds E and E(B) under one constant while B's two bindings keep their own, so
    # the RE binding square ties B's two dissociation constants. The occupied binding's
    # constant is the dependent one, and v reads it.
    folded = @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
            E + B ⇌ E(B)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    @test occursin("# Wegscheider constraints:\nK_EAB_to_EA_B = K_EB_to_E_B\n",
                   rate_equation_string(folded))
    # S binds E and Estar under one constant, so the square of the two isomerizations
    # ties their equilibrium constants; the isomerization between occupied forms is the
    # dependent one.
    iso_tie = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E + S ⇌ E(S), Estar + S ⇌ Estar(S))
            E ⇌ Estar
            E(S) ⇌ Estar(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    @test occursin("\nK_ES_to_EstarS = K_E_to_Estar\n", rate_equation_string(iso_tie))
    # I binding E and I binding Estar share one group and run opposite ways around the
    # cycle, so the Haldane relation fixes the Theorell–Chance constant at Keq; v reads it.
    keq_tie = @enzyme_mechanism begin
        substrates: S
        products: P
        regulators: I
        steps: begin
            E + S ⇌ Estar + P
            (E + I <--> E(I), Estar + I <--> E(I))
        end
    end
    @test occursin("\nK_E_S_to_Estar_P = Keq\n", rate_equation_string(keq_tie))
end
```

   - In "rate_equation_string prints residual names bare in constraint lines", replace
     the regex test with
     `    @test occursin(r"\nk_EQ_to_EB_res_\+A_-P = .* \* k_EA_to_EP_res_\+A_-P \* ", s)`
     (the Haldane dependent is now the uptake of Q into the covalent form, E7).
   - Replace the testset `"_cat_params and _kreg_params"` (inside "Parameter-struct
     allosteric helpers", keeping `cm_src`) with:

```julia
    @testset "_step_names and _kreg_params" begin
        names_of(am, params) = [EnzymeRates.name(p, am) for p in params]
        flat_names(am, state) = reduce(vcat, EnzymeRates._step_names(am, state))
        # :NonequalAI cat groups + :NonequalAI reg ligand → both conformations named.
        aem = _testhelper_allo_from_source(
            cm_src, (2, (:NonequalAI, :EqualAI, :NonequalAI)),
            (((:R,), 1, (:NonequalAI,)),),
        )
        am = EnzymeRates.AllostericMechanism(aem)
        # Catalytic side: every step of the inactive law takes its inactive name, the
        # :EqualAI iso its one shared name.
        @test Set(flat_names(am, :I)) ==
              Set([:K_I_ES_to_E_S, :k_ES_to_EP, :k_EP_to_ES, :K_I_EP_to_E_P])
        @test Set(flat_names(am, :A)) ==
              Set([:K_A_ES_to_E_S, :k_ES_to_EP, :k_EP_to_ES, :K_A_EP_to_E_P])
        @test names_of(am, EnzymeRates._kreg_params(am, :A)) == [:K_A_Rreg]
        @test names_of(am, EnzymeRates._kreg_params(am, :I)) == [:K_I_Rreg]

        # Balanced :OnlyA bindings (both S and P) keep the :NonequalAI iso
        # Haldane-valid. The inactive conformation keeps no form but E, so its law has
        # no step and the :OnlyA reg ligand no inactive constant.
        aem_skip = _testhelper_allo_from_source(
            cm_src, (2, (:OnlyA, :NonequalAI, :OnlyA)),
            (((:R,), 1, (:OnlyA,)),),
        )
        am_skip = EnzymeRates.AllostericMechanism(aem_skip)
        @test isempty(EnzymeRates._step_names(am_skip, :I))
        @test isempty(EnzymeRates._kreg_params(am_skip, :I))
        @test names_of(am_skip, EnzymeRates._kreg_params(am_skip, :A)) == [:K_A_Rreg]

        # An :OnlyI reg ligand has an inactive constant and no active one; an :EqualAI
        # ligand one bare constant, listed with the active ones.
        am_onlyi = EnzymeRates.AllostericMechanism(_testhelper_allo_from_source(
            cm_src, (2, (:NonequalAI, :EqualAI, :NonequalAI)),
            (((:R,), 1, (:OnlyI,)),),
        ))
        @test isempty(EnzymeRates._kreg_params(am_onlyi, :A))
        @test names_of(am_onlyi, EnzymeRates._kreg_params(am_onlyi, :I)) == [:K_I_Rreg]
        am_equal = EnzymeRates.AllostericMechanism(_testhelper_allo_from_source(
            cm_src, (2, (:NonequalAI, :EqualAI, :NonequalAI)),
            (((:R,), 1, (:EqualAI,)),),
        ))
        @test names_of(am_equal, EnzymeRates._kreg_params(am_equal, :A)) == [:K_Rreg]
        @test isempty(EnzymeRates._kreg_params(am_equal, :I))
    end
```

   - Replace the testset "allosteric single-symbol Wegscheider tie is absorbed in both
     states" (keep its two fixtures, `tie` and `grouped`, verbatim) with a testset named
     `"an allosteric Wegscheider tie leaves the grouped mechanism's model"`, whose
     comment reads

```julia
    # A binds E and E(B) under one shared K while B's two bindings keep their own, so
    # the RE binding square ties B's two K's (K_EAB_to_EA_B = K_EB_to_E_B) in each
    # conformation. The tied constants are dependents, so the mechanism fits the same
    # constants as the one that groups B's two bindings outright, and its law takes the
    # same values.
```

     and whose assertions after the two fixtures read

```julia
    @test EnzymeRates.fitted_params(tie) == EnzymeRates.fitted_params(grouped) ==
        (:K_A_EA_to_E_A, :K_A_EB_to_E_B, :K_A_EP_to_E_P, :K_I_EA_to_E_A,
         :K_I_EB_to_E_B, :K_I_EP_to_E_P, :k_A_EAB_to_EP, :k_I_EAB_to_EP, :L)
    rng = Random.MersenneTwister(4)
    for _ in 1:5
        params = _testhelper_random_reduced_params(tie; rng)
        concs = (A = 0.5 + rand(rng), B = 0.5 + rand(rng), P = 0.5 + rand(rng))
        @test rate_equation(tie, concs, params) ≈ rate_equation(grouped, concs, params) rtol = 1e-12
    end
```

   - In the testset with the comment "This degenerate over-parametrized Wegscheider
     cycle has an ambiguous independent-parameter basis", replace those five comment
     lines with

```julia
    # Which koff is Wegscheider-dependent follows the solve's priority. The values
    # come from the mechanism's own fitted_params, so the Haldane property (zero net
    # rate at equilibrium) is verified whatever the priority picks.
```

3. `test/test_allosteric_collapse.jl`, testset "binding-Wegscheider single inner edge ->
   full collapse": replace

```julia
        s=replace(ER.rate_equation_string(am)," "=>"")
        @test occursin("K_I_EAB_to_EB_A=K_A_EAB_to_EB_A", s)  # explicit mirror
```

   with

```julia
        dep,_=ER._dependent_param_exprs(am)
        @test dep[:K_I_EAB_to_EB_A] == dep[:K_A_EAB_to_EB_A]  # one value in both states
```

   (both copies of the forbidden split are now dependents with the same right side,
   `K_EAB_to_EA_B * K_EA_to_E_A * (1 / K_EB_to_E_B)`; the testset's `abs(veq) < 1e-8`
   check already holds the value).
4. `test/allosteric_ground_truth.jl`, gate ":NonequalAI ping-pong MWC derivation matches
   mass-action ground truth": replace the `fp` tuple with

```julia
    @test fp == (Symbol("k_A_EA_to_E_res_+A_-P_P"), Symbol("k_A_EB_res_+A_-P_to_E_Q"),
                 Symbol("k_A_E_Q_to_EB_res_+A_-P"), :k_EA_to_E_A,
                 Symbol("k_EB_res_+A_-P_to_E_res_+A_-P_B"), :k_E_A_to_EA,
                 Symbol("k_E_res_+A_-P_B_to_EB_res_+A_-P"),
                 Symbol("k_I_EA_to_E_res_+A_-P_P"), Symbol("k_I_EB_res_+A_-P_to_E_Q"),
                 Symbol("k_I_E_Q_to_EB_res_+A_-P"), :L)
```

   replace the four comment lines from
   ``        # `k_EA_to_E_A` and `k_I_E_Q_to_EB_res_+A_-P` are absent from`` through
   ``        # k2r_I is computed here the same way, and the derivation rederives k4r_I.``
   with

```julia
        # `k_A_E_res_+A_-P_P_to_EA` and `k_I_E_res_+A_-P_P_to_EA`, the P rebinding
        # rates, are absent from fitted_params: each conformation's Haldane makes the
        # constant that runs against the reaction dependent. The oracle takes k2r_A and
        # k4r_I and derives k1r and k2r_I by the Haldane relations, so both are
        # computed here the same way, and the derivation rederives k2r_A and k2r_I.
```

   and add `                 :k_EA_to_E_A => k1r,` as the second entry of `d` (after
   `:k_E_A_to_EA => k1f,`). The gate's 1e-10 comparison with the hand-built network is
   the value check.
5. `test/test_fitting_scale.jl`:
   - "a law that never saturates peaks at the upper edge and is flagged": the fused
     release's dissociation constant is now the Haldane dependent (other chemistry is
     eliminated before a plain binding). Append to the comment
     `        # release's dissociation constant K_ES_to_E_P = Keq·k_ES_to_E_S/k_E_S_to_ES is`
     `        # 1.5 here.` (after its `at zero product v = k_E_S_to_ES·S.` line, starting
     that sentence with `The`), and replace
     `        p = (K_ES_to_E_P = 1.5, k_E_S_to_ES = 2.0, Keq = 2.0)` with
     `        p = (k_ES_to_E_S = 1.5, k_E_S_to_ES = 2.0, Keq = 2.0)` — the same constants.
   - "a product can raise the peak above the product-free turnover": set the slow step
     by role, so the pins `0.249991171308`, `125.62056353` and `14120.94` stay:

```julia
        pn = EnzymeRates.fitted_params(rob)
        rob_mech = EnzymeRates.Mechanism(rob)
        release_p = only(s for (s, _) in EnzymeRates._flat_steps(rob_mech)
                         if EnzymeRates.consumed(s) == [EnzymeRates.Product(:P)] &&
                            isempty(EnzymeRates.bound(EnzymeRates.from_species(s))))
        slow = Set(EnzymeRates.name(c, rob_mech)
                   for c in EnzymeRates._step_constants(release_p, :None))
        p = merge(NamedTuple{pn}(Tuple(k in slow ? 1e-3 : 1.0 for k in pn)),
                  (Keq = 1.0,))
```

     (replacing the two lines that build `p` from `k === :k_E_P_to_EP`). With every
     other constant 1 and `Keq = 1`, each dependent comes out 1 or 1e-3, so the step
     constants are those the comment describes.
6. `test/test_mechanism_enumeration.jl`:
   - Testset "AllostericMechanism — a split whose count needs a law over the term cap":
     rename it `"AllostericMechanism — the split count of a law over the term cap"`,
     replace its comment with

```julia
    # Merged random-order ter-bi with the abortive complexes E(B, P) and E(A, B, P). Every
    # steady-state binding sits alone in its group except the A bindings at E(B) and
    # E(B, P), and the active state has more denominator terms than
    # MAX_RATE_EQUATION_TERMS. Splitting the P group by A leaves E(A, B) + P alone, and the
    # square E(B), E(A, B), E(A, B, P), E(B, P), closed by the shared A group, ties its K
    # to the other part's: the joint solve counts that child, which derives no rate
    # polynomial, at the parent's count, so the split gains nothing and is not emitted.
    # The move emits the two splits that gain, the P group by B and the A group by P, and
    # expand_mechanisms keeps them beside the other moves' children.
```

     and replace
     `    @test_throws "polynomial terms" EnzymeRates._independent_param_count(p_by_a)`
     with

```julia
    @test EnzymeRates._independent_param_count(p_by_a) ==
          EnzymeRates._independent_param_count(am)
```

     (its children assertions stay: still exactly `p_by_b` and `a_by_p`).
   - Replace the testset "AllostericMechanism — a parent whose count needs a law over
     the term cap" (keep the closing `end` of the enclosing testset after it) with:

```julia
@testset "AllostericMechanism — a parent over the term cap splits by its count" begin
    # Merged random-order ter-bi over MAX_RATE_EQUATION_TERMS whose P group is split by A:
    # E(A, B) + P sits alone, its K tied to the other part's by the square E(B), E(A, B),
    # E(A, B, P), E(B, P) that the shared A group closes. The joint solve counts the parent
    # without deriving its law, so the move weighs every split against that count: the A
    # group split by P and the remaining P group split by B gain, and expand_mechanisms
    # keeps them beside the other moves' children.
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q
        catalytic_steps: begin
            E + A <--> E(A)                                         :: OnlyA
            (E(B) + A <--> E(A, B), E(B, P) + A <--> E(A, B, P))    :: OnlyA
            E(D) + A <--> E(A, D)                                   :: OnlyA
            E(B, D) + A <--> E(P, Q)                                :: OnlyA
            E + B <--> E(B)                                         :: OnlyA
            E(A) + B <--> E(A, B)                                   :: OnlyA
            E(D) + B <--> E(B, D)                                   :: OnlyA
            E(A, D) + B <--> E(P, Q)                                :: OnlyA
            E + D <--> E(D)                                         :: OnlyA
            E(A) + D <--> E(A, D)                                   :: OnlyA
            E(B) + D <--> E(B, D)                                   :: OnlyA
            E(A, B) + D <--> E(P, Q)                                :: OnlyA
            E(A, B) + P ⇌ E(A, B, P)                                :: OnlyA
            (E(B) + P ⇌ E(B, P), E(Q) + P ⇌ E(P, Q))                :: OnlyA
            E + Q <--> E(Q)                                         :: EqualAI
        end
    end)
    a_by_p = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q
        catalytic_steps: begin
            E + A <--> E(A)                                         :: OnlyA
            E(B) + A <--> E(A, B)                                   :: OnlyA
            E(B, P) + A <--> E(A, B, P)                             :: OnlyA
            E(D) + A <--> E(A, D)                                   :: OnlyA
            E(B, D) + A <--> E(P, Q)                                :: OnlyA
            E + B <--> E(B)                                         :: OnlyA
            E(A) + B <--> E(A, B)                                   :: OnlyA
            E(D) + B <--> E(B, D)                                   :: OnlyA
            E(A, D) + B <--> E(P, Q)                                :: OnlyA
            E + D <--> E(D)                                         :: OnlyA
            E(A) + D <--> E(A, D)                                   :: OnlyA
            E(B) + D <--> E(B, D)                                   :: OnlyA
            E(A, B) + D <--> E(P, Q)                                :: OnlyA
            E(A, B) + P ⇌ E(A, B, P)                                :: OnlyA
            (E(B) + P ⇌ E(B, P), E(Q) + P ⇌ E(P, Q))                :: OnlyA
            E + Q <--> E(Q)                                         :: EqualAI
        end
    end)
    p_by_b = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q
        catalytic_steps: begin
            E + A <--> E(A)                                         :: OnlyA
            (E(B) + A <--> E(A, B), E(B, P) + A <--> E(A, B, P))    :: OnlyA
            E(D) + A <--> E(A, D)                                   :: OnlyA
            E(B, D) + A <--> E(P, Q)                                :: OnlyA
            E + B <--> E(B)                                         :: OnlyA
            E(A) + B <--> E(A, B)                                   :: OnlyA
            E(D) + B <--> E(B, D)                                   :: OnlyA
            E(A, D) + B <--> E(P, Q)                                :: OnlyA
            E + D <--> E(D)                                         :: OnlyA
            E(A) + D <--> E(A, D)                                   :: OnlyA
            E(B) + D <--> E(B, D)                                   :: OnlyA
            E(A, B) + D <--> E(P, Q)                                :: OnlyA
            E(A, B) + P ⇌ E(A, B, P)                                :: OnlyA
            E(B) + P ⇌ E(B, P)                                      :: OnlyA
            E(Q) + P ⇌ E(P, Q)                                      :: OnlyA
            E + Q <--> E(Q)                                         :: EqualAI
        end
    end)
    @test EnzymeRates._eq_complexity(am) > EnzymeRates.MAX_RATE_EQUATION_TERMS
    @test EnzymeRates._independent_param_count(am) == 22
    @test EnzymeRates._independent_param_count(a_by_p) == 24
    @test EnzymeRates._independent_param_count(p_by_b) == 23
    kids = EnzymeRates._expand_split_kinetic_group(am)
    @test length(kids) == 2
    @test Set(kids) == Set([a_by_p, p_by_b])
    others = [EnzymeRates._expand_re_to_ss(am); EnzymeRates._expand_change_allo_state(am)]
    @test Set(EnzymeRates.expand_mechanisms([am], EnzymeRates.reaction(am))) ==
          Set([others; a_by_p; p_by_b])
end
```

7. `test/test_identify_rate_equation.jl`:
   - In testset "_rate_eq_dedup_key", replace the comment and the four lines that build
     and compare `a` and `b` with

```julia
        # differs only in a comment header line:
        a = replace(base, "# Haldane constraints:" => "# Haldane constraint:")
        @test EnzymeRates._rate_eq_dedup_key(a) ==
              EnzymeRates._rate_eq_dedup_key(base)
        # differs in a Wegscheider tie -> different key:
        b = replace(base, "# Haldane constraints:" =>
                          "# Wegscheider constraints:\nK_x = K_a\n# Haldane constraints:")
        @test EnzymeRates._rate_eq_dedup_key(base) !=
              EnzymeRates._rate_eq_dedup_key(b)
```

   - Replace `_testhelper_dedup_twins` and its comment (the regrouping pair no longer
     prints one law) with:

```julia
# Two distinct bi-bi mechanisms that render one rate equation: the second declares a
# competitive inhibitor that no step binds, which the compiled mechanism leaves out.
# Written out directly so the tests are self-contained and cheap. The tests re-verify the
# collision at run time.
function _testhelper_dedup_twins()
    m1 = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products:   P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
        end
    end)
    rxn = EnzymeRates.reaction(m1)
    with_inhibitor = EnzymeRates.EnzymeReaction(EnzymeRates.reactants(rxn),
        [EnzymeRates.RegulatorMults(EnzymeRates.CompetitiveInhibitor(:I), [1])],
        EnzymeRates.allowed_catalytic_multiplicities(rxn))
    m2 = EnzymeRates.Mechanism(with_inhibitor, EnzymeRates.steps(m1))
    return m1, m2
end
```

   - Rename the testset "allosteric split/merge pair: the split form carries one more
     param" to `"allosteric split/merge pair: the split form is the merged form's model"`
     and replace its last comment and two assertions with

```julia
    # The split constant is tied back: the square E, E(Lactate), E(Lactate, NADH),
    # E(NADH), closed by the shared NADH group, makes K_ELactateNADH_to_ENADH_Lactate equal
    # the free-enzyme Lactate constant in the active conformation, and the inactive one
    # binds no NADH. Both forms fit the same constants.
    @test EnzymeRates._independent_param_count(am1) ==
          EnzymeRates._independent_param_count(am2) == 7
    @test EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(am1)) ==
          EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(am2))
```

     Today the split form counts 8: the active conformation folds the free-enzyme
     Lactate constant into the split one while the inactive law still reads it (known
     issue 2), which the old finite-difference rank of 9 against 8 measured as a real
     parameter. The joint solve makes it the merged form's model.
8. Regenerate the golden reference after every other edit (five lines change: PFK-1,
   HK and m_all print `K_Pireg`/`K_R2reg` bare with no mirror line; the two LDH 6-group
   specs name a catalytic class after its kept member):

```bash
$JL -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); src = read("test/test_allosteric_golden.jl", String); a = findfirst("\"\"\"Canonical serialization", src); b = findfirst("@testset", src); include_string(Main, src[first(a):first(b)-1]); path = "test/reference/allosteric_golden_reference.txt"; write(path, join(_testhelper_allosteric_golden_lines(), "\n") * "\n"); println("golden lines ", length(readlines(path)))'
git diff --stat test/reference/allosteric_golden_reference.txt
```

   Expected: `golden lines 22`, `5 insertions(+), 5 deletions(-)`. Read each changed
   line: the names may change, the structure may not (each spec's `Rate Equation String`
   test evaluates the string against `rate_equation`).

- [ ] **Step 10: Run the focused tests**

Run the focused commands for `test_oracles`, `test_joint_solve`, `test_types`,
`test_rate_eq_derivation`, `test_allosteric_collapse`, `test_allosteric_golden`,
`allosteric_ground_truth`, `test_fitting_scale`, `test_identify_rate_equation` and
`test_mechanism_enumeration`. Expected: all pass; `Laws over a regenerated population`
passes all its tests (about 30 s). A failure in a pin not listed above means the code
differs from the prototype: find the cause with superpowers:systematic-debugging before
touching the pin.

- [ ] **Step 11: Confirm nothing reads a removed name**

Run `grep -rn "ANNOTATION_SUBSTITUTED\|_build_wegscheider_rename_map\|_leaf_names\|_step_parameters\|_free_enz_set\|_step_priority\|_group_rep\|_rep_sides\|_enumerate_parameters_full\|_param_columns\|_assemble_constraints\|_solve_dependent_set\|_state_allo_mechanism\|_state_parts\|_state_mechanism\|_cat_params" src test`.
Expected: no output.

- [ ] **Step 12: Run the full suite, then commit**

Run the full suite. Expected: all pass, about 33,300 tests in about 7 minutes; the three
`rate_equation` info lines stay near 77, 117 and 103 ns. If JET reports a method error
on a `Union` argument, give the function one method per type (as `_state_groups` has).

```bash
git status
git add src test
git commit -m "Derive every law from one joint thermodynamic solve" \
  -m "One solve over the constants of both conformations picks the dependents by the fixed
priority, names each kinetic class after the member it keeps, and drops the second rename
pass. The inactive conformation's law is derived on the one graph with the mechanism's
own orientation and names, so an :EqualAI constant can no longer take a second fitted
name or fold one way in one conformation and the other way in the other (ROADMAP known
issues 1 and 2). :EqualAI regulator ligands print bare. Every law still equals the exact
mass-action oracle; the re-recorded pins change names and the choice of dependents only." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---
### Task 8: Fold the guard's cycle helpers into `src/thermo.jl`

**Files:**
- Modify: `src/thermo.jl` (append what remains of the old file)
- Delete: `src/thermodynamic_constr_for_rate_eq_derivation.jl`
- Modify: `src/EnzymeRates.jl`, `ROADMAP.md`

**Interfaces:**
- Consumes: Task 7's state of both files.
- Produces: one file for the thermodynamics: `_rref_partition`, `_rational_nullspace`,
  `_integer_nullspace`, `_thermodynamic_constraints`, `_count_kind`, `_add_step_column!`,
  `_partition_independent_count`, `_has_strict_positive_combination`,
  `_onlya_violation_message` and `_onlya_haldane_violation` move unchanged.

A pure move (spec section 3 proposes `thermo.jl` for the solve): no function changes.
The `:OnlyA` guard keeps its `Rational{BigInt}` null space until phase 4 replaces it with
admission.

- [ ] **Step 1: Move the code**

1. Append to `src/thermo.jl` the line
   `# ─── Cycle-space helpers of the :OnlyA guard and the split move's counter ─────`
   followed by a blank line and every line of
   `src/thermodynamic_constr_for_rate_eq_derivation.jl` after its two ABOUTME lines,
   unchanged.
2. Replace the two ABOUTME lines of `src/thermo.jl` with

```julia
# ABOUTME: The joint thermodynamic solve over the constants of both conformations, and the
# ABOUTME: cycle-space helpers of the :OnlyA guard and of the split move's counter.
```

3. `git rm src/thermodynamic_constr_for_rate_eq_derivation.jl`, and in
   `src/EnzymeRates.jl` delete the line
   `include("thermodynamic_constr_for_rate_eq_derivation.jl")`.
4. `ROADMAP.md`: in the entries that cite it, replace
   `src/thermodynamic_constr_for_rate_eq_derivation.jl` with `src/thermo.jl`.

- [ ] **Step 2: Verify**

Run `grep -rn "thermodynamic_constr_for_rate_eq_derivation" src test docs/src ROADMAP.md README.md .claude`.
Expected: no output. Run the focused commands for `test_types` (the guard and the
chokepoint AST walk, which reads every `src/*.jl`) and `test_joint_solve`. Expected:
all pass.

- [ ] **Step 3: Run the full suite, then commit**

Run the full suite. Expected: all pass.

```bash
git status
git add src/thermo.jl src/EnzymeRates.jl ROADMAP.md
git commit -m "Keep the thermodynamics in one file" \
  -m "The cycle-space helpers of the :OnlyA guard and of the split move's counter move
unchanged into src/thermo.jl beside the joint solve." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

(`git rm` already staged the deletion.)

---

### Task 9: Rewrite the textbook bridge by role

**Files:**
- Modify: `test/test_rate_eq_derivation.jl` (`_testhelper_positional_flat_idx`,
  `_testhelper_positional_params`, `_testhelper_analytical_oracle_params` replaced;
  four call sites)
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (two comments)

**Interfaces:**
- Consumes: `all_parameters` (Task 3), the E3 chokepoint (Task 7),
  `_testhelper_positional_site_idx` (kept), `EnzymeRates._reg_K`.
- Produces: `_testhelper_stored_params(m, all_params::NamedTuple) -> NamedTuple`
  (`K{i}`, `k{i}f`, `k{i}r` in stored flat order and direction, for the QSSA and ODE
  cross-checks) and `_testhelper_written_params(m, all_params::NamedTuple;
  source_steps, source_reg_sites = nothing) -> NamedTuple` (the textbook names).

Spec section 14 (C5): each textbook test reads its written steps and takes the package
constant of the matching stored step in the written direction, through
`all_parameters`. The constructor stores every binding with its metabolite taken up, so
a written binding's direction comes from its metabolite's role: a product's binding is
written as its release (every spec writes product steps that way). Every other written
step keeps its direction in `source_steps`; on the 36 specs that carry them, each
matches its stored step in the stored direction, but the bridge inverts an RE
equilibrium constant and swaps a rate pair when one does not. Measured on the Task 7
state: all 34 analytical rate laws and 11 analytical kcats pass. The by-role bridge
reads the constants through the chokepoint, so it needs no group representative and no
rename; it follows the switch because before it, three plain specs had a rename-absorbed
binding constant that `all_parameters` did not list.

- [ ] **Step 1: Replace the positional helpers**

In `test/test_rate_eq_derivation.jl`:

1. Delete `_testhelper_positional_flat_idx` with its docstring (the block starting
   `Per-step flat index for positional oracle naming.`). Keep
   `_testhelper_positional_site_idx`.
2. Replace `_testhelper_positional_params` and `_testhelper_analytical_oracle_params`,
   with their docstrings (from the docstring starting `Re-key a structural-named
   parameter NamedTuple to the **per-step** positional` through the `end` of
   `_testhelper_analytical_oracle_params`), with:

```julia
"""
The constants of the plain mechanism `m` (values by name in `all_params`, from
`all_parameters`) under the positional names the QSSA and ODE cross-checks read: stored
flat step `i` gives `K{i}` at rapid equilibrium and `k{i}f`, `k{i}r` at steady state, in
its stored direction. Every entry of `all_params` passes through.
"""
function _testhelper_stored_params(m, all_params::NamedTuple)
    mech = _testhelper_lift(m)
    out = Pair{Symbol, Any}[]
    for (i, (s, _)) in enumerate(EnzymeRates._flat_steps(EnzymeRates.steps(mech)))
        if EnzymeRates.is_equilibrium(s)
            push!(out, Symbol("K", i) =>
                all_params[EnzymeRates.name(EnzymeRates.Krapid(s, :None), mech)])
        else
            push!(out,
                Symbol("k", i, "f") =>
                    all_params[EnzymeRates.name(EnzymeRates.Kfor(s, :None), mech)],
                Symbol("k", i, "r") =>
                    all_params[EnzymeRates.name(EnzymeRates.Krev(s, :None), mech)])
        end
    end
    merge((; out...), all_params)
end

"""
The constants of `m` (values by name in `all_params`, from `all_parameters`) under the
positional names a textbook rate law reads, by role. Written step `i` (in the order of
the flattened `source_steps`) takes the constants of the stored step it matches, in the
written direction: `k{i}f` and `k{i}r` at steady state, and at rapid equilibrium `K{i}`,
a binding's dissociation constant or, for any other step, the equilibrium constant of the
written direction. The inactive conformation of a `:NonequalAI` group adds `k{i}f_T`,
`k{i}r_T` or `K{i}_T`. The constructor stores every binding with its metabolite taken up,
so a written binding's direction is read from its metabolite's role: a product's binding
is written as its release. A ligand of regulatory site `j` (the order of
`source_reg_sites`) gives `K_<ligand>_reg{j}` and, when its inactive conformation has a
constant of its own, `K_<ligand>_T_reg{j}`. Every entry of `all_params` passes through.
"""
function _testhelper_written_params(m, all_params::NamedTuple; source_steps,
                                    source_reg_sites = nothing)
    mech = _testhelper_lift(m)
    allo = mech isa EnzymeRates.AllostericMechanism
    stored = EnzymeRates._flat_steps(EnzymeRates.steps(mech))
    out = Pair{Symbol, Any}[]
    for (i, w) in enumerate(Iterators.flatten(source_steps))
        j = findfirst(x -> first(x) == w || first(x) == EnzymeRates._reversed(w), stored)
        j === nothing && error("written step $i matches no stored step")
        s, g = stored[j]
        forward = EnzymeRates.is_binding(w) ?
                  !(EnzymeRates.bound_metabolite(w) isa EnzymeRates.Product) : s == w
        tag = allo ? EnzymeRates.cat_allo_state(mech, g) : :None
        states = tag === :NonequalAI ? ((:A, ""), (:I, "_T")) :
                 ((tag === :OnlyA ? :A : tag, ""),)
        for (state, suffix) in states
            if EnzymeRates.is_equilibrium(s)
                key = EnzymeRates.name(EnzymeRates.Krapid(s, state), mech)
                state === :I && !haskey(all_params, key) && continue
                K = all_params[key]
                push!(out, Symbol("K", i, suffix) =>
                    EnzymeRates.is_binding(s) || forward ? K : 1 / K)
            else
                kf_key = EnzymeRates.name(EnzymeRates.Kfor(s, state), mech)
                state === :I && !haskey(all_params, kf_key) && continue
                kf = all_params[kf_key]
                kr = all_params[EnzymeRates.name(EnzymeRates.Krev(s, state), mech)]
                f, r = forward ? (kf, kr) : (kr, kf)
                push!(out, Symbol("k", i, "f", suffix) => f, Symbol("k", i, "r", suffix) => r)
            end
        end
    end
    if allo
        for (pos, site) in enumerate(EnzymeRates.regulatory_sites(mech))
            j = _testhelper_positional_site_idx(site, pos, source_reg_sites)
            for (lig, tag) in zip(EnzymeRates.ligands(site), EnzymeRates.allo_states(site))
                for (inactive, suffix) in ((false, ""), (true, "_T"))
                    tag === :EqualAI && inactive && continue
                    K = EnzymeRates._reg_K(mech, site, lig, tag, inactive)
                    K === nothing ||
                        push!(out, Symbol("K_", EnzymeRates.name(lig), suffix, "_reg", j) =>
                                   all_params[K])
                end
            end
        end
    end
    merge((; out...), all_params)
end
```

3. Call sites:
   - `_testhelper_test_reference_qssa` and `_testhelper_test_ode_steadystate`: replace
     `_testhelper_positional_params(m, all_params)` with
     `_testhelper_stored_params(m, all_params)`.
   - `_testhelper_test_analytical_rate` and `_testhelper_test_analytical_kcat`: replace
     `_testhelper_analytical_oracle_params(` with `_testhelper_written_params(` (the
     keyword arguments stay).
   - Update the comment in `_testhelper_test_reference_qssa`,
     `# _testhelper_reference_qssa uses positional step-indexed k$(i)f/k$(i)r keys`, to
     `# _testhelper_reference_qssa reads stored-step k$(i)f/k$(i)r keys`; and in
     `_testhelper_test_analytical_kcat` replace the four comment lines from
     `        # The oracle's positional formula may reference a forward rate that is` to
     `        # included) to positional names.` with
     `        # The textbook formula may read a dependent rate constant, so the bridge`
     `        # reads every constant through all_parameters.`

In the fixtures file:

4. In the struct `_testhelper_MechanismTestSpec`, replace the comment line
   ``    # k1,k2,… (and reg1,reg2,…) by source order; the bridge in`` and the next line
   ``    # `_testhelper_analytical_oracle_params` uses these to remap onto canonical stored``
   with ``    # k1,k2,… (and reg1,reg2,…) by source order; `_testhelper_written_params` reads``
   and ``    # them to take each written step's constants from its stored step, by role``,
   keeping the line that follows.
5. In the MWC tetramer spec, replace the comment lines
   `    # Param naming uses kinetic-group representative-step indices:` through
   `    # k13f/k13r=catalysis SS (rep step 13).` with

```julia
    # Each written step i gives K{i} or k{i}f/k{i}r; the law reads the first step of
    # each group: K1 = S1 binding (step 1), K4 = P1 binding (step 4), K7 = S2 binding
    # (step 7), K10 = P2 binding (step 10), k13f/k13r = catalysis (step 13).
```

- [ ] **Step 2: Run the derivation tests**

Run the focused command for `test_rate_eq_derivation`.
Expected: all pass, the 34 `Analytical Rate` and 11 `Peak turnover approaches the
analytical kcat` testsets among them. Then
`grep -rn "_testhelper_positional_params\|_testhelper_analytical_oracle_params\|_testhelper_positional_flat_idx" test`
prints nothing.

- [ ] **Step 3: Run the full suite, then commit**

Run the full suite. Expected: all pass.

```bash
git status
git add test/test_rate_eq_derivation.jl test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "Map textbook constants to the package's by role" \
  -m "Each textbook test takes the constants of the stored step that matches a written
step, in the written direction, through all_parameters and the naming chokepoint. A
written binding's direction comes from its metabolite's role. The bridge needs no group
representative, so it survives any naming rule." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---

### Task 10: Document the joint solve; update CLAUDE.md and the ROADMAP

**Files:**
- Modify: `docs/src/deriving/thermodynamic_constraints.md`,
  `docs/src/deriving/kinetic_groups.md`, `docs/src/deriving/mwc_allostery.md`,
  `docs/src/developer.md`, `.claude/CLAUDE.md`, `ROADMAP.md`

**Interfaces:**
- Consumes: the behavior of Tasks 3–9. No code changes.

Spec section 16, phase-2 rows: the fixed priority and C8's note in
`thermodynamic_constraints.md`; the "Canonical Step Form" section of CLAUDE.md (only that
section changes in CLAUDE.md); the ROADMAP entries phase 2 retires or rewrites; every
page that states a changed name. Verified against the Task 7 state: the README
quickstart, `getting_started.md`, `textbooks.md`, `fitting/tutorial.md`,
`normalized_vs_absolute.md`, `identify/tutorial.md`, `iso.md`, `regulators.md`,
`re_vs_ss.md` and `cha_king_altman.md` name only constants that do not change; the
`@example` blocks regenerate their output.

- [ ] **Step 1: `docs/src/deriving/thermodynamic_constraints.md`**

1. Replace the section `## How the constraints are found` (its two paragraphs) with:

```markdown
## How the constraints are found

The package takes the fundamental cycles of a spanning forest of the enzyme-form graph,
one per step outside the forest, so the cycle rows come out **linearly independent by
construction**. Each cycle is classified by its net metabolite change: proportional to
the overall reaction makes it a Haldane cycle, zero makes it a Wegscheider cycle. An
allosteric mechanism writes every cycle once per conformation, and the steps of a kinetic
group add rows that tie their constants together. One exact rational elimination over
all of these rows chooses the dependent constants (integer arithmetic that switches to
big integers if a value grows too large), so the same constraints come out every time.
If a mechanism is thermodynamically infeasible — its rows combine to
`0 = c · log(Keq)` with `c ≠ 0` — the package raises an error rather than emitting a
silent or wrong equation.
```

2. Replace the section `## The choice of independent vs dependent parameters` (its
   numbered list and the paragraph after it) with:

```markdown
## The choice of independent vs dependent parameters

Each independent cycle costs one constant. The package decides which constant becomes
dependent by a fixed priority, read in turn; the earlier a constant ranks, the sooner it
becomes dependent:

1. **Conformation.** In an allosteric mechanism, the inactive conformation's own
   constants go before the active ones.
2. **Kind of step.** Isomerizations (steps that exchange no metabolite, on-enzyme
   chemistry included) go first, then other chemistry (a binding that also converts, a
   Theorell–Chance step), then plain bindings.
3. **Occupancy.** Steps whose two enzyme forms both hold a metabolite go before steps
   that touch free E or a form that holds only a covalent residual.
4. **Speed.** Steady-state constants go before rapid-equilibrium ones.
5. **Direction.** Within a steady-state step, the rate constant that runs against the
   reaction goes first: the reverse rate, except that for the binding of a product the
   uptake runs against the reaction. So in `E(S) <--> E + P`, stored as the binding
   `E + P → E(S)`, the product uptake is dependent and the catalytic rate stays fitted.

The free-enzyme binding constants a biochemist measures therefore stay independent
whenever a cycle leaves them free. Ties go to a content key built from the forms and
metabolites of each step, never from the order in which steps were written, so two
mechanisms that differ only in step order, step direction or how the kinetic groups were
labelled have the same dependents. A kinetic group's constant takes the name of the step
the solve keeps, so splitting a group whose parts the cycles already tie changes no
fitted name. The key reads metabolite names, though: renaming a metabolite can change
which of two tied constants is dependent.
```

- [ ] **Step 2: `docs/src/deriving/kinetic_groups.md`**

Replace the paragraph after the `ungrouped` example (from
`There are nine independent constants, with form-specific names:` through
`closes a thermodynamic cycle.`) with:

```markdown
There are nine independent constants, with form-specific names: `K_EA_to_E_A` for `A`
on free E, `K_EAB_to_EA_B` and `K_EAQ_to_EA_Q` for `B` and `Q` on `E(A)`, and similar
families for the others. Four more bindings are not fit at all but fixed by Wegscheider
relations (`K_EAB_to_EB_A`, `K_EAQ_to_EQ_A`, `K_EBP_to_EP_B`, `K_EPQ_to_EQ_P`), since
each catalytic and abortive loop closes a thermodynamic cycle.
```

- [ ] **Step 3: `docs/src/deriving/mwc_allostery.md`**

1. Replace `such equations, one per conformation, adding constraints neither has alone.`
   (and the sentence's start on the line before, `An MWC equation stacks two`) with
   `An MWC equation solves the relations of both conformations together, adding
   constraints neither has alone.` (re-wrap).
2. Replace the paragraph starting
   ``Now `K_I_EP_to_E_P = K_A_EP_to_E_P · K_I_ES_to_E_S / K_A_ES_to_E_S`, i.e.`` (two
   lines, through ``K_P^I / K_S^I = K_P^A / K_S^A`:``) with

```markdown
Now `K_I_ES_to_E_S = K_A_ES_to_E_S · K_I_EP_to_E_P / K_A_EP_to_E_P`, i.e.
`K_S^I / K_S^A = K_P^I / K_P^A`:
```

   keeping the sentence that follows.
3. In the paragraph that starts ``Often a K-system`` nothing changes. Confirm the
   statements ``The equation reports `K_I_ES_to_E_S = K_A_ES_to_E_S`.`` and
   ``(`k_I_ES_to_E_S = k_A_ES_to_E_S · k_I_E_S_to_ES / k_A_E_S_to_ES`)`` against the
   example output in the docs build (Step 7): both hold.

- [ ] **Step 4: `docs/src/developer.md`**

1. Replace ``free-enzyme forms and each group's naming representative and rendered sides), which``
   with ``rendered sides of each kinetic class's kept member), which``.
2. Replace the lines

```text
`_onlya_haldane_violation`, the steady-state pivot tie-break in
`_assemble_constraints` and `_fused_substrate_binding` read it to tell catalysis from
binding.
```

   with

```text
`_onlya_haldane_violation`, the joint solve's priority (`_priority_key`) and
`_fused_substrate_binding` read it to tell catalysis from binding.
```

3. Replace the sentences

```text
(`_partition_independent_count`). For an `AllostericMechanism` it is taken on the built
child and, when a Wegscheider tie folds a constant, reads both conformations' rate
polynomials, so a candidate whose law exceeds `MAX_RATE_EQUATION_TERMS` has no count and
gains nothing (`_split_gain_test`). The RE→SS move flips whole groups and accepts a
```

   with

```text
(`_partition_independent_count`). For an `AllostericMechanism` it is the joint solve's
count on the built child, which derives no rate polynomial (`_split_gain_test`). The
RE→SS move flips whole groups and accepts a
```

   then re-wrap the paragraph to 92 characters.
4. After the paragraph that ends ``…compile once per optimizer, not once per mechanism.``
   (the one Task 5 edited; it ends with the `_peak_turnover` sentence), add a
   subsection:

```markdown
### The joint thermodynamic solve

`_dependent_param_exprs` (`src/thermo.jl`) builds one integer system per mechanism
(`_joint_system`): a column per constant of every step (`_Column`), with an inactive
copy for each step of a `:NonequalAI` or `:OnlyA` group; a row per fundamental cycle
and conformation, and a row tying each further member of a kinetic class to its first.
It orders the columns by `_priority_key` and reduces the system exactly
(`_reduce_joint_system`: `Rational{Int}`, redone in `Rational{BigInt}` on overflow).
Pivot columns are dependents. Columns outside the inactive conformation's law (an
`:OnlyA` step's inactive copy, a step on a form the inactive conformation cannot reach)
come first, so no kept dependent reads them. Each class is named after its member that
comes last in that order (`_kept_members`), through the `name(p, m)` chokepoint, and
the engine derives each conformation's polynomials from `(reaction, groups, names)`:
the inactive one on `_inactive_groups(am)`, with the full mechanism's own steps.
```

- [ ] **Step 5: `.claude/CLAUDE.md`**

Replace the paragraph under `### Canonical Step Form (load-bearing guard)` with:

```markdown
Step direction, step order, and group order are canonicalized in the `Step` and `Mechanism`/`AllostericMechanism` constructors. This is **load-bearing, not cosmetic**: structural deduplication (`unique!`) works by pure `==`/`hash` only because construction is canonical, and every constant is named after the stored sides of its class's kept member. The joint solve (`src/thermo.jl`) picks dependent parameters by a fixed content priority (`_priority_key`), so `fitted_params` and the rate equation do not depend on step order, step direction or group labels; renaming a metabolite can still change which tied constant is dependent. Do not relax or reorder this without reading the Developer page in the docs.
```

Nothing else in CLAUDE.md changes.

- [ ] **Step 6: `ROADMAP.md`**

1. Delete the known issues **The inactive state can give an `:EqualAI` constant a
   second, fitted name.** and **The two conformations fold Wegscheider ties
   separately.**
2. Replace the known issue **Splitting a tied kinetic group changes fitted names and
   `eq_hash` without changing the model.** with:

```markdown
- **Splitting a tied kinetic group changes `eq_hash` without changing the model.**
  When a square of rapid-equilibrium bindings already ties two binding constants,
  splitting their kinetic group leaves the rate function and the fitted names unchanged,
  but the split part's constant becomes a dependent with its own line under
  `# Wegscheider constraints:`, and the law names it, so the printed equation and the
  dedup key (`_rate_eq_dedup_key`, `src/identify_rate_equation.jl`) differ, and identify
  fits both. The rename-proof key under Planned removes this.
```

3. In the Planned entry **Rewrite the derivation around the reaction graph.**, replace
   `Four PRs remain: one thermodynamic solve that` / `derives both conformations from one graph; the graph mechanism type` with
   `Three PRs remain: the graph mechanism type` (re-wrap).
4. Replace the Planned entry **Derive both conformations from one catalytic rate
   equation.** with:

```markdown
- **Take `:OnlyA` steps to their limit in the inactive conformation.**
  Both conformations come from one joint solve on one graph, but the inactive law still
  deletes every `:OnlyA` group and prunes the forms that deletion strands. An `:OnlyA`
  binding should lose its affinity for the inactive conformation instead (`K → ∞`, or
  `kon → 0` with `koff` fixed for a steady-state binding), and an `:OnlyA` chemistry step
  be blocked (both rate constants → 0, their ratio set by thermodynamics); each step of an
  `:OnlyA` group already keeps its own class. A tag set is admitted when an exact
  feasibility test passes and every form the free enzyme cannot reach in the inactive
  conformation lies behind a binding, on its bound side; the limit then equals deletion
  and pruning, and enumeration stops producing the `:OnlyA` ping-pong mechanisms this
  rejects. The conformations will flip at free E, which every mechanism must then
  contain; until then the constructor admits a second free form only in rapid
  equilibrium with the first (`_assert_one_flip_segment`), where the choice only rescales
  `L`. The rewrite's admission PR brings the limits, admission and the flip at free E;
  `test/allosteric_ground_truth.jl` is the oracle. Priority: high.
```

5. Delete the Planned entry **Build the Wegscheider rename into the constraint
   columns.**

- [ ] **Step 7: Build the docs**

The docs build runs every `@example` block and doctest and checks that every export has
a docstring (`docs/make.jl`); the test suite does not. Run, under the lock with
`timeout -s KILL 3600`:

```bash
timeout -s KILL 3600 flock -w 14400 /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock julia --project=docs --heap-size-hint=2500M -e 'using Pkg; Pkg.develop(PackageSpec(path = pwd())); Pkg.instantiate()'
timeout -s KILL 3600 flock -w 14400 /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock julia --project=docs --heap-size-hint=2500M docs/make.jl 2>&1 | tee /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase2/logs/docs.log | tail -40
```

Expected: `makedocs` finishes without an error (`deploydocs` skips outside CI). Read the
rendered output of `thermodynamic_constraints.md`'s example and `mwc_allostery.md`'s
`coupled` and `ss` examples in `docs/build/` and confirm the prose of Steps 1 and 3
matches them. Do not commit `docs/build/` or any `docs/Manifest.toml` change the
`Pkg.develop` call makes: `docs/Manifest.toml` is untracked, and if `git status` shows
`docs/Project.toml` modified, restore it with `git checkout docs/Project.toml`.

- [ ] **Step 8: Commit**

```bash
git status
git add docs/src/deriving/thermodynamic_constraints.md docs/src/deriving/kinetic_groups.md \
        docs/src/deriving/mwc_allostery.md docs/src/developer.md .claude/CLAUDE.md ROADMAP.md
git commit -m "Document the joint thermodynamic solve and retire its ROADMAP entries" \
  -m "The thermodynamic-constraints page states the fixed priority and that renaming a
metabolite can change which tied constant is dependent; the kinetic-groups and MWC pages
name the dependents the solve now picks; the developer page describes the solve; the
Canonical Step Form note no longer says the reduction follows step order. ROADMAP known
issues 1 and 2 and the rename-pass entry are retired, and the conformations entry keeps
its :OnlyA-limit half." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---

### Task 11: Re-measure the GC live-byte growth

**Files:**
- Modify: `ROADMAP.md` (the known issue **Repeated constraint solves inflate the GC's
  live-byte count.**)
- Create (scratch, not committed):
  `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase2/probes/gc_live_bytes_task11.jl`

**Interfaces:**
- Consumes: the Task 7 solve.

Spec section 18 item 11 and controller ruling 10: re-measure with the entry's own loop;
retire the entry only if the growth is gone, otherwise update its numbers. Measured on
the prototype of Task 7 (aarch64, Julia 1.12.5), 100,000 calls each: a plain random-order
bi-bi grows by nothing (the reading oscillates between 2 and 60 MB; no full sweep; 5.1 s,
against 4.85 GB, 123 full sweeps and 33 s before the switch); an allosteric bi-bi with
`:OnlyA` groups still grows by 1.5 GB with 17 full sweeps in 34 s (7.6 GB, 219 sweeps
and 76 s before), because lifting the compiled type constructs an `AllostericMechanism`,
whose `:OnlyA` guard runs its null space in `Rational{BigInt}`. So the entry is updated,
not retired.

- [ ] **Step 1: Write the measurement script**

Create the scratch file above with:

```julia
# The ROADMAP's GC live-byte loop: repeated _dependent_param_exprs calls on one mechanism,
# printing Base.gc_live_bytes() growth and full sweeps every 10,000 calls.
using EnzymeRates
const ER = EnzymeRates
function measure(label, f)
    f()
    GC.gc(true)
    live0, sweeps0, t0 = Base.gc_live_bytes(), Base.gc_num().full_sweep, time()
    for block in 1:10
        for _ in 1:10_000
            f()
        end
        println(label, " calls=", block * 10_000, " gc_live_MB=",
                round((Base.gc_live_bytes() - live0) / 2^20; digits = 1),
                " full_sweeps=", Base.gc_num().full_sweep - sweeps0,
                " wall_s=", round(time() - t0; digits = 1))
    end
end
plain = ER.Mechanism(@enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)
        E + B <--> E(B)
        E(A) + B <--> E(A, B)
        E(B) + A <--> E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) <--> E(P) + Q
        E(P, Q) <--> E(Q) + P
        E(P) <--> E + P
        E(Q) <--> E + Q
    end
end)
allo = typeof(@allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    catalytic_multiplicity: 2
    catalytic_steps: begin
        E + A <--> E(A)        :: OnlyA
        E + B <--> E(B)        :: EqualAI
        E(A) + B <--> E(A, B)  :: EqualAI
        E(B) + A <--> E(A, B)  :: OnlyA
        E(A, B) <--> E(P, Q)   :: OnlyA
        E(P, Q) <--> E(P) + Q  :: EqualAI
        E(P, Q) <--> E(Q) + P  :: EqualAI
        E(P) <--> E + P        :: EqualAI
        E(Q) <--> E + Q        :: EqualAI
    end
end)
measure("plain", () -> ER._dependent_param_exprs(plain))
measure("allosteric", () -> ER._dependent_param_exprs(allo))
```

- [ ] **Step 2: Run it**

`timeout -s KILL 1200 flock -w 14400 /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock julia --project=. --heap-size-hint=2500M /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase2/probes/gc_live_bytes_task11.jl`

Expected: the `plain` lines show no growth trend and `full_sweeps=0`; the `allosteric`
lines grow by roughly 15 KB per call with a few full sweeps. If the allosteric growth is
gone too, delete the entry instead of Step 3's rewrite and say so in the commit message.

- [ ] **Step 3: Update the entry**

Replace the known issue **Repeated constraint solves inflate the GC's live-byte count.**
with (put your measured numbers in place of these if they differ):

```markdown
- **The `:OnlyA` guard inflates the GC's live-byte count.**
  Each `AllostericMechanism` construction with an `:OnlyA` group runs
  `_onlya_haldane_violation` (`src/thermo.jl`), whose null space is computed in
  `Rational{BigInt}`, and so does every lift of a compiled allosteric type. 100,000
  `_dependent_param_exprs` calls on an allosteric random-order bi-bi raise
  `Base.gc_live_bytes` by about 1.5 GB that no live object holds, with 17 full sweeps; the
  joint solve itself runs in `Int64` and leaves no growth on a plain mechanism. Admission,
  which replaces the guard, runs in integers. A loop of `_dependent_param_exprs` calls on
  one allosteric mechanism that prints `Base.gc_live_bytes()` and
  `Base.gc_num().full_sweep` shows it.
```

- [ ] **Step 4: Commit**

```bash
git status
git add ROADMAP.md
git commit -m "Re-measure the GC live-byte growth after the joint solve" \
  -m "100,000 solves of a plain bi-bi no longer grow the live-byte count (4.85 GB and 123
full sweeps before); an allosteric one still grows by 1.5 GB (7.6 GB before), from the
:OnlyA guard's Rational{BigInt} null space, so the entry now names the guard." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---

### Task 12: Propose v0.11.0, verify, open the PR

**Files:**
- Modify: `Project.toml` (`version = "0.10.0"` → `version = "0.11.0"`)

Fitted names change (E3, E7, B10), `Full`/`Reduced` go and `all_parameters` arrives:
a breaking release before 1.0, for Denis to decide.

- [ ] **Step 1: Bump the version**

Edit `Project.toml`, then run the full suite. Expected: all pass (about 33,300 tests).

- [ ] **Step 2: Commit and push**

```bash
git status
git add Project.toml
git commit -m "Propose v0.11.0: the joint solve renames fitted constants" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
git push -u origin joint-solve
```

- [ ] **Step 3: Open the PR (draft, stacked)**

Open a draft PR with base `main` and head `joint-solve`, titled "Derive every law from
one joint thermodynamic solve (v0.11.0)". Body sections: **Why** (two rename passes
and a rebuilt inactive mechanism made fitted names depend on grouping and gave
`:EqualAI` constants phantom or folded names; spec sections 6 and 7), **What changes**
(one bullet per Task 1–11), **Breaking changes** (fitted names: dependents follow the
fixed priority, so many mechanisms fit a different but equivalent set, e.g. Segel's
ordered bi-bi fits `k_EA_to_E_A` and derives `k_EQ_P_to_EAB`; classes are named after
their kept member; `:EqualAI` regulator ligands are bare `K_<ligand>reg`; `Full`,
`Reduced` and the mode arguments are removed; `all_parameters` is exported; identify's
CSV parameter columns and `eq_hash` values change, and a same-model regrouping no longer
shares an `eq_hash`), **Known issues fixed** (ROADMAP known issues 1 and 2: 71 of the
2,344 mechanisms of the one-round ping-pong-capable bi-bi population lost a phantom or
folded parameter), **Stacking** (based on `fitting-scale`, PR #91, itself stacked on
#85–#89; the commits of this PR are those after `e5fea31`; it will be rebased onto main
once they merge), **Version** (v0.11.0 proposed for Denis to decide), **Tests** (the
exact oracle and its tolerances, the populations, the counts) and **Verification**
(full-suite count and time, the three `rate_equation` timings, CI results). End with
`🤖 Generated with [Claude Code](https://claude.com/claude-code)`; no session link.

- [ ] **Step 4: Wait for CI and report**

Poll `gh pr checks <n>` in a bounded loop (sleep 60, at most 60 minutes). Expected:
every job passes, the docs job included (it runs the `@example` blocks and doctests). If
a job fails, find the root cause with superpowers:systematic-debugging, fix it
test-first, push, and poll again.
