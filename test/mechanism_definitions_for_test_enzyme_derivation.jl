# ABOUTME: Shared mechanism test specifications used across the derivation tests.
# ABOUTME: Each mechanism is defined inline with its expected properties.

using LinearAlgebra
using Random

# ── _testhelper_MechanismTestSpec struct ────────────────────────────────────

"""
Data-driven test specification for enzyme mechanisms.
Contains all expected properties for comprehensive testing.
"""
Base.@kwdef struct _testhelper_MechanismTestSpec
    # Core data
    name::String                          # Human-readable name for test labels
    mechanism::Any                        # EnzymeMechanism or AllostericEnzymeMechanism
    metabolite_names::Vector{Symbol}      # metabolites(mechanism) order; for params/concs

    # Structural expectations
    expected_n_states::Int
    expected_n_steps::Int

    # Constraint expectations
    expected_n_haldane_constraints::Int       # RHS references Keq (catalytic-cycle closure)
    expected_n_mirror_constraints::Int = 0    # RHS is a single Symbol: a tie K_a = K_b
    expected_n_wegscheider_constraints::Int = 0   # RHS Expr without Keq (multi-cycle futile-cycle closure)
    expected_n_independent_params::Int        # 2*n_steps - n_constraints

    # Test configuration (optional)
    run_ode_test::Bool = true
    reference_rtol::Float64 = 1e-8
    ode_rtol::Float64 = 1e-6

    # Optional analytical rate function for extra validation (unicyclic mechanisms)
    # Signature: (all_params::NamedTuple, concs::NamedTuple) -> Float64
    analytical_rate_fn::Union{Function,Nothing} = nothing

    # Optional analytical peak turnover, which the peak search approaches on a wide range:
    # the kcat (saturating substrates, zero products), or, for a spec whose product raises
    # the turnover (m_OnlyA_prod), the supremum over product levels
    # Signature: (params::NamedTuple) -> Float64
    # params holds `_testhelper_written_params`' constants (the textbook names `k1f`,
    # `K2`, ... beside every name `all_parameters` gives) and `Et`
    analytical_kcat_fn::Union{Function,Nothing} = nothing

    # Factored form validation (optional — for mechanisms with known factored forms)
    # Expected numerator/denominator from rate_equation_string output
    expected_factored_num::Union{String,Nothing} = nothing
    expected_factored_denom::Union{String,Nothing} = nothing

    # The mechanism as written, the DSL walker's raw value: its steps and, for an
    # allosteric mechanism, its regulatory sites in written order, before the
    # constructors canonicalize them. Positional textbook oracles number k1,k2,… (and
    # reg1,reg2,…) by written order; `_testhelper_written_params` reads it to take each
    # written step's constants from its stored step, by role. The reaction-graph tests
    # build each spec's graph from it. Every spec fills it, through
    # `@_testhelper_with_raw`.
    raw::Union{EnzymeRates.RawMechanism,EnzymeRates.RawAllostericMechanism,Nothing} =
        nothing
end

"""
    @_testhelper_with_raw call

`(mechanism, raw)` for `call`, an `@enzyme_mechanism` or `@allosteric_mechanism` call:
the compiled mechanism the call builds, and the raw value the package's DSL walker
reads from its block, which keeps the steps and regulatory sites as written. It fills
`_testhelper_MechanismTestSpec.raw`.
"""
macro _testhelper_with_raw(call)
    allosteric = call isa Expr && call.head === :macrocall &&
                 call.args[1] === Symbol("@allosteric_mechanism")
    allosteric || call isa Expr && call.head === :macrocall &&
                  call.args[1] === Symbol("@enzyme_mechanism") ||
        error("@_testhelper_with_raw wraps an @enzyme_mechanism or " *
              "@allosteric_mechanism call; got $call")
    esc(:(($call, $(EnzymeRates._raw_mechanism(call.args[end], allosteric)))))
end

# Every declared substrate and product appears in some step: some step takes it up or
# gives it off. Regulators are exempt: init_mechanisms declares a dead-end inhibitor that
# no step binds yet (expand_mechanisms binds it later; _sig_of leaves it out of the
# compiled type). Substrates and products are never dropped.
function _testhelper_assert_mechanism_invariants(m)
    g = _testhelper_graph(m)
    appearing = Set(x for s in g.steps for x in [s.consumed; s.released])
    for x in [g.substrates; g.products]
        x in appearing || error("declared substrate/product $x appears in no step")
    end
    nothing
end

# The Mechanism lifted from an EnzymeMechanism, or the AllostericMechanism lifted from an
# AllostericEnzymeMechanism; a concrete mechanism is returned as is.
_testhelper_lift(m::EnzymeRates.EnzymeMechanism) = EnzymeRates.Mechanism(m)
_testhelper_lift(m::EnzymeRates.AllostericEnzymeMechanism) =
    EnzymeRates.AllostericMechanism(m)
_testhelper_lift(m::Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}) = m

"""
    _testhelper_graph(m) -> NamedTuple

`m` (concrete or compiled) as the test helpers read it, a `NamedTuple`: `forms`, each
form's name (`Vector{String}`, stored order); `free`, whether each form holds no
metabolite and no residual; `steps`, one `NamedTuple` per step in stored order and
direction with `from` and `to` (form positions), `consumed` and `released` (metabolite
names, a copy by its reactant's name), `consumed_modifier` and `released_modifier` (which
of them bind as modifiers), `rapid_equilibrium`, `group`, and `names` and
`inactive_names`, the package names of its constants in the active (or only) and the
inactive conformation; `tags`; `inactive`, which steps the inactive law keeps;
`catalytic_multiplicity`; `sites`, each `(; multiplicity, ligands, allo_states)`;
`substrates`, `products` and `modifiers`. A plain mechanism has no tags, no inactive
conformation (`inactive` and every `inactive_names` empty), no sites and
`catalytic_multiplicity` `nothing`. Forms come in the canonical form order: free E first,
then by content.
"""
function _testhelper_graph(m)
    mech = _testhelper_lift(m)
    allo = mech isa EnzymeRates.AllostericMechanism
    groups = EnzymeRates.steps(mech)
    forms = EnzymeRates._forms(groups)
    position = Dict(f => i for (i, f) in enumerate(forms))
    sides = EnzymeRates._class_sides(mech)
    tags = allo ? copy(EnzymeRates.cat_allo_states(mech)) : Symbol[]
    flat = EnzymeRates._flat_steps(groups)
    _testhelper_names(s, state) = Symbol[EnzymeRates.name(p, mech, sides)
                                         for p in EnzymeRates._step_constants(s, state)]
    steps = map(flat) do (s, g)
        shared = allo && tags[g] === :EqualAI
        (from = position[EnzymeRates.from_species(s)],
         to = position[EnzymeRates.to_species(s)],
         consumed = Symbol[EnzymeRates.name(x) for x in EnzymeRates.consumed(s)],
         released = Symbol[EnzymeRates.name(x) for x in EnzymeRates.released(s)],
         consumed_modifier = Bool[x isa EnzymeRates.CompetitiveInhibitor
                                  for x in EnzymeRates.consumed(s)],
         released_modifier = Bool[x isa EnzymeRates.CompetitiveInhibitor
                                  for x in EnzymeRates.released(s)],
         rapid_equilibrium = EnzymeRates.is_equilibrium(s), group = g,
         names = _testhelper_names(s, !allo ? :None : shared ? :EqualAI : :A),
         inactive_names = allo ? _testhelper_names(s, shared ? :EqualAI : :I) : Symbol[])
    end
    kept = allo ? EnzymeRates._inactive_groups(mech) : nothing
    sites = [(multiplicity = EnzymeRates.multiplicity(site),
              ligands = copy(EnzymeRates.ligands(site)),
              allo_states = copy(EnzymeRates.allo_states(site)))
             for site in (allo ? EnzymeRates.regulatory_sites(mech) :
                                 EnzymeRates.RegulatorySite[])]
    (forms = String[EnzymeRates._form_label(f, mech) for f in forms],
     free = Bool[isempty(EnzymeRates.bound(f)) && isempty(EnzymeRates.residual(f))
                 for f in forms],
     steps, tags,
     inactive = Int[j for (j, (s, g)) in enumerate(flat) if allo && s in kept[g]],
     catalytic_multiplicity = allo ? EnzymeRates.catalytic_multiplicity(mech) : nothing,
     sites, substrates = EnzymeRates.substrates(mech),
     products = EnzymeRates.products(mech), modifiers = EnzymeRates.modifiers(mech))
end

"""The name of every step constant of the plain mechanism `m`, fitted or dependent: each
step's constants as the naming chokepoint renders them."""
_testhelper_step_constant_names(m) =
    Set(x for s in _testhelper_graph(m).steps for x in s.names)

# The names of a compiled mechanism's substrates and products, and of its regulators: the
# modifiers an EnzymeMechanism's steps bind, an AllostericEnzymeMechanism's
# regulatory-site ligands.
_testhelper_substrates(m) = _testhelper_graph(m).substrates
_testhelper_products(m) = _testhelper_graph(m).products
_testhelper_regulators(m::EnzymeRates.EnzymeMechanism) = _testhelper_graph(m).modifiers
_testhelper_regulators(m::EnzymeRates.AllostericEnzymeMechanism) =
    unique(Symbol[l for site in _testhelper_graph(m).sites for l in site.ligands])

"""Generate random independent (fitted) params + Keq + E_total for a mechanism."""
function _testhelper_random_independent_params(m; rng=Random.default_rng())
    fp = EnzymeRates.fitted_params(m)
    vals = Tuple(0.1 + 9.9 * rand(rng) for _ in fp)
    Keq_val = 0.1 + 9.9 * rand(rng)
    E_total_val = 0.1 + 9.9 * rand(rng)
    keys_out = (fp..., :Keq, :E_total)
    vals_out = (vals..., Keq_val, E_total_val)
    NamedTuple{keys_out}(vals_out)
end

"""
A reaction over the names of `m`: its substrates and products with placeholder atoms that
balance for any counts (each substrate as many carbons as `m` has products, each product
as many as it has substrates), and its modifiers as competitive inhibitors. The
constructors over Step groups take a reaction for its names, so a helper that rebuilds a
mechanism passes this one.
"""
function _testhelper_placeholder_reaction(m)
    g = _testhelper_graph(m)
    reactants = EnzymeRates.ReactantAtoms[
        [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(s), [:C => length(g.products)])
         for s in g.substrates];
        [EnzymeRates.ReactantAtoms(EnzymeRates.Product(p), [:C => length(g.substrates)])
         for p in g.products]]
    EnzymeRates.EnzymeReaction(reactants, EnzymeRates.RegulatorMults[
        EnzymeRates.RegulatorMults(EnzymeRates.CompetitiveInhibitor(x), [1])
        for x in g.modifiers], [1])
end

"""
Test that `m` stores its steps in the direction the orientation tiers give: `m` written
back as raw input (`_raw`), each step as stored, has no step that canonical form turns
(`_canonicalize`), group by group, and writing every step the other way round, in the
reverse order, builds the same mechanism.
"""
function _testhelper_check_stored_orientation(m)
    mech = _testhelper_lift(m)
    raw = EnzymeRates._raw(mech)
    allo = raw isa EnzymeRates.RawAllostericMechanism
    cat = allo ? raw.catalytic : raw
    # `_raw` writes edge `k` as step `k`.
    flipped = EnzymeRates._canonicalize(cat).flipped
    for members in EnzymeRates._group_members(EnzymeRates.catalytic(mech))
        @test !any(flipped[members])
    end
    written = EnzymeRates.RawMechanism(cat.substrates, cat.products, cat.modifiers,
        reverse([EnzymeRates.RawStep(s.to, s.from, s.released, s.consumed,
                                     s.rapid_equilibrium) for s in cat.steps]),
        reverse(cat.group))
    rebuilt = !allo ? EnzymeRates.Mechanism(written) :
        EnzymeRates.AllostericMechanism(EnzymeRates.RawAllostericMechanism(written,
            raw.cat_allo_states, raw.catalytic_multiplicity, raw.sites))
    @test rebuilt == mech
end

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
        _testhelper_placeholder_reaction(am), EnzymeRates.steps(am),
        EnzymeRates.cat_allo_states(am), 1, EnzymeRates.regulatory_sites(am)))
end

"""
The free form at which the two conformations of an allosteric mechanism interconvert in
its law, as its position in `g` (`_testhelper_graph`): the free form (no bound
metabolite, no residual) first in the canonical form order (free E first), which roots
the rapid-equilibrium segment that holds every free form.
"""
_testhelper_flip_form(g) = findfirst(g.free)

"""Whether `s`, a step of `_testhelper_graph`, is a binding: it takes up one metabolite and
gives off none."""
_testhelper_is_binding(s) = length(s.consumed) == 1 && isempty(s.released)

"""Whether `s`, a step of `_testhelper_graph`, is an isomerization: it takes up and gives
off nothing."""
_testhelper_is_iso(s) = isempty(s.consumed) && isempty(s.released)

"""The first step of each kinetic group of `g` (`_testhelper_graph`), in group order."""
_testhelper_first_steps(g) = [g.steps[findfirst(s -> s.group == k, g.steps)]
                              for k in unique(s.group for s in g.steps)]

"""The text of `s`, a step of `g` (`_testhelper_graph`), as `_step_text` writes it, each
side its form's name and its metabolites' names: `"E + S ⇌ ES"`, `"ES <--> EP"`."""
_testhelper_step_text(g, s) =
    join([g.forms[s.from]; String.(s.consumed)], " + ") *
    (s.rapid_equilibrium ? " ⇌ " : " <--> ") *
    join([g.forms[s.to]; String.(s.released)], " + ")

"""
    _testhelper_consistent_point(m, rng) -> (; vals, Keq, active, inactive, L)

A thermodynamically consistent parameter point of `m`, built from free energies and so
independent of the package's thermodynamic solve. Each enzyme form has a free energy in
each conformation and each metabolite one; a step's association constant is
`exp(g_from + Σ g_consumed − g_to − Σ g_released)`, so every cycle multiplies to `Keq^n`.
The free energies are projected so that the steps of a kinetic group share their
constants (their metabolites cancel) and, in an allosteric mechanism, an `:EqualAI` step
has one constant in both conformations. A steady-state class draws one forward rate.

`vals` maps the name of every step constant of the law (`names`, and `inactive_names` of
the inactive law's steps, in `_testhelper_graph`), every regulator constant and `L` to its
value. `active[j]` (and `inactive[j]` for the inactive law's steps) holds a step of
`_testhelper_graph` as `(step, Ka, kf)`, `kf === nothing` for a rapid-equilibrium step.
"""
function _testhelper_consistent_point(m, rng)
    mech = _testhelper_lift(m)
    g = _testhelper_graph(mech)
    allo = g.catalytic_multiplicity !== nothing
    mets = unique!([g.substrates; g.products; g.modifiers])
    nf, nconf = length(g.forms), allo ? 2 : 1
    mi = Dict(x => nconf * nf + i for (i, x) in enumerate(mets))
    n = nconf * nf + length(mets)
    # The form part of step s's free-energy drop in conformation c (1 or 2).
    function _testhelper_drop(c, s)
        v = zeros(BigFloat, n)
        v[(c - 1) * nf + s.from] += 1
        v[(c - 1) * nf + s.to] -= 1
        v
    end
    tag(k) = allo ? g.tags[k] : :None
    ties = Vector{BigFloat}[]
    for (j, s) in enumerate(g.steps)
        j1 = findfirst(t -> t.group == s.group, g.steps)
        s1 = g.steps[j1]
        j == j1 || push!(ties, _testhelper_drop(1, s) - _testhelper_drop(1, s1))
        allo || continue
        tag(s.group) === :NonequalAI && j != j1 &&
            push!(ties, _testhelper_drop(2, s) - _testhelper_drop(2, s1))
        tag(s.group) === :EqualAI &&
            push!(ties, _testhelper_drop(2, s) - _testhelper_drop(1, s))
    end
    basis = Vector{BigFloat}[]
    for v in ties
        for b in basis; v = v - (b' * v) * b; end
        nv = sqrt(v' * v)
        nv > big(10)^-30 && push!(basis, v / nv)
    end
    u = BigFloat[4 * big(rand(rng)) - 2 for _ in 1:n]
    for b in basis; u = u - (b' * u) * b; end
    gm(x) = u[mi[x]]
    Keq = exp(sum(gm, g.substrates) - sum(gm, g.products))
    vals = Dict{Symbol, BigFloat}()
    kf_of = Dict{Symbol, BigFloat}()
    function _testhelper_conformation_steps(c, inactive)
        out = []
        for (j, s) in enumerate(g.steps)
            inactive && j ∉ g.inactive && continue
            Ka = exp(u[(c - 1) * nf + s.from] + sum(gm, s.consumed; init = big(0)) -
                     u[(c - 1) * nf + s.to] - sum(gm, s.released; init = big(0)))
            names = inactive ? s.inactive_names : s.names
            if s.rapid_equilibrium
                vals[only(names)] = _testhelper_is_binding(s) ? 1 / Ka : Ka
                push!(out, (s, Ka, nothing))
            else
                kf = get!(() -> exp(4 * big(rand(rng)) - 2), kf_of, names[1])
                vals[names[1]] = kf
                vals[names[2]] = kf / Ka
                push!(out, (s, Ka, kf))
            end
        end
        out
    end
    active = _testhelper_conformation_steps(1, false)
    inactive = allo ? _testhelper_conformation_steps(2, true) : []
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
        g = _testhelper_graph(mech)
        node = Dict{Tuple{Symbol, Int}, Int}()
        id(c, f) = get!(node, (c, f), length(node) + 1)
        x1 = first(g.substrates)
        edges = Tuple{Int, Int, BigFloat, BigFloat, Int}[]
        for (c, list) in ((:A, point.active), (:I, point.inactive)), (s, Ka, kf) in list
            f, r = kf === nothing ? (Λ, Λ / Ka) : (kf, kf / Ka)
            rf = f * prod((concs[x] for x in s.consumed); init = big(1))
            rr = r * prod((concs[x] for x in s.released); init = big(1))
            ω = count(==(x1), s.consumed) - count(==(x1), s.released)
            push!(edges, (id(c, s.from), id(c, s.to), rf, rr, ω))
        end
        if g.catalytic_multiplicity !== nothing
            # Regulators bind each conformation at rapid equilibrium, independently of
            # the catalytic form, so summed over regulator occupancy the flip ratio is
            # L times the ratio of the conformations' regulatory-site factors.
            function _testhelper_site_factor(inactive)
                prod((1 + sum((concs[lig] / point.vals[K]
                               for (lig, tag) in zip(site.ligands, site.allo_states)
                               for K in (EnzymeRates._reg_K(mech, lig, tag, inactive),)
                               if K !== nothing); init = big(0)))^site.multiplicity
                     for site in g.sites; init = big(1))
            end
            flip = _testhelper_flip_form(g)
            push!(edges, (id(:A, flip), id(:I, flip),
                          point.L * _testhelper_site_factor(true) /
                          _testhelper_site_factor(false), big(1), 0))
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

# ── Detailed balance, declared dimension and populations ────────────────────

"""The value of the polynomial `p` at the values `env` (by symbol), in BigFloat."""
_testhelper_poly_value(p, env) =
    sum((c * prod((env[s]^e for (s, e) in mono); init = big(1)) for (mono, c) in p);
        init = big(0))

"""
    _testhelper_detailed_balance(m, rng) -> Vector{BigFloat}

Detailed-balance residuals of `m`'s rate law at a random parameter point, in 256-bit
BigFloat: the law's value at equilibrium concentrations (the first product scaled so
that Π[P]/Π[S] = Keq) over its value with that product ten times lower, then, for an
allosteric mechanism, the same ratio for each conformation's numerator that is not
identically zero (`_state_rate_polys`). Each is zero in exact arithmetic when the law
satisfies every Haldane and Wegscheider relation.
"""
function _testhelper_detailed_balance(m, rng)
    setprecision(BigFloat, 256) do
        em = m isa EnzymeRates.AbstractEnzymeMechanism ? m :
             EnzymeRates.compile_mechanism(m)
        mech = _testhelper_lift(em)
        fitted = EnzymeRates._dependent_param_exprs(typeof(em))[2]
        params = Dict{Symbol, BigFloat}(k => exp(4 * big(rand(rng)) - 2) for k in fitted)
        params[:Keq] = exp(4 * big(rand(rng)) - 2)
        g = _testhelper_graph(mech)
        subs, prods = g.substrates, g.products
        concs = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1)
                                       for x in EnzymeRates.metabolites(em))
        p1 = first(prods)
        concs[p1] = params[:Keq] * prod(concs[x] for x in subs) /
                    prod((concs[x] for x in prods if x != p1); init = big(1))
        off = copy(concs)
        off[p1] /= 10
        ratios = BigFloat[abs(_testhelper_law_value(em, params, concs)) /
                          abs(_testhelper_law_value(em, params, off))]
        g.catalytic_multiplicity === nothing && return ratios
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

"""A basis of the cycle space of the step graph of `g` (`_testhelper_graph`): one integer
vector over its steps (each counted in its stored direction) per independent cycle, from
the null space of the form-by-step incidence matrix."""
function _testhelper_cycle_basis(g)
    R = zeros(Rational{BigInt}, length(g.forms), length(g.steps))
    for (j, s) in enumerate(g.steps)
        R[s.from, j] -= 1
        R[s.to, j] += 1
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
    g = _testhelper_graph(m)
    allo = g.catalytic_multiplicity !== nothing
    tag(k) = allo ? g.tags[k] : :None
    index = Dict{Any, Int}()
    in_law = Bool[]
    function _testhelper_class(key, inlaw)
        i = get!(index, key) do
            push!(in_law, false)
            length(in_law)
        end
        in_law[i] |= inlaw
        i
    end
    function _testhelper_classes(j, inactive)
        s = g.steps[j]
        slots = s.rapid_equilibrium ? (1,) : (1, 2)
        !inactive || tag(s.group) === :EqualAI ?
            [_testhelper_class((:A, s.group, k), true) for k in slots] :
        tag(s.group) === :OnlyA ? [_testhelper_class((:I, j, k), false) for k in slots] :
            [_testhelper_class((:I, s.group, k), j in g.inactive) for k in slots]
    end
    foreach(j -> _testhelper_classes(j, false), eachindex(g.steps))
    allo && foreach(j -> _testhelper_classes(j, true), eachindex(g.steps))
    rows = Vector{Int}[]
    for z in _testhelper_cycle_basis(g), inactive in (allo ? (false, true) : (false,))
        row = zeros(Int, length(in_law))
        for (j, s) in enumerate(g.steps)
            z[j] == 0 && continue
            cs = _testhelper_classes(j, inactive)
            if s.rapid_equilibrium
                row[cs[1]] += _testhelper_is_binding(s) ? -z[j] : z[j]
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
    n_reg = sum((t === :NonequalAI ? 2 : 1 for site in g.sites for t in site.allo_states);
                init = 0)
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

# ── Mechanism test specifications ───────────────────────────────────────────

const MECHANISM_TEST_SPECS = _testhelper_MechanismTestSpec[]

# 1. Uni-Uni (simplest): E + S ⇌ ES ⇌ E + P
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end
    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Uni-Uni",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P],
        expected_n_states=2,
        expected_n_steps=2,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=3,
        analytical_rate_fn=(p, c) -> begin
            num = p.k1f * p.k2f * c.S - p.k1r * p.k2r * c.P
            denom = p.k1r + p.k2f + p.k1f * c.S + p.k2r * c.P
            p.Et * num / denom
        end
    ))
end

# 2. Segel Uni Uni (replaces Three-Step Iso): E + A ⇌ EA ⇌ EP ⇌ E + P
#    Reference: Segel, Enzyme Kinetics, Eq. IX-8
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> E + P
        end
    end

    # Segel Eq. IX-8: Uni Uni steady-state rate
    function _testhelper_rate_uni_uni(params, concs)
        (; k1f, k1r, k2f, k2r, k3f, k3r, Etotal) = params
        (; A, P) = concs
        num = k1f * k2f * k3f * A - k1r * k2r * k3r * P
        denom = (k1r * k2r + k1r * k3f + k2f * k3f) +
                k1f * (k2r + k2f + k3f) * A +
                k3r * (k1r + k2r + k2f) * P
        return Etotal * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Segel Uni Uni",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :P],
        expected_n_states=3,
        expected_n_steps=3,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=5,
        analytical_rate_fn=(p, c) -> _testhelper_rate_uni_uni(merge(p, (Etotal=p.Et,)), c)
    ))
end

# 3. Segel Iso Uni Uni (new): E + A ⇌ EA ⇌ EP ⇌ F + P, F ⇌ E
#    Reference: Segel, Enzyme Kinetics, Eq. IX-45
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
        end
    end

    # Segel Eq. IX-45: Iso Uni Uni steady-state rate
    function _testhelper_rate_iso_uni_uni(params, concs)
        (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, Etotal) = params
        (; A, P) = concs
        num = k1f * k2f * k3f * k4f * A - k1r * k2r * k3r * k4r * P
        denom = (k4f + k4r) * (k1r * k3f + k1r * k2r + k2f * k3f) +
                k1f * (k2f * k3f + k2f * k4f + k2r * k4f + k3f * k4f) * A +
                k3r * (k1r * k2r + k1r * k4r + k2f * k4r + k2r * k4r) * P +
                k1f * k3r * (k2f + k2r) * A * P
        return Etotal * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Segel Iso Uni Uni",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :P],
        expected_n_states=4,
        expected_n_steps=4,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=7,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_iso_uni_uni(merge(p, (Etotal=p.Et,)), c)
    ))
end

# 4. Segel Ordered Uni Bi (replaces Seq Uni-Bi): E + A ⇌ (EA≡EPQ) ⇌ EQ + P ⇌ E + Q
#    Reference: Segel, Enzyme Kinetics, Eq. IX-60
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end

    # Segel Eq. IX-60: Ordered Uni Bi steady-state rate
    function _testhelper_rate_ordered_uni_bi(params, concs)
        (; k1f, k1r, k2f, k2r, k3f, k3r, Etotal) = params
        (; A, P, Q) = concs
        num = k1f * k2f * k3f * A - k1r * k2r * k3r * P * Q
        denom = k3f * (k2f + k1r) +
                k1f * (k2f + k3f) * A +
                k1r * k2r * P +
                k3r * (k2f + k1r) * Q +
                k1f * k2r * A * P +
                k2r * k3r * P * Q
        return Etotal * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Segel Ordered Uni Bi",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :P, :Q],
        expected_n_states=3,
        expected_n_steps=3,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=5,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_ordered_uni_bi(merge(p, (Etotal=p.Et,)), c)
    ))
end

# 5. Segel Ordered Bi Bi (replaces Seq Bi-Bi):
#    E + A ⇌ EA + B ⇌ (EAB≡EPQ) ⇌ EQ + P ⇌ E + Q
#    Reference: Segel, Enzyme Kinetics, Eq. IX-87
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end

    # Segel Eq. IX-87: Ordered Bi Bi steady-state rate
    function _testhelper_rate_ordered_bi_bi(params, concs)
        (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, Etotal) = params
        (; A, B, P, Q) = concs
        num = k1f * k2f * k3f * k4f * A * B - k1r * k2r * k3r * k4r * P * Q
        denom = k1r * k4f * (k2r + k3f) +
                k1f * k4f * (k2r + k3f) * A +
                k2f * k3f * k4f * B +
                k1r * k2r * k3r * P +
                k1r * k4r * (k2r + k3f) * Q +
                k1f * k2f * (k3f + k4f) * A * B +
                k1f * k2r * k3r * A * P +
                k2f * k3f * k4r * B * Q +
                k3r * k4r * (k1r + k2r) * P * Q +
                k1f * k2f * k3r * A * B * P +
                k2f * k3r * k4r * B * P * Q
        return Etotal * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Segel Ordered Bi Bi",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :B, :P, :Q],
        expected_n_states=4,
        expected_n_steps=4,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=7,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_ordered_bi_bi(merge(p, (Etotal=p.Et,)), c),
        analytical_kcat_fn=p -> p.k3f * p.k4f / (p.k3f + p.k4f),
    ))
end

# 6. Segel Theorell-Chance Bi Bi: E + A ⇌ EA, EA + B ⇌ EQ + P, EQ ⇌ E + Q
#    B binds and P leaves in one step, which carries metabolites on both sides.
#    Reference: Segel, Enzyme Kinetics, Eq. IX-122
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end

    # Segel Eq. IX-122: Theorell-Chance Bi Bi steady-state rate
    function _testhelper_rate_theorell_chance_bi_bi(params, concs)
        (; k1f, k1r, k2f, k2r, k3f, k3r, Etotal) = params
        (; A, B, P, Q) = concs
        num = k1f * k2f * k3f * A * B - k1r * k2r * k3r * P * Q
        denom = k1r * k3f +
                k1f * k3f * A +
                k2f * k3f * B +
                k1r * k2r * P +
                k1r * k3r * Q +
                k1f * k2f * A * B +
                k1f * k2r * A * P +
                k2f * k3r * B * Q +
                k2r * k3r * P * Q
        return Etotal * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Segel Theorell-Chance Bi Bi",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :B, :P, :Q],
        expected_n_states=3,
        expected_n_steps=3,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=5,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_theorell_chance_bi_bi(merge(p, (Etotal=p.Et,)), c),
        analytical_kcat_fn=p -> p.k3f,
    ))
end

# 7. Segel Ping Pong Bi Bi (replaces Ping-Pong Bi-Bi):
#    E + A ⇌ (EA≡FP) ⇌ F + P, F + B ⇌ (FB≡EQ) ⇌ E + Q
#    Reference: Segel, Enzyme Kinetics, Eq. IX-140
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> F + P
            F + B <--> F(B)
            F(B) <--> E + Q
        end
    end

    # Segel Eq. IX-140: Ping Pong Bi Bi steady-state rate
    function _testhelper_rate_ping_pong_bi_bi(params, concs)
        (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, Etotal) = params
        (; A, B, P, Q) = concs
        num = k1f * k2f * k3f * k4f * A * B - k1r * k2r * k3r * k4r * P * Q
        denom = k1f * k2f * (k3r + k4f) * A +
                k3f * k4f * (k1r + k2f) * B +
                k1r * k2r * (k3r + k4f) * P +
                k3r * k4r * (k1r + k2f) * Q +
                k1f * k3f * (k2f + k4f) * A * B +
                k1f * k2r * (k3r + k4f) * A * P +
                k3f * k4r * (k1r + k2f) * B * Q +
                k2r * k4r * (k1r + k3r) * P * Q
        return Etotal * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Segel Ping Pong Bi Bi",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :B, :P, :Q],
        expected_n_states=4,
        expected_n_steps=4,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=7,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_ping_pong_bi_bi(
                merge(p, (Etotal=p.Et,)), c),
        analytical_kcat_fn=p -> p.k2f * p.k4f / (p.k2f + p.k4f),
    ))
end

# 8. Segel Ordered Ter Bi (replaces Seq Ter-Bi):
#     E + A ⇌ EA + B ⇌ EAB + C ⇌ (EABC≡EPQ) ⇌ EQ + P ⇌ E + Q
#     Reference: Segel, Enzyme Kinetics, Eq. IX-195
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) + C <--> E(A, B, C)
            E(A, B, C) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end

    # Segel Eq. IX-195: Ordered Ter Bi steady-state rate
    function _testhelper_rate_ordered_ter_bi(params, concs)
        (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, Etotal) = params
        (; A, B, C, P, Q) = concs
        num = k1f * k2f * k3f * k4f * k5f * A * B * C -
              k1r * k2r * k3r * k4r * k5r * P * Q
        denom = # constant
            k1r * k2r * k5f * (k3r + k4f) +
            # single substrate/product
            k1f * k2r * k5f * (k3r + k4f) * A +
            k1r * k3f * k4f * k5f * C +
            k1r * k2r * k3r * k4r * P +
            k1r * k2r * k5r * (k3r + k4f) * Q +
            # two substrates/products
            k1f * k2f * k5f * (k3r + k4f) * A * B +
            k1f * k3f * k4f * k5f * A * C +
            k2f * k3f * k4f * k5f * B * C +
            k1f * k2r * k3r * k4r * A * P +
            k1r * k3f * k4f * k5r * C * Q +
            k4r * k5r * (k1r * k2r + k1r * k3r + k2r * k3r) * P * Q +
            # three substrates/products
            k1f * k2f * k3f * (k4f + k5f) * A * B * C +
            k1f * k2f * k3r * k4r * A * B * P +
            k2f * k3f * k4f * k5r * B * C * Q +
            k1r * k3f * k4r * k5r * C * P * Q +
            k2f * k3r * k4r * k5r * B * P * Q +
            # four substrates/products
            k1f * k2f * k3f * k4r * A * B * C * P +
            k2f * k3f * k4r * k5r * B * C * P * Q
        return Etotal * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Segel Ordered Ter Bi",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :B, :C, :P, :Q],
        expected_n_states=5,
        expected_n_steps=5,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=9,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_ordered_ter_bi(merge(p, (Etotal=p.Et,)), c)
    ))
end

# 9. Segel Ordered Ter Ter (replaces Seq Ter-Ter):
#     E + A ⇌ EA + B ⇌ EAB + C ⇌ (EABC≡EPQR) ⇌ EQR + P ⇌ ER + Q ⇌ E + R
#     Reference: Segel, Enzyme Kinetics, Eq. IX-261
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q, R
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) + C <--> E(A, B, C)
            E(A, B, C) <--> E(Q, R) + P
            E(Q, R) <--> E(R) + Q
            E(R) <--> E + R
        end
    end

    # Segel Eq. IX-261: Ordered Ter Ter steady-state rate
    function _testhelper_rate_ordered_ter_ter(params, concs)
        (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal) = params
        (; A, B, C, P, Q, R) = concs
        num = k1f * k2f * k3f * k4f * k5f * k6f * A * B * C -
              k1r * k2r * k3r * k4r * k5r * k6r * P * Q * R
        denom = # constant
            k1r * k2r * k5f * k6f * (k3r + k4f) +
            # single substrate/product
            k1f * k2r * k5f * k6f * (k3r + k4f) * A +
            k1r * k3f * k4f * k5f * k6f * C +
            k1r * k2r * k3r * k4r * k6f * P +
            k1r * k2r * k5f * k6r * (k3r + k4f) * R +
            # two substrates/products
            k1f * k2f * k5f * k6f * (k3r + k4f) * A * B +
            k1f * k3f * k4f * k5f * k6f * A * C +
            k1f * k2r * k3r * k4r * k6f * A * P +
            k2f * k3f * k4f * k5f * k6f * B * C +
            k1r * k3f * k4f * k5f * k6r * C * R +
            k1r * k2r * k3r * k4r * k5r * P * Q +
            k1r * k2r * k3r * k4r * k6r * P * R +
            k1r * k2r * k5r * k6r * (k3r + k4f) * Q * R +
            # three substrates/products
            k1f * k2f * k3f * (k4f * k5f + k4f * k6f + k5f * k6f) * A * B * C +
            k1f * k2f * k3r * k4r * k6f * A * B * P +
            k1f * k2r * k3r * k4r * k5r * A * P * Q +
            k2f * k3f * k4f * k5f * k6r * B * C * R +
            k1r * k3f * k4f * k5r * k6r * C * Q * R +
            k4r * k5r * k6r * (k1r * k2r + k1r * k3r + k2r * k3r) * P * Q * R +
            # four substrates/products
            k1f * k2f * k3f * k4r * k6f * A * B * C * P +
            k1f * k2f * k3f * k4f * k5r * A * B * C * Q +
            k1f * k2f * k3r * k4r * k5r * A * B * P * Q +
            k2f * k3f * k4f * k5r * k6r * B * C * Q * R +
            k2f * k3r * k4r * k5r * k6r * B * P * Q * R +
            k1r * k3f * k4r * k5r * k6r * C * P * Q * R +
            # five substrates/products
            k1f * k2f * k3f * k4r * k5r * A * B * C * P * Q +
            k2f * k3f * k4r * k5r * k6r * B * C * P * Q * R
        return Etotal * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Segel Ordered Ter Ter",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :B, :C, :P, :Q, :R],
        expected_n_states=6,
        expected_n_steps=6,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=11,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_ordered_ter_ter(
                merge(p, (Etotal=p.Et,)), c)
    ))
end

# 10. Segel Bi Uni Uni Uni Ping Pong Ter Bi (new):
#     E + A ⇌ EA + B ⇌ (EAB≡FP) ⇌ F + P, F + C ⇌ (FC≡EQ) ⇌ E + Q
#     Reference: Segel, Enzyme Kinetics, Eq. IX-228
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> F + P
            F + C <--> F(C)
            F(C) <--> E + Q
        end
    end

    # Segel Eq. IX-228: Bi Uni Uni Uni Ping Pong Ter Bi steady-state rate
    function _testhelper_rate_bi_uni_uni_uni_ping_pong_ter_bi(params, concs)
        (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, Etotal) = params
        (; A, B, C, P, Q) = concs
        num = k1f * k2f * k3f * k4f * k5f * A * B * C -
              k1r * k2r * k3r * k4r * k5r * P * Q
        denom = # single substrate/product
            k1r * k4f * k5f * (k2r + k3f) * C +
            k1r * k2r * k3r * (k4r + k5f) * P +
            k1r * k4r * k5r * (k2r + k3f) * Q +
            # two substrates/products
            k1f * k2f * k3f * (k4r + k5f) * A * B +
            k1f * k4f * k5f * (k2r + k3f) * A * C +
            k1f * k2r * k3r * (k4r + k5f) * A * P +
            k2f * k3f * k4f * k5f * B * C +
            k2f * k3f * k4r * k5r * B * Q +
            k1r * k4f * k5r * (k2r + k3f) * C * Q +
            k3r * k5r * (k1r * k2r + k1r * k4r + k2r * k4r) * P * Q +
            # three substrates/products
            k1f * k2f * k4f * (k3f + k5f) * A * B * C +
            k1f * k2f * k3r * (k4r + k5f) * A * B * P +
            k2f * k3f * k4f * k5r * B * C * Q +
            k2f * k3r * k4r * k5r * B * P * Q
        return Etotal * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Segel Bi Uni Uni Uni PP Ter Bi",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :B, :C, :P, :Q],
        expected_n_states=5,
        expected_n_steps=5,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=9,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_bi_uni_uni_uni_ping_pong_ter_bi(
                merge(p, (Etotal=p.Et,)), c)
    ))
end

# 11. Random-order Bi-Bi (branched): Two substrate binding orders converge
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E + B <--> E(B)
            E(A) + B <--> E(A, B)
            E(B) + A <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end
    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Random-order Bi-Bi",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :B, :P, :Q],
        expected_n_states=6,
        expected_n_steps=7,
        expected_n_haldane_constraints=1,
        expected_n_wegscheider_constraints=1,
        expected_n_independent_params=12,
    ))
end

# 11b. Random-order, chemistry-SS Bi-Bi: two substrate binding orders,
#      a steady-state chemistry step, single product. ODE is ground truth.
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S1, S2
        products: P
        steps: begin
            E + S1 ⇌ E(S1)
            E + S2 ⇌ E(S2)
            E(S1) + S2 <--> E(S1, S2)
            E(S2) + S1 <--> E(S1, S2)
            E(S1, S2) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Numerator: random chem-SS (RE/SS)",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S1, :S2, :P],
        expected_n_states=5,
        expected_n_steps=6,
        expected_n_haldane_constraints=1,
        expected_n_wegscheider_constraints=1,
        expected_n_independent_params=7,
    ))
end

# 11c. Catalytic isomerization RE, binding+release SS: random-order binding
#      and product release around a rapid-equilibrium chemistry step.
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S1, S2
        products: P1, P2
        steps: begin
            E + S1 ⇌ E(S1)
            E + S2 ⇌ E(S2)
            E(S1) + S2 <--> E(S1, S2)
            E(S2) + S1 <--> E(S1, S2)
            E(S1, S2) ⇌ E(P1, P2)
            E(P1, P2) <--> E(P1) + P2
            E(P1, P2) <--> E(P2) + P1
            E(P1) ⇌ E + P1
            E(P2) ⇌ E + P2
        end
    end
    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Numerator: RE-chemistry",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S1, :S2, :P1, :P2],
        expected_n_states=7,
        expected_n_steps=9,
        expected_n_haldane_constraints=1,
        expected_n_wegscheider_constraints=2,
        expected_n_independent_params=10,
    ))
end

# 11d. Single-RE-segment with a redundant SS binding branch (one branch
#      RE, one SS). ODE is ground truth.
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S1, S2
        products: P
        steps: begin
            E + S1 <--> E(S1)
            E + S2 ⇌ E(S2)
            E(S1) + S2 ⇌ E(S1, S2)
            E(S2) + S1 ⇌ E(S1, S2)
            E(S1, S2) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Numerator: redundant SS-bind",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S1, :S2, :P],
        expected_n_states=5,
        expected_n_steps=6,
        expected_n_haldane_constraints=1,
        expected_n_wegscheider_constraints=1,
        expected_n_independent_params=6,
    ))
end

# 12. Segel Bi Uni Uni Bi Ping Pong Ter Ter (new):
#     E + A ⇌ EA + B ⇌ (EAB≡FP) ⇌ F + P, F + C ⇌ (FC≡EQR) ⇌ ER + Q ⇌ E + R
#     Reference: Segel, Enzyme Kinetics, Eq. IX-278
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q, R
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> F + P
            F + C <--> F(C)
            F(C) <--> E(R) + Q
            E(R) <--> E + R
        end
    end

    # Segel Eq. IX-278: Bi Uni Uni Bi Ping Pong Ter Ter steady-state rate
    function _testhelper_rate_bi_uni_uni_bi_ping_pong_ter_ter(params, concs)
        (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal) = params
        (; A, B, C, P, Q, R) = concs
        num = k1f * k2f * k3f * k4f * k5f * k6f * A * B * C -
              k1r * k2r * k3r * k4r * k5r * k6r * P * Q * R
        denom = # single substrate/product
            k1r * k4f * k5f * k6f * (k2r + k3f) * C +
            k1r * k2r * k3r * k6f * (k4r + k5f) * P +
            # two substrates/products
            k1f * k2f * k3f * k6f * (k4r + k5f) * A * B +
            k1f * k4f * k5f * k6f * (k2r + k3f) * A * C +
            k1f * k2r * k3r * k6f * (k4r + k5f) * A * P +
            k2f * k3f * k4f * k5f * k6f * B * C +
            k1r * k4f * k5f * k6r * (k2r + k3f) * C * R +
            k1r * k2r * k3r * k4r * k5r * P * Q +
            k1r * k2r * k3r * k6r * (k4r + k5f) * P * R +
            k1r * k4r * k5r * k6r * (k2r + k3f) * Q * R +
            # three substrates/products
            k1f * k2f * k4f * (k3f * k5f + k3f * k6f + k5f * k6f) * A * B * C +
            k1f * k2f * k3r * k6f * (k4r + k5f) * A * B * P +
            k1f * k2f * k3f * k4r * k5r * A * B * Q +
            k1f * k2r * k3r * k4r * k5r * A * P * Q +
            k2f * k3f * k4f * k5f * k6r * B * C * R +
            k2f * k3f * k4r * k5r * k6r * B * Q * R +
            k1r * k4f * k5r * k6r * (k2r + k3f) * C * Q * R +
            k3r * k5r * k6r * (k1r * k2r + k1r * k4r + k2r * k4r) * P * Q * R +
            # four substrates/products
            k1f * k2f * k3f * k4f * k5r * A * B * C * Q +
            k1f * k2f * k3r * k4r * k5r * A * B * P * Q +
            k2f * k3f * k4f * k5r * k6r * B * C * Q * R +
            k2f * k3r * k4r * k5r * k6r * B * P * Q * R
        return Etotal * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Segel Bi Uni Uni Bi PP Ter Ter",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :B, :C, :P, :Q, :R],
        expected_n_states=6,
        expected_n_steps=6,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=11,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_bi_uni_uni_bi_ping_pong_ter_ter(
                merge(p, (Etotal=p.Et,)), c)
    ))
end

# 13. Segel Bi Bi Uni Uni Ping Pong Ter Ter (new):
#     E + A ⇌ EA + B ⇌ (EAB≡FPQ) ⇌ FQ + P, FQ ⇌ F + Q, F + C ⇌ (FC≡ER) ⇌ E + R
#     Reference: Segel, Enzyme Kinetics, Eq. IX-288
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q, R
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> F(Q) + P
            F(Q) <--> F + Q
            F + C <--> F(C)
            F(C) <--> E + R
        end
    end

    # Segel Eq. IX-288: Bi Bi Uni Uni Ping Pong Ter Ter steady-state rate
    function _testhelper_rate_bi_bi_uni_uni_ping_pong_ter_ter(params, concs)
        (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal) = params
        (; A, B, C, P, Q, R) = concs
        num = k1f * k2f * k3f * k4f * k5f * k6f * A * B * C -
              k1r * k2r * k3r * k4r * k5r * k6r * P * Q * R
        denom = # single substrate/product
            k1r * k4f * k5f * k6f * (k2r + k3f) * C +
            k1r * k4f * k5r * k6r * (k2r + k3f) * R +
            # two substrates/products
            k1f * k2f * k3f * k4f * (k5r + k6f) * A * B +
            k1f * k4f * k5f * k6f * (k2r + k3f) * A * C +
            k2f * k3f * k4f * k5f * k6f * B * C +
            k2f * k3f * k4f * k5r * k6r * B * R +
            k1r * k2r * k3r * k5f * k6f * C * P +
            k1r * k4f * k5f * k6r * (k2r + k3f) * C * R +
            k1r * k2r * k3r * k4r * (k5r + k6f) * P * Q +
            k1r * k2r * k3r * k5r * k6r * P * R +
            k1r * k4r * k5r * k6r * (k2r + k3f) * Q * R +
            # three substrates/products
            k1f * k2f * k5f * (k3f * k4f + k3f * k6f + k4f * k6f) * A * B * C +
            k1f * k2f * k3f * k4r * (k5r + k6f) * A * B * Q +
            k1f * k2r * k3r * k5f * k6f * A * C * P +
            k1f * k2r * k3r * k4r * (k5r + k6f) * A * P * Q +
            k2f * k3f * k4f * k5f * k6r * B * C * R +
            k2f * k3f * k4r * k5r * k6r * B * Q * R +
            k1r * k2r * k3r * k5f * k6r * C * P * R +
            k3r * k4r * k6r * (k1r * k2r + k1r * k5r + k2r * k5r) * P * Q * R +
            # four substrates/products
            k1f * k2f * k3r * k5f * k6f * A * B * C * P +
            k1f * k2f * k3r * k4r * (k5r + k6f) * A * B * P * Q +
            k2f * k3r * k4r * k5r * k6r * B * P * Q * R
        return Etotal * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Segel Bi Bi Uni Uni PP Ter Ter",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :B, :C, :P, :Q, :R],
        expected_n_states=6,
        expected_n_steps=6,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=11,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_bi_bi_uni_uni_ping_pong_ter_ter(
                merge(p, (Etotal=p.Et,)), c)
    ))
end

# 14. Segel Hexa Uni Ping Pong (new):
#     E + A ⇌ (EA≡FP) ⇌ F + P, F + B ⇌ (FB≡GQ) ⇌ G + Q, G + C ⇌ (GC≡ER) ⇌ E + R
#     Reference: Segel, Enzyme Kinetics, Eq. IX-308
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q, R
        steps: begin
            E + A <--> E(A)
            E(A) <--> F + P
            F + B <--> F(B)
            F(B) <--> G + Q
            G + C <--> G(C)
            G(C) <--> E + R
        end
    end

    # Segel Eq. IX-308: Hexa Uni Ping Pong steady-state rate
    function _testhelper_rate_hexa_uni_ping_pong(params, concs)
        (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal) = params
        (; A, B, C, P, Q, R) = concs
        num = k1f * k2f * k3f * k4f * k5f * k6f * A * B * C -
              k1r * k2r * k3r * k4r * k5r * k6r * P * Q * R
        denom = # two substrates/products (9 classes, each with 2 terms → factored)
            k1f * k2f * k3f * k4f * (k5r + k6f) * A * B +
            k1f * k2f * k5f * k6f * (k3r + k4f) * A * C +
            k1f * k2f * k3r * k4r * (k5r + k6f) * A * Q +
            k3f * k4f * k5f * k6f * (k1r + k2f) * B * C +
            k3f * k4f * k5r * k6r * (k1r + k2f) * B * R +
            k1r * k2r * k5f * k6f * (k3r + k4f) * C * P +
            k1r * k2r * k3r * k4r * (k5r + k6f) * P * Q +
            k1r * k2r * k5r * k6r * (k3r + k4f) * P * R +
            k3r * k4r * k5r * k6r * (k1r + k2f) * Q * R +
            # three substrates/products (8 classes)
            k1f * k3f * k5f * (k2f * k4f + k2f * k6f + k4f * k6f) * A * B * C +
            k1f * k2f * k3f * k4r * (k5r + k6f) * A * B * Q +
            k1f * k2r * k5f * k6f * (k3r + k4f) * A * C * P +
            k1f * k2r * k3r * k4r * (k5r + k6f) * A * P * Q +
            k3f * k4f * k5f * k6r * (k1r + k2f) * B * C * R +
            k3f * k4r * k5r * k6r * (k1r + k2f) * B * Q * R +
            k1r * k2r * k5f * k6r * (k3r + k4f) * C * P * R +
            k2r * k4r * k6r * (k1r * k3r + k1r * k5r + k3r * k5r) * P * Q * R
        return Etotal * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Segel Hexa Uni Ping Pong",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :B, :C, :P, :Q, :R],
        expected_n_states=6,
        expected_n_steps=6,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=11,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_hexa_uni_ping_pong(
                merge(p, (Etotal=p.Et,)), c)
    ))
end

# ── Rapid-Equilibrium (RE) Mechanisms ──────────────────────────────────────

# 15. RE Uni-Uni: E + A ⇌_RE EA <-->_SS E + P
#     Rapid-equilibrium substrate binding, steady-state catalysis
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E + P
        end
    end

    # Rapid-equilibrium Michaelis-Menten (K1 = Kd = [E][A]/[EA]):
    # rate = E_t * (k2f * A/K1 - k2r * P) / (1 + A/K1)
    function _testhelper_rate_re_uni_uni(params, concs)
        (; K1, k2f, k2r, Et) = params
        (; A, P) = concs
        num = k2f * A / K1 - k2r * P
        denom = 1.0 + A / K1
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="RE Uni-Uni",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :P],
        expected_n_states=2,
        expected_n_steps=2,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=2,
        analytical_rate_fn=(p, c) -> _testhelper_rate_re_uni_uni(merge(p, (Et=p.Et,)), c),
        analytical_kcat_fn=p -> p.k2f,
    ))
end

# 16. RE Ordered Bi-Bi: substrate binding is RE, catalysis and product release are SS
#     E + A ⇌_RE EA, EA + B ⇌_RE EAB, (EAB≡EPQ) <-->_SS EQ + P, EQ <-->_SS E + Q
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end

    # Groups: {E, EA, EAB} (RE group) and {EQ} (singleton)
    # K1, K2 are Kd (dissociation constants):
    #   σ = 1 + A/K1 + A*B/(K1*K2)
    # num = k3f*k4f*A*B/(K1*K2) - k3r*k4r*P*Q
    # denom = (1+A/K1+A*B/(K1*K2))*(k3r*P+k4f) + k3f*A*B/(K1*K2) + k4r*Q
    function _testhelper_rate_re_ordered_bi_bi(params, concs)
        (; K1, K2, k3f, k3r, k4f, k4r, Et) = params
        (; A, B, P, Q) = concs
        num = k3f * k4f * A * B / (K1 * K2) - k3r * k4r * P * Q
        sigma1 = 1.0 + A / K1 + A * B / (K1 * K2)
        R12 = k3f * A * B / (K1 * K2) + k4r * Q
        R21 = k3r * P + k4f
        denom = sigma1 * R21 + R12
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="RE Ordered Bi-Bi",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :B, :P, :Q],
        expected_n_states=4,
        expected_n_steps=4,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=5,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_re_ordered_bi_bi(merge(p, (Et=p.Et,)), c)
    ))
end

# 17. RE Random Bi-Bi: binding to free enzyme is RE, other steps are SS
#     E + A ⇌_RE EA, E + B ⇌_RE EB, EA + B <-->_SS EAB, EB + A <-->_SS EAB,
#     EAB <-->_SS EPQ, EPQ <-->_SS EQ + P, EQ <-->_SS E + Q
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E + B ⇌ E(B)
            E(A) + B <--> E(A, B)
            E(B) + A <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end
    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="RE Random Bi-Bi",
        mechanism=m,
        raw=raw,
        metabolite_names=[:A, :B, :P, :Q],
        expected_n_states=6,
        expected_n_steps=7,
        expected_n_haldane_constraints=1,
        expected_n_wegscheider_constraints=1,
        expected_n_independent_params=10,
    ))
end

# ── Classical Inhibitor/Activator Mechanisms (factored form tests) ────────

# 18. Competitive inhibitor: E + S ⇌ ES, ES ⇌ EP (SS), EP ⇌ E + P, E + R ⇌ ER
#     Dead-end inhibitor R binds free enzyme only.
#     No Cartesian product structure → flat sum denominator.
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: R
        steps: begin
            E + S ⇌ E(S)      # K1
            E(S) <--> E(P)    # k2f, k2r (SS)
            E(P) ⇌ E + P      # K3
            E + R ⇌ E(R)      # K4
        end
    end

    # v = Et * (k2f*S/K1 - k2r*P/K3) / (1 + S/K1 + P/K3 + R/K4)
    function _testhelper_rate_competitive_inh(params, concs)
        (; K1, k2f, k2r, K3, K4, Et) = params
        (; S, P, R) = concs
        num = k2f * S / K1 - k2r * P / K3
        denom = 1.0 + S / K1 + P / K3 + R / K4
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Competitive Inhibitor",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P, :R],
        expected_n_states=4,
        expected_n_steps=4,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=4,
        analytical_rate_fn=(p, c) -> _testhelper_rate_competitive_inh(
            merge(p, (Et=p.Et,)), c),
        analytical_kcat_fn=p -> p.k2f,
        # Textbook: flat sum denominator (no Cartesian product structure)
        expected_factored_num=
        "k_ES_to_EP * S / K_ES_to_E_S - k_EP_to_ES * P / K_EP_to_E_P",
        expected_factored_denom=
        "1 + P / K_EP_to_E_P + R / K_ER_to_E_R + S / K_ES_to_E_S",
    ))
end

# 19. Non-competitive inhibitor: R binds both free E and ES with same K
#     Forms: E, E_S, E_P, E_R, E_S_R; SS: E_S↔E_P; K5=K4
#     Denom factors as (1+R/K4)*(1+S/K1) + P/K3
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: R
        steps: begin
            E + S ⇌ E(S)        # K1
            E(S) <--> E(P)      # k2f, k2r (SS)
            E + P ⇌ E(P)        # K3
            (E + R ⇌ E(R), E(S) + R ⇌ E(S, R))  # K4 = K5 (R binding shared)
            E(R) + S ⇌ E(S, R)  # K6
        end
    end

    # v = Et * (k2f*S/K1 - k2r*P/K3) / ((1+R/K4)*(1+S/K1) + P/K3)
    function _testhelper_rate_noncompetitive_inh(params, concs)
        (; K1, k2f, k2r, K3, K4, Et) = params
        (; S, P, R) = concs
        num = k2f * S / K1 - k2r * P / K3
        denom = (1.0 + R / K4) * (1.0 + S / K1) + P / K3
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Non-competitive Inhibitor",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P, :R],
        expected_n_states=5,
        expected_n_steps=6,
        expected_n_haldane_constraints=1,
        expected_n_mirror_constraints=1,
        expected_n_independent_params=4,
        analytical_rate_fn=(p, c) -> _testhelper_rate_noncompetitive_inh(
            merge(p, (Et=p.Et,)), c),
        expected_factored_num=
        "k_ES_to_EP * S / K_ES_to_E_S - k_EP_to_ES * P / K_EP_to_E_P",
        expected_factored_denom=
        "1 + P / K_EP_to_E_P + R / K_ER_to_E_R + S / K_ES_to_E_S" *
        " + R * S / (K_ERS_to_ER_S * K_ER_to_E_R)",
    ))
end

# 20. Uncompetitive inhibitor: R binds ES only (not free E)
#     Forms: E, E_S, E_P, E_S_R; SS: E_S↔E_P; No extra constraints
#     Denom: 1 + P/K3 + S/K1*(1+R/K4)
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: R
        steps: begin
            E + S ⇌ E(S)        # K1
            E(S) <--> E(P)      # k2f, k2r (SS)
            E + P ⇌ E(P)        # K3
            E(S) + R ⇌ E(S, R)  # K4
        end
    end

    # v = Et * (k2f*S/K1 - k2r*P/K3) / (1 + P/K3 + S/K1*(1+R/K4))
    function _testhelper_rate_uncompetitive_inh(params, concs)
        (; K1, k2f, k2r, K3, K4, Et) = params
        (; S, P, R) = concs
        num = k2f * S / K1 - k2r * P / K3
        denom = 1.0 + P / K3 + (S / K1) * (1.0 + R / K4)
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Uncompetitive Inhibitor",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P, :R],
        expected_n_states=4,
        expected_n_steps=4,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=4,
        analytical_rate_fn=(p, c) -> _testhelper_rate_uncompetitive_inh(
            merge(p, (Et=p.Et,)), c),
        expected_factored_num=
        "k_ES_to_EP * S / K_ES_to_E_S - k_EP_to_ES * P / K_EP_to_E_P",
        expected_factored_denom=
        "1 + P / K_EP_to_E_P + S / K_ES_to_E_S" *
        " + R * S / (K_ERS_to_ES_R * K_ES_to_E_S)",
    ))
end

# 21. Essential activator: R must bind before S can bind
#     Forms: E, E_R, E_S_R, E_P_R; SS: E_S_R↔E_P_R
#     Num: R/K4 * (k2f*S/K1 - k2r*P/K3)
#     Denom: 1 + R/K4*(1+S/K1+P/K3)
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: R
        steps: begin
            E(R) + S ⇌ E(S, R)    # K1
            E(S, R) <--> E(P, R)  # k2f, k2r (SS)
            E(R) + P ⇌ E(P, R)    # K3
            E + R ⇌ E(R)          # K4
        end
    end

    # v = Et * R/K4 * (k2f*S/K1 - k2r*P/K3) / (1 + R/K4*(1+S/K1+P/K3))
    function _testhelper_rate_essential_activator(params, concs)
        (; K1, k2f, k2r, K3, K4, Et) = params
        (; S, P, R) = concs
        num = (R / K4) * (k2f * S / K1 - k2r * P / K3)
        denom = 1.0 + (R / K4) * (1.0 + S / K1 + P / K3)
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Essential Activator",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P, :R],
        expected_n_states=4,
        expected_n_steps=4,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=4,
        analytical_rate_fn=(p, c) -> _testhelper_rate_essential_activator(
            merge(p, (Et=p.Et,)), c),
        analytical_kcat_fn=p -> p.k2f,
        expected_factored_num=
        "k_ERS_to_EPR * R * S / (K_ERS_to_ER_S * K_ER_to_E_R) - k_EPR_to_ERS * P * R / (K_EPR_to_ER_P * K_ER_to_E_R)",
        expected_factored_denom=
        "1 + R / K_ER_to_E_R + P * R / (K_EPR_to_ER_P * K_ER_to_E_R) + R * S / (K_ERS_to_ER_S * K_ER_to_E_R)",
    ))
end

# 22. Non-essential activator (general modifier):
#     R modifies catalysis but isn't required. Two parallel SS cycles.
#     Forms: E, E_S, E_P, E_R, E_S_R, E_P_R; SS: E_S↔E_P, E_S_R↔E_P_R
#     K8=K7, K9=K7 (R binding independent of S/P)
#     K4=K1 and K6=K3 are implied by Wegscheider relations
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: R
        steps: begin
            E + S ⇌ E(S)            # K1
            E(S) <--> E(P)          # k2f, k2r (SS)
            E + P ⇌ E(P)            # K3
            E(R) + S ⇌ E(S, R)      # K4
            E(S, R) <--> E(P, R)    # k5f, k5r (SS)
            E(R) + P ⇌ E(P, R)      # K6
            (E + R ⇌ E(R),          # K7 = K8 = K9 (R binding shared)
             E(S) + R ⇌ E(S, R),
             E(P) + R ⇌ E(P, R))
        end
    end

    # v = Et * ((k2f*S/K1-k2r*P/K3) + R/K7*(k5f*S/K1-k5r*P/K3))
    #         / ((1+S/K1+P/K3)*(1+R/K7))
    function _testhelper_rate_nonessential_activator(params, concs)
        (; K1, k2f, k2r, K3, k5f, k5r, K7, Et) = params
        (; S, P, R) = concs
        num = (k2f * S / K1 - k2r * P / K3) +
              (R / K7) * (k5f * S / K1 - k5r * P / K3)
        denom = (1.0 + S / K1 + P / K3) * (1.0 + R / K7)
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Non-essential Activator",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P, :R],
        expected_n_states=6,
        expected_n_steps=9,
        expected_n_haldane_constraints=2,
        expected_n_mirror_constraints=2,
        expected_n_independent_params=5,
        analytical_rate_fn=(p, c) -> _testhelper_rate_nonessential_activator(
            merge(p, (Et=p.Et,)), c),
        analytical_kcat_fn=p -> max(p.k2f, p.k5f),
        expected_factored_num=
        "k_ES_to_EP * S / K_ES_to_E_S + k_ERS_to_EPR * R * S / (K_ERS_to_ER_S * K_ER_to_E_R) - (k_EP_to_ES * P / K_EP_to_E_P + k_EPR_to_ERS * P * R / (K_EP_to_E_P * K_ER_to_E_R))",
        expected_factored_denom=
        "1 + P / K_EP_to_E_P + R / K_ER_to_E_R + S / K_ES_to_E_S + P * R / (K_EP_to_E_P * K_ER_to_E_R) + R * S / (K_ERS_to_ER_S * K_ER_to_E_R)",
    ))
end

# 23. Non-essential activator + competitive inhibitor:
#     Combines non-essential activation with competitive inhibition.
#     A modifies catalysis but isn't required (binds E, E_S, E_P with same K).
#     I binds only free E (competitive dead-end).
#     Forms: E, E_S, E_P, E_A, E_S_A, E_P_A, E_I
#     SS: E_S↔E_P, E_S_A↔E_P_A
#     Constraints: K8=K7, K9=K7 (A binding K independent of S/P)
#     Wegscheider gives K4=K1, K6=K3
#     Denom: (1+S/K1+P/K3)*(1+A/K7) + I/K10
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: A, I
        steps: begin
            E + S ⇌ E(S)            # K1
            E(S) <--> E(P)          # k2f, k2r (SS)
            E + P ⇌ E(P)            # K3
            E(A) + S ⇌ E(A, S)      # K4
            E(A, S) <--> E(A, P)    # k5f, k5r (SS)
            E(A) + P ⇌ E(A, P)      # K6
            (E + A ⇌ E(A),          # K7 = K8 = K9 (A binding shared)
             E(S) + A ⇌ E(A, S),
             E(P) + A ⇌ E(A, P))
            E + I ⇌ E(I)            # K10
        end
    end

    # v = Et * [(k2f*S/K1 - k2r*P/K3) + (A/K7)*(k5f*S/K1 - k5r*P/K3)]
    #         / [(1+S/K1+P/K3)*(1+A/K7) + I/K10]
    function _testhelper_rate_activator_inhibitor(params, concs)
        (; K1, k2f, k2r, K3, k5f, k5r, K7, K10, Et) = params
        (; S, P, A, I) = concs
        num = (k2f * S / K1 - k2r * P / K3) +
              (A / K7) * (k5f * S / K1 - k5r * P / K3)
        denom = (1.0 + S / K1 + P / K3) * (1.0 + A / K7) +
                I / K10
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Activator + Competitive Inhibitor",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P, :A, :I],
        expected_n_states=7,
        expected_n_steps=10,
        expected_n_haldane_constraints=2,
        expected_n_mirror_constraints=2,
        expected_n_independent_params=6,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_activator_inhibitor(
                merge(p, (Et=p.Et,)), c),
        # Denom has both multiplicative (activator) and additive
        # (inhibitor) structure
        expected_factored_num=
        "k_ES_to_EP * S / K_ES_to_E_S + k_EAS_to_EAP * A * S / (K_EAS_to_EA_S * K_EA_to_E_A) - (k_EP_to_ES * P / K_EP_to_E_P + k_EAP_to_EAS * A * P / (K_EAP_to_EA_P * K_EA_to_E_A))",
        expected_factored_denom=
        "1 + A / K_EA_to_E_A + I / K_EI_to_E_I + P / K_EP_to_E_P + S / K_ES_to_E_S + A * P / (K_EAP_to_EA_P * K_EA_to_E_A) + A * S / (K_EAS_to_EA_S * K_EA_to_E_A)",
    ))
end

# 24. MWC (Monod-Wyman-Changeux) dimer with R/T conformational states
#     Two conformational states: R (relaxed, active) and T (tense, inactive).
#     Two identical substrate binding sites per state. Only R state is
#     catalytically active (SS isomerization S→P at each site).
#     R ⇌ T isomerization at base form only (classical MWC).
#     R forms: R_00, R_S0, R_0S, R_SS, R_P0, R_0P, R_SP, R_PS, R_PP (9)
#     T forms: T_00, T_S0, T_0S, T_SS, T_P0, T_0P, T_SP, T_PS, T_PP (9)
#     Total: 18 forms, 23 steps
#     Denominator: (1 + S/Ks_R + P/Kp_R)^2 + L*(1 + S/Ks_T + P/Kp_T)^2
#     Numerator: (
#                   (k_cat_R*S/Ks_R - k_rev_R*P/Kp_R)*(1 + S/Ks_R + P/Kp_R) +
#                   L*(k_cat_T*S/Ks_T - k_rev_T*P/Kp_T)*(1 + S/Ks_T + P/Kp_T)
#                 )
# The explicit flat homodimer fixture is intentionally omitted; the
# AllostericEnzymeMechanism form below is the canonical representation.

# 24B. AllostericEnzymeMechanism MWC Dimer
#      2 catalytic sites × 2 conformations (R/T). No explicit equality constraints
#      needed — symmetric subunits are captured by the site multiplicity.
#      Conformational equilibrium: L (= K37 in the EnzymeMechanism above).
let
    m, raw = @_testhelper_with_raw @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)    :: NonequalAI
            E + P ⇌ E(P)    :: NonequalAI
            E(S) <--> E(P)  :: NonequalAI
        end
    end

    function _testhelper_rate_mwc_dimer_oligo(params, concs)
        (; K1, K2, k3f, k3r, K1_T, K2_T, k3f_T, k3r_T, L, Et) = params
        (; S, P) = concs
        r_flux   = k3f * S / K1 - k3r * P / K2
        t_flux   = k3f_T * S / K1_T - k3r_T * P / K2_T
        r_factor = 1.0 + S / K1 + P / K2
        t_factor = 1.0 + S / K1_T + P / K2_T
        return Et * (r_flux * r_factor + L * t_flux * t_factor) /
                   (r_factor^2 + L * t_factor^2)
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="MWC Dimer [AllostericEnzymeMechanism]",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P],
        expected_n_states=3,          # catalytic subunit: E, E(S), E(P)
        expected_n_steps=3,           # 2 RE + 1 SS per subunit
        expected_n_haldane_constraints=2,         # k3r per conformation × 2
        expected_n_independent_params=7,
        run_ode_test=false,
        analytical_rate_fn=_testhelper_rate_mwc_dimer_oligo,
    ))
end

# ── Edge-case factoring tests ─────────────────────────────────────────────
# These mechanisms test factoring patterns not covered by the classical
# inhibitor/activator mechanisms above.

# The explicit flat homodimer inhibitor fixture is intentionally omitted;
# the AllostericEnzymeMechanism form below is canonical.

# 25B. AllostericEnzymeMechanism Homodimer + Non-competitive Inhibitor
#      I binds all enzyme forms independently with the same Ki (enzyme-level).
#      sigma = Q_cat^2 * (1 + I/K_I_reg1)  (multiplicative factor).
#      The Wegscheider closure of the R_00I ⇌ T_00I cycle is automatic in the
#      conformational assembly formula; the DSL needs no explicit constraint.
let
    m, raw = @_testhelper_with_raw @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: I
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)    :: NonequalAI
            E + P ⇌ E(P)    :: NonequalAI
            E(S) <--> E(P)  :: NonequalAI
        end
        regulatory_site(1): I::NonequalAI
    end

    function _testhelper_rate_homodimer_noncomp_inh_oligo(params, concs)
        (; K1, K2, k3f, k3r, K1_T, K2_T, k3f_T, k3r_T,
           L, K_I_reg1, K_I_T_reg1, Et) = params
        (; S, P, I) = concs
        r_flux   = k3f * S / K1 - k3r * P / K2
        t_flux   = k3f_T * S / K1_T - k3r_T * P / K2_T
        r_factor = 1.0 + S / K1 + P / K2
        t_factor = 1.0 + S / K1_T + P / K2_T
        num   = r_flux * r_factor * (1.0 + I / K_I_reg1) +
                L * t_flux * t_factor * (1.0 + I / K_I_T_reg1)
        denom = r_factor^2 * (1.0 + I / K_I_reg1) +
                L * t_factor^2 * (1.0 + I / K_I_T_reg1)
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Homodimer + Non-competitive Inhibitor [AllostericEnzymeMechanism]",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P, :I],
        expected_n_states=3,
        expected_n_steps=3,
        expected_n_haldane_constraints=2,
        expected_n_independent_params=9,
        run_ode_test=false,
        analytical_rate_fn=_testhelper_rate_homodimer_noncomp_inh_oligo,
    ))
end

# 27. Two Competitive Inhibitors (monomer)
#     Tests: multiple additive dead-end terms (flat sum denominator).
#     I1 and I2 both bind only free E (competitive with S/P and each other).
#     Denom: 1 + S/K1 + P/K3 + I1/K4 + I2/K5
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I1, I2
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
            E + I1 ⇌ E(I1)
            E + I2 ⇌ E(I2)
        end
    end

    function _testhelper_rate_two_comp_inh(p, c)
        (; K1, k2f, k2r, K3, K4, K5, Et) = p
        (; S, P, I1, I2) = c
        num = k2f * S / K1 - k2r * P / K3
        denom = 1.0 + S / K1 + P / K3 + I1 / K4 + I2 / K5
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Two Competitive Inhibitors",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P, :I1, :I2],
        expected_n_states=5,
        expected_n_steps=5,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=5,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_two_comp_inh(merge(p, (Et=p.Et,)), c),
        expected_factored_num=
        "k_ES_to_EP * S / K_ES_to_E_S - k_EP_to_ES * P / K_EP_to_E_P",
        expected_factored_denom=
        "1 + I1 / K_EI1_to_E_I1 + I2 / K_EI2_to_E_I2 + P / K_EP_to_E_P" *
        " + S / K_ES_to_E_S",
    ))
end

# 28. Two Non-competitive Inhibitors (monomer)
#     Tests: triple multiplicative product factoring.
#     I1 and I2 bind independently to all forms (E, ES, EP) at separate sites.
#     Denom: (1+S/K1+P/K3) * (1+I1/K4) * (1+I2/K7)
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I1, I2
        steps: begin
            # S binding (shared K1)
            (E + S ⇌ E(S),
             E(I1) + S ⇌ E(S, I1),
             E(I2) + S ⇌ E(S, I2),
             E(I1, I2) + S ⇌ E(S, I1, I2))
            E(S) <--> E(P)
            # P binding (shared K3)
            (E + P ⇌ E(P),
             E(I1) + P ⇌ E(P, I1),
             E(I2) + P ⇌ E(P, I2),
             E(I1, I2) + P ⇌ E(P, I1, I2))
            # I1 binding (shared K4)
            (E + I1 ⇌ E(I1),
             E(S) + I1 ⇌ E(S, I1),
             E(P) + I1 ⇌ E(P, I1),
             E(I2) + I1 ⇌ E(I1, I2),
             E(S, I2) + I1 ⇌ E(S, I1, I2),
             E(P, I2) + I1 ⇌ E(P, I1, I2))
            # I2 binding (shared K7)
            (E + I2 ⇌ E(I2),
             E(S) + I2 ⇌ E(S, I2),
             E(P) + I2 ⇌ E(P, I2),
             E(I1) + I2 ⇌ E(I1, I2),
             E(S, I1) + I2 ⇌ E(S, I1, I2),
             E(P, I1) + I2 ⇌ E(P, I1, I2))
        end
    end

    # Each written step i gives K{i} or k{i}f/k{i}r; the law reads:
    # K1=S-binding, k5f=iso, K6=P-binding, K10=I1-binding, K16=I2-binding.
    function _testhelper_rate_two_noncomp_inh(p, c)
        (; K1, k5f, k5r, K6, K10, K16, Et) = p
        (; S, P, I1, I2) = c
        num = k5f * S / K1 - k5r * P / K6
        denom = (1.0 + S / K1 + P / K6) *
                (1.0 + I1 / K10) * (1.0 + I2 / K16)
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Two Non-competitive Inhibitors",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P, :I1, :I2],
        expected_n_states=12,
        expected_n_steps=21,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=5,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_two_noncomp_inh(merge(p, (Et=p.Et,)), c),
        expected_factored_num=
        "k_ES_to_EP * S / K_ES_to_E_S - k_EP_to_ES * P / K_EP_to_E_P",
        expected_factored_denom=
        "1 + I1 / K_EI1_to_E_I1 + I2 / K_EI2_to_E_I2 + P / K_EP_to_E_P + S / K_ES_to_E_S + I1 * I2 / (K_EI1_to_E_I1 * K_EI2_to_E_I2) + I1 * P / (K_EI1_to_E_I1 * K_EP_to_E_P) + I1 * S / (K_EI1_to_E_I1 * K_ES_to_E_S) + I2 * P / (K_EI2_to_E_I2 * K_EP_to_E_P) + I2 * S / (K_EI2_to_E_I2 * K_ES_to_E_S) + I1 * I2 * P / (K_EI1_to_E_I1 * K_EI2_to_E_I2 * K_EP_to_E_P) + I1 * I2 * S / (K_EI1_to_E_I1 * K_EI2_to_E_I2 * K_ES_to_E_S)",
    ))
end

# 29. Non-competitive + Competitive Inhibitor (monomer)
#     Tests: multiplicative product + additive dead-end term.
#     I1 binds all forms independently (non-competitive).
#     I2 binds only free E (competitive).
#     Denom: (1+S/K1+P/K3)*(1+I1/K4) + I2/K9
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I1, I2
        steps: begin
            # S binding (shared K1)
            (E + S ⇌ E(S), E(I1) + S ⇌ E(S, I1))
            E(S) <--> E(P)
            # P binding (shared K3)
            (E + P ⇌ E(P), E(I1) + P ⇌ E(P, I1))
            # I1 binding (shared K4)
            (E + I1 ⇌ E(I1),
             E(S) + I1 ⇌ E(S, I1),
             E(P) + I1 ⇌ E(P, I1))
            E + I2 ⇌ E(I2)
        end
    end

    # Each written step i gives K{i} or k{i}f/k{i}r; the law reads:
    # K1=S, k3f=iso, K4=P, K6=I1-binding, K9=I2-dead-end.
    function _testhelper_rate_noncomp_comp_inh(p, c)
        (; K1, k3f, k3r, K4, K6, K9, Et) = p
        (; S, P, I1, I2) = c
        num = k3f * S / K1 - k3r * P / K4
        denom = (1.0 + S / K1 + P / K4) *
                (1.0 + I1 / K6) + I2 / K9
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Non-competitive + Competitive Inhibitor",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P, :I1, :I2],
        expected_n_states=7,
        expected_n_steps=9,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=5,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_noncomp_comp_inh(merge(p, (Et=p.Et,)), c),
        expected_factored_num=
        "k_ES_to_EP * S / K_ES_to_E_S - k_EP_to_ES * P / K_EP_to_E_P",
        expected_factored_denom=
        "1 + I1 / K_EI1_to_E_I1 + I2 / K_EI2_to_E_I2 + P / K_EP_to_E_P + S / K_ES_to_E_S + I1 * P / (K_EI1_to_E_I1 * K_EP_to_E_P) + I1 * S / (K_EI1_to_E_I1 * K_ES_to_E_S)",
    ))
end

# 30. Uncompetitive + Competitive Inhibitor (monomer)
#     Tests: mixed additive structure with nested multiplicative term.
#     I1 binds only ES (uncompetitive). I2 binds only free E (competitive).
#     Denom: 1 + I2/K5 + P/K3 + (S/K1)*(1+I1/K4)
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I1, I2
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
            E(S) + I1 ⇌ E(S, I1)
            E + I2 ⇌ E(I2)
        end
    end

    function _testhelper_rate_uncomp_comp_inh(p, c)
        (; K1, k2f, k2r, K3, K4, K5, Et) = p
        (; S, P, I1, I2) = c
        num = k2f * S / K1 - k2r * P / K3
        denom = 1.0 + I2 / K5 + P / K3 +
                (S / K1) * (1.0 + I1 / K4)
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Uncompetitive + Competitive Inhibitor",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P, :I1, :I2],
        expected_n_states=5,
        expected_n_steps=5,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=5,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_uncomp_comp_inh(merge(p, (Et=p.Et,)), c),
        expected_factored_num=
        "k_ES_to_EP * S / K_ES_to_E_S - k_EP_to_ES * P / K_EP_to_E_P",
        expected_factored_denom=
        "1 + I2 / K_EI2_to_E_I2 + P / K_EP_to_E_P + S / K_ES_to_E_S" *
        " + I1 * S / (K_EI1S_to_ES_I1 * K_ES_to_E_S)",
    ))
end

# 31. Two Same-site (Competing) Non-competitive Inhibitors (monomer)
#     Tests: multiplicative product with flat-sum allosteric factor.
#     I1 and I2 compete for the same allosteric site X (mutually exclusive).
#     Both bind independently of S/P (non-competitive).
#     Denom: (1+S/K1+P/K3) * (1+I1/K4+I2/K9)
let
    m, raw = @_testhelper_with_raw @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I1, I2
        steps: begin
            # S binding (shared K1)
            (E + S ⇌ E(S),
             E(I1) + S ⇌ E(S, I1),
             E(I2) + S ⇌ E(S, I2))
            E(S) <--> E(P)
            # P binding (shared K3)
            (E + P ⇌ E(P),
             E(I1) + P ⇌ E(P, I1),
             E(I2) + P ⇌ E(P, I2))
            # I1 binding (shared K4)
            (E + I1 ⇌ E(I1),
             E(S) + I1 ⇌ E(S, I1),
             E(P) + I1 ⇌ E(P, I1))
            # I2 binding (shared K9)
            (E + I2 ⇌ E(I2),
             E(S) + I2 ⇌ E(S, I2),
             E(P) + I2 ⇌ E(P, I2))
        end
    end

    # Each written step i gives K{i} or k{i}f/k{i}r; the law reads:
    # K1=S, k4f=iso, K5=P, K8=I1, K11=I2.
    function _testhelper_rate_two_samesite_inh(p, c)
        (; K1, k4f, k4r, K5, K8, K11, Et) = p
        (; S, P, I1, I2) = c
        num = k4f * S / K1 - k4r * P / K5
        denom = (1.0 + S / K1 + P / K5) *
                (1.0 + I1 / K8 + I2 / K11)
        return Et * num / denom
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="Two Same-site Inhibitors",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P, :I1, :I2],
        expected_n_states=9,
        expected_n_steps=13,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=5,
        analytical_rate_fn=(p, c) ->
            _testhelper_rate_two_samesite_inh(merge(p, (Et=p.Et,)), c),
        expected_factored_num=
        "k_ES_to_EP * S / K_ES_to_E_S - k_EP_to_ES * P / K_EP_to_E_P",
        expected_factored_denom=
        "1 + I1 / K_EI1_to_E_I1 + I2 / K_EI2_to_E_I2 + P / K_EP_to_E_P + S / K_ES_to_E_S + I1 * P / (K_EI1_to_E_I1 * K_EP_to_E_P) + I1 * S / (K_EI1_to_E_I1 * K_ES_to_E_S) + I2 * P / (K_EI2_to_E_I2 * K_EP_to_E_P) + I2 * S / (K_EI2_to_E_I2 * K_ES_to_E_S)",
    ))
end

# 32. MWC Tetramer — Random-Order Bi-Bi RE + Two Allosteric Sites
#
# Reaction: S1 + S2 ⇌ P1 + P2  (Bi-Bi)
#
# Catalytic site (×4 subunits, 2 conformations):
#   Site A: S1 and P1 bind competitively (Kd = K1, K2 in R state)
#   Site B: S2 and P2 bind competitively (Kd = K3, K4 in R state)
#   Sites A and B are independent — cross-occupancy (S1+P2, P1+S2) is allowed.
#   All binding steps are RE; only the ternary complex interconversion is SS:
#     E_S1S2 <--> E_P1P2  (k13f, k13r; k13r derived from shared Keq via Haldane)
#   Wegscheider constraints enforce site independence:
#     K5=K3, K6=K1  (two routes to E_S1S2)
#     K7=K4, K8=K1  (two routes to E_S1P2)
#     K9=K3, K10=K2 (two routes to E_P1S2)
#     K11=K4, K12=K2 (two routes to E_P1P2)
#
# Reg site type 1 (×4 per tetramer): R1 and R2 bind competitively
# Reg site type 2 (×4 per tetramer): R3 binds exclusively
#
# Analytical rate equation (all Kd = dissociation constants):
#
#   R-state catalytic site:
#     Q_A_R   = 1 + S1/K1 + P1/K2
#     Q_B_R   = 1 + S2/K3 + P2/K4
#     Q_cat_R = Q_A_R * Q_B_R          (site A ⊗ site B independence)
#     N_cat_R = k13f*S1*S2/(K1*K3) - k13r*P1*P2/(K2*K4)
#
#   T-state: same with _T suffix on every K and k parameter.
#
#   Regulatory partition functions (star topology → direct sum):
#     Q_reg1_R = 1 + R1/K_R1_reg1 + R2/K_R2_reg1
#     Q_reg2_R = 1 + R3/K_R3_reg2
#     (T-state: K_R1_T_reg1, K_R2_T_reg1, K_R3_T_reg2)
#
#   MWC assembly (n_cat=4, n_reg1=4, n_reg2=4):
#     Q_R = Q_cat_R^4 * Q_reg1_R^4 * Q_reg2_R^4
#     Q_T = Q_cat_T^4 * Q_reg1_T^4 * Q_reg2_T^4
#     Z   = Q_R + L * Q_T
#     Rate = Et * (N_cat_R * Q_cat_R^3 * Q_reg1_R^4 * Q_reg2_R^4
#                + L * N_cat_T * Q_cat_T^3 * Q_reg1_T^4 * Q_reg2_T^4) / Z
#
# Independent parameters (17, excluding Keq and Et):
#   R-state catalytic: K1, K2, K3, K4, k13f
#   T-state catalytic: K1_T, K2_T, K3_T, K4_T, k13f_T
#   Reg site 1 R: K_R1_reg1, K_R2_reg1
#   Reg site 1 T: K_R1_T_reg1, K_R2_T_reg1
#   Reg site 2 R: K_R3_reg2
#   Reg site 2 T: K_R3_T_reg2
#   Conformational equilibrium: L
let
    m, raw = @_testhelper_with_raw @allosteric_mechanism begin
        substrates: S1, S2
        products: P1, P2
        allosteric_regulators: R1, R2, R3
        catalytic_multiplicity: 4
        steps: begin
            # S1 binding (shared K)
            (E + S1 ⇌ E(S1),
             E(S2) + S1 ⇌ E(S1, S2),
             E(P2) + S1 ⇌ E(S1, P2)) :: NonequalAI
            # P1 binding (shared K)
            (E + P1 ⇌ E(P1),
             E(S2) + P1 ⇌ E(P1, S2),
             E(P2) + P1 ⇌ E(P1, P2)) :: NonequalAI
            # S2 binding (shared K)
            (E + S2 ⇌ E(S2),
             E(S1) + S2 ⇌ E(S1, S2),
             E(P1) + S2 ⇌ E(P1, S2)) :: NonequalAI
            # P2 binding (shared K)
            (E + P2 ⇌ E(P2),
             E(S1) + P2 ⇌ E(S1, P2),
             E(P1) + P2 ⇌ E(P1, P2)) :: NonequalAI
            E(S1, S2) <--> E(P1, P2) :: NonequalAI
        end
        regulatory_site(4): R1::NonequalAI, R2::NonequalAI
        regulatory_site(4): R3::NonequalAI
    end

    # Parameter naming convention:
    #   R-state catalytic Kd: K1=Kd(S1), K4=Kd(P1), K7=Kd(S2), K10=Kd(P2)
    #   R-state SS rate:      k13f (k13r is Haldane-derived)
    #   T-state catalytic:    K1_T, K4_T, K7_T, K10_T, k13f_T, k13r_T
    #   Reg site 1, R state:  K_R1_reg1, K_R2_reg1
    #   Reg site 1, T state:  K_R1_T_reg1, K_R2_T_reg1
    #   Reg site 2, R state:  K_R3_reg2
    #   Reg site 2, T state:  K_R3_T_reg2
    #   Shared: L (= [E_T]/[E_R] for bare enzyme), Keq, Et
    #
    # Each written step i gives K{i} or k{i}f/k{i}r; the law reads the first step of
    # each group: K1 = S1 binding (step 1), K4 = P1 binding (step 4), K7 = S2 binding
    # (step 7), K10 = P2 binding (step 10), k13f/k13r = catalysis (step 13).
    function _testhelper_rate_mwc_tetramer_bi_bi(params, concs)
        (; K1, K4, K7, K10, k13f, k13r,
           K1_T, K4_T, K7_T, K10_T, k13f_T, k13r_T,
           K_R1_reg1, K_R2_reg1, K_R1_T_reg1, K_R2_T_reg1,
           K_R3_reg2, K_R3_T_reg2,
           L, Et) = params
        (; S1, S2, P1, P2, R1, R2, R3) = concs

        # R-state: catalytic site factors by site-A ⊗ site-B independence
        Q_A_R   = 1.0 + S1 / K1 + P1 / K4
        Q_B_R   = 1.0 + S2 / K7 + P2 / K10
        Q_cat_R = Q_A_R * Q_B_R
        N_cat_R = k13f * S1 * S2 / (K1 * K7) - k13r * P1 * P2 / (K4 * K10)

        # R-state: regulatory site partition functions (star topology → direct sum)
        Q_reg1_R = 1.0 + R1 / K_R1_reg1 + R2 / K_R2_reg1
        Q_reg2_R = 1.0 + R3 / K_R3_reg2

        # T-state: same structure with _T parameters
        Q_A_T   = 1.0 + S1 / K1_T + P1 / K4_T
        Q_B_T   = 1.0 + S2 / K7_T + P2 / K10_T
        Q_cat_T = Q_A_T * Q_B_T
        N_cat_T = k13f_T * S1 * S2 / (K1_T * K7_T) - k13r_T * P1 * P2 / (K4_T * K10_T)

        Q_reg1_T = 1.0 + R1 / K_R1_T_reg1 + R2 / K_R2_T_reg1
        Q_reg2_T = 1.0 + R3 / K_R3_T_reg2

        # MWC assembly (n_cat=4, n_reg1=4, n_reg2=4)
        Q_R = Q_cat_R^4 * Q_reg1_R^4 * Q_reg2_R^4
        Q_T = Q_cat_T^4 * Q_reg1_T^4 * Q_reg2_T^4
        Z   = Q_R + L * Q_T

        num = N_cat_R * Q_cat_R^3 * Q_reg1_R^4 * Q_reg2_R^4 +
              L * N_cat_T * Q_cat_T^3 * Q_reg1_T^4 * Q_reg2_T^4

        return Et * num / Z
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S1, :S2, :P1, :P2, :R1, :R2, :R3],
        expected_n_states=9,           # catalytic subunit states
        expected_n_steps=13,           # catalytic subunit steps
        expected_n_haldane_constraints=2,          # one k13r per conformation (R and T)
        expected_n_independent_params=17,
        run_ode_test=false,
        analytical_rate_fn=_testhelper_rate_mwc_tetramer_bi_bi,
    ))
end

# ── PFK-1 hand-verified mechanism ───────────────────────────────────────
# Reaction: F6P + ATP ⇌ F16BP + ADP, 4 catalytic subunits, 2 conformations.
# F6P binding is :OnlyA — T-state can't bind F6P, so the T-state cycle is
# broken in both directions and N_cat_T = 0. ATP appears as both
# substrate and allosteric regulator (different tags per context).
let
    m, raw = @_testhelper_with_raw @allosteric_mechanism begin
        substrates: F6P, ATP
        products:   F16BP, ADP
        allosteric_regulators: Pi, ATP, ADP, Citrate, F26BP

        catalytic_multiplicity: 4
        steps: begin
            (E + F6P ⇌ E(F6P), E(ATP) + F6P ⇌ E(F6P, ATP))           :: OnlyA
            (E + ATP ⇌ E(ATP), E(F6P) + ATP ⇌ E(F6P, ATP))           :: EqualAI
            E(F6P, ATP) <--> E(F16BP, ADP)                            :: OnlyA
            (E(F16BP, ADP) ⇌ E(ADP) + F16BP, E(F16BP) ⇌ E + F16BP)   :: EqualAI
            (E(F16BP, ADP) ⇌ E(F16BP) + ADP, E(ADP) ⇌ E + ADP)       :: EqualAI
        end

        regulatory_site(4): Pi::EqualAI, ATP::OnlyI
        regulatory_site(4): ADP::OnlyA
        regulatory_site(4): Citrate::OnlyI
        regulatory_site(4): F26BP::NonequalAI
    end

    function _testhelper_pfk_rate_analytical(params, concs)
        (; K1, K3, k5f, K6, K8,
           K_Pi_reg1, K_ATP_T_reg1,
           K_ADP_reg2, K_Citrate_T_reg3,
           K_F26BP_reg4, K_F26BP_T_reg4,
           L, Keq, Et) = params
        (; F6P, ATP, F16BP, ADP, Pi, Citrate, F26BP) = concs
        k5r = k5f * K6 * K8 / (Keq * K1 * K3)

        Q_cat_R = 1 + F6P/K1 + ATP/K3 + F6P*ATP/(K1*K3) +
                  F16BP/K6 + ADP/K8 + F16BP*ADP/(K6*K8)
        Q_cat_T = 1 + ATP/K3 + F16BP/K6 + ADP/K8 + F16BP*ADP/(K6*K8)

        # N_cat_T = 0: T-state cycle is broken (F6P binding :OnlyA), so
        # the Cha-nominal reverse term -k5r*F16BP*ADP/(K6*K8) is
        # non-physical at steady state.
        N_cat_R = k5f * F6P * ATP / (K1 * K3) - k5r * F16BP * ADP / (K6 * K8)
        N_cat_T = 0.0

        Q_reg1_R = 1 + Pi / K_Pi_reg1
        Q_reg1_T = 1 + Pi / K_Pi_reg1 + ATP / K_ATP_T_reg1
        Q_reg2_R = 1 + ADP / K_ADP_reg2
        Q_reg2_T = 1
        Q_reg3_R = 1
        Q_reg3_T = 1 + Citrate / K_Citrate_T_reg3
        Q_reg4_R = 1 + F26BP / K_F26BP_reg4
        Q_reg4_T = 1 + F26BP / K_F26BP_T_reg4

        Q_R = Q_cat_R^4 * Q_reg1_R^4 * Q_reg2_R^4 * Q_reg3_R^4 * Q_reg4_R^4
        Q_T = Q_cat_T^4 * Q_reg1_T^4 * Q_reg2_T^4 * Q_reg3_T^4 * Q_reg4_T^4
        num = N_cat_R * Q_cat_R^3 * Q_reg1_R^4 * Q_reg2_R^4 * Q_reg3_R^4 * Q_reg4_R^4 +
              L * N_cat_T * Q_cat_T^3 * Q_reg1_T^4 * Q_reg2_T^4 * Q_reg3_T^4 * Q_reg4_T^4
        return Et * num / (Q_R + L * Q_T)
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="PFK-1",
        mechanism=m,
        raw=raw,
        metabolite_names=[:ATP, :F6P, :ADP, :F16BP, :Pi, :Citrate, :F26BP],
        expected_n_states=7,
        expected_n_steps=9,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=12,
        run_ode_test=false,
        analytical_rate_fn=_testhelper_pfk_rate_analytical,
        # F6P binding (group 1) is :OnlyA → the F6P·ATP saturating
        # pattern is unreachable in T-state, so kcat = k5f
        # for every regulator corner. Guards the peak turnover for a saturating
        # pattern absent from the I-state.
        analytical_kcat_fn = p -> p.k5f,
    ))
end

# ── Hexokinase hand-verified mechanism ──────────────────────────────────
#
# G6P has THREE distinct binding sites on HK:
#
#   1. Catalytic site (glucose pocket): G6P competes with Glucose, but
#      binds together with ATP at the nucleotide pocket. This is the
#      product-release equilibrium — group 4 (K7).
#
#   2. Inhibitory site (= nucleotide pocket): G6P competes with ATP/ADP
#      at the same pocket, but can co-bind with Glucose or with G6P at
#      the catalytic site — group 6 (K12). G6P here EXCLUDES ATP/ADP at
#      the same site, so there are no E_ATP_G6Pi or E_ADP_G6Pi forms.
#
#   3. Allosteric site: G6P competes with Pi at a separate regulatory
#      site — regulatory_site(2) with G6P::OnlyI and Pi::EqualAI.
#
# ATP binding (group 2) is :OnlyA — T-state can't bind ATP, so the
# catalytic cycle is broken and N_cat_T = 0.
let
    m, raw = @_testhelper_with_raw @allosteric_mechanism begin
        substrates: Glucose, ATP
        products:   G6P, ADP
        allosteric_regulators: G6P, Pi
        modifiers:  G6P

        catalytic_multiplicity: 2
        steps: begin
            # Group 1 (Glucose binding at catalytic site, EqualAI)
            (E + Glucose ⇌ E(Glucose),
             E(ATP) + Glucose ⇌ E(Glucose, ATP),
             E(G6P::Mod) + Glucose ⇌ E(Glucose, G6P::Mod))    :: EqualAI
            # Group 2 (ATP binding at nucleotide pocket, OnlyA —
            # T-state can't bind ATP)
            (E + ATP ⇌ E(ATP),
             E(Glucose) + ATP ⇌ E(Glucose, ATP))              :: OnlyA
            # Group 3 (catalysis SS, OnlyA)
            E(Glucose, ATP) <--> E(G6P, ADP)                  :: OnlyA
            # Group 4 (G6P binding/release at catalytic site, EqualAI)
            (E(G6P, ADP) ⇌ E(ADP) + G6P,
             E(G6P) ⇌ E + G6P,
             E(G6P, G6P::Mod) ⇌ E(G6P::Mod) + G6P)            :: EqualAI
            # Group 5 (ADP release, EqualAI)
            (E(G6P, ADP) ⇌ E(G6P) + ADP,
             E(ADP) ⇌ E + ADP)                                :: EqualAI
            # Group 6 (G6P binding at INHIBITORY site, EqualAI) —
            # G6P at site 2 competes with ATP/ADP, can co-bind with
            # Glucose at site 1 or G6P at site 1.
            (E + G6P::Mod ⇌ E(G6P::Mod),
             E(Glucose) + G6P::Mod ⇌ E(Glucose, G6P::Mod),
             E(G6P) + G6P::Mod ⇌ E(G6P, G6P::Mod))            :: EqualAI
        end

        regulatory_site(2): G6P::OnlyI, Pi::EqualAI
    end

    # Each written step i gives K{i} or k{i}f/k{i}r; the law reads:
    #   K1  (Glucose binding, group 1, step 1)
    #   K4  (ATP binding, group 2, step 4)
    #   k6f (catalysis, group 3, step 6)
    #   K7  (G6P at catalytic site, group 4, step 7)
    #   K10 (ADP release, group 5, step 10)
    #   K12 (G6P at inhibitory site, group 6, step 12) — single K
    #        for all three E_G6Pi-form bindings.
    function _testhelper_hk_rate_analytical(params, concs)
        (; K1, K4, k6f, K7, K10, K12,
           K_Pi_reg1, K_G6P_T_reg1,
           L, Keq, Et) = params
        (; Glucose, ATP, G6P, ADP, Pi) = concs
        k6r = k6f * K7 * K10 / (Keq * K1 * K4)

        # R-state catalytic partition function (10 enzyme forms).
        Q_cat_R = 1 +
                  Glucose / K1 +
                  ATP / K4 +
                  Glucose * ATP / (K1 * K4) +
                  G6P * ADP / (K7 * K10) +
                  G6P / K7 +
                  ADP / K10 +
                  G6P / K12 +
                  Glucose * G6P / (K1 * K12) +
                  G6P^2 / (K7 * K12)
        # T-state: ATP group :OnlyA → zero ATP terms.
        Q_cat_T = 1 +
                  Glucose / K1 +
                  G6P * ADP / (K7 * K10) +
                  G6P / K7 +
                  ADP / K10 +
                  G6P / K12 +
                  Glucose * G6P / (K1 * K12) +
                  G6P^2 / (K7 * K12)

        # N_cat_T = 0: T-state cycle is broken (ATP binding :OnlyA).
        N_cat_R = k6f * Glucose * ATP / (K1 * K4) -
                  k6r * G6P * ADP / (K7 * K10)
        N_cat_T = 0.0

        Q_reg1_R = 1 + Pi / K_Pi_reg1
        Q_reg1_T = 1 + Pi / K_Pi_reg1 + G6P / K_G6P_T_reg1

        Q_R = Q_cat_R^2 * Q_reg1_R^2
        Q_T = Q_cat_T^2 * Q_reg1_T^2
        num = N_cat_R * Q_cat_R * Q_reg1_R^2 +
              L * N_cat_T * Q_cat_T * Q_reg1_T^2
        return Et * num / (Q_R + L * Q_T)
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="HK",
        mechanism=m,
        raw=raw,
        metabolite_names=[:ATP, :Glucose, :ADP, :G6P, :Pi],
        expected_n_states=10,    # 7 catalytic-cycle + 3 G6Pi dead-end
        expected_n_steps=14,     # 3+2+1+3+2+3
        expected_n_haldane_constraints=1,
        expected_n_independent_params=9,
        run_ode_test=false,
        analytical_rate_fn=_testhelper_hk_rate_analytical,
        # ATP binding (group 2) is :OnlyA → the Glucose·ATP saturating
        # pattern is unreachable in T-state, so kcat = k6f
        # for every regulator corner. Guards the peak turnover for a saturating
        # pattern absent from the I-state.
        analytical_kcat_fn = p -> p.k6f,
    ))
end

# ── Pyruvate kinase (PK) hand-verified mechanism ──────────────────────────
# Reaction: PEP + ADP ⇌ Pyruvate + ATP, 4 catalytic subunits.
# PEP binding and catalysis are :OnlyA, so the T-state cycle is dead; the
# ADP binding and both releases are :EqualAI. Reg sites have MISMATCHED
# multiplicities:
#   ATP::OnlyI at mult 2
#   F16BP::OnlyA at mult 4 (matches catalytic mult)
# This exercises the symmetric all-reg-sites contribution to both numerator
# and denominator.
# Independent parameters (8): K1, K3, k5f, K6, K8, K_ATP_T_reg1, K_F16BP_reg2, L
let
    m, raw = @_testhelper_with_raw @allosteric_mechanism begin
        substrates: PEP, ADP
        products:   Pyruvate, ATP
        allosteric_regulators: ATP, F16BP

        catalytic_multiplicity: 4
        steps: begin
            (E + PEP ⇌ E(PEP),
             E(ADP) + PEP ⇌ E(PEP, ADP))                          :: OnlyA
            (E + ADP ⇌ E(ADP),
             E(PEP) + ADP ⇌ E(PEP, ADP))                          :: EqualAI
            E(PEP, ADP) <--> E(Pyruvate, ATP)                     :: OnlyA
            (E(Pyruvate, ATP) ⇌ E(ATP) + Pyruvate,
             E(Pyruvate) ⇌ E + Pyruvate)                          :: EqualAI
            (E(Pyruvate, ATP) ⇌ E(Pyruvate) + ATP,
             E(ATP) ⇌ E + ATP)                                    :: EqualAI
        end

        regulatory_site(2): ATP::OnlyI
        regulatory_site(4): F16BP::OnlyA
    end

    # Param mapping:
    #   K1  : PEP binding (group 1, :OnlyA — R-state only)
    #   K3  : ADP binding (group 2, :EqualAI)
    #   k5f : catalysis SS forward rate (group 3, :OnlyA)
    #   K6  : Pyruvate release (group 4, :EqualAI)
    #   K8  : ATP release (group 5, :EqualAI)
    #
    # PEP binding AND catalysis are :OnlyA (the exclusive-binding K-system):
    # the T-state cannot bind PEP, so it cannot catalyze — the T-catalytic
    # cycle is dead (N_T = 0) and the T-state is a clean binding partition with
    # no PEP term. k5r derives via the R-state Haldane k5r = k5f·K6·K8/(Keq·K1·K3).
    # At saturation the R-state dominates, so forward kcat = k5f.
    function _testhelper_pk_rate_analytical(params, concs)
        (; K1, K3, k5f, K6, K8,
           K_ATP_T_reg1, K_F16BP_reg2,
           L, Keq, Et) = params
        (; PEP, ADP, Pyruvate, ATP, F16BP) = concs

        k5r = k5f * K6 * K8 / (Keq * K1 * K3)

        Q_cat_R = 1 + PEP/K1 + ADP/K3 + PEP*ADP/(K1 * K3) +
                  Pyruvate/K6 + ATP/K8 + Pyruvate*ATP/(K6 * K8)
        Q_cat_T = 1 + ADP/K3 + Pyruvate/K6 + ATP/K8 +
                  Pyruvate*ATP/(K6 * K8)                    # no PEP (:OnlyA)

        N_R = k5f * PEP * ADP / (K1 * K3) - k5r * Pyruvate * ATP / (K6 * K8)
        # N_T = 0: catalysis :OnlyA — the dead T-state contributes no flux.

        Q_reg1_R = 1                                     # ATP::OnlyI, no R term
        Q_reg1_T = 1 + ATP / K_ATP_T_reg1
        Q_reg2_R = 1 + F16BP / K_F16BP_reg2               # F16BP::OnlyA, no T term
        Q_reg2_T = 1

        num = N_R * Q_cat_R^3 * Q_reg1_R^2 * Q_reg2_R^4
        den = Q_cat_R^4 * Q_reg1_R^2 * Q_reg2_R^4 +
              L * Q_cat_T^4 * Q_reg1_T^2 * Q_reg2_T^4
        return Et * num / den
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="PK",
        mechanism=m,
        raw=raw,
        metabolite_names=[:ADP, :PEP, :ATP, :Pyruvate, :F16BP],
        expected_n_states=7,           # E, E_PEP, E_ADP, E_PEP_ADP, E_Pyr_ATP, E_Pyr, E_ATP
        expected_n_steps=9,
        expected_n_haldane_constraints=1,
        # PEP binding and catalysis are :OnlyA (the T-state is pruned), so
        # there is no mirror/collapse — the T-catalytic cycle is simply dead
        # (N_T = 0). One R-state Haldane derives k5r; no Wegscheider tie.
        expected_n_independent_params=8,
        run_ode_test=false,
        analytical_rate_fn=_testhelper_pk_rate_analytical,
        analytical_kcat_fn = p -> p.k5f,
    ))
end

# ── m_all: NonequalAI coverage on substrate + catalysis + product ───────
# Two-substrate two-product reaction with explicit :NonequalAI on
# S1 binding (group 1), catalysis (group 3), and P1 release (group 4),
# and a 2-ligand mixed-state reg site (R1::NonequalAI + R2::EqualAI).
# Catalysis :NonequalAI combined with :NonequalAI substrate yields
# independent R-state and T-state Haldanes (one per state).
# Independent parameters (12): K1, K1_T, K3, k5f, k5f_T, K6, K6_T, K8,
# K_R1_reg1, K_R1_T_reg1, K_R2_reg1, L
let
    m, raw = @_testhelper_with_raw @allosteric_mechanism begin
        substrates: S1, S2
        products:   P1, P2
        allosteric_regulators: R1, R2

        catalytic_multiplicity: 2
        steps: begin
            (E + S1 ⇌ E(S1),
             E(S2) + S1 ⇌ E(S1, S2))    :: NonequalAI
            (E + S2 ⇌ E(S2),
             E(S1) + S2 ⇌ E(S1, S2))    :: EqualAI
            E(S1, S2) <--> E(P1, P2)    :: NonequalAI
            (E(P1, P2) ⇌ E(P2) + P1,
             E(P1) ⇌ E + P1)            :: NonequalAI
            (E(P1, P2) ⇌ E(P1) + P2,
             E(P2) ⇌ E + P2)            :: EqualAI
        end

        regulatory_site(2): R1::NonequalAI, R2::EqualAI
    end

    # Each written step i gives K{i} or k{i}f/k{i}r; the law reads:
    #   K1, K1_T : S1 binding (group 1, NonequalAI)
    #   K3       : S2 binding (group 2, EqualAI)
    #   k5f, k5f_T : catalysis SS (group 3, NonequalAI)
    #   K6, K6_T : P1 release (group 4, NonequalAI)
    #   K8       : P2 release (group 5, EqualAI)
    function _testhelper_m_all_rate_analytical(params, concs)
        (; K1, K1_T, K3, k5f, k5f_T, K6, K6_T, K8,
           K_R1_reg1, K_R1_T_reg1, K_R2_reg1,
           L, Keq, Et) = params
        (; S1, S2, P1, P2, R1, R2) = concs

        k5r   = k5f   * K6   * K8 / (Keq * K1   * K3)
        k5r_T = k5f_T * K6_T * K8 / (Keq * K1_T * K3)

        Q_cat_R = 1 + S1/K1   + S2/K3 + S1*S2/(K1   * K3) +
                  P1/K6   + P2/K8 + P1*P2/(K6   * K8)
        Q_cat_T = 1 + S1/K1_T + S2/K3 + S1*S2/(K1_T * K3) +
                  P1/K6_T + P2/K8 + P1*P2/(K6_T * K8)

        N_R = k5f   * S1 * S2 / (K1   * K3) - k5r   * P1 * P2 / (K6   * K8)
        N_T = k5f_T * S1 * S2 / (K1_T * K3) - k5r_T * P1 * P2 / (K6_T * K8)

        Q_reg1_R = 1 + R1/K_R1_reg1   + R2/K_R2_reg1
        Q_reg1_T = 1 + R1/K_R1_T_reg1 + R2/K_R2_reg1

        num = N_R * Q_cat_R   * Q_reg1_R^2 + L * N_T * Q_cat_T   * Q_reg1_T^2
        den = Q_cat_R^2 * Q_reg1_R^2       + L *      Q_cat_T^2 * Q_reg1_T^2

        return Et * num / den
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="m_all",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S1, :S2, :P1, :P2, :R1, :R2],
        expected_n_states=7,
        expected_n_steps=9,
        expected_n_haldane_constraints=2,
        # :EqualAI catalytic groups and regulator ligands each take one bare name: no mirror
        expected_n_independent_params=12,
        run_ode_test=false,
        analytical_rate_fn=_testhelper_m_all_rate_analytical,
        analytical_kcat_fn=nothing,      # cat is :NonequalAI → kcat L-dependent
    ))
end

# ── m_OnlyA_prod: catalysis :OnlyA → T-state truly inactive ─────────────
# Exercises the zero I-state numerator when the catalytic conversion is
# :OnlyA: the T-state cannot run the reaction (forward or reverse), so its
# cycle is dead — N_T is 0 and the L*num_T branch is dropped. Analytical
# kcat = k2f/(1+L) — L-dependent because the saturating R-state
# pattern (S only) IS reachable in T-state (Q_cat_T at sat S = S/K1,
# same as Q_cat_R), so B_T ≠ 0 and the 1/(1+L) factor appears.
let
    m, raw = @_testhelper_with_raw @allosteric_mechanism begin
        substrates: S
        products:   P

        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)    :: EqualAI      # group 1, K1
            E(S) <--> E(P)  :: OnlyA        # group 2, k2f catalysis :OnlyA
            E + P ⇌ E(P)    :: OnlyA        # group 3, K3 (P binding)
        end
    end

    # Param mapping:
    #   K1   : S binding (group 1, EqualAI)
    #   k2f  : catalysis SS (group 2, OnlyA, k2r derived via Haldane)
    #   K3   : P release (group 3, OnlyA)
    function _testhelper_m_OnlyA_prod_rate_analytical(params, concs)
        (; K1, k2f, K3, L, Keq, Et) = params
        (; S, P) = concs

        k2r = k2f * K3 / (Keq * K1)

        Q_cat_R = 1 + S/K1 + P/K3
        Q_cat_T = 1 + S/K1                    # E(P) unreachable in T-state (catalysis + P-binding :OnlyA)

        N_R = k2f * S/K1 - k2r * P/K3
        # N_T = 0: the T-state cycle is dead (catalysis group 2 is :OnlyA)

        num = N_R * Q_cat_R                   # L*N_T*Q_cat_T term elided
        den = Q_cat_R^2 + L * Q_cat_T^2

        return Et * num / den
    end

    # The peak turnover over S and P, from the written constants `p`.
    # kcat at saturating S, zero P:
    #   A_R = k2f/K1², B_R = 1/K1², B_T = 1/K1² (T-state pattern same as R)
    #   kcat = A_R / (B_R + L · B_T) = k2f / (1 + L)
    # P binds the R state only, so it can draw the enzyme out of T and raise the
    # turnover above kcat. At any fixed P/S the turnover rises with S, so the peak
    # lies at saturating S, where, with P/K3 = c·S/K1, it is
    #   g(c) = (k2f − k2r·c)(1 + c) / ((1 + c)² + L).
    # g rises at c = 0 when k2f·(L − 1) > k2r·(1 + L). Its peak over c ≥ 0 sits at
    # 1 + c = u, the positive root of (k2f + k2r)·u² + 2·k2r·L·u − (k2f + k2r)·L = 0,
    # or at c = 0 when u ≤ 1.
    function _testhelper_m_OnlyA_prod_peak_turnover(p)
        g(c) = (p.k2f - p.k2r * c) * (1 + c) / ((1 + c)^2 + p.L)
        s = p.k2f + p.k2r
        u = (sqrt((p.k2r * p.L)^2 + s^2 * p.L) - p.k2r * p.L) / s
        g(max(0.0, u - 1))
    end

    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="m_OnlyA_prod",
        mechanism=m,
        raw=raw,
        metabolite_names=[:S, :P],
        expected_n_states=3,                  # E, E_S, E_P
        expected_n_steps=3,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=4,
        run_ode_test=false,
        analytical_rate_fn=_testhelper_m_OnlyA_prod_rate_analytical,
        analytical_kcat_fn=_testhelper_m_OnlyA_prod_peak_turnover,
    ))
end

# The three LDH i-state mechanisms that exposed the Bug-2 fitted_params
# leak (a Haldane-dependent reverse rate that landed in the independent
# set). Written as @allosteric_mechanism for readability; RE-containing
# with a products=0 boundary, so the ODE cross-check is skipped. The parameter
# headers of their golden RATE_EQUATION_STRING lines confirm the dependent reverse
# rate is excluded.
let
    m, raw = @_testhelper_with_raw @allosteric_mechanism begin
        substrates: NADH, Pyruvate
        products: Lactate, NAD
        catalytic_multiplicity: 4
        steps: begin
            (E + Lactate <--> E(Lactate),
             E(NADH) + Lactate <--> E(Lactate, NADH)) :: NonequalAI
            (E + NAD ⇌ E(NAD),
             E(Lactate) + NAD ⇌ E(Lactate, NAD),
             E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate)) :: EqualAI
            (E + NADH <--> E(NADH),
             E(Lactate) + NADH <--> E(Lactate, NADH),
             E(Pyruvate) + NADH <--> E(NADH, Pyruvate)) :: OnlyA
            (E + Pyruvate ⇌ E(Pyruvate),
             E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
             E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate)) :: EqualAI
            E(NAD) + Lactate ⇌ E(Lactate, NAD) :: EqualAI
            E(NADH, Pyruvate) <--> E(Lactate, NAD) :: OnlyA
        end
    end
    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="LDH i-state NonequalAI 6-group",
        mechanism=m,
        raw=raw,
        metabolite_names=[:NADH, :Pyruvate, :Lactate, :NAD],
        expected_n_states=9, expected_n_steps=13,
        expected_n_haldane_constraints=1,
        expected_n_wegscheider_constraints=2, expected_n_independent_params=9,
        run_ode_test=false))
end

let
    m, raw = @_testhelper_with_raw @allosteric_mechanism begin
        substrates: NADH, Pyruvate
        products: Lactate, NAD
        catalytic_multiplicity: 4
        steps: begin
            (E + Lactate ⇌ E(Lactate),
             E(NAD) + Lactate ⇌ E(Lactate, NAD),
             E(NADH) + Lactate ⇌ E(Lactate, NADH)) :: NonequalAI
            (E + NAD ⇌ E(NAD),
             E(Lactate) + NAD ⇌ E(Lactate, NAD)) :: EqualAI
            (E + NADH <--> E(NADH),
             E(Lactate) + NADH <--> E(Lactate, NADH)) :: OnlyA
            E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate) :: EqualAI
            E(NADH, Pyruvate) <--> E(Lactate, NAD) :: OnlyA
        end
    end
    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="LDH i-state NonequalAI 5-group",
        mechanism=m,
        raw=raw,
        metabolite_names=[:NADH, :Pyruvate, :Lactate, :NAD],
        expected_n_states=7, expected_n_steps=9,
        expected_n_haldane_constraints=1,
        expected_n_independent_params=8,
        run_ode_test=false))
end

let
    m, raw = @_testhelper_with_raw @allosteric_mechanism begin
        substrates: NADH, Pyruvate
        products: Lactate, NAD
        catalytic_multiplicity: 4
        steps: begin
            E + NAD ⇌ E(NAD) :: EqualAI
            (E + NADH ⇌ E(NADH),
             E(Pyruvate) + NADH ⇌ E(NADH, Pyruvate)) :: EqualAI
            (E + Pyruvate ⇌ E(Pyruvate),
             E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
             E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate)) :: EqualAI
            (E(NAD) + Lactate ⇌ E(Lactate, NAD),
             E(NADH) + Lactate ⇌ E(Lactate, NADH)) :: EqualAI
            E(NADH, Pyruvate) <--> E(Lactate, NAD) :: EqualAI
            E(Pyruvate) + NAD <--> E(NAD, Pyruvate) :: NonequalAI
        end
    end
    push!(MECHANISM_TEST_SPECS, _testhelper_MechanismTestSpec(
        name="LDH i-state EqualAI-NonequalAI 6-group",
        mechanism=m,
        raw=raw,
        metabolite_names=[:NADH, :Pyruvate, :Lactate, :NAD],
        expected_n_states=8, expected_n_steps=10,
        expected_n_haldane_constraints=1,
        expected_n_wegscheider_constraints=2, expected_n_independent_params=8,
        run_ode_test=false))
end

# The PFK-1, HK and PK specs above pair one :OnlyA binding with :OnlyA catalysis. A
# mechanism that binds both substrates :OnlyA has no spec here; its derivation is
# checked against the n=1 mass-action ground truth in allosteric_ground_truth.jl
# (testset "multi-OnlyA MWC derivation matches mass-action ground truth").
# `_expand_to_allosteric` reaches it by tagging the two substrate binding groups
# :OnlyA with the chemistry.
